use uuid::Uuid;

use deadpool_postgres::{Pool, Manager};

use crate::api::query::Query;

use std::str::FromStr;

pub struct SiteDB {
    pool: Pool,
}

impl SiteDB {
    pub async fn new(connection_string: &str) -> Result<SiteDB, anyhow::Error> {
        let config = tokio_postgres::Config::from_str(connection_string)?;

        let manager = Manager::new(config, tokio_postgres::NoTls);
        let pool = Pool::new(manager, 16);

        Ok(SiteDB {
            pool,
        })
    }

    pub async fn query_by_id(&self, uuid: &Uuid) -> Option<Query> {
        let client =
            match self.pool.get().await {
                Ok(client) => client,
                Err(err) => {
                    eprintln!("error querying by ID: {}", err);
                    return None;
                },
            };

        let rs =
            match client.query("SELECT query_json FROM query WHERE id = $1",
                               &[uuid]).await {
                Ok(rs) => rs,
                Err(err) => {
                    eprintln!("error querying by ID: {}", err);
                    return None;
                },
            };

        if !rs.is_empty() {
            let query_value: Option<String> = rs[0].get(0);

            let query: Query = match serde_json::from_str(&query_value.unwrap()) {
                Ok(v) => v,
                Err(_) => return None,
            };

            return Some(query);
        }

        None
    }

    pub async fn id_from_query(&self, query: &Query) -> Option<Uuid> {
        let query_value: String = match serde_json::value::to_value(query) {
            Ok(v) => v.to_string(),
            Err(err) => {
                eprintln!("error converting query to string: {:?}", err);
                return None;
            }
        };

        let client =
            match self.pool.get().await {
                Ok(client) => client,
                Err(err) => {
                    eprintln!("can't get client: {:?}", err);
                    return None;
                }
            };

        let rs =
            match client.query("SELECT id::uuid FROM query WHERE digest(query_json, 'sha256') = digest($1, 'sha256');",
                               &[&query_value]).await
        {
            Ok(rs) => rs,
            Err(err) => {
                eprintln!("error querying in id_from_query(): {:?}", err);
                return None;
            }
        };

        if !rs.is_empty() {
            let id: Uuid = rs[0].get("id");
            Some(id)
        } else {
            None
        }
    }

    pub async fn save_query(&self, uuid: &Uuid, query: &Query) -> Result<(), String> {
        let mut client =
            match self.pool.get().await {
                Ok(client) => client,
                Err(err) => {
                    return Err(format!("error saving query: {}", err));
                },
            };

        let trans = match client.transaction().await {
            Ok(t) => t,
            Err(e) => return Err(format!("couldn't begin transaction: {}", e)),
        };

        let serde_value: String = match serde_json::value::to_value(&query) {
            Ok(v) => v.to_string(),
            Err(e) => return Err(format!("serde error: {}", e)),
        };

        match trans.execute("INSERT INTO query(id, query_json) values ($1, $2)",
                            &[uuid, &serde_value]).await {
            Ok(_) => (),
            Err(e) => return Err(format!("error executing query: {}", e)),
        };

        match trans.commit().await {
            Ok(_) => Ok(()),
            Err(e) => Err(format!("error saving query: {}", e)),
        }
    }
}
