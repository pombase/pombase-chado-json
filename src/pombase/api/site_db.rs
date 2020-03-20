use uuid::Uuid;
use postgres::{Connection, TlsMode};

use crate::api::query::Query;

pub struct SiteDB {
    conn: Connection,
}

impl SiteDB {
    pub fn new(connection_string: &str) -> SiteDB {
        let conn = match Connection::connect(connection_string, TlsMode::None) {
            Ok(conn) => conn,
            Err(err) => panic!("failed to connect using: {}, err: {}",
                               connection_string, err)
        };

        SiteDB {
            conn,
        }
    }

    pub fn query_by_id(&self, uuid: &Uuid) -> Option<Query> {
        let rs = match self.conn.query("SELECT query_json FROM query WHERE id = $1",
                                       &[uuid])
        {
            Ok(rs) => rs,
            Err(_) => return None,
        };

        if rs.len() > 0 {
            let query_value: Option<serde_json::Value> = rs.get(0).get(0);

            let query: Query = match serde_json::value::from_value(query_value.unwrap()) {
                Ok(v) => v,
                Err(_) => return None,
            };

            return Some(query);
        }

        None
    }

    pub fn id_from_query(&self, query: &Query) -> Option<Uuid> {
        let query_value: String = match serde_json::value::to_value(query) {
            Ok(v) => v.to_string(),
            Err(err) => {
                eprint!("error converting query to string: {:?}\n", err);
                return None;
            }
        };

        let rs = match self.conn.query("SELECT id::uuid FROM query WHERE digest(query_json, 'sha256') = digest($1, 'sha256');",
                                       &[&query_value])
        {
            Ok(rs) => rs,
            Err(err) => {
                eprint!("error querying in id_from_query(): {:?}\n", err);
                return None;
            }
        };

        if rs.len() > 0 {
            let id: Uuid = rs.get(0).get("id");
            Some(id)
        } else {
            eprint!("query not in DB: {}\n", query_value);
            None
        }
    }

    pub fn save_query(&self, uuid: &Uuid, query: &Query) -> Result<(), String> {
        let trans = match self.conn.transaction() {
            Ok(t) => t,
            Err(e) => return Err(format!("couldn't begin transaction: {}", e)),
        };

        let serde_value: String = match serde_json::value::to_value(&query) {
            Ok(v) => v.to_string(),
            Err(e) => return Err(format!("serde error: {}", e)),
        };

        match trans.execute("INSERT INTO query(id, query_json) values ($1, $2)",
                            &[uuid, &serde_value]) {
            Ok(_) => (),
            Err(e) => return Err(format!("Error executing query: {}", e)),
        };

        match trans.commit() {
            Ok(_) => Ok(()),
            Err(e) => Err(format!("{}", e)),
        }
    }
}
