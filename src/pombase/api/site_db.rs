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
        let query_value = match serde_json::value::to_value(query) {
            Ok(v) => v,
            Err(_) => return None,
        };

        let rs = match self.conn.query("SELECT id FROM query WHERE query_json = $1",
                                       &[&query_value])
        {
            Ok(rs) => rs,
            Err(_) => return None,
        };

        if rs.len() > 0 {
            Some(rs.get(0).get(0))
        } else {
            None
        }
    }

    pub fn save_query(&self, uuid: &Uuid, query: &Query) -> Result<(), String> {
        let trans = match self.conn.transaction() {
            Ok(t) => t,
            Err(e) => return Err(format!("couldn't begin transaction: {}", e)),
        };

        let serde_value = match serde_json::value::to_value(&query) {
            Ok(v) => v,
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
