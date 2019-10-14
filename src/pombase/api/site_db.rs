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

    pub fn save_query(&self, uuid: &Uuid, query: &Query) -> Result<(), String> {
        let trans = match self.conn.transaction() {
            Ok(t) => t,
            Err(e) => return Err(format!("{}", e)),
        };

        let serde_value = match serde_json::value::to_value(&query) {
            Ok(v) => v,
            Err(e) => return Err(format!("{}", e)),
        };

        match trans.execute("INSERT INTO query(id, query_json) values ($1, $2)",
                            &[uuid, &serde_value]) {
            Ok(_) => (),
            Err(e) => return Err(format!("{}", e)),
        };

        match trans.commit() {
            Ok(_) => Ok(()),
            Err(e) => Err(format!("{}", e)),
        }
    }
}
