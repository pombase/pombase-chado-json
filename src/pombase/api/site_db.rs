use postgres::{Connection, TlsMode};

use crate::web::config::SiteDBConfig;

pub struct SiteDB {
    conn: Connection,
}

impl SiteDB {
    pub fn new(conf: &SiteDBConfig) -> SiteDB {
        let conn_string = conf.connection_string.as_str();
        let conn = match Connection::connect(conn_string, TlsMode::None) {
            Ok(conn) => conn,
            Err(err) => panic!("failed to connect using: {}, err: {}", conn_string, err)
        };

        SiteDB {
            conn,
        }
    }

    pub fn connection(&self) -> &Connection {
        &self.conn
    }
}
