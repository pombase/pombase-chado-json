use uuid::Uuid;

use crate::api::result::*;
use crate::api::query::*;
use crate::api::server_data::ServerData;

use pombase_rc_string::RcString;

pub struct QueryExec {
    server_data: ServerData,
}

impl QueryExec {
    pub fn new(server_data: ServerData) -> QueryExec {
        QueryExec {
            server_data,
        }
    }

    // add missing term definitions
    fn fill_query_details(&self, query: &Query) -> Query {
        query.clone()
    }

    // execute a Query
    pub fn exec(&self, query: &Query) -> QueryAPIResult {
        let filled_query = self.fill_query_details(query);
        let rows_result = query.exec(&self.server_data);

        match rows_result {
            Ok(rows) => QueryAPIResult {
                query: filled_query,
                id: RcString::from(&Uuid::new_v4().to_hyphenated().to_string()),
                status: RcString::from("ok"),
                rows,
            },
            Err(mess) => QueryAPIResult {
                query: filled_query,
                id: RcString::from("error"),
                status: mess,
                rows: vec![],
            },
        }
    }

    pub fn get_server_data(&self) -> &ServerData {
        &self.server_data
    }

    pub fn reload(&mut self) {
        self.server_data.reload();
    }
}
