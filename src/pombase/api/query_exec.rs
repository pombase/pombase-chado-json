use api::result::*;
use api::query::*;
use api::server_data::ServerData;



pub struct QueryExec {
    server_data: ServerData,
}

impl QueryExec {
    pub fn new(server_data: ServerData) -> QueryExec {
        QueryExec {
            server_data: server_data,
        }
    }

    // add missing term definitions
    fn fill_query_details(&self, query: &Query) -> Query {
        query.clone()
    }

    pub fn exec(&self, query: &Query) -> QueryAPIResult {
        let filled_query = self.fill_query_details(query);
        let rows_result = query.exec(&self.server_data);

        match rows_result {
            Ok(rows) => QueryAPIResult {
                query: filled_query,
                status: "ok".to_owned(),
                rows: rows,
            },
            Err(mess) => QueryAPIResult {
                query: filled_query,
                status: mess,
                rows: vec![],
            },
        }
    }

    pub fn reload(&mut self) {
        self.server_data.reload();
    }
}
