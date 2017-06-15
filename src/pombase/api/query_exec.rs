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

    pub fn exec(&self, query: &Query) -> Result {
        Result {
            status: ResultStatus::OK,
            rows: vec![ResultRow{ gene_uniquename: "SPAC27E2.05".into(),
            }]
        }
    }

    pub fn reload(&mut self) {
        self.server_data.reload();
    }
}
