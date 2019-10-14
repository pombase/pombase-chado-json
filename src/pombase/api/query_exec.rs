use uuid::Uuid;

use crate::api::result::*;
use crate::api::query::*;
use crate::api::site_db::SiteDB;
use crate::api::server_data::ServerData;

use pombase_rc_string::RcString;

pub struct QueryExec {
    server_data: ServerData,
    site_db: Option<SiteDB>,
}

impl QueryExec {
    pub fn new(server_data: ServerData, site_db: Option<SiteDB>) -> QueryExec {
        QueryExec {
            server_data,
            site_db,
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
            Ok(rows) => {
                let uuid = &Uuid::new_v4();
                let id = RcString::from(&uuid.hyphenated().to_string());

                if let Some(site_db) = &self.site_db {
                    match site_db.save_query(uuid, &filled_query) {
                        Ok(_) => (),
                        Err(mess) =>
                            return QueryAPIResult::new_error(&filled_query, &mess),
                    }
                }

                QueryAPIResult {
                    query: filled_query,
                    id,
                    status: RcString::from("ok"),
                    rows,
                }
            },
            Err(mess) => QueryAPIResult::new_error(&filled_query, &mess),
        }
    }

    pub fn get_server_data(&self) -> &ServerData {
        &self.server_data
    }

    pub fn reload(&mut self) {
        self.server_data.reload();
    }
}
