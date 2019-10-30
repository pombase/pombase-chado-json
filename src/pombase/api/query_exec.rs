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

    // execute a Query
    pub fn exec(&self, query: &Query) -> QueryAPIResult {
        let query =
            if let Some(ref site_db) = self.site_db {
                if let QueryNode::QueryId { id } = query.get_constraints() {
                    if let Some(existing_query) = site_db.query_by_id(id) {
                        existing_query
                    } else {
                        let message = RcString::from(&format!("can't find query for ID {}", id));
                        return QueryAPIResult::new_error(&query, &message);

                    }
                } else {
                    query.clone()
                }
            } else {
                query.clone()
            };

        let uuid =
            if let Some(ref site_db) = self.site_db {
                if let Some(existing_uuid) = site_db.id_from_query(&query) {
                    existing_uuid
                } else {
                    let new_uuid = Uuid::new_v4();
                    match site_db.save_query(&new_uuid, &query) {
                        Ok(_) => (),
                        Err(mess) =>
                            return QueryAPIResult::new_error(&query, &mess),
                    }
                    new_uuid
                }
            } else {
                Uuid::new_v4()
            };

        let id = RcString::from(&uuid.hyphenated().to_string());

        let rows_result = query.exec(&self.server_data, &self.site_db);

        match rows_result {
            Ok(rows) => {
                QueryAPIResult {
                    query,
                    id,
                    status: RcString::from("ok"),
                    rows,
                }
            },
            Err(mess) => QueryAPIResult::new_error(&query, &mess),
        }
    }

    pub fn get_server_data(&self) -> &ServerData {
        &self.server_data
    }

    pub fn reload(&mut self) {
        self.server_data.reload();
    }
}
