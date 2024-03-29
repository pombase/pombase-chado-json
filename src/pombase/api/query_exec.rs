use uuid::Uuid;

use crate::api::result::*;
use crate::api::query::*;
use crate::api::site_db::SiteDB;
use crate::api_data::APIData;

use flexstr::{ToSharedStr, shared_str as flex_str, shared_fmt as flex_fmt};

pub struct QueryExec {
    api_data: APIData,
    site_db: Option<SiteDB>,
}

impl QueryExec {
    pub fn new(api_data: APIData, site_db: Option<SiteDB>) -> QueryExec {
        QueryExec {
            api_data,
            site_db,
        }
    }

    // execute a Query
    pub async fn exec(&self, query: &Query) -> QueryAPIResult {
        let query =
            if let Some(ref site_db) = self.site_db {
                if let Some(ref query_id) = query.get_constraints().get_query_id() {
                    if let Some(existing_query) = site_db.query_by_id(query_id).await {
                        existing_query
                    } else {
                        let message =
                            flex_fmt!("can't find query for ID {}", query_id);
                        return QueryAPIResult::new_error(query, message);
                    }
                } else {
                    query.clone()
                }
            } else {
                query.clone()
            };

        let uuid =
            if let Some(ref site_db) = self.site_db {
                if let Some(existing_uuid) = site_db.id_from_query(&query).await {
                    existing_uuid
                } else {
                    let new_uuid = Uuid::new_v4();
                    match site_db.save_query(&new_uuid, &query).await {
                        Ok(_) => (),
                        Err(mess) =>
                            return QueryAPIResult::new_error(&query, mess.to_shared_str()),
                    }
                    new_uuid
                }
            } else {
                Uuid::new_v4()
            };

        let id = uuid.as_hyphenated().to_string().to_shared_str();

        let rows_result = query.exec(&self.api_data, &self.site_db).await;

        match rows_result {
            Ok(rows) => {
                QueryAPIResult {
                    query,
                    id,
                    status: flex_str!("ok"),
                    rows,
                }
            },
            Err(mess) => QueryAPIResult::new_error(&query, mess),
        }
    }

    pub fn get_api_data(&self) -> &APIData {
        &self.api_data
    }
}
