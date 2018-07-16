use api::query::Query;

use web::data::{DeletionViability, GeneQueryTermData};

#[derive(Serialize, Deserialize, Debug, PartialEq)]
pub struct ResultRow {
    pub gene_uniquename: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub deletion_viability: Option<DeletionViability>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub go_component: Option<GeneQueryTermData>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub sequence: Option<String>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct QueryAPIResult {
    pub query: Query,
    pub status: String,
    pub rows: Vec<ResultRow>,
}
