use api::query::Query;

#[derive(Serialize, Deserialize, Debug)]
pub struct ResultRow {
    pub gene_uniquename: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub sequence: Option<String>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct QueryAPIResult {
    pub query: Query,
    pub status: String,
    pub rows: Vec<ResultRow>,
}
