use api::query::Query;

#[derive(Serialize, Deserialize, Debug)]
pub struct ResultRow {
    pub gene_uniquename: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct QueryAPIResult {
    pub query: Query,
    pub status: String,
    pub rows: Vec<ResultRow>,
}
