#[derive(Serialize, Deserialize)]
pub struct ResultRow {
    pub gene_uniquename: String,
}

#[derive(Serialize, Deserialize)]
pub enum ResultStatus {
    OK,
    Error(String),
}

#[derive(Serialize, Deserialize)]
pub struct Result {
    pub status: ResultStatus,
    pub rows: Vec<ResultRow>,
}
