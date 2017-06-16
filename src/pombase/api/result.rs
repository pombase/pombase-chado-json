#[derive(Serialize, Deserialize, Debug)]
pub struct ResultRow {
    pub gene_uniquename: String,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug)]
pub enum ResultStatus {
    Ok,
    Error(String),
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Result {
    pub status: ResultStatus,
    pub rows: Vec<ResultRow>,
}
