#[derive(Serialize, Deserialize)]
pub enum QueryPart {
#[serde(rename = "or")]
    Or(Vec<QueryPart>),
#[serde(rename = "and")]
    And(Vec<QueryPart>),
#[serde(rename = "termid")]
    TermId(String),
}


#[derive(Serialize, Deserialize)]
pub struct Query {
    pub query: QueryPart,
}
