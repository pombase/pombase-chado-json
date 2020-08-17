use std::collections::HashSet;

use crate::api::query::Query;
use crate::data_types::{DeletionViability, PresentAbsent, GeneQueryTermData, GeneQueryAttrName};
use crate::types::{TermId, GeneUniquename};

use pombase_rc_string::RcString;

#[derive(Serialize, Deserialize, Debug, PartialEq)]
pub struct ResultRow {
    pub gene_uniquename: RcString,
    #[serde(skip_serializing_if="Option::is_none")]
    pub deletion_viability: Option<DeletionViability>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub go_component: Option<GeneQueryTermData>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub go_process_superslim: Option<GeneQueryTermData>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub go_function: Option<GeneQueryTermData>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub characterisation_status: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub taxonomic_distribution: Option<RcString>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub ortholog_taxonids: HashSet<u32>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub physical_interactors: HashSet<GeneUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub tmm: Option<PresentAbsent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub molecular_weight: Option<f32>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub protein_length: Option<usize>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub protein_length_bin: Option<GeneQueryAttrName>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub sequence: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gaf_lines: Option<String>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub subsets: HashSet<TermId>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct QueryAPIResult {
    pub query: Query,
    pub id: RcString,
    pub status: RcString,
    pub rows: Vec<ResultRow>,
}

impl QueryAPIResult {
    pub fn new_error(query: &Query, error_message: &str) -> QueryAPIResult {
        QueryAPIResult {
            query: query.clone(),
            id: RcString::from("error"),
            status: RcString::from(error_message),
            rows: vec![],
        }
    }
}
