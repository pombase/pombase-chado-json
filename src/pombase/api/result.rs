use std::collections::HashSet;

use crate::api::query::Query;
use crate::data_types::{DeletionViability, GeneQueryAttrName, GeneQueryTermData, GoCamId, PresentAbsent};
use crate::types::{TermId, GeneUniquename, ReferenceUniquename, PdbId};

use flexstr::{SharedStr as FlexStr, shared_str as flex_str};

#[derive(Serialize, Deserialize, Debug, PartialEq)]
pub struct GeneExValue {
    pub dataset_name: FlexStr,
    pub value: FlexStr,
}

#[derive(Serialize, Deserialize, Debug, PartialEq)]
pub struct ResultRow {
    pub gene_uniquename: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub deletion_viability: Option<DeletionViability>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub go_component: Option<GeneQueryTermData>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub go_process_superslim: Option<GeneQueryTermData>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub go_function: Option<GeneQueryTermData>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub characterisation_status: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub taxonomic_distribution: Option<FlexStr>,
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
    pub spliced_rna_length: Option<usize>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub unspliced_rna_length: Option<usize>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub sequence: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gaf_lines: Option<String>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub reference_uniquenames: HashSet<ReferenceUniquename>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub pdb_ids: HashSet<PdbId>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub gocam_ids: HashSet<GoCamId>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub subsets: HashSet<TermId>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gene_expression: Vec<GeneExValue>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct QueryAPIResult {
    pub query: Query,
    pub id: FlexStr,
    pub status: FlexStr,
    pub rows: Vec<ResultRow>,
}

impl QueryAPIResult {
    pub fn new_error(query: &Query, error_message: FlexStr) -> QueryAPIResult {
        QueryAPIResult {
            query: query.clone(),
            id: flex_str!("error"),
            status: error_message,
            rows: vec![],
        }
    }
}
