use hashbrown::HashSet;

use crate::api::query::Query;

use crate::web::data::{DeletionViability, PresentAbsent, GeneQueryTermData};

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
    pub tmm: Option<PresentAbsent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub sequence: Option<RcString>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct QueryAPIResult {
    pub query: Query,
    pub status: RcString,
    pub rows: Vec<ResultRow>,
}
