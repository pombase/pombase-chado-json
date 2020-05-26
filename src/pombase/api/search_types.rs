use std::collections::HashMap;

pub type SolrMatchId = String;
pub type SolrFieldName = String;
pub type SolrHighlightDetails = String;

pub type SolrMatchHighlight = HashMap<SolrFieldName, Vec<SolrHighlightDetails>>;
