// configuration for extension display names and for the "Target of" section
#[derive(Deserialize, Clone, Debug)]
pub struct ExtensionDisplayNames {
    pub rel_name: String, // name of extension relation
    pub display_name: String, // text to display
    pub if_descendent_of: Option<String>, // None if applies to any extension
    pub reciprocal_display: Option<String>, // None if reciprocal shouldn't be displayed
}

// "interesting parents" are those stored in the JSON in the TermShort structs
#[derive(Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct InterestingParent {
    pub termid: String,
    pub rel_name: String,
}

// the order of relations within an extension:
#[derive(Deserialize, Clone, Debug)]
pub struct RelationOrder {
    // put the relations in this order in the displayed extensions:
    pub relation_order: Vec<String>,
    // except for these reactions which should always come last:
    pub always_last: Vec<String>,
}


#[derive(Deserialize, Clone, Debug)]
pub struct TermFilterCategory {
    display_name: String,
    // this category matches these terms and their descendants
    ancestors: Vec<TermId>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct FilterConfig {
    filter_name: String,
    display_name: String,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    term_categories: Vec<TermFilterCategory>
}

#[derive(Deserialize, Clone, Debug)]
pub struct CvConfig {
    pub feature_type: String,
    // filtering configured per CV
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub filters: Vec<FilterConfig>,
    // relations to not show in the summary
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub summary_relations_to_hide: Vec<String>,
    // relations where the range is a gene ID to display like:
    //   has substrate pom1, cdc1 involved in negative regulation of ...
    // rather than as two lines
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub summary_gene_relations_to_collect: Vec<String>,
}

pub type ShortEvidenceCode = String;
pub type LongEvidenceCode = String;

#[derive(Deserialize, Clone, Debug)]
pub struct ConfigOrganism {
    pub genus: String,
    pub species: String,
}

impl ConfigOrganism {
    pub fn full_name(&self) -> String {
        self.genus.clone() + "_" + self.species.as_str()
    }
}

#[derive(Deserialize, Clone, Debug)]
pub struct Config {
    pub load_organism: ConfigOrganism,
    pub extension_display_names: Vec<ExtensionDisplayNames>,
    pub extension_relation_order: RelationOrder,
    pub evidence_types: HashMap<ShortEvidenceCode, LongEvidenceCode>,
    pub cv_config: HashMap<CvName, CvConfig>,
// when creating a TermShort struct, for each of these termids if the term has
// an "interesting parent" using the given rel_name, we store it in the
// interesting_parents field of the TermShort
    pub interesting_parents: Vec<InterestingParent>,
}

impl Config {
    pub fn cv_config_by_name(&self, cv_name: &str) -> CvConfig {
        if let Some(config) = self.cv_config.get(cv_name) {
            config.clone()
        } else {
            CvConfig {
                feature_type: "gene".into(),
                filters: vec![],
                summary_relations_to_hide: vec![],
                summary_gene_relations_to_collect: vec![],
            }
        }
    }
}
