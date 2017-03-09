// configuration for extension display names and for the "Target of" section
#[derive(Deserialize, Clone, Debug)]
pub struct ExtensionConfig {
    pub rel_name: String, // name of extension relation
    pub display_name: String, // text to display
    pub if_descendent_of: Option<String>, // None if applies to any extension
    pub reciprocal_display: Option<String>, // None if reciprocal shouldn't be displayed
}

#[derive(Deserialize, Clone, Debug)]
pub struct CvConfig {
    pub feature_type: String,
    // relations to not show in the summary
    pub summary_relations_to_hide: Vec<String>,
    // relations where the range is a gene ID to display like:
    //   has substrate pom1, cdc1 involved in negative regulation of ...
    // rather than as two lines
    pub summary_gene_relations_to_collect: Vec<String>,
}

pub type ShortEvidenceCode = String;
pub type LongEvidenceCode = String;

#[derive(Deserialize, Clone, Debug)]
pub struct Config {
    pub extensions: Vec<ExtensionConfig>,
    pub evidence_types: HashMap<ShortEvidenceCode, LongEvidenceCode>,
    pub cv_config: HashMap<CvName, CvConfig>,
}

impl Config {
    pub fn cv_config_by_name(&self, cv_name: &str) -> CvConfig {
        if let Some(config) = self.cv_config.get(cv_name) {
            config.clone()
        } else {
            CvConfig {
                feature_type: "gene".into(),
                summary_relations_to_hide: vec![],
                summary_gene_relations_to_collect: vec![],
            }
        }
    }
}
