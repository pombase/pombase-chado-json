pub const POMBASE_ANN_EXT_TERM_CV_NAME: &'static str = "PomBase annotation extension terms";
pub const ANNOTATION_EXT_REL_PREFIX: &'static str = "annotation_extension_relation-";

pub const DB_NAME: &'static str = "PomBase";

pub enum FeatureRelAnnotationType {
    Interaction,
    Ortholog,
    Paralog,
}
pub struct FeatureRelConfig {
    pub rel_type_name: &'static str,
    pub annotation_type: FeatureRelAnnotationType,
}
pub const FEATURE_REL_CONFIGS: [FeatureRelConfig; 4] =
    [
        FeatureRelConfig {
            rel_type_name: "interacts_physically",
            annotation_type: FeatureRelAnnotationType::Interaction,
        },
        FeatureRelConfig {
            rel_type_name: "interacts_genetically",
            annotation_type: FeatureRelAnnotationType::Interaction,
        },
        FeatureRelConfig {
            rel_type_name: "orthologous_to",
            annotation_type: FeatureRelAnnotationType::Ortholog,
        },
        FeatureRelConfig {
            rel_type_name: "paralogous_to",
            annotation_type: FeatureRelAnnotationType::Paralog,
        },
    ];

// relations to use when copy annotation to parents (ie. adding the
// annotation of child terms to parents)
pub const DESCENDANT_REL_NAMES: [&'static str; 3] = ["is_a", "part_of", "regulates"];

// number of genes before (and after) to add to the gene_neighbourhood field
pub const GENE_NEIGHBOURHOOD_DISTANCE: usize = 5;

