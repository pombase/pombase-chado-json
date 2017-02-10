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
// relationships where the more specific term is the object, used when reading cvtermpath
pub const INVERSE_REL_NAMES: [&'static str; 1] = ["has_part"];
// only consider inverse relations for these ontologies:
pub const INVERSE_REL_CV_NAMES: [&'static str; 1] = ["fission_yeast_phenotype"];

// number of genes before (and after) to add to the gene_neighbourhood field
pub const GENE_NEIGHBOURHOOD_DISTANCE: usize = 5;

pub const TRANSCRIPT_FEATURE_TYPES: [&'static str; 6] =
    ["snRNA", "rRNA", "mRNA", "snoRNA", "ncRNA", "tRNA"];

// "interesting parents" are those stored in the JSON in the TermShort structs
#[derive(Clone)]
pub struct InterestingParent {
    pub termid: &'static str,
    pub rel_name: &'static str,
}

// when creating a TermShort struct, for each of these termids if the term has
// an "interesting parent" using the given rel_name, we store it in the
// interesting_parents field of the TermShort
pub const INTERESTING_PARENTS: [InterestingParent; 4] = [
    InterestingParent {
        termid: "FYPO:0000002",
        rel_name: "is_a",
    },
    InterestingParent {
        termid: "FYPO:0000003",
        rel_name: "is_a",
    },
    InterestingParent {
        termid: "FYPO:0000652",
        rel_name: "is_a",
    },

    // a with property on a protein binding (GO:0005515) is displayed
    // as a binds extension - configure this to allow easy checking for that in
    // process_feature_cvterms()
    // see: https://github.com/pombase/website/issues/108
    InterestingParent {
        termid: "GO:0005515",
        rel_name: "is_a",
    },
];
