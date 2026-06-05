pub const API_MAPS_SQLITE3_FILE_NAME: &str = "api_maps.sqlite3";

pub const API_MAPS_TABLE_NAMES: &[&str; 7] =
    &["terms", "genes", "alleles", "refs", "genotypes",
      "annotation_detail", "termid_genotype_annotations"];

pub const FYPO_ROOT_TERM_ID: &str = "FYPO:0000001";

pub const MOLECULAR_FUNCTION_ROOT_TERM_ID: &str = "GO:0003674";
pub const CELLULAR_COMPONENT_ROOT_TERM_ID: &str = "GO:0005575";
pub const BIOLOGICAL_PROCESS_ROOT_TERM_ID: &str = "GO:0008150";
