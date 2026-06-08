pub const API_MAPS_SQLITE3_FILE_NAME: &str = "api_maps.sqlite3";

pub const API_MAPS_TABLE_NAMES: &[&str; 7] =
    &["terms", "genes", "alleles", "refs", "genotypes",
      "annotation_detail", "termid_genotype_annotations"];

pub const FYPO_ROOT_TERM_ID: &str = "FYPO:0000001";
pub const FYPO_CV_NAME: &str = "fission_yeast_phenotype";

pub const MOLECULAR_FUNCTION_ROOT_TERM_ID: &str = "GO:0003674";
pub const MOLECULAR_FUNCTION_ROOT_TERM_NAME: &str = "molecular_function";
pub const CELLULAR_COMPONENT_ROOT_TERM_ID: &str = "GO:0005575";
pub const CELLULAR_COMPONENT_ROOT_TERM_NAME: &str = "cellular_component";
pub const BIOLOGICAL_PROCESS_ROOT_TERM_ID: &str = "GO:0008150";
pub const BIOLOGICAL_PROCESS_ROOT_TERM_NAME: &str = "biological_process";

pub fn is_go_root_name(term_name: &str) -> bool {
    term_name == MOLECULAR_FUNCTION_ROOT_TERM_NAME ||
        term_name == CELLULAR_COMPONENT_ROOT_TERM_NAME ||
        term_name == BIOLOGICAL_PROCESS_ROOT_TERM_NAME
}
