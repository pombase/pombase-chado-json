use std::collections::hash_map::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;
use serde_json;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Location {
    start: usize,
    end: usize,
    score: f32,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct InterProMatch {
    id: String,
    dbname: String,
    name: String,
    evidence: String,
    interpro_id: String,
    interpro_name: String,
    interpro_type: String,
    locations: Vec<Location>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct UniprotResult {
    pub interpro_matches: Vec<InterProMatch>,
}

pub fn parse_interpro(file_name: &str) -> HashMap<String, UniprotResult> {
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            print!("Failed to read {}: {}\n", file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    match serde_json::from_reader(reader) {
        Ok(uniprot_results) => uniprot_results,
        Err(err) => {
            print!("failed to parse {}: {}", file_name, err);
            process::exit(1);
        },
    }
}
