use hashbrown::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;
use serde_json;

use crate::web::config::Config;

use pombase_rc_string::RcString;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Location {
    pub start: usize,
    pub end: usize,
    pub score: f32,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct InterProMatch {
    pub id: RcString,
    pub dbname: RcString,
    pub name: RcString,
    pub model: Option<RcString>,
    pub evidence: RcString,
    pub interpro_id: RcString,
    pub interpro_name: RcString,
    pub interpro_type: RcString,
    pub locations: Vec<Location>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct TMMatch {
    pub start: usize,
    pub end: usize,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct UniprotResult {
    pub interpro_matches: Vec<InterProMatch>,
    pub tmhmm_matches: Vec<TMMatch>,
}

pub fn parse_interpro(config: &Config, file_name: &str) -> HashMap<RcString, UniprotResult> {
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            print!("Failed to read {}: {}\n", file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    let uniprot_results: HashMap<RcString, UniprotResult> =
        match serde_json::from_reader(reader) {
            Ok(uniprot_results) => uniprot_results,
            Err(err) => {
                print!("failed to parse {}: {}", file_name, err);
                process::exit(1);
            },
        };

    let mut filtered_results = HashMap::new();

    for (uniprot_identifier, mut results) in uniprot_results {
        let new_interpro_matches =
            results.interpro_matches.into_iter()
            .filter(|interpro_match|
                    !config.interpro.dbnames_to_filter.contains(&interpro_match.dbname))
            .collect();

        results.interpro_matches = new_interpro_matches;
        filtered_results.insert(uniprot_identifier, results);
    }

    filtered_results
}
