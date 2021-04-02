use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;
use serde_json;

use pombase_rc_string::RcString;

#[derive(Debug, Deserialize)]
pub struct PfamMotifMetadata {
    pub database: RcString,
}

#[derive(Debug, Deserialize)]
pub struct PfamMotif {
    #[serde(rename = "type")]
    pub motif_type: RcString,
    pub metadata: PfamMotifMetadata,
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Deserialize)]
pub struct PfamProteinDetails {
    pub motifs: Vec<PfamMotif>,
}

pub fn parse_pfam(file_name: &str) -> HashMap<RcString, PfamProteinDetails> {
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
