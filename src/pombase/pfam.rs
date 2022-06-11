use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;
use serde_json;

use flexstr::SharedStr as FlexStr;

#[derive(Debug, Deserialize)]
pub struct PfamMotifMetadata {
    pub database: FlexStr,
}

#[derive(Debug, Deserialize)]
pub struct PfamMotif {
    #[serde(rename = "type")]
    pub motif_type: FlexStr,
    pub metadata: PfamMotifMetadata,
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Deserialize)]
pub struct PfamProteinDetails {
    pub motifs: Vec<PfamMotif>,
}

pub fn parse_pfam(file_name: &str) -> HashMap<FlexStr, PfamProteinDetails> {
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            println!("Failed to read {}: {}", file_name, err);
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
