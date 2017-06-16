use std::process;
use std::fs::File;
use std::io::BufReader;

use serde_json;

use web::data::SearchAPIMaps;
use types::GeneUniquename;


pub struct ServerData {
    search_api_maps: SearchAPIMaps,
    file_name: String,
}

fn load(file_name: &str) -> SearchAPIMaps {
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            print!("Failed to read {}: {}\n", file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    let search_api_maps: SearchAPIMaps =
        match serde_json::from_reader(reader) {
            Ok(results) => results,
            Err(err) => {
                print!("failed to parse {}: {}", file_name, err);
                process::exit(1);
            },
        };

    search_api_maps
}


impl ServerData {
    pub fn new(file_name: &str) -> ServerData {
        ServerData {
            file_name: file_name.into(),
            search_api_maps: load(file_name),
        }
    }

    pub fn genes_of_termid(&self, term_id: &str) -> Vec<GeneUniquename> {
        match self.search_api_maps.termid_genes.get(term_id) {
            Some(gene_uniquenames) => {
                gene_uniquenames.iter().cloned().collect::<Vec<_>>()
            },
            None => vec![],
        }
    }

    pub fn reload(&mut self) {
        self.search_api_maps = load(&self.file_name);
    }
}

