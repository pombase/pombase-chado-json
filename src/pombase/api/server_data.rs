use std::process;
use std::fs::File;
use std::io::BufReader;

use serde_json;
use std::collections::HashMap;

use web::data::{SearchAPIMaps, IdGeneSubsetMap};
use web::config::Config;
use types::GeneUniquename;


pub struct ServerData {
    config_file_name: String,
    maps: SearchAPIMaps,
    search_maps_file_name: String,
    gene_subsets: IdGeneSubsetMap,
    gene_subsets_file_name: String,
}


fn load(config: &Config, search_maps_file_name: &str, gene_subsets_file_name: &str)
        -> (SearchAPIMaps, IdGeneSubsetMap)
{
    let file = match File::open(search_maps_file_name) {
        Ok(file) => file,
        Err(err) => {
            print!("Failed to read {}: {}\n", search_maps_file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    let search_api_maps: SearchAPIMaps =
        match serde_json::from_reader(reader) {
            Ok(results) => results,
            Err(err) => {
                print!("failed to parse {}: {}", search_maps_file_name, err);
                process::exit(1);
            },
        };

    let file = match File::open(gene_subsets_file_name) {
        Ok(file) => file,
        Err(err) => {
            print!("Failed to read {}: {}\n", gene_subsets_file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    let mut gene_subsets: IdGeneSubsetMap =
        match serde_json::from_reader(reader) {
            Ok(results) => results,
            Err(err) => {
                print!("failed to parse {}: {}", gene_subsets_file_name, err);
                process::exit(1);
            },
        };

    let mut new_entries: IdGeneSubsetMap = HashMap::new();

    let prefixes_to_remove: Vec<String> =
        config.server.subsets.prefixes_to_remove
        .iter().map(|prefix| prefix.clone() + ":").collect();

    // strip prefixes and add to map
    for (subset_name, subset_details) in gene_subsets.iter() {
        for prefix in &prefixes_to_remove {
            if subset_name.starts_with(prefix) {
                let new_subset_name = &subset_name[prefix.len()..];
                new_entries.insert(String::from(new_subset_name), subset_details.clone());
            }
        }
    }

    gene_subsets.extend(new_entries);

    (search_api_maps, gene_subsets)
}

impl ServerData {
    pub fn new(config_file_name: &str, search_maps_file_name: &str,
               gene_subsets_file_name: &str)
               -> ServerData
    {
        let config = Config::read(config_file_name);

        let (maps, gene_subsets) =
            load(&config, search_maps_file_name, gene_subsets_file_name);
        ServerData {
            config_file_name: config_file_name.into(),
            search_maps_file_name: search_maps_file_name.into(),
            maps: maps,
            gene_subsets_file_name: gene_subsets_file_name.into(),
            gene_subsets: gene_subsets,
        }
    }

    pub fn genes_of_termid(&self, term_id: &str) -> Vec<GeneUniquename> {
        match self.maps.termid_genotype_genes.get(term_id) {
            Some(genotype_genes) => {
                genotype_genes.single_allele.iter().cloned().collect::<Vec<_>>()
            },
            None =>
                match self.maps.termid_genes.get(term_id) {
                    Some(gene_uniquenames) => {
                        gene_uniquenames.iter().cloned().collect::<Vec<_>>()
                    },
                    None => vec![],
                },
        }
    }

    pub fn genes_of_subset(&self, subset_name: &str) -> Vec<GeneUniquename> {
        match self.gene_subsets.get(subset_name) {
            Some(subset_details) => {
                subset_details.elements.iter().cloned().collect::<Vec<_>>()
            },
            None => vec![],
        }
    }

    pub fn reload(&mut self) {
        let config = Config::read(&self.config_file_name);
        let (maps, gene_subsets) =
            load(&config, &self.search_maps_file_name,
                 &self.gene_subsets_file_name);
        self.maps = maps;
        self.gene_subsets = gene_subsets;
    }
}
