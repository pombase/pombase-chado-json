use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;

use crate::web::config::Config;

use flexstr::SharedStr as FlexStr;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Location {
    pub start: usize,
    pub end: usize,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct InterProMatch {
    pub id: FlexStr,
    pub dbname: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub description: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub interpro_id: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub interpro_name: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub interpro_description: Option<FlexStr>,
    pub match_start: usize,
    pub match_end: usize,
    pub locations: Vec<Location>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct TMMatch {
    pub start: usize,
    pub end: usize,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct GeneMatches {
    pub interpro_matches: Vec<InterProMatch>,
    pub segmasker_matches: Vec<Location>,
    pub tmhmm_matches: Vec<TMMatch>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct DomainData {
    pub interproscan_version: FlexStr,
    pub domains_by_id: HashMap<FlexStr, GeneMatches>,
}

fn sort_matches(config: &Config, matches: &mut [InterProMatch]) {
    let db_order = &config.protein_feature_db_order;

    let order_pos = move |m: &InterProMatch| {
        #[allow(clippy::needless_range_loop)]
        for pos in 0..db_order.len() {
            if m.dbname.as_str().starts_with(db_order[pos].as_str()) {
                return pos;
            }
        }

        db_order.len()
    };


    let sort_fn = |a: &InterProMatch, b: &InterProMatch| {
        let a_pos = order_pos(a);
        let b_pos = order_pos(b);

        a_pos.cmp(&b_pos)
    };

    matches.sort_by(sort_fn);
}

pub fn parse_interpro(config: &Config, file_name: &str) -> DomainData {
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            println!("Failed to read {}: {}", file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    let domain_data: DomainData =
        match serde_json::from_reader(reader) {
            Ok(uniprot_results) => uniprot_results,
            Err(err) => {
                print!("failed to parse {}: {}", file_name, err);
                process::exit(1);
            },
        };

    let mut filtered_domains = HashMap::new();

    for (gene_uniquename, mut results) in domain_data.domains_by_id {
        let mut new_interpro_matches: Vec<InterProMatch> =
            results.interpro_matches.into_iter()
            .filter(|interpro_match|
                    !config.interpro.dbnames_to_filter.contains(&interpro_match.dbname))
            .collect();

        sort_matches(config, &mut new_interpro_matches);

        results.interpro_matches = new_interpro_matches;
        filtered_domains.insert(gene_uniquename, results);
    }

    DomainData {
        interproscan_version: domain_data.interproscan_version,
        domains_by_id: filtered_domains,
    }
}
