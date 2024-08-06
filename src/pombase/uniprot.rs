use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;

use flexstr::{SharedStr as FlexStr, shared_str as flex_str};
use regex::{Captures, Regex};

use crate::data_types::{ActiveSite, BetaStrand, BindingSite, Chain,
                        Helix, PeptideRange, Propeptide, SignalPeptide,
                        TransitPeptide, Turn};
use crate::types::GeneUniquename;

#[derive(Debug, Serialize, Deserialize)]
pub struct UniProtDataEntry {
    pub gene_uniquename: GeneUniquename,
    pub sequence: FlexStr,
    pub signal_peptide: Option<SignalPeptide>,
    pub transit_peptide: Option<TransitPeptide>,
    pub binding_sites: Vec<BindingSite>,
    pub active_sites: Vec<ActiveSite>,
    pub beta_strands: Vec<BetaStrand>,
    pub helices: Vec<Helix>,
    pub turns: Vec<Turn>,
    pub propeptides: Vec<Propeptide>,
    pub chains: Vec<Chain>,
}
pub type UniProtDataMap = HashMap<GeneUniquename, UniProtDataEntry>;

// columns names come from the UniProt advanced search download
#[derive(Debug, Deserialize)]
struct UniProtDataRecord {
//    #[serde(rename = "Entry")]
//    uniprot_accession: String,
    #[serde(rename = "PomBase")]
    gene_uniquename: String,
    #[serde(rename = "Signal peptide")]
    signal_peptide: String,
    #[serde(rename = "Transit peptide")]
    transit_peptide: String,
    #[serde(rename = "Binding site")]
    binding_sites: String,
    #[serde(rename = "Active site")]
    active_sites: String,
    #[serde(rename = "Beta strand")]
    beta_strands: String,
    #[serde(rename = "Helix")]
    helices: String,
    #[serde(rename = "Turn")]
    turns: String,
    #[serde(rename = "Chain")]
    chains: String,
    #[serde(rename = "Propeptide")]
    propeptides: String,

    #[serde(rename = "Sequence")]
    sequence: String,
/*
 Catalytic activity
 Gene Names (synonym)
 Post-translational modification
 Modified residue
 Cofactor
 Kinetics
 */
}

lazy_static! {
    static ref SPLIT_RE: Regex = Regex::new(r"(SIGNAL|TRANSIT|BINDING|ACT_SITE|STRAND|HELIX|TURN|PROPEP|CHAIN) ").unwrap();
    static ref RANGE_RE: Regex = Regex::new(r"^(\?|\d+)(?:\.\.([\?\d]+))?.*?(;.*)?").unwrap();
    static ref LIGAND_RE: Regex = Regex::new(r#"/ligand="([^"]+)""#).unwrap();
}

fn get_range(cap: Captures) -> Option<PeptideRange> {
    let start: usize = cap.get(1).unwrap().as_str().parse().ok()?;
    if let Some(cap2) = cap.get(2) {
        let end: usize = cap2.as_str().parse().ok()?;
        Some(PeptideRange {
            start,
            end,
        })
    } else {
        // single base like:
        // BINDING 158; /ligand="(1,3-beta-D-glucosyl)n"
        Some(PeptideRange {
            start,
            end: start,
        })
    }
}

fn get_ligand(rest: &str) -> FlexStr {
    if let Some(ligand_capture) = LIGAND_RE.captures_iter(rest).next() {
        if let Some(ligand_match) = ligand_capture.get(1) {
            return ligand_match.as_str().into()
        }
    }

    flex_str!("unknown ligand")
}

fn first_field_part(field: &str) -> Option<String> {
    let mut parts_iter = SPLIT_RE.split(field);
    parts_iter.next();
    parts_iter.next().map(|s| s.to_owned())
}

fn get_propeptides(uniprot_record: &UniProtDataRecord) -> Vec<Propeptide> {
    let mut propeptides_parts_iter = SPLIT_RE.split(&uniprot_record.propeptides);
    propeptides_parts_iter.next();  // remove blank

    propeptides_parts_iter.filter_map(|field_part| {
        let cap = RANGE_RE.captures_iter(field_part).next()?;
        let range = get_range(cap)?;
        Some(Propeptide {
            range,
        })
    })
    .collect()
}

fn get_chains(uniprot_record: &UniProtDataRecord) -> Vec<Chain> {
    let mut chains_parts_iter = SPLIT_RE.split(&uniprot_record.chains);
    chains_parts_iter.next();  // remove blank

    chains_parts_iter.filter_map(|field_part| {
        let cap = RANGE_RE.captures_iter(field_part).next()?;
        let range = get_range(cap)?;
        Some(Chain {
                range,
        })
    })
    .collect()
}

fn process_record(uniprot_record: UniProtDataRecord) -> UniProtDataEntry {
    let gene_uniquename =
        if uniprot_record.gene_uniquename.ends_with(";") {
            let mut uniquename = uniprot_record.gene_uniquename.to_string();
            uniquename.pop();
            uniquename
        } else {
            uniprot_record.gene_uniquename.to_string()
        };

    let mut signal_peptide = None;
    if let Some(field_part) = first_field_part(&uniprot_record.signal_peptide) {
        if let Some(cap) = RANGE_RE.captures_iter(&field_part).next() {
            if let Some(range) = get_range(cap) {
                signal_peptide = Some(SignalPeptide {
                    range,
                })
            }
        }
    }

    let mut transit_peptide = None;
    if let Some(field_part) = first_field_part(&uniprot_record.transit_peptide) {
        if let Some(cap) = RANGE_RE.captures_iter(&field_part).next() {
            if let Some(range) = get_range(cap) {
                transit_peptide = Some(TransitPeptide {
                    range,
                });
            }
        }
    }

    let mut binding_sites_parts_iter =
        SPLIT_RE.split(&uniprot_record.binding_sites);
    binding_sites_parts_iter.next();  // remove blank

    let binding_sites =
        binding_sites_parts_iter.filter_map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                let ligand =
                    if let Some(rest) = cap.get(3) {
                        get_ligand(rest.as_str())
                    } else {
                        flex_str!("unknown ligand")
                    };
                let range = get_range(cap)?;
                Some(BindingSite {
                    ligand,
                    range,
                })
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    let mut active_sites_parts_iter =
        SPLIT_RE.split(&uniprot_record.active_sites);
    active_sites_parts_iter.next();  // remove blank

    let active_sites =
        active_sites_parts_iter.filter_map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                let range = get_range(cap)?;
                Some(ActiveSite {
                    range,
                })
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    let mut beta_strands_parts_iter =
        SPLIT_RE.split(&uniprot_record.beta_strands);
    beta_strands_parts_iter.next();  // remove blank

    let beta_strands =
        beta_strands_parts_iter.filter_map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                let range = get_range(cap)?;
                Some(BetaStrand {
                    range,
                })
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    let mut helices_parts_iter = SPLIT_RE.split(&uniprot_record.helices);
    helices_parts_iter.next();  // remove blank

    let helices =
        helices_parts_iter.filter_map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                let range = get_range(cap)?;
                Some(Helix {
                    range,
                })
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    let mut turns_parts_iter = SPLIT_RE.split(&uniprot_record.turns);
    turns_parts_iter.next();  // remove blank

    let turns =
        turns_parts_iter.filter_map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                let range = get_range(cap)?;
                Some(Turn {
                    range,
                })
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    let propeptides = get_propeptides(&uniprot_record);
    let chains = get_chains(&uniprot_record);

    UniProtDataEntry {
        gene_uniquename: gene_uniquename.into(),
        sequence: uniprot_record.sequence.into(),
        signal_peptide,
        transit_peptide,
        binding_sites,
        active_sites,
        beta_strands,
        helices,
        turns,
        propeptides,
        chains,
    }
}

pub fn parse_uniprot(file_name: &str) -> UniProtDataMap {
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            println!("Failed to open {}: {}", file_name, err);
            process::exit(1);
        }
    };

    let reader = BufReader::new(file);

    let mut csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(reader);

    let mut res: UniProtDataMap = HashMap::new();

    for result in csv_reader.deserialize() {
        let record: UniProtDataRecord =
            result.unwrap_or_else(|e| {
                panic!("failed to read gene history CSV file: {}", e);
            });

        let entry = process_record(record);

        res.insert(entry.gene_uniquename.clone(), entry);
    }

    res
}

