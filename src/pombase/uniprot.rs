use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;

use flexstr::{SharedStr as FlexStr, shared_str as flex_str};
use regex::{Captures, Regex};

use crate::data_types::{ActiveSite, BetaStrand, BindingSite, Helix,
                        PeptideRange, SignalPeptide, TransitPeptide};
use crate::types::GeneUniquename;

#[derive(Debug, Serialize, Deserialize)]
pub struct UniProtDataEntry {
    pub gene_uniquename: GeneUniquename,
    pub signal_peptide: Option<SignalPeptide>,
    pub transit_peptide: Option<TransitPeptide>,
    pub binding_sites: Vec<BindingSite>,
    pub active_sites: Vec<ActiveSite>,
    pub beta_strands: Vec<BetaStrand>,
    pub helices: Vec<Helix>,
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
    static ref SPLIT_RE: Regex = Regex::new(r"(SIGNAL|TRANSIT|BINDING|ACT_SITE|STRAND|HELIX) ").unwrap();
    static ref RANGE_RE: Regex = Regex::new(r"^(\d+)(?:\.\.([\?\d]+))?.*?(;.*)?").unwrap();
    static ref LIGAND_RE: Regex = Regex::new(r#"/ligand="([^"]+)""#).unwrap();
}

fn get_range(cap: Captures) -> PeptideRange {
  let start: usize = cap.get(1).unwrap().as_str().parse().unwrap();
  if let Some(end_cap) = cap.get(2) {
     if let Ok(end) = end_cap.as_str().parse() {
         PeptideRange {
            start,
            end,
        }
     } else {
        // matched 1..?:
        PeptideRange {
            start,
            end: start,
        }
     }
  } else {
     PeptideRange {
        start,
        end: start,
     }
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

fn process_record(uniprot_record: UniProtDataRecord) -> UniProtDataEntry {
    let gene_uniquename =
        if uniprot_record.gene_uniquename.ends_with(";") {
            let mut uniquename = uniprot_record.gene_uniquename.to_string();
            uniquename.pop();
            uniquename
        } else {
            uniprot_record.gene_uniquename.to_string()
        };

    let signal_peptide =
        if let Some(field_part) = first_field_part(&uniprot_record.signal_peptide) {
            RANGE_RE.captures_iter(&field_part).next()
                .map(|cap| SignalPeptide {
                    range: get_range(cap),
                })
        } else {
            None
        };

    let transit_peptide =
        if let Some(field_part) = first_field_part(&uniprot_record.transit_peptide) {
            RANGE_RE.captures_iter(&field_part).next()
                .map(|cap| TransitPeptide {
                    range: get_range(cap),
                })
        } else {
            None
        };

    let mut binding_sites_parts_iter =
        SPLIT_RE.split(&uniprot_record.binding_sites);
    binding_sites_parts_iter.next();  // remove blank

    let binding_sites =
        binding_sites_parts_iter.map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                let ligand =
                    if let Some(rest) = cap.get(3) {
                        get_ligand(rest.as_str())
                    } else {
                        flex_str!("unknown ligand")
                    };
                BindingSite {
                    ligand,
                    range: get_range(cap),
                }
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    let mut active_sites_parts_iter =
        SPLIT_RE.split(&uniprot_record.active_sites);
    active_sites_parts_iter.next();  // remove blank

    let active_sites =
        active_sites_parts_iter.map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                ActiveSite {
                    range: get_range(cap),
                }
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    let mut beta_strands_parts_iter =
        SPLIT_RE.split(&uniprot_record.beta_strands);
    beta_strands_parts_iter.next();  // remove blank

    let beta_strands =
        beta_strands_parts_iter.map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                BetaStrand {
                    range: get_range(cap),
                }
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    let mut helices_parts_iter = SPLIT_RE.split(&uniprot_record.helices);
    helices_parts_iter.next();  // remove blank

    let helices =
        helices_parts_iter.map(|field_part| {
            if let Some(cap) = RANGE_RE.captures_iter(field_part).next() {
                Helix {
                    range: get_range(cap),
                }
            } else {
                panic!("failed to parse UniProt data file, no range in {}", field_part);
            }
        })
        .collect();

    UniProtDataEntry {
        gene_uniquename: gene_uniquename.into(),
        signal_peptide,
        transit_peptide,
        binding_sites,
        active_sites,
        beta_strands,
        helices,
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

