use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;

use regex::Regex;

use crate::data_types::{PeptideRange, SignalPeptide, TransitPeptide};
use crate::types::GeneUniquename;

#[derive(Debug, Serialize, Deserialize)]
pub struct UniProtDataEntry {
    pub gene_uniquename: GeneUniquename,
    pub signal_peptide: Option<SignalPeptide>,
    pub transit_peptide: Option<TransitPeptide>,
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

/*
 Binding site
 Active site
 Catalytic activity
 Gene Names (synonym)
 Post-translational modification
 Modified residue
 Cofactor
 Kinetics
 */
}

lazy_static! {
    static ref SIGNAL_PEPTIDE_RE: Regex = Regex::new(r"^SIGNAL (\d+)\.\.(\d+)").unwrap();
    static ref TRANSIT_PEPTIDE_RE: Regex = Regex::new(r"^TRANSIT (\d+)\.\.(\d+)").unwrap();
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

    if let Some(cap) = SIGNAL_PEPTIDE_RE.captures_iter(&uniprot_record.signal_peptide).next() {
        let start: usize = cap.get(1).unwrap().as_str().parse().unwrap();
        let end: usize = cap.get(2).unwrap().as_str().parse().unwrap();

        signal_peptide = Some(SignalPeptide {
            range: PeptideRange {
                start,
                end,
            }
        });
    };

    let mut transit_peptide = None;

    if let Some(cap) = TRANSIT_PEPTIDE_RE.captures_iter(&uniprot_record.transit_peptide).next() {
        let start: usize = cap.get(1).unwrap().as_str().parse().unwrap();
        let end: usize = cap.get(2).unwrap().as_str().parse().unwrap();

        transit_peptide = Some(TransitPeptide {
            range: PeptideRange {
                start,
                end,
            }
        });
    };

    UniProtDataEntry {
        gene_uniquename: gene_uniquename.into(),
        signal_peptide,
        transit_peptide,
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

