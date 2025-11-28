use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;

use crate::types::GeneUniquename;

use crate::data_types::GeneHistoryEntry;

#[derive(Debug, Deserialize)]
struct GeneHistoryFileRecord {
     date: String,
    systematic_id: String,
    feature_type: String,
    value: String,
    db_xref_added: Option<String>,
    pombase_reference_added: Option<String>,
    pombase_comments_added: Option<String>,
    genome_snapshot: Option<String>,
}

const FEATURE_TYPES: [&str; 8] = ["CDS", "lncRNA", "rRNA", "sncRNA", "snoRNA", "snRNA", "ncRNA", "tRNA"];

pub type GeneHistoryMap = HashMap<GeneUniquename, Vec<GeneHistoryEntry>>;

fn add_entry(map: &mut GeneHistoryMap, record: GeneHistoryFileRecord) {

    let mut references = vec![];

    if let Some(ref reference) = record.pombase_reference_added {
        references.push(reference.into())
    }

    if let Some(db_xref) = record.db_xref_added {
        if let Some(reference) = record.pombase_reference_added {
            if !reference.contains(&db_xref) {
                references.push(db_xref.into())
            }
        } else {
            references.push(db_xref.into())
        }
    }

    let entry = GeneHistoryEntry {
        previous_coords: record.value,
        date: record.date.into(),
        references,
        comments: record.pombase_comments_added,
        genome_snapshot_link: record.genome_snapshot,
    };

    map.entry(record.systematic_id.into())
        .or_default()
        .push(entry);

}

pub fn parse_gene_history(file_name: &str) -> GeneHistoryMap
{
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            println!("Failed to read {}: {}", file_name, err);
            process::exit(1);
        }
    };

    let reader = BufReader::new(file);

    let mut csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(reader);

    let mut res: GeneHistoryMap = HashMap::new();

    for result in csv_reader.deserialize() {
        let record: GeneHistoryFileRecord =
            result.unwrap_or_else(|e| {
                panic!("failed to read gene history CSV file: {}", e);
            });

        if !FEATURE_TYPES.contains(&record.feature_type.as_ref()) {
            continue;
        }

        add_entry(&mut res, record);
    }


    for entries in res.values_mut() {
        entries.sort_by(|a, b| a.date.cmp(&b.date))
    }

    res
}
