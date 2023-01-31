use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;

use crate::{data_types::{GeneHistoryEntry, GeneHistoryEntryType}, types::GeneUniquename};

#[derive(Debug, Deserialize)]
struct GeneHistoryFileRecord {
    revision: String,
    date: String,
    systematic_id: String,
    feature_type: String,
    added_or_removed: GeneHistoryEntryType,
    db_xref: Option<String>,
    pombase_reference: Option<String>,
    pombase_comments: Option<String>,
}


pub fn parse_gene_history(file_name: &str) -> HashMap<GeneUniquename, Vec<GeneHistoryEntry>>
{
    let mut res: HashMap<GeneUniquename, Vec<GeneHistoryEntry>> = HashMap::new();

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

    for result in csv_reader.deserialize() {
        let record: GeneHistoryFileRecord =
            result.unwrap_or_else(|e| {
                panic!("failed to read gene history CSV file: {}", e);
            });

        if record.feature_type != "CDS" {
            continue;
        }

        let systematic_id: &str = record.systematic_id.as_ref();

        let entries = res.get_mut(systematic_id);

        let existing_entry =
            if let Some(entries) = entries {
                if entries.iter().filter(|entry_item| { entry_item.revision == record.revision }).count() > 1 {
                    panic!("internal error, more than one item at revision {} for {}",
                           record.revision, record.systematic_id);
                }

                entries.iter_mut().find(|entry_item| {
                    entry_item.revision == record.revision
                })
            } else {
                None
            };

        if let Some(existing_entry) = existing_entry {
            if record.added_or_removed == GeneHistoryEntryType::Added {
                if existing_entry.entry_type == GeneHistoryEntryType::Removed {
                    existing_entry.entry_type = GeneHistoryEntryType::Changed
                } else {
                    eprintln!("gene {} added twice in revision {}",
                              record.systematic_id, record.revision);
                }
            } else if record.added_or_removed == GeneHistoryEntryType::Removed {
                if existing_entry.entry_type == GeneHistoryEntryType::Added {
                    existing_entry.entry_type = GeneHistoryEntryType::Changed;
                } else {
                    eprintln!("gene {} removed twice in revision {}",
                              record.systematic_id, record.revision);
                }
            }
        } else {

            let mut references = vec![];

            if let Some(reference) = record.pombase_reference {
                references.push(reference.into())
            }

            if let Some(db_xref) = record.db_xref {
                references.push(db_xref.into())
            }

            let entry = GeneHistoryEntry {
                revision: record.revision,
                date: record.date.into(),
                entry_type: record.added_or_removed,
                references,
                comment: record.pombase_comments.clone(),
            };

            res.entry(record.systematic_id.into())
                .or_insert_with(Vec::new)
                .push(entry);
        };
    }

    res
}