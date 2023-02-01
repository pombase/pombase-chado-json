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
    value: String,
    db_xref: Option<String>,
    pombase_reference: Option<String>,
    pombase_comments: Option<String>,
}

const FEATURE_TYPES: [&str; 8] = ["CDS", "lncRNA", "rRNA", "sncRNA", "snoRNA", "snRNA", "ncRNA", "tRNA"];

pub type GeneHistoryMap = HashMap<GeneUniquename, Vec<GeneHistoryEntry>>;

fn add_entry(map: &mut GeneHistoryMap, record: GeneHistoryFileRecord) {
    let (old_coords, new_coords) =
        if record.added_or_removed == GeneHistoryEntryType::Added {
            (None, Some(record.value))
        } else {
            (Some(record.value), None)
        };

    let mut references = vec![];

    if let Some(reference) = record.pombase_reference {
        references.push(reference.into())
    }

    if let Some(db_xref) = record.db_xref {
        references.push(db_xref.into())
    }

    let entry = GeneHistoryEntry {
        revision: record.revision,
        old_coords,
        new_coords,
        date: record.date.into(),
        entry_type: record.added_or_removed,
        references,
        comment: record.pombase_comments.clone(),
    };

    map.entry(record.systematic_id.into())
        .or_insert_with(Vec::new)
        .push(entry);

}

pub fn parse_gene_history(file_name: &str) -> GeneHistoryMap
{
    let mut res: GeneHistoryMap = HashMap::new();

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

    let mut removed_records = vec![];

    for result in csv_reader.deserialize() {
        let record: GeneHistoryFileRecord =
            result.unwrap_or_else(|e| {
                panic!("failed to read gene history CSV file: {}", e);
            });

        if !FEATURE_TYPES.contains(&record.feature_type.as_ref()) {
            continue;
        }

        if record.added_or_removed == GeneHistoryEntryType::Removed {
            removed_records.push(record);
            continue;
        }

        if record.added_or_removed != GeneHistoryEntryType::Added {
            panic!("expected record type: {:?}", record.added_or_removed);
        }

        add_entry(&mut res, record);
    }

    for removed_record in removed_records {
        let systematic_id: &str = removed_record.systematic_id.as_ref();

        let entries = res.get_mut(systematic_id);

        let existing_entry =
            if let Some(entries) = entries {
                entries.iter_mut().find(|entry_item| {
                    entry_item.entry_type == GeneHistoryEntryType::Added &&
                        entry_item.revision == removed_record.revision
                })
            } else {
                None
            };

        if let Some(existing_entry) = existing_entry {
            existing_entry.entry_type = GeneHistoryEntryType::Changed;
            existing_entry.old_coords = Some(removed_record.value);
        } else {
            add_entry(&mut res, removed_record);
        }
    }

    for entries in res.values_mut() {
        entries.sort_by(|a, b| a.date.cmp(&b.date))
    }

    res
}
