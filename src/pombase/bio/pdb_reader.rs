use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;
use crate::data_types::{PdbId, PDBEntry};

pub type PDBEntryMap = HashMap<PdbId, PDBEntry>;

pub fn read_pdb_data(file_name: &str) -> PDBEntryMap {
  let mut map = HashMap::new();

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
     let record: PDBEntry =
         result.unwrap_or_else(|e| {
             panic!("failed to read PDB data file: {}", e);
         });

     map.insert(record.pdb_id.clone(), record);
  }

  map
}