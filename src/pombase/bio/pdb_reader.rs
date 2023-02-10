use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;
use crate::types::{GeneUniquename, ReferenceUniquename, PdbId};
use crate::data_types::{PDBGeneEntry, PDBRefEntry, PDBGeneChain};

pub type PDBGeneEntryMap = HashMap<GeneUniquename, Vec<PDBGeneEntry>>;
pub type PDBRefEntryMap = HashMap<ReferenceUniquename, HashMap<PdbId, PDBRefEntry>>;

fn make_ref_entry(record: &PDBGeneEntry) -> PDBRefEntry {
    PDBRefEntry {
        pdb_id: record.pdb_id.clone(),
        gene_chains: vec![],
        title: record.title.clone(),
        entry_authors: record.entry_authors.clone(),
        entry_authors_abbrev: record.entry_authors_abbrev.clone(),
        reference: record.reference.clone(),
        experimental_method: record.experimental_method.clone(),
        resolution: record.experimental_method.clone(),
    }
}

fn make_gene_chain(record: &PDBGeneEntry) -> PDBGeneChain {
    PDBGeneChain {
        gene_uniquename: record.gene_uniquename.clone(),
        chain: record.chain.clone(),
        position: record.position.clone(),
    }
}

pub fn read_pdb_data(file_name: &str)
 -> (PDBGeneEntryMap, PDBRefEntryMap)
{
  let mut gene_entry_map = HashMap::new();
  let mut ref_entry_map = HashMap::new();

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
     let record: PDBGeneEntry =
         result.unwrap_or_else(|e| {
             panic!("failed to read PDB data file: {}", e);
         });

     let gene_uniquename = record.gene_uniquename.clone();

     gene_entry_map.entry(gene_uniquename)
        .or_insert_with(Vec::new)
        .push(record.clone());

     let ref_uniquename = record.reference.clone();

     ref_entry_map.entry(ref_uniquename)
        .or_insert_with(HashMap::new)
        .entry(record.pdb_id.clone())
        .or_insert_with(|| make_ref_entry(&record))
        .gene_chains.push(make_gene_chain(&record));
  }

  (gene_entry_map, ref_entry_map)
}