use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::process;
use crate::types::{GeneUniquename, ReferenceUniquename, PdbId};
use crate::data_types::{PDBEntry, PDBGeneChain};


#[derive(Serialize, Deserialize, Clone, Debug)]
struct PDBSourceEntry {
    pub gene_uniquename: GeneUniquename,
    pub pdb_id: PdbId,
    pub title: String,
    pub entry_authors: String,
    pub entry_authors_abbrev: String,
    pub reference: ReferenceUniquename,
    pub experimental_method: String,
    pub resolution: String,
    pub chain: String,
    pub position: String,
}


pub type PDBGeneEntryMap = HashMap<GeneUniquename, Vec<PDBEntry>>;
pub type PDBRefEntryMap = HashMap<ReferenceUniquename, HashMap<PdbId, PDBEntry>>;

fn make_entry(record: &PDBSourceEntry) -> PDBEntry {
    let gene_chain = make_gene_chain(record);

    PDBEntry {
        pdb_id: record.pdb_id.clone(),
        gene_chains: vec![gene_chain],
        title: record.title.clone(),
        entry_authors: record.entry_authors.clone(),
        entry_authors_abbrev: record.entry_authors_abbrev.clone(),
        reference: record.reference.clone(),
        experimental_method: record.experimental_method.clone(),
        resolution: record.resolution.clone(),
    }
}

fn make_gene_chain(record: &PDBSourceEntry) -> PDBGeneChain {
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
     let record: PDBSourceEntry =
         result.unwrap_or_else(|e| {
             panic!("failed to read PDB data file: {}", e);
         });

     let gene_uniquename = record.gene_uniquename.clone();

     gene_entry_map.entry(gene_uniquename)
        .or_insert_with(Vec::new)
        .push(make_entry(&record));

     let ref_uniquename = record.reference.clone();

     ref_entry_map.entry(ref_uniquename)
        .or_insert_with(HashMap::new)
        .entry(record.pdb_id.clone())
        .or_insert_with(|| make_entry(&record))
        .gene_chains.push(make_gene_chain(&record));
  }

  (gene_entry_map, ref_entry_map)
}