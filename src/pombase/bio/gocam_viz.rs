use std::collections::{HashMap, HashSet};

use crate::{data_types::{GoCamDetails, GoCamId, TermIdDetailsMap, UniquenameGeneMap}, types::GeneUniquename, web::config::TermAndName};

pub fn make_gocam_data_by_gene(gene_details_maps: &UniquenameGeneMap)
  -> HashMap<GeneUniquename, HashSet<GoCamId>>
{
  let mut ret = HashMap::new();

  for gene in gene_details_maps.values() {
    if gene.gocams.len() > 0 {
      let gocam_ids = gene.gocams.iter().map(|gocam| gocam.gocam_id.clone()).collect();
      ret.insert(gene.uniquename.clone(), gocam_ids);
    }
  }

  ret
}

pub fn make_gocam_data_by_id(gene_details_map: &UniquenameGeneMap,
                             term_details_map: &TermIdDetailsMap)
  -> HashMap<GoCamId, GoCamDetails>
{
  let mut ret = HashMap::new();

  for gene in gene_details_map.values() {
    for gocam in &gene.gocams {
      ret.entry(gocam.gocam_id.clone())
         .or_insert_with(|| GoCamDetails::new(&gocam.gocam_id, gocam.title.clone()))
         .genes.insert(gene.uniquename.clone());
    }
  }

  for term in term_details_map.values() {
    for gocam in &term.gocams {
      ret.entry(gocam.gocam_id.clone())
         .or_insert_with(|| GoCamDetails::new(&gocam.gocam_id, gocam.title.clone()))
         .terms.insert(TermAndName {
                         termid: term.termid.clone(),
                         name: term.name.clone(),
                       });
    }
  }

  ret
}