use std::collections::{HashMap, HashSet};

use crate::{data_types::{GoCamDetails, GoCamId, TermIdDetailsMap, UniquenameGeneMap}, types::GeneUniquename, web::config::TermAndName};

pub fn make_gocam_data_by_gene(gene_details_maps: &UniquenameGeneMap)
  -> HashMap<GeneUniquename, HashSet<GoCamId>>
{
  let mut ret = HashMap::new();

  for gene in gene_details_maps.values() {
    if gene.gocam_ids.len() > 0 {
      ret.insert(gene.uniquename.clone(), gene.gocam_ids.clone());
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
    for gocam_id in &gene.gocam_ids {
      ret.entry(gocam_id.clone())
         .or_insert_with(|| GoCamDetails::new(&gocam_id))
         .genes.insert(gene.uniquename.clone());
    }
  }

  for term in term_details_map.values() {
    for gocam_id in &term.gocam_ids {
      ret.entry(gocam_id.clone())
         .or_insert_with(|| GoCamDetails::new(&gocam_id))
         .terms.insert(TermAndName {
                         termid: term.termid.clone(),
                         name: term.name.clone(),
                       });
    }
  }

  ret
}