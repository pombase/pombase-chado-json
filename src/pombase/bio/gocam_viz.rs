use std::collections::{HashMap, HashSet};

use crate::{data_types::{GoCamId, UniquenameGeneMap}, types::GeneUniquename};

pub fn make_gocam_data(gene_details_maps: &UniquenameGeneMap)
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