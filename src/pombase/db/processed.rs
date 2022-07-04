use std::collections::HashMap;

use flexstr::SharedStr as FlexStr;

use crate::types::{OrganismTaxonId, GeneUniquename, GeneName, AlleleUniquename, AlleleName, ReferenceUniquename};
use crate::web::config::ConfigOrganism;

use super::Raw;

type Organism = ConfigOrganism;
type OrganismTaxonIdMap = HashMap<OrganismTaxonId, Organism>;

fn make_organism_map(raw: &Raw) -> OrganismTaxonIdMap {
    let mut organism_map = HashMap::new();

    for chado_organism in &raw.organisms {
        let mut maybe_taxonid: Option<u32> = None;
        for prop in chado_organism.organismprops.borrow().iter() {
            if prop.prop_type.name == "taxon_id" {
                maybe_taxonid = Some(prop.value.parse().unwrap());
            }
        }
        if let Some(taxonid) = maybe_taxonid {
            let org =
                ConfigOrganism::new(taxonid, chado_organism.genus.clone(),
                                    chado_organism.species.clone(),
                                    vec![], None);
            organism_map.insert(taxonid, org);
        }
    }

    organism_map
}

pub struct Processed {
  organism_map: OrganismTaxonIdMap,

  gene_map: GeneMap,
  allele_map: AlleleMap,

  reference_map: ReferenceMap,

//  term_map: TermMap,

}

/*
type TermId = FlexStr;
type TermName = FlexStr;

pub struct Term {
  id: TermId,
  name: TermName,
}

type TermMap = HashMap<TermId, Term>;

fn make_term_map(raw: &Raw) -> TermMap {
  HashMap::new()
}
 */

pub struct Gene {
  uniquename: GeneUniquename,
  name: Option<GeneName>,
  feature_type: FlexStr,
}

impl Gene {
  pub fn new(uniquename: &GeneUniquename, name: &Option<GeneName>, feature_type: &FlexStr) -> Gene {
    Gene {
      uniquename: uniquename.clone(),
      name: name.clone(),
      feature_type: feature_type.clone(),
    }
  }

  pub fn uniquename(&self) -> GeneUniquename {
    self.uniquename.clone()
  }

  pub fn name(&self) -> Option<GeneName> {
    self.name.clone()
  }

  pub fn feature_type(&self) -> FlexStr {
    self.feature_type.clone()
  }
}

type GeneMap = HashMap<GeneUniquename, Gene>;

pub struct Allele {
  uniquename: AlleleUniquename,
  name: Option<AlleleName>,
}

impl Allele {
  pub fn new(uniquename: &AlleleUniquename, name: &Option<AlleleName>) -> Allele {
    Allele {
      uniquename: uniquename.clone(),
      name: name.clone(),
    }
  }

  pub fn uniquename(&self) -> AlleleUniquename {
    self.uniquename.clone()
  }

  pub fn name(&self) -> Option<AlleleName> {
    self.name.clone()
  }
}

type AlleleMap = HashMap<AlleleUniquename, Allele>;

fn make_feature_maps(raw: &Raw) -> (GeneMap, AlleleMap)
{
  let mut gene_map = HashMap::new();
  let mut allele_map = HashMap::new();

  for rc_feature in &raw.features {
    let feature_type_name = &rc_feature.feat_type.name;

    match feature_type_name.as_str() {
      "gene"|"pseudogene" => {
        gene_map.insert(rc_feature.uniquename.clone(),
                        Gene::new(&rc_feature.uniquename, &rc_feature.name,
                        &rc_feature.feat_type.name));
      },
      "allele" => {
        allele_map.insert(rc_feature.uniquename.clone(),
                          Allele::new(&rc_feature.uniquename, &rc_feature.name));
      },
      _ => (),
    }
  }

  (gene_map, allele_map)
}

pub struct Reference {
  uniquename: ReferenceUniquename,
}

impl Reference {
  pub fn uniquename(&self) -> ReferenceUniquename {
    self.uniquename.clone()
  }
}

type ReferenceMap = HashMap<ReferenceUniquename, Reference>;

fn make_reference_map(raw: &Raw) -> ReferenceMap {
  let mut ref_map = HashMap::new();

  for rc_pub in &raw.publications {
    let reference = Reference {
      uniquename: rc_pub.uniquename.clone(),
    };
    ref_map.insert(rc_pub.uniquename.clone(), reference);
  }

  ref_map
}

impl<'a> Processed {
  pub fn new(raw: Raw) -> Processed {
//    let mut term_map = make_term_map(&raw);

    let organism_map = make_organism_map(&raw);
    let (gene_map, allele_map) = make_feature_maps(&raw);

    let reference_map = make_reference_map(&raw);

    Processed {
        organism_map,

        gene_map,
        allele_map,

        reference_map,

      // term_map,
    }
  }
/*
  pub fn terms(&'a self) -> &'a TermMap {
    &self.term_map
  }
   */

   pub fn organism_by_taxonid(&'a self, taxonid: OrganismTaxonId) -> Option<&'a Organism> {
     self.organism_map.get(&taxonid)
   }

   pub fn gene_by_uniquename(&'a self, uniquename: &GeneUniquename) -> Option<&'a Gene> {
     self.gene_map.get(uniquename)
   }

   pub fn allele_by_uniquename(&'a self, uniquename: &AlleleUniquename) -> Option<&'a Allele> {
     self.allele_map.get(uniquename)
   }

   pub fn reference_by_uniquename(&'a self, uniquename: &ReferenceUniquename) -> Option<&'a Reference> {
    self.reference_map.get(uniquename)
   }
}
