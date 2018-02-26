use std::process;
use std::fs::File;
use std::io::{Read, BufReader};

use serde_json;
use std::collections::HashMap;
use std::collections::HashSet;

use web::data::{APIMaps, IdGeneSubsetMap, APIGeneSummary, APIAlleleDetails,
                GeneDetails, TermDetails, GenotypeDetails, ReferenceDetails,
                InteractionType, OntAnnotationMap, IdOntAnnotationDetailMap,
                TermShort, TermShortOptionMap,
                ReferenceShort, ReferenceShortOptionMap,
                GeneShort, GeneShortOptionMap};
use web::config::Config;
use api::query::{SingleOrMultiAllele, QueryExpressionFilter};

use flate2::read::GzDecoder;

use types::GeneUniquename;

pub struct ServerData {
    config_file_name: String,
    maps: APIMaps,
    search_maps_file_name: String,
    gene_subsets: IdGeneSubsetMap,
    gene_subsets_file_name: String,
}


fn load(config: &Config, search_maps_file_name: &str, gene_subsets_file_name: &str)
        -> (APIMaps, IdGeneSubsetMap)
{
    let file = match File::open(search_maps_file_name) {
        Ok(file) => file,
        Err(err) => {
            print!("Failed to read {}: {}\n", search_maps_file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    let mut decoder = GzDecoder::new(reader).unwrap();
    let mut decoded_json = String::new();
    decoder.read_to_string(&mut decoded_json).unwrap();

    let query_api_maps: APIMaps =
        match serde_json::from_str(&decoded_json) {
            Ok(results) => results,
            Err(err) => {
                print!("failed to parse {}: {}", search_maps_file_name, err);
                process::exit(1);
            },
        };

    let file = match File::open(gene_subsets_file_name) {
        Ok(file) => file,
        Err(err) => {
            print!("Failed to read {}: {}\n", gene_subsets_file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    let mut gene_subsets: IdGeneSubsetMap =
        match serde_json::from_reader(reader) {
            Ok(results) => results,
            Err(err) => {
                print!("failed to parse {}: {}", gene_subsets_file_name, err);
                process::exit(1);
            },
        };

    let mut new_entries: IdGeneSubsetMap = HashMap::new();

    let prefixes_to_remove: Vec<String> =
        config.server.subsets.prefixes_to_remove
        .iter().map(|prefix| prefix.clone() + ":").collect();

    // strip prefixes and add to map
    for (subset_name, subset_details) in gene_subsets.iter() {
        for prefix in &prefixes_to_remove {
            if subset_name.starts_with(prefix) {
                let new_subset_name = &subset_name[prefix.len()..];
                new_entries.insert(String::from(new_subset_name), subset_details.clone());
            }
        }
    }

    gene_subsets.extend(new_entries);

    (query_api_maps, gene_subsets)
}

impl ServerData {
    pub fn new(config_file_name: &str, search_maps_file_name: &str,
               gene_subsets_file_name: &str)
               -> ServerData
    {
        let config = Config::read(config_file_name);

        let (maps, gene_subsets) =
            load(&config, search_maps_file_name, gene_subsets_file_name);
        ServerData {
            config_file_name: config_file_name.into(),
            search_maps_file_name: search_maps_file_name.into(),
            maps: maps,
            gene_subsets_file_name: gene_subsets_file_name.into(),
            gene_subsets: gene_subsets,
        }
    }

    pub fn gene_uniquename_of_id(&self, id: &str) -> Option<GeneUniquename> {
        if self.maps.gene_summaries.contains_key(id) {
            Some(id.to_owned())
        } else {
            if let Some(gene_uniquename) =
                self.maps.gene_name_gene_map.get(id) {
                    Some(gene_uniquename.clone())
                } else {
                    let lower_case_id = id.to_lowercase();
                    for gene_uniquename in self.maps.gene_summaries.keys() {
                        if gene_uniquename.to_lowercase() == lower_case_id {
                            return Some(gene_uniquename.clone());
                        }
                    }
                    None
                }
        }
    }

    pub fn get_gene_summary(&self, gene_uniquename: &str) -> Option<&APIGeneSummary> {
        self.maps.gene_summaries.get(gene_uniquename)
    }

    pub fn genes_of_termid(&self, term_id: &str) -> Vec<GeneUniquename> {
        match self.maps.termid_genes.get(term_id) {
            Some(gene_uniquenames) => {
                gene_uniquenames.iter().cloned().collect::<Vec<_>>()
            },
            None => vec![],
        }
    }

    pub fn genes_of_subset(&self, search_name: &str) -> Vec<GeneUniquename> {
        if search_name.starts_with('!') || search_name.ends_with('*') {
            let mut trimmed_search_name = search_name.to_owned();
            let invert_search = search_name.starts_with('!');
            let wildcard = search_name.ends_with('*');
            if invert_search {
                trimmed_search_name.remove(0);
            }
            if wildcard {
                trimmed_search_name.pop();
            }
            let mut genes = HashSet::new();
            for (subset_name, subset_details) in &self.gene_subsets {
                let name_matches =
                    wildcard && subset_name.starts_with(&trimmed_search_name) ||
                    trimmed_search_name.eq(subset_name);

                if !invert_search && name_matches || invert_search && !name_matches {
                    genes.extend(subset_details.elements.iter().cloned());
                }
            }
            genes.iter().cloned().collect::<Vec<_>>()
        } else {
            match self.gene_subsets.get(search_name) {
                Some(subset_details) => {
                    subset_details.elements.iter().cloned().collect::<Vec<_>>()
                },
                None => vec![],
            }
        }
    }

    pub fn genes_of_genotypes(&self, term_id: &str,
                              single_or_multi_allele: &SingleOrMultiAllele,
                              expression_filter: &Option<QueryExpressionFilter>)
                              -> Vec<GeneUniquename>
    {
        if let Some(annotations) = self.maps.termid_genotype_annotation.get(term_id) {
            let mut genes: HashSet<GeneUniquename> = HashSet::new();
            for annotation in annotations {

                let mut add_single = false;
                let mut add_multi = false;

                match *single_or_multi_allele {
                    SingleOrMultiAllele::Single => add_single = true,
                    SingleOrMultiAllele::Multi => add_multi = true,
                    SingleOrMultiAllele::Both => {
                        add_single = true;
                        add_multi = true;
                    },
                }

                let expression_matches =
                    |allele_details: &APIAlleleDetails| {
                        let expression = &allele_details.expression;
                        let allele_type = &allele_details.allele_type;
                        if let Some(ref expression_filter) = *expression_filter {
                            if *expression_filter == QueryExpressionFilter::Any {
                                return true;
                            }
                            if *expression_filter == QueryExpressionFilter::Null {
                                if allele_type == "deletion" {
                                    return true;
                                }
                                if let Some(ref expression) = *expression {
                                    expression == "Null"
                                } else {
                                    false
                                }
                            } else {
                                if *expression_filter == QueryExpressionFilter::WtOverexpressed {
                                    if let Some(ref expression) = *expression {
                                        expression == "Overexpression" &&
                                            allele_type == "wild_type"
                                    } else {
                                        false
                                    }
                                } else {
                                    false
                                }
                            }
                        } else {
                            true
                        }
                    };

                if annotation.is_multi && add_multi ||
                    !annotation.is_multi && add_single &&
                    expression_matches(&annotation.alleles[0]) {
                        for allele_details in &annotation.alleles {
                            genes.insert(allele_details.gene.clone());
                        }
                }
            }
            genes.into_iter().collect::<Vec<_>>()
        } else {
            vec![]
        }
    }

    pub fn filter_genes(&self, p: &Fn(&APIGeneSummary) -> bool)
                        -> Vec<GeneUniquename>
    {
        self.maps.gene_summaries.values()
            .filter(|ref summ| p(summ))
            .map(|ref summ| summ.uniquename.clone())
            .collect()
    }

    fn detail_map_of_cv_annotations(&self, ont_annotation_map: &OntAnnotationMap)
                                    -> IdOntAnnotationDetailMap
    {
        let mut details_map = HashMap::new();

        for (_, term_annotations) in ont_annotation_map {
            for term_annotation in term_annotations {
                for annotation_detail_id in &term_annotation.annotations {
                    let details =
                        self.maps.annotation_details.get(annotation_detail_id).unwrap();
                    details_map.insert(*annotation_detail_id, details.clone());
                }
            }
        }

        details_map
    }

    fn fill_term_map(&self, term_map: &TermShortOptionMap) -> TermShortOptionMap {
        let mut ret = term_map.clone();
        for (termid, _) in term_map {
            let term_details = self.maps.terms.get(termid).unwrap();
            let term_short = TermShort::from_term_details(term_details);
            ret.insert(termid.clone(), Some(term_short));
        }
        ret
    }

    fn fill_gene_map(&self, gene_map: &GeneShortOptionMap) -> GeneShortOptionMap {
        let mut ret = gene_map.clone();
        for (gene_uniquename, _) in gene_map {
            let gene_details = self.maps.genes.get(gene_uniquename).unwrap();
            let gene_short = GeneShort::from_gene_details(gene_details);
            ret.insert(gene_uniquename.clone(), Some(gene_short));
        }
        ret
    }

    fn fill_reference_map(&self, reference_map: &ReferenceShortOptionMap) -> ReferenceShortOptionMap {
        let mut ret = reference_map.clone();
        for (reference_uniquename, _) in reference_map {
            let reference_details = self.maps.references.get(reference_uniquename).unwrap();
            let reference_short = ReferenceShort::from_reference_details(reference_details);
            ret.insert(reference_uniquename.clone(), Some(reference_short));
        }
        ret
    }

    pub fn get_gene_details(&self, gene_uniquename: &str) -> Option<GeneDetails> {
        if let Some(gene_ref) = self.maps.genes.get(gene_uniquename) {
            let mut gene = gene_ref.clone();
            let details_map = self.detail_map_of_cv_annotations(&gene.cv_annotations);
            gene.annotation_details = details_map;
            gene.terms_by_termid = self.fill_term_map(&gene.terms_by_termid);
            gene.genes_by_uniquename = self.fill_gene_map(&gene.genes_by_uniquename);
            gene.references_by_uniquename =
                self.fill_reference_map(&gene.references_by_uniquename);
            Some(gene)
        } else {
            None
        }
    }

    pub fn get_genotype_details(&self, genotype_uniquename: &str) -> Option<GenotypeDetails> {
        if let Some(genotype_ref) = self.maps.genotypes.get(genotype_uniquename) {
            let mut genotype = genotype_ref.clone();
            let details_map = self.detail_map_of_cv_annotations(&genotype.cv_annotations);
            genotype.annotation_details = details_map;
            genotype.terms_by_termid = self.fill_term_map(&genotype.terms_by_termid);
            genotype.genes_by_uniquename = self.fill_gene_map(&genotype.genes_by_uniquename);
            genotype.references_by_uniquename =
                self.fill_reference_map(&genotype.references_by_uniquename);
            Some(genotype)
        } else {
            None
        }
    }

    pub fn get_term_details(&self, termid: &str) -> Option<TermDetails> {
        if let Some(term_ref) = self.maps.terms.get(termid) {
            let mut term = term_ref.clone();
            let details_map = self.detail_map_of_cv_annotations(&term.cv_annotations);
            term.annotation_details = details_map;
            term.terms_by_termid = self.fill_term_map(&term.terms_by_termid);
            term.genes_by_uniquename = self.fill_gene_map(&term.genes_by_uniquename);
            term.references_by_uniquename =
                self.fill_reference_map(&term.references_by_uniquename);
            Some(term)
        } else {
            None
        }
    }

    pub fn get_reference_details(&self, reference_uniquename: &str) -> Option<ReferenceDetails> {
        if let Some(reference_ref) = self.maps.references.get(reference_uniquename) {
            let mut reference = reference_ref.clone();
            let details_map = self.detail_map_of_cv_annotations(&reference.cv_annotations);
            reference.annotation_details = details_map;
            reference.terms_by_termid = self.fill_term_map(&reference.terms_by_termid);
            reference.genes_by_uniquename = self.fill_gene_map(&reference.genes_by_uniquename);
            Some(reference)
        } else {
            None
        }
    }

    pub fn interactors_of_genes(&self, gene_uniquename: &GeneUniquename,
                                interaction_type: InteractionType)
                                -> Vec<GeneUniquename> {
        if let Some(interactors) = self.maps.interactors_of_genes.get(gene_uniquename) {
            interactors.iter()
                .filter(|interactor| {
                    interactor.interaction_type == interaction_type
                })
                .map(|interactor| interactor.interactor_uniquename.clone())
                .collect::<Vec<_>>()
        } else {
            vec![]
        }
    }

    pub fn reload(&mut self) {
        let config = Config::read(&self.config_file_name);
        let (maps, gene_subsets) =
            load(&config, &self.search_maps_file_name,
                 &self.gene_subsets_file_name);
        self.maps = maps;
        self.gene_subsets = gene_subsets;
    }
}
