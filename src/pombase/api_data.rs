use std::process;
use std::fs::File;
use std::io::{Read, BufReader};

use serde_json;
use std::collections::{HashMap, HashSet};

use crate::data_types::{APIMaps, IdGeneSubsetMap, APIGeneSummary, APIAlleleDetails,
                GeneDetails, TermDetails, GenotypeDetails, ReferenceDetails,
                InteractionType, OntAnnotationMap, IdOntAnnotationDetailMap,
                OntAnnotationDetail,
                TermShort, TermShortOptionMap, ChromosomeDetails,
                ReferenceShort, ReferenceShortOptionMap,
                FeatureShort,
                GeneShort, GeneShortOptionMap, GeneQueryData,
                ExtPart, ExtRange, GeneAndGeneProduct, WithFromValue,
                Ploidiness};
use crate::web::config::{Config, TermAndName};
use crate::api::query::{SingleOrMultiLocus, QueryExpressionFilter};
use crate::web::cv_summary::make_cv_summaries;

use flate2::read::GzDecoder;

use crate::types::{TermId, GeneUniquename};

use pombase_rc_string::RcString;

pub struct APIData {
    config: Config,
    maps: APIMaps,
    secondary_identifiers_map: HashMap<TermId, TermId>,
}


fn make_secondary_identifiers_map(terms: &HashMap<TermId, TermDetails>)
                                  -> HashMap<TermId, TermId>
{
    let mut ret_map = HashMap::new();

    for term_details in terms.values() {
        for secondary_identifier in &term_details.secondary_identifiers {
            ret_map.insert(secondary_identifier.clone(),
                           term_details.termid.clone());
        }
    }

    ret_map
}


pub fn api_maps_from_file(search_maps_file_name: &str) -> APIMaps
{
    let file = match File::open(search_maps_file_name) {
        Ok(file) => file,
        Err(err) => {
            eprint!("Failed to read {}: {}\n", search_maps_file_name, err);
            process::exit(1);
        }
    };

    let reader = BufReader::new(file);
    let mut decoder = GzDecoder::new(reader);

    //  this uses less peak memory but is 4X slower
    //  See: https://github.com/serde-rs/json/issues/160
    //        let serde_result = serde_json::de::from_reader(&mut decoder);

    let mut decoded_json = String::new();
    decoder.read_to_string(&mut decoded_json).unwrap();
    let serde_result = serde_json::from_str(&decoded_json);

    match serde_result {
        Ok(results) => results,
        Err(err) => {
            eprint!("failed to parse {}: {}", search_maps_file_name, err);
            process::exit(1);
        },
    }
}

impl APIData {
    pub fn new(config: &Config, maps: &APIMaps)
               ->APIData
    {
        let mut maps = maps.clone();
        let mut new_entries: IdGeneSubsetMap = HashMap::new();

        let prefixes_to_remove: Vec<String> =
            config.server.subsets.prefixes_to_remove
            .iter().map(|prefix| prefix.clone() + ":").collect();

        // strip prefixes and add to map
        for (subset_name, subset_details) in &maps.gene_subsets {
            for prefix in &prefixes_to_remove {
                if subset_name.starts_with(prefix) {
                    let new_subset_name = &subset_name[prefix.len()..];
                    new_entries.insert(RcString::from(new_subset_name), subset_details.clone());
                }
            }
        }

        maps.gene_subsets.extend(new_entries);

        let secondary_identifiers_map =
            make_secondary_identifiers_map(&maps.terms);

        APIData {
            config: config.clone(),
            secondary_identifiers_map,
            maps
        }
    }

    pub fn get_config(&self) -> &Config {
        &self.config
    }

    pub fn get_maps(&self) -> &APIMaps {
        &self.maps
    }

    pub fn gene_uniquename_of_id(&self, id: &RcString) -> Option<GeneUniquename> {
        if self.maps.gene_summaries.contains_key(id) {
            Some(id.clone())
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

    pub fn get_gene_details(&self, gene_uniquename: &str) -> Option<&GeneDetails> {
        self.maps.genes.get(gene_uniquename)
    }

    pub fn get_gene_query_data(&self, gene_uniquename: &str) -> Option<&GeneQueryData> {
        self.maps.gene_query_data_map.get(gene_uniquename)
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
            for (subset_name, subset_details) in &self.maps.gene_subsets {
                let name_matches =
                    wildcard && subset_name.starts_with(&trimmed_search_name) ||
                    trimmed_search_name.eq(subset_name);

                if !invert_search && name_matches || invert_search && !name_matches {
                    genes.extend(subset_details.elements.iter().cloned());
                }
            }
            genes.iter().cloned().collect::<Vec<_>>()
        } else {
            match self.maps.gene_subsets.get(search_name) {
                Some(subset_details) => {
                    subset_details.elements.iter().cloned().collect::<Vec<_>>()
                },
                None => vec![],
            }
        }
    }

    pub fn genes_of_genotypes(&self, term_id: &str,
                              single_or_multi_locus: &SingleOrMultiLocus,
                              query_ploidiness: &Ploidiness,
                              expression_filter: &Option<QueryExpressionFilter>,
                              conditions_filter: &HashSet<TermAndName>,
                              excluded_conditions_filter: &HashSet<TermAndName>)
                              -> Vec<GeneUniquename>
    {
        if let Some(annotations) = self.maps.termid_genotype_annotation.get(term_id) {
            let mut genes: HashSet<GeneUniquename> = HashSet::new();
            for annotation in annotations {

                let mut add_single = false;
                let mut add_multi = false;

                match *single_or_multi_locus {
                    SingleOrMultiLocus::Single => add_single = true,
                    SingleOrMultiLocus::Multi => add_multi = true,
                    SingleOrMultiLocus::Both => {
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

                let condition_matches =
                    |conditions: &HashSet<TermAndName>, test_conditions_filter: &HashSet<TermAndName>| {
                        if test_conditions_filter.len() == 0 {
                            // don't filter in this case
                            return true;
                        }

                        for filter_condition in test_conditions_filter {
                            if conditions.iter().find(
                                |cond| {
                                    filter_condition.termid == cond.termid
                                })
                                .is_some() {
                                    return true
                                }
                        }

                        false
                    };

                if !condition_matches(&annotation.conditions, conditions_filter) ||
                    excluded_conditions_filter.len() > 0 &&
                    condition_matches(&annotation.conditions, excluded_conditions_filter) {
                        continue;
                    }


                let mut add_genotype_genes = false;

                if annotation.is_multi && add_multi {
                    if *query_ploidiness == annotation.ploidiness ||
                        *query_ploidiness == Ploidiness::Any
                    {
                        add_genotype_genes = true;
                    }
                }

                if !annotation.is_multi && add_single {
                    if (*query_ploidiness != Ploidiness::Haploid &&
                        annotation.ploidiness == Ploidiness::Diploid) ||
                        (annotation.ploidiness == Ploidiness::Haploid &&
                         (*query_ploidiness == Ploidiness::Any ||
                          (*query_ploidiness == Ploidiness::Haploid &&
                           expression_matches(&annotation.alleles[0]))))
                        {
                            add_genotype_genes = true;
                        }
                }

                if add_genotype_genes {
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

    pub fn filter_genes(&self, p: & dyn Fn(&APIGeneSummary) -> bool)
                        -> Vec<GeneUniquename>
    {
        self.maps.gene_summaries.values()
            .filter(|summ| p(summ))
            .map(|summ| summ.uniquename.clone())
            .collect()
    }


    fn strip_db_prefix(&self, uniquename: &RcString) -> RcString {
        if uniquename.starts_with("PomBase:SP") {
            RcString::from(&uniquename[8..])
        } else {
            uniquename.clone()
        }
    }

    // return a fake extension for "with" properties on protein binding annotations
    fn get_with_extension(&self, with_value: WithFromValue) -> ExtPart {
        let ext_range =
            match with_value {
                WithFromValue::Gene(gene_short) => {
                    let uniquename = &gene_short.uniquename;
                    ExtRange::Gene(self.strip_db_prefix(&uniquename))
                },
                _ => panic!("unexpected WithFromValue variant: {:#?}", with_value),
            };

        // a with property on a protein binding (GO:0005515) is
        // displayed as a binds extension
        // https://github.com/pombase/website/issues/108
        ExtPart {
            rel_type_id: None,
            rel_type_name: "binds".into(),
            rel_type_display_name: "binds".into(),
            ext_range,
        }
    }

    // add the with value as a fake extension if the cvterm is_a protein binding and
    // remove the with value
    fn maybe_move_with(&self, term_details: &TermDetails,
                       annotation: &mut OntAnnotationDetail) {
        if let Some(ref evidence_code) = annotation.evidence {
            if evidence_code == "IPI" &&
                (term_details.termid == "GO:0005515" ||
                 term_details.interesting_parent_ids.contains("GO:0005515") ||
                 term_details.termid == "GO:0003723" ||
                 term_details.interesting_parent_ids.contains("GO:0003723"))
            {
                let mut first_with = None;

                for with_value in annotation.withs.drain() {
                    if first_with.is_none() {
                        first_with = Some(with_value);
                    } else {
                        annotation.extension.insert(0, self.get_with_extension(with_value));
                    }
                }

                // if there is an ExtRange that is a PRO ID (ExtRange::GeneProduct),
                // and there is a with value, combine them into an
                // ExtRange::GeneAndGeneProduct
                if let Some(with_value) = first_with.clone() {
                    annotation.extension.iter_mut().for_each(|mut ext_part| {
                        if ext_part.rel_type_name == "has_direct_input" {
                            match &mut ext_part.ext_range {
                                ExtRange::GeneProduct(range_termid) => {
                                    match &with_value {
                                        WithFromValue::Gene(gene_short) => {
                                            let gene_uniquename =
                                                self.strip_db_prefix(&gene_short.uniquename);
                                            let val = GeneAndGeneProduct {
                                                product: range_termid.clone(),
                                                gene_uniquename,
                                            };
                                            ext_part.ext_range = ExtRange::GeneAndGeneProduct(val);
                                            first_with = None;
                                        },
                                        _ => (),
                                    }
                                },
                                _ => ()
                            }
                        }
                    });
                }

                if let Some(with_value) = first_with {
                    annotation.extension.insert(0, self.get_with_extension(with_value));
                }
            }
        }
    }


    fn detail_map_of_cv_annotations(&self, ont_annotation_map: &OntAnnotationMap)
                                    -> IdOntAnnotationDetailMap
    {
        let mut details_map = HashMap::new();

        for term_annotations in ont_annotation_map.values() {
            for term_annotation in term_annotations {
                let termid = &term_annotation.term;
                if let Some(ref term_details) = self.maps.terms.get(termid) {
                    for annotation_detail_id in &term_annotation.annotations {
                        let mut annotation =
                            self.maps.annotation_details[annotation_detail_id].clone();
                        self.maybe_move_with(term_details, &mut annotation);
                        details_map.insert(*annotation_detail_id, annotation);
                    }
                } else {
                    panic!("failed to find TermDetails for {}", termid);
                }
            }
        }

        details_map
    }

    fn fill_term_map(&self, term_map: &TermShortOptionMap) -> TermShortOptionMap {
        let mut ret = term_map.clone();
        for termid in term_map.keys() {
            if let Some(term_details) = self.maps.terms.get(termid) {
                let term_short = TermShort::from_term_details(term_details);
                ret.insert(termid.clone(), Some(term_short));
            } else {
                eprint!("WARNING missing short info for: {}\n", termid);
                ret.insert(termid.clone(), None);
            }
        }
        ret
    }

    fn fill_gene_map(&self, gene_map: &GeneShortOptionMap) -> GeneShortOptionMap {
        let mut ret = gene_map.clone();
        for gene_uniquename in gene_map.keys() {
            let gene_details = &self.maps.genes[gene_uniquename];
            let gene_short = GeneShort::from_gene_details(gene_details);
            ret.insert(gene_uniquename.clone(), Some(gene_short));
        }
        ret
    }

    fn fill_reference_map(&self, reference_map: &ReferenceShortOptionMap) -> ReferenceShortOptionMap {
        let mut ret = reference_map.clone();
        for reference_uniquename in reference_map.keys() {
            let reference_details = &self.maps.references[reference_uniquename];
            let reference_short = ReferenceShort::from_reference_details(reference_details);
            ret.insert(reference_uniquename.clone(), Some(reference_short));
        }
        ret
    }

    // return a GeneDetails object with the term, genes and references maps filled in
    pub fn get_full_gene_details(&self, gene_uniquename: &str) -> Option<GeneDetails> {
        if let Some(gene_ref) = self.maps.genes.get(gene_uniquename) {
            let mut gene = gene_ref.clone();
            let details_map = self.detail_map_of_cv_annotations(&gene.cv_annotations);
            gene.terms_by_termid = self.fill_term_map(&gene.terms_by_termid);
            gene.genes_by_uniquename = self.fill_gene_map(&gene.genes_by_uniquename);
            gene.references_by_uniquename =
                self.fill_reference_map(&gene.references_by_uniquename);
            make_cv_summaries(&mut gene, &self.config, &self.maps.children_by_termid,
                              false, true, &self.maps.genes, &self.maps.genotypes,
                              &details_map);

            gene.annotation_details = details_map;
            Some(gene)
        } else {
            None
        }
    }

    pub fn get_genotype_details(&self, genotype_uniquename: &str) -> Option<GenotypeDetails> {
        if let Some(genotype_ref) = self.maps.genotypes.get(genotype_uniquename) {
            let mut genotype = genotype_ref.clone();
            let details_map = self.detail_map_of_cv_annotations(&genotype.cv_annotations);
            genotype.terms_by_termid = self.fill_term_map(&genotype.terms_by_termid);
            genotype.genes_by_uniquename = self.fill_gene_map(&genotype.genes_by_uniquename);
            genotype.references_by_uniquename =
                self.fill_reference_map(&genotype.references_by_uniquename);
            make_cv_summaries(&mut genotype, &self.config, &self.maps.children_by_termid,
                              false, false, &self.maps.genes, &self.maps.genotypes,
                              &details_map);
            genotype.annotation_details = details_map;
            Some(genotype)
        } else {
            None
        }
    }

    fn term_details_helper(&self, term_ref: &TermDetails) -> Option<TermDetails> {
        let mut term = term_ref.clone();
        let details_map = self.detail_map_of_cv_annotations(&term.cv_annotations);
        term.terms_by_termid = self.fill_term_map(&term.terms_by_termid);
        term.genes_by_uniquename = self.fill_gene_map(&term.genes_by_uniquename);
        term.references_by_uniquename =
            self.fill_reference_map(&term.references_by_uniquename);
        make_cv_summaries(&mut term, &self.config, &self.maps.children_by_termid,
                          true, true, &self.maps.genes, &self.maps.genotypes,
                          &details_map);
        term.annotation_details = details_map;
        Some(term)
    }

    pub fn get_term_details(&self, termid: &str) -> Option<TermDetails> {
        let termid =
            if let Some(real_termid) = self.secondary_identifiers_map.get(termid) {
                real_termid
            } else {
                termid
            };
        if let Some(term_ref) = self.maps.terms.get(termid) {
            self.term_details_helper(term_ref)
        } else {
            None
        }
    }

    pub fn get_reference_details(&self, reference_uniquename: &str) -> Option<ReferenceDetails> {
        if let Some(reference_ref) = self.maps.references.get(reference_uniquename) {
            let mut reference = reference_ref.clone();
            let details_map = self.detail_map_of_cv_annotations(&reference.cv_annotations);
            reference.terms_by_termid = self.fill_term_map(&reference.terms_by_termid);
            reference.genes_by_uniquename = self.fill_gene_map(&reference.genes_by_uniquename);
            make_cv_summaries(&mut reference, &self.config, &self.maps.children_by_termid,
                              true, true, &self.maps.genes, &self.maps.genotypes,
                              &details_map);
            reference.annotation_details = details_map;
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

    pub fn get_chr_details(&self, chr_name: &str) -> Option<&ChromosomeDetails> {
        self.maps.chromosomes.get(chr_name)
    }

    pub fn seq_feature_page_features(&self) -> Vec<FeatureShort> {
        self.maps.seq_feature_page_features.clone()
    }
}
