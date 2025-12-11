
use std::process;
use std::fs::File;
use std::io::{Read, BufReader};

use std::sync::{Mutex, Arc, RwLock};

use zstd::stream::Decoder;
use rusqlite::Connection;

use std::collections::{BTreeSet, HashMap, HashSet};

use crate::data_types::{APIAlleleDetails, APIGeneSummary, APIGenotypeAnnotation, APIMaps, AlleleDetails, AlleleShort, ChromosomeDetails, DataLookup, ExtPart, ExtRange, FeatureShort, GeneAndGeneProduct, GeneDetails, GeneQueryData, GeneShort, GeneShortOptionMap, GenotypeDetails, GoCamSummary, GoCamId, IdGeneSubsetMap, IdOntAnnotationDetailMap, InteractionType, OntAnnotationDetail, OntAnnotationId, OntAnnotationMap, Ploidiness, ProteinViewData, ProteinViewType, ReferenceDetails, ReferenceShort, ReferenceShortOptionMap, TermDetails, TermShort, TermShortOptionMap, Throughput, TranscriptDetailsOptionMap, WithFromValue};

use crate::sort_annotations::sort_cv_annotation_details;
use crate::web::config::{Config, TermAndName};
use crate::api::query::{QueryExpressionFilter, SingleOrMultiLocus, TargetOfType};
use crate::web::cv_summary::make_cv_summaries;

use crate::types::{TermId, GeneUniquename, AlleleUniquename,
                   GenotypeDisplayUniquename, GenotypeUniquename, ReferenceUniquename};

use flexstr::{SharedStr as FlexStr, shared_str as flex_str, shared_fmt as flex_fmt, ToSharedStr};

pub struct APIData {
    config: Config,
    maps: APIMaps,
    maps_database: APIMapsDatabase,
}

impl APIData {
    fn get_termid_genotype_annotation(&self, termid: &TermId)
           -> Option<Arc<Vec<APIGenotypeAnnotation>>>
    {
        self.maps_database.get_termid_genotype_annotation(termid)
    }
}

impl DataLookup for APIData {
    fn get_term(&self, termid: &TermId) -> Option<Arc<TermDetails>> {
        self.maps_database.get_term(termid)
    }

    fn get_gene(&self, gene_uniquename: &GeneUniquename) -> Option<Arc<GeneDetails>> {
        self.maps_database.get_gene(gene_uniquename)
    }

    fn get_allele(&self, allele_uniquename: &AlleleUniquename) -> Option<Arc<AlleleDetails>> {
        self.maps_database.get_allele(allele_uniquename)
    }

    fn get_reference(&self, reference_uniquename: &str)
           -> Option<Arc<ReferenceDetails>>
    {
        self.maps_database.get_reference(reference_uniquename)
    }

    fn get_genotype(&self, genotype_display_uniquename: &GenotypeDisplayUniquename)
           -> Option<Arc<GenotypeDetails>>
    {
        self.maps_database.get_genotype(genotype_display_uniquename)
    }

    fn get_annotation_detail(&self, annotation_id: OntAnnotationId)
           -> Option<Arc<OntAnnotationDetail>> {
        self.maps_database.get_annotation_detail(annotation_id)
    }
}

pub struct APIMapsDatabase {
    api_maps_database_conn: Arc<Mutex<Connection>>,
    genotype_cache: RwLock<HashMap<GenotypeUniquename, Arc<GenotypeDetails>>>,
    term_cache: RwLock<HashMap<TermId, Arc<TermDetails>>>,
    reference_cache: RwLock<HashMap<ReferenceUniquename, Arc<ReferenceDetails>>>,
    gene_cache: RwLock<HashMap<GeneUniquename, Arc<GeneDetails>>>,
    allele_cache: RwLock<HashMap<AlleleUniquename, Arc<AlleleDetails>>>,
    termid_genotype_annotation_cache: RwLock<HashMap<TermId, Arc<Vec<APIGenotypeAnnotation>>>>,
    annotation_detail_cache: RwLock<HashMap<OntAnnotationId, Arc<OntAnnotationDetail>>>,
}

impl APIMapsDatabase {
    pub fn new(api_maps_database_conn: Connection) -> APIMapsDatabase {
        APIMapsDatabase {
            api_maps_database_conn: Arc::new(Mutex::new(api_maps_database_conn)),
            genotype_cache: RwLock::new(HashMap::new()),
            term_cache: RwLock::new(HashMap::new()),
            reference_cache: RwLock::new(HashMap::new()),
            gene_cache: RwLock::new(HashMap::new()),
            allele_cache: RwLock::new(HashMap::new()),
            termid_genotype_annotation_cache: RwLock::new(HashMap::new()),
            annotation_detail_cache: RwLock::new(HashMap::new()),
        }
    }

    pub fn get_genotype(&self, genotype_display_uniquename: &GenotypeDisplayUniquename)
           -> Option<Arc<GenotypeDetails>>
    {
        let mut cache = self.genotype_cache.write().unwrap();

        if !cache.contains_key(genotype_display_uniquename) {
            let conn = self.api_maps_database_conn.lock().unwrap();

            let mut stmt = conn.prepare("SELECT data FROM genotypes WHERE id = :id").unwrap();

            let mut genotypes =
                stmt.query_map(&[(":id", genotype_display_uniquename.as_ref())],
                               |row| {
                                   let json: String = row.get(0)?;
                                   let genotype_details: GenotypeDetails =
                                       serde_json::from_str(&json).unwrap();
                                   Ok(genotype_details)
                               }).unwrap();

            let result_genotype = genotypes.next();

            let maybe_genotype = result_genotype.map(|g| g.unwrap());

            if let Some(genotype_details) = maybe_genotype {
                cache.insert(genotype_display_uniquename.clone(), Arc::new(genotype_details));
            }
        }

        let genotype_value = cache.get(genotype_display_uniquename);

        genotype_value.map(|g| g.to_owned())
    }

    pub fn get_term(&self, termid: &TermId)
           -> Option<Arc<TermDetails>>
    {
        let mut cache = self.term_cache.write().unwrap();

        if !cache.contains_key(termid) {
            let conn = self.api_maps_database_conn.lock().unwrap();

            let mut stmt = conn.prepare("SELECT data FROM terms WHERE id = :id").unwrap();

            let mut terms =
                stmt.query_map(&[(":id", termid.as_ref())],
                               |row| {
                                   let json: String = row.get(0)?;
                                   let term_details: TermDetails =
                                       serde_json::from_str(&json).unwrap();
                                   Ok(term_details)
                               }).unwrap();

            let result_term = terms.next();

            let maybe_term = result_term.map(|g| g.unwrap());

            if let Some(term_details) = maybe_term {
                cache.insert(termid.clone(), Arc::new(term_details));
            }
        }

        let term_value = cache.get(termid);

        term_value.map(|t| t.to_owned())
    }

    pub fn get_reference(&self, reference_uniquename: &str)
           -> Option<Arc<ReferenceDetails>>
    {
        let mut cache = self.reference_cache.write().unwrap();

        if !cache.contains_key(reference_uniquename) {
            let conn = self.api_maps_database_conn.lock().unwrap();

            let mut stmt = conn.prepare("SELECT data FROM refs WHERE id = :id").unwrap();

            let mut references =
                stmt.query_map(&[(":id", reference_uniquename)],
                               |row| {
                                   let json: String = row.get(0)?;
                                   let reference_details: ReferenceDetails =
                                       serde_json::from_str(&json).unwrap();
                                   Ok(reference_details)
                               }).unwrap();

            let result_reference = references.next();

            let maybe_reference = result_reference.map(|g| g.unwrap());

            if let Some(reference_details) = maybe_reference {
                cache.insert(reference_uniquename.to_shared_str(),
                             Arc::new(reference_details));
            }
        }

        let reference_value = cache.get(reference_uniquename);

        reference_value.map(|t| t.to_owned())
    }

    pub fn get_gene(&self, gene_uniquename: &GeneUniquename)
           -> Option<Arc<GeneDetails>>
    {
        let mut cache = self.gene_cache.write().unwrap();

        if !cache.contains_key(gene_uniquename) {
            let conn = self.api_maps_database_conn.lock().unwrap();

            let mut stmt = conn.prepare("SELECT data FROM genes WHERE id = :id").unwrap();

            let mut genes =
                stmt.query_map(&[(":id", gene_uniquename.as_ref())],
                               |row| {
                                   let json: String = row.get(0)?;
                                   let gene_details: GeneDetails =
                                       serde_json::from_str(&json).unwrap();
                                   Ok(gene_details)
                               }).unwrap();

            let result_gene = genes.next();

            let maybe_gene = result_gene.map(|g| g.unwrap());

            if let Some(gene_details) = maybe_gene {
                cache.insert(gene_uniquename.clone(), Arc::new(gene_details));
            }
        }

        let gene_value = cache.get(gene_uniquename);

        gene_value.map(|t| t.to_owned())
    }

    pub fn get_allele(&self, allele_uniquename: &AlleleUniquename)
           -> Option<Arc<AlleleDetails>>
    {
        let mut cache = self.allele_cache.write().unwrap();

        if !cache.contains_key(allele_uniquename) {
            let conn = self.api_maps_database_conn.lock().unwrap();

            let mut stmt = conn.prepare("SELECT data FROM alleles WHERE id = :id").unwrap();

            let mut alleles =
                stmt.query_map(&[(":id", allele_uniquename.as_ref())],
                               |row| {
                                   let json: String = row.get(0)?;
                                   let allele_details: AlleleDetails =
                                       serde_json::from_str(&json).unwrap();
                                   Ok(allele_details)
                               }).unwrap();

            let result_allele = alleles.next();

            let maybe_allele = result_allele.map(|g| g.unwrap());

            if let Some(allele_details) = maybe_allele {
                cache.insert(allele_uniquename.clone(), Arc::new(allele_details));
            }
        }

        let allele_value = cache.get(allele_uniquename);

        allele_value.map(|t| t.to_owned())
    }

    pub fn get_termid_genotype_annotation(&self, termid: &TermId)
           -> Option<Arc<Vec<APIGenotypeAnnotation>>>
    {
        let mut cache = self.termid_genotype_annotation_cache.write().unwrap();

        if !cache.contains_key(termid) {
            let conn = self.api_maps_database_conn.lock().unwrap();

            let mut stmt = conn.prepare("SELECT data FROM termid_genotype_annotations WHERE id = :id").unwrap();

            let mut termid_genotype_annotations =
                stmt.query_map(&[(":id", termid.as_ref())],
                               |row| {
                                   let json: String = row.get(0)?;
                                   let termid_genotype_annotations: Vec<APIGenotypeAnnotation> =
                                       serde_json::from_str(&json).unwrap();
                                   Ok(termid_genotype_annotations)
                               }).unwrap();

            let result_termid_genotype_annotations = termid_genotype_annotations.next();

            let maybe_termid_genotype_annotations =
                result_termid_genotype_annotations.map(|v| v.unwrap());

            if let Some(termid_genotype_annotations) = maybe_termid_genotype_annotations {
                cache.insert(termid.clone(), Arc::new(termid_genotype_annotations));
            }
        }

        let term_value = cache.get(termid);

        term_value.map(|t| t.to_owned())
    }

    pub fn get_annotation_detail(&self, annotation_id: OntAnnotationId)
           -> Option<Arc<OntAnnotationDetail>>
    {
        let mut cache = self.annotation_detail_cache.write().unwrap();

        if let std::collections::hash_map::Entry::Vacant(e) = cache.entry(annotation_id) {
            let conn = self.api_maps_database_conn.lock().unwrap();

            let mut stmt = conn.prepare("SELECT data FROM annotation_detail WHERE id = :id").unwrap();

            let mut annotation_detail =
                stmt.query_map(&[(":id", &annotation_id)],
                               |row| {
                                   let json: String = row.get(0)?;
                                   let annotation_details: OntAnnotationDetail =
                                       serde_json::from_str(&json).unwrap();
                                   Ok(annotation_details)
                               }).unwrap();

            let result_annotation_detail = annotation_detail.next();

            let maybe_annotation_detail =
                result_annotation_detail.map(|v| v.unwrap());

            if let Some(annotation_detail) = maybe_annotation_detail {
                e.insert(Arc::new(annotation_detail));
            }
        }

        let term_value = cache.get(&annotation_id);

        term_value.map(|t| t.to_owned())
    }

}

pub fn api_maps_from_file(search_maps_file_name: &str) -> APIMaps
{
    let file = match File::open(search_maps_file_name) {
        Ok(file) => file,
        Err(err) => {
            eprintln!("Failed to read {}: {}", search_maps_file_name, err);
            process::exit(1);
        }
    };

    let reader = BufReader::new(file);
    let mut decoder = Decoder::new(reader).unwrap();

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
    pub fn new(config: &Config, maps_database_conn: Connection, mut maps: APIMaps)
               ->APIData
    {
        let mut new_entries: IdGeneSubsetMap = HashMap::new();

        let prefixes_to_remove: Vec<FlexStr> =
            config.server.subsets.prefixes_to_remove
            .iter().map(|prefix| prefix.clone() + ":").collect();

        // strip prefixes and add to map
        for (subset_name, subset_details) in &maps.gene_subsets {
            for prefix in &prefixes_to_remove {
                if subset_name.starts_with(prefix.as_ref()) {
                    let new_subset_name = subset_name[prefix.len()..].to_shared_str();
                    new_entries.insert(new_subset_name, subset_details.clone());
                }
            }
        }

        maps.gene_subsets.extend(new_entries);

        let maps_database = APIMapsDatabase::new(maps_database_conn);

        APIData {
            config: config.clone(),
            maps,
            maps_database,
        }
    }

    pub fn get_config(&self) -> &Config {
        &self.config
    }

    pub fn get_maps(&self) -> &APIMaps {
        &self.maps
    }

    pub fn gene_uniquename_of_id(&self, id: &FlexStr) -> Option<GeneUniquename> {
        if self.maps.gene_summaries.contains_key(id) {
            Some(id.clone())
        } else if let Some(gene_uniquename) =
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

    pub fn get_gene_summary(&self, gene_uniquename: &FlexStr) -> Option<&APIGeneSummary> {
        self.maps.gene_summaries.get(gene_uniquename)
    }

    pub fn get_gene_details(&self, gene_uniquename: &FlexStr) -> Option<Arc<GeneDetails>> {
        self.get_gene(gene_uniquename)
    }

    pub fn get_gene_query_data(&self, gene_uniquename: &FlexStr) -> Option<&GeneQueryData> {
        self.maps.gene_query_data_map.get(gene_uniquename)
    }

    pub fn genes_of_termid(&self, term_id: &FlexStr) -> Vec<GeneUniquename> {
        match self.maps.termid_genes.get(term_id) {
            Some(gene_uniquenames) => {
                gene_uniquenames.iter().cloned().collect::<Vec<_>>()
            },
            None => vec![],
        }
    }

    pub fn genes_of_subset(&self, search_name: &FlexStr) -> Vec<GeneUniquename> {
        if search_name.starts_with('!') || search_name.ends_with('*') {
            let mut trimmed_search_name = search_name.to_string();
            let invert_search = search_name.starts_with('!');
            let wildcard = search_name.ends_with('*');
            if invert_search {
                trimmed_search_name.remove(0);
            }
            if wildcard {
                trimmed_search_name.pop();
            }
            let trimmed_search_name = trimmed_search_name.to_shared_str();
            let mut genes = HashSet::new();
            for (subset_name, subset_details) in &self.maps.gene_subsets {
                let name_matches =
                    wildcard && subset_name.starts_with(trimmed_search_name.as_ref()) ||
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

    pub fn genes_of_genotypes(&self, term_id: &FlexStr,
                              single_or_multi_locus: &SingleOrMultiLocus,
                              query_ploidiness: &Ploidiness,
                              expression_filter: &Option<QueryExpressionFilter>,
                              conditions_filter: &HashSet<TermAndName>,
                              excluded_conditions_filter: &HashSet<TermAndName>)
                              -> Vec<GeneUniquename>
    {
        if let Some(annotations) = self.get_termid_genotype_annotation(term_id) {
            let mut genes: HashSet<GeneUniquename> = HashSet::new();
            for annotation in annotations.as_ref() {

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
                            } else if *expression_filter == QueryExpressionFilter::WtOverexpressed {
                                if let Some(ref expression) = *expression {
                                    expression == "Overexpression" &&
                                        allele_type == "wild_type"
                                } else {
                                    false
                                }
                            } else {
                                false
                            }
                        } else {
                            true
                        }
                    };

                let condition_matches =
                    |conditions: &HashSet<TermAndName>, test_conditions_filter: &HashSet<TermAndName>| {
                        if test_conditions_filter.is_empty() {
                            // don't filter in this case
                            return true;
                        }

                        for filter_condition in test_conditions_filter {
                            if !conditions.contains(filter_condition) {
                                return false
                            }
                        }

                        true
                    };

                if !condition_matches(&annotation.conditions, conditions_filter) ||
                    !excluded_conditions_filter.is_empty() &&
                    condition_matches(&annotation.conditions, excluded_conditions_filter) {
                        continue;
                    }


                let mut add_genotype_genes = false;

                if annotation.is_multi && add_multi && (*query_ploidiness == annotation.ploidiness || *query_ploidiness == Ploidiness::Any) {
                    add_genotype_genes = true;
                }

                if !annotation.is_multi && add_single &&
                    ((*query_ploidiness != Ploidiness::Haploid &&
                        annotation.ploidiness == Ploidiness::Diploid) || (annotation.ploidiness == Ploidiness::Haploid &&
                         (*query_ploidiness == Ploidiness::Any ||
                          (*query_ploidiness == Ploidiness::Haploid &&
                           expression_matches(&annotation.alleles[0]))))) {
                    add_genotype_genes = true;
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

    // return all genes where the predicate returns true
    pub fn filter_genes(&self, p: & dyn Fn(&APIGeneSummary) -> bool)
                        -> Vec<GeneUniquename>
    {
        self.maps.gene_summaries.values()
            .filter(|summ| p(summ))
            .map(|summ| summ.uniquename.clone())
            .collect()
    }

    // return protein coding genes where the predicate returns true
    pub fn filter_proteins(&self, p: & dyn Fn(&APIGeneSummary) -> bool)
                           -> Vec<GeneUniquename>
    {
        self.maps.gene_summaries.values()
            .filter(|summ| summ.feature_type == "mRNA gene" && p(summ))
            .map(|summ| summ.uniquename.clone())
            .collect()
    }

    fn strip_db_prefix(&self, uniquename: &FlexStr) -> FlexStr {
        if uniquename.starts_with("PomBase:SP") {
            uniquename[8..].to_shared_str()
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
                    ExtRange::Gene(self.strip_db_prefix(uniquename))
                },
                WithFromValue::Transcript(transcript_uniquename) => {
                    ExtRange::Transcript(transcript_uniquename)
                }
                WithFromValue::Identifier(identifier) => {
                    ExtRange::Misc(identifier)
                }
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
        if let Some(ref evidence_code) = annotation.evidence
            && evidence_code == "IPI" &&
                (term_details.termid == "GO:0005515" ||
                 term_details.interesting_parent_ids.contains(&flex_str!("GO:0005515")) ||
                 term_details.termid == "GO:0003723" ||
                 term_details.interesting_parent_ids.contains(&flex_str!("GO:0003723")))
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
                // if there is ExtRange::GeneAndGeneProduct that mentions the with gene,
                // drop the with for this annotation
                if let Some(with_value) = first_with.clone() {
                    annotation.extension.iter_mut().for_each(|ext_part| {
                        if ext_part.rel_type_name == "has_input" {
                            match &mut ext_part.ext_range {
                                ExtRange::GeneProduct(range_termid) => {
                                    if let WithFromValue::Gene(gene_short) = &with_value {
                                        let gene_uniquename =
                                            self.strip_db_prefix(&gene_short.uniquename);
                                        let val = GeneAndGeneProduct {
                                            product: range_termid.clone(),
                                            gene_uniquename,
                                        };
                                        ext_part.ext_range = ExtRange::GeneAndGeneProduct(val);
                                        first_with = None;
                                    }
                                },
                                ExtRange::GeneAndGeneProduct(GeneAndGeneProduct { gene_uniquename, product: _ }) => {
                                    if let WithFromValue::Gene(ref with_gene) = with_value
                                        && *gene_uniquename == with_gene.uniquename {
                                           first_with = None;
                                        }
                                },
                                _ => (),
                            }
                        }
                    });
                }

                if let Some(with_value) = first_with {
                    annotation.extension.insert(0, self.get_with_extension(with_value));
                }
            }
    }

    fn get_gene_prod_extension(&self, prod_value: &FlexStr) -> ExtPart {
      let ext_range =
        if let Some(term_details) = self.get_term(prod_value) {
          if let Some(ref pombase_gene_id) = term_details.pombase_gene_id {
            let gene_and_product = GeneAndGeneProduct {
              gene_uniquename: pombase_gene_id.clone(),
              product: prod_value.to_owned(),
            };
            ExtRange::GeneAndGeneProduct(gene_and_product)
          } else {
            ExtRange::GeneProduct(prod_value.to_owned())
          }
        } else {
          ExtRange::GeneProduct(prod_value.to_owned())
        };
      ExtPart {
        rel_type_id: None,
        rel_type_name: "active_form".into(),
        rel_type_display_name: "active form".into(),
        ext_range,
      }
    }

    fn detail_map_of_cv_annotations(&self, ont_annotation_map: &OntAnnotationMap)
                                    -> IdOntAnnotationDetailMap
    {
        let mut details_map = HashMap::new();

        for term_annotations in ont_annotation_map.values() {
            for term_annotation in term_annotations {
                let termid = &term_annotation.term;
                if let Some(term_details) = self.get_term(termid) {
                    let annotations: Vec<OntAnnotationDetail> =
                        term_annotation.annotations
                        .iter()
                        .map(|annotation_detail_id| {
                            let mut annotation =
                                self.get_annotation_detail(*annotation_detail_id)
                                    .unwrap().as_ref().to_owned();

                           self.maybe_move_with(&term_details, &mut annotation);

                           if let Some(ref gene_product_form_id) = annotation.gene_product_form_id {
                               let gene_prod_extension = self.get_gene_prod_extension(gene_product_form_id);
                               annotation.extension.insert(0, gene_prod_extension);
                           }

                           annotation
                        })
                        .collect();

                    for annotation in annotations {
                        let annotation_detail_id = annotation.id;
                        details_map.insert(annotation_detail_id, annotation);
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
            if let Some(term_details) = self.get_term(termid) {
                let term_short = TermShort::from_term_details(&term_details);
                ret.insert(termid.clone(), Some(term_short));
            } else {
                eprintln!("WARNING missing short info for: {}", termid);
                ret.insert(termid.clone(), None);
            }
        }
        ret
    }

    fn fill_gene_map(&self, gene_map: &GeneShortOptionMap) -> GeneShortOptionMap {
        let mut ret = gene_map.clone();
        for gene_uniquename in gene_map.keys() {
            let gene_details_arc = self.get_gene(gene_uniquename).unwrap();
            let gene_details = gene_details_arc.as_ref();
            let gene_short = GeneShort::from_gene_details(gene_details);
            ret.insert(gene_uniquename.clone(), Some(gene_short));
        }
        ret
    }

    fn fill_transcript_map(&self, transcript_map: &TranscriptDetailsOptionMap) -> TranscriptDetailsOptionMap {
        let mut ret = transcript_map.clone();

        for transcript_uniquename in transcript_map.keys() {
            let transcript_details = &self.maps.transcripts[transcript_uniquename];
            ret.insert(transcript_uniquename.clone(), Some(transcript_details.to_owned()));
        }
        ret
    }

    fn fill_reference_map(&self, reference_map: &ReferenceShortOptionMap) -> ReferenceShortOptionMap {
        let mut ret = reference_map.clone();
        for reference_uniquename in reference_map.keys() {
            let reference_arc = &self.get_reference(reference_uniquename).unwrap();
            let reference_details = reference_arc.as_ref();
            let reference_short = ReferenceShort::from_reference_details(reference_details);
            ret.insert(reference_uniquename.clone(), Some(reference_short));
        }
        ret
    }

    // return a GeneDetails object with the term, genes and references maps filled in
    pub fn get_full_gene_details(&self, gene_uniquename: &str) -> Option<GeneDetails> {
        let gene_uniquename = gene_uniquename.to_shared_str();
        if let Some(gene_arc) = self.get_gene(&gene_uniquename) {
            let mut gene = gene_arc.as_ref().to_owned();
            let details_map = self.detail_map_of_cv_annotations(&gene.cv_annotations);
            gene.terms_by_termid = self.fill_term_map(&gene.terms_by_termid);
            gene.genes_by_uniquename = self.fill_gene_map(&gene.genes_by_uniquename);
            gene.transcripts_by_uniquename = self.fill_transcript_map(&gene.transcripts_by_uniquename);
            gene.references_by_uniquename =
                self.fill_reference_map(&gene.references_by_uniquename);
            sort_cv_annotation_details(&mut gene, &self.config,
                                       self,
                                       &details_map);
            make_cv_summaries(&mut gene, &self.config, &self.maps.children_by_termid,
                              false, true, self,
                              &details_map);

            gene.annotation_details = details_map;
            Some(gene)
        } else {
            None
        }
    }

    pub fn get_genotype_details(&self, genotype_uniquename: &str) -> Option<GenotypeDetails> {
        let genotype_uniquename = genotype_uniquename.to_shared_str();
        if let Some(genotype_ref) = self.get_genotype(&genotype_uniquename) {
            let mut genotype = genotype_ref.as_ref().to_owned();
            let details_map = self.detail_map_of_cv_annotations(&genotype.cv_annotations);
            genotype.terms_by_termid = self.fill_term_map(&genotype.terms_by_termid);
            genotype.genes_by_uniquename = self.fill_gene_map(&genotype.genes_by_uniquename);
            genotype.transcripts_by_uniquename = self.fill_transcript_map(&genotype.transcripts_by_uniquename);
            genotype.references_by_uniquename =
                self.fill_reference_map(&genotype.references_by_uniquename);
            sort_cv_annotation_details(&mut genotype, &self.config,
                                       self,
                                       &details_map);
            make_cv_summaries(&mut genotype, &self.config, &self.maps.children_by_termid,
                              false, false, self,
                              &details_map);
            genotype.annotation_details = details_map;
            Some(genotype)
        } else {
            None
        }
    }

    pub fn get_allele_details(&self, allele_uniquename: &str) -> Option<AlleleDetails> {
        let allele_uniquename = allele_uniquename.to_shared_str();

        if let Some(allele_arc) = self.get_allele(&allele_uniquename) {
            let mut allele_details = allele_arc.as_ref().to_owned();
            let mut alleles_by_uniquename = HashMap::new();

            for genotype_short in &allele_details.genotypes {
                for locus in &genotype_short.loci {
                    for expressed_allele in &locus.expressed_alleles {
                        let allele_uniquename = &expressed_allele.allele_uniquename;
                        let allele_short: AlleleShort =
                            self.get_allele(allele_uniquename).unwrap().as_ref().into();
                        alleles_by_uniquename.insert(allele_uniquename.clone(), allele_short);
                    }
                }
            }

            allele_details.alleles_by_uniquename = alleles_by_uniquename;

            Some(allele_details)
        } else {
            None
        }
    }

    fn term_details_helper(&self, term_ref: &TermDetails) -> Option<TermDetails> {
        let mut term = term_ref.clone();
        let details_map = self.detail_map_of_cv_annotations(&term.cv_annotations);
        term.terms_by_termid = self.fill_term_map(&term.terms_by_termid);
        term.genes_by_uniquename = self.fill_gene_map(&term.genes_by_uniquename);
        term.transcripts_by_uniquename = self.fill_transcript_map(&term.transcripts_by_uniquename);
        term.references_by_uniquename =
            self.fill_reference_map(&term.references_by_uniquename);
        sort_cv_annotation_details(&mut term, &self.config,
                                   self,
                                   &details_map);
        make_cv_summaries(&mut term, &self.config, &self.maps.children_by_termid,
                          true, true, self,
                          &details_map);
        term.annotation_details = details_map;
        Some(term)
    }

    pub fn get_term_details(&self, termid: &str) -> Option<TermDetails> {
        let termid = termid.to_shared_str();
        let termid =
            if let Some(real_termid) =
               self.maps.secondary_identifiers_map.get(&termid) {
                real_termid
            } else {
                &termid
            };
        if let Some(term_ref) = self.get_term(termid) {
            self.term_details_helper(&term_ref)
        } else {
            None
        }
    }

    pub fn get_reference_details(&self, reference_uniquename: &str) -> Option<ReferenceDetails> {
        let reference_uniquename = reference_uniquename.to_shared_str();
        if let Some(reference_ref) = self.get_reference(&reference_uniquename) {
            let mut reference = reference_ref.as_ref().to_owned();
            let details_map = self.detail_map_of_cv_annotations(&reference.cv_annotations);
            reference.terms_by_termid = self.fill_term_map(&reference.terms_by_termid);
            reference.genes_by_uniquename = self.fill_gene_map(&reference.genes_by_uniquename);
            reference.transcripts_by_uniquename = self.fill_transcript_map(&reference.transcripts_by_uniquename);
            sort_cv_annotation_details(&mut reference, &self.config,
                                       self,
                                       &details_map);
            make_cv_summaries(&mut reference, &self.config, &self.maps.children_by_termid,
                              true, true, self,
                              &details_map);
            reference.annotation_details = details_map;
            Some(reference)
        } else {
            None
        }
    }

    pub fn interactors_of_genes(&self, gene_uniquename: &GeneUniquename,
                                interaction_type_constraint: InteractionType,
                                throughput_constraint: &Option<Throughput>,
                                evidence_type_constraint: &Option<String>)
                                -> Vec<GeneUniquename> {
        if let Some(interactors) = self.maps.interactors_of_genes.get(gene_uniquename) {
            interactors.iter()
                .filter(|interactor| {
                    if interactor.interaction_type != interaction_type_constraint {
                        return false;
                    }
                    if let Some(throughput_constraint) = throughput_constraint {
                        if let Some(ref interactor_throughput) = interactor.throughput {
                            if interactor_throughput != throughput_constraint {
                                return false;
                            }
                        } else {
                            return false;
                        }
                    }
                    if let Some(evidence_type) = evidence_type_constraint
                        && interactor.evidence_type != *evidence_type {
                            return false;
                        }
                    true
                })
                .map(|interactor| interactor.interactor_uniquename.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect()
        } else {
            vec![]
        }
    }

    pub fn downstream_genes(&self, cv_name: &str, gene_uniquename: &GeneUniquename,
                             phase_term: &Option<TermId>)
       -> HashSet<GeneUniquename>
    {
        let key = flex_fmt!("{}--{}", cv_name, gene_uniquename);
        let Some(downstream_genes_by_term) = self.maps.downstream_genes.get(&key)
        else {
            return HashSet::new();
        };

        if let Some(phase_term) = phase_term {
            if let Some(downstream_genes) = downstream_genes_by_term.get(phase_term) {
                downstream_genes.to_owned()
            } else {
                HashSet::new()
            }
        } else {
            let mut ret_downstream_genes = HashSet::new();
            for downstream_genes in downstream_genes_by_term.values() {
                ret_downstream_genes.extend(downstream_genes.to_owned());
            }
            ret_downstream_genes
        }
    }

    pub fn genes_targeting(&self, gene_uniquename: &GeneUniquename, target_of_type: TargetOfType)
      -> HashSet<GeneUniquename>
    {
        let Some(gene_arc) = self.get_gene(gene_uniquename)
        else {
            return HashSet::new();
        };

        let gene_details = gene_arc.as_ref();

        let mut ret_hashset = HashSet::new();

        let target_of_config = &self.config.target_of_config;
        let unknown = &flex_str!("unknown");

        for target_of_annotation in &gene_details.target_of_annotations {
            let ontology_label =
                target_of_config.ontology_labels.get(&target_of_annotation.ontology_name)
                .unwrap_or(unknown);

            if (target_of_type == TargetOfType::GO || target_of_type == TargetOfType::All) &&
                ontology_label == "GO" {
                ret_hashset.insert(target_of_annotation.gene.clone());
            }
            if (target_of_type == TargetOfType::Phenotype || target_of_type == TargetOfType::All) &&
                ontology_label != "GO" {
                ret_hashset.insert(target_of_annotation.gene.clone());
            }
        }

        ret_hashset
    }

    pub fn get_chr_details(&self, chr_name: &FlexStr) -> Option<&ChromosomeDetails> {
        self.maps.chromosomes.get(chr_name)
    }

    pub fn seq_feature_page_features(&self) -> Vec<FeatureShort> {
        self.maps.seq_feature_page_features.clone()
    }

    pub fn get_protein_features_of_gene(&self, scope: ProteinViewType,
                                        gene_uniquename: &str)
        -> Option<ProteinViewData>
    {
        let gene_uniquename = FlexStr::from(gene_uniquename);

        let data = self.maps.protein_view_data.get(&gene_uniquename)?;

        let filtered_tracks = data.tracks
            .iter()
            .filter(|track| {
                let prot_feat_conf = &self.config.protein_feature_view;
                match scope {
                    ProteinViewType::Widget =>
                        prot_feat_conf.widget_tracks.contains(&track.name),
                    ProteinViewType::DomainsAndFeatures => {
                        let lc_track_name = track.name.to_ascii_lowercase();
                        (lc_track_name.starts_with("pfam") &&
                         lc_track_name != "pfam domains") ||
                        prot_feat_conf.domains_and_features_tracks.contains(&track.name)
                    },
                    ProteinViewType::Modifications =>
                        prot_feat_conf.modification_section_tracks.contains(&track.name),
                    ProteinViewType::Full =>
                        !prot_feat_conf.full_display_excluded.contains(&track.name),

                }
            })
            .map(Clone::clone)
            .collect();

        Some(ProteinViewData {
            sequence: data.sequence.clone(),
            tracks: filtered_tracks
        })
    }

    pub fn get_gocam_data_of_gene(&self, gene_uniquename: &str)
        -> Option<HashSet<GoCamId>>
    {
        let gene_uniquename = FlexStr::from(gene_uniquename);

        self.maps.gocam_data_by_gene.get(&gene_uniquename).map(|d| d.to_owned())
    }

    pub fn get_gocam_details_by_id(&self, gocam_id: &str)
        -> Option<GoCamSummary>
    {
        self.maps.gocam_data_by_gocam_id.get(gocam_id).map(|d| d.to_owned())
    }

    pub fn get_all_gocam_data(&self)
        -> HashMap<GoCamId, GoCamSummary>
    {
        self.maps.gocam_data_by_gocam_id.clone()
    }
}
