use std::rc::Rc;
use std::collections::hash_map::HashMap;
use std::collections::HashSet;
use std::borrow::Borrow;
use std::cmp::Ordering;

use regex::Regex;

use db::*;
use web::data::*;
use web::config::*;

include!(concat!(env!("OUT_DIR"), "/config_serde.rs"));

fn make_organism_short(rc_organism: &Rc<Organism>) -> OrganismShort {
    OrganismShort {
        genus: rc_organism.genus.clone(),
        species: rc_organism.species.clone(),
    }
}

#[derive(Clone)]
pub struct AlleleAndExpression {
    allele_uniquename: String,
    expression: Option<String>,
}

pub struct WebDataBuild<'a> {
    raw: &'a Raw,
    config: &'a Config,

    genes: UniquenameGeneMap,
    transcripts: UniquenameTranscriptMap,
    genotypes: UniquenameGenotypeMap,
    alleles: UniquenameAlleleMap,
    terms: HashMap<TermId, TermDetails>,
    references: IdReferenceMap,
    all_ont_annotations: HashMap<TermId, Vec<Rc<OntAnnotationDetail>>>,
    all_not_ont_annotations: HashMap<TermId, Vec<Rc<OntAnnotationDetail>>>,

    genes_of_transcripts: HashMap<String, String>,
    transcripts_of_polypeptides: HashMap<String, String>,
    genes_of_alleles: HashMap<String, String>,
    alleles_of_genotypes: HashMap<String, Vec<AlleleAndExpression>>,
    // gene_uniquename vs transcript_type_name:
    transcript_type_of_genes: HashMap<String, String>,

    // a map from IDs of terms from the "PomBase annotation extension terms" cv
    // to a Vec of the details of each of the extension
    parts_of_extensions: HashMap<TermId, Vec<ExtPart>>,

    base_term_of_extensions: HashMap<TermId, TermId>,
}

fn get_maps() ->
    (HashMap<String, ReferenceShortMap>,
     HashMap<String, GeneShortMap>,
     HashMap<String, GenotypeShortMap>,
     HashMap<String, AlleleShortMap>,
     HashMap<GeneUniquename, TermShortMap>)
{
    (HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new())
}

fn get_feat_rel_expression(feature_relationship: &FeatureRelationship) -> Option<String> {
    for prop in feature_relationship.feature_relationshipprops.borrow().iter() {
        if prop.prop_type.name == "expression" {
            return prop.value.clone();
        }
    }

    None
}

fn is_gene_type(feature_type_name: &str) -> bool {
    feature_type_name == "gene" || feature_type_name == "pseudogene"
}

impl <'a> WebDataBuild<'a> {
    pub fn new(raw: &'a Raw, config: &'a Config) -> WebDataBuild<'a> {
        WebDataBuild {
            raw: raw,
            config: config,

            genes: HashMap::new(),
            transcripts: HashMap::new(),
            genotypes: HashMap::new(),
            alleles: HashMap::new(),
            terms: HashMap::new(),
            references: HashMap::new(),
            all_ont_annotations: HashMap::new(),
            all_not_ont_annotations: HashMap::new(),

            genes_of_transcripts: HashMap::new(),
            transcripts_of_polypeptides: HashMap::new(),
            genes_of_alleles: HashMap::new(),
            alleles_of_genotypes: HashMap::new(),
            transcript_type_of_genes: HashMap::new(),

            parts_of_extensions: HashMap::new(),

            base_term_of_extensions: HashMap::new(),
        }
    }

    fn add_ref_to_hash(&self,
                       seen_references: &mut HashMap<String, ReferenceShortMap>,
                       identifier: String,
                       maybe_reference_uniquename: Option<ReferenceUniquename>) {
        if let Some(reference_uniquename) = maybe_reference_uniquename {
            if let Some(reference_short) = self.make_reference_short(&reference_uniquename) {
                seen_references
                    .entry(identifier.clone())
                    .or_insert(HashMap::new())
                    .insert(reference_uniquename.clone(),
                            reference_short);
            }
        }
    }

    fn add_gene_to_hash(&self,
                        seen_genes: &mut HashMap<String, GeneShortMap>,
                        identifier: String,
                        other_gene_uniquename: GeneUniquename) {
        seen_genes
            .entry(identifier)
            .or_insert(HashMap::new())
            .insert(other_gene_uniquename.clone(),
                    self.make_gene_short(&other_gene_uniquename));
    }

    fn add_genotype_to_hash(&self,
                            seen_genotypes: &mut HashMap<String, GenotypeShortMap>,
                            identifier: String,
                            genotype_uniquename: &GenotypeUniquename) {
        seen_genotypes
            .entry(identifier)
            .or_insert(HashMap::new())
            .insert(genotype_uniquename.clone(),
                    self.make_genotype_short(genotype_uniquename));
    }

    fn add_allele_to_hash(&self,
                          seen_alleles: &mut HashMap<String, AlleleShortMap>,
                          identifier: String,
                          allele_uniquename: AlleleUniquename) -> AlleleShort {
        let allele_short = self.make_allele_short(&allele_uniquename);
        seen_alleles
            .entry(identifier)
            .or_insert(HashMap::new())
            .insert(allele_uniquename, allele_short.clone());
        allele_short
    }

    fn add_term_to_hash(&self,
                        seen_terms: &mut HashMap<TermId, TermShortMap>,
                        identifier: String,
                        other_termid: TermId) {
        seen_terms
            .entry(identifier)
            .or_insert(HashMap::new())
            .insert(other_termid.clone(),
                    self.make_term_short(&other_termid));
    }

    fn get_gene<'b>(&'b self, gene_uniquename: &'b str) -> &'b GeneDetails {
        if let Some(gene_details) = self.genes.get(gene_uniquename) {
            gene_details
        } else {
            panic!("can't find GeneDetails for gene uniquename {}", gene_uniquename)
        }
    }

    fn get_gene_mut<'b>(&'b mut self, gene_uniquename: &'b str) -> &'b mut GeneDetails {
        if let Some(gene_details) = self.genes.get_mut(gene_uniquename) {
            gene_details
        } else {
            panic!("can't find GeneDetails for gene uniquename {}", gene_uniquename)
        }
    }

    fn make_gene_short(&self, gene_uniquename: &str) -> GeneShort {
        let gene_details = self.get_gene(&gene_uniquename);
        GeneShort {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
        }
    }

    fn make_gene_summary(&self, gene_uniquename: &str) -> GeneSummary {
        let gene_details = self.get_gene(&gene_uniquename);
        GeneSummary {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
            synonyms: gene_details.synonyms.clone(),
            feature_type: gene_details.feature_type.clone(),
            organism: gene_details.organism.clone(),
            location: gene_details.location.clone(),
        }
    }

    fn make_reference_short(&self, reference_uniquename: &str) -> Option<ReferenceShort> {
        if reference_uniquename == "null" {
            None
        } else {
            let reference_details = self.references.get(reference_uniquename).unwrap();

            let reference_short =
                ReferenceShort {
                    uniquename: String::from(reference_uniquename),
                    title: reference_details.title.clone(),
                    citation: reference_details.citation.clone(),
                    publication_year: reference_details.publication_year.clone(),
                    authors: reference_details.authors.clone(),
                    authors_abbrev: reference_details.authors_abbrev.clone(),
                    gene_count: reference_details.genes_by_uniquename.keys().len(),
                    genotype_count: reference_details.genotypes_by_uniquename.keys().len(),
                };

            Some(reference_short)
        }
    }

    fn make_term_short(&self, termid: &str) -> TermShort {
        if let Some(term_details) = self.terms.get(termid) {
            TermShort {
                name: term_details.name.clone(),
                cv_name: term_details.cv_name.clone(),
                interesting_parents: term_details.interesting_parents.clone(),
                termid: term_details.termid.clone(),
                is_obsolete: term_details.is_obsolete,
                gene_count: term_details.genes_by_uniquename.keys().len(),
                genotype_count: term_details.genotypes_by_uniquename.keys().len(),
            }
        } else {
            panic!("can't find TermDetails for termid: {}", termid)
        }
    }

    fn add_characterisation_status(&mut self, gene_uniquename: &String, cvterm_name: &String) {
        let mut gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        gene_details.characterisation_status = Some(cvterm_name.clone());
    }

    fn add_gene_product(&mut self, gene_uniquename: &String, product: &String) {
        let mut gene_details = self.get_gene_mut(gene_uniquename);
        gene_details.product = Some(product.clone());
    }

    fn add_name_description(&mut self, gene_uniquename: &str, name_description: &str) {
        let mut gene_details = self.get_gene_mut(gene_uniquename);
        gene_details.name_descriptions.push(name_description.into());
    }

    fn add_annotation(&mut self, cvterm: &Cvterm, is_not: bool,
                      annotation_template: OntAnnotationDetail) {
        let termid =
            match self.base_term_of_extensions.get(&cvterm.termid()) {
                Some(base_termid) => base_termid.clone(),
                None => cvterm.termid(),
            };

        let extension_parts =
            match self.parts_of_extensions.get(&cvterm.termid()) {
                Some(parts) => parts.clone(),
                None => vec![],
            };

        let mut new_extension = extension_parts.clone();
        let mut existing_extensions = annotation_template.extension.clone();
        new_extension.append(&mut existing_extensions);

        let ont_annotation_detail =
            OntAnnotationDetail {
                extension: new_extension,
                .. annotation_template
            };

        let annotation_map = if is_not {
            &mut self.all_not_ont_annotations
        } else {
            &mut self.all_ont_annotations
        };

        let entry = annotation_map.entry(termid.clone());
        entry.or_insert(
            vec![]
        ).push(Rc::new(ont_annotation_detail));
    }

    fn process_references(&mut self) {
        for rc_publication in &self.raw.publications {
            let reference_uniquename = &rc_publication.uniquename;

            let mut pubmed_authors: Option<String> = None;
            let mut pubmed_publication_date: Option<String> = None;
            let mut pubmed_abstract: Option<String> = None;

            for prop in rc_publication.publicationprops.borrow().iter() {
                match &prop.prop_type.name as &str {
                    "pubmed_publication_date" =>
                        pubmed_publication_date = Some(prop.value.clone()),
                    "pubmed_authors" =>
                        pubmed_authors = Some(prop.value.clone()),
                    "pubmed_abstract" =>
                        pubmed_abstract = Some(prop.value.clone()),
                    _ => ()
                }
            }

            let mut authors_abbrev = None;
            let mut publication_year = None;

            if let Some(authors) = pubmed_authors.clone() {
                if authors.contains(",") {
                    let author_re = Regex::new(r"^(?P<f>[^,]+),.*$").unwrap();
                    authors_abbrev = Some(author_re.replace_all(&authors, "$f et al."));
                } else {
                    authors_abbrev = Some(authors.clone());
                }
            }

            if let Some(publication_date) = pubmed_publication_date.clone() {
                let date_re = Regex::new(r"^(.* )?(?P<y>\d\d\d\d)$").unwrap();
                publication_year = Some(date_re.replace_all(&publication_date, "$y"));
            }

            self.references.insert(reference_uniquename.clone(),
                                   ReferenceDetails {
                                       uniquename: reference_uniquename.clone(),
                                       title: rc_publication.title.clone(),
                                       citation: rc_publication.miniref.clone(),
                                       pubmed_abstract: pubmed_abstract.clone(),
                                       authors: pubmed_authors.clone(),
                                       authors_abbrev: authors_abbrev,
                                       pubmed_publication_date: pubmed_publication_date.clone(),
                                       publication_year: publication_year,
                                       cv_annotations: HashMap::new(),
                                       physical_interactions: vec![],
                                       genetic_interactions: vec![],
                                       ortholog_annotations: vec![],
                                       paralog_annotations: vec![],
                                       genes_by_uniquename: HashMap::new(),
                                       genotypes_by_uniquename: HashMap::new(),
                                       alleles_by_uniquename: HashMap::new(),
                                       terms_by_termid: HashMap::new(),
                                   });
        }
    }

    fn make_feature_rel_maps(&mut self) {
        for feature_rel in self.raw.feature_relationships.iter() {
            let subject_type_name = &feature_rel.subject.feat_type.name;
            let rel_name = &feature_rel.rel_type.name;
            let object_type_name = &feature_rel.object.feat_type.name;
            let subject_uniquename = &feature_rel.subject.uniquename;
            let object_uniquename = &feature_rel.object.uniquename;

            if TRANSCRIPT_FEATURE_TYPES.contains(&subject_type_name.as_str()) &&
                rel_name == "part_of" &&
                (object_type_name == "gene" || object_type_name == "pseudogene") {
                    self.genes_of_transcripts.insert(subject_uniquename.clone(),
                                                     object_uniquename.clone());
                    self.transcript_type_of_genes.insert(object_uniquename.clone(),
                                                         subject_type_name.clone());
                    continue;
                }
            if subject_type_name == "polypeptide" &&
                rel_name == "derives_from" &&
                object_type_name == "mRNA" {
                    self.transcripts_of_polypeptides.insert(subject_uniquename.clone(),
                                                            object_uniquename.clone());
                    continue;
                }
            if subject_type_name == "allele" {
                if feature_rel.rel_type.name == "instance_of" &&
                    (object_type_name == "gene" || object_type_name == "pseudogene") {
                        self.genes_of_alleles.insert(subject_uniquename.clone(),
                                                     object_uniquename.clone());
                        continue;
                    }
                if feature_rel.rel_type.name == "part_of" &&
                    object_type_name == "genotype" {
                        let allele_and_expression =
                            AlleleAndExpression {
                                allele_uniquename: subject_uniquename.clone(),
                                expression: get_feat_rel_expression(feature_rel),
                            };
                        let entry = self.alleles_of_genotypes.entry(object_uniquename.clone());
                        entry.or_insert(Vec::new()).push(allele_and_expression);
                        continue;
                    }
            }
        }
    }

    fn make_location(&self, feat: &Feature) -> Option<ChromosomeLocation> {
        let feature_locs = feat.featurelocs.borrow();
        match feature_locs.get(0) {
            Some(feature_loc) => {
                let start_pos =
                    if feature_loc.fmin + 1 >= 1 {
                        (feature_loc.fmin + 1) as u32
                    } else {
                        panic!("start_pos less than 1");
                    };
                let end_pos =
                    if feature_loc.fmax >= 1 {
                        feature_loc.fmax as u32
                    } else {
                        panic!("start_end less than 1");
                    };
                Some(ChromosomeLocation {
                    chromosome_name: feature_loc.srcfeature.uniquename.clone(),
                    start_pos: start_pos,
                    end_pos: end_pos,
                    strand: match feature_loc.strand {
                        1 => Strand::Forward,
                        -1 => Strand::Reverse,
                        _ => panic!(),
                    },
                })
            },
            None => None,
        }
    }

    fn store_gene_details(&mut self, feat: &Feature) {
        let location = self.make_location(&feat);

        let organism = make_organism_short(&feat.organism);


        let feature_type =
            if let Some(transcript_type) =
                self.transcript_type_of_genes.get(&feat.uniquename) {
                    transcript_type.clone() + " " + &feat.feat_type.name
                } else {
                    feat.feat_type.name.clone()
                };

        let gene_feature = GeneDetails {
            uniquename: feat.uniquename.clone(),
            name: feat.name.clone(),
            organism: organism,
            product: None,
            name_descriptions: vec![],
            synonyms: vec![],
            feature_type: feature_type,
            characterisation_status: None,
            location: location,
            gene_neighbourhood: vec![],
            cds_location: None,
            cv_annotations: HashMap::new(),
            physical_interactions: vec![],
            genetic_interactions: vec![],
            ortholog_annotations: vec![],
            paralog_annotations: vec![],
            target_of_annotations: vec![],
            transcripts: vec![],
            genes_by_uniquename: HashMap::new(),
            genotypes_by_uniquename: HashMap::new(),
            alleles_by_uniquename: HashMap::new(),
            references_by_uniquename: HashMap::new(),
            terms_by_termid: HashMap::new(),
        };

        self.genes.insert(feat.uniquename.clone(), gene_feature);
    }

    fn store_genotype_details(&mut self, feat: &Feature) {
        let mut background = None;

        for prop in feat.featureprops.borrow().iter() {
            if prop.prop_type.name == "genotype_background" {
                background = prop.value.clone()
            }
        }

        self.genotypes.insert(feat.uniquename.clone(),
                              GenotypeDetails {
                                  uniquename: feat.uniquename.clone(),
                                  name: feat.name.clone(),
                                  background: background,
                                  expressed_alleles: vec![],
                                  cv_annotations: HashMap::new(),
                                  genes_by_uniquename: HashMap::new(),
                                  alleles_by_uniquename: HashMap::new(),
                                  references_by_uniquename: HashMap::new(),
                                  terms_by_termid: HashMap::new(),
                              });
    }

    fn store_allele_details(&mut self, feat: &Feature) {
        let mut allele_type = None;
        let mut description = None;

        for prop in feat.featureprops.borrow().iter() {
            match &prop.prop_type.name as &str {
                "allele_type" =>
                    allele_type = prop.value.clone(),
                "description" =>
                    description = prop.value.clone(),
                _ => ()
            }
        }

        if allele_type.is_none() {
            panic!("no allele_type cvtermprop for {}", &feat.uniquename);
        }

        let gene_uniquename =
            self.genes_of_alleles.get(&feat.uniquename).unwrap();
        let allele_details = AlleleShort {
            uniquename: feat.uniquename.clone(),
            name: feat.name.clone(),
            gene_uniquename: gene_uniquename.clone(),
            allele_type: allele_type.unwrap(),
            description: description,
        };
        self.alleles.insert(feat.uniquename.clone(), allele_details);
    }

    fn process_feature(&mut self, feat: &Feature) {
        match &feat.feat_type.name as &str {
            "gene" | "pseudogene" =>
                self.store_gene_details(feat),

            _ => {
                if TRANSCRIPT_FEATURE_TYPES.contains(&feat.feat_type.name.as_str()) {
                    self.transcripts.insert(feat.uniquename.clone(),
                                            TranscriptDetails {
                                                uniquename: feat.uniquename.clone(),
                                                name: feat.name.clone(),
                                            });
                }
            }
        }
    }

    fn process_features(&mut self) {
        for feat in &self.raw.features {
            if feat.feat_type.name != "genotype" && feat.feat_type.name != "allele" {
                self.process_feature(&feat);
            }
        }
    }

    fn add_interesting_parents(&mut self) {
        let mut interesting_parents_by_termid: HashMap<String, HashSet<String>> =
            HashMap::new();

        for cvtermpath in &self.raw.cvtermpaths {
            let subject_term = &cvtermpath.subject;
            let subject_termid = subject_term.termid();
            let object_term = &cvtermpath.object;
            let object_termid = object_term.termid();

            let rel_termid =
                match cvtermpath.rel_type {
                    Some(ref rel_type) => {
                        rel_type.termid()
                    },
                    None => panic!("no relation type for {} <-> {}\n",
                                   &subject_term.name, &object_term.name)
                };
            let rel_term_name =
                self.make_term_short(&rel_termid).name;

            if self.is_interesting_parent(&object_termid, &rel_term_name) {
                interesting_parents_by_termid
                    .entry(subject_termid.clone())
                    .or_insert(HashSet::new())
                    .insert(object_termid.into());
            };
        }

        for (termid, interesting_parents) in interesting_parents_by_termid {
            let mut term_details = self.terms.get_mut(&termid).unwrap();
            term_details.interesting_parents = interesting_parents;
        }
    }

    fn process_allele_features(&mut self) {
        for feat in &self.raw.features {
            if feat.feat_type.name == "allele" {
                self.store_allele_details(&feat);
            }
        }
    }

    fn process_genotype_features(&mut self) {
        for feat in &self.raw.features {
            if feat.feat_type.name == "genotype" {
                self.store_genotype_details(&feat);
            }
        }
    }

    fn add_gene_neighbourhoods(&mut self) {
        struct GeneAndLoc {
            gene_uniquename: String,
            loc: ChromosomeLocation,
        };

        let mut genes_and_locs: Vec<GeneAndLoc> = vec![];

        for gene_details in self.genes.values() {
            if let Some(ref location) = gene_details.location {
                genes_and_locs.push(GeneAndLoc {
                    gene_uniquename: gene_details.uniquename.clone(),
                    loc: location.clone(),
                });
            }
        }

        let cmp = |a: &GeneAndLoc, b: &GeneAndLoc| {
            let order = a.loc.chromosome_name.cmp(&b.loc.chromosome_name);
            if order == Ordering::Equal {
                a.loc.start_pos.cmp(&b.loc.start_pos)
            } else {
                order
            }
        };

        genes_and_locs.sort_by(cmp);

        for (i, this_gene_and_loc) in genes_and_locs.iter().enumerate() {
            let mut nearby_genes: Vec<GeneShort> = vec![];
            if i > 0 {
                let start_index =
                    if i > GENE_NEIGHBOURHOOD_DISTANCE {
                        i - GENE_NEIGHBOURHOOD_DISTANCE
                    } else {
                        0
                    };

                for back_index in (start_index..i).rev() {
                    let back_gene_and_loc = &genes_and_locs[back_index];

                    if back_gene_and_loc.loc.chromosome_name !=
                        this_gene_and_loc.loc.chromosome_name {
                            break;
                        }
                    let back_gene_short = self.make_gene_short(&back_gene_and_loc.gene_uniquename);
                    nearby_genes.insert(0, back_gene_short);
                }
            }

            let gene_short = self.make_gene_short(&this_gene_and_loc.gene_uniquename);
            nearby_genes.push(gene_short);

            if i < genes_and_locs.len() - 1 {
                let end_index =
                    if i + GENE_NEIGHBOURHOOD_DISTANCE >= genes_and_locs.len() {
                        genes_and_locs.len()
                    } else {
                        i + GENE_NEIGHBOURHOOD_DISTANCE + 1
                    };

                for forward_index in i+1..end_index {
                    let forward_gene_and_loc = &genes_and_locs[forward_index];

                    if forward_gene_and_loc.loc.chromosome_name !=
                        this_gene_and_loc.loc.chromosome_name {
                            break;
                        }

                    let forward_gene_short = self.make_gene_short(&forward_gene_and_loc.gene_uniquename);
                    nearby_genes.push(forward_gene_short);
                }
            }

            let mut this_gene_details =
                self.genes.get_mut(&this_gene_and_loc.gene_uniquename).unwrap();

            this_gene_details.gene_neighbourhood.append(&mut nearby_genes);
        }
    }

    fn add_alleles_to_genotypes(&mut self) {
        let mut alleles_to_add: HashMap<String, Vec<ExpressedAllele>> = HashMap::new();

        for genotype_uniquename in self.genotypes.keys() {
            let allele_uniquenames: Vec<AlleleAndExpression> =
                self.alleles_of_genotypes.get(genotype_uniquename).unwrap().clone();
            let expressed_allele_vec: Vec<ExpressedAllele> =
                allele_uniquenames.iter()
                .map(|ref allele_and_expression| {
                    ExpressedAllele {
                        allele_uniquename: allele_and_expression.allele_uniquename.clone(),
                        expression: allele_and_expression.expression.clone(),
                    }
                })
                .collect();

            alleles_to_add.insert(genotype_uniquename.clone(), expressed_allele_vec);
        }

        for (genotype_uniquename, genotype_details) in &mut self.genotypes {
            genotype_details.expressed_alleles =
                alleles_to_add.remove(genotype_uniquename).unwrap();
        }
    }

    // add interaction, ortholog and paralog annotations
    fn process_annotation_feature_rels(&mut self) {
        for feature_rel in self.raw.feature_relationships.iter() {
            let rel_name = &feature_rel.rel_type.name;
            let subject_uniquename = &feature_rel.subject.uniquename;
            let object_uniquename = &feature_rel.object.uniquename;

            for rel_config in FEATURE_REL_CONFIGS.iter() {
                if rel_name == rel_config.rel_type_name &&
                    is_gene_type(&feature_rel.subject.feat_type.name) &&
                    is_gene_type(&feature_rel.object.feat_type.name) {
                        let mut evidence: Option<Evidence> = None;

                        let borrowed_publications = feature_rel.publications.borrow();
                        let maybe_publication = borrowed_publications.get(0).clone();
                        let maybe_reference_uniquename =
                            match maybe_publication {
                                Some(publication) => Some(publication.uniquename.clone()),
                                None => None,
                            };


                        for prop in feature_rel.feature_relationshipprops.borrow().iter() {
                            if prop.prop_type.name == "evidence" {
                                if let Some(evidence_long) = prop.value.clone() {
                                    if let Some(code) = self.config.evidence_types.get(&evidence_long) {
                                        evidence = Some(code.clone());
                                    } else {
                                        evidence = Some(evidence_long);
                                    }
                                }
                            }
                        }

                        let evidence_clone = evidence.clone();

                        let gene_uniquename = subject_uniquename;
                        let gene_organism_short = {
                            self.genes.get(subject_uniquename).unwrap().organism.clone()
                        };
                        let other_gene_uniquename = object_uniquename;
                        let other_gene_organism_short = {
                             self.genes.get(object_uniquename).unwrap().organism.clone()
                        };
                        {
                            let mut gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                            match rel_config.annotation_type {
                                FeatureRelAnnotationType::Interaction => {
                                    let interaction_annotation =
                                        InteractionAnnotation {
                                            gene_uniquename: gene_uniquename.clone(),
                                            interactor_uniquename: other_gene_uniquename.clone(),
                                            evidence: evidence,
                                            reference_uniquename: maybe_reference_uniquename.clone(),
                                        };
                                    if rel_name == "interacts_physically" {
                                        gene_details.physical_interactions.push(interaction_annotation.clone());
                                    } else {
                                        if rel_name == "interacts_genetically" {
                                            gene_details.genetic_interactions.push(interaction_annotation.clone());
                                        } else {
                                            panic!("unknown interaction type: {}", rel_name);
                                        }
                                    };

                                    if let Some(ref_details) =
                                        if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                            self.references.get_mut(reference_uniquename)
                                        } else {
                                            None
                                        }
                                    {
                                        if rel_name == "interacts_physically" {
                                            ref_details.physical_interactions.push(interaction_annotation.clone());
                                        } else {
                                            if rel_name == "interacts_genetically" {
                                                ref_details.genetic_interactions.push(interaction_annotation.clone());
                                            } else {
                                                panic!("unknown interaction type: {}", rel_name);
                                            }
                                        };
                                    }
                                },
                                FeatureRelAnnotationType::Ortholog => {
                                    let ortholog_annotation =
                                        OrthologAnnotation {
                                            gene_uniquename: gene_uniquename.clone(),
                                            ortholog_uniquename: other_gene_uniquename.clone(),
                                            ortholog_organism: other_gene_organism_short,
                                            evidence: evidence,
                                            reference_uniquename: maybe_reference_uniquename.clone(),
                                        };
                                    gene_details.ortholog_annotations.push(ortholog_annotation.clone());
                                    if let Some(ref_details) =
                                        if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                            self.references.get_mut(reference_uniquename)
                                        } else {
                                            None
                                        }
                                    {
                                        ref_details.ortholog_annotations.push(ortholog_annotation);
                                    }
                                },
                                FeatureRelAnnotationType::Paralog => {
                                    let paralog_annotation =
                                        ParalogAnnotation {
                                            gene_uniquename: gene_uniquename.clone(),
                                            paralog_uniquename: other_gene_uniquename.clone(),
                                            evidence: evidence,
                                            reference_uniquename: maybe_reference_uniquename.clone(),
                                        };
                                    gene_details.paralog_annotations.push(paralog_annotation.clone());
                                    if let Some(ref_details) =
                                        if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                            self.references.get_mut(reference_uniquename)
                                        } else {
                                            None
                                        }
                                    {
                                        ref_details.paralog_annotations.push(paralog_annotation);
                                    }
                                }
                            }
                        }
                        {
                            // for orthologs and paralogs, store the reverses annotation too
                            let mut other_gene_details = self.genes.get_mut(object_uniquename).unwrap();
                            match rel_config.annotation_type {
                                FeatureRelAnnotationType::Interaction => {},
                                FeatureRelAnnotationType::Ortholog =>
                                    other_gene_details.ortholog_annotations.push(
                                        OrthologAnnotation {
                                            gene_uniquename: other_gene_uniquename.clone(),
                                            ortholog_uniquename: gene_uniquename.clone(),
                                            ortholog_organism: gene_organism_short,
                                            evidence: evidence_clone,
                                            reference_uniquename: maybe_reference_uniquename.clone(),
                                        }),
                                FeatureRelAnnotationType::Paralog =>
                                    other_gene_details.paralog_annotations.push(
                                        ParalogAnnotation {
                                            gene_uniquename: other_gene_uniquename.clone(),
                                            paralog_uniquename: gene_uniquename.clone(),
                                            evidence: evidence_clone,
                                            reference_uniquename: maybe_reference_uniquename
                                        }),
                            }
                        }
                    }
            }
        }

        for (_, ref_details) in &mut self.references {
            ref_details.physical_interactions.sort();
            ref_details.genetic_interactions.sort();
            ref_details.ortholog_annotations.sort();
            ref_details.paralog_annotations.sort();
        }

        for (_, gene_details) in &mut self.genes {
            gene_details.physical_interactions.sort();
            gene_details.genetic_interactions.sort();
            gene_details.ortholog_annotations.sort();
            gene_details.paralog_annotations.sort();
        }
    }

    fn matching_ext_config(&self, annotation_termid: &str,
                           rel_type_name: &str) -> Option<ExtensionConfig> {
        let ext_configs = &self.config.extensions;

        if let Some(annotation_term_details) = self.terms.get(annotation_termid) {
            for ext_config in ext_configs {
                if ext_config.rel_name == rel_type_name {
                    if let Some(if_descendent_of) = ext_config.if_descendent_of.clone() {
                        if annotation_term_details.interesting_parents.contains(&if_descendent_of) {
                            return Some((*ext_config).clone());
                        }
                    } else {
                        return Some((*ext_config).clone());
                    }
                }
            }
        } else {
            panic!("can't find details for term: {}\n", annotation_termid);
        }

        None
    }

    // create and returns any TargetOfAnnotations implied by the extension
    fn make_target_of_for_ext(&self, cv_name: &String,
                              gene_uniquename: &String,
                              reference_uniquename: &Option<String>,
                              annotation_termid: &String,
                              extension: &Vec<ExtPart>) -> Vec<(GeneUniquename, TargetOfAnnotation)> {
        let mut ret_vec = vec![];

        for ext_part in extension {
            let maybe_ext_config =
                self.matching_ext_config(annotation_termid, &ext_part.rel_type_name);
            if let ExtRange::Gene(ref target_gene_uniquename) = ext_part.ext_range {
                if let Some(ext_config) = maybe_ext_config {
                    if let Some(reciprocal_display_name) =
                        ext_config.reciprocal_display {
                            ret_vec.push(((*target_gene_uniquename).clone(),
                                          TargetOfAnnotation {
                                              ontology_name: cv_name.clone(),
                                              ext_rel_display_name: reciprocal_display_name,
                                              gene_uniquename: gene_uniquename.clone(),
                                              reference_uniquename: reference_uniquename.clone(),
                                          }));
                        }
                }
            }
        }

        ret_vec
    }

    fn add_target_of_annotations(&mut self) {
        let mut target_of_annotations: HashMap<GeneUniquename, HashSet<TargetOfAnnotation>> =
            HashMap::new();

        for (gene_uniquename, gene_details) in &self.genes {
            for (cv_name, feat_annotations) in &gene_details.cv_annotations {
                for feat_annotation in feat_annotations.iter() {
                    for detail in &feat_annotation.annotations {
                        let new_annotations =
                            self.make_target_of_for_ext(&cv_name, &gene_uniquename,
                                                        &detail.reference_uniquename,
                                                        &feat_annotation.term.termid, &detail.extension);
                        for (target_gene_uniquename, new_annotation) in new_annotations {
                            target_of_annotations
                                .entry(target_gene_uniquename.clone())
                                .or_insert(HashSet::new())
                                .insert(new_annotation);
                        }
                    }
                }
            }
        }

        for (gene_uniquename, mut target_of_annotations) in target_of_annotations {
            let mut gene_details = self.genes.get_mut(&gene_uniquename).unwrap();
            gene_details.target_of_annotations = target_of_annotations.drain().collect();
        }
    }

    fn process_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name != POMBASE_ANN_EXT_TERM_CV_NAME {
                let cv_config =
                    self.config.cv_config_by_name(&cvterm.cv.name);
                let annotation_feature_type =
                    cv_config.feature_type.clone();
                self.terms.insert(cvterm.termid(),
                                  TermDetails {
                                      name: cvterm.name.clone(),
                                      cv_name: cvterm.cv.name.clone(),
                                      annotation_feature_type: annotation_feature_type,
                                      interesting_parents: HashSet::new(),
                                      termid: cvterm.termid(),
                                      definition: cvterm.definition.clone(),
                                      is_obsolete: cvterm.is_obsolete,
                                      single_allele_genotype_uniquenames: HashSet::new(),
                                      rel_annotations: vec![],
                                      not_rel_annotations: vec![],
                                      genes_by_uniquename: HashMap::new(),
                                      genotypes_by_uniquename: HashMap::new(),
                                      alleles_by_uniquename: HashMap::new(),
                                      references_by_uniquename: HashMap::new(),
                                      terms_by_termid: HashMap::new(),
                                  });
            }
        }
    }

    fn get_ext_rel_display_name(&self, annotation_termid: &String,
                                ext_rel_name: &String) -> String {
        if let Some(ext_conf) = self.matching_ext_config(annotation_termid, ext_rel_name) {
            ext_conf.display_name.clone()
        } else {
            let re = Regex::new("_").unwrap();
            re.replace_all(&ext_rel_name, " ")
        }
    }

    fn process_extension_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                    if (*cvtermprop).prop_type.name.starts_with(ANNOTATION_EXT_REL_PREFIX) {
                        let ext_rel_name_str =
                            &(*cvtermprop).prop_type.name[ANNOTATION_EXT_REL_PREFIX.len()..];
                        let ext_rel_name = String::from(ext_rel_name_str);
                        let ext_range = (*cvtermprop).value.clone();
                        let range: ExtRange = if ext_range.starts_with("SP") {
                            ExtRange::Gene(ext_range)
                        } else {
                            ExtRange::Misc(ext_range)
                        };

                        if let Some(base_termid) =
                            self.base_term_of_extensions.get(&cvterm.termid()) {
                                let rel_type_display_name =
                                    self.get_ext_rel_display_name(&base_termid, &ext_rel_name);

                                self.parts_of_extensions.entry(cvterm.termid())
                                    .or_insert(Vec::new()).push(ExtPart {
                                        rel_type_name: String::from(ext_rel_name),
                                        rel_type_display_name: rel_type_display_name,
                                        ext_range: range,
                                    });
                            } else {
                                panic!("can't find details for term: {}\n", cvterm.termid());
                            }
                    }
                }
            }
        }
    }

    fn process_cvterm_rels(&mut self) {
        for cvterm_rel in &self.raw.cvterm_relationships {
            let subject_term = &cvterm_rel.subject;
            let object_term = &cvterm_rel.object;
            let rel_type = &cvterm_rel.rel_type;

            if subject_term.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                let subject_termid = subject_term.termid();
                if rel_type.name == "is_a" {
                    self.base_term_of_extensions.insert(subject_termid.clone(),
                                                        object_term.termid().clone());
                }
            }
        }

        for cvterm_rel in &self.raw.cvterm_relationships {
            let subject_term = &cvterm_rel.subject;
            let object_term = &cvterm_rel.object;
            let rel_type = &cvterm_rel.rel_type;

            if subject_term.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                let subject_termid = subject_term.termid();
                if rel_type.name != "is_a" {
                    if let Some(base_termid) =
                        self.base_term_of_extensions.get(&subject_term.termid()) {
                            let rel_type_display_name =
                                self.get_ext_rel_display_name(base_termid, &rel_type.name);

                            self.parts_of_extensions.entry(subject_termid)
                                .or_insert(Vec::new()).push(ExtPart {
                                    rel_type_name: rel_type.name.clone(),
                                    rel_type_display_name: rel_type_display_name,
                                    ext_range: ExtRange::Term(object_term.termid().clone()),
                                });
                        } else {
                            panic!("can't find details for {}\n", object_term.termid());
                        }
                }
            }
        }
    }

    fn process_feature_synonyms(&mut self) {
        for feature_synonym in self.raw.feature_synonyms.iter() {
            let feature = &feature_synonym.feature;
            let synonym = &feature_synonym.synonym;

            if let Some(ref mut gene_details) = self.genes.get_mut(&feature.uniquename) {
                gene_details.synonyms.push(SynonymDetails {
                    name: synonym.name.clone(),
                    synonym_type: synonym.synonym_type.name.clone()
                });
            }
        }
    }

    fn make_genotype_short(&self, genotype_uniquename: &str) -> GenotypeShort {
        let details = self.genotypes.get(genotype_uniquename).unwrap().clone();

        GenotypeShort {
            uniquename: details.uniquename,
            name: details.name,
            background: details.background,
            expressed_alleles: details.expressed_alleles,
        }
    }

    fn make_allele_short(&self, allele_uniquename: &str) -> AlleleShort {
        self.alleles.get(allele_uniquename).unwrap().clone()
    }

    // process feature properties stored as cvterms,
    // eg. characterisation_status and product
    fn process_props_from_feature_cvterms(&mut self) {
        for feature_cvterm in self.raw.feature_cvterms.iter() {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;

            let gene_uniquenames_vec: Vec<GeneUniquename> =
                if cvterm.cv.name == "PomBase gene products" {
                    if feature.feat_type.name == "polypeptide" {
                        if let Some(transcript_uniquename) =
                            self.transcripts_of_polypeptides.get(&feature.uniquename) {
                                if let Some(gene_uniquename) =
                                    self.genes_of_transcripts.get(transcript_uniquename) {
                                        vec![gene_uniquename.clone()]
                                    } else {
                                        vec![]
                                    }
                            } else {
                                vec![]
                            }
                    } else {
                        if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                            if let Some(gene_uniquename) =
                                self.genes_of_transcripts.get(&feature.uniquename) {
                                    vec![gene_uniquename.clone()]
                                } else {
                                    vec![]
                                }
                        } else {
                            vec![]
                        }
                    }
                } else {
                    vec![]
                };
            for gene_uniquename in &gene_uniquenames_vec {
                self.add_gene_product(&gene_uniquename, &cvterm.name);
            }
            if feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene" {
                if cvterm.cv.name == "PomBase gene characterisation status" {
                    self.add_characterisation_status(&feature.uniquename, &cvterm.name);
                } else {
                    if cvterm.cv.name == "name_description" {
                        self.add_name_description(&feature.uniquename, &cvterm.name);
                    }
                }
            }
        }
    }

    fn get_gene_prod_extension(&self, prod_value: &String) -> ExtPart {
        let ext_range =
            if prod_value.starts_with("PR:") {
                ExtRange::GeneProduct(prod_value.clone())
            } else {
                ExtRange::Misc(prod_value.clone())
            };

        ExtPart {
            rel_type_name: "active_form".into(),
            rel_type_display_name: "active form".into(),
            ext_range: ext_range,
        }
    }

    // return a fake extension for "with" properties on protein binding annotations
    fn get_with_extension(&self, with_value: &String) -> ExtPart {
        let ext_range =
            if with_value.starts_with("SP%") {
                ExtRange::Gene(with_value.clone())
            } else {
                if with_value.starts_with("PomBase:SP") {
                    let gene_uniquename =
                        String::from(&with_value[8..]);
                    ExtRange::Gene(gene_uniquename)
                } else {
                    if with_value.to_lowercase().starts_with("pfam:") {
                        ExtRange::Domain(with_value.clone())
                    } else {
                        ExtRange::Misc(with_value.clone())
                    }
                }
            };

        // a with property on a protein binding (GO:0005515) is
        // displayed as a binds extension
        // https://github.com/pombase/website/issues/108
        ExtPart {
            rel_type_name: "binds".into(),
            rel_type_display_name: "binds".into(),
            ext_range: ext_range,
        }
    }

    fn make_with_or_from_value(&self, with_or_from_value: String) -> WithFromValue {
        let db_prefix_patt = String::from("^") + DB_NAME + ":";
        let re = Regex::new(&db_prefix_patt).unwrap();
        let gene_uniquename = re.replace_all(&with_or_from_value, "");
        if self.genes.contains_key(&gene_uniquename) {
            let gene_short = self.make_gene_short(&gene_uniquename);
            WithFromValue::Gene(gene_short)
        } else {
            if self.terms.get(&with_or_from_value).is_some() {
                WithFromValue::Term(self.make_term_short(&with_or_from_value))
            } else {
                WithFromValue::Identifier(with_or_from_value)
            }
        }
    }

    // add the with value as a fake extension if the cvterm is_a protein binding,
    // otherwise return the value
    fn make_with_extension(&self, termid: &String, evidence_code: Option<String>,
                           extension: &mut Vec<ExtPart>,
                           with_value: String) -> WithFromValue {
        let base_termid =
            match self.base_term_of_extensions.get(termid) {
                Some(base_termid) => base_termid.clone(),
                None => termid.clone(),
            };

        let base_term_short = self.make_term_short(&base_termid);

        if evidence_code.is_some() &&
            evidence_code.unwrap() == "IPI" &&
            (base_term_short.termid == "GO:0005515" ||
            base_term_short.interesting_parents
            .contains("GO:0005515")) {
                extension.push(self.get_with_extension(&with_value));
            } else {
                return self.make_with_or_from_value(with_value);
            }
        WithFromValue::None
    }

    // process annotation
    fn process_feature_cvterms(&mut self) {
        for feature_cvterm in self.raw.feature_cvterms.iter() {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;

            let mut extension = vec![];

            if cvterm.cv.name == "PomBase gene characterisation status" ||
                cvterm.cv.name == "PomBase gene products" ||
                cvterm.cv.name == "name_description" {
                    continue;
                }

            let publication = &feature_cvterm.publication;
            let mut extra_props: HashMap<String, String> = HashMap::new();
            let mut conditions: Vec<TermId> = vec![];
            let mut with: WithFromValue = WithFromValue::None;
            let mut from: WithFromValue = WithFromValue::None;
            let mut qualifiers: Vec<Qualifier> = vec![];
            let mut evidence: Option<String> = None;
            let mut raw_with_value: Option<String> = None;
            for ref prop in feature_cvterm.feature_cvtermprops.borrow().iter() {
                match &prop.type_name() as &str {
                    "residue" | "scale" |
                    "quant_gene_ex_copies_per_cell" |
                    "quant_gene_ex_avg_copies_per_cell" => {
                        if let Some(value) = prop.value.clone() {
                            extra_props.insert(prop.type_name().clone(), value);
                        }
                    },
                    "evidence" =>
                        if let Some(evidence_long) = prop.value.clone() {
                            if let Some(code) = self.config.evidence_types.get(&evidence_long) {
                                evidence = Some(code.clone());
                            } else {
                                evidence = Some(evidence_long);
                            }
                        },
                    "condition" =>
                        if let Some(value) = prop.value.clone() {
                            conditions.push(value.clone());
                        },
                    "qualifier" =>
                        if let Some(value) = prop.value.clone() {
                            qualifiers.push(value);
                        },
                    "with" => {
                        raw_with_value = prop.value.clone();
                    },
                    "from" => {
                        if let Some(value) = prop.value.clone() {
                            from = self.make_with_or_from_value(value);
                        }
                    },
                    "gene_product_form_id" => {
                        if let Some(value) = prop.value.clone() {
                            extension.push(self.get_gene_prod_extension(&value));
                        }
                    },
                    _ => ()
                }
            }

            if let Some(value) = raw_with_value {
                let with_gene_short =
                    self.make_with_extension(&cvterm.termid(), evidence.clone(),
                                             &mut extension, value);
                if with_gene_short.is_some() {
                    with = with_gene_short;
                }
            }

            let mut maybe_genotype_uniquename = None;
            let mut gene_uniquenames_vec: Vec<GeneUniquename> =
                match &feature.feat_type.name as &str {
                    "polypeptide" => {
                        if let Some(transcript_uniquename) =
                            self.transcripts_of_polypeptides.get(&feature.uniquename) {
                                if let Some(gene_uniquename) =
                                    self.genes_of_transcripts.get(transcript_uniquename) {
                                        vec![gene_uniquename.clone()]
                                    } else {
                                        vec![]
                                    }
                            } else {
                                vec![]
                            }
                    },
                    "genotype" => {
                        let genotype_short = self.make_genotype_short(&feature.uniquename);
                        maybe_genotype_uniquename = Some(genotype_short.uniquename.clone());
                        genotype_short.expressed_alleles.iter()
                            .map(|expressed_allele| {
                                let allele_short =
                                    self.make_allele_short(&expressed_allele.allele_uniquename);
                                allele_short.gene_uniquename.clone()
                            })
                            .collect()
                    },
                    _ => {
                        if feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene" {
                            vec![feature.uniquename.clone()]
                        } else {
                            if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                                if let Some(gene_uniquename) =
                                    self.genes_of_transcripts.get(&feature.uniquename) {
                                        vec![gene_uniquename.clone()]
                                    } else {
                                        vec![]
                                    }
                            } else {
                                vec![]
                            }
                        }
                    }
                };

            gene_uniquenames_vec.dedup();

            let reference_uniquename =
                if publication.uniquename == "null" {
                    None
                } else {
                    Some(publication.uniquename.clone())
                };

            for gene_uniquename in &gene_uniquenames_vec {
                let mut extra_props_clone = extra_props.clone();
                let copies_per_cell = extra_props_clone.remove("quant_gene_ex_copies_per_cell");
                let avg_copies_per_cell = extra_props_clone.remove("quant_gene_ex_avg_copies_per_cell");
                let scale = extra_props_clone.remove("scale");
                let gene_ex_props =
                    if copies_per_cell.is_some() || avg_copies_per_cell.is_some() {
                        Some(GeneExProps {
                            copies_per_cell: copies_per_cell,
                            avg_copies_per_cell: avg_copies_per_cell,
                            scale: scale,
                        })
                    } else {
                        None
                    };
                let annotation = OntAnnotationDetail {
                    id: feature_cvterm.feature_cvterm_id,
                    gene_uniquename: gene_uniquename.clone(),
                    reference_uniquename: reference_uniquename.clone(),
                    genotype_uniquename: maybe_genotype_uniquename.clone(),
                    with: with.clone(),
                    from: from.clone(),
                    residue: extra_props_clone.remove("residue"),
                    gene_ex_props: gene_ex_props,
                    qualifiers: qualifiers.clone(),
                    evidence: evidence.clone(),
                    conditions: conditions.clone(),
                    extension: extension.clone(),
                };

                self.add_annotation(cvterm.borrow(), feature_cvterm.is_not,
                                    annotation);
            }
        }
    }

    fn make_term_annotations(&self, termid: &str, details: &Vec<Rc<OntAnnotationDetail>>,
                             is_not: bool)
                       -> Vec<(CvName, OntTermAnnotations)> {
        let term_short = self.make_term_short(termid);

        let cv_name = term_short.cv_name.clone();

        match cv_name.as_ref() {
            "gene_ex" => {
                if is_not {
                    panic!("gene_ex annotations can't be NOT annotations");
                }
                let mut qual_annotations =
                    OntTermAnnotations {
                        term: term_short.clone(),
                        is_not: false,
                        annotations: vec![],
                    };
                let mut quant_annotations =
                    OntTermAnnotations {
                        term: term_short.clone(),
                        is_not: false,
                        annotations: vec![],
                    };
                for detail in details {
                    if detail.gene_ex_props.is_some() {
                        quant_annotations.annotations.push(detail.clone())
                    } else {
                        qual_annotations.annotations.push(detail.clone())
                    }
                }

                let mut return_vec = vec![];

                if qual_annotations.annotations.len() > 0 {
                    return_vec.push((String::from("qualitative_gene_expression"),
                                     qual_annotations));
                }

                if quant_annotations.annotations.len() > 0 {
                    return_vec.push((String::from("quantitative_gene_expression"),
                                     quant_annotations));
                }

                return_vec
            },
            "fission_yeast_phenotype" => {
                let mut single_allele =
                    OntTermAnnotations {
                        term: term_short.clone(),
                        is_not: is_not,
                        annotations: vec![],
                    };
                let mut multi_allele =
                    OntTermAnnotations {
                        term: term_short.clone(),
                        is_not: is_not,
                        annotations: vec![],
                    };

                for detail in details {
                    let genotype_uniquename = detail.genotype_uniquename.clone().unwrap();
                    if let Some(genotype_details) = self.genotypes.get(&genotype_uniquename) {
                        if genotype_details.expressed_alleles.len() == 1 {
                            single_allele.annotations.push(detail.clone())
                        } else {
                            multi_allele.annotations.push(detail.clone())
                        }
                    } else {
                        panic!("can't find genotype details for {}\n", genotype_uniquename);
                    }
                }

                let mut return_vec = vec![];

                if single_allele.annotations.len() > 0 {
                    return_vec.push((String::from("single_allele_phenotype"),
                                     single_allele));
                }

                if multi_allele.annotations.len() > 0 {
                    return_vec.push((String::from("multi_allele_phenotype"),
                                     multi_allele));
                }

                return_vec
            },
            _ => {
                vec![(cv_name,
                      OntTermAnnotations {
                          term: term_short.clone(),
                          is_not: is_not,
                          annotations: details.clone(),
                      })]
            }
        }
    }

    // store the OntTermAnnotations in the TermDetails, GeneDetails,
    // GenotypeDetails and ReferenceDetails
    fn store_ont_annotations(&mut self, is_not: bool) {
        let ont_annotations = if is_not {
            &self.all_not_ont_annotations
        } else {
            &self.all_ont_annotations
        };

        let mut gene_annotation_by_term: HashMap<GeneUniquename, HashMap<TermId, Vec<Rc<OntAnnotationDetail>>>> =
            HashMap::new();
        let mut genotype_annotation_by_term: HashMap<GenotypeUniquename, HashMap<TermId, Vec<Rc<OntAnnotationDetail>>>> =
            HashMap::new();
        let mut ref_annotation_by_term: HashMap<String, HashMap<TermId, Vec<Rc<OntAnnotationDetail>>>> =
            HashMap::new();

        for (termid, annotations) in ont_annotations {
            let term_short = self.make_term_short(termid);

            if let Some(ref mut term_details) = self.terms.get_mut(termid) {
                let new_rel_ont_annotation = RelOntAnnotation {
                    rel_names: HashSet::new(),
                    term: term_short.clone(),
                    annotations: annotations.clone(),
                };
                if is_not {
                    term_details.not_rel_annotations.push(new_rel_ont_annotation);
                } else {
                    term_details.rel_annotations.push(new_rel_ont_annotation);
                }
            } else {
                panic!("missing termid: {}\n", termid);
            }

            for detail in annotations {
                gene_annotation_by_term.entry(detail.gene_uniquename.clone())
                    .or_insert(HashMap::new())
                    .entry(termid.clone())
                    .or_insert(vec![])
                    .push(detail.clone());

                if let Some(ref genotype_uniquename) = detail.genotype_uniquename {
                    let mut existing =
                        genotype_annotation_by_term.entry(genotype_uniquename.clone())
                        .or_insert(HashMap::new())
                        .entry(termid.clone())
                        .or_insert(vec![]);
                    if !existing.contains(detail) {
                        existing.push(detail.clone());
                    }
                }

                if let Some(reference_uniquename) = detail.reference_uniquename.clone() {
                    ref_annotation_by_term.entry(reference_uniquename)
                        .or_insert(HashMap::new())
                        .entry(termid.clone())
                        .or_insert(vec![])
                        .push(detail.clone());
                }

                for condition_termid in &detail.conditions {
                    let condition_term_short = {
                        self.make_term_short(&condition_termid)
                    };
                    if let Some(ref mut condition_term_details) =
                        self.terms.get_mut(&condition_termid.clone())
                    {
                        if condition_term_details.rel_annotations.len() == 0 {
                            condition_term_details.rel_annotations.push(
                                RelOntAnnotation {
                                    term: condition_term_short,
                                    rel_names: HashSet::new(),
                                    annotations: vec![],
                                });
                        }
                        if let Some(rel_annotation) = condition_term_details.rel_annotations.get_mut(0) {
                            rel_annotation.annotations.push(detail.clone())
                        }
                    }
                }
            }
        }

        for (gene_uniquename, term_annotation_map) in &gene_annotation_by_term {
            for (termid, details) in term_annotation_map {
                let new_annotations =
                    self.make_term_annotations(&termid, &details, is_not);

                let mut gene_details = self.genes.get_mut(gene_uniquename).unwrap();

                for (cv_name, new_annotation) in new_annotations {
                    gene_details.cv_annotations.entry(cv_name.clone())
                        .or_insert(Vec::new())
                        .push(new_annotation);
                }
            }

            let mut gene_details = self.genes.get_mut(gene_uniquename).unwrap();
            for (_, mut cv_annotations) in &mut gene_details.cv_annotations {
                cv_annotations.sort()
            }
        }

        for (genotype_uniquename, term_annotation_map) in &genotype_annotation_by_term {
            for (termid, details) in term_annotation_map {
                let new_annotations =
                    self.make_term_annotations(&termid, &details, is_not);

                let mut details = self.genotypes.get_mut(genotype_uniquename).unwrap();

                for (cv_name, new_annotation) in new_annotations {
                    details.cv_annotations.entry(cv_name.clone())
                        .or_insert(Vec::new())
                        .push(new_annotation);
                }
            }

            let mut details = self.genotypes.get_mut(genotype_uniquename).unwrap();
            for (_, mut cv_annotations) in &mut details.cv_annotations {
                cv_annotations.sort()
            }
        }

        for (reference_uniquename, ref_annotation_map) in &ref_annotation_by_term {
            for (termid, details) in ref_annotation_map {
                let new_annotations =
                    self.make_term_annotations(&termid, &details, is_not);

                let mut ref_details = self.references.get_mut(reference_uniquename).unwrap();

                for (cv_name, new_annotation) in new_annotations {
                    ref_details.cv_annotations.entry(cv_name).or_insert(Vec::new())
                        .push(new_annotation.clone());
                }
            }

            let mut ref_details = self.references.get_mut(reference_uniquename).unwrap();
            for (_, mut term_annotations) in &mut ref_details.cv_annotations {
                term_annotations.sort()
            }
        }
    }

    fn is_interesting_parent(&self, termid: &str, rel_name: &str) -> bool {
        for parent_conf in INTERESTING_PARENTS.iter() {
            if parent_conf.termid == termid &&
                parent_conf.rel_name == rel_name {
                    return true;
                }
        }

        for ext_conf in &self.config.extensions {
            if let Some(ref conf_termid) = ext_conf.if_descendent_of {
                if conf_termid == termid && rel_name == "is_a" {
                    return true;
                }
            }
        }

        false
    }

    fn process_cvtermpath(&mut self) {
        let mut annotation_by_id: HashMap<i32, Rc<OntAnnotationDetail>> = HashMap::new();
        let mut new_annotations: HashMap<TermId, HashMap<TermId, HashMap<i32, HashSet<RelName>>>> =
            HashMap::new();

        for cvtermpath in &self.raw.cvtermpaths {
            let subject_term = &cvtermpath.subject;
            let subject_termid = subject_term.termid();
            let object_term = &cvtermpath.object;
            let object_termid = object_term.termid();

            if let Some(subject_term_details) = self.terms.get(&subject_termid) {
                let rel_termid =
                    match cvtermpath.rel_type {
                        Some(ref rel_type) => {
                            rel_type.termid()
                        },
                        None => panic!("no relation type for {} <-> {}\n",
                                       &subject_term.name, &object_term.name)
                    };

                let rel_term_name =
                    self.make_term_short(&rel_termid).name;

                if rel_term_name == "has_part" &&
                    !HAS_PART_CV_NAMES.contains(&subject_term_details.cv_name.as_str()) {
                    continue;
                }

                if !DESCENDANT_REL_NAMES.contains(&rel_term_name.as_str()) {
                    continue;
                }

                let annotations = &subject_term_details.rel_annotations;
                for rel_annotation in annotations {
                    let RelOntAnnotation {
                        rel_names: _,
                        term: _,
                        annotations: existing_details
                    } = rel_annotation.clone();

                    for detail in &existing_details {
                        if !annotation_by_id.contains_key(&detail.id) {
                            annotation_by_id.insert(detail.id, detail.clone());
                        }
                        let (dest_termid, source_termid) =
                            (object_termid.clone(), subject_termid.clone());
                        new_annotations.entry(dest_termid)
                            .or_insert(HashMap::new())
                            .entry(source_termid)
                            .or_insert(HashMap::new())
                            .entry(detail.id)
                            .or_insert(HashSet::new())
                            .insert(rel_term_name.clone());
                    }
                }
            } else {
                panic!("TermDetails not found for {}", &subject_termid);
            }
        }

        for (dest_termid, dest_annotations_map) in new_annotations.drain() {
            for (source_termid, source_annotations_map) in dest_annotations_map {
                let mut new_details: Vec<Rc<OntAnnotationDetail>> = vec![];
                let mut all_rel_names: HashSet<String> = HashSet::new();
                for (id, rel_names) in source_annotations_map {
                    let detail = annotation_by_id.get(&id).unwrap().clone();
                    new_details.push(detail);
                    for rel_name in rel_names {
                        all_rel_names.insert(rel_name);
                    }
                }
                let source_term_short = self.make_term_short(&source_termid);
                let mut dest_term_details = {
                    self.terms.get_mut(&dest_termid).unwrap()
                };
                dest_term_details.rel_annotations.push(RelOntAnnotation {
                    rel_names: all_rel_names,
                    term: source_term_short.clone(),
                    annotations: new_details,
                });
            }
        }
    }

    fn make_metadata(&mut self) -> Metadata {
        let mut db_creation_datetime = None;

        for chadoprop in &self.raw.chadoprops {
            if chadoprop.prop_type.name == "db_creation_datetime" {
                db_creation_datetime = chadoprop.value.clone();
            }
        }

        const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
        const VERSION: &'static str = env!("CARGO_PKG_VERSION");

        Metadata {
            export_prog_name: String::from(PKG_NAME),
            export_prog_version: String::from(VERSION),
            db_creation_datetime: db_creation_datetime.unwrap(),
            gene_count: self.genes.len(),
            term_count: self.terms.len(),
        }
    }

    pub fn make_search_api_maps(&self) -> SearchAPIMaps {
        let mut gene_summaries: Vec<GeneSummary> = vec![];

        let gene_uniquenames: Vec<String> =
            self.genes.keys().map(|uniquename| uniquename.clone()).collect();

        for gene_uniquename in gene_uniquenames {
            gene_summaries.push(self.make_gene_summary(&gene_uniquename));
        }

        let mut term_summaries: HashSet<TermShort> = HashSet::new();

        let mut termid_genes: HashMap<TermId, HashSet<GeneUniquename>> = HashMap::new();
        let mut term_name_genes: HashMap<TermName, HashSet<GeneUniquename>> = HashMap::new();

        for (termid, term_details) in &self.terms {
            term_summaries.insert(self.make_term_short(&termid));
            for gene_uniquename in term_details.genes_by_uniquename.keys() {
                termid_genes.entry(termid.clone())
                    .or_insert(HashSet::new())
                    .insert(gene_uniquename.clone());
                term_name_genes.entry(term_details.name.clone())
                    .or_insert(HashSet::new())
                    .insert(gene_uniquename.clone());
            }
        }

        SearchAPIMaps {
            gene_summaries: gene_summaries,
            termid_genes: termid_genes,
            term_name_genes: term_name_genes,
            term_summaries: term_summaries,
        }
    }

    fn add_cv_annotations_to_maps(&self,
                                  identifier: &String,
                                  cv_annotations: &OntAnnotationMap,
                                  seen_references: &mut HashMap<String, ReferenceShortMap>,
                                  seen_genes: &mut HashMap<String, GeneShortMap>,
                                  seen_genotypes: &mut HashMap<String, GenotypeShortMap>,
                                  seen_alleles: &mut HashMap<String, AlleleShortMap>,
                                  seen_terms: &mut HashMap<String, TermShortMap>) {
        for (_, feat_annotations) in cv_annotations {
            for feat_annotation in feat_annotations.iter() {
                for detail in &feat_annotation.annotations {
                    self.add_ref_to_hash(seen_references,
                                         identifier.clone(), detail.reference_uniquename.clone());
                    for condition_termid in &detail.conditions {
                        self.add_term_to_hash(seen_terms,
                                              identifier.clone(), condition_termid.clone());
                    }
                    for ext_part in &detail.extension {
                        match ext_part.ext_range {
                            ExtRange::Term(ref range_termid) =>
                                self.add_term_to_hash(seen_terms, identifier.clone(), range_termid.clone()),
                            ExtRange::Gene(ref allele_gene_uniquename) =>
                                self.add_gene_to_hash(seen_genes, identifier.clone(),
                                                      allele_gene_uniquename.clone()),
                            _ => {},
                        }
                    }
                    if let Some(ref genotype_uniquename) = detail.genotype_uniquename {
                        let genotype = self.make_genotype_short(genotype_uniquename);
                        for expressed_allele in &genotype.expressed_alleles {
                            let allele_short =
                                self.add_allele_to_hash(seen_alleles, identifier.clone(),
                                                        expressed_allele.allele_uniquename.clone());
                            let allele_gene_uniquename =
                                allele_short.gene_uniquename.clone();
                            self.add_gene_to_hash(seen_genes, identifier.clone(), allele_gene_uniquename);
                        }
                        self.add_genotype_to_hash(seen_genotypes, identifier.clone(),
                                                  &genotype.uniquename);
                    }
                }
            }
        }
    }

    fn set_term_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_terms) = get_maps();

        for (termid, term_details) in &self.terms {
            for rel_annotation in &term_details.rel_annotations {
                for detail in &rel_annotation.annotations {
                    let gene_uniquename = detail.gene_uniquename.clone();
                    self.add_gene_to_hash(&mut seen_genes, termid.clone(), gene_uniquename.clone());
                    self.add_ref_to_hash(&mut seen_references, termid.clone(), detail.reference_uniquename.clone());
                    for condition_termid in &detail.conditions {
                        self.add_term_to_hash(&mut seen_terms, termid.clone(), condition_termid.clone());
                    }
                    for ext_part in &detail.extension {
                        match ext_part.ext_range {
                            ExtRange::Term(ref range_termid) =>
                                self.add_term_to_hash(&mut seen_terms, termid.clone(), range_termid.clone()),
                            ExtRange::Gene(ref allele_gene_uniquename) =>
                                self.add_gene_to_hash(&mut seen_genes, termid.clone(),
                                                 allele_gene_uniquename.clone()),
                            _ => {},
                        }
                    }
                    if let Some(ref genotype_uniquename) = detail.genotype_uniquename {
                        let genotype = self.make_genotype_short(genotype_uniquename);
                        for expressed_allele in &genotype.expressed_alleles {
                            let allele_short =
                                self.add_allele_to_hash(&mut seen_alleles, termid.clone(),
                                                        expressed_allele.allele_uniquename.clone());
                            let allele_gene_uniquename =
                                allele_short.gene_uniquename.clone();
                            self.add_gene_to_hash(&mut seen_genes, termid.clone(), allele_gene_uniquename);
                        }
                        self.add_genotype_to_hash(&mut seen_genotypes, termid.clone(),
                                                  &genotype.uniquename);
                    }
                }
            }
        }

        for (termid, term_details) in &mut self.terms {
            if let Some(genes) = seen_genes.remove(termid) {
                term_details.genes_by_uniquename = genes;
            }
            if let Some(genotypes) = seen_genotypes.remove(termid) {
                term_details.genotypes_by_uniquename = genotypes;
            }
            if let Some(alleles) = seen_alleles.remove(termid) {
                term_details.alleles_by_uniquename = alleles;
            }
            if let Some(references) = seen_references.remove(termid) {
                term_details.references_by_uniquename = references;
            }
             if let Some(terms) = seen_terms.remove(termid) {
                term_details.terms_by_termid = terms;
            }
       }
    }

    fn set_gene_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_terms) = get_maps();

        {
            for (gene_uniquename, gene_details) in &self.genes {
                self.add_cv_annotations_to_maps(&gene_uniquename,
                                                &gene_details.cv_annotations,
                                                &mut seen_references,
                                                &mut seen_genes,
                                                &mut seen_genotypes,
                                                &mut seen_alleles,
                                                &mut seen_terms);

                let interaction_iter =
                    gene_details.physical_interactions.iter().chain(&gene_details.genetic_interactions);
                for interaction in interaction_iter {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename.clone(), interaction.reference_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), interaction.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), interaction.interactor_uniquename.clone());
                }

                for ortholog_annotation in &gene_details.ortholog_annotations {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename.clone(), ortholog_annotation.reference_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), ortholog_annotation.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), ortholog_annotation.ortholog_uniquename.clone());
                }
                for paralog_annotation in &gene_details.paralog_annotations {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename.clone(), paralog_annotation.reference_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), paralog_annotation.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), paralog_annotation.paralog_uniquename.clone());
                }
                for target_of_annotation in &gene_details.target_of_annotations {
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(),
                                     target_of_annotation.gene_uniquename.clone());
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename.clone(),
                                    target_of_annotation.reference_uniquename.clone());
                }
            }
        }

        for (gene_uniquename, gene_details) in &mut self.genes {
            if let Some(references) = seen_references.remove(gene_uniquename) {
                gene_details.references_by_uniquename = references;
            }
            if let Some(alleles) = seen_alleles.remove(gene_uniquename) {
                gene_details.alleles_by_uniquename = alleles;
            }
            if let Some(genes) = seen_genes.remove(gene_uniquename) {
                gene_details.genes_by_uniquename = genes;
            }
            if let Some(genotypes) = seen_genotypes.remove(gene_uniquename) {
                gene_details.genotypes_by_uniquename = genotypes;
            }
            if let Some(terms) = seen_terms.remove(gene_uniquename) {
                gene_details.terms_by_termid = terms;
            }
        }
    }

    fn set_genotype_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_terms) = get_maps();

        for (genotype_uniquename, genotype_details) in &self.genotypes {
            self.add_cv_annotations_to_maps(&genotype_uniquename,
                                            &genotype_details.cv_annotations,
                                            &mut seen_references,
                                            &mut seen_genes,
                                            &mut seen_genotypes,
                                            &mut seen_alleles,
                                            &mut seen_terms);
        }

        for (genotype_uniquename, genotype_details) in &mut self.genotypes {
            if let Some(references) = seen_references.remove(genotype_uniquename) {
                genotype_details.references_by_uniquename = references;
            }
            if let Some(alleles) = seen_alleles.remove(genotype_uniquename) {
                genotype_details.alleles_by_uniquename = alleles;
            }
            if let Some(genotypes) = seen_genes.remove(genotype_uniquename) {
                genotype_details.genes_by_uniquename = genotypes;
            }
            if let Some(terms) = seen_terms.remove(genotype_uniquename) {
                genotype_details.terms_by_termid = terms;
            }
        }
    }

    fn set_reference_details_maps(&mut self) {
        type GeneShortMap = HashMap<GeneUniquename, GeneShort>;
        let mut seen_genes: HashMap<String, GeneShortMap> = HashMap::new();

        type GenotypeShortMap = HashMap<GenotypeUniquename, GenotypeShort>;
        let mut seen_genotypes: HashMap<ReferenceUniquename, GenotypeShortMap> = HashMap::new();

        type AlleleShortMap = HashMap<AlleleUniquename, AlleleShort>;
        let mut seen_alleles: HashMap<TermId, AlleleShortMap> = HashMap::new();

        type TermShortMap = HashMap<TermId, TermShort>;
        let mut seen_terms: HashMap<GeneUniquename, TermShortMap> = HashMap::new();

        {
            let mut add_gene_to_hash =
                |reference_uniquename: ReferenceUniquename, gene_uniquename: GeneUniquename| {
                    seen_genes
                        .entry(reference_uniquename.clone())
                        .or_insert(HashMap::new())
                        .insert(gene_uniquename.clone(),
                                self.make_gene_short(&gene_uniquename));
                };

            let mut add_allele_to_hash =
                |reference_uniquename: ReferenceUniquename, allele_uniquename: AlleleUniquename| {
                    let allele_short = self.make_allele_short(&allele_uniquename);
                    seen_alleles
                        .entry(reference_uniquename.clone())
                        .or_insert(HashMap::new())
                        .insert(allele_uniquename.clone(), allele_short.clone());
                    allele_short
                };

            let mut add_term_to_hash =
                |reference_uniquename: ReferenceUniquename, other_termid: TermId| {
                    seen_terms
                        .entry(reference_uniquename.clone())
                        .or_insert(HashMap::new())
                        .insert(other_termid.clone(),
                                self.make_term_short(&other_termid));
                };

            for (reference_uniquename, reference_details) in &self.references {
                for (_, feat_annotations) in &reference_details.cv_annotations {
                    for feat_annotation in feat_annotations.iter() {
                        for detail in &feat_annotation.annotations {
                            add_gene_to_hash(reference_uniquename.clone(),
                                             detail.gene_uniquename.clone());
                            for condition_termid in &detail.conditions {
                                add_term_to_hash(reference_uniquename.clone(), condition_termid.clone());
                            }
                            for ext_part in &detail.extension {
                                match ext_part.ext_range {
                                    ExtRange::Term(ref range_termid) =>
                                        add_term_to_hash(reference_uniquename.clone(), range_termid.clone()),
                                    ExtRange::Gene(ref allele_gene_uniquename) =>
                                        add_gene_to_hash(reference_uniquename.clone(),
                                                         allele_gene_uniquename.clone()),
                                    _ => {},
                                }
                            }
                            if let Some(ref genotype_uniquename) = detail.genotype_uniquename {
                                let genotype = self.make_genotype_short(genotype_uniquename);
                                for expressed_allele in &genotype.expressed_alleles {
                                    let allele_short =
                                        add_allele_to_hash(reference_uniquename.clone(),
                                                           expressed_allele.allele_uniquename.clone());
                                    let allele_gene_uniquename =
                                        allele_short.gene_uniquename.clone();
                                    add_gene_to_hash(reference_uniquename.clone(), allele_gene_uniquename);
                                }
                                self.add_genotype_to_hash(&mut seen_genotypes, reference_uniquename.clone(),
                                                          &genotype.uniquename);
                           }
                        }
                   }
                }

                let interaction_iter =
                    reference_details.physical_interactions.iter().chain(&reference_details.genetic_interactions);
                for interaction in interaction_iter {
                    add_gene_to_hash(reference_uniquename.clone(), interaction.gene_uniquename.clone());
                    add_gene_to_hash(reference_uniquename.clone(), interaction.interactor_uniquename.clone());
                }

                for ortholog_annotation in &reference_details.ortholog_annotations {
                    add_gene_to_hash(reference_uniquename.clone(), ortholog_annotation.gene_uniquename.clone());
                    add_gene_to_hash(reference_uniquename.clone(), ortholog_annotation.ortholog_uniquename.clone());
                }
                for paralog_annotation in &reference_details.paralog_annotations {
                    add_gene_to_hash(reference_uniquename.clone(), paralog_annotation.gene_uniquename.clone());
                    add_gene_to_hash(reference_uniquename.clone(), paralog_annotation.paralog_uniquename.clone());
                }

            }
        }

        for (reference_uniquename, reference_details) in &mut self.references {
            if let Some(genes) = seen_genes.remove(reference_uniquename) {
                reference_details.genes_by_uniquename = genes;
            }
            if let Some(genotypes) = seen_genotypes.remove(reference_uniquename) {
                reference_details.genotypes_by_uniquename = genotypes;
            }
            if let Some(alleles) = seen_alleles.remove(reference_uniquename) {
                reference_details.alleles_by_uniquename = alleles;
            }
            if let Some(terms) = seen_terms.remove(reference_uniquename) {
                reference_details.terms_by_termid = terms;
            }
        }
    }

    pub fn set_counts(&mut self) {
        let mut term_seen_genes: HashMap<TermId, HashSet<GeneUniquename>> = HashMap::new();
        let mut term_seen_genotypes: HashMap<TermId, HashSet<GenotypeUniquename>> = HashMap::new();
        let mut term_seen_single_allele_genotypes: HashMap<TermId, HashSet<GenotypeUniquename>> = HashMap::new();
        let mut ref_seen_genes: HashMap<ReferenceUniquename, HashSet<GeneUniquename>> = HashMap::new();

        for (termid, term_details) in &self.terms {
            let mut seen_genes: HashSet<GeneUniquename> = HashSet::new();
            let mut seen_genotypes: HashSet<GenotypeUniquename> = HashSet::new();
            let mut seen_single_allele_genotypes: HashSet<GenotypeUniquename> = HashSet::new();
            for rel_annotation in &term_details.rel_annotations {
                for annotation in &rel_annotation.annotations {
                    seen_genes.insert(annotation.gene_uniquename.clone());
                    if let Some(ref genotype_uniquename) = annotation.genotype_uniquename {
                        seen_genotypes.insert(genotype_uniquename.clone());
                        let genotype = self.genotypes.get(genotype_uniquename).unwrap();
                        if genotype.expressed_alleles.len() == 1 {
                            seen_single_allele_genotypes.insert(genotype_uniquename.clone());
                        }
                    }
                }
            }
            term_seen_genes.insert(termid.clone(), seen_genes);
            term_seen_genotypes.insert(termid.clone(), seen_genotypes);
            term_seen_single_allele_genotypes.insert(termid.clone(), seen_single_allele_genotypes);
        }

        for (reference_uniquename, reference_details) in &self.references {
            let mut seen_genes: HashSet<GeneUniquename> = HashSet::new();
            for (_, rel_annotations) in &reference_details.cv_annotations {
                for rel_annotation in rel_annotations {
                    for annotation in &rel_annotation.annotations {
                        if !rel_annotation.is_not {
                            seen_genes.insert(annotation.gene_uniquename.clone());
                        }
                    }
                }
            }
            let interaction_iter =
                reference_details.physical_interactions.iter().chain(&reference_details.genetic_interactions);
            for interaction in interaction_iter {
                seen_genes.insert(interaction.gene_uniquename.clone());
                seen_genes.insert(interaction.interactor_uniquename.clone());
            }
            for ortholog_annotation in &reference_details.ortholog_annotations {
                seen_genes.insert(ortholog_annotation.gene_uniquename.clone());
            }
            ref_seen_genes.insert(reference_uniquename.clone(), seen_genes);
        }

        for (_, gene_details) in &mut self.genes {
            for (_, feat_annotations) in &mut gene_details.cv_annotations {
                for mut feat_annotation in feat_annotations.iter_mut() {
                    feat_annotation.term.gene_count =
                        term_seen_genes.get(&feat_annotation.term.termid).unwrap().len();
                    feat_annotation.term.genotype_count =
                        term_seen_genotypes.get(&feat_annotation.term.termid).unwrap().len();
                }
            }

            for (reference_uniquename, reference_short) in
                &mut gene_details.references_by_uniquename {
                    reference_short.gene_count =
                        ref_seen_genes.get(reference_uniquename).unwrap().len();
                }
        }

        for (_, genotype_details) in &mut self.genotypes {
            for (_, feat_annotations) in &mut genotype_details.cv_annotations {
                for mut feat_annotation in feat_annotations.iter_mut() {
                    feat_annotation.term.genotype_count =
                        term_seen_genotypes.get(&feat_annotation.term.termid).unwrap().len();
                }
            }
        }

        for (_, ref_details) in &mut self.references {
            for (_, ref_annotations) in &mut ref_details.cv_annotations {
                for ref_annotation in ref_annotations {
                    ref_annotation.term.gene_count =
                        term_seen_genes.get(&ref_annotation.term.termid).unwrap().len();
                    ref_annotation.term.genotype_count =
                        term_seen_genotypes.get(&ref_annotation.term.termid).unwrap().len();
                }
            }
        }

        for (_, term_details) in &mut self.terms {
            for rel_annotation in &mut term_details.rel_annotations {
                rel_annotation.term.gene_count =
                    term_seen_genes.get(&rel_annotation.term.termid).unwrap().len();
                rel_annotation.term.genotype_count =
                    term_seen_genotypes.get(&rel_annotation.term.termid).unwrap().len();
            }

            for (reference_uniquename, reference_short) in
                &mut term_details.references_by_uniquename {
                    reference_short.gene_count =
                        ref_seen_genes.get(reference_uniquename).unwrap().len();
                }

            term_details.single_allele_genotype_uniquenames =
                term_seen_single_allele_genotypes.remove(&term_details.termid).unwrap();
        }
    }

    pub fn get_web_data(&mut self) -> WebData {
        self.process_references();
        self.make_feature_rel_maps();
        self.process_features();
        self.add_gene_neighbourhoods();
        self.process_props_from_feature_cvterms();
        self.process_allele_features();
        self.process_genotype_features();
        self.add_alleles_to_genotypes();
        self.process_cvterms();
        self.add_interesting_parents();
        self.process_cvterm_rels();
        self.process_extension_cvterms();
        self.process_feature_synonyms();
        self.process_feature_cvterms();
        self.store_ont_annotations(false);
        self.store_ont_annotations(true);
        self.process_cvtermpath();
        self.process_annotation_feature_rels();
        self.add_target_of_annotations();
        self.set_term_details_maps();
        self.set_gene_details_maps();
        self.set_genotype_details_maps();
        self.set_reference_details_maps();
        self.set_counts();

        let mut web_data_terms: IdTermDetailsMap = HashMap::new();

        let search_api_maps = self.make_search_api_maps();

        for (termid, term_details) in self.terms.drain() {
            web_data_terms.insert(termid.clone(), Rc::new(term_details));
        }

        self.terms = HashMap::new();

        let mut used_terms: IdTermDetailsMap = HashMap::new();

        // remove terms with no annotation
        for (termid, term_details) in &web_data_terms {
            if term_details.rel_annotations.len() > 0 {
                used_terms.insert(termid.clone(), term_details.clone());
            }
        }

        let metadata = self.make_metadata();

        WebData {
            genes: self.genes.clone(),
            genotypes: self.genotypes.clone(),
            terms: web_data_terms,
            used_terms: used_terms,
            metadata: metadata,
            references: self.references.clone(),

            search_api_maps: search_api_maps,
        }
    }
}
