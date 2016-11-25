use std::rc::Rc;
use std::collections::hash_map::HashMap;
use std::collections::HashSet;
use std::borrow::Borrow;

use regex::Regex;

use db::*;
use web::data::*;
use web::config::*;

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

    genes: UniquenameGeneMap,
    transcripts: UniquenameTranscriptMap,
    genotypes: UniquenameGenotypeMap,
    alleles: UniquenameAlleleMap,
    terms: HashMap<TermId, TermDetails>,
    references: IdReferenceMap,
    all_ont_annotations: Vec<OntAnnotation>,

    genes_of_transcripts: HashMap<String, String>,
    transcripts_of_polypeptides: HashMap<String, String>,
    genes_of_alleles: HashMap<String, String>,
    alleles_of_genotypes: HashMap<String, Vec<AlleleAndExpression>>,

    // a map from IDs of terms from the "PomBase annotation extension terms" cv
    // to a Vec of the details of each of the extension
    parts_of_extensions: HashMap<String, Vec<ExtPart>>,

    base_term_of_extensions: HashMap<TermId, TermId>,
}

fn get_feat_rel_expression(feature_relationship: &FeatureRelationship) -> Option<String> {
    for prop in feature_relationship.feature_relationshipprops.borrow().iter() {
        if prop.prop_type.name == "expression" {
            return prop.value.clone();
        }
    }

    None
}

impl <'a> WebDataBuild<'a> {
    pub fn new(raw: &'a Raw) -> WebDataBuild<'a> {
        WebDataBuild {
            raw: raw,

            genes: HashMap::new(),
            transcripts: HashMap::new(),
            genotypes: HashMap::new(),
            alleles: HashMap::new(),
            terms: HashMap::new(),
            references: HashMap::new(),
            all_ont_annotations: vec![],

            genes_of_transcripts: HashMap::new(),
            transcripts_of_polypeptides: HashMap::new(),
            genes_of_alleles: HashMap::new(),
            alleles_of_genotypes: HashMap::new(),

            parts_of_extensions: HashMap::new(),

            base_term_of_extensions: HashMap::new(),
        }
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
            synonyms: gene_details.synonyms.clone(),
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
                    authors_abbrev: reference_details.authors_abbrev.clone(),
                };

            Some(reference_short)
        }
    }

    fn make_term_short(&self, termid: &str) -> TermShort {
        if let Some(term_details) = self.terms.get(termid) {
            TermShort {
                name: term_details.name.clone(),
                cv_name: term_details.cv_name.clone(),
                termid: term_details.termid.clone(),
                is_obsolete: term_details.is_obsolete,
                gene_count: None,
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

    fn add_annotation(&mut self, cvterm: &Cvterm,
                      annotation_template: OntAnnotation) {
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

        let term_short = self.make_term_short(&termid).clone();
        let ont_annotation =
            OntAnnotation {
                term: Some(term_short),
                extension: extension_parts,
                .. annotation_template
            };

        self.all_ont_annotations.push(ont_annotation);
   }

    fn process_references(&mut self) {
        for rc_publication in &self.raw.publications {
            let reference_uniquename = &rc_publication.uniquename;

            let mut pubmed_authors: Option<String> = None;
            let mut pubmed_publication_date: Option<String> = None;

            for prop in rc_publication.publicationprops.borrow().iter() {
                match &prop.prop_type.name as &str {
                    "pubmed_publication_date" =>
                        pubmed_publication_date = Some(prop.value.clone()),
                    "pubmed_authors" =>
                        pubmed_authors = Some(prop.value.clone()),
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
                                       authors: pubmed_authors.clone(),
                                       authors_abbrev: authors_abbrev,
                                       pubmed_publication_date: pubmed_publication_date.clone(),
                                       publication_year: publication_year,
                                       annotations: HashMap::new(),
                                       interaction_annotations: HashMap::new(),
                                       ortholog_annotations: vec![],
                                       paralog_annotations: vec![],
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

            if subject_type_name == "mRNA" &&
                rel_name == "part_of" &&
                (object_type_name == "gene" || object_type_name == "pseudogene") {
                    self.genes_of_transcripts.insert(subject_uniquename.clone(),
                                                     object_uniquename.clone());
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

        let gene_feature = GeneDetails {
            uniquename: feat.uniquename.clone(),
            name: feat.name.clone(),
            organism: organism,
            product: None,
            synonyms: vec![],
            feature_type: feat.feat_type.name.clone(),
            characterisation_status: None,
            location: location,
            cds_location: None,
            annotations: HashMap::new(),
            interaction_annotations: HashMap::new(),
            ortholog_annotations: vec![],
            paralog_annotations: vec![],
            transcripts: vec![],
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
                              GenotypeShort {
                                  uniquename: feat.uniquename.clone(),
                                  name: feat.name.clone(),
                                  background: background,
                                  expressed_alleles: vec![],
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
            gene: self.make_gene_short(&gene_uniquename),
            allele_type: allele_type.unwrap(),
            description: description,
        };
        self.alleles.insert(feat.uniquename.clone(), allele_details);
    }

    fn process_feature(&mut self, feat: &Feature) {
        match &feat.feat_type.name as &str {
            "gene" | "pseudogene" =>
                self.store_gene_details(feat),

            // TODO: mRNA isn't the only transcript type
            "mRNA" => {
                self.transcripts.insert(feat.uniquename.clone(),
                                        TranscriptDetails {
                                            uniquename: feat.uniquename.clone(),
                                            name: feat.name.clone(),
                                        });
            },
            _ => (),
        }
    }

    fn process_features(&mut self) {
        for feat in &self.raw.features {
            if feat.feat_type.name != "genotype" && feat.feat_type.name != "allele" {
                self.process_feature(&feat);
            }
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

    fn add_alleles_to_genotypes(&mut self) {
        let mut alleles_to_add: HashMap<String, Vec<ExpressedAllele>> = HashMap::new();

        for genotype_uniquename in self.genotypes.keys() {
            let allele_uniquenames: Vec<AlleleAndExpression> =
                self.alleles_of_genotypes.get(genotype_uniquename).unwrap().clone();
            let expressed_allele_vec: Vec<ExpressedAllele> =
                allele_uniquenames.iter()
                .map(|ref allele_and_expression| {
                    ExpressedAllele {
                        allele: self.make_allele_short(&allele_and_expression.allele_uniquename),
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
                    (feature_rel.subject.feat_type.name == "gene" ||
                     feature_rel.subject.feat_type.name == "pseudogene") &&
                    (feature_rel.object.feat_type.name == "gene" ||
                     feature_rel.object.feat_type.name == "pseudogene") {
                        let mut evidence: Option<Evidence> = None;

                        let borrowed_publications = feature_rel.publications.borrow();
                        let maybe_publication = borrowed_publications.get(0).clone();
                        let maybe_reference_short = match maybe_publication {
                            Some(publication) => self.make_reference_short(&publication.uniquename),
                            None => None,
                        };

                        for prop in feature_rel.feature_relationshipprops.borrow().iter() {
                            if prop.prop_type.name == "evidence" {
                                evidence = prop.value.clone();
                            }
                        }

                        let evidence_clone = evidence.clone();

                        let (gene, gene_organism_short) = {
                            (self.make_gene_short(&subject_uniquename).clone(),
                             self.genes.get(subject_uniquename).unwrap().organism.clone())
                        };
                        let gene_clone = gene.clone();
                        let (other_gene, other_gene_organism_short) = {
                            (self.make_gene_short(&object_uniquename).clone(),
                             self.genes.get(object_uniquename).unwrap().organism.clone())
                        };
                        {
                            let mut gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                            match rel_config.annotation_type {
                                FeatureRelAnnotationType::Interaction => {
                                    let interaction_annotation =
                                        InteractionAnnotation {
                                            gene: gene.clone(),
                                            interactor: other_gene.clone(),
                                            evidence: evidence,
                                            reference: maybe_reference_short.clone(),
                                        };
                                    gene_details.interaction_annotations
                                        .entry(rel_name.clone()).or_insert(Vec::new()).push(interaction_annotation.clone());
                                    if let Some(ref_details) =
                                        if let Some(ref reference_short) = maybe_reference_short {
                                            self.references.get_mut(&reference_short.uniquename)
                                        } else {
                                            None
                                        }
                                    {
                                        ref_details.interaction_annotations
                                            .entry(rel_name.clone()).or_insert(Vec::new()).push(interaction_annotation);
                                    }
                                },
                                FeatureRelAnnotationType::Ortholog => {
                                    let ortholog_annotation =
                                        OrthologAnnotation {
                                            gene: gene.clone(),
                                            ortholog: other_gene.clone(),
                                            ortholog_organism: other_gene_organism_short,
                                            evidence: evidence,
                                            reference: maybe_reference_short.clone(),
                                        };
                                    gene_details.ortholog_annotations.push(ortholog_annotation.clone());
                                    if let Some(ref_details) =
                                        if let Some(ref reference_short) = maybe_reference_short {
                                            self.references.get_mut(&reference_short.uniquename)
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
                                            gene: gene.clone(),
                                            paralog: other_gene.clone(),
                                            evidence: evidence,
                                            reference: maybe_reference_short.clone(),
                                        };
                                    gene_details.paralog_annotations.push(paralog_annotation.clone());
                                    if let Some(ref_details) =
                                        if let Some(ref reference_short) = maybe_reference_short {
                                            self.references.get_mut(&reference_short.uniquename)
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
                            let mut other_gene_details = self.genes.get_mut(object_uniquename).unwrap();
                            match rel_config.annotation_type {
                                FeatureRelAnnotationType::Interaction => {},
                                FeatureRelAnnotationType::Ortholog =>
                                    other_gene_details.ortholog_annotations.push(
                                        OrthologAnnotation {
                                            gene: other_gene,
                                            ortholog: gene_clone,
                                            ortholog_organism: gene_organism_short,
                                            evidence: evidence_clone,
                                            reference: maybe_reference_short.clone(),
                                        }),
                                FeatureRelAnnotationType::Paralog =>
                                    other_gene_details.paralog_annotations.push(
                                        ParalogAnnotation {
                                            gene: other_gene,
                                            paralog: gene_clone,
                                            evidence: evidence_clone,
                                            reference: maybe_reference_short,
                                        }),
                            }
                        }
                    }
            }
        }
    }

    fn process_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name != POMBASE_ANN_EXT_TERM_CV_NAME {
                self.terms.insert(cvterm.termid(),
                                  TermDetails {
                                      name: cvterm.name.clone(),
                                      cv_name: cvterm.cv.name.clone(),
                                      termid: cvterm.termid(),
                                      definition: cvterm.definition.clone(),
                                      is_obsolete: cvterm.is_obsolete,
                                      genes: vec![],
                                      annotations: vec![],
                                  });
            }
        }
    }

    fn process_extension_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                    if (*cvtermprop).prop_type.name.starts_with(ANNOTATION_EXT_REL_PREFIX) {
                        let ext_rel_name: &str = &(*cvtermprop).prop_type.name[ANNOTATION_EXT_REL_PREFIX.len()..];
                        let ext_range = (*cvtermprop).value.clone();
                        let range: ExtRange = if ext_range.starts_with("SP") {
                            ExtRange::Gene(self.make_gene_short(&ext_range))
                        } else {
                            ExtRange::Misc(ext_range)
                        };

                        self.parts_of_extensions.entry(cvterm.termid())
                            .or_insert(Vec::new()).push(ExtPart {
                                rel_type_name: String::from(ext_rel_name),
                                ext_range: range,
                            });
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
                } else {
                    let object_term_short = self.make_term_short(&object_term.termid());

                    self.parts_of_extensions.entry(subject_termid)
                        .or_insert(Vec::new()).push(ExtPart {
                            rel_type_name: rel_type.name.clone(),
                            ext_range: ExtRange::Term(object_term_short.clone()),
                        });
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
        self.genotypes.get(genotype_uniquename).unwrap().clone()
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
                if feature.feat_type.name == "polypeptide" && cvterm.cv.name == "PomBase gene products" {
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
                    vec![]
                };
            for gene_uniquename in &gene_uniquenames_vec {
                self.add_gene_product(&gene_uniquename, &cvterm.name);
            }
            if(feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene")
                && cvterm.cv.name == "PomBase gene characterisation status" {
                    self.add_characterisation_status(&feature.uniquename, &cvterm.name);
                }
        }
    }

    // process annotation
    fn process_feature_cvterms(&mut self) {
        let db_prefix_patt = String::from("^") + DB_NAME + ":";
        for feature_cvterm in self.raw.feature_cvterms.iter() {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;
            let publication = &feature_cvterm.publication;
            let mut extra_props: HashMap<String, String> = HashMap::new();
            let mut conditions: Vec<TermShort> = vec![];
            let mut qualifiers: Vec<Qualifier> = vec![];
            for ref prop in feature_cvterm.feature_cvtermprops.borrow().iter() {
                match &prop.type_name() as &str {
                    "evidence" | "residue" | "scale" |
                    "quant_gene_ex_copies_per_cell" |
                    "quant_gene_ex_avg_copies_per_cell" => {
                        if let Some(value) = prop.value.clone() {
                            extra_props.insert(prop.type_name().clone(), value);
                        }
                    },
                    "condition" =>
                        if let Some(value) = prop.value.clone() {
                            conditions.push(self.make_term_short(&value));
                        },
                    "qualifier" =>
                        if let Some(value) = prop.value.clone() {
                            qualifiers.push(value);
                        },
                    "with" | "from" => {
                        if let Some(value) = prop.value.clone() {
                        let re = Regex::new(&db_prefix_patt).unwrap();
                            extra_props.insert(prop.type_name().clone(),
                                               re.replace_all(&value, ""));
                        }
                    },
                    _ => ()
                }
            }
            let mut maybe_genotype_short = None;
            let gene_uniquenames_vec: Vec<GeneUniquename> =
                match &feature.feat_type.name as &str {
                    "mRNA" => {
                        if let Some(gene_uniquename) =
                            self.genes_of_transcripts.get(&feature.uniquename) {
                                vec![gene_uniquename.clone()]
                            } else {
                                vec![]
                            }
                    },
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
                        maybe_genotype_short = Some(genotype_short.clone());
                        genotype_short.expressed_alleles.iter()
                            .map(|expressed_allele| {
                                expressed_allele.allele.gene.uniquename.clone()
                            })
                            .collect()
                    },
                    _ => {
                        if feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene" {
                            vec![feature.uniquename.clone()]
                        } else {
                            vec![]
                        }
                    }
                };

            let reference_short = self.make_reference_short(&publication.uniquename);

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
                let annotation_template = OntAnnotation {
                    id: feature_cvterm.feature_cvterm_id,
                    gene: Some(self.make_gene_short(&gene_uniquename)),
                    term: None,
                    reference: reference_short.clone(),
                    genotype: maybe_genotype_short.clone(),
                    with: extra_props_clone.remove("with"),
                    residue: extra_props_clone.remove("residue"),
                    gene_ex_props: gene_ex_props,
                    qualifiers: qualifiers.clone(),
                    evidence: extra_props_clone.remove("evidence"),
                    conditions: conditions.clone(),
                    extension: vec![],
                    is_not: feature_cvterm.is_not,
                };

                if cvterm.cv.name != "PomBase gene characterisation status" &&
                    cvterm.cv.name != "PomBase gene products" {
                        self.add_annotation(cvterm.borrow(), annotation_template)
                    }
            }
        }
    }

    fn store_ont_annotations(&mut self) {
        let mut term_gene_sets: HashMap<TermId, HashSet<GeneShort>> =
            HashMap::new();

        for ont_annotation in &self.all_ont_annotations {
            let term_short = ont_annotation.term.clone().unwrap();
            let gene_short = ont_annotation.gene.clone().unwrap();

            term_gene_sets.entry(term_short.termid.clone())
                .or_insert(HashSet::new())
                .insert(gene_short);
        }

        for ont_annotation in &mut self.all_ont_annotations {
            let mut new_term_short = ont_annotation.term.clone().unwrap();
            let term_gene_uniquenames = term_gene_sets.get(&new_term_short.termid).unwrap();
            new_term_short.gene_count = Some(term_gene_uniquenames.len());
            ont_annotation.term = Some(new_term_short);
        }

        self.all_ont_annotations.sort();

        for ont_annotation in &self.all_ont_annotations {
            let gene_uniquename = ont_annotation.gene.clone().unwrap().uniquename;
            let term = ont_annotation.term.clone().unwrap();
            let termid = term.termid;
            let cv_name = term.cv_name.clone();

            let rc_annotation = Rc::new(ont_annotation.clone());

            let mut gene_details = self.genes.get_mut(&gene_uniquename).unwrap();
            gene_details.annotations.entry(cv_name.clone())
                .or_insert(Vec::new()).push(rc_annotation.clone());

            if let Some(ref mut term_details) = self.terms.get_mut(&termid) {
                term_details.annotations.push(RelOntAnnotation {
                    rel_names: HashSet::new(),
                    annotation: rc_annotation.clone(),
                });

                if let Some(mut gene_short_vec) = term_gene_sets.remove(&termid) {
                    term_details.genes.append(&mut gene_short_vec.drain().collect());
                    term_details.genes.sort();
                }
            } else {
                panic!("missing termid: {}\n", &termid);
            }

            for condition_term_short in &rc_annotation.conditions {
                if let Some(ref mut condition_term_details) =
                    self.terms.get_mut(&condition_term_short.termid) {
                        condition_term_details.annotations.
                            push(RelOntAnnotation {
                                rel_names: HashSet::new(),
                                annotation: rc_annotation.clone(),
                            });
                    }
            }

            if let Some(reference) = ont_annotation.reference.clone() {
                let reference_uniquename = reference.uniquename.clone();
                let mut ref_details =
                    self.references.get_mut(&reference_uniquename).unwrap();
                ref_details.annotations.entry(cv_name).or_insert(Vec::new())
                    .push(rc_annotation.clone());
            }
        }
    }

    fn process_cvtermpath(&mut self) {
        let mut annotation_by_id: HashMap<i32, Rc<OntAnnotation>> = HashMap::new();
        let mut new_annotations: HashMap<TermId, HashMap<i32, HashSet<RelName>>> =
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

            if !DESCENDANT_REL_NAMES.contains(&&rel_term_name[0..]) {
                continue;
            }

            if let Some(term_details) = self.terms.get(&subject_termid) {
                let subject_annotations = &term_details.annotations;
                for rel_annotation in subject_annotations {
                    let RelOntAnnotation {
                        rel_names: _,
                        annotation: existing_annotation
                    } = rel_annotation.clone();

                    if !annotation_by_id.contains_key(&existing_annotation.id) {
                        annotation_by_id.insert(existing_annotation.id,
                                                existing_annotation.clone());
                    }
                    new_annotations.entry(object_termid.clone())
                        .or_insert(HashMap::new())
                        .entry(existing_annotation.id)
                        .or_insert(HashSet::new())
                        .insert(rel_term_name.clone());
                }
            } else {
                panic!("TermDetails not found for {}", &subject_termid);
            }
        }
        for (termid, annotations_map) in new_annotations.drain() {
            let mut term_details = self.terms.get_mut(&termid).unwrap();
            for (id, rel_names) in annotations_map {
                let annotation = annotation_by_id.get(&id).unwrap().clone();
                term_details.annotations.push(RelOntAnnotation {
                    rel_names: rel_names,
                    annotation: annotation
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
        let mut gene_summaries: HashSet<GeneShort> = HashSet::new();

        let gene_uniquenames: Vec<String> =
            self.genes.keys().map(|uniquename| uniquename.clone()).collect();

        for gene_uniquename in gene_uniquenames {
            gene_summaries.insert(self.make_gene_short(&gene_uniquename));
        }

        let mut term_summaries: HashSet<TermShort> = HashSet::new();

        let mut termid_genes: HashMap<TermId, HashSet<GeneUniquename>> = HashMap::new();
        let mut term_name_genes: HashMap<TermName, HashSet<GeneUniquename>> = HashMap::new();

        for (termid, term_details) in &self.terms {
            term_summaries.insert(self.make_term_short(&termid));
            for gene_short in &term_details.genes {
                termid_genes.entry(termid.clone())
                    .or_insert(HashSet::new())
                    .insert(gene_short.uniquename.clone());
                term_name_genes.entry(term_details.name.clone())
                    .or_insert(HashSet::new())
                    .insert(gene_short.uniquename.clone());
            }
        }

        SearchAPIMaps {
            gene_summaries: gene_summaries,
            termid_genes: termid_genes,
            term_name_genes: term_name_genes,
            term_summaries: term_summaries,
        }
    }

    pub fn get_web_data(&mut self) -> WebData {
        self.process_references();
        self.make_feature_rel_maps();
        self.process_features();
        self.process_props_from_feature_cvterms();
        self.process_allele_features();
        self.process_genotype_features();
        self.add_alleles_to_genotypes();
        self.process_cvterms();
        self.process_extension_cvterms();
        self.process_cvterm_rels();
        self.process_feature_synonyms();
        self.process_feature_cvterms();
        self.store_ont_annotations();
        self.process_cvtermpath();
        self.process_annotation_feature_rels();
//        self.set_term_short_use_counts();

        let mut web_data_terms: IdTermDetailsMap = HashMap::new();

        let search_api_maps = self.make_search_api_maps();

        for (termid, term_details) in self.terms.drain() {
            web_data_terms.insert(termid.clone(), Rc::new(term_details));
        }

        self.terms = HashMap::new();

        let mut used_terms: IdTermDetailsMap = HashMap::new();

        // remove terms with no annotation
        for (termid, term_details) in &web_data_terms {
            if term_details.annotations.len() > 0 {
                used_terms.insert(termid.clone(), term_details.clone());
            }
        }

        let metadata = self.make_metadata();

        WebData {
            genes: self.genes.clone(),
            terms: web_data_terms,
            used_terms: used_terms,
            metadata: metadata,
            references: self.references.clone(),

            search_api_maps: search_api_maps,
        }
    }
}
