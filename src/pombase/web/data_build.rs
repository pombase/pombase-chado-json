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

type Termid = String;

struct TermidCvname {
    termid: Termid,
    cv_name: CvName,
}

pub struct WebDataBuild<'a> {
    raw: &'a Raw,

    genes: UniquenameGeneMap,
    transcripts: UniquenameTranscriptMap,
    genotypes: UniquenameGenotypeMap,
    alleles: UniquenameAlleleShortMap,
    terms: HashMap<TermId, TermDetails>,
    references: IdReferenceMap,

    genes_of_transcripts: HashMap<String, String>,
    transcripts_of_polypeptides: HashMap<String, String>,
    genes_of_alleles: HashMap<String, String>,
    alleles_of_genotypes: HashMap<String, Vec<String>>,

    // a map from IDs of terms from the "PomBase annotation extension terms" cv
    // to a Vec of the details of each of the extension
    parts_of_extensions: HashMap<String, Vec<ExtPart>>,

    base_term_of_extensions: HashMap<String, TermidCvname>,
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

            genes_of_transcripts: HashMap::new(),
            transcripts_of_polypeptides: HashMap::new(),
            genes_of_alleles: HashMap::new(),
            alleles_of_genotypes: HashMap::new(),

            parts_of_extensions: HashMap::new(),

            base_term_of_extensions: HashMap::new(),
        }
    }

    fn make_gene_short(&self, gene_uniquename: &str) -> GeneShort {
        if let Some(gene_details) = self.genes.get(gene_uniquename) {
            GeneShort {
                uniquename: gene_details.uniquename.clone(),
                name: gene_details.name.clone(),
                product: gene_details.product.clone(),
                synonyms: gene_details.synonyms.clone(),
            }
        } else {
            panic!("can't find GeneShort for gene uniquename {}", gene_uniquename);
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
        let mut gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        gene_details.product = Some(product.clone());
    }

    fn add_annotation(&mut self, gene_uniquename: &String,
                      cvterm: &Cvterm, evidence: &Option<Evidence>,
                      reference_opt: &Option<ReferenceShort>,
                      genotype_and_alleles: &Option<GenotypeAndAlleles>) {
        let (termid, cv_name) =
            match self.base_term_of_extensions.get(&cvterm.termid()) {
                Some(termid_cv_name) => (termid_cv_name.termid.clone(),
                                         termid_cv_name.cv_name.clone()),
                None => (cvterm.termid(), cvterm.cv.name.clone()),
            };

        let extension_parts =
            match self.parts_of_extensions.get(&cvterm.termid()) {
                Some(parts) => parts.clone(),
                None => vec![],
            };

        let term_short = self.make_term_short(&termid).clone();
        let gene_short = self.make_gene_short(&gene_uniquename).clone();

        let feature_annotation =
            FeatureAnnotation {
                term: term_short.clone(),
                extension: extension_parts.clone(),
                evidence: evidence.clone(),
                reference: reference_opt.clone(),
                genotype: genotype_and_alleles.clone(),
            };
        let mut gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        gene_details.annotations.entry(cv_name.clone())
            .or_insert(Vec::new()).push(feature_annotation);

        let term_annotation =
            TermAnnotation {
                term: term_short.clone(),
                gene: gene_short.clone(),
                evidence: evidence.clone(),
                reference: reference_opt.clone(),
                extension: extension_parts.clone(),
            };
        if let Some(ref mut term_details) = self.terms.get_mut(&termid) {
            term_details.annotations.entry(String::from("direct")).or_insert(Vec::new()).push(Rc::new(term_annotation));
        } else {
            panic!("missing termid: {}\n", &termid);
        }

        if let Some(reference) = reference_opt.clone() {
            let ref_annotation =
                ReferenceAnnotation {
                    gene: gene_short.clone(),
                    term: term_short.clone(),
                    evidence: evidence.clone(),
                    extension: extension_parts.clone(),
                    genotype: genotype_and_alleles.clone(),
                };
            let mut ref_details = self.references.get_mut(&reference.uniquename).unwrap();
            ref_details.annotations.entry(cv_name).or_insert(Vec::new()).push(ref_annotation);
        }
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
                        let entry = self.alleles_of_genotypes.entry(object_uniquename.clone());
                        entry.or_insert(Vec::new()).push(subject_uniquename.clone());
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

    fn process_feature(&mut self, feat: &Feature) {
        match &feat.feat_type.name as &str {
            "gene" | "pseudogene" => {
                let location = self.make_location(&feat);

                let organism = make_organism_short(&feat.organism);
                self.genes.insert(feat.uniquename.clone(),
                                  GeneDetails {
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
                                  });

            },
            // TODO: mRNA isn't the only transcript type
            "mRNA" => {
                self.transcripts.insert(feat.uniquename.clone(),
                                        TranscriptDetails {
                                            uniquename: feat.uniquename.clone(),
                                            name: feat.name.clone(),
                                        });
            },
            "genotype" => {
                self.genotypes.insert(feat.uniquename.clone(),
                                      GenotypeDetails {
                                          uniquename: feat.uniquename.clone(),
                                          name: feat.name.clone(),
                                          annotations: HashMap::new(),
                                      });
            },
            "allele" => {
                let gene_uniquename =
                    self.genes_of_alleles.get(&feat.uniquename).unwrap();
                self.alleles.insert(feat.uniquename.clone(),
                                    AlleleShort {
                                        uniquename: feat.uniquename.clone(),
                                        name: feat.name.clone(),
                                        gene_uniquename: gene_uniquename.clone(),
                                    });
            },
            _ => (),
        }
    }

    fn process_features(&mut self) {
        for feat in &self.raw.features {
            self.process_feature(&feat);
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
                                      annotations: HashMap::new(),
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
                    let base = TermidCvname {
                        termid: object_term.termid().clone(),
                        cv_name: object_term.cv.name.clone()
                    };
                    self.base_term_of_extensions.insert(subject_termid.clone(), base);
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

    fn process_feature_cvterms(&mut self) {
        for feature_cvterm in self.raw.feature_cvterms.iter() {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;
            let mut evidence: Option<String> = None;
            for prop in feature_cvterm.feature_cvtermprops.borrow().iter() {
                if (*prop).type_name() == "evidence" {
                    evidence = prop.value.clone();
                }
            }
            let reference_short =
                self.make_reference_short(&feature_cvterm.publication.uniquename);
            let mut genotype_alleles = None;
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
                        if let Some(allele_uniquenames_vec) =
                            self.alleles_of_genotypes.get(&feature.uniquename) {
                                let genotype_alleles_ret = GenotypeAndAlleles {
                                    alleles: allele_uniquenames_vec.iter()
                                        .map(|allele_uniquename| {
                                            self.alleles.get(allele_uniquename).unwrap().clone()
                                        })
                                        .collect::<Vec<_>>()
                                };
                                genotype_alleles = Some(genotype_alleles_ret);
                                allele_uniquenames_vec.iter()
                                    .map(|allele_uniquename|
                                         {
                                             let gene_uniquename =
                                                 self.genes_of_alleles.get(allele_uniquename).unwrap();
                                             gene_uniquename.clone()
                                         })
                                    .collect()
                            } else {
                                vec![]
                            }
                    },
                    _ => {
                        if feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene" {
                            vec![feature.uniquename.clone()]
                        } else {
                            vec![]
                        }
                    }
                };

            for gene_uniquename in &gene_uniquenames_vec {
                match &cvterm.cv.name as &str {
                    "PomBase gene characterisation status" => self.add_characterisation_status(&gene_uniquename, &cvterm.name),
                    "PomBase gene products" => self.add_gene_product(&gene_uniquename, &cvterm.name),
                    _ => self.add_annotation(&gene_uniquename, cvterm.borrow(),
                                             &evidence, &reference_short, &genotype_alleles)
                }
            }
        }
    }

    fn process_cvtermpath(&mut self) {
        let mut new_annotations: HashMap<Termid, TermAnnotationMap> = HashMap::new();

        for cvtermpath in &self.raw.cvtermpaths {
            let subject_term = &cvtermpath.subject;
            let subject_termid = subject_term.termid();
            let object_term = &cvtermpath.object;
            let object_termid = object_term.termid();

            // cope with multiple paths between terms
            let mut seen_pairs: HashSet<(Termid, Termid)> = HashSet::new();

            let key = (subject_termid.clone(), object_termid.clone());
            if let Some(_) = seen_pairs.get(&key) {
                continue;
            }

            let distance = cvtermpath.pathdistance.unwrap();

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

            seen_pairs.insert((subject_termid.clone(), object_termid.clone()));
            if let Some(term_details) = self.terms.get(&subject_termid) {
                let subject_annotations = &term_details.annotations;
                for (_, annotations) in subject_annotations {
                    for annotation in annotations {
                        let new_annotation = annotation.clone();
                        let key = &format!("{}::{}", rel_term_name, distance);
                        new_annotations.entry(object_termid.clone())
                            .or_insert(HashMap::new())
                            .entry(key.clone())
                            .or_insert(Vec::new()).push(new_annotation);
                    }
                }
            } else {
                panic!("TermDetails not found for {}", &subject_termid);
            }
        }
        for (termid, mut annotations) in new_annotations.drain() {
            let mut term_details = self.terms.get_mut(&termid).unwrap();
            for (key, annotation_vec) in annotations.drain() {
                let mut vec_clone = annotation_vec.clone();
                term_details.annotations.entry(key)
                    .or_insert(annotation_vec).append(&mut vec_clone);
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

    fn set_term_short_use_counts(&mut self) {
        let mut seen_genes: HashMap<Termid, HashSet<String>> = HashMap::new();

        for (termid, term_details) in &self.terms {
            let mut term_seen_genes: HashSet<String> = HashSet::new();
            for (_, annotation_vec) in &term_details.annotations {
                for annotation in annotation_vec {
                    term_seen_genes.insert(annotation.gene.uniquename.clone());
                }
            }
            seen_genes.insert(termid.clone(), term_seen_genes);
        }

        for (_, gene_details) in &mut self.genes {
            for (_, feat_annotations) in &mut gene_details.annotations {
                for mut feat_annotation in feat_annotations.iter_mut() {
                    feat_annotation.term.gene_count =
                        Some(seen_genes.get(&feat_annotation.term.termid).unwrap().len());
                }
            }
        }

        for (_, ref_details) in &mut self.references {
            for (_, ref_annotations) in &mut ref_details.annotations {
                for ref_annotation in ref_annotations {
                    ref_annotation.term.gene_count =
                        Some(seen_genes.get(&ref_annotation.term.termid).unwrap().len());
                }
            }
        }

        let mut gene_short_map: HashMap<Termid, Vec<GeneShort>> = HashMap::new();

        for (termid, term_seen_genes) in seen_genes {
            for gene_uniquename in term_seen_genes {
                let gene_short = self.make_gene_short(&gene_uniquename);
                gene_short_map.entry(termid.clone())
                    .or_insert(Vec::new()).push(gene_short);
            }
        }

        for (ref termid, ref mut gene_short_vec) in gene_short_map.drain() {
            let mut term_details = self.terms.get_mut(termid).unwrap();
            term_details.genes.append(gene_short_vec);
        }
    }

    pub fn get_web_data(&mut self) -> WebData {
        self.process_references();
        self.make_feature_rel_maps();
        self.process_features();
        self.process_cvterms();
        self.process_extension_cvterms();
        self.process_cvterm_rels();
        self.process_feature_synonyms();
        self.process_feature_cvterms();
        self.process_cvtermpath();
        self.process_annotation_feature_rels();
        self.set_term_short_use_counts();

        let mut web_data_terms: IdTermDetailsMap = HashMap::new();

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

        let mut gene_summaries: IdGeneShortMap = HashMap::new();

        let gene_uniquenames: Vec<String> =
            self.genes.keys().map(|uniquename| uniquename.clone()).collect();

        for gene_uniquename in gene_uniquenames {
            gene_summaries.insert(gene_uniquename.clone(), self.make_gene_short(&gene_uniquename));
        }

        WebData {
            genes: self.genes.clone(),
            gene_summaries: gene_summaries,
            terms: web_data_terms,
            used_terms: used_terms,
            metadata: metadata,
            references: self.references.clone(),
        }
    }
}
