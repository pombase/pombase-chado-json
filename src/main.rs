#![feature(plugin, proc_macro)]
extern crate postgres;

#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_json;

extern crate getopts;

use postgres::{Connection, TlsMode};

use std::env;
use getopts::Options;
use std::fs;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::collections::HashMap;
use std::rc::Rc;
use std::borrow::Borrow;

mod pombase;

use pombase::web::data::*;
use pombase::db::*;

const POMBASE_ANN_EXT_TERM_CV_NAME: &'static str = "PomBase annotation extension terms";
const ANNOTATION_EXT_REL_PREFIX: &'static str = "annotation_extension_relation-";
enum FeatureRelAnnotationType {
    Interaction,
    Ortholog,
    Paralog,
}
struct FeatureRelConfig {
    rel_type_name: &'static str,
    annotation_type: FeatureRelAnnotationType,
}
const FEATURE_REL_CONFIGS: [FeatureRelConfig; 4] =
    [
        FeatureRelConfig {
            rel_type_name: "interacts_physically",
            annotation_type: FeatureRelAnnotationType::Interaction,
        },
        FeatureRelConfig {
            rel_type_name: "interacts_genetically",
            annotation_type: FeatureRelAnnotationType::Interaction,
        },
        FeatureRelConfig {
            rel_type_name: "orthologous_to",
            annotation_type: FeatureRelAnnotationType::Ortholog,
        },
        FeatureRelConfig {
            rel_type_name: "paralogous_to",
            annotation_type: FeatureRelAnnotationType::Paralog,
        },
    ];


fn make_publication_short(rc_publication: Rc<pombase::db::Publication>) -> Option<PublicationShort> {
    if rc_publication.uniquename == "null" {
        None
    } else {
        Some(PublicationShort {
            uniquename: rc_publication.uniquename.clone(),
            title: rc_publication.title.clone(),
            citation: rc_publication.miniref.clone(),
        })
    }
}

#[derive(Serialize)]
struct WebData {
    genes: IdGeneMap,
    terms: IdTermMap,
    used_terms: IdTermMap,
    metadata: Metadata,
}

struct WebDataBuild<'a> {
    raw: &'a Raw,

    genes: UniquenameGeneMap,
    transcripts: UniquenameTranscriptMap,
    genotypes: UniquenameGenotypeMap,
    alleles: UniquenameAlleleShortMap,
    terms: IdTermMap,

    genes_of_transcripts: HashMap<String, String>,
    transcripts_of_polypeptides: HashMap<String, String>,
    genes_of_alleles: HashMap<String, String>,
    alleles_of_genotypes: HashMap<String, Vec<String>>,

    // a map from IDs of terms from the "PomBase annotation extension terms" cv
    // to a Vec of the details of each of the extension
    parts_of_extensions: HashMap<String, Vec<ExtPart>>
}

impl <'a> WebDataBuild<'a> {
    fn new(raw: &'a Raw) -> WebDataBuild<'a> {
        WebDataBuild {
            raw: raw,

            genes: HashMap::new(),
            transcripts: HashMap::new(),
            genotypes: HashMap::new(),
            alleles: HashMap::new(),
            terms: HashMap::new(),

            genes_of_transcripts: HashMap::new(),
            transcripts_of_polypeptides: HashMap::new(),
            genes_of_alleles: HashMap::new(),
            alleles_of_genotypes: HashMap::new(),

            parts_of_extensions: HashMap::new(),
        }
    }

    fn make_term_short(&self, cvterm: &Cvterm) -> (TermShort, CvName, Vec<ExtPart>) {
        let mut extension: Vec<ExtPart> = vec![];

        let term =
            if cvterm.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                if let Some(ext_parts) = self.parts_of_extensions.get(&cvterm.termid()) {
                    let mut base_term_details_opt = None;

                    for ext_part in ext_parts {
                        if ext_part.rel_type_name == "is_a" {
                            base_term_details_opt = self.terms.get(&ext_part.ext_range);
                        } else {
                            extension.push(ext_part.clone());
                        }
                    }

                    if let Some(term_details) = base_term_details_opt {
                        term_details
                    } else {
                        self.terms.get(&cvterm.termid()).unwrap()
                    }
                } else {
                    self.terms.get(&cvterm.termid()).unwrap()
                }
            } else {
                self.terms.get(&cvterm.termid()).unwrap()
            };

        (
            TermShort {
                name: term.name.clone(),
                termid: term.termid.clone(),
                is_obsolete: term.is_obsolete
            },
            term.cv_name.clone(),
            extension
        )
    }

    fn add_characterisation_status(&mut self, gene_uniquename: &String, cvterm_name: &String) {
        let mut gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        gene_details.characterisation_status = Some(cvterm_name.clone());
    }

    fn add_gene_product(&mut self, gene_uniquename: &String, product: &String) {
        let mut gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        gene_details.product = Some(product.clone());
    }

    fn add_annotation_to_gene(&mut self, gene_uniquename: &String,
                              cvterm: &Cvterm, evidence: &Option<Evidence>,
                              publication: &Option<PublicationShort>,
                              genotype_and_alleles: &Option<GenotypeAndAlleles>) {
        let (term, cv_name, extension) = self.make_term_short(&cvterm);
        let mut gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        let feature_annotation =
            FeatureAnnotation {
                term: term,
                extension: extension.clone(),
                evidence: evidence.clone(),
                publication: publication.clone(),
                genotype: genotype_and_alleles.clone(),
            };
        gene_details.annotations.entry(cv_name).or_insert(Vec::new()).push(feature_annotation);
        let term_annotation =
            TermAnnotation {
                gene: make_gene_short(&gene_details),
                evidence: evidence.clone(),
                publication: publication.clone(),
                extension: extension.clone(),
            };
        let termid = cvterm.termid();
        let term_details = self.terms.get_mut(&termid).unwrap();
        term_details.annotations.push(term_annotation);
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

    fn process_features(&mut self) {
        for feat in &self.raw.features {
            match &feat.feat_type.name as &str {
                "gene" | "pseudogene" => {
                    let feature_locs = feat.featurelocs.borrow();
                    let location =
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
                        };
                    let organism = OrganismDetails {
                        genus: feat.organism.genus.clone(),
                        species: feat.organism.species.clone(),
                    };
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
    }

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

                        for prop in feature_rel.feature_relationshipprops.borrow().iter() {
                            if prop.prop_type.name == "evidence" {
                                evidence = prop.value.clone();
                            }
                        }

                        let evidence_clone = evidence.clone();

                        let gene = {
                            let gene_details = self.genes.get(subject_uniquename).unwrap();
                            make_gene_short(gene_details)
                        };
                        let gene_clone = gene.clone();
                        let other_gene = {
                            let other_gene_details = self.genes.get(object_uniquename).unwrap();
                            make_gene_short(other_gene_details)
                        };
                        {
                            let mut gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                            match rel_config.annotation_type {
                                FeatureRelAnnotationType::Interaction =>
                                    gene_details.interaction_annotations.entry(rel_name.clone()).or_insert(Vec::new()).push(
                                        InteractionAnnotation {
                                            gene: gene,
                                            interactor: other_gene,
                                            evidence: evidence,
                                            publication: None, // FIXME
                                        }),
                                FeatureRelAnnotationType::Ortholog =>
                                    gene_details.ortholog_annotations.push(
                                        OrthologAnnotation {
                                            ortholog: other_gene,
                                            evidence: evidence,
                                            publication: None, // FIXME
                                        }),
                                FeatureRelAnnotationType::Paralog =>
                                    gene_details.paralog_annotations.push(
                                        ParalogAnnotation {
                                            paralog: other_gene,
                                            evidence: evidence,
                                            publication: None, // FIXME
                                        }),
                            }
                        }
                        {
                            let mut other_gene_details = self.genes.get_mut(object_uniquename).unwrap();
                            match rel_config.annotation_type {
                                FeatureRelAnnotationType::Interaction => {},
                                FeatureRelAnnotationType::Ortholog =>
                                    other_gene_details.ortholog_annotations.push(
                                        OrthologAnnotation {
                                            ortholog: gene_clone,
                                            evidence: evidence_clone,
                                            publication: None, // FIXME
                                        }),
                                FeatureRelAnnotationType::Paralog =>
                                    other_gene_details.paralog_annotations.push(
                                        ParalogAnnotation {
                                            paralog: gene_clone,
                                            evidence: evidence_clone,
                                            publication: None, // FIXME
                                        }),
                            }
                        }
                    }
            }
        }
    }

    fn process_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                    if (*cvtermprop).prop_type.name.starts_with(ANNOTATION_EXT_REL_PREFIX) {
                        let ext_rel_name: &str = &(*cvtermprop).prop_type.name[ANNOTATION_EXT_REL_PREFIX.len()..];
                        let ext_range = (*cvtermprop).value.clone();
                        let range_type = if ext_range.starts_with("SP") {
                            ExtRangeType::Gene
                        } else {
                            ExtRangeType::Misc
                        };

                        self.parts_of_extensions.entry(cvterm.termid())
                            .or_insert(Vec::new()).push(ExtPart {
                                rel_type_name: String::from(ext_rel_name),
                                range_type: range_type,
                                ext_range: ext_range,
                            });
                    }
                }
            }

            self.terms.insert(cvterm.termid(),
                              TermDetails {
                                  name: cvterm.name.clone(),
                                  cv_name: cvterm.cv.name.clone(),
                                  termid: cvterm.termid(),
                                  definition: cvterm.definition.clone(),
                                  is_obsolete: cvterm.is_obsolete,
                                  annotations: vec![],
                              });
        }
    }

    fn process_cvterm_rels(&mut self) {
        for cvterm_rel in &self.raw.cvterm_relationships {
            let subject_term = &cvterm_rel.subject;
            let object_term = &cvterm_rel.object;
            let rel_type = &cvterm_rel.rel_type;

            if subject_term.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                self.parts_of_extensions.entry(subject_term.termid())
                    .or_insert(Vec::new()).push(ExtPart {
                        rel_type_name: rel_type.name.clone(),
                        range_type: ExtRangeType::Term,
                        ext_range: object_term.termid(),
                    });
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
            let publication = make_publication_short(feature_cvterm.publication.clone());
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
                    _ => self.add_annotation_to_gene(&gene_uniquename, cvterm.borrow(),
                                                     &evidence, &publication, &genotype_alleles)
                }
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

        Metadata {
            db_creation_datetime: db_creation_datetime.unwrap(),
            gene_count: self.genes.len(),
            term_count: self.terms.len(),
        }
    }

    fn get_web_data(&mut self) -> WebData {
        self.make_feature_rel_maps();
        self.process_features();
        self.process_cvterms();
        self.process_cvterm_rels();
        self.process_feature_synonyms();
        self.process_feature_cvterms();
        self.process_annotation_feature_rels();

        // remove terms with no annotation
        let used_terms = self.terms.clone().into_iter()
            .filter(|&(_, ref t)| t.annotations.len() > 0)
            .collect();

        let metadata = self.make_metadata();

        WebData {
            genes: self.genes.clone(),
            terms: self.terms.clone(),
            used_terms: used_terms,
            metadata: metadata,
        }
    }
}

fn write_gene_details(output_dir: &str, genes: &IdGeneMap) {
    let s = serde_json::to_string(&genes).unwrap();
    let file_name = String::new() + &output_dir + "/genes.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    for (gene_uniquename, gene_details) in genes {
        let s = serde_json::to_string(&gene_details).unwrap();
        let file_name = String::new() + &output_dir + "/gene/" + &gene_uniquename + ".json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }
}

fn make_gene_short(gene_details: &GeneDetails) -> GeneShort {
    GeneShort {
        uniquename: gene_details.uniquename.clone(),
        name: gene_details.name.clone(),
        product: gene_details.product.clone(),
        synonyms: gene_details.synonyms.clone(),
    }
}

fn write_gene_summary(output_dir: &str, genes: &IdGeneMap, organism_genus_species: &str) {

    let gene_summaries =
        genes.values()
        .filter(|gene_details| {
            let feature_org_genus_species = String::new() +
                &gene_details.organism.genus + "_" + &gene_details.organism.species;
            feature_org_genus_species == organism_genus_species
        })
        .map(|gene_details| make_gene_short(&gene_details))
        .collect::<Vec<GeneShort>>();
    let s = serde_json::to_string(&gene_summaries).unwrap();
    let file_name = String::new() + &output_dir + "/gene_summaries.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");
}

fn write_terms(output_dir: &str, terms: &IdTermMap) {
    let s = serde_json::to_string(&terms).unwrap();
    let file_name = String::new() + &output_dir + "/terms.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    for (termid, term_details) in terms {
        let s = serde_json::to_string(&term_details).unwrap();
        let file_name = String::new() + &output_dir + "/term/" + &termid + ".json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }
}

fn write_metadata(output_dir: &str, metadata: &Metadata) {
    let s = serde_json::to_string(&metadata).unwrap();
    let file_name = String::new() + &output_dir + "/metadata.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");
}

fn write_web_data(output_dir: &str, web_data: &WebData, organism_genus_species: &str) {
    let s = serde_json::to_string(&web_data).unwrap();
    let file_name = String::new() + output_dir + "/all.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    write_gene_details(output_dir, &web_data.genes);
    write_gene_summary(output_dir, &web_data.genes, organism_genus_species);
    write_terms(output_dir, &web_data.terms);
    write_metadata(output_dir, &web_data.metadata);

    println!("wrote {} genes", web_data.genes.len());
    println!("wrote {} terms", web_data.terms.len());
}

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const VERSION: &'static str = env!("CARGO_PKG_VERSION");

fn main() {
    print!("{} v{}\n", PKG_NAME, VERSION);

    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();

    opts.reqopt("c", "connection-string",
                "PostgresSQL connection string like: postgres://user:pass@host/db_name",
                "CONN_STR");
    opts.reqopt("d", "output-directory",
                "Destination directory for JSON output", "DIR");
    opts.reqopt("O", "organism",
                "Only output genes from this organism (eg. 'Schizosaccharomyces_pombe')", "GENUS_SPECIES");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("Invalid options\n{}", f)
    };

    let connection_string = matches.opt_str("c").unwrap();
    let output_dir = matches.opt_str("d").unwrap();
    let organism_genus_species = matches.opt_str("O").unwrap();

    let subdirs = vec!["gene", "term"];

    for subdir in &subdirs {
        let dir = String::new() + &output_dir + "/" + subdir;
        fs::create_dir_all(&dir).unwrap_or_else(|why| {
            println!("Creating output directory failed: {:?}", why.kind());
        });
    }

    let conn = Connection::connect(connection_string.as_str(), TlsMode::None).unwrap();
    let raw = Raw::new(&conn);
    let mut web_data_build = WebDataBuild::new(&raw);
    let web_data = web_data_build.get_web_data();

    write_web_data(&output_dir, &web_data, &organism_genus_species);
}

#[allow(dead_code)]
fn make_test_cvterm_dbxref(cvterms: &mut Vec<Rc<Cvterm>>, dbxrefs: &mut Vec<Rc<Dbxref>>,
                           cv: &Rc<Cv>, db: &Rc<Db>, cvterm_name: &str, accession: &str)
                           -> Rc<Cvterm> {
    use std::cell::RefCell;
    let dbxref = Rc::new(Dbxref {
        accession: String::from(accession),
        db: db.clone(),
    });
    let cvterm = Rc::new(Cvterm {
        name: String::from(cvterm_name),
        definition: None,
        cv: cv.clone(),
        dbxref: dbxref.clone(),
        is_obsolete: false,
        cvtermprops: RefCell::new(vec![]),
    });
    cvterms.push(cvterm.clone());
    dbxrefs.push(dbxref.clone());
    cvterm
}

#[allow(dead_code)]
fn make_test_cv(cvs: &mut Vec<Rc<Cv>>, cv_name: &str) -> Rc<Cv> {
    let cv = Rc::new(Cv {
        name: String::from(cv_name),
    });
    cvs.push(cv.clone());
    cv
}

#[allow(dead_code)]
fn make_test_db(dbs: &mut Vec<Rc<Db>>, db_name: &str) -> Rc<Db> {
    let db = Rc::new(Db {
        name: String::from(db_name),
    });
    dbs.push(db.clone());
    db
}

#[allow(dead_code)]
fn make_test_feature(features: &mut Vec<Rc<Feature>>, organism: &Rc<Organism>,
                     type_cvterm: &Rc<Cvterm>, uniquename: &str, name: Option<String>)
                     -> Rc<Feature> {
    use std::cell::RefCell;
    let feature = Rc::new(Feature {
        organism: organism.clone(),
        uniquename: String::from(uniquename),
        name: name,
        feat_type: type_cvterm.clone(),
        featurelocs: RefCell::new(vec![]),
        featureprops: RefCell::new(vec![]),
    });
    features.push(feature.clone());
    feature
}

#[allow(dead_code)]
fn make_test_feature_cvterm(feature_cvterms: &mut Vec<Rc<FeatureCvterm>>,
                            feature: &Rc<Feature>, cvterm: &Rc<Cvterm>,
                            publication: &Rc<Publication>)
                     -> Rc<FeatureCvterm> {
    use std::cell::RefCell;
    let feature_cvterm = Rc::new(FeatureCvterm {
        feature_cvterm_id: 0,
        feature: feature.clone(),
        cvterm: cvterm.clone(),
        publication: publication.clone(),
        is_not: false,
        feature_cvtermprops: RefCell::new(vec![]),
    });
    feature_cvterms.push(feature_cvterm.clone());
    feature_cvterm
}

#[allow(dead_code)]
fn make_test_feature_rel(feature_relationships: &mut Vec<Rc<FeatureRelationship>>,
                         subject: &Rc<Feature>, rel: &Rc<Cvterm>, object: &Rc<Feature>) {
    use std::cell::RefCell;
    let rel = Rc::new(FeatureRelationship {
        subject: subject.clone(),
        rel_type: rel.clone(),
        object: object.clone(),
        feature_relationshipprops: RefCell::new(vec![]),
    });
    feature_relationships.push(rel);
}

#[allow(dead_code)]
fn make_test_cvterm_rel(cvterm_relationships: &mut Vec<Rc<CvtermRelationship>>,
                        subject: &Rc<Cvterm>, rel_type: &Rc<Cvterm>, object: &Rc<Cvterm>) {
    let rel = Rc::new(CvtermRelationship {
        subject: subject.clone(),
        rel_type: rel_type.clone(),
        object: object.clone(),
    });
    cvterm_relationships.push(rel);
}

#[allow(dead_code)]
fn get_test_raw() -> Raw {
    use std::cell::RefCell;
    let mut feature_relationships: Vec<Rc<FeatureRelationship>> = vec![];
    let mut cvterm_relationships: Vec<Rc<CvtermRelationship>> = vec![];
    let mut cvs: Vec<Rc<Cv>> = vec![];
    let mut dbs: Vec<Rc<Db>> = vec![];
    let mut cvterms: Vec<Rc<Cvterm>> = vec![];
    let mut dbxrefs: Vec<Rc<Dbxref>> = vec![];
    let mut features: Vec<Rc<Feature>> = vec![];
    let mut feature_cvterms: Vec<Rc<FeatureCvterm>> = vec![];

    let pombe_organism =
        Rc::new(Organism{
            genus: String::from("Schizosaccharomyces"),
            species: String::from("pombe"),
            abbreviation: String::from("Spombe"),
            common_name: String::from("pombe"),
        });

    let bp_cv = make_test_cv(&mut cvs, "biological_process");
    let fypo_cv = make_test_cv(&mut cvs, "fission_yeast_phenotype");
    let extension_cv = make_test_cv(&mut cvs, POMBASE_ANN_EXT_TERM_CV_NAME);
    let relations_cv = make_test_cv(&mut cvs, "relations");
    let pombase_relations_cv = make_test_cv(&mut cvs, "pombase_relations");
    let fypo_ext_relations_cv = make_test_cv(&mut cvs, "fypo_extension_relations");
    let sequence_cv = make_test_cv(&mut cvs, "sequence");
    let chadoprop_types_cv = make_test_cv(&mut cvs, "PomBase chadoprop types");
    let pub_type_cv = make_test_cv(&mut cvs, "PomBase publication types");

    let pbo_db = make_test_db(&mut dbs, "PBO");
    let go_db = make_test_db(&mut dbs, "GO");
    let fypo_db = make_test_db(&mut dbs, "FYPO");
    let fypo_ext_db = make_test_db(&mut dbs, "FYPO_EXT");
    let obo_rel_db = make_test_db(&mut dbs, "OBO_REL");
    let so_db = make_test_db(&mut dbs, "SO");

    let db_creation_datetime_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &chadoprop_types_cv, &pbo_db,
                                "db_creation_datetime", "0017491");
    let is_a_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &relations_cv, &obo_rel_db,
                                "is_a", "is_a");
    let part_of_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &relations_cv, &obo_rel_db,
                                "part_of", "0000050");
    let instance_of_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &pombase_relations_cv, &pbo_db,
                                "instance_of", "999990000");
    let derives_from_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &relations_cv, &obo_rel_db,
                                "derives_from", "0001000");
    let has_expressivity_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &pombase_relations_cv, &pbo_db,
                                "has_expressivity", "999990001");
    let medium_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &fypo_ext_relations_cv, &fypo_ext_db,
                                "medium", "0000002");
    let gene_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &sequence_cv, &so_db,
                                "gene", "0000704");
    let mrna_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &sequence_cv, &so_db,
                                "mRNA", "0000234");
    let genotype_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &sequence_cv, &so_db,
                                "genotype", "0001027");
    let allele_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &sequence_cv, &so_db,
                                "allele", "0001023");
    let polypeptide_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &sequence_cv, &so_db,
                                "polypeptide", "0000104");
    let chromosome_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &sequence_cv, &so_db,
                                "chromosome", "0000340");
    let paper_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &pub_type_cv, &pbo_db,
                                "paper", "0000044");

    let go0031030_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &bp_cv, &go_db, "negative regulation of septation initiation signaling",
                                "0031030");
    let pbo0022440_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &extension_cv, &pbo_db,
                                "decreased cell population growth at high temperature [has_expressivity] medium",
                                "0000082");
    make_test_cvterm_rel(&mut cvterm_relationships,
                         &pbo0022440_cvterm, &has_expressivity_cvterm, &medium_cvterm);

    let fypo0000082_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &fypo_cv, &fypo_db,
                                "decreased cell population growth at high temperature",
                                "0000082");
    make_test_cvterm_rel(&mut cvterm_relationships, &pbo0022440_cvterm, &is_a_cvterm, &fypo0000082_cvterm);

    let chadoprops = vec![Rc::new(Chadoprop {
        prop_type: db_creation_datetime_cvterm,
        value: Some(String::from("2016-10-17 03:41:56")),
    })];

    let chr_1 = make_test_feature(&mut features, &pombe_organism,
                                  &chromosome_cvterm, "chromosome_1", None);
    let chr_3 = make_test_feature(&mut features, &pombe_organism,
                                  &chromosome_cvterm, "chromosome_3", None);

    let cdc16_gene = make_test_feature(&mut features, &pombe_organism, &gene_cvterm,
                                      "SPAC6F6.08c", Some(String::from("cdc16")));
    let cdc16_allele1 = make_test_feature(&mut features, &pombe_organism, &gene_cvterm,
                                          "SPAC6F6.08c-allele1", Some(String::from("cdc16-116")));
    make_test_feature_rel(&mut feature_relationships,
                          &cdc16_allele1, &instance_of_cvterm, &cdc16_gene);

    let par1_gene = make_test_feature(&mut features, &pombe_organism, &gene_cvterm,
                                      "SPCC188.02", Some(String::from("par1")));
    let par1_delta_allele = make_test_feature(&mut features, &pombe_organism, &allele_cvterm,
                                              "SPCC188.02-allele1", Some(String::from("par1delta")));
    make_test_feature_rel(&mut feature_relationships,
                          &par1_delta_allele, &instance_of_cvterm, &par1_gene);

    let genotype1 = make_test_feature(&mut features, &pombe_organism,
                                      &genotype_cvterm, "test-genotype1", None);
    make_test_feature_rel(&mut feature_relationships,
                          &par1_delta_allele, &part_of_cvterm, &genotype1);
    make_test_feature_rel(&mut feature_relationships,
                          &cdc16_allele1, &part_of_cvterm, &genotype1);

    let par1_mrna = make_test_feature(&mut features, &pombe_organism,
                                      &mrna_cvterm, "SPCC188.02.1", None);
    make_test_feature_rel(&mut feature_relationships, &par1_mrna, &part_of_cvterm, &par1_gene);

    let par1_polypeptide = make_test_feature(&mut features, &pombe_organism, &polypeptide_cvterm,
                                             "SPCC188.02:pep", None);
    make_test_feature_rel(&mut feature_relationships,
                          &par1_polypeptide, &derives_from_cvterm, &par1_mrna);

    let publication = Rc::new(Publication {
        uniquename: String::from("PMID:11707284"),
        title: Some(String::from("The protein phosphatase 2A B'-regulatory subunit par1p is implicated in regulation of the S. pombe septation initiation network.")),
        pub_type: paper_cvterm,
        miniref: Some(String::from("FEBS Lett. 2001 Nov 9;508(1):136-42")),
        publicationprops: RefCell::new(vec![]),
    });

    make_test_feature_cvterm(&mut feature_cvterms, &par1_mrna, &go0031030_cvterm, &publication);
    make_test_feature_cvterm(&mut feature_cvterms, &genotype1, &pbo0022440_cvterm, &publication);

    par1_gene.featurelocs.borrow_mut().push(Rc::new(Featureloc {
        feature: par1_gene.clone(),
        fmin: 1475493,
        fmax: 1478103,
        strand: 1,
        srcfeature: chr_3.clone(),
    }));

    cdc16_gene.featurelocs.borrow_mut().push(Rc::new(Featureloc {
        feature: par1_gene.clone(),
        fmin: 2746666,
        fmax: 2748180,
        strand: -1,
        srcfeature: chr_1.clone(),
    }));

    Raw {
        organisms: vec![pombe_organism],
        cvs: cvs,
        dbs: dbs,
        dbxrefs: dbxrefs,
        cvterms: cvterms,
        cvtermprops: vec![],
        cvtermpaths: vec![],
        cvterm_relationships: cvterm_relationships,
        publications: vec![publication],
        synonyms: vec![],
        features: features,
        feature_synonyms: vec![],
        featureprops: vec![],
        featurelocs: vec![],
        feature_cvterms: feature_cvterms,
        feature_cvtermprops: vec![],
        feature_relationships: feature_relationships,
        chadoprops: chadoprops,
    }
}

#[allow(dead_code)]
fn get_test_web_data() -> WebData {
    let raw = get_test_raw();

    let mut web_data_build = WebDataBuild::new(&raw);
    web_data_build.get_web_data()
}

#[test]
fn test_gene_details() {
    let web_data = get_test_web_data();

    assert_eq!(web_data.genes.len(), 3);
    let par1_gene = web_data.genes.get("SPCC188.02").unwrap().clone();
    assert_eq!(par1_gene.uniquename, "SPCC188.02");
    assert_eq!(par1_gene.name.unwrap(), "par1");
    assert_eq!(par1_gene.annotations.len(), 2);

    if par1_gene.annotations.get(POMBASE_ANN_EXT_TERM_CV_NAME).is_some() {
        panic!("extension cv shouldn't be in the annotations");
    }
}
