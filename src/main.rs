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
const FEATURE_REL_ANNOTATIONS: [&'static str; 2] =
    ["interacts_physically", "interacts_genetically"];


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

struct WebDataBuild {
    genes: UniquenameGeneMap,
    transcripts: UniquenameTranscriptMap,
    genotypes: UniquenameGenotypeMap,
    alleles: UniquenameAlleleShortMap,
    terms: IdTermMap,

    genes_of_transcripts: HashMap<String, String>,
    genes_of_alleles: HashMap<String, String>,
    alleles_of_genotypes: HashMap<String, Vec<String>>,

    parts_of_extensions: HashMap<String, Vec<ExtPart>>
}

impl WebDataBuild {
    fn new() -> WebDataBuild {
        WebDataBuild {
            genes: HashMap::new(),
            transcripts: HashMap::new(),
            genotypes: HashMap::new(),
            alleles: HashMap::new(),
            terms: HashMap::new(),

            genes_of_transcripts: HashMap::new(),
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

    fn get_web_data(&mut self, raw: &Raw, organism_genus_species: &String) -> WebData {

        for feature_rel in raw.feature_relationships.iter() {
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

        //    type IdAnnotationMap = HashMap<i32, (FeatureAnnotation, TermAnnotation)>;

        //    let annotation_map: IdAnnotationMap = HashMap::new();

        for feat in raw.features.iter().filter(|&f| {
            let feature_org_genus_species = String::new() +
                &f.organism.genus + "_" + &f.organism.species;
            feature_org_genus_species == *organism_genus_species
        } ) {

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
                    self.genes.insert(feat.uniquename.clone(),
                                      GeneDetails {
                                          uniquename: feat.uniquename.clone(),
                                          name: feat.name.clone(),
                                          product: None,
                                          synonyms: vec![],
                                          feature_type: feat.feat_type.name.clone(),
                                          characterisation_status: None,
                                          location: location,
                                          cds_location: None,
                                          annotations: HashMap::new(),
                                          interaction_annotations: HashMap::new(),
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

        for feature_rel in raw.feature_relationships.iter() {
            let rel_name = &feature_rel.rel_type.name;
            let subject_uniquename = &feature_rel.subject.uniquename;
            let object_uniquename = &feature_rel.object.uniquename;

            for annotation_rel_type_name in FEATURE_REL_ANNOTATIONS.iter() {
                if rel_name == annotation_rel_type_name &&
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

                        let gene = {
                            let gene_details = self.genes.get(subject_uniquename).unwrap();
                            make_gene_short(gene_details)
                        };
                        let interactor = {
                            let interactor_details = self.genes.get(object_uniquename).unwrap();
                            make_gene_short(interactor_details)
                        };
                        let mut gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                        gene_details.interaction_annotations.entry(rel_name.clone()).or_insert(Vec::new()).push(
                            InteractionAnnotation {
                                gene: gene,
                                interactor: interactor,
                                evidence: evidence,
                                publication: None, // FIXME
                            });
                    }
            }
        }

        // a map from IDs of terms from the "PomBase annotation extension terms" cv
        // to a Vec of the details of each of the extension
        let mut parts_of_extensions: HashMap<String, Vec<ExtPart>> = HashMap::new();

        for cvterm in &raw.cvterms {
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

                        parts_of_extensions.entry(cvterm.termid())
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

        for cvterm_rel in &raw.cvterm_relationships {
            let subject_term = &cvterm_rel.subject;
            let object_term = &cvterm_rel.object;
            let rel_type = &cvterm_rel.rel_type;

            if subject_term.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                parts_of_extensions.entry(subject_term.termid())
                    .or_insert(Vec::new()).push(ExtPart {
                        rel_type_name: rel_type.name.clone(),
                        range_type: ExtRangeType::Term,
                        ext_range: object_term.termid(),
                    });
            }
        }

        for feature_synonym in raw.feature_synonyms.iter() {
            let feature = &feature_synonym.feature;
            let synonym = &feature_synonym.synonym;

            if let Some(ref mut gene_details) = self.genes.get_mut(&feature.uniquename) {
                gene_details.synonyms.push(SynonymDetails {
                    name: synonym.name.clone(),
                    synonym_type: synonym.synonym_type.name.clone()
                });
            }
        }

        for feature_cvterm in raw.feature_cvterms.iter() {
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
                if feature.feat_type.name == "mRNA" {
                    if let Some(gene_uniquename) =
                        self.genes_of_transcripts.get(&feature.uniquename) {
                            vec![gene_uniquename.clone()]
                        } else {
                            vec![]
                        }
                } else {
                    if feature.feat_type.name == "genotype" {
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
                    } else {
                        if feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene" {
                            vec![feature.uniquename.clone()]
                        } else {
                            vec![]
                        }
                    }
                };


            for gene_uniquename in &gene_uniquenames_vec {
                if cvterm.cv.name == "PomBase gene characterisation status" {
                    self.add_characterisation_status(&gene_uniquename, &cvterm.name);
                } else {
                    self.add_annotation_to_gene(&gene_uniquename, cvterm.borrow(),
                                                &evidence, &publication, &genotype_alleles);
                }
            }
        }

        // remove terms with no annotation
        let used_terms = self.terms.clone().into_iter()
            .filter(|&(_, ref t)| t.annotations.len() > 0)
            .collect();

        let mut db_creation_datetime = None;

        for chadoprop in &raw.chadoprops {
            if chadoprop.prop_type.name == "db_creation_datetime" {
                db_creation_datetime = chadoprop.value.clone();
            }
        }

        let metadata = Metadata {
            db_creation_datetime: db_creation_datetime.unwrap(),
            gene_count: self.genes.len(),
            term_count: self.terms.len(),
        };

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

fn write_gene_summary(output_dir: &str, genes: &IdGeneMap) {
    let gene_summaries =
        genes.values()
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

fn write_web_data(output_dir: &str, web_data: &WebData) {
    let s = serde_json::to_string(&web_data).unwrap();
    let file_name = String::new() + output_dir + "/all.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    write_gene_details(output_dir, &web_data.genes);
    write_gene_summary(output_dir, &web_data.genes);
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
    let mut web_data_build = WebDataBuild::new();
    let web_data = web_data_build.get_web_data(&raw, &organism_genus_species);

    write_web_data(&output_dir, &web_data);
}
