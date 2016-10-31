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

mod pombase;

use pombase::web::data::*;
use pombase::db::*;

fn make_term_short(rc_cvterm: Rc<pombase::db::Cvterm>) -> TermShort {
    TermShort {
        name: rc_cvterm.name.clone(),
        termid: rc_cvterm.termid(),
        is_obsolete: rc_cvterm.is_obsolete
    }
}

fn make_gene_short(gene_details: &GeneDetails) -> GeneShort {
    GeneShort {
        uniquename: gene_details.uniquename.clone(),
        name: gene_details.name.clone(),
    }
}

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

fn get_web_data(raw: &Raw, organism_genus_species: &String) -> WebData {
    let mut genes: UniquenameGeneMap = HashMap::new();
    let mut transcripts: UniquenameTranscriptMap = HashMap::new();
    let mut genotypes: UniquenameGenotypeMap = HashMap::new();
    let mut terms: IdTermMap = HashMap::new();

//    type IdAnnotationMap = HashMap<i32, (FeatureAnnotation, TermAnnotation)>;

//    let annotation_map: IdAnnotationMap = HashMap::new();

    for feat in raw.features.iter().filter(|&f| {
        let feature_org_genus_species = String::new() +
            &f.organism.genus + "_" + &f.organism.species;
        feature_org_genus_species == *organism_genus_species
    } ) {

        match &feat.feat_type.name as &str {
            "gene" => {
                genes.insert(feat.uniquename.clone(),
                             GeneDetails {
                                 uniquename: feat.uniquename.clone(),
                                 name: feat.name.clone(),
                                 annotations: HashMap::new(),
                                 transcripts: vec![],
                             });
                ()},
            // TODO: mRNA isn't the only transcript type
            "mRNA" => {
                transcripts.insert(feat.uniquename.clone(),
                                   TranscriptDetails {
                                       uniquename: feat.uniquename.clone(),
                                       name: feat.name.clone(),
//                                       annotations: HashMap::new(),
                                   });
                ()},
            "genotype" => {
                genotypes.insert(feat.uniquename.clone(),
                                 GenotypeDetails {
                                     uniquename: feat.uniquename.clone(),
                                     name: feat.name.clone(),
                                     annotations: HashMap::new(),
                                 });
                ()},
            _ => (),
        }
    }

    for cvterm in &raw.cvterms {
        terms.insert(cvterm.termid(),
                     TermDetails {
                         name: cvterm.name.clone(),
                         termid: cvterm.termid(),
                         definition: cvterm.definition.clone(),
                         is_obsolete: cvterm.is_obsolete,
                         annotations: vec![],
                     });
    }

    let mut feature_parents: HashMap<String, String> = HashMap::new();

    for feature_rel in raw.feature_relationships.iter() {
        if feature_rel.subject.feat_type.name == "mRNA" &&
            feature_rel.rel_type.name == "part_of" &&
            feature_rel.object.feat_type.name == "gene" {
                feature_parents.insert(feature_rel.subject.uniquename.clone(),
                                       feature_rel.object.uniquename.clone());
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
        let feature_annotation =
            FeatureAnnotation {
                term: make_term_short(cvterm.clone()),
                evidence: evidence.clone(),
                publication: publication.clone(),
            };
        let cv_name = cvterm.cv.name.clone();
        let gene_details_opt =
            if feature.feat_type.name == "mRNA" {
                if let Some(gene_uniquename) = feature_parents.get(&feature.uniquename) {
                    genes.get_mut(gene_uniquename)
                } else {
                    None
                }
            } else {
                genes.get_mut(&feature.uniquename)
            };

        if let Some(gene_details) = gene_details_opt {
            gene_details.annotations.entry(cv_name).or_insert(Vec::new()).push(feature_annotation);
            let term_annotation =
                TermAnnotation {
                    gene: make_gene_short(&gene_details),
                    evidence: evidence,
                    publication: publication,
                };
            let termid = cvterm.termid();
            let term_details = terms.get_mut(&termid).unwrap();
            term_details.annotations.push(term_annotation);
        }
    }

    // remove terms with no annotation
    let used_terms = terms.clone().into_iter()
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
        gene_count: genes.len(),
        term_count: terms.len(),
    };

    WebData {
        genes: genes,
        terms: terms,
        used_terms: used_terms,
        metadata: metadata,
    }
}

fn write_genes(output_dir: &str, genes: &IdGeneMap) {
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

    write_genes(output_dir, &web_data.genes);
    write_terms(output_dir, &web_data.terms);
    write_metadata(output_dir, &web_data.metadata);

    println!("wrote {} genes", web_data.genes.len());
    println!("wrote {} terms", web_data.terms.len());
}

fn main() {
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
    let web_data = get_web_data(&raw, &organism_genus_species);

    write_web_data(&output_dir, &web_data);
}
