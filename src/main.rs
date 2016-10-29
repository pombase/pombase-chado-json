#![feature(proc_macro)]

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

mod pombase {
    use std::rc::Rc;
    use std::cell::RefCell;
    use std::collections::HashMap;
    use postgres::Connection;

    use pombase::db::*;

    pub struct Raw {
        pub organisms: Vec<Rc<Organism>>,
        pub cvs: Vec<Rc<Cv>>,
        pub dbs: Vec<Rc<Db>>,
        pub dbxrefs: Vec<Rc<Dbxref>>,
        pub cvterms: Vec<Rc<Cvterm>>,
        pub cvtermprops: Vec<Rc<Cvtermprop>>,
        pub cvtermpaths: Vec<Rc<Cvtermpath>>,
        pub cvterm_relationships: Vec<Rc<CvtermRelationship>>,
        pub publications: Vec<Rc<Publication>>,
        pub features: Vec<Rc<Feature>>,
        pub featureprops: Vec<Rc<Featureprop>>,
        pub feature_cvterms: Vec<Rc<FeatureCvterm>>,
        pub feature_cvtermprops: Vec<Rc<FeatureCvtermprop>>,
        pub feature_relationships: Vec<Rc<FeatureRelationship>>,
        pub chadoprops: Vec<Rc<Chadoprop>>,
    }

    impl Raw {

        pub fn new(conn: &Connection) -> Raw {
            let mut ret = Raw {
                organisms: vec![],
                cvs: vec![], dbs: vec![], dbxrefs: vec![], cvterms: vec![],
                cvtermprops: vec![],
                cvtermpaths: vec![], cvterm_relationships: vec![],
                publications: vec![], features: vec![], featureprops: vec![],
                feature_cvterms: vec![], feature_cvtermprops: vec![],
                feature_relationships: vec![], chadoprops: vec![],
            };

            let mut organism_map: HashMap<i32, Rc<Organism>> = HashMap::new();
            let mut cv_map: HashMap<i32, Rc<Cv>> = HashMap::new();
            let mut db_map: HashMap<i32, Rc<Db>> = HashMap::new();
            let mut dbxref_map: HashMap<i32, Rc<Dbxref>> = HashMap::new();
            let mut cvterm_map: HashMap<i32, Rc<Cvterm>> = HashMap::new();
            let mut feature_map: HashMap<i32, Rc<Feature>> = HashMap::new();
            let mut feature_cvterm_map: HashMap<i32, Rc<FeatureCvterm>> = HashMap::new();
            let mut publication_map: HashMap<i32, Rc<Publication>> = HashMap::new();

            fn get_db(db_map: &mut HashMap<i32, Rc<Db>>, db_id: i32) -> Rc<Db> {
                db_map.get(&db_id).unwrap().clone()
            };

            fn get_dbxref(dbxref_map: &mut HashMap<i32, Rc<Dbxref>>, dbxref_id: i32) -> Rc<Dbxref> {
                dbxref_map.get(&dbxref_id).unwrap().clone()
            };

            fn get_cv(cv_map: &mut HashMap<i32, Rc<Cv>>, cv_id: i32) -> Rc<Cv> {
                cv_map.get(&cv_id).unwrap().clone()
            };

            fn get_cvterm(cvterm_map: &mut HashMap<i32, Rc<Cvterm>>, cvterm_id: i32) -> Rc<Cvterm> {
                cvterm_map.get(&cvterm_id).unwrap().clone()
            };

            fn get_feature(feature_map: &mut HashMap<i32, Rc<Feature>>, feature_id: i32) -> Rc<Feature> {
                feature_map.get(&feature_id).unwrap().clone()
            };

            for row in &conn.query("SELECT organism_id, genus, species, abbreviation, common_name FROM organism", &[]).unwrap() {
                let organism = Organism {
                    genus: row.get(1),
                    species: row.get(2),
                    abbreviation: row.get(3),
                    common_name: row.get(4),
                };
                let rc_organism = Rc::new(organism);
                ret.organisms.push(rc_organism.clone());
                organism_map.insert(row.get(0), rc_organism);
            }

            for row in &conn.query("SELECT cv_id, name FROM cv", &[]).unwrap() {
                let cv = Cv {
                    name: row.get(1),
                };
                let rc_cv = Rc::new(cv);
                ret.cvs.push(rc_cv.clone());
                cv_map.insert(row.get(0), rc_cv);
            }

            for row in &conn.query("SELECT db_id, name FROM db", &[]).unwrap() {
                let db_id = row.get(0);
                let db = Db {
                    name: row.get(1),
                };
                let rc_db = Rc::new(db);
                ret.dbs.push(rc_db.clone());
                db_map.insert(db_id, rc_db);
            }

            for row in &conn.query("SELECT dbxref_id, db_id, accession FROM dbxref", &[]).unwrap() {
                let dbxref_id: i32 = row.get(0);
                let db_id: i32 = row.get(1);
                let dbxref = Dbxref {
                    db: get_db(&mut db_map, db_id),
                    accession: row.get(2),
                };
                let rc_dbxref = Rc::new(dbxref);
                ret.dbxrefs.push(rc_dbxref.clone());
                dbxref_map.insert(dbxref_id, rc_dbxref);
            }

            for row in &conn.query("SELECT cvterm_id, cv_id, dbxref_id, name, definition, is_obsolete FROM cvterm", &[]).unwrap() {
                let cvterm_id: i32 = row.get(0);
                let cv_id: i32 = row.get(1);
                let dbxref_id: i32 = row.get(2);
                let name: String = row.get(3);
                let definition: Option<String> = row.get(4);
                let is_obsolete: i32 = row.get(5);
                let cvterm = Cvterm {
                    cv: get_cv(&mut cv_map, cv_id),
                    dbxref: get_dbxref(&mut dbxref_map, dbxref_id),
                    name: name,
                    is_obsolete: is_obsolete == 1,
                    definition: definition,
                    cvtermprops: RefCell::new(vec![]),
                };
                let rc_cvterm = Rc::new(cvterm);
                ret.cvterms.push(rc_cvterm.clone());
                cvterm_map.insert(cvterm_id, rc_cvterm);
            }
            for row in &conn.query("SELECT cvterm_id, type_id, value FROM cvtermprop", &[]).unwrap() {
                let cvterm_id: i32 = row.get(0);
                let cvterm = get_cvterm(&mut cvterm_map, cvterm_id);
                let type_id: i32 = row.get(1);
                let value: String = row.get(2);
                let cvtermprop = Cvtermprop {
                    cvterm: cvterm.clone(),
                    prop_type: get_cvterm(&mut cvterm_map, type_id),
                    value: value,
                };
                let rc_cvtermprop = Rc::new(cvtermprop);
                ret.cvtermprops.push(rc_cvtermprop.clone());
                cvterm.cvtermprops.borrow_mut().push(rc_cvtermprop.clone());
            }

            for row in &conn.query("SELECT pub_id, uniquename, type_id, title, miniref FROM pub", &[]).unwrap() {
                let pub_id: i32 = row.get(0);
                let uniquename: String = row.get(1);
                let type_id: i32 = row.get(2);
                let title: Option<String> = row.get(3);
                let miniref: Option<String> = row.get(4);
                let publication = Publication {
                    uniquename: uniquename,
                    pub_type: get_cvterm(&mut cvterm_map, type_id),
                    title: title,
                    miniref: miniref,
                    publicationprops: RefCell::new(vec![]),
                };
                let rc_publication = Rc::new(publication);
                ret.publications.push(rc_publication.clone());
                publication_map.insert(pub_id, rc_publication);
            }

            for row in &conn.query("SELECT feature_id, uniquename, name, type_id, organism_id FROM feature", &[]).unwrap() {
                let feature_id = row.get(0);
                let type_id: i32 = row.get(3);
                let organism_id: i32 = row.get(4);
                let feature = Feature {
                    uniquename: row.get(1),
                    name: row.get(2),
                    feat_type: get_cvterm(&mut cvterm_map, type_id),
                    organism: organism_map.get(&organism_id).unwrap().clone(),
                    featureprops: RefCell::new(vec![]),
                };
                let rc_feature = Rc::new(feature);
                ret.features.push(rc_feature.clone());
                feature_map.insert(feature_id, rc_feature);
            }

            for row in &conn.query("SELECT feature_id, type_id, value FROM featureprop", &[]).unwrap() {
                let feature_id: i32 = row.get(0);
                let feature = get_feature(&mut feature_map, feature_id);
                let type_id: i32 = row.get(1);
                let value: Option<String> = row.get(2);
                let featureprop = Featureprop {
                    feature: feature.clone(),
                    prop_type: get_cvterm(&mut cvterm_map, type_id),
                    value: value,
                };
                let rc_featureprop = Rc::new(featureprop);
                ret.featureprops.push(rc_featureprop.clone());
                feature.featureprops.borrow_mut().push(rc_featureprop.clone());
            }

            for row in &conn.query("SELECT feature_cvterm_id, feature_id, cvterm_id, pub_id, is_not FROM feature_cvterm", &[]).unwrap() {
                let feature_cvterm_id = row.get(0);
                let feature_id = row.get(1);
                let cvterm_id: i32 = row.get(2);
                let pub_id: i32 = row.get(3);
                let is_not: bool = row.get(4);
                let feature_cvterm = FeatureCvterm {
                    feature_cvterm_id: feature_cvterm_id,
                    feature: feature_map.get(&feature_id).unwrap().clone(),
                    cvterm: cvterm_map.get(&cvterm_id).unwrap().clone(),
                    publication: publication_map.get(&pub_id).unwrap().clone(),
                    is_not: is_not,
                    feature_cvtermprops: RefCell::new(vec![]),
                };
                let rc_feature_cvterm = Rc::new(feature_cvterm);
                ret.feature_cvterms.push(rc_feature_cvterm.clone());
                feature_cvterm_map.insert(feature_cvterm_id, rc_feature_cvterm);
            }

            for row in &conn.query("SELECT feature_cvterm_id, type_id, value FROM feature_cvtermprop", &[]).unwrap() {
                let feature_cvterm_id: i32 = row.get(0);
                let feature_cvterm = feature_cvterm_map.get(&feature_cvterm_id).unwrap().clone();
                let type_id: i32 = row.get(1);
                let value: Option<String> = row.get(2);
                let feature_cvtermprop = FeatureCvtermprop {
                    feature_cvterm: feature_cvterm.clone(),
                    prop_type: get_cvterm(&mut cvterm_map, type_id),
                    value: value,
                };
                let rc_feature_cvtermprop = Rc::new(feature_cvtermprop);
                ret.feature_cvtermprops.push(rc_feature_cvtermprop.clone());
                feature_cvterm.feature_cvtermprops.borrow_mut().push(rc_feature_cvtermprop);
            }

            for row in &conn.query("SELECT subject_id, object_id, type_id FROM feature_relationship", &[]).unwrap() {
                let subject_id = row.get(0);
                let object_id: i32 = row.get(1);
                let type_id: i32 = row.get(2);
                let feature_relationship = FeatureRelationship {
                    subject: feature_map.get(&subject_id).unwrap().clone(),
                    object: feature_map.get(&object_id).unwrap().clone(),
                    rel_type: get_cvterm(&mut cvterm_map, type_id),
                    feature_relationshipprops: RefCell::new(vec![]),
                };
                let rc_feature_relationship = Rc::new(feature_relationship);
                ret.feature_relationships.push(rc_feature_relationship.clone());
            }

            for row in &conn.query("SELECT object_id, subject_id, type_id FROM cvterm_relationship", &[]).unwrap() {
                let subject_id = row.get(0);
                let object_id: i32 = row.get(1);
                let type_id: i32 = row.get(2);
                let cvterm_relationship = CvtermRelationship {
                    subject: cvterm_map.get(&subject_id).unwrap().clone(),
                    object: cvterm_map.get(&object_id).unwrap().clone(),
                    rel_type: get_cvterm(&mut cvterm_map, type_id),
                };
                let rc_cvterm_relationship = Rc::new(cvterm_relationship);
                ret.cvterm_relationships.push(rc_cvterm_relationship.clone());
            }

            for row in &conn.query("SELECT object_id, subject_id, type_id, pathdistance FROM cvtermpath", &[]).unwrap() {
                let subject_id = row.get(0);
                let object_id: i32 = row.get(1);
                let type_id: Option<i32> = row.get(2);
                let rel_type: Option<Rc<Cvterm>> = match type_id {
                    Some(cvterm_id) => Some(get_cvterm(&mut cvterm_map, cvterm_id)),
                    None => None
                };
                let pathdistance: Option<i32> = row.get(3);
                let cvtermpath = Cvtermpath {
                    subject: cvterm_map.get(&subject_id).unwrap().clone(),
                    object: cvterm_map.get(&object_id).unwrap().clone(),
                    rel_type: rel_type,
                    pathdistance: pathdistance,
                };
                let rc_cvtermpath = Rc::new(cvtermpath);
                ret.cvtermpaths.push(rc_cvtermpath.clone());
            }

            for row in &conn.query("SELECT type_id, value FROM chadoprop", &[]).unwrap() {
                let type_id: i32 = row.get(0);
                let value = row.get(1);
                let chadoprop = Chadoprop {
                    prop_type: cvterm_map.get(&type_id).unwrap().clone(),
                    value: value,
                };
                let rc_chadoprop = Rc::new(chadoprop);
                ret.chadoprops.push(rc_chadoprop.clone());
            }

            ret
        }
    }

    pub mod db {
        use std::rc::Rc;
        use std::cell::RefCell;

        pub trait Prop {
            fn type_name(&self) -> String;
            fn value(&self) -> Option<String>;
        }

        pub struct Chadoprop {
            pub prop_type: Rc<Cvterm>,
            pub value: Option<String>,
        }

        pub struct Organism {
            pub genus: String,
            pub species: String,
            pub abbreviation: String,
            pub common_name: String,
        }
        pub struct Cv {
            pub name: String,
        }
        pub struct Db {
            pub name: String,
        }
        pub struct Dbxref {
            pub accession: String,
            pub db: Rc<Db>,
        }
        pub struct Cvterm {
            pub name: String,
            pub cv: Rc<Cv>,
            pub dbxref: Rc<Dbxref>,
            pub definition: Option<String>,
            pub is_obsolete: bool,
            pub cvtermprops: RefCell<Vec<Rc<Cvtermprop>>>,
        }
        impl Cvterm {
            pub fn termid(&self) -> String {
                String::new() + &self.dbxref.db.name + ":" + &self.dbxref.accession
            }
        }
        pub struct Cvtermprop {
            pub cvterm: Rc<Cvterm>,
            pub prop_type: Rc<Cvterm>,
            pub value: String,
        }
        pub struct Publication {
            pub uniquename: String,
            pub pub_type: Rc<Cvterm>,
            pub title: Option<String>,
            pub miniref: Option<String>,
            pub publicationprops: RefCell<Vec<Rc<Publicationprop>>>,
        }
        pub struct Publicationprop {
            pub publication: Rc<Publication>,
            pub prop_type: Rc<Cvterm>,
            pub value: Option<String>,
        }
        pub struct CvtermRelationship {
            pub subject: Rc<Cvterm>,
            pub object: Rc<Cvterm>,
            pub rel_type: Rc<Cvterm>,
        }
        pub struct Cvtermpath {
            pub subject: Rc<Cvterm>,
            pub object: Rc<Cvterm>,
            pub rel_type: Option<Rc<Cvterm>>,
            pub pathdistance: Option<i32>,
        }
        pub struct Feature {
            pub uniquename: String,
            pub name: Option<String>,
            pub feat_type: Rc<Cvterm>,
            pub organism: Rc<Organism>,
            pub featureprops: RefCell<Vec<Rc<Featureprop>>>,
        }
        pub struct Featureprop {
            pub feature: Rc<Feature>,
            pub prop_type: Rc<Cvterm>,
            pub value: Option<String>,
        }
        pub struct FeatureCvterm {
            pub feature_cvterm_id: i32,
            pub feature: Rc<Feature>,
            pub cvterm: Rc<Cvterm>,
            pub publication: Rc<Publication>,
            pub is_not: bool,
            pub feature_cvtermprops: RefCell<Vec<Rc<FeatureCvtermprop>>>,
        }
        pub struct FeatureCvtermprop {
            pub feature_cvterm: Rc<FeatureCvterm>,
            pub prop_type: Rc<Cvterm>,
            pub value: Option<String>,
        }
        impl Prop for FeatureCvtermprop {
            fn type_name(&self) -> String {
                self.prop_type.name.clone()
            }
            fn value(&self) -> Option<String> {
                self.value.clone()
            }
        }

        pub struct FeatureRelationship {
            pub subject: Rc<Feature>,
            pub object: Rc<Feature>,
            pub rel_type: Rc<Cvterm>,
            pub feature_relationshipprops: RefCell<Vec<Rc<FeatureRelationshipprop>>>,
        }
        pub struct FeatureRelationshipprop {
            pub feature_relationship: Rc<FeatureRelationship>,
            pub prop_type: Rc<Cvterm>,
            pub value: Option<String>,
        }
    }

    pub mod web {
        pub mod data {
            use std::collections::HashMap;

            pub type GeneUniquename = String;
            pub type GeneName = String;
            pub type TypeName = String;

            #[derive(Serialize, Clone)]
            pub struct GeneShort {
                pub uniquename: GeneUniquename,
                pub name: Option<GeneName>,
            }

            pub type TranscriptUniquename = String;

            #[derive(Serialize)]
            pub struct TranscriptShort {
                pub uniquename: TranscriptUniquename,
//                pub exons: Vec<ExonShort>,
//                pub utrs: Vec<UTRShort>,
            }

            pub type TermName = String;
            pub type TermId = String;
            pub type TermDef = String;

            #[derive(Serialize)]
            pub struct TermShort {
                pub name: TermName,
                pub termid: TermId,
                pub is_obsolete: bool,
            }

            #[derive(Serialize, Clone)]
            pub struct PublicationShort {
                pub uniquename: String,
                pub title: Option<String>,
                pub citation: Option<String>,
            }

            type Evidence = String;

            #[derive(Serialize)]
            pub struct FeatureAnnotation {
                pub term: TermShort,
                pub evidence: Option<Evidence>,
                pub publication: Option<PublicationShort>,
            }

            pub type TypeFeatureAnnotationMap =
                HashMap<TypeName, Vec<FeatureAnnotation>>;

            #[derive(Serialize)]
            pub struct GeneDetails {
                pub uniquename: GeneUniquename,
                pub name: Option<String>,
                pub transcripts: Vec<TranscriptShort>,
                pub annotations: TypeFeatureAnnotationMap,
            }

            pub type UniquenameGeneMap =
                HashMap<GeneUniquename, GeneDetails>;

            #[derive(Serialize)]
            pub struct TranscriptDetails {
                pub uniquename: TranscriptUniquename,
                pub name: Option<String>,
                pub annotations: TypeFeatureAnnotationMap,
            }

            pub type UniquenameTranscriptMap =
                HashMap<TranscriptUniquename, TranscriptDetails>;

            #[derive(Serialize, Clone)]
            pub struct TermAnnotation {
                pub gene: GeneShort,
                pub evidence: Option<Evidence>,
                pub publication: Option<PublicationShort>,
            }

            #[derive(Serialize)]
            pub struct TermDetails {
                pub name: TermName,
                pub termid: TermId,
                pub definition: Option<TermDef>,
                pub is_obsolete: bool,
                pub annotations: Vec<TermAnnotation>,
            }

            impl Clone for TermDetails {
                fn clone(&self) -> TermDetails {
                    TermDetails {
                        name: self.name.clone(),
                        termid: self.termid.clone(),
                        definition: self.definition.clone(),
                        is_obsolete: self.is_obsolete,
                        annotations: self.annotations.clone(),
                    }
                }
            }

            pub type IdTermMap =
                HashMap<TermId, TermDetails>;

            #[derive(Serialize)]
            pub struct Metadata {
                pub db_creation_datetime: String,
                pub gene_count: usize,
                pub term_count: usize,
            }
        }
    }
}

use pombase::Raw;
use pombase::web::data::*;
use pombase::db::Prop;

use std::rc::Rc;

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
    genes: HashMap<GeneUniquename, GeneDetails>,
    terms: HashMap<TermId, TermDetails>,
    used_terms: HashMap<TermId, TermDetails>,
    metadata: Metadata,
}

fn get_web_data(raw: &Raw, organism_genus_species: &String) -> WebData {
    let mut genes: UniquenameGeneMap = HashMap::new();
    let mut transcripts: UniquenameTranscriptMap = HashMap::new();
    let mut terms: IdTermMap = HashMap::new();

//    type IdAnnotationMap = HashMap<i32, (FeatureAnnotation, TermAnnotation)>;

//    let annotation_map: IdAnnotationMap = HashMap::new();

    for feat in raw.features.iter().filter(|&f| {
        let feature_org_genus_species = String::new() +
            &f.organism.genus + "_" + &f.organism.species;
        feature_org_genus_species == *organism_genus_species
    } ) {

        if feat.feat_type.name == "gene" {
            genes.insert(feat.uniquename.clone(),
                         GeneDetails {
                             uniquename: feat.uniquename.clone(),
                             name: feat.name.clone(),
                             annotations: HashMap::new(),
                             transcripts: vec![],
                         });
        } else {
            // TODO: mRNA isn't the only transcript type
            if feat.feat_type.name == "mRNA" {
                transcripts.insert(feat.uniquename.clone(),
                             TranscriptDetails {
                                 uniquename: feat.uniquename.clone(),
                                 name: feat.name.clone(),
                                 annotations: HashMap::new(),
                             });
            }
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

    let WebData {
        genes,
        terms,
        used_terms: _used_terms,
        metadata
    } = web_data;

    let s = serde_json::to_string(&genes).unwrap();
    let file_name = String::new() + &output_dir + "/genes.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    for (gene_uniquename, gene_details) in &genes {
        let s = serde_json::to_string(&gene_details).unwrap();
        let file_name = String::new() + &output_dir + "/gene/" + &gene_uniquename + ".json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }

    let s = serde_json::to_string(&terms).unwrap();
    let file_name = String::new() + &output_dir + "/terms.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    for (termid, term_details) in &terms {
        let s = serde_json::to_string(&term_details).unwrap();
        let file_name = String::new() + &output_dir + "/term/" + &termid + ".json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }

    let s = serde_json::to_string(&metadata).unwrap();
    let file_name = String::new() + &output_dir + "/metadata.json";
    let f = File::create(file_name).expect("Unable to open file");
    let mut writer = BufWriter::new(&f);
    writer.write_all(s.as_bytes()).expect("Unable to write!");

    println!("wrote {} genes", genes.len());
    println!("wrote {} terms", terms.len());
}
