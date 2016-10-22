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

mod pombase {
    use std::rc::Rc;
    use std::collections::HashMap;
    use postgres::Connection;

    use pombase::db::*;

    pub struct Raw {
        pub cvs: Vec<Rc<Cv>>,
        pub dbs: Vec<Rc<Db>>,
        pub dbxrefs: Vec<Rc<Dbxref>>,
        pub cvterms: Vec<Rc<Cvterm>>,
        pub cvtermpaths: Vec<Rc<Cvtermpath>>,
        pub cvterm_relationships: Vec<Rc<CvtermRelationship>>,
        pub publications: Vec<Rc<Publication>>,
        pub features: Vec<Rc<Feature>>,
    }

    impl Raw {
        pub fn new(conn: &Connection) -> Raw {
            let mut ret = Raw {
                cvs: vec![], dbs: vec![], dbxrefs: vec![], cvterms: vec![],
                cvtermpaths: vec![], cvterm_relationships: vec![],
                publications: vec![], features: vec![],
            };

            let mut cv_map: HashMap<i32, Rc<Cv>> = HashMap::new();
            let mut db_map: HashMap<i32, Rc<Db>> = HashMap::new();
            let mut dbxref_map: HashMap<i32, Rc<Dbxref>> = HashMap::new();
            let mut cvterm_map: HashMap<i32, Rc<Cvterm>> = HashMap::new();
            let mut feature_map: HashMap<i32, Rc<Feature>> = HashMap::new();
            let mut publication_map: HashMap<i32, Rc<Publication>> = HashMap::new();

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
                    db: db_map.get(&db_id).unwrap().clone(),
                    accession: row.get(2),
                };
                let rc_dbxref = Rc::new(dbxref);
                ret.dbxrefs.push(rc_dbxref.clone());
                dbxref_map.insert(dbxref_id, rc_dbxref);
            }

            for row in &conn.query("SELECT cvterm_id, cv_id, dbxref_id, name FROM cvterm", &[]).unwrap() {
                let cvterm_id: i32 = row.get(0);
                let cv_id: i32 = row.get(1);
                let dbxref_id: i32 = row.get(2);
                let cvterm = Cvterm {
                    cv: cv_map.get(&cv_id).unwrap().clone(),
                    dbxref: dbxref_map.get(&dbxref_id).unwrap().clone(),
                    name: row.get(3),
                };
                let rc_cvterm = Rc::new(cvterm);
                ret.cvterms.push(rc_cvterm.clone());
                cvterm_map.insert(cvterm_id, rc_cvterm);
            }

            for row in &conn.query("SELECT pub_id, uniquename, type_id, title FROM pub", &[]).unwrap() {
                let pub_id: i32 = row.get(0);
                let uniquename: String = row.get(1);
                let type_id: i32 = row.get(2);
                let title: Option<String> = row.get(3);
                let publication = Publication {
                    uniquename: uniquename,
                    pub_type: cvterm_map.get(&type_id).unwrap().clone(),
                    title: title
                };
                let rc_publication = Rc::new(publication);
                ret.publications.push(rc_publication.clone());
                publication_map.insert(pub_id, rc_publication);
            }

            for row in &conn.query("SELECT feature_id, uniquename, name, type_id FROM feature", &[]).unwrap() {
                let feature_id = row.get(0);
                let type_id: i32 = row.get(3);
                let feature = Feature {
                    uniquename: row.get(1),
                    name: row.get(2),
                    feat_type: cvterm_map.get(&type_id).unwrap().clone(),
                };
                let rc_feature = Rc::new(feature);
                ret.features.push(rc_feature.clone());
                feature_map.insert(feature_id, rc_feature);
            }

            ret
        }
    }

    pub mod db {
        use std::rc::Rc;

        pub struct Feature {
            pub uniquename: String,
            pub name: Option<String>,
            pub feat_type: Rc<Cvterm>,
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
        }
        pub struct Publication {
            pub uniquename: String,
            pub pub_type: Rc<Cvterm>,
            pub title: Option<String>,
        }
        pub struct CvtermRelationship {
            pub subject: Rc<Cvterm>,
            pub object: Rc<Cvterm>,
            pub reltype: Rc<Cvterm>,
        }
        pub struct Cvtermpath {
            pub subject: Rc<Cvterm>,
            pub object: Rc<Cvterm>,
            pub rel_type: Rc<Cvterm>,
            pub pathdistance: Option<i32>,
        }
    }

    pub mod web {
        pub mod data {
            use std::collections::HashMap;

            type GeneUniquename = String;
            type GeneName = String;

            #[derive(Serialize)]
            pub struct GeneShort {
                pub uniquename: GeneUniquename,
                pub name: Option<GeneName>,
            }

            type TermName = String;
            type TermId = String;
            type TermDef = String;

            #[derive(Serialize)]
            pub struct TermShort {
                pub name: TermName,
                pub termid: TermId,
                pub definition: TermDef,
                pub is_obsolete: bool,
            }

            type Evidence = String;

            #[derive(Serialize)]
            pub struct FeatureAnnotation {
                pub term: TermShort,
                pub evidence: Evidence,
            }

            type TypeName = String;
            pub type TypeFeatureAnnotationMap =
                HashMap<TypeName, Vec<FeatureAnnotation>>;

            #[derive(Serialize)]
            pub struct GeneDetails {
                pub uniquename: GeneUniquename,
                pub name: Option<String>,
                pub annotations: TypeFeatureAnnotationMap,
            }

            type UniquenameGeneMap =
                HashMap<GeneUniquename, Vec<GeneDetails>>;

            #[derive(Serialize)]
            pub struct TermAnnotation {
                pub gene: GeneShort,
                pub evidence: Evidence,
            }

            #[derive(Serialize)]
            pub struct TermDetails {
                pub name: TermName,
                pub termid: TermId,
                pub definition: TermDef,
                pub is_obsolete: bool,
                pub annotations: Vec<TermAnnotation>,
            }
        }
    }
}

use pombase::Raw;
use pombase::web::data::*;

fn get_web_data(raw: &Raw) -> (Vec<GeneDetails>, Vec<TermDetails>) {
    let mut genes: Vec<GeneDetails> = vec![];
    let terms: Vec<TermDetails> = vec![];

    for feat in raw.features.iter().filter(|&f| f.feat_type.name == "gene") {
        genes.push(GeneDetails {
            uniquename: feat.uniquename.clone(),
            name: feat.name.clone(),
            annotations: TypeFeatureAnnotationMap::new(),
        });
    }

    (genes, terms)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();

    opts.reqopt("c", "connection-string",
                "PostgresSQL connection string like: postgres://user:pass@host/db_name",
                "CONN_STR");
    opts.reqopt("d", "output-directory",
                "Destination directory for JSON output", "DIR");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("Invalid options\n{}", f)
    };

    let connection_string = matches.opt_str("c").unwrap();
    let output_dir = matches.opt_str("d").unwrap();

    fs::create_dir_all(&output_dir).unwrap_or_else(|why| {
        println!("Creating output directory failed: {:?}", why.kind());
    });

    let conn = Connection::connect(connection_string.as_str(), TlsMode::None).unwrap();

    let raw = Raw::new(&conn);

    let (genes, terms) = get_web_data(&raw);

    println!("{}", genes.len() + terms.len());
    println!("{}", genes.get(0).unwrap().uniquename);

    for gene in &genes {
        let s = serde_json::to_string(&gene).unwrap();
        let file_name = String::new() + &output_dir + "/" + &gene.uniquename + ".json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }

}
