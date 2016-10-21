extern crate postgres;

use postgres::{Connection, TlsMode};

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
        pub features: Vec<Rc<Feature>>,
    }

    impl Raw {
        pub fn new(conn: &Connection) -> Raw {
            let mut ret = Raw {
                cvs: vec![], dbs: vec![], dbxrefs: vec![], cvterms: vec![],
                features: vec![],
            };

            let mut cv_map: HashMap<i32, Rc<Cv>> = HashMap::new();
            let mut db_map: HashMap<i32, Rc<Db>> = HashMap::new();
            let mut dbxref_map: HashMap<i32, Rc<Dbxref>> = HashMap::new();
            let mut cvterm_map: HashMap<i32, Rc<Cvterm>> = HashMap::new();
            let mut feature_map: HashMap<i32, Rc<Feature>> = HashMap::new();

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

            for row in &conn.query("SELECT feature_id, uniquename, name FROM feature", &[]).unwrap() {
                let feature_id = row.get(0);
                let feature = Feature {
                    uniquename: row.get(1),
                    name: row.get(2),
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
            pub name: Option<String>
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
    }
}

use pombase::Raw;

fn main() {
    let conn = Connection::connect("postgres://kmr44:kmr44@localhost/pombase-build-2016-09-20-v1", TlsMode::None).unwrap();

    let raw = Raw::new(&conn);

    println!("{}", raw.features.len());
}
