extern crate postgres;

use std::rc::Rc;
use std::cell::RefCell;
use std::collections::HashMap;
use self::postgres::Connection;

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
    pub publicationprops: Vec<Rc<Publicationprop>>,
    pub synonyms: Vec<Rc<Synonym>>,
    pub features: Vec<Rc<Feature>>,
    pub feature_synonyms: Vec<Rc<FeatureSynonym>>,
    pub featureprops: Vec<Rc<Featureprop>>,
    pub featurelocs: Vec<Rc<Featureloc>>,
    pub feature_cvterms: Vec<Rc<FeatureCvterm>>,
    pub feature_cvtermprops: Vec<Rc<FeatureCvtermprop>>,
    pub feature_relationships: Vec<Rc<FeatureRelationship>>,
    pub chadoprops: Vec<Rc<Chadoprop>>,
}

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
pub type CvName = String;
pub struct Cv {
    pub name: CvName,
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
    pub value: String,
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
    pub featurelocs: RefCell<Vec<Rc<Featureloc>>>,
}
pub struct Featureloc {
    pub feature: Rc<Feature>,
    pub srcfeature: Rc<Feature>,
    pub fmin: i32,
    pub fmax: i32,
    pub strand: i16,
//    pub phase: i32,
}
pub struct Featureprop {
    pub feature: Rc<Feature>,
    pub prop_type: Rc<Cvterm>,
    pub value: Option<String>,
}
pub struct Synonym {
    pub name: String,
    pub synonym_type: Rc<Cvterm>
}
pub struct FeatureSynonym {
    pub feature: Rc<Feature>,
    pub synonym: Rc<Synonym>,
    pub publication: Rc<Publication>,
    pub is_current: bool,
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
    pub publications: RefCell<Vec<Rc<Publication>>>,
}
pub struct FeatureRelationshipprop {
    pub feature_relationship: Rc<FeatureRelationship>,
    pub prop_type: Rc<Cvterm>,
    pub value: Option<String>,
}

impl Raw {
    pub fn new(conn: &Connection) -> Raw {
        let mut ret = Raw {
            organisms: vec![],
            cvs: vec![], dbs: vec![], dbxrefs: vec![], cvterms: vec![],
            synonyms: vec![], cvtermprops: vec![],
            cvtermpaths: vec![], cvterm_relationships: vec![],
            publications: vec![], publicationprops: vec![],
            features: vec![], featureprops: vec![],
            featurelocs: vec![], feature_synonyms: vec![],
            feature_cvterms: vec![], feature_cvtermprops: vec![],
            feature_relationships: vec![], chadoprops: vec![],
        };

        let mut organism_map: HashMap<i32, Rc<Organism>> = HashMap::new();
        let mut cv_map: HashMap<i32, Rc<Cv>> = HashMap::new();
        let mut db_map: HashMap<i32, Rc<Db>> = HashMap::new();
        let mut dbxref_map: HashMap<i32, Rc<Dbxref>> = HashMap::new();
        let mut cvterm_map: HashMap<i32, Rc<Cvterm>> = HashMap::new();
        let mut synonym_map: HashMap<i32, Rc<Synonym>> = HashMap::new();
        let mut feature_map: HashMap<i32, Rc<Feature>> = HashMap::new();
        let mut feature_cvterm_map: HashMap<i32, Rc<FeatureCvterm>> = HashMap::new();
        let mut feature_relationship_map: HashMap<i32, Rc<FeatureRelationship>> = HashMap::new();
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

        fn get_synonym(synonym_map: &mut HashMap<i32, Rc<Synonym>>, synonym_id: i32) -> Rc<Synonym> {
            synonym_map.get(&synonym_id).unwrap().clone()
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

        for row in &conn.query("SELECT pub_id, type_id, value FROM pubprop", &[]).unwrap() {
            let pub_id: i32 = row.get(0);
            let type_id: i32 = row.get(1);
            let value: String = row.get(2);
            let publication = publication_map.get(&pub_id).unwrap().clone();
            let publicationprop = Publicationprop {
                publication: publication.clone(),
                prop_type: get_cvterm(&mut cvterm_map, type_id),
                value: value,
            };
            let rc_publicationprop = Rc::new(publicationprop);
            ret.publicationprops.push(rc_publicationprop.clone());
            publication.publicationprops.borrow_mut().push(rc_publicationprop.clone());
        }

        for row in &conn.query("SELECT synonym_id, name, type_id FROM synonym", &[]).unwrap() {
            let synonym_id = row.get(0);
            let name = row.get(1);
            let type_id: i32 = row.get(2);
            let synonym = Synonym {
                name: name,
                synonym_type: get_cvterm(&mut cvterm_map, type_id),
            };
            let rc_synonym = Rc::new(synonym);
            ret.synonyms.push(rc_synonym.clone());
            synonym_map.insert(synonym_id, rc_synonym);
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
                featurelocs: RefCell::new(vec![]),
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

        for row in &conn.query("SELECT feature_id, srcfeature_id, fmin, fmax, strand FROM featureloc", &[]).unwrap() {
            let feature_id: i32 = row.get(0);
            let feature = get_feature(&mut feature_map, feature_id);
            let srcfeature_id: i32 = row.get(1);
            let srcfeature = get_feature(&mut feature_map, srcfeature_id);
            let fmin: i32 = row.get(2);
            let fmax: i32 = row.get(3);
            let strand: i16 = row.get(4);
//            let phase: i32 = row.get(5);
            let featureloc = Featureloc {
                feature: feature.clone(),
                srcfeature: srcfeature.clone(),
                fmin: fmin,
                fmax: fmax,
                strand: strand,
//                phase: phase,
            };
            let rc_featureloc = Rc::new(featureloc);
            ret.featurelocs.push(rc_featureloc.clone());
            feature.featurelocs.borrow_mut().push(rc_featureloc.clone());
        }

        for row in &conn.query("SELECT feature_id, synonym_id, pub_id, is_current FROM feature_synonym", &[]).unwrap() {
            let feature_id: i32 = row.get(0);
            let synonym_id: i32 = row.get(1);
            let pub_id: i32 = row.get(2);
            let is_current: bool = row.get(3);
            let feature_synonym = FeatureSynonym {
                feature: feature_map.get(&feature_id).unwrap().clone(),
                synonym: get_synonym(&mut synonym_map, synonym_id),
                publication: publication_map.get(&pub_id).unwrap().clone(),
                is_current: is_current
            };
            let rc_feature_synonym = Rc::new(feature_synonym);
            ret.feature_synonyms.push(rc_feature_synonym.clone());
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

        for row in &conn.query("SELECT feature_relationship_id, subject_id, object_id, type_id FROM feature_relationship", &[]).unwrap() {
            let feature_relationship_id = row.get(0);
            let subject_id = row.get(1);
            let object_id: i32 = row.get(2);
            let type_id: i32 = row.get(3);
            let feature_relationship = FeatureRelationship {
                subject: feature_map.get(&subject_id).unwrap().clone(),
                object: feature_map.get(&object_id).unwrap().clone(),
                rel_type: get_cvterm(&mut cvterm_map, type_id),
                feature_relationshipprops: RefCell::new(vec![]),
                publications: RefCell::new(vec![]),
            };
            let rc_feature_relationship = Rc::new(feature_relationship);
            ret.feature_relationships.push(rc_feature_relationship.clone());
            feature_relationship_map.insert(feature_relationship_id, rc_feature_relationship);
        }

        for row in &conn.query("SELECT feature_relationship_id, type_id, value FROM feature_relationshipprop", &[]).unwrap() {
            let feature_relationship_id: i32 = row.get(0);
            let feature_relationship =
                feature_relationship_map.get(&feature_relationship_id).unwrap().clone();
            let type_id: i32 = row.get(1);
            let value: Option<String> = row.get(2);
            let feature_relationshipprop = FeatureRelationshipprop {
                feature_relationship: feature_relationship.clone(),
                prop_type: get_cvterm(&mut cvterm_map, type_id),
                value: value,
            };
            let rc_feature_relationshipprop = Rc::new(feature_relationshipprop);
            feature_relationship.feature_relationshipprops.borrow_mut().push(rc_feature_relationshipprop.clone());
        }

        for row in &conn.query("SELECT feature_relationship_id, pub_id FROM feature_relationship_pub", &[]).unwrap() {
            let feature_relationship_id: i32 = row.get(0);
            let feature_relationship =
                feature_relationship_map.get(&feature_relationship_id).unwrap().clone();
            let pub_id: i32 = row.get(1);
            let publication = publication_map.get(&pub_id).unwrap().clone();
            feature_relationship.publications.borrow_mut().push(publication.clone());
        }

        for row in &conn.query("SELECT object_id, subject_id, type_id FROM cvterm_relationship", &[]).unwrap() {
            let object_id = row.get(0);
            let subject_id: i32 = row.get(1);
            let type_id: i32 = row.get(2);
            let cvterm_relationship = CvtermRelationship {
                subject: cvterm_map.get(&subject_id).unwrap().clone(),
                object: cvterm_map.get(&object_id).unwrap().clone(),
                rel_type: get_cvterm(&mut cvterm_map, type_id),
            };
            let rc_cvterm_relationship = Rc::new(cvterm_relationship);
            ret.cvterm_relationships.push(rc_cvterm_relationship.clone());
        }

        for row in &conn.query("SELECT subject_id, object_id, type_id, pathdistance FROM cvtermpath", &[]).unwrap() {
            let subject_id: i32 = row.get(0);
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
