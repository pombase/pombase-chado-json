use std::rc::Rc;
use std::cell::RefCell;
use std::collections::HashMap;

extern crate pombase;

use pombase::db::*;
use pombase::web::config::*;
use pombase::web::data::*;
use pombase::web::data_build::*;

fn make_test_cvterm_dbxref(cvterms: &mut Vec<Rc<Cvterm>>, dbxrefs: &mut Vec<Rc<Dbxref>>,
                           cv: &Rc<Cv>, db: &Rc<Db>, cvterm_name: &str, accession: &str)
                           -> Rc<Cvterm> {
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
        is_relationshiptype: false,
        cvtermprops: RefCell::new(vec![]),
    });
    cvterms.push(cvterm.clone());
    dbxrefs.push(dbxref.clone());
    cvterm
}

fn make_test_cv(cvs: &mut Vec<Rc<Cv>>, cv_name: &str) -> Rc<Cv> {
    let cv = Rc::new(Cv {
        name: String::from(cv_name),
    });
    cvs.push(cv.clone());
    cv
}

fn make_test_db(dbs: &mut Vec<Rc<Db>>, db_name: &str) -> Rc<Db> {
    let db = Rc::new(Db {
        name: String::from(db_name),
    });
    dbs.push(db.clone());
    db
}

fn make_test_feature(features: &mut Vec<Rc<Feature>>, organism: &Rc<Organism>,
                     type_cvterm: &Rc<Cvterm>, uniquename: &str, name: Option<String>)
                     -> Rc<Feature> {
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

fn make_test_featureprop(featureprops: &mut Vec<Rc<Featureprop>>, feature: &Rc<Feature>,
                         type_cvterm: &Rc<Cvterm>, value: Option<String>) -> Rc<Featureprop> {
    let featureprop = Rc::new(Featureprop {
        feature: feature.clone(),
        prop_type: type_cvterm.clone(),
        value: value,
    });
    feature.featureprops.borrow_mut().push(featureprop.clone());
    featureprops.push(featureprop.clone());
    featureprop
}

fn make_test_feature_cvterm(feature_cvterms: &mut Vec<Rc<FeatureCvterm>>,
                            feature: &Rc<Feature>, cvterm: &Rc<Cvterm>,
                            publication: &Rc<Publication>)
                     -> Rc<FeatureCvterm> {
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

fn make_test_feature_cvtermprop(feature_cvtermprops: &mut Vec<Rc<FeatureCvtermprop>>,
                                feature_cvterm: &Rc<FeatureCvterm>, prop_type: &Rc<Cvterm>,
                                value: &str) -> Rc<FeatureCvtermprop> {
    let feature_cvtermprop = Rc::new(FeatureCvtermprop {
        feature_cvterm: feature_cvterm.clone(),
        prop_type: prop_type.clone(),
        value: Some(value.into()),
    });
    feature_cvtermprops.push(feature_cvtermprop.clone());
    feature_cvterm.feature_cvtermprops.borrow_mut().push(feature_cvtermprop.clone());
    feature_cvtermprop
}

fn make_test_feature_rel(feature_relationships: &mut Vec<Rc<FeatureRelationship>>,
                         publication: &Rc<Publication>,
                         subject: &Rc<Feature>, rel: &Rc<Cvterm>, object: &Rc<Feature>) {
    let rel = Rc::new(FeatureRelationship {
        subject: subject.clone(),
        rel_type: rel.clone(),
        object: object.clone(),
        feature_relationshipprops: RefCell::new(vec![]),
        publications: RefCell::new(vec![publication.clone()]),
    });
    feature_relationships.push(rel);
}

fn make_test_cvterm_rel(cvterm_relationships: &mut Vec<Rc<CvtermRelationship>>,
                        subject: &Rc<Cvterm>, rel_type: &Rc<Cvterm>, object: &Rc<Cvterm>) {
    let rel = Rc::new(CvtermRelationship {
        subject: subject.clone(),
        rel_type: rel_type.clone(),
        object: object.clone(),
    });
    cvterm_relationships.push(rel);
}

fn make_test_cvtermpath(cvtermpaths: &mut Vec<Rc<Cvtermpath>>,
                        subject: &Rc<Cvterm>, rel_type: &Rc<Cvterm>, object: &Rc<Cvterm>,
                        pathdistance: i32) {
    let rel = Rc::new(Cvtermpath {
        subject: subject.clone(),
        rel_type: Some(rel_type.clone()),
        object: object.clone(),
        pathdistance: Some(pathdistance),
    });
    cvtermpaths.push(rel);
}

fn get_test_raw() -> Raw {
    let mut feature_relationships: Vec<Rc<FeatureRelationship>> = vec![];
    let mut cvterm_relationships: Vec<Rc<CvtermRelationship>> = vec![];
    let mut cvtermpaths: Vec<Rc<Cvtermpath>> = vec![];
    let mut cvs: Vec<Rc<Cv>> = vec![];
    let mut dbs: Vec<Rc<Db>> = vec![];
    let mut cvterms: Vec<Rc<Cvterm>> = vec![];
    let mut dbxrefs: Vec<Rc<Dbxref>> = vec![];
    let mut features: Vec<Rc<Feature>> = vec![];
    let mut featureprops: Vec<Rc<Featureprop>> = vec![];
    let mut feature_cvterms: Vec<Rc<FeatureCvterm>> = vec![];
    let mut feature_cvtermprops: Vec<Rc<FeatureCvtermprop>> = vec![];

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
    let featureprop_types_cv = make_test_cv(&mut cvs, "PomBase feature property types");
    let pub_type_cv = make_test_cv(&mut cvs, "PomBase publication types");
    let pubprop_type_cv = make_test_cv(&mut cvs, "pubprop_type");
    let feature_cvtermprop_type_cv = make_test_cv(&mut cvs, "feature_cvtermprop_type");

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
    let allele_type_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &featureprop_types_cv, &pbo_db,
                                "allele_type", "0011954");
    let description_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &featureprop_types_cv, &pbo_db,
                                "description", "0011059");
    let polypeptide_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &sequence_cv, &so_db,
                                "polypeptide", "0000104");
    let chromosome_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &sequence_cv, &so_db,
                                "chromosome", "0000340");
    let paper_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &pub_type_cv, &pbo_db,
                                "paper", "0000044");
    let pubmed_publication_date_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &pubprop_type_cv, &pbo_db,
                                "pubmed_publication_date", "0034034");
    let pubmed_authors_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &pubprop_type_cv, &pbo_db,
                                "pubmed_authors", "0034035");
    let with_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &feature_cvtermprop_type_cv, &pbo_db,
                                "with", "0000098");

    let publication = Rc::new(Publication {
        uniquename: String::from("PMID:11707284"),
        title: Some(String::from("The protein phosphatase 2A B'-regulatory subunit par1p is implicated in regulation of the S. pombe septation initiation network.")),
        pub_type: paper_cvterm,
        miniref: Some(String::from("FEBS Lett. 2001 Nov 9;508(1):136-42")),
        publicationprops: RefCell::new(vec![]),
    });

    let publication_pub_date = Rc::new(Publicationprop {
        publication: publication.clone(),
        prop_type: pubmed_publication_date_cvterm,
        value: String::from("9 Nov 2001"),
    });
    let publication_authors = Rc::new(Publicationprop {
        publication: publication.clone(),
        prop_type: pubmed_authors_cvterm,
        value: String::from("Le Goff X, Buvelot S, Salimova E, Guerry F, Schmidt S, Cueille N, Cano E, Simanis V"),
    });

    publication.publicationprops.borrow_mut().push(publication_pub_date);
    publication.publicationprops.borrow_mut().push(publication_authors);

    let bp_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &bp_cv, &go_db, "biological_process",
                                "0008150");
    let go0031030_cvterm =
        make_test_cvterm_dbxref(&mut cvterms, &mut dbxrefs, &bp_cv, &go_db, "negative regulation of septation initiation signaling",
                                "0031030");

    make_test_cvtermpath(&mut cvtermpaths, &go0031030_cvterm, &is_a_cvterm, &bp_cvterm,
                         11);
    make_test_cvtermpath(&mut cvtermpaths, &bp_cvterm, &is_a_cvterm, &go0031030_cvterm,
                         -11);

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

    let pom1_gene = make_test_feature(&mut features, &pombe_organism, &gene_cvterm,
                                      "SPAC2F7.03c", Some(String::from("pom1")));

    let cdc16_gene = make_test_feature(&mut features, &pombe_organism, &gene_cvterm,
                                      "SPAC6F6.08c", Some(String::from("cdc16")));
    let cdc16_allele1 = make_test_feature(&mut features, &pombe_organism, &allele_cvterm,
                                          "SPAC6F6.08c-allele1", Some(String::from("cdc16::ura4+")));
    make_test_featureprop(&mut featureprops, &cdc16_allele1, &allele_type_cvterm,
                          Some(String::from("disruption")));
    make_test_feature_rel(&mut feature_relationships, &publication,
                          &cdc16_allele1, &instance_of_cvterm, &cdc16_gene);

    let cdc16_delta_allele = make_test_feature(&mut features, &pombe_organism, &allele_cvterm,
                                               "SPAC6F6.08c-allele2", Some(String::from("cdc16delta")));
    make_test_featureprop(&mut featureprops, &cdc16_delta_allele, &allele_type_cvterm,
                          Some(String::from("deletion")));
    make_test_feature_rel(&mut feature_relationships, &publication,
                          &cdc16_delta_allele, &instance_of_cvterm, &cdc16_gene);

    let par1_gene = make_test_feature(&mut features, &pombe_organism, &gene_cvterm,
                                      "SPCC188.02", Some(String::from("par1")));
    let par1_delta_allele = make_test_feature(&mut features, &pombe_organism, &allele_cvterm,
                                              "SPCC188.02-allele1", Some(String::from("par1delta")));
    make_test_featureprop(&mut featureprops, &par1_delta_allele, &allele_type_cvterm,
                          Some(String::from("deletion")));
    make_test_featureprop(&mut featureprops, &par1_delta_allele, &description_cvterm,
                          Some(String::from("deletion")));
    make_test_feature_rel(&mut feature_relationships, &publication,
                          &par1_delta_allele, &instance_of_cvterm, &par1_gene);

    let genotype1 = make_test_feature(&mut features, &pombe_organism,
                                      &genotype_cvterm, "test-genotype1", None);
    make_test_feature_rel(&mut feature_relationships, &publication,
                          &par1_delta_allele, &part_of_cvterm, &genotype1);
    make_test_feature_rel(&mut feature_relationships, &publication,
                          &cdc16_allele1, &part_of_cvterm, &genotype1);
    make_test_feature_rel(&mut feature_relationships, &publication,
                          &cdc16_delta_allele, &part_of_cvterm, &genotype1);

    let par1_mrna = make_test_feature(&mut features, &pombe_organism,
                                      &mrna_cvterm, "SPCC188.02.1", None);
    make_test_feature_rel(&mut feature_relationships,  &publication,
                          &par1_mrna, &part_of_cvterm, &par1_gene);

    let par1_polypeptide = make_test_feature(&mut features, &pombe_organism, &polypeptide_cvterm,
                                             "SPCC188.02:pep", None);
    make_test_feature_rel(&mut feature_relationships, &publication,
                          &par1_polypeptide, &derives_from_cvterm, &par1_mrna);

    let par1_go0031030_fc =
        make_test_feature_cvterm(&mut feature_cvterms, &par1_mrna, &go0031030_cvterm, &publication);
    make_test_feature_cvtermprop(&mut feature_cvtermprops, &par1_go0031030_fc, &with_cvterm, "SPAC6F6.08c");

    make_test_feature_cvterm(&mut feature_cvterms, &genotype1, &pbo0022440_cvterm, &publication);

    par1_gene.featurelocs.borrow_mut().push(Rc::new(Featureloc {
        feature: par1_gene.clone(),
        fmin: 1475493,
        fmax: 1478103,
        strand: 1,
        srcfeature: chr_3.clone(),
    }));

    cdc16_gene.featurelocs.borrow_mut().push(Rc::new(Featureloc {
        feature: cdc16_gene.clone(),
        fmin: 2746666,
        fmax: 2748180,
        strand: -1,
        srcfeature: chr_1.clone(),
    }));

    pom1_gene.featurelocs.borrow_mut().push(Rc::new(Featureloc {
        feature: pom1_gene.clone(),
        fmin: 534119,
        fmax: 537869,
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
        publicationprops: vec![],
        synonyms: vec![],
        features: features,
        feature_synonyms: vec![],
        featureprops: featureprops,
        featurelocs: vec![],
        feature_cvterms: feature_cvterms,
        feature_cvtermprops: vec![],
        feature_relationships: feature_relationships,
        chadoprops: chadoprops,
    }
}

fn get_test_web_data() -> WebData {
    let raw = get_test_raw();
    let config = Config {
        extensions: vec![],
        evidence_types: HashMap::new(),
        cv_config: HashMap::new(),
    };
    let mut web_data_build = WebDataBuild::new(&raw, &config);
    web_data_build.get_web_data()
}

#[test]
fn test_gene_details() {
    let web_data = get_test_web_data();

    assert_eq!(web_data.genes.len(), 3);
    let par1_gene = web_data.genes.get("SPCC188.02").unwrap().clone();
    assert_eq!(par1_gene.uniquename, "SPCC188.02");
    assert_eq!(par1_gene.name.unwrap(), "par1");
    assert_eq!(par1_gene.cv_annotations.len(), 2);

    if par1_gene.cv_annotations.get(POMBASE_ANN_EXT_TERM_CV_NAME).is_some() {
        panic!("extension cv shouldn't be in the annotations");
    }
}

#[test]
fn test_make_publication_short() {
    let web_data = get_test_web_data();

    let pmid = "PMID:11707284";
    let ref_short = web_data.references.get(pmid).unwrap();

    assert_eq!(ref_short.uniquename, pmid);

    assert_eq!(ref_short.authors_abbrev.clone().unwrap(), "Le Goff X et al.");
    assert_eq!(ref_short.publication_year.clone().unwrap(), "2001");
}

#[test]
fn test_term_gene_count() {
    let web_data = get_test_web_data();
    let par1_gene = web_data.genes.get("SPCC188.02").unwrap().clone();
    let cv_annotations = par1_gene.cv_annotations;
    let biological_process_annotations = cv_annotations.get("biological_process").unwrap();
    assert_eq!(biological_process_annotations.len(), 1);
    let first_annotation = &biological_process_annotations[0];
    let actual_count = first_annotation.term.clone().gene_count;
    assert_eq!(actual_count, 1);
}

#[test]
fn test_gene_with() {
    let web_data = get_test_web_data();
    let par1_gene = web_data.genes.get("SPCC188.02").unwrap().clone();
    let cv_annotations = par1_gene.cv_annotations;
    let biological_process_annotations = cv_annotations.get("biological_process").unwrap();
    assert_eq!(biological_process_annotations.len(), 1);
    let first_annotation = &biological_process_annotations[0].annotations[0];
    match first_annotation.with {
        WithFromValue::Gene(ref gene_short) => {
            assert_eq!(gene_short.clone().uniquename, "SPAC6F6.08c");
            assert_eq!(gene_short.clone().name.unwrap(), "cdc16");
        },
        _ => {
            panic!();
        }
    }
}

#[test]
fn test_genotype_annotation() {
    // make sure that if a genotype has two allele from the same gene,
    // we get only one annotation
    let web_data = get_test_web_data();
    let cdc16_gene = web_data.genes.get("SPAC6F6.08c").unwrap().clone();
    let fypo_annotations = cdc16_gene.cv_annotations.get("multi_allele_phenotype").unwrap();

    assert_eq!(fypo_annotations.len(), 1);
}

#[test]
fn test_terms() {
    let web_data = get_test_web_data();

    let go0031030_cvterm = web_data.terms.get("GO:0031030").unwrap();
    assert_eq!("negative regulation of septation initiation signaling",
               go0031030_cvterm.name);
    assert_eq!(go0031030_cvterm.rel_annotations.len(), 1);

    let bp_cvterm = web_data.terms.get("GO:0008150").unwrap();
    assert_eq!("biological_process", bp_cvterm.name);
    assert_eq!(go0031030_cvterm.rel_annotations.len(), 1);
}

#[test]
fn test_locations() {
    let web_data = get_test_web_data();

    let pom1_gene = web_data.genes.get("SPAC2F7.03c").unwrap().clone();
    let pom1_loc = pom1_gene.location.unwrap().clone();

    assert_eq!(&pom1_loc.chromosome_name, "chromosome_1");
    assert_eq!(pom1_loc.start_pos, 534120);
    assert_eq!(pom1_loc.end_pos, 537869);
    assert_eq!(pom1_loc.strand, Strand::Reverse);

    assert_eq!(pom1_gene.gene_neighbourhood.len(), 2);
    assert_eq!(&pom1_gene.gene_neighbourhood[0].uniquename, "SPAC2F7.03c");
    assert_eq!(&pom1_gene.gene_neighbourhood[1].uniquename, "SPAC6F6.08c");

    let cdc16_gene = web_data.genes.get("SPAC6F6.08c").unwrap().clone();
    let cdc16_loc = cdc16_gene.location.unwrap().clone();

    assert_eq!(&cdc16_loc.chromosome_name, "chromosome_1");
    assert_eq!(cdc16_loc.start_pos, 2746667);
    assert_eq!(cdc16_loc.end_pos, 2748180);
    assert_eq!(cdc16_loc.strand, Strand::Reverse);

    assert_eq!(cdc16_gene.gene_neighbourhood.len(), 2);
    assert_eq!(&cdc16_gene.gene_neighbourhood[0].uniquename, "SPAC2F7.03c");
    assert_eq!(&cdc16_gene.gene_neighbourhood[1].uniquename, "SPAC6F6.08c");
}

