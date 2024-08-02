use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use flexstr::ToSharedStr;

use pombase::api_data::{APIData, api_maps_from_file};
use pombase::data_types::{APIGenotypeAnnotation, AlleleDetails, DeletionViability, DisplayUniquenameGenotypeMap, ExpressedAllele, ExtPart, ExtRange, GeneDetails, GenotypeDetails, GenotypeLocus, IdGenotypeMap, IdOntAnnotationDetailMap, OntAnnotationDetail, Ploidiness, TermDetails, TermIdDetailsMap, Throughput, UniquenameAlleleDetailsMap, UniquenameAlleleMap, UniquenameGeneMap, UniquenameReferenceMap};
use pombase::types::TermId;
use pombase::utils::{make_maps_database_tables, store_maps_into_database};
use pombase::web::config::Config;
use rusqlite::Connection;

#[allow(dead_code)]
pub fn setup_test_maps_database(mut conn: &mut Connection,
                                terms: &TermIdDetailsMap,
                                genes: &UniquenameGeneMap,
                                alleles: &UniquenameAlleleMap,
                                references: &UniquenameReferenceMap,
                                genotypes: &IdGenotypeMap,
                                annotation_details: &IdOntAnnotationDetailMap,
                                termid_genotype_annotation: &HashMap<TermId, Vec<APIGenotypeAnnotation>>) {
    make_maps_database_tables(conn).unwrap();

    store_maps_into_database(&mut conn, terms, genes, alleles, references,
                             genotypes, annotation_details, termid_genotype_annotation).unwrap();
}

#[allow(dead_code)]
pub fn get_api_data() -> APIData {
    use std::path::PathBuf;
    let mut search_maps_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    search_maps_path.push("tests/test_search_data.json.zst");
    let mut config_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    config_path.push("tests/test_config.json");
    let config = Config::read(config_path.to_str().expect("config"));
    let api_maps = api_maps_from_file(search_maps_path.to_str().expect("search maps"));
    let mut maps_db_conn = Connection::open_in_memory().unwrap();
    let genes = get_test_genes_map();
    let genotypes = get_test_genotypes_map();
    let alleles = get_test_alleles_map();
    let terms = get_test_terms_map();
    let references = get_test_references_map();
    let annotation_details_maps = get_test_annotation_details_map();
    let termid_genotype_annotation = HashMap::new();
    setup_test_maps_database(&mut maps_db_conn, &terms, &genes, &alleles, &references, &genotypes,
                             &annotation_details_maps, &termid_genotype_annotation);
    APIData::new(&config, maps_db_conn, api_maps)
}

#[warn(dead_code)]
pub fn make_one_genotype(display_uniquename: &str, name: Option<&str>,
                     loci: Vec<GenotypeLocus>) -> GenotypeDetails {
    GenotypeDetails {
        display_uniquename: display_uniquename.into(),
        display_name: display_uniquename.into(),
        name: name.map(|s| s.to_shared_str()),
        taxonid: 4896,
        loci: loci,
        ploidiness: Ploidiness::Haploid,
        comment: None,
        cv_annotations: HashMap::new(),
        double_mutant_genetic_interactions: HashMap::new(),
        rescue_genetic_interactions: HashMap::new(),
        references_by_uniquename: HashMap::new(),
        genes_by_uniquename: HashMap::new(),
        alleles_by_uniquename: HashMap::new(),
        transcripts_by_uniquename: HashMap::new(),
        genotypes_by_uniquename: HashMap::new(),
        terms_by_termid: HashMap::new(),
        annotation_details: HashMap::new(),
        annotation_count: 0,
    }
}

#[warn(dead_code)]
pub fn make_test_gene(uniquename: &str, name: Option<&str>) -> GeneDetails {
    GeneDetails {
        uniquename: uniquename.into(),
        name: name.map(|s| s.to_shared_str()),
        taxonid: 4896,
        product: None,
        deletion_viability: DeletionViability::Unknown,
        uniprot_identifier: None,
        secondary_identifier: None,
        biogrid_interactor_id: None,
        rnacentral_urs_identifier: None,
        rnacentral_2d_structure_id: None,
        pdb_entries: vec![],
        interpro_matches: vec![],
        tm_domain_coords: vec![],
        disordered_region_coords: vec![],
        low_complexity_region_coords: vec![],
        coiled_coil_coords: vec![],
        signal_peptide: None,
        transit_peptide: None,
        binding_sites: vec![],
        active_sites: vec![],
        beta_strands: vec![],
        has_protein_features: false,
        rfam_annotations: vec![],
        orfeome_identifier: None,
        pombephosphoproteomics_unige_ch_starvation_mating_gene: None,
        pombephosphoproteomics_unige_ch_fusion_gene: None,
        name_descriptions: vec![],
        synonyms: vec![],
        dbxrefs: HashSet::new(),
        gocams: HashSet::new(),
        flags: HashSet::new(),
        feature_type: "gene".to_shared_str(),
        feature_so_termid: "SO:0000704".to_shared_str(),
        transcript_so_termid: Some("SO:0001217".to_shared_str()),
        characterisation_status: None,
        taxonomic_distribution: None,
        location: None,
        gene_neighbourhood: vec![],
        transcripts: vec![],
        transcripts_by_uniquename: HashMap::new(),
        cv_annotations: HashMap::new(),
        physical_interactions: vec![],
        genetic_interactions: HashMap::new(),
        ortholog_annotations: vec![],
        paralog_annotations: vec![],
        target_of_annotations: vec![],
        references_by_uniquename: HashMap::new(),
        genes_by_uniquename: HashMap::new(),
        genotypes_by_uniquename: HashMap::new(),
        alleles_by_uniquename: HashMap::new(),
        terms_by_termid: HashMap::new(),
        annotation_details: HashMap::new(),
        feature_publications: HashSet::new(),
        subset_termids: HashSet::new(),
        gene_history: vec![],
    }
}

#[cfg(test)]
pub fn get_test_genes_map() -> UniquenameGeneMap {
    use pombase::data_types::OntTermAnnotations;

    let mut ret = BTreeMap::new();

    let gene_data =
        vec![("SPBC11B10.09", Some("cdc2")),
             ("SPAC144.13c", Some("srw1")),
             ("SPAC25G10.07c", Some("cut7")),
             ("SPBC146.03c", Some("cut3")),
             ("SPBC32F12.09", None),
             ("SPBC646.13", Some("sds23")),
             ("SPBC6B1.04", Some("mde4")),
             ("SPBC776.02c", Some("dis2")),
             ("SPAC3G9.09c", Some("tif211")),
             ("SPAC977.09c", Some("plb4")),
             ("SPAC16.01", Some("rho2")),
             ("SPAC24C9.02c", Some("cyt2")),
             ("SPAC24B11.06c", None),
             ("SPAC19G12.03", None),
             ("SPAC2F3.09", None)];

    for (uniquename, name) in gene_data {
        ret.insert(uniquename.into(), make_test_gene(uniquename, name));
    }

    let gene_spac27e2_05 = make_test_gene("SPAC27E2.05", None);

    let mut cv_annotations = HashMap::new();
    cv_annotations.insert("biological_process".into(),
                          vec![OntTermAnnotations {
                              term: "GO:0044237".into(),
                              is_not: false,
                              rel_names: HashSet::new(),
                              annotations: vec![10000],
                              summary: Some(vec![]),
                          }]);
    let mut transcripts_by_uniquename = HashMap::new();
    transcripts_by_uniquename.insert("SPBC11C11.05.1".into(), None);

    ret.insert("SPAC27E2.05".into(),
               GeneDetails {
                 deletion_viability: DeletionViability::DependsOnConditions,
                 cv_annotations,
                 transcripts: vec!["SPBC11C11.05.1".into()],
                 transcripts_by_uniquename,
                 ..gene_spac27e2_05
               });

    ret
}

pub fn get_test_references_map() -> UniquenameReferenceMap {
    HashMap::new()
}

#[cfg(test)]
pub fn get_test_genotypes_map() -> DisplayUniquenameGenotypeMap {
    let mut ret = HashMap::new();

    ret.insert("e674fe7ceba478aa-genotype-2".to_shared_str(),
               make_one_genotype(
                   "G799D(G799D)",
                   Some("test genotype name"),
                   vec![
                       GenotypeLocus {
                           expressed_alleles: vec![
                               ExpressedAllele {
                                   expression: Some("Not assayed".to_shared_str()),
                                   allele_uniquename: "SPBC16A3.11:allele-7".to_shared_str(),
                                   promoter_gene: None,
                               }
                           ]
                       }
                   ]
               ));

    ret.insert("d6c914796c35e3b5-genotype-4".to_shared_str(),
               make_one_genotype(
                   "C-terminal truncation 940-1516(940-1516)",
                   None,
                   vec![
                       GenotypeLocus {
                           expressed_alleles: vec![
                               ExpressedAllele {
                                   expression: Some("Not assayed".to_shared_str()),
                                   allele_uniquename: "SPCC1919.10c:allele-5".to_shared_str(),
                                   promoter_gene: None,
                               }
                           ]
                       }
                   ]
               ));

    ret.insert("65c76fa511461156-genotype-3".to_shared_str(),
               make_one_genotype(
                   "cdc25-22(c532y)",
                   None,
                   vec![
                       GenotypeLocus {
                           expressed_alleles: vec![
                               ExpressedAllele {
                                   expression: Some("Not assayed".to_shared_str()),
                                   allele_uniquename: "SPAC24H6.05:allele-3".to_shared_str(),
                                   promoter_gene: None,
                                }
                           ]
                       }
                   ]
               ));

    ret.insert("d6c914796c35e3b5-genotype-2".to_shared_str(),
               make_one_genotype(
                   "ATPase dead mutant(unknown)",
                   Some("ZZ-name"),
                   vec![
                       GenotypeLocus {
                           expressed_alleles: vec![
                               ExpressedAllele {
                                   expression: Some("Not assayed".to_shared_str()),
                                   allele_uniquename: "SPCC1919.10c:allele-4".to_shared_str(),
                                   promoter_gene: None,
                               }
                           ]
                       }
                   ]
               ));

    ret.insert("d6c914796c35e3b5-genotype-3".to_shared_str(),
               make_one_genotype(
                   "C-terminal truncation(1320-1516)",
                   None,
                   vec![
                       GenotypeLocus {
                           expressed_alleles: vec![
                               ExpressedAllele {
                                   expression: Some("Not assayed".to_shared_str()),
                                   allele_uniquename: "SPCC1919.10c:allele-6".to_shared_str(),
                                   promoter_gene: None,
                               }
                           ]
                       }
                   ]
               ));

    ret.insert("fd4f3f52f1d38106-genotype-4".to_shared_str(),
               make_one_genotype(
                   "K418R(K418R)",
                   None,
                   vec![
                       GenotypeLocus {
                           expressed_alleles: vec![
                               ExpressedAllele {
                                   expression: Some("Wild type product level".to_shared_str()),
                                   allele_uniquename: "SPAC25A8.01c:allele-5".to_shared_str(),
                                   promoter_gene: None,
                              }
                           ]
                       }
                   ]
               ));

    ret.insert("a6d8f45c20c2227d-genotype-9".to_shared_str(),
               make_one_genotype(
                   "UBS-I&II(F18A,F21A,W26A,L40A,W41A,W45A)",
                   None,
                   vec![
                       GenotypeLocus {
                           expressed_alleles: vec![
                               ExpressedAllele {
                                   expression: Some("Not assayed".to_shared_str()),
                                   allele_uniquename: "SPAC3G6.02:allele-7".to_shared_str(),
                                   promoter_gene: None,
                               }
                           ]
                       }
                   ]
               ));

    ret
}

pub fn make_one_allele_details(uniquename: &str, name: &str, allele_type: &str,
                           description: Option<&str>, gene_uniquename: &str) -> AlleleDetails {
     let ref gene = make_test_gene(gene_uniquename, None);
     AlleleDetails::new(
        uniquename.into(),
        &Some(name.into()),
        allele_type,
        &description.map(|s| s.to_shared_str()),
        &vec![],
        false,
        gene.into(),
    )
}

pub fn get_test_alleles_map() -> UniquenameAlleleDetailsMap {
    let mut ret = BTreeMap::new();

    ret.insert("SPCC1919.10c:allele-4".to_shared_str(),
               make_one_allele_details("SPCC1919.10c:allele-4", "ATPase dead mutant", "unknown", None, "SPCC1919.10c"));

    ret.insert("SPCC1919.10c:allele-5".to_shared_str(),
               make_one_allele_details("SPCC1919.10c:allele-5", "C-terminal truncation 940-1516", "partial_amino_acid_deletion",
                                     Some("940-1516"), "SPCC1919.10c"));

    ret.insert("SPCC1919.10c:allele-6".to_shared_str(),
               make_one_allele_details("SPCC1919.10c:allele-6", "C-terminal truncation", "partial_amino_acid_deletion", Some("1320-1516"),
                                     "SPCC1919.10c"));

    ret.insert("SPBC16A3.11:allele-7".to_shared_str(),
               make_one_allele_details("SPBC16A3.11:allele-7", "G799D", "amino_acid_mutation", Some("G799D"), "SPBC16A3.11"));


    ret.insert("SPAC25A8.01c:allele-5".to_shared_str(),
               make_one_allele_details("SPAC25A8.01c:allele-5", "K418R", "amino_acid_mutation", Some("K418R"), "SPAC25A8.01c"));

    ret.insert("SPAC3G6.02:allele-7".to_shared_str(),
               make_one_allele_details("SPAC3G6.02:allele-7", "UBS-I&II", "amino_acid_mutation", Some("F18A,F21A,W26A,L40A,W41A,W45A"), "SPAC3G6.02"));

    ret.insert("SPAC24H6.05:allele-3".to_shared_str(),
               make_one_allele_details("SPAC24H6.05:allele-3", "cdc25-22", "amino_acid_mutation", Some("C532Y"), "SPAC24H6.05"));

    ret
}

pub fn make_test_term_details(id: &str, name: &str, cv_name: &str) -> TermDetails {
    TermDetails {
        termid: id.into(),
        name: name.into(),
        cv_name: cv_name.into(),
        annotation_feature_type: "gene".to_shared_str(),
        interesting_parent_ids: HashSet::new(),
        interesting_parent_details: HashSet::new(),
        in_subsets: HashSet::from_iter(vec!["goslim_pombe".to_shared_str()]),
        synonyms: vec![],
        definition: None,
        direct_ancestors: vec![],
        definition_xrefs: HashSet::new(),
        secondary_identifiers: HashSet::new(),
        annotated_genes: HashSet::new(),
        single_locus_annotated_genes: HashSet::new(),
        multi_locus_annotated_genes: HashSet::new(),
        is_obsolete: false,
        single_locus_genotype_uniquenames: HashSet::new(),
        cv_annotations: HashMap::new(),
        genes_by_uniquename: HashMap::new(),
        genotypes_by_uniquename: HashMap::new(),
        alleles_by_uniquename: HashMap::new(),
        transcripts_by_uniquename: HashMap::new(),
        references_by_uniquename: HashMap::new(),
        terms_by_termid: HashMap::new(),
        annotation_details: HashMap::new(),
        double_mutant_genetic_interactions: HashMap::new(),
        single_allele_genetic_interactions: HashMap::new(),
        gene_count: 0,
        genotype_count: 0,
        xrefs: HashMap::new(),
        pombase_gene_id: None,
        gocams: HashSet::new(),
    }
}

pub fn get_test_terms_map() -> TermIdDetailsMap {
    let mut ret = HashMap::new();

    let term_data = vec![
        ("GO:0022403", "cell cycle phase", "biological_process"),
        ("GO:0051318", "G1 phase", "biological_process"),
        ("GO:0000080", "mitotic G1 phase", "biological_process"),
        ("GO:0088888", "fake child of mitotic G1 phase", "biological_process"),
        ("GO:0000089", "mitotic metaphase", "biological_process"),
        ("GO:0099999", "fake child of  mitotic metaphase term", "biological_process"),
        ("GO:0042307", "positive regulation of protein import into nucleus", "biological_process"),
        ("GO:0098783", "correction of merotelic kinetochore attachment, mitotic", "biological_process"),
        ("GO:1902845", "negative regulation of mitotic spindle elongation", "biological_process"),
        ("GO:1903380", "positive regulation of mitotic chromosome condensation", "biological_process"),
        ("GO:1903693", "regulation of mitotic G1 cell cycle arrest in response to nitrogen starvation", "biological_process"),
        ("GO:1905785", "negative regulation of anaphase-promoting complex-dependent catabolic process", "biological_process"),
    ];

    for (id, name, cv_name) in term_data {
        ret.insert(id.into(), make_test_term_details(id, name, cv_name));
    }

    // this term is used to test querying:
    let term_0044237 = make_test_term_details("GO:0044237", "cellular metabolic process",
                                              "biological_process");

    ret.insert("GO:0044237".into(), TermDetails {
      is_obsolete: false,
      gene_count: 10,
      genotype_count: 0,
      annotation_feature_type: "gene".into(),
      ..term_0044237
    });

    ret
}

#[allow(dead_code)]
fn make_one_detail(id: i32, gene_uniquename: &str, reference_uniquename: &str,
                   maybe_genotype_uniquename: Option<&str>, evidence: &str,
                   extension: Vec<ExtPart>,
                   conditions: HashSet<TermId>) -> OntAnnotationDetail {
    let condition_details: BTreeSet<_> =
        conditions.iter().map(|cond| (cond.clone(), None)).collect();

    OntAnnotationDetail {
        id: id,
        genes: vec![gene_uniquename.into()],
        transcript_uniquenames: vec![],
        genotype: maybe_genotype_uniquename.map(|s| s.to_shared_str()),
        genotype_background: None,
        reference: Some(reference_uniquename.into()),
        evidence: Some(evidence.into()),
        eco_evidence: Some(evidence.into()),
        annotation_phenotype_score: None,
        withs: HashSet::new(),
        froms: HashSet::new(),
        date: Some("2009-02-13".to_shared_str()),
        residue: None,
        gene_product_form_id: None,
        allele_promoters: vec![],
        qualifiers: vec![],
        extension: extension,
        gene_ex_props: None,
        conditions,
        condition_details,
        assigned_by: Some("PomBase".to_shared_str()),
        throughput: Some(Throughput::HighThroughput),
        curator: None,
    }
}

#[allow(dead_code)]
pub fn make_test_ext_part(rel_type_name: &str, rel_type_display_name: &str,
                          ext_range: ExtRange) -> ExtPart {
    ExtPart {
        rel_type_id: Some("RO:0000000".to_shared_str()),
        rel_type_name: rel_type_name.into(),
        rel_type_display_name: rel_type_display_name.into(),
        ext_range: ext_range,
    }
}

#[allow(dead_code)]
pub fn get_test_annotation_details_map() -> IdOntAnnotationDetailMap {
    let mut map = HashMap::new();
    map.insert(188_448,
               make_one_detail(188_448, "SPBC11B10.09", "PMID:3322810", None,
                               "IDA", vec![], HashSet::new()));
    map.insert(202_017,
               make_one_detail(202_017,"SPBC11B10.09", "PMID:2665944", None,
                               "IDA", vec![], HashSet::new()));

    map.insert(10000, make_one_detail(10000, "SPAC27E2.05",
                                     "PB_REF:0000001", None,
                                     "ISS", vec![], HashSet::new()));

    let mut test_conditions = HashSet::new();
    test_conditions.insert("FYECO:0000103".to_shared_str());
    test_conditions.insert("FYECO:0000137".to_shared_str());

    let mut fypo_details = vec![
        make_one_detail(41_717, "SPBC11B10.09", "PMID:9242669", None,
                        "IDA",vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC646.13".to_shared_str())), //  sds23
                        ], HashSet::new()),
        make_one_detail(41_718, "SPBC11B10.09", "PMID:11937031", None,
                        "IDA", vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC32F12.09".to_shared_str())), // no name
                        ], HashSet::new()),
        make_one_detail(187_893, "SPBC11B10.09", "PMID:19523829", None, "IMP",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC6B1.04".to_shared_str())), //  mde4
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:1902845".to_shared_str())),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term("GO:0000089".to_shared_str())),
                        ],
                        HashSet::new()),
        make_one_detail(193_221, "SPBC11B10.09", "PMID:10921876", None, "IMP",
                        vec![
                            make_test_ext_part("directly_negatively_regulates", "directly inhibits",
                                               ExtRange::Gene("SPAC144.13c".to_shared_str())), //  srw1
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:1903693".to_shared_str())),
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:1905785".to_shared_str())),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term("GO:0000080".to_shared_str())),
                        ],
                        HashSet::new()),
        make_one_detail(194_213, "SPBC11B10.09", "PMID:7957097", None, "IDA",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC776.02c".to_shared_str())),  // dis2
                        ],
                        HashSet::new()),
        make_one_detail(194_661, "SPBC11B10.09", "PMID:10485849", None, "IMP",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC146.03c".to_shared_str())), //  cut3
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:1903380".to_shared_str())),
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:0042307".to_shared_str())),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term("GO:0000089".to_shared_str())),
                        ],
                        HashSet::new()),
        make_one_detail(223_656,
                        "SPBC16A3.11",
                        "PMID:23050226",
                        Some("e674fe7ceba478aa-genotype-2"),
                        "Cell growth assay",
                        vec![],
                        test_conditions.clone()),
        make_one_detail(201_099,
                        "SPCC1919.10c",
                        "PMID:16421926",
                        Some("d6c914796c35e3b5-genotype-4"),
                        "Cell growth assay",
                        vec![],
                        HashSet::new()),
        make_one_detail(201_095,
                        "SPCC1919.10c",
                        "PMID:16421926",
                        Some("d6c914796c35e3b5-genotype-3"),
                        "Cell growth assay",
                        vec![],
                        HashSet::new()),
        make_one_detail(204_063,
                        "SPAC25A8.01c",
                        "PMID:25798942",
                        Some("fd4f3f52f1d38106-genotype-4"),
                        "Cell growth assay",
                        vec![],
                        test_conditions.clone()),
        make_one_detail(227_452,
                        "SPAC3G6.02",
                        "PMID:25306921",
                        Some("a6d8f45c20c2227d-genotype-9"),
                        "Cell growth assay",
                        vec![],
                        HashSet::new()),
        make_one_detail(201_094,
                        "SPCC1919.10c",
                        "PMID:16421926",
                        Some("d6c914796c35e3b5-genotype-2"),
                        "Cell growth assay",
                        vec![],
                        HashSet::new()),
        make_one_detail(186_589,
                        "SPAC24H6.05",
                        "PMID:1464319",
                        Some("65c76fa511461156-genotype-3"),
                        "Cell growth assay",
                        vec![],
                        test_conditions)];

    for detail in fypo_details.drain(0..) {
        map.insert(detail.id, detail);
    }

    map
}
