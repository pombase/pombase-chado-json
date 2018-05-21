use std::collections::*;

use std::cmp::Ordering;

extern crate pombase;

use self::pombase::types::*;
use self::pombase::web::data::*;
use self::pombase::web::config::*;
use self::pombase::web::cv_summary::*;

#[allow(dead_code)]
fn get_test_config() -> Config {
    let mut config = Config {
        database_name: "PomBase".into(),
        load_organism_taxonid: 4896,
        organisms: vec![
            ConfigOrganism {
                taxonid: 4896,
                genus: "Schizosaccharomyces".into(),
                species: "pombe".into(),
            },
            ConfigOrganism {
                taxonid: 9606,
                genus: "Homo".into(),
                species: "sapiens".into(),
            },
            ConfigOrganism {
                taxonid: 4932,
                genus: "Saccharomyces".into(),
                species: "cerevisiae".into(),
            }
        ],
        api_seq_chunk_sizes: vec![10_000, 200_000],
        extension_display_names: vec![],
        extension_relation_order: RelationOrder{
            relation_order: vec![
                String::from("directly_positively_regulates"),
                String::from("has_direct_input"),
                String::from("involved_in"),
                String::from("occurs_at"),
                String::from("occurs_in"),
                String::from("added_by"),
                String::from("added_during"),
                String::from("has_penetrance"),
            ],
            always_last: vec![String::from("happens_during"),
                              String::from("exists_during")],
        },
        evidence_types: HashMap::new(),
        cv_config: HashMap::new(),
        interesting_parents: vec![],
        viability_terms: ViabilityTerms {
            viable: "FYPO:0002058".into(),
            inviable: "FYPO:0002059".into(),
        },
        go_slim_terms: vec![],
        interpro: InterPro {
            dbnames_to_filter: vec![],
        },
        server: ServerConfig {
            subsets: ServerSubsetConfig {
                prefixes_to_remove: vec![],
            },
            solr_url: "http://localhost:8983/solr".to_owned(),
            close_synonym_boost: 0.6,
            distant_synonym_boost: 0.3,
        },
        extra_database_aliases: HashMap::new(),
        chromosomes: HashMap::new(),
    };

    config.cv_config.insert(String::from("molecular_function"),
                            CvConfig {
                                feature_type: String::from("Gene"),
                                filters: vec![],
                                split_by_parents: vec![],
                                summary_relations_to_hide: vec![],
                                summary_relation_ranges_to_collect: vec![String::from("has_substrate")],
                                sort_details_by: None,
                            });

    config
}


#[test]
fn test_compare_ext_part_with_config() {
    let config = get_test_config();
    let mut ext_part1 = ExtPart {
        rel_type_name: String::from("has_direct_input"),
        rel_type_display_name: String::from("NA"),
        ext_range: ExtRange::Misc(String::from("misc_ext_part_1")),
    };
    let mut ext_part2 = ExtPart {
        rel_type_name: String::from("has_direct_input"),
        rel_type_display_name: String::from("NA"),
        ext_range: ExtRange::Misc(String::from("misc_ext_part_2")),
    };
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Equal);

    ext_part1.rel_type_name = "directly_positively_regulates".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "has_direct_input".into();
    ext_part2.rel_type_name = "directly_positively_regulates".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part2.rel_type_name = "absent_during".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part2.rel_type_name = "misc_rel".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "other_misc_rel".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = "other_misc_rel".into();
    ext_part2.rel_type_name = "other_misc_rel".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Equal);

    ext_part2.rel_type_name = "happens_during".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "happens_during".into();
    ext_part2.rel_type_name = "misc_rel".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = "has_direct_input".into();
    ext_part2.rel_type_name = "happens_during".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "happens_during".into();
    ext_part2.rel_type_name = "has_direct_input".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = "happens_during".into();
    ext_part2.rel_type_name = "exists_during".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "happens_during".into();
    ext_part2.rel_type_name = "happens_during".into();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Equal);
}

#[allow(dead_code)]
fn make_test_ext_part(rel_type_name: &str, rel_type_display_name: &str,
                      ext_range: ExtRange) -> ExtPart {
    ExtPart {
        rel_type_name: rel_type_name.into(),
        rel_type_display_name: rel_type_display_name.into(),
        ext_range: ext_range,
    }
}

#[allow(dead_code)]
fn get_test_annotation_details_map() -> IdOntAnnotationDetailMap {
    let mut map = HashMap::new();
    map.insert(188_448,
               make_one_detail(188_448, "SPBC11B10.09", "PMID:3322810", None,
                               "IDA", vec![], HashSet::new()));
    map.insert(202_017,
               make_one_detail(202_017,"SPBC11B10.09", "PMID:2665944", None,
                               "IDA", vec![], HashSet::new()));


    let mut test_conditions = HashSet::new();
    test_conditions.insert("PECO:0000103".into());
    test_conditions.insert("PECO:0000137".into());

    let mut fypo_details = vec![
        make_one_detail(41_717, "SPBC11B10.09", "PMID:9242669", None,
                        "IDA",vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC646.13".into())), //  sds23
                        ], HashSet::new()),
        make_one_detail(41_718, "SPBC11B10.09", "PMID:11937031", None,
                        "IDA", vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC32F12.09".into())), // no name
                        ], HashSet::new()),
        make_one_detail(187_893, "SPBC11B10.09", "PMID:19523829", None, "IMP",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC6B1.04".into())), //  mde4
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:1902845".into())),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term("GO:0000089".into())),
                        ],
                        HashSet::new()),
        make_one_detail(193_221, "SPBC11B10.09", "PMID:10921876", None, "IMP",
                        vec![
                            make_test_ext_part("directly_negatively_regulates", "directly inhibits",
                                               ExtRange::Gene("SPAC144.13c".into())), //  srw1
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:1903693".into())),
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:1905785".into())),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term("GO:0000080".into())),
                        ],
                        HashSet::new()),
        make_one_detail(194_213, "SPBC11B10.09", "PMID:7957097", None, "IDA",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC776.02c".into())),  // dis2
                        ],
                        HashSet::new()),
        make_one_detail(194_661, "SPBC11B10.09", "PMID:10485849", None, "IMP",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene("SPBC146.03c".into())), //  cut3
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:1903380".into())),
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term("GO:0042307".into())),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term("GO:0000089".into())),
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

#[allow(dead_code)]
fn get_test_annotations() -> Vec<OntTermAnnotations> {
    let annotations1 = vec![188_448, 202_017];
    let ont_term1 = OntTermAnnotations {
        term: "GO:0097472".into(),
        is_not: false,
        rel_names: HashSet::new(),
        annotations: annotations1,
        summary: None,
    };

    let annotations2 =
        vec![41_717, 41_718, 187_893, 193_221, 194_213, 194_661];

    let ont_term2 = OntTermAnnotations {
        term: "GO:0004693".into(),
        is_not: false,
        rel_names: HashSet::new(),
        annotations: annotations2,
        summary: None,
    };

    vec![ont_term1, ont_term2]
}

#[allow(dead_code)]
fn make_one_detail(id: i32, gene_uniquename: &str, reference_uniquename: &str,
                   maybe_genotype_uniquename: Option<&str>, evidence: &str,
                   extension: Vec<ExtPart>,
                   conditions: HashSet<TermId>) -> OntAnnotationDetail {
    OntAnnotationDetail {
        id: id,
        genes: vec![gene_uniquename.into()],
        genotype: maybe_genotype_uniquename.map(str::to_string),
        reference: Some(reference_uniquename.into()),
        evidence: Some(evidence.into()),
        withs: vec![],
        froms: vec![],
        residue: None,
        qualifiers: vec![],
        extension: extension,
        gene_ex_props: None,
        conditions: conditions,
    }
}

#[allow(dead_code)]
fn get_test_fypo_term_details() -> Vec<i32> {
    vec![223_656, 201_099, 201_095, 204_063, 227_452, 201_094, 186_589]
}

#[allow(dead_code)]
fn make_one_genotype(display_uniquename: &str, name: Option<&str>,
                     expressed_alleles: Vec<ExpressedAllele>) -> GenotypeDetails {
    GenotypeDetails {
        display_uniquename: display_uniquename.into(),
        name: name.map(str::to_string),
        expressed_alleles: expressed_alleles,
        cv_annotations: HashMap::new(),
        references_by_uniquename: HashMap::new(),
        genes_by_uniquename: HashMap::new(),
        alleles_by_uniquename: HashMap::new(),
        terms_by_termid: HashMap::new(),
        annotation_details: HashMap::new(),
    }
}

#[allow(dead_code)]
fn make_test_gene(uniquename: &str, name: Option<&str>) -> GeneDetails {
    GeneDetails {
        uniquename: uniquename.into(),
        name: name.map(str::to_string),
        taxonid: 4896,
        product: None,
        deletion_viability: DeletionViability::Unknown,
        uniprot_identifier: None,
        interpro_matches: vec![],
        tm_domain_coords: vec![],
        orfeome_identifier: None,
        name_descriptions: vec![],
        synonyms: vec![],
        dbxrefs: HashSet::new(),
        feature_type: "gene".into(),
        characterisation_status: None,
        location: None,
        gene_neighbourhood: vec![],
        transcripts: vec![],
        cv_annotations: HashMap::new(),
        physical_interactions: vec![],
        genetic_interactions: vec![],
        ortholog_annotations: vec![],
        paralog_annotations: vec![],
        target_of_annotations: vec![],
        references_by_uniquename: HashMap::new(),
        genes_by_uniquename: HashMap::new(),
        genotypes_by_uniquename: HashMap::new(),
        alleles_by_uniquename: HashMap::new(),
        terms_by_termid: HashMap::new(),
        annotation_details: HashMap::new(),
    }
}

#[allow(dead_code)]
fn get_test_genes_map() -> UniquenameGeneMap {
    let mut ret = BTreeMap::new();

    let gene_data =
        vec![("SPBC11B10.09", Some("cdc2")),
             ("SPAC144.13c", Some("srw1")),
             ("SPAC25G10.07c", Some("cut7")),
             ("SPBC146.03c", Some("cut3")),
             ("SPBC32F12.09", None),
             ("SPBC646.13", Some("sds23")),
             ("SPBC6B1.04", Some("mde4")),
             ("SPBC776.02c", Some("dis2"))];

    for (uniquename, name) in gene_data {
        ret.insert(uniquename.into(), make_test_gene(uniquename, name));
    }

    ret
}

#[allow(dead_code)]
fn get_test_genotypes_map() -> UniquenameGenotypeMap {
    let mut ret = HashMap::new();

    ret.insert(String::from("e674fe7ceba478aa-genotype-2"),
               make_one_genotype(
                   "G799D(G799D)",
                   Some("test genotype name"),
                   vec![
                       ExpressedAllele {
                           expression: Some("Not assayed".into()),
                           allele_uniquename: "SPBC16A3.11:allele-7".into(),
                       }
                   ]
               ));

    ret.insert(String::from("d6c914796c35e3b5-genotype-4"),
               make_one_genotype(
                   "C-terminal truncation 940-1516(940-1516)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some("Not assayed".into()),
                           allele_uniquename: "SPCC1919.10c:allele-5".into(),
                       }
                   ]
               ));

    ret.insert(String::from("65c76fa511461156-genotype-3"),
               make_one_genotype(
                   "cdc25-22(c532y)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some("Not assayed".into()),
                           allele_uniquename: "SPAC24H6.05:allele-3".into(),
                       }
                   ]
               ));

    ret.insert(String::from("d6c914796c35e3b5-genotype-2"),
               make_one_genotype(
                   "ATPase dead mutant(unknown)",
                   Some("ZZ-name"),
                   vec![
                       ExpressedAllele {
                           expression: Some("Not assayed".into()),
                           allele_uniquename: "SPCC1919.10c:allele-4".into(),
                       }
                   ]
               ));

    ret.insert(String::from("d6c914796c35e3b5-genotype-3"),
               make_one_genotype(
                   "C-terminal truncation(1320-1516)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some("Not assayed".into()),
                           allele_uniquename: "SPCC1919.10c:allele-6".into(),
                       }
                   ]
               ));

    ret.insert(String::from("fd4f3f52f1d38106-genotype-4"),
               make_one_genotype(
                   "K418R(K418R)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some("Wild type product level".into()),
                           allele_uniquename: "SPAC25A8.01c:allele-5".into(),
                       }
                   ]
               ));

    ret.insert(String::from("a6d8f45c20c2227d-genotype-9"),
               make_one_genotype(
                   "UBS-I&II(F18A,F21A,W26A,L40A,W41A,W45A)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some("Not assayed".into()),
                           allele_uniquename: "SPAC3G6.02:allele-7".into(),
                       }
                   ]
               ));

    ret
}

#[allow(dead_code)]
fn make_one_allele_short(uniquename: &str, name: &str, allele_type: &str,
                         description: Option<&str>, gene_uniquename: &str) -> AlleleShort {
    AlleleShort {
        uniquename: uniquename.into(),
        description: description.map(str::to_string),
        name: Some(name.into()),
        allele_type: allele_type.into(),
        gene_uniquename: gene_uniquename.into(),
    }
}

#[allow(dead_code)]
fn get_test_alleles_map() -> UniquenameAlleleMap {
    let mut ret = HashMap::new();

    ret.insert(String::from("SPCC1919.10c:allele-4"),
               make_one_allele_short("SPCC1919.10c:allele-4", "ATPase dead mutant", "unknown", None, "SPCC1919.10c"));

    ret.insert(String::from("SPCC1919.10c:allele-5"),
               make_one_allele_short("SPCC1919.10c:allele-5", "C-terminal truncation 940-1516", "partial_amino_acid_deletion",
                                     Some("940-1516"), "SPCC1919.10c"));

    ret.insert(String::from("SPCC1919.10c:allele-6"),
               make_one_allele_short("SPCC1919.10c:allele-6", "C-terminal truncation", "partial_amino_acid_deletion", Some("1320-1516"),
                                     "SPCC1919.10c"));

    ret.insert(String::from("SPBC16A3.11:allele-7"),
               make_one_allele_short("SPBC16A3.11:allele-7", "G799D", "amino_acid_mutation", Some("G799D"), "SPBC16A3.11"));


    ret.insert(String::from("SPAC25A8.01c:allele-5"),
               make_one_allele_short("SPAC25A8.01c:allele-5", "K418R", "amino_acid_mutation", Some("K418R"), "SPAC25A8.01c"));

    ret.insert(String::from("SPAC3G6.02:allele-7"),
               make_one_allele_short("SPAC3G6.02:allele-7", "UBS-I&II", "amino_acid_mutation", Some("F18A,F21A,W26A,L40A,W41A,W45A"), "SPAC3G6.02"));

    ret.insert(String::from("SPAC24H6.05:allele-3"),
               make_one_allele_short("SPAC24H6.05:allele-3", "cdc25-22", "amino_acid_mutation", Some("C532Y"), "SPAC24H6.05"));

    ret
}

#[allow(dead_code)]
fn make_test_term_details(id: &str, name: &str, cv_name: &str) -> TermDetails {
    TermDetails {
        termid: id.into(),
        name: name.into(),
        cv_name: cv_name.into(),
        annotation_feature_type: "gene".into(),
        interesting_parents: HashSet::new(),
        subsets: vec!["goslim_pombe".into()],
        synonyms: vec![],
        definition: None,
        direct_ancestors: vec![],
        genes_annotated_with: HashSet::new(),
        is_obsolete: false,
        single_allele_genotype_uniquenames: HashSet::new(),
        cv_annotations: HashMap::new(),
        genes_by_uniquename: HashMap::new(),
        genotypes_by_uniquename: HashMap::new(),
        alleles_by_uniquename: HashMap::new(),
        references_by_uniquename: HashMap::new(),
        terms_by_termid: HashMap::new(),
        annotation_details: HashMap::new(),
        gene_count: 0,
        genotype_count: 0,
    }
}

#[allow(dead_code)]
fn get_test_terms_map() -> TermIdDetailsMap {
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

    ret
}

#[test]
fn test_cmp_ont_annotation_detail() {
    let mut details_vec = get_test_fypo_term_details();
    let genes = get_test_genes_map();
    let genotypes = get_test_genotypes_map();
    let alleles = get_test_alleles_map();
    let terms = get_test_terms_map();
    let annotation_details_maps = get_test_annotation_details_map();

    let cmp_detail_with_genotypes =
        |id1: &i32, id2: &i32| {
            let annotation1 = annotation_details_maps.get(id1).expect(&format!("{}", id1));
            let annotation2 = annotation_details_maps.get(id2).expect(&format!("{}", id2));
            let cv_config = &get_test_config().cv_config_by_name("molecular_function");
            pombase::web::data_build::cmp_ont_annotation_detail(cv_config,
                                      annotation1, annotation2, &genes,
                                      &genotypes, &terms).unwrap()
        };

    details_vec.sort_by(&cmp_detail_with_genotypes);

    let expected: Vec<String> =
        vec!["atpase_dead_mutant-unknown-unknown-expression-not_assayed",
        "c-terminal_truncation-1320-1516-partial_amino_acid_deletion-expression-not_assayed",
        "c-terminal_truncation_940-1516-940-1516-partial_amino_acid_deletion-expression-not_assayed",
        "cdc25-22-c532y-amino_acid_mutation-expression-not_assayed",
        "g799d-g799d-amino_acid_mutation-expression-not_assayed",
        "k418r-k418r-amino_acid_mutation-expression-wild_type_product_level",
        "ubs-i&ii-f18a,f21a,w26a,l40a,w41a,w45a-amino_acid_mutation-expression-not_assayed"]
        .iter().map(|s| str::to_string(s)).collect();

    let mut result_genotype_display_names =
        details_vec.drain(0..)
        .map(|detail_id| {
            let detail = annotation_details_maps.get(&detail_id).unwrap();
            let genotype_uniquename = detail.clone().genotype.unwrap();
            let genotype = genotypes.get(&genotype_uniquename).unwrap();
            pombase::web::data_build::make_genotype_display_name(&genotype.expressed_alleles, &alleles)
                .to_lowercase()
        }).collect::<Vec<String>>();

    result_genotype_display_names.sort();

    assert_eq!(result_genotype_display_names, expected);

    let test_term_annotations = get_test_annotations();
    let mut extension_details_vec = test_term_annotations[1].annotations.clone();

    extension_details_vec.sort_by(&cmp_detail_with_genotypes);

    let annotation_sort_results: Vec<(String, String)> =
        extension_details_vec.iter().map(|detail_id| {
            let detail = annotation_details_maps.get(detail_id).unwrap();
            ((*detail).genes[0].clone(),
             (*detail).reference.clone().unwrap())
        }).collect();

    let expected_annotation_sort: Vec<(String, String)> =
        vec![("SPBC11B10.09", "PMID:10921876"),
             ("SPBC11B10.09", "PMID:10485849" /* has_direct_input(cut3) */),
             ("SPBC11B10.09", "PMID:7957097" /* has_direct_input(dis2) */),
             ("SPBC11B10.09", "PMID:19523829" /* has_direct_input(mde4), part_of(...), happens_during(...) */),
             ("SPBC11B10.09", "PMID:9242669" /* has_direct_input(sds23) */),
             ("SPBC11B10.09", "PMID:11937031" /* has_direct_input(SPBC32F12.09) */),
        ]
        .iter()
        .map(|&(gene, reference)|
             (gene.into(), reference.into())).collect();

    assert_eq![annotation_sort_results, expected_annotation_sort];
}

#[allow(dead_code)]
fn make_test_summary(termid: &str, rows: Vec<TermSummaryRow>) -> OntTermAnnotations {
    OntTermAnnotations {
        term: termid.to_owned(),
        is_not: false,
        annotations: vec![],
        rel_names: HashSet::new(),
        summary: Some(rows),
    }
}


#[test]
fn test_summary_row_equals() {
    let r1 = TermSummaryRow {
        gene_uniquenames: vec!["SPAPB1A10.09".into()],
        genotype_uniquenames: vec![],
        extension: vec![],
    };
    let r2 = TermSummaryRow {
        gene_uniquenames: vec!["SPAPB1A10.09".into()],
        genotype_uniquenames: vec![],
        extension: vec![],
    };
    assert!(r1 == r2);
}

#[allow(dead_code)]
fn get_test_summaries() -> Vec<OntTermAnnotations> {
    let mut summaries = vec![];

    let ext = make_test_ext_part("part_of", "involved in",
                                 ExtRange::Term("GO:1905785".into()));
    let ext2 = make_test_ext_part("some_rel", "some_rel_display_name",
                                  ExtRange::Term("GO:1234567".into()));

    summaries.push(make_test_summary("GO:0022403", vec![]));
    summaries.push(make_test_summary("GO:0051318",
                                     vec![TermSummaryRow {
                                         gene_uniquenames: vec![],
                                         genotype_uniquenames: vec![],
                                         extension: vec![ext.clone()],
                                     }]));
    summaries.push(make_test_summary("GO:0000080", vec![]));
    summaries.push(make_test_summary("GO:0000089",
                                     vec![
                                         TermSummaryRow {
                                             gene_uniquenames: vec![],
                                             genotype_uniquenames: vec![],
                                             extension: vec![ext.clone()],
                                         }
                                     ]));
    summaries.push(make_test_summary("GO:0099999",
                                     vec![
                                         TermSummaryRow {
                                             gene_uniquenames: vec![],
                                             genotype_uniquenames: vec![],
                                             extension: vec![ext.clone(), ext2],
                                         }
                                     ]));

    summaries
}

#[allow(dead_code)]
fn get_test_children_by_termid() -> HashMap<TermId, HashSet<TermId>> {
    let mut children_by_termid = HashMap::new();

    let mut children_of_0022403 = HashSet::new();
    children_of_0022403.insert("GO:0051318".into());
    children_of_0022403.insert("GO:0000080".into());
    children_of_0022403.insert("GO:0088888".into());
    children_of_0022403.insert("GO:0000089".into());
    children_of_0022403.insert("GO:0099999".into());

    let mut children_of_0051318 = HashSet::new();
    children_of_0051318.insert("GO:0000080".into());
    children_of_0051318.insert("GO:0088888".into());
    children_of_0051318.insert("GO:0000089".into());
    children_of_0051318.insert("GO:0099999".into());

    let mut children_of_0000080 = HashSet::new();
    children_of_0000080.insert("GO:0088888".into());

    let mut children_of_0000089 = HashSet::new();
    children_of_0000089.insert("GO:0099999".into());

    children_by_termid.insert("GO:0022403".into(), children_of_0022403);
    children_by_termid.insert("GO:0051318".into(), children_of_0051318);
    children_by_termid.insert("GO:0000080".into(), children_of_0000080);
    children_by_termid.insert("GO:0000089".into(), children_of_0000089);

    children_by_termid
}

#[test]
fn test_remove_redundant_summaries() {
    let mut term_annotations: Vec<OntTermAnnotations> = get_test_summaries();

    let children_by_termid = get_test_children_by_termid();
    assert_eq!(term_annotations.len(), 5);

    pombase::web::cv_summary::remove_redundant_summaries(&children_by_termid, &mut term_annotations);

    assert_eq!(term_annotations.iter().filter(|term_annotation| {
        term_annotation.summary.is_some()
    }).collect::<Vec<&OntTermAnnotations>>().len(), 3);
}


fn get_test_gene_short_map() -> IdGeneShortMap {
    let mut ret_map = HashMap::new();

    ret_map.insert("SPAC977.09c".into(),
                   GeneShort {
                       uniquename: "SPAC977.09c".into(),
                       name: None,
                       product: Some("phospholipase (predicted)".into()),
                   });
    ret_map.insert("SPAC3G9.09c".into(),
                   GeneShort {
                       uniquename: "SPAC3G9.09c".into(),
                       name: Some("tif211".into()),
                       product: Some("translation initiation factor eIF2 alpha subunit".into()),
                   });
    ret_map.insert("SPAC16.01".into(),
                   GeneShort {
                       uniquename: "SPAC16.01".into(),
                       name: Some("rho2".into()),
                       product: Some("Rho family GTPase Rho2".into()),
                   });
    ret_map.insert("SPAC24C9.02c".into(),
                   GeneShort { 
                       uniquename: "SPAC24C9.02c".into(),
                       name: Some("cyt2".into()),
                       product: Some("cytochrome c1 heme lyase Cyt2 (predicted)".into()),
                   });

    ret_map
}




#[test]
fn test_merge_ext_part_ranges() {
    let ext_part1 = ExtPart {
        rel_type_name: "has_substrate".into(),
        rel_type_display_name: "has substrate".into(),
        ext_range: ExtRange::SummaryGenes(vec![vec!["SPAC977.09c".into()]]),
    };
    let ext_part2 = ExtPart {
        rel_type_name: "has_substrate".into(),
        rel_type_display_name: "has substrate".into(),
        ext_range: ExtRange::SummaryGenes(vec![vec!["SPAC24C9.02c".into()]]),
    };

    let gene_short_map = get_test_gene_short_map();
    let res = merge_ext_part_ranges(&ext_part1, &ext_part2, &gene_short_map);

    assert_eq!(res.ext_range,
               ExtRange::SummaryGenes(vec![vec!["SPAC24C9.02c".into()],
                                           vec!["SPAC977.09c".into()]]));
}

fn get_test_summary_rows() -> Vec<TermSummaryRow> {
    let mut rows = vec![];

    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC25B2.01c")],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC25B2.01c")],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: String::from("some_good_rel"),
                rel_type_display_name: String::from("some rel"),
                ext_range: ExtRange::Term(String::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: String::from("has_substrate"),
                rel_type_display_name: String::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![String::from("SPAC3G9.09c")]]),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: String::from("has_substrate"),
                rel_type_display_name: String::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![String::from("SPAC3G9.09c")]]),
            },
            ExtPart {
                rel_type_name: String::from("during"),
                rel_type_display_name: String::from("during"),
                ext_range: ExtRange::Term(String::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: String::from("has_substrate"),
                rel_type_display_name: String::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![String::from("SPAC3G9.09c")]]),
            },
            ExtPart {
                rel_type_name: String::from("during"),
                rel_type_display_name: String::from("during"),
                ext_range: ExtRange::Term(String::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC1786.03")], // change annotated gene
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: String::from("has_substrate"),
                rel_type_display_name: String::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![String::from("SPAC3G9.09c")]]),
            },
            ExtPart {
                rel_type_name: String::from("during"),
                rel_type_display_name: String::from("during"),
                ext_range: ExtRange::Term(String::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC1786.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: String::from("has_substrate"),
                rel_type_display_name: String::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![String::from("SPAC16.01")]]),   // change substrate
            },
            ExtPart {
                rel_type_name: String::from("during"),
                rel_type_display_name: String::from("during"),
                ext_range: ExtRange::Term(String::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![String::from("SPAC1786.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: String::from("has_substrate"),
                rel_type_display_name: String::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![String::from("SPAC16.01")]]),
            },
            ExtPart {
                rel_type_name: String::from("during"),
                rel_type_display_name: String::from("during"),
                ext_range: ExtRange::Term(String::from("GO:0071472")), // change during term
            }],
    });

    rows
}

#[test]
fn test_collect_ext_summary_genes() {
    let config = get_test_config();

    let mut rows = get_test_summary_rows();
    assert_eq!(rows.len(), 10);

    let gene_short_map = get_test_gene_short_map();

    pombase::web::cv_summary::collect_ext_summary_genes(&config.cv_config_by_name("molecular_function"),
                                                        &mut rows, &gene_short_map);
    assert_eq!(rows.len(), 8);

    let collected_ext = rows.get(6).unwrap();
    let collected_ext_ext_part_1 = collected_ext.extension.get(0).unwrap();
    let summary_genes_vec = vec![vec![String::from("SPAC16.01")],
                                 vec![String::from("SPAC3G9.09c")]];
    assert_eq!(collected_ext_ext_part_1.ext_range,
               ExtRange::SummaryGenes(summary_genes_vec));
}

#[test]
fn test_remove_redundant_summary_rows() {
    let mut rows = get_test_summary_rows();
    assert_eq!(rows.len(), 10);

    pombase::web::cv_summary::remove_redundant_summary_rows(&mut rows);
    assert_eq!(rows.len(), 6);
}

