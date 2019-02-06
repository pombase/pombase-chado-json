use std::collections::BTreeMap;
use hashbrown::{HashMap, HashSet};
use std::iter::FromIterator;

use std::cmp::Ordering;

extern crate pombase;
extern crate pombase_rc_string;

use self::pombase::types::*;
use self::pombase::web::data::*;
use self::pombase::web::config::*;
use self::pombase::web::cv_summary::*;

use pombase_rc_string::RcString;

#[allow(dead_code)]
fn get_test_config() -> Config {
    let mut config = Config {
        database_name: "PomBase".into(),
        database_citation: RcString::from("PMID:22039153"),
        load_organism_taxonid: Some(4896),
        base_url: RcString::from("https://www.pombase.org"),
        organisms: vec![
            ConfigOrganism {
                taxonid: 4896,
                genus: RcString::from("Schizosaccharomyces"),
                species: RcString::from("pombe"),
                assembly_version: Some(RcString::from("ASM294v2")),
            },
            ConfigOrganism {
                taxonid: 9606,
                genus: RcString::from("Homo"),
                species: RcString::from("sapiens"),
                assembly_version: None,
            },
            ConfigOrganism {
                taxonid: 4932,
                genus: RcString::from("Saccharomyces"),
                species: RcString::from("cerevisiae"),
                assembly_version: None,
            }
        ],
        api_seq_chunk_sizes: vec![10_000, 200_000],
        extension_display_names: vec![],
        extension_relation_order: RelationOrder{
            relation_order: vec![
                RcString::from("directly_positively_regulates"),
                RcString::from("has_direct_input"),
                RcString::from("involved_in"),
                RcString::from("occurs_at"),
                RcString::from("occurs_in"),
                RcString::from("added_by"),
                RcString::from("added_during"),
                RcString::from("has_penetrance"),
            ],
            always_last: vec![RcString::from("happens_during"),
                              RcString::from("exists_during")],
        },
        evidence_types: HashMap::new(),
        cv_config: HashMap::new(),
        interesting_parents: vec![],
        viability_terms: ViabilityTerms {
            viable: RcString::from("FYPO:0002058"),
            inviable: RcString::from("FYPO:0002059"),
        },
        reference_page_config: ReferencePageConfig {
            triage_status_to_ignore: vec![],
        },
        slims: HashMap::new(),
        interpro: InterPro {
            dbnames_to_filter: vec![],
        },
        server: ServerConfig {
            subsets: ServerSubsetConfig {
                prefixes_to_remove: vec![],
            },
            solr_url: String::from("http://localhost:8983/solr"),
            close_synonym_boost: 0.6,
            distant_synonym_boost: 0.3,
            django_url: String::from("http://localhost:8999"),
        },
        extra_database_aliases: HashMap::new(),
        chromosomes: HashMap::new(),
        gene_results: GeneResultsConfig {
            visualisation: GeneResultVisConfig {
                columns: vec![],
            }
        },
        query_data_config: QueryDataConfig {
            ortholog_presence_taxonids: HashSet::from_iter(vec![9606, 4932]),
        },
        file_exports: FileExportConfig {
            macromolecular_complexes: None,
            rnacentral: None,
        },
    };

    config.cv_config.insert(RcString::from("molecular_function"),
                            CvConfig {
                                feature_type: RcString::from("Gene"),
                                filters: vec![],
                                split_by_parents: vec![],
                                summary_relations_to_hide: vec![],
                                summary_relation_ranges_to_collect: vec![RcString::from("has_substrate")],
                                sort_details_by: None,
                                source_config: HashMap::new(),
                            });

    config
}


#[test]
fn test_compare_ext_part_with_config() {
    let config = get_test_config();
    let mut ext_part1 = ExtPart {
        rel_type_name: RcString::from("has_direct_input"),
        rel_type_display_name: RcString::from("NA"),
        ext_range: ExtRange::Misc(RcString::from("misc_ext_part_1")),
    };
    let mut ext_part2 = ExtPart {
        rel_type_name: RcString::from("has_direct_input"),
        rel_type_display_name: RcString::from("NA"),
        ext_range: ExtRange::Misc(RcString::from("misc_ext_part_2")),
    };
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Equal);

    ext_part1.rel_type_name = RcString::from("directly_positively_regulates");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = RcString::from("has_direct_input");
    ext_part2.rel_type_name = RcString::from("directly_positively_regulates");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part2.rel_type_name = RcString::from("absent_during");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part2.rel_type_name = RcString::from("misc_rel");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = RcString::from("other_misc_rel");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = RcString::from("other_misc_rel");
    ext_part2.rel_type_name = RcString::from("other_misc_rel");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Equal);

    ext_part2.rel_type_name = RcString::from("happens_during");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = RcString::from("happens_during");
    ext_part2.rel_type_name = RcString::from("misc_rel");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = RcString::from("has_direct_input");
    ext_part2.rel_type_name = RcString::from("happens_during");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = RcString::from("happens_during");
    ext_part2.rel_type_name = RcString::from("has_direct_input");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = RcString::from("happens_during");
    ext_part2.rel_type_name = RcString::from("exists_during");
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = RcString::from("happens_during");
    ext_part2.rel_type_name = RcString::from("happens_during");
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
    test_conditions.insert(RcString::from("PECO:0000103"));
    test_conditions.insert(RcString::from("PECO:0000137"));

    let mut fypo_details = vec![
        make_one_detail(41_717, "SPBC11B10.09", "PMID:9242669", None,
                        "IDA",vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene(RcString::from("SPBC646.13"))), //  sds23
                        ], HashSet::new()),
        make_one_detail(41_718, "SPBC11B10.09", "PMID:11937031", None,
                        "IDA", vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene(RcString::from("SPBC32F12.09"))), // no name
                        ], HashSet::new()),
        make_one_detail(187_893, "SPBC11B10.09", "PMID:19523829", None, "IMP",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene(RcString::from("SPBC6B1.04"))), //  mde4
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term(RcString::from("GO:1902845"))),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term(RcString::from("GO:0000089"))),
                        ],
                        HashSet::new()),
        make_one_detail(193_221, "SPBC11B10.09", "PMID:10921876", None, "IMP",
                        vec![
                            make_test_ext_part("directly_negatively_regulates", "directly inhibits",
                                               ExtRange::Gene(RcString::from("SPAC144.13c"))), //  srw1
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term(RcString::from("GO:1903693"))),
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term(RcString::from("GO:1905785"))),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term(RcString::from("GO:0000080"))),
                        ],
                        HashSet::new()),
        make_one_detail(194_213, "SPBC11B10.09", "PMID:7957097", None, "IDA",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene(RcString::from("SPBC776.02c"))),  // dis2
                        ],
                        HashSet::new()),
        make_one_detail(194_661, "SPBC11B10.09", "PMID:10485849", None, "IMP",
                        vec![
                            make_test_ext_part("has_direct_input", "has substrate",
                                               ExtRange::Gene(RcString::from("SPBC146.03c"))), //  cut3
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term(RcString::from("GO:1903380"))),
                            make_test_ext_part("part_of", "involved in",
                                               ExtRange::Term(RcString::from("GO:0042307"))),
                            make_test_ext_part("happens_during", "during",
                                               ExtRange::Term(RcString::from("GO:0000089"))),
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
        term: RcString::from("GO:0097472"),
        is_not: false,
        rel_names: HashSet::new(),
        annotations: annotations1,
        summary: None,
    };

    let annotations2 =
        vec![41_717, 41_718, 187_893, 193_221, 194_213, 194_661];

    let ont_term2 = OntTermAnnotations {
        term: RcString::from("GO:0004693"),
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
        genotype: maybe_genotype_uniquename.map(RcString::from),
        genotype_background: None,
        reference: Some(reference_uniquename.into()),
        evidence: Some(evidence.into()),
        withs: HashSet::new(),
        froms: HashSet::new(),
        residue: None,
        qualifiers: vec![],
        extension: extension,
        gene_ex_props: None,
        conditions: conditions,
        assigned_by: Some(RcString::from("PomBase")),
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
        name: name.map(RcString::from),
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
        name: name.map(RcString::from),
        taxonid: 4896,
        product: None,
        deletion_viability: DeletionViability::Unknown,
        uniprot_identifier: None,
        biogrid_interactor_id: None,
        interpro_matches: vec![],
        tm_domain_coords: vec![],
        orfeome_identifier: None,
        name_descriptions: vec![],
        synonyms: vec![],
        dbxrefs: HashSet::new(),
        feature_type: RcString::from("gene"),
        transcript_so_termid: RcString::from("SO:0001217"),
        characterisation_status: None,
        taxonomic_distribution: None,
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
        feature_publications: HashSet::new(),
        subset_termids: HashSet::new(),
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

    ret.insert(RcString::from("e674fe7ceba478aa-genotype-2"),
               make_one_genotype(
                   "G799D(G799D)",
                   Some("test genotype name"),
                   vec![
                       ExpressedAllele {
                           expression: Some(RcString::from("Not assayed")),
                           allele_uniquename: RcString::from("SPBC16A3.11:allele-7"),
                       }
                   ]
               ));

    ret.insert(RcString::from("d6c914796c35e3b5-genotype-4"),
               make_one_genotype(
                   "C-terminal truncation 940-1516(940-1516)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some(RcString::from("Not assayed")),
                           allele_uniquename: RcString::from("SPCC1919.10c:allele-5"),
                       }
                   ]
               ));

    ret.insert(RcString::from("65c76fa511461156-genotype-3"),
               make_one_genotype(
                   "cdc25-22(c532y)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some(RcString::from("Not assayed")),
                           allele_uniquename: RcString::from("SPAC24H6.05:allele-3"),
                       }
                   ]
               ));

    ret.insert(RcString::from("d6c914796c35e3b5-genotype-2"),
               make_one_genotype(
                   "ATPase dead mutant(unknown)",
                   Some("ZZ-name"),
                   vec![
                       ExpressedAllele {
                           expression: Some(RcString::from("Not assayed")),
                           allele_uniquename: RcString::from("SPCC1919.10c:allele-4"),
                       }
                   ]
               ));

    ret.insert(RcString::from("d6c914796c35e3b5-genotype-3"),
               make_one_genotype(
                   "C-terminal truncation(1320-1516)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some(RcString::from("Not assayed")),
                           allele_uniquename: RcString::from("SPCC1919.10c:allele-6"),
                       }
                   ]
               ));

    ret.insert(RcString::from("fd4f3f52f1d38106-genotype-4"),
               make_one_genotype(
                   "K418R(K418R)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some(RcString::from("Wild type product level")),
                           allele_uniquename: RcString::from("SPAC25A8.01c:allele-5"),
                       }
                   ]
               ));

    ret.insert(RcString::from("a6d8f45c20c2227d-genotype-9"),
               make_one_genotype(
                   "UBS-I&II(F18A,F21A,W26A,L40A,W41A,W45A)",
                   None,
                   vec![
                       ExpressedAllele {
                           expression: Some(RcString::from("Not assayed")),
                           allele_uniquename: RcString::from("SPAC3G6.02:allele-7"),
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
        description: description.map(RcString::from),
        name: Some(name.into()),
        allele_type: allele_type.into(),
        gene_uniquename: gene_uniquename.into(),
    }
}

#[allow(dead_code)]
fn get_test_alleles_map() -> UniquenameAlleleMap {
    let mut ret = HashMap::new();

    ret.insert(RcString::from("SPCC1919.10c:allele-4"),
               make_one_allele_short("SPCC1919.10c:allele-4", "ATPase dead mutant", "unknown", None, "SPCC1919.10c"));

    ret.insert(RcString::from("SPCC1919.10c:allele-5"),
               make_one_allele_short("SPCC1919.10c:allele-5", "C-terminal truncation 940-1516", "partial_amino_acid_deletion",
                                     Some("940-1516"), "SPCC1919.10c"));

    ret.insert(RcString::from("SPCC1919.10c:allele-6"),
               make_one_allele_short("SPCC1919.10c:allele-6", "C-terminal truncation", "partial_amino_acid_deletion", Some("1320-1516"),
                                     "SPCC1919.10c"));

    ret.insert(RcString::from("SPBC16A3.11:allele-7"),
               make_one_allele_short("SPBC16A3.11:allele-7", "G799D", "amino_acid_mutation", Some("G799D"), "SPBC16A3.11"));


    ret.insert(RcString::from("SPAC25A8.01c:allele-5"),
               make_one_allele_short("SPAC25A8.01c:allele-5", "K418R", "amino_acid_mutation", Some("K418R"), "SPAC25A8.01c"));

    ret.insert(RcString::from("SPAC3G6.02:allele-7"),
               make_one_allele_short("SPAC3G6.02:allele-7", "UBS-I&II", "amino_acid_mutation", Some("F18A,F21A,W26A,L40A,W41A,W45A"), "SPAC3G6.02"));

    ret.insert(RcString::from("SPAC24H6.05:allele-3"),
               make_one_allele_short("SPAC24H6.05:allele-3", "cdc25-22", "amino_acid_mutation", Some("C532Y"), "SPAC24H6.05"));

    ret
}

#[allow(dead_code)]
fn make_test_term_details(id: &str, name: &str, cv_name: &str) -> TermDetails {
    TermDetails {
        termid: id.into(),
        name: name.into(),
        cv_name: cv_name.into(),
        annotation_feature_type: RcString::from("gene"),
        interesting_parents: HashSet::new(),
        in_subsets: HashSet::from_iter(vec![RcString::from("goslim_pombe")]),
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
        xrefs: HashMap::new(),
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
        "ubs-i_ii-f18a,f21a,w26a,l40a,w41a,w45a-amino_acid_mutation-expression-not_assayed"]
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
            ((*detail).genes[0].to_string(),
             (*detail).reference.clone().unwrap().to_string())
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
        term: RcString::from(termid),
        is_not: false,
        annotations: vec![],
        rel_names: HashSet::new(),
        summary: Some(rows),
    }
}


#[test]
fn test_summary_row_equals() {
    let r1 = TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAPB1A10.09")],
        genotype_uniquenames: vec![],
        extension: vec![],
    };
    let r2 = TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAPB1A10.09")],
        genotype_uniquenames: vec![],
        extension: vec![],
    };
    assert!(r1 == r2);
}

#[allow(dead_code)]
fn get_test_summaries() -> Vec<OntTermAnnotations> {
    let mut summaries = vec![];

    let ext = make_test_ext_part("part_of", "involved in",
                                 ExtRange::Term(RcString::from("GO:1905785")));
    let ext2 = make_test_ext_part("some_rel", "some_rel_display_name",
                                  ExtRange::Term(RcString::from("GO:1234567")));

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
    children_of_0022403.insert(RcString::from("GO:0051318"));
    children_of_0022403.insert(RcString::from("GO:0000080"));
    children_of_0022403.insert(RcString::from("GO:0088888"));
    children_of_0022403.insert(RcString::from("GO:0000089"));
    children_of_0022403.insert(RcString::from("GO:0099999"));

    let mut children_of_0051318 = HashSet::new();
    children_of_0051318.insert(RcString::from("GO:0000080"));
    children_of_0051318.insert(RcString::from("GO:0088888"));
    children_of_0051318.insert(RcString::from("GO:0000089"));
    children_of_0051318.insert(RcString::from("GO:0099999"));

    let mut children_of_0000080 = HashSet::new();
    children_of_0000080.insert(RcString::from("GO:0088888"));

    let mut children_of_0000089 = HashSet::new();
    children_of_0000089.insert(RcString::from("GO:0099999"));

    children_by_termid.insert(RcString::from("GO:0022403"), children_of_0022403);
    children_by_termid.insert(RcString::from("GO:0051318"), children_of_0051318);
    children_by_termid.insert(RcString::from("GO:0000080"), children_of_0000080);
    children_by_termid.insert(RcString::from("GO:0000089"), children_of_0000089);

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

    ret_map.insert(RcString::from("SPAC977.09c"),
                   GeneShort {
                       uniquename: RcString::from("SPAC977.09c"),
                       name: None,
                       product: Some(RcString::from("phospholipase (predicted)")),
                   });
    ret_map.insert(RcString::from("SPAC3G9.09c"),
                   GeneShort {
                       uniquename: RcString::from("SPAC3G9.09c"),
                       name: Some(RcString::from("tif211")),
                       product: Some(RcString::from("translation initiation factor eIF2 alpha subunit")),
                   });
    ret_map.insert(RcString::from("SPAC16.01"),
                   GeneShort {
                       uniquename: RcString::from("SPAC16.01"),
                       name: Some(RcString::from("rho2")),
                       product: Some(RcString::from("Rho family GTPase Rho2")),
                   });
    ret_map.insert(RcString::from("SPAC24C9.02c"),
                   GeneShort { 
                       uniquename: RcString::from("SPAC24C9.02c"),
                       name: Some(RcString::from("cyt2")),
                       product: Some(RcString::from("cytochrome c1 heme lyase Cyt2 (predicted)")),
                   });

    ret_map
}




#[test]
fn test_merge_ext_part_ranges() {
    let ext_part1 = ExtPart {
        rel_type_name: RcString::from("has_substrate"),
        rel_type_display_name: RcString::from("has substrate"),
        ext_range: ExtRange::SummaryGenes(vec![vec![RcString::from("SPAC977.09c")]]),
    };
    let ext_part2 = ExtPart {
        rel_type_name: RcString::from("has_substrate"),
        rel_type_display_name: RcString::from("has substrate"),
        ext_range: ExtRange::SummaryGenes(vec![vec![RcString::from("SPAC24C9.02c")]]),
    };

    let gene_short_map = get_test_gene_short_map();
    let res = merge_ext_part_ranges(&ext_part1, &ext_part2, &gene_short_map);

    assert_eq!(res.ext_range,
               ExtRange::SummaryGenes(vec![vec![RcString::from("SPAC24C9.02c")],
                                           vec![RcString::from("SPAC977.09c")]]));
}

fn get_test_summary_rows() -> Vec<TermSummaryRow> {
    let mut rows = vec![];

    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC25B2.01c")],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC25B2.01c")],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("some_good_rel"),
                rel_type_display_name: RcString::from("some rel"),
                ext_range: ExtRange::Term(RcString::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("has_substrate"),
                rel_type_display_name: RcString::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![RcString::from("SPAC3G9.09c")]]),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("has_substrate"),
                rel_type_display_name: RcString::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![RcString::from("SPAC3G9.09c")]]),
            },
            ExtPart {
                rel_type_name: RcString::from("during"),
                rel_type_display_name: RcString::from("during"),
                ext_range: ExtRange::Term(RcString::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC25B8.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("has_substrate"),
                rel_type_display_name: RcString::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![RcString::from("SPAC3G9.09c")]]),
            },
            ExtPart {
                rel_type_name: RcString::from("during"),
                rel_type_display_name: RcString::from("during"),
                ext_range: ExtRange::Term(RcString::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC1786.03")], // change annotated gene
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("has_substrate"),
                rel_type_display_name: RcString::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![RcString::from("SPAC3G9.09c")]]),
            },
            ExtPart {
                rel_type_name: RcString::from("during"),
                rel_type_display_name: RcString::from("during"),
                ext_range: ExtRange::Term(RcString::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC1786.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("has_substrate"),
                rel_type_display_name: RcString::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![RcString::from("SPAC16.01")]]),   // change substrate
            },
            ExtPart {
                rel_type_name: RcString::from("during"),
                rel_type_display_name: RcString::from("during"),
                ext_range: ExtRange::Term(RcString::from("GO:0070301")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC1786.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("has_substrate"),
                rel_type_display_name: RcString::from("has substrate"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![RcString::from("SPAC16.01")]]),
            },
            ExtPart {
                rel_type_name: RcString::from("during"),
                rel_type_display_name: RcString::from("during"),
                ext_range: ExtRange::Term(RcString::from("GO:0071472")), // change during term
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC222.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("binds"),
                rel_type_display_name: RcString::from("binds"),
                ext_range: ExtRange::Term(RcString::from("PR:000027629")),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC222.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("binds"),
                rel_type_display_name: RcString::from("binds"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![RcString::from("SPAC16.01")]]),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec![RcString::from("SPAC222.03")],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_name: RcString::from("binds"),
                rel_type_display_name: RcString::from("binds"),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec![RcString::from("SPAC16.01")]]),
            },
            ExtPart {
                rel_type_name: RcString::from("binds"),
                rel_type_display_name: RcString::from("binds"),
                ext_range: ExtRange::Term(RcString::from("PR:000027629")),
            },
],
    });

    rows
}

#[test]
fn test_collect_ext_summary_genes() {
    let config = get_test_config();

    let mut rows = get_test_summary_rows();
    assert_eq!(rows.len(), 13);

    let gene_short_map = get_test_gene_short_map();

    pombase::web::cv_summary::collect_ext_summary_genes(&config.cv_config_by_name("molecular_function"),
                                                        &mut rows, &gene_short_map);
    assert_eq!(rows.len(), 11);

    let collected_ext = rows.get(6).unwrap();
    let collected_ext_ext_part_1 = collected_ext.extension.get(0).unwrap();
    let summary_genes_vec = vec![vec![RcString::from("SPAC16.01")],
                                 vec![RcString::from("SPAC3G9.09c")]];
    assert_eq!(collected_ext_ext_part_1.ext_range,
               ExtRange::SummaryGenes(summary_genes_vec));
}

#[test]
fn test_remove_redundant_summary_rows() {
    let mut rows = get_test_summary_rows();
    assert_eq!(rows.len(), 13);

    pombase::web::cv_summary::remove_redundant_summary_rows(&mut rows);
    assert_eq!(rows.len(), 7);
}

