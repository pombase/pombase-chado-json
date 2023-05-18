use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;

use std::cmp::Ordering;

extern crate pombase;
extern crate flexstr;

mod util;

use crate::util::{setup_test_maps_database, get_test_genes_map, get_test_genotypes_map, get_test_alleles_map, get_test_terms_map, get_test_references_map, get_api_data};

use self::pombase::types::*;
use self::pombase::data_types::*;
use self::pombase::web::config::*;
use self::pombase::web::cv_summary::*;

use flexstr::{ToSharedStr, shared_str as flex_str};

use rusqlite::Connection;

#[allow(dead_code)]
fn get_test_config() -> Config {
    let mut config = Config {
        database_name: "PomBase".into(),
        database_long_name: "PomBase".to_shared_str(),
        database_citation: "PMID:22039153".to_shared_str(),
        funder: "Wellcome Trust".to_shared_str(),
        site_description: "PomBase".to_shared_str(),
        load_organism_taxonid: Some(4896),
        base_url: "https://www.pombase.org".to_shared_str(),
        helpdesk_address: "helpdesk@some.domain".to_shared_str(),
        doc_page_aliases: HashMap::new(),
        sequence_feature_page: SeqFeaturePageConfig {
            so_types_to_show: vec![flex_str!("regional_centromere")],
        },
        organisms: vec![
            ConfigOrganism {
                taxonid: 4896,
                genus: "Schizosaccharomyces".to_shared_str(),
                species: "pombe".to_shared_str(),
                alternative_names: vec![],
                assembly_version: Some("ASM294v2".to_shared_str()),
            },
            ConfigOrganism {
                taxonid: 9606,
                genus: "Homo".to_shared_str(),
                species: "sapiens".to_shared_str(),
                alternative_names: vec![],
                assembly_version: None,
            },
            ConfigOrganism {
                taxonid: 4932,
                genus: "Saccharomyces".to_shared_str(),
                species: "cerevisiae".to_shared_str(),
                alternative_names: vec![],
                assembly_version: None,
            }
        ],
        api_seq_chunk_sizes: vec![10_000, 200_000],
        extension_display_names: vec![],
        extension_relation_order: RelationOrder{
            relation_order: vec![
                "directly_positively_regulates".to_shared_str(),
                "has_direct_input".to_shared_str(),
                "involved_in".to_shared_str(),
                "occurs_at".to_shared_str(),
                "occurs_in".to_shared_str(),
                "added_by".to_shared_str(),
                "added_during".to_shared_str(),
                "has_penetrance".to_shared_str(),
            ],
            always_last: vec!["happens_during".to_shared_str(),
                              "exists_during".to_shared_str()],
        },
        evidence_types: HashMap::new(),
        cv_config: HashMap::new(),
        target_of_config: TargetOfConfig {
            relation_priority: HashMap::new(),
        },
        interesting_parents: vec![],
        viability_terms: ViabilityTerms {
            viable: "FYPO:0002058".to_shared_str(),
            inviable: "FYPO:0002059".to_shared_str(),
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
            term_definition_boost: 0.1,
            django_url: String::from("http://localhost:8999"),
            cv_name_for_terms_search: String::from(""),
            gene_uniquename_re: String::from("^SP([ABCN][CTP]|RR|SN|MIT|MT)[\\d\\w]*\\.\\d\\d\\d?\\d?c?$"),
        },
        extra_database_aliases: HashMap::new(),
        chromosomes: vec![],
        gene_results: GeneResultsConfig {
            field_config: HashMap::new(),
            visualisation_field_names: vec![]
        },
        ortholog_taxonids: HashSet::from_iter(vec![9606, 4932]),
        file_exports: FileExportConfig {
            site_map_reference_prefixes: vec![],
            site_map_term_prefixes: vec![],
            macromolecular_complexes: None,
            rnacentral: None,
            annotation_subsets: vec![],
            gpad_gpi: GpadGpiConfig {
                extension_relation_mappings: HashMap::new(),
                transcript_gene_so_term_map: HashMap::new(),
                go_aspect_terms: HashMap::new(),
            },
            nd_reference: String::from("GO_REF:0000015"),
            phaf_cv_name: String::from("single_locus_phenotype"),
            phaf_parental_strain: HashMap::new(),
        },
        gene_expression: GeneExpressionConfig {
            datasets: vec![],
        },
        feature_sub_groups: HashMap::new(),
    };

    config.file_exports.gpad_gpi.go_aspect_terms.insert("molecular_function".to_shared_str(),
                                                        "GO:0003674".to_shared_str());
    config.file_exports.gpad_gpi.go_aspect_terms.insert("cellular_component".to_shared_str(),
                                                        "GO:0005575".to_shared_str());
    config.file_exports.gpad_gpi.go_aspect_terms.insert("biological_process".to_shared_str(),
                                                        "GO:0008150".to_shared_str());

    config.cv_config.insert("molecular_function".to_shared_str(),
                            CvConfig {
                                feature_type: "Gene".to_shared_str(),
                                display_name: Some("molecular function".to_shared_str()),
                                single_or_multi_locus: SingleOrMultiLocusConfig::NotApplicable,
                                filters: vec![],
                                split_by_parents: vec![],
                                summary_relations_to_hide: vec![],
                                summary_relation_ranges_to_collect: vec!["has_substrate".to_shared_str()],
                                sort_details_by: None,
                                source_config: HashMap::new(),
                            });

    config
}


#[test]
fn test_compare_ext_part_with_config() {
    let config = get_test_config();
    let mut ext_part1 = ExtPart {
        rel_type_id: Some("RO:0002400".to_shared_str()),
        rel_type_name: "has_direct_input".to_shared_str(),
        rel_type_display_name: "NA".to_shared_str(),
        ext_range: ExtRange::Misc("misc_ext_part_1".to_shared_str()),
    };
    let mut ext_part2 = ExtPart {
        rel_type_id: Some("RO:0002400".to_shared_str()),
        rel_type_name: "has_direct_input".to_shared_str(),
        rel_type_display_name: "NA".to_shared_str(),
        ext_range: ExtRange::Misc("misc_ext_part_2".to_shared_str()),
    };
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Equal);

    ext_part1.rel_type_name = "directly_positively_regulates".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "has_direct_input".to_shared_str();
    ext_part2.rel_type_name = "directly_positively_regulates".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part2.rel_type_name = "absent_during".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part2.rel_type_name = "misc_rel".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "other_misc_rel".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = "other_misc_rel".to_shared_str();
    ext_part2.rel_type_name = "other_misc_rel".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Equal);

    ext_part2.rel_type_name = "happens_during".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "happens_during".to_shared_str();
    ext_part2.rel_type_name = "misc_rel".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = "has_direct_input".to_shared_str();
    ext_part2.rel_type_name = "happens_during".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "happens_during".to_shared_str();
    ext_part2.rel_type_name = "has_direct_input".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Greater);

    ext_part1.rel_type_name = "happens_during".to_shared_str();
    ext_part2.rel_type_name = "exists_during".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Less);

    ext_part1.rel_type_name = "happens_during".to_shared_str();
    ext_part2.rel_type_name = "happens_during".to_shared_str();
    assert_eq!(pombase::web::data_build::compare_ext_part_with_config(&config.extension_relation_order, &ext_part1, &ext_part2),
               Ordering::Equal);
}

#[allow(dead_code)]
fn make_test_ext_part(rel_type_name: &str, rel_type_display_name: &str,
                      ext_range: ExtRange) -> ExtPart {
    ExtPart {
        rel_type_id: Some("RO:0000000".to_shared_str()),
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

#[allow(dead_code)]
fn get_test_annotations() -> Vec<OntTermAnnotations> {
    let annotations1 = vec![188_448, 202_017];
    let ont_term1 = OntTermAnnotations {
        term: "GO:0097472".to_shared_str(),
        is_not: false,
        rel_names: HashSet::new(),
        annotations: annotations1,
        summary: None,
    };

    let annotations2 =
        vec![41_717, 41_718, 187_893, 193_221, 194_213, 194_661];

    let ont_term2 = OntTermAnnotations {
        term: "GO:0004693".to_shared_str(),
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
        transcript_uniquenames: vec![],
        genotype: maybe_genotype_uniquename.map(|s| s.to_shared_str()),
        genotype_background: None,
        reference: Some(reference_uniquename.into()),
        evidence: Some(evidence.into()),
        eco_evidence: Some(evidence.into()),
        withs: HashSet::new(),
        froms: HashSet::new(),
        date: None,
        residue: None,
        gene_product_form_id: None,
        qualifiers: vec![],
        extension: extension,
        gene_ex_props: None,
        conditions: conditions,
        assigned_by: Some("PomBase".to_shared_str()),
        throughput: Some(Throughput::HighThroughput),
    }
}

#[allow(dead_code)]
fn get_test_fypo_term_details() -> Vec<i32> {
    vec![223_656, 201_099, 201_095, 204_063, 227_452, 201_094, 186_589]
}

#[test]
fn test_cmp_ont_annotation_detail() {
    let mut details_vec = get_test_fypo_term_details();
    let genes = get_test_genes_map();
    let genotypes = get_test_genotypes_map();
    let alleles = get_test_alleles_map();
    let terms = get_test_terms_map();
    let references = get_test_references_map();
    let annotation_details_maps = get_test_annotation_details_map();
    let mut maps_db_conn = Connection::open_in_memory().unwrap();
    setup_test_maps_database(&mut maps_db_conn, &terms, &genes, &references,
                             &genotypes);

    let config = get_test_config();
    let api_data = get_api_data();

    let cmp_detail_with_genotypes =
        |id1: &i32, id2: &i32| {
            let annotation1 = annotation_details_maps.get(id1).expect(&format!("{}", id1));
            let annotation2 = annotation_details_maps.get(id2).expect(&format!("{}", id2));
            let cv_config = &config.cv_config_by_name(&flex_str!("molecular_function"));
            pombase::sort_annotations::cmp_ont_annotation_detail(cv_config,
                                      annotation1, annotation2,
                                      &api_data).unwrap()
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

    let mut result_genotype_display_uniquenames =
        details_vec.drain(0..)
        .map(|detail_id| {
            let detail = annotation_details_maps.get(&detail_id).unwrap();
            let genotype_uniquename = detail.clone().genotype.unwrap();
            let genotype = genotypes.get(&genotype_uniquename).unwrap();
            pombase::web::data_build::make_genotype_display_uniquename(&genotype.loci, &alleles)
                .to_lowercase()
        }).collect::<Vec<String>>();

    result_genotype_display_uniquenames.sort();

    assert_eq!(result_genotype_display_uniquenames, expected);

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
        term: termid.to_shared_str(),
        is_not: false,
        annotations: vec![],
        rel_names: HashSet::new(),
        summary: Some(rows),
    }
}


#[test]
fn test_summary_row_equals() {
    let r1 = TermSummaryRow {
        gene_uniquenames: vec!["SPAPB1A10.09".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![],
    };
    let r2 = TermSummaryRow {
        gene_uniquenames: vec!["SPAPB1A10.09".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![],
    };
    assert!(r1 == r2);
}

#[allow(dead_code)]
fn get_test_summaries() -> Vec<OntTermAnnotations> {
    let mut summaries = vec![];

    let ext = make_test_ext_part("part_of", "involved in",
                                 ExtRange::Term("GO:1905785".to_shared_str()));
    let ext2 = make_test_ext_part("some_rel", "some_rel_display_name",
                                  ExtRange::Term("GO:1234567".to_shared_str()));

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
    children_of_0022403.insert("GO:0051318".to_shared_str());
    children_of_0022403.insert("GO:0000080".to_shared_str());
    children_of_0022403.insert("GO:0088888".to_shared_str());
    children_of_0022403.insert("GO:0000089".to_shared_str());
    children_of_0022403.insert("GO:0099999".to_shared_str());

    let mut children_of_0051318 = HashSet::new();
    children_of_0051318.insert("GO:0000080".to_shared_str());
    children_of_0051318.insert("GO:0088888".to_shared_str());
    children_of_0051318.insert("GO:0000089".to_shared_str());
    children_of_0051318.insert("GO:0099999".to_shared_str());

    let mut children_of_0000080 = HashSet::new();
    children_of_0000080.insert("GO:0088888".to_shared_str());

    let mut children_of_0000089 = HashSet::new();
    children_of_0000089.insert("GO:0099999".to_shared_str());

    children_by_termid.insert("GO:0022403".to_shared_str(), children_of_0022403);
    children_by_termid.insert("GO:0051318".to_shared_str(), children_of_0051318);
    children_by_termid.insert("GO:0000080".to_shared_str(), children_of_0000080);
    children_by_termid.insert("GO:0000089".to_shared_str(), children_of_0000089);

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


#[test]
fn test_merge_ext_part_ranges() {
    let ext_part1 = ExtPart {
        rel_type_id: Some("RO:0000000".to_shared_str()),
        rel_type_name: "has_substrate".to_shared_str(),
        rel_type_display_name: "has substrate".to_shared_str(),
        ext_range: ExtRange::SummaryGenes(vec![vec!["SPAC977.09c".to_shared_str()]]),
    };
    let ext_part2 = ExtPart {
        rel_type_id: Some("RO:0000000".to_shared_str()),
        rel_type_name: "has_substrate".to_shared_str(),
        rel_type_display_name: "has substrate".to_shared_str(),
        ext_range: ExtRange::SummaryGenes(vec![vec!["SPAC24C9.02c".to_shared_str()]]),
    };

    let api_data = get_api_data();

    let res = merge_ext_part_ranges(&ext_part1, &ext_part2, &api_data);

    assert_eq!(res.ext_range,
               ExtRange::SummaryGenes(vec![vec!["SPAC24C9.02c".to_shared_str()],
                                           vec!["SPAC977.09c".to_shared_str()]]));
}

fn get_test_summary_rows() -> Vec<TermSummaryRow> {
    let mut rows = vec![];

    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC25B2.01c".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC25B2.01c".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC25B8.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC25B8.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "some_good_rel".to_shared_str(),
                rel_type_display_name: "some rel".to_shared_str(),
                ext_range: ExtRange::Term("GO:0070301".to_shared_str()),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC25B8.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "has_substrate".to_shared_str(),
                rel_type_display_name: "has substrate".to_shared_str(),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec!["SPAC3G9.09c".to_shared_str()]]),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC25B8.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "has_substrate".to_shared_str(),
                rel_type_display_name: "has substrate".to_shared_str(),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec!["SPAC3G9.09c".to_shared_str()]]),
            },
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "during".to_shared_str(),
                rel_type_display_name: "during".to_shared_str(),
                ext_range: ExtRange::Term("GO:0070301".to_shared_str()),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC25B8.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "has_substrate".to_shared_str(),
                rel_type_display_name: "has substrate".to_shared_str(),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec!["SPAC3G9.09c".to_shared_str()]]),
            },
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "during".to_shared_str(),
                rel_type_display_name: "during".to_shared_str(),
                ext_range: ExtRange::Term("GO:0070301".to_shared_str()),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC1786.03".to_shared_str()], // change annotated gene
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "has_substrate".to_shared_str(),
                rel_type_display_name: "has substrate".to_shared_str(),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec!["SPAC3G9.09c".to_shared_str()]]),
            },
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "during".to_shared_str(),
                rel_type_display_name: "during".to_shared_str(),
                ext_range: ExtRange::Term("GO:0070301".to_shared_str()),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC1786.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "has_substrate".to_shared_str(),
                rel_type_display_name: "has substrate".to_shared_str(),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec!["SPAC16.01".to_shared_str()]]),   // change substrate
            },
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "during".to_shared_str(),
                rel_type_display_name: "during".to_shared_str(),
                ext_range: ExtRange::Term("GO:0070301".to_shared_str()),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC1786.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "has_substrate".to_shared_str(),
                rel_type_display_name: "has substrate".to_shared_str(),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec!["SPAC16.01".to_shared_str()]]),
            },
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "during".to_shared_str(),
                rel_type_display_name: "during".to_shared_str(),
                ext_range: ExtRange::Term("GO:0071472".to_shared_str()), // change during term
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC222.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "binds".to_shared_str(),
                rel_type_display_name: "binds".to_shared_str(),
                ext_range: ExtRange::Term("PR:000027629".to_shared_str()),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC222.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "binds".to_shared_str(),
                rel_type_display_name: "binds".to_shared_str(),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec!["SPAC16.01".to_shared_str()]]),
            }],
    });
    rows.push(TermSummaryRow {
        gene_uniquenames: vec!["SPAC222.03".to_shared_str()],
        genotype_uniquenames: vec![],
        extension: vec![
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "binds".to_shared_str(),
                rel_type_display_name: "binds".to_shared_str(),
                ext_range: ExtRange::SummaryGenes(
                    vec![vec!["SPAC16.01".to_shared_str()]]),
            },
            ExtPart {
                rel_type_id: Some("RO:0000000".to_shared_str()),
                rel_type_name: "binds".to_shared_str(),
                rel_type_display_name: "binds".to_shared_str(),
                ext_range: ExtRange::Term("PR:000027629".to_shared_str()),
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

    let api_data = get_api_data();

    pombase::web::cv_summary::collect_ext_summary_genes(&config.cv_config_by_name(&flex_str!("molecular_function")),
                                                        &mut rows, &api_data);
    assert_eq!(rows.len(), 11);

    let collected_ext = rows.get(6).unwrap();
    let collected_ext_ext_part_1 = collected_ext.extension.get(0).unwrap();
    let summary_genes_vec = vec![vec!["SPAC16.01".to_shared_str()],
                                 vec!["SPAC3G9.09c".to_shared_str()]];
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

