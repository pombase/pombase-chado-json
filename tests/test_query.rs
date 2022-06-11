extern crate pombase;
extern crate flexstr;

use std::iter::Iterator;

use std::collections::HashSet;

use self::pombase::api::query::*;
use self::pombase::api::result::*;
use self::pombase::api::query_exec::*;
use self::pombase::web::config::{Config, TermAndName};
use self::pombase::data_types::{GeneShort, DeletionViability, GeneQueryTermData};
use self::pombase::api_data::*;
use self::pombase::bio::go_format_writer::GO_ASPECT_NAMES;

use flexstr::{ToSharedStr};

fn get_api_data() -> APIData {
    use std::path::PathBuf;
    let mut search_maps_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    search_maps_path.push("tests/test_search_data.json.zst");
    let mut config_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    config_path.push("tests/test_config.json");
    let config = Config::read(&config_path.to_str().expect("config"));
    let api_maps = api_maps_from_file(search_maps_path.to_str().expect("search maps"));
    APIData::new(&config, &api_maps)
}

async fn check_gene_result(query: &Query, genes: Vec<&str>) {
    let api_data = get_api_data();
    let query_exec = QueryExec::new(api_data, None);
    let result = query_exec.exec(&query).await;

    let result_genes_iter =
        result.rows.into_iter()
        .map(|row: ResultRow| row.gene_uniquename)
        .collect::<HashSet<_>>();
    let expect_genes_iter =
        genes.iter().cloned()
        .map(|s: &str| s.to_shared_str())
        .collect::<HashSet<_>>();
    assert_eq!(result_genes_iter,
               expect_genes_iter);
}

async fn check_gene_result_with_viability(query: &Query,
                                    expected_results: &Vec<ResultRow>) {
    let api_data = get_api_data();
    let query_exec = QueryExec::new(api_data, None);
    let results = query_exec.exec(&query).await;

    assert_eq!(expected_results, &results.rows);
}

fn make_genes(ids: Vec<&str>) -> Vec<GeneShort> {
    let mut ret = vec![];
    for id in ids {
        ret.push(GeneShort { uniquename: id.into(), name: None, product: None, transcript_count: 1, });
    }
    ret
}

#[tokio::test]
async fn test_and_or_not() {
    let qp1 = QueryNode {
        gene_list: Some(GeneListNode {
            genes: make_genes(vec!["SPAC19G12.04", "SPAC1805.15c", "SPAC27E2.05"])
        }),
        .. QueryNode::template_node()
    };
    let qp2 = QueryNode {
        gene_list: Some(GeneListNode {
            genes: make_genes(vec!["SPAC1805.15c", "SPAC27E2.05", "SPAC2F3.09"])
        }),
        .. QueryNode::template_node()
    };
    let qp3 = QueryNode {
        gene_list: Some(GeneListNode {
            genes: make_genes(vec!["SPAC19G12.04", "SPAC1805.15c", "SPAC27E2.05", "SPAC2F3.09"]),
        }),
        .. QueryNode::template_node()
    };
    let qp4 = QueryNode {
        gene_list: Some(GeneListNode {
            genes: make_genes(vec!["SPAC1805.15c", "SPAC27E2.05", "SPAC2F3.09", "SPRRNA.26"]),
        }),
        .. QueryNode::template_node()
    };
    let qp5 = QueryNode {
        gene_list: Some(GeneListNode {
            genes: make_genes(vec!["SPAC27E2.05", "SPAC27E2.05", "SPRRNA.26"])
        }),
        .. QueryNode::template_node()
    };
    let and_query_node_1 = QueryNode {
        and: Some(vec![qp1.clone(), qp2.clone()]),
        .. QueryNode::template_node()
    };
    let and_query_node_2 = QueryNode {
        and: Some(vec![qp3.clone(), qp4.clone()]),
        .. QueryNode::template_node()
    };

    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_shared_str()],
        sequence: SeqType::None,
        gaf_options: None,
        ancestor_terms: HashSet::new(),
        flags: HashSet::new(),
    };

    let and_query_node =
        QueryNode {
            and: Some(vec![and_query_node_1.clone(), and_query_node_2.clone()]),
            .. QueryNode::template_node()
        };
    let and_query = Query::new(and_query_node, opts.clone());

    check_gene_result(&and_query, vec!["SPAC1805.15c", "SPAC27E2.05"]).await;

    let or_query_node = QueryNode {
        or: Some(vec![and_query_node_1, and_query_node_2]),
        .. QueryNode::template_node()
    };
    let or_query = Query::new(or_query_node, opts.clone());

    check_gene_result(&or_query, vec!["SPAC1805.15c", "SPAC27E2.05", "SPAC2F3.09"]).await;

    let not_query_node =
        QueryNode {
            not: Some(NotNode {
                node_a: Box::new(qp1.clone()),
                node_b: Box::new(qp2.clone())
            }),
            .. QueryNode::template_node()
        };
    let not_query = Query::new(not_query_node, opts.clone());
    check_gene_result(&not_query, vec!["SPAC19G12.04"]).await;

    let not_query_2_node = QueryNode {
        not: Some(NotNode {
            node_a: Box::new(qp1.clone()),
            node_b: Box::new(qp5),
        }),
        .. QueryNode::template_node()
    };
    let not_query_2 = Query::new(not_query_2_node, opts.clone());
    check_gene_result(&not_query_2, vec!["SPAC19G12.04", "SPAC1805.15c"]).await;
}

#[tokio::test]
async fn test_output_options() {

    // try extra output options
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_shared_str(),
                          "deletion_viability".to_shared_str(),
                          "go_component".to_shared_str()],
        sequence: SeqType::None,
        gaf_options: None,
        ancestor_terms: HashSet::new(),
        flags: HashSet::new(),
    };

    let expected_results =
        vec![
            ResultRow {
                gene_uniquename: "SPAC19G12.04".to_shared_str(),
                deletion_viability: Some(DeletionViability::Inviable),
                go_component: None,
                go_process_superslim: None,
                go_function: None,
                characterisation_status: None,
                taxonomic_distribution: None,
                ortholog_taxonids: HashSet::new(),
                physical_interactors: HashSet::new(),
                tmm: None,
                molecular_weight: None,
                protein_length: None,
                protein_length_bin: None,
                sequence: None,
                gaf_lines: None,
                subsets: HashSet::new(),
                gene_expression: vec![],
            },
            ResultRow {
                gene_uniquename: "SPAC1805.15c".to_shared_str(),
                deletion_viability: Some(DeletionViability::Viable),
                go_component: Some(GeneQueryTermData::Other),
                go_process_superslim: None,
                go_function: None,
                characterisation_status: None,
                taxonomic_distribution: None,
                ortholog_taxonids: HashSet::new(),
                physical_interactors: HashSet::new(),
                molecular_weight: None,
                protein_length: None,
                protein_length_bin: None,
                tmm: None,
                sequence: None,
                gaf_lines: None,
                subsets: HashSet::new(),
                gene_expression: vec![],
            },
            ResultRow {
                gene_uniquename: "SPAC27E2.05".to_shared_str(),
                deletion_viability: Some(DeletionViability::DependsOnConditions),
                go_component: Some(GeneQueryTermData::Term(TermAndName {
                    termid: "GO:0005634".to_shared_str(),
                    name: "nucleus".to_shared_str(),
                })),
                go_process_superslim: None,
                go_function: None,
                characterisation_status: None,
                taxonomic_distribution: None,
                ortholog_taxonids: HashSet::new(),
                physical_interactors: HashSet::new(),
                molecular_weight: None,
                protein_length: None,
                protein_length_bin: None,
                tmm: None,
                sequence: None,
                gaf_lines: None,
                subsets: HashSet::new(),
                gene_expression: vec![],
            }];

    let qp1 = QueryNode {
        gene_list: Some(GeneListNode {
            genes: make_genes(vec!["SPAC19G12.04", "SPAC1805.15c", "SPAC27E2.05"])
        }),
        .. QueryNode::template_node()
    };

    let query = Query::new(qp1, opts);

    check_gene_result_with_viability(&query, &expected_results).await;
}

#[tokio::test]
async fn test_termid() {
    let qp1 = QueryNode {
        term: Some(TermNode {
            termid: "GO:0044237".into(),
            name: None,
            single_or_multi_locus: None,
            ploidiness: None,
            expression: None,
            conditions: HashSet::new(),
            excluded_conditions: HashSet::new(),
        }),
        .. QueryNode::template_node()
    };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_shared_str()],
        sequence: SeqType::None,
        gaf_options: None,
        ancestor_terms: HashSet::new(),
        flags: HashSet::new(),
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1, vec!["SPAC24B11.06c", "SPAC19G12.03", "SPAC2F3.09", "SPAC27E2.05"]).await;
}

#[tokio::test]
async fn test_gene_subset() {
    let qp1 = QueryNode {
        subset: Some(SubsetNode {
            subset_name: "PTHR17490:SF9".into(),
        }),
        .. QueryNode::template_node()
    };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_shared_str()],
        sequence: SeqType::None,
        gaf_options: None,
        ancestor_terms: HashSet::new(),
        flags: HashSet::new(),
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1, vec!["SPCC895.03c"]).await;
}

#[tokio::test]
async fn test_gene_subset_invert() {
    let qp1 = QueryNode {
        subset: Some(SubsetNode {
            subset_name: "!interpro:IPR002906".into(),
        }),
        .. QueryNode::template_node()
    };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_shared_str()],
        sequence: SeqType::None,
        gaf_options: None,
        ancestor_terms: HashSet::new(),
        flags: HashSet::new(),
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1,
                      vec!["SPAC589.10c", "SPAC17H9.16",
                           "SPBC460.01c", "SPBC1652.02", "SPCC736.11",
                           "SPBC19F8.06c", "SPAC6G10.11c",
                           "SPBC359.01", "SPAC1039.09", "SPBC359.03c",
                           "SPAC7D4.10", "SPCC895.03c"]).await;
}

#[tokio::test]
async fn test_gene_subset_wildcard() {
    let qp1 = QueryNode {
        subset: Some(SubsetNode {
            subset_name: "interpro:IPR*".into(),
        }),
        .. QueryNode::template_node()
    };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_shared_str()],
        sequence: SeqType::None,
        gaf_options: None,
        ancestor_terms: HashSet::new(),
        flags: HashSet::new(),
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1,
                      vec!["SPAC589.10c", "SPCC736.11", "SPAC6G10.11c"]).await;
}

#[tokio::test]
async fn test_gene_subset_not_wildcard() {
    let qp1 = QueryNode {
        subset: Some(SubsetNode {
            subset_name: "!interpro:IPR*".into(),
        }),
        .. QueryNode::template_node()
    };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_shared_str()],
        sequence: SeqType::None,
        gaf_options: None,
        ancestor_terms: HashSet::new(),
        flags: HashSet::new(),
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1,
                      vec!["SPAC589.10c", "SPAC7D4.10",
                           "SPBC359.01", "SPAC17H9.16", "SPAC1039.09",
                           "SPBC359.03c", "SPCC736.11", "SPBC1652.02",
                           "SPBC460.01c", "SPAC6G10.11c", "SPCC895.03c",
                           "SPBC19F8.06c"]).await;
}

#[tokio::test]
async fn test_gene_gaf() {
    let qp1 = QueryNode {
        term: Some(TermNode {
            termid: "GO:0044237".into(),
            name: None,
            single_or_multi_locus: None,
            ploidiness: None,
            expression: None,
            conditions: HashSet::new(),
            excluded_conditions: HashSet::new(),
        }),
        .. QueryNode::template_node()
    };
    let mut aspects = HashSet::new();
    for aspect in GO_ASPECT_NAMES.iter() {
        aspects.insert(aspect.clone());
    }

    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_shared_str()],
        sequence: SeqType::None,
        gaf_options: Some(GAFOptions {
            aspects,
        }),
        ancestor_terms: HashSet::new(),
        flags: HashSet::new(),
    };
    let query = Query::new(qp1, opts);

    let api_data = get_api_data();
    let query_exec = QueryExec::new(api_data, None);
    let result = query_exec.exec(&query).await;

    let mut gaf_lines = String::new();

    for row in result.rows.into_iter() {
        if let Some(ref lines) = row.gaf_lines {
            gaf_lines += lines;
        }
    }

    assert_eq!(gaf_lines, "PomBase\tSPAC27E2.05\tSPAC27E2.05\t\tGO:0044237\tPB_REF:0000001\tISS\t\tP\t\t\tprotein\ttaxon:4893\t20090213\tPomBase\t\t\n");
}
