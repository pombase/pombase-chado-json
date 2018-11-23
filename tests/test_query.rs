extern crate pombase;
extern crate pombase_rc_string;

use std::collections::hash_set::HashSet;
use std::iter::Iterator;

use self::pombase::api::query::*;
use self::pombase::api::result::*;
use self::pombase::api::server_data::*;
use self::pombase::api::query_exec::*;
use self::pombase::web::config::TermAndName;
use self::pombase::web::data::{GeneShort, DeletionViability, GeneQueryTermData};

use self::pombase_rc_string::RcString;

fn get_server_data() -> ServerData {
    use std::path::PathBuf;
    let mut search_maps_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    search_maps_path.push("tests/test_search_data.json.gz");
    let mut gene_subsets_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    gene_subsets_path.push("tests/test_gene_subsets.json");
    let mut config_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    config_path.push("tests/test_config.json");
    ServerData::new(config_path.to_str().expect("config"),
                    search_maps_path.to_str().expect("search maps"),
                    gene_subsets_path.to_str().expect("subsets"))
}

fn check_gene_result(query: &Query, genes: Vec<&str>) {
    let server_data = get_server_data();
    let query_exec = QueryExec::new(server_data);
    let result = query_exec.exec(&query);

    let result_genes_iter =
        result.rows.into_iter()
        .map(|row: ResultRow| row.gene_uniquename)
        .collect::<HashSet<_>>();
    let expect_genes_iter =
        genes.iter().cloned()
        .map(|s: &str| RcString::from(s))
        .collect::<HashSet<_>>();
    assert_eq!(result_genes_iter,
               expect_genes_iter);
}

fn check_gene_result_with_viability(query: &Query,
                                    expected_results: &Vec<ResultRow>) {
    let server_data = get_server_data();
    let query_exec = QueryExec::new(server_data);
    let results = query_exec.exec(&query);

    assert_eq!(expected_results, &results.rows);
}

fn make_genes(ids: Vec<&str>) -> Vec<GeneShort> {
    let mut ret = vec![];
    for id in ids {
        ret.push(GeneShort { uniquename: id.into(), name: None, product: None });
    }
    ret
}

#[test]
fn test_and_or_not() {
    let qp1 = QueryNode::GeneList {
        genes: make_genes(vec!["SPAC19G12.04", "SPAC1805.15c", "SPAC27E2.05"])
    };
    let qp2 = QueryNode::GeneList {
        genes: make_genes(vec!["SPAC1805.15c", "SPAC27E2.05", "SPAC2F3.09"])
    };
    let qp3 = QueryNode::GeneList {
        genes: make_genes(vec!["SPAC19G12.04", "SPAC1805.15c", "SPAC27E2.05", "SPAC2F3.09"])
    };
    let qp4 = QueryNode::GeneList {
        genes: make_genes(vec!["SPAC1805.15c", "SPAC27E2.05", "SPAC2F3.09", "SPRRNA.26"])
    };
    let qp5 = QueryNode::GeneList {
        genes: make_genes(vec!["SPAC27E2.05", "SPAC27E2.05", "SPRRNA.26"])
    };
    let and_query_node_1 = QueryNode::And(vec![qp1.clone(), qp2.clone()]);
    let and_query_node_2 = QueryNode::And(vec![qp3.clone(), qp4.clone()]);

    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_owned()],
        sequence: SeqType::None,
    };

    let and_query =
        Query::new(QueryNode::And(vec![and_query_node_1.clone(), and_query_node_2.clone()]),
                   opts.clone());

    check_gene_result(&and_query, vec!["SPAC1805.15c", "SPAC27E2.05"]);

    let or_query =
        Query::new(QueryNode::Or(vec![and_query_node_1, and_query_node_2]), opts.clone());

    check_gene_result(&or_query, vec!["SPAC1805.15c", "SPAC27E2.05", "SPAC2F3.09"]);

    let not_query_node =
        QueryNode::Not { node_a: Box::new(qp1.clone()), node_b: Box::new(qp2.clone()) };
    let not_query = Query::new(not_query_node, opts.clone());
    check_gene_result(&not_query, vec!["SPAC19G12.04"]);

    let not_query_2_node = QueryNode::Not { node_a: Box::new(qp1.clone()), node_b: Box::new(qp5) };
    let not_query_2 = Query::new(not_query_2_node, opts.clone());
    check_gene_result(&not_query_2, vec!["SPAC19G12.04", "SPAC1805.15c"]);
}

#[test]
fn test_output_options() {

    // try extra output options
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_owned(),
                          "deletion_viability".to_owned(),
                          "go_component".to_owned()],
        sequence: SeqType::None,
    };

    let expected_results =
        vec![
            ResultRow {
                gene_uniquename: RcString::from("SPAC19G12.04"),
                deletion_viability: Some(DeletionViability::Inviable),
                go_component: None,
                go_process_superslim: None,
                go_function: None,
                ortholog_taxonids: HashSet::new(),
                tmm: None,
                sequence: None,
            },
            ResultRow {
                gene_uniquename: RcString::from("SPAC1805.15c"),
                deletion_viability: Some(DeletionViability::Viable),
                go_component: Some(GeneQueryTermData::Other),
                go_process_superslim: None,
                go_function: None,
                ortholog_taxonids: HashSet::new(),
                tmm: None,
                sequence: None,
            },
            ResultRow {
                gene_uniquename: RcString::from("SPAC27E2.05"),
                deletion_viability: Some(DeletionViability::DependsOnConditions),
                go_component: Some(GeneQueryTermData::Term(TermAndName {
                    termid: RcString::from("GO:0005634"),
                    name: RcString::from("nucleus"),
                })),
                go_process_superslim: None,
                go_function: None,
                ortholog_taxonids: HashSet::new(),
                tmm: None,
                sequence: None,
            }];

    let qp1 = QueryNode::GeneList {
        genes: make_genes(vec!["SPAC19G12.04", "SPAC1805.15c", "SPAC27E2.05"])
    };

    let query = Query::new(qp1, opts);

    check_gene_result_with_viability(&query, &expected_results);
}

#[test]
fn test_termid() {
    let qp1 = QueryNode::Term { termid: "GO:0044237".into(), name: None,
                                single_or_multi_allele: None, expression: None };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_owned()],
        sequence: SeqType::None,
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1, vec!["SPAC24B11.06c", "SPAC19G12.03", "SPAC2F3.09", "SPAC27E2.05"]);
}

#[test]
fn test_gene_subset() {
    let qp1 = QueryNode::Subset { subset_name: "PTHR17490:SF9".into() };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_owned()],
        sequence: SeqType::None,
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1, vec!["SPCC895.03c"]);
}

#[test]
fn test_gene_subset_invert() {
    let qp1 = QueryNode::Subset { subset_name: "!interpro:IPR002906".into() };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_owned()],
        sequence: SeqType::None,
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1,
                      vec!["SPAC589.10c", "SPAC17H9.16",
                           "SPBC460.01c", "SPBC1652.02", "SPCC736.11",
                           "SPBC19F8.06c", "SPAC6G10.11c",
                           "SPBC359.01", "SPAC1039.09", "SPBC359.03c",
                           "SPAC7D4.10", "SPCC895.03c"]);
}

#[test]
fn test_gene_subset_wildcard() {
    let qp1 = QueryNode::Subset { subset_name: "interpro:IPR*".into() };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_owned()],
        sequence: SeqType::None,
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1,
                      vec!["SPAC589.10c", "SPCC736.11", "SPAC6G10.11c"]);
}

#[test]
fn test_gene_subset_not_wildcard() {
    let qp1 = QueryNode::Subset { subset_name: "!interpro:IPR*".into() };
    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_owned()],
        sequence: SeqType::None,
    };
    let q1 = Query::new(qp1, opts);

    check_gene_result(&q1,
                      vec!["SPAC589.10c", "SPAC7D4.10",
                           "SPBC359.01", "SPAC17H9.16", "SPAC1039.09",
                           "SPBC359.03c", "SPCC736.11", "SPBC1652.02",
                           "SPBC460.01c", "SPAC6G10.11c", "SPCC895.03c",
                           "SPBC19F8.06c"]);
}
