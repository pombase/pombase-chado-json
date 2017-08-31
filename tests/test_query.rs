extern crate pombase;

use std::collections::hash_set::HashSet;
use std::iter::Iterator;

use self::pombase::api::query::*;
use self::pombase::api::result::*;
use self::pombase::api::server_data::*;
use self::pombase::api::query_exec::*;

fn get_server_data() -> ServerData {
    use std::path::PathBuf;
    let mut search_maps_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    search_maps_path.push("tests/test_search_data.json.xz");
    let mut gene_subsets_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    gene_subsets_path.push("tests/test_gene_subsets.json");
    let mut config_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    config_path.push("tests/test_config.json");
    ServerData::new(config_path.to_str().unwrap(), search_maps_path.to_str().unwrap(),
                    gene_subsets_path.to_str().unwrap())
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
        .map(|s: &str| String::from(s))
        .collect::<HashSet<_>>();
    assert_eq!(result_genes_iter,
               expect_genes_iter);
}

#[test]
fn test_and_or_not() {
    let qp1 = QueryNode::GeneList { ids: vec!["id_one".into(), "id_two".into(), "id_three".into()] };
    let qp2 = QueryNode::GeneList { ids: vec!["id_two".into(), "id_three".into(), "id_four".into()] };
    let qp3 = QueryNode::GeneList { ids: vec!["id_one".into(), "id_two".into(),
                                              "id_three".into(), "id_four".into()] };
    let qp4 = QueryNode::GeneList { ids: vec!["id_two".into(), "id_three".into(),
                                              "id_four".into(), "id_five".into()] };
    let qp5 = QueryNode::GeneList { ids: vec!["id_three".into(), "id_four".into(),
                                              "id_five".into()] };

    let and_query_node_1 = QueryNode::And(vec![qp1.clone(), qp2.clone()]);
    let and_query_node_2 = QueryNode::And(vec![qp3.clone(), qp4.clone()]);

    let opts = QueryOutputOptions {
        field_names: vec!["gene_uniquename".to_owned()],
        sequence: SeqType::None,
    };

    let and_query =
        Query::new(QueryNode::And(vec![and_query_node_1.clone(), and_query_node_2.clone()]),
                   opts.clone());

    check_gene_result(&and_query, vec!["id_two", "id_three"]);

    let or_query =
        Query::new(QueryNode::Or(vec![and_query_node_1, and_query_node_2]), opts.clone());

    check_gene_result(&or_query, vec!["id_two", "id_three", "id_four"]);

    let not_query_node =
        QueryNode::Not { node_a: Box::new(qp1.clone()), node_b: Box::new(qp2.clone()) };
    let not_query = Query::new(not_query_node, opts.clone());
    check_gene_result(&not_query, vec!["id_one"]);

    let not_query_2_node = QueryNode::Not { node_a: Box::new(qp1), node_b: Box::new(qp5) };
    let not_query_2 = Query::new(not_query_2_node, opts.clone());
    check_gene_result(&not_query_2, vec!["id_one", "id_two"]);
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
