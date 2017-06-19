extern crate pombase;

use std::collections::hash_set::HashSet;
use std::iter::Iterator;

use self::pombase::api::query::*;
use self::pombase::api::result::*;
use self::pombase::api::server_data::*;
use self::pombase::api::query_exec::*;

fn get_server_data() -> ServerData {
    use std::path::PathBuf;
    let mut data_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    data_path.push("tests/test_search_data.json");
    ServerData::new(data_path.to_str().unwrap())
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
fn test_and_or() {
    let qp1 = QueryNode::GeneList(vec!["id_one".into(), "id_two".into(), "id_three".into()]);
    let qp2 = QueryNode::GeneList(vec!["id_two".into(), "id_three".into(), "id_four".into()]);
    let qp3 = QueryNode::GeneList(vec!["id_one".into(), "id_two".into(),
                                       "id_three".into(), "id_four".into()]);
    let qp4 = QueryNode::GeneList(vec!["id_two".into(), "id_three".into(),
                                       "id_four".into(), "id_five".into()]);

    let and_query_node_1 = QueryNode::And(vec![qp1, qp2]);
    let and_query_node_2 = QueryNode::And(vec![qp3, qp4]);

    let and_query =
        Query::from_node(QueryNode::And(vec![and_query_node_1.clone(), and_query_node_2.clone()]));

    check_gene_result(&and_query, vec!["id_two", "id_three"]);

    let or_query =
        Query::from_node(QueryNode::Or(vec![and_query_node_1, and_query_node_2]));

    check_gene_result(&or_query, vec!["id_two", "id_three", "id_four"]);
}

#[test]
fn test_termid() {
    let qp1 = QueryNode::TermId("GO:0044237".into());
    let q1 = Query::from_node(qp1);

    check_gene_result(&q1, vec!["SPAC24B11.06c", "SPAC19G12.03", "SPAC2F3.09", "SPAC27E2.05"]);
}
