use std::collections::hash_set::HashSet;
use std::iter::FromIterator;

use api::server_data::ServerData;
use api::result::*;
use types::GeneUniquename;

#[derive(Serialize, Deserialize, Eq, PartialEq, Clone)]
pub enum QueryNode {
#[serde(rename = "or")]
    Or(Vec<QueryNode>),
#[serde(rename = "and")]
    And(Vec<QueryNode>),
#[serde(rename = "termid")]
    TermId(String),
#[serde(rename = "genelist")]
    GeneList(Vec<GeneUniquename>),
}

fn exec_or(server_data: &ServerData, nodes: &Vec<QueryNode>) -> Result {
    if nodes.len() == 0 {
        return Result {
            status: ResultStatus::Error("illegal query: OR operator has no nodes".into()),
            rows: vec![],
        }
    }

    let mut seen_genes = HashSet::new();
    let mut or_rows = vec![];

    for node in nodes {
        let node_result = node.exec(server_data);
        if node_result.status != ResultStatus::Ok {
            return Result {
                status: node_result.status,
                rows: vec![],
            }
        }
        for node_row in &node_result.rows {
            if !seen_genes.contains(&node_row.gene_uniquename) {
                or_rows.push(ResultRow {
                    gene_uniquename: node_row.gene_uniquename.clone(),
                });
                seen_genes.insert(node_row.gene_uniquename.clone());
            }
        }
    }
    Result {
        status: ResultStatus::Ok,
        rows: or_rows
    }
}

fn exec_and(server_data: &ServerData, nodes: &Vec<QueryNode>) -> Result {
    if nodes.len() == 0 {
        return Result {
            status: ResultStatus::Error("illegal query: AND operator has no nodes".into()),
            rows: vec![],
        };
    }

    let first_node_result = nodes[0].exec(server_data);

    if first_node_result.status != ResultStatus::Ok {
        return Result {
            status: first_node_result.status,
            rows: vec![],
        };
    }

    let current_genes =
        first_node_result.rows.into_iter().map(|row| row.gene_uniquename);

    let mut current_gene_set =
        HashSet::from_iter(current_genes);

    for node in nodes[1..].iter() {
        let node_result = node.exec(server_data);
        if node_result.status != ResultStatus::Ok {
            return Result {
                status: node_result.status,
                rows: vec![],
            }
        }

        let node_genes =
            node_result.rows.into_iter().map(|row| row.gene_uniquename)
            .collect::<HashSet<_>>();

        current_gene_set = current_gene_set.intersection(&node_genes).cloned().collect();
    }
    Result {
        status: ResultStatus::Ok,
        rows: current_gene_set.iter().map(|gene_uniquename| ResultRow {
            gene_uniquename: gene_uniquename.clone()
        }).collect(),
    }
}

fn exec_termid(server_data: &ServerData, term_id: &str) -> Result {
    let rows = server_data.genes_of_termid(term_id).iter()
        .map(|gene_uniquename| ResultRow { gene_uniquename: gene_uniquename.clone() })
        .collect::<Vec<_>>();

    Result {
        status: ResultStatus::Ok,
        rows: rows,
    }
}

fn exec_gene_list(gene_uniquenames: &Vec<GeneUniquename>) -> Result {
    Result {
        status: ResultStatus::Ok,
        rows: gene_uniquenames.iter()
            .map(|gene_uniquename| ResultRow {
                gene_uniquename: gene_uniquename.clone(),
            }).collect(),
    }
}

impl QueryNode {
    pub fn exec(&self, server_data: &ServerData) -> Result {
        use self::QueryNode::*;
        match *self {
            Or(ref nodes) => exec_or(server_data, nodes),
            And(ref nodes) => exec_and(server_data, nodes),
            TermId(ref term_id) => exec_termid(server_data, term_id),
            GeneList(ref gene_list) => exec_gene_list(gene_list),
        }
    }
}


#[derive(Serialize, Deserialize)]
pub struct Query {
    constraints: QueryNode,
}

impl Query {
    pub fn from_node(node: QueryNode) -> Query {
        Query {
            constraints: node
        }
    }

    pub fn exec(&self, server_data: &ServerData) -> Result {
        self.constraints.exec(server_data)
    }
}
