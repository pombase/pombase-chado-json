use std::collections::hash_set::HashSet;
use std::iter::FromIterator;

use api::server_data::ServerData;
use api::result::*;
use web::data::APIGeneSummary;

use types::GeneUniquename;

#[derive(Serialize, Deserialize, Eq, PartialEq, Clone)]
pub enum IntRangeType {
#[serde(rename = "genome_range_contains")]
    GenomeRangeContains,
#[serde(rename = "protein_length")]
    ProteinLength,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Clone)]
pub enum FloatRangeType {
#[serde(rename = "protein_mol_weight")]
    ProteinMolWeight,
}

#[derive(Serialize, Deserialize, Clone)]
pub enum QueryNode {
#[serde(rename = "or")]
    Or(Vec<QueryNode>),
#[serde(rename = "and")]
    And(Vec<QueryNode>),
#[serde(rename = "not")]
    Not(Box<QueryNode>, Box<QueryNode>),
#[serde(rename = "termid")]
    TermId(String),
#[serde(rename = "subset")]
    Subset(String),
#[serde(rename = "gene_list")]
    GeneList(Vec<GeneUniquename>),
#[serde(rename = "int_range")]
    IntRange(IntRangeType, Option<u64>, Option<u64>),
#[serde(rename = "float_range")]
    FloatRange(FloatRangeType, Option<f64>, Option<f64>),
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

fn exec_not(server_data: &ServerData, node_a: &QueryNode, node_b: &QueryNode) -> Result {
    let node_b_result = node_b.exec(server_data);
    if node_b_result.status != ResultStatus::Ok {
        return Result {
            status: node_b_result.status,
            rows: vec![],
        }
    }

    let node_b_gene_set: HashSet<GeneUniquename> =
        HashSet::from_iter(node_b_result.rows.into_iter().map(|row| row.gene_uniquename));

    let node_a_result = node_a.exec(server_data);
    if node_a_result.status != ResultStatus::Ok {
        return Result {
            status: node_a_result.status,
            rows: vec![],
        }
    }

    let mut not_rows = vec![];

    for row in &node_a_result.rows {
        if !node_b_gene_set.contains(&row.gene_uniquename) {
            not_rows.push(ResultRow {
                gene_uniquename: row.gene_uniquename.clone(),
            });
        }
    }
    Result {
        status: ResultStatus::Ok,
        rows: not_rows,
    }
}

fn results_from_gene_vec(genes: Vec<GeneUniquename>) -> Result {
    Result {
        status: ResultStatus::Ok,
        rows: genes.into_iter()
            .map(|gene_uniquename| ResultRow {
                gene_uniquename: gene_uniquename,
            }).collect::<Vec<_>>()
    }
}

fn exec_termid(server_data: &ServerData, term_id: &str) -> Result {
    results_from_gene_vec(server_data.genes_of_termid(term_id))
}

fn exec_subset(server_data: &ServerData, subset_name: &str) -> Result {
    results_from_gene_vec(server_data.genes_of_subset(subset_name))
}

fn exec_gene_list(gene_uniquenames: &Vec<GeneUniquename>) -> Result {
    results_from_gene_vec(gene_uniquenames.clone())
}

fn exec_genome_range_overlaps(server_data: &ServerData,
                              range_start: Option<u64>, range_end: Option<u64>)
                              -> Result
{
    let gene_uniquenames =
        server_data.filter_genes(&|gene: &APIGeneSummary| {
            if let Some(ref location) = gene.location {
                (range_end.is_none() || location.start_pos as u64 <= range_end.unwrap()) &&
                (range_start.is_none() || location.end_pos as u64 >= range_start.unwrap())
            } else {
                false
            }
        });
    results_from_gene_vec(gene_uniquenames)
}

fn exec_protein_length_range(server_data: &ServerData,
                             range_start: Option<u64>, range_end: Option<u64>)
                             -> Result
{
    let gene_uniquenames =
        server_data.filter_genes(&|gene: &APIGeneSummary| {
            if gene.transcripts.len() > 0 {
                if let Some(ref protein) = gene.transcripts[0].protein {
                    (range_start.is_none() || protein.sequence.len() as u64 >= range_start.unwrap()) &&
                    (range_end.is_none() || protein.sequence.len() as u64 <= range_end.unwrap())
                } else {
                    false
                }
            } else {
                false
            }
        });
    results_from_gene_vec(gene_uniquenames)
}

fn exec_int_range(server_data: &ServerData, range_type: &IntRangeType,
                  start: Option<u64>, end: Option<u64>) -> Result {
    match *range_type {
        IntRangeType::GenomeRangeContains => exec_genome_range_overlaps(server_data, start, end),
        IntRangeType::ProteinLength => exec_protein_length_range(server_data, start, end),
    }
}

fn exec_mol_weight_range(server_data: &ServerData, range_start: Option<f64>, range_end: Option<f64>)
                         -> Result
{
    let gene_uniquenames =
        server_data.filter_genes(&|gene: &APIGeneSummary| {
            if gene.transcripts.len() > 0 {
                if let Some(ref protein) = gene.transcripts[0].protein {
                    (range_start.is_none() ||
                        protein.molecular_weight as f64 >= range_start.unwrap()) &&
                    (range_end.is_none() ||
                        protein.molecular_weight as f64 <= range_end.unwrap())
                } else {
                    false
                }
            } else {
                false
            }
        });
    results_from_gene_vec(gene_uniquenames)
}

fn exec_float_range(server_data: &ServerData, range_type: &FloatRangeType,
                    start: Option<f64>, end: Option<f64>) -> Result {
    match *range_type {
        FloatRangeType::ProteinMolWeight => exec_mol_weight_range(server_data, start, end)
    }
}

impl QueryNode {
    pub fn exec(&self, server_data: &ServerData) -> Result {
        use self::QueryNode::*;
        match *self {
            Or(ref nodes) => exec_or(server_data, nodes),
            And(ref nodes) => exec_and(server_data, nodes),
            Not(ref node_a, ref node_b) => exec_not(server_data, node_a, node_b),
            TermId(ref term_id) => exec_termid(server_data, term_id),
            Subset(ref subset_name) => exec_subset(server_data, subset_name),
            GeneList(ref gene_list) => exec_gene_list(gene_list),
            IntRange(ref range_type, start, end) =>
                exec_int_range(server_data, range_type, start, end),
            FloatRange(ref range_type, start, end) =>
                exec_float_range(server_data, range_type, start, end),
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
