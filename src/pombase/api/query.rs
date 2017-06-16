use std::collections::hash_set::HashSet;
use std::iter::FromIterator;

use api::server_data::ServerData;
use api::result::*;
use types::GeneUniquename;

#[derive(Serialize, Deserialize, Eq, PartialEq, Clone)]
pub enum QueryPart {
#[serde(rename = "or")]
    Or(Vec<QueryPart>),
#[serde(rename = "and")]
    And(Vec<QueryPart>),
#[serde(rename = "termid")]
    TermId(String),
#[serde(rename = "gene_list")]
    GeneList(Vec<GeneUniquename>),
}

fn exec_or(server_data: &ServerData, parts: &Vec<QueryPart>) -> Result {
    if parts.len() == 0 {
        return Result {
            status: ResultStatus::Error("illegal query: OR operator has no parts".into()),
            rows: vec![],
        }
    }

    let mut seen_genes = HashSet::new();
    let mut or_rows = vec![];

    for part in parts {
        let part_result = part.exec(server_data);
        if part_result.status != ResultStatus::Ok {
            return Result {
                status: part_result.status,
                rows: vec![],
            }
        }
        for part_row in &part_result.rows {
            if !seen_genes.contains(&part_row.gene_uniquename) {
                or_rows.push(ResultRow {
                    gene_uniquename: part_row.gene_uniquename.clone(),
                });
                seen_genes.insert(part_row.gene_uniquename.clone());
            }
        }
    }
    Result {
        status: ResultStatus::Ok,
        rows: or_rows
    }
}

fn exec_and(server_data: &ServerData, parts: &Vec<QueryPart>) -> Result {
    if parts.len() == 0 {
        return Result {
            status: ResultStatus::Error("illegal query: AND operator has no parts".into()),
            rows: vec![],
        };
    }

    let first_part_result = parts[0].exec(server_data);

    if first_part_result.status != ResultStatus::Ok {
        return Result {
            status: first_part_result.status,
            rows: vec![],
        };
    }

    let current_genes =
        first_part_result.rows.into_iter().map(|row| row.gene_uniquename);

    let mut current_gene_set =
        HashSet::from_iter(current_genes);

    for part in parts[1..].iter() {
        let part_result = part.exec(server_data);
        if part_result.status != ResultStatus::Ok {
            return Result {
                status: part_result.status,
                rows: vec![],
            }
        }

        let part_genes =
            part_result.rows.into_iter().map(|row| row.gene_uniquename)
            .collect::<HashSet<_>>();

        current_gene_set = current_gene_set.intersection(&part_genes).cloned().collect();
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

impl QueryPart {
    pub fn exec(&self, server_data: &ServerData) -> Result {
        use self::QueryPart::*;
        match *self {
            Or(ref parts) => exec_or(server_data, parts),
            And(ref parts) => exec_and(server_data, parts),
            TermId(ref term_id) => exec_termid(server_data, term_id),
            GeneList(ref gene_list) => exec_gene_list(gene_list),
        }
    }
}


#[derive(Serialize, Deserialize)]
pub struct Query {
    part: QueryPart,
}

impl Query {
    pub fn from_part(part: QueryPart) -> Query {
        Query {
            part: part
        }
    }

    pub fn exec(&self, server_data: &ServerData) -> Result {
        self.part.exec(server_data)
    }
}
