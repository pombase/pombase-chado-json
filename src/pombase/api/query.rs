use std::iter::FromIterator;
use std::cmp;

use std::collections::HashSet;

use uuid::Uuid;

use crate::api::server_data::ServerData;
use crate::api::result::*;
use crate::web::data::{APIGeneSummary, TranscriptDetails, FeatureType, GeneShort, InteractionType,
                       ChromosomeDetails, Strand};

use crate::bio::util::rev_comp;

use crate::types::GeneUniquename;

use pombase_rc_string::RcString;

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum IntRangeType {
#[serde(rename = "protein_length")]
    ProteinLength,
#[serde(rename = "tm_domain_count")]
    TMDomainCount,
#[serde(rename = "exon_count")]
    ExonCount,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum FloatRangeType {
#[serde(rename = "protein_mol_weight")]
    ProteinMolWeight,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum SingleOrMultiAllele {
#[serde(rename = "single")]
    Single,
#[serde(rename = "multi")]
    Multi,
#[serde(rename = "both")]
    Both,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum QueryExpressionFilter {
#[serde(rename = "any")]
    Any,
#[serde(rename = "null")]
    Null,
#[serde(rename = "wt-overexpressed")]
    WtOverexpressed,
}

type TermName = String;
type QueryRowsResult = Result<Vec<ResultRow>, RcString>;
type GeneUniquenameVecResult = Result<Vec<GeneUniquename>, RcString>;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum QueryNode {
#[serde(rename = "or")]
    Or(Vec<QueryNode>),
#[serde(rename = "and")]
    And(Vec<QueryNode>),
#[serde(rename = "not")]
    Not { node_a: Box<QueryNode>, node_b: Box<QueryNode> },
#[serde(rename = "term")]
    Term {
        termid: String,
        name: Option<TermName>,
        single_or_multi_allele: Option<SingleOrMultiAllele>,
        expression: Option<QueryExpressionFilter>,
    },
#[serde(rename = "subset")]
    Subset { subset_name: String },
#[serde(rename = "gene_list")]
    GeneList { genes: Vec<GeneShort> },
#[serde(rename = "int_range")]
    IntRange { range_type: IntRangeType, start: Option<usize>, end: Option<usize> },
#[serde(rename = "float_range")]
    FloatRange { range_type: FloatRangeType, start: Option<f64>, end: Option<f64> },
#[serde(rename = "genome_range")]
    GenomeRange { start: Option<usize>, end: Option<usize>, chromosome_name: String, },
#[serde(rename = "interactors")]
    Interactors { gene_uniquename: GeneUniquename, interaction_type: String },
#[serde(rename = "query_id")]
    QueryId { id: Uuid },
}

fn exec_or(server_data: &ServerData, nodes: &[QueryNode]) -> GeneUniquenameVecResult {
    if nodes.is_empty() {
        return Err(RcString::from("illegal query: OR operator has no nodes"));
    }

    let mut seen_genes = HashSet::new();
    let mut or_rows = vec![];

    for node in nodes {
        let exec_rows = node.exec(server_data)?;

        for row_gene_uniquename in &exec_rows {
            if !seen_genes.contains(row_gene_uniquename) {
                or_rows.push(row_gene_uniquename.clone());
                seen_genes.insert(row_gene_uniquename.clone());
            }
        }
    }

    Ok(or_rows)
}

fn exec_and(server_data: &ServerData, nodes: &[QueryNode]) -> GeneUniquenameVecResult {
    if nodes.is_empty() {
        return Err("illegal query: AND operator has no nodes".into());
    }

    let first_node_genes = nodes[0].exec(server_data)?;

    let current_genes = first_node_genes;

    let mut current_gene_set = HashSet::from_iter(current_genes);

    for node in nodes[1..].iter() {
        let node_result_rows = node.exec(server_data)?;
        let node_genes = node_result_rows.into_iter().collect::<HashSet<_>>();

        current_gene_set = current_gene_set.intersection(&node_genes).cloned().collect();
    }

    Ok(current_gene_set.into_iter().collect())
}

fn exec_not(server_data: &ServerData, node_a: &QueryNode, node_b: &QueryNode)
             -> GeneUniquenameVecResult
{
    let node_b_result = node_b.exec(server_data)?;

    let node_b_gene_set: HashSet<GeneUniquename> =
        HashSet::from_iter(node_b_result.into_iter());

    let node_a_result = node_a.exec(server_data)?;

    let mut not_rows = vec![];

    for row_gene_uniquename in &node_a_result {
        if !node_b_gene_set.contains(row_gene_uniquename) {
            not_rows.push(row_gene_uniquename.clone());
        }
    }

    Ok(not_rows)
}

fn exec_termid(server_data: &ServerData, term_id: &str,
               maybe_single_or_multi_allele: &Option<SingleOrMultiAllele>,
               expression: &Option<QueryExpressionFilter>)  -> GeneUniquenameVecResult {
    if let Some(ref single_or_multi_allele) = *maybe_single_or_multi_allele {
        let genes = server_data.genes_of_genotypes(term_id, single_or_multi_allele,
                                                   expression);
        Ok(genes)
    } else {
        Ok(server_data.genes_of_termid(term_id))
    }
}

fn exec_subset(server_data: &ServerData, subset_name: &str)  -> GeneUniquenameVecResult {
    Ok(server_data.genes_of_subset(subset_name))
}

fn exec_gene_list(genes: &[GeneShort])
                  -> GeneUniquenameVecResult
{
    let mut ret = vec![];

    for gene in genes {
        ret.push(gene.uniquename.clone());
    }

    Ok(ret)
}

fn exec_genome_range_overlaps(server_data: &ServerData, start: Option<usize>, end: Option<usize>,
                              chromosome_name: &str) -> GeneUniquenameVecResult {
    let gene_uniquenames =
        server_data.filter_genes(&|gene: &APIGeneSummary| {
            if let Some(ref location) = gene.location {
                (end.is_none() || usize::from(location.start_pos) <= end.unwrap()) &&
                    (start.is_none() || usize::from(location.end_pos) >= start.unwrap()) &&
                    location.chromosome_name == chromosome_name
            } else {
                false
            }
        });
    Ok(gene_uniquenames)
}

fn exec_protein_length_range(server_data: &ServerData,
                             range_start: Option<usize>, range_end: Option<usize>)
                              -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        server_data.filter_genes(&|gene: &APIGeneSummary| {
            if !gene.transcripts.is_empty() {
                if let Some(ref protein) = gene.transcripts[0].protein {
                    (range_start.is_none() || protein.sequence.len() >= range_start.unwrap()) &&
                    (range_end.is_none() || protein.sequence.len() <= range_end.unwrap())
                } else {
                    false
                }
            } else {
                false
            }
        });
    Ok(gene_uniquenames)
}

fn exec_tm_domain_count_range(server_data: &ServerData,
                              range_start: Option<usize>, range_end: Option<usize>)
                               -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        server_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.tm_domain_count >= range_start.unwrap()) &&
            (range_end.is_none() || gene.tm_domain_count <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_exon_count_range(server_data: &ServerData,
                         range_start: Option<usize>, range_end: Option<usize>)
                         -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        server_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.exon_count >= range_start.unwrap()) &&
            (range_end.is_none() || gene.exon_count <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_int_range(server_data: &ServerData, range_type: &IntRangeType,
                  start: Option<usize>, end: Option<usize>) -> GeneUniquenameVecResult {
    match *range_type {
        IntRangeType::ProteinLength => exec_protein_length_range(server_data, start, end),
        IntRangeType::TMDomainCount => exec_tm_domain_count_range(server_data, start, end),
        IntRangeType::ExonCount => exec_exon_count_range(server_data, start, end),
    }
}

fn exec_mol_weight_range(server_data: &ServerData, range_start: Option<f64>, range_end: Option<f64>)
                         -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        server_data.filter_genes(&|gene: &APIGeneSummary| {
            if !gene.transcripts.is_empty() {
                if let Some(ref protein) = gene.transcripts[0].protein {
                    (range_start.is_none() ||
                        f64::from(protein.molecular_weight) >= range_start.unwrap()) &&
                    (range_end.is_none() ||
                        f64::from(protein.molecular_weight) <= range_end.unwrap())
                } else {
                    false
                }
            } else {
                false
            }
        });
    Ok(gene_uniquenames)
}

fn exec_float_range(server_data: &ServerData, range_type: &FloatRangeType,
                    start: Option<f64>, end: Option<f64>) -> GeneUniquenameVecResult {
    match *range_type {
        FloatRangeType::ProteinMolWeight => exec_mol_weight_range(server_data, start, end)
    }
}

fn exec_interactors_of_gene(server_data: &ServerData, gene_uniquename: &GeneUniquename,
                            interaction_type: InteractionType) -> GeneUniquenameVecResult {
    Ok(server_data.interactors_of_genes(gene_uniquename, interaction_type))
}

impl QueryNode {
    pub fn exec(&self, server_data: &ServerData) -> GeneUniquenameVecResult {
        use self::QueryNode::*;
        match *self {
            Or(ref nodes) => exec_or(server_data, nodes),
            And(ref nodes) => exec_and(server_data, nodes),
            Not { ref node_a, ref node_b } => exec_not(server_data, node_a, node_b),
            Term {
                ref termid,
                ref single_or_multi_allele,
                ref expression,
                ..
            } => exec_termid(server_data, termid, single_or_multi_allele, expression),
            Subset { ref subset_name } => exec_subset(server_data, subset_name),
            GeneList { ref genes } => exec_gene_list(genes),
            GenomeRange { start, end, ref chromosome_name } =>
                exec_genome_range_overlaps(server_data, start, end, chromosome_name),
            Interactors { ref gene_uniquename, ref interaction_type } => {
                match &interaction_type as &str {
                    "physical" =>
                        exec_interactors_of_gene(server_data, gene_uniquename,
                                                 InteractionType::Physical),
                    "genetic" =>
                        exec_interactors_of_gene(server_data, gene_uniquename,
                                                 InteractionType::Genetic),
                    _ => Err(RcString::from(&format!("No such interaction type: {}",
                                                     interaction_type)))
                }
            },
            IntRange { ref range_type, start, end } =>
                exec_int_range(server_data, range_type, start, end),
            FloatRange { ref range_type, start, end } =>
                exec_float_range(server_data, range_type, start, end),
            QueryId { id } => panic!("This case shouldn't be reached - id: {}", id),
       }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct NucleotideDownloadOptions {
    include_introns: bool,
    include_5_prime_utr: bool,
    include_3_prime_utr: bool,
    upstream_bases: Option<usize>,
    downstream_bases: Option<usize>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum SeqType {
#[serde(rename = "protein")]
    Protein,
#[serde(rename = "nucleotide")]
    Nucleotide(NucleotideDownloadOptions),
#[serde(rename = "none")]
    None,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct QueryOutputOptions {
    pub sequence: SeqType,
    pub field_names: Vec<String>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub flags: HashSet<String>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Query {
    output_options: QueryOutputOptions,
    constraints: QueryNode,
}

#[derive(PartialEq)]
enum BeforeOrAfter {
    Before,
    After,
}

fn get_chr_range(chr_residues: &RcString, feature_edge: usize, base_count: usize,
                 before_or_after: BeforeOrAfter) -> RcString
{
    print!("{:?} {:?}\n", feature_edge, base_count);

    let (start_pos, end_pos) =
        if before_or_after == BeforeOrAfter::Before {
            let start_pos =
                if feature_edge < 1 + base_count {
                    // prevent addition underflowing
                    0
                } else {
                    cmp::max(0, feature_edge - 1 - base_count)
                };
            (start_pos, feature_edge - 1)
        } else {
            let end_pos =
                if base_count > chr_residues.len() {
                    // prevent addition overflowing
                    chr_residues.len()
                } else {
                    cmp::min(chr_residues.len(), feature_edge + base_count)
                };
            (feature_edge, end_pos)
        };

    print!("{:?} {:?}\n", start_pos, end_pos);

    RcString::from(&chr_residues.as_str()[start_pos..end_pos])
}

impl Query {
    pub fn new(constraints: QueryNode, output_options: QueryOutputOptions) -> Query {
        Query {
            output_options,
            constraints
        }
    }

    fn make_nucl_seq(&self, transcript: &TranscriptDetails,
                     maybe_chr_details: Option<&ChromosomeDetails>,
                     options: &NucleotideDownloadOptions) -> RcString {
        let mut seq = String::from("");

        let loc = &transcript.location;
        let cds_loc = &transcript.cds_location;

        if let (Some(chr_details), Some(upstream_bases)) =
            (maybe_chr_details, options.upstream_bases) {
                let (before_or_after, feature_edge) =
                    if options.include_5_prime_utr || cds_loc.is_none() {
                        if loc.strand == Strand::Forward {
                            (BeforeOrAfter::Before, loc.start_pos)
                        } else {
                            (BeforeOrAfter::After, loc.end_pos)
                        }
                    } else {
                        if loc.strand == Strand::Forward {
                            (BeforeOrAfter::Before, cds_loc.as_ref().unwrap().start_pos)
                        } else {
                            (BeforeOrAfter::After, cds_loc.as_ref().unwrap().end_pos)
                        }
                    };

                let range_seq = get_chr_range(&chr_details.residues, feature_edge,
                                              upstream_bases, before_or_after);

                if loc.strand == Strand::Forward {
                    seq += range_seq.as_str();
                } else {
                    seq += rev_comp(range_seq.as_str()).as_str();
                }
            }

        for part in &transcript.parts {
            if part.feature_type == FeatureType::Exon ||
                options.include_introns && part.feature_type == FeatureType::CdsIntron ||
                options.include_5_prime_utr &&
                (part.feature_type == FeatureType::FivePrimeUtr ||
                 options.include_introns && part.feature_type == FeatureType::FivePrimeUtrIntron) ||
                options.include_3_prime_utr &&
                (part.feature_type == FeatureType::ThreePrimeUtr ||
                 options.include_introns && part.feature_type == FeatureType::ThreePrimeUtrIntron) {
                    seq += &part.residues;
                }
        }

        if let (Some(chr_details), Some(downstream_bases)) =
            (maybe_chr_details, options.downstream_bases) {
                let (before_or_after, feature_edge) =
                    if options.include_3_prime_utr || cds_loc.is_none() {
                        if loc.strand == Strand::Forward {
                            (BeforeOrAfter::After, loc.end_pos)
                        } else {
                            (BeforeOrAfter::Before, loc.start_pos)
                        }
                    } else {
                        if loc.strand == Strand::Forward {
                            (BeforeOrAfter::After, cds_loc.as_ref().unwrap().end_pos)
                        } else {
                            (BeforeOrAfter::Before, cds_loc.as_ref().unwrap().start_pos)
                        }
                    };

                let range_seq = get_chr_range(&chr_details.residues, feature_edge,
                                              downstream_bases, before_or_after);

                if loc.strand == Strand::Forward {
                    seq += range_seq.as_str();
                } else {
                    seq += rev_comp(range_seq.as_str()).as_str();
                }
            }

        RcString::from(&seq)
    }

    fn make_sequence(&self, server_data: &ServerData,
                     gene_uniquename: &RcString) -> Option<RcString> {
        let maybe_gene_summary = server_data.get_gene_summary(gene_uniquename);

        if let Some(gene_summary) = maybe_gene_summary {
            let maybe_transcript = gene_summary.transcripts.get(0);
            if let Some(transcript) = maybe_transcript {
                match self.output_options.sequence {
                    SeqType::Protein =>
                        if let Some(ref protein) = transcript.protein {
                            return Some(protein.sequence.clone());
                        }
                    SeqType::Nucleotide(ref options) => {
                        let chr = gene_summary.location.as_ref()
                            .map(|ref l| -> &RcString { &l.chromosome_name } )
                            .map_or(None, |ref name| server_data.get_chr_details(name.as_str()));
                        let seq = self.make_nucl_seq(transcript, chr, options);
                        return Some(seq)
                    }
                    SeqType::None => (),
                }
            }
        }

        None
    }

    fn make_result_rows(&self, server_data: &ServerData,
                        genes: Vec<RcString>) -> QueryRowsResult {
        Ok(genes.into_iter()
           .map(|gene_uniquename| {
               let mut deletion_viability = None;
               let mut go_component = None;
               let mut go_process_superslim = None;
               let mut go_function = None;
               let mut characterisation_status = None;
               let mut taxonomic_distribution = None;
               let mut ortholog_taxonids = HashSet::new();
               let mut tmm = None;
               let mut protein_length_bin = None;
               let mut subsets = HashSet::new();

               let maybe_gene_data = server_data.get_gene_query_data(&gene_uniquename);

               if let Some(gene_data) = maybe_gene_data {
                   for field_name in &self.output_options.field_names {
                       match field_name as &str {
                           "deletion_viability" =>
                                deletion_viability = Some(gene_data.deletion_viability.clone()),
                           "ortholog_taxonids" =>
                                ortholog_taxonids = gene_data.ortholog_taxonids.clone(),
                           "go_component" =>
                                go_component = gene_data.go_component.clone(),
                           "go_process_superslim" =>
                                go_process_superslim = gene_data.go_process_superslim.clone(),
                           "characterisation_status" =>
                               characterisation_status = gene_data.characterisation_status.clone(),
                           "taxonomic_distribution" =>
                               taxonomic_distribution = gene_data.taxonomic_distribution.clone(),
                           "go_function" =>
                                go_function = gene_data.go_function.clone(),
                           "tmm" => tmm = gene_data.tmm.clone(),
                           "protein_length_bin" =>
                               protein_length_bin = gene_data.protein_length_bin.clone(),
                           "gene_uniquename" => (),
                           _ => eprintln!("warning - no such option field: {}", field_name),
                       }
                   }


                   subsets =
                       if self.output_options.flags.contains("include_gene_subsets") {
                           gene_data.subset_termids.clone()
                       } else {
                           HashSet::new()
                       };
               }

               let sequence = self.make_sequence(server_data, &gene_uniquename);
               ResultRow {
                   sequence,
                   deletion_viability,
                   go_component,
                   go_process_superslim,
                   go_function,
                   characterisation_status,
                   taxonomic_distribution,
                   ortholog_taxonids,
                   tmm,
                   gene_uniquename,
                   protein_length_bin,
                   subsets,
               }
           }).collect::<Vec<_>>())
    }

    pub fn exec(&self, server_data: &ServerData) -> QueryRowsResult {
        let genes_result = self.constraints.exec(server_data);

        match genes_result {
            Ok(genes) => self.make_result_rows(server_data, genes),
            Err(err) => Err(err)
        }
    }
}
