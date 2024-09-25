use std::cmp;
use std::str;

use std::collections::HashSet;

use async_recursion::async_recursion;

use itertools::Itertools;
use uuid::Uuid;

use crate::api_data::APIData;
use crate::api::site_db::SiteDB;
use crate::api::result::*;
use crate::data_types::DataLookup;
use crate::data_types::{APIGeneSummary, TranscriptDetails, FeatureType, GeneShort, InteractionType,
                       ChromosomeDetails, Strand, Ploidiness, GeneQueryPropFlag};
use crate::types::CvName;
use crate::types::TermId;
use crate::web::config::TermAndName;

use crate::bio::util::rev_comp;
use crate::bio::go_format_writer::{write_go_annotation_format, GpadGafWriteMode};


use crate::types::GeneUniquename;

use flexstr::{SharedStr as FlexStr, ToSharedStr, shared_str as flex_str, shared_fmt as flex_fmt};

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum IntRangeType {
#[serde(rename = "spliced_rna_length")]
    SplicedRnaLength,
#[serde(rename = "unspliced_rna_length")]
    UnsplicedRnaLength,
#[serde(rename = "protein_length")]
    ProteinLength,
#[serde(rename = "tm_domain_count")]
    TMDomainCount,
#[serde(rename = "coiled_coils_count")]
    CoiledCoilsCount,
#[serde(rename = "disordered_regions_count")]
    DisorderedRegionsCount,
#[serde(rename = "disordered_percent")]
    DisorderedPercent,
#[serde(rename = "low_complexity_regions_count")]
    LowComplexityRegionsCount,
#[serde(rename = "low_complexity_percent")]
    LowComplexityPercent,
#[serde(rename = "exon_count")]
    ExonCount,
#[serde(rename = "transcript_count")]
    TranscriptCount,
#[serde(rename = "pdb_structure_count")]
    PDBStructureCount,
#[serde(rename = "gocam_model_count")]
    GoCamModelCount,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum FloatRangeType {
#[serde(rename = "protein_mol_weight")]
    ProteinMolWeight,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Copy, Clone)]
pub enum TargetOfType {
#[serde(rename = "go")]
    GO,
#[serde(rename = "phenotype")]
    Phenotype,
#[serde(rename = "all")]
    All,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum SingleOrMultiLocus {
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

type TermName = FlexStr;
type QueryRowsResult = Result<Vec<ResultRow>, FlexStr>;
type GeneUniquenameVecResult = Result<Vec<GeneUniquename>, FlexStr>;

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct NotNode {
    pub node_a: Box<QueryNode>,
    pub node_b: Box<QueryNode>,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct TermNode {
    pub termid: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<TermName>,
    pub single_or_multi_locus: Option<SingleOrMultiLocus>,
    pub ploidiness: Option<Ploidiness>,
    pub expression: Option<QueryExpressionFilter>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub conditions: HashSet<TermAndName>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub excluded_conditions: HashSet<TermAndName>,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct RefGenesNode {
    pub reference_uniquename: FlexStr,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct SubsetNode {
    pub subset_name: FlexStr,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct GeneListNode {
    pub genes: Vec<GeneShort>,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct TargetOfNode {
    pub genes: Vec<GeneShort>,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct HasOrthologNode {
    pub taxonid: u32,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct IntRangeNode {
    pub range_type: IntRangeType,
    pub start: Option<usize>,
    pub end: Option<usize>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub options: HashSet<String>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct FloatRangeNode {
    pub range_type: FloatRangeType,
    pub start: Option<f64>,
    pub end: Option<f64>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub options: HashSet<String>,
}

impl Eq for FloatRangeNode {
}
impl PartialEq for FloatRangeNode {
    fn eq(&self, other: &Self) -> bool {
        if self.range_type != other.range_type {
            return false;
        }
        if self.start.is_none() != other.start.is_none() {
            return false;
        }
        if self.start.is_some() &&
            (self.start.unwrap() - other.start.unwrap()).abs() > 1e-8 {
                return false;
            }
        if self.end.is_some() &&
            (self.end.unwrap() - other.end.unwrap()).abs() > 1e-8 {
                return false;
            }

        true
    }
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct GenomeRangeNode {
    #[serde(skip_serializing_if="Option::is_none")]
    pub start: Option<usize>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub end: Option<usize>,
    pub chromosome_name: String,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct InteractorsNode {
    pub gene_uniquename: GeneUniquename,
    pub interaction_type: String
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct SubstratesNode {
    pub gene_uniquename: GeneUniquename,
    pub phase_term: Option<TermId>,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct DownstreamGenesNode {
    pub gene_uniquename: GeneUniquename,
    pub cv_name: CvName,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct GenesTargetingNode {
    pub gene_uniquename: GeneUniquename,
    pub target_of_type: TargetOfType,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct QueryGenePropNode {
    pub property_flags: HashSet<GeneQueryPropFlag>,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct QueryIdNode {
    pub id: Uuid,
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub struct QueryNode {
    #[serde(skip_serializing_if="Option::is_none")]
    pub node_name: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub or: Option<Vec<QueryNode>>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub and: Option<Vec<QueryNode>>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub not: Option<NotNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub term: Option<TermNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub ref_genes: Option<RefGenesNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub subset: Option<SubsetNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_list: Option<GeneListNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub target_of: Option<TargetOfNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub has_ortholog: Option<HasOrthologNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub int_range: Option<IntRangeNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub float_range: Option<FloatRangeNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genome_range: Option<GenomeRangeNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub interactors: Option<InteractorsNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub substrates: Option<SubstratesNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub downstream_genes: Option<DownstreamGenesNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genes_targeting: Option<GenesTargetingNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_properties: Option<QueryGenePropNode>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub query_id: Option<QueryIdNode>,
}

#[async_recursion]
async fn exec_or(api_data: &APIData, site_db: &Option<SiteDB>,
           nodes: &[QueryNode]) -> GeneUniquenameVecResult {
    if nodes.is_empty() {
        return Err(flex_str!("illegal query: OR operator has no nodes"));
    }

    let mut seen_genes = HashSet::new();
    let mut or_rows = vec![];

    for node in nodes {
        let exec_rows = node.exec(api_data, site_db).await?;

        for row_gene_uniquename in &exec_rows {
            if !seen_genes.contains(row_gene_uniquename) {
                or_rows.push(row_gene_uniquename.clone());
                seen_genes.insert(row_gene_uniquename.clone());
            }
        }
    }

    Ok(or_rows)
}

#[async_recursion]
async fn exec_and(api_data: &APIData, site_db: &Option<SiteDB>,
                  nodes: &[QueryNode]) -> GeneUniquenameVecResult {
    if nodes.is_empty() {
        return Err("illegal query: AND operator has no nodes".into());
    }

    let first_node_genes = nodes[0].exec(api_data, site_db).await?;

    let current_genes = first_node_genes;

    let mut current_gene_set = HashSet::from_iter(current_genes);

    for node in nodes[1..].iter() {
        let node_result_rows = node.exec(api_data, site_db).await?;
        let node_genes = node_result_rows.into_iter().collect::<HashSet<_>>();

        current_gene_set = current_gene_set.intersection(&node_genes).cloned().collect();
    }

    Ok(current_gene_set.into_iter().collect())
}

#[async_recursion]
async fn exec_not(api_data: &APIData, site_db: &Option<SiteDB>,
                  node_a: &QueryNode, node_b: &QueryNode)
                  -> GeneUniquenameVecResult
{
    let node_b_result = node_b.exec(api_data, site_db).await?;

    let node_b_gene_set: HashSet<GeneUniquename> =
        HashSet::from_iter(node_b_result);

    let node_a_result = node_a.exec(api_data, site_db).await?;

    let mut not_rows = vec![];

    for row_gene_uniquename in &node_a_result {
        if !node_b_gene_set.contains(row_gene_uniquename) {
            not_rows.push(row_gene_uniquename.clone());
        }
    }

    Ok(not_rows)
}

fn exec_termid(api_data: &APIData, term_id: &FlexStr,
               maybe_single_or_multi_locus: &Option<SingleOrMultiLocus>,
               maybe_ploidiness: &Option<Ploidiness>,
               expression: &Option<QueryExpressionFilter>,
               conditions: &HashSet<TermAndName>,
               excluded_conditions: &HashSet<TermAndName>)  -> GeneUniquenameVecResult {
    if let Some(ref single_or_multi_locus) = *maybe_single_or_multi_locus {
        let ploidiness = maybe_ploidiness.clone().unwrap_or(Ploidiness::Any);
        let genes = api_data.genes_of_genotypes(term_id, single_or_multi_locus,
                                                &ploidiness,
                                                expression, conditions, excluded_conditions);
        Ok(genes)
    } else {
        Ok(api_data.genes_of_termid(term_id))
    }
}

fn exec_ref_genes(api_data: &APIData, reference_uniquename: &FlexStr) -> GeneUniquenameVecResult {
    if let Some(reference) = api_data.get_reference(reference_uniquename) {
        Ok(reference.genes_by_uniquename.keys().cloned().collect_vec())
    } else {
        Err(format!("no such reference: {}", reference_uniquename).into())
    }
}

fn exec_subset(api_data: &APIData, subset_name: &FlexStr)  -> GeneUniquenameVecResult {
    Ok(api_data.genes_of_subset(subset_name))
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


fn exec_target_of(api_data: &APIData, genes: &[GeneShort])
                  -> GeneUniquenameVecResult
{
    let mut ret = vec![];

    for gene in genes {
        let gene_details = api_data.get_gene_details(&gene.uniquename)
            .unwrap_or_else(|| panic!("failed to find gene_details for {}", gene.uniquename));
        for target_of_annotation in &gene_details.target_of_annotations {
            ret.push(target_of_annotation.gene.clone());
        }
    }

    Ok(ret)
}

fn exec_has_ortholog(api_data: &APIData, taxonid: u32)
                     -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            gene.ortholog_taxonids.contains(&taxonid)
        });

    Ok(gene_uniquenames)
}

fn exec_genome_range_overlaps(api_data: &APIData, start: Option<usize>, end: Option<usize>,
                              chromosome_name: &str) -> GeneUniquenameVecResult {
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            if let Some(ref location) = gene.location {
                (end.is_none() || location.start_pos <= end.unwrap()) &&
                    (start.is_none() || location.end_pos >= start.unwrap()) &&
                    location.chromosome_name == chromosome_name
            } else {
                false
            }
        });
    Ok(gene_uniquenames)
}

fn exec_protein_length_range(api_data: &APIData,
                             range_start: Option<usize>, range_end: Option<usize>)
                              -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            if !gene.transcripts.is_empty() {
                if let Some(ref protein) = gene.transcripts[0].protein {
                    (range_start.is_none() || protein.number_of_residues >= range_start.unwrap()) &&
                    (range_end.is_none() || protein.number_of_residues <= range_end.unwrap())
                } else {
                    false
                }
            } else {
                false
            }
        });
    Ok(gene_uniquenames)
}

fn exec_spliced_rna_length_range(api_data: &APIData,
                                 range_start: Option<usize>, range_end: Option<usize>)
                                 -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            if let Some(first_transcript) = gene.transcripts.get(0) {
                if let Some(spliced_rna_len) = first_transcript.rna_seq_length_spliced {
                    (range_start.is_none() || spliced_rna_len.get() >= range_start.unwrap()) &&
                    (range_end.is_none() || spliced_rna_len.get() <= range_end.unwrap())
                } else {
                    false
                }
            } else {
                false
            }
        });
    Ok(gene_uniquenames)
}

fn exec_unspliced_rna_length_range(api_data: &APIData,
                                 range_start: Option<usize>, range_end: Option<usize>)
                                 -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            if let Some(first_transcript) = gene.transcripts.get(0) {
                if let Some(unspliced_rna_len) = first_transcript.rna_seq_length_unspliced {
                    (range_start.is_none() || unspliced_rna_len.get() >= range_start.unwrap()) &&
                    (range_end.is_none() || unspliced_rna_len.get() <= range_end.unwrap())
                } else {
                    false
                }
            } else {
                false
            }
        });
    Ok(gene_uniquenames)
}

fn exec_tm_domain_count_range(api_data: &APIData,
                              range_start: Option<usize>, range_end: Option<usize>)
                               -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.tm_domain_count >= range_start.unwrap()) &&
                (range_end.is_none() || gene.tm_domain_count <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_coiled_coils_count_range(api_data: &APIData,
                                 range_start: Option<usize>, range_end: Option<usize>)
                                 -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.coiled_coil_count >= range_start.unwrap()) &&
                (range_end.is_none() || gene.coiled_coil_count <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_disordered_regions_count_range(api_data: &APIData,
                                       range_start: Option<usize>, range_end: Option<usize>)
                                       -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.disordered_regions_count >= range_start.unwrap()) &&
                (range_end.is_none() || gene.disordered_regions_count <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_disordered_percent_range(api_data: &APIData,
                                         range_start: Option<usize>, range_end: Option<usize>)
                                         -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.disordered_percent >= range_start.unwrap()) &&
                (range_end.is_none() || gene.disordered_percent <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_low_complexity_regions_count_range(api_data: &APIData,
                                           range_start: Option<usize>, range_end: Option<usize>)
                                           -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.low_complexity_regions_count >= range_start.unwrap()) &&
                (range_end.is_none() || gene.low_complexity_regions_count <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_low_complexity_percent_range(api_data: &APIData,
                                             range_start: Option<usize>, range_end: Option<usize>)
                                             -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.low_complexity_percent >= range_start.unwrap()) &&
                (range_end.is_none() || gene.low_complexity_percent <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_exon_count_range(api_data: &APIData,
                         range_start: Option<usize>, range_end: Option<usize>,
                         options: &HashSet<String>)
                         -> GeneUniquenameVecResult
{

    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            let count =
                if options.contains("five-prime-utr-exons") {
                    gene.five_prime_exon_count
                } else {
                    if options.contains("three-prime-utr-exons") {
                        gene.three_prime_exon_count
                    } else {
                        gene.coding_exon_count
                    }
                };
            (range_start.is_none() || count >= range_start.unwrap()) &&
            (range_end.is_none() || count <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_transcript_count_range(api_data: &APIData,
                               range_start: Option<usize>, range_end: Option<usize>)
                               -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.transcript_count >= range_start.unwrap()) &&
            (range_end.is_none() || gene.transcript_count <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_pdb_id_count_range(api_data: &APIData,
                           range_start: Option<usize>, range_end: Option<usize>)
                       -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.pdb_ids.len() >= range_start.unwrap()) &&
            (range_end.is_none() || gene.pdb_ids.len() <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_gocam_model_count_range(api_data: &APIData,
                                 range_start: Option<usize>, range_end: Option<usize>)
                       -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
            (range_start.is_none() || gene.gocam_ids.len() >= range_start.unwrap()) &&
            (range_end.is_none() || gene.gocam_ids.len() <= range_end.unwrap())
        });
    Ok(gene_uniquenames)
}

fn exec_int_range(api_data: &APIData, range_type: &IntRangeType,
                  start: Option<usize>, end: Option<usize>,
                  options: &HashSet<String>) -> GeneUniquenameVecResult {
    match *range_type {
        IntRangeType::ProteinLength => exec_protein_length_range(api_data, start, end),
        IntRangeType::SplicedRnaLength => exec_spliced_rna_length_range(api_data, start, end),
        IntRangeType::UnsplicedRnaLength => exec_unspliced_rna_length_range(api_data, start, end),
        IntRangeType::TMDomainCount => exec_tm_domain_count_range(api_data, start, end),
        IntRangeType::CoiledCoilsCount => exec_coiled_coils_count_range(api_data, start, end),
        IntRangeType::DisorderedRegionsCount => exec_disordered_regions_count_range(api_data, start, end),
        IntRangeType::DisorderedPercent => exec_disordered_percent_range(api_data, start, end),
        IntRangeType::LowComplexityRegionsCount => exec_low_complexity_regions_count_range(api_data, start, end),
        IntRangeType::LowComplexityPercent => exec_low_complexity_percent_range(api_data, start, end),
        IntRangeType::ExonCount => exec_exon_count_range(api_data, start, end, options),
        IntRangeType::TranscriptCount => exec_transcript_count_range(api_data, start, end),
        IntRangeType::PDBStructureCount => exec_pdb_id_count_range(api_data, start, end),
        IntRangeType::GoCamModelCount => exec_gocam_model_count_range(api_data, start, end),
    }
}

fn exec_mol_weight_range(api_data: &APIData, range_start: Option<f64>, range_end: Option<f64>)
                         -> GeneUniquenameVecResult
{
    let gene_uniquenames =
        api_data.filter_genes(&|gene: &APIGeneSummary| {
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

fn exec_float_range(api_data: &APIData, range_type: &FloatRangeType,
                    start: Option<f64>, end: Option<f64>,
                    _: &HashSet<String>) -> GeneUniquenameVecResult {
    match *range_type {
        FloatRangeType::ProteinMolWeight => exec_mol_weight_range(api_data, start, end)
    }
}

fn exec_interactors_of_gene(api_data: &APIData, gene_uniquename: &GeneUniquename,
                            interaction_type: InteractionType) -> GeneUniquenameVecResult {
    Ok(api_data.interactors_of_genes(gene_uniquename, interaction_type))
}

fn exec_substrates_of_gene(api_data: &APIData, gene_uniquename: &GeneUniquename,
                           phase_term: &Option<TermId>)
    -> GeneUniquenameVecResult
{
    Ok(api_data.downstream_genes("molecular_function", gene_uniquename, phase_term)
        .into_iter().collect())
}

fn exec_downstream_genes(api_data: &APIData, gene_uniquename: &GeneUniquename, cv_name: &str)
    -> GeneUniquenameVecResult
{
    Ok(api_data.downstream_genes(cv_name, gene_uniquename, &None).into_iter().collect())
}

fn exec_genes_targeting(api_data: &APIData, gene_uniquename: &GeneUniquename,
                        target_of_type: TargetOfType)
    -> GeneUniquenameVecResult
{
    let mut genes_targeting = api_data.genes_targeting(gene_uniquename, target_of_type);

    Ok(genes_targeting.drain().collect())
}

fn exec_gene_properties(api_data: &APIData, query_properties: &HashSet<GeneQueryPropFlag>)
    -> GeneUniquenameVecResult
{
    let mut gene_list = vec![];

    'GENE:
    for gene_query_data in api_data.get_maps().gene_query_data_map.values() {
        for query_prop in query_properties.iter() {
            if !gene_query_data.property_flags.contains(query_prop) {
                continue 'GENE;
            }
        }

        gene_list.push(gene_query_data.gene_uniquename.to_owned());
    }

    Ok(gene_list)
}

async fn exec_query_id(api_data: &APIData,
                       maybe_site_db: &Option<SiteDB>, id: &Uuid)
                       -> GeneUniquenameVecResult
{
    if let Some(site_db) = maybe_site_db {
        if let Some(query) = site_db.query_by_id(id).await {
            match query.exec(api_data, maybe_site_db).await {
                Ok(res) => {
                    Ok(res.iter().map(|row| { row.gene_uniquename.clone() }).collect())
                },
                Err(err) => {
                    Err(err)
                }
            }
        } else {
            Err(flex_fmt!("can't find query for ID {}", id))
        }
    } else {
        Err(flex_fmt!("can't find query for ID {} - no database", id))
    }
}

impl QueryNode {
    pub fn template_node() -> QueryNode {
        QueryNode {
            node_name: None,
            or: None,
            and: None,
            not: None,
            term: None,
            ref_genes: None,
            subset: None,
            gene_list: None,
            target_of: None,
            has_ortholog: None,
            int_range: None,
            float_range: None,
            genome_range: None,
            interactors: None,
            substrates: None,
            downstream_genes: None,
            genes_targeting: None,
            gene_properties: None,
            query_id: None,
        }
    }

    pub fn get_query_id(&self) -> Option<Uuid> {
        self.query_id.as_ref().map(|query_id_node| query_id_node.id)
    }

    #[async_recursion]
    pub async fn exec<'a>(&'a self, api_data: &'a APIData,
                      site_db: &'a Option<SiteDB>) -> GeneUniquenameVecResult {
        if let Some(ref nodes) = self.or {
            return exec_or(api_data, site_db, nodes).await;
        }
        if let Some(ref nodes) = self.and {
            return exec_and(api_data, site_db, nodes).await;
        }
        if let Some(ref not_node) = self.not {
            return exec_not(api_data, site_db, &not_node.node_a, &not_node.node_b).await;
        }
        if let Some(ref term) = self.term {
            return exec_termid(api_data, &term.termid, &term.single_or_multi_locus,
                               &term.ploidiness,
                               &term.expression, &term.conditions,
                               &term.excluded_conditions);
        }
        if let Some(ref ref_genes) = self.ref_genes {
            return exec_ref_genes(api_data, &ref_genes.reference_uniquename);
        }
        if let Some(ref subset_node) = self.subset {
            return exec_subset(api_data, &subset_node.subset_name);
        }
        if let Some(ref gene_list_node) = self.gene_list {
            return exec_gene_list(&gene_list_node.genes);
        }
        if let Some(ref target_of_node) = self.target_of {
            return exec_target_of(api_data, &target_of_node.genes);
        }
        if let Some(ref has_ortholog_node) = self.has_ortholog {
            return exec_has_ortholog(api_data, has_ortholog_node.taxonid);
        }
        if let Some(ref genome_range_node) = self.genome_range {
            return exec_genome_range_overlaps(api_data, genome_range_node.start,
                                              genome_range_node.end,
                                              &genome_range_node.chromosome_name);
        }
        if let Some(ref interactors_node) = self.interactors {
            return match &interactors_node.interaction_type as &str {
                "physical" =>
                    exec_interactors_of_gene(api_data, &interactors_node.gene_uniquename,
                                             InteractionType::Physical),
                "genetic" =>
                    exec_interactors_of_gene(api_data, &interactors_node.gene_uniquename,
                                             InteractionType::Genetic),
                _ => Err(flex_fmt!("No such interaction type: {}",
                                   interactors_node.interaction_type))
            };
        }
        if let Some(ref substrates_node) = self.substrates {
            return exec_substrates_of_gene(api_data, &substrates_node.gene_uniquename,
                                           &substrates_node.phase_term);
        }
        if let Some(ref downstream_genes_node) = self.downstream_genes {
            return exec_downstream_genes(api_data, &downstream_genes_node.gene_uniquename,
                                         &downstream_genes_node.cv_name);
        }
        if let Some(ref genes_targeting_node) = self.genes_targeting {
            return exec_genes_targeting(api_data, &genes_targeting_node.gene_uniquename,
                                        genes_targeting_node.target_of_type);
        }
        if let Some(ref int_range_node) = self.int_range {
            return exec_int_range(api_data, &int_range_node.range_type,
                                  int_range_node.start, int_range_node.end,
                                  &int_range_node.options);
        }
        if let Some(ref float_range_node) = self.float_range {
            return exec_float_range(api_data, &float_range_node.range_type,
                                    float_range_node.start, float_range_node.end,
                                    &float_range_node.options);
        }
        if let Some(ref gene_properties_node) = self.gene_properties {
            return exec_gene_properties(api_data, &gene_properties_node.property_flags);
        }
        if let Some(ref query_id_node) = self.query_id {
            return exec_query_id(api_data, site_db, &query_id_node.id).await;
        }

        // fall through:
        panic!("unsupported query node: {:?}", self);
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct NucleotideDownloadOptions {
    include_introns: bool,
    include_exons: bool,
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
pub struct GAFOptions {
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub aspects: HashSet<FlexStr>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct QueryOutputOptions {
    pub sequence: SeqType,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gaf_options: Option<GAFOptions>,
    pub field_names: Vec<FlexStr>,

    // If a gene in the results is directly or indirectly annotation with one of
    // the ancestor_terms, that term will be added to the "subsets" field of the
    // ResultRow
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub ancestor_terms: HashSet<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub flags: HashSet<FlexStr>,
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

fn get_chr_range(chr_residues: &FlexStr, feature_edge: usize, base_count: usize,
                 before_or_after: BeforeOrAfter) -> FlexStr
{
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

    chr_residues.as_str()[start_pos..end_pos].to_shared_str()
}

impl Query {
    pub fn new(constraints: QueryNode,
               output_options: QueryOutputOptions) -> Query {
        Query {
            output_options,
            constraints
        }
    }

    fn make_nucl_seq(&self, transcript: &TranscriptDetails,
                     maybe_chr_details: Option<&ChromosomeDetails>,
                     options: &NucleotideDownloadOptions) -> String {
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
            if options.include_exons && part.feature_type == FeatureType::Exon ||
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

        seq
    }

    fn make_sequence(&self, api_data: &APIData,
                     gene_uniquename: &FlexStr) -> Option<String> {
        let maybe_gene_summary = api_data.get_gene_summary(gene_uniquename);

        if let Some(gene_summary) = maybe_gene_summary {
            let maybe_transcript = gene_summary.transcripts.get(0);
            if let Some(transcript) = maybe_transcript {
                match self.output_options.sequence {
                    SeqType::Protein =>
                        if let Some(ref protein) = transcript.protein {
                            return Some(protein.sequence.to_string());
                        }
                    SeqType::Nucleotide(ref options) => {
                        let chr = gene_summary.location.as_ref()
                            .map(|l| -> &FlexStr { &l.chromosome_name } )
                            .and_then(|name| api_data.get_chr_details(name));
                        let seq = self.make_nucl_seq(transcript, chr, options);
                        return Some(seq)
                    }
                    SeqType::None => (),
                }
            }
        }

        None
    }

    fn make_gaf_lines(&self, api_data: &APIData,
                      gene_uniquename: &FlexStr) -> Option<String> {

        let gaf_options = &self.output_options.gaf_options;

        if gaf_options.is_none() {
            return None;
        }

        let aspects = &gaf_options.as_ref().unwrap().aspects;

        let maybe_gene_details = api_data.get_gene_details(gene_uniquename);

        let mut gaf_bytes: Vec<u8> = vec![];

        if let Some(gene_details) = maybe_gene_details {
            for aspect in aspects.iter() {
                let result = write_go_annotation_format(&mut gaf_bytes, api_data.get_config(),
                                                        api_data,
                                                        GpadGafWriteMode::PomBaseGaf,
                                                        gene_details.as_ref(),
                                                        &api_data.get_maps().transcripts,
                                                        aspect);

                if result.is_err() {
                    return None;
                }
            }
        }

        if let Ok(gaf_str) = str::from_utf8(&gaf_bytes) {
            if !gaf_str.is_empty() {
                return Some(gaf_str.to_owned())
            }
        }
        None
    }

    fn get_gene_ex_for_results(&self, api_data: &APIData, dataset_name: &FlexStr,
                               gene_uniquename: &FlexStr)
        -> Vec<GeneExValue>
    {
        let mut result = vec![];

        let maybe_datasets =
            api_data.get_maps().gene_expression_measurements.get(gene_uniquename);

        if let Some(datasets) = maybe_datasets {
            if let Some(measurement) = datasets.get(dataset_name) {
                if let Some(ref avg_copies_per_cell) = measurement.avg_copies_per_cell {
                    let new_value = GeneExValue {
                        dataset_name: dataset_name.clone(),
                        value: avg_copies_per_cell.clone(),
                    };
                    result.push(new_value);
                }
            }
        }

        result
    }

    fn add_ancestor_terms(&self, api_data: &APIData,  gene_uniquename: &FlexStr,
                          subsets: &mut HashSet<FlexStr>) {

        let needed_ancestor_terms = &self.output_options.ancestor_terms;
        for ancestor_term in needed_ancestor_terms.iter() {
            if let Some(genes_of_ancestor_term) =
                api_data.get_maps().termid_genes.get(ancestor_term)
            {
                if genes_of_ancestor_term.contains(gene_uniquename) {
                    subsets.insert(ancestor_term.clone());
                }
            }
        }
    }

    fn make_result_rows(&self, api_data: &APIData,
                        genes: Vec<FlexStr>) -> QueryRowsResult {
        Ok(genes.into_iter()
           .map(|gene_uniquename| {
               let mut deletion_viability = None;
               let mut go_component = None;
               let mut go_process_superslim = None;
               let mut go_function = None;
               let mut characterisation_status = None;
               let mut taxonomic_distribution = None;
               let mut ortholog_taxonids = HashSet::new();
               let mut physical_interactors = HashSet::new();
               let mut reference_uniquenames = HashSet::new();
               let mut tmm = None;
               let mut pdb_ids = HashSet::new();
               let mut rnacentral_id = None;
               let mut gocam_ids = HashSet::new();
               let mut paralogs = HashSet::new();
               let mut molecular_weight = None;
               let mut protein_length = None;
               let mut protein_length_bin = None;
               let mut spliced_rna_length = None;
               let mut unspliced_rna_length = None;
               let mut subsets = HashSet::new();
               let mut gene_expression = vec![];

               let maybe_gene_data = api_data.get_gene_query_data(&gene_uniquename);

               if let Some(gene_data) = maybe_gene_data {
                   for field_name in &self.output_options.field_names {
                       match field_name.as_ref() {
                           "deletion_viability" =>
                                deletion_viability = Some(gene_data.deletion_viability.clone()),
                           "ortholog_taxonids" =>
                                ortholog_taxonids = gene_data.ortholog_taxonids.clone(),
                           "physical_interactors" =>
                                physical_interactors = gene_data.physical_interactors.clone(),
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
                           "pdb_ids" => pdb_ids = gene_data.pdb_ids.clone(),
                           "rnacentral_id" => rnacentral_id = gene_data.rnacentral_urs_identifier.clone(),
                           "gocam_ids" => gocam_ids = gene_data.gocam_ids.clone(),
                           "paralogs" => paralogs = gene_data.paralogs.clone(),
                           "molecular_weight" =>
                               molecular_weight = gene_data.molecular_weight,
                           "protein_length" =>
                               protein_length = gene_data.protein_length,
                           "protein_length_bin" =>
                               protein_length_bin = gene_data.protein_length_bin.clone(),
                           "spliced_rna_length" =>
                               spliced_rna_length = gene_data.spliced_rna_length,
                           "unspliced_rna_length" =>
                               unspliced_rna_length = gene_data.unspliced_rna_length,
                           "reference_uniquenames" =>
                               reference_uniquenames = gene_data.reference_uniquenames.clone(),
                           "gene_uniquename" => (),
                           _ => {
                               if let Some(without_prefix) = field_name.strip_prefix("gene_ex_avg_copies_per_cell:") {
                                   let dataset_name = without_prefix.to_shared_str();
                                   let results_measurements =
                                       self.get_gene_ex_for_results(api_data, &dataset_name,
                                                                    &gene_uniquename);
                                   gene_expression.extend(results_measurements);
                               } else {
                                   eprintln!("warning - no such field: {}", field_name);
                               }
                           }
                       }
                   }


                   subsets =
                       if self.output_options.flags.contains(&flex_str!("include_gene_subsets")) {
                           gene_data.subset_termids.clone()
                       } else {
                           HashSet::new()
                       };

                   self.add_ancestor_terms(api_data, &gene_uniquename, &mut subsets);
               }

               let sequence = self.make_sequence(api_data, &gene_uniquename);
               let gaf_lines = self.make_gaf_lines(api_data, &gene_uniquename);

               ResultRow {
                   sequence,
                   gaf_lines,
                   deletion_viability,
                   go_component,
                   go_process_superslim,
                   go_function,
                   characterisation_status,
                   taxonomic_distribution,
                   ortholog_taxonids,
                   physical_interactors,
                   tmm,
                   gene_uniquename,
                   molecular_weight,
                   protein_length,
                   protein_length_bin,
                   spliced_rna_length,
                   unspliced_rna_length,
                   reference_uniquenames,
                   pdb_ids,
                   rnacentral_id,
                   gocam_ids,
                   paralogs,
                   subsets,
                   gene_expression,
               }
           }).collect::<Vec<_>>())
    }

    pub async fn exec(&self, api_data: &APIData,
                      site_db: &Option<SiteDB>)
                -> QueryRowsResult
    {
        let genes_result = self.constraints.exec(api_data, site_db).await;

        match genes_result {
            Ok(genes) => self.make_result_rows(api_data, genes),
            Err(err) => Err(err)
        }
    }

    pub fn get_constraints(&self) -> &QueryNode {
        &self.constraints
    }
}
