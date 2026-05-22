use std::collections::{BTreeSet, HashSet};
use std::num::NonZeroUsize;

use flexstr::{SharedStr as FlexStr, ToFlex};


use crate::data_types::{ActiveSite, AssignedByPeptideRange, BasicProteinFeature, BetaStrand, BindingSite, Chain, ChromosomeLocation, DeletionViability, DisulfideBond, FeatureType, GeneHistoryEntry, GlycosylationSite, GoCamIdAndTitle, Helix, LipidationSite, PDBEntry, ReferenceAndSource, Residues, SynonymDetails, Turn};
use crate::interpro::InterProMatch;
use crate::types::{GeneName, GeneUniquename, ProteinUniquename, RnaUrsId, TermId, TranscriptUniquename};

use crate::api_data::APIData;

pub struct RestExec {
}

impl RestExec {
    pub fn new() -> RestExec {
        RestExec { }
    }

    pub async fn genes_by_id(&self, api_data: &APIData, gene_ids: &[&str])
        -> PublicAPIGeneLookupResponse
    {
        eprintln!("{:?}", gene_ids);
        PublicAPIGeneLookupResponse {
            found: vec![],
            missing: gene_ids.iter().map(|s| s.to_flex()).collect(),
        }
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PubAPIFeaturePart {
    pub feature_type: FeatureType,
    pub systematic_id: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub location: ChromosomeLocation,
    pub residues: Residues,
    pub comment: Option<FlexStr>,
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIProteinDetails {
    pub systematic_id: ProteinUniquename,
    pub sequence: FlexStr,
    pub number_of_residues: usize,  // residue count not including stop codon
    pub product: Option<FlexStr>,
    pub molecular_weight: f32,
    pub average_residue_weight: f32,
    pub charge_at_ph7: f32,
    pub isoelectric_point: f32,
    pub codon_adaptation_index: f32,
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPITranscriptDetails {
    pub systematic_id: TranscriptUniquename,
    pub name: Option<GeneName>,
    pub location: ChromosomeLocation,
    pub parts: Vec<PubAPIFeaturePart>,
    pub transcript_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub protein: Option<PublicAPIProteinDetails>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<ChromosomeLocation>,

    // the CDS length or RNA length without introns - sum of lengths of exons
    pub rna_seq_length_spliced: Option<NonZeroUsize>,
    // the CDS length (protein coding) or RNA length (non-coding) including introns
    pub rna_seq_length_unspliced: Option<NonZeroUsize>,
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIGeneDetails {
    pub systematic_id: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub taxonid: u32,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<FlexStr>,
    pub deletion_viability: DeletionViability,
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub secondary_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub agr_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub biogrid_interactor_id: Option<u32>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub rnacentral_urs_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub interpro_matches: Vec<InterProMatch>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub tm_domain_coords: Vec<AssignedByPeptideRange>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub disordered_region_coords: Vec<AssignedByPeptideRange>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub low_complexity_region_coords: Vec<AssignedByPeptideRange>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub coiled_coil_coords: Vec<AssignedByPeptideRange>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub signal_peptide: Option<BasicProteinFeature>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub transit_peptide: Option<BasicProteinFeature>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub binding_sites: Vec<BindingSite>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub active_sites: Vec<ActiveSite>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub beta_strands: Vec<BetaStrand>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub helices: Vec<Helix>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub turns: Vec<Turn>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub propeptides: Vec<BasicProteinFeature>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub chains: Vec<Chain>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub glycosylation_sites: Vec<GlycosylationSite>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub disulfide_bonds: Vec<DisulfideBond>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub lipidation_sites: Vec<LipidationSite>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub pdb_entries: Vec<PDBEntry>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub orfeome_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub schizosaccharomyces_orthogroup: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub tfexplorer_chipseq_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub tfexplorer_ipms_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pombephosphoproteomics_unige_ch_starvation_mating_gene: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pombephosphoproteomics_unige_ch_fusion_gene: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub name_descriptions: Vec<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<SynonymDetails>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub dbxrefs: HashSet<FlexStr>,

    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    // possible values: "is_histone"
    pub flags: HashSet<FlexStr>,

    pub feature_type: FlexStr,
    pub feature_so_termid: FlexStr,
    pub transcript_so_termid: Option<TermId>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub characterisation_status: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub taxonomic_distribution: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,

    #[serde(skip_serializing_if="BTreeSet::is_empty", default)]
    // transcripts of this gene
    pub transcripts: BTreeSet<PublicAPITranscriptDetails>,

    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub feature_publications: HashSet<ReferenceAndSource>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gocams: Vec<GoCamIdAndTitle>,

    #[serde(skip_serializing_if="Option::is_none")]
    pub rnacentral_2d_structure_id: Option<RnaUrsId>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gene_history: Vec<GeneHistoryEntry>,
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIGeneLookupResponse {
    found: Vec<PublicAPIGeneDetails>,
    missing: Vec<FlexStr>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct RestIdRequest {
    pub ids: String,
}
