use std::collections::HashSet;
use std::num::NonZeroUsize;

use flexstr::{SharedStr as FlexStr, ToFlex, shared_str as flex_str};
use itertools::{Either, Itertools};

use crate::data_types::{ActiveSite, AssignedByPeptideRange, BasicProteinFeature, BetaStrand,
                        BindingSite, Chain, ChromosomeLocation, DeletionViability,
                        DisulfideBond, FeatureShort, FeatureType, GeneDetails, GeneHistoryEntry,
                        GlycosylationSite, GoCamIdAndTitle, Helix, LipidationSite, PDBEntry,
                        ProteinDetails, Residues, SynonymDetails, TranscriptDetails, Turn,
                        OrthologAnnotation};
use crate::interpro::InterProMatch;
use crate::types::{GeneName, GeneUniquename, ProteinUniquename, RnaUrsId, TermId,
                   TranscriptUniquename, Evidence, ReferenceUniquename};

use crate::api_data::APIData;

pub struct RestExec {
}

impl RestExec {
    pub fn new() -> RestExec {
        RestExec { }
    }

    pub async fn gene_by_id(&self, api_data: &APIData, gene_id: &str)
        -> Option<PublicAPIGeneDetails>
    {
        api_data.get_full_gene_details(gene_id).as_ref()
            .map(|gene_details| (gene_details as &GeneDetails).into())
    }

    pub async fn genes_by_id(&self, api_data: &APIData, gene_ids: &[&str])
        -> PublicAPIGeneLookupResponse
    {
        let (found, missing) = gene_ids.iter()
            .partition_map(|id| {
                let id = id.to_flex();
                if let Some(ref gene_details) = api_data.get_full_gene_details(&id) {
                    let gd: PublicAPIGeneDetails = (gene_details as &GeneDetails).into();
                    Either::Left(gd)
                } else {
                    Either::Right(id)
                }
            });

        PublicAPIGeneLookupResponse {
            found,
            missing,
        }
    }
}

impl Default for RestExec {
    fn default() -> Self {
        RestExec::new()
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIFeaturePart {
    pub feature_type: FeatureType,
    pub systematic_id: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub location: ChromosomeLocation,
    pub residues: Residues,
}

impl From<&FeatureShort> for PublicAPIFeaturePart {
    fn from(feat: &FeatureShort) -> Self {
        PublicAPIFeaturePart {
            feature_type: feat.feature_type.clone(),
            systematic_id: feat.uniquename.clone(),
            name: feat.name.clone(),
            location: feat.location.clone(),
            residues: feat.residues.clone(),
        }
    }
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

impl From<&ProteinDetails> for PublicAPIProteinDetails {
    fn from(prot: &ProteinDetails) -> Self {
        PublicAPIProteinDetails {
            systematic_id: prot.uniquename.clone(),
            sequence: prot.sequence.clone(),
            number_of_residues: prot.number_of_residues,
            product: prot.product.clone(),
            molecular_weight: prot.molecular_weight,
            average_residue_weight: prot.average_residue_weight,
            charge_at_ph7: prot.charge_at_ph7,
            isoelectric_point: prot.isoelectric_point,
            codon_adaptation_index: prot.codon_adaptation_index,
        }
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPITranscriptDetails {
    pub systematic_id: TranscriptUniquename,
    pub name: Option<GeneName>,
    pub location: ChromosomeLocation,
    pub parts: Vec<PublicAPIFeaturePart>,
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

impl From<&TranscriptDetails> for PublicAPITranscriptDetails {
    fn from(tr: &TranscriptDetails) -> Self {
        let protein = tr.protein.as_ref().map(|p| p.into());
        PublicAPITranscriptDetails {
            systematic_id: tr.uniquename.clone(),
            name: tr.name.clone(),
            location: tr.location.clone(),
            parts: tr.parts.iter().map(|p| p.into()).collect(),
            transcript_type: tr.transcript_type.clone(),
            protein,
            cds_location: tr.cds_location.clone(),
            rna_seq_length_spliced: tr.rna_seq_length_spliced,
            rna_seq_length_unspliced: tr.rna_seq_length_unspliced,
        }
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIOrthologAnnotation {
    pub gene_systematic_id: GeneUniquename,
    pub ortholog_taxonid: u32,
    pub ortholog_systematic_id: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference: Option<ReferenceUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub qualifier: Option<FlexStr>,
}

impl From<&OrthologAnnotation> for PublicAPIOrthologAnnotation {
    fn from(ortholog_annotation: &OrthologAnnotation) -> Self {
        PublicAPIOrthologAnnotation {
            gene_systematic_id: ortholog_annotation.gene_uniquename.clone(),
            ortholog_taxonid: ortholog_annotation.ortholog_taxonid,
            ortholog_systematic_id: ortholog_annotation.ortholog_uniquename.clone(),
            evidence: ortholog_annotation.evidence.clone(),
            reference: ortholog_annotation.reference_uniquename.clone(),
            qualifier: ortholog_annotation.qualifier.clone(),
        }
    }
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

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub ortholog_annotations: Vec<PublicAPIOrthologAnnotation>,

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

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcripts: Vec<PublicAPITranscriptDetails>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gocams: Vec<GoCamIdAndTitle>,

    #[serde(skip_serializing_if="Option::is_none")]
    pub rnacentral_2d_structure_id: Option<RnaUrsId>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gene_history: Vec<GeneHistoryEntry>,
}

impl From<&GeneDetails> for PublicAPIGeneDetails {
    fn from(gene: &GeneDetails) -> Self {
        let transcripts = gene.transcripts.iter()
            .filter_map(|transcript_uniquename| {
                let maybe_maybe_details = gene.transcripts_by_uniquename.get(transcript_uniquename);

                let maybe_details = maybe_maybe_details?;

                let Some(details) = maybe_details
                else {
                    return None
                };

                let api_details: PublicAPITranscriptDetails = details.into();
                Some(api_details)
            })
            .collect();

        let feature_type =
            if gene.feature_type == "mRNA gene" {
                flex_str!("protein")
            } else {
                gene.feature_type.clone()
            };

        let ortholog_annotations = gene.ortholog_annotations.iter()
            .map(|orth| orth.into()).collect();

        PublicAPIGeneDetails {
            systematic_id: gene.uniquename.clone(),
            name: gene.name.clone(),
            taxonid: gene.taxonid,
            product: gene.product.clone(),
            deletion_viability: gene.deletion_viability.clone(),
            uniprot_identifier: gene.uniprot_identifier.clone(),
            secondary_identifier: gene.secondary_identifier.clone(),
            agr_identifier: gene.agr_identifier.clone(),
            biogrid_interactor_id: gene.biogrid_interactor_id,
            rnacentral_urs_identifier: gene.rnacentral_urs_identifier.clone(),
            interpro_matches: gene.interpro_matches.clone(),
            tm_domain_coords: gene.tm_domain_coords.clone(),
            disordered_region_coords: gene.disordered_region_coords.clone(),
            low_complexity_region_coords: gene.low_complexity_region_coords.clone(),
            coiled_coil_coords: gene.coiled_coil_coords.clone(),
            signal_peptide: gene.signal_peptide.clone(),
            transit_peptide: gene.transit_peptide.clone(),
            binding_sites: gene.binding_sites.clone(),
            active_sites: gene.active_sites.clone(),
            beta_strands: gene.beta_strands.clone(),
            helices: gene.helices.clone(),
            turns: gene.turns.clone(),
            propeptides: gene.propeptides.clone(),
            chains: gene.chains.clone(),
            glycosylation_sites: gene.glycosylation_sites.clone(),
            disulfide_bonds: gene.disulfide_bonds.clone(),
            lipidation_sites: gene.lipidation_sites.clone(),
            pdb_entries: gene.pdb_entries.clone(),
            orfeome_identifier: gene.orfeome_identifier.clone(),
            schizosaccharomyces_orthogroup: gene.schizosaccharomyces_orthogroup.clone(),
            tfexplorer_chipseq_identifier: gene.tfexplorer_chipseq_identifier.clone(),
            tfexplorer_ipms_identifier: gene.tfexplorer_ipms_identifier.clone(),
            pombephosphoproteomics_unige_ch_starvation_mating_gene: gene.pombephosphoproteomics_unige_ch_starvation_mating_gene.clone(),
            pombephosphoproteomics_unige_ch_fusion_gene: gene.pombephosphoproteomics_unige_ch_fusion_gene.clone(),
            name_descriptions: gene.name_descriptions.clone(),
            synonyms: gene.synonyms.clone(),
            dbxrefs: gene.dbxrefs.clone(),
            flags: gene.flags.clone(),
            feature_type,
            feature_so_termid: gene.feature_so_termid.clone(),
            transcript_so_termid: gene.transcript_so_termid.clone(),
            characterisation_status: gene.characterisation_status.clone(),
            taxonomic_distribution: gene.taxonomic_distribution.clone(),
            location: gene.location.clone(),
            transcripts,
            ortholog_annotations,
            gocams: gene.gocams.clone(),
            rnacentral_2d_structure_id: gene.rnacentral_2d_structure_id.clone(),
            gene_history: gene.gene_history.clone(),
        }
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIGeneLookupResponse {
    found: Vec<PublicAPIGeneDetails>,
    missing: Vec<FlexStr>,
}
