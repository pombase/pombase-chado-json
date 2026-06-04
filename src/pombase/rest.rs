use std::collections::{BTreeSet, HashSet};
use std::num::NonZeroUsize;

use flexstr::{SharedStr as FlexStr, ToFlex, shared_str as flex_str};
use itertools::{Either, Itertools};

use crate::bio::ExportCommentsMode;
use crate::bio::go_format_writer::{GpadGafWriteMode, make_gaf_line};
use crate::bio::phenotype_format_writer::{FypoEvidenceType, make_phenotype_line_parts};
use crate::data_types::{ActiveSite, AssignedByPeptideRange, BasicProteinFeature, BetaStrand, BindingSite, Chain, ChromosomeLocation, DeletionViability, DisulfideBond, ExpressedAllele, Expression, FeatureShort, FeatureType, GeneDetails, GeneHistoryEntry, GeneShort, GenotypeDetails, GenotypeLocus, GlycosylationSite, GoCamIdAndTitle, Helix, LipidationSite, OntAnnotationDetail, OrthologAnnotation, PDBEntry, ProteinDetails, Residues, SynonymDetails, Throughput, TranscriptDetails, Turn};
use crate::interpro::InterProMatch;
use crate::types::{AlleleUniquename, Evidence, GeneName, GeneProduct, GeneUniquename, GenotypeDisplayName, GenotypeDisplayUniquename, ProteinUniquename, ReferenceUniquename, RnaUrsId, TermId, TermName, TranscriptUniquename};

use crate::api_data::APIData;
use crate::data_types::DataLookup;
use crate::web::config::Config;

pub struct RestExec {
}

#[derive(Deserialize, Clone, Debug)]
#[serde(rename_all = "lowercase")]
pub enum PublicAPIPhenotypeOutputType {
    PHAF,
    JSON,
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

    pub async fn gene_by_uniprot_accession(&self, api_data: &APIData, uniprot_accession: &str)
        -> Option<PublicAPIGeneDetails>
    {
        api_data.get_gene_details_by_uniprot_accession(uniprot_accession).as_ref()
            .map(|gene_details| (gene_details as &GeneDetails).into())
    }

    pub async fn genes_by_uniprot_accession(&self, api_data: &APIData, uniprot_accessions: &[&str])
        -> PublicAPIGeneLookupResponse
    {
        let (found, missing) = uniprot_accessions.iter()
            .partition_map(|id| {
                let id = id.to_flex();
                if let Some(ref details) = api_data.get_gene_details_by_uniprot_accession(&id) {
                    let gd: PublicAPIGeneDetails = (details as &GeneDetails).into();
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

    // return None if the termid doesn't exist or doesn't have annotations
    pub async fn phenotype_annotation_by_termid(&self, config: &Config,
                                                api_data: &APIData, ancestor_termid: &str,
                                                output_type: PublicAPIPhenotypeOutputType)
        -> Option<String>
    {
        let ancestor_termid = ancestor_termid.to_flex();
        let mut lines = vec![];
        let mut annotations = vec![];

        let ancestor_term_details = api_data.get_term(&ancestor_termid)?;

        for term_annotations in ancestor_term_details.cv_annotations.values() {
            for term_annotation in term_annotations {
                let termid = &term_annotation.term;
                let term_name = &api_data.get_term(termid).unwrap().name;

                for annotation_id in &term_annotation.annotations {
                    let annotation_details = api_data.get_annotation_detail(*annotation_id).unwrap();
                    let annotation_details = annotation_details.as_ref();
                    let genotype_uniquename = annotation_details.genotype.as_ref().unwrap();
                    let genotype_details =
                        api_data.get_genotype(genotype_uniquename).unwrap();

                    match output_type {
                        PublicAPIPhenotypeOutputType::PHAF => {
                            if let Some(line) =
                                make_phenotype_line_parts(config, api_data, termid,
                                                          annotation_details,
                                                          &genotype_details,
                                                          FypoEvidenceType::PomBase) {
                                    lines.push(line);
                                }
                        },
                        PublicAPIPhenotypeOutputType::JSON => {
                            let phenotype_annotation =
                                make_phenotype_annotation(api_data,
                                                          termid.clone(), term_name.clone(),
                                                          annotation_details, &genotype_details);
                            annotations.push(phenotype_annotation);
                        }
                    }
                }
            }
        }

        match output_type {
            PublicAPIPhenotypeOutputType::PHAF =>
                Some(lines.iter().map(|line| line.join("\t")).join("\n")),
            PublicAPIPhenotypeOutputType::JSON => Some(serde_json::to_string(&annotations).unwrap())
        }
    }


    // return None if the termid doesn't exist or doesn't have annotations
    pub async fn go_annotation_by_termid(&self, config: &Config,
                                         api_data: &APIData, ancestor_termid: &str)
        -> Option<String>
    {
        let ancestor_termid = ancestor_termid.to_flex();
        let mut lines = vec![];

        let ancestor_term_details = api_data.get_term(&ancestor_termid)?;

        for (cv_name, term_annotations) in ancestor_term_details.cv_annotations.iter() {
            for term_annotation in term_annotations {
                let termid = &term_annotation.term;
                if term_annotation.is_not {
                    continue;
                }

                for annotation_id in &term_annotation.annotations {
                    let annotation_details = api_data.get_annotation_detail(*annotation_id).unwrap();
                    let annotation_details = annotation_details.as_ref();
                    let gene_uniquename = &annotation_details.genes[0];
                    let gene_details = api_data.get_gene(gene_uniquename).unwrap();
                    let gene_details = gene_details.as_ref();

                    let Some(line) =
                        make_gaf_line(config, api_data, GpadGafWriteMode::GafForRest,
                                      ExportCommentsMode::NoExport, gene_details,
                                      annotation_details, termid, false, cv_name)
                    else {
                        continue;
                    };

                    lines.push(line);
                }
            }
        }
        Some(lines.iter().map(|line| line.join("\t")).join("\n"))
    }
}


impl Default for RestExec {
    fn default() -> Self {
        RestExec::new()
    }
}

fn make_phenotype_annotation(api_data: &APIData,
                             termid: TermId, term_name: TermName,
                             annotation_details: &OntAnnotationDetail,
                             genotype_details: &GenotypeDetails)
    -> PublicAPIPhenotypeAnnotation
 {
        PublicAPIPhenotypeAnnotation {
             genotype: make_genotype(api_data, genotype_details),
             termid,
             term_name,
             conditions: annotation_details.condition_details.clone(),
             date: annotation_details.date.clone(),
             throughput: annotation_details.throughput,
             evidence: annotation_details.evidence.clone(),
             eco_evidence: annotation_details.eco_evidence.clone(),
             reference: annotation_details.reference.clone(),
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
    pub taxonid: u32,
    pub systematic_id: GeneUniquename,
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
            taxonid: ortholog_annotation.ortholog_taxonid,
            systematic_id: ortholog_annotation.ortholog_uniquename.clone(),
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
    pub orthologs: Vec<PublicAPIOrthologAnnotation>,

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

        let orthologs = gene.ortholog_annotations.iter().map(|orth| orth.into()).collect();

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
            orthologs,
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

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIGeneShort {
    pub systematic_id: GeneUniquename,
    pub name: Option<GeneName>,
    pub product: Option<GeneProduct>,
}

impl From<&GeneShort> for PublicAPIGeneShort {
    fn from(gene_short: &GeneShort) -> Self {
        PublicAPIGeneShort {
            systematic_id: gene_short.uniquename.clone(),
            name: gene_short.name.clone(),
            product: gene_short.product.clone(),
        }
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIAllele {
    pub name: Option<FlexStr>,
    pub allele_type: FlexStr,
    pub description: Option<FlexStr>,
    pub gene: PublicAPIGeneShort,
    pub synonyms: Vec<SynonymDetails>,
}

fn make_allele(api_data: &APIData, allele_uniquename: &AlleleUniquename)
    -> PublicAPIAllele
{
    let allele = api_data.get_allele(allele_uniquename).unwrap();
    let gene = &allele.gene;

    PublicAPIAllele {
        name: allele.name.clone(),
        allele_type: allele.allele_type.clone(),
        description: allele.description.clone(),
        gene: gene.into(),
        synonyms: allele.synonyms.clone(),
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIExpressedAllele {
    pub expression: Option<Expression>,
    pub promoter_gene: Option<FlexStr>,
    pub allele: PublicAPIAllele,
}

fn make_expressed_allele(api_data: &APIData, ea: &ExpressedAllele)
    -> PublicAPIExpressedAllele
{
    PublicAPIExpressedAllele {
        expression: ea.expression.clone(),
        promoter_gene: ea.promoter_gene.clone(),
        allele: make_allele(api_data, &ea.allele_uniquename),
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIGenotypeLocus {
    pub expressed_alleles: Vec<PublicAPIExpressedAllele>,
}

fn make_locus(api_data: &APIData, locus: &GenotypeLocus)
    -> PublicAPIGenotypeLocus
{
   let expressed_alleles =
       locus.expressed_alleles.iter().map(|ea| make_expressed_allele(api_data, ea)).collect();
   PublicAPIGenotypeLocus {
       expressed_alleles,
   }
}


#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIGenotype {
    pub display_uniquename: GenotypeDisplayUniquename,
    pub display_name: GenotypeDisplayName,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub loci: Vec<PublicAPIGenotypeLocus>,
}

fn make_genotype(api_data: &APIData, gd: &GenotypeDetails)
    -> PublicAPIGenotype
{
    PublicAPIGenotype {
        display_uniquename: gd.display_uniquename.clone(),
        display_name: gd.display_name.clone(),
        name: gd.name.clone(),
        loci: gd.loci.iter().map(|l| make_locus(api_data, l)).collect(),
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIPhenotypeAnnotation {
    pub genotype: PublicAPIGenotype,
    pub termid: TermId,
    pub term_name: TermName,
    pub conditions: BTreeSet<(TermId, Option<String>)>,
    pub date: Option<FlexStr>,
    pub throughput: Option<Throughput>,
    pub evidence: Option<Evidence>,
    pub eco_evidence: Option<Evidence>,
    pub reference: Option<ReferenceUniquename>,
}

