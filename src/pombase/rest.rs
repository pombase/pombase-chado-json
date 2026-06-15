use std::collections::{BTreeSet, HashSet, HashMap};
use std::num::NonZeroUsize;

use flexstr::{SharedStr as FlexStr, ToFlex, shared_str as flex_str};
use itertools::{Either, Itertools};

use crate::bio::ExportCommentsMode;
use crate::bio::go_format_writer::{GpadGafWriteMode, make_gaf_line};
use crate::bio::phenotype_format_writer::{FypoEvidenceType, make_phenotype_line_parts};
use crate::constants::{FYPO_CV_NAME, is_go_root_name};

use crate::data_types::{ActiveSite, AnnotationCurator, AnnotationExtension, AssignedByPeptideRange, BasicProteinFeature, BetaStrand, BindingSite, ChromosomeLocation, DeletionViability, DisulfideBond, ExpressedAllele, Expression, FeatureShort, FeatureType, GeneDetails, GeneHistoryEntry, GeneShort, GenotypeDetails, GenotypeLocus, GlycosylationSite, GoCamIdAndTitle, Helix, LipidationSite, OntAnnotationDetail, OrthologAnnotation, PDBEntry, Phase, ProteinDetails, ReferenceDetails, Residues, RheaId, Strand, SynonymDetails, TermAndRelation, TermDetails, TermXref, Throughput, TranscriptDetails, Turn, WithFromValue};
use crate::interpro::InterProMatch;
use crate::types::{AlleleUniquename, CvName, Evidence, GeneName, GeneProduct, GeneUniquename, GenotypeDisplayName, GenotypeDisplayUniquename, ProteinUniquename, Qualifier, ReferenceUniquename, RnaUrsId, TermDef, TermId, TermName, TranscriptUniquename};

use crate::api_data::APIData;
use crate::data_types::DataLookup;
use crate::web::config::{Config, InterestingParent};

pub struct RestExec {
}

#[derive(Deserialize, Clone, Copy, Debug)]
#[serde(rename_all = "lowercase")]
pub enum PublicAPIOutputType {
    TSV,
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
            .map(|gene_details| gene_details.into())
    }

    pub async fn genes_by_id(&self, api_data: &APIData, gene_ids: &[&str])
        -> PublicAPIGeneLookupResponse
    {
        let (found, missing) = gene_ids.iter()
            .partition_map(|id| {
                let id = id.to_flex();
                if let Some(ref gene_details) = api_data.get_full_gene_details(&id) {
                    let gd: PublicAPIGeneDetails = gene_details.into();
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
            .map(|gene_details| gene_details.into())
    }

    pub async fn genes_by_uniprot_accession(&self, api_data: &APIData, uniprot_accessions: &[&str])
        -> PublicAPIGeneLookupResponse
    {
        let (found, missing) = uniprot_accessions.iter()
            .partition_map(|id| {
                let id = id.to_flex();
                if let Some(ref details) = api_data.get_gene_details_by_uniprot_accession(&id) {
                    let gd: PublicAPIGeneDetails = details.into();
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
    pub fn phenotype_annotation_by_termid(&self, config: &Config,
                                          api_data: &dyn DataLookup,
                                          ancestor_termids: &[&str],
                                          output_type: PublicAPIOutputType)
        -> Option<String>
    {
        let mut lines = vec![];
        let mut annotations = HashSet::new();

        for ancestor_termid in ancestor_termids {
            let ancestor_termid = ancestor_termid.to_flex();
            let ancestor_term_details = api_data.get_term(&ancestor_termid)?;

            if ancestor_term_details.cv_name != FYPO_CV_NAME {
                continue;
            }

            for term_annotations in ancestor_term_details.cv_annotations.values() {
                for term_annotation in term_annotations {
                    let termid = &term_annotation.term;
                    let term_name = &api_data.get_term(termid).unwrap().name;

                    for annotation_id in &term_annotation.annotations {
                        let Some(annotation_details) = api_data.get_annotation_detail(*annotation_id)
                        else {
                            continue;
                        };
                        let annotation_details = annotation_details.as_ref();
                        let Some(genotype_uniquename) = annotation_details.genotype.as_ref()
                        else {
                            continue;
                        };
                        let genotype_details =
                            api_data.get_genotype(genotype_uniquename).unwrap();

                        match output_type {
                            PublicAPIOutputType::TSV => {
                                if let Some(line) =
                                    make_phenotype_line_parts(config, api_data, termid,
                                                              annotation_details,
                                                              &genotype_details,
                                                              FypoEvidenceType::PomBase) {
                                        lines.push(line);
                                    }
                            },
                            PublicAPIOutputType::JSON => {
                                let phenotype_annotation =
                                    make_phenotype_annotation(api_data,
                                                              termid.clone(), term_name.clone(),
                                                              annotation_details, &genotype_details);
                                annotations.insert(phenotype_annotation);
                            },

                        }
                    }
                }
            }
        }

        match output_type {
            PublicAPIOutputType::TSV =>
                Some(lines.iter().map(|line| line.join("\t")).join("\n")),
            PublicAPIOutputType::JSON => Some(serde_json::to_string(&annotations).unwrap())
        }
    }


    // return None if the termid doesn't exist or doesn't have annotations
    pub fn go_annotation_by_termid(&self, config: &Config,
                                   api_data: &dyn DataLookup,
                                   ancestor_termids: &[&str],
                                   output_type: PublicAPIOutputType)
        -> Option<String>
    {
        let mut lines = vec![];
        let mut annotations = vec![];

        for ancestor_termid in ancestor_termids {
            let ancestor_termid = ancestor_termid.to_flex();
            let ancestor_term_details = api_data.get_term(&ancestor_termid)?;

            if !is_go_root_name(&ancestor_term_details.cv_name) {
                continue;
            }

            for (cv_name, term_annotations) in ancestor_term_details.cv_annotations.iter() {
                for term_annotation in term_annotations {
                    let termid = &term_annotation.term;

                    if term_annotation.is_not {
                        continue;
                    }

                    let term_name = &api_data.get_term(termid).unwrap().name;

                    for annotation_id in &term_annotation.annotations {
                        let annotation_details = api_data.get_annotation_detail(*annotation_id).unwrap();
                        let annotation_details = annotation_details.as_ref();
                        let gene_uniquename = &annotation_details.genes[0];
                        let gene_details = api_data.get_gene(gene_uniquename).unwrap();
                        let gene_details = gene_details.as_ref();

                        if let Some(ref characterisation_status) = gene_details.characterisation_status
                            && (characterisation_status == "dubious" || characterisation_status == "transposon") {
                                continue;
                            }

                        if gene_details.feature_type == "pseudogene" {
                            continue;
                        }

                        match output_type {
                            PublicAPIOutputType::TSV => {
                                let Some(line) =
                                    make_gaf_line(config, api_data, GpadGafWriteMode::GafForRest,
                                                  ExportCommentsMode::NoExport, gene_details,
                                                  annotation_details, termid, false, cv_name)
                                else {
                                    continue;
                                };

                                lines.push(line);
                            },
                            PublicAPIOutputType::JSON => {
                                let go_annotation =
                                    make_go_annotation(config, api_data,
                                                       termid.clone(), term_name.clone(),
                                                       annotation_details);
                                annotations.push(go_annotation);
                            }
                        }
                    }
                }
            }
        }

        match output_type {
            PublicAPIOutputType::TSV =>
                Some(lines.iter().map(|line| line.join("\t")).join("\n")),
            PublicAPIOutputType::JSON =>
                Some(serde_json::to_string(&annotations).unwrap())
        }
    }
}


impl Default for RestExec {
    fn default() -> Self {
        RestExec::new()
    }
}

fn make_go_annotation(config: &Config, api_data: &dyn DataLookup,
                      termid: TermId, term_name: TermName,
                      annotation_details: &OntAnnotationDetail)
    -> PublicAPIGOAnnotation
{
    let gene_uniquename = &annotation_details.genes[0];
    let gene = api_data.get_gene(gene_uniquename).unwrap();

    let db_prefix = &config.database_name;

    let add_db_prefix = |v: &WithFromValue| {
        let id = v.id();
        if id.contains(":") {
            id
        } else {
            format!("{}:{}", db_prefix, id).to_flex()
        }
    };

    let with = annotation_details.withs.iter()
        .map(add_db_prefix).collect();
    let from = annotation_details.froms.iter()
        .map(add_db_prefix).collect();

    PublicAPIGOAnnotation {
        gene_systematic_id: gene_uniquename.to_owned(),
        gene_name: gene.name.clone(),
        product: gene.product.clone(),
        taxonid: gene.taxonid,
        feature_type: gene.feature_type.clone(),
        feature_so_termid: gene.feature_so_termid.clone(),
        transcript_so_termid: gene.transcript_so_termid.clone(),
        termid,
        term_name,
        annotation_extension: annotation_details.extension.clone(),
        with,
        from,
        gene_product_form_id: annotation_details.gene_product_form_id.clone(),
        qualifiers: annotation_details.qualifiers.clone(),
        date: annotation_details.date.clone(),
        throughput: annotation_details.throughput,
        evidence: annotation_details.evidence.clone(),
        eco_evidence: annotation_details.eco_evidence.clone(),
        reference: annotation_details.reference.clone(),
    }
}

fn make_phenotype_annotation(api_data: &dyn DataLookup,
                             termid: TermId, term_name: TermName,
                             annotation_details: &OntAnnotationDetail,
                             gd: &GenotypeDetails)
    -> PublicAPIPhenotypeAnnotation
{
    let conditions = annotation_details.condition_details.iter()
        .map(|(termid, _)| {
            let term_details = api_data.get_term(termid).unwrap();
            let name = term_details.name.to_std_string();
            PublicAPICondition {
                termid: termid.clone(),
                name: Some(name),
            }
        })
        .collect();

    PublicAPIPhenotypeAnnotation {
        genotype_display_uniquename: gd.display_uniquename.clone(),
        genotype_display_name: gd.display_name.clone(),
        genotype_name: gd.name.clone(),
        genotype_loci: gd.loci.iter().map(|l| make_locus(api_data, l)).collect(),
        termid,
        term_name,
        conditions,
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
    pub chromosome_location: PublicAPIChromosomeLocation,
    pub residues: Residues,
}

impl From<&FeatureShort> for PublicAPIFeaturePart {
    fn from(feat: &FeatureShort) -> Self {
        PublicAPIFeaturePart {
            feature_type: feat.feature_type.clone(),
            systematic_id: feat.uniquename.clone(),
            name: feat.name.clone(),
            chromosome_location: (&feat.location).into(),
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
    pub chromosome_location: PublicAPIChromosomeLocation,
    pub parts: Vec<PublicAPIFeaturePart>,
    pub transcript_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub protein: Option<PublicAPIProteinDetails>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<PublicAPIChromosomeLocation>,

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
            chromosome_location: (&tr.location).into(),
            parts: tr.parts.iter().map(|p| p.into()).collect(),
            transcript_type: tr.transcript_type.clone(),
            protein,
            cds_location: tr.cds_location.as_ref().map(|l| l.into()),
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
    pub symbol: Option<FlexStr>,
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

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub name_descriptions: Vec<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<SynonymDetails>,

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
    pub chromosome_location: Option<PublicAPIChromosomeLocation>,

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
            symbol: gene.name.clone(),
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
            glycosylation_sites: gene.glycosylation_sites.clone(),
            disulfide_bonds: gene.disulfide_bonds.clone(),
            lipidation_sites: gene.lipidation_sites.clone(),
            pdb_entries: gene.pdb_entries.clone(),
            orfeome_identifier: gene.orfeome_identifier.clone(),
            schizosaccharomyces_orthogroup: gene.schizosaccharomyces_orthogroup.clone(),
            tfexplorer_chipseq_identifier: gene.tfexplorer_chipseq_identifier.clone(),
            tfexplorer_ipms_identifier: gene.tfexplorer_ipms_identifier.clone(),
            name_descriptions: gene.name_descriptions.clone(),
            synonyms: gene.synonyms.clone(),
            flags: gene.flags.clone(),
            feature_type,
            feature_so_termid: gene.feature_so_termid.clone(),
            transcript_so_termid: gene.transcript_so_termid.clone(),
            characterisation_status: gene.characterisation_status.clone(),
            taxonomic_distribution: gene.taxonomic_distribution.clone(),
            chromosome_location: gene.location.as_ref().map(|l| l.into()),
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

#[derive(Serialize, Clone, Debug, Hash, Eq, PartialEq)]
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

#[derive(Serialize, Clone, Debug, Hash, Eq, PartialEq)]
pub struct PublicAPIAllele {
    pub name: Option<FlexStr>,
    pub allele_type: FlexStr,
    pub description: Option<FlexStr>,
    pub gene: PublicAPIGeneShort,
    pub synonyms: Vec<SynonymDetails>,
}

fn make_allele(api_data: &dyn DataLookup, allele_uniquename: &AlleleUniquename)
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

#[derive(Serialize, Clone, Debug, Hash, Eq, PartialEq)]
pub struct PublicAPIExpressedAllele {
    pub expression: Option<Expression>,
    pub promoter_gene: Option<FlexStr>,
    pub allele: PublicAPIAllele,
}

fn make_expressed_allele(api_data: &dyn DataLookup, ea: &ExpressedAllele)
    -> PublicAPIExpressedAllele
{
    PublicAPIExpressedAllele {
        expression: ea.expression.clone(),
        promoter_gene: ea.promoter_gene.clone(),
        allele: make_allele(api_data, &ea.allele_uniquename),
    }
}

#[derive(Serialize, Clone, Debug, Hash, Eq, PartialEq)]
pub struct PublicAPIGenotypeLocus {
    pub expressed_alleles: Vec<PublicAPIExpressedAllele>,
}

fn make_locus(api_data: &dyn DataLookup, locus: &GenotypeLocus)
    -> PublicAPIGenotypeLocus
{
   let expressed_alleles =
       locus.expressed_alleles.iter().map(|ea| make_expressed_allele(api_data, ea)).collect();
   PublicAPIGenotypeLocus {
       expressed_alleles,
   }
}

#[derive(Serialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PublicAPICondition {
    pub termid: FlexStr,
    pub name: Option<String>,
}

#[derive(Serialize, Clone, Debug, Hash, Eq, PartialEq)]
pub struct PublicAPIPhenotypeAnnotation {
    pub genotype_display_uniquename: GenotypeDisplayUniquename,
    pub genotype_display_name: GenotypeDisplayName,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_name: Option<FlexStr>,
    pub genotype_loci: Vec<PublicAPIGenotypeLocus>,
    pub termid: TermId,
    pub term_name: TermName,
    pub conditions: BTreeSet<PublicAPICondition>,
    pub date: Option<FlexStr>,
    pub throughput: Option<Throughput>,
    pub evidence: Option<Evidence>,
    pub eco_evidence: Option<Evidence>,
    pub reference: Option<ReferenceUniquename>,
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIGOAnnotation {
    pub gene_systematic_id: GeneUniquename,
    pub gene_name: Option<GeneName>,
    pub product: Option<GeneProduct>,
    pub feature_type: FlexStr,
    pub taxonid: u32,
    pub feature_so_termid: FlexStr,
    pub transcript_so_termid: Option<TermId>,
    pub termid: TermId,
    pub term_name: TermName,
    pub annotation_extension: AnnotationExtension,
    pub with: BTreeSet<FlexStr>,
    pub from: BTreeSet<FlexStr>,
    pub gene_product_form_id: Option<FlexStr>,
    pub qualifiers: Vec<Qualifier>,
    pub date: Option<FlexStr>,
    pub throughput: Option<Throughput>,
    pub evidence: Option<Evidence>,
    pub eco_evidence: Option<Evidence>,
    pub reference: Option<ReferenceUniquename>,
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIChromosomeLocation {
    pub chromosome_name: FlexStr,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    #[serde(skip_serializing_if="Option::is_none")]
    pub phase: Option<Phase>,
}

impl From<&ChromosomeLocation> for PublicAPIChromosomeLocation {
    fn from(loc: &ChromosomeLocation) -> Self {
        PublicAPIChromosomeLocation {
            chromosome_name: loc.chromosome_name.clone(),
            start: loc.start_pos,
            end: loc.end_pos,
            strand: loc.strand,
            phase: loc.phase,
        }
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct PublicAPIReferenceDetails {
    pub uniquename: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none", rename = "abstract")]
    pub pubmed_abstract: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none", rename = "doi")]
    pub pubmed_doi: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_publication_date: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_entrez_date: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_session_key: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_annotation_status: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_triage_status: Option<FlexStr>,
    pub canto_curator_role: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_name: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_first_approved_date: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_approved_date: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_approver_orcid: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_session_submitted_date: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_added_date: Option<FlexStr>,

    // the curators of the annotations from Canto, may be different from the canto_curator_name
    pub annotation_curators: Vec<AnnotationCurator>,

    #[serde(skip_serializing_if="Option::is_none")]
    pub file_curator_name: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub file_curator_role: Option<FlexStr>,

    pub annotation_file_curators: Vec<AnnotationCurator>,

    pub genes: Vec<FlexStr>,
    pub gene_count: usize,

    // count of genes annotated in LTP experiments
    pub ltp_gene_count: usize,

    // This is set to the year part of canto_first_approved_date if it is
    // not None, otherwise set to the year part of canto_approved_date, otherwise
    // canto_session_submitted_date
    #[serde(skip_serializing_if="Option::is_none")]
    pub approved_date: Option<FlexStr>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub pdb_entries: Vec<PDBEntry>,
}

impl From<&ReferenceDetails> for PublicAPIReferenceDetails {
    fn from(ref_details: &ReferenceDetails) -> Self {
        let genes = ref_details.genes_by_uniquename.keys().cloned().collect();
        PublicAPIReferenceDetails {
            uniquename: ref_details.uniquename.clone(),
            title: ref_details.title.clone(),
            citation: ref_details.citation.clone(),
            pubmed_abstract: ref_details.pubmed_abstract.clone(),
            pubmed_doi: ref_details.pubmed_doi.clone(),
            authors: ref_details.authors.clone(),
            authors_abbrev: ref_details.authors_abbrev.clone(),
            pubmed_publication_date: ref_details.pubmed_publication_date.clone(),
            pubmed_entrez_date: ref_details.pubmed_entrez_date.clone(),
            publication_year: ref_details.publication_year.clone(),
            canto_session_key: ref_details.canto_session_key.clone(),
            canto_annotation_status: ref_details.canto_annotation_status.clone(),
            canto_triage_status: ref_details.canto_triage_status.clone(),
            canto_curator_role: ref_details.canto_curator_role.clone(),
            canto_curator_name: ref_details.canto_curator_name.clone(),
            canto_first_approved_date: ref_details.canto_first_approved_date.clone(),
            canto_approved_date: ref_details.canto_approved_date.clone(),
            canto_approver_orcid: ref_details.canto_approver_orcid.clone(),
            canto_session_submitted_date: ref_details.canto_session_submitted_date.clone(),
            canto_added_date: ref_details.canto_added_date.clone(),
            annotation_curators: ref_details.annotation_curators.clone(),
            file_curator_name: ref_details.file_curator_name.clone(),
            file_curator_role: ref_details.file_curator_role.clone(),
            annotation_file_curators: ref_details.annotation_file_curators.clone(),
            genes,
            gene_count: ref_details.gene_count,
            ltp_gene_count: ref_details.ltp_gene_count,
            approved_date: ref_details.approved_date.clone(),
            pdb_entries: ref_details.pdb_entries.clone(),
        }
    }
}



#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct PublicAPITermDetails {
    pub termid: TermId,
    pub name: TermName,
    pub cv_name: CvName,
    pub annotation_feature_type: FlexStr,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parent_ids: HashSet<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parent_details: HashSet<InterestingParent>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub in_subsets: HashSet<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<SynonymDetails>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub definition: Option<TermDef>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub direct_ancestors: Vec<TermAndRelation>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub definition_xrefs: HashSet<TermId>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub secondary_identifiers: HashSet<TermId>,

    // genes annotated with this term, except if the term is a phenotype
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub annotated_genes: HashSet<GeneUniquename>,

    // genes in single locus genotypes annotated with this term
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub single_locus_annotated_genes: HashSet<GeneUniquename>,

    // genes in multi locus genotypes annotated with this term
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub multi_locus_annotated_genes: HashSet<GeneUniquename>,

    pub is_obsolete: bool,

    // count of genes annotated with this term, excluding NOT annotations
    pub gene_count: usize,

    pub genotype_count: usize,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub xrefs: HashMap<FlexStr, TermXref>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pombase_gene_id: Option<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub gocams: HashSet<GoCamIdAndTitle>,

    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    // for ChEBI terms, the IDs of reactions containing this chemical:
    pub rhea_reaction_ids: HashSet<RheaId>,
}

impl From<&TermDetails> for PublicAPITermDetails {
    fn from(td: &TermDetails) -> Self {
        PublicAPITermDetails {
            termid: td.termid.clone(),
            name: td.name.clone(),
            cv_name: td.cv_name.clone(),
            annotation_feature_type: td.annotation_feature_type.clone(),
            interesting_parent_ids: td.interesting_parent_ids.clone(),
            interesting_parent_details: td.interesting_parent_details.clone(),
            in_subsets: td.in_subsets.clone(),
            synonyms: td.synonyms.clone(),
            definition: td.definition.clone(),
            direct_ancestors: td.direct_ancestors.clone(),
            definition_xrefs: td.definition_xrefs.clone(),
            secondary_identifiers: td.secondary_identifiers.clone(),
            annotated_genes: td.annotated_genes.clone(),
            single_locus_annotated_genes: td.single_locus_annotated_genes.clone(),
            multi_locus_annotated_genes: td.multi_locus_annotated_genes.clone(),
            is_obsolete: td.is_obsolete,
            gene_count: td.gene_count,
            genotype_count: td.genotype_count,
            xrefs: td.xrefs.clone(),
            pombase_gene_id: td.pombase_gene_id.clone(),
            gocams: td.gocams.clone(),
            rhea_reaction_ids: td.rhea_reaction_ids.clone(),
        }
    }
}
