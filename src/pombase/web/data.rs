extern crate serde_json;
extern crate postgres;

use std::cmp::min;
use std::fs::{File, create_dir_all};
use std::io::{Write, BufWriter};
use std::collections::HashMap;
use std::fmt::Display;
use std::fmt;

use self::postgres::Connection;

type CvName = String;

pub type TypeInteractionAnnotationMap =
    HashMap<TypeName, Vec<InteractionAnnotation>>;
pub type UniquenameGeneMap =
    HashMap<GeneUniquename, GeneDetails>;
pub type UniquenameTranscriptMap =
    HashMap<TranscriptUniquename, TranscriptDetails>;
pub type UniquenameProteinMap =
    HashMap<ProteinUniquename, ProteinDetails>;
pub type UniquenameReferenceMap =
    HashMap<TermId, ReferenceDetails>;

pub type UniquenameAlleleMap = HashMap<AlleleUniquename, AlleleShort>;
pub type UniquenameGenotypeMap = HashMap<GenotypeUniquename, GenotypeDetails>;
pub type TermIdDetailsMap = HashMap<TermId, TermDetails>;
pub type ChrNameDetailsMap = HashMap<ChromosomeName, ChromosomeDetails>;

pub type IdGenotypeMap = HashMap<GenotypeUniquename, GenotypeDetails>;
pub type IdGeneShortMap = HashMap<GeneUniquename, GeneShort>;
pub type IdRcTermShortMap = HashMap<TermId, Rc<TermShort>>;
pub type IdRcTermDetailsMap = HashMap<TermId, Rc<TermDetails>>;

pub type ReferenceShortMap = HashMap<ReferenceUniquename, ReferenceShort>;
pub type GeneShortMap = HashMap<GeneUniquename, GeneShort>;
pub type GenotypeShortMap = HashMap<GeneUniquename, GenotypeShort>;
pub type AlleleShortMap = HashMap<AlleleUniquename, AlleleShort>;
pub type TermShortMap = HashMap<TermId, TermShort>;

use std::rc::Rc;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::collections::HashSet;

use web::config::*;
use types::*;
use interpro::InterProMatch;

#[derive(Serialize, Clone, Debug, PartialEq, Eq, PartialOrd, Hash)]
pub enum ExtRange {
#[serde(rename = "gene_uniquename")]
    Gene(GeneUniquename),
#[serde(rename = "summary_gene_uniquenames")]
    // the inner Vec length will be > 1 for cases like "binds abc1 and def2, cdc2"
    SummaryGenes(Vec<Vec<String>>),
#[serde(rename = "termid")]
    Term(TermId),
#[serde(rename = "summary_termids")]
    // See: merge_ext_part_ranges()
    SummaryTerms(Vec<TermId>),
#[serde(rename = "misc")]
    Misc(String),
#[serde(rename = "domain")]
    Domain(String),
#[serde(rename = "gene_product")]
    GeneProduct(String),
}

impl ExtRange {
    pub fn is_gene(&self) -> bool {
        match *self {
            ExtRange::Gene(_) => true,
            _ => false,
        }
    }
}

// A single part of an extension.
#[derive(Serialize, Clone, Debug, PartialEq, Eq)]
pub struct ExtPart {
    pub rel_type_name: String,
    pub rel_type_display_name: String,
    pub ext_range: ExtRange,
}
impl Hash for ExtPart {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.rel_type_name.hash(state);
        self.ext_range.hash(state);
    }
}

// minimal information about a gene used in other objects
#[derive(Serialize, Clone, Debug)]
pub struct GeneShort {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,
}

impl PartialEq for GeneShort {
    fn eq(&self, other: &GeneShort) -> bool {
        self.uniquename == other.uniquename
    }
}
impl Eq for GeneShort { }
impl Ord for GeneShort {
    fn cmp(&self, other: &GeneShort) -> Ordering {
        if self.name.is_some() {
            if other.name.is_some() {
                self.name.cmp(&other.name)
            } else { Ordering::Less }
        } else {
            if other.name.is_some() {
                Ordering::Greater
            } else { self.uniquename.cmp(&other.uniquename) }
        }
    }
}
impl PartialOrd for GeneShort {
    fn partial_cmp(&self, other: &GeneShort) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Hash for GeneShort {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.uniquename.hash(state);
    }
}

// a gene uniquename and an organism ID
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct IdAndOrganism {
    pub identifier: String,
    pub taxonid: u32,
}

// identifiers used for autocomplete in the search box
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneSummary {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,
    pub taxonid: OrganismTaxonId,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,
    pub synonyms: Vec<String>,
    pub orthologs: Vec<IdAndOrganism>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    pub feature_type: String,
}

// minimal information about a terms used in other objects
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermShort {
    pub name: TermName,
    pub cv_name: String,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parents: HashSet<String>,
    pub termid: TermId,
    pub is_obsolete: bool,
    pub gene_count: usize,
    pub genotype_count: usize,
}

impl PartialEq for TermShort {
    fn eq(&self, other: &TermShort) -> bool {
        self.termid == other.termid
    }
}
impl Eq for TermShort { }
impl Ord for TermShort {
    fn cmp(&self, other: &TermShort) -> Ordering {
        let order = self.name.cmp(&other.name);
        if order == Ordering::Equal {
            self.termid.cmp(&other.termid)
        } else {
            order
        }
    }
}
impl PartialOrd for TermShort {
    fn partial_cmp(&self, other: &TermShort) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Hash for TermShort {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.termid.hash(state);
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct ChromosomeDetails {
    pub name: String,
    pub residues: String,
    pub ena_identifier: String,
}

#[derive(Serialize, Clone, Debug)]
pub struct ReferenceShort {
    pub uniquename: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<String>,
    pub gene_count: usize,
    pub genotype_count: usize,
}

#[derive(Serialize, Clone, Debug)]
pub struct ReferenceDetails {
    pub uniquename: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<String>,
    #[serde(skip_serializing_if="Option::is_none", rename = "abstract")]
    pub pubmed_abstract: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_publication_date: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_role: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_name: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_approved_date: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_session_submitted_date: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_added_date: Option<String>,
    pub cv_annotations: OntAnnotationMap,
    pub physical_interactions: Vec<InteractionAnnotation>,
    pub genetic_interactions: Vec<InteractionAnnotation>,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
    pub genes_by_uniquename: HashMap<GeneUniquename, GeneShort>,
    pub genotypes_by_uniquename: HashMap<GenotypeUniquename, GenotypeShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub terms_by_termid: HashMap<TermId, TermShort>,
}

// the GO with/from
#[derive(Serialize, Clone, Debug)]
pub enum WithFromValue {
#[serde(rename = "gene")]
    Gene(GeneShort),
#[serde(rename = "term")]
    Term(TermShort),
#[serde(rename = "identifier")]
    Identifier(String)
}

#[derive(Serialize, Clone, Debug)]
pub struct OntAnnotationDetail {
    pub id: i32,
    pub genes: Vec<GeneUniquename>,
    pub reference: Option<ReferenceUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    pub extension: Vec<ExtPart>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub withs: Vec<WithFromValue>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub froms: Vec<WithFromValue>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub residue: Option<Residue>,
    pub qualifiers: Vec<Qualifier>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_ex_props: Option<GeneExProps>,
    // only for genotype/phenotype annotation:
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype: Option<GenotypeUniquename>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub conditions: HashSet<TermId>,
}

impl PartialEq for OntAnnotationDetail {
    fn eq(&self, other: &OntAnnotationDetail) -> bool {
        self.id == other.id
    }
}
impl Eq for OntAnnotationDetail { }

#[derive(Serialize, Clone, Debug)]
pub struct OntTermAnnotations {
    pub term: TermShort,
    pub is_not: bool,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub rel_names: HashSet<RelName>,
    pub annotations: Vec<Rc<OntAnnotationDetail>>,
    pub summary: Option<Vec<TermSummaryRow>>,
}

impl PartialEq for OntTermAnnotations {
    fn eq(&self, other: &OntTermAnnotations) -> bool {
        self.term.termid == other.term.termid
    }
}
impl Eq for OntTermAnnotations { }
impl Ord for OntTermAnnotations {
    fn cmp(&self, other: &OntTermAnnotations) -> Ordering {
        if !self.is_not && other.is_not {
            return Ordering::Less;
        }
        if self.is_not && !other.is_not {
            return Ordering::Greater;
        }
        self.term.name.to_lowercase().cmp(&other.term.name.to_lowercase())
    }
}
impl PartialOrd for OntTermAnnotations {
    fn partial_cmp(&self, other: &OntTermAnnotations) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Hash for OntTermAnnotations {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.term.termid.hash(state);
        self.is_not.hash(state);
    }
}

#[derive(Serialize, Clone, Debug, Eq, PartialEq)]
pub struct TermSummaryRow {
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gene_uniquenames: Vec<GeneUniquename>, // for term and ref pages
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub genotype_uniquenames: Vec<GenotypeUniquename>, // for term pages
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub extension: Vec<ExtPart>,
}

impl Hash for TermSummaryRow {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.gene_uniquenames.hash(state);
        self.genotype_uniquenames.hash(state);
        for ext_part in &self.extension {
            ext_part.hash(state);
        }
    }
}

#[derive(Serialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct TargetOfAnnotation {
    pub ontology_name: String,
    pub ext_rel_display_name: String,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub genes: Vec<GeneUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_uniquename: Option<GenotypeUniquename>,
    pub reference_uniquename: Option<ReferenceUniquename>,
}

#[derive(Serialize, Clone, Debug)]
pub struct SynonymDetails {
    pub name: String,
    #[serde(rename = "type")]
    pub synonym_type: String
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
pub enum Strand {
    #[serde(rename="forward")]
    Forward = 1,
    #[serde(rename="reverse")]
    Reverse = -1,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChromosomeShort {
    pub name: String,
    pub length: usize,
    pub ena_identifier: String,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChromosomeLocation {
    pub chromosome: ChromosomeShort,
    pub start_pos: u32,
    pub end_pos: u32,
    pub strand: Strand,
}

#[derive(Serialize, Clone, Debug)]
pub enum DeletionViability {
    #[serde(rename="viable")]
    Viable,
    #[serde(rename="inviable")]
    Inviable,
    #[serde(rename="depends_on_conditions")]
    DependsOnConditions,
    #[serde(rename="unknown")]
    Unknown,
}

#[derive(Serialize, Clone, Debug)]
pub struct GeneDetails {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    pub taxonid: u32,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<String>,
    pub deletion_viability: DeletionViability,
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<String>,
    pub interpro_matches: Vec<InterProMatch>,
    // non-InterPro domains:
    pub tm_domain_coords: Vec<(usize, usize) >,
    pub orfeome_identifier: Option<String>,
    pub name_descriptions: Vec<String>,
    pub synonyms: Vec<SynonymDetails>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub dbxrefs: HashSet<String>,
    pub feature_type: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub characterisation_status: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<ChromosomeLocation>,
    pub gene_neighbourhood: Vec<GeneShort>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcripts: Vec<TranscriptDetails>,
    pub cv_annotations: OntAnnotationMap,
    pub physical_interactions: Vec<InteractionAnnotation>,
    pub genetic_interactions: Vec<InteractionAnnotation>,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
    pub target_of_annotations: Vec<TargetOfAnnotation>,
    pub references_by_uniquename: HashMap<ReferenceUniquename, ReferenceShort>,
    pub genes_by_uniquename: HashMap<GeneUniquename, GeneShort>,
    pub genotypes_by_uniquename: HashMap<GenotypeUniquename, GenotypeShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub terms_by_termid: HashMap<TermId, TermShort>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinDetails {
    pub uniquename: TranscriptUniquename,
    pub sequence: String,
    pub molecular_weight: f32,
    pub average_residue_weight: f32,
    pub charge_at_ph7: f32,
    pub isoelectric_point: f32,
    pub codon_adaptation_index: f32,
}

pub type Residues = String;

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum FeatureType {
#[serde(rename = "five_prime_utr")]
    FivePrimeUtr,
#[serde(rename = "five_prime_utr_intron")]
    FivePrimeUtrIntron,
#[serde(rename = "exon")]
    Exon,
#[serde(rename = "cds_intron")]
    // type for introns between exons
    CdsIntron,
#[serde(rename = "three_prime_utr")]
    ThreePrimeUtr,
#[serde(rename = "three_prime_utr_intron")]
    ThreePrimeUtrIntron,
}

impl Display for FeatureType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(match *self {
            FeatureType::FivePrimeUtr => "5'UTR",
            FeatureType::FivePrimeUtrIntron => "5'UTR_intron",
            FeatureType::Exon => "exon",
            FeatureType::CdsIntron => "cds_intron",
            FeatureType::ThreePrimeUtr => "3'UTR",
            FeatureType::ThreePrimeUtrIntron => "3'UTR_intron",
        })
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct FeatureShort {
    pub feature_type: FeatureType,
    pub uniquename: String,
    pub location: ChromosomeLocation,
    pub residues: Residues,
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TranscriptDetails {
    pub uniquename: TranscriptUniquename,
    pub parts: Vec<FeatureShort>,
    pub transcript_type: String,
    pub protein: Option<ProteinDetails>,
}

#[derive(Serialize, Clone, Debug)]
pub struct GenotypeShort {
    pub uniquename: GenotypeUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub background: Option<String>,
    pub expressed_alleles: Vec<ExpressedAllele>,
}

#[derive(Serialize, Clone, Debug)]
pub struct GenotypeDetails {
    pub uniquename: GenotypeUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub background: Option<String>,
    pub expressed_alleles: Vec<ExpressedAllele>,
    pub cv_annotations: OntAnnotationMap,
    pub references_by_uniquename: HashMap<ReferenceUniquename, ReferenceShort>,
    pub genes_by_uniquename: HashMap<GeneUniquename, GeneShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub terms_by_termid: HashMap<TermId, TermShort>,
}

#[derive(Serialize, Clone, Debug)]
pub struct ExpressedAllele {
    #[serde(skip_serializing_if="Option::is_none")]
    pub expression: Option<String>,
    pub allele_uniquename: AlleleUniquename,
}

#[derive(Serialize, Clone, Debug)]
pub struct AlleleShort {
    pub uniquename: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    pub allele_type: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub description: Option<String>,
    pub gene_uniquename: GeneUniquename,
}

pub type RelName = String;

#[derive(Serialize, Clone, Debug)]
pub struct GeneExProps {
    #[serde(skip_serializing_if="Option::is_none")]
    pub copies_per_cell: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub avg_copies_per_cell: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub scale: Option<String>,
}

pub type OntName = String;
pub type OntAnnotationMap = HashMap<OntName, Vec<OntTermAnnotations>>;

#[derive(Serialize, Clone, Debug)]
pub struct TermAndRelation {
    pub termid: TermId,
    pub term_name: TermName,
    pub relation_name: RelName,
}

#[derive(Serialize, Clone, Debug)]
pub struct TermDetails {
    pub name: TermName,
    pub cv_name: CvName,
    pub annotation_feature_type: String,
    pub interesting_parents: HashSet<String>,
    pub subsets: Vec<String>,
    pub termid: TermId,
    #[serde(skip_serializing_if="Option::is_none")]
    pub definition: Option<TermDef>,
    pub direct_ancestors: Vec<TermAndRelation>,
    pub genes_annotated_with: HashSet<GeneUniquename>,
    pub is_obsolete: bool,
    pub single_allele_genotype_uniquenames: HashSet<String>,
    pub cv_annotations: OntAnnotationMap,
    pub genes_by_uniquename: HashMap<GeneUniquename, GeneShort>,
    pub genotypes_by_uniquename: HashMap<GenotypeUniquename, GenotypeShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub references_by_uniquename: HashMap<ReferenceUniquename, ReferenceShort>,
    pub terms_by_termid: HashMap<TermId, TermShort>,
}

#[derive(Serialize, Clone, Debug)]
pub struct InteractionAnnotation {
    pub gene_uniquename: GeneUniquename,
    pub interactor_uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference_uniquename: Option<ReferenceUniquename>,
}
impl PartialEq for InteractionAnnotation {
    fn eq(&self, other: &Self) -> bool {
        if let Some(ref evidence) = self.evidence {
            if let Some(ref other_evidence) = other.evidence {
                return evidence == other_evidence;
            }
        }
        (&self.gene_uniquename, &self.interactor_uniquename) ==
            (&other.gene_uniquename, &other.interactor_uniquename)
    }
}
impl Eq for InteractionAnnotation { }
impl Ord for InteractionAnnotation {
    fn cmp(&self, other: &Self) -> Ordering {
        if let Some(ref evidence) = self.evidence {
            if let Some(ref other_evidence) = other.evidence {
                let order = evidence.cmp(&other_evidence);
                if order != Ordering::Equal {
                    return order;
                }
            }
        }
        (&self.gene_uniquename, &self.interactor_uniquename)
            .cmp(&(&other.gene_uniquename, &other.interactor_uniquename))
    }
}
impl PartialOrd for InteractionAnnotation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct OrthologAnnotation {
    pub gene_uniquename: GeneUniquename,
    pub ortholog_taxonid: u32,
    pub ortholog_uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference_uniquename: Option<ReferenceUniquename>,
}
impl PartialEq for OrthologAnnotation {
    fn eq(&self, other: &Self) -> bool {
        (&self.gene_uniquename, &self.ortholog_uniquename) ==
            (&other.gene_uniquename, &other.ortholog_uniquename)
    }
}
impl Eq for OrthologAnnotation { }
impl Ord for OrthologAnnotation {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.gene_uniquename, &self.ortholog_uniquename)
            .cmp(&(&other.gene_uniquename, &other.ortholog_uniquename))
    }
}
impl PartialOrd for OrthologAnnotation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct ParalogAnnotation {
    pub gene_uniquename: GeneUniquename,
    pub paralog_uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference_uniquename: Option<ReferenceUniquename>,
}
impl PartialEq for ParalogAnnotation {
    fn eq(&self, other: &Self) -> bool {
        (&self.gene_uniquename, &self.paralog_uniquename) ==
            (&other.gene_uniquename, &other.paralog_uniquename)
    }
}
impl Eq for ParalogAnnotation { }
impl Ord for ParalogAnnotation {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.gene_uniquename, &self.paralog_uniquename)
            .cmp(&(&other.gene_uniquename, &other.paralog_uniquename))
    }
}
impl PartialOrd for ParalogAnnotation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Metadata {
    pub db_creation_datetime: String,
    pub export_prog_name: String,
    pub export_prog_version: String,
    pub gene_count: usize,
    pub term_count: usize,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIAlleleDetails {
    pub gene: GeneUniquename,
    pub allele_type: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub expression: Option<String>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIGenotypeAnnotation {
    pub is_multi: bool,
    pub alleles: Vec<APIAlleleDetails>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIGeneSummary {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<String>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub exact_synonyms: Vec<String>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub dbxrefs: HashSet<String>,
    pub location: Option<ChromosomeLocation>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<ChromosomeLocation>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcripts: Vec<TranscriptDetails>,
    pub tm_domain_count: usize,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIMaps {
    pub termid_genes: HashMap<TermId, HashSet<GeneUniquename>>,
    pub termid_genotype_annotation: HashMap<TermId, Vec<APIGenotypeAnnotation>>,
    pub gene_summaries: Vec<APIGeneSummary>,
    pub term_summaries: HashSet<TermShort>,
}

#[derive(Serialize, Clone, Debug)]
pub struct RecentReferences {
    // most recent papers from PubMed
    pub pubmed: Vec<ReferenceShort>,
    // most recent admin curated papers
    pub admin_curated: Vec<ReferenceShort>,
    // most recent community curated
    pub community_curated: Vec<ReferenceShort>,
}

#[derive(Serialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct TermSubsetElement {
    pub name: String,
    pub termid: TermId,
    pub gene_count: usize,
}

#[derive(Serialize, Clone, Debug)]
pub struct TermSubsetDetails {
    pub name: String,
    pub total_gene_count: usize, // total unique genes in all subsets
    pub elements: HashSet<TermSubsetElement>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneSubsetDetails {
    pub name: String,
    pub display_name: String,
    pub elements: HashSet<GeneUniquename>,
}

pub type IdTermSubsetMap = HashMap<String, TermSubsetDetails>;
pub type IdGeneSubsetMap = HashMap<String, GeneSubsetDetails>;

#[derive(Serialize, Clone, Debug)]
pub struct WebData {
    pub genes: UniquenameGeneMap,
    pub genotypes: IdGenotypeMap,
    pub terms: IdRcTermDetailsMap,
    pub used_terms: IdRcTermDetailsMap,
    pub metadata: Metadata,
    pub chromosomes: ChrNameDetailsMap,
    pub references: UniquenameReferenceMap,
    pub recent_references: RecentReferences,
    pub api_maps: APIMaps,
    pub search_gene_summaries: Vec<GeneSummary>,
    pub term_subsets: IdTermSubsetMap,
    pub gene_subsets: IdGeneSubsetMap,
}

impl WebData {
    fn get_genes(&self) -> &UniquenameGeneMap {
        &self.genes
    }
    fn get_genotypes(&self) -> &IdGenotypeMap {
        &self.genotypes
    }
    fn get_references(&self) -> &UniquenameReferenceMap {
        &self.references
    }
    fn get_chromosomes(&self) -> &ChrNameDetailsMap {
        &self.chromosomes
    }
    fn get_terms(&self) -> &IdRcTermDetailsMap {
        &self.terms
    }

    fn create_dir(&self, output_dir: &str, dir_name: &str) -> String {
        let path = String::new() + output_dir + "/" + dir_name;
        create_dir_all(&path).unwrap_or_else(|why| {
            println!("Creating output directory failed: {:?}", why.kind());
        });
        path
    }

    fn write_chromosome_seq_chunks(&self, output_dir: &str, chunk_sizes: &Vec<usize>) {
        for chunk_size in chunk_sizes {
        for (chromosome_uniquename, chromosome_details) in &self.chromosomes {
            let new_path_part = &format!("{}/sequence/{}", chromosome_uniquename, chunk_size);
            let chr_path = self.create_dir(output_dir, new_path_part);
            let mut index = 0;
            let max_index = chromosome_details.residues.len() / chunk_size;
            while index <= max_index {
                let start_pos = index*chunk_size;
                let end_pos = min(start_pos+chunk_size, chromosome_details.residues.len());
                let chunk: String = chromosome_details.residues[start_pos..end_pos].into();
                let file_name = format!("{}/chunk_{}", chr_path, index);
                let f = File::create(file_name).expect("Unable to open file");
                let mut writer = BufWriter::new(&f);
                writer.write_all(chunk.as_bytes()).expect("Unable to write chromosome chunk");
                index += 1;
            }
        }
        }
    }

    fn write_chromosomes(&self, config: &Config, output_dir: &str) {
        let new_path = self.create_dir(output_dir, "chromosome");
        for (chromosome_uniquename, chromosome_details) in &self.chromosomes {
            let s = serde_json::to_string(&chromosome_details).unwrap();
            let file_name = format!("{}/{}.json", new_path, &chromosome_uniquename);
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write chromosome JSON");
        }
        self.write_chromosome_seq_chunks(&new_path, &config.api_seq_chunk_sizes);
    }

    fn write_reference_details(&self, output_dir: &str) {
        let new_path = self.create_dir(output_dir, "reference");
        for (reference_uniquename, reference_details) in &self.references {
            let s = serde_json::to_string(&reference_details).unwrap();
            let file_name = format!("{}/{}.json", new_path, &reference_uniquename);
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write reference JSON");
        }
    }

    fn write_gene_details(&self, output_dir: &str) {
        let new_path = self.create_dir(output_dir, "gene");
        for (gene_uniquename, gene_details) in &self.genes {
            let s = serde_json::to_string(&gene_details).unwrap();
            let file_name = format!("{}/{}.json", new_path, &gene_uniquename);
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write gene JSON");
        }
    }

    fn write_genotype_details(&self, output_dir: &str) {
        let new_path = self.create_dir(output_dir, "genotype");
        for (genotype_uniquename, genotype_details) in &self.genotypes {
            let s = serde_json::to_string(&genotype_details).unwrap();
            let file_name = format!("{}/{}.json", new_path, &genotype_uniquename);
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write genotype JSON");
        }
    }

    fn write_gene_summaries(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.search_gene_summaries).unwrap();
        let file_name = String::new() + &output_dir + "/gene_summaries.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write gene_summaries.json");
    }

    fn write_terms(&self, output_dir: &str) {
        let new_path = self.create_dir(output_dir, "term");
        for (termid, term_details) in &self.terms {
            let s = serde_json::to_string(&term_details).unwrap();
            let file_name = format!("{}/{}.json", new_path, &termid);
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write term JSON");
        }
    }

    fn write_metadata(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.metadata).unwrap();
        let file_name = String::new() + &output_dir + "/metadata.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write metadata.json");
    }

    fn write_recent_references(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.recent_references).unwrap();
        let file_name = String::new() + &output_dir + "/recent_references.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write recent references JSON");
    }

    fn write_api_maps(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.api_maps).unwrap();
        let file_name = String::new() + &output_dir + "/api_maps.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }

    fn write_subsets(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.term_subsets).unwrap();
        let file_name = String::new() + &output_dir + "/term_subsets.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");

        let s = serde_json::to_string(&self.gene_subsets).unwrap();
        let file_name = String::new() + &output_dir + "/gene_subsets.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }

    pub fn write(&self, config: &Config, output_dir: &str) {
        self.write_chromosomes(config, output_dir);
        println!("wrote {} chromosomes", self.get_chromosomes().len());
        self.write_reference_details(output_dir);
        println!("wrote {} references", self.get_references().len());
        self.write_gene_details(output_dir);
        println!("wrote {} genes", self.get_genes().len());
        self.write_genotype_details(output_dir);
        println!("wrote {} genotypes", self.get_genotypes().len());
        self.write_gene_summaries(output_dir);
        println!("wrote gene summaries");
        self.write_terms(output_dir);
        println!("wrote {} terms", self.get_terms().len());
        self.write_metadata(output_dir);
        println!("wrote metadata");
        self.write_recent_references(output_dir);
        println!("wrote recent references");
        self.write_api_maps(output_dir);
        println!("wrote search data");
        self.write_subsets(output_dir);
        println!("wrote subsets");
    }

    pub fn store_jsonb(&self, conn: &Connection) {
        let trans = conn.transaction().unwrap();

        for (uniquename, gene_details) in &self.genes {
            let serde_value = serde_json::value::to_value(&gene_details).unwrap();
            trans.execute("INSERT INTO web_json.gene (uniquename, data) values ($1, $2)",
                          &[&uniquename, &serde_value]).unwrap();
        }
        for (uniquename, ref_details) in &self.references {
            let serde_value = serde_json::value::to_value(&ref_details).unwrap();
            trans.execute("INSERT INTO web_json.reference (uniquename, data) values ($1, $2)",
                          &[&uniquename, &serde_value]).unwrap();
        }
        for (termid, term_details) in &self.terms {
            let serde_value = serde_json::value::to_value(&term_details).unwrap();
            trans.execute("INSERT INTO web_json.term (termid, data) values ($1, $2)",
                         &[&termid, &serde_value]).unwrap();
        }

        trans.execute("CREATE INDEX gene_jsonb_idx ON web_json.gene USING gin (data jsonb_path_ops)", &[]).unwrap();
        trans.execute("CREATE INDEX gene_jsonb_name_idx ON web_json.gene USING gin ((data->>'name') gin_trgm_ops);", &[]).unwrap();
        trans.execute("CREATE INDEX term_jsonb_idx ON web_json.term USING gin (data jsonb_path_ops)", &[]).unwrap();
        trans.execute("CREATE INDEX term_jsonb_name_idx ON web_json.term USING gin ((data->>'name') gin_trgm_ops);", &[]).unwrap();
        trans.execute("CREATE INDEX reference_jsonb_idx ON web_json.reference USING gin (data jsonb_path_ops)", &[]).unwrap();
        trans.execute("CREATE INDEX reference_jsonb_title_idx ON web_json.reference USING gin ((data->>'title') gin_trgm_ops);", &[]).unwrap();

        trans.commit().unwrap();
    }
}
