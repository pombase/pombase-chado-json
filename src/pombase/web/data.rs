extern crate serde_json;
extern crate postgres;

use std::cmp::min;
use std::fs::{File, create_dir_all};
use std::io::{Write, BufWriter};
use std::io;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Display;
use std::fmt;

use regex::Regex;

use bio::util::format_fasta;

use flate2::Compression;
use flate2::write::GzEncoder;

use self::postgres::Connection;

type CvName = String;

pub type TypeInteractionAnnotationMap =
    HashMap<TypeName, Vec<InteractionAnnotation>>;
pub type UniquenameGeneMap =
    BTreeMap<GeneUniquename, GeneDetails>;
pub type UniquenameTranscriptMap =
    HashMap<TranscriptUniquename, TranscriptDetails>;
pub type UniquenameProteinMap =
    HashMap<ProteinUniquename, ProteinDetails>;
pub type UniquenameReferenceMap =
    HashMap<TermId, ReferenceDetails>;

pub type UniquenameAlleleMap = HashMap<AlleleUniquename, AlleleShort>;
pub type UniquenameGenotypeMap = HashMap<GenotypeUniquename, GenotypeDetails>;
pub type TermIdDetailsMap = HashMap<TermId, TermDetails>;
pub type ChrNameDetailsMap = BTreeMap<ChromosomeName, ChromosomeDetails>;

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

const FASTA_SEQ_COLUMNS: usize = 60;

fn write_as_fasta(writer: &mut Write, id: &str, desc: Option<String>, seq: &str) {
    let fasta = format_fasta(id, desc, &seq, FASTA_SEQ_COLUMNS);
    writer.write(fasta.as_bytes()).unwrap();
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, PartialOrd, Hash)]
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
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
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
#[derive(Serialize, Deserialize, Clone, Debug)]
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
    pub uniprot_identifier: Option<String>,
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

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChromosomeDetails {
    pub name: String,
    pub residues: String,
    pub ena_identifier: String,
    pub gene_uniquenames: Vec<String>,
    pub taxonid: OrganismTaxonId,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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
    pub canto_triage_status: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_role: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_name: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_first_approved_date: Option<String>,
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
#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum WithFromValue {
#[serde(rename = "gene")]
    Gene(GeneShort),
#[serde(rename = "term")]
    Term(TermShort),
#[serde(rename = "identifier")]
    Identifier(String)
}

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct OntTermAnnotations {
    pub term: TermShort,
    pub is_not: bool,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub rel_names: HashSet<RelName>,
    pub annotations: Vec<OntAnnotationDetail>,
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

#[derive(Serialize, Deserialize, Clone, Debug, Eq, PartialEq)]
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

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct TargetOfAnnotation {
    pub ontology_name: String,
    pub ext_rel_display_name: String,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub genes: Vec<GeneUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_uniquename: Option<GenotypeUniquename>,
    pub reference_uniquename: Option<ReferenceUniquename>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<ChromosomeLocation>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GenotypeShort {
    pub uniquename: GenotypeUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub background: Option<String>,
    pub expressed_alleles: Vec<ExpressedAllele>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ExpressedAllele {
    #[serde(skip_serializing_if="Option::is_none")]
    pub expression: Option<String>,
    pub allele_uniquename: AlleleUniquename,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermAndRelation {
    pub termid: TermId,
    pub term_name: TermName,
    pub relation_name: RelName,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermDetails {
    pub name: TermName,
    pub cv_name: CvName,
    pub annotation_feature_type: String,
    pub interesting_parents: HashSet<String>,
    pub subsets: Vec<String>,
    pub termid: TermId,
    pub synonyms: Vec<SynonymDetails>,
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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
                let order = evidence.cmp(other_evidence);
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
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
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcripts: Vec<TranscriptDetails>,
    pub tm_domain_count: usize,
    pub exon_count: usize,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIMaps {
    pub termid_genes: HashMap<TermId, HashSet<GeneUniquename>>,
    pub termid_genotype_annotation: HashMap<TermId, Vec<APIGenotypeAnnotation>>,
    pub gene_summaries: HashMap<GeneUniquename, APIGeneSummary>,
    pub term_summaries: HashSet<TermShort>,
    pub genes: UniquenameGeneMap,
    pub gene_name_gene_map: HashMap<String, GeneUniquename>,
    pub genotypes: IdGenotypeMap,
    pub terms: HashMap<TermId, TermDetails>,
    pub references: UniquenameReferenceMap,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrTermSummary {
    pub id: TermId,
    pub name: TermName,
    pub cv_name: CvName,
    #[serde(skip_serializing_if="Option::is_none")]
    pub definition: Option<TermDef>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub close_synonyms: Vec<String>,   // exact and narrow
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub distant_synonyms: Vec<String>, // broad and related
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parents: HashSet<String>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrData {
    pub term_summaries: HashMap<TermId, SolrTermSummary>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct RecentReferences {
    // most recent papers from PubMed
    pub pubmed: Vec<ReferenceShort>,
    // most recent admin curated papers
    pub admin_curated: Vec<ReferenceShort>,
    // most recent community curated
    pub community_curated: Vec<ReferenceShort>,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct TermSubsetElement {
    pub name: String,
    pub termid: TermId,
    pub gene_count: usize,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
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

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct WebData {
    pub metadata: Metadata,
    pub chromosomes: ChrNameDetailsMap,
    pub recent_references: RecentReferences,
    pub api_maps: APIMaps,
    pub solr_data: SolrData,
    pub search_gene_summaries: Vec<GeneSummary>,
    pub term_subsets: IdTermSubsetMap,
    pub gene_subsets: IdGeneSubsetMap,
}

impl WebData {
    fn get_chromosomes(&self) -> &ChrNameDetailsMap {
        &self.chromosomes
    }

    fn create_dir(&self, output_dir: &str, dir_name: &str) -> String {
        let path = String::new() + output_dir + "/" + dir_name;
        create_dir_all(&path).unwrap_or_else(|why| {
            println!("Creating output directory failed: {:?}", why.kind());
        });
        path
    }

    fn write_chromosome_seq_chunks(&self, output_dir: &str, chunk_sizes: &[usize]) {
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

    fn write_chromosome_json(&self, config: &Config, output_dir: &str) {
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

    fn write_gene_summaries(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.search_gene_summaries).unwrap();
        let file_name = String::new() + output_dir + "/gene_summaries.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write gene_summaries.json");
    }

    fn write_metadata(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.metadata).unwrap();
        let file_name = String::new() + output_dir + "/metadata.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write metadata.json");
    }

    fn write_recent_references(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.recent_references).unwrap();
        let file_name = String::new() + output_dir + "/recent_references.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write recent references JSON");
    }

    fn write_api_maps(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.api_maps).unwrap();
        let file_name = String::new() + output_dir + "/api_maps.json.gz";
        let f = File::create(file_name).expect("Unable to open file");

        let mut compressor = GzEncoder::new(f, Compression::Default);
        compressor.write_all(s.as_bytes()).unwrap();
        compressor.finish().unwrap();
    }

    fn write_solr_data(&self, output_dir: &str) {
        let new_path = self.create_dir(output_dir, "solr_data/terms");
        for (termid, term_summary) in &self.solr_data.term_summaries {
            let s = serde_json::to_string(&term_summary).unwrap();
            let file_name = format!("{}/{}.json", new_path, &termid);
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write term JSON");
        }
    }

    fn write_subsets(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.term_subsets).unwrap();
        let file_name = String::new() + output_dir + "/term_subsets.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write");

        let s = serde_json::to_string(&self.gene_subsets).unwrap();
        let file_name = String::new() + output_dir + "/gene_subsets.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write");
    }

    fn write_feature_sequences(&self, output_dir: &str) {
        let make_seq_writer = |name: &str| {
            let file_name = String::new() + output_dir + "/" + name;
            let file = File::create(file_name).expect("Unable to open file");
            BufWriter::new(file)
        };

        let mut cds_writer = make_seq_writer("cds.fa");
        let mut cds_introns_writer = make_seq_writer("cds+introns.fa");
        let mut cds_introns_utrs_writer = make_seq_writer("cds+introns+utrs.fa");
        let mut peptide_writer = make_seq_writer("peptide.fa");

        for (gene_uniquename, gene_details) in &self.api_maps.genes {
            if let Some(transcript) = gene_details.transcripts.get(0) {
                let mut cds_seq = String::new();
                let mut cds_introns_seq = String::new();
                let mut cds_introns_utrs_seq = String::new();
                for part in &transcript.parts {
                    if part.feature_type == FeatureType::Exon {
                        cds_seq += &part.residues;
                        cds_introns_seq += &part.residues;
                    }
                    if part.feature_type == FeatureType::CdsIntron {
                        cds_introns_seq += &part.residues;
                    }
                    cds_introns_utrs_seq += &part.residues;
                }

                write_as_fasta(&mut cds_writer, gene_uniquename, None, &cds_seq);
                write_as_fasta(&mut cds_introns_writer, gene_uniquename, None, &cds_introns_seq);
                write_as_fasta(&mut cds_introns_utrs_writer,
                               gene_uniquename, None, &cds_introns_utrs_seq);
                if let Some(ref protein) = transcript.protein {
                    write_as_fasta(&mut peptide_writer, &(gene_uniquename.to_owned() + ":pep"),
                                   None, &protein.sequence);
                }
            }
        }

        cds_writer.flush().unwrap();
        cds_introns_utrs_writer.flush().unwrap();
        peptide_writer.flush().unwrap();
    }

    pub fn write_chromosome_sequences(&self, config: &Config, output_dir: &str) {
        let make_seq_writer = |name: &str| {
            let file_name = String::new() + output_dir + "/" + name;
            let file = File::create(file_name).expect("Unable to open file");
            BufWriter::new(file)
        };

        let load_org_name = config.load_organism().full_name();
        let chromosomes_file_name = load_org_name.clone() + "_all_chromosomes.fa";
        let mut chromosomes_writer = make_seq_writer(&chromosomes_file_name);

        for (uniquename, details) in &self.chromosomes {
            write_as_fasta(&mut chromosomes_writer, &uniquename, Some(load_org_name.clone()),
                           &details.residues);

            let this_chr_file_name = load_org_name.clone() + "_chromosome_" +
                uniquename + ".fa";
            let mut this_chr_writer = make_seq_writer(&this_chr_file_name);
            write_as_fasta(&mut this_chr_writer, uniquename, Some(load_org_name.clone()),
                           &details.residues);
            this_chr_writer.flush().unwrap();

        }

        chromosomes_writer.flush().unwrap();
    }

    fn write_gene_id_table(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        let load_org_taxonid = config.load_organism_taxonid;

        let gene_file_name = output_dir.to_owned() + "/sysID2product.tsv";
        let rna_file_name = output_dir.to_owned() + "/sysID2product.rna.tsv";
        let pseudogenes_file_name = output_dir.to_owned() + "/pseudogeneIDs.tsv";
        let all_names_file_name = output_dir.to_owned() + "/allNames.tsv";

        let gene_file = File::create(gene_file_name).expect("Unable to open file");
        let rna_file = File::create(rna_file_name).expect("Unable to open file");
        let pseudogenes_file = File::create(pseudogenes_file_name).expect("Unable to open file");
        let all_names_file = File::create(all_names_file_name).expect("Unable to open file");

        let mut gene_writer = BufWriter::new(&gene_file);
        let mut rna_writer = BufWriter::new(&rna_file);
        let mut pseudogenes_writer = BufWriter::new(&pseudogenes_file);
        let mut all_names_writer = BufWriter::new(&all_names_file);

        let db_version = format!("# Chado database date: {}\n", self.metadata.db_creation_datetime);
        gene_writer.write(db_version.as_bytes())?;
        rna_writer.write(db_version.as_bytes())?;
        pseudogenes_writer.write(db_version.as_bytes())?;
        all_names_writer.write(db_version.as_bytes())?;

        let rna_re = Regex::new(r"RNA").unwrap();

        for (_, gene_details) in &self.api_maps.genes {
            if gene_details.taxonid != load_org_taxonid {
                continue;
            }

            let synonyms =
                gene_details.synonyms.iter().filter(|synonym| {
                    synonym.synonym_type == "exact"
                })
                .map(|synonym| synonym.name.clone())
                .collect::<Vec<_>>()
                .join(",");

            let line = format!("{}\t{}\t{}\n",
                               gene_details.uniquename,
                               gene_details.name.clone().unwrap_or("".to_owned()),
                               synonyms);

            let line_with_product = format!("{}\t{}\t{}\t{}\n",
                                            gene_details.uniquename,
                                            gene_details.name.clone().unwrap_or("".to_owned()),
                                            synonyms,
                                            gene_details.product.clone().unwrap_or("".to_owned()));

            all_names_writer.write(line.as_bytes())?;

            if gene_details.feature_type == "pseudogene" {
                pseudogenes_writer.write(line.as_bytes())?;
            } else {
                if gene_details.feature_type == "mRNA gene" {
                    gene_writer.write(line_with_product.as_bytes())?;
                } else {
                    if rna_re.is_match(&gene_details.feature_type) {
                        rna_writer.write(line_with_product.as_bytes())?;
                    }
                }
            }
        }

        gene_writer.flush()?;
        rna_writer.flush()?;
        pseudogenes_writer.flush()?;
        all_names_writer.flush()?;

        Ok(())
    }

    fn write_protein_features(&self, output_dir: &str) -> Result<(), io::Error> {
        let peptide_stats_name = format!("{}/PeptideStats.tsv", output_dir);
        let peptide_stats_file = File::create(peptide_stats_name).expect("Unable to open file");
        let mut peptide_stats_writer = BufWriter::new(&peptide_stats_file);

        let peptide_stats_header = "Systematic_ID\tMass (kDa)\tpI\tCharge\tResidues\tCAI\n";
        peptide_stats_writer.write(peptide_stats_header.as_bytes())?;
 
        let protein_features_name = format!("{}/ProteinFeatures.tsv", output_dir);
        let protein_features_file = File::create(protein_features_name).expect("Unable to open file");
        let mut protein_features_writer = BufWriter::new(&protein_features_file);

        let protein_features_header =
            "systematic_id\tgene_name\tpeptide_id\tdomain_id\tdatabase\tseq_start\tseq_end\n";
        protein_features_writer.write(protein_features_header.as_bytes())?;

        for (gene_uniquename, gene_details) in &self.api_maps.genes {
            if let Some(transcript) = gene_details.transcripts.get(0) {
                if let Some(ref protein) = transcript.protein {
                    let line = format!("{}\t{:.2}\t{}\t{}\t{}\t{}\n",
                                       gene_uniquename, protein.molecular_weight,
                                       protein.isoelectric_point,
                                       protein.charge_at_ph7,
                                       protein.sequence.len() - 1,
                                       protein.codon_adaptation_index);
                    peptide_stats_writer.write(line.as_bytes())?;

                    let gene_name = gene_details.name.clone().unwrap_or("".to_owned());
                    for interpro_match in &gene_details.interpro_matches {
                        let line_start = format!("{}\t{}\t{}\t{}\t{}",
                                                 gene_uniquename, gene_name,
                                                 protein.uniquename, interpro_match.id,
                                                 interpro_match.dbname);
                        for location in &interpro_match.locations {
                            let line = format!("{}\t{}\t{}\n", line_start,
                                               location.start, location.end);
                            protein_features_writer.write(line.as_bytes())?;
                        }
                    }
                }
            }

        }

        peptide_stats_writer.flush()?;

        Ok(())
    }

    fn write_feature_coords(&self, config: &Config, output_dir: &str)
                            -> Result<(), io::Error>
    {
        let load_org_taxonid = config.load_organism_taxonid;

        let write_line =
            |uniquename: &str, location: &ChromosomeLocation,
             writer: &mut BufWriter<&File>| {
                let display_strand =
                    if location.strand == Strand::Forward {1} else {-1};
                let line = format!("{}\t{}\t{}\t{}\n",
                                   uniquename, location.start_pos,
                                   location.end_pos, display_strand);
                writer.write(line.as_bytes())
        };

        for (chr_uniquename, chr_details) in &self.chromosomes {

            if chr_details.taxonid != load_org_taxonid {
                continue;
            }

            let gene_file_name = format!("{}/{}.gene.coords.tsv", output_dir, chr_uniquename);
            let cds_file_name = format!("{}/{}.cds.coords.tsv", output_dir, chr_uniquename);
            let exon_file_name = format!("{}/{}.exon.coords.tsv", output_dir, chr_uniquename);

            let gene_file = File::create(gene_file_name).expect("Unable to open file");
            let cds_file = File::create(cds_file_name).expect("Unable to open file");
            let exon_file = File::create(exon_file_name).expect("Unable to open file");

            let mut gene_writer = BufWriter::new(&gene_file);
            let mut cds_writer = BufWriter::new(&cds_file);
            let mut exon_writer = BufWriter::new(&exon_file);

            for gene_uniquename in &chr_details.gene_uniquenames {
                let gene = self.api_maps.genes.get(gene_uniquename).unwrap();
                if let Some(ref gene_location) = gene.location {
                    write_line(gene_uniquename, gene_location, &mut gene_writer)?;

                    for transcript in &gene.transcripts {
                        if let Some(ref cds_location) = transcript.cds_location {
                            write_line(gene_uniquename, cds_location, &mut cds_writer)?;
                        }

                        let is_forward =
                            transcript.parts[0].location.strand == Strand::Forward;

                        if is_forward {
                            for part in &transcript.parts {
                                if part.feature_type == FeatureType::Exon {
                                    write_line(gene_uniquename, &part.location, &mut exon_writer)?;
                                }
                            }
                        } else {
                            for part in transcript.parts.iter().rev() {
                                if part.feature_type == FeatureType::Exon {
                                    write_line(gene_uniquename, &part.location, &mut exon_writer)?;
                                }
                            }
                        }
                    }
                }
            }

            gene_writer.flush()?;
            cds_writer.flush()?;
            exon_writer.flush()?;
        }

        Ok(())
    }

    pub fn write(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        let web_json_path = self.create_dir(output_dir, "web-json");

        self.write_chromosome_json(config, &web_json_path);
        println!("wrote {} chromosomes", self.get_chromosomes().len());
        self.write_gene_summaries(&web_json_path);
        println!("wrote gene summaries");
        self.write_metadata(&web_json_path);
        println!("wrote metadata");
        self.write_recent_references(&web_json_path);
        println!("wrote recent references");
        self.write_api_maps(&web_json_path);
        self.write_solr_data(&web_json_path);
        println!("wrote search data");
        self.write_subsets(&web_json_path);
        println!("wrote subsets");

        let fasta_path = self.create_dir(output_dir, "fasta");
        let feature_sequences_path = self.create_dir(&fasta_path, "feature_sequences");
        self.write_feature_sequences(&feature_sequences_path);
        let chromosomes_path = self.create_dir(&fasta_path, "chromosomes");
        self.write_chromosome_sequences(config, &chromosomes_path);

        let misc_path = self.create_dir(output_dir, "misc");
        self.write_gene_id_table(&config, &misc_path)?;
        self.write_protein_features(&misc_path)?;
        self.write_feature_coords(&config, &misc_path)?;

        Ok(())
    }

    pub fn store_jsonb(&self, conn: &Connection) {
        let trans = conn.transaction().unwrap();

        for (uniquename, gene_details) in &self.api_maps.genes {
            let serde_value = serde_json::value::to_value(&gene_details).unwrap();
            trans.execute("INSERT INTO web_json.gene (uniquename, data) values ($1, $2)",
                          &[&uniquename, &serde_value]).unwrap();
        }
        for (uniquename, ref_details) in &self.api_maps.references {
            let serde_value = serde_json::value::to_value(&ref_details).unwrap();
            trans.execute("INSERT INTO web_json.reference (uniquename, data) values ($1, $2)",
                          &[&uniquename, &serde_value]).unwrap();
        }
        for (termid, term_details) in &self.api_maps.terms {
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
