use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};

use std::fmt::Display;
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::Arc;

use regex::Regex;

use flexstr::{SharedStr as FlexStr, shared_str as flex_str, ToSharedStr, shared_fmt as flex_fmt};
use serde_with::serde_as;

pub type TypeInteractionAnnotationMap =
    HashMap<TypeName, Vec<InteractionAnnotation>>;
pub type UniquenameGeneMap =
    BTreeMap<GeneUniquename, GeneDetails>;
pub type UniquenameAlleleMap =
    BTreeMap<AlleleUniquename, AlleleDetails>;
pub type UniquenameTranscriptMap =
    HashMap<TranscriptUniquename, TranscriptDetails>;
pub type UniquenameProteinMap =
    HashMap<ProteinUniquename, ProteinDetails>;
pub type UniquenameReferenceMap =
    HashMap<ReferenceUniquename, ReferenceDetails>;

pub type UniquenameAlleleDetailsMap = BTreeMap<AlleleUniquename, AlleleDetails>;
pub type DisplayUniquenameGenotypeMap = HashMap<GenotypeDisplayUniquename, GenotypeDetails>;
pub type UniquenameFeatureShortMap = HashMap<FlexStr, FeatureShort>;
pub type TermIdDetailsMap = HashMap<TermId, TermDetails>;
pub type ChrNameDetailsMap = BTreeMap<ChromosomeName, ChromosomeDetails>;

pub type IdGenotypeMap = HashMap<GenotypeUniquename, GenotypeDetails>;
pub type IdGeneShortMap = HashMap<GeneUniquename, GeneShort>;
pub type IdRcTermShortMap = HashMap<TermId, Rc<TermShort>>;
pub type IdRcTermDetailsMap = HashMap<TermId, Rc<TermDetails>>;

pub type GeneShortMap = HashMap<GeneUniquename, GeneShort>;
pub type GenotypeShortMap = HashMap<GeneUniquename, GenotypeShort>;
pub type AlleleShortMap = HashMap<AlleleUniquename, AlleleShort>;
pub type TermShortMap = HashMap<TermId, TermShort>;

pub type OntAnnotationId = i32;
pub type IdOntAnnotationDetailMap = HashMap<OntAnnotationId, OntAnnotationDetail>;

pub type TermShortOptionMap = HashMap<TermId, Option<TermShort>>;
pub type GeneShortOptionMap = HashMap<GeneUniquename, Option<GeneShort>>;
pub type ReferenceShortOptionMap = HashMap<ReferenceUniquename, Option<ReferenceShort>>;
pub type TranscriptDetailsOptionMap =
    HashMap<TranscriptUniquename, Option<TranscriptDetails>>;

pub type GeneticInteractionMap = HashMap<GeneticInteractionKey, Vec<GeneticInteractionDetail>>;

pub type ProteinComplexMap = HashMap<ProteinComplexUniquename, ProteinComplexDetails>;

use std::rc::Rc;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

use crate::api::search_types::SolrMatchHighlight;
use crate::db::chado_queries::CommunityResponseRate;
use crate::web::config::*;
use crate::types::*;
use crate::interpro::InterProMatch;
use crate::rnacentral::RfamAnnotation;

pub trait DataLookup {
    fn get_term(&self, termid: &TermId) -> Option<Arc<TermDetails>>;

    fn get_gene(&self, gene_uniquename: &GeneUniquename) -> Option<Arc<GeneDetails>>;

    fn get_allele(&self, allele_unquename: &AlleleUniquename) -> Option<Arc<AlleleDetails>>;

    fn get_reference(&self, reference_uniquename: &ReferenceUniquename) -> Option<Arc<ReferenceDetails>>;

    fn get_genotype(&self, genotype_display_uniquename: &GenotypeDisplayUniquename)
           -> Option<Arc<GenotypeDetails>>;

    fn get_annotation_detail(&self, annotation_id: OntAnnotationId)
           -> Option<Arc<OntAnnotationDetail>>;
}

pub type RNAcentralAnnotations = HashMap<FlexStr, Vec<RfamAnnotation>>;

fn empty_string()-> FlexStr {
    flex_str!("")
}

fn is_default<T: Default + PartialEq>(t: &T) -> bool {
    t == &T::default()
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Debug, Clone)]
pub enum Ploidiness {
#[serde(rename = "haploid")]
    Haploid,
#[serde(rename = "diploid")]
    Diploid,
#[serde(rename = "any")]
    Any,
}

#[derive(Serialize, Deserialize, Clone, Debug, Hash, PartialEq, Eq)]
pub enum Throughput {
#[serde(rename = "high")]
    HighThroughput,
#[serde(rename = "low")]
    LowThroughput,
#[serde(rename = "non-experimental")]
    NonExperimental,
}

impl Display for Throughput {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = match self {
            Self::LowThroughput => "low",
            Self::HighThroughput => "high",
            Self::NonExperimental => "non-experimental",
        };

        f.write_str(s)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, PartialOrd, Hash)]
pub struct GeneAndGeneProduct {
    pub gene_uniquename: GeneUniquename,
    pub product: TermId, // eg.  PR:000027705
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, PartialOrd, Hash)]
pub enum ExtRange {
#[serde(rename = "gene_uniquename")]
    Gene(GeneUniquename),
#[serde(rename = "transcript_uniquename")]
    Transcript(TranscriptUniquename),
#[serde(rename = "promoter_gene_uniquename")]
    Promoter(FlexStr),
#[serde(rename = "summary_gene_uniquenames")]
    // the inner Vec length will be > 1 for cases like "binds abc1 and def2, cdc2"
    SummaryGenes(Vec<Vec<FlexStr>>),
#[serde(rename = "summary_transcript_uniquenames")]
    SummaryTranscripts(Vec<Vec<TranscriptUniquename>>),
#[serde(rename = "termid")]
    Term(TermId),
#[serde(rename = "summary_termids")]
    // See: merge_ext_part_ranges()
    SummaryTerms(Vec<TermId>),
#[serde(rename = "misc")]
    Misc(FlexStr),
#[serde(rename = "domain")]
    Domain(FlexStr),
#[serde(rename = "gene_product")]
    GeneProduct(TermId),  // eg.  PR:000027705
#[serde(rename = "gene_and_gene_product")]
// see: https://github.com/pombase/website/issues/1625
    GeneAndGeneProduct(GeneAndGeneProduct),
#[serde(rename = "summary_residues")]
    ModifiedResidues(Vec<Residue>),
}

impl ExtRange {
    pub fn is_gene(&self) -> bool {
        matches!(*self, ExtRange::Gene(_))
    }
}

impl fmt::Display for ExtRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExtRange::Gene(ref gene_uniquename) | ExtRange::Promoter(ref gene_uniquename) =>
                write!(f, "{}", gene_uniquename),
            ExtRange::Transcript(ref transcript_uniquename) =>
                write!(f, "{}", transcript_uniquename),
            ExtRange::SummaryGenes(_) => panic!("can't handle SummaryGenes\n"),
            ExtRange::SummaryTranscripts(_) => panic!("can't handle SummaryTranscripts\n"),
            ExtRange::Term(ref termid) => write!(f, "{}", termid),
            ExtRange::ModifiedResidues(ref residue) =>
                write!(f, "{}", residue.iter().map(FlexStr::to_string).collect::<Vec<_>>().join(",")),
            ExtRange::SummaryTerms(_) => panic!("can't handle SummaryGenes\n"),
            ExtRange::Misc(ref misc) => write!(f, "{}", misc),
            ExtRange::Domain(ref domain) => write!(f, "{}", domain),
            ExtRange::GeneProduct(ref gene_product) => write!(f, "{}", gene_product),
            ExtRange::GeneAndGeneProduct(ref gene_and_gene_product) =>
                write!(f, "{} ({})", gene_and_gene_product.gene_uniquename,
                       gene_and_gene_product.product),
        }
    }
}

// A single part of an extension.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ExtPart {
    pub rel_type_name: FlexStr,
    pub rel_type_display_name: FlexStr,
    pub rel_type_id: Option<TermId>,
    pub ext_range: ExtRange,
}
impl PartialEq for ExtPart {
    fn eq(&self, other: &Self) -> bool {
        self.rel_type_name == other.rel_type_name &&
            self.ext_range == other.ext_range
    }
}
impl Eq for ExtPart {
}
impl Hash for ExtPart {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.rel_type_name.hash(state);
        self.ext_range.hash(state);
    }
}

/* needs work and isn't currently needed:
impl fmt::Display for ExtPart {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
         write!(f, "{}/{} ({})", self.rel_type_display_name,
                self.rel_type_name, self.ext_range)
    }
}
*/


// used to parse Canto format extension stored as strings in prop values
#[derive(Deserialize, Debug)]
#[serde(rename_all = "camelCase")]
pub struct CantoExtPart {
    pub range_display_name: Option<FlexStr>,
    pub range_type: Option<FlexStr>,
    pub range_value: FlexStr,
    pub relation: FlexStr,
}

// minimal information about a gene used in other objects
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneShort {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,
    #[serde(skip_serializing_if="is_one", default = "one")]
    pub transcript_count: usize,

    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub flags: HashSet<FlexStr>,
}

impl GeneShort {
    pub fn from_gene_details(gene_details: &GeneDetails) -> Self {
        GeneShort {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
            transcript_count: gene_details.transcripts.len(),
            flags: gene_details.flags.clone(),
        }
    }

    pub fn display_name(&self) -> String {
        if let Some(ref name) = self.name {
            format!("{} ({})", name, self.uniquename)
        } else {
            self.uniquename.to_string()
        }
    }
}

impl From<&GeneDetails> for GeneShort {
    fn from(details: &GeneDetails) -> GeneShort {
        GeneShort {
            uniquename: details.uniquename.clone(),
            name: details.name.clone(),
            product: details.product.clone(),
            transcript_count: details.transcripts.len(),
            flags: details.flags.clone(),
        }
    }
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
        } else if other.name.is_some() {
            Ordering::Greater
        } else {
            self.uniquename.cmp(&other.uniquename)
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
    pub identifier: FlexStr,
    pub taxonid: u32,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct IdNameAndOrganism {
    pub identifier: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub secondary_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub taxonid: u32,
}

// used in serialisation
#[allow(clippy::trivially_copy_pass_by_ref)]
fn is_one(num: &usize) -> bool {
    *num == 1
}

fn one() -> usize {
    1
}

#[allow(clippy::trivially_copy_pass_by_ref)]
fn is_true(b: &bool) -> bool {
    *b
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
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub secondary_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub orthologs: Vec<IdNameAndOrganism>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    #[serde(skip_serializing_if="is_one")]
    pub transcript_count: usize,
    pub feature_type: FlexStr,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermXref {
    pub xref_id: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub xref_display_name: Option<FlexStr>,
}

// minimal information about a terms used in other objects
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermShort {
    pub name: TermName,
    pub cv_name: FlexStr,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parent_ids: HashSet<FlexStr>,
    pub termid: TermId,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub secondary_identifiers: HashSet<TermId>,
    #[serde(skip_serializing_if="is_true", default)]
    pub is_obsolete: bool,
    pub gene_count: usize,
    pub genotype_count: usize,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub xrefs: HashMap<FlexStr, TermXref>,
}

impl TermShort {
    pub fn from_term_details(term_details: &TermDetails) -> Self {
        TermShort {
            name: term_details.name.clone(),
            cv_name: term_details.cv_name.clone(),
            interesting_parent_ids: term_details.interesting_parent_ids.clone(),
            termid: term_details.termid.clone(),
            is_obsolete: term_details.is_obsolete,
            gene_count: term_details.gene_count,
            genotype_count: term_details.genotype_count,
            xrefs: term_details.xrefs.clone(),
            secondary_identifiers: term_details.secondary_identifiers.clone(),
        }
    }
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

impl From<&TermDetails> for TermShort {
    fn from(details: &TermDetails) -> TermShort {
        TermShort {
            name: details.name.clone(),
            termid: details.termid.clone(),
            cv_name: details.cv_name.clone(),
            secondary_identifiers: details.secondary_identifiers.clone(),
            interesting_parent_ids: details.interesting_parent_ids.clone(),
            is_obsolete: details.is_obsolete,
            gene_count: details.gene_count,
            genotype_count: details.genotype_count,
            xrefs: details.xrefs.clone(),
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChromosomeDetails {
    pub name: FlexStr,
    pub residues: FlexStr,
    pub ena_identifier: FlexStr,
    pub gene_uniquenames: Vec<FlexStr>,
    pub taxonid: OrganismTaxonId,
    pub gene_count: usize,
    pub coding_gene_count: usize,
}

impl ChromosomeDetails {
    pub fn make_chromosome_short(&self) -> ChromosomeShort {
        ChromosomeShort {
            name: self.name.clone(),
            length: self.residues.len(),
            ena_identifier: self.ena_identifier.clone(),
            gene_count: self.gene_count,
            coding_gene_count: self.coding_gene_count,
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ReferenceShort {
    pub uniquename: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub approved_date: Option<FlexStr>,
    pub gene_count: usize,
    pub genotype_count: usize,
}

impl ReferenceShort {
    pub fn from_reference_details(reference_details: &ReferenceDetails) -> ReferenceShort {
        ReferenceShort {
            uniquename: reference_details.uniquename.clone(),
            title: reference_details.title.clone(),
            citation: reference_details.citation.clone(),
            publication_year: reference_details.publication_year.clone(),
            authors_abbrev: reference_details.authors_abbrev.clone(),
            approved_date: reference_details.approved_date.clone(),
            gene_count: reference_details.gene_count,
            genotype_count: reference_details.genotypes_by_uniquename.keys().len(),
        }
    }
}

#[derive(PartialEq)]
pub enum ContainerType {
    Gene,
    Term,
    Reference,
    Genotype,
}

pub trait Container {
    fn container_type(&self) -> ContainerType;
}

pub trait AnnotationContainer: Container {
    fn cv_annotations(&self) -> &OntAnnotationMap;
    fn cv_annotations_mut(&mut self) -> &mut OntAnnotationMap;
    fn annotation_details(&self) -> &IdOntAnnotationDetailMap;
    fn terms_by_termid(&self) -> &TermShortOptionMap;
    fn genes_by_uniquename(&self) -> &GeneShortOptionMap;
    fn genotypes_by_uniquename(&self) -> &HashMap<GenotypeUniquename, GenotypeShort>;

    fn annotation_count(&self) -> usize {
        let mut annotation_ids = HashSet::new();

        for term_annotations in self.cv_annotations().values() {
            for term_annotation in term_annotations {
                for annotation_detail_id in &term_annotation.annotations {
                    annotation_ids.insert(annotation_detail_id);
                }
            }
        }

        annotation_ids.len()
    }
}

pub trait OrthologAnnotationContainer: AnnotationContainer {
    fn ortholog_annotations(&self) -> &Vec<OrthologAnnotation>;
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct AnnotationCurator {
    pub name: FlexStr,
    pub community_curator: bool,
    pub annotation_count: usize,
    pub orcid: Option<FlexStr>,
    pub file_type: Option<FlexStr>,
    pub file_name: Option<FlexStr>,
}

#[serde_as]
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ReferenceDetails {
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

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub pubmed_keyword_genes: Vec<FlexStr>,

    // count of genes from the main organism of the site (ie. only pombe)
    pub gene_count: usize,

    // count of genes annotated in LTP experiments
    pub ltp_gene_count: usize,

    // This is set to the year part of canto_first_approved_date if it is
    // not None, otherwise set to the year part of canto_approved_date, otherwise
    // canto_session_submitted_date
    #[serde(skip_serializing_if="Option::is_none")]
    pub approved_date: Option<FlexStr>,
    pub cv_annotations: OntAnnotationMap,
    pub physical_interactions: Vec<InteractionAnnotation>,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    #[serde_as(as = "Vec<(_, _)>")]
    pub genetic_interactions: GeneticInteractionMap,

    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
    pub genes_by_uniquename: GeneShortOptionMap,
    pub genotypes_by_uniquename: HashMap<GenotypeUniquename, GenotypeShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub references_by_uniquename: ReferenceShortOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub transcripts_by_uniquename: TranscriptDetailsOptionMap,
    pub terms_by_termid: TermShortOptionMap,
    pub annotation_details: IdOntAnnotationDetailMap,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub pdb_entries: Vec<PDBEntry>,
}

impl Container for ReferenceDetails {
    fn container_type(&self) -> ContainerType {
        ContainerType::Reference
    }
}

impl AnnotationContainer for ReferenceDetails {
    fn cv_annotations(&self) -> &OntAnnotationMap {
        &self.cv_annotations
    }
    fn cv_annotations_mut(&mut self) -> &mut OntAnnotationMap {
        &mut self.cv_annotations
    }
    fn annotation_details(&self) -> &IdOntAnnotationDetailMap {
        &self.annotation_details
    }
    fn terms_by_termid(&self) -> &TermShortOptionMap {
        &self.terms_by_termid
    }
    fn genes_by_uniquename(&self) -> &GeneShortOptionMap {
        &self.genes_by_uniquename
    }
    fn genotypes_by_uniquename(&self) -> &HashMap<GenotypeUniquename, GenotypeShort> {
        &self.genotypes_by_uniquename
    }
}

impl OrthologAnnotationContainer for ReferenceDetails {
    fn ortholog_annotations(&self) -> &Vec<OrthologAnnotation> {
        &self.ortholog_annotations
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct IdentifierAndName {
    pub identifier: FlexStr,
    pub name: FlexStr,
}

// the GO with/from
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum WithFromValue {
    Gene(GeneShort),
    Transcript(TranscriptUniquename),
    Term(TermShort),
    Identifier(FlexStr),
    IdentifierAndName(IdentifierAndName),
}

impl WithFromValue {
    pub fn id(&self) -> FlexStr {
        match self {
            WithFromValue::Gene(gene) => gene.uniquename.clone(),
            WithFromValue::Term(term) => term.termid.clone(),
            WithFromValue::Transcript(transcript_uniquename) => transcript_uniquename.clone(),
            WithFromValue::Identifier(id) => id.clone(),
            WithFromValue::IdentifierAndName(id_and_name) =>
                id_and_name.identifier.clone(),
        }
    }
}

impl From<WithFromValue> for FlexStr {
    fn from(with_from: WithFromValue) -> Self {
        match with_from {
            WithFromValue::Gene(gene) => gene.uniquename,
            WithFromValue::Transcript(transcript_uniquename) => transcript_uniquename,
            WithFromValue::Term(term) => term.termid,
            WithFromValue::Identifier(id) => id,
            WithFromValue::IdentifierAndName(id_and_name) =>
                flex_fmt!("{} ({})", id_and_name.name, id_and_name.identifier),
        }
    }
}

impl Ord for WithFromValue {
    fn cmp(& self, other: & WithFromValue) -> Ordering {
        let self_string =
            match self {
                WithFromValue::Gene(ref gene) => &gene.uniquename,
                WithFromValue::Transcript(ref transcript_uniquename) => transcript_uniquename,
                WithFromValue::Term(ref term) => &term.name,
                WithFromValue::Identifier(ref id) => id,
                WithFromValue::IdentifierAndName(ref id_and_name) => &id_and_name.identifier,
            };
        let other_string =
            match other {
                WithFromValue::Gene(ref gene) => &gene.uniquename,
                WithFromValue::Transcript(ref transcript_uniquename) => transcript_uniquename,
                WithFromValue::Term(ref term) => &term.name,
                WithFromValue::Identifier(ref id) => id,
                WithFromValue::IdentifierAndName(ref id_and_name) => &id_and_name.identifier,
            };

        self_string.cmp(other_string)
    }
}
impl PartialOrd for WithFromValue {
    fn partial_cmp(&self, other: &WithFromValue) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub type CuratorOrcid = FlexStr;

pub type Promoter = FlexStr;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct AnnotationPromoter {
    pub allele_uniquename: FlexStr,
    pub allele_display_name: FlexStr,
    pub allele_expression: Option<FlexStr>,
    pub allele_gene: GeneShort,
    pub promoter: Promoter,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct OntAnnotationDetail {
    pub id: i32,
    pub genes: Vec<GeneUniquename>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcript_uniquenames: Vec<TranscriptUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference: Option<ReferenceUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub eco_evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub annotation_phenotype_score: Option<FlexStr>,
    pub extension: Vec<ExtPart>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub withs: HashSet<WithFromValue>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub froms: HashSet<WithFromValue>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub residue: Option<Residue>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_product_form_id: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub qualifiers: Vec<Qualifier>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_ex_props: Option<GeneExProps>,
    // only for genotype/phenotype annotation:
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype: Option<GenotypeUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_background: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub allele_promoters: Vec<AnnotationPromoter>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub conditions: HashSet<TermId>,
    #[serde(skip_serializing_if="BTreeSet::is_empty", default)]
    pub condition_details: BTreeSet<(TermId, Option<String>)>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub date: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub assigned_by: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub throughput: Option<Throughput>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub curator: Option<CuratorOrcid>,
}

impl PartialEq for OntAnnotationDetail {
    fn eq(&self, other: &OntAnnotationDetail) -> bool {
        self.genes == other.genes &&
            self.reference == other.reference &&
            self.evidence == other.evidence &&
            self.extension == other.extension &&
            self.withs == other.withs &&
            self.froms == other.froms &&
            self.residue == other.residue &&
            self.gene_product_form_id == other.gene_product_form_id &&
            self.qualifiers == other.qualifiers &&
            self.gene_ex_props == other.gene_ex_props &&
            self.genotype == other.genotype &&
            self.genotype_background == other.genotype_background &&
            self.conditions == other.conditions &&
            self.date == other.date &&
            self.assigned_by == other.assigned_by &&
            self.throughput == other.throughput
    }
}
impl Eq for OntAnnotationDetail { }

fn hash_a_hashset<T: Hash+Ord, H: Hasher>(hash_set: &HashSet<T>, state: &mut H) {
    let mut v: Vec<_> = hash_set.iter().collect();
    v.sort();
    for el in v.iter() {
        (*el).hash(state);
    }
}

impl Hash for OntAnnotationDetail {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.genes.hash(state);
        self.reference.hash(state);
        self.evidence.hash(state);
        self.extension.hash(state);
        hash_a_hashset(&self.withs, state);
        hash_a_hashset(&self.froms, state);
        self.residue.hash(state);
        self.gene_product_form_id.hash(state);
        self.qualifiers.hash(state);
        self.gene_ex_props.hash(state);
        self.genotype.hash(state);
        self.genotype_background.hash(state);
        hash_a_hashset(&self.conditions, state);
        self.date.hash(state);
        self.assigned_by.hash(state);
        self.throughput.hash(state);
    }
}

impl PartialOrd for OntAnnotationDetail {
    fn partial_cmp(&self, other: &OntAnnotationDetail) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OntAnnotationDetail {
    fn cmp(&self, other: &OntAnnotationDetail) -> Ordering {
        if self == other {
            Ordering::Equal
        } else {
            self.id.cmp(&other.id)
        }
    }
}
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct OntTermAnnotations {
    pub term: TermId,
    pub is_not: bool,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub rel_names: HashSet<RelName>,
    pub annotations: Vec<OntAnnotationId>,
    pub summary: Option<Vec<TermSummaryRow>>,
}

impl PartialEq for OntTermAnnotations {
    fn eq(&self, other: &OntTermAnnotations) -> bool {
        self.term == other.term
    }
}
impl Eq for OntTermAnnotations { }
impl Hash for OntTermAnnotations {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.term.hash(state);
        self.is_not.hash(state);
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct OntAnnotation {
    pub term_short: TermShort,
    pub id: i32,
    pub genes: HashSet<GeneShort>,
    pub reference_short: Option<ReferenceShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    pub extension: Vec<ExtPart>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub withs: HashSet<WithFromValue>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub froms: HashSet<WithFromValue>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub residue: Option<Residue>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub qualifiers: Vec<Qualifier>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_ex_props: Option<GeneExProps>,
    // only for genotype/phenotype annotation:
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_short: Option<GenotypeShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_background: Option<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub conditions: HashSet<TermShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub assigned_by: Option<FlexStr>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermSummaryRow {
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gene_uniquenames: Vec<GeneUniquename>, // for term and ref pages
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub genotype_uniquenames: Vec<GenotypeUniquename>, // for term pages
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub extension: Vec<ExtPart>,
}
impl PartialEq for TermSummaryRow {
    fn eq(&self, other: &TermSummaryRow) -> bool {
        self.gene_uniquenames == other.gene_uniquenames &&
            self.genotype_uniquenames == other.genotype_uniquenames &&
            self.extension == other.extension
    }
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
    pub show_in_summary: bool,
    pub ontology_name: FlexStr,
    pub ext_rel_display_name: FlexStr,
    pub gene: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none", default)]
    pub genotype_uniquename: Option<GenotypeUniquename>,
    #[serde(skip_serializing_if="Option::is_none", default)]
    pub reference_uniquename: Option<ReferenceUniquename>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SynonymDetails {
    pub name: FlexStr,
    #[serde(rename = "type")]
    pub synonym_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none", default)]
    pub reference: Option<ReferenceUniquename>,
}

impl PartialEq for SynonymDetails {
    fn eq(&self, other: &SynonymDetails) -> bool {
        self.name == other.name && self.synonym_type == other.synonym_type
    }
}
impl Eq for SynonymDetails { }
impl Ord for SynonymDetails {
    fn cmp(&self, other: &SynonymDetails) -> Ordering {
        let type_cmp = self.synonym_type.cmp(&other.synonym_type);
        if type_cmp == Ordering::Equal {
            self.name.cmp(&other.name)
        } else {
            type_cmp
        }
    }
}
impl PartialOrd for SynonymDetails {
    fn partial_cmp(&self, other: &SynonymDetails) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Hash for SynonymDetails {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
        self.synonym_type.hash(state);
    }
}


#[derive(Serialize, Deserialize, Copy, Clone, Debug, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward = 1,
    Reverse = -1,
    Unstranded = 0,
}

impl Strand {
    pub fn to_gff_str(&self) -> &'static str {
        match *self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unstranded => ".",
        }
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.to_gff_str())
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChromosomeShort {
    pub name: FlexStr,
    pub length: usize,
    pub ena_identifier: FlexStr,
    pub gene_count: usize,
    pub coding_gene_count: usize,
}

#[derive(Serialize, Deserialize, PartialEq, Clone, Debug)]
pub enum Phase {
#[serde(rename = "0")]
    Zero,
#[serde(rename = "1")]
    One,
#[serde(rename = "2")]
    Two,
}

impl Phase {
    pub fn to_gff_str(&self) -> &'static str {
        match *self {
            Phase::Zero => "0",
            Phase::One => "1",
            Phase::Two => "2",
        }
    }
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChromosomeLocation {
    pub chromosome_name: FlexStr,
    pub start_pos: usize,
    pub end_pos: usize,
    pub strand: Strand,
    #[serde(skip_serializing_if="Option::is_none")]
    pub phase: Option<Phase>,
}

impl ChromosomeLocation {
    pub fn len(&self) -> usize {
        (self.end_pos + 1).saturating_sub(self.start_pos)
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum DeletionViability {
    Viable,
    Inviable,
    DependsOnConditions,
    Unknown,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum PresentAbsent {
    Present,
    Absent,
    NotApplicable,
    Unknown,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneHistoryEntry {
    pub previous_coords: String,
    pub date: FlexStr,
    pub references: Vec<FlexStr>,
    pub comments: Option<String>,
    pub genome_snapshot_link: Option<String>,
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct PDBGeneChain {
    pub gene_uniquename: GeneUniquename,
    pub chain: String,
    pub position: String,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct PDBEntry {
    pub pdb_id: PdbId,
    pub gene_chains: Vec<PDBGeneChain>,
    pub title: String,
    pub entry_authors: String,
    pub entry_authors_abbrev: String,
    pub reference_uniquename: Option<ReferenceUniquename>,
    pub experimental_method: String,
    pub resolution: String,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ReferenceAndSource {
    pub reference_uniquename: FlexStr,
    pub source: FlexStr,
}

#[serde_as]
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneDetails {
    pub uniquename: GeneUniquename,
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
    pub biogrid_interactor_id: Option<u32>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub rnacentral_urs_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub interpro_matches: Vec<InterProMatch>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub tm_domain_coords: Vec<(usize, usize)>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub disordered_region_coords: Vec<(usize, usize)>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub low_complexity_region_coords: Vec<(usize, usize)>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub coiled_coil_coords: Vec<(usize, usize)>,
    pub has_protein_features: bool,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub rfam_annotations: Vec<RfamAnnotation>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub pdb_entries: Vec<PDBEntry>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub orfeome_identifier: Option<FlexStr>,
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
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gene_neighbourhood: Vec<GeneShort>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    // transcripts of this gene
    pub transcripts: Vec<TranscriptUniquename>,
    pub cv_annotations: OntAnnotationMap,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub physical_interactions: Vec<InteractionAnnotation>,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    #[serde_as(as = "Vec<(_, _)>")]
    pub genetic_interactions: GeneticInteractionMap,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub paralog_annotations: Vec<ParalogAnnotation>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub target_of_annotations: Vec<TargetOfAnnotation>,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub references_by_uniquename: ReferenceShortOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub genes_by_uniquename: GeneShortOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub genotypes_by_uniquename: HashMap<GenotypeUniquename, GenotypeShort>,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    // transcripts referenced by this gene
    pub transcripts_by_uniquename: TranscriptDetailsOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub terms_by_termid: TermShortOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub annotation_details: IdOntAnnotationDetailMap,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub feature_publications: HashSet<ReferenceAndSource>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    // A Vec of the term IDs of subsets for this gene.  Any useful subset
    // that contains any term for any annotation in the gene is included.
    // "useful" means that the front end might need it, eg. slim term IDs
    pub subset_termids: HashSet<TermId>,

    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub gocams: HashSet<GoCamIdAndTitle>,

    #[serde(skip_serializing_if="Option::is_none")]
    pub rnacentral_2d_structure_id: Option<RnaUrsId>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub gene_history: Vec<GeneHistoryEntry>,
}

impl GeneDetails {
    pub fn display_name(&self) -> FlexStr {
        if let Some(ref name) = self.name {
            flex_fmt!("{} ({})", name, self.uniquename)
        } else {
            self.uniquename.clone()
        }
    }
}

impl PartialEq for GeneDetails {
    fn eq(&self, other: &GeneDetails) -> bool {
        self.uniquename == other.uniquename
    }
}
impl Eq for GeneDetails { }
impl Ord for GeneDetails {
    fn cmp(&self, other: &GeneDetails) -> Ordering {
        if self.name.is_some() {
            if other.name.is_some() {
                self.name.cmp(&other.name)
            } else { Ordering::Less }
        } else if other.name.is_some() {
            Ordering::Greater
        } else { self.uniquename.cmp(&other.uniquename) }
    }
}
impl PartialOrd for GeneDetails {
    fn partial_cmp(&self, other: &GeneDetails) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Hash for GeneDetails {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.uniquename.hash(state);
    }
}

impl Container for GeneDetails {
    fn container_type(&self) -> ContainerType {
        ContainerType::Gene
    }
}

impl AnnotationContainer for GeneDetails {
    fn cv_annotations(&self) -> &OntAnnotationMap {
        &self.cv_annotations
    }
    fn cv_annotations_mut(&mut self) -> &mut OntAnnotationMap {
        &mut self.cv_annotations
    }
    fn annotation_details(&self) -> &IdOntAnnotationDetailMap {
        &self.annotation_details
    }
    fn terms_by_termid(&self) -> &TermShortOptionMap {
        &self.terms_by_termid
    }
    fn genes_by_uniquename(&self) -> &GeneShortOptionMap {
        &self.genes_by_uniquename
    }
    fn genotypes_by_uniquename(&self) -> &HashMap<GenotypeUniquename, GenotypeShort> {
        &self.genotypes_by_uniquename
    }
}

impl OrthologAnnotationContainer for GeneDetails {
    fn ortholog_annotations(&self) -> &Vec<OrthologAnnotation> {
        &self.ortholog_annotations
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinDetails {
    pub uniquename: ProteinUniquename,
    pub sequence: FlexStr,
    pub number_of_residues: usize,  // residue count not including stop codon
    pub product: Option<FlexStr>,
    pub molecular_weight: f32,
    pub average_residue_weight: f32,
    pub charge_at_ph7: f32,
    pub isoelectric_point: f32,
    pub codon_adaptation_index: f32,
}

impl ProteinDetails {
    pub fn sequence_length(&self) -> usize {
        if self.sequence.ends_with("*") {
            self.sequence.len() - 1
        } else {
            self.sequence.len()
        }
    }
}

pub type Residues = FlexStr;

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
#[serde(rename = "dg_repeat")]
    DGRepeat,
#[serde(rename = "dh_repeat")]
    DHRepeat,
#[serde(rename = "gap")]
    Gap,
#[serde(rename = "gene_group")]
    GeneGroup,
#[serde(rename = "long_terminal_repeat")]
    LongTerminalRepeat,
#[serde(rename = "low_complexity_region")]
    LowComplexityRegion,
#[serde(rename = "LTR_retrotransposon")]
    LTRRetrotransposon,
#[serde(rename = "mating_type_region")]
    MatingTypeRegion,
#[serde(rename = "nuclear_mt_pseudogene")]
    NuclearMtPseudogene,
#[serde(rename = "origin_of_replication")]
    OriginOfReplication,
#[serde(rename = "polyA_signal_sequence")]
    PolyASignalSequence,
#[serde(rename = "polyA_site")]
    PolyASite,
#[serde(rename = "promoter")]
    Promoter,
#[serde(rename = "region")]
    Region,
#[serde(rename = "regional_centromere")]
    RegionalCentromere,
#[serde(rename = "regional_centromere_central_core")]
    RegionalCentromereCentralCore,
#[serde(rename = "regional_centromere_inner_repeat_region")]
    RegionalCentromereInnerRepeatRegion,
#[serde(rename = "repeat_region")]
    RepeatRegion,
#[serde(rename = "TR_box")]
    TRBox,
#[serde(rename = "antisense_RNA")]
    AntisenseRNA,
#[serde(rename = "sncRNA")]
    SncRNA,
#[serde(rename = "lncRNA")]
    LncRNA,
#[serde(rename = "guide_RNA")]
    GuideRNA,
    SNP,
}

impl FeatureType {
    pub fn is_any_intron_type(&self) -> bool {
        *self == FeatureType::CdsIntron || *self == FeatureType::FivePrimeUtrIntron ||
        *self == FeatureType::ThreePrimeUtrIntron
    }
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
            FeatureType::DGRepeat => "dg_repeat",
            FeatureType::DHRepeat => "dh_repeat",
            FeatureType::Gap => "gap",
            FeatureType::GeneGroup => "gene_group",
            FeatureType::LongTerminalRepeat => "long_terminal_repeat",
            FeatureType::LowComplexityRegion => "low_complexity_region",
            FeatureType::LTRRetrotransposon => "LTR_retrotransposon",
            FeatureType::MatingTypeRegion => "mating_type_region",
            FeatureType::NuclearMtPseudogene => "nuclear_mt_pseudogene",
            FeatureType::OriginOfReplication => "origin_of_replication",
            FeatureType::PolyASignalSequence => "polyA_signal_sequence",
            FeatureType::PolyASite => "polyA_site",
            FeatureType::Promoter => "promoter",
            FeatureType::Region => "region",
            FeatureType::RegionalCentromere => "regional_centromere",
            FeatureType::RegionalCentromereCentralCore => "regional_centromere_central_core",
            FeatureType::RegionalCentromereInnerRepeatRegion => "regional_centromere_inner_repeat_region",
            FeatureType::RepeatRegion => "repeat_region",
            FeatureType::TRBox => "TR_box",
            FeatureType::AntisenseRNA => "antisense_RNA",
            FeatureType::SncRNA => "sncRNA",
            FeatureType::LncRNA => "lncRNA",
            FeatureType::GuideRNA => "guide_RNA",
            FeatureType::SNP => "SNP",
        })
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct FeatureShort {
    pub feature_type: FeatureType,
    pub uniquename: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub location: ChromosomeLocation,
    pub residues: Residues,
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TranscriptDetails {
    pub uniquename: TranscriptUniquename,
    pub name: Option<GeneName>,
    pub location: ChromosomeLocation,
    pub parts: Vec<FeatureShort>,
    pub transcript_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub protein: Option<ProteinDetails>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<ChromosomeLocation>,
    pub gene_uniquename: FlexStr,
    // the CDS length or RNA length without introns - sum of lengths of exons
    pub rna_seq_length_spliced: Option<NonZeroUsize>,
    // the CDS length (protein coding) or RNA length (non-coding) including introns
    pub rna_seq_length_unspliced: Option<NonZeroUsize>,
}

impl TranscriptDetails {
    pub fn spliced_transcript_sequence(&self) -> FlexStr {
        let mut seq = String::new();

        for part in &self.parts {
            if part.feature_type == FeatureType::Exon {
                seq += &part.residues;
            }
        }

        seq.to_shared_str()
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct GenotypeLocus {
    pub expressed_alleles: Vec<ExpressedAllele>,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct GenotypeShort {
    pub display_uniquename: GenotypeDisplayUniquename,
    pub display_name: GenotypeDisplayName,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub loci: Vec<GenotypeLocus>,
}

impl GenotypeShort {
    pub fn ploidiness(&self) -> Ploidiness {
        for locus in &self.loci {
            if locus.expressed_alleles.len() > 1 {
                return Ploidiness::Diploid;
            }
        }
        Ploidiness::Haploid
    }
}

#[serde_as]
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GenotypeDetails {
    pub display_uniquename: GenotypeUniquename,
    pub display_name: GenotypeDisplayName,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub taxonid: u32,
    pub loci: Vec<GenotypeLocus>,
    pub comment: Option<FlexStr>,
    pub ploidiness: Ploidiness,
    pub cv_annotations: OntAnnotationMap,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    #[serde_as(as = "Vec<(_, _)>")]
    pub double_mutant_genetic_interactions: GeneticInteractionMap,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    #[serde_as(as = "Vec<(_, _)>")]
    pub rescue_genetic_interactions: GeneticInteractionMap,

    pub references_by_uniquename: ReferenceShortOptionMap,
    pub genes_by_uniquename: GeneShortOptionMap,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    // transcripts referenced by this genotype
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub transcripts_by_uniquename: TranscriptDetailsOptionMap,
    pub genotypes_by_uniquename: HashMap<GenotypeUniquename, GenotypeShort>,
    pub terms_by_termid: TermShortOptionMap,
    pub annotation_details: IdOntAnnotationDetailMap,
    pub annotation_count: usize,
}

impl GenotypeDetails {
    pub fn ploidiness(&self) -> Ploidiness {
        for locus in &self.loci {
            if locus.expressed_alleles.len() > 1 {
                return Ploidiness::Diploid;
            }
        }
        Ploidiness::Haploid
    }

    pub fn get_gene_short(&self, gene_uniquename: &GeneUniquename) -> Option<GeneShort> {
        self.genes_by_uniquename().get(gene_uniquename)
           .map(|opt_gene| { opt_gene.to_owned().unwrap() })
    }

    pub fn get_allele(&self, allele_uniquename: &AlleleUniquename) -> Option<&AlleleShort> {
        self.alleles_by_uniquename.get(allele_uniquename)
    }
}

impl From<&GenotypeDetails> for GenotypeShort {
    fn from(details: &GenotypeDetails) -> GenotypeShort {
        GenotypeShort {
            display_uniquename: details.display_uniquename.clone(),
            display_name: details.display_name.clone(),
            name: details.name.clone(),
            loci: details.loci.clone(),
        }
    }
}

impl Container for GenotypeDetails {
    fn container_type(&self) -> ContainerType {
        ContainerType::Genotype
    }
}

impl AnnotationContainer for GenotypeDetails {
    fn cv_annotations(&self) -> &OntAnnotationMap {
        &self.cv_annotations
    }
    fn cv_annotations_mut(&mut self) -> &mut OntAnnotationMap {
        &mut self.cv_annotations
    }
    fn annotation_details(&self) -> &IdOntAnnotationDetailMap {
        &self.annotation_details
    }
    fn terms_by_termid(&self) -> &TermShortOptionMap {
        &self.terms_by_termid
    }
    fn genes_by_uniquename(&self) -> &GeneShortOptionMap {
        &self.genes_by_uniquename
    }
    fn genotypes_by_uniquename(&self) -> &HashMap<GenotypeUniquename, GenotypeShort> {
        &self.genotypes_by_uniquename
    }
}

type Expression = FlexStr;

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct ExpressedAllele {
    #[serde(skip_serializing_if="Option::is_none")]
    pub expression: Option<Expression>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub promoter_gene: Option<FlexStr>,
    pub allele_uniquename: AlleleUniquename,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct AlleleShort {
    pub uniquename: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub allele_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub description: Option<FlexStr>,
    pub gene_uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<SynonymDetails>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub comments: Vec<CommentAndReference>,
}

 lazy_static! {
     static ref MUTATION_DESC_RE: Regex = Regex::new(r"\w+\d+\w+$").unwrap();
     static ref ENDS_WITH_INTEGER_RE: Regex = Regex::new(r"\d+$").unwrap();
 }

fn description_with_residue_type(allele: &AlleleShort) -> FlexStr {
    let description = allele.description.clone();
    let allele_type = &allele.allele_type;

    if let Some(description) = description {

        if description.is_empty() {
            return allele_type.clone();
        }

        if allele_type.ends_with("mutation") && MUTATION_DESC_RE.is_match(description.as_str()) {
            if allele_type.contains("amino_acid") {
                return description + " aa";
            } else if allele_type.contains("nucleotide") {
                return description + " nt";
            }
        }

        description
    } else {
        allele_type.clone()
    }
}

impl AlleleShort {
   fn display_name_helper(&self, compact: bool) -> FlexStr {
        let name = self.name.clone().unwrap_or_else(|| flex_str!("unnamed"));
        let allele_type = &self.allele_type;
        let mut description =
            self.description.as_ref().unwrap_or(allele_type).clone();

        if allele_type == "deletion" &&
            (name.ends_with('Δ') || name.ends_with("delta")) ||
            allele_type == "wild_type" && name.ends_with('+') {
                if name.ends_with(description.as_str()) || allele_type == description ||
                    compact && allele_type.as_str() == &description.replace(" ", "_") {
                    return name;
                } else {
                    return flex_fmt!("{}({})", name, description);
                }
            }

        if name.contains(description.as_str()) {
            if allele_type.contains("amino_acid") || allele_type.contains("amino acid") {
                return name + "(aa)";
            } else if allele_type.contains("nucleotide") {
                return name + "(nt)";
            } else {
                return name;
            }
        }

        description = description_with_residue_type(self);

        if allele_type.starts_with("partial") &&
            ENDS_WITH_INTEGER_RE.is_match(description.as_str()) {
                if allele_type == "partial_amino_acid_deletion" {
                    description = description + " Δaa";
                } else if allele_type == "partial_nucleotide_deletion" {
                    description = description + " Δnt";
                }
            }

        flex_fmt!("{}({})", name, description)
    }

    pub fn short_display_name(&self) -> FlexStr {
        self.display_name_helper(true)
    }

    pub fn display_name(&self) -> FlexStr {
        self.display_name_helper(false)
    }
}

fn allele_encoded_name_and_type(allele_name: &Option<FlexStr>, allele_type: &str,
                                allele_description: &Option<FlexStr>) -> FlexStr {
    let name = allele_name.clone().unwrap_or_else(|| flex_str!("unnamed"));
    let allele_type = allele_type.to_owned();
    let description =
        allele_description.clone().unwrap_or_else(|| allele_type.to_shared_str());

    if allele_type == "deletion" && name.ends_with("delta") ||
        allele_type.starts_with("wild_type") && name.ends_with('+') {
            let normalised_description = description.replace("[\\s_]+", "");
            let normalised_allele_type = allele_type.replace("[\\s_]+", "");
            if normalised_description != normalised_allele_type {
                return flex_fmt!("{}({})", name, description);
            } else {
                return name;
            }
        }

    let display_name =
        if allele_type == "deletion" {
            name + "-" + description.as_str()
        } else {
            name + "-" + description.as_str() + "-" + &allele_type
        };
    display_name
}

impl From<&AlleleDetails> for AlleleShort {
    fn from(details: &AlleleDetails) -> AlleleShort {
        AlleleShort {
           uniquename: details.uniquename.clone(),
           name: details.name.clone(),
           allele_type: details.allele_type.clone(),
           description: details.description.clone(),
           gene_uniquename: details.gene.uniquename.clone(),
           synonyms: details.synonyms.clone(),
           comments: details.comments.clone(),
       }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct CommentAndReference {
    pub comment: FlexStr,
    pub reference: Option<FlexStr>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct AlleleDetails {
    pub uniquename: FlexStr,
    // don't serialise since the web code doesn't use this field:
    #[serde(skip, default="empty_string")]
    pub encoded_name_and_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub allele_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub description: Option<FlexStr>,
    pub gene: GeneShort,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<SynonymDetails>,
    // genotypes containing this allele:
    pub genotypes: Vec<GenotypeShort>,
    // the phenotypes of those genotypes
    pub phenotypes: Vec<TermShort>,

    #[serde(default, skip_serializing_if = "is_default")]
    pub is_obsolete: bool,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub comments: Vec<CommentAndReference>,

    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub references_by_uniquename: ReferenceShortOptionMap,
    pub genes_by_uniquename: GeneShortOptionMap,
}

impl AlleleDetails {
    pub fn new(uniquename: &str,
               name: &Option<FlexStr>,
               allele_type: &str,
               description: &Option<FlexStr>,
               comments: &[CommentAndReference],
               is_obsolete: bool,
               gene: GeneShort) -> AlleleDetails {
        let encoded_name_and_type =
            allele_encoded_name_and_type(name, allele_type, description);
        AlleleDetails {
            uniquename: uniquename.to_shared_str(),
            encoded_name_and_type,
            name: name.clone(),
            allele_type: allele_type.to_shared_str(),
            description: description.clone(),
            gene,
            synonyms: vec![],
            genotypes: vec![],
            phenotypes: vec![],
            is_obsolete,
            comments: comments.to_owned(),
            alleles_by_uniquename: HashMap::new(),
            references_by_uniquename: HashMap::new(),
            genes_by_uniquename: HashMap::new(),
        }
    }
}


pub type RelName = FlexStr;

#[derive(Serialize, Deserialize, Clone, Debug, Hash, PartialEq, Eq)]
pub struct GeneExProps {
    #[serde(skip_serializing_if="Option::is_none")]
    pub copies_per_cell: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub avg_copies_per_cell: Option<FlexStr>,
    pub scale: FlexStr,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GeneExMeasurement {
    pub reference_uniquename: FlexStr,
    pub level_type_termid: FlexStr,
    pub during_termid: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub copies_per_cell: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub avg_copies_per_cell: Option<FlexStr>,
    pub scale: FlexStr,
}

pub type OntName = FlexStr;
pub type OntAnnotationMap = HashMap<OntName, Vec<OntTermAnnotations>>;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermAndRelation {
    pub termid: TermId,
    pub term_name: TermName,
    pub relation_name: RelName,
}

#[serde_as]
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermDetails {
    pub name: TermName,
    pub cv_name: CvName,
    pub annotation_feature_type: FlexStr,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parent_ids: HashSet<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parent_details: HashSet<InterestingParent>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub in_subsets: HashSet<FlexStr>,
    pub termid: TermId,
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

    #[serde(skip_serializing_if="is_true", default)]
    pub is_obsolete: bool,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub single_locus_genotype_uniquenames: HashSet<FlexStr>,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub cv_annotations: OntAnnotationMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub genes_by_uniquename: GeneShortOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub genotypes_by_uniquename: HashMap<GenotypeUniquename, GenotypeShort>,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub references_by_uniquename: ReferenceShortOptionMap,
    // transcripts referenced by this term
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub transcripts_by_uniquename: TranscriptDetailsOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub terms_by_termid: TermShortOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub annotation_details: IdOntAnnotationDetailMap,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    #[serde_as(as = "Vec<(_, _)>")]
    pub double_mutant_genetic_interactions: GeneticInteractionMap,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    #[serde_as(as = "Vec<(_, _)>")]
    pub single_allele_genetic_interactions: GeneticInteractionMap,

    pub gene_count: usize,
    pub genotype_count: usize,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub xrefs: HashMap<FlexStr, TermXref>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pombase_gene_id: Option<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub gocams: HashSet<GoCamIdAndTitle>,
}

impl Container for TermDetails {
    fn container_type(&self) -> ContainerType {
        ContainerType::Term
    }
}

impl AnnotationContainer for TermDetails {
    fn cv_annotations(&self) -> &OntAnnotationMap {
        &self.cv_annotations
    }
    fn cv_annotations_mut(&mut self) -> &mut OntAnnotationMap {
        &mut self.cv_annotations
    }
    fn annotation_details(&self) -> &IdOntAnnotationDetailMap {
        &self.annotation_details
    }
    fn terms_by_termid(&self) -> &TermShortOptionMap {
        &self.terms_by_termid
    }
    fn genes_by_uniquename(&self) -> &GeneShortOptionMap {
        &self.genes_by_uniquename
    }
    fn genotypes_by_uniquename(&self) -> &HashMap<GenotypeUniquename, GenotypeShort> {
       &self.genotypes_by_uniquename
    }
}

// pub type GeneticInteractionKey = (GeneUniquename, Evidence, GeneUniquename);

#[derive(Serialize, Deserialize, Clone, Debug, Eq, PartialEq)]
pub struct GeneticInteractionKey {
    pub gene_a_uniquename: GeneUniquename,
    pub gene_b_uniquename: GeneUniquename,
    pub interaction_type: Evidence,
}
impl Ord for GeneticInteractionKey {
    fn cmp(&self, other: &Self) -> Ordering {
        let order = self.interaction_type.cmp(&other.interaction_type);
        if order != Ordering::Equal {
            return order;
        }

        (&self.gene_a_uniquename, &self.interaction_type, &self.gene_b_uniquename)
            .cmp(&(&other.gene_a_uniquename, &other.interaction_type, &other.gene_b_uniquename))
    }
}
impl PartialOrd for GeneticInteractionKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Hash for GeneticInteractionKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.gene_a_uniquename.hash(state);
        self.interaction_type.hash(state);
        self.gene_b_uniquename.hash(state);
    }
}



#[derive(Serialize, Deserialize, Clone, Debug, Eq, PartialEq)]
pub struct GeneticInteractionDetail {
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_a_uniquename: Option<GenotypeDisplayUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_b_uniquename: Option<GenotypeDisplayUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub double_mutant_phenotype_termid: Option<TermId>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub double_mutant_extension: Vec<ExtPart>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub double_mutant_genotype_display_uniquename: Option<GenotypeDisplayUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub rescued_phenotype_termid: Option<TermId>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub rescued_phenotype_extension: Vec<ExtPart>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference_uniquename: Option<ReferenceUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub throughput: Option<Throughput>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub interaction_note: Option<FlexStr>,
}
impl Ord for GeneticInteractionDetail {
    fn cmp(&self, other: &Self) -> Ordering {
        let order = self.reference_uniquename.cmp(&other.reference_uniquename);
        if order != Ordering::Equal {
            return order;
        }

        (&self.genotype_a_uniquename, &self.genotype_b_uniquename)
            .cmp(&(&other.genotype_a_uniquename, &other.genotype_b_uniquename))
    }
}
impl PartialOrd for GeneticInteractionDetail {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct InteractionAnnotation {
    pub gene_uniquename: GeneUniquename,
    pub interactor_uniquename: GeneUniquename,
    pub evidence: Evidence,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference_uniquename: Option<ReferenceUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub throughput: Option<Throughput>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub interaction_note: Option<FlexStr>,
}
impl PartialEq for InteractionAnnotation {
    fn eq(&self, other: &Self) -> bool {
        (&self.gene_uniquename, &self.interactor_uniquename,
         &self.evidence) ==
            (&other.gene_uniquename, &other.interactor_uniquename,
             &other.evidence)
    }
}
impl Eq for InteractionAnnotation { }
impl Ord for InteractionAnnotation {
    fn cmp(&self, other: &Self) -> Ordering {
        let order = self.evidence.cmp(&other.evidence);
        if order != Ordering::Equal {
            return order;
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
    #[serde(skip_serializing_if="Option::is_none")]
    pub qualifier: Option<FlexStr>,
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

type DataSourceName = FlexStr;
type DataSourceVersion = FlexStr;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Metadata {
    pub db_creation_datetime: FlexStr,
    pub date_version: FlexStr,
    pub export_prog_name: FlexStr,
    pub export_prog_version: FlexStr,
    pub data_source_versions: HashMap<DataSourceName, DataSourceVersion>,
    pub gene_count: usize,
    pub term_count: usize,
    pub cv_versions: HashMap<FlexStr, FlexStr>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIAlleleDetails {
    pub gene: GeneUniquename,
    pub allele_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub expression: Option<FlexStr>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIGenotypeAnnotation {
    pub is_multi: bool,
    pub ploidiness: Ploidiness,
    pub conditions: HashSet<TermAndName>,
    pub alleles: Vec<APIAlleleDetails>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIGeneSummary {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub exact_synonyms: Vec<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub dbxrefs: HashSet<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub pdb_ids: HashSet<PdbId>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub gocam_ids: HashSet<GoCamId>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcripts: Vec<TranscriptDetails>,
    pub tm_domain_count: usize,
    pub coiled_coil_count: usize,
    pub disordered_regions_count: usize,
    pub low_complexity_regions_count: usize,
    pub coding_exon_count: usize,
    pub five_prime_exon_count: usize,
    pub three_prime_exon_count: usize,
    pub transcript_count: usize,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub ortholog_taxonids: HashSet<u32>,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum GeneQueryTermData {
    Term(TermAndName),
    Other,
}

pub type GeneQueryAttrName = FlexStr;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneQueryData {
    pub gene_uniquename: GeneUniquename,
    pub deletion_viability: DeletionViability,
#[serde(skip_serializing_if="Option::is_none")]
    pub go_component: Option<GeneQueryTermData>,
#[serde(skip_serializing_if="Option::is_none")]
    pub go_process_superslim: Option<GeneQueryTermData>,
#[serde(skip_serializing_if="Option::is_none")]
    pub go_function: Option<GeneQueryTermData>,
#[serde(skip_serializing_if="Option::is_none")]
    pub characterisation_status: Option<FlexStr>,
#[serde(skip_serializing_if="Option::is_none")]
    pub taxonomic_distribution: Option<FlexStr>,
#[serde(skip_serializing_if="Option::is_none")]
    pub tmm: Option<PresentAbsent>,
#[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub ortholog_taxonids: HashSet<u32>,
#[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub physical_interactors: HashSet<GeneUniquename>,
    // None for RNA genes:
#[serde(skip_serializing_if="Option::is_none")]
    pub molecular_weight: Option<f32>,
#[serde(skip_serializing_if="Option::is_none")]
    pub protein_length: Option<usize>,
#[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub reference_uniquenames: HashSet<ReferenceUniquename>,
    // bin is None for RNA genes:
#[serde(skip_serializing_if="Option::is_none")]
    pub protein_length_bin: Option<GeneQueryAttrName>,
#[serde(skip_serializing_if="Option::is_none")]
    pub spliced_rna_length: Option<usize>,
#[serde(skip_serializing_if="Option::is_none")]
    pub unspliced_rna_length: Option<usize>,
#[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub subset_termids: HashSet<TermId>,
#[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub pdb_ids: HashSet<PdbId>,
#[serde(skip_serializing_if="Option::is_none")]
    pub rnacentral_urs_identifier: Option<FlexStr>,
#[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub gocam_ids: HashSet<GoCamId>,
}

#[derive(Serialize, Deserialize, Copy, Clone, Debug, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum InteractionType {
    Physical,
    Genetic
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
pub struct APIInteractor {
    pub interaction_type: InteractionType,
    pub interactor_uniquename: GeneUniquename,
}

pub type GeneExDataSetName = FlexStr;

pub type GeneExDataSetMeasurements =
    HashMap<GeneUniquename, HashMap<GeneExDataSetName, GeneExMeasurement>>;

#[derive(PartialEq)]
pub enum ProteinViewType {
    Full,
    Widget,
}

impl TryFrom<&str> for ProteinViewType {
    type Error = String;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "full" => Ok(ProteinViewType::Full),
            "widget" => Ok(ProteinViewType::Widget),
            _ => Err(format!("unknown protein view type: {}", value)),
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct TermNameAndId {
    pub name: TermName,
    pub id: TermId,
}

impl Ord for TermNameAndId {
    fn cmp(&self, other: &TermNameAndId) -> Ordering {
        self.name.cmp(&other.name)
    }
}
impl PartialOrd for TermNameAndId {
    fn partial_cmp(&self, other: &TermNameAndId) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub type ProteinViewFeaturePos = (FlexStr, usize, usize);

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinViewFeature {
    pub id: FlexStr,

    #[serde(skip_serializing_if="Option::is_none")]
    pub display_name: Option<FlexStr>,

    #[serde(skip_serializing_if="BTreeSet::is_empty", default)]
    pub annotated_terms: BTreeSet<TermNameAndId>,

    #[serde(skip_serializing_if="Option::is_none")]
    // used for grouping modifications by type so that the different types
    // can be coloured
    pub feature_group: Option<ProteinFeatureViewModGroup>,

    #[serde(skip_serializing_if="BTreeSet::is_empty", default)]
    pub display_extension: BTreeSet<FlexStr>,

    // start, end pairs:
    pub positions: Vec<ProteinViewFeaturePos>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinViewTrack {
    pub name: FlexStr,
    pub display_type: FlexStr,
    pub features: Vec<ProteinViewFeature>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinViewData {
    pub sequence: FlexStr,
    pub tracks: Vec<ProteinViewTrack>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinComplexDetails {
    pub complex_uniquename: FlexStr,
    pub complex_name: Option<FlexStr>,
    pub genes: HashSet<GeneUniquename>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinComplexGene {
    pub gene_short: GeneShort,
    pub annotation_details: HashSet<(Option<ReferenceUniquename>,
                                 Option<AssignedBy>, Evidence)>
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinComplexTerm {
    pub term_name: TermName,
    pub complex_genes: BTreeMap<GeneUniquename, ProteinComplexGene>
}

pub type ProteinComplexData =
    HashMap<TermId, ProteinComplexTerm>;

pub type GoCamId = FlexStr;
pub type GoCamTitle = FlexStr;

#[derive(Serialize, Deserialize, Clone, Debug, Eq, PartialEq, Hash)]
pub struct GoCamIdAndTitle {
    pub gocam_id: GoCamId,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<GoCamTitle>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamDetails {
    pub gocam_id: GoCamId,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<FlexStr>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub genes: HashSet<GeneUniquename>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub terms: HashSet<TermAndName>,
}

impl GoCamDetails {
    pub fn new(gocam_id: &GoCamId, gocam_title: Option<GoCamTitle>) -> GoCamDetails {
        GoCamDetails {
            gocam_id: gocam_id.to_owned(),
            title: gocam_title,
            genes: HashSet::new(),
            terms: HashSet::new(),
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIMaps {
    pub termid_genes: HashMap<TermId, HashSet<GeneUniquename>>,
    pub gene_summaries: HashMap<GeneUniquename, APIGeneSummary>,
    pub gene_query_data_map: HashMap<GeneUniquename, GeneQueryData>,
    pub transcripts: UniquenameTranscriptMap,
    pub gene_name_gene_map: HashMap<FlexStr, GeneUniquename>,
    pub interactors_of_genes: HashMap<GeneUniquename, Vec<APIInteractor>>,
    pub downstream_genes: HashMap<FlexStr, HashMap<TermId, HashSet<GeneUniquename>>>,
    pub other_features: UniquenameFeatureShortMap,
    pub seq_feature_page_features: Vec<FeatureShort>,
    pub chromosomes: ChrNameDetailsMap,
    pub term_subsets: IdTermSubsetMap,
    pub gene_subsets: IdGeneSubsetMap,
    pub children_by_termid: HashMap<TermId, HashSet<TermId>>,
    pub gene_expression_measurements: GeneExDataSetMeasurements,
    pub secondary_identifiers_map: HashMap<TermId, TermId>,
    pub protein_view_data: HashMap<GeneUniquename, ProteinViewData>,
    pub gocam_data_by_gene: HashMap<GeneUniquename, HashSet<GoCamId>>,
    pub gocam_data_by_gocam_id: HashMap<GoCamId, GoCamDetails>,

    // based on annotations to "protein-containing complex" (GO:0032991):
    pub protein_complex_data: ProteinComplexData,

    // data from Complex Portal:
    pub protein_complexes: ProteinComplexMap,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrGeneSummary {
    pub id: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,
    pub taxonid: OrganismTaxonId,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<FlexStr>,
    pub feature_type: FlexStr,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrAlleleSummary {
    pub id: AlleleUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<FlexStr>,
    pub allele_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub description: Option<FlexStr>,
    pub gene_uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_name: Option<FlexStr>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<FlexStr>,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub highlighting: SolrMatchHighlight,
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrTermSummary {
    pub id: TermId,
    pub name: TermName,
    pub cv_name: CvName,
    #[serde(skip_serializing_if="Option::is_none")]
    pub definition: Option<TermDef>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub close_synonyms: Vec<FlexStr>,   // exact and narrow
    // a uniquified list of the words in all close synonyms
    pub close_synonym_words: FlexStr,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub distant_synonyms: Vec<FlexStr>, // broad and related
    // a uniquified list of the words in all distant synonyms
    pub distant_synonym_words: FlexStr,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parent_ids: HashSet<FlexStr>,

    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub definition_xrefs: HashSet<FlexStr>,

    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub secondary_identifiers: HashSet<TermId>,

    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub gocam_ids: HashSet<GoCamId>,

    pub annotation_count: usize,
    pub gene_count: usize,
    pub genotype_count: usize,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub highlighting: SolrMatchHighlight,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrReferenceSummary {
    pub id: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_abstract: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_publication_date: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_entrez_date: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<u32>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub approved_date: Option<FlexStr>,
    pub gene_count: usize,
    pub genotype_count: usize,
    pub annotation_count: usize,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_annotation_status: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_name: Option<FlexStr>,
    pub canto_curator_role: FlexStr,

    #[serde(skip_serializing_if="Option::is_none")]
    pub file_curator_name: Option<FlexStr>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub file_curator_role: Option<FlexStr>,

    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub highlighting: SolrMatchHighlight,
}

impl SolrReferenceSummary {
    pub fn from_reference_details(reference_details: &ReferenceDetails) -> SolrReferenceSummary {
        let pub_year_as_int: Option<u32> =
            if let Some(ref pub_year) = reference_details.publication_year {
                pub_year.parse().ok()
            } else {
                None
            };

        let annotation_count = reference_details.annotation_details.len() +
            reference_details.genetic_interactions.len() +
            reference_details.physical_interactions.len() +
            reference_details.ortholog_annotations.len() +
            reference_details.paralog_annotations.len();

        SolrReferenceSummary {
            id: reference_details.uniquename.clone(),
            title: reference_details.title.clone(),
            pubmed_abstract: reference_details.pubmed_abstract.clone(),
            citation: reference_details.citation.clone(),
            publication_year: pub_year_as_int,
            pubmed_publication_date: reference_details.pubmed_publication_date.clone(),
            pubmed_entrez_date: reference_details.pubmed_entrez_date.clone(),
            authors: reference_details.authors.clone(),
            authors_abbrev: reference_details.authors_abbrev.clone(),
            approved_date: reference_details.approved_date.clone(),
            gene_count: reference_details.genes_by_uniquename.keys().len(),
            genotype_count: reference_details.genotypes_by_uniquename.keys().len(),
            annotation_count,
            canto_annotation_status: reference_details.canto_annotation_status.clone(),
            canto_curator_name: reference_details.canto_curator_name.clone(),
            canto_curator_role: reference_details.canto_curator_role.clone(),

            file_curator_name: reference_details.file_curator_name.clone(),
            file_curator_role: reference_details.file_curator_role.clone(),

            highlighting: HashMap::new(),
        }
    }
}



#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrData {
    pub term_summaries: Vec<SolrTermSummary>,
    pub gene_summaries: Vec<SolrGeneSummary>,
    pub allele_summaries: Vec<SolrAlleleSummary>,
    pub reference_summaries: Vec<SolrReferenceSummary>,
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct InterMineGeneDetails {
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<SynonymDetails>,

    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,

    pub systematic_id: GeneUniquename,

    pub feature_type: FlexStr,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcripts: Vec<TranscriptDetails>,


    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<FlexStr>,
    pub taxonid: OrganismTaxonId,

    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub references: Vec<ReferenceUniquename>,
}

impl InterMineGeneDetails {
    pub fn from_gene_details(transcripts_by_uniquename: &UniquenameTranscriptMap,
                             gene_details: &GeneDetails) -> InterMineGeneDetails {
        let transcripts = gene_details.transcripts.iter()
            .map(|transcript_uniquename|
                 transcripts_by_uniquename.get(transcript_uniquename)
                 .unwrap_or_else(|| panic!("internal error, failed to find transcript: {}",
                                  transcript_uniquename))
                 .to_owned())
            .collect::<Vec<_>>();

        let references =
            gene_details.references_by_uniquename.keys()
            .map(|uniquename_ref| uniquename_ref.to_owned()).collect();

        InterMineGeneDetails {
            name: gene_details.name.clone(),
            synonyms: gene_details.synonyms.clone(),
            product: gene_details.product.clone(),
            systematic_id: gene_details.uniquename.clone(),
            feature_type: gene_details.feature_type.clone(),
            location: gene_details.location.clone(),
            transcripts,
            uniprot_identifier: gene_details.uniprot_identifier.clone(),
            taxonid: gene_details.taxonid,
            references,
        }
    }
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
    pub name: FlexStr,
    pub gene_count: usize,

    // used for displaying the FYPO slim
    pub single_locus_gene_count: usize,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermSubsetDetails {
    pub name: FlexStr,
    pub total_gene_count: usize, // total unique genes in all subsets

    // total in single locus genotypes in all subsets, used by the FYPO slim
    pub total_single_locus_gene_count: usize,

    pub elements: HashMap<TermId, TermSubsetElement>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneSubsetDetails {
    pub name: FlexStr,
    pub display_name: FlexStr,
    pub elements: HashSet<GeneUniquename>,
}

pub type IdTermSubsetMap = HashMap<FlexStr, TermSubsetDetails>;
pub type IdGeneSubsetMap = HashMap<FlexStr, GeneSubsetDetails>;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct StatCountsByTaxon {
    pub genes: usize,
    pub annotations: usize,
    #[serde(skip_serializing_if="Option::is_none")]
    pub annotation_type_counts_by_year: Option<StatsIntegerTable>,
}

// could be a year, a year-month "2021-02" or a date "2021-02-14"
pub type DateString = String;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Stats {
    pub by_taxon: HashMap<OrganismTaxonId, StatCountsByTaxon>,
    pub community_pubs_count: usize,
    pub non_community_pubs_count: usize,
}


pub type StatsIntegerTableRow = (DateString, Vec<usize>);

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct StatsIntegerTable {
    pub header: Vec<String>,
    pub data: Vec<StatsIntegerTableRow>,
}

pub type StatsFloatTableRow = (DateString, f32);

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct StatsFloatTable {
    pub header: Vec<String>,
    pub data: Vec<StatsFloatTableRow>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct DetailedStats {
    pub curated_by_month: StatsIntegerTable,
    pub curated_by_year: StatsIntegerTable,
    pub cumulative_curated_by_month: StatsIntegerTable,
    pub cumulative_curated_by_year: StatsIntegerTable,
    pub ltp_genes_per_pub_per_year_range: StatsFloatTable,
    pub ltp_annotations_per_pub_per_year_range: StatsFloatTable,
    pub htp_annotations_per_pub_per_year_range: StatsFloatTable,
    pub community_response_rates: Vec<CommunityResponseRate>,
    pub annotation_type_counts_by_year: StatsIntegerTable,
    pub cumulative_annotation_type_counts_by_year: StatsIntegerTable,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ApicuronReport {
    pub activity_term: FlexStr,
    pub timestamp: FlexStr,
    pub curator_orcid: FlexStr,
    pub entity_uri: FlexStr,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ApicuronData {
    pub resource_id: FlexStr,
    pub reports: Vec<ApicuronReport>,
}
