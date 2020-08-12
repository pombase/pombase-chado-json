use std::collections::{HashMap, HashSet, BTreeMap};

use std::fmt::Display;
use std::fmt;

use pombase_rc_string::RcString;

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
pub type UniquenameFeatureShortMap = HashMap<RcString, FeatureShort>;
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

use std::rc::Rc;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

use crate::web::config::*;
use crate::types::*;
use crate::interpro::InterProMatch;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum Throughput {
#[serde(rename = "high")]
    HighThroughput,
#[serde(rename = "low")]
    LowThroughput,
#[serde(rename = "non-experimental")]
    NonExperimental,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, PartialOrd, Hash)]
pub enum ExtRange {
#[serde(rename = "gene_uniquename")]
    Gene(GeneUniquename),
#[serde(rename = "promoter_gene_uniquename")]
    Promoter(RcString),
#[serde(rename = "summary_gene_uniquenames")]
    // the inner Vec length will be > 1 for cases like "binds abc1 and def2, cdc2"
    SummaryGenes(Vec<Vec<RcString>>),
#[serde(rename = "termid")]
    Term(TermId),
#[serde(rename = "summary_termids")]
    // See: merge_ext_part_ranges()
    SummaryTerms(Vec<TermId>),
#[serde(rename = "misc")]
    Misc(RcString),
#[serde(rename = "domain")]
    Domain(RcString),
#[serde(rename = "gene_product")]
    GeneProduct(TermId),  // eg.  PR:000027705
#[serde(rename = "summary_residues")]
    SummaryModifiedResidues(Vec<Residue>),
}

impl ExtRange {
    pub fn is_gene(&self) -> bool {
        match *self {
            ExtRange::Gene(_) => true,
            _ => false,
        }
    }
}

impl fmt::Display for ExtRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExtRange::Gene(ref gene_uniquename) | ExtRange::Promoter(ref gene_uniquename) =>
                write!(f, "{}", gene_uniquename),
            ExtRange::SummaryGenes(_) => panic!("can't handle SummaryGenes\n"),
            ExtRange::Term(ref termid) => write!(f, "{}", termid),
            ExtRange::SummaryModifiedResidues(ref residue) =>
                write!(f, "{}", &residue.join(",")),
            ExtRange::SummaryTerms(_) => panic!("can't handle SummaryGenes\n"),
            ExtRange::Misc(ref misc) => write!(f, "{}", misc),
            ExtRange::Domain(ref domain) => write!(f, "{}", domain),
            ExtRange::GeneProduct(ref gene_product) => write!(f, "{}", gene_product),
        }
    }
}

// A single part of an extension.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ExtPart {
    pub rel_type_name: RcString,
    pub rel_type_display_name: RcString,
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

// minimal information about a gene used in other objects
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneShort {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,
}

impl GeneShort {
    pub fn from_gene_details(gene_details: &GeneDetails) -> Self {
        GeneShort {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
        }
    }

    pub fn display_name(&self) -> String {
        if let Some(ref name) = self.name {
            format!("{} ({})", name, self.uniquename)
        } else {
            String::from(&self.uniquename)
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
    pub identifier: RcString,
    pub taxonid: u32,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct IdNameAndOrganism {
    pub identifier: RcString,
    pub name: Option<RcString>,
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
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<RcString>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<RcString>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub orthologs: Vec<IdNameAndOrganism>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    pub feature_type: RcString,
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermXref {
    pub xref_id: RcString,
    #[serde(skip_serializing_if="Option::is_none")]
    pub xref_display_name: Option<RcString>,
}

// minimal information about a terms used in other objects
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermShort {
    pub name: TermName,
    pub cv_name: RcString,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parents: HashSet<RcString>,
    pub termid: TermId,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub secondary_identifiers: HashSet<TermId>,
    pub is_obsolete: bool,
    pub gene_count: usize,
    pub genotype_count: usize,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub xrefs: HashMap<RcString, TermXref>,
}

impl TermShort {
    pub fn from_term_details(term_details: &TermDetails) -> Self {
        TermShort {
            name: term_details.name.clone(),
            cv_name: term_details.cv_name.clone(),
            interesting_parents: term_details.interesting_parents.clone(),
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

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChromosomeDetails {
    pub name: RcString,
    pub residues: RcString,
    pub ena_identifier: RcString,
    pub gene_uniquenames: Vec<RcString>,
    pub taxonid: OrganismTaxonId,
}

impl ChromosomeDetails {
    pub fn make_chromosome_short(&self) -> ChromosomeShort {
        ChromosomeShort {
            name: self.name.clone(),
            length: self.residues.len(),
            ena_identifier: self.ena_identifier.clone(),
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ReferenceShort {
    pub uniquename: RcString,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub approved_date: Option<RcString>,
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
            gene_count: reference_details.genes_by_uniquename.keys().len(),
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
    fn genotypes_by_uniquename(&self) -> Option<&HashMap<GenotypeUniquename, GenotypeShort>>;
}

pub trait OrthologAnnotationContainer: AnnotationContainer {
    fn ortholog_annotations(&self) -> &Vec<OrthologAnnotation>;
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ReferenceDetails {
    pub uniquename: RcString,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none", rename = "abstract")]
    pub pubmed_abstract: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none", rename = "doi")]
    pub pubmed_doi: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_publication_date: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_annotation_status: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_triage_status: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_role: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_name: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_first_approved_date: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_approved_date: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_session_submitted_date: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_added_date: Option<RcString>,

    // This is set to the year part of canto_first_approved_date if it is
    // not None, otherwise set to the year part of canto_approved_date, otherwise
    // canto_session_submitted_date
    #[serde(skip_serializing_if="Option::is_none")]
    pub approved_date: Option<RcString>,
    pub cv_annotations: OntAnnotationMap,
    pub physical_interactions: Vec<InteractionAnnotation>,
    pub genetic_interactions: Vec<InteractionAnnotation>,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
    pub genes_by_uniquename: GeneShortOptionMap,
    pub genotypes_by_uniquename: HashMap<GenotypeUniquename, GenotypeShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub terms_by_termid: TermShortOptionMap,
    pub annotation_details: IdOntAnnotationDetailMap,
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
    fn genotypes_by_uniquename(&self) -> Option<&HashMap<GenotypeUniquename, GenotypeShort>> {
        Some(&self.genotypes_by_uniquename)
    }
}

impl OrthologAnnotationContainer for ReferenceDetails {
    fn ortholog_annotations(&self) -> &Vec<OrthologAnnotation> {
        &self.ortholog_annotations
    }
}

// the GO with/from
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum WithFromValue {
    Gene(GeneShort),
    Term(TermShort),
    Identifier(RcString)
}

impl From<WithFromValue> for RcString {
    fn from(with_from: WithFromValue) -> Self {
        match with_from {
            WithFromValue::Gene(gene) => gene.uniquename,
            WithFromValue::Term(term) => term.termid,
            WithFromValue::Identifier(id) => id,
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct OntAnnotationDetail {
    pub id: i32,
    pub genes: Vec<GeneUniquename>,
    pub reference: Option<ReferenceUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    pub extension: Vec<ExtPart>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub withs: HashSet<WithFromValue>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub froms: HashSet<WithFromValue>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub residue: Option<Residue>,
    pub qualifiers: Vec<Qualifier>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_ex_props: Option<GeneExProps>,
    // only for genotype/phenotype annotation:
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype: Option<GenotypeUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_background: Option<RcString>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub conditions: HashSet<TermId>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub date: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub assigned_by: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub throughput: Option<Throughput>,
}

impl PartialEq for OntAnnotationDetail {
    fn eq(&self, other: &OntAnnotationDetail) -> bool {
        self.id == other.id
    }
}
impl Eq for OntAnnotationDetail { }

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
    pub qualifiers: Vec<Qualifier>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_ex_props: Option<GeneExProps>,
    // only for genotype/phenotype annotation:
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_short: Option<GenotypeShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype_background: Option<RcString>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub conditions: HashSet<TermShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub assigned_by: Option<RcString>,
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
    pub ontology_name: RcString,
    pub ext_rel_display_name: RcString,
    pub gene: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none", default)]
    pub genotype_uniquename: Option<GenotypeUniquename>,
    pub reference_uniquename: Option<ReferenceUniquename>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SynonymDetails {
    pub name: RcString,
    #[serde(rename = "type")]
    pub synonym_type: RcString
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
    pub name: RcString,
    pub length: usize,
    pub ena_identifier: RcString,
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
    pub chromosome_name: RcString,
    pub start_pos: usize,
    pub end_pos: usize,
    pub strand: Strand,
    #[serde(skip_serializing_if="Option::is_none")]
    pub phase: Option<Phase>,
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
pub struct GeneDetails {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<RcString>,
    pub taxonid: u32,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<RcString>,
    pub deletion_viability: DeletionViability,
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<RcString>,
    pub biogrid_interactor_id: Option<u32>,
    pub interpro_matches: Vec<InterProMatch>,
    // non-InterPro domains:
    pub tm_domain_coords: Vec<(usize, usize) >,
    pub orfeome_identifier: Option<RcString>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub name_descriptions: Vec<RcString>,
    pub synonyms: Vec<SynonymDetails>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub dbxrefs: HashSet<RcString>,
    pub feature_type: RcString,
    pub feature_so_termid: RcString,
    pub transcript_so_termid: TermId,
    #[serde(skip_serializing_if="Option::is_none")]
    pub characterisation_status: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub taxonomic_distribution: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    pub gene_neighbourhood: Vec<GeneShort>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcripts: Vec<TranscriptDetails>,
    pub cv_annotations: OntAnnotationMap,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub physical_interactions: Vec<InteractionAnnotation>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub genetic_interactions: Vec<InteractionAnnotation>,
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
    pub terms_by_termid: TermShortOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub annotation_details: IdOntAnnotationDetailMap,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub feature_publications: HashSet<ReferenceUniquename>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    // A Vec of the term IDs of subsets for this gene.  Any useful subset
    // that contains any term for any annotation in the gene is included.
    // "useful" means that the front end might need it, eg. slim term IDs
    pub subset_termids: HashSet<TermId>,
}

impl GeneDetails {
    pub fn spliced_transcript_sequence(&self) -> Option<RcString> {
        if self.transcripts.len() > 1 {
            panic!("no support for multi-transcript genes");
        }

        if let Some(transcript) = self.transcripts.get(0) {
            let mut seq = String::new();

            for part in &transcript.parts {
                if part.feature_type == FeatureType::Exon {
                    seq += &part.residues;
                }
            }

            Some(RcString::from(&seq))
        } else {
            None
        }
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
    fn genotypes_by_uniquename(&self) -> Option<&HashMap<GenotypeUniquename, GenotypeShort>> {
        Some(&self.genotypes_by_uniquename)
    }
}

impl OrthologAnnotationContainer for GeneDetails {
    fn ortholog_annotations(&self) -> &Vec<OrthologAnnotation> {
        &self.ortholog_annotations
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ProteinDetails {
    pub uniquename: TranscriptUniquename,
    pub sequence: RcString,
    pub molecular_weight: f32,
    pub average_residue_weight: f32,
    pub charge_at_ph7: f32,
    pub isoelectric_point: f32,
    pub codon_adaptation_index: f32,
}

pub type Residues = RcString;

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
    SNP,
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
            FeatureType::SNP => "SNP",
        })
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct FeatureShort {
    pub feature_type: FeatureType,
    pub uniquename: RcString,
    pub location: ChromosomeLocation,
    pub residues: Residues,
}


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TranscriptDetails {
    pub uniquename: TranscriptUniquename,
    pub location: ChromosomeLocation,
    pub parts: Vec<FeatureShort>,
    pub transcript_type: RcString,
    pub protein: Option<ProteinDetails>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<ChromosomeLocation>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GenotypeShort {
    pub display_uniquename: GenotypeUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<RcString>,
    pub expressed_alleles: Vec<ExpressedAllele>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GenotypeDetails {
    pub display_uniquename: GenotypeUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<RcString>,
    pub expressed_alleles: Vec<ExpressedAllele>,
    pub cv_annotations: OntAnnotationMap,
    pub references_by_uniquename: ReferenceShortOptionMap,
    pub genes_by_uniquename: GeneShortOptionMap,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub terms_by_termid: TermShortOptionMap,
    pub annotation_details: IdOntAnnotationDetailMap,
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
    fn genotypes_by_uniquename(&self) -> Option<&HashMap<GenotypeUniquename, GenotypeShort>> {
        None
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ExpressedAllele {
    #[serde(skip_serializing_if="Option::is_none")]
    pub expression: Option<RcString>,
    pub allele_uniquename: AlleleUniquename,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct AlleleShort {
    pub uniquename: RcString,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<RcString>,
    pub allele_type: RcString,
    #[serde(skip_serializing_if="Option::is_none")]
    pub description: Option<RcString>,
    pub gene_uniquename: GeneUniquename,
}

pub type RelName = RcString;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneExProps {
    #[serde(skip_serializing_if="Option::is_none")]
    pub copies_per_cell: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub avg_copies_per_cell: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub scale: Option<RcString>,
}

pub type OntName = RcString;
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
    pub annotation_feature_type: RcString,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parents: HashSet<RcString>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub in_subsets: HashSet<RcString>,
    pub termid: TermId,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<SynonymDetails>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub definition: Option<TermDef>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub direct_ancestors: Vec<TermAndRelation>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub secondary_identifiers: HashSet<TermId>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub genes_annotated_with: HashSet<GeneUniquename>,
    pub is_obsolete: bool,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub single_allele_genotype_uniquenames: HashSet<RcString>,
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
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub terms_by_termid: TermShortOptionMap,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub annotation_details: IdOntAnnotationDetailMap,
    pub gene_count: usize,
    pub genotype_count: usize,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub xrefs: HashMap<RcString, TermXref>,
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
    fn genotypes_by_uniquename(&self) -> Option<&HashMap<GenotypeUniquename, GenotypeShort>> {
        Some(&self.genotypes_by_uniquename)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct InteractionAnnotation {
    pub gene_uniquename: GeneUniquename,
    pub interactor_uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference_uniquename: Option<ReferenceUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub throughput: Option<Throughput>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub interaction_note: Option<RcString>,
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
    pub db_creation_datetime: RcString,
    pub export_prog_name: RcString,
    pub export_prog_version: RcString,
    pub gene_count: usize,
    pub term_count: usize,
    pub cv_versions: HashMap<RcString, RcString>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIAlleleDetails {
    pub gene: GeneUniquename,
    pub allele_type: RcString,
    #[serde(skip_serializing_if="Option::is_none")]
    pub expression: Option<RcString>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIGenotypeAnnotation {
    pub is_multi: bool,
    pub conditions: HashSet<TermAndName>,
    pub alleles: Vec<APIAlleleDetails>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIGeneSummary {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub uniprot_identifier: Option<RcString>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub exact_synonyms: Vec<RcString>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub dbxrefs: HashSet<RcString>,
    pub location: Option<ChromosomeLocation>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub transcripts: Vec<TranscriptDetails>,
    pub tm_domain_count: usize,
    pub exon_count: usize,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum GeneQueryTermData {
    Term(TermAndName),
    Other,
}

pub type GeneQueryAttrName = RcString;

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
    pub characterisation_status: Option<RcString>,
#[serde(skip_serializing_if="Option::is_none")]
    pub taxonomic_distribution: Option<RcString>,
#[serde(skip_serializing_if="Option::is_none")]
    pub tmm: Option<PresentAbsent>,
    pub ortholog_taxonids: HashSet<u32>,
    pub physical_interactors: HashSet<GeneUniquename>,
    // None for RNA genes:
    pub molecular_weight: Option<f32>,
    pub protein_length: Option<usize>,
    // bin is None for RNA genes:
    pub protein_length_bin: Option<GeneQueryAttrName>,
#[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub subset_termids: HashSet<TermId>,
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

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct APIMaps {
    pub termid_genes: HashMap<TermId, HashSet<GeneUniquename>>,
    pub termid_genotype_annotation: HashMap<TermId, Vec<APIGenotypeAnnotation>>,
    pub gene_summaries: HashMap<GeneUniquename, APIGeneSummary>,
    pub gene_query_data_map: HashMap<GeneUniquename, GeneQueryData>,
    pub term_summaries: HashSet<TermShort>,
    pub genes: UniquenameGeneMap,
    pub gene_name_gene_map: HashMap<RcString, GeneUniquename>,
    pub genotypes: IdGenotypeMap,
    pub terms: HashMap<TermId, TermDetails>,
    pub interactors_of_genes: HashMap<GeneUniquename, Vec<APIInteractor>>,
    pub references: UniquenameReferenceMap,
    pub other_features: UniquenameFeatureShortMap,
    pub annotation_details: IdOntAnnotationDetailMap,
    pub chromosomes: ChrNameDetailsMap,
    pub term_subsets: IdTermSubsetMap,
    pub gene_subsets: IdGeneSubsetMap,
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
    pub uniprot_identifier: Option<RcString>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub synonyms: Vec<RcString>,
    pub feature_type: RcString,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrTermSummary {
    pub id: TermId,
    pub name: TermName,
    pub cv_name: CvName,
    #[serde(skip_serializing_if="Option::is_none")]
    pub definition: Option<TermDef>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub close_synonyms: Vec<RcString>,   // exact and narrow
    // a uniquified list of the words in all close synonyms
    pub close_synonym_words: RcString,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub distant_synonyms: Vec<RcString>, // broad and related
    // a uniquified list of the words in all distant synonyms
    pub distant_synonym_words: RcString,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parents: HashSet<RcString>,
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub secondary_identifiers: HashSet<TermId>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrReferenceSummary {
    pub id: RcString,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_abstract: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub pubmed_publication_date: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<u32>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub approved_date: Option<RcString>,
    pub gene_count: usize,
    pub genotype_count: usize,
    pub annotation_count: usize,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_annotation_status: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_name: Option<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub canto_curator_role: Option<RcString>,
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
            authors: reference_details.authors.clone(),
            authors_abbrev: reference_details.authors_abbrev.clone(),
            approved_date: reference_details.approved_date.clone(),
            gene_count: reference_details.genes_by_uniquename.keys().len(),
            genotype_count: reference_details.genotypes_by_uniquename.keys().len(),
            annotation_count,
            canto_annotation_status: reference_details.canto_annotation_status.clone(),
            canto_curator_name: reference_details.canto_curator_name.clone(),
            canto_curator_role: reference_details.canto_curator_role.clone(),
        }
    }
}



#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SolrData {
    pub term_summaries: Vec<SolrTermSummary>,
    pub gene_summaries: Vec<SolrGeneSummary>,
    pub reference_summaries: Vec<SolrReferenceSummary>,
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
    pub name: RcString,
    pub gene_count: usize,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct TermSubsetDetails {
    pub name: RcString,
    pub total_gene_count: usize, // total unique genes in all subsets
    pub elements: HashMap<TermId, TermSubsetElement>,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GeneSubsetDetails {
    pub name: RcString,
    pub display_name: RcString,
    pub elements: HashSet<GeneUniquename>,
}

pub type IdTermSubsetMap = HashMap<RcString, TermSubsetDetails>;
pub type IdGeneSubsetMap = HashMap<RcString, GeneSubsetDetails>;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct StatCountsByTaxon {
    pub genes: usize,
    pub annotations: usize,
}

impl StatCountsByTaxon {
    pub fn empty() -> Self {
        StatCountsByTaxon {
            genes: 0,
            annotations: 0,
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Stats {
    pub by_taxon: HashMap<OrganismTaxonId, StatCountsByTaxon>,
    pub community_pubs_count: usize,
    pub non_community_pubs_count: usize,
}
