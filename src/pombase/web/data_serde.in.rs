use std::rc::Rc;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::collections::HashSet;


#[derive(Serialize, Clone)]
pub enum ExtRange {
#[serde(rename = "gene")]
    Gene(GeneShort),
#[serde(rename = "term")]
    Term(TermShort),
#[serde(rename = "misc")]
    Misc(String),
}

#[derive(Serialize, Clone)]
pub struct ExtPart {
    pub rel_type_name: String,
    pub ext_range: ExtRange,
}

#[derive(Serialize, Clone)]
pub struct GeneShort {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,
    pub synonyms: Vec<SynonymDetails>,
}

pub struct TermIdPair {
    pub firstid: TermId,
    pub secondid: TermId,
}

impl TermIdPair {
    pub fn new(firstid: &TermId, secondid: &TermId) -> TermIdPair {
        TermIdPair {
            firstid: firstid.clone(),
            secondid: secondid.clone(),
        }
    }
}

impl PartialEq for TermIdPair {
    fn eq(&self, other: &TermIdPair) -> bool { 
        self.firstid == other.firstid &&
            self.secondid == other.secondid
    }
}
impl Eq for TermIdPair { }
impl Hash for TermIdPair {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.firstid.hash(state);
        self.secondid.hash(state);
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

#[derive(Serialize, Clone)]
pub struct TranscriptShort {
    pub uniquename: TranscriptUniquename,
    //                pub exons: Vec<ExonShort>,
    //                pub utrs: Vec<UTRShort>,
}

#[derive(Serialize, Clone)]
pub struct TermShort {
    pub name: TermName,
    pub cv_name: String,
    pub termid: TermId,
    pub is_obsolete: bool,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_count: Option<usize>,
}

#[derive(Serialize, Clone)]
pub struct ReferenceShort {
    pub uniquename: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub title: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub citation: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub authors_abbrev: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<String>,
}

#[derive(Serialize, Clone)]
pub struct ReferenceDetails {
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
    pub pubmed_publication_date: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub publication_year: Option<String>,
    pub annotations: OntAnnotationMap,
    pub interaction_annotations: TypeInteractionAnnotationMap,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
}

#[derive(Serialize, Clone)]
pub struct OntAnnotation {
    pub id: i32,
    #[serde(skip_serializing_if="Option::is_none")]
    pub term: Option<TermShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene: Option<GeneShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference: Option<ReferenceShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    pub extension: Vec<ExtPart>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub with: Option<With>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub residue: Option<Residue>,
    pub qualifiers: Vec<Qualifier>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_ex_props: Option<GeneExProps>,
    // only for genotype/phenotype annotation:
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype: Option<GenotypeShort>,
    pub conditions: Vec<TermShort>,
    pub is_not: bool,
}

#[derive(Serialize, Clone)]
pub struct SynonymDetails {
    pub name: String,
    pub synonym_type: String
}

#[derive(Serialize, Clone)]
pub enum Strand {
    #[serde(rename="forward")]
    Forward = 1,
    #[serde(rename="reverse")]
    Reverse = -1,
}

#[derive(Serialize, Clone)]
pub struct ChromosomeLocation {
    pub chromosome_name: String,
    pub start_pos: u32,
    pub end_pos: u32,
    pub strand: Strand,
}

#[derive(Serialize, Clone)]
pub struct OrganismShort {
    pub genus: String,
    pub species: String,
}

#[derive(Serialize, Clone)]
pub struct GeneDetails {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    pub organism: OrganismShort,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<String>,
    pub synonyms: Vec<SynonymDetails>,
    pub feature_type: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub characterisation_status: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<ChromosomeLocation>,
    pub transcripts: Vec<TranscriptShort>,
    pub annotations: OntAnnotationMap,
    pub interaction_annotations: TypeInteractionAnnotationMap,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
}

#[derive(Serialize, Clone)]
pub struct TranscriptDetails {
    pub uniquename: TranscriptUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
//    pub annotations: TypeFeatureAnnotationMap,
}

#[derive(Serialize, Clone)]
pub struct GenotypeShort {
    pub uniquename: GenotypeUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub background: Option<String>,
    pub expressed_alleles: Vec<ExpressedAllele>,
//    pub annotations: TypeFeatureAnnotationMap,
}

#[derive(Serialize, Clone)]
pub struct ExpressedAllele {
    #[serde(skip_serializing_if="Option::is_none")]
    pub expression: Option<String>,
    pub allele: AlleleShort,
}

#[derive(Serialize, Clone)]
pub struct AlleleShort {
    pub uniquename: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<String>,
    pub allele_type: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub description: Option<String>,
    pub gene: GeneShort,
}

pub type RelName = String;

#[derive(Serialize, Clone)]
pub struct GeneExProps {
    #[serde(skip_serializing_if="Option::is_none")]
    pub copies_per_cell: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub avg_copies_per_cell: Option<String>,
}

pub type OntName = String;
pub type OntAnnotationMap = HashMap<OntName, Vec<Rc<OntAnnotation>>>;

#[derive(Serialize, Clone)]
pub struct RelOntAnnotation {
    pub rel_names: HashSet<RelName>,
    pub annotation: Rc<OntAnnotation>,
}

#[derive(Serialize, Clone)]
pub struct TermDetails {
    pub name: TermName,
    pub cv_name: CvName,
    pub termid: TermId,
    #[serde(skip_serializing_if="Option::is_none")]
    pub definition: Option<TermDef>,
    pub is_obsolete: bool,
    pub genes: Vec<GeneShort>,
    pub annotations: Vec<RelOntAnnotation>,
}

#[derive(Serialize, Clone)]
pub struct InteractionAnnotation {
    pub gene: GeneShort,
    pub interactor: GeneShort,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference: Option<ReferenceShort>,
}

#[derive(Serialize, Clone)]
pub struct OrthologAnnotation {
    pub gene: GeneShort,
    pub ortholog_organism: OrganismShort,
    pub ortholog: GeneShort,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference: Option<ReferenceShort>,
}

#[derive(Serialize, Clone)]
pub struct ParalogAnnotation {
    pub gene: GeneShort,
    pub paralog: GeneShort,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference: Option<ReferenceShort>,
}

#[derive(Serialize, Clone)]
pub struct Metadata {
    pub db_creation_datetime: String,
    pub export_prog_name: String,
    pub export_prog_version: String,
    pub gene_count: usize,
    pub term_count: usize,
}

#[derive(Serialize, Clone)]
pub struct WebData {
    pub genes: IdGeneMap,
    pub gene_summaries: IdGeneShortMap,
    pub terms: IdTermDetailsMap,
    pub used_terms: IdTermDetailsMap,
    pub metadata: Metadata,
    pub references: IdReferenceMap,
}
