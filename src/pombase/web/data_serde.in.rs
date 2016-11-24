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
    pub name: Option<GeneName>,
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
    pub gene_count: Option<usize>,
}

#[derive(Serialize, Clone)]
pub struct ReferenceShort {
    pub uniquename: String,
    pub title: Option<String>,
    pub citation: Option<String>,
    pub authors_abbrev: Option<String>,
    pub publication_year: Option<String>,
}

#[derive(Serialize, Clone)]
pub struct ReferenceDetails {
    pub uniquename: String,
    pub title: Option<String>,
    pub citation: Option<String>,
    pub authors: Option<String>,
    pub authors_abbrev: Option<String>,
    pub pubmed_publication_date: Option<String>,
    pub publication_year: Option<String>,
    pub annotations: OntAnnotationMap,
    pub interaction_annotations: TypeInteractionAnnotationMap,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
}

#[derive(Serialize, Clone)]
pub struct OntAnnotation {
    pub id: i32,
    pub term: Option<TermShort>,
    pub gene: Option<GeneShort>,
    pub reference: Option<ReferenceShort>,
    pub evidence: Option<Evidence>,
    pub extension: Vec<ExtPart>,
    pub with: Option<With>,
    pub residue: Option<Residue>,
    pub qualifiers: Vec<Qualifier>,
    // only for genotype/phenotype annotation:
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
    pub name: Option<String>,
    pub organism: OrganismShort,
    pub product: Option<String>,
    pub synonyms: Vec<SynonymDetails>,
    pub feature_type: String,
    pub characterisation_status: Option<String>,
    pub location: Option<ChromosomeLocation>,
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
    pub name: Option<String>,
//    pub annotations: TypeFeatureAnnotationMap,
}

#[derive(Serialize, Clone)]
pub struct GenotypeShort {
    pub uniquename: GenotypeUniquename,
    pub name: Option<String>,
    pub background: Option<String>,
    pub expressed_alleles: Vec<ExpressedAllele>,
//    pub annotations: TypeFeatureAnnotationMap,
}

#[derive(Serialize, Clone)]
pub struct ExpressedAllele {
    pub expression: Option<String>,
    pub allele: AlleleShort,
}

#[derive(Serialize, Clone)]
pub struct AlleleShort {
    pub uniquename: String,
    pub name: Option<String>,
    pub allele_type: String,
    pub description: Option<String>,
    pub gene: GeneShort,
}

pub type RelName = String;

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
    pub definition: Option<TermDef>,
    pub is_obsolete: bool,
    pub genes: Vec<GeneShort>,
    pub annotations: Vec<RelOntAnnotation>,
}

#[derive(Serialize, Clone)]
pub struct InteractionAnnotation {
    pub gene: GeneShort,
    pub interactor: GeneShort,
    pub evidence: Option<Evidence>,
    pub reference: Option<ReferenceShort>,
}

#[derive(Serialize, Clone)]
pub struct OrthologAnnotation {
    pub gene: GeneShort,
    pub ortholog_organism: OrganismShort,
    pub ortholog: GeneShort,
    pub evidence: Option<Evidence>,
    pub reference: Option<ReferenceShort>,
}

#[derive(Serialize, Clone)]
pub struct ParalogAnnotation {
    pub gene: GeneShort,
    pub paralog: GeneShort,
    pub evidence: Option<Evidence>,
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
