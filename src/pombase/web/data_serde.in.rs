use std::rc::Rc;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::collections::HashSet;


#[derive(Serialize, Clone)]
pub enum ExtRange {
#[serde(rename = "gene_uniquename")]
    Gene(GeneUniquename),
#[serde(rename = "termid")]
    Term(TermId),
#[serde(rename = "misc")]
    Misc(String),
}

#[derive(Serialize, Clone)]
pub struct ExtPart {
    pub rel_type_name: String,
    pub rel_type_display_name: String,
    pub ext_range: ExtRange,
}

#[derive(Serialize, Clone)]
pub struct GeneShort {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,
}

#[derive(Serialize, Clone)]
pub struct GeneSummary {
    pub uniquename: GeneUniquename,
    #[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<GeneName>,
    pub organism: OrganismShort,
    #[serde(skip_serializing_if="Option::is_none")]
    pub product: Option<GeneProduct>,
    pub synonyms: Vec<SynonymDetails>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    pub feature_type: String,
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
    #[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub interesting_parents: HashSet<String>,
    pub termid: TermId,
    pub is_obsolete: bool,
    pub gene_count: usize,
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
    pub cv_annotations: OntAnnotationMap,
    pub physical_interactions: Vec<InteractionAnnotation>,
    pub genetic_interactions: Vec<InteractionAnnotation>,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
    pub genes_by_uniquename: HashMap<GeneUniquename, GeneShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub terms_by_termid: HashMap<TermId, TermShort>,
}

#[derive(Serialize, Clone)]
pub struct OntTermAnnotations {
    pub term: TermShort,
    pub annotations: Vec<Rc<OntAnnotationDetail>>,
}

#[derive(Serialize, Clone)]
pub struct OntAnnotationDetail {
    pub id: i32,
    pub gene_uniquename: GeneUniquename,
    pub reference_uniquename: Option<ReferenceUniquename>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub evidence: Option<Evidence>,
    pub extension: Vec<ExtPart>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub with: Option<GeneShort>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub residue: Option<Residue>,
    pub qualifiers: Vec<Qualifier>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub gene_ex_props: Option<GeneExProps>,
    // only for genotype/phenotype annotation:
    #[serde(skip_serializing_if="Option::is_none")]
    pub genotype: Option<GenotypeShort>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub conditions: Vec<TermId>,
    pub is_not: bool,
}

impl PartialEq for OntTermAnnotations {
    fn eq(&self, other: &OntTermAnnotations) -> bool {
        self.term.termid == other.term.termid
    }
}
impl Eq for OntTermAnnotations { }
impl Ord for OntTermAnnotations {
    fn cmp(&self, other: &OntTermAnnotations) -> Ordering {
        return self.term.name.cmp(&other.term.name);
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
    }
}

#[derive(Serialize, Clone, PartialEq, Eq, Hash)]
pub struct TargetOfAnnotation {
    pub ontology_name: String,
    pub ext_rel_display_name: String,
    pub gene_uniquename: GeneUniquename,
    pub reference_uniquename: Option<ReferenceUniquename>,
}

#[derive(Serialize, Clone)]
pub struct SynonymDetails {
    pub name: String,
    #[serde(rename = "type")]
    pub synonym_type: String
}

#[derive(Serialize, Clone, Debug, PartialEq)]
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
    pub name_descriptions: Vec<String>,
    pub synonyms: Vec<SynonymDetails>,
    pub feature_type: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub characterisation_status: Option<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub location: Option<ChromosomeLocation>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub cds_location: Option<ChromosomeLocation>,
    pub gene_neighbourhood: Vec<GeneShort>,
    pub transcripts: Vec<TranscriptShort>,
    pub cv_annotations: OntAnnotationMap,
    pub physical_interactions: Vec<InteractionAnnotation>,
    pub genetic_interactions: Vec<InteractionAnnotation>,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
    pub target_of_annotations: Vec<TargetOfAnnotation>,
    pub references_by_uniquename: HashMap<ReferenceUniquename, ReferenceShort>,
    // genes mentioned in orthologs, paralogs and interactions
    pub genes_by_uniquename: HashMap<GeneUniquename, GeneShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub terms_by_termid: HashMap<TermId, TermShort>,
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
    pub allele_uniquename: AlleleUniquename,
}

#[derive(Serialize, Clone)]
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

#[derive(Serialize, Clone)]
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

#[derive(Serialize, Clone)]
pub struct RelOntAnnotation {
    pub term: TermShort,
    pub rel_names: HashSet<RelName>,
    pub annotations: Vec<Rc<OntAnnotationDetail>>,
}

#[derive(Serialize, Clone)]
pub struct TermDetails {
    pub name: TermName,
    pub cv_name: CvName,
    pub interesting_parents: HashSet<String>,
    pub termid: TermId,
    #[serde(skip_serializing_if="Option::is_none")]
    pub definition: Option<TermDef>,
    pub is_obsolete: bool,
    pub rel_annotations: Vec<RelOntAnnotation>,
    pub genes_by_uniquename: HashMap<GeneUniquename, GeneShort>,
    pub alleles_by_uniquename: HashMap<AlleleUniquename, AlleleShort>,
    pub references_by_uniquename: HashMap<ReferenceUniquename, ReferenceShort>,
    pub terms_by_termid: HashMap<TermId, TermShort>,
}

#[derive(Serialize, Clone)]
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

#[derive(Serialize, Clone)]
pub struct OrthologAnnotation {
    pub gene_uniquename: GeneUniquename,
    pub ortholog_organism: OrganismShort,
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

#[derive(Serialize, Clone)]
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

#[derive(Serialize, Clone)]
pub struct Metadata {
    pub db_creation_datetime: String,
    pub export_prog_name: String,
    pub export_prog_version: String,
    pub gene_count: usize,
    pub term_count: usize,
}

#[derive(Serialize, Clone)]
pub struct SearchAPIMaps {
    pub termid_genes: HashMap<TermId, HashSet<GeneUniquename>>,
    pub term_name_genes: HashMap<TermName, HashSet<GeneUniquename>>,
    pub gene_summaries: Vec<GeneSummary>,
    pub term_summaries: HashSet<TermShort>,
}

#[derive(Serialize, Clone)]
pub struct WebData {
    pub genes: IdGeneMap,
    pub terms: IdTermDetailsMap,
    pub used_terms: IdTermDetailsMap,
    pub metadata: Metadata,
    pub references: IdReferenceMap,

    pub search_api_maps: SearchAPIMaps,
}
