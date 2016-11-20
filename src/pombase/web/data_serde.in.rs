use std::rc::Rc;

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

#[derive(Serialize, Clone)]
pub struct TranscriptShort {
    pub uniquename: TranscriptUniquename,
    //                pub exons: Vec<ExonShort>,
    //                pub utrs: Vec<UTRShort>,
}

#[derive(Serialize, Clone)]
pub struct TermShort {
    pub name: TermName,
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
    pub annotations: TypeReferenceAnnotationMap,
    pub interaction_annotations: TypeInteractionAnnotationMap,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
}

#[derive(Serialize, Clone)]
pub struct ReferenceAnnotation {
    pub gene: GeneShort,
    pub term: TermShort,
    pub evidence: Option<Evidence>,
    pub extension: Vec<ExtPart>,
    // only for genotype/phenotype annotation:
    pub genotype: Option<GenotypeAndAlleles>,
    pub is_not: bool,
}

#[derive(Serialize, Clone)]
pub struct FeatureAnnotation {
    pub term: TermShort,
    pub extension: Vec<ExtPart>,
    pub evidence: Option<Evidence>,
    pub reference: Option<ReferenceShort>,
    // only for genotype/phenotype annotation:
    pub genotype: Option<GenotypeAndAlleles>,
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
    pub annotations: TypeFeatureAnnotationMap,
    pub interaction_annotations: TypeInteractionAnnotationMap,
    pub ortholog_annotations: Vec<OrthologAnnotation>,
    pub paralog_annotations: Vec<ParalogAnnotation>,
}

#[derive(Serialize)]
pub struct TranscriptDetails {
    pub uniquename: TranscriptUniquename,
    pub name: Option<String>,
//    pub annotations: TypeFeatureAnnotationMap,
}

#[derive(Serialize)]
pub struct GenotypeDetails {
    pub uniquename: GenotypeUniquename,
    pub name: Option<String>,
    pub background: Option<String>,
    pub annotations: TypeFeatureAnnotationMap,
}

#[derive(Serialize, Clone)]
pub struct GenotypeShort {
    pub uniquename: GenotypeUniquename,
    pub name: Option<String>,
    pub background: Option<String>,
}

#[derive(Serialize, Clone)]
pub struct GenotypeAndAlleles {
    pub genotype: GenotypeShort,
    pub alleles: Vec<AlleleShort>,
}

#[derive(Serialize, Clone)]
pub struct AlleleShort {
    pub uniquename: String,
    pub name: Option<String>,
    pub gene_uniquename: String,
}

#[derive(Serialize, Clone)]
pub struct TermAnnotation {
    pub term: TermShort,
    pub gene: GeneShort,
    pub extension: Vec<ExtPart>,
    pub evidence: Option<Evidence>,
    pub reference: Option<ReferenceShort>,
    pub is_not: bool,
}

pub type TermAnnotationKey = String;

pub type TermAnnotationMap = HashMap<TermAnnotationKey, Vec<Rc<TermAnnotation>>>;

#[derive(Serialize, Clone)]
pub struct TermDetails {
    pub name: TermName,
    pub cv_name: CvName,
    pub termid: TermId,
    pub definition: Option<TermDef>,
    pub is_obsolete: bool,
    pub genes: Vec<GeneShort>,
    pub annotations: TermAnnotationMap,
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

#[derive(Serialize)]
pub struct Metadata {
    pub db_creation_datetime: String,
    pub export_prog_name: String,
    pub export_prog_version: String,
    pub gene_count: usize,
    pub term_count: usize,
}

#[derive(Serialize)]
pub struct WebData {
    pub genes: IdGeneMap,
    pub gene_summaries: IdGeneShortMap,
    pub terms: IdTermDetailsMap,
    pub used_terms: IdTermDetailsMap,
    pub metadata: Metadata,
    pub references: IdReferenceMap,
}
