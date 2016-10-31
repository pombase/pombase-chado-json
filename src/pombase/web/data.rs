use std::collections::HashMap;

pub type GeneUniquename = String;
pub type GeneName = String;
pub type TypeName = String;

#[derive(Serialize, Clone)]
pub struct GeneShort {
    pub uniquename: GeneUniquename,
    pub name: Option<GeneName>,
}

#[derive(Serialize, Clone)]
pub struct TranscriptShort {
    pub uniquename: TranscriptUniquename,
    //                pub exons: Vec<ExonShort>,
    //                pub utrs: Vec<UTRShort>,
}

pub type TermName = String;
pub type TermId = String;
pub type TermDef = String;

#[derive(Serialize, Clone)]
pub struct TermShort {
    pub name: TermName,
    pub termid: TermId,
    pub is_obsolete: bool,
}

#[derive(Serialize, Clone)]
pub struct PublicationShort {
    pub uniquename: String,
    pub title: Option<String>,
    pub citation: Option<String>,
}

type Evidence = String;

#[derive(Serialize, Clone)]
pub struct FeatureAnnotation {
    pub term: TermShort,
    pub evidence: Option<Evidence>,
    pub publication: Option<PublicationShort>,
    // only for genotype/phenotype annotation:
    pub genotype: Option<GenotypeAndAlleles>,
}

pub type TypeFeatureAnnotationMap =
    HashMap<TypeName, Vec<FeatureAnnotation>>;

#[derive(Serialize, Clone)]
pub struct GeneDetails {
    pub uniquename: GeneUniquename,
    pub name: Option<String>,
    pub transcripts: Vec<TranscriptShort>,
    pub annotations: TypeFeatureAnnotationMap,
}

pub type UniquenameGeneMap =
    HashMap<GeneUniquename, GeneDetails>;

pub type TranscriptUniquename = String;

#[derive(Serialize)]
pub struct TranscriptDetails {
    pub uniquename: TranscriptUniquename,
    pub name: Option<String>,
//    pub annotations: TypeFeatureAnnotationMap,
}

pub type UniquenameTranscriptMap =
    HashMap<TranscriptUniquename, TranscriptDetails>;

pub type GenotypeUniquename = String;

#[derive(Serialize)]
pub struct GenotypeDetails {
    pub uniquename: GenotypeUniquename,
    pub name: Option<String>,
    pub annotations: TypeFeatureAnnotationMap,
}

pub type UniquenameGenotypeMap =
    HashMap<GenotypeUniquename, GenotypeDetails>;

#[derive(Serialize, Clone)]
pub struct GenotypeAndAlleles {
    pub alleles: Vec<AlleleShort>,
}

pub type AlleleUniquename = String;

#[derive(Serialize, Clone)]
pub struct AlleleShort {
    pub uniquename: String,
    pub name: Option<String>,
    pub gene_uniquename: String,
}

pub type UniquenameAlleleShortMap =
    HashMap<AlleleUniquename, AlleleShort>;

#[derive(Serialize, Clone)]
pub struct TermAnnotation {
    pub gene: GeneShort,
    pub evidence: Option<Evidence>,
    pub publication: Option<PublicationShort>,
}

#[derive(Serialize)]
pub struct TermDetails {
    pub name: TermName,
    pub termid: TermId,
    pub definition: Option<TermDef>,
    pub is_obsolete: bool,
    pub annotations: Vec<TermAnnotation>,
}

impl Clone for TermDetails {
    fn clone(&self) -> TermDetails {
        TermDetails {
            name: self.name.clone(),
            termid: self.termid.clone(),
            definition: self.definition.clone(),
            is_obsolete: self.is_obsolete,
            annotations: self.annotations.clone(),
        }
    }
}

pub type IdGeneMap = HashMap<GeneUniquename, GeneDetails>;
pub type IdTermMap = HashMap<TermId, TermDetails>;

#[derive(Serialize)]
pub struct Metadata {
    pub db_creation_datetime: String,
    pub gene_count: usize,
    pub term_count: usize,
}
