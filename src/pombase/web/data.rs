use std::collections::HashMap;

type CvName = String;

pub type MiscExtRange = String;

pub type GeneUniquename = String;
pub type GeneName = String;
pub type TypeName = String;
pub type GeneProduct = String;

pub type TermName = String;
pub type TermId = String;
pub type TermDef = String;

pub type Evidence = String;

pub type TypeFeatureAnnotationMap =
    HashMap<TypeName, Vec<FeatureAnnotation>>;
pub type TypeReferenceAnnotationMap =
    HashMap<TypeName, Vec<ReferenceAnnotation>>;
pub type TypeInteractionAnnotationMap =
    HashMap<TypeName, Vec<InteractionAnnotation>>;

pub type UniquenameGeneMap =
    HashMap<GeneUniquename, GeneDetails>;

pub type TranscriptUniquename = String;

pub type UniquenameTranscriptMap =
    HashMap<TranscriptUniquename, TranscriptDetails>;

pub type GenotypeUniquename = String;

pub type UniquenameGenotypeMap =
    HashMap<GenotypeUniquename, GenotypeDetails>;

pub type AlleleUniquename = String;

pub type UniquenameAlleleShortMap =
    HashMap<AlleleUniquename, AlleleShort>;

pub type IdGeneMap = HashMap<GeneUniquename, GeneDetails>;
pub type IdGeneShortMap = HashMap<GeneUniquename, GeneShort>;
pub type IdTermMap = HashMap<TermId, TermDetails>;
pub type IdReferenceMap = HashMap<TermId, ReferenceDetails>;

include!(concat!(env!("OUT_DIR"), "/data_serde.rs"));

impl WebData {
    pub fn get_gene_summaries(&self) -> &IdGeneShortMap {
        &self.gene_summaries
    }
    pub fn get_genes(&self) -> &IdGeneMap {
        &self.genes
    }
    pub fn get_terms(&self) -> &IdTermMap {
        &self.terms
    }
    pub fn get_metadata(&self) -> &Metadata {
        &self.metadata
    }
    pub fn get_references(&self) -> &IdReferenceMap {
        &self.references
    }

}
