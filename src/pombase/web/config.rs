use std::collections::{HashMap, HashSet, BTreeMap};
use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;

use crate::types::*;
use serde_json;

use pombase_rc_string::RcString;

// configuration for extension display names and for the "Target of" section
#[derive(Deserialize, Clone, Debug)]
pub struct ExtensionDisplayNames {
    pub rel_name: RcString, // name of extension relation
    pub display_name: RcString, // text to display
    pub if_descendant_of: Option<RcString>, // None if applies to any extension
    pub reciprocal_display: Option<RcString>, // None if reciprocal shouldn't be displayed
}

// "interesting parents" are those stored in the JSON in the TermShort structs
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct InterestingParent {
    pub termid: RcString,
    pub rel_name: RcString,
}

// the order of relations within an extension:
#[derive(Deserialize, Clone, Debug)]
pub struct RelationOrder {
    // put the relations in this order in the displayed extensions:
    pub relation_order: Vec<RcString>,
    // except for these reactions which should always come last:
    pub always_last: Vec<RcString>,
}


#[derive(Deserialize, Clone, Debug)]
pub struct AncestorFilterCategory {
    pub display_name: RcString,
    // this category matches these terms and their descendants
    pub ancestors: Vec<RcString>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct FilterConfig {
    pub filter_name: String,
    pub display_name: String,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub term_categories: Vec<AncestorFilterCategory>,
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub extension_categories: Vec<AncestorFilterCategory>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct SplitByParentsConfig {
    pub termids: Vec<RcString>,
    pub display_name: RcString,
}

#[derive(Deserialize, Clone, Debug)]
pub struct ChromosomeConfig {
    pub name: RcString,
    // string to use for this chromosome in a file name, eg. "chromosome_II"
    // or "mitochondrial_chromosome"
    pub export_file_id: RcString,
    // string to use within files, eg. "II" or "mitochondrial"
    pub export_id: RcString,
    // eg. "Chromosome II" or "Mitochondrial chromosome"
    pub long_display_name: RcString,
    // eg. "II" or "Mitochondrial"
    pub short_display_name: RcString,
}

#[derive(Deserialize, Clone, Debug)]
pub struct CvSourceConfig {
    // a type name for the cvtermprop to display to the user
    pub display_name_prop: Option<RcString>,
    // the cvtermprop type name for the ID used for linking
    // or "ACCESSION" if the accession ID of the term should be used
    pub id_source: Option<RcString>,
}

pub type TargetRelationName = String;

#[derive(Deserialize, Clone, Debug)]
pub struct TargetOfConfig {
    // these priorities are used to order the list in the "Target of" section
    // and to filter the "Target of" summary
    // https://github.com/pombase/website/issues/299
    pub relation_priority: HashMap<TargetRelationName, u32>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct CvConfig {
    pub feature_type: RcString,
    pub display_name: Option<RcString>,
    // filtering configured per CV
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub filters: Vec<FilterConfig>,
    // config for splitting cv annotation tables into sub-sections
    // based on ancestry
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub split_by_parents: Vec<SplitByParentsConfig>,
    // relations to not show in the summary
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub summary_relations_to_hide: Vec<RcString>,
    // relations where the range is a gene ID to display like:
    //   has substrate pom1, cdc1 involved in negative regulation of ...
    // rather than as two lines
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub summary_relation_ranges_to_collect: Vec<RcString>,

    #[serde(default="SingleOrMultiLocusConfig::not_applicable")]
    pub single_or_multi_locus: SingleOrMultiLocusConfig,

    // the field to sort by
    #[serde(skip_serializing_if="Option::is_none")]
    pub sort_details_by: Option<Vec<RcString>>,

    // This is the configuration for the "Source" column, a map from
    // source name to config
    // See Disease association for an example.  If there is no config
    // there will be no Source column will be displayed
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    pub source_config: HashMap<RcString, CvSourceConfig>,
}

pub type ShortEvidenceCode = RcString;
pub type LongEvidenceCode = RcString;

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct ConfigOrganism {
    pub taxonid: OrganismTaxonId,
    pub genus: RcString,
    pub species: RcString,
    pub alternative_names: Vec<RcString>,
    pub assembly_version: Option<RcString>,
}

impl ConfigOrganism {
    pub fn full_name(&self) -> String {
        self.genus.clone() + "_" + self.species.as_str()
    }

    pub fn scientific_name(&self) -> String {
        self.genus.clone() + " " + self.species.as_str()
    }
}

#[derive(Deserialize, Clone, Debug)]
pub struct ViabilityTerms {
    pub viable: RcString,
    pub inviable: RcString,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct TermAndName {
    pub termid: RcString,
    pub name: RcString,
}

#[derive(Deserialize, Clone, Debug)]
pub struct ReferencePageConfig {
    pub triage_status_to_ignore: Vec<String>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct InterPro {
    pub dbnames_to_filter: Vec<RcString>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct ServerSubsetConfig {
    pub prefixes_to_remove: Vec<RcString>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct ServerConfig {
    pub subsets: ServerSubsetConfig,
    pub solr_url: String,
    pub close_synonym_boost: f32,
    pub distant_synonym_boost: f32,
    pub term_definition_boost: f32,
    pub django_url: String,
    pub cv_name_for_terms_search: String,
    pub gene_uniquename_re: String,
}

#[derive(Deserialize, Clone, Debug)]
pub struct EvidenceDetails {
    pub long: LongEvidenceCode,
    pub link: Option<RcString>,
}

pub type DatabaseName = RcString;
pub type DatabaseAliases = HashMap<DatabaseName, DatabaseName>;

#[derive(Deserialize, Clone, Debug)]
pub struct MacromolecularComplexesConfig {
    pub parent_complex_termid: RcString,
    pub excluded_terms: HashSet<RcString>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct RNAcentralConfig {
    // SO termids of RNA features to export
    pub export_so_ids: HashSet<RcString>,
}

#[derive(Deserialize, Clone, Debug, PartialEq)]
pub enum SingleOrMultiLocusConfig {
#[serde(rename = "single")]
    Single,
#[serde(rename = "multi")]
    Multi,
#[serde(rename = "na")]
    NotApplicable
}

impl SingleOrMultiLocusConfig {
    pub fn not_applicable() -> SingleOrMultiLocusConfig {
        SingleOrMultiLocusConfig::NotApplicable
    }
}

#[derive(Deserialize, Clone, Debug)]
pub struct ExportColumnConfig {
    pub name: RcString,
    pub display_name: RcString
}

#[derive(Deserialize, Clone, Debug)]
pub struct AnnotationSubsetConfig {
    pub term_ids: Vec<TermId>,
    pub file_name: RcString,
    pub columns: Vec<ExportColumnConfig>,
    #[serde(default="SingleOrMultiLocusConfig::not_applicable")]
    pub single_or_multi_locus: SingleOrMultiLocusConfig,
}

#[derive(Deserialize, Clone, Debug)]
pub struct GpadGpiConfig {
    // the term IDs of the three GO aspects
    pub go_aspect_terms: HashMap<String, TermId>,
    // Map a relation term name to a term ID, unless the term ID is None in
    // which case we skip writing this extension part
    pub extension_relation_mappings: HashMap<String, Option<TermId>>,
    // A map from the SO type of a transcript to the SO type of the gene is
    // derives from
    pub transcript_gene_so_term_map: HashMap<String, String>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct FileExportConfig {
    pub site_map_term_prefixes: Vec<RcString>,
    pub site_map_reference_prefixes: Vec<RcString>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub macromolecular_complexes: Option<MacromolecularComplexesConfig>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub rnacentral: Option<RNAcentralConfig>,
    pub annotation_subsets: Vec<AnnotationSubsetConfig>,
    pub gpad_gpi: GpadGpiConfig,
    // the reference to use for ND lines in GPAD/GAF output
    pub nd_reference: String,
}

#[derive(Deserialize, Clone, Debug)]
pub struct GeneResultVisAttrValueConfig {
    pub termid: Option<RcString>,
    pub name: RcString,
    pub bin_start: Option<usize>,
    pub bin_end: Option<usize>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct GeneResultVisColumnConfig {
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub attr_values: Vec<GeneResultVisAttrValueConfig>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct GeneResultsConfig {
    pub field_config: HashMap<RcString, GeneResultVisColumnConfig>,
    pub visualisation_field_names: Vec<RcString>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct SlimConfig {
    pub slim_display_name: RcString,
    pub cv_name: RcString,
    pub terms: Vec<TermAndName>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct SeqFeaturePageConfig {
    pub so_types_to_show: Vec<String>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct GeneExDatasetConfig {
    pub name: RcString,
    pub pubmed_id: RcString,
    pub level_type_termid: RcString,
    pub during_termid: RcString,
    pub scale: RcString,
}

#[derive(Deserialize, Clone, Debug)]
pub struct GeneExpressionConfig {
    pub datasets: Vec<GeneExDatasetConfig>,
}

#[derive(Deserialize, Clone, Debug)]
pub struct Config {
    pub database_name: RcString,
    pub database_long_name: RcString,
    pub database_citation: RcString,
    pub funder: RcString,
    pub site_description: RcString,
    pub load_organism_taxonid: Option<OrganismTaxonId>,
    pub base_url: RcString,
    pub doc_page_aliases: HashMap<String, String>,
    pub organisms: Vec<ConfigOrganism>,
    pub api_seq_chunk_sizes: Vec<usize>,
    pub sequence_feature_page: SeqFeaturePageConfig,
    pub extension_display_names: Vec<ExtensionDisplayNames>,
    pub extension_relation_order: RelationOrder,
    pub evidence_types: HashMap<ShortEvidenceCode, EvidenceDetails>,
    pub cv_config: HashMap<CvName, CvConfig>,
    pub target_of_config: TargetOfConfig,
// when creating a TermShort struct, for each of these termids if the term has
// an "interesting parent" using the given rel_name, we store it in the
// interesting_parents field of the TermShort
    pub interesting_parents: Vec<InterestingParent>,
    pub viability_terms: ViabilityTerms,
    // slim sets by slim name:
    pub slims: HashMap<RcString, SlimConfig>,
    pub reference_page_config: ReferencePageConfig,
    pub interpro: InterPro,
    pub server: ServerConfig,
    pub extra_database_aliases: DatabaseAliases,
    pub chromosomes: Vec<ChromosomeConfig>,
    pub gene_results: GeneResultsConfig,
    pub ortholog_taxonids: HashSet<u32>,
    pub file_exports: FileExportConfig,
    pub gene_expression: GeneExpressionConfig,
}

impl Config {

    pub fn read(config_file_name: &str) -> Config {
        let file = match File::open(config_file_name) {
            Ok(file) => file,
            Err(err) => {
                panic!("Failed to read {}: {}\n", config_file_name, err)
            }
        };
        let reader = BufReader::new(file);

        match serde_json::from_reader(reader) {
            Ok(config) => config,
            Err(err) => {
                panic!("failed to parse {}: {}", config_file_name, err)
            },
        }
    }

    pub fn cv_config_by_name(&self, cv_name: &str) -> CvConfig {
        if let Some(config) = self.cv_config.get(cv_name) {
            config.clone()
        } else {
            let empty_cv_config = 
                CvConfig {
                    feature_type: "".into(),
                    display_name: Some("".into()),
                    single_or_multi_locus: SingleOrMultiLocusConfig::NotApplicable,
                    filters: vec![],
                    split_by_parents: vec![],
                    summary_relations_to_hide: vec![],
                    summary_relation_ranges_to_collect: vec![],
                    sort_details_by: None,
                    source_config: HashMap::new(),
                };

            if cv_name.starts_with("extension:") {
                if cv_name.ends_with(":gene") {
                    CvConfig {
                        feature_type: "gene".into(),
                        ..empty_cv_config
                    }
                } else {
                    CvConfig {
                        feature_type: "genotype".into(),
                        ..empty_cv_config
                    }
                }
            } else {
                CvConfig {
                    feature_type: "gene".into(),
                    ..empty_cv_config
                }
            }
        }
    }

    pub fn organism_by_taxonid(&self, lookup_taxonid: u32) -> Option<ConfigOrganism> {
        for org in &self.organisms {
            if org.taxonid == lookup_taxonid {
                return Some(org.clone());
            }
        }

        None
    }

    pub fn load_organism(&self) -> Option<ConfigOrganism> {
        if let Some(load_organism_taxonid) = self.load_organism_taxonid {
            let org = self.organism_by_taxonid(load_organism_taxonid);
            if org.is_none() {
                panic!("can't find configuration for load_organism_taxonid: {}",
                       load_organism_taxonid);
            }
            org
        } else {
            None
        }
    }

    pub fn find_chromosome_config<'a>(&'a self, chromosome_name: &str)
                                      -> &'a ChromosomeConfig
    {
        for chr_config in &self.chromosomes {
            if chr_config.name == chromosome_name {
                return chr_config;
            }
        }
        panic!("can't find chromosome configuration for {}", &chromosome_name);
    }
}

pub const POMBASE_ANN_EXT_TERM_CV_NAME: &str = "PomBase annotation extension terms";
pub const ANNOTATION_EXT_REL_PREFIX: &str = "annotation_extension_relation-";

pub enum FeatureRelAnnotationType {
    Interaction,
    Ortholog,
    Paralog,
}
pub struct FeatureRelConfig {
    pub rel_type_name: &'static str,
    pub annotation_type: FeatureRelAnnotationType,
}
pub const FEATURE_REL_CONFIGS: [FeatureRelConfig; 4] =
    [
        FeatureRelConfig {
            rel_type_name: "interacts_physically",
            annotation_type: FeatureRelAnnotationType::Interaction,
        },
        FeatureRelConfig {
            rel_type_name: "interacts_genetically",
            annotation_type: FeatureRelAnnotationType::Interaction,
        },
        FeatureRelConfig {
            rel_type_name: "orthologous_to",
            annotation_type: FeatureRelAnnotationType::Ortholog,
        },
        FeatureRelConfig {
            rel_type_name: "paralogous_to",
            annotation_type: FeatureRelAnnotationType::Paralog,
        },
    ];

// relations to use when copy annotation to parents (ie. adding the
// annotation of child terms to parents)
pub const DESCENDANT_REL_NAMES: [&str; 7] =
    ["is_a", "part_of", "regulates", "positively_regulates", "negatively_regulates",
     "has_part", "output_of"];
// only consider has_part relations for these ontologies:
pub const HAS_PART_CV_NAMES: [&str; 1] = ["fission_yeast_phenotype"];

// number of genes before (and after) to add to the gene_neighbourhood field
pub const GENE_NEIGHBOURHOOD_DISTANCE: usize = 5;

pub const TRANSCRIPT_FEATURE_TYPES: [&str; 8] =
    ["snRNA", "rRNA", "mRNA", "snoRNA", "ncRNA", "tRNA", "pseudogenic_transcript",
     "transcript"];
pub const TRANSCRIPT_PART_TYPES: [&str; 4] =
    ["five_prime_UTR", "exon", "pseudogenic_exon", "three_prime_UTR"];
// any feature with a type not in this list or in the two TRANSCRIPT lists above
// will be stored in the other_features map
pub const HANDLED_FEATURE_TYPES: [&str; 7] =
    ["gene", "pseudogene", "intron", "genotype", "allele", "chromosome", "polypeptide"];


#[derive(Deserialize, Clone, Debug)]
pub struct DocConfig {
    pub pages: BTreeMap<RcString, RcString>,
}

impl DocConfig {

    pub fn read(doc_config_file_name: &str) -> DocConfig {
        let file = match File::open(doc_config_file_name) {
            Ok(file) => file,
            Err(err) => {
                panic!("Failed to read {}: {}\n", doc_config_file_name, err)
            }
        };
        let reader = BufReader::new(file);

        match serde_json::from_reader(reader) {
            Ok(config) => config,
            Err(err) => {
                panic!("failed to parse {}: {}", doc_config_file_name, err)
            },
        }
    }
}


pub struct GoEcoMapping {
    mapping: HashMap<(String, String), String>,
}

impl GoEcoMapping {
    pub fn read(file_name: &str) -> Result<GoEcoMapping, std::io::Error> {
        let file = match File::open(file_name) {
            Ok(file) => file,
            Err(err) => {
                panic!("Failed to read {}: {}\n", file_name, err)
            }
        };
        let reader = BufReader::new(file);

        let mut mapping = HashMap::new();

        for line_result in reader.lines() {
            match line_result {
                Ok(line) => {
                    if line.starts_with("#") {
                        continue;
                    }
                    let parts: Vec<&str> = line.split('\t').collect();
                    mapping.insert((String::from(parts[0]), String::from(parts[1])),
                                   String::from(parts[2]));
                },
                Err(err) => return Err(err)
            };
        }

        Ok(GoEcoMapping {
            mapping
        })
    }

    pub fn lookup_default(&self, go_evidence_code: &str) -> Option<String> {
        self.mapping.get(&(String::from(go_evidence_code), String::from("Default")))
            .map(|s| String::from(s))
    }

    pub fn lookup_with_go_ref(&self, go_evidence_code: &str, go_ref: &str)
                              -> Option<String>
    {
        self.mapping.get(&(String::from(go_evidence_code), String::from(go_ref)))
            .map(|s| String::from(s))
    }
}
