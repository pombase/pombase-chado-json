use core::panic;
use std::hash::Hash;
use std::num::NonZeroUsize;
use std::rc::Rc;
use std::collections::{BTreeMap, BTreeSet};
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::sync::{Arc, RwLock};

use pombase_gocam::overlaps::{GoCamNodeOverlap, find_activity_overlaps, find_chemical_overlaps};
use pombase_gocam_process::find_holes;
use regex::Regex;

use std::collections::{HashMap, HashSet};

use pombase_gocam::{GoCamModel, GoCamNode};

use crate::bio::pdb_reader::{PDBGeneEntryMap, PDBRefEntryMap};
use crate::bio::protein_view::make_protein_view_data_map;

use crate::db::{raw::*, ChadoQueries};

use crate::gene_history::GeneHistoryMap;
use crate::types::*;
use crate::data_types::*;
use crate::uniprot::UniProtDataMap;
use crate::web::data::*;
use crate::web::config::*;
use crate::web::util::cmp_str_dates;
use crate::utils::join;

use crate::bio::util::{compare_ext_part_with_config, rev_comp};

use flexstr::{SharedStr as FlexStr, shared_str as flex_str, ToSharedStr, shared_fmt as flex_fmt};

use crate::interpro::DomainData;

lazy_static! {
    static ref ISO_DATE_RE: Regex =
        Regex::new(r"^(\d\d\d\d)-(\d\d)-(\d\d)($|\s)").unwrap();
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum AddToHashFlag {
    AlleleCommentRefs,
    AlleleSynonymRefs,
}

fn make_organism(rc_organism: &Rc<ChadoOrganism>) -> ConfigOrganism {
    let mut maybe_taxonid: Option<u32> = None;
    for prop in rc_organism.organismprops.borrow().iter() {
        if prop.prop_type.name == "taxon_id" {
            maybe_taxonid = Some(prop.value.parse().unwrap());
        }
    }
    ConfigOrganism {
        taxonid: maybe_taxonid.unwrap(),
        genus: rc_organism.genus.clone(),
        species: rc_organism.species.clone(),
        alternative_names: vec![],
        assembly_version: None,
    }
}

type TermShortOptionMap = HashMap<TermId, Option<TermShort>>;

type GenotypeInteractionUniquename = FlexStr;

pub struct WebDataBuild<'a> {
    raw: &'a Raw,
    domain_data: DomainData,

    uniprot_data: Option<UniProtDataMap>,
    rnacentral_data: Option<RNAcentralAnnotations>,
    all_gene_history: Option<GeneHistoryMap>,
    pdb_gene_entry_map: Option<PDBGeneEntryMap>,
    pdb_ref_entry_map: Option<PDBRefEntryMap>,
    chado_queries: ChadoQueries,
    orcid_name_map: HashMap<CuratorOrcid, FlexStr>,
    gocam_models: Vec<GoCamModel>,
    config: &'a Config,

    genes: UniquenameGeneMap,
    genotypes: DisplayUniquenameGenotypeMap,
    genotypes_of_alleles: HashMap<AlleleUniquename, HashSet<GenotypeUniquename>>,
    genotype_backgrounds: HashMap<GenotypeUniquename, FlexStr>,
    alleles: UniquenameAlleleDetailsMap,
    transcripts: UniquenameTranscriptMap,
    other_features: UniquenameFeatureShortMap,
    terms: TermIdDetailsMap,
    chromosomes: ChrNameDetailsMap,
    references: UniquenameReferenceMap,
    all_ont_annotations: HashMap<TermId, Vec<OntAnnotationId>>,
    all_not_ont_annotations: HashMap<TermId, Vec<OntAnnotationId>>,

    protein_complexes: HashMap<ProteinComplexUniquename, ProteinComplexDetails>,
    genes_of_complexes: HashMap<ProteinComplexUniquename, HashSet<GeneUniquename>>,

    // map from term name to term ID (ie "nucleus" -> "GO:0005634")
    term_ids_by_name: HashMap<FlexStr, TermId>,

    genes_of_transcripts: HashMap<FlexStr, FlexStr>,
    transcripts_of_polypeptides: HashMap<FlexStr, FlexStr>,
    parts_of_transcripts: HashMap<FlexStr, Vec<FeatureShort>>,
    genes_of_alleles: HashMap<FlexStr, FlexStr>,
    loci_of_genotypes: HashMap<FlexStr, HashMap<FlexStr, GenotypeLocus>>,
    genotype_display_uniquenames: HashMap<GenotypeUniquename, GenotypeDisplayUniquename>,

    // maps used to collect the genotype interaction part from the feature_relationship table
    genotype_interaction_double_mutant: HashMap<GenotypeInteractionUniquename, GenotypeUniquename>,
    genotype_interaction_genotype_a: HashMap<GenotypeInteractionUniquename, GenotypeUniquename>,
    genotype_interaction_genotype_b: HashMap<GenotypeInteractionUniquename, GenotypeUniquename>,

    // a map from IDs of terms from the "PomBase annotation extension terms" cv
    // to a Vec of the details of each of the extension
    parts_of_extensions: HashMap<TermId, Vec<ExtPart>>,

    base_term_of_extensions: HashMap<TermId, TermId>,

    // a set of child terms for each term from the cvtermpath table
    children_by_termid: HashMap<TermId, HashSet<TermId>>,
    dbxrefs_of_features: HashMap<FlexStr, HashSet<FlexStr>>,

    possible_interesting_parents: HashSet<InterestingParent>,

    recent_references: RecentReferences,
    all_community_curated: Vec<ReferenceShort>,
    all_admin_curated: Vec<ReferenceShort>,

    gene_expression_measurements: GeneExDataSetMeasurements,

    term_subsets: IdTermSubsetMap,
    gene_subsets: IdGeneSubsetMap,

    annotation_details: IdOntAnnotationDetailMap,

    ont_annotations: Vec<OntAnnotation>,

    physical_interaction_annotations: HashSet<InteractionAnnotation>,
    genetic_interaction_annotations: HashMap<GeneticInteractionKey, Vec<GeneticInteractionDetail>>,

    gocam_summaries: GoCamSummaryMap,

    protein_complex_data: ProteinComplexData,

    gocam_overlaps: Vec<GoCamNodeOverlap>,
    gocam_overlaps_merge_by_chemical: Vec<GoCamNodeOverlap>,

    gocam_holes: Vec<GoCamNode>,

    // transcripts with overlapping exons
    transcript_frameshifts_to_check: Vec<(TranscriptUniquename, usize)>,

    // for passing to GoCamModel::genes_in_model()
    pro_term_to_gene: HashMap<String, String>,
}

#[allow(clippy::type_complexity)]
fn get_maps() ->
    (HashMap<FlexStr, ReferenceShortOptionMap>,
     HashMap<FlexStr, GeneShortOptionMap>,
     HashMap<FlexStr, GenotypeShortMap>,
     HashMap<FlexStr, AlleleShortMap>,
     HashMap<FlexStr, TranscriptDetailsOptionMap>,
     HashMap<GeneUniquename, TermShortOptionMap>)
{
    (HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new(),
     HashMap::new(), HashMap::new())
}

// returns (Option<expression>, Option<promoter gene>, Option<exogenous promoter>)
fn get_feat_rel_expression_and_promoter(feature: &Feature,
                                        feature_relationship: &FeatureRelationship)
   -> (Option<FlexStr>, Option<FlexStr>)
{
    for feature_prop in feature.featureprops.borrow().iter() {
        if feature_prop.prop_type.name == "allele_type"
            && let Some(ref value) = feature_prop.value
                && value == "deletion" {
                    return (Some("Null".into()), None);
                }
    }

    let mut maybe_expression = None;
    let mut maybe_promoter_gene = None;

    for rel_prop in feature_relationship.feature_relationshipprops.borrow().iter() {
        match rel_prop.prop_type.name.as_str() {
            "expression" => maybe_expression = rel_prop.value.clone(),
            "promoter_gene" => maybe_promoter_gene = rel_prop.value.clone(),
            _ => (),
        }
    }

    (maybe_expression, maybe_promoter_gene)
}

fn get_feat_rel_prop_value(prop_name: &FlexStr,
                           feature_relationship: &FeatureRelationship) -> Option<FlexStr> {
    for rel_prop in feature_relationship.feature_relationshipprops.borrow().iter() {
        if rel_prop.prop_type.name == *prop_name {
            return rel_prop.value.clone();
        }
    }

    None
}

fn reference_has_annotation(reference_details: &ReferenceDetails) -> bool {
    !reference_details.cv_annotations.is_empty() ||
        !reference_details.physical_interactions.is_empty() ||
        !reference_details.genetic_interactions.is_empty() ||
        !reference_details.ortholog_annotations.is_empty() ||
        !reference_details.paralog_annotations.is_empty()
}

fn is_gene_type(feature_type_name: &FlexStr) -> bool {
    feature_type_name == "gene" || feature_type_name == "pseudogene"
}

lazy_static! {
    static ref BAD_GENOTYPE_NAME_CHARS_RE: Regex =
        Regex::new(r"[% /&;?\\]").unwrap();
}

pub fn make_genotype_display_uniquename(loci: &[GenotypeLocus],
                                  allele_map: &UniquenameAlleleDetailsMap) -> FlexStr {
    let mut locus_display_names: Vec<String> =
        loci.iter().map(|locus| {
            let mut allele_display_names: Vec<String> =
                locus.expressed_alleles.iter().map(|expressed_allele| {
                    let allele_details = allele_map.get(&expressed_allele.allele_uniquename).unwrap();
                    let mut encoded_name_and_type =
                        allele_details.encoded_name_and_type.to_string();
                    if allele_details.allele_type != "deletion" {
                        if encoded_name_and_type == "unnamed-unrecorded-unrecorded" {
                            encoded_name_and_type =
                                format!("{}-{}", allele_details.gene.uniquename,
                                        encoded_name_and_type);
                        }
                        if let Some(ref expression) = expressed_allele.expression {
                            encoded_name_and_type += &format!("-expression-{}", expression.to_lowercase());
                        }
                    }
                    encoded_name_and_type
                }).collect();
            allele_display_names.sort();
            allele_display_names.join("/")
        }).collect();

    locus_display_names.sort();

    let joined_alleles = locus_display_names.join("  ");

    let clean_display_name =
        BAD_GENOTYPE_NAME_CHARS_RE.replace_all(&joined_alleles, "_");
    clean_display_name.to_shared_str()
}

fn make_genotype_display_name(loci: &[GenotypeLocus],
                              allele_map: &UniquenameAlleleDetailsMap) -> FlexStr {
    let mut locus_display_names: Vec<String> =
        loci.iter().map(|locus| {
            let mut allele_display_names: Vec<String> =
                locus.expressed_alleles.iter().map(|expressed_allele| {
                    let allele_details = allele_map.get(&expressed_allele.allele_uniquename).unwrap();
                    let allele_short: AlleleShort = allele_details.into();
                    let allele_display_name = allele_short.display_name();

                    if let Some(ref expression) = expressed_allele.expression {
                        if expression.starts_with("Not assayed") ||
                            allele_details.allele_type == "deletion" &&
                            expression == "Null"
                        {
                           allele_display_name.to_string()
                        } else {
                           format!("{}[{}]", allele_display_name, expression)
                        }

                    } else {
                        allele_display_name.to_string()
                    }
                }).collect();
            allele_display_names.sort();
            allele_display_names.join("/")
        }).collect();

    locus_display_names.sort();

    let joined_alleles = locus_display_names.join(" ");

    joined_alleles.to_shared_str()
}

fn make_phase(feature_loc: &Featureloc) -> Option<Phase> {
    if let Some(phase) = feature_loc.phase {
        match phase {
            0 => Some(Phase::Zero),
            1 => Some(Phase::One),
            2 => Some(Phase::Two),
            _ => panic!(),
        }
    } else {
        None
    }
}

fn make_location(chromosome_map: &ChrNameDetailsMap,
                 feat: &Feature) -> Option<ChromosomeLocation> {
    let feature_locs = feat.featurelocs.borrow();
    match feature_locs.first() {
        Some(feature_loc) => {
            let start_pos =
                if feature_loc.fmin + 1 >= 1 {
                    (feature_loc.fmin + 1) as usize
                } else {
                    panic!("start_pos less than 1");
                };
            let end_pos =
                if feature_loc.fmax >= 1 {
                    feature_loc.fmax as usize
                } else {
                    panic!("start_end less than 1");
                };
            let feature_uniquename = &feature_loc.srcfeature.uniquename;
            let chr_short = make_chromosome_short(chromosome_map, feature_uniquename);
            Some(ChromosomeLocation {
                chromosome_name: chr_short.name,
                start_pos,
                end_pos,
                strand: match feature_loc.strand {
                    1 => Strand::Forward,
                    -1 => Strand::Reverse,
                    _ => panic!(),
                },
                phase: make_phase(feature_loc),
            })
        },
        None => None,
    }
}

fn get_loc_residues(chr: &ChromosomeDetails,
                    loc: &ChromosomeLocation) -> Residues {
    let start = loc.start_pos - 1;
    let end = loc.end_pos;
    let residues: Residues = chr.residues[start..end].into();
    if loc.strand == Strand::Forward {
        residues
    } else {
        rev_comp(&residues)
    }
}

fn make_feature_short(chromosome_map: &ChrNameDetailsMap, feat: &Feature) -> FeatureShort {
    let maybe_loc = make_location(chromosome_map, feat);
    if let Some(loc) = maybe_loc {
        if let Some(chr) = chromosome_map.get(&loc.chromosome_name) {
            let residues = get_loc_residues(chr, &loc);
            let feature_type = match &feat.feat_type.name as &str {
                "five_prime_UTR" => FeatureType::FivePrimeUtr,
                "pseudogenic_exon" | "exon" => FeatureType::Exon,
                "three_prime_UTR" => FeatureType::ThreePrimeUtr,
                "dg_repeat" => FeatureType::DGRepeat,
                "dh_repeat" => FeatureType::DHRepeat,
                "gap" => FeatureType::Gap,
                "gene_group" => FeatureType::GeneGroup,
                "long_terminal_repeat" => FeatureType::LongTerminalRepeat,
                "low_complexity_region" => FeatureType::LowComplexityRegion,
                "LTR_retrotransposon" => FeatureType::LTRRetrotransposon,
                "mating_type_region" => FeatureType::MatingTypeRegion,
                "nuclear_mt_pseudogene" => FeatureType::NuclearMtPseudogene,
                "origin_of_replication" => FeatureType::OriginOfReplication,
                "polyA_signal_sequence" => FeatureType::PolyASignalSequence,
                "polyA_site" => FeatureType::PolyASite,
                "promoter" => FeatureType::Promoter,
                "region" => FeatureType::Region,
                "regional_centromere" => FeatureType::RegionalCentromere,
                "regional_centromere_central_core" => FeatureType::RegionalCentromereCentralCore,
                "regional_centromere_inner_repeat_region" => FeatureType::RegionalCentromereInnerRepeatRegion,
                "repeat_region" => FeatureType::RepeatRegion,
                "TR_box" => FeatureType::TRBox,
                "antisense_RNA" => FeatureType::AntisenseRNA,
                "sncRNA" => FeatureType::AntisenseRNA,
                "lncRNA" => FeatureType::LncRNA,
                "guide_RNA" => FeatureType::GuideRNA,
                "SNP" => FeatureType::SNP,
                _ => panic!("can't handle feature type: {}", feat.feat_type.name),
            };
            FeatureShort {
                feature_type,
                uniquename: feat.uniquename.clone(),
                name: feat.name.clone(),
                location: loc,
                residues,
            }
        } else {
            panic!("can't find chromosome {}", loc.chromosome_name);
        }
    } else {
        panic!("{} has no featureloc", feat.uniquename);
    }
}

pub fn make_chromosome_short<'a>(chromosome_map: &'a ChrNameDetailsMap,
                                 chromosome_name: &'a FlexStr) -> ChromosomeShort {
    if let Some(chr) = chromosome_map.get(chromosome_name) {
        chr.make_chromosome_short()
    } else {
        panic!("can't find chromosome: {}", chromosome_name);
    }
}


fn make_reference_short(reference_map: &UniquenameReferenceMap,
                        reference_uniquename: &FlexStr) -> Option<ReferenceShort> {
    if reference_uniquename == "null" {
        None
    } else {
        let reference_details = reference_map.get(reference_uniquename)
            .unwrap_or_else(|| panic!("missing reference in make_reference_short(): {}",
                            reference_uniquename));

        let reference_short = ReferenceShort::from_reference_details(reference_details);

        Some(reference_short)
    }
}

lazy_static! {
    static ref PROMOTER_RE: Regex = Regex::new(r"^(?P<gene>.*)-promoter$").unwrap();
    static ref PREFIX_AND_ID_RE: Regex =
        Regex::new(r"^(?P<prefix>\S+):(?P<id>\S+)$").unwrap();
    static ref TRANSCRIPT_ID_RE: Regex = Regex::new(r"^(?P<gene>.*)\.(?P<suffix>\d+)$").unwrap();
}

// Some ancestor terms are useful in the web code.  This function uses the Config and returns
// the terms that might be useful.
fn get_possible_interesting_parents(config: &Config) -> HashSet<InterestingParent> {
    let mut ret = HashSet::new();

    for parent_conf in &config.interesting_parents {
        ret.insert(parent_conf.clone());
    }

    for ext_conf in &config.extension_display_names {
        if let Some(ref conf_termid) = ext_conf.if_descendant_of {
            ret.insert(InterestingParent {
                termid: conf_termid.clone(),
                rel_name: flex_str!("is_a"),
            });
        }
    }

    let add_to_set = |set: &mut HashSet<_>, termid: FlexStr| {
        for rel_name in &DESCENDANT_REL_NAMES {
            set.insert(InterestingParent {
                termid: termid.to_owned(),
                rel_name: flex_str!(*rel_name),
            });
        }
    };

    for slim_config in config.slims.values() {
        for go_slim_conf in &slim_config.terms {
            add_to_set(&mut ret, go_slim_conf.termid.clone());
        }
    }

    for field_name in &config.gene_results.visualisation_field_names {
        if let Some(column_conf) = config.gene_results.field_config.get(field_name) {

            for attr_value_conf in &column_conf.attr_values {
                if let Some(ref termid) = attr_value_conf.termid {
                    add_to_set(&mut ret, termid.clone());
                }
            }
        } else {
            panic!["can't find field configuration for {}", field_name];
        }
    }

    ret.insert(InterestingParent {
        termid: config.viability_terms.viable.clone(),
        rel_name: flex_str!("is_a"),
    });
    ret.insert(InterestingParent {
        termid: config.viability_terms.inviable.clone(),
        rel_name: flex_str!("is_a"),
    });

    let add_filter_ancestor =
        |set: &mut HashSet<_>, category: &AncestorFilterCategory, cv_name: &FlexStr| {
            for ancestor in &category.ancestors {
                for config_rel_name in &DESCENDANT_REL_NAMES {
                    if *config_rel_name == "has_part" &&
                        !HAS_PART_CV_NAMES.contains(cv_name) {
                            continue;
                        }
                    set.insert(InterestingParent {
                        termid: ancestor.clone(),
                        rel_name: flex_str!(*config_rel_name),
                    });
                }
            }
        };

    for (cv_name, conf) in &config.cv_config {
        for filter in &conf.filters {
            for category in &filter.term_categories {
                add_filter_ancestor(&mut ret, category, cv_name);
            }
            for category_configs in config.extension_categories.values() {
                for cat_config in category_configs {
                    add_filter_ancestor(&mut ret, cat_config, cv_name);
                }
            }
        }

        for split_by_parent_config in &conf.split_by_parents {
            for ancestor in &split_by_parent_config.termids {
                let ancestor_termid =
                    if let Some(without_prefix) = ancestor.strip_prefix("NOT ") {
                        without_prefix.to_shared_str()
                    } else {
                        ancestor.clone()
                    };
                ret.insert(InterestingParent {
                    termid: ancestor_termid,
                    rel_name: "is_a".into(),
                });
            }
        }
    }

    for mod_group in &config.protein_feature_view.modification_groups {
        ret.insert(InterestingParent {
            termid: mod_group.termid.clone(),
            rel_name: "is_a".into(),
        });
    }

    ret
}

const MAX_RECENT_REFS: usize = 20;

fn make_recently_added(references_map: &UniquenameReferenceMap,
                       all_ref_uniquenames: &[FlexStr]) -> Vec<ReferenceShort> {
    let mut date_sorted_pub_uniquenames = all_ref_uniquenames.to_owned();

    {
        let ref_added_date_cmp =
            |ref_uniquename1: &ReferenceUniquename, ref_uniquename2: &ReferenceUniquename| {
                let ref1 = references_map.get(ref_uniquename1).unwrap();
                let ref2 = references_map.get(ref_uniquename2).unwrap();

                if let Some(ref ref1_added_date) = ref1.canto_added_date {
                    if let Some(ref ref2_added_date) = ref2.canto_added_date {
                        cmp_str_dates(ref1_added_date, ref2_added_date).reverse()
                    } else {
                        Ordering::Less
                    }
                } else if ref2.canto_added_date.is_some() {
                    Ordering::Greater
                } else {
                    Ordering::Equal
                }
            };

        date_sorted_pub_uniquenames.sort_by(ref_added_date_cmp);
    }

    let recently_added_iter =
        date_sorted_pub_uniquenames.iter().take(MAX_RECENT_REFS);

    let mut recently_added: Vec<ReferenceShort> = vec![];

    for ref_uniquename in recently_added_iter {
        let ref_short_maybe = make_reference_short(references_map, ref_uniquename);
        if let Some(ref_short) = ref_short_maybe {
            recently_added.push(ref_short);
        }
    }

    recently_added
}

fn make_canto_curated(references_map: &UniquenameReferenceMap,
                      all_ref_uniquenames: &[FlexStr])
                      -> (Vec<ReferenceShort>, Vec<ReferenceShort>, Vec<ReferenceShort>,
                          Vec<ReferenceShort>) {
    let mut sorted_pub_uniquenames: Vec<ReferenceUniquename> =
        all_ref_uniquenames.iter()
        .filter(|ref_uniquename| {
            let reference = references_map.get(*ref_uniquename).unwrap();
            reference.canto_approved_date.is_some()
        })
        .cloned()
        .collect();

    {
        let pub_date_cmp =
            |ref_uniquename1: &ReferenceUniquename, ref_uniquename2: &ReferenceUniquename| {
                let ref1 = references_map.get(ref_uniquename1).unwrap();
                let ref2 = references_map.get(ref_uniquename2).unwrap();

                // use first approval date, but fall back to the most recent approval date

                let ref1_date =
                    ref1.canto_first_approved_date.as_ref()
                    .unwrap_or_else(|| ref1.canto_session_submitted_date.as_ref().unwrap());
                let ref2_date = ref2.canto_first_approved_date.as_ref()
                    .unwrap_or_else(|| ref2.canto_session_submitted_date.as_ref().unwrap());

                cmp_str_dates(ref2_date, ref1_date)
            };

        sorted_pub_uniquenames.sort_by(pub_date_cmp);
    }

    let mut recent_admin_curated = vec![];
    let mut recent_community_curated = vec![];
    let mut all_community_curated = vec![];
    let mut all_admin_curated = vec![];

    let ref_uniquename_iter = sorted_pub_uniquenames.iter();

    for ref_uniquename in ref_uniquename_iter {
        let reference = references_map.get(ref_uniquename).unwrap();

        let ref_short = make_reference_short(references_map, ref_uniquename).unwrap();
        if reference.canto_curator_role.to_lowercase() == "community" {
            all_community_curated.push(ref_short.clone());
            if recent_community_curated.len() <= MAX_RECENT_REFS {
                recent_community_curated.push(ref_short);
            }
        } else {
            all_admin_curated.push(ref_short.clone());
            if recent_admin_curated.len() <= MAX_RECENT_REFS {
                recent_admin_curated.push(ref_short);
            }
        }
    }

    (recent_admin_curated, recent_community_curated, all_community_curated,
     all_admin_curated)
}
// return true iff the gene has a frameshifted warning annotation
fn has_annotated_frame_shift(gene_details: &GeneDetails) -> bool {
    let Some(term_annotations) = gene_details.cv_annotations.get("warning")
    else {
        return false;
    };

    for term_annotation in term_annotations {
        if term_annotation.term == "PBO:0000108" {
            return true;
        }
    }

    false
}

// returns a list of IDs for transcripts that have potential frameshifts, with the position
// in the transcript
fn add_introns_to_transcript(chromosome: &ChromosomeDetails,
                             transcript_uniquename: &FlexStr,
                             strand: Strand,
                             parts: &mut Vec<FeatureShort>)
   -> Option<usize>
{
    let mut new_parts: Vec<FeatureShort> = vec![];

    let mut maybe_frameshift_pos = None;

    let mut intron_count = 0;

    for part in parts.drain(0..) {
        let mut maybe_new_intron = None;

        if let Some(prev_part) = new_parts.last() {
            let intron_start = prev_part.location.end_pos + 1;
            let intron_end = part.location.start_pos - 1;

            let overlap = part.location.start_pos as i64 - prev_part.location.end_pos as i64 - 1;

            if overlap <= 0 {
                if overlap == 0 && (prev_part.feature_type != FeatureType::Exon || part.feature_type != FeatureType::Exon) {
                    // no gap between coding exon and 5'/3' UTR
                } else if overlap == -1 || overlap == -2 || overlap == -3 {
                    // Probably a overlap that represents a frameshift in the reference.
                    // We'll check later for "warning, frameshifted".
                    // See:
                    // https://github.com/pombase/curation/issues/1453#issuecomment-303214177
                    maybe_frameshift_pos = Some(intron_start);
                } else {
                    println!("no gap between exons at {}..{} in {}", intron_start, intron_end,
                             transcript_uniquename);
                }
            } else {

                let new_intron_loc = ChromosomeLocation {
                    chromosome_name: prev_part.location.chromosome_name.clone(),
                    start_pos: intron_start,
                    end_pos: intron_end,
                    strand: prev_part.location.strand,
                    phase: None,
                };

                let intron_residues = get_loc_residues(chromosome, &new_intron_loc);

                let intron_type =
                    if prev_part.feature_type == FeatureType::Exon &&
                    part.feature_type == FeatureType::Exon {
                        FeatureType::CdsIntron
                    } else if prev_part.feature_type == FeatureType::FivePrimeUtr {
                        FeatureType::FivePrimeUtrIntron
                    } else {
                        FeatureType::ThreePrimeUtrIntron
                    };
                maybe_new_intron = Some(FeatureShort {
                    feature_type: intron_type,
                    uniquename: flex_str!("placeholder"), // we set this later
                    name: None,
                    location: new_intron_loc,
                    residues: intron_residues,
                });

                intron_count += 1;
            }
        }

        if let Some(new_intron) = maybe_new_intron {
            new_parts.push(new_intron);
        }

        new_parts.push(part);
    }

    let mut count_for_uniquename =
        if strand == Strand::Forward {
            1
        } else {
            intron_count
        };

    for part in &mut new_parts {
        if part.feature_type.is_any_intron_type() {
           let intron_uniquename =
              format!("{}:intron:{}", transcript_uniquename, count_for_uniquename);

           part.uniquename = intron_uniquename.into();

           if strand == Strand::Forward {
               count_for_uniquename += 1;
           } else {
               count_for_uniquename -= 1;
           }
        }
    }

    *parts = new_parts;

    maybe_frameshift_pos
}

fn validate_transcript_parts(transcript_uniquename: &FlexStr, parts: &[FeatureShort]) {
    let mut seen_exon = false;
    for part in parts {
        if part.feature_type == FeatureType::Exon {
            seen_exon = true;
            break;
        }
    }
    if !seen_exon {
        panic!("transcript has no exons: {}", transcript_uniquename);
    }

    if parts[0].feature_type != FeatureType::Exon {
        for i in 1..parts.len() {
            let part = &parts[i];
            if part.feature_type == FeatureType::Exon {
                let last_utr_before_exons = &parts[i-1];

                let first_exon = &parts[i];

                if last_utr_before_exons.location.end_pos + 1 != first_exon.location.start_pos {
                    println!("{} and exon don't meet up: {} at pos {}",
                             last_utr_before_exons.feature_type, transcript_uniquename,
                             last_utr_before_exons.location.end_pos);
                }

                break;
            } else if part.location.strand == Strand::Forward {
                if part.feature_type != FeatureType::FivePrimeUtr {
                    println!("{:?}", parts);
                    panic!("wrong feature type '{}' before exons in {}",
                           part.feature_type, transcript_uniquename);
                }
            } else if part.feature_type != FeatureType::ThreePrimeUtr {
                println!("{:?}", parts);
                panic!("wrong feature type '{}' after exons in {}",
                       part.feature_type, transcript_uniquename);
            }
        }
    }

    let last_part = parts.last().unwrap();

    if last_part.feature_type != FeatureType::Exon {
        for i in (0..parts.len()-1).rev() {
            let part = &parts[i];
            if part.feature_type == FeatureType::Exon {
                let first_utr_after_exons = &parts[i+1];

                let last_exon = &parts[i];

                if last_exon.location.end_pos + 1 != first_utr_after_exons.location.start_pos {
                    println!("{} and exon don't meet up: {} at pos {}",
                             first_utr_after_exons.feature_type, transcript_uniquename,
                             first_utr_after_exons.location.end_pos);
                }

                break;
            } else if part.location.strand == Strand::Forward {
                if part.feature_type != FeatureType::ThreePrimeUtr {
                    panic!("wrong feature type '{}' before exons in {}",
                           part.feature_type, transcript_uniquename);
                }
            } else if part.feature_type != FeatureType::FivePrimeUtr {
                panic!("wrong feature type '{}' after exons in {}",
                       part.feature_type, transcript_uniquename);
            }
        }
    }
}

fn set_has_protein_features(genes: &mut UniquenameGeneMap, protein_view_data: &HashMap<GeneUniquename, ProteinViewData>) {
    for (gene_uniquename, protein_feature_data) in protein_view_data {
        if let Some(gene_details) = genes.get_mut(gene_uniquename) {
            let mut has_protein_features = false;

            for track in &protein_feature_data.tracks {
                if !track.features.is_empty() {
                    has_protein_features = true;
                    break;
                }
            }

            gene_details.has_protein_features = has_protein_features;
        }
    }
}


lazy_static! {
    static ref CONDITION_DETAIL_RE: Regex =
        Regex::new(r"^([A-Z][\w]+:\d\d\d+)(?:\((.+)\))?$").unwrap();
}

// parse something like "FYECO:0000005(32C)" or "FYECO:0000329(2% (v/v))" into a term ID
// and a String containers the details from the brackets
// parse "FYECO:0000005" into (TermId, None)
fn parse_condition_with_detail(condition_string: &str) -> (TermId, Option<String>) {
    let captures = CONDITION_DETAIL_RE.captures(condition_string).unwrap();

    let term_id = captures.get(1).as_ref().unwrap().as_str().into();

    if let Some(detail) = captures.get(2) {
        (term_id, Some(detail.as_str().into()))
    } else {
        (term_id, None)
    }
}

fn get_cumulative_annotation_type_counts(annotation_type_counts: StatsIntegerTable)
     -> StatsIntegerTable
{
    let mut cumulative_counts = annotation_type_counts;

    for row_idx in 1..cumulative_counts.data.len() {
        let prev = cumulative_counts.data[row_idx-1].1.clone();
        let current = &mut cumulative_counts.data[row_idx].1;
        for col_idx in 0..cumulative_counts.header.len() {
            current[col_idx] += prev[col_idx];
        }
    }

    cumulative_counts
}

type CuratedStats =(Vec<(DateString, Vec<usize>)>, Vec<(DateString, Vec<usize>)>);

impl <'a> WebDataBuild<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(raw: &'a Raw,
               domain_data: DomainData,
               uniprot_data: Option<UniProtDataMap>,
               rnacentral_data: Option<RNAcentralAnnotations>,
               all_gene_history: Option<GeneHistoryMap>,
               pdb_gene_entry_map: Option<PDBGeneEntryMap>,
               pdb_ref_entry_map: Option<PDBRefEntryMap>,
               chado_queries: ChadoQueries,
               orcid_name_map: HashMap<CuratorOrcid, FlexStr>,
               gocam_models: Vec<GoCamModel>,
               config: &'a Config) -> WebDataBuild<'a>
    {
        WebDataBuild {
            raw,
            domain_data,

            uniprot_data,
            rnacentral_data,
            all_gene_history,
            pdb_gene_entry_map,
            pdb_ref_entry_map,
            chado_queries,
            orcid_name_map,
            gocam_models,
            config,

            genes: BTreeMap::new(),
            genotypes: HashMap::new(),
            genotypes_of_alleles: HashMap::new(),
            genotype_backgrounds: HashMap::new(),
            genotype_interaction_genotype_a: HashMap::new(),
            genotype_interaction_genotype_b: HashMap::new(),
            genotype_interaction_double_mutant: HashMap::new(),
            alleles: BTreeMap::new(),
            transcripts: HashMap::new(),
            other_features: HashMap::new(),

            protein_complexes: HashMap::new(),
            genes_of_complexes: HashMap::new(),
            terms: HashMap::new(),
            chromosomes: BTreeMap::new(),
            references: HashMap::new(),
            all_ont_annotations: HashMap::new(),
            all_not_ont_annotations: HashMap::new(),
            physical_interaction_annotations: HashSet::new(),
            genetic_interaction_annotations: HashMap::new(),
            recent_references: RecentReferences {
                admin_curated: vec![],
                community_curated: vec![],
                pubmed: vec![],
            },
            all_community_curated: vec![],
            all_admin_curated: vec![],

            term_ids_by_name: HashMap::new(),

            genes_of_transcripts: HashMap::new(),
            transcripts_of_polypeptides: HashMap::new(),
            parts_of_transcripts: HashMap::new(),
            genes_of_alleles: HashMap::new(),
            loci_of_genotypes: HashMap::new(),
            genotype_display_uniquenames: HashMap::new(),

            parts_of_extensions: HashMap::new(),

            base_term_of_extensions: HashMap::new(),

            children_by_termid: HashMap::new(),
            dbxrefs_of_features: HashMap::new(),

            possible_interesting_parents: get_possible_interesting_parents(config),

            term_subsets: HashMap::new(),
            gene_subsets: HashMap::new(),

            annotation_details: HashMap::new(),

            ont_annotations: vec![],

            gene_expression_measurements: HashMap::new(),

            gocam_summaries: HashMap::new(),

            gocam_overlaps: vec![],
            gocam_overlaps_merge_by_chemical: vec![],
            gocam_holes: vec![],

            protein_complex_data: HashMap::new(),

            transcript_frameshifts_to_check: vec![],

            pro_term_to_gene: HashMap::new(),
       }
    }

    fn add_ref_to_hash(&self,
                       seen_references: &mut HashMap<FlexStr, ReferenceShortOptionMap>,
                       identifier: &FlexStr,
                       maybe_reference_uniquename: &Option<ReferenceUniquename>) {
        if let Some(reference_uniquename) = maybe_reference_uniquename
            && reference_uniquename != "null" {
                seen_references
                    .entry(identifier.into())
                    .or_default()
                    .insert(reference_uniquename.clone(), None);
            }
    }


    fn add_gene_to_hash(&self,
                        seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                        identifier: &FlexStr,
                        other_gene_uniquename: &GeneUniquename) {
        if !self.genes.contains_key(other_gene_uniquename) {
            panic!("{}", other_gene_uniquename);
        }
        seen_genes
            .entry(identifier.clone())
            .or_default()
            .insert(other_gene_uniquename.clone(), None);
    }

    fn add_genotype_to_hash(&self,
                            seen_genotypes: &mut HashMap<FlexStr, GenotypeShortMap>,
                            seen_alleles: &mut HashMap<FlexStr, AlleleShortMap>,
                            seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                            seen_references: &mut HashMap<FlexStr, ReferenceShortOptionMap>,
                            identifier: &FlexStr,
                            genotype_uniquename: &FlexStr) {
        let genotype_short = self.make_genotype_short(genotype_uniquename);
        for locus in &genotype_short.loci {
            for expressed_allele in &locus.expressed_alleles {
                self.add_allele_to_hash(HashSet::new(), seen_alleles, seen_genes,
                                        seen_references, identifier,
                                        &expressed_allele.allele_uniquename);
            }
        }

        seen_genotypes
            .entry(identifier.clone())
            .or_default()
            .insert(genotype_uniquename.clone(), genotype_short);
    }

    fn add_allele_to_hash(&self,
                          flags: HashSet<AddToHashFlag>,
                          seen_alleles: &mut HashMap<FlexStr, AlleleShortMap>,
                          seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                          seen_references: &mut HashMap<FlexStr, ReferenceShortOptionMap>,
                          identifier: &FlexStr,
                          allele_uniquename: &AlleleUniquename) -> AlleleShort {
        let allele_short = self.make_allele_short(allele_uniquename);
        {
            if flags.contains(&AddToHashFlag::AlleleCommentRefs) {
                for comment_details in &allele_short.comments {
                    self.add_ref_to_hash(seen_references, identifier,
                                         &comment_details.reference);
                }
            }
            if flags.contains(&AddToHashFlag::AlleleSynonymRefs) {
                for synonym_details in &allele_short.synonyms {
                    self.add_ref_to_hash(seen_references, identifier,
                                         &synonym_details.reference);
                }
            }

            let allele_gene_uniquename = &allele_short.gene_uniquename;
            self.add_gene_to_hash(seen_genes, identifier, allele_gene_uniquename);
            seen_alleles
                .entry(identifier.clone())
                .or_default()
                .insert(allele_uniquename.clone(), allele_short.clone());
        }
        allele_short
    }

    fn add_transcript_to_hashes(&self,
                              seen_transcripts: &mut HashMap<FlexStr, TranscriptDetailsOptionMap>,
                              seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                              identifier: &FlexStr,
                              transcript_uniquename: &TranscriptUniquename) {
        if let Some(transcript_details) = self.transcripts.get(transcript_uniquename) {
            seen_transcripts
                .entry(identifier.clone())
                .or_default()
                .insert(transcript_uniquename.clone(), None);
            self.add_gene_to_hash(seen_genes, identifier,
                                  &transcript_details.gene_uniquename);
        } else {
            panic!("internal error, can't find transcript {}",
                   transcript_uniquename);
        }
    }

    fn add_term_to_hash(&self,
                        seen_terms: &mut HashMap<TermId, TermShortOptionMap>,
                        identifier: &FlexStr,
                        other_termid: &TermId) {
        seen_terms
            .entry(identifier.clone())
            .or_default()
            .insert(other_termid.clone(), None);
    }

    fn get_gene<'b>(&'b self, gene_uniquename: &'b FlexStr) -> &'b GeneDetails {
        if let Some(gene_details) = self.genes.get(gene_uniquename) {
            gene_details
        } else {
            panic!("can't find GeneDetails for gene uniquename {}", gene_uniquename)
        }
    }

    fn get_gene_mut<'b>(&'b mut self, gene_uniquename: &'b FlexStr) -> &'b mut GeneDetails {
        if let Some(gene_details) = self.genes.get_mut(gene_uniquename) {
            gene_details
        } else {
            panic!("can't find GeneDetails for gene uniquename {}", gene_uniquename)
        }
    }

    fn make_gene_short(&self, gene_uniquename: &FlexStr) -> GeneShort {
        self.get_gene(gene_uniquename).into()
    }

    fn make_gene_summary(&self, gene_uniquename: &FlexStr) -> GeneSummary {
        let gene_details = self.get_gene(gene_uniquename);
        let synonyms =
            gene_details.synonyms.iter()
            .filter(|synonym| synonym.synonym_type == "exact")
            .map(|synonym| synonym.name.clone())
            .collect::<Vec<FlexStr>>();
        let ortholog_ids =
            gene_details.ortholog_annotations.iter()
            .map(|ortholog_annotation| {
                let ortholog_uniquename = ortholog_annotation.ortholog_uniquename.clone();
                let ortholog_gene_summary = &self.genes.get(&ortholog_uniquename).unwrap();
                let maybe_secondary_identifier = ortholog_gene_summary.secondary_identifier.clone();
                let maybe_agr_identifier = ortholog_gene_summary.agr_identifier.clone();
                let maybe_ortholog_name = ortholog_gene_summary.name.clone();

                IdNameAndOrganism {
                    identifier: ortholog_uniquename,
                    secondary_identifier: maybe_secondary_identifier,
                    agr_identifier: maybe_agr_identifier,
                    name: maybe_ortholog_name,
                    taxonid: ortholog_annotation.ortholog_taxonid,
                }
            })
            .collect::<Vec<IdNameAndOrganism>>();

        GeneSummary {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
            uniprot_identifier: gene_details.uniprot_identifier.clone(),
            secondary_identifier: gene_details.secondary_identifier.clone(),
            synonyms,
            orthologs: ortholog_ids,
            feature_type: gene_details.feature_type.clone(),
            taxonid: gene_details.taxonid,
            transcript_count: gene_details.transcripts.len(),
            location: gene_details.location.clone(),
        }
    }

    fn make_api_gene_summary(&self, gene_uniquename: &FlexStr) -> APIGeneSummary {
        let gene_details = self.get_gene(gene_uniquename);
        let synonyms =
            gene_details.synonyms.iter()
            .filter(|synonym| synonym.synonym_type == "exact")
            .map(|synonym| synonym.name.clone())
            .collect::<Vec<FlexStr>>();
        let mut coding_exon_count = 0;
        let mut five_prime_exon_count = 0;
        let mut three_prime_exon_count = 0;

        if let Some(transcript_uniquename) = gene_details.transcripts.first() {
            let transcript = self.transcripts
                .get(transcript_uniquename)
                .unwrap_or_else(|| panic!("internal error, can't find transcript details for {}",
                                transcript_uniquename));

            for part in &transcript.parts {
                match part.feature_type {
                    FeatureType::Exon => coding_exon_count += 1,
                    FeatureType::FivePrimeUtr => five_prime_exon_count += 1,
                    FeatureType::ThreePrimeUtr => three_prime_exon_count += 1,
                    _ => (),
                }
            }
        }

        let mut ortholog_taxonids = HashSet::new();
        for ortholog_annotation in &gene_details.ortholog_annotations {
            ortholog_taxonids.insert(ortholog_annotation.ortholog_taxonid);
        }

        let transcript_details = gene_details.transcripts
            .iter()
            .map(|transcript_uniquename| {
                self.transcripts.get(transcript_uniquename)
                    .unwrap_or_else(|| panic!("internal error, failed to find transcript: {}",
                                    transcript_uniquename))
                    .clone()
            }).collect::<Vec<_>>();

        let pdb_ids = gene_details.pdb_entries.iter()
            .map(|pdb_entry| pdb_entry.pdb_id.clone())
            .collect();

        let mut gocam_ids = HashSet::new();
        let mut enables_gocam_activity_ids = HashSet::new();

        let uniquename_with_prefix =
            format!("{}:{}", self.config.database_name, gene_uniquename);

        for model in &self.gocam_models {
            if model.genes_in_model().contains(&uniquename_with_prefix) {
                gocam_ids.insert(model.id().into());
            }
            if model.genes_enabling_activities().contains_key(&uniquename_with_prefix) {
                enables_gocam_activity_ids.insert(model.id().into());
            }
        }

        let disordered_regions_count = gene_details.disordered_region_coords.len();
        let low_complexity_regions_count = gene_details.low_complexity_region_coords.len();
        let mut coiled_coil_percent = 0;
        let mut disordered_percent = 0;
        let mut low_complexity_percent = 0;

        if let Some(transcript) = transcript_details.first()
            && let Some(ref protein) = transcript.protein {
                coiled_coil_percent =
                    100usize * gene_details.coiled_coil_aa_count() / protein.sequence_length();
                disordered_percent =
                    100usize * gene_details.disordered_aa_count() / protein.sequence_length();
                low_complexity_percent =
                    100usize * gene_details.low_complexity_aa_count() / protein.sequence_length();
            };

        APIGeneSummary {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
            feature_type: gene_details.feature_type.clone(),
            uniprot_identifier: gene_details.uniprot_identifier.clone(),
            exact_synonyms: synonyms,
            dbxrefs: gene_details.dbxrefs.clone(),
            pdb_ids,
            gocam_ids,
            enables_gocam_activity_ids,
            location: gene_details.location.clone(),
            transcripts: transcript_details,
            tm_domain_count: gene_details.tm_domain_coords.len(),
            coiled_coil_count: gene_details.coiled_coil_coords.len(),
            coiled_coil_percent,
            disordered_regions_count,
            disordered_percent,
            low_complexity_regions_count,
            low_complexity_percent,
            coding_exon_count,
            five_prime_exon_count,
            three_prime_exon_count,
            transcript_count: gene_details.transcripts.len(),
            ortholog_taxonids,
        }
    }

    fn make_term_short(&self, termid: &FlexStr) -> TermShort {
        if let Some(term_details) = self.terms.get(termid) {
            TermShort::from_term_details(term_details)
        } else {
            panic!("can't find TermDetails for termid: {}", termid)
        }
    }

    fn add_characterisation_status(&mut self, gene_uniquename: &FlexStr,
                                   cvterm_name: &FlexStr) {
        let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        gene_details.characterisation_status = Some(cvterm_name.clone());
    }

    fn add_gene_product(&mut self, gene_uniquename: &FlexStr, new_product: &FlexStr) {
        let gene_details = self.get_gene_mut(gene_uniquename);
        if let Some(ref product) = gene_details.product {
            gene_details.product = Some(flex_fmt!("{} OR {}", product, new_product));
        } else {
            gene_details.product = Some(new_product.clone());
        }
    }

    fn add_name_description(&mut self, gene_uniquename: &FlexStr, name_description: &FlexStr) {
        let gene_details = self.get_gene_mut(gene_uniquename);
        if !gene_details.name_descriptions.contains(name_description) {
          gene_details.name_descriptions.push(name_description.into());
        }
    }

    fn annotation_from_template(&self,
                                extension_relation_order: &RelationOrder,
                                cvterm: &Cvterm,
                                annotation_template: OntAnnotationDetail)
        -> OntAnnotationDetail
    {

        let extension_parts =
            match self.parts_of_extensions.get(&cvterm.termid()) {
                Some(parts) => parts.clone(),
                None => vec![],
            };

        let mut new_extension = extension_parts;

        let mut existing_extensions = annotation_template.extension.clone();
        new_extension.append(&mut existing_extensions);

        let compare_ext_part_func =
            |e1: &ExtPart, e2: &ExtPart| {
                compare_ext_part_with_config(extension_relation_order, e1, e2)
            };

        new_extension.sort_by(compare_ext_part_func);

        OntAnnotationDetail {
            extension: new_extension,
            .. annotation_template
        }
    }

    // return the gene of a single locus
    fn gene_from_genotype(&self, genotype_uniquename: &GenotypeUniquename)
       -> GeneUniquename
    {
        let genotype_display_uniquename =
            self.genotype_display_uniquenames.get(genotype_uniquename).unwrap();

        let genotype =
            self.genotypes.get(genotype_display_uniquename)
            .unwrap_or_else(|| panic!("internal error: can't find genotype {}", genotype_uniquename));

        let loci = &genotype.loci;

        if loci.len() > 1 {
            panic!("genotype has multiple loci");
        }

        let expressed_alleles = &loci[0].expressed_alleles;

        let allele_uniquename = &expressed_alleles[0].allele_uniquename;

        let allele = self.alleles.get(allele_uniquename).unwrap();

        allele.gene.uniquename.clone()
    }

    fn parse_extension_prop(&self, annotation_termid: &TermId, extension_string: &str) -> Vec<ExtPart> {
        let ext: Vec<Vec<CantoExtPart>> =
            serde_json::from_str(extension_string)
            .unwrap_or_else(|_| panic!("failed to parse Canto extension from property: {}", extension_string));

        if ext.len() > 1 {
            eprintln!("\
currently we can't handle extensions with multiple parts for single allele
phenotypes, so just the first part of this extension will be used:
{}", extension_string);
        }

        if ext.is_empty() {
            return vec![];
        }

        let inner = ext.first().unwrap();

        if inner.is_empty() {
            return vec![];
        }

        let to_parsed = |part: &CantoExtPart| {
            let ext_range =
                if let Some(ref range_type) = part.range_type {
                    match range_type.as_str() {
                       "Ontology" => ExtRange::Term(part.range_value.clone()),
                       "Gene" => ExtRange::Gene(part.range_value.clone()),
                       _ => ExtRange::Misc(part.range_value.clone()),
                    }
                } else {
                    let value = part.range_value.clone();
                    if part.range_value.contains(":") {
                        ExtRange::Term(value)
                    } else {
                        ExtRange::Misc(value)
                    }
                };

            let rel_type_display_name =
                self.get_ext_rel_display_name(annotation_termid, &part.relation);

            ExtPart {
                rel_type_name: part.relation.clone(),
                rel_type_display_name,
                rel_type_id: None,
                ext_range,
            }
        };

        inner.iter().map(|part: &CantoExtPart| to_parsed(part)).collect()
    }

    fn add_genetic_interaction(&mut self, genotype_interaction_feature: &Feature,
                               extension_relation_order: &RelationOrder,
                               cvterm: &Cvterm,
                               annotation_template: OntAnnotationDetail) {
        let genotype_interaction_uniquename = &genotype_interaction_feature.uniquename;
        let genotype_a_uniquename =
            self.genotype_interaction_genotype_a.get(genotype_interaction_uniquename)
                .unwrap_or_else(|| panic!("can't find genotype_a of {}",
                                 genotype_interaction_uniquename));
        let genotype_a_display_uniquename =
            self.genotype_display_uniquenames.get(genotype_a_uniquename).unwrap().clone();

        let genotype_b_uniquename =
            self.genotype_interaction_genotype_b.get(genotype_interaction_uniquename)
                .unwrap_or_else(|| panic!("can't find genotype_b of {}",
                                 genotype_interaction_uniquename));
        let genotype_b_display_uniquename =
            self.genotype_display_uniquenames.get(genotype_b_uniquename).unwrap().clone();

        let double_mutant_uniquename =
            self.genotype_interaction_double_mutant.get(genotype_interaction_uniquename);
        let double_mutant_genotype_display_uniquename =
            if let Some(double_genotype_uniquename) = double_mutant_uniquename {
                let double_genotype_display_uniquename =
                    self.genotype_display_uniquenames.get(double_genotype_uniquename).unwrap().clone();
                Some(double_genotype_display_uniquename)
            } else {
                None
            };

        let ont_annotation_detail =
            self.annotation_from_template(extension_relation_order,
                                          cvterm, annotation_template);

        let mut interaction_type = None;
        let mut interaction_note = None;
        let mut rescued_phenotype_termid = None;
        let mut rescued_phenotype_extension_value = None;
        let mut annotation_date = None;

        for prop in genotype_interaction_feature.featureprops.borrow().iter() {
            match prop.prop_type.name.as_str() {
                "interaction_type" => interaction_type = prop.value.clone(),
                "interaction_note" => interaction_note = prop.value.clone(),
                "interaction_rescued_phenotype_id" => rescued_phenotype_termid = prop.value.clone(),
                "interaction_rescued_phenotype_extension" =>
                    rescued_phenotype_extension_value = prop.value.clone(),
                "annotation_date" => annotation_date = prop.value.clone(),
                _ => (),
            }
        }

        let termid = cvterm.termid();
        let base_termid = self.base_term_of_extensions.get(&termid);

        let double_mutant_phenotype_termid =
            if let Some(base_termid) = base_termid {
                Some(base_termid.clone())
            } else {
                Some(termid.clone())
            };

        let rescued_phenotype_termid =
            if let Some(rescued_phenotype_termid) = rescued_phenotype_termid {
                if let Some(base_rescued_phenotype_termid) =
                    self.base_term_of_extensions.get(&rescued_phenotype_termid)
                {
                    Some(base_rescued_phenotype_termid.clone())
                } else {
                    Some(rescued_phenotype_termid)
                }
            } else {
                None
            };

        let ext_filter = |ext_part: &ExtPart| {
            ext_part.rel_type_name == "has_penetrance" ||
            ext_part.rel_type_name == "has_severity"
        };

        let filter_interaction_extension = |mut ext: Vec<ExtPart>| {
            ext.drain(0..).filter(ext_filter)
                .collect()
        };

        let rescued_phenotype_extension =
            if let (Some(termid), Some(ext)) =
                (&rescued_phenotype_termid, &rescued_phenotype_extension_value)
            {
                let ext = self.parse_extension_prop(termid, ext);
                // See: https://github.com/pombase/pombase-chado/issues/1114
                filter_interaction_extension(ext)
            } else {
                vec![]
            };

        let interaction_type =
            interaction_type.unwrap_or_else(|| panic!("interaction_type missing for {}",
                                             genotype_interaction_uniquename));

        let gene_a_uniquename = self.gene_from_genotype(genotype_a_uniquename);
        let gene_b_uniquename = self.gene_from_genotype(genotype_b_uniquename);

        let interaction_key =
            GeneticInteractionKey {
                gene_a_uniquename,
                interaction_type,
                gene_b_uniquename,
            };

        let double_mutant_extension =
            filter_interaction_extension(ont_annotation_detail.extension);

        let interaction_annotation =
            GeneticInteractionDetail {
                genotype_a_uniquename: Some(genotype_a_display_uniquename),
                genotype_b_uniquename: Some(genotype_b_display_uniquename),
                reference_uniquename: ont_annotation_detail.reference,
                double_mutant_phenotype_termid,
                double_mutant_extension,
                double_mutant_genotype_display_uniquename,
                rescued_phenotype_termid,
                rescued_phenotype_extension,
                interaction_note,
                throughput: ont_annotation_detail.throughput,
                source_database: Some(self.config.database_name.clone()),
                annotation_date,
            };

        self.genetic_interaction_annotations
            .entry(interaction_key)
            .or_default()
            .push(interaction_annotation);
    }

    fn add_annotation(&mut self, extension_relation_order: &RelationOrder,
                      cvterm: &Cvterm, is_not: bool,
                      annotation_template: OntAnnotationDetail) {
        let termid =
            match self.base_term_of_extensions.get(&cvterm.termid()) {
                Some(base_termid) => base_termid.clone(),
                None => cvterm.termid(),
            };

        let ont_annotation_detail =
            self.annotation_from_template(extension_relation_order,
                                          cvterm, annotation_template);

        let annotation_map = if is_not {
            &mut self.all_not_ont_annotations
        } else {
            &mut self.all_ont_annotations
        };

        let entry = annotation_map.entry(termid);
        entry.or_insert_with(Vec::new).push(ont_annotation_detail.id);

        self.annotation_details.insert(ont_annotation_detail.id,
                                          ont_annotation_detail);
    }

    fn process_dbxrefs(&mut self) {
        let mut map = HashMap::new();

        for feature_dbxref in &self.raw.feature_dbxrefs {
            let feature = &feature_dbxref.feature;
            let dbxref = &feature_dbxref.dbxref;

            map.entry(feature.uniquename.clone())
                .or_insert_with(HashSet::new)
                .insert(dbxref.identifier());
        }

        self.dbxrefs_of_features = map;
    }

    fn get_extra_gene_pubs(&self, publication: &Publication)
         -> (Vec<GeneUniquename>, Vec<GeneUniquename>)
    {
        let mut pubmed_keyword_genes = vec![];
        let mut extra_genes = vec![];
        for feat_pub in publication.feature_publications.borrow().iter() {
            for prop in feat_pub.feature_pubprops.borrow().iter() {
                if prop.prop_type.name == "feature_pub_source"
                    && let Some(ref prop_value) = prop.value {
                        let feature = &feat_pub.feature;
                        if feature.feat_type.name == "gene" {
                            let gene_uniquename = feature.uniquename.clone();
                            if prop_value == "pubmed_keyword" {
                                pubmed_keyword_genes.push(gene_uniquename);
                            } else {
                                extra_genes.push(gene_uniquename);
                            }
                        }
                    }
            }
        }
        (pubmed_keyword_genes, extra_genes)
    }

    fn process_references(&mut self) {
        let parse_annotation_curator =
            |value: &str, prop_type: &str| {
                let curator: AnnotationCurator = serde_json::from_str(value)
                    .unwrap_or_else(|_| panic!("failed to parse {} pupprop: {}", prop_type, value));
                curator
            };

        let author_re = Regex::new(r"^(?P<f>[^,]+),.*$").unwrap();
        let publication_date_re = Regex::new(r"^(.* )?(?P<y>\d\d\d\d)$").unwrap();
        let approved_date_re = Regex::new(r"^(?P<date>\d\d\d\d-\d\d-\d\d).*").unwrap();

        for rc_publication in &self.raw.publications {
            let reference_uniquename = &rc_publication.uniquename;

            if reference_uniquename.to_lowercase() == "null" {
                continue;
            }

            let mut pubmed_authors: Option<FlexStr> = None;
            let mut pubmed_publication_date: Option<FlexStr> = None;
            let mut pubmed_entrez_date: Option<FlexStr> = None;
            let mut pubmed_abstract: Option<FlexStr> = None;
            let mut pubmed_doi: Option<FlexStr> = None;
            let mut non_pubmed_authors: Option<FlexStr> = None;
            let mut non_pubmed_abstract: Option<FlexStr> = None;
            let mut canto_session_key: Option<FlexStr> = None;
            let mut canto_annotation_status: Option<FlexStr> = None;
            let mut canto_triage_status: Option<FlexStr> = None;
            let mut canto_curator_role = self.config.database_name.clone();
            let mut canto_curator_name: Option<FlexStr> = None;
            let mut canto_first_approved_date: Option<FlexStr> = None;
            let mut canto_approved_date: Option<FlexStr> = None;
            let mut canto_approver_orcid: Option<FlexStr> = None;
            let mut canto_added_date: Option<FlexStr> = None;
            let mut canto_session_submitted_date: Option<FlexStr> = None;
            let mut annotation_curators = vec![];

            let mut file_curator_role = None;
            let mut file_curator_name: Option<FlexStr> = None;
            let mut annotation_file_curators = vec![];

            for prop in rc_publication.publicationprops.borrow().iter() {
                match &prop.prop_type.name as &str {
                    "pubmed_publication_date" =>
                        pubmed_publication_date = Some(prop.value.clone()),
                    "pubmed_entrez_date" =>
                        pubmed_entrez_date = Some(prop.value.clone()),
                    "pubmed_authors" =>
                        pubmed_authors = Some(prop.value.clone()),
                    "pubmed_abstract" =>
                        pubmed_abstract = Some(prop.value.clone()),
                    "pubmed_doi" =>
                        pubmed_doi = Some(prop.value.clone()),
                    "authors" =>
                        non_pubmed_authors = Some(prop.value.clone()),
                    "abstract" =>
                        non_pubmed_abstract = Some(prop.value.clone()),
                    "canto_session" =>
                        canto_session_key = Some(prop.value.clone()),
                    "canto_annotation_status" =>
                        canto_annotation_status = Some(prop.value.clone()),
                    "canto_triage_status" =>
                        canto_triage_status = Some(prop.value.clone()),
                    "canto_first_approved_date" =>
                        canto_first_approved_date = Some(prop.value.clone()),
                    "canto_approved_date" =>
                        canto_approved_date = Some(prop.value.clone()),
                    "canto_approver_orcid" =>
                        canto_approver_orcid = Some(prop.value.clone()),
                    "canto_added_date" =>
                        canto_added_date = Some(prop.value.clone()),
                    "canto_session_submitted_date" =>
                        canto_session_submitted_date = Some(prop.value.clone()),
                    "annotation_curator" => {
                        let curator = parse_annotation_curator(&prop.value, &prop.prop_type.name);
                        annotation_curators.push(curator)
                    },
                    "annotation_file_curator" => {
                        let curator = parse_annotation_curator(&prop.value, &prop.prop_type.name);
                        annotation_file_curators.push(curator);
                    },
                    _ => ()
                }
            }

            if let Some(ref canto_triage_status) = canto_triage_status {
                let triage_status_to_ignore =
                    &self.config.reference_page_config.triage_status_to_ignore;
                if triage_status_to_ignore.contains(canto_triage_status) {
                    continue;
                }
            }

            let pdb_entries =
                if let Some(ref pdb_ref_entry_map) = self.pdb_ref_entry_map {
                    if let Some(pdbid_entry_map) = pdb_ref_entry_map.get(reference_uniquename.as_ref()) {
                        pdbid_entry_map.values().cloned().collect()
                    } else {
                        vec![]
                    }
                } else {
                    vec![]
                };

            let mut authors_abbrev = None;
            let mut publication_year = None;

            if let Some(authors) = pubmed_authors.clone() {
                if authors.contains(',') {
                    let replaced: String =
                        author_re.replace_all(&authors, "$f et al.").into();
                    authors_abbrev = Some(replaced.to_shared_str());
                } else {
                    authors_abbrev = Some(authors.clone());
                }
            }

            if let Some(publication_date) = pubmed_publication_date.clone() {
                publication_year = Some(publication_date_re.replace_all(&publication_date, "$y").to_shared_str());
            }

            let mut approved_date = canto_first_approved_date.clone();

            if approved_date.is_none() {
                approved_date = canto_session_submitted_date.clone();
            }

            approved_date =
                approved_date.map(|date| approved_date_re.replace_all(&date, "$date").to_shared_str());

            if let Some(ref canto_annotation_status) = canto_annotation_status
                && canto_annotation_status != "APPROVED" {
                    approved_date = None;
                }

            let authors = pubmed_authors.or(non_pubmed_authors);

            for annotation_curator in &annotation_curators {
              if annotation_curator.community_curator {
                canto_curator_name = Some(annotation_curator.name.clone());
                canto_curator_role = flex_str!("community");
              } else if canto_curator_name.is_none() {
                canto_curator_name = Some(annotation_curator.name.clone());
              }
            }

            for annotation_curator in &annotation_file_curators {
              if annotation_curator.community_curator {
                file_curator_name = Some(annotation_curator.name.clone());
                file_curator_role = Some(flex_str!("community"));
              }
            }

            if let Some(first_curator) = annotation_file_curators.first()
              && file_curator_name.is_none() && file_curator_role.is_none() {
                file_curator_name = Some(first_curator.name.clone());
                file_curator_role = Some(self.config.database_name.clone());
              }

            let (pubmed_keyword_genes, extra_genes) =
                self.get_extra_gene_pubs(rc_publication);

            self.references.insert(reference_uniquename.clone(),
                                   ReferenceDetails {
                                       uniquename: reference_uniquename.clone(),
                                       title: rc_publication.title.clone(),
                                       citation: rc_publication.miniref.clone(),
                                       pubmed_abstract: pubmed_abstract.or(non_pubmed_abstract),
                                       pubmed_doi,
                                       authors,
                                       authors_abbrev,
                                       pubmed_publication_date,
                                       pubmed_entrez_date,
                                       canto_session_key,
                                       canto_annotation_status,
                                       canto_triage_status,
                                       canto_curator_role,
                                       canto_curator_name,
                                       canto_first_approved_date,
                                       canto_approved_date,
                                       canto_approver_orcid,
                                       canto_session_submitted_date,
                                       canto_added_date,
                                       annotation_curators,
                                       file_curator_role,
                                       file_curator_name,
                                       annotation_file_curators,
                                       pubmed_keyword_genes,
                                       extra_genes,
                                       approved_date,
                                       publication_year,
                                       cv_annotations: HashMap::new(),
                                       physical_interactions: vec![],
                                       genetic_interactions: HashMap::new(),
                                       ortholog_annotations: vec![],
                                       paralog_annotations: vec![],
                                       genes_by_uniquename: HashMap::new(),
                                       genotypes_by_uniquename: HashMap::new(),
                                       alleles_by_uniquename: HashMap::new(),
                                       references_by_uniquename: HashMap::new(),
                                       transcripts_by_uniquename: HashMap::new(),
                                       terms_by_termid: HashMap::new(),
                                       annotation_details: HashMap::new(),
                                       gene_count: 0,
                                       ltp_gene_count: 0,

                                       pdb_entries,
                                   });
        }
    }

    // make maps from genes to transcript, transcripts to polypeptide,
    // genotypes to genotype interactions, exon, intron, UTRs
    fn make_feature_rel_maps(&mut self) {
        for feature_rel in &self.raw.feature_relationships {
            let subject_type_name = &feature_rel.subject.feat_type.name;
            let rel_name = &feature_rel.rel_type.name;
            let object_type_name = &feature_rel.object.feat_type.name;
            let subject_uniquename = &feature_rel.subject.uniquename;
            let object_uniquename = &feature_rel.object.uniquename;

            if TRANSCRIPT_FEATURE_TYPES.contains(&subject_type_name.as_str()) &&
                rel_name == "part_of" && is_gene_type(object_type_name) {
                    self.genes_of_transcripts.insert(subject_uniquename.clone(),
                                                     object_uniquename.clone());
                    continue;
                }
            if subject_type_name == "polypeptide" &&
                rel_name == "derives_from" &&
                object_type_name == "mRNA" {
                    self.transcripts_of_polypeptides.insert(subject_uniquename.clone(),
                                                            object_uniquename.clone());
                    continue;
                }
            if subject_type_name == "allele" {
                if feature_rel.rel_type.name == "instance_of" &&
                    (object_type_name == "gene" || object_type_name == "pseudogene") {
                        self.genes_of_alleles.insert(subject_uniquename.clone(),
                                                     object_uniquename.clone());
                        continue;
                    }
                if feature_rel.rel_type.name == "part_of" &&
                    object_type_name == "genotype" {
                        let (expression, promoter_gene) =
                            get_feat_rel_expression_and_promoter(&feature_rel.subject, feature_rel);
                        let genotype_locus_identifier =
                            get_feat_rel_prop_value(&flex_str!("genotype_locus"), feature_rel)
                            .unwrap_or_else(|| {
                                flex_fmt!("{}-{}", feature_rel.object.uniquename,
                                          feature_rel.feature_relationship_id)
                            });

                        let allele_uniquename = subject_uniquename;
                        let allele_and_expression =
                            ExpressedAllele {
                                allele_uniquename: allele_uniquename.clone(),
                                expression,
                                promoter_gene,
                            };

                        let genotype_uniquename = object_uniquename;
                        let genotype_entry = self.loci_of_genotypes.entry(genotype_uniquename.clone());
                        let locus_map = genotype_entry.or_default();

                        let genotype_locus =
                            locus_map.entry(genotype_locus_identifier)
                            .or_insert_with(|| GenotypeLocus {
                                expressed_alleles: vec![]
                            });

                        genotype_locus.expressed_alleles.push(allele_and_expression);

                        self.genotypes_of_alleles.entry(allele_uniquename.clone())
                           .or_default()
                           .insert(genotype_uniquename.clone());

                        continue;
                    }
            }
            if TRANSCRIPT_PART_TYPES.contains(&subject_type_name.as_str()) {
                let entry = self.parts_of_transcripts.entry(object_uniquename.clone());
                let part = make_feature_short(&self.chromosomes, &feature_rel.subject);
                entry.or_default().push(part);
            }

            if subject_type_name == "gene" && object_type_name == "protein-containing complex" {
                let entry = self.genes_of_complexes.entry(object_uniquename.clone());
                entry.or_default().insert(subject_uniquename.clone());
            }

            if object_type_name == "genotype_interaction" {
                let genotype_interaction_uniquename = object_uniquename.clone();


                match rel_name.as_str() {
                    "interaction_genotype_a" => {
                        let genotype_a_uniquename = subject_uniquename.clone();
                        self.genotype_interaction_genotype_a.insert(genotype_interaction_uniquename,
                                                                    genotype_a_uniquename);
                    },
                    "interaction_genotype_b" => {
                        let genotype_b_uniquename = subject_uniquename.clone();
                        self.genotype_interaction_genotype_b.insert(genotype_interaction_uniquename,
                                                                    genotype_b_uniquename);
                    },
                    "interaction_double_mutant_genotype" => {
                        let double_mutant_genotype = subject_uniquename.clone();
                        self.genotype_interaction_double_mutant.insert(genotype_interaction_uniquename,
                                                                       double_mutant_genotype);
                    },
                    _ => {
                        panic!("unknown relation type {} for interaction {}", rel_name,
                               object_uniquename);
                    },
                }
            }
        }
    }

    fn get_feature_dbxrefs(&self, feature: &Feature) -> HashSet<FlexStr> {
        if let Some(dbxrefs) = self.dbxrefs_of_features.get(&feature.uniquename) {
            dbxrefs.clone()
        } else {
            HashSet::new()
        }
    }

    fn store_gene_details(&mut self, feat: &Feature) {
        let gene_uniquename = feat.uniquename.clone();

        let maybe_location = make_location(&self.chromosomes, feat);

        if let Some(ref location) = maybe_location
            && let Some(ref mut chr) = self.chromosomes.get_mut(&location.chromosome_name) {
                chr.gene_uniquenames.push(feat.uniquename.clone());
            }

        let organism = make_organism(&feat.organism);
        let dbxrefs = self.get_feature_dbxrefs(feat);

        let mut orfeome_identifier = None;
        for dbxref in &dbxrefs {
            if let Some(without_prefix) = dbxref.strip_prefix("SPD:") {
                orfeome_identifier = Some(without_prefix.to_shared_str());
            }
        }

        let mut uniprot_identifier = None;
        let mut secondary_identifier = None;
        let mut agr_identifier = None;
        let mut biogrid_interactor_id: Option<u32> = None;
        let mut rnacentral_urs_identifier = None;
        let mut rnacentral_2d_structure_id = None;
        let pdb_entries =
            if let Some(ref pdb_entry_map) = self.pdb_gene_entry_map {
                if let Some(pdb_entries) = pdb_entry_map.get(&gene_uniquename) {
                    pdb_entries.clone()
                } else {
                    vec![]
                }
            } else {
                vec![]
            };

        let mut tfexplorer_chipseq_identifier = None;
        let mut tfexplorer_ipms_identifier = None;
        let mut pombephosphoproteomics_unige_ch_starvation_mating_gene: Option<FlexStr> = None;
        let mut pombephosphoproteomics_unige_ch_fusion_gene: Option<FlexStr> = None;

        for prop in feat.featureprops.borrow().iter() {
            match prop.prop_type.name.as_str() {
                "uniprot_identifier" => uniprot_identifier = prop.value.clone(),
                "sgd_identifier" => secondary_identifier = prop.value.clone(),
                "agr_identifier" => agr_identifier = prop.value.clone(),
                "biogrid_interactor_id" => {
                    if let Some(ref chado_biogrid_id) = prop.value {
                        biogrid_interactor_id = match chado_biogrid_id.parse::<u32>() {
                            Ok(val) => Some(val),
                            Err(err) =>
                                panic!("error parsing BioGRID interactor ID from Chado: {}", err),
                        }
                    }
                },
                "pombephosphoproteomics_unige_ch_starvation_mating_gene" =>
                    pombephosphoproteomics_unige_ch_starvation_mating_gene = prop.value.clone(),
                "pombephosphoproteomics_unige_ch_fusion_gene" =>
                    pombephosphoproteomics_unige_ch_fusion_gene = prop.value.clone(),
                "rnacentral_identifier" => rnacentral_urs_identifier = prop.value.clone(),
                "rnacentral_2d_structure_id" => rnacentral_2d_structure_id = prop.value.clone(),
                "tfexplorer_chipseq_identifier" =>
                    tfexplorer_chipseq_identifier = prop.value.clone(),
                "tfexplorer_ipms_identifier" =>
                    tfexplorer_ipms_identifier = prop.value.clone(),
                _ => (),
            }
        }

        let (interpro_matches, tm_domain_coords, coiled_coil_coords,
             disordered_region_coords, low_complexity_region_coords) =
            if let Some(result) = self.domain_data.domains_by_id.get(gene_uniquename.as_str()) {
                let tm_domain_matches = result.tmhmm_matches.iter()
                    .map(|tm_match| AssignedByPeptideRange {
                        range: PeptideRange {
                            start: tm_match.start,
                            end: tm_match.end,
                        },
                        assigned_by: Some(flex_str!("TMHMM")),
                    })
                    .collect::<Vec<_>>();
                let mut coiled_coil_coords = vec![];
                let mut disordered_region_coords = vec![];

                for interpro_match in result.interpro_matches.iter() {
                    let dbname = interpro_match.dbname.as_str();
                    if dbname == "COILS" || dbname == "MOBIDB-Disorder" {
                        for loc in interpro_match.locations.iter() {
                            let range = AssignedByPeptideRange {
                                range: PeptideRange {
                                    start: loc.start,
                                    end: loc.end,
                                },
                                // adding the DB to every range would be excessive
                                assigned_by: None,
                            };

                            if dbname == "COILS" {
                                coiled_coil_coords.push(range);
                            } else {
                                disordered_region_coords.push(range);
                            }
                        }
                    }
                }

                let low_complexity_region_coords = result.segmasker_matches
                    .iter()
                    .map(|loc| {
                        AssignedByPeptideRange {
                            range: PeptideRange {
                                start: loc.start,
                                end: loc.end,
                            },
                            // adding "segmasker" to every struct would be excessive
                            assigned_by: None,
                        }
                    })
                    .collect();

                (result.interpro_matches.clone(), tm_domain_matches,
                 coiled_coil_coords, disordered_region_coords, low_complexity_region_coords)
            } else {
                (vec![], vec![], vec![], vec![], vec![])
            };

        let rfam_annotations =
            if let Some(ref rnacentral_data) = self.rnacentral_data {
                if let Some(ref rnacentral_urs_identifier) = rnacentral_urs_identifier {
                    if let Some(result) = rnacentral_data.get(rnacentral_urs_identifier) {
                        result.clone()
                    } else {
                        vec![]
                    }
                } else {
                    vec![]
                }
            } else {
                vec![]
            };

        let gene_history =
            if let Some(ref all_gene_history) = self.all_gene_history {
                if let Some(gene_history) = all_gene_history.get(&gene_uniquename) {
                    gene_history.clone()
                } else {
                    vec![]
                }
            } else {
                vec![]
            };

        let gene_feature = GeneDetails {
            uniquename: gene_uniquename,
            name: feat.name.clone(),
            taxonid: organism.taxonid,
            product: None,
            deletion_viability: DeletionViability::Unknown,
            uniprot_identifier,
            secondary_identifier,
            agr_identifier,
            biogrid_interactor_id,
            rnacentral_urs_identifier,
            rnacentral_2d_structure_id,
            pdb_entries,
            interpro_matches,
            tm_domain_coords,
            disordered_region_coords,
            low_complexity_region_coords,
            coiled_coil_coords,
            signal_peptide: None,
            transit_peptide: None,
            binding_sites: vec![],
            active_sites: vec![],
            beta_strands: vec![],
            helices: vec![],
            turns: vec![],
            propeptides: vec![],
            chains: vec![],
            glycosylation_sites: vec![],
            disulfide_bonds: vec![],
            lipidation_sites: vec![],
            has_protein_features: false, // is set later
            rfam_annotations,
            orfeome_identifier,
            tfexplorer_chipseq_identifier,
            tfexplorer_ipms_identifier,
            pombephosphoproteomics_unige_ch_starvation_mating_gene,
            pombephosphoproteomics_unige_ch_fusion_gene,
            gocams: HashSet::new(),
            name_descriptions: vec![],
            synonyms: vec![],
            dbxrefs,
            flags: HashSet::new(),
            feature_type: feat.feat_type.name.clone(),
            feature_so_termid: feat.feat_type.termid(),
            transcript_so_termid: None,
            characterisation_status: None,
            taxonomic_distribution: None,
            location: maybe_location,
            gene_neighbourhood: vec![],
            cv_annotations: HashMap::new(),
            physical_interactions: vec![],
            genetic_interactions: HashMap::new(),
            ortholog_annotations: vec![],
            paralog_annotations: vec![],
            target_of_annotations: vec![],
            transcripts: BTreeSet::new(),
            transcripts_by_uniquename: HashMap::new(),
            genes_by_uniquename: HashMap::new(),
            genotypes_by_uniquename: HashMap::new(),
            alleles_by_uniquename: HashMap::new(),
            references_by_uniquename: HashMap::new(),
            terms_by_termid: HashMap::new(),
            annotation_details: HashMap::new(),
            feature_publications: HashSet::new(),
            subset_termids: HashSet::new(),
            split_by_parent_groups: HashMap::new(),

            gene_history,
        };

        self.genes.insert(feat.uniquename.clone(), gene_feature);
    }

    fn get_transcript_parts(&mut self, transcript_uniquename: &FlexStr) -> Vec<FeatureShort> {
        if let Some(mut parts) = self.parts_of_transcripts.remove(transcript_uniquename) {
            if parts.is_empty() {
                panic!("transcript has no parts: {}", transcript_uniquename);
            }

            let strand = parts[0].location.strand;

            let part_cmp = |a: &FeatureShort, b: &FeatureShort| {
                a.location.start_pos.cmp(&b.location.start_pos)
            };

            parts.sort_by(&part_cmp);

            validate_transcript_parts(transcript_uniquename, &parts);

            let chr_name = &parts[0].location.chromosome_name.clone();
            if let Some(chromosome) = self.chromosomes.get(chr_name) {
                let maybe_frameshift_pos =
                    add_introns_to_transcript(chromosome, transcript_uniquename, strand, &mut parts);
                if let Some(frameshift_pos) = maybe_frameshift_pos {
                    let frame_shift_detail = (transcript_uniquename.to_owned(), frameshift_pos);
                    self.transcript_frameshifts_to_check.push(frame_shift_detail);
                }
            } else {
                panic!("can't find chromosome details for: {}", chr_name);
            }

            if parts[0].location.strand == Strand::Reverse {
                parts.reverse();
            }

            parts
        } else {
            vec![]
        }
    }

    fn store_transcript_details(&mut self, feat: &Feature) {
        let transcript_uniquename = feat.uniquename.clone();

        let parts = self.get_transcript_parts(&transcript_uniquename);

        if parts.is_empty() {
            return;
        }

        let mut transcript_start = usize::MAX;
        let mut transcript_end = 0;
        let mut rna_seq_length_spliced = 0;

        for part in &parts {
            if part.location.start_pos < transcript_start {
                transcript_start = part.location.start_pos;
            }
            if part.location.end_pos > transcript_end {
                transcript_end = part.location.end_pos;
            }

            if part.feature_type == FeatureType::Exon {
                rna_seq_length_spliced += part.location.len();
            }
        }

        let rna_seq_length_spliced = NonZeroUsize::new(rna_seq_length_spliced);

        // use the first part as a template to get the chromosome details
        let transcript_location =
            ChromosomeLocation {
                start_pos: transcript_start,
                end_pos: transcript_end,
                phase: None,
                .. parts[0].location.clone()
            };

        let (rna_seq_length_unspliced, maybe_cds_location) =
            if feat.feat_type.name == "mRNA" {
                let mut cds_start = usize::MAX;
                let mut cds_end = 0;

                for part in &parts {
                    if part.feature_type == FeatureType::Exon {
                        if part.location.start_pos < cds_start {
                            cds_start = part.location.start_pos;
                        }
                        if part.location.end_pos > cds_end {
                            cds_end = part.location.end_pos;
                        }
                    }
                }

                if cds_end == 0 {
                    (None, None)
                } else if let Some(mrna_location) = feat.featurelocs.borrow().first() {
                    let first_part_loc = &parts[0].location;

                    (NonZeroUsize::new((cds_end + 1).saturating_sub(cds_start)),
                     Some(ChromosomeLocation {
                          chromosome_name: first_part_loc.chromosome_name.clone(),
                          start_pos: cds_start,
                          end_pos: cds_end,
                          strand: first_part_loc.strand,
                          phase: make_phase(mrna_location),
                      }))
                } else {
                    (None, None)
                }
            } else {
                let rna_length =
                    (transcript_location.end_pos + 1).saturating_sub(transcript_location.start_pos);
                let rna_length = NonZeroUsize::new(rna_length);
                (rna_length, None)
            };

        if let Some(gene_uniquename) =
            self.genes_of_transcripts.get(&transcript_uniquename) {
                let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
                let transcript_type = feat.feat_type.name.clone();
                if gene_details.feature_type == "gene" {
                    let feature_type = format!("{} {}", transcript_type, gene_details.feature_type);
                    gene_details.feature_type = feature_type.to_shared_str();
                }
                let name =
                    if let Some(ref gene_name) = gene_details.name {
                        if let Some(captures) = TRANSCRIPT_ID_RE.captures(transcript_uniquename.as_ref()) {
                            if &captures["gene"] == gene_uniquename.as_ref() {
                                Some(flex_fmt!("{}.{}", gene_name, &captures["suffix"]))
                            } else {
                                panic!("transcript uniquename ({}) doesn't start with the gene uniquename ({})",
                                       transcript_uniquename, gene_uniquename);
                            }
                        } else {
                            panic!("ID doesn't like like a transcript ID: {}", transcript_uniquename);
                        }
                    } else {
                        None
                    };
                let transcript = TranscriptDetails {
                    uniquename: transcript_uniquename.clone(),
                    name,
                    location: transcript_location,
                    transcript_type,
                    parts,
                    protein: None,
                    cds_location: maybe_cds_location,
                    gene_uniquename: gene_uniquename.to_owned(),
                    rna_seq_length_spliced,
                    rna_seq_length_unspliced,
                };

                self.transcripts.insert(transcript_uniquename.clone(), transcript);

                gene_details.transcripts.insert(transcript_uniquename);
                gene_details.transcript_so_termid = Some(feat.feat_type.termid());
            } else {
                panic!("can't find gene for transcript: {}", transcript_uniquename);
            }
    }

    fn store_protein_details(&mut self, feat: &Feature) {
        if let Some(residues) = feat.residues.clone() {
            let protein_uniquename = feat.uniquename.clone();

            let mut molecular_weight = None;
            let mut average_residue_weight = None;
            let mut charge_at_ph7    = None;
            let mut isoelectric_point = None;
            let mut codon_adaptation_index = None;

            let parse_prop_as_f32 = |p: &Option<FlexStr>| {
                if let Some(prop_value) = p {
                    let maybe_value = prop_value.parse();
                    if let Ok(parsed_prop) = maybe_value {
                        Some(parsed_prop)
                    } else {
                        println!("{}: couldn't parse {} as f32",
                                 feat.uniquename, &prop_value);
                        None
                    }
                } else {
                    None
                }
            };

            for prop in feat.featureprops.borrow().iter() {
                if prop.prop_type.name == "molecular_weight"
                    && let Some(value) = parse_prop_as_f32(&prop.value) {
                        molecular_weight = Some(value / 1000.0);
                    }
                if prop.prop_type.name == "average_residue_weight"
                    && let Some(value) = parse_prop_as_f32(&prop.value) {
                        average_residue_weight = Some(value / 1000.0);
                    }
                if prop.prop_type.name == "charge_at_ph7" {
                    charge_at_ph7 = parse_prop_as_f32(&prop.value);
                }
                if prop.prop_type.name == "isoelectric_point" {
                    isoelectric_point = parse_prop_as_f32(&prop.value);
                }
                if prop.prop_type.name == "codon_adaptation_index" {
                    codon_adaptation_index = parse_prop_as_f32(&prop.value);
                }
            }

            if molecular_weight.is_none() {
                panic!("{} has no molecular_weight", feat.uniquename)
            }

            let number_of_residues =
                if residues.ends_with("*") {
                    residues.len() - 1
                } else {
                    residues.len()
                };

            let protein = ProteinDetails {
                uniquename: feat.uniquename.clone(),
                sequence: residues.to_shared_str(),
                number_of_residues,
                product: None,
                molecular_weight: molecular_weight.unwrap(),
                average_residue_weight: average_residue_weight.unwrap(),
                charge_at_ph7: charge_at_ph7.unwrap(),
                isoelectric_point: isoelectric_point.unwrap(),
                codon_adaptation_index: codon_adaptation_index.unwrap(),
            };

            if let Some(transcript_uniquename) =
                self.transcripts_of_polypeptides.get(&protein_uniquename) {
                    self.transcripts.get_mut(transcript_uniquename)
                        .unwrap_or_else(|| panic!("internal error, failed to find transcript: {}",
                                         transcript_uniquename))
                        .protein = Some(protein);
                } else {
                    panic!("can't find transcript of polypeptide: {}", protein_uniquename)
                }
        } else {
            panic!("no residues for protein: {}", feat.uniquename);
        }
    }

    fn store_protein_complex(&mut self, feat: &Feature) {
        let complex_uniquename = feat.uniquename.clone();

        let details = ProteinComplexDetails {
            complex_uniquename: complex_uniquename.clone(),
            complex_name: feat.name.clone(),
            genes: self.genes_of_complexes.get(&complex_uniquename)
              .map_or_else(HashSet::new, HashSet::clone),
        };

        self.protein_complexes.insert(complex_uniquename, details);
    }

    fn store_chromosome_details(&mut self, feat: &Feature) {
        let mut ena_identifier = None;

        for prop in feat.featureprops.borrow().iter() {
            if prop.prop_type.name == "ena_id" {
                ena_identifier = prop.value.clone()
            }
        }

        if feat.residues.is_none() {
            panic!("{:?}", feat.uniquename);
        }

        let org = make_organism(&feat.organism);

        let residues = feat.residues.clone().unwrap();

        if !residues.is_ascii() {
            panic!("sequence for chromosome {} contains non-ascii characters",
                   feat.uniquename);
        }

        let chr = ChromosomeDetails {
            name: feat.uniquename.clone(),
            residues: residues.to_shared_str(),
            ena_identifier: ena_identifier.unwrap().to_shared_str(),
            gene_uniquenames: vec![],
            taxonid: org.taxonid,
            gene_count: 0,    // we'll update the counts once the genes are processed
            coding_gene_count: 0,
        };

        self.chromosomes.insert(feat.uniquename.clone(), chr);
    }

    fn store_genotype_details(&mut self, feat: &Feature) {
        let genotype_uniquename = &feat.uniquename;
        let mut loci: Vec<_> =
            self.loci_of_genotypes[genotype_uniquename]
            .values().cloned().collect();

        // example: "aps1delta__asp1-H397A-H397A-amino_acid_mutation-expression-not_assayed"
        let genotype_display_uniquename =
            make_genotype_display_uniquename(&loci, &self.alleles);

        // example: "aps1delta asp1-H397A(aa)"
        let display_name =
            make_genotype_display_name(&loci, &self.alleles);

        self.genotype_display_uniquenames.insert(genotype_uniquename.clone(),
                                           genotype_display_uniquename.clone());

        let mut ploidiness = Ploidiness::Haploid;
        let mut comment: Option<FlexStr> = None;

        let organism = make_organism(&feat.organism);

        for locus in &loci {
            if locus.expressed_alleles.len() > 1 {
                ploidiness = Ploidiness::Diploid;
                break;
            }
        }

        {
            let loci_cmp = |locus1: &GenotypeLocus, locus2: &GenotypeLocus| {
                let locus1_display_name =
                    &self.alleles[&locus1.expressed_alleles[0].allele_uniquename]
                    .encoded_name_and_type;
                let locus2_display_name =
                    &self.alleles[&locus2.expressed_alleles[0].allele_uniquename]
                    .encoded_name_and_type;
                locus1_display_name.cmp(locus2_display_name)
            };

            loci.sort_by(&loci_cmp);
        }

        for prop in feat.featureprops.borrow().iter() {
            if prop.prop_type.name == "genotype_background" {
                if let Some(ref background) = prop.value {
                    self.genotype_backgrounds.insert(feat.uniquename.clone(),
                                                     background.clone());
                }
            } else if prop.prop_type.name == "genotype_comment"
            && let Some(ref comment_ref) = prop.value {
                comment = Some(comment_ref.to_shared_str());
            }
        }

        let display_uniquename = genotype_display_uniquename.to_shared_str();

        self.genotypes.insert(display_uniquename.clone(),
                              GenotypeDetails {
                                  display_uniquename,
                                  display_name,
                                  name: feat.name.as_ref().map(|s| s.to_shared_str()),
                                  taxonid: organism.taxonid,
                                  loci,
                                  ploidiness,
                                  comment,
                                  cv_annotations: HashMap::new(),
                                  double_mutant_genetic_interactions: HashMap::new(),
                                  rescue_genetic_interactions: HashMap::new(),
                                  genes_by_uniquename: HashMap::new(),
                                  alleles_by_uniquename: HashMap::new(),
                                  references_by_uniquename: HashMap::new(),
                                  transcripts_by_uniquename: HashMap::new(),
                                  genotypes_by_uniquename: HashMap::new(),
                                  terms_by_termid: HashMap::new(),
                                  annotation_details: HashMap::new(),
                                  annotation_count: 0,
                              });
    }

    fn store_allele_details(&mut self, feat: &Feature) {
        let mut allele_type = None;
        let mut description = None;
        let mut comments = vec![];

        for prop in feat.featureprops.borrow().iter() {
            match &prop.prop_type.name as &str {
                "allele_type" =>
                    allele_type = prop.value.clone(),
                "description" =>
                    description = prop.value.clone(),
                "comment" => if let Some(ref comment) = prop.value {
                    let comment = CommentAndReference {
                      comment: comment.clone(),
                      reference:
                         prop.featureprop_pubs.borrow().first()
                             .map(|publication| publication.uniquename.clone()),
                    };
                    comments.push(comment);
                },
                _ => ()
            }
        }

        if let Some(allele_type) = allele_type {
            let gene_uniquename = &self.genes_of_alleles[&feat.uniquename];
            let gene_details = &self.genes[gene_uniquename];
            let allele_details = AlleleDetails::new(&feat.uniquename,
                                                    &feat.name,
                                                    &allele_type,
                                                    &description,
                                                    &comments,
                                                    feat.is_obsolete,
                                                    gene_details.into());
            self.alleles.insert(feat.uniquename.clone(), allele_details);
        } else {
            panic!("no allele_type cvtermprop for {}", &feat.uniquename);
        }
    }

    fn process_chromosome_features(&mut self) {
        // we need to process all chromosomes before other featuers
        for feat in &self.raw.features {
            if feat.feat_type.name == "chromosome" {
                self.store_chromosome_details(feat);
            }
        }
    }

    fn process_features(&mut self) {
        // we need to process all genes before transcripts
        for feat in &self.raw.features {
            if feat.feat_type.name == "gene" || feat.feat_type.name == "pseudogene" {
                self.store_gene_details(feat);
            }
        }

        for feat in &self.raw.features {
            if TRANSCRIPT_FEATURE_TYPES.contains(&feat.feat_type.name.as_str()) {
                self.store_transcript_details(feat)
            }
        }

        for feat in &self.raw.features {
            if feat.feat_type.name == "polypeptide"{
                self.store_protein_details(feat);
            }
        }

        for feat in &self.raw.features {
            if feat.feat_type.name == "protein-containing complex" {
                self.store_protein_complex(feat);
            }

            if !TRANSCRIPT_FEATURE_TYPES.contains(&feat.feat_type.name.as_str()) &&
                !TRANSCRIPT_PART_TYPES.contains(&feat.feat_type.name.as_str()) &&
                !HANDLED_FEATURE_TYPES.contains(&feat.feat_type.name.as_str())
            {
                // for now, ignore features without locations
                if !feat.featurelocs.borrow().is_empty() {
                    let feature_short = make_feature_short(&self.chromosomes, feat);
                    self.other_features.insert(feat.uniquename.clone(), feature_short);
                }
            }
        }
    }

    fn add_interesting_parents(&mut self) {
        let mut interesting_parents_by_termid: HashMap<FlexStr, HashSet<InterestingParent>> =
            HashMap::new();

        for cvtermpath in &self.raw.cvtermpaths {
            let subject_term = &cvtermpath.subject;
            let subject_termid = subject_term.termid();
            let object_term = &cvtermpath.object;
            let object_termid = object_term.termid();

            let rel_termid =
                match cvtermpath.rel_type {
                    Some(ref rel_type) => {
                        rel_type.termid()
                    },
                    None => panic!("no relation type for {} <-> {}\n",
                                   &subject_term.name, &object_term.name)
                };
            let term_details = self.terms.get(&rel_termid);
            let rel_term_name = term_details.unwrap().name.clone();

            if self.is_interesting_parent(&object_termid, &rel_term_name) {
                interesting_parents_by_termid
                    .entry(subject_termid.clone())
                    .or_default()
                    .insert(InterestingParent {
                        termid: object_termid,
                        rel_name: rel_term_name,
                    });
            };
        }

        for (termid, interesting_parents) in interesting_parents_by_termid {
            let term_details = self.terms.get_mut(&termid).unwrap();
            let interesting_parent_ids = interesting_parents.iter()
                .map(|p| p.termid.clone())
                .collect::<HashSet<_>>();
            term_details.interesting_parent_ids = interesting_parent_ids;
            term_details.interesting_parent_details = interesting_parents;
        }
    }

    fn process_allele_features(&mut self) {
        for feat in &self.raw.features {
            if feat.feat_type.name == "allele" {
                self.store_allele_details(feat);
            }
        }
    }

    fn process_genotype_features(&mut self) {
        for feat in &self.raw.features {
            if feat.feat_type.name == "genotype" {
                self.store_genotype_details(feat);
            }
        }
    }

    fn add_genotypes_to_allele_details(&mut self) {
        let genotype_cmp = |geno1: &GenotypeShort, geno2: &GenotypeShort| {
            let geno1_ploidiness = geno1.ploidiness();
            let geno2_ploidiness = geno2.ploidiness();

            if geno1_ploidiness == Ploidiness::Haploid &&
               geno2_ploidiness == Ploidiness::Diploid {
                return Ordering::Less;
            }

            if geno1_ploidiness == Ploidiness::Diploid &&
               geno2_ploidiness == Ploidiness::Haploid {
                return Ordering::Greater;
            }

            geno1.display_uniquename.to_ascii_lowercase()
                 .cmp(&geno2.display_uniquename.to_ascii_lowercase())
        };
        let phenotype_cmp = |pheno1: &TermShort, pheno2: &TermShort| {
            pheno1.name.to_ascii_lowercase().cmp(&pheno2.name.to_ascii_lowercase())
        };

        for allele_details in self.alleles.values_mut() {
            let mut genotypes: HashSet<GenotypeShort> = HashSet::new();
            let mut phenotypes: HashSet<TermShort> = HashSet::new();

            if let Some(genotype_uniquenames) = self.genotypes_of_alleles.get(&allele_details.uniquename) {

                for genotype_uniquename in genotype_uniquenames {

                    let genotype_display_uniquename =
                        self.genotype_display_uniquenames.get(genotype_uniquename).unwrap();

                    if let Some(genotype_details) = self.genotypes.get(genotype_display_uniquename) {
                        for term_annotations in genotype_details.cv_annotations.values() {
                            for term_annotation in term_annotations {
                                let termid = &term_annotation.term;
                                let term_details =
                                     self.terms.get(termid).unwrap_or_else(|| panic!("missing termid {}", termid));
                                let term_short: TermShort = term_details.into();
                                phenotypes.insert(term_short);
                            }
                        }
                        if genotype_details.annotation_count > 0 {
                            genotypes.insert(genotype_details.into());
                        }
                    } else {
                        panic!("can't find GenotypeDetails for {}", genotype_display_uniquename);
                    }
                }
            }

            allele_details.genotypes.extend(genotypes.drain());
            allele_details.phenotypes.extend(phenotypes.drain());

            allele_details.genotypes.sort_by(genotype_cmp);
            allele_details.phenotypes.sort_by(phenotype_cmp);
        }
    }

    fn add_gene_neighbourhoods(&mut self) {
        struct GeneAndLoc {
            gene_uniquename: FlexStr,
            loc: ChromosomeLocation,
        }

        let mut genes_and_locs: Vec<GeneAndLoc> = vec![];

        for gene_details in self.genes.values() {
            // mRNA or pseudogenic_transcript
            if let Some(ref transcript_so_termid) = gene_details.transcript_so_termid
                && (transcript_so_termid == "SO:0000234" ||
                    transcript_so_termid == "SO:0000516")
                    && let Some(ref location) = gene_details.location {
                        genes_and_locs.push(GeneAndLoc {
                            gene_uniquename: gene_details.uniquename.clone(),
                            loc: location.clone(),
                        });
                    }
        }

        let cmp = |a: &GeneAndLoc, b: &GeneAndLoc| {
            let order = a.loc.chromosome_name.cmp(&b.loc.chromosome_name);
            if order == Ordering::Equal {
                a.loc.start_pos.cmp(&b.loc.start_pos)
            } else {
                order
            }
        };

        genes_and_locs.sort_by(cmp);

        for (i, this_gene_and_loc) in genes_and_locs.iter().enumerate() {
            let mut nearby_genes: Vec<GeneShort> = vec![];
            if i > 0 {
                let start_index =
                    i.saturating_sub(GENE_NEIGHBOURHOOD_DISTANCE);

                for back_index in (start_index..i).rev() {
                    let back_gene_and_loc = &genes_and_locs[back_index];

                    if back_gene_and_loc.loc.chromosome_name !=
                        this_gene_and_loc.loc.chromosome_name {
                            break;
                        }
                    let back_gene_short = self.make_gene_short(&back_gene_and_loc.gene_uniquename);
                    nearby_genes.insert(0, back_gene_short);
                }
            }

            let gene_short = self.make_gene_short(&this_gene_and_loc.gene_uniquename);
            nearby_genes.push(gene_short);

            if i < genes_and_locs.len() - 1 {
                let end_index =
                    if i + GENE_NEIGHBOURHOOD_DISTANCE >= genes_and_locs.len() {
                        genes_and_locs.len()
                    } else {
                        i + GENE_NEIGHBOURHOOD_DISTANCE + 1
                    };

                #[allow(clippy::needless_range_loop)]
                for forward_index in i+1..end_index {
                    let forward_gene_and_loc = &genes_and_locs[forward_index];

                    if forward_gene_and_loc.loc.chromosome_name !=
                        this_gene_and_loc.loc.chromosome_name {
                            break;
                        }

                    let forward_gene_short = self.make_gene_short(&forward_gene_and_loc.gene_uniquename);
                    nearby_genes.push(forward_gene_short);
                }
            }

            let this_gene_details =
                self.genes.get_mut(&this_gene_and_loc.gene_uniquename).unwrap();

            this_gene_details.gene_neighbourhood.append(&mut nearby_genes);
        }
    }

    fn save_genetic_interaction(&mut self, interaction: InteractionAnnotation) {
        let interaction_key =
            GeneticInteractionKey {
                gene_a_uniquename: interaction.gene_uniquename.clone(),
                interaction_type: interaction.evidence.clone(),
                gene_b_uniquename: interaction.interactor_uniquename.clone(),
            };

        let interaction_annotation =
            GeneticInteractionDetail {
                genotype_a_uniquename: None,
                genotype_b_uniquename: None,
                reference_uniquename: interaction.reference_uniquename,
                double_mutant_phenotype_termid: None,
                double_mutant_extension: vec![],
                double_mutant_genotype_display_uniquename: None,
                rescued_phenotype_termid: None,
                rescued_phenotype_extension: vec![],
                interaction_note: interaction.interaction_note,
                throughput: interaction.throughput,
                source_database: interaction.source_database,
                annotation_date: interaction.annotation_date,
            };

        self.genetic_interaction_annotations
            .entry(interaction_key)
            .or_default()
            .push(interaction_annotation);
    }

    // add physical interaction, legacy genetic interaction, ortholog, paralog annotations
    // and GO-CAM genes
    fn process_feature_rels(&mut self) {
        for feature_rel in &self.raw.feature_relationships {
            let rel_name = &feature_rel.rel_type.name;
            let subject_uniquename = &feature_rel.subject.uniquename;
            let object_uniquename = &feature_rel.object.uniquename;

            for rel_config in &FEATURE_REL_CONFIGS {
                if rel_name == rel_config.rel_type_name &&
                    is_gene_type(&feature_rel.subject.feat_type.name) &&
                    is_gene_type(&feature_rel.object.feat_type.name) {
                        let mut evidence: Option<Evidence> = None;
                        let mut throughput: Option<Throughput> = None;
                        let mut is_inferred_interaction: bool = false;
                        let mut interaction_note: Option<FlexStr> = None;
                        let mut source_database = None;
                        let mut ortholog_qualifier = None;
                        let mut annotation_date = None;

                        let borrowed_publications = feature_rel.publications.borrow();
                        let maybe_publication = borrowed_publications.first();
                        let maybe_reference_uniquename =
                            match maybe_publication {
                                Some(publication) =>
                                    if publication.uniquename == "null" {
                                        None
                                    } else {
                                        Some(publication.uniquename.clone())
                                    },
                                None => None,
                            };

                        for prop in feature_rel.feature_relationshipprops.borrow().iter() {
                            if prop.prop_type.name == "evidence"
                                && let Some(ref evidence_long) = prop.value {
                                    for (evidence_code, ev_details) in &self.config.evidence_types {
                                        if &ev_details.long == evidence_long {
                                            evidence = Some(evidence_code.clone());
                                        }
                                    }
                                    if evidence.is_none() {
                                        evidence = Some(evidence_long.clone());
                                    }
                                }
                            if prop.prop_type.name == "is_inferred"
                                && let Some(is_inferred_value) = prop.value.clone()
                                    && is_inferred_value == "yes" {
                                        is_inferred_interaction = true;
                                    }
                            if prop.prop_type.name == "annotation_throughput_type"
                                && let Some(throughput_type) = prop.value.clone() {
                                    throughput = Some(match throughput_type.as_ref() {
                                        "low throughput" => Throughput::LowThroughput,
                                        "high throughput" => Throughput::HighThroughput,
                                        "non-experimental" => Throughput::NonExperimental,
                                        _ => {
                                            panic!("unknown throughput type: {}",
                                                   throughput_type);
                                        }
                                    });
                                }
                            if prop.prop_type.name == "interaction_note"
                                && let Some(interaction_note_value) = prop.value.clone() {
                                    interaction_note = Some(interaction_note_value);
                                }
                            if prop.prop_type.name == "ortholog qualifier" ||
                                prop.prop_type.name == "ortholog_qualifier" {
                                ortholog_qualifier = prop.value.clone()
                            }
                            if prop.prop_type.name == "source_database"
                                && let Some(source_database_value) = prop.value.clone() {
                                    source_database = Some(source_database_value);
                                }
                            if prop.prop_type.name == "date"
                                && let Some(date_value) = prop.value.clone() {
                                    annotation_date = Some(date_value);
                                }
                        }

                        let evidence_clone = evidence.clone();

                        let gene_uniquename = subject_uniquename;
                        let gene_organism_taxonid = {
                            self.genes[subject_uniquename].taxonid
                        };
                        let other_gene_uniquename = object_uniquename;
                        let other_gene_organism_taxonid = {
                            self.genes[object_uniquename].taxonid
                        };
                        match rel_config.annotation_type {
                            FeatureRelAnnotationType::Interaction =>
                                if !is_inferred_interaction {
                                    let evidence =
                                        evidence.unwrap_or_else(|| panic!("evidence missing for feature_relationship_id: {}",
                                                                 feature_rel.feature_relationship_id));
                                    let interaction_annotation =
                                        InteractionAnnotation {
                                            gene_uniquename: gene_uniquename.clone(),
                                            interactor_uniquename: other_gene_uniquename.clone(),
                                            evidence,
                                            reference_uniquename: maybe_reference_uniquename.clone(),
                                            throughput,
                                            interaction_note,
                                            source_database,
                                            annotation_date,
                                        };

                                    if rel_name == "interacts_genetically" {
                                        self.save_genetic_interaction(interaction_annotation);
                                    } else {
                                    {
                                        let gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                                        if rel_name == "interacts_physically" {
                                            gene_details.physical_interactions.push(interaction_annotation.clone());
                                        } else {
                                            panic!("unknown interaction type: {}", rel_name);
                                        };
                                    }
                                    if gene_uniquename != other_gene_uniquename {
                                        let other_gene_details = self.genes.get_mut(object_uniquename).unwrap();
                                        if rel_name == "interacts_physically" {
                                            other_gene_details.physical_interactions.push(interaction_annotation.clone());
                                        } else {
                                            panic!("unknown interaction type: {}", rel_name);
                                        };
                                    }

                                    if let Some(ref_details) =
                                        if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                            self.references.get_mut(reference_uniquename)
                                        } else {
                                            None
                                        }
                                    {
                                        if rel_name == "interacts_physically" {
                                            ref_details.physical_interactions.push(interaction_annotation.clone());
                                        } else {
                                            panic!("unknown interaction type: {}", rel_name);
                                        };
                                    }
                                    self.physical_interaction_annotations.insert(interaction_annotation);
                                    }
                                },
                            FeatureRelAnnotationType::Ortholog => {
                                let ortholog_annotation =
                                    OrthologAnnotation {
                                        gene_uniquename: gene_uniquename.clone(),
                                        ortholog_uniquename: other_gene_uniquename.clone(),
                                        ortholog_taxonid: other_gene_organism_taxonid,
                                        evidence,
                                        reference_uniquename: maybe_reference_uniquename.clone(),
                                        qualifier: ortholog_qualifier.clone(),
                                    };
                                if Some(gene_organism_taxonid) == self.config.load_organism_taxonid
                                    && let Some(ref_details) =
                                        if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                            self.references.get_mut(reference_uniquename)
                                        } else {
                                            None
                                        }
                                    {
                                        ref_details.ortholog_annotations.push(ortholog_annotation.clone());
                                    }
                                let gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                                gene_details.ortholog_annotations.push(ortholog_annotation);
                            },
                            FeatureRelAnnotationType::Paralog => {
                                let paralog_annotation =
                                    ParalogAnnotation {
                                        gene_uniquename: gene_uniquename.clone(),
                                        paralog_uniquename: other_gene_uniquename.clone(),
                                        evidence,
                                        reference_uniquename: maybe_reference_uniquename.clone(),
                                    };
                                let gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                                gene_details.paralog_annotations.push(paralog_annotation.clone());
                                if let Some(ref_details) =
                                    if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                        self.references.get_mut(reference_uniquename)
                                    } else {
                                        None
                                    }
                                    && (self.config.load_organism_taxonid.is_some() &&
                                        self.config.load_organism_taxonid.unwrap() == gene_details.taxonid ||
                                        gene_organism_taxonid < other_gene_organism_taxonid)
                                    {
                                        ref_details.paralog_annotations.push(paralog_annotation);
                                    }
                            }
                        }

                        // for orthologs and paralogs, store the reverse annotation too
                        let other_gene_details = self.genes.get_mut(object_uniquename).unwrap();
                        match rel_config.annotation_type {
                            FeatureRelAnnotationType::Interaction => {},
                            FeatureRelAnnotationType::Ortholog => {
                                let ortholog_annotation =
                                    OrthologAnnotation {
                                        gene_uniquename: other_gene_uniquename.clone(),
                                        ortholog_uniquename: gene_uniquename.clone(),
                                        ortholog_taxonid: gene_organism_taxonid,
                                        evidence: evidence_clone,
                                        reference_uniquename: maybe_reference_uniquename.clone(),
                                        qualifier: ortholog_qualifier,
                                    };
                                if Some(other_gene_organism_taxonid) == self.config.load_organism_taxonid
                                    && let Some(ref_details) =
                                        if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                            self.references.get_mut(reference_uniquename)
                                        } else {
                                            None
                                        }
                                    {
                                        ref_details.ortholog_annotations.push(ortholog_annotation.clone());
                                    }
                                other_gene_details.ortholog_annotations.push(ortholog_annotation);
                            },
                            FeatureRelAnnotationType::Paralog => {
                                let paralog_annotation =
                                    ParalogAnnotation {
                                        gene_uniquename: other_gene_uniquename.clone(),
                                        paralog_uniquename: gene_uniquename.clone(),
                                        evidence: evidence_clone,
                                        reference_uniquename: maybe_reference_uniquename.clone(),
                                    };
                                other_gene_details.paralog_annotations.push(paralog_annotation.clone());
                                if let Some(ref_details) =
                                    if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                        self.references.get_mut(reference_uniquename)
                                    } else {
                                        None
                                    }
                                    && (self.config.load_organism_taxonid.is_some() &&
                                        self.config.load_organism_taxonid.unwrap() == other_gene_details.taxonid ||
                                        gene_organism_taxonid > other_gene_organism_taxonid)
                                    {
                                        ref_details.paralog_annotations.push(paralog_annotation);
                                    }
                            },
                        }
                    }
            }
        }

        for ref_details in self.references.values_mut() {
            ref_details.physical_interactions.sort();
            ref_details.ortholog_annotations.sort();
            ref_details.paralog_annotations.sort();
        }

        for gene_details in self.genes.values_mut() {
            gene_details.physical_interactions.sort();
            gene_details.ortholog_annotations.sort();
            gene_details.paralog_annotations.sort();
        }
    }

    // find the extension_display_names config for the given termid and relation type name
    fn matching_ext_config(&self, annotation_termid: &FlexStr,
                           rel_type_name: &FlexStr) -> Option<ExtensionDisplayNames> {
        let ext_configs = &self.config.extension_display_names;

        if let Some(annotation_term_details) = self.terms.get(annotation_termid) {
            for ext_config in ext_configs {
                if ext_config.rel_name == *rel_type_name {
                    if let Some(ref if_descendant_of) = ext_config.if_descendant_of {
                        if annotation_termid == if_descendant_of.as_str() ||
                            annotation_term_details.interesting_parent_ids.contains(if_descendant_of) {
                            return Some((*ext_config).clone());
                        }
                    } else {
                        return Some((*ext_config).clone());
                    }
                }
            }
        }

        None
    }

    // create and returns any TargetOfAnnotations implied by the extension
    fn make_target_of_for_ext(&self, cv_name: &FlexStr,
                              annotation_termid: &FlexStr,
                              annotation: &OntAnnotationDetail)
        -> Vec<(GeneUniquename, TargetOfAnnotation)>
    {
        let genes = &annotation.genes;
        if genes.len() != 1 {
            panic!("expected an annotation with one gene for {}, got: {:?}",
                   annotation_termid, genes);
        }
        let gene = &genes[0];

        let genotype_uniquename = &annotation.genotype;
        let reference_uniquename = &annotation.reference;
        let assigned_by = &annotation.assigned_by;
        let evidence = &annotation.evidence;
        let extension = &annotation.extension;
        let date = &annotation.date;

        let mut ret_vec = vec![];

        for ext_part in extension {
            let maybe_ext_config =
                self.matching_ext_config(annotation_termid, &ext_part.rel_type_name);
            match ext_part.ext_range {
                ExtRange::Gene(ref target_gene_uniquename) |
                ExtRange::GeneAndGeneProduct(GeneAndGeneProduct {
                     gene_uniquename: ref target_gene_uniquename, product: _
                }) =>
                if let Some(ext_config) = maybe_ext_config
                    && let Some(reciprocal_display_name) =
                        ext_config.reciprocal_display {
                            let (annotation_gene_uniquename, annotation_genotype_uniquename) =
                                if genotype_uniquename.is_some() {
                                    (gene.clone(), genotype_uniquename.clone())
                                } else {
                                    (gene.clone(), None)
                                };
                            ret_vec.push(((*target_gene_uniquename).clone(),
                                          TargetOfAnnotation {
                                              show_in_summary: true,  // set this later
                                              ontology_name: cv_name.into(),
                                              ext_rel_display_name: reciprocal_display_name,
                                              gene: annotation_gene_uniquename,
                                              genotype_uniquename: annotation_genotype_uniquename,
                                              reference_uniquename: reference_uniquename.clone(),
                                              assigned_by: assigned_by.clone(),
                                              date: date.clone(),
                                              evidence: evidence.clone(),
                                          }));
                        },
                _ => (),
            }
        }

        ret_vec
    }

    // return an ordered vector of annotations, setting the show_in_summary flag
    // see: https://github.com/pombase/website/issues/299
    fn process_target_of_annotations(&self, gene_details: &GeneDetails,
                                     annotations: &mut HashSet<TargetOfAnnotation>)
                                     -> Vec<TargetOfAnnotation>
    {
        let mut processed_annotations = annotations.drain().collect::<Vec<_>>();

        let target_of_config = &self.config.target_of_config;
        let priority_config = &target_of_config.relation_priority;

        for annotation in &processed_annotations {
            if priority_config.get(annotation.ext_rel_display_name.as_str()).is_none() {
                eprintln!(r#"No priority configured for "{}" (from {})"#,
                          annotation.ext_rel_display_name, gene_details.uniquename);
            }
        }

        let cmp_fn = |a: &TargetOfAnnotation, b: &TargetOfAnnotation| {
            let a_rel_name = a.ext_rel_display_name.as_str();
            let a_pri = priority_config.get(a_rel_name).unwrap_or(&0);
            let b_rel_name = b.ext_rel_display_name.as_str();
            let b_pri = priority_config.get(b_rel_name).unwrap_or(&0);

            let pri_order = b_pri.cmp(a_pri);

            if pri_order == Ordering::Equal {
                let rel_name_order = a_rel_name.cmp(b_rel_name);
                if rel_name_order == Ordering::Equal {
                    let a_gene_details = self.genes.get(&a.gene).unwrap();
                    let b_gene_details = self.genes.get(&b.gene).unwrap();

                    if let (Some(a_name), Some(b_name)) =
                        (&a_gene_details.name, &b_gene_details.name)
                    {
                        a_name.cmp(b_name)
                    } else if a_gene_details.name.is_some() {
                        Ordering::Less
                    } else if b_gene_details.name.is_some() {
                        Ordering::Greater
                    } else {
                        a_gene_details.uniquename.cmp(&b_gene_details.uniquename)
                    }
                } else {
                    rel_name_order
                }
            } else {
                pri_order
            }
        };

        processed_annotations.sort_by(cmp_fn);

        let mut seen_gene_rels = HashMap::new();

        for annotation in processed_annotations.iter_mut() {
            let rel_priority = priority_config.get(annotation.ext_rel_display_name.as_str())
                .unwrap_or(&0);

            let key = flex_fmt!("{}-{}", annotation.ontology_name, annotation.gene);
            let existing_rel = seen_gene_rels.get(&key);

            if let Some(existing_rel) = existing_rel
                && *existing_rel > *rel_priority {
                    annotation.show_in_summary = false;
                    continue;
                }
            seen_gene_rels.insert(key, *rel_priority);
        }

        processed_annotations
    }

    fn add_target_of_annotations(&mut self) {
        let mut target_of_annotations: HashMap<GeneUniquename, HashSet<TargetOfAnnotation>> =
            HashMap::new();

        for term_details in self.terms.values() {
            for term_annotations in term_details.cv_annotations.values() {
                for term_annotation in term_annotations {
                    'ANNOTATION: for annotation_id in &term_annotation.annotations {
                        let annotation = self.annotation_details
                            .get(annotation_id).expect("can't find OntAnnotationDetail");
                        if let Some(ref genotype_uniquename) = annotation.genotype {
                            let genotype = &self.genotypes[genotype_uniquename];

                            if genotype.loci.len() > 1 ||
                                genotype.loci[0].expressed_alleles.len() > 1 {
                                break 'ANNOTATION;
                            }
                        }

                        let new_annotations =
                            self.make_target_of_for_ext(&term_details.cv_name,
                                                        &term_details.termid,
                                                        annotation);
                        for (target_gene_uniquename, new_annotation) in new_annotations {
                           if self.genes.contains_key(&target_gene_uniquename) {
                               if target_gene_uniquename != new_annotation.gene {
                                   target_of_annotations
                                       .entry(target_gene_uniquename.clone())
                                       .or_default()
                                       .insert(new_annotation);
                               }
                           } else {
                               eprintln!("can't find gene {} in extension for {}",
                                         target_gene_uniquename, term_details.termid);
                               for annotation_gene in &annotation.genes {
                                   eprintln!("  in annotation of {}", annotation_gene);
                               }
                           }
                        }
                    }
                }
            }
        }

        for (gene_uniquename, mut target_of_annotations) in target_of_annotations {
            let gene_details = self.genes.get(&gene_uniquename).unwrap();
            let processed_target_of_annotations =
                self.process_target_of_annotations(gene_details, &mut target_of_annotations);
            let gene_details = self.genes.get_mut(&gene_uniquename).unwrap();
            gene_details.target_of_annotations = processed_target_of_annotations;
        }
    }

    fn set_deletion_viability(&mut self) {
        let some_null = Some(flex_str!("Null"));

        let mut gene_statuses = HashMap::new();

        let condition_string =
            |condition_ids: HashSet<FlexStr>| {
                let mut ids_vec: Vec<FlexStr> = condition_ids.iter().cloned().collect();
                ids_vec.sort();
                join(&ids_vec, " ")
            };

        let viable_termid = &self.config.viability_terms.viable;
        let inviable_termid = &self.config.viability_terms.inviable;

        for (gene_uniquename, gene_details) in &mut self.genes {
            let mut new_status = DeletionViability::Unknown;

            if let Some(single_locus_term_annotations) =
                gene_details.cv_annotations.get(&flex_str!("single_locus_phenotype")) {
                    let mut viable_conditions: HashMap<FlexStr, TermId> = HashMap::new();
                    let mut inviable_conditions: HashMap<FlexStr, TermId> = HashMap::new();

                    for term_annotation in single_locus_term_annotations {
                        'ANNOTATION: for annotation_id in &term_annotation.annotations {
                            let annotation = self.annotation_details
                                .get(annotation_id).expect("can't find OntAnnotationDetail");

                            let genotype_uniquename = annotation.genotype.as_ref().unwrap();

                            let genotype = &self.genotypes[genotype_uniquename];
                            if genotype.loci[0].expressed_alleles.len() > 1 {
                                // diploid locus
                                continue 'ANNOTATION;
                            }
                            let expressed_allele = &genotype.loci[0].expressed_alleles[0];
                            let allele = &self.alleles[&expressed_allele.allele_uniquename];
                            if allele.allele_type != "deletion" &&
                                expressed_allele.expression != some_null {
                                continue 'ANNOTATION;
                            }

                            let term = &self.terms[&term_annotation.term];
                            let interesting_parent_ids = &term.interesting_parent_ids;
                            let conditions_as_string =
                                condition_string(annotation.conditions.clone());
                            if interesting_parent_ids.contains(viable_termid) ||
                                *viable_termid == term_annotation.term {
                                    viable_conditions.insert(conditions_as_string,
                                                             term_annotation.term.clone());
                                } else if interesting_parent_ids.contains(inviable_termid) ||
                                *inviable_termid == term_annotation.term {
                                    inviable_conditions.insert(conditions_as_string,
                                                               term_annotation.term.clone());
                                }
                        }
                    }

                    if viable_conditions.is_empty() {
                        if !inviable_conditions.is_empty() {
                            new_status = DeletionViability::Inviable;
                        }
                    } else if inviable_conditions.is_empty() {
                        new_status = DeletionViability::Viable;
                    } else {
                        new_status = DeletionViability::DependsOnConditions;

                        let viable_conditions_set: HashSet<FlexStr> =
                            viable_conditions.keys().cloned().collect();
                        let inviable_conditions_set: HashSet<FlexStr> =
                            inviable_conditions.keys().cloned().collect();

                        let intersecting_conditions =
                            viable_conditions_set.intersection(&inviable_conditions_set);
                        if intersecting_conditions.clone().count() > 0 {
                            println!("{} is viable and inviable with", gene_uniquename);
                            for cond in intersecting_conditions {
                                if cond.is_empty() {
                                    println!("  no conditions");
                                } else {
                                    println!("  conditions: {}", cond);
                                }
                                println!("   viable term: {}",
                                         viable_conditions[cond]);
                                println!("   inviable term: {}",
                                         inviable_conditions[cond]);
                            }
                        }
                    }
                }
            gene_statuses.insert(gene_uniquename.clone(), new_status);
        }

        for (gene_uniquename, status) in &gene_statuses {
            if let Some(ref mut gene_details) = self.genes.get_mut(gene_uniquename) {
                gene_details.deletion_viability = status.clone();
            }
        }
    }

    fn set_gene_flags(&mut self) {
        let mut nucleosome_genes = HashSet::new();

        for gene_details in self.genes.values() {
            if let Some(cc_annotations) = gene_details.cv_annotations.get("cellular_component") {
                for cc_term_annotation in cc_annotations {
                    let termid = &cc_term_annotation.term;
                    let term_details = self.terms.get(termid).expect("term missing from terms map");

                    if termid == "GO:0000786" ||
                        term_details.interesting_parent_ids.contains("GO:0000786")
                    {
                        nucleosome_genes.insert(gene_details.uniquename.clone());
                    }
                }
            }
        }

        for gene_details in self.genes.values_mut() {
            if nucleosome_genes.contains(&gene_details.uniquename) {
                gene_details.flags.insert(flex_str!("is_histone"));
            }
            if let Some(load_organism_taxonid) = self.config.load_organism_taxonid
                && gene_details.taxonid != load_organism_taxonid {
                    gene_details.flags.insert(flex_str!("not_load_organism"));
                }
        }
    }

    fn set_term_details_subsets(&mut self) {
        let mut subsets_by_termid = HashMap::new();
        for (slim_name, slim_config) in self.config.slims.iter() {
            for term_and_name in &slim_config.terms {
                subsets_by_termid
                    .entry(term_and_name.termid.clone())
                    .or_insert_with(HashSet::new)
                    .insert(slim_name.clone());
            }
        }


        for term_details in self.terms.values_mut() {
            if let Some(subsets) = subsets_by_termid.remove(&term_details.termid) {
                term_details.in_subsets = subsets;
            }
        }
    }

    // On each GeneDetails, add a set of the term IDs of subsets for
    // this gene.  Any useful subset that contains any term for any
    // annotation in the gene is included.  "useful" means that the
    // front end might need it, eg. slim term IDs
    fn set_gene_details_subset_termids(&mut self) {
        let is_subset_member =
            |subset_termid: &FlexStr, test_termid: &FlexStr| {
                if subset_termid == test_termid {
                    return true;
                }
                if let Some(children) = self.children_by_termid.get(subset_termid) {
                    children.contains(test_termid)
                } else {
                    false
                }
            };

        let mut subsets_by_gene = HashMap::new();
        for slim_config in self.config.slims.values() {
            for term_and_name in &slim_config.terms {
                for gene_details in self.genes.values() {
                    for term_annotations in gene_details.cv_annotations.values() {
                        for term_annotation in term_annotations {
                            if term_annotation.is_not {
                                // NOT annotations aren't included in subsets
                                continue;
                            }
                            let gene_termid = &term_annotation.term;
                            if is_subset_member(&term_and_name.termid, gene_termid) {
                                subsets_by_gene
                                    .entry(gene_details.uniquename.clone())
                                    .or_insert_with(HashSet::new)
                                    .insert(term_and_name.termid.clone());
                            }
                        }
                    }
                }
            }
        }

        for gene_details in self.genes.values_mut() {
            if let Some(subset_termids) = subsets_by_gene.remove(&gene_details.uniquename) {
                gene_details.subset_termids = subset_termids;
            }
        }
    }

    fn set_taxonomic_distributions(&mut self) {
        let mut term_name_map = HashMap::new();

        let in_archaea = "conserved in archaea".to_shared_str();
        let in_bacteria = "conserved in bacteria".to_shared_str();
        let in_fungi_only = "conserved in fungi only".to_shared_str();
        let in_metazoa = "conserved in metazoa".to_shared_str();
        let pombe_specific = "Schizosaccharomyces pombe specific".to_shared_str();
        let schizo_specific = "Schizosaccharomyces specific".to_shared_str();

        let names = vec![in_archaea.clone(), in_bacteria.clone(),
                         in_fungi_only.clone(), in_metazoa.clone(),
                         pombe_specific.clone(), schizo_specific.clone()];

        for name in names {
            if let Some(termid) = self.term_ids_by_name.get(&name) {
                term_name_map.insert(termid.clone(), name.to_owned());
            } else {
                eprintln!("configuration error: can't find {} in term_ids_by_name map", name);
                eprintln!("skipping taxonomic distribution");
                return;
            }
        }

        'GENE:
        for gene_details in self.genes.values_mut() {
            let mut dist_names = HashSet::new();

            if let Some(species_dists) = gene_details.cv_annotations.get(&flex_str!("species_dist")) {
                for ont_term_annotations in species_dists {
                    let term = &ont_term_annotations.term;
                    if let Some(term_name) = term_name_map.get(term) {
                        dist_names.insert(term_name.to_owned());
                    }
                }
            }

            if (dist_names.contains(&in_archaea) || dist_names.contains(&in_bacteria))
                && !dist_names.contains(&in_metazoa) {
                    gene_details.taxonomic_distribution =
                        Some(flex_str!("fungi and prokaryotes"));
                    continue 'GENE;
                }
            if dist_names.contains(&in_metazoa) &&
                !((dist_names.contains(&in_archaea) || dist_names.contains(&in_bacteria))
                  && dist_names.contains(&in_metazoa)) {
                    gene_details.taxonomic_distribution =
                        Some(flex_str!("eukaryotes only, fungi and metazoa"));
                    continue 'GENE;
                }


            if (dist_names.contains(&in_archaea) || dist_names.contains(&in_bacteria)) &&
                dist_names.contains(&in_metazoa) {
                    gene_details.taxonomic_distribution =
                        Some(flex_str!("eukaryotes and prokaryotes"));
                    continue 'GENE;
                }

            if dist_names.contains(&in_fungi_only) {
                gene_details.taxonomic_distribution = Some(flex_str!("fungi only"));
                continue 'GENE;
            }

            if dist_names.contains(&pombe_specific) {
                gene_details.taxonomic_distribution = Some(flex_str!("S. pombe specific"));
                continue 'GENE;
            }

            if dist_names.contains(&schizo_specific) {
                gene_details.taxonomic_distribution = Some(flex_str!("Schizos. specific"));
                continue 'GENE;
            }

            if let Some(ref characterisation_status) = gene_details.characterisation_status
                && characterisation_status == "dubious" {
                    gene_details.taxonomic_distribution = Some(flex_str!("dubious"));
                    continue 'GENE;
                }

            if gene_details.feature_type != "mRNA gene" {
                gene_details.taxonomic_distribution = Some(flex_str!("not curated"));
                continue 'GENE;
            }

            gene_details.taxonomic_distribution = Some(flex_str!("other"));
        }
    }

    fn process_cvterms(&mut self) {
        let mut pro_term_to_gene = HashMap::new();

        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name != POMBASE_ANN_EXT_TERM_CV_NAME {
                let cv_config = self.config.cv_config_by_name_with_default(&cvterm.cv.name);
                let annotation_feature_type = cv_config.feature_type.clone();

                let mut pombase_gene_id = None;
                for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                    if cvtermprop.prop_type.name.as_str() == "pombase_gene_id" {
                        pombase_gene_id = Some(cvtermprop.value.clone());
                        let gene_for_map = format!("{}:{}", self.config.database_name,
                                                   cvtermprop.value);
                        pro_term_to_gene.insert(cvterm.termid().to_string(), gene_for_map);
                    }
                }

                let mut xrefs = HashMap::new();

                for (source_name, source_config) in cv_config.source_config {
                    let mut maybe_xref_id = None;
                    if let Some(ref term_xref_id_prop) = source_config.id_source {
                        if let Some(term_xref_id_prop) = term_xref_id_prop.strip_prefix("prop_name:") {
                            for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                                if cvtermprop.prop_type.name == *term_xref_id_prop {
                                    maybe_xref_id = Some(cvtermprop.value.clone());
                                    break;
                                }
                            }
                        } else if term_xref_id_prop == "ACCESSION" {
                            let dbxref: &Dbxref = cvterm.dbxref.borrow();
                            maybe_xref_id = Some(dbxref.accession.clone());
                        }
                    }
                    let mut maybe_xref_display_name = None;
                    if let Some(ref xref_display_name_prop) = source_config.display_name_prop {
                        for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                            if cvtermprop.prop_type.name == *xref_display_name_prop {
                                maybe_xref_display_name = Some(cvtermprop.value.clone());
                            }
                        }
                    }
                    if let Some(xref_id) = maybe_xref_id {
                        let term_xref = TermXref {
                            xref_id,
                            xref_display_name: maybe_xref_display_name,
                        };

                        xrefs.insert(source_name.clone(), term_xref);
                    }
                }

                let synonyms =
                    cvterm.cvtermsynonyms.borrow().iter().map(|syn| {
                        SynonymDetails {
                            synonym_type: syn.synonym_type.name.clone(),
                            name: syn.name.clone(),
                            reference: None,
                        }
                    }).collect::<Vec<_>>();

                let definition_xrefs =
                    cvterm.definition_xrefs.borrow().iter()
                    .map(|dbxref| {
                        dbxref.identifier()
                    }).collect::<HashSet<_>>();

                let secondary_identifiers =
                    cvterm.other_dbxrefs.borrow().iter()
                    .map(|dbxref| {
                        dbxref.identifier()
                    }).collect::<HashSet<_>>();

                self.terms.insert(cvterm.termid(),
                                  TermDetails {
                                      name: cvterm.name.clone(),
                                      cv_name: cvterm.cv.name.clone(),
                                      annotation_feature_type,
                                      interesting_parent_ids: HashSet::new(),
                                      interesting_parent_details: HashSet::new(),
                                      in_subsets: HashSet::new(),
                                      termid: cvterm.termid(),
                                      synonyms,
                                      definition: cvterm.definition.clone(),
                                      direct_ancestors: vec![],
                                      definition_xrefs,
                                      secondary_identifiers,
                                      annotated_genes: HashSet::new(),
                                      single_locus_annotated_genes: HashSet::new(),
                                      multi_locus_annotated_genes: HashSet::new(),
                                      is_obsolete: cvterm.is_obsolete,
                                      single_locus_genotype_uniquenames: HashSet::new(),
                                      cv_annotations: HashMap::new(),
                                      genes_by_uniquename: HashMap::new(),
                                      genotypes_by_uniquename: HashMap::new(),
                                      alleles_by_uniquename: HashMap::new(),
                                      transcripts_by_uniquename: HashMap::new(),
                                      references_by_uniquename: HashMap::new(),
                                      terms_by_termid: HashMap::new(),
                                      annotation_details: HashMap::new(),
                                      double_mutant_genetic_interactions: HashMap::new(),
                                      single_allele_genetic_interactions: HashMap::new(),
                                      gene_count: 0,
                                      genotype_count: 0,
                                      xrefs,
                                      pombase_gene_id,
                                      gocams: HashSet::new(),
                                  });
                self.term_ids_by_name.insert(cvterm.name.clone(), cvterm.termid());
            }
        }

        self.pro_term_to_gene = pro_term_to_gene;
    }

    fn get_ext_rel_display_name(&self, annotation_termid: &FlexStr,
                                ext_rel_name: &FlexStr) -> FlexStr {
        if let Some(ext_conf) = self.matching_ext_config(annotation_termid, ext_rel_name) {
            ext_conf.display_name
        } else {
            str::replace(ext_rel_name, "_", " ").to_shared_str()
        }
    }

    fn process_extension_cvterms(&mut self) {
        let db_prefix = format!("{}:", self.config.database_name);

        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                    if cvtermprop.prop_type.name.starts_with(ANNOTATION_EXT_REL_PREFIX) {
                        let ext_rel_name_str =
                            &cvtermprop.prop_type.name[ANNOTATION_EXT_REL_PREFIX.len()..];
                        let ext_rel_name = ext_rel_name_str.to_shared_str();
                        let ext_range = cvtermprop.value.clone();
                        let range: ExtRange = if ext_range.starts_with(&db_prefix) {
                            let db_feature_uniquename = ext_range[db_prefix.len()..].to_shared_str();
                            if let Some(captures) = PROMOTER_RE.captures(&db_feature_uniquename) {
                                let gene_uniquename = captures["gene"].to_shared_str();
                                if self.genes.contains_key(&gene_uniquename) {
                                    ExtRange::Promoter(gene_uniquename)
                                } else {
                                    panic!("unknown gene in promoter: {}", db_feature_uniquename);
                                }
                            } else if self.genes.contains_key(&db_feature_uniquename) {
                                ExtRange::Gene(db_feature_uniquename.clone())
                            } else if let Some(captures) = TRANSCRIPT_ID_RE.captures(db_feature_uniquename.as_ref()) {
                                if self.genes.contains_key(&captures["gene"].to_shared_str()) {
                                    ExtRange::Transcript(db_feature_uniquename.clone())
                                } else {
                                    panic!("unknown gene for transcript: {}", db_feature_uniquename);
                                }
                            } else {
                                panic!("can't find gene or transcript for: {}", db_feature_uniquename);
                            }
                        } else {
                            ExtRange::Misc(ext_range)
                        };

                        if let Some(base_termid) =
                            self.base_term_of_extensions.get(&cvterm.termid()) {
                                let rel_type_display_name =
                                    self.get_ext_rel_display_name(base_termid, &ext_rel_name);

                                let rel_type_id =
                                    self.term_ids_by_name.get(&ext_rel_name).cloned();

                                self.parts_of_extensions.entry(cvterm.termid())
                                    .or_default().push(ExtPart {
                                        rel_type_id,
                                        rel_type_name: ext_rel_name,
                                        rel_type_display_name,
                                        ext_range: range,
                                    });
                            } else {
                                panic!("can't find details for term: {}\n", cvterm.termid());
                            }
                    }
                }
            }
        }
    }

    fn make_term_ext_range(&self, term_id: &FlexStr) -> ExtRange {
       if term_id.starts_with("PR:") {
         if let Some(term_details) = self.terms.get(term_id) {
           if let Some(ref pombase_gene_id) = term_details.pombase_gene_id {
             let gene_and_product = GeneAndGeneProduct {
                gene_uniquename: pombase_gene_id.clone(),
                product: term_id.clone(),
             };
             ExtRange::GeneAndGeneProduct(gene_and_product)
           } else {
             ExtRange::GeneProduct(term_id.clone())
           }
         } else {
           ExtRange::GeneProduct(term_id.clone())
         }
       } else {
         ExtRange::Term(term_id.clone())
       }
    }

    fn process_cvterm_rels(&mut self) {
        for cvterm_rel in &self.raw.cvterm_relationships {
            let subject_term = &cvterm_rel.subject;
            let object_term = &cvterm_rel.object;
            let rel_type = &cvterm_rel.rel_type;

            if subject_term.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                let subject_termid = subject_term.termid();
                if rel_type.name == "is_a" {
                    self.base_term_of_extensions.insert(subject_termid.clone(),
                                                        object_term.termid().clone());
                }
            } else {
                let object_term_short =
                    self.make_term_short(&object_term.termid());
                if let Some(ref mut subject_term_details) = self.terms.get_mut(&subject_term.termid()) {
                    subject_term_details.direct_ancestors.push(TermAndRelation {
                        termid: object_term_short.termid.clone(),
                        term_name: object_term_short.name.clone(),
                        relation_name: rel_type.name.clone(),
                    });
                }
            }
        }

        for cvterm_rel in &self.raw.cvterm_relationships {
            let subject_term = &cvterm_rel.subject;
            let object_term = &cvterm_rel.object;
            let rel_type = &cvterm_rel.rel_type;

            if subject_term.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                let subject_termid = subject_term.termid();
                if rel_type.name != "is_a" {
                    let object_termid = object_term.termid();
                    if let Some(base_termid) =
                        self.base_term_of_extensions.get(&subject_term.termid()) {
                            let rel_type_display_name =
                                self.get_ext_rel_display_name(base_termid, &rel_type.name);

                            let ext_range =
                               self.make_term_ext_range(&object_termid);

                            self.parts_of_extensions.entry(subject_termid)
                                .or_default().push(ExtPart {
                                    rel_type_id: Some(rel_type.termid()),
                                    rel_type_name: rel_type.name.clone(),
                                    rel_type_display_name,
                                    ext_range,
                                });
                        } else {
                            panic!("can't find details for {}\n", object_termid);
                        }
                }
            }
        }
    }

    fn process_feature_synonyms(&mut self) {
        for feature_synonym in &self.raw.feature_synonyms {
            let feature = &feature_synonym.feature;
            let synonym = &feature_synonym.synonym;

            let reference =
                if feature_synonym.publication.uniquename == "null" {
                    None
                } else {
                    Some(feature_synonym.publication.uniquename.clone())
                };
            let make_synonym = || {
                SynonymDetails {
                    name: synonym.name.clone(),
                    synonym_type: synonym.synonym_type.name.clone(),
                    reference,
                }
            };

            if let Some(ref mut gene_details) = self.genes.get_mut(&feature.uniquename) {
                gene_details.synonyms.push(make_synonym());
            } else if let Some(ref mut allele) = self.alleles.get_mut(&feature.uniquename) {
                let synonym = make_synonym();
                if let Err(insert_pos) = allele.synonyms.binary_search(&synonym) {
                    // keep synonym list ordered
                    allele.synonyms.insert(insert_pos, synonym);
                }
            }
        }
    }

    fn process_feature_publications(&mut self) {
        for feature_pub in &self.raw.feature_pubs {
            let feature = &feature_pub.feature;
            let publication = &feature_pub.publication;

            if publication.uniquename.starts_with("PMID:") {
                let Some(source) = feature_pub.feature_pubprops.borrow()
                    .iter().filter(|prop| prop.prop_type.name == "feature_pub_source")
                    .map(|prop| prop.value.clone().unwrap())
                    .next()
                else {
                    continue;
                };

                if let Some(ref mut gene_details) = self.genes.get_mut(&feature.uniquename) {
                    let ref_and_source = ReferenceAndSource {
                        reference_uniquename: publication.uniquename.clone(),
                        source: source.clone(),
                    };
                    gene_details.feature_publications.insert(ref_and_source);
                }
            }
        }
    }

    fn make_genotype_short(&self, genotype_display_uniquename: &FlexStr) -> GenotypeShort {
        if let Some(details) = self.genotypes.get(genotype_display_uniquename) {
            details.into()
        } else {
            panic!("can't find genotype {}", genotype_display_uniquename);
        }
    }

    fn make_allele_short(&self, allele_uniquename: &FlexStr) -> AlleleShort {
        if let Some(details) = self.alleles.get(allele_uniquename) {
            details.into()
        } else {
            panic!("can't find allele for {}", allele_uniquename);
        }
    }

    fn add_product_to_protein(&mut self, transcript_uniquename: &FlexStr,
                              product: FlexStr) {
        if let Some(transcript_details) =
            self.transcripts.get_mut(transcript_uniquename)
            && let Some(ref mut protein) = transcript_details
                .protein
            {
                protein.product = Some(product);
            }
    }

    // process feature properties stored as cvterms,
    // eg. characterisation_status and product
    fn process_props_from_feature_cvterms(&mut self) {
        for feature_cvterm in &self.raw.feature_cvterms {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;

            let (maybe_gene_uniquename, maybe_transcript_uniquename) =
                if cvterm.cv.name == "PomBase gene products" {
                    if feature.feat_type.name == "polypeptide" {
                        if let Some(transcript_uniquename) =
                            self.transcripts_of_polypeptides.get(&feature.uniquename) {
                                if let Some(gene_uniquename) =
                                    self.genes_of_transcripts.get(transcript_uniquename) {
                                        (Some(gene_uniquename.clone()),
                                         Some(transcript_uniquename.clone()))
                                    } else {
                                        (None, None)
                                    }
                            } else {
                                (None, None)
                            }
                    } else if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                        if let Some(gene_uniquename) =
                            self.genes_of_transcripts.get(&feature.uniquename) {
                                (Some(gene_uniquename.clone()), Some(feature.uniquename.clone()))
                            } else {
                                (None, None)
                            }
                    } else if feature.feat_type.name == "gene" {
                        (Some(feature.uniquename.clone()), None)
                    } else {
                        (None, None)
                    }
                } else {
                    (None, None)
                };

            if let Some(gene_uniquename) = maybe_gene_uniquename
                && let Some(transcript_uniquename) = maybe_transcript_uniquename {
                    self.add_gene_product(&gene_uniquename, &cvterm.name);

                    self.add_product_to_protein(&transcript_uniquename,
                                                cvterm.name.clone());
                }

            if feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene" {
                if cvterm.cv.name == "PomBase gene characterisation status" {
                    self.add_characterisation_status(&feature.uniquename, &cvterm.name);
                } else if cvterm.cv.name == "name_description" {
                    self.add_name_description(&feature.uniquename, &cvterm.name);
                }
            }
        }
    }

    fn make_with_or_from_value(&self, with_or_from_value: &FlexStr) -> WithFromValue {
        if let Some(captures) = PREFIX_AND_ID_RE.captures(with_or_from_value) {
            let prefix = captures["prefix"].to_shared_str();
            let id = captures["id"].to_shared_str();

            if self.genes.contains_key(&id) {
                let gene_short = self.make_gene_short(&id);
                if self.config.database_name == prefix {
                    // a gene from the main organism
                    return WithFromValue::Gene(gene_short);
                } else if let Some(name) = &gene_short.name {
                    return WithFromValue::IdentifierAndName({
                        IdentifierAndName {
                            identifier: with_or_from_value.clone(),
                            name: name.clone(),
                        }
                    });
                }
            } else if self.transcripts.contains_key(&id)
            && self.config.database_name == prefix {
                return WithFromValue::Transcript(id);
            }
        } else if self.genes.contains_key(with_or_from_value) {
            let gene_short = self.make_gene_short(with_or_from_value);
            // a gene from the main organism
            return WithFromValue::Gene(gene_short);
        } else if self.transcripts.contains_key(with_or_from_value) {
            return WithFromValue::Transcript(with_or_from_value.clone());
        }

        if self.terms.contains_key(with_or_from_value) {
            return WithFromValue::Term(self.make_term_short(with_or_from_value))
        }

        WithFromValue::Identifier(with_or_from_value.clone())
    }

    // process annotation
    fn process_feature_cvterms(&mut self) {
        let rel_order = self.config.extension_relation_order.clone();

        'FEATURE_CVTERM:
        for feature_cvterm in &self.raw.feature_cvterms {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;

            if feature.type_name() == "gocam_model" {
                continue;
            }

            let termid =
                match self.base_term_of_extensions.get(&cvterm.termid()) {
                    Some(base_termid) => base_termid.clone(),
                    None => cvterm.termid(),
                };

            let mut transcript_uniquenames = vec![];

            let mut extension = vec![];

            if cvterm.cv.name == "PomBase gene characterisation status" ||
                cvterm.cv.name == "PomBase gene products" ||
                cvterm.cv.name == "name_description" {
                    continue;
                }

            let publication = &feature_cvterm.publication;
            let mut extra_props: HashMap<FlexStr, FlexStr> = HashMap::new();
            let mut conditions: HashSet<TermId> = HashSet::new();
            let mut condition_details: BTreeSet<(TermId, Option<String>)> = BTreeSet::new();
            let mut withs: HashSet<WithFromValue> = HashSet::new();
            let mut froms: HashSet<WithFromValue> = HashSet::new();
            let mut qualifiers: Vec<Qualifier> = vec![];
            let mut date: Option<FlexStr> = None;
            let mut assigned_by: Option<FlexStr> = None;
            let mut evidence: Option<FlexStr> = None;
            let mut eco_evidence: Option<FlexStr> = None;
            let mut annotation_phenotype_score: Option<FlexStr> = None;
            let mut genotype_background: Option<FlexStr> = None;
            let mut allele_promoters = vec![];
            let mut throughput: Option<Throughput> = None;
            let mut curator_orcid: Option<CuratorOrcid> = None;
            let mut submitter_comment: Option<FlexStr> = None;
            let mut curation_session: Option<CurationSessionKey> = None;

            // need to get evidence first as it's used later
            // See: https://github.com/pombase/website/issues/455
            for prop in feature_cvterm.feature_cvtermprops.borrow().iter() {
                if &prop.type_name() == "evidence"
                    && let Some(ref evidence_long) = prop.value {
                        for (evidence_code, ev_details) in &self.config.evidence_types {
                            if &ev_details.long == evidence_long {
                                evidence = Some(evidence_code.clone());
                            }
                        }
                        if evidence.is_none() {
                            evidence = Some(evidence_long.clone());
                        }
                    }
            }

            for prop in feature_cvterm.feature_cvtermprops.borrow().iter() {
                match &prop.type_name() as &str {
                    "residue" | "scale" | "gene_product_form_id" |
                    "quant_gene_ex_copies_per_cell" |
                    "quant_gene_ex_avg_copies_per_cell" => {
                        if let Some(value) = prop.value.clone() {
                            if prop.type_name() == "residue" &&
                                &cvterm.cv.name != "sequence"
                            {
                                let residue = value.clone();
                                let rel_type_name = flex_str!("residue");
                                let display_name =
                                    self.get_ext_rel_display_name(&termid, &rel_type_name);

                                let residue_range_part = ExtPart {
                                    rel_type_id: None,
                                    rel_type_name,
                                    rel_type_display_name: display_name,
                                    ext_range: ExtRange::ModifiedResidues(vec![residue]),
                                };
                                extension.insert(0, residue_range_part);
                            }

                            extra_props.insert(prop.type_name().clone(), value);
                        }
                    },
                    "condition" =>
                        if let Some(value) = prop.value.clone() {
                            if value.contains(':') {
                                conditions.insert(value.clone());
                            } else {
                                eprintln!(r#"ignoring condition that isn't a term ID "{}" (from annotation of {} with {})"#,
                                          value, feature.uniquename, termid);
                            }
                        },
                    "condition_detail" =>
                        if let Some(ref value) = prop.value {
                            if value.contains(':') {
                                let parsed_value = parse_condition_with_detail(value);
                                condition_details.insert(parsed_value);
                            } else {
                                eprintln!(r#"ignoring condition that doesn't contain a term ID "{}" (from annotation of {} with {})"#,
                                          value, feature.uniquename, termid);
                            }
                        },
                    "qualifier" =>
                        if let Some(value) = prop.value.clone() {
                            qualifiers.push(value);
                        },
                    "assigned_by" =>
                        if let Some(value) = prop.value.clone() {
                            assigned_by = Some(value);
                        },
                    "date" => {
                        if let Some(value) = prop.value.clone() {
                            date = Some(value);
                        }
                    },
                    "eco_evidence" => {
                        eco_evidence = prop.value.clone();
                    },
                    "annotation_phenotype_score" => {
                        annotation_phenotype_score = prop.value.clone();
                    }
                    "with" => {
                        if let Some(value) = prop.value.clone() {
                            withs.insert(self.make_with_or_from_value(&value));
                        }
                    },
                    "from" => {
                        if let Some(value) = prop.value.clone() {
                            froms.insert(self.make_with_or_from_value(&value));
                        }
                    },
                    "curator_orcid" => {
                        curator_orcid = prop.value.clone();
                    },
                    "submitter_comment" => {
                        submitter_comment = prop.value.clone();
                    },
                    "canto_session" => {
                        curation_session = prop.value.clone();
                    },
                    "annotation_throughput_type" => {
                        if let Some(throughput_type) = prop.value.clone() {
                            throughput = Some(match throughput_type.as_ref() {
                                "low throughput" => Throughput::LowThroughput,
                                "high throughput" => Throughput::HighThroughput,
                                "non-experimental" => Throughput::NonExperimental,
                                _ => {
                                    panic!("unknown throughput type: {}",
                                           throughput_type);
                                }
                            });
                        }
                    },
                    _ => ()
                }
            }

            let mut maybe_genotype_uniquename = None;
            let mut gene_uniquenames_vec: Vec<GeneUniquename> =
                match &feature.feat_type.name as &str {
                    "polypeptide" => {
                        if let Some(transcript_uniquename) =
                            self.transcripts_of_polypeptides.get(&feature.uniquename) {
                                if let Some(gene_uniquename) =
                                    self.genes_of_transcripts.get(transcript_uniquename) {
                                        vec![gene_uniquename.clone()]
                                    } else {
                                        vec![]
                                    }
                            } else {
                                vec![]
                            }
                    },
                    "genotype" => {
                        let loci: Vec<_> =
                            self.loci_of_genotypes[&feature.uniquename]
                            .values().cloned().collect();
                        let genotype_display_uniquename =
                            make_genotype_display_uniquename(&loci, &self.alleles);
                        maybe_genotype_uniquename = Some(genotype_display_uniquename);
                        genotype_background =
                            self.genotype_backgrounds.get(&feature.uniquename)
                            .cloned();

                        loci.iter()
                            .map(|locus| {
                                locus.expressed_alleles.iter()
                                    .map(|expressed_allele| {
                                        let allele_details =
                                            self.alleles.get(&expressed_allele.allele_uniquename).unwrap();
                                        let allele_short: AlleleShort = allele_details.into();
                                        let allele_display_name = allele_short.short_display_name();
                                        let allele_uniquename =
                                            expressed_allele.allele_uniquename.clone();
                                        let allele_expression = expressed_allele.expression.clone();
                                        if let Some(ref promoter_gene) = expressed_allele.promoter_gene {
                                            allele_promoters.push(AnnotationPromoter {
                                                allele_uniquename,
                                                allele_display_name,
                                                allele_expression,
                                                allele_gene: allele_details.gene.clone(),
                                                promoter: promoter_gene.into(),
                                            })
                                        }
                                        allele_details.gene.uniquename.clone()
                                    })
                                    .collect()
                            })
                            .collect::<Vec<Vec<_>>>()
                            .concat()
                    },
                    "gene" | "pseudogene" => {
                        vec![feature.uniquename.clone()]
                    },
                    "genotype_interaction" => {
                        vec![]
                    },
                    "gocam_model" => {
                        vec![]
                    },
                    _ =>
                        if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                            if let Some(gene_uniquename) =
                                self.genes_of_transcripts.get(&feature.uniquename) {
                                    if let Some(gene_details) = self.genes.get(gene_uniquename)
                                        && gene_details.transcripts.len() > 1 {
                                            // only bother to record the specific transcript if
                                            // there is more than one
                                            transcript_uniquenames.push(feature.uniquename.clone());
                                        }
                                    vec![gene_uniquename.clone()]
                                } else {
                                    vec![]
                                }
                        } else {
                            eprintln!("can't handle annotation on {} {}",
                                      &feature.feat_type.name, &feature.uniquename);
                            continue 'FEATURE_CVTERM;
                        }
                };

            gene_uniquenames_vec.dedup();

            gene_uniquenames_vec =
                gene_uniquenames_vec.iter().map(|gene_uniquename: &FlexStr| {
                    self.make_gene_short(gene_uniquename).uniquename
                }).collect();

            let reference_uniquename =
                if publication.uniquename == "null" {
                    None
                } else {
                    Some(publication.uniquename.clone())
                };

            let mut extra_props_clone = extra_props.clone();
            let copies_per_cell = extra_props_clone.remove(&flex_str!("quant_gene_ex_copies_per_cell"));
            let avg_copies_per_cell = extra_props_clone.remove(&flex_str!("quant_gene_ex_avg_copies_per_cell"));
            let gene_ex_props =
                if copies_per_cell.is_some() || avg_copies_per_cell.is_some() {
                    let scale = extra_props_clone.remove(&flex_str!("scale"))
                        .expect("gene ex scale missing");
                    Some(GeneExProps {
                        copies_per_cell,
                        avg_copies_per_cell,
                        scale,
                    })
                } else {
                    None
                };

            if gene_uniquenames_vec.len() > 1 && maybe_genotype_uniquename.is_none() {
                panic!("non-genotype annotation has more than one gene");
            }

            let annotation_detail = OntAnnotationDetail {
                id: feature_cvterm.feature_cvterm_id,
                genes: gene_uniquenames_vec,
                transcript_uniquenames,
                reference: reference_uniquename,
                curation_session,
                genotype: maybe_genotype_uniquename,
                genotype_background,
                allele_promoters,
                withs,
                froms,
                residue: extra_props_clone.remove(&flex_str!("residue")),
                gene_product_form_id: extra_props_clone.remove(&flex_str!("gene_product_form_id")),
                gene_ex_props,
                qualifiers,
                evidence,
                eco_evidence,
                annotation_phenotype_score,
                conditions,
                condition_details,
                extension,
                date,
                assigned_by,
                throughput,
                curator: curator_orcid,
                submitter_comment,
            };

            if &feature.feat_type.name == "genotype_interaction" {
                self.add_genetic_interaction(feature, &rel_order,
                                             cvterm.borrow(), annotation_detail);
            } else {
                self.add_annotation(&rel_order,
                                    cvterm.borrow(), feature_cvterm.is_not,
                                    annotation_detail);
            }
        }
    }

    fn make_term_annotations(&self, termid: &FlexStr, detail_ids: &[OntAnnotationId],
                             is_not: bool)
                       -> Vec<(CvName, OntTermAnnotations)> {
        let term_details = &self.terms[termid];

        let cv_name = term_details.cv_name.clone();

        match cv_name.as_ref() {
            "gene_ex" | "PomGeneExRNA" | "PomGeneExProt" | "PomGeneExRD" => {
                if is_not {
                    panic!("gene_ex annotations can't be NOT annotations");
                }
                let mut qual_annotations =
                    OntTermAnnotations {
                        term: termid.clone(),
                        is_not: false,
                        rel_names: HashSet::new(),
                        annotations: vec![],
                        summary: None,
                    };
                let mut quant_annotations =
                    OntTermAnnotations {
                        term: termid.clone(),
                        is_not: false,
                        rel_names: HashSet::new(),
                        annotations: vec![],
                        summary: None,
                    };
                for annotation_id in detail_ids {
                    let annotation = self.annotation_details.
                        get(annotation_id).expect("can't find OntAnnotationDetail");

                    if annotation.gene_ex_props.is_some() {
                        quant_annotations.annotations.push(*annotation_id)
                    } else {
                        qual_annotations.annotations.push(*annotation_id)
                    }
                }

                let mut return_vec = vec![];

                if !qual_annotations.annotations.is_empty() {
                    return_vec.push((flex_str!("qualitative_gene_expression"),
                                     qual_annotations));
                }

                if !quant_annotations.annotations.is_empty() {
                    return_vec.push((flex_str!("quantitative_gene_expression"),
                                     quant_annotations));
                }

                return_vec
            },
            "fission_yeast_phenotype" => {
                let mut single_locus =
                    OntTermAnnotations {
                        term: termid.clone(),
                        is_not,
                        rel_names: HashSet::new(),
                        annotations: vec![],
                        summary: None,
                    };
                let mut multi_locus =
                    OntTermAnnotations {
                        term: termid.clone(),
                        is_not,
                        rel_names: HashSet::new(),
                        annotations: vec![],
                        summary: None,
                    };

                for annotation_id in detail_ids {
                    let annotation = self.annotation_details.
                        get(annotation_id).expect("can't find OntAnnotationDetail");

                    let genotype_uniquename = annotation.genotype.as_ref().unwrap();

                    if let Some(genotype_details) = self.genotypes.get(genotype_uniquename) {
                        if genotype_details.loci.len() == 1 {
                            single_locus.annotations.push(*annotation_id);
                        } else if !multi_locus.annotations.contains(annotation_id) {
                            multi_locus.annotations.push(*annotation_id);
                        }
                    } else {
                        panic!("can't find genotype details for {}\n", genotype_uniquename);
                    }
                }

                let mut return_vec = vec![];

                if !single_locus.annotations.is_empty() {
                    return_vec.push((flex_str!("single_locus_phenotype"),
                                     single_locus));
                }

                if !multi_locus.annotations.is_empty() {
                    return_vec.push((flex_str!("multi_locus_phenotype"),
                                     multi_locus));
                }

                return_vec
            },
            _ => {
                vec![(cv_name,
                      OntTermAnnotations {
                          term: termid.clone(),
                          is_not,
                          rel_names: HashSet::new(),
                          annotations: detail_ids.to_owned(),
                          summary: None,
                      })]
            }
        }
    }


    fn remove_duplicate_transcript_annotation(&mut self) {
        let ont_annotation_map = &mut self.all_ont_annotations;

        for annotations in ont_annotation_map.values_mut() {
            let (no_transcript_annotations, mut has_transcript_annotations): (Vec<i32>, Vec<i32>) =
                annotations
                .iter()
                .partition(|&annotation_id| {
                    if let Some(ont_annotation_detail) =
                        self.annotation_details.get(annotation_id) {
                            ont_annotation_detail.transcript_uniquenames.is_empty()
                        } else {
                            panic!("can't find annotation details for {}", annotation_id);
                        }
                });

            *annotations = no_transcript_annotations;

            if has_transcript_annotations.len() >= 2 {
                // merge annotations that differ only by transcript ID

                has_transcript_annotations.sort();

                let mut prev_annotation_id = has_transcript_annotations.remove(0);

                for current_annotation_id in has_transcript_annotations.drain(0..) {
                    let (annotations_equal, current_transcript_uniquename) = {
                        let prev_annotation =
                            self.annotation_details.get(&prev_annotation_id).unwrap();
                        let current_annotation =
                            self.annotation_details.get(&current_annotation_id).unwrap();
                        (prev_annotation == current_annotation,
                         current_annotation.transcript_uniquenames[0].clone())
                    };
                    if annotations_equal {
                        if let Some(annotation_details) = self.annotation_details.get(&prev_annotation_id)
                            && !annotation_details.transcript_uniquenames.contains(&current_transcript_uniquename) {
                                self.annotation_details.get_mut(&prev_annotation_id).unwrap()
                                    .transcript_uniquenames.push(current_transcript_uniquename);
                            }
                    } else {
                        annotations.push(prev_annotation_id);
                        prev_annotation_id = current_annotation_id;
                    }
                }
                annotations.push(prev_annotation_id);
            } else {
                annotations.extend(has_transcript_annotations.iter());
            }
        }
    }

    // store the OntTermAnnotations in the TermDetails, GeneDetails,
    // GenotypeDetails and ReferenceDetails
    fn store_ont_annotations(&mut self, is_not: bool) {
        let ont_annotation_map = if is_not {
            &self.all_not_ont_annotations
        } else {
            &self.all_ont_annotations
        };

        let mut gene_annotation_by_term: HashMap<GeneUniquename, HashMap<TermId, Vec<OntAnnotationId>>> =
            HashMap::new();
        let mut genotype_annotation_by_term: HashMap<GenotypeUniquename, HashMap<TermId, Vec<OntAnnotationId>>> =
            HashMap::new();
        let mut ref_annotation_by_term: HashMap<FlexStr, HashMap<TermId, Vec<OntAnnotationId>>> =
            HashMap::new();

        let mut ont_annotations = vec![];

        for (termid, annotations) in ont_annotation_map {
                let new_annotations =
                    self.make_term_annotations(termid, annotations, is_not);

                if let Some(ref mut term_details) = self.terms.get_mut(termid) {
                    for (cv_name, new_annotation) in new_annotations {
                        term_details.cv_annotations.entry(cv_name.clone())
                            .or_insert_with(Vec::new)
                            .push(new_annotation);
                    }
                } else {
                    panic!("missing termid: {}\n", termid);
                }

            for annotation_id in annotations {
                let annotation = self.annotation_details.
                    get(annotation_id).expect("can't find OntAnnotationDetail");

                for gene_uniquename in &annotation.genes {
                    gene_annotation_by_term.entry(gene_uniquename.clone())
                        .or_default()
                        .entry(termid.clone())
                        .or_default()
                        .push(*annotation_id);
                }

                if let Some(ref genotype_uniquename) = annotation.genotype {
                    let existing =
                        genotype_annotation_by_term.entry(genotype_uniquename.clone())
                        .or_default()
                        .entry(termid.clone())
                        .or_default();
                    if !existing.contains(annotation_id) {
                        existing.push(*annotation_id);
                    }
                }

                if let Some(reference_uniquename) = annotation.reference.clone() {
                    ref_annotation_by_term.entry(reference_uniquename)
                        .or_default()
                        .entry(termid.clone())
                        .or_default()
                        .push(*annotation_id);
                }

                for condition_termid in &annotation.conditions {
                    let cv_name =
                        if let Some(term_details) = self.terms.get(condition_termid) {
                            term_details.cv_name.clone()
                        } else {
                            panic!("can't find term details for {}", condition_termid);
                        };

                    if let Some(ref mut condition_term_details) =
                        self.terms.get_mut(&condition_termid.clone())
                    {
                        condition_term_details.cv_annotations
                            .entry(cv_name.clone())
                            .or_insert({
                                let mut new_vec = Vec::new();
                                let new_term_annotation =
                                    OntTermAnnotations {
                                        term: condition_termid.clone(),
                                        is_not,
                                        rel_names: HashSet::new(),
                                        annotations: vec![],
                                        summary: None,
                                    };
                                new_vec.push(new_term_annotation);
                                new_vec
                            });
                        condition_term_details.cv_annotations.get_mut(&cv_name)
                            .unwrap()[0]
                            .annotations.push(*annotation_id);
                    }
                }

                /*

                Remove for now because it's messing with the gene counts.
                See: https://github.com/pombase/website/issues/1705

                // Add annotations to terms referred to in extensions.  They
                // are added to fake CV that have a name starting with
                // "extension:".  The CV name will end with ":genotype" if the
                // annotation is a phentoype/genotype, and will end with ":end"
                // otherwise.  The middle of the fake CV name is the display
                // name for the extension relation.
                // eg. "extension:directly activates:gene"
                for ext_part in &annotation.extension {
                    if let ExtRange::Term(ref part_termid) = ext_part.ext_range {
                        let cv_name = "extension:".to_owned() + &ext_part.rel_type_display_name;

                        if let Some(ref mut part_term_details) =
                            self.terms.get_mut(part_termid)
                        {
                            let extension_cv_name =
                                if annotation.genotype.is_some() {
                                    cv_name.clone() + ":genotype"
                                } else {
                                    cv_name.clone() + ":gene"
                                };

                            part_term_details.cv_annotations
                                .entry(flex_str!(&extension_cv_name))
                                .or_insert({
                                    let mut new_vec = Vec::new();
                                    let new_term_annotation =
                                        OntTermAnnotations {
                                            term: part_termid.to_owned(),
                                            is_not,
                                            rel_names: HashSet::new(),
                                            annotations: vec![],
                                            summary: None,
                                        };
                                    new_vec.push(new_term_annotation);
                                    new_vec
                                });
                            part_term_details.cv_annotations.get_mut(&extension_cv_name)
                                .unwrap()[0]
                                .annotations.push(annotation_id);
                        }
                    }
                }

                */

                let gene_short_list =
                    annotation.genes.iter().map(|uniquename: &FlexStr| {
                        self.make_gene_short(uniquename)
                    }).collect::<HashSet<_>>();

                let reference_short =
                    annotation.reference.as_ref().and_then(|uniquename: &FlexStr| {
                        make_reference_short(&self.references, uniquename)
                    });

                let genotype_short =
                    annotation.genotype.as_ref().map(|uniquename: &FlexStr| {
                        self.make_genotype_short(uniquename)
                    });

                let conditions =
                    annotation.conditions.iter().map(|termid| {
                        self.make_term_short(termid)
                    }).collect::<HashSet<_>>();

                if gene_short_list.is_empty() {
                    panic!("no genes for {:?}", &annotation);
                }

                let ont_annotation = OntAnnotation {
                    term_short: self.make_term_short(termid),
                    id: annotation.id,
                    genes: gene_short_list,
                    reference_short,
                    genotype_short,
                    genotype_background: annotation.genotype_background.clone(),
                    withs: annotation.withs.clone(),
                    froms: annotation.froms.clone(),
                    residue: annotation.residue.clone(),
                    gene_ex_props: annotation.gene_ex_props.clone(),
                    qualifiers: annotation.qualifiers.clone(),
                    evidence: annotation.evidence.clone(),
                    conditions,
                    extension: annotation.extension.clone(),
                    assigned_by: annotation.assigned_by.clone(),
                };

                ont_annotations.push(ont_annotation);
            }
        }

        let mut term_names = HashMap::new();
        for (termid, term_details) in &self.terms {
            term_names.insert(termid.clone(), term_details.name.to_lowercase());
        }

        let ont_term_cmp = |ont_term_1: &OntTermAnnotations, ont_term_2: &OntTermAnnotations| {
            if !ont_term_1.is_not && ont_term_2.is_not {
                return Ordering::Less;
            }
            if ont_term_1.is_not && !ont_term_2.is_not {
                return Ordering::Greater;
            }
            let term1 = &term_names[&ont_term_1.term];
            let term2 = &term_names[&ont_term_2.term];

            term1.cmp(term2)
        };

        for (gene_uniquename, term_annotation_map) in &gene_annotation_by_term {
            for (termid, details) in term_annotation_map {
                let new_annotations =
                    self.make_term_annotations(termid, details, is_not);

                let gene_details = self.genes.get_mut(gene_uniquename).unwrap();

                for (cv_name, new_annotation) in new_annotations {
                    gene_details.cv_annotations.entry(cv_name.clone())
                        .or_default()
                        .push(new_annotation);
                }
            }

            let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
            for cv_annotations in gene_details.cv_annotations.values_mut() {
                cv_annotations.sort_by(&ont_term_cmp)
            }
        }

        for (genotype_uniquename, term_annotation_map) in &genotype_annotation_by_term {
            for (termid, details) in term_annotation_map {
                let new_annotations =
                    self.make_term_annotations(termid, details, is_not);

                let details = self.genotypes.get_mut(genotype_uniquename).unwrap();

                for (cv_name, new_annotation) in new_annotations {
                    details.cv_annotations.entry(cv_name.clone())
                        .or_default()
                        .push(new_annotation);
                }
            }

            let details = self.genotypes.get_mut(genotype_uniquename).unwrap();
            for cv_annotations in details.cv_annotations.values_mut() {
                cv_annotations.sort_by(&ont_term_cmp)
            }
        }

        for (reference_uniquename, ref_annotation_map) in &ref_annotation_by_term {
            for (termid, details) in ref_annotation_map {
                let new_annotations =
                    self.make_term_annotations(termid, details, is_not);

                let ref_details = self.references.get_mut(reference_uniquename).unwrap();

                for (cv_name, new_annotation) in new_annotations {
                    ref_details.cv_annotations.entry(cv_name).or_default()
                        .push(new_annotation.clone());
                }
            }

            let ref_details = self.references.get_mut(reference_uniquename).unwrap();
            for cv_annotations in ref_details.cv_annotations.values_mut() {
                cv_annotations.sort_by(&ont_term_cmp)
            }
        }

        for ont_annotation in ont_annotations.drain(0..) {
            self.ont_annotations.push(ont_annotation);
        }
    }

    // return true if the term could or should appear in the interesting_parent_details
    // field of the TermDetails and TermShort structs
    fn is_interesting_parent(&self, termid: &FlexStr, rel_name: &FlexStr) -> bool {
        self.possible_interesting_parents.contains(&InterestingParent {
            termid: termid.into(),
            rel_name: rel_name.into(),
        })
    }

    fn process_cvtermpath(&mut self) {
        let mut slim_termids = HashSet::new();
        for slim_config in self.config.slims.values() {
            for term_and_name in &slim_config.terms {
                slim_termids.insert(term_and_name.termid.clone());
            }
        }

        #[allow(clippy::type_complexity)]
        let mut new_annotations: HashMap<(CvName, TermId), HashMap<TermId, HashMap<i32, HashSet<RelName>>>> =
            HashMap::new();

        let mut children_by_termid: HashMap<TermId, HashSet<TermId>> = HashMap::new();

        for cvtermpath in &self.raw.cvtermpaths {
            let subject_term = &cvtermpath.subject;
            let subject_termid = subject_term.termid();
            let object_term = &cvtermpath.object;
            let object_termid = object_term.termid();

            if let Some(subject_term_details) = self.terms.get(&subject_termid) {
                let rel_termid =
                    match cvtermpath.rel_type {
                        Some(ref rel_type) => {
                            rel_type.termid()
                        },
                        None => panic!("no relation type for {} <-> {}\n",
                                       &subject_term.name, &object_term.name)
                    };

                let rel_term_name =
                    self.make_term_short(&rel_termid).name;

                if rel_term_name == "has_part" &&
                    !HAS_PART_CV_NAMES.contains(&subject_term_details.cv_name) {
                        continue;
                    }

                if !DESCENDANT_REL_NAMES.contains(&rel_term_name.as_str()) {
                    continue;
                }

                if subject_term_details.cv_annotations.keys().len() > 0 ||
                    slim_termids.contains(&object_termid)
                {
                    children_by_termid
                        .entry(object_termid.clone())
                        .or_default()
                        .insert(subject_termid.clone());
                }

                for (cv_name, term_annotations) in &subject_term_details.cv_annotations {
                    for term_annotation in term_annotations {
                        for annotation_id in &term_annotation.annotations {
                            let dest_termid = object_termid.clone();
                            let source_termid = subject_termid.clone();

                            if !term_annotation.is_not {
                                new_annotations.entry((cv_name.clone(), dest_termid))
                                    .or_default()
                                    .entry(source_termid)
                                    .or_default()
                                    .entry(*annotation_id)
                                    .or_default()
                                    .insert(rel_term_name.clone());
                            }
                        }
                    }
                }
            } else {
                panic!("TermDetails not found for {}", &subject_termid);
            }
        }

        for ((dest_cv_name, dest_termid), dest_annotations_map) in new_annotations.drain() {
            for (source_termid, source_annotations_map) in dest_annotations_map {
                let mut new_annotations: Vec<OntAnnotationId> = vec![];
                let mut all_rel_names: HashSet<FlexStr> = HashSet::new();
                for (annotation_id, rel_names) in source_annotations_map {
                    new_annotations.push(annotation_id);
                    for rel_name in rel_names {
                        all_rel_names.insert(rel_name);
                    }
                }

                let new_annotations =
                    self.make_term_annotations(&source_termid, &new_annotations, false);

                let dest_term_details = {
                    self.terms.get_mut(&dest_termid).unwrap()
                };

                for (_, new_annotation) in new_annotations {
                    let mut new_annotation_clone = new_annotation.clone();

                    new_annotation_clone.rel_names.extend(all_rel_names.clone());

                    dest_term_details.cv_annotations
                        .entry(dest_cv_name.clone())
                        .or_insert_with(Vec::new)
                        .push(new_annotation_clone);
                }
            }
        }


        let mut term_names = HashMap::new();
        for (termid, term_details) in &self.terms {
            term_names.insert(termid.clone(), term_details.name.to_lowercase());
        }

        for term_details in self.terms.values_mut() {
            let term_details_termid = &term_details.termid;
            for term_annotations in term_details.cv_annotations.values_mut() {
                let ont_term_cmp = |ont_term_1: &OntTermAnnotations, ont_term_2: &OntTermAnnotations| {
                    if ont_term_1.term == ont_term_2.term {
                        return Ordering::Equal;
                    }

                    // put direct annotation first on page
                    if ont_term_1.term == *term_details_termid {
                        return Ordering::Less;
                    }
                    if ont_term_2.term == *term_details_termid {
                        return Ordering::Greater;
                    }

                    if !ont_term_1.is_not && ont_term_2.is_not {
                        return Ordering::Less;
                    }
                    if ont_term_1.is_not && !ont_term_2.is_not {
                        return Ordering::Greater;
                    }
                    let term1 = &term_names[&ont_term_1.term];
                    let term2 = &term_names[&ont_term_2.term];

                    term1.cmp(term2)
                };

                term_annotations.sort_by(&ont_term_cmp);
            }
        }

        self.children_by_termid = children_by_termid;
    }

    fn make_metadata(&mut self) -> Metadata {
        let mut db_creation_datetime = None;
        let mut date_version = None;

        let mut data_source_versions = HashMap::new();

        for chadoprop in &self.raw.chadoprops {
            if chadoprop.prop_type.name == "db_creation_datetime" {
                db_creation_datetime = chadoprop.value.clone();
                continue;
            }
            if chadoprop.prop_type.name == "date_version" {
                date_version = chadoprop.value.clone();
                continue;
            }
            if chadoprop.prop_type.name == "db_date_version" {
                continue;
            }
            if chadoprop.prop_type.name.ends_with("_version")
                && let Some(ref value) = chadoprop.value {
                    let trimmed_type =
                        chadoprop.prop_type.name.trim_end_matches("_version").to_shared_str();
                    data_source_versions.insert(trimmed_type,
                                                value.to_owned());
                }
        }

        let mut cv_versions = HashMap::new();

        for cvprop in &self.raw.cvprops {
            if cvprop.prop_type.name == "cv_version" {
                cv_versions.insert(cvprop.cv.name.clone(), cvprop.value.clone());
            }
        }

        const PKG_NAME: &str = env!("CARGO_PKG_NAME");
        const VERSION: &str = env!("CARGO_PKG_VERSION");

        data_source_versions.insert(flex_str!("InterPro"),
                                    self.domain_data.interproscan_version.clone());

        Metadata {
            export_prog_name: flex_str!(PKG_NAME),
            export_prog_version: flex_str!(VERSION),
            db_creation_datetime: db_creation_datetime.unwrap(),
            date_version: date_version.unwrap(),
            data_source_versions,
            gene_count: self.genes.len(),
            term_count: self.terms.len(),
            cv_versions,
        }
    }

    pub fn get_api_genotype_annotation(&self) -> HashMap<TermId, Vec<APIGenotypeAnnotation>>
    {
        let mut app_genotype_annotation = HashMap::new();

        for term_details in self.terms.values() {
            for annotations_vec in term_details.cv_annotations.values() {
                for ont_term_annotations in annotations_vec {
                    'DETAILS: for annotation_id in &ont_term_annotations.annotations {
                        let annotation_details = self.annotation_details.
                            get(annotation_id).expect("can't find OntAnnotationDetail");

                        if annotation_details.genotype.is_none() {
                            continue 'DETAILS;
                        }
                        let genotype_uniquename = annotation_details.genotype.clone().unwrap();
                        let genotype =
                            &term_details.genotypes_by_uniquename[&genotype_uniquename];
                        let conditions =
                            annotation_details.conditions.iter()
                            .map(|cond_termid| {
                                let cond_term = self.terms.get(cond_termid).unwrap();
                                TermAndName {
                                    termid: cond_term.termid.clone(),
                                    name: cond_term.name.clone(),
                                }
                            })
                            .collect::<HashSet<_>>();
                        let mut api_annotation = APIGenotypeAnnotation {
                            is_multi: genotype.loci.len() > 1,
                            ploidiness: genotype.ploidiness(),
                            conditions,
                            alleles: vec![],
                        };
                        for locus in &genotype.loci {
                            for allele in &locus.expressed_alleles {
                                let allele_uniquename = &allele.allele_uniquename;
                                let allele_short =
                                    self.alleles.get(allele_uniquename).expect("Can't find allele");
                                let allele_gene_uniquename =
                                    allele_short.gene.uniquename.clone();
                                let allele_details = APIAlleleDetails {
                                    gene: allele_gene_uniquename,
                                    allele_type: allele_short.allele_type.clone(),
                                    expression: allele.expression.clone(),
                                };
                                api_annotation.alleles.push(allele_details);
                            }
                        }
                        app_genotype_annotation
                            .entry(term_details.termid.clone())
                            .or_insert_with(Vec::new)
                            .push(api_annotation);
                    }
                }
            }
        }

        app_genotype_annotation
    }

    fn make_protein_data(&self, gene_details: &GeneDetails)
                         -> (Option<f32>, Option<usize>, Option<GeneQueryAttrName>)
    {
        let mut molecular_weight = None;
        let mut protein_length = None;

        for transcript_uniquename in &gene_details.transcripts {
            if let Some(transcript) = self.transcripts.get(transcript_uniquename)
                && let Some(ref protein) = transcript.protein {
                    molecular_weight = Some((100.0 * protein.molecular_weight).round() / 100.0);
                    if protein.sequence.ends_with('*') {
                        protein_length = Some(protein.sequence.len() - 1);
                    } else {
                        protein_length = Some(protein.sequence.len());
                    }
                    break;
                }
        }

        for field_name in &self.config.gene_results.visualisation_field_names {
            let column_conf = &self.config.gene_results.field_config[field_name];
            for attr_value_conf in &column_conf.attr_values {
                if let (Some(ref bin_start), Some(ref bin_end)) =
                    (attr_value_conf.bin_start, attr_value_conf.bin_end)
                        && let Some(prot_len) = protein_length
                            && *bin_start <= prot_len && *bin_end >= prot_len {
                                return (molecular_weight,
                                        Some(prot_len),
                                        Some(attr_value_conf.name.clone()));
                            }
            }
        }

        (None, None, None)
    }

    fn make_gene_query_go_data(&self, gene_details: &GeneDetails, term_config: &[TermId],
                               cv_name: &FlexStr) -> Option<GeneQueryTermData>
    {
        let component_term_annotations =
            gene_details.cv_annotations.get(cv_name)?;

        let in_component = |check_termid: &FlexStr| {
            for term_annotation in component_term_annotations {
                let maybe_term_details = self.terms.get(&term_annotation.term);

                let term_details =
                    maybe_term_details .unwrap_or_else(|| {
                        panic!("can't find TermDetails for {}", &term_annotation.term)
                    });

                let interesting_parent_ids = &term_details.interesting_parent_ids;

                if !term_annotation.is_not &&
                    (term_annotation.term == *check_termid ||
                     interesting_parent_ids.contains(check_termid))
                {
                    return true;
                }
            }
            false
        };

        for go_component_termid in term_config {
            if in_component(go_component_termid) {
                return Some(GeneQueryTermData::Term(TermAndName {
                    termid: go_component_termid.to_owned(),
                    name: self.terms.get(go_component_termid).unwrap().name.clone(),
                }));
            }
        }

        Some(GeneQueryTermData::Other)
    }

    fn get_ortholog_taxonids(&self, gene_details: &GeneDetails)
                             -> HashSet<u32>
    {
        let mut return_set = HashSet::new();

        for ortholog_annotation in &gene_details.ortholog_annotations {
            return_set.insert(ortholog_annotation.ortholog_taxonid);
        }

        return_set
    }

    fn get_physical_interactors(&self, gene_details: &GeneDetails)
                                -> HashSet<GeneUniquename>
    {
        let mut return_set = HashSet::new();

        for physical_interaction in &gene_details.physical_interactions {
            if gene_details.uniquename == physical_interaction.gene_uniquename {
                return_set.insert(physical_interaction.interactor_uniquename.clone());
            } else {
                // gene is the prey for this interaction
                return_set.insert(physical_interaction.gene_uniquename.clone());
            }
        }

        return_set
    }

    fn get_rna_lengths(&self, gene_details: &GeneDetails) -> (Option<usize>, Option<usize>) {
        if let Some(transcript_uniquename) = gene_details.transcripts.first() {
            if let Some(transcript) =
                self.transcripts.get(transcript_uniquename) {
                    let spliced_length = transcript.rna_seq_length_spliced.map(|i| i.get());
                    let unspliced_length = transcript.rna_seq_length_unspliced.map(|i| i.get());

                    (spliced_length, unspliced_length)
                } else {
                    (None, None)
                }
        } else {
            (None, None)
        }
    }


    fn make_gene_query_data_map(&self) -> HashMap<GeneUniquename, GeneQueryData> {
        let mut gene_query_data_map = HashMap::new();

        for gene_details in self.genes.values() {
            if self.config.load_organism_taxonid.is_some() &&
                self.config.load_organism_taxonid.unwrap() != gene_details.taxonid
            {
                continue;
            }

            let gene_uniquename = &gene_details.uniquename;
            let ortholog_taxonids = self.get_ortholog_taxonids(gene_details);
            let physical_interactors = self.get_physical_interactors(gene_details);
            let reference_uniquenames =
                 gene_details.references_by_uniquename.keys().cloned().collect();

            let mut cc_terms = vec![];
            let mut process_terms = vec![];
            let mut function_terms = vec![];

            for field_name in &self.config.gene_results.visualisation_field_names {
                let column_conf = &self.config.gene_results.field_config[field_name];

                for attr_value_conf in &column_conf.attr_values {
                    if let Some(ref termid) = attr_value_conf.termid {
                        match field_name.as_ref() {
                            "go_component" => cc_terms.push(termid.clone()),
                            "go_process_superslim" => process_terms.push(termid.clone()),
                            "go_function" => function_terms.push(termid.clone()),
                            _ => (),
                        }
                    }
                }
            }


            let go_component =
                self.make_gene_query_go_data(gene_details, &cc_terms,
                    &flex_str!("cellular_component"));
            let go_process_superslim =
                self.make_gene_query_go_data(gene_details, &process_terms,
                    &flex_str!("biological_process"));
            let go_function =
                self.make_gene_query_go_data(gene_details, &function_terms,
                    &flex_str!("molecular_function"));

            let tmm =
                if gene_details.feature_type == "mRNA gene" {
                    if gene_details.tm_domain_coords.is_empty() {
                        Some(PresentAbsent::Absent)
                    } else {
                        Some(PresentAbsent::Present)
                    }
                } else {
                    Some(PresentAbsent::NotApplicable)
                };

            let (molecular_weight, protein_length, protein_length_bin) =
                self.make_protein_data(gene_details);

            let (spliced_rna_length, unspliced_rna_length) = self.get_rna_lengths(gene_details);
            let pdb_ids = gene_details.pdb_entries.iter()
                .map(|pdb_entry| pdb_entry.pdb_id.clone())
                .collect();
            let rnacentral_urs_identifier = gene_details.rnacentral_urs_identifier.clone();

            let mut gocam_ids = HashSet::new();
            let mut enables_gocam_activity_ids = HashSet::new();

            let uniquename_with_prefix =
                if gene_uniquename.contains(':') {
                    gene_uniquename.to_std_string()
                } else {
                    format!("{}:{}", self.config.database_name, gene_uniquename)
                };

            for model in &self.gocam_models {
                if model.genes_in_model().contains(&uniquename_with_prefix) {
                    gocam_ids.insert(model.id().into());
                }
                if model.genes_enabling_activities().contains_key(&uniquename_with_prefix) {
                    enables_gocam_activity_ids.insert(model.id().into());
                }
            }

            let paralogs =
                gene_details.paralog_annotations.iter()
                .map(|paralog_annotation| {
                    paralog_annotation.paralog_uniquename.clone()
                })
                .collect();

            let mut property_flags = HashSet::new();
            if gene_details.rnacentral_2d_structure_id.is_some() {
                property_flags.insert(GeneQueryPropFlag::Rnacentral2DStructure);
            }
            if !gene_details.paralog_annotations.is_empty() {
                property_flags.insert(GeneQueryPropFlag::HasParalog);
            }

            let gene_query_data = GeneQueryData {
                gene_uniquename: gene_details.uniquename.clone(),
                deletion_viability: gene_details.deletion_viability.clone(),
                go_component,
                go_process_superslim,
                go_function,
                characterisation_status: gene_details.characterisation_status.clone(),
                taxonomic_distribution: gene_details.taxonomic_distribution.clone(),
                tmm,
                ortholog_taxonids,
                physical_interactors,
                molecular_weight,
                protein_length,
                protein_length_bin,
                spliced_rna_length,
                unspliced_rna_length,
                reference_uniquenames,
                pdb_ids,
                rnacentral_urs_identifier,
                gocam_ids,
                enables_gocam_activity_ids,
                paralogs,
                subset_termids: gene_details.subset_termids.clone(),

                property_flags,
            };

            gene_query_data_map.insert(gene_details.uniquename.clone(), gene_query_data);
        }

        gene_query_data_map
    }


    fn make_secondary_identifiers_map(&self)
                                  -> HashMap<TermId, TermId>
    {
        let mut ret_map = HashMap::new();

        for term_details in self.terms.values() {
            for secondary_identifier in &term_details.secondary_identifiers {
                ret_map.insert(secondary_identifier.clone(),
                               term_details.termid.clone());
            }
        }

        ret_map
    }

    fn make_api_maps_downstream_genes(&self)
        -> HashMap<GeneUniquename, HashMap<TermId, HashSet<GeneUniquename>>>
    {
        let mut downstream_genes = HashMap::new();

        let Some(during_config) = self.config.extension_categories.get("during")
        else {
            return downstream_genes;
        };

        let mut possible_phases = HashSet::new();

        for during_conf in during_config {
            possible_phases.extend(during_conf.ancestors.iter());
        }

        for gene_details in self.genes.values() {
           for (cv_name, term_annotations) in gene_details.cv_annotations.iter() {
              let Some(cv_config) = self.config.cv_config.get(cv_name)
              else {
                continue;
              };

              let downstream_relations = &cv_config.downstream_relations;
              if downstream_relations.is_empty() {
                continue;
              }

              for term_annotation in term_annotations {
                 for annotation in &term_annotation.annotations {
                    let annotation_details = self.annotation_details.get(annotation).unwrap();
                    let mut maybe_downstream_gene: Option<GeneUniquename> = None;
                    let mut downstream_gene_phase: Option<TermId> = None;

                    for ext_part in &annotation_details.extension {
                        let Some(ref rel_type_id) = ext_part.rel_type_id
                        else {
                            continue;
                        };
                        if rel_type_id == "RO:0002092" {
                            if let ExtRange::Term(ref phase_termid) = ext_part.ext_range {
                                downstream_gene_phase = Some(phase_termid.into());
                            }
                        } else if (downstream_relations.contains(rel_type_id) ||
                        downstream_relations.contains(&ext_part.rel_type_name))
                         && let ExtRange::Gene(ref target_gene_uniquename) = ext_part.ext_range {
                             maybe_downstream_gene = Some(target_gene_uniquename.into());
                         }
                    }

                    if let Some(ref downstream_gene_uniquename) = maybe_downstream_gene {
                        let gene_uniquename = &gene_details.uniquename;

                        let mut downstream_gene_phase_and_parents = vec![];

                        if let Some(ref downstream_gene_phase) = downstream_gene_phase {
                            let phase_term_details = self.terms.get(downstream_gene_phase)
                                .unwrap_or_else(|| panic!("internal error: failed to find term {}", downstream_gene_phase));

                            for parent_id in phase_term_details.interesting_parent_ids.iter() {
                                if possible_phases.contains(parent_id) {
                                    downstream_gene_phase_and_parents.push(parent_id.to_owned());
                                }
                            }

                            downstream_gene_phase_and_parents.push(downstream_gene_phase.to_owned());
                        } else {
                            downstream_gene_phase_and_parents.push(flex_str!(""));
                        }

                        for downstream_gene_phase_key in downstream_gene_phase_and_parents {
                            let key = flex_fmt!("{}--{}", cv_name, gene_uniquename);
                            downstream_genes
                                .entry(key)
                                .or_insert_with(HashMap::new)
                                .entry(downstream_gene_phase_key.clone())
                                .or_insert_with(HashSet::new)
                                .insert(downstream_gene_uniquename.clone());
                        }
                    }
                 }
              }
            }
        }

        downstream_genes
    }

    fn make_gocam_data_by_gene(&self) -> HashMap<GeneUniquename, HashSet<GoCamId>> {
        let mut ret = HashMap::new();

        for (gocam_id, gocam_summary) in &self.gocam_summaries {
            for gene_uniquename in &gocam_summary.activity_enabling_genes {
                ret.entry(gene_uniquename.clone())
                   .or_insert_with(HashSet::new)
                   .insert(gocam_id.clone());
            }
        }

        ret
    }

    pub fn make_api_maps(mut self) -> APIMaps {
        let mut gene_summaries: HashMap<GeneUniquename, APIGeneSummary> = HashMap::new();
        let mut gene_name_gene_map = HashMap::new();
        let mut interactors_of_genes = HashMap::new();
        let downstream_genes = self.make_api_maps_downstream_genes();

        for (gene_uniquename, gene_details) in &self.genes {
            if self.config.load_organism_taxonid.is_none() ||
                self.config.load_organism_taxonid.unwrap() == gene_details.taxonid {
                let gene_summary = self.make_api_gene_summary(gene_uniquename);
                if let Some(ref gene_name) = gene_summary.name {
                    gene_name_gene_map.insert(gene_name.clone(), gene_uniquename.clone());
                }
                gene_summaries.insert(gene_uniquename.clone(), gene_summary);

                let mut interactors = HashSet::new();

                for interaction_annotation in &gene_details.physical_interactions {
                    let interactor_uniquename =
                        if gene_uniquename == interaction_annotation.gene_uniquename {
                            interaction_annotation.interactor_uniquename.clone()
                        } else {
                            interaction_annotation.gene_uniquename.clone()
                        };
                    let interactor = APIInteractor {
                        interaction_type: InteractionType::Physical,
                        interactor_uniquename,
                        throughput: interaction_annotation.throughput.clone(),
                        evidence_type: interaction_annotation.evidence.clone(),
                    };
                    interactors.insert(interactor);
                }
                for (interaction_key, annotations) in gene_details.genetic_interactions.iter() {
                    let interactor_uniquename =
                        if gene_uniquename == interaction_key.gene_a_uniquename {
                            interaction_key.gene_b_uniquename.clone()
                        } else {
                            interaction_key.gene_a_uniquename.clone()
                        };
                    for annotation in annotations {
                        let interactor = APIInteractor {
                            interaction_type: InteractionType::Genetic,
                            interactor_uniquename: interactor_uniquename.clone(),
                            throughput: annotation.throughput.clone(),
                            evidence_type: interaction_key.interaction_type.clone(),
                        };
                        interactors.insert(interactor);
                    }
                }
                interactors_of_genes.insert(gene_uniquename.clone(), interactors);
            }
        }

        let protein_view_data =
            make_protein_view_data_map(&self.genes,
                                       &self.terms,
                                       &self.annotation_details,
                                       &self.genotypes, &self.alleles,
                                       &self.transcripts, &self.references,
                                       self.config);

        let gocam_data_by_gene = self.make_gocam_data_by_gene();

        let gene_query_data_map = self.make_gene_query_data_map();

        let mut termid_genes: HashMap<TermId, HashSet<GeneUniquename>> = HashMap::new();

        for (termid, term_details) in self.terms.drain() {
            let cv_config = &self.config.cv_config;
            if let Some(term_config) = cv_config.get(&term_details.cv_name)
                && term_config.feature_type == "gene" {
                    termid_genes.insert(termid.clone(),
                                        term_details.annotated_genes.clone());
                }
        }

        let seq_feature_page_features: Vec<FeatureShort> =
            self.other_features.values()
            .filter(|feature_short| {
                let so_types_to_show =
                    &self.config.sequence_feature_page.so_types_to_show;
                let feature_type_string =
                    feature_short.feature_type.to_string().to_shared_str();
                so_types_to_show.contains(&feature_type_string)
            })
            .map(|feature_short| {
                let mut new_feature = feature_short.clone();
                // we don't need the residues for the seq feature page
                new_feature.residues = flex_str!("");
                new_feature
            }).collect();

        let secondary_identifiers_map = self.make_secondary_identifiers_map();

        // avoid clone()
        let mut term_subsets = HashMap::new();
        std::mem::swap(&mut term_subsets, &mut self.term_subsets);

        let mut gene_subsets = HashMap::new();
        std::mem::swap(&mut gene_subsets, &mut self.gene_subsets);

        let mut children_by_termid = HashMap::new();
        std::mem::swap(&mut children_by_termid, &mut self.children_by_termid);

        let mut gene_expression_measurements = HashMap::new();
        std::mem::swap(&mut gene_expression_measurements,
                       &mut self.gene_expression_measurements);

        APIMaps {
            gene_summaries,
            gene_query_data_map,
            termid_genes,
            gene_name_gene_map,
            transcripts: self.transcripts,
            interactors_of_genes,
            downstream_genes,
            other_features: self.other_features,
            seq_feature_page_features,
            chromosomes: self.chromosomes,
            term_subsets,
            gene_subsets,
            children_by_termid,
            gene_expression_measurements,
            secondary_identifiers_map,
            protein_view_data,
            gocam_data_by_gene,
            gocam_data_by_gocam_id: self.gocam_summaries,
            gocam_overlaps: self.gocam_overlaps,
            gocam_overlaps_merge_by_chemical: self.gocam_overlaps_merge_by_chemical,
            gocam_holes: self.gocam_holes,
            pro_term_to_gene_map: self.pro_term_to_gene,

            protein_complex_data: self.protein_complex_data,
            protein_complexes: self.protein_complexes,
       }
    }

    fn gene_uniquenames_from_genotype(&self, genotype_details: &GenotypeDetails) -> Vec<FlexStr> {
        let mut ret = vec![];

        for locus in &genotype_details.loci {
            for expressed_allele in &locus.expressed_alleles {
                let allele_uniquename = &expressed_allele.allele_uniquename;
                let gene_uniquename =
                    self.alleles.get(allele_uniquename)
                    .unwrap().gene.uniquename.clone();
                ret.push(gene_uniquename);
            }
        }

        ret
    }

    #[allow(clippy::too_many_arguments)]
    fn add_cv_annotations_to_maps(&self,
                                  identifier: &FlexStr,
                                  cv_annotations: &OntAnnotationMap,
                                  seen_references: &mut HashMap<FlexStr, ReferenceShortOptionMap>,
                                  seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                                  seen_genotypes: &mut HashMap<FlexStr, GenotypeShortMap>,
                                  seen_alleles: &mut HashMap<FlexStr, AlleleShortMap>,
                                  seen_transcripts: &mut HashMap<FlexStr, TranscriptDetailsOptionMap>,
                                  seen_terms: &mut HashMap<FlexStr, TermShortOptionMap>) {
        for feat_annotations in cv_annotations.values() {
            for feat_annotation in feat_annotations.iter() {
                self.add_term_to_hash(seen_terms, identifier,
                                      &feat_annotation.term);

                for annotation_detail_id in &feat_annotation.annotations {
                    let annotation_detail = self.annotation_details.
                        get(annotation_detail_id).expect("can't find OntAnnotationDetail");

                    for transcript_uniquename in &annotation_detail.transcript_uniquenames {
                        self.add_transcript_to_hashes(seen_transcripts,
                                                      seen_genes,
                                                      identifier,
                                                      transcript_uniquename);
                    }

                    self.add_ref_to_hash(seen_references, identifier,
                                         &annotation_detail.reference);
                    for condition_termid in &annotation_detail.conditions {
                        self.add_term_to_hash(seen_terms, identifier,
                                              condition_termid);
                    }
                    if let Some(ref gene_product_form_id) =
                        annotation_detail.gene_product_form_id
                            && gene_product_form_id.starts_with("PR:") {
                                self.add_term_to_hash(seen_terms, identifier,
                                                      gene_product_form_id);
                            }
                    self.add_extension_to_maps(&annotation_detail.extension, seen_genes,
                                               seen_transcripts, seen_terms, identifier);
                    if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  seen_references,
                                                  identifier, genotype_uniquename);
                    }

                    let with_from_iter = annotation_detail.withs
                        .iter()
                        .chain(annotation_detail.froms.iter());

                    for with_from_value in with_from_iter {
                        match with_from_value {
                            WithFromValue::Gene(gene_short) => {
                                self.add_gene_to_hash(seen_genes, identifier,
                                                      &gene_short.uniquename)
                            },
                            WithFromValue::Transcript(transcript_uniquename) => {
                                self.add_transcript_to_hashes(seen_transcripts, seen_genes,
                                                              identifier,  transcript_uniquename);
                            },
                            _ => (),
                        }
                    }
                }
            }
        }
    }

    fn add_extension_to_maps(&self, extension: &Vec<ExtPart>,
                             seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                             seen_transcripts: &mut HashMap<FlexStr, TranscriptDetailsOptionMap>,
                             seen_terms: &mut HashMap<FlexStr, TermShortOptionMap>,
                             map_key: &FlexStr) {
        for ext_part in extension {
            match ext_part.ext_range {
                ExtRange::Term(ref range_termid) |
                ExtRange::GeneProduct(ref range_termid) =>
                    self.add_term_to_hash(seen_terms, map_key,
                                          range_termid),
                ExtRange::GeneAndGeneProduct(GeneAndGeneProduct { ref gene_uniquename, ref product }) => {
                    self.add_gene_to_hash(seen_genes, map_key, gene_uniquename);
                    self.add_term_to_hash(seen_terms, map_key, product);
                },
                ExtRange::Gene(ref gene_uniquename) |
                ExtRange::Promoter(ref gene_uniquename) =>
                    self.add_gene_to_hash(seen_genes, map_key,
                                          gene_uniquename),
                ExtRange::Transcript(ref transcript_uniquename) =>
                    self.add_transcript_to_hashes(seen_transcripts, seen_genes,
                                                map_key,
                                                transcript_uniquename),
                _ => {},
            }
        }
    }

    fn set_term_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_transcripts, mut seen_terms) = get_maps();

        let mut annotated_genes_map: HashMap<TermId, HashSet<GeneUniquename>> =
            HashMap::new();
        let mut single_locus_annotated_genes_map: HashMap<TermId, HashSet<GeneUniquename>> =
            HashMap::new();
        let mut multi_locus_annotated_genes_map: HashMap<TermId, HashSet<GeneUniquename>> =
            HashMap::new();

        let set_interaction_maps = |genetic_interactions: &GeneticInteractionMap,
                                    termid: &FlexStr,
                                    seen_references: &mut HashMap<FlexStr, ReferenceShortOptionMap>,
                                    seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                                    seen_genotypes: &mut HashMap<FlexStr, GenotypeShortMap>,
                                    seen_alleles: &mut HashMap<FlexStr, AlleleShortMap>,
                                    seen_transcripts: &mut HashMap<FlexStr, TranscriptDetailsOptionMap>,
                                    seen_terms: &mut HashMap<GeneUniquename, TermShortOptionMap>| {

            let interaction_iter = genetic_interactions.iter();
            for (interaction_key, interaction_details) in interaction_iter {
                self.add_gene_to_hash(seen_genes, termid,
                                      &interaction_key.gene_a_uniquename);
                self.add_gene_to_hash(seen_genes, termid,
                                      &interaction_key.gene_b_uniquename);

                for interaction_detail in interaction_details {

                    self.add_ref_to_hash(seen_references, termid,
                                         &interaction_detail.reference_uniquename);

                    if let Some(ref genotype_a_uniquename) = interaction_detail.genotype_a_uniquename {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  seen_references,
                                                  termid,
                                                  genotype_a_uniquename);
                    }
                    if let Some(ref genotype_b_uniquename) = interaction_detail.genotype_b_uniquename {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  seen_references,
                                                  termid,
                                                  genotype_b_uniquename);
                    }
                    if let Some(ref double_mutant_genotype_display_uniquename) = interaction_detail.double_mutant_genotype_display_uniquename {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  seen_references,
                                                  termid,
                                                  double_mutant_genotype_display_uniquename);
                    }

                    if let Some(ref double_mutant_phenotype_termid) = interaction_detail.double_mutant_phenotype_termid {
                        self.add_term_to_hash(seen_terms, termid,
                                              double_mutant_phenotype_termid);
                    }
                    if let Some(ref rescued_phenotype_termid) = interaction_detail.rescued_phenotype_termid {
                        self.add_term_to_hash(seen_terms, termid,
                                              rescued_phenotype_termid);
                    }

                    self.add_extension_to_maps(&interaction_detail.double_mutant_extension,
                                               seen_genes, seen_transcripts, seen_terms,
                                               termid);

                    self.add_extension_to_maps(&interaction_detail.rescued_phenotype_extension,
                                               seen_genes, seen_transcripts, seen_terms,
                                               termid);
                }
            }
        };

        for (termid, term_details) in &self.terms {
            for xref in &term_details.definition_xrefs {
                if xref.starts_with("PMID:") && self.references.contains_key(xref) {
                    self.add_ref_to_hash(&mut seen_references, termid, &Some(xref.clone()));
                }
            }

            for (cv_name, term_annotations) in &term_details.cv_annotations {
                for term_annotation in term_annotations {
                    self.add_term_to_hash(&mut seen_terms, termid,
                                          &term_annotation.term);
                    for annotation_detail_id in &term_annotation.annotations {
                        let annotation_detail = self.annotation_details
                            .get(annotation_detail_id).expect("can't find OntAnnotationDetail");

                        // prevent extension annotations from appearing
                        // in the normal query builder searches
                        if !cv_name.starts_with("extension:") &&
                            annotation_detail.genotype.is_none() {
                                for gene_uniquename in &annotation_detail.genes {
                                    self.add_gene_to_hash(&mut seen_genes, termid,
                                                          gene_uniquename);

                                    if !term_annotation.is_not {
                                    // prevent NOT annotation from appearing in the
                                    // counts on term pages and in the query builder
                                    annotated_genes_map
                                        .entry(termid.clone()).or_default()
                                        .insert(gene_uniquename.clone());
                                    }
                                }
                            }

                        if let Some(genotype_uniquename) = annotation_detail.genotype.as_ref() {
                            let genotype_details = self.genotypes.get(genotype_uniquename).unwrap();
                            let gene_uniquenames = self.gene_uniquenames_from_genotype(genotype_details);
                            if genotype_details.loci.len() == 1 {
                                single_locus_annotated_genes_map
                                    .entry(termid.clone()).or_default()
                                    .extend(gene_uniquenames.iter().cloned());
                            } else {
                                multi_locus_annotated_genes_map
                                    .entry(termid.clone()).or_default()
                                    .extend(gene_uniquenames.iter().cloned());
                            }
                        }

                        for transcript_uniquename in &annotation_detail.transcript_uniquenames {
                            self.add_transcript_to_hashes(&mut seen_transcripts,
                                                          &mut seen_genes,
                                                          termid, transcript_uniquename);
                        }

                        self.add_ref_to_hash(&mut seen_references, termid,
                                             &annotation_detail.reference);
                        for condition_termid in &annotation_detail.conditions {
                            self.add_term_to_hash(&mut seen_terms, termid,
                                                  condition_termid);
                        }
                        if let Some(ref gene_product_form_id) =
                            annotation_detail.gene_product_form_id
                                && gene_product_form_id.starts_with("PR:") {
                                    self.add_term_to_hash(&mut seen_terms, termid,
                                                          gene_product_form_id);
                                }
                        for ext_part in &annotation_detail.extension {
                            match ext_part.ext_range {
                                ExtRange::Term(ref range_termid) |
                                ExtRange::GeneProduct(ref range_termid) =>
                                    self.add_term_to_hash(&mut seen_terms, termid,
                                                          range_termid),
                                ExtRange::GeneAndGeneProduct(GeneAndGeneProduct { ref gene_uniquename, ref product }) => {
                                    self.add_gene_to_hash(&mut seen_genes, termid, gene_uniquename);
                                    self.add_term_to_hash(&mut seen_terms, termid, product);
                                },
                                ExtRange::Gene(ref gene_uniquename) |
                                ExtRange::Promoter(ref gene_uniquename) =>
                                    self.add_gene_to_hash(&mut seen_genes, termid,
                                                          gene_uniquename),
                                ExtRange::Transcript(ref transcript_uniquename) =>
                                    self.add_transcript_to_hashes(&mut seen_transcripts,
                                                                &mut seen_genes,
                                                                termid, transcript_uniquename),
                                _ => {},
                            }
                        }
                        if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                            self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                      &mut seen_genes,
                                                      &mut seen_references, termid,
                                                      genotype_uniquename);
                        }

                        let with_from_iter = annotation_detail.withs
                            .iter()
                            .chain(annotation_detail.froms.iter());

                        for with_from_value in with_from_iter {
                            match with_from_value {
                                WithFromValue::Gene(gene_short) => {
                                    self.add_gene_to_hash(&mut seen_genes, termid,
                                                          &gene_short.uniquename)
                                },
                                WithFromValue::Transcript(transcript_uniquename) => {
                                    self.add_transcript_to_hashes(&mut seen_transcripts, &mut seen_genes,
                                                                  termid,  transcript_uniquename);
                                },
                                _ => (),
                            }
                        }
                    }
                }
            }


            set_interaction_maps(&term_details.double_mutant_genetic_interactions,
                                 termid,
                                 &mut seen_references,
                                 &mut seen_genes,
                                 &mut seen_genotypes,
                                 &mut seen_alleles,
                                 &mut seen_transcripts,
                                 &mut seen_terms);

            set_interaction_maps(&term_details.single_allele_genetic_interactions,
                                 termid,
                                 &mut seen_references,
                                 &mut seen_genes,
                                 &mut seen_genotypes,
                                 &mut seen_alleles,
                                 &mut seen_transcripts,
                                 &mut seen_terms);
        }

        for (termid, term_details) in &mut self.terms {
            if let Some(genes) = seen_genes.remove(termid) {
                term_details.genes_by_uniquename = genes;
            }
            if let Some(genotypes) = seen_genotypes.remove(termid) {
                term_details.genotypes_by_uniquename = genotypes;
            }
            if let Some(alleles) = seen_alleles.remove(termid) {
                term_details.alleles_by_uniquename = alleles;
            }
            if let Some(references) = seen_references.remove(termid) {
                term_details.references_by_uniquename = references;
            }
            if let Some(transcripts) = seen_transcripts.remove(termid) {
                term_details.transcripts_by_uniquename = transcripts;
            }
            if let Some(terms) = seen_terms.remove(termid) {
                term_details.terms_by_termid = terms;
            }
            if let Some(gene_uniquename_set) = annotated_genes_map.remove(termid) {
                term_details.annotated_genes = gene_uniquename_set;
            }
            if let Some(gene_uniquename_set) = single_locus_annotated_genes_map.remove(termid) {
                term_details.single_locus_annotated_genes = gene_uniquename_set;
            }
            if let Some(gene_uniquename_set) = multi_locus_annotated_genes_map.remove(termid) {
                term_details.multi_locus_annotated_genes = gene_uniquename_set;
            }
        }
    }

    fn set_gene_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_transcripts, mut seen_terms) = get_maps();

        {
            for (gene_uniquename, gene_details) in &self.genes {
                self.add_cv_annotations_to_maps(gene_uniquename,
                                                &gene_details.cv_annotations,
                                                &mut seen_references,
                                                &mut seen_genes,
                                                &mut seen_genotypes,
                                                &mut seen_alleles,
                                                &mut seen_transcripts,
                                                &mut seen_terms);

                for transcript_uniquename in &gene_details.transcripts {
                    self.add_transcript_to_hashes(&mut seen_transcripts,
                                                &mut seen_genes,
                                                gene_uniquename,
                                                transcript_uniquename);
                }

                let interaction_iter = gene_details.physical_interactions.iter();
                for interaction in interaction_iter {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &interaction.reference_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &interaction.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &interaction.interactor_uniquename);
                }


                let interaction_iter = gene_details.genetic_interactions.iter();

                for (interaction_key, interaction_details) in interaction_iter {
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &interaction_key.gene_a_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &interaction_key.gene_b_uniquename);

                    for interaction_detail in interaction_details {

                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &interaction_detail.reference_uniquename);

                    if let Some(ref genotype_a_uniquename) = interaction_detail.genotype_a_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                  &mut seen_genes,
                                                  &mut seen_references,
                                                  gene_uniquename,
                                                  genotype_a_uniquename);
                    }
                    if let Some(ref genotype_b_uniquename) = interaction_detail.genotype_b_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                  &mut seen_genes,
                                                  &mut seen_references,
                                                  gene_uniquename,
                                                  genotype_b_uniquename);
                    }
                    if let Some(ref double_mutant_genotype_display_uniquename) = interaction_detail.double_mutant_genotype_display_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                  &mut seen_genes,
                                                  &mut seen_references,
                                                  gene_uniquename,
                                                  double_mutant_genotype_display_uniquename);
                    }

                    if let Some(ref double_mutant_phenotype_termid) = interaction_detail.double_mutant_phenotype_termid {
                        self.add_term_to_hash(&mut seen_terms, gene_uniquename,
                                              double_mutant_phenotype_termid);
                    }
                    if let Some(ref rescued_phenotype_termid) = interaction_detail.rescued_phenotype_termid {
                        self.add_term_to_hash(&mut seen_terms, gene_uniquename,
                                              rescued_phenotype_termid);
                    }
                    self.add_extension_to_maps(&interaction_detail.rescued_phenotype_extension,
                                               &mut seen_genes, &mut seen_transcripts, &mut seen_terms,
                                               gene_uniquename)

                    }
                }

                for ortholog_annotation in &gene_details.ortholog_annotations {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &ortholog_annotation.reference_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &ortholog_annotation.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &ortholog_annotation.ortholog_uniquename);
                }
                for paralog_annotation in &gene_details.paralog_annotations {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &paralog_annotation.reference_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &paralog_annotation.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &paralog_annotation.paralog_uniquename);
                }
                for target_of_annotation in &gene_details.target_of_annotations {
                    let target_of_gene = &target_of_annotation.gene;
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename, target_of_gene);
                    if let Some(ref annotation_genotype_uniquename) = target_of_annotation.genotype_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                  &mut seen_genes,
                                                  &mut seen_references,
                                                  gene_uniquename,
                                                  annotation_genotype_uniquename)
                    }
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &target_of_annotation.reference_uniquename);
                }

                for pub_and_source in &gene_details.feature_publications {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &Some(pub_and_source.reference_uniquename.clone()));
                }

                for pdb_ref_entry in &gene_details.pdb_entries {
                    if let Some(ref reference_uniquename) = pdb_ref_entry.reference_uniquename {
                         for chain in &pdb_ref_entry.gene_chains {
                           self.add_gene_to_hash(&mut seen_genes, &gene_details.uniquename,
                                                 &chain.gene_uniquename);
                         }

                        if !self.references.contains_key(reference_uniquename) {
                          // we don't have publication pages for all PDB references
                          continue;
                      }

                      self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                           &Some(reference_uniquename.clone()));
                      }
                }
            }
        }

        for (gene_uniquename, gene_details) in &mut self.genes {
            if let Some(references) = seen_references.remove(gene_uniquename) {
                gene_details.references_by_uniquename = references;
            }
            if let Some(alleles) = seen_alleles.remove(gene_uniquename) {
                gene_details.alleles_by_uniquename = alleles;
            }
            if let Some(genes) = seen_genes.remove(gene_uniquename) {
                gene_details.genes_by_uniquename = genes;
            }
            if let Some(genotypes) = seen_genotypes.remove(gene_uniquename) {
                gene_details.genotypes_by_uniquename = genotypes;
            }
            if let Some(transcripts) = seen_transcripts.remove(gene_uniquename) {
                gene_details.transcripts_by_uniquename = transcripts;
            }
            if let Some(terms) = seen_terms.remove(gene_uniquename) {
                gene_details.terms_by_termid = terms;
            }
        }
    }

    fn set_allele_details_maps(&mut self) {
        let mut seen_alleles = HashMap::new();
        let mut seen_genes = HashMap::new();
        let mut seen_references = HashMap::new();

        for (allele_uniquename, allele_details) in &self.alleles {
            for genotype_short in &allele_details.genotypes {
                for locus in &genotype_short.loci {
                    for expressed_allele in &locus.expressed_alleles {
                        let mut flags = HashSet::new();
                        flags.insert(AddToHashFlag::AlleleCommentRefs);
                        flags.insert(AddToHashFlag::AlleleSynonymRefs);
                        self.add_allele_to_hash(flags, &mut seen_alleles, &mut seen_genes,
                                                &mut seen_references, allele_uniquename,
                                                &expressed_allele.allele_uniquename);
                    }
                }
            }
        }

        for (allele_uniquename, allele_details) in &mut self.alleles {
            if let Some(references) = seen_references.remove(allele_uniquename) {
                allele_details.references_by_uniquename = references;
            }
            if let Some(alleles) = seen_alleles.remove(allele_uniquename) {
                allele_details.alleles_by_uniquename = alleles;
            }
            if let Some(genes) = seen_genes.remove(allele_uniquename) {
                allele_details.genes_by_uniquename = genes;
            }
        }
    }

    fn set_genotype_details_maps(&mut self) {
        let set_interaction_maps = |genetic_interactions: &GeneticInteractionMap,
                                    genotype_uniquename: &FlexStr,
                                    seen_references: &mut HashMap<FlexStr, ReferenceShortOptionMap>,
                                    seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                                    seen_genotypes: &mut HashMap<FlexStr, GenotypeShortMap>,
                                    seen_alleles: &mut HashMap<FlexStr, AlleleShortMap>,
                                    seen_transcripts: &mut HashMap<FlexStr, TranscriptDetailsOptionMap>,
                                    seen_terms: &mut HashMap<GeneUniquename, TermShortOptionMap>| {

            let interaction_iter = genetic_interactions.iter();

            for (interaction_key, interaction_details) in interaction_iter {
                self.add_gene_to_hash(seen_genes, genotype_uniquename,
                                      &interaction_key.gene_a_uniquename);
                self.add_gene_to_hash(seen_genes, genotype_uniquename,
                                      &interaction_key.gene_b_uniquename);

                for interaction_detail in interaction_details {

                    self.add_ref_to_hash(seen_references, genotype_uniquename,
                                         &interaction_detail.reference_uniquename);

                    if let Some(ref genotype_a_uniquename) = interaction_detail.genotype_a_uniquename {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  seen_references,
                                                  genotype_uniquename,
                                                  genotype_a_uniquename);
                    }
                    if let Some(ref genotype_b_uniquename) = interaction_detail.genotype_b_uniquename {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  seen_references,
                                                  genotype_uniquename,
                                                  genotype_b_uniquename);
                    }
                    if let Some(ref double_mutant_genotype_display_uniquename) = interaction_detail.double_mutant_genotype_display_uniquename {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  seen_references,
                                                  genotype_uniquename,
                                                  double_mutant_genotype_display_uniquename);
                    }

                    if let Some(ref double_mutant_phenotype_termid) = interaction_detail.double_mutant_phenotype_termid {
                        self.add_term_to_hash(seen_terms, genotype_uniquename,
                                              double_mutant_phenotype_termid);
                    }
                    if let Some(ref rescued_phenotype_termid) = interaction_detail.rescued_phenotype_termid {
                        self.add_term_to_hash(seen_terms, genotype_uniquename,
                                              rescued_phenotype_termid);
                    }
                    self.add_extension_to_maps(&interaction_detail.rescued_phenotype_extension,
                                               seen_genes, seen_transcripts, seen_terms,
                                               genotype_uniquename)
                }
            }
        };

        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_transcripts, mut seen_terms) = get_maps();


        for (genotype_uniquename, genotype_details) in &self.genotypes {
            self.add_cv_annotations_to_maps(genotype_uniquename,
                                            &genotype_details.cv_annotations,
                                            &mut seen_references,
                                            &mut seen_genes,
                                            &mut seen_genotypes,
                                            &mut seen_alleles,
                                            &mut seen_transcripts,
                                            &mut seen_terms);

            set_interaction_maps(&genotype_details.double_mutant_genetic_interactions,
                                 genotype_uniquename,
                                 &mut seen_references,
                                 &mut seen_genes,
                                 &mut seen_genotypes,
                                 &mut seen_alleles,
                                 &mut seen_transcripts,
                                 &mut seen_terms);

            set_interaction_maps(&genotype_details.rescue_genetic_interactions,
                                 genotype_uniquename,
                                 &mut seen_references,
                                 &mut seen_genes,
                                 &mut seen_genotypes,
                                 &mut seen_alleles,
                                 &mut seen_transcripts,
                                 &mut seen_terms);
        }

        for (genotype_uniquename, genotype_details) in &mut self.genotypes {
            if let Some(references) = seen_references.remove(genotype_uniquename) {
                genotype_details.references_by_uniquename = references;
            }
            if let Some(alleles) = seen_alleles.remove(genotype_uniquename) {
                genotype_details.alleles_by_uniquename = alleles;
            }
            if let Some(genes) = seen_genes.remove(genotype_uniquename) {
                genotype_details.genes_by_uniquename = genes;
            }
            if let Some(genotypes) = seen_genotypes.remove(genotype_uniquename) {
                genotype_details.genotypes_by_uniquename = genotypes;
            }
            if let Some(transcripts) = seen_transcripts.remove(genotype_uniquename) {
                genotype_details.transcripts_by_uniquename = transcripts;
            }
            if let Some(terms) = seen_terms.remove(genotype_uniquename) {
                genotype_details.terms_by_termid = terms;
            }
        }
    }

    fn set_reference_details_maps(&mut self) {
        // for calculating the gene_count field, we don't incude non-pombe genes
        let mut gene_count_hash: HashMap<FlexStr, GeneShortOptionMap> =
            HashMap::new();

        // counts of directly annotated genes from LTP annotations
        let mut ltp_gene_count_hash: HashMap<FlexStr, GeneShortOptionMap> =
            HashMap::new();

        let mut maybe_add_to_gene_count_hash =
            |reference_uniquename: &FlexStr, gene_uniquename: &GeneUniquename, ltp_only: bool| {
                let Some(load_org_taxonid) = self.config.load_organism_taxonid
                else {
                    return;
                };

                let Some(gene_details) = self.genes.get(gene_uniquename)
                else {
                    return;
                };

                if gene_details.taxonid == load_org_taxonid {
                    let hash = if ltp_only {
                        &mut ltp_gene_count_hash
                    } else {
                        &mut gene_count_hash
                    };
                    self.add_gene_to_hash(hash,
                                          reference_uniquename,
                                          gene_uniquename);
                }
            };


        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_transcripts, mut seen_terms) = get_maps();

        {
            for (reference_uniquename, reference_details) in &self.references {
                for feat_annotations in reference_details.cv_annotations.values() {
                    for feat_annotation in feat_annotations.iter() {
                        self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                              &feat_annotation.term);

                       for annotation_detail_id in &feat_annotation.annotations {
                            let annotation_detail = self.annotation_details
                                .get(annotation_detail_id).expect("can't find OntAnnotationDetail");

                           for transcript_uniquename in &annotation_detail.transcript_uniquenames {
                               self.add_transcript_to_hashes(&mut seen_transcripts,
                                                             &mut seen_genes,
                                                             reference_uniquename, transcript_uniquename);
                           }

                            for gene_uniquename in &annotation_detail.genes {
                                self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                                      gene_uniquename);
                                maybe_add_to_gene_count_hash(reference_uniquename,
                                                             gene_uniquename, false);
                                if annotation_detail.throughput == Some(Throughput::LowThroughput) {
                                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                             gene_uniquename, true);
                                }
                            }
                            for condition_termid in &annotation_detail.conditions {
                                self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                                      condition_termid);
                            }
                           if let Some(ref gene_product_form_id) =
                               annotation_detail.gene_product_form_id
                                   && gene_product_form_id.starts_with("PR:") {
                                       self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                                             gene_product_form_id);
                                   }
                            for ext_part in &annotation_detail.extension {
                                match ext_part.ext_range {
                                    ExtRange::Term(ref range_termid) |
                                    ExtRange::GeneProduct(ref range_termid) =>
                                        self.add_term_to_hash(&mut seen_terms,
                                                              reference_uniquename,
                                                              range_termid),
                                    ExtRange::GeneAndGeneProduct(GeneAndGeneProduct { ref gene_uniquename, ref product }) => {
                                        self.add_gene_to_hash(&mut seen_genes, reference_uniquename, gene_uniquename);
                                        self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                                              product);
                                    },
                                    ExtRange::Gene(ref gene_uniquename) |
                                    ExtRange::Promoter(ref gene_uniquename) => {
                                        self.add_gene_to_hash(&mut seen_genes,
                                                              reference_uniquename,
                                                              gene_uniquename);
                                        maybe_add_to_gene_count_hash(reference_uniquename,
                                                                     gene_uniquename, false);
                                    },
                                    ExtRange::Transcript(ref transcript_uniquename) =>
                                        self.add_transcript_to_hashes(&mut seen_transcripts,
                                                                      &mut seen_genes,
                                                                      reference_uniquename,
                                                                      transcript_uniquename),
                                    _ => {},
                                }
                            }

                            let with_from_iter = annotation_detail.withs
                                .iter()
                                .chain(annotation_detail.froms.iter());

                            for with_from_value in with_from_iter {
                                match with_from_value {
                                    WithFromValue::Gene(gene_short) => {
                                        self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                                              &gene_short.uniquename);
                                        maybe_add_to_gene_count_hash(reference_uniquename,
                                                                     &gene_short.uniquename, false);
                                    },
                                    WithFromValue::Transcript(transcript_uniquename) => {
                                        self.add_transcript_to_hashes(&mut seen_transcripts, &mut seen_genes,
                                                                      reference_uniquename,
                                                                      transcript_uniquename);
                                    },
                                    _ => (),
                                }
                            }

                            if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                                let genotype = self.make_genotype_short(genotype_uniquename);
                                self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                          &mut seen_genes,
                                                          &mut seen_references,
                                                          reference_uniquename,
                                                          &genotype.display_uniquename);
                            }
                        }
                    }
                }

                let interaction_iter =
                    reference_details.physical_interactions.iter();
                for interaction in interaction_iter {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction.gene_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &interaction.gene_uniquename, false);
                    if interaction.throughput == Some(Throughput::LowThroughput) {
                        maybe_add_to_gene_count_hash(reference_uniquename,
                                                     &interaction.gene_uniquename, true);
                    }
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction.interactor_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &interaction.interactor_uniquename, false);
                    if interaction.throughput == Some(Throughput::LowThroughput) {
                        maybe_add_to_gene_count_hash(reference_uniquename,
                                                     &interaction.interactor_uniquename, true);
                    }
                }

                let interaction_iter =
                    reference_details.genetic_interactions.iter();
                for (interaction_key, interaction_details) in interaction_iter {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction_key.gene_a_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &interaction_key.gene_a_uniquename, false);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction_key.gene_b_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &interaction_key.gene_b_uniquename, false);

                    for interaction_detail in interaction_details {
                    if let Some(ref genotype_a_uniquename) = interaction_detail.genotype_a_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                  &mut seen_genes,
                                                  &mut seen_references,
                                                  reference_uniquename,
                                                  genotype_a_uniquename);
                    }
                    if let Some(ref genotype_b_uniquename) = interaction_detail.genotype_b_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                  &mut seen_genes,
                                                  &mut seen_references,
                                                  reference_uniquename,
                                                  genotype_b_uniquename);
                    }
                    if let Some(ref double_mutant_genotype_display_uniquename) = interaction_detail.double_mutant_genotype_display_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                  &mut seen_genes,
                                                  &mut seen_references,
                                                  reference_uniquename,
                                                  double_mutant_genotype_display_uniquename);
                    }

                    if let Some(ref double_mutant_phenotype_termid) = interaction_detail.double_mutant_phenotype_termid {
                        self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                              double_mutant_phenotype_termid);
                    }
                    if let Some(ref rescued_phenotype_termid) = interaction_detail.rescued_phenotype_termid {
                        self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                              rescued_phenotype_termid);
                    }

                        if interaction_detail.throughput == Some(Throughput::LowThroughput) {
                            maybe_add_to_gene_count_hash(reference_uniquename,
                                                         &interaction_key.gene_a_uniquename, true);
                        }
                        if interaction_detail.throughput == Some(Throughput::LowThroughput) {
                            maybe_add_to_gene_count_hash(reference_uniquename,
                                                         &interaction_key.gene_b_uniquename, true);
                        }
                    }
                }

                for ortholog_annotation in &reference_details.ortholog_annotations {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &ortholog_annotation.gene_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &ortholog_annotation.gene_uniquename, false);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &ortholog_annotation.ortholog_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &ortholog_annotation.ortholog_uniquename, false);
                }
                for paralog_annotation in &reference_details.paralog_annotations {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &paralog_annotation.gene_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &paralog_annotation.gene_uniquename, false);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &paralog_annotation.paralog_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &paralog_annotation.paralog_uniquename, false);
                }

                for pdb_ref_entry in &reference_details.pdb_entries {
                    for pdb_gene_chain in &pdb_ref_entry.gene_chains {
                        let chain_gene_uniquename = &pdb_gene_chain.gene_uniquename;
                        self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                              chain_gene_uniquename);
                        maybe_add_to_gene_count_hash(&reference_details.uniquename,
                                                     chain_gene_uniquename, false);
                        maybe_add_to_gene_count_hash(&reference_details.uniquename,
                                                     chain_gene_uniquename, true);
                    }
                }

                for pubmed_keyword_gene in &reference_details.pubmed_keyword_genes {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          pubmed_keyword_gene);
                    maybe_add_to_gene_count_hash(&reference_details.uniquename,
                                                 pubmed_keyword_gene, false);
                }

                for extra_gene in &reference_details.extra_genes {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          extra_gene);
                    maybe_add_to_gene_count_hash(&reference_details.uniquename,
                                                 extra_gene, false);
                }
            }
        }

        for (reference_uniquename, reference_details) in &mut self.references {
            if let Some(genes) = seen_genes.remove(reference_uniquename) {
                reference_details.genes_by_uniquename = genes;
            }
            if let Some(genotypes) = seen_genotypes.remove(reference_uniquename) {
                reference_details.genotypes_by_uniquename = genotypes;
            }
            if let Some(alleles) = seen_alleles.remove(reference_uniquename) {
                reference_details.alleles_by_uniquename = alleles;
            }
            if let Some(terms) = seen_terms.remove(reference_uniquename) {
                reference_details.terms_by_termid = terms;
            }
            if let Some(references) = seen_references.remove(reference_uniquename) {
                reference_details.references_by_uniquename = references;
            }
            if let Some(transcripts) = seen_transcripts.remove(reference_uniquename) {
                reference_details.transcripts_by_uniquename = transcripts;
            }
            if let Some(gene_count_genes) = gene_count_hash.remove(reference_uniquename) {
                reference_details.gene_count = gene_count_genes.len();
            }
            if let Some(ltp_gene_count_genes) = ltp_gene_count_hash.remove(reference_uniquename) {
                reference_details.ltp_gene_count = ltp_gene_count_genes.len();
            }
        }
    }

    pub fn set_counts(&mut self) {
        let mut term_seen_genes: HashMap<TermId, HashSet<GeneUniquename>> = HashMap::new();
        let mut term_seen_genotypes: HashMap<TermId, HashSet<GenotypeUniquename>> = HashMap::new();
        let mut term_seen_single_locus_genotypes: HashMap<TermId, HashSet<GenotypeUniquename>> = HashMap::new();

        for (termid, term_details) in &self.terms {
            let mut seen_genes: HashSet<GeneUniquename> = HashSet::new();
            let mut seen_genotypes: HashSet<GenotypeUniquename> = HashSet::new();
            let mut seen_single_locus_genotypes: HashSet<GenotypeUniquename> = HashSet::new();
            for term_annotations in term_details.cv_annotations.values() {
                for term_annotation in term_annotations {
                    for annotation_detail_id in &term_annotation.annotations {
                        let annotation_detail = self.annotation_details
                            .get(annotation_detail_id).expect("can't find OntAnnotationDetail");
                        for gene_uniquename in &annotation_detail.genes {
                            if !term_annotation.is_not {
                                seen_genes.insert(gene_uniquename.clone());
                            }
                        }
                        if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                            seen_genotypes.insert(genotype_uniquename.clone());
                            let genotype = &self.genotypes[genotype_uniquename];
                            if genotype.loci.len() == 1 {
                                seen_single_locus_genotypes.insert(genotype_uniquename.clone());
                            }
                        }
                    }
                }
            }
            term_seen_genes.insert(termid.clone(), seen_genes);
            term_seen_genotypes.insert(termid.clone(), seen_genotypes);
            term_seen_single_locus_genotypes.insert(termid.clone(), seen_single_locus_genotypes);
        }

        let mut all_published_uniquenames = vec![];

        for (reference_uniquename, reference_details) in &self.references {
            if reference_details.pubmed_publication_date.is_some() {
                all_published_uniquenames.push(reference_uniquename.clone());
            }
        }

        let (recent_admin_curated, recent_community_curated,
             all_community_curated, all_admin_curated) =
            make_canto_curated(&self.references, &all_published_uniquenames);

        let recent_references = RecentReferences {
            pubmed: make_recently_added(&self.references, &all_published_uniquenames),
            admin_curated: recent_admin_curated,
            community_curated: recent_community_curated,
        };

        self.recent_references = recent_references;

        self.all_community_curated = all_community_curated;
        self.all_admin_curated = all_admin_curated;


        for term_details in self.terms.values_mut() {
            term_details.single_locus_genotype_uniquenames =
                term_seen_single_locus_genotypes.remove(&term_details.termid).unwrap();

            term_details.gene_count =
                term_seen_genes[&term_details.termid].len();
            term_details.genotype_count =
                term_seen_genotypes[&term_details.termid].len();
        }

        for genotype in self.genotypes.values_mut() {
            let mut annotation_count = 0;

            for term_annotations in genotype.cv_annotations().values() {
                for term_annotation in term_annotations {
                    annotation_count += term_annotation.annotations.len()
                }
            }

            genotype.annotation_count = annotation_count;
        }
    }

    fn set_split_by_parent_sets(&mut self) {
        for gene_details in self.genes.values_mut() {
            for (cv_name, term_annotations) in gene_details.cv_annotations.iter() {
                let Some(cv_config) = self.config.cv_config_by_name(cv_name)
                else {
                    continue;
                };

                let split_by_parents = &cv_config.split_by_parents;
                let mut split_by_parent_groups = HashMap::new();

                for term_annotation in term_annotations {
                    let term_details = self.terms.get(&term_annotation.term).unwrap();
                    let interesting_parent_ids = &term_details.interesting_parent_ids;

                    for split_by_config in split_by_parents {
                        for split_by_termid in &split_by_config.termids {
                            let (split_by_termid, not_flag) =
                                if let Some(without_prefex) = split_by_termid.strip_prefix("NOT ") {
                                    (without_prefex, true)
                                } else {
                                    (split_by_termid.as_str(), false)
                                };

                            let is_in_this_split =
                                term_annotation.term == split_by_termid ||
                                interesting_parent_ids.contains(split_by_termid);

                            if not_flag && !is_in_this_split || !not_flag && is_in_this_split {
                                split_by_parent_groups
                                    .entry(split_by_config.config_name.clone())
                                    .or_insert_with(HashSet::new)
                                    .insert(term_annotation.term.clone());
                            }
                        }
                    }
                }

                if !split_by_parent_groups.is_empty() {
                    gene_details.split_by_parent_groups
                        .insert(cv_name.clone(), split_by_parent_groups);
                }
            }
        }
    }

    // make gene subsets for genes the are not in a slim category
    fn make_non_slim_subset(&self, cv_name: &FlexStr, slim_subset: &TermSubsetDetails)
                            -> IdGeneSubsetMap
    {
        let slim_termid_set: HashSet<FlexStr> =
            slim_subset.elements.keys().cloned().collect();

        let mut non_slim_with_bp_annotation = HashSet::new();
        let mut non_slim_without_bp_annotation = HashSet::new();

        let has_parent_in_slim = |term_annotations: &Vec<OntTermAnnotations>| {
            for term_annotation in term_annotations {
                let interesting_parent_ids =
                    &self.terms[&term_annotation.term].interesting_parent_ids;
                if !term_annotation.is_not &&
                    (slim_termid_set.contains(&term_annotation.term) ||
                     interesting_parent_ids.intersection(&slim_termid_set).count() > 0)
                {
                    return true;
                }
            }
            false
        };

        for gene_details in self.genes.values() {
            if let Some(load_organism_taxonid) = self.config.load_organism_taxonid
                && load_organism_taxonid != gene_details.taxonid {
                    continue;
                }

            if gene_details.feature_type != "mRNA gene" {
                continue;
            }

            if gene_details.characterisation_status == Some(flex_str!("transposon")) ||
                gene_details.characterisation_status == Some(flex_str!("dubious"))
            {
                continue;
            }

            let mut bp_count = 0;

            if let Some(annotations) =
                gene_details.cv_annotations.get(cv_name) {
                    if has_parent_in_slim(annotations) {
                        continue
                    }
                    bp_count = annotations.len();
                }

            if bp_count == 0 {
                non_slim_without_bp_annotation.insert(gene_details.uniquename.clone());
            } else {
                non_slim_with_bp_annotation.insert(gene_details.uniquename.clone());
            }
        }

        let mut return_map = HashMap::new();

        let cv_display_name = str::replace(cv_name, "_", " ");

        let with_annotation_display_name =
            String::from("Gene products with ") + &cv_display_name +
            " annotation that are not in a slim category";
        let name = flex_fmt!("non_slim_with_{}_annotation", cv_name);
        return_map.insert(name.clone(),
                          GeneSubsetDetails {
                              name,
                              display_name: with_annotation_display_name.to_shared_str(),
                              elements: non_slim_with_bp_annotation,
                          });
        let without_annotation_display_name =
            String::from("Gene products with no ") + &cv_display_name +
            " annotation and are not in a slim category";
        let name = flex_fmt!("non_slim_without_{}_annotation", cv_name);
        return_map.insert(name.clone(),
                          GeneSubsetDetails {
                              name,
                              display_name: without_annotation_display_name.to_shared_str(),
                              elements: non_slim_without_bp_annotation,
                          });
        return_map
    }

    fn make_slim_subset(&self, slim_name: &FlexStr) -> TermSubsetDetails {
        let mut all_genes = HashSet::new();
        let mut all_single_locus_genes = HashSet::new();
        let mut slim_subset: HashMap<TermId, TermSubsetElement> = HashMap::new();
        let slim_config = self.config.slims.get(slim_name)
            .unwrap_or_else(|| panic!("no slim config for {}", slim_name));
        for slim_conf in &slim_config.terms {
            let slim_termid = &slim_conf.termid;
            let term_details = self.terms.get(slim_termid)
                .unwrap_or_else(|| panic!("can't find TermDetails for {}", slim_termid));

            let subset_element = TermSubsetElement {
                name: term_details.name.clone(),
                gene_count: term_details.annotated_genes.len(),
                single_locus_gene_count: term_details.single_locus_annotated_genes.len(),
            };

            for gene in &term_details.annotated_genes {
                all_genes.insert(gene);
            }

            if term_details.annotation_feature_type == "gene" {
                for gene in &term_details.single_locus_annotated_genes {
                    all_single_locus_genes.insert(gene);
                }
            }

            slim_subset.insert(slim_termid.clone(), subset_element);
        }

        TermSubsetDetails {
            name: slim_name.clone(),
            total_gene_count: all_genes.len(),
            total_single_locus_gene_count: all_genes.len(),
            elements: slim_subset,
        }
    }

    fn make_feature_type_subsets(&self, subsets: &mut IdGeneSubsetMap) {
        let mut add_to_subset = |subset_name: &str, gene_details: &GeneDetails| {
            let re = Regex::new(r"[\s,:]+").unwrap();
            let subset_name_no_spaces = re.replace_all(subset_name, "_").to_shared_str();
            subsets.entry(subset_name_no_spaces.clone())
                .or_insert(GeneSubsetDetails {
                    name: subset_name_no_spaces,
                    display_name: subset_name.into(),
                    elements: HashSet::new()
                })
                .elements.insert(gene_details.uniquename.clone());
        };
        for gene_details in self.genes.values() {
            if let Some(load_organism_taxonid) = self.config.load_organism_taxonid
                && load_organism_taxonid != gene_details.taxonid {
                    continue;
                }
            let subset_name =
                format!("feature_type:{}", gene_details.feature_type);
            add_to_subset(&subset_name, gene_details);

            if let Some(parent_types) =
               self.config.feature_sub_groups.get(gene_details.feature_type.as_ref())
            {
                for parent_type in parent_types.iter() {
                    let subset_name = format!("feature_type:{}", parent_type);
                    add_to_subset(&subset_name, gene_details);
                }
            }
        }
    }

    // make subsets using the characterisation_status field of GeneDetails
    fn make_characterisation_status_subsets(&self, subsets: &mut IdGeneSubsetMap) {
        let status_fix_re = Regex::new(r"[\s,:]+").unwrap();
        for gene_details in self.genes.values() {
            if let Some(load_organism_taxonid) = self.config.load_organism_taxonid
                && load_organism_taxonid != gene_details.taxonid {
                    continue;
                }

            if gene_details.feature_type != "mRNA gene" {
                continue;
            }

            if let Some(ref characterisation_status) = gene_details.characterisation_status {
                let subset_name =
                    flex_str!("characterisation_status:") + characterisation_status;
                let subset_name_no_spaces = status_fix_re.replace_all(&subset_name, "_").to_shared_str();
                subsets.entry(subset_name_no_spaces.clone())
                    .or_insert(GeneSubsetDetails {
                        name: subset_name_no_spaces,
                        display_name: subset_name.to_shared_str(),
                        elements: HashSet::new()
                    })
                    .elements.insert(gene_details.uniquename.clone());
            }
        }
    }

    // make InterPro subsets using the interpro_matches field of GeneDetails
    fn make_interpro_subsets(&mut self, subsets: &mut IdGeneSubsetMap) {
        for (gene_uniquename, gene_details) in &self.genes {
            if self.config.load_organism_taxonid.is_none() ||
                self.config.load_organism_taxonid.unwrap() != gene_details.taxonid {
                    continue;
                }

            for interpro_match in &gene_details.interpro_matches {

                let mut new_subset_names = vec![];

                if let (Some(interpro_id), Some(interpro_name)) =
                    (&interpro_match.interpro_id, &interpro_match.interpro_name)
                    && !interpro_id.is_empty() {
                        let subset_name =
                            String::from("interpro:") + interpro_id;
                        new_subset_names.push((subset_name.to_shared_str(),
                                               interpro_name.clone()));
                    }

                let subset_name = format!("interpro:{}:{}",
                                          interpro_match.dbname,
                                          interpro_match.id);
                let subset_display_name = interpro_match.name.as_ref()
                    .or(interpro_match.description.as_ref())
                    .map(|s| s.as_str())
                    .unwrap_or(subset_name.as_str())
                    .into();
                new_subset_names.push((subset_name.to_shared_str(), subset_display_name));

                for (subset_name, display_name) in new_subset_names {
                    subsets.entry(subset_name.clone())
                        .or_insert(GeneSubsetDetails {
                            name: subset_name,
                            display_name,
                            elements: HashSet::new(),
                        })
                        .elements.insert(gene_uniquename.clone());
                }
            }
        }
    }

    // populate the subsets HashMap
    fn make_subsets(&mut self) {
        let mut gene_subsets: IdGeneSubsetMap = HashMap::new();

        for (slim_name, slim_config) in &self.config.slims {
            let slim_subset = self.make_slim_subset(slim_name);
            let gene_subset = self.make_non_slim_subset(&slim_config.cv_name, &slim_subset);
            gene_subsets.extend(gene_subset);
            self.term_subsets.insert(slim_name.clone(), slim_subset);
        }

        self.make_feature_type_subsets(&mut gene_subsets);
        self.make_characterisation_status_subsets(&mut gene_subsets);
        self.make_interpro_subsets(&mut gene_subsets);

        self.gene_subsets = gene_subsets;
    }

    // sort the list of genes in the ChromosomeDetails by start_pos
    pub fn sort_chromosome_genes(&mut self) {
        let mut genes_to_sort: HashMap<ChromosomeName, Vec<GeneUniquename>> =
            HashMap::new();

        {
            let sorter = |uniquename1: &GeneUniquename, uniquename2: &GeneUniquename| {
                let gene1 = &self.genes[uniquename1];
                let gene2 = &self.genes[uniquename2];

                if let Some(ref gene1_loc) = gene1.location
                    && let Some(ref gene2_loc) = gene2.location {
                        let cmp = gene1_loc.start_pos.cmp(&gene2_loc.start_pos);
                        if cmp != Ordering::Equal {
                            return cmp;
                        }
                    }
                if gene1.name.is_some() {
                    if gene2.name.is_some() {
                        gene1.name.cmp(&gene2.name)
                    } else {
                        Ordering::Less
                    }
                } else if gene2.name.is_some() {
                    Ordering::Greater
                } else {
                    gene1.uniquename.cmp(&gene2.uniquename)
                }
            };

            for (chr_uniquename, chr_details) in &self.chromosomes {
                genes_to_sort.insert(chr_uniquename.clone(),
                                     chr_details.gene_uniquenames.clone());
            }

            for gene_uniquenames in genes_to_sort.values_mut() {
                gene_uniquenames.sort_by(&sorter);
            }
        }

        for (chr_uniquename, gene_uniquenames) in genes_to_sort {
            self.chromosomes.get_mut(&chr_uniquename).unwrap().gene_uniquenames =
                gene_uniquenames;
        }
    }

    fn get_dataset_name_for_measurement(&self, reference_uniquename: &FlexStr,
                                        level_type_termid: &FlexStr,
                                        during_termid: &FlexStr, scale: &FlexStr)
                                        -> Option<FlexStr>
    {
        for conf in &self.config.gene_expression.datasets {
            if conf.pubmed_id == *reference_uniquename &&
                conf.level_type_termid == *level_type_termid &&
                conf.during_termid == *during_termid &&
                conf.scale == *scale {
                    return Some(conf.name.clone());
                }
        }

        None
    }

    fn set_gene_expression_measurements(&mut self) {
       let mut measurements = HashMap::new();

        for annotation in &self.ont_annotations {
            if &annotation.term_short.cv_name != "gene_ex" {
                continue;
            }

            let gene_uniquename =
                if let Some(gene_short) = annotation.genes.iter().next() {
                    gene_short.uniquename.clone()
                } else {
                    continue;
                };

            let level_type_termid = annotation.term_short.termid.clone();

            let reference_uniquename =
                if let Some(ref_short) = &annotation.reference_short {
                    ref_short.uniquename.clone()
                } else {
                    continue;
                };

            let mut during_ext = None;

            for extpart in &annotation.extension {
                if extpart.rel_type_name == "during" {
                    during_ext = Some(&extpart.ext_range);
                }
            }

            let during_termid =
                if let Some(ExtRange::Term(termid)) = during_ext {
                    termid.clone()
                } else {
                    continue;
                };

            let gene_ex_props =
                if let Some(ref props) = annotation.gene_ex_props {
                    props
                } else {
                    continue;
                };

            let scale = gene_ex_props.scale.clone();

            let copies_per_cell =
                gene_ex_props.copies_per_cell.as_ref().cloned();

            let avg_copies_per_cell =
                gene_ex_props.avg_copies_per_cell.as_ref().cloned();

            if let Some(dataset_name) =
                self.get_dataset_name_for_measurement(&reference_uniquename,
                                                      &level_type_termid,
                                                      &during_termid, &scale) {
                    measurements
                        .entry(gene_uniquename)
                        .or_insert_with(HashMap::new)
                        .insert(dataset_name,
                                GeneExMeasurement {
                                    reference_uniquename,
                                    level_type_termid,
                                    during_termid,
                                    copies_per_cell,
                                    avg_copies_per_cell,
                                    scale
                                });
                }
        }

        self.gene_expression_measurements = measurements;
    }

    fn set_chromosome_gene_counts(&mut self) {
        let mut counts = HashMap::new();
        let mut coding_counts = HashMap::new();

        for gene_details in self.genes.values() {
            if let Some(ref loc) = gene_details.location {
                *counts
                    .entry(&loc.chromosome_name)
                    .or_insert(0) += 1;
            }

            if gene_details.feature_type == "mRNA gene"
                && let Some(ref loc) = gene_details.location {
                    *coding_counts
                        .entry(&loc.chromosome_name)
                        .or_insert(0) += 1;
                }
        }

        for chromosome_detail in self.chromosomes.values_mut() {
            if let Some(count) = counts.get(&chromosome_detail.name) {
                chromosome_detail.gene_count = *count;
            }

            if let Some(count) = coding_counts.get(&chromosome_detail.name) {
                chromosome_detail.coding_gene_count = *count;
            }
        }
    }

    fn macromolecular_complex_data(&self, config: &Config)
                                   -> ProteinComplexData
    {
        let mut complex_data = HashMap::new();

        let no_evidence = flex_str!("NO_EVIDENCE");

        let Some(ref complexes_config) = config.file_exports.macromolecular_complexes
        else {
            return HashMap::new();
        };

        let check_parent_term = |termid: &FlexStr| {
            *termid == complexes_config.parent_complex_termid
        };

        for (gene_uniquename, gene_details) in &self.genes {
            let uniprot_identifier = &gene_details.uniprot_identifier;
            for term_annotations in gene_details.cv_annotations.values() {
                'TERM_ANNOTATION: for term_annotation in term_annotations {
                    if term_annotation.is_not {
                        continue 'TERM_ANNOTATION;
                    }
                    let termid = &term_annotation.term;
                    'TERM: for annotation_detail_id in &term_annotation.annotations {
                        let annotation_detail = self.annotation_details
                            .get(annotation_detail_id).expect("can't find OntAnnotationDetail");

                        let term_details = self.terms.get(termid).unwrap();
                        let term_name = &term_details.name;

                        let reference_uniquename =
                            annotation_detail.reference.clone();

                        if complexes_config.excluded_terms.contains(termid) {
                            continue 'TERM;
                        }

                        if !term_details.interesting_parent_ids.iter().any(check_parent_term) {
                            continue 'TERM;
                        }

                        let evidence = annotation_detail.evidence.clone()
                            .unwrap_or_else(|| no_evidence.clone());
                        let gene_short = gene_details.into();

                        complex_data.entry(termid.clone())
                            .or_insert_with(|| ProteinComplexTerm {
                                term_name: term_name.to_owned(),
                                complex_genes: BTreeMap::new(),
                            })
                            .complex_genes
                            .entry(gene_uniquename.clone())
                            .or_insert_with(|| ProteinComplexGene {
                                gene_short,
                                uniprot_identifier: uniprot_identifier.to_owned(),
                                annotation_details: HashSet::new(),
                            })
                            .annotation_details
                            .insert((reference_uniquename, annotation_detail.assigned_by.clone(),
                                     evidence));
                    }
                }
            }
        }

        complex_data
    }

    // remove some of the refs that have no annotations.
    // See: https://github.com/pombase/website/issues/628
    fn remove_non_curatable_refs(&mut self) {
        let filtered_refs = self.references.drain()
            .filter(|(_, reference_details)| {
                if reference_has_annotation(reference_details) {
                    return true;
                }
                if let Some(ref triage_status) = reference_details.canto_triage_status
                    && (triage_status == "New" || triage_status == "Wrong organism" && triage_status == "Loaded in error"){
                        return false;
                    }
                // default to true because there are references that
                // haven't or shouldn't be triaged, eg. GO_REF:...
                true
            })
            .collect();

        self.references = filtered_refs;
    }

    fn make_gocam_summaries(&mut self) {
        self.gocam_summaries = self.gocam_models.iter()
            .map(|m| {
                let summ = GoCamSummary::new_from_model(m, &self.orcid_name_map,
                                                        &self.terms, &self.children_by_termid);
                (summ.gocam_id.clone(), summ)
            })
            .collect()
    }

    fn set_gene_and_term_gocams(&mut self) {
        let mut gocams_of_genes = HashMap::new();
        let mut gocams_of_terms = HashMap::new();

        for gocam_details in self.gocam_summaries.values() {
            for gene_uniquename in &gocam_details.activity_enabling_genes {
                gocams_of_genes.entry(gene_uniquename.clone())
                   .or_insert_with(HashSet::new)
                   .insert(gocam_details.into());
            }

            for cvterm in &gocam_details.title_child_process_terms {
                gocams_of_terms.entry(cvterm.termid.clone())
                   .or_insert_with(HashSet::new)
                   .insert(gocam_details.into());
            }
        }

        for (gene_uniquename, gocam_id_and_title) in gocams_of_genes.drain() {
            if let Some(ref mut gene_details) = self.genes.get_mut(&gene_uniquename) {
                gene_details.gocams = gocam_id_and_title;
            }
        }

        for (termid, gocam_id_and_title) in gocams_of_terms.drain() {
            self.terms.get_mut(&termid).unwrap().gocams = gocam_id_and_title;
        }

    }

    fn store_genetic_interactions(&mut self) {
        fn get_interaction_details(interaction_annotations: &mut GeneticInteractionMap,
                               interaction_key: GeneticInteractionKey)
             -> &mut Vec<GeneticInteractionDetail>
        {
            interaction_annotations.entry(interaction_key).or_default()
        }

        for (interaction_key, interaction_annotation_details) in self.genetic_interaction_annotations.drain() {
            let gene_a_uniquename = &interaction_key.gene_a_uniquename;
            let gene_a = self.genes.get_mut(gene_a_uniquename).unwrap();
            get_interaction_details(&mut gene_a.genetic_interactions, interaction_key.clone())
                 .extend(interaction_annotation_details.clone());

            let gene_b_uniquename = &interaction_key.gene_b_uniquename;
            let gene_b = self.genes.get_mut(gene_b_uniquename).unwrap();
            get_interaction_details(&mut gene_b.genetic_interactions, interaction_key.clone())
                 .extend(interaction_annotation_details.clone());


            for interaction_detail in interaction_annotation_details {
                if let Some(ref reference_uniquename) =
                    interaction_detail.reference_uniquename
                {
                    let reference = self.references.get_mut(reference_uniquename).unwrap();
                    get_interaction_details(&mut reference.genetic_interactions, interaction_key.clone())
                        .push(interaction_detail.clone());
                }

                if let Some(ref double_mutant_phenotype_termid) =
                    interaction_detail.double_mutant_phenotype_termid
                {
                    let term_details = self.terms.get_mut(double_mutant_phenotype_termid).unwrap();
                    get_interaction_details(&mut term_details.double_mutant_genetic_interactions, interaction_key.clone())
                        .push(interaction_detail.clone());
                }

                if let Some(ref rescued_phenotype_termid) =
                    interaction_detail.rescued_phenotype_termid
                {
                    let term_details = self.terms.get_mut(rescued_phenotype_termid).unwrap();
                    get_interaction_details(&mut term_details.single_allele_genetic_interactions, interaction_key.clone())
                        .push(interaction_detail.clone());
                }

                if let Some(ref double_mutant_genotype_uniquename) =
                     interaction_detail.double_mutant_genotype_display_uniquename
                {
                    let genotype_details = self.genotypes.get_mut(double_mutant_genotype_uniquename).unwrap();
                    get_interaction_details(&mut genotype_details.double_mutant_genetic_interactions,
                                            interaction_key.clone())
                                            .push(interaction_detail.clone());
                }

                if let Some(ref rescue_genotype_uniquename) =
                     interaction_detail.genotype_a_uniquename
                {
                    let genotype_details = self.genotypes.get_mut(rescue_genotype_uniquename).unwrap();
                    get_interaction_details(&mut genotype_details.rescue_genetic_interactions,
                                            interaction_key.clone())
                                            .push(interaction_detail.clone());
                }
            }
        }
    }

    fn make_gene_summaries(&mut self) -> (Vec<GeneSummary>, Vec<SolrGeneSummary>) {
        let mut gene_summaries = vec![];
        let mut solr_gene_summaries = vec![];

        for (gene_uniquename, gene_details) in &self.genes {
            if self.config.load_organism_taxonid.is_none() ||
                self.config.load_organism_taxonid.unwrap() == gene_details.taxonid {
                    let gene_summary = self.make_gene_summary(gene_uniquename);
                    let solr_gene_summary =
                        SolrGeneSummary {
                            id: gene_summary.uniquename.clone(),
                            name: gene_summary.name.clone(),
                            taxonid: gene_summary.taxonid,
                            product: gene_summary.product.clone(),
                            uniprot_identifier: gene_summary.uniprot_identifier.clone(),
                            synonyms: gene_summary.synonyms.clone(),
                            feature_type: gene_summary.feature_type.clone(),
                        };
                    gene_summaries.push(gene_summary);
                    solr_gene_summaries.push(solr_gene_summary);
                }
        }

        (gene_summaries, solr_gene_summaries)
    }

    fn make_solr_allele_summaries(&mut self) -> Vec<SolrAlleleSummary> {
        self.alleles.values().map(|details| {
            SolrAlleleSummary {
                id: details.uniquename.clone(),
                name: details.name.clone(),
                allele_type: details.allele_type.clone(),
                description: details.description.clone(),
                gene_uniquename: details.gene.uniquename.clone(),
                gene_name: details.gene.name.clone(),
                synonyms: details.synonyms.iter().map(|s| s.name.clone()).collect(),
                highlighting: HashMap::new(),
            }
        })
        .collect()
    }

    fn make_solr_term_summaries(&mut self) -> Vec<SolrTermSummary> {
        let mut return_summaries = vec![];

        let term_name_split_re = Regex::new(r"\W+").unwrap();

        for (termid, term_details) in &self.terms {
            if term_details.is_obsolete {
                continue;
            }

            let trimmable_p = |c: char| {
                c.is_whitespace() || c == ',' || c == ':'
                    || c == ';' || c == '.' || c == '\''
            };

            let term_name_words =
                term_name_split_re.split(&term_details.name)
                .map(|s: &str| {
                    s.trim_matches(&trimmable_p).to_shared_str()
                }).collect::<Vec<_>>();

            let mut exact_synonyms = vec![];
            let mut exact_synonym_words_vec: Vec<FlexStr> = vec![];
            let mut narrow_synonyms = vec![];
            let mut narrow_synonym_words_vec: Vec<FlexStr> = vec![];
            let mut distant_synonyms = vec![];
            let mut distant_synonym_words_vec: Vec<FlexStr> = vec![];

            let add_to_words_vec = |synonym: &FlexStr, words_vec: &mut Vec<FlexStr>| {
                let synonym_words = term_name_split_re.split(synonym);
                for word in synonym_words {
                    let word_string = word.trim_matches(&trimmable_p).to_shared_str();
                    if !words_vec.contains(&word_string) &&
                        !term_name_words.contains(&word_string) {
                            words_vec.push(word_string);
                        }
                }
            };

            for synonym in &term_details.synonyms {
                match synonym.synonym_type.as_str() {
                    "exact" => {
                        add_to_words_vec(&synonym.name, &mut exact_synonym_words_vec);
                        exact_synonyms.push(synonym.name.clone());
                    },
                    "narrow" => {
                        add_to_words_vec(&synonym.name, &mut narrow_synonym_words_vec);
                        narrow_synonyms.push(synonym.name.clone());
                    },
                    _ => {
                        add_to_words_vec(&synonym.name, &mut distant_synonym_words_vec);
                        distant_synonyms.push(synonym.name.clone());
                    }
                }
            }

            distant_synonyms = distant_synonyms.into_iter()
                .filter(|synonym| {
                    !exact_synonyms.contains(synonym)
                })
                .collect::<Vec<_>>();

            let annotation_count = term_details.annotation_count();
            let interesting_parent_ids_for_solr =
                term_details.interesting_parent_ids.clone();

            let gocam_ids =
                term_details.gocams.iter().map(|gocam| gocam.gocam_id.clone()).collect();

            let term_summ = SolrTermSummary {
                id: termid.clone(),
                cv_name: term_details.cv_name.clone(),
                name: term_details.name.clone(),
                name_str_field: term_details.name.clone(),
                definition: term_details.definition.clone(),
                exact_synonyms_str_field: exact_synonyms.clone(),
                exact_synonyms,
                exact_synonym_words: join(&exact_synonym_words_vec," "),
                narrow_synonyms,
                narrow_synonym_words: join(&narrow_synonym_words_vec," "),
                distant_synonyms,
                distant_synonym_words: join(&distant_synonym_words_vec, " "),
                interesting_parent_ids: interesting_parent_ids_for_solr,
                definition_xrefs: term_details.definition_xrefs.clone(),
                secondary_identifiers: term_details.secondary_identifiers.clone(),
                gocam_ids,
                annotation_count,
                gene_count: term_details.gene_count,
                genotype_count: term_details.genotype_count,
                highlighting: HashMap::new(),
            };
            return_summaries.push(term_summ);
        }

        return_summaries
    }

    fn make_solr_reference_summaries(&mut self) -> Vec<SolrReferenceSummary> {
        let mut return_summaries = vec![];

        for reference_details in self.references.values() {
            return_summaries.push(SolrReferenceSummary::from_reference_details(reference_details));
        }

        return_summaries
    }

    fn get_stats(&self) -> Stats {
        let mut by_taxon_collector = HashMap::new();

        for gene_details in self.genes.values() {
            let taxonid = gene_details.taxonid;

            let count = by_taxon_collector
                .entry(taxonid)
                .or_insert(0);
            *count += 1usize;
        }

        let mut by_taxon = HashMap::new();

        for (taxonid, counts) in by_taxon_collector.iter() {
            let mut total_annotation_count = 0;

            let annotation_type_counts_by_year =
                if self.config.load_organism_taxonid == Some(*taxonid) {
                    let annotation_counts = self.chado_queries.annotation_type_counts_by_year.clone();
                    for (_, row_data) in &annotation_counts.data {
                        for count in row_data {
                            total_annotation_count += count;
                        }
                    }
                    Some(annotation_counts)
                } else {
                    None
                };


            by_taxon.insert(*taxonid, StatCountsByTaxon {
                genes: *counts,
                annotations: total_annotation_count,
                annotation_type_counts_by_year,
            });
        }

        Stats {
            by_taxon,
            community_pubs_count: self.all_community_curated.len(),
            non_community_pubs_count: self.all_admin_curated.len(),
        }
    }

    fn get_pub_curated_stats(&self) -> CuratedStats
    {
        let mut year_month_map = HashMap::new();
        let mut year_map = HashMap::new();

        let mut lowest_year = "9999".to_owned();
        let mut highest_year = "0000".to_owned();

        for ref_details in self.references.values() {
            let Some(ref canto_triage_status) = ref_details.canto_triage_status
            else {
                continue;
            };

            let is_curatable = canto_triage_status == "Curatable" ||
                ref_details.annotation_count() > 0;

            let is_community_curated = ref_details.canto_curator_role == "community";

            let Some(ref pubmed_entrez_date) = ref_details.pubmed_entrez_date
            else {
                continue;
            };

            let Some(captures) = ISO_DATE_RE.captures(pubmed_entrez_date)
            else {
                continue;
            };

            let (Some(entrez_date_year), Some(entrez_date_month)) =
                (captures.get(1), captures.get(2))
            else {
                continue;
            };

            let entrez_date_year = entrez_date_year.as_str();
            if entrez_date_year.cmp(&lowest_year) == Ordering::Less {
                lowest_year = entrez_date_year.to_owned();
            }
            if entrez_date_year.cmp(&highest_year) == Ordering::Greater {
                highest_year = entrez_date_year.to_owned();
            }

            let entrez_date_year_month =
                format!("{}-{}", entrez_date_year, entrez_date_month.as_str());

            if is_curatable {

            year_month_map.entry(entrez_date_year_month)
                .or_insert((0, 0, 0, 0))
                .0 += 1;

            year_map.entry(entrez_date_year.to_owned())
                .or_insert((0, 0, 0, 0))
                .0 += 1;

            } else {

            year_month_map.entry(entrez_date_year_month)
                .or_insert((0, 0, 0, 0))
                .3 += 1;

            year_map.entry(entrez_date_year.to_owned())
                .or_insert((0, 0, 0, 0))
                .3 += 1;

            }
            let Some(ref submitted_date) = ref_details.canto_session_submitted_date
            else {
                continue;
            };

            let Some(captures) = ISO_DATE_RE.captures(submitted_date)
            else {
                continue;
            };

            let (Some(submitted_year), Some(submitted_month)) =
                (captures.get(1), captures.get(2))
            else {
                continue;
            };

            let submitted_year = submitted_year.as_str();

            if submitted_year.cmp(&lowest_year) == Ordering::Less {
                lowest_year = submitted_year.to_owned();
            }
            if submitted_year.cmp(&highest_year) == Ordering::Greater {
                highest_year = submitted_year.to_owned();
            }

            let submitted_year_month =
                format!("{}-{}", submitted_year, submitted_month.as_str());

            let value = year_month_map.entry(submitted_year_month)
                .or_insert((0, 0, 0, 0));
            if is_community_curated {
                value.1 += 1;
            } else {
                value.2 += 1;
            }

            let value = year_map.entry(submitted_year.to_owned())
                .or_insert((0, 0, 0, 0));
            if is_community_curated {
                value.1 += 1;
            } else {
                value.2 += 1;
            }
        }

        let mut year_months = vec![];

        let lowest_year: usize = lowest_year.parse().unwrap();
        let highest_year: usize = highest_year.parse().unwrap();
        for year in lowest_year..=highest_year {
           for month in 1..=12 {
              year_months.push(format!("{}-{:02}", year, month));
           }
        }

        let mut year_month_return = vec![];

        for year_month in &year_months {
            let (curatable, community_curated, admin_curated, uncuratable) =
                year_month_map.get(year_month).unwrap_or(&(0,0,0,0));

            year_month_return.push((year_month.to_owned(), vec![*curatable, *community_curated, *admin_curated, *uncuratable]));
        }

        let mut year_return = vec![];

        for year in lowest_year..=highest_year {
            let year = year.to_string();

            let (curatable, community_curated, admin_curated, uncuratable) =
                year_map.get(&year).unwrap_or(&(0,0,0,0));

            year_return.push((year.to_owned(), vec![*curatable, *community_curated, *admin_curated, *uncuratable]));
        }

        (year_month_return, year_return)
    }

    fn get_cumulative_curated_stats(&self, stats: &[(DateString, Vec<usize>)])
      -> Vec<(DateString, Vec<usize>)>
    {
        if stats.is_empty() {
            return vec![];
        }

        let mut return_vec = vec![];

        let (first_date_string, sums) = stats[0].to_owned();
        let (mut curatable_sum, mut community_curated_sum, mut admin_curated_sum, mut uncuratable_sum) =
            (sums[0], sums[1], sums[2], sums[3]);

        return_vec.push((first_date_string,
                         vec![curatable_sum, community_curated_sum, admin_curated_sum]));

        for (date_string, counts) in &stats[1..] {
            if let &[this_curatable, this_community_curated, this_admin_curated, this_uncuratable] =
                   &counts[0..] {
              curatable_sum += this_curatable;
              community_curated_sum += this_community_curated;
              admin_curated_sum += this_admin_curated;
              uncuratable_sum += this_uncuratable;
              return_vec.push((date_string.to_owned(),
                               vec![curatable_sum, community_curated_sum,
                                    admin_curated_sum, uncuratable_sum]));
            }
        }

        return_vec
    }

    fn annotations_per_pub(&self)
         -> (Vec<StatsFloatTableRow>, Vec<StatsFloatTableRow>,
             Vec<StatsFloatTableRow>)
    {
        let mut ltp_gene_map = HashMap::new();
        let mut ltp_map = HashMap::new();
        let mut htp_map = HashMap::new();

        let range_of_year = |year| {
            let range_start = (year-1)/5*5+1;
            let range_end = range_start + 4;
            format!("{range_start}-{range_end}")
        };

        for reference in self.references.values() {
            let Some(ref canto_annotation_status) = reference.canto_annotation_status
            else {
                continue;
            };

            if canto_annotation_status != "APPROVED" {
                continue;
            }

            let Some(ref pub_year) = reference.publication_year
            else {
                continue;
            };

            let Ok(pub_year) = pub_year.parse::<usize>()
            else {
                continue;
            };


            let year_range = range_of_year(pub_year);

            ltp_gene_map.entry(year_range.clone())
                .or_insert_with(HashMap::new)
                .insert(reference.uniquename.clone(), reference.ltp_gene_count);

            let mut ltp_count = 0usize;
            let mut htp_count = 0usize;

            let mut add_count = |throughput: &Throughput| {
                match throughput {
                    Throughput::LowThroughput => ltp_count += 1,
                    Throughput::HighThroughput => htp_count += 1,
                    _ => (),
                }
            };

            for term_annotations_vec in reference.cv_annotations.values() {
                for term_annotations in term_annotations_vec {
                    for annotation_id in &term_annotations.annotations {
                        let annotation = self.annotation_details.get(annotation_id).unwrap();

                        let Some(ref throughput) = annotation.throughput
                        else {
                            continue;
                        };

                        add_count(throughput);
                    }
                }
            }

            for annotation in &reference.physical_interactions {
                let Some(ref throughput) = annotation.throughput
                else {
                    continue;
                };

                add_count(throughput);
            }

            for annotations in reference.genetic_interactions.values() {
                for annotation in annotations {
                   let Some(ref throughput) = annotation.throughput
                    else {
                        continue;
                    };

                    add_count(throughput);
                }
            }

            if ltp_count > 0 {
                ltp_map.entry(year_range.clone())
                       .or_insert_with(HashMap::new)
                       .insert(reference.uniquename.clone(), ltp_count);
            }

            if htp_count > 0 {
                htp_map.entry(year_range.clone())
                       .or_insert_with(HashMap::new)
                       .insert(reference.uniquename.clone(), htp_count);
            }
        }

        let mut ranges: BTreeSet<String> = BTreeSet::new();

        ranges.extend(ltp_gene_map.keys().cloned());
        ranges.extend(ltp_map.keys().cloned());
        ranges.extend(htp_map.keys().cloned());

        let mut genes_data: Vec<StatsFloatTableRow> = vec![];
        let mut ltp_annotations_data: Vec<StatsFloatTableRow> = vec![];
        let mut htp_annotations_data: Vec<StatsFloatTableRow> = vec![];

        for range in ranges.iter() {
            for (map, data) in [(&mut ltp_gene_map, &mut genes_data),
                                (&mut ltp_map, &mut ltp_annotations_data),
                                (&mut htp_map, &mut htp_annotations_data)] {
                if let Some(gene_ref_counts) = map.get(range) {
                    let non_zero_values: Vec<f32> = gene_ref_counts.values()
                       .map(|v| *v as f32)
                       .collect();
                    let sum: f32 = non_zero_values.iter().sum();
                    let average = sum / non_zero_values.len() as f32;
                    data.push((range.to_owned(), average));
                } else {
                    data.push((range.to_owned(), 0.0f32))
                }
            }
        }

        (genes_data, ltp_annotations_data, htp_annotations_data)
    }

    fn get_micropublications_by_year(&self)
        -> (StatsIntegerTable, StatsIntegerTable)
    {
        let mut by_year_map = BTreeMap::new();

        for reference in self.references.values() {
            let Some(ref pub_year) = reference.publication_year
            else {
                continue;
            };

            let Some(ref citation) = reference.citation
            else {
                continue;
            };

            if !citation.contains("MicroPubl Biol") {
                continue;
            }

            *by_year_map.entry(pub_year.clone())
                .or_insert(0) += 1;
        }

        let header = vec!["date".to_owned(), "count".to_owned()];
        let mut data = vec![];
        let mut cumulative_data = vec![];
        let mut total = 0;

        for (year, count) in by_year_map.iter() {
            data.push((year.to_std_string(), vec![*count]));
            total += *count;
            cumulative_data.push((year.to_std_string(), vec![total]));
        }

        let by_year = StatsIntegerTable {
            header: header.clone(),
            data,
        };

        let cumulative_by_year = StatsIntegerTable {
            header: header.clone(),
            data: cumulative_data,
        };

        (by_year, cumulative_by_year)
    }

    fn get_detailed_stats(&self) -> DetailedStats {
        let pub_stats_header = vec!["date".to_owned(), "curatable".to_owned(),
                                    "community_curated".to_owned(),
                                    "admin_curated".to_owned(),
                                    "uncuratable".to_owned()];

        let (curated_by_month, curated_by_year) = self.get_pub_curated_stats();
        let cumulative_curated_by_month =
            self.get_cumulative_curated_stats(&curated_by_month);
        let cumulative_curated_by_year =
            self.get_cumulative_curated_stats(&curated_by_year);

        let annotations_per_year_header =
            vec!["year_range".to_owned(), "average_annotations_per_publication".to_owned()];
        let (ltp_genes_per_pub_per_year_range,
             ltp_annotations_per_year_range, htp_annotations_per_year_range) =
            self.annotations_per_pub();

        let annotation_type_counts_by_year =
            self.chado_queries.annotation_type_counts_by_year.clone();
        let cumulative_annotation_type_counts_by_year =
            get_cumulative_annotation_type_counts(annotation_type_counts_by_year.clone());
        let (micropublications_by_year, cumulative_micropublications_by_year) =
            self.get_micropublications_by_year();

        DetailedStats {
            curated_by_month: StatsIntegerTable {
                 header: pub_stats_header.clone(),
                 data: curated_by_month,
            },
            cumulative_curated_by_month: StatsIntegerTable {
                 header: pub_stats_header.clone(),
                 data: cumulative_curated_by_month,
            },
            curated_by_year: StatsIntegerTable {
                 header: pub_stats_header.clone(),
                 data: curated_by_year,
            },
            cumulative_curated_by_year: StatsIntegerTable {
                 header: pub_stats_header.clone(),
                 data: cumulative_curated_by_year,
            },
            ltp_genes_per_pub_per_year_range: StatsFloatTable {
                 header: vec!["year_range".to_owned(), "genes_per_pub".to_owned()],
                 data: ltp_genes_per_pub_per_year_range,
            },
            ltp_annotations_per_pub_per_year_range: StatsFloatTable {
                 header: annotations_per_year_header.clone(),
                 data: ltp_annotations_per_year_range,
            },
            htp_annotations_per_pub_per_year_range: StatsFloatTable {
                 header: annotations_per_year_header,
                 data: htp_annotations_per_year_range,
            },
            community_response_rates: self.chado_queries.community_response_rates.clone(),
            annotation_type_counts_by_year,
            cumulative_annotation_type_counts_by_year,
            micropublications_by_year,
            cumulative_micropublications_by_year,
        }
    }

    fn set_protein_feature_fields(&mut self) {

        for (gene_uniquename, gene_details) in &mut self.genes {
            for transcript_uniquename in &gene_details.transcripts {
                let Some(transcript) = self.transcripts.get(transcript_uniquename)
                else {
                    continue;
                };
                let Some(ref protein) = transcript.protein
                else {
                    continue;
                };

                let Some(ref uniprot_data) = self.uniprot_data
                else {
                    continue;
                };

                let Some(uniprot_data_entry) = uniprot_data.get(gene_uniquename)
                else {
                    continue;
                };

                let uniprot_sequence = &uniprot_data_entry.sequence;

                gene_details.signal_peptide = uniprot_data_entry.signal_peptide.clone();
                gene_details.transit_peptide = uniprot_data_entry.transit_peptide.clone();

                if &protein.sequence != uniprot_sequence &&
                    !(protein.sequence.ends_with("*") &&
                      &protein.sequence[0..protein.sequence.len() - 1] ==
                      uniprot_sequence.as_str()) {
                    eprintln!("Sequence doesn't match Uniprot for: {}",
                              gene_uniquename);
                    continue;
                }

                gene_details.binding_sites = uniprot_data_entry.binding_sites.clone();
                gene_details.active_sites = uniprot_data_entry.active_sites.clone();
                gene_details.beta_strands = uniprot_data_entry.beta_strands.clone();
                gene_details.helices = uniprot_data_entry.helices.clone();
                gene_details.turns = uniprot_data_entry.turns.clone();
                gene_details.propeptides = uniprot_data_entry.propeptides.clone();
                gene_details.chains = uniprot_data_entry.chains.clone();
                gene_details.glycosylation_sites = uniprot_data_entry.glycosylation_sites.clone();
                gene_details.disulfide_bonds = uniprot_data_entry.disulfide_bonds.clone();
                gene_details.lipidation_sites = uniprot_data_entry.lipidation_sites.clone();
            }
        }
    }

    fn interactions_for_export(&self)
       -> (Vec<InteractionAnnotation>, Vec<InteractionAnnotation>)
    {
        let physical_interaction_annotations: Vec<_> =
            self.physical_interaction_annotations.iter().cloned().collect();
        let mut genetic_interaction_annotations = vec![];

        for (key, interaction_detail_vec) in self.genetic_interaction_annotations.iter() {
            for interaction_detail in interaction_detail_vec.iter() {
                let GeneticInteractionKey {
                    gene_a_uniquename,
                    gene_b_uniquename,
                    interaction_type,
                } = key.clone();

                let interaction = InteractionAnnotation {
                    gene_uniquename: gene_a_uniquename,
                    interactor_uniquename: gene_b_uniquename,
                    evidence: interaction_type,
                    reference_uniquename: interaction_detail.reference_uniquename.clone(),
                    throughput: interaction_detail.throughput.clone(),
                    interaction_note: interaction_detail.interaction_note.clone(),
                    source_database: interaction_detail.source_database.clone(),
                    annotation_date: interaction_detail.annotation_date.clone(),
                };
                genetic_interaction_annotations.push(interaction);
            }
        }

        (physical_interaction_annotations, genetic_interaction_annotations)
    }

    fn report_missing_frameshift_quals(&self) {
        for (transcript_uniquename, pos) in &self.transcript_frameshifts_to_check {
            let Some(gene_uniquename) = self.genes_of_transcripts.get(transcript_uniquename)
            else {
                panic!("transcript has no gene: {}", transcript_uniquename);
            };

            let gene_details = self.genes.get(gene_uniquename).unwrap();

            if !has_annotated_frame_shift(gene_details) {
                eprintln!("no gap between exons at {} in {}", pos, transcript_uniquename);
            }
        }
    }

    pub fn get_web_data(mut self) -> WebData {
        self.process_dbxrefs();
        self.process_references();
        self.process_chromosome_features();
        self.make_feature_rel_maps();
        self.process_features();
        self.add_gene_neighbourhoods();
        self.process_props_from_feature_cvterms();
        self.process_allele_features();
        self.process_genotype_features();
        self.process_cvterms();
        self.add_interesting_parents();
        self.process_cvterm_rels();
        self.process_extension_cvterms();
        self.process_feature_synonyms();
        self.process_feature_publications();
        self.process_feature_cvterms();
        self.remove_duplicate_transcript_annotation();
        self.store_ont_annotations(false);
        self.store_ont_annotations(true);
        self.process_cvtermpath();
        self.process_feature_rels();
        self.add_target_of_annotations();
        self.set_deletion_viability();
        self.set_gene_flags();
        self.set_term_details_subsets();
        self.set_taxonomic_distributions();
        self.remove_non_curatable_refs();

        self.make_gocam_summaries();
        self.set_gene_and_term_gocams();

        let (physical_interaction_annotations,
            genetic_interaction_annotations) = self.interactions_for_export();
        self.store_genetic_interactions();

        self.set_term_details_maps();
        self.set_gene_details_maps();
        self.set_gene_details_subset_termids();
        self.set_allele_details_maps();
        self.set_genotype_details_maps();
        self.set_reference_details_maps();
        self.set_chromosome_gene_counts();
        self.set_counts();
        self.set_split_by_parent_sets();
        self.set_protein_feature_fields();
        self.add_genotypes_to_allele_details();
        self.make_subsets();
        self.sort_chromosome_genes();
        self.set_gene_expression_measurements();

        self.report_missing_frameshift_quals();

        let stats = self.get_stats();

        let detailed_stats = self.get_detailed_stats();

        let metadata = self.make_metadata();

        let (gene_summaries, solr_gene_summaries) = self.make_gene_summaries();

        let solr_allele_summaries = self.make_solr_allele_summaries();
        let solr_term_summaries = self.make_solr_term_summaries();
        let solr_reference_summaries = self.make_solr_reference_summaries();

        let gocam_models = self.gocam_models.clone();

        self.gocam_overlaps = find_activity_overlaps(&gocam_models);
        self.gocam_overlaps_merge_by_chemical =
            find_chemical_overlaps(&gocam_models);
        self.gocam_holes = gocam_models.iter()
            .flat_map(find_holes).collect();

        let solr_data = SolrData {
            term_summaries: solr_term_summaries,
            gene_summaries: solr_gene_summaries,
            allele_summaries: solr_allele_summaries,
            reference_summaries: solr_reference_summaries,
        };

        let chromosomes = self.chromosomes.clone();
        let mut chromosome_summaries = vec![];

        for chr_details in self.chromosomes.values() {
            chromosome_summaries.push(chr_details.make_chromosome_short());
        }

        let recent_references = self.recent_references.clone();
        let all_community_curated = self.all_community_curated.clone();
        let all_admin_curated = self.all_admin_curated.clone();
        let ont_annotations = self.ont_annotations.clone();

        self.protein_complex_data = self.macromolecular_complex_data(self.config);

        let mut terms_for_api: HashMap<TermId, TermDetails> = HashMap::new();

        for (termid, term_details) in &self.terms {
            terms_for_api.insert(termid.clone(), term_details.clone());
        }

        let mut genes = self.genes.clone();
        let alleles = self.alleles.clone();
        let genotypes = self.genotypes.clone();
        let references = self.references.clone();
        let annotation_details = self.annotation_details.clone();

        let termid_genotype_annotation = self.get_api_genotype_annotation();

        for model in &mut self.gocam_models {
            model.add_pro_term_to_gene_map(&self.pro_term_to_gene);
        }

        let api_maps = self.make_api_maps();

        set_has_protein_features(&mut genes, &api_maps.protein_view_data);

        WebData {
            metadata,
            chromosomes,
            chromosome_summaries,
            recent_references,
            all_community_curated,
            all_admin_curated,
            api_maps,
            terms: terms_for_api,
            genes,
            alleles,
            genotypes,
            references,
            annotation_details,
            termid_genotype_annotation,
            search_gene_summaries: gene_summaries,
            solr_data,
            ont_annotations,
            stats,
            detailed_stats,
            gocam_models,

            physical_interaction_annotations,
            genetic_interaction_annotations,

            // used to implement the DataLookup trait:
            arc_terms: Arc::new(RwLock::new(HashMap::new())),
            arc_genes: Arc::new(RwLock::new(HashMap::new())),
            arc_alleles: Arc::new(RwLock::new(HashMap::new())),
            arc_references: Arc::new(RwLock::new(HashMap::new())),
            arc_genotypes: Arc::new(RwLock::new(HashMap::new())),
            arc_annotation_details: Arc::new(RwLock::new(HashMap::new())),
        }
    }
}
