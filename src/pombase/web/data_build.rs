use std::rc::Rc;
use std::collections::BTreeMap;
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::usize;
use std::string::ToString;

use regex::Regex;

use std::collections::{HashMap, HashSet};

use crate::db::raw::*;

use crate::types::*;
use crate::data_types::*;
use crate::web::data::*;
use crate::web::config::*;
use crate::web::util::cmp_str_dates;
use crate::utils::join;

use crate::bio::util::rev_comp;

use flexstr::{SharedStr as FlexStr, shared_str as flex_str, ToSharedStr, shared_fmt as flex_fmt};

use crate::interpro::UniprotResult;
use crate::pfam::PfamProteinDetails;

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

type UniprotIdentifier = FlexStr;

type GenotypeInteractionUniquename = FlexStr;

pub struct WebDataBuild<'a> {
    raw: &'a Raw,
    domain_data: &'a HashMap<UniprotIdentifier, UniprotResult>,
    pfam_data: &'a Option<HashMap<UniprotIdentifier, PfamProteinDetails>>,
    rnacentral_data: &'a Option<RNAcentralAnnotations>,
    config: &'a Config,

    genes: UniquenameGeneMap,
    genotypes: UniquenameGenotypeMap,
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

    // map from term name to term ID (ie "nucleus" -> "GO:0005634")
    term_ids_by_name: HashMap<FlexStr, TermId>,

    genes_of_transcripts: HashMap<FlexStr, FlexStr>,
    transcripts_of_polypeptides: HashMap<FlexStr, FlexStr>,
    parts_of_transcripts: HashMap<FlexStr, Vec<FeatureShort>>,
    genes_of_alleles: HashMap<FlexStr, FlexStr>,
    loci_of_genotypes: HashMap<FlexStr, HashMap<FlexStr, GenotypeLocus>>,
    genotype_display_names: HashMap<GenotypeUniquename, FlexStr>,

    // maps used to collect the genotype interaction part from the feature_relationship table
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

    genetic_interaction_annotations: Vec<GeneticInteractionAnnotation>,
}

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

fn get_feat_rel_expression(feature: &Feature,
                           feature_relationship: &FeatureRelationship) -> Option<FlexStr> {
    for feature_prop in feature.featureprops.borrow().iter() {
        if feature_prop.prop_type.name == "allele_type" {
            if let Some(ref value) = feature_prop.value {
                if value == "deletion" {
                    return Some("Null".into());
                }
            }
        }
    }

    for rel_prop in feature_relationship.feature_relationshipprops.borrow().iter() {
        if rel_prop.prop_type.name == "expression" {
            return rel_prop.value.clone();
        }
    }

    None
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

pub fn compare_ext_part_with_config(extension_relation_order: &RelationOrder,
                                    ep1: &ExtPart, ep2: &ExtPart) -> Ordering {
    let rel_order_conf = extension_relation_order;
    let order_conf = &rel_order_conf.relation_order;
    let always_last_conf = &rel_order_conf.always_last;

    let maybe_ep1_index = order_conf.iter().position(|r| *r == ep1.rel_type_name);
    let maybe_ep2_index = order_conf.iter().position(|r| *r == ep2.rel_type_name);

    if let Some(ep1_index) = maybe_ep1_index {
        if let Some(ep2_index) = maybe_ep2_index {
            ep1_index.cmp(&ep2_index)
        } else {
            Ordering::Less
        }
    } else {
        if maybe_ep2_index.is_some() {
            Ordering::Greater
        } else {
            let maybe_ep1_last_index = always_last_conf.iter().position(|r| *r == ep1.rel_type_name);
            let maybe_ep2_last_index = always_last_conf.iter().position(|r| *r == ep2.rel_type_name);

            if let Some(ep1_last_index) = maybe_ep1_last_index {
                if let Some(ep2_last_index) = maybe_ep2_last_index {
                    ep1_last_index.cmp(&ep2_last_index)
                } else {
                    Ordering::Greater
                }
            } else {
                if maybe_ep2_last_index.is_some() {
                    Ordering::Less
                } else {
                    let name_cmp = ep1.rel_type_name.cmp(&ep2.rel_type_name);

                    if name_cmp == Ordering::Equal {
                        if ep1.ext_range.is_gene() && !ep2.ext_range.is_gene() {
                            Ordering::Less
                        } else {
                            if !ep1.ext_range.is_gene() && ep2.ext_range.is_gene() {
                                Ordering::Greater
                            } else {
                                Ordering::Equal
                            }
                        }
                    } else {
                        name_cmp
                    }
                }
            }
        }
    }
}

lazy_static! {
    static ref BAD_GENOTYPE_NAME_CHARS_RE: Regex =
        Regex::new(r"[% /&;?]").unwrap();
}

pub fn make_genotype_display_name(loci: &[GenotypeLocus],
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
    match feature_locs.get(0) {
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
    let start = (loc.start_pos - 1) as usize;
    let end = loc.end_pos as usize;
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
            for category in &filter.extension_categories {
                add_filter_ancestor(&mut ret, category, cv_name);
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
                } else {
                    if ref2.canto_added_date.is_some() {
                        Ordering::Greater
                    } else {
                        Ordering::Equal
                    }
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
            reference.canto_approved_date.is_some() &&
                reference.canto_curator_role.is_some()
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
        if reference.canto_curator_role == Some("community".into()) {
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

fn add_introns_to_transcript(chromosome: &ChromosomeDetails,
                             transcript_uniquename: &FlexStr,
                             strand: Strand,
                             parts: &mut Vec<FeatureShort>) {
    let mut new_parts: Vec<FeatureShort> = vec![];

    let mut intron_count = 0;

    for part in parts.drain(0..) {
        let mut maybe_new_intron = None;

        if let Some(prev_part) = new_parts.last() {
            let intron_start = prev_part.location.end_pos + 1;
            let intron_end = part.location.start_pos - 1;

            if intron_start > intron_end {
                if intron_start > intron_end + 1 {
                    println!("no gap between exons at {}..{} in {}", intron_start, intron_end,
                             transcript_uniquename);
                }
                // if intron_start == intron_end-1 then it is a one base overlap that
                // represents a frameshift in the reference See:
                // https://github.com/pombase/curation/issues/1453#issuecomment-303214177
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
                    } else {
                        if prev_part.feature_type == FeatureType::FivePrimeUtr {
                            FeatureType::FivePrimeUtrIntron
                        } else {
                            FeatureType::ThreePrimeUtrIntron
                        }
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
            } else {
                if part.location.strand == Strand::Forward {
                    if part.feature_type != FeatureType::FivePrimeUtr {
                        println!("{:?}", parts);
                        panic!("wrong feature type '{}' before exons in {}",
                               part.feature_type, transcript_uniquename);
                    }
                } else {
                    if part.feature_type != FeatureType::ThreePrimeUtr {
                        println!("{:?}", parts);
                        panic!("wrong feature type '{}' after exons in {}",
                               part.feature_type, transcript_uniquename);
                    }
                }
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
            } else {
                if part.location.strand == Strand::Forward {
                    if part.feature_type != FeatureType::ThreePrimeUtr {
                        panic!("wrong feature type '{}' before exons in {}",
                               part.feature_type, transcript_uniquename);
                    }
                } else {
                    if part.feature_type != FeatureType::FivePrimeUtr {
                        panic!("wrong feature type '{}' after exons in {}",
                               part.feature_type, transcript_uniquename);
                    }
                }
            }
        }
    }
}


impl <'a> WebDataBuild<'a> {
    pub fn new(raw: &'a Raw,
               domain_data: &'a HashMap<UniprotIdentifier, UniprotResult>,
               pfam_data: &'a Option<HashMap<UniprotIdentifier, PfamProteinDetails>>,
               rnacentral_data: &'a Option<RNAcentralAnnotations>,
               config: &'a Config) -> WebDataBuild<'a>
    {
        WebDataBuild {
            raw,
            domain_data,
            pfam_data,
            rnacentral_data,
            config,

            genes: BTreeMap::new(),
            genotypes: HashMap::new(),
            genotypes_of_alleles: HashMap::new(),
            genotype_backgrounds: HashMap::new(),
            genotype_interaction_genotype_a: HashMap::new(),
            genotype_interaction_genotype_b: HashMap::new(),
            alleles: HashMap::new(),
            transcripts: HashMap::new(),
            other_features: HashMap::new(),
            terms: HashMap::new(),
            chromosomes: BTreeMap::new(),
            references: HashMap::new(),
            all_ont_annotations: HashMap::new(),
            all_not_ont_annotations: HashMap::new(),
            genetic_interaction_annotations: vec![],
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
            genotype_display_names: HashMap::new(),

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
       }
    }

    fn add_ref_to_hash(&self,
                       seen_references: &mut HashMap<FlexStr, ReferenceShortOptionMap>,
                       identifier: &FlexStr,
                       maybe_reference_uniquename: &Option<ReferenceUniquename>) {
        if let Some(reference_uniquename) = maybe_reference_uniquename {
            if reference_uniquename != "null" {
                seen_references
                    .entry(identifier.into())
                    .or_insert_with(HashMap::new)
                    .insert(reference_uniquename.clone(), None);
            }
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
            .or_insert_with(HashMap::new)
            .insert(other_gene_uniquename.clone(), None);
    }

    fn add_genotype_to_hash(&self,
                            seen_genotypes: &mut HashMap<FlexStr, GenotypeShortMap>,
                            seen_alleles: &mut HashMap<FlexStr, AlleleShortMap>,
                            seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                            identifier: &FlexStr,
                            genotype_uniquename: &FlexStr) {
        let genotype_short = self.make_genotype_short(genotype_uniquename);
        for locus in &genotype_short.loci {
            for expressed_allele in &locus.expressed_alleles {
                self.add_allele_to_hash(seen_alleles, seen_genes, identifier,
                                        &expressed_allele.allele_uniquename);
            }
        }

        seen_genotypes
            .entry(identifier.clone())
            .or_insert_with(HashMap::new)
            .insert(genotype_uniquename.clone(), genotype_short);
    }

    fn add_allele_to_hash(&self,
                          seen_alleles: &mut HashMap<FlexStr, AlleleShortMap>,
                          seen_genes: &mut HashMap<FlexStr, GeneShortOptionMap>,
                          identifier: &FlexStr,
                          allele_uniquename: &AlleleUniquename) -> AlleleShort {
        let allele_short = self.make_allele_short(allele_uniquename);
        {
            let allele_gene_uniquename = &allele_short.gene_uniquename;
            self.add_gene_to_hash(seen_genes, identifier, allele_gene_uniquename);
            seen_alleles
                .entry(identifier.clone())
                .or_insert_with(HashMap::new)
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
                .or_insert_with(HashMap::new)
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
            .or_insert_with(HashMap::new)
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
                let maybe_ortholog_name = ortholog_gene_summary.name.clone();

                IdNameAndOrganism {
                    identifier: ortholog_uniquename,
                    secondary_identifier: maybe_secondary_identifier,
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
        let exon_count =
            if let Some(transcript_uniquename) = gene_details.transcripts.get(0) {

                let transcript = self.transcripts
                    .get(transcript_uniquename)
                    .expect(&format!("internal error, can't find transcript details for {}",
                                     transcript_uniquename));

                let mut count = 0;
                for part in &transcript.parts {
                    if part.feature_type == FeatureType::Exon {
                        count += 1;
                    }
                }
                count
            } else {
                0
            };
        let mut ortholog_taxonids = HashSet::new();
        for ortholog_annotation in &gene_details.ortholog_annotations {
            ortholog_taxonids.insert(ortholog_annotation.ortholog_taxonid);
        }

        let transcript_details = gene_details.transcripts
            .iter()
            .map(|transcript_uniquename| {
                self.transcripts.get(transcript_uniquename)
                    .expect(&format!("internal error, failed to find transcript: {}",
                                    transcript_uniquename))
                    .clone()
            }).collect::<Vec<_>>();

        APIGeneSummary {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
            uniprot_identifier: gene_details.uniprot_identifier.clone(),
            exact_synonyms: synonyms,
            dbxrefs: gene_details.dbxrefs.clone(),
            location: gene_details.location.clone(),
            transcripts: transcript_details,
            tm_domain_count: gene_details.tm_domain_coords.len(),
            coiled_coil_count: gene_details.coiled_coil_coords.len(),
            disordered_regions_count: gene_details.disordered_region_coords.len(),
            low_complexity_regions_count: gene_details.low_complexity_region_coords.len(),
            exon_count,
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

    fn add_gene_product(&mut self, gene_uniquename: &FlexStr, product: &FlexStr) {
        let gene_details = self.get_gene_mut(gene_uniquename);
        gene_details.product = Some(product.clone());
    }

    fn add_name_description(&mut self, gene_uniquename: &FlexStr, name_description: &FlexStr) {
        let gene_details = self.get_gene_mut(gene_uniquename);
        gene_details.name_descriptions.push(name_description.into());
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

    // return the gene of a single allele
    fn gene_from_genotype(&self, genotype_uniquename: &GenotypeUniquename)
       -> GeneUniquename
    {
        let genotype_display_name =
            self.genotype_display_names.get(genotype_uniquename).unwrap();

        let genotype =
            self.genotypes.get(genotype_display_name)
            .expect(&format!("internal error: can't find genotype {}", genotype_uniquename));

        let loci = &genotype.loci;

        if loci.len() > 1 {
            panic!("genotype has multiple loci");
        }

        let expressed_alleles = &loci[0].expressed_alleles;
        if expressed_alleles.len() > 1 {
            panic!("genotype is a diploid");
        }

        let allele_uniquename = &expressed_alleles[0].allele_uniquename;

        let allele = self.alleles.get(allele_uniquename).unwrap();

        allele.gene.uniquename.clone()
    }

    fn add_genetic_interaction(&mut self, genotype_interaction_feature: &Feature,
                               extension_relation_order: &RelationOrder,
                               cvterm: &Cvterm,
                               annotation_template: OntAnnotationDetail) {
        let genotype_interaction_uniquename = &genotype_interaction_feature.uniquename;
        let genotype_a_uniquename =
            self.genotype_interaction_genotype_a.get(genotype_interaction_uniquename)
                .expect(&format!("can't find genotype_a of {}",
                                 genotype_interaction_uniquename));
        let genotype_a_display_name =
            self.genotype_display_names.get(genotype_a_uniquename).unwrap().clone();

        let genotype_b_uniquename =
            self.genotype_interaction_genotype_b.get(genotype_interaction_uniquename)
                .expect(&format!("can't find genotype_b of {}",
                                 genotype_interaction_uniquename));
        let genotype_b_display_name =
            self.genotype_display_names.get(genotype_b_uniquename).unwrap().clone();

        let ont_annotation_detail =
            self.annotation_from_template(extension_relation_order,
                                          cvterm, annotation_template);

        let mut interaction_type = None;
        let mut interaction_note = None;
        let mut rescued_phenotype_termid = None;

        for prop in genotype_interaction_feature.featureprops.borrow().iter() {
            match prop.prop_type.name.as_str() {
                "interaction_type" => interaction_type = prop.value.clone(),
                "interaction_note" => interaction_note = prop.value.clone(),
                "interaction_rescued_phenotype_id" => rescued_phenotype_termid = prop.value.clone(),
                _ => (),
            }
        }

        let double_mutant_phenotype_termid  = Some(cvterm.termid());

        let interaction_type =
            interaction_type.expect(&format!("interaction_type missing for {}",
                                             genotype_interaction_uniquename));

        let interaction_annotation =
            GeneticInteractionAnnotation {
                gene_a_uniquename: self.gene_from_genotype(genotype_a_uniquename),
                gene_b_uniquename: self.gene_from_genotype(genotype_b_uniquename),
                genotype_a_uniquename: Some(genotype_a_display_name),
                genotype_b_uniquename: Some(genotype_b_display_name),
                double_mutant_phenotype_termid,
                rescued_phenotype_termid,
                reference_uniquename: ont_annotation_detail.reference,
                throughput: ont_annotation_detail.throughput,
                interaction_type,
                interaction_note,
            };

        self.genetic_interaction_annotations.push(interaction_annotation);
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

    fn process_references(&mut self) {
        for rc_publication in &self.raw.publications {
            let reference_uniquename = &rc_publication.uniquename;

            if reference_uniquename.to_lowercase() == "null" {
                continue;
            }

            let mut pubmed_authors: Option<FlexStr> = None;
            let mut pubmed_publication_date: Option<FlexStr> = None;
            let mut pubmed_abstract: Option<FlexStr> = None;
            let mut pubmed_doi: Option<FlexStr> = None;
            let mut non_pubmed_authors: Option<FlexStr> = None;
            let mut non_pubmed_abstract: Option<FlexStr> = None;
            let mut canto_annotation_status: Option<FlexStr> = None;
            let mut canto_triage_status: Option<FlexStr> = None;
            let mut canto_curator_role: Option<FlexStr> = None;
            let mut canto_curator_name: Option<FlexStr> = None;
            let mut canto_first_approved_date: Option<FlexStr> = None;
            let mut canto_approved_date: Option<FlexStr> = None;
            let mut canto_added_date: Option<FlexStr> = None;
            let mut canto_session_submitted_date: Option<FlexStr> = None;
            let mut annotation_curators = vec![];

            for prop in rc_publication.publicationprops.borrow().iter() {
                match &prop.prop_type.name as &str {
                    "pubmed_publication_date" =>
                        pubmed_publication_date = Some(prop.value.clone()),
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
                    "canto_annotation_status" =>
                        canto_annotation_status = Some(prop.value.clone()),
                    "canto_triage_status" =>
                        canto_triage_status = Some(prop.value.clone()),
                    "canto_first_approved_date" =>
                        canto_first_approved_date = Some(prop.value.clone()),
                    "canto_approved_date" =>
                        canto_approved_date = Some(prop.value.clone()),
                    "canto_added_date" =>
                        canto_added_date = Some(prop.value.clone()),
                    "canto_session_submitted_date" =>
                        canto_session_submitted_date = Some(prop.value.clone()),
                    "annotation_curator" => {
                        let curator: AnnotationCurator = serde_json::from_str(&prop.value)
                            .expect(&format!("failed to parse annotation_curators pupprop: {}", prop.value));
                        annotation_curators.push(curator);
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

            let mut authors_abbrev = None;
            let mut publication_year = None;

            if let Some(authors) = pubmed_authors.clone() {
                if authors.contains(',') {
                    let author_re = Regex::new(r"^(?P<f>[^,]+),.*$").unwrap();
                    let replaced: String =
                        author_re.replace_all(&authors, "$f et al.").into();
                    authors_abbrev = Some(replaced.to_shared_str());
                } else {
                    authors_abbrev = Some(authors.clone());
                }
            }

            if let Some(publication_date) = pubmed_publication_date.clone() {
                let date_re = Regex::new(r"^(.* )?(?P<y>\d\d\d\d)$").unwrap();
                publication_year = Some(date_re.replace_all(&publication_date, "$y").to_shared_str());
            }

            let mut approved_date = canto_first_approved_date.clone();

            if approved_date.is_none() {
                approved_date = canto_session_submitted_date.clone();
            }

            approved_date =
                if let Some(date) = approved_date {
                    let re = Regex::new(r"^(?P<date>\d\d\d\d-\d\d-\d\d).*").unwrap();
                    Some(re.replace_all(&date, "$date").to_shared_str())
                } else {
                    None
                };

            if let Some(ref canto_annotation_status) = canto_annotation_status {
                if canto_annotation_status != "APPROVED" {
                    approved_date = None;
                }
            }

            let authors = pubmed_authors.or(non_pubmed_authors);

            for annotation_curator in &annotation_curators {
              if annotation_curator.community_curator {
                canto_curator_name = Some(annotation_curator.name.clone());
                canto_curator_role = Some(flex_str!("community"));
              } else {
                if canto_curator_role.is_none() {
                  canto_curator_role = Some(self.config.database_name.clone());
                  canto_curator_name = Some(annotation_curator.name.clone());
                }
              }
            }

            self.references.insert(reference_uniquename.clone(),
                                   ReferenceDetails {
                                       uniquename: reference_uniquename.clone(),
                                       title: rc_publication.title.clone(),
                                       citation: rc_publication.miniref.clone(),
                                       pubmed_abstract: pubmed_abstract.or(non_pubmed_abstract),
                                       pubmed_doi,
                                       authors,
                                       authors_abbrev,
                                       pubmed_publication_date: pubmed_publication_date.clone(),
                                       canto_annotation_status,
                                       canto_triage_status,
                                       canto_curator_role,
                                       canto_curator_name,
                                       canto_first_approved_date,
                                       canto_approved_date,
                                       canto_session_submitted_date,
                                       canto_added_date,
                                       annotation_curators,
                                       approved_date,
                                       publication_year,
                                       cv_annotations: HashMap::new(),
                                       physical_interactions: vec![],
                                       genetic_interactions: vec![],
                                       ortholog_annotations: vec![],
                                       paralog_annotations: vec![],
                                       genes_by_uniquename: HashMap::new(),
                                       genotypes_by_uniquename: HashMap::new(),
                                       alleles_by_uniquename: HashMap::new(),
                                       transcripts_by_uniquename: HashMap::new(),
                                       terms_by_termid: HashMap::new(),
                                       annotation_details: HashMap::new(),
                                       gene_count: 0,
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
                        let expression = get_feat_rel_expression(&feature_rel.subject, feature_rel);
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
                            };

                        let genotype_uniquename = object_uniquename;
                        let genotype_entry = self.loci_of_genotypes.entry(genotype_uniquename.clone());
                        let locus_map = genotype_entry.or_insert_with(HashMap::new);

                        let genotype_locus =
                            locus_map.entry(genotype_locus_identifier)
                            .or_insert_with(|| GenotypeLocus {
                                expressed_alleles: vec![]
                            });

                        genotype_locus.expressed_alleles.push(allele_and_expression);

                        self.genotypes_of_alleles.entry(allele_uniquename.clone())
                           .or_insert_with(HashSet::new)
                           .insert(genotype_uniquename.clone());

                        continue;
                    }
            }
            if TRANSCRIPT_PART_TYPES.contains(&subject_type_name.as_str()) {
                let entry = self.parts_of_transcripts.entry(object_uniquename.clone());
                let part = make_feature_short(&self.chromosomes, &feature_rel.subject);
                entry.or_insert_with(Vec::new).push(part);
            }

            if object_type_name == "genotype_interaction" {
                let genotype_interaction_uniquename = object_uniquename.clone();

                if rel_name == "interaction_genotype_a" {
                    let genotype_a_uniquename = subject_uniquename.clone();
                    self.genotype_interaction_genotype_a.insert(genotype_interaction_uniquename,
                                                                genotype_a_uniquename);
                } else {
                    if rel_name == "interaction_genotype_b" {
                        let genotype_b_uniquename = subject_uniquename.clone();
                        self.genotype_interaction_genotype_b.insert(genotype_interaction_uniquename,
                                                                    genotype_b_uniquename);
                    } else {
                        panic!("unknown relation type {} for interaction {}", rel_name,
                               object_uniquename);
                    }
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
        let maybe_location = make_location(&self.chromosomes, feat);

        if let Some(ref location) = maybe_location {
            if let Some(ref mut chr) = self.chromosomes.get_mut(&location.chromosome_name) {
                chr.gene_uniquenames.push(feat.uniquename.clone());
            }
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
        let mut biogrid_interactor_id: Option<u32> = None;
        let mut rnacentral_urs_identifier = None;

        for prop in feat.featureprops.borrow().iter() {
            match prop.prop_type.name.as_str() {
                "uniprot_identifier" => uniprot_identifier = prop.value.clone(),
                "sgd_identifier" => secondary_identifier = prop.value.clone(),
                "biogrid_interactor_id" => {
                    if let Some(ref chado_biogrid_id) = prop.value {
                        biogrid_interactor_id = match chado_biogrid_id.parse::<u32>() {
                            Ok(val) => Some(val),
                            Err(err) =>
                                panic!("error parsing BioGRID interactor ID from Chado: {}", err),
                        }
                    }
                },
                "rnacentral_identifier" => rnacentral_urs_identifier = prop.value.clone(),
                _ => (),
            }
        }

        let (interpro_matches, tm_domain_coords) =
            if let Some(ref uniprot_identifier) = uniprot_identifier {
                if let Some(result) = self.domain_data.get(uniprot_identifier) {
                    let tm_domain_matches = result.tmhmm_matches.iter()
                        .map(|tm_match| (tm_match.start, tm_match.end))
                        .collect::<Vec<_>>();
                    (result.interpro_matches.clone(), tm_domain_matches)
                } else {
                    (vec![], vec![])
                }
            } else {
                (vec![], vec![])
            };

        let (disordered_region_coords, low_complexity_region_coords, coiled_coil_coords) =
            if let Some(pfam_data) = self.pfam_data {
                if let Some(ref uniprot_identifier) = uniprot_identifier {
                    if let Some(result) = pfam_data.get(uniprot_identifier) {
                        let mut disordered_region_coords = vec![];
                        let mut low_complexity_region_coords = vec![];
                        let mut coiled_coil_coords = vec![];
                        for motif in &result.motifs {
                            match &motif.motif_type as &str {
                                "disorder" =>
                                    disordered_region_coords.push((motif.start, motif.end)),
                                "low_complexity" =>
                                    low_complexity_region_coords.push((motif.start, motif.end)),
                                "coiled_coil" =>
                                    coiled_coil_coords.push((motif.start, motif.end)),
                                _ => (),
                            }
                        }
                        (disordered_region_coords, low_complexity_region_coords, coiled_coil_coords)
                    } else {
                        (vec![], vec![], vec![])
                    }
                } else {
                    (vec![], vec![], vec![])
                }
            } else {
                (vec![], vec![], vec![])
            };

        let rfam_annotations =
            if let Some(rnacentral_data) = self.rnacentral_data {
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

        let gene_feature = GeneDetails {
            uniquename: feat.uniquename.clone(),
            name: feat.name.clone(),
            taxonid: organism.taxonid,
            product: None,
            deletion_viability: DeletionViability::Unknown,
            uniprot_identifier,
            secondary_identifier,
            biogrid_interactor_id,
            rnacentral_urs_identifier,
            interpro_matches,
            tm_domain_coords,
            disordered_region_coords,
            low_complexity_region_coords,
            coiled_coil_coords,
            rfam_annotations,
            orfeome_identifier,
            name_descriptions: vec![],
            synonyms: vec![],
            dbxrefs,
            feature_type: feat.feat_type.name.clone(),
            feature_so_termid: feat.feat_type.termid(),
            transcript_so_termid: feat.feat_type.termid(),
            characterisation_status: None,
            taxonomic_distribution: None,
            location: maybe_location,
            gene_neighbourhood: vec![],
            cv_annotations: HashMap::new(),
            physical_interactions: vec![],
            genetic_interactions: vec![],
            ortholog_annotations: vec![],
            paralog_annotations: vec![],
            target_of_annotations: vec![],
            transcripts: vec![],
            transcripts_by_uniquename: HashMap::new(),
            genes_by_uniquename: HashMap::new(),
            genotypes_by_uniquename: HashMap::new(),
            alleles_by_uniquename: HashMap::new(),
            references_by_uniquename: HashMap::new(),
            terms_by_termid: HashMap::new(),
            annotation_details: HashMap::new(),
            feature_publications: HashSet::new(),
            subset_termids: HashSet::new(),
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
                add_introns_to_transcript(chromosome, transcript_uniquename, strand, &mut parts);
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

        for part in &parts {
            if part.location.start_pos < transcript_start {
                transcript_start = part.location.start_pos;
            }
            if part.location.end_pos > transcript_end {
                transcript_end = part.location.end_pos;
            }
        }

        // use the first part as a template to get the chromosome details
        let transcript_location =
            ChromosomeLocation {
                start_pos: transcript_start,
                end_pos: transcript_end,
                phase: None,
                .. parts[0].location.clone()
            };

        let maybe_cds_location =
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
                    None
                } else {
                    if let Some(mrna_location) = feat.featurelocs.borrow().get(0) {
                        let first_part_loc = &parts[0].location;

                        Some(ChromosomeLocation {
                            chromosome_name: first_part_loc.chromosome_name.clone(),
                            start_pos: cds_start,
                            end_pos: cds_end,
                            strand: first_part_loc.strand,
                            phase: make_phase(mrna_location),
                        })
                    } else {
                        None
                    }
                }
            } else {
                None
            };

        if let Some(gene_uniquename) =
            self.genes_of_transcripts.get(&transcript_uniquename) {
                let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
                let transcript_type = feat.feat_type.name.clone();
                if gene_details.feature_type == "gene" {
                    let feature_type = format!("{} {}", transcript_type, gene_details.feature_type);
                    gene_details.feature_type = feature_type.to_shared_str();
                }
                let transcript = TranscriptDetails {
                    uniquename: transcript_uniquename.clone(),
                    location: transcript_location,
                    transcript_type,
                    parts,
                    protein: None,
                    cds_location: maybe_cds_location,
                    gene_uniquename: gene_uniquename.to_owned(),
                };

                self.transcripts.insert(transcript_uniquename.clone(), transcript.clone());

                gene_details.transcripts.push(transcript_uniquename);
                gene_details.transcript_so_termid = feat.feat_type.termid();
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
                if let Some(ref prop_value) = p {
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
                if prop.prop_type.name == "molecular_weight" {
                    if let Some(value) = parse_prop_as_f32(&prop.value) {
                        molecular_weight = Some(value / 1000.0);
                    }
                }
                if prop.prop_type.name == "average_residue_weight" {
                    if let Some(value) = parse_prop_as_f32(&prop.value) {
                        average_residue_weight = Some(value / 1000.0);
                    }
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

            let protein = ProteinDetails {
                uniquename: feat.uniquename.clone(),
                sequence: residues.to_shared_str(),
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
                        .expect(&format!("internal error, failed to find transcript: {}",
                                         transcript_uniquename))
                        .protein = Some(protein);
                } else {
                    panic!("can't find transcript of polypeptide: {}", protein_uniquename)
                }
        } else {
            panic!("no residues for protein: {}", feat.uniquename);
        }
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
        let genotype_display_uniquename =
            make_genotype_display_name(&loci, &self.alleles);

        self.genotype_display_names.insert(genotype_uniquename.clone(),
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
            } else {
                if prop.prop_type.name == "genotype_comment" {
                    if let Some(ref comment_ref) = prop.value {
                        comment = Some(comment_ref.to_shared_str());
                    }
                }
            }
        }

        let rc_display_name = genotype_display_uniquename.to_shared_str();

        self.genotypes.insert(rc_display_name.clone(),
                              GenotypeDetails {
                                  display_uniquename: rc_display_name,
                                  name: feat.name.as_ref().map(|s| s.to_shared_str()),
                                  taxonid: organism.taxonid,
                                  loci,
                                  ploidiness,
                                  comment,
                                  cv_annotations: HashMap::new(),
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
                      reference_uniquename:
                         prop.featureprop_pubs.borrow().get(0)
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
            if !TRANSCRIPT_FEATURE_TYPES.contains(&feat.feat_type.name.as_str()) &&
                !TRANSCRIPT_PART_TYPES.contains(&feat.feat_type.name.as_str()) &&
                !HANDLED_FEATURE_TYPES.contains(&feat.feat_type.name.as_str())
            {
                // for now, ignore features without locations
                if feat.featurelocs.borrow().len() > 0 {
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
            let rel_term_name =
                self.make_term_short(&rel_termid).name;

            if self.is_interesting_parent(&object_termid, &rel_term_name) {
                interesting_parents_by_termid
                    .entry(subject_termid.clone())
                    .or_insert_with(HashSet::new)
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
                        self.genotype_display_names.get(genotype_uniquename).unwrap();

                    if let Some(genotype_details) = self.genotypes.get(genotype_display_uniquename) {
                        for term_annotations in genotype_details.cv_annotations.values() {
                            for term_annotation in term_annotations {
                                let termid = &term_annotation.term;
                                let term_details =
                                     self.terms.get(termid).expect(&format!("missing termid {}", termid));
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
            if let Some(ref location) = gene_details.location {
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
                    if i > GENE_NEIGHBOURHOOD_DISTANCE {
                        i - GENE_NEIGHBOURHOOD_DISTANCE
                    } else {
                        0
                    };

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


    // add physical interaction, legacy genetic interaction, ortholog and paralog annotations
    fn process_annotation_feature_rels(&mut self) {
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
                        let mut ortholog_qualifier = None;

                        let borrowed_publications = feature_rel.publications.borrow();
                        let maybe_publication = borrowed_publications.get(0);
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
                            if prop.prop_type.name == "evidence" {
                                if let Some(ref evidence_long) = prop.value {
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
                            if prop.prop_type.name == "is_inferred" {
                                if let Some(is_inferred_value) = prop.value.clone() {
                                    if is_inferred_value == "yes" {
                                        is_inferred_interaction = true;
                                    }
                                }
                            }
                            if prop.prop_type.name == "annotation_throughput_type" {
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
                            }
                            if prop.prop_type.name == "interaction_note" {
                                if let Some(interaction_note_value) = prop.value.clone() {
                                    interaction_note = Some(interaction_note_value);
                                }
                            }
                            if prop.prop_type.name == "ortholog qualifier" ||
                                prop.prop_type.name == "ortholog_qualifier" {
                                ortholog_qualifier = prop.value.clone()
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
                                        evidence.expect(&format!("evidence missing for feature_relationship_id: {}",
                                                                 feature_rel.feature_relationship_id));
                                    let interaction_annotation =
                                        InteractionAnnotation {
                                            gene_uniquename: gene_uniquename.clone(),
                                            interactor_uniquename: other_gene_uniquename.clone(),
                                            evidence,
                                            reference_uniquename: maybe_reference_uniquename.clone(),
                                            throughput,
                                            interaction_note,
                                        };
                                    if rel_name == "interacts_genetically" {
                                        let interaction_annotation = &interaction_annotation;
                                        self.genetic_interaction_annotations.push(interaction_annotation.into());
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
                                let gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                                gene_details.ortholog_annotations.push(ortholog_annotation.clone());
                                if let Some(ref_details) =
                                    if let Some(ref reference_uniquename) = maybe_reference_uniquename {
                                        self.references.get_mut(reference_uniquename)
                                    } else {
                                        None
                                    }
                                {
                                    ref_details.ortholog_annotations.push(ortholog_annotation);
                                }
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
                                {
                                    if self.config.load_organism_taxonid.is_some() &&
                                        self.config.load_organism_taxonid.unwrap() == gene_details.taxonid ||
                                        gene_organism_taxonid < other_gene_organism_taxonid
                                    {
                                        ref_details.paralog_annotations.push(paralog_annotation);
                                    }
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
                                {
                                    if self.config.load_organism_taxonid.is_some() &&
                                        self.config.load_organism_taxonid.unwrap() == other_gene_details.taxonid ||
                                        gene_organism_taxonid > other_gene_organism_taxonid
                                    {
                                        ref_details.paralog_annotations.push(paralog_annotation);
                                    }
                                }
                            },
                        }
                    }
            }
        }

        for ref_details in self.references.values_mut() {
            ref_details.physical_interactions.sort();
            ref_details.genetic_interactions.sort();
            ref_details.ortholog_annotations.sort();
            ref_details.paralog_annotations.sort();
        }

        for gene_details in self.genes.values_mut() {
            gene_details.physical_interactions.sort();
            gene_details.genetic_interactions.sort();
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
                              genes: &[FlexStr],
                              maybe_genotype_uniquename: &Option<FlexStr>,
                              reference_uniquename: &Option<FlexStr>,
                              annotation_termid: &FlexStr,
                              extension: &[ExtPart]) -> Vec<(GeneUniquename, TargetOfAnnotation)> {
        if genes.len() != 1 {
            panic!("expected an annotation with one gene for {}, got: {:?}",
                   annotation_termid, genes);
        }
        let gene = &genes[0];
        let mut ret_vec = vec![];

        for ext_part in extension {
            let maybe_ext_config =
                self.matching_ext_config(annotation_termid, &ext_part.rel_type_name);
            if let ExtRange::Gene(ref target_gene_uniquename) = ext_part.ext_range {
                if let Some(ext_config) = maybe_ext_config {
                    if let Some(reciprocal_display_name) =
                        ext_config.reciprocal_display {
                            let (annotation_gene_uniquename, annotation_genotype_uniquename) =
                                if maybe_genotype_uniquename.is_some() {
                                    (gene.clone(), maybe_genotype_uniquename.clone())
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
                                          }));
                        }
                }
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

            let existing_rel = seen_gene_rels.get(&annotation.gene);

            if let Some(existing_rel) = existing_rel {
                if *existing_rel > rel_priority {
                    annotation.show_in_summary = false;
                    continue;
                }
            }
            seen_gene_rels.insert(annotation.gene.clone(), rel_priority);
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
                                                        &annotation.genes,
                                                        &annotation.genotype,
                                                        &annotation.reference,
                                                        &term_details.termid,
                                                        &annotation.extension);
                        for (target_gene_uniquename, new_annotation) in new_annotations {
                           if self.genes.get(&target_gene_uniquename).is_some() {
                               target_of_annotations
                                   .entry(target_gene_uniquename.clone())
                                   .or_insert_with(HashSet::new)
                                   .insert(new_annotation);
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
                                } else {
                                    if interesting_parent_ids.contains(inviable_termid) ||
                                        *inviable_termid == term_annotation.term {
                                            inviable_conditions.insert(conditions_as_string,
                                                                       term_annotation.term.clone());
                                        }
                                }
                        }
                    }

                    if viable_conditions.is_empty() {
                        if !inviable_conditions.is_empty() {
                            new_status = DeletionViability::Inviable;
                        }
                    } else {
                        if inviable_conditions.is_empty() {
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
                }
            gene_statuses.insert(gene_uniquename.clone(), new_status);
        }

        for (gene_uniquename, status) in &gene_statuses {
            if let Some(ref mut gene_details) = self.genes.get_mut(gene_uniquename) {
                gene_details.deletion_viability = status.clone();
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

            if let Some(ref characterisation_status) = gene_details.characterisation_status {
                if characterisation_status == "dubious" {
                    gene_details.taxonomic_distribution = Some(flex_str!("dubious"));
                    continue 'GENE;
                }
            }

            if gene_details.feature_type != "mRNA gene" {
                gene_details.taxonomic_distribution = Some(flex_str!("not curated"));
                continue 'GENE;
            }

            gene_details.taxonomic_distribution = Some(flex_str!("other"));
        }
    }

    fn process_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name != POMBASE_ANN_EXT_TERM_CV_NAME {
                let cv_config = self.config.cv_config_by_name(&cvterm.cv.name);
                let annotation_feature_type = cv_config.feature_type.clone();

                let mut maybe_pombase_gene_id = None;

                for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                  if cvtermprop.prop_type.name == "pombase_gene_id" {
                    maybe_pombase_gene_id = Some(cvtermprop.value.clone());
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
                        } else {
                            if term_xref_id_prop == "ACCESSION" {
                                let dbxref: &Dbxref = cvterm.dbxref.borrow();
                                maybe_xref_id = Some(dbxref.accession.clone());
                            }
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
                            synonym_type: (*syn).synonym_type.name.clone(),
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
                                      gene_count: 0,
                                      genotype_count: 0,
                                      xrefs,
                                      pombase_gene_id: maybe_pombase_gene_id,
                                  });
                self.term_ids_by_name.insert(cvterm.name.clone(), cvterm.termid());
            }
        }
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
                    if (*cvtermprop).prop_type.name.starts_with(ANNOTATION_EXT_REL_PREFIX) {
                        let ext_rel_name_str =
                            &(*cvtermprop).prop_type.name[ANNOTATION_EXT_REL_PREFIX.len()..];
                        let ext_rel_name = ext_rel_name_str.to_shared_str();
                        let ext_range = (*cvtermprop).value.clone();
                        let range: ExtRange = if ext_range.starts_with(&db_prefix) {
                            let db_feature_uniquename = ext_range[db_prefix.len()..].to_shared_str();
                            if let Some(captures) = PROMOTER_RE.captures(&db_feature_uniquename) {
                                let gene_uniquename = captures["gene"].to_shared_str();
                                if self.genes.contains_key(&gene_uniquename) {
                                    ExtRange::Promoter(gene_uniquename)
                                } else {
                                    panic!("unknown gene in promoter: {}", db_feature_uniquename);
                                }
                            } else {
                                if self.genes.contains_key(&db_feature_uniquename) {
                                    ExtRange::Gene(db_feature_uniquename.clone())
                                } else {
                                    if let Some(captures) = TRANSCRIPT_ID_RE.captures(db_feature_uniquename.as_ref()) {
                                        if self.genes.contains_key(&captures["gene"].to_shared_str()) {
                                            ExtRange::Transcript(db_feature_uniquename.clone())
                                        } else {
                                            panic!("unknown gene for transcript: {}", db_feature_uniquename);
                                        }
                                    } else {
                                        panic!("can't find gene or transcript for: {}", db_feature_uniquename);
                                    }
                                }
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
                                    .or_insert_with(Vec::new).push(ExtPart {
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
                                .or_insert_with(Vec::new).push(ExtPart {
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
            } else {
                if let Some(ref mut allele) = self.alleles.get_mut(&feature.uniquename) {
                    let synonym = make_synonym();
                    if let Err(insert_pos) = allele.synonyms.binary_search(&synonym) {
                        // keep synonym list ordered
                        allele.synonyms.insert(insert_pos, synonym);
                    }
                }
            }
        }
    }

    fn process_feature_publications(&mut self) {
        for feature_pub in &self.raw.feature_pubs {
            let feature = &feature_pub.feature;
            let publication = &feature_pub.publication;

            if publication.uniquename.starts_with("PMID:") {
                if let Some(ref mut gene_details) = self.genes.get_mut(&feature.uniquename) {
                    gene_details.feature_publications.insert(publication.uniquename.clone());
                }
            }
        }
    }

    fn make_genotype_short(&self, genotype_display_name: &FlexStr) -> GenotypeShort {
        if let Some(details) = self.genotypes.get(genotype_display_name) {
            GenotypeShort {
                display_uniquename: details.display_uniquename.clone(),
                name: details.name.clone(),
                loci: details.loci.clone(),
            }
        } else {
            panic!("can't find genotype {}", genotype_display_name);
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
        {
            if let Some(ref mut protein) = transcript_details
                .protein
            {
                protein.product = Some(product);
            }
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
                                         Some(transcript_uniquename.clone())
)                                    } else {
                                        (None, None)
                                    }
                            } else {
                                (None, None)
                            }
                    } else {
                        if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                            if let Some(gene_uniquename) =
                                self.genes_of_transcripts.get(&feature.uniquename) {
                                    (Some(gene_uniquename.clone()), Some(feature.uniquename.clone()))
                                } else {
                                    (None, None)
                                }
                        } else {
                            if feature.feat_type.name == "gene" {
                                (Some(feature.uniquename.clone()), None)
                            } else {
                                (None, None)
                            }
                        }
                    }
                } else {
                    (None, None)
                };

            if let Some(gene_uniquename) = maybe_gene_uniquename {
                if let Some(transcript_uniquename) = maybe_transcript_uniquename {
                    if transcript_uniquename.ends_with(".1") {
                        // for multi-transcript genes, use the product
                        // from the first transcript
                        self.add_gene_product(&gene_uniquename, &cvterm.name);
                    }

                    self.add_product_to_protein(&transcript_uniquename,
                                                cvterm.name.clone());
                }
            }

            if feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene" {
                if cvterm.cv.name == "PomBase gene characterisation status" {
                    self.add_characterisation_status(&feature.uniquename, &cvterm.name);
                } else {
                    if cvterm.cv.name == "name_description" {
                        self.add_name_description(&feature.uniquename, &cvterm.name);
                    }
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
                } else {
                    if let Some(name) = &gene_short.name {
                        return WithFromValue::IdentifierAndName({
                            IdentifierAndName {
                                identifier: with_or_from_value.clone(),
                                name: name.clone(),
                            }
                        });
                    }
                }
            } else {
                if self.transcripts.contains_key(&id) {
                    if self.config.database_name == prefix {
                        return WithFromValue::Transcript(id);
                    }
                }
            }
        } else {
            if self.genes.contains_key(with_or_from_value) {
                let gene_short = self.make_gene_short(with_or_from_value);
                // a gene from the main organism
                return WithFromValue::Gene(gene_short);
            } else {
                if self.transcripts.contains_key(with_or_from_value) {
                    return WithFromValue::Transcript(with_or_from_value.clone());
                }
            }
        }

        if self.terms.get(with_or_from_value).is_some() {
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

            let termid = cvterm.termid();

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
            let mut withs: HashSet<WithFromValue> = HashSet::new();
            let mut froms: HashSet<WithFromValue> = HashSet::new();
            let mut qualifiers: Vec<Qualifier> = vec![];
            let mut date: Option<FlexStr> = None;
            let mut assigned_by: Option<FlexStr> = None;
            let mut evidence: Option<FlexStr> = None;
            let mut eco_evidence: Option<FlexStr> = None;
            let mut genotype_background: Option<FlexStr> = None;
            let mut throughput: Option<Throughput> = None;

            // need to get evidence first as it's used later
            // See: https://github.com/pombase/website/issues/455
            for prop in feature_cvterm.feature_cvtermprops.borrow().iter() {
                if &prop.type_name() == "evidence" {
                    if let Some(ref evidence_long) = prop.value {
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
                                let display_name =
                                    self.get_ext_rel_display_name(&termid,
                                                                  &flex_str!("modified residue"));

                                let residue_range_part = ExtPart {
                                    rel_type_id: None,
                                    rel_type_name: display_name.clone(),
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
                        let genotype_display_name =
                            make_genotype_display_name(&loci, &self.alleles);
                        maybe_genotype_uniquename = Some(genotype_display_name);
                        genotype_background =
                            self.genotype_backgrounds.get(&feature.uniquename)
                            .cloned();
                        loci.iter()
                            .map(|locus| {
                                locus.expressed_alleles.iter()
                                    .map(|expressed_allele| {
                                        let allele_short =
                                            self.make_allele_short(&expressed_allele.allele_uniquename);
                                        allele_short.gene_uniquename
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
                    _ =>
                        if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                            if let Some(gene_uniquename) =
                                self.genes_of_transcripts.get(&feature.uniquename) {
                                    if let Some(gene_details) = self.genes.get(gene_uniquename) {
                                        if gene_details.transcripts.len() > 1 {
                                            // only bother to record the specific transcript if
                                            // there is more than one
                                            transcript_uniquenames.push(feature.uniquename.clone());
                                        }
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
                genotype: maybe_genotype_uniquename,
                genotype_background,
                withs,
                froms,
                residue: extra_props_clone.remove(&flex_str!("residue")),
                gene_product_form_id: extra_props_clone.remove(&flex_str!("gene_product_form_id")),
                gene_ex_props,
                qualifiers,
                evidence,
                eco_evidence,
                conditions,
                extension,
                date,
                assigned_by,
                throughput,
            };

            if &feature.feat_type.name == "genotype_interaction" {
                self.add_genetic_interaction(&feature, &rel_order,
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
                        } else {
                            if !multi_locus.annotations.contains(annotation_id) {
                                multi_locus.annotations.push(*annotation_id);
                            }
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

        for (_, annotations) in ont_annotation_map {
            let (no_transcript_annotations, mut has_transcript_annotations): (Vec<i32>, Vec<i32>) =
                annotations
                .iter()
                .partition(|&annotation_id| {
                    if let Some(ont_annotation_detail) =
                        self.annotation_details.get(annotation_id) {
                            ont_annotation_detail.transcript_uniquenames.len() == 0
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
                        if let Some(ref annotation_details) = self.annotation_details.get(&prev_annotation_id) {
                            if !annotation_details.transcript_uniquenames.contains(&current_transcript_uniquename) {
                                self.annotation_details.get_mut(&prev_annotation_id).unwrap()
                                    .transcript_uniquenames.push(current_transcript_uniquename);
                            }
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
                    self.make_term_annotations(termid, &annotations, is_not);

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
                    get(&annotation_id).expect("can't find OntAnnotationDetail");

                for gene_uniquename in &annotation.genes {
                    gene_annotation_by_term.entry(gene_uniquename.clone())
                        .or_insert_with(HashMap::new)
                        .entry(termid.clone())
                        .or_insert_with(Vec::new)
                        .push(*annotation_id);
                }

                if let Some(ref genotype_uniquename) = annotation.genotype {
                    let existing =
                        genotype_annotation_by_term.entry(genotype_uniquename.clone())
                        .or_insert_with(HashMap::new)
                        .entry(termid.clone())
                        .or_insert_with(Vec::new);
                    if !existing.contains(&annotation_id) {
                        existing.push(*annotation_id);
                    }
                }

                if let Some(reference_uniquename) = annotation.reference.clone() {
                    ref_annotation_by_term.entry(reference_uniquename)
                        .or_insert_with(HashMap::new)
                        .entry(termid.clone())
                        .or_insert_with(Vec::new)
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
                        .or_insert_with(Vec::new)
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
                        .or_insert_with(Vec::new)
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
                    ref_details.cv_annotations.entry(cv_name).or_insert_with(Vec::new)
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
                        .or_insert_with(HashSet::new)
                        .insert(subject_termid.clone());
                }

                for (cv_name, term_annotations) in &subject_term_details.cv_annotations {
                    for term_annotation in term_annotations {
                        for annotation_id in &term_annotation.annotations {
                            let dest_termid = object_termid.clone();
                            let source_termid = subject_termid.clone();

                            if !term_annotation.is_not {
                                new_annotations.entry((cv_name.clone(), dest_termid))
                                    .or_insert_with(HashMap::new)
                                    .entry(source_termid)
                                    .or_insert_with(HashMap::new)
                                    .entry(*annotation_id)
                                    .or_insert_with(HashSet::new)
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

        for chadoprop in &self.raw.chadoprops {
            if chadoprop.prop_type.name == "db_creation_datetime" {
                db_creation_datetime = chadoprop.value.clone();
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

        Metadata {
            export_prog_name: flex_str!(PKG_NAME),
            export_prog_version: flex_str!(VERSION),
            db_creation_datetime: db_creation_datetime.unwrap(),
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
            if let Some(transcript) = self.transcripts.get(transcript_uniquename) {
                if let Some(ref protein) = transcript.protein {
                    molecular_weight = Some((100.0 * protein.molecular_weight).round() / 100.0);
                    if protein.sequence.ends_with('*') {
                        protein_length = Some(protein.sequence.len() - 1);
                    } else {
                        protein_length = Some(protein.sequence.len());
                    }
                    break;
                }
            }
        }

        for field_name in &self.config.gene_results.visualisation_field_names {
            let column_conf = &self.config.gene_results.field_config[field_name];
            for attr_value_conf in &column_conf.attr_values {
                if let (Some(ref bin_start), Some(ref bin_end)) =
                    (attr_value_conf.bin_start, attr_value_conf.bin_end) {
                        if let Some(prot_len) = protein_length {
                            if *bin_start <= prot_len && *bin_end >= prot_len {
                                return (molecular_weight,
                                        Some(prot_len),
                                        Some(attr_value_conf.name.clone()));
                            }
                        }
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

    fn make_gene_query_data_map(&self) -> HashMap<GeneUniquename, GeneQueryData> {
        let mut gene_query_data_map = HashMap::new();

        for gene_details in self.genes.values() {
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
                reference_uniquenames,
                subset_termids: gene_details.subset_termids.clone(),
            };

            gene_query_data_map.insert(gene_details.uniquename.clone(), gene_query_data);
        }

        gene_query_data_map
    }

    pub fn make_api_maps(mut self) -> APIMaps {
        let mut gene_summaries: HashMap<GeneUniquename, APIGeneSummary> = HashMap::new();
        let mut gene_name_gene_map = HashMap::new();
        let mut interactors_of_genes = HashMap::new();

        for (gene_uniquename, gene_details) in &self.genes {
            if self.config.load_organism_taxonid.is_none() ||
                self.config.load_organism_taxonid.unwrap() == gene_details.taxonid {
                let gene_summary = self.make_api_gene_summary(gene_uniquename);
                if let Some(ref gene_name) = gene_summary.name {
                    gene_name_gene_map.insert(gene_name.clone(), gene_uniquename.clone());
                }
                gene_summaries.insert(gene_uniquename.clone(), gene_summary);

                let mut interactors = vec![];

                for interaction_annotation in &gene_details.physical_interactions {
                    let interactor_uniquename =
                        if gene_uniquename == &interaction_annotation.gene_uniquename {
                            interaction_annotation.interactor_uniquename.clone()
                        } else {
                            interaction_annotation.gene_uniquename.clone()
                        };
                    let interactor = APIInteractor {
                        interaction_type: InteractionType::Physical,
                        interactor_uniquename,
                    };
                    if !interactors.contains(&interactor) {
                        interactors.push(interactor);
                    }
                }
                for interaction_annotation in &gene_details.genetic_interactions {
                    let interactor_uniquename =
                        if gene_uniquename == &interaction_annotation.gene_a_uniquename {
                            interaction_annotation.gene_b_uniquename.clone()
                        } else {
                            interaction_annotation.gene_a_uniquename.clone()
                        };
                    let interactor = APIInteractor {
                        interaction_type: InteractionType::Genetic,
                        interactor_uniquename,
                    };
                    if !interactors.contains(&interactor) {
                        interactors.push(interactor);
                    }
                }
                interactors_of_genes.insert(gene_uniquename.clone(), interactors);
            }
        }

        let gene_query_data_map = self.make_gene_query_data_map();

        let mut term_summaries: HashSet<TermShort> = HashSet::new();
        let mut termid_genes: HashMap<TermId, HashSet<GeneUniquename>> = HashMap::new();

        let mut terms_for_api: HashMap<TermId, TermDetails> = HashMap::new();

        for termid in self.terms.keys() {
            term_summaries.insert(self.make_term_short(termid));
        }

        let termid_genotype_annotation: HashMap<TermId, Vec<APIGenotypeAnnotation>> =
            self.get_api_genotype_annotation();

        for (termid, term_details) in self.terms.drain() {
            let cv_config = &self.config.cv_config;
            if let Some(term_config) = cv_config.get(&term_details.cv_name) {
                if term_config.feature_type == "gene" {
                    termid_genes.insert(termid.clone(),
                                        term_details.annotated_genes.clone());
                }
            }

            terms_for_api.insert(termid.clone(), term_details);
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
            termid_genotype_annotation,
            term_summaries,
            genes: self.genes,
            gene_name_gene_map,
            transcripts: self.transcripts,
            alleles: self.alleles,
            genotypes: self.genotypes,
            terms: terms_for_api,
            interactors_of_genes,
            references: self.references,
            other_features: self.other_features,
            seq_feature_page_features,
            annotation_details: self.annotation_details,
            chromosomes: self.chromosomes,
            term_subsets,
            gene_subsets,
            children_by_termid,
            gene_expression_measurements,
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
                        annotation_detail.gene_product_form_id {
                            if gene_product_form_id.starts_with("PR:") {
                                self.add_term_to_hash(seen_terms, identifier,
                                                      gene_product_form_id);
                            }
                        }
                    for ext_part in &annotation_detail.extension {
                        match ext_part.ext_range {
                            ExtRange::Term(ref range_termid) |
                            ExtRange::GeneProduct(ref range_termid) =>
                                self.add_term_to_hash(seen_terms, identifier,
                                                      range_termid),
                            ExtRange::GeneAndGeneProduct(GeneAndGeneProduct { ref gene_uniquename, ref product }) => {
                                self.add_gene_to_hash(seen_genes, identifier, gene_uniquename);
                                self.add_term_to_hash(seen_terms, identifier, product);
                            },
                            ExtRange::Gene(ref gene_uniquename) |
                            ExtRange::Promoter(ref gene_uniquename) =>
                                self.add_gene_to_hash(seen_genes, identifier,
                                                      gene_uniquename),
                            ExtRange::Transcript(ref transcript_uniquename) =>
                                self.add_transcript_to_hashes(seen_transcripts, seen_genes,
                                                            identifier,
                                                            transcript_uniquename),
                            _ => {},
                        }
                    }
                    if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  identifier, genotype_uniquename);
                    }

                    let with_from_iter = annotation_detail.withs
                        .iter()
                        .chain(annotation_detail.froms.iter());

                    for with_from_value in with_from_iter {
                        match with_from_value {
                            WithFromValue::Gene(ref gene_short) => {
                                self.add_gene_to_hash(seen_genes, identifier,
                                                      &gene_short.uniquename)
                            },
                            &WithFromValue::Transcript(ref transcript_uniquename) => {
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

    fn set_term_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_transcripts, mut seen_terms) = get_maps();

        let mut annotated_genes_map: HashMap<TermId, HashSet<GeneUniquename>> =
            HashMap::new();
        let mut single_locus_annotated_genes_map: HashMap<TermId, HashSet<GeneUniquename>> =
            HashMap::new();
        let mut multi_locus_annotated_genes_map: HashMap<TermId, HashSet<GeneUniquename>> =
            HashMap::new();

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
                                        .entry(termid.clone()).or_insert_with(HashSet::new)
                                        .insert(gene_uniquename.clone());
                                    }
                                }
                            }

                        if let Some(genotype_uniquename) = annotation_detail.genotype.as_ref() {
                            let genotype_details = self.genotypes.get(genotype_uniquename).unwrap();
                            let gene_uniquenames = self.gene_uniquenames_from_genotype(genotype_details);
                            if genotype_details.loci.len() == 1 {
                                single_locus_annotated_genes_map
                                    .entry(termid.clone()).or_insert_with(HashSet::new)
                                    .extend(gene_uniquenames.iter().cloned());
                            } else {
                                multi_locus_annotated_genes_map
                                    .entry(termid.clone()).or_insert_with(HashSet::new)
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
                            annotation_detail.gene_product_form_id {
                                if gene_product_form_id.starts_with("PR:") {
                                    self.add_term_to_hash(&mut seen_terms, termid,
                                                          gene_product_form_id);
                                }
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
                                                      &mut seen_genes, termid,
                                                      genotype_uniquename);
                        }

                        let with_from_iter = annotation_detail.withs
                            .iter()
                            .chain(annotation_detail.froms.iter());

                        for with_from_value in with_from_iter {
                            match with_from_value {
                                WithFromValue::Gene(ref gene_short) => {
                                    self.add_gene_to_hash(&mut seen_genes, termid,
                                                          &gene_short.uniquename)
                                },
                                &WithFromValue::Transcript(ref transcript_uniquename) => {
                                    self.add_transcript_to_hashes(&mut seen_transcripts, &mut seen_genes,
                                                                  termid,  transcript_uniquename);
                                },
                                _ => (),
                            }
                        }
                    }
                }
            }
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
                for interaction in interaction_iter {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &interaction.reference_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &interaction.gene_a_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &interaction.gene_b_uniquename);
                    if let Some(ref genotype_a_uniquename) = interaction.genotype_a_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
                                                  gene_uniquename,
                                                  genotype_a_uniquename);
                    }
                    if let Some(ref genotype_b_uniquename) = interaction.genotype_b_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
                                                  gene_uniquename,
                                                  genotype_b_uniquename);
                    }

                    if let Some(ref double_mutant_phenotype_termid) = interaction.double_mutant_phenotype_termid {
                        self.add_term_to_hash(&mut seen_terms, gene_uniquename,
                                              double_mutant_phenotype_termid);
                    }
                    if let Some(ref rescued_phenotype_termid) = interaction.rescued_phenotype_termid {
                        self.add_term_to_hash(&mut seen_terms, gene_uniquename,
                                              rescued_phenotype_termid);
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
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
                                                  gene_uniquename,
                                                  annotation_genotype_uniquename)
                    }
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &target_of_annotation.reference_uniquename);
                }

                for publication in &gene_details.feature_publications {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &Some(publication.clone()));
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

    fn set_genotype_details_maps(&mut self) {
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

        let mut maybe_add_to_gene_count_hash =
            |reference_uniquename: &FlexStr, gene_uniquename: &GeneUniquename| {
                if let Some(load_org_taxonid) = self.config.load_organism_taxonid {
                    if let Some(gene_details) = self.genes.get(gene_uniquename) {
                        if gene_details.taxonid == load_org_taxonid {
                            self.add_gene_to_hash(&mut gene_count_hash,
                                                  reference_uniquename,
                                                  gene_uniquename);
                        }
                    }
                }
            };


        let (_, mut seen_genes, mut seen_genotypes,
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
                                                             gene_uniquename);
                            }
                            for condition_termid in &annotation_detail.conditions {
                                self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                                      condition_termid);
                            }
                           if let Some(ref gene_product_form_id) =
                               annotation_detail.gene_product_form_id {
                                   if gene_product_form_id.starts_with("PR:") {
                                       self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                                             gene_product_form_id);
                                   }
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
                                                                     gene_uniquename);
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
                                    WithFromValue::Gene(ref gene_short) => {
                                        self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                                              &gene_short.uniquename);
                                        maybe_add_to_gene_count_hash(reference_uniquename,
                                                                     &gene_short.uniquename);
                                    },
                                    WithFromValue::Transcript(ref transcript_uniquename) => {
                                        self.add_transcript_to_hashes(&mut seen_transcripts, &mut seen_genes,
                                                                      reference_uniquename,
                                                                      transcript_uniquename);
                                    },
                                    _ => (),
                                }
                            }

                            if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                                let genotype = self.make_genotype_short(genotype_uniquename);
                                self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
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
                                                 &interaction.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction.interactor_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &interaction.interactor_uniquename);
                }

                let interaction_iter =
                    reference_details.genetic_interactions.iter();
                for interaction in interaction_iter {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction.gene_a_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &interaction.gene_a_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction.gene_b_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &interaction.gene_b_uniquename);

                    if let Some(ref genotype_a_uniquename) = interaction.genotype_a_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
                                                  reference_uniquename,
                                                  genotype_a_uniquename);
                    }
                    if let Some(ref genotype_b_uniquename) = interaction.genotype_b_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
                                                  reference_uniquename,
                                                  genotype_b_uniquename);
                    }

                    if let Some(ref double_mutant_phenotype_termid) = interaction.double_mutant_phenotype_termid {
                        self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                              double_mutant_phenotype_termid);
                    }
                    if let Some(ref rescued_phenotype_termid) = interaction.rescued_phenotype_termid {
                        self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                              rescued_phenotype_termid);
                    }
                }

                for ortholog_annotation in &reference_details.ortholog_annotations {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &ortholog_annotation.gene_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &ortholog_annotation.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &ortholog_annotation.ortholog_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &ortholog_annotation.ortholog_uniquename);

                }
                for paralog_annotation in &reference_details.paralog_annotations {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &paralog_annotation.gene_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &paralog_annotation.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &paralog_annotation.paralog_uniquename);
                    maybe_add_to_gene_count_hash(reference_uniquename,
                                                 &paralog_annotation.paralog_uniquename);
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
            if let Some(transcripts) = seen_transcripts.remove(reference_uniquename) {
                reference_details.transcripts_by_uniquename = transcripts;
            }
            if let Some(gene_count_genes) = gene_count_hash.remove(reference_uniquename) {
                reference_details.gene_count = gene_count_genes.len();
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
                            seen_genes.insert(gene_uniquename.clone());
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

            for (_, term_annotations) in genotype.cv_annotations() {
                for term_annotation in term_annotations {
                    annotation_count += term_annotation.annotations.len()
                }
            }

            genotype.annotation_count = annotation_count;
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
            if let Some(load_organism_taxonid) = self.config.load_organism_taxonid {
                if load_organism_taxonid != gene_details.taxonid {
                    continue;
                }
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
        for gene_details in self.genes.values() {
            if let Some(load_organism_taxonid) = self.config.load_organism_taxonid {
                if load_organism_taxonid != gene_details.taxonid {
                    continue;
                }
            }

            let subset_name =
                flex_str!("feature_type:") + &gene_details.feature_type;
            let re = Regex::new(r"[\s,:]+").unwrap();
            let subset_name_no_spaces = re.replace_all(&subset_name, "_").to_shared_str();
            subsets.entry(subset_name_no_spaces.clone())
                .or_insert(GeneSubsetDetails {
                    name: subset_name_no_spaces,
                    display_name: subset_name.clone(),
                    elements: HashSet::new()
                })
                .elements.insert(gene_details.uniquename.clone());
        }
    }

    // make subsets using the characterisation_status field of GeneDetails
    fn make_characterisation_status_subsets(&self, subsets: &mut IdGeneSubsetMap) {
        for gene_details in self.genes.values() {
            if let Some(load_organism_taxonid) = self.config.load_organism_taxonid {
                if load_organism_taxonid != gene_details.taxonid {
                    continue;
                }
            }

            if gene_details.feature_type != "mRNA gene" {
                continue;
            }

            if let Some(ref characterisation_status) = gene_details.characterisation_status {
                let subset_name =
                    flex_str!("characterisation_status:") + characterisation_status;
                let re = Regex::new(r"[\s,:]+").unwrap();
                let subset_name_no_spaces = re.replace_all(&subset_name, "_").to_shared_str();
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

                if !interpro_match.interpro_id.is_empty() {
                    let subset_name =
                        String::from("interpro:") + &interpro_match.interpro_id;
                    new_subset_names.push((subset_name.to_shared_str(),
                                           interpro_match.interpro_name.clone()));
                }

                let subset_name = String::from("interpro:") +
                     &interpro_match.dbname.clone() + ":" + &interpro_match.id;
                new_subset_names.push((subset_name.to_shared_str(), interpro_match.name.clone()));

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

    // populated the subsets HashMap
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

                if let Some(ref gene1_loc) = gene1.location {
                    if let Some(ref gene2_loc) = gene2.location {
                        let cmp = gene1_loc.start_pos.cmp(&gene2_loc.start_pos);
                        if cmp != Ordering::Equal {
                            return cmp;
                        }
                    }
                }
                if gene1.name.is_some() {
                    if gene2.name.is_some() {
                        gene1.name.cmp(&gene2.name)
                    } else {
                        Ordering::Less
                    }
                } else {
                    if gene2.name.is_some() {
                        Ordering::Greater
                    } else {
                        gene1.uniquename.cmp(&gene2.uniquename)
                    }
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

            if gene_details.feature_type == "mRNA gene" {
                if let Some(ref loc) = gene_details.location {
                    *coding_counts
                        .entry(&loc.chromosome_name)
                        .or_insert(0) += 1;
                }
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

    // remove some of the refs that have no annotations.
    // See: https://github.com/pombase/website/issues/628
    fn remove_non_curatable_refs(&mut self) {
        let filtered_refs = self.references.drain()
            .filter(|&(_, ref reference_details)| {
                if reference_has_annotation(reference_details) {
                    return true;
                }
                if let Some(ref triage_status) = reference_details.canto_triage_status {
                    if triage_status == "New" || triage_status == "Wrong organism" && triage_status == "Loaded in error"{
                        return false;
                    }
                }
                // default to true because there are references that
                // haven't or shouldn't be triaged, eg. GO_REF:...
                true
            })
            .collect();

        self.references = filtered_refs;
    }

    fn store_genetic_interactions(&mut self) {
        for interaction_annotation in self.genetic_interaction_annotations.drain(0..) {
            let gene_a_uniquename = &interaction_annotation.gene_a_uniquename;
            let gene_a = self.genes.get_mut(gene_a_uniquename).unwrap();
            gene_a.genetic_interactions.push(interaction_annotation.clone());

            let gene_b_uniquename = &interaction_annotation.gene_b_uniquename;
            let gene_b = self.genes.get_mut(gene_b_uniquename).unwrap();
            gene_b.genetic_interactions.push(interaction_annotation.clone());


            if let Some(ref reference_uniquename) =
                interaction_annotation.reference_uniquename
            {
                let reference = self.references.get_mut(reference_uniquename).unwrap();
                reference.genetic_interactions.push(interaction_annotation);

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

            let mut close_synonyms = vec![];
            let mut close_synonym_words_vec: Vec<FlexStr> = vec![];
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
                if synonym.synonym_type == "exact" || synonym.synonym_type == "narrow" {
                    add_to_words_vec(&synonym.name, &mut close_synonym_words_vec);
                    close_synonyms.push(synonym.name.clone());
                } else {
                    add_to_words_vec(&synonym.name, &mut distant_synonym_words_vec);
                    distant_synonyms.push(synonym.name.clone());
                }
            }

            distant_synonyms = distant_synonyms.into_iter()
                .filter(|synonym| {
                    !close_synonyms.contains(synonym)
                })
                .collect::<Vec<_>>();

            let annotation_count = term_details.annotation_count();
            let interesting_parent_ids_for_solr =
                term_details.interesting_parent_ids.clone();
            let term_summ = SolrTermSummary {
                id: termid.clone(),
                cv_name: term_details.cv_name.clone(),
                name: term_details.name.clone(),
                definition: term_details.definition.clone(),
                close_synonyms,
                close_synonym_words: join(&close_synonym_words_vec," "),
                distant_synonyms,
                distant_synonym_words: join(&distant_synonym_words_vec, " "),
                interesting_parent_ids: interesting_parent_ids_for_solr,
                secondary_identifiers: term_details.secondary_identifiers.clone(),
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
        let mut by_taxon = HashMap::new();

        for gene_details in self.genes.values() {
            let taxonid = gene_details.taxonid;

            by_taxon
                .entry(taxonid)
                .or_insert_with(StatCountsByTaxon::empty)
                .genes += 1;

            let mut annotation_count = 0;

            for term_annotations in gene_details.cv_annotations.values() {
                for term_annotation in term_annotations {
                    annotation_count += term_annotation.annotations.len();
                }
            }

            by_taxon
                .entry(taxonid)
                .or_insert_with(StatCountsByTaxon::empty)
                .annotations += annotation_count;
        }

        Stats {
            by_taxon,
            community_pubs_count: self.all_community_curated.len(),
            non_community_pubs_count: self.all_admin_curated.len(),
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
        self.process_annotation_feature_rels();
        self.add_target_of_annotations();
        self.set_deletion_viability();
        self.set_term_details_subsets();
        self.set_taxonomic_distributions();
        self.remove_non_curatable_refs();
        self.store_genetic_interactions();        self.set_term_details_maps();
        self.set_gene_details_maps();
        self.set_gene_details_subset_termids();
        self.set_genotype_details_maps();
        self.set_reference_details_maps();
        self.set_chromosome_gene_counts();
        self.set_counts();
        self.add_genotypes_to_allele_details();
        self.make_subsets();
        self.sort_chromosome_genes();
        self.set_gene_expression_measurements();

        let stats = self.get_stats();

        let metadata = self.make_metadata();

        let (gene_summaries, solr_gene_summaries) = self.make_gene_summaries();

        let solr_allele_summaries = self.make_solr_allele_summaries();
        let solr_term_summaries = self.make_solr_term_summaries();
        let solr_reference_summaries = self.make_solr_reference_summaries();

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

        WebData {
            metadata,
            chromosomes,
            chromosome_summaries,
            recent_references,
            all_community_curated,
            all_admin_curated,
            api_maps: self.make_api_maps(),
            search_gene_summaries: gene_summaries,
            solr_data,
            ont_annotations,
            stats,
        }
    }
}
