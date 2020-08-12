use std::rc::Rc;
use std::collections::BTreeMap;
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::usize;

use regex::Regex;

use std::collections::{HashMap, HashSet};

use crate::db::*;
use crate::types::*;
use crate::data_types::*;
use crate::web::data::*;
use crate::web::config::*;
use crate::web::cv_summary::make_cv_summaries;
use crate::web::cmp_utils::cmp_residues;
use crate::web::util::cmp_str_dates;

use crate::bio::util::rev_comp;

use pombase_rc_string::RcString;

use crate::interpro::UniprotResult;

fn make_organism(rc_organism: &Rc<Organism>) -> ConfigOrganism {
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

type UniprotIdentifier = RcString;

pub struct WebDataBuild<'a> {
    raw: &'a Raw,
    domain_data: &'a HashMap<UniprotIdentifier, UniprotResult>,
    config: &'a Config,

    genes: UniquenameGeneMap,
    genotypes: UniquenameGenotypeMap,
    genotype_backgrounds: HashMap<GenotypeUniquename, RcString>,
    alleles: UniquenameAlleleMap,
    other_features: UniquenameFeatureShortMap,
    terms: TermIdDetailsMap,
    chromosomes: ChrNameDetailsMap,
    references: UniquenameReferenceMap,
    all_ont_annotations: HashMap<TermId, Vec<OntAnnotationId>>,
    all_not_ont_annotations: HashMap<TermId, Vec<OntAnnotationId>>,

    // map from term name to term ID (ie "nucleus" -> "GO:0005634")
    term_ids_by_name: HashMap<RcString, TermId>,

    genes_of_transcripts: HashMap<RcString, RcString>,
    transcripts_of_polypeptides: HashMap<RcString, RcString>,
    parts_of_transcripts: HashMap<RcString, Vec<FeatureShort>>,
    genes_of_alleles: HashMap<RcString, RcString>,
    alleles_of_genotypes: HashMap<RcString, Vec<ExpressedAllele>>,

    // a map from IDs of terms from the "PomBase annotation extension terms" cv
    // to a Vec of the details of each of the extension
    parts_of_extensions: HashMap<TermId, Vec<ExtPart>>,

    base_term_of_extensions: HashMap<TermId, TermId>,

    // a set of child terms for each term from the cvtermpath table
    children_by_termid: HashMap<TermId, HashSet<TermId>>,
    dbxrefs_of_features: HashMap<RcString, HashSet<RcString>>,

    possible_interesting_parents: HashSet<InterestingParent>,

    recent_references: RecentReferences,
    all_community_curated: Vec<ReferenceShort>,
    all_admin_curated: Vec<ReferenceShort>,

    term_subsets: IdTermSubsetMap,
    gene_subsets: IdGeneSubsetMap,

    annotation_details: IdOntAnnotationDetailMap,

    ont_annotations: Vec<OntAnnotation>,
}

fn get_maps() ->
    (HashMap<RcString, ReferenceShortOptionMap>,
     HashMap<RcString, GeneShortOptionMap>,
     HashMap<RcString, GenotypeShortMap>,
     HashMap<RcString, AlleleShortMap>,
     HashMap<GeneUniquename, TermShortOptionMap>)
{
    (HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new())
}

fn get_feat_rel_expression(feature: &Feature,
                           feature_relationship: &FeatureRelationship) -> Option<RcString> {
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

fn reference_has_annotation(reference_details: &ReferenceDetails) -> bool {
    !reference_details.cv_annotations.is_empty() ||
        !reference_details.physical_interactions.is_empty() ||
        !reference_details.genetic_interactions.is_empty() ||
        !reference_details.ortholog_annotations.is_empty() ||
        !reference_details.paralog_annotations.is_empty()
}

fn is_gene_type(feature_type_name: &str) -> bool {
    feature_type_name == "gene" || feature_type_name == "pseudogene"
}

pub fn compare_ext_part_with_config(config: &Config, ep1: &ExtPart, ep2: &ExtPart) -> Ordering {
    let rel_order_conf = &config.extension_relation_order;
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

fn string_from_ext_range(ext_range: &ExtRange,
                         genes: &UniquenameGeneMap, terms: &TermIdDetailsMap) -> RcString {
    match *ext_range {
        ExtRange::Gene(ref gene_uniquename) | ExtRange::Promoter(ref gene_uniquename) => {
            let gene = genes.get(gene_uniquename)
                .unwrap_or_else(|| panic!("can't find gene: {}", gene_uniquename));
            gene_display_name(gene)
        },
        ExtRange::SummaryGenes(_) => panic!("can't handle SummaryGenes\n"),
        ExtRange::Term(ref termid) => RcString::from(&terms.get(termid).unwrap().name),
        ExtRange::SummaryModifiedResidues(ref residue) => RcString::from(&residue.join(",")),
        ExtRange::SummaryTerms(_) => panic!("can't handle SummaryGenes\n"),
        ExtRange::Misc(ref misc) => misc.clone(),
        ExtRange::Domain(ref domain) => domain.clone(),
        ExtRange::GeneProduct(ref gene_product) => gene_product.clone(),
    }
}

fn cmp_ext_part(ext_part1: &ExtPart, ext_part2: &ExtPart,
                genes: &UniquenameGeneMap,
                terms: &TermIdDetailsMap) -> Ordering {
    let ord = ext_part1.rel_type_display_name.cmp(&ext_part2.rel_type_display_name);

    if ord == Ordering::Equal {
        let ext_part1_str = string_from_ext_range(&ext_part1.ext_range, genes, terms);
        let ext_part2_str = string_from_ext_range(&ext_part2.ext_range, genes, terms);

        ext_part1_str.to_lowercase().cmp(&ext_part2_str.to_lowercase())
    } else {
        ord
    }
}

// compare the extension up to the last common index
fn cmp_extension_prefix(cv_config: &CvConfig, ext1: &[ExtPart], ext2: &[ExtPart],
                        genes: &UniquenameGeneMap,
                        terms: &TermIdDetailsMap) -> Ordering {
    let conf_rel_ranges = &cv_config.summary_relation_ranges_to_collect;

    let is_grouping_rel_name =
        |ext: &ExtPart| !conf_rel_ranges.contains(&ext.rel_type_name);

    // put the extension that will be grouped in the summary at the end
    // See: https://github.com/pombase/pombase-chado/issues/636
    let (mut ext1_for_cmp, ext1_rest): (Vec<ExtPart>, Vec<ExtPart>) =
        ext1.to_vec().into_iter().partition(&is_grouping_rel_name);
    ext1_for_cmp.extend(ext1_rest.into_iter());

    let (mut ext2_for_cmp, ext2_rest): (Vec<ExtPart>, Vec<ExtPart>) =
        ext2.to_vec().into_iter().partition(&is_grouping_rel_name);
    ext2_for_cmp.extend(ext2_rest.into_iter());

    let iter = ext1_for_cmp.iter().zip(&ext2_for_cmp).enumerate();
    for (_, (ext1_part, ext2_part)) in iter {
        let ord = cmp_ext_part(ext1_part, ext2_part, genes, terms);

        if ord != Ordering::Equal {
            return ord
        }
    }

    Ordering::Equal
}

fn cmp_extension(cv_config: &CvConfig, ext1: &[ExtPart], ext2: &[ExtPart],
                 genes: &UniquenameGeneMap,
                 terms: &TermIdDetailsMap) -> Ordering {
    let cmp = cmp_extension_prefix(cv_config, ext1, ext2, genes, terms);

    if cmp == Ordering::Equal {
        ext1.len().cmp(&ext2.len())
    } else {
        cmp
    }
}

fn cmp_genotypes(genotype1: &GenotypeDetails, genotype2: &GenotypeDetails) -> Ordering {
    genotype1.display_uniquename.to_lowercase().cmp(&genotype2.display_uniquename.to_lowercase())
}


fn allele_display_name(allele: &AlleleShort) -> RcString {
    let name = allele.name.clone().unwrap_or_else(|| RcString::from("unnamed"));
    let allele_type = allele.allele_type.clone();
    let description = allele.description.clone().unwrap_or_else(|| allele_type.clone());

    if allele_type == "deletion" && name.ends_with("delta") ||
        allele_type.starts_with("wild_type") && name.ends_with('+') {
            let normalised_description = description.replace("[\\s_]+", "");
            let normalised_allele_type = allele_type.replace("[\\s_]+", "");
            if normalised_description != normalised_allele_type {
                return RcString::from(&(name + "(" + description.as_str() + ")"));
            } else {
                return name;
            }
        }

    let display_name =
        if allele_type == "deletion" {
            name + "-" + description.as_str()
        } else {
            name + "-" + description.as_str() + "-" + &allele.allele_type
        };
    RcString::from(&display_name)
}


fn gene_display_name(gene: &GeneDetails) -> RcString {
    if let Some(ref name) = gene.name {
        name.clone()
    } else {
        gene.uniquename.clone()
    }
}

lazy_static! {
    static ref BAD_GENOTYPE_NAME_CHARS_RE: Regex =
        Regex::new(r"[% /&;]").unwrap();
}

pub fn make_genotype_display_name(genotype_expressed_alleles: &[ExpressedAllele],
                                  allele_map: &UniquenameAlleleMap) -> RcString {
    let mut allele_display_names: Vec<String> =
        genotype_expressed_alleles.iter().map(|expressed_allele| {
            let allele_short = allele_map.get(&expressed_allele.allele_uniquename).unwrap();
            let mut display_name = allele_display_name(allele_short).to_string();
            if allele_short.allele_type != "deletion" {
                if display_name == "unnamed-unrecorded-unrecorded" {
                    display_name = format!("{}-{}", allele_short.gene_uniquename,
                                            display_name);
                }
                if let Some(ref expression) = expressed_allele.expression {
                    display_name += &format!("-expression-{}", expression.to_lowercase());
                }
            }
            display_name
        }).collect();

    allele_display_names.sort();

    let joined_alleles = allele_display_names.join("  ");

    let clean_display_name =
        BAD_GENOTYPE_NAME_CHARS_RE.replace_all(&joined_alleles, "_");
    RcString::from(&clean_display_name)
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
                phase: make_phase(&feature_loc),
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
                                 chromosome_name: &'a str) -> ChromosomeShort {
    if let Some(chr) = chromosome_map.get(chromosome_name) {
        chr.make_chromosome_short()
    } else {
        panic!("can't find chromosome: {}", chromosome_name);
    }
}

fn make_gene_short<'b>(gene_map: &'b UniquenameGeneMap,
                       gene_uniquename: &'b str) -> GeneShort {
    if let Some(gene_details) = gene_map.get(gene_uniquename) {
        GeneShort {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
        }
    } else {
        panic!("can't find GeneDetails for gene uniquename {}", gene_uniquename)
    }
}

fn make_reference_short<'a>(reference_map: &'a UniquenameReferenceMap,
                            reference_uniquename: &str) -> Option<ReferenceShort> {
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

// compare two gene vectors which must be ordered vecs
fn cmp_gene_vec(genes: &UniquenameGeneMap,
                gene_vec1: &[GeneUniquename],
                gene_vec2: &[GeneUniquename]) -> Ordering {

    let gene_short_vec1: Vec<GeneShort> =
        gene_vec1.iter().map(|gene_uniquename: &RcString| {
            make_gene_short(genes, gene_uniquename)
        }).collect();
    let gene_short_vec2: Vec<GeneShort> =
        gene_vec2.iter().map(|gene_uniquename: &RcString| {
            make_gene_short(genes, gene_uniquename)
        }).collect();

    gene_short_vec1.cmp(&gene_short_vec2)
}

lazy_static! {
    static ref PROMOTER_RE: Regex = Regex::new(r"^(?P<gene>.*)-promoter$").unwrap();
}

pub fn cmp_ont_annotation_detail(cv_config: &CvConfig,
                             detail1: &OntAnnotationDetail,
                             detail2: &OntAnnotationDetail,
                             genes: &UniquenameGeneMap,
                             genotypes: &UniquenameGenotypeMap,
                             terms: &TermIdDetailsMap) -> Result<Ordering, String> {
    if let Some(ref detail1_genotype_uniquename) = detail1.genotype {
        if let Some(ref detail2_genotype_uniquename) = detail2.genotype {
            let genotype1 = genotypes.get(detail1_genotype_uniquename).unwrap();
            let genotype2 = genotypes.get(detail2_genotype_uniquename).unwrap();

            let ord = cmp_genotypes(genotype1, genotype2);

            if ord == Ordering::Equal {
                Ok(cmp_extension(cv_config, &detail1.extension, &detail2.extension,
                                 genes, terms))
            } else {
                Ok(ord)
            }
        } else {
            Err(format!("comparing two OntAnnotationDetail but one has a genotype and
 one a gene:\n{:?}\n{:?}\n", detail1, detail2))
        }
    } else {
        if detail2.genotype.is_some() {
            Err(format!("comparing two OntAnnotationDetail but one has a genotype and
 one a gene:\n{:?}\n{:?}\n", detail1, detail2))
        } else {
            let ord = cmp_gene_vec(genes, &detail1.genes, &detail2.genes);

            if ord == Ordering::Equal {
                if let Some(ref sort_details_by) = cv_config.sort_details_by {
                    for sort_type in sort_details_by {
                        if sort_type == "modification" {
                            let res = cmp_residues(&detail1.residue, &detail2.residue);
                            if res != Ordering::Equal {
                                return Ok(res);
                            }
                        } else {
                            let res = cmp_extension(cv_config, &detail1.extension,
                                                    &detail2.extension, genes, terms);
                            if res != Ordering::Equal {
                                return Ok(res);
                            }
                        }
                    }
                    Ok(Ordering::Equal)
                } else {
                    Ok(cmp_extension(cv_config, &detail1.extension, &detail2.extension,
                                     genes, terms))
                }
            } else {
                Ok(ord)
            }
        }
    }
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
                rel_name: RcString::from("is_a"),
            });
        }
    }

    let add_to_set = |set: &mut HashSet<_>, termid: RcString| {
        for rel_name in &DESCENDANT_REL_NAMES {
            set.insert(InterestingParent {
                termid: termid.to_owned(),
                rel_name: RcString::from(rel_name),
            });
        }
    };

    for slim_config in config.slims.values() {
        for go_slim_conf in &slim_config.terms {
            add_to_set(&mut ret, go_slim_conf.termid.clone());
        }
    }

    for field_name in &config.gene_results.visualisation_field_names {
        let column_conf = &config.gene_results.field_config[field_name];

        for attr_value_conf in &column_conf.attr_values {
            if let Some(ref termid) = attr_value_conf.termid {
                add_to_set(&mut ret, termid.clone());
            }
        }
    }

    ret.insert(InterestingParent {
        termid: config.viability_terms.viable.clone(),
        rel_name: RcString::from("is_a"),
    });
    ret.insert(InterestingParent {
        termid: config.viability_terms.inviable.clone(),
        rel_name: RcString::from("is_a"),
    });

    let add_filter_ancestor =
        |set: &mut HashSet<_>, category: &AncestorFilterCategory, cv_name: &str| {
            for ancestor in &category.ancestors {
                for config_rel_name in &DESCENDANT_REL_NAMES {
                    if *config_rel_name == "has_part" &&
                        !HAS_PART_CV_NAMES.contains(&cv_name) {
                            continue;
                        }
                    set.insert(InterestingParent {
                        termid: ancestor.clone(),
                        rel_name: RcString::from(*config_rel_name),
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
                    if ancestor.starts_with("NOT ") {
                        RcString::from(&ancestor[4..])
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
                       all_ref_uniquenames: &[RcString]) -> Vec<ReferenceShort> {
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
                      all_ref_uniquenames: &[RcString])
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
                             transcript_uniquename: &str, parts: &mut Vec<FeatureShort>) {
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

                intron_count += 1;

                let new_intron_loc = ChromosomeLocation {
                    chromosome_name: prev_part.location.chromosome_name.clone(),
                    start_pos: intron_start,
                    end_pos: intron_end,
                    strand: prev_part.location.strand.clone(),
                    phase: None,
                };

                let intron_uniquename =
                    format!("{}:intron:{}", transcript_uniquename, intron_count);
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
                    uniquename: RcString::from(&intron_uniquename),
                    location: new_intron_loc,
                    residues: intron_residues,
                });
            }
        }

        if let Some(new_intron) = maybe_new_intron {
            new_parts.push(new_intron);
        }

        new_parts.push(part);
    }

    *parts = new_parts;
}

fn validate_transcript_parts(transcript_uniquename: &str, parts: &[FeatureShort]) {
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
    pub fn new(raw: &'a Raw, domain_data: &'a HashMap<UniprotIdentifier, UniprotResult>,
               config: &'a Config) -> WebDataBuild<'a>
    {
        WebDataBuild {
            raw,
            domain_data,
            config,

            genes: BTreeMap::new(),
            genotypes: HashMap::new(),
            genotype_backgrounds: HashMap::new(),
            alleles: HashMap::new(),
            other_features: HashMap::new(),
            terms: HashMap::new(),
            chromosomes: BTreeMap::new(),
            references: HashMap::new(),
            all_ont_annotations: HashMap::new(),
            all_not_ont_annotations: HashMap::new(),
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
            alleles_of_genotypes: HashMap::new(),

            parts_of_extensions: HashMap::new(),

            base_term_of_extensions: HashMap::new(),

            children_by_termid: HashMap::new(),
            dbxrefs_of_features: HashMap::new(),

            possible_interesting_parents: get_possible_interesting_parents(config),

            term_subsets: HashMap::new(),
            gene_subsets: HashMap::new(),

            annotation_details: HashMap::new(),

            ont_annotations: vec![],
       }
    }

    fn add_ref_to_hash(&self,
                       seen_references: &mut HashMap<RcString, ReferenceShortOptionMap>,
                       identifier: &str,
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
                        seen_genes: &mut HashMap<RcString, GeneShortOptionMap>,
                        identifier: &RcString,
                        other_gene_uniquename: &GeneUniquename) {
        seen_genes
            .entry(identifier.clone())
            .or_insert_with(HashMap::new)
            .insert(other_gene_uniquename.clone(), None);
    }

    fn add_genotype_to_hash(&self,
                            seen_genotypes: &mut HashMap<RcString, GenotypeShortMap>,
                            seen_alleles: &mut HashMap<RcString, AlleleShortMap>,
                            seen_genes: &mut HashMap<RcString, GeneShortOptionMap>,
                            identifier: &RcString,
                            genotype_uniquename: &RcString) {
        let genotype_short = self.make_genotype_short(&genotype_uniquename);
        for expressed_allele in &genotype_short.expressed_alleles {
            self.add_allele_to_hash(seen_alleles, seen_genes, identifier,
                                    &expressed_allele.allele_uniquename);
        }

        seen_genotypes
            .entry(identifier.clone())
            .or_insert_with(HashMap::new)
            .insert(genotype_uniquename.clone(), genotype_short);
    }

    fn add_allele_to_hash(&self,
                          seen_alleles: &mut HashMap<RcString, AlleleShortMap>,
                          seen_genes: &mut HashMap<RcString, GeneShortOptionMap>,
                          identifier: &RcString,
                          allele_uniquename: &AlleleUniquename) -> AlleleShort {
        let allele_short = self.make_allele_short(&allele_uniquename);
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

    fn add_term_to_hash(&self,
                        seen_terms: &mut HashMap<TermId, TermShortOptionMap>,
                        identifier: &RcString,
                        other_termid: &TermId) {
        seen_terms
            .entry(identifier.clone())
            .or_insert_with(HashMap::new)
            .insert(other_termid.clone(), None);
    }

    fn get_gene<'b>(&'b self, gene_uniquename: &'b str) -> &'b GeneDetails {
        if let Some(gene_details) = self.genes.get(gene_uniquename) {
            gene_details
        } else {
            panic!("can't find GeneDetails for gene uniquename {}", gene_uniquename)
        }
    }

    fn get_gene_mut<'b>(&'b mut self, gene_uniquename: &'b str) -> &'b mut GeneDetails {
        if let Some(gene_details) = self.genes.get_mut(gene_uniquename) {
            gene_details
        } else {
            panic!("can't find GeneDetails for gene uniquename {}", gene_uniquename)
        }
    }

    fn make_gene_short(&self, gene_uniquename: &str) -> GeneShort {
        let gene_details = self.get_gene(gene_uniquename);
        GeneShort {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
        }
    }

    fn make_gene_summary(&self, gene_uniquename: &str) -> GeneSummary {
        let gene_details = self.get_gene(gene_uniquename);
        let synonyms =
            gene_details.synonyms.iter()
            .filter(|synonym| synonym.synonym_type == "exact")
            .map(|synonym| synonym.name.clone())
            .collect::<Vec<RcString>>();
        let ortholog_ids =
            gene_details.ortholog_annotations.iter()
            .map(|ortholog_annotation| {
                let ortholog_uniquename = ortholog_annotation.ortholog_uniquename.clone();
                let maybe_ortholog_name =
                    self.genes.get(&ortholog_uniquename).unwrap().name.clone();

                IdNameAndOrganism {
                    identifier: ortholog_uniquename,
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
            synonyms,
            orthologs: ortholog_ids,
            feature_type: gene_details.feature_type.clone(),
            taxonid: gene_details.taxonid,
            location: gene_details.location.clone(),
        }
    }

    fn make_api_gene_summary(&self, gene_uniquename: &str) -> APIGeneSummary {
        let gene_details = self.get_gene(gene_uniquename);
        let synonyms =
            gene_details.synonyms.iter()
            .filter(|synonym| synonym.synonym_type == "exact")
            .map(|synonym| synonym.name.clone())
            .collect::<Vec<RcString>>();
        let exon_count =
            if let Some(transcript) = gene_details.transcripts.get(0) {
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
        APIGeneSummary {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
            uniprot_identifier: gene_details.uniprot_identifier.clone(),
            exact_synonyms: synonyms,
            dbxrefs: gene_details.dbxrefs.clone(),
            location: gene_details.location.clone(),
            transcripts: gene_details.transcripts.clone(),
            tm_domain_count: gene_details.tm_domain_coords.len(),
            exon_count,
        }
    }

    fn make_term_short(&self, termid: &str) -> TermShort {
        if let Some(term_details) = self.terms.get(termid) {
            TermShort::from_term_details(&term_details)
        } else {
            panic!("can't find TermDetails for termid: {}", termid)
        }
    }

    fn add_characterisation_status(&mut self, gene_uniquename: &str,
                                   cvterm_name: &RcString) {
        let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        gene_details.characterisation_status = Some(cvterm_name.clone());
    }

    fn add_gene_product(&mut self, gene_uniquename: &str, product: &RcString) {
        let gene_details = self.get_gene_mut(gene_uniquename);
        gene_details.product = Some(product.clone());
    }

    fn add_name_description(&mut self, gene_uniquename: &str, name_description: &str) {
        let gene_details = self.get_gene_mut(gene_uniquename);
        gene_details.name_descriptions.push(name_description.into());
    }

    fn add_annotation(&mut self, cvterm: &Cvterm, is_not: bool,
                      annotation_template: OntAnnotationDetail) {
        let termid =
            match self.base_term_of_extensions.get(&cvterm.termid()) {
                Some(base_termid) => base_termid.clone(),
                None => cvterm.termid(),
            };

        let extension_parts =
            match self.parts_of_extensions.get(&cvterm.termid()) {
                Some(parts) => parts.clone(),
                None => vec![],
            };

        let mut new_extension = extension_parts.clone();

        let mut existing_extensions = annotation_template.extension.clone();
        new_extension.append(&mut existing_extensions);

        {
            let compare_ext_part_func =
                |e1: &ExtPart, e2: &ExtPart| compare_ext_part_with_config(self.config, e1, e2);

            new_extension.sort_by(compare_ext_part_func);
        };

        let ont_annotation_detail =
            OntAnnotationDetail {
                extension: new_extension,
                .. annotation_template
            };

        let annotation_map = if is_not {
            &mut self.all_not_ont_annotations
        } else {
            &mut self.all_ont_annotations
        };

        let entry = annotation_map.entry(termid.clone());
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
        let mut all_uniquenames = vec![];

        for rc_publication in &self.raw.publications {
            let reference_uniquename = &rc_publication.uniquename;

            if reference_uniquename.to_lowercase() == "null" {
                continue;
            }

            let mut pubmed_authors: Option<RcString> = None;
            let mut pubmed_publication_date: Option<RcString> = None;
            let mut pubmed_abstract: Option<RcString> = None;
            let mut pubmed_doi: Option<RcString> = None;
            let mut canto_annotation_status: Option<RcString> = None;
            let mut canto_triage_status: Option<RcString> = None;
            let mut canto_curator_role: Option<RcString> = None;
            let mut canto_curator_name: Option<RcString> = None;
            let mut canto_first_approved_date: Option<RcString> = None;
            let mut canto_approved_date: Option<RcString> = None;
            let mut canto_added_date: Option<RcString> = None;
            let mut canto_session_submitted_date: Option<RcString> = None;

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
                    "canto_annotation_status" =>
                        canto_annotation_status = Some(prop.value.clone()),
                    "canto_triage_status" =>
                        canto_triage_status = Some(prop.value.clone()),
                    "canto_curator_role" =>
                        canto_curator_role = Some(prop.value.clone()),
                    "canto_curator_name" =>
                        canto_curator_name = Some(prop.value.clone()),
                    "canto_first_approved_date" =>
                        canto_first_approved_date = Some(prop.value.clone()),
                    "canto_approved_date" =>
                        canto_approved_date = Some(prop.value.clone()),
                    "canto_added_date" =>
                        canto_added_date = Some(prop.value.clone()),
                    "canto_session_submitted_date" =>
                        canto_session_submitted_date = Some(prop.value.clone()),
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
                    authors_abbrev = Some(RcString::from(&replaced));
                } else {
                    authors_abbrev = Some(authors.clone());
                }
            }

            if let Some(publication_date) = pubmed_publication_date.clone() {
                let date_re = Regex::new(r"^(.* )?(?P<y>\d\d\d\d)$").unwrap();
                publication_year = Some(RcString::from(&date_re.replace_all(&publication_date, "$y")));
            }

            let mut approved_date = canto_first_approved_date.clone();

            if approved_date.is_none() {
                approved_date = canto_session_submitted_date.clone();
            }

            approved_date =
                if let Some(date) = approved_date {
                    let re = Regex::new(r"^(?P<date>\d\d\d\d-\d\d-\d\d).*").unwrap();
                    Some(RcString::from(&re.replace_all(&date, "$date")))
                } else {
                    None
                };

            if let Some(ref canto_annotation_status) = canto_annotation_status {
                if canto_annotation_status != "APPROVED" {
                    approved_date = None;
                }
            }

            self.references.insert(reference_uniquename.clone(),
                                   ReferenceDetails {
                                       uniquename: reference_uniquename.clone(),
                                       title: rc_publication.title.clone(),
                                       citation: rc_publication.miniref.clone(),
                                       pubmed_abstract,
                                       pubmed_doi,
                                       authors: pubmed_authors.clone(),
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
                                       terms_by_termid: HashMap::new(),
                                       annotation_details: HashMap::new(),
                                   });

            if pubmed_publication_date.is_some() {
                all_uniquenames.push(reference_uniquename.clone());
            }
        }

        let (recent_admin_curated, recent_community_curated,
             all_community_curated, all_admin_curated) =
            make_canto_curated(&self.references, &all_uniquenames);

        let recent_references = RecentReferences {
            pubmed: make_recently_added(&self.references, &all_uniquenames),
            admin_curated: recent_admin_curated,
            community_curated: recent_community_curated,
        };

        self.recent_references = recent_references;

        self.all_community_curated = all_community_curated;
        self.all_admin_curated = all_admin_curated;
    }

    // make maps from genes to transcript, transcripts to polypeptide,
    // exon, intron, UTRs
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
                        let allele_and_expression =
                            ExpressedAllele {
                                allele_uniquename: subject_uniquename.clone(),
                                expression,
                            };
                        let entry = self.alleles_of_genotypes.entry(object_uniquename.clone());
                        entry.or_insert_with(Vec::new).push(allele_and_expression);
                        continue;
                    }
            }
            if TRANSCRIPT_PART_TYPES.contains(&subject_type_name.as_str()) {
                let entry = self.parts_of_transcripts.entry(object_uniquename.clone());
                let part = make_feature_short(&self.chromosomes, &feature_rel.subject);
                entry.or_insert_with(Vec::new).push(part);
            }
        }
    }

    fn get_feature_dbxrefs(&self, feature: &Feature) -> HashSet<RcString> {
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
            if dbxref.starts_with("SPD:") {
                orfeome_identifier = Some(RcString::from(&dbxref[4..]));
            }
        }

        let mut uniprot_identifier = None;
        let mut biogrid_interactor_id: Option<u32> = None;

        for prop in feat.featureprops.borrow().iter() {
            if prop.prop_type.name == "uniprot_identifier" {
                uniprot_identifier = prop.value.clone();
            } else {
                if prop.prop_type.name == "biogrid_interactor_id" {
                    if let Some(ref chado_biogrid_id) = prop.value {
                        biogrid_interactor_id = match chado_biogrid_id.parse::<u32>() {
                            Ok(val) => Some(val),
                            Err(err) =>
                                panic!("error parsing BioGRID interactor ID from Chado: {}", err),
                        }
                    }
                }
            }
        }

        let (interpro_matches, tm_domain_coords) =
            if let Some(ref uniprot_identifier) = uniprot_identifier {
                if let Some(result) = self.domain_data.get(uniprot_identifier as &str) {
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

        let gene_feature = GeneDetails {
            uniquename: feat.uniquename.clone(),
            name: feat.name.clone(),
            taxonid: organism.taxonid,
            product: None,
            deletion_viability: DeletionViability::Unknown,
            uniprot_identifier,
            biogrid_interactor_id,
            interpro_matches,
            tm_domain_coords,
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

    fn get_transcript_parts(&mut self, transcript_uniquename: &str) -> Vec<FeatureShort> {
        if let Some(mut parts) = self.parts_of_transcripts.remove(transcript_uniquename) {
            if parts.is_empty() {
                panic!("transcript has no parts: {}", transcript_uniquename);
            }

            let part_cmp = |a: &FeatureShort, b: &FeatureShort| {
                a.location.start_pos.cmp(&b.location.start_pos)
            };

            parts.sort_by(&part_cmp);

            validate_transcript_parts(transcript_uniquename, &parts);

            let chr_name = &parts[0].location.chromosome_name.clone();
            if let Some(chromosome) = self.chromosomes.get(chr_name) {
                add_introns_to_transcript(chromosome, transcript_uniquename, &mut parts);
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
                            strand: first_part_loc.strand.clone(),
                            phase: make_phase(&mrna_location),
                        })
                    } else {
                        None
                    }
                }
            } else {
                None
            };

        let transcript = TranscriptDetails {
            uniquename: transcript_uniquename.clone(),
            location: transcript_location,
            transcript_type: feat.feat_type.name.clone(),
            parts,
            protein: None,
            cds_location: maybe_cds_location,
        };

        if let Some(gene_uniquename) =
            self.genes_of_transcripts.get(&transcript_uniquename) {
                let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
                if gene_details.feature_type != "pseudogene" {
                    let feature_type =
                        transcript.transcript_type.clone() + " " + &gene_details.feature_type;
                    gene_details.feature_type = RcString::from(&feature_type);
                }
                gene_details.transcripts.push(transcript);
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

            let parse_prop_as_f32 = |p: &Option<RcString>| {
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
                sequence: RcString::from(&residues),
                molecular_weight: molecular_weight.unwrap(),
                average_residue_weight: average_residue_weight.unwrap(),
                charge_at_ph7: charge_at_ph7.unwrap(),
                isoelectric_point: isoelectric_point.unwrap(),
                codon_adaptation_index: codon_adaptation_index.unwrap(),
            };

            if let Some(transcript_uniquename) =
                self.transcripts_of_polypeptides.get(&protein_uniquename) {
                    if let Some(gene_uniquename) =
                        self.genes_of_transcripts.get(transcript_uniquename) {
                            let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
                            if gene_details.transcripts.len() > 1 {
                                panic!("unimplemented - can't handle multiple transcripts for: {}",
                                       gene_uniquename);
                            } else {
                                if gene_details.transcripts.is_empty() {
                                    panic!("gene has no transcript: {}", gene_uniquename);
                                } else {
                                    gene_details.transcripts[0].protein = Some(protein);
                                }
                            }
                        } else {
                            panic!("can't find gene for transcript: {}", transcript_uniquename);
                        }
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
            residues: RcString::from(&residues),
            ena_identifier: RcString::from(&ena_identifier.unwrap()),
            gene_uniquenames: vec![],
            taxonid: org.taxonid,
        };

        self.chromosomes.insert(feat.uniquename.clone(), chr);
    }

    fn store_genotype_details(&mut self, feat: &Feature) {
        let mut expressed_alleles =
            self.alleles_of_genotypes[&feat.uniquename].clone();
        let genotype_display_uniquename =
            make_genotype_display_name(&expressed_alleles, &self.alleles);

        {
            let allele_cmp = |allele1: &ExpressedAllele, allele2: &ExpressedAllele| {
                let allele1_display_name =
                    allele_display_name(&self.alleles[&allele1.allele_uniquename]);
                let allele2_display_name =
                    allele_display_name(&self.alleles[&allele2.allele_uniquename]);
                allele1_display_name.cmp(&allele2_display_name)
            };

            expressed_alleles.sort_by(&allele_cmp);
        }

        for prop in feat.featureprops.borrow().iter() {
            if prop.prop_type.name == "genotype_background" {
                if let Some(ref background) = prop.value {
                    self.genotype_backgrounds.insert(feat.uniquename.clone(),
                                                     background.clone());
                }
            }
        }

        let rc_display_name = RcString::from(&genotype_display_uniquename);

        self.genotypes.insert(rc_display_name.clone(),
                              GenotypeDetails {
                                  display_uniquename: rc_display_name,
                                  name: feat.name.as_ref().map(|s| RcString::from(s)),
                                  expressed_alleles,
                                  cv_annotations: HashMap::new(),
                                  genes_by_uniquename: HashMap::new(),
                                  alleles_by_uniquename: HashMap::new(),
                                  references_by_uniquename: HashMap::new(),
                                  terms_by_termid: HashMap::new(),
                                  annotation_details: HashMap::new(),
                              });
    }

    fn store_allele_details(&mut self, feat: &Feature) {
        let mut allele_type = None;
        let mut description = None;

        for prop in feat.featureprops.borrow().iter() {
            match &prop.prop_type.name as &str {
                "allele_type" =>
                    allele_type = prop.value.clone(),
                "description" =>
                    description = prop.value.clone(),
                _ => ()
            }
        }

        if allele_type.is_none() {
            panic!("no allele_type cvtermprop for {}", &feat.uniquename);
        }

        let gene_uniquename =
            self.genes_of_alleles[&feat.uniquename].clone();
        let allele_details = AlleleShort {
            uniquename: feat.uniquename.clone(),
            name: feat.name.clone(),
            gene_uniquename,
            allele_type: allele_type.unwrap(),
            description,
        };
        self.alleles.insert(feat.uniquename.clone(), allele_details);
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
                    let feature_short = make_feature_short(&self.chromosomes, &feat);
                    self.other_features.insert(feat.uniquename.clone(), feature_short);
                }
            }
        }
    }

    fn add_interesting_parents(&mut self) {
        let mut interesting_parents_by_termid: HashMap<RcString, HashSet<RcString>> =
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
                    .insert(object_termid);
            };
        }

        for (termid, interesting_parents) in interesting_parents_by_termid {
            let term_details = self.terms.get_mut(&termid).unwrap();
            term_details.interesting_parents = interesting_parents;
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

    fn add_gene_neighbourhoods(&mut self) {
        struct GeneAndLoc {
            gene_uniquename: RcString,
            loc: ChromosomeLocation,
        };

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


    // add interaction, ortholog and paralog annotations
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
                        let mut interaction_note: Option<RcString> = None;

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
                                    let interaction_annotation =
                                        InteractionAnnotation {
                                            gene_uniquename: gene_uniquename.clone(),
                                            interactor_uniquename: other_gene_uniquename.clone(),
                                            evidence,
                                            reference_uniquename: maybe_reference_uniquename.clone(),
                                            throughput,
                                            interaction_note,
                                        };
                                    {
                                        let gene_details = self.genes.get_mut(subject_uniquename).unwrap();
                                        if rel_name == "interacts_physically" {
                                            gene_details.physical_interactions.push(interaction_annotation.clone());
                                        } else {
                                            if rel_name == "interacts_genetically" {
                                                gene_details.genetic_interactions.push(interaction_annotation.clone());
                                            } else {
                                                panic!("unknown interaction type: {}", rel_name);
                                            }
                                        };
                                    }
                                    if gene_uniquename != other_gene_uniquename {
                                        let other_gene_details = self.genes.get_mut(object_uniquename).unwrap();
                                        if rel_name == "interacts_physically" {
                                            other_gene_details.physical_interactions.push(interaction_annotation.clone());
                                        } else {
                                            if rel_name == "interacts_genetically" {
                                                other_gene_details.genetic_interactions.push(interaction_annotation.clone());
                                            } else {
                                                panic!("unknown interaction type: {}", rel_name);
                                            }
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
                                            if rel_name == "interacts_genetically" {
                                                ref_details.genetic_interactions.push(interaction_annotation.clone());
                                            } else {
                                                panic!("unknown interaction type: {}", rel_name);
                                            }
                                        };
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
                                    // avoid duplicates in the reference pages
                                    if self.config.load_organism_taxonid.is_some() &&
                                        self.config.load_organism_taxonid.unwrap() == gene_details.taxonid ||
                                        gene_organism_taxonid < other_gene_organism_taxonid
                                    {
                                        ref_details.ortholog_annotations.push(ortholog_annotation);
                                    }
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
                                    };
                                other_gene_details.ortholog_annotations.push(ortholog_annotation.clone());
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
                                        ref_details.ortholog_annotations.push(ortholog_annotation);
                                    }
                                }
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
    fn matching_ext_config(&self, annotation_termid: &str,
                           rel_type_name: &str) -> Option<ExtensionDisplayNames> {
        let ext_configs = &self.config.extension_display_names;

        if let Some(annotation_term_details) = self.terms.get(annotation_termid) {
            for ext_config in ext_configs {
                if ext_config.rel_name == rel_type_name {
                    if let Some(if_descendant_of) = ext_config.if_descendant_of.clone() {
                        if annotation_term_details.interesting_parents.contains(&if_descendant_of) {
                            return Some((*ext_config).clone());
                        }
                    } else {
                        return Some((*ext_config).clone());
                    }
                }
            }
        } else {
            panic!("can't find details for term: {}\n", annotation_termid);
        }

        None
    }

    // create and returns any TargetOfAnnotations implied by the extension
    fn make_target_of_for_ext(&self, cv_name: &str,
                              genes: &[RcString],
                              maybe_genotype_uniquename: &Option<RcString>,
                              reference_uniquename: &Option<RcString>,
                              annotation_termid: &str,
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
            if priority_config.get(annotation.ext_rel_display_name.as_ref()).is_none() {
                eprintln!(r#"No priority configured for "{}" (from {})"#,
                          annotation.ext_rel_display_name, gene_details.uniquename);
            }
        }

        let cmp_fn = |a: &TargetOfAnnotation, b: &TargetOfAnnotation| {
            let a_rel_name = a.ext_rel_display_name.as_ref();
            let a_pri = priority_config.get(a_rel_name).unwrap_or(&0);
            let b_rel_name = b.ext_rel_display_name.as_ref();
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
                            .get(&annotation_id).expect("can't find OntAnnotationDetail");
                        if let Some(ref genotype_uniquename) = annotation.genotype {
                            let genotype = &self.genotypes[genotype_uniquename];

                            if genotype.expressed_alleles.len() > 1 {
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
                            target_of_annotations
                                .entry(target_gene_uniquename.clone())
                                .or_insert_with(HashSet::new)
                                .insert(new_annotation);
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
        let some_null = Some(RcString::from("Null"));

        let mut gene_statuses = HashMap::new();

        let condition_string =
            |condition_ids: HashSet<RcString>| {
                let mut ids_vec: Vec<RcString> = condition_ids.iter().cloned().collect();
                ids_vec.sort();
                RcString::from(&ids_vec.join(" "))
            };

        let viable_termid = &self.config.viability_terms.viable;
        let inviable_termid = &self.config.viability_terms.inviable;

        for (gene_uniquename, gene_details) in &mut self.genes {
            let mut new_status = DeletionViability::Unknown;

            if let Some(single_allele_term_annotations) =
                gene_details.cv_annotations.get("single_allele_phenotype") {
                    let mut viable_conditions: HashMap<RcString, TermId> = HashMap::new();
                    let mut inviable_conditions: HashMap<RcString, TermId> = HashMap::new();

                    for term_annotation in single_allele_term_annotations {
                        'ANNOTATION: for annotation_id in &term_annotation.annotations {
                            let annotation = self.annotation_details
                                .get(&annotation_id).expect("can't find OntAnnotationDetail");

                            let genotype_uniquename = annotation.genotype.as_ref().unwrap();

                            let genotype = &self.genotypes[genotype_uniquename];
                            let expressed_allele = &genotype.expressed_alleles[0];
                            let allele = &self.alleles[&expressed_allele.allele_uniquename];
                            if allele.allele_type != "deletion" &&
                                expressed_allele.expression != some_null {
                                continue 'ANNOTATION;
                            }

                            let term = &self.terms[&term_annotation.term];
                            let interesting_parents = &term.interesting_parents;
                            let conditions_as_string =
                                condition_string(annotation.conditions.clone());
                            if interesting_parents.contains(viable_termid) ||
                                *viable_termid == term_annotation.term {
                                    viable_conditions.insert(conditions_as_string,
                                                             term_annotation.term.clone());
                                } else {
                                    if interesting_parents.contains(inviable_termid) ||
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

                            let viable_conditions_set: HashSet<RcString> =
                                viable_conditions.keys().cloned().collect();
                            let inviable_conditions_set: HashSet<RcString> =
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
            |subset_termid: &str, test_termid: &str| {
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

        let in_archaea = "conserved in archaea";
        let in_bacteria = "conserved in bacteria";
        let in_fungi_only = "conserved in fungi only";
        let in_metazoa = "conserved in metazoa";
        let pombe_specific = "Schizosaccharomyces pombe specific";
        let schizo_specific = "Schizosaccharomyces specific";

        let names = vec![in_archaea, in_bacteria, in_fungi_only, in_metazoa,
                         pombe_specific, schizo_specific];

        for name in names {
            if let Some(termid) = self.term_ids_by_name.get(name) {
                term_name_map.insert(termid.clone(), name.clone());
            } else {
                eprintln!("configuration error: can't find {} in term_ids_by_name map", name);
                eprintln!("skipping taxonomic distribution");
                return;
            }
        }

        'GENE:
        for gene_details in self.genes.values_mut() {
            let mut dist_names = HashSet::new();

            if let Some(species_dists) = gene_details.cv_annotations.get("species_dist") {
                for ont_term_annotations in species_dists {
                    let ref term = ont_term_annotations.term;
                    if let Some(term_name) = term_name_map.get(term) {
                        dist_names.insert(term_name.clone());
                    }
                }
            }

            if (dist_names.contains(in_archaea) || dist_names.contains(in_bacteria))
                && !dist_names.contains(in_metazoa) {
                    gene_details.taxonomic_distribution =
                        Some(RcString::from("fungi and prokaryotes"));
                    continue 'GENE;
                }
            if dist_names.contains(in_metazoa) &&
                !((dist_names.contains(in_archaea) || dist_names.contains(in_bacteria))
                  && dist_names.contains(in_metazoa)) {
                    gene_details.taxonomic_distribution =
                        Some(RcString::from("eukaryotes only, fungi and metazoa"));
                    continue 'GENE;
                }


            if (dist_names.contains(in_archaea) || dist_names.contains(in_bacteria)) &&
                dist_names.contains(in_metazoa) {
                    gene_details.taxonomic_distribution =
                        Some(RcString::from("eukaryotes and prokaryotes"));
                    continue 'GENE;
                }

            if dist_names.contains(in_fungi_only) {
                gene_details.taxonomic_distribution = Some(RcString::from("fungi only"));
                continue 'GENE;
            }

            if dist_names.contains(pombe_specific) {
                gene_details.taxonomic_distribution = Some(RcString::from("S. pombe specific"));
                continue 'GENE;
            }

            if dist_names.contains(schizo_specific) {
                gene_details.taxonomic_distribution = Some(RcString::from("Schizos. specific"));
                continue 'GENE;
            }

            if let Some(ref characterisation_status) = gene_details.characterisation_status {
                if characterisation_status == "dubious" {
                    gene_details.taxonomic_distribution = Some(RcString::from("dubious"));
                    continue 'GENE;
                }
            }

            if gene_details.feature_type != "mRNA gene" {
                gene_details.taxonomic_distribution = Some(RcString::from("not curated"));
                continue 'GENE;
            }

            gene_details.taxonomic_distribution = Some(RcString::from("other"));
        }
    }

    fn make_gene_short_map(&self) -> IdGeneShortMap {
        let mut ret_map = HashMap::new();

        for gene_uniquename in self.genes.keys() {
            ret_map.insert(gene_uniquename.clone(),
                           make_gene_short(&self.genes, &gene_uniquename));
        }

        ret_map
    }

    fn make_residue_extensions(&mut self) {
        for annotation in self.annotation_details.values_mut() {
            if let Some(ref residue) = annotation.residue {
                let display_name = RcString::from("modified residue");
                let residue_range_part = ExtPart {
                    rel_type_id: None,
                    rel_type_name: display_name.clone(),
                    rel_type_display_name: display_name,
                    ext_range: ExtRange::SummaryModifiedResidues(vec![residue.clone()]),
                };
                annotation.extension.insert(0, residue_range_part);
            }
        }
    }

    fn make_all_cv_summaries(&mut self) {
        let gene_short_map = self.make_gene_short_map();

        for term_details in self.terms.values_mut() {
            make_cv_summaries(term_details, self.config, &self.children_by_termid,
                              true, true, &gene_short_map, &self.annotation_details);
        }

        for gene_details in self.genes.values_mut() {
            make_cv_summaries(gene_details, &self.config, &self.children_by_termid,
                              false, true, &gene_short_map, &self.annotation_details);
        }

        for genotype_details in self.genotypes.values_mut() {
            make_cv_summaries(genotype_details, &self.config, &self.children_by_termid,
                              false, false, &gene_short_map, &self.annotation_details);
        }

        for reference_details in self.references.values_mut() {
            make_cv_summaries( reference_details, &self.config, &self.children_by_termid,
                              true, true, &gene_short_map, &self.annotation_details);
        }
    }

    fn process_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name != POMBASE_ANN_EXT_TERM_CV_NAME {
                let cv_config = self.config.cv_config_by_name(&cvterm.cv.name);
                let annotation_feature_type = cv_config.feature_type.clone();

                let mut xrefs = HashMap::new();

                for (source_name, source_config) in cv_config.source_config {
                    let mut maybe_xref_id = None;
                    if let Some(ref term_xref_id_prop) = source_config.id_prop {
                        for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                            if cvtermprop.prop_type.name == *term_xref_id_prop {
                                maybe_xref_id = Some(cvtermprop.value.clone());
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
                        }
                    }).collect::<Vec<_>>();

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
                                      interesting_parents: HashSet::new(),
                                      in_subsets: HashSet::new(),
                                      termid: cvterm.termid(),
                                      synonyms,
                                      definition: cvterm.definition.clone(),
                                      direct_ancestors: vec![],
                                      secondary_identifiers,
                                      genes_annotated_with: HashSet::new(),
                                      is_obsolete: cvterm.is_obsolete,
                                      single_allele_genotype_uniquenames: HashSet::new(),
                                      cv_annotations: HashMap::new(),
                                      genes_by_uniquename: HashMap::new(),
                                      genotypes_by_uniquename: HashMap::new(),
                                      alleles_by_uniquename: HashMap::new(),
                                      references_by_uniquename: HashMap::new(),
                                      terms_by_termid: HashMap::new(),
                                      annotation_details: HashMap::new(),
                                      gene_count: 0,
                                      genotype_count: 0,
                                      xrefs,
                                  });
                self.term_ids_by_name.insert(cvterm.name.clone(), cvterm.termid());
            }
        }
    }

    fn get_ext_rel_display_name(&self, annotation_termid: &str,
                                ext_rel_name: &str) -> RcString {
        if let Some(ext_conf) = self.matching_ext_config(annotation_termid, ext_rel_name) {
            ext_conf.display_name.clone()
        } else {
            RcString::from(&str::replace(&ext_rel_name, "_", " "))
        }
    }

    fn process_extension_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                    if (*cvtermprop).prop_type.name.starts_with(ANNOTATION_EXT_REL_PREFIX) {
                        let ext_rel_name_str =
                            &(*cvtermprop).prop_type.name[ANNOTATION_EXT_REL_PREFIX.len()..];
                        let ext_rel_name = RcString::from(ext_rel_name_str);
                        let ext_range = (*cvtermprop).value.clone();
                        let range: ExtRange = if ext_range.starts_with("SP") {
                            if let Some(captures) = PROMOTER_RE.captures(&ext_range) {
                                let gene_uniquename = RcString::from(&captures["gene"]);
                                ExtRange::Promoter(gene_uniquename)
                            } else {
                                ExtRange::Gene(ext_range.clone())
                            }
                        } else {
                            ExtRange::Misc(ext_range)
                        };

                        if let Some(base_termid) =
                            self.base_term_of_extensions.get(&cvterm.termid()) {
                                let rel_type_display_name =
                                    self.get_ext_rel_display_name(&base_termid, &ext_rel_name);

                                let rel_type_id =
                                    self.term_ids_by_name.get(&ext_rel_name)
                                    .map(|t| t.clone());

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
                                if object_termid.starts_with("PR:") {
                                    ExtRange::GeneProduct(object_termid)
                                } else {
                                    ExtRange::Term(object_termid)
                                };

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

            if let Some(ref mut gene_details) = self.genes.get_mut(&feature.uniquename) {
                gene_details.synonyms.push(SynonymDetails {
                    name: synonym.name.clone(),
                    synonym_type: synonym.synonym_type.name.clone()
                });
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

    fn make_genotype_short(&self, genotype_display_name: &str) -> GenotypeShort {
        if let Some(ref details) = self.genotypes.get(genotype_display_name) {
            GenotypeShort {
                display_uniquename: details.display_uniquename.clone(),
                name: details.name.clone(),
                expressed_alleles: details.expressed_alleles.clone(),
            }
        } else {
            panic!("can't find genotype {}", genotype_display_name);
        }
    }

    fn make_allele_short(&self, allele_uniquename: &str) -> AlleleShort {
        self.alleles[allele_uniquename].clone()
    }

    // process feature properties stored as cvterms,
    // eg. characterisation_status and product
    fn process_props_from_feature_cvterms(&mut self) {
        for feature_cvterm in &self.raw.feature_cvterms {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;

            let gene_uniquenames_vec: Vec<GeneUniquename> =
                if cvterm.cv.name == "PomBase gene products" {
                    if feature.feat_type.name == "polypeptide" {
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
                    } else {
                        if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                            if let Some(gene_uniquename) =
                                self.genes_of_transcripts.get(&feature.uniquename) {
                                    vec![gene_uniquename.clone()]
                                } else {
                                    vec![]
                                }
                        } else {
                            if feature.feat_type.name == "gene" {
                                vec![feature.uniquename.clone()]
                            } else {
                                vec![]
                            }
                        }
                    }
                } else {
                    vec![]
                };
            for gene_uniquename in &gene_uniquenames_vec {
                self.add_gene_product(&gene_uniquename, &cvterm.name);
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

    fn get_gene_prod_extension(&self, prod_value: RcString) -> ExtPart {
        let ext_range =
            if prod_value.starts_with("PR:") {
                ExtRange::GeneProduct(prod_value)
            } else {
                ExtRange::Misc(prod_value)
            };

        let active_form = "active_form";

        ExtPart {
            rel_type_id: None,
            rel_type_name: active_form.into(),
            rel_type_display_name: "active form".into(),
            ext_range,
        }
    }

    // return a fake extension for "with" properties on protein binding annotations
    fn get_with_extension(&self, with_value: &RcString) -> ExtPart {
        let ext_range =
            if with_value.starts_with("SP%") {
                ExtRange::Gene(with_value.clone())
            } else {
                if with_value.starts_with("PomBase:SP") {
                    let gene_uniquename =
                        RcString::from(&with_value[8..]);
                    ExtRange::Gene(gene_uniquename)
                } else {
                    if with_value.to_lowercase().starts_with("pfam:") {
                        ExtRange::Domain(with_value.clone())
                    } else {
                        ExtRange::Misc(with_value.clone())
                    }
                }
            };

        // a with property on a protein binding (GO:0005515) is
        // displayed as a binds extension
        // https://github.com/pombase/website/issues/108
        ExtPart {
            rel_type_id: None,
            rel_type_name: "binds".into(),
            rel_type_display_name: "binds".into(),
            ext_range,
        }
    }

    fn make_with_or_from_value(&self, with_or_from_value: &RcString) -> WithFromValue {
        let db_prefix_patt = RcString::from("^") + DB_NAME + ":";
        let re = Regex::new(&db_prefix_patt).unwrap();
        let gene_uniquename = RcString::from(&re.replace_all(&with_or_from_value, ""));
        if self.genes.contains_key(&gene_uniquename) {
            let gene_short = self.make_gene_short(&gene_uniquename);
            WithFromValue::Gene(gene_short)
        } else {
            if self.terms.get(with_or_from_value).is_some() {
                WithFromValue::Term(self.make_term_short(with_or_from_value))
            } else {
                WithFromValue::Identifier(with_or_from_value.clone())
            }
        }
    }

    // add the with value as a fake extension if the cvterm is_a protein binding,
    // otherwise return the value
    fn make_with_extension(&self, termid: &RcString, evidence_code: Option<RcString>,
                           extension: &mut Vec<ExtPart>,
                           with_value: &RcString) -> Option<WithFromValue> {
        let base_termid =
            match self.base_term_of_extensions.get(termid) {
                Some(base_termid) => base_termid.clone(),
                None => termid.clone(),
            };

        let base_term_short = self.make_term_short(&base_termid);

        if evidence_code.is_some() &&
            evidence_code.unwrap() == "IPI" &&
            // add new IDs to the interesting_parents config
            (base_term_short.termid == "GO:0005515" ||
             base_term_short.interesting_parents.contains("GO:0005515") ||
             base_term_short.termid == "GO:0003723" ||
             base_term_short.interesting_parents.contains("GO:0003723")) {
                extension.push(self.get_with_extension(with_value));
            } else {
                return Some(self.make_with_or_from_value(with_value));
            }
        None
    }

    // process annotation
    fn process_feature_cvterms(&mut self) {
        for feature_cvterm in &self.raw.feature_cvterms {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;

            let termid = cvterm.termid();

            let mut extension = vec![];

            if cvterm.cv.name == "PomBase gene characterisation status" ||
                cvterm.cv.name == "PomBase gene products" ||
                cvterm.cv.name == "name_description" {
                    continue;
                }

            let publication = &feature_cvterm.publication;
            let mut extra_props: HashMap<RcString, RcString> = HashMap::new();
            let mut conditions: HashSet<TermId> = HashSet::new();
            let mut withs: HashSet<WithFromValue> = HashSet::new();
            let mut froms: HashSet<WithFromValue> = HashSet::new();
            let mut qualifiers: Vec<Qualifier> = vec![];
            let mut date: Option<RcString> = None;
            let mut assigned_by: Option<RcString> = None;
            let mut evidence: Option<RcString> = None;
            let mut genotype_background: Option<RcString> = None;
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
                    "residue" | "scale" |
                    "quant_gene_ex_copies_per_cell" |
                    "quant_gene_ex_avg_copies_per_cell" => {
                        if let Some(value) = prop.value.clone() {
                            extra_props.insert(prop.type_name().clone(), value);
                        }
                    },
                    "condition" =>
                        if let Some(value) = prop.value.clone() {
                            if value.contains(":") {
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
                    "with" => {
                        if let Some(ref with_value) = prop.value.clone() {
                            if let Some(with_gene_short) =
                                self.make_with_extension(&termid, evidence.clone(),
                                                         &mut extension, with_value) {
                                    withs.insert(with_gene_short);
                                }
                        }
                    },
                    "from" => {
                        if let Some(value) = prop.value.clone() {
                            froms.insert(self.make_with_or_from_value(&value));
                        }
                    },
                    "gene_product_form_id" => {
                        if let Some(value) = prop.value.clone() {
                            extension.push(self.get_gene_prod_extension(value));
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
                        let expressed_alleles =
                            &self.alleles_of_genotypes[&feature.uniquename];
                        let genotype_display_name =
                            make_genotype_display_name(&expressed_alleles, &self.alleles);
                        maybe_genotype_uniquename = Some(genotype_display_name.clone());
                        genotype_background =
                            self.genotype_backgrounds.get(&feature.uniquename)
                            .map(|s| s.clone());
                        expressed_alleles.iter()
                            .map(|expressed_allele| {
                                let allele_short =
                                    self.make_allele_short(&expressed_allele.allele_uniquename);
                                allele_short.gene_uniquename.clone()
                            })
                            .collect()
                    },
                    "gene" | "pseudogene" => {
                        vec![feature.uniquename.clone()]
                    },
                    _ =>
                        if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                            if let Some(gene_uniquename) =
                                self.genes_of_transcripts.get(&feature.uniquename) {
                                    vec![gene_uniquename.clone()]
                                } else {
                                    vec![]
                                }
                        } else {
                            panic!("can't handle annotation on feature type: {:?} {}",
                                   &feature.feat_type.name, &feature.uniquename);
                        }
                };

            gene_uniquenames_vec.dedup();

            gene_uniquenames_vec =
                gene_uniquenames_vec.iter().map(|gene_uniquename: &RcString| {
                    self.make_gene_short(&gene_uniquename).uniquename
                }).collect();

            let reference_uniquename =
                if publication.uniquename == "null" {
                    None
                } else {
                    Some(publication.uniquename.clone())
                };

            let mut extra_props_clone = extra_props.clone();
            let copies_per_cell = extra_props_clone.remove("quant_gene_ex_copies_per_cell");
            let avg_copies_per_cell = extra_props_clone.remove("quant_gene_ex_avg_copies_per_cell");
            let scale = extra_props_clone.remove("scale");
            let gene_ex_props =
                if copies_per_cell.is_some() || avg_copies_per_cell.is_some() {
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
                reference: reference_uniquename,
                genotype: maybe_genotype_uniquename,
                genotype_background,
                withs,
                froms,
                residue: extra_props_clone.remove("residue"),
                gene_ex_props,
                qualifiers,
                evidence,
                conditions,
                extension,
                date,
                assigned_by,
                throughput,
            };

            self.add_annotation(cvterm.borrow(), feature_cvterm.is_not,
                                annotation_detail);
        }
    }

    fn make_term_annotations(&self, termid: &RcString, detail_ids: &[OntAnnotationId],
                             is_not: bool)
                       -> Vec<(CvName, OntTermAnnotations)> {
        let term_details = &self.terms[termid];

        let cv_name = term_details.cv_name.clone();

        match cv_name.as_ref() {
            "gene_ex" => {
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
                        get(&annotation_id).expect("can't find OntAnnotationDetail");

                    if annotation.gene_ex_props.is_some() {
                        quant_annotations.annotations.push(*annotation_id)
                    } else {
                        qual_annotations.annotations.push(*annotation_id)
                    }
                }

                let mut return_vec = vec![];

                if !qual_annotations.annotations.is_empty() {
                    return_vec.push((RcString::from("qualitative_gene_expression"),
                                     qual_annotations));
                }

                if !quant_annotations.annotations.is_empty() {
                    return_vec.push((RcString::from("quantitative_gene_expression"),
                                     quant_annotations));
                }

                return_vec
            },
            "fission_yeast_phenotype" => {
                let mut single_allele =
                    OntTermAnnotations {
                        term: termid.clone(),
                        is_not,
                        rel_names: HashSet::new(),
                        annotations: vec![],
                        summary: None,
                    };
                let mut multi_allele =
                    OntTermAnnotations {
                        term: termid.clone(),
                        is_not,
                        rel_names: HashSet::new(),
                        annotations: vec![],
                        summary: None,
                    };

                for annotation_id in detail_ids {
                    let annotation = self.annotation_details.
                        get(&annotation_id).expect("can't find OntAnnotationDetail");

                    let genotype_uniquename = annotation.genotype.as_ref().unwrap();

                    if let Some(genotype_details) = self.genotypes.get(genotype_uniquename) {
                        if genotype_details.expressed_alleles.len() == 1 {
                            single_allele.annotations.push(*annotation_id);
                        } else {
                            multi_allele.annotations.push(*annotation_id);
                        }
                    } else {
                        panic!("can't find genotype details for {}\n", genotype_uniquename);
                    }
                }

                let mut return_vec = vec![];

                if !single_allele.annotations.is_empty() {
                    return_vec.push((RcString::from("single_allele_phenotype"),
                                     single_allele));
                }

                if !multi_allele.annotations.is_empty() {
                    return_vec.push((RcString::from("multi_allele_phenotype"),
                                     multi_allele));
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
        let mut ref_annotation_by_term: HashMap<RcString, HashMap<TermId, Vec<OntAnnotationId>>> =
            HashMap::new();

        let mut ont_annotations = vec![];

        for (termid, annotations) in ont_annotation_map {
            let mut sorted_annotations = annotations.clone();

            if !is_not {
                let cv_config = {
                    let term = &self.terms[termid];
                    &self.config.cv_config_by_name(&term.cv_name)
                };

                {
                    let cmp_detail_with_maps =
                        |id1: &i32, id2: &i32| {
                            let annotation1 = self.annotation_details.
                                get(&id1).expect("can't find OntAnnotationDetail");
                            let annotation2 = self.annotation_details.
                                get(&id2).expect("can't find OntAnnotationDetail");

                            let result =
                                cmp_ont_annotation_detail(cv_config,
                                                          annotation1, annotation2, &self.genes,
                                                          &self.genotypes,
                                                          &self.terms);
                            result.unwrap_or_else(|err| panic!("error from cmp_ont_annotation_detail: {}", err))
                        };

                    sorted_annotations.sort_by(cmp_detail_with_maps);
                }

                let new_annotations =
                    self.make_term_annotations(&termid, &sorted_annotations, is_not);

                if let Some(ref mut term_details) = self.terms.get_mut(termid) {
                    for (cv_name, new_annotation) in new_annotations {
                        term_details.cv_annotations.entry(cv_name.clone())
                            .or_insert_with(Vec::new)
                            .push(new_annotation);
                    }
                } else {
                    panic!("missing termid: {}\n", termid);
                }
            }

            for annotation_id in sorted_annotations {
                let annotation = self.annotation_details.
                    get(&annotation_id).expect("can't find OntAnnotationDetail");

                for gene_uniquename in &annotation.genes {
                    gene_annotation_by_term.entry(gene_uniquename.clone())
                        .or_insert_with(HashMap::new)
                        .entry(termid.clone())
                        .or_insert_with(|| vec![])
                        .push(annotation_id);
                }

                if let Some(ref genotype_uniquename) = annotation.genotype {
                    let existing =
                        genotype_annotation_by_term.entry(genotype_uniquename.clone())
                        .or_insert_with(HashMap::new)
                        .entry(termid.clone())
                        .or_insert_with(|| vec![]);
                    if !existing.contains(&annotation_id) {
                        existing.push(annotation_id);
                    }
                }

                if let Some(reference_uniquename) = annotation.reference.clone() {
                    ref_annotation_by_term.entry(reference_uniquename)
                        .or_insert_with(HashMap::new)
                        .entry(termid.clone())
                        .or_insert_with(|| vec![])
                        .push(annotation_id);
                }

                for condition_termid in &annotation.conditions {
                    let cv_name =
                        if let Some(ref term_details) = self.terms.get(condition_termid) {
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
                            .annotations.push(annotation_id);
                    }
                }

                /*

                disable for now:

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
                                .entry(RcString::from(&extension_cv_name))
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
                    annotation.genes.iter().map(|uniquename: &RcString| {
                        self.make_gene_short(uniquename)
                    }).collect::<HashSet<_>>();

                let reference_short =
                    annotation.reference.as_ref().map_or(None, |uniquename: &RcString| {
                        make_reference_short(&self.references, uniquename)
                    });

                let genotype_short =
                    annotation.genotype.as_ref().map(|uniquename: &RcString| {
                        self.make_genotype_short(uniquename)
                    });

                let conditions =
                    annotation.conditions.iter().map(|termid| {
                        self.make_term_short(termid)
                    }).collect::<HashSet<_>>();

                if gene_short_list.len() == 0 {
                    panic!("no genes for {:?}", &annotation);
                }

                let ont_annotation = OntAnnotation {
                    term_short: self.make_term_short(&termid),
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

            term1.cmp(&term2)
        };

        for (gene_uniquename, term_annotation_map) in &gene_annotation_by_term {
            for (termid, details) in term_annotation_map {
                let new_annotations =
                    self.make_term_annotations(&termid, &details, is_not);

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
                    self.make_term_annotations(&termid, &details, is_not);

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
                    self.make_term_annotations(&termid, &details, is_not);

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

    // return true if the term could or should appear in the interesting_parents
    // field of the TermDetails and TermShort structs
    fn is_interesting_parent(&self, termid: &str, rel_name: &str) -> bool {
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
                    !HAS_PART_CV_NAMES.contains(&subject_term_details.cv_name.as_str()) {
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
                let mut all_rel_names: HashSet<RcString> = HashSet::new();
                for (annotation_id, rel_names) in source_annotations_map {
                    new_annotations.push(annotation_id);
                    for rel_name in rel_names {
                        all_rel_names.insert(rel_name);
                    }
                }

                let dest_cv_config = &self.config.cv_config_by_name(&dest_cv_name);

                {
                    let cmp_detail_with_genotypes =
                        |id1: &i32, id2: &i32| {
                            let annotation1 = self.annotation_details.
                                get(&id1).expect("can't find OntAnnotationDetail");
                            let annotation2 = self.annotation_details.
                                get(&id2).expect("can't find OntAnnotationDetail");
                            let result =
                                cmp_ont_annotation_detail(dest_cv_config,
                                                          annotation1, annotation2, &self.genes,
                                                          &self.genotypes, &self.terms);
                            result.unwrap_or_else(|err| {
                                panic!("error from cmp_ont_annotation_detail: {} with terms: {} and {}",
                                       err, source_termid, dest_termid)
                            })
                        };

                    new_annotations.sort_by(cmp_detail_with_genotypes);
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

        let ont_term_cmp = |ont_term_1: &OntTermAnnotations, ont_term_2: &OntTermAnnotations| {
            if !ont_term_1.is_not && ont_term_2.is_not {
                return Ordering::Less;
            }
            if ont_term_1.is_not && !ont_term_2.is_not {
                return Ordering::Greater;
            }
            let term1 = &term_names[&ont_term_1.term];
            let term2 = &term_names[&ont_term_2.term];

            term1.cmp(&term2)
        };

        for term_details in self.terms.values_mut() {
            for term_annotations in term_details.cv_annotations.values_mut() {
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
            export_prog_name: RcString::from(PKG_NAME),
            export_prog_version: RcString::from(VERSION),
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
                            get(&annotation_id).expect("can't find OntAnnotationDetail");

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
                            is_multi: genotype.expressed_alleles.len() > 1,
                            conditions,
                            alleles: vec![],
                        };
                        for allele in &genotype.expressed_alleles {
                            let allele_uniquename = &allele.allele_uniquename;
                            let allele_short =
                                self.alleles.get(allele_uniquename).expect("Can't find allele");
                            let allele_gene_uniquename =
                                allele_short.gene_uniquename.clone();
                            let allele_details = APIAlleleDetails {
                                gene: allele_gene_uniquename,
                                allele_type: allele_short.allele_type.clone(),
                                expression: allele.expression.clone(),
                            };
                            api_annotation.alleles.push(allele_details);
                        }
                        app_genotype_annotation
                            .entry(term_details.termid.clone())
                            .or_insert_with(|| vec![])
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

        for transcript in &gene_details.transcripts {
            if let Some(ref protein) = transcript.protein {
                molecular_weight = Some((100.0 * protein.molecular_weight).round() / 100.0);
                if protein.sequence.ends_with("*") {
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

    fn make_gene_query_go_data(&self, gene_details: &GeneDetails, term_config: &Vec<TermId>,
                               cv_name: &str) -> Option<GeneQueryTermData>
    {
        let component_term_annotations =
            gene_details.cv_annotations.get(cv_name);

        if component_term_annotations.is_none() {
            return None;
        }

        let in_component = |check_termid: &str| {
            for term_annotation in component_term_annotations.unwrap() {
                let maybe_term_details = self.terms.get(&term_annotation.term);

                let term_details =
                    maybe_term_details .unwrap_or_else(|| {
                        panic!("can't find TermDetails for {}", &term_annotation.term)
                    });

                let interesting_parents = &term_details.interesting_parents;

                if !term_annotation.is_not &&
                    (term_annotation.term == check_termid ||
                     interesting_parents.contains(check_termid))
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
                                             "cellular_component");
            let go_process_superslim =
                self.make_gene_query_go_data(gene_details, &process_terms,
                                             "biological_process");
            let go_function =
                self.make_gene_query_go_data(gene_details, &function_terms,
                                             "molecular_function");

            let tmm =
                if gene_details.feature_type == "mRNA gene" {
                    if gene_details.tm_domain_coords.len() == 0 {
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
                let gene_summary = self.make_api_gene_summary(&gene_uniquename);
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
                        if gene_uniquename == &interaction_annotation.gene_uniquename {
                            interaction_annotation.interactor_uniquename.clone()
                        } else {
                            interaction_annotation.gene_uniquename.clone()
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
                                        term_details.genes_annotated_with.clone());
                }
            }

            terms_for_api.insert(termid.clone(), term_details);
        }

        let term_subsets = self.term_subsets.clone();
        let gene_subsets = self.gene_subsets.clone();

        APIMaps {
            gene_summaries,
            gene_query_data_map,
            termid_genes,
            termid_genotype_annotation,
            term_summaries,
            genes: self.genes,
            gene_name_gene_map,
            genotypes: self.genotypes,
            terms: terms_for_api,
            interactors_of_genes,
            references: self.references,
            other_features: self.other_features,
            annotation_details: self.annotation_details,
            chromosomes: self.chromosomes,
            term_subsets,
            gene_subsets,
       }
    }

    fn add_cv_annotations_to_maps(&self,
                                  identifier: &RcString,
                                  cv_annotations: &OntAnnotationMap,
                                  seen_references: &mut HashMap<RcString, ReferenceShortOptionMap>,
                                  seen_genes: &mut HashMap<RcString, GeneShortOptionMap>,
                                  seen_genotypes: &mut HashMap<RcString, GenotypeShortMap>,
                                  seen_alleles: &mut HashMap<RcString, AlleleShortMap>,
                                  seen_terms: &mut HashMap<RcString, TermShortOptionMap>) {
        for feat_annotations in cv_annotations.values() {
            for feat_annotation in feat_annotations.iter() {
                self.add_term_to_hash(seen_terms, identifier,
                                      &feat_annotation.term);

                for annotation_detail_id in &feat_annotation.annotations {
                    let annotation_detail = self.annotation_details.
                        get(&annotation_detail_id).expect("can't find OntAnnotationDetail");

                    self.add_ref_to_hash(seen_references, identifier,
                                         &annotation_detail.reference);
                    for condition_termid in &annotation_detail.conditions {
                        self.add_term_to_hash(seen_terms, identifier,
                                              &condition_termid);
                    }
                    for ext_part in &annotation_detail.extension {
                        match ext_part.ext_range {
                            ExtRange::Term(ref range_termid) |
                            ExtRange::GeneProduct(ref range_termid) =>
                                self.add_term_to_hash(seen_terms, identifier,
                                                      range_termid),
                            ExtRange::Gene(ref gene_uniquename) |
                            ExtRange::Promoter(ref gene_uniquename) =>
                                self.add_gene_to_hash(seen_genes, identifier,
                                                      &gene_uniquename),
                            _ => {},
                        }
                    }
                    if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  identifier, genotype_uniquename);
                    }
                }
            }
        }
    }

    fn set_term_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_terms) = get_maps();

        let mut genes_annotated_with_map: HashMap<TermId, HashSet<GeneUniquename>> =
            HashMap::new();

        for (termid, term_details) in &self.terms {
            for (cv_name, term_annotations) in &term_details.cv_annotations {
                for term_annotation in term_annotations {
                    self.add_term_to_hash(&mut seen_terms, termid,
                                          &term_annotation.term);
                    for annotation_detail_id in &term_annotation.annotations {
                        let annotation_detail = self.annotation_details
                            .get(&annotation_detail_id).expect("can't find OntAnnotationDetail");

                        for gene_uniquename in &annotation_detail.genes {
                            self.add_gene_to_hash(&mut seen_genes, termid,
                                                  gene_uniquename);
                            if !cv_name.starts_with("extension:") {
                                // prevent extension annotations from appearing
                                // in the normal query builder searches
                                genes_annotated_with_map
                                    .entry(termid.clone()).or_insert_with(HashSet::new)
                                    .insert(gene_uniquename.clone());
                            }
                        }
                        self.add_ref_to_hash(&mut seen_references, termid,
                                             &annotation_detail.reference);
                        for condition_termid in &annotation_detail.conditions {
                            self.add_term_to_hash(&mut seen_terms, termid,
                                                  condition_termid);
                        }
                        for ext_part in &annotation_detail.extension {
                            match ext_part.ext_range {
                                ExtRange::Term(ref range_termid) |
                                ExtRange::GeneProduct(ref range_termid) =>
                                    self.add_term_to_hash(&mut seen_terms, termid,
                                                          range_termid),
                                ExtRange::Gene(ref gene_uniquename) |
                                ExtRange::Promoter(ref gene_uniquename) =>
                                    self.add_gene_to_hash(&mut seen_genes, termid,
                                                          gene_uniquename),
                                _ => {},
                            }
                        }
                        if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                            self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                      &mut seen_genes, &termid,
                                                      genotype_uniquename);
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
            if let Some(terms) = seen_terms.remove(termid) {
                term_details.terms_by_termid = terms;
            }
            if let Some(gene_uniquename_set) = genes_annotated_with_map.remove(termid) {
                term_details.genes_annotated_with = gene_uniquename_set;
            }
        }
    }

    fn set_gene_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_terms) = get_maps();

        {
            for (gene_uniquename, gene_details) in &self.genes {
                self.add_cv_annotations_to_maps(&gene_uniquename,
                                                &gene_details.cv_annotations,
                                                &mut seen_references,
                                                &mut seen_genes,
                                                &mut seen_genotypes,
                                                &mut seen_alleles,
                                                &mut seen_terms);

                let interaction_iter =
                    gene_details.physical_interactions.iter().chain(&gene_details.genetic_interactions);
                for interaction in interaction_iter {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename,
                                         &interaction.reference_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &interaction.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename,
                                          &interaction.interactor_uniquename);
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
                    let ref target_of_gene = target_of_annotation.gene;
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename, target_of_gene);
                    if let Some(ref annotation_genotype_uniquename) = target_of_annotation.genotype_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
                                                  gene_uniquename,
                                                  &annotation_genotype_uniquename)
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
            if let Some(terms) = seen_terms.remove(gene_uniquename) {
                gene_details.terms_by_termid = terms;
            }
        }
    }

    fn set_genotype_details_maps(&mut self) {
        let (mut seen_references, mut seen_genes, mut seen_genotypes,
             mut seen_alleles, mut seen_terms) = get_maps();

        for (genotype_uniquename, genotype_details) in &self.genotypes {
            self.add_cv_annotations_to_maps(&genotype_uniquename,
                                            &genotype_details.cv_annotations,
                                            &mut seen_references,
                                            &mut seen_genes,
                                            &mut seen_genotypes,
                                            &mut seen_alleles,
                                            &mut seen_terms);
        }

        for (genotype_uniquename, genotype_details) in &mut self.genotypes {
            if let Some(references) = seen_references.remove(genotype_uniquename) {
                genotype_details.references_by_uniquename = references;
            }
            if let Some(alleles) = seen_alleles.remove(genotype_uniquename) {
                genotype_details.alleles_by_uniquename = alleles;
            }
            if let Some(genotypes) = seen_genes.remove(genotype_uniquename) {
                genotype_details.genes_by_uniquename = genotypes;
            }
            if let Some(terms) = seen_terms.remove(genotype_uniquename) {
                genotype_details.terms_by_termid = terms;
            }
        }
    }

    fn set_reference_details_maps(&mut self) {
        let mut seen_genes: HashMap<RcString, GeneShortOptionMap> = HashMap::new();

        type GenotypeShortMap = HashMap<GenotypeUniquename, GenotypeShort>;
        let mut seen_genotypes: HashMap<ReferenceUniquename, GenotypeShortMap> = HashMap::new();

        type AlleleShortMap = HashMap<AlleleUniquename, AlleleShort>;
        let mut seen_alleles: HashMap<TermId, AlleleShortMap> = HashMap::new();

        let mut seen_terms: HashMap<GeneUniquename, TermShortOptionMap> = HashMap::new();

        {
            for (reference_uniquename, reference_details) in &self.references {
                for feat_annotations in reference_details.cv_annotations.values() {
                    for feat_annotation in feat_annotations.iter() {
                        self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                              &feat_annotation.term);

                       for annotation_detail_id in &feat_annotation.annotations {
                            let annotation_detail = self.annotation_details
                                .get(&annotation_detail_id).expect("can't find OntAnnotationDetail");
                            for gene_uniquename in &annotation_detail.genes {
                                self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                                      gene_uniquename)
                            }
                            for condition_termid in &annotation_detail.conditions {
                                self.add_term_to_hash(&mut seen_terms, reference_uniquename,
                                                      condition_termid);
                            }
                            for ext_part in &annotation_detail.extension {
                                match ext_part.ext_range {
                                    ExtRange::Term(ref range_termid) |
                                    ExtRange::GeneProduct(ref range_termid) =>
                                        self.add_term_to_hash(&mut seen_terms,
                                                              reference_uniquename,
                                                              range_termid),
                                    ExtRange::Gene(ref gene_uniquename) |
                                    ExtRange::Promoter(ref gene_uniquename) =>
                                        self.add_gene_to_hash(&mut seen_genes,
                                                              reference_uniquename,
                                                              gene_uniquename),
                                    _ => {},
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
                    reference_details.physical_interactions.iter()
                    .chain(&reference_details.genetic_interactions);
                for interaction in interaction_iter {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &interaction.interactor_uniquename);
                }

                for ortholog_annotation in &reference_details.ortholog_annotations {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &ortholog_annotation.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &ortholog_annotation.ortholog_uniquename);
                }
                for paralog_annotation in &reference_details.paralog_annotations {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
                                          &paralog_annotation.gene_uniquename);
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename,
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
        }
    }

    pub fn set_counts(&mut self) {
        let mut term_seen_genes: HashMap<TermId, HashSet<GeneUniquename>> = HashMap::new();
        let mut term_seen_genotypes: HashMap<TermId, HashSet<GenotypeUniquename>> = HashMap::new();
        let mut term_seen_single_allele_genotypes: HashMap<TermId, HashSet<GenotypeUniquename>> = HashMap::new();
        let mut ref_seen_genes: HashMap<ReferenceUniquename, HashSet<GeneUniquename>> = HashMap::new();

        for (termid, term_details) in &self.terms {
            let mut seen_genes: HashSet<GeneUniquename> = HashSet::new();
            let mut seen_genotypes: HashSet<GenotypeUniquename> = HashSet::new();
            let mut seen_single_allele_genotypes: HashSet<GenotypeUniquename> = HashSet::new();
            for term_annotations in term_details.cv_annotations.values() {
                for term_annotation in term_annotations {
                    for annotation_detail_id in &term_annotation.annotations {
                        let annotation_detail = self.annotation_details
                            .get(&annotation_detail_id).expect("can't find OntAnnotationDetail");
                        for gene_uniquename in &annotation_detail.genes {
                            seen_genes.insert(gene_uniquename.clone());
                        }
                        if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                            seen_genotypes.insert(genotype_uniquename.clone());
                            let genotype = &self.genotypes[genotype_uniquename];
                            if genotype.expressed_alleles.len() == 1 {
                                seen_single_allele_genotypes.insert(genotype_uniquename.clone());
                            }
                        }
                    }
                }
            }
            term_seen_genes.insert(termid.clone(), seen_genes);
            term_seen_genotypes.insert(termid.clone(), seen_genotypes);
            term_seen_single_allele_genotypes.insert(termid.clone(), seen_single_allele_genotypes);
        }

        for (reference_uniquename, reference_details) in &self.references {
            let mut seen_genes: HashSet<GeneUniquename> = HashSet::new();
            for rel_annotations in reference_details.cv_annotations.values() {
                for rel_annotation in rel_annotations {
                    for annotation_detail_id in &rel_annotation.annotations {
                        let annotation_detail = self.annotation_details
                            .get(&annotation_detail_id).expect("can't find OntAnnotationDetail");
                        if !rel_annotation.is_not {
                            for gene_uniquename in &annotation_detail.genes {
                                seen_genes.insert(gene_uniquename.clone());
                            }
                        }
                    }
                }
            }
            let interaction_iter =
                reference_details.physical_interactions.iter().chain(&reference_details.genetic_interactions);
            for interaction in interaction_iter {
                seen_genes.insert(interaction.gene_uniquename.clone());
                seen_genes.insert(interaction.interactor_uniquename.clone());
            }
            for ortholog_annotation in &reference_details.ortholog_annotations {
                seen_genes.insert(ortholog_annotation.gene_uniquename.clone());
            }
            ref_seen_genes.insert(reference_uniquename.clone(), seen_genes);
        }

        for term_details in self.terms.values_mut() {
            term_details.single_allele_genotype_uniquenames =
                term_seen_single_allele_genotypes.remove(&term_details.termid).unwrap();

            term_details.gene_count =
                term_seen_genes[&term_details.termid].len();
            term_details.genotype_count =
                term_seen_genotypes[&term_details.termid].len();
        }
    }

    // make gene subsets for genes the are not in a slim category
    fn make_non_slim_subset(&self, cv_name: &str, slim_subset: &TermSubsetDetails)
                            -> IdGeneSubsetMap
    {
        let slim_termid_set: HashSet<RcString> =
            slim_subset.elements.keys().map(|termid| termid.clone()).collect();

        let mut non_slim_with_bp_annotation = HashSet::new();
        let mut non_slim_without_bp_annotation = HashSet::new();

        let has_parent_in_slim = |term_annotations: &Vec<OntTermAnnotations>| {
            for term_annotation in term_annotations {
                let interesting_parents =
                    &self.terms[&term_annotation.term].interesting_parents;
                if !term_annotation.is_not &&
                    (slim_termid_set.contains(&term_annotation.term) ||
                     interesting_parents.intersection(&slim_termid_set).count() > 0)
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

            if gene_details.characterisation_status == Some(RcString::from("transposon")) ||
                gene_details.characterisation_status == Some(RcString::from("dubious"))
            {
                continue;
            }

            let mut bp_count = 0;

            if let Some(annotations) =
                gene_details.cv_annotations.get(cv_name) {
                    if has_parent_in_slim(&annotations) {
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
        let name = RcString::from(&format!("non_slim_with_{}_annotation", cv_name));
        return_map.insert(name.clone(),
                          GeneSubsetDetails {
                              name,
                              display_name: RcString::from(&with_annotation_display_name),
                              elements: non_slim_with_bp_annotation,
                          });
        let without_annotation_display_name =
            String::from("Gene products with no ") + &cv_display_name +
            " annotation and are not in a slim category";
        let name = RcString::from(&format!("non_slim_without_{}_annotation", cv_name));
        return_map.insert(name.clone(),
                          GeneSubsetDetails {
                              name,
                              display_name: RcString::from(&without_annotation_display_name),
                              elements: non_slim_without_bp_annotation,
                          });
        return_map
    }

    fn make_slim_subset(&self, slim_name: &str) -> TermSubsetDetails {
        let mut all_genes = HashSet::new();
        let mut slim_subset: HashMap<TermId, TermSubsetElement> = HashMap::new();
        let slim_config = self.config.slims.get(slim_name)
            .expect(&format!("no slim config for {}", slim_name));
        for slim_conf in slim_config.terms.clone() {
            let slim_termid = slim_conf.termid;
            let term_details = self.terms.get(&slim_termid)
                .unwrap_or_else(|| panic!("can't find TermDetails for {}", &slim_termid));

            let subset_element = TermSubsetElement {
                name: term_details.name.clone(),
                gene_count: term_details.genes_annotated_with.len(),
            };

            for gene in &term_details.genes_annotated_with {
                all_genes.insert(gene);
            }
            slim_subset.insert(slim_termid.clone(), subset_element);
        }

        TermSubsetDetails {
            name: RcString::from(slim_name),
            total_gene_count: all_genes.len(),
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
                RcString::from("feature_type:") + &gene_details.feature_type;
            let re = Regex::new(r"[\s,:]+").unwrap();
            let subset_name_no_spaces = RcString::from(&re.replace_all(&subset_name, "_"));
            subsets.entry(subset_name_no_spaces.clone())
                .or_insert(GeneSubsetDetails {
                    name: subset_name_no_spaces,
                    display_name: RcString::from(&subset_name),
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
                    RcString::from("characterisation_status:") + &characterisation_status;
                let re = Regex::new(r"[\s,:]+").unwrap();
                let subset_name_no_spaces = RcString::from(&re.replace_all(&subset_name, "_"));
                subsets.entry(subset_name_no_spaces.clone())
                    .or_insert(GeneSubsetDetails {
                        name: subset_name_no_spaces,
                        display_name: RcString::from(&subset_name),
                        elements: HashSet::new()
                    })
                    .elements.insert(gene_details.uniquename.clone());
            }
        }
    }

    // make InterPro subsets using the interpro_matches field of GeneDetails
    fn make_interpro_subsets(&mut self, subsets: &mut IdGeneSubsetMap) {
        for (gene_uniquename, gene_details) in &self.genes {
            for interpro_match in &gene_details.interpro_matches {

                let mut new_subset_names = vec![];

                if !interpro_match.interpro_id.is_empty() {
                    let subset_name =
                        String::from("interpro:") + &interpro_match.interpro_id;
                    new_subset_names.push((RcString::from(&subset_name),
                                           interpro_match.interpro_name.clone()));
                }

                let subset_name = String::from("interpro:") +
                     &interpro_match.dbname.clone() + ":" + &interpro_match.id;
                new_subset_names.push((RcString::from(&subset_name), interpro_match.name.clone()));

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
            let slim_subset = self.make_slim_subset(&slim_name);
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

    // remove some of the refs that have no annotations.
    // See: https://github.com/pombase/website/issues/628
    fn remove_non_curatable_refs(&mut self) {
        let filtered_refs = self.references.drain()
            .filter(|&(_, ref reference_details)| {
                if reference_has_annotation(reference_details) {
                    return true;
                }
                if let Some(ref canto_triage_status) = reference_details.canto_triage_status {
                    if canto_triage_status == "New" {
                        return false;
                    }
                } else {
                    if reference_details.uniquename.starts_with("PMID:") {
                        print!("reference {} has no canto_triage_status\n", reference_details.uniquename);
                    }
                }
                if let Some (ref triage_status) = reference_details.canto_triage_status {
                    return triage_status != "Wrong organism" && triage_status != "Loaded in error";
                }
                // default to true because there are references that
                // haven't or shouldn't be triaged, eg. GO_REF:...
                true
            })
            .into_iter().collect();

        self.references = filtered_refs;
    }

    fn make_solr_term_summaries(&mut self) -> Vec<SolrTermSummary> {
        let mut return_summaries = vec![];

        let term_name_split_re = Regex::new(r"\W+").unwrap();

        for (termid, term_details) in &self.terms {
            if term_details.cv_annotations.keys()
                .filter(|cv_name| !cv_name.starts_with("extension:"))
                .next().is_none() {
                continue;
            }

            let trimmable_p = |c: char| {
                c.is_whitespace() || c == ',' || c == ':'
                    || c == ';' || c == '.' || c == '\''
            };

            let term_name_words =
                term_name_split_re.split(&term_details.name)
                .map(|s: &str| {
                    s.trim_matches(&trimmable_p).to_owned()
                }).collect::<Vec<String>>();

            let mut close_synonyms = vec![];
            let mut close_synonym_words_vec: Vec<RcString> = vec![];
            let mut distant_synonyms = vec![];
            let mut distant_synonym_words_vec: Vec<RcString> = vec![];

            let add_to_words_vec = |synonym: &str, words_vec: &mut Vec<RcString>| {
                let synonym_words = term_name_split_re.split(&synonym);
                for word in synonym_words {
                    let word_string = RcString::from(word.trim_matches(&trimmable_p));
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
                    !close_synonyms.contains(&synonym)
                })
                .collect::<Vec<_>>();

            let interesting_parents_for_solr =
                term_details.interesting_parents.clone();
            let term_summ = SolrTermSummary {
                id: termid.clone(),
                cv_name: term_details.cv_name.clone(),
                name: term_details.name.clone(),
                definition: term_details.definition.clone(),
                close_synonyms,
                close_synonym_words: RcString::from(&close_synonym_words_vec.join(" ")),
                distant_synonyms,
                distant_synonym_words: RcString::from(&distant_synonym_words_vec.join(" ")),
                interesting_parents: interesting_parents_for_solr,
                secondary_identifiers: term_details.secondary_identifiers.clone(),
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
        self.store_ont_annotations(false);
        self.store_ont_annotations(true);
        self.process_cvtermpath();
        self.process_annotation_feature_rels();
        self.add_target_of_annotations();
        self.set_deletion_viability();
        self.set_term_details_subsets();
        self.set_taxonomic_distributions();
        self.make_residue_extensions();
        self.make_all_cv_summaries();
        self.remove_non_curatable_refs();
        self.set_term_details_maps();
        self.set_gene_details_maps();
        self.set_gene_details_subset_termids();
        self.set_genotype_details_maps();
        self.set_reference_details_maps();
        self.set_counts();
        self.make_subsets();
        self.sort_chromosome_genes();

        let stats = self.get_stats();

        let metadata = self.make_metadata();

        let mut gene_summaries: Vec<GeneSummary> = vec![];
        let mut solr_gene_summaries: Vec<SolrGeneSummary> = vec![];

        for (gene_uniquename, gene_details) in &self.genes {
            if self.config.load_organism_taxonid.is_none() ||
                self.config.load_organism_taxonid.unwrap() == gene_details.taxonid {
                    let gene_summary = self.make_gene_summary(&gene_uniquename);
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

        let solr_term_summaries = self.make_solr_term_summaries();
        let solr_reference_summaries = self.make_solr_reference_summaries();

        let solr_data = SolrData {
            term_summaries: solr_term_summaries,
            gene_summaries: solr_gene_summaries,
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

