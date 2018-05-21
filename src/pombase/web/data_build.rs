use std::rc::Rc;
use std::collections::{BTreeMap, HashMap};
use std::collections::HashSet;
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::u32;

use regex::Regex;

use db::*;
use types::*;
use web::data::*;
use web::config::*;
use web::cv_summary::make_cv_summaries;
use web::util::cmp_str_dates;

use interpro::UniprotResult;

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
    }
}

#[derive(Clone)]
pub struct AlleleAndExpression {
    allele_uniquename: String,
    expression: Option<String>,
}

type TermShortOptionMap = HashMap<TermId, Option<TermShort>>;

type UniprotIdentifier = String;

pub struct WebDataBuild<'a> {
    raw: &'a Raw,
    domain_data: &'a HashMap<UniprotIdentifier, UniprotResult>,
    config: &'a Config,

    genes: UniquenameGeneMap,
    genotypes: UniquenameGenotypeMap,
    alleles: UniquenameAlleleMap,
    other_features: UniquenameFeatureShortMap,
    terms: TermIdDetailsMap,
    chromosomes: ChrNameDetailsMap,
    references: UniquenameReferenceMap,
    all_ont_annotations: HashMap<TermId, Vec<OntAnnotationId>>,
    all_not_ont_annotations: HashMap<TermId, Vec<OntAnnotationId>>,

    genes_of_transcripts: HashMap<String, String>,
    transcripts_of_polypeptides: HashMap<String, String>,
    parts_of_transcripts: HashMap<String, Vec<FeatureShort>>,
    genes_of_alleles: HashMap<String, String>,
    alleles_of_genotypes: HashMap<String, Vec<AlleleAndExpression>>,

    // a map from IDs of terms from the "PomBase annotation extension terms" cv
    // to a Vec of the details of each of the extension
    parts_of_extensions: HashMap<TermId, Vec<ExtPart>>,

    base_term_of_extensions: HashMap<TermId, TermId>,

    children_by_termid: HashMap<TermId, HashSet<TermId>>,
    dbxrefs_of_features: HashMap<String, HashSet<String>>,

    possible_interesting_parents: HashSet<InterestingParent>,

    recent_references: RecentReferences,
    all_community_curated: Vec<ReferenceShort>,

    term_subsets: IdTermSubsetMap,
    gene_subsets: IdGeneSubsetMap,

    annotation_details: IdOntAnnotationDetailMap,
}

fn get_maps() ->
    (HashMap<String, ReferenceShortOptionMap>,
     HashMap<String, GeneShortOptionMap>,
     HashMap<String, GenotypeShortMap>,
     HashMap<String, AlleleShortMap>,
     HashMap<GeneUniquename, TermShortOptionMap>)
{
    (HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new(), HashMap::new())
}

fn get_feat_rel_expression(feature: &Feature,
                           feature_relationship: &FeatureRelationship) -> Option<String> {
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
                         genes: &UniquenameGeneMap, terms: &TermIdDetailsMap) -> String {
    match *ext_range {
        ExtRange::Gene(ref gene_uniquename) => {
            let gene = genes.get(gene_uniquename)
                .unwrap_or_else(|| panic!("can't find gene: {}", gene_uniquename));
            gene_display_name(gene)
        },
        ExtRange::Promoter(ref promoter_uniquename) => promoter_uniquename.clone(),
        ExtRange::SummaryGenes(_) => panic!("can't handle SummaryGenes\n"),
        ExtRange::Term(ref termid) => terms.get(termid).unwrap().name.clone(),
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

fn cmp_genotypes(genotype1: &GenotypeDetails, genotype2: &GenotypeDetails,
                 alleles: &UniquenameAlleleMap) -> Ordering {
    let name1 = genotype_display_name(genotype1, alleles);
    let name2 = genotype_display_name(genotype2, alleles);
    name1.to_lowercase().cmp(&name2.to_lowercase())
}


fn allele_display_name(allele: &AlleleShort) -> String {
    let name = allele.name.clone().unwrap_or_else(|| "unnamed".into());
    let allele_type = allele.allele_type.clone();
    let description = allele.description.clone().unwrap_or_else(|| allele_type.clone());

    if allele_type == "deletion" && name.ends_with("delta") ||
        allele_type.starts_with("wild_type") && name.ends_with('+') {
            let normalised_description = description.replace("[\\s_]+", "");
            let normalised_allele_type = allele_type.replace("[\\s_]+", "");
            if normalised_description != normalised_allele_type {
                return name + "(" + description.as_str() + ")";
            } else {
                return name;
            }
        }

    name + "(" + description.as_str() + ")"
}


fn gene_display_name(gene: &GeneDetails) -> String {
    if let Some(name) = gene.name.clone() {
        name
    } else {
        gene.uniquename.clone()
    }
}

pub fn genotype_display_name(genotype: &GenotypeDetails,
                             alleles: &UniquenameAlleleMap) -> String {
    if let Some(ref name) = genotype.name {
        name.clone()
    } else {
        let allele_display_names: Vec<String> =
            genotype.expressed_alleles.iter().map(|expressed_allele| {
                let allele_short = alleles.get(&expressed_allele.allele_uniquename).unwrap();
                allele_display_name(allele_short)
            }).collect();

        allele_display_names.join(" ")
    }
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
                    (feature_loc.fmin + 1) as u32
                } else {
                    panic!("start_pos less than 1");
                };
            let end_pos =
                if feature_loc.fmax >= 1 {
                    feature_loc.fmax as u32
                } else {
                    panic!("start_end less than 1");
                };
            let feature_uniquename = &feature_loc.srcfeature.uniquename;
            let chr_short = make_chromosome_short(chromosome_map, feature_uniquename);
            Some(ChromosomeLocation {
                chromosome_name: chr_short.name,
                start_pos: start_pos,
                end_pos: end_pos,
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

fn complement_char(base: char) -> char {
    match base {
        'a' => 't',
        'A' => 'T',
        't' => 'a',
        'T' => 'A',
        'g' => 'c',
        'G' => 'C',
        'c' => 'g',
        'C' => 'G',
        _ => 'n',
    }
}

fn rev_comp(residues: &str) -> Residues {
    residues.chars()
        .rev().map(complement_char)
        .collect()
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
                _ => panic!("can't handle feature type: {}", feat.feat_type.name),
            };
            FeatureShort {
                feature_type: feature_type,
                uniquename: feat.uniquename.clone(),
                location: loc,
                residues: residues,
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
            .expect(&format!("missing reference in make_reference_short(): {}",
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
        gene_vec1.iter().map(|gene_uniquename: &String| {
            make_gene_short(genes, gene_uniquename)
        }).collect();
    let gene_short_vec2: Vec<GeneShort> =
        gene_vec2.iter().map(|gene_uniquename: &String| {
            make_gene_short(genes, gene_uniquename)
        }).collect();

    gene_short_vec1.cmp(&gene_short_vec2)
}

lazy_static
! {
    static ref MODIFICATION_RE: Regex = Regex::new(r"^(?P<aa>[A-Z])(?P<pos>\d+)$").unwrap();
}

fn cmp_residues(residue1: &Option<Residue>, residue2: &Option<Residue>) -> Ordering {
    if let Some(ref res1) = *residue1 {
        if let Some(ref res2) = *residue2 {
            if let (Some(res1_captures), Some(res2_captures)) =
                (MODIFICATION_RE.captures(res1), MODIFICATION_RE.captures(res2))
            {
                let res1_aa = res1_captures.name("aa").unwrap().as_str();
                let res2_aa = res2_captures.name("aa").unwrap().as_str();
                let aa_order = res1_aa.cmp(&res2_aa);
                if aa_order == Ordering::Equal {
                    let res1_pos =
                        res1_captures.name("pos").unwrap().as_str().parse::<i32>().unwrap();
                    let res2_pos =
                        res2_captures.name("pos").unwrap().as_str().parse::<i32>().unwrap();
                    res1_pos.cmp(&res2_pos)
                } else {
                    aa_order
                }
            } else {
                res1.cmp(&res2)
            }
        } else {
            Ordering::Less
        }
    } else {
        if let Some(_) = *residue2 {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

pub fn cmp_ont_annotation_detail(cv_config: &CvConfig,
                             detail1: &OntAnnotationDetail,
                             detail2: &OntAnnotationDetail,
                             genes: &UniquenameGeneMap,
                             genotypes: &UniquenameGenotypeMap,
                             alleles: &UniquenameAlleleMap,
                             terms: &TermIdDetailsMap) -> Result<Ordering, String> {
    if let Some(ref detail1_genotype_uniquename) = detail1.genotype {
        if let Some(ref detail2_genotype_uniquename) = detail2.genotype {
            let genotype1 = genotypes.get(detail1_genotype_uniquename).unwrap();
            let genotype2 = genotypes.get(detail2_genotype_uniquename).unwrap();

            let ord = cmp_genotypes(genotype1, genotype2, alleles);

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
        if let Some(ref conf_termid) = ext_conf.if_descendent_of {
            ret.insert(InterestingParent {
                termid: conf_termid.clone(),
                rel_name: "is_a".into(),
            });
        }
    }

    for go_slim_conf in &config.go_slim_terms {
        for rel_name in DESCENDANT_REL_NAMES.iter() {
            ret.insert(InterestingParent {
                termid: go_slim_conf.termid.clone(),
                rel_name: (*rel_name).to_owned(),
            });
        }
    }

    ret.insert(InterestingParent {
        termid: config.viability_terms.viable.clone(),
        rel_name: "is_a".into(),
    });
    ret.insert(InterestingParent {
        termid: config.viability_terms.inviable.clone(),
        rel_name: "is_a".into(),
    });

    for (cv_name, conf) in &config.cv_config {
        for filter in &conf.filters {
            for category in &filter.term_categories {
                for ancestor in &category.ancestors {
                    for config_rel_name in &DESCENDANT_REL_NAMES {
                        if *config_rel_name == "has_part" &&
                            !HAS_PART_CV_NAMES.contains(&cv_name.as_str()) {
                                continue;
                            }
                        ret.insert(InterestingParent {
                            termid: ancestor.clone(),
                            rel_name: String::from(*config_rel_name),
                        });
                    }
                }
            }
        }

        for split_by_parent_config in &conf.split_by_parents {
            for ancestor in &split_by_parent_config.termids {
                let ancestor_termid =
                    if ancestor.starts_with("NOT ") {
                        ancestor[4..].to_owned()
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
                       all_ref_uniquenames: &[String]) -> Vec<ReferenceShort> {
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
                      all_ref_uniquenames: &[String])
                      -> (Vec<ReferenceShort>, Vec<ReferenceShort>, Vec<ReferenceShort>) {
    let mut sorted_pub_uniquenames: Vec<ReferenceUniquename> =
        all_ref_uniquenames.iter()
        .filter(|ref_uniquename| {
            let reference = references_map.get(*ref_uniquename).unwrap();
            (reference.canto_first_approved_date.is_some() ||
             reference.canto_session_submitted_date.is_some()) &&
                reference.canto_curator_role.is_some()
        })
        .cloned()
        .collect();

    {
        let submitted_date_cmp =
            |ref_uniquename1: &ReferenceUniquename, ref_uniquename2: &ReferenceUniquename| {
                let ref1 = references_map.get(ref_uniquename1).unwrap();
                let ref2 = references_map.get(ref_uniquename2).unwrap();

                // use first approval date, but fall back to the most recent approval date

                if let Some(ref ref1_date) = ref1.canto_first_approved_date {
                    if let Some(ref ref2_date) = ref2.canto_first_approved_date {
                        return cmp_str_dates(ref1_date, ref2_date).reverse();
                    }
                }

                if let Some(ref ref1_date) = ref1.canto_session_submitted_date {
                    if let Some(ref ref2_date) = ref2.canto_first_approved_date {
                        return cmp_str_dates(ref1_date, ref2_date).reverse();
                    }
                }

                if let Some(ref ref1_date) = ref1.canto_first_approved_date {
                    if let Some(ref ref2_date) = ref2.canto_session_submitted_date {
                        return cmp_str_dates(ref1_date, ref2_date).reverse();
                    }
                }

                if let Some(ref ref1_date) = ref1.canto_session_submitted_date {
                    if let Some(ref ref2_date) = ref2.canto_session_submitted_date {
                        return cmp_str_dates(ref1_date, ref2_date).reverse();
                    }
                }

                panic!();
            };

        sorted_pub_uniquenames.sort_by(submitted_date_cmp);
    }

    let mut recent_admin_curated = vec![];
    let mut recent_community_curated = vec![];
    let mut all_community_curated = vec![];

    let ref_uniquename_iter = sorted_pub_uniquenames.iter();

    for ref_uniquename in ref_uniquename_iter {
        let reference = references_map.get(ref_uniquename).unwrap();

        if reference.canto_curator_role == Some("community".into()) {
            let ref_short = make_reference_short(references_map, ref_uniquename).unwrap();
            all_community_curated.push(ref_short.clone());
            if recent_community_curated.len() <= MAX_RECENT_REFS {
                recent_community_curated.push(ref_short);
            }
        } else {
            if recent_admin_curated.len() <= MAX_RECENT_REFS {
                let ref_short = make_reference_short(references_map, ref_uniquename).unwrap();
                recent_admin_curated.push(ref_short);
            }
        }
    }

    (recent_admin_curated, recent_community_curated, all_community_curated)
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
                    uniquename: intron_uniquename,
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
            raw: raw,
            domain_data: domain_data,
            config: config,

            genes: BTreeMap::new(),
            genotypes: HashMap::new(),
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
       }
    }

    fn add_ref_to_hash(&self,
                       seen_references: &mut HashMap<String, ReferenceShortOptionMap>,
                       identifier: String,
                       maybe_reference_uniquename: Option<ReferenceUniquename>) {
        if let Some(reference_uniquename) = maybe_reference_uniquename {
            seen_references
                .entry(identifier.clone())
                .or_insert_with(HashMap::new)
                .insert(reference_uniquename.clone(), None);
        }
    }


    fn add_gene_to_hash(&self,
                        seen_genes: &mut HashMap<String, GeneShortOptionMap>,
                        identifier: String,
                        other_gene_uniquename: GeneUniquename) {
        seen_genes
            .entry(identifier)
            .or_insert_with(HashMap::new)
            .insert(other_gene_uniquename.clone(), None);
    }

    fn add_genotype_to_hash(&self,
                            seen_genotypes: &mut HashMap<String, GenotypeShortMap>,
                            seen_alleles: &mut HashMap<String, AlleleShortMap>,
                            seen_genes: &mut HashMap<String, GeneShortOptionMap>,
                            identifier: String,
                            genotype_uniquename: &str) {
        let genotype = self.make_genotype_short(genotype_uniquename);
        for expressed_allele in &genotype.expressed_alleles {
            self.add_allele_to_hash(seen_alleles, seen_genes, identifier.clone(),
                                    expressed_allele.allele_uniquename.clone());
        }

        seen_genotypes
            .entry(identifier)
            .or_insert_with(HashMap::new)
            .insert(genotype_uniquename.to_owned(),
                    self.make_genotype_short(genotype_uniquename));
    }

    fn add_allele_to_hash(&self,
                          seen_alleles: &mut HashMap<String, AlleleShortMap>,
                          seen_genes: &mut HashMap<String, GeneShortOptionMap>,
                          identifier: String,
                          allele_uniquename: AlleleUniquename) -> AlleleShort {
        let allele_short = self.make_allele_short(&allele_uniquename);
        let allele_gene_uniquename =
            allele_short.gene_uniquename.clone();
        self.add_gene_to_hash(seen_genes, identifier.clone(), allele_gene_uniquename);
        seen_alleles
            .entry(identifier)
            .or_insert_with(HashMap::new)
            .insert(allele_uniquename, allele_short.clone());
        allele_short
    }

    fn add_term_to_hash(&self,
                        seen_terms: &mut HashMap<TermId, TermShortOptionMap>,
                        identifier: String,
                        other_termid: TermId) {
        seen_terms
            .entry(identifier)
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
            .collect::<Vec<String>>();
        let mut ortholog_ids =
            gene_details.ortholog_annotations.iter()
            .map(|ortholog_annotation| {
                IdAndOrganism {
                    identifier: ortholog_annotation.ortholog_uniquename.clone(),
                    taxonid: ortholog_annotation.ortholog_taxonid,
                }
            })
            .collect::<Vec<IdAndOrganism>>();

        for ortholog_annotation in &gene_details.ortholog_annotations {
            let orth_uniquename = &ortholog_annotation.ortholog_uniquename;
            if let Some(orth_gene) =
                self.genes.get(orth_uniquename) {
                    if let Some(ref orth_name) = orth_gene.name {
                        let id_and_org =IdAndOrganism {
                            identifier: orth_name.clone(),
                            taxonid: ortholog_annotation.ortholog_taxonid,
                        };
                        ortholog_ids.push(id_and_org);
                    }
                } else {
                    panic!("missing GeneShort for: {:?}", orth_uniquename);
                }
        }
        GeneSummary {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
            uniprot_identifier: gene_details.uniprot_identifier.clone(),
            synonyms: synonyms,
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
            .collect::<Vec<String>>();
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
            exon_count: exon_count,
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
                                   cvterm_name: &String) {
        let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
        gene_details.characterisation_status = Some(cvterm_name.clone());
    }

    fn add_gene_product(&mut self, gene_uniquename: &str, product: &str) {
        let gene_details = self.get_gene_mut(gene_uniquename);
        gene_details.product = Some(product.to_owned());
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

            let mut pubmed_authors: Option<String> = None;
            let mut pubmed_publication_date: Option<String> = None;
            let mut pubmed_abstract: Option<String> = None;
            let mut canto_triage_status: Option<String> = None;
            let mut canto_curator_role: Option<String> = None;
            let mut canto_curator_name: Option<String> = None;
            let mut canto_first_approved_date: Option<String> = None;
            let mut canto_approved_date: Option<String> = None;
            let mut canto_added_date: Option<String> = None;
            let mut canto_session_submitted_date: Option<String> = None;

            for prop in rc_publication.publicationprops.borrow().iter() {
                match &prop.prop_type.name as &str {
                    "pubmed_publication_date" =>
                        pubmed_publication_date = Some(prop.value.clone()),
                    "pubmed_authors" =>
                        pubmed_authors = Some(prop.value.clone()),
                    "pubmed_abstract" =>
                        pubmed_abstract = Some(prop.value.clone()),
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

            let mut authors_abbrev = None;
            let mut publication_year = None;

            if let Some(authors) = pubmed_authors.clone() {
                if authors.contains(',') {
                    let author_re = Regex::new(r"^(?P<f>[^,]+),.*$").unwrap();
                    let replaced: String =
                        author_re.replace_all(&authors, "$f et al.").into();
                    authors_abbrev = Some(replaced);
                } else {
                    authors_abbrev = Some(authors.clone());
                }
            }

            if let Some(publication_date) = pubmed_publication_date.clone() {
                let date_re = Regex::new(r"^(.* )?(?P<y>\d\d\d\d)$").unwrap();
                publication_year = Some(date_re.replace_all(&publication_date, "$y").into());
            }

            let mut approved_date = canto_first_approved_date.clone();

            if approved_date.is_none() {
                approved_date = canto_session_submitted_date.clone();
            }
            approved_date =
                if let Some(date) = approved_date {
                    let re = Regex::new(r"^(?P<date>\d\d\d\d-\d\d-\d\d).*").unwrap();
                    Some(re.replace_all(&date, "$date").into())
                } else {
                    None
                };

            self.references.insert(reference_uniquename.clone(),
                                   ReferenceDetails {
                                       uniquename: reference_uniquename.clone(),
                                       title: rc_publication.title.clone(),
                                       citation: rc_publication.miniref.clone(),
                                       pubmed_abstract: pubmed_abstract.clone(),
                                       authors: pubmed_authors.clone(),
                                       authors_abbrev: authors_abbrev,
                                       pubmed_publication_date: pubmed_publication_date.clone(),
                                       canto_triage_status: canto_triage_status,
                                       canto_curator_role: canto_curator_role,
                                       canto_curator_name: canto_curator_name,
                                       canto_first_approved_date: canto_first_approved_date,
                                       canto_approved_date: canto_approved_date,
                                       canto_session_submitted_date: canto_session_submitted_date,
                                       canto_added_date: canto_added_date,
                                       approved_date: approved_date,
                                       publication_year: publication_year,
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
             all_community_curated) =
            make_canto_curated(&self.references, &all_uniquenames);

        let recent_references = RecentReferences {
            pubmed: make_recently_added(&self.references, &all_uniquenames),
            admin_curated: recent_admin_curated,
            community_curated: recent_community_curated,
        };

        self.recent_references = recent_references;

        self.all_community_curated = all_community_curated;
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
                            AlleleAndExpression {
                                allele_uniquename: subject_uniquename.clone(),
                                expression: expression,
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

    fn get_feature_dbxrefs(&self, feature: &Feature) -> HashSet<String> {
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
                orfeome_identifier = Some(String::from(&dbxref[4..]));
            }
        }

        let mut uniprot_identifier = None;
        for prop in feat.featureprops.borrow().iter() {
            if prop.prop_type.name == "uniprot_identifier" {
                uniprot_identifier = prop.value.clone();
                break;
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

        let gene_feature = GeneDetails {
            uniquename: feat.uniquename.clone(),
            name: feat.name.clone(),
            taxonid: organism.taxonid,
            product: None,
            deletion_viability: DeletionViability::Unknown,
            uniprot_identifier: uniprot_identifier,
            interpro_matches: interpro_matches,
            tm_domain_coords: tm_domain_coords,
            orfeome_identifier: orfeome_identifier,
            name_descriptions: vec![],
            synonyms: vec![],
            dbxrefs: dbxrefs,
            feature_type: feat.feat_type.name.clone(),
            characterisation_status: None,
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
        };

        self.genes.insert(feat.uniquename.clone(), gene_feature);
    }

    fn get_transcript_parts(&mut self, transcript_uniquename: &str) -> Vec<FeatureShort> {
        let mut parts = self.parts_of_transcripts.remove(transcript_uniquename)
            .expect("can't find transcript");

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
    }

    fn store_transcript_details(&mut self, feat: &Feature) {
        let transcript_uniquename = feat.uniquename.clone();

        let parts = self.get_transcript_parts(&transcript_uniquename);

        if parts.len() == 0 {
            panic!("transcript has no parts");
        }

        let mut transcript_start = u32::MAX;
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
                let mut cds_start = u32::MAX;
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
            parts: parts,
            protein: None,
            cds_location: maybe_cds_location,
        };

        if let Some(gene_uniquename) =
            self.genes_of_transcripts.get(&transcript_uniquename) {
                let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
                if gene_details.feature_type != "pseudogene" {
                    gene_details.feature_type =
                        transcript.transcript_type.clone() + " " + &gene_details.feature_type;
                }
                gene_details.transcripts.push(transcript);
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

            let parse_prop_as_f32 = |p: &Option<String>| {
                if let Some(prop_value) = p.clone() {
                    let maybe_value = prop_value.parse();
                    if let Ok(parsed_prop) = maybe_value {
                        Some(parsed_prop)
                    } else {
                        println!("{}: couldn't parse {} as f32",
                                 feat.uniquename, prop_value);
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
                sequence: residues,
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

        let chr = ChromosomeDetails {
            name: feat.uniquename.clone(),
            residues: feat.residues.clone().unwrap(),
            ena_identifier: ena_identifier.unwrap(),
            gene_uniquenames: vec![],
            taxonid: org.taxonid,
        };

        self.chromosomes.insert(feat.uniquename.clone(), chr);
    }

    fn store_genotype_details(&mut self, feat: &Feature) {
        let mut background = None;

        for prop in feat.featureprops.borrow().iter() {
            if prop.prop_type.name == "genotype_background" {
                background = prop.value.clone()
            }
        }

        self.genotypes.insert(feat.uniquename.clone(),
                              GenotypeDetails {
                                  uniquename: feat.uniquename.clone(),
                                  name: feat.name.clone(),
                                  background: background,
                                  expressed_alleles: vec![],
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
            gene_uniquename: gene_uniquename,
            allele_type: allele_type.unwrap(),
            description: description,
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
        let mut interesting_parents_by_termid: HashMap<String, HashSet<String>> =
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
                    .insert(object_termid.into());
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
            gene_uniquename: String,
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

    fn add_alleles_to_genotypes(&mut self) {
        let mut alleles_to_add: HashMap<String, Vec<ExpressedAllele>> = HashMap::new();

        for genotype_uniquename in self.genotypes.keys() {
            let allele_uniquenames: Vec<AlleleAndExpression> =
                self.alleles_of_genotypes[genotype_uniquename].clone();
            let expressed_allele_vec: Vec<ExpressedAllele> =
                allele_uniquenames.iter()
                .map(|allele_and_expression| {
                    ExpressedAllele {
                        allele_uniquename: allele_and_expression.allele_uniquename.clone(),
                        expression: allele_and_expression.expression.clone(),
                    }
                })
                .collect();

            alleles_to_add.insert(genotype_uniquename.clone(), expressed_allele_vec);
        }

        {
            let allele_cmp = |allele1: &ExpressedAllele, allele2: &ExpressedAllele| {
                let allele1_display_name =
                    allele_display_name(&self.alleles.get(&allele1.allele_uniquename).unwrap());
                let allele2_display_name =
                    allele_display_name(&self.alleles.get(&allele2.allele_uniquename).unwrap());
                allele1_display_name.cmp(&allele2_display_name)
            };

            for (_, alleles) in &mut alleles_to_add {
                alleles.sort_by(&allele_cmp);
            }
        }

        for (genotype_uniquename, genotype_details) in &mut self.genotypes {
            genotype_details.expressed_alleles =
                alleles_to_add.remove(genotype_uniquename).unwrap();
        }
    }

    // add interaction, ortholog and paralog annotations
    fn process_annotation_feature_rels(&mut self) {
        for feature_rel in self.raw.feature_relationships.iter() {
            let rel_name = &feature_rel.rel_type.name;
            let subject_uniquename = &feature_rel.subject.uniquename;
            let object_uniquename = &feature_rel.object.uniquename;

            for rel_config in FEATURE_REL_CONFIGS.iter() {
                if rel_name == rel_config.rel_type_name &&
                    is_gene_type(&feature_rel.subject.feat_type.name) &&
                    is_gene_type(&feature_rel.object.feat_type.name) {
                        let mut evidence: Option<Evidence> = None;
                        let mut is_inferred_interaction: bool = false;

                        let borrowed_publications = feature_rel.publications.borrow();
                        let maybe_publication = borrowed_publications.get(0).clone();
                        let maybe_reference_uniquename =
                            match maybe_publication {
                                Some(publication) => Some(publication.uniquename.clone()),
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
                        }

                        let evidence_clone = evidence.clone();

                        let gene_uniquename = subject_uniquename;
                        let gene_organism_taxonid = {
                            self.genes.get(subject_uniquename).unwrap().taxonid.clone()
                        };
                        let other_gene_uniquename = object_uniquename;
                        let other_gene_organism_taxonid = {
                            self.genes.get(object_uniquename).unwrap().taxonid.clone()
                        };
                        match rel_config.annotation_type {
                            FeatureRelAnnotationType::Interaction =>
                                if !is_inferred_interaction {
                                    let interaction_annotation =
                                        InteractionAnnotation {
                                            gene_uniquename: gene_uniquename.clone(),
                                            interactor_uniquename: other_gene_uniquename.clone(),
                                            evidence: evidence,
                                            reference_uniquename: maybe_reference_uniquename.clone(),
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
                                        evidence: evidence,
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
                                    if self.config.load_organism_taxonid == gene_details.taxonid {
                                        ref_details.ortholog_annotations.push(ortholog_annotation);
                                    }
                                }
                            },
                            FeatureRelAnnotationType::Paralog => {
                                let paralog_annotation =
                                    ParalogAnnotation {
                                        gene_uniquename: gene_uniquename.clone(),
                                        paralog_uniquename: other_gene_uniquename.clone(),
                                        evidence: evidence,
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
                                    if self.config.load_organism_taxonid == gene_details.taxonid {
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
                                    if self.config.load_organism_taxonid == other_gene_details.taxonid {
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
                                    if self.config.load_organism_taxonid == other_gene_details.taxonid {
                                        ref_details.paralog_annotations.push(paralog_annotation);
                                    }
                                }
                            },
                        }
                    }
            }
        }

        for (_, ref_details) in &mut self.references {
            ref_details.physical_interactions.sort();
            ref_details.genetic_interactions.sort();
            ref_details.ortholog_annotations.sort();
            ref_details.paralog_annotations.sort();
        }

        for (_, gene_details) in &mut self.genes {
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
                    if let Some(if_descendent_of) = ext_config.if_descendent_of.clone() {
                        if annotation_term_details.interesting_parents.contains(&if_descendent_of) {
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
    fn make_target_of_for_ext(&self, cv_name: &String,
                              genes: &Vec<String>,
                              maybe_genotype_uniquename: &Option<String>,
                              reference_uniquename: &Option<String>,
                              annotation_termid: &String,
                              extension: &Vec<ExtPart>) -> Vec<(GeneUniquename, TargetOfAnnotation)> {
        let mut ret_vec = vec![];

        for ext_part in extension {
            let maybe_ext_config =
                self.matching_ext_config(annotation_termid, &ext_part.rel_type_name);
            if let ExtRange::Gene(ref target_gene_uniquename) = ext_part.ext_range {
                if let Some(ext_config) = maybe_ext_config {
                    if let Some(reciprocal_display_name) =
                        ext_config.reciprocal_display {
                            let (annotation_gene_uniquenames, annotation_genotype_uniquename) =
                                if maybe_genotype_uniquename.is_some() {
                                    (genes.clone(), maybe_genotype_uniquename.clone())
                                } else {
                                    (genes.clone(), None)
                                };
                            ret_vec.push(((*target_gene_uniquename).clone(),
                                          TargetOfAnnotation {
                                              ontology_name: cv_name.clone(),
                                              ext_rel_display_name: reciprocal_display_name,
                                              genes: annotation_gene_uniquenames,
                                              genotype_uniquename: annotation_genotype_uniquename,
                                              reference_uniquename: reference_uniquename.clone(),
                                          }));
                        }
                }
            }
        }

        ret_vec
    }

    fn add_target_of_annotations(&mut self) {
        let mut target_of_annotations: HashMap<GeneUniquename, HashSet<TargetOfAnnotation>> =
            HashMap::new();

        for (_, term_details) in &self.terms {
            for (_, term_annotations) in &term_details.cv_annotations {
                for term_annotation in term_annotations {
                    'ANNOTATION: for annotation_id in &term_annotation.annotations {
                        let annotation = self.annotation_details
                            .get(&annotation_id).expect("can't find OntAnnotationDetail");
                        if let Some(ref genotype_uniquename) = annotation.genotype {
                            let genotype = self.genotypes.get(genotype_uniquename).unwrap();

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
                                .or_insert(HashSet::new())
                                .insert(new_annotation);
                        }
                    }
                }
            }
        }

        for (gene_uniquename, mut target_of_annotations) in target_of_annotations {
            let gene_details = self.genes.get_mut(&gene_uniquename).unwrap();
            gene_details.target_of_annotations = target_of_annotations.drain().collect();
        }
    }

    fn set_deletion_viability(&mut self) {
        let mut gene_statuses = HashMap::new();

        let condition_string =
            |condition_ids: HashSet<String>| {
                let mut ids_vec: Vec<String> = condition_ids.iter().cloned().collect();
                ids_vec.sort();
                ids_vec.join(" ")
            };

        let viable_termid = &self.config.viability_terms.viable;
        let inviable_termid = &self.config.viability_terms.inviable;

        for (gene_uniquename, gene_details) in &mut self.genes {
            let mut new_status = DeletionViability::Unknown;

            if let Some(single_allele_term_annotations) =
                gene_details.cv_annotations.get("single_allele_phenotype") {
                    let mut viable_conditions: HashMap<String, TermId> = HashMap::new();
                    let mut inviable_conditions: HashMap<String, TermId> = HashMap::new();

                    for term_annotation in single_allele_term_annotations {
                        'ANNOTATION: for annotation_id in &term_annotation.annotations {
                            let annotation = self.annotation_details
                                .get(&annotation_id).expect("can't find OntAnnotationDetail");

                            let genotype_uniquename = annotation.genotype.as_ref().unwrap();

                            let genotype = self.genotypes.get(genotype_uniquename).unwrap();
                            let allele_uniquename =
                                genotype.expressed_alleles[0].allele_uniquename.clone();
                            let allele = self.alleles.get(&allele_uniquename).unwrap();
                            if allele.allele_type != "deletion" {
                                continue 'ANNOTATION;
                            }

                            let term = self.terms.get(&term_annotation.term).unwrap();
                            let ref interesting_parents = term.interesting_parents;
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

                    if viable_conditions.len() == 0 {
                        if inviable_conditions.len() > 0 {
                            new_status = DeletionViability::Inviable;
                        }
                    } else {
                        if inviable_conditions.len() == 0 {
                            new_status = DeletionViability::Viable;
                        } else {
                            new_status = DeletionViability::DependsOnConditions;

                            let viable_conditions_set: HashSet<String> =
                                viable_conditions.keys().cloned().collect();
                            let inviable_conditions_set: HashSet<String> =
                                inviable_conditions.keys().cloned().collect();

                            let intersecting_conditions =
                                viable_conditions_set.intersection(&inviable_conditions_set);
                            if intersecting_conditions.clone().count() > 0 {
                                println!("{} is viable and inviable with", gene_uniquename);
                                for cond in intersecting_conditions {
                                    if cond.len() == 0 {
                                        println!("  no conditions");
                                    } else {
                                        println!("  conditions: {}", cond);
                                    }
                                    println!("   viable term: {}",
                                             viable_conditions.get(cond).unwrap());
                                    println!("   inviable term: {}",
                                             inviable_conditions.get(cond).unwrap());
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
        'TERM: for go_slim_conf in self.config.go_slim_terms.clone() {
            let slim_termid = &go_slim_conf.termid;
            for (_, term_details) in &mut self.terms {
                if term_details.termid == *slim_termid {
                    term_details.subsets.push("goslim_pombe".into());
                    break 'TERM;
                }
            }
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

    fn make_all_cv_summaries(&mut self) {
        let gene_short_map = self.make_gene_short_map();

        for (_, term_details) in &mut self.terms {
            make_cv_summaries(term_details, self.config, &self.children_by_termid,
                              true, true, &gene_short_map, &self.annotation_details);
        }

        for (_, gene_details) in &mut self.genes {
            make_cv_summaries(gene_details, &self.config, &self.children_by_termid,
                              false, true, &gene_short_map, &self.annotation_details);
        }

        for (_, genotype_details) in &mut self.genotypes {
            make_cv_summaries(genotype_details, &self.config, &self.children_by_termid,
                              false, false, &gene_short_map, &self.annotation_details);
        }

        for (_, reference_details) in &mut self.references {
            make_cv_summaries( reference_details, &self.config, &self.children_by_termid,
                              true, true, &gene_short_map, &self.annotation_details);
        }
    }

    fn process_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name != POMBASE_ANN_EXT_TERM_CV_NAME {
                let cv_config =
                    self.config.cv_config_by_name(&cvterm.cv.name);
                let annotation_feature_type =
                    cv_config.feature_type.clone();
                let synonyms =
                    cvterm.cvtermsynonyms.borrow().iter().map(|syn| {
                        SynonymDetails {
                            synonym_type: (*syn).synonym_type.name.clone(),
                            name: syn.name.clone(),
                        }
                    }).collect::<Vec<_>>();
                self.terms.insert(cvterm.termid(),
                                  TermDetails {
                                      name: cvterm.name.clone(),
                                      cv_name: cvterm.cv.name.clone(),
                                      annotation_feature_type: annotation_feature_type,
                                      interesting_parents: HashSet::new(),
                                      subsets: vec![],
                                      termid: cvterm.termid(),
                                      synonyms: synonyms,
                                      definition: cvterm.definition.clone(),
                                      direct_ancestors: vec![],
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
                                  });
            }
        }
    }

    fn get_ext_rel_display_name(&self, annotation_termid: &String,
                                ext_rel_name: &String) -> String {
        if let Some(ext_conf) = self.matching_ext_config(annotation_termid, ext_rel_name) {
            ext_conf.display_name.clone()
        } else {
            let re = Regex::new("_").unwrap();
            re.replace_all(&ext_rel_name, " ").into()
        }
    }

    fn process_extension_cvterms(&mut self) {
        for cvterm in &self.raw.cvterms {
            if cvterm.cv.name == POMBASE_ANN_EXT_TERM_CV_NAME {
                for cvtermprop in cvterm.cvtermprops.borrow().iter() {
                    if (*cvtermprop).prop_type.name.starts_with(ANNOTATION_EXT_REL_PREFIX) {
                        let ext_rel_name_str =
                            &(*cvtermprop).prop_type.name[ANNOTATION_EXT_REL_PREFIX.len()..];
                        let ext_rel_name = String::from(ext_rel_name_str);
                        let ext_range = (*cvtermprop).value.clone();
                        let range: ExtRange = if ext_range.starts_with("SP") {
                            if ext_range.ends_with("-promoter") {
                                ExtRange::Promoter(ext_range)
                            } else {
                                ExtRange::Gene(ext_range)
                            }
                        } else {
                            ExtRange::Misc(ext_range)
                        };

                        if let Some(base_termid) =
                            self.base_term_of_extensions.get(&cvterm.termid()) {
                                let rel_type_display_name =
                                    self.get_ext_rel_display_name(&base_termid, &ext_rel_name);

                                self.parts_of_extensions.entry(cvterm.termid())
                                    .or_insert(Vec::new()).push(ExtPart {
                                        rel_type_name: String::from(ext_rel_name),
                                        rel_type_display_name: rel_type_display_name,
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
                    if let Some(base_termid) =
                        self.base_term_of_extensions.get(&subject_term.termid()) {
                            let rel_type_display_name =
                                self.get_ext_rel_display_name(base_termid, &rel_type.name);

                            self.parts_of_extensions.entry(subject_termid)
                                .or_insert(Vec::new()).push(ExtPart {
                                    rel_type_name: rel_type.name.clone(),
                                    rel_type_display_name: rel_type_display_name,
                                    ext_range: ExtRange::Term(object_term.termid().clone()),
                                });
                        } else {
                            panic!("can't find details for {}\n", object_term.termid());
                        }
                }
            }
        }
    }

    fn process_feature_synonyms(&mut self) {
        for feature_synonym in self.raw.feature_synonyms.iter() {
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

    fn make_genotype_short(&self, genotype_uniquename: &str) -> GenotypeShort {
        if let Some(ref details) = self.genotypes.get(genotype_uniquename) {
            GenotypeShort {
                uniquename: details.uniquename.clone(),
                name: details.name.clone(),
                background: details.background.clone(),
                expressed_alleles: details.expressed_alleles.clone(),
            }
        } else {
            panic!("can't find genotype");
        }
    }

    fn make_allele_short(&self, allele_uniquename: &str) -> AlleleShort {
        self.alleles.get(allele_uniquename).unwrap().clone()
    }

    // process feature properties stored as cvterms,
    // eg. characterisation_status and product
    fn process_props_from_feature_cvterms(&mut self) {
        for feature_cvterm in self.raw.feature_cvterms.iter() {
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

    fn get_gene_prod_extension(&self, prod_value: &String) -> ExtPart {
        let ext_range =
            if prod_value.starts_with("PR:") {
                ExtRange::GeneProduct(prod_value.clone())
            } else {
                ExtRange::Misc(prod_value.clone())
            };

        ExtPart {
            rel_type_name: "active_form".into(),
            rel_type_display_name: "active form".into(),
            ext_range: ext_range,
        }
    }

    // return a fake extension for "with" properties on protein binding annotations
    fn get_with_extension(&self, with_value: &String) -> ExtPart {
        let ext_range =
            if with_value.starts_with("SP%") {
                ExtRange::Gene(with_value.clone())
            } else {
                if with_value.starts_with("PomBase:SP") {
                    let gene_uniquename =
                        String::from(&with_value[8..]);
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
            rel_type_name: "binds".into(),
            rel_type_display_name: "binds".into(),
            ext_range: ext_range,
        }
    }

    fn make_with_or_from_value(&self, with_or_from_value: String) -> WithFromValue {
        let db_prefix_patt = String::from("^") + DB_NAME + ":";
        let re = Regex::new(&db_prefix_patt).unwrap();
        let gene_uniquename: String = re.replace_all(&with_or_from_value, "").into();
        if self.genes.contains_key(&gene_uniquename) {
            let gene_short = self.make_gene_short(&gene_uniquename);
            WithFromValue::Gene(gene_short)
        } else {
            if self.terms.get(&with_or_from_value).is_some() {
                WithFromValue::Term(self.make_term_short(&with_or_from_value))
            } else {
                WithFromValue::Identifier(with_or_from_value)
            }
        }
    }

    // add the with value as a fake extension if the cvterm is_a protein binding,
    // otherwise return the value
    fn make_with_extension(&self, termid: &String, evidence_code: Option<String>,
                           extension: &mut Vec<ExtPart>,
                           with_value: String) -> Option<WithFromValue> {
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
                extension.push(self.get_with_extension(&with_value));
            } else {
                return Some(self.make_with_or_from_value(with_value));
            }
        None
    }

    // process annotation
    fn process_feature_cvterms(&mut self) {
        for feature_cvterm in self.raw.feature_cvterms.iter() {
            let feature = &feature_cvterm.feature;
            let cvterm = &feature_cvterm.cvterm;

            let mut extension = vec![];

            if cvterm.cv.name == "PomBase gene characterisation status" ||
                cvterm.cv.name == "PomBase gene products" ||
                cvterm.cv.name == "name_description" {
                    continue;
                }

            let publication = &feature_cvterm.publication;
            let mut extra_props: HashMap<String, String> = HashMap::new();
            let mut conditions: HashSet<TermId> = HashSet::new();
            let mut withs: Vec<WithFromValue> = vec![];
            let mut froms: Vec<WithFromValue> = vec![];
            let mut qualifiers: Vec<Qualifier> = vec![];
            let mut evidence: Option<String> = None;

            // need to get evidence first as it's used later
            // See: https://github.com/pombase/website/issues/455
            for ref prop in feature_cvterm.feature_cvtermprops.borrow().iter() {
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

            for ref prop in feature_cvterm.feature_cvtermprops.borrow().iter() {
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
                            conditions.insert(value.clone());
                        },
                    "qualifier" =>
                        if let Some(value) = prop.value.clone() {
                            qualifiers.push(value);
                        },
                    "with" => {
                        if let Some(with_value) = prop.value.clone() {
                            if let Some(with_gene_short) =
                                self.make_with_extension(&cvterm.termid(), evidence.clone(),
                                                         &mut extension, with_value) {
                                    withs.push(with_gene_short);
                                }
                        }
                    },
                    "from" => {
                        if let Some(value) = prop.value.clone() {
                            froms.push(self.make_with_or_from_value(value));
                        }
                    },
                    "gene_product_form_id" => {
                        if let Some(value) = prop.value.clone() {
                            extension.push(self.get_gene_prod_extension(&value));
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
                        let genotype_short = self.make_genotype_short(&feature.uniquename);
                        maybe_genotype_uniquename = Some(genotype_short.uniquename.clone());
                        genotype_short.expressed_alleles.iter()
                            .map(|expressed_allele| {
                                let allele_short =
                                    self.make_allele_short(&expressed_allele.allele_uniquename);
                                allele_short.gene_uniquename.clone()
                            })
                            .collect()
                    },
                    _ => {
                        if feature.feat_type.name == "gene" || feature.feat_type.name == "pseudogene" {
                            vec![feature.uniquename.clone()]
                        } else {
                            if TRANSCRIPT_FEATURE_TYPES.contains(&feature.feat_type.name.as_str()) {
                                if let Some(gene_uniquename) =
                                    self.genes_of_transcripts.get(&feature.uniquename) {
                                        vec![gene_uniquename.clone()]
                                    } else {
                                        vec![]
                                    }
                            } else {
                                vec![]
                            }
                        }
                    }
                };

            gene_uniquenames_vec.dedup();

            gene_uniquenames_vec =
                gene_uniquenames_vec.iter().map(|gene_uniquename: &String| {
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
                        copies_per_cell: copies_per_cell,
                        avg_copies_per_cell: avg_copies_per_cell,
                        scale: scale,
                    })
                } else {
                    None
                };

            if gene_uniquenames_vec.len() > 1 && maybe_genotype_uniquename.is_none() {
                panic!("non-genotype annotation has more than one gene");
            }

            let annotation = OntAnnotationDetail {
                id: feature_cvterm.feature_cvterm_id,
                genes: gene_uniquenames_vec,
                reference: reference_uniquename.clone(),
                genotype: maybe_genotype_uniquename.clone(),
                withs: withs.clone(),
                froms: froms.clone(),
                residue: extra_props_clone.remove("residue"),
                gene_ex_props: gene_ex_props,
                qualifiers: qualifiers.clone(),
                evidence: evidence.clone(),
                conditions: conditions.clone(),
                extension: extension.clone(),
            };

            self.add_annotation(cvterm.borrow(), feature_cvterm.is_not,
                                annotation);
        }
    }

    fn make_term_annotations(&self, termid: &str, detail_ids: &Vec<OntAnnotationId>,
                             is_not: bool)
                       -> Vec<(CvName, OntTermAnnotations)> {
        let term_details = self.terms.get(termid).unwrap();

        let cv_name = term_details.cv_name.clone();

        match cv_name.as_ref() {
            "gene_ex" => {
                if is_not {
                    panic!("gene_ex annotations can't be NOT annotations");
                }
                let mut qual_annotations =
                    OntTermAnnotations {
                        term: termid.to_owned(),
                        is_not: false,
                        rel_names: HashSet::new(),
                        annotations: vec![],
                        summary: None,
                    };
                let mut quant_annotations =
                    OntTermAnnotations {
                        term: termid.to_owned(),
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

                if qual_annotations.annotations.len() > 0 {
                    return_vec.push((String::from("qualitative_gene_expression"),
                                     qual_annotations));
                }

                if quant_annotations.annotations.len() > 0 {
                    return_vec.push((String::from("quantitative_gene_expression"),
                                     quant_annotations));
                }

                return_vec
            },
            "fission_yeast_phenotype" => {
                let mut single_allele =
                    OntTermAnnotations {
                        term: termid.to_owned(),
                        is_not: is_not,
                        rel_names: HashSet::new(),
                        annotations: vec![],
                        summary: None,
                    };
                let mut multi_allele =
                    OntTermAnnotations {
                        term: termid.to_owned(),
                        is_not: is_not,
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

                if single_allele.annotations.len() > 0 {
                    return_vec.push((String::from("single_allele_phenotype"),
                                     single_allele));
                }

                if multi_allele.annotations.len() > 0 {
                    return_vec.push((String::from("multi_allele_phenotype"),
                                     multi_allele));
                }

                return_vec
            },
            _ => {
                vec![(cv_name,
                      OntTermAnnotations {
                          term: termid.to_owned(),
                          is_not: is_not,
                          rel_names: HashSet::new(),
                          annotations: detail_ids.clone(),
                          summary: None,
                      })]
            }
        }
    }

    // store the OntTermAnnotations in the TermDetails, GeneDetails,
    // GenotypeDetails and ReferenceDetails
    fn store_ont_annotations(&mut self, is_not: bool) {
        let ont_annotations = if is_not {
            &self.all_not_ont_annotations
        } else {
            &self.all_ont_annotations
        };

        let mut gene_annotation_by_term: HashMap<GeneUniquename, HashMap<TermId, Vec<OntAnnotationId>>> =
            HashMap::new();
        let mut genotype_annotation_by_term: HashMap<GenotypeUniquename, HashMap<TermId, Vec<OntAnnotationId>>> =
            HashMap::new();
        let mut ref_annotation_by_term: HashMap<String, HashMap<TermId, Vec<OntAnnotationId>>> =
            HashMap::new();

        for (termid, annotations) in ont_annotations {
            let mut sorted_annotations = annotations.clone();

            if !is_not {
                let cv_config = {
                    let term = self.terms.get(termid).unwrap();
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
                                                          &self.genotypes, &self.alleles,
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
                            .or_insert(Vec::new())
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
                        .or_insert(HashMap::new())
                        .entry(termid.clone())
                        .or_insert(vec![])
                        .push(annotation_id);
                }

                if let Some(ref genotype_uniquename) = annotation.genotype {
                    let existing =
                        genotype_annotation_by_term.entry(genotype_uniquename.clone())
                        .or_insert(HashMap::new())
                        .entry(termid.clone())
                        .or_insert(vec![]);
                    if !existing.contains(&annotation_id) {
                        existing.push(annotation_id);
                    }
                }

                if let Some(reference_uniquename) = annotation.reference.clone() {
                    ref_annotation_by_term.entry(reference_uniquename)
                        .or_insert(HashMap::new())
                        .entry(termid.clone())
                        .or_insert(vec![])
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
                                        is_not: is_not,
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
                                .entry(extension_cv_name.clone())
                                .or_insert({
                                    let mut new_vec = Vec::new();
                                    let new_term_annotation =
                                        OntTermAnnotations {
                                            term: part_termid.to_owned(),
                                            is_not: is_not,
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
            let term1 = term_names.get(&ont_term_1.term).unwrap();
            let term2 = term_names.get(&ont_term_2.term).unwrap();

            term1.cmp(&term2)
        };

        for (gene_uniquename, term_annotation_map) in &gene_annotation_by_term {
            for (termid, details) in term_annotation_map {
                let new_annotations =
                    self.make_term_annotations(&termid, &details, is_not);

                let gene_details = self.genes.get_mut(gene_uniquename).unwrap();

                for (cv_name, new_annotation) in new_annotations {
                    gene_details.cv_annotations.entry(cv_name.clone())
                        .or_insert(Vec::new())
                        .push(new_annotation);
                }
            }

            let gene_details = self.genes.get_mut(gene_uniquename).unwrap();
            for (_, cv_annotations) in &mut gene_details.cv_annotations {
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
                        .or_insert(Vec::new())
                        .push(new_annotation);
                }
            }

            let details = self.genotypes.get_mut(genotype_uniquename).unwrap();
            for (_, cv_annotations) in &mut details.cv_annotations {
                cv_annotations.sort_by(&ont_term_cmp)
            }
        }

        for (reference_uniquename, ref_annotation_map) in &ref_annotation_by_term {
            for (termid, details) in ref_annotation_map {
                let new_annotations =
                    self.make_term_annotations(&termid, &details, is_not);

                let ref_details = self.references.get_mut(reference_uniquename).unwrap();

                for (cv_name, new_annotation) in new_annotations {
                    ref_details.cv_annotations.entry(cv_name).or_insert(Vec::new())
                        .push(new_annotation.clone());
                }
            }

            let ref_details = self.references.get_mut(reference_uniquename).unwrap();
            for (_, cv_annotations) in &mut ref_details.cv_annotations {
                cv_annotations.sort_by(&ont_term_cmp)
            }
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

                if subject_term_details.cv_annotations.keys().len() > 0 {
                    if let Some(object_term_details) = self.terms.get(&object_termid) {
                        if object_term_details.cv_annotations.keys().len() > 0 {
                            children_by_termid
                                .entry(object_termid.clone())
                                .or_insert(HashSet::new())
                                .insert(subject_termid.clone());
                        }
                    }
                }

                for (cv_name, term_annotations) in &subject_term_details.cv_annotations {
                    for term_annotation in term_annotations {
                        for annotation_id in &term_annotation.annotations {
                            let dest_termid = object_termid.clone();
                            let source_termid = subject_termid.clone();

                            if !term_annotation.is_not {
                                new_annotations.entry((cv_name.clone(), dest_termid))
                                    .or_insert(HashMap::new())
                                    .entry(source_termid)
                                    .or_insert(HashMap::new())
                                    .entry(*annotation_id)
                                    .or_insert(HashSet::new())
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
                let mut all_rel_names: HashSet<String> = HashSet::new();
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
                                                          &self.genotypes, &self.alleles, &self.terms);
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
                        .or_insert(Vec::new())
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
            let term1 = term_names.get(&ont_term_1.term).unwrap();
            let term2 = term_names.get(&ont_term_2.term).unwrap();

            term1.cmp(&term2)
        };

        for (_, term_details) in &mut self.terms {
            for (_, term_annotations) in &mut term_details.cv_annotations {
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
            export_prog_name: String::from(PKG_NAME),
            export_prog_version: String::from(VERSION),
            db_creation_datetime: db_creation_datetime.unwrap(),
            gene_count: self.genes.len(),
            term_count: self.terms.len(),
            cv_versions: cv_versions,
        }
    }

    pub fn get_api_genotype_annotation(&self) -> HashMap<TermId, Vec<APIGenotypeAnnotation>>
    {
        let mut app_genotype_annotation = HashMap::new();

        'TERM: for (_, term_details) in &self.terms {
            for (_, annotations_vec) in &term_details.cv_annotations {
                for ont_term_annotations in annotations_vec {
                    'DETAILS: for annotation_id in ont_term_annotations.annotations.iter() {
                        let annotation_details = self.annotation_details.
                            get(&annotation_id).expect("can't find OntAnnotationDetail");

                        if annotation_details.genotype.is_none() {
                            continue 'DETAILS;
                        }
                        let genotype_uniquename = annotation_details.genotype.clone().unwrap();
                        let genotype =
                            term_details.genotypes_by_uniquename.get(&genotype_uniquename).unwrap();
                        let mut api_annotation = APIGenotypeAnnotation {
                            is_multi: genotype.expressed_alleles.len() > 1,
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
                            .or_insert(vec![])
                            .push(api_annotation);
                    }
                }
            }
        }

        app_genotype_annotation
    }

    pub fn make_api_maps(mut self) -> APIMaps {
        let mut gene_summaries: HashMap<GeneUniquename, APIGeneSummary> = HashMap::new();
        let mut gene_name_gene_map = HashMap::new();
        let mut interactors_of_genes = HashMap::new();

        for (gene_uniquename, gene_details) in &self.genes {
            if self.config.load_organism_taxonid == gene_details.taxonid {
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
                        interactor_uniquename: interactor_uniquename,
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
                        interactor_uniquename: interactor_uniquename,
                    };
                    if !interactors.contains(&interactor) {
                        interactors.push(interactor);
                    }
                }
                interactors_of_genes.insert(gene_uniquename.clone(), interactors);
            }
        }

        let mut term_summaries: HashSet<TermShort> = HashSet::new();
        let mut termid_genes: HashMap<TermId, HashSet<GeneUniquename>> = HashMap::new();

        let mut terms_for_api: HashMap<TermId, TermDetails> = HashMap::new();

        for ref termid in self.terms.keys() {
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

            terms_for_api.insert(termid.clone(), term_details.clone());
        }

        APIMaps {
            gene_summaries: gene_summaries,
            termid_genes: termid_genes,
            termid_genotype_annotation: termid_genotype_annotation,
            term_summaries: term_summaries,
            genes: self.genes.clone(),
            gene_name_gene_map: gene_name_gene_map,
            genotypes: self.genotypes.clone(),
            terms: terms_for_api.clone(),
            interactors_of_genes: interactors_of_genes,
            references: self.references.clone(),
            other_features: self.other_features.clone(),
            annotation_details: self.annotation_details.clone(),
       }
    }

    fn add_cv_annotations_to_maps(&self,
                                  identifier: &String,
                                  cv_annotations: &OntAnnotationMap,
                                  seen_references: &mut HashMap<String, ReferenceShortOptionMap>,
                                  seen_genes: &mut HashMap<String, GeneShortOptionMap>,
                                  seen_genotypes: &mut HashMap<String, GenotypeShortMap>,
                                  seen_alleles: &mut HashMap<String, AlleleShortMap>,
                                  seen_terms: &mut HashMap<String, TermShortOptionMap>) {
        for (_, feat_annotations) in cv_annotations {
            for feat_annotation in feat_annotations.iter() {
                self.add_term_to_hash(seen_terms,
                                      identifier.clone(), feat_annotation.term.clone());

                for annotation_detail_id in &feat_annotation.annotations {
                    let annotation_detail = self.annotation_details.
                        get(&annotation_detail_id).expect("can't find OntAnnotationDetail");

                    self.add_ref_to_hash(seen_references,
                                         identifier.clone(), annotation_detail.reference.clone());
                    for condition_termid in &annotation_detail.conditions {
                        self.add_term_to_hash(seen_terms,
                                              identifier.clone(), condition_termid.clone());
                    }
                    for ext_part in &annotation_detail.extension {
                        match ext_part.ext_range {
                            ExtRange::Term(ref range_termid) =>
                                self.add_term_to_hash(seen_terms, identifier.clone(), range_termid.clone()),
                            ExtRange::Gene(ref allele_gene_uniquename) =>
                                self.add_gene_to_hash(seen_genes, identifier.clone(),
                                                      allele_gene_uniquename.clone()),
                            _ => {},
                        }
                    }
                    if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                        self.add_genotype_to_hash(seen_genotypes, seen_alleles, seen_genes,
                                                  identifier.clone(),
                                                  &genotype_uniquename);
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
                    self.add_term_to_hash(&mut seen_terms,
                                          termid.clone(), term_annotation.term.clone());
                    for annotation_detail_id in &term_annotation.annotations {
                        let annotation_detail = self.annotation_details
                            .get(&annotation_detail_id).expect("can't find OntAnnotationDetail");

                        for gene_uniquename in &annotation_detail.genes {
                            self.add_gene_to_hash(&mut seen_genes, termid.clone(), gene_uniquename.clone());
                            if !cv_name.starts_with("extension:") {
                                // prevent extension annotations from appears
                                // in the normal query builder searches
                                genes_annotated_with_map
                                    .entry(termid.clone()).or_insert(HashSet::new())
                                    .insert(gene_uniquename.clone());
                            }
                        }
                        self.add_ref_to_hash(&mut seen_references, termid.clone(),
                                             annotation_detail.reference.clone());
                        for condition_termid in &annotation_detail.conditions {
                            self.add_term_to_hash(&mut seen_terms, termid.clone(), condition_termid.clone());
                        }
                        for ext_part in &annotation_detail.extension {
                            match ext_part.ext_range {
                                ExtRange::Term(ref range_termid) =>
                                    self.add_term_to_hash(&mut seen_terms, termid.clone(), range_termid.clone()),
                                ExtRange::Gene(ref ext_gene_uniquename) =>
                                    self.add_gene_to_hash(&mut seen_genes, termid.clone(),
                                                          ext_gene_uniquename.clone()),
                                _ => {},
                            }
                        }
                        if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                            self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles,
                                                      &mut seen_genes, termid.clone(),
                                                      &genotype_uniquename);
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
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename.clone(), interaction.reference_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), interaction.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), interaction.interactor_uniquename.clone());
                }

                for ortholog_annotation in &gene_details.ortholog_annotations {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename.clone(), ortholog_annotation.reference_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), ortholog_annotation.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), ortholog_annotation.ortholog_uniquename.clone());
                }
                for paralog_annotation in &gene_details.paralog_annotations {
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename.clone(), paralog_annotation.reference_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), paralog_annotation.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(), paralog_annotation.paralog_uniquename.clone());
                }
                for target_of_annotation in &gene_details.target_of_annotations {
                    for annotation_gene_uniquename in &target_of_annotation.genes {
                        self.add_gene_to_hash(&mut seen_genes, gene_uniquename.clone(),
                                              annotation_gene_uniquename.clone());
                    }
                    if let Some(ref annotation_genotype_uniquename) = target_of_annotation.genotype_uniquename {
                        self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
                                                  gene_uniquename.clone(),
                                                  &annotation_genotype_uniquename.clone())
                    }
                    self.add_ref_to_hash(&mut seen_references, gene_uniquename.clone(),
                                         target_of_annotation.reference_uniquename.clone());
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
        let mut seen_genes: HashMap<String, GeneShortOptionMap> = HashMap::new();

        type GenotypeShortMap = HashMap<GenotypeUniquename, GenotypeShort>;
        let mut seen_genotypes: HashMap<ReferenceUniquename, GenotypeShortMap> = HashMap::new();

        type AlleleShortMap = HashMap<AlleleUniquename, AlleleShort>;
        let mut seen_alleles: HashMap<TermId, AlleleShortMap> = HashMap::new();

        let mut seen_terms: HashMap<GeneUniquename, TermShortOptionMap> = HashMap::new();

        {
            for (reference_uniquename, reference_details) in &self.references {
                for (_, feat_annotations) in &reference_details.cv_annotations {
                    for feat_annotation in feat_annotations.iter() {
                        self.add_term_to_hash(&mut seen_terms, reference_uniquename.clone(),
                                              feat_annotation.term.clone());

                       for annotation_detail_id in &feat_annotation.annotations {
                            let annotation_detail = self.annotation_details
                                .get(&annotation_detail_id).expect("can't find OntAnnotationDetail");
                            for gene_uniquename in &annotation_detail.genes {
                                self.add_gene_to_hash(&mut seen_genes, reference_uniquename.clone(),
                                                      gene_uniquename.clone())
                            }
                            for condition_termid in &annotation_detail.conditions {
                                self.add_term_to_hash(&mut seen_terms, reference_uniquename.clone(), condition_termid.clone());
                            }
                            for ext_part in &annotation_detail.extension {
                                match ext_part.ext_range {
                                    ExtRange::Term(ref range_termid) =>
                                        self.add_term_to_hash(&mut seen_terms, reference_uniquename.clone(), range_termid.clone()),
                                    ExtRange::Gene(ref allele_gene_uniquename) =>
                                        self.add_gene_to_hash(&mut seen_genes, reference_uniquename.clone(),
                                                              allele_gene_uniquename.clone()),
                                    _ => {},
                                }
                            }
                            if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                                let genotype = self.make_genotype_short(genotype_uniquename);
                                self.add_genotype_to_hash(&mut seen_genotypes, &mut seen_alleles, &mut seen_genes,
                                                          reference_uniquename.clone(),
                                                          &genotype.uniquename);
                            }
                        }
                    }
                }

                let interaction_iter =
                    reference_details.physical_interactions.iter().chain(&reference_details.genetic_interactions);
                for interaction in interaction_iter {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename.clone(), interaction.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename.clone(), interaction.interactor_uniquename.clone());
                }

                for ortholog_annotation in &reference_details.ortholog_annotations {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename.clone(), ortholog_annotation.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename.clone(), ortholog_annotation.ortholog_uniquename.clone());
                }
                for paralog_annotation in &reference_details.paralog_annotations {
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename.clone(), paralog_annotation.gene_uniquename.clone());
                    self.add_gene_to_hash(&mut seen_genes, reference_uniquename.clone(), paralog_annotation.paralog_uniquename.clone());
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
            for (_, term_annotations) in &term_details.cv_annotations {
                for term_annotation in term_annotations {
                    for annotation_detail_id in &term_annotation.annotations {
                        let annotation_detail = self.annotation_details
                            .get(&annotation_detail_id).expect("can't find OntAnnotationDetail");
                        for gene_uniquename in &annotation_detail.genes {
                            seen_genes.insert(gene_uniquename.clone());
                        }
                        if let Some(ref genotype_uniquename) = annotation_detail.genotype {
                            seen_genotypes.insert(genotype_uniquename.clone());
                            let genotype = self.genotypes.get(genotype_uniquename).unwrap();
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
            for (_, rel_annotations) in &reference_details.cv_annotations {
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

        for (_, term_details) in &mut self.terms {
            term_details.single_allele_genotype_uniquenames =
                term_seen_single_allele_genotypes.remove(&term_details.termid).unwrap();

            term_details.gene_count =
                term_seen_genes.get(&term_details.termid).unwrap().len();
            term_details.genotype_count =
                term_seen_genotypes.get(&term_details.termid).unwrap().len();
        }
    }

    fn make_non_bp_slim_gene_subset(&self, go_slim_subset: &TermSubsetDetails)
                                    -> IdGeneSubsetMap
    {
        let slim_termid_set: HashSet<String> =
            go_slim_subset.elements
            .iter().map(|element| element.termid.clone()).collect();

        let mut non_slim_with_bp_annotation = HashSet::new();
        let mut non_slim_without_bp_annotation = HashSet::new();

        let has_parent_in_slim = |term_annotations: &Vec<OntTermAnnotations>| {
            for term_annotation in term_annotations {
                let ref interesting_parents =
                    self.terms.get(&term_annotation.term).unwrap().interesting_parents;
                if !term_annotation.is_not &&
                    (slim_termid_set.contains(&term_annotation.term) ||
                     interesting_parents.intersection(&slim_termid_set).count() > 0)
                {
                    return true;
                }
            }
            false
        };

        'GENE: for (_, gene_details) in &self.genes {
            if self.config.load_organism_taxonid != gene_details.taxonid {
                continue;
            }

            if gene_details.feature_type != "mRNA gene" {
                continue;
            }

            if gene_details.characterisation_status == Some("transposon".into()) ||
                gene_details.characterisation_status == Some("dubious".into())
            {
                continue;
            }

            let mut bp_count = 0;

            if let Some(annotations) =
                gene_details.cv_annotations.get("biological_process") {
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

        return_map.insert("non_go_slim_with_bp_annotation".into(),
                          GeneSubsetDetails {
                              name: "non_go_slim_with_bp_annotation".into(),
                              display_name: String::from("Proteins with biological process ") +
                                  "annotation that are not in a slim category",
                              elements: non_slim_with_bp_annotation,
                          });
        return_map.insert("non_go_slim_without_bp_annotation".into(),
                          GeneSubsetDetails {
                              name: "non_go_slim_without_bp_annotation".into(),
                              display_name: String::from("Proteins with no biological process ") +
                                  "annotation and are not in a slim category",
                              elements: non_slim_without_bp_annotation,
                          });
        return_map
    }

    fn make_bp_go_slim_subset(&self) -> TermSubsetDetails {
        let mut all_genes = HashSet::new();
        let mut go_slim_subset: HashSet<TermSubsetElement> = HashSet::new();
        'TERM: for go_slim_conf in self.config.go_slim_terms.clone() {
            let slim_termid = go_slim_conf.termid;
            let term_details = self.terms.get(&slim_termid)
                .expect(&format!("can't find TermDetails for {}", &slim_termid));

            let subset_element = TermSubsetElement {
                name: term_details.name.clone(),
                termid: slim_termid.clone(),
                gene_count: term_details.genes_annotated_with.len(),
            };

            for gene in &term_details.genes_annotated_with {
                all_genes.insert(gene);
            }
            go_slim_subset.insert(subset_element);
        }

        TermSubsetDetails {
            name: "goslim_pombe".into(),
            total_gene_count: all_genes.len(),
            elements: go_slim_subset,
        }
    }

    fn make_feature_type_subsets(&self, subsets: &mut IdGeneSubsetMap) {
        for (_, gene_details) in &self.genes {
            if self.config.load_organism_taxonid != gene_details.taxonid {
                continue;
            }

            let subset_name =
                String::from("feature_type:") + &gene_details.feature_type;
            let re = Regex::new(r"[\s,:]+").unwrap();
            let subset_name_no_spaces: String = re.replace_all(&subset_name, "_").into();
            subsets.entry(subset_name_no_spaces.clone())
                .or_insert(GeneSubsetDetails {
                    name: subset_name_no_spaces,
                    display_name: subset_name,
                    elements: HashSet::new()
                })
                .elements.insert(gene_details.uniquename.clone());
        }
    }

    // make subsets using the characterisation_status field of GeneDetails
    fn make_characterisation_status_subsets(&self, subsets: &mut IdGeneSubsetMap) {
        for (_, gene_details) in &self.genes {
            if self.config.load_organism_taxonid != gene_details.taxonid {
                continue;
            }

            if gene_details.feature_type != "mRNA gene" {
                continue;
            }

            if let Some(ref characterisation_status) = gene_details.characterisation_status {
                let subset_name =
                    String::from("characterisation_status:") + &characterisation_status;
                let re = Regex::new(r"[\s,:]+").unwrap();
                let subset_name_no_spaces: String = re.replace_all(&subset_name, "_").into();
                subsets.entry(subset_name_no_spaces.clone())
                    .or_insert(GeneSubsetDetails {
                        name: subset_name_no_spaces,
                        display_name: subset_name,
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

                if interpro_match.interpro_id.len() > 0 {
                    let subset_name =
                        String::from("interpro:") + &interpro_match.interpro_id;
                    new_subset_names.push((subset_name,
                                           interpro_match.interpro_name.clone()));
                }

                let subset_name = String::from("interpro:") +
                     &interpro_match.dbname.clone() + ":" + &interpro_match.id;
                new_subset_names.push((subset_name, interpro_match.name.clone()));

                for (subset_name, display_name) in new_subset_names {
                    subsets.entry(subset_name.clone())
                        .or_insert(GeneSubsetDetails {
                            name: subset_name,
                            display_name: display_name,
                            elements: HashSet::new(),
                        })
                        .elements.insert(gene_uniquename.clone());
                }
            }
        }
    }

    // populated the subsets HashMap
    fn make_subsets(&mut self) {
        let bp_go_slim_subset = self.make_bp_go_slim_subset();
        let mut gene_subsets =
            self.make_non_bp_slim_gene_subset(&bp_go_slim_subset);

        self.term_subsets.insert("bp_goslim_pombe".into(), bp_go_slim_subset);

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
                let gene1 = self.genes.get(uniquename1).unwrap();
                let gene2 = self.genes.get(uniquename2).unwrap();

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

            for (_, gene_uniquenames) in &mut genes_to_sort {
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
        self.add_alleles_to_genotypes();
        self.process_cvterms();
        self.add_interesting_parents();
        self.process_cvterm_rels();
        self.process_extension_cvterms();
        self.process_feature_synonyms();
        self.process_feature_cvterms();
        self.store_ont_annotations(false);
        self.store_ont_annotations(true);
        self.process_cvtermpath();
        self.process_annotation_feature_rels();
        self.add_target_of_annotations();
        self.set_deletion_viability();
        self.set_term_details_subsets();
        self.make_all_cv_summaries();
        self.remove_non_curatable_refs();
        self.set_term_details_maps();
        self.set_gene_details_maps();
        self.set_genotype_details_maps();
        self.set_reference_details_maps();
        self.set_counts();
        self.make_subsets();
        self.sort_chromosome_genes();

        let metadata = self.make_metadata();

        let mut gene_summaries: Vec<GeneSummary> = vec![];

        for (gene_uniquename, gene_details) in &self.genes {
            if self.config.load_organism_taxonid == gene_details.taxonid {
                gene_summaries.push(self.make_gene_summary(&gene_uniquename));
            }
        }

        let mut solr_term_summaries = vec![];

        let term_name_split_re = Regex::new(r"\W+").unwrap();

        for (termid, term_details) in self.terms.iter() {
            if term_details.cv_annotations.keys().len() == 0 {
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
            let mut close_synonym_words_vec: Vec<String> = vec![];
            let mut distant_synonyms = vec![];
            let mut distant_synonym_words_vec: Vec<String> = vec![];

            let add_to_words_vec = |synonym: &str, words_vec: &mut Vec<String>| {
                let synonym_words = term_name_split_re.split(&synonym);
                for word in synonym_words {
                    let word_string = word.trim_matches(&trimmable_p).to_owned();
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
                close_synonyms: close_synonyms,
                close_synonym_words: close_synonym_words_vec.join(" "),
                distant_synonyms: distant_synonyms,
                distant_synonym_words: distant_synonym_words_vec.join(" "),
                interesting_parents: interesting_parents_for_solr,
            };
            solr_term_summaries.push(term_summ);
        }

        let solr_data = SolrData {
            term_summaries: solr_term_summaries,
            gene_summaries: gene_summaries.clone(),
        };

        let chromosomes = self.chromosomes.clone();
        let mut chromosome_summaries = vec![];

        for (_, chr_details) in &self.chromosomes {
            chromosome_summaries.push(chr_details.make_chromosome_short());
        }

        let term_subsets = self.term_subsets.clone();
        let gene_subsets = self.gene_subsets.clone();
        let recent_references = self.recent_references.clone();
        let all_community_curated = self.all_community_curated.clone();

        WebData {
            metadata: metadata,
            chromosomes: chromosomes,
            chromosome_summaries: chromosome_summaries,
            recent_references: recent_references,
            all_community_curated: all_community_curated,
            api_maps: self.make_api_maps(),
            search_gene_summaries: gene_summaries,
            term_subsets: term_subsets,
            gene_subsets: gene_subsets,
            solr_data: solr_data,
        }
    }
}

