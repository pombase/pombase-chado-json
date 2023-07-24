use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use regex::Regex;


use crate::types::GeneUniquename;
use crate::data_types::{ProteinViewData, UniquenameGeneMap,
                        IdOntAnnotationDetailMap, DisplayUniquenameGenotypeMap,
                        UniquenameAlleleDetailsMap, AlleleShort,
                        ProteinViewFeature, ProteinViewTrack,
                        UniquenameTranscriptMap, GeneDetails,
                        TermIdDetailsMap, AlleleDetails, ProteinDetails,
                        ProteinViewFeaturePos, TermNameAndId, ExtPart};
use crate::web::config::{Config, CvConfig};

use flexstr::{shared_str as flex_str, SharedStr as FlexStr, shared_fmt as flex_fmt};

#[derive(PartialEq)]
enum TrackType {
    PartialDeletions,
    AminoAcidMutations,
}

lazy_static! {
    static ref MUTATION_DESC_RE: Regex =
        Regex::new(r"^([ARNDCQEGHILKMFPOSUTWYVBZXJ]+)-?(\d+)-?[ARNDCQEGHILKMFPOSUTWYVBZXJ]+$").unwrap();
}

fn parse_mutation_postion(desc_part: &str) -> Option<(FlexStr, usize, usize)> {
    let Some(captures) = MUTATION_DESC_RE.captures(desc_part)
    else {
        return None;
    };

    let (Some(amino_acids), Some(pos_match)) = (captures.get(1), captures.get(2))
    else {
        return None;
    };

    let Ok(pos) = pos_match.as_str().parse::<usize>()
    else {
        return None;
    };

    let end_pos = pos + amino_acids.len() - 1;
    let pos_name = flex_fmt!("{}-{}..{}", desc_part, pos, end_pos);

    Some((pos_name, pos, end_pos))
}

lazy_static! {
    static ref DELETION_DESC_RE: Regex =
       Regex::new(r"^(\d+)-(\d+)$").unwrap();
}

fn parse_deletion_postion(desc_part: &str)
        -> Option<(FlexStr, usize, usize)>
{
    let Some(captures) = DELETION_DESC_RE.captures(desc_part)
    else {
        return None;
    };

    let (Some(start), Some(end)) = (captures.get(1), captures.get(2))
    else {
        return None;
    };

    let Ok(start) = start.as_str().parse::<usize>()
    else {
        return None;
    };

    let Ok(end) = end.as_str().parse::<usize>()
    else {
        return None;
    };

    let pos_name = flex_fmt!("{}-{}..{}", desc_part, start, end);

    Some((pos_name, start, end))
}

lazy_static! {
    static ref TRUNCATION_DESC_RE: Regex =
       Regex::new(r"^[ARNDCQEGHILKMFPOSUTWYVBZXJ](\d+)\*$").unwrap();
}

fn parse_truncation_postion(desc_part: &str, seq_lenth: usize)
        -> Option<(FlexStr, usize, usize)>
{
    let Some(captures) = TRUNCATION_DESC_RE.captures(desc_part)
    else {
        return None;
    };

    let Some(start) = captures.get(1)
    else {
        return None;
    };

    let Ok(start) = start.as_str().parse::<usize>()
    else {
        return None;
    };

    let pos_name = flex_fmt!("{}-{}..{}", desc_part, start, seq_lenth);

    Some((pos_name, start, seq_lenth))
}

fn feature_from_allele(allele_details: &AlleleDetails, seq_length: usize)
       -> Option<ProteinViewFeature>
{
    let allele: AlleleShort = allele_details.into();

    let Some(ref allele_description) = allele.description
    else {
        return None;
    };

    let mut positions = vec![];

    for desc_part in allele_description.split(",") {
        if let Some(mutation_position) = parse_mutation_postion(desc_part) {
            positions.push(mutation_position);
            continue;
        }
        if let Some(deletion_position) = parse_deletion_postion(desc_part) {
            positions.push(deletion_position);
            continue;
        } else {
            if let Some(deletion_position) = parse_truncation_postion(desc_part, seq_length) {
                positions.push(deletion_position);
                continue;
            }
        }
    }

    if positions.len() > 0 {
        Some(ProteinViewFeature {
            id: allele.uniquename.clone(),
            display_name: Some(allele.display_name()),
            annotated_terms: HashSet::new(),
            display_extension: vec![],
            positions,
        })
    } else {
        None
    }
}

fn make_mutant_summary(mutants_track: &ProteinViewTrack) -> ProteinViewTrack {
    let mut summary_features = HashMap::new();

    for mutant_feat in &mutants_track.features {
        for (_, feat_start, feat_end) in &mutant_feat.positions {
            let key = (*feat_start, *feat_end);
            summary_features.entry(key)
               .or_insert_with(|| {});
        }
    }

    let features = summary_features.keys()
        .map(|(start, end)| {
            let pos_str = if start == end {
                flex_fmt!("{}", start)
            } else {
                flex_fmt!("{}..{}", start, end)
            };
            let pos_display_name = flex_fmt!("AA {}", pos_str);
            ProteinViewFeature {
               id: pos_str.clone(),
               display_name: Some(pos_display_name),
               annotated_terms: HashSet::new(),
               display_extension: vec![],
               positions: vec![(pos_str, *start, *end)],
            }
         })
         .collect();

    ProteinViewTrack {
        name: flex_str!("Mutant positions"),
        display_type: flex_str!("block"),
        features,
    }
}

fn make_mutants_track(gene_details: &GeneDetails,
                      pheno_type_config: &CvConfig,
                      protein: &ProteinDetails,
                      track_type: TrackType,
                      annotation_details_maps: &IdOntAnnotationDetailMap,
                      term_details_map: &TermIdDetailsMap,
                      genotypes: &DisplayUniquenameGenotypeMap,
                      alleles: &UniquenameAlleleDetailsMap) -> ProteinViewTrack {
    let mut track = ProteinViewTrack {
        name: if track_type == TrackType::AminoAcidMutations {
            flex_str!("Mutants")
        } else {
            flex_str!("Partial deletions")
        },
        display_type: flex_str!("block"),
        features: vec![],
    };

    let mut all_categories = HashSet::new();

    if let Some(term_filter) = pheno_type_config.filters.iter().find(|cat| cat.filter_name == "term") {
        for term_category in &term_filter.term_categories {
            for ancestor in &term_category.ancestors {
                all_categories.insert(ancestor.clone());
            }
        }
    }

    let mut allele_terms = HashMap::new();

    if let Some(term_annotations) = gene_details.cv_annotations.get("single_locus_phenotype") {

        for term_annotation in term_annotations {
            let Some(term_details) = term_details_map.get(&term_annotation.term)
            else {
                continue;
            };

            let mut term_category_ids = HashSet::new();

            for category_termid in all_categories.iter() {
                if term_details.interesting_parent_ids.contains(category_termid) {
                    term_category_ids.insert(category_termid);
                }
            }

            for annotation_id in &term_annotation.annotations {
                let annotation_detail = annotation_details_maps.get(&annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                let Some(ref genotype_uniquename) = annotation_detail.genotype
                else {
                    continue;
                };

                let Some(genotype) = genotypes.get(genotype_uniquename)
                else {
                    continue;
                };

                for locus in &genotype.loci {
                    for expressed_allele in &locus.expressed_alleles {
                        if allele_terms.contains_key(&expressed_allele.allele_uniquename) {
                            continue;
                        }

                        let allele_categories =
                            allele_terms
                            .entry(expressed_allele.allele_uniquename.clone())
                            .or_insert_with(HashSet::new);

                        allele_categories.extend(term_category_ids.iter());
                    }
                }
            }
        }
    }



    for (allele_uniquename, allele_categories) in &allele_terms {

        let Some(allele_details) = alleles.get(allele_uniquename)
        else {
            continue;
        };

        match allele_details.allele_type.as_str() {
            "amino_acid_mutation" => {
                if track_type != TrackType::AminoAcidMutations {
                    continue;
                }
            },
            "partial_amino_acid_deletion" => {
                if track_type != TrackType::PartialDeletions {
                    continue;
                }
            },
            _ => {
                continue;
            }
        }

        let sequence_length = protein.sequence_length();

        if let Some(mut feature) = feature_from_allele(&allele_details, sequence_length) {
            let categories_with_names = allele_categories
                .iter()
                .filter_map(|cat_term_id: &&FlexStr| {
                    let term_details = term_details_map.get(*cat_term_id)?;

                    Some(TermNameAndId {
                        id: (*cat_term_id).clone(),
                        name: term_details.name.clone(),
                    })
                });

            feature.annotated_terms.extend(categories_with_names);

            track.features.push(feature);
        }
    }

    track
}

lazy_static! {
    static ref MODIFICATION_RESIDUE_RE: Regex =
        Regex::new(r"^([ARNDCQEGHILKMFPOSUTWYVBZXJ]+(\d+))$").unwrap();
}


fn make_modification_track(gene_details: &GeneDetails,
                           config: &Config,
                           term_details_map: &TermIdDetailsMap,
                           annotation_details_map: &IdOntAnnotationDetailMap) -> ProteinViewTrack {
    let mut features = vec![];

    let ext_rel_types = &config.protein_feature_view.modification_extension_rel_types;

    if let Some(term_annotations) = gene_details.cv_annotations.get("PSI-MOD") {
        let mut seen_modifications = HashSet::new();

        for term_annotation in term_annotations {
            let termid = &term_annotation.term;
            let ref term_name = term_details_map
                .get(termid)
                .expect(&format!("term: {}", termid)).name;
            let annotations = term_annotation.annotations.clone();

            for annotation_id in annotations {
                let annotation_detail = annotation_details_map.get(&annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                if let Some(ref residue_str) = annotation_detail.residue {
                    let Some(captures) = MODIFICATION_RESIDUE_RE.captures(residue_str.as_str())
                    else {
                        continue;
                    };

                    let (Some(residue_match), Some(pos_match)) = (captures.get(1), captures.get(2))
                    else {
                        continue;
                    };

                    let Ok(residue_pos) = pos_match.as_str().parse::<usize>()
                    else {
                        continue;
                    };

                    let description = flex_fmt!("{}: {}", term_name, residue_match.as_str());

                    if !seen_modifications.contains(&description) {
                        seen_modifications.insert(description.clone());

                        let mut annotated_terms = HashSet::new();
                        let name_and_id = TermNameAndId {
                            name: term_name.clone(),
                            id: termid.clone(),
                        };
                        annotated_terms.insert(name_and_id);

                        let display_extension: Vec<ExtPart> = annotation_detail.extension
                            .iter()
                            .filter(|ext_part| {
                                ext_rel_types.contains(&ext_part.rel_type_name)
                            })
                            .map(|ext_part| ext_part.clone())
                            .collect();

                        let feature = ProteinViewFeature {
                            id: description.clone(),
                            display_name: Some(description.clone()),
                            annotated_terms,
                            display_extension,
                            positions: vec![(description, residue_pos, residue_pos)],
                        };

                        features.push(feature);
                    }
                }
            }
        }
    }

    ProteinViewTrack {
        name: flex_str!("Modifications"),
        display_type: flex_str!("pin"),
        features,
    }
}

fn make_pfam_track(gene_details: &GeneDetails) -> ProteinViewTrack {
    let mut features = vec![];

    for interpro_match in &gene_details.interpro_matches {
        if interpro_match.dbname == "PFAM" {
            let positions = interpro_match
                .locations
                .iter()
                .map(|loc| (flex_fmt!("{}-{}..{}", interpro_match.name, loc.start, loc.end), loc.start, loc.end))
                .collect();

            let display_name =
                flex_fmt!("{}: {}", interpro_match.id, interpro_match.name);

            let feature = ProteinViewFeature {
                id: interpro_match.id.clone(),
                display_name: Some(display_name),
                annotated_terms: HashSet::new(),
                display_extension: vec![],
                positions,
            };

            features.push(feature);
        }
    }

    ProteinViewTrack {
        name: flex_str!("Pfam families"),
        display_type: flex_str!("block"),
        features,
    }
}

fn make_generic_track(track_name: FlexStr, feature_coords: &Vec<(usize, usize)>)
    -> ProteinViewTrack
{
    let features =
        feature_coords.iter().map(|(start, end)| {
            let feature_name = flex_fmt!("{} {}..{}", track_name, start, end);
            let feature_pos_name = flex_fmt!("{}-{}..{}", track_name, start, end);
            ProteinViewFeature {
                id: feature_name.clone(),
                display_name: None,
                annotated_terms: HashSet::new(),
                display_extension: vec![],
                positions: vec![(feature_pos_name, *start, *end)],
            }
        })
        .collect();

    ProteinViewTrack {
        name: track_name,
        display_type: flex_str!("block"),
        features,
    }
}

fn sort_deletions(deletions_track: &mut ProteinViewTrack) {
    let sort_helper = |(_, p1_start, p1_end): &ProteinViewFeaturePos,
                       (_, p2_start, p2_end): &ProteinViewFeaturePos| {
        let first_pos_cmp = p1_start.cmp(&p2_start);
        if first_pos_cmp == Ordering::Equal {
            p1_end.cmp(&p2_end)
        } else {
            first_pos_cmp
        }
    };

    let sorter = |f1: &ProteinViewFeature, f2: &ProteinViewFeature| {
        for (p1, p2) in f1.positions.iter().zip(f2.positions.iter()) {
            let cmp = sort_helper(p1, p2);

            if cmp != Ordering::Equal {
                return cmp;
            }
        }

        Ordering::Equal
    };

    deletions_track.features.sort_by(sorter);
}

pub fn make_protein_view_data_map(gene_details_maps: &UniquenameGeneMap,
                                  term_details_map: &TermIdDetailsMap,
                                  annotation_details_map: &IdOntAnnotationDetailMap,
                                  genotypes: &DisplayUniquenameGenotypeMap,
                                  alleles: &UniquenameAlleleDetailsMap,
                                  transcripts: &UniquenameTranscriptMap,
                                  config: &Config)
                                  -> HashMap<GeneUniquename, ProteinViewData>
{
    let mut gene_map = HashMap::new();

    let Some(load_org_taxonid) = config.load_organism_taxonid
    else {
        return gene_map;
    };

    let Some(phenotype_config) = config.cv_config.get("fission_yeast_phenotype")
    else {
        return gene_map;
    };

    for gene_details in gene_details_maps.values() {
        if gene_details.taxonid != load_org_taxonid {
            continue;
        }

        let Some(transcript_uniquename) = gene_details.transcripts.get(0)
        else {
            continue;
        };

        let transcript = transcripts.get(transcript_uniquename).unwrap();

        let Some(ref protein) = transcript.protein
        else {
            continue;
        };

        let mutants_track = make_mutants_track(gene_details, phenotype_config, protein,
                                               TrackType::AminoAcidMutations,
                                               annotation_details_map,
                                               term_details_map,
                                               genotypes, alleles);

        let mutant_summary_track = make_mutant_summary(&mutants_track);

        let mut deletions_track = make_mutants_track(gene_details, phenotype_config, protein,
                                                 TrackType::PartialDeletions,
                                                 annotation_details_map,
                                                 term_details_map,
                                                 genotypes, alleles);
        sort_deletions(&mut deletions_track);

        let modification_track =
            make_modification_track(gene_details, config, term_details_map, annotation_details_map);

        let pfam_track = make_pfam_track(gene_details);

        let tm_domains_track =
            make_generic_track(flex_str!("TM domains"),
                               &gene_details.tm_domain_coords);
        let disordered_regions_track =
            make_generic_track(flex_str!("Disordered regions"),
                               &gene_details.disordered_region_coords);
        let low_complexity_regions_track =
            make_generic_track(flex_str!("Low complexity"),
                               &gene_details.low_complexity_region_coords);
        let coiled_coil_coords =
            make_generic_track(flex_str!("Coiled coils"), &gene_details.coiled_coil_coords);

        let protein_view_data = ProteinViewData {
            sequence: protein.sequence.clone(),
            tracks: vec![mutant_summary_track, mutants_track, deletions_track,
                         modification_track, pfam_track,
                         tm_domains_track, disordered_regions_track,
                         low_complexity_regions_track, coiled_coil_coords],
        };

        gene_map.insert(gene_details.uniquename.clone(), protein_view_data);
    }

    gene_map
}
