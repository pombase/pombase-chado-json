use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap, HashSet};

use regex::Regex;


use crate::types::GeneUniquename;
use crate::data_types::{AlleleDetails, AlleleShort, DisplayUniquenameGenotypeMap,
                        ExtPart, ExtRange, GeneDetails, IdOntAnnotationDetailMap,
                        ProteinDetails, ProteinViewData, ProteinViewFeature,
                        ProteinViewFeaturePos, ProteinViewTrack, TermIdDetailsMap,
                        TermNameAndId, UniquenameAlleleDetailsMap, UniquenameGeneMap,
                        UniquenameReferenceMap, UniquenameTranscriptMap};
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
            annotated_terms: BTreeSet::new(),
            feature_group: None,
            display_extension: BTreeSet::new(),
            assigned_by: None,
            author_and_year: None,
            evidence: None,
            inviable_or_abnormal: None,
            positions,
        })
    } else {
        None
    }
}

fn make_mutant_summary(mutants_track: &ProteinViewTrack) -> ProteinViewTrack {
    let mut summary_features = HashSet::new();

    for mutant_feat in &mutants_track.features {
        for (id, feat_start, feat_end) in &mutant_feat.positions {
            for (residue, pos) in id.chars().zip(*feat_start..*feat_end+1) {
                let key = (flex_fmt!("{}{}", residue, pos), pos);

                summary_features.insert(key);
            }
        }
    }

    let features: Vec<_> = summary_features.into_iter()
        .map(|(residue_and_pos, pos)| {
            ProteinViewFeature {
               id: residue_and_pos.clone(),
               display_name: Some(residue_and_pos.clone()),
               annotated_terms: BTreeSet::new(),
               feature_group: None,
               display_extension: BTreeSet::new(),
               assigned_by: None,
               author_and_year: None,
               evidence: None,
               inviable_or_abnormal: None,
               positions: vec![(residue_and_pos, pos, pos)],
            }
         })
         .collect();

    ProteinViewTrack {
        name: flex_str!("AA substitution positions"),
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
            flex_str!("AA substitution alleles")
        } else {
            flex_str!("Partial deletions")
        },
        display_type: flex_str!("block"),
        features: vec![],
    };

    let mut all_categories = HashSet::new();

    if let Some(term_filter) = pheno_type_config.filters.iter().find(|cat| cat.filter_type == "term") {
        for term_category in &term_filter.term_categories {
            for ancestor in &term_category.ancestors {
                all_categories.insert(ancestor.clone());
            }
        }
    }

    let mut inviable_or_abnormal_termids = HashSet::new();
    inviable_or_abnormal_termids.extend(["FYPO:0001985".into(), // inviable cell
                                         "FYPO:0002059".into(), // inviable cell population
                                         "FYPO:0000049".into(), // abnormal phenotype
                                        ].into_iter());

    all_categories.extend(inviable_or_abnormal_termids.clone());

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
            let mut inviable_or_abnormal = false;

            let mut categories_with_names: Vec<_> = allele_categories
                .iter()
                .filter_map(|cat_term_id: &&FlexStr| {
                    let term_details = term_details_map.get(*cat_term_id)?;

                    Some(TermNameAndId {
                        id: (*cat_term_id).clone(),
                        name: term_details.name.clone(),
                    })
                })
                .collect();

            for cat_with_name in &categories_with_names {
                if inviable_or_abnormal_termids.contains(&cat_with_name.id) {
                    inviable_or_abnormal = true;
                }
            }

            categories_with_names.sort_by(|c1, c2| c1.cmp(c2));

            feature.annotated_terms.extend(categories_with_names);

            feature.inviable_or_abnormal = Some(inviable_or_abnormal);

            track.features.push(feature);
        }
    }

    track
}

fn make_mod_extension(extension: &Vec<ExtPart>,
                      ext_rel_types: &HashSet<FlexStr>,
                      gene_details_maps: &UniquenameGeneMap,
                      term_details_map: &TermIdDetailsMap) -> Vec<FlexStr> {

    let ext_display_range = |ext_part: &ExtPart| {
        match &ext_part.ext_range {
            ExtRange::GeneAndGeneProduct(gene_and_prod) => {
                if let Some(gene_details) = gene_details_maps.get(&gene_and_prod.gene_uniquename) {
                    let gene_display_name = gene_details.display_name();
                    flex_fmt!("{} / {}", gene_display_name, gene_and_prod.product)
                } else {
                    flex_fmt!("{} / {}", gene_and_prod.gene_uniquename, gene_and_prod.product)
                }
            },
            ExtRange::Term(termid) => {
                let term_details = term_details_map.get(termid).unwrap();
                flex_fmt!("{} ({})", term_details.name, termid)
            },
            ExtRange::Gene(gene_uniquename) => {
                if let Some(gene_details) = gene_details_maps.get(gene_uniquename) {
                    gene_details.display_name()
                } else {
                    gene_uniquename.clone()
                }
            },
            _ => { panic!("unknown relation: {:#?}", ext_part)}
        }
    };

    let mut display_extension: Vec<FlexStr> = extension
        .iter()
        .filter(|ext_part| {
            ext_rel_types.contains(&ext_part.rel_type_name)
        })
        .map(|ext_part| {
            let range_str = ext_display_range(ext_part);
            flex_fmt!("{} {}", ext_part.rel_type_display_name, range_str)
        })
        .collect();

    display_extension.sort_by(|s1, s2| s1.cmp(s2));

    display_extension
}

lazy_static! {
    static ref MODIFICATION_RESIDUE_RE: Regex =
        Regex::new(r"^([ARNDCQEGHILKMFPOSUTWYVBZXJ]*(\d+))$").unwrap();
}

fn make_modification_track(gene_details: &GeneDetails,
                           config: &Config,
                           gene_details_maps: &UniquenameGeneMap,
                           term_details_map: &TermIdDetailsMap,
                           references_map: &UniquenameReferenceMap,
                           annotation_details_map: &IdOntAnnotationDetailMap) -> ProteinViewTrack {
    let ext_rel_types = &config.protein_feature_view.modification_extension_rel_types;

    let mut seen_modifications = HashMap::new();

    if let Some(term_annotations) = gene_details.cv_annotations.get("PSI-MOD") {

        for term_annotation in term_annotations {
            let termid = &term_annotation.term;
            let ref term_details = term_details_map
                .get(termid)
                .expect(&format!("internal error, can't find term: {}", termid));
            let term_name = &term_details.name;
            let annotations = term_annotation.annotations.clone();

            for annotation_id in annotations {
                let annotation_detail = annotation_details_map.get(&annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                let assigned_by = &annotation_detail.assigned_by;
                let mut author_and_year = None;

                if let Some(ref reference_uniquename) = annotation_detail.reference {
                    if let Some(ref ref_details) = references_map.get(reference_uniquename) {
                        if let Some(ref authors_abbrev) = ref_details.authors_abbrev {
                            if let Some(ref year) = ref_details.publication_year {
                                author_and_year = Some(flex_fmt!("{} ({})", authors_abbrev, year));
                            } else {
                                author_and_year = Some(authors_abbrev.clone());
                            }
                        }
                    }
                }

                let evidence = &annotation_detail.evidence;

                let mut annotation_residues = vec![];

                if let Some(ref residue_str) = annotation_detail.residue {
                    if let Some(psi_mod_config) = config.cv_config.get("PSI-MOD") {
                        let mod_abbrev_conf = &psi_mod_config.modification_abbreviations;
                        let gene_uniquename = &gene_details.uniquename;

                        if let Some(res_config) = mod_abbrev_conf.get(gene_uniquename) {
                            if let Some(expanded_residues) = res_config.get(residue_str) {
                                annotation_residues.extend(expanded_residues.split(","));
                            }
                        }
                    }

                    if annotation_residues.len() == 0 {
                        annotation_residues.push(residue_str.as_str());
                    }
                };

                for residue_str in &annotation_residues {
                    let Some(captures) = MODIFICATION_RESIDUE_RE.captures(residue_str)
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

                    let feature = seen_modifications
                        .entry(description.clone())
                        .or_insert_with(|| {

                            let mut feature_group = None;

                            'GROUP: for mod_group in &config.protein_feature_view.modification_groups {
                                if &mod_group.termid == termid {
                                    feature_group = Some(mod_group.clone());
                                    break 'GROUP;
                                } else {
                                    for interesting_parent in &term_details.interesting_parent_ids {
                                        if interesting_parent == mod_group.termid {
                                            feature_group = Some(mod_group.clone());
                                            break 'GROUP;
                                        }
                                    }
                                }
                            }

                            ProteinViewFeature {
                                id: description.clone(),
                                display_name: Some(description.clone()),
                                annotated_terms: BTreeSet::new(),
                                feature_group,
                                display_extension: BTreeSet::new(),
                                assigned_by: assigned_by.clone(),
                                author_and_year: author_and_year.clone(),
                                evidence: evidence.clone(),
                                inviable_or_abnormal: None,
                                positions: vec![(description, residue_pos, residue_pos)],
                            }
                        });

                    let name_and_id = TermNameAndId {
                        name: term_name.clone(),
                        id: termid.clone(),
                    };

                    feature.annotated_terms.insert(name_and_id);

                    let display_extension =
                        make_mod_extension(&annotation_detail.extension, ext_rel_types,
                                           gene_details_maps, term_details_map);

                    feature.display_extension.extend(display_extension);
                }
            }
        }
    }

    let features = seen_modifications.into_iter()
                      .map(|(_, feature)| feature)
                      .collect();

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
                annotated_terms: BTreeSet::new(),
                feature_group: None,
                assigned_by: Some(flex_str!["InterPro"]),
                evidence: None,
                author_and_year: None,
                display_extension: BTreeSet::new(),
                inviable_or_abnormal: None,
                positions,
            };

            features.push(feature);
        }
    }

    ProteinViewTrack {
        name: flex_str!("Pfam domains"),
        display_type: flex_str!("block"),
        features,
    }
}

fn make_generic_track(track_name: FlexStr, feature_coords: &Vec<(usize, usize)>,
                      split_start_and_end: bool)
    -> ProteinViewTrack
{
    let features =
        feature_coords.iter().map(|(start, end)| {
            let feature_name = flex_fmt!("{} {}..{}", track_name, start, end);
            let feature_pos_name = flex_fmt!("{}-{}..{}", track_name, start, end);
            let positions =
                if split_start_and_end {
                    vec![(feature_pos_name.clone(), *start, *start),
                         (feature_pos_name, *end, *end)]
                } else {
                    vec![(feature_pos_name, *start, *end)]
                };
            ProteinViewFeature {
                id: feature_name.clone(),
                display_name: Some(feature_name.clone()),
                annotated_terms: BTreeSet::new(),
                feature_group: None,
                display_extension: BTreeSet::new(),
                assigned_by: None,
                author_and_year: None,
                evidence: None,
                inviable_or_abnormal: None,
                positions,
            }
        })
        .collect();

    ProteinViewTrack {
        name: track_name,
        display_type: flex_str!("block"),
        features,
    }
}

fn make_binding_sites_track(gene_details: &GeneDetails) -> ProteinViewTrack {
    let track_name = flex_str!("Binding sites");
    let features =
        gene_details.binding_sites.iter().map(|binding_site| {
            let ligand = &binding_site.ligand;
            let start = binding_site.range.start;
            let end = binding_site.range.end;
            let feature_name = flex_fmt!("Binding site, ligand: {}", ligand);
            let feature_id = flex_fmt!("Binding site {}..{}, ligand: {}", start, end, ligand);
            ProteinViewFeature {
                id: feature_id,
                display_name: Some(feature_name.clone()),
                annotated_terms: BTreeSet::new(),
                feature_group: None,
                display_extension: BTreeSet::new(),
                assigned_by: None,
                author_and_year: None,
                evidence: None,
                inviable_or_abnormal: None,
                positions: vec![(feature_name.clone(), start, end)],
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
        p1_start.cmp(&p2_start)
                .then_with(|| p1_end.cmp(&p2_end))
    };

    let sorter = |f1: &ProteinViewFeature, f2: &ProteinViewFeature| {
        for (p1, p2) in f1.positions.iter().zip(f2.positions.iter()) {
            let cmp = sort_helper(p1, p2);

            if cmp != Ordering::Equal {
                return cmp;
            }
        }

        f1.positions.len().cmp(&f2.positions.len())
    };

    deletions_track.features.sort_by(sorter);
}

pub fn make_protein_view_data_map(gene_details_maps: &UniquenameGeneMap,
                                  term_details_map: &TermIdDetailsMap,
                                  annotation_details_map: &IdOntAnnotationDetailMap,
                                  genotypes: &DisplayUniquenameGenotypeMap,
                                  alleles: &UniquenameAlleleDetailsMap,
                                  transcripts: &UniquenameTranscriptMap,
                                  references: &UniquenameReferenceMap,
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
            make_modification_track(gene_details, config, gene_details_maps, term_details_map,
                                    references, annotation_details_map);

        let pfam_track = make_pfam_track(gene_details);

        let tm_domains_track =
            make_generic_track(flex_str!("TM domains"),
                               &gene_details.tm_domain_coords, false);
        let disordered_regions_track =
            make_generic_track(flex_str!("Disordered regions"),
                               &gene_details.disordered_region_coords, false);
        let low_complexity_regions_track =
            make_generic_track(flex_str!("Low complexity"),
                               &gene_details.low_complexity_region_coords, false);
        let coiled_coil_coords =
            make_generic_track(flex_str!("Coiled coils"), &gene_details.coiled_coil_coords, false);

        let mut signal_peptide_coords = vec![];

        if let Some(ref signal_peptide) = gene_details.signal_peptide {
            signal_peptide_coords.push((signal_peptide.range.start,
                                        signal_peptide.range.end));
        }

        let signal_peptide_track =
            make_generic_track(flex_str!("Signal peptide"), &signal_peptide_coords, false);

        let mut transit_peptide_coords = vec![];

        if let Some(ref transit_peptide) = gene_details.transit_peptide {
            transit_peptide_coords.push((transit_peptide.range.start,
                                        transit_peptide.range.end));
        }

        let transit_peptide_track =
            make_generic_track(flex_str!("Transit peptide"), &transit_peptide_coords, false);

        let binding_sites_track = make_binding_sites_track(gene_details);

        let active_sites_coords =
            gene_details.active_sites.iter()
                .map(|site| (site.range.start, site.range.end))
                .collect();
        let active_sites_track =
            make_generic_track(flex_str!("Active sites"), &active_sites_coords, false);

        let beta_strand_coords =
            gene_details.beta_strands.iter()
                .map(|strand| (strand.range.start, strand.range.end))
                .collect();
        let beta_strands_track =
            make_generic_track(flex_str!("Beta strands"), &beta_strand_coords, false);

        let helix_coords =
            gene_details.helices.iter()
                .map(|helix| (helix.range.start, helix.range.end))
                .collect();
        let helix_track =
            make_generic_track(flex_str!("Helices"), &helix_coords, false);

        let turn_coords =
            gene_details.turns.iter()
                .map(|turn| (turn.range.start, turn.range.end))
                .collect();
        let turns_track =
            make_generic_track(flex_str!("Turns"), &turn_coords, false);

        let propeptide_coords =
            gene_details.propeptides.iter()
                .map(|propeptide| (propeptide.range.start, propeptide.range.end))
                .collect();
        let propeptides_track =
            make_generic_track(flex_str!("Propeptides"), &propeptide_coords, false);

        let chain_coords =
            gene_details.chains.iter()
                .map(|chain| (chain.range.start, chain.range.end))
                .collect();
        let chains_track =
            make_generic_track(flex_str!("Chains"), &chain_coords, false);

        let disulfide_bond_coords =
            gene_details.disulfide_bonds.iter()
                .map(|disulfide_bond| (disulfide_bond.range.start, disulfide_bond.range.end))
                .collect();
        let disulfide_bonds_track =
            make_generic_track(flex_str!("Disulfide bonds"), &disulfide_bond_coords, true);

        let protein_view_data = ProteinViewData {
            sequence: protein.sequence.clone(),
            tracks: vec![mutant_summary_track, mutants_track, deletions_track,
                         modification_track,
                         disulfide_bonds_track, pfam_track,
                         tm_domains_track, disordered_regions_track,
                         low_complexity_regions_track, coiled_coil_coords,
                         signal_peptide_track, transit_peptide_track,
                         binding_sites_track, active_sites_track,
                         beta_strands_track, helix_track, turns_track,
                         propeptides_track, chains_track],
        };

        gene_map.insert(gene_details.uniquename.clone(), protein_view_data);
    }

    gene_map
}
