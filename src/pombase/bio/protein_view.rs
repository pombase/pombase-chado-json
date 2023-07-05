use std::collections::{HashMap, HashSet};

use regex::Regex;


use crate::types::GeneUniquename;
use crate::data_types::{ProteinViewData, UniquenameGeneMap,
                        IdOntAnnotationDetailMap, DisplayUniquenameGenotypeMap,
                        UniquenameAlleleDetailsMap, AlleleShort,
                        ProteinViewFeature, ProteinViewTrack,
                        UniquenameTranscriptMap};
use crate::web::config::Config;

use flexstr::{shared_str as flex_str, SharedStr as FlexStr, shared_fmt as flex_fmt};


lazy_static! {
    static ref MUTATION_DESC_RE: Regex =
       Regex::new(r"^([ARNDCQEGHILKMFPOSUTWYVBZXJ]+)(\d+)[ARNDCQEGHILKMFPOSUTWYVBZXJ]+$").unwrap();
}


fn variant_from_allele(allele: &AlleleShort) -> Option<ProteinViewFeature> {
    let Some(ref allele_description) = allele.description
    else {
        return None;
    };

    let Some(captures) = MUTATION_DESC_RE.captures(allele_description.as_str())
    else {
        return None;
    };

    let (Some(amino_acids), Some(pos_match)) = (captures.get(1), captures.get(2))
    else {
        return None;
    };

    let mut positions = vec![];

    let Ok(pos) = pos_match.as_str().parse::<usize>()
    else {
        return None;
    };

    let end_pos = pos + amino_acids.len() - 1;
    positions.push((pos, end_pos));

    let mutation_details = ProteinViewFeature {
         id: allele.uniquename.clone(),
         display_name: Some(allele.display_name()),
         positions,
    };

    Some(mutation_details)
}

fn make_generic_track(track_name: FlexStr, feature_coords: &Vec<(usize, usize)>)
    -> ProteinViewTrack
{
    let features =
        feature_coords.iter().map(|(start, end)| {
            let feature_name = flex_fmt!("{} {}..{}", track_name, start, end);
            ProteinViewFeature {
                 id: feature_name.clone(),
                 display_name: None,
                 positions: vec![(*start, *end)],
            }
        })
        .collect();

    ProteinViewTrack {
        name: track_name + "s",
        display_type: flex_str!("block"),
        features,
    }
}


pub fn make_protein_view_data_map(gene_details_maps: &UniquenameGeneMap,
                                  annotation_details_maps: &IdOntAnnotationDetailMap,
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

    let variants_name = flex_str!("Variants");

    for gene_details in gene_details_maps.values() {
        if gene_details.taxonid != load_org_taxonid {
            continue;
        }

        let Some(term_annotations) = gene_details.cv_annotations.get("single_locus_phenotype")
        else {
            continue;
        };

        let Some(transcript_uniquename) = gene_details.transcripts.get(0)
        else {
            continue;
        };

        let transcript = transcripts.get(transcript_uniquename).unwrap();

        let Some(ref protein) = transcript.protein
        else {
            continue;
        };

        let mut variant_track = ProteinViewTrack {
            name: variants_name.clone(),
            display_type: flex_str!("block"),
            features: vec![],
        };

        let mut seen_alleles = HashSet::new();

        for term_annotation in term_annotations {
            let annotations = term_annotation.annotations.clone();

            for annotation_id in annotations {
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
                        if seen_alleles.contains(&expressed_allele.allele_uniquename) {
                            continue;
                        }

                        seen_alleles.insert(expressed_allele.allele_uniquename.clone());

                        let Some(allele_details) = alleles.get(&expressed_allele.allele_uniquename)
                        else {
                            continue;
                        };

                        let allele_short: AlleleShort = allele_details.into();

                        let Some(variant_features) = variant_from_allele(&allele_short)
                        else {
                            continue;
                        };

                        variant_track.features.push(variant_features);
                    }
                }
            }
        }

        let tm_domains_track =
            make_generic_track(flex_str!("TM domain"),
                               &gene_details.tm_domain_coords);
        let disordered_regions_track =
            make_generic_track(flex_str!("Disordered region"),
                               &gene_details.disordered_region_coords);
        let low_complexity_regions_track =
            make_generic_track(flex_str!("Low complexity region"),
                               &gene_details.low_complexity_region_coords);
        let coiled_coil_coords =
            make_generic_track(flex_str!("Coiled coil"), &gene_details.coiled_coil_coords);

        let protein_view_data = ProteinViewData {
            sequence: protein.sequence.clone(),
            tracks: vec![tm_domains_track, disordered_regions_track,
                         low_complexity_regions_track, coiled_coil_coords,
                         variant_track],
        };

        gene_map.insert(gene_details.uniquename.clone(), protein_view_data);
    }

    gene_map
}
