use std::cmp::Ordering;

use crate::types::*;
use crate::data_types::*;
use crate::utils::join;
use crate::web::config::*;
use crate::web::cmp_utils::cmp_residues;
use crate::web::util::*;

use flexstr::{SharedStr as FlexStr, shared_fmt as flex_fmt};

fn gene_display_name(gene: &GeneDetails) -> FlexStr {
    if let Some(ref name) = gene.name {
        name.clone()
    } else {
        gene.uniquename.clone()
    }
}

fn string_from_ext_range(ext_range: &ExtRange,
                         data_lookup: &dyn DataLookup) -> FlexStr {

    let gene_name = |gene_uniquename: &FlexStr| {
        let gene = data_lookup.get_gene(gene_uniquename)
           .unwrap_or_else(|| panic!("can't find gene: {}", gene_uniquename));
        gene_display_name(gene.as_ref())
    };

    match *ext_range {
        ExtRange::Gene(ref gene_uniquename) | ExtRange::Promoter(ref gene_uniquename) => {
            gene_name(gene_uniquename)
        },
        ExtRange::Transcript(ref transcript_uniquename) => transcript_uniquename.clone(),
        ExtRange::SummaryGenes(ref summ_genes) => {
            gene_name(&summ_genes[0][0])
        },
        ExtRange::SummaryTranscripts(ref summ_transcripts) => summ_transcripts[0][0].clone(),
        ExtRange::Term(ref termid) => data_lookup.get_term(termid).unwrap().name.clone(),
        ExtRange::ModifiedResidues(ref residue) => join(residue, ","),
        ExtRange::SummaryTerms(ref summ_terms) => {
            data_lookup.get_term(&summ_terms[0]).unwrap().name.clone()
        },
        ExtRange::Misc(ref misc) => misc.clone(),
        ExtRange::Domain(ref domain) => domain.clone(),
        ExtRange::GeneProduct(ref gene_product) => gene_product.clone(),
        ExtRange::GeneAndGeneProduct(ref gene_and_gene_product) =>
            flex_fmt!("{} ({})", gene_and_gene_product.gene_uniquename,
                      gene_and_gene_product.product),
    }
}

fn cmp_ext_part(ext_part1: &ExtPart, ext_part2: &ExtPart,
                data_lookup: &dyn DataLookup) -> Ordering {

    let ord = ext_part1.rel_type_display_name.cmp(&ext_part2.rel_type_display_name);

    if ord == Ordering::Equal {
        let ext_part1_str = string_from_ext_range(&ext_part1.ext_range, data_lookup);
        let ext_part2_str = string_from_ext_range(&ext_part2.ext_range, data_lookup);

        ext_part1_str.to_lowercase().cmp(&ext_part2_str.to_lowercase())
    } else {
        ord
    }
}

// compare the extension up to the last common index
fn cmp_extension_prefix(conf_rel_ranges: &Vec<FlexStr>, ext1: &[ExtPart], ext2: &[ExtPart],
                        data_lookup: &dyn DataLookup) -> Ordering {
    let mut ext1_for_cmp = ext1.to_owned();
    let mut ext2_for_cmp = ext2.to_owned();

    let is_grouping_rel_name =
        |ext: &ExtPart| conf_rel_ranges.contains(&ext.rel_type_name);

    // put the extension that will be grouped in the summary at the end
    // See: https://github.com/pombase/pombase-chado/issues/636

    let first_ext1_collect_rel_idx = ext1.iter().position(is_grouping_rel_name);

    if let Some(ext1_collect_rel_idx) = first_ext1_collect_rel_idx {
        let ext1_cellect_rel = ext1_for_cmp.remove(ext1_collect_rel_idx);
        ext1_for_cmp.push(ext1_cellect_rel);
    }

    let first_ext2_collect_rel_idx = ext2.iter().position(is_grouping_rel_name);

    if let Some(ext2_collect_rel_idx) = first_ext2_collect_rel_idx {
        let ext2_cellect_rel = ext2_for_cmp.remove(ext2_collect_rel_idx);
        ext2_for_cmp.push(ext2_cellect_rel);
    }

    let iter = ext1_for_cmp.iter().zip(&ext2_for_cmp).enumerate();
    for (_, (ext1_part, ext2_part)) in iter {
        let ord = cmp_ext_part(ext1_part, ext2_part, data_lookup);

        if ord != Ordering::Equal {
            return ord
        }
    }

    Ordering::Equal
}

pub fn cmp_extension(conf_rel_ranges: &Vec<FlexStr>, ext1: &[ExtPart], ext2: &[ExtPart],
                     data_lookup: &dyn DataLookup) -> Ordering {
    let cmp = cmp_extension_prefix(conf_rel_ranges, ext1, ext2, data_lookup);
    if cmp == Ordering::Equal {
        ext1.len().cmp(&ext2.len())
    } else {
        cmp
    }
}

fn cmp_genotypes(genotype1: &GenotypeDetails, genotype2: &GenotypeDetails) -> Ordering {
    genotype1.display_uniquename.to_lowercase().cmp(&genotype2.display_uniquename.to_lowercase())
}


// compare two gene vectors which must be ordered vecs
fn cmp_gene_vec(data_lookup: &dyn DataLookup,
                gene_vec1: &[GeneUniquename],
                gene_vec2: &[GeneUniquename]) -> Ordering {

    let gene_short_vec1: Vec<GeneShort> =
        gene_vec1.iter().map(|gene_uniquename: &FlexStr| {
            make_gene_short(data_lookup, gene_uniquename)
        }).collect();
    let gene_short_vec2: Vec<GeneShort> =
        gene_vec2.iter().map(|gene_uniquename: &FlexStr| {
            make_gene_short(data_lookup, gene_uniquename)
        }).collect();

    gene_short_vec1.cmp(&gene_short_vec2)
}


pub fn cmp_ont_annotation_detail(cv_config: &CvConfig,
                                 detail1: &OntAnnotationDetail,
                                 detail2: &OntAnnotationDetail,
                                 data_lookup: &dyn DataLookup) -> Result<Ordering, String> {
    let conf_rel_ranges = &cv_config.summary_relation_ranges_to_collect;
    if let Some(ref detail1_genotype_uniquename) = detail1.genotype {
        if let Some(ref detail2_genotype_uniquename) = detail2.genotype {
            let genotype1 = data_lookup.get_genotype(detail1_genotype_uniquename).unwrap();
            let genotype2 = data_lookup.get_genotype(detail2_genotype_uniquename).unwrap();

            let ord = cmp_genotypes(&genotype1, &genotype2);

            if ord == Ordering::Equal {
                Ok(cmp_extension(conf_rel_ranges, &detail1.extension, &detail2.extension,
                                 data_lookup))
            } else {
                Ok(ord)
            }
        } else {
            // needed to display the FYECO term pages
            Ok(Ordering::Less)
        }
    } else if detail2.genotype.is_some() {
        // needed to display the FYECO term pages
        Ok(Ordering::Greater)
    } else {
        let ord = cmp_gene_vec(data_lookup, &detail1.genes, &detail2.genes);

        if ord == Ordering::Equal {
            if let Some(ref sort_details_by) = cv_config.sort_details_by {
                for sort_type in sort_details_by {
                    if sort_type == "modification" {
                        let res = cmp_residues(&detail1.residue, &detail2.residue);
                        if res != Ordering::Equal {
                            return Ok(res);
                        }
                    } else {
                        let res = cmp_extension(conf_rel_ranges, &detail1.extension,
                                                &detail2.extension,
                                                data_lookup);
                        if res != Ordering::Equal {
                            return Ok(res);
                        }
                    }
                }
                Ok(Ordering::Equal)
            } else {
                Ok(cmp_extension(conf_rel_ranges, &detail1.extension, &detail2.extension,
                                 data_lookup))
            }
        } else {
            Ok(ord)
        }
    }
}

pub fn sort_cv_annotation_details<T: AnnotationContainer>
    (container: &mut T,
     config: &Config,
     data_lookup: &dyn DataLookup,
     annotation_details_map: &IdOntAnnotationDetailMap)
{

    let detail_cmp_using_cv_name = |cv_name: &FlexStr| {
        let cv_config = config.cv_config.get(cv_name);

        move |id1: &i32, id2: &i32| {
            let annotation1: &OntAnnotationDetail =
                annotation_details_map.get(id1).expect("can't find OntAnnotationDetail");
            let annotation2: &OntAnnotationDetail =
                annotation_details_map.get(id2).expect("can't find OntAnnotationDetail");


            if let Some(cv_config) = cv_config {
                let result =
                    cmp_ont_annotation_detail(cv_config,
                                              annotation1, annotation2,
                                              data_lookup);
                result.unwrap_or_else(|err| panic!("error from cmp_ont_annotation_detail: {}", err))
            } else {
                Ordering::Equal
            }
        }
    };

    for term_annotations in container.cv_annotations_mut().values_mut() {
        for term_annotation in term_annotations {
            let termid = &term_annotation.term;
            if let Some(term_details) = data_lookup.get_term(termid) {
                let cmp = detail_cmp_using_cv_name(&term_details.cv_name);
                term_annotation.annotations.sort_by(cmp);
            } else {
                panic!("failed to find TermDetails for {}", termid);
            }
        }
    }
}
