use std::cmp::Ordering;

use crate::types::*;
use crate::data_types::*;
use crate::web::config::*;
use crate::web::cmp_utils::cmp_residues;
use crate::web::util::*;

use pombase_rc_string::RcString;

fn gene_display_name(gene: &GeneDetails) -> RcString {
    if let Some(ref name) = gene.name {
        name.clone()
    } else {
        gene.uniquename.clone()
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
        ExtRange::Transcript(ref transcript_uniquename) => transcript_uniquename.clone(),
        ExtRange::SummaryGenes(_) => panic!("can't handle SummaryGenes\n"),
        ExtRange::Term(ref termid) => RcString::from(&terms.get(termid).unwrap().name),
        ExtRange::SummaryModifiedResidues(ref residue) => RcString::from(&residue.join(",")),
        ExtRange::SummaryTerms(_) => panic!("can't handle SummaryGenes\n"),
        ExtRange::Misc(ref misc) => misc.clone(),
        ExtRange::Domain(ref domain) => domain.clone(),
        ExtRange::GeneProduct(ref gene_product) => gene_product.clone(),
        ExtRange::GeneAndGeneProduct(ref gene_and_gene_product) =>
            RcString::from(&format!("{} ({})", gene_and_gene_product.gene_uniquename,
                                    gene_and_gene_product.product)),
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
    ext1_for_cmp.extend(ext1_rest);

    let (mut ext2_for_cmp, ext2_rest): (Vec<ExtPart>, Vec<ExtPart>) =
        ext2.to_vec().into_iter().partition(&is_grouping_rel_name);
    ext2_for_cmp.extend(ext2_rest);

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

pub fn sort_cv_annotation_details<T: AnnotationContainer>
    (container: &mut T,
     config: &Config,
     gene_map: &UniquenameGeneMap,
     genotype_map: &UniquenameGenotypeMap,
     term_map: &TermIdDetailsMap,
     annotation_details_map: &IdOntAnnotationDetailMap)
{

    let detail_cmp_using_cv_name = |cv_name: &str| {
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
                                              gene_map,
                                              genotype_map,
                                              term_map);
                result.unwrap_or_else(|err| panic!("error from cmp_ont_annotation_detail: {}", err))
            } else {
                Ordering::Equal
            }
        }
    };

    for (_, term_annotations) in container.cv_annotations_mut() {
        for term_annotation in term_annotations {
            let termid = &term_annotation.term;
            if let Some(term_details) = term_map.get(termid) {
                let cmp = detail_cmp_using_cv_name(&term_details.cv_name);
                term_annotation.annotations.sort_by(cmp);
            } else {
                panic!("failed to find TermDetails for {}", termid);
            }
        }
    }
}
