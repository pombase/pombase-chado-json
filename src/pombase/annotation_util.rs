use std::collections::{HashSet, HashMap};

use flexstr::AFlexStr as FlexStr;

use crate::web::config::{CvConfig, AnnotationSubsetConfig};
use crate::types::{CvName};
use crate::data_types::APIMaps;
use crate::utils::join;


pub fn table_for_export(api_maps: &APIMaps, cv_config_map: &HashMap<CvName, CvConfig>,
                        subset_config: &AnnotationSubsetConfig)
    -> Vec<Vec<FlexStr>>
{
    let mut seen: HashSet<Vec<FlexStr>> = HashSet::new();

    let mut result: Vec<Vec<FlexStr>> = vec![];

    for termid in &subset_config.term_ids {
        let term_details = api_maps.terms.get(&termid)
            .unwrap_or_else(|| panic!("no term details found for {} for config file", termid));

        for (cv_name, term_annotations) in &term_details.cv_annotations {
            if let Some(cv_config) = cv_config_map.get(cv_name) {
                if subset_config.single_or_multi_locus != cv_config.single_or_multi_locus {
                    continue;
                }
            } else {
                panic!("can't find configuration for CV: {}", cv_name);
            }

            for term_annotation in term_annotations {
                let termid = &term_annotation.term;

                let annotation_term_details = api_maps.terms.get(termid).unwrap();

                if term_annotation.is_not {
                    continue;
                }

                for annotation_id in &term_annotation.annotations {
                    let mut row = vec![];
                    let annotation_details = api_maps.annotation_details
                        .get(annotation_id).expect("can't find OntAnnotationDetail");

                    let gene_uniquenames = &annotation_details.genes;

                    let maybe_genotype_short =
                        if let Some(ref genotype_uniquename) = annotation_details.genotype {
                            api_maps.genotypes.get(genotype_uniquename)
                        } else {
                            None
                        };

                    for column_config in &subset_config.columns {
                        if column_config.name == "cv_name" {
                            row.push(cv_name.clone());
                        }
                        if column_config.name == "termid" {
                            row.push(termid.clone());
                        }
                        if column_config.name == "term_name" {
                            let term_name = annotation_term_details.name.clone();
                            row.push(term_name);
                        }
                        if column_config.name == "allele" {
                            if let Some(genotype_short) = maybe_genotype_short {
                                row.push(genotype_short.display_uniquename.clone());
                            }
                        }
                        if column_config.name == "gene_uniquename" {
                            let gene_uniquenames_string = join(&gene_uniquenames, ",");
                            row.push(gene_uniquenames_string);
                        }
                        if column_config.name == "gene_name" {
                            
                            let filtered_gene_uniquenames = gene_uniquenames.iter()
                                .filter_map(|uniquename| {
                                   api_maps.genes.get(uniquename).unwrap().name.clone()
                                }).collect::<Vec<_>>();
                       
                            row.push(join(&filtered_gene_uniquenames, ","));
                        }
                    }

                    if !seen.contains(&row) {
                        result.push(row.clone());
                        seen.insert(row.clone());
                    }
                }
            }
        }
    }

    result
}
