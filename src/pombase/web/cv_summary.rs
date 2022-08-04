use std::mem;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

use std::collections::{HashSet, HashMap};

use crate::types::*;
use crate::data_types::*;
use crate::web::config::*;
use crate::web::vec_set::*;
use crate::web::util::*;

use flexstr::SharedStr as FlexStr;

pub fn make_cv_summaries<T: AnnotationContainer>
    (container: &mut T,
     config: &Config,
     children_by_termid: &HashMap<TermId, HashSet<TermId>>,
     include_gene: bool, include_genotype: bool,
     gene_map: &UniquenameGeneMap,
     genotype_map: &UniquenameGenotypeMap,
     annotation_details: &IdOntAnnotationDetailMap) {
    for (cv_name, mut term_annotations) in container.cv_annotations_mut() {
        let cv_config = config.cv_config_by_name(cv_name);
        make_cv_summary(&cv_config, children_by_termid,
                        include_gene, include_genotype, &mut term_annotations,
                        gene_map, genotype_map, annotation_details);
    }
}

// merge two ExtPart objects into one by merging ranges
pub fn merge_ext_part_ranges(ext_part1: &ExtPart, ext_part2: &ExtPart,
                             genes: &UniquenameGeneMap) -> ExtPart {
    if ext_part1.rel_type_display_name == ext_part2.rel_type_display_name {
        match ext_part1.ext_range {
            ExtRange::SummaryGenes(ref part1_summ_genes) => {
                if let ExtRange::SummaryGenes(ref part2_summ_genes) = ext_part2.ext_range {
                    let mut ret_ext_part = ext_part1.clone();
                    let mut new_gene_uniquenames = [part1_summ_genes.clone(), part2_summ_genes.clone()].concat();
                    let cmp =
                        |vec1: &Vec<FlexStr>, vec2: &Vec<FlexStr>| {
                            let gene1 = &genes[&vec1[0]];
                            let gene2 = &genes[&vec2[0]];
                            gene1.cmp(gene2)
                        };
                    new_gene_uniquenames.sort_by(cmp);
                    new_gene_uniquenames.dedup();
                    ret_ext_part.ext_range = ExtRange::SummaryGenes(new_gene_uniquenames);
                    return ret_ext_part
                }
            },
            ExtRange::SummaryTranscripts(ref part1_summ_transcripts) => {
                if let ExtRange::SummaryTranscripts(ref part2_summ_transcripts) = ext_part2.ext_range {
                    let mut ret_ext_part = ext_part1.clone();
                    let mut new_transcript_uniquenames =
                        [part1_summ_transcripts.clone(), part2_summ_transcripts.clone()].concat();
                    let cmp =
                        |vec1: &Vec<FlexStr>, vec2: &Vec<FlexStr>| {
                            let transcript_uniquename_1 = &vec1[0];
                            let transcript_uniquename_2 = &vec2[0];
                            transcript_uniquename_1.cmp(transcript_uniquename_2)
                        };
                    new_transcript_uniquenames.sort_by(cmp);
                    new_transcript_uniquenames.dedup();
                    ret_ext_part.ext_range = ExtRange::SummaryTranscripts(new_transcript_uniquenames);
                    return ret_ext_part
                }
            },
            ExtRange::SummaryTerms(ref part1_summ_termids) => {
                if let ExtRange::SummaryTerms(ref part2_summ_termids) = ext_part2.ext_range {
                    let mut ret_ext_part = ext_part1.clone();
                    let mut new_terms =
                        [part1_summ_termids.clone(), part2_summ_termids.clone()].concat();
                    new_terms.sort();
                    new_terms.dedup();
                    ret_ext_part.ext_range = ExtRange::SummaryTerms(new_terms);
                    return ret_ext_part
                }
            },
            ExtRange::ModifiedResidues(ref part1_summ_residues) => {
                if let ExtRange::ModifiedResidues(ref part2_summ_residues) = ext_part2.ext_range
                {
                    let mut ret_ext_part = ext_part1.clone();
                    let mut new_terms =
                        [part1_summ_residues.clone(), part2_summ_residues.clone()].concat();
                    new_terms.sort();
                    new_terms.dedup();
                    ret_ext_part.ext_range = ExtRange::ModifiedResidues(new_terms);
                    return ret_ext_part
                }
            },
            _ => () // fall through and panic
        }
        panic!("passed ExtPart objects that have ranges that aren't genes or terms
to merge_ext_part_ranges(): {:?} {:?}", ext_part1, ext_part2);
    } else {
        panic!("passed ExtPart objects with mismatched relations to merge_ext_part_ranges():
  {} {}\n", ext_part1.rel_type_name, ext_part2.rel_type_name);
    }
}

// turn "has_substrate(gene1)" "has_substrate(gene2)" into "has_substrate(gene1,gene2)"
pub fn collect_ext_summary_genes(cv_config: &CvConfig, rows: &mut Vec<TermSummaryRow>,
                                 genes: &UniquenameGeneMap) {

    let conf_rel_ranges = &cv_config.summary_relation_ranges_to_collect;
    let merge_range_rel_p =
        |ext_part: &ExtPart| {
            match ext_part.ext_range {
                ExtRange::SummaryGenes(_) | ExtRange::SummaryTranscripts(_) | ExtRange::SummaryTerms(_) |
                ExtRange::ModifiedResidues(_) =>
                    conf_rel_ranges.contains(&ext_part.rel_type_name),
                _ =>false
            }
        };
    let mut ret_rows = vec![];

    {
        let mut row_iter = rows.iter().cloned();

        if let Some(mut prev_row) = row_iter.next() {
            for current_row in row_iter {
                if prev_row.gene_uniquenames != current_row.gene_uniquenames ||
                    prev_row.genotype_uniquenames != current_row.genotype_uniquenames {
                        ret_rows.push(prev_row);
                        prev_row = current_row;
                        continue;
                    }

                let mut prev_row_extension = prev_row.extension.clone();
                let prev_matching_ext_part_result =
                    remove_first_with_index(&mut prev_row_extension, &merge_range_rel_p);
                let mut current_row_extension = current_row.extension.clone();
                let current_matching_ext_part =
                    remove_first(&mut current_row_extension, &merge_range_rel_p);

                if let (Some((prev_ext_part, prev_removed_index)), Some(current_ext_part)) =
                    (prev_matching_ext_part_result, current_matching_ext_part) {

                        if mem::discriminant(&prev_ext_part.ext_range) ==
                            mem::discriminant(&current_ext_part.ext_range) &&
                            current_row_extension == prev_row_extension &&
                            prev_ext_part.rel_type_display_name == current_ext_part.rel_type_display_name {
                                let merged_ext_parts =
                                    merge_ext_part_ranges(&prev_ext_part,
                                                          &current_ext_part,
                                                          genes);
                                let mut new_ext = prev_row_extension.clone();
                                new_ext.insert(prev_removed_index, merged_ext_parts);
                                prev_row.extension = new_ext;
                                continue;
                            }
                    }

                ret_rows.push(prev_row);
                prev_row = current_row;
            }

            ret_rows.push(prev_row);
        }
    }

    *rows = ret_rows;
}

fn sort_genotype_uniquenames(genotypes: &IdGenotypeMap,
                             genotype_uniquenames: &mut Vec<FlexStr>) {
    let cmp_genotype_ploidiness =
        |genotype_a_uniquename: &FlexStr, genotype_b_uniquename: &FlexStr| {
            let genotype_a = genotypes.get(genotype_a_uniquename)
                .unwrap_or_else(|| panic!("missing genotype {}", genotype_a_uniquename));
            let genotype_b = genotypes.get(genotype_b_uniquename)
                .unwrap_or_else(|| panic!("missing genotype {}", genotype_b_uniquename));

            let genotype_a_ploidiness = genotype_a.ploidiness();
            let genotype_b_ploidiness = genotype_b.ploidiness();

            if genotype_a_ploidiness == genotype_b_ploidiness {
                return Ordering::Equal;
            }

            if genotype_a.ploidiness() == Ploidiness::Haploid {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        };

    genotype_uniquenames.sort_by(cmp_genotype_ploidiness);
}

// combine rows that have a gene or genotype but no extension into one row
fn collect_summary_rows(genes: &UniquenameGeneMap, genotypes: &IdGenotypeMap,
                        rows: &mut Vec<TermSummaryRow>) {
    let mut no_ext_rows = vec![];
    let mut other_rows = vec![];

    for row in rows.drain(0..) {
        if (!row.gene_uniquenames.is_empty() || !row.genotype_uniquenames.is_empty())
            && row.extension.is_empty() {
                if row.gene_uniquenames.len() > 1 {
                    panic!("row has more than one gene\n");
                }
                if row.genotype_uniquenames.len() > 1 {
                    panic!("row has more than one genotype\n");
                }
                no_ext_rows.push(row);
            } else {
                other_rows.push(row);
            }
    }

    let mut gene_uniquenames: Vec<FlexStr> =
        no_ext_rows.iter().filter(|row| !row.gene_uniquenames.is_empty())
        .map(|row| row.gene_uniquenames[0].clone())
        .collect();

    let gene_cmp =
        |uniquename1: &FlexStr, uniquename2: &FlexStr| {
            let gene1 = genes.get(uniquename1).unwrap();
            let gene2 = genes.get(uniquename2).unwrap();
            gene1.cmp(gene2)
        };
    gene_uniquenames.sort_by(gene_cmp);
    gene_uniquenames.dedup();

    let mut genotype_uniquenames: Vec<FlexStr> =
        no_ext_rows.iter().filter(|row| !row.genotype_uniquenames.is_empty())
        .map(|row| row.genotype_uniquenames[0].clone())
        .collect();

    sort_genotype_uniquenames(genotypes, &mut genotype_uniquenames);

    rows.clear();

    if !gene_uniquenames.is_empty() || !genotype_uniquenames.is_empty() {
        let genes_row = TermSummaryRow {
            gene_uniquenames,
            genotype_uniquenames,
            extension: vec![],
        };

        rows.push(genes_row);
    }
    rows.append(&mut other_rows);
}

// used to compare ExtParts by rel_type_display_name instead of rel_type_name

#[derive(Clone)]
struct ExtPartWrapper(ExtPart);

impl PartialEq for ExtPartWrapper {
    fn eq(&self, other: &Self) -> bool {
        self.0.rel_type_display_name == other.0.rel_type_display_name &&
            self.0.ext_range == other.0.ext_range
    }
}
impl Eq for ExtPartWrapper {
}
impl Hash for ExtPartWrapper {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.rel_type_display_name.hash(state);
        self.0.ext_range.hash(state);
    }
}

// Remove annotations from the summary where there is another more
// specific annotation.  ie. the same annotation but with extra part(s) in the
// extension.
// See: https://github.com/pombase/website/issues/185
pub fn remove_redundant_summary_rows(rows: &mut Vec<TermSummaryRow>) {
    let mut results = vec![];

    if rows.len() <= 1 {
        return;
    }

    let ext_parts_wrap = |ext_parts: &[ExtPart]| {
        ext_parts
            .iter()
            .map(|part| {
                ExtPartWrapper(part.clone())
            })
            .collect::<Vec<_>>()
    };

    let process_feature_rows =
        |results: &mut Vec<TermSummaryRow>, vec_set: &mut VecSet<ExtPartWrapper>,
         current_rows: &mut Vec<TermSummaryRow>| {

            let mut seen = HashSet::new();

            for current in current_rows.drain(0..) {
                if !vec_set.contains_proper_superset(&ext_parts_wrap(&current.extension)) {
                    if !seen.contains(&current.extension) {
                        // no extension is more specific (ie. has more parts), so include in the summary:
                        results.push(current.clone());
                        seen.insert(current.extension);
                    }
                }
            }
        };

    let mut current_rows = vec![];
    let mut vec_set = VecSet::new();

    let mut prev = rows.remove(0);

    vec_set.insert(ext_parts_wrap(&prev.extension));

    current_rows.push(prev.clone());

    for current in rows.drain(0..) {
        if current.gene_uniquenames.len() > 1 {
            panic!("remove_redundant_summary_rows() failed: num genes > 1\n");
        }

        if (&prev.gene_uniquenames, &prev.genotype_uniquenames) !=
            (&current.gene_uniquenames, &current.genotype_uniquenames) {
                process_feature_rows(&mut results, &mut vec_set, &mut current_rows);
                vec_set = VecSet::new();
                current_rows = vec![];
            }

        vec_set.insert(ext_parts_wrap(&current.extension));
        current_rows.push(current.clone());

        prev = current;
    }

    process_feature_rows(&mut results, &mut vec_set, &mut current_rows);

    *rows = results;
}

// Remove annotation from the summary if there is a more specific annotation
pub fn remove_redundant_summaries(children_by_termid: &HashMap<TermId, HashSet<TermId>>,
                                  term_annotations: &mut Vec<OntTermAnnotations>) {
    let mut term_annotations_by_termid = HashMap::new();

    for term_annotation in &*term_annotations {
        if !term_annotation.is_not {
            // NOT annotation don't count as more specific
            term_annotations_by_termid.insert(term_annotation.term.clone(),
                                              term_annotation.clone());
        }
    }

    for term_annotation in term_annotations.iter_mut() {
        if let Some(child_termids) = children_by_termid.get(&term_annotation.term) {
            if term_annotation.summary.as_ref().unwrap().is_empty() {
                let mut found_child_match = false;
                for child_termid in child_termids {
                    if term_annotations_by_termid.get(child_termid).is_some() {
                        found_child_match = true;
                    }
                }
                if found_child_match {
                    term_annotation.summary = None;
                }
            } else {
                let mut filtered_rows: Vec<TermSummaryRow> = vec![];

                for row in term_annotation.summary.as_mut().unwrap() {
                    let mut found_child_match = false;
                    for child_termid in child_termids {
                        if let Some(ref mut child_term_annotation) =
                            term_annotations_by_termid.get(child_termid) {
                                for child_row in child_term_annotation.summary.as_ref().unwrap() {
                                    if *row == *child_row {
                                        found_child_match = true;
                                        break;
                                    }
                                }
                            }
                    }
                    if !found_child_match {
                        filtered_rows.push(row.clone());
                    }
                }

                if filtered_rows.is_empty() {
                    term_annotation.summary = None;
                } else {
                    term_annotation.summary = Some(filtered_rows);
                }
            }
        }
    }
}


// turns binds([[gene1]]),binds([[gene2]]),other_rel(...) into:
// binds([[gene1, gene2]]),other_rel(...)
pub fn collect_duplicated_relations(ext: &mut Vec<ExtPart>) {
    let mut result: Vec<ExtPart> = vec![];

    {
        let mut iter = ext.iter().cloned();

        if let Some(mut prev) = iter.next() {
            for current in iter {
                if prev.rel_type_name != current.rel_type_name {
                    result.push(prev);
                    prev = current;
                    continue;
                }

                if let ExtRange::SummaryGenes(ref current_summ_genes) = current.ext_range {
                    if let ExtRange::SummaryGenes(ref mut prev_summ_genes) = prev.ext_range {
                        let mut current_genes = current_summ_genes[0].clone();
                        prev_summ_genes[0].append(& mut current_genes);

                        continue;
                    }
                }

                if let ExtRange::SummaryTranscripts(ref current_summ_transcripts) = current.ext_range {
                    if let ExtRange::SummaryTranscripts(ref mut prev_summ_transcripts) = prev.ext_range {
                        let mut current_transcripts = current_summ_transcripts[0].clone();
                        prev_summ_transcripts[0].append(& mut current_transcripts);

                        continue;
                    }
                }

                result.push(prev);
                prev = current;
            }

            result.push(prev);
        }
    }

    ext.clear();
    ext.append(&mut result);
}

fn make_cv_summary(cv_config: &CvConfig,
                   children_by_termid: &HashMap<TermId, HashSet<TermId>>,
                   include_gene: bool, include_genotype: bool,
                   term_and_annotations_vec: &mut Vec<OntTermAnnotations>,
                   genes: &UniquenameGeneMap,
                   genotypes: &UniquenameGenotypeMap,
                   annotation_details: &IdOntAnnotationDetailMap) {
    for term_and_annotations in term_and_annotations_vec.iter_mut() {
        let mut rows = vec![];

        let mut summary_sorted_annotations = term_and_annotations.annotations.clone();

        // in the summary, sort by extension type and length to fix:
        // https://github.com/pombase/website/issues/228
        let ext_comp = |id1: &i32, id2: &i32| {
            let a1: &OntAnnotationDetail =
                annotation_details.get(id1).expect("can't find OntAnnotationDetail");
            let a2: &OntAnnotationDetail =
                annotation_details.get(id2).expect("can't find OntAnnotationDetail");

            if a1.extension.is_empty() || a2.extension.is_empty() {
                // sort is stable (hopefully) so this will keep the ordering the same
                Ordering::Equal
            } else {
                if a1.extension[0].rel_type_display_name == a2.extension[0].rel_type_display_name {
                    a1.extension.len().cmp(&a2.extension.len())
                } else {
                    Ordering::Equal
                }
            }
        };
        summary_sorted_annotations.sort_by(ext_comp);

        for annotation_id in &summary_sorted_annotations {
            let annotation =
                annotation_details.get(annotation_id).expect("can't find OntAnnotationDetail");
            let gene_uniquenames =
                if include_gene && cv_config.feature_type == "gene" {
                    annotation.genes.clone()
                } else {
                    vec![]
                };

            let mut genotype_uniquenames =
                if include_genotype && cv_config.feature_type == "genotype" {
                    if let Some(ref genotype_uniquename) = annotation.genotype {
                        vec![genotype_uniquename.clone()]
                    } else {
                        vec![]
                    }
                } else {
                    vec![]
                };

            sort_genotype_uniquenames(genotypes, &mut genotype_uniquenames);

            let summary_relations_to_hide = &cv_config.summary_relations_to_hide;

            let mut summary_extension = annotation.extension.iter().cloned()
                .filter(|ext_part| {
                    !summary_relations_to_hide.contains(&ext_part.rel_type_name) &&
                        *summary_relations_to_hide != vec!["ALL"]
                })
                .map(move |mut ext_part| {
                    match ext_part.ext_range.clone() {
                        ExtRange::Gene(gene_uniquename) => {
                            let summ_genes = vec![gene_uniquename];
                            ext_part.ext_range = ExtRange::SummaryGenes(vec![summ_genes]);
                        },
                        ExtRange::Transcript(transcript_uniquename) => {
                            let summ_transcripts = vec![transcript_uniquename];
                            ext_part.ext_range = ExtRange::SummaryTranscripts(vec![summ_transcripts]);
                        },
                        ExtRange::Term(termid) => {
                            ext_part.ext_range = ExtRange::SummaryTerms(vec![termid]);
                        },
                        _ => (),
                    }
                    ext_part })
                .collect::<Vec<ExtPart>>();

            if gene_uniquenames.is_empty() &&
                genotype_uniquenames.is_empty() &&
                summary_extension.is_empty() {
                    continue;
                }

            collect_duplicated_relations(&mut summary_extension);

            let row = TermSummaryRow {
                gene_uniquenames,
                genotype_uniquenames,
                extension: summary_extension,
            };

            rows.push(row);
        }

        remove_redundant_summary_rows(&mut rows);

        term_and_annotations.summary = Some(rows);
    }

    remove_redundant_summaries(children_by_termid, term_and_annotations_vec);

    for term_and_annotations in &mut term_and_annotations_vec.iter_mut() {
        if let Some(ref mut summary) = term_and_annotations.summary {
            collect_summary_rows(genes, genotypes, summary);
            collect_ext_summary_genes(cv_config, summary, genes);
        }
    }
}
