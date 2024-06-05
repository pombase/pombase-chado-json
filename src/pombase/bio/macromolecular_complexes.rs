use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{self, Write, BufWriter};
use std::fs::File;

use flexstr::{SharedStr as FlexStr, shared_str as flex_str};

use crate::data_types::{OntAnnotation, ProteinComplexData, ProteinComplexTerm,
                        ProteinComplexGene};

use crate::web::config::Config;
use crate::utils::join;

pub fn macromolecular_complex_data(ont_annotations: &Vec<OntAnnotation>, config: &Config)
    -> ProteinComplexData
{
    let mut complex_data = HashMap::new();

    let no_evidence = flex_str!("NO_EVIDENCE");

    let annotation_details = |annotation: &OntAnnotation| {
        let evidence = annotation.evidence.clone().unwrap_or_else(|| no_evidence.clone());
        (annotation.term_short.clone(), annotation.genes.iter().next().unwrap().clone(),
         evidence)
    };

    if let Some(ref complexes_config) = config.file_exports.macromolecular_complexes {
        let check_parent_term = |el: &FlexStr| {
            *el == complexes_config.parent_complex_termid
        };
        'TERM: for annotation in ont_annotations {
            let term_short = &annotation.term_short;
            let termid = &term_short.termid;

            let reference_uniquename =
                if let Some(ref reference_short) = annotation.reference_short {
                    Some(reference_short.uniquename.clone())
                } else {
                    None
                };

            if complexes_config.excluded_terms.contains(termid) {
                continue 'TERM;
            }
            if !term_short.interesting_parent_ids.iter().any(check_parent_term) {
                continue 'TERM;
            }

            let (term_short, gene_short, evidence) = annotation_details(annotation);

            complex_data.entry(term_short.termid.clone())
                .or_insert_with(|| ProteinComplexTerm {
                    term_short,
                    complex_genes: BTreeMap::new(),
                })
                .complex_genes
                .entry(gene_short.uniquename.clone())
                .or_insert_with(|| ProteinComplexGene {
                    gene_short,
                    annotation_details: HashSet::new(),
                })
                .annotation_details
                .insert((reference_uniquename, annotation.assigned_by.clone(),
                          evidence));
        }
    }

    complex_data
}

pub fn write_macromolecular_complexes(complex_data: &ProteinComplexData,
                                      output_dir: &str)
    -> Result<(), io::Error>
{
    let complexes_file_name = format!("{}/Complex_annotation.tsv", output_dir);
    let complexes_file = File::create(complexes_file_name).expect("Unable to open file");
    let mut complexes_writer = BufWriter::new(&complexes_file);

    let header = "acc\tGO_name\tsystematic_id\tsymbol\tgene_product_description\tevidence_code\tsource\tassigned_by";

    complexes_writer.write_all(header.as_bytes())?;
    complexes_writer.write_all(b"\n")?;

    let mut lines = vec![];

    for (_, term_details) in complex_data.iter() {

        for (_, gene_details) in &term_details.complex_genes {
            let mut evidence_set = BTreeSet::new();
            let mut evidence_refs = HashMap::new();
            let mut evidence_assigned_by = HashMap::new();

            for (maybe_ref_uniquename, maybe_assigned_by, evidence) in &gene_details.annotation_details {
                evidence_set.insert(evidence.clone());
                if let Some(ref_uniquename) = maybe_ref_uniquename {
                    evidence_refs.entry(evidence.clone())
                    .or_insert_with(HashSet::new)
                    .insert(ref_uniquename.to_owned());
                }
                if let Some(assigned_by) = maybe_assigned_by {
                    evidence_assigned_by.entry(evidence.clone())
                    .or_insert_with(HashSet::new)
                    .insert(assigned_by.to_owned());
                }
            }

            for evidence in evidence_set.iter() {

                let mut refs_vec = evidence_refs.remove(evidence)
                    .unwrap_or_else(HashSet::new)
                    .into_iter().collect::<Vec<_>>();
                refs_vec.sort();
                let mut assigned_bys_vec = evidence_assigned_by.remove(evidence)
                    .unwrap_or_else(HashSet::new)
                    .into_iter().collect::<Vec<_>>();
                assigned_bys_vec.sort();

                let refs_string = join(&refs_vec, ",");
                let assigned_by_string = join(&assigned_bys_vec, ",");
                let gene_display_name = gene_details.gene_short.name.as_ref().map(FlexStr::as_str)
                    .unwrap_or_else(|| gene_details.gene_short.uniquename.as_str());

                let line_bits = vec![term_details.term_short.termid.as_str(),
                                     term_details.term_short.name.as_str(),
                                     gene_details.gene_short.uniquename.as_str(),
                                     gene_display_name,
                                     gene_details.gene_short.product.as_ref().map(FlexStr::as_str).unwrap_or_else(|| ""),
                                     evidence.as_str(), refs_string.as_str(),
                                     assigned_by_string.as_str()];

                lines.push(line_bits.join("\t"));
            }
        }
    }

   lines.sort();

    for line in lines.drain(0..) {
        complexes_writer.write_all(line.as_bytes())?;
        complexes_writer.write_all(b"\n")?
    }

    Ok(())
}
