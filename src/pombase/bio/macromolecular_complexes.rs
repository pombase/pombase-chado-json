use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::{self, Write, BufWriter};
use std::fs::File;

use flexstr::SharedStr as FlexStr;

use crate::data_types::ProteinComplexData;

use crate::utils::join;

pub fn write_macromolecular_complexes(complex_data: &ProteinComplexData,
                                      output_dir: &str)
    -> Result<(), io::Error>
{
    let complexes_file_name = format!("{}/Complex_annotation.tsv", output_dir);
    let complexes_file = File::create(complexes_file_name).expect("Unable to open file");
    let mut complexes_writer = BufWriter::new(&complexes_file);

    let header = "complex_term_id\tGO_term_name\tsystematic_id\tsymbol\tgene_product_description\tUniProt_ID\tevidence_code\tsource\tassigned_by";

    complexes_writer.write_all(header.as_bytes())?;
    complexes_writer.write_all(b"\n")?;

    let mut lines = vec![];

    for (termid, term_details) in complex_data.iter() {

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

                let line_bits = vec![termid.as_str(),
                                     term_details.term_name.as_str(),
                                     gene_details.gene_short.uniquename.as_str(),
                                     gene_display_name,
                                     gene_details.gene_short.product.as_ref().map(FlexStr::as_str).unwrap_or_default(),
                                     gene_details.uniprot_identifier.as_ref().map(FlexStr::as_str).unwrap_or_default(),
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
