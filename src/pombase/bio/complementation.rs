use std::{fs::File, io::{self, BufWriter, Write}};

use flexstr::shared_str as flex_str;
use itertools::Itertools;

use crate::data_types::{DataLookup, UniquenameGeneMap};

pub fn write_complementation(data_lookup: &dyn DataLookup,
                             genes: &UniquenameGeneMap,
                             output_dir: &str)
    -> Result<(), io::Error>
{

    let complemented_by_file_name = format!("{}/complemented_by_annotation.tsv", output_dir);
    let complemented_by_file = File::create(complemented_by_file_name)
        .expect("Unable to open file");
    let mut complemented_by_writer = BufWriter::new(complemented_by_file);

    let complements_file_name = format!("{}/complements_annotation.tsv", output_dir);
    let complements_file = File::create(complements_file_name)
        .expect("Unable to open file");
    let mut complements_writer = BufWriter::new(complements_file);

    let not_complemented_by_file_name = format!("{}/not_complemented_by_annotation.tsv", output_dir);
    let not_complemented_by_file = File::create(not_complemented_by_file_name)
        .expect("Unable to open file");
    let mut not_complemented_by_writer = BufWriter::new(not_complemented_by_file);

    let does_not_complement_file_name = format!("{}/does_not_complement_annotation.tsv", output_dir);
    let does_not_complement_file = File::create(does_not_complement_file_name)
        .expect("Unable to open file");
    let mut does_not_complement_writer = BufWriter::new(does_not_complement_file);

    let empty_string = flex_str!("");

    for gene_details in genes.values() {
        let Some(term_annotations) = gene_details
            .cv_annotations.get("complementation")
        else {
            continue;
        };

        for term_annotation in term_annotations {
            let term =
                data_lookup.get_term(&term_annotation.term)
                .unwrap_or_else(|| panic!("failed to find term summary for {}",
                                          term_annotation.term));

            for annotation_id in &term_annotation.annotations {
                let annotation_detail = data_lookup.get_annotation_detail(*annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                let gene_name = gene_details.name
                    .as_ref().unwrap_or_else(|| &empty_string);
                let reference = annotation_detail.reference
                    .as_ref().unwrap_or_else(|| &empty_string);

                let qualifiers = annotation_detail.qualifiers.iter().join(",");

                let writer =
                    if term.name.starts_with("functionally complemented by") {
                        &mut complemented_by_writer
                    } else {
                        if term.name.starts_with("functionally complements") {
                            &mut complements_writer
                        } else {
                            if term.name.starts_with("is not functionally complemented by") {
                                &mut not_complemented_by_writer
                            } else {
                                if term.name.starts_with("not functionally complemented") {
                                    &mut does_not_complement_writer
                                } else {
                                    eprintln!("can't parse complementation term name {} {}: {}",
                                              gene_details.uniquename, reference, term.name);
                                    continue;
                                }
                            }
                        }
                    };
                writer.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\n",
                                              gene_details.uniquename,
                                              gene_name,
                                              term.name, qualifiers, reference))?;
            }
        }
    }

    Ok(())
}
