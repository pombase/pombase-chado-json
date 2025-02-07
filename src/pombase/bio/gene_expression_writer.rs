use std::io;
use std::io::Write;

use flexstr::shared_str as flex_str;

use crate::{data_types::{ExtPart, ExtRange, GeneDetails, OntAnnotationDetail, TermIdDetailsMap}, types::TermId};

fn extension_to_string(extension: &[ExtPart]) -> String {
    let mut part_strings = vec![];

    for ext_part in extension {
        let ExtRange::Term(ref range_string) = ext_part.ext_range
        else {
           panic!("unexpect range in gene expression extension: {}", ext_part.ext_range);
        };

        let part_string = format!("{}({})", ext_part.rel_type_display_name, range_string);

        part_strings.push(part_string);
    }

    part_strings.join(",")
}

pub fn write_quantitative_expression_row(writer: &mut dyn Write,
                                         terms: &TermIdDetailsMap,
                                         gene_details: &GeneDetails,
                                         annotation_termid: &TermId,
                                         annotation_detail: &OntAnnotationDetail)
   -> Result<(), io::Error>
{
    let empty_string = flex_str!("");

    let gene_uniquename = &gene_details.uniquename;
    let gene_name = gene_details.name.as_ref().unwrap_or(&empty_string);

    let Some(annotation_term) = terms.get(annotation_termid)
    else {
        panic!("can't find term details for: {}", annotation_termid);
    };

    let annotation_type =
        if annotation_term.name == "RNA level" {
            "RNA"
        } else {
            if annotation_term.name == "protein level" {
                "protein"
            } else {
                panic!("unexpected gene expression term name: {}",
                       annotation_term.name);
            }
        };

    let extension_string = extension_to_string(&annotation_detail.extension);
    let reference = annotation_detail.reference.as_ref().unwrap_or(&empty_string);

    let Some(ref gene_ex_props) = &annotation_detail.gene_ex_props
    else {
        if reference.len() == 0 {
            panic!("no gene_ex_props for gene expression");
        } else {
            panic!("no gene_ex_props for gene expression for reference {}", reference);
        }
    };

    let copies_per_cell =
        if let Some(ref copies_per_cell) = gene_ex_props.copies_per_cell {
            if copies_per_cell == "ND" {
                &empty_string
            } else {
                copies_per_cell
            }
        } else {
            &empty_string
        };
    let avg_copies_per_cell = gene_ex_props.avg_copies_per_cell.as_ref().unwrap_or(&empty_string);
    let scale = &gene_ex_props.scale;

    let evidence = annotation_detail.evidence.as_ref().unwrap_or(&empty_string);

    let date = annotation_detail.date.as_ref().unwrap_or(&empty_string);

    let mut condition_names =
        annotation_detail.conditions.iter().map(|termid| {
            terms.get(termid).expect(&format!("failed to find term details for {}", termid)).name.clone()
        })
        .collect::<Vec<_>>();

    condition_names.sort();

    let conditions = condition_names.join(",");

    let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                       gene_uniquename, gene_name, annotation_type,
                       extension_string, avg_copies_per_cell, copies_per_cell,
                       evidence, scale, conditions,
                       reference, gene_details.taxonid, date);

    writer.write_all(line.as_bytes())?;

    Ok(())
}

const VALID_QUALIFIERS: [&str; 8] =
    ["qualifier", "increased", "decreased", "present", "unchanged", "absent",
     "constant", "fluctuates"];

pub fn write_qualitative_expression_row(writer: &mut dyn Write,
                                        terms: &TermIdDetailsMap,
                                        gene_details: &GeneDetails,
                                        annotation_termid: &TermId,
                                        annotation_detail: &OntAnnotationDetail)
   -> Result<(), io::Error>
{
    let empty_string = flex_str!("");

    let gene_uniquename = &gene_details.uniquename;
    let gene_name = gene_details.name.as_ref().unwrap_or(&empty_string);

    let Some(annotation_term) = terms.get(annotation_termid)
    else {
        panic!("can't find term details for: {}", annotation_termid);
    };

    let Some((annotation_type, qualifier)) = annotation_term.name.rsplit_once(" ")
    else {
        panic!("can't parse term name {} {}", annotation_term.termid, annotation_term.name);
    };

    if !VALID_QUALIFIERS.contains(&qualifier) {
        eprintln!(r#"invalid qualifier "{}" from term name "{}""#,
                  qualifier, annotation_term.name);
        return Ok(());
    }

    let extension_string = extension_to_string(&annotation_detail.extension);
    let reference = annotation_detail.reference.as_ref().unwrap_or(&empty_string);

    let evidence = annotation_detail.evidence.as_ref().unwrap_or(&empty_string);

    let date = annotation_detail.date.as_ref().unwrap_or(&empty_string);

    let line =
        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                gene_uniquename, gene_name, annotation_type,
                evidence, qualifier, extension_string,
                reference, gene_details.taxonid, date);

    writer.write_all(line.as_bytes())?;

    Ok(())
}
