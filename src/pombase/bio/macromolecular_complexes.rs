use std::collections::{HashMap, HashSet};
use std::io::{self, Write, BufWriter};
use std::fs::File;

use flexstr::{SharedStr as FlexStr, shared_str as flex_str};

use crate::data_types::{GeneShort, OntAnnotation, TermShort};
use crate::web::config::Config;
use crate::utils::join;


pub fn write_macromolecular_complexes(ont_annotations: &Vec<OntAnnotation>,
                                      config: &Config, output_dir: &str)
                                      -> Result<(), io::Error>
{
    let mut complex_data: HashMap<(TermShort, GeneShort, FlexStr), _> = HashMap::new();

    let no_evidence = flex_str!("NO_EVIDENCE");

    let make_key = |annotation: &OntAnnotation| {
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

            if complexes_config.excluded_terms.contains(termid) {
                continue 'TERM;
            }
            if !term_short.interesting_parent_ids.iter().any(check_parent_term) {
                continue 'TERM;
            }

            let key: (TermShort, GeneShort, FlexStr) = make_key(annotation);
            complex_data.entry(key)
                .or_insert_with(Vec::new)
                .push((annotation.reference_short.clone(), annotation.assigned_by.clone()));
        }
    }

    let complexes_file_name = format!("{}/Complex_annotation.tsv", output_dir);
    let complexes_file = File::create(complexes_file_name).expect("Unable to open file");
    let mut complexes_writer = BufWriter::new(&complexes_file);

    let header = "acc\tGO_name\tsystematic_id\tsymbol\tgene_product_description\tevidence_code\tsource\tassigned_by";

    complexes_writer.write_all(header.as_bytes())?;
    complexes_writer.write_all(b"\n")?;

    let mut lines = vec![];

    for (key, values) in complex_data.drain() {
        let (term_short, gene_short, evidence) = key;
        let mut refs = HashSet::new();
        let mut assigned_bys = HashSet::new();
        for (maybe_ref_short, maybe_assigned_by) in values {
            if let Some(ref_short) = maybe_ref_short {
                refs.insert(ref_short.uniquename);
            }
            if let Some(assigned_by) = maybe_assigned_by {
                assigned_bys.insert(assigned_by);
            }
        }

        let mut refs_vec = refs.into_iter().collect::<Vec<_>>();
        refs_vec.sort();
        let mut assigned_bys_vec = assigned_bys.into_iter().collect::<Vec<_>>();
        assigned_bys_vec.sort();

        let refs_string = join(&refs_vec, ",");
        let assigned_by_string = join(&assigned_bys_vec, ",");

        let line_bits = vec![term_short.termid.as_str(), term_short.name.as_str(),
                             gene_short.uniquename.as_str(),
                             gene_short.name.as_ref().map(FlexStr::as_str)
                             .unwrap_or_else(|| gene_short.uniquename.as_str()),
                             gene_short.product.as_ref().map(FlexStr::as_str).unwrap_or_else(|| ""),
                             evidence.as_str(), refs_string.as_str(),
                             assigned_by_string.as_str()];

        lines.push(line_bits.join("\t"));
    }

    lines.sort();

    for line in lines.drain(0..) {
        complexes_writer.write_all(line.as_bytes())?;
        complexes_writer.write_all(b"\n")?
    }

    Ok(())
}
