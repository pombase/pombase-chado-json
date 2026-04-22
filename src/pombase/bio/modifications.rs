use std::io::{self, Write, BufWriter};
use std::fs::File;

use flexstr::{shared_str as flex_str, ToSharedStr};

use crate::bio::util::process_modification_ext;
use crate::data_types::{DataLookup, UniquenameGeneMap};
use crate::web::config::Config;


pub fn write_modifications(config: &Config, genes: &UniquenameGeneMap,
                           data_lookup: &dyn DataLookup, output_dir: &str)
         -> Result<(), io::Error>
{
    let load_org_taxonid =
        if let Some(load_org_taxonid) = config.load_organism_taxonid {
            load_org_taxonid
        } else {
            return Ok(())
        };

    let file_name = format!("{}/modifications.tsv", output_dir);
    let file = File::create(file_name)?;
    let mut writer = BufWriter::new(&file);

    let file_name_with_comments = format!("{}/modifications_with_comments.tsv", output_dir);
    let file_with_comments = File::create(file_name_with_comments)?;
    let mut writer_with_comments = BufWriter::new(&file_with_comments);

    let header = "#gene_systematic_id\tgene_name\tmodification_term_id\tevidence\tmodification\textension\treference\ttaxon_id\tdate\tassigned_by";
    writeln!(writer, "{}", header)?;
    writeln!(writer_with_comments, "{}\tcomments", header)?;

    for gene_details in genes.values() {
        if gene_details.taxonid != load_org_taxonid {
            continue;
        }

        if let Some(term_annotations) = gene_details.cv_annotations.get(&flex_str!("PSI-MOD")) {
            for term_annotation in term_annotations {
                if term_annotation.is_not {
                    continue;
                }
                for annotation_id in &term_annotation.annotations {
                    let annotation_detail = data_lookup.get_annotation_detail(*annotation_id)
                        .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                    let gene_name = gene_details.name.as_deref().unwrap_or("");
                    let mut maybe_evidence = annotation_detail.evidence.clone();
                    if let Some(ref evidence) = maybe_evidence
                        && let Some(ev_config) = config.evidence_types.get(evidence) {
                            maybe_evidence = Some(ev_config.long.to_shared_str());
                        }
                    let (modification, extension) =
                        process_modification_ext(config, data_lookup, &gene_details.uniquename,
                                                 &annotation_detail.extension);

                    let reference =
                        annotation_detail.reference.as_deref().unwrap_or("");
                    let date = annotation_detail.date.as_deref().unwrap_or("");
                    let assigned_by = annotation_detail.assigned_by.as_deref().unwrap_or("");
                    let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                       gene_details.uniquename,
                                       gene_name,
                                       term_annotation.term,
                                       maybe_evidence.unwrap_or_else(|| "".into()),
                                       modification,
                                       extension,
                                       reference,
                                       load_org_taxonid,
                                       date,
                                       assigned_by);

                    writeln!(writer, "{}", line)?;

                    if let Some(ref comment) = annotation_detail.submitter_comment {
                        writeln!(writer_with_comments, "{}\t{}", line, comment)?;
                    } else {
                        writeln!(writer_with_comments, "{}", line)?;
                    }
                }
            }
        }
    }

    Ok(())
}
