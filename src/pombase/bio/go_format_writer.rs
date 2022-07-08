use std::io::{self, BufWriter};

use std::cmp::Ordering;
use std::collections::{HashSet, HashMap};

use std::io::Write;
use std::fs::File;

use chrono::prelude::{Local, DateTime};

use flexstr::{SharedStr as FlexStr, shared_str as flex_str, shared_fmt as flex_fmt};

use crate::utils::join;
use crate::web::config::*;
use crate::data_types::*;
use crate::types::TermId;

use super::util::make_extension_string;

#[derive(Clone, PartialEq, Eq)]
pub enum GpadGafWriteMode {
  PomBaseGaf,
  StandardGaf,
  Gpad,
}


pub const GO_ASPECT_NAMES: [FlexStr; 3] =
    [flex_str!("cellular_component"), flex_str!("biological_process"),
     flex_str!("molecular_function")];


pub fn write_go_annotation_files(api_maps: &APIMaps, config: &Config,
                                 db_creation_datetime: &FlexStr,
                                 go_eco_mappping: &GoEcoMapping,
                                 output_dir: &str)
        -> Result<(), io::Error>
{
    let load_org_taxonid =
        if let Some(load_org_taxonid) = config.load_organism_taxonid {
            load_org_taxonid
        } else {
            return Ok(())
        };

    let database_name = &config.database_name;

    let gpi_file_name =
        format!("{}/gene_product_information_taxonid_{}.tsv", output_dir,
                load_org_taxonid);
    let gpad_file_name =
        format!("{}/gene_product_annotation_data_taxonid_{}.tsv", output_dir,
                load_org_taxonid);
    let pombase_gaf_file_name = format!("{}/pombase_style_gaf.tsv", output_dir);
    let standard_gaf_file_name = format!("{}/go_style_gaf.tsv", output_dir);

    let gpi_file = File::create(gpi_file_name).expect("Unable to open file");
    let gpad_file = File::create(gpad_file_name).expect("Unable to open file");
    let pombase_gaf_file =
        File::create(pombase_gaf_file_name).expect("Unable to open file");
    let standard_gaf_file =
        File::create(standard_gaf_file_name).expect("Unable to open file");
    let mut gpi_writer = BufWriter::new(&gpi_file);
    let mut gpad_writer = BufWriter::new(&gpad_file);
    let mut pombase_gaf_writer = BufWriter::new(&pombase_gaf_file);
    let mut standard_gaf_writer = BufWriter::new(&standard_gaf_file);

    let generated_by = format!("!generated-by: {}\n", database_name);
    let iso_date = db_creation_datetime.replace(" ", "T");
    let date_generated = format!("!date-generated: {}\n", &iso_date);
    let url_header = format!("!URL: {}\n", &config.base_url);
    let funding_header = format!("!funding: {}\n", &config.funder);

    gpi_writer.write_all("!gpi-version: 2.0\n".as_bytes())?;
    gpi_writer.write_all(format!("!namespace: {}\n", database_name).as_bytes())?;
    gpi_writer.write_all(generated_by.as_bytes())?;
    gpi_writer.write_all(date_generated.as_bytes())?;
    gpi_writer.write_all(url_header.as_bytes())?;
    gpi_writer.write_all(funding_header.as_bytes())?;

    gpad_writer.write_all("!gpa-version: 2.0\n".as_bytes())?;
    gpad_writer.write_all(generated_by.as_bytes())?;
    gpad_writer.write_all(date_generated.as_bytes())?;
    gpad_writer.write_all(url_header.as_bytes())?;
    gpad_writer.write_all(funding_header.as_bytes())?;

    standard_gaf_writer.write_all("!gaf-version: 2.2\n".as_bytes())?;
    standard_gaf_writer.write_all(generated_by.as_bytes())?;
    standard_gaf_writer.write_all(date_generated.as_bytes())?;
    standard_gaf_writer.write_all(url_header.as_bytes())?;
    let contact = format!("!contact: {}\n", &config.helpdesk_address);
    standard_gaf_writer.write_all(contact.as_bytes())?;

    for gene_details in api_maps.genes.values() {
        if gene_details.taxonid != load_org_taxonid {
            continue;
        }

        if gene_details.feature_type == "ncRNA gene" {
            let mut found = false;
            for aspect in GO_ASPECT_NAMES.iter() {
                let term_annotations = gene_details.cv_annotations.get(aspect);
                if term_annotations.is_some() {
                    found = true;
                    break;
                }
            }
            if !found {
                // special case for ncRNA genes: only write ND lines for an aspect
                // if there are some annotations for the other aspects
                continue;
            }
        }

        write_gene_to_gpi(&mut gpi_writer, config, api_maps, gene_details)?;

        write_gene_product_annotation(&mut gpad_writer, go_eco_mappping, config,
                                      api_maps, gene_details)?;

        for aspect_name in &GO_ASPECT_NAMES {
            write_go_annotation_format(&mut pombase_gaf_writer, config,
                                       GpadGafWriteMode::PomBaseGaf,
                                       api_maps, gene_details,
                                       aspect_name)?;
            write_go_annotation_format(&mut standard_gaf_writer, config,
                                       GpadGafWriteMode::StandardGaf,
                                       api_maps, gene_details,
                                       aspect_name)?;
        }
    }

    Ok(())
}

fn eco_evidence_from_annotation(mapping: &GoEcoMapping,
                                annotation: &OntAnnotationDetail)
                                -> FlexStr
{
    if let Some(ref go_evidence) = annotation.
        evidence {
        if let Some(ref reference) = annotation.reference {
            mapping.lookup_with_go_ref(go_evidence, reference)
                .unwrap_or_else(|| mapping.lookup_default(go_evidence)
                                .unwrap_or_else(|| panic!("failed to find {} for {} in GO mapping",
                                                          go_evidence, reference)))
        } else {
            mapping.lookup_default(go_evidence)
                .unwrap_or_else(|| panic!("failed to find {} in GO mapping",
                                          go_evidence))
        }
    } else {
        flex_str!("")
    }
}

fn get_gpad_relation_of(term_details: &TermDetails,
                        annotation_detail: &OntAnnotationDetail) -> FlexStr {
    match term_details.cv_name.as_str() {
        "molecular_function" => {
            if annotation_detail.qualifiers.iter()
                .any(|q| q.as_str() == "contributes_to")
            {
                flex_str!("RO:0002326")
            } else {
                flex_str!("RO:0002327")
            }
        },
        "biological_process" => flex_str!("RO:0002331"),
        "cellular_component" => {
            if term_details.termid == flex_str!("GO:0032991") ||
                term_details.interesting_parent_details.iter()
                .any(|p| p.termid == flex_str!("GO:0032991"))
            {
                flex_str!("BFO:0000050")
            } else {
                flex_str!("RO:0002432")
            }
        },
        _ => panic!("unknown cv_name in GPAD export: {}", term_details.cv_name)
    }
}

fn get_gpad_relation_name_of(term_map: &HashMap<TermId, TermDetails>,
                             term_details: &TermDetails,
                             annotation_detail: &OntAnnotationDetail) -> FlexStr {
    let rel_id = get_gpad_relation_of(term_details, annotation_detail);
    term_map.get(&rel_id)
            .unwrap_or_else(|| panic!("internal error, can't find term {}", rel_id))
            .name.clone()
}

fn get_gpad_nd_relation_of(aspect: &FlexStr) -> FlexStr {
    match aspect.as_ref() {
        "molecular_function" => flex_str!("RO:0002327"),
        "biological_process" => flex_str!("RO:0002331"),
        "cellular_component" => flex_str!("RO:0002432"),
        _ => panic!("unknown aspect in GPAD export: {}", aspect)
    }
}

fn get_gpad_nd_relation_name_of(term_map: &HashMap<TermId, TermDetails>,
                                aspect: &FlexStr) -> FlexStr {
    let rel_id = get_gpad_nd_relation_of(aspect);
    term_map.get(&rel_id)
            .unwrap_or_else(|| panic!("internal error, can't find term {}", rel_id))
            .name.clone()
}

fn compare_withs(withs1: &HashSet<WithFromValue>,
                 withs2: &HashSet<WithFromValue>)
                 -> Ordering
{
    let mut withs1_vec: Vec<_> = withs1.iter().collect();
    let mut withs2_vec: Vec<_> = withs2.iter().collect();

    withs1_vec.sort();
    withs2_vec.sort();

    withs1_vec.cmp(&withs2_vec)
}

pub fn write_gene_to_gpi(gpi_writer: &mut dyn Write, config: &Config, api_maps: &APIMaps,
                         gene_details: &GeneDetails)
                         -> Result<(), io::Error>
{
    let database_name = &config.database_name;

    let db_object_id = format!("{}:{}", database_name, gene_details.uniquename);
    let db_object_symbol =
        gene_details.product.clone().unwrap_or(flex_str!(""));
    let db_object_name =
        gene_details.name.clone().unwrap_or(flex_str!(""));

    let mut db_object_synonyms =
        gene_details.synonyms.iter().filter(|synonym| {
            synonym.synonym_type == "exact"
        })
        .map(|synonym| synonym.name.to_string())
        .collect::<Vec<String>>();

    db_object_synonyms.sort();

    let db_object_synonyms_string = db_object_synonyms.join("|");

    let db_object_type = config.file_exports.gpad_gpi.transcript_gene_so_term_map
        .get(&gene_details.transcript_so_termid)
        .unwrap_or_else(|| {
            panic!("failed for find configuration for {} in transcript_gene_so_term_map",
                   &gene_details.transcript_so_termid);
        });

    let db_object_taxon =
        format!("NCBITaxon:{}",
                config.load_organism_taxonid.expect("internal error, no load_organism_taxonid"));
    let db_xrefs =
        if let Some(ref uniprot_id) = gene_details.uniprot_identifier {
            flex_fmt!("UniProtKB:{}", uniprot_id)
        } else {
            flex_str!("")
        };

    let gpi_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t\t\t\t{}\tgo-annotation-summary={}\n",
                           db_object_id,
                           db_object_name,
                           db_object_name,
                           db_object_synonyms_string,
                           db_object_type,
                           db_object_taxon,
                           db_xrefs,
                           db_object_symbol);
    gpi_writer.write_all(gpi_line.as_bytes())?;

    let db_protein_id =
        if let Some(transcript_uniquename) = gene_details.transcripts.get(0) {
            if let Some(maybe_transcript_details) =
                gene_details.transcripts_by_uniquename.get(transcript_uniquename)
            {
                if let Some(transcript_details) = maybe_transcript_details {
                    if let Some(ref protein) = transcript_details.protein {
                        format!("{}:{}", database_name, protein.uniquename)
                    } else {
                        return Ok(())
                    }
                } else {
                        return Ok(())
                }
            } else {
                return Ok(())
            }
        } else {
            return Ok(())
        };

    let mut pr_ids_seen = HashSet::new();

    for aspect in GO_ASPECT_NAMES.iter() {
        let term_annotations =
        match gene_details.cv_annotations.get(aspect) {
            Some(term_annotations) => term_annotations.clone(),
            None => continue,
        };

        for term_annotation in term_annotations {
            for annotation_id in &term_annotation.annotations {
                let annotation_detail = api_maps.annotation_details
                    .get(annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                if let Some(ref gene_product_form_id) = annotation_detail.gene_product_form_id {
                    if pr_ids_seen.contains(gene_product_form_id) {
                        continue;
                    }

                    pr_ids_seen.insert(gene_product_form_id);

                    let gpi_line =
                        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\tgo-annotation-summary={}\n",
                                gene_product_form_id,
                                db_object_name,
                                db_object_name,
                                db_object_synonyms_string,
                                db_object_type,
                                db_object_taxon,
                                db_object_id,
                                db_protein_id,
                                db_xrefs,
                                db_object_symbol);
                    gpi_writer.write_all(gpi_line.as_bytes())?;
                }
            }
        }
    }

    Ok(())
}

pub fn write_gene_product_annotation(gpad_writer: &mut dyn io::Write,
                                     go_eco_mappping: &GoEcoMapping, config: &Config,
                                     api_maps: &APIMaps, gene_details: &GeneDetails)
                                     -> Result<(), io::Error>
{
    let database_name = &config.database_name;
    let db_object_id = flex_fmt!("{}:{}", database_name, gene_details.uniquename);
    let local: DateTime<Local> = Local::now();
    let local_iso_date = local.format("%F");
    let assigned_by = &config.database_name;

    if let Some(ref product) = gene_details.product {
        if product == "dubious" {
            return Ok(());
        }
    }

    if gene_details.feature_type == "pseudogene" {
        return Ok(());
    }

    for aspect in GO_ASPECT_NAMES.iter() {
        let term_annotations = gene_details.cv_annotations.get(aspect);
        if term_annotations.is_none() {
            let relation = get_gpad_nd_relation_of(aspect);
            let go_aspect_termid =
                config.file_exports.gpad_gpi.go_aspect_terms.get(aspect).unwrap();
            let nd_ref = &config.file_exports.nd_reference;
            let line = format!("{}\t\t{}\t{}\t{}\tECO:0000307\t\t\t{}\t{}\t\t\n",
                               db_object_id,
                               relation,
                               &go_aspect_termid,
                               nd_ref,
                               local_iso_date, assigned_by);
            gpad_writer.write_all(line.as_bytes())?;
        }
    }

    for aspect in GO_ASPECT_NAMES.iter() {
        let mut sorted_term_annotations =
            match gene_details.cv_annotations.get(aspect) {
                Some(term_annotations) => term_annotations.clone(),
                None => continue,
            };

        sorted_term_annotations.sort_by(|ta1, ta2| {
            ta1.term.cmp(&ta2.term)
        });

        for term_annotation in sorted_term_annotations {
            let term = api_maps.terms.get(&term_annotation.term)
                .unwrap_or_else(|| panic!("failed to find term summary for {}",
                                 term_annotation.term));

            let mut sorted_annotations = term_annotation.annotations.clone();

            sorted_annotations.sort_by(|a1, a2| {
                let detail1 =
                    api_maps.annotation_details.get(a1)
                    .unwrap_or_else(|| panic!("can't find annotation {}", a1));
                let detail2 =
                    api_maps.annotation_details.get(a2)
                    .unwrap_or_else(|| panic!("can't find annotation {}", a1));

                let ev_res = detail1.evidence.cmp(&detail2.evidence);
                if ev_res != Ordering::Equal {
                    return ev_res;
                }

                let ref_res = detail1.reference.cmp(&detail2.reference);
                if ref_res != Ordering::Equal {
                    return ref_res;
                }

                let date_res = detail1.date.cmp(&detail2.date);
                if date_res != Ordering::Equal {
                    return date_res
                }

                compare_withs(&detail1.withs, &detail2.withs)
            });

            for annotation_id in &sorted_annotations {
                let annotation_detail = api_maps.annotation_details
                    .get(annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));
                let not = if term_annotation.is_not {
                    "NOT"
                } else {
                    ""
                };

                let relation = get_gpad_relation_of(term, annotation_detail);

                let ontology_class_id = &term_annotation.term;
                let reference_uniquename =
                    annotation_detail.reference.clone()
                    .unwrap_or_else(|| flex_str!(""));
                let evidence_type =
                    eco_evidence_from_annotation(go_eco_mappping, annotation_detail);

                let with_iter = annotation_detail.withs.iter();
                let from_iter = annotation_detail.froms.iter();
                let mut with_or_from_parts =
                    with_iter.chain(from_iter).map(|s| {
                        let with_from_id: FlexStr = s.id();
                        if with_from_id.contains(':') {
                            with_from_id
                        } else {
                            flex_fmt!("{}:{}", assigned_by, with_from_id)
                        }
                    })
                    .collect::<Vec<FlexStr>>();
                with_or_from_parts.sort_unstable();
                let with_or_from = join(&with_or_from_parts, ",");

                if let Some(ref date) = annotation_detail.date {
                    let annotation_extensions =
                        make_extension_string(config, &api_maps.terms,
                                              &GpadGafWriteMode::Gpad,
                                              &annotation_detail.extension);
                    let db_object_id =
                        if let Some(ref gene_product_form_id) = annotation_detail.gene_product_form_id {
                          &gene_product_form_id
                        } else {
                          &db_object_id
                        };
                    let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t\n",
                                       db_object_id,
                                       not, relation, ontology_class_id,
                                       reference_uniquename, evidence_type,
                                       with_or_from, date, assigned_by,
                                       annotation_extensions);
                    gpad_writer.write_all(line.as_bytes())?;
                } else {
                    println!("WARNING while writing GAF file: \
                              date missing from annotation {}", annotation_id);
                }
            }
        }
    }

    Ok(())
}

pub fn write_go_annotation_format(writer: &mut dyn io::Write, config: &Config,
                                  write_mode: GpadGafWriteMode,
                                  api_maps: &APIMaps, gene_details: &GeneDetails,
                                  cv_name: &FlexStr)
                                  -> Result<(), io::Error>
{
    if !GO_ASPECT_NAMES.contains(&cv_name) {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
                                  format!("unknown CV: {}", cv_name)));
    }

    let database_name = &config.database_name;
    let db_object_id = &gene_details.uniquename;
    let db_object_symbol =
        if let Some(ref gene_name) = gene_details.name {
            gene_name.as_str()
        } else {
            &db_object_id
        };
    let mut db_object_synonyms =
        gene_details.synonyms.iter().filter(|synonym| {
            synonym.synonym_type == "exact"
        })
        .map(|synonym| synonym.name.to_string())
        .collect::<Vec<String>>();

    db_object_synonyms.sort();

    let db_object_synonyms_string = db_object_synonyms.join("|");

        let db_object_name =
        if let Some(ref product) = gene_details.product {
            product.as_str()
        } else {
            ""
        };

    let db_object_type =
        if let Some(transcript_uniquename) = gene_details.transcripts.get(0) {
            let transcript_details = api_maps
                .transcripts.get(transcript_uniquename)
                .expect(&format!("internal error, failed to find transcript: {}",
                                 transcript_uniquename));

            let transcript_type = transcript_details.transcript_type.as_str();

            if transcript_type == "mRNA" {
                "protein"
            } else {
                transcript_type
            }
        } else {
            return Ok(());
        };

    let single_letter_aspect =
        if cv_name.starts_with('b') {
            "P"
        } else {
            if cv_name.starts_with('c') {
                "C"
            } else {
                "F"
            }
        };

    if let Some(term_annotations) = gene_details.cv_annotations.get(cv_name) {
        for term_annotation in term_annotations {
            let go_term = api_maps.terms.get(&term_annotation.term)
                .unwrap_or_else(|| panic!("failed to find term summary for {}",
                                          term_annotation.term));

            for annotation_id in &term_annotation.annotations {
                let annotation_detail = api_maps.annotation_details
                    .get(annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                let go_id = &term_annotation.term;

                let qualifiers = if write_mode == GpadGafWriteMode::PomBaseGaf {
                    let mut qualifier_parts = vec![];
                    if term_annotation.is_not {
                        qualifier_parts.push("NOT");
                    };

                    qualifier_parts.extend(annotation_detail.qualifiers.iter()
                                   .map(|s| s.as_str()));

                    qualifier_parts.sort();

                    qualifier_parts.join("|")
                } else {
                    let relation_name =
                        get_gpad_relation_name_of(&api_maps.terms,
                                                  go_term, annotation_detail)
                                                  .to_string();

                    if term_annotation.is_not {
                        format!("NOT|{}", relation_name)
                    } else {
                        relation_name
                    }
                };

                let reference_uniquename =
                    annotation_detail.reference.clone()
                    .unwrap_or_else(|| flex_str!(""));
                let assigned_by =
                    if let Some(ref assigned_by) = annotation_detail.assigned_by {
                        assigned_by.as_str()
                    } else {
                        database_name.as_str()
                    };

                let with_iter = annotation_detail.withs.iter();
                let from_iter = annotation_detail.froms.iter();
                let mut with_or_from_parts =
                    with_iter.chain(from_iter).map(|s| {
                        let with_from_id: FlexStr = s.id();
                        if with_from_id.contains(':') {
                            with_from_id
                        } else {
                            flex_fmt!("{}:{}", database_name, with_from_id)
                        }
                    })
                    .collect::<Vec<FlexStr>>();
                with_or_from_parts.sort_unstable();
                let with_or_from = join(&with_or_from_parts,",");

                let evidence_code =
                    annotation_detail.evidence
                    .as_ref()
                    .map(|s| s.as_str())
                    .unwrap_or_else(|| "");

                let date =
                    if let Some(ref raw_date) = annotation_detail.date {
                        raw_date.replace("-", "")
                    } else {
                        println!("WARNING while writing GPAD file: \
                                  date missing from annotation {}", annotation_id);
                        continue;
                    };

                let annotation_extensions =
                    make_extension_string(config, &api_maps.terms,
                                          &write_mode, &annotation_detail.extension);

                let gene_product_form_id =
                    if let Some(ref gene_product_form_id) = annotation_detail.gene_product_form_id
                {
                    gene_product_form_id.as_str()
                } else {
                    ""
                };

                let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\ttaxon:{}\t{}\t{}\t{}\t{}\n",
                                   database_name,
                                   db_object_id,
                                   db_object_symbol,
                                   qualifiers,
                                   go_id,
                                   reference_uniquename,
                                   evidence_code,
                                   with_or_from,
                                   single_letter_aspect,
                                   db_object_name,
                                   db_object_synonyms_string,
                                   db_object_type,
                                   gene_details.taxonid,
                                   date,
                                   assigned_by,
                                   annotation_extensions,
                                   gene_product_form_id);
                writer.write_all(line.as_bytes())?;
            }
        }

    } else {
        if write_mode == GpadGafWriteMode::StandardGaf && db_object_type == "protein" {
            let local: DateTime<Local> = Local::now();
            let date = local.format("%Y%m%d");
            let relation_name =
                get_gpad_nd_relation_name_of(&api_maps.terms, cv_name);
            let go_aspect_termid =
                config.file_exports.gpad_gpi.go_aspect_terms.get(cv_name).unwrap();
            let nd_ref = &config.file_exports.nd_reference;
            let assigned_by = database_name.as_str();

            let line = format!("{}\t{}\t{}\t{}\t{}\t{}\tND\t\t{}\t{}\t{}\t{}\ttaxon:{}\t{}\t{}\t\t\n",
                               database_name,
                               db_object_id,
                               db_object_symbol,
                               relation_name,
                               go_aspect_termid,
                               nd_ref,
                               single_letter_aspect,
                               db_object_name,
                               db_object_synonyms_string,
                               db_object_type,
                               gene_details.taxonid,
                               date,
                               assigned_by);

            writer.write_all(line.as_bytes())?;
        }
    }

    Ok(())
}

