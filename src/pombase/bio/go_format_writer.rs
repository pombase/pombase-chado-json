use std::io;

use std::fmt::Write;

use pombase_rc_string::RcString;

use crate::web::config::*;
use crate::data_types::*;

pub const GO_ASPECT_NAMES: [&str; 3] =
    ["cellular_component", "biological_process", "molecular_function"];


fn eco_evidence_from_annotation(mapping: &GoEcoMapping,
                                annotation: &OntAnnotationDetail)
                                -> String
{
    if let Some(ref go_evidence) = annotation.evidence {
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
        String::from("")
    }
}

fn get_gpad_relation_of(term_details: &TermDetails,
                        annotation_detail: &OntAnnotationDetail) -> String {
    match term_details.cv_name.as_str() {
        "molecular_function" => {
            if annotation_detail.qualifiers.iter()
                .find(|q| q.as_str() == "contributes_to").is_some()
            {
                String::from("RO:0002326")
            } else {
                String::from("RO:0002327")
            }
        },
        "biological_process" => String::from("RO:0002331"),
        "cellular_component" => {
            if term_details.interesting_parent_details.iter()
                .any(|p| p.termid == RcString::from("GO:0032991"))
            {
                String::from("BFO:0000050")
            } else {
                String::from("RO:0001025")
            }
        },
        _ => panic!("unknown cv_name in GPAD export: {}", term_details.cv_name)
    }
}

fn make_gpad_extension_string(config: &Config, extension: &Vec<ExtPart>) -> String {
    let rel_mapping = &config.file_exports.gpad_gpi.extension_relation_mappings;
    let get_rel_termid = |ext_part: &ExtPart| {
        if let Some(map_termid) = rel_mapping.get(ext_part.rel_type_name.as_ref()) {
            map_termid.clone().unwrap()
        } else {
            ext_part.rel_type_id.clone().unwrap()
        }
    };

    let get_range = |ext_part: &ExtPart| {
        let mut range_copy = ext_part.ext_range.clone();

        if let ExtRange::Gene(ref mut gene_uniquename) = range_copy {
            if !gene_uniquename.contains(":") {
                let new_uniquename =
                    RcString::from(&format!("{}:{}", config.database_name,
                                            gene_uniquename));
                *gene_uniquename = new_uniquename;
            }
        }
        range_copy
    };

    extension.iter()
        .filter(|ext_part| {
            if let Some(map_termid) = rel_mapping.get(ext_part.rel_type_name.as_ref()) {
                map_termid.is_some()
            } else {
                ext_part.rel_type_id.is_some()
            }
        })
        .map(|ext_part| format!("{}({})", get_rel_termid(ext_part),
                                get_range(ext_part)))
        .collect::<Vec<_>>().join(",")
}


pub fn write_gene_product_annotation(gpad_writer: &mut dyn io::Write,
                                     go_eco_mappping: &GoEcoMapping, config: &Config,
                                     api_maps: &APIMaps, gene_details: &GeneDetails)
                                     -> Result<(), io::Error>
{
    let database_name = &config.database_name;
    let db_object_id = format!("{}:{}", database_name, gene_details.uniquename);

    for (cv_name, term_annotations) in &gene_details.cv_annotations {
        if !GO_ASPECT_NAMES.contains(&cv_name.as_ref()) {
            continue;
        }

        for term_annotation in term_annotations {
            let term = api_maps.terms.get(&term_annotation.term)
                .expect(&format!("failed to find term summary for {}",
                                 term_annotation.term));

            for annotation_id in &term_annotation.annotations {
                let annotation_detail = api_maps.annotation_details
                    .get(&annotation_id)
                    .expect(&format!("can't find annotation {}", annotation_id));
                let not = if term_annotation.is_not {
                    "NOT"
                } else {
                    ""
                };

                let relation = get_gpad_relation_of(term, annotation_detail);

                let ontology_class_id = &term_annotation.term;
                let reference_uniquename =
                    annotation_detail.reference.clone()
                    .unwrap_or_else(|| RcString::from(""));
                let evidence_type =
                    eco_evidence_from_annotation(go_eco_mappping, annotation_detail);
                let assigned_by = &config.database_name;

                let with_iter = annotation_detail.withs.iter();
                let from_iter = annotation_detail.froms.iter();
                let mut with_or_from_parts =
                    with_iter.chain(from_iter).map(|s| {
                        let s: RcString = s.clone().into();
                        if s.contains(":") {
                            s
                        } else {
                            RcString::from(&format!("{}:{}", assigned_by, s))
                        }
                    })
                    .collect::<Vec<RcString>>();
                with_or_from_parts.sort_unstable();
                let with_or_from = with_or_from_parts.join(",");

                if let Some(ref date) = annotation_detail.date {
                    let annotation_extensions =
                        make_gpad_extension_string(config, &annotation_detail.extension);
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


fn make_gaf_extension_string(config: &Config, extension: &Vec<ExtPart>) -> String {
    let mut ret_string = String::new();

    for (idx, ext_part) in extension.iter().enumerate() {
        ret_string += &ext_part.rel_type_name;
        ret_string.push('(');

        if let ExtRange::Gene(ref gene_uniquename) = ext_part.ext_range {
            if gene_uniquename.contains(":") {
                ret_string += gene_uniquename;
            } else {
                ret_string += &config.database_name;
                ret_string.push(':');
                ret_string += &gene_uniquename;
            }
        } else {
            write!(ret_string, "{}", &ext_part.ext_range).unwrap();
        }
        ret_string.push(')');

        if idx < extension.len() - 1 {
            ret_string.push(',');
        }
    }

    ret_string
}


pub fn write_go_annotation_format(writer: &mut dyn io::Write, config: &Config,
                                  api_maps: &APIMaps, gene_details: &GeneDetails,
                                  cv_name: &str)
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
    let db_object_synonyms =
        gene_details.synonyms.iter().filter(|synonym| {
            synonym.synonym_type == "exact"
        })
        .map(|synonym| synonym.name.to_string())
        .collect::<Vec<String>>()
        .join("|");

    let single_letter_aspect =
        if cv_name.starts_with("b") {
            "P"
        } else {
            if cv_name.starts_with("c") {
                "C"
            } else {
                "F"
            }
        };

    if let Some(term_annotations) = gene_details.cv_annotations.get(cv_name) {
        for term_annotation in term_annotations {
            for annotation_id in &term_annotation.annotations {
                let annotation_detail = api_maps.annotation_details
                    .get(&annotation_id)
                    .expect(&format!("can't find annotation {}", annotation_id));

                let mut qualifier_parts = vec![];
                if term_annotation.is_not {
                    qualifier_parts.push("NOT");
                };

                qualifier_parts.extend(annotation_detail.qualifiers.iter()
                                       .map(|s| s.as_str()));

                let go_id = &term_annotation.term;
                let reference_uniquename =
                    annotation_detail.reference.clone()
                    .unwrap_or_else(|| RcString::from(""));
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
                        let s: RcString = s.clone().into();
                        if s.contains(":") {
                            s
                        } else {
                            RcString::from(&format!("{}:{}", assigned_by, s))
                        }
                    })
                    .collect::<Vec<RcString>>();
                with_or_from_parts.sort_unstable();
                let with_or_from = with_or_from_parts.join(",");

                let evidence_code =
                    annotation_detail.evidence
                    .as_ref()
                    .map(|s| s.as_str())
                    .unwrap_or_else(|| "");

                let db_object_name =
                    if let Some(ref product) = gene_details.product {
                        product.as_str()
                    } else {
                        ""
                    };

                let transcript_type_str =
                    gene_details.transcripts[0].transcript_type.as_str();
                let db_object_type =
                    if transcript_type_str == "mRNA" {
                        "protein"
                    } else {
                        transcript_type_str
                    };
                let date =
                    if let Some(ref raw_date) = annotation_detail.date {
                        raw_date.replace("-", "")
                    } else {
                        println!("WARNING while writing GPAD file: \
                                  date missing from annotation {}", annotation_id);
                        continue;
                    };

                let annotation_extensions =
                    make_gaf_extension_string(config, &annotation_detail.extension);

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
                                   qualifier_parts.join("|"),
                                   go_id,
                                   reference_uniquename,
                                   evidence_code,
                                   with_or_from,
                                   single_letter_aspect,
                                   db_object_name,
                                   db_object_synonyms,
                                   db_object_type,
                                   gene_details.taxonid,
                                   date,
                                   assigned_by,
                                   annotation_extensions,
                                   gene_product_form_id);
                writer.write_all(line.as_bytes())?;
            }
        }
    }

    Ok(())
}

