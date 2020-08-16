use std::io::Write;
use std::io;

use pombase_rc_string::RcString;

use crate::web::config::*;
use crate::data_types::*;

const GO_ASPECT_NAMES: [&str; 3] =
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

fn make_extension_string(config: &Config, extension: &Vec<ExtPart>) -> String {
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


pub fn write_gene_product_annotation(gpad_writer: &mut dyn Write,
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

                let mut with_or_from_parts =
                    annotation_detail.withs.iter().map(|s| {
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

                let date = annotation_detail.date.clone()
                    .expect(&format!("date missing from annotation with ID {}",
                                     annotation_id));
                let annotation_extensions =
                    make_extension_string(config, &annotation_detail.extension);
                let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t\n",
                                   db_object_id,
                                   not, relation, ontology_class_id,
                                   reference_uniquename, evidence_type,
                                   with_or_from, date, assigned_by,
                                   annotation_extensions);
                gpad_writer.write_all(line.as_bytes())?;
            }
        }
    }

    Ok(())
}

