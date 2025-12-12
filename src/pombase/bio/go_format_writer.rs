use std::io::{self, BufWriter};

use std::cmp::Ordering;
use std::collections::HashSet;

use std::io::Write;
use std::fs::File;
use std::sync::Arc;

use arrow::datatypes::Date32Type;
use chrono::NaiveDate;

use arrow_array::{Array, Date32Array, LargeStringArray, RecordBatch};
use arrow_schema::{DataType, Field, Schema};
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;
use parquet::basic::{Compression, ZstdLevel};

use chrono::prelude::{Local, DateTime};

use flexstr::{SharedStr as FlexStr, shared_str as flex_str, shared_fmt as flex_fmt};
use regex::Regex;

use crate::bio::{get_submitter_comment, ExportComments};
use crate::utils::join;
use crate::web::config::*;
use crate::data_types::*;

use crate::bio::util::make_extension_string;

#[derive(Clone, PartialEq, Eq)]
pub enum GpadGafWriteMode {
  PomBaseGaf,
  ExtendedPomBaseGaf,
  StandardGaf,
  Gpad,
}

pub const GO_ASPECT_NAMES: [FlexStr; 3] =
    [flex_str!("cellular_component"), flex_str!("biological_process"),
     flex_str!("molecular_function")];

fn abbreviation_of_go_aspect(cv_name: &str)
    -> &'static str
{
    if cv_name.starts_with('b') {
        "P"
    } else if cv_name.starts_with('c') {
        "C"
    } else {
        "F"
    }
}

fn write_tsv_lines(writer: &mut dyn io::Write, lines: &Vec<Vec<String>>)
  -> Result<(), io::Error>
{
    for line_parts in lines {
        let line = line_parts.join("\t");
        writeln!(writer, "{}", line)?;
    }

    Ok(())
}

fn write_parquet(writer: &mut BufWriter<&File>,
                 header_parts: &[&str],
                 lines: &[Vec<String>])
  -> Result<(), io::Error>
{
    let mut gaf_schema_fields = vec![];

    let mut date_idx = 0;

    for (idx, header_part) in header_parts.iter().enumerate() {
        let field = if *header_part == "date" {
            date_idx = idx;
            Field::new(header_part.to_owned(), DataType::Date32, false)
        } else {
            Field::new(header_part.to_owned(), DataType::LargeUtf8, false)
        };
        gaf_schema_fields.push(field);
    }
    let schema = Arc::new(Schema::new(gaf_schema_fields));

    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::try_new(16).unwrap()))
        .build();

    let mut arrow_writer =
        ArrowWriter::try_new(writer, schema.clone(), Some(props))?;

    for chunk in lines.chunks(1000) {
        let chunk_lines: Vec<_> = chunk.iter().collect();

        let mut columns: Vec<Vec<String>> = vec![];
        for _ in 0..chunk_lines[0].len() {
            columns.push(vec![]);
        }

        for chunk_line in chunk_lines.into_iter() {
            for (idx, value) in chunk_line.iter().enumerate() {
                columns[idx].push(value.to_owned());
            }
        }

        let str_vec_vec: Vec<Arc<dyn Array>> =
            columns.iter().enumerate()
            .map(|(idx, v)| {
                if idx == date_idx {
                    let date_iter = v.iter()
                        .map(|d| {
                            let naive_date = NaiveDate::parse_from_str(d, "%Y%m%d").unwrap();
                            Date32Type::from_naive_date(naive_date)
                        });
                    Arc::new(Date32Array::from_iter_values(date_iter)) as _
                } else {
                    Arc::new(LargeStringArray::from_iter_values(v)) as _
                }
            })
            .collect();

        let record_batch = RecordBatch::try_new(
            schema.clone(),
            str_vec_vec,
        ).unwrap();

        arrow_writer.write(&record_batch)?;
    }

    arrow_writer.close()?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn write_go_annotation_files(api_maps: &APIMaps, config: &Config,
                                 data_lookup: &dyn DataLookup,
                                 db_creation_datetime: &FlexStr,
                                 go_eco_mappping: &GoEcoMapping,
                                 genes: &UniquenameGeneMap,
                                 transcripts: &UniquenameTranscriptMap,
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
    let pombase_gaf_parquet_file_name = format!("{}/pombase_style_gaf.parquet", output_dir);
    let extended_pombase_gaf_file_name = format!("{}/extended_pombase_style_gaf.tsv", output_dir);
    let extended_pombase_gaf_parquet_file_name = format!("{}/extended_pombase_style_gaf.parquet", output_dir);
    let standard_gaf_file_name = format!("{}/go_style_gaf.tsv", output_dir);
    let standard_gaf_parquet_file_name = format!("{}/go_style_gaf.parquet", output_dir);
    let comments_gaf_file_name = format!("{}/canto_go_annotations_with_comments.tsv", output_dir);
    let comments_gaf_parquet_file_name = format!("{}/canto_go_annotations_with_comments.parquet", output_dir);

    let gpi_file = File::create(gpi_file_name).expect("Unable to open file");
    let gpad_file = File::create(gpad_file_name).expect("Unable to open file");
    let pombase_gaf_file =
        File::create(pombase_gaf_file_name).expect("Unable to open file");
    let pombase_gaf_parquet_file =
        File::create(pombase_gaf_parquet_file_name).expect("Unable to open file");
    let extended_pombase_gaf_file =
        File::create(extended_pombase_gaf_file_name).expect("Unable to open file");
    let extended_pombase_gaf_parquet_file =
        File::create(extended_pombase_gaf_parquet_file_name).expect("Unable to open file");
    let standard_gaf_file =
        File::create(standard_gaf_file_name).expect("Unable to open file");
    let standard_gaf_parquet_file =
        File::create(standard_gaf_parquet_file_name).expect("Unable to open file");
    let comments_gaf_file =
        File::create(comments_gaf_file_name).expect("Unable to open file");
    let comments_gaf_parquet_file =
        File::create(comments_gaf_parquet_file_name).expect("Unable to open file");

    let mut gpi_writer = BufWriter::new(&gpi_file);
    let mut gpad_writer = BufWriter::new(&gpad_file);
    let mut extended_pombase_gaf_writer = BufWriter::new(&extended_pombase_gaf_file);
    let mut extended_pombase_gaf_parquet_writer = BufWriter::new(&extended_pombase_gaf_parquet_file);
    let mut pombase_gaf_writer = BufWriter::new(&pombase_gaf_file);
    let mut pombase_gaf_parquet_writer = BufWriter::new(&pombase_gaf_parquet_file);
    let mut standard_gaf_writer = BufWriter::new(&standard_gaf_file);
    let mut standard_gaf_parquet_writer = BufWriter::new(&standard_gaf_parquet_file);
    let mut comments_gaf_writer = BufWriter::new(&comments_gaf_file);
    let mut comments_gaf_parquet_writer = BufWriter::new(&comments_gaf_parquet_file);

    let generated_by = format!("!generated-by: {}\n", database_name);
    let iso_date = db_creation_datetime.replace(' ', "T");
    let date_generated = format!("!date-generated: {}\n", &iso_date);
    let url_header = format!("!URL: {}\n", &config.base_url);
    let funding_header = format!("!funding: {}\n", &config.funder);

    gpi_writer.write_all("!gpi-version: 2.0\n".as_bytes())?;
    gpi_writer.write_all(format!("!namespace: {}\n", database_name).as_bytes())?;
    gpi_writer.write_all(generated_by.as_bytes())?;
    gpi_writer.write_all(date_generated.as_bytes())?;
    gpi_writer.write_all(url_header.as_bytes())?;
    gpi_writer.write_all(funding_header.as_bytes())?;

    gpad_writer.write_all("!gpad-version: 2.0\n".as_bytes())?;
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

    let extended_pombase_gaf_header_parts =
        vec!["db", "db_object_id", "db_object_symbol", "db_object_description",
             "qualifier", "go_id", "go_term_name", "db:reference", "evidence_code",
             "with_or_from", "aspect", "db_object_name", "db_object_synonym",
             "db_object_type", "taxon", "date", "assigned_by", "annotation_extension",
             "gene_product_form_id"];
    let standard_gaf_header_parts =
        vec!["db", "db_object_id", "db_object_symbol", "qualifier", "go_id", "db:reference",
             "evidence_code", "with_or_from", "aspect", "db_object_name", "db_object_synonym",
             "db_object_type", "taxon", "date", "assigned_by", "annotation_extension",
             "gene_product_form_id"];
    let mut comments_gaf_header_parts = standard_gaf_header_parts.clone();
    comments_gaf_header_parts.push("comment_or_text_span");

    writeln!(comments_gaf_writer, "{}\n", comments_gaf_header_parts.join("\t"))?;

    let mut pombase_gaf_lines = vec![];
    let mut extended_pombase_gaf_lines = vec![];
    let mut standard_gaf_lines = vec![];
    let mut comments_gaf_lines = vec![];

    for gene_details in genes.values() {
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

        if let Some(ref characterisation_status) = gene_details.characterisation_status
            && (characterisation_status == "dubious" || characterisation_status == "transposon") {
                continue;
            }

        if gene_details.feature_type == "pseudogene" {
            continue;
        }

        write_to_gpi(&mut gpi_writer, config, api_maps, data_lookup, gene_details)?;

        write_gene_product_annotation(&mut gpad_writer, data_lookup, go_eco_mappping, config,
                                      gene_details)?;

        for aspect_name in &GO_ASPECT_NAMES {
            make_gene_association_lines(config,
                                        data_lookup,
                                        GpadGafWriteMode::PomBaseGaf,
                                        ExportComments::NoExport,
                                        gene_details, transcripts,
                                        aspect_name,
                                        &mut pombase_gaf_lines);
            make_gene_association_lines(config,
                                        data_lookup,
                                        GpadGafWriteMode::ExtendedPomBaseGaf,
                                        ExportComments::NoExport,
                                        gene_details, transcripts,
                                        aspect_name,
                                        &mut extended_pombase_gaf_lines);
            make_gene_association_lines(config,
                                        data_lookup,
                                        GpadGafWriteMode::StandardGaf,
                                        ExportComments::NoExport,
                                        gene_details, transcripts,
                                        aspect_name,
                                        &mut standard_gaf_lines);
            make_gene_association_lines(config,
                                        data_lookup,
                                        GpadGafWriteMode::StandardGaf,
                                        ExportComments::Export,
                                        gene_details, transcripts,
                                        aspect_name,
                                        &mut comments_gaf_lines);
        }
    }

    write_tsv_lines(&mut pombase_gaf_writer, &pombase_gaf_lines)?;
    write_tsv_lines(&mut extended_pombase_gaf_writer, &extended_pombase_gaf_lines)?;
    write_tsv_lines(&mut standard_gaf_writer, &standard_gaf_lines)?;
    write_tsv_lines(&mut comments_gaf_writer, &comments_gaf_lines)?;

    write_parquet(&mut pombase_gaf_parquet_writer, &standard_gaf_header_parts,
                  &pombase_gaf_lines)?;
    write_parquet(&mut extended_pombase_gaf_parquet_writer, &extended_pombase_gaf_header_parts,
                  &extended_pombase_gaf_lines)?;
    write_parquet(&mut standard_gaf_parquet_writer, &standard_gaf_header_parts,
                  &standard_gaf_lines)?;
    write_parquet(&mut comments_gaf_parquet_writer, &comments_gaf_header_parts,
                  &comments_gaf_lines)?;

    Ok(())
}

fn needs_nd_annotation(term_annotations: Option<&Vec<OntTermAnnotations>>) -> bool
{
    if let Some(annotations) = term_annotations {
        for annotation in annotations {
            if !annotation.is_not {
                return false;
            }
        }
    }

    true
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
                .any(|p| p.termid == flex_str!("GO:0032991") && p.rel_name == "is_a")
            {
                flex_str!("BFO:0000050")
            } else {
                flex_str!("RO:0002432")
            }
        },
        _ => panic!("unknown cv_name in GPAD export: {}", term_details.cv_name)
    }
}

fn get_gpad_relation_name_of(data_lookup: &dyn DataLookup,
                             term_details: &TermDetails,
                             annotation_detail: &OntAnnotationDetail) -> FlexStr {
    let rel_id = get_gpad_relation_of(term_details, annotation_detail);
    data_lookup.get_term(&rel_id)
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

fn get_gpad_nd_relation_name_of(data_lookup: &dyn DataLookup,
                                aspect: &FlexStr) -> FlexStr {
    let rel_id = get_gpad_nd_relation_of(aspect);
    data_lookup.get_term(&rel_id)
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

lazy_static! {
    static ref POMBASE_ORIG_TERM_NAME_RE: Regex = Regex::new(r#"PomBase_original_PRO_term_name: (.*)"#).unwrap();
}

fn get_pr_term_name(data_lookup: &dyn DataLookup, gene_product_form_id: &FlexStr)
    -> (FlexStr, FlexStr)
{
    let term_details = data_lookup.get_term(gene_product_form_id)
        .unwrap_or_else(|| panic!("can't find details for term {}", gene_product_form_id));

    let mut orig_term_name = term_details.name.clone();

    for syn in &term_details.synonyms {
        if syn.synonym_type == "exact" {
            let Some(captures) = POMBASE_ORIG_TERM_NAME_RE.captures(&syn.name)
            else {
                continue;
            };

            orig_term_name = captures.get(1).unwrap().as_str().into();
        }
    }

    (term_details.name.clone(), orig_term_name)
}

pub fn write_to_gpi(gpi_writer: &mut dyn Write, config: &Config, api_maps: &APIMaps,
                         data_lookup: &dyn DataLookup,
                         gene_details: &GeneDetails)
                         -> Result<(), io::Error>
{
    let database_name = &config.database_name;

    let db_object_id = format!("{}:{}", database_name, gene_details.uniquename);
    let db_object_symbol =
        gene_details.name.clone().unwrap_or(gene_details.uniquename.clone());

    let product =
        gene_details.product.clone().unwrap_or(flex_str!("unknown product"));
    let db_object_name =
        gene_details.product.clone().unwrap_or(db_object_symbol.clone());

    let mut db_object_synonyms =
        gene_details.synonyms.iter().filter(|synonym| {
            synonym.synonym_type == "exact"
        })
        .map(|synonym| synonym.name.to_string())
        .collect::<Vec<String>>();

    db_object_synonyms.sort();

    let db_object_synonyms_string = db_object_synonyms.join("|");

    let Some(transcript_so_termid) =
        gene_details.transcript_so_termid.as_ref() else {
           eprintln!("no transcript_so_termid for gene: {}", db_object_id);
           return Ok(());
        };

    let Some(db_object_type_map_val) = config.file_exports.gpad_gpi.transcript_gene_so_term_map
        .get(transcript_so_termid)
        else {
            eprintln!("failed for find configuration for {} in transcript_gene_so_term_map from {}",
                      transcript_so_termid, db_object_id);
            return Ok(());
        };

    let Some(db_object_type) = db_object_type_map_val else {
        // there is no available gene SO term for this transcript's SO type
        return Ok(());
    };

    let db_object_taxon = make_ncbi_taxon_id(config);

    let mut db_xrefs_vec = vec![];

    if let Some(ref uniprot_id) = gene_details.uniprot_identifier {
        db_xrefs_vec.push(flex_fmt!("UniProtKB:{}", uniprot_id))
    };

    if let Some(ref rnacentral_urs_identifier) = gene_details.rnacentral_urs_identifier {
        db_xrefs_vec.push(flex_fmt!("RNAcentral:{}_{}", rnacentral_urs_identifier,
                                    gene_details.taxonid))
    }

    let db_xrefs = db_xrefs_vec.join("|");

    let gpi_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t\t{}\tgo-annotation-summary={}\n",
                           db_object_id,
                           db_object_symbol,
                           db_object_name,
                           db_object_synonyms_string,
                           db_object_type,
                           db_object_taxon,
                           db_object_id,
                           db_xrefs,
                           product);
    gpi_writer.write_all(gpi_line.as_bytes())?;

    if gene_details.feature_type == "mRNA gene" {
        for (idx, transcript_id) in gene_details.transcripts.iter().enumerate() {
            let transcript_number = idx + 1;
            let transcript = api_maps.transcripts.get(transcript_id).unwrap();
            let transcript_name = transcript.name.as_ref().unwrap_or(transcript_id);

            let gpi_line = format!("{}.{}\t{}\t{}\t\tSO:0000234\t{}\t{}\t{}\t\t\tgo-annotation-summary={}\n",
                                   db_object_id, transcript_number,
                                   transcript_name,
                                   db_object_name,
                                   db_object_taxon,
                                   db_object_id,
                                   db_object_id,
                                   product);
            gpi_writer.write_all(gpi_line.as_bytes())?;
        }
    }

    let mut pr_ids_seen = HashSet::new();

    for aspect in GO_ASPECT_NAMES.iter() {
        let term_annotations =
        match gene_details.cv_annotations.get(aspect) {
            Some(term_annotations) => term_annotations.clone(),
            None => continue,
        };

        for term_annotation in term_annotations {
            for annotation_id in term_annotation.annotations {
                let maybe_annotation_detail = data_lookup.get_annotation_detail(annotation_id);
                let Some(annotation_detail) = maybe_annotation_detail.clone()
                else {
                    panic!("can't find annotation {}", annotation_id);
                };

                if let Some(ref reference) = annotation_detail.reference
                    && config.file_exports.exclude_references.contains(reference) {
                        continue;
                    }

                let gene_product_form_id = annotation_detail.gene_product_form_id.to_owned();

                if let Some(gene_product_form_id) = gene_product_form_id {
                    if pr_ids_seen.contains(&gene_product_form_id) {
                        continue;
                    }

                    pr_ids_seen.insert(gene_product_form_id.clone());

                    let (pr_term_symbol, pr_term_name) =
                        get_pr_term_name(data_lookup, &gene_product_form_id);

                    let gpi_line =
                        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t\tgo-annotation-summary={}\n",
                                gene_product_form_id,
                                pr_term_symbol,
                                pr_term_name,
                                db_object_synonyms_string,
                                db_object_type,
                                db_object_taxon,
                                db_object_id,
                                db_object_id,
                                product);
                    gpi_writer.write_all(gpi_line.as_bytes())?;
                }
            }
        }
    }

    Ok(())
}

pub fn write_gene_product_annotation(gpad_writer: &mut dyn io::Write,
                                     data_lookup: &dyn DataLookup,
                                     go_eco_mappping: &GoEcoMapping, config: &Config,
                                     gene_details: &GeneDetails)
                                     -> Result<(), io::Error>
{
    let database_name = &config.database_name;
    let db_object_id = flex_fmt!("{}:{}", database_name, gene_details.uniquename);
    let local: DateTime<Local> = Local::now();
    let local_iso_date = local.format("%F");
    let assigned_by = &config.database_name;

    if config.file_exports.include_nd_lines {
        for aspect in GO_ASPECT_NAMES.iter() {
            let term_annotations = gene_details.cv_annotations.get(aspect);
            if needs_nd_annotation(term_annotations) {
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
            let maybe_term = data_lookup.get_term(&term_annotation.term);
            let term = maybe_term
                .unwrap_or_else(|| panic!("failed to find term summary for {}",
                                 term_annotation.term));

            let mut sorted_annotations = term_annotation.annotations.clone();

            sorted_annotations.sort_by(|a1, a2| {
                let detail1 =
                    data_lookup.get_annotation_detail(*a1)
                    .unwrap_or_else(|| panic!("can't find annotation {}", a1));
                let detail2 =
                    data_lookup.get_annotation_detail(*a2)
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

            for annotation_id in sorted_annotations {
                let annotation_detail = data_lookup.get_annotation_detail(annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));
                let not = if term_annotation.is_not {
                    "NOT"
                } else {
                    ""
                };

                let relation = get_gpad_relation_of(term.as_ref(), annotation_detail.as_ref());

                let ontology_class_id = &term_annotation.term;
                let reference_uniquename =
                    annotation_detail.reference.clone()
                    .unwrap_or_else(|| flex_str!(""));
                let evidence_type =
                    eco_evidence_from_annotation(go_eco_mappping, annotation_detail.as_ref());

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
                        make_extension_string(config, data_lookup,
                                              &GpadGafWriteMode::Gpad,
                                              &annotation_detail.extension);
                    let db_object_id =
                        if let Some(ref gene_product_form_id) = annotation_detail.gene_product_form_id {
                          gene_product_form_id
                        } else {
                          &db_object_id
                        };
                    let annotation_props =
                        if let Some(ref curator_orcid) = annotation_detail.curator {
                            flex_fmt!("contributor-id=orcid:{}", curator_orcid)
                        } else {
                            flex_str!("")
                        };

                    let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\n",
                                       db_object_id,
                                       not, relation, ontology_class_id,
                                       reference_uniquename, evidence_type,
                                       with_or_from, date, assigned_by,
                                       annotation_extensions, annotation_props);
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

#[allow(clippy::too_many_arguments)]
pub fn make_gene_association_lines(config: &Config,
                                   data_lookup: &dyn DataLookup,
                                   write_mode: GpadGafWriteMode,
                                   export_comments: ExportComments,
                                   gene_details: &GeneDetails,
                                   transcripts: &UniquenameTranscriptMap,
                                   cv_name: &FlexStr,
                                   lines: &mut Vec<Vec<String>>) {
    let database_name = config.database_name.to_std_string();
    let db_object_id = gene_details.uniquename.to_std_string();
    let db_object_symbol =
        if let Some(ref gene_name) = gene_details.name {
            gene_name.to_std_string()
        } else {
            db_object_id.clone()
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
            product.to_std_string()
        } else {
            "".to_owned()
        };

    let db_object_type =
        if let Some(transcript_uniquename) = gene_details.transcripts.first() {
            let transcript_details =
                transcripts.get(transcript_uniquename)
                .unwrap_or_else(|| panic!("internal error, failed to find transcript: {}",
                                          transcript_uniquename));

            let transcript_type = transcript_details.transcript_type.as_str();

            if transcript_type == "mRNA" {
                "protein".to_owned()
            } else {
                transcript_type.to_owned()
            }
        } else {
            return;
        };

    let single_letter_aspect = abbreviation_of_go_aspect(cv_name).to_owned();

    let mut positive_annotation_count = 0;

    if let Some(term_annotations) = gene_details.cv_annotations.get(cv_name) {
        for term_annotation in term_annotations {
            let go_term = data_lookup.get_term(&term_annotation.term)
                .unwrap_or_else(|| panic!("failed to find term summary for {}",
                                          term_annotation.term));

            for annotation_id in &term_annotation.annotations {
                let annotation_detail = data_lookup.get_annotation_detail(*annotation_id)
                    .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));
                if let Some(ref reference) = annotation_detail.reference
                    && config.file_exports.exclude_references.contains(reference) {
                        continue;
                    }

                if export_comments == ExportComments::Export &&
                    let Some(ref reference_uniquename) = annotation_detail.reference {
                        if let Some(reference_details) = data_lookup.get_reference(reference_uniquename) {
                            if !reference_details.is_canto_curated() {
                                continue;
                            }
                        } else {
                            continue;
                        }
                    }

                let go_id = term_annotation.term.to_std_string();

                let term_details_arc = data_lookup.get_term(&term_annotation.term);
                let term_details_ref = term_details_arc.as_deref().unwrap().to_owned();

                let qualifiers = if write_mode == GpadGafWriteMode::PomBaseGaf ||
                    write_mode == GpadGafWriteMode::ExtendedPomBaseGaf
                {
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
                        get_gpad_relation_name_of(data_lookup,
                                                  &go_term, annotation_detail.as_ref())
                        .to_string();

                    if term_annotation.is_not {
                        format!("NOT|{}", relation_name)
                    } else {
                        relation_name
                    }
                };

                let assigned_by =
                    if let Some(ref assigned_by) = annotation_detail.assigned_by {
                        assigned_by.to_std_string()
                    } else {
                        database_name.clone()
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
                let with_or_from = itertools::join(&with_or_from_parts,",");

                let evidence_code =
                    annotation_detail.evidence
                    .as_ref()
                    .map(|s| s.as_str())
                    .unwrap_or_else(|| "")
                    .to_string();

                let date =
                    if let Some(ref raw_date) = annotation_detail.date {
                        raw_date.replace('-', "").to_string()
                    } else {
                        println!("WARNING while writing GPAD file: \
                                  date missing from annotation {}", annotation_id);
                        continue;
                    };

                let annotation_extensions =
                    make_extension_string(config, data_lookup,
                                          &write_mode, &annotation_detail.extension);

                let gene_product_form_id =
                    if let Some(ref gene_product_form_id) = annotation_detail.gene_product_form_id
                {
                    gene_product_form_id.to_std_string()
                } else {
                    "".to_owned()
                };

                let taxonid = format!("taxon:{}", gene_details.taxonid);

                if !term_annotation.is_not {
                    positive_annotation_count += 1;
                }

                let mut line_parts: Vec<String> =
                    vec![database_name.clone(),
                         db_object_id.clone(),
                         db_object_symbol.clone()];

                if write_mode == GpadGafWriteMode::ExtendedPomBaseGaf {
                    if let Some(ref product) = gene_details.product {
                        line_parts.push(product.to_std_string())
                    } else {
                        line_parts.push("".to_owned());
                    }
                }

                line_parts.push(qualifiers);
                line_parts.push(go_id);
                if write_mode == GpadGafWriteMode::ExtendedPomBaseGaf {
                    line_parts.push(term_details_ref.name.to_std_string());
                }

                if let Some(ref reference_uniquename) =
                    annotation_detail.reference {
                        line_parts.push(reference_uniquename.to_std_string());
                    } else {
                        line_parts.push("".to_owned());
                    };
                line_parts.push(evidence_code);
                line_parts.push(with_or_from);
                line_parts.push(single_letter_aspect.clone());
                line_parts.push(db_object_name.clone());
                line_parts.push(db_object_synonyms_string.clone());
                line_parts.push(db_object_type.clone());
                line_parts.push(taxonid);
                line_parts.push(date);
                line_parts.push(assigned_by);
                line_parts.push(annotation_extensions);
                line_parts.push(gene_product_form_id);
                if export_comments == ExportComments::Export {
                    if let Some(ref tmp_submitter_comment) = get_submitter_comment(annotation_detail.as_ref()) {
                        line_parts.push(tmp_submitter_comment.to_owned());
                    } else {
                        line_parts.push("".to_owned());
                    }
                }

                lines.push(line_parts);
            }
        }
    }

    if positive_annotation_count == 0 && config.file_exports.include_nd_lines &&
        write_mode == GpadGafWriteMode::StandardGaf && db_object_type == "protein" &&
        export_comments == ExportComments::NoExport
    {
        let local: DateTime<Local> = Local::now();
        let date = local.format("%Y%m%d").to_string();
        let relation_name =
            get_gpad_nd_relation_name_of(data_lookup, cv_name).to_std_string();
        let go_aspect_termid =
            config.file_exports.gpad_gpi.go_aspect_terms.get(cv_name).unwrap().to_std_string();
        let nd_ref = config.file_exports.nd_reference.clone();
        let taxonid = format!("taxon:{}", gene_details.taxonid);
        let assigned_by = database_name.clone();

        let line_parts = vec![database_name.clone(),
                              db_object_id,
                              db_object_symbol,
                              relation_name,
                              go_aspect_termid,
                              nd_ref,
                              "ND".to_owned(),
                              "".to_owned(),
                              single_letter_aspect,
                              db_object_name,
                              db_object_synonyms_string,
                              db_object_type,
                              taxonid,
                              date,
                              assigned_by,
                              "".to_owned(),
                              "".to_owned()];

        lines.push(line_parts);
    }
}

#[allow(clippy::too_many_arguments)]
pub fn write_go_annotation_format(tsv_writer: &mut dyn io::Write,
                                  config: &Config,
                                  data_lookup: &dyn DataLookup,
                                  write_mode: GpadGafWriteMode,
                                  export_comments: ExportComments,
                                  gene_details: &GeneDetails,
                                  transcripts: &UniquenameTranscriptMap,
                                  cv_name: &FlexStr)
         -> Result<(), io::Error>
{
    if !GO_ASPECT_NAMES.contains(cv_name) {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
                                  format!("unknown CV: {}", cv_name)));
    }

    let mut lines = vec![];

    make_gene_association_lines(config, data_lookup, write_mode, export_comments,
                                gene_details, transcripts, cv_name, &mut lines);

    write_tsv_lines(tsv_writer, &lines)?;

    Ok(())
}

fn make_ncbi_taxon_id(config: &Config) -> FlexStr {
    let taxon_str = format!("NCBITaxon:{}",
                            config.load_organism_taxonid.expect("internal error, no load_organism_taxonid"));

    taxon_str.into()
}
