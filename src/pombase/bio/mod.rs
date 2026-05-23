use std::fs::File;
use std::sync::Arc;
use std::io::{self, BufWriter};

use arrow::datatypes::Date32Type;
use chrono::NaiveDate;

use arrow_array::{Array, Date32Array, LargeStringArray, RecordBatch};
use arrow_schema::{DataType, Field, Schema};
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;
use parquet::basic::{Compression, ZstdLevel};

use crate::{bio::util::{ANNOTATION_COMMENT_NESTED_BRACKETS_RE, ANNOTATION_COMMENT_RE}, data_types::OntAnnotationDetail};

pub mod util;
pub mod go_format_writer;
pub mod phenotype_format_writer;
pub mod pdb_reader;
pub mod protein_view;
pub mod macromolecular_complexes;
pub mod generic_annotation_writer;
pub mod gene_expression_writer;
pub mod gocam_model_process;
pub mod complementation;
pub mod modifications;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum ExportCommentsMode {
    Export,
    NoExport,
}

fn get_submitter_comment(annotation_detail: &OntAnnotationDetail)
   -> Option<String>
{
    let submitter_comment = annotation_detail.submitter_comment.as_ref()?;

    // remove comments while handling nested brackets (in a hacky way)
    let submitter_comment =
        ANNOTATION_COMMENT_NESTED_BRACKETS_RE.replace_all(submitter_comment, "");
    let submitter_comment =
        ANNOTATION_COMMENT_RE.replace_all(&submitter_comment, "");
    let submitter_comment = submitter_comment.trim();

    if submitter_comment.is_empty() {
        None
    } else {
        Some(submitter_comment.to_owned())
    }
}

pub fn write_parquet(writer: &mut BufWriter<&File>,
                     header_parts: &[&str],
                     lines: &[Vec<String>])
  -> Result<(), io::Error>
{
    let mut schema_fields = vec![];

    let mut date_idx: i64 = -1;

    for (idx, header_part) in header_parts.iter().enumerate() {
        let field = if header_part.to_lowercase() == "date" {
            date_idx = idx as i64;
            Field::new(header_part.to_owned(), DataType::Date32, false)
        } else {
            Field::new(header_part.to_owned(), DataType::LargeUtf8, false)
        };
        schema_fields.push(field);
    }
    let schema = Arc::new(Schema::new(schema_fields));

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
                if idx as i64 == date_idx {
                    let date_iter = v.iter()
                        .map(|d| {
                            let naive_date =
                                if d.contains("-") {
                                    NaiveDate::parse_from_str(d, "%Y-%m-%d")
                                } else {
                                    NaiveDate::parse_from_str(d, "%Y%m%d")
                                }
                                .unwrap_or_else(|_| panic!("failed to parse date: {}", d));
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
