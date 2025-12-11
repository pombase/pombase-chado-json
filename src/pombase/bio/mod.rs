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


#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum ExportComments {
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

