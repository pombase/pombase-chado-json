use std::collections::HashMap;

use flexstr::SharedStr as FlexStr;
use rusqlite::Connection;

use crate::{constants::API_MAPS_TABLE_NAMES, data_types::{TermIdDetailsMap, UniquenameGeneMap, UniquenameReferenceMap, IdGenotypeMap, APIGenotypeAnnotation}, types::TermId};

pub fn join(v: &[FlexStr], connector: &str) -> FlexStr {
    let result = itertools::join(v.iter().map(FlexStr::as_ref), connector);
    result.into()
}

pub fn make_maps_database_tables(conn: &mut Connection) -> rusqlite::Result<()> {
    let tx = conn.transaction()?;

    for table_name in API_MAPS_TABLE_NAMES.iter() {
        tx.execute(
            &format!("CREATE TABLE {} (
                        id    TEXT PRIMARY KEY,
                        data  TEXT NOT NULL
                     )",
                     table_name),
            (),
        )?;
    }

    tx.commit()?;

    Ok(())
}


pub fn store_maps_into_database(conn: &mut Connection, terms: &TermIdDetailsMap,
                                genes: &UniquenameGeneMap, references: &UniquenameReferenceMap,
                                genotypes: &IdGenotypeMap,
                                termid_genotype_annotation: &HashMap<TermId, Vec<APIGenotypeAnnotation>>) -> anyhow::Result<()> {
    let tx = conn.transaction()?;

    for (termid, term_details) in terms {
        let json = serde_json::value::to_value(term_details)?;

        tx.execute("INSERT INTO terms (id, data) VALUES (?1, ?2)",
                   (termid.as_ref(), &json))?;
    }
    for (gene_uniquename, gene_details) in genes {
        let json = serde_json::value::to_value(gene_details)?;

        tx.execute("INSERT INTO genes (id, data) VALUES (?1, ?2)",
                   (gene_uniquename.as_ref(), &json))?;
    }
    for (ref_uniquename, reference_details) in references {
        let json = serde_json::value::to_value(reference_details)?;

        tx.execute("INSERT INTO refs (id, data) VALUES (?1, ?2)",
                   (ref_uniquename.as_ref(), &json))?;
    }
    for (genotype_display_uniquename, genotype_details) in genotypes {
        let json = serde_json::value::to_value(genotype_details)?;

        tx.execute("INSERT INTO genotypes (id, data) VALUES (?1, ?2)",
                   (genotype_display_uniquename.as_ref(), &json))?;
    }
    for (termid, genotype_annotations) in termid_genotype_annotation {
        let json = serde_json::value::to_value(genotype_annotations)?;

        tx.execute("INSERT INTO termid_genotype_annotations (id, data) VALUES (?1, ?2)",
                   (termid.as_ref(), &json))?;
    }

    tx.commit()?;

    Ok(())
}
