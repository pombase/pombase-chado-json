use std::error::Error;

use deadpool_postgres::Client;

use crate::data_types::AlleleShortMap;
use crate::db::Processed;
use crate::types::OrganismTaxonId;

pub struct Loader {
    conn: Client,
    taxonid: OrganismTaxonId,
    processed: Processed,
}

impl Loader {
    pub fn new(conn: Client, taxonid: OrganismTaxonId, processed: Processed) -> Loader {
        Loader {
            conn,
            taxonid,
            processed,
        }
    }

    pub async fn load_alleles(&mut self, alleles: &AlleleShortMap) -> Result<(), Box<dyn Error>>
    {
        let organism = self.processed.organism_by_taxonid(self.taxonid)
            .expect(&format!("load failed, organism not in database with taxon ID: {}", self.taxonid));

        let trans = self.conn.transaction().await?;

        let allele_insert_sql = "INSERT INTO feature(uniquename, name, type_id, organism_id)
VALUES ($1, $2,
        (SELECT cvterm_id FROM cvterm t join cv ON cv.cv_id = t.cv_id WHERE t.name = 'allele' AND cv.name = 'sequence' LIMIT 1),
        (SELECT organism_id FROM organism WHERE genus = $3 AND species = $4 LIMIT 1))";

        let allele_smt = trans.prepare(allele_insert_sql).await?;

        let rel_sql = "INSERT INTO feature_relationship(subject_id, object_id, type_id)
VALUES ((SELECT feature_id FROM feature f join cvterm t ON f.type_id = t.cvterm_id WHERE t.name = 'allele' AND f.uniquename = $1),
        (SELECT feature_id FROM feature f join cvterm t ON f.type_id = t.cvterm_id WHERE (t.name = 'gene' OR t.name = 'pseudogene') AND f.uniquename = $2),
        (SELECT cvterm_id FROM cvterm t join cv ON cv.cv_id = t.cv_id WHERE t.name = 'instance_of' AND cv.name = 'pombase_relations' LIMIT 1))";

        let rel_smt = trans.prepare(rel_sql).await?;

        let add_feat_prop_sql = "INSERT INTO featureprop(feature_id, type_id, value)
VALUES ((SELECT feature_id FROM feature f join cvterm t ON f.type_id = t.cvterm_id WHERE t.name = 'allele' AND f.uniquename = $1),
        (SELECT cvterm_id FROM cvterm t join cv ON cv.cv_id = t.cv_id WHERE t.name = $2 AND cv.name = 'PomBase feature property types' LIMIT 1),
        $3)";

        let add_feat_prop_smt = trans.prepare(add_feat_prop_sql).await?;

        let add_synonym_sql =
            "INSERT INTO synonym (name, type_id, synonym_sgml)
                                   VALUES ($1, (SELECT cvterm_id FROM cvterm t JOIN cv ON cv.cv_id = t.cv_id
                                                 WHERE t.name = 'exact' AND cv.name = 'synonym_type'), '')
                              ON CONFLICT DO NOTHING";

        let add_synonym_smt = trans.prepare(add_synonym_sql).await?;

        let add_pub_sql =
            "INSERT INTO pub (uniquename, type_id)
                         VALUES ($1, (SELECT cvterm_id FROM cvterm WHERE name = 'paper'))
                        ON CONFLICT DO NOTHING";

        let add_pub_smt = trans.prepare(add_pub_sql).await?;

        let add_feat_synonym_sql =
            "INSERT INTO feature_synonym (feature_id, synonym_id, pub_id)
                 SELECT (SELECT feature_id FROM feature f WHERE f.uniquename = $1),
                        (SELECT synonym_id FROM synonym
                              JOIN cvterm t ON t.cvterm_id = synonym.type_id
                              JOIN cv ON cv.cv_id = t.cv_id
                             WHERE synonym.name = $2 AND t.name = 'exact' AND cv.name = 'synonym_type'),
                        (SELECT pub_id FROM pub WHERE uniquename = $3)";

        let add_feat_synonym_smt = trans.prepare(add_feat_synonym_sql).await?;

        for allele in alleles.values() {
            if self.processed.allele_by_uniquename(&allele.uniquename).is_some() {
                panic!("allele loading failed: can't find {} in the database", allele.uniquename);
            }

            if self.processed.gene_by_uniquename(&allele.gene_uniquename).is_none() {
                panic!("allele loading failed: can't find {} in the database", allele.gene_uniquename);
            }

            let allele_uniquename = allele.uniquename.as_str();
            let allele_name = allele.name.as_ref().map(|s| s.as_str());
            let genus = organism.genus.as_str();
            let species = organism.species.as_str();
            let gene_uniquename = allele.gene_uniquename.as_str();

            trans.query(&allele_smt,
                        &[&allele_uniquename, &allele_name, &genus, &species]).await?;

            trans.query(&rel_smt, &[&allele_uniquename, &gene_uniquename]).await?;

            let allele_type = allele.allele_type.as_str();
            trans.query(&add_feat_prop_smt, &[&allele_uniquename, &"allele_type", &allele_type]).await?;

            let allele_description = allele.description.as_ref().map(|s| s.as_str());
            trans.query(&add_feat_prop_smt, &[&allele_uniquename, &"description", &allele_description]).await?;

            for synonym in &allele.synonyms {
              let synonym_name = synonym.name.as_str();

              trans.query(&add_synonym_smt, &[&synonym_name]).await?;

              let pub_uniquename =
                  if let Some(ref uniquename) = synonym.reference {
                    uniquename.as_str()
                  } else {
                    "null"
                  };

              trans.query(&add_pub_smt, &[&pub_uniquename]).await?;

              trans.query(&add_feat_synonym_smt, &[&allele_uniquename, &synonym_name, &pub_uniquename]).await?;
            }
        }

        trans.commit().await?;

        Ok(())
    }
}
