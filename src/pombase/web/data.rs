extern crate serde_json;
extern crate postgres;

use std::fs::File;
use std::io::{Write, BufWriter};
use std::collections::HashMap;

use self::postgres::Connection;

type CvName = String;

pub type MiscExtRange = String;

pub type GeneUniquename = String;
pub type TermUniquename = String;
pub type GeneName = String;
pub type TypeName = String;
pub type GeneProduct = String;

pub type TermName = String;
pub type TermId = String;
pub type TermDef = String;

pub type Evidence = String;
pub type With = String;
pub type Residue = String;
pub type Condition = String;
pub type Qualifier = String;

pub type TypeInteractionAnnotationMap =
    HashMap<TypeName, Vec<InteractionAnnotation>>;

pub type UniquenameGeneMap =
    HashMap<GeneUniquename, GeneDetails>;

pub type TranscriptUniquename = String;

pub type UniquenameTranscriptMap =
    HashMap<TranscriptUniquename, TranscriptDetails>;

pub type GenotypeUniquename = String;
pub type AlleleUniquename = String;

pub type UniquenameAlleleMap = HashMap<AlleleUniquename, AlleleShort>;
pub type UniquenameGenotypeMap = HashMap<GenotypeUniquename, GenotypeShort>;

pub type IdGeneMap = HashMap<GeneUniquename, GeneDetails>;
pub type IdGeneShortMap = HashMap<GeneUniquename, GeneShort>;
pub type IdTermDetailsMap = HashMap<TermId, Rc<TermDetails>>;
pub type IdReferenceMap = HashMap<TermId, ReferenceDetails>;

include!(concat!(env!("OUT_DIR"), "/data_serde.rs"));

impl WebData {
    fn get_genes(&self) -> &IdGeneMap {
        &self.genes
    }
    fn get_terms(&self) -> &IdTermDetailsMap {
        &self.terms
    }

    fn write_reference_details(&self, output_dir: &str) {
        for (reference_uniquename, reference_details) in &self.references {
            let s = serde_json::to_string(&reference_details).unwrap();
            let file_name = String::new() + &output_dir + "/reference/" + &reference_uniquename + ".json";
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write!");
        }
    }

    fn write_gene_details(&self, output_dir: &str) {
        for (gene_uniquename, gene_details) in &self.genes {
            let s = serde_json::to_string(&gene_details).unwrap();
            let file_name = String::new() + &output_dir + "/gene/" + &gene_uniquename + ".json";
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write!");
        }
    }

    fn write_gene_summary(&self, output_dir: &str, organism_genus_species: &str) {
        let gene_summaries =
            self.search_api_maps.gene_summaries.values()
            .filter(|gene_short| {
                let gene_uniquename = &gene_short.uniquename;
                let gene_details = self.get_genes().get(gene_uniquename).unwrap();
                let feature_org_genus_species = String::new() +
                    &gene_details.organism.genus + "_" + &gene_details.organism.species;
                feature_org_genus_species == organism_genus_species
            })
            .map(|gene_short| gene_short.clone())
            .collect::<Vec<GeneShort>>();
        let s = serde_json::to_string(&gene_summaries).unwrap();
        let file_name = String::new() + &output_dir + "/gene_summaries.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }

    fn write_terms(&self, output_dir: &str) {
        for (termid, term_details) in &self.terms {
            let s = serde_json::to_string(&term_details).unwrap();
            let file_name = String::new() + &output_dir + "/term/" + &termid + ".json";
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write!");
        }
    }

    fn write_metadata(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.metadata).unwrap();
        let file_name = String::new() + &output_dir + "/metadata.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }

    fn write_search_api_maps(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.search_api_maps).unwrap();
        let file_name = String::new() + &output_dir + "/search_api_maps.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write!");
    }

    pub fn write(&self, output_dir: &str, organism_genus_species: &str) {
        self.write_reference_details(output_dir);
        self.write_gene_details(output_dir);
        self.write_gene_summary(output_dir, organism_genus_species);
        self.write_terms(output_dir);
        self.write_metadata(output_dir);
        self.write_search_api_maps(output_dir);

        println!("wrote {} genes", self.get_genes().len());
        println!("wrote {} terms", self.get_terms().len());
    }

    pub fn store_jsonb(&self, conn: &Connection) {
        let trans = conn.transaction().unwrap();

        for (uniquename, gene_details) in &self.genes {
            let serde_value = serde_json::value::to_value(&gene_details);
            trans.execute("INSERT INTO web_json.gene (uniquename, data) values ($1, $2)",
                          &[&uniquename, &serde_value]).unwrap();
        }
        for (uniquename, ref_details) in &self.references {
            let serde_value = serde_json::value::to_value(&ref_details);
            trans.execute("INSERT INTO web_json.reference (uniquename, data) values ($1, $2)",
                          &[&uniquename, &serde_value]).unwrap();
        }
        for (termid, term_details) in &self.terms {
            let serde_value = serde_json::value::to_value(&term_details);
            trans.execute("INSERT INTO web_json.term (termid, data) values ($1, $2)",
                         &[&termid, &serde_value]).unwrap();
        }

        trans.execute("CREATE INDEX gene_jsonb_idx ON web_json.gene USING gin (data jsonb_path_ops)", &[]).unwrap();
        trans.execute("CREATE INDEX gene_jsonb_name_idx ON web_json.gene USING gin ((data->>'name') gin_trgm_ops);", &[]).unwrap();
        trans.execute("CREATE INDEX term_jsonb_idx ON web_json.term USING gin (data jsonb_path_ops)", &[]).unwrap();
        trans.execute("CREATE INDEX term_jsonb_name_idx ON web_json.term USING gin ((data->>'name') gin_trgm_ops);", &[]).unwrap();
        trans.execute("CREATE INDEX reference_jsonb_idx ON web_json.reference USING gin (data jsonb_path_ops)", &[]).unwrap();
        trans.execute("CREATE INDEX reference_jsonb_title_idx ON web_json.reference USING gin ((data->>'title') gin_trgm_ops);", &[]).unwrap();

        trans.commit().unwrap();
    }
}
