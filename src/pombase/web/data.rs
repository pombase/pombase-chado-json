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
pub type ReferenceUniquename = String;
pub type GeneName = String;
pub type TypeName = String;
pub type GeneProduct = String;

pub type TermName = String;
pub type TermId = String;
pub type TermDef = String;

pub type Evidence = String;
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
pub type UniquenameGenotypeMap = HashMap<GenotypeUniquename, GenotypeDetails>;
pub type TermIdDetailsMap = HashMap<TermId, TermDetails>;

pub type IdGeneMap = HashMap<GeneUniquename, GeneDetails>;
pub type IdGenotypeMap = HashMap<GenotypeUniquename, GenotypeDetails>;
pub type IdGeneShortMap = HashMap<GeneUniquename, GeneShort>;
pub type IdRcTermShortMap = HashMap<TermId, Rc<TermShort>>;
pub type IdRcTermDetailsMap = HashMap<TermId, Rc<TermDetails>>;
pub type IdReferenceMap = HashMap<TermId, ReferenceDetails>;

pub type ReferenceShortMap = HashMap<ReferenceUniquename, ReferenceShort>;
pub type GeneShortMap = HashMap<GeneUniquename, GeneShort>;
pub type GenotypeShortMap = HashMap<GeneUniquename, GenotypeShort>;
pub type AlleleShortMap = HashMap<AlleleUniquename, AlleleShort>;
pub type TermShortMap = HashMap<TermId, TermShort>;

include!(concat!(env!("OUT_DIR"), "/data_serde.rs"));

impl WebData {
    fn get_genes(&self) -> &IdGeneMap {
        &self.genes
    }
    fn get_genotypes(&self) -> &IdGenotypeMap {
        &self.genotypes
    }
    fn get_references(&self) -> &IdReferenceMap {
        &self.references
    }
    fn get_terms(&self) -> &IdRcTermDetailsMap {
        &self.terms
    }

    fn write_reference_details(&self, output_dir: &str) {
        for (reference_uniquename, reference_details) in &self.references {
            let s = serde_json::to_string(&reference_details).unwrap();
            let file_name = String::new() + &output_dir + "/reference/" + &reference_uniquename + ".json";
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write reference JSON");
        }
    }

    fn write_gene_details(&self, output_dir: &str) {
        for (gene_uniquename, gene_details) in &self.genes {
            let s = serde_json::to_string(&gene_details).unwrap();
            let file_name = String::new() + &output_dir + "/gene/" + &gene_uniquename + ".json";
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write gene JSON");
        }
    }

    fn write_genotype_details(&self, output_dir: &str) {
        for (genotype_uniquename, genotype_details) in &self.genotypes {
            let s = serde_json::to_string(&genotype_details).unwrap();
            let file_name = String::new() + &output_dir + "/genotype/" + &genotype_uniquename + ".json";
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write genotype JSON");
        }
    }

    fn write_gene_summaries(&self, output_dir: &str, organism_genus_species: &str) {
        let filtered_gene_summaries =
            self.search_api_maps.gene_summaries.iter()
            .filter(|gene_summary| {
                let feature_org_genus_species = String::new() +
                    &gene_summary.organism.genus + "_" + &gene_summary.organism.species;
                feature_org_genus_species == organism_genus_species
            })
            .map(|gene_summary| gene_summary.clone())
            .collect::<Vec<GeneSummary>>();
        let s = serde_json::to_string(&filtered_gene_summaries).unwrap();
        let file_name = String::new() + &output_dir + "/gene_summaries.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write gene_summaries.json");
    }

    fn write_terms(&self, output_dir: &str) {
        for (termid, term_details) in &self.terms {
            let s = serde_json::to_string(&term_details).unwrap();
            let file_name = String::new() + &output_dir + "/term/" + &termid + ".json";
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write term JSON");
        }
    }

    fn write_metadata(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.metadata).unwrap();
        let file_name = String::new() + &output_dir + "/metadata.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write metadata.json");
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
        println!("wrote {} references", self.get_references().len());
        self.write_gene_details(output_dir);
        println!("wrote {} genes", self.get_genes().len());
        self.write_genotype_details(output_dir);
        println!("wrote {} genotypes", self.get_genotypes().len());
        self.write_gene_summaries(output_dir, organism_genus_species);
        println!("wrote gene summaries");
        self.write_terms(output_dir);
        println!("wrote {} terms", self.get_terms().len());
        self.write_metadata(output_dir);
        println!("wrote metadata");
        self.write_search_api_maps(output_dir);
        println!("wrote search data");
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
