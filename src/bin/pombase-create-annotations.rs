extern crate pombase;

use std::collections::HashMap;
use std::error::Error;
use std::env;

use std::fs::File;
use std::io;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::process;

use chrono::DateTime;
use chrono::Local;
use flexstr::shared_fmt as flex_fmt;
use getopts::Options;

use getopts::ParsingStyle;

use pombase::bio::generic_annotation_writer::write_from_uniprot_map;
use pombase::bio::generic_annotation_writer::UniProtTermidMap;
use pombase::bio::generic_annotation_writer::write_generic_annotation;
use pombase::bio::util::read_fasta;
use pombase::interpro::DomainData;
use pombase::interpro::InterProMatch;
use pombase::uniprot::filter_uniprot_map;
use pombase::uniprot::parse_uniprot;

fn usage_message(program: &str) -> String {
    format!("Usage: {} [input_file_type] [file_name]

Reads a data file and writes an annotation format lines to STDOUT.
The output is suitable for piping to pombase-import.pl with the generic-annotation
option.
",
            program)
}

fn print_usage(program: &str, opts: &Options) {
    let message = usage_message(program);
    print!("{}", opts.usage(&message));
}

fn eprint_usage(program: &str, opts: &Options) {
    let message = usage_message(program);
    eprint!("{}", opts.usage(&message));
}


fn main() -> Result<(), Box<dyn Error>> {
    let mut args: Vec<String> = env::args().collect();
    let mut opts = Options::new();
    let opts = opts.parsing_style(ParsingStyle::StopAtFirstFree);

    opts.optflag("h", "help", "print this help message");

    let program = args.remove(0);

    let matches = match opts.parse(&args) {
        Ok(m) => m,
        Err(e) => {
            eprint_usage(&program, opts);
            println!("\noption error: {}", e);
            process::exit(1);
        }
    };

    if matches.opt_present("help") {
        print_usage(&program, opts);
        process::exit(0);
    }

    if args.len() < 2 {
        println!("needs [input_file_type] and [file_name] arguments");
        eprint_usage(&program, opts);
        process::exit(1);
    }

    let input_file_type = args.remove(0);

    match input_file_type.as_str() {
        "uniprot-data-tsv" => {
            process_uniprot(&program, &args)?;
        },
        "interpro-domain-json" => {
            process_interpro_json(&program, &args)?;
        },
        _ => {
            eprintln!("unknown input file type: {}", input_file_type);
            eprint_usage(&program, opts);
            process::exit(1);
        }
    }

    Ok(())
}

fn process_uniprot(program: &str, args: &[String]) -> Result<(), Box<dyn Error>> {
    let mut sub_opts = Options::new();
    let sub_opts = sub_opts.parsing_style(ParsingStyle::StopAtFirstFree);

    sub_opts.reqopt("", "n-glycsylated-residue-termid",
                    "The term ID for N-glycsylated residue", "TERMID");

    sub_opts.reqopt("", "glycosylation-site-termid",
                    "The term ID to use in glycosylation modification annotations",
                    "TERMID");

    sub_opts.reqopt("", "disulphide-bond-termid",
                    "The term ID to use in disulphide bond annotations",
                    "TERMID");

    sub_opts.reqopt("", "peptide-fasta",
                    "A file of peptides to compare to UniProt sequences", "FASTA_FILE");

    sub_opts.reqopt("", "uniprot-reference",
                    "The PubMed ID of the most recent UniProt publication to add to the annotations",
                    "PMID");

    sub_opts.reqopt("", "assigned-by",
                    "Value for the assigned_by column", "SOURCE");

    sub_opts.optopt("", "filter-references",
                    "A comma separated list of PMIDs to ignore", "PMIDs");

    let sub_matches = match sub_opts.parse(args) {
        Ok(m) => m,
        Err(e) => {
            eprint_usage(program, sub_opts);
            println!("\nerror: {}", e);
            process::exit(1);
        },
    };

    let termid_map = UniProtTermidMap {
        n_glycsylated_residue_termid: sub_matches.opt_str("n-glycsylated-residue-termid").unwrap(),
        glycosylation_site_termid: sub_matches.opt_str("glycosylation-site-termid").unwrap(),
        disulphide_bond_termid: sub_matches.opt_str("disulphide-bond-termid").unwrap(),
    };

    let uniprot_pmid = sub_matches.opt_str("uniprot-reference").unwrap();
    let assigned_by = sub_matches.opt_str("assigned-by").unwrap();
    let peptide_filename = sub_matches.opt_str("peptide-fasta").unwrap();
    let mut peptide_file = File::open(peptide_filename)?;
    let peptides = read_fasta(&mut peptide_file)?;
    let filter_references: Vec<_> = sub_matches.opt_str("filter-references")
        .unwrap_or_default()
        .trim()
        .split(",")
        .map(|s| flex_fmt!("PMID:{}", s.trim_start_matches("PMID:")))
        .collect();

    let file_name_pos = 6 + if !filter_references.is_empty() { 1 } else { 0 };
    let file_name = &args[file_name_pos];

    let uniprot_data_map =
        filter_uniprot_map(parse_uniprot(file_name, &filter_references),
                           &peptides);

    let mut stdout = io::stdout().lock();
    write_from_uniprot_map(&uniprot_data_map, &uniprot_pmid,
                           &termid_map, &assigned_by,
                           &mut stdout)?;

    Ok(())
}

fn process_interpro_json(program: &str, args: &[String]) -> Result<(), Box<dyn Error>> {
    let mut sub_opts = Options::new();
    let sub_opts = sub_opts.parsing_style(ParsingStyle::StopAtFirstFree);

    sub_opts.reqopt("", "tm-assigned-by",
                    "Assigned by field value for transmembrane domains", "ASSIGNED_BY");


    let sub_matches = match sub_opts.parse(args) {
        Ok(m) => m,
        Err(e) => {
            eprint_usage(program, sub_opts);
            println!("\nerror: {}", e);
            process::exit(1);
        },
    };

    let tm_assigned_by = sub_matches.opt_str("tm-assigned-by").unwrap();

    let file_name_pos = 1;
    let file_name = &args[file_name_pos];

    let interpro_data = parse_interpro(file_name);

    let mut stdout = io::stdout().lock();
    write_interpro_domains(interpro_data, &tm_assigned_by,
                           &mut stdout)?;


    Ok(())
}

pub fn parse_interpro(file_name: &str) -> DomainData {
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            println!("Failed to read {}: {}", file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    let domain_data: DomainData =
        match serde_json::from_reader(reader) {
            Ok(uniprot_results) => uniprot_results,
            Err(err) => {
                print!("failed to parse {}: {}", file_name, err);
                process::exit(1);
            },
        };

    let mut filtered_domains = HashMap::new();

    for (gene_uniquename, mut results) in domain_data.domains_by_id {
        let new_interpro_matches: Vec<InterProMatch> =
            results.interpro_matches.into_iter()
            .filter(|interpro_match| interpro_match.dbname == "DeepTMHMM-Signal-Peptide")
            .collect();

        results.interpro_matches = new_interpro_matches;
        filtered_domains.insert(gene_uniquename, results);
    }

    DomainData {
        interproscan_version: domain_data.interproscan_version,
        domains_by_id: filtered_domains,
    }
}

fn write_interpro_domains(domains: DomainData,
                          tm_assigned_by: &str,
                          out: &mut dyn Write)
                          -> Result<(), io::Error>
{
    let mut writer = BufWriter::new(out);

    let local: DateTime<Local> = Local::now();
    let date = local.format("%F").to_string();

    for (gene_uniquename, entry) in domains.domains_by_id.into_iter() {
        for domain in &entry.interpro_matches {
            for loc in &domain.locations {
                let residue = if loc.start == loc.end {
                    format!("residue({})", loc.start)
                } else {
                    format!("residue({}-{})", loc.start, loc.end)
                };

                write_generic_annotation(&mut writer,
                                         &format!("{}.1", gene_uniquename),
                                         "",
                                         "SO:0000418",
                                         "ECO:0000255",
                                         "",
                                         &date,
                                         "",
                                         &residue,
                                         tm_assigned_by)?;
            }
        }
    }
    Ok(())
}
