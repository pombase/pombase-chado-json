extern crate pombase;

use std::error::Error;
use std::env;

use std::fs::File;
use std::io;
use std::process;

use flexstr::shared_fmt as flex_fmt;
use getopts::Options;

use getopts::ParsingStyle;

use pombase::bio::generic_annotation_writer::write_from_uniprot_map;
use pombase::bio::generic_annotation_writer::UniProtTermidMap;
use pombase::bio::util::read_fasta;
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

    if input_file_type != "uniprot-data-tsv" {
        eprintln!("unknown input file type {}\n", input_file_type);
        eprint_usage(&program, opts);
        process::exit(1);
    }

    if args.is_empty() {
        eprintln!(r#"missing [file_name] argument after "{}"\n"#, input_file_type);
        eprint_usage(&program, opts);
        process::exit(1);
    }

    if input_file_type == "uniprot-data-tsv" {
        let mut sub_opts = Options::new();
        let sub_opts = sub_opts.parsing_style(ParsingStyle::StopAtFirstFree);

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

        let sub_matches = match sub_opts.parse(&args) {
            Ok(m) => m,
            Err(e) => {
                eprint_usage(&program, sub_opts);
                println!("\nerror: {}", e);
                process::exit(1);
            },
        };

        let termid_map = UniProtTermidMap {
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

        let file_name_pos = 5 + if filter_references.len() > 0 { 1 } else { 0 };
        let file_name = &args[file_name_pos];

        let uniprot_data_map =
            filter_uniprot_map(parse_uniprot(file_name, &filter_references),
                               &peptides);

        let mut stdout = io::stdout().lock();
        write_from_uniprot_map(&uniprot_data_map, &uniprot_pmid,
                               &termid_map, &assigned_by,
                               &mut stdout)?;
    } else {
        eprintln!("unknown input file type: {}", input_file_type);
        eprint_usage(&program, opts);
        process::exit(1);
    }

    Ok(())
}
