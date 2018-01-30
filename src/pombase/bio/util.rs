use std::collections::{HashSet, HashMap};
use web::data::{GeneDetails, ChromosomeLocation, DeletionViability, Strand,
                ChromosomeShort, TranscriptDetails, FeatureShort,
                ProteinDetails, FeatureType};

pub fn format_fasta(id: &str, maybe_desc: Option<String>,
                    seq: &str, width: usize) -> String {
    let mut ret = ">".to_owned() + id;

    if let Some(desc) = maybe_desc {
        ret.push(' ');
        ret.push_str(&desc);
    }

    ret.push('\n');

    if seq.len() == 0 {
        ret.push('\n');
    } else {
        let mut count = 0;
        for c in seq.chars() {
            ret.push(c);
            count += 1;
            if count % width == 0 {
                ret.push('\n');
            }
        }

        if count % width != 0 {
            ret.push('\n');
        }
    }

    ret
}

fn to_gff(source: &str, feat_type: &str,
          location: &ChromosomeLocation, maybe_parent: &Option<String>) -> String {
    let phase_char =
        if let Some(ref phase) = location.phase {
            phase.to_gff_str()
        } else {
            "."
        };
    let start_pos_str = format!("{}", location.start_pos);
    let end_pos_str = format!("{}", location.end_pos);
    let mut ret_val = String::new();
    let bits: Vec<&str> =
        vec![&location.chromosome.name, source, feat_type,
             &start_pos_str, &end_pos_str,
             ".", location.strand.to_gff_str(), phase_char];
    for s in bits {
        ret_val.push_str(s);
        ret_val.push_str("\t");
    }

    if let &Some(ref parent_id) = maybe_parent {
        ret_val.push_str("ID=");
        ret_val.push_str(&parent_id);
    };

    ret_val
}


pub fn format_gene_gff(source: &str, gene: &GeneDetails) -> Option<String> {
    if let Some (ref gene_loc) = gene.location {
        Some(to_gff(source, "gene", &gene_loc, &None))
    } else {
        None
    }
}

#[test]
fn test_format_fasta() {
    assert!(format_fasta("id1", None, "", 8) == ">id1\n\n");
    assert!(format_fasta("id2", Some("desc1".to_owned()), "atgc", 4) ==
            ">id2 desc1\natgc\n");
    assert!(format_fasta("id3", None, "atgc", 4) == ">id3\natgc\n");
    let input1 = "acgtacgattattaccggttacgcatccgtgtaaca";
    let expected1 = ">id4 desc4\nacgtacga\nttattacc\nggttacgc\natccgtgt\naaca\n";
    assert!(format_fasta("id4", Some("desc4".to_owned()), input1, 8) == expected1);
    let input2 = "acgtacgattattaccggttacgcatccgtgt";
    let expected2 = ">id5\nacgtacga\nttattacc\nggttacgc\natccgtgt\n";
    assert!(format_fasta("id5", None, input2, 8) == expected2);
}

#[allow(dead_code)]
fn make_test_gene() -> GeneDetails {
    GeneDetails {
        uniquename: "SPCC18B5.06".to_owned(),
        name: Some("dom34".to_owned()),
        taxonid: 4896,
        product: Some("Dom34-Hbs1 translation release factor complex subunit, peloto ortholog Dom34".to_owned()),
        deletion_viability: DeletionViability::Viable,
        uniprot_identifier: Some("Q9USL5".to_owned()),
        interpro_matches: vec![],
        tm_domain_coords: vec![],
        orfeome_identifier: None,
        name_descriptions: vec![],
        synonyms: vec![],
        dbxrefs: HashSet::new(),
        feature_type: "mRNA gene".to_owned(),
        characterisation_status: Some("published".to_owned()),
        location: Some(ChromosomeLocation {
            chromosome: ChromosomeShort {
                name: "chromosome_3".to_owned(),
                length: 2452883,
                ena_identifier: "CU329672.1".to_owned(),
            },
            start_pos: 729054,
            end_pos: 730829,
            strand: Strand::Forward,
            phase: None,
        }),
        gene_neighbourhood: vec![],
        transcripts: vec![
            TranscriptDetails {
                uniquename:"SPCC18B5.06.1".to_owned(),
                parts: vec![
                    FeatureShort {
                        feature_type: FeatureType::FivePrimeUtr,
                        uniquename: "SPCC18B5.06.1:five_prime_UTR:1".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729054,
                            end_pos: 729132,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "AAAATTTCGAGATATTATGTCAGTCAAAATAGCCTAAAAAATTCCTGTTCACTTAAATTCTTCGTCAACCACATTCAAT".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:1".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729133,
                            end_pos: 729212,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "ATGAAGTTGATTCAAAAAAACATCGAAAAAAATGGCTCCGGATGGATAACCATGTGCCCTGAAGAGCCAGAAGATATGTG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: "SPCC18B5.06.1:intron:1".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729213,
                            end_pos: 729265,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTATGCCTGGTTGCTGTATATCTGCAATCAATGCACAAACTAACTTAGTTTAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:2".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729266,
                            end_pos: 729319,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GCATTTGTATAATATTCTTCAAGTTGGAGATCAGCTGAAGGCTTCTACAGTTCG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: "SPCC18B5.06.1:intron:2".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729320,
                            end_pos: 729379,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTATGAAATAAAGTGTATCTTCAATGTTCGATAATAACACCTGTTTTTTCTGTTGTTTAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:3".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729380,
                            end_pos: 729678,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "TCGTGTAGTGAAAGTTGGCGCTACAGGAAGTACGTCAGGTTCAAGAGTTGTGATGAAACTACGTATTTTAGTTGAGAATATGGACTTTGATACAAAGGCTGCTCAATTGCACATCAAAGGACGGACAACTGAATACCATCCTGAAGTTAAGATGGGATCCTACCATACCTTGGACTTAGAACTACATCGCAATTTTACTCTATATAAAAATGAATGGGATGCATTTGCATTGGACCGTGTAGATGCTGCTTGTAATCCTTCAAGAAATGCTGAAATAGGTGCTGTGGTTTTAGATGAAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: "SPCC18B5.06.1:intron:3".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729679,
                            end_pos: 729734,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTAAATTTCTTGGCTTACATCGCTTTATTAGTTCGTGCATGTTTACTGACTTGCAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:4".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729735,
                            end_pos: 729841,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTCTTGCCAATATTTGTCTCATTACAGATTATATGACCATCCTGCGTCAAAGGATTGATCAAGTGATTCCAAGGAAACGGAGAGGGGACAGCAGCGCTTACCAAAAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: "SPCC18B5.06.1:intron:4".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729842,
                            end_pos: 729889,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTAAATTTTTAGACTTTGATTTTTCGTCCGTACTGATCAATTTTTTAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:5".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729890,
                            end_pos: 730522,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GGCCTTGATAAATTTTATGACTCTGTTTTTCAATCCATTAACTCAGAATTTGATTTTGATAAATTGAAAGTCGTTATTCTTGCTTCACCAGGATTTGTGGCTCGAGGCTTGTATGACTACATATTCAGCATGGCCGTGAAGTTAGACTTGAAACAAATTGTTAAATCAAAGAATAAATTTGTCATCCTTCATTCTAGCACTGGTCATATTCATTCCCTTAATGAAATTTTGAAGGACCCTGCTGTTGAATCAAAACTAGCCGACACAAAATACGTACAAGAAATTCGCGTTCTGAATAAATTTTACGATGTCATGAATGAAGATGATAGAAAGGCATGGTATGGTCCAAATCATGTTTTGAAGGCTTTTGAACTTGGCGCGATCGGAGAACTTCTGATTAGCGATTCTCTGTTCAGGAGTTCTGACATTGCTACTAGAAAAAAATGGGTTTCATTAGTAGAAGGTGTTAAGGAGATTAACTGTCCTGTTTATATTTTCAGTAGTTTGCATGAGTCAGGGAAGCAGCTGGATCTGTTGTCAGGTATTGCCGCCATTCTCACTTACCCAGTCGATGAAGAGGATATATCAGAAGATGAAGAGGATGAGGAATCCCAAAATTTTGAACATAGTTAA".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::ThreePrimeUtr,
                        uniquename: "SPCC18B5.06.1:three_prime_UTR:1".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 730523,
                            end_pos: 730829,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "AGTTCATCAGTATCCGAATTGTCATGAATCTAATTATTGCTAAGCCAATATTTCATACTTTAAGCTCGGTTAGAACAATTTGTTTCATTCTCTTAAAAAATTTATTTATGGGCTCGTTTTGGTAGTCATTATTTATGCTTTTACTTGGATGTTTTAGGGTATTTACTATGATAAACATGCAAAAAATTAGGTGTTAGAATGGTCAAAAATTGATACCCTAAAATGATTTATACTTATCGATTATAACTCTTAACTTGTAAAATTAAGCTGTTAATTATAGCCGGTCCAATACAGTATTCAATTACAG".to_owned(),
                    }
                ],
                transcript_type: "mRNA".to_owned(),
                protein: Some(ProteinDetails {
                    uniquename: "SPCC18B5.06.1:pep".to_owned(),
                    sequence: "MKLIQKNIEKNGSGWITMCPEEPEDMWHLYNILQVGDQLKASTVRRVVKVGATGSTSGSRVVMKLRILVENMDFDTKAAQLHIKGRTTEYHPEVKMGSYHTLDLELHRNFTLYKNEWDAFALDRVDAACNPSRNAEIGAVVLDEGLANICLITDYMTILRQRIDQVIPRKRRGDSSAYQKGLDKFYDSVFQSINSEFDFDKLKVVILASPGFVARGLYDYIFSMAVKLDLKQIVKSKNKFVILHSSTGHIHSLNEILKDPAVESKLADTKYVQEIRVLNKFYDVMNEDDRKAWYGPNHVLKAFELGAIGELLISDSLFRSSDIATRKKWVSLVEGVKEINCPVYIFSSLHESGKQLDLLSGIAAILTYPVDEEDISEDEEDEESQNFEHS*".to_owned(),
                    molecular_weight: 44.3361,
                    average_residue_weight: 0.11339,
                    charge_at_ph7: -8.45,
                    isoelectric_point: 5.64,
                    codon_adaptation_index: 0.59
                }),
                cds_location: Some(ChromosomeLocation {
                    chromosome: ChromosomeShort {
                        name: "chromosome_3".to_owned(),
                        length: 2452883,
                        ena_identifier: "CU329672.1".to_owned(),
                    },
                    start_pos: 729133,
                    end_pos: 730522,
                    strand: Strand::Forward,
                    phase: None,
                }),
            },
        ],
        cv_annotations: HashMap::new(),
        physical_interactions: vec![],
        genetic_interactions: vec![],
        ortholog_annotations: vec![],
        paralog_annotations: vec![],
        target_of_annotations: vec![],
        references_by_uniquename: HashMap::new(),
        genes_by_uniquename: HashMap::new(),
        genotypes_by_uniquename: HashMap::new(),
        alleles_by_uniquename: HashMap::new(),
        terms_by_termid: HashMap::new(),
    }
}

#[test]
fn test_format_gff() {
    let gene = make_test_gene();
    assert!(format_gene_gff("PomBase", &gene).unwrap() == "chromosome_3\tPomBase\tgene\t729054\t730829\t.\t+\t.\t");
}
