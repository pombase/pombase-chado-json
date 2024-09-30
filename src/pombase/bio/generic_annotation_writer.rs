use std::io::{self, BufWriter, Write};

use crate::{data_types::PeptideRange,
            uniprot::{UniProtDataMap, UniProtDataEntry}};

use chrono::prelude::{Local, DateTime};

pub struct UniProtTermidMap {
  pub glycosylation_site_termid: String,
  pub disulphide_bond_termid: String,
}

fn get_aa_at_pos(sequence: &str, range: &PeptideRange)
    -> char
{
    if range.len() != 1 {
        panic!("get_aa_at_pos() needs a position, not: {}..{}",
               range.start, range.end);
    }
    let pos = range.start - 1;
    sequence.as_bytes()[pos] as char
}

fn aa_for_residue(sequence: &str, termid: &str, range: &PeptideRange)
    -> Option<char>
{
    match termid {
        "MOD:00042"
            => Some('D'),
        "MOD:00047"
            => Some('T'),
        "MOD:00048" | "MOD:00156"
            => Some('Y'),
        "MOD:00050"
            => Some('A'),
        "MOD:00046" | "MOD:00060" | "MOD:00159" | "MOD:00171"
            => Some('S'),
        "MOD:00080"
            => Some('Q'),
        "MOD:00111" | "MOD:00114" | "MOD:00115" | "MOD:00234" |
        "MOD:00257" | "MOD:00441" | "MOD:00689"
            => Some('C'),
        "MOD:00064" | "MOD:00083" | "MOD:00084" | "MOD:00085" |
        "MOD:00123" | "MOD:00125" | "MOD:00126" | "MOD:00127" |
        "MOD:00128"
            => Some('K'),
        "MOD:00153" | "MOD:00226" | "MOD:00890"
            => Some('H'),
        "MOD:00304"
            => Some('L'),
        "MOD:00310"
            => Some('R'),
        "MOD:00068" | "MOD:00351" | "MOD:01625" | "MOD:01982"
            => Some('G'),
        "MOD:00818" => {
            let char_at_pos = get_aa_at_pos(sequence, range);
            match char_at_pos {
                'G'|'D'|'S'|'N'|'C' => Some(char_at_pos),
                _ => {
                    eprintln!("unexpected char {} as {}", char_at_pos, range);
                    Some('?')
                }
            }
        },
        "MOD:01154" => {
            let char_at_pos = get_aa_at_pos(sequence, range);
            match char_at_pos {
                'S'|'T'|'K' => Some(char_at_pos),
                _ => {
                    eprintln!("unexpected char {} as {}", char_at_pos, range);
                    Some('?')
                },
            }
        },
        _ => None,
    }
}

fn get_residue_extension(uniprot_data: &UniProtDataEntry,
                         termid: &str, range: &PeptideRange)
    -> String
{
    let residue_abbrev = aa_for_residue(uniprot_data.sequence.as_str(),
                                        termid, range);
    if let Some(residue_abbrev) = residue_abbrev {
        format!("residue({}{})", residue_abbrev, range)
    } else {
        "".to_owned()
    }
}

pub fn write_from_uniprot_map(uniprot_data_map: &UniProtDataMap,
                              uniprot_pmid: &str,
                              termid_map: &UniProtTermidMap,
                              assigned_by: &str,
                              out: &mut dyn Write)
   -> Result<(), io::Error>
{
  let mut writer = BufWriter::new(out);

  let local: DateTime<Local> = Local::now();
  let date = local.format("%F").to_string();

  for uniprot_data in uniprot_data_map.values() {
    for site in &uniprot_data.glycosylation_sites {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or(uniprot_pmid);
      let termid = &termid_map.glycosylation_site_termid;
      let residue_extension = format!("residue(N{})", site.range);
      write_generic_annotation(&mut writer,
                               &uniprot_data.gene_uniquename,
                               "",
                               termid,
                               evidence,
                               reference,
                               &date,
                               "",
                               &residue_extension,
                               assigned_by)?;
    }
    for site in &uniprot_data.disulfide_bonds {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or(uniprot_pmid);
      let termid = &termid_map.disulphide_bond_termid;
      let residue_extension = get_residue_extension(uniprot_data, termid,
                                                    &site.range);
      write_generic_annotation(&mut writer,
                               &uniprot_data.gene_uniquename,
                               "",
                               termid,
                               evidence,
                               reference,
                               &date,
                               "",
                               &residue_extension,
                               assigned_by)?;
    }
    for site in &uniprot_data.lipidation_sites {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or(uniprot_pmid);
      let termid = &site.termid;
      let residue_extension = get_residue_extension(uniprot_data, termid,
                                                    &site.range);
      write_generic_annotation(&mut writer,
                               &uniprot_data.gene_uniquename,
                               "",
                               termid,
                               evidence,
                               reference,
                               &date,
                               "",
                               &residue_extension,
                               assigned_by)?;
    }
    for site in &uniprot_data.modified_residues {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or(uniprot_pmid);
      let termid = &site.termid;
      let residue_extension = get_residue_extension(uniprot_data, termid,
                                                    &site.range);
      write_generic_annotation(&mut writer,
                               &uniprot_data.gene_uniquename,
                               "",
                               &site.termid,
                               evidence,
                               reference,
                               &date,
                               "",
                               &residue_extension,
                               assigned_by)?;
    }
  }

  Ok(())
}

fn write_generic_annotation(out: &mut BufWriter<&mut dyn Write>,
                            uniquename: &str, name: &str, termid: &str,
                            evidence: &str,
                            reference_uniquename: &str,
                            date: &str,
                            qualifiers: &str,
                            extension: &str,
                            assigned_by: &str)
         -> Result<(), io::Error>
{
  writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
           uniquename, name, termid, evidence,
           reference_uniquename, date, qualifiers, extension, assigned_by)?;
  Ok(())
}
