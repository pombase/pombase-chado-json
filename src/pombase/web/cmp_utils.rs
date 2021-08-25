use regex::Regex;
use std::cmp::Ordering;

use crate::types::*;

lazy_static! {
    static ref MODIFICATION_RE: Regex = Regex::new(r"^(?P<aa>[A-Z])(?P<pos>\d+)$").unwrap();
}

pub fn cmp_residues(residue1: &Option<Residue>, residue2: &Option<Residue>) -> Ordering {
    if let Some(ref res1) = *residue1 {
        if let Some(ref res2) = *residue2 {
            if let (Some(res1_captures), Some(res2_captures)) =
                (MODIFICATION_RE.captures(res1), MODIFICATION_RE.captures(res2))
            {
                let res1_aa = res1_captures.name("aa").unwrap().as_str();
                let res2_aa = res2_captures.name("aa").unwrap().as_str();
                let aa_order = res1_aa.cmp(res2_aa);
                if aa_order == Ordering::Equal {
                    let res1_pos =
                        res1_captures.name("pos").unwrap().as_str().parse::<i32>().unwrap();
                    let res2_pos =
                        res2_captures.name("pos").unwrap().as_str().parse::<i32>().unwrap();
                    res1_pos.cmp(&res2_pos)
                } else {
                    aa_order
                }
            } else {
                res1.cmp(res2)
            }
        } else {
            Ordering::Less
        }
    } else {
        if residue2.is_some() {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}
