use std::cmp::Ordering;

use chrono::NaiveDate;

use crate::data_types::*;

use flexstr::SharedStr as FlexStr;

 // Parse two date strings and compare them.  If both can't be parsed, return Equal.
pub fn cmp_str_dates(date_str1: &str, date_str2: &str) -> Ordering {
    let datetime1_res = NaiveDate::parse_from_str(date_str1, "%Y-%m-%d %H:%M:%S");
    let datetime2_res = NaiveDate::parse_from_str(date_str2, "%Y-%m-%d %H:%M:%S");

    match datetime1_res {
        Ok(datetime1) => {
            match datetime2_res {
                Ok(datetime2) => datetime1.cmp(&datetime2),
                Err(_) => Ordering::Greater
            }
        },
        Err(_) => match datetime2_res {
            Ok(_) => Ordering::Less,
            Err(_) => Ordering::Equal
        }
    }
}

pub fn remove_first<T, P>(vec: &mut Vec<T>, predicate: P) -> Option<T>
    where P: FnMut(&T) -> bool {
    if let Some(pos) = vec.iter().position(predicate) {
        return Some(vec.remove(pos));
    }

    None
}

pub fn remove_first_with_index<T, P>(vec: &mut Vec<T>, predicate: P) -> Option<(T, usize)>
    where P: FnMut(&T) -> bool {
    if let Some(pos) = vec.iter().position(predicate) {
        return Some((vec.remove(pos), pos));
    }

    None
}

#[test]
fn test_remove_first_string() {
    let mut arr =  vec!["foo", "bar", "ZZZ"];

    assert!(remove_first(&mut arr, |x| x.starts_with("_DUMMY_")).is_none());
    assert_eq!(arr.len(), 3);
    assert_eq!(remove_first(&mut arr, |x| x.to_lowercase() == "zzz").unwrap(), "ZZZ");
    assert_eq!(arr.len(), 2);
    assert_eq!(remove_first(&mut arr, |x| x.starts_with("foo")).unwrap(), "foo");
    assert_eq!(arr.len(), 1);
    assert_eq!(remove_first(&mut arr, |x| x.starts_with("b")).unwrap(), "bar");
    assert_eq!(arr.len(), 0);
}

pub fn make_gene_short(data_lookup: &dyn DataLookup,
                       gene_uniquename: &FlexStr) -> GeneShort {
    if let Some(gene_details) = data_lookup.get_gene(gene_uniquename) {
        GeneShort {
            uniquename: gene_details.uniquename.clone(),
            name: gene_details.name.clone(),
            product: gene_details.product.clone(),
            transcript_count: gene_details.transcripts.len(),
            flags: gene_details.flags.clone(),
        }
    } else {
        panic!("can't find GeneDetails for gene uniquename {}", gene_uniquename)
    }
}
