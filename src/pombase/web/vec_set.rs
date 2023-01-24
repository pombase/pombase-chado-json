use std::hash::Hash;
use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::collections::HashMap;
use bit_set::BitSet;

// A data structure for lookup to find supersets of sets.
// iff a superset of the argument was inserted earlier.

#[derive(Debug, Default)]
pub struct VecSet<T: Eq + Hash + Clone> {
    bit_sets: Vec<BitSet>,
    data: HashMap<T, usize>,
    curr_index: usize,
}

impl<T: Eq + Hash + Clone> VecSet<T> {
    pub fn new() -> VecSet<T> {
        VecSet {
            bit_sets: vec![],
            data: HashMap::new(),
            curr_index: 0,
        }
    }

    fn make_bit_set(&mut self, tvec: &[T]) -> BitSet {
        let mut new_bitset = BitSet::new();
        for t in tvec {
            let mut curr = self.curr_index;
            let index = match self.data.entry(t.clone()) {
                Vacant(entry) => {
                    curr += 1;
                    entry.insert(curr);
                    curr
                }
                Occupied(entry) => *entry.get()
            };

            self.curr_index = curr;

            new_bitset.insert(index);
        }
        new_bitset
    }

    // Insert the argument set
    pub fn insert(&mut self, tvec: &[T]) {
        let new_bitset = self.make_bit_set(tvec);
        self.bit_sets.push(new_bitset);
    }

    // Return true iff a superset of the argument was inserted
    // with insert() earlier
    pub fn contains_superset(&mut self, tvec: &[T]) -> bool {
        let new_bitset = self.make_bit_set(tvec);
        self.bit_sets.iter().any(|s| new_bitset.is_subset(s))
    }

    // Return true iff a proper superset of the argument was inserted
    // with insert() earlier
    pub fn contains_proper_superset(&mut self, tvec: &[T]) -> bool {
        let new_bitset = self.make_bit_set(tvec);
        self.bit_sets.iter().any(|s| new_bitset.is_subset(s) && &new_bitset != s)
    }
}

#[test]
fn test_insert() {
    let mut vec_set: VecSet<String> = VecSet::new();

    let check = |vec_set: &VecSet<String>, d, b, n| {
        assert_eq!(vec_set.data.keys().len(), d);
        assert_eq!(vec_set.bit_sets.len(), b);
        assert_eq!(vec_set.curr_index, n);
    };

    let str_data =
        |s_vec: Vec<&str>| s_vec.iter()
            .map(|s| String::from(*s)).collect::<Vec<String>>();

    let s1 = str_data(vec!["one","two","three","four"]);
    vec_set.insert(&s1);
    check(&vec_set, 4,1,4);

    let s2 = str_data(vec!["one","two","four","five"]);
    vec_set.insert(&s2);
    check(&vec_set, 5,2,5);

    let s4 = str_data(vec!["one","two","three","five"]);
    vec_set.insert(&s4);
    check(&vec_set, 5,3,5);

    let s6 = str_data(vec!["one","four"]);
    vec_set.insert(&s6);
    check(&vec_set, 5,4,5);

    let s7 = str_data(vec!["one","six"]);
    vec_set.insert(&s7);
    check(&vec_set, 6,5,6);

    let s8 = str_data(vec!["one"]);
    vec_set.insert(&s8);
    check(&vec_set, 6,6,6);
}

#[test]
fn test_contains_superset() {
    let mut vec_set: VecSet<String> = VecSet::new();

    let str_data =
        |s_vec: Vec<&str>| s_vec.iter()
            .map(|s| String::from(*s)).collect::<Vec<String>>();

    let s0 = str_data(vec!["one","seven","three","four"]);
    vec_set.insert(&s0);
    let s1 = str_data(vec!["one","two","three","four"]);
    let s2 = str_data(vec!["one","two","three"]);
    let s2_2 = str_data(vec!["one","two","five"]);
    let s3 = str_data(vec!["one"]);
    vec_set.insert(&s1);
    assert!(vec_set.contains_proper_superset(&s2));
    vec_set.insert(&s2);
    assert!(vec_set.contains_proper_superset(&s2));
    vec_set.insert(&s2_2);
    assert!(!vec_set.contains_proper_superset(&s2_2));
    assert!(vec_set.contains_superset(&s2_2));
    vec_set.insert(&s3);
    assert!(vec_set.contains_proper_superset(&s3));

    let s4 = str_data(vec!["one","two"]);
    assert!(vec_set.contains_proper_superset(&s4));

    let s5 = str_data(vec!["five"]);
    assert!(vec_set.contains_proper_superset(&s5));

    let s6 = str_data(vec!["six"]);
    assert!(!vec_set.contains_proper_superset(&s6));
}
