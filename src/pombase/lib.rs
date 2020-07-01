#![feature(nll)]

extern crate jemallocator;

#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

extern crate regex;
extern crate bit_set;
extern crate chrono;
extern crate serde_json;
extern crate reqwest;
extern crate flate2;
#[macro_use] extern crate lazy_static;
#[macro_use] extern crate serde_derive;
extern crate pombase_rc_string;
extern crate uuid;
extern crate postgres;

pub mod db;
pub mod web;
pub mod types;
pub mod interpro;
pub mod api;
pub mod bio;
pub mod rnacentral;
pub mod annotation_util;
pub mod data_types;
pub mod api_data;
