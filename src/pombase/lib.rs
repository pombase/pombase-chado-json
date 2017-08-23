#![feature(discriminant_value)]

extern crate regex;
extern crate bit_set;
extern crate chrono;
extern crate serde_json;
extern crate reqwest;

#[macro_use] extern crate serde_derive;

pub mod db;
pub mod web;
pub mod types;
pub mod interpro;
pub mod api;
