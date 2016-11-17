extern crate serde_codegen;

use std::env;
use std::path::Path;

fn main() {
    let out_dir = env::var_os("OUT_DIR").unwrap();

    let src = Path::new("src/pombase/web/data_serde.in.rs");
    let dst = Path::new(&out_dir).join("data_serde.rs");

    serde_codegen::expand(&src, &dst).unwrap();
}
