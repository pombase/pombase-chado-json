use std::fs::{self, File};

use pombase_gocam::{gocam_parse, GoCamRawModel};

use anyhow::Result;

pub fn read_gocam_models(model_dir: &str)
    -> Result<Vec<GoCamRawModel>>
{
    let mut ret = vec![];

    let entries = fs::read_dir(model_dir)?;

    for entry in entries {
        let file_name = entry?.path();
        if file_name.to_string_lossy().ends_with(".json") {
            let mut file = File::open(file_name)?;
            let model = gocam_parse(&mut file)?;
            ret.push(model);
        }
    }

    Ok(ret)
}
