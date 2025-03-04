use std::fs::{self, File};

use anyhow::Result;

use pombase_gocam::{make_gocam_model, GoCamModel};

pub fn read_gocam_models(model_dir: &str)
    -> Result<Vec<GoCamModel>>
{
    let mut ret = vec![];

    let entries = fs::read_dir(model_dir)?;

    for entry in entries {
        let file_name = entry?.path();
        if file_name.to_string_lossy().ends_with(".json") {
            let mut file = File::open(file_name)?;
            let model = make_gocam_model(&mut file)?;
            ret.push(model);
        }
    }

    Ok(ret)
}
