use std::{collections::HashSet, fs::{self, File}, io::Cursor, vec};

use anyhow::Result;

use pombase_gocam::{parse_gocam_model, GoCamModel, GoCamNodeOverlap, GoCamNodeType};
use tokio::io::AsyncReadExt as _;

use crate::data_types::GoCamDetails;

pub fn read_gocam_models(model_dir: &str)
    -> Result<Vec<GoCamModel>>
{
    let mut ret = vec![];

    let entries = fs::read_dir(model_dir)?;

    for entry in entries {
        let file_name = entry?.path();
        if file_name.to_string_lossy().ends_with(".json") {
            let mut file = File::open(file_name)?;
            let model = parse_gocam_model(&mut file)?;
            ret.push(model);
        }
    }

    Ok(ret)
}

pub async fn read_gocam_model(web_root_dir: &str, gocam_id: &str) -> anyhow::Result<GoCamModel> {
    let file_name = format!("{}/web-json/go-cam/gomodel:{}.json", web_root_dir, gocam_id);
    let mut source = tokio::fs::File::open(file_name).await?;
    let mut contents = vec![];
    source.read_to_end(&mut contents).await?;
    let mut cursor = Cursor::new(contents);
    parse_gocam_model(&mut cursor).into()
}

pub async fn read_merged_gocam_model(web_root_dir: &str, all_gocam_data: &Vec<GoCamDetails>)
    -> anyhow::Result<GoCamModel>
{
    let mut models = vec![];

    for detail in all_gocam_data {
        let gocam_id = &detail.gocam_id;
        let model = read_gocam_model(web_root_dir, gocam_id).await?;
        models.push(model);
    }

    let merge_res = GoCamModel::merge_models("merged", "merged models", &models);

    merge_res
}

pub async fn read_connected_gocam_models(web_root_dir: &str,
                                         all_gocam_data: &Vec<GoCamDetails>,
                                         overlaps: &Vec<GoCamNodeOverlap>)
    -> anyhow::Result<GoCamModel>
{
    let mut overlapping_gocam_ids = HashSet::new();
    for overlap in overlaps {
        if let GoCamNodeType::Activity(_) = overlap.node_type {
            for (model_id, _) in overlap.models.iter() {
                overlapping_gocam_ids.insert(model_id.replace("gomodel:", ""));
            }
        }
    }

    let mut models = vec![];

    for detail in all_gocam_data {
        let gocam_id = &detail.gocam_id;

        if !overlapping_gocam_ids.contains(gocam_id.as_str()) {
            continue;
        }

        let model = read_gocam_model(web_root_dir, gocam_id).await?;
        models.push(model);
    }

    let merge_res = GoCamModel::merge_models("merged", "merged models", &models);

    merge_res
}
