use std::{collections::{BTreeSet, HashMap, HashSet}, fs::{self, File}, io::Cursor, vec};

use anyhow::Result;

use pombase_gocam::{overlaps::GoCamNodeOverlap, parse_gocam_model, GoCamMergeAlgorithm,
                    GoCamGeneIdentifier, GoCamModel, RemoveType};
use tokio::io::AsyncReadExt as _;

use crate::data_types::{GoCamId, GoCamSummary};

pub fn read_gocam_models_from_dir(model_dir: &str)
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

pub async fn read_all_gocam_models(web_root_dir: &str,
                                   all_gocam_data: &HashMap<GoCamId, GoCamSummary>)
    -> anyhow::Result<Vec<GoCamModel>>
{
    let mut models = vec![];

    let mut flags = HashSet::new();

    flags.insert("with_chemicals".to_owned());
    flags.insert("with_inputs".to_owned());

    for gocam_id in all_gocam_data.keys() {
        let model = read_gocam_model(web_root_dir, gocam_id, &flags).await?;
        models.push(model);
    }

    Ok(models)
}

pub async fn read_gocam_model(web_root_dir: &str, gocam_id: &str, flags: &HashSet<String>)
    -> anyhow::Result<GoCamModel>
{
    let file_name = format!("{}/web-json/go-cam/gomodel:{}.json", web_root_dir, gocam_id);
    let mut source = tokio::fs::File::open(file_name).await?;
    let mut contents = vec![];
    source.read_to_end(&mut contents).await?;
    let mut cursor = Cursor::new(contents);
    let model_res = parse_gocam_model(&mut cursor);

    let mut remove_types = HashSet::new();

    if flags.contains("no_inputs") {
        remove_types.insert(RemoveType::Targets);
    }
    if flags.contains("no_chemicals") {
        remove_types.insert(RemoveType::Chemicals);
    }

    if remove_types.is_empty() {
        model_res
    } else {
        model_res.map(|model| model.remove_nodes(remove_types))
    }
}

pub async fn read_merged_gocam_model(web_root_dir: &str, all_gocam_data: &HashMap<GoCamId, GoCamSummary>,
                                     flags: &HashSet<String>,
                                     gene_list: &BTreeSet<GoCamGeneIdentifier>)
                                     -> anyhow::Result<GoCamModel>
{
    let models: Vec<_> = read_all_gocam_models(web_root_dir, all_gocam_data).await?
        .into_iter()
        .filter(|model| {
            if flags.contains("trim_models") {
                model.model_activity_enabled_by(gene_list)
            } else {
                true
            }
        })
        .collect();

    let merge_res = GoCamModel::merge_models("merged", "merged models", &models,
                                             GoCamMergeAlgorithm::Activity);

    let mut remove_types = HashSet::new();

    if flags.contains("no_chemicals") {
        remove_types.insert(RemoveType::Chemicals);
    }
    if flags.contains("no_inputs") {
        remove_types.insert(RemoveType::Targets);
    }

    if remove_types.is_empty() {
        merge_res
    } else {
        merge_res.map(|model| model.remove_nodes(remove_types))
    }
}

pub async fn read_connected_gocam_models(web_root_dir: &str,
                                         all_gocam_data: &HashMap<GoCamId, GoCamSummary>,
                                         overlaps: &Vec<GoCamNodeOverlap>,
                                         flags: &HashSet<String>)
    -> anyhow::Result<GoCamModel>
{
    let mut overlapping_gocam_ids = HashSet::new();
    for overlap in overlaps {
        for (model_id, _, _) in overlap.models.iter() {
            overlapping_gocam_ids.insert(model_id.replace("gomodel:", ""));
        }
    }

    let mut models = vec![];

    for gocam_id in all_gocam_data.keys() {
        if !overlapping_gocam_ids.contains(gocam_id.as_str()) {
            continue;
        }

        let flags = HashSet::new();
        let model = read_gocam_model(web_root_dir, gocam_id, &flags).await?;
        models.push(model);
    }

    let merge_algorithm = if flags.contains("merge_by_chemical") {
        GoCamMergeAlgorithm::Chemical
    } else {
        GoCamMergeAlgorithm::Activity
    };

    let merge_res = GoCamModel::merge_models("merged", "merged models", &models,
                                             merge_algorithm);
    let mut remove_types = HashSet::new();

    if flags.contains("no_chemicals") {
        remove_types.insert(RemoveType::Chemicals);
    }
    if flags.contains("no_inputs") {
        remove_types.insert(RemoveType::Targets);
    }

    if remove_types.is_empty() {
        merge_res
    } else {
        merge_res.map(|model| model.remove_nodes(remove_types)
                         .retain_largest_subgraph())
    }
}
