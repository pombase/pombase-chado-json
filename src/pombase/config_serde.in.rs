// configuration for extension display names and for the "Target of" section
#[derive(Deserialize, Clone, Debug)]
pub struct ExtensionConfig {
    pub rel_name: String, // name of extension relation
    pub display_name: String, // text to display
    pub if_descendent_of: Option<String>, // None if applies to any extension
    pub reciprocal_display: Option<String>, // None if reciprocal shouldn't be displayed
}

#[derive(Deserialize, Clone, Debug)]
pub struct Config {
    pub extensions: Vec<ExtensionConfig>,
}
