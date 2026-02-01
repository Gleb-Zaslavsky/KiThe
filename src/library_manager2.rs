use crate::library_manager::{LibraryManager, LibraryType};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

impl LibraryManager {
    /// Recreates the AllKeysSubstance file from the SubstanceBase file.
    ///
    /// Reads the SubstanceBase file (structure: {"database_name":{"substance_name":{...}}})
    /// and creates AllKeysSubstance file (structure: [["database_name", "substance_name"], ...])
    ///
    /// # Returns
    /// * `Ok(())` - If recreation was successful
    /// * `Err(Box<dyn std::error::Error>)` - If file operations failed
    pub fn recreate_all_keys_substance(&self) -> Result<(), Box<dyn std::error::Error>> {
        let substance_base_path = self.substance_base_path();
        let all_keys_path = self.all_keys_substance_path();

        if !Path::new(substance_base_path).exists() {
            return Err(
                format!("SubstanceBase file does not exist: {}", substance_base_path).into(),
            );
        }

        let content = fs::read_to_string(substance_base_path)?;
        let substance_base: serde_json::Value = serde_json::from_str(&content)?;

        let mut all_keys = Vec::new();

        if let Some(databases) = substance_base.as_object() {
            for (database_name, substances) in databases {
                if let Some(substance_map) = substances.as_object() {
                    for substance_name in substance_map.keys() {
                        all_keys.push(vec![database_name.clone(), substance_name.clone()]);
                    }
                }
            }
        }

        let all_keys_json = serde_json::to_string_pretty(&all_keys)?;
        fs::write(all_keys_path, all_keys_json)?;

        Ok(())
    }

    /// Creates the elements base file from the SubstanceBase file.
    ///
    /// Reads the SubstanceBase file and extracts composition data to create
    /// an elements file with structure: {"element_name": [["substance_name", "database_name"], ...]}
    ///
    /// # Returns
    /// * `Ok(())` - If creation was successful
    /// * `Err(Box<dyn std::error::Error>)` - If file operations failed
    pub fn create_elements_base(&self) -> Result<(), Box<dyn std::error::Error>> {
        let substance_base_path = self.substance_base_path();
        let elements_path = self.elements_path();

        if !Path::new(substance_base_path).exists() {
            return Err(
                format!("SubstanceBase file does not exist: {}", substance_base_path).into(),
            );
        }

        let content = fs::read_to_string(substance_base_path)?;
        let substance_base: serde_json::Value = serde_json::from_str(&content)?;

        let mut elements_map: HashMap<String, Vec<Vec<String>>> = HashMap::new();

        if let Some(databases) = substance_base.as_object() {
            for (database_name, substances) in databases {
                if let Some(substance_map) = substances.as_object() {
                    for (substance_name, substance_data) in substance_map {
                        if let Some(composition) = substance_data.get("composition") {
                            if let Some(comp_str) = composition.as_str() {
                                // Parse composition string like "{Ti: 1, O: 1, Cl: 2}"
                                let comp_str = comp_str.trim_matches(|c| c == '{' || c == '}');
                                for element_pair in comp_str.split(',') {
                                    let element_pair = element_pair.trim();
                                    if let Some(colon_pos) = element_pair.find(':') {
                                        let element_name = element_pair[..colon_pos].trim();
                                        elements_map
                                            .entry(element_name.to_string())
                                            .or_insert_with(Vec::new)
                                            .push(vec![
                                                substance_name.clone(),
                                                database_name.clone(),
                                            ]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        let elements_json = serde_json::to_string_pretty(&elements_map)?;
        fs::write(elements_path, elements_json)?;

        Ok(())
    }
}
