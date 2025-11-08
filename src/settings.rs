//! # Settings Module
//!
//! ## Purpose
//! Provides a user-friendly interface for managing JSON library file versions in KiThe.
//! This module acts as a high-level wrapper around the LibraryManager, offering
//! GUI-friendly display names and simplified operations for library version management.
//!
//! ## Key Features
//! - **Display Name Mapping**: Maps user-friendly names to internal library keys
//! - **Batch Operations**: Supports updating multiple libraries simultaneously
//! - **Configuration Sync**: Automatically syncs with the underlying LibraryManager
//! - **Error Handling**: Provides string-based error messages suitable for GUI display
//!
//! ## Usage Pattern
//! ```rust
//! use KiThe::settings::Settings;
//!
//! let mut settings = Settings::new();
//! settings.set_library_version("Substance Base", "substance_base_v3.json")?;
//! settings.reset_to_defaults()?;
//! ```
//!
//! ## Library Mappings
//! | Display Name | Internal Key | Default File |
//! |--------------|--------------|-------------|
//! | "Substance Base" | "substance_base" | substance_base_v2.json |
//! | "All Keys Substance" | "all_keys_substance" | all_keys_substance.json |
//! | "Reactbase" | "reactbase" | Reactbase.json |
//! | "Dict Reaction" | "dict_reaction" | dict_reaction.json |

use crate::library_manager::{LibraryConfig, with_library_manager, with_library_manager_mut};
use std::collections::HashMap;

/// User-friendly settings interface for managing library versions.
///
/// This struct provides a simplified API for GUI applications to manage
/// JSON library file versions without dealing with internal library manager complexity.
/// It maintains a local cache of library versions mapped to display names.
#[derive(Debug, Clone, Default)]
pub struct Settings {
    /// Cache of library versions using display names as keys
    /// Maps "Substance Base" -> "substance_base_v2.json", etc.
    pub library_versions: HashMap<String, String>,
}

impl Settings {
    /// Creates a new Settings instance by loading current configuration from LibraryManager.
    ///
    /// This constructor fetches the current library configuration and populates
    /// the internal cache with display name mappings.
    ///
    /// # Returns
    /// A new Settings instance with current library versions loaded
    ///
    /// # Example
    /// ```rust
    /// let settings = Settings::new();
    /// assert_eq!(settings.library_versions.len(), 4);
    /// ```
    pub fn new() -> Self {
        let config = with_library_manager(|manager| manager.get_config().clone());

        let mut library_versions = HashMap::new();
        library_versions.insert("Substance Base".to_string(), config.substance_base.clone());
        library_versions.insert(
            "All Keys Substance".to_string(),
            config.all_keys_substance.clone(),
        );
        library_versions.insert("Reactbase".to_string(), config.reactbase.clone());
        library_versions.insert("Dict Reaction".to_string(), config.dict_reaction.clone());

        Self { library_versions }
    }

    /// Retrieves the current file path for a given library display name.
    ///
    /// # Arguments
    /// * `library_name` - Display name of the library (e.g., "Substance Base")
    ///
    /// # Returns
    /// * `Some(&String)` - Current file path if library exists
    /// * `None` - If library name is not recognized
    ///
    /// # Example
    /// ```rust
    /// let settings = Settings::new();
    /// let path = settings.get_library_version("Substance Base");
    /// assert_eq!(path, Some(&"substance_base_v2.json".to_string()));
    /// ```
    pub fn get_library_version(&self, library_name: &str) -> Option<&String> {
        self.library_versions.get(library_name)
    }

    /// Updates a library to use a new file version.
    ///
    /// This method validates that the new file exists, updates the underlying
    /// LibraryManager configuration, and syncs the local cache.
    ///
    /// # Arguments
    /// * `library_name` - Display name of the library (e.g., "Substance Base")
    /// * `version` - New file path to use for this library
    ///
    /// # Returns
    /// * `Ok(())` - If update was successful
    /// * `Err(String)` - If library name is unknown or file doesn't exist
    ///
    /// # Example
    /// ```rust
    /// let mut settings = Settings::new();
    /// settings.set_library_version("Substance Base", "substance_base_v3.json")?;
    /// ```
    pub fn set_library_version(&mut self, library_name: &str, version: &str) -> Result<(), String> {
        let result = with_library_manager_mut(|manager| match library_name {
            "Substance Base" => manager.set_substance_base(version),
            "All Keys Substance" => manager.set_all_keys_substance(version),
            "Reactbase" => manager.set_reactbase(version),
            "Dict Reaction" => manager.set_dict_reaction(version),
            _ => return Err(format!("Unknown library: {}", library_name).into()),
        });

        match result {
            Ok(_) => {
                self.library_versions
                    .insert(library_name.to_string(), version.to_string());
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        }
    }

    /// Returns a list of all available library display names.
    ///
    /// # Returns
    /// Vector of references to library display names
    ///
    /// # Example
    /// ```rust
    /// let settings = Settings::new();
    /// let libraries = settings.get_available_libraries();
    /// assert!(libraries.contains(&&"Substance Base".to_string()));
    /// ```
    pub fn get_available_libraries(&self) -> Vec<&String> {
        self.library_versions.keys().collect()
    }

    /// Resets all library versions to their default values.
    ///
    /// This method restores the default configuration (substance_base_v2.json, etc.)
    /// and updates both the LibraryManager and local cache.
    ///
    /// # Returns
    /// * `Ok(())` - If reset was successful
    /// * `Err(String)` - If reset operation failed
    ///
    /// # Example
    /// ```rust
    /// let mut settings = Settings::new();
    /// settings.reset_to_defaults()?;
    /// ```
    pub fn reset_to_defaults(&mut self) -> Result<(), String> {
        with_library_manager_mut(|manager| match manager.reset_to_defaults() {
            Ok(_) => {
                let config = manager.get_config();
                self.library_versions.clear();
                self.library_versions
                    .insert("Substance Base".to_string(), config.substance_base.clone());
                self.library_versions.insert(
                    "All Keys Substance".to_string(),
                    config.all_keys_substance.clone(),
                );
                self.library_versions
                    .insert("Reactbase".to_string(), config.reactbase.clone());
                self.library_versions
                    .insert("Dict Reaction".to_string(), config.dict_reaction.clone());
                Ok(())
            }
            Err(e) => Err(e.to_string()),
        })
    }

    /// Updates multiple libraries simultaneously in a single atomic operation.
    ///
    /// This method is optimized for programmatic use where multiple library
    /// versions need to be changed at once. All files are validated before
    /// any changes are made.
    ///
    /// # Arguments
    /// * `updates` - HashMap mapping display names to new file paths
    ///
    /// # Returns
    /// * `Ok(())` - If all updates were successful
    /// * `Err(String)` - If any file doesn't exist or library name is unknown
    ///
    /// # Example
    /// ```rust, ignore
    /// let mut settings = Settings::new();
    /// let mut updates = HashMap::new();
    /// updates.insert("Substance Base", "substance_base_v3.json");
    /// updates.insert("Reactbase", "Reactbase_v2.json");
    /// settings.update_multiple_libraries(updates)?;
    /// ```
    pub fn update_multiple_libraries(
        &mut self,
        updates: HashMap<&str, &str>,
    ) -> Result<(), String> {
        // Convert display names to internal keys
        let mut internal_updates = HashMap::new();
        for (display_name, path) in &updates {
            let internal_key = match *display_name {
                "Substance Base" => "substance_base",
                "All Keys Substance" => "all_keys_substance",
                "Reactbase" => "reactbase",
                "Dict Reaction" => "dict_reaction",
                _ => return Err(format!("Unknown library: {}", display_name)),
            };
            internal_updates.insert(internal_key, *path);
        }

        with_library_manager_mut(|manager| {
            match manager.update_libraries(internal_updates) {
                Ok(_) => {
                    // Update local cache
                    for (display_name, path) in updates {
                        self.library_versions
                            .insert(display_name.to_string(), path.to_string());
                    }
                    Ok(())
                }
                Err(e) => Err(e.to_string()),
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::library_manager::LibraryManager;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_manager() -> LibraryManager {
        let (substance_file, keys_file, react_file, dict_file) = create_test_files();
        let mut config_file = NamedTempFile::new().unwrap();

        let config = LibraryConfig {
            substance_base: substance_file.path().to_str().unwrap().to_string(),
            all_keys_substance: keys_file.path().to_str().unwrap().to_string(),
            reactbase: react_file.path().to_str().unwrap().to_string(),
            dict_reaction: dict_file.path().to_str().unwrap().to_string(),
        };

        let config_json = serde_json::to_string_pretty(&config).unwrap();
        config_file.write_all(config_json.as_bytes()).unwrap();

        LibraryManager::with_config_file(config_file.path().to_str().unwrap())
    }

    fn create_test_files() -> (NamedTempFile, NamedTempFile, NamedTempFile, NamedTempFile) {
        let mut substance_file = NamedTempFile::new().unwrap();
        let mut keys_file = NamedTempFile::new().unwrap();
        let mut react_file = NamedTempFile::new().unwrap();
        let mut dict_file = NamedTempFile::new().unwrap();

        substance_file.write_all(b"{}").unwrap();
        keys_file.write_all(b"[]").unwrap();
        react_file.write_all(b"{}").unwrap();
        dict_file.write_all(b"{}").unwrap();

        (substance_file, keys_file, react_file, dict_file)
    }

    #[test]
    fn test_settings_new() {
        crate::library_manager::set_test_manager(create_test_manager());
        let settings = Settings::new();
        assert_eq!(settings.library_versions.len(), 4);
        assert!(settings.library_versions.contains_key("Substance Base"));
        assert!(settings.library_versions.contains_key("All Keys Substance"));
        assert!(settings.library_versions.contains_key("Reactbase"));
        assert!(settings.library_versions.contains_key("Dict Reaction"));
        crate::library_manager::clear_test_manager();
    }

    #[test]
    fn test_get_library_version() {
        crate::library_manager::set_test_manager(create_test_manager());
        let settings = Settings::new();
        let version = settings.get_library_version("Substance Base");
        assert!(version.is_some());
        crate::library_manager::clear_test_manager();
    }

    #[test]
    fn test_set_library_version() {
        crate::library_manager::set_test_manager(create_test_manager());
        let mut temp_file = NamedTempFile::new().unwrap();
        temp_file.write_all(b"{}").unwrap();

        let mut settings = Settings::new();
        let result =
            settings.set_library_version("Substance Base", temp_file.path().to_str().unwrap());

        assert!(result.is_ok());
        assert_eq!(
            settings.get_library_version("Substance Base").unwrap(),
            temp_file.path().to_str().unwrap()
        );
        crate::library_manager::clear_test_manager();
    }

    #[test]
    fn test_update_multiple_libraries() {
        crate::library_manager::set_test_manager(create_test_manager());
        let mut temp_file1 = NamedTempFile::new().unwrap();
        let mut temp_file2 = NamedTempFile::new().unwrap();
        temp_file1.write_all(b"{}").unwrap();
        temp_file2.write_all(b"[]").unwrap();

        let mut settings = Settings::new();
        let mut updates = HashMap::new();
        updates.insert("Substance Base", temp_file1.path().to_str().unwrap());
        updates.insert("All Keys Substance", temp_file2.path().to_str().unwrap());

        let result = settings.update_multiple_libraries(updates);
        assert!(result.is_ok());

        assert_eq!(
            settings.get_library_version("Substance Base").unwrap(),
            temp_file1.path().to_str().unwrap()
        );
        assert_eq!(
            settings.get_library_version("All Keys Substance").unwrap(),
            temp_file2.path().to_str().unwrap()
        );
        crate::library_manager::clear_test_manager();
    }
}
