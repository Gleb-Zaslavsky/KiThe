//! # Library Manager Module
//!
//! ## Purpose
//! Provides centralized management of JSON library file paths for the KiThe chemical
//! thermodynamics and kinetics crate. This module eliminates hardcoded file paths
//! throughout the codebase and enables dynamic library version switching.
//!
//! ## Architecture
//! - **LibraryConfig**: Serializable configuration structure
//! - **LibraryManager**: Core manager with file validation and persistence
//! - **Global Access**: Thread-safe singleton pattern with test isolation
//! - **Configuration File**: JSON-based persistent storage (library_config.json)
//!
//! ## Key Features
//! - **Centralized Paths**: All JSON library paths managed in one location
//! - **Version Switching**: Easy switching between library file versions
//! - **File Validation**: Ensures files exist before updating configuration
//! - **Thread Safety**: Safe concurrent access using Mutex and OnceLock
//! - **Test Isolation**: Separate test manager to prevent test interference
//! - **Persistence**: Automatic saving/loading of configuration
//!
//! ## Configuration Format
//! ```json
//! {
//!   "substance_base": "substance_base_v2.json",
//!   "all_keys_substance": "all_keys_substance.json",
//!   "reactbase": "Reactbase.json",
//!   "dict_reaction": "dict_reaction.json"
//! }
//! ```
//!
//! ## Usage Patterns
//!
//! ### Read-only Access
//! ```rust
//! use KiThe::library_manager::with_library_manager;
//!
//! let path = with_library_manager(|manager| {
//!     manager.substance_base_path().to_string()
//! });
//! ```
//!
//! ### Mutable Access
//! ```rust
//! use KiThe::library_manager::with_library_manager_mut;
//!
//! with_library_manager_mut(|manager| {
//!     manager.set_substance_base("substance_base_v3.json")
//! })?;
//! ```

use chrono::{DateTime, Local};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::sync::{Mutex, OnceLock};

/// Enum representing the different types of libraries in the system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LibraryType {
    SubstanceBase,
    AllKeysSubstance,
    Elements,
    Reactbase,
    DictReaction,
}

impl LibraryType {
    /// Returns the pseudonym string for the library type.
    pub fn as_str(&self) -> &'static str {
        match self {
            LibraryType::SubstanceBase => "substance_base",
            LibraryType::AllKeysSubstance => "all_keys_substance",
            LibraryType::Elements => "elements",
            LibraryType::Reactbase => "reactbase",
            LibraryType::DictReaction => "dict_reaction",
        }
    }
}

/// Configuration structure for JSON library file paths.
///
/// This struct defines the mapping between internal library keys and their
/// corresponding JSON file paths. It supports serialization for persistent
/// storage in library_config.json.
///
/// # Fields
/// * `substance_base` - Path to main substance thermodynamics database
/// * `all_keys_substance` - Path to substance key mapping file  
/// * `reactbase` - Path to chemical reaction kinetics database
/// * `dict_reaction` - Path to reaction dictionary/mapping file
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibraryConfig {
    pub substance_base: String,
    pub all_keys_substance: String,
    pub elements: String,
    pub reactbase: String,
    pub dict_reaction: String,
    pub problems_folder: String,
}

impl Default for LibraryConfig {
    /// Creates default configuration with standard KiThe library file names.
    ///
    /// # Returns
    /// LibraryConfig with default file paths:
    /// - substance_base_v2.json
    /// - all_keys_substance.json  
    /// - Reactbase.json
    /// - dict_reaction.json
    fn default() -> Self {
        Self {
            substance_base: "substance_base_v2.json".to_string(),
            all_keys_substance: "all_keys_substance.json".to_string(),
            elements: "elements.json".to_string(),
            reactbase: "Reactbase.json".to_string(),
            dict_reaction: "dict_reaction.json".to_string(),
            problems_folder: "problems".to_string(),
        }
    }
}

/// Core library manager responsible for JSON file path management.
///
/// This struct handles loading, saving, and updating library configuration.
/// It maintains both the current configuration and the path to the config file
/// for persistence operations.
///
/// # Fields
/// * `config` - Current library configuration
/// * `config_file` - Path to the configuration file for persistence
#[derive(Debug, Clone)]
pub struct LibraryManager {
    config: LibraryConfig,
    config_file: String,
}
impl Default for LibraryManager {
    fn default() -> Self {
        Self::new()
    }
}
impl LibraryManager {
    /// Creates a new LibraryManager with default configuration file.
    ///
    /// Attempts to load configuration from "library_config.json" in the current
    /// directory. If the file doesn't exist or is invalid, uses default configuration.
    ///
    /// # Returns
    /// New LibraryManager instance with loaded or default configuration
    pub fn new() -> Self {
        let config_file = "library_config.json".to_string();
        let config = Self::load_config(&config_file).unwrap_or_default();

        let manager = Self {
            config,
            config_file,
        };

        // Ensure problems folder exists
        let _ = manager.ensure_problems_folder();

        manager
    }

    /// Creates a new LibraryManager with a custom configuration file path.
    ///
    /// This constructor is primarily used for testing or when a non-standard
    /// configuration file location is required.
    ///
    /// # Arguments
    /// * `config_file` - Path to the configuration file
    ///
    /// # Returns
    /// New LibraryManager instance with configuration loaded from specified file
    pub fn with_config_file(config_file: &str) -> Self {
        let config = Self::load_config(config_file).unwrap_or_default();

        Self {
            config,
            config_file: config_file.to_string(),
        }
    }

    /// Loads configuration from a JSON file.
    ///
    /// Attempts to read and parse the configuration file. If the file doesn't exist
    /// or parsing fails, returns the default configuration instead of an error.
    ///
    /// # Arguments
    /// * `config_file` - Path to the configuration file
    ///
    /// # Returns
    /// * `Ok(LibraryConfig)` - Loaded or default configuration
    /// * `Err(Box<dyn std::error::Error>)` - Only on severe I/O errors
    fn load_config(config_file: &str) -> Result<LibraryConfig, Box<dyn std::error::Error>> {
        if Path::new(config_file).exists() {
            let content = fs::read_to_string(config_file)?;
            let config: LibraryConfig = serde_json::from_str(&content)?;
            Ok(config)
        } else {
            Ok(LibraryConfig::default())
        }
    }

    /// Saves current configuration to the config file.
    ///
    /// Serializes the current configuration to JSON and writes it to the config file.
    /// During tests, this method does nothing to prevent pollution of the real config file.
    ///
    /// # Returns
    /// * `Ok(())` - If save was successful or during tests
    /// * `Err(Box<dyn std::error::Error>)` - If file write or JSON serialization failed
    pub fn save_config(&self) -> Result<(), Box<dyn std::error::Error>> {
        #[cfg(test)]
        {
            // Don't save config during tests to avoid polluting the real config file
            return Ok(());
        }

        #[cfg(not(test))]
        {
            let content = serde_json::to_string_pretty(&self.config)?;
            fs::write(&self.config_file, content)?;
            Ok(())
        }
    }

    /// Returns the current path to the substance base JSON file.
    ///
    /// # Returns
    /// String slice containing the file path (e.g., "substance_base_v2.json")
    pub fn substance_base_path(&self) -> &str {
        &self.config.substance_base
    }

    /// Returns the current path to the substance keys mapping JSON file.
    ///
    /// # Returns  
    /// String slice containing the file path (e.g., "all_keys_substance.json")
    pub fn all_keys_substance_path(&self) -> &str {
        &self.config.all_keys_substance
    }

    /// Returns the current path to the reaction base JSON file.
    ///
    /// # Returns
    /// String slice containing the file path (e.g., "Reactbase.json")
    pub fn reactbase_path(&self) -> &str {
        &self.config.reactbase
    }

    /// Returns the current path to the reaction dictionary JSON file.
    ///
    /// # Returns
    /// String slice containing the file path (e.g., "dict_reaction.json")
    pub fn dict_reaction_path(&self) -> &str {
        &self.config.dict_reaction
    }

    /// Returns the current path to the elements JSON file.
    ///
    /// # Returns
    /// String slice containing the file path (e.g., "elements.json")
    pub fn elements_path(&self) -> &str {
        &self.config.elements
    }

    /// Returns the current path to the problems folder.
    ///
    /// # Returns
    /// String slice containing the folder path (e.g., "problems")
    pub fn problems_folder_path(&self) -> &str {
        &self.config.problems_folder
    }

    /// Returns the actual filename for a given library type.
    ///
    /// # Arguments
    /// * `library_type` - The type of library to get the filename for
    ///
    /// # Returns
    /// String slice containing the actual filename
    pub fn get_library_filename(&self, library_type: LibraryType) -> &str {
        match library_type {
            LibraryType::SubstanceBase => &self.config.substance_base,
            LibraryType::AllKeysSubstance => &self.config.all_keys_substance,
            LibraryType::Elements => &self.config.elements,
            LibraryType::Reactbase => &self.config.reactbase,
            LibraryType::DictReaction => &self.config.dict_reaction,
        }
    }

    /// Updates the substance base library file path.
    ///
    /// Validates that the new file exists before updating the configuration.
    /// Automatically saves the configuration after successful update.
    ///
    /// # Arguments
    /// * `path` - New file path for the substance base library
    ///
    /// # Returns
    /// * `Ok(())` - If file exists and update was successful
    /// * `Err(Box<dyn std::error::Error>)` - If file doesn't exist or save failed
    pub fn set_substance_base(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        if Path::new(path).exists() {
            self.config.substance_base = path.to_string();
            self.save_config()?;
            Ok(())
        } else {
            Err(format!("File does not exist: {}", path).into())
        }
    }

    /// Updates the substance keys mapping file path.
    ///
    /// Validates that the new file exists before updating the configuration.
    /// Automatically saves the configuration after successful update.
    ///
    /// # Arguments
    /// * `path` - New file path for the substance keys library
    ///
    /// # Returns
    /// * `Ok(())` - If file exists and update was successful
    /// * `Err(Box<dyn std::error::Error>)` - If file doesn't exist or save failed
    pub fn set_all_keys_substance(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        if Path::new(path).exists() {
            self.config.all_keys_substance = path.to_string();
            self.save_config()?;
            Ok(())
        } else {
            Err(format!("File does not exist: {}", path).into())
        }
    }

    /// Updates the reaction base library file path.
    ///
    /// Validates that the new file exists before updating the configuration.
    /// Automatically saves the configuration after successful update.
    ///
    /// # Arguments
    /// * `path` - New file path for the reaction base library
    ///
    /// # Returns
    /// * `Ok(())` - If file exists and update was successful
    /// * `Err(Box<dyn std::error::Error>)` - If file doesn't exist or save failed
    pub fn set_reactbase(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        if Path::new(path).exists() {
            self.config.reactbase = path.to_string();
            self.save_config()?;
            Ok(())
        } else {
            Err(format!("File does not exist: {}", path).into())
        }
    }

    /// Updates the elements library file path.
    ///
    /// Validates that the new file exists before updating the configuration.
    /// Automatically saves the configuration after successful update.
    ///
    /// # Arguments
    /// * `path` - New file path for the elements library
    ///
    /// # Returns
    /// * `Ok(())` - If file exists and update was successful
    /// * `Err(Box<dyn std::error::Error>)` - If file doesn't exist or save failed
    pub fn set_elements(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        if Path::new(path).exists() {
            self.config.elements = path.to_string();
            self.save_config()?;
            Ok(())
        } else {
            Err(format!("File does not exist: {}", path).into())
        }
    }

    /// Updates the reaction dictionary file path.
    ///
    /// Validates that the new file exists before updating the configuration.
    /// Automatically saves the configuration after successful update.
    ///
    /// # Arguments
    /// * `path` - New file path for the reaction dictionary
    ///
    /// # Returns
    /// * `Ok(())` - If file exists and update was successful
    /// * `Err(Box<dyn std::error::Error>)` - If file doesn't exist or save failed
    pub fn set_dict_reaction(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        if Path::new(path).exists() {
            self.config.dict_reaction = path.to_string();
            self.save_config()?;
            Ok(())
        } else {
            Err(format!("File does not exist: {}", path).into())
        }
    }

    /// Updates multiple library paths in a single atomic operation.
    ///
    /// This method validates all files exist before making any changes, ensuring
    /// either all updates succeed or none are applied. More efficient than multiple
    /// individual setter calls.
    ///
    /// # Arguments
    /// * `updates` - HashMap mapping internal keys to new file paths
    ///   - "substance_base" -> path to substance base file
    ///   - "all_keys_substance" -> path to substance keys file
    ///   - "reactbase" -> path to reaction base file  
    ///   - "dict_reaction" -> path to reaction dictionary file
    ///
    /// # Returns
    /// * `Ok(())` - If all files exist and updates were successful
    /// * `Err(Box<dyn std::error::Error>)` - If any file doesn't exist or save failed
    ///
    /// # Example
    /// ```rust, ignore
    /// let mut updates = HashMap::new();
    /// updates.insert("substance_base", "substance_base_v3.json");
    /// updates.insert("reactbase", "Reactbase_v2.json");
    /// manager.update_libraries(updates)?;
    /// ```
    pub fn update_libraries(
        &mut self,
        updates: HashMap<&str, &str>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        for (key, path) in &updates {
            if !Path::new(path).exists() {
                return Err(format!("File does not exist: {}", path).into());
            }
        }

        for (key, path) in updates.clone() {
            match key {
                "substance_base" => self.config.substance_base = path.to_string(),
                "all_keys_substance" => self.config.all_keys_substance = path.to_string(),
                "elements" => self.config.elements = path.to_string(),
                "reactbase" => self.config.reactbase = path.to_string(),
                "dict_reaction" => self.config.dict_reaction = path.to_string(),
                "problems_folder" => {
                    if !Path::new(path).exists() {
                        fs::create_dir_all(path)?;
                    }
                    self.config.problems_folder = path.to_string();
                }
                _ => return Err(format!("Unknown library key: {}", key).into()),
            }
        }

        self.save_config()?;
        Ok(())
    }

    /// Returns a reference to the current library configuration.
    ///
    /// # Returns
    /// Reference to the LibraryConfig containing all current file paths
    pub fn get_config(&self) -> &LibraryConfig {
        &self.config
    }

    /// Updates the problems folder path.
    ///
    /// Creates the folder if it doesn't exist. Automatically saves the configuration after successful update.
    ///
    /// # Arguments
    /// * `path` - New folder path for the problems directory
    ///
    /// # Returns
    /// * `Ok(())` - If folder was created/exists and update was successful
    /// * `Err(Box<dyn std::error::Error>)` - If folder creation or save failed
    pub fn set_problems_folder(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        if !Path::new(path).exists() {
            fs::create_dir_all(path)?;
        }
        self.config.problems_folder = path.to_string();
        self.save_config()?;
        Ok(())
    }

    /// Ensures the problems folder exists, creating it if necessary.
    ///
    /// # Returns
    /// * `Ok(())` - If folder exists or was created successfully
    /// * `Err(Box<dyn std::error::Error>)` - If folder creation failed
    pub fn ensure_problems_folder(&self) -> Result<(), Box<dyn std::error::Error>> {
        if !Path::new(&self.config.problems_folder).exists() {
            fs::create_dir_all(&self.config.problems_folder)?;
        }
        Ok(())
    }

    /// Resets all library paths to their default values.
    ///
    /// Restores the configuration to default file names and saves the changes.
    ///
    /// # Returns
    /// * `Ok(())` - If reset and save were successful
    /// * `Err(Box<dyn std::error::Error>)` - If save operation failed
    pub fn reset_to_defaults(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        self.config = LibraryConfig::default();
        self.save_config()?;
        Ok(())
    }

    /// Logs library update operations to a log file.
    ///
    /// This helper function records when data is saved to library files,
    /// creating a persistent log of all library modifications.
    ///
    /// # Arguments
    /// * `library_type` - The type of library being updated
    /// * `content` - Description of what was saved
    ///
    /// # Returns
    /// * `Ok(())` - If log entry was written successfully
    /// * `Err(Box<dyn std::error::Error>)` - If log file write failed
    pub fn log_library_update(
        &self,
        library_type: LibraryType,
        content: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let log_file = "library_log.txt";
        let timestamp: DateTime<Local> = Local::now();
        let date_str = timestamp.format("%d.%m.%Y").to_string();

        let actual_filename = match library_type {
            LibraryType::SubstanceBase => &self.config.substance_base,
            LibraryType::AllKeysSubstance => &self.config.all_keys_substance,
            LibraryType::Elements => &self.config.elements,
            LibraryType::Reactbase => &self.config.reactbase,
            LibraryType::DictReaction => &self.config.dict_reaction,
        };

        let log_entry = format!(
            "{} in {} \"{}\"
",
            date_str, actual_filename, content
        );

        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(log_file)?;

        file.write_all(log_entry.as_bytes())?;
        Ok(())
    }

    /// Adds a new top-level key with empty value to a JSON library file.
    ///
    /// # Arguments
    /// * `library_type` - The type of library to modify
    /// * `key` - The top-level key to add
    ///
    /// # Returns
    /// * `Ok(())` - If key was added successfully
    /// * `Err(Box<dyn std::error::Error>)` - If file operation failed
    pub fn add_library_key(
        &self,
        library_type: LibraryType,
        key: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let filename = self.get_library_filename(library_type);

        let mut json_data: serde_json::Value = if Path::new(filename).exists() {
            let content = fs::read_to_string(filename)?;
            serde_json::from_str(&content)?
        } else {
            serde_json::json!({})
        };

        if let Some(obj) = json_data.as_object_mut() {
            obj.insert(key.to_string(), serde_json::json!({}));
        }

        let content = serde_json::to_string_pretty(&json_data)?;
        fs::write(filename, content)?;

        Ok(())
    }

    /// Adds a secondary key-value pair to an existing primary key in a JSON library file.
    ///
    /// # Arguments
    /// * `library_type` - The type of library to modify
    /// * `primary_key` - The top-level key to find
    /// * `secondary_key` - The key to add inside the primary key
    /// * `value` - The value to associate with the secondary key
    ///
    /// # Returns
    /// * `Ok(())` - If key-value pair was added successfully
    /// * `Err(Box<dyn std::error::Error>)` - If file operation failed or primary key not found
    pub fn add_secondary_key(
        &self,
        library_type: LibraryType,
        primary_key: &str,
        secondary_key: &str,
        value: serde_json::Value,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let filename = self.get_library_filename(library_type);

        let mut json_data: serde_json::Value = if Path::new(filename).exists() {
            let content = fs::read_to_string(filename)?;
            serde_json::from_str(&content)?
        } else {
            return Err(format!("Library file does not exist: {}", filename).into());
        };

        if let Some(obj) = json_data.as_object_mut() {
            if let Some(primary_obj) = obj.get_mut(primary_key) {
                if let Some(primary_map) = primary_obj.as_object_mut() {
                    primary_map.insert(secondary_key.to_string(), value);
                } else {
                    return Err(format!("Primary key '{}' is not an object", primary_key).into());
                }
            } else {
                return Err(format!("Primary key '{}' not found", primary_key).into());
            }
        }

        let content = serde_json::to_string_pretty(&json_data)?;
        fs::write(filename, content)?;

        Ok(())
    }
}

/// Global singleton instance of LibraryManager using thread-safe OnceLock pattern
static GLOBAL_LIBRARY_MANAGER: OnceLock<Mutex<LibraryManager>> = OnceLock::new();

/// Test-specific manager instance to isolate tests from global state
#[cfg(test)]
static TEST_MANAGER: std::sync::Mutex<Option<LibraryManager>> = std::sync::Mutex::new(None);

/// Sets a test-specific manager instance for test isolation.
///
/// This function is only available during testing and allows tests to use
/// their own LibraryManager instance without affecting the global state.
///
/// # Arguments
/// * `manager` - LibraryManager instance to use for tests
#[cfg(test)]
pub fn set_test_manager(manager: LibraryManager) {
    *TEST_MANAGER.lock().unwrap() = Some(manager);
}

/// Clears the test-specific manager instance.
///
/// This function should be called at the end of each test to ensure
/// test isolation and prevent state leakage between tests.
#[cfg(test)]
pub fn clear_test_manager() {
    *TEST_MANAGER.lock().unwrap() = None;
}

/// Returns a mutex guard to the global LibraryManager instance.
///
/// This function provides thread-safe access to the singleton LibraryManager.
/// During tests, it uses the test-specific manager if one has been set.
///
/// # Returns
/// MutexGuard providing exclusive access to the LibraryManager
///
/// # Panics
/// Panics if the mutex is poisoned (should not happen in normal operation)
pub fn get_library_manager() -> std::sync::MutexGuard<'static, LibraryManager> {
    #[cfg(test)]
    {
        if let Some(ref manager) = *TEST_MANAGER.lock().unwrap() {
            // For tests, we need to return a different type, so we'll use the global one but reset it
            let _ = GLOBAL_LIBRARY_MANAGER.set(Mutex::new(manager.clone()));
        }
    }

    GLOBAL_LIBRARY_MANAGER
        .get_or_init(|| Mutex::new(LibraryManager::new()))
        .lock()
        .unwrap()
}

/// Executes a closure with read-only access to the LibraryManager.
///
/// This function provides a convenient way to access the LibraryManager
/// for read operations without needing to manage the mutex guard directly.
///
/// # Arguments
/// * `f` - Closure that receives a reference to LibraryManager
///
/// # Returns
/// The result of the closure execution
///
/// # Example
/// ```rust
/// let path = with_library_manager(|manager| {
///     manager.substance_base_path().to_string()
/// });
/// ```
pub fn with_library_manager<F, R>(f: F) -> R
where
    F: FnOnce(&LibraryManager) -> R,
{
    let manager = get_library_manager();
    f(&*manager)
}

/// Executes a closure with mutable access to the LibraryManager.
///
/// This function provides a convenient way to access the LibraryManager
/// for write operations without needing to manage the mutex guard directly.
///
/// # Arguments
/// * `f` - Closure that receives a mutable reference to LibraryManager
///
/// # Returns
/// The result of the closure execution
///
/// # Example
/// ```rust
/// let result = with_library_manager_mut(|manager| {
///     manager.set_substance_base("substance_base_v3.json")
/// });
/// ```
pub fn with_library_manager_mut<F, R>(f: F) -> R
where
    F: FnOnce(&mut LibraryManager) -> R,
{
    let mut manager = get_library_manager();
    f(&mut *manager)
}

/// Convenience function to log library updates.
///
/// This function provides easy access to the logging functionality
/// for library update operations.
///
/// # Arguments
/// * `library_type` - The type of library being updated
/// * `content` - Description of what was saved
///
/// # Returns
/// * `Ok(())` - If log entry was written successfully
/// * `Err(Box<dyn std::error::Error>)` - If log file write failed
///
/// # Example
/// ```rust
/// log_library_update(LibraryType::SubstanceBase, "Added CH4 thermodynamic data")?;
/// ```
pub fn log_library_update(
    library_type: LibraryType,
    content: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    with_library_manager(|manager| manager.log_library_update(library_type, content))
}

/// Returns the actual filename for a given library type.
///
/// # Arguments
/// * `library_type` - The type of library to get the filename for
///
/// # Returns
/// String containing the actual filename
///
/// # Example
/// ```rust
/// let filename = get_library_filename(LibraryType::SubstanceBase);
/// // Returns "substance_base_v2.json"
/// ```
pub fn get_library_filename(library_type: LibraryType) -> String {
    with_library_manager(|manager| manager.get_library_filename(library_type).to_string())
}

/// Adds a new top-level key with empty value to a JSON library file.
///
/// # Arguments
/// * `library_type` - The type of library to modify
/// * `key` - The top-level key to add
///
/// # Returns
/// * `Ok(())` - If key was added successfully
/// * `Err(Box<dyn std::error::Error>)` - If file operation failed
///
/// # Example
/// ```rust
/// add_library_key(LibraryType::SubstanceBase, "CH4")?;
/// // Adds "CH4": {} to the substance base JSON file
/// ```
pub fn add_library_key(
    library_type: LibraryType,
    key: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    with_library_manager(|manager| manager.add_library_key(library_type, key))
}

/// Adds a secondary key-value pair to an existing primary key in a JSON library file.
///
/// # Arguments
/// * `library_type` - The type of library to modify
/// * `primary_key` - The top-level key to find
/// * `secondary_key` - The key to add inside the primary key
/// * `value` - The value to associate with the secondary key
///
/// # Returns
/// * `Ok(())` - If key-value pair was added successfully
/// * `Err(Box<dyn std::error::Error>)` - If file operation failed or primary key not found
///
/// # Example
/// ```rust
/// add_secondary_key(LibraryType::SubstanceBase, "CH4", "formula", serde_json::json!("CH4"))?;
/// // Adds "formula": "CH4" inside the CH4 object
/// ```
pub fn add_secondary_key(
    library_type: LibraryType,
    primary_key: &str,
    secondary_key: &str,
    value: serde_json::Value,
) -> Result<(), Box<dyn std::error::Error>> {
    with_library_manager(|manager| {
        manager.add_secondary_key(library_type, primary_key, secondary_key, value)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn setup_test_files() -> (NamedTempFile, NamedTempFile, NamedTempFile, NamedTempFile) {
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
    fn test_library_manager_new() {
        let manager = LibraryManager::new();
        assert_eq!(manager.all_keys_substance_path(), "all_keys_substance.json");
        assert_eq!(manager.substance_base_path(), "substance_base_v2.json");

        assert_eq!(manager.reactbase_path(), "Reactbase.json");
        assert_eq!(manager.dict_reaction_path(), "dict_reaction.json");
    }

    #[test]
    fn test_library_manager_with_config() {
        let mut temp_config = NamedTempFile::new().unwrap();
        let mut temp_substance = NamedTempFile::new().unwrap();
        let mut temp_keys = NamedTempFile::new().unwrap();
        let mut temp_react = NamedTempFile::new().unwrap();
        let mut temp_dict = NamedTempFile::new().unwrap();

        temp_substance.write_all(b"{}").unwrap();
        temp_keys.write_all(b"[]").unwrap();
        temp_react.write_all(b"{}").unwrap();
        temp_dict.write_all(b"{}").unwrap();

        let config = LibraryConfig {
            substance_base: temp_substance.path().to_str().unwrap().to_string(),
            all_keys_substance: temp_keys.path().to_str().unwrap().to_string(),
            elements: temp_keys.path().to_str().unwrap().to_string(),
            reactbase: temp_react.path().to_str().unwrap().to_string(),
            dict_reaction: temp_dict.path().to_str().unwrap().to_string(),
            problems_folder: "problems".to_string(),
        };

        let config_json = serde_json::to_string_pretty(&config).unwrap();
        temp_config.write_all(config_json.as_bytes()).unwrap();

        let manager = LibraryManager::with_config_file(temp_config.path().to_str().unwrap());
        assert_eq!(
            manager.substance_base_path(),
            temp_substance.path().to_str().unwrap()
        );
        assert_eq!(
            manager.all_keys_substance_path(),
            temp_keys.path().to_str().unwrap()
        );
    }

    #[test]
    fn test_update_libraries() {
        let mut temp_config = NamedTempFile::new().unwrap();
        let mut temp_substance = NamedTempFile::new().unwrap();
        let mut temp_keys = NamedTempFile::new().unwrap();

        // Create dummy files
        temp_substance.write_all(b"{}").unwrap();
        temp_keys.write_all(b"[]").unwrap();

        let mut manager = LibraryManager::with_config_file(temp_config.path().to_str().unwrap());

        let mut updates = HashMap::new();
        updates.insert("substance_base", temp_substance.path().to_str().unwrap());
        updates.insert("all_keys_substance", temp_keys.path().to_str().unwrap());

        let result = manager.update_libraries(updates);
        assert!(result.is_ok());

        assert_eq!(
            manager.substance_base_path(),
            temp_substance.path().to_str().unwrap()
        );
        assert_eq!(
            manager.all_keys_substance_path(),
            temp_keys.path().to_str().unwrap()
        );
    }
}
