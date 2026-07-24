use crate::library_manager::LibraryManager;
use chrono::Local;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io;
use std::path::Path;
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ThermoLibraryDropReport {
    pub library: String,
    pub substance: String,
    pub removed_from_substance_base: bool,
    pub removed_from_all_keys: usize,
    pub removed_from_elements: usize,
}
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ThermoLibraryWriteReport {
    pub library: String,
    pub substance: String,
    pub replaced_existing: bool,
    pub refreshed_all_keys: usize,
    pub refreshed_elements: usize,
}
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ThermoLibraryWriteOperation {
    Upsert,
    Drop,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[allow(dead_code)]
pub(crate) struct ThermoLibraryJournalEntry {
    pub timestamp: String,
    #[serde(default)]
    pub source_format: String,
    pub operation: ThermoLibraryWriteOperation,
    pub library: String,
    pub substance: String,
    pub canonical_names: Vec<String>,
    pub files: Vec<String>,
}

/// Typed front-door description for writing one substance into a thermo library.
///
/// The payload stays structured so the higher-level API does not need to pass raw
/// JSON blobs around, while the lower-level writer still stores the final record
/// as JSON in the canonical substance base file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[allow(dead_code)]
pub(crate) struct ThermoSubstanceWriteSpec {
    pub library: String,
    pub substance: String,
    pub source_format: ThermoSubstanceSourceFormat,
    pub composition: String,
    #[serde(default)]
    pub fields: BTreeMap<String, serde_json::Value>,
}

impl ThermoSubstanceWriteSpec {
    /// Converts the structured request into the JSON payload stored in the
    /// canonical substance base file.
    #[allow(dead_code)]
    fn into_record(self) -> serde_json::Value {
        let mut payload = serde_json::Map::new();
        payload.insert(
            "composition".to_string(),
            serde_json::Value::String(self.composition),
        );
        for (key, value) in self.fields {
            payload.insert(key, value);
        }
        serde_json::Value::Object(payload)
    }
}

/// Declares the broad source format for a typed thermo write request.
///
/// This is intentionally narrow: the user-facing write workflow should choose
/// from the known families instead of inventing new spellings.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ThermoSubstanceSourceFormat {
    Nasa,
    Nist,
    Aramco,
    Cea,
    Custom(String),
}

impl ThermoSubstanceSourceFormat {
    #[allow(dead_code)]
    fn as_str(&self) -> &str {
        match self {
            Self::Nasa => "NASA",
            Self::Nist => "NIST",
            Self::Aramco => "Aramco",
            Self::Cea => "CEA",
            Self::Custom(value) => value.as_str(),
        }
    }
    #[allow(dead_code)]
    fn expected_library_prefix(&self) -> Option<&str> {
        match self {
            Self::Nasa => Some("NASA"),
            Self::Nist => Some("NIST"),
            Self::Aramco => Some("Aramco"),
            Self::Cea => Some("CEA"),
            Self::Custom(_) => None,
        }
    }
}
#[allow(dead_code)]
#[derive(Debug)]
pub(crate) enum ThermoLibraryWriteError {
    MissingLibraryFile {
        path: String,
        source: io::Error,
    },
    InvalidLibraryJson {
        path: String,
        source: serde_json::Error,
    },
    MissingLibraryEntry {
        library: String,
    },
    MissingSubstanceEntry {
        library: String,
        substance: String,
    },
    LibraryFormatMismatch {
        library: String,
        source_format: String,
    },
    WriteFailed {
        path: String,
        source: io::Error,
    },
}

impl std::fmt::Display for ThermoLibraryWriteError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingLibraryFile { path, source } => {
                write!(
                    f,
                    "failed to open thermo library file '{}': {}",
                    path, source
                )
            }
            Self::InvalidLibraryJson { path, source } => {
                write!(
                    f,
                    "failed to parse thermo library file '{}': {}",
                    path, source
                )
            }
            Self::MissingLibraryEntry { library } => {
                write!(f, "library '{}' was not found in substance base", library)
            }
            Self::MissingSubstanceEntry { library, substance } => {
                write!(
                    f,
                    "substance '{}' was not found in library '{}'",
                    substance, library
                )
            }
            Self::LibraryFormatMismatch {
                library,
                source_format,
            } => {
                write!(
                    f,
                    "library '{}' does not match typed source format '{}'",
                    library, source_format
                )
            }
            Self::WriteFailed { path, source } => {
                write!(
                    f,
                    "failed to write thermo library file '{}': {}",
                    path, source
                )
            }
        }
    }
}

impl std::error::Error for ThermoLibraryWriteError {}
#[allow(dead_code)]
fn read_json_file<T: serde::de::DeserializeOwned>(
    path: &str,
) -> Result<T, ThermoLibraryWriteError> {
    let content =
        fs::read_to_string(path).map_err(|source| ThermoLibraryWriteError::MissingLibraryFile {
            path: path.to_string(),
            source,
        })?;
    serde_json::from_str(&content).map_err(|source| ThermoLibraryWriteError::InvalidLibraryJson {
        path: path.to_string(),
        source,
    })
}
#[allow(dead_code)]
fn write_json_file<T: serde::Serialize>(
    path: &str,
    value: &T,
) -> Result<(), ThermoLibraryWriteError> {
    let json = serde_json::to_string_pretty(value).map_err(|source| {
        ThermoLibraryWriteError::InvalidLibraryJson {
            path: path.to_string(),
            source,
        }
    })?;
    fs::write(path, json).map_err(|source| ThermoLibraryWriteError::WriteFailed {
        path: path.to_string(),
        source,
    })
}
#[allow(dead_code)]
fn append_journal_entry(
    path: &str,
    entry: ThermoLibraryJournalEntry,
) -> Result<(), ThermoLibraryWriteError> {
    let mut entries: Vec<ThermoLibraryJournalEntry> = if Path::new(path).exists() {
        read_json_file(path)?
    } else {
        Vec::new()
    };
    entries.push(entry);
    write_json_file(path, &entries)
}

impl LibraryManager {
    /// Inserts or replaces one substance in one concrete library and keeps the
    /// derived catalog files synchronized.
    ///
    /// This method is crate-visible so external consumers cannot mutate the
    /// on-disk library graph in an ad hoc way.
    #[allow(dead_code)]
    pub(crate) fn upsert_substance_in_library(
        &self,
        library: &str,
        substance: &str,
        substance_data: serde_json::Value,
    ) -> Result<ThermoLibraryWriteReport, ThermoLibraryWriteError> {
        let substance_base_path = self.substance_base_path().to_string();
        let all_keys_path = self.all_keys_substance_path().to_string();
        let elements_path = self.elements_path().to_string();
        let journal_path = self.thermo_write_journal_path().to_string();

        let mut substance_base: HashMap<String, HashMap<String, serde_json::Value>> =
            read_json_file(&substance_base_path)?;
        let library_bucket = substance_base
            .entry(library.to_string())
            .or_insert_with(HashMap::new);
        let replaced_existing = library_bucket
            .insert(substance.to_string(), substance_data)
            .is_some();
        write_json_file(&substance_base_path, &substance_base)?;

        // Regenerate the derived catalogs from the canonical substance base so
        // the address and element layers remain consistent after every write.
        self.recreate_all_keys_substance().map_err(|source| {
            ThermoLibraryWriteError::WriteFailed {
                path: all_keys_path.clone(),
                source: io::Error::new(io::ErrorKind::Other, source.to_string()),
            }
        })?;
        self.create_elements_base()
            .map_err(|source| ThermoLibraryWriteError::WriteFailed {
                path: elements_path.clone(),
                source: io::Error::new(io::ErrorKind::Other, source.to_string()),
            })?;

        append_journal_entry(
            &journal_path,
            ThermoLibraryJournalEntry {
                timestamp: Local::now().to_rfc3339(),
                operation: ThermoLibraryWriteOperation::Upsert,
                source_format: "manual".to_string(),
                library: library.to_string(),
                substance: substance.to_string(),
                canonical_names: vec![substance.to_string()],
                files: vec![
                    substance_base_path.clone(),
                    all_keys_path.clone(),
                    elements_path.clone(),
                ],
            },
        )?;

        Ok(ThermoLibraryWriteReport {
            library: library.to_string(),
            substance: substance.to_string(),
            replaced_existing,
            refreshed_all_keys: 1,
            refreshed_elements: 1,
        })
    }

    /// Writes one substance through the typed thermo write workflow.
    ///
    /// This is the higher-level entry point intended for future UI and script
    /// callers. It keeps the request structured until the last moment and only
    /// materializes the JSON payload right before persistence.
    #[allow(dead_code)]
    pub(crate) fn write_substance_from_spec(
        &self,
        spec: ThermoSubstanceWriteSpec,
    ) -> Result<ThermoLibraryWriteReport, ThermoLibraryWriteError> {
        if let Some(expected_prefix) = spec.source_format.expected_library_prefix() {
            if !spec.library.starts_with(expected_prefix) {
                return Err(ThermoLibraryWriteError::LibraryFormatMismatch {
                    library: spec.library,
                    source_format: spec.source_format.as_str().to_string(),
                });
            }
        }

        let record = spec.clone().into_record();
        let library = spec.library;
        let substance = spec.substance;
        let source_format = spec.source_format.as_str().to_string();
        let report = self.upsert_substance_in_library(&library, &substance, record)?;

        // Update the last journal entry so typed writes retain the format family
        // that was used to materialize the JSON record.
        let journal_path = self.thermo_write_journal_path().to_string();
        if Path::new(&journal_path).exists() {
            let mut entries: Vec<ThermoLibraryJournalEntry> = read_json_file(&journal_path)?;
            if let Some(last) = entries.last_mut() {
                last.source_format = source_format;
                write_json_file(&journal_path, &entries)?;
            }
        }

        Ok(report)
    }

    /// Removes one substance from one concrete library and keeps derived catalog files in sync.
    ///
    /// This is crate-visible on purpose: dropping library data is a structural
    /// mutation and should not be available to external consumers by accident.
    #[allow(dead_code)]
    pub(crate) fn drop_substance_from_library(
        &self,
        library: &str,
        substance: &str,
    ) -> Result<ThermoLibraryDropReport, ThermoLibraryWriteError> {
        let substance_base_path = self.substance_base_path().to_string();
        let all_keys_path = self.all_keys_substance_path().to_string();
        let elements_path = self.elements_path().to_string();
        let journal_path = self.thermo_write_journal_path().to_string();

        let mut substance_base: HashMap<String, HashMap<String, serde_json::Value>> =
            read_json_file(&substance_base_path)?;
        let Some(library_bucket) = substance_base.get_mut(library) else {
            return Err(ThermoLibraryWriteError::MissingLibraryEntry {
                library: library.to_string(),
            });
        };
        if library_bucket.remove(substance).is_none() {
            return Err(ThermoLibraryWriteError::MissingSubstanceEntry {
                library: library.to_string(),
                substance: substance.to_string(),
            });
        }
        if library_bucket.is_empty() {
            substance_base.remove(library);
        }
        write_json_file(&substance_base_path, &substance_base)?;

        let mut all_keys: Vec<(String, String)> = read_json_file(&all_keys_path)?;
        let removed_from_all_keys_before = all_keys.len();
        all_keys.retain(|(lib, sub)| !(lib == library && sub == substance));
        let removed_from_all_keys = removed_from_all_keys_before.saturating_sub(all_keys.len());
        write_json_file(&all_keys_path, &all_keys)?;

        let mut elements: HashMap<String, Vec<Vec<String>>> = read_json_file(&elements_path)?;
        let mut removed_from_elements = 0usize;
        let element_keys: Vec<String> = elements.keys().cloned().collect();
        for element in element_keys {
            if let Some(entries) = elements.get_mut(&element) {
                let before = entries.len();
                entries.retain(|pair| {
                    !(pair.len() >= 2 && pair[0] == substance && pair[1] == library)
                });
                removed_from_elements += before.saturating_sub(entries.len());
                if entries.is_empty() {
                    elements.remove(&element);
                }
            }
        }
        write_json_file(&elements_path, &elements)?;

        // Regenerate derived catalogs from the updated substance base so the
        // address and element libraries always reflect the canonical source.
        self.recreate_all_keys_substance().map_err(|source| {
            ThermoLibraryWriteError::WriteFailed {
                path: all_keys_path.clone(),
                source: io::Error::new(io::ErrorKind::Other, source.to_string()),
            }
        })?;
        self.create_elements_base()
            .map_err(|source| ThermoLibraryWriteError::WriteFailed {
                path: elements_path.clone(),
                source: io::Error::new(io::ErrorKind::Other, source.to_string()),
            })?;

        append_journal_entry(
            &journal_path,
            ThermoLibraryJournalEntry {
                timestamp: Local::now().to_rfc3339(),
                source_format: "manual".to_string(),
                operation: ThermoLibraryWriteOperation::Drop,
                library: library.to_string(),
                substance: substance.to_string(),
                canonical_names: vec![substance.to_string()],
                files: vec![
                    substance_base_path.clone(),
                    all_keys_path.clone(),
                    elements_path.clone(),
                ],
            },
        )?;

        Ok(ThermoLibraryDropReport {
            library: library.to_string(),
            substance: substance.to_string(),
            removed_from_substance_base: true,
            removed_from_all_keys,
            removed_from_elements,
        })
    }

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

        all_keys.sort();

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

        for entries in elements_map.values_mut() {
            entries.sort();
        }

        let elements_json = serde_json::to_string_pretty(&elements_map)?;
        fs::write(elements_path, elements_json)?;

        Ok(())
    }
}

#[cfg(test)]
mod write_workflow_tests {
    use super::*;
    use serde_json::json;
    use tempfile::tempdir;

    #[test]
    fn drop_substance_updates_all_catalog_files() {
        let dir = tempdir().expect("tempdir should be created");
        let config_path = dir.path().join("library_config.json");
        let substance_base_path = dir.path().join("substance_base.json");
        let all_keys_path = dir.path().join("all_keys_substance.json");
        let elements_path = dir.path().join("elements.json");
        let journal_path = dir.path().join("thermo_write_journal.json");

        let substance_base = json!({
            "NASA_gas": {
                "CO": { "composition": "{C:1, O:1}", "Cp": [1.0] },
                "H2O": { "composition": "{H:2, O:1}", "Cp": [2.0] }
            },
            "NIST": {
                "CO": { "composition": "{C:1, O:1}", "Cp": [3.0] }
            }
        });
        let all_keys = json!([["NASA_gas", "CO"], ["NASA_gas", "H2O"], ["NIST", "CO"]]);
        let elements = json!({
            "C": [["CO", "NASA_gas"], ["CO", "NIST"]],
            "H": [["H2O", "NASA_gas"]],
            "O": [["CO", "NASA_gas"], ["H2O", "NASA_gas"], ["CO", "NIST"]]
        });

        fs::write(
            &substance_base_path,
            serde_json::to_string_pretty(&substance_base).unwrap(),
        )
        .unwrap();
        fs::write(
            &all_keys_path,
            serde_json::to_string_pretty(&all_keys).unwrap(),
        )
        .unwrap();
        fs::write(
            &elements_path,
            serde_json::to_string_pretty(&elements).unwrap(),
        )
        .unwrap();
        fs::write(&journal_path, "[]").unwrap();

        let mut manager = LibraryManager::with_config_file(config_path.to_str().unwrap());
        manager
            .set_substance_base(substance_base_path.to_str().unwrap())
            .unwrap();
        manager
            .set_all_keys_substance(all_keys_path.to_str().unwrap())
            .unwrap();
        manager
            .set_elements(elements_path.to_str().unwrap())
            .unwrap();
        manager
            .set_thermo_write_journal(journal_path.to_str().unwrap())
            .unwrap();

        let report = manager
            .drop_substance_from_library("NASA_gas", "CO")
            .expect("drop should succeed");

        assert_eq!(report.library, "NASA_gas");
        assert_eq!(report.substance, "CO");
        assert!(report.removed_from_substance_base);
        assert_eq!(report.removed_from_all_keys, 1);
        assert_eq!(report.removed_from_elements, 2);

        let updated_base: HashMap<String, HashMap<String, serde_json::Value>> =
            serde_json::from_str(&fs::read_to_string(&substance_base_path).unwrap()).unwrap();
        assert!(!updated_base.get("NASA_gas").unwrap().contains_key("CO"));
        assert!(updated_base.get("NASA_gas").unwrap().contains_key("H2O"));
        assert!(updated_base.get("NIST").unwrap().contains_key("CO"));

        let updated_all_keys: Vec<(String, String)> =
            serde_json::from_str(&fs::read_to_string(&all_keys_path).unwrap()).unwrap();
        assert_eq!(
            updated_all_keys,
            vec![
                ("NASA_gas".to_string(), "H2O".to_string()),
                ("NIST".to_string(), "CO".to_string())
            ]
        );

        let updated_elements: HashMap<String, Vec<Vec<String>>> =
            serde_json::from_str(&fs::read_to_string(&elements_path).unwrap()).unwrap();
        assert_eq!(
            updated_elements.get("C").unwrap(),
            &vec![vec!["CO".to_string(), "NIST".to_string()]]
        );
        assert_eq!(
            updated_elements.get("H").unwrap(),
            &vec![vec!["H2O".to_string(), "NASA_gas".to_string()]]
        );
        assert_eq!(
            updated_elements.get("O").unwrap(),
            &vec![
                vec!["CO".to_string(), "NIST".to_string()],
                vec!["H2O".to_string(), "NASA_gas".to_string()]
            ]
        );

        let journal: Vec<ThermoLibraryJournalEntry> =
            serde_json::from_str(&fs::read_to_string(&journal_path).unwrap()).unwrap();
        assert_eq!(journal.len(), 1);
        assert_eq!(journal[0].operation, ThermoLibraryWriteOperation::Drop);
        assert_eq!(journal[0].library, "NASA_gas");
        assert_eq!(journal[0].substance, "CO");
        assert_eq!(journal[0].canonical_names, vec!["CO".to_string()]);
        assert_eq!(
            journal[0].files,
            vec![
                substance_base_path.to_str().unwrap().to_string(),
                all_keys_path.to_str().unwrap().to_string(),
                elements_path.to_str().unwrap().to_string()
            ]
        );
    }

    #[test]
    fn drop_substance_rejects_missing_pair() {
        let dir = tempdir().expect("tempdir should be created");
        let config_path = dir.path().join("library_config.json");
        let substance_base_path = dir.path().join("substance_base.json");
        let all_keys_path = dir.path().join("all_keys_substance.json");
        let elements_path = dir.path().join("elements.json");
        let journal_path = dir.path().join("thermo_write_journal.json");

        let substance_base = json!({
            "NASA_gas": {
                "H2O": { "composition": "{H:2, O:1}", "Cp": [2.0] }
            }
        });
        let all_keys = json!([["NASA_gas", "H2O"]]);
        let elements = json!({
            "H": [["H2O", "NASA_gas"]],
            "O": [["H2O", "NASA_gas"]]
        });

        fs::write(
            &substance_base_path,
            serde_json::to_string_pretty(&substance_base).unwrap(),
        )
        .unwrap();
        fs::write(
            &all_keys_path,
            serde_json::to_string_pretty(&all_keys).unwrap(),
        )
        .unwrap();
        fs::write(
            &elements_path,
            serde_json::to_string_pretty(&elements).unwrap(),
        )
        .unwrap();
        fs::write(&journal_path, "[]").unwrap();

        let mut manager = LibraryManager::with_config_file(config_path.to_str().unwrap());
        manager
            .set_substance_base(substance_base_path.to_str().unwrap())
            .unwrap();
        manager
            .set_all_keys_substance(all_keys_path.to_str().unwrap())
            .unwrap();
        manager
            .set_elements(elements_path.to_str().unwrap())
            .unwrap();
        manager
            .set_thermo_write_journal(journal_path.to_str().unwrap())
            .unwrap();

        let err = manager
            .drop_substance_from_library("NASA_gas", "CO")
            .expect_err("missing substance should fail");
        let message = err.to_string();
        assert!(message.contains("substance 'CO'"));
        assert!(message.contains("NASA_gas"));

        let journal: Vec<ThermoLibraryJournalEntry> =
            serde_json::from_str(&fs::read_to_string(&journal_path).unwrap()).unwrap();
        assert!(journal.is_empty());
    }

    #[test]
    fn upsert_substance_updates_all_catalog_files_without_touching_real_paths() {
        let dir = tempdir().expect("tempdir should be created");
        let config_path = dir.path().join("library_config.json");
        let substance_base_path = dir.path().join("substance_base.json");
        let all_keys_path = dir.path().join("all_keys_substance.json");
        let elements_path = dir.path().join("elements.json");
        let journal_path = dir.path().join("thermo_write_journal.json");

        let substance_base = json!({
            "NASA_gas": {
                "H2O": { "composition": "{H:2, O:1}", "Cp": [2.0] }
            },
            "NIST": {
                "CO": { "composition": "{C:1, O:1}", "Cp": [3.0] }
            }
        });
        let all_keys = json!([["NASA_gas", "H2O"], ["NIST", "CO"]]);
        let elements = json!({
            "C": [["CO", "NIST"]],
            "H": [["H2O", "NASA_gas"]],
            "O": [["CO", "NIST"], ["H2O", "NASA_gas"]]
        });

        fs::write(
            &substance_base_path,
            serde_json::to_string_pretty(&substance_base).unwrap(),
        )
        .unwrap();
        fs::write(
            &all_keys_path,
            serde_json::to_string_pretty(&all_keys).unwrap(),
        )
        .unwrap();
        fs::write(
            &elements_path,
            serde_json::to_string_pretty(&elements).unwrap(),
        )
        .unwrap();
        fs::write(&journal_path, "[]").unwrap();

        let mut manager = LibraryManager::with_config_file(config_path.to_str().unwrap());
        manager
            .set_substance_base(substance_base_path.to_str().unwrap())
            .unwrap();
        manager
            .set_all_keys_substance(all_keys_path.to_str().unwrap())
            .unwrap();
        manager
            .set_elements(elements_path.to_str().unwrap())
            .unwrap();
        manager
            .set_thermo_write_journal(journal_path.to_str().unwrap())
            .unwrap();

        let report = manager
            .write_substance_from_spec(ThermoSubstanceWriteSpec {
                library: "NASA_gas".to_string(),
                substance: "CH4".to_string(),
                source_format: ThermoSubstanceSourceFormat::Nasa,
                composition: "{C:1, H:4}".to_string(),
                fields: BTreeMap::from([("Cp".to_string(), json!([4.0]))]),
            })
            .expect("upsert should succeed");

        assert_eq!(report.library, "NASA_gas");
        assert_eq!(report.substance, "CH4");
        assert!(!report.replaced_existing);
        assert_eq!(report.refreshed_all_keys, 1);
        assert_eq!(report.refreshed_elements, 1);

        let updated_base: HashMap<String, HashMap<String, serde_json::Value>> =
            serde_json::from_str(&fs::read_to_string(&substance_base_path).unwrap()).unwrap();
        assert!(updated_base.get("NASA_gas").unwrap().contains_key("CH4"));
        assert!(updated_base.get("NASA_gas").unwrap().contains_key("H2O"));
        assert!(updated_base.get("NIST").unwrap().contains_key("CO"));

        let updated_all_keys: Vec<(String, String)> =
            serde_json::from_str(&fs::read_to_string(&all_keys_path).unwrap()).unwrap();
        assert_eq!(
            updated_all_keys,
            vec![
                ("NASA_gas".to_string(), "CH4".to_string()),
                ("NASA_gas".to_string(), "H2O".to_string()),
                ("NIST".to_string(), "CO".to_string())
            ]
        );

        let updated_elements: HashMap<String, Vec<Vec<String>>> =
            serde_json::from_str(&fs::read_to_string(&elements_path).unwrap()).unwrap();
        assert_eq!(
            updated_elements.get("C").unwrap(),
            &vec![
                vec!["CH4".to_string(), "NASA_gas".to_string()],
                vec!["CO".to_string(), "NIST".to_string()]
            ]
        );
        assert_eq!(
            updated_elements.get("H").unwrap(),
            &vec![
                vec!["CH4".to_string(), "NASA_gas".to_string()],
                vec!["H2O".to_string(), "NASA_gas".to_string()]
            ]
        );
        assert_eq!(
            updated_elements.get("O").unwrap(),
            &vec![
                vec!["CO".to_string(), "NIST".to_string()],
                vec!["H2O".to_string(), "NASA_gas".to_string()]
            ]
        );

        let journal: Vec<ThermoLibraryJournalEntry> =
            serde_json::from_str(&fs::read_to_string(&journal_path).unwrap()).unwrap();
        assert_eq!(journal.len(), 1);
        assert_eq!(journal[0].source_format, "NASA");
        assert_eq!(journal[0].operation, ThermoLibraryWriteOperation::Upsert);
        assert_eq!(journal[0].library, "NASA_gas");
        assert_eq!(journal[0].substance, "CH4");
        assert_eq!(journal[0].canonical_names, vec!["CH4".to_string()]);
    }

    #[test]
    fn typed_write_rejects_library_format_mismatch() {
        let dir = tempdir().expect("tempdir should be created");
        let config_path = dir.path().join("library_config.json");
        let substance_base_path = dir.path().join("substance_base.json");
        let all_keys_path = dir.path().join("all_keys_substance.json");
        let elements_path = dir.path().join("elements.json");
        let journal_path = dir.path().join("thermo_write_journal.json");

        let substance_base = json!({
            "NASA_gas": {}
        });

        fs::write(
            &substance_base_path,
            serde_json::to_string_pretty(&substance_base).unwrap(),
        )
        .unwrap();
        fs::write(&all_keys_path, "[]").unwrap();
        fs::write(&elements_path, "{}").unwrap();
        fs::write(&journal_path, "[]").unwrap();

        let mut manager = LibraryManager::with_config_file(config_path.to_str().unwrap());
        manager
            .set_substance_base(substance_base_path.to_str().unwrap())
            .unwrap();
        manager
            .set_all_keys_substance(all_keys_path.to_str().unwrap())
            .unwrap();
        manager
            .set_elements(elements_path.to_str().unwrap())
            .unwrap();
        manager
            .set_thermo_write_journal(journal_path.to_str().unwrap())
            .unwrap();

        let err = manager
            .write_substance_from_spec(ThermoSubstanceWriteSpec {
                library: "NASA_gas".to_string(),
                substance: "CO".to_string(),
                source_format: ThermoSubstanceSourceFormat::Nist,
                composition: "{C:1, O:1}".to_string(),
                fields: BTreeMap::from([("Cp".to_string(), json!([1.0]))]),
            })
            .expect_err("mismatched format should fail");

        let message = err.to_string();
        assert!(message.contains("NASA_gas"));
        assert!(message.contains("NIST"));

        let journal: Vec<ThermoLibraryJournalEntry> =
            serde_json::from_str(&fs::read_to_string(&journal_path).unwrap()).unwrap();
        assert!(journal.is_empty());
    }

    #[test]
    fn upsert_then_drop_restores_original_catalog_snapshot() {
        let dir = tempdir().expect("tempdir should be created");
        let config_path = dir.path().join("library_config.json");
        let substance_base_path = dir.path().join("substance_base.json");
        let all_keys_path = dir.path().join("all_keys_substance.json");
        let elements_path = dir.path().join("elements.json");
        let journal_path = dir.path().join("thermo_write_journal.json");

        let substance_base = json!({
            "NASA_gas": {
                "H2O": { "composition": "{H:2, O:1}", "Cp": [2.0] }
            }
        });
        let all_keys = json!([["NASA_gas", "H2O"]]);
        let elements = json!({
            "H": [["H2O", "NASA_gas"]],
            "O": [["H2O", "NASA_gas"]]
        });

        let original_base: serde_json::Value = substance_base.clone();
        let original_keys: serde_json::Value = all_keys.clone();
        let original_elements: serde_json::Value = elements.clone();

        fs::write(
            &substance_base_path,
            serde_json::to_string_pretty(&original_base).unwrap(),
        )
        .unwrap();
        fs::write(
            &all_keys_path,
            serde_json::to_string_pretty(&original_keys).unwrap(),
        )
        .unwrap();
        fs::write(
            &elements_path,
            serde_json::to_string_pretty(&original_elements).unwrap(),
        )
        .unwrap();
        fs::write(&journal_path, "[]").unwrap();

        let mut manager = LibraryManager::with_config_file(config_path.to_str().unwrap());
        manager
            .set_substance_base(substance_base_path.to_str().unwrap())
            .unwrap();
        manager
            .set_all_keys_substance(all_keys_path.to_str().unwrap())
            .unwrap();
        manager
            .set_elements(elements_path.to_str().unwrap())
            .unwrap();
        manager
            .set_thermo_write_journal(journal_path.to_str().unwrap())
            .unwrap();

        manager
            .write_substance_from_spec(ThermoSubstanceWriteSpec {
                library: "NASA_gas".to_string(),
                substance: "CH4".to_string(),
                source_format: ThermoSubstanceSourceFormat::Nasa,
                composition: "{C:1, H:4}".to_string(),
                fields: BTreeMap::from([("Cp".to_string(), json!([4.0]))]),
            })
            .unwrap();
        manager
            .drop_substance_from_library("NASA_gas", "CH4")
            .unwrap();

        let final_base: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(&substance_base_path).unwrap()).unwrap();
        let final_keys: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(&all_keys_path).unwrap()).unwrap();
        let final_elements: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(&elements_path).unwrap()).unwrap();

        assert_eq!(final_base, original_base);
        assert_eq!(final_keys, original_keys);
        assert_eq!(final_elements, original_elements);

        let journal: Vec<ThermoLibraryJournalEntry> =
            serde_json::from_str(&fs::read_to_string(&journal_path).unwrap()).unwrap();
        assert_eq!(journal.len(), 2);
        assert_eq!(journal[0].source_format, "NASA");
        assert_eq!(journal[0].operation, ThermoLibraryWriteOperation::Upsert);
        assert_eq!(journal[1].source_format, "manual");
        assert_eq!(journal[1].operation, ThermoLibraryWriteOperation::Drop);
    }
}
