//! # Kinetics Library API Module
//!
//! ## Purpose
//! This module provides a comprehensive API for accessing and searching chemical reaction databases.
//! It enables users to load kinetic libraries, search for reactions by various criteria, and retrieve
//! detailed reaction data including Arrhenius parameters and reaction mechanisms.
//!
//! ## Main Data Structures
//! - `KineticData`: Core structure that manages reaction libraries and provides search functionality
//!   - `LibKineticData`: HashMap storing reaction ID -> reaction data mappings
//!   - `HashMapOfReactantsAndProducts`: Nested HashMap for substance-based searches
//!   - `AllLibraries`: Vector of available reaction library names
//!   - Search result vectors for found reactions by products/reagents
//!
//! ## Key Logic Implementation
//! 1. **Library Loading**: Reads JSON files containing reaction databases (Reactbase.json, dict_reaction.json)
//! 2. **Multi-criteria Search**: Supports searching by reaction ID, equation, substances, or field values
//! 3. **Substance Matching**: Uses HashSet operations for efficient subset matching of reagents/products
//! 4. **Data Retrieval**: Returns structured reaction data as serde_json::Value objects
//!
//! ## Usage Pattern
//! ```rust
//! use KiThe::Kinetics::kinetics_lib_api::KineticData;
//! let mut kin_instance = KineticData::new();
//! kin_instance.open_json_files("NUIG");  // Load library
//! kin_instance.print_all_reactions();    // Get all reactions
//! kin_instance.search_reaction_by_reagents_and_products(vec!["CO".to_string()]);
//! println!("Found reactions: {:?}", kin_instance.FoundReactionsByProducts);
//! ```
//!
//! ## Interesting Features
//! - **Dual JSON Structure**: Uses separate files for kinetic parameters and substance mappings
//! - **Flexible Search**: Supports both exact matches and subset searches for reaction participants
//! - **Library Agnostic**: Can work with multiple reaction databases (NUIG, Cantera, Aramco, etc.)
//! - **Efficient Indexing**: Pre-builds equation-to-ID mappings for fast lookups
//! - **Reusable loaded view**: Open a library once, run queries, and call
//!   `close()` to release instance-owned data. The parsed process cache remains
//!   available for other callers and later opens.
//! - **Opt-in load diagnostics**: `open_json_files_with_diagnostics()` reports
//!   cache, file I/O, JSON parsing, and library materialization stages for a
//!   higher-level telemetry collector. The ordinary open path does not time them.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, Instant, SystemTime};

use crate::Kinetics::error::{KineticsError, KineticsResult};
use crate::library_manager::with_library_manager;
use serde::de::DeserializeOwned;
use serde_json::Value;
use sha2::{Digest, Sha256};

type LibraryReactbase = HashMap<String, HashMap<String, Value>>;
type LibraryReactionDb = HashMap<String, HashMap<String, HashMap<String, HashSet<String>>>>;

#[derive(Debug, Clone, PartialEq, Eq)]
struct FileRevision {
    path: PathBuf,
    length: u64,
    modified: Option<SystemTime>,
}

struct CachedLibrary<T> {
    revision: FileRevision,
    content_sha256: String,
    value: Arc<T>,
}

type ReactbaseCache = Mutex<Option<CachedLibrary<LibraryReactbase>>>;
type ReactionDbCache = Mutex<Option<CachedLibrary<LibraryReactionDb>>>;

static REACTBASE_CACHE: OnceLock<ReactbaseCache> = OnceLock::new();
static REACTION_DB_CACHE: OnceLock<ReactionDbCache> = OnceLock::new();

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Result of resolving a catalog against the process-wide parsed cache.
pub enum CatalogCacheStatus {
    Hit,
    Miss,
    RevisionChanged,
    ConcurrentLoadReused,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Catalog stage active when a load failed.
pub enum CatalogLoadStage {
    RevisionCheck,
    CacheLookup,
    FileOpen,
    FileRead,
    JsonDeserialize,
    LibrarySelection,
    LibraryMaterialization,
}

#[derive(Debug, Clone, Default)]
/// Per-file observations. Durations are present only on the explicit diagnostic path.
pub struct CatalogFileDiagnostics {
    pub path: PathBuf,
    pub bytes: Option<u64>,
    pub content_sha256: Option<String>,
    pub cache_status: Option<CatalogCacheStatus>,
    pub revision_check: Option<Duration>,
    pub revision_confirmation: Option<Duration>,
    pub cache_lookup: Option<Duration>,
    pub file_open: Option<Duration>,
    pub file_read: Option<Duration>,
    pub json_deserialize: Option<Duration>,
    pub failure_stage: Option<CatalogLoadStage>,
}

#[derive(Debug, Clone, Default)]
/// Optional load report intended to be aggregated by the calling application.
pub struct CatalogLoadDiagnostics {
    pub library: String,
    pub configuration_snapshot: Option<Duration>,
    pub library_selection: Option<Duration>,
    pub library_materialization: Option<Duration>,
    pub reactbase: CatalogFileDiagnostics,
    pub reaction_database: CatalogFileDiagnostics,
}

#[derive(Debug)]
/// A catalog-loading error with diagnostics for stages completed before failure.
pub struct CatalogLoadFailure {
    pub diagnostics: CatalogLoadDiagnostics,
    pub source: KineticsError,
}

impl std::fmt::Display for CatalogLoadFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "catalog load failed: {}", self.source)
    }
}

impl std::error::Error for CatalogLoadFailure {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

fn cache_lock_error(cache_name: &'static str) -> KineticsError {
    KineticsError::CacheLockPoisoned { cache: cache_name }
}

fn accumulate_duration(slot: &mut Option<Duration>, elapsed: Duration) {
    *slot = Some(slot.unwrap_or_default().saturating_add(elapsed));
}

#[cfg(test)]
fn read_json_library<T: DeserializeOwned>(path: &str, library: &'static str) -> KineticsResult<T> {
    read_json_library_with_hash(path, library).map(|(value, _)| value)
}

fn read_json_library_with_hash<T: DeserializeOwned>(
    path: &str,
    library: &'static str,
) -> KineticsResult<(T, String)> {
    let mut file = File::open(path).map_err(|error| {
        KineticsError::Io(error).with_context(format!("opening {library} at `{path}`"))
    })?;
    let mut contents = String::new();
    file.read_to_string(&mut contents).map_err(|error| {
        KineticsError::Io(error).with_context(format!("reading {library} at `{path}`"))
    })?;
    let digest = format!("{:x}", Sha256::digest(contents.as_bytes()));
    let parsed = serde_json::from_str(&contents).map_err(|error| {
        KineticsError::Json(error).with_context(format!("parsing {library} at `{path}`"))
    })?;
    Ok((parsed, digest))
}

fn read_json_library_with_diagnostics<T: DeserializeOwned>(
    path: &str,
    library: &'static str,
    diagnostics: &mut CatalogFileDiagnostics,
) -> KineticsResult<(T, String)> {
    let started = Instant::now();
    let mut file = match File::open(path) {
        Ok(file) => file,
        Err(error) => {
            diagnostics.file_open = Some(started.elapsed());
            diagnostics.failure_stage = Some(CatalogLoadStage::FileOpen);
            return Err(
                KineticsError::Io(error).with_context(format!("opening {library} at `{path}`"))
            );
        }
    };
    diagnostics.file_open = Some(started.elapsed());

    let started = Instant::now();
    let mut contents = String::new();
    if let Err(error) = file.read_to_string(&mut contents) {
        diagnostics.file_read = Some(started.elapsed());
        diagnostics.failure_stage = Some(CatalogLoadStage::FileRead);
        return Err(KineticsError::Io(error).with_context(format!("reading {library} at `{path}`")));
    }
    diagnostics.file_read = Some(started.elapsed());
    diagnostics.bytes = Some(contents.len() as u64);
    let digest = format!("{:x}", Sha256::digest(contents.as_bytes()));
    diagnostics.content_sha256 = Some(digest.clone());

    let started = Instant::now();
    match serde_json::from_str(&contents).map_err(|error| {
        KineticsError::Json(error).with_context(format!("parsing {library} at `{path}`"))
    }) {
        Ok(parsed) => {
            diagnostics.json_deserialize = Some(started.elapsed());
            Ok((parsed, digest))
        }
        Err(error) => {
            diagnostics.json_deserialize = Some(started.elapsed());
            diagnostics.failure_stage = Some(CatalogLoadStage::JsonDeserialize);
            Err(error)
        }
    }
}

fn file_revision(path: &str, library: &'static str) -> KineticsResult<FileRevision> {
    let metadata = std::fs::metadata(path).map_err(|error| {
        KineticsError::Io(error).with_context(format!("checking {library} at `{path}`"))
    })?;
    Ok(FileRevision {
        path: PathBuf::from(path),
        length: metadata.len(),
        modified: metadata.modified().ok(),
    })
}

fn cached_json_library<T: DeserializeOwned + Send + Sync + 'static>(
    cache: &Mutex<Option<CachedLibrary<T>>>,
    cache_name: &'static str,
    path: &str,
    library: &'static str,
) -> KineticsResult<Arc<T>> {
    cached_json_library_inner(cache, cache_name, path, library, None)
}

fn cached_json_library_observed<T: DeserializeOwned + Send + Sync + 'static>(
    cache: &Mutex<Option<CachedLibrary<T>>>,
    cache_name: &'static str,
    path: &str,
    library: &'static str,
    diagnostics: &mut CatalogFileDiagnostics,
) -> KineticsResult<Arc<T>> {
    diagnostics.path = PathBuf::from(path);
    cached_json_library_inner(cache, cache_name, path, library, Some(diagnostics))
}

fn cached_json_library_inner<T: DeserializeOwned + Send + Sync + 'static>(
    cache: &Mutex<Option<CachedLibrary<T>>>,
    cache_name: &'static str,
    path: &str,
    library: &'static str,
    mut diagnostics: Option<&mut CatalogFileDiagnostics>,
) -> KineticsResult<Arc<T>> {
    for _ in 0..2 {
        let revision_started = diagnostics.as_ref().map(|_| Instant::now());
        let revision_result = file_revision(path, library);
        if let (Some(diagnostics), Some(started)) = (diagnostics.as_deref_mut(), revision_started) {
            accumulate_duration(&mut diagnostics.revision_check, started.elapsed());
        }
        let revision = match revision_result {
            Ok(revision) => revision,
            Err(error) => {
                if let Some(diagnostics) = diagnostics.as_deref_mut() {
                    diagnostics.failure_stage = Some(CatalogLoadStage::RevisionCheck);
                }
                return Err(error);
            }
        };
        if let Some(diagnostics) = diagnostics.as_deref_mut() {
            diagnostics.bytes = Some(revision.length);
        }
        {
            let lookup_started = diagnostics.as_ref().map(|_| Instant::now());
            let guard_result = cache.lock();
            let guard = match guard_result {
                Ok(guard) => guard,
                Err(_) => {
                    if let Some(diagnostics) = diagnostics.as_deref_mut() {
                        if let Some(started) = lookup_started {
                            accumulate_duration(&mut diagnostics.cache_lookup, started.elapsed());
                        }
                        diagnostics.failure_stage = Some(CatalogLoadStage::CacheLookup);
                    }
                    return Err(cache_lock_error(cache_name));
                }
            };
            let cached = guard
                .as_ref()
                .filter(|entry| entry.revision == revision)
                .map(|entry| (entry.value.clone(), entry.content_sha256.clone()));
            if let Some(diagnostics) = diagnostics.as_deref_mut() {
                if let Some(started) = lookup_started {
                    accumulate_duration(&mut diagnostics.cache_lookup, started.elapsed());
                }
                if cached.is_some() {
                    diagnostics.cache_status = Some(CatalogCacheStatus::Hit);
                    diagnostics.content_sha256 = cached.as_ref().map(|(_, digest)| digest.clone());
                } else {
                    diagnostics.cache_status = Some(if guard.is_some() {
                        CatalogCacheStatus::RevisionChanged
                    } else {
                        CatalogCacheStatus::Miss
                    });
                }
            }
            if let Some(cached) = cached {
                return Ok(cached.0);
            }
        }

        let (parsed_value, content_sha256) = match diagnostics.as_deref_mut() {
            Some(diagnostics) => read_json_library_with_diagnostics(path, library, diagnostics)?,
            None => read_json_library_with_hash(path, library)?,
        };
        let parsed: Arc<T> = Arc::new(parsed_value);
        let confirmation_started = diagnostics.as_ref().map(|_| Instant::now());
        let confirmation_result = file_revision(path, library);
        if let (Some(diagnostics), Some(started)) =
            (diagnostics.as_deref_mut(), confirmation_started)
        {
            accumulate_duration(&mut diagnostics.revision_confirmation, started.elapsed());
        }
        let confirmed_revision = match confirmation_result {
            Ok(revision) => revision,
            Err(error) => {
                if let Some(diagnostics) = diagnostics.as_deref_mut() {
                    diagnostics.failure_stage = Some(CatalogLoadStage::RevisionCheck);
                }
                return Err(error);
            }
        };
        if confirmed_revision != revision {
            continue;
        }

        let lookup_started = diagnostics.as_ref().map(|_| Instant::now());
        let guard_result = cache.lock();
        let mut guard = match guard_result {
            Ok(guard) => guard,
            Err(_) => {
                if let Some(diagnostics) = diagnostics.as_deref_mut() {
                    if let Some(started) = lookup_started {
                        accumulate_duration(&mut diagnostics.cache_lookup, started.elapsed());
                    }
                    diagnostics.failure_stage = Some(CatalogLoadStage::CacheLookup);
                }
                return Err(cache_lock_error(cache_name));
            }
        };
        let cached = guard
            .as_ref()
            .filter(|entry| entry.revision == revision)
            .map(|entry| (entry.value.clone(), entry.content_sha256.clone()));
        if let Some(diagnostics) = diagnostics.as_deref_mut() {
            if let Some(started) = lookup_started {
                accumulate_duration(&mut diagnostics.cache_lookup, started.elapsed());
            }
            if cached.is_some() {
                diagnostics.cache_status = Some(CatalogCacheStatus::ConcurrentLoadReused);
                diagnostics.content_sha256 = cached.as_ref().map(|(_, digest)| digest.clone());
            }
        }
        if let Some(cached) = cached {
            return Ok(cached.0);
        }
        *guard = Some(CachedLibrary {
            revision,
            content_sha256,
            value: parsed.clone(),
        });
        return Ok(parsed);
    }

    if let Some(diagnostics) = diagnostics.as_deref_mut() {
        diagnostics.failure_stage = Some(CatalogLoadStage::RevisionCheck);
    }
    Err(KineticsError::CatalogChangedDuringLoad {
        path: path.to_string(),
    })
}

fn configured_catalog_paths() -> (String, String) {
    with_library_manager(|manager| {
        (
            manager.reactbase_path().to_owned(),
            manager.dict_reaction_path().to_owned(),
        )
    })
}

fn cached_reactbase_at(path: &str) -> KineticsResult<Arc<LibraryReactbase>> {
    cached_json_library(
        REACTBASE_CACHE.get_or_init(|| Mutex::new(None)),
        "Reactbase",
        path,
        "reaction parameter library",
    )
}

fn cached_reaction_db_at(path: &str) -> KineticsResult<Arc<LibraryReactionDb>> {
    cached_json_library(
        REACTION_DB_CACHE.get_or_init(|| Mutex::new(None)),
        "reaction database",
        path,
        "reaction relationship database",
    )
}

/// Load and cache the reaction parameter library once per process.
#[cfg(test)]
pub(crate) fn cached_reactbase() -> KineticsResult<Arc<LibraryReactbase>> {
    let (path, _) = configured_catalog_paths();
    cached_reactbase_at(&path)
}

/// Load and cache the reagent/product relationship database once per process.
#[cfg(test)]
pub(crate) fn cached_reaction_db() -> KineticsResult<Arc<LibraryReactionDb>> {
    let (_, path) = configured_catalog_paths();
    cached_reaction_db_at(&path)
}

pub(crate) fn cached_catalog_pair()
-> KineticsResult<(Arc<LibraryReactbase>, Arc<LibraryReactionDb>)> {
    let (reactbase_path, reaction_db_path) = configured_catalog_paths();
    Ok((
        cached_reactbase_at(&reactbase_path)?,
        cached_reaction_db_at(&reaction_db_path)?,
    ))
}

fn cached_catalog_pair_observed(
    diagnostics: &mut CatalogLoadDiagnostics,
) -> KineticsResult<(Arc<LibraryReactbase>, Arc<LibraryReactionDb>)> {
    let started = Instant::now();
    let (reactbase_path, reaction_db_path) = configured_catalog_paths();
    diagnostics.configuration_snapshot = Some(started.elapsed());
    diagnostics.reactbase.path = PathBuf::from(&reactbase_path);
    diagnostics.reaction_database.path = PathBuf::from(&reaction_db_path);
    let reactbase = cached_json_library_observed(
        REACTBASE_CACHE.get_or_init(|| Mutex::new(None)),
        "Reactbase",
        &reactbase_path,
        "reaction parameter library",
        &mut diagnostics.reactbase,
    )?;
    let reaction_db = cached_json_library_observed(
        REACTION_DB_CACHE.get_or_init(|| Mutex::new(None)),
        "reaction database",
        &reaction_db_path,
        "reaction relationship database",
        &mut diagnostics.reaction_database,
    )?;
    Ok((reactbase, reaction_db))
}
#[derive(Debug, Clone)]
pub struct KineticData {
    pub HashMapOfReactantsAndProducts: HashMap<String, HashMap<String, HashSet<String>>>, // {'reaction ID':{'reagent'/"product": HashSet[substance]}}
    pub lib_name: String,
    pub LibKineticData: HashMap<String, Value>, // {'reaction ID':{data of reaction}}
    pub AllLibraries: Vec<String>,              // all reaction libraries
    pub AllEquations: Vec<String>,              // all reaction equations in library
    pub EquationReactionMap: HashMap<String, Vec<String>>, // equation -> all matching reaction IDs
    pub UserEquationReactionMap: HashMap<String, String>,
    pub FoundReactionsByProducts: Vec<String>,
    pub FoundReactionsByReagents: Vec<String>,
    pub FoundReactionDatasByIDs: Vec<Value>,
}
impl KineticData {
    pub fn new() -> Self {
        Self {
            HashMapOfReactantsAndProducts: HashMap::new(),
            lib_name: String::new(),
            LibKineticData: HashMap::new(),
            AllLibraries: Vec::new(),
            AllEquations: Vec::new(),
            EquationReactionMap: HashMap::new(),
            UserEquationReactionMap: HashMap::new(),
            FoundReactionsByProducts: Vec::new(),
            FoundReactionsByReagents: Vec::new(),
            FoundReactionDatasByIDs: Vec::new(),
        }
    }
    /// Load the reaction library from a JSON file by the name of the reaction library
    pub fn open_json_files(&mut self, big_mech: &str) -> KineticsResult<()> {
        self.open_json_files_inner(big_mech, None)
    }

    /// Load a library and return opt-in diagnostics for the caller's telemetry.
    ///
    /// Unlike `open_json_files`, this method measures catalog loading stages.
    /// On failure, partial diagnostics are returned alongside the typed cause.
    pub fn open_json_files_with_diagnostics(
        &mut self,
        big_mech: &str,
    ) -> Result<CatalogLoadDiagnostics, CatalogLoadFailure> {
        let mut diagnostics = CatalogLoadDiagnostics {
            library: big_mech.to_string(),
            ..CatalogLoadDiagnostics::default()
        };
        match self.open_json_files_inner(big_mech, Some(&mut diagnostics)) {
            Ok(()) => Ok(diagnostics),
            Err(source) => Err(CatalogLoadFailure {
                diagnostics,
                source,
            }),
        }
    }

    fn open_json_files_inner(
        &mut self,
        big_mech: &str,
        mut diagnostics: Option<&mut CatalogLoadDiagnostics>,
    ) -> KineticsResult<()> {
        let (reactlibrary, reaction_db) = match diagnostics.as_deref_mut() {
            Some(diagnostics) => cached_catalog_pair_observed(diagnostics)?,
            None => cached_catalog_pair()?,
        };
        let all_libraries = {
            let mut libraries: Vec<_> = reactlibrary.keys().cloned().collect();
            libraries.sort();
            libraries
        };
        let selection_started = diagnostics.as_ref().map(|_| Instant::now());
        let library_of_kinetic_parameters = match reactlibrary.get(big_mech) {
            Some(library) => library,
            None => {
                if let Some(diagnostics) = diagnostics.as_deref_mut() {
                    diagnostics.library_selection = selection_started.map(|s| s.elapsed());
                    diagnostics.reactbase.failure_stage = Some(CatalogLoadStage::LibrarySelection);
                }
                return Err(KineticsError::MissingLibrary(big_mech.to_string()));
            }
        };
        let library_of_reagents_and_products: &HashMap<String, HashMap<String, HashSet<String>>> =
            match reaction_db.get(big_mech) {
                Some(library) => library,
                None => {
                    if let Some(diagnostics) = diagnostics.as_deref_mut() {
                        diagnostics.library_selection = selection_started.map(|s| s.elapsed());
                        diagnostics.reaction_database.failure_stage =
                            Some(CatalogLoadStage::LibrarySelection);
                    }
                    return Err(KineticsError::MissingLibrary(big_mech.to_string()));
                }
            };
        if let Some(diagnostics) = diagnostics.as_deref_mut() {
            diagnostics.library_selection = selection_started.map(|s| s.elapsed());
        }

        let materialization_started = diagnostics.as_ref().map(|_| Instant::now());
        let kinetic_data = library_of_kinetic_parameters.to_owned();
        let reaction_parts = library_of_reagents_and_products.to_owned();
        if let Some(diagnostics) = diagnostics.as_deref_mut() {
            diagnostics.library_materialization = materialization_started.map(|s| s.elapsed());
        }
        self.lib_name = big_mech.to_owned();
        self.AllLibraries = all_libraries;
        self.LibKineticData = kinetic_data;
        self.HashMapOfReactantsAndProducts = reaction_parts;
        self.AllEquations.clear();
        self.EquationReactionMap.clear();
        self.UserEquationReactionMap.clear();
        self.FoundReactionsByProducts.clear();
        self.FoundReactionsByReagents.clear();
        self.FoundReactionDatasByIDs.clear();
        Ok(())
    }

    /// Whether this instance currently owns an opened library view.
    pub fn is_open(&self) -> bool {
        !self.lib_name.is_empty()
    }

    /// Release this instance's selected-library view and derived query state.
    ///
    /// The shared process-wide parsed catalog cache is intentionally retained;
    /// other users may still hold references to it, and subsequent opens can
    /// reuse it without reparsing the JSON files.
    pub fn close(&mut self) {
        *self = Self::new();
    }

    /// returns vector of all reaction equations and HashMap {reaction equation : reaction ID}
    pub fn print_all_reactions(&mut self) -> KineticsResult<()> {
        let mut all_reactions: Vec<String> = self.LibKineticData.keys().cloned().collect();
        all_reactions.sort();
        let all_equations: Vec<String> = all_reactions
            .iter()
            .map(|k| {
                self.LibKineticData
                    .get(k)
                    .and_then(|reaction| reaction.get("eq"))
                    .and_then(|eq| eq.as_str())
                    .map(|eq| eq.to_string())
                    .ok_or_else(|| KineticsError::MissingReaction(k.clone()))
            })
            .collect::<KineticsResult<Vec<String>>>()?;
        self.AllEquations = all_equations.clone();
        let mut equation_reaction_ids: HashMap<String, Vec<String>> = HashMap::new();
        for (equation, reaction_id) in all_equations.iter().zip(all_reactions.iter()) {
            equation_reaction_ids
                .entry(equation.clone())
                .or_default()
                .push(reaction_id.clone());
        }
        for reaction_ids in equation_reaction_ids.values_mut() {
            reaction_ids.sort();
        }
        self.EquationReactionMap = equation_reaction_ids;
        Ok(())
    } //end print_all_reactions
    /// Return the unique reaction for an equation, rejecting duplicate matches.
    pub fn search_reaction_by_equation(
        &mut self,
        equation: &str,
    ) -> KineticsResult<(String, Value)> {
        let reaction_ids = self
            .EquationReactionMap
            .get(equation)
            .ok_or_else(|| KineticsError::MissingReaction(equation.to_string()))?;
        if reaction_ids.len() != 1 {
            return Err(KineticsError::AmbiguousReactionEquation {
                equation: equation.to_string(),
                reaction_ids: reaction_ids.clone(),
            });
        }
        let reaction_id = &reaction_ids[0];
        let reaction = self
            .LibKineticData
            .get(reaction_id)
            .ok_or_else(|| KineticsError::MissingReaction(reaction_id.clone()))?;
        Ok((reaction_id.clone(), reaction.clone()))
    }

    /// Return every reaction whose catalog equation exactly matches `equation`.
    pub fn search_reactions_by_equation(
        &self,
        equation: &str,
    ) -> KineticsResult<Vec<(String, Value)>> {
        let reaction_ids = self
            .EquationReactionMap
            .get(equation)
            .ok_or_else(|| KineticsError::MissingReaction(equation.to_string()))?;
        reaction_ids
            .iter()
            .map(|reaction_id| {
                let reaction = self
                    .LibKineticData
                    .get(reaction_id)
                    .ok_or_else(|| KineticsError::MissingReaction(reaction_id.clone()))?;
                Ok((reaction_id.clone(), reaction.clone()))
            })
            .collect()
    }

    // /returns reaction ID and reaction data (parsed from json) for given reagents/products
    fn search_reaction_by_substances(
        &mut self,
        substances: Vec<String>,
        SubstanceType: &str,
    ) -> KineticsResult<Vec<String>> {
        let substances: HashSet<String> = substances.iter().cloned().collect();
        let mut found_reactions: Vec<String> = Vec::new();
        for (r_id, reaction_parts) in &self.HashMapOfReactantsAndProducts {
            let Some(base_substances) = reaction_parts.get(SubstanceType) else {
                return Err(KineticsError::MissingReaction(r_id.to_string()));
            };
            if base_substances
                .iter()
                .filter(|subs| !matches!(subs.as_str(), "M" | "M)" | "(M" | "(M)"))
                .all(|subs| substances.contains(subs))
            {
                found_reactions.push(r_id.to_string());
            }
        }
        found_reactions.sort();
        Ok(found_reactions)
    }
    /// search reactions where these substances are reagents/products
    pub fn search_reaction_by_reagents_and_products(
        &mut self,
        reagents: Vec<String>,
    ) -> KineticsResult<()> {
        let found_reactions_by_reagents =
            self.search_reaction_by_substances(reagents.clone(), "reagents")?;
        self.FoundReactionsByReagents = found_reactions_by_reagents;
        let found_reactions_by_products =
            self.search_reaction_by_substances(reagents, "products")?;
        self.FoundReactionsByProducts = found_reactions_by_products;
        Ok(())
    }
    /// get concrete reaction data
    pub fn search_reactdata_by_reaction_id(&mut self, reaction_id: &str) -> KineticsResult<Value> {
        let reaction = self
            .LibKineticData
            .get(reaction_id)
            .ok_or_else(|| KineticsError::MissingReaction(reaction_id.to_string()))?;
        return Ok(reaction.clone());
    }
    /// get concrete reaction data for vector of IDs
    pub fn search_reactdata_for_vector_of_IDs(
        &mut self,
        reaction_ids: Vec<String>,
    ) -> KineticsResult<Vec<Value>> {
        let mut vec_of_reactions: Vec<Value> = Vec::new();
        for r_id in reaction_ids {
            let reaction = self
                .LibKineticData
                .get(&r_id)
                .ok_or_else(|| KineticsError::MissingReaction(r_id.clone()))?;
            vec_of_reactions.push(reaction.clone());
        }
        self.FoundReactionDatasByIDs = vec_of_reactions.clone();
        Ok(vec_of_reactions)
    }
    //TODO! test
    /// get reaction data for a certain value of json field
    pub fn get_reactions_by_field(&self, field: &str, field_value: &str) -> Vec<Value> {
        self.LibKineticData
            .values()
            .filter(|reaction| {
                reaction
                    .get(field)
                    .and_then(|v| v.as_str())
                    .map_or(false, |v| v == field_value)
            })
            .cloned()
            .collect()
    }
    #[allow(dead_code)]
    /// form user chosen data
    pub fn form_user_chosen_data(
        &self,
        what_to_get: SearchVariants,
        mech_name: Option<String>,
    ) -> KineticsResult<HashMap<String, HashMap<String, Value>>> {
        let big_mech = if let Some(mech_name) = { mech_name } {
            mech_name
        } else {
            self.lib_name.clone()
        };
        let reaction_ids = self.get_search_results(what_to_get);
        if reaction_ids.len() == 0 {
            return Ok(HashMap::new());
        }
        let mut hash_map_id_data = HashMap::with_capacity(reaction_ids.len());
        for r_id in reaction_ids {
            let reaction = self
                .LibKineticData
                .get(&r_id)
                .ok_or_else(|| KineticsError::MissingReaction(r_id.clone()))?;
            hash_map_id_data.insert(r_id, reaction.clone());
        }
        let mut hash_map = HashMap::with_capacity(1);
        hash_map.insert(big_mech, hash_map_id_data);
        Ok(hash_map)
    }
    ///
    #[allow(dead_code)]
    pub fn get_search_results(&self, what_to_get: SearchVariants) -> Vec<String> {
        match what_to_get {
            SearchVariants::FoundReactionsByProducts => self.FoundReactionsByProducts.clone(),
            SearchVariants::FoundReactionsByReagents => self.FoundReactionsByReagents.clone(),
        }
    }
}

pub enum SearchVariants {
    FoundReactionsByProducts,
    FoundReactionsByReagents,
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TESTS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn poisoned_cache_lock_maps_to_a_typed_error() {
        assert!(matches!(
            cache_lock_error("Reactbase"),
            KineticsError::CacheLockPoisoned { cache: "Reactbase" }
        ));
    }

    #[test]
    fn library_json_errors_keep_stage_path_and_typed_source() {
        let temp = tempfile::tempdir().expect("temporary fixture directory");
        let missing_path = temp.path().join("missing.json");
        let missing_path = missing_path.to_string_lossy().into_owned();
        let open_error = read_json_library::<Value>(&missing_path, "test library").unwrap_err();
        let KineticsError::Context { context, source } = &open_error else {
            panic!("missing open-stage context: {open_error}");
        };
        assert!(context.contains("opening test library"));
        assert!(context.contains(&missing_path));
        assert!(matches!(**source, KineticsError::Io(_)));

        let invalid_path = temp.path().join("invalid.json");
        std::fs::write(&invalid_path, "{").expect("write malformed JSON fixture");
        let invalid_path = invalid_path.to_string_lossy().into_owned();
        let parse_error = read_json_library::<Value>(&invalid_path, "test library").unwrap_err();
        let KineticsError::Context { context, source } = &parse_error else {
            panic!("missing parse-stage context: {parse_error}");
        };
        assert!(context.contains("parsing test library"));
        assert!(context.contains(&invalid_path));
        assert!(matches!(**source, KineticsError::Json(_)));
    }

    #[test]
    fn catalog_cache_tracks_path_and_file_revision_and_retries_after_failure() {
        let temp = tempfile::tempdir().expect("temporary fixture directory");
        let path_a = temp.path().join("a.json");
        let path_b = temp.path().join("b.json");
        std::fs::write(&path_a, r#"{"value":1}"#).unwrap();
        std::fs::write(&path_b, r#"{"value":333}"#).unwrap();
        let path_a = path_a.to_string_lossy().into_owned();
        let path_b = path_b.to_string_lossy().into_owned();
        let cache = Mutex::new(None);

        let first = cached_json_library::<Value>(&cache, "test", &path_a, "test fixture").unwrap();
        let reused = cached_json_library::<Value>(&cache, "test", &path_a, "test fixture").unwrap();
        assert!(Arc::ptr_eq(&first, &reused));

        std::fs::write(&path_a, r#"{"value":22}"#).unwrap();
        let revised =
            cached_json_library::<Value>(&cache, "test", &path_a, "test fixture").unwrap();
        assert_eq!(revised["value"], 22);
        assert!(!Arc::ptr_eq(&first, &revised));

        let switched =
            cached_json_library::<Value>(&cache, "test", &path_b, "test fixture").unwrap();
        assert_eq!(switched["value"], 333);
        assert!(!Arc::ptr_eq(&revised, &switched));

        let retry_path = temp
            .path()
            .join("retry.json")
            .to_string_lossy()
            .into_owned();
        assert!(cached_json_library::<Value>(&cache, "test", &retry_path, "test fixture").is_err());
        std::fs::write(&retry_path, r#"{"ready":true}"#).unwrap();
        let retry = cached_json_library::<Value>(&cache, "test", &retry_path, "test fixture")
            .expect("failed first load must not poison the cache entry");
        assert_eq!(retry["ready"], true);
    }

    #[test]
    fn catalog_diagnostics_distinguish_cold_parse_from_cache_hit() {
        let temp = tempfile::tempdir().expect("temporary fixture directory");
        let path = temp.path().join("observed.json");
        std::fs::write(&path, r#"{"value":7}"#).unwrap();
        let path = path.to_string_lossy().into_owned();
        let cache = Mutex::new(None);
        let mut cold = CatalogFileDiagnostics::default();

        cached_json_library_observed::<Value>(&cache, "test", &path, "test fixture", &mut cold)
            .unwrap();
        assert_eq!(cold.cache_status, Some(CatalogCacheStatus::Miss));
        assert_eq!(cold.bytes, Some(r#"{"value":7}"#.len() as u64));
        assert!(cold.file_open.is_some());
        assert!(cold.file_read.is_some());
        assert!(cold.json_deserialize.is_some());
        assert!(cold.content_sha256.is_some());
        assert_eq!(cold.failure_stage, None);

        let mut warm = CatalogFileDiagnostics::default();
        cached_json_library_observed::<Value>(&cache, "test", &path, "test fixture", &mut warm)
            .unwrap();
        assert_eq!(warm.cache_status, Some(CatalogCacheStatus::Hit));
        assert_eq!(warm.content_sha256, cold.content_sha256);
        assert!(warm.file_open.is_none());
        assert!(warm.file_read.is_none());
        assert!(warm.json_deserialize.is_none());

        std::fs::write(&path, r#"{"value":77}"#).unwrap();
        let mut revised = CatalogFileDiagnostics::default();
        cached_json_library_observed::<Value>(&cache, "test", &path, "test fixture", &mut revised)
            .unwrap();
        assert_eq!(
            revised.cache_status,
            Some(CatalogCacheStatus::RevisionChanged)
        );
        assert_ne!(revised.content_sha256, cold.content_sha256);
    }

    #[test]
    fn catalog_diagnostics_identify_json_parse_failure_stage() {
        let temp = tempfile::tempdir().expect("temporary fixture directory");
        let path = temp.path().join("invalid.json");
        std::fs::write(&path, "{").unwrap();
        let path = path.to_string_lossy().into_owned();
        let mut diagnostics = CatalogFileDiagnostics::default();

        let err =
            read_json_library_with_diagnostics::<Value>(&path, "test fixture", &mut diagnostics)
                .unwrap_err();
        assert!(matches!(err, KineticsError::Context { .. }));
        assert_eq!(
            diagnostics.failure_stage,
            Some(CatalogLoadStage::JsonDeserialize)
        );
        assert!(diagnostics.json_deserialize.is_some());
    }

    #[test]
    fn concurrent_catalog_initialization_returns_one_shared_cached_value() {
        let temp = tempfile::tempdir().expect("temporary fixture directory");
        let path = temp.path().join("concurrent.json");
        std::fs::write(&path, r#"{"value":7}"#).unwrap();
        let path = Arc::new(path.to_string_lossy().into_owned());
        let cache = Arc::new(Mutex::new(None));
        let workers: Vec<_> = (0..8)
            .map(|_| {
                let path = Arc::clone(&path);
                let cache = Arc::clone(&cache);
                std::thread::spawn(move || {
                    cached_json_library::<Value>(&cache, "test", &path, "test fixture").unwrap()
                })
            })
            .collect();
        let loaded: Vec<_> = workers
            .into_iter()
            .map(|worker| worker.join().unwrap())
            .collect();
        assert!(loaded.iter().all(|entry| Arc::ptr_eq(&loaded[0], entry)));
    }

    #[test]
    fn test_kinetic_data_new() {
        let mut kin_instance = KineticData::new();
        let lib = "NUIG";
        let LibVec: HashSet<String> = vec![
            "Beckstead".to_string(),
            "NUIG".to_string(),
            "Cantera".to_string(),
            "Beckstead_c".to_string(),
            "Aramco".to_string(),
        ]
        .into_iter()
        .collect();

        // collecting reaction data for library name lib
        let _ = kin_instance.open_json_files(lib);
        // veiew all reactions in library
        let _ = kin_instance.print_all_reactions();
        assert_eq!(kin_instance.HashMapOfReactantsAndProducts.is_empty(), false);
        assert_eq!(kin_instance.LibKineticData.is_empty(), false);
        let all_lib: HashSet<String> = kin_instance.AllLibraries.clone().into_iter().collect();
        assert_eq!(all_lib, LibVec);
        assert_eq!(kin_instance.AllEquations.is_empty(), false);

        // search reaction by equation
        let equation = kin_instance
            .LibKineticData
            .get("1")
            .unwrap()
            .get("eq")
            .unwrap()
            .as_str()
            .unwrap()
            .to_string();
        // returns reaction ID and reaction data (parsed from json) for given reaction equation
        let (reaction_id, reaction) = kin_instance.search_reaction_by_equation(&equation).unwrap();
        assert!(kin_instance.LibKineticData.contains_key(&reaction_id));
        assert_eq!(reaction.get("eq").unwrap().as_str().unwrap(), equation);
        let substances: Vec<String> = vec!["CO".to_string(), "O".to_string()];
        // search reactions by substances
        kin_instance
            .search_reaction_by_reagents_and_products(substances)
            .unwrap();
        assert_eq!(kin_instance.FoundReactionsByProducts.is_empty(), false);
        assert_eq!(kin_instance.FoundReactionsByReagents.is_empty(), false);
        kin_instance
            .search_reactdata_for_vector_of_IDs(vec!["1".to_string(), "2".to_string()])
            .unwrap();
        assert_eq!(kin_instance.FoundReactionDatasByIDs.is_empty(), false);
    }

    #[test]
    fn test_open_json_files_rejects_unknown_library() {
        let mut kin_instance = KineticData::new();
        kin_instance.lib_name = "existing".into();
        kin_instance
            .LibKineticData
            .insert("old".into(), serde_json::json!({"eq":"A=B"}));
        kin_instance.HashMapOfReactantsAndProducts.insert(
            "old".into(),
            HashMap::from([("reagents".into(), HashSet::from(["A".into()]))]),
        );
        kin_instance.AllEquations = vec!["A=B".into()];
        kin_instance
            .EquationReactionMap
            .insert("A=B".into(), vec!["old".into()]);
        kin_instance.FoundReactionsByProducts = vec!["old".into()];
        let err = kin_instance.open_json_files("NotARealLibrary").unwrap_err();

        assert!(matches!(
            err,
            KineticsError::MissingLibrary(library) if library == "NotARealLibrary"
        ));
        assert_eq!(kin_instance.lib_name, "existing");
        assert!(kin_instance.LibKineticData.contains_key("old"));
        assert!(
            kin_instance
                .HashMapOfReactantsAndProducts
                .contains_key("old")
        );
        assert_eq!(kin_instance.AllEquations, vec!["A=B"]);
        assert_eq!(
            kin_instance.EquationReactionMap.get("A=B"),
            Some(&vec!["old".into()])
        );
        assert_eq!(kin_instance.FoundReactionsByProducts, vec!["old"]);
    }

    #[test]
    fn duplicate_equations_retain_all_ids_and_singular_lookup_is_typed_ambiguous() {
        let mut data = KineticData::new();
        data.LibKineticData.insert(
            "R2".into(),
            serde_json::json!({"eq": "A => B", "type": "elem"}),
        );
        data.LibKineticData.insert(
            "R1".into(),
            serde_json::json!({"eq": "A => B", "type": "falloff"}),
        );
        data.print_all_reactions().unwrap();

        assert_eq!(
            data.EquationReactionMap.get("A => B"),
            Some(&vec!["R1".to_string(), "R2".to_string()])
        );
        let err = data.search_reaction_by_equation("A => B").unwrap_err();
        assert!(matches!(
            err,
            KineticsError::AmbiguousReactionEquation { equation, reaction_ids }
                if equation == "A => B" && reaction_ids == vec!["R1", "R2"]
        ));
        let all = data.search_reactions_by_equation("A => B").unwrap();
        assert_eq!(
            all.iter().map(|(id, _)| id.as_str()).collect::<Vec<_>>(),
            vec!["R1", "R2"]
        );
    }

    #[test]
    fn close_releases_instance_view_and_keeps_shared_catalog_cache_reusable() {
        let mut kinetic_data = KineticData::new();
        assert!(!kinetic_data.is_open());

        let first = kinetic_data
            .open_json_files_with_diagnostics("NUIG")
            .expect("catalog should open");
        assert_eq!(first.library, "NUIG");
        assert!(first.configuration_snapshot.is_some());
        assert!(first.library_selection.is_some());
        assert!(first.library_materialization.is_some());
        assert!(first.reactbase.path.exists());
        assert!(first.reaction_database.path.exists());
        assert!(first.reactbase.cache_status.is_some());
        assert!(first.reaction_database.cache_status.is_some());
        assert!(kinetic_data.is_open());
        kinetic_data.print_all_reactions().unwrap();
        assert!(!kinetic_data.AllEquations.is_empty());

        kinetic_data.close();
        assert!(!kinetic_data.is_open());
        assert!(kinetic_data.LibKineticData.is_empty());
        assert!(kinetic_data.HashMapOfReactantsAndProducts.is_empty());
        assert!(kinetic_data.AllLibraries.is_empty());
        assert!(kinetic_data.AllEquations.is_empty());
        assert!(kinetic_data.EquationReactionMap.is_empty());
        assert!(kinetic_data.UserEquationReactionMap.is_empty());
        assert!(kinetic_data.FoundReactionsByProducts.is_empty());
        assert!(kinetic_data.FoundReactionsByReagents.is_empty());
        assert!(kinetic_data.FoundReactionDatasByIDs.is_empty());

        let second = kinetic_data
            .open_json_files_with_diagnostics("NUIG")
            .expect("cached catalog should reopen");
        assert_eq!(second.reactbase.cache_status, Some(CatalogCacheStatus::Hit));
        assert_eq!(
            second.reaction_database.cache_status,
            Some(CatalogCacheStatus::Hit)
        );
    }

    #[test]
    fn catalog_diagnostics_include_partial_report_for_selection_failure() {
        let mut kinetic_data = KineticData::new();
        let failure = kinetic_data
            .open_json_files_with_diagnostics("NotARealLibrary")
            .unwrap_err();

        assert!(matches!(failure.source, KineticsError::MissingLibrary(_)));
        assert_eq!(failure.diagnostics.library, "NotARealLibrary");
        assert_eq!(
            failure.diagnostics.reactbase.failure_stage,
            Some(CatalogLoadStage::LibrarySelection)
        );
        assert!(failure.diagnostics.reactbase.cache_status.is_some());
        assert!(!kinetic_data.is_open());
    }

    #[test]
    fn test_search_reactdata_by_reaction_id_rejects_missing_reaction() {
        let mut kin_instance = KineticData::new();
        kin_instance.open_json_files("NUIG").unwrap();

        let err = kin_instance
            .search_reactdata_by_reaction_id("definitely_missing_id")
            .unwrap_err();

        assert!(matches!(
            err,
            KineticsError::MissingReaction(reaction_id)
                if reaction_id == "definitely_missing_id"
        ));
    }

    #[test]
    fn test_search_reaction_by_equation_rejects_missing_equation() {
        let mut kin_instance = KineticData::new();
        kin_instance.open_json_files("NUIG").unwrap();
        kin_instance.print_all_reactions().unwrap();

        let err = kin_instance
            .search_reaction_by_equation("A + B => C")
            .unwrap_err();

        assert!(matches!(
            err,
            KineticsError::MissingReaction(equation) if equation == "A + B => C"
        ));
    }

    #[test]
    fn test_search_reactdata_for_vector_of_ids_rejects_missing_reaction() {
        let mut kin_instance = KineticData::new();
        kin_instance.open_json_files("NUIG").unwrap();

        let err = kin_instance
            .search_reactdata_for_vector_of_IDs(vec!["definitely_missing_id".to_string()])
            .unwrap_err();

        assert!(matches!(
            err,
            KineticsError::MissingReaction(reaction_id)
                if reaction_id == "definitely_missing_id"
        ));
    }

    #[test]
    fn test_search_value() {
        let mut kin_instance = KineticData::new();
        let lib = "NUIG";

        // collecting reaction data for library name lib
        let _ = kin_instance.open_json_files(lib);
        // veiew all reactions in library
        let _ = kin_instance.print_all_reactions();

        // search reaction by equation
        let data_of_reacion = kin_instance.LibKineticData.get("1").unwrap();
        println!(" data_of_reacion {:?}", data_of_reacion);
    }

    #[test]
    fn test_search_data_by_field() {
        let mut kin_instance = KineticData::new();
        let _ = kin_instance.open_json_files("NUIG");

        let reactions_by_field = kin_instance.get_reactions_by_field("type", "falloff");
        println!(" reactions_by_field {:?}", reactions_by_field);
        assert_eq!(reactions_by_field.len(), 109);
    }

    #[test]
    fn test_form_user_chosen_data() {
        let mut kin_instance = KineticData::new();
        let _ = kin_instance.open_json_files("NUIG");

        // Prepare some test data
        kin_instance.FoundReactionsByProducts = vec!["1".to_string(), "2".to_string()];
        kin_instance.FoundReactionsByReagents = vec!["3".to_string(), "4".to_string()];

        // Test with FoundReactionsByProducts and default mechanism name
        let result = kin_instance
            .form_user_chosen_data(SearchVariants::FoundReactionsByProducts, None)
            .unwrap();
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("NUIG"));
        let nuig_data = result.get("NUIG").unwrap();
        assert_eq!(nuig_data.len(), 2);
        assert!(nuig_data.contains_key("1"));
        assert!(nuig_data.contains_key("2"));

        // Test with FoundReactionsByReagents and specified mechanism name
        let result = kin_instance
            .form_user_chosen_data(
                SearchVariants::FoundReactionsByReagents,
                Some("CustomMech".to_string()),
            )
            .unwrap();
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("CustomMech"));
        let custom_mech_data = result.get("CustomMech").unwrap();
        assert_eq!(custom_mech_data.len(), 2);
        assert!(custom_mech_data.contains_key("3"));
        assert!(custom_mech_data.contains_key("4"));

        // Test with empty results
        kin_instance.FoundReactionsByProducts = vec![];
        let result = kin_instance
            .form_user_chosen_data(SearchVariants::FoundReactionsByProducts, None)
            .unwrap();
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_form_user_chosen_data_integration() {
        let mut kin_instance = KineticData::new();
        let _ = kin_instance.open_json_files("NUIG");

        // Perform a search to populate FoundReactionsByProducts
        kin_instance
            .search_reaction_by_reagents_and_products(vec!["O".to_string(), "H".to_string()])
            .unwrap();

        // Test with actual search results
        let result = kin_instance
            .form_user_chosen_data(SearchVariants::FoundReactionsByProducts, None)
            .unwrap();
        assert!(!result.is_empty());
        assert!(result.contains_key("NUIG"));
        let nuig_data = result.get("NUIG").unwrap();
        assert!(!nuig_data.is_empty());

        // Verify that the returned data matches the LibKineticData
        for (reaction_id, reaction_data) in nuig_data {
            assert_eq!(
                reaction_data,
                kin_instance.LibKineticData.get(reaction_id).unwrap()
            );
        }
    }

    #[test]
    fn test_cached_reactbase_reuses_same_allocation() {
        let first = cached_reactbase().unwrap();
        let second = cached_reactbase().unwrap();

        assert!(Arc::ptr_eq(&first, &second));
    }

    #[test]
    fn test_cached_reaction_db_reuses_same_allocation() {
        let first = cached_reaction_db().unwrap();
        let second = cached_reaction_db().unwrap();

        assert!(Arc::ptr_eq(&first, &second));
    }

    #[test]
    #[ignore = "full bundled-library stress story; included in the Kinetics campaign"]
    fn story_cantera_full_library_deserializes_into_typed_reactions() {
        use crate::Kinetics::mechfinder_api::parse_kinetic_data_vec;

        let mut catalog = KineticData::new();
        catalog
            .open_json_files("Cantera")
            .expect("bundled Cantera catalog should load");
        assert!(!catalog.LibKineticData.is_empty());

        let mut reaction_ids: Vec<_> = catalog.LibKineticData.keys().cloned().collect();
        reaction_ids.sort();
        let payloads: Vec<_> = reaction_ids
            .iter()
            .map(|reaction_id| {
                catalog
                    .LibKineticData
                    .get(reaction_id)
                    .expect("sorted reaction id came from this library")
                    .clone()
            })
            .collect();
        let raw_count = payloads.len();
        let (reactions, equations) = parse_kinetic_data_vec(payloads)
            .expect("every Cantera reaction payload should deserialize and validate");

        assert_eq!(reactions.len(), raw_count);
        assert_eq!(equations.len(), raw_count);
        println!(
            "Cantera full-library parse: payloads={raw_count}, typed={}, equations={}",
            reactions.len(),
            equations.len()
        );
    }
    /*
    #[test]
    fn test_kinetic_data_open_json_files() {
        let mut kinetic_data = KineticData::new();

        // Replace "big_mech" with the actual file name
        kinetic_data.open_json_files("big_mech");

        // Add assertions to verify the content of kinetic_data.HashMapOfReactantsAndProducts,
        // kinetic_data.LibKineticData, kinetic_data.AllLibraries, etc.
    }



    */
    // Add more tests for other methods as needed
}
