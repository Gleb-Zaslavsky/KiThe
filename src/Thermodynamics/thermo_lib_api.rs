use crate::Thermodynamics::physical_state::{
    PhysicalState, PhysicalStateEvidence, ThermoRecordQuery,
};
use crate::library_manager::with_library_manager;
use prettytable::{Table, row};
use serde::de::DeserializeOwned;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;
use std::sync::{Arc, OnceLock};

#[derive(Debug)]
pub enum ThermoLibraryError {
    MissingLibraryFile {
        path: String,
        source: std::io::Error,
    },
    InvalidLibraryJson {
        path: String,
        source: serde_json::Error,
    },
    EmptyLibraryCatalog {
        path: String,
    },
    InvalidLibrarySchema {
        library: String,
        reason: String,
    },
    AmbiguousPhysicalState {
        library: String,
        substance: String,
        requested: PhysicalState,
        candidates: Vec<String>,
    },
    ConflictingPhysicalStateRequest {
        substance: String,
        requested: PhysicalState,
        encoded: PhysicalState,
    },
}

impl std::fmt::Display for ThermoLibraryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ThermoLibraryError::MissingLibraryFile { path, source } => {
                write!(
                    f,
                    "failed to open thermo library file '{}': {}",
                    path, source
                )
            }
            ThermoLibraryError::InvalidLibraryJson { path, source } => {
                write!(
                    f,
                    "failed to parse thermo library file '{}': {}",
                    path, source
                )
            }
            ThermoLibraryError::EmptyLibraryCatalog { path } => {
                write!(f, "thermo library catalog '{}' is empty", path)
            }
            ThermoLibraryError::InvalidLibrarySchema { library, reason } => {
                write!(
                    f,
                    "invalid thermo library schema for '{}': {}",
                    library, reason
                )
            }
            ThermoLibraryError::AmbiguousPhysicalState {
                library,
                substance,
                requested,
                candidates,
            } => write!(
                f,
                "physical-state lookup for '{}' in '{}' is ambiguous for {:?}: {}",
                substance,
                library,
                requested,
                candidates.join(", ")
            ),
            ThermoLibraryError::ConflictingPhysicalStateRequest {
                substance,
                requested,
                encoded,
            } => write!(
                f,
                "physical-state request for '{}' conflicts: requested {:?}, key encodes {:?}",
                substance, requested, encoded
            ),
        }
    }
}

impl Error for ThermoLibraryError {}

/// Coarse-grained capability of a library record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LibraryCapability {
    Thermo,
    Transport,
}

/// Canonical identifiers for the library records we know how to classify.
///
/// The enum deliberately covers the repository keys and their accepted aliases
/// so the rest of the codebase can stop repeating string matching logic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LibraryId {
    CanteraNasaBaseGas,
    CanteraNasaBaseCond,
    Nasa,
    Nasa7,
    NasaGas,
    NasaCond,
    NuigThermo,
    Cea,
    Nist,
    Nist9,
    AramcoTransport,
}

impl LibraryId {
    /// Returns the canonical storage key used by the catalog.
    pub const fn canonical_name(self) -> &'static str {
        match self {
            Self::CanteraNasaBaseGas => "NASA_gas",
            Self::CanteraNasaBaseCond => "NASA_cond",
            Self::Nasa => "NASA",
            Self::Nasa7 => "NASA7",
            Self::NasaGas => "NASA_gas",
            Self::NasaCond => "NASA_cond",
            Self::NuigThermo => "nuig_thermo",
            Self::Cea => "CEA",
            Self::Nist => "NIST",
            Self::Nist9 => "NIST9",
            Self::AramcoTransport => "Aramco_transport",
        }
    }

    /// Returns the semantic capability of the library.
    pub const fn capability(self) -> LibraryCapability {
        match self {
            Self::Cea | Self::AramcoTransport => LibraryCapability::Transport,
            _ => LibraryCapability::Thermo,
        }
    }

    /// Returns the family-specific backend name used by the calculator factory.
    pub const fn handler_name(self) -> &'static str {
        match self.capability() {
            LibraryCapability::Thermo => match self {
                Self::Nist | Self::Nist9 => "NIST",
                _ => "NASA",
            },
            LibraryCapability::Transport => match self {
                Self::Cea => "CEA",
                Self::AramcoTransport => "transport",
                _ => "transport",
            },
        }
    }

    /// Resolves accepted aliases to a typed library identifier.
    pub fn resolve(lib_name: &str) -> Option<Self> {
        match lib_name {
            "Cantera_nasa_base_gas" => Some(Self::CanteraNasaBaseGas),
            "Cantera_nasa_base_cond" => Some(Self::CanteraNasaBaseCond),
            "NASA" => Some(Self::Nasa),
            "NASA7" => Some(Self::Nasa7),
            "NASA_gas" => Some(Self::NasaGas),
            "NASA_cond" => Some(Self::NasaCond),
            "nuig_thermo" | "NUIG_thermo" => Some(Self::NuigThermo),
            "CEA" => Some(Self::Cea),
            "NIST" => Some(Self::Nist),
            "NIST9" => Some(Self::Nist9),
            "Aramco_transpot" | "Aramco_transport" => Some(Self::AramcoTransport),
            _ => None,
        }
    }
}

/// Immutable thermodynamic catalog shared between `ThermoData` instances.
///
/// The repository owns the loaded library metadata and database payloads.
/// `ThermoData` keeps a handle to this catalog and maintains only query-side
/// mutable state locally.
#[derive(Debug)]
pub struct ThermoRepository {
    pub VecOfSubsAdresses: Arc<Vec<(String, String)>>,
    pub LibThermoData: Arc<HashMap<String, HashMap<String, Value>>>,
    pub ElementsData: Arc<HashMap<String, Vec<Vec<String>>>>,
    pub AllLibraries: Arc<Vec<String>>,
    pub library_aliases: Arc<HashMap<String, String>>,
    pub library_handlers: Arc<HashMap<String, String>>,
    pub thermo_libs: Arc<Vec<String>>,
    pub transport_libs: Arc<Vec<String>>,
}

/// A record selected from a catalog with explicit lookup provenance.
#[derive(Debug, Clone)]
pub struct ResolvedThermoRecord {
    pub record_key: String,
    pub data: Value,
    pub physical_state: Option<PhysicalState>,
    pub state_evidence: Option<PhysicalStateEvidence>,
}

impl ThermoRepository {
    fn empty() -> Self {
        Self {
            VecOfSubsAdresses: Arc::new(Vec::new()),
            LibThermoData: Arc::new(HashMap::new()),
            ElementsData: Arc::new(HashMap::new()),
            AllLibraries: Arc::new(Vec::new()),
            library_aliases: Arc::new(HashMap::new()),
            library_handlers: Arc::new(HashMap::new()),
            thermo_libs: Arc::new(Vec::new()),
            transport_libs: Arc::new(Vec::new()),
        }
    }

    pub(crate) fn from_parts(
        VecOfSubsAdresses: Vec<(String, String)>,
        LibThermoData: HashMap<String, HashMap<String, Value>>,
        ElementsData: HashMap<String, Vec<Vec<String>>>,
        AllLibraries: Vec<String>,
        library_aliases: HashMap<String, String>,
        library_handlers: HashMap<String, String>,
        thermo_libs: Vec<String>,
        transport_libs: Vec<String>,
    ) -> Self {
        Self {
            VecOfSubsAdresses: Arc::new(VecOfSubsAdresses),
            LibThermoData: Arc::new(LibThermoData),
            ElementsData: Arc::new(ElementsData),
            AllLibraries: Arc::new(AllLibraries),
            library_aliases: Arc::new(library_aliases),
            library_handlers: Arc::new(library_handlers),
            thermo_libs: Arc::new(thermo_libs),
            transport_libs: Arc::new(transport_libs),
        }
    }

    /// Resolves a library record while respecting a requested physical state.
    ///
    /// State suffixes are interpreted only for library families that document
    /// them. In particular, NUIG labels such as `CH2(S)` are left literal;
    /// treating every parenthesized suffix as a solid state would corrupt
    /// chemistry identifiers.
    pub fn resolve_thermo_record(
        &self,
        library: &str,
        query: &ThermoRecordQuery,
    ) -> Result<Option<ResolvedThermoRecord>, ThermoLibraryError> {
        resolve_thermo_record_from_catalog(self.LibThermoData.as_ref(), library, query)
    }
}

fn resolve_thermo_record_from_catalog(
    catalog: &HashMap<String, HashMap<String, Value>>,
    library: &str,
    query: &ThermoRecordQuery,
) -> Result<Option<ResolvedThermoRecord>, ThermoLibraryError> {
    let canonical_library = ThermoData::canonical_library_name(library);
    let Some(records) = catalog.get(&canonical_library) else {
        return Ok(None);
    };
    let library_id = LibraryId::resolve(&canonical_library);
    let (query_base, encoded_state) = split_query_state(library_id, query.substance());
    if let (Some(requested), Some(encoded)) = (query.physical_state(), encoded_state) {
        if !requested.accepts(encoded) {
            return Err(ThermoLibraryError::ConflictingPhysicalStateRequest {
                substance: query.substance().to_string(),
                requested,
                encoded,
            });
        }
    }
    let requested_state = query.physical_state().or(encoded_state);

    if requested_state.is_none() {
        if let Some(data) = records.get(query.substance()) {
            let (state, evidence) = record_state(library_id, query.substance());
            return Ok(Some(ResolvedThermoRecord {
                record_key: query.substance().to_string(),
                data: data.clone(),
                physical_state: state,
                state_evidence: evidence,
            }));
        }
        return Ok(None);
    }

    let Some(requested_state) = requested_state else {
        return Ok(None);
    };
    let mut candidates = Vec::new();
    for (record_key, data) in records.iter() {
        let (record_base, _) = split_query_state(library_id, record_key);
        let (state, evidence) = record_state(library_id, record_key);
        if record_base == query_base
            && state.is_some_and(|observed| requested_state.accepts(observed))
        {
            candidates.push(ResolvedThermoRecord {
                record_key: record_key.clone(),
                data: data.clone(),
                physical_state: state,
                state_evidence: evidence,
            });
        }
    }
    match candidates.len() {
        0 => Ok(None),
        1 => Ok(candidates.pop()),
        _ => {
            let mut candidate_keys: Vec<_> = candidates
                .into_iter()
                .map(|candidate| candidate.record_key)
                .collect();
            candidate_keys.sort();
            Err(ThermoLibraryError::AmbiguousPhysicalState {
                library: canonical_library,
                substance: query_base.to_string(),
                requested: requested_state,
                candidates: candidate_keys,
            })
        }
    }
}

fn split_query_state<'a>(
    library: Option<LibraryId>,
    key: &'a str,
) -> (&'a str, Option<PhysicalState>) {
    let Some((base, state)) = state_suffix(library, key) else {
        return (key, None);
    };
    (base, Some(state))
}

fn record_state(
    library: Option<LibraryId>,
    key: &str,
) -> (Option<PhysicalState>, Option<PhysicalStateEvidence>) {
    if let Some((_, state)) = state_suffix(library, key) {
        return (Some(state), Some(PhysicalStateEvidence::KeyConvention));
    }
    match library {
        Some(LibraryId::NasaGas | LibraryId::CanteraNasaBaseGas) => (
            Some(PhysicalState::Gas),
            Some(PhysicalStateEvidence::LibraryDefault),
        ),
        Some(LibraryId::NasaCond | LibraryId::CanteraNasaBaseCond) => (
            Some(PhysicalState::Condensed),
            Some(PhysicalStateEvidence::LibraryDefault),
        ),
        _ => (None, None),
    }
}

/// Applies only documented, family-specific state conventions.
fn state_suffix(library: Option<LibraryId>, key: &str) -> Option<(&str, PhysicalState)> {
    let supports_state_suffix = matches!(
        library,
        Some(
            LibraryId::NasaCond
                | LibraryId::CanteraNasaBaseCond
                | LibraryId::Nist
                | LibraryId::Nist9
        )
    );
    if !supports_state_suffix || !key.ends_with(')') {
        return None;
    }
    let open = key.rfind('(')?;
    let base = &key[..open];
    let tag = &key[open + 1..key.len() - 1];
    let state = match tag {
        "g" | "G" => PhysicalState::Gas,
        "l" | "L" => PhysicalState::Liquid,
        "s" | "S" => PhysicalState::Solid,
        // NASA-condensed uses alphabetical allotrope labels for records such
        // as Fe(a), Fe(c), and Fe(d). They are all solids but the exact
        // allotrope must remain explicit when more than one exists.
        "a" | "b" | "c" | "d" => PhysicalState::Solid,
        _ => return None,
    };
    Some((base, state))
}

static DEFAULT_THERMO_REPOSITORY: OnceLock<Arc<ThermoRepository>> = OnceLock::new();

#[derive(Debug, Clone)]
pub struct ThermoData {
    /// Shared immutable catalog handle for all read-only library data.
    pub repository: Arc<ThermoRepository>,
    pub VecOfSubsAdresses: Arc<Vec<(String, String)>>,
    /// adresses of all substances (lib:sustance)
    pub LibThermoData: Arc<HashMap<String, HashMap<String, Value>>>,
    /// library of thermodymical and heat mass transfer data  lib:{ substance :{ data }}
    pub ElementsData: Arc<HashMap<String, Vec<Vec<String>>>>,
    /// elements base data {"element_name": [["substance_name", "database_name"], ...]}
    pub AllLibraries: Arc<Vec<String>>,
    pub library_handlers: Arc<HashMap<String, String>>,
    /// all sub-libraries in big library
    pub subs_to_search: Vec<String>,
    /// all substances to search in libraries
    pub hashmap_of_thermo_data: HashMap<String, HashMap<String, Value>>,
    /// search result {substance: {lib: data}}
    pub all_substances: Vec<String>,
    pub thermo_libs: Arc<Vec<String>>,
    pub transport_libs: Arc<Vec<String>>,
}
impl ThermoData {
    fn default_library_aliases() -> Arc<HashMap<String, String>> {
        static DEFAULT_LIBRARY_ALIASES: OnceLock<Arc<HashMap<String, String>>> = OnceLock::new();
        if let Some(registry) = DEFAULT_LIBRARY_ALIASES.get() {
            return Arc::clone(registry);
        }

        let registry = Arc::new(HashMap::from([
            ("Cantera_nasa_base_gas".to_string(), "NASA_gas".to_string()),
            (
                "Cantera_nasa_base_cond".to_string(),
                "NASA_cond".to_string(),
            ),
        ]));

        let _ = DEFAULT_LIBRARY_ALIASES.set(Arc::clone(&registry));
        registry
    }

    fn default_library_handlers() -> Arc<HashMap<String, String>> {
        static DEFAULT_LIBRARY_HANDLERS: OnceLock<Arc<HashMap<String, String>>> = OnceLock::new();
        if let Some(registry) = DEFAULT_LIBRARY_HANDLERS.get() {
            return Arc::clone(registry);
        }

        let registry = Arc::new(HashMap::from([
            ("Cantera_nasa_base_gas".to_string(), "NASA".to_string()),
            ("Cantera_nasa_base_cond".to_string(), "NASA".to_string()),
        ]));

        let _ = DEFAULT_LIBRARY_HANDLERS.set(Arc::clone(&registry));
        registry
    }

    pub fn canonical_library_name(lib_name: &str) -> String {
        LibraryId::resolve(lib_name)
            .map(|id| id.canonical_name().to_string())
            .or_else(|| Self::default_library_aliases().get(lib_name).cloned())
            .unwrap_or_else(|| lib_name.to_string())
    }

    /// Resolves against this query object's current catalog view.
    ///
    /// `ThermoData` normally shares immutable data with a repository, but a
    /// caller may intentionally detach its `Arc` for an isolated transaction
    /// or test fixture. Lookup must honor that local view.
    pub fn resolve_thermo_record(
        &self,
        library: &str,
        query: &ThermoRecordQuery,
    ) -> Result<Option<ResolvedThermoRecord>, ThermoLibraryError> {
        resolve_thermo_record_from_catalog(self.LibThermoData.as_ref(), library, query)
    }

    /// Returns the typed capability for a library record if it is recognized.
    pub fn library_capability(lib_name: &str) -> Option<LibraryCapability> {
        LibraryId::resolve(lib_name)
            .map(|id| id.capability())
            .or_else(|| {
                match Self::default_library_handlers()
                    .get(lib_name)
                    .map(|s| s.as_str())
                {
                    Some("NASA") => Some(LibraryCapability::Thermo),
                    Some("NIST") => Some(LibraryCapability::Thermo),
                    Some("CEA") => Some(LibraryCapability::Transport),
                    Some("transport") => Some(LibraryCapability::Transport),
                    _ => None,
                }
            })
    }

    fn classify_library(lib_name: &str) -> String {
        LibraryId::resolve(lib_name)
            .map(|id| id.handler_name().to_string())
            .or_else(|| Self::default_library_handlers().get(lib_name).cloned())
            .unwrap_or_else(|| lib_name.to_string())
    }

    fn split_libraries_by_handler(
        lib_thermo_data: &HashMap<String, HashMap<String, Value>>,
    ) -> (Vec<String>, Vec<String>) {
        let mut thermo_libs = Vec::new();
        let mut transport_libs = Vec::new();

        for lib_name in lib_thermo_data.keys() {
            match Self::library_capability(lib_name) {
                Some(LibraryCapability::Thermo) => thermo_libs.push(lib_name.clone()),
                Some(LibraryCapability::Transport) => transport_libs.push(lib_name.clone()),
                None => {}
            }
        }

        thermo_libs.sort();
        transport_libs.sort();
        (thermo_libs, transport_libs)
    }

    fn validate_loaded_catalog(
        all_keys_path: &str,
        substance_base_path: &str,
        elements_path: &str,
        vec_of_subs_adresses: &[(String, String)],
        lib_thermo_data: &HashMap<String, HashMap<String, Value>>,
        elements_data: &HashMap<String, Vec<Vec<String>>>,
    ) -> Result<(), ThermoLibraryError> {
        if vec_of_subs_adresses.is_empty() {
            return Err(ThermoLibraryError::EmptyLibraryCatalog {
                path: all_keys_path.to_string(),
            });
        }
        if lib_thermo_data.is_empty() {
            return Err(ThermoLibraryError::EmptyLibraryCatalog {
                path: substance_base_path.to_string(),
            });
        }
        if elements_data.is_empty() {
            return Err(ThermoLibraryError::EmptyLibraryCatalog {
                path: elements_path.to_string(),
            });
        }

        Ok(())
    }

    fn load_json_from_path<T: DeserializeOwned>(path: &str) -> Result<T, ThermoLibraryError> {
        let mut file =
            File::open(path).map_err(|source| ThermoLibraryError::MissingLibraryFile {
                path: path.to_string(),
                source,
            })?;
        let mut file_contents = String::new();
        file.read_to_string(&mut file_contents).map_err(|source| {
            ThermoLibraryError::MissingLibraryFile {
                path: path.to_string(),
                source,
            }
        })?;
        serde_json::from_str(&file_contents).map_err(|source| {
            ThermoLibraryError::InvalidLibraryJson {
                path: path.to_string(),
                source,
            }
        })
    }

    fn try_load_repository_from_default_paths() -> Result<ThermoRepository, ThermoLibraryError> {
        let all_keys_path =
            with_library_manager(|manager| manager.all_keys_substance_path().to_string());
        let substance_base_path =
            with_library_manager(|manager| manager.substance_base_path().to_string());
        let elements_path = with_library_manager(|manager| manager.elements_path().to_string());

        let lib_sub_pairs_vec: Vec<(String, String)> = Self::load_json_from_path(&all_keys_path)?;
        let lib_thermo_data: HashMap<String, HashMap<String, Value>> =
            Self::load_json_from_path(&substance_base_path)?;
        let elements_data: HashMap<String, Vec<Vec<String>>> =
            Self::load_json_from_path(&elements_path)?;

        Self::validate_loaded_catalog(
            &all_keys_path,
            &substance_base_path,
            &elements_path,
            &lib_sub_pairs_vec,
            &lib_thermo_data,
            &elements_data,
        )?;

        let all_libraries: Vec<String> = lib_thermo_data.keys().cloned().collect();
        let (thermo_libs, transport_libs) = Self::split_libraries_by_handler(&lib_thermo_data);

        Ok(ThermoRepository::from_parts(
            lib_sub_pairs_vec,
            lib_thermo_data,
            elements_data,
            all_libraries,
            (*Self::default_library_aliases()).clone(),
            (*Self::default_library_handlers()).clone(),
            thermo_libs,
            transport_libs,
        ))
    }

    pub fn try_new_from_paths(
        all_keys_path: &str,
        substance_base_path: &str,
        elements_path: &str,
    ) -> Result<Self, ThermoLibraryError> {
        let vec_of_subs_adresses: Vec<(String, String)> = Self::load_json_from_path(all_keys_path)?;
        let lib_thermo_data: HashMap<String, HashMap<String, Value>> =
            Self::load_json_from_path(substance_base_path)?;
        let elements_data: HashMap<String, Vec<Vec<String>>> =
            Self::load_json_from_path(elements_path)?;
        Self::validate_loaded_catalog(
            all_keys_path,
            substance_base_path,
            elements_path,
            &vec_of_subs_adresses,
            &lib_thermo_data,
            &elements_data,
        )?;
        let all_libraries: Vec<String> = lib_thermo_data.keys().cloned().collect();
        let (thermo_libs, transport_libs) = Self::split_libraries_by_handler(&lib_thermo_data);

        let repository = ThermoRepository::from_parts(
            vec_of_subs_adresses,
            lib_thermo_data,
            elements_data,
            all_libraries,
            (*Self::default_library_aliases()).clone(),
            (*Self::default_library_handlers()).clone(),
            thermo_libs,
            transport_libs,
        );
        Ok(Self::from_repository(Arc::new(repository)))
    }

    fn default_repository() -> Result<Arc<ThermoRepository>, ThermoLibraryError> {
        if let Some(repository) = DEFAULT_THERMO_REPOSITORY.get() {
            return Ok(Arc::clone(repository));
        }

        let repository = Arc::new(Self::try_load_repository_from_default_paths()?);
        let _ = DEFAULT_THERMO_REPOSITORY.set(Arc::clone(&repository));
        Ok(repository)
    }

    /// Returns the process-wide immutable catalog used by default lookups.
    ///
    /// Callers that resolve a compound object, such as a multi-phase system,
    /// can pass this handle through explicitly so every child query is tied to
    /// the same repository lifecycle and no hidden catalog load is involved.
    pub fn try_default_repository() -> Result<Arc<ThermoRepository>, ThermoLibraryError> {
        Self::default_repository()
    }

    pub fn try_new() -> Result<Self, ThermoLibraryError> {
        let repository = Self::default_repository()?;
        Ok(Self::from_repository(repository))
    }

    pub fn try_new_fresh() -> Result<Self, ThermoLibraryError> {
        let repository = Arc::new(Self::try_load_repository_from_default_paths()?);
        Ok(Self::from_repository(repository))
    }

    pub fn from_repository(repository: Arc<ThermoRepository>) -> Self {
        Self {
            repository: Arc::clone(&repository),
            VecOfSubsAdresses: Arc::clone(&repository.VecOfSubsAdresses),
            LibThermoData: Arc::clone(&repository.LibThermoData),
            ElementsData: Arc::clone(&repository.ElementsData),
            AllLibraries: Arc::clone(&repository.AllLibraries),
            library_handlers: Arc::clone(&repository.library_handlers),
            all_substances: Vec::new(),
            subs_to_search: Vec::new(),
            hashmap_of_thermo_data: HashMap::new(),
            thermo_libs: Arc::clone(&repository.thermo_libs),
            transport_libs: Arc::clone(&repository.transport_libs),
        }
    }

    /// Compatibility constructor kept for existing callers.
    /// Prefer `try_new()` when the caller wants explicit file-loading errors.
    pub fn new() -> Self {
        Self::try_new().expect(
            "ThermoData::new() requires the bundled thermo catalog; use try_new() or try_new_from_paths() when loading may fail",
        )
    }

    pub fn empty() -> Self {
        Self::from_repository(Arc::new(ThermoRepository::empty()))
    }
    //////////////////////////////////////////////DATA MANIPULATION//////////////////////////////////////////////////////////
    /// add substances in the form {substanse:{type_of_data*:{data} }}
    /// *type of data may be 1) NASA format of thermodynamic data 2) NIST format of thermodynamic data 3) "transport" format of heat mass transfer data
    pub fn append_substances_data(
        &mut self,
        hashmap_of_thermo_data: HashMap<String, HashMap<String, Value>>,
    ) {
        let new_keys: Vec<String> = hashmap_of_thermo_data.keys().cloned().collect();
        self.hashmap_of_thermo_data.extend(hashmap_of_thermo_data);
        self.subs_to_search.extend(new_keys);
    }
    pub fn subs_of_this_lib(&mut self, lib: &str) -> Vec<String> {
        let subs_of_this_lib: Vec<String> = self
            .VecOfSubsAdresses
            .iter()
            .filter(|(k, _)| k == lib)
            .map(|(_, v)| v.clone())
            .collect();
        subs_of_this_lib
    }

    ///  remove substance from hashmap and vector of substances
    pub fn remove_by_name(&mut self, name: String) -> bool {
        let removed_from_search =
            if let Some(position) = self.subs_to_search.iter().position(|x| x == &name) {
                self.subs_to_search.remove(position);
                true
            } else {
                false
            };
        let removed_from_map = self.hashmap_of_thermo_data.remove(&name).is_some();
        removed_from_search || removed_from_map
    }
    /////////////////////////////////PARSE DATA FROM LIBRARIES////////////////////////////////////////////////
    /// Search for a specific substance in all reaction libraries and return all libraries where it is found
    pub fn search_libs(
        &mut self,
        substance: &str,
        allowed_libs: Option<Vec<String>>,
    ) -> Vec<String> {
        let allowed_libs = allowed_libs.map(|libs| libs.into_iter().collect::<HashSet<_>>());
        self.search_libs_with_allowed_set(substance, allowed_libs.as_ref())
    }

    fn search_libs_with_allowed_set(
        &self,
        substance: &str,
        allowed_libs: Option<&HashSet<String>>,
    ) -> Vec<String> {
        self.VecOfSubsAdresses
            .iter()
            .filter(|(_, subs)| subs == substance)
            .filter(|(library, _)| {
                allowed_libs
                    .map(|allowed| allowed.contains(library.as_str()))
                    .unwrap_or(true)
            })
            .map(|(library, _)| library.clone())
            .collect()
    }
    /// search a substances from vector in all libraries and making a hashmap of information in the form {substanse:{type_of_data*:{data} }}
    /// *type of data may be 1) NASA format of thermodynamic data 2) NIST format of thermodynamic data 3) "transport" format of heat mass transfer data
    pub fn search_libs_for_subs(
        &mut self,
        substances: Vec<String>,
        allowed_libs: Option<Vec<String>>,
    ) -> Vec<String> {
        let allowed_libs = allowed_libs.map(|libs| libs.into_iter().collect::<HashSet<_>>());
        let mut hashmap_of_thermo_data: HashMap<String, HashMap<String, Value>> = HashMap::new();
        let mut libs_for_this_subs = Vec::new();

        for sub in substances {
            // for this sub let us find all libraries where this sub is present
            let libs_for_this_sub: Vec<String> =
                self.search_libs_with_allowed_set(&sub, allowed_libs.as_ref());
            let hashmap_lib_data: HashMap<String, Value> = libs_for_this_sub // iterating through all libraries where this sub is present and finding data for this sub
                .iter()
                .filter_map(|lib| {
                    // search data for this sub in this libraries
                    let canonical_lib = Self::canonical_library_name(lib);
                    self.LibThermoData
                        .get(&canonical_lib) // get data for this sub from this library and save it in hashmap {lib:substance}
                        .and_then(|lib_data| lib_data.get(&sub).map(|v| ((lib.clone()), v.clone())))
                })
                .collect();
            hashmap_of_thermo_data.insert(sub, hashmap_lib_data);
            libs_for_this_subs.extend(libs_for_this_sub);
        }

        self.hashmap_of_thermo_data = hashmap_of_thermo_data;
        libs_for_this_subs
    }
    pub fn rename(lib_name: &str) -> &str {
        match lib_name {
            "Cantera_nasa_base_gas" => "NASA_gas",
            "Cantera_nasa_base_cond" => "NASA_cond",
            _ => lib_name,
        }
    }

    /// Get elements data
    pub fn get_elements_data(&self) -> &HashMap<String, Vec<Vec<String>>> {
        self.ElementsData.as_ref()
    }

    /// Set elements data
    pub fn set_elements_data(&mut self, elements_data: HashMap<String, Vec<Vec<String>>>) {
        self.ElementsData = Arc::new(elements_data);
    }

    /// Search substances by elements and populate hashmap_of_thermo_data
    pub fn search_by_elements(&mut self, elements: Vec<String>) -> Vec<String> {
        let mut found_substances = Vec::new();
        let mut seen_substances = HashSet::new();
        let mut hashmap_of_thermo_data: HashMap<String, HashMap<String, Value>> = HashMap::new();

        for element in elements {
            if let Some(substance_pairs) = self.ElementsData.get(&element) {
                for pair in substance_pairs {
                    if pair.len() >= 2 {
                        let substance_name = &pair[0];
                        let database_name = &pair[1];

                        if let Some(lib_data) = self.LibThermoData.get(database_name) {
                            if let Some(substance_data) = lib_data.get(substance_name) {
                                hashmap_of_thermo_data
                                    .entry(substance_name.clone())
                                    .or_insert_with(HashMap::new)
                                    .insert(database_name.clone(), substance_data.clone());

                                if seen_substances.insert(substance_name.clone()) {
                                    found_substances.push(substance_name.clone());
                                }
                            }
                        }
                    }
                }
            }
        }

        self.hashmap_of_thermo_data.extend(hashmap_of_thermo_data);
        found_substances
    }

    /// Search substances containing only the specified elements or subsets
    pub fn search_by_elements_only(&mut self, elements: Vec<String>) -> Vec<String> {
        let mut found_substances = Vec::new();
        let mut seen_substances = HashSet::new();
        let mut hashmap_of_thermo_data: HashMap<String, HashMap<String, Value>> = HashMap::new();
        let allowed_elements: HashSet<String> = elements.into_iter().collect();

        // Collect all substances and their elements
        let mut substance_elements: HashMap<String, (HashSet<String>, HashSet<String>)> =
            HashMap::new();

        let mut element_keys: Vec<String> = self.ElementsData.keys().cloned().collect();
        element_keys.sort();

        for element in element_keys {
            let Some(substance_pairs) = self.ElementsData.get(&element) else {
                continue;
            };
            for pair in substance_pairs {
                if pair.len() >= 2 {
                    let substance_name = &pair[0];
                    let database_name = &pair[1];

                    substance_elements
                        .entry(substance_name.clone())
                        .or_insert_with(|| (HashSet::new(), HashSet::new()))
                        .0
                        .insert(database_name.clone());
                    substance_elements
                        .entry(substance_name.clone())
                        .or_insert_with(|| (HashSet::new(), HashSet::new()))
                        .1
                        .insert(element.clone());
                }
            }
        }

        // Filter substances that contain only allowed elements (or subsets)
        let mut substance_names: Vec<String> = substance_elements.keys().cloned().collect();
        substance_names.sort();

        for substance_name in substance_names {
            let Some((database_names, substance_element_set)) =
                substance_elements.remove(&substance_name)
            else {
                continue;
            };
            if substance_element_set.is_subset(&allowed_elements) {
                let mut libraries: Vec<String> = database_names.into_iter().collect();
                libraries.sort();
                if let Some(database_name) = libraries.first() {
                    if let Some(lib_data) = self.LibThermoData.get(database_name) {
                        if let Some(substance_data) = lib_data.get(&substance_name) {
                            hashmap_of_thermo_data
                                .entry(substance_name.clone())
                                .or_insert_with(HashMap::new)
                                .insert(database_name.clone(), substance_data.clone());

                            if seen_substances.insert(substance_name.clone()) {
                                found_substances.push(substance_name);
                            }
                        }
                    }
                }
            }
        }

        self.hashmap_of_thermo_data.extend(hashmap_of_thermo_data);
        found_substances
    }
    ///////////////////INPUT/OUTPUT/////////////////////////////////////////////////
    ///
    pub fn save_thermo_data_to_json<P: AsRef<Path>>(
        &self,
        filename: P,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let filename = filename.as_ref();
        let json = serde_json::to_string(&self.hashmap_of_thermo_data)?;
        let mut file = File::create(filename)?;
        file.write_all(json.as_bytes())?;
        Ok(())
    }

    pub fn load_thermo_data_from_json<P: AsRef<Path>>(
        &mut self,
        filename: P,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let filename = filename.as_ref();
        let mut file = File::open(filename)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        self.hashmap_of_thermo_data = serde_json::from_str(&contents)?;
        self.subs_to_search = self
            .hashmap_of_thermo_data
            .keys()
            .map(|k| k.to_string())
            .collect();
        Ok(())
    }
    pub fn pretty_print_thermo_data(&self) {
        let data = self.hashmap_of_thermo_data.clone();
        let mut table = Table::new();
        table.add_row(row!["Substance", "Library", "Data",]);
        for (sub, hashmap_lib_data) in data.iter() {
            for (lib, data) in hashmap_lib_data.iter() {
                table.add_row(row![sub, lib, data]);
            }
        }
        table.printstd();
    }

    pub fn create_substance_document<P: AsRef<Path>>(
        &self,
        file_name: P,
    ) -> std::io::Result<HashMap<String, HashMap<String, Value>>> {
        let file_name = file_name.as_ref();
        let path = file_name;
        let file_exists = path.exists();

        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .append(true)
            .open(file_name)?;

        if !file_exists {
            writeln!(file, "SUBSTANCES DATA")?;
        } else {
            let reader = BufReader::new(File::open(file_name)?);
            let has_header = reader.lines().any(|line| {
                line.map(|line| {
                    let normalized = line.trim().to_uppercase();
                    normalized == "SUBSTANCES DATA" || normalized == "SUBS DATA"
                })
                .unwrap_or(false)
            });

            if !has_header {
                writeln!(file, "\nSUBSTANCES DATA")?;
            }
        }
        let hashmap_to_save: HashMap<String, HashMap<String, Value>> =
            self.hashmap_of_thermo_data.clone();
        if hashmap_to_save.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "hashmap_to_save is empty",
            ));
        }

        writeln!(file, "{}", serde_json::to_string_pretty(&hashmap_to_save)?)?;

        Ok(hashmap_to_save)
    }
    pub fn what_handler_to_use(lib_name: &str) -> String {
        Self::classify_library(lib_name)
    }
    ///////////////////////////////////////////////////////////////////////////
    pub fn compare_2_libs(
        &self,
    ) -> (
        Vec<(Option<String>, Option<String>)>,
        Vec<(Option<String>, Option<String>)>,
    ) {
        let vec_lib_subs = self.VecOfSubsAdresses.as_ref();
        let LibThermoData = self.LibThermoData.as_ref();
        let vec_lib_subs_set: HashSet<(&str, &str)> = vec_lib_subs
            .iter()
            .map(|(library, substance)| (library.as_str(), substance.as_str()))
            .collect();
        let mut present_in_lib_sub_pairs_but_not_in_big_lib: Vec<(Option<String>, Option<String>)> =
            Vec::new();
        for (library, substance) in vec_lib_subs.iter() {
            match LibThermoData.get(library) {
                Some(lib_data) => {
                    if !lib_data.contains_key(substance) {
                        present_in_lib_sub_pairs_but_not_in_big_lib
                            .push((Some(library.clone()), Some(substance.clone())));
                    }
                }
                None => {
                    present_in_lib_sub_pairs_but_not_in_big_lib.push((Some(library.clone()), None));
                }
            }
        }
        let mut present_in_big_lib_but_not_in_lib_sub_pairs: Vec<(Option<String>, Option<String>)> =
            Vec::new();
        for (library, lib_data) in LibThermoData.iter() {
            for substance in lib_data.keys() {
                if !vec_lib_subs_set.contains(&(library.as_str(), substance.as_str())) {
                    present_in_big_lib_but_not_in_lib_sub_pairs
                        .push((Some(library.clone()), Some(substance.clone())));
                }
            }
        }
        // Display results in tables
        let mut table1 = Table::new();
        table1.add_row(row!["Library", "Substance"]);
        for (lib, sub) in &present_in_lib_sub_pairs_but_not_in_big_lib {
            table1.add_row(row![
                lib.as_ref().unwrap_or(&"N/A".to_string()),
                sub.as_ref().unwrap_or(&"N/A".to_string())
            ]);
        }
        println!("Present in lib_sub_pairs but not in big_lib:");
        table1.printstd();

        let mut table2 = Table::new();
        table2.add_row(row!["Library", "Substance"]);
        for (lib, sub) in &present_in_big_lib_but_not_in_lib_sub_pairs {
            table2.add_row(row![
                lib.as_ref().unwrap_or(&"N/A".to_string()),
                sub.as_ref().unwrap_or(&"N/A".to_string())
            ]);
        }
        println!("\nPresent in big_lib but not in lib_sub_pairs:");
        table2.printstd();
        (
            present_in_lib_sub_pairs_but_not_in_big_lib.clone(),
            present_in_big_lib_but_not_in_lib_sub_pairs.clone(),
        )
    }

    pub fn compare_2_libs2(&self) -> (HashSet<String>, HashSet<String>) {
        let libs: HashSet<String> = self
            .VecOfSubsAdresses
            .iter()
            .map(|(lib, _subs)| lib.clone())
            .collect();
        let libs2: HashSet<String> = self.LibThermoData.keys().cloned().collect();
        (libs, libs2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use std::sync::Arc;

    use tempfile::{NamedTempFile, tempdir};
    #[test]
    fn test_load_json_from_path_rejects_invalid_json() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "{{ not valid json").unwrap();
        let path = temp_file.path().to_str().unwrap().to_string();

        let err = ThermoData::load_json_from_path::<Vec<(String, String)>>(&path).unwrap_err();
        assert!(matches!(err, ThermoLibraryError::InvalidLibraryJson { .. }));
    }

    #[test]
    fn test_load_json_from_path_rejects_missing_file() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path().to_str().unwrap().to_string();
        temp_file.close().unwrap();

        let err = ThermoData::load_json_from_path::<Vec<(String, String)>>(&path).unwrap_err();
        assert!(matches!(err, ThermoLibraryError::MissingLibraryFile { .. }));
    }

    #[test]
    fn test_remove_by_name_is_recoverable() {
        let mut thermo_data = ThermoData::empty();

        assert!(!thermo_data.remove_by_name("missing".to_string()));
        assert!(thermo_data.subs_to_search.is_empty());
    }

    #[test]
    fn test_search_libs_filters_disallowed_entries() {
        let mut thermo_data = ThermoData::empty();
        thermo_data.VecOfSubsAdresses = Arc::new(vec![
            ("NASA_gas".to_string(), "CO2".to_string()),
            ("NIST".to_string(), "CO2".to_string()),
            ("CEA".to_string(), "H2O".to_string()),
        ]);

        let libs = thermo_data.search_libs("CO2", Some(vec!["NASA_gas".to_string()]));
        assert_eq!(libs, vec!["NASA_gas".to_string()]);
    }

    #[test]
    fn test_search_libs() {
        let thermo_data = ThermoData::new();
        let (
            _present_in_lib_sub_pairs_but_not_in_big_lib,
            present_in_big_lib_but_not_in_lib_sub_pairs,
        ) = thermo_data.compare_2_libs();
        println!(
            "present_in_lib_sub_pairs_but_not_in_big_lib: {:?}",
            present_in_big_lib_but_not_in_lib_sub_pairs
        );
        let (lib, lib2) = thermo_data.compare_2_libs2();
        println!("lib: {:?}", lib);
        println!("lib2: {:?}", lib2);
    }
    #[test]
    fn test_create_substance_document_new_file() {
        let mut thermo_data = ThermoData::new();
        let mut test_data = HashMap::new();
        let mut substance_data = HashMap::new();
        substance_data.insert(
            "test_lib".to_string(),
            serde_json::json!({"test_key": "test_value"}),
        );
        test_data.insert("test_substance".to_string(), substance_data);
        thermo_data.hashmap_of_thermo_data = test_data;

        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path();

        let result = thermo_data.create_substance_document(file_path);
        println!("result: {:#?}", result);
        assert!(result.is_ok());

        let mut file_content = String::new();
        temp_file
            .as_file()
            .read_to_string(&mut file_content)
            .unwrap();

        assert!(file_content.contains("SUBSTANCES DATA"));
        assert!(file_content.contains("test_substance"));
        assert!(file_content.contains("test_lib"));
        assert!(file_content.contains("test_key"));
        assert!(file_content.contains("test_value"));
    }
    #[test]
    fn NamedTempFiletest() {
        let text = "Brian was here. Briefly.";

        // Create a file inside of `env::temp_dir()`.
        let mut named_temp_file = NamedTempFile::new().unwrap();

        // Re-open it.
        let mut file2 = named_temp_file.reopen().unwrap();

        // Write some test data to the first handle.
        named_temp_file.write_all(text.as_bytes()).unwrap();

        // Read the test data using the second handle.
        let mut buf = String::new();
        file2.read_to_string(&mut buf).unwrap();
        assert_eq!(buf, text);
    }
    #[test]
    fn test_create_substance_document_existing_file() {
        let mut thermo_data = ThermoData::new();
        let mut test_data = HashMap::new();
        let mut substance_data = HashMap::new();
        substance_data.insert(
            "test_lib".to_string(),
            serde_json::json!({"test_key": "test_value"}),
        );
        test_data.insert("test_substance".to_string(), substance_data);
        thermo_data.hashmap_of_thermo_data = test_data;

        let mut temp_file = NamedTempFile::new().unwrap();
        // Re-open it.
        let mut file2 = temp_file.reopen().unwrap();
        writeln!(temp_file, "Existing content").unwrap();
        temp_file.flush().unwrap();
        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        println!("content \n {:#?}", file_content);
        assert!(file_content.contains("Existing content"));

        let file_path = temp_file.path();
        let result = thermo_data.create_substance_document(file_path);
        assert!(result.is_ok());

        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        println!("content \n \n {:#?}", file_content);
        //    assert!(file_content.contains("Existing content\n"));
        assert!(file_content.contains("SUBSTANCES DATA"));
        assert!(file_content.contains("test_substance"));
        assert!(file_content.contains("test_lib"));
        assert!(file_content.contains("test_key"));
        assert!(file_content.contains("test_value"));
    }

    #[test]
    fn test_create_substance_document_empty_data() {
        let thermo_data = ThermoData::new();
        let temp_file = NamedTempFile::new().unwrap();
        let file_path = temp_file.path();

        let result = thermo_data.create_substance_document(file_path);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().kind(), std::io::ErrorKind::InvalidInput);
    }
    #[test]
    fn test_with_real_data() {
        let mut td = ThermoData::new();
        let subdata = td.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = subdata.get("CO").unwrap();
        println!("CO_data: {}", CO_data);
        td.search_libs_for_subs(
            vec!["CO".to_owned()],
            Some(vec!["Cantera_nasa_base_gas".to_owned()]),
        );
        println!(
            "td.hashmap_of_thermo_data: {:#?}",
            td.hashmap_of_thermo_data
        );
        let mut temp_file = NamedTempFile::new().unwrap();
        // Re-open it.
        let mut file2 = temp_file.reopen().unwrap();
        writeln!(temp_file, "Existing content").unwrap();
        let file_path = temp_file.path();
        let result = td.create_substance_document(file_path);
        assert!(result.is_ok());
        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        println!("content \n \n {:#?}", file_content);
        assert!(file_content.contains("Existing content\n"));
        assert!(file_content.contains("SUBSTANCES DATA"));
    }
    #[test]
    fn test_with_real_data2() {
        let mut td = ThermoData::new();
        let subdata = td.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = subdata.get("CO").unwrap();
        println!("CO_data: {}", CO_data);
        td.search_libs_for_subs(
            vec!["CO".to_owned(), "O2".to_owned()],
            Some(vec![
                "Cantera_nasa_base_gas".to_owned(),
                "Aramco_transport".to_owned(),
            ]),
        );
        println!(
            "td.hashmap_of_thermo_data: {:#?}",
            td.hashmap_of_thermo_data
        );
        let mut temp_file = NamedTempFile::new().unwrap();
        // Re-open it.
        let mut file2 = temp_file.reopen().unwrap();
        writeln!(temp_file, "Existing content").unwrap();
        let file_path = temp_file.path();
        let result = td.create_substance_document(file_path);
        assert!(result.is_ok());
        let mut file_content = String::new();
        file2.read_to_string(&mut file_content).unwrap();
        // println!("content \n \n {:#?}", file_content);
        assert!(file_content.contains("Existing content\n"));
        assert!(file_content.contains("SUBSTANCES DATA"));
    }

    #[test]
    fn test_save_and_load_thermo_data_with_pathbuf() {
        let mut thermo_data = ThermoData::new();
        let mut test_data = HashMap::new();
        let mut substance_data = HashMap::new();
        substance_data.insert(
            "test_lib".to_string(),
            serde_json::json!({"test_key": "test_value"}),
        );
        test_data.insert("test_substance".to_string(), substance_data);
        thermo_data.hashmap_of_thermo_data = test_data;

        let temp_file = NamedTempFile::new().unwrap();
        let path_buf = temp_file.path().to_path_buf();

        thermo_data
            .save_thermo_data_to_json(&path_buf)
            .expect("save should accept path-like values");

        let mut loaded = ThermoData::new();
        loaded
            .load_thermo_data_from_json(path_buf)
            .expect("load should accept path-like values");

        assert_eq!(
            loaded.hashmap_of_thermo_data,
            thermo_data.hashmap_of_thermo_data
        );
    }

    #[test]
    fn test_iterate_all_libraries() {
        let thermo_data = ThermoData::new();

        for library in thermo_data.AllLibraries.as_ref().iter() {
            println!("Library: {}", library);

            if let Some(lib_data) = thermo_data.LibThermoData.get(library) {
                if let Some((substance, data)) = lib_data.iter().next() {
                    println!("  Sample substance: {}", substance);
                    println!("  Data: {}", data);
                }
            }
            println!();
        }
    }

    #[test]
    fn test_library_classification() {
        let thermo_data = ThermoData::new();

        println!("Thermo libraries: {:?}", thermo_data.thermo_libs);
        println!("Transport libraries: {:?}", thermo_data.transport_libs);
        println!("All libraries: {:?}", thermo_data.AllLibraries);
        assert!(!thermo_data.thermo_libs.is_empty() || !thermo_data.transport_libs.is_empty());

        for lib in thermo_data.thermo_libs.as_ref().iter() {
            let handler = ThermoData::what_handler_to_use(lib);
            assert!(
                matches!(handler.as_str(), "NASA" | "NIST"),
                "Thermo lib {} resolved to unexpected handler {}",
                lib,
                handler
            );
        }

        for lib in thermo_data.transport_libs.as_ref().iter() {
            let handler = ThermoData::what_handler_to_use(lib);
            assert!(
                matches!(handler.as_str(), "CEA" | "transport"),
                "Transport lib {} resolved to unexpected handler {}",
                lib,
                handler
            );
        }
    }

    #[test]
    fn test_handler_registry_and_alias_resolution() {
        assert_eq!(
            ThermoData::what_handler_to_use("Cantera_nasa_base_gas"),
            "NASA"
        );
        assert_eq!(ThermoData::what_handler_to_use("NASA7"), "NASA");
        assert_eq!(ThermoData::what_handler_to_use("NIST9"), "NIST");
        assert_eq!(
            ThermoData::what_handler_to_use("Aramco_transpot"),
            "transport"
        );
        assert_eq!(ThermoData::rename("Cantera_nasa_base_cond"), "NASA_cond");
        assert_eq!(
            ThermoData::canonical_library_name("Aramco_transpot"),
            "Aramco_transport"
        );
        assert_eq!(
            ThermoData::canonical_library_name("NUIG_thermo"),
            "nuig_thermo"
        );
        assert_eq!(
            ThermoData::library_capability("Aramco_transpot"),
            Some(LibraryCapability::Transport)
        );
        assert_eq!(
            ThermoData::library_capability("NASA_gas"),
            Some(LibraryCapability::Thermo)
        );
    }

    #[test]
    fn test_elements_data_loading() {
        let thermo_data = ThermoData::new();
        println!(
            "Elements data loaded: {} elements",
            thermo_data.ElementsData.len()
        );
    }

    #[test]
    fn test_search_by_elements() {
        let mut thermo_data = ThermoData::new();

        // Test with common elements
        let elements = vec!["C".to_string(), "O".to_string()];
        let found_substances = thermo_data.search_by_elements(elements);

        println!("Found substances: {:?}", found_substances);
        println!(
            "Thermo data entries: {}",
            thermo_data.hashmap_of_thermo_data.len()
        );

        // Verify that substances were found and data was populated
        if !found_substances.is_empty() {
            assert!(!thermo_data.hashmap_of_thermo_data.is_empty());
        }
    }

    #[test]
    fn test_elements_getters_setters() {
        let mut thermo_data = ThermoData::new();

        let _original_size = thermo_data.get_elements_data().len();

        let mut test_elements = HashMap::new();
        test_elements.insert(
            "H".to_string(),
            vec![vec!["H2O".to_string(), "test_lib".to_string()]],
        );

        thermo_data.set_elements_data(test_elements);

        assert_eq!(thermo_data.get_elements_data().len(), 1);
        assert!(thermo_data.get_elements_data().contains_key("H"));
    }

    #[test]
    fn test_search_by_elements_only() {
        let mut thermo_data = ThermoData::new();

        // Test with elements that should form compounds together
        let elements = vec!["H".to_string(), "O".to_string()];
        let found_substances = thermo_data.search_by_elements_only(elements);

        println!("Found substances with only H and O: {:?}", found_substances);
        println!(
            "Thermo data entries: {}",
            thermo_data.hashmap_of_thermo_data.len()
        );

        // Verify that substances were found and data was populated
        if !found_substances.is_empty() {
            assert!(!thermo_data.hashmap_of_thermo_data.is_empty());
        }
    }
    #[test]
    fn test_search_by_elements_only2() {
        let mut thermo_data = ThermoData::new();

        // Test with elements that should form compounds together
        let elements = vec!["N".to_string(), "O".to_string()];
        let found_substances = thermo_data.search_by_elements_only(elements);

        println!("Found substances with only H and O: {:?}", found_substances);
        println!(
            "Thermo data entries: {}",
            thermo_data.hashmap_of_thermo_data.len()
        );

        // Verify that substances were found and data was populated
        if !found_substances.is_empty() {
            assert!(!thermo_data.hashmap_of_thermo_data.is_empty());
        }
    }
    #[test]
    fn test_search_by_elements_only_with_mock_data() {
        let mut thermo_data = ThermoData::new();

        // Create mock elements data
        let mut mock_elements = HashMap::new();
        mock_elements.insert(
            "H".to_string(),
            vec![
                vec!["H2O".to_string(), "test_lib".to_string()],
                vec!["H2".to_string(), "test_lib".to_string()],
                vec!["CH4".to_string(), "test_lib".to_string()],
            ],
        );
        mock_elements.insert(
            "O".to_string(),
            vec![
                vec!["H2O".to_string(), "test_lib".to_string()],
                vec!["O2".to_string(), "test_lib".to_string()],
            ],
        );
        mock_elements.insert(
            "C".to_string(),
            vec![vec!["CH4".to_string(), "test_lib".to_string()]],
        );

        thermo_data.set_elements_data(mock_elements);

        // Search for substances with only H and O
        let elements = vec!["H".to_string(), "O".to_string()];
        let found_substances = thermo_data.search_by_elements_only(elements);

        // Should find H2O but not CH4 (contains C) or H2/O2 (missing elements)
        println!("Mock test - Found substances: {:?}", found_substances);
        assert!(found_substances.contains(&"H2O".to_string()) || found_substances.is_empty());
    }

    #[test]
    fn test_search_by_elements_only_picks_deterministic_library_for_duplicate_matches() {
        let mut thermo_data = ThermoData::new();

        Arc::make_mut(&mut thermo_data.LibThermoData).insert(
            "beta_lib".to_string(),
            HashMap::from([("H2O".to_string(), Value::Null)]),
        );
        Arc::make_mut(&mut thermo_data.LibThermoData).insert(
            "alpha_lib".to_string(),
            HashMap::from([("H2O".to_string(), Value::Null)]),
        );

        let mut mock_elements = HashMap::new();
        mock_elements.insert(
            "O".to_string(),
            vec![
                vec!["H2O".to_string(), "beta_lib".to_string()],
                vec!["H2O".to_string(), "alpha_lib".to_string()],
            ],
        );
        mock_elements.insert(
            "H".to_string(),
            vec![
                vec!["H2O".to_string(), "beta_lib".to_string()],
                vec!["H2O".to_string(), "alpha_lib".to_string()],
            ],
        );
        thermo_data.set_elements_data(mock_elements);

        let found_substances =
            thermo_data.search_by_elements_only(vec!["H".to_string(), "O".to_string()]);

        assert_eq!(found_substances, vec!["H2O".to_string()]);
        let chosen_library = thermo_data
            .hashmap_of_thermo_data
            .get("H2O")
            .and_then(|map| map.keys().next())
            .cloned();
        assert_eq!(chosen_library, Some("alpha_lib".to_string()));
    }

    #[test]
    fn test_search_by_elements_deduplicates_substances() {
        let mut thermo_data = ThermoData::new();

        let mut mock_elements = HashMap::new();
        mock_elements.insert(
            "H".to_string(),
            vec![
                vec!["H2O".to_string(), "test_lib".to_string()],
                vec!["H2O".to_string(), "other_lib".to_string()],
            ],
        );
        thermo_data.set_elements_data(mock_elements);

        Arc::make_mut(&mut thermo_data.LibThermoData).insert(
            "test_lib".to_string(),
            HashMap::from([("H2O".to_string(), Value::Null)]),
        );
        Arc::make_mut(&mut thermo_data.LibThermoData).insert(
            "other_lib".to_string(),
            HashMap::from([("H2O".to_string(), Value::Null)]),
        );

        let found = thermo_data.search_by_elements(vec!["H".to_string(), "H".to_string()]);
        assert_eq!(found, vec!["H2O".to_string()]);
    }

    #[test]
    fn test_default_repository_handle_is_shared() {
        let repo_a = ThermoData::default_repository().unwrap();
        let repo_b = ThermoData::default_repository().unwrap();

        assert!(Arc::ptr_eq(&repo_a, &repo_b));
    }

    #[test]
    fn test_two_aggregators_share_repository_but_not_mutable_state() {
        let mut left = ThermoData::try_new().unwrap();
        let right = ThermoData::try_new().unwrap();

        assert!(Arc::ptr_eq(&left.repository, &right.repository));
        assert!(Arc::ptr_eq(&left.LibThermoData, &right.LibThermoData));

        left.subs_to_search.push("synthetic_substance".to_string());
        assert!(right.subs_to_search.is_empty());
    }

    #[test]
    fn test_try_new_fresh_uses_a_new_repository_handle() {
        let cached_repository = ThermoData::default_repository().unwrap();
        let fresh = ThermoData::try_new_fresh().unwrap();

        assert!(!Arc::ptr_eq(&cached_repository, &fresh.repository));
    }

    #[test]
    fn test_custom_repository_loads_from_explicit_paths() {
        let mut all_keys = NamedTempFile::new().unwrap();
        let mut substance_base = NamedTempFile::new().unwrap();
        let mut elements = NamedTempFile::new().unwrap();

        writeln!(all_keys, "{}", serde_json::json!([["custom_lib", "H2O"]])).unwrap();
        writeln!(
            substance_base,
            "{}",
            serde_json::json!({
                "custom_lib": {
                    "H2O": {
                        "Cp": [200.0, 1000.0, 6000.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        "composition": "{H: 2, O: 1}",
                        "model": "NASA7"
                    }
                }
            })
        )
        .unwrap();
        writeln!(
            elements,
            "{}",
            serde_json::json!({
                "H": [["H2O", "custom_lib"]],
                "O": [["H2O", "custom_lib"]]
            })
        )
        .unwrap();

        let repository = ThermoData::try_new_from_paths(
            all_keys.path().to_str().unwrap(),
            substance_base.path().to_str().unwrap(),
            elements.path().to_str().unwrap(),
        )
        .unwrap();

        assert!(repository.LibThermoData.contains_key("custom_lib"));
        assert!(repository.ElementsData.contains_key("H"));
        assert_eq!(
            repository.AllLibraries.as_ref(),
            &vec!["custom_lib".to_string()]
        );
    }

    #[test]
    fn test_try_new_from_paths_rejects_empty_catalog() {
        let mut all_keys = NamedTempFile::new().unwrap();
        let mut substance_base = NamedTempFile::new().unwrap();
        let mut elements = NamedTempFile::new().unwrap();

        writeln!(all_keys, "[]").unwrap();
        writeln!(substance_base, "{{}}").unwrap();
        writeln!(elements, "{{}}").unwrap();

        let err = ThermoData::try_new_from_paths(
            all_keys.path().to_str().unwrap(),
            substance_base.path().to_str().unwrap(),
            elements.path().to_str().unwrap(),
        )
        .unwrap_err();

        assert!(matches!(
            err,
            ThermoLibraryError::EmptyLibraryCatalog { .. }
        ));
    }

    #[test]
    fn test_repository_story_loads_from_temporary_fixture_directory() {
        let temp_dir = tempdir().unwrap();
        let all_keys_path = temp_dir.path().join("all_keys.json");
        let substance_base_path = temp_dir.path().join("substance_base.json");
        let elements_path = temp_dir.path().join("elements.json");

        std::fs::write(
            &all_keys_path,
            serde_json::json!([["story_lib", "H2O"]]).to_string(),
        )
        .unwrap();
        std::fs::write(
            &substance_base_path,
            serde_json::json!({
                "story_lib": {
                    "H2O": {
                        "Cp": [200.0, 1000.0, 6000.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        "composition": "{H: 2, O: 1}",
                        "model": "NASA7"
                    }
                }
            })
            .to_string(),
        )
        .unwrap();
        std::fs::write(
            &elements_path,
            serde_json::json!({
                "H": [["H2O", "story_lib"]],
                "O": [["H2O", "story_lib"]]
            })
            .to_string(),
        )
        .unwrap();

        let repository = ThermoData::try_new_from_paths(
            all_keys_path.to_str().unwrap(),
            substance_base_path.to_str().unwrap(),
            elements_path.to_str().unwrap(),
        )
        .unwrap();

        assert!(repository.LibThermoData.contains_key("story_lib"));
        assert!(repository.ElementsData.contains_key("H"));
        assert_eq!(
            repository.AllLibraries.as_ref(),
            &vec!["story_lib".to_string()]
        );
    }

    fn physical_state_fixture_repository() -> ThermoRepository {
        let catalog = HashMap::from([
            (
                "NASA_cond".to_string(),
                HashMap::from([
                    ("Fe(L)".to_string(), serde_json::json!({"record": "liquid"})),
                    (
                        "H2O".to_string(),
                        serde_json::json!({"record": "condensed-default"}),
                    ),
                    ("Fe(a)".to_string(), serde_json::json!({"record": "alpha"})),
                    ("Fe(c)".to_string(), serde_json::json!({"record": "gamma"})),
                ]),
            ),
            (
                "NASA_gas".to_string(),
                HashMap::from([("Fe".to_string(), serde_json::json!({"record": "gas"}))]),
            ),
            (
                "nuig_thermo".to_string(),
                HashMap::from([(
                    "CH2(S)".to_string(),
                    serde_json::json!({"record": "electronic-state-label"}),
                )]),
            ),
        ]);
        ThermoRepository::from_parts(
            Vec::new(),
            catalog,
            HashMap::new(),
            vec!["NASA_cond".into(), "NASA_gas".into(), "nuig_thermo".into()],
            HashMap::new(),
            HashMap::new(),
            vec!["NASA_cond".into(), "NASA_gas".into(), "nuig_thermo".into()],
            Vec::new(),
        )
    }

    #[test]
    fn physical_state_lookup_selects_nasa_condensed_liquid_record() {
        let repository = physical_state_fixture_repository();
        let query = ThermoRecordQuery::new("Fe").with_physical_state(PhysicalState::Liquid);

        let record = repository
            .resolve_thermo_record("NASA_cond", &query)
            .unwrap()
            .unwrap();

        assert_eq!(record.record_key, "Fe(L)");
        assert_eq!(record.physical_state, Some(PhysicalState::Liquid));
        assert_eq!(
            record.state_evidence,
            Some(PhysicalStateEvidence::KeyConvention)
        );
    }

    #[test]
    fn physical_state_lookup_refuses_to_guess_between_solid_polymorphs() {
        let repository = physical_state_fixture_repository();
        let query = ThermoRecordQuery::new("Fe").with_physical_state(PhysicalState::Solid);

        let error = repository
            .resolve_thermo_record("NASA_cond", &query)
            .unwrap_err();

        let ThermoLibraryError::AmbiguousPhysicalState { mut candidates, .. } = error else {
            panic!("solid Fe must require an explicit NASA-condensed polymorph key");
        };
        candidates.sort();
        assert_eq!(candidates, vec!["Fe(a)".to_string(), "Fe(c)".to_string()]);
    }

    #[test]
    fn physical_state_lookup_keeps_nuig_parenthesized_labels_literal() {
        let repository = physical_state_fixture_repository();
        let record = repository
            .resolve_thermo_record("nuig_thermo", &ThermoRecordQuery::new("CH2(S)"))
            .unwrap()
            .unwrap();

        assert_eq!(record.record_key, "CH2(S)");
        assert_eq!(record.physical_state, None);
    }

    #[test]
    fn physical_state_lookup_uses_nasa_gas_default_without_key_suffix() {
        let repository = physical_state_fixture_repository();
        let query = ThermoRecordQuery::new("Fe").with_physical_state(PhysicalState::Gas);
        let record = repository
            .resolve_thermo_record("NASA_gas", &query)
            .unwrap()
            .unwrap();

        assert_eq!(record.record_key, "Fe");
        assert_eq!(record.physical_state, Some(PhysicalState::Gas));
        assert_eq!(
            record.state_evidence,
            Some(PhysicalStateEvidence::LibraryDefault)
        );
    }

    #[test]
    fn physical_state_lookup_accepts_an_untagged_nasa_condensed_record() {
        let repository = physical_state_fixture_repository();
        let query = ThermoRecordQuery::new("H2O").with_physical_state(PhysicalState::Liquid);
        let record = repository
            .resolve_thermo_record("NASA_cond", &query)
            .unwrap()
            .unwrap();

        assert_eq!(record.record_key, "H2O");
        assert_eq!(record.physical_state, Some(PhysicalState::Condensed));
    }
}
