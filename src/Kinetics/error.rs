use thiserror::Error;

pub type KineticsResult<T> = Result<T, KineticsError>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum KinDataStateRequirement {
    #[error("shortcut metadata or shortcut reaction IDs")]
    ShortcutMetadataOrIds,
    #[error("library pairs or library reaction IDs")]
    LibraryPairsOrIds,
    #[error("raw reaction values or document IDs")]
    RawValuesOrDocumentIds,
    #[error("parsed reaction data")]
    ParsedReactionData,
    #[error("reaction data from mechanism construction")]
    MechanismReactionData,
    #[error("parsed data, equations, or analyzed artifacts")]
    AnalyzedDataOrArtifacts,
    #[error("a coherent normalized reaction view")]
    NormalizedReactionView,
}

#[derive(Debug, Error)]
pub enum KineticsError {
    #[error("{context}: {source}")]
    Context {
        context: String,
        #[source]
        source: Box<KineticsError>,
    },
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("library `{0}` was not found")]
    MissingLibrary(String),
    #[error("library `{library}` was not found in catalog `{catalog}`")]
    MissingCatalogLibrary {
        catalog: &'static str,
        library: String,
    },
    #[error("catalog file `{path}` changed repeatedly while it was being loaded")]
    CatalogChangedDuringLoad { path: String },
    #[error("the `{cache}` kinetics cache lock is poisoned")]
    CacheLockPoisoned { cache: &'static str },
    #[error("reaction `{0}` was not found")]
    MissingReaction(String),
    #[error("equation `{equation}` maps to multiple reactions: {reaction_ids:?}")]
    AmbiguousReactionEquation {
        equation: String,
        reaction_ids: Vec<String>,
    },
    #[error("invalid shortcut range `{0}`")]
    InvalidShortcutRange(String),
    #[error("invalid reaction shortcut `{input}`: {reason}")]
    InvalidReactionShortcut { input: String, reason: String },
    #[error("mechanism `{mechanism}` contains duplicate reaction ID `{reaction_id}`")]
    DuplicateReactionId {
        mechanism: String,
        reaction_id: String,
    },
    #[error("mechanism `{mechanism}` has non-numeric reaction ID `{reaction_id}`")]
    InvalidMechanismReactionId {
        mechanism: String,
        reaction_id: String,
    },
    #[error("reaction data validation failed: {0}")]
    InvalidReactionData(String),
    #[error("reaction `{reaction}` has invalid payload field `{field}`: {reason}")]
    InvalidReactionPayload {
        reaction: String,
        field: String,
        reason: String,
    },
    #[error("reaction `{reaction}` received invalid `{field}` input: {reason}")]
    InvalidKineticInput {
        reaction: String,
        field: &'static str,
        reason: String,
    },
    #[error("reaction `{reaction}` requires `{field}` input")]
    MissingKineticInput {
        reaction: String,
        field: &'static str,
    },
    #[error("invalid reaction syntax in `{input}`: {reason}")]
    ReactionSyntax { input: String, reason: String },
    #[error("invalid formula `{formula}`: {reason}")]
    FormulaSyntax { formula: String, reason: String },
    #[error("formula `{formula}` contains unknown element `{symbol}`")]
    UnknownElement { formula: String, symbol: String },
    #[error("custom group `{group}` is invalid: {reason}")]
    InvalidCustomGroup { group: String, reason: String },
    #[error("formula {resource} limit exceeded: requested {requested}, limit {limit}")]
    FormulaLimitExceeded {
        resource: &'static str,
        requested: usize,
        limit: usize,
    },
    #[error(
        "reaction `{reaction}` references species `{species}` missing from the supplied species list"
    )]
    UnknownReactionSpecies { reaction: String, species: String },
    #[error("failed to compile {context} regex: {source}")]
    RegexCompile {
        context: &'static str,
        #[source]
        source: regex::Error,
    },
    #[error("unsupported reaction model: {0}")]
    UnsupportedReactionModel(String),
    #[error("rate evaluation for reaction `{reaction}` produced a non-finite value: {value}")]
    NonFiniteRate { reaction: String, value: String },
    #[error("{resource} request has {requested} items; configured limit is {limit}")]
    ResourceLimitExceeded {
        resource: &'static str,
        requested: usize,
        limit: usize,
    },
    #[error("length mismatch in {context}: left={left}, right={right}")]
    LengthMismatch {
        context: &'static str,
        left: usize,
        right: usize,
    },
    #[error("invalid KinData state: {0}")]
    InvalidState(String),
    #[error("operation `{operation}` is not allowed in KinData state `{state}`")]
    OperationNotAllowedInState {
        operation: &'static str,
        state: &'static str,
    },
    #[error("KinData state `{state}` does not satisfy its requirement: {requirement}")]
    StateContractMismatch {
        state: &'static str,
        requirement: KinDataStateRequirement,
    },
    #[error("solver error: {0}")]
    Solver(String),
}

impl KineticsError {
    /// Add operation context without discarding the typed cause.
    pub fn with_context(self, context: impl Into<String>) -> Self {
        Self::Context {
            context: context.into(),
            source: Box::new(self),
        }
    }
}
