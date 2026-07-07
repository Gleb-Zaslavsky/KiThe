use thiserror::Error;

pub type KineticsResult<T> = Result<T, KineticsError>;

#[derive(Debug, Error)]
pub enum KineticsError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("library `{0}` was not found")]
    MissingLibrary(String),
    #[error("reaction `{0}` was not found")]
    MissingReaction(String),
    #[error("invalid shortcut range `{0}`")]
    InvalidShortcutRange(String),
    #[error("reaction data validation failed: {0}")]
    InvalidReactionData(String),
    #[error("length mismatch in {context}: left={left}, right={right}")]
    LengthMismatch {
        context: &'static str,
        left: usize,
        right: usize,
    },
    #[error("invalid KinData state: {0}")]
    InvalidState(String),
    #[error("solver error: {0}")]
    Solver(String),
}
