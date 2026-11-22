use RustedSciThe::numerical::Radau::api::RadauError;
use RustedSciThe::numerical::StiffAPI::StiffSolverKind;
use RustedSciThe::numerical::StiffAPI::{
    StiffApiError, StiffCancellationError, StiffPostprocessError, StiffRestartError,
    StiffSolveError, StiffTelemetryError,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum NetworkIvpError {
    #[error("invalid network field '{field}': {reason}")]
    InvalidInput { field: String, reason: &'static str },
    #[error("duplicate species ID '{id}'")]
    DuplicateSpecies { id: String },
    #[error("duplicate reaction ID '{id}'")]
    DuplicateReaction { id: String },
    #[error("reaction '{reaction_id}' references unknown species '{species_id}'")]
    UnknownSpecies {
        reaction_id: String,
        species_id: String,
    },
    #[error("reaction '{reaction_id}' has unsupported symbol '{symbol}' in its rate law")]
    RateSymbol { reaction_id: String, symbol: String },
    #[error("SolidState stage '{reaction_id}' has invalid {field}: {reason}")]
    SolidStateInput {
        reaction_id: String,
        field: &'static str,
        reason: &'static str,
    },
    #[error(
        "SolidState stages '{reaction_id}' and '{other_reaction_id}' both consume '{species_id}'"
    )]
    SharedSolidStateResource {
        reaction_id: String,
        other_reaction_id: String,
        species_id: String,
    },
    #[error(
        "SolidState stage '{reaction_id}' needs {required} mol/m^3 of '{species_id}', but only {available} mol/m^3 is initially available"
    )]
    SolidStateResourceShortfall {
        reaction_id: String,
        species_id: String,
        required: f64,
        available: f64,
    },
    #[error("SolidState stage '{reaction_id}' model {model} is invalid at alpha={alpha}: {reason}")]
    SolidStateDomain {
        reaction_id: String,
        model: String,
        alpha: f64,
        reason: &'static str,
    },
    #[error("{family:?} does not support {capability}")]
    UnsupportedCapability {
        family: StiffSolverKind,
        capability: &'static str,
    },
    #[error("temperature program contains unsupported symbol '{symbol}'")]
    TemperatureSymbol { symbol: String },
    #[error("temperature is invalid at t={time}: {value} K")]
    TemperatureDomain { time: f64, value: f64 },
    #[error("temperature rebind expects {expected} values, got {actual}")]
    TemperatureBindingMismatch {
        expected: &'static str,
        actual: &'static str,
    },
    #[error("rate coefficient binding expects {expected} reactions, got {actual}")]
    RateBindingCountMismatch { expected: usize, actual: usize },
    #[error("rate coefficient binding contains duplicate reaction ID '{reaction_id}'")]
    DuplicateRateBinding { reaction_id: String },
    #[error("reaction '{reaction_id}' has no rebindable rate coefficient")]
    UnknownRateBinding { reaction_id: String },
    #[error(
        "reaction '{reaction_id}' rate binding kind mismatch: expected {expected}, got {actual}"
    )]
    RateBindingKindMismatch {
        reaction_id: String,
        expected: &'static str,
        actual: &'static str,
    },
    #[error("arbitrary symbolic temperature programs cannot be rebound numerically")]
    UnbindableTemperatureProgram,
    #[error(
        "reaction-network RHS for species '{species_id}' (component {component}) is non-finite at t={time}"
    )]
    RateDomain {
        time: f64,
        component: usize,
        species_id: String,
    },
    #[error("RST preparation failed: {0}")]
    Preparation(#[source] StiffApiError),
    #[error("RST integration failed: {0}")]
    Integration(#[source] StiffSolveError),
    #[error("RST restart failed: {0}")]
    Restart(#[source] StiffRestartError),
    #[error("RST cancellation-aware integration failed: {0}")]
    CancellationAwareIntegration(#[source] StiffCancellationError),
    #[error("RST result extraction failed: {0}")]
    Postprocessing(#[source] StiffPostprocessError),
    #[error("RST telemetry snapshot failed: {0}")]
    Telemetry(#[source] StiffTelemetryError),
    #[error("Radau result extraction failed: {0}")]
    RadauExtraction(#[source] RadauError),
    #[error("solver output is invalid: {reason}")]
    InvalidOutput { reason: &'static str },
    #[error(
        "solver output concentration for species '{species_id}' is invalid at t={time}: {value} mol/m^3"
    )]
    InvalidConcentration {
        species_id: String,
        time: f64,
        value: f64,
    },
    #[error(
        "solver output conversion for stage '{reaction_id}' is invalid at t={time}: alpha={value}"
    )]
    InvalidConversion {
        reaction_id: String,
        time: f64,
        value: f64,
    },
    #[error("solver output time sample {index} is invalid: t={time}, previous={previous_time:?}")]
    InvalidTimeSample {
        index: usize,
        time: f64,
        previous_time: Option<f64>,
    },
}

pub(crate) fn invalid(field: impl Into<String>, reason: &'static str) -> NetworkIvpError {
    NetworkIvpError::InvalidInput {
        field: field.into(),
        reason,
    }
}
