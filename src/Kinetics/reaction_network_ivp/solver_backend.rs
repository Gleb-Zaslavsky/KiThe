use RustedSciThe::numerical::BDF::BDF_api::BdfSolverOptions;
use RustedSciThe::numerical::BE::BeSolverOptions;
use RustedSciThe::numerical::LSODE2::Lsode2ProblemConfig;
use RustedSciThe::numerical::Radau::api::{RadauConfig, RadauOutputPolicy, RadauProblem};
use RustedSciThe::numerical::SolverTelemetryMode;
use RustedSciThe::numerical::StiffAPI::{
    PreparedStiffSolver, StiffCallbackPolicy, StiffFrontend, StiffMatrixLayout, StiffPreset,
    StiffSolverKind, StiffSolverRequest, StiffTelemetrySnapshot,
};
use RustedSciThe::numerical::logging::SolverLogLevel;
use nalgebra::DVector;

use super::assembly::AssembledNetwork;
use super::error::{NetworkIvpError, invalid};
use super::problem::NetworkProblem;
use super::rate_law::RateCoefficientBinding;

/// Conservative vector-IVP controls; unsupported native capability combinations
/// are rejected by StiffAPI during preparation, never silently substituted.
#[derive(Debug, Clone)]
pub struct NetworkIvpOptions {
    pub family: StiffSolverKind,
    pub frontend: StiffFrontend,
    pub matrix_layout: StiffMatrixLayout,
    /// None selects sequential callbacks where the native family exposes them.
    pub callback_policy: Option<StiffCallbackPolicy>,
    pub rtol: f64,
    pub atol: f64,
    pub max_step: f64,
    pub telemetry: SolverTelemetryMode,
    pub logging: SolverLogLevel,
}

impl Default for NetworkIvpOptions {
    fn default() -> Self {
        Self {
            family: StiffSolverKind::Bdf,
            frontend: StiffFrontend::ExprLegacy,
            matrix_layout: StiffMatrixLayout::Dense,
            callback_policy: None,
            rtol: 1e-8,
            atol: 1e-10,
            max_step: f64::INFINITY,
            telemetry: SolverTelemetryMode::Off,
            logging: SolverLogLevel::Off,
        }
    }
}

impl NetworkIvpOptions {
    pub(crate) fn validate(&self) -> Result<(), NetworkIvpError> {
        if !self.rtol.is_finite() || self.rtol <= 0.0 {
            return Err(invalid("options.rtol", "must be finite and positive"));
        }
        if !self.atol.is_finite() || self.atol <= 0.0 {
            return Err(invalid("options.atol", "must be finite and positive"));
        }
        if self.max_step.is_nan() || self.max_step <= 0.0 {
            return Err(invalid("options.max_step", "must be positive"));
        }
        if self.family == StiffSolverKind::BackwardEuler
            && (self.rtol != 1e-8 || self.atol != 1e-10)
        {
            return Err(NetworkIvpError::UnsupportedCapability {
                family: self.family,
                capability: "rtol/atol integration controls; Backward Euler uses a separate Newton tolerance",
            });
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct NetworkIvpSolution {
    /// Exact inputs accepted for this run; retained so stored trajectories can
    /// be traced back to their prepared-run bindings.
    pub run_metadata: NetworkIvpRunMetadata,
    pub species_ids: Vec<String>,
    pub solid_state_reaction_ids: Vec<String>,
    pub time_seconds: Vec<f64>,
    /// Row-major physical concentrations in mol/m^3, in `species_ids` order.
    pub concentrations: Vec<Vec<f64>>,
    /// Row-major conversion extents in `solid_state_reaction_ids` order.
    pub conversions: Vec<Vec<f64>>,
    pub temperature_kelvin: Vec<f64>,
    /// None means telemetry was not requested, not that the backend did no work.
    pub telemetry: Option<StiffTelemetrySnapshot>,
    /// KiThe-owned cumulative stage data for the prepared session up to this
    /// accepted result. Native RST diagnostics remain in `telemetry`.
    pub stage_diagnostics: Option<super::diagnostics::NetworkIvpStageDiagnostics>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NetworkIvpRunMetadata {
    pub initial_concentrations: Vec<f64>,
    pub time_interval_seconds: (f64, f64),
    pub temperature_binding: Option<super::temperature::TemperatureBinding>,
    pub rate_coefficients: Vec<RateCoefficientBinding>,
}

/// Bounded raw native samples observed after a failed integration attempt.
/// These values are diagnostic only, lack an RST solve ID, and have not passed
/// physical acceptance.
#[derive(Debug, Clone)]
pub struct NetworkIvpPartialTrajectory {
    pub time_seconds: Vec<f64>,
    pub states: Vec<Vec<f64>>,
    pub omitted_samples: usize,
}

pub(crate) fn prepare_native(
    problem: &NetworkProblem,
    assembled: AssembledNetwork,
    options: &NetworkIvpOptions,
) -> Result<PreparedStiffSolver, NetworkIvpError> {
    let (t0, tf) = problem.time_interval();
    let initial = problem.initial_state(&problem.initial_concentrations());
    let request = match options.family {
        StiffSolverKind::Bdf => {
            let native = BdfSolverOptions::for_bdf(
                assembled.equations,
                assembled.variables,
                "t".to_owned(),
                t0,
                DVector::from_vec(initial.clone()),
                tf,
                options.max_step,
                options.rtol,
                options.atol,
                None,
                false,
                None,
            );
            let native = if assembled.parameter_names.is_empty() {
                native
            } else {
                native
                    .with_equation_parameters(assembled.parameter_names)
                    .with_equation_parameter_values(DVector::from_vec(assembled.parameter_values))
            };
            StiffSolverRequest::Bdf { options: native }
        }
        StiffSolverKind::Radau => {
            let problem = RadauProblem::new(assembled.equations, assembled.variables, "t")
                .with_parameter_names(assembled.parameter_names);
            let config = RadauConfig {
                t0,
                t_bound: tf,
                rtol: options.rtol,
                atol: options.atol,
                max_step: options.max_step,
                output: RadauOutputPolicy::Dense,
                ..RadauConfig::default()
            };
            StiffSolverRequest::Radau {
                problem,
                config,
                initial_state: initial.clone(),
                parameters: assembled.parameter_values,
            }
        }
        StiffSolverKind::BackwardEuler => {
            let native = BeSolverOptions::new(
                assembled.equations,
                assembled.variables,
                "t".to_owned(),
                1e-10,
                30,
                options.max_step.is_finite().then_some(options.max_step),
                t0,
                tf,
                DVector::from_vec(initial.clone()),
            );
            StiffSolverRequest::BackwardEuler {
                options: native,
                equation_parameters: assembled.parameter_names,
                parameter_values: assembled.parameter_values,
            }
        }
        StiffSolverKind::Lsode2 => {
            let native = Lsode2ProblemConfig::new(
                assembled.equations,
                assembled.variables,
                "t".to_owned(),
                t0,
                DVector::from_vec(initial.clone()),
                tf,
                options.max_step,
                options.rtol,
                options.atol,
            );
            let native = if assembled.parameter_names.is_empty() {
                native
            } else {
                native
                    .with_equation_parameters(assembled.parameter_names)
                    .with_equation_parameter_values(DVector::from_vec(assembled.parameter_values))
            };
            StiffSolverRequest::Lsode2 {
                config: native,
                logging_level: None,
            }
        }
    };
    let mut preset = StiffPreset::default()
        .with_frontend(options.frontend)
        .with_matrix_layout(options.matrix_layout)
        .with_telemetry(options.telemetry)
        .with_logging(options.logging);
    if let Some(policy) = options.callback_policy {
        preset = preset.with_callback_policy(policy);
    } else if options.family != StiffSolverKind::BackwardEuler {
        preset = preset.with_callback_policy(StiffCallbackPolicy::Sequential);
    }
    request
        .prepare(Some(&preset))
        .map_err(NetworkIvpError::Preparation)
}
