//! Solver-backend facade for condensed-phase IVP tasks.
//!
//! The IVP reactor workflow starts from symbolic equations and typed initial
//! conditions, so the user-facing facade selects a solver family and keeps
//! each family's native controls behind one typed adapter.

use RustedSciThe::command_interpreter::task_parser_ivp_hybrid::HybridStiffSettings;
use RustedSciThe::numerical::BDF::BDF_api::{BdfSolverOptions, BdfTelemetryMode};
use RustedSciThe::numerical::LSODE2::{
    Lsode2AotProfile, Lsode2AotToolchain, Lsode2ControllerConfig, Lsode2LinearSolverPolicy,
    Lsode2LinearSystemStructure, Lsode2Method, Lsode2NativeExecutionConfig, Lsode2ProblemConfig,
    Lsode2ResidualJacobianSource, Lsode2SymbolicAssemblyBackend, Lsode2SymbolicExecutionMode,
};
use RustedSciThe::numerical::Radau::api::{
    RadauConfig, RadauFrontend, RadauMatrixLayout, RadauOutputPolicy, RadauProblem,
    RadauTelemetryMode,
};
use RustedSciThe::numerical::SolverTelemetryMode;
use RustedSciThe::numerical::StiffAPI::{
    PreparedStiffSolver, StiffApiError, StiffCallbackPolicy, StiffExecution, StiffFrontend,
    StiffMatrixLayout, StiffPreset, StiffSolverKind, StiffSolverRequest,
};
use RustedSciThe::numerical::logging::SolverLogLevel;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DVector;
use std::fmt;

const REACTOR_IVP_SPECIES_STATE_OFFSET: usize = 2;

/// Typed failure while mapping and preparing a condensed IVP solver request.
#[derive(Debug)]
pub enum ReactorIvpBackendError {
    InvalidConfiguration(String),
    UnsupportedCapability {
        family: ReactorIvpSolverFamily,
        capability: &'static str,
    },
    StiffApiPreparation(StiffApiError),
}

impl fmt::Display for ReactorIvpBackendError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfiguration(msg) => write!(f, "Invalid configuration: {msg}"),
            Self::UnsupportedCapability { family, capability } => {
                write!(f, "{family:?} does not support {capability}")
            }
            Self::StiffApiPreparation(error) => {
                write!(f, "StiffAPI request preparation failed: {error}")
            }
        }
    }
}

impl std::error::Error for ReactorIvpBackendError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::StiffApiPreparation(error) => Some(error),
            Self::InvalidConfiguration(_) | Self::UnsupportedCapability { .. } => None,
        }
    }
}

/// Time-integration family exposed by the condensed-phase IVP facade.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactorIvpMethod {
    Auto,
    Bdf,
    Adams,
}

/// Integration family. LSODE2's internal Auto/BDF/Adams controller remains a
/// separate setting and is meaningful only when `family == Lsode2`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ReactorIvpSolverFamily {
    #[default]
    Lsode2,
    Bdf,
    Radau,
}

impl ReactorIvpSolverFamily {
    pub const fn to_stiff_kind(self) -> StiffSolverKind {
        match self {
            Self::Lsode2 => StiffSolverKind::Lsode2,
            Self::Bdf => StiffSolverKind::Bdf,
            Self::Radau => StiffSolverKind::Radau,
        }
    }
}

impl Default for ReactorIvpMethod {
    fn default() -> Self {
        Self::Auto
    }
}

impl ReactorIvpMethod {
    pub fn to_controller(self) -> Lsode2ControllerConfig {
        match self {
            Self::Auto => Lsode2ControllerConfig::automatic_adams_bdf(),
            Self::Bdf => Lsode2ControllerConfig::bdf_only(),
            Self::Adams => Lsode2ControllerConfig::adams_only(),
        }
    }

    /// LSODE2 currently exposes the fixed BDF method at the problem layer.
    ///
    /// The Auto/BDF/Adams choice is carried by the controller instead of by a
    /// separate lower-level method enum.
    pub fn to_lsode2_method(self) -> Lsode2Method {
        Lsode2Method::Bdf
    }
}

/// Symbolic assembly route used before Lambdify/AOT lowering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactorIvpSymbolicBackend {
    AtomView,
    ExprLegacy,
}

impl Default for ReactorIvpSymbolicBackend {
    fn default() -> Self {
        Self::AtomView
    }
}

impl ReactorIvpSymbolicBackend {
    pub fn to_rusted(self) -> Lsode2SymbolicAssemblyBackend {
        match self {
            Self::AtomView => Lsode2SymbolicAssemblyBackend::AtomView,
            Self::ExprLegacy => Lsode2SymbolicAssemblyBackend::ExprLegacy,
        }
    }
}

/// Matrix representation used by generated IVP callbacks.
///
/// The public facade intentionally keeps `Banded` parameter-free. LSODE2 can
/// infer the effective band width from the symbolic Jacobian route, so the
/// user-facing contract should not ask for `kl`/`ku` manually.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactorIvpMatrixBackend {
    Dense,
    Sparse,
    Banded,
}

impl Default for ReactorIvpMatrixBackend {
    fn default() -> Self {
        Self::Sparse
    }
}

impl ReactorIvpMatrixBackend {
    /// Validate the matrix route for a problem with `variable_count` states.
    pub fn validate_for_dimension(
        self,
        variable_count: usize,
    ) -> Result<(), ReactorIvpBackendError> {
        match self {
            Self::Dense => Ok(()),
            Self::Sparse => Ok(()),
            Self::Banded => {
                if variable_count == 0 {
                    return Err(ReactorIvpBackendError::InvalidConfiguration(
                        "banded matrix backend requires at least one variable".to_string(),
                    ));
                }
                Ok(())
            }
        }
    }

    pub fn to_rusted(self) -> Lsode2LinearSystemStructure {
        match self {
            Self::Dense => Lsode2LinearSystemStructure::Dense,
            Self::Sparse => Lsode2LinearSystemStructure::Sparse,
            Self::Banded => Lsode2LinearSystemStructure::BandedAuto,
        }
    }

    pub fn is_banded(self) -> bool {
        matches!(self, Self::Banded)
    }
}

/// Execution route for symbolic residual/Jacobian lowering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorIvpExecutionBackend {
    Lambdify,
    Aot {
        toolchain: Lsode2AotToolchain,
        profile: Lsode2AotProfile,
    },
}

/// Scheduling policy for residual/Jacobian callbacks on supported Lambdify routes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorIvpCallbackPolicy {
    Sequential,
    Parallel { min_work: usize },
    Auto { min_work: usize },
}

impl ReactorIvpCallbackPolicy {
    const fn to_stiff(self) -> StiffCallbackPolicy {
        match self {
            Self::Sequential => StiffCallbackPolicy::Sequential,
            Self::Parallel { min_work } => StiffCallbackPolicy::Parallel { min_work },
            Self::Auto { min_work } => StiffCallbackPolicy::Auto { min_work },
        }
    }
}

impl Default for ReactorIvpExecutionBackend {
    fn default() -> Self {
        Self::Lambdify
    }
}

/// Family-owned iteration controls. A request can no longer carry a Radau
/// Newton cap while selecting LSODE2, or LSODE2 native execution on BDF.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ReactorIvpFamilyOptions {
    Lsode2 {
        controller: ReactorIvpMethod,
        native_execution: Lsode2NativeExecutionConfig,
    },
    Bdf {
        max_steps: usize,
    },
    Radau {
        max_steps: usize,
        max_newton_iterations: usize,
    },
}

impl Default for ReactorIvpFamilyOptions {
    fn default() -> Self {
        Self::Lsode2 {
            controller: ReactorIvpMethod::default(),
            native_execution: Lsode2NativeExecutionConfig::faithful_bdf_solve(200_000, 200_000),
        }
    }
}

impl ReactorIvpFamilyOptions {
    pub const fn family(self) -> ReactorIvpSolverFamily {
        match self {
            Self::Lsode2 { .. } => ReactorIvpSolverFamily::Lsode2,
            Self::Bdf { .. } => ReactorIvpSolverFamily::Bdf,
            Self::Radau { .. } => ReactorIvpSolverFamily::Radau,
        }
    }

    pub const fn lsode2_controller(self) -> Option<ReactorIvpMethod> {
        match self {
            Self::Lsode2 { controller, .. } => Some(controller),
            Self::Bdf { .. } | Self::Radau { .. } => None,
        }
    }
}

/// Early-termination rule for the condensed IVP solver.
///
/// The reactor UI exposes this as a selected species concentration plus a
/// lower-or-equal threshold. We keep the backend representation index-based so
/// the typed config stays `Copy` and can be cheaply moved through the task.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactorIvpStopCondition {
    pub species_index: usize,
    pub threshold: f64,
}

impl ReactorIvpStopCondition {
    pub fn new(species_index: usize, threshold: f64) -> Self {
        Self {
            species_index,
            threshold,
        }
    }

    /// Validate the stop rule against the current variable list.
    pub fn validate(self, species_count: usize) -> Result<(), ReactorIvpBackendError> {
        if self.species_index >= species_count {
            return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                "stop condition species index {} is out of range for {} species",
                self.species_index, species_count
            )));
        }
        if !self.threshold.is_finite() {
            return Err(ReactorIvpBackendError::InvalidConfiguration(
                "stop condition threshold must be finite".to_string(),
            ));
        }
        if self.threshold < 0.0 {
            return Err(ReactorIvpBackendError::InvalidConfiguration(
                "stop condition threshold must not be negative".to_string(),
            ));
        }
        Ok(())
    }
}

/// Typed condensed-phase IVP solver configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactorIvpSolverConfig {
    pub family_options: ReactorIvpFamilyOptions,
    pub symbolic_backend: ReactorIvpSymbolicBackend,
    pub matrix_backend: ReactorIvpMatrixBackend,
    pub execution_backend: ReactorIvpExecutionBackend,
    /// Optional override; `None` retains the selected StiffAPI route default.
    pub callback_policy: Option<ReactorIvpCallbackPolicy>,
    pub stop_condition: Option<ReactorIvpStopCondition>,
    /// Lower dimensionless spatial coordinate `z=x/L` (legacy field name).
    pub x0: f64,
    /// Upper dimensionless spatial coordinate `z=x/L` (legacy field name).
    pub x_bound: f64,
    pub first_step: Option<f64>,
    pub max_step: f64,
    pub rtol: f64,
    pub atol: f64,
    pub telemetry: SolverTelemetryMode,
    pub logging: SolverLogLevel,
}

impl Default for ReactorIvpSolverConfig {
    fn default() -> Self {
        Self {
            family_options: ReactorIvpFamilyOptions::default(),
            symbolic_backend: ReactorIvpSymbolicBackend::default(),
            matrix_backend: ReactorIvpMatrixBackend::default(),
            execution_backend: ReactorIvpExecutionBackend::default(),
            callback_policy: None,
            stop_condition: None,
            x0: 0.0,
            x_bound: 1.0,
            first_step: None,
            max_step: 1.0,
            rtol: 1e-6,
            atol: 1e-8,
            telemetry: SolverTelemetryMode::Off,
            logging: SolverLogLevel::Off,
        }
    }
}

impl ReactorIvpSolverConfig {
    pub fn default_lambdify() -> Self {
        Self::default()
    }

    pub fn with_lsode2_controller(mut self, controller: ReactorIvpMethod) -> Self {
        let native_execution = match self.family_options {
            ReactorIvpFamilyOptions::Lsode2 {
                native_execution, ..
            } => native_execution,
            _ => match ReactorIvpFamilyOptions::default() {
                ReactorIvpFamilyOptions::Lsode2 {
                    native_execution, ..
                } => native_execution,
                _ => unreachable!(),
            },
        };
        self.family_options = ReactorIvpFamilyOptions::Lsode2 {
            controller,
            native_execution,
        };
        self
    }

    pub fn with_family(mut self, family: ReactorIvpSolverFamily) -> Self {
        self.family_options = match family {
            ReactorIvpSolverFamily::Lsode2 => ReactorIvpFamilyOptions::default(),
            ReactorIvpSolverFamily::Bdf => ReactorIvpFamilyOptions::Bdf { max_steps: 200_000 },
            ReactorIvpSolverFamily::Radau => ReactorIvpFamilyOptions::Radau {
                max_steps: 200_000,
                max_newton_iterations: 6,
            },
        };
        if family != ReactorIvpSolverFamily::Lsode2 {
            self.matrix_backend = ReactorIvpMatrixBackend::Dense;
            self.execution_backend = ReactorIvpExecutionBackend::Lambdify;
            self.stop_condition = None;
        }
        self
    }

    pub const fn family(self) -> ReactorIvpSolverFamily {
        self.family_options.family()
    }

    pub const fn lsode2_controller(self) -> Option<ReactorIvpMethod> {
        self.family_options.lsode2_controller()
    }

    pub fn with_telemetry(mut self, telemetry: SolverTelemetryMode) -> Self {
        self.telemetry = telemetry;
        self
    }

    pub fn with_logging(mut self, logging: SolverLogLevel) -> Self {
        self.logging = logging;
        self
    }

    pub fn with_symbolic_backend(mut self, backend: ReactorIvpSymbolicBackend) -> Self {
        self.symbolic_backend = backend;
        self
    }

    pub fn with_matrix_backend(mut self, backend: ReactorIvpMatrixBackend) -> Self {
        self.matrix_backend = backend;
        self
    }

    pub fn with_execution_backend(mut self, backend: ReactorIvpExecutionBackend) -> Self {
        self.execution_backend = backend;
        self
    }

    /// Override residual/Jacobian callback scheduling on supported routes.
    pub fn with_callback_policy(mut self, policy: ReactorIvpCallbackPolicy) -> Self {
        self.callback_policy = Some(policy);
        self
    }

    /// Attach a stop condition that ends the solve when the selected species
    /// concentration reaches the supplied lower bound.
    pub fn with_stop_condition_le(mut self, species_index: usize, threshold: f64) -> Self {
        self.stop_condition = Some(ReactorIvpStopCondition::new(species_index, threshold));
        self
    }

    /// Clear any configured stop condition.
    pub fn without_stop_condition(mut self) -> Self {
        self.stop_condition = None;
        self
    }

    /// Convenience constructor for the explicit AOT route.
    pub fn with_aot_backend(
        mut self,
        toolchain: Lsode2AotToolchain,
        profile: Lsode2AotProfile,
    ) -> Self {
        self.execution_backend = ReactorIvpExecutionBackend::Aot { toolchain, profile };
        self
    }

    pub fn with_integration_domain(mut self, x0: f64, x_bound: f64) -> Self {
        self.x0 = x0;
        self.x_bound = x_bound;
        self
    }

    pub fn with_first_step(mut self, first_step: Option<f64>) -> Self {
        self.first_step = first_step;
        self
    }

    pub fn with_max_step(mut self, max_step: f64) -> Self {
        self.max_step = max_step;
        self
    }

    pub fn with_rtol(mut self, rtol: f64) -> Self {
        self.rtol = rtol;
        self
    }

    pub fn with_atol(mut self, atol: f64) -> Self {
        self.atol = atol;
        self
    }

    pub fn with_bdf_max_steps(mut self, max_steps: usize) -> Self {
        if self.family() != ReactorIvpSolverFamily::Bdf {
            self = self.with_family(ReactorIvpSolverFamily::Bdf);
        }
        self.family_options = ReactorIvpFamilyOptions::Bdf { max_steps };
        self
    }

    pub fn with_radau_iteration_limits(
        mut self,
        max_steps: usize,
        max_newton_iterations: usize,
    ) -> Self {
        if self.family() != ReactorIvpSolverFamily::Radau {
            self = self.with_family(ReactorIvpSolverFamily::Radau);
        }
        self.family_options = ReactorIvpFamilyOptions::Radau {
            max_steps,
            max_newton_iterations,
        };
        self
    }

    pub fn with_native_execution(mut self, native_execution: Lsode2NativeExecutionConfig) -> Self {
        let controller = self
            .lsode2_controller()
            .unwrap_or_else(ReactorIvpMethod::default);
        self.family_options = ReactorIvpFamilyOptions::Lsode2 {
            controller,
            native_execution,
        };
        self
    }

    pub fn is_aot(self) -> bool {
        matches!(
            self.execution_backend,
            ReactorIvpExecutionBackend::Aot { .. }
        )
    }

    /// Validate the route-level backend contract before building LSODE2 config.
    pub fn validate_backend_contract(
        self,
        variable_count: usize,
    ) -> Result<(), ReactorIvpBackendError> {
        self.matrix_backend.validate_for_dimension(variable_count)?;
        Ok(())
    }

    pub fn to_rusted_linear_policy(self) -> Lsode2LinearSolverPolicy {
        Lsode2LinearSolverPolicy::Auto
    }

    pub fn to_rusted_residual_source(self) -> Lsode2ResidualJacobianSource {
        let execution = match self.execution_backend {
            ReactorIvpExecutionBackend::Lambdify => Lsode2SymbolicExecutionMode::LambdifyExpr,
            ReactorIvpExecutionBackend::Aot { toolchain, profile } => {
                Lsode2SymbolicExecutionMode::Aot { toolchain, profile }
            }
        };
        Lsode2ResidualJacobianSource::Symbolic {
            assembly: self.symbolic_backend.to_rusted(),
            execution,
        }
    }

    pub fn to_rusted_problem_config(
        self,
        eq_system: Vec<Expr>,
        values: Vec<String>,
        arg: String,
        y0: DVector<f64>,
    ) -> Result<Lsode2ProblemConfig, ReactorIvpBackendError> {
        let (controller, native_execution) = match self.family_options {
            ReactorIvpFamilyOptions::Lsode2 {
                controller,
                native_execution,
            } => (controller, native_execution),
            _ => {
                return Err(ReactorIvpBackendError::InvalidConfiguration(
                    "LSODE2 problem config requires LSODE2 family options".to_string(),
                ));
            }
        };
        if eq_system.len() != values.len() {
            return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                "equation count {} must match variable count {}",
                eq_system.len(),
                values.len()
            )));
        }
        if y0.len() != values.len() {
            return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                "initial-state length {} must match variable count {}",
                y0.len(),
                values.len()
            )));
        }
        self.validate_backend_contract(values.len())?;
        if !self.x0.is_finite() || !self.x_bound.is_finite() {
            return Err(ReactorIvpBackendError::InvalidConfiguration(
                "integration bounds must be finite".to_string(),
            ));
        }
        if self.x_bound <= self.x0 {
            return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                "integration bound must be larger than x0 (x0={}, x_bound={})",
                self.x0, self.x_bound
            )));
        }
        if !self.max_step.is_finite() || self.max_step <= 0.0 {
            return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                "max_step must be finite and positive, got {}",
                self.max_step
            )));
        }
        if !self.rtol.is_finite() || self.rtol <= 0.0 {
            return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                "rtol must be finite and positive, got {}",
                self.rtol
            )));
        }
        if !self.atol.is_finite() || self.atol <= 0.0 {
            return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                "atol must be finite and positive, got {}",
                self.atol
            )));
        }
        if let Some(first_step) = self.first_step {
            if !first_step.is_finite() || first_step <= 0.0 {
                return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                    "first_step must be finite and positive when set, got {}",
                    first_step
                )));
            }
        }
        let species_count = values
            .len()
            .checked_sub(REACTOR_IVP_SPECIES_STATE_OFFSET)
            .ok_or_else(|| {
                ReactorIvpBackendError::InvalidConfiguration(
                    "condensed IVP stop conditions require at least two thermal states".to_string(),
                )
            })?;
        if let Some(stop_condition) = self.stop_condition {
            stop_condition.validate(species_count)?;
        }

        let mut config = Lsode2ProblemConfig::new(
            eq_system,
            values,
            arg,
            self.x0,
            y0,
            self.x_bound,
            self.max_step,
            self.rtol,
            self.atol,
        );
        config.method = controller.to_lsode2_method();
        config.first_step = self.first_step;
        config.controller = controller.to_controller();
        config.residual_jacobian_source = self.to_rusted_residual_source();
        config.linear_system_structure = self.matrix_backend.to_rusted();
        config.linear_solver_policy = self.to_rusted_linear_policy();
        config.native_execution = native_execution;
        if let Some(stop_condition) = self.stop_condition {
            let species_name = config
                .values
                .get(REACTOR_IVP_SPECIES_STATE_OFFSET + stop_condition.species_index)
                .cloned()
                .ok_or_else(|| {
                    ReactorIvpBackendError::InvalidConfiguration(format!(
                        "stop condition species index {} is out of range",
                        stop_condition.species_index
                    ))
                })?;
            config = config.with_stop_condition_le(species_name, stop_condition.threshold);
        }
        Ok(config)
    }

    /// Build a native StiffAPI request while preserving family-specific
    /// controls and rejecting options that cannot be honored by that family.
    pub fn to_stiff_request(
        self,
        eq_system: Vec<Expr>,
        values: Vec<String>,
        arg: String,
        y0: DVector<f64>,
    ) -> Result<StiffSolverRequest, ReactorIvpBackendError> {
        self.validate_common_inputs(&eq_system, &values, &y0)?;
        match self.family() {
            ReactorIvpSolverFamily::Lsode2 => {
                let config = self.to_rusted_problem_config(eq_system, values, arg, y0)?;
                Ok(StiffSolverRequest::Lsode2 {
                    config,
                    logging_level: Some(self.logging),
                })
            }
            ReactorIvpSolverFamily::Bdf => {
                if !matches!(self.execution_backend, ReactorIvpExecutionBackend::Lambdify) {
                    return Err(ReactorIvpBackendError::UnsupportedCapability {
                        family: self.family(),
                        capability: "AOT execution mapping for standalone BDF",
                    });
                }
                if self.matrix_backend != ReactorIvpMatrixBackend::Dense {
                    return Err(ReactorIvpBackendError::UnsupportedCapability {
                        family: self.family(),
                        capability: "non-Dense linear algebra",
                    });
                }
                let max_steps = match self.family_options {
                    ReactorIvpFamilyOptions::Bdf { max_steps } => max_steps,
                    _ => unreachable!("family() and family options are derived together"),
                };
                Ok(StiffSolverRequest::Bdf {
                    options: BdfSolverOptions::for_bdf(
                        eq_system,
                        values,
                        arg,
                        self.x0,
                        y0,
                        self.x_bound,
                        self.max_step,
                        self.rtol,
                        self.atol,
                        None,
                        false,
                        self.first_step,
                    )
                    .with_max_steps(max_steps)
                    .with_telemetry_mode(match self.telemetry {
                        SolverTelemetryMode::Off => BdfTelemetryMode::Off,
                        SolverTelemetryMode::Counters => BdfTelemetryMode::Counters,
                        SolverTelemetryMode::Timings => BdfTelemetryMode::Timings,
                        SolverTelemetryMode::Detailed => BdfTelemetryMode::Detailed,
                    })
                    .with_logging_level(self.logging),
                })
            }
            ReactorIvpSolverFamily::Radau => {
                if self.stop_condition.is_some() {
                    return Err(ReactorIvpBackendError::UnsupportedCapability {
                        family: self.family(),
                        capability: "species stop condition",
                    });
                }
                if !matches!(self.execution_backend, ReactorIvpExecutionBackend::Lambdify) {
                    return Err(ReactorIvpBackendError::UnsupportedCapability {
                        family: self.family(),
                        capability: "AOT execution mapping for Radau",
                    });
                }
                let (max_steps, max_newton_iterations) = match self.family_options {
                    ReactorIvpFamilyOptions::Radau {
                        max_steps,
                        max_newton_iterations,
                    } => (max_steps, max_newton_iterations),
                    _ => unreachable!("family() and family options are derived together"),
                };
                let problem = RadauProblem::new(eq_system, values, arg);
                let config = RadauConfig {
                    t0: self.x0,
                    t_bound: self.x_bound,
                    rtol: self.rtol,
                    atol: self.atol,
                    first_step: self.first_step,
                    max_step: self.max_step,
                    max_steps,
                    max_newton_iterations,
                    frontend: match self.symbolic_backend {
                        ReactorIvpSymbolicBackend::AtomView => RadauFrontend::AtomViewNative,
                        ReactorIvpSymbolicBackend::ExprLegacy => RadauFrontend::ExprLegacy,
                    },
                    matrix_layout: match self.matrix_backend {
                        ReactorIvpMatrixBackend::Dense => RadauMatrixLayout::Dense,
                        ReactorIvpMatrixBackend::Sparse => RadauMatrixLayout::Sparse,
                        ReactorIvpMatrixBackend::Banded => RadauMatrixLayout::BandedAuto,
                    },
                    telemetry: match self.telemetry {
                        SolverTelemetryMode::Off => RadauTelemetryMode::Off,
                        SolverTelemetryMode::Counters => RadauTelemetryMode::Counters,
                        SolverTelemetryMode::Timings => RadauTelemetryMode::Timings,
                        SolverTelemetryMode::Detailed => RadauTelemetryMode::Detailed,
                    },
                    logging_level: self.logging,
                    output: RadauOutputPolicy::Dense,
                    ..RadauConfig::default()
                };
                Ok(StiffSolverRequest::Radau {
                    problem,
                    config,
                    initial_state: y0.iter().copied().collect(),
                    parameters: Vec::new(),
                })
            }
        }
    }

    /// Prepare the typed request through RustedSciThe's shared stiff-solver
    /// dispatch boundary.
    pub fn prepare_stiff_solver(
        self,
        eq_system: Vec<Expr>,
        values: Vec<String>,
        arg: String,
        y0: DVector<f64>,
    ) -> Result<PreparedStiffSolver, ReactorIvpBackendError> {
        self.prepare_stiff_solver_with_hybrid_settings(eq_system, values, arg, y0, None)
    }

    pub fn prepare_stiff_solver_with_hybrid_settings(
        self,
        eq_system: Vec<Expr>,
        values: Vec<String>,
        arg: String,
        y0: DVector<f64>,
        hybrid_settings: Option<&HybridStiffSettings>,
    ) -> Result<PreparedStiffSolver, ReactorIvpBackendError> {
        let request = self.to_stiff_request(eq_system, values.clone(), arg, y0)?;
        self.prepare_stiff_request(request, &values, hybrid_settings)
    }

    /// Prepare an already-built request through the shared StiffAPI boundary.
    /// This lets callers measure request construction separately from RST setup.
    pub fn prepare_stiff_request(
        self,
        request: StiffSolverRequest,
        values: &[String],
        hybrid_settings: Option<&HybridStiffSettings>,
    ) -> Result<PreparedStiffSolver, ReactorIvpBackendError> {
        let mut prepared = if let Some(settings) = hybrid_settings {
            settings.clone().prepare(request)
        } else {
            request.prepare(Some(&self.stiff_preset()?))
        }
        .map_err(ReactorIvpBackendError::StiffApiPreparation)?;
        if let Some(stop_condition) = self.stop_condition {
            let mut condition = std::collections::HashMap::new();
            condition.insert(
                values[REACTOR_IVP_SPECIES_STATE_OFFSET + stop_condition.species_index].clone(),
                stop_condition.threshold,
            );
            if let Some(solver) = prepared.as_bdf_mut() {
                solver.try_set_stop_condition(condition).map_err(|error| {
                    ReactorIvpBackendError::InvalidConfiguration(error.to_string())
                })?;
            }
        }
        Ok(prepared)
    }

    fn validate_common_inputs(
        self,
        eq_system: &[Expr],
        values: &[String],
        y0: &DVector<f64>,
    ) -> Result<(), ReactorIvpBackendError> {
        if eq_system.len() != values.len() || y0.len() != values.len() {
            return Err(ReactorIvpBackendError::InvalidConfiguration(format!(
                "equation, variable, and initial-state dimensions differ (equations={}, variables={}, initial={})",
                eq_system.len(),
                values.len(),
                y0.len()
            )));
        }
        if self.x_bound <= self.x0 || !self.x0.is_finite() || !self.x_bound.is_finite() {
            return Err(ReactorIvpBackendError::InvalidConfiguration(
                "integration bounds must be finite and x_bound must exceed x0".to_string(),
            ));
        }
        if self.matrix_backend == ReactorIvpMatrixBackend::Banded {
            self.matrix_backend.validate_for_dimension(values.len())?;
        }
        if let Some(condition) = self.stop_condition {
            condition.validate(
                values
                    .len()
                    .checked_sub(REACTOR_IVP_SPECIES_STATE_OFFSET)
                    .ok_or_else(|| {
                        ReactorIvpBackendError::InvalidConfiguration(
                            "condensed IVP requires two thermal states before species".to_string(),
                        )
                    })?,
            )?;
        }
        Ok(())
    }

    fn stiff_preset(self) -> Result<StiffPreset, ReactorIvpBackendError> {
        let mut preset = StiffPreset::default()
            .with_frontend(match self.symbolic_backend {
                ReactorIvpSymbolicBackend::AtomView => StiffFrontend::AtomViewNative,
                ReactorIvpSymbolicBackend::ExprLegacy => StiffFrontend::ExprLegacy,
            })
            .with_telemetry(self.telemetry)
            .with_logging(self.logging);
        preset = match self.matrix_backend {
            ReactorIvpMatrixBackend::Dense => preset.with_matrix_layout(StiffMatrixLayout::Dense),
            ReactorIvpMatrixBackend::Sparse => preset.with_matrix_layout(StiffMatrixLayout::Sparse),
            ReactorIvpMatrixBackend::Banded => {
                preset.with_matrix_layout(StiffMatrixLayout::BandedAuto)
            }
        };
        if matches!(self.execution_backend, ReactorIvpExecutionBackend::Lambdify) {
            preset = preset.with_execution(StiffExecution::Lambdify);
        }
        if let Some(policy) = self.callback_policy {
            preset = preset.with_callback_policy(policy.to_stiff());
        }
        Ok(preset)
    }
}
