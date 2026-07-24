//! Solver-backend facade for condensed-phase IVP tasks.
//!
//! The IVP reactor workflow starts from symbolic equations and typed initial
//! conditions, so the user-facing facade mirrors only the meaningful LSODE2
//! routes: method family, symbolic execution backend, and matrix structure.

use RustedSciThe::numerical::LSODE2::{
    Lsode2AotProfile, Lsode2AotToolchain, Lsode2ControllerConfig, Lsode2LinearSolverPolicy,
    Lsode2LinearSystemStructure, Lsode2Method, Lsode2NativeExecutionConfig, Lsode2ProblemConfig,
    Lsode2ResidualJacobianSource, Lsode2SymbolicAssemblyBackend, Lsode2SymbolicExecutionMode,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DVector;
use std::fmt;

const REACTOR_IVP_SPECIES_STATE_OFFSET: usize = 2;

/// Typed solver-backend error used while building the LSODE2 handoff.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReactorIvpBackendError {
    InvalidConfiguration(String),
}

impl fmt::Display for ReactorIvpBackendError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfiguration(msg) => write!(f, "Invalid configuration: {msg}"),
        }
    }
}

impl std::error::Error for ReactorIvpBackendError {}

/// Time-integration family exposed by the condensed-phase IVP facade.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactorIvpMethod {
    Auto,
    Bdf,
    Adams,
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
            Self::Sparse => Lsode2LinearSystemStructure::Sparse,
            Self::Banded => Lsode2LinearSystemStructure::Banded { kl: 0, ku: 0 },
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

impl Default for ReactorIvpExecutionBackend {
    fn default() -> Self {
        Self::Lambdify
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
    pub method: ReactorIvpMethod,
    pub symbolic_backend: ReactorIvpSymbolicBackend,
    pub matrix_backend: ReactorIvpMatrixBackend,
    pub execution_backend: ReactorIvpExecutionBackend,
    pub stop_condition: Option<ReactorIvpStopCondition>,
    pub x0: f64,
    pub x_bound: f64,
    pub first_step: Option<f64>,
    pub max_step: f64,
    pub rtol: f64,
    pub atol: f64,
    pub native_execution: Lsode2NativeExecutionConfig,
}

impl Default for ReactorIvpSolverConfig {
    fn default() -> Self {
        Self {
            method: ReactorIvpMethod::default(),
            symbolic_backend: ReactorIvpSymbolicBackend::default(),
            matrix_backend: ReactorIvpMatrixBackend::default(),
            execution_backend: ReactorIvpExecutionBackend::default(),
            stop_condition: None,
            x0: 0.0,
            x_bound: 1.0,
            first_step: None,
            max_step: 1.0,
            rtol: 1e-6,
            atol: 1e-8,
            native_execution: Lsode2NativeExecutionConfig::faithful_bdf_solve(200_000, 200_000),
        }
    }
}

impl ReactorIvpSolverConfig {
    pub fn default_lambdify() -> Self {
        Self::default()
    }

    pub fn with_method(mut self, method: ReactorIvpMethod) -> Self {
        self.method = method;
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

    pub fn with_native_execution(mut self, native_execution: Lsode2NativeExecutionConfig) -> Self {
        self.native_execution = native_execution;
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
        config.method = self.method.to_lsode2_method();
        config.first_step = self.first_step;
        config.controller = self.method.to_controller();
        config.residual_jacobian_source = self.to_rusted_residual_source();
        config.linear_system_structure = self.matrix_backend.to_rusted();
        config.linear_solver_policy = self.to_rusted_linear_policy();
        config.native_execution = self.native_execution;
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
}
