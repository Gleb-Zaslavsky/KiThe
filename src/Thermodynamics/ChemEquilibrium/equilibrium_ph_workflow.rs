//! Fixed-pressure, fixed-total-enthalpy workflow.
//!
//! The canonical resolved-data implementation is the coupled monolithic
//! solve. The outer safeguarded scalar solve around
//! [`super::phase_equilibrium_workflow::solve_resolved_pt`] remains an
//! independent reference and recovery path behind an explicit solve-mode
//! boundary.
//!
//! # Supported physical contract
//!
//! This workflow currently supports closed, fixed-pressure systems made from
//! the ideal-gas and pure-condensed phase models accepted by the canonical
//! `P,T` engine. Total enthalpy is the additive standard-state quantity
//! `H = sum_i n_i h_i(T)` in joules. Consequently, a request is physically
//! meaningful only when every selected component record uses a compatible
//! thermochemical reference convention and the supplied `H_target` was built
//! from that same convention.
//!
//! In this scope, a `P,H` request represents an adiabatic isobaric equilibrium
//! calculation only when heat transfer, shaft/electrical work, and kinetic or
//! potential-energy terms are either absent or already included in the target.
//! Non-ideal solution excess enthalpy, pressure-dependent real-fluid
//! enthalpy, and unimplemented phase-model contributions are intentionally not
//! inferred from this additive model.
//!
//! The workflow is **thermochemical-format agnostic**. It never branches on
//! NASA, NIST, or a library name: phase resolution selects records first, then
//! [`ResolvedThermochemistry`] asks each selected `ThermoCalculator` for the
//! same `g0(T)`, `h(T)`, optional `Cp(T)`, and valid-temperature capabilities.
//! Provenance remains in the result precisely so format/library differences
//! stay observable without leaking into numerical residual construction.
//!
//! # Validation scope
//!
//! Optional equilibrium-constant validation remains attached to an inner
//! fixed-`P,T` chemical candidate. It is evidence about that candidate only;
//! outer `P,H` acceptance is independently determined from the typed enthalpy
//! residual, scalar bracket, conservation, and normal `P,T` acceptance gate.

use std::cell::RefCell;
use std::rc::Rc;
use std::time::{Duration, Instant};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    additive_total_enthalpy, EnthalpyScale, EquilibriumConstraint, TemperatureBounds,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::{
    EquilibriumExecutionControl, EquilibriumProgressEvent, EquilibriumProgressStage,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionExtentError, ReactionExtentErrorKind,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_formulation::PreparedPhFormulation;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_monolithic::{
    solve_monolithic_active_set_candidate, PreparedMonolithicPhRunner,
};
#[cfg(test)]
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_nested::safeguarded_interpolation_step;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_nested::{
    solve_bracketed_temperature as solve_nested_bracketed_temperature, NestedBracketOptions,
};
pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_nested::{
    PhTemperatureStepKind, PhTemperatureTrial, PhTrialInnerEvidence, PhTrialPhaseState,
    PhTrialPreparation, PhTrialTimingReport,
};
pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_options::{
    PhAcceptanceOptions, PhMonolithicOptions, PhMonotonicityPolicy, PhNestedOptions,
};
pub use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::{
    MolarEnthalpyFunction, MolarThermoFunction, ResolvedThermochemistry, ThermochemistryProvenance,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::LogMolesInitialGuess;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    prepare_rst_symbolic_ph_problem, RstPreparedProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, MultiStartSolveReport,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::{
    EquilibriumTimingMode, EquilibriumTimingReport,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    MultiphaseAcceptanceReport, PhaseControlledSolveReport,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    build_phase_equilibrium_problem_with_timing, PhaseEquilibriumBuildRequest,
    PreparedPhaseEquilibriumTemplate, SupportedPhaseModelPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    solve_resolved_pt, EquilibriumSolveOptions, PhaseControlPolicy, PhaseEquilibriumSolveMode,
    ResolvedPhaseEquilibriumRequest,
};
use crate::Thermodynamics::User_PhaseOrSolution::ResolvedPhaseSystem;

impl ResolvedThermochemistry {
    /// Converts the mandatory enthalpy column to the compatibility model used
    /// by the independent nested P,H route. The capability bundle remains the
    /// source of truth; this is only an adapter for that older input shape.
    pub fn enthalpy_model(&self) -> EnthalpyModel<'static> {
        EnthalpyModel {
            functions: self.enthalpy_functions(),
            heat_capacity: self.heat_capacity_functions(),
        }
    }
}

/// Component-aligned molar enthalpy capabilities for the additive model.
///
/// The component order must be the same as the resolved `SystemLayout` and
/// the initial composition. The resolved-data bridge is
/// `ResolvedThermochemistry`; this smaller type remains useful for analytic
/// and adapter-level callers.
#[derive(Clone)]
pub struct EnthalpyModel<'a> {
    functions: Vec<MolarEnthalpyFunction<'a>>,
    heat_capacity: Vec<Option<MolarEnthalpyFunction<'a>>>,
}

/// Pure additive enthalpy evidence at one temperature.
///
/// `partial_temperature_derivative` deliberately means `dH/dT` at fixed
/// composition. It is not the full outer derivative of
/// `H(solution_of_P_T(T), T)`, which also contains the implicit composition
/// term `sum(h_i * d n_i / dT)`.
#[derive(Debug, Clone, PartialEq)]
pub struct EnthalpyEvaluation {
    molar_enthalpies: Vec<f64>,
    heat_capacities: Vec<f64>,
    total_enthalpy: f64,
    partial_temperature_derivative: f64,
}

/// Numerical formulation selected for one fixed-pressure, fixed-enthalpy solve.
///
/// The coupled monolithic path is the canonical default for requests built
/// from resolved thermochemistry. The safeguarded nested path remains an
/// explicit independent reference/fallback route, while `Auto` records a
/// monolithic attempt and recovers through nested solving only after a
/// classified formulation failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhSolveMode {
    /// Solve composition and bounded temperature in one nonlinear system.
    Monolithic,
    /// Safeguarded scalar temperature solve around accepted fixed-P,T solves.
    NestedTemperature,
    /// Monolithic first, then a nested fallback only for classified numerical
    /// failures.
    Auto,
}

impl Default for PhSolveMode {
    fn default() -> Self {
        Self::Monolithic
    }
}

/// Numerical route that produced an accepted `P,H` result.
///
/// This is result evidence, not a caller preference. An `Auto` request may
/// legitimately publish a nested recovery result after a classified
/// monolithic failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhSolvePath {
    /// One coupled nonlinear solve over log-moles and bounded temperature for
    /// a fixed declared active phase set.
    MonolithicFixedActiveSet,
    /// Coupled log-mole/temperature candidates executed by the shared
    /// transactional bounded phase-control lifecycle.
    MonolithicPhaseControl,
    /// Safeguarded scalar temperature search containing accepted P,T solves.
    NestedTemperature,
}

/// Why an automatic P,H request left the monolithic formulation.
///
/// The reason is retained in the immutable result rather than being emitted
/// only to logs. This makes `Auto` auditable and prevents a numerical recovery
/// from looking like an ordinary nested solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhFallbackReason {
    error_kind: ReactionExtentErrorKind,
    message: String,
}

impl PhFallbackReason {
    /// Classifies a [`ReactionExtentError`] into a typed fallback reason.
    ///
    /// Captures the high-level error kind and the full display message so the
    /// monolithic-to-nested recovery path can report *why* the coupled solve
    /// failed without retaining the original error's dynamic type or backtrace.
    /// The classification is deliberately lossy: only the kind and message
    /// survive, which is sufficient for diagnostics and summary rows.
    fn from_error(error: &ReactionExtentError) -> Self {
        Self {
            error_kind: error.kind(),
            message: error.to_string(),
        }
    }

    /// High-level classification of the rejected monolithic attempt.
    pub fn error_kind(&self) -> ReactionExtentErrorKind {
        self.error_kind
    }

    /// Human-readable typed error rendered at the fallback boundary.
    pub fn message(&self) -> &str {
        &self.message
    }
}

impl EnthalpyEvaluation {
    /// Component-aligned molar enthalpies in J/mol.
    pub fn molar_enthalpies(&self) -> &[f64] {
        &self.molar_enthalpies
    }

    /// Component-aligned heat capacities in J/(mol K).
    pub fn heat_capacities(&self) -> &[f64] {
        &self.heat_capacities
    }

    /// Additive total enthalpy in J.
    pub fn total_enthalpy(&self) -> f64 {
        self.total_enthalpy
    }

    /// `dH/dT` while holding the supplied component moles fixed.
    pub fn partial_temperature_derivative(&self) -> f64 {
        self.partial_temperature_derivative
    }
}

impl<'a> EnthalpyModel<'a> {
    /// Creates a model from functions returning J/mol at temperature in K.
    pub fn from_functions(
        functions: Vec<MolarEnthalpyFunction<'a>>,
    ) -> Result<Self, ReactionExtentError> {
        if functions.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_functions",
                message: "at least one molar enthalpy function is required".to_string(),
            });
        }
        Ok(Self {
            heat_capacity: vec![None; functions.len()],
            functions,
        })
    }

    /// Creates an enthalpy model with component-aligned heat-capacity
    /// capabilities for derivative-ready callers.
    ///
    /// A `None` entry is retained intentionally: it means the source record
    /// has no usable `Cp(T)` capability. The ordinary P,H bracketed workflow
    /// may still use such a model, while derivative evaluation rejects it at
    /// the point where the missing capability matters.
    pub fn from_functions_with_heat_capacity(
        functions: Vec<MolarEnthalpyFunction<'a>>,
        heat_capacity: Vec<Option<MolarEnthalpyFunction<'a>>>,
    ) -> Result<Self, ReactionExtentError> {
        let model = Self::from_functions(functions)?;
        if model.functions.len() != heat_capacity.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "enthalpy model has {} functions but {} heat-capacity functions",
                model.functions.len(),
                heat_capacity.len()
            )));
        }
        Ok(Self {
            heat_capacity,
            ..model
        })
    }

    /// Builds component-aligned `dH(T)` capabilities from resolved phase data.
    ///
    /// This compatibility constructor delegates to the complete bundle so it
    /// cannot accidentally trust a stale zero-valued closure cache. The
    /// returned model retains the bundle's temperature-aware calculator path.
    pub fn from_resolved_system(
        resolved: &'a ResolvedPhaseSystem,
    ) -> Result<Self, ReactionExtentError> {
        Ok(ResolvedThermochemistry::from_resolved_system(resolved)?.enthalpy_model())
    }

    /// Number of component-aligned capabilities.
    pub fn len(&self) -> usize {
        self.functions.len()
    }

    /// Whether this model contains no capabilities.
    pub fn is_empty(&self) -> bool {
        self.functions.is_empty()
    }

    /// Evaluates all molar enthalpies at one temperature in K.
    pub fn evaluate_molar(&self, temperature: f64) -> Result<Vec<f64>, ReactionExtentError> {
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature",
                message: "enthalpy temperature must be finite and positive".to_string(),
            });
        }

        self.functions
            .iter()
            .enumerate()
            .map(|(index, function)| {
                let value = function(temperature)?;
                if !value.is_finite() {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "molar_enthalpy",
                        message: format!("enthalpy function {index} returned a non-finite value"),
                    });
                }
                Ok(value)
            })
            .collect()
    }

    /// Evaluates the additive total enthalpy in J for a physical composition.
    pub fn evaluate_total(
        &self,
        moles: &[f64],
        temperature: f64,
    ) -> Result<f64, ReactionExtentError> {
        let molar = self.evaluate_molar(temperature)?;
        additive_total_enthalpy(moles, &molar)
    }

    /// Evaluates additive enthalpy and its explicit temperature derivative.
    ///
    /// This method does not solve the chemical equilibrium sensitivity
    /// problem. It is therefore safe to use as a building block for a future
    /// implicit derivative, where the missing term is obtained from the
    /// inner residual Jacobian and its temperature column.
    pub fn evaluate_total_with_partial_temperature_derivative(
        &self,
        moles: &[f64],
        temperature: f64,
    ) -> Result<EnthalpyEvaluation, ReactionExtentError> {
        let molar_enthalpies = self.evaluate_molar(temperature)?;
        let total_enthalpy = additive_total_enthalpy(moles, &molar_enthalpies)?;
        let mut heat_capacities = Vec::with_capacity(self.heat_capacity.len());
        for (index, function) in self.heat_capacity.iter().enumerate() {
            let Some(function) = function else {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "heat_capacity",
                    message: format!("heat capacity capability {index} is unavailable"),
                });
            };
            let value = function(temperature)?;
            if !value.is_finite() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "heat_capacity",
                    message: format!(
                        "heat capacity capability {index} returned a non-finite value"
                    ),
                });
            }
            heat_capacities.push(value);
        }
        let partial_temperature_derivative = moles
            .iter()
            .zip(&heat_capacities)
            .map(|(&moles_i, &cp_i)| moles_i * cp_i)
            .sum::<f64>();
        if !partial_temperature_derivative.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_temperature_derivative",
                message: "partial enthalpy temperature derivative is not finite".to_string(),
            });
        }
        Ok(EnthalpyEvaluation {
            molar_enthalpies,
            heat_capacities,
            total_enthalpy,
            partial_temperature_derivative,
        })
    }
}

/// Scalar solver controls for the outer `P,H` temperature problem.
#[derive(Debug, Clone)]
pub struct PhTemperatureSolveOptions {
    /// Dimensionless enthalpy residual tolerance after scaling.
    pub scaled_enthalpy_tolerance: f64,
    /// Absolute enthalpy tolerance in joules. The effective acceptance limit
    /// is the larger of this floor and the scale-aware relative tolerance.
    pub absolute_enthalpy_tolerance_joules: f64,
    /// Temperature interval tolerance in K.
    pub temperature_tolerance: f64,
    /// Maximum number of interior temperature trials.
    pub max_iterations: usize,
    /// Global maximum number of scalar temperature evaluations, including
    /// both bracket endpoints. This is not multiplied by inner backend work.
    pub max_temperature_evaluations: usize,
    /// Policy for sampled reversals of the outer `H_eq(P,T)` branch.
    pub monotonicity_policy: PhMonotonicityPolicy,
    /// Optional global budget for started inner nonlinear backend attempts.
    ///
    /// This is deliberately separate from the scalar-evaluation budget: one
    /// temperature trial may use several fallback backends, so a bounded
    /// outer solve must be able to cap that nested work explicitly.
    pub max_inner_backend_attempts: Option<usize>,
    /// Optional global budget for nonlinear iterations across all inner P,T
    /// solves and all continuation seeds.
    pub max_inner_nonlinear_iterations: Option<usize>,
    /// Optional global budget for accepted phase transitions across all inner
    /// bounded phase-control solves.
    pub max_phase_control_transitions: Option<usize>,
    /// Optional wall-clock budget for the complete outer transaction.
    pub max_wall_time: Option<Duration>,
    /// Optional cancellation/progress control shared with every inner P,T
    /// solve in this outer transaction.
    execution_control: Option<EquilibriumExecutionControl>,
}

impl Default for PhTemperatureSolveOptions {
    fn default() -> Self {
        Self {
            scaled_enthalpy_tolerance: 1.0e-8,
            absolute_enthalpy_tolerance_joules: 1.0e-6,
            temperature_tolerance: 1.0e-8,
            max_iterations: 80,
            max_temperature_evaluations: 82,
            monotonicity_policy: PhMonotonicityPolicy::default(),
            max_inner_backend_attempts: None,
            max_inner_nonlinear_iterations: None,
            max_phase_control_transitions: None,
            max_wall_time: None,
            execution_control: None,
        }
    }
}

impl PhTemperatureSolveOptions {
    /// Limits the total number of started inner backend attempts across all
    /// scalar temperature trials.
    ///
    /// The limit is checked after each completed inner solve. A failed inner
    /// solve is still accounted for by its own backend policy and remains
    /// visible as the original typed error.
    pub fn with_max_inner_backend_attempts(
        mut self,
        limit: usize,
    ) -> Result<Self, ReactionExtentError> {
        if limit == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_inner_backend_attempts",
                message: "inner backend-attempt budget must be positive".to_string(),
            });
        }
        self.max_inner_backend_attempts = Some(limit);
        Ok(self)
    }

    /// Removes the optional global inner backend-attempt limit.
    pub fn without_inner_backend_attempt_limit(mut self) -> Self {
        self.max_inner_backend_attempts = None;
        self
    }

    /// Limits total nonlinear iterations across the complete outer solve.
    pub fn with_max_inner_nonlinear_iterations(
        mut self,
        limit: usize,
    ) -> Result<Self, ReactionExtentError> {
        if limit == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_inner_nonlinear_iterations",
                message: "inner nonlinear-iteration budget must be positive".to_string(),
            });
        }
        self.max_inner_nonlinear_iterations = Some(limit);
        Ok(self)
    }

    /// Removes the optional global nonlinear-iteration limit.
    pub fn without_inner_nonlinear_iteration_limit(mut self) -> Self {
        self.max_inner_nonlinear_iterations = None;
        self
    }

    /// Limits accepted phase transitions across the complete outer solve.
    pub fn with_max_phase_control_transitions(
        mut self,
        limit: usize,
    ) -> Result<Self, ReactionExtentError> {
        if limit == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_phase_control_transitions",
                message: "phase-transition budget must be positive".to_string(),
            });
        }
        self.max_phase_control_transitions = Some(limit);
        Ok(self)
    }

    /// Removes the optional global phase-transition limit.
    pub fn without_phase_control_transition_limit(mut self) -> Self {
        self.max_phase_control_transitions = None;
        self
    }

    /// Sets the explicit policy for sampled outer-branch reversals.
    pub fn with_monotonicity_policy(mut self, policy: PhMonotonicityPolicy) -> Self {
        self.monotonicity_policy = policy;
        self
    }

    /// Validates scalar solver controls before any inner equilibrium solve.
    pub fn validate(&self) -> Result<(), ReactionExtentError> {
        if !self.scaled_enthalpy_tolerance.is_finite() || self.scaled_enthalpy_tolerance <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "scaled_enthalpy_tolerance",
                message: "tolerance must be finite and positive".to_string(),
            });
        }
        if !self.absolute_enthalpy_tolerance_joules.is_finite()
            || self.absolute_enthalpy_tolerance_joules <= 0.0
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "absolute_enthalpy_tolerance_joules",
                message: "tolerance must be finite and positive".to_string(),
            });
        }
        if !self.temperature_tolerance.is_finite() || self.temperature_tolerance <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_tolerance",
                message: "tolerance must be finite and positive".to_string(),
            });
        }
        if self.max_iterations == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_iterations",
                message: "maximum iterations must be greater than zero".to_string(),
            });
        }
        if self.max_temperature_evaluations < 2 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_temperature_evaluations",
                message: "at least two evaluations are required to test a bracket".to_string(),
            });
        }
        if matches!(self.max_inner_backend_attempts, Some(limit) if limit == 0) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_inner_backend_attempts",
                message: "inner backend-attempt budget must be positive when provided".to_string(),
            });
        }
        if matches!(self.max_inner_nonlinear_iterations, Some(limit) if limit == 0) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_inner_nonlinear_iterations",
                message: "inner nonlinear-iteration budget must be positive when provided"
                    .to_string(),
            });
        }
        if matches!(self.max_phase_control_transitions, Some(limit) if limit == 0) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_phase_control_transitions",
                message: "phase-transition budget must be positive when provided".to_string(),
            });
        }
        if matches!(self.max_wall_time, Some(limit) if limit.is_zero()) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_wall_time",
                message: "wall-time budget must be positive when provided".to_string(),
            });
        }
        Ok(())
    }

    /// Returns the absolute joule limit corresponding to a validated scale.
    ///
    /// Using the maximum makes the contract meaningful both for ordinary
    /// extensive enthalpies and for targets close to zero. A relative-only
    /// test would otherwise demand an unrealistically tiny absolute residual
    /// whenever the reference scale is small.
    pub fn accepted_enthalpy_error_limit_joules(&self, scale: EnthalpyScale) -> f64 {
        self.absolute_enthalpy_tolerance_joules
            .max(self.scaled_enthalpy_tolerance * scale.joules())
    }

    /// Tests the complete absolute-plus-relative enthalpy acceptance contract.
    pub fn accepts_enthalpy_error(&self, error_joules: f64, scale: EnthalpyScale) -> bool {
        error_joules.is_finite()
            && error_joules.abs() <= self.accepted_enthalpy_error_limit_joules(scale)
    }

    /// Installs cooperative cancellation and progress reporting for the
    /// complete outer transaction. The same handle is forwarded to inner P,T
    /// attempts, so cancellation cannot leave an inner solve running after
    /// the outer caller has requested shutdown.
    pub fn with_execution_control(mut self, control: EquilibriumExecutionControl) -> Self {
        self.execution_control = Some(control);
        self
    }

    /// Returns the shared execution control, if configured.
    pub fn execution_control(&self) -> Option<&EquilibriumExecutionControl> {
        self.execution_control.as_ref()
    }

    /// Returns the common enthalpy acceptance contract shared by both routes.
    ///
    /// The fields remain public for the existing builder surface, so this
    /// method validates them again before exposing the typed value object.
    pub fn acceptance_options(&self) -> Result<PhAcceptanceOptions, ReactionExtentError> {
        PhAcceptanceOptions::new(
            self.scaled_enthalpy_tolerance,
            self.absolute_enthalpy_tolerance_joules,
        )
    }

    /// Returns controls consumed only by the coupled monolithic runner.
    pub fn monolithic_options(&self) -> Result<PhMonolithicOptions, ReactionExtentError> {
        Ok(PhMonolithicOptions::new(self.acceptance_options()?))
    }

    /// Returns controls consumed only by the nested scalar/bracket runner.
    pub fn nested_options(&self) -> Result<PhNestedOptions, ReactionExtentError> {
        PhNestedOptions::new(
            self.acceptance_options()?,
            self.temperature_tolerance,
            self.max_iterations,
            self.max_temperature_evaluations,
            self.monotonicity_policy,
            self.max_wall_time,
        )
    }
}

/// Evidence from the outer scalar temperature solve.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PhTemperatureTimingReport {
    enabled: bool,
    total: Duration,
    scalar_orchestration: Duration,
    enthalpy_evaluation: Duration,
}

impl PhTemperatureTimingReport {
    /// Whether this report contains measured outer-stage timings.
    pub fn enabled(&self) -> bool {
        self.enabled
    }

    /// Wall-clock duration of the complete outer P,H transaction.
    pub fn total(&self) -> Duration {
        self.total
    }

    /// Time spent in bracket orchestration, including nested P,T work.
    pub fn scalar_orchestration(&self) -> Duration {
        self.scalar_orchestration
    }

    /// Aggregate time spent evaluating additive total enthalpy at trials.
    pub fn enthalpy_evaluation(&self) -> Duration {
        self.enthalpy_evaluation
    }
}

/// Backend and active-set evidence for a coupled monolithic P,H solve.
///
/// This is deliberately separate from [`PhTemperatureTrial`]. A monolithic
/// solve has no outer scalar temperature trials: its Newton/LM/TR iterations
/// operate on one coupled `[log-moles, temperature]` system. Nested scalar
/// diagnostics therefore must not be used as a container for this evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct PhMonolithicEvidence {
    solve_report: EquilibriumSolveReport,
    multi_start_report: Option<MultiStartSolveReport>,
    phase_control_report: Option<PhaseControlledSolveReport>,
    acceptance_report: Option<MultiphaseAcceptanceReport>,
    inner_timing: EquilibriumTimingReport,
}

impl PhMonolithicEvidence {
    /// Ordered backend attempts for the coupled monolithic candidate.
    pub fn solve_report(&self) -> &EquilibriumSolveReport {
        &self.solve_report
    }

    /// Seed comparison evidence, when the inner policy used multi-start.
    pub fn multi_start_report(&self) -> Option<&MultiStartSolveReport> {
        self.multi_start_report.as_ref()
    }

    /// Transactional phase-control evidence for bounded monolithic mode.
    pub fn phase_control_report(&self) -> Option<&PhaseControlledSolveReport> {
        self.phase_control_report.as_ref()
    }

    /// Final fixed-P,T acceptance evidence carried by the monolithic solve.
    pub fn acceptance_report(&self) -> Option<&MultiphaseAcceptanceReport> {
        self.acceptance_report.as_ref()
    }

    /// Timing report for the coupled inner equilibrium solve.
    pub fn inner_timing(&self) -> EquilibriumTimingReport {
        self.inner_timing
    }

    /// Residual evaluations reported by the coupled backend cascade.
    pub fn residual_evaluations(&self) -> usize {
        self.solve_report
            .attempts
            .iter()
            .filter_map(|attempt| attempt.metrics.as_ref())
            .map(|metrics| metrics.residual_evaluations)
            .sum()
    }

    /// Jacobian evaluations reported by the coupled backend cascade.
    pub fn jacobian_evaluations(&self) -> usize {
        self.solve_report
            .attempts
            .iter()
            .filter_map(|attempt| attempt.metrics.as_ref())
            .map(|metrics| metrics.jacobian_evaluations)
            .sum()
    }
}

/// Evidence from the outer scalar temperature solve.
#[derive(Debug, Clone, PartialEq)]
pub struct PhTemperatureSolveReport {
    solve_path: PhSolvePath,
    fallback_reason: Option<PhFallbackReason>,
    target_enthalpy: f64,
    enthalpy_scale_joules: f64,
    absolute_enthalpy_tolerance_joules: f64,
    scaled_enthalpy_tolerance: f64,
    max_temperature_evaluations: usize,
    initial_temperature_seed: f64,
    monotonicity_policy: PhMonotonicityPolicy,
    max_inner_backend_attempts: Option<usize>,
    max_inner_nonlinear_iterations: Option<usize>,
    max_phase_control_transitions: Option<usize>,
    max_wall_time: Option<Duration>,
    inner_backend_attempts: usize,
    inner_nonlinear_iterations: usize,
    phase_control_transitions: usize,
    fixed_formulation_builds: usize,
    fixed_formulation_reuses: usize,
    inner_timing: EquilibriumTimingReport,
    timing: PhTemperatureTimingReport,
    iterations: usize,
    trials: Vec<PhTemperatureTrial>,
    monolithic_evidence: Option<PhMonolithicEvidence>,
}

impl PhTemperatureSolveReport {
    /// Numerical route that accepted this immutable result.
    pub fn solve_path(&self) -> PhSolvePath {
        self.solve_path
    }

    /// Numerical reason that caused `Auto` to recover through the nested
    /// formulation, if recovery occurred.
    pub fn fallback_reason(&self) -> Option<&PhFallbackReason> {
        self.fallback_reason.as_ref()
    }

    /// Target total enthalpy in J.
    pub fn target_enthalpy(&self) -> f64 {
        self.target_enthalpy
    }

    /// Scale used to make the energy residual dimensionless, in J.
    pub fn enthalpy_scale_joules(&self) -> f64 {
        self.enthalpy_scale_joules
    }

    /// Absolute tolerance floor used by the accepted result contract.
    pub fn absolute_enthalpy_tolerance_joules(&self) -> f64 {
        self.absolute_enthalpy_tolerance_joules
    }

    /// Scale-aware relative tolerance used by the accepted result contract.
    pub fn scaled_enthalpy_tolerance(&self) -> f64 {
        self.scaled_enthalpy_tolerance
    }

    /// Global scalar evaluation budget configured for this solve.
    pub fn max_temperature_evaluations(&self) -> usize {
        self.max_temperature_evaluations
    }

    /// User-supplied P,H temperature seed evaluated inside the explicit
    /// thermochemistry bounds. It is evidence only and never widens bounds.
    pub fn initial_temperature_seed(&self) -> f64 {
        self.initial_temperature_seed
    }

    /// Outer-branch policy used while sampling the temperature bracket.
    pub fn monotonicity_policy(&self) -> PhMonotonicityPolicy {
        self.monotonicity_policy
    }

    /// Optional global budget for started inner backend attempts.
    pub fn max_inner_backend_attempts(&self) -> Option<usize> {
        self.max_inner_backend_attempts
    }

    /// Optional global nonlinear-iteration budget configured for this solve.
    pub fn max_inner_nonlinear_iterations(&self) -> Option<usize> {
        self.max_inner_nonlinear_iterations
    }

    /// Optional global phase-transition budget configured for this solve.
    pub fn max_phase_control_transitions(&self) -> Option<usize> {
        self.max_phase_control_transitions
    }

    /// Optional global wall-time budget configured for this solve.
    pub fn max_wall_time(&self) -> Option<Duration> {
        self.max_wall_time
    }

    /// Total inner backend attempts across all successful temperature trials.
    pub fn inner_backend_attempts(&self) -> usize {
        self.inner_backend_attempts
    }

    /// Total nonlinear iterations across all successful temperature trials.
    pub fn inner_nonlinear_iterations(&self) -> usize {
        self.inner_nonlinear_iterations
    }

    /// Total accepted phase-control transitions across all successful trials.
    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_transitions
    }

    /// Number of fixed-phase structural formulations built by this outer
    /// transaction. Bounded phase control reports zero because its active set
    /// has a separate lifecycle and is never reused across arbitrary bracket
    /// trial ordering.
    pub fn fixed_formulation_builds(&self) -> usize {
        self.fixed_formulation_builds
    }

    /// Number of accepted P,H trials that reused the fixed declared-phase
    /// formulation after its first local construction.
    pub fn fixed_formulation_reuses(&self) -> usize {
        self.fixed_formulation_reuses
    }

    /// Aggregate inner P,T timing across all successful temperature trials.
    pub fn inner_timing(&self) -> EquilibriumTimingReport {
        self.inner_timing
    }

    /// Outer scalar and enthalpy timing, when explicitly enabled.
    pub fn timing(&self) -> PhTemperatureTimingReport {
        self.timing
    }

    /// Effective absolute acceptance limit in joules.
    pub fn accepted_enthalpy_error_limit_joules(&self) -> f64 {
        self.absolute_enthalpy_tolerance_joules
            .max(self.scaled_enthalpy_tolerance * self.enthalpy_scale_joules)
    }

    /// Number of interior bisection trials.
    pub fn iterations(&self) -> usize {
        self.iterations
    }

    /// All endpoint and interior scalar trial records in evaluation order.
    ///
    /// This is empty for monolithic mode because a coupled Newton/LM/TR solve
    /// does not perform an outer scalar temperature search. Use
    /// [`Self::monolithic_evidence`] for that route's backend diagnostics.
    pub fn trials(&self) -> &[PhTemperatureTrial] {
        &self.trials
    }

    /// Coupled backend and phase-control evidence for monolithic mode.
    pub fn monolithic_evidence(&self) -> Option<&PhMonolithicEvidence> {
        self.monolithic_evidence.as_ref()
    }
}

/// Immutable result of one accepted `P,H` solve.
#[derive(Debug, Clone, PartialEq)]
pub struct FixedPressureEnthalpySolution {
    solution: MultiphaseEquilibriumSolution,
    calculated_enthalpy: f64,
    target_enthalpy: f64,
    report: PhTemperatureSolveReport,
    thermochemistry: Option<ResolvedThermochemistry>,
}

impl FixedPressureEnthalpySolution {
    /// Accepted fixed-`P,T` inner solution.
    pub fn equilibrium(&self) -> &MultiphaseEquilibriumSolution {
        &self.solution
    }

    /// Solved equilibrium temperature in K.
    pub fn temperature(&self) -> f64 {
        self.solution.conditions().temperature()
    }

    /// Accepted pressure in Pa.
    pub fn pressure(&self) -> f64 {
        self.solution.conditions().pressure()
    }

    /// Calculated total enthalpy in J.
    pub fn calculated_enthalpy(&self) -> f64 {
        self.calculated_enthalpy
    }

    /// Requested total enthalpy in J.
    pub fn target_enthalpy(&self) -> f64 {
        self.target_enthalpy
    }

    /// Raw energy error in J.
    pub fn enthalpy_error(&self) -> f64 {
        self.calculated_enthalpy - self.target_enthalpy
    }

    /// Scale-aware energy error used by the outer acceptance contract.
    pub fn scaled_enthalpy_error(&self) -> f64 {
        self.enthalpy_error() / self.report.enthalpy_scale_joules()
    }

    /// Effective absolute energy tolerance used to accept this result.
    pub fn enthalpy_error_limit_joules(&self) -> f64 {
        self.report.accepted_enthalpy_error_limit_joules()
    }

    /// Outer scalar-solver diagnostics.
    pub fn report(&self) -> &PhTemperatureSolveReport {
        &self.report
    }

    /// Thermochemistry bundle used for the accepted energy evaluation, when
    /// the request was built from resolved thermochemical records.
    pub fn thermochemistry(&self) -> Option<&ResolvedThermochemistry> {
        self.thermochemistry.as_ref()
    }
}

/// Complete request for one typed `P,H` workflow.
#[derive(Clone)]
pub struct ResolvedPhaseEnthalpyRequest<'a> {
    /// Immutable snapshot captured when the request is built. The solve does
    /// not borrow the caller's mutable lifecycle or depend on its storage.
    resolved: ResolvedPhaseSystem,
    initial_composition: MultiphaseInitialComposition,
    constraint: EquilibriumConstraint,
    temperature_bounds: TemperatureBounds,
    enthalpy: EnthalpyModel<'a>,
    solve_options: EquilibriumSolveOptions,
    solve_mode: PhaseEquilibriumSolveMode,
    ph_solve_mode: PhSolveMode,
    temperature_options: PhTemperatureSolveOptions,
    thermochemistry: Option<ResolvedThermochemistry>,
}

impl<'a> ResolvedPhaseEnthalpyRequest<'a> {
    /// Creates the legacy/reference request from arbitrary enthalpy closures.
    ///
    /// This constructor intentionally selects [`PhSolveMode::NestedTemperature`]
    /// and cannot construct the coupled monolithic formulation because it has
    /// no resolved Gibbs/Cp bundle. New production callers must use
    /// [`Self::from_resolved_thermochemistry`] instead. The constructor is
    /// retained only for compatibility and for the independent nested
    /// reference route.
    #[deprecated(
        note = "use ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry for the canonical production API"
    )]
    pub fn new<'resolved>(
        resolved: &'resolved ResolvedPhaseSystem,
        initial_composition: MultiphaseInitialComposition,
        constraint: EquilibriumConstraint,
        temperature_bounds: TemperatureBounds,
        enthalpy: EnthalpyModel<'a>,
    ) -> Result<Self, ReactionExtentError> {
        Self::from_legacy_parts(
            resolved,
            initial_composition,
            constraint,
            temperature_bounds,
            enthalpy,
        )
    }

    /// Shared construction helper used by both the deprecated [`new`] and the
    /// canonical [`from_resolved_thermochemistry`] constructors.
    ///
    /// Validates the component count against the resolved layout, checks that
    /// the enthalpy model dimension matches, and stores the immutable request
    /// fields without touching the repository or evaluating thermochemistry.
    /// The caller is responsible for attaching the optional
    /// [`ResolvedThermochemistry`] bundle when the thermochemistry-aware path
    /// is used.
    fn from_legacy_parts<'resolved>(
        resolved: &'resolved ResolvedPhaseSystem,
        initial_composition: MultiphaseInitialComposition,
        constraint: EquilibriumConstraint,
        temperature_bounds: TemperatureBounds,
        enthalpy: EnthalpyModel<'a>,
    ) -> Result<Self, ReactionExtentError> {
        let request = Self {
            resolved: resolved.clone(),
            initial_composition,
            constraint,
            temperature_bounds,
            enthalpy,
            solve_options: EquilibriumSolveOptions::default(),
            solve_mode: PhaseEquilibriumSolveMode::fixed_declared_phases(),
            // This constructor accepts arbitrary enthalpy closures and does
            // not have the Gibbs/Cp bundle required by the coupled system.
            // Keep it on the independent reference route; resolved-data
            // callers are promoted to the monolithic default below.
            ph_solve_mode: PhSolveMode::NestedTemperature,
            temperature_options: PhTemperatureSolveOptions::default(),
            thermochemistry: None,
        };
        request.validate()?;
        Ok(request)
    }

    /// Creates a P,H request from the resolved thermochemistry bundle.
    ///
    /// The user may narrow the bundle interval, but cannot extend it. This
    /// keeps scalar trials inside the domain supported by every selected
    /// thermochemical record.
    pub fn from_resolved_thermochemistry<'resolved>(
        resolved: &'resolved ResolvedPhaseSystem,
        initial_composition: MultiphaseInitialComposition,
        constraint: EquilibriumConstraint,
        temperature_bounds: TemperatureBounds,
        thermochemistry: ResolvedThermochemistry,
    ) -> Result<Self, ReactionExtentError> {
        let bundle_bounds = thermochemistry.temperature_bounds();
        if temperature_bounds.lower() < bundle_bounds.lower()
            || temperature_bounds.upper() > bundle_bounds.upper()
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_bounds",
                message: "P,H bounds extend beyond the selected thermochemistry domain".to_string(),
            });
        }
        let enthalpy = thermochemistry.enthalpy_model();
        let mut request = Self::from_legacy_parts(
            resolved,
            initial_composition,
            constraint,
            temperature_bounds,
            enthalpy,
        )?;
        request.thermochemistry = Some(thermochemistry);
        request.ph_solve_mode = PhSolveMode::default();
        Ok(request)
    }

    /// Validates the complete `P,H` request before any inner equilibrium solve.
    ///
    /// Checks that the constraint is a valid `P,H` specification (not `P,T`),
    /// that the initial temperature seed lies inside the declared thermochemistry
    /// bounds, and that the initial composition dimension matches the resolved
    /// layout. Returns a typed error with the offending field name on the first
    /// violation so callers receive a diagnostic message rather than a silent
    /// assertion failure during the scalar search.
    fn validate(&self) -> Result<(), ReactionExtentError> {
        let (target_enthalpy, initial_temperature) = validated_ph_parameters(self.constraint)?;
        self.constraint
            .validate_temperature(initial_temperature, self.temperature_bounds)?;
        if self.initial_composition.moles().len() != self.resolved.layout().component_count() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial composition has {} entries but resolved layout has {} components",
                self.initial_composition.moles().len(),
                self.resolved.layout().component_count()
            )));
        }
        if self.enthalpy.len() != self.initial_composition.moles().len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "enthalpy model has {} functions but resolved layout has {} components",
                self.enthalpy.len(),
                self.initial_composition.moles().len()
            )));
        }
        self.temperature_options.validate()?;
        let initial_molar = self.enthalpy.evaluate_molar(initial_temperature)?;
        EnthalpyScale::from_magnitudes(
            target_enthalpy,
            self.initial_composition.moles(),
            &initial_molar,
        )?;
        Ok(())
    }

    /// Replaces inner fixed-`P,T` solver controls.
    pub fn with_solve_options(mut self, options: EquilibriumSolveOptions) -> Self {
        self.solve_options = options;
        self
    }

    /// Uses bounded active-set phase control for each temperature trial.
    pub fn with_phase_control_policy(mut self, policy: PhaseControlPolicy) -> Self {
        self.solve_mode = PhaseEquilibriumSolveMode::bounded_phase_control(policy);
        self
    }

    /// Selects the P,H numerical formulation explicitly.
    ///
    /// `Monolithic` requires a complete [`ResolvedThermochemistry`] bundle.
    /// With bounded phase control, the coupled candidate is still executed by
    /// the shared transactional active-set lifecycle. `Auto` tries this
    /// monolithic route first and falls back only for retryable numerical
    /// failures; its report retains the classified reason.
    pub fn with_ph_solve_mode(mut self, mode: PhSolveMode) -> Self {
        self.ph_solve_mode = mode;
        self
    }

    /// Replaces outer scalar solver controls after validation.
    pub fn with_temperature_options(
        mut self,
        options: PhTemperatureSolveOptions,
    ) -> Result<Self, ReactionExtentError> {
        options.validate()?;
        self.temperature_options = options;
        Ok(self)
    }

    /// Replaces the physical composition used to seed the next `P,H` solve.
    ///
    /// This is the continuation boundary for a target-enthalpy sweep. The
    /// composition is validated against the immutable resolved layout before
    /// the request is returned, so a point from another system cannot be
    /// accidentally used as a seed.
    pub fn with_initial_composition(
        mut self,
        composition: MultiphaseInitialComposition,
    ) -> Result<Self, ReactionExtentError> {
        let layout = MultiphaseEquilibriumLayout::new(self.resolved.phase_specs().to_vec())?;
        composition.validate_for(&layout)?;
        self.initial_composition = composition;
        self.validate()?;
        Ok(self)
    }

    /// Replaces the target enthalpy and continuation temperature seed.
    ///
    /// The temperature is only an initial estimate; the accepted solution
    /// still has to satisfy the full `P,H` solve and its bounds. Keeping this
    /// operation on the typed request prevents a batch layer from mutating
    /// solver state between points.
    pub fn with_target_enthalpy_and_seed(
        mut self,
        target_enthalpy: crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::
            TotalEnthalpyJoules,
        initial_temperature: f64,
    ) -> Result<Self, ReactionExtentError> {
        self.constraint = EquilibriumConstraint::ph_joules(
            self.constraint.pressure(),
            self.constraint.reference_pressure(),
            target_enthalpy,
            initial_temperature,
        )?;
        self.validate()?;
        Ok(self)
    }

    /// Returns the immutable resolved system borrowed by this request.
    pub fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    /// Returns the physical composition used to initialize this request.
    ///
    /// Batch continuation uses this value only as the first-point seed; all
    /// later seeds come from the preceding accepted solution.
    pub fn initial_composition(&self) -> &MultiphaseInitialComposition {
        &self.initial_composition
    }

    /// Returns the `P,H` constraint.
    pub fn constraint(&self) -> EquilibriumConstraint {
        self.constraint
    }

    /// Selected P,H numerical formulation.
    pub fn ph_solve_mode(&self) -> PhSolveMode {
        self.ph_solve_mode
    }

    /// Creates a shallow clone of this request with a different `P,H` solve mode.
    ///
    /// The resolved system, initial composition, constraint, temperature bounds,
    /// enthalpy model, and solve options are shared via cloning the `Arc`-wrapped
    /// resolved system and copying the remaining value-type fields. Only the
    /// [`PhSolveMode`] is replaced, which lets the monolithic-to-nested fallback
    /// path retry with a different numerical route without rebuilding the request.
    fn cloned_with_ph_solve_mode(&self, mode: PhSolveMode) -> Self {
        Self {
            resolved: self.resolved.clone(),
            initial_composition: self.initial_composition.clone(),
            constraint: self.constraint,
            temperature_bounds: self.temperature_bounds,
            enthalpy: self.enthalpy.clone(),
            solve_options: self.solve_options.clone(),
            solve_mode: self.solve_mode.clone(),
            ph_solve_mode: mode,
            temperature_options: self.temperature_options.clone(),
            thermochemistry: self.thermochemistry.clone(),
        }
    }

    /// Borrows the inner fixed-`P,T` solver options so a frontend can attach
    /// cancellation, progress, or diagnostic policy without opening the
    /// request's mutable internals.
    pub fn solve_options(&self) -> &EquilibriumSolveOptions {
        &self.solve_options
    }

    /// Whether this request owns a fixed declared phase layout rather than a
    /// phase-control lifecycle. Prepared P,H range state is valid only for
    /// this topology-stable route.
    pub(crate) fn uses_fixed_declared_phases(&self) -> bool {
        matches!(
            self.solve_mode,
            PhaseEquilibriumSolveMode::FixedDeclaredPhases
        )
    }
}

/// Prepared fixed-active P,H state shared by a target-enthalpy range.
///
/// The structural problem and symbolic residual are immutable across points;
/// only the accepted composition/temperature seed and the two scalar energy
/// parameters are retargeted. This state is crate-private so the public range
/// facade remains the sole orchestration boundary.
pub(crate) struct PreparedPhContinuationState {
    prepared: crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
        PreparedFixedActiveProblem,
    formulation: PreparedPhFormulation,
    thermochemistry: ResolvedThermochemistry,
    rst_problem: Option<RstPreparedProblem>,
}

/// Shared fixed-`P,T` template for a nested `P,H` target range.
///
/// The nested route still owns an independent scalar bracket for every
/// target. Only the topology-stable inner formulation is shared here. In
/// particular, this type never carries a bracket, a rejected trial, or a
/// bounded phase-control active set from one target into the next.
pub(crate) struct PreparedNestedPhContinuationState {
    fixed_template: Rc<RefCell<Option<PreparedPhaseEquilibriumTemplate>>>,
}

impl PreparedNestedPhContinuationState {
    /// Creates an empty shared template for a fixed-declared nested range.
    pub(crate) fn new(
        request: &ResolvedPhaseEnthalpyRequest<'_>,
    ) -> Result<Self, ReactionExtentError> {
        request.validate()?;
        if !matches!(request.ph_solve_mode, PhSolveMode::NestedTemperature)
            || !request.uses_fixed_declared_phases()
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_range_nested_template",
                message: "nested P,H template requires fixed declared phases".into(),
            });
        }
        Ok(Self {
            fixed_template: Rc::new(RefCell::new(None)),
        })
    }

    /// Solves one scalar target while sharing only the prepared inner
    /// fixed-`P,T` formulation with earlier targets.
    pub(crate) fn solve(
        &self,
        request: ResolvedPhaseEnthalpyRequest<'_>,
    ) -> Result<FixedPressureEnthalpySolution, ReactionExtentError> {
        solve_resolved_ph_nested_with_template(request, Some(Rc::clone(&self.fixed_template)))
    }
}

impl PreparedPhContinuationState {
    /// Builds structural P,H state from the first range point.
    pub(crate) fn new(
        request: &ResolvedPhaseEnthalpyRequest<'_>,
    ) -> Result<Self, ReactionExtentError> {
        request.validate()?;
        if !matches!(request.ph_solve_mode, PhSolveMode::Monolithic)
            || !matches!(
                request.solve_mode,
                PhaseEquilibriumSolveMode::FixedDeclaredPhases
            )
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_range_prepared_state",
                message: "prepared monolithic P,H state requires fixed declared phases".into(),
            });
        }
        let thermochemistry =
            request
                .thermochemistry
                .clone()
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "ph_monolithic_thermochemistry",
                    message: "prepared monolithic P,H state requires resolved thermochemistry"
                        .into(),
                })?;
        let (target, initial_temperature) = validated_ph_parameters(request.constraint)?;
        let initial_molar_enthalpies = thermochemistry.evaluate_enthalpy(initial_temperature)?;
        let scale = EnthalpyScale::from_magnitudes(
            target,
            request.initial_composition.moles(),
            &initial_molar_enthalpies,
        )?;
        let conditions = request.constraint.conditions_at(initial_temperature)?;
        let bundle = build_phase_equilibrium_problem_with_timing(
            PhaseEquilibriumBuildRequest::new(
                &request.resolved,
                conditions,
                request.initial_composition.clone(),
                request.solve_options.trace_seed_policy(),
                SupportedPhaseModelPolicy::default(),
            )?,
            request.solve_options.timing_mode(),
        )?;
        let prepared = bundle.into_prepared_fixed_active_problem()?;
        let formulation = PreparedPhFormulation::new(
            prepared.prepared.clone(),
            thermochemistry.clone(),
            request.temperature_bounds,
            target,
            scale,
        )?;
        let rst_problem = request
            .solve_options
            .prepares_rst_backend()
            .then(|| prepare_rst_symbolic_ph_problem(&formulation))
            .transpose()?;
        Ok(Self {
            prepared,
            formulation,
            thermochemistry,
            rst_problem,
        })
    }

    /// Solves one accepted target using the shared prepared state.
    pub(crate) fn solve(
        &mut self,
        request: &ResolvedPhaseEnthalpyRequest<'_>,
        formulation_reused: bool,
    ) -> Result<FixedPressureEnthalpySolution, ReactionExtentError> {
        let started = Instant::now();
        request.validate()?;
        let (target, initial_temperature) = validated_ph_parameters(request.constraint)?;
        let initial_molar_enthalpies = self
            .thermochemistry
            .evaluate_enthalpy(initial_temperature)?;
        let scale = EnthalpyScale::from_magnitudes(
            target,
            request.initial_composition.moles(),
            &initial_molar_enthalpies,
        )?;
        let conditions = request.constraint.conditions_at(initial_temperature)?;
        let log_seed =
            LogMolesInitialGuess::from_initial_moles(request.initial_composition.moles())?;
        let gibbs = self
            .thermochemistry
            .gibbs_snapshot_for_legacy_boundary(initial_temperature)?;
        let prepared =
            self.prepared
                .prepared
                .retarget_with_gibbs(conditions, log_seed.clone(), gibbs)?;
        self.prepared.prepared = prepared.clone();
        self.formulation = self.formulation.retarget(prepared, target, scale)?;
        if let Some(rst_problem) = self.rst_problem.as_mut() {
            rst_problem.set_ph_parameters(target, scale.joules())?;
        }

        let options = request.temperature_options.clone();
        let acceptance_options = options.monolithic_options()?.acceptance();
        let runner = PreparedMonolithicPhRunner::new(
            self.formulation.clone(),
            request.solve_options.clone().into_settings(),
            acceptance_options,
        )?;
        if let Some(control) = options.execution_control() {
            control.check_cancelled()?;
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::FormulationPreparation,
                Some(0),
                Some(1),
                Some(initial_temperature),
            ));
        }
        let outcome = runner.solve_from_log_moles_and_temperature_seed(
            &log_seed,
            initial_temperature,
            self.rst_problem.as_ref(),
        )?;
        let accepted_conditions = request
            .constraint
            .conditions_at(outcome.snapshot.temperature)?;
        let accepted_solution = self.prepared.prepared.accepted_solution_at_conditions(
            outcome.snapshot.log_moles.clone(),
            outcome.pt_validation.clone(),
            accepted_conditions,
        )?;
        let accepted_gibbs = self
            .thermochemistry
            .evaluate_gibbs(outcome.snapshot.temperature)?
            .into_iter()
            .map(|value| Rc::new(move |_| value) as GibbsFn)
            .collect::<Vec<_>>();
        let accepted_report = self
            .prepared
            .report
            .at_conditions(accepted_conditions, &accepted_gibbs)?;
        let solution = MultiphaseEquilibriumSolution::from_fixed_active_parts(
            self.prepared.metadata.clone(),
            accepted_report,
            accepted_solution,
            outcome.solve_report.clone(),
            self.prepared.timing,
        )?;
        let timing_enabled = request.solve_options.timing_mode() == EquilibriumTimingMode::Enabled;
        let report = PhTemperatureSolveReport {
            solve_path: PhSolvePath::MonolithicFixedActiveSet,
            fallback_reason: None,
            target_enthalpy: target,
            enthalpy_scale_joules: scale.joules(),
            absolute_enthalpy_tolerance_joules: options.absolute_enthalpy_tolerance_joules,
            scaled_enthalpy_tolerance: options.scaled_enthalpy_tolerance,
            max_temperature_evaluations: options.max_temperature_evaluations,
            initial_temperature_seed: initial_temperature,
            monotonicity_policy: options.monotonicity_policy,
            max_inner_backend_attempts: options.max_inner_backend_attempts,
            max_inner_nonlinear_iterations: options.max_inner_nonlinear_iterations,
            max_phase_control_transitions: options.max_phase_control_transitions,
            max_wall_time: options.max_wall_time,
            inner_backend_attempts: outcome.solve_report.started_attempt_count(),
            inner_nonlinear_iterations: outcome.solve_report.nonlinear_iterations(),
            phase_control_transitions: 0,
            fixed_formulation_builds: usize::from(!formulation_reused),
            fixed_formulation_reuses: usize::from(formulation_reused),
            inner_timing: *solution.timing_report(),
            timing: PhTemperatureTimingReport {
                enabled: timing_enabled,
                total: if timing_enabled {
                    started.elapsed()
                } else {
                    Duration::ZERO
                },
                scalar_orchestration: Duration::ZERO,
                enthalpy_evaluation: Duration::ZERO,
            },
            iterations: 0,
            trials: Vec::new(),
            monolithic_evidence: Some(PhMonolithicEvidence {
                solve_report: solution.solve_report().clone(),
                multi_start_report: solution.multi_start_report().cloned(),
                phase_control_report: None,
                acceptance_report: solution.acceptance_report().cloned(),
                inner_timing: *solution.timing_report(),
            }),
        };
        Ok(FixedPressureEnthalpySolution {
            solution,
            calculated_enthalpy: outcome.snapshot.total_enthalpy,
            target_enthalpy: target,
            report,
            thermochemistry: Some(self.thermochemistry.clone()),
        })
    }
}

/// Solves one resolved system at fixed pressure and total enthalpy.
///
/// The selected [`PhSolveMode`] determines the numerical formulation. The
/// canonical resolved-data request defaults to the coupled monolithic path.
/// The nested path keeps every scalar trial local until its bracket accepts
/// one temperature; no failed trial is published. `Auto` attempts monolithic
/// solving first and falls back only for classified numerical failures.
pub fn solve_resolved_ph(
    request: ResolvedPhaseEnthalpyRequest<'_>,
) -> Result<FixedPressureEnthalpySolution, ReactionExtentError> {
    match request.ph_solve_mode {
        PhSolveMode::NestedTemperature => solve_resolved_ph_nested(request),
        PhSolveMode::Monolithic => match request.solve_mode {
            PhaseEquilibriumSolveMode::FixedDeclaredPhases => solve_resolved_ph_monolithic(request),
            PhaseEquilibriumSolveMode::BoundedPhaseControl(_) => {
                solve_resolved_ph_monolithic_phase_control(request)
            }
        },
        PhSolveMode::Auto => {
            let monolithic_request = request.cloned_with_ph_solve_mode(PhSolveMode::Monolithic);
            match solve_resolved_ph(monolithic_request) {
                Ok(solution) => Ok(solution),
                Err(error) if error.is_retryable_formulation_failure() => {
                    let fallback_reason = PhFallbackReason::from_error(&error);
                    match solve_resolved_ph(
                        request.cloned_with_ph_solve_mode(PhSolveMode::NestedTemperature),
                    ) {
                        Ok(mut nested) => {
                            nested.report.fallback_reason = Some(fallback_reason);
                            Ok(nested)
                        }
                        Err(nested) => Err(ReactionExtentError::PhAutoFallbackFailed {
                            monolithic: Box::new(error),
                            nested: Box::new(nested),
                        }),
                    }
                }
                Err(error) => Err(error),
            }
        }
    }
}

/// Runs monolithic P,H candidates through the canonical transactional phase
/// lifecycle.
///
/// The nonlinear formulation is rebuilt only for the active mask selected by
/// the phase controller. Activation, destruction, hysteresis, cycle checks,
/// and final publication remain owned by `PreparedPhaseControlRunner`; this
/// function only adapts one fixed active mask to the coupled P,H equations.
fn solve_resolved_ph_monolithic_phase_control(
    request: ResolvedPhaseEnthalpyRequest<'_>,
) -> Result<FixedPressureEnthalpySolution, ReactionExtentError> {
    let started = Instant::now();
    request.validate()?;
    let thermochemistry = request.thermochemistry.clone().ok_or_else(|| {
        ReactionExtentError::InvalidProblem {
            field: "ph_monolithic_thermochemistry",
            message: "monolithic P,H requires a ResolvedThermochemistry bundle with Gibbs, enthalpy, and Cp capabilities".to_string(),
        }
    })?;
    let (target, initial_temperature) = validated_ph_parameters(request.constraint)?;
    let initial_molar_enthalpies = thermochemistry.evaluate_enthalpy(initial_temperature)?;
    let scale = EnthalpyScale::from_magnitudes(
        target,
        request.initial_composition.moles(),
        &initial_molar_enthalpies,
    )?;
    let initial_conditions = request.constraint.conditions_at(initial_temperature)?;
    let timing_mode = request.solve_options.timing_mode();
    let bundle = build_phase_equilibrium_problem_with_timing(
        PhaseEquilibriumBuildRequest::new(
            &request.resolved,
            initial_conditions,
            request.initial_composition.clone(),
            request.solve_options.trace_seed_policy(),
            SupportedPhaseModelPolicy::default(),
        )?,
        timing_mode,
    )?;
    let policy = match request.solve_mode {
        PhaseEquilibriumSolveMode::BoundedPhaseControl(policy) => policy,
        PhaseEquilibriumSolveMode::FixedDeclaredPhases => {
            return Err(ReactionExtentError::InvalidProblem {
                field: "ph_monolithic_phase_control",
                message: "the bounded monolithic adapter requires a phase-control policy"
                    .to_string(),
            });
        }
    };
    let template = bundle.into_phase_control_template(|manager| {
        *manager = policy.into_phase_manager();
    })?;
    let settings = request.solve_options.clone().into_settings();
    let temperature_options = request.temperature_options.clone();
    let acceptance_options = temperature_options.monolithic_options()?.acceptance();
    let temperature_bounds = request.temperature_bounds;
    let mut temperature_seed = initial_temperature;
    let solution = template.solve_with_fixed_active_solver(
        settings,
        |runner, active, seed, species_phase, full_element_totals| {
            solve_monolithic_active_set_candidate(
                runner,
                active,
                seed,
                species_phase,
                full_element_totals,
                &thermochemistry,
                temperature_bounds,
                target,
                scale,
                &mut temperature_seed,
                &acceptance_options,
            )
        },
    )?;
    let phase_control_report = solution.phase_control_report().cloned().ok_or_else(|| {
        ReactionExtentError::InvalidCandidate {
            field: "ph_monolithic_phase_control_report",
            message: "bounded monolithic solve published no phase-control report".to_string(),
        }
    })?;
    let calculated_enthalpy = thermochemistry.enthalpy_model().evaluate_total(
        solution.component_moles(),
        solution.conditions().temperature(),
    )?;
    let inner_backend_attempts = phase_control_report
        .nonlinear_reports
        .iter()
        .map(|report| report.started_attempt_count())
        .sum();
    let inner_nonlinear_iterations = phase_control_report
        .nonlinear_reports
        .iter()
        .map(|report| report.nonlinear_iterations())
        .sum();
    let timing_enabled = timing_mode == EquilibriumTimingMode::Enabled;
    let total = started.elapsed();
    let report = PhTemperatureSolveReport {
        solve_path: PhSolvePath::MonolithicPhaseControl,
        fallback_reason: None,
        target_enthalpy: target,
        enthalpy_scale_joules: scale.joules(),
        absolute_enthalpy_tolerance_joules: temperature_options.absolute_enthalpy_tolerance_joules,
        scaled_enthalpy_tolerance: temperature_options.scaled_enthalpy_tolerance,
        max_temperature_evaluations: temperature_options.max_temperature_evaluations,
        initial_temperature_seed: initial_temperature,
        monotonicity_policy: temperature_options.monotonicity_policy,
        max_inner_backend_attempts: temperature_options.max_inner_backend_attempts,
        max_inner_nonlinear_iterations: temperature_options.max_inner_nonlinear_iterations,
        max_phase_control_transitions: temperature_options.max_phase_control_transitions,
        max_wall_time: temperature_options.max_wall_time,
        inner_backend_attempts,
        inner_nonlinear_iterations,
        phase_control_transitions: phase_control_report.transitions.len(),
        // These counters describe reusable scalar P,H formulations, not the
        // per-active-set candidates owned by phase control. Their lifecycle
        // evidence is retained in `PhaseControlledSolveReport` instead.
        fixed_formulation_builds: 0,
        fixed_formulation_reuses: 0,
        inner_timing: *solution.timing_report(),
        timing: PhTemperatureTimingReport {
            enabled: timing_enabled,
            total: timing_enabled.then_some(total).unwrap_or(Duration::ZERO),
            scalar_orchestration: Duration::ZERO,
            enthalpy_evaluation: Duration::ZERO,
        },
        iterations: phase_control_report.iterations,
        trials: Vec::new(),
        monolithic_evidence: Some(PhMonolithicEvidence {
            solve_report: solution.solve_report().clone(),
            multi_start_report: solution.multi_start_report().cloned(),
            phase_control_report: Some(phase_control_report),
            acceptance_report: solution.acceptance_report().cloned(),
            inner_timing: *solution.timing_report(),
        }),
    };
    let result = FixedPressureEnthalpySolution {
        solution,
        calculated_enthalpy,
        target_enthalpy: target,
        report,
        thermochemistry: Some(thermochemistry),
    };
    if let Some(control) = temperature_options.execution_control() {
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::PublicationCompleted,
            Some(1),
            Some(1),
            Some(result.temperature()),
        ));
        control.check_cancelled()?;
    }
    Ok(result)
}

/// Runs the coupled P,H equations for one declared, fixed active phase set.
///
/// This private helper is the fixed-phase branch of the public monolithic
/// facade. Bounded phase control uses the sibling adapter above and never
/// treats a changing active set as one smooth residual system.
fn solve_resolved_ph_monolithic(
    request: ResolvedPhaseEnthalpyRequest<'_>,
) -> Result<FixedPressureEnthalpySolution, ReactionExtentError> {
    let started = Instant::now();
    request.validate()?;
    let thermochemistry = request.thermochemistry.clone().ok_or_else(|| {
        ReactionExtentError::InvalidProblem {
            field: "ph_monolithic_thermochemistry",
            message: "monolithic P,H requires a ResolvedThermochemistry bundle with Gibbs, enthalpy, and Cp capabilities".to_string(),
        }
    })?;
    let (target, initial_temperature) = validated_ph_parameters(request.constraint)?;
    let initial_molar_enthalpies = thermochemistry.evaluate_enthalpy(initial_temperature)?;
    let scale = EnthalpyScale::from_magnitudes(
        target,
        request.initial_composition.moles(),
        &initial_molar_enthalpies,
    )?;
    let initial_conditions = request.constraint.conditions_at(initial_temperature)?;
    let timing_mode = request.solve_options.timing_mode();
    let bundle = build_phase_equilibrium_problem_with_timing(
        PhaseEquilibriumBuildRequest::new(
            &request.resolved,
            initial_conditions,
            request.initial_composition.clone(),
            request.solve_options.trace_seed_policy(),
            SupportedPhaseModelPolicy::default(),
        )?,
        timing_mode,
    )?;
    let prepared = bundle.into_prepared_fixed_active_problem()?;
    let formulation = PreparedPhFormulation::new(
        prepared.prepared.clone(),
        thermochemistry.clone(),
        request.temperature_bounds,
        target,
        scale,
    )?;
    let options = request.temperature_options.clone();
    let acceptance_options = options.monolithic_options()?.acceptance();
    let runner = PreparedMonolithicPhRunner::new(
        formulation,
        request.solve_options.clone().into_settings(),
        acceptance_options,
    )?;
    if let Some(control) = options.execution_control() {
        control.check_cancelled()?;
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::FormulationPreparation,
            Some(0),
            Some(1),
            Some(initial_temperature),
        ));
    }
    let outcome = runner.solve_from_temperature_seed(initial_temperature)?;
    let conditions = request
        .constraint
        .conditions_at(outcome.snapshot.temperature)?;
    let accepted_solution = prepared.prepared.accepted_solution_at_conditions(
        outcome.snapshot.log_moles.clone(),
        outcome.pt_validation.clone(),
        conditions,
    )?;
    // The structural bridge was prepared at the temperature seed, while the
    // coupled solve owns the accepted temperature. Rebuild only the immutable
    // report view at that accepted condition before publication; otherwise the
    // result gate would correctly reject a solution/report temperature
    // mismatch even though the numerical candidate itself is valid.
    let accepted_gibbs = thermochemistry
        .evaluate_gibbs(outcome.snapshot.temperature)?
        .into_iter()
        .map(|value| Rc::new(move |_| value) as GibbsFn)
        .collect::<Vec<_>>();
    let accepted_report = prepared.report.at_conditions(conditions, &accepted_gibbs)?;
    let solution = MultiphaseEquilibriumSolution::from_fixed_active_parts(
        prepared.metadata,
        accepted_report,
        accepted_solution,
        outcome.solve_report.clone(),
        prepared.timing,
    )?;
    let timing_enabled = timing_mode == EquilibriumTimingMode::Enabled;
    let total = started.elapsed();
    let report = PhTemperatureSolveReport {
        solve_path: PhSolvePath::MonolithicFixedActiveSet,
        fallback_reason: None,
        target_enthalpy: target,
        enthalpy_scale_joules: scale.joules(),
        absolute_enthalpy_tolerance_joules: options.absolute_enthalpy_tolerance_joules,
        scaled_enthalpy_tolerance: options.scaled_enthalpy_tolerance,
        max_temperature_evaluations: options.max_temperature_evaluations,
        initial_temperature_seed: initial_temperature,
        monotonicity_policy: options.monotonicity_policy,
        max_inner_backend_attempts: options.max_inner_backend_attempts,
        max_inner_nonlinear_iterations: options.max_inner_nonlinear_iterations,
        max_phase_control_transitions: options.max_phase_control_transitions,
        max_wall_time: options.max_wall_time,
        inner_backend_attempts: outcome.solve_report.started_attempt_count(),
        inner_nonlinear_iterations: outcome.solve_report.nonlinear_iterations(),
        phase_control_transitions: 0,
        fixed_formulation_builds: 1,
        fixed_formulation_reuses: 0,
        inner_timing: *solution.timing_report(),
        timing: PhTemperatureTimingReport {
            enabled: timing_enabled,
            total: timing_enabled.then_some(total).unwrap_or(Duration::ZERO),
            scalar_orchestration: Duration::ZERO,
            enthalpy_evaluation: Duration::ZERO,
        },
        iterations: 0,
        trials: Vec::new(),
        monolithic_evidence: Some(PhMonolithicEvidence {
            solve_report: solution.solve_report().clone(),
            multi_start_report: solution.multi_start_report().cloned(),
            phase_control_report: None,
            acceptance_report: solution.acceptance_report().cloned(),
            inner_timing: *solution.timing_report(),
        }),
    };
    let result = FixedPressureEnthalpySolution {
        solution,
        calculated_enthalpy: outcome.snapshot.total_enthalpy,
        target_enthalpy: target,
        report,
        thermochemistry: Some(thermochemistry),
    };
    if let Some(control) = options.execution_control() {
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::PublicationCompleted,
            Some(1),
            Some(1),
            Some(result.temperature()),
        ));
        control.check_cancelled()?;
    }
    Ok(result)
}

/// Retained safeguarded reference implementation.
fn solve_resolved_ph_nested(
    request: ResolvedPhaseEnthalpyRequest<'_>,
) -> Result<FixedPressureEnthalpySolution, ReactionExtentError> {
    solve_resolved_ph_nested_with_template(request, None)
}

fn solve_resolved_ph_nested_with_template(
    request: ResolvedPhaseEnthalpyRequest<'_>,
    shared_template: Option<Rc<RefCell<Option<PreparedPhaseEquilibriumTemplate>>>>,
) -> Result<FixedPressureEnthalpySolution, ReactionExtentError> {
    let outer_started = Instant::now();
    request.validate()?;
    let constraint = request.constraint;
    let (target, initial_temperature) = validated_ph_parameters(constraint)?;
    let initial_molar = request.enthalpy.evaluate_molar(initial_temperature)?;
    let scale = EnthalpyScale::from_magnitudes(
        target,
        request.initial_composition.moles(),
        &initial_molar,
    )?;

    let bounds = request.temperature_bounds;
    let options = request.temperature_options;
    let resolved = &request.resolved;
    let initial_composition = request.initial_composition.clone();
    let enthalpy = request.enthalpy.clone();
    let solve_options = request.solve_options.clone();
    let solve_mode = request.solve_mode.clone();
    let thermochemistry = request.thermochemistry.clone();
    let timing_enabled = solve_options.timing_mode() == EquilibriumTimingMode::Enabled;
    let execution_control = options.execution_control().cloned();
    // The evaluator owns one clone; the outer publication boundary keeps a
    // second handle so cancellation is still checked after the bracket has
    // returned but before the immutable result is exposed.
    let publication_execution_control = execution_control.clone();
    let evaluator_options = options.clone();
    if let Some(control) = &execution_control {
        control.check_cancelled()?;
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::FormulationPreparation,
            None,
            None,
            None,
        ));
    }
    // A fixed declared phase layout has no active-set history. Its structural
    // preparation can therefore be safely retained while the P,H bracket
    // samples temperatures in non-monotone order. Bounded phase control is
    // deliberately excluded: its active set has lifecycle semantics and must
    // not inherit state from an arbitrary rejected scalar trial.
    let fixed_template = shared_template.or_else(|| {
        matches!(&solve_mode, PhaseEquilibriumSolveMode::FixedDeclaredPhases)
            .then(|| Rc::new(RefCell::new(None)))
    });
    let fixed_template_for_report = fixed_template.clone();
    let template_was_initialized = fixed_template_for_report
        .as_ref()
        .is_some_and(|template| template.borrow().is_some());
    // Continuation is safe for fixed declared phases. Bounded phase control
    // keeps the original inventory seed because its active set can change.
    let use_continuation = matches!(&solve_mode, PhaseEquilibriumSolveMode::FixedDeclaredPhases);
    // A bounded active set may change between two bracket samples. Secant
    // interpolation would then act as if their enthalpy values came from one
    // smooth branch, so only the topology-safe bisection step is allowed.
    let allow_interpolation = use_continuation;
    let mut continuation_seed: Option<LogMolesInitialGuess> = None;
    let mut temperature_evaluations = 0usize;
    let mut inner_backend_attempts_total = 0usize;
    let mut inner_nonlinear_iterations_total = 0usize;
    let mut phase_control_transitions_total = 0usize;
    let nested_evidence = Rc::new(RefCell::new(Vec::<PhNestedTrialEvidence>::new()));
    let nested_evidence_for_trial = Rc::clone(&nested_evidence);
    let enthalpy_timing = Rc::new(RefCell::new(Duration::ZERO));
    let enthalpy_timing_for_trial = Rc::clone(&enthalpy_timing);
    let evaluator = move |temperature: f64| {
        if temperature_evaluations >= evaluator_options.max_temperature_evaluations {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_budget",
                message: format!(
                    "outer P,H solve exceeded the global temperature evaluation budget of {}",
                    evaluator_options.max_temperature_evaluations
                ),
            });
        }
        if let Some(control) = &execution_control {
            control.check_cancelled()?;
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::TemperatureTrialStarted,
                Some(temperature_evaluations),
                Some(evaluator_options.max_temperature_evaluations),
                Some(temperature),
            ));
        }
        temperature_evaluations += 1;
        if let Some(control) = &execution_control {
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::InnerSolveStarted,
                Some(temperature_evaluations - 1),
                Some(evaluator_options.max_temperature_evaluations),
                Some(temperature),
            ));
            control.check_cancelled()?;
        }
        let trial_evaluator = PhTrialEvaluator {
            resolved,
            initial_composition: &initial_composition,
            constraint,
            enthalpy: &enthalpy,
            solve_options: &solve_options,
            solve_mode: &solve_mode,
            fixed_template: fixed_template.as_deref(),
            execution_control: execution_control.as_ref(),
            timing_enabled,
        };
        let trial_continuation_seed = use_continuation.then(|| continuation_seed.take()).flatten();
        let outcome = match trial_evaluator.evaluate(temperature, trial_continuation_seed) {
            Ok(outcome) => outcome,
            Err(error) => {
                if !matches!(&error, ReactionExtentError::Cancelled) {
                    if let Some(control) = &execution_control {
                        control.report(EquilibriumProgressEvent::new(
                            EquilibriumProgressStage::TemperatureTrialRejected,
                            Some(temperature_evaluations - 1),
                            Some(evaluator_options.max_temperature_evaluations),
                            Some(temperature),
                        ));
                    }
                }
                return Err(error);
            }
        };
        if let Some(control) = &execution_control {
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::InnerSolveCompleted,
                Some(temperature_evaluations - 1),
                Some(evaluator_options.max_temperature_evaluations),
                Some(temperature),
            ));
            for _ in 0..outcome.evidence.phase_control_transitions {
                control.report(EquilibriumProgressEvent::new(
                    EquilibriumProgressStage::PhaseTransitionAccepted,
                    Some(temperature_evaluations - 1),
                    Some(evaluator_options.max_temperature_evaluations),
                    Some(temperature),
                ));
            }
            control.check_cancelled()?;
        }
        let attempts = outcome.evidence.backend_attempts;
        inner_backend_attempts_total = inner_backend_attempts_total.saturating_add(attempts);
        if let Some(limit) = evaluator_options.max_inner_backend_attempts {
            if inner_backend_attempts_total > limit {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "inner_backend_budget",
                    message: format!(
                        "outer P,H solve exceeded the global inner backend-attempt budget of {limit}"
                    ),
                });
            }
        }
        let nonlinear_iterations = outcome.evidence.nonlinear_iterations;
        inner_nonlinear_iterations_total =
            inner_nonlinear_iterations_total.saturating_add(nonlinear_iterations);
        if let Some(limit) = evaluator_options.max_inner_nonlinear_iterations {
            if inner_nonlinear_iterations_total > limit {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "inner_nonlinear_iteration_budget",
                    message: format!(
                        "outer P,H solve exceeded the global inner nonlinear-iteration budget of {limit}"
                    ),
                });
            }
        }
        phase_control_transitions_total = phase_control_transitions_total
            .saturating_add(outcome.evidence.phase_control_transitions);
        if let Some(limit) = evaluator_options.max_phase_control_transitions {
            if phase_control_transitions_total > limit {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "phase_control_transition_budget",
                    message: format!(
                        "outer P,H solve exceeded the global phase-transition budget of {limit}"
                    ),
                });
            }
        }
        if use_continuation {
            continuation_seed = Some(LogMolesInitialGuess::from_initial_moles(
                outcome.solution.component_moles(),
            )?);
        }
        *enthalpy_timing_for_trial.borrow_mut() += outcome.enthalpy_evaluation;
        nested_evidence_for_trial
            .borrow_mut()
            .push(outcome.evidence.clone());
        if let Some(control) = &execution_control {
            control.check_cancelled()?;
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::TemperatureTrialAccepted,
                Some(temperature_evaluations - 1),
                Some(evaluator_options.max_temperature_evaluations),
                Some(temperature),
            ));
        }
        Ok((outcome.solution, outcome.total_enthalpy))
    };

    let report_options = options.clone();
    let scalar_started = Instant::now();
    let (solution, total_enthalpy, iterations, mut trials) = solve_bracketed_temperature_from(
        bounds,
        scale,
        target,
        Some(initial_temperature),
        allow_interpolation,
        options,
        outer_started,
        evaluator,
    )?;
    if let Some(control) = &publication_execution_control {
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::PublicationStarted,
            None,
            None,
            Some(solution.conditions().temperature()),
        ));
        control.check_cancelled()?;
    }
    let timing = PhTemperatureTimingReport {
        enabled: timing_enabled,
        total: if timing_enabled {
            outer_started.elapsed()
        } else {
            Duration::ZERO
        },
        scalar_orchestration: if timing_enabled {
            scalar_started.elapsed()
        } else {
            Duration::ZERO
        },
        enthalpy_evaluation: if timing_enabled {
            *enthalpy_timing.borrow()
        } else {
            Duration::ZERO
        },
    };
    let (
        inner_backend_attempts,
        inner_nonlinear_iterations,
        phase_control_transitions,
        inner_timing,
    ) = attach_nested_evidence(&mut trials, &nested_evidence.borrow())?;
    let template_is_initialized = fixed_template_for_report
        .as_ref()
        .is_some_and(|template| template.borrow().is_some());
    let fixed_formulation_builds =
        usize::from(!template_was_initialized && template_is_initialized);
    let fixed_formulation_reuses = if !template_is_initialized {
        0
    } else if template_was_initialized {
        trials.len()
    } else {
        trials.len().saturating_sub(1)
    };
    let report = PhTemperatureSolveReport {
        solve_path: PhSolvePath::NestedTemperature,
        fallback_reason: None,
        target_enthalpy: target,
        enthalpy_scale_joules: scale.joules(),
        absolute_enthalpy_tolerance_joules: report_options.absolute_enthalpy_tolerance_joules,
        scaled_enthalpy_tolerance: report_options.scaled_enthalpy_tolerance,
        max_temperature_evaluations: report_options.max_temperature_evaluations,
        initial_temperature_seed: initial_temperature,
        monotonicity_policy: report_options.monotonicity_policy,
        max_inner_backend_attempts: report_options.max_inner_backend_attempts,
        max_inner_nonlinear_iterations: report_options.max_inner_nonlinear_iterations,
        max_phase_control_transitions: report_options.max_phase_control_transitions,
        max_wall_time: report_options.max_wall_time,
        inner_backend_attempts,
        inner_nonlinear_iterations,
        phase_control_transitions,
        fixed_formulation_builds,
        fixed_formulation_reuses,
        inner_timing,
        timing,
        iterations,
        trials,
        monolithic_evidence: None,
    };
    let result = FixedPressureEnthalpySolution {
        solution,
        calculated_enthalpy: total_enthalpy,
        target_enthalpy: target,
        report,
        thermochemistry,
    };
    if let Some(control) = &publication_execution_control {
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::PublicationCompleted,
            None,
            None,
            Some(result.temperature()),
        ));
        control.check_cancelled()?;
    }
    Ok(result)
}

#[cfg(test)]
fn solve_bracketed_temperature<T, F>(
    bounds: TemperatureBounds,
    scale: EnthalpyScale,
    target: f64,
    options: PhTemperatureSolveOptions,
    evaluate: F,
) -> Result<(T, f64, usize, Vec<PhTemperatureTrial>), ReactionExtentError>
where
    F: FnMut(f64) -> Result<(T, f64), ReactionExtentError>,
{
    solve_bracketed_temperature_from(
        bounds,
        scale,
        target,
        None,
        true,
        options,
        Instant::now(),
        evaluate,
    )
}

fn solve_bracketed_temperature_from<T, F>(
    bounds: TemperatureBounds,
    scale: EnthalpyScale,
    target: f64,
    seed_temperature: Option<f64>,
    allow_interpolation: bool,
    options: PhTemperatureSolveOptions,
    started: Instant,
    evaluate: F,
) -> Result<(T, f64, usize, Vec<PhTemperatureTrial>), ReactionExtentError>
where
    F: FnMut(f64) -> Result<(T, f64), ReactionExtentError>,
{
    options.validate()?;
    let route_options = options.nested_options()?;
    let acceptance = route_options.acceptance();
    let nested_options = NestedBracketOptions {
        scaled_enthalpy_tolerance: acceptance.scaled_enthalpy_tolerance(),
        absolute_enthalpy_tolerance_joules: acceptance.absolute_enthalpy_tolerance_joules(),
        temperature_tolerance: route_options.temperature_tolerance(),
        max_iterations: route_options.max_iterations(),
        max_temperature_evaluations: route_options.max_temperature_evaluations(),
        allow_interpolation,
        allow_bracketed_sign_search: matches!(
            route_options.monotonicity_policy(),
            PhMonotonicityPolicy::AllowBracketedSignSearch
        ),
        max_wall_time: route_options.max_wall_time(),
    };
    let result = solve_nested_bracketed_temperature(
        bounds,
        scale,
        target,
        seed_temperature,
        nested_options,
        started,
        evaluate,
    )?;
    let trials = result
        .trials
        .into_iter()
        .map(PhTemperatureTrial::from)
        .collect();
    Ok((
        result.value,
        result.total_enthalpy,
        result.iterations,
        trials,
    ))
}
/// Publishes inner execution evidence onto the corresponding outer trials.
///
/// The bracket helper and the evaluator are intentionally separate layers, so
/// their record counts are checked at this boundary. A silent `zip` would
/// make an incomplete report look valid and could hide a publication bug.
#[derive(Debug, Clone)]
struct PhNestedTrialEvidence {
    backend_attempts: usize,
    nonlinear_iterations: usize,
    phase_control_transitions: usize,
    preparation: PhTrialPreparation,
    phase_states: Vec<PhTrialPhaseState>,
    timing: EquilibriumTimingReport,
    trial_timing: PhTrialTimingReport,
    inner_evidence: Option<PhTrialInnerEvidence>,
}

/// Immutable outcome of one local `P,T -> H` trial. It is deliberately not
/// published until the outer scalar workflow accepts its temperature.
struct PhTrialOutcome {
    solution: MultiphaseEquilibriumSolution,
    total_enthalpy: f64,
    enthalpy_evaluation: Duration,
    evidence: PhNestedTrialEvidence,
}

/// Builds one P,H scalar trial from immutable request data.
///
/// This object owns no global counters, continuation state, or publication
/// state. The outer solve is therefore the sole authority for budgets,
/// progress notifications, and transactional result publication.
struct PhTrialEvaluator<'request, 'enthalpy> {
    resolved: &'request ResolvedPhaseSystem,
    initial_composition: &'request MultiphaseInitialComposition,
    constraint: EquilibriumConstraint,
    enthalpy: &'request EnthalpyModel<'enthalpy>,
    solve_options: &'request EquilibriumSolveOptions,
    solve_mode: &'request PhaseEquilibriumSolveMode,
    fixed_template: Option<&'request RefCell<Option<PreparedPhaseEquilibriumTemplate>>>,
    execution_control: Option<&'request EquilibriumExecutionControl>,
    timing_enabled: bool,
}

impl<'request, 'enthalpy> PhTrialEvaluator<'request, 'enthalpy> {
    /// Evaluates `H(solution_of_P_T(T), T)` without mutating the request or
    /// publishing a partial result.
    ///
    /// Solves the inner fixed-`P,T` equilibrium at the supplied temperature
    /// using the configured solve mode (fixed-declared-phases or bounded phase
    /// control), then evaluates the additive total enthalpy from the accepted
    /// composition. The optional continuation seed is forwarded to the inner
    /// solver when the caller wants to warm-start from a previous point.
    /// Timing is collected when the request has timing enabled.
    fn evaluate(
        &self,
        temperature: f64,
        continuation_seed: Option<LogMolesInitialGuess>,
    ) -> Result<PhTrialOutcome, ReactionExtentError> {
        let trial_started = self.timing_enabled.then(Instant::now);
        let conditions = self.constraint.conditions_at(temperature)?;
        let trial_solve_options = if let Some(control) = self.execution_control {
            self.solve_options
                .clone()
                .with_execution_control(control.clone())
        } else {
            self.solve_options.clone()
        };
        let (solution, preparation) = match self.solve_mode {
            PhaseEquilibriumSolveMode::FixedDeclaredPhases => {
                let mut seeds = vec![LogMolesInitialGuess::from_initial_moles(
                    self.initial_composition.moles(),
                )?];
                if let Some(seed) = continuation_seed {
                    seeds.push(seed);
                }
                let template =
                    self.fixed_template
                        .ok_or_else(|| ReactionExtentError::InvalidProblem {
                            field: "ph_fixed_template",
                            message: "fixed P,H trial is missing its prepared formulation"
                                .to_string(),
                        })?;
                let mut template = template.borrow_mut();
                let reused = template.is_some();
                if template.is_none() {
                    let build_request = PhaseEquilibriumBuildRequest::new(
                        self.resolved,
                        conditions,
                        self.initial_composition.clone(),
                        self.solve_options.trace_seed_policy(),
                        SupportedPhaseModelPolicy::default(),
                    )?;
                    let bundle = build_phase_equilibrium_problem_with_timing(
                        build_request,
                        self.solve_options.timing_mode(),
                    )?;
                    *template = Some(
                        bundle.into_temperature_template(
                            self.solve_options.prepares_rst_backend(),
                            self.solve_options.timing_mode(),
                        )?,
                    );
                }
                let solution = template
                    .as_mut()
                    .ok_or_else(|| ReactionExtentError::InvalidProblem {
                        field: "ph_fixed_template",
                        message: "fixed P,H formulation was not initialized".to_string(),
                    })?
                    .solve_at_with_initial_guesses(
                        conditions,
                        seeds,
                        trial_solve_options.into_settings(),
                        self.solve_options.timing_mode(),
                    )?;
                (
                    solution,
                    if reused {
                        PhTrialPreparation::FixedFormulationReused
                    } else {
                        PhTrialPreparation::FixedFormulationInitial
                    },
                )
            }
            PhaseEquilibriumSolveMode::BoundedPhaseControl(policy) => {
                let inner = ResolvedPhaseEquilibriumRequest::new(
                    self.resolved,
                    conditions,
                    self.initial_composition.clone(),
                )
                .with_solve_options(trial_solve_options)
                .with_phase_control_policy(policy.clone());
                (
                    solve_resolved_pt(inner)?,
                    PhTrialPreparation::BoundedPhaseControlIsolated,
                )
            }
        };
        let enthalpy_started = self.timing_enabled.then(Instant::now);
        let total_enthalpy = self
            .enthalpy
            .evaluate_total(solution.component_moles(), temperature)?;
        let enthalpy_evaluation =
            enthalpy_started.map_or(Duration::ZERO, |started| started.elapsed());
        let trial_timing = PhTrialTimingReport {
            enabled: self.timing_enabled,
            total: trial_started.map_or(Duration::ZERO, |started| started.elapsed()),
            inner_equilibrium: if self.timing_enabled {
                solution.timing_report().total()
            } else {
                Duration::ZERO
            },
            enthalpy_evaluation,
        };
        let evidence = PhNestedTrialEvidence {
            backend_attempts: solution.started_backend_attempts(),
            nonlinear_iterations: solution.nonlinear_iterations(),
            phase_control_transitions: solution.phase_control_transitions(),
            preparation,
            phase_states: phase_states_from_solution(&solution)?,
            timing: *solution.timing_report(),
            trial_timing,
            inner_evidence: Some(PhTrialInnerEvidence {
                solve_report: solution.solve_report().clone(),
                multi_start_report: solution.multi_start_report().cloned(),
                phase_control_report: solution.phase_control_report().cloned(),
                acceptance_report: solution.acceptance_report().cloned(),
            }),
        };
        Ok(PhTrialOutcome {
            solution,
            total_enthalpy,
            enthalpy_evaluation,
            evidence,
        })
    }
}

fn attach_nested_evidence(
    trials: &mut [PhTemperatureTrial],
    evidence: &[PhNestedTrialEvidence],
) -> Result<(usize, usize, usize, EquilibriumTimingReport), ReactionExtentError> {
    if trials.len() != evidence.len() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "ph_nested_evidence",
            message: format!(
                "outer trial count {} does not match nested evidence count {}",
                trials.len(),
                evidence.len()
            ),
        });
    }

    let mut inner_backend_attempts = 0usize;
    let mut inner_nonlinear_iterations = 0usize;
    let mut phase_control_transitions = 0usize;
    let mut inner_timing = EquilibriumTimingReport::default();
    for (trial, evidence) in trials.iter_mut().zip(evidence) {
        trial.inner_backend_attempts = evidence.backend_attempts;
        trial.inner_nonlinear_iterations = evidence.nonlinear_iterations;
        trial.phase_control_transitions = evidence.phase_control_transitions;
        trial.preparation = evidence.preparation;
        trial.phase_states = evidence.phase_states.clone();
        trial.inner_timing = evidence.timing;
        trial.timing = evidence.trial_timing;
        trial.inner_evidence = evidence.inner_evidence.clone();
        inner_backend_attempts = inner_backend_attempts.saturating_add(evidence.backend_attempts);
        inner_nonlinear_iterations =
            inner_nonlinear_iterations.saturating_add(evidence.nonlinear_iterations);
        phase_control_transitions =
            phase_control_transitions.saturating_add(evidence.phase_control_transitions);
        inner_timing.accumulate(evidence.timing);
    }
    Ok((
        inner_backend_attempts,
        inner_nonlinear_iterations,
        phase_control_transitions,
        inner_timing,
    ))
}

/// Captures the phase lifecycle states that identify an accepted P,T branch.
fn phase_states_from_solution(
    solution: &MultiphaseEquilibriumSolution,
) -> Result<Vec<PhTrialPhaseState>, ReactionExtentError> {
    solution
        .phases()
        .iter()
        .map(|descriptor| {
            let phase = descriptor.id().clone();
            let status = solution.phase_status(&phase).ok_or_else(|| {
                ReactionExtentError::InvalidProblem {
                    field: "phase_status",
                    message: format!(
                        "accepted equilibrium solution omitted the status for phase '{phase:?}'"
                    ),
                }
            })?;
            Ok(PhTrialPhaseState { phase, status })
        })
        .collect()
}

/// Extracts the scalar parameters required by the P,H workflow without
/// relying on an `expect` after validation. This keeps the public solve path
/// typed-error-only even if a future request mutation violates its invariant.
fn validated_ph_parameters(
    constraint: EquilibriumConstraint,
) -> Result<(f64, f64), ReactionExtentError> {
    match constraint {
        EquilibriumConstraint::PH {
            target_enthalpy,
            initial_temperature,
            ..
        } => Ok((target_enthalpy.joules(), initial_temperature)),
        EquilibriumConstraint::PT { .. } => Err(ReactionExtentError::InvalidProblem {
            field: "constraint",
            message: "resolved phase enthalpy workflow requires a PH constraint".to_string(),
        }),
    }
}

#[cfg(test)]
#[allow(deprecated)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::sync::{Arc, Mutex};

    use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::MultiphaseEquilibriumLayout;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::EquilibriumSolveReport;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus;
    use crate::Thermodynamics::User_PhaseOrSolution::PhaseSpec;
    use crate::Thermodynamics::User_substances::{DataType, SubsData};

    fn linear_enthalpy_model() -> EnthalpyModel<'static> {
        EnthalpyModel::from_functions(vec![Arc::new(|temperature| Ok(temperature * 10.0))]).unwrap()
    }

    fn synthetic_inner_evidence() -> PhTrialInnerEvidence {
        use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
        use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
            SolverAttemptOutcome, SolverAttemptReport, SolverBackend, SolverPolicy,
        };

        let backend = SolverBackend::Legacy(Solvers::NR);
        PhTrialInnerEvidence {
            solve_report: EquilibriumSolveReport {
                policy: SolverPolicy::Single(backend),
                attempts: vec![SolverAttemptReport {
                    backend,
                    outcome: SolverAttemptOutcome::Accepted,
                    metrics: None,
                }],
                accepted_backend: backend,
            },
            multi_start_report: None,
            phase_control_report: None,
            acceptance_report: None,
        }
    }

    #[test]
    fn monolithic_mode_requires_resolved_thermochemistry_capabilities() {
        assert_eq!(PhSolveMode::default(), PhSolveMode::Monolithic);

        let mut data = SubsData::new();
        data.set_substances(vec!["A".to_string()]);
        let phase = PhaseSpec::ideal_gas(PhaseId::new(None), vec!["A".to_string()]).unwrap();
        let resolved =
            ResolvedPhaseSystem::new(vec![phase.clone()], HashMap::from([(None, data)])).unwrap();
        let layout = MultiphaseEquilibriumLayout::new(vec![phase]).unwrap();
        let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1.0]).unwrap();
        let request = ResolvedPhaseEnthalpyRequest::new(
            &resolved,
            composition,
            EquilibriumConstraint::ph(101_325.0, 101_325.0, 0.0, 700.0).unwrap(),
            TemperatureBounds::new(300.0, 1_000.0).unwrap(),
            EnthalpyModel::from_functions(vec![Arc::new(|_| Ok(0.0))]).unwrap(),
        )
        .unwrap()
        .with_ph_solve_mode(PhSolveMode::Monolithic);

        assert_eq!(request.ph_solve_mode(), PhSolveMode::Monolithic);
        let stages = Arc::new(Mutex::new(Vec::new()));
        let stages_for_sink = Arc::clone(&stages);
        let control = EquilibriumExecutionControl::new().with_progress_sink(move |event| {
            stages_for_sink.lock().unwrap().push(event.stage());
        });
        let auto_request = request
            .cloned_with_ph_solve_mode(PhSolveMode::Auto)
            .with_temperature_options(
                PhTemperatureSolveOptions::default().with_execution_control(control),
            )
            .unwrap();
        assert!(matches!(
            solve_resolved_ph(request),
            Err(ReactionExtentError::InvalidProblem {
                field: "ph_monolithic_thermochemistry",
                ..
            })
        ));
        assert!(matches!(
            solve_resolved_ph(auto_request),
            Err(ReactionExtentError::InvalidProblem {
                field: "ph_monolithic_thermochemistry",
                ..
            })
        ));
        assert!(!stages
            .lock()
            .unwrap()
            .contains(&EquilibriumProgressStage::TemperatureTrialStarted));
    }

    #[test]
    fn auto_fallback_reason_preserves_retryable_backend_classification() {
        let error = ReactionExtentError::BackendFailure {
            backend: "monolithic".to_string(),
            kind: crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::
                BackendFailureKind::NumericalBreakdown,
            message: "synthetic breakdown".to_string(),
        };
        let reason = PhFallbackReason::from_error(&error);

        assert_eq!(reason.error_kind(), ReactionExtentErrorKind::BackendFailure);
        assert!(reason.message().contains("monolithic"));
        assert!(error.is_retryable_backend_failure());
    }

    #[test]
    fn validated_ph_parameter_extraction_rejects_pt_without_panicking() {
        let constraint = EquilibriumConstraint::pt(
            EquilibriumConditions::new(700.0, 101_325.0, 101_325.0).unwrap(),
        );

        let error = validated_ph_parameters(constraint).unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "constraint",
                ..
            }
        ));
    }

    #[test]
    fn additive_enthalpy_derivative_reports_only_the_explicit_temperature_term() {
        let model = EnthalpyModel::from_functions_with_heat_capacity(
            vec![Arc::new(|temperature| Ok(100.0 + 2.0 * temperature))],
            vec![Some(Arc::new(|_| Ok(2.0)))],
        )
        .unwrap();

        let evaluation = model
            .evaluate_total_with_partial_temperature_derivative(&[3.0], 400.0)
            .unwrap();

        assert_eq!(evaluation.molar_enthalpies(), &[900.0]);
        assert_eq!(evaluation.heat_capacities(), &[2.0]);
        assert_eq!(evaluation.total_enthalpy(), 2_700.0);
        assert_eq!(evaluation.partial_temperature_derivative(), 6.0);
    }

    #[test]
    fn additive_enthalpy_derivative_rejects_missing_heat_capacity() {
        let model = linear_enthalpy_model();
        let error = model
            .evaluate_total_with_partial_temperature_derivative(&[1.0], 400.0)
            .unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "heat_capacity",
                ..
            }
        ));
    }

    #[test]
    fn sampled_non_monotone_branch_is_rejected_by_default_but_explicit_compatibility_allows_it() {
        let branch = |temperature: f64| {
            if temperature <= 850.0 {
                -1.0 - 99.0 * (temperature - 500.0) / 350.0
            } else {
                -100.0 + 101.0 * (temperature - 850.0) / 350.0
            }
        };
        let bounds = TemperatureBounds::new(500.0, 1_200.0).unwrap();
        let scale = EnthalpyScale::new(1.0).unwrap();

        let error = solve_bracketed_temperature(
            bounds,
            scale,
            0.0,
            PhTemperatureSolveOptions::default(),
            |temperature| Ok((temperature, branch(temperature))),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "enthalpy_non_monotone",
                ..
            }
        ));

        let compatibility = PhTemperatureSolveOptions::default()
            .with_monotonicity_policy(PhMonotonicityPolicy::AllowBracketedSignSearch);
        let (temperature, enthalpy, _, _) =
            solve_bracketed_temperature(bounds, scale, 0.0, compatibility, |temperature| {
                Ok((temperature, branch(temperature)))
            })
            .unwrap();
        assert!((temperature - (850.0 + 100.0 * 350.0 / 101.0)).abs() < 1.0e-5);
        assert!(enthalpy.abs() < 1.0e-6);
    }

    #[test]
    fn ph_temperature_seed_is_a_real_trial_and_can_accept_the_root() {
        let (temperature, enthalpy, iterations, trials) = solve_bracketed_temperature_from(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            Some(750.0),
            true,
            PhTemperatureSolveOptions::default(),
            Instant::now(),
            |temperature| Ok((temperature, temperature - 750.0)),
        )
        .unwrap();

        assert_eq!(iterations, 0);
        assert_eq!(temperature, 750.0);
        assert_eq!(enthalpy, 0.0);
        assert_eq!(trials.len(), 3);
        assert_eq!(trials[2].temperature(), 750.0);
    }

    #[test]
    fn ph_seed_reports_observed_multiple_brackets_instead_of_selecting_one() {
        let error = solve_bracketed_temperature_from(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            Some(750.0),
            true,
            PhTemperatureSolveOptions::default(),
            Instant::now(),
            |temperature| Ok((temperature, (temperature - 650.0) * (temperature - 850.0))),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "enthalpy_multiple_brackets",
                ..
            }
        ));
    }

    #[test]
    fn safeguarded_scalar_solver_uses_an_interior_secant_step_for_a_linear_branch() {
        let (temperature, enthalpy, iterations, trials) = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            PhTemperatureSolveOptions::default(),
            |temperature| Ok((temperature, temperature - 775.0)),
        )
        .unwrap();

        assert_eq!(iterations, 1);
        assert_eq!(temperature, 775.0);
        assert_eq!(enthalpy, 0.0);
        assert_eq!(trials[2].step_kind(), PhTemperatureStepKind::Interpolation);
    }

    #[test]
    fn phase_control_scalar_path_uses_bisection_without_branch_interpolation() {
        let (temperature, enthalpy, iterations, trials) = solve_bracketed_temperature_from(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            None,
            false,
            PhTemperatureSolveOptions::default(),
            Instant::now(),
            |temperature| Ok((temperature, temperature - 750.0)),
        )
        .unwrap();

        assert_eq!(iterations, 1);
        assert_eq!(temperature, 750.0);
        assert_eq!(enthalpy, 0.0);
        assert_eq!(trials[2].step_kind(), PhTemperatureStepKind::Bisection);
    }

    #[test]
    fn interpolation_step_falls_back_when_secant_hugs_a_bracket_endpoint() {
        assert_eq!(
            safeguarded_interpolation_step(0.0, -1.0, 100.0, 100.0),
            None
        );
        assert_eq!(
            safeguarded_interpolation_step(0.0, -1.0, 100.0, 1.0),
            Some(50.0)
        );
    }

    #[test]
    fn derivative_ready_enthalpy_model_rejects_misaligned_capabilities() {
        let error = match EnthalpyModel::from_functions_with_heat_capacity(
            vec![Arc::new(|temperature| Ok(temperature))],
            Vec::new(),
        ) {
            Ok(_) => panic!("misaligned heat-capacity capabilities must be rejected"),
            Err(error) => error,
        };

        assert!(matches!(
            error,
            ReactionExtentError::DimensionMismatch(message)
                if message.contains("heat-capacity functions")
        ));
    }

    #[test]
    fn bracketed_solver_recovers_known_temperature_without_a_chemical_solver() {
        let model = linear_enthalpy_model();
        let target = 7_500.0;
        let scale = EnthalpyScale::from_magnitudes(target, &[1.0], &[7_000.0]).unwrap();
        let options = PhTemperatureSolveOptions {
            scaled_enthalpy_tolerance: 1.0e-10,
            absolute_enthalpy_tolerance_joules: 1.0e-6,
            temperature_tolerance: 1.0e-10,
            max_iterations: 80,
            max_temperature_evaluations: 82,
            monotonicity_policy: PhMonotonicityPolicy::default(),
            max_inner_backend_attempts: None,
            max_inner_nonlinear_iterations: None,
            max_phase_control_transitions: None,
            max_wall_time: None,
            execution_control: None,
        };
        let (value, enthalpy, _, _) = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            scale,
            target,
            options,
            |temperature| Ok((temperature, model.evaluate_total(&[1.0], temperature)?)),
        )
        .unwrap();

        assert!((value - 750.0).abs() < 1.0e-8);
        assert!((enthalpy - target).abs() < 1.0e-7);
    }

    #[test]
    fn nonreacting_sensible_heat_fixture_recovers_the_analytic_root() {
        // H(T) = 2 * (100 + 4T) + 3 * (-50 + 2T) = 50 + 14T.
        // This deliberately contains no chemistry, so it isolates the outer
        // energy equation from all fixed-P,T equilibrium behavior.
        let model = EnthalpyModel::from_functions(vec![
            Arc::new(|temperature| Ok(100.0 + 4.0 * temperature)),
            Arc::new(|temperature| Ok(-50.0 + 2.0 * temperature)),
        ])
        .unwrap();
        let moles = [2.0, 3.0];
        let target_temperature = 725.0;
        let target = model.evaluate_total(&moles, target_temperature).unwrap();
        let initial_molar = model.evaluate_molar(500.0).unwrap();
        let scale = EnthalpyScale::from_magnitudes(target, &moles, &initial_molar).unwrap();

        let (temperature, enthalpy, _, _) = solve_bracketed_temperature(
            TemperatureBounds::new(300.0, 1_000.0).unwrap(),
            scale,
            target,
            PhTemperatureSolveOptions::default(),
            |temperature| Ok((temperature, model.evaluate_total(&moles, temperature)?)),
        )
        .unwrap();

        assert!((temperature - target_temperature).abs() < 1.0e-5);
        assert!((enthalpy - target).abs() < 1.0e-4);
    }

    #[test]
    fn inventory_scaling_preserves_temperature_and_mole_fractions() {
        let model = EnthalpyModel::from_functions(vec![
            Arc::new(|temperature| Ok(300.0 + 3.0 * temperature)),
            Arc::new(|temperature| Ok(-100.0 + 5.0 * temperature)),
        ])
        .unwrap();
        let base_moles = [1.5, 2.5];
        let scaled_moles = [1_500.0, 2_500.0];
        let target_temperature = 810.0;
        let base_target = model
            .evaluate_total(&base_moles, target_temperature)
            .unwrap();
        let scaled_target = model
            .evaluate_total(&scaled_moles, target_temperature)
            .unwrap();
        let base_scale = EnthalpyScale::from_magnitudes(
            base_target,
            &base_moles,
            &model.evaluate_molar(400.0).unwrap(),
        )
        .unwrap();
        let scaled_scale = EnthalpyScale::from_magnitudes(
            scaled_target,
            &scaled_moles,
            &model.evaluate_molar(400.0).unwrap(),
        )
        .unwrap();
        let options = PhTemperatureSolveOptions::default();

        let (base_temperature, _, _, _) = solve_bracketed_temperature(
            TemperatureBounds::new(300.0, 1_000.0).unwrap(),
            base_scale,
            base_target,
            options.clone(),
            |temperature| Ok((temperature, model.evaluate_total(&base_moles, temperature)?)),
        )
        .unwrap();
        let (scaled_temperature, _, _, _) = solve_bracketed_temperature(
            TemperatureBounds::new(300.0, 1_000.0).unwrap(),
            scaled_scale,
            scaled_target,
            options,
            |temperature| {
                Ok((
                    temperature,
                    model.evaluate_total(&scaled_moles, temperature)?,
                ))
            },
        )
        .unwrap();

        assert!((base_temperature - target_temperature).abs() < 1.0e-5);
        assert!((scaled_temperature - base_temperature).abs() < 1.0e-5);
        let base_fraction = base_moles[0] / base_moles.iter().sum::<f64>();
        let scaled_fraction = scaled_moles[0] / scaled_moles.iter().sum::<f64>();
        assert!((base_fraction - scaled_fraction).abs() < 1.0e-15);
    }

    #[test]
    fn ph_constraint_rejects_nonfinite_target_before_solving() {
        let error = EquilibriumConstraint::ph(101_325.0, 101_325.0, f64::NAN, 700.0).unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "target_enthalpy",
                ..
            }
        ));
    }

    #[test]
    fn acceptance_contract_keeps_an_absolute_floor_for_zero_target() {
        let scale = EnthalpyScale::new(1.0e-12).unwrap();
        let options = PhTemperatureSolveOptions {
            scaled_enthalpy_tolerance: 1.0e-12,
            absolute_enthalpy_tolerance_joules: 1.0e-4,
            temperature_tolerance: 1.0e-8,
            max_iterations: 4,
            max_temperature_evaluations: 6,
            monotonicity_policy: PhMonotonicityPolicy::default(),
            max_inner_backend_attempts: None,
            max_inner_nonlinear_iterations: None,
            max_phase_control_transitions: None,
            max_wall_time: None,
            execution_control: None,
        };

        assert_eq!(options.accepted_enthalpy_error_limit_joules(scale), 1.0e-4);
        assert!(options.accepts_enthalpy_error(5.0e-5, scale));
        assert!(!options.accepts_enthalpy_error(2.0e-4, scale));
    }

    #[test]
    fn inner_backend_budget_is_typed_and_rejects_zero() {
        let invalid = PhTemperatureSolveOptions {
            max_inner_backend_attempts: Some(0),
            ..PhTemperatureSolveOptions::default()
        };
        assert!(matches!(
            invalid.validate(),
            Err(ReactionExtentError::InvalidProblem {
                field: "max_inner_backend_attempts",
                ..
            })
        ));

        let limited = PhTemperatureSolveOptions::default()
            .with_max_inner_backend_attempts(3)
            .unwrap();
        assert_eq!(limited.max_inner_backend_attempts, Some(3));
        assert_eq!(
            limited
                .clone()
                .without_inner_backend_attempt_limit()
                .max_inner_backend_attempts,
            None
        );
        assert!(limited.with_max_inner_backend_attempts(0).is_err());

        let iteration_limited = PhTemperatureSolveOptions::default()
            .with_max_inner_nonlinear_iterations(7)
            .unwrap();
        assert_eq!(iteration_limited.max_inner_nonlinear_iterations, Some(7));
        assert_eq!(
            iteration_limited
                .clone()
                .without_inner_nonlinear_iteration_limit()
                .max_inner_nonlinear_iterations,
            None
        );
        assert!(iteration_limited
            .with_max_inner_nonlinear_iterations(0)
            .is_err());

        let transition_limited = PhTemperatureSolveOptions::default()
            .with_max_phase_control_transitions(4)
            .unwrap();
        assert_eq!(transition_limited.max_phase_control_transitions, Some(4));
        assert_eq!(
            transition_limited
                .clone()
                .without_phase_control_transition_limit()
                .max_phase_control_transitions,
            None
        );
        assert!(transition_limited
            .with_max_phase_control_transitions(0)
            .is_err());
    }

    #[test]
    fn bracket_solver_does_not_accept_temperature_width_without_energy_accuracy() {
        let scale = EnthalpyScale::new(1.0).unwrap();
        let options = PhTemperatureSolveOptions {
            scaled_enthalpy_tolerance: 1.0e-12,
            absolute_enthalpy_tolerance_joules: 1.0e-12,
            temperature_tolerance: 10.0,
            max_iterations: 16,
            max_temperature_evaluations: 18,
            monotonicity_policy: PhMonotonicityPolicy::default(),
            max_inner_backend_attempts: None,
            max_inner_nonlinear_iterations: None,
            max_phase_control_transitions: None,
            max_wall_time: None,
            execution_control: None,
        };
        let error = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            scale,
            0.0,
            options,
            |temperature| Ok(((), if temperature < 750.0 { -1.0 } else { 1.0 })),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "enthalpy_acceptance",
                ..
            }
        ));
    }

    #[test]
    fn outer_temperature_budget_is_global_across_bracket_and_midpoints() {
        let scale = EnthalpyScale::new(1.0).unwrap();
        let options = PhTemperatureSolveOptions {
            scaled_enthalpy_tolerance: 1.0e-12,
            absolute_enthalpy_tolerance_joules: 1.0e-12,
            temperature_tolerance: 1.0e-12,
            max_iterations: 80,
            max_temperature_evaluations: 2,
            monotonicity_policy: PhMonotonicityPolicy::default(),
            max_inner_backend_attempts: None,
            max_inner_nonlinear_iterations: None,
            max_phase_control_transitions: None,
            max_wall_time: None,
            execution_control: None,
        };
        let evaluations = std::cell::Cell::new(0usize);
        let error = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            scale,
            0.0,
            options,
            |temperature| {
                evaluations.set(evaluations.get() + 1);
                Ok((temperature, temperature - 800.0))
            },
        )
        .unwrap_err();

        assert_eq!(evaluations.get(), 2);
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "temperature_budget",
                ..
            }
        ));
    }

    #[test]
    fn nested_evidence_is_attached_to_trials_without_silent_truncation() {
        let mut trials = vec![PhTemperatureTrial {
            temperature: 700.0,
            step_kind: PhTemperatureStepKind::Bisection,
            total_enthalpy: 10.0,
            enthalpy_error_joules: 0.0,
            scaled_error: 0.0,
            inner_backend_attempts: 0,
            inner_nonlinear_iterations: 0,
            phase_control_transitions: 0,
            preparation: PhTrialPreparation::Unspecified,
            phase_states: Vec::new(),
            inner_timing: EquilibriumTimingReport::default(),
            timing: PhTrialTimingReport::default(),
            inner_evidence: None,
        }];
        let evidence = vec![PhNestedTrialEvidence {
            backend_attempts: 3,
            nonlinear_iterations: 11,
            phase_control_transitions: 2,
            preparation: PhTrialPreparation::BoundedPhaseControlIsolated,
            phase_states: vec![PhTrialPhaseState {
                phase: PhaseId::new(Some("condensed".to_string())),
                status: PhaseStatus::Appeared,
            }],
            timing: EquilibriumTimingReport::default(),
            trial_timing: PhTrialTimingReport {
                enabled: true,
                total: Duration::from_millis(3),
                inner_equilibrium: Duration::from_millis(2),
                enthalpy_evaluation: Duration::from_micros(10),
            },
            inner_evidence: Some(synthetic_inner_evidence()),
        }];

        let totals = attach_nested_evidence(&mut trials, &evidence).unwrap();
        assert_eq!(totals.0, 3);
        assert_eq!(totals.1, 11);
        assert_eq!(totals.2, 2);
        assert_eq!(trials[0].inner_backend_attempts(), 3);
        assert_eq!(trials[0].inner_nonlinear_iterations(), 11);
        assert_eq!(trials[0].phase_control_transitions(), 2);
        assert_eq!(
            trials[0].preparation(),
            PhTrialPreparation::BoundedPhaseControlIsolated
        );
        assert_eq!(trials[0].phase_states().len(), 1);
        assert_eq!(trials[0].phase_states()[0].status(), PhaseStatus::Appeared);
        assert!(trials[0].timing().enabled());
        assert_eq!(trials[0].timing().total(), Duration::from_millis(3));
        assert_eq!(
            trials[0].timing().inner_equilibrium(),
            Duration::from_millis(2)
        );
        let inner = trials[0]
            .inner_evidence()
            .expect("accepted nested trial must retain its backend trace");
        assert_eq!(inner.solve_report().started_attempt_count(), 1);
        assert!(inner.multi_start_report().is_none());
        assert!(inner.phase_control_report().is_none());

        let mismatch = attach_nested_evidence(&mut trials, &[]).unwrap_err();
        assert!(matches!(
            mismatch,
            ReactionExtentError::InvalidProblem {
                field: "ph_nested_evidence",
                ..
            }
        ));
    }

    #[test]
    fn wall_time_budget_is_validated_and_enforced_by_the_bracket_helper() {
        let invalid = PhTemperatureSolveOptions {
            max_wall_time: Some(Duration::ZERO),
            ..PhTemperatureSolveOptions::default()
        };
        assert!(matches!(
            invalid.validate(),
            Err(ReactionExtentError::InvalidProblem {
                field: "max_wall_time",
                ..
            })
        ));

        let options = PhTemperatureSolveOptions {
            max_wall_time: Some(Duration::from_nanos(1)),
            ..PhTemperatureSolveOptions::default()
        };
        let error = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            options,
            |temperature| Ok((temperature, temperature - 800.0)),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "wall_time_budget",
                ..
            }
        ));
    }

    #[test]
    fn wall_time_budget_rejects_a_trial_that_finishes_after_the_deadline() {
        let options = PhTemperatureSolveOptions {
            max_wall_time: Some(Duration::from_millis(1)),
            ..PhTemperatureSolveOptions::default()
        };
        let error = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            options,
            |temperature| {
                std::thread::sleep(Duration::from_millis(5));
                Ok((temperature, temperature - 800.0))
            },
        )
        .unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "wall_time_budget",
                ..
            }
        ));
    }

    #[test]
    fn bracket_solver_reports_endpoint_temperature_on_inner_failure() {
        let error = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            PhTemperatureSolveOptions::default(),
            |_| {
                Err::<((), f64), _>(ReactionExtentError::InvalidProblem {
                    field: "inner",
                    message: "synthetic endpoint failure".to_string(),
                })
            },
        )
        .unwrap_err();

        match error {
            ReactionExtentError::TemperatureTrialFailed {
                trial_index,
                temperature,
                cause,
            } => {
                assert_eq!(trial_index, 0);
                assert_eq!(temperature, 500.0);
                assert!(matches!(
                    *cause,
                    ReactionExtentError::InvalidProblem { field: "inner", .. }
                ));
            }
            other => panic!("unexpected endpoint error: {other}"),
        }
    }

    #[test]
    fn bracket_solver_reports_midpoint_temperature_on_inner_failure() {
        let error = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            PhTemperatureSolveOptions::default(),
            |temperature| {
                if temperature == 750.0 {
                    Err(ReactionExtentError::InvalidProblem {
                        field: "inner",
                        message: "synthetic midpoint failure".to_string(),
                    })
                } else if temperature < 750.0 {
                    // The secant point would hug the lower endpoint, so the
                    // safeguarded hybrid must deliberately fall back to 750 K.
                    Ok((temperature, -1.0))
                } else {
                    Ok((temperature, 100.0))
                }
            },
        )
        .unwrap_err();

        match error {
            ReactionExtentError::TemperatureTrialFailed {
                trial_index,
                temperature,
                cause,
            } => {
                assert_eq!(trial_index, 2);
                assert_eq!(temperature, 750.0);
                assert!(matches!(
                    *cause,
                    ReactionExtentError::InvalidProblem { field: "inner", .. }
                ));
            }
            other => panic!("unexpected midpoint error: {other}"),
        }
    }

    #[test]
    fn bracket_solver_keeps_inner_cancellation_as_a_top_level_cancellation() {
        let error = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            EnthalpyScale::new(1.0).unwrap(),
            0.0,
            PhTemperatureSolveOptions::default(),
            |_| Err::<((), f64), _>(ReactionExtentError::Cancelled),
        )
        .unwrap_err();

        assert!(matches!(error, ReactionExtentError::Cancelled));
    }

    #[test]
    fn cancelled_ph_request_stops_before_an_inner_equilibrium_attempt() {
        let mut data = SubsData::new();
        data.set_substances(vec!["A".to_string()]);
        let phase = PhaseSpec::ideal_gas(PhaseId::new(None), vec!["A".to_string()]).unwrap();
        let resolved =
            ResolvedPhaseSystem::new(vec![phase], HashMap::from([(None, data)])).unwrap();
        let layout = MultiphaseEquilibriumLayout::new(vec![PhaseSpec::ideal_gas(
            PhaseId::new(None),
            vec!["A".to_string()],
        )
        .unwrap()])
        .unwrap();
        let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1.0]).unwrap();
        let constraint = EquilibriumConstraint::ph(101_325.0, 101_325.0, 0.0, 700.0).unwrap();
        let control = EquilibriumExecutionControl::new();
        control.request_cancel();
        let options = PhTemperatureSolveOptions::default().with_execution_control(control);
        let request = ResolvedPhaseEnthalpyRequest::new(
            &resolved,
            composition,
            constraint,
            TemperatureBounds::new(300.0, 1_000.0).unwrap(),
            EnthalpyModel::from_functions(vec![Arc::new(|_| Ok(0.0))]).unwrap(),
        )
        .unwrap()
        .with_temperature_options(options)
        .unwrap();
        // The request owns its resolved snapshot; the caller's source object
        // may be released before the expensive solve begins.
        drop(resolved);

        assert!(matches!(
            solve_resolved_ph(request),
            Err(ReactionExtentError::Cancelled)
        ));
    }

    #[test]
    fn ph_cancellation_after_inner_start_never_publishes_a_result() {
        let mut data = SubsData::new();
        data.set_substances(vec!["A".to_string()]);
        let phase = PhaseSpec::ideal_gas(PhaseId::new(None), vec!["A".to_string()]).unwrap();
        let resolved =
            ResolvedPhaseSystem::new(vec![phase], HashMap::from([(None, data)])).unwrap();
        let layout = MultiphaseEquilibriumLayout::new(vec![PhaseSpec::ideal_gas(
            PhaseId::new(None),
            vec!["A".to_string()],
        )
        .unwrap()])
        .unwrap();
        let composition = MultiphaseInitialComposition::from_dense(&layout, vec![1.0]).unwrap();
        let constraint = EquilibriumConstraint::ph(101_325.0, 101_325.0, 0.0, 700.0).unwrap();
        let stages = Arc::new(Mutex::new(Vec::new()));
        let stages_for_sink = Arc::clone(&stages);
        let control = EquilibriumExecutionControl::new();
        let cancellation = control.clone();
        let control = control.with_progress_sink(move |event| {
            stages_for_sink
                .lock()
                .expect("progress stage lock is not poisoned")
                .push(event.stage());
            if event.stage() == EquilibriumProgressStage::InnerSolveStarted {
                cancellation.request_cancel();
            }
        });
        let options = PhTemperatureSolveOptions::default().with_execution_control(control);
        let request = ResolvedPhaseEnthalpyRequest::new(
            &resolved,
            composition,
            constraint,
            TemperatureBounds::new(300.0, 1_000.0).unwrap(),
            EnthalpyModel::from_functions(vec![Arc::new(|_| Ok(0.0))]).unwrap(),
        )
        .unwrap()
        .with_temperature_options(options)
        .unwrap();

        assert!(matches!(
            solve_resolved_ph(request),
            Err(ReactionExtentError::Cancelled)
        ));
        let stages = stages.lock().unwrap().clone();
        assert!(stages.contains(&EquilibriumProgressStage::InnerSolveStarted));
        assert!(!stages.contains(&EquilibriumProgressStage::PublicationCompleted));
    }

    #[test]
    fn bracketed_solver_rejects_unreachable_target() {
        let model = linear_enthalpy_model();
        let scale = EnthalpyScale::new(10_000.0).unwrap();
        let error = solve_bracketed_temperature(
            TemperatureBounds::new(500.0, 1_000.0).unwrap(),
            scale,
            20_000.0,
            PhTemperatureSolveOptions::default(),
            |temperature| Ok((temperature, model.evaluate_total(&[1.0], temperature)?)),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "enthalpy_bracket",
                ..
            }
        ));
    }

    #[test]
    fn enthalpy_model_rejects_wrong_component_count() {
        let error = match EnthalpyModel::from_functions(Vec::new()) {
            Ok(_) => panic!("an empty enthalpy model must be rejected"),
            Err(error) => error,
        };
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "enthalpy_functions",
                ..
            }
        ));
    }

    #[test]
    fn resolved_bridge_does_not_trust_unprovenanced_dh_cache() {
        let mut data = SubsData::new();
        data.set_substances(vec!["A".to_string()]);
        data.therm_map_of_fun.insert(
            "A".to_string(),
            HashMap::from([(
                DataType::dH_fun,
                Some(Box::new(|temperature: f64| temperature * 4.0)
                    as Box<dyn Fn(f64) -> f64 + Send + Sync>),
            )]),
        );
        let phase = PhaseSpec::ideal_gas(PhaseId::new(None), vec!["A".to_string()]).unwrap();
        let resolved =
            ResolvedPhaseSystem::new(vec![phase], HashMap::from([(None, data)])).unwrap();

        let error = match EnthalpyModel::from_resolved_system(&resolved) {
            Ok(_) => panic!("an unresolved cache must not bypass thermochemistry provenance"),
            Err(error) => error,
        };
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "thermochemistry_provenance",
                ..
            }
        ));
        assert!(resolved
            .phase_data()
            .get(&None)
            .unwrap()
            .get_thermo_function("A", DataType::dH_fun)
            .is_some());
    }

    fn bundle_provenance(label: &str) -> ThermochemistryProvenance {
        ThermochemistryProvenance::new(
            PhaseComponentId::new(PhaseId::new(None), label),
            "synthetic",
            label,
            "gas",
        )
    }

    #[test]
    fn thermochemistry_bundle_keeps_order_and_optional_cp_explicit() {
        let bounds = TemperatureBounds::new(300.0, 1_000.0).unwrap();
        let bundle = ResolvedThermochemistry::from_functions(
            vec![bundle_provenance("A"), bundle_provenance("B")],
            bounds,
            vec![
                Arc::new(|temperature| Ok(temperature + 1.0)),
                Arc::new(|temperature| Ok(temperature + 2.0)),
            ],
            vec![
                Arc::new(|temperature| Ok(2.0 * temperature)),
                Arc::new(|temperature| Ok(3.0 * temperature)),
            ],
            vec![Some(Arc::new(|_| Ok(10.0))), None],
        )
        .unwrap();

        assert_eq!(
            bundle
                .provenance()
                .iter()
                .map(|row| row.component().label())
                .collect::<Vec<_>>(),
            vec!["A", "B"]
        );
        assert_eq!(bundle.evaluate_gibbs(500.0).unwrap(), vec![501.0, 502.0]);
        assert_eq!(
            bundle.evaluate_enthalpy(500.0).unwrap(),
            vec![1_000.0, 1_500.0]
        );
        assert_eq!(
            bundle.evaluate_heat_capacity(500.0).unwrap(),
            vec![Some(10.0), None]
        );
    }

    #[test]
    fn thermochemistry_bundle_rejects_dimension_mismatch_and_domain_escape() {
        let error = match ResolvedThermochemistry::from_functions(
            vec![bundle_provenance("A")],
            TemperatureBounds::new(300.0, 1_000.0).unwrap(),
            vec![Arc::new(|_| Ok(0.0)), Arc::new(|_| Ok(0.0))],
            vec![Arc::new(|_| Ok(0.0))],
            vec![None],
        ) {
            Ok(_) => panic!("a dimension-mismatched bundle must be rejected"),
            Err(error) => error,
        };
        assert!(matches!(error, ReactionExtentError::DimensionMismatch(_)));

        let bundle = ResolvedThermochemistry::from_functions(
            vec![bundle_provenance("A")],
            TemperatureBounds::new(300.0, 1_000.0).unwrap(),
            vec![Arc::new(|_| Ok(0.0))],
            vec![Arc::new(|_| Ok(0.0))],
            vec![None],
        )
        .unwrap();
        let error = bundle.evaluate_enthalpy(1_001.0).unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "temperature",
                ..
            }
        ));
    }
}
