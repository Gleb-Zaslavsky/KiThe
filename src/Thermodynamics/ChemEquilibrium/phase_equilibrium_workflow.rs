//! Narrow public fixed-`P,T` phase-equilibrium workflow.
//!
//! This facade owns orchestration only: it joins validated resolved data,
//! physical inventory, numerical settings, bridge construction, and immutable
//! result publication. It does not duplicate residual construction or expose
//! the historical mutable solver as an alternative public engine.

use serde::{Deserialize, Serialize};
use std::fmt;
use std::sync::Arc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_candidate_selection::{
    EquilibriumCandidatePhasePlan, EquilibriumCandidateSelectionReport,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::EquilibriumDiagnosticsOptions;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::{
    EquilibriumExecutionControl, EquilibriumProgressEvent, EquilibriumProgressStage,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumSolverSettings, Solvers,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, LogMolesInitialGuess, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverBackend, SolverCascadeBudget, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::{
    TemperatureGrid, TemperatureRangeRequest, TemperatureRangeSolution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::{
    EquilibriumTimingMode, EquilibriumTimingReport, EquilibriumTimingStage,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::InitialPhaseSet;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseManager;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::PhaseEquilibriumBuildRequest;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::SupportedPhaseModelPolicy;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::User_PhaseOrSolution::{
    ResolvedPhaseSystem, ResolvedPhaseSystemReport, SubstanceSystemFactory,
    SubstanceSystemFactoryError, SubstanceSystemSpec,
};
use crate::Thermodynamics::phase_layout::PhaseComponentId;
use crate::Thermodynamics::thermo_lib_api::ThermoRepository;

/// Current solve mode exposed by the typed public facade.
///
/// Fixed declared phases and bounded phase control share the same resolved
/// data bridge and immutable result model. The mutable solver is only the
/// retained compatibility implementation; bounded mode uses the prepared
/// immutable runner and callers do not receive an alternative mutable
/// production API.
#[derive(Debug, Clone)]
pub enum PhaseEquilibriumSolveMode {
    /// Solve all declared phases as one immutable fixed active set.
    FixedDeclaredPhases,
    /// Run the bounded active-set phase-control algorithm with an explicit
    /// hysteresis and initial-phase policy.
    BoundedPhaseControl(PhaseControlPolicy),
}

impl Default for PhaseEquilibriumSolveMode {
    fn default() -> Self {
        Self::FixedDeclaredPhases
    }
}

impl PhaseEquilibriumSolveMode {
    /// Explicit safe default for the canonical production path.
    pub fn fixed_declared_phases() -> Self {
        Self::FixedDeclaredPhases
    }

    /// Explicit constructor for bounded phase-control mode.
    pub fn bounded_phase_control(policy: PhaseControlPolicy) -> Self {
        Self::BoundedPhaseControl(policy)
    }
}

/// Narrow validated wrapper for solver policy knobs exposed at the facade boundary.
///
/// The inner `EquilibriumSolverSettings` still exists as the canonical backend
/// structure, but production callers should interact with this wrapper so the
/// public workflow can evolve without exposing every internal field.
#[derive(Clone)]
pub struct EquilibriumSolveOptions {
    settings: EquilibriumSolverSettings,
    timing_mode: EquilibriumTimingMode,
    diagnostics: EquilibriumDiagnosticsOptions,
    execution_control: Option<EquilibriumExecutionControl>,
}

/// Serializable, read-only record of effective numerical controls.
///
/// This deliberately records the resolved backend order rather than only the
/// caller's optional `SolverPolicy`: an absent explicit policy still expands
/// to the production cascade at solve time. It excludes the live execution
/// handle because cancellation/progress callbacks are process-local, not part
/// of a reproducible numerical contract.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EquilibriumSolveOptionsSnapshot {
    pub preferred_legacy_backend: String,
    pub effective_backend_order: Vec<String>,
    pub max_iterations: usize,
    pub tolerance: f64,
    pub scaling_enabled: bool,
    pub trace_seed_policy: String,
    pub continuation_seed_policy: String,
    pub equilibrium_constant_validation_mode: String,
    pub timing_mode: String,
    pub diagnostics_mode: String,
    pub diagnostics_max_events: usize,
    pub diagnostics_range_policy: String,
    pub solver_budget: Option<EquilibriumSolverBudgetSnapshot>,
    pub execution_control_attached: bool,
}

/// Serializable resource limits attached to a solver-cascade snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumSolverBudgetSnapshot {
    pub max_attempts: usize,
    pub max_iterations_per_attempt: usize,
    pub max_total_iterations: usize,
}

impl fmt::Debug for EquilibriumSolveOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EquilibriumSolveOptions")
            .finish_non_exhaustive()
    }
}

impl Default for EquilibriumSolveOptions {
    fn default() -> Self {
        Self {
            settings: EquilibriumSolverSettings::default(),
            timing_mode: EquilibriumTimingMode::default(),
            diagnostics: EquilibriumDiagnosticsOptions::default(),
            execution_control: None,
        }
    }
}

impl EquilibriumSolveOptions {
    /// Creates one wrapper around validated backend solver settings.
    pub fn new() -> Self {
        Self::default()
    }

    /// Wraps an existing backend settings value after validating it.
    #[cfg(test)]
    pub(crate) fn from_settings(
        settings: EquilibriumSolverSettings,
    ) -> Result<Self, ReactionExtentError> {
        settings.validate()?;
        Ok(Self {
            settings,
            timing_mode: EquilibriumTimingMode::default(),
            diagnostics: EquilibriumDiagnosticsOptions::default(),
            execution_control: None,
        })
    }

    /// Consumes the wrapper for the crate-internal immutable runner.
    pub(crate) fn into_settings(self) -> EquilibriumSolverSettings {
        let mut settings = self.settings;
        settings.execution_control = self.execution_control;
        settings
    }

    /// Installs cooperative cancellation and progress reporting for this
    /// transaction. The handle is shared with the worker-facing caller.
    pub fn with_execution_control(mut self, control: EquilibriumExecutionControl) -> Self {
        self.execution_control = Some(control);
        self
    }

    /// Returns the execution handle, if this request is externally controlled.
    pub fn execution_control(&self) -> Option<&EquilibriumExecutionControl> {
        self.execution_control.as_ref()
    }

    /// Uses the standard RST-first production cascade with the configured
    /// legacy backend family retained as a fallback.
    pub fn with_production_cascade(mut self) -> Self {
        self.settings.solver_policy = Some(SolverPolicy::production_default(self.settings.solver));
        self
    }

    /// Changes the nonlinear iteration budget without exposing raw settings.
    pub fn with_max_iterations(mut self, max_iter: usize) -> Result<Self, ReactionExtentError> {
        self.settings.solver_params.max_iter = max_iter;
        self.settings.validate()?;
        Ok(self)
    }

    /// Installs explicit resource limits for the ordered backend cascade.
    ///
    /// The GUI and other typed frontends use this boundary instead of reaching
    /// into the historical mutable settings object. Leaving it unset keeps the
    /// policy-derived production budget intact.
    pub fn with_solver_budget(
        mut self,
        budget: SolverCascadeBudget,
    ) -> Result<Self, ReactionExtentError> {
        self.settings.solver_budget = Some(budget);
        self.settings.validate()?;
        Ok(self)
    }

    /// Changes the common residual tolerance without exposing raw settings.
    pub fn with_tolerance(mut self, tolerance: f64) -> Result<Self, ReactionExtentError> {
        self.settings.solver_params.tol = tolerance;
        self.settings.validate()?;
        Ok(self)
    }

    /// Enables or disables the canonical residual/Jacobian scaling contract.
    pub fn with_scaling(mut self, enabled: bool) -> Self {
        self.settings.scaling_flag = enabled;
        self
    }

    /// Replaces the preferred backend ordering used by the solver cascade.
    pub fn with_solver_backend(mut self, solver: Solvers) -> Self {
        self.settings.solver = solver;
        // `solver_policy` is an explicit ordered cascade and therefore wins
        // over the legacy preferred-solver field. Clear it here so this
        // convenience method cannot silently become a no-op.
        self.settings.solver_policy = None;
        self
    }

    /// Installs an explicit ordered backend policy after validating it.
    pub fn with_solver_policy(mut self, policy: SolverPolicy) -> Result<Self, ReactionExtentError> {
        self.settings.solver_policy = Some(policy);
        self.settings.validate()?;
        Ok(self)
    }

    /// Replaces the policy used to seed positive log-mole coordinates.
    pub fn with_trace_seed_policy(mut self, policy: TraceSpeciesSeedPolicy) -> Self {
        self.settings.trace_seed_policy = policy;
        self
    }

    /// Returns the single trace-seed policy carried by this options bundle.
    pub fn trace_seed_policy(&self) -> TraceSpeciesSeedPolicy {
        self.settings.trace_seed_policy
    }

    /// Requests independent equilibrium-constant validation when the resolved
    /// system belongs to the validator's supported domain.
    ///
    /// This leaves the numerical solve unchanged. It only controls whether the
    /// accepted candidate carries a secondary `K_eq` comparison report, or is
    /// rejected when validation was explicitly required but unavailable.
    pub fn with_keq_validation_mode(mut self, mode: EquilibriumConstantValidationMode) -> Self {
        self.settings.keq_validation_mode = mode;
        self
    }

    /// Enables or disables stage timing in the immutable result report.
    ///
    /// Timing is off by default so ordinary production solves do not pay for
    /// repeated clock reads. Enable it for characterization and bottleneck
    /// investigations.
    pub fn with_timing_mode(mut self, mode: EquilibriumTimingMode) -> Self {
        self.timing_mode = mode;
        self
    }

    /// Enables bounded typed diagnostics for this solve transaction.
    ///
    /// Diagnostics remain observational: they never change the numerical
    /// formulation, phase-control policy, or acceptance contract.
    pub fn with_diagnostics(mut self, diagnostics: EquilibriumDiagnosticsOptions) -> Self {
        self.diagnostics = diagnostics;
        self
    }

    /// Returns the explicit diagnostics policy for this request.
    pub fn diagnostics_options(&self) -> &EquilibriumDiagnosticsOptions {
        &self.diagnostics
    }

    /// Returns the timing policy carried by this solve request.
    pub fn timing_mode(&self) -> EquilibriumTimingMode {
        self.timing_mode
    }

    /// Captures the effective immutable numerical contract for reports or
    /// reproducibility records. Calling this never changes the options or the
    /// attached execution handle.
    pub fn reproducibility_snapshot(&self) -> EquilibriumSolveOptionsSnapshot {
        let policy = self
            .settings
            .solver_policy
            .clone()
            .unwrap_or_else(|| SolverPolicy::production_default(self.settings.solver));
        EquilibriumSolveOptionsSnapshot {
            preferred_legacy_backend: format!("{:?}", self.settings.solver),
            effective_backend_order: policy
                .ordered_backends()
                .into_iter()
                .map(|backend| format!("{backend:?}"))
                .collect(),
            max_iterations: self.settings.solver_params.max_iter,
            tolerance: self.settings.solver_params.tol,
            scaling_enabled: self.settings.scaling_flag,
            trace_seed_policy: format!("{:?}", self.settings.trace_seed_policy),
            continuation_seed_policy: format!("{:?}", self.settings.continuation_seed_policy),
            equilibrium_constant_validation_mode: format!(
                "{:?}",
                self.settings.keq_validation_mode
            ),
            timing_mode: format!("{:?}", self.timing_mode),
            diagnostics_mode: format!("{:?}", self.diagnostics.mode()),
            diagnostics_max_events: self.diagnostics.max_events(),
            diagnostics_range_policy: format!("{:?}", self.diagnostics.range_policy()),
            solver_budget: self.settings.solver_budget.map(|budget| {
                EquilibriumSolverBudgetSnapshot {
                    max_attempts: budget.max_attempts,
                    max_iterations_per_attempt: budget.max_iterations_per_attempt,
                    max_total_iterations: budget.max_total_iterations,
                }
            }),
            execution_control_attached: self.execution_control.is_some(),
        }
    }

    /// Returns whether this policy can execute an RST symbolic backend.
    ///
    /// Resolved thermochemical problems have symbolic expressions by
    /// construction, so the implicit production policy is RST-first. The
    /// range facade uses this bit to decide whether symbolic preparation can
    /// be retained and parameter-updated across points.
    pub(crate) fn prepares_rst_backend(&self) -> bool {
        self.settings
            .solver_policy
            .as_ref()
            .map(|policy| {
                policy
                    .ordered_backends()
                    .iter()
                    .any(|backend| matches!(backend, SolverBackend::RustedSciThe(_)))
            })
            .unwrap_or(true)
    }
}

/// Narrow validated wrapper for phase-control policy exposed at the facade boundary.
///
/// The inner `PhaseManager` remains the canonical backend controller, but the
/// public workflow should traffic in this policy type so callers do not need to
/// know about mutable solver internals or historical phase-control helpers.
#[derive(Clone)]
pub struct PhaseControlPolicy {
    phase_manager: PhaseManager,
}

impl fmt::Debug for PhaseControlPolicy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PhaseControlPolicy").finish_non_exhaustive()
    }
}

impl Default for PhaseControlPolicy {
    fn default() -> Self {
        Self {
            phase_manager: PhaseManager::default(),
        }
    }
}

impl PhaseControlPolicy {
    /// Creates one policy wrapper from the backend controller.
    pub(crate) fn new(phase_manager: PhaseManager) -> Result<Self, ReactionExtentError> {
        phase_manager.validate()?;
        Ok(Self { phase_manager })
    }

    /// Convenience constructor for explicit phase hysteresis thresholds.
    pub fn with_explicit_hysteresis(
        phase_eps: f64,
        dg_create: f64,
        dg_keep: f64,
    ) -> Result<Self, ReactionExtentError> {
        let phase_manager = PhaseManager::new(phase_eps, dg_create, dg_keep);
        Self::new(phase_manager)
    }

    /// Convenience constructor for temperature-scaled hysteresis thresholds.
    pub fn with_temperature_scaled_hysteresis(
        phase_eps: f64,
        create_rt_factor: f64,
        keep_rt_factor: f64,
    ) -> Result<Self, ReactionExtentError> {
        let phase_manager = PhaseManager::with_temperature_scaled_hysteresis(
            phase_eps,
            create_rt_factor,
            keep_rt_factor,
        );
        Self::new(phase_manager)
    }

    /// Consumes the wrapper for the crate-internal immutable phase runner.
    pub(crate) fn into_phase_manager(self) -> PhaseManager {
        self.phase_manager
    }

    /// Changes the phase destruction threshold through the typed policy.
    pub fn with_phase_epsilon(mut self, phase_eps: f64) -> Result<Self, ReactionExtentError> {
        self.phase_manager.phase_eps = phase_eps;
        self.phase_manager.validate()?;
        Ok(self)
    }

    /// Changes the bounded outer-loop budget through the typed policy.
    pub fn with_max_phase_iterations(
        mut self,
        iterations: usize,
    ) -> Result<Self, ReactionExtentError> {
        self.phase_manager.max_phase_iterations = iterations;
        self.phase_manager.validate()?;
        Ok(self)
    }

    /// Chooses the initial active-phase policy through the typed boundary.
    pub fn with_initial_phase_set(
        mut self,
        initial_phase_set: InitialPhaseSet,
    ) -> Result<Self, ReactionExtentError> {
        self.phase_manager.initial_phase_set = initial_phase_set;
        self.phase_manager.validate()?;
        Ok(self)
    }
}

/// High-level production pipeline error.
#[derive(Debug)]
pub enum PhaseEquilibriumPipelineError {
    /// The phase/specification layer rejected the requested system.
    Resolve(SubstanceSystemFactoryError),
    /// The bridge or solver rejected the resolved problem.
    Solve(ReactionExtentError),
}

impl fmt::Display for PhaseEquilibriumPipelineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Resolve(error) => write!(f, "phase-system resolution failed: {error}"),
            Self::Solve(error) => write!(f, "phase-equilibrium solve failed: {error}"),
        }
    }
}

impl std::error::Error for PhaseEquilibriumPipelineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Resolve(error) => Some(error),
            Self::Solve(error) => Some(error),
        }
    }
}

impl From<SubstanceSystemFactoryError> for PhaseEquilibriumPipelineError {
    fn from(value: SubstanceSystemFactoryError) -> Self {
        Self::Resolve(value)
    }
}

impl From<ReactionExtentError> for PhaseEquilibriumPipelineError {
    fn from(value: ReactionExtentError) -> Self {
        Self::Solve(value)
    }
}

#[derive(Debug, Clone)]
enum PipelineInitialComposition {
    Dense(Vec<f64>),
    Sparse(Vec<(PhaseComponentId, f64)>),
}

/// High-level one-shot request that resolves phase specs and then solves them.
///
/// The caller supplies phase declarations, lookup policy, the initial
/// composition in canonical layout order, and fixed-`P,T` solve controls. The
/// builder performs resolve -> build -> solve as one transactional pipeline.
#[derive(Clone)]
pub struct PhaseEquilibriumPipelineRequest {
    spec: SubstanceSystemSpec,
    initial_composition: PipelineInitialComposition,
    conditions: EquilibriumConditions,
    model_policy: SupportedPhaseModelPolicy,
    solve_options: EquilibriumSolveOptions,
    solve_mode: PhaseEquilibriumSolveMode,
    multi_start_seeds: Vec<LogMolesInitialGuess>,
    repository: Option<Arc<ThermoRepository>>,
}

impl PhaseEquilibriumPipelineRequest {
    /// Creates one pipeline request with default numerical and model policies.
    pub fn new(
        spec: SubstanceSystemSpec,
        initial_moles: Vec<f64>,
        conditions: EquilibriumConditions,
    ) -> Self {
        Self {
            spec,
            initial_composition: PipelineInitialComposition::Dense(initial_moles),
            conditions,
            model_policy: SupportedPhaseModelPolicy::default(),
            solve_options: EquilibriumSolveOptions::default(),
            solve_mode: PhaseEquilibriumSolveMode::fixed_declared_phases(),
            multi_start_seeds: Vec::new(),
            repository: None,
        }
    }

    /// Creates a production pipeline from an element-search report and an
    /// explicit physical phase plan.
    ///
    /// The selection report supplies exact record keys and library
    /// provenance; the phase plan supplies the activity model.  Keeping both
    /// inputs mandatory prevents the top-level API from silently turning an
    /// arbitrary element search into an ideal-gas or condensed calculation.
    pub fn from_candidate_selection(
        repository: Arc<ThermoRepository>,
        selection: &EquilibriumCandidateSelectionReport,
        phase_plan: &EquilibriumCandidatePhasePlan,
        initial_moles: Vec<f64>,
        conditions: EquilibriumConditions,
    ) -> Result<Self, SubstanceSystemFactoryError> {
        let spec = phase_plan.build_spec(selection)?;
        Ok(Self::new(spec, initial_moles, conditions).with_repository(repository))
    }

    /// Creates a pipeline request with phase-qualified initial inventory.
    ///
    /// Unlike the dense constructor, this form does not make callers predict
    /// the post-resolution vector order. Missing declared components receive
    /// zero physical moles; duplicate or unknown component identities are
    /// rejected transactionally when the resolved layout is built.
    pub fn new_with_sparse_initial_composition(
        spec: SubstanceSystemSpec,
        entries: Vec<(PhaseComponentId, f64)>,
        conditions: EquilibriumConditions,
    ) -> Self {
        Self {
            spec,
            initial_composition: PipelineInitialComposition::Sparse(entries),
            conditions,
            model_policy: SupportedPhaseModelPolicy::default(),
            solve_options: EquilibriumSolveOptions::default(),
            solve_mode: PhaseEquilibriumSolveMode::fixed_declared_phases(),
            multi_start_seeds: Vec::new(),
            repository: None,
        }
    }

    /// Uses an explicit immutable repository instead of the default search path.
    pub fn with_repository(mut self, repository: Arc<ThermoRepository>) -> Self {
        self.repository = Some(repository);
        self
    }

    /// Replaces only the numerical trace-coordinate policy.
    pub fn with_trace_seed_policy(mut self, policy: TraceSpeciesSeedPolicy) -> Self {
        self.solve_options = self.solve_options.with_trace_seed_policy(policy);
        self
    }

    /// Replaces only the supported physical-model policy.
    pub fn with_model_policy(mut self, policy: SupportedPhaseModelPolicy) -> Self {
        self.model_policy = policy;
        self
    }

    /// Replaces only the phase-control policy used by bounded outer loops.
    pub fn with_phase_control_policy(mut self, policy: PhaseControlPolicy) -> Self {
        self.solve_mode = PhaseEquilibriumSolveMode::BoundedPhaseControl(policy);
        self
    }

    /// Explicitly returns the request to the safe fixed-declared-phases path.
    pub fn with_fixed_declared_phases(mut self) -> Self {
        self.solve_mode = PhaseEquilibriumSolveMode::fixed_declared_phases();
        self
    }

    /// Replaces only numerical backend and acceptance settings.
    pub fn with_solve_options(mut self, options: EquilibriumSolveOptions) -> Self {
        self.solve_options = options;
        self
    }

    /// Borrows the typed numerical options so an execution handle can be
    /// attached by a frontend without exposing the request fields.
    pub fn solve_options(&self) -> &EquilibriumSolveOptions {
        &self.solve_options
    }

    /// Requests an explicit deterministic multi-start solve for fixed declared
    /// phases.
    ///
    /// Every seed must use the canonical component order of this request and
    /// contain log-moles, not ordinary mole numbers. The solver reuses the
    /// prepared formulation and publishes the per-seed comparison report. This
    /// recovery policy is intentionally unavailable for bounded phase control:
    /// retrying an entire active-set lifecycle needs its own transition policy.
    pub fn with_multi_start_seeds(mut self, seeds: Vec<LogMolesInitialGuess>) -> Self {
        self.multi_start_seeds = seeds;
        self
    }

    /// Legacy compatibility shim for callers that still select the whole mode enum directly.
    #[deprecated(note = "use with_fixed_declared_phases() or with_phase_control_policy() instead")]
    pub fn with_solve_mode(mut self, mode: PhaseEquilibriumSolveMode) -> Self {
        self.solve_mode = mode;
        self
    }

    /// Resolves phase specifications into immutable lookup data.
    pub fn resolve(self) -> Result<ResolvedPhaseSystem, PhaseEquilibriumPipelineError> {
        let Self {
            spec, repository, ..
        } = self;
        let resolved = match repository {
            Some(repository) => {
                SubstanceSystemFactory::resolve_phase_system_with_repository(spec, repository)?
            }
            None => SubstanceSystemFactory::resolve_phase_system(spec)?,
        };
        Ok(resolved)
    }

    /// Resolves and solves the complete pipeline in one transactional pass.
    pub fn solve(self) -> Result<ResolvedPhaseEquilibriumOutcome, PhaseEquilibriumPipelineError> {
        let started = std::time::Instant::now();
        let Self {
            spec,
            initial_composition,
            conditions,
            model_policy,
            solve_options,
            solve_mode,
            multi_start_seeds,
            repository,
        } = self;
        let execution_control = solve_options.execution_control().cloned();
        if let Some(control) = &execution_control {
            control.check_cancelled()?;
        }

        let lookup_started = std::time::Instant::now();
        let resolved = match repository {
            Some(repository) => {
                SubstanceSystemFactory::resolve_phase_system_with_repository(spec, repository)?
            }
            None => SubstanceSystemFactory::resolve_phase_system(spec)?,
        };
        let layout = crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::MultiphaseEquilibriumLayout::new(
            resolved.phase_specs().to_vec(),
        )?;
        if let Some(control) = &execution_control {
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::RepositoryLookup,
                None,
                None,
                Some(conditions.temperature()),
            ));
            control.check_cancelled()?;
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::FormulationPreparation,
                None,
                None,
                Some(conditions.temperature()),
            ));
        }
        let initial_composition = match initial_composition {
            PipelineInitialComposition::Dense(initial_moles) => {
                MultiphaseInitialComposition::from_dense(&layout, initial_moles)?
            }
            PipelineInitialComposition::Sparse(entries) => {
                MultiphaseInitialComposition::from_sparse(&layout, entries)?
            }
        };
        let lookup_elapsed = lookup_started.elapsed();
        let request =
            ResolvedPhaseEquilibriumRequest::new(&resolved, conditions, initial_composition)
                .with_model_policy(model_policy)
                .with_solve_options(solve_options)
                .with_multi_start_seeds(multi_start_seeds)
                .with_fixed_declared_phases();
        let request = match solve_mode {
            PhaseEquilibriumSolveMode::FixedDeclaredPhases => request,
            PhaseEquilibriumSolveMode::BoundedPhaseControl(policy) => {
                request.with_phase_control_policy(policy)
            }
        };
        let solution = solve_resolved_pt(request)?
            .with_timing_stage(EquilibriumTimingStage::RepositoryLookup, lookup_elapsed)
            .with_timing_total(started.elapsed());

        Ok(ResolvedPhaseEquilibriumOutcome { resolved, solution })
    }

    /// Resolves once and runs a typed temperature range.
    ///
    /// Fixed mode reuses one active formulation. Bounded mode additionally
    /// carries the accepted phase set between points and records any active-set
    /// transition that requires a new projection. In both modes the resolved
    /// layout and conserved elemental inventory remain immutable.
    pub fn solve_temperature_range(
        self,
        temperatures: TemperatureGrid,
    ) -> Result<TemperatureRangeSolution, PhaseEquilibriumPipelineError> {
        let Self {
            spec,
            initial_composition,
            conditions,
            model_policy,
            solve_options,
            solve_mode,
            multi_start_seeds,
            repository,
        } = self;
        if let Some(control) = solve_options.execution_control() {
            control.check_cancelled()?;
        }
        if !multi_start_seeds.is_empty() {
            return Err(PhaseEquilibriumPipelineError::Solve(
                ReactionExtentError::InvalidProblem {
                    field: "multi_start_temperature_range",
                    message: "explicit multi-start is supported only for one fixed-P,T solve"
                        .to_string(),
                },
            ));
        }
        let phase_control_policy = match solve_mode {
            PhaseEquilibriumSolveMode::FixedDeclaredPhases => None,
            PhaseEquilibriumSolveMode::BoundedPhaseControl(policy) => Some(policy),
        };

        let resolved = match repository {
            Some(repository) => {
                SubstanceSystemFactory::resolve_phase_system_with_repository(spec, repository)?
            }
            None => SubstanceSystemFactory::resolve_phase_system(spec)?,
        };
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
        let initial_composition = match initial_composition {
            PipelineInitialComposition::Dense(initial_moles) => {
                MultiphaseInitialComposition::from_dense(&layout, initial_moles)?
            }
            PipelineInitialComposition::Sparse(entries) => {
                MultiphaseInitialComposition::from_sparse(&layout, entries)?
            }
        };
        let request = TemperatureRangeRequest::new(
            &resolved,
            initial_composition,
            conditions.pressure(),
            conditions.reference_pressure(),
            temperatures,
        )
        .map(|request| {
            request
                .with_model_policy(model_policy)
                .with_solve_options(solve_options)
        })?;
        let request = match phase_control_policy {
            Some(policy) => request.with_phase_control_policy(policy),
            None => request,
        };

        request.solve().map_err(PhaseEquilibriumPipelineError::from)
    }
}

/// Fully typed high-level pipeline outcome.
#[derive(Debug, Clone)]
pub struct ResolvedPhaseEquilibriumOutcome {
    resolved: ResolvedPhaseSystem,
    solution: MultiphaseEquilibriumSolution,
}

impl ResolvedPhaseEquilibriumOutcome {
    /// Immutable resolved phase system and its provenance.
    pub fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    /// Immutable lookup provenance forwarded from the resolved phase system.
    pub fn lookup_report(&self) -> &ResolvedPhaseSystemReport {
        self.resolved.report()
    }

    /// Immutable accepted phase-equilibrium result.
    pub fn solution(&self) -> &MultiphaseEquilibriumSolution {
        &self.solution
    }

    /// Immutable stage timing evidence for this pipeline solve.
    pub fn timing_report(&self) -> &EquilibriumTimingReport {
        self.solution.timing_report()
    }

    /// Consumes the outcome and keeps only the accepted solution.
    pub fn into_solution(self) -> MultiphaseEquilibriumSolution {
        self.solution
    }
}

/// Complete public input for one resolved fixed-pressure, fixed-temperature
/// equilibrium solve.
#[derive(Clone)]
pub struct ResolvedPhaseEquilibriumRequest<'a> {
    resolved: &'a ResolvedPhaseSystem,
    conditions: EquilibriumConditions,
    initial_composition: MultiphaseInitialComposition,
    model_policy: SupportedPhaseModelPolicy,
    solve_options: EquilibriumSolveOptions,
    solve_mode: PhaseEquilibriumSolveMode,
    multi_start_seeds: Vec<LogMolesInitialGuess>,
}

impl<'a> ResolvedPhaseEquilibriumRequest<'a> {
    /// Creates one explicit solver request. The resolved system remains
    /// borrowed and immutable for the full build/solve transaction.
    pub fn new(
        resolved: &'a ResolvedPhaseSystem,
        conditions: EquilibriumConditions,
        initial_composition: MultiphaseInitialComposition,
    ) -> Self {
        Self {
            resolved,
            conditions,
            initial_composition,
            model_policy: SupportedPhaseModelPolicy::default(),
            solve_options: EquilibriumSolveOptions::default(),
            solve_mode: PhaseEquilibriumSolveMode::fixed_declared_phases(),
            multi_start_seeds: Vec::new(),
        }
    }

    /// Replaces only the numerical trace-coordinate policy.
    pub fn with_trace_seed_policy(mut self, policy: TraceSpeciesSeedPolicy) -> Self {
        self.solve_options = self.solve_options.with_trace_seed_policy(policy);
        self
    }

    /// Replaces only the supported physical-model policy.
    pub fn with_model_policy(mut self, policy: SupportedPhaseModelPolicy) -> Self {
        self.model_policy = policy;
        self
    }

    /// Replaces only the phase-control policy used by bounded outer loops.
    pub fn with_phase_control_policy(mut self, policy: PhaseControlPolicy) -> Self {
        self.solve_mode = PhaseEquilibriumSolveMode::BoundedPhaseControl(policy);
        self
    }

    /// Explicitly returns the request to the safe fixed-declared-phases path.
    pub fn with_fixed_declared_phases(mut self) -> Self {
        self.solve_mode = PhaseEquilibriumSolveMode::fixed_declared_phases();
        self
    }

    /// Replaces only numerical backend and acceptance settings.
    pub fn with_solve_options(mut self, options: EquilibriumSolveOptions) -> Self {
        self.solve_options = options;
        self
    }

    /// Requests an explicit deterministic multi-start solve for fixed declared
    /// phases. Seeds use the canonical component order and log-mole units.
    /// Bounded phase control rejects this option because retrying a complete
    /// active-set lifecycle needs a separate transition policy.
    pub fn with_multi_start_seeds(mut self, seeds: Vec<LogMolesInitialGuess>) -> Self {
        self.multi_start_seeds = seeds;
        self
    }

    /// Legacy compatibility shim for callers that still select the whole mode enum directly.
    #[deprecated(note = "use with_fixed_declared_phases() or with_phase_control_policy() instead")]
    pub fn with_solve_mode(mut self, mode: PhaseEquilibriumSolveMode) -> Self {
        self.solve_mode = mode;
        self
    }

    /// Borrow the immutable resolved system.
    pub fn resolved(&self) -> &'a ResolvedPhaseSystem {
        self.resolved
    }

    /// Fixed thermodynamic conditions.
    pub fn conditions(&self) -> EquilibriumConditions {
        self.conditions
    }

    /// Physical initial inventory in canonical layout order.
    pub fn initial_composition(&self) -> &MultiphaseInitialComposition {
        &self.initial_composition
    }

    /// Numerical policy to be applied after the physical inventory is validated.
    pub fn solve_options(&self) -> &EquilibriumSolveOptions {
        &self.solve_options
    }
}

/// Builds and solves one resolved phase system through the canonical bridge.
///
/// All preparation stays transactional: a failed lookup, validation, or
/// backend attempt returns an error without modifying the supplied resolved
/// phase system or publishing a partial result.
pub fn solve_resolved_pt(
    request: ResolvedPhaseEquilibriumRequest<'_>,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    let timing_mode = request.solve_options.timing_mode();
    let diagnostics = request.solve_options.diagnostics_options().clone();
    let started = std::time::Instant::now();
    let execution_control = request.solve_options.execution_control().cloned();
    if let Some(control) = &execution_control {
        control.check_cancelled()?;
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::PointStarted,
            Some(0),
            Some(1),
            Some(request.conditions.temperature()),
        ));
    }
    let solved = match request.solve_mode {
        PhaseEquilibriumSolveMode::FixedDeclaredPhases => {
            let trace_seed_policy = request.solve_options.trace_seed_policy();
            let bundle = crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
                build_phase_equilibrium_problem_with_timing(PhaseEquilibriumBuildRequest::new(
                request.resolved,
                request.conditions,
                request.initial_composition,
                trace_seed_policy,
                request.model_policy,
            )?, timing_mode)?;
            let settings = request.solve_options.into_settings();
            let solved = if request.multi_start_seeds.is_empty() {
                bundle.solve_with(|configured| *configured = settings)
            } else {
                bundle.solve_with_initial_guesses(request.multi_start_seeds, |configured| {
                    *configured = settings
                })
            };
            solved
                .and_then(|bundle| bundle.into_multiphase_solution())
                .map(|solution| solution.with_timing_total(started.elapsed()))
        }
        PhaseEquilibriumSolveMode::BoundedPhaseControl(phase_control_policy) => {
            if !request.multi_start_seeds.is_empty() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "multi_start_phase_control",
                    message: "explicit multi-start is supported only for fixed declared phases"
                        .to_string(),
                });
            }
            let trace_seed_policy = request.solve_options.trace_seed_policy();
            let bundle = crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
                build_phase_equilibrium_problem_with_timing(PhaseEquilibriumBuildRequest::new(
                request.resolved,
                request.conditions,
                request.initial_composition,
                trace_seed_policy,
                request.model_policy,
            )?, timing_mode)?;
            let settings = request.solve_options.into_settings();
            bundle
                .solve_with_bounded_phase_control_with_diagnostics(
                    |configured| *configured = settings,
                    |configured| *configured = phase_control_policy.into_phase_manager(),
                    diagnostics,
                )
                .map(|solution| solution.with_timing_total(started.elapsed()))
        }
    };
    if let Some(control) = &execution_control {
        control.check_cancelled()?;
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::PointAccepted,
            Some(0),
            Some(1),
            Some(request.conditions.temperature()),
        ));
    }
    solved
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
        EquilibriumDiagnosticEvent, EquilibriumDiagnosticsMode, EquilibriumDiagnosticsOptions,
    };
    use crate::Thermodynamics::ChemEquilibrium::prelude::{
        LegacyEquilibriumSolver, LogMolesInitialGuess, RustedSciTheSolver, SolverBackend,
        SolverPolicy,
    };
    use crate::Thermodynamics::User_PhaseOrSolution::{
        PhaseModel, PhaseSpec, SubstanceSystemFactory, SubstanceSystemSpec,
        SubstanceSystemSpecBuilder, SubstancesContainer,
    };
    use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
    use crate::Thermodynamics::physical_state::PhysicalState;
    use crate::Thermodynamics::thermo_lib_api::ThermoData;

    #[test]
    fn solve_options_round_trip_explicit_cascade_budget() {
        let budget = SolverCascadeBudget::new(2, 10, 1);
        let options = EquilibriumSolveOptions::new()
            .with_solver_budget(budget)
            .expect("positive cascade budget must validate");
        assert_eq!(options.into_settings().solver_budget, Some(budget));

        let error = EquilibriumSolveOptions::new()
            .with_solver_budget(SolverCascadeBudget::new(0, 10, 10))
            .expect_err("zero attempt budget must be rejected");
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "solver_budget",
                ..
            }
        ));
    }

    #[test]
    fn cancelled_pipeline_stops_before_repository_resolution_for_point_and_range() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let conditions = EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap();
        let control =
            crate::Thermodynamics::ChemEquilibrium::prelude::EquilibriumExecutionControl::new();
        control.request_cancel();
        let options = EquilibriumSolveOptions::new().with_execution_control(control);

        let point_error = PhaseEquilibriumPipelineRequest::new(spec.clone(), vec![1.0], conditions)
            .with_solve_options(options.clone())
            .solve()
            .expect_err("cancelled point must not enter repository lookup");
        assert!(matches!(
            point_error,
            PhaseEquilibriumPipelineError::Solve(ReactionExtentError::Cancelled)
        ));

        let range_error = PhaseEquilibriumPipelineRequest::new(spec, vec![1.0], conditions)
            .with_solve_options(options)
            .solve_temperature_range(TemperatureGrid::new(vec![500.0, 600.0]).unwrap())
            .expect_err("cancelled range must not enter repository lookup");
        assert!(matches!(
            range_error,
            PhaseEquilibriumPipelineError::Solve(ReactionExtentError::Cancelled)
        ));
    }

    #[test]
    fn production_prelude_exposes_complete_backend_policy_contract() {
        let rst = SolverPolicy::Single(SolverBackend::RustedSciThe(
            RustedSciTheSolver::MinpackLevenbergMarquardt,
        ));
        let legacy = SolverPolicy::Single(SolverBackend::Legacy(LegacyEquilibriumSolver::NR));

        EquilibriumSolveOptions::new()
            .with_solver_policy(rst)
            .expect("RST policy exported by the production prelude must validate");
        EquilibriumSolveOptions::new()
            .with_solver_policy(legacy)
            .expect("legacy fallback policy exported by the production prelude must validate");
    }

    #[test]
    fn typed_ideal_solution_spec_resolves_through_the_repository_boundary() {
        // This intentionally bypasses the compatibility `MultiPhase` builder:
        // its historical phase-nature enum cannot express a mixing model.
        // Production callers and the GUI use `from_phases`, which must retain
        // the explicit semantic model through ordinary local lookup.
        let phase = PhaseSpec::ideal_solution(
            PhaseId::new(Some("oxide_solution".to_string())),
            vec!["AL(cr)".to_string(), "AL2O3(a)".to_string()],
            PhysicalState::Solid,
        )
        .expect("the typed ideal-solution declaration must be valid");
        let spec = SubstanceSystemSpec::from_phases(vec![phase])
            .expect("typed phase specification must validate")
            .with_lookup_policy(vec!["NASA_cond".to_string()], Vec::new(), None, false);
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            spec,
            ThermoData::try_default_repository().expect("bundled repository must load"),
        )
        .expect("bundled condensed records must resolve through the typed specification");

        assert_eq!(resolved.phase_specs().len(), 1);
        assert_eq!(resolved.phase_specs()[0].model(), PhaseModel::IdealSolution);
        assert_eq!(
            resolved.phase_specs()[0].components(),
            ["AL(cr)", "AL2O3(a)"]
        );
        let data = resolved
            .phase_data()
            .get(&Some("oxide_solution".to_string()))
            .expect("resolved phase payload must retain its declared identity");
        assert_eq!(data.substances, vec!["AL(cr)", "AL2O3(a)"]);
        assert_eq!(
            resolved.report().nist_fallback_enabled(),
            false,
            "the contract fixture must remain fully local and deterministic"
        );
    }

    #[test]
    fn enabled_timing_is_published_on_pipeline_solution() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let options =
            EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled);
        let outcome = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_solve_options(options)
        .solve()
        .unwrap();

        let timing = outcome.timing_report();
        assert!(timing.enabled());
        assert!(timing.total() > std::time::Duration::ZERO);
        assert!(timing.thermochemistry_preparation() > std::time::Duration::ZERO);
        assert!(timing.numeric_closure_construction() > std::time::Duration::ZERO);
        assert!(timing.symbolic_construction() > std::time::Duration::ZERO);
        assert!(timing.equation_construction() > std::time::Duration::ZERO);
        assert!(timing.nonlinear_solve() > std::time::Duration::ZERO);
        assert!(timing.validation() >= std::time::Duration::ZERO);
        assert!(timing.postprocessing() > std::time::Duration::ZERO);
    }

    #[test]
    fn bounded_phase_control_publishes_opt_in_lifecycle_diagnostics() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let options = EquilibriumSolveOptions::new()
            .with_timing_mode(EquilibriumTimingMode::Enabled)
            .with_diagnostics(EquilibriumDiagnosticsOptions::enabled(
                EquilibriumDiagnosticsMode::PhaseLifecycle,
            ));
        let solution = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_solve_options(options)
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve()
        .unwrap()
        .into_solution();

        let diagnostics = solution
            .diagnostics_report()
            .expect("enabled diagnostics must be attached to the accepted solution");
        assert!(matches!(
            diagnostics.events().first(),
            Some(EquilibriumDiagnosticEvent::SolveStarted { .. })
        ));
        assert!(
            diagnostics.events().iter().any(|event| matches!(
                event,
                EquilibriumDiagnosticEvent::StabilityEvaluated { .. }
            ))
        );
        assert!(matches!(
            diagnostics.events().last(),
            Some(EquilibriumDiagnosticEvent::SolveAccepted { .. })
        ));
        let rendered = crate::Thermodynamics::ChemEquilibrium::
            equilibrium_diagnostics_display::format_solution_diagnostics(&solution)
            .expect("enabled diagnostics must render through the presentation adapter");
        assert!(rendered.contains("solve started"));
        assert!(rendered.contains("stability iteration"));
        assert!(rendered.contains("solve accepted"));
        assert!(rendered.contains("timing ms: total="));
    }

    #[test]
    fn diagnostics_remain_absent_when_the_default_options_are_used() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let solution = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default())
        .solve()
        .unwrap()
        .into_solution();

        assert!(solution.diagnostics_report().is_none());
    }

    #[test]
    fn resolve_and_solve_pipeline_round_trips_one_gas_phase_with_two_species() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(true)
        .build()
        .unwrap();
        let request = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        );

        let outcome = request.solve().unwrap();

        assert_eq!(outcome.resolved().phase_specs().len(), 1);
        assert_eq!(outcome.solution().component_moles().len(), 2);
        assert_eq!(outcome.solution().build_report().components().len(), 2);
    }

    #[test]
    fn pipeline_and_resolved_pt_facades_publish_the_same_fixed_solution() {
        // The pipeline is intentionally exercised against the already
        // resolved facade with one explicit backend. This isolates the
        // resolve/build orchestration difference from backend cascade policy
        // and makes parity a durable API contract rather than a coincidence
        // of two independent defaults.
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let conditions = EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap();
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
            .expect("single legacy backend policy must validate");
        let pipeline = PhaseEquilibriumPipelineRequest::new(spec, vec![0.79, 0.21], conditions)
            .with_solve_options(options.clone());
        let resolved = pipeline
            .clone()
            .resolve()
            .expect("pipeline resolution must succeed");
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
        let composition = MultiphaseInitialComposition::from_dense(&layout, vec![0.79, 0.21])
            .expect("direct facade composition must match the resolved layout");

        let pipeline_solution = pipeline
            .solve()
            .expect("pipeline facade must solve")
            .into_solution();
        let direct_solution = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(&resolved, conditions, composition)
                .with_solve_options(options),
        )
        .expect("resolved facade must solve");

        assert_eq!(
            pipeline_solution.component_moles().len(),
            direct_solution.component_moles().len()
        );
        for (pipeline_moles, direct_moles) in pipeline_solution
            .component_moles()
            .iter()
            .zip(direct_solution.component_moles())
        {
            assert!((pipeline_moles - direct_moles).abs() <= 1e-8);
        }
        assert!(
            (pipeline_solution
                .accepted_solution()
                .validation()
                .residual_l2_norm
                - direct_solution
                    .accepted_solution()
                    .validation()
                    .residual_l2_norm)
                .abs()
                <= 1e-12
        );
        assert_eq!(
            pipeline_solution.build_report().layout_fingerprint(),
            direct_solution.build_report().layout_fingerprint()
        );
    }

    #[test]
    fn fixed_pipeline_publishes_explicit_multi_start_evidence() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let conditions = EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap();
        let seed = LogMolesInitialGuess::from_initial_moles(&[0.79, 0.21]).unwrap();

        let outcome = PhaseEquilibriumPipelineRequest::new(spec, vec![0.79, 0.21], conditions)
            .with_multi_start_seeds(vec![seed.clone(), seed])
            .solve()
            .unwrap();

        let report = outcome
            .solution()
            .multi_start_report()
            .expect("fixed multi-start must publish its comparison report");
        assert_eq!(report.attempts.len(), 2);
        assert!(report.attempts.iter().all(|attempt| attempt.accepted));
        assert!(report.selected_start < report.attempts.len());
    }

    #[test]
    fn multi_start_is_rejected_for_bounded_phase_control() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let seed = LogMolesInitialGuess::from_initial_moles(&[0.79, 0.21]).unwrap();

        let error = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_multi_start_seeds(vec![seed])
        .solve()
        .unwrap_err();

        assert!(matches!(
            error,
            PhaseEquilibriumPipelineError::Solve(ReactionExtentError::InvalidProblem {
                field: "multi_start_phase_control",
                ..
            })
        ));
    }

    #[test]
    fn multi_start_rejects_a_seed_with_the_wrong_component_dimension() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let bad_seed = LogMolesInitialGuess::from_initial_moles(&[1.0]).unwrap();

        let error = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_multi_start_seeds(vec![bad_seed])
        .solve()
        .unwrap_err();

        assert!(matches!(
            error,
            PhaseEquilibriumPipelineError::Solve(ReactionExtentError::InvalidProblem {
                field: "multi_start_seed_dimension",
                ..
            })
        ));
    }

    #[test]
    fn sparse_pipeline_inventory_matches_dense_ordered_inventory() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(true)
        .build()
        .unwrap();
        let conditions = EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap();
        let dense =
            PhaseEquilibriumPipelineRequest::new(spec.clone(), vec![0.79, 0.21], conditions)
                .solve()
                .unwrap();
        let sparse = PhaseEquilibriumPipelineRequest::new_with_sparse_initial_composition(
            spec,
            vec![
                (PhaseComponentId::new(PhaseId::new(None), "N2"), 0.79),
                (PhaseComponentId::new(PhaseId::new(None), "O2"), 0.21),
            ],
            conditions,
        )
        .solve()
        .unwrap();

        assert_eq!(
            dense.solution().component_moles(),
            sparse.solution().component_moles()
        );
        assert_eq!(
            dense.solution().layout_fingerprint(),
            sparse.solution().layout_fingerprint()
        );
    }

    #[test]
    fn sparse_pipeline_inventory_rejects_duplicate_phase_components() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(true)
        .build()
        .unwrap();
        let duplicate = PhaseComponentId::new(PhaseId::new(None), "N2");
        let error = PhaseEquilibriumPipelineRequest::new_with_sparse_initial_composition(
            spec,
            vec![(duplicate.clone(), 0.79), (duplicate, 0.21)],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .solve()
        .unwrap_err();

        assert!(matches!(
            error,
            PhaseEquilibriumPipelineError::Solve(ReactionExtentError::InvalidProblem {
                field: "initial_composition",
                ..
            })
        ));
    }

    #[test]
    fn phase_control_policy_rejects_an_invalid_hysteresis_order() {
        assert!(PhaseControlPolicy::with_explicit_hysteresis(1e-6, 1.0, 0.5).is_err());
    }

    #[test]
    fn phase_control_policy_rejects_invalid_scalar_limits_at_construction() {
        let mut manager = PhaseManager::default();
        manager.phase_eps = 0.0;
        assert!(PhaseControlPolicy::new(manager).is_err());

        let mut manager = PhaseManager::default();
        manager.max_phase_iterations = 0;
        assert!(PhaseControlPolicy::new(manager).is_err());
    }

    #[test]
    fn typed_policy_builders_reject_invalid_limits_before_solving() {
        assert!(
            EquilibriumSolveOptions::default()
                .with_max_iterations(0)
                .is_err()
        );
        assert!(
            EquilibriumSolveOptions::default()
                .with_tolerance(0.0)
                .is_err()
        );
        assert!(
            PhaseControlPolicy::default()
                .with_phase_epsilon(0.0)
                .is_err()
        );
        assert!(
            PhaseControlPolicy::default()
                .with_max_phase_iterations(0)
                .is_err()
        );
    }

    #[test]
    fn trace_seed_policy_has_one_typed_source_of_truth() {
        let policy = TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
            fraction: 1e-8,
            minimum_floor: 1e-16,
        };
        let options = EquilibriumSolveOptions::default().with_trace_seed_policy(policy);
        assert_eq!(options.trace_seed_policy(), policy);
    }

    #[test]
    fn phase_control_policy_rejects_indices_after_resolved_phase_count_is_known() {
        let manager = PhaseManager {
            initial_phase_set: crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::
                InitialPhaseSet::Explicit {
                active: vec![crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex::new(2, 3).unwrap()],
                excluded: Vec::new(),
            },
            ..PhaseManager::default()
        };

        PhaseControlPolicy::new(manager)
            .unwrap()
            .into_phase_manager()
            .validate_for_phase_count(2)
            .expect_err("resolved phase bounds must be checked before active-set construction");
    }

    #[test]
    fn selecting_a_preferred_backend_clears_an_explicit_policy() {
        let settings = EquilibriumSolverSettings {
            solver_policy: Some(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR))),
            ..EquilibriumSolverSettings::default()
        };
        let options = EquilibriumSolveOptions::from_settings(settings)
            .unwrap()
            .with_solver_backend(Solvers::TR);
        let settings = options.into_settings();

        assert_eq!(settings.solver, Solvers::TR);
        assert!(settings.solver_policy.is_none());
    }

    #[test]
    fn pipeline_outcome_exposes_the_original_lookup_provenance() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(true)
        .build()
        .unwrap();
        let request = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        );

        let outcome = request.solve().unwrap();

        assert_eq!(outcome.lookup_report(), outcome.resolved().report());
        assert_eq!(
            outcome.solution().build_report().lookup_report(),
            outcome.lookup_report()
        );
    }

    #[test]
    fn solve_mode_defaults_to_the_explicit_fixed_declared_phase_path() {
        assert!(matches!(
            PhaseEquilibriumSolveMode::default(),
            PhaseEquilibriumSolveMode::FixedDeclaredPhases
        ));
        assert!(matches!(
            PhaseEquilibriumSolveMode::fixed_declared_phases(),
            PhaseEquilibriumSolveMode::FixedDeclaredPhases
        ));
    }

    #[test]
    fn request_can_explicitly_reset_to_the_fixed_declared_phase_path() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(true)
        .build()
        .unwrap();
        let request = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_fixed_declared_phases();

        let outcome = request.solve().unwrap();

        assert_eq!(outcome.resolved().phase_specs().len(), 1);
        assert_eq!(outcome.solution().component_moles().len(), 2);
    }
}
