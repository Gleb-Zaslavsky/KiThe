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
    ElementInventoryCandidateSelection, EquilibriumCandidatePhasePlan,
    EquilibriumCandidateSelectionReport,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
    EquilibriumDiagnosticEvent, EquilibriumDiagnosticsCollector, EquilibriumDiagnosticsMode,
    EquilibriumDiagnosticsOptions,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_element_inventory::ElementInventory;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::{
    EquilibriumExecutionControl, EquilibriumProgressEvent, EquilibriumProgressStage,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_extensive_normalization::ExtensiveNormalization;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumSolverSettings, Solvers,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
#[cfg(test)]
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    BackendFailureKind, SolveError,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionExtentError, ReactionExtentErrorKind,
};
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
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    InitialPhaseSet, PhaseManager, PhaseSet,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::PhaseEquilibriumBuildRequest;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::SupportedPhaseModelPolicy;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::ExtensiveNormalizationRecoveryEvidence;
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

/// Policy for recovering a numerically ill-conditioned extensive inventory.
///
/// Recovery never changes the published physical problem. It first solves one
/// exactly normalized equivalent and then prefers a physical-coordinate retry.
/// If that retry remains ill-conditioned, the accepted normalized snapshot is
/// passed through the audited physical publication boundary; provenance says
/// explicitly which route produced the public result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExtensiveNormalizationPolicy {
    /// Do not attempt extensive-coordinate recovery.
    Disabled,
    /// Retry only after a classified numerical failure of the physical solve.
    OnNumericalFailure,
}

impl Default for ExtensiveNormalizationPolicy {
    fn default() -> Self {
        Self::OnNumericalFailure
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
    extensive_normalization: ExtensiveNormalizationPolicy,
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
    pub extensive_normalization_policy: String,
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
            extensive_normalization: ExtensiveNormalizationPolicy::default(),
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
            extensive_normalization: ExtensiveNormalizationPolicy::default(),
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

    /// Selects whether a classified extensive-conditioning failure may use
    /// normalized basin discovery and the audited physical publication path.
    pub fn with_extensive_normalization_policy(
        mut self,
        policy: ExtensiveNormalizationPolicy,
    ) -> Self {
        self.extensive_normalization = policy;
        self
    }

    /// Returns the extensive-conditioning recovery policy.
    pub fn extensive_normalization_policy(&self) -> ExtensiveNormalizationPolicy {
        self.extensive_normalization
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
            extensive_normalization_policy: format!("{:?}", self.extensive_normalization),
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

    /// Converts only the absolute phase-mole threshold for normalized space.
    ///
    /// The returned policy keeps its active-set budget, initial phase set, and
    /// intensive TPD hysteresis. Pair it with a request whose composition and
    /// other extensive inputs were normalized by the same
    /// [`ExtensiveNormalization`]. This is an explicit representation mapping,
    /// not an automatic retry or default solver route.
    pub fn normalized_for_extensive_representation(
        mut self,
        normalization: ExtensiveNormalization,
    ) -> Result<Self, ReactionExtentError> {
        self.phase_manager.phase_eps =
            normalization.normalize_physical_mole_threshold(self.phase_manager.phase_eps)?;
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
    ElementInventory(ElementInventory),
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
    candidate_selection: Option<EquilibriumCandidateSelectionReport>,
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
            candidate_selection: None,
            conditions,
            model_policy: SupportedPhaseModelPolicy::default(),
            solve_options: EquilibriumSolveOptions::default(),
            solve_mode: PhaseEquilibriumSolveMode::fixed_declared_phases(),
            multi_start_seeds: Vec::new(),
            repository: None,
        }
    }

    /// Creates a P,T pipeline from a closed elemental inventory and an
    /// explicit real phase/species universe.
    ///
    /// The inventory is not expanded into an artificial species. After lookup
    /// it enters the canonical bridge as `b`, while a seed is built from the
    /// resolved real components only.
    pub fn from_element_inventory(
        spec: SubstanceSystemSpec,
        element_inventory: ElementInventory,
        conditions: EquilibriumConditions,
    ) -> Self {
        Self {
            spec,
            initial_composition: PipelineInitialComposition::ElementInventory(element_inventory),
            candidate_selection: None,
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
        Ok(Self::new(spec, initial_moles, conditions)
            .with_repository(repository)
            .with_candidate_selection(selection.clone()))
    }

    /// Creates an element-defined pipeline from a deterministic catalog
    /// selection and explicit physical phase plan.
    pub fn from_candidate_selection_with_element_inventory(
        repository: Arc<ThermoRepository>,
        selection: &EquilibriumCandidateSelectionReport,
        phase_plan: &EquilibriumCandidatePhasePlan,
        element_inventory: ElementInventory,
        conditions: EquilibriumConditions,
    ) -> Result<Self, SubstanceSystemFactoryError> {
        let spec = phase_plan.build_spec(selection)?;
        Ok(
            Self::from_element_inventory(spec, element_inventory, conditions)
                .with_repository(repository)
                .with_candidate_selection(selection.clone()),
        )
    }

    /// Creates a reproducible element-defined pipeline directly from typed
    /// selection evidence. The inventory remains the physical conservation
    /// vector; the nested report explains the selected real record universe.
    pub fn from_element_inventory_candidate_selection(
        repository: Arc<ThermoRepository>,
        selection: &ElementInventoryCandidateSelection,
        phase_plan: &EquilibriumCandidatePhasePlan,
        conditions: EquilibriumConditions,
    ) -> Result<Self, SubstanceSystemFactoryError> {
        Self::from_candidate_selection_with_element_inventory(
            repository,
            selection.report(),
            phase_plan,
            selection.inventory().clone(),
            conditions,
        )
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
            candidate_selection: None,
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

    /// Retains the immutable catalog-selection transaction in the eventual
    /// bridge build report.
    pub fn with_candidate_selection(
        mut self,
        selection: EquilibriumCandidateSelectionReport,
    ) -> Self {
        self.candidate_selection = Some(selection);
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
            candidate_selection,
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
            PipelineInitialComposition::ElementInventory(element_inventory) => {
                if !multi_start_seeds.is_empty() {
                    return Err(PhaseEquilibriumPipelineError::Solve(
                        ReactionExtentError::InvalidProblem {
                            field: "multi_start_element_inventory",
                            message: "explicit multi-start is not yet available for element-defined P,T pipelines".to_string(),
                        },
                    ));
                }
                let phase_control = match solve_mode {
                    PhaseEquilibriumSolveMode::FixedDeclaredPhases => None,
                    PhaseEquilibriumSolveMode::BoundedPhaseControl(policy) => Some(policy),
                };
                let solution = solve_resolved_pt_from_element_inventory_with_selection(
                    &resolved,
                    conditions,
                    element_inventory,
                    None,
                    candidate_selection,
                    solve_options,
                    phase_control,
                )?
                .with_timing_stage(
                    EquilibriumTimingStage::RepositoryLookup,
                    lookup_started.elapsed(),
                )
                .with_timing_total(started.elapsed());
                return Ok(ResolvedPhaseEquilibriumOutcome { resolved, solution });
            }
        };
        let lookup_elapsed = lookup_started.elapsed();
        let request =
            ResolvedPhaseEquilibriumRequest::new(&resolved, conditions, initial_composition)
                .with_model_policy(model_policy)
                .with_solve_options(solve_options)
                .with_multi_start_seeds(multi_start_seeds)
                .with_fixed_declared_phases()
                .with_candidate_selection(candidate_selection);
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
            candidate_selection: _,
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
            PipelineInitialComposition::ElementInventory(_) => {
                return Err(PhaseEquilibriumPipelineError::Solve(
                    ReactionExtentError::InvalidProblem {
                        field: "element_inventory_temperature_range",
                        message: "element-defined temperature ranges require the prepared range migration and are not available yet".to_string(),
                    },
                ));
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
    /// Rebuilds the transport wrapper after an integration layer has consumed
    /// a calculator-facade point. No numerical payload is cloned.
    pub fn from_parts(
        resolved: ResolvedPhaseSystem,
        solution: MultiphaseEquilibriumSolution,
    ) -> Self {
        Self { resolved, solution }
    }

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
    physical_input: ResolvedPhaseEquilibriumPhysicalInput,
    candidate_selection: Option<EquilibriumCandidateSelectionReport>,
    model_policy: SupportedPhaseModelPolicy,
    solve_options: EquilibriumSolveOptions,
    solve_mode: PhaseEquilibriumSolveMode,
    multi_start_seeds: Vec<LogMolesInitialGuess>,
    /// Internal numerical seed used by extensive-conditioning recovery. It is
    /// never interpreted as accepted phase-set continuation.
    recovery_basin_seed: Option<LogMolesInitialGuess>,
    /// Accepted active set carried only by a prepared continuation owner.
    /// Keeping this crate-internal prevents arbitrary callers from claiming
    /// unvalidated phase history.
    continuation_phase_set: Option<PhaseSet>,
}

#[derive(Clone)]
enum ResolvedPhaseEquilibriumPhysicalInput {
    Composition,
    ElementInventory(ElementInventory),
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
            physical_input: ResolvedPhaseEquilibriumPhysicalInput::Composition,
            candidate_selection: None,
            model_policy: SupportedPhaseModelPolicy::default(),
            solve_options: EquilibriumSolveOptions::default(),
            solve_mode: PhaseEquilibriumSolveMode::fixed_declared_phases(),
            multi_start_seeds: Vec::new(),
            recovery_basin_seed: None,
            continuation_phase_set: None,
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

    /// Carries a previously accepted numerical state into an internal range
    /// recovery without exposing mutable phase-control state publicly.
    pub(crate) fn with_continuation_state(
        mut self,
        seed: LogMolesInitialGuess,
        phase_set: Option<PhaseSet>,
    ) -> Self {
        self.recovery_basin_seed = Some(seed);
        self.continuation_phase_set = phase_set;
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

    /// Replaces only the physical conservation source with a validated
    /// elemental inventory. The composition remains a numerical seed and is
    /// checked against this inventory during bridge construction.
    pub fn with_element_inventory(mut self, inventory: ElementInventory) -> Self {
        self.physical_input = ResolvedPhaseEquilibriumPhysicalInput::ElementInventory(inventory);
        self
    }

    /// Carries immutable candidate-selection provenance into the build report.
    pub fn with_candidate_selection(
        mut self,
        selection: Option<EquilibriumCandidateSelectionReport>,
    ) -> Self {
        self.candidate_selection = selection;
        self
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
    let solved = solve_resolved_pt_transaction(request)?.with_timing_total(started.elapsed());
    if let Some(control) = &execution_control {
        control.check_cancelled()?;
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::PointAccepted,
            Some(0),
            Some(1),
            Some(solved.conditions().temperature()),
        ));
    }
    Ok(solved)
}

/// Solves one P,T point from a closed elemental inventory over an already
/// resolved real phase/species universe.
///
/// This is a typed sibling of [`solve_resolved_pt`] rather than a second
/// solver request hierarchy. The inventory is consumed by the bridge to
/// produce `A`, `b`, and a real-component numerical seed; all backend,
/// acceptance, and optional phase-control behavior remains identical.
pub fn solve_resolved_pt_from_element_inventory(
    resolved: &ResolvedPhaseSystem,
    conditions: EquilibriumConditions,
    element_inventory: ElementInventory,
    solve_options: EquilibriumSolveOptions,
    phase_control: Option<PhaseControlPolicy>,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    solve_resolved_pt_from_element_inventory_with_numerical_seed(
        resolved,
        conditions,
        element_inventory,
        None,
        solve_options,
        phase_control,
    )
}

/// Internal continuation-capable form of
/// [`solve_resolved_pt_from_element_inventory`].
///
/// `numerical_seed` affects only log-coordinate initialization. The bridge
/// still derives conserved totals exclusively from `element_inventory`.
pub(crate) fn solve_resolved_pt_from_element_inventory_with_numerical_seed(
    resolved: &ResolvedPhaseSystem,
    conditions: EquilibriumConditions,
    element_inventory: ElementInventory,
    numerical_seed: Option<MultiphaseInitialComposition>,
    solve_options: EquilibriumSolveOptions,
    phase_control: Option<PhaseControlPolicy>,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    solve_resolved_pt_from_element_inventory_with_selection(
        resolved,
        conditions,
        element_inventory,
        numerical_seed,
        None,
        solve_options,
        phase_control,
    )
}

fn solve_resolved_pt_from_element_inventory_with_selection(
    resolved: &ResolvedPhaseSystem,
    conditions: EquilibriumConditions,
    element_inventory: ElementInventory,
    numerical_seed: Option<MultiphaseInitialComposition>,
    candidate_selection: Option<EquilibriumCandidateSelectionReport>,
    solve_options: EquilibriumSolveOptions,
    phase_control: Option<PhaseControlPolicy>,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    let trace_seed_policy = solve_options.trace_seed_policy();
    let build_request = PhaseEquilibriumBuildRequest::from_element_inventory(
        resolved,
        conditions,
        element_inventory,
        trace_seed_policy,
        SupportedPhaseModelPolicy::default(),
    )?;
    let build_request = match candidate_selection {
        Some(selection) => build_request.with_candidate_selection(selection),
        None => build_request,
    };
    let build_request = match numerical_seed {
        Some(seed) => build_request.with_numerical_seed(seed)?,
        None => build_request,
    };
    match phase_control {
        Some(policy) => crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
            build_phase_equilibrium_problem(build_request)?
            .solve_with_bounded_phase_control_with_diagnostics(
            |settings| *settings = solve_options.clone().into_settings(),
            |manager| *manager = policy.into_phase_manager(),
            solve_options.diagnostics_options().clone(),
        ),
        None => crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
            build_phase_equilibrium_problem(build_request)?
            .solve_with(|settings| *settings = solve_options.into_settings())?
            .into_multiphase_solution(),
    }
}

/// Runs one complete point transaction without emitting range-level progress.
///
/// Prepared continuation drivers use this only after their cached fast path
/// fails. Progress publication stays with the outer owner, so a recovered
/// range point cannot masquerade as a nested one-point request.
pub(crate) fn solve_resolved_pt_transaction(
    request: ResolvedPhaseEquilibriumRequest<'_>,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    match solve_resolved_pt_once(request.clone()) {
        Ok(solution) => Ok(solution),
        Err(physical_error) => recover_resolved_pt_after_numerical_failure(request, physical_error),
    }
}

/// Applies canonical extensive recovery to an already observed point error.
///
/// This avoids repeating the same failed physical attempt when a prepared
/// range formulation has already produced complete typed failure evidence.
pub(crate) fn recover_resolved_pt_after_numerical_failure<'a>(
    request: ResolvedPhaseEquilibriumRequest<'a>,
    physical_error: ReactionExtentError,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    if request.solve_options.extensive_normalization_policy()
        == ExtensiveNormalizationPolicy::OnNumericalFailure
        && is_extensive_normalization_recovery_candidate(&physical_error)
    {
        recover_resolved_pt_by_extensive_normalization(request, physical_error)
    } else {
        Err(physical_error)
    }
}

fn build_phase_equilibrium_request<'a>(
    request: &ResolvedPhaseEquilibriumRequest<'a>,
    numerical_seed: MultiphaseInitialComposition,
) -> Result<PhaseEquilibriumBuildRequest<'a>, ReactionExtentError> {
    let trace_seed_policy = request.solve_options.trace_seed_policy();
    let build_request = match &request.physical_input {
        ResolvedPhaseEquilibriumPhysicalInput::Composition => PhaseEquilibriumBuildRequest::new(
            request.resolved,
            request.conditions,
            numerical_seed,
            trace_seed_policy,
            request.model_policy,
        ),
        ResolvedPhaseEquilibriumPhysicalInput::ElementInventory(inventory) => {
            PhaseEquilibriumBuildRequest::from_element_inventory(
                request.resolved,
                request.conditions,
                inventory.clone(),
                trace_seed_policy,
                request.model_policy,
            )?
            .with_numerical_seed(numerical_seed)
        }
    }?;
    Ok(match &request.candidate_selection {
        Some(selection) => build_request.with_candidate_selection(selection.clone()),
        None => build_request,
    })
}

/// Executes exactly one physical or explicitly normalized request.
///
/// Keeping retry policy outside this function prevents recursive recovery and
/// makes every individual attempt retain the ordinary solver contract.
fn solve_resolved_pt_once(
    request: ResolvedPhaseEquilibriumRequest<'_>,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    let timing_mode = request.solve_options.timing_mode();
    let diagnostics = request.solve_options.diagnostics_options().clone();
    let started = std::time::Instant::now();
    match request.solve_mode.clone() {
        PhaseEquilibriumSolveMode::FixedDeclaredPhases => {
            let bundle = crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
                build_phase_equilibrium_problem_with_timing(
                    build_phase_equilibrium_request(&request, request.initial_composition.clone())?,
                    timing_mode,
                )?;
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
            let bundle = crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
                build_phase_equilibrium_problem_with_timing(
                    build_phase_equilibrium_request(&request, request.initial_composition.clone())?,
                    timing_mode,
                )?;
            let settings = request.solve_options.into_settings();
            bundle
                .solve_with_bounded_phase_control_with_diagnostics_and_seed(
                    |configured| *configured = settings,
                    |configured| *configured = phase_control_policy.into_phase_manager(),
                    diagnostics,
                    request.recovery_basin_seed,
                    request.continuation_phase_set,
                )
                .map(|solution| solution.with_timing_total(started.elapsed()))
        }
    }
}

/// The recovery route is reserved for genuinely extensive conditioning.
///
/// A factor in `[0.1, 10]` is deliberately treated as near-unit: applying a
/// second formulation there would turn an ordinary numerical failure into an
/// implicit multi-start policy and would make the failure provenance harder to
/// interpret. The bound is expressed in the physical scale domain so the
/// exact boundary is testable without constructing a complete solve request.
const EXTENSIVE_RECOVERY_MIN_SCALE_FACTOR: f64 = 10.0;

// The normalized recovery state is used as a cross-representation witness.
// A tighter floor prevents a scale-dependent phase-control amount from being
// published when the ordinary production tolerance is sufficient for solve
// acceptance but not for strict extensive equivalence.
const EXTENSIVE_RECOVERY_MAX_SOLVER_TOLERANCE: f64 = 1.0e-12;

fn is_materially_extensive_scaled(scale: f64) -> bool {
    scale.is_finite() && scale > 0.0 && scale.ln().abs() > EXTENSIVE_RECOVERY_MIN_SCALE_FACTOR.ln()
}

fn is_extensive_normalization_recovery_candidate(error: &ReactionExtentError) -> bool {
    matches!(
        error.kind(),
        ReactionExtentErrorKind::AllBackendsFailed
            | ReactionExtentErrorKind::NonConvergence
            | ReactionExtentErrorKind::BackendFailure
    )
}

/// Finds a basin in an exactly normalized equivalent, then proves that basin
/// against the original physical extensive formulation before publication.
fn recover_resolved_pt_by_extensive_normalization<'a>(
    request: ResolvedPhaseEquilibriumRequest<'a>,
    physical_error: ReactionExtentError,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    let normalization = match &request.physical_input {
        ResolvedPhaseEquilibriumPhysicalInput::Composition => {
            ExtensiveNormalization::from_physical_moles(request.initial_composition.moles())?
        }
        ResolvedPhaseEquilibriumPhysicalInput::ElementInventory(inventory) => {
            ExtensiveNormalization::from_element_inventory(inventory)?
        }
    };
    // A near-unit inventory is not an extensive-conditioning problem. Avoid
    // turning this route into an implicit general-purpose multi-start policy.
    if !is_materially_extensive_scaled(normalization.physical_inventory_scale()) {
        return Err(physical_error);
    }

    let trigger_kind = physical_error.kind();
    let trigger_message = physical_error.to_string();
    let layout = MultiphaseEquilibriumLayout::new(request.resolved.phase_specs().to_vec())?;
    let normalized_composition =
        normalization.normalize_initial_composition(&layout, &request.initial_composition)?;
    let normalized_seeds = request
        .multi_start_seeds
        .iter()
        .map(|seed| normalization.normalize_log_mole_guess(seed))
        .collect::<Result<Vec<_>, _>>()?;
    let normalized_basin_seed = request
        .recovery_basin_seed
        .as_ref()
        .map(|seed| normalization.normalize_log_mole_guess(seed))
        .transpose()?;
    let normalized_physical_input = match &request.physical_input {
        ResolvedPhaseEquilibriumPhysicalInput::Composition => {
            ResolvedPhaseEquilibriumPhysicalInput::Composition
        }
        ResolvedPhaseEquilibriumPhysicalInput::ElementInventory(inventory) => {
            ResolvedPhaseEquilibriumPhysicalInput::ElementInventory(
                normalization.normalize_element_inventory(inventory)?,
            )
        }
    };
    let normalized_mode = match request.solve_mode.clone() {
        PhaseEquilibriumSolveMode::FixedDeclaredPhases => {
            PhaseEquilibriumSolveMode::FixedDeclaredPhases
        }
        PhaseEquilibriumSolveMode::BoundedPhaseControl(policy) => {
            PhaseEquilibriumSolveMode::BoundedPhaseControl(
                policy.normalized_for_extensive_representation(normalization)?,
            )
        }
    };
    let normalized_options = request
        .solve_options
        .clone()
        .with_extensive_normalization_policy(ExtensiveNormalizationPolicy::Disabled)
        // Discovery is an internal numerical operation. Streaming its events
        // would expose normalized mole units through a physical API. A
        // successful physical retry retains ordinary diagnostics; an explicit
        // reconstruction retains typed recovery provenance instead.
        .with_diagnostics(EquilibriumDiagnosticsOptions::disabled())
        .with_trace_seed_policy(
            normalization.normalize_trace_seed_policy(request.solve_options.trace_seed_policy())?,
        );
    let recovery_tolerance = normalized_options
        .reproducibility_snapshot()
        .tolerance
        .min(EXTENSIVE_RECOVERY_MAX_SOLVER_TOLERANCE);
    let normalized_options = normalized_options.with_tolerance(recovery_tolerance)?;
    let normalized_request = ResolvedPhaseEquilibriumRequest {
        resolved: request.resolved,
        conditions: request.conditions,
        initial_composition: normalized_composition,
        physical_input: normalized_physical_input,
        candidate_selection: request.candidate_selection.clone(),
        model_policy: request.model_policy,
        solve_options: normalized_options,
        solve_mode: normalized_mode,
        multi_start_seeds: normalized_seeds,
        recovery_basin_seed: normalized_basin_seed,
        continuation_phase_set: request.continuation_phase_set.clone(),
    };
    let normalized_solution = match solve_resolved_pt_once(normalized_request) {
        Ok(solution) => solution,
        Err(recovery) => {
            return Err(ReactionExtentError::ExtensiveNormalizationRecoveryFailed {
                stage: "normalized basin discovery",
                physical: Box::new(physical_error),
                recovery: Box::new(recovery),
            });
        }
    };
    let discovery_backend = normalized_solution.solve_report().accepted_backend.clone();

    let mut physical_retry = request.clone();
    physical_retry.solve_options = physical_retry
        .solve_options
        .with_extensive_normalization_policy(ExtensiveNormalizationPolicy::Disabled);
    match &physical_retry.solve_mode {
        PhaseEquilibriumSolveMode::FixedDeclaredPhases => {
            let recovered_seed = normalization
                .denormalize_log_moles(normalized_solution.accepted_solution().log_moles())?;
            physical_retry.multi_start_seeds.insert(0, recovered_seed);
        }
        PhaseEquilibriumSolveMode::BoundedPhaseControl(_) => {
            // For an activation lifecycle, `restart_seed` retains the accepted
            // pre-activation coordinates for every previously active species
            // and only adds the incipient phase seed. The original physical
            // active set therefore sees the basin it needs while still
            // reproducing the activation transition itself.
            let normalized_seed = normalized_solution
                .phase_control_report()
                .and_then(|report| report.transitions.first())
                .filter(|transition| {
                    !transition.activated.is_empty() && transition.deactivated.is_empty()
                })
                .map(|transition| transition.restart_seed.as_slice())
                .unwrap_or_else(|| normalized_solution.accepted_solution().log_moles());
            physical_retry.recovery_basin_seed =
                Some(normalization.denormalize_log_moles(normalized_seed)?);
        }
    }

    match solve_resolved_pt_once(physical_retry) {
        Ok(physical_solution) => {
            let physical_retry_backend = physical_solution.solve_report().accepted_backend.clone();
            Ok(attach_extensive_normalization_recovery(
                physical_solution,
                ExtensiveNormalizationRecoveryEvidence {
                    physical_inventory_scale: normalization.physical_inventory_scale(),
                    trigger_kind,
                    trigger_message,
                    discovery_backend,
                    physical_retry_backend: Some(physical_retry_backend),
                    physical_retry_failure: None,
                    reconstructed_physical_boundary: false,
                },
                request.solve_options.diagnostics_options(),
            ))
        }
        Err(physical_retry_error) => {
            let retry_message = physical_retry_error.to_string();
            let physical_publication = (|| {
                let physical_build_report = crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
                    build_phase_equilibrium_problem_with_timing(
                        build_phase_equilibrium_request(
                            &request,
                            request.initial_composition.clone(),
                        )?,
                        request.solve_options.timing_mode(),
                    )?
                    .report()
                    .clone();
                normalized_solution
                    .reconstruct_physical_extensive_boundary(normalization, physical_build_report)
            })();
            match physical_publication {
                Ok(physical_solution) => Ok(attach_extensive_normalization_recovery(
                    physical_solution,
                    ExtensiveNormalizationRecoveryEvidence {
                        physical_inventory_scale: normalization.physical_inventory_scale(),
                        trigger_kind,
                        trigger_message,
                        discovery_backend,
                        physical_retry_backend: None,
                        physical_retry_failure: Some(retry_message),
                        reconstructed_physical_boundary: true,
                    },
                    request.solve_options.diagnostics_options(),
                )),
                Err(recovery) => Err(ReactionExtentError::ExtensiveNormalizationRecoveryFailed {
                    stage: "physical publication reconstruction",
                    physical: Box::new(physical_error),
                    recovery: Box::new(recovery),
                }),
            }
        }
    }
}

/// Publishes one recovery fact through the same opt-in retained/sink
/// diagnostics contract as the rest of the canonical workflow.
fn attach_extensive_normalization_recovery(
    solution: MultiphaseEquilibriumSolution,
    evidence: ExtensiveNormalizationRecoveryEvidence,
    diagnostics: &EquilibriumDiagnosticsOptions,
) -> MultiphaseEquilibriumSolution {
    let event = EquilibriumDiagnosticEvent::ExtensiveNormalizationRecoveryAccepted {
        physical_inventory_scale: evidence.physical_inventory_scale,
        trigger_kind: evidence.trigger_kind,
        discovery_backend: format!("{:?}", evidence.discovery_backend),
        physical_retry_backend: evidence
            .physical_retry_backend
            .as_ref()
            .map(|backend| format!("{backend:?}")),
        reconstructed_physical_boundary: evidence.reconstructed_physical_boundary,
    };
    let mut collector = solution.diagnostics_report().map_or_else(
        || EquilibriumDiagnosticsCollector::new(diagnostics.clone()),
        |report| EquilibriumDiagnosticsCollector::from_report(diagnostics.clone(), report),
    );
    collector.record(EquilibriumDiagnosticsMode::Summary, event);
    solution
        .with_diagnostics(collector.finish())
        .with_extensive_normalization_recovery(evidence)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_candidate_selection::{
        EquilibriumCandidatePhaseAssignment, EquilibriumCandidatePolicy,
        EquilibriumCandidateSelector,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
        EquilibriumDiagnosticEvent, EquilibriumDiagnosticsMode, EquilibriumDiagnosticsOptions,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
        PhaseEquilibriumInputKind, PhaseEquilibriumSeedSource,
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
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
    use std::collections::HashMap;

    fn real_h_o_multiphase_repository() -> Arc<ThermoRepository> {
        let source = ThermoData::try_default_repository()
            .expect("bundled repository must load for the real multiphase fixture");
        let records = [
            ("NASA_gas", "O2"),
            ("NASA_gas", "H2O"),
            ("NASA_cond", "H2O(L)"),
        ];
        let mut payloads = HashMap::new();
        for (library, record_key) in records {
            payloads
                .entry(library.to_string())
                .or_insert_with(HashMap::new)
                .insert(
                    record_key.to_string(),
                    source
                        .LibThermoData
                        .get(library)
                        .and_then(|records| records.get(record_key))
                        .cloned()
                        .expect("fixture record must exist in bundled JSON"),
                );
        }

        let mut elements = HashMap::new();
        for element in ["H", "O"] {
            let rows = source
                .ElementsData
                .get(element)
                .into_iter()
                .flatten()
                .filter(|row| {
                    row.len() >= 2
                        && records.iter().any(|(library, record_key)| {
                            row[0] == *record_key && row[1] == *library
                        })
                })
                .cloned()
                .collect();
            elements.insert(element.to_string(), rows);
        }

        Arc::new(ThermoRepository::from_parts(
            records
                .iter()
                .map(|(library, record_key)| (library.to_string(), record_key.to_string()))
                .collect(),
            payloads,
            elements,
            vec!["NASA_gas".to_string(), "NASA_cond".to_string()],
            HashMap::new(),
            HashMap::from([
                ("NASA_gas".to_string(), "NASA".to_string()),
                ("NASA_cond".to_string(), "NASA".to_string()),
            ]),
            vec!["NASA_gas".to_string(), "NASA_cond".to_string()],
            Vec::new(),
        ))
    }

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
        let summary = crate::Thermodynamics::ChemEquilibrium::
            equilibrium_diagnostics_display::format_solution_execution_summary(&solution);
        assert!(summary.contains("equilibrium execution:"));
        assert!(summary.contains("phase control:"));
        assert!(summary.contains("verbose diagnostics: disabled"));
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
    fn element_inventory_pipeline_solves_without_dense_component_moles() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "H2".to_string(),
            "O2".to_string(),
            "H2O".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let outcome = PhaseEquilibriumPipelineRequest::from_element_inventory(
            spec,
            ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)]).unwrap(),
            EquilibriumConditions::new(1_200.0, 101_325.0, 101_325.0).unwrap(),
        )
        .solve()
        .expect("element-defined pipeline must resolve and solve locally");

        assert_eq!(outcome.resolved().phase_specs()[0].components().len(), 3);
        assert!(
            outcome
                .solution()
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite() && *moles >= 0.0)
        );
    }

    #[test]
    fn explicit_species_universe_changes_the_problem_without_catalog_supplementation() {
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let conditions = EquilibriumConditions::new(1_200.0, 101_325.0, 101_325.0)
            .expect("conditions must validate");
        let full_spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "H2".to_string(),
            "O2".to_string(),
            "H2O".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("full explicit universe must validate");
        let restricted_spec =
            SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
                "H2".to_string(),
                "O2".to_string(),
            ]))
            .with_library_priorities(vec!["NASA_gas".to_string()])
            .with_search_in_nist(false)
            .build()
            .expect("restricted explicit universe must validate");

        let full = PhaseEquilibriumPipelineRequest::from_element_inventory(
            full_spec,
            inventory.clone(),
            conditions,
        )
        .solve()
        .expect("full explicit universe must solve");
        let restricted = PhaseEquilibriumPipelineRequest::from_element_inventory(
            restricted_spec,
            inventory,
            conditions,
        )
        .solve()
        .expect("restricted explicit universe must solve");

        assert_eq!(
            full.resolved().phase_specs()[0].components(),
            ["H2", "O2", "H2O"]
        );
        assert_eq!(
            restricted.resolved().phase_specs()[0].components(),
            ["H2", "O2"]
        );
        assert_eq!(restricted.solution().component_moles().len(), 2);
        assert!(
            restricted
                .solution()
                .build_report()
                .components()
                .iter()
                .all(|component| component.component().label() != "H2O")
        );
        assert!(full.solution().component_moles()[2] > 1.0e-12);
        assert!(
            restricted
                .solution()
                .accepted_solution()
                .validation()
                .max_abs_element_balance_error
                < 1.0e-8
        );
    }

    #[test]
    fn truncated_catalog_selection_cannot_bypass_inventory_representability() {
        let repository = ThermoData::try_default_repository()
            .expect("bundled repository must load for truncation validation");
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let policy = EquilibriumCandidatePolicy::default()
            .with_library_preference(vec!["NASA_gas".to_string()])
            .with_temperature_range(1_200.0, 1_200.0)
            .expect("candidate temperature range must validate")
            .with_max_candidates(1)
            .expect("candidate limit must validate");
        let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
            .select_inventory(&inventory, policy)
            .expect("truncated H/O selection must still produce evidence");
        assert!(selection.is_truncated());
        assert_eq!(selection.selected().len(), 1);

        let phase_plan = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::ideal_gas(
                PhaseId::new(Some("gas".to_string())),
                selection
                    .selected()
                    .iter()
                    .map(|candidate| candidate.record_key().to_string())
                    .collect(),
            ),
        ]);
        let error =
            PhaseEquilibriumPipelineRequest::from_candidate_selection_with_element_inventory(
                repository,
                &selection,
                &phase_plan,
                inventory,
                EquilibriumConditions::new(1_200.0, 101_325.0, 101_325.0).unwrap(),
            )
            .expect("truncated selection must still produce a typed phase plan")
            .solve()
            .expect_err("a truncated universe that cannot represent b must be rejected");

        assert!(matches!(
            error,
            PhaseEquilibriumPipelineError::Solve(ReactionExtentError::Preparation(_))
        ));
    }

    #[test]
    fn candidate_selection_provenance_survives_bridge_and_reproducibility() {
        let repository = ThermoData::try_default_repository()
            .expect("bundled repository must load for selection provenance");
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let policy = EquilibriumCandidatePolicy::default()
            .with_library_preference(vec!["NASA_gas".to_string()])
            .with_temperature_range(1_200.0, 1_200.0)
            .expect("temperature selection range must validate")
            .with_max_candidates(3)
            .expect("candidate cap must validate");
        let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
            .select_inventory(&inventory, policy)
            .expect("local H/O selection must succeed");
        assert_eq!(selection.selected().len(), 3);
        let phase_plan = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::ideal_gas(
                PhaseId::new(Some("gas".to_string())),
                selection
                    .selected()
                    .iter()
                    .map(|candidate| candidate.record_key().to_string())
                    .collect(),
            ),
        ]);
        let outcome =
            PhaseEquilibriumPipelineRequest::from_candidate_selection_with_element_inventory(
                repository,
                &selection,
                &phase_plan,
                inventory,
                EquilibriumConditions::new(1_200.0, 101_325.0, 101_325.0).unwrap(),
            )
            .expect("candidate phase plan must build")
            .with_solve_options(
                EquilibriumSolveOptions::new()
                    .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)))
                    .expect("single local backend policy must validate"),
            )
            .solve()
            .expect("selected local H/O universe must solve");

        assert_eq!(
            outcome.solution().build_report().candidate_selection(),
            Some(&selection)
        );
        let capsule = crate::Thermodynamics::ChemEquilibrium::equilibrium_reproducibility::
            EquilibriumReproducibilityCapsule::from_outcome(
                &outcome,
                &EquilibriumSolveOptions::new(),
            );
        assert!(capsule.candidate_selection.is_some());
        assert_eq!(
            capsule
                .candidate_selection
                .as_ref()
                .unwrap()
                .selected_records
                .len(),
            selection.selected().len()
        );
        assert_eq!(
            capsule
                .candidate_selection
                .as_ref()
                .unwrap()
                .rejected_records
                .len(),
            selection.rejected().len()
        );
        assert!(capsule.candidate_selection.as_ref().unwrap().truncated);
    }

    #[test]
    fn real_offline_multiphase_fixture_tracks_phase_lifecycle_inventory_and_universe_provenance() {
        let source = ThermoData::try_default_repository()
            .expect("bundled repository must load for the real multiphase fixture");
        let source_payloads = source.LibThermoData.as_ref().clone();
        let source_elements = source.ElementsData.as_ref().clone();
        let repository = real_h_o_multiphase_repository();
        let inventory = ElementInventory::from_amounts([("H", 1.0), ("O", 1.0)])
            .expect("closed H/O inventory must validate");
        let policy = EquilibriumCandidatePolicy::default()
            .with_library_preference(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
            .with_temperature_range(300.0, 3_000.0)
            .expect("real fixture range must validate")
            .with_max_candidates(4)
            .expect("real fixture candidate cap must validate");
        let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
            .select_inventory(&inventory, policy)
            .expect("real H/O fixture selection must succeed");

        assert_eq!(
            selection
                .selected()
                .iter()
                .map(|candidate| candidate.record_key())
                .collect::<Vec<_>>(),
            vec!["H2O", "O2", "H2O(L)"]
        );
        assert_eq!(
            selection.selected()[2].physical_state(),
            Some(PhysicalState::Liquid)
        );

        let phase_plan = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::ideal_gas(
                PhaseId::new(Some("gas".to_string())),
                vec!["H2O".to_string(), "O2".to_string()],
            ),
            EquilibriumCandidatePhaseAssignment::pure_condensed(
                PhaseId::new(Some("liquid".to_string())),
                PhysicalState::Liquid,
                vec!["H2O(L)".to_string()],
            ),
        ]);

        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            phase_plan
                .build_spec(&selection)
                .expect("real phase plan must build a typed specification"),
            Arc::clone(&repository),
        )
        .expect("real offline phase universe must resolve");
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
            .expect("real offline phase layout must validate");
        let solve_automatic = |temperature, initial_moles| {
            let initial_composition =
                MultiphaseInitialComposition::from_dense(&layout, initial_moles)
                    .expect("element-conserving numerical seed must match the layout");
            solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    &resolved,
                    EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
                    initial_composition,
                )
                .with_element_inventory(inventory.clone())
                .with_candidate_selection(Some(selection.clone()))
                .with_solve_options(
                    EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
                )
                .with_phase_control_policy(PhaseControlPolicy::default()),
            )
            .expect("real offline multiphase fixture must solve")
        };
        let low_temperature = solve_automatic(350.0, vec![0.5, 0.25, 0.0]);
        let high_temperature = solve_automatic(550.0, vec![0.25, 0.25, 0.25]);

        for outcome in [&low_temperature, &high_temperature] {
            let report = outcome.build_report();
            assert_eq!(
                report.input_kind(),
                PhaseEquilibriumInputKind::ElementInventory
            );
            assert_eq!(
                report.seed_evidence().source(),
                PhaseEquilibriumSeedSource::SuppliedNumericalSeed
            );
            assert_eq!(report.candidate_selection(), Some(&selection));
            let totals = report
                .element_labels()
                .iter()
                .zip(report.element_totals().iter().copied())
                .map(|(label, total)| (label.as_str(), total))
                .collect::<HashMap<_, _>>();
            assert_eq!(totals.get("H"), Some(&1.0));
            assert_eq!(totals.get("O"), Some(&1.0));
            assert!(
                outcome
                    .acceptance_report()
                    .expect("bounded solve must publish acceptance evidence")
                    .final_validation
                    .max_abs_element_balance_error
                    < 1e-6
            );
        }

        let low_liquid_status = low_temperature
            .phase_status(&PhaseId::new(Some("liquid".to_string())))
            .expect("liquid status must be published");
        assert!(matches!(
            low_liquid_status,
            PhaseStatus::Appeared | PhaseStatus::Active
        ));
        assert!(
            low_temperature
                .phase_control_report()
                .expect("low-temperature phase control report must be published")
                .transitions
                .iter()
                .any(|transition| transition.activated.iter().any(|phase| phase.index() == 1))
        );

        let high_liquid_status = high_temperature
            .phase_status(&PhaseId::new(Some("liquid".to_string())))
            .expect("liquid status must be published");
        assert!(matches!(
            high_liquid_status,
            PhaseStatus::Inactive | PhaseStatus::Disappeared
        ));
        assert!(
            high_temperature
                .phase_control_report()
                .expect("high-temperature phase control report must be published")
                .transitions
                .iter()
                .any(|transition| transition
                    .deactivated
                    .iter()
                    .any(|phase| phase.index() == 1))
        );

        let explicit_initial =
            MultiphaseInitialComposition::from_dense(&layout, vec![0.5, 0.25, 0.0])
                .expect("explicit numerical seed must match the resolved layout");
        let explicit = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                &resolved,
                EquilibriumConditions::new(350.0, 101_325.0, 101_325.0).unwrap(),
                explicit_initial,
            )
            .with_element_inventory(inventory)
            .with_solve_options(
                EquilibriumSolveOptions::new().with_timing_mode(EquilibriumTimingMode::Enabled),
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )
        .expect("explicit real phase universe must solve");
        assert!(explicit.build_report().candidate_selection().is_none());
        assert_eq!(
            explicit.build_report().input_kind(),
            PhaseEquilibriumInputKind::ElementInventory
        );
        assert_eq!(
            explicit.phases(),
            low_temperature.phases(),
            "automatic selection and explicit species declaration must resolve the same phases"
        );
        for (automatic, explicit) in low_temperature
            .component_moles()
            .iter()
            .zip(explicit.component_moles())
        {
            assert!((automatic - explicit).abs() <= 1e-8);
        }

        let build_report = low_temperature.build_report();
        let acceptance = low_temperature
            .acceptance_report()
            .expect("release story must retain acceptance evidence");
        println!(
            "offline automatic multiphase release story\n  b={:?}\n  selected={:?}\n  rejected={:?}\n  matrix={}x{} rank={} reaction_dimension={}\n  seed={:?} backend={:?} iterations={} residual={:.3e} balance={:.3e}\n  transitions={} provenance={:?}\n  timing={:?}\n  explicit_backend={:?} explicit_iterations={}",
            build_report.element_totals(),
            selection
                .selected()
                .iter()
                .map(|candidate| format!("{}:{}", candidate.library(), candidate.record_key()))
                .collect::<Vec<_>>(),
            selection
                .rejected()
                .iter()
                .map(|rejection| format!("{}:{}", rejection.library(), rejection.substance()))
                .collect::<Vec<_>>(),
            build_report.components().len(),
            build_report.element_labels().len(),
            build_report.solver_element_labels().len(),
            build_report
                .components()
                .len()
                .saturating_sub(build_report.solver_element_labels().len()),
            build_report.seed_evidence().source(),
            low_temperature.solve_report().accepted_backend,
            low_temperature.nonlinear_iterations(),
            acceptance.final_validation.residual_l2_norm,
            acceptance.final_validation.max_abs_element_balance_error,
            low_temperature.phase_control_transitions(),
            build_report
                .components()
                .iter()
                .map(|component| {
                    format!(
                        "{}<-{}:{}",
                        component.component().label(),
                        component.thermo_source().library(),
                        component.thermo_source().record_key()
                    )
                })
                .collect::<Vec<_>>(),
            low_temperature.timing_report(),
            explicit.solve_report().accepted_backend,
            explicit.nonlinear_iterations(),
        );

        assert_eq!(source.LibThermoData.as_ref(), &source_payloads);
        assert_eq!(source.ElementsData.as_ref(), &source_elements);
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
    fn normalized_phase_policy_scales_only_the_absolute_phase_epsilon() {
        let policy = PhaseControlPolicy::with_explicit_hysteresis(1.0e-12, -2.0, -1.0)
            .unwrap()
            .normalized_for_extensive_representation(
                ExtensiveNormalization::from_physical_inventory_scale(1.0e6).unwrap(),
            )
            .unwrap();
        let manager = policy.into_phase_manager();
        assert!((manager.phase_eps - 1.0e-18).abs() <= 1.0e-30);
        assert_eq!(
            manager.phase_hysteresis,
            crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::
                PhaseHysteresisPolicy::Explicit {
                dg_create: -2.0,
                dg_keep: -1.0,
            }
        );
    }

    #[test]
    fn extensive_normalization_policy_is_explicit_and_reproducible() {
        let defaults = EquilibriumSolveOptions::default();
        assert_eq!(
            defaults.extensive_normalization_policy(),
            ExtensiveNormalizationPolicy::OnNumericalFailure
        );
        assert_eq!(
            defaults
                .reproducibility_snapshot()
                .extensive_normalization_policy,
            "OnNumericalFailure"
        );

        let disabled =
            defaults.with_extensive_normalization_policy(ExtensiveNormalizationPolicy::Disabled);
        assert_eq!(
            disabled.extensive_normalization_policy(),
            ExtensiveNormalizationPolicy::Disabled
        );
        assert_eq!(
            disabled
                .reproducibility_snapshot()
                .extensive_normalization_policy,
            "Disabled"
        );
    }

    #[test]
    fn extensive_recovery_guard_is_symmetric_and_has_an_explicit_boundary() {
        assert!(!is_materially_extensive_scaled(0.1));
        assert!(!is_materially_extensive_scaled(1.0));
        assert!(!is_materially_extensive_scaled(10.0));
        assert!(is_materially_extensive_scaled(0.099));
        assert!(is_materially_extensive_scaled(10.001));
        assert!(!is_materially_extensive_scaled(0.0));
        assert!(!is_materially_extensive_scaled(f64::NAN));
        assert!(!is_materially_extensive_scaled(f64::INFINITY));
    }

    #[test]
    fn extensive_recovery_authorizes_only_retryable_numerical_failures() {
        let retryable = [
            ReactionExtentError::AllBackendsFailed { attempts: vec![] },
            ReactionExtentError::SolveError(SolveError::MaxIterations),
            ReactionExtentError::SolveError(SolveError::SingularMatrix),
            ReactionExtentError::BackendFailure {
                backend: "test".to_string(),
                kind: BackendFailureKind::NumericalBreakdown,
                message: "test numerical breakdown".to_string(),
            },
        ];
        assert!(
            retryable
                .iter()
                .all(is_extensive_normalization_recovery_candidate)
        );

        let non_retryable = [
            ReactionExtentError::Cancelled,
            ReactionExtentError::InvalidProblem {
                field: "test",
                message: "invalid".to_string(),
            },
            ReactionExtentError::InvalidConditions {
                parameter: "temperature",
                value: 0.0,
            },
            ReactionExtentError::InvalidCandidate {
                field: "residual",
                message: "rejected".to_string(),
            },
            ReactionExtentError::ResidualEvaluation("test residual".to_string()),
            ReactionExtentError::JacobianEvaluation("test jacobian".to_string()),
            ReactionExtentError::UnsupportedBackendCapability {
                backend: "test".to_string(),
                capability: "test",
                alternatives: "none".to_string(),
            },
            ReactionExtentError::CascadeAborted {
                attempts: vec![],
                cause: Box::new(ReactionExtentError::InvalidProblem {
                    field: "test",
                    message: "abort".to_string(),
                }),
            },
        ];
        assert!(
            non_retryable
                .iter()
                .all(|error| !is_extensive_normalization_recovery_candidate(error))
        );
    }

    #[test]
    fn options_snapshot_round_trips_the_extensive_recovery_policy() {
        for policy in [
            ExtensiveNormalizationPolicy::Disabled,
            ExtensiveNormalizationPolicy::OnNumericalFailure,
        ] {
            let snapshot = EquilibriumSolveOptions::new()
                .with_extensive_normalization_policy(policy)
                .reproducibility_snapshot();
            let json = serde_json::to_string(&snapshot).unwrap();
            let decoded: EquilibriumSolveOptionsSnapshot = serde_json::from_str(&json).unwrap();
            assert_eq!(decoded, snapshot);
            assert_eq!(
                decoded.extensive_normalization_policy,
                format!("{policy:?}")
            );
        }
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
