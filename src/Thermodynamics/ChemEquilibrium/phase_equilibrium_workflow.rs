//! Narrow public fixed-`P,T` phase-equilibrium workflow.
//!
//! This facade owns orchestration only: it joins validated resolved data,
//! physical inventory, numerical settings, bridge construction, and immutable
//! result publication. It does not duplicate residual construction or expose
//! the historical mutable solver as an alternative public engine.

use std::fmt;
use std::sync::Arc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_candidate_selection::{
    EquilibriumCandidatePhasePlan, EquilibriumCandidateSelectionReport,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumSolverSettings, Solvers,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverBackend, SolverPolicy,
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
        })
    }

    /// Consumes the wrapper for the crate-internal immutable runner.
    pub(crate) fn into_settings(self) -> EquilibriumSolverSettings {
        self.settings
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

    /// Returns the timing policy carried by this solve request.
    pub fn timing_mode(&self) -> EquilibriumTimingMode {
        self.timing_mode
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
            repository,
        } = self;

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
            repository,
        } = self;
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
    let started = std::time::Instant::now();
    match request.solve_mode {
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
            bundle
                .solve_with(|configured| *configured = settings)
                .and_then(|bundle| bundle.into_multiphase_solution())
                .map(|solution| solution.with_timing_total(started.elapsed()))
        }
        PhaseEquilibriumSolveMode::BoundedPhaseControl(phase_control_policy) => {
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
                .solve_with_bounded_phase_control(
                    |configured| *configured = settings,
                    |configured| *configured = phase_control_policy.into_phase_manager(),
                )
                .map(|solution| solution.with_timing_total(started.elapsed()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::prelude::{
        LegacyEquilibriumSolver, RustedSciTheSolver, SolverBackend, SolverPolicy,
    };
    use crate::Thermodynamics::User_PhaseOrSolution::{
        SubstanceSystemSpecBuilder, SubstancesContainer,
    };
    use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};

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
