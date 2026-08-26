//! Conversion boundary from the equilibrium GUI model to the production API.
//!
//! View code must stop at [`ValidatedEquilibriumGuiConfig`]. This module is the
//! only place that knows how GUI enums map to `PhaseSpec`, `PhaseComponentId`,
//! solver policies, and the fixed-pressure equilibrium workflow.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::EquilibriumDiagnosticEvent;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::EquilibriumExecutionControl;
use crate::Thermodynamics::ChemEquilibrium::prelude::{
    CandidateSelectionError, ElementSearchMode, EquilibriumCandidatePolicy,
    EquilibriumCandidateSelectionReport, EquilibriumCandidateSelector, EquilibriumConditions,
    EquilibriumConstantValidationMode, EquilibriumConstraint, EquilibriumDiagnosticsMode,
    EquilibriumDiagnosticsOptions, EquilibriumRangeDiagnosticsPolicy, EquilibriumSolveOptions,
    EquilibriumTimingMode, FixedPressureEnthalpySolution, LegacyEquilibriumSolver,
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition, PhSolveMode, PhaseComponentId,
    PhaseControlPolicy, PhaseEquilibriumPipelineError, PhaseEquilibriumPipelineRequest, PhaseId,
    PhaseModel, PhaseSpec, ResolvedPhaseEnthalpyRequest, ResolvedPhaseEquilibriumOutcome,
    ResolvedThermochemistry, RustedSciTheSolver, SolverBackend, SolverCascadeBudget, SolverPolicy,
    SubstanceSystemFactory, SubstanceSystemFactoryError, SubstanceSystemSpec, TemperatureBounds,
    TemperatureGrid, TemperatureRangeSolution, ThermoRepository, TraceSpeciesSeedPolicy,
};
use crate::gui::equilibrium_gui_model::{
    EquilibriumSolverDraft, GuiKeqValidationMode, GuiPhaseLifecycleTrace, GuiPhaseModel,
    GuiPhysicalState, GuiRangeLifecycleTrace, GuiSolverBackend, ValidatedCandidatePolicy,
    ValidatedEquilibriumGuiConfig, ValidatedInventory, ValidatedLookup, ValidatedPhaseMode,
    ValidatedProblem, ValidatedSolver, ValidatedTemperature, ValidatedTraceSeedPolicy,
};
use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

/// A validated request ready for the production point or range facade.
pub enum EquilibriumGuiSolveRequest {
    /// One fixed-pressure, fixed-temperature solve.
    Point(PhaseEquilibriumPipelineRequest),
    /// One resolved policy plus a typed continuation grid.
    Range {
        request: PhaseEquilibriumPipelineRequest,
        temperatures: TemperatureGrid,
    },
    /// Fixed-pressure, fixed-total-enthalpy request. Resolution is delayed
    /// until the worker executes it, just like the ordinary PT request.
    Ph(EquilibriumGuiPhRequest),
}

/// GUI-owned P,H request payload.
///
/// Keeping the unresolved `SubstanceSystemSpec` here is important: preparing
/// a document remains a pure editor operation, while the worker performs the
/// repository lookup and captures one immutable thermochemistry bundle.
pub struct EquilibriumGuiPhRequest {
    spec: SubstanceSystemSpec,
    sparse_initial_moles: Vec<(PhaseComponentId, f64)>,
    pressure_pa: f64,
    reference_pressure_pa: f64,
    target_enthalpy_j: f64,
    temperature_bounds: TemperatureBounds,
    seed_temperature_k: f64,
    options: EquilibriumSolveOptions,
    phase_control_policy: Option<PhaseControlPolicy>,
    repository: Option<Arc<ThermoRepository>>,
}

impl EquilibriumGuiPhRequest {
    fn with_execution_control(mut self, control: EquilibriumExecutionControl) -> Self {
        self.options = self.options.with_execution_control(control);
        self
    }

    fn with_diagnostic_sink<F>(mut self, sink: F) -> Self
    where
        F: Fn(EquilibriumDiagnosticEvent) + Send + Sync + 'static,
    {
        let diagnostics = self.options.diagnostics_options().clone().with_sink(sink);
        self.options = self.options.with_diagnostics(diagnostics);
        self
    }

    fn solve(self) -> Result<FixedPressureEnthalpySolution, PhaseEquilibriumPipelineError> {
        let Self {
            spec,
            sparse_initial_moles,
            pressure_pa,
            reference_pressure_pa,
            target_enthalpy_j,
            temperature_bounds,
            seed_temperature_k,
            options,
            phase_control_policy,
            repository,
        } = self;
        let resolved = match repository {
            Some(repository) => {
                SubstanceSystemFactory::resolve_phase_system_with_repository(spec, repository)?
            }
            None => SubstanceSystemFactory::resolve_phase_system(spec)?,
        };
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
        let composition = MultiphaseInitialComposition::from_sparse(&layout, sparse_initial_moles)?;
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
        let constraint = EquilibriumConstraint::ph(
            pressure_pa,
            reference_pressure_pa,
            target_enthalpy_j,
            seed_temperature_k,
        )?;
        let ph_solve_mode = match thermochemistry.symbolic_for_bounds(temperature_bounds) {
            Ok(Some(_)) => PhSolveMode::Monolithic,
            // A broad bracket may cross a native NASA/NIST coefficient switch.
            // The GUI must keep such a valid numeric request usable, while
            // never pretending that RST can differentiate a piecewise source
            // as one smooth symbolic expression.
            Ok(None) | Err(_) => PhSolveMode::NestedTemperature,
        };
        let mut request = ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            composition,
            constraint,
            temperature_bounds,
            thermochemistry,
        )?
        // Prefer the coupled route when the resolved data expose one exact
        // symbolic coefficient interval. A bracket crossing a native source
        // boundary remains a valid numeric request, but uses the safeguarded
        // nested reference route instead of hiding a symbolic capability
        // error from the GUI user.
        .with_ph_solve_mode(ph_solve_mode)
        .with_solve_options(options);
        if let Some(policy) = phase_control_policy {
            request = request.with_phase_control_policy(policy);
        }
        crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::solve_resolved_ph(request)
            .map_err(Into::into)
    }
}

impl EquilibriumGuiSolveRequest {
    /// Attaches a worker-owned execution handle without changing the request
    /// payload or its document fingerprint.
    pub fn with_execution_control(self, control: EquilibriumExecutionControl) -> Self {
        match self {
            Self::Point(request) => {
                let options = request
                    .solve_options()
                    .clone()
                    .with_execution_control(control);
                Self::Point(request.with_solve_options(options))
            }
            Self::Range {
                request,
                temperatures,
            } => {
                let options = request
                    .solve_options()
                    .clone()
                    .with_execution_control(control);
                Self::Range {
                    request: request.with_solve_options(options),
                    temperatures,
                }
            }
            Self::Ph(request) => Self::Ph(request.with_execution_control(control)),
        }
    }

    /// Attaches an observational lifecycle stream owned by the GUI worker.
    ///
    /// The engine still decides whether diagnostics are enabled from the
    /// validated document. Supplying a sink alone cannot turn tracing on or
    /// alter the numerical transaction.
    pub fn with_diagnostic_sink<F>(self, sink: F) -> Self
    where
        F: Fn(EquilibriumDiagnosticEvent) + Send + Sync + 'static,
    {
        let sink: Arc<dyn Fn(EquilibriumDiagnosticEvent) + Send + Sync> = Arc::new(sink);
        match self {
            Self::Point(request) => {
                let diagnostics = request
                    .solve_options()
                    .diagnostics_options()
                    .clone()
                    .with_sink({
                        let sink = Arc::clone(&sink);
                        move |event| sink(event)
                    });
                let options = request
                    .solve_options()
                    .clone()
                    .with_diagnostics(diagnostics);
                Self::Point(request.with_solve_options(options))
            }
            Self::Range {
                request,
                temperatures,
            } => {
                let diagnostics = request
                    .solve_options()
                    .diagnostics_options()
                    .clone()
                    .with_sink({
                        let sink = Arc::clone(&sink);
                        move |event| sink(event)
                    });
                let options = request
                    .solve_options()
                    .clone()
                    .with_diagnostics(diagnostics);
                Self::Range {
                    request: request.with_solve_options(options),
                    temperatures,
                }
            }
            Self::Ph(request) => Self::Ph(request.with_diagnostic_sink(move |event| sink(event))),
        }
    }

    /// Executes through the canonical production orchestration boundary.
    pub fn solve(self) -> Result<EquilibriumGuiSolveOutcome, PhaseEquilibriumPipelineError> {
        match self {
            Self::Point(request) => request.solve().map(EquilibriumGuiSolveOutcome::Point),
            Self::Range {
                request,
                temperatures,
            } => request
                .solve_temperature_range(temperatures)
                .map(EquilibriumGuiSolveOutcome::Range),
            Self::Ph(request) => request.solve().map(EquilibriumGuiSolveOutcome::Ph),
        }
    }

    /// Whether this request is a temperature continuation.
    pub const fn is_range(&self) -> bool {
        matches!(self, Self::Range { .. })
    }
}

/// Immutable production outcome owned by the GUI worker/result layer.
#[derive(Debug, Clone)]
pub enum EquilibriumGuiSolveOutcome {
    Point(ResolvedPhaseEquilibriumOutcome),
    Range(TemperatureRangeSolution),
    Ph(FixedPressureEnthalpySolution),
}

/// Errors specific to converting an already validated GUI model into an
/// engine request. User-editable errors belong to the model validation report.
#[derive(Debug)]
pub enum EquilibriumGuiRequestError {
    /// Kept as an explicit guard if a future caller bypasses run validation.
    FeatureUnavailable(&'static str),
    /// The phase declaration was incompatible with the production phase API.
    PhaseSpecification(SubstanceSystemFactoryError),
    /// The production request rejected a typed condition or policy.
    EngineInput(crate::Thermodynamics::ChemEquilibrium::prelude::ReactionExtentError),
    /// Candidate discovery failed before a solver request could be built.
    CandidateSelection(CandidateSelectionError),
}

impl fmt::Display for EquilibriumGuiRequestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FeatureUnavailable(feature) => write!(f, "feature is unavailable: {feature}"),
            Self::PhaseSpecification(error) => write!(f, "phase specification failed: {error}"),
            Self::EngineInput(error) => write!(f, "equilibrium request failed validation: {error}"),
            Self::CandidateSelection(error) => write!(f, "candidate selection failed: {error}"),
        }
    }
}

impl std::error::Error for EquilibriumGuiRequestError {}

impl From<SubstanceSystemFactoryError> for EquilibriumGuiRequestError {
    fn from(value: SubstanceSystemFactoryError) -> Self {
        Self::PhaseSpecification(value)
    }
}

impl From<crate::Thermodynamics::ChemEquilibrium::prelude::ReactionExtentError>
    for EquilibriumGuiRequestError
{
    fn from(value: crate::Thermodynamics::ChemEquilibrium::prelude::ReactionExtentError) -> Self {
        Self::EngineInput(value)
    }
}

impl From<CandidateSelectionError> for EquilibriumGuiRequestError {
    fn from(value: CandidateSelectionError) -> Self {
        Self::CandidateSelection(value)
    }
}

/// Runs the engine's deterministic element candidate selector for a validated
/// GUI document. The returned report is immutable derived state and does not
/// mutate the repository or create a solver request.
pub fn select_equilibrium_candidates(
    config: &ValidatedEquilibriumGuiConfig,
    repository: Arc<ThermoRepository>,
) -> Result<EquilibriumCandidateSelectionReport, EquilibriumGuiRequestError> {
    let ValidatedInventory::ElementCandidates {
        elements,
        candidate_policy,
        ..
    } = &config.inventory
    else {
        return Err(EquilibriumGuiRequestError::FeatureUnavailable(
            "candidate preview requires element inventory mode",
        ));
    };
    let policy = build_candidate_policy(candidate_policy, &config.lookup)?;
    EquilibriumCandidateSelector::new(repository)
        .select(elements, policy)
        .map_err(Into::into)
}

fn build_candidate_policy(
    candidate: &ValidatedCandidatePolicy,
    lookup: &ValidatedLookup,
) -> Result<EquilibriumCandidatePolicy, EquilibriumGuiRequestError> {
    let mode = match candidate.element_mode {
        crate::gui::equilibrium_gui_model::GuiElementSearchMode::Exact => {
            ElementSearchMode::ExactSet
        }
        crate::gui::equilibrium_gui_model::GuiElementSearchMode::SubsetOf => {
            ElementSearchMode::SubsetOf
        }
    };
    let mut policy = EquilibriumCandidatePolicy::new(mode);
    if !candidate.physical_states.is_empty() {
        policy = policy.with_physical_states(
            candidate
                .physical_states
                .iter()
                .copied()
                .map(map_physical_state)
                .collect(),
        );
    }
    policy = policy
        .with_temperature_range(candidate.temperature_lower_k, candidate.temperature_upper_k)?
        .with_max_candidates(candidate.max_candidates)?;
    if let ValidatedLookup::Explicit {
        priority_libraries,
        permitted_libraries,
        ..
    } = lookup
    {
        // Candidate selection has one ordered preference list, while the GUI
        // exposes the engine's two-level lookup policy.  When permitted
        // libraries are present they form the closed candidate universe;
        // priorities only control its order. This prevents a preview from
        // displaying records that the subsequent resolver is forbidden to
        // use.
        let mut libraries = Vec::new();
        if permitted_libraries.is_empty() {
            libraries.extend(priority_libraries.iter().cloned());
        } else {
            libraries.extend(
                priority_libraries
                    .iter()
                    .filter(|library| permitted_libraries.contains(*library))
                    .cloned(),
            );
            for library in permitted_libraries {
                if !libraries.contains(library) {
                    libraries.push(library.clone());
                }
            }
        }
        policy = policy.with_library_preference(libraries);
    }
    Ok(policy)
}

/// Builds the canonical production request from a validated GUI config.
///
/// No repository lookup occurs here. The request retains an optional shared
/// repository handle and resolves transactionally only when the worker calls
/// `solve()`.
pub fn build_equilibrium_request(
    config: ValidatedEquilibriumGuiConfig,
    repository: Option<Arc<ThermoRepository>>,
) -> Result<EquilibriumGuiSolveRequest, EquilibriumGuiRequestError> {
    let ValidatedEquilibriumGuiConfig {
        problem,
        inventory,
        lookup,
        phase_mode,
        solver,
        diagnostics,
        ..
    } = config;

    let (phases, sparse_initial_moles, component_libraries) = phases_and_inventory(&inventory)?;
    let mut spec = SubstanceSystemSpec::from_phases(phases)?;
    apply_lookup_policy(&mut spec, lookup, component_libraries);

    if let ValidatedProblem::FixedPh {
        pressure_pa,
        reference_pressure_pa,
        target_enthalpy_j,
        temperature_bounds,
    } = problem
    {
        let policy = match phase_mode {
            ValidatedPhaseMode::FixedDeclared => None,
            ValidatedPhaseMode::Bounded {
                phase_epsilon,
                dg_create,
                dg_keep,
                max_phase_iterations,
            } => Some(
                PhaseControlPolicy::with_explicit_hysteresis(phase_epsilon, dg_create, dg_keep)?
                    .with_max_phase_iterations(max_phase_iterations)?,
            ),
        };
        return Ok(EquilibriumGuiSolveRequest::Ph(EquilibriumGuiPhRequest {
            spec,
            sparse_initial_moles,
            pressure_pa,
            reference_pressure_pa,
            target_enthalpy_j,
            temperature_bounds: TemperatureBounds::new(
                temperature_bounds.lower_k,
                temperature_bounds.upper_k,
            )?,
            seed_temperature_k: temperature_bounds.seed_k,
            options: build_solve_options(&solver, &diagnostics)?,
            phase_control_policy: policy,
            repository,
        }));
    }

    let ValidatedProblem::FixedPt {
        pressure_pa,
        reference_pressure_pa,
        temperature,
    } = problem
    else {
        unreachable!("validated equilibrium problem has an unsupported variant")
    };

    let first_temperature = match &temperature {
        ValidatedTemperature::Point(value) => *value,
        ValidatedTemperature::Range { start_k, .. } => *start_k,
    };
    let conditions =
        EquilibriumConditions::new(first_temperature, pressure_pa, reference_pressure_pa)?;
    let mut request = PhaseEquilibriumPipelineRequest::new_with_sparse_initial_composition(
        spec,
        sparse_initial_moles,
        conditions,
    );
    if let Some(repository) = repository {
        request = request.with_repository(repository);
    }

    let options = build_solve_options(&solver, &diagnostics)?;
    request = request.with_solve_options(options);
    request = apply_phase_mode(request, phase_mode)?;

    match temperature {
        ValidatedTemperature::Point(_) => Ok(EquilibriumGuiSolveRequest::Point(request)),
        ValidatedTemperature::Range {
            start_k,
            end_k,
            point_count,
        } => {
            let values = linspace(start_k, end_k, point_count);
            let grid = TemperatureGrid::new(values)?;
            Ok(EquilibriumGuiSolveRequest::Range {
                request,
                temperatures: grid,
            })
        }
    }
}

fn phases_and_inventory(
    inventory: &ValidatedInventory,
) -> Result<
    (
        Vec<PhaseSpec>,
        Vec<(PhaseComponentId, f64)>,
        HashMap<String, String>,
    ),
    EquilibriumGuiRequestError,
> {
    let phases = match inventory {
        ValidatedInventory::ExplicitPhases(phases) => phases,
        ValidatedInventory::ElementCandidates { assignments, .. } => assignments,
    };
    if phases.is_empty() {
        return Err(EquilibriumGuiRequestError::FeatureUnavailable(
            "element candidate preview has no confirmed phase assignments",
        ));
    }

    let mut specs = Vec::with_capacity(phases.len());
    let mut initial = Vec::new();
    let mut component_libraries = HashMap::new();
    for phase in phases {
        let phase_id = PhaseId::new(phase.id.clone());
        let components = phase
            .components
            .iter()
            .map(|component| component.substance.clone())
            .collect::<Vec<_>>();
        let spec = PhaseSpec::new(
            phase_id.clone(),
            components,
            map_physical_state(phase.physical_state),
            map_phase_model(phase.model),
        )?;
        for component in &phase.components {
            initial.push((
                PhaseComponentId::new(phase_id.clone(), component.substance.clone()),
                component.initial_moles,
            ));
            if let Some(source_library) = &component.source_library {
                component_libraries.insert(component.substance.clone(), source_library.clone());
            }
        }
        specs.push(spec);
    }
    Ok((specs, initial, component_libraries))
}

fn apply_lookup_policy(
    spec: &mut SubstanceSystemSpec,
    lookup: ValidatedLookup,
    component_libraries: HashMap<String, String>,
) {
    match lookup {
        ValidatedLookup::Default => {
            if !component_libraries.is_empty() {
                *spec = spec.clone().with_lookup_policy(
                    Vec::new(),
                    Vec::new(),
                    Some(component_libraries),
                    false,
                );
            }
        }
        ValidatedLookup::Explicit {
            priority_libraries,
            permitted_libraries,
            explicit_search_instructions,
            search_in_nist,
        } => {
            let mut instructions = explicit_search_instructions
                .into_iter()
                .collect::<HashMap<_, _>>();
            // A candidate preview is an explicit user confirmation. Its
            // provenance therefore takes precedence over a broad policy.
            instructions.extend(component_libraries);
            let instructions = (!instructions.is_empty()).then_some(instructions);
            let phases = spec.phases().to_vec();
            *spec = SubstanceSystemSpec::from_phases(phases)
                .expect("validated phase specification remains valid")
                .with_lookup_policy(
                    priority_libraries,
                    permitted_libraries,
                    instructions,
                    search_in_nist,
                );
        }
    }
}

fn build_solve_options(
    solver: &ValidatedSolver,
    diagnostics: &crate::gui::equilibrium_gui_model::ValidatedDiagnostics,
) -> Result<EquilibriumSolveOptions, EquilibriumGuiRequestError> {
    let mut options = EquilibriumSolveOptions::new();
    match &solver.selection {
        EquilibriumSolverDraft::ProductionDefault => {
            options = options.with_production_cascade();
        }
        EquilibriumSolverDraft::SingleBackend { backend } => {
            options =
                options.with_solver_policy(SolverPolicy::Single(map_solver_backend(*backend)))?;
        }
        EquilibriumSolverDraft::CustomCascade { backends } => {
            options = options.with_solver_policy(SolverPolicy::Cascade(
                backends.iter().copied().map(map_solver_backend).collect(),
            ))?;
        }
    }
    if let Some(tolerance) = solver.tolerance {
        options = options.with_tolerance(tolerance)?;
    }
    if let Some(max_iterations) = solver.max_iterations {
        options = options.with_max_iterations(max_iterations)?;
    }
    if let Some(budget) = solver.cascade_budget {
        options = options.with_solver_budget(SolverCascadeBudget::new(
            budget.max_attempts,
            budget.max_iterations_per_attempt,
            budget.max_total_iterations,
        ))?;
    }
    if let Some(policy) = &solver.trace_seed_policy {
        options = options.with_trace_seed_policy(match policy {
            ValidatedTraceSeedPolicy::Absolute { floor } => {
                TraceSpeciesSeedPolicy::Absolute { floor: *floor }
            }
            ValidatedTraceSeedPolicy::RelativeToLargestInitialMole {
                fraction,
                minimum_floor,
            } => TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
                fraction: *fraction,
                minimum_floor: *minimum_floor,
            },
        });
    }
    options = options.with_scaling(solver.scaling_enabled);
    options = options.with_timing_mode(if diagnostics.collect_timing {
        EquilibriumTimingMode::Enabled
    } else {
        EquilibriumTimingMode::Disabled
    });
    let diagnostics_mode = match diagnostics.phase_lifecycle_trace {
        GuiPhaseLifecycleTrace::Off => EquilibriumDiagnosticsMode::Disabled,
        GuiPhaseLifecycleTrace::Summary => EquilibriumDiagnosticsMode::Summary,
        GuiPhaseLifecycleTrace::PhaseLifecycle => EquilibriumDiagnosticsMode::PhaseLifecycle,
        GuiPhaseLifecycleTrace::Detailed => EquilibriumDiagnosticsMode::Detailed,
    };
    let range_policy = match diagnostics.range_lifecycle_trace {
        GuiRangeLifecycleTrace::Endpoints => EquilibriumRangeDiagnosticsPolicy::Endpoints,
        GuiRangeLifecycleTrace::TransitionsOnly => {
            EquilibriumRangeDiagnosticsPolicy::TransitionsOnly
        }
        GuiRangeLifecycleTrace::EveryPoint => EquilibriumRangeDiagnosticsPolicy::EveryPoint,
    };
    options = options.with_diagnostics(
        EquilibriumDiagnosticsOptions::enabled(diagnostics_mode).with_range_policy(range_policy),
    );
    options = options.with_keq_validation_mode(match diagnostics.keq_validation {
        GuiKeqValidationMode::Off => EquilibriumConstantValidationMode::Off,
        GuiKeqValidationMode::WhenApplicable => EquilibriumConstantValidationMode::WhenApplicable,
        GuiKeqValidationMode::Required => EquilibriumConstantValidationMode::Required,
    });
    Ok(options)
}

fn apply_phase_mode(
    request: PhaseEquilibriumPipelineRequest,
    phase_mode: ValidatedPhaseMode,
) -> Result<PhaseEquilibriumPipelineRequest, EquilibriumGuiRequestError> {
    match phase_mode {
        ValidatedPhaseMode::FixedDeclared => Ok(request.with_fixed_declared_phases()),
        ValidatedPhaseMode::Bounded {
            phase_epsilon,
            dg_create,
            dg_keep,
            max_phase_iterations,
        } => {
            let policy =
                PhaseControlPolicy::with_explicit_hysteresis(phase_epsilon, dg_create, dg_keep)?
                    .with_max_phase_iterations(max_phase_iterations)?;
            Ok(request.with_phase_control_policy(policy))
        }
    }
}

fn map_physical_state(
    state: GuiPhysicalState,
) -> crate::Thermodynamics::physical_state::PhysicalState {
    match state {
        GuiPhysicalState::Gas => crate::Thermodynamics::physical_state::PhysicalState::Gas,
        GuiPhysicalState::Liquid => crate::Thermodynamics::physical_state::PhysicalState::Liquid,
        GuiPhysicalState::Solid => crate::Thermodynamics::physical_state::PhysicalState::Solid,
        GuiPhysicalState::Condensed => {
            crate::Thermodynamics::physical_state::PhysicalState::Condensed
        }
    }
}

fn map_phase_model(model: GuiPhaseModel) -> PhaseModel {
    match model {
        GuiPhaseModel::IdealGas => PhaseModel::IdealGas,
        GuiPhaseModel::IdealSolution => PhaseModel::IdealSolution,
        GuiPhaseModel::PureCondensed => PhaseModel::PureCondensed,
    }
}

fn map_solver_backend(backend: GuiSolverBackend) -> SolverBackend {
    match backend {
        GuiSolverBackend::RstLm => {
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt)
        }
        GuiSolverBackend::RstMinpackLm => {
            SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt)
        }
        GuiSolverBackend::RstNielsenLm => {
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt)
        }
        GuiSolverBackend::RstTrustRegionLm => {
            SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt)
        }
        GuiSolverBackend::RstPowellDogleg => {
            SolverBackend::RustedSciThe(RustedSciTheSolver::PowellDogleg)
        }
        GuiSolverBackend::RstDampedNewton => {
            SolverBackend::RustedSciThe(RustedSciTheSolver::DampedNewton)
        }
        GuiSolverBackend::LegacyLm => SolverBackend::Legacy(LegacyEquilibriumSolver::LM),
        GuiSolverBackend::LegacyNr => SolverBackend::Legacy(LegacyEquilibriumSolver::NR),
        GuiSolverBackend::LegacyTr => SolverBackend::Legacy(LegacyEquilibriumSolver::TR),
    }
}

fn linspace(start: f64, end: f64, count: usize) -> Vec<f64> {
    if count == 1 {
        return vec![start];
    }
    let denominator = (count - 1) as f64;
    (0..count)
        .map(|index| start + (end - start) * index as f64 / denominator)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gui::equilibrium_gui_model::{EquilibriumGuiDocument, EquilibriumProblemDraft};

    #[test]
    fn default_document_builds_a_point_request_without_opening_the_repository() {
        let validated = EquilibriumGuiDocument::new()
            .validate_for_run()
            .expect("default document is runnable");
        let request = build_equilibrium_request(validated, None).expect("request builds");
        assert!(!request.is_range());
    }

    #[test]
    fn descending_range_becomes_a_typed_range_request() {
        let mut document = EquilibriumGuiDocument::new();
        let EquilibriumProblemDraft::FixedPt { temperature, .. } = &mut document.config.problem
        else {
            panic!("default problem must be P,T");
        };
        *temperature = crate::gui::equilibrium_gui_model::TemperatureDraft::Range {
            start_k: "1200".into(),
            end_k: "300".into(),
            point_count: "4".into(),
        };
        let validated = document.validate_for_run().expect("range validates");
        let request = build_equilibrium_request(validated, None).expect("range request builds");
        assert!(request.is_range());
    }

    #[test]
    fn fixed_ph_crosses_the_request_boundary_as_a_typed_request() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.problem = EquilibriumProblemDraft::FixedPh {
            pressure_pa: "101325".into(),
            reference_pressure_pa: "101325".into(),
            target_enthalpy_j: "1000".into(),
            temperature_bounds:
                crate::gui::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
        };
        let validated = document.validate_for_run().expect("P,H validates");
        let request = build_equilibrium_request(validated, None).expect("P,H request builds");
        assert!(matches!(request, EquilibriumGuiSolveRequest::Ph(_)));
    }

    #[test]
    fn solver_backend_mapping_covers_rst_and_legacy_families() {
        assert!(matches!(
            map_solver_backend(GuiSolverBackend::RstLm),
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt)
        ));
        assert!(matches!(
            map_solver_backend(GuiSolverBackend::LegacyNr),
            SolverBackend::Legacy(LegacyEquilibriumSolver::NR)
        ));
    }

    #[test]
    fn lifecycle_diagnostics_map_to_the_canonical_engine_policy() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.diagnostics.phase_lifecycle_trace =
            crate::gui::equilibrium_gui_model::GuiPhaseLifecycleTrace::Detailed;
        document.config.diagnostics.range_lifecycle_trace =
            crate::gui::equilibrium_gui_model::GuiRangeLifecycleTrace::TransitionsOnly;
        let validated = document.validate_for_run().expect("document validates");
        let options = build_solve_options(&validated.solver, &validated.diagnostics)
            .expect("diagnostic policy maps");
        assert_eq!(
            options.diagnostics_options().mode(),
            EquilibriumDiagnosticsMode::Detailed
        );
        assert_eq!(
            options.diagnostics_options().range_policy(),
            EquilibriumRangeDiagnosticsPolicy::TransitionsOnly
        );
    }

    #[test]
    fn permitted_libraries_bound_candidate_preview_and_preserve_priority_order() {
        let candidate = ValidatedCandidatePolicy {
            element_mode: crate::gui::equilibrium_gui_model::GuiElementSearchMode::Exact,
            physical_states: Vec::new(),
            temperature_lower_k: 300.0,
            temperature_upper_k: 1000.0,
            max_candidates: 10,
        };
        let lookup = ValidatedLookup::Explicit {
            priority_libraries: vec!["NASA_gas".into(), "NASA_cond".into()],
            permitted_libraries: vec!["NASA_cond".into(), "NIST".into()],
            explicit_search_instructions: Default::default(),
            search_in_nist: false,
        };
        let policy = build_candidate_policy(&candidate, &lookup).expect("policy validates");
        assert_eq!(
            policy.library_preference(),
            &["NASA_cond".to_string(), "NIST".to_string()]
        );
    }
}
