//! Conversion boundary from the equilibrium GUI model to the production API.
//!
//! View code must stop at [`ValidatedEquilibriumGuiConfig`]. This module is the
//! only place that knows how GUI enums map to `PhaseSpec`, `PhaseComponentId`,
//! solver policies, and the fixed-pressure equilibrium workflow.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::EquilibriumDiagnosticEvent;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::EquilibriumExecutionControl;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::PhRangeSolution;
use crate::Thermodynamics::ChemEquilibrium::prelude::{
    CandidateSelectionError, ElementSearchMode, EquilibriumCalculator,
    EquilibriumCalculatorBuilder, EquilibriumCalculatorError, EquilibriumCandidatePolicy,
    EquilibriumCandidateSelectionReport, EquilibriumCandidateSelector,
    EquilibriumConstantValidationMode, EquilibriumDiagnosticsMode, EquilibriumDiagnosticsOptions,
    EquilibriumRangeDiagnosticsPolicy, EquilibriumSolveOptions, EquilibriumTimingMode,
    FixedPressureEnthalpySolution, InitialPhaseSet, LegacyEquilibriumSolver, PhSolveMode,
    PhaseComponentId, PhaseControlPolicy, PhaseId, PhaseModel, PhaseSpec,
    ResolvedPhaseEquilibriumOutcome, RustedSciTheSolver, SolverBackend, SolverCascadeBudget,
    SolverPolicy, SubstanceSystemFactoryError, TemperatureBounds, TemperatureRangeSolution,
    ThermoRepository, TraceSpeciesSeedPolicy,
};
use crate::gui::equilibrium_gui_model::{
    EquilibriumSolverDraft, GuiInitialPhasePolicyDraft, GuiKeqValidationMode, GuiPhSolveMode,
    GuiPhaseLifecycleTrace, GuiPhaseModel, GuiPhysicalState, GuiRangeLifecycleTrace,
    GuiSolverBackend, ValidatedCandidatePolicy, ValidatedEquilibriumGuiConfig, ValidatedInventory,
    ValidatedLookup, ValidatedPhaseMode, ValidatedProblem, ValidatedSolver, ValidatedTemperature,
    ValidatedTraceSeedPolicy,
};
use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

/// A validated request ready for the production point or range facade.
pub(crate) enum EquilibriumGuiSolveRequest {
    /// Application-facade request used by the migrated worker path.
    Facade(EquilibriumCalculatorBuilder),
}

/// Builds the new application facade from validated GUI state.
///
/// This is the only production conversion path from validated GUI state to
/// the equilibrium engine.
pub fn build_equilibrium_calculator(
    config: ValidatedEquilibriumGuiConfig,
    repository: Option<Arc<ThermoRepository>>,
) -> Result<EquilibriumCalculatorBuilder, EquilibriumGuiRequestError> {
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
    let mut builder =
        EquilibriumCalculator::from_phases(phases).initial_phase_moles(sparse_initial_moles);
    if let Some(repository) = repository {
        builder = builder.with_repository(repository);
    }
    match lookup {
        ValidatedLookup::Default => {
            if !component_libraries.is_empty() {
                builder = builder.component_library_instructions(component_libraries);
            }
        }
        ValidatedLookup::Explicit {
            priority_libraries,
            permitted_libraries,
            explicit_search_instructions,
            search_in_nist,
        } => {
            let mut instructions = explicit_search_instructions;
            instructions.extend(component_libraries);
            builder = builder
                .prefer_libraries(priority_libraries)
                .permit_libraries(permitted_libraries);
            if !instructions.is_empty() {
                builder = builder.component_library_instructions(instructions);
            }
            if search_in_nist {
                builder = builder.with_exact_state_nist_fallback();
            } else {
                builder = builder.offline_only();
            }
        }
    }
    builder = builder
        .solve_options(build_solve_options(&solver, &diagnostics)?)
        .ph_solve_mode(map_ph_solve_mode(solver.ph_solve_mode));
    if let ValidatedPhaseMode::Bounded {
        phase_epsilon,
        dg_create,
        dg_keep,
        max_phase_iterations,
        initial_phase_policy,
    } = phase_mode
    {
        let policy =
            PhaseControlPolicy::with_explicit_hysteresis(phase_epsilon, dg_create, dg_keep)?
                .with_max_phase_iterations(max_phase_iterations)?
                .with_initial_phase_set(match initial_phase_policy {
                    GuiInitialPhasePolicyDraft::FromInitialMoles => {
                        InitialPhaseSet::FromInitialMoles
                    }
                    GuiInitialPhasePolicyDraft::AllDeclaredCandidates => {
                        InitialPhaseSet::AllCandidatePhases
                    }
                })?;
        builder = builder.phase_control(policy);
    }
    match problem {
        ValidatedProblem::FixedPt {
            pressure_pa,
            reference_pressure_pa,
            temperature,
        } => {
            builder = builder
                .pressure_pa(pressure_pa)
                .reference_pressure_pa(reference_pressure_pa);
            match temperature {
                ValidatedTemperature::Point(value) => Ok(builder.at_temperature(value)),
                ValidatedTemperature::Range {
                    start_k,
                    end_k,
                    point_count,
                } => Ok(builder.over_temperature_range(linspace(start_k, end_k, point_count))?),
            }
        }
        ValidatedProblem::FixedPh {
            pressure_pa,
            reference_pressure_pa,
            target_enthalpy_j,
            temperature_bounds,
        } => Ok(builder
            .pressure_pa(pressure_pa)
            .reference_pressure_pa(reference_pressure_pa)
            .at_total_enthalpy(
                crate::Thermodynamics::ChemEquilibrium::prelude::TotalEnthalpyJoules::new(
                    target_enthalpy_j,
                )?,
                temperature_bounds.seed_k,
                TemperatureBounds::new(temperature_bounds.lower_k, temperature_bounds.upper_k)?,
            )),
    }
}

/// Wraps the migrated builder in the worker request boundary.
pub(crate) fn build_equilibrium_facade_request(
    config: ValidatedEquilibriumGuiConfig,
    repository: Option<Arc<ThermoRepository>>,
) -> Result<EquilibriumGuiSolveRequest, EquilibriumGuiRequestError> {
    Ok(EquilibriumGuiSolveRequest::Facade(
        build_equilibrium_calculator(config, repository)?,
    ))
}

impl EquilibriumGuiSolveRequest {
    /// Attaches a worker-owned execution handle without changing the request
    /// payload or its document fingerprint.
    pub fn with_execution_control(self, control: EquilibriumExecutionControl) -> Self {
        match self {
            Self::Facade(builder) => Self::Facade(builder.execution_control(control)),
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
            Self::Facade(builder) => {
                Self::Facade(builder.diagnostic_sink(move |event| sink(event)))
            }
        }
    }

    /// Executes through the canonical production orchestration boundary.
    pub fn solve(self) -> Result<EquilibriumGuiSolveOutcome, EquilibriumGuiRequestError> {
        match self {
            Self::Facade(builder) => {
                let outcome = builder.solve().map_err(EquilibriumGuiRequestError::from)?;
                EquilibriumGuiSolveOutcome::from_calculator_outcome(outcome)
                    .map_err(EquilibriumGuiRequestError::FeatureUnavailable)
            }
        }
    }
}

/// Immutable production outcome owned by the GUI worker/result layer.
#[derive(Debug, Clone)]
pub enum EquilibriumGuiSolveOutcome {
    Point(ResolvedPhaseEquilibriumOutcome),
    Range(TemperatureRangeSolution),
    Ph(FixedPressureEnthalpySolution),
    PhRange(PhRangeSolution),
}

impl EquilibriumGuiSolveOutcome {
    /// Converts a facade outcome at the GUI boundary without rebuilding or
    /// copying accepted numerical results. P,T postprocessing is presentation
    /// data and remains owned by the facade result until the GUI range snapshot
    /// consumes a dedicated presentation adapter.
    pub fn from_calculator_outcome(
        outcome: crate::Thermodynamics::ChemEquilibrium::prelude::EquilibriumCalculatorOutcome,
    ) -> Result<Self, &'static str> {
        use crate::Thermodynamics::ChemEquilibrium::prelude::EquilibriumCalculatorOutcome;
        match outcome {
            EquilibriumCalculatorOutcome::PtPoint(point) => {
                let (resolved, solution) = point.into_parts();
                Ok(Self::Point(
                    crate::Thermodynamics::ChemEquilibrium::prelude::ResolvedPhaseEquilibriumOutcome::from_parts(
                        resolved, solution,
                    ),
                ))
            }
            EquilibriumCalculatorOutcome::PtRange(range) => {
                let (_, solution, _) = range.into_parts();
                Ok(Self::Range(solution))
            }
            EquilibriumCalculatorOutcome::PhPoint(point) => {
                let (_, solution) = point.into_parts();
                Ok(Self::Ph(solution))
            }
            EquilibriumCalculatorOutcome::PhRange(range) => {
                let (_, solution) = range.into_parts();
                Ok(Self::PhRange(solution))
            }
        }
    }
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
    /// The application facade rejected a typed builder input.
    Calculator(EquilibriumCalculatorError),
    /// Candidate discovery failed before a solver request could be built.
    CandidateSelection(CandidateSelectionError),
}

impl fmt::Display for EquilibriumGuiRequestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FeatureUnavailable(feature) => write!(f, "feature is unavailable: {feature}"),
            Self::PhaseSpecification(error) => write!(f, "phase specification failed: {error}"),
            Self::EngineInput(error) => write!(f, "equilibrium request failed validation: {error}"),
            Self::Calculator(error) => write!(f, "calculator request failed validation: {error}"),
            Self::CandidateSelection(error) => write!(f, "candidate selection failed: {error}"),
        }
    }
}

impl std::error::Error for EquilibriumGuiRequestError {}

impl From<EquilibriumCalculatorError> for EquilibriumGuiRequestError {
    fn from(value: EquilibriumCalculatorError) -> Self {
        Self::Calculator(value)
    }
}

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
    let mut diagnostic_options =
        EquilibriumDiagnosticsOptions::enabled(diagnostics_mode).with_range_policy(range_policy);
    if let Some(max_events) = diagnostics.max_lifecycle_events {
        diagnostic_options = diagnostic_options.with_max_events(max_events);
    }
    options = options.with_diagnostics(diagnostic_options);
    options = options.with_keq_validation_mode(match diagnostics.keq_validation {
        GuiKeqValidationMode::Off => EquilibriumConstantValidationMode::Off,
        GuiKeqValidationMode::WhenApplicable => EquilibriumConstantValidationMode::WhenApplicable,
        GuiKeqValidationMode::Required => EquilibriumConstantValidationMode::Required,
    });
    Ok(options)
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

fn map_ph_solve_mode(mode: GuiPhSolveMode) -> PhSolveMode {
    match mode {
        GuiPhSolveMode::Auto => PhSolveMode::Auto,
        GuiPhSolveMode::Monolithic => PhSolveMode::Monolithic,
        GuiPhSolveMode::NestedTemperature => PhSolveMode::NestedTemperature,
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
    use crate::gui::equilibrium_gui_model::{
        EquilibriumGuiDocument, EquilibriumLookupDraft, EquilibriumProblemDraft,
    };

    #[test]
    fn default_document_builds_a_facade_request_without_opening_the_repository() {
        let validated = EquilibriumGuiDocument::new()
            .validate_for_run()
            .expect("default document is runnable");
        let request =
            build_equilibrium_facade_request(validated, None).expect("facade request builds");
        assert!(matches!(request, EquilibriumGuiSolveRequest::Facade(_)));
    }

    #[test]
    fn validated_document_also_builds_the_application_facade_request() {
        let validated = EquilibriumGuiDocument::new()
            .validate_for_run()
            .expect("default document is runnable");
        let _builder = build_equilibrium_calculator(validated, None)
            .expect("validated GUI state must map to the calculator facade");
    }

    #[test]
    fn facade_request_executes_and_returns_the_existing_gui_outcome_shape() {
        // A solve needs a concrete local-library policy. `Default` deliberately
        // preserves repository defaults and is not itself a promise that every
        // installation has a searchable catalog configured.
        let mut document = EquilibriumGuiDocument::new();
        document.config.lookup = EquilibriumLookupDraft::Explicit {
            priority_libraries: vec!["NASA_gas".into()],
            permitted_libraries: vec!["NASA_gas".into()],
            explicit_search_instructions: Default::default(),
            search_in_nist: false,
        };
        let validated = document
            .validate_for_run()
            .expect("default document is runnable");
        let request = build_equilibrium_facade_request(validated, None)
            .expect("validated GUI state must build a facade request");
        let outcome = request
            .solve()
            .expect("facade-backed GUI request must solve");
        assert!(matches!(outcome, EquilibriumGuiSolveOutcome::Point(_)));
    }

    #[test]
    fn descending_range_becomes_a_facade_request() {
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
        let request =
            build_equilibrium_facade_request(validated, None).expect("range facade request builds");
        assert!(matches!(request, EquilibriumGuiSolveRequest::Facade(_)));
    }

    #[test]
    fn fixed_ph_crosses_the_facade_request_boundary() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.problem = EquilibriumProblemDraft::FixedPh {
            pressure_pa: "101325".into(),
            reference_pressure_pa: "101325".into(),
            target_enthalpy_j: "1000".into(),
            temperature_bounds:
                crate::gui::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
        };
        let validated = document.validate_for_run().expect("P,H validates");
        let request =
            build_equilibrium_facade_request(validated, None).expect("P,H facade request builds");
        assert!(matches!(request, EquilibriumGuiSolveRequest::Facade(_)));
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
        document.config.diagnostics.max_lifecycle_events = "17".into();
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
        assert_eq!(options.diagnostics_options().max_events(), 17);
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
