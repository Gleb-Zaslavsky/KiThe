//! High-level application facade for the canonical equilibrium workflows.
//!
//! The lower-level `ResolvedPhase*Request` types deliberately expose every
//! physical and numerical boundary. They are appropriate for integrations and
//! diagnostics, but unnecessarily verbose for an application that begins with
//! phase declarations, an inventory, and one requested thermodynamic mode.
//! This module owns only request assembly:
//!
//! ```text
//! phase declarations + lookup policy + inventory + conditions/options
//!                              |
//!                              v
//!              resolve -> canonical P,T/P,H workflow -> typed outcome
//! ```
//!
//! It never implements a solver, fallback, or cache. The returned outcomes
//! retain the resolved system, lookup provenance, numerical reports, and all
//! phase-control evidence from the canonical workflows.

use std::collections::HashMap;
use std::sync::Arc;

use thiserror::Error;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_candidate_selection::{
    CandidateSelectionError, CandidateTemperatureRange, EquilibriumCandidatePhaseAssignment,
    EquilibriumCandidatePhasePlan, EquilibriumCandidatePolicy, EquilibriumCandidateSelectionReport,
    EquilibriumCandidateSelector,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EquilibriumConstraint, TemperatureBounds, TotalEnthalpyJoules,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
    EquilibriumDiagnosticsMode, EquilibriumDiagnosticsOptions, EquilibriumRangeDiagnosticsPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_element_inventory::ElementInventory;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::EquilibriumExecutionControl;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::{
    PhEnthalpyGrid, PhRangeError, PhRangeRequest, PhRangeSolution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
    FixedPressureEnthalpySolution, PhSolveMode, PhTemperatureSolveOptions,
    ResolvedPhaseEnthalpyRequest, ResolvedThermochemistry, solve_resolved_ph,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
    TemperaturePostprocessingPolicy, TemperaturePostprocessingResult,
    postprocess_temperature_range_solution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::{
    TemperatureGrid, TemperatureRangeRequest, TemperatureRangeSolution,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    SupportedPhaseModelPolicy, build_phase_equilibrium_problem,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    EquilibriumSolveOptions, PhaseControlPolicy, ResolvedPhaseEquilibriumRequest,
    solve_resolved_pt, solve_resolved_pt_transaction,
};
use crate::Thermodynamics::User_PhaseOrSolution::{
    PhaseSpec, ResolvedPhaseSystem, SubstanceSystemFactoryError, SubstanceSystemSpec,
};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::{NistFallbackPolicy, PhysicalState};
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use crate::Thermodynamics::thermo_lib_api::ThermoRepository;

/// Creates an application-facing equilibrium request.
pub struct EquilibriumCalculator;

impl EquilibriumCalculator {
    /// Starts an empty builder. A calculation must declare at least one phase,
    /// either molecular initial amounts or an elemental inventory,
    /// pressure/reference pressure, and one mode.
    pub fn builder() -> EquilibriumCalculatorBuilder {
        EquilibriumCalculatorBuilder::default()
    }

    /// Starts a builder from already validated phase declarations.
    ///
    /// This is the intended bridge for editors such as the GUI: validation of
    /// editable fields stays in the editor, while lookup and all numerical
    /// validation remain owned by the calculator facade.
    pub fn from_phases<I>(phases: I) -> EquilibriumCalculatorBuilder
    where
        I: IntoIterator<Item = PhaseSpec>,
    {
        EquilibriumCalculatorBuilder {
            phases: phases.into_iter().collect(),
            ..EquilibriumCalculatorBuilder::default()
        }
    }
}

/// A requested terminal calculation mode.
#[derive(Debug, Clone)]
pub enum EquilibriumCalculationMode {
    /// One fixed-pressure, fixed-temperature equilibrium point.
    PtPoint { temperature_kelvin: f64 },
    /// A fixed-pressure temperature continuation with accepted-state reuse.
    PtRange { temperatures: TemperatureGrid },
    /// One fixed-pressure, fixed-total-enthalpy equilibrium point.
    PhPoint {
        target_enthalpy: TotalEnthalpyJoules,
        initial_temperature_kelvin: f64,
        temperature_bounds: TemperatureBounds,
    },
    /// A fixed-pressure continuation over total-enthalpy targets.
    PhRange {
        target_enthalpies: PhEnthalpyGrid,
        initial_temperature_kelvin: f64,
        temperature_bounds: TemperatureBounds,
    },
}

/// Defines how the real species/phase universe is assembled at the facade
/// boundary. Once resolved, downstream workflows receive only the ordered
/// [`ResolvedPhaseSystem`] and do not branch on this policy.
#[derive(Debug, Clone)]
pub enum EquilibriumSpeciesUniversePolicy {
    /// Use the phase declarations already supplied to the builder.
    Explicit,
    /// Select records from the closed elemental inventory, then apply the
    /// reviewed phase plan without inferring a phase model from chemistry.
    FromElements {
        candidate_policy: EquilibriumCandidatePolicy,
        phase_plan: EquilibriumCandidatePhasePlan,
    },
}

/// Named application policies for observation of an equilibrium run.
///
/// Presets configure existing timing and diagnostic value objects only. They
/// do not alter equations, tolerances, backend order, phase policy, or the
/// accepted physical result. Individual builder methods can override a preset
/// afterwards.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumCalculatorPreset {
    /// No timing and no retained lifecycle events.
    Silent,
    /// Measure stage timing while retaining no lifecycle event trace.
    Timed,
    /// Retain bounded phase-control/TPD lifecycle evidence for one solve.
    LifecycleDiagnostics,
    /// Keep timing and every point's bounded summary for range
    /// characterization runs.
    RangeCharacterization,
}

/// Fluent assembly boundary for the production equilibrium calculator.
///
/// It keeps physical input and policy choices explicit while hiding repetitive
/// layout construction, thermochemistry preparation, and hand-off between the
/// resolved P,T/P,H request types.
#[derive(Clone)]
pub struct EquilibriumCalculatorBuilder {
    phases: Vec<PhaseSpec>,
    species_universe_policy: EquilibriumSpeciesUniversePolicy,
    initial_moles: Option<Vec<f64>>,
    sparse_initial_moles: Option<Vec<(PhaseComponentId, f64)>>,
    element_inventory: Option<ElementInventory>,
    element_numerical_seed: Option<MultiphaseInitialComposition>,
    library_priorities: Vec<String>,
    permitted_libraries: Vec<String>,
    component_library_instructions: Option<HashMap<String, String>>,
    nist_fallback_policy: NistFallbackPolicy,
    repository: Option<Arc<ThermoRepository>>,
    pressure_pa: Option<f64>,
    reference_pressure_pa: Option<f64>,
    mode: Option<EquilibriumCalculationMode>,
    solve_options: EquilibriumSolveOptions,
    phase_control_policy: Option<PhaseControlPolicy>,
    ph_solve_mode: PhSolveMode,
    ph_temperature_options: PhTemperatureSolveOptions,
    temperature_postprocessing: Option<TemperaturePostprocessingPolicy>,
}

impl Default for EquilibriumCalculatorBuilder {
    fn default() -> Self {
        Self {
            phases: Vec::new(),
            species_universe_policy: EquilibriumSpeciesUniversePolicy::Explicit,
            initial_moles: None,
            sparse_initial_moles: None,
            element_inventory: None,
            element_numerical_seed: None,
            library_priorities: Vec::new(),
            permitted_libraries: Vec::new(),
            component_library_instructions: None,
            nist_fallback_policy: NistFallbackPolicy::Disabled,
            repository: None,
            pressure_pa: None,
            reference_pressure_pa: None,
            mode: None,
            solve_options: EquilibriumSolveOptions::default(),
            phase_control_policy: None,
            ph_solve_mode: PhSolveMode::default(),
            ph_temperature_options: PhTemperatureSolveOptions::default(),
            temperature_postprocessing: None,
        }
    }
}

impl EquilibriumCalculatorBuilder {
    /// Applies a named observation preset. Numerical and physical policies are
    /// intentionally left untouched.
    pub fn preset(mut self, preset: EquilibriumCalculatorPreset) -> Self {
        match preset {
            EquilibriumCalculatorPreset::Silent => {
                self.solve_options = self
                    .solve_options
                    .with_timing_mode(EquilibriumTimingMode::Disabled)
                    .with_diagnostics(EquilibriumDiagnosticsOptions::disabled());
            }
            EquilibriumCalculatorPreset::Timed => {
                self.solve_options = self
                    .solve_options
                    .with_timing_mode(EquilibriumTimingMode::Enabled)
                    .with_diagnostics(EquilibriumDiagnosticsOptions::disabled());
            }
            EquilibriumCalculatorPreset::LifecycleDiagnostics => {
                self.solve_options = self
                    .solve_options
                    .with_timing_mode(EquilibriumTimingMode::Enabled)
                    .with_diagnostics(EquilibriumDiagnosticsOptions::enabled(
                        EquilibriumDiagnosticsMode::PhaseLifecycle,
                    ));
            }
            EquilibriumCalculatorPreset::RangeCharacterization => {
                self.solve_options = self
                    .solve_options
                    .with_timing_mode(EquilibriumTimingMode::Enabled)
                    .with_diagnostics(
                        EquilibriumDiagnosticsOptions::enabled(EquilibriumDiagnosticsMode::Summary)
                            .with_range_policy(EquilibriumRangeDiagnosticsPolicy::EveryPoint),
                    );
            }
        }
        self
    }

    /// Adds an anonymous ideal-gas phase, suitable for ordinary one-gas-phase
    /// equilibrium calculations.
    pub fn single_ideal_gas<I, S>(
        mut self,
        components: I,
    ) -> Result<Self, EquilibriumCalculatorError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.phases.push(PhaseSpec::ideal_gas(
            PhaseId::new(None),
            components.into_iter().map(Into::into).collect(),
        )?);
        Ok(self)
    }

    /// Adds one named ideal-gas phase to a multiphase declaration.
    pub fn ideal_gas_phase<I, S>(
        mut self,
        name: impl Into<String>,
        components: I,
    ) -> Result<Self, EquilibriumCalculatorError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.phases.push(PhaseSpec::ideal_gas(
            PhaseId::new(Some(name.into())),
            components.into_iter().map(Into::into).collect(),
        )?);
        Ok(self)
    }

    /// Adds a named multicomponent ideal condensed solution phase.
    pub fn ideal_solution_phase<I, S>(
        mut self,
        name: impl Into<String>,
        components: I,
        physical_state: PhysicalState,
    ) -> Result<Self, EquilibriumCalculatorError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.phases.push(PhaseSpec::ideal_solution(
            PhaseId::new(Some(name.into())),
            components.into_iter().map(Into::into).collect(),
            physical_state,
        )?);
        Ok(self)
    }

    /// Adds a one-component pure condensed phase without exposing the legacy
    /// container/phase-nature compatibility API.
    pub fn pure_condensed_phase(
        mut self,
        name: impl Into<String>,
        component: impl Into<String>,
        physical_state: PhysicalState,
    ) -> Result<Self, EquilibriumCalculatorError> {
        self.phases.push(PhaseSpec::pure_condensed(
            PhaseId::new(Some(name.into())),
            vec![component.into()],
            physical_state,
        )?);
        Ok(self)
    }

    /// Adds an already validated phase declaration for advanced activity-model
    /// integrations.
    pub fn phase(mut self, phase: PhaseSpec) -> Self {
        self.phases.push(phase);
        self
    }

    /// Chooses the species-universe assembly policy for this request.
    ///
    /// `FromElements` requires [`Self::element_inventory`] and uses the
    /// candidate selector plus the supplied physical phase plan. `Explicit`
    /// keeps the builder's declared phases as the complete universe.
    pub fn species_universe_policy(mut self, policy: EquilibriumSpeciesUniversePolicy) -> Self {
        self.species_universe_policy = policy;
        self
    }

    /// Supplies moles in the canonical declared component order.
    ///
    /// For a single phase this is exactly the order passed to
    /// [`Self::single_ideal_gas`]. Named multi-phase systems retain canonical
    /// phase ordering; use explicit [`PhaseSpec`] declarations when that order
    /// needs to be reviewed alongside the application model.
    pub fn initial_moles(mut self, moles: impl Into<Vec<f64>>) -> Self {
        self.initial_moles = Some(moles.into());
        self
    }

    /// Supplies a closed elemental inventory. The selected phase declarations
    /// remain the complete real species universe; the inventory is not turned
    /// into a synthetic species or an implicit catalog entry.
    pub fn element_inventory(mut self, inventory: ElementInventory) -> Self {
        self.element_inventory = Some(inventory);
        self
    }

    /// Supplies a numerical seed for an element-defined request.
    ///
    /// The seed is only a solver starting point. The physical conservation
    /// vector remains the separately supplied [`ElementInventory`], and the
    /// accepted build report therefore continues to identify the request as
    /// element-defined. This is useful for multiphase systems where an
    /// answer-independent feasible projection is too far from the desired
    /// phase basin.
    pub fn element_numerical_seed(mut self, seed: MultiphaseInitialComposition) -> Self {
        self.element_numerical_seed = Some(seed);
        self
    }

    /// Supplies phase-qualified initial amounts without relying on dense
    /// canonical layout order. Components omitted from the list are physical
    /// zeroes; duplicate and unknown identities are rejected during solve.
    /// This method is mutually exclusive with [`Self::initial_moles`].
    pub fn initial_phase_moles<I>(mut self, entries: I) -> Self
    where
        I: IntoIterator<Item = (PhaseComponentId, f64)>,
    {
        self.sparse_initial_moles = Some(entries.into_iter().collect());
        self
    }

    /// Selects preferred local libraries in lookup order.
    pub fn prefer_libraries<I, S>(mut self, libraries: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.library_priorities = libraries.into_iter().map(Into::into).collect();
        self
    }

    /// Restricts resolution to the supplied local libraries when non-empty.
    pub fn permit_libraries<I, S>(mut self, libraries: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.permitted_libraries = libraries.into_iter().map(Into::into).collect();
        self
    }

    /// Supplies component-specific library choices produced by a candidate
    /// selector or an editor. These instructions are lookup metadata only;
    /// they do not bypass phase-state validation or provenance reporting.
    pub fn component_library_instructions(
        mut self,
        instructions: impl IntoIterator<Item = (String, String)>,
    ) -> Self {
        self.component_library_instructions = Some(instructions.into_iter().collect());
        self
    }

    /// Uses only resolved local records and never initiates a NIST request.
    pub fn offline_only(mut self) -> Self {
        self.nist_fallback_policy = NistFallbackPolicy::Disabled;
        self
    }

    /// Enables exact-state NIST fallback after local lookup failure.
    pub fn with_exact_state_nist_fallback(mut self) -> Self {
        self.nist_fallback_policy = NistFallbackPolicy::ExactRequestedState;
        self
    }

    /// Uses an application-owned immutable repository for every phase lookup.
    pub fn with_repository(mut self, repository: Arc<ThermoRepository>) -> Self {
        self.repository = Some(repository);
        self
    }

    /// Sets the system pressure in Pa.
    pub fn pressure_pa(mut self, pressure_pa: f64) -> Self {
        self.pressure_pa = Some(pressure_pa);
        self
    }

    /// Sets the activity standard-state pressure in Pa.
    ///
    /// This is intentionally distinct from system pressure. The facade never
    /// silently sets `P0 = P`, which would erase ideal-gas pressure dependence.
    pub fn reference_pressure_pa(mut self, reference_pressure_pa: f64) -> Self {
        self.reference_pressure_pa = Some(reference_pressure_pa);
        self
    }

    /// Selects one fixed-temperature calculation.
    pub fn at_temperature(mut self, temperature_kelvin: f64) -> Self {
        self.mode = Some(EquilibriumCalculationMode::PtPoint { temperature_kelvin });
        self
    }

    /// Selects a fixed-pressure temperature continuation.
    pub fn over_temperature_range(
        mut self,
        temperatures: impl Into<Vec<f64>>,
    ) -> Result<Self, EquilibriumCalculatorError> {
        self.mode = Some(EquilibriumCalculationMode::PtRange {
            temperatures: TemperatureGrid::new(temperatures.into())?,
        });
        Ok(self)
    }

    /// Selects one fixed-pressure, fixed-total-enthalpy calculation.
    pub fn at_total_enthalpy(
        mut self,
        target_enthalpy: TotalEnthalpyJoules,
        initial_temperature_kelvin: f64,
        temperature_bounds: TemperatureBounds,
    ) -> Self {
        self.mode = Some(EquilibriumCalculationMode::PhPoint {
            target_enthalpy,
            initial_temperature_kelvin,
            temperature_bounds,
        });
        self
    }

    /// Selects a fixed-pressure continuation over total-enthalpy targets.
    pub fn over_enthalpy_range(
        mut self,
        target_enthalpies: impl Into<Vec<f64>>,
        initial_temperature_kelvin: f64,
        temperature_bounds: TemperatureBounds,
    ) -> Result<Self, EquilibriumCalculatorError> {
        self.mode = Some(EquilibriumCalculationMode::PhRange {
            target_enthalpies: PhEnthalpyGrid::new(target_enthalpies.into())?,
            initial_temperature_kelvin,
            temperature_bounds,
        });
        Ok(self)
    }

    /// Uses the default RST-first production cascade with retained legacy
    /// numerical fallbacks.
    pub fn production_cascade(mut self) -> Self {
        self.solve_options = self.solve_options.with_production_cascade();
        self
    }

    /// Selects an explicit solver policy for every inner nonlinear solve.
    pub fn solver_policy(
        mut self,
        policy: SolverPolicy,
    ) -> Result<Self, EquilibriumCalculatorError> {
        self.solve_options = self.solve_options.with_solver_policy(policy)?;
        Ok(self)
    }

    /// Replaces the shared canonical solve options.
    pub fn solve_options(mut self, options: EquilibriumSolveOptions) -> Self {
        self.solve_options = options;
        self
    }

    /// Enables or disables measured workflow timing.
    pub fn timing(mut self, mode: EquilibriumTimingMode) -> Self {
        self.solve_options = self.solve_options.with_timing_mode(mode);
        self
    }

    /// Attaches cooperative cancellation and progress reporting to the
    /// eventual worker transaction.
    pub fn execution_control(mut self, control: EquilibriumExecutionControl) -> Self {
        self.solve_options = self.solve_options.with_execution_control(control);
        self
    }

    /// Attaches a live diagnostic sink while retaining the immutable bounded
    /// report configured by the selected preset/options.
    pub fn diagnostic_sink<F>(mut self, sink: F) -> Self
    where
        F: Fn(crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::EquilibriumDiagnosticEvent)
            + Send
            + Sync
            + 'static,
    {
        let diagnostics = self
            .solve_options
            .diagnostics_options()
            .clone()
            .with_sink(sink);
        self.solve_options = self.solve_options.with_diagnostics(diagnostics);
        self
    }

    /// Enables bounded typed lifecycle diagnostics without changing numerical
    /// settings or accepted-state publication.
    pub fn diagnostics(mut self, options: EquilibriumDiagnosticsOptions) -> Self {
        self.solve_options = self.solve_options.with_diagnostics(options);
        self
    }

    /// Enables bounded phase creation/deactivation for all relevant point or
    /// continuation solves.
    pub fn phase_control(mut self, policy: PhaseControlPolicy) -> Self {
        self.phase_control_policy = Some(policy);
        self
    }

    /// Selects the P,H route. `Auto`, `Monolithic`, and `NestedTemperature`
    /// preserve their existing canonical semantics.
    pub fn ph_solve_mode(mut self, mode: PhSolveMode) -> Self {
        self.ph_solve_mode = mode;
        self
    }

    /// Replaces outer P,H scalar/coupled-solver controls.
    pub fn ph_temperature_options(mut self, options: PhTemperatureSolveOptions) -> Self {
        self.ph_temperature_options = options;
        self
    }

    /// Requests presentation-only PCHIP postprocessing after a successful
    /// P,T temperature range. Resampling across an accepted phase transition
    /// remains rejected by the existing postprocessing contract.
    pub fn temperature_postprocessing(mut self, policy: TemperaturePostprocessingPolicy) -> Self {
        self.temperature_postprocessing = Some(policy);
        self
    }

    /// Resolves phase data once and executes the selected canonical workflow.
    pub fn solve(self) -> Result<EquilibriumCalculatorOutcome, EquilibriumCalculatorError> {
        let mode = self
            .mode
            .clone()
            .ok_or(EquilibriumCalculatorError::MissingMode)?;
        if self.temperature_postprocessing.is_some()
            && !matches!(mode, EquilibriumCalculationMode::PtRange { .. })
        {
            return Err(EquilibriumCalculatorError::PostprocessingRequiresPtRange);
        }
        let pressure_pa = self
            .pressure_pa
            .ok_or(EquilibriumCalculatorError::MissingPressure)?;
        let reference_pressure_pa = self
            .reference_pressure_pa
            .ok_or(EquilibriumCalculatorError::MissingReferencePressure)?;
        if self.element_inventory.is_some()
            && (self.initial_moles.is_some() || self.sparse_initial_moles.is_some())
        {
            return Err(EquilibriumCalculatorError::ConflictingPhysicalInputs);
        }
        if self.element_numerical_seed.is_some() && self.element_inventory.is_none() {
            return Err(EquilibriumCalculatorError::MissingElementInventoryForNumericalSeed);
        }
        if let Some(inventory) = self.element_inventory.clone() {
            let (resolved, candidate_selection) =
                self.resolve_system_with_selection(Some(&mode))?;
            match mode {
                EquilibriumCalculationMode::PtPoint { temperature_kelvin } => {
                    let conditions = EquilibriumConditions::new(
                        temperature_kelvin,
                        pressure_pa,
                        reference_pressure_pa,
                    )?;
                    let build_request = crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::
                        PhaseEquilibriumBuildRequest::from_element_inventory(
                            &resolved,
                            conditions,
                            inventory.clone(),
                            self.solve_options.trace_seed_policy(),
                            SupportedPhaseModelPolicy::default(),
                        )?;
                    let build_request = match self.element_numerical_seed.clone() {
                        Some(seed) => build_request.with_numerical_seed(seed)?,
                        None => build_request,
                    };
                    let build_request = match candidate_selection.clone() {
                        Some(selection) => build_request.with_candidate_selection(selection),
                        None => build_request,
                    };
                    let bundle = build_phase_equilibrium_problem(build_request)?;
                    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
                    let numerical_seed = MultiphaseInitialComposition::from_dense(
                        &layout,
                        bundle.problem().initial_moles().to_vec(),
                    )?;
                    let request =
                        ResolvedPhaseEquilibriumRequest::new(&resolved, conditions, numerical_seed)
                            .with_element_inventory(inventory)
                            .with_candidate_selection(candidate_selection)
                            .with_solve_options(self.solve_options);
                    let request = match self.phase_control_policy {
                        Some(policy) => request.with_phase_control_policy(policy),
                        None => request,
                    };
                    let solution = solve_resolved_pt_transaction(request)?;
                    return Ok(EquilibriumCalculatorOutcome::PtPoint(
                        EquilibriumCalculatorPoint { resolved, solution },
                    ));
                }
                EquilibriumCalculationMode::PhPoint {
                    target_enthalpy,
                    initial_temperature_kelvin,
                    temperature_bounds,
                } => {
                    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
                    let request = ResolvedPhaseEnthalpyRequest::from_element_inventory(
                        &resolved,
                        inventory,
                        EquilibriumConstraint::ph_joules(
                            pressure_pa,
                            reference_pressure_pa,
                            target_enthalpy,
                            initial_temperature_kelvin,
                        )?,
                        temperature_bounds,
                        thermochemistry,
                    )?;
                    let request = match self.element_numerical_seed.clone() {
                        Some(seed) => request.with_initial_composition(seed)?,
                        None => request,
                    }
                    .with_solve_options(self.solve_options)
                    .with_ph_solve_mode(self.ph_solve_mode)
                    .with_temperature_options(self.ph_temperature_options)?;
                    let request = match candidate_selection.clone() {
                        Some(selection) => request.with_candidate_selection(selection),
                        None => request,
                    };
                    let request = match self.phase_control_policy {
                        Some(policy) => request.with_phase_control_policy(policy),
                        None => request,
                    };
                    let solution = solve_resolved_ph(request)?;
                    return Ok(EquilibriumCalculatorOutcome::PhPoint(
                        EquilibriumCalculatorPhPoint { resolved, solution },
                    ));
                }
                EquilibriumCalculationMode::PhRange {
                    target_enthalpies,
                    initial_temperature_kelvin,
                    temperature_bounds,
                } => {
                    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
                    let request = PhRangeRequest::from_element_inventory(
                        &resolved,
                        inventory,
                        pressure_pa,
                        reference_pressure_pa,
                        target_enthalpies,
                        temperature_bounds,
                        initial_temperature_kelvin,
                        thermochemistry,
                    )?;
                    let request = match self.element_numerical_seed.clone() {
                        Some(seed) => request.with_initial_composition(seed)?,
                        None => request,
                    }
                    .with_solve_options(self.solve_options)
                    .with_ph_solve_mode(self.ph_solve_mode)
                    .with_temperature_options(self.ph_temperature_options)?;
                    let request = match candidate_selection.clone() {
                        Some(selection) => request.with_candidate_selection(selection),
                        None => request,
                    };
                    let request = match self.phase_control_policy {
                        Some(policy) => request.with_phase_control_policy(policy),
                        None => request,
                    };
                    let solution = request.solve()?;
                    return Ok(EquilibriumCalculatorOutcome::PhRange(
                        EquilibriumCalculatorPhRange { resolved, solution },
                    ));
                }
                EquilibriumCalculationMode::PtRange { temperatures } => {
                    let request = TemperatureRangeRequest::from_element_inventory(
                        &resolved,
                        inventory,
                        pressure_pa,
                        reference_pressure_pa,
                        temperatures,
                    )?;
                    let request = match self.element_numerical_seed.clone() {
                        Some(seed) => request.with_initial_composition(seed)?,
                        None => request,
                    }
                    .with_solve_options(self.solve_options);
                    let request = match candidate_selection.clone() {
                        Some(selection) => request.with_candidate_selection(selection),
                        None => request,
                    };
                    let request = match self.phase_control_policy {
                        Some(policy) => request.with_phase_control_policy(policy),
                        None => request,
                    };
                    let solution = request.solve()?;
                    let postprocessing = self
                        .temperature_postprocessing
                        .as_ref()
                        .map(|policy| postprocess_temperature_range_solution(&solution, policy))
                        .transpose()?;
                    return Ok(EquilibriumCalculatorOutcome::PtRange(
                        EquilibriumCalculatorTemperatureRange {
                            resolved,
                            solution,
                            postprocessing,
                        },
                    ));
                }
            }
        }
        let (resolved, initial) = self.resolve_and_prepare_initial()?;

        match mode {
            EquilibriumCalculationMode::PtPoint { temperature_kelvin } => {
                let conditions = EquilibriumConditions::new(
                    temperature_kelvin,
                    pressure_pa,
                    reference_pressure_pa,
                )?;
                let request = self.pt_request(&resolved, conditions, initial);
                let solution = solve_resolved_pt(request)?;
                Ok(EquilibriumCalculatorOutcome::PtPoint(
                    EquilibriumCalculatorPoint { resolved, solution },
                ))
            }
            EquilibriumCalculationMode::PtRange { temperatures } => {
                let request = TemperatureRangeRequest::new(
                    &resolved,
                    initial,
                    pressure_pa,
                    reference_pressure_pa,
                    temperatures,
                )?
                .with_solve_options(self.solve_options);
                let request = match self.phase_control_policy {
                    Some(policy) => request.with_phase_control_policy(policy),
                    None => request,
                };
                let solution = request.solve()?;
                let postprocessing = self
                    .temperature_postprocessing
                    .as_ref()
                    .map(|policy| postprocess_temperature_range_solution(&solution, policy))
                    .transpose()?;
                Ok(EquilibriumCalculatorOutcome::PtRange(
                    EquilibriumCalculatorTemperatureRange {
                        resolved,
                        solution,
                        postprocessing,
                    },
                ))
            }
            EquilibriumCalculationMode::PhPoint {
                target_enthalpy,
                initial_temperature_kelvin,
                temperature_bounds,
            } => {
                let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
                let request = ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
                    &resolved,
                    initial,
                    EquilibriumConstraint::ph_joules(
                        pressure_pa,
                        reference_pressure_pa,
                        target_enthalpy,
                        initial_temperature_kelvin,
                    )?,
                    temperature_bounds,
                    thermochemistry,
                )?
                .with_solve_options(self.solve_options)
                .with_ph_solve_mode(self.ph_solve_mode)
                .with_temperature_options(self.ph_temperature_options)?;
                let request = match self.phase_control_policy {
                    Some(policy) => request.with_phase_control_policy(policy),
                    None => request,
                };
                let solution = solve_resolved_ph(request)?;
                Ok(EquilibriumCalculatorOutcome::PhPoint(
                    EquilibriumCalculatorPhPoint { resolved, solution },
                ))
            }
            EquilibriumCalculationMode::PhRange {
                target_enthalpies,
                initial_temperature_kelvin,
                temperature_bounds,
            } => {
                let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
                let request = PhRangeRequest::from_resolved_thermochemistry(
                    &resolved,
                    initial,
                    pressure_pa,
                    reference_pressure_pa,
                    target_enthalpies,
                    temperature_bounds,
                    initial_temperature_kelvin,
                    thermochemistry,
                )?
                .with_solve_options(self.solve_options)
                .with_ph_solve_mode(self.ph_solve_mode)
                .with_temperature_options(self.ph_temperature_options)?;
                let request = match self.phase_control_policy {
                    Some(policy) => request.with_phase_control_policy(policy),
                    None => request,
                };
                let solution = request.solve()?;
                Ok(EquilibriumCalculatorOutcome::PhRange(
                    EquilibriumCalculatorPhRange { resolved, solution },
                ))
            }
        }
    }

    fn resolve_and_prepare_initial(
        &self,
    ) -> Result<(ResolvedPhaseSystem, MultiphaseInitialComposition), EquilibriumCalculatorError>
    {
        let resolved = self.resolve_system()?;
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
        if self.initial_moles.is_none() && self.sparse_initial_moles.is_none() {
            return Err(EquilibriumCalculatorError::MissingInitialMoles);
        }
        let dense_initial_moles = self.initial_moles.clone();
        let sparse_initial_moles = self.sparse_initial_moles.clone();
        let initial = match (dense_initial_moles, sparse_initial_moles) {
            (Some(moles), None) => MultiphaseInitialComposition::from_dense(&layout, moles)?,
            (None, Some(entries)) => MultiphaseInitialComposition::from_sparse(&layout, entries)?,
            _ => unreachable!("initial composition conflict is checked above"),
        };
        Ok((resolved, initial))
    }

    fn resolve_system(&self) -> Result<ResolvedPhaseSystem, EquilibriumCalculatorError> {
        self.resolve_system_with_selection(None)
            .map(|(resolved, _)| resolved)
    }

    fn resolve_system_with_selection(
        &self,
        mode: Option<&EquilibriumCalculationMode>,
    ) -> Result<
        (
            ResolvedPhaseSystem,
            Option<EquilibriumCandidateSelectionReport>,
        ),
        EquilibriumCalculatorError,
    > {
        if let EquilibriumSpeciesUniversePolicy::FromElements {
            candidate_policy,
            phase_plan,
        } = &self.species_universe_policy
        {
            if !self.phases.is_empty() {
                return Err(EquilibriumCalculatorError::ConflictingSpeciesUniverseDeclarations);
            }
            let inventory = self
                .element_inventory
                .as_ref()
                .ok_or(EquilibriumCalculatorError::MissingElementInventoryForUniverseSelection)?;
            let repository = match self.repository.clone() {
                Some(repository) => repository,
                None => ThermoData::try_default_repository()
                    .map_err(SubstanceSystemFactoryError::Repository)?,
            };
            let candidate_policy = Self::candidate_policy_for_mode(
                candidate_policy,
                mode.ok_or(EquilibriumCalculatorError::MissingMode)?,
            )?;
            let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
                .select_inventory(inventory, candidate_policy)?;
            let spec = phase_plan.build_spec(&selection)?;
            let resolved = spec.resolve_with_repository(repository)?;
            return Ok((resolved, Some(selection)));
        }
        if self.phases.is_empty() {
            return Err(EquilibriumCalculatorError::MissingPhases);
        }
        if self.initial_moles.is_some() && self.sparse_initial_moles.is_some() {
            return Err(EquilibriumCalculatorError::ConflictingInitialComposition);
        }
        let mut spec = SubstanceSystemSpec::from_phases(self.phases.clone())?;
        // An empty builder lookup section means "use the repository default";
        // it must not replace the default catalog policy with empty lists.
        if !self.library_priorities.is_empty()
            || !self.permitted_libraries.is_empty()
            || self.component_library_instructions.is_some()
        {
            spec = spec.with_lookup_policy(
                self.library_priorities.clone(),
                self.permitted_libraries.clone(),
                self.component_library_instructions.clone(),
                false,
            );
        }
        let spec = spec.with_nist_fallback_policy(self.nist_fallback_policy);
        let resolved = match self.repository.clone() {
            Some(repository) => spec.resolve_with_repository(repository)?,
            None => spec.resolve()?,
        };
        Ok((resolved, None))
    }

    fn candidate_policy_for_mode(
        policy: &EquilibriumCandidatePolicy,
        mode: &EquilibriumCalculationMode,
    ) -> Result<EquilibriumCandidatePolicy, CandidateSelectionError> {
        let (required_lower, required_upper) = match mode {
            EquilibriumCalculationMode::PtPoint { temperature_kelvin } => {
                (*temperature_kelvin, *temperature_kelvin)
            }
            EquilibriumCalculationMode::PtRange { temperatures } => (
                temperatures
                    .values()
                    .iter()
                    .copied()
                    .fold(f64::INFINITY, f64::min),
                temperatures
                    .values()
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max),
            ),
            EquilibriumCalculationMode::PhPoint {
                temperature_bounds, ..
            }
            | EquilibriumCalculationMode::PhRange {
                temperature_bounds, ..
            } => (temperature_bounds.lower(), temperature_bounds.upper()),
        };
        let required = CandidateTemperatureRange::new(required_lower, required_upper)?;
        match policy.temperature_range() {
            Some(configured)
                if configured.lower() <= required.lower()
                    && configured.upper() >= required.upper() =>
            {
                Ok(policy.clone())
            }
            Some(configured) => Err(CandidateSelectionError::InvalidPolicy(format!(
                "candidate temperature range [{:.6}, {:.6}] does not cover calculation range [{:.6}, {:.6}]",
                configured.lower(),
                configured.upper(),
                required.lower(),
                required.upper(),
            ))),
            None => policy
                .clone()
                .with_temperature_range(required.lower(), required.upper()),
        }
    }

    fn pt_request<'a>(
        &self,
        resolved: &'a ResolvedPhaseSystem,
        conditions: EquilibriumConditions,
        initial: MultiphaseInitialComposition,
    ) -> ResolvedPhaseEquilibriumRequest<'a> {
        let request = ResolvedPhaseEquilibriumRequest::new(resolved, conditions, initial)
            .with_solve_options(self.solve_options.clone());
        match self.phase_control_policy.clone() {
            Some(policy) => request.with_phase_control_policy(policy),
            None => request,
        }
    }
}

/// One accepted P,T point paired with the immutable resolved lookup report.
#[derive(Debug, Clone)]
pub struct EquilibriumCalculatorPoint {
    resolved: ResolvedPhaseSystem,
    solution: MultiphaseEquilibriumSolution,
}

impl EquilibriumCalculatorPoint {
    pub fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    pub fn solution(&self) -> &MultiphaseEquilibriumSolution {
        &self.solution
    }

    /// Moves the resolved system and accepted solution to an integration
    /// adapter without cloning either payload.
    pub fn into_parts(self) -> (ResolvedPhaseSystem, MultiphaseEquilibriumSolution) {
        (self.resolved, self.solution)
    }
}

/// One accepted P,H point paired with the immutable resolved lookup report.
#[derive(Debug, Clone)]
pub struct EquilibriumCalculatorPhPoint {
    resolved: ResolvedPhaseSystem,
    solution: FixedPressureEnthalpySolution,
}

impl EquilibriumCalculatorPhPoint {
    pub fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    pub fn solution(&self) -> &FixedPressureEnthalpySolution {
        &self.solution
    }

    /// Moves the resolved system and accepted P,H solution to an adapter.
    pub fn into_parts(self) -> (ResolvedPhaseSystem, FixedPressureEnthalpySolution) {
        (self.resolved, self.solution)
    }
}

/// One accepted P,T range plus optional display-only postprocessing.
#[derive(Debug, Clone)]
pub struct EquilibriumCalculatorTemperatureRange {
    resolved: ResolvedPhaseSystem,
    solution: TemperatureRangeSolution,
    postprocessing: Option<TemperaturePostprocessingResult>,
}

impl EquilibriumCalculatorTemperatureRange {
    pub fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    pub fn solution(&self) -> &TemperatureRangeSolution {
        &self.solution
    }

    pub fn postprocessing(&self) -> Option<&TemperaturePostprocessingResult> {
        self.postprocessing.as_ref()
    }

    /// Moves a range result to a presentation adapter without rebuilding raw
    /// points or the optional postprocessed series.
    pub fn into_parts(
        self,
    ) -> (
        ResolvedPhaseSystem,
        TemperatureRangeSolution,
        Option<TemperaturePostprocessingResult>,
    ) {
        (self.resolved, self.solution, self.postprocessing)
    }
}

/// One accepted P,H target range paired with its resolved lookup report.
#[derive(Debug, Clone)]
pub struct EquilibriumCalculatorPhRange {
    resolved: ResolvedPhaseSystem,
    solution: PhRangeSolution,
}

impl EquilibriumCalculatorPhRange {
    pub fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    pub fn solution(&self) -> &PhRangeSolution {
        &self.solution
    }

    /// Moves a P,H range result to a presentation adapter without copying it.
    pub fn into_parts(self) -> (ResolvedPhaseSystem, PhRangeSolution) {
        (self.resolved, self.solution)
    }
}

/// Terminal typed outcome from [`EquilibriumCalculatorBuilder::solve`].
#[derive(Debug, Clone)]
pub enum EquilibriumCalculatorOutcome {
    PtPoint(EquilibriumCalculatorPoint),
    PtRange(EquilibriumCalculatorTemperatureRange),
    PhPoint(EquilibriumCalculatorPhPoint),
    PhRange(EquilibriumCalculatorPhRange),
}

/// Typed facade assembly and execution error.
#[derive(Debug, Error)]
pub enum EquilibriumCalculatorError {
    #[error("equilibrium calculator requires at least one phase")]
    MissingPhases,
    #[error("FromElements species-universe policy requires an elemental inventory")]
    MissingElementInventoryForUniverseSelection,
    #[error("an elemental numerical seed requires an elemental inventory")]
    MissingElementInventoryForNumericalSeed,
    #[error("FromElements species-universe policy cannot be combined with explicit phases")]
    ConflictingSpeciesUniverseDeclarations,
    #[error("equilibrium calculator requires initial moles")]
    MissingInitialMoles,
    #[error(
        "equilibrium calculator accepts either dense or phase-qualified initial composition, not both"
    )]
    ConflictingInitialComposition,
    #[error(
        "equilibrium calculator accepts either an elemental inventory or molecular initial amounts, not both"
    )]
    ConflictingPhysicalInputs,
    #[error("equilibrium calculator requires a P,T or P,H mode")]
    MissingMode,
    #[error("equilibrium calculator requires system pressure in Pa")]
    MissingPressure,
    #[error("equilibrium calculator requires standard-state reference pressure in Pa")]
    MissingReferencePressure,
    #[error("temperature postprocessing is available only for a P,T temperature range")]
    PostprocessingRequiresPtRange,
    #[error(transparent)]
    Resolve(#[from] SubstanceSystemFactoryError),
    #[error(transparent)]
    CandidateSelection(#[from] CandidateSelectionError),
    #[error(transparent)]
    Preparation(ReactionExtentError),
    #[error(transparent)]
    Solve(ReactionExtentError),
    #[error(transparent)]
    PhRange(PhRangeError),
}

impl From<ReactionExtentError> for EquilibriumCalculatorError {
    fn from(value: ReactionExtentError) -> Self {
        if matches!(&value, ReactionExtentError::Preparation(_)) {
            Self::Preparation(value)
        } else {
            Self::Solve(value)
        }
    }
}

impl From<PhRangeError> for EquilibriumCalculatorError {
    fn from(value: PhRangeError) -> Self {
        match value {
            PhRangeError::InvalidProblem(error) => Self::from(error),
            PhRangeError::Point(error) => Self::PhRange(PhRangeError::Point(error)),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_candidate_selection::{
    CandidateSelectionError, CandidateTemperatureRange, EquilibriumCandidatePhaseAssignment,
    EquilibriumCandidatePhasePlan, EquilibriumCandidatePolicy, 
    EquilibriumCandidateSelector,
};
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::TemperatureResamplingGrid;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::PhaseEquilibriumInputKind;

    const PRESSURE_PA: f64 = 101_325.0;
    const REFERENCE_PRESSURE_PA: f64 = 100_000.0;

    fn gas_point_builder() -> EquilibriumCalculatorBuilder {
        EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .expect("small gas phase declaration must validate")
            .initial_moles(vec![0.1, 0.05, 1.9])
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .production_cascade()
    }

    #[test]
    fn facade_solves_pt_then_ph_without_manual_resolved_request_assembly() {
        let point = match gas_point_builder()
            .at_temperature(2_500.0)
            .solve()
            .expect("facade P,T point must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected P,T point, got {other:?}"),
        };
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(point.resolved())
            .expect("resolved gas phase must provide P,H thermochemistry");
        let target_enthalpy = thermochemistry
            .enthalpy_model()
            .evaluate_total(point.solution().component_moles(), 2_500.0)
            .expect("accepted P,T state must evaluate total enthalpy");

        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let elemental_point = match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(inventory.clone())
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_total_enthalpy(
                TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
                2_300.0,
                TemperatureBounds::new(2_100.0, 2_900.0).unwrap(),
            )
            .ph_solve_mode(PhSolveMode::NestedTemperature)
            .preset(EquilibriumCalculatorPreset::Timed)
            .solve()
            .expect("element-defined P,H point must solve")
        {
            EquilibriumCalculatorOutcome::PhPoint(point) => point,
            other => panic!("expected element-defined P,H point, got {other:?}"),
        };
        assert!((elemental_point.solution().temperature() - 2_500.0).abs() < 1.0e-4);
        assert!(
            elemental_point
                .solution()
                .equilibrium()
                .accepted_solution()
                .validation()
                .max_abs_element_balance_error
                < 1.0e-7
        );
        assert!(elemental_point.solution().report().timing().enabled());
        assert!(elemental_point.solution().report().inner_timing().enabled());

        let elemental_range = match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(inventory)
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .over_enthalpy_range(
                vec![target_enthalpy, target_enthalpy + 1_000.0],
                2_500.0,
                TemperatureBounds::new(2_100.0, 2_900.0).unwrap(),
            )
            .unwrap()
            .ph_solve_mode(PhSolveMode::NestedTemperature)
            .preset(EquilibriumCalculatorPreset::Timed)
            .solve()
            .expect("element-defined P,H range must solve")
        {
            EquilibriumCalculatorOutcome::PhRange(range) => range,
            other => panic!("expected element-defined P,H range, got {other:?}"),
        };
        assert_eq!(elemental_range.solution().points().len(), 2);
        assert_eq!(elemental_range.solution().report().continuation_points(), 1);
        assert_eq!(elemental_range.solution().report().formulation_builds(), 1);
        assert!(elemental_range.solution().report().formulation_reuses() > 0);
        assert!(
            elemental_range.solution().report().point_timing().total() > std::time::Duration::ZERO
        );
        for (point, expected_target) in elemental_range
            .solution()
            .points()
            .iter()
            .zip([target_enthalpy, target_enthalpy + 1_000.0])
        {
            let solution = point.solution();
            assert_eq!(
                solution.equilibrium().build_report().input_kind(),
                PhaseEquilibriumInputKind::ElementInventory
            );
            assert_eq!(
                solution.equilibrium().build_report().element_totals(),
                [4.0, 2.0]
            );
            assert_eq!(solution.target_enthalpy(), expected_target);
            assert!(solution.enthalpy_error().abs() <= solution.enthalpy_error_limit_joules());
        }

        let recovered = match gas_point_builder()
            .at_total_enthalpy(
                TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
                2_300.0,
                TemperatureBounds::new(2_100.0, 2_900.0).unwrap(),
            )
            .ph_solve_mode(PhSolveMode::NestedTemperature)
            .solve()
            .expect("facade P,H point must solve")
        {
            EquilibriumCalculatorOutcome::PhPoint(point) => point,
            other => panic!("expected P,H point, got {other:?}"),
        };

        assert!(
            (recovered.solution().temperature() - 2_500.0).abs() < 1.0e-4,
            "P,H facade must recover the P,T reference temperature"
        );
        assert!(
            recovered
                .solution()
                .equilibrium()
                .accepted_solution()
                .validation()
                .max_abs_element_balance_error
                < 1.0e-7
        );
        assert!(!recovered.resolved().report().nist_fallback_enabled());
    }

    #[test]
    fn elemental_ph_range_preserves_transactional_rollback_after_continuation_failure() {
        let reference = match gas_point_builder()
            .at_temperature(2_500.0)
            .solve()
            .expect("reference P,T point must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected reference P,T point, got {other:?}"),
        };
        let thermochemistry =
            ResolvedThermochemistry::from_resolved_system(reference.resolved()).unwrap();
        let target_enthalpy = thermochemistry
            .enthalpy_model()
            .evaluate_total(reference.solution().component_moles(), 2_500.0)
            .unwrap();
        let unreachable_target = target_enthalpy + 1.0e9;
        let error = EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)]).unwrap())
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .over_enthalpy_range(
                vec![target_enthalpy, unreachable_target],
                2_500.0,
                TemperatureBounds::new(2_100.0, 2_900.0).unwrap(),
            )
            .unwrap()
            .ph_solve_mode(PhSolveMode::NestedTemperature)
            .solve()
            .expect_err("an unreachable second target must abort the complete range");

        match error {
            EquilibriumCalculatorError::PhRange(PhRangeError::Point(point)) => {
                assert_eq!(point.index(), 1);
                assert_eq!(point.target_enthalpy_joules(), unreachable_target);
            }
            other => panic!("expected point-indexed P,H range rollback error, got {other:?}"),
        }
    }

    #[test]
    fn facade_accepts_element_inventory_without_a_molecular_seed() {
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let result = EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(inventory)
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(2_500.0)
            .solve()
            .expect("element-defined P,T point must solve");

        let point = match result {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected element-defined P,T point, got {other:?}"),
        };
        assert_eq!(point.solution().component_moles().len(), 3);
        assert!(
            point
                .solution()
                .component_moles()
                .iter()
                .all(|moles| moles.is_finite() && *moles >= 0.0)
        );
    }

    #[test]
    fn from_elements_policy_resolves_and_publishes_selection_provenance() {
        let repository = ThermoData::try_default_repository()
            .expect("bundled repository must load for automatic universe selection");
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let candidate_policy = EquilibriumCandidatePolicy::default()
            .with_library_preference(vec!["NASA_gas".into()])
            .with_temperature_range(300.0, 1_200.0)
            .unwrap()
            .with_max_candidates(3)
            .unwrap();
        let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
            .select_inventory(&inventory, candidate_policy.clone())
            .expect("automatic candidate selection must succeed");
        let phase_plan = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::ideal_gas(
                PhaseId::new(Some("gas".into())),
                selection
                    .selected()
                    .iter()
                    .map(|candidate| candidate.record_key().to_string())
                    .collect(),
            ),
        ]);

        let result = EquilibriumCalculator::builder()
            .species_universe_policy(EquilibriumSpeciesUniversePolicy::FromElements {
                candidate_policy,
                phase_plan,
            })
            .element_inventory(inventory)
            .with_repository(repository)
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(1_200.0)
            .solve()
            .expect("automatic element-defined P,T point must solve");
        let EquilibriumCalculatorOutcome::PtPoint(point) = result else {
            panic!("expected an automatic element-defined P,T point");
        };
        let report = point
            .solution()
            .build_report()
            .candidate_selection()
            .expect("selection provenance must reach the accepted build report");
        assert_eq!(report.selected().len(), 3);
        assert_eq!(point.resolved().phase_specs().len(), 1);
    }

    #[test]
    fn from_elements_facade_derives_point_temperature_screening_when_policy_omits_range() {
        let repository = ThermoData::try_default_repository()
            .expect("bundled repository must load for automatic universe selection");
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let candidate_policy = EquilibriumCandidatePolicy::default()
            .with_library_preference(vec!["NASA_gas".into()])
            .with_max_candidates(3)
            .unwrap();
        let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
            .select_inventory(&inventory, candidate_policy.clone())
            .expect("automatic candidate selection must succeed");
        let phase_plan = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::ideal_gas(
                PhaseId::new(Some("gas".into())),
                selection
                    .selected()
                    .iter()
                    .map(|candidate| candidate.record_key().to_string())
                    .collect(),
            ),
        ]);

        let result = EquilibriumCalculator::builder()
            .species_universe_policy(EquilibriumSpeciesUniversePolicy::FromElements {
                candidate_policy,
                phase_plan,
            })
            .element_inventory(inventory)
            .with_repository(repository)
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(1_200.0)
            .solve()
            .expect("facade must derive the point screening range");
        assert!(matches!(result, EquilibriumCalculatorOutcome::PtPoint(_)));
    }

    #[test]
    fn automatic_candidate_temperature_policy_covers_ranges_and_rejects_narrow_ph_bounds() {
        let policy = EquilibriumCandidatePolicy::default();
        let pt_range = EquilibriumCalculationMode::PtRange {
            temperatures: TemperatureGrid::new(vec![900.0, 1_000.0, 1_200.0]).unwrap(),
        };
        let derived = EquilibriumCalculatorBuilder::candidate_policy_for_mode(&policy, &pt_range)
            .expect("missing policy range must be derived from the PT grid");
        assert_eq!(
            derived.temperature_range().unwrap(),
            CandidateTemperatureRange::new(900.0, 1_200.0).unwrap()
        );

        let configured = policy.with_temperature_range(300.0, 1_000.0).unwrap();
        let ph_point = EquilibriumCalculationMode::PhPoint {
            target_enthalpy: TotalEnthalpyJoules::new(0.0).unwrap(),
            initial_temperature_kelvin: 950.0,
            temperature_bounds: TemperatureBounds::new(900.0, 1_200.0).unwrap(),
        };
        let error = EquilibriumCalculatorBuilder::candidate_policy_for_mode(&configured, &ph_point)
            .expect_err("a policy narrower than PH bounds must be rejected");
        assert!(matches!(
            error,
            CandidateSelectionError::InvalidPolicy(message)
                if message.contains("does not cover calculation range")
        ));
    }

    #[test]
    fn molecular_and_elemental_facade_inputs_match_at_pressure_distinct_from_reference_pressure() {
        assert_ne!(PRESSURE_PA, REFERENCE_PRESSURE_PA);
        let molecular_pt = match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .initial_moles(vec![2.0, 1.0, 0.0])
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(2_500.0)
            .production_cascade()
            .solve()
            .expect("molecular P,T facade input must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected molecular P,T point, got {other:?}"),
        };
        let elemental_pt = match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)]).unwrap())
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(2_500.0)
            .production_cascade()
            .solve()
            .expect("elemental P,T facade input must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected elemental P,T point, got {other:?}"),
        };

        assert_eq!(
            molecular_pt.resolved().phase_specs(),
            elemental_pt.resolved().phase_specs()
        );
        assert_eq!(
            molecular_pt.solution().build_report().element_totals(),
            elemental_pt.solution().build_report().element_totals()
        );
        for (index, (molecular_moles, elemental_moles)) in molecular_pt
            .solution()
            .component_moles()
            .iter()
            .zip(elemental_pt.solution().component_moles())
            .enumerate()
        {
            assert!(
                (molecular_moles - elemental_moles).abs() < 1.0e-6,
                "PT component {index}: molecular={molecular_moles:?}, elemental={elemental_moles:?}"
            );
        }

        let thermochemistry =
            ResolvedThermochemistry::from_resolved_system(molecular_pt.resolved()).unwrap();
        let target_enthalpy = thermochemistry
            .enthalpy_model()
            .evaluate_total(molecular_pt.solution().component_moles(), 2_500.0)
            .unwrap();
        let bounds = TemperatureBounds::new(2_100.0, 2_900.0).unwrap();
        let mut ph_temperature_options = PhTemperatureSolveOptions::default();
        ph_temperature_options.scaled_enthalpy_tolerance = 1.0e-7;
        let molecular_ph = match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .initial_moles(vec![2.0, 1.0, 0.0])
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_total_enthalpy(
                TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
                2_300.0,
                bounds,
            )
            .ph_solve_mode(PhSolveMode::NestedTemperature)
            .ph_temperature_options(ph_temperature_options.clone())
            .production_cascade()
            .solve()
            .expect("molecular P,H facade input must solve")
        {
            EquilibriumCalculatorOutcome::PhPoint(point) => point,
            other => panic!("expected molecular P,H point, got {other:?}"),
        };
        let elemental_ph = match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)]).unwrap())
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_total_enthalpy(
                TotalEnthalpyJoules::new(target_enthalpy).unwrap(),
                2_300.0,
                bounds,
            )
            .ph_solve_mode(PhSolveMode::NestedTemperature)
            .ph_temperature_options(ph_temperature_options)
            .production_cascade()
            .solve()
            .expect("elemental P,H facade input must solve")
        {
            EquilibriumCalculatorOutcome::PhPoint(point) => point,
            other => panic!("expected elemental P,H point, got {other:?}"),
        };

        assert_eq!(
            molecular_ph.resolved().phase_specs(),
            elemental_ph.resolved().phase_specs()
        );
        assert_eq!(
            molecular_ph
                .solution()
                .equilibrium()
                .build_report()
                .element_totals(),
            elemental_ph
                .solution()
                .equilibrium()
                .build_report()
                .element_totals()
        );
        assert!(
            (molecular_ph.solution().temperature() - elemental_ph.solution().temperature()).abs()
                < 1.0e-4,
            "PH temperature: molecular={}, elemental={}",
            molecular_ph.solution().temperature(),
            elemental_ph.solution().temperature()
        );
        for (index, (molecular_moles, elemental_moles)) in molecular_ph
            .solution()
            .equilibrium()
            .component_moles()
            .iter()
            .zip(elemental_ph.solution().equilibrium().component_moles())
            .enumerate()
        {
            assert!(
                (molecular_moles - elemental_moles).abs() < 1.0e-6,
                "PH component {index}: molecular={molecular_moles:?}, elemental={elemental_moles:?}"
            );
        }
    }

    #[test]
    fn calculator_three_way_pt_metamorphic_inputs_with_same_b_match() {
        let solve_molecular = |initial_moles| match gas_point_builder()
            .initial_moles(initial_moles)
            .at_temperature(2_500.0)
            .solve()
            .expect("molecular P,T feed must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected molecular P,T point, got {other:?}"),
        };
        let molecular_reactants = solve_molecular(vec![2.0, 1.0, 0.0]);
        let molecular_mixed = solve_molecular(vec![1.0, 0.5, 1.0]);
        let elemental = match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)]).unwrap())
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(2_500.0)
            .production_cascade()
            .solve()
            .expect("elemental P,T feed must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected elemental P,T point, got {other:?}"),
        };

        let reference = molecular_reactants.solution();
        for (name, point) in [
            ("molecular-mixed", &molecular_mixed),
            ("elemental", &elemental),
        ] {
            let solution = point.solution();
            assert_eq!(
                point.resolved().phase_specs(),
                molecular_reactants.resolved().phase_specs()
            );
            assert_eq!(
                solution.metadata().components(),
                reference.metadata().components()
            );
            assert_eq!(solution.conditions(), reference.conditions());
            assert_eq!(solution.build_report().element_labels(), ["H", "O"]);
            assert_eq!(solution.build_report().element_totals(), [4.0, 2.0]);

            for component in reference.metadata().components() {
                let expected_fraction = reference
                    .mole_fraction_for(component.id())
                    .expect("reference component fraction must exist");
                let actual_fraction = solution
                    .mole_fraction_for(component.id())
                    .expect("candidate component fraction must exist");
                assert!(
                    (actual_fraction - expected_fraction).abs() < 1.0e-8,
                    "{name} fraction for {} differs: actual={actual_fraction:?}, expected={expected_fraction:?}",
                    component.label()
                );
            }
            for (index, (actual, expected)) in solution
                .component_moles()
                .iter()
                .zip(reference.component_moles())
                .enumerate()
            {
                assert!(
                    (actual - expected).abs() < 1.0e-6,
                    "{name} component {index}: actual={actual:?}, expected={expected:?}"
                );
            }

            let actual_validation = solution.accepted_solution().validation();
            let reference_validation = reference.accepted_solution().validation();
            assert!(
                actual_validation.residual_l2_norm < 1.0e-5
                    && reference_validation.residual_l2_norm < 1.0e-5
                    && (actual_validation.residual_l2_norm - reference_validation.residual_l2_norm)
                        .abs()
                        < 1.0e-5,
                "{name} residual quality differs too much: actual={}, expected={}",
                actual_validation.residual_l2_norm,
                reference_validation.residual_l2_norm
            );
            assert!(actual_validation.max_abs_element_balance_error / 4.0 < 1.0e-7);
        }
    }

    #[test]
    fn elemental_calculator_ph_distinct_targets_recover_distinct_temperatures() {
        let solve_pt = |temperature| match gas_point_builder()
            .at_temperature(temperature)
            .solve()
            .expect("reference P,T point must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected reference P,T point, got {other:?}"),
        };
        let low_reference = solve_pt(2_200.0);
        let high_reference = solve_pt(2_600.0);
        let thermochemistry =
            ResolvedThermochemistry::from_resolved_system(low_reference.resolved()).unwrap();
        let low_target = thermochemistry
            .enthalpy_model()
            .evaluate_total(low_reference.solution().component_moles(), 2_200.0)
            .unwrap();
        let high_target = thermochemistry
            .enthalpy_model()
            .evaluate_total(high_reference.solution().component_moles(), 2_600.0)
            .unwrap();
        assert!(high_target > low_target);

        let mut temperature_options = PhTemperatureSolveOptions::default();
        temperature_options.scaled_enthalpy_tolerance = 1.0e-7;
        let solve_elemental_ph =
            |target, initial_temperature| match EquilibriumCalculator::builder()
                .single_ideal_gas(["H2", "O2", "H2O"])
                .unwrap()
                .element_inventory(
                    ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)]).unwrap(),
                )
                .prefer_libraries(["NASA_gas"])
                .offline_only()
                .pressure_pa(PRESSURE_PA)
                .reference_pressure_pa(REFERENCE_PRESSURE_PA)
                .at_total_enthalpy(
                    TotalEnthalpyJoules::new(target).unwrap(),
                    initial_temperature,
                    TemperatureBounds::new(1_900.0, 2_900.0).unwrap(),
                )
                .ph_solve_mode(PhSolveMode::NestedTemperature)
                .ph_temperature_options(temperature_options.clone())
                .production_cascade()
                .solve()
                .expect("elemental P,H target must solve")
            {
                EquilibriumCalculatorOutcome::PhPoint(point) => point,
                other => panic!("expected elemental P,H point, got {other:?}"),
            };
        let low = solve_elemental_ph(low_target, 2_100.0);
        let high = solve_elemental_ph(high_target, 2_500.0);

        assert!(high.solution().temperature() > low.solution().temperature());
        assert!((high.solution().temperature() - low.solution().temperature()).abs() > 100.0);
        for point in [&low, &high] {
            assert_eq!(
                point.solution().equilibrium().build_report().input_kind(),
                PhaseEquilibriumInputKind::ElementInventory
            );
            assert_eq!(
                point
                    .solution()
                    .equilibrium()
                    .build_report()
                    .element_totals(),
                [4.0, 2.0]
            );
            assert!(
                point.solution().enthalpy_error().abs()
                    <= point.solution().enthalpy_error_limit_joules()
            );
            assert!(point.solution().scaled_enthalpy_error().abs() <= 1.0e-7);
        }
    }

    #[test]
    fn elemental_calculator_pt_is_extensive_under_inventory_scaling() {
        let solve = |scale| match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(
                ElementInventory::from_amounts([("H", 4.0 * scale), ("O", 2.0 * scale)]).unwrap(),
            )
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(2_500.0)
            .production_cascade()
            .solve()
            .expect("scaled elemental P,T feed must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected elemental P,T point, got {other:?}"),
        };
        let baseline = solve(1.0);
        let scaled = solve(10.0);

        assert_eq!(
            baseline.resolved().phase_specs(),
            scaled.resolved().phase_specs()
        );
        assert_eq!(baseline.solution().phases(), scaled.solution().phases());
        assert_eq!(
            baseline.solution().metadata().components(),
            scaled.solution().metadata().components()
        );
        assert_eq!(
            baseline.solution().build_report().element_totals(),
            [4.0, 2.0]
        );
        assert_eq!(
            scaled.solution().build_report().element_totals(),
            [40.0, 20.0]
        );
        for component in baseline.solution().metadata().components() {
            let baseline_fraction = baseline
                .solution()
                .mole_fraction_for(component.id())
                .unwrap();
            let scaled_fraction = scaled.solution().mole_fraction_for(component.id()).unwrap();
            assert!((scaled_fraction - baseline_fraction).abs() < 1.0e-8);
        }
        for (index, (baseline_moles, scaled_moles)) in baseline
            .solution()
            .component_moles()
            .iter()
            .zip(scaled.solution().component_moles())
            .enumerate()
        {
            assert!(
                (scaled_moles / 10.0 - baseline_moles).abs() < 1.0e-6,
                "component {index}: baseline={baseline_moles:?}, scaled={scaled_moles:?}"
            );
        }
        let baseline_balance = baseline
            .solution()
            .accepted_solution()
            .validation()
            .max_abs_element_balance_error
            / 4.0;
        let scaled_balance = scaled
            .solution()
            .accepted_solution()
            .validation()
            .max_abs_element_balance_error
            / 40.0;
        assert!(baseline_balance < 1.0e-7);
        assert!(scaled_balance < 1.0e-7);
        assert!((scaled_balance - baseline_balance).abs() < 1.0e-7);
    }

    #[test]
    fn elemental_calculator_ph_is_extensive_under_inventory_and_enthalpy_scaling() {
        let reference = match gas_point_builder()
            .at_temperature(2_500.0)
            .solve()
            .expect("reference P,T point must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected reference P,T point, got {other:?}"),
        };
        let thermochemistry =
            ResolvedThermochemistry::from_resolved_system(reference.resolved()).unwrap();
        let target_enthalpy = thermochemistry
            .enthalpy_model()
            .evaluate_total(reference.solution().component_moles(), 2_500.0)
            .unwrap();
        let mut temperature_options = PhTemperatureSolveOptions::default();
        temperature_options.scaled_enthalpy_tolerance = 1.0e-7;
        let solve = |scale| match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(
                ElementInventory::from_amounts([("H", 4.0 * scale), ("O", 2.0 * scale)]).unwrap(),
            )
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_total_enthalpy(
                TotalEnthalpyJoules::new(target_enthalpy * scale).unwrap(),
                2_300.0,
                TemperatureBounds::new(2_100.0, 2_900.0).unwrap(),
            )
            .ph_solve_mode(PhSolveMode::NestedTemperature)
            .ph_temperature_options(temperature_options.clone())
            .production_cascade()
            .solve()
            .expect("scaled elemental P,H feed must solve")
        {
            EquilibriumCalculatorOutcome::PhPoint(point) => point,
            other => panic!("expected elemental P,H point, got {other:?}"),
        };
        let baseline = solve(1.0);
        let scaled = solve(10.0);

        assert!(
            (scaled.solution().temperature() - baseline.solution().temperature()).abs() < 1.0e-4
        );
        assert_eq!(
            baseline.solution().equilibrium().phases(),
            scaled.solution().equilibrium().phases()
        );
        assert_eq!(
            baseline
                .solution()
                .equilibrium()
                .build_report()
                .element_totals(),
            [4.0, 2.0]
        );
        assert_eq!(
            scaled
                .solution()
                .equilibrium()
                .build_report()
                .element_totals(),
            [40.0, 20.0]
        );
        for component in baseline.solution().equilibrium().metadata().components() {
            let baseline_fraction = baseline
                .solution()
                .equilibrium()
                .mole_fraction_for(component.id())
                .unwrap();
            let scaled_fraction = scaled
                .solution()
                .equilibrium()
                .mole_fraction_for(component.id())
                .unwrap();
            assert!((scaled_fraction - baseline_fraction).abs() < 1.0e-8);
        }
        for (index, (baseline_moles, scaled_moles)) in baseline
            .solution()
            .equilibrium()
            .component_moles()
            .iter()
            .zip(scaled.solution().equilibrium().component_moles())
            .enumerate()
        {
            assert!(
                (scaled_moles / 10.0 - baseline_moles).abs() < 1.0e-6,
                "component {index}: baseline={baseline_moles:?}, scaled={scaled_moles:?}"
            );
        }
        for point in [&baseline, &scaled] {
            assert!(
                point.solution().enthalpy_error().abs()
                    <= point.solution().enthalpy_error_limit_joules()
            );
            assert!(point.solution().scaled_enthalpy_error().abs() <= 1.0e-7);
        }
    }

    #[test]
    fn explicit_and_from_elements_policies_have_different_universes_when_explicit_is_restricted() {
        let repository = ThermoData::try_default_repository()
            .expect("bundled repository must load for automatic universe selection");
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let candidate_policy = EquilibriumCandidatePolicy::default()
            .with_library_preference(vec!["NASA_gas".into()])
            .with_temperature_range(300.0, 1_200.0)
            .unwrap()
            .with_max_candidates(3)
            .unwrap();
        let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
            .select_inventory(&inventory, candidate_policy.clone())
            .expect("automatic candidate selection must succeed");
        assert!(
            selection
                .selected()
                .iter()
                .any(|candidate| candidate.record_key() == "H2O")
        );
        let phase_plan = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::ideal_gas(
                PhaseId::new(Some("gas".into())),
                selection
                    .selected()
                    .iter()
                    .map(|candidate| candidate.record_key().to_string())
                    .collect(),
            ),
        ]);

        let automatic = match EquilibriumCalculator::builder()
            .species_universe_policy(EquilibriumSpeciesUniversePolicy::FromElements {
                candidate_policy,
                phase_plan,
            })
            .element_inventory(inventory.clone())
            .with_repository(Arc::clone(&repository))
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(1_200.0)
            .solve()
            .expect("automatic universe must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected an automatic P,T point, got {other:?}"),
        };
        let explicit = match EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2"])
            .expect("restricted explicit phase must validate")
            .element_inventory(inventory)
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(1_200.0)
            .solve()
            .expect("restricted explicit universe must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected an explicit P,T point, got {other:?}"),
        };

        let automatic_components = automatic.resolved().phase_specs()[0].components();
        let explicit_components = explicit.resolved().phase_specs()[0].components();
        assert!(
            automatic_components
                .iter()
                .any(|component| component == "H2O")
        );
        assert_eq!(explicit_components, ["H2", "O2"]);
        assert_ne!(automatic_components, explicit_components);
        let water_index = automatic_components
            .iter()
            .position(|component| component == "H2O")
            .expect("automatic universe must retain H2O");
        assert!(automatic.solution().component_moles()[water_index] > 0.0);
        assert_eq!(explicit.solution().component_moles().len(), 2);
    }

    #[test]
    fn explicit_and_from_elements_policies_match_when_the_resolved_universe_matches() {
        let repository = ThermoData::try_default_repository()
            .expect("bundled repository must load for automatic universe selection");
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let candidate_policy = EquilibriumCandidatePolicy::default()
            .with_library_preference(vec!["NASA_gas".into()])
            .with_temperature_range(300.0, 1_200.0)
            .unwrap()
            .with_max_candidates(3)
            .unwrap();
        let selection = EquilibriumCandidateSelector::new(Arc::clone(&repository))
            .select_inventory(&inventory, candidate_policy.clone())
            .expect("automatic candidate selection must succeed");
        let selected_components = selection
            .selected()
            .iter()
            .map(|candidate| candidate.record_key().to_string())
            .collect::<Vec<_>>();
        let phase_plan = EquilibriumCandidatePhasePlan::new(vec![
            EquilibriumCandidatePhaseAssignment::ideal_gas(
                PhaseId::new(Some("gas".into())),
                selected_components.clone(),
            ),
        ]);

        let automatic = match EquilibriumCalculator::builder()
            .species_universe_policy(EquilibriumSpeciesUniversePolicy::FromElements {
                candidate_policy,
                phase_plan,
            })
            .element_inventory(inventory.clone())
            .with_repository(Arc::clone(&repository))
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(1_200.0)
            .solve()
            .expect("automatic universe must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected an automatic P,T point, got {other:?}"),
        };
        let explicit = match EquilibriumCalculator::builder()
            .ideal_gas_phase("gas", selected_components)
            .expect("matching explicit phase must validate")
            .element_inventory(inventory)
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(1_200.0)
            .solve()
            .expect("matching explicit universe must solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected an explicit P,T point, got {other:?}"),
        };

        assert_eq!(
            automatic.resolved().phase_specs(),
            explicit.resolved().phase_specs()
        );
        assert_eq!(
            automatic.solution().build_report().element_totals(),
            explicit.solution().build_report().element_totals()
        );
        assert_eq!(
            automatic.solution().component_moles().len(),
            explicit.solution().component_moles().len()
        );
        for (automatic_moles, explicit_moles) in automatic
            .solution()
            .component_moles()
            .iter()
            .zip(explicit.solution().component_moles())
        {
            assert!((automatic_moles - explicit_moles).abs() < 1.0e-8);
        }
        assert!(
            (automatic.solution().conditions().temperature()
                - explicit.solution().conditions().temperature())
            .abs()
                < 1.0e-10
        );
    }

    #[test]
    fn from_elements_policy_rejects_a_second_explicit_species_declaration() {
        let policy = EquilibriumCandidatePolicy::default();
        let phase_plan = EquilibriumCandidatePhasePlan::new(Vec::new());
        let error = EquilibriumCalculator::builder()
            .single_ideal_gas(["H2"])
            .expect("phase declaration must validate")
            .species_universe_policy(EquilibriumSpeciesUniversePolicy::FromElements {
                candidate_policy: policy,
                phase_plan,
            })
            .element_inventory(
                ElementInventory::from_amounts([("H", 2.0)])
                    .expect("element inventory must validate"),
            )
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(1_200.0)
            .solve()
            .expect_err("mixed universe declarations must be rejected");
        assert!(matches!(
            error,
            EquilibriumCalculatorError::ConflictingSpeciesUniverseDeclarations
        ));
    }

    #[test]
    fn facade_keeps_element_inventory_preparation_failures_out_of_solve_errors() {
        let inventory = ElementInventory::from_amounts([("C", 1.0), ("H", 2.0)])
            .expect("element inventory must validate");
        let result = EquilibriumCalculator::builder()
            .single_ideal_gas(["H2"])
            .expect("phase declaration must validate")
            .element_inventory(inventory)
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .at_temperature(1_200.0)
            .solve();

        assert!(matches!(
            result,
            Err(EquilibriumCalculatorError::Preparation(
                ReactionExtentError::Preparation(_)
            ))
        ));
    }

    #[test]
    fn calculator_error_types_keep_selection_preparation_and_solve_boundaries() {
        let selection = EquilibriumCalculatorError::from(CandidateSelectionError::UnknownLibrary(
            "missing".into(),
        ));
        assert!(matches!(
            selection,
            EquilibriumCalculatorError::CandidateSelection(_)
        ));

        let preparation = EquilibriumCalculatorError::from(PhRangeError::InvalidProblem(
            ReactionExtentError::Preparation(
                crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::
                    EquilibriumPreparationError::EmptySpeciesUniverse,
            ),
        ));
        assert!(matches!(
            preparation,
            EquilibriumCalculatorError::Preparation(ReactionExtentError::Preparation(_))
        ));

        let solve = EquilibriumCalculatorError::from(ReactionExtentError::InvalidCandidate {
            field: "candidate",
            message: "nonlinear candidate rejected".into(),
        });
        assert!(matches!(solve, EquilibriumCalculatorError::Solve(_)));
    }

    #[test]
    fn facade_accepts_element_inventory_for_pt_temperature_range() {
        let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])
            .expect("element inventory must validate");
        let result = EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .unwrap()
            .element_inventory(inventory)
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .over_temperature_range(vec![2_400.0, 2_500.0, 2_600.0])
            .unwrap()
            .solve()
            .expect("element-defined P,T range must solve");

        let EquilibriumCalculatorOutcome::PtRange(range) = result else {
            panic!("expected element-defined P,T range");
        };
        assert_eq!(range.solution().points().len(), 3);
        let balance_errors = range
            .solution()
            .points()
            .iter()
            .map(|point| {
                point
                    .solution()
                    .accepted_solution()
                    .validation()
                    .max_abs_element_balance_error
            })
            .collect::<Vec<_>>();
        let inventory_scale = 4.0_f64;
        assert!(
            balance_errors
                .iter()
                .all(|error| *error / inventory_scale < 1.0e-7),
            "unexpected range conservation residuals: {balance_errors:?}"
        );
        for point in range.solution().points() {
            assert_eq!(
                point.solution().build_report().input_kind(),
                PhaseEquilibriumInputKind::ElementInventory
            );
            assert_eq!(point.solution().build_report().element_labels(), ["H", "O"]);
            assert_eq!(point.solution().build_report().element_totals(), [4.0, 2.0]);
        }
        assert!(
            range.solution().points()[1]
                .report()
                .used_continuation_seed()
        );
    }

    #[test]
    fn facade_rejects_temperature_postprocessing_for_a_non_range_mode() {
        let policy = TemperaturePostprocessingPolicy::default();
        let error = gas_point_builder()
            .at_temperature(2_500.0)
            .temperature_postprocessing(policy)
            .solve()
            .expect_err("point solve must not silently ignore a range postprocessor");
        assert!(matches!(
            error,
            EquilibriumCalculatorError::PostprocessingRequiresPtRange
        ));
    }

    #[test]
    fn facade_returns_postprocessed_pt_range_without_mutating_raw_points() {
        let policy = TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::Uniform { points: 5 },
            ..TemperaturePostprocessingPolicy::default()
        };
        let range = match gas_point_builder()
            .over_temperature_range(vec![2_300.0, 2_500.0, 2_700.0])
            .expect("temperature grid must validate")
            .temperature_postprocessing(policy)
            .solve()
            .expect("facade P,T range must solve")
        {
            EquilibriumCalculatorOutcome::PtRange(range) => range,
            other => panic!("expected P,T range, got {other:?}"),
        };
        assert_eq!(range.solution().points().len(), 3);
        let postprocessed = range
            .postprocessing()
            .expect("requested postprocessing must be retained with the outcome");
        assert_eq!(postprocessed.raw.point_count(), 3);
        assert_eq!(
            postprocessed
                .resampled
                .as_ref()
                .expect("uniform policy must produce a separate resampled series")
                .point_count(),
            5
        );
    }

    #[test]
    fn facade_presets_change_observation_only() {
        let silent = match gas_point_builder()
            .preset(EquilibriumCalculatorPreset::Silent)
            .at_temperature(2_500.0)
            .solve()
            .expect("silent preset must preserve the solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected P,T point, got {other:?}"),
        };
        assert!(!silent.solution().timing_report().enabled());

        let timed = match gas_point_builder()
            .preset(EquilibriumCalculatorPreset::Timed)
            .at_temperature(2_500.0)
            .solve()
            .expect("timed preset must preserve the solve")
        {
            EquilibriumCalculatorOutcome::PtPoint(point) => point,
            other => panic!("expected P,T point, got {other:?}"),
        };
        assert!(timed.solution().timing_report().enabled());
        assert_eq!(
            silent.solution().component_moles().len(),
            timed.solution().component_moles().len()
        );
    }

    #[test]
    fn facade_accepts_phase_qualified_inventory_without_dense_order() {
        let phase = PhaseId::new(None);
        let outcome = gas_point_builder_without_inventory()
            .initial_phase_moles([
                (PhaseComponentId::new(phase.clone(), "O2"), 0.05),
                (PhaseComponentId::new(phase, "H2"), 0.1),
            ])
            .at_temperature(2_500.0)
            .solve()
            .expect("phase-qualified inventory must solve");
        let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
            panic!("expected P,T point");
        };
        assert_eq!(point.solution().component_moles().len(), 3);
        assert!(
            point
                .solution()
                .component_moles()
                .iter()
                .all(|mole| mole.is_finite())
        );
    }

    fn gas_point_builder_without_inventory() -> EquilibriumCalculatorBuilder {
        EquilibriumCalculator::builder()
            .single_ideal_gas(["H2", "O2", "H2O"])
            .expect("small gas phase declaration must validate")
            .prefer_libraries(["NASA_gas"])
            .offline_only()
            .pressure_pa(PRESSURE_PA)
            .reference_pressure_pa(REFERENCE_PRESSURE_PA)
            .production_cascade()
    }
}
