//! Typed fixed-pressure/temperature continuation over resolved phase data.
//!
//! This module is the production range facade. It resolves and prepares one
//! immutable phase layout, then refreshes only temperature-dependent
//! thermochemistry while using the previous accepted log-mole vector as the
//! next initial guess. The result is published only after every requested
//! point succeeds.

use std::time::{Duration, Instant};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_candidate_selection::EquilibriumCandidateSelectionReport;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::EquilibriumRangeDiagnosticsPolicy;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_element_inventory::ElementInventory;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::{
    EquilibriumExecutionControl, EquilibriumProgressEvent, EquilibriumProgressStage,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    DEFAULT_TRACE_MOLE_FLOOR, EquilibriumConditions, LogMolesInitialGuess, TraceSpeciesSeedPolicy,
};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::{
    EquilibriumTimingMode, EquilibriumTimingReport,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseSet;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    PhaseEquilibriumBuildRequest, SupportedPhaseModelPolicy,
    build_phase_equilibrium_problem_with_timing,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    EquilibriumSolveOptions, PhaseControlPolicy, ResolvedPhaseEquilibriumRequest,
    recover_resolved_pt_after_numerical_failure,
};
use crate::Thermodynamics::User_PhaseOrSolution::ResolvedPhaseSystem;

/// Direction of a validated temperature grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TemperatureRangeDirection {
    /// Temperatures increase from the first point to the last point.
    Ascending,
    /// Temperatures decrease from the first point to the last point.
    Descending,
}

/// Explicit temperature grid for one continuation request.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureGrid {
    values: Vec<f64>,
    direction: TemperatureRangeDirection,
}

impl TemperatureGrid {
    /// Validates a non-empty, strictly monotone finite grid.
    pub fn new(values: Vec<f64>) -> Result<Self, ReactionExtentError> {
        if values.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_grid",
                message: "temperature grid must contain at least one point".to_string(),
            });
        }
        if values
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_grid",
                message: "temperature grid values must be finite and strictly positive".to_string(),
            });
        }

        let direction = if values.len() == 1 {
            TemperatureRangeDirection::Ascending
        } else if values.windows(2).all(|pair| pair[1] > pair[0]) {
            TemperatureRangeDirection::Ascending
        } else if values.windows(2).all(|pair| pair[1] < pair[0]) {
            TemperatureRangeDirection::Descending
        } else {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_grid",
                message: "temperature grid must be strictly ascending or descending".to_string(),
            });
        };

        Ok(Self { values, direction })
    }

    /// Temperature values in the requested continuation order.
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Direction of the validated grid.
    pub fn direction(&self) -> TemperatureRangeDirection {
        self.direction
    }
}

/// Why a range point was prepared or reused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TemperatureRangePointPreparation {
    /// The first point created the shared formulation.
    InitialFormulation,
    /// The formulation was reused and only temperature-dependent state was
    /// refreshed before solving this point.
    ReusedFormulation,
    /// The cached physical formulation failed and the point was accepted by
    /// the canonical extensive-normalization recovery transaction.
    RecoveryFormulation,
}

/// Timing snapshot for one reduced active-set formulation retained by a
/// phase-control temperature-range solve.
///
/// The active mask is expressed in canonical phase order. A point report
/// contains the complete deterministic cache snapshot, rather than only the
/// aggregate time spent by that point, so a transition that creates several
/// formulations remains auditable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TemperatureRangeFormulationCacheTiming {
    active_mask: Vec<bool>,
    formulation_build: Duration,
}

impl TemperatureRangeFormulationCacheTiming {
    pub(crate) fn new(active_mask: Vec<bool>, formulation_build: Duration) -> Self {
        Self {
            active_mask,
            formulation_build,
        }
    }

    /// Active phases represented by this reduced formulation.
    pub fn active_mask(&self) -> &[bool] {
        &self.active_mask
    }

    /// Cumulative build time for this cache entry.
    pub fn formulation_build(&self) -> Duration {
        self.formulation_build
    }
}

/// Immutable timing and continuation evidence for one accepted point.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TemperatureRangePointReport {
    temperature_bits: u64,
    preparation: TemperatureRangePointPreparation,
    continuation: bool,
    thermochemistry_refreshed: bool,
    symbolic_parameter_reused: bool,
    phase_control_transitions: usize,
    phase_control_iterations: usize,
    phase_set_reused: bool,
    formulation_build: Duration,
    formulation_cache_timings: Vec<TemperatureRangeFormulationCacheTiming>,
    timing: EquilibriumTimingReport,
}

impl TemperatureRangePointReport {
    /// Temperature represented by this point report.
    pub fn temperature(&self) -> f64 {
        f64::from_bits(self.temperature_bits)
    }

    /// Preparation/reuse classification.
    pub fn preparation(&self) -> TemperatureRangePointPreparation {
        self.preparation
    }

    /// Whether the seed came from the previous accepted point.
    pub fn used_continuation_seed(&self) -> bool {
        self.continuation
    }

    /// Whether temperature-dependent thermochemistry was evaluated for this
    /// point's build report.
    pub fn thermochemistry_refreshed(&self) -> bool {
        self.thermochemistry_refreshed
    }

    /// Whether the RST symbolic problem was retained and its `T` parameter
    /// updated instead of rebuilding symbolic closures.
    pub fn symbolic_parameter_reused(&self) -> bool {
        self.symbolic_parameter_reused
    }

    /// Number of active-set transitions accepted at this temperature.
    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_transitions
    }

    /// Number of bounded outer-loop iterations used at this temperature.
    pub fn phase_control_iterations(&self) -> usize {
        self.phase_control_iterations
    }

    /// Whether the previous point's accepted phase set was used as the seed.
    pub fn phase_set_reused(&self) -> bool {
        self.phase_set_reused
    }

    /// Wall time spent rebuilding a symbolic/reduced formulation for this
    /// point. The initial formulation is reported by the range-level build
    /// timing; a zero value here means the point reused all prepared equation
    /// objects or timing was disabled.
    pub fn formulation_build(&self) -> Duration {
        self.formulation_build
    }

    /// Complete deterministic timing snapshot for reduced cache entries.
    pub fn formulation_cache_timings(&self) -> &[TemperatureRangeFormulationCacheTiming] {
        &self.formulation_cache_timings
    }

    /// Timing attached to this accepted point.
    pub fn timing(&self) -> &EquilibriumTimingReport {
        &self.timing
    }
}

/// Robust summary of accepted-point wall-clock durations.
///
/// Detailed stage evidence remains attached to every point. This value object
/// gives release characterization a stable total/mean/median/worst contract
/// without requiring console-log parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TemperatureRangeDurationSummary {
    total: Duration,
    mean: Duration,
    median: Duration,
    worst: Duration,
}

impl TemperatureRangeDurationSummary {
    pub fn total(&self) -> Duration {
        self.total
    }

    pub fn mean(&self) -> Duration {
        self.mean
    }

    pub fn median(&self) -> Duration {
        self.median
    }

    pub fn worst(&self) -> Duration {
        self.worst
    }
}

/// Immutable range-level characterization report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TemperatureRangeSolveReport {
    direction: TemperatureRangeDirection,
    point_count: usize,
    formulation_builds: usize,
    formulation_reuses: usize,
    symbolic_parameter_updates: usize,
    symbolic_problem_reused: bool,
    phase_control_enabled: bool,
    phase_projection_cache_entries: usize,
    phase_prepared_cache_entries: usize,
    phase_rst_cache_entries: usize,
    phase_stability_geometry_cache_entries: usize,
    phase_stability_geometry_cache_builds: usize,
    phase_stability_geometry_cache_reuses: usize,
    phase_control_transitions: usize,
    initial_formulation_timing: EquilibriumTimingReport,
    point_timing: TemperatureRangeDurationSummary,
    total: Duration,
}

impl TemperatureRangeSolveReport {
    /// Requested grid direction.
    pub fn direction(&self) -> TemperatureRangeDirection {
        self.direction
    }

    /// Number of transactionally accepted points.
    pub fn point_count(&self) -> usize {
        self.point_count
    }

    /// Number of full formulation builds.
    pub fn formulation_builds(&self) -> usize {
        self.formulation_builds
    }

    /// Number of points that reused the prepared formulation.
    pub fn formulation_reuses(&self) -> usize {
        self.formulation_reuses
    }

    /// Number of in-place symbolic temperature updates.
    pub fn symbolic_parameter_updates(&self) -> usize {
        self.symbolic_parameter_updates
    }

    /// Whether one RST symbolic problem was retained for the sweep.
    pub fn symbolic_problem_reused(&self) -> bool {
        self.symbolic_problem_reused
    }

    /// Whether the bounded phase-control outer loop was used.
    pub fn phase_control_enabled(&self) -> bool {
        self.phase_control_enabled
    }

    /// Number of distinct active-set projections retained by the sweep.
    pub fn phase_projection_cache_entries(&self) -> usize {
        self.phase_projection_cache_entries
    }

    /// Number of reduced prepared formulations retained across the sweep.
    pub fn phase_prepared_cache_entries(&self) -> usize {
        self.phase_prepared_cache_entries
    }

    /// Number of active-set entries retaining a prepared RST symbolic
    /// problem. This is separate from numeric preparation because symbolic
    /// construction is the expensive part of the default RST path.
    pub fn phase_rst_cache_entries(&self) -> usize {
        self.phase_rst_cache_entries
    }

    /// Number of retained immutable SVD/range/null-space geometries used by
    /// TPD phase-stability probes during the sweep.
    pub fn phase_stability_geometry_cache_entries(&self) -> usize {
        self.phase_stability_geometry_cache_entries
    }

    /// Number of times a TPD elemental geometry was constructed.
    pub fn phase_stability_geometry_cache_builds(&self) -> usize {
        self.phase_stability_geometry_cache_builds
    }

    /// Number of TPD probes that reused a previously constructed elemental
    /// geometry. This is diagnostic evidence, not a numerical iteration count.
    pub fn phase_stability_geometry_cache_reuses(&self) -> usize {
        self.phase_stability_geometry_cache_reuses
    }

    /// Total number of accepted phase transitions across all points.
    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_transitions
    }

    /// Timing accumulated while the first formulation was built.
    pub fn initial_formulation_timing(&self) -> &EquilibriumTimingReport {
        &self.initial_formulation_timing
    }

    /// Total/mean/median/worst wall-clock duration of accepted points.
    pub fn point_timing(&self) -> TemperatureRangeDurationSummary {
        self.point_timing
    }

    /// Wall-clock time for the complete transactional range operation.
    pub fn total(&self) -> Duration {
        self.total
    }
}

/// One accepted temperature point and its immutable phase-aware result.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureRangePoint {
    solution: MultiphaseEquilibriumSolution,
    report: TemperatureRangePointReport,
}

impl TemperatureRangePoint {
    /// Accepted phase-aware equilibrium result.
    pub fn solution(&self) -> &MultiphaseEquilibriumSolution {
        &self.solution
    }

    /// Continuation and timing evidence for this point.
    pub fn report(&self) -> &TemperatureRangePointReport {
        &self.report
    }
}

/// Transactionally published result of one typed temperature-range solve.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureRangeSolution {
    points: Vec<TemperatureRangePoint>,
    report: TemperatureRangeSolveReport,
}

impl TemperatureRangeSolution {
    /// Accepted points in exactly the requested grid order.
    pub fn points(&self) -> &[TemperatureRangePoint] {
        &self.points
    }

    /// Range-level timing and reuse evidence.
    pub fn report(&self) -> &TemperatureRangeSolveReport {
        &self.report
    }
}

/// Production request for a fixed-pressure, fixed-reference-pressure
/// temperature continuation over an already resolved phase system.
pub struct TemperatureRangeRequest<'a> {
    resolved: &'a ResolvedPhaseSystem,
    initial_composition: MultiphaseInitialComposition,
    element_inventory: Option<ElementInventory>,
    pressure: f64,
    reference_pressure: f64,
    temperatures: TemperatureGrid,
    model_policy: SupportedPhaseModelPolicy,
    solve_options: EquilibriumSolveOptions,
    phase_control_policy: Option<PhaseControlPolicy>,
    candidate_selection: Option<EquilibriumCandidateSelectionReport>,
}

impl<'a> TemperatureRangeRequest<'a> {
    /// Creates a validated range request without touching the repository.
    pub fn new(
        resolved: &'a ResolvedPhaseSystem,
        initial_composition: MultiphaseInitialComposition,
        pressure: f64,
        reference_pressure: f64,
        temperatures: TemperatureGrid,
    ) -> Result<Self, ReactionExtentError> {
        let first_conditions =
            EquilibriumConditions::new(temperatures.values()[0], pressure, reference_pressure)?;
        PhaseEquilibriumBuildRequest::new(
            resolved,
            first_conditions,
            initial_composition.clone(),
            TraceSpeciesSeedPolicy::Absolute {
                floor: DEFAULT_TRACE_MOLE_FLOOR,
            },
            SupportedPhaseModelPolicy::default(),
        )?;

        Ok(Self {
            resolved,
            initial_composition,
            element_inventory: None,
            pressure,
            reference_pressure,
            temperatures,
            model_policy: SupportedPhaseModelPolicy::default(),
            solve_options: EquilibriumSolveOptions::default(),
            phase_control_policy: None,
            candidate_selection: None,
        })
    }

    /// Creates a temperature range request from a closed elemental inventory.
    ///
    /// The first real-component seed is built once during request assembly;
    /// every prepared formulation and recovery attempt retains the original
    /// inventory as the physical conservation source.
    pub fn from_element_inventory(
        resolved: &'a ResolvedPhaseSystem,
        element_inventory: ElementInventory,
        pressure: f64,
        reference_pressure: f64,
        temperatures: TemperatureGrid,
    ) -> Result<Self, ReactionExtentError> {
        let first_conditions =
            EquilibriumConditions::new(temperatures.values()[0], pressure, reference_pressure)?;
        let bundle = build_phase_equilibrium_problem_with_timing(
            PhaseEquilibriumBuildRequest::from_element_inventory(
                resolved,
                first_conditions,
                element_inventory.clone(),
                TraceSpeciesSeedPolicy::Absolute {
                    floor: DEFAULT_TRACE_MOLE_FLOOR,
                },
                SupportedPhaseModelPolicy::default(),
            )?,
            EquilibriumTimingMode::Disabled,
        )?;
        let layout = crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::
            MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
        let initial_composition = MultiphaseInitialComposition::from_dense(
            &layout,
            bundle.problem().initial_moles().to_vec(),
        )?;

        Ok(Self {
            resolved,
            initial_composition,
            element_inventory: Some(element_inventory),
            pressure,
            reference_pressure,
            temperatures,
            model_policy: SupportedPhaseModelPolicy::default(),
            solve_options: EquilibriumSolveOptions::default(),
            phase_control_policy: None,
            candidate_selection: None,
        })
    }

    /// Replaces only the first numerical seed of an element-defined range.
    ///
    /// The closed elemental inventory remains the physical conservation
    /// source; the composition is retained solely as the first solver seed.
    pub fn with_initial_composition(
        mut self,
        composition: MultiphaseInitialComposition,
    ) -> Result<Self, ReactionExtentError> {
        let layout = MultiphaseEquilibriumLayout::new(self.resolved.phase_specs().to_vec())?;
        composition.validate_for(&layout)?;
        self.initial_composition = composition;
        Ok(self)
    }

    /// Replaces the supported activity-model policy.
    pub fn with_model_policy(mut self, policy: SupportedPhaseModelPolicy) -> Self {
        self.model_policy = policy;
        self
    }

    /// Replaces numerical and timing options for every point.
    pub fn with_solve_options(mut self, options: EquilibriumSolveOptions) -> Self {
        self.solve_options = options;
        self
    }

    /// Retains the immutable catalog-selection transaction in every range
    /// point's bridge report.
    pub fn with_candidate_selection(
        mut self,
        selection: EquilibriumCandidateSelectionReport,
    ) -> Self {
        self.candidate_selection = Some(selection);
        self
    }

    /// Enables bounded phase-control continuation for the range.
    ///
    /// The resolved layout and elemental inventory remain fixed. Only the
    /// accepted active phase set may change between points; such a transition
    /// is recorded and causes the corresponding projection to be prepared.
    pub fn with_phase_control_policy(mut self, policy: PhaseControlPolicy) -> Self {
        self.phase_control_policy = Some(policy);
        self
    }

    fn recover_failed_point(
        &self,
        conditions: EquilibriumConditions,
        seed: LogMolesInitialGuess,
        phase_set: Option<PhaseSet>,
        physical_error: ReactionExtentError,
    ) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
        let request = ResolvedPhaseEquilibriumRequest::new(
            self.resolved,
            conditions,
            self.initial_composition.clone(),
        )
        .with_model_policy(self.model_policy)
        .with_solve_options(self.solve_options.clone());
        let request = match &self.element_inventory {
            Some(inventory) => request.with_element_inventory(inventory.clone()),
            None => request,
        };
        let request = match &self.candidate_selection {
            Some(selection) => request.with_candidate_selection(Some(selection.clone())),
            None => request,
        };
        let request = match self.phase_control_policy.clone() {
            Some(policy) => request
                .with_phase_control_policy(policy)
                .with_continuation_state(seed, phase_set),
            None => request.with_multi_start_seeds(vec![seed]),
        };
        recover_resolved_pt_after_numerical_failure(request, physical_error)
    }

    /// Solves every point transactionally with continuation from the previous
    /// accepted point. A failure returns an error and publishes no range.
    pub fn solve(self) -> Result<TemperatureRangeSolution, ReactionExtentError> {
        let started = Instant::now();
        let execution_control: Option<EquilibriumExecutionControl> =
            self.solve_options.execution_control().cloned();
        if let Some(control) = &execution_control {
            control.check_cancelled()?;
        }
        let timing_mode = self.solve_options.timing_mode();
        let trace_policy = self.solve_options.trace_seed_policy();
        let first_conditions = EquilibriumConditions::new(
            self.temperatures.values()[0],
            self.pressure,
            self.reference_pressure,
        )?;
        let build_request = match &self.element_inventory {
            Some(inventory) => PhaseEquilibriumBuildRequest::from_element_inventory(
                self.resolved,
                first_conditions,
                inventory.clone(),
                trace_policy,
                self.model_policy,
            )?
            .with_numerical_seed(self.initial_composition.clone())?,
            None => PhaseEquilibriumBuildRequest::new(
                self.resolved,
                first_conditions,
                self.initial_composition.clone(),
                trace_policy,
                self.model_policy,
            )?,
        };
        let build_request = match &self.candidate_selection {
            Some(selection) => build_request.with_candidate_selection(selection.clone()),
            None => build_request,
        };
        let bundle = build_phase_equilibrium_problem_with_timing(build_request, timing_mode)?;
        if let Some(control) = &execution_control {
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::FormulationPreparation,
                None,
                Some(self.temperatures.values().len()),
                None,
            ));
            control.check_cancelled()?;
        }
        if let Some(phase_control_policy) = self.phase_control_policy.clone() {
            return self.solve_phase_control_range(
                bundle,
                phase_control_policy,
                started,
                timing_mode,
                trace_policy,
                execution_control,
            );
        }
        let mut template = bundle
            .into_temperature_template(self.solve_options.prepares_rst_backend(), timing_mode)?;
        let initial_formulation_timing = template.build_timing();
        let symbolic_problem_reused = template.symbolic_problem_reused();
        let mut seed = LogMolesInitialGuess::from_moles_with_policy(
            self.initial_composition.moles(),
            trace_policy,
        )?;
        let mut points = Vec::with_capacity(self.temperatures.values().len());

        for (index, &temperature) in self.temperatures.values().iter().enumerate() {
            if let Some(control) = &execution_control {
                control.check_cancelled()?;
                control.report(EquilibriumProgressEvent::new(
                    EquilibriumProgressStage::PointStarted,
                    Some(index),
                    Some(self.temperatures.values().len()),
                    Some(temperature),
                ));
            }
            let conditions =
                EquilibriumConditions::new(temperature, self.pressure, self.reference_pressure)?;
            let continuation = index > 0;
            let point_seed = seed.clone();
            let point_started = Instant::now();
            let (solution, recovered) = match template.solve_at(
                conditions,
                point_seed.clone(),
                self.solve_options.clone().into_settings(),
                timing_mode,
            ) {
                Ok(solution) => (solution, false),
                Err(error) => (
                    self.recover_failed_point(conditions, point_seed, None, error)
                        .map_err(|error| range_point_error(index, temperature, error))?,
                    true,
                ),
            };
            let solution = if recovered {
                solution.with_timing_total(point_started.elapsed())
            } else {
                solution
            };
            seed = LogMolesInitialGuess::new(solution.accepted_solution().log_moles().to_vec())?;
            points.push(TemperatureRangePoint {
                report: TemperatureRangePointReport {
                    temperature_bits: temperature.to_bits(),
                    preparation: if recovered {
                        TemperatureRangePointPreparation::RecoveryFormulation
                    } else if continuation {
                        TemperatureRangePointPreparation::ReusedFormulation
                    } else {
                        TemperatureRangePointPreparation::InitialFormulation
                    },
                    continuation,
                    thermochemistry_refreshed: true,
                    symbolic_parameter_reused: !recovered
                        && template.last_symbolic_parameter_reused(),
                    phase_control_transitions: 0,
                    phase_control_iterations: 0,
                    phase_set_reused: false,
                    formulation_build: if recovered {
                        recovery_formulation_duration(solution.timing_report())
                    } else {
                        template.last_formulation_build()
                    },
                    formulation_cache_timings: Vec::new(),
                    timing: *solution.timing_report(),
                },
                solution,
            });
            if let Some(control) = &execution_control {
                control.report(EquilibriumProgressEvent::new(
                    EquilibriumProgressStage::PointAccepted,
                    Some(index),
                    Some(self.temperatures.values().len()),
                    Some(temperature),
                ));
            }
        }

        let point_timing = summarize_point_timing(&points);

        Ok(TemperatureRangeSolution {
            report: TemperatureRangeSolveReport {
                direction: self.temperatures.direction(),
                point_count: points.len(),
                formulation_builds: 1 + points
                    .iter()
                    .filter(|point| {
                        point.report.preparation()
                            == TemperatureRangePointPreparation::RecoveryFormulation
                    })
                    .count(),
                formulation_reuses: points
                    .iter()
                    .filter(|point| {
                        point.report.preparation()
                            == TemperatureRangePointPreparation::ReusedFormulation
                    })
                    .count(),
                symbolic_parameter_updates: points
                    .iter()
                    .filter(|point| point.report.symbolic_parameter_reused())
                    .count(),
                symbolic_problem_reused,
                phase_control_enabled: false,
                phase_projection_cache_entries: 0,
                phase_prepared_cache_entries: 0,
                phase_rst_cache_entries: 0,
                phase_stability_geometry_cache_entries: 0,
                phase_stability_geometry_cache_builds: 0,
                phase_stability_geometry_cache_reuses: 0,
                phase_control_transitions: 0,
                initial_formulation_timing,
                point_timing,
                total: started.elapsed(),
            },
            points,
        })
    }

    /// Solves a temperature sweep with bounded active-set phase control.
    ///
    /// Each temperature point uses the [`PreparedPhaseControlRunner`] instead
    /// of the fixed-active-set path. The phase-control policy, trace seed
    /// policy, and execution control are forwarded to every point. Timing is
    /// collected when the timing mode is enabled. A single point failure
    /// returns an error and publishes no partial range result.
    fn solve_phase_control_range(
        self,
        bundle: crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::PhaseEquilibriumProblemBundle,
        phase_control_policy: PhaseControlPolicy,
        started: Instant,
        timing_mode: crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode,
        trace_policy: TraceSpeciesSeedPolicy,
        execution_control: Option<EquilibriumExecutionControl>,
    ) -> Result<TemperatureRangeSolution, ReactionExtentError> {
        let mut template = bundle.into_phase_control_template_with_diagnostics(
            |configured| *configured = phase_control_policy.into_phase_manager(),
            self.solve_options.diagnostics_options().clone(),
        )?;
        let initial_formulation_timing = template.build_timing();
        let mut seed = LogMolesInitialGuess::from_moles_with_policy(
            self.initial_composition.moles(),
            trace_policy,
        )?;
        let mut points = Vec::with_capacity(self.temperatures.values().len());
        let mut transitions = 0usize;

        for (index, &temperature) in self.temperatures.values().iter().enumerate() {
            if let Some(control) = &execution_control {
                control.check_cancelled()?;
                control.report(EquilibriumProgressEvent::new(
                    EquilibriumProgressStage::PointStarted,
                    Some(index),
                    Some(self.temperatures.values().len()),
                    Some(temperature),
                ));
            }
            let conditions =
                EquilibriumConditions::new(temperature, self.pressure, self.reference_pressure)?;
            let point_diagnostics = self
                .solve_options
                .diagnostics_options()
                .for_range_point(index, self.temperatures.values().len());
            template.set_diagnostics_options(point_diagnostics.clone());
            let continuation = index > 0;
            let phase_set = points.last().and_then(|point: &TemperatureRangePoint| {
                point
                    .solution()
                    .phase_control_report()
                    .map(|report| report.final_phase_set.clone())
            });
            let point_seed = seed.clone();
            let point_started = Instant::now();
            let (mut solution, recovered) = match template.solve_at(
                conditions,
                point_seed.clone(),
                self.solve_options.clone().into_settings(),
                timing_mode,
                phase_set.clone(),
            ) {
                Ok(solution) => (solution, false),
                Err(error) => (
                    self.recover_failed_point(conditions, point_seed, phase_set, error)
                        .map_err(|error| range_point_error(index, temperature, error))?,
                    true,
                ),
            };
            if recovered {
                solution = solution.with_timing_total(point_started.elapsed());
            }
            let (point_transitions, phase_control_iterations) = solution
                .phase_control_report()
                .map(|report| (report.transitions.len(), report.iterations))
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "temperature_range_phase_control",
                    message: "bounded range point did not publish phase-control evidence"
                        .to_string(),
                })?;
            if self.solve_options.diagnostics_options().range_policy()
                == EquilibriumRangeDiagnosticsPolicy::TransitionsOnly
            {
                if point_transitions > 0 {
                    if let Some(report) = solution.diagnostics_report() {
                        point_diagnostics.replay_retained_events(report);
                    }
                } else {
                    solution = solution.without_diagnostics();
                }
            }
            transitions += point_transitions;
            let formulation_cache_timings = template
                .formulation_cache_timings()
                .into_iter()
                .map(|(active_mask, duration)| {
                    TemperatureRangeFormulationCacheTiming::new(active_mask, duration)
                })
                .collect();
            seed = LogMolesInitialGuess::new(solution.accepted_solution().log_moles().to_vec())?;
            points.push(TemperatureRangePoint {
                report: TemperatureRangePointReport {
                    temperature_bits: temperature.to_bits(),
                    preparation: if recovered {
                        TemperatureRangePointPreparation::RecoveryFormulation
                    } else if continuation {
                        TemperatureRangePointPreparation::ReusedFormulation
                    } else {
                        TemperatureRangePointPreparation::InitialFormulation
                    },
                    continuation,
                    thermochemistry_refreshed: true,
                    symbolic_parameter_reused: !recovered && template.last_rst_symbolic_reused(),
                    phase_control_transitions: point_transitions,
                    phase_control_iterations,
                    phase_set_reused: continuation,
                    formulation_build: if recovered {
                        recovery_formulation_duration(solution.timing_report())
                    } else {
                        template.last_formulation_build()
                    },
                    formulation_cache_timings,
                    timing: *solution.timing_report(),
                },
                solution,
            });
            if let Some(control) = &execution_control {
                control.report(EquilibriumProgressEvent::new(
                    EquilibriumProgressStage::PointAccepted,
                    Some(index),
                    Some(self.temperatures.values().len()),
                    Some(temperature),
                ));
            }
        }

        let phase_stability_geometry_cache = template.phase_stability_geometry_cache_statistics();
        Ok(TemperatureRangeSolution {
            report: TemperatureRangeSolveReport {
                direction: self.temperatures.direction(),
                point_count: points.len(),
                formulation_builds: 1 + points
                    .iter()
                    .filter(|point| {
                        point.report.preparation()
                            == TemperatureRangePointPreparation::RecoveryFormulation
                    })
                    .count(),
                formulation_reuses: points
                    .iter()
                    .filter(|point| {
                        point.report.preparation()
                            == TemperatureRangePointPreparation::ReusedFormulation
                    })
                    .count(),
                symbolic_parameter_updates: points
                    .iter()
                    .filter(|point| point.report.symbolic_parameter_reused())
                    .count(),
                symbolic_problem_reused: template.rst_cache_size() > 0,
                phase_control_enabled: true,
                phase_projection_cache_entries: template.projection_cache_size(),
                phase_prepared_cache_entries: template.prepared_cache_size(),
                phase_rst_cache_entries: template.rst_cache_size(),
                phase_stability_geometry_cache_entries: phase_stability_geometry_cache.entries,
                phase_stability_geometry_cache_builds: phase_stability_geometry_cache.builds,
                phase_stability_geometry_cache_reuses: phase_stability_geometry_cache.reuses,
                phase_control_transitions: transitions,
                initial_formulation_timing,
                point_timing: summarize_point_timing(&points),
                total: started.elapsed(),
            },
            points,
        })
    }
}

fn summarize_point_timing(points: &[TemperatureRangePoint]) -> TemperatureRangeDurationSummary {
    if points.is_empty() {
        return TemperatureRangeDurationSummary::default();
    }
    let durations: Vec<_> = points
        .iter()
        .map(|point| point.report.timing.total())
        .collect();
    summarize_durations(&durations)
}

fn recovery_formulation_duration(report: &EquilibriumTimingReport) -> Duration {
    report.thermochemistry_preparation()
        + report.numeric_closure_construction()
        + report.symbolic_construction()
        + report.equation_construction()
        + report.numerical_problem_preparation()
        + report.projection_build()
}

fn summarize_durations(input: &[Duration]) -> TemperatureRangeDurationSummary {
    if input.is_empty() {
        return TemperatureRangeDurationSummary::default();
    }
    let mut durations = input.to_vec();
    durations.sort_unstable();
    let total_nanos: u128 = durations.iter().map(Duration::as_nanos).sum();
    let mean_nanos = total_nanos / durations.len() as u128;
    let median_nanos = if durations.len() % 2 == 1 {
        durations[durations.len() / 2].as_nanos()
    } else {
        let upper = durations.len() / 2;
        (durations[upper - 1].as_nanos() + durations[upper].as_nanos()) / 2
    };
    let to_duration = |nanos: u128| Duration::from_nanos(nanos.min(u64::MAX as u128) as u64);
    TemperatureRangeDurationSummary {
        total: to_duration(total_nanos),
        mean: to_duration(mean_nanos),
        median: to_duration(median_nanos),
        worst: durations.last().copied().unwrap_or_default(),
    }
}

fn range_point_error(
    index: usize,
    temperature: f64,
    error: ReactionExtentError,
) -> ReactionExtentError {
    // A release diagnostic needs the real attempt reports, not a string that
    // merely mentions them. The wrapper adds range coordinates without
    // destroying backend failure, termination, or timing evidence.
    ReactionExtentError::TemperatureRangePointFailed {
        point_index: index,
        temperature,
        cause: Box::new(error),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_preserves_ascending_and_descending_order() {
        let ascending = TemperatureGrid::new(vec![300.0, 500.0, 900.0]).unwrap();
        assert_eq!(ascending.direction(), TemperatureRangeDirection::Ascending);
        assert_eq!(ascending.values(), &[300.0, 500.0, 900.0]);

        let descending = TemperatureGrid::new(vec![900.0, 500.0, 300.0]).unwrap();
        assert_eq!(
            descending.direction(),
            TemperatureRangeDirection::Descending
        );
        assert_eq!(descending.values(), &[900.0, 500.0, 300.0]);
    }

    #[test]
    fn grid_rejects_duplicates_and_non_monotone_values() {
        assert!(TemperatureGrid::new(vec![300.0, 300.0]).is_err());
        assert!(TemperatureGrid::new(vec![300.0, 500.0, 400.0]).is_err());
        assert!(TemperatureGrid::new(Vec::new()).is_err());
    }

    #[test]
    fn duration_summary_reports_mean_median_and_worst() {
        let summary = summarize_durations(&[
            Duration::from_millis(1),
            Duration::from_millis(3),
            Duration::from_millis(8),
            Duration::from_millis(10),
        ]);
        assert_eq!(summary.total(), Duration::from_millis(22));
        assert_eq!(summary.mean(), Duration::from_micros(5500));
        assert_eq!(summary.median(), Duration::from_micros(5500));
        assert_eq!(summary.worst(), Duration::from_millis(10));
    }

    #[test]
    fn range_point_error_keeps_backend_attempt_diagnostics() {
        use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptReport;
        let error = range_point_error(
            0,
            1_000.0,
            ReactionExtentError::AllBackendsFailed {
                attempts: vec![SolverAttemptReport {
                    backend: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverBackend::RustedSciThe(
                        crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver::LevenbergMarquardt,
                    ),
                    outcome: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptOutcome::Failed {
                        kind: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptFailureKind::Solver,
                        reason: "step too small".to_string(),
                    },
                    metrics: None,
                }],
            },
        );

        let ReactionExtentError::TemperatureRangePointFailed {
            point_index,
            temperature,
            cause,
        } = error
        else {
            panic!("range point failure must preserve its typed source");
        };
        assert_eq!(point_index, 0);
        assert_eq!(temperature, 1_000.0);
        let ReactionExtentError::AllBackendsFailed { attempts } = cause.as_ref() else {
            panic!("range point source must retain the backend cascade");
        };
        assert_eq!(attempts.len(), 1);
        assert!(attempts[0].summary().contains("outcome=failed"));
        assert!(attempts[0].summary().contains("step too small"));
    }
}
