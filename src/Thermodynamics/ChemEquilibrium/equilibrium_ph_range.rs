//! Transactional fixed-pressure, fixed-enthalpy continuation.
//!
//! A `P,H` sweep is ordered by target enthalpy, not by temperature. Each
//! accepted point supplies its physical composition and solved temperature as
//! the seed for the next point. The underlying single-point workflow remains
//! the source of truth; this module owns only batch ordering, continuation,
//! timing, and publication semantics.

use std::error::Error;
use std::fmt;
use std::time::{Duration, Instant};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EquilibriumConstraint, TemperatureBounds, TotalEnthalpyJoules,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::EquilibriumRangeDiagnosticsPolicy;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
    FixedPressureEnthalpySolution, PhFallbackReason, PhSolveMode, PhSolvePath,
    PhTemperatureSolveOptions, PreparedNestedPhContinuationState, PreparedPhContinuationState,
    ResolvedPhaseEnthalpyRequest, ResolvedThermochemistry, solve_resolved_ph,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    EquilibriumSolveOptions, PhaseControlPolicy,
};
use crate::Thermodynamics::User_PhaseOrSolution::ResolvedPhaseSystem;

/// Direction of a validated target-enthalpy grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhRangeDirection {
    /// Targets increase in the requested continuation order.
    Ascending,
    /// Targets decrease in the requested continuation order.
    Descending,
}

/// Strictly monotone extensive enthalpy targets in joules.
#[derive(Debug, Clone, PartialEq)]
pub struct PhEnthalpyGrid {
    values: Vec<TotalEnthalpyJoules>,
    direction: PhRangeDirection,
}

impl PhEnthalpyGrid {
    /// Validates a non-empty, strictly monotone target grid.
    pub fn new(values: Vec<f64>) -> Result<Self, ReactionExtentError> {
        let mut typed = Vec::with_capacity(values.len());
        for value in values {
            typed.push(TotalEnthalpyJoules::new(value)?);
        }
        Self::from_joules(typed)
    }

    /// Validates an already unit-typed target grid.
    pub fn from_joules(values: Vec<TotalEnthalpyJoules>) -> Result<Self, ReactionExtentError> {
        if values.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_grid",
                message: "enthalpy grid must contain at least one point".into(),
            });
        }
        let direction = if values.len() == 1 {
            PhRangeDirection::Ascending
        } else if values
            .windows(2)
            .all(|pair| pair[1].joules() > pair[0].joules())
        {
            PhRangeDirection::Ascending
        } else if values
            .windows(2)
            .all(|pair| pair[1].joules() < pair[0].joules())
        {
            PhRangeDirection::Descending
        } else {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_grid",
                message: "enthalpy grid must be strictly ascending or descending".into(),
            });
        };
        Ok(Self { values, direction })
    }

    /// Targets in the requested continuation order.
    pub fn values(&self) -> &[TotalEnthalpyJoules] {
        &self.values
    }

    /// Numeric target values in joules.
    pub fn joules(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.values.iter().map(|value| value.joules())
    }

    /// Direction of the validated grid.
    pub fn direction(&self) -> PhRangeDirection {
        self.direction
    }
}

/// Why a batch point was prepared.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhRangePointPreparation {
    /// The first target used the caller-provided initial composition/seed.
    Initial,
    /// The previous accepted point supplied the continuation seed.
    Continued,
}

/// Immutable evidence attached to one accepted target-enthalpy point.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhRangePointReport {
    index: usize,
    target_enthalpy_bits: u64,
    seed_temperature_bits: u64,
    solved_temperature_bits: u64,
    preparation: PhRangePointPreparation,
    elapsed: Duration,
    phase_control_transitions: usize,
    formulation_builds: usize,
    formulation_reuses: usize,
    solve_path: PhSolvePath,
    fallback_reason: Option<PhFallbackReason>,
}

impl PhRangePointReport {
    /// Zero-based point index in the requested target order.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Target enthalpy in joules.
    pub fn target_enthalpy_joules(&self) -> f64 {
        f64::from_bits(self.target_enthalpy_bits)
    }

    /// Temperature used to initialize this point's solve.
    pub fn seed_temperature(&self) -> f64 {
        f64::from_bits(self.seed_temperature_bits)
    }

    /// Accepted equilibrium temperature.
    pub fn solved_temperature(&self) -> f64 {
        f64::from_bits(self.solved_temperature_bits)
    }

    /// Initial versus continuation preparation.
    pub fn preparation(&self) -> PhRangePointPreparation {
        self.preparation
    }

    /// Whether this point used the previous accepted point as its seed.
    pub fn used_continuation(&self) -> bool {
        self.preparation == PhRangePointPreparation::Continued
    }

    /// Wall time spent solving this point.
    pub fn elapsed(&self) -> Duration {
        self.elapsed
    }

    /// Accepted phase-control transitions at this target.
    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_transitions
    }

    /// Number of prepared formulation builds attributed to this point.
    pub fn formulation_builds(&self) -> usize {
        self.formulation_builds
    }

    /// Number of prepared formulation reuses attributed to this point.
    pub fn formulation_reuses(&self) -> usize {
        self.formulation_reuses
    }

    /// Numerical route that accepted this point.
    pub fn solve_path(&self) -> PhSolvePath {
        self.solve_path
    }

    /// Monolithic-to-nested fallback evidence, when `Auto` recovered.
    pub fn fallback_reason(&self) -> Option<&PhFallbackReason> {
        self.fallback_reason.as_ref()
    }
}

/// One accepted solution in a target-enthalpy sweep.
#[derive(Debug, Clone, PartialEq)]
pub struct PhRangePoint {
    solution: FixedPressureEnthalpySolution,
    report: PhRangePointReport,
}

impl PhRangePoint {
    /// Accepted immutable single-point result.
    pub fn solution(&self) -> &FixedPressureEnthalpySolution {
        &self.solution
    }

    /// Continuation and route evidence for this point.
    pub fn report(&self) -> &PhRangePointReport {
        &self.report
    }
}

/// Stable total/mean/median/worst timing summary for a successful sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PhRangeDurationSummary {
    total: Duration,
    mean: Duration,
    median: Duration,
    worst: Duration,
}

impl PhRangeDurationSummary {
    /// Total wall time spent in all accepted point solves.
    pub fn total(&self) -> Duration {
        self.total
    }

    /// Arithmetic mean of accepted point durations.
    pub fn mean(&self) -> Duration {
        self.mean
    }

    /// Median accepted point duration.
    pub fn median(&self) -> Duration {
        self.median
    }

    /// Slowest accepted point duration.
    pub fn worst(&self) -> Duration {
        self.worst
    }
}

/// Range-level continuation and timing evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhRangeSolveReport {
    direction: PhRangeDirection,
    point_count: usize,
    continuation_points: usize,
    phase_control_transitions: usize,
    formulation_builds: usize,
    formulation_reuses: usize,
    point_timing: PhRangeDurationSummary,
    total: Duration,
}

impl PhRangeSolveReport {
    /// Direction in which target enthalpies were solved.
    pub fn direction(&self) -> PhRangeDirection {
        self.direction
    }

    /// Number of accepted points published in the batch.
    pub fn point_count(&self) -> usize {
        self.point_count
    }

    /// Number of points seeded from a previous accepted solution.
    pub fn continuation_points(&self) -> usize {
        self.continuation_points
    }

    /// Total accepted phase-control transitions across the batch.
    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_transitions
    }

    /// Total prepared formulation builds across the range.
    pub fn formulation_builds(&self) -> usize {
        self.formulation_builds
    }

    /// Total prepared formulation reuses across the range.
    pub fn formulation_reuses(&self) -> usize {
        self.formulation_reuses
    }

    /// Aggregate point timing summary.
    pub fn point_timing(&self) -> PhRangeDurationSummary {
        self.point_timing
    }

    /// Wall time for the complete successful batch.
    pub fn total(&self) -> Duration {
        self.total
    }
}

/// Transactionally published target-enthalpy sweep.
#[derive(Debug, Clone, PartialEq)]
pub struct PhRangeSolution {
    points: Vec<PhRangePoint>,
    report: PhRangeSolveReport,
}

impl PhRangeSolution {
    /// Accepted points in the requested enthalpy order.
    pub fn points(&self) -> &[PhRangePoint] {
        &self.points
    }

    /// Immutable continuation and timing evidence for the batch.
    pub fn report(&self) -> &PhRangeSolveReport {
        &self.report
    }
}

/// A point-indexed failure from a transactional `P,H` sweep.
#[derive(Debug)]
pub struct PhRangePointError {
    index: usize,
    target_enthalpy: f64,
    source: ReactionExtentError,
}

impl PhRangePointError {
    /// Zero-based index of the failed target.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Failed target enthalpy in joules.
    pub fn target_enthalpy_joules(&self) -> f64 {
        self.target_enthalpy
    }

    /// Original typed error from the single-point P,H workflow.
    pub fn source(&self) -> &ReactionExtentError {
        &self.source
    }
}

impl fmt::Display for PhRangePointError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "P,H target point {} ({:.6e} J) failed: {}",
            self.index, self.target_enthalpy, self.source
        )
    }
}

impl Error for PhRangePointError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.source)
    }
}

/// Errors returned by a target-enthalpy sweep.
#[derive(Debug)]
pub enum PhRangeError {
    /// The target grid or request contract is invalid before solving starts.
    InvalidProblem(ReactionExtentError),
    /// A point failed; no partial `PhRangeSolution` is returned.
    Point(PhRangePointError),
}

impl fmt::Display for PhRangeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidProblem(error) => write!(formatter, "invalid P,H range: {error}"),
            Self::Point(error) => error.fmt(formatter),
        }
    }
}

impl Error for PhRangeError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidProblem(error) => Some(error),
            Self::Point(error) => Some(error),
        }
    }
}

/// Typed request for a fixed-pressure sweep over target enthalpies.
pub struct PhRangeRequest<'a> {
    prototype: ResolvedPhaseEnthalpyRequest<'a>,
    targets: PhEnthalpyGrid,
}

impl<'a> PhRangeRequest<'a> {
    /// Builds a production range request from one resolved thermochemistry
    /// bundle. The first target uses `initial_temperature`; later points use
    /// the previous accepted temperature and physical composition.
    pub fn from_resolved_thermochemistry(
        resolved: &'a ResolvedPhaseSystem,
        initial_composition: MultiphaseInitialComposition,
        pressure: f64,
        reference_pressure: f64,
        targets: PhEnthalpyGrid,
        temperature_bounds: TemperatureBounds,
        initial_temperature: f64,
        thermochemistry: ResolvedThermochemistry,
    ) -> Result<Self, PhRangeError> {
        let first_target = targets.values()[0];
        let constraint = EquilibriumConstraint::ph_joules(
            pressure,
            reference_pressure,
            first_target,
            initial_temperature,
        )
        .map_err(PhRangeError::InvalidProblem)?;
        let prototype = ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            resolved,
            initial_composition,
            constraint,
            temperature_bounds,
            thermochemistry,
        )
        .map_err(PhRangeError::InvalidProblem)?;
        Ok(Self { prototype, targets })
    }

    /// Selects the inner fixed-`P,T` backend policy for every point.
    pub fn with_solve_options(mut self, options: EquilibriumSolveOptions) -> Self {
        self.prototype = self.prototype.with_solve_options(options);
        self
    }

    /// Selects bounded phase control for every point.
    pub fn with_phase_control_policy(mut self, policy: PhaseControlPolicy) -> Self {
        self.prototype = self.prototype.with_phase_control_policy(policy);
        self
    }

    /// Selects monolithic, nested, or classified-auto P,H solving for every
    /// point.
    pub fn with_ph_solve_mode(mut self, mode: PhSolveMode) -> Self {
        self.prototype = self.prototype.with_ph_solve_mode(mode);
        self
    }

    /// Replaces scalar P,H controls shared by every target point.
    pub fn with_temperature_options(
        mut self,
        options: PhTemperatureSolveOptions,
    ) -> Result<Self, PhRangeError> {
        self.prototype = self
            .prototype
            .with_temperature_options(options)
            .map_err(PhRangeError::InvalidProblem)?;
        Ok(self)
    }

    /// Requested target grid.
    pub fn targets(&self) -> &PhEnthalpyGrid {
        &self.targets
    }

    /// Solves all targets transactionally with continuation.
    pub fn solve(self) -> Result<PhRangeSolution, PhRangeError> {
        let started = Instant::now();
        let layout =
            MultiphaseEquilibriumLayout::new(self.prototype.resolved().phase_specs().to_vec())
                .map_err(PhRangeError::InvalidProblem)?;
        let mut previous_composition = self.prototype.initial_composition().clone();
        let initial_seed = self
            .prototype
            .constraint()
            .initial_temperature()
            .ok_or_else(|| {
                PhRangeError::InvalidProblem(ReactionExtentError::InvalidProblem {
                    field: "constraint",
                    message: "P,H range requires a PH constraint".into(),
                })
            })?;
        let mut previous_temperature = initial_seed;
        let mut points = Vec::with_capacity(self.targets.values().len());
        let mut durations = Vec::with_capacity(self.targets.values().len());
        let mut transitions = 0usize;
        let mut formulation_builds = 0usize;
        let mut formulation_reuses = 0usize;
        let mut prepared_monolithic = if self.prototype.ph_solve_mode() == PhSolveMode::Monolithic
            && self.prototype.uses_fixed_declared_phases()
        {
            Some(
                PreparedPhContinuationState::new(&self.prototype)
                    .map_err(PhRangeError::InvalidProblem)?,
            )
        } else {
            None
        };
        let prepared_nested = if self.prototype.ph_solve_mode() == PhSolveMode::NestedTemperature
            && self.prototype.uses_fixed_declared_phases()
        {
            Some(
                PreparedNestedPhContinuationState::new(&self.prototype)
                    .map_err(PhRangeError::InvalidProblem)?,
            )
        } else {
            None
        };

        for (index, target) in self.targets.values().iter().copied().enumerate() {
            let seed_temperature = if index == 0 {
                initial_seed
            } else {
                previous_temperature
            };
            let preparation = if index == 0 {
                PhRangePointPreparation::Initial
            } else {
                PhRangePointPreparation::Continued
            };
            let mut request = self
                .prototype
                .clone()
                .with_target_enthalpy_and_seed(target, seed_temperature)
                .map_err(PhRangeError::InvalidProblem)?;
            // Diagnostics are intentionally sampled per target just like the
            // P,T range workflow. A live sink for a long P,H continuation
            // therefore observes only the caller-selected points instead of
            // receiving hundreds of otherwise identical lifecycle traces.
            let point_diagnostics = self
                .prototype
                .solve_options()
                .diagnostics_options()
                .for_range_point(index, self.targets.values().len());
            let point_options = request
                .solve_options()
                .clone()
                .with_diagnostics(point_diagnostics.clone());
            request = request.with_solve_options(point_options);
            if index > 0 {
                request = request
                    .with_initial_composition(previous_composition.clone())
                    .map_err(PhRangeError::InvalidProblem)?;
            }
            let point_started = Instant::now();
            let formulation_reused = prepared_monolithic.is_some() && index > 0;
            let mut solution = if let Some(state) = prepared_monolithic.as_mut() {
                state.solve(&request, formulation_reused)
            } else if let Some(state) = prepared_nested.as_ref() {
                state.solve(request)
            } else {
                solve_resolved_ph(request)
            }
            .map_err(|source| {
                PhRangeError::Point(PhRangePointError {
                    index,
                    target_enthalpy: target.joules(),
                    source,
                })
            })?;
            let elapsed = point_started.elapsed();
            let point_transitions = solution.report().phase_control_transitions();
            if self
                .prototype
                .solve_options()
                .diagnostics_options()
                .range_policy()
                == EquilibriumRangeDiagnosticsPolicy::TransitionsOnly
            {
                if point_transitions > 0 {
                    if let Some(report) = solution.equilibrium().diagnostics_report() {
                        point_diagnostics.replay_retained_events(report);
                    }
                } else {
                    solution = solution.without_equilibrium_diagnostics();
                }
            }
            previous_temperature = solution.temperature();
            previous_composition = MultiphaseInitialComposition::from_dense(
                &layout,
                solution.equilibrium().component_moles().to_vec(),
            )
            .map_err(PhRangeError::InvalidProblem)?;
            transitions += point_transitions;
            formulation_builds += solution.report().fixed_formulation_builds();
            formulation_reuses += solution.report().fixed_formulation_reuses();
            durations.push(elapsed);
            points.push(PhRangePoint {
                report: PhRangePointReport {
                    index,
                    target_enthalpy_bits: target.joules().to_bits(),
                    seed_temperature_bits: seed_temperature.to_bits(),
                    solved_temperature_bits: solution.temperature().to_bits(),
                    preparation,
                    elapsed,
                    phase_control_transitions: point_transitions,
                    formulation_builds: solution.report().fixed_formulation_builds(),
                    formulation_reuses: solution.report().fixed_formulation_reuses(),
                    solve_path: solution.report().solve_path(),
                    fallback_reason: solution.report().fallback_reason().cloned(),
                },
                solution,
            });
        }

        let report = PhRangeSolveReport {
            direction: self.targets.direction(),
            point_count: points.len(),
            continuation_points: points
                .iter()
                .filter(|point| point.report.used_continuation())
                .count(),
            phase_control_transitions: transitions,
            formulation_builds,
            formulation_reuses,
            point_timing: summarize_durations(&durations),
            total: started.elapsed(),
        };
        Ok(PhRangeSolution { points, report })
    }
}

fn summarize_durations(values: &[Duration]) -> PhRangeDurationSummary {
    if values.is_empty() {
        return PhRangeDurationSummary::default();
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let total = values.iter().copied().sum();
    let mean = total / values.len() as u32;
    let median = sorted[sorted.len() / 2];
    let worst = sorted[sorted.len() - 1];
    PhRangeDurationSummary {
        total,
        mean,
        median,
        worst,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enthalpy_grid_accepts_both_continuation_directions() {
        let ascending = PhEnthalpyGrid::new(vec![-10.0, 0.0, 25.0]).unwrap();
        assert_eq!(ascending.direction(), PhRangeDirection::Ascending);
        assert_eq!(
            ascending.joules().collect::<Vec<_>>(),
            vec![-10.0, 0.0, 25.0]
        );

        let descending = PhEnthalpyGrid::new(vec![25.0, 0.0, -10.0]).unwrap();
        assert_eq!(descending.direction(), PhRangeDirection::Descending);
    }

    #[test]
    fn enthalpy_grid_rejects_empty_duplicate_and_non_monotone_targets() {
        for values in [vec![], vec![1.0, 1.0], vec![1.0, 3.0, 2.0]] {
            assert!(matches!(
                PhEnthalpyGrid::new(values),
                Err(ReactionExtentError::InvalidProblem {
                    field: "enthalpy_grid",
                    ..
                })
            ));
        }
        assert!(PhEnthalpyGrid::new(vec![0.0, f64::NAN]).is_err());
    }

    #[test]
    fn point_error_keeps_target_index_and_typed_source() {
        let source = ReactionExtentError::InvalidProblem {
            field: "temperature",
            message: "outside bracket".into(),
        };
        let error = PhRangePointError {
            index: 2,
            target_enthalpy: 42.0,
            source,
        };
        assert_eq!(error.index(), 2);
        assert_eq!(error.target_enthalpy_joules(), 42.0);
        assert!(error.to_string().contains("point 2"));
        assert!(error.to_string().contains("4.200000e1"));
    }

    #[test]
    fn duration_summary_is_deterministic() {
        let summary = summarize_durations(&[
            Duration::from_millis(30),
            Duration::from_millis(10),
            Duration::from_millis(20),
        ]);
        assert_eq!(summary.total(), Duration::from_millis(60));
        assert_eq!(summary.mean(), Duration::from_millis(20));
        assert_eq!(summary.median(), Duration::from_millis(20));
        assert_eq!(summary.worst(), Duration::from_millis(30));
    }
}
