//! Stateless numerical helpers for the nested fixed-pressure/fixed-enthalpy path.
//!
//! This module deliberately contains no phase lifecycle, repository lookup,
//! `SubsData`, or result publication. The workflow owns those concerns and
//! supplies only sampled scalar evidence. Keeping these helpers independent
//! makes the bracketing contract testable without constructing a chemical
//! equilibrium problem. The workflow facade re-exports the public trial
//! snapshots, while the scalar engine and nested-specific records remain
//! owned here.

use std::time::{Duration, Instant};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EnthalpyScale, TemperatureBounds,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, MultiStartSolveReport,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingReport;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    MultiphaseAcceptanceReport, PhaseControlledSolveReport, PhaseStatus,
};
use crate::Thermodynamics::phase_layout::PhaseId;

/// Final lifecycle state of one semantic phase at an accepted temperature
/// trial. The phase id remains explicit because the same substance can occur
/// in distinct gas, liquid, or solid phases.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhTrialPhaseState {
    pub(crate) phase: PhaseId,
    pub(crate) status: PhaseStatus,
}

/// How the scalar outer solver selected one temperature trial.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhTemperatureStepKind {
    /// Lower user-provided temperature bound.
    LowerBound,
    /// Upper user-provided temperature bound.
    UpperBound,
    /// Explicit P,H seed inside the bounds.
    Seed,
    /// Accepted safeguarded secant interpolation inside the active bracket.
    Interpolation,
    /// Guaranteed-progress midpoint fallback.
    Bisection,
}

/// Structural preparation used by one accepted `P,H` temperature trial.
///
/// The value records a contract decision, not merely a performance counter:
/// bounded phase control is isolated because scalar bracket samples have no
/// monotone active-set history to carry safely between them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhTrialPreparation {
    /// Internal scalar-helper fixture with no nested equilibrium evidence.
    Unspecified,
    /// The first fixed-declared-phase trial built the reusable formulation.
    FixedFormulationInitial,
    /// A later fixed-declared-phase trial retargeted that formulation.
    FixedFormulationReused,
    /// Bounded phase control solved this trial as an isolated active-set
    /// lifecycle.
    BoundedPhaseControlIsolated,
}

impl PhTrialPhaseState {
    /// Semantic phase identity from the immutable resolved layout.
    pub fn phase(&self) -> &PhaseId {
        &self.phase
    }

    /// Lifecycle status accepted by the inner fixed-`P,T` solve.
    pub fn status(&self) -> PhaseStatus {
        self.status
    }
}

/// Complete immutable evidence from the accepted inner `P,T` solve for one
/// outer `P,H` temperature trial.
///
/// The outer report keeps its compact counters for tabular consumers, while
/// this snapshot keeps the detailed backend and active-set provenance needed
/// to explain a particular temperature point. `None` is reserved for pure
/// scalar test fixtures that do not execute an inner equilibrium solve.
#[derive(Debug, Clone, PartialEq)]
pub struct PhTrialInnerEvidence {
    pub(crate) solve_report: EquilibriumSolveReport,
    pub(crate) multi_start_report: Option<MultiStartSolveReport>,
    pub(crate) phase_control_report: Option<PhaseControlledSolveReport>,
    pub(crate) acceptance_report: Option<MultiphaseAcceptanceReport>,
}

impl PhTrialInnerEvidence {
    /// Ordered backend cascade which accepted this trial's final candidate.
    pub fn solve_report(&self) -> &EquilibriumSolveReport {
        &self.solve_report
    }

    /// Seed-level comparison evidence when continuation multi-start was used.
    pub fn multi_start_report(&self) -> Option<&MultiStartSolveReport> {
        self.multi_start_report.as_ref()
    }

    /// Transactional active-set lifecycle evidence for bounded phase control.
    pub fn phase_control_report(&self) -> Option<&PhaseControlledSolveReport> {
        self.phase_control_report.as_ref()
    }

    /// Final complementarity and numerical acceptance evidence, when present.
    pub fn acceptance_report(&self) -> Option<&MultiphaseAcceptanceReport> {
        self.acceptance_report.as_ref()
    }

    /// Sum of residual evaluations reported by all inner backend attempts.
    pub fn residual_evaluations(&self) -> usize {
        self.solve_report
            .attempts
            .iter()
            .filter_map(|attempt| attempt.metrics.as_ref())
            .map(|metrics| metrics.residual_evaluations)
            .sum()
    }

    /// Sum of Jacobian evaluations reported by all inner backend attempts.
    pub fn jacobian_evaluations(&self) -> usize {
        self.solve_report
            .attempts
            .iter()
            .filter_map(|attempt| attempt.metrics.as_ref())
            .map(|metrics| metrics.jacobian_evaluations)
            .sum()
    }
}

/// One accepted or rejected outer temperature trial.
#[derive(Debug, Clone, PartialEq)]
pub struct PhTemperatureTrial {
    pub(crate) temperature: f64,
    pub(crate) step_kind: PhTemperatureStepKind,
    pub(crate) total_enthalpy: f64,
    pub(crate) enthalpy_error_joules: f64,
    pub(crate) scaled_error: f64,
    pub(crate) inner_backend_attempts: usize,
    pub(crate) inner_nonlinear_iterations: usize,
    pub(crate) phase_control_transitions: usize,
    pub(crate) preparation: PhTrialPreparation,
    pub(crate) phase_states: Vec<PhTrialPhaseState>,
    pub(crate) inner_timing: EquilibriumTimingReport,
    pub(crate) timing: PhTrialTimingReport,
    pub(crate) inner_evidence: Option<PhTrialInnerEvidence>,
}

impl PhTemperatureTrial {
    /// Trial temperature in K.
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// Scalar-search step that selected this temperature.
    pub fn step_kind(&self) -> PhTemperatureStepKind {
        self.step_kind
    }

    /// Calculated total enthalpy in J.
    pub fn total_enthalpy(&self) -> f64 {
        self.total_enthalpy
    }

    /// Raw `H - H_target` in joules at this trial.
    pub fn enthalpy_error_joules(&self) -> f64 {
        self.enthalpy_error_joules
    }

    /// Dimensionless scaled `H - H_target`.
    pub fn scaled_error(&self) -> f64 {
        self.scaled_error
    }

    /// Number of inner nonlinear backend attempts used by this trial.
    pub fn inner_backend_attempts(&self) -> usize {
        self.inner_backend_attempts
    }

    /// Total nonlinear iterations spent by all inner backend attempts for
    /// this temperature trial, including continuation multi-start work.
    pub fn inner_nonlinear_iterations(&self) -> usize {
        self.inner_nonlinear_iterations
    }

    /// Number of accepted phase-control transitions in this trial.
    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_transitions
    }

    /// Structural preparation contract used by this trial.
    pub fn preparation(&self) -> PhTrialPreparation {
        self.preparation
    }

    /// Accepted lifecycle status of every semantic phase at this trial.
    pub fn phase_states(&self) -> &[PhTrialPhaseState] {
        &self.phase_states
    }

    /// Timing evidence published by the inner P,T solve.
    pub fn inner_timing(&self) -> EquilibriumTimingReport {
        self.inner_timing
    }

    /// Timing evidence for the complete scalar trial.
    ///
    /// This separates the outer trial wall time from the nested fixed-`P,T`
    /// solve and from the additive enthalpy evaluation. The report is zeroed
    /// when timing is disabled, so ordinary production solves do not acquire
    /// a hidden measurement contract.
    pub fn timing(&self) -> PhTrialTimingReport {
        self.timing
    }

    /// Complete backend and active-set trace for this trial's accepted inner
    /// equilibrium calculation. Pure scalar test fixtures have no such trace.
    pub fn inner_evidence(&self) -> Option<&PhTrialInnerEvidence> {
        self.inner_evidence.as_ref()
    }
}

/// Per-temperature timing evidence published with one `P,H` trial.
///
/// `inner_equilibrium` is the immutable timing report returned by the nested
/// fixed-`P,T` solve. `total` includes scalar-trial preparation, that solve,
/// and the enthalpy evaluation; the fields are intentionally not additive
/// wall-clock partitions because nested reports may include inclusive stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PhTrialTimingReport {
    pub(crate) enabled: bool,
    pub(crate) total: Duration,
    pub(crate) inner_equilibrium: Duration,
    pub(crate) enthalpy_evaluation: Duration,
}

impl PhTrialTimingReport {
    /// Whether this trial contains measured durations.
    pub fn enabled(&self) -> bool {
        self.enabled
    }

    /// Wall-clock duration of the complete trial.
    pub fn total(&self) -> Duration {
        self.total
    }

    /// Wall-clock duration reported by the nested fixed-`P,T` solve.
    pub fn inner_equilibrium(&self) -> Duration {
        self.inner_equilibrium
    }

    /// Time spent evaluating the additive enthalpy of the accepted inner
    /// composition.
    pub fn enthalpy_evaluation(&self) -> Duration {
        self.enthalpy_evaluation
    }
}



/// Scalar controls consumed by the nested bracket engine.
///
/// This is intentionally smaller than the public workflow options. Inner
/// backend budgets and phase-control policies are enforced by the trial
/// evaluator; the bracket engine sees only scalar acceptance and termination
/// rules.
#[derive(Debug, Clone, Copy)]
pub(crate) struct NestedBracketOptions {
    pub(crate) scaled_enthalpy_tolerance: f64,
    pub(crate) absolute_enthalpy_tolerance_joules: f64,
    pub(crate) temperature_tolerance: f64,
    pub(crate) max_iterations: usize,
    pub(crate) max_temperature_evaluations: usize,
    pub(crate) allow_interpolation: bool,
    pub(crate) allow_bracketed_sign_search: bool,
    pub(crate) max_wall_time: Option<Duration>,
}

impl NestedBracketOptions {
    /// Validates all scalar bracket controls before the first temperature trial.
    ///
    /// Checks that the enthalpy tolerance, absolute joule floor, evaluation
    /// budget, and optional wall-time limit are finite and positive. Returns
    /// a typed [`ReactionExtentError`] with the offending field name when any
    /// constraint is violated, so callers receive a diagnostic message instead
    /// of a silent default or a panic during the scalar search.
    fn validate(self) -> Result<(), ReactionExtentError> {
        if !self.scaled_enthalpy_tolerance.is_finite()
            || self.scaled_enthalpy_tolerance <= 0.0
        {
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
        if self.max_wall_time.is_some_and(|limit| limit.is_zero()) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "max_wall_time",
                message: "wall-time budget must be positive when provided".to_string(),
            });
        }
        Ok(())
    }

    /// Tests whether a raw enthalpy error satisfies the absolute-plus-relative
    /// acceptance contract used by the nested scalar-search path.
    ///
    /// The error is accepted when its absolute value does not exceed the larger
    /// of the absolute joule floor and the scale-aware relative tolerance.
    /// Non-finite errors are always rejected regardless of magnitude.
    fn accepts(&self, error_joules: f64, scale: EnthalpyScale) -> bool {
        error_joules.is_finite()
            && error_joules.abs()
                <= self
                    .absolute_enthalpy_tolerance_joules
                    .max(self.scaled_enthalpy_tolerance * scale.joules())
    }
}

/// Numerical step selected by the nested bracket engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NestedStepKind {
    LowerBound,
    UpperBound,
    Seed,
    Interpolation,
    Bisection,
}

/// Scalar evidence produced for one nested temperature evaluation.
#[derive(Debug, Clone, Copy)]
pub(crate) struct NestedTrial {
    pub(crate) temperature: f64,
    pub(crate) step_kind: NestedStepKind,
    pub(crate) total_enthalpy: f64,
    pub(crate) enthalpy_error_joules: f64,
    pub(crate) scaled_error: f64,
}

impl From<NestedTrial> for PhTemperatureTrial {
    fn from(trial: NestedTrial) -> Self {
        let step_kind = match trial.step_kind {
            NestedStepKind::LowerBound => PhTemperatureStepKind::LowerBound,
            NestedStepKind::UpperBound => PhTemperatureStepKind::UpperBound,
            NestedStepKind::Seed => PhTemperatureStepKind::Seed,
            NestedStepKind::Interpolation => PhTemperatureStepKind::Interpolation,
            NestedStepKind::Bisection => PhTemperatureStepKind::Bisection,
        };
        Self {
            temperature: trial.temperature,
            step_kind,
            total_enthalpy: trial.total_enthalpy,
            enthalpy_error_joules: trial.enthalpy_error_joules,
            scaled_error: trial.scaled_error,
            inner_backend_attempts: 0,
            inner_nonlinear_iterations: 0,
            phase_control_transitions: 0,
            preparation: PhTrialPreparation::Unspecified,
            phase_states: Vec::new(),
            inner_timing: EquilibriumTimingReport::default(),
            timing: PhTrialTimingReport::default(),
            inner_evidence: None,
        }
    }
}

/// Accepted value and scalar evidence returned by the nested engine.
pub(crate) struct NestedBracketResult<T> {
    pub(crate) value: T,
    pub(crate) total_enthalpy: f64,
    pub(crate) iterations: usize,
    pub(crate) trials: Vec<NestedTrial>,
}

/// Runs a safeguarded scalar temperature search around an arbitrary trial
/// evaluator. The evaluator remains responsible for the inner P,T solve;
/// this function owns only bracket invariants and scalar budgets.
pub(crate) fn solve_bracketed_temperature<T, F>(
    bounds: TemperatureBounds,
    scale: EnthalpyScale,
    target: f64,
    seed_temperature: Option<f64>,
    options: NestedBracketOptions,
    started: Instant,
    mut evaluate: F,
) -> Result<NestedBracketResult<T>, ReactionExtentError>
where
    F: FnMut(f64) -> Result<(T, f64), ReactionExtentError>,
{
    options.validate()?;
    let mut trials = Vec::new();
    let mut sampled_errors = Vec::<(f64, f64)>::new();
    let ensure_temperature_budget = |evaluations: usize| {
        if evaluations >= options.max_temperature_evaluations {
            Err(ReactionExtentError::InvalidProblem {
                field: "temperature_budget",
                message: format!(
                    "outer P,H solve exceeded the global temperature evaluation budget of {}",
                    options.max_temperature_evaluations
                ),
            })
        } else {
            ensure_wall_time_budget(started, options.max_wall_time)
        }
    };
    let trial_error = |index: usize, temperature: f64, cause: ReactionExtentError| {
        if matches!(&cause, ReactionExtentError::Cancelled) {
            cause
        } else {
            ReactionExtentError::TemperatureTrialFailed {
                trial_index: index,
                temperature,
                cause: Box::new(cause),
            }
        }
    };
    let make_trial = |temperature, step_kind, enthalpy| {
        let enthalpy_error_joules = enthalpy - target;
        let scaled_error = scale.scale_error(enthalpy_error_joules)?;
        Ok::<NestedTrial, ReactionExtentError>(NestedTrial {
            temperature,
            step_kind,
            total_enthalpy: enthalpy,
            enthalpy_error_joules,
            scaled_error,
        })
    };

    ensure_temperature_budget(trials.len())?;
    let lower_temperature = bounds.lower();
    let (lower_value, lower_enthalpy) = evaluate(lower_temperature)
        .map_err(|cause| trial_error(0, lower_temperature, cause))?;
    ensure_wall_time_budget(started, options.max_wall_time)?;
    let lower_trial = make_trial(
        lower_temperature,
        NestedStepKind::LowerBound,
        lower_enthalpy,
    )?;
    let lower_error = lower_trial.scaled_error;
    sampled_errors.push((lower_temperature, lower_error));
    trials.push(lower_trial);
    if options.accepts(lower_trial.enthalpy_error_joules, scale) {
        return Ok(NestedBracketResult {
            value: lower_value,
            total_enthalpy: lower_enthalpy,
            iterations: 0,
            trials,
        });
    }

    ensure_temperature_budget(trials.len())?;
    let upper_temperature = bounds.upper();
    let (upper_value, upper_enthalpy) = evaluate(upper_temperature)
        .map_err(|cause| trial_error(1, upper_temperature, cause))?;
    ensure_wall_time_budget(started, options.max_wall_time)?;
    let upper_trial = make_trial(
        upper_temperature,
        NestedStepKind::UpperBound,
        upper_enthalpy,
    )?;
    let upper_error = upper_trial.scaled_error;
    sampled_errors.push((upper_temperature, upper_error));
    trials.push(upper_trial);
    if options.accepts(upper_trial.enthalpy_error_joules, scale) {
        return Ok(NestedBracketResult {
            value: upper_value,
            total_enthalpy: upper_enthalpy,
            iterations: 0,
            trials,
        });
    }

    let mut lower = lower_temperature;
    let mut upper = upper_temperature;
    let mut lower_error = lower_error;
    let mut upper_error = upper_error;
    if let Some(seed) = seed_temperature
        .filter(|seed| *seed > bounds.lower() && *seed < bounds.upper() && seed.is_finite())
    {
        ensure_temperature_budget(trials.len())?;
        let trial_index = trials.len();
        let (seed_value, seed_enthalpy) =
            evaluate(seed).map_err(|cause| trial_error(trial_index, seed, cause))?;
        ensure_wall_time_budget(started, options.max_wall_time)?;
        let seed_trial = make_trial(seed, NestedStepKind::Seed, seed_enthalpy)?;
        let seed_error = seed_trial.scaled_error;
        sampled_errors.push((seed, seed_error));
        trials.push(seed_trial);
        if options.accepts(seed_trial.enthalpy_error_joules, scale) {
            return Ok(NestedBracketResult {
                value: seed_value,
                total_enthalpy: seed_enthalpy,
                iterations: 0,
                trials,
            });
        }

        let endpoints_have_same_sign = lower_error.signum() == upper_error.signum();
        let seed_differs_from_lower = seed_error.signum() != lower_error.signum();
        let seed_differs_from_upper = seed_error.signum() != upper_error.signum();
        if endpoints_have_same_sign && seed_differs_from_lower {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_multiple_brackets",
                message: format!(
                    "seed {seed} K creates two sign-change brackets inside [{}, {}] K",
                    bounds.lower(),
                    bounds.upper()
                ),
            });
        }
        if !endpoints_have_same_sign {
            if seed_differs_from_lower {
                upper = seed;
                upper_error = seed_error;
            } else if seed_differs_from_upper {
                lower = seed;
                lower_error = seed_error;
            }
        }
        validate_sampled_monotonicity(&sampled_errors, options.allow_bracketed_sign_search)?;
    }

    if lower_error.signum() == upper_error.signum() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "enthalpy_bracket",
            message: format!(
                "target enthalpy is not bracketed: scaled errors are {lower_error:e} and {upper_error:e}"
            ),
        });
    }

    for iteration in 1..=options.max_iterations {
        let midpoint = lower + (upper - lower) * 0.5;
        let (trial_temperature, step_kind) = if options.allow_interpolation {
            safeguarded_interpolation_step(lower, lower_error, upper, upper_error).map_or(
                (midpoint, NestedStepKind::Bisection),
                |temperature| (temperature, NestedStepKind::Interpolation),
            )
        } else {
            (midpoint, NestedStepKind::Bisection)
        };
        ensure_temperature_budget(trials.len())?;
        let trial_index = trials.len();
        let (value, enthalpy) = evaluate(trial_temperature)
            .map_err(|cause| trial_error(trial_index, trial_temperature, cause))?;
        ensure_wall_time_budget(started, options.max_wall_time)?;
        let trial = make_trial(trial_temperature, step_kind, enthalpy)?;
        let error = trial.scaled_error;
        sampled_errors.push((trial_temperature, error));
        validate_sampled_monotonicity(&sampled_errors, options.allow_bracketed_sign_search)?;
        let error_joules = trial.enthalpy_error_joules;
        trials.push(trial);
        if options.accepts(error_joules, scale) {
            return Ok(NestedBracketResult {
                value,
                total_enthalpy: enthalpy,
                iterations: iteration,
                trials,
            });
        }
        if (upper - lower).abs() <= options.temperature_tolerance {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_acceptance",
                message: format!(
                    "temperature bracket reached tolerance but enthalpy error {error_joules:e} J exceeds the accepted limit {} J",
                    options
                        .absolute_enthalpy_tolerance_joules
                        .max(options.scaled_enthalpy_tolerance * scale.joules())
                ),
            });
        }
        if error.signum() == lower_error.signum() {
            lower = trial_temperature;
            lower_error = error;
        } else {
            upper = trial_temperature;
            upper_error = error;
        }
    }

    Err(ReactionExtentError::InvalidProblem {
        field: "enthalpy_solver",
        message: format!(
            "temperature bracket did not converge in {} iterations",
            options.max_iterations
        ),
    })
}

/// Proposes a secant point only when it remains comfortably inside the active
/// sign bracket. Otherwise the caller must use bisection, preserving progress
/// near flat, discontinuous, or ill-scaled branches.
pub(crate) fn safeguarded_interpolation_step(
    lower_temperature: f64,
    lower_error: f64,
    upper_temperature: f64,
    upper_error: f64,
) -> Option<f64> {
    let width = upper_temperature - lower_temperature;
    let denominator = upper_error - lower_error;
    if !width.is_finite() || width <= 0.0 || !denominator.is_finite() || denominator == 0.0 {
        return None;
    }
    let candidate = upper_temperature - upper_error * width / denominator;
    let guard = width * 0.1;
    (candidate.is_finite()
        && candidate > lower_temperature + guard
        && candidate < upper_temperature - guard)
        .then_some(candidate)
}

/// Enforces the outer wall-time contract after an evaluator returns and before
/// another nested trial starts. An inner P,T solve may cross the deadline while
/// it is running, so checking only before evaluation is insufficient.
pub(crate) fn ensure_wall_time_budget(
    started: Instant,
    limit: Option<Duration>,
) -> Result<(), ReactionExtentError> {
    if limit.is_some_and(|limit| started.elapsed() >= limit) {
        Err(ReactionExtentError::InvalidProblem {
            field: "wall_time_budget",
            message: "outer P,H solve exceeded its wall-time budget".to_string(),
        })
    } else {
        Ok(())
    }
}

/// Checks only the scalar evidence sampled by the safeguarded search.
///
/// This cannot prove global monotonicity. It does prevent the solver from
/// silently treating an observed reversal as a smooth single-root branch.
pub(crate) fn validate_sampled_monotonicity(
    samples: &[(f64, f64)],
    allow_bracketed_sign_search: bool,
) -> Result<(), ReactionExtentError> {
    if allow_bracketed_sign_search || samples.len() < 3 {
        return Ok(());
    }

    let mut ordered = samples.to_vec();
    ordered.sort_by(|left, right| left.0.total_cmp(&right.0));
    let mut previous_slope_sign = None;
    for window in ordered.windows(2) {
        let delta_temperature = window[1].0 - window[0].0;
        let delta_error = window[1].1 - window[0].1;
        if delta_temperature <= 0.0 || delta_error == 0.0 {
            continue;
        }
        let slope_sign = delta_error.signum();
        if previous_slope_sign.is_some_and(|previous| previous != slope_sign) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "enthalpy_non_monotone",
                message: format!(
                    "sampled outer enthalpy errors reverse direction near {} K",
                    window[0].0
                ),
            });
        }
        previous_slope_sign = Some(slope_sign);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolation_stays_inside_a_guarded_sign_bracket() {
        let candidate = safeguarded_interpolation_step(300.0, -2.0, 700.0, 2.0);
        assert_eq!(candidate, Some(500.0));
        assert!(safeguarded_interpolation_step(300.0, 0.0, 700.0, 0.0).is_none());
    }

    #[test]
    fn monotonicity_guard_rejects_a_sampled_reversal_only_when_requested() {
        let samples = [(300.0, -1.0), (500.0, 1.0), (700.0, -1.0)];
        assert!(validate_sampled_monotonicity(&samples, false).is_err());
        assert!(validate_sampled_monotonicity(&samples, true).is_ok());
    }

    #[test]
    fn wall_time_guard_reports_a_typed_budget_error() {
        let error = ensure_wall_time_budget(Instant::now(), Some(Duration::ZERO)).unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "wall_time_budget",
                ..
            }
        ));
    }

    #[test]
    fn scalar_trial_conversion_keeps_step_and_defers_inner_evidence() {
        let trial = PhTemperatureTrial::from(NestedTrial {
            temperature: 812.5,
            step_kind: NestedStepKind::Interpolation,
            total_enthalpy: 42.0,
            enthalpy_error_joules: -0.25,
            scaled_error: -0.005,
        });

        assert_eq!(trial.temperature(), 812.5);
        assert_eq!(trial.step_kind(), PhTemperatureStepKind::Interpolation);
        assert_eq!(trial.total_enthalpy(), 42.0);
        assert_eq!(trial.enthalpy_error_joules(), -0.25);
        assert_eq!(trial.scaled_error(), -0.005);
        assert_eq!(trial.preparation(), PhTrialPreparation::Unspecified);
        assert!(trial.phase_states().is_empty());
        assert!(trial.inner_evidence().is_none());
    }
}
