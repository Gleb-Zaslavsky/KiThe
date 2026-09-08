//! Read-only presentation views for fixed-pressure temperature sweeps.
//!
//! This module turns the transactional `P,T` range result into stable labelled
//! series and point diagnostics. It deliberately preserves raw accepted points
//! and phase-transition boundaries; interpolation remains a presentation-only
//! operation and is rejected across a reported transition.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
    TemperaturePostprocessingPolicy, TemperaturePostprocessingResult, TemperatureResamplingGrid,
    TemperatureSweepSeries,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::TemperatureRangeSolution;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use std::time::Duration;

/// Presentation row for one accepted `P,T` range point.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureRangePresentationPointRow {
    /// Original point position in the caller-requested temperature grid.
    pub point_index: usize,
    /// Accepted temperature in K.
    pub temperature_kelvin: f64,
    /// Whether the previous accepted point provided the composition seed.
    pub used_continuation_seed: bool,
    /// Whether temperature-dependent thermochemistry was refreshed.
    pub thermochemistry_refreshed: bool,
    /// Whether the RST symbolic problem was retained and retargeted.
    pub symbolic_parameter_reused: bool,
    /// Number of accepted phase lifecycle transitions at this point.
    pub phase_control_transitions: usize,
    /// Number of phase-control outer-loop iterations at this point.
    pub phase_control_iterations: usize,
    /// Whether the previous accepted active phase set was reused.
    pub phase_set_reused: bool,
    /// Time attributed to rebuilding a reduced formulation for this point.
    pub formulation_build: Duration,
    /// Complete per-point timing snapshot when timing is enabled.
    pub total_duration: Duration,
    /// Backend that accepted the point.
    pub accepted_backend: String,
    /// Accepted scaled residual L2 norm.
    pub residual_l2_norm: f64,
    /// Maximum absolute elemental-balance error for the physical result.
    pub max_abs_element_balance_error: f64,
}

/// Boundary between adjacent ordered temperatures that must not be smoothed.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureRangeTransitionBoundary {
    /// Original point index below the boundary in temperature order.
    pub lower_point_index: usize,
    /// Original point index above the boundary in temperature order.
    pub upper_point_index: usize,
    /// Lower accepted temperature in K.
    pub lower_temperature_kelvin: f64,
    /// Upper accepted temperature in K.
    pub upper_temperature_kelvin: f64,
    /// Explicit phase-control transitions recorded at the upper point.
    pub transitions_at_upper_point: usize,
    /// Whether at least one explicit phase lifecycle status changed.
    pub phase_status_changed: bool,
}

/// Complete presentation projection of one accepted fixed-pressure temperature range.
///
/// Point diagnostics and series are ordered by increasing temperature, while
/// `point_index` retains the caller's original grid order. This keeps plotting
/// deterministic even for descending continuation runs without hiding which
/// accepted point supplied a seed or triggered a phase transition.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureRangePresentationReport {
    /// Point-level solver, continuation, and validation evidence.
    pub points: Vec<TemperatureRangePresentationPointRow>,
    /// Physical component mole amounts in canonical component order.
    pub component_moles: TemperatureSweepSeries,
    /// Local phase mole fractions in the same component order as `component_moles`.
    pub component_mole_fractions: TemperatureSweepSeries,
    /// Physical phase totals in canonical phase order.
    pub phase_totals: TemperatureSweepSeries,
    /// Residual, elemental-balance, and wall-time evidence at raw accepted points.
    /// Columns are `residual_l2_norm`, `max_abs_element_balance_error`, and
    /// `point_elapsed_ms` in that order.
    pub solver_metrics: TemperatureSweepSeries,
    /// Temperature boundaries where continuous interpolation is not valid.
    pub transition_boundaries: Vec<TemperatureRangeTransitionBoundary>,
}

impl TemperatureRangePresentationReport {
    /// Builds a pure presentation projection of one accepted `P,T` range.
    pub fn from_solution(solution: &TemperatureRangeSolution) -> Result<Self, ReactionExtentError> {
        if solution.points().is_empty() {
            return Err(invalid_range_presentation(
                "temperature range must contain at least one accepted point",
            ));
        }

        let mut ordered = solution.points().iter().enumerate().collect::<Vec<_>>();
        ordered.sort_by(|(_, left), (_, right)| {
            left.solution()
                .conditions()
                .temperature()
                .total_cmp(&right.solution().conditions().temperature())
        });

        let first_solution = ordered[0].1.solution();
        let expected_component_labels = component_labels(first_solution);
        let expected_phase_labels = phase_labels(first_solution);
        let mut point_rows = Vec::with_capacity(ordered.len());
        let mut component_rows = Vec::with_capacity(ordered.len());
        let mut component_fraction_rows = Vec::with_capacity(ordered.len());
        let mut phase_total_rows = Vec::with_capacity(ordered.len());
        let mut solver_metric_rows = Vec::with_capacity(ordered.len());
        let mut phase_statuses_by_point = Vec::with_capacity(ordered.len());

        for (point_index, point) in &ordered {
            let accepted = point.solution();
            if expected_component_labels != component_labels(accepted)
                || expected_phase_labels != phase_labels(accepted)
            {
                return Err(invalid_range_presentation(
                    "accepted range points do not share one canonical phase/component layout",
                ));
            }

            let report = point.report();
            let validation = accepted.accepted_solution().validation();
            let temperature = accepted.conditions().temperature();
            point_rows.push(TemperatureRangePresentationPointRow {
                point_index: *point_index,
                temperature_kelvin: temperature,
                used_continuation_seed: report.used_continuation_seed(),
                thermochemistry_refreshed: report.thermochemistry_refreshed(),
                symbolic_parameter_reused: report.symbolic_parameter_reused(),
                phase_control_transitions: report.phase_control_transitions(),
                phase_control_iterations: report.phase_control_iterations(),
                phase_set_reused: report.phase_set_reused(),
                formulation_build: report.formulation_build(),
                total_duration: report.timing().total(),
                accepted_backend: format!("{:?}", accepted.solve_report().accepted_backend),
                residual_l2_norm: validation.residual_l2_norm,
                max_abs_element_balance_error: validation.max_abs_element_balance_error,
            });
            component_rows.push((temperature, accepted.component_moles().to_vec()));
            component_fraction_rows.push((
                temperature,
                accepted
                    .metadata()
                    .components()
                    .iter()
                    .map(|component| accepted.mole_fraction_for(component.id()).unwrap_or(0.0))
                    .collect(),
            ));
            phase_total_rows.push((
                temperature,
                accepted
                    .phases()
                    .iter()
                    .map(|phase| accepted.phase_total(phase.id()).unwrap_or(0.0))
                    .collect(),
            ));
            solver_metric_rows.push((
                temperature,
                vec![
                    validation.residual_l2_norm,
                    validation.max_abs_element_balance_error,
                    report.timing().total().as_secs_f64() * 1_000.0,
                ],
            ));
            phase_statuses_by_point.push(phase_statuses(accepted));
        }

        let transition_boundaries = point_rows
            .windows(2)
            .zip(phase_statuses_by_point.windows(2))
            .filter_map(|(points, statuses)| {
                let upper = &points[1];
                let phase_status_changed = statuses[0] != statuses[1];
                (upper.phase_control_transitions > 0 || phase_status_changed).then(|| {
                    TemperatureRangeTransitionBoundary {
                        lower_point_index: points[0].point_index,
                        upper_point_index: upper.point_index,
                        lower_temperature_kelvin: points[0].temperature_kelvin,
                        upper_temperature_kelvin: upper.temperature_kelvin,
                        transitions_at_upper_point: upper.phase_control_transitions,
                        phase_status_changed,
                    }
                })
            })
            .collect();

        Ok(Self {
            points: point_rows,
            component_moles: TemperatureSweepSeries::from_rows(
                expected_component_labels.clone(),
                &component_rows,
            )?,
            component_mole_fractions: TemperatureSweepSeries::from_rows(
                expected_component_labels,
                &component_fraction_rows,
            )?,
            phase_totals: TemperatureSweepSeries::from_rows(
                expected_phase_labels,
                &phase_total_rows,
            )?,
            solver_metrics: TemperatureSweepSeries::from_rows(
                vec![
                    "residual_l2_norm".to_string(),
                    "max_abs_element_balance_error".to_string(),
                    "point_elapsed_ms".to_string(),
                ],
                &solver_metric_rows,
            )?,
            transition_boundaries,
        })
    }

    /// Returns `true` when the complete plotted range has no transition boundary.
    pub fn is_phase_stable(&self) -> bool {
        self.transition_boundaries.is_empty()
    }

    /// Postprocesses component moles without allowing a false smooth curve across transitions.
    pub fn postprocess_component_moles(
        &self,
        policy: &TemperaturePostprocessingPolicy,
    ) -> Result<TemperaturePostprocessingResult, ReactionExtentError> {
        postprocess_series_phase_safely(&self.component_moles, &self.transition_boundaries, policy)
    }

    /// Postprocesses phase totals without allowing a false smooth curve across transitions.
    pub fn postprocess_phase_totals(
        &self,
        policy: &TemperaturePostprocessingPolicy,
    ) -> Result<TemperaturePostprocessingResult, ReactionExtentError> {
        postprocess_series_phase_safely(&self.phase_totals, &self.transition_boundaries, policy)
    }
}

/// Extracts phase-qualified component labels in canonical solver order.
///
/// These labels form the column headers of a temperature-range series and must
/// stay stable and deterministic for plotting and export.
fn component_labels(solution: &MultiphaseEquilibriumSolution) -> Vec<String> {
    solution
        .metadata()
        .components()
        .iter()
        .map(|component| component.label())
        .collect()
}

/// Extracts `label::total_moles` series headers for every phase.
///
/// Uses the semantic phase name (falling back to `single`) so phase totals can
/// be projected as distinct labelled columns in a temperature-range series.
fn phase_labels(solution: &MultiphaseEquilibriumSolution) -> Vec<String> {
    solution
        .phases()
        .iter()
        .map(|phase| {
            let label = phase
                .id()
                .as_option()
                .clone()
                .unwrap_or_else(|| "single".to_string());
            format!("{label}::total_moles")
        })
        .collect()
}

/// Projects the lifecycle status of every phase into a stable string series.
///
/// Renders the accepted status (`Active`, `Inactive`, ...) per phase, using
/// `Unknown` when the phase is absent from the solution, so point diagnostics
/// reveal where phase transitions occurred across the sweep.
fn phase_statuses(solution: &MultiphaseEquilibriumSolution) -> Vec<String> {
    solution
        .phases()
        .iter()
        .map(|phase| {
            solution
                .phase_status(phase.id())
                .map(|status| format!("{status:?}"))
                .unwrap_or_else(|| "Unknown".to_string())
        })
        .collect()
}

/// Post-processes a temperature series, forbidding resampling across reported
/// phase-transition boundaries.
///
/// Interpolating over a phase transition would fabricate compositions that
/// never existed. When any transition boundary is present and the policy asks
/// for resampling (anything other than `RawOnly`), this helper returns a typed
/// error asking the caller to keep raw points or split the range into
/// phase-stable segments.
fn postprocess_series_phase_safely(
    series: &TemperatureSweepSeries,
    transition_boundaries: &[TemperatureRangeTransitionBoundary],
    policy: &TemperaturePostprocessingPolicy,
) -> Result<TemperaturePostprocessingResult, ReactionExtentError> {
    if !transition_boundaries.is_empty()
        && !matches!(policy.grid, TemperatureResamplingGrid::RawOnly)
    {
        return Err(invalid_range_presentation(
            "range resampling is forbidden across reported phase-transition boundaries; keep raw points or split the range into phase-stable segments",
        ));
    }
    Ok(TemperaturePostprocessingResult {
        raw: series.clone(),
        resampled: series.resample(policy)?,
    })
}

/// Constructs an invalid-problem error scoped to range presentation.
///
/// Centralizes the diagnostic field label so every presentation rejection
/// points at the same namespace and carries a consistent user-facing message.
fn invalid_range_presentation(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidProblem {
        field: "temperature_range_presentation",
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
        TemperatureInterpolationPolicy, TemperatureResamplingGrid,
    };
    use crate::Thermodynamics::ChemEquilibrium::prelude::{
        EquilibriumConditions, EquilibriumSolveOptions, LegacyEquilibriumSolver,
        PhaseEquilibriumPipelineRequest, SolverBackend, SolverPolicy, SubstanceSystemSpecBuilder,
        SubstancesContainer, TemperatureGrid,
    };

    fn series() -> TemperatureSweepSeries {
        TemperatureSweepSeries::from_rows(
            vec!["gas::A".to_string()],
            &[(300.0, vec![1.0]), (400.0, vec![2.0]), (500.0, vec![3.0])],
        )
        .unwrap()
    }

    fn boundary() -> TemperatureRangeTransitionBoundary {
        TemperatureRangeTransitionBoundary {
            lower_point_index: 0,
            upper_point_index: 1,
            lower_temperature_kelvin: 300.0,
            upper_temperature_kelvin: 400.0,
            transitions_at_upper_point: 1,
            phase_status_changed: true,
        }
    }

    #[test]
    fn phase_stable_series_can_resample_but_transition_series_cannot() {
        let policy = TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::Uniform { points: 5 },
            interpolation: TemperatureInterpolationPolicy::default(),
        };
        let stable = postprocess_series_phase_safely(&series(), &[], &policy).unwrap();
        assert_eq!(stable.resampled.unwrap().point_count(), 5);

        let error = postprocess_series_phase_safely(&series(), &[boundary()], &policy)
            .expect_err("a phase-transition boundary must forbid smoothing");
        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "temperature_range_presentation",
                ..
            }
        ));
    }

    #[test]
    fn transition_series_still_exposes_raw_points() {
        let raw_policy = TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::RawOnly,
            interpolation: TemperatureInterpolationPolicy::default(),
        };
        let result =
            postprocess_series_phase_safely(&series(), &[boundary()], &raw_policy).unwrap();
        assert_eq!(result.raw.point_count(), 3);
        assert!(result.resampled.is_none());
    }

    #[test]
    fn offline_local_range_projects_canonical_series_and_continuation() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("offline local-NASA specification must validate");
        let options = EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(
                LegacyEquilibriumSolver::NR,
            )))
            .expect("single legacy-NR policy must validate");
        let range = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(400.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_solve_options(options)
        .solve_temperature_range(TemperatureGrid::new(vec![400.0, 500.0, 600.0]).unwrap())
        .expect("offline local-NASA range must solve");

        let presentation = TemperatureRangePresentationReport::from_solution(&range)
            .expect("accepted range must produce a presentation report");
        assert_eq!(presentation.points.len(), 3);
        assert_eq!(presentation.component_moles.labels(), &["N2", "O2"]);
        assert_eq!(
            presentation.component_mole_fractions.labels(),
            &["N2", "O2"]
        );
        assert!(
            presentation
                .component_mole_fractions
                .rows()
                .iter()
                .all(|row| (row.iter().sum::<f64>() - 1.0).abs() < 1e-12)
        );
        assert_eq!(presentation.phase_totals.labels(), &["single::total_moles"]);
        assert_eq!(
            presentation.solver_metrics.labels(),
            &[
                "residual_l2_norm".to_string(),
                "max_abs_element_balance_error".to_string(),
                "point_elapsed_ms".to_string(),
            ]
        );
        assert_eq!(presentation.solver_metrics.rows().len(), 3);
        assert!(presentation.is_phase_stable());
        assert!(presentation.points[1].used_continuation_seed);
        assert_eq!(
            presentation.component_moles.temperatures(),
            &[400.0, 500.0, 600.0]
        );

        let postprocessed = presentation
            .postprocess_component_moles(&TemperaturePostprocessingPolicy {
                grid: TemperatureResamplingGrid::Uniform { points: 5 },
                interpolation: TemperatureInterpolationPolicy::default(),
            })
            .expect("phase-stable local range may be resampled");
        assert_eq!(
            postprocessed
                .resampled
                .expect("uniform grid requested")
                .point_count(),
            5
        );
    }
}
