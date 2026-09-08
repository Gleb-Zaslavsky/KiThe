//! Read-only presentation of fixed-pressure, target-enthalpy sweeps.
//!
//! P,H uses target enthalpy as its independent coordinate and solved
//! temperature as output. Its raw series therefore intentionally do not reuse
//! a temperature-grid type.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::PhRangeSolution;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use std::time::Duration;

/// Labelled values sampled on a strictly increasing target-enthalpy axis.
#[derive(Debug, Clone, PartialEq)]
pub struct EnthalpySweepSeries {
    labels: Vec<String>,
    target_enthalpies_joules: Vec<f64>,
    rows: Vec<Vec<f64>>,
}

impl EnthalpySweepSeries {
    /// Validates and stores raw range rows in increasing target-enthalpy order.
    pub fn from_rows(
        labels: Vec<String>,
        rows: &[(f64, Vec<f64>)],
    ) -> Result<Self, ReactionExtentError> {
        if labels.is_empty()
            || labels.iter().any(|label| label.trim().is_empty())
            || rows.is_empty()
        {
            return Err(invalid(
                "P,H presentation requires non-empty labels and rows",
            ));
        }
        let mut targets = Vec::with_capacity(rows.len());
        let mut previous = f64::NEG_INFINITY;
        for (index, (target, values)) in rows.iter().enumerate() {
            if !target.is_finite()
                || *target <= previous
                || values.len() != labels.len()
                || values.iter().any(|value| !value.is_finite())
            {
                return Err(invalid(format!("invalid P,H presentation row {index}")));
            }
            targets.push(*target);
            previous = *target;
        }
        Ok(Self {
            labels,
            target_enthalpies_joules: targets,
            rows: rows.iter().map(|(_, row)| row.clone()).collect(),
        })
    }

    /// Canonical component or phase labels.
    pub fn labels(&self) -> &[String] {
        &self.labels
    }
    /// Sorted target enthalpies in J.
    pub fn target_enthalpies_joules(&self) -> &[f64] {
        &self.target_enthalpies_joules
    }
    /// Row-major values aligned with labels and targets.
    pub fn rows(&self) -> &[Vec<f64>] {
        &self.rows
    }
}

/// One accepted P,H target and its route/validation evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct PhRangePresentationPointRow {
    /// Zero-based position in the caller-requested target-enthalpy grid.
    pub point_index: usize,
    /// Target total enthalpy for this point, in J.
    pub target_enthalpy_joules: f64,
    /// Temperature used to initialize this point's solve, in K.
    pub seed_temperature_kelvin: f64,
    /// Accepted equilibrium temperature for this point, in K.
    pub solved_temperature_kelvin: f64,
    /// Raw `H(T_solved) - H_target` at this point, in J.
    pub enthalpy_error_joules: f64,
    /// Dimensionless scale-aware enthalpy error used by acceptance.
    pub scaled_enthalpy_error: f64,
    /// Whether this point was warm-started from the previous accepted point.
    pub used_continuation: bool,
    /// Wall time spent solving this point.
    pub elapsed: Duration,
    /// Number of accepted phase-control transitions at this target.
    pub phase_control_transitions: usize,
    /// Number of prepared formulation builds attributed to this point.
    pub formulation_builds: usize,
    /// Number of prepared formulation reuses attributed to this point.
    pub formulation_reuses: usize,
    /// Numerical route that accepted this point (`Monolithic`, `Nested`, `Auto`).
    pub solve_path: String,
    /// Monolithic-to-nested fallback reason, when `Auto` recovered.
    pub fallback_reason: Option<String>,
    /// Backend that accepted the final candidate.
    pub accepted_backend: String,
    /// Residual L2 norm of the accepted chemical candidate.
    pub residual_l2_norm: f64,
    /// Largest absolute elemental-balance error of the accepted candidate.
    pub max_abs_element_balance_error: f64,
}

/// Raw `P,H` plotting and diagnostic data, sorted by target enthalpy.
#[derive(Debug, Clone, PartialEq)]
pub struct PhRangePresentationReport {
    /// Per-point diagnostics in the requested target-enthalpy order.
    pub points: Vec<PhRangePresentationPointRow>,
    /// Solved temperature series ordered by target enthalpy.
    pub solved_temperature: EnthalpySweepSeries,
    /// Phase-qualified component moles in canonical solver order.
    pub component_moles: EnthalpySweepSeries,
    /// Local phase mole fractions in the same component order as `component_moles`.
    pub component_mole_fractions: EnthalpySweepSeries,
    /// Accepted total mole amount per phase across the sweep.
    pub phase_totals: EnthalpySweepSeries,
    /// Residual, elemental-balance, and per-target elapsed-time evidence.
    /// Columns are `residual_l2_norm`, `max_abs_element_balance_error`, and
    /// `point_elapsed_ms` in that order.
    pub solver_metrics: EnthalpySweepSeries,
}

impl PhRangePresentationReport {
    /// Projects an accepted target-enthalpy range without reopening data or solving again.
    pub fn from_solution(solution: &PhRangeSolution) -> Result<Self, ReactionExtentError> {
        if solution.points().is_empty() {
            return Err(invalid("P,H range contains no accepted points"));
        }
        let mut ordered = solution.points().iter().enumerate().collect::<Vec<_>>();
        ordered.sort_by(|(_, left), (_, right)| {
            left.report()
                .target_enthalpy_joules()
                .total_cmp(&right.report().target_enthalpy_joules())
        });
        let first = ordered[0].1.solution().equilibrium();
        let expected_component_labels = component_labels(first);
        let expected_phase_labels = phase_labels(first);
        let mut points = Vec::with_capacity(ordered.len());
        let mut temperatures = Vec::with_capacity(ordered.len());
        let mut components = Vec::with_capacity(ordered.len());
        let mut component_fractions = Vec::with_capacity(ordered.len());
        let mut phases = Vec::with_capacity(ordered.len());
        let mut metrics = Vec::with_capacity(ordered.len());
        for (point_index, point) in ordered {
            let fixed = point.solution();
            let equilibrium = fixed.equilibrium();
            if expected_component_labels != component_labels(equilibrium)
                || expected_phase_labels != phase_labels(equilibrium)
            {
                return Err(invalid(
                    "P,H range points do not share one canonical phase/component layout",
                ));
            }
            let report = point.report();
            let validation = equilibrium.accepted_solution().validation();
            let target = report.target_enthalpy_joules();
            points.push(PhRangePresentationPointRow {
                point_index,
                target_enthalpy_joules: target,
                seed_temperature_kelvin: report.seed_temperature(),
                solved_temperature_kelvin: report.solved_temperature(),
                enthalpy_error_joules: fixed.enthalpy_error(),
                scaled_enthalpy_error: fixed.scaled_enthalpy_error(),
                used_continuation: report.used_continuation(),
                elapsed: report.elapsed(),
                phase_control_transitions: report.phase_control_transitions(),
                formulation_builds: report.formulation_builds(),
                formulation_reuses: report.formulation_reuses(),
                solve_path: format!("{:?}", report.solve_path()),
                fallback_reason: report
                    .fallback_reason()
                    .map(|reason| reason.message().to_string()),
                accepted_backend: format!("{:?}", equilibrium.solve_report().accepted_backend),
                residual_l2_norm: validation.residual_l2_norm,
                max_abs_element_balance_error: validation.max_abs_element_balance_error,
            });
            temperatures.push((target, vec![fixed.temperature()]));
            components.push((target, equilibrium.component_moles().to_vec()));
            component_fractions.push((
                target,
                equilibrium
                    .metadata()
                    .components()
                    .iter()
                    .map(|component| equilibrium.mole_fraction_for(component.id()).unwrap_or(0.0))
                    .collect(),
            ));
            phases.push((
                target,
                equilibrium
                    .phases()
                    .iter()
                    .map(|phase| equilibrium.phase_total(phase.id()).unwrap_or(0.0))
                    .collect(),
            ));
            metrics.push((
                target,
                vec![
                    validation.residual_l2_norm,
                    validation.max_abs_element_balance_error,
                    report.elapsed().as_secs_f64() * 1_000.0,
                ],
            ));
        }
        Ok(Self {
            points,
            solved_temperature: EnthalpySweepSeries::from_rows(
                vec!["solved_temperature_K".to_string()],
                &temperatures,
            )?,
            component_moles: EnthalpySweepSeries::from_rows(
                expected_component_labels.clone(),
                &components,
            )?,
            component_mole_fractions: EnthalpySweepSeries::from_rows(
                expected_component_labels,
                &component_fractions,
            )?,
            phase_totals: EnthalpySweepSeries::from_rows(expected_phase_labels, &phases)?,
            solver_metrics: EnthalpySweepSeries::from_rows(
                vec![
                    "residual_l2_norm".to_string(),
                    "max_abs_element_balance_error".to_string(),
                    "point_elapsed_ms".to_string(),
                ],
                &metrics,
            )?,
        })
    }
}

/// Extracts phase-qualified component labels in canonical solver order.
///
/// These labels form the column headers of a target-enthalpy series and must
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
/// be projected as distinct labelled columns in a `P,H` sweep series.
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

/// Constructs an invalid-problem error scoped to `P,H` range presentation.
///
/// Centralizes the diagnostic field label so every presentation rejection
/// points at the same namespace and carries a consistent user-facing message.
fn invalid(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidProblem {
        field: "ph_range_presentation",
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TemperatureBounds;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
        MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::{
        PhEnthalpyGrid, PhRangeRequest,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::ResolvedThermochemistry;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
        EquilibriumSolveOptions, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
    };
    use crate::Thermodynamics::User_PhaseOrSolution::{
        SubstanceSystemFactory, SubstanceSystemSpecBuilder, SubstancesContainer,
    };

    #[test]
    fn enthalpy_series_preserves_rows_and_rejects_non_monotone_targets() {
        let series = EnthalpySweepSeries::from_rows(
            vec!["T".to_string()],
            &[(-2.0, vec![300.0]), (3.0, vec![500.0])],
        )
        .unwrap();
        assert_eq!(series.target_enthalpies_joules(), &[-2.0, 3.0]);
        assert_eq!(series.rows(), &[vec![300.0], vec![500.0]]);
        assert!(
            EnthalpySweepSeries::from_rows(
                vec!["T".to_string()],
                &[(2.0, vec![300.0]), (2.0, vec![500.0])],
            )
            .is_err()
        );
    }

    #[test]
    fn local_nasa_ph_range_projects_raw_enthalpy_axis_and_continuation_evidence() {
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "H2".to_string(),
            "O2".to_string(),
            "H2O".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .unwrap();
        let resolved = SubstanceSystemFactory::resolve_phase_system(spec).unwrap();
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).unwrap();
        let initial =
            MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9]).unwrap();
        let pressure = 101_325.0;
        // Presentation must be independent of the numerical backend that the
        // canonical production cascade selects for a target.
        let options = EquilibriumSolveOptions::new().with_production_cascade();
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved).unwrap();
        let mut targets = [2_300.0, 2_500.0, 2_700.0]
            .into_iter()
            .map(|temperature| {
                let solution = solve_resolved_pt(
                    ResolvedPhaseEquilibriumRequest::new(
                        &resolved,
                        EquilibriumConditions::new(temperature, pressure, pressure).unwrap(),
                        initial.clone(),
                    )
                    .with_solve_options(options.clone()),
                )
                .unwrap();
                thermochemistry
                    .enthalpy_model()
                    .evaluate_total(solution.component_moles(), temperature)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        targets.sort_by(f64::total_cmp);
        let solution = PhRangeRequest::from_resolved_thermochemistry(
            &resolved,
            initial,
            pressure,
            pressure,
            PhEnthalpyGrid::new(targets.clone()).unwrap(),
            TemperatureBounds::new(2_100.0, 2_900.0).unwrap(),
            2_300.0,
            thermochemistry,
        )
        .unwrap()
        .with_solve_options(options)
        .solve()
        .unwrap();

        let presentation = PhRangePresentationReport::from_solution(&solution).unwrap();
        assert_eq!(presentation.points.len(), 3);
        assert_eq!(presentation.component_moles.labels().len(), 3);
        assert_eq!(presentation.phase_totals.labels(), &["single::total_moles"]);
        assert_eq!(presentation.points[0].point_index, 0);
        assert!(presentation.points[1].used_continuation);
        assert!(
            presentation
                .solved_temperature
                .target_enthalpies_joules()
                .windows(2)
                .all(|pair| pair[0] < pair[1])
        );
        assert_eq!(presentation.component_moles.rows().len(), 3);
        assert_eq!(presentation.component_mole_fractions.rows().len(), 3);
        assert!(
            presentation
                .component_mole_fractions
                .rows()
                .iter()
                .all(|row| (row.iter().sum::<f64>() - 1.0).abs() < 1e-10)
        );
        assert_eq!(presentation.phase_totals.rows().len(), 3);
        assert_eq!(presentation.solver_metrics.rows().len(), 3);
        assert_eq!(
            presentation.solver_metrics.labels(),
            &[
                "residual_l2_norm".to_string(),
                "max_abs_element_balance_error".to_string(),
                "point_elapsed_ms".to_string(),
            ]
        );
    }
}
