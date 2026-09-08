//! Read-only comparisons between two accepted equilibrium solutions.
//!
//! The comparison boundary is intentionally strict: phase-qualified component
//! identities and semantic phase order must match exactly. This prevents the
//! common reporting bug of comparing two unrelated solver vectors merely
//! because they happen to have the same length.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;

/// Summary-level deltas between two accepted solutions with a shared layout.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumComparisonSummary {
    /// Shared process-independent layout fingerprint.
    pub layout_fingerprint: u64,
    /// Left-hand temperature in K.
    pub left_temperature_kelvin: f64,
    /// Right-hand temperature in K.
    pub right_temperature_kelvin: f64,
    /// Right minus left temperature in K.
    pub temperature_delta_kelvin: f64,
    /// Left-hand pressure in Pa.
    pub left_pressure_pa: f64,
    /// Right-hand pressure in Pa.
    pub right_pressure_pa: f64,
    /// Right minus left pressure in Pa.
    pub pressure_delta_pa: f64,
    /// Backend that accepted the left solution.
    pub left_backend: String,
    /// Backend that accepted the right solution.
    pub right_backend: String,
    /// Left accepted residual norm.
    pub left_residual_l2_norm: f64,
    /// Right accepted residual norm.
    pub right_residual_l2_norm: f64,
    /// Left maximum absolute elemental-balance error.
    pub left_max_abs_element_balance_error: f64,
    /// Right maximum absolute elemental-balance error.
    pub right_max_abs_element_balance_error: f64,
}

/// One phase-qualified component delta, in canonical component order.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumComponentComparisonRow {
    /// Collision-free `phase::substance` identity.
    pub component: String,
    /// Left published physical amount in mol.
    pub left_moles: f64,
    /// Right published physical amount in mol.
    pub right_moles: f64,
    /// Right minus left amount in mol.
    pub delta_moles: f64,
    /// Scale-aware absolute relative delta using `max(|left|, |right|, 1e-30)`.
    pub relative_delta: f64,
    /// Left selected source library.
    pub left_library: String,
    /// Right selected source library.
    pub right_library: String,
    /// Left selected source record key.
    pub left_record_key: String,
    /// Right selected source record key.
    pub right_record_key: String,
    /// Whether library or record selection differs between both solutions.
    pub provenance_changed: bool,
}

/// One semantic phase-total delta, in canonical phase order.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumPhaseComparisonRow {
    /// Semantic phase label; anonymous one-phase systems use `single`.
    pub phase: String,
    /// Left lifecycle state.
    pub left_status: String,
    /// Right lifecycle state.
    pub right_status: String,
    /// Left published total in mol.
    pub left_total_moles: f64,
    /// Right published total in mol.
    pub right_total_moles: f64,
    /// Right minus left total in mol.
    pub delta_moles: f64,
    /// Whether the explicit phase lifecycle state changed.
    pub status_changed: bool,
}

/// Full strict comparison report between two accepted equilibrium results.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumComparisonReport {
    /// Conditions, backend, and validation evidence.
    pub summary: EquilibriumComparisonSummary,
    /// Per-component physical/provenance deltas.
    pub components: Vec<EquilibriumComponentComparisonRow>,
    /// Per-phase physical/lifecycle deltas.
    pub phases: Vec<EquilibriumPhaseComparisonRow>,
}

impl EquilibriumComparisonReport {
    /// Compares two accepted solutions with exactly the same canonical layout.
    ///
    /// This method never aligns by vector position alone. A layout fingerprint
    /// mismatch or any phase/component identity mismatch is rejected before
    /// values or provenance are read.
    pub fn between(
        left: &MultiphaseEquilibriumSolution,
        right: &MultiphaseEquilibriumSolution,
    ) -> Result<Self, ReactionExtentError> {
        ensure_compatible_layouts(left, right)?;

        let left_validation = left.accepted_solution().validation();
        let right_validation = right.accepted_solution().validation();
        let left_conditions = left.conditions();
        let right_conditions = right.conditions();
        let summary = EquilibriumComparisonSummary {
            layout_fingerprint: left.metadata().layout_fingerprint(),
            left_temperature_kelvin: left_conditions.temperature(),
            right_temperature_kelvin: right_conditions.temperature(),
            temperature_delta_kelvin: right_conditions.temperature()
                - left_conditions.temperature(),
            left_pressure_pa: left_conditions.pressure(),
            right_pressure_pa: right_conditions.pressure(),
            pressure_delta_pa: right_conditions.pressure() - left_conditions.pressure(),
            left_backend: format!("{:?}", left.solve_report().accepted_backend),
            right_backend: format!("{:?}", right.solve_report().accepted_backend),
            left_residual_l2_norm: left_validation.residual_l2_norm,
            right_residual_l2_norm: right_validation.residual_l2_norm,
            left_max_abs_element_balance_error: left_validation.max_abs_element_balance_error,
            right_max_abs_element_balance_error: right_validation.max_abs_element_balance_error,
        };

        let components = left
            .metadata()
            .components()
            .iter()
            .zip(left.build_report().components())
            .zip(right.build_report().components())
            .enumerate()
            .map(|(index, ((component, left_source), right_source))| {
                let left_moles = left.component_moles()[index];
                let right_moles = right.component_moles()[index];
                let delta_moles = right_moles - left_moles;
                let left_lookup = left_source.thermo_source();
                let right_lookup = right_source.thermo_source();
                let left_library = left_lookup.library().to_string();
                let right_library = right_lookup.library().to_string();
                let left_record_key = left_lookup.record_key().to_string();
                let right_record_key = right_lookup.record_key().to_string();
                EquilibriumComponentComparisonRow {
                    component: component.label(),
                    left_moles,
                    right_moles,
                    delta_moles,
                    relative_delta: relative_delta(left_moles, right_moles),
                    provenance_changed: left_library != right_library
                        || left_record_key != right_record_key,
                    left_library,
                    right_library,
                    left_record_key,
                    right_record_key,
                }
            })
            .collect();

        let phases = left
            .phases()
            .iter()
            .zip(right.phases())
            .map(|(left_phase, right_phase)| {
                let left_status = left
                    .phase_status(left_phase.id())
                    .map(|status| format!("{status:?}"))
                    .unwrap_or_else(|| "Unknown".to_string());
                let right_status = right
                    .phase_status(right_phase.id())
                    .map(|status| format!("{status:?}"))
                    .unwrap_or_else(|| "Unknown".to_string());
                let left_total_moles = left.phase_total(left_phase.id()).unwrap_or(0.0);
                let right_total_moles = right.phase_total(right_phase.id()).unwrap_or(0.0);
                EquilibriumPhaseComparisonRow {
                    phase: phase_label(left_phase.id().as_option()),
                    status_changed: left_status != right_status,
                    left_status,
                    right_status,
                    left_total_moles,
                    right_total_moles,
                    delta_moles: right_total_moles - left_total_moles,
                }
            })
            .collect();

        Ok(Self {
            summary,
            components,
            phases,
        })
    }
}

/// Verifies that two accepted solutions share an identical phase layout.
///
/// Comparison is only meaningful when both solutions use the same
/// phase-qualified component identities in the same order. This guard rejects
/// a vector-by-position comparison of solutions that merely happen to have the
/// same length, returning an invalid-problem error naming the layout field.
fn ensure_compatible_layouts(
    left: &MultiphaseEquilibriumSolution,
    right: &MultiphaseEquilibriumSolution,
) -> Result<(), ReactionExtentError> {
    if left.metadata().layout_fingerprint() != right.metadata().layout_fingerprint()
        || left.metadata().components() != right.metadata().components()
        || left.phases() != right.phases()
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "equilibrium_comparison_layout",
            message: "accepted solutions use different phase-qualified layouts and cannot be compared by vector position".to_string(),
        });
    }
    Ok(())
}

/// Computes a scale-invariant relative difference between two values.
///
/// Normalizes the absolute difference by the larger magnitude, floored at
/// `1e-30` so a comparison between two zero values returns `0.0` instead of
/// dividing by zero. Used for mole-fraction and energy deltas where absolute
/// differences alone would be misleading across vastly different scales.
fn relative_delta(left: f64, right: f64) -> f64 {
    (right - left).abs() / left.abs().max(right.abs()).max(1e-30)
}

/// Renders an optional semantic phase name into a stable display label.
///
/// Returns the stored name when present, otherwise the placeholder
/// `"single"` for the common single-phase case. This keeps comparison rows and
/// report labels deterministic for unnamed phases.
fn phase_label(phase: &Option<String>) -> String {
    phase.clone().unwrap_or_else(|| "single".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::prelude::{
        EquilibriumConditions, EquilibriumSolveOptions, LegacyEquilibriumSolver,
        PhaseEquilibriumPipelineRequest, SolverBackend, SolverPolicy, SubstanceSystemSpecBuilder,
        SubstancesContainer,
    };

    fn local_solution(temperature: f64) -> MultiphaseEquilibriumSolution {
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
        PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
        )
        .with_solve_options(options)
        .solve()
        .expect("offline local-NASA point must solve")
        .into_solution()
    }

    #[test]
    fn comparison_aligns_local_solutions_by_qualified_layout_and_provenance() {
        let report =
            EquilibriumComparisonReport::between(&local_solution(400.0), &local_solution(600.0))
                .expect("same offline layout must compare");

        assert_eq!(report.summary.temperature_delta_kelvin, 200.0);
        assert_eq!(report.components.len(), 2);
        assert_eq!(report.components[0].component, "N2");
        assert!(
            report
                .components
                .iter()
                .all(|component| component.left_library == "NASA_gas")
        );
        assert!(
            report
                .components
                .iter()
                .all(|component| !component.provenance_changed)
        );
        assert_eq!(report.phases.len(), 1);
        assert_eq!(report.phases[0].phase, "single");
    }

    #[test]
    fn relative_delta_handles_trace_scale_without_division_by_zero() {
        assert_eq!(relative_delta(0.0, 0.0), 0.0);
        assert_eq!(relative_delta(0.0, 1e-31), 0.1);
        assert_eq!(relative_delta(2.0, 3.0), 1.0 / 3.0);
    }
}
