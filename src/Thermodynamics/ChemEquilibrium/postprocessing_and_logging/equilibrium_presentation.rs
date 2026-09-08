//! Read-only presentation views for accepted chemical-equilibrium solutions.
//!
//! The solver result already owns the canonical physical state, lookup
//! provenance, validation evidence, backend trace, and optional timing. This
//! module projects that immutable evidence into deterministic row-oriented
//! data suitable for a GUI, CLI table, CSV writer, or application log. It
//! deliberately does not own cache state and never re-evaluates thermochemistry.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverAttemptOutcome, SolverAttemptReport,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use std::fmt::Write;
use std::time::Duration;

/// Stable summary of one accepted equilibrium result.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumPresentationSummary {
    /// Solved temperature in K.
    pub temperature_kelvin: f64,
    /// Solved pressure in Pa.
    pub pressure_pa: f64,
    /// Stable fingerprint of the phase/component layout.
    pub layout_fingerprint: u64,
    /// Backend whose candidate passed common acceptance checks.
    pub accepted_backend: String,
    /// Number of backends that actually started.
    pub started_backend_attempts: usize,
    /// Total nonlinear iterations represented by the accepted solve.
    pub nonlinear_iterations: usize,
    /// Number of accepted phase lifecycle transitions.
    pub phase_transitions: usize,
    /// Norm of the residual used by the common acceptance gate.
    pub residual_l2_norm: f64,
    /// Maximum absolute elemental-balance error of the accepted physical state.
    pub max_abs_element_balance_error: f64,
}

/// One phase row in canonical semantic phase order.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumPhasePresentationRow {
    /// User-facing semantic phase name; anonymous one-phase systems use `single`.
    pub phase: String,
    /// Lifecycle state such as `Active` or `Inactive`.
    pub status: String,
    /// Published physical total in mol.
    pub physical_total_moles: f64,
    /// Positive numerical total retained for trace-floor diagnostics in mol.
    pub numerical_total_moles: f64,
}

/// One component row in exact solver-vector order.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumComponentPresentationRow {
    /// Collision-free `phase::substance` label.
    pub component: String,
    /// Semantic phase owning the component.
    pub phase: String,
    /// Bare canonical substance name.
    pub substance: String,
    /// Published physical amount in mol.
    pub physical_moles: f64,
    /// Positive numerical amount retained for diagnostics in mol.
    pub numerical_moles: f64,
    /// Local phase mole fraction.
    pub mole_fraction: f64,
    /// Physical initial amount before trace seeding in mol.
    pub initial_moles: f64,
    /// Standard-state Gibbs energy evaluated at the solved conditions in J/mol.
    pub standard_gibbs_j_per_mol: f64,
    /// Library that supplied the thermochemical record.
    pub library: String,
    /// Exact source-record key selected by lookup.
    pub record_key: String,
    /// Lookup tier that selected the record.
    pub lookup_priority: String,
}

/// One solver-cascade attempt, retained in execution order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumBackendAttemptPresentationRow {
    /// Zero-based position in the ordered cascade.
    pub attempt_index: usize,
    /// Backend label.
    pub backend: String,
    /// `accepted`, `failed`, `rejected_candidate`, or `skipped`.
    pub outcome: String,
    /// Optional typed backend failure category.
    pub failure_kind: Option<String>,
    /// Optional backend stopping condition.
    pub termination: Option<String>,
    /// Backend-reported convergence flag, when metrics exist.
    pub backend_converged: Option<bool>,
    /// Completed nonlinear iterations, when metrics exist.
    pub iterations: Option<usize>,
    /// Residual callback evaluations, when metrics exist.
    pub residual_evaluations: Option<usize>,
    /// Jacobian callback evaluations, when metrics exist.
    pub jacobian_evaluations: Option<usize>,
    /// Linear subproblems solved, when metrics exist.
    pub linear_solves: Option<usize>,
    /// Backend wall-clock time, when metrics exist.
    pub elapsed_millis: Option<u128>,
    /// Residual-callback time, when the backend provides a split.
    pub residual_evaluation_micros: Option<u128>,
    /// Jacobian-callback time, when the backend provides a split.
    pub jacobian_evaluation_micros: Option<u128>,
    /// Solver-internal time outside callbacks, when the backend provides a split.
    pub solver_overhead_micros: Option<u128>,
    /// Stable human-readable detail for a failed, rejected, or skipped attempt.
    pub detail: Option<String>,
}

/// One optional workflow timing stage.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumTimingPresentationRow {
    /// Stable stage key shared by terminal and GUI consumers.
    pub stage: &'static str,
    /// Measured duration for this stage.
    pub duration: Duration,
}

/// Complete read-only presentation projection of one accepted solution.
///
/// Its vectors preserve the solution's canonical ordering. Consumers can sort
/// a copy for display, but must retain this order when correlating rows with
/// solver vectors or phase-local composition.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumPresentationReport {
    /// High-level accepted-state and validation evidence.
    pub summary: EquilibriumPresentationSummary,
    /// Phase totals and lifecycle status in canonical phase order.
    pub phases: Vec<EquilibriumPhasePresentationRow>,
    /// Component amounts and lookup provenance in canonical component order.
    pub components: Vec<EquilibriumComponentPresentationRow>,
    /// Full solver-cascade trace, including failed and skipped backends.
    pub backend_attempts: Vec<EquilibriumBackendAttemptPresentationRow>,
    /// Optional workflow timing stages. Empty when timing was disabled.
    pub timing: Vec<EquilibriumTimingPresentationRow>,
}

impl EquilibriumPresentationReport {
    /// Builds a pure projection of an already accepted equilibrium solution.
    pub fn from_solution(solution: &MultiphaseEquilibriumSolution) -> Self {
        let conditions = solution.conditions();
        let validation = solution.accepted_solution().validation();
        let summary = EquilibriumPresentationSummary {
            temperature_kelvin: conditions.temperature(),
            pressure_pa: conditions.pressure(),
            layout_fingerprint: solution.metadata().layout_fingerprint(),
            accepted_backend: format!("{:?}", solution.solve_report().accepted_backend),
            started_backend_attempts: solution.started_backend_attempts(),
            nonlinear_iterations: solution.nonlinear_iterations(),
            phase_transitions: solution.phase_control_transitions(),
            residual_l2_norm: validation.residual_l2_norm,
            max_abs_element_balance_error: validation.max_abs_element_balance_error,
        };

        let phases = solution
            .phases()
            .iter()
            .map(|phase| EquilibriumPhasePresentationRow {
                phase: phase_label(phase.id().as_option()),
                status: solution
                    .phase_status(phase.id())
                    .map(|status| format!("{status:?}"))
                    .unwrap_or_else(|| "Unknown".to_string()),
                physical_total_moles: solution.phase_total(phase.id()).unwrap_or(0.0),
                numerical_total_moles: solution.numerical_phase_total(phase.id()).unwrap_or(0.0),
            })
            .collect();

        let components = solution
            .metadata()
            .components()
            .iter()
            .zip(solution.build_report().components())
            .enumerate()
            .map(|(index, (component, preparation))| {
                let source = preparation.thermo_source();
                EquilibriumComponentPresentationRow {
                    component: component.label(),
                    phase: phase_label(component.id().phase.as_option()),
                    substance: component.substance().to_string(),
                    physical_moles: solution.component_moles()[index],
                    numerical_moles: solution.numerical_component_moles()[index],
                    mole_fraction: solution.mole_fraction_for(component.id()).unwrap_or(0.0),
                    initial_moles: preparation.initial_moles(),
                    standard_gibbs_j_per_mol: preparation.standard_gibbs_at_conditions(),
                    library: source.library().to_string(),
                    record_key: source.record_key().to_string(),
                    lookup_priority: source.priority().to_string(),
                }
            })
            .collect();

        Self {
            summary,
            phases,
            components,
            backend_attempts: backend_attempt_rows(solution.solve_report()),
            timing: timing_rows(solution.timing_report().enabled(), solution),
        }
    }

    /// Renders a compact, ASCII-only diagnostic view for CLI logs and examples.
    ///
    /// Structured rows remain the canonical API. This renderer intentionally
    /// keeps only the fields most useful when a solve is inspected manually.
    pub fn render_compact(&self) -> String {
        let mut output = String::new();
        let _ = writeln!(
            output,
            "equilibrium: T={:.6} K, P={:.6e} Pa, backend={}, residual={:.3e}, balance={:.3e}",
            self.summary.temperature_kelvin,
            self.summary.pressure_pa,
            self.summary.accepted_backend,
            self.summary.residual_l2_norm,
            self.summary.max_abs_element_balance_error,
        );
        let _ = writeln!(
            output,
            "phases: phase | status | physical mol | numerical mol"
        );
        for phase in &self.phases {
            let _ = writeln!(
                output,
                "  {} | {} | {:.6e} | {:.6e}",
                phase.phase, phase.status, phase.physical_total_moles, phase.numerical_total_moles,
            );
        }
        let _ = writeln!(output, "components: component | mol | x | library | record");
        for component in &self.components {
            let _ = writeln!(
                output,
                "  {} | {:.6e} | {:.6e} | {} | {}",
                component.component,
                component.physical_moles,
                component.mole_fraction,
                component.library,
                component.record_key,
            );
        }
        let _ = writeln!(
            output,
            "backends: index | backend | outcome | iterations | elapsed ms"
        );
        for attempt in &self.backend_attempts {
            let _ = writeln!(
                output,
                "  {} | {} | {} | {} | {}",
                attempt.attempt_index,
                attempt.backend,
                attempt.outcome,
                optional_display(attempt.iterations),
                optional_display(attempt.elapsed_millis),
            );
        }
        if !self.timing.is_empty() {
            let _ = writeln!(output, "timing: stage | ms");
            for timing in &self.timing {
                let _ = writeln!(
                    output,
                    "  {} | {:.3}",
                    timing.stage,
                    timing.duration.as_secs_f64() * 1_000.0,
                );
            }
        }
        output
    }
}

/// Projects a complete backend cascade without requiring a full solution.
///
/// This is useful for typed failures where no accepted solution exists, and
/// makes failure diagnostics available to GUI worker code without fabricated
/// solution rows.
pub fn backend_attempt_rows(
    report: &EquilibriumSolveReport,
) -> Vec<EquilibriumBackendAttemptPresentationRow> {
    report
        .attempts
        .iter()
        .enumerate()
        .map(|(attempt_index, attempt)| backend_attempt_row(attempt_index, attempt))
        .collect()
}

/// Extracts backend-attempt diagnostics from a typed failed solve when present.
///
/// Nested P,H seed recovery, automatic route fallback, and range/trial
/// wrappers are traversed recursively. The returned rows therefore retain
/// the backend order even when no accepted solution exists at the top level.
pub fn backend_attempt_rows_from_error(
    error: &ReactionExtentError,
) -> Option<Vec<EquilibriumBackendAttemptPresentationRow>> {
    let mut attempts = Vec::new();
    collect_backend_attempts(error, &mut attempts);
    (!attempts.is_empty()).then(|| {
        attempts
            .into_iter()
            .enumerate()
            .map(|(attempt_index, attempt)| backend_attempt_row(attempt_index, attempt))
            .collect()
    })
}

fn collect_backend_attempts<'a>(
    error: &'a ReactionExtentError,
    output: &mut Vec<&'a SolverAttemptReport>,
) {
    match error {
        ReactionExtentError::AllBackendsFailed { attempts }
        | ReactionExtentError::CascadeAborted { attempts, .. } => output.extend(attempts),
        ReactionExtentError::PhMonolithicSeedRecoveryFailed { attempts } => {
            for attempt in attempts {
                collect_backend_attempts(attempt.cause(), output);
            }
        }
        ReactionExtentError::PhAutoFallbackFailed { monolithic, nested } => {
            collect_backend_attempts(monolithic, output);
            collect_backend_attempts(nested, output);
        }
        ReactionExtentError::TemperatureTrialFailed { cause, .. }
        | ReactionExtentError::TemperatureRangePointFailed { cause, .. } => {
            collect_backend_attempts(cause, output);
        }
        _ => {}
    }
}

fn backend_attempt_row(
    attempt_index: usize,
    attempt: &SolverAttemptReport,
) -> EquilibriumBackendAttemptPresentationRow {
    let (outcome, detail) = match &attempt.outcome {
        SolverAttemptOutcome::Accepted => ("accepted".to_string(), None),
        SolverAttemptOutcome::Failed { reason, .. } => ("failed".to_string(), Some(reason.clone())),
        SolverAttemptOutcome::RejectedCandidate { reason } => {
            ("rejected_candidate".to_string(), Some(reason.clone()))
        }
        SolverAttemptOutcome::Skipped { reason } => ("skipped".to_string(), Some(reason.clone())),
    };
    let metrics = attempt.metrics.as_ref();
    let timing = metrics.and_then(|metrics| metrics.evaluation_timing.as_ref());
    EquilibriumBackendAttemptPresentationRow {
        attempt_index,
        backend: format!("{:?}", attempt.backend),
        outcome,
        failure_kind: attempt.failure_kind().map(|kind| format!("{kind:?}")),
        termination: metrics.map(|metrics| format!("{:?}", metrics.termination)),
        backend_converged: metrics.map(|metrics| metrics.backend_converged),
        iterations: metrics.map(|metrics| metrics.iterations),
        residual_evaluations: metrics.map(|metrics| metrics.residual_evaluations),
        jacobian_evaluations: metrics.map(|metrics| metrics.jacobian_evaluations),
        linear_solves: metrics.map(|metrics| metrics.linear_solves),
        elapsed_millis: metrics.map(|metrics| metrics.elapsed_millis),
        residual_evaluation_micros: timing.map(|timing| timing.residual_evaluation_micros),
        jacobian_evaluation_micros: timing.map(|timing| timing.jacobian_evaluation_micros),
        solver_overhead_micros: timing.map(|timing| timing.solver_overhead_micros),
        detail,
    }
}

/// Projects stage timings into presentation rows when timing is enabled.
///
/// Returns an empty vector when timing was disabled so the presentation never
/// shows fabricated zero-duration stages. Otherwise it maps every measured
/// stage of the accepted solution's timing report to a labelled row.
fn timing_rows(
    enabled: bool,
    solution: &MultiphaseEquilibriumSolution,
) -> Vec<EquilibriumTimingPresentationRow> {
    if !enabled {
        return Vec::new();
    }
    let timing = solution.timing_report();
    [
        ("total", timing.total()),
        ("repository_lookup", timing.repository_lookup()),
        (
            "thermochemistry_preparation",
            timing.thermochemistry_preparation(),
        ),
        (
            "numeric_closure_construction",
            timing.numeric_closure_construction(),
        ),
        ("symbolic_construction", timing.symbolic_construction()),
        ("equation_construction", timing.equation_construction()),
        (
            "numerical_problem_preparation",
            timing.numerical_problem_preparation(),
        ),
        ("projection_build", timing.projection_build()),
        ("nonlinear_solve", timing.nonlinear_solve()),
        ("phase_control", timing.phase_control()),
        ("validation", timing.validation()),
        ("postprocessing", timing.postprocessing()),
    ]
    .into_iter()
    .map(|(stage, duration)| EquilibriumTimingPresentationRow { stage, duration })
    .collect()
}

/// Renders an optional semantic phase name into a stable display label.
///
/// Returns the stored name when present, otherwise the placeholder `"single"`
/// for the common unnamed single-phase case.
fn phase_label(phase: &Option<String>) -> String {
    phase.clone().unwrap_or_else(|| "single".to_string())
}

/// Renders an optional value as its display text, or `-` when absent.
///
/// Used to keep presentation rows column-aligned when a solver metric or
/// evidence field was not produced.
fn optional_display<T: std::fmt::Display>(value: Option<T>) -> String {
    value.map_or_else(|| "-".to_string(), |value| value.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::PhMonolithicSeedFailure;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        SolverAttemptFailureKind, SolverAttemptMetrics, SolverBackend, SolverEvaluationTiming,
        SolverPolicy, SolverTermination,
    };
    use crate::Thermodynamics::ChemEquilibrium::prelude::{
        EquilibriumConditions, EquilibriumSolveOptions, LegacyEquilibriumSolver,
        PhaseEquilibriumPipelineRequest, SubstanceSystemSpecBuilder, SubstancesContainer,
    };

    #[test]
    fn backend_rows_preserve_outcomes_metrics_and_diagnostics() {
        let legacy_lm = SolverBackend::Legacy(Solvers::LM);
        let legacy_nr = SolverBackend::Legacy(Solvers::NR);
        let legacy_tr = SolverBackend::Legacy(Solvers::TR);
        let report = EquilibriumSolveReport {
            policy: SolverPolicy::Cascade(vec![legacy_lm, legacy_nr, legacy_tr]),
            attempts: vec![
                SolverAttemptReport {
                    backend: legacy_lm,
                    outcome: SolverAttemptOutcome::Failed {
                        kind: SolverAttemptFailureKind::Backend,
                        reason: "step rejected".to_string(),
                    },
                    metrics: Some(SolverAttemptMetrics {
                        termination: SolverTermination::Stagnation,
                        backend_converged: false,
                        iterations: 7,
                        residual_evaluations: 8,
                        jacobian_evaluations: 7,
                        linear_solves: 7,
                        elapsed_millis: 12,
                        evaluation_timing: Some(SolverEvaluationTiming {
                            residual_evaluation_micros: 20,
                            jacobian_evaluation_micros: 30,
                            solver_overhead_micros: 40,
                        }),
                    }),
                },
                SolverAttemptReport {
                    backend: legacy_nr,
                    outcome: SolverAttemptOutcome::Accepted,
                    metrics: None,
                },
                SolverAttemptReport {
                    backend: legacy_tr,
                    outcome: SolverAttemptOutcome::Skipped {
                        reason: "cascade budget exhausted".to_string(),
                    },
                    metrics: None,
                },
            ],
            accepted_backend: legacy_nr,
        };

        let rows = backend_attempt_rows(&report);
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].outcome, "failed");
        assert_eq!(rows[0].failure_kind.as_deref(), Some("Backend"));
        assert_eq!(rows[0].iterations, Some(7));
        assert_eq!(rows[0].jacobian_evaluation_micros, Some(30));
        assert_eq!(rows[1].outcome, "accepted");
        assert_eq!(rows[1].detail, None);
        assert_eq!(rows[2].outcome, "skipped");
        assert_eq!(rows[2].detail.as_deref(), Some("cascade budget exhausted"));
    }

    #[test]
    fn failed_cascade_projects_the_same_attempt_trace_without_a_solution() {
        let error = ReactionExtentError::AllBackendsFailed {
            attempts: vec![SolverAttemptReport {
                backend: SolverBackend::Legacy(Solvers::TR),
                outcome: SolverAttemptOutcome::RejectedCandidate {
                    reason: "physical balances did not pass acceptance".to_string(),
                },
                metrics: None,
            }],
        };

        let rows = backend_attempt_rows_from_error(&error)
            .expect("all-backends failure must retain attempt rows");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].outcome, "rejected_candidate");
        assert_eq!(
            rows[0].detail.as_deref(),
            Some("physical balances did not pass acceptance")
        );
        assert!(backend_attempt_rows_from_error(&ReactionExtentError::Cancelled).is_none());
    }

    #[test]
    fn monolithic_seed_failure_flattens_every_preserved_backend_trace() {
        let failure = |solver| ReactionExtentError::AllBackendsFailed {
            attempts: vec![SolverAttemptReport {
                backend: SolverBackend::Legacy(solver),
                outcome: SolverAttemptOutcome::Failed {
                    kind: SolverAttemptFailureKind::Solver,
                    reason: "iteration limit".to_string(),
                },
                metrics: None,
            }],
        };
        let error = ReactionExtentError::PhMonolithicSeedRecoveryFailed {
            attempts: vec![
                PhMonolithicSeedFailure::new(600.0, failure(Solvers::LM)),
                PhMonolithicSeedFailure::new(1_000.0, failure(Solvers::NR)),
            ],
        };

        let rows = backend_attempt_rows_from_error(&error)
            .expect("seed recovery must expose its preserved backend traces");
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].attempt_index, 0);
        assert_eq!(rows[0].backend, "Legacy(LM)");
        assert_eq!(rows[1].attempt_index, 1);
        assert_eq!(rows[1].backend, "Legacy(NR)");
    }

    #[test]
    fn compact_renderer_keeps_missing_backend_metrics_explicit() {
        let report = EquilibriumPresentationReport {
            summary: EquilibriumPresentationSummary {
                temperature_kelvin: 900.0,
                pressure_pa: 101_325.0,
                layout_fingerprint: 7,
                accepted_backend: "Legacy(NR)".to_string(),
                started_backend_attempts: 1,
                nonlinear_iterations: 4,
                phase_transitions: 0,
                residual_l2_norm: 1.0e-9,
                max_abs_element_balance_error: 2.0e-12,
            },
            phases: Vec::new(),
            components: Vec::new(),
            backend_attempts: vec![EquilibriumBackendAttemptPresentationRow {
                attempt_index: 0,
                backend: "Legacy(NR)".to_string(),
                outcome: "accepted".to_string(),
                failure_kind: None,
                termination: None,
                backend_converged: None,
                iterations: None,
                residual_evaluations: None,
                jacobian_evaluations: None,
                linear_solves: None,
                elapsed_millis: None,
                residual_evaluation_micros: None,
                jacobian_evaluation_micros: None,
                solver_overhead_micros: None,
                detail: None,
            }],
            timing: Vec::new(),
        };

        let rendered = report.render_compact();
        assert!(rendered.contains("equilibrium: T=900.000000 K"));
        assert!(rendered.contains("0 | Legacy(NR) | accepted | - | -"));
    }

    #[test]
    fn accepted_local_solution_projects_physical_rows_and_lookup_provenance() {
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
        let outcome = PhaseEquilibriumPipelineRequest::new(
            spec,
            vec![0.79, 0.21],
            EquilibriumConditions::new(500.0, 101_325.0, 101_325.0).unwrap(),
        )
        .with_solve_options(options)
        .solve()
        .expect("offline local-NASA point must solve");

        let presentation = EquilibriumPresentationReport::from_solution(outcome.solution());
        assert_eq!(presentation.phases.len(), 1);
        assert_eq!(presentation.phases[0].phase, "single");
        assert_eq!(presentation.components.len(), 2);
        assert_eq!(presentation.components[0].component, "N2");
        assert_eq!(presentation.components[1].component, "O2");
        assert!(
            presentation
                .components
                .iter()
                .all(|component| component.library == "NASA_gas")
        );
        assert!(
            presentation
                .components
                .iter()
                .all(|component| component.physical_moles.is_finite())
        );
        assert!(
            presentation
                .render_compact()
                .contains("components: component | mol | x | library | record")
        );
    }
}
