//! Human-readable rendering for structured equilibrium diagnostics.
//!
//! This module is deliberately outside the solver and phase-stability code.
//! It translates immutable typed evidence into text only after a caller has
//! explicitly enabled diagnostics and chosen to inspect or log the result.

use std::fmt::Write;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
    EquilibriumDiagnosticEvent, EquilibriumDiagnosticsReport, PhaseStabilityDiagnostic,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseSet;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;

/// Formats the diagnostic evidence attached to an accepted public solution.
///
/// Returns `None` when the request left diagnostics disabled. The phase-aware
/// solution supplies stable names for otherwise compact phase-index events.
pub fn format_solution_diagnostics(solution: &MultiphaseEquilibriumSolution) -> Option<String> {
    solution
        .diagnostics_report()
        .map(|report| format_diagnostics(report, solution))
}

/// Emits the formatted diagnostic tree through the application's configured
/// [`log`] facade. This is an explicit presentation action, never an engine
/// side effect.
pub fn log_solution_diagnostics(solution: &MultiphaseEquilibriumSolution) {
    if let Some(text) = format_solution_diagnostics(solution) {
        log::info!("{text}");
    }
}

/// Formats one immutable diagnostics report against its phase-qualified
/// solution metadata.
pub fn format_diagnostics(
    report: &EquilibriumDiagnosticsReport,
    solution: &MultiphaseEquilibriumSolution,
) -> String {
    let mut text = String::new();
    writeln!(text, "equilibrium diagnostics: mode={:?}", report.mode()).ok();
    for event in report.events() {
        match event {
            EquilibriumDiagnosticEvent::SolveStarted {
                conditions,
                initial_phase_set,
            } => {
                writeln!(
                    text,
                    "solve started: T={:.6} K, P={:.6} Pa, initial=[{}]",
                    conditions.temperature(),
                    conditions.pressure(),
                    format_phase_set(initial_phase_set, solution),
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::OuterIterationStarted {
                iteration,
                active_phase_set,
            } => {
                writeln!(
                    text,
                    "  outer iteration {iteration}: active=[{}]",
                    format_phase_set(active_phase_set, solution),
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::ActiveSetCandidateAccepted {
                iteration,
                validation,
                backend_summary,
                ..
            } => {
                writeln!(
                    text,
                    "    candidate accepted at iteration {iteration}: residual={:.3e}, balance={:.3e}, {backend_summary}",
                    validation.residual_l2_norm,
                    validation.max_abs_element_balance_error,
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::ActiveSetCandidateRejected {
                iteration,
                active_phase_set,
                message,
            } => {
                writeln!(
                    text,
                    "    candidate rejected at iteration {iteration}: active=[{}], reason={message}",
                    format_phase_set(active_phase_set, solution),
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::StabilityEvaluated {
                iteration,
                dg_create_j_per_mol,
                dg_keep_j_per_mol,
                phases,
            } => {
                writeln!(
                    text,
                    "    stability iteration {iteration}: dg_create={dg_create_j_per_mol:.3e} J/mol, dg_keep={dg_keep_j_per_mol:.3e} J/mol",
                )
                .ok();
                for phase in phases {
                    format_stability_line(&mut text, phase, solution);
                }
            }
            EquilibriumDiagnosticEvent::TransitionAccepted {
                iteration,
                activated_phase_indices,
                deactivated_phase_indices,
                reason,
                previous_phase_set,
                new_phase_set,
                incipient_composition,
            } => {
                writeln!(
                    text,
                    "    transition accepted at iteration {iteration}: [{}] -> [{}], activated=[{}], deactivated=[{}], reason={reason:?}",
                    format_phase_set(previous_phase_set, solution),
                    format_phase_set(new_phase_set, solution),
                    format_phase_indices(activated_phase_indices, solution),
                    format_phase_indices(deactivated_phase_indices, solution),
                )
                .ok();
                if let (Some(&phase), Some(composition)) = (
                    activated_phase_indices.first(),
                    incipient_composition.as_deref(),
                ) {
                    writeln!(
                        text,
                        "      restart seed from TPD minimizer: {}",
                        format_composition(phase, composition, solution),
                    )
                    .ok();
                }
            }
            EquilibriumDiagnosticEvent::TransitionHeldByHysteresis {
                iteration,
                phase_index,
            } => {
                writeln!(
                    text,
                    "    hysteresis retained phase {} at iteration {iteration}",
                    phase_label(*phase_index, solution),
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::RecoveryProbeAccepted {
                iteration,
                phase_index,
                previous_phase_set,
                new_phase_set,
            } => {
                writeln!(
                    text,
                    "    boundary recovery accepted at iteration {iteration}: removed {}, [{}] -> [{}]",
                    phase_label(*phase_index, solution),
                    format_phase_set(previous_phase_set, solution),
                    format_phase_set(new_phase_set, solution),
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::RecoveryProbeStarted {
                iteration,
                removed_phase_index,
                attempted_phase_set,
            } => {
                writeln!(
                    text,
                    "    boundary recovery probe at iteration {iteration}: remove {}, active=[{}]",
                    phase_label(*removed_phase_index, solution),
                    format_phase_set(attempted_phase_set, solution),
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::RecoveryProbeRejected {
                iteration,
                removed_phase_index,
                message,
            } => {
                writeln!(
                    text,
                    "    boundary recovery rejected at iteration {iteration}: remove {}, reason={message}",
                    phase_label(*removed_phase_index, solution),
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::ContinuationRestored {
                retained_phase_set,
                retained_seed,
            } => {
                let phases = retained_phase_set
                    .as_ref()
                    .map(|set| format_phase_set(set, solution))
                    .unwrap_or_else(|| "none".to_string());
                writeln!(
                    text,
                    "continuation restored: seed_retained={retained_seed}, phase_set=[{phases}]",
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::PhaseControlBudgetExhausted {
                max_outer_iterations,
            } => {
                writeln!(
                    text,
                    "phase-control budget exhausted after {max_outer_iterations} outer iterations",
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::PhaseControlCycleDetected {
                iteration,
                repeated_phase_set,
            } => {
                writeln!(
                    text,
                    "phase-control cycle detected at iteration {iteration}: repeated active=[{}]",
                    format_phase_set(repeated_phase_set, solution),
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::PhRouteFallback {
                from_route,
                to_route,
                error_kind,
                message,
            } => {
                writeln!(
                    text,
                    "P,H route fallback: {from_route:?} -> {to_route:?}, error_kind={error_kind:?}, reason={message}",
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::PhRouteFailed {
                route,
                error_kind,
                message,
            } => {
                writeln!(
                    text,
                    "P,H route failed: {route:?}, error_kind={error_kind:?}, reason={message}",
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::SolveFailed {
                message,
                continuation_restored,
            } => {
                writeln!(
                    text,
                    "solve failed: continuation_restored={continuation_restored}, reason={message}",
                )
                .ok();
            }
            EquilibriumDiagnosticEvent::SolveAccepted {
                final_phase_set,
                outer_iterations,
                transition_count,
            } => {
                writeln!(
                    text,
                    "solve accepted: active=[{}], outer_iterations={outer_iterations}, transitions={transition_count}",
                    format_phase_set(final_phase_set, solution),
                )
                .ok();
            }
        }
    }
    if report.dropped_events() > 0 {
        writeln!(
            text,
            "diagnostics truncated: {} retained events omitted",
            report.dropped_events()
        )
        .ok();
    }
    format_timing_summary(&mut text, solution);
    text
}

/// Appends already-measured stage timing without making diagnostics itself a
/// second timer. Disabled timing deliberately stays absent from the trace.
fn format_timing_summary(text: &mut String, solution: &MultiphaseEquilibriumSolution) {
    let timing = solution.timing_report();
    if !timing.enabled() {
        return;
    }
    let milliseconds = |duration: std::time::Duration| duration.as_secs_f64() * 1_000.0;
    writeln!(
        text,
        "timing ms: total={:.3}, lookup={:.3}, thermo={:.3}, closures={:.3}, symbolic={:.3}, projection={:.3}, nonlinear={:.3}, phase_control={:.3}, validation={:.3}, post={:.3}",
        milliseconds(timing.total()),
        milliseconds(timing.repository_lookup()),
        milliseconds(timing.thermochemistry_preparation()),
        milliseconds(timing.numeric_closure_construction()),
        milliseconds(timing.symbolic_construction()),
        milliseconds(timing.projection_build()),
        milliseconds(timing.nonlinear_solve()),
        milliseconds(timing.phase_control()),
        milliseconds(timing.validation()),
        milliseconds(timing.postprocessing()),
    )
    .ok();
}

fn format_stability_line(
    text: &mut String,
    phase: &PhaseStabilityDiagnostic,
    solution: &MultiphaseEquilibriumSolution,
) {
    let tpd = phase
        .minimum_tpd_j_per_mol
        .map(|value| format!("{value:.3e} J/mol"))
        .unwrap_or_else(|| "not evaluated".to_string());
    writeln!(
        text,
        "      {}: {}, total={:.3e} mol, minimum_tpd={tpd}",
        phase_label(phase.phase_index, solution),
        if phase.active { "active" } else { "inactive" },
        phase.phase_total_moles,
    )
    .ok();
    if let Some(composition) = phase.incipient_composition.as_deref() {
        writeln!(
            text,
            "        TPD minimizer: {}",
            format_composition(phase.phase_index, composition, solution),
        )
        .ok();
    }
}

fn format_phase_set(set: &PhaseSet, solution: &MultiphaseEquilibriumSolution) -> String {
    set.active_mask()
        .iter()
        .enumerate()
        .filter_map(|(index, active)| active.then(|| phase_label(index, solution)))
        .collect::<Vec<_>>()
        .join(", ")
}

fn format_phase_indices(indices: &[usize], solution: &MultiphaseEquilibriumSolution) -> String {
    indices
        .iter()
        .map(|&index| phase_label(index, solution))
        .collect::<Vec<_>>()
        .join(", ")
}

fn phase_label(index: usize, solution: &MultiphaseEquilibriumSolution) -> String {
    solution
        .phases()
        .get(index)
        .and_then(|phase| phase.id().as_option().clone())
        .unwrap_or_else(|| format!("phase_{index}"))
}

fn format_composition(
    phase_index: usize,
    composition: &[f64],
    solution: &MultiphaseEquilibriumSolution,
) -> String {
    let Some(phase) = solution.phases().get(phase_index) else {
        return format!("{:?}", composition);
    };
    solution.metadata().components()[phase.component_range()]
        .iter()
        .zip(composition)
        .map(|(component, &fraction)| format!("{}={fraction:.4}", component.label()))
        .collect::<Vec<_>>()
        .join(", ")
}
