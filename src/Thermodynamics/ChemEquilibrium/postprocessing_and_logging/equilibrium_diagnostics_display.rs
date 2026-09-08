//! Human-readable rendering for structured equilibrium diagnostics.
//!
//! This module is deliberately outside the solver and phase-stability code.
//! It translates immutable typed evidence into text only after a caller has
//! explicitly enabled diagnostics and chosen to inspect or log the result.

use std::fmt::Write;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
    EquilibriumDiagnosticEvent, EquilibriumDiagnosticsReport, PhaseStabilityDiagnostic,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::FixedPressureEnthalpySolution;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseSet;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::ExtensiveNormalizationRecoveryEvidence;
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

/// Formats compact, always-available execution evidence for an accepted P,T
/// solution.
///
/// Unlike [`format_solution_diagnostics`], this summary does not require the
/// optional chronological event stream to have been enabled. It exposes the
/// immutable publication evidence that matters after a long run: accepted
/// backend and validation, committed phase-set transitions, extensive
/// normalization recovery, and opt-in timing. It never reconstructs a solve
/// history from mutable solver state.
pub fn format_solution_execution_summary(solution: &MultiphaseEquilibriumSolution) -> String {
    let validation = solution.accepted_solution().validation();
    let mut text = String::new();
    writeln!(
        text,
        "equilibrium execution: T={:.6} K, P={:.6e} Pa, backend={:?}, residual={:.3e}, balance={:.3e}",
        solution.conditions().temperature(),
        solution.conditions().pressure(),
        solution.solve_report().accepted_backend,
        validation.residual_l2_norm,
        validation.max_abs_element_balance_error,
    )
    .ok();
    format_phase_control_summary(&mut text, solution);
    if let Some(evidence) = solution.extensive_normalization_recovery() {
        writeln!(text, "{}", format_extensive_recovery_detail(evidence)).ok();
    }
    if let Some(report) = solution.diagnostics_report() {
        writeln!(
            text,
            "verbose diagnostics: retained_events={}, dropped_events={}",
            report.events().len(),
            report.dropped_events(),
        )
        .ok();
    } else {
        writeln!(text, "verbose diagnostics: disabled").ok();
    }
    format_timing_summary(&mut text, solution);
    text
}

/// Emits [`format_solution_execution_summary`] through the application's
/// configured [`log`] facade. This is an explicit caller action and remains
/// silent unless the caller invokes it.
pub fn log_solution_execution_summary(solution: &MultiphaseEquilibriumSolution) {
    log::info!("{}", format_solution_execution_summary(solution));
}

/// Formats the P,H route decision together with the accepted inner P,T
/// execution evidence.
///
/// A monolithic P,H solve has no scalar temperature trials, while a nested
/// solve does. The summary deliberately reports route-level counters instead
/// of fabricating trial rows, and makes an `Auto` fallback visible even when
/// detailed diagnostics were disabled.
pub fn format_ph_solution_execution_summary(solution: &FixedPressureEnthalpySolution) -> String {
    let report = solution.report();
    let mut text = String::new();
    writeln!(
        text,
        "P,H execution: route={:?}, T={:.6} K, P={:.6e} Pa, H_target={:.6e} J, H_calculated={:.6e} J, scaled_error={:.3e}, limit={:.3e} J",
        report.solve_path(),
        solution.temperature(),
        solution.pressure(),
        solution.target_enthalpy(),
        solution.calculated_enthalpy(),
        solution.scaled_enthalpy_error(),
        solution.enthalpy_error_limit_joules(),
    )
    .ok();
    if let Some(reason) = report.fallback_reason() {
        writeln!(text, "P,H Auto fallback accepted: {reason:?}").ok();
    }
    for decision in report.route_decisions() {
        writeln!(
            text,
            "P,H route decision: {:?} -> {:?}, reason={:?}",
            decision.from_route(),
            decision.to_route(),
            decision.reason(),
        )
        .ok();
    }
    writeln!(
        text,
        "P,H work: trials={}, inner_attempts={}, inner_iterations={}, phase_transitions={}, formulation_builds={}, formulation_reuses={}",
        report.trials().len(),
        report.inner_backend_attempts(),
        report.inner_nonlinear_iterations(),
        report.phase_control_transitions(),
        report.fixed_formulation_builds(),
        report.fixed_formulation_reuses(),
    )
    .ok();
    if let Some(monolithic) = report.monolithic_evidence() {
        writeln!(
            text,
            "P,H monolithic evidence: backend={:?}, attempts={}, residual_evaluations={}, jacobian_evaluations={}, multi_start={}, temperature_seed_recovery={}",
            monolithic.solve_report().accepted_backend,
            monolithic.solve_report().attempts.len(),
            monolithic.residual_evaluations(),
            monolithic.jacobian_evaluations(),
            monolithic.multi_start_report().is_some(),
            monolithic.temperature_seed_report().is_some(),
        )
        .ok();
    }
    let timing = report.timing();
    if timing.enabled() {
        writeln!(
            text,
            "P,H timing ms: total={:.3}, scalar_orchestration={:.3}, enthalpy_evaluation={:.3}",
            timing.total().as_secs_f64() * 1_000.0,
            timing.scalar_orchestration().as_secs_f64() * 1_000.0,
            timing.enthalpy_evaluation().as_secs_f64() * 1_000.0,
        )
        .ok();
    }
    text.push_str(&format_solution_execution_summary(solution.equilibrium()));
    text
}

/// Emits the route-aware P,H execution summary through [`log`].
pub fn log_ph_solution_execution_summary(solution: &FixedPressureEnthalpySolution) {
    log::info!("{}", format_ph_solution_execution_summary(solution));
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
            EquilibriumDiagnosticEvent::ExtensiveNormalizationRecoveryAccepted {
                physical_inventory_scale,
                trigger_kind,
                discovery_backend,
                physical_retry_backend,
                reconstructed_physical_boundary,
            } => {
                writeln!(
                    text,
                    "extensive normalization accepted: scale={physical_inventory_scale:.6e} mol, trigger={trigger_kind:?}, discovery={discovery_backend}, physical_retry={}, reconstructed={reconstructed_physical_boundary}",
                    physical_retry_backend.as_deref().unwrap_or("none"),
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
    if let Some(evidence) = solution.extensive_normalization_recovery() {
        writeln!(text, "{}", format_extensive_recovery_detail(evidence)).ok();
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

/// Explains the recovery route without presenting normalized coordinates as
/// physical result data. The immutable evidence retains the original failure,
/// so a successful recovery remains auditable rather than looking like a
/// backend silently succeeded on a different problem.
fn format_extensive_recovery_detail(evidence: &ExtensiveNormalizationRecoveryEvidence) -> String {
    format!(
        "recovery detail: physical solve failed numerically; normalized basin accepted by {:?}; publication={}; physical_failure={}; physical_retry_failure={}",
        evidence.discovery_backend,
        if evidence.reconstructed_physical_boundary {
            "reconstruction"
        } else {
            "physical retry"
        },
        evidence.trigger_message,
        evidence.physical_retry_failure.as_deref().unwrap_or("none"),
    )
}

/// Renders the committed active-set history retained by the accepted solution.
///
/// This is intentionally derived from `PhaseControlledSolveReport`, not from
/// the optional bounded diagnostics stream: a quiet production request still
/// deserves an auditable account of phase activation, deactivation, and the
/// accepted fixed-set candidate that caused each published transition.
fn format_phase_control_summary(text: &mut String, solution: &MultiphaseEquilibriumSolution) {
    let Some(report) = solution.phase_control_report() else {
        writeln!(text, "phase control: not used").ok();
        return;
    };
    writeln!(
        text,
        "phase control: iterations={}, transitions={}, initial=[{}], final=[{}], final_residual={:.3e}, final_balance={:.3e}",
        report.iterations,
        report.transitions.len(),
        format_phase_set(&report.initial_phase_set, solution),
        format_phase_set(&report.final_phase_set, solution),
        report.final_validation.residual_l2_norm,
        report.final_validation.max_abs_element_balance_error,
    )
    .ok();

    for transition in &report.transitions {
        let minimum_tpds = transition
            .minimum_tpds
            .iter()
            .enumerate()
            .filter_map(|(index, value)| {
                value.map(|tpd| format!("{}={tpd:.3e}", phase_label(index, solution)))
            })
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(
            text,
            "  transition {}: [{}] -> [{}], activated=[{}], deactivated=[{}], reason={:?}, backend={:?}, residual={:.3e}, balance={:.3e}, duration_ms={:.3}, TPD=[{}]",
            transition.iteration,
            format_phase_set(&transition.previous_phase_set, solution),
            format_phase_set(&transition.new_phase_set, solution),
            format_phase_index_ids(&transition.activated, solution),
            format_phase_index_ids(&transition.deactivated, solution),
            transition.reason,
            transition.nonlinear_report.accepted_backend,
            transition.candidate_validation.residual_l2_norm,
            transition.candidate_validation.max_abs_element_balance_error,
            transition.transition_duration.as_secs_f64() * 1_000.0,
            minimum_tpds,
        )
        .ok();
        if let (Some(&phase), Some(composition)) = (
            transition.activated.first(),
            transition.incipient_composition.as_deref(),
        ) {
            writeln!(
                text,
                "    activation seed: {}",
                format_composition(phase.into(), composition, solution),
            )
            .ok();
        }
    }
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

/// Appends one per-phase stability line to a human-readable diagnostics block.
///
/// Renders the phase label, active/inactive state, total mole amount, and the
/// minimum TPD value (or `not evaluated` when the stability analysis was not
/// run). When an incipient composition is available it is formatted on a
/// following indented line.
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

/// Renders the lifecycle status of every phase in a set as a compact string.
///
/// Produces a deterministic `label=status, ...` sequence in canonical phase
/// order so a diagnostics trace is reproducible and easy to scan.
fn format_phase_set(set: &PhaseSet, solution: &MultiphaseEquilibriumSolution) -> String {
    set.active_mask()
        .iter()
        .enumerate()
        .filter_map(|(index, active)| active.then(|| phase_label(index, solution)))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Joins a list of dense phase indices into a labelled, stable string.
///
/// Resolves each index to its semantic phase label and comma-joins the results.
/// Used to summarize which phases were activated or deactivated in a transition.
fn format_phase_indices(indices: &[usize], solution: &MultiphaseEquilibriumSolution) -> String {
    indices
        .iter()
        .map(|&index| phase_label(index, solution))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Converts the strongly typed phase indices retained by phase-control reports
/// only at the text-rendering boundary.
fn format_phase_index_ids(
    indices: &[crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex],
    solution: &MultiphaseEquilibriumSolution,
) -> String {
    indices
        .iter()
        .map(|&index| phase_label(index.into(), solution))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Resolves a dense phase index into its semantic display label.
///
/// Delegates to the accepted solution's phase metadata, falling back to a
/// numeric placeholder when the index is out of bounds so diagnostics rendering
/// never panics on malformed evidence.
fn phase_label(index: usize, solution: &MultiphaseEquilibriumSolution) -> String {
    solution
        .phases()
        .get(index)
        .and_then(|phase| phase.id().as_option().clone())
        .unwrap_or_else(|| format!("phase_{index}"))
}

/// Formats an incipient TPD composition as `label=x, ...` in component order.
///
/// Pairs each component label from the named phase with its mole fraction in
/// the minimizer composition, producing a deterministic human-readable row.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentErrorKind;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverBackend;

    #[test]
    fn recovery_detail_explains_failure_and_publication_route() {
        let evidence = ExtensiveNormalizationRecoveryEvidence {
            physical_inventory_scale: 1.0e4,
            trigger_kind: ReactionExtentErrorKind::AllBackendsFailed,
            trigger_message: "all physical backends failed".to_string(),
            discovery_backend: SolverBackend::Legacy(Solvers::NR),
            physical_retry_backend: None,
            physical_retry_failure: Some("physical retry remained ill-conditioned".to_string()),
            reconstructed_physical_boundary: true,
        };

        let rendered = format_extensive_recovery_detail(&evidence);
        assert!(rendered.contains("physical solve failed numerically"));
        assert!(rendered.contains("normalized basin accepted"));
        assert!(rendered.contains("publication=reconstruction"));
        assert!(rendered.contains("all physical backends failed"));
        assert!(rendered.contains("physical retry remained ill-conditioned"));
        assert!(!rendered.contains("normalized_moles"));
        assert!(!rendered.contains("normalized_total"));
    }
}
