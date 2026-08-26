//! Immutable UI projection of accepted chemical-equilibrium solutions.
//!
//! The solver result already owns authoritative reports and provenance. This
//! module adds one compact, phase-qualified, aligned projection for tables and
//! plotting while retaining the source solutions through `Arc`.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_diagnostics::{
    EquilibriumDiagnosticEvent, PhaseStabilityDiagnostic,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
    FixedPressureEnthalpySolution, PhRouteDecision,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_range::TemperatureRangeSolveReport;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::gui::equilibrium_gui_request::EquilibriumGuiSolveOutcome;
use std::sync::Arc;

/// One immutable point in the UI result model.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumGuiPointSnapshot {
    temperature_k: f64,
    component_moles: Vec<f64>,
    mole_fractions: Vec<f64>,
    phase_totals: Vec<f64>,
    phase_statuses: Vec<String>,
    component_active: Vec<bool>,
    phase_active: Vec<bool>,
    lifecycle_trace: Option<EquilibriumGuiLifecycleTraceSnapshot>,
    source: Arc<MultiphaseEquilibriumSolution>,
}

/// One readable fact in the accepted phase-control decision tree.
///
/// The source event remains owned by the engine report. This small immutable
/// projection resolves phase indices into the labels already displayed by the
/// result table, so the GUI never has to reproduce lifecycle semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumGuiLifecycleEventSnapshot {
    title: String,
    detail: String,
}

impl EquilibriumGuiLifecycleEventSnapshot {
    pub fn title(&self) -> &str {
        &self.title
    }

    pub fn detail(&self) -> &str {
        &self.detail
    }
}

/// Bounded diagnostics attached to one accepted point.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumGuiLifecycleTraceSnapshot {
    mode: String,
    events: Vec<EquilibriumGuiLifecycleEventSnapshot>,
    dropped_events: usize,
}

impl EquilibriumGuiLifecycleTraceSnapshot {
    fn from_solution(solution: &MultiphaseEquilibriumSolution) -> Option<Self> {
        let report = solution.diagnostics_report()?;
        if !report.enabled() {
            return None;
        }
        Some(Self {
            mode: format!("{:?}", report.mode()),
            events: report
                .events()
                .iter()
                .map(|event| lifecycle_event_snapshot(event, solution))
                .collect(),
            dropped_events: report.dropped_events(),
        })
    }

    pub fn mode(&self) -> &str {
        &self.mode
    }

    pub fn events(&self) -> &[EquilibriumGuiLifecycleEventSnapshot] {
        &self.events
    }

    pub fn dropped_events(&self) -> usize {
        self.dropped_events
    }
}

/// Scalar energy evidence attached to a fixed-`P,H` result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumGuiEnthalpySnapshot {
    target_enthalpy_j: f64,
    calculated_enthalpy_j: f64,
    enthalpy_error_j: f64,
    relative_enthalpy_error: f64,
}

impl EquilibriumGuiEnthalpySnapshot {
    pub fn target_enthalpy_j(self) -> f64 {
        self.target_enthalpy_j
    }

    pub fn calculated_enthalpy_j(self) -> f64 {
        self.calculated_enthalpy_j
    }

    pub fn enthalpy_error_j(self) -> f64 {
        self.enthalpy_error_j
    }

    pub fn relative_enthalpy_error(self) -> f64 {
        self.relative_enthalpy_error
    }
}

/// Immutable outer-solver evidence attached to a fixed-`P,H` result.
///
/// The GUI deliberately copies the small, user-relevant part of the engine
/// report instead of retaining a mutable workflow object. This keeps the
/// result safe to display after the worker has exited and makes fallback or
/// budget behavior visible rather than silently collapsing it into a scalar
/// energy error.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumGuiPhDiagnosticsSnapshot {
    solved_temperature_k: f64,
    solve_path: String,
    fallback_reason: Option<String>,
    trial_count: usize,
    iterations: usize,
    inner_backend_attempts: usize,
    inner_nonlinear_iterations: usize,
    phase_control_transitions: usize,
    fixed_formulation_builds: usize,
    fixed_formulation_reuses: usize,
    accepted_enthalpy_error_limit_j: f64,
    timing_enabled: bool,
    timing_total_ms: f64,
    route_decisions: Vec<EquilibriumGuiPhRouteDecisionSnapshot>,
    trials: Vec<EquilibriumGuiPhTrialSnapshot>,
    monolithic: Option<EquilibriumGuiPhMonolithicSnapshot>,
}

/// One immutable high-level P,H routing decision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumGuiPhRouteDecisionSnapshot {
    from_route: String,
    to_route: String,
    reason: String,
}

impl EquilibriumGuiPhRouteDecisionSnapshot {
    pub fn from_route(&self) -> &str {
        &self.from_route
    }

    pub fn to_route(&self) -> &str {
        &self.to_route
    }

    pub fn reason(&self) -> &str {
        &self.reason
    }
}

/// Immutable GUI projection of the coupled monolithic P,H evidence.
///
/// Monolithic solves do not have outer scalar-temperature trials. Keeping this
/// evidence in a separate route-specific snapshot prevents the UI from
/// inventing trial rows merely to make two numerically different workflows
/// look alike.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumGuiPhMonolithicSnapshot {
    backend_attempts: Vec<String>,
    accepted_backend: String,
    acceptance_rows: Vec<(String, String)>,
    phase_control_rows: Vec<(String, String)>,
    residual_evaluations: usize,
    jacobian_evaluations: usize,
    inner_timing_ms: f64,
}

impl EquilibriumGuiPhMonolithicSnapshot {
    pub fn backend_attempts(&self) -> &[String] {
        &self.backend_attempts
    }

    pub fn accepted_backend(&self) -> &str {
        &self.accepted_backend
    }

    pub fn acceptance_rows(&self) -> &[(String, String)] {
        &self.acceptance_rows
    }

    pub fn phase_control_rows(&self) -> &[(String, String)] {
        &self.phase_control_rows
    }

    pub fn residual_evaluations(&self) -> usize {
        self.residual_evaluations
    }

    pub fn jacobian_evaluations(&self) -> usize {
        self.jacobian_evaluations
    }

    pub fn inner_timing_ms(&self) -> f64 {
        self.inner_timing_ms
    }
}

/// One immutable row of the outer P,H temperature search.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumGuiPhTrialSnapshot {
    temperature_k: f64,
    step_kind: String,
    total_enthalpy_j: f64,
    enthalpy_error_j: f64,
    scaled_error: f64,
    inner_backend_attempts: usize,
    phase_control_transitions: usize,
    total_ms: f64,
}

impl EquilibriumGuiPhTrialSnapshot {
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    pub fn step_kind(&self) -> &str {
        &self.step_kind
    }

    pub fn total_enthalpy_j(&self) -> f64 {
        self.total_enthalpy_j
    }

    pub fn enthalpy_error_j(&self) -> f64 {
        self.enthalpy_error_j
    }

    pub fn scaled_error(&self) -> f64 {
        self.scaled_error
    }

    pub fn inner_backend_attempts(&self) -> usize {
        self.inner_backend_attempts
    }

    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_transitions
    }

    pub fn total_ms(&self) -> f64 {
        self.total_ms
    }
}

impl EquilibriumGuiPhDiagnosticsSnapshot {
    fn from_solution(solution: &FixedPressureEnthalpySolution) -> Self {
        let report = solution.report();
        let timing = report.timing();
        let trials = report
            .trials()
            .iter()
            .map(|trial| EquilibriumGuiPhTrialSnapshot {
                temperature_k: trial.temperature(),
                step_kind: format!("{:?}", trial.step_kind()),
                total_enthalpy_j: trial.total_enthalpy(),
                enthalpy_error_j: trial.enthalpy_error_joules(),
                scaled_error: trial.scaled_error(),
                inner_backend_attempts: trial.inner_backend_attempts(),
                phase_control_transitions: trial.phase_control_transitions(),
                total_ms: trial.timing().total().as_secs_f64() * 1_000.0,
            })
            .collect();
        let monolithic = report.monolithic_evidence().map(|evidence| {
            let acceptance_rows = evidence
                .acceptance_report()
                .map(|report| {
                    report
                        .summary_rows()
                        .into_iter()
                        .map(|row| (row.label, row.value))
                        .collect()
                })
                .unwrap_or_default();
            let phase_control_rows = evidence
                .phase_control_report()
                .map(|report| {
                    report
                        .summary_rows()
                        .into_iter()
                        .map(|row| (row.label, row.value))
                        .collect()
                })
                .unwrap_or_default();
            EquilibriumGuiPhMonolithicSnapshot {
                backend_attempts: evidence
                    .solve_report()
                    .attempts
                    .iter()
                    .map(|attempt| attempt.summary())
                    .collect(),
                accepted_backend: format!("{:?}", evidence.solve_report().accepted_backend),
                acceptance_rows,
                phase_control_rows,
                residual_evaluations: evidence.residual_evaluations(),
                jacobian_evaluations: evidence.jacobian_evaluations(),
                inner_timing_ms: evidence.inner_timing().total().as_secs_f64() * 1_000.0,
            }
        });
        let route_decisions = report
            .route_decisions()
            .iter()
            .map(EquilibriumGuiPhRouteDecisionSnapshot::from_engine)
            .collect();
        Self {
            solved_temperature_k: solution.temperature(),
            solve_path: format!("{:?}", report.solve_path()),
            fallback_reason: report
                .fallback_reason()
                .map(|reason| format!("{:?}: {}", reason.error_kind(), reason.message())),
            trial_count: report.trials().len(),
            iterations: report.iterations(),
            inner_backend_attempts: report.inner_backend_attempts(),
            inner_nonlinear_iterations: report.inner_nonlinear_iterations(),
            phase_control_transitions: report.phase_control_transitions(),
            fixed_formulation_builds: report.fixed_formulation_builds(),
            fixed_formulation_reuses: report.fixed_formulation_reuses(),
            accepted_enthalpy_error_limit_j: solution.enthalpy_error_limit_joules(),
            timing_enabled: timing.enabled(),
            timing_total_ms: timing.total().as_secs_f64() * 1_000.0,
            route_decisions,
            trials,
            monolithic,
        }
    }

    pub fn solved_temperature_k(&self) -> f64 {
        self.solved_temperature_k
    }

    pub fn solve_path(&self) -> &str {
        &self.solve_path
    }

    pub fn fallback_reason(&self) -> Option<&str> {
        self.fallback_reason.as_deref()
    }

    pub fn trial_count(&self) -> usize {
        self.trial_count
    }

    pub fn iterations(&self) -> usize {
        self.iterations
    }

    pub fn inner_backend_attempts(&self) -> usize {
        self.inner_backend_attempts
    }

    pub fn inner_nonlinear_iterations(&self) -> usize {
        self.inner_nonlinear_iterations
    }

    pub fn phase_control_transitions(&self) -> usize {
        self.phase_control_transitions
    }

    pub fn fixed_formulation_builds(&self) -> usize {
        self.fixed_formulation_builds
    }

    pub fn fixed_formulation_reuses(&self) -> usize {
        self.fixed_formulation_reuses
    }

    pub fn accepted_enthalpy_error_limit_j(&self) -> f64 {
        self.accepted_enthalpy_error_limit_j
    }

    pub fn timing_enabled(&self) -> bool {
        self.timing_enabled
    }

    pub fn timing_total_ms(&self) -> f64 {
        self.timing_total_ms
    }

    pub fn trials(&self) -> &[EquilibriumGuiPhTrialSnapshot] {
        &self.trials
    }

    /// High-level route evidence is independent from scalar temperature trials.
    pub fn route_decisions(&self) -> &[EquilibriumGuiPhRouteDecisionSnapshot] {
        &self.route_decisions
    }

    /// Coupled evidence, present only for monolithic P,H routes.
    pub fn monolithic(&self) -> Option<&EquilibriumGuiPhMonolithicSnapshot> {
        self.monolithic.as_ref()
    }
}

impl EquilibriumGuiPhRouteDecisionSnapshot {
    fn from_engine(decision: &PhRouteDecision) -> Self {
        Self {
            from_route: format!("{:?}", decision.from_route()),
            to_route: format!("{:?}", decision.to_route()),
            reason: format!(
                "{:?}: {}",
                decision.reason().error_kind(),
                decision.reason().message()
            ),
        }
    }
}

fn lifecycle_event_snapshot(
    event: &EquilibriumDiagnosticEvent,
    solution: &MultiphaseEquilibriumSolution,
) -> EquilibriumGuiLifecycleEventSnapshot {
    use EquilibriumDiagnosticEvent as Event;

    let (title, detail) = match event {
        Event::SolveStarted {
            conditions,
            initial_phase_set,
        } => (
            "Solve started".to_string(),
            format!(
                "T={:.6} K, P={:.6} Pa, active=[{}]",
                conditions.temperature(),
                conditions.pressure(),
                lifecycle_phase_set(initial_phase_set.active_mask(), solution),
            ),
        ),
        Event::OuterIterationStarted {
            iteration,
            active_phase_set,
        } => (
            format!("Outer iteration {iteration}"),
            format!(
                "active=[{}]",
                lifecycle_phase_set(active_phase_set.active_mask(), solution)
            ),
        ),
        Event::ActiveSetCandidateAccepted {
            iteration,
            validation,
            backend_summary,
            ..
        } => (
            format!("Candidate accepted at iteration {iteration}"),
            format!(
                "residual={:.3e}, balance={:.3e}; {backend_summary}",
                validation.residual_l2_norm, validation.max_abs_element_balance_error
            ),
        ),
        Event::ActiveSetCandidateRejected {
            iteration,
            active_phase_set,
            message,
        } => (
            format!("Candidate rejected at iteration {iteration}"),
            format!(
                "active=[{}]; {message}",
                lifecycle_phase_set(active_phase_set.active_mask(), solution)
            ),
        ),
        Event::StabilityEvaluated {
            iteration,
            dg_create_j_per_mol,
            dg_keep_j_per_mol,
            phases,
        } => (
            format!("TPD stability at iteration {iteration}"),
            format!(
                "dg_create={dg_create_j_per_mol:.3e} J/mol, dg_keep={dg_keep_j_per_mol:.3e} J/mol; {}",
                lifecycle_stability_summary(phases, solution)
            ),
        ),
        Event::TransitionAccepted {
            iteration,
            activated_phase_indices,
            deactivated_phase_indices,
            reason,
            previous_phase_set,
            new_phase_set,
            ..
        } => (
            format!("Phase transition accepted at iteration {iteration}"),
            format!(
                "[{}] -> [{}]; activated=[{}], deactivated=[{}], reason={reason:?}",
                lifecycle_phase_set(previous_phase_set.active_mask(), solution),
                lifecycle_phase_set(new_phase_set.active_mask(), solution),
                lifecycle_phase_indices(activated_phase_indices, solution),
                lifecycle_phase_indices(deactivated_phase_indices, solution),
            ),
        ),
        Event::TransitionHeldByHysteresis {
            iteration,
            phase_index,
        } => (
            format!("Hysteresis hold at iteration {iteration}"),
            format!("retained {}", lifecycle_phase_label(*phase_index, solution)),
        ),
        Event::RecoveryProbeStarted {
            iteration,
            removed_phase_index,
            attempted_phase_set,
        } => (
            format!("Boundary recovery started at iteration {iteration}"),
            format!(
                "remove {}; active=[{}]",
                lifecycle_phase_label(*removed_phase_index, solution),
                lifecycle_phase_set(attempted_phase_set.active_mask(), solution),
            ),
        ),
        Event::RecoveryProbeAccepted {
            iteration,
            phase_index,
            previous_phase_set,
            new_phase_set,
        } => (
            format!("Boundary recovery accepted at iteration {iteration}"),
            format!(
                "removed {}; [{}] -> [{}]",
                lifecycle_phase_label(*phase_index, solution),
                lifecycle_phase_set(previous_phase_set.active_mask(), solution),
                lifecycle_phase_set(new_phase_set.active_mask(), solution),
            ),
        ),
        Event::RecoveryProbeRejected {
            iteration,
            removed_phase_index,
            message,
        } => (
            format!("Boundary recovery rejected at iteration {iteration}"),
            format!(
                "remove {}; {message}",
                lifecycle_phase_label(*removed_phase_index, solution)
            ),
        ),
        Event::ContinuationRestored {
            retained_phase_set,
            retained_seed,
        } => (
            "Continuation restored".to_string(),
            format!(
                "seed_retained={retained_seed}, active=[{}]",
                retained_phase_set
                    .as_ref()
                    .map(|set| lifecycle_phase_set(set.active_mask(), solution))
                    .unwrap_or_else(|| "none".into())
            ),
        ),
        Event::PhaseControlBudgetExhausted {
            max_outer_iterations,
        } => (
            "Phase-control budget exhausted".to_string(),
            format!("max_outer_iterations={max_outer_iterations}"),
        ),
        Event::PhaseControlCycleDetected {
            iteration,
            repeated_phase_set,
        } => (
            format!("Phase-control cycle at iteration {iteration}"),
            format!(
                "repeated active=[{}]",
                lifecycle_phase_set(repeated_phase_set.active_mask(), solution)
            ),
        ),
        Event::PhRouteFallback {
            from_route,
            to_route,
            error_kind,
            message,
        } => (
            "P,H route fallback".to_string(),
            format!("{from_route:?} -> {to_route:?}; {error_kind:?}: {message}"),
        ),
        Event::PhRouteFailed {
            route,
            error_kind,
            message,
        } => (
            "P,H route failed".to_string(),
            format!("{route:?}; {error_kind:?}: {message}"),
        ),
        Event::SolveFailed {
            message,
            continuation_restored,
        } => (
            "Solve failed".to_string(),
            format!("continuation_restored={continuation_restored}; {message}"),
        ),
        Event::SolveAccepted {
            final_phase_set,
            outer_iterations,
            transition_count,
        } => (
            "Solve accepted".to_string(),
            format!(
                "active=[{}], outer_iterations={outer_iterations}, transitions={transition_count}",
                lifecycle_phase_set(final_phase_set.active_mask(), solution)
            ),
        ),
    };
    EquilibriumGuiLifecycleEventSnapshot { title, detail }
}

fn lifecycle_stability_summary(
    phases: &[PhaseStabilityDiagnostic],
    solution: &MultiphaseEquilibriumSolution,
) -> String {
    phases
        .iter()
        .map(|phase| {
            let tpd = phase
                .minimum_tpd_j_per_mol
                .map(|value| format!("{value:.3e} J/mol"))
                .unwrap_or_else(|| "not evaluated".into());
            format!(
                "{}: {}, total={:.3e} mol, minimum_tpd={tpd}",
                lifecycle_phase_label(phase.phase_index, solution),
                if phase.active { "active" } else { "inactive" },
                phase.phase_total_moles,
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn lifecycle_phase_set(active_mask: Vec<bool>, solution: &MultiphaseEquilibriumSolution) -> String {
    active_mask
        .iter()
        .enumerate()
        .filter_map(|(index, active)| active.then(|| lifecycle_phase_label(index, solution)))
        .collect::<Vec<_>>()
        .join(", ")
}

fn lifecycle_phase_indices(
    phase_indices: &[usize],
    solution: &MultiphaseEquilibriumSolution,
) -> String {
    phase_indices
        .iter()
        .map(|&index| lifecycle_phase_label(index, solution))
        .collect::<Vec<_>>()
        .join(", ")
}

fn lifecycle_phase_label(index: usize, solution: &MultiphaseEquilibriumSolution) -> String {
    solution
        .phases()
        .get(index)
        .and_then(|phase| phase.id().as_option().clone())
        .unwrap_or_else(|| format!("phase_{index}"))
}

impl EquilibriumGuiPointSnapshot {
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    pub fn component_moles(&self) -> &[f64] {
        &self.component_moles
    }

    pub fn mole_fractions(&self) -> &[f64] {
        &self.mole_fractions
    }

    pub fn phase_totals(&self) -> &[f64] {
        &self.phase_totals
    }

    pub fn phase_statuses(&self) -> &[String] {
        &self.phase_statuses
    }

    /// Active mask aligned with the component order. This is lifecycle
    /// metadata, not a filter over physical zero values.
    pub fn component_active(&self) -> &[bool] {
        &self.component_active
    }

    /// Active mask aligned with the phase order.
    pub fn phase_active(&self) -> &[bool] {
        &self.phase_active
    }

    /// Optional bounded lifecycle trace requested by the GUI document.
    pub fn lifecycle_trace(&self) -> Option<&EquilibriumGuiLifecycleTraceSnapshot> {
        self.lifecycle_trace.as_ref()
    }

    /// Authoritative solver reports and lookup provenance for this point.
    pub fn source(&self) -> &MultiphaseEquilibriumSolution {
        self.source.as_ref()
    }
}

/// Transactionally accepted UI result for one point or a temperature range.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumGuiResultSnapshot {
    component_labels: Vec<String>,
    phase_labels: Vec<String>,
    points: Vec<EquilibriumGuiPointSnapshot>,
    range_report: Option<TemperatureRangeSolveReport>,
    enthalpy: Option<EquilibriumGuiEnthalpySnapshot>,
    ph_diagnostics: Option<EquilibriumGuiPhDiagnosticsSnapshot>,
}

impl EquilibriumGuiResultSnapshot {
    /// Creates a UI snapshot only from an already accepted engine outcome.
    ///
    /// The conversion is fallible because a range result must be structurally
    /// homogeneous: every point must have the same phase-qualified layout and
    /// component ordering. A malformed worker payload must become a normal GUI
    /// failure, never a panic and never a partially published table.
    pub fn from_outcome(outcome: EquilibriumGuiSolveOutcome) -> Result<Self, String> {
        match outcome {
            EquilibriumGuiSolveOutcome::Point(solution) => {
                let source = Arc::new(solution.solution().clone());
                Self::from_sources(vec![source], None, None, None)
            }
            EquilibriumGuiSolveOutcome::Range(range) => {
                let sources = range
                    .points()
                    .iter()
                    .map(|point| Arc::new(point.solution().clone()))
                    .collect();
                Self::from_sources(sources, Some(range.report().clone()), None, None)
            }
            EquilibriumGuiSolveOutcome::Ph(solution) => {
                let source = Arc::new(solution.equilibrium().clone());
                let ph_diagnostics = EquilibriumGuiPhDiagnosticsSnapshot::from_solution(&solution);
                let enthalpy = EquilibriumGuiEnthalpySnapshot {
                    target_enthalpy_j: solution.target_enthalpy(),
                    calculated_enthalpy_j: solution.calculated_enthalpy(),
                    enthalpy_error_j: solution.enthalpy_error(),
                    relative_enthalpy_error: solution.enthalpy_error().abs()
                        / solution.target_enthalpy().abs().max(1.0),
                };
                Self::from_sources(vec![source], None, Some(enthalpy), Some(ph_diagnostics))
            }
        }
    }

    fn from_sources(
        sources: Vec<Arc<MultiphaseEquilibriumSolution>>,
        range_report: Option<TemperatureRangeSolveReport>,
        enthalpy: Option<EquilibriumGuiEnthalpySnapshot>,
        ph_diagnostics: Option<EquilibriumGuiPhDiagnosticsSnapshot>,
    ) -> Result<Self, String> {
        let first = sources
            .first()
            .ok_or_else(|| "accepted equilibrium outcome contains no points".to_string())?;
        validate_accepted_candidate(first, 0)?;
        let component_labels = first
            .metadata()
            .components()
            .iter()
            .map(|component| component.label())
            .collect::<Vec<_>>();
        let phase_labels = first
            .metadata()
            .phases()
            .iter()
            .map(|phase| {
                phase
                    .id()
                    .as_option()
                    .clone()
                    .unwrap_or_else(|| "<anonymous>".into())
            })
            .collect::<Vec<_>>();
        let layout_fingerprint = first.layout_fingerprint();
        for (point_index, source) in sources.iter().enumerate().skip(1) {
            validate_accepted_candidate(source, point_index)?;
            if source.layout_fingerprint() != layout_fingerprint {
                return Err(format!(
                    "range result point {point_index} has a different layout fingerprint"
                ));
            }
            let point_components = source
                .metadata()
                .components()
                .iter()
                .map(|component| component.label())
                .collect::<Vec<_>>();
            if point_components != component_labels {
                return Err(format!(
                    "range result point {point_index} has a different component ordering"
                ));
            }
            let point_phases = source
                .metadata()
                .phases()
                .iter()
                .map(|phase| {
                    phase
                        .id()
                        .as_option()
                        .clone()
                        .unwrap_or_else(|| "<anonymous>".into())
                })
                .collect::<Vec<_>>();
            if point_phases != phase_labels {
                return Err(format!(
                    "range result point {point_index} has a different phase ordering"
                ));
            }
        }
        let points = sources
            .into_iter()
            .map(|source| {
                let component_ids = source
                    .metadata()
                    .components()
                    .iter()
                    .map(|component| component.id().clone())
                    .collect::<Vec<_>>();
                let phase_ids = source
                    .metadata()
                    .phases()
                    .iter()
                    .map(|phase| phase.id().clone())
                    .collect::<Vec<_>>();
                let phase_statuses = phase_ids
                    .iter()
                    .map(|phase| {
                        source
                            .phase_status(phase)
                            .map(|status| format!("{status:?}"))
                            .unwrap_or_else(|| "Unknown".into())
                    })
                    .collect();
                let phase_active = phase_ids
                    .iter()
                    .map(|phase| {
                        matches!(
                            source.phase_status(phase),
                            Some(PhaseStatus::Active | PhaseStatus::Appeared)
                        )
                    })
                    .collect::<Vec<_>>();
                let component_active = component_ids
                    .iter()
                    .map(|component| {
                        matches!(
                            source.phase_status(&component.phase),
                            Some(PhaseStatus::Active | PhaseStatus::Appeared)
                        )
                    })
                    .collect();
                let component_moles = source.component_moles().to_vec();
                let mole_fractions = component_ids
                    .iter()
                    .map(|component| source.mole_fraction_for(component).unwrap_or(0.0))
                    .collect();
                let phase_totals = phase_ids
                    .iter()
                    .map(|phase| source.phase_total(phase).unwrap_or(0.0))
                    .collect();
                let lifecycle_trace = EquilibriumGuiLifecycleTraceSnapshot::from_solution(&source);
                EquilibriumGuiPointSnapshot {
                    temperature_k: source.conditions().temperature(),
                    component_moles,
                    mole_fractions,
                    phase_totals,
                    phase_statuses,
                    component_active,
                    phase_active,
                    lifecycle_trace,
                    source,
                }
            })
            .collect();

        Ok(Self {
            component_labels,
            phase_labels,
            points,
            range_report,
            enthalpy,
            ph_diagnostics,
        })
    }

    pub fn component_labels(&self) -> &[String] {
        &self.component_labels
    }

    pub fn phase_labels(&self) -> &[String] {
        &self.phase_labels
    }

    /// Points in exact solver-grid order. A point solve contains one point.
    pub fn points(&self) -> &[EquilibriumGuiPointSnapshot] {
        &self.points
    }

    pub fn range_report(&self) -> Option<&TemperatureRangeSolveReport> {
        self.range_report.as_ref()
    }

    pub fn is_range(&self) -> bool {
        self.range_report.is_some()
    }

    /// Returns the energy contract for a `P,H` result, if this snapshot came
    /// from that problem family.
    pub fn enthalpy(&self) -> Option<EquilibriumGuiEnthalpySnapshot> {
        self.enthalpy
    }

    /// Returns outer-solver evidence for a fixed-`P,H` result.
    pub fn ph_diagnostics(&self) -> Option<EquilibriumGuiPhDiagnosticsSnapshot> {
        self.ph_diagnostics.clone()
    }

    /// Returns one aligned series for plotting or table rendering.
    pub fn component_mole_series(&self, component_index: usize) -> Option<Vec<f64>> {
        if component_index >= self.component_labels.len() {
            return None;
        }
        Some(
            self.points
                .iter()
                .map(|point| point.component_moles[component_index])
                .collect(),
        )
    }

    pub fn component_fraction_series(&self, component_index: usize) -> Option<Vec<f64>> {
        if component_index >= self.component_labels.len() {
            return None;
        }
        Some(
            self.points
                .iter()
                .map(|point| point.mole_fractions[component_index])
                .collect(),
        )
    }

    pub fn phase_total_series(&self, phase_index: usize) -> Option<Vec<f64>> {
        if phase_index >= self.phase_labels.len() {
            return None;
        }
        Some(
            self.points
                .iter()
                .map(|point| point.phase_totals[phase_index])
                .collect(),
        )
    }

    pub fn temperatures(&self) -> Vec<f64> {
        self.points
            .iter()
            .map(|point| point.temperature_k)
            .collect()
    }
}

/// Validates the small set of invariants the GUI relies on when displaying an
/// accepted engine result.
///
/// The solver owns the authoritative acceptance decision. This is a second,
/// cheap publication-boundary check: it catches a malformed worker payload in
/// which the candidate is marked accepted but the copied validation report no
/// longer describes the physical moles. The GUI must reject that payload
/// transactionally instead of rendering a trustworthy-looking table.
fn validate_accepted_candidate(
    source: &MultiphaseEquilibriumSolution,
    point_index: usize,
) -> Result<(), String> {
    let validation = source.accepted_solution().validation();
    let metrics = [
        validation.residual_l2_norm,
        validation.residual_rms,
        validation.max_abs_residual,
        validation.raw_residual_l2_norm,
        validation.raw_residual_rms,
        validation.raw_max_abs_residual,
        validation.max_abs_element_balance_error,
        validation.reaction_affinity_l2_norm,
        validation.max_abs_reaction_affinity,
        validation.min_moles,
    ];
    if metrics
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(format!(
            "accepted result point {point_index} has non-finite or negative validation evidence"
        ));
    }

    let physical_min_moles = source
        .accepted_solution()
        .moles()
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let tolerance = 1.0e-12 * physical_min_moles.abs().max(1.0);
    if (validation.min_moles - physical_min_moles).abs() > tolerance {
        return Err(format!(
            "accepted result point {point_index} has validation evidence inconsistent with physical moles"
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::PhaseEquilibriumPipelineRequest;
    use crate::Thermodynamics::ChemEquilibrium::prelude::{
        SubstanceSystemFactory, SubstanceSystemSpecBuilder, SubstancesContainer,
    };
    use crate::Thermodynamics::thermo_lib_api::ThermoData;

    fn accepted_real_solution(
        names: &[&str],
        initial_moles: Vec<f64>,
    ) -> Arc<MultiphaseEquilibriumSolution> {
        let repository = ThermoData::try_default_repository()
            .expect("the bundled thermochemistry repository must be available");
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(
            names.iter().map(|name| (*name).to_string()).collect(),
        ))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("the real GUI mismatch fixture must build");
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            spec.clone(),
            repository.clone(),
        )
        .expect("the real GUI mismatch fixture must resolve");
        let outcome = PhaseEquilibriumPipelineRequest::new(
            spec,
            initial_moles,
            EquilibriumConditions::new(2_500.0, 101_325.0, 101_325.0)
                .expect("fixture conditions must validate"),
        )
        .with_repository(repository)
        .solve()
        .expect("the real GUI mismatch fixture must solve");
        assert_eq!(outcome.resolved().phase_specs(), resolved.phase_specs());
        Arc::new(outcome.into_solution())
    }

    #[test]
    fn source_snapshot_requires_an_accepted_point() {
        let error = EquilibriumGuiResultSnapshot::from_sources(Vec::new(), None, None, None)
            .expect_err("an empty worker outcome must be rejected");
        assert!(error.contains("contains no points"));
    }

    #[test]
    #[ignore = "requires the bundled NASA gas catalog"]
    fn source_snapshot_rejects_an_accepted_range_with_mismatched_layout() {
        let first = accepted_real_solution(&["H2", "O2", "H2O"], vec![0.1, 0.05, 1.9]);
        let second = accepted_real_solution(&["H2", "O2"], vec![0.1, 0.05]);

        let error =
            EquilibriumGuiResultSnapshot::from_sources(vec![first, second], None, None, None)
                .expect_err("a range with incompatible accepted layouts must be rejected");
        assert!(
            error.contains("different layout fingerprint"),
            "unexpected mismatch error: {error}"
        );
    }

    #[test]
    #[ignore = "requires the bundled NASA gas catalog"]
    fn source_snapshot_rejects_an_accepted_candidate_with_mismatched_validation() {
        let source = accepted_real_solution(&["H2", "O2", "H2O"], vec![0.1, 0.05, 1.9]);
        let mut validation = source.accepted_solution().validation().clone();
        validation.min_moles *= 2.0;
        let malformed = Arc::new(source.as_ref().clone().with_validation_for_test(validation));

        let error = EquilibriumGuiResultSnapshot::from_sources(vec![malformed], None, None, None)
            .expect_err("a candidate with mismatched validation must not publish");
        assert!(
            error.contains("validation evidence inconsistent with physical moles"),
            "unexpected validation mismatch error: {error}"
        );
    }
}
