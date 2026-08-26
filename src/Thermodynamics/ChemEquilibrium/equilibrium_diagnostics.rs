//! Optional structured diagnostics for canonical equilibrium workflows.
//!
//! The nonlinear and phase-control layers emit typed facts, not formatted
//! strings or direct `log` calls.  A caller may retain the immutable report,
//! attach a live sink for a CLI or GUI, or render it afterwards.  Diagnostics
//! are disabled by default so ordinary solves do not allocate event payloads.

use std::fmt;
use std::sync::Arc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentErrorKind;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    PhaseSet, PhaseTransitionReason,
};

/// Requested amount of diagnostic evidence for one equilibrium transaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EquilibriumDiagnosticsMode {
    /// Emit and retain no diagnostic events.
    #[default]
    Disabled,
    /// Retain transaction boundaries and the final decision only.
    Summary,
    /// Include the phase-control lifecycle, TPD evaluations, and transitions.
    PhaseLifecycle,
    /// Include accepted fixed-set candidates and recovery probes as well.
    Detailed,
}

impl EquilibriumDiagnosticsMode {
    pub(crate) fn captures_lifecycle(self) -> bool {
        matches!(self, Self::PhaseLifecycle | Self::Detailed)
    }
}

/// Selects which points of a temperature range emit diagnostic events.
///
/// Point-level diagnostics can be verbose even when each individual trace is
/// bounded. The default records only the two endpoints; callers must opt in
/// explicitly before streaming every point to a live sink.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EquilibriumRangeDiagnosticsPolicy {
    /// Retain and stream a point trace only after that accepted point changed
    /// the active phase set. The range workflow defers its live sink until it
    /// knows whether a transition occurred.
    TransitionsOnly,
    /// Keep diagnostics for the first and last point only.
    #[default]
    Endpoints,
    /// Keep diagnostics for point zero and every `stride`-th point.
    EveryNth { stride: usize },
    /// Keep diagnostics for every solved point.
    EveryPoint,
}

/// Controls bounded structured diagnostic collection.
///
/// `max_events` applies only to the immutable report. A live sink still sees
/// events after the retained report reaches its limit, which lets a caller
/// stream a long temperature range without retaining an unbounded trace.
#[derive(Clone)]
pub struct EquilibriumDiagnosticsOptions {
    mode: EquilibriumDiagnosticsMode,
    max_events: usize,
    range_policy: EquilibriumRangeDiagnosticsPolicy,
    sink: Option<EquilibriumDiagnosticSink>,
    defer_live_sink: bool,
}

impl fmt::Debug for EquilibriumDiagnosticsOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EquilibriumDiagnosticsOptions")
            .field("mode", &self.mode)
            .field("max_events", &self.max_events)
            .field("range_policy", &self.range_policy)
            .field("has_sink", &self.sink.is_some())
            .field("defer_live_sink", &self.defer_live_sink)
            .finish()
    }
}

impl Default for EquilibriumDiagnosticsOptions {
    fn default() -> Self {
        Self {
            mode: EquilibriumDiagnosticsMode::Disabled,
            max_events: 256,
            range_policy: EquilibriumRangeDiagnosticsPolicy::default(),
            sink: None,
            defer_live_sink: false,
        }
    }
}

impl EquilibriumDiagnosticsOptions {
    /// Starts with diagnostics disabled.
    pub fn disabled() -> Self {
        Self::default()
    }

    /// Enables one diagnostic detail level with a bounded in-result trace.
    pub fn enabled(mode: EquilibriumDiagnosticsMode) -> Self {
        Self {
            mode,
            ..Self::default()
        }
    }

    /// Selects the retained diagnostic detail level.
    pub fn with_mode(mut self, mode: EquilibriumDiagnosticsMode) -> Self {
        self.mode = mode;
        self
    }

    /// Limits event retention in the immutable result report.
    pub fn with_max_events(mut self, max_events: usize) -> Self {
        self.max_events = max_events;
        self
    }

    /// Controls which points in a temperature range emit diagnostics.
    pub fn with_range_policy(mut self, policy: EquilibriumRangeDiagnosticsPolicy) -> Self {
        self.range_policy = policy;
        self
    }

    /// Receives events as they are produced by the canonical workflow.
    ///
    /// The sink is observational: it must not mutate solver state or depend
    /// on a partially published numerical result.
    pub fn with_sink<F>(mut self, sink: F) -> Self
    where
        F: Fn(EquilibriumDiagnosticEvent) + Send + Sync + 'static,
    {
        self.sink = Some(Arc::new(sink));
        self
    }

    pub fn mode(&self) -> EquilibriumDiagnosticsMode {
        self.mode
    }

    pub fn max_events(&self) -> usize {
        self.max_events
    }

    pub fn range_policy(&self) -> EquilibriumRangeDiagnosticsPolicy {
        self.range_policy
    }

    pub fn has_sink(&self) -> bool {
        self.sink.is_some()
    }

    /// Derives the policy for one zero-based point in a finite range.
    ///
    /// Disabled points suppress the live sink too, rather than merely dropping
    /// report entries after verbose logging has already happened.
    pub(crate) fn for_range_point(&self, index: usize, point_count: usize) -> Self {
        if self.mode == EquilibriumDiagnosticsMode::Disabled {
            return self.clone();
        }
        let selected = match self.range_policy {
            EquilibriumRangeDiagnosticsPolicy::TransitionsOnly => {
                let mut deferred = self.clone();
                deferred.defer_live_sink = true;
                return deferred;
            }
            EquilibriumRangeDiagnosticsPolicy::Endpoints => {
                index == 0 || index.saturating_add(1) == point_count
            }
            EquilibriumRangeDiagnosticsPolicy::EveryNth { stride } => {
                stride > 0 && (index == 0 || index % stride == 0)
            }
            EquilibriumRangeDiagnosticsPolicy::EveryPoint => true,
        };
        if selected {
            self.clone()
        } else {
            Self::disabled()
        }
    }

    /// Streams a workflow-level observation that has no accepted solution to
    /// carry an immutable report yet, such as an `Auto` P,H route fallback.
    pub(crate) fn emit(
        &self,
        minimum: EquilibriumDiagnosticsMode,
        event: EquilibriumDiagnosticEvent,
    ) {
        if self.mode == EquilibriumDiagnosticsMode::Disabled || !mode_includes(self.mode, minimum) {
            return;
        }
        if !self.defer_live_sink {
            if let Some(sink) = &self.sink {
                sink(event);
            }
        }
    }

    /// Replays one retained accepted-point report after a range workflow
    /// confirms that the point satisfies its deferred publication rule.
    pub(crate) fn replay_retained_events(&self, report: &EquilibriumDiagnosticsReport) {
        if self.mode == EquilibriumDiagnosticsMode::Disabled || !self.defer_live_sink {
            return;
        }
        if let Some(sink) = &self.sink {
            for event in report.events() {
                sink(event.clone());
            }
        }
    }
}

/// Live observer for one diagnostic event.
pub type EquilibriumDiagnosticSink = Arc<dyn Fn(EquilibriumDiagnosticEvent) + Send + Sync>;

/// Compact TPD evidence for one declared phase.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseStabilityDiagnostic {
    pub phase_index: usize,
    pub active: bool,
    pub phase_total_moles: f64,
    pub minimum_tpd_j_per_mol: Option<f64>,
    pub incipient_composition: Option<Vec<f64>>,
}

/// Chronological fact emitted by the canonical phase-control workflow.
#[derive(Debug, Clone, PartialEq)]
pub enum EquilibriumDiagnosticEvent {
    SolveStarted {
        conditions: EquilibriumConditions,
        initial_phase_set: PhaseSet,
    },
    OuterIterationStarted {
        iteration: usize,
        active_phase_set: PhaseSet,
    },
    ActiveSetCandidateAccepted {
        iteration: usize,
        active_phase_set: PhaseSet,
        validation: EquilibriumCandidateReport,
        backend_summary: String,
    },
    /// A fixed active-set solve was rejected before it could contribute
    /// stability evidence. The lifecycle may still attempt a boundary
    /// recovery with a smaller active set.
    ActiveSetCandidateRejected {
        iteration: usize,
        active_phase_set: PhaseSet,
        message: String,
    },
    StabilityEvaluated {
        iteration: usize,
        dg_create_j_per_mol: f64,
        dg_keep_j_per_mol: f64,
        phases: Vec<PhaseStabilityDiagnostic>,
    },
    TransitionAccepted {
        iteration: usize,
        activated_phase_indices: Vec<usize>,
        deactivated_phase_indices: Vec<usize>,
        reason: PhaseTransitionReason,
        previous_phase_set: PhaseSet,
        new_phase_set: PhaseSet,
        incipient_composition: Option<Vec<f64>>,
    },
    TransitionHeldByHysteresis {
        iteration: usize,
        phase_index: usize,
    },
    RecoveryProbeAccepted {
        iteration: usize,
        phase_index: usize,
        previous_phase_set: PhaseSet,
        new_phase_set: PhaseSet,
    },
    /// A boundary recovery is trying a smaller active phase set after the
    /// primary fixed-set solve failed.
    RecoveryProbeStarted {
        iteration: usize,
        removed_phase_index: usize,
        attempted_phase_set: PhaseSet,
    },
    /// A boundary recovery probe did not establish a valid physical boundary.
    RecoveryProbeRejected {
        iteration: usize,
        removed_phase_index: usize,
        message: String,
    },
    /// A failed transaction restored the last accepted continuation state.
    ContinuationRestored {
        retained_phase_set: Option<PhaseSet>,
        retained_seed: bool,
    },
    /// The bounded outer loop exhausted its transition budget.
    PhaseControlBudgetExhausted { max_outer_iterations: usize },
    /// A proposed transition revisited a settled active phase set, so the
    /// outer loop stopped rather than oscillating indefinitely.
    PhaseControlCycleDetected {
        iteration: usize,
        repeated_phase_set: PhaseSet,
    },
    /// `Auto` P,H abandoned its coupled formulation after a classified error
    /// and began the independent nested-temperature route.
    PhRouteFallback {
        from_route: PhDiagnosticRoute,
        to_route: PhDiagnosticRoute,
        error_kind: ReactionExtentErrorKind,
        message: String,
    },
    /// A P,H route ended before it could publish an accepted solution.
    PhRouteFailed {
        route: PhDiagnosticRoute,
        error_kind: ReactionExtentErrorKind,
        message: String,
    },
    /// The transaction ended without publishing a new accepted solution.
    ///
    /// This event is mainly useful to a live sink because failed transactions
    /// intentionally do not return an immutable solution report.
    SolveFailed {
        message: String,
        continuation_restored: bool,
    },
    SolveAccepted {
        final_phase_set: PhaseSet,
        outer_iterations: usize,
        transition_count: usize,
    },
}

/// Coarse P,H route identity used in diagnostics without exposing workflow
/// implementation types to the generic phase-control trace.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhDiagnosticRoute {
    Monolithic,
    NestedTemperature,
}

/// Immutable bounded diagnostic evidence attached to an accepted result.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumDiagnosticsReport {
    mode: EquilibriumDiagnosticsMode,
    events: Vec<EquilibriumDiagnosticEvent>,
    dropped_events: usize,
}

impl EquilibriumDiagnosticsReport {
    pub fn mode(&self) -> EquilibriumDiagnosticsMode {
        self.mode
    }

    pub fn enabled(&self) -> bool {
        self.mode != EquilibriumDiagnosticsMode::Disabled
    }

    pub fn events(&self) -> &[EquilibriumDiagnosticEvent] {
        &self.events
    }

    pub fn dropped_events(&self) -> usize {
        self.dropped_events
    }
}

/// Mutable collector kept inside one phase-control transaction.
#[derive(Debug, Clone)]
pub(crate) struct EquilibriumDiagnosticsCollector {
    options: EquilibriumDiagnosticsOptions,
    events: Vec<EquilibriumDiagnosticEvent>,
    dropped_events: usize,
}

impl EquilibriumDiagnosticsCollector {
    pub(crate) fn new(options: EquilibriumDiagnosticsOptions) -> Self {
        Self {
            options,
            events: Vec::new(),
            dropped_events: 0,
        }
    }

    pub(crate) fn set_options(&mut self, options: EquilibriumDiagnosticsOptions) {
        self.options = options;
        self.events.clear();
        self.dropped_events = 0;
    }

    pub(crate) fn mode(&self) -> EquilibriumDiagnosticsMode {
        self.options.mode
    }

    pub(crate) fn record(
        &mut self,
        minimum: EquilibriumDiagnosticsMode,
        event: EquilibriumDiagnosticEvent,
    ) {
        if self.options.mode == EquilibriumDiagnosticsMode::Disabled
            || !mode_includes(self.options.mode, minimum)
        {
            return;
        }
        self.options.emit(minimum, event.clone());
        if self.events.len() < self.options.max_events {
            self.events.push(event);
        } else {
            self.dropped_events += 1;
        }
    }

    pub(crate) fn finish(&self) -> EquilibriumDiagnosticsReport {
        EquilibriumDiagnosticsReport {
            mode: self.options.mode,
            events: self.events.clone(),
            dropped_events: self.dropped_events,
        }
    }
}

fn mode_includes(mode: EquilibriumDiagnosticsMode, required: EquilibriumDiagnosticsMode) -> bool {
    fn rank(mode: EquilibriumDiagnosticsMode) -> u8 {
        match mode {
            EquilibriumDiagnosticsMode::Disabled => 0,
            EquilibriumDiagnosticsMode::Summary => 1,
            EquilibriumDiagnosticsMode::PhaseLifecycle => 2,
            EquilibriumDiagnosticsMode::Detailed => 3,
        }
    }
    rank(mode) >= rank(required)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::InitialPhaseSet;
    use std::sync::mpsc;

    #[test]
    fn disabled_collector_keeps_no_events() {
        let mut collector =
            EquilibriumDiagnosticsCollector::new(EquilibriumDiagnosticsOptions::disabled());
        collector.record(
            EquilibriumDiagnosticsMode::Summary,
            EquilibriumDiagnosticEvent::SolveAccepted {
                final_phase_set: PhaseSet::from_policy(
                    &InitialPhaseSet::AllCandidatePhases,
                    &[true],
                )
                .unwrap(),
                outer_iterations: 1,
                transition_count: 0,
            },
        );
        assert!(!collector.finish().enabled());
        assert!(collector.finish().events().is_empty());
    }

    #[test]
    fn bounded_report_drops_only_retained_events() {
        let mut collector = EquilibriumDiagnosticsCollector::new(
            EquilibriumDiagnosticsOptions::enabled(EquilibriumDiagnosticsMode::Summary)
                .with_max_events(1),
        );
        let phase_set =
            PhaseSet::from_policy(&InitialPhaseSet::AllCandidatePhases, &[true]).unwrap();
        for _ in 0..2 {
            collector.record(
                EquilibriumDiagnosticsMode::Summary,
                EquilibriumDiagnosticEvent::SolveAccepted {
                    final_phase_set: phase_set.clone(),
                    outer_iterations: 1,
                    transition_count: 0,
                },
            );
        }
        let report = collector.finish();
        assert_eq!(report.events().len(), 1);
        assert_eq!(report.dropped_events(), 1);
    }

    #[test]
    fn enabled_sink_receives_events_even_when_retention_is_disabled() {
        let (sender, receiver) = mpsc::channel();
        let options = EquilibriumDiagnosticsOptions::enabled(EquilibriumDiagnosticsMode::Summary)
            .with_max_events(0)
            .with_sink(move |event| sender.send(event).unwrap());
        let mut collector = EquilibriumDiagnosticsCollector::new(options);
        let phase_set =
            PhaseSet::from_policy(&InitialPhaseSet::AllCandidatePhases, &[true]).unwrap();
        collector.record(
            EquilibriumDiagnosticsMode::Summary,
            EquilibriumDiagnosticEvent::SolveAccepted {
                final_phase_set: phase_set,
                outer_iterations: 1,
                transition_count: 0,
            },
        );
        assert!(matches!(
            receiver.recv().unwrap(),
            EquilibriumDiagnosticEvent::SolveAccepted { .. }
        ));
        assert!(collector.finish().events().is_empty());
    }

    #[test]
    fn endpoint_range_policy_suppresses_middle_point_diagnostics() {
        let options = EquilibriumDiagnosticsOptions::enabled(EquilibriumDiagnosticsMode::Summary);
        assert_ne!(
            options.for_range_point(0, 5).mode(),
            EquilibriumDiagnosticsMode::Disabled
        );
        assert_eq!(
            options.for_range_point(2, 5).mode(),
            EquilibriumDiagnosticsMode::Disabled
        );
        assert_ne!(
            options.for_range_point(4, 5).mode(),
            EquilibriumDiagnosticsMode::Disabled
        );
    }

    #[test]
    fn transitions_only_policy_defers_then_replays_accepted_point_evidence() {
        let (sender, receiver) = mpsc::channel();
        let options = EquilibriumDiagnosticsOptions::enabled(EquilibriumDiagnosticsMode::Summary)
            .with_range_policy(EquilibriumRangeDiagnosticsPolicy::TransitionsOnly)
            .with_sink(move |event| sender.send(event).unwrap());
        let point_options = options.for_range_point(4, 10);
        let mut collector = EquilibriumDiagnosticsCollector::new(point_options.clone());
        let phase_set =
            PhaseSet::from_policy(&InitialPhaseSet::AllCandidatePhases, &[true]).unwrap();
        collector.record(
            EquilibriumDiagnosticsMode::Summary,
            EquilibriumDiagnosticEvent::SolveAccepted {
                final_phase_set: phase_set,
                outer_iterations: 1,
                transition_count: 1,
            },
        );
        let report = collector.finish();
        assert!(receiver.try_recv().is_err());

        point_options.replay_retained_events(&report);
        assert!(matches!(
            receiver.recv().unwrap(),
            EquilibriumDiagnosticEvent::SolveAccepted {
                transition_count: 1,
                ..
            }
        ));
    }
}
