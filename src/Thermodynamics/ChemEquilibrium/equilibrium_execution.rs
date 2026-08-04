//! Cooperative execution control for typed equilibrium workflows.
//!
//! The control handle is deliberately small and cloneable.  The solver checks
//! it at safe transaction boundaries and reports progress without owning GUI
//! state, so a CLI, service, or egui frontend can use the same contract.

use std::fmt;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;

/// Coarse-grained stage reported by a typed equilibrium execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumProgressStage {
    RepositoryLookup,
    FormulationPreparation,
    PointStarted,
    PointAccepted,
    /// One scalar temperature trial of the outer fixed-`P,H` solve started.
    TemperatureTrialStarted,
    /// One scalar temperature trial produced an accepted inner candidate.
    TemperatureTrialAccepted,
    /// One scalar temperature trial could not produce an accepted candidate.
    /// The typed solve error carries the detailed cause and location.
    TemperatureTrialRejected,
    /// The outer P,H workflow is about to start its inner fixed-P,T solve.
    InnerSolveStarted,
    /// The inner fixed-P,T solve produced a candidate for outer enthalpy
    /// evaluation. The candidate is not yet a published P,H result.
    InnerSolveCompleted,
    /// One phase-control transition was accepted inside an inner solve.
    PhaseTransitionAccepted,
    /// The outer P,H bracket has produced a final candidate and publication
    /// validation is beginning.
    PublicationStarted,
    /// The immutable P,H result has passed final validation and is returned.
    PublicationCompleted,
    /// One numerical backend attempt is about to start.
    InnerBackendAttemptStarted,
    /// One numerical backend attempt has returned and its candidate/error is
    /// about to be classified by the cascade.
    InnerBackendAttemptFinished,
}

/// Immutable progress event suitable for a UI status line or a structured log.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumProgressEvent {
    stage: EquilibriumProgressStage,
    point_index: Option<usize>,
    point_count: Option<usize>,
    temperature: Option<f64>,
}

impl EquilibriumProgressEvent {
    pub fn new(
        stage: EquilibriumProgressStage,
        point_index: Option<usize>,
        point_count: Option<usize>,
        temperature: Option<f64>,
    ) -> Self {
        Self {
            stage,
            point_index,
            point_count,
            temperature,
        }
    }

    pub fn stage(self) -> EquilibriumProgressStage {
        self.stage
    }
    pub fn point_index(self) -> Option<usize> {
        self.point_index
    }
    pub fn point_count(self) -> Option<usize> {
        self.point_count
    }
    pub fn temperature(self) -> Option<f64> {
        self.temperature
    }
}

type ProgressSink = Arc<dyn Fn(EquilibriumProgressEvent) + Send + Sync + 'static>;

/// Shared cancellation and progress handle for one solve transaction.
#[derive(Clone, Default)]
pub struct EquilibriumExecutionControl {
    cancelled: Arc<AtomicBool>,
    progress_sink: Option<ProgressSink>,
}

impl fmt::Debug for EquilibriumExecutionControl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EquilibriumExecutionControl")
            .field("cancelled", &self.is_cancel_requested())
            .field("has_progress_sink", &self.progress_sink.is_some())
            .finish()
    }
}

impl EquilibriumExecutionControl {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_progress_sink<F>(mut self, sink: F) -> Self
    where
        F: Fn(EquilibriumProgressEvent) + Send + Sync + 'static,
    {
        self.progress_sink = Some(Arc::new(sink));
        self
    }

    pub fn request_cancel(&self) {
        self.cancelled.store(true, Ordering::Release);
    }

    pub fn is_cancel_requested(&self) -> bool {
        self.cancelled.load(Ordering::Acquire)
    }

    pub(crate) fn check_cancelled(&self) -> Result<(), ReactionExtentError> {
        if self.is_cancel_requested() {
            Err(ReactionExtentError::Cancelled)
        } else {
            Ok(())
        }
    }

    pub(crate) fn report(&self, event: EquilibriumProgressEvent) {
        if let Some(sink) = &self.progress_sink {
            sink(event);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;

    #[test]
    fn cloned_controls_share_cancellation() {
        let first = EquilibriumExecutionControl::new();
        let second = first.clone();
        second.request_cancel();
        assert!(first.is_cancel_requested());
        assert!(matches!(
            first.check_cancelled(),
            Err(ReactionExtentError::Cancelled)
        ));
    }

    #[test]
    fn progress_sink_receives_typed_events() {
        let (sender, receiver) = mpsc::channel();
        let control = EquilibriumExecutionControl::new().with_progress_sink(move |event| {
            sender.send(event).expect("progress receiver remains alive");
        });
        control.report(EquilibriumProgressEvent::new(
            EquilibriumProgressStage::PointAccepted,
            Some(2),
            Some(5),
            Some(800.0),
        ));
        let event = receiver.recv().expect("event must be delivered");
        assert_eq!(event.stage(), EquilibriumProgressStage::PointAccepted);
        assert_eq!(event.point_index(), Some(2));
        assert_eq!(event.point_count(), Some(5));
        assert_eq!(event.temperature(), Some(800.0));
    }

    #[test]
    fn rejected_temperature_trial_is_a_distinct_progress_stage() {
        let event = EquilibriumProgressEvent::new(
            EquilibriumProgressStage::TemperatureTrialRejected,
            Some(4),
            Some(9),
            Some(1_250.0),
        );

        assert_eq!(
            event.stage(),
            EquilibriumProgressStage::TemperatureTrialRejected
        );
        assert_eq!(event.point_index(), Some(4));
        assert_eq!(event.point_count(), Some(9));
        assert_eq!(event.temperature(), Some(1_250.0));
    }
}
