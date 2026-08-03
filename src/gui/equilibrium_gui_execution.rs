//! Transactional lifecycle gate for chemical-equilibrium GUI workers.
//!
//! The gate is intentionally independent of egui and thread creation. A
//! caller can run the production request on any executor, then hand its
//! completion back here. Only the currently active `(run_id, fingerprint)`
//! pair may publish a result.

use std::marker::PhantomData;

/// Monotonic identity attached to one worker request.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EquilibriumGuiRunTicket {
    run_id: u64,
    fingerprint: u64,
}

impl EquilibriumGuiRunTicket {
    pub fn run_id(self) -> u64 {
        self.run_id
    }

    pub fn fingerprint(self) -> u64 {
        self.fingerprint
    }
}

/// Visible lifecycle state of one equilibrium editor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumGuiRunState {
    Idle,
    Resolving { run_id: u64 },
    Solving { run_id: u64 },
    Cancelling { run_id: u64 },
    Completed { run_id: u64 },
    Failed { run_id: u64 },
}

impl Default for EquilibriumGuiRunState {
    fn default() -> Self {
        Self::Idle
    }
}

/// Completion message delivered by a repository/solver worker.
#[derive(Debug, Clone, PartialEq)]
pub enum EquilibriumGuiWorkerMessage<R> {
    Completed {
        ticket: EquilibriumGuiRunTicket,
        result: R,
    },
    Failed {
        ticket: EquilibriumGuiRunTicket,
        error: String,
    },
}

/// Result of attempting to publish a worker message.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumGuiPublication {
    Accepted,
    Stale,
}

/// Generic lifecycle owner. The result type is generic so the gate can be
/// tested with a tiny fake and later reused for the real immutable equilibrium
/// snapshot without coupling lifecycle tests to database fixtures.
#[derive(Debug)]
pub struct EquilibriumGuiExecution<R> {
    next_run_id: u64,
    active: Option<EquilibriumGuiRunTicket>,
    state: EquilibriumGuiRunState,
    accepted_result: Option<R>,
    accepted_fingerprint: Option<u64>,
    last_error: Option<String>,
    _result_marker: PhantomData<fn() -> R>,
}

impl<R> Default for EquilibriumGuiExecution<R> {
    fn default() -> Self {
        Self {
            next_run_id: 0,
            active: None,
            state: EquilibriumGuiRunState::Idle,
            accepted_result: None,
            accepted_fingerprint: None,
            last_error: None,
            _result_marker: PhantomData,
        }
    }
}

impl<R> EquilibriumGuiExecution<R> {
    /// Clears every accepted and in-flight result while retaining the
    /// monotonic run-id source. Document loading uses this as a hard lifecycle
    /// boundary: no result from the previous document may be reused.
    pub fn reset(&mut self) {
        self.active = None;
        self.state = EquilibriumGuiRunState::Idle;
        self.accepted_result = None;
        self.accepted_fingerprint = None;
        self.last_error = None;
    }

    /// Starts a new resolve/solve transaction. Any previous worker becomes
    /// stale immediately, while the last accepted result remains available as
    /// an explicitly old snapshot until a new result is accepted.
    pub fn begin(&mut self, fingerprint: u64) -> EquilibriumGuiRunTicket {
        self.next_run_id = self.next_run_id.wrapping_add(1).max(1);
        let ticket = EquilibriumGuiRunTicket {
            run_id: self.next_run_id,
            fingerprint,
        };
        self.active = Some(ticket);
        self.state = EquilibriumGuiRunState::Resolving {
            run_id: ticket.run_id,
        };
        self.last_error = None;
        ticket
    }

    /// Marks the active ticket as solver-running.
    pub fn mark_solving(&mut self, ticket: EquilibriumGuiRunTicket) -> bool {
        if self.active == Some(ticket) {
            self.state = EquilibriumGuiRunState::Solving {
                run_id: ticket.run_id,
            };
            true
        } else {
            false
        }
    }

    /// Requests cancellation. Completion from this ticket is stale even if a
    /// backend finishes later because the current engine has no cooperative
    /// cancellation token yet.
    pub fn cancel(&mut self, ticket: EquilibriumGuiRunTicket) -> bool {
        if self.active == Some(ticket) {
            self.active = None;
            self.state = EquilibriumGuiRunState::Cancelling {
                run_id: ticket.run_id,
            };
            true
        } else {
            false
        }
    }

    /// Invalidates a worker ticket without pretending that the worker failed.
    /// This is used when the owning document changes while a worker is still
    /// running; its eventual message must be stale even if it reports success.
    pub fn discard(&mut self, ticket: EquilibriumGuiRunTicket) -> bool {
        self.cancel(ticket)
    }

    /// Applies a worker message only when its ticket is still current.
    pub fn publish(
        &mut self,
        message: EquilibriumGuiWorkerMessage<R>,
    ) -> EquilibriumGuiPublication {
        let ticket = match &message {
            EquilibriumGuiWorkerMessage::Completed { ticket, .. }
            | EquilibriumGuiWorkerMessage::Failed { ticket, .. } => *ticket,
        };
        if self.active != Some(ticket) {
            return EquilibriumGuiPublication::Stale;
        }

        self.active = None;
        match message {
            EquilibriumGuiWorkerMessage::Completed { result, .. } => {
                self.accepted_result = Some(result);
                self.accepted_fingerprint = Some(ticket.fingerprint);
                self.last_error = None;
                self.state = EquilibriumGuiRunState::Completed {
                    run_id: ticket.run_id,
                };
            }
            EquilibriumGuiWorkerMessage::Failed { error, .. } => {
                self.last_error = Some(error);
                self.state = EquilibriumGuiRunState::Failed {
                    run_id: ticket.run_id,
                };
            }
        }
        EquilibriumGuiPublication::Accepted
    }

    pub fn state(&self) -> EquilibriumGuiRunState {
        self.state
    }

    pub fn active_ticket(&self) -> Option<EquilibriumGuiRunTicket> {
        self.active
    }

    /// Last successful immutable result, even while a newer run is pending.
    pub fn accepted_result(&self) -> Option<&R> {
        self.accepted_result.as_ref()
    }

    pub fn accepted_fingerprint(&self) -> Option<u64> {
        self.accepted_fingerprint
    }

    pub fn result_matches(&self, fingerprint: u64) -> bool {
        self.accepted_fingerprint == Some(fingerprint)
    }

    pub fn last_error(&self) -> Option<&str> {
        self.last_error.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newer_editor_run_makes_old_completion_stale() {
        let mut execution = EquilibriumGuiExecution::<String>::default();
        let first = execution.begin(10);
        execution.mark_solving(first);
        let second = execution.begin(20);
        assert_eq!(
            execution.state(),
            EquilibriumGuiRunState::Resolving {
                run_id: second.run_id()
            }
        );

        let stale = execution.publish(EquilibriumGuiWorkerMessage::Completed {
            ticket: first,
            result: "old".into(),
        });
        assert_eq!(stale, EquilibriumGuiPublication::Stale);
        assert!(execution.accepted_result().is_none());
        assert_eq!(execution.active_ticket(), Some(second));

        let accepted = execution.publish(EquilibriumGuiWorkerMessage::Completed {
            ticket: second,
            result: "new".into(),
        });
        assert_eq!(accepted, EquilibriumGuiPublication::Accepted);
        assert_eq!(execution.accepted_result(), Some(&"new".to_string()));
        assert!(execution.result_matches(20));
    }

    #[test]
    fn cancellation_discards_late_completion() {
        let mut execution = EquilibriumGuiExecution::<u32>::default();
        let ticket = execution.begin(7);
        assert!(execution.mark_solving(ticket));
        assert!(execution.cancel(ticket));
        assert_eq!(
            execution.state(),
            EquilibriumGuiRunState::Cancelling {
                run_id: ticket.run_id()
            }
        );
        assert_eq!(
            execution.publish(EquilibriumGuiWorkerMessage::Completed { ticket, result: 42 }),
            EquilibriumGuiPublication::Stale
        );
        assert!(execution.accepted_result().is_none());
    }

    #[test]
    fn failed_run_keeps_previous_success_but_marks_new_error() {
        let mut execution = EquilibriumGuiExecution::<&'static str>::default();
        let first = execution.begin(1);
        execution.publish(EquilibriumGuiWorkerMessage::Completed {
            ticket: first,
            result: "accepted",
        });
        let second = execution.begin(2);
        assert_eq!(
            execution.publish(EquilibriumGuiWorkerMessage::Failed {
                ticket: second,
                error: "backend failed".into(),
            }),
            EquilibriumGuiPublication::Accepted
        );
        assert_eq!(execution.accepted_result(), Some(&"accepted"));
        assert!(!execution.result_matches(2));
        assert_eq!(execution.last_error(), Some("backend failed"));
        assert_eq!(
            execution.state(),
            EquilibriumGuiRunState::Failed {
                run_id: second.run_id()
            }
        );
    }

    #[test]
    fn reset_discards_accepted_result_and_returns_to_idle() {
        let mut execution = EquilibriumGuiExecution::<String>::default();
        let ticket = execution.begin(42);
        execution.mark_solving(ticket);
        assert_eq!(
            execution.publish(EquilibriumGuiWorkerMessage::Completed {
                ticket,
                result: "accepted".into(),
            }),
            EquilibriumGuiPublication::Accepted
        );
        execution.reset();
        assert_eq!(execution.state(), EquilibriumGuiRunState::Idle);
        assert!(execution.active_ticket().is_none());
        assert!(execution.accepted_result().is_none());
        assert!(execution.last_error().is_none());
    }
}
