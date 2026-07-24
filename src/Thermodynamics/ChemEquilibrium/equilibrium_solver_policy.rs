//! Backend selection policy, cascade budget, and diagnostic solve reports.
//!
//! # Purpose
//!
//! This module defines **how** the equilibrium solver selects and sequences
//! numerical backends. The policy is deliberately independent of residual
//! construction and solver implementations. It records what was attempted so
//! that fallback behavior is a **testable contract** rather than a sequence
//! visible only in logs.
//!
//! # Key Concepts
//!
//! ## Solver Backend
//!
//! A [`SolverBackend`] is one concrete numerical method that can solve the
//! equilibrium system. Currently supported:
//!
//! - **Legacy backends**: [`Solvers::LM`](crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers),
//!   [`Solvers::NR`](crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers),
//!   [`Solvers::TR`](crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers)
//!   — hand-written implementations in [`equilibrium_nonlinear`](super::equilibrium_nonlinear).
//! - **RustedSciThe backends**: symbolic computation engine with automatic
//!   differentiation via [`equilibrium_rst_backend`](super::equilibrium_rst_backend).
//!
//! ## Solver Policy
//!
//! [`SolverPolicy`] determines the execution strategy:
//!
//! - `Single(backend)` — run exactly one backend, no fallback.
//! - `Cascade(vec![backends])` — run backends in declared order, skipping
//!   duplicates, until one succeeds or the budget is exhausted.
//!
//! ## Cascade Budget
//!
//! [`SolverCascadeBudget`] imposes deterministic resource limits:
//!
//! - `max_attempts` — maximum number of backends that may start.
//! - `max_iterations_per_attempt` — nonlinear iterations per backend.
//! - `max_total_iterations` — total iterations across all backends.
//!
//! ## Solve Report
//!
//! [`EquilibriumSolveReport`] records the full cascade execution trace:
//! which backends were attempted, their outcomes, and which (if any) was
//! accepted. This enables deterministic testing of fallback logic.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`SolverBackend`] | Enum: Legacy(LM/NR/TR) or RustedSciThe |
//! | [`SolverPolicy`] | Enum: Single(backend) or Cascade(vec) |
//! | [`SolverCascadeBudget`] | Resource limits for cascade execution |
//! | [`EquilibriumSolveReport`] | Full cascade execution trace |
//! | [`SolverAttemptReport`] | One backend attempt record |
//! | [`SolverAttemptOutcome`] | Enum: Skipped, Accepted, Rejected, Failed, BudgetExhausted |
//! | [`SolverAttemptMetrics`] | Iteration count and residual history |
//!
//! # Dataflow
//!
//! ```text
//!   EquilibriumLogMores::solver_impl()
//!     │
//!     ├── Resolve SolverPolicy into ordered list of backends
//!     ├── Create SolverCascadeBudget from settings
//!     │
//!     └── solve_backend_cascade(backends, budget, request)
//!           │
//!           ├── For each backend in order:
//!           │     ├── Check budget (max_attempts, max_total_iterations)
//!           │     ├── Run backend.solve(request)
//!           │     ├── Record SolverAttemptReport
//!           │     ├── If accepted → return candidate
//!           │     └── If failed → check retryability
//!           │
//!           └── Return EquilibriumSolveReport with all attempts
//! ```
//!
//! # Examples
//!
//! ```rust, ignore
//! use equilibrium_solver_policy::*;
//!
//! // Try RST symbolic first, fall back to legacy LM
//! let policy = SolverPolicy::Cascade(vec![
//!     SolverBackend::RustedSciThe(RustedSciTheSolver::LM),
//!     SolverBackend::Legacy(Solvers::LM),
//! ]);
//!
//! let budget = SolverCascadeBudget {
//!     max_attempts: 2,
//!     max_iterations_per_attempt: 100,
//!     max_total_iterations: 200,
//! };
//! ```
//!
//! # Non-obvious Details
//!
//! - The cascade **skips duplicate backends** — if the same backend appears
//!   twice in the list, the second occurrence is recorded as `Skipped`.
//! - `SolverAttemptOutcome::BudgetExhausted` is recorded when the global
//!   iteration budget runs out before a backend completes, ensuring the
//!   report accurately reflects why later backends were not tried.
//! - `SolverAttemptFailureKind` classifies failures as retryable (numerical
//!   issues) or non-retryable (input errors), controlling cascade continuation.
//!
//! # Related Modules
//!
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — orchestrator that uses the policy
//! - [`equilibrium_backend_adapter`](super::equilibrium_backend_adapter) — adapter dispatching to backends
//! - [`equilibrium_nonlinear`](super::equilibrium_nonlinear) — error types used in outcomes
//!

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
use std::fmt;

/// One concrete numerical backend accepted by an equilibrium solve policy.
///
/// Legacy hand-written methods remain visible only as an explicit fallback
/// family. New policies should prefer the RST symbolic backends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverBackend {
    /// Historical KiThe implementation using a hand-written Jacobian.
    Legacy(Solvers),
    /// RST symbolic backend that owns residual/Jacobian preparation.
    RustedSciThe(RustedSciTheSolver),
}

/// Ordered selection of numerical backends for one equilibrium solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolverPolicy {
    /// Run exactly one backend and do not fall back.
    Single(SolverBackend),
    /// Run backends in the declared order, skipping duplicate entries.
    Cascade(Vec<SolverBackend>),
}

/// Deterministic resource limits for one solver cascade.
///
/// The limits apply to the control plane before a backend is started. They do
/// not claim to measure every internal evaluation of legacy methods, but they
/// guarantee that no attempt receives more iterations than allocated and that
/// later backends are skipped in a reproducible order once the budget is gone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SolverCascadeBudget {
    /// Maximum number of backends that may actually start.
    pub max_attempts: usize,
    /// Maximum nonlinear iterations granted to one started backend.
    pub max_iterations_per_attempt: usize,
    /// Total nonlinear-iteration allocation across all started backends.
    pub max_total_iterations: usize,
}

impl SolverCascadeBudget {
    /// Creates an explicit cascade budget. Validation happens at the solver
    /// boundary so deserialized or programmatic configurations share one path.
    pub const fn new(
        max_attempts: usize,
        max_iterations_per_attempt: usize,
        max_total_iterations: usize,
    ) -> Self {
        Self {
            max_attempts,
            max_iterations_per_attempt,
            max_total_iterations,
        }
    }
}

impl SolverPolicy {
    /// Builds the historical fallback order while making it explicit.
    pub fn legacy_default(preferred: Solvers) -> Self {
        Self::Cascade(vec![
            SolverBackend::Legacy(preferred),
            SolverBackend::Legacy(Solvers::LM),
            SolverBackend::Legacy(Solvers::NR),
            SolverBackend::Legacy(Solvers::TR),
        ])
    }

    /// Builds the preferred RST-only default for new equilibrium workflows.
    pub fn rusted_scithe_default() -> Self {
        Self::Cascade(
            RustedSciTheSolver::recommended_cascade()
                .into_iter()
                .map(SolverBackend::RustedSciThe)
                .collect(),
        )
    }

    /// Returns the deterministic, duplicate-free backend order.
    pub fn ordered_backends(&self) -> Vec<SolverBackend> {
        let requested = match self {
            Self::Single(backend) => std::slice::from_ref(backend),
            Self::Cascade(backends) => backends.as_slice(),
        };

        let mut unique = Vec::with_capacity(requested.len());
        for &backend in requested {
            if !unique.contains(&backend) {
                unique.push(backend);
            }
        }
        unique
    }
}

/// Machine-readable class of a recoverable failed backend attempt.
///
/// This deliberately describes the equilibrium adapter boundary rather than a
/// particular library's private error enum. Callers can therefore react to a
/// failed cascade without parsing diagnostic text.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverAttemptFailureKind {
    /// The nonlinear method terminated without a usable candidate.
    Solver,
    /// A named numerical backend reported a typed internal failure.
    Backend,
    /// Residual evaluation failed during a numerical iteration.
    ResidualEvaluation,
    /// Jacobian evaluation failed during a numerical iteration.
    JacobianEvaluation,
}

/// Observable outcome of one backend attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolverAttemptOutcome {
    /// The backend produced the published candidate.
    Accepted,
    /// The backend returned a numerical failure before producing a candidate.
    Failed {
        /// Typed recoverable failure category.
        kind: SolverAttemptFailureKind,
        /// Stable diagnostic detail for a human or log consumer.
        reason: String,
    },
    /// A candidate was returned but rejected by common acceptance checks.
    RejectedCandidate { reason: String },
    /// The policy did not start this backend because the cascade budget was
    /// exhausted before its turn.
    Skipped { reason: String },
}

impl SolverAttemptOutcome {
    /// Returns `true` when the backend produced a candidate that passed the
    /// common acceptance gate.
    pub fn is_accepted(&self) -> bool {
        matches!(self, Self::Accepted)
    }

    /// Returns `true` when the backend was not started because of budget
    /// exhaustion.
    pub fn is_skipped(&self) -> bool {
        matches!(self, Self::Skipped { .. })
    }

    /// Returns `true` when the backend started and then failed or was rejected.
    pub fn is_started(&self) -> bool {
        !self.is_skipped()
    }

    /// Returns the stable human-readable reason if the outcome carries one.
    pub fn reason(&self) -> Option<&str> {
        match self {
            Self::Accepted => None,
            Self::Failed { reason, .. }
            | Self::RejectedCandidate { reason }
            | Self::Skipped { reason } => Some(reason.as_str()),
        }
    }

    /// Returns the failure class for a recoverable failed backend attempt.
    pub fn failure_kind(&self) -> Option<SolverAttemptFailureKind> {
        match self {
            Self::Failed { kind, .. } => Some(*kind),
            _ => None,
        }
    }
}

/// Backend stopping reason in a solver-independent report.
///
/// RST provides these values directly. Temporary legacy adapters leave
/// [`SolverAttemptReport::metrics`] empty until they expose equivalent data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverTermination {
    /// The backend satisfied its convergence criterion.
    Converged,
    /// The backend consumed its allowed nonlinear iterations.
    MaxIterations,
    /// The trial step became too small to make further progress.
    StepTooSmall,
    /// The method detected stalled progress.
    Stagnation,
    /// The backend exhausted its rejected-step allowance.
    RejectedStepLimit,
}

/// Backend-provided diagnostics for a candidate-producing solve attempt.
///
/// Counters are optional because temporary legacy methods do not expose a
/// comparable engine-level trace. RST adapters fill every field they receive
/// from `SolveResult`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolverAttemptMetrics {
    /// Typed stopping classification reported by the backend.
    pub termination: SolverTermination,
    /// Whether the backend itself reported its canonical convergence state.
    /// A `false` value does not by itself discard a candidate: the common
    /// equilibrium acceptance gate remains the final authority.
    pub backend_converged: bool,
    /// Completed nonlinear iterations.
    pub iterations: usize,
    /// Number of residual evaluations performed by the backend.
    pub residual_evaluations: usize,
    /// Number of Jacobian evaluations performed by the backend.
    pub jacobian_evaluations: usize,
    /// Number of linear subproblems solved by the backend.
    pub linear_solves: usize,
    /// Wall-clock time spent in the backend call.
    pub elapsed_millis: u128,
}

/// One entry in a deterministic solver-cascade trace.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolverAttemptReport {
    /// Numerical backend selected for this attempt.
    pub backend: SolverBackend,
    /// Terminal outcome of this attempt.
    pub outcome: SolverAttemptOutcome,
    /// Numerical counters when the backend exposes them.
    pub metrics: Option<SolverAttemptMetrics>,
}

impl SolverAttemptReport {
    /// Returns `true` when this backend attempt produced the accepted candidate.
    pub fn is_accepted(&self) -> bool {
        self.outcome.is_accepted()
    }

    /// Returns `true` when this backend attempt was skipped by policy.
    pub fn is_skipped(&self) -> bool {
        self.outcome.is_skipped()
    }

    /// Returns `true` when this backend attempt started execution.
    pub fn is_started(&self) -> bool {
        self.outcome.is_started()
    }

    /// Returns the failure kind when the attempt ended in a recoverable failure.
    pub fn failure_kind(&self) -> Option<SolverAttemptFailureKind> {
        self.outcome.failure_kind()
    }

    /// Returns a compact human-readable summary for logging or UI.
    pub fn summary(&self) -> String {
        match &self.outcome {
            SolverAttemptOutcome::Accepted => {
                format!("backend={:?}, outcome=accepted", self.backend)
            }
            SolverAttemptOutcome::Failed { kind, reason } => {
                format!(
                    "backend={:?}, outcome=failed({kind:?}): {reason}",
                    self.backend
                )
            }
            SolverAttemptOutcome::RejectedCandidate { reason } => {
                format!(
                    "backend={:?}, outcome=rejected_candidate: {reason}",
                    self.backend
                )
            }
            SolverAttemptOutcome::Skipped { reason } => {
                format!("backend={:?}, outcome=skipped: {reason}", self.backend)
            }
        }
    }
}

impl fmt::Display for SolverAttemptReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.summary())
    }
}

/// Diagnostics for the backend selection phase of one accepted solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquilibriumSolveReport {
    /// Effective policy used after legacy defaults or explicit overrides.
    pub policy: SolverPolicy,
    /// Attempts in execution order.
    pub attempts: Vec<SolverAttemptReport>,
    /// Backend whose candidate passed acceptance checks.
    pub accepted_backend: SolverBackend,
}

impl EquilibriumSolveReport {
    /// Total number of backend entries recorded in this solve.
    pub fn attempt_count(&self) -> usize {
        self.attempts.len()
    }

    /// Number of started attempts, excluding those skipped by the budget.
    pub fn started_attempt_count(&self) -> usize {
        self.attempts
            .iter()
            .filter(|attempt| !matches!(attempt.outcome, SolverAttemptOutcome::Skipped { .. }))
            .count()
    }

    /// Number of fallback attempts after the first recorded backend entry.
    pub fn fallback_attempt_count(&self) -> usize {
        self.started_attempt_count().saturating_sub(1)
    }

    /// Number of skipped attempts recorded because the budget was exhausted.
    pub fn skipped_attempt_count(&self) -> usize {
        self.attempts
            .iter()
            .filter(|attempt| matches!(attempt.outcome, SolverAttemptOutcome::Skipped { .. }))
            .count()
    }

    /// Returns the index of the accepted backend in the recorded attempt list.
    pub fn accepted_attempt_index(&self) -> Option<usize> {
        self.attempts
            .iter()
            .position(|attempt| attempt.backend == self.accepted_backend)
    }

    /// Returns the full report entry for the accepted backend, if it exists.
    pub fn accepted_attempt(&self) -> Option<&SolverAttemptReport> {
        self.accepted_attempt_index()
            .and_then(|index| self.attempts.get(index))
    }

    /// Returns one attempt by index if it exists.
    pub fn attempt(&self, index: usize) -> Option<&SolverAttemptReport> {
        self.attempts.get(index)
    }

    /// Returns the ordered attempt list as a pure view.
    pub fn attempt_backends(&self) -> Vec<SolverBackend> {
        self.attempts
            .iter()
            .map(|attempt| attempt.backend)
            .collect()
    }

    /// Returns `true` if the solve reached acceptance after at least one fallback.
    pub fn accepted_after_fallback(&self) -> bool {
        self.started_attempt_count() > 1
    }

    /// Returns a compact one-line summary suitable for logs or UI badges.
    pub fn summary(&self) -> String {
        format!(
            "policy={:?}, attempts={}, started={}, skipped={}, fallback={}, accepted_backend={:?}",
            self.policy,
            self.attempt_count(),
            self.started_attempt_count(),
            self.skipped_attempt_count(),
            self.fallback_attempt_count(),
            self.accepted_backend
        )
    }
}

impl fmt::Display for EquilibriumSolveReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.summary())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_policy_keeps_preferred_backend_first_and_removes_duplicates() {
        let policy = SolverPolicy::legacy_default(Solvers::NR);
        assert_eq!(
            policy.ordered_backends(),
            vec![
                SolverBackend::Legacy(Solvers::NR),
                SolverBackend::Legacy(Solvers::LM),
                SolverBackend::Legacy(Solvers::TR),
            ]
        );
    }

    #[test]
    fn single_policy_never_adds_a_fallback_backend() {
        assert_eq!(
            SolverPolicy::Single(SolverBackend::Legacy(Solvers::TR)).ordered_backends(),
            vec![SolverBackend::Legacy(Solvers::TR)]
        );
    }

    #[test]
    fn solve_report_summarizes_attempts_without_manual_vector_traversal() {
        let report = EquilibriumSolveReport {
            policy: SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM)),
            attempts: vec![
                SolverAttemptReport {
                    backend: SolverBackend::Legacy(Solvers::LM),
                    outcome: SolverAttemptOutcome::Failed {
                        kind: SolverAttemptFailureKind::Solver,
                        reason: "first".to_string(),
                    },
                    metrics: None,
                },
                SolverAttemptReport {
                    backend: SolverBackend::Legacy(Solvers::NR),
                    outcome: SolverAttemptOutcome::Accepted,
                    metrics: None,
                },
                SolverAttemptReport {
                    backend: SolverBackend::Legacy(Solvers::TR),
                    outcome: SolverAttemptOutcome::Skipped {
                        reason: "budget exhausted".to_string(),
                    },
                    metrics: None,
                },
            ],
            accepted_backend: SolverBackend::Legacy(Solvers::NR),
        };

        assert_eq!(report.attempt_count(), 3);
        assert_eq!(report.started_attempt_count(), 2);
        assert_eq!(report.fallback_attempt_count(), 1);
        assert_eq!(report.skipped_attempt_count(), 1);
        assert_eq!(report.accepted_attempt_index(), Some(1));
        assert!(report.accepted_attempt().unwrap().is_accepted());
        assert!(report.attempt(2).unwrap().is_skipped());
        assert_eq!(
            report.attempt_backends(),
            vec![
                SolverBackend::Legacy(Solvers::LM),
                SolverBackend::Legacy(Solvers::NR),
                SolverBackend::Legacy(Solvers::TR),
            ]
        );
        assert!(report.accepted_after_fallback());
        assert!(report.summary().contains("attempts=3"));
        assert!(report.to_string().contains("accepted_backend"));
        assert!(report.attempt(0).unwrap().to_string().contains("failed"));
    }

    #[test]
    fn attempt_outcome_exposes_reason_and_failure_kind() {
        let outcome = SolverAttemptOutcome::Failed {
            kind: SolverAttemptFailureKind::ResidualEvaluation,
            reason: "temporary residual failure".to_string(),
        };
        assert!(outcome.is_started());
        assert!(!outcome.is_skipped());
        assert_eq!(
            outcome.failure_kind(),
            Some(SolverAttemptFailureKind::ResidualEvaluation)
        );
        assert_eq!(outcome.reason(), Some("temporary residual failure"));

        let skipped = SolverAttemptOutcome::Skipped {
            reason: "budget exhausted".to_string(),
        };
        assert!(skipped.is_skipped());
        assert_eq!(skipped.reason(), Some("budget exhausted"));
    }
}
