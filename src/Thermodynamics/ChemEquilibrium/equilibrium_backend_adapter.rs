//! Backend adapter trait and uniform request/result types for solver dispatch.
//!
//! # Purpose
//!
//! This module defines the **adapter contract** between the solver policy
//! (which decides which backend to run) and the concrete backend implementations
//! (legacy LM/NR/TR and RustedSciThe). It provides:
//!
//! - [`BackendSolveRequest`] — uniform input to any backend.
//! - [`BackendSolveResult`] — uniform output from any backend.
//! - [`EquilibriumNonlinearBackend`] — trait that all backends implement.
//!
//! The module is deliberately **tiny and mechanical**. It isolates the policy
//! from implementation details and provides a place for deterministic fake
//! backends used in cascade tests.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`BackendSolveRequest`] | Bundles initial guess, residual, Jacobian, feasibility, params |
//! | [`BackendSolveResult`] | Solution vector + optional metrics |
//! | [`EquilibriumNonlinearBackend`] | Trait: `backend()` + `solve()` |
//!
//! # Dataflow
//!
//! ```text
//!   solve_backend_cascade_with_control()
//!     │
//!     ├── Builds BackendSolveRequest from solver state
//!     │     ├── initial_guess: Vec<f64>
//!     │     ├── residual: &dyn Fn(&[f64]) -> Result<Vec<f64>>
//!     │     ├── jacobian: Option<&dyn Fn(&[f64]) -> Result<DMatrix<f64>>>
//!     │     ├── feasible: &dyn Fn(&[f64]) -> bool
//!     │     ├── params: &SolverParams
//!     │     ├── initial_moles: &[f64]
//!     │     ├── reactions: &DMatrix<f64>
//!     │     ├── max_iterations: usize
//!     │     └── rst_problem: Option<&RstPreparedProblem>
//!     │
//!     ├── For each backend in policy:
//!     │     ├── backend.solve(request) -> Result<BackendSolveResult>
//!     │     └── On success: validate candidate
//!     │
//!     └── Return accepted candidate or cascade failure
//! ```
//!
//! # Backend Implementations
//!
//! | Backend | Module | Description |
//! |---------|--------|-------------|
//! | Legacy(LM/NR/TR) | [`equilibrium_legacy_backend`](super::equilibrium_legacy_backend) | Hand-written solvers |
//! | RustedSciThe | [`equilibrium_rst_backend`](super::equilibrium_rst_backend) | Symbolic computation engine |
//! | FakeBackend | (test-only) | Deterministic mock for cascade tests |
//!
//! # Non-obvious Details
//!
//! - The `jacobian` field is `Option` because some backends (RustedSciThe) compute
//!   their own Jacobian and don't need the analytical one.
//! - `BackendSolveResult.metrics` is `Option<SolverAttemptMetrics>` because some
//!   backends may not report iteration counts.
//! - The trait is `pub(crate)` — it is an internal implementation detail, not
//!   part of the public API.
//!
//! # Related Modules
//!
//! - [`equilibrium_solver_policy`](super::equilibrium_solver_policy) — backend selection
//! - [`equilibrium_legacy_backend`](super::equilibrium_legacy_backend) — legacy adapter
//! - [`equilibrium_rst_backend`](super::equilibrium_rst_backend) — RST adapter
//!

use crate::Thermodynamics::ChemEquilibrium::equilibrium_execution::{
    EquilibriumExecutionControl, EquilibriumProgressEvent, EquilibriumProgressStage,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_legacy_backend::solve_legacy_backend;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    recoverable_backend_failure_kind, SolverParams,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, RustedSciTheSolveContract,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverAttemptMetrics, SolverAttemptOutcome, SolverAttemptReport,
    SolverBackend, SolverCascadeBudget, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
use log::info;
use nalgebra::DMatrix;

type ResidualFn<'a> = dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError> + 'a;
type JacobianFn<'a> = dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError> + 'a;
type FeasibilityFn<'a> = dyn Fn(&[f64]) -> bool + 'a;

/// Domain-independent prepared nonlinear-system capabilities.
///
/// The coordinate meaning is intentionally absent. A fixed-P,T problem may
/// use log-moles only; a future P,H problem may append a temperature
/// coordinate. The adapter only needs the dimension and callable residual,
/// Jacobian, and feasibility contracts. The optional RST payload is an
/// implementation detail of the symbolic backend, not a thermodynamic model.
pub(crate) struct PreparedNonlinearSystem<'a> {
    /// Number of unknowns expected by the prepared capabilities.
    pub dimension: usize,
    /// Residual capability for the prepared coordinate vector.
    pub residual: &'a ResidualFn<'a>,
    /// Optional analytic Jacobian capability.
    pub jacobian: Option<&'a JacobianFn<'a>>,
    /// Domain-side feasibility predicate.
    pub feasible: &'a FeasibilityFn<'a>,
    /// Optional symbolic payload consumed only by an RST backend.
    pub rst_problem: Option<&'a RstPreparedProblem>,
}

impl<'a> PreparedNonlinearSystem<'a> {
    pub(crate) fn new(
        dimension: usize,
        residual: &'a ResidualFn<'a>,
        jacobian: Option<&'a JacobianFn<'a>>,
        feasible: &'a FeasibilityFn<'a>,
        rst_problem: Option<&'a RstPreparedProblem>,
    ) -> Result<Self, ReactionExtentError> {
        if dimension == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "nonlinear_system.dimension",
                message: "prepared nonlinear system must contain at least one unknown".to_string(),
            });
        }
        Ok(Self {
            dimension,
            residual,
            jacobian,
            feasible,
            rst_problem,
        })
    }
}

/// Mechanical request passed to one backend attempt.
pub(crate) struct BackendSolveRequest<'a> {
    pub initial_guess: Vec<f64>,
    pub system: PreparedNonlinearSystem<'a>,
    pub params: &'a SolverParams,
    pub max_iterations: usize,
}

/// Uniform output of one backend attempt.
#[derive(Debug, Clone)]
pub(crate) struct BackendSolveResult {
    pub solution: Vec<f64>,
    pub metrics: Option<SolverAttemptMetrics>,
}

/// Internal adapter contract for one numerical backend.
pub(crate) trait EquilibriumNonlinearBackend {
    fn backend(&self) -> SolverBackend;

    fn solve(
        &self,
        request: BackendSolveRequest<'_>,
    ) -> Result<BackendSolveResult, ReactionExtentError>;
}

impl EquilibriumNonlinearBackend for SolverBackend {
    fn backend(&self) -> SolverBackend {
        *self
    }

    fn solve(
        &self,
        request: BackendSolveRequest<'_>,
    ) -> Result<BackendSolveResult, ReactionExtentError> {
        if request.initial_guess.len() != request.system.dimension {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "nonlinear initial guess has {} entries, prepared system requires {}",
                request.initial_guess.len(),
                request.system.dimension
            )));
        }
        match *self {
            SolverBackend::Legacy(legacy_backend) => solve_legacy_backend(
                legacy_backend,
                request.initial_guess,
                request.system.residual,
                request.system.jacobian,
                request.system.feasible,
                request.params,
                request.max_iterations,
            )
            .map(|solution| BackendSolveResult {
                solution,
                metrics: None,
            }),
            SolverBackend::RustedSciThe(rst_solver) => match request.system.rst_problem {
                Some(problem) => rst_solver
                    .solve(
                        problem,
                        &request.initial_guess,
                        &RustedSciTheSolveContract::new(
                            request.params.tol,
                            request.max_iterations,
                        )?,
                    )
                    .map(|outcome| BackendSolveResult {
                        solution: outcome.solution,
                        metrics: Some(outcome.metrics),
                    }),
                None => Err(ReactionExtentError::InvalidProblem {
                    field: "rst_symbolic_problem",
                    message: "RST backend selected without a prepared symbolic problem".to_string(),
                }),
            },
        }
    }
}

pub(crate) fn solve_backend_cascade_with_control(
    backends: &[&dyn EquilibriumNonlinearBackend],
    initial_guess: Vec<f64>,
    f: &dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
    j: Option<&dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>,
    feasible: &dyn Fn(&[f64]) -> bool,
    validate_candidate: &dyn Fn(&[f64]) -> Result<EquilibriumCandidateReport, ReactionExtentError>,
    policy: SolverPolicy,
    budget: SolverCascadeBudget,
    solver_params: &SolverParams,
    rst_problem: Option<&RstPreparedProblem>,
    execution_control: Option<&EquilibriumExecutionControl>,
) -> Result<(Vec<f64>, EquilibriumCandidateReport, EquilibriumSolveReport), ReactionExtentError> {
    if backends.is_empty() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "solver_policy",
            message: "a solver policy must contain at least one backend".to_string(),
        });
    }

    if budget.max_attempts == 0
        || budget.max_iterations_per_attempt == 0
        || budget.max_total_iterations == 0
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "solver_budget",
            message: "attempt, per-attempt, and total iteration limits must be positive"
                .to_string(),
        });
    }

    let system = PreparedNonlinearSystem::new(initial_guess.len(), f, j, feasible, rst_problem)?;

    let mut attempts = Vec::with_capacity(backends.len());
    let mut started_attempts = 0usize;
    let mut remaining_iterations = budget.max_total_iterations;
    for &backend in backends {
        let backend_id = backend.backend();
        if started_attempts >= budget.max_attempts || remaining_iterations == 0 {
            let reason = if started_attempts >= budget.max_attempts {
                format!(
                    "cascade attempt budget exhausted after {} started backend(s)",
                    budget.max_attempts
                )
            } else {
                "cascade total iteration budget exhausted".to_string()
            };
            attempts.push(SolverAttemptReport {
                backend: backend_id,
                outcome: SolverAttemptOutcome::Skipped { reason },
                metrics: None,
            });
            continue;
        }
        let max_iterations = budget.max_iterations_per_attempt.min(remaining_iterations);
        started_attempts += 1;
        remaining_iterations -= max_iterations;

        // Every backend receives this same validated seed. A previous
        // rejected candidate is evidence for diagnostics, never an
        // implicit warm start.
        if let Some(control) = execution_control {
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::InnerBackendAttemptStarted,
                None,
                None,
                None,
            ));
            control.check_cancelled()?;
        }
        let result = backend.solve(BackendSolveRequest {
            initial_guess: initial_guess.clone(),
            system: PreparedNonlinearSystem {
                dimension: system.dimension,
                residual: system.residual,
                jacobian: system.jacobian,
                feasible: system.feasible,
                rst_problem: system.rst_problem,
            },
            params: solver_params,
            max_iterations,
        });
        if let Some(control) = execution_control {
            control.report(EquilibriumProgressEvent::new(
                EquilibriumProgressStage::InnerBackendAttemptFinished,
                None,
                None,
                None,
            ));
            control.check_cancelled()?;
        }

        match result {
            Ok(candidate) => match validate_candidate(&candidate.solution) {
                Ok(validation) => {
                    attempts.push(SolverAttemptReport {
                        backend: backend_id,
                        outcome: SolverAttemptOutcome::Accepted,
                        metrics: candidate.metrics,
                    });
                    let report = EquilibriumSolveReport {
                        policy,
                        attempts,
                        accepted_backend: backend_id,
                    };
                    if report.accepted_after_fallback() {
                        info!("Solver accepted after fallback: {}", report);
                    }
                    return Ok((candidate.solution, validation, report));
                }
                Err(error) => attempts.push(SolverAttemptReport {
                    backend: backend_id,
                    outcome: SolverAttemptOutcome::RejectedCandidate {
                        reason: error.to_string(),
                    },
                    metrics: candidate.metrics,
                }),
            },
            Err(e) => {
                let Some(kind) = recoverable_backend_failure_kind(&e) else {
                    return if attempts.is_empty() {
                        Err(e)
                    } else {
                        Err(ReactionExtentError::CascadeAborted {
                            attempts,
                            cause: Box::new(e),
                        })
                    };
                };
                attempts.push(SolverAttemptReport {
                    backend: backend_id,
                    outcome: SolverAttemptOutcome::Failed {
                        kind,
                        // Attempt reports are part of the public solve
                        // contract, so keep their wording independent of
                        // a private Debug implementation.
                        reason: e.to_string(),
                    },
                    metrics: None,
                });
            }
        }
    }

    Err(ReactionExtentError::AllBackendsFailed { attempts })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::BackendFailureKind;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverBackend;

    #[derive(Clone, Copy)]
    enum FakeBackendBehavior {
        Succeed,
        Fail,
    }

    struct FakeBackend {
        backend: SolverBackend,
        behavior: FakeBackendBehavior,
    }

    impl EquilibriumNonlinearBackend for FakeBackend {
        fn backend(&self) -> SolverBackend {
            self.backend
        }

        fn solve(
            &self,
            request: BackendSolveRequest<'_>,
        ) -> Result<BackendSolveResult, ReactionExtentError> {
            match self.behavior {
                FakeBackendBehavior::Succeed => Ok(BackendSolveResult {
                    solution: request.initial_guess,
                    metrics: Some(SolverAttemptMetrics {
                        termination: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverTermination::Converged,
                        backend_converged: true,
                        iterations: 1,
                        residual_evaluations: 1,
                        jacobian_evaluations: 1,
                        linear_solves: 1,
                        elapsed_millis: 0,
                        evaluation_timing: None,
                    }),
                }),
                FakeBackendBehavior::Fail => Err(ReactionExtentError::BackendFailure {
                    backend: "fake".to_string(),
                    kind: BackendFailureKind::NumericalBreakdown,
                    message: "forced fake failure".to_string(),
                }),
            }
        }
    }

    #[test]
    fn fake_backend_can_be_used_as_a_deterministic_test_double() {
        let residual = |values: &[f64]| Ok(values.to_vec());
        let feasible = |_values: &[f64]| true;
        let request = BackendSolveRequest {
            initial_guess: vec![1.0],
            system: PreparedNonlinearSystem::new(1, &residual, None, &feasible, None).unwrap(),
            params: &SolverParams::default(),
            max_iterations: 1,
        };

        let backend = FakeBackend {
            backend: SolverBackend::Legacy(Solvers::LM),
            behavior: FakeBackendBehavior::Succeed,
        };
        let result = backend.solve(request).unwrap();
        assert_eq!(result.solution, vec![1.0]);
        assert!(result.metrics.as_ref().unwrap().backend_converged);
    }

    #[test]
    fn fake_backend_can_force_a_repeatable_failure() {
        let residual = |values: &[f64]| Ok(values.to_vec());
        let feasible = |_values: &[f64]| true;
        let request = BackendSolveRequest {
            initial_guess: vec![1.0],
            system: PreparedNonlinearSystem::new(1, &residual, None, &feasible, None).unwrap(),
            params: &SolverParams::default(),
            max_iterations: 1,
        };

        let backend = FakeBackend {
            backend: SolverBackend::Legacy(Solvers::LM),
            behavior: FakeBackendBehavior::Fail,
        };
        assert!(matches!(
            backend.solve(request),
            Err(ReactionExtentError::BackendFailure { .. })
        ));
    }

    #[test]
    fn prepared_nonlinear_system_rejects_an_empty_coordinate_contract() {
        let residual = |_values: &[f64]| Ok(Vec::new());
        let feasible = |_values: &[f64]| true;

        assert!(matches!(
            PreparedNonlinearSystem::new(0, &residual, None, &feasible, None),
            Err(ReactionExtentError::InvalidProblem {
                field: "nonlinear_system.dimension",
                ..
            })
        ));
    }

    #[test]
    fn backend_rejects_a_seed_with_the_wrong_prepared_dimension() {
        let residual = |values: &[f64]| Ok(values.to_vec());
        let feasible = |_values: &[f64]| true;
        let request = BackendSolveRequest {
            initial_guess: vec![1.0],
            system: PreparedNonlinearSystem::new(2, &residual, None, &feasible, None).unwrap(),
            params: &SolverParams::default(),
            max_iterations: 1,
        };
        let backend = SolverBackend::Legacy(Solvers::LM);

        assert!(matches!(
            backend.solve(request),
            Err(ReactionExtentError::DimensionMismatch(message))
                if message.contains("prepared system requires 2")
        ));
    }
}
