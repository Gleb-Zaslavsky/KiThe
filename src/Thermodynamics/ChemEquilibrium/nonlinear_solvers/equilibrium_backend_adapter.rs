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
//! `BackendSolveError` carries the typed backend error together with any
//! metrics collected before the backend stopped.
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
use crate::Thermodynamics::ChemEquilibrium::equilibrium_instrumentation::{
    LegacySolverDiagnosticsSnapshot, begin_legacy_solver_diagnostics,
    finish_legacy_solver_diagnostics,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_legacy_backend::solve_legacy_backend;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    SolverParams, recoverable_backend_failure_kind,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, RustedSciTheSolveContract,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverAttemptMetrics, SolverAttemptOutcome, SolverAttemptReport,
    SolverBackend, SolverCascadeBudget, SolverEvaluationTiming, SolverPolicy, SolverTermination,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
use RustedSciThe::numerical::SolverTelemetryMode;
use log::info;
use nalgebra::DMatrix;
use std::time::Instant;

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
    pub rst_telemetry_mode: SolverTelemetryMode,
}

/// Uniform output of one backend attempt.
#[derive(Debug, Clone)]
pub(crate) struct BackendSolveResult {
    pub solution: Vec<f64>,
    pub metrics: Option<SolverAttemptMetrics>,
}

/// Backend error plus any diagnostics collected before the backend stopped.
/// Keeping these together prevents failed Legacy attempts from becoming a
/// wall-only row in the cascade ledger.
#[derive(Debug)]
pub(crate) struct BackendSolveError {
    pub error: ReactionExtentError,
    pub metrics: Option<SolverAttemptMetrics>,
}

/// Internal adapter contract for one numerical backend.
pub(crate) trait EquilibriumNonlinearBackend {
    fn backend(&self) -> SolverBackend;

    fn solve(
        &self,
        request: BackendSolveRequest<'_>,
    ) -> Result<BackendSolveResult, BackendSolveError>;
}

fn legacy_solver_metrics(
    diagnostics: LegacySolverDiagnosticsSnapshot,
    elapsed: std::time::Duration,
    termination: SolverTermination,
    backend_converged: bool,
) -> SolverAttemptMetrics {
    let callback_micros = diagnostics
        .residual_micros
        .saturating_add(diagnostics.jacobian_micros);
    let accounted_micros = callback_micros.saturating_add(diagnostics.linear_system_micros);
    SolverAttemptMetrics {
        termination,
        backend_converged,
        iterations: diagnostics.jacobian_calls,
        residual_evaluations: diagnostics.residual_calls,
        jacobian_evaluations: diagnostics.jacobian_calls,
        linear_solves: diagnostics.linear_system_solves,
        counters_collected: true,
        timings_collected: true,
        elapsed_millis: elapsed.as_millis(),
        elapsed_micros: elapsed.as_micros(),
        evaluation_timing: Some(SolverEvaluationTiming {
            residual_evaluation_micros: diagnostics.residual_micros,
            jacobian_evaluation_micros: diagnostics.jacobian_micros,
            legacy_same_iterate_callback_pairs: Some(diagnostics.same_iterate_callback_pairs),
            linear_system_micros: Some(diagnostics.linear_system_micros),
            solver_overhead_micros: elapsed.as_micros().saturating_sub(accounted_micros),
            matrix_operations_micros: diagnostics.matrix_operations_micros,
            norm_micros: diagnostics.norm_micros,
            trial_vector_micros: diagnostics.trial_vector_micros,
            feasibility_micros: diagnostics.feasibility_micros,
            globalization_micros: diagnostics.globalization_micros,
            state_residual_evaluations: None,
            trial_residual_evaluations: None,
            state_jacobian_evaluations: None,
            trial_jacobian_evaluations: None,
            linear_factorization_micros: None,
            linear_system_solve_micros: None,
            linear_factorizations: None,
            accepted_steps: None,
            rejected_steps: None,
            reusable_trial_points: None,
            iteration_telemetry: None,
        }),
    }
}

fn legacy_failure_termination(error: &ReactionExtentError) -> SolverTermination {
    match error {
        ReactionExtentError::SolveError(crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::SolveError::MaxIterations) => {
            SolverTermination::MaxIterations
        }
        _ => SolverTermination::Stagnation,
    }
}

impl EquilibriumNonlinearBackend for SolverBackend {
    fn backend(&self) -> SolverBackend {
        *self
    }

    fn solve(
        &self,
        request: BackendSolveRequest<'_>,
    ) -> Result<BackendSolveResult, BackendSolveError> {
        if request.initial_guess.len() != request.system.dimension {
            return Err(BackendSolveError {
                error: ReactionExtentError::DimensionMismatch(format!(
                    "nonlinear initial guess has {} entries, prepared system requires {}",
                    request.initial_guess.len(),
                    request.system.dimension
                )),
                metrics: None,
            });
        }
        match *self {
            SolverBackend::Legacy(legacy_backend) => {
                let Some(diagnostics_session) = begin_legacy_solver_diagnostics() else {
                    return solve_legacy_backend(
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
                    })
                    .map_err(|error| BackendSolveError {
                        error,
                        metrics: None,
                    });
                };
                let started = Instant::now();
                let solution_result = solve_legacy_backend(
                    legacy_backend,
                    request.initial_guess,
                    request.system.residual,
                    request.system.jacobian,
                    request.system.feasible,
                    request.params,
                    request.max_iterations,
                );
                let diagnostics = finish_legacy_solver_diagnostics(diagnostics_session);
                let elapsed = started.elapsed();
                match solution_result {
                    Ok(solution) => Ok(BackendSolveResult {
                        solution,
                        metrics: Some(legacy_solver_metrics(
                            diagnostics,
                            elapsed,
                            SolverTermination::Converged,
                            true,
                        )),
                    }),
                    Err(error) => Err(BackendSolveError {
                        metrics: Some(legacy_solver_metrics(
                            diagnostics,
                            elapsed,
                            legacy_failure_termination(&error),
                            false,
                        )),
                        error,
                    }),
                }
            }
            SolverBackend::RustedSciThe(rst_solver) => match request.system.rst_problem {
                Some(problem) => rst_solver
                    .solve(
                        problem,
                        &request.initial_guess,
                        &RustedSciTheSolveContract::new(request.params.tol, request.max_iterations)
                            .map(|contract| {
                                contract.with_telemetry_mode(request.rst_telemetry_mode)
                            })
                            .map_err(|error| BackendSolveError {
                                error,
                                metrics: None,
                            })?,
                    )
                    .map(|outcome| BackendSolveResult {
                        solution: outcome.solution,
                        metrics: Some(outcome.metrics),
                    })
                    .map_err(|error| BackendSolveError {
                        error,
                        metrics: None,
                    }),
                None => Err(BackendSolveError {
                    error: ReactionExtentError::InvalidProblem {
                        field: "rst_symbolic_problem",
                        message: "RST backend selected without a prepared symbolic problem"
                            .to_string(),
                    },
                    metrics: None,
                }),
            },
        }
    }
}

#[cfg(test)]
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
    solve_backend_cascade_with_control_and_telemetry(
        backends,
        initial_guess,
        f,
        j,
        feasible,
        validate_candidate,
        policy,
        budget,
        solver_params,
        rst_problem,
        execution_control,
        SolverTelemetryMode::Off,
    )
}

pub(crate) fn solve_backend_cascade_with_control_and_telemetry(
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
    rst_telemetry_mode: SolverTelemetryMode,
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
        let backend_started = Instant::now();
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
            rst_telemetry_mode,
        });
        let backend_elapsed_micros = backend_started.elapsed().as_micros();
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
            Err(backend_error) => {
                let BackendSolveError { error: e, metrics } = backend_error;
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
                        elapsed_micros: backend_elapsed_micros,
                    },
                    metrics,
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
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
        RstPreparedProblem, RstSymbolicProblem, RustedSciTheSolver,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverBackend;
    use RustedSciThe::symbolic::symbolic_engine::Expr;

    #[derive(Clone, Copy)]
    enum FakeBackendBehavior {
        Succeed,
        Fail,
        FailAfterDelay,
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
        ) -> Result<BackendSolveResult, BackendSolveError> {
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
                        counters_collected: true,
                        timings_collected: false,
                        elapsed_millis: 0,
                        elapsed_micros: 0,
                        evaluation_timing: None,
                    }),
                }),
                FakeBackendBehavior::Fail => Err(BackendSolveError {
                    error: ReactionExtentError::BackendFailure {
                        backend: "fake".to_string(),
                        kind: BackendFailureKind::NumericalBreakdown,
                        message: "forced fake failure".to_string(),
                    },
                    metrics: None,
                }),
                FakeBackendBehavior::FailAfterDelay => {
                    std::thread::sleep(std::time::Duration::from_millis(2));
                    Err(BackendSolveError {
                        error: ReactionExtentError::BackendFailure {
                            backend: "fake".to_string(),
                            kind: BackendFailureKind::NumericalBreakdown,
                            message: "forced delayed fake failure".to_string(),
                        },
                        metrics: None,
                    })
                }
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
            rst_telemetry_mode: SolverTelemetryMode::Off,
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
    fn calculator_telemetry_mode_reaches_the_rst_backend_attempt() {
        let rst_problem = RstPreparedProblem::new(
            RstSymbolicProblem::from_expressions(
                vec![Expr::parse_expression("x - 1.0")],
                vec!["x".to_string()],
                None,
                None,
            )
            .unwrap(),
        );
        let residual = |values: &[f64]| Ok(vec![values[0] - 1.0]);
        let feasible = |_values: &[f64]| true;
        let reject = |_values: &[f64]| {
            Err(ReactionExtentError::ResidualEvaluation(
                "intentional telemetry probe rejection".to_string(),
            ))
        };
        let backend = SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt);
        let backends: [&dyn EquilibriumNonlinearBackend; 1] = [&backend];
        let error = solve_backend_cascade_with_control_and_telemetry(
            &backends,
            vec![0.0],
            &residual,
            None,
            &feasible,
            &reject,
            SolverPolicy::Single(backend),
            SolverCascadeBudget::new(1, 20, 20),
            &SolverParams::default(),
            Some(&rst_problem),
            None,
            SolverTelemetryMode::Detailed,
        )
        .expect_err("the candidate is deliberately rejected after solving");
        let reports = error.backend_attempt_reports();
        let metrics = reports[0]
            .metrics
            .as_ref()
            .expect("RST attempt telemetry must survive candidate rejection");

        assert!(metrics.counters_collected);
        assert!(metrics.timings_collected);
        assert!(metrics.residual_evaluations > 0);
        assert!(
            metrics
                .evaluation_timing
                .as_ref()
                .unwrap()
                .iteration_telemetry
                .is_some()
        );
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
            rst_telemetry_mode: SolverTelemetryMode::Off,
        };

        let backend = FakeBackend {
            backend: SolverBackend::Legacy(Solvers::LM),
            behavior: FakeBackendBehavior::Fail,
        };
        assert!(matches!(
            backend.solve(request),
            Err(BackendSolveError {
                error: ReactionExtentError::BackendFailure { .. },
                ..
            })
        ));
    }

    #[test]
    fn cascade_retains_wall_time_for_failed_attempt_without_metrics() {
        let residual = |values: &[f64]| Ok(values.to_vec());
        let feasible = |_values: &[f64]| true;
        let validate =
            |_values: &[f64]| -> Result<EquilibriumCandidateReport, ReactionExtentError> {
                panic!("failed backend must not reach candidate validation")
            };
        let backend = FakeBackend {
            backend: SolverBackend::Legacy(Solvers::LM),
            behavior: FakeBackendBehavior::FailAfterDelay,
        };
        let backends: [&dyn EquilibriumNonlinearBackend; 1] = [&backend];
        let error = solve_backend_cascade_with_control(
            &backends,
            vec![1.0],
            &residual,
            None,
            &feasible,
            &validate,
            SolverPolicy::Single(backend.backend),
            SolverCascadeBudget::new(1, 4, 4),
            &SolverParams::default(),
            None,
            None,
        )
        .expect_err("forced backend failure must reject the cascade");
        let ReactionExtentError::AllBackendsFailed { attempts } = error else {
            panic!("expected an all-backends failure, got {error:?}")
        };
        let SolverAttemptOutcome::Failed { elapsed_micros, .. } = attempts[0].outcome else {
            panic!("expected a failed attempt, got {:?}", attempts[0].outcome)
        };
        assert!(
            elapsed_micros >= 1_000,
            "failed attempt wall time must include backend execution: {elapsed_micros} us"
        );
    }

    #[test]
    fn failed_legacy_attempt_preserves_callback_metrics() {
        use crate::Thermodynamics::ChemEquilibrium::equilibrium_instrumentation::capture;

        let residual = |_values: &[f64]| {
            Err(ReactionExtentError::ResidualEvaluation(
                "forced legacy residual failure".to_string(),
            ))
        };
        let jacobian = |_values: &[f64]| Ok(DMatrix::identity(1, 1));
        let feasible = |_values: &[f64]| true;
        let backend = SolverBackend::Legacy(Solvers::LM);
        let (result, _) = capture(|| {
            backend.solve(BackendSolveRequest {
                initial_guess: vec![0.0],
                system: PreparedNonlinearSystem::new(
                    1,
                    &residual,
                    Some(&jacobian),
                    &feasible,
                    None,
                )
                .unwrap(),
                params: &SolverParams::default(),
                max_iterations: 4,
                rst_telemetry_mode: SolverTelemetryMode::Off,
            })
        });
        let error = result.expect_err("forced residual failure must be returned");
        let metrics = error
            .metrics
            .expect("failed Legacy attempt must retain diagnostics");
        assert!(!metrics.backend_converged);
        assert_eq!(metrics.residual_evaluations, 1);
        assert!(metrics.evaluation_timing.is_some());
    }

    #[test]
    fn legacy_solver_diagnostics_are_scoped_to_calculator_capture() {
        use crate::Thermodynamics::ChemEquilibrium::equilibrium_instrumentation::capture;

        let residual = |values: &[f64]| Ok(vec![values[0] - 1.0]);
        let jacobian = |_values: &[f64]| Ok(DMatrix::identity(1, 1));
        let feasible = |_values: &[f64]| true;
        for backend in [
            SolverBackend::Legacy(Solvers::LM),
            SolverBackend::Legacy(Solvers::NR),
            SolverBackend::Legacy(Solvers::TR),
        ] {
            let solve = || {
                backend.solve(BackendSolveRequest {
                    initial_guess: vec![0.0],
                    system: PreparedNonlinearSystem::new(
                        1,
                        &residual,
                        Some(&jacobian),
                        &feasible,
                        None,
                    )
                    .unwrap(),
                    params: &SolverParams::default(),
                    max_iterations: 20,
                    rst_telemetry_mode: SolverTelemetryMode::Off,
                })
            };

            let (captured, _) = capture(solve);
            let metrics = captured
                .unwrap()
                .metrics
                .expect("capture must attach legacy solver diagnostics");
            assert!(metrics.residual_evaluations > 0);
            assert!(metrics.jacobian_evaluations > 0);
            if !matches!(backend, SolverBackend::Legacy(Solvers::TR)) {
                assert!(metrics.linear_solves > 0);
            }
            let timing = metrics
                .evaluation_timing
                .expect("capture must attach legacy timing split");
            assert!(timing.linear_system_micros.is_some());
        }

        let uncaptured = SolverBackend::Legacy(Solvers::LM)
            .solve(BackendSolveRequest {
                initial_guess: vec![0.0],
                system: PreparedNonlinearSystem::new(
                    1,
                    &residual,
                    Some(&jacobian),
                    &feasible,
                    None,
                )
                .unwrap(),
                params: &SolverParams::default(),
                max_iterations: 20,
                rst_telemetry_mode: SolverTelemetryMode::Off,
            })
            .unwrap();
        assert!(
            uncaptured.metrics.is_none(),
            "legacy diagnostics must stay disabled outside calculator capture"
        );
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
            rst_telemetry_mode: SolverTelemetryMode::Off,
        };
        let backend = SolverBackend::Legacy(Solvers::LM);

        assert!(matches!(
            backend.solve(request),
            Err(BackendSolveError {
                error: ReactionExtentError::DimensionMismatch(message),
                ..
            }) if message.contains("prepared system requires 2")
        ));
    }
}
