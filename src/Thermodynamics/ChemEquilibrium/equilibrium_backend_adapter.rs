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
//!   EquilibriumLogMoles::solve_backend_cascade()
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

use crate::Thermodynamics::ChemEquilibrium::equilibrium_legacy_backend::solve_legacy_backend;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::SolverParams;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, RustedSciTheSolveContract,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverAttemptMetrics, SolverBackend,
};
use nalgebra::DMatrix;

type ResidualFn<'a> = dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError> + 'a;
type JacobianFn<'a> = dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError> + 'a;
type FeasibilityFn<'a> = dyn Fn(&[f64]) -> bool + 'a;

/// Mechanical request passed to one backend attempt.
pub(crate) struct BackendSolveRequest<'a> {
    pub initial_guess: Vec<f64>,
    pub residual: &'a ResidualFn<'a>,
    pub jacobian: Option<&'a JacobianFn<'a>>,
    pub feasible: &'a FeasibilityFn<'a>,
    pub params: &'a SolverParams,
    pub initial_moles: &'a [f64],
    pub reactions: &'a DMatrix<f64>,
    pub max_iterations: usize,
    pub rst_problem: Option<&'a RstPreparedProblem>,
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
        match *self {
            SolverBackend::Legacy(legacy_backend) => solve_legacy_backend(
                legacy_backend,
                request.initial_guess,
                request.residual,
                request.jacobian,
                request.feasible,
                request.params,
                request.initial_moles,
                request.reactions,
                request.max_iterations,
            )
            .map(|solution| BackendSolveResult {
                solution,
                metrics: None,
            }),
            SolverBackend::RustedSciThe(rst_solver) => match request.rst_problem {
                Some(problem) => rst_solver
                    .solve(
                        problem,
                        &request.initial_guess,
                        RustedSciTheSolveContract::new(request.params.tol, request.max_iterations)?,
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
            residual: &residual,
            jacobian: None,
            feasible: &feasible,
            params: &SolverParams::default(),
            initial_moles: &[1.0],
            reactions: &DMatrix::zeros(1, 1),
            max_iterations: 1,
            rst_problem: None,
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
            residual: &residual,
            jacobian: None,
            feasible: &feasible,
            params: &SolverParams::default(),
            initial_moles: &[1.0],
            reactions: &DMatrix::zeros(1, 1),
            max_iterations: 1,
            rst_problem: None,
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
}
