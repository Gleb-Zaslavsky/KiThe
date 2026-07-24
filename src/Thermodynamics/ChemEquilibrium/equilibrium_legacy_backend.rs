//! Adapter for the historical hand-written LM/NR/TR nonlinear solvers.
//!
//! # Purpose
//!
//! This module provides the **mechanical mapping** from a selected legacy
//! backend (LM, NR, or TR) to its concrete solver implementation in
//! [`equilibrium_nonlinear`](super::equilibrium_nonlinear). It is a thin
//! adapter that:
//!
//! 1. Receives a [`BackendSolveRequest`](super::equilibrium_backend_adapter::BackendSolveRequest).
//! 2. Constructs the appropriate solver (`LMSolver`, `NRSolver`, or `TrustRegionSolver`).
//! 3. Calls `solve()` and returns the result.
//!
//! # Status: DEPRECATED
//!
//! This module is retained for **regression testing and backward compatibility**
//! only. It can be deleted as one unit once the RustedSciThe symbolic backends
//! cover the retained regression matrix. The canonical engine owns policy
//! selection, budgets, candidate validation, and publication — this module
//! owns only the mechanical dispatch.
//!
//! # Key Functions
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`solve_legacy_backend`] | Dispatches to LM, NR, or TR solver |
//!
//! # Dataflow
//!
//! ```text
//!   EquilibriumLogMoles::solve_backend_cascade()
//!     │
//!     ├── Selects SolverBackend::Legacy(Solvers::LM/NR/TR)
//!     │
//!     v
//!   solve_legacy_backend(backend, initial_guess, residual, jacobian,
//!                        feasible, params, initial_moles, reactions, max_iter)
//!     │
//!     ├── Solvers::LM → LMSolver { f, jacobian, feasible, lambda, tol, max_iter, alpha_min }
//!     │                    └── solver.solve(initial_guess)
//!     │
//!     ├── Solvers::NR → NRSolver { f, jacobian, feasible, tol, max_iter, alpha_min }
//!     │                    └── solver.solve(initial_guess)
//!     │
//!     └── Solvers::TR → TrustRegionSolver { f, jacobian, feasible, tol, max_iter, delta_init, delta_max, eta }
//!                          └── solver.solve(initial_guess)
//!     │
//!     v
//!   Returns solution Vec<f64> or ReactionExtentError
//! ```
//!
//! # Non-obvious Details
//!
//! - The legacy backends **require an analytical Jacobian**. If `jacobian` is
//!   `None`, the function returns `ReactionExtentError::InvalidProblem`.
//! - The adapter has **no fallback logic** — it runs exactly one backend and
//!   returns. Cascade logic is in [`equilibrium_log_moles`](super::equilibrium_log_moles).
//! - The adapter does **not mutate** any solver instance — it is a pure function.
//!
//! # Related Modules
//!
//! - [`equilibrium_nonlinear`](super::equilibrium_nonlinear) — solver implementations
//! - [`equilibrium_backend_adapter`](super::equilibrium_backend_adapter) — adapter trait
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — orchestrator
//!

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{SolverParams, Solvers};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    LMSolver, NRSolver, ReactionExtentError, TrustRegionSolver,
};
use nalgebra::DMatrix;

type ResidualFn<'a> = dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError> + 'a;
type JacobianFn<'a> = dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError> + 'a;
type FeasibilityFn<'a> = dyn Fn(&[f64]) -> bool + 'a;

/// Runs exactly one requested legacy backend with the caller-provided budget.
///
/// The adapter deliberately has no fallback logic and never mutates a solver
/// instance. A missing Jacobian is a configuration error because all retained
/// legacy algorithms require one.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_legacy_backend(
    backend: Solvers,
    initial_guess: Vec<f64>,
    residual: &ResidualFn<'_>,
    jacobian: Option<&JacobianFn<'_>>,
    feasible: &FeasibilityFn<'_>,
    params: &SolverParams,
    initial_moles: &[f64],
    reactions: &DMatrix<f64>,
    max_iterations: usize,
) -> Result<Vec<f64>, ReactionExtentError> {
    let jacobian = jacobian.ok_or_else(|| ReactionExtentError::InvalidProblem {
        field: "legacy_jacobian",
        message: "a legacy backend requires the analytical Jacobian".to_string(),
    })?;

    match backend {
        Solvers::LM => {
            let mut solver = LMSolver {
                f: residual,
                jacobian,
                feasible,
                lambda: params.lambda,
                tol: 1e-12,
                max_iter: max_iterations,
                alpha_min: params.alpha_min,
            };
            solver
                .solve(initial_guess)
                .map_err(ReactionExtentError::SolveError)
        }
        Solvers::NR => {
            let mut solver = NRSolver {
                f: residual,
                jacobian,
                feasible,
                n0: initial_moles.to_vec(),
                reactions: reactions.clone(),
                tol: params.tol,
                max_iter: max_iterations,
                alpha_min: params.alpha_min,
            };
            solver
                .solve(initial_guess)
                .map_err(ReactionExtentError::SolveError)
        }
        Solvers::TR => {
            let solver = TrustRegionSolver {
                f: residual,
                jacobian,
                feasible,
                tol: params.tol,
                max_iter: max_iterations,
                delta_init: params.delta_init,
                delta_max: params.delta_max,
                eta: params.eta,
            };
            solver
                .solve(initial_guess)
                .map_err(ReactionExtentError::SolveError)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::solve_legacy_backend;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{SolverParams, Solvers};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use nalgebra::DMatrix;

    #[test]
    fn missing_legacy_jacobian_is_a_typed_configuration_error() {
        let residual = |_values: &[f64]| Ok(vec![0.0]);
        let feasible = |_values: &[f64]| true;

        assert!(matches!(
            solve_legacy_backend(
                Solvers::LM,
                vec![0.0],
                &residual,
                None,
                &feasible,
                &SolverParams::default(),
                &[1.0],
                &DMatrix::zeros(1, 0),
                1,
            ),
            Err(ReactionExtentError::InvalidProblem {
                field: "legacy_jacobian",
                ..
            })
        ));
    }
}
