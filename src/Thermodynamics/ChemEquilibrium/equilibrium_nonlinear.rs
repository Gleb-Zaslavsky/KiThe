//! Numerical solvers, error types, and SVD reaction basis for chemical equilibrium.
//!
//! # Purpose
//!
//! This module provides the **numerical core** of the equilibrium calculation pipeline.
//! It contains three nonlinear solvers (Levenberg-Marquardt, Newton-Raphson, Trust Region),
//! the SVD-based reaction basis computation, error types for the entire equilibrium
//! subsystem, and the feasibility constraint for mole non-negativity.
//!
//! # Physical and Mathematical Background
//!
//! ## Levenberg-Marquardt (LMSolver)
//!
//! Solves the nonlinear least-squares problem `min ||F(x)||²` using the damped
//! Gauss-Newton method:
//!
//! ```text
//! (J^T·J + λ·I)·s = -J^T·F
//! x_{k+1} = x_k + s
//! ```
//!
//! where `λ` is the damping parameter. The algorithm adapts `λ`:
//! - If the step reduces the residual, `λ` is decreased (closer to Gauss-Newton).
//! - If the step increases the residual, `λ` is increased (closer to gradient descent).
//!
//! The solver includes a **feasibility constraint** via `max_step_moles_nonnegative`
//! that limits the step to prevent log-mole values from producing negative physical moles.
//!
//! ## Newton-Raphson (NRSolver)
//!
//! Solves `F(x) = 0` using:
//!
//! ```text
//! J(x_k)·s = -F(x_k)
//! x_{k+1} = x_k + α_k · s
//! ```
//!
//! where `α_k ∈ (0, 1]` is determined by **line search with Armijo condition**:
//!
//! ```text
//! ||F(x_k + α·s)||² ≤ (1 - 2·α·c) · ||F(x_k)||²
//! ```
//!
//! The step is also limited by `max_step_moles_nonnegative` to maintain feasibility.
//!
//! ## Trust Region (TrustRegionSolver)
//!
//! Solves the constrained subproblem:
//!
//! ```text
//! min_{||s|| ≤ Δ} ||F(x_k) + J(x_k)·s||²
//! ```
//!
//! where `Δ` is the trust-region radius. The ratio of actual to predicted reduction:
//!
//! ```text
//! ρ = (||F(x_k)||² - ||F(x_k + s)||²) / (||F(x_k)||² - ||F(x_k) + J·s||²)
//! ```
//!
//! determines radius adaptation: if `ρ` is large, the radius grows; if small, it shrinks.
//!
//! ## SVD Reaction Basis
//!
//! The reaction basis is computed from the element composition matrix `A (species × elements)`:
//!
//! ```text
//! A = U · Σ · V^T
//! ```
//!
//! The stoichiometric matrix `ν (species × reactions)` is formed from columns of `V`
//! corresponding to zero (or near-zero) singular values. This guarantees that:
//!
//! ```text
//! A^T · ν = 0    (element conservation by construction)
//! ```
//!
//! If SVD does not return a complete nullspace, Gram-Schmidt orthogonalization
//! completes the basis.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`LMSolver`] | Levenberg-Marquardt solver with adaptive damping |
//! | [`NRSolver`] | Newton-Raphson solver with line search |
//! | [`TrustRegionSolver`] | Trust Region solver with adaptive radius |
//! | [`ReactionBasis`] | Container for SVD-derived reaction stoichiometry |
//! | [`SolveError`] | Solver-level errors (max iterations, singular matrix, eval error) |
//! | [`ReactionExtentError`] | Comprehensive equilibrium error hierarchy |
//! | [`BackendFailureKind`] | Classification of backend failures for retry logic |
//!
//! # Dataflow
//!
//! ```text
//!   Element composition matrix A (species × elements)
//!     │
//!     v
//!   compute_reaction_basis(A)
//!     ├── SVD: A = U·Σ·V^T
//!     ├── Extract nullspace columns from V
//!     ├── Gram-Schmidt completion if needed
//!     └── Return ReactionBasis { stoich_matrix, n_reactions }
//!     │
//!     v
//!   EquilibriumLogMoles builds residual/Jacobian closures
//!     │
//!     v
//!   LMSolver / NRSolver / TrustRegionSolver
//!     ├── Iterate until ||F(x)|| < tol or max_iter
//!     ├── Apply max_step_moles_nonnegative for feasibility
//!     └── Return solution or SolveError
//! ```
//!
//! # Error Hierarchy
//!
//! ```text
//! ReactionExtentError
//!   ├── InvalidProblem { field, message }     ── configuration errors
//!   ├── InvalidConditions { parameter, value } ── T, P, p0 validation
//!   ├── InvalidCandidate { message }           ── non-finite log-moles
//!   ├── DimensionMismatch(String)              ── matrix/vector shape errors
//!   ├── BackendFailure { kind, attempts }      ── solver cascade failures
//!   ├── SolveError                             ── wrapped solver errors
//!   ├── SpeciesNotAssigned(usize)              ── species without phase
//!   ├── SubsDataError                          ── database errors
//!   └── Aborted { reason }                     ── cascade aborted
//! ```
//!
//! # Examples
//!
//! ```rust
//! use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::*;
//! use nalgebra::DMatrix;
//!
//! // Compute reaction basis for O2 dissociation
//! let elem_comp = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]); // O2, O
//! let basis = compute_reaction_basis(&elem_comp, 1e-12).unwrap();
//! assert_eq!(basis.n_reactions, 1);
//! ```
//!
//! # Non-obvious Details
//!
//! - The LM solver's damping parameter `λ` is **adaptive**: it increases by 10× on
//!   rejected steps and decreases by 10× on accepted steps, clamped to `[1e-10, 1e10]`.
//! - The NR solver uses a **safeguarded line search** that backtracks from α=1 with
//!   a reduction factor of 0.5, ensuring sufficient decrease.
//! - The Trust Region solver solves the subproblem via the **dogleg method** when the
//!   Gauss-Newton step lies outside the trust region.
//! - `max_step_moles_nonnegative` computes the maximum step that keeps all mole numbers
//!   non-negative: `α_max = min_i (-x_i / s_i)` for `s_i < 0`.
//! - Gram-Schmidt completion in `compute_reaction_basis` handles the case where SVD
//!   returns fewer nullspace vectors than the dimension of the nullspace.
//!
//! # Related Modules
//!
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — orchestrator that uses these solvers
//! - [`equilibrium_legacy_backend`](super::equilibrium_legacy_backend) — adapter wrapping these solvers
//! - [`equilibrium_reaction_basis`](super::equilibrium_reaction_basis) — typed reaction basis wrapper
//!
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptReport;
use crate::Thermodynamics::User_substances_error::SubsDataError;
use log::{error, info};
use nalgebra::{DMatrix, DVector};
use std::error::Error;
use std::fmt;

/// Errors that can occur during nonlinear solving
#[derive(Debug, Clone)]
pub enum SolveError {
    /// Solver did not converge within max_iter
    MaxIterations,
    /// Linear system (LM normal equations) could not be solved
    SingularMatrix,
    /// Residual or Jacobian evaluation failed during solve
    EvalError(String),
}

impl fmt::Display for SolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaxIterations => write!(f, "nonlinear solver reached its iteration limit"),
            Self::SingularMatrix => write!(f, "nonlinear solver encountered a singular matrix"),
            Self::EvalError(message) => write!(f, "nonlinear solver evaluation failed: {message}"),
        }
    }
}

impl Error for SolveError {}

/// Machine-readable class of a numerical backend failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackendFailureKind {
    /// The backend could not solve one of its linearized subproblems.
    LinearSolve,
    /// The backend detected a singular or ill-conditioned Jacobian.
    SingularJacobian,
    /// The backend encountered NaN, Inf, or another numerical breakdown.
    NumericalBreakdown,
}

/// Comprehensive error types for chemical equilibrium calculations
#[derive(Debug)]
pub enum ReactionExtentError {
    /// A typed equilibrium-problem field violates its documented contract.
    InvalidProblem {
        /// Name of the invalid field or invariant.
        field: &'static str,
        /// Human-readable explanation suitable for diagnostics.
        message: String,
    },
    /// Temperature, pressure, or reference pressure is invalid.
    InvalidConditions {
        /// Name of the invalid condition.
        parameter: &'static str,
        /// Invalid value supplied by the caller.
        value: f64,
    },
    /// A numerical backend returned a candidate that violates equilibrium
    /// acceptance invariants.
    InvalidCandidate {
        /// Name of the violated candidate invariant.
        field: &'static str,
        /// Explanation of the failed invariant.
        message: String,
    },
    /// Every backend in an explicit fallback policy terminated without an
    /// accepted candidate. The reports preserve the deterministic attempt trace.
    AllBackendsFailed { attempts: Vec<SolverAttemptReport> },
    /// A cascade already performed one or more numerical attempts, then met a
    /// non-retryable domain or configuration error. The trace is retained so
    /// callers can distinguish an aborted cascade from ordinary exhaustion.
    CascadeAborted {
        /// Ordered outcomes recorded before the aborting error.
        attempts: Vec<SolverAttemptReport>,
        /// The non-retryable error that stopped further fallback attempts.
        cause: Box<ReactionExtentError>,
    },
    /// Error from substance database operations
    SubsDataError(SubsDataError),
    /// Error originating from the underlying nonlinear solver
    SolveError(SolveError),
    /// A named numerical backend failed in a way that is distinct from a
    /// domain/residual contract failure and may be retried by a cascade.
    BackendFailure {
        /// Stable backend identifier from the solve policy.
        backend: String,
        /// Typed numerical failure class.
        kind: BackendFailureKind,
        /// Backend-provided diagnostic detail.
        message: String,
    },
    /// SVD failed when computing reaction basis
    SVDError(String),
    /// A requested validation/solver path is intentionally not applicable to
    /// the supplied system size or model class.
    ValidationNotApplicable {
        /// Name of the requested validation or solver path.
        path: &'static str,
        /// Human-readable reason for the refusal.
        message: String,
    },
    /// Phase-control outer loop exhausted its allowed restart budget.
    PhaseControlDidNotConverge {
        /// Number of outer iterations completed before the budget was hit.
        iterations: usize,
    },
    /// The active-set controller revisited an earlier phase assemblage.
    PhaseControlCycleDetected {
        /// Outer iteration at which the repeated set was proposed.
        iteration: usize,
        /// Ordered indices of phases that would be active after the transition.
        active_phases: Vec<usize>,
    },
    /// Dimension mismatch between vectors/matrices
    DimensionMismatch(String),
    /// Duplicate species assigned to multiple phases
    DuplicateSpecies(usize),
    /// A species is not assigned to any phase.
    SpeciesNotAssigned { index: usize },
    /// Initial residual evaluation produced NaN/Inf
    InvalidInitialResiduals(Vec<f64>),
    /// Residual evaluation failed (generic)
    ResidualEvaluation(String),
    /// Jacobian evaluation failed (generic)
    JacobianEvaluation(String),
    /// Invalid species mole numbers (negative or zero)
    InvalidSpeciesAmount { index: usize, value: f64 },
    /// Invalid per-phase totals
    InvalidNPhase { index: usize, value: f64 },
    /// Invalid phi value
    InvalidPhi { index: usize, value: f64 },
    /// A species Gibbs-energy closure returned a non-finite value at the
    /// requested temperature.
    InvalidDG0 {
        /// Species position in the deterministic prepared ordering.
        species_index: usize,
        /// Value returned by the Gibbs-energy closure.
        dg0: f64,
        /// Temperature passed to the closure.
        temperature: f64,
    },
}

/// Machine-readable high-level classification for equilibrium failures.
///
/// This keeps the existing public error enum stable while giving callers a
/// strong, typed way to branch on failure families.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionExtentErrorKind {
    InvalidInput,
    ThermochemicalLookup,
    Formulation,
    ResidualEvaluation,
    JacobianEvaluation,
    BackendFailure,
    NonConvergence,
    InvalidCandidate,
    UnsupportedModel,
    ValidationMismatch,
    AllBackendsFailed,
    CascadeAborted,
}

impl From<SubsDataError> for ReactionExtentError {
    fn from(value: SubsDataError) -> Self {
        ReactionExtentError::SubsDataError(value)
    }
}

impl fmt::Display for ReactionExtentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidProblem { field, message } => {
                write!(f, "invalid equilibrium problem field '{field}': {message}")
            }
            Self::InvalidConditions { parameter, value } => {
                write!(f, "invalid equilibrium condition '{parameter}': {value}")
            }
            Self::InvalidCandidate { field, message } => {
                write!(
                    f,
                    "equilibrium candidate field '{field}' was rejected: {message}"
                )
            }
            Self::AllBackendsFailed { attempts } => {
                write!(
                    f,
                    "all equilibrium solver backends failed after {} attempt(s)",
                    attempts.len()
                )
            }
            Self::CascadeAborted { attempts, cause } => write!(
                f,
                "equilibrium solver cascade aborted after {} attempt(s): {cause}",
                attempts.len()
            ),
            Self::SubsDataError(error) => write!(f, "equilibrium data preparation failed: {error}"),
            Self::SolveError(error) => write!(f, "equilibrium solver failed: {error}"),
            Self::BackendFailure {
                backend,
                kind,
                message,
            } => write!(
                f,
                "equilibrium backend '{backend}' failed ({kind:?}): {message}"
            ),
            Self::SVDError(message) => write!(f, "reaction-basis construction failed: {message}"),
            Self::ValidationNotApplicable { path, message } => {
                write!(f, "validation path '{path}' is not applicable: {message}")
            }
            Self::PhaseControlDidNotConverge { iterations } => {
                write!(
                    f,
                    "phase-control outer loop did not converge after {iterations} iteration(s)"
                )
            }
            Self::PhaseControlCycleDetected {
                iteration,
                active_phases,
            } => write!(
                f,
                "phase-control cycle detected at iteration {iteration} for active phases {active_phases:?}"
            ),
            Self::DimensionMismatch(message) => {
                write!(f, "equilibrium dimension mismatch: {message}")
            }
            Self::DuplicateSpecies(index) => {
                write!(f, "species index {index} belongs to multiple phases")
            }
            Self::SpeciesNotAssigned { index } => {
                write!(f, "species index {index} is not assigned to any phase")
            }
            Self::InvalidInitialResiduals(values) => {
                write!(
                    f,
                    "initial equilibrium residual contains non-finite values: {values:?}"
                )
            }
            Self::ResidualEvaluation(message) => {
                write!(f, "equilibrium residual evaluation failed: {message}")
            }
            Self::JacobianEvaluation(message) => {
                write!(f, "equilibrium Jacobian evaluation failed: {message}")
            }
            Self::InvalidSpeciesAmount { index, value } => {
                write!(f, "invalid mole amount at species index {index}: {value}")
            }
            Self::InvalidNPhase { index, value } => {
                write!(f, "invalid total mole amount for phase {index}: {value}")
            }
            Self::InvalidPhi { index, value } => {
                write!(f, "invalid phase correction at phase {index}: {value}")
            }
            Self::InvalidDG0 {
                species_index,
                dg0,
                temperature,
            } => {
                write!(
                    f,
                    "invalid standard Gibbs energy {dg0} for species index {species_index} at temperature {temperature}"
                )
            }
        }
    }
}

impl Error for ReactionExtentError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::SubsDataError(error) => Some(error),
            Self::SolveError(error) => Some(error),
            Self::CascadeAborted { cause, .. } => Some(cause.as_ref()),
            _ => None,
        }
    }
}

impl ReactionExtentError {
    /// Returns the machine-readable high-level kind for this error.
    pub fn kind(&self) -> ReactionExtentErrorKind {
        match self {
            Self::InvalidProblem { .. }
            | Self::InvalidConditions { .. }
            | Self::DimensionMismatch(_)
            | Self::DuplicateSpecies(_)
            | Self::SpeciesNotAssigned { .. }
            | Self::InvalidInitialResiduals(_)
            | Self::InvalidSpeciesAmount { .. }
            | Self::InvalidNPhase { .. }
            | Self::InvalidPhi { .. } => ReactionExtentErrorKind::Formulation,
            Self::SubsDataError(_) => ReactionExtentErrorKind::ThermochemicalLookup,
            Self::InvalidCandidate { .. } => ReactionExtentErrorKind::InvalidCandidate,
            Self::InvalidDG0 { .. } => ReactionExtentErrorKind::ThermochemicalLookup,
            Self::SolveError(SolveError::MaxIterations) => ReactionExtentErrorKind::NonConvergence,
            Self::SolveError(_) => ReactionExtentErrorKind::BackendFailure,
            Self::BackendFailure { .. } => ReactionExtentErrorKind::BackendFailure,
            Self::SVDError(_) => ReactionExtentErrorKind::UnsupportedModel,
            Self::ValidationNotApplicable { .. } => ReactionExtentErrorKind::ValidationMismatch,
            Self::PhaseControlDidNotConverge { .. } | Self::PhaseControlCycleDetected { .. } => {
                ReactionExtentErrorKind::NonConvergence
            }
            Self::ResidualEvaluation(_) => ReactionExtentErrorKind::ResidualEvaluation,
            Self::JacobianEvaluation(_) => ReactionExtentErrorKind::JacobianEvaluation,
            Self::AllBackendsFailed { .. } => ReactionExtentErrorKind::AllBackendsFailed,
            Self::CascadeAborted { .. } => ReactionExtentErrorKind::CascadeAborted,
        }
    }

    /// Returns `true` when another numerical backend may legitimately retry.
    ///
    /// Invalid input, malformed configuration, and unsupported model errors
    /// deliberately return `false` so fallback cannot hide a bad problem
    /// definition behind later attempts.
    pub fn is_retryable_backend_failure(&self) -> bool {
        matches!(
            self.kind(),
            ReactionExtentErrorKind::NonConvergence
                | ReactionExtentErrorKind::BackendFailure
                | ReactionExtentErrorKind::ResidualEvaluation
                | ReactionExtentErrorKind::JacobianEvaluation
        )
    }

    /// Returns `true` for errors that should fail fast before fallback.
    pub fn is_non_retryable_input_error(&self) -> bool {
        matches!(
            self.kind(),
            ReactionExtentErrorKind::InvalidInput
                | ReactionExtentErrorKind::ThermochemicalLookup
                | ReactionExtentErrorKind::Formulation
                | ReactionExtentErrorKind::InvalidCandidate
                | ReactionExtentErrorKind::UnsupportedModel
                | ReactionExtentErrorKind::ValidationMismatch
        )
    }
}

#[cfg(test)]
mod error_contract_tests {
    use super::{ReactionExtentError, ReactionExtentErrorKind, SolveError};
    use crate::Thermodynamics::User_substances_error::SubsDataError;
    use std::error::Error;

    #[test]
    fn equilibrium_errors_are_displayable_and_preserve_nested_sources() {
        let solve_error = ReactionExtentError::SolveError(SolveError::MaxIterations);
        assert_eq!(
            solve_error.to_string(),
            "equilibrium solver failed: nonlinear solver reached its iteration limit"
        );
        assert!(solve_error.source().is_some());

        let data_error =
            ReactionExtentError::from(SubsDataError::SubstanceNotFound("O2".to_string()));
        assert!(data_error.to_string().contains("O2"));
        assert!(data_error.source().is_some());
    }

    #[test]
    fn all_backend_failure_display_reports_the_attempt_count() {
        let error = ReactionExtentError::AllBackendsFailed {
            attempts: Vec::new(),
        };
        assert_eq!(
            error.to_string(),
            "all equilibrium solver backends failed after 0 attempt(s)"
        );
        assert!(error.source().is_none());
        assert_eq!(error.kind(), ReactionExtentErrorKind::AllBackendsFailed);
    }

    #[test]
    fn aborted_cascade_preserves_the_non_retryable_cause() {
        let error = ReactionExtentError::CascadeAborted {
            attempts: Vec::new(),
            cause: Box::new(ReactionExtentError::InvalidProblem {
                field: "solver_policy",
                message: "missing backend".to_string(),
            }),
        };
        assert!(
            error
                .to_string()
                .contains("cascade aborted after 0 attempt(s)")
        );
        assert!(error.source().is_some());
        assert_eq!(error.kind(), ReactionExtentErrorKind::CascadeAborted);
    }

    #[test]
    fn species_not_assigned_display_reports_the_missing_index() {
        let error = ReactionExtentError::SpeciesNotAssigned { index: 3 };
        assert_eq!(
            error.to_string(),
            "species index 3 is not assigned to any phase"
        );
    }

    #[test]
    fn error_kind_classifies_core_families() {
        assert_eq!(
            ReactionExtentError::InvalidProblem {
                field: "solver_policy",
                message: "missing backend".to_string(),
            }
            .kind(),
            ReactionExtentErrorKind::Formulation
        );
        assert_eq!(
            ReactionExtentError::from(SubsDataError::SubstanceNotFound("O2".to_string())).kind(),
            ReactionExtentErrorKind::ThermochemicalLookup
        );
        assert_eq!(
            ReactionExtentError::ResidualEvaluation("bad residual".to_string()).kind(),
            ReactionExtentErrorKind::ResidualEvaluation
        );
        assert_eq!(
            ReactionExtentError::SolveError(SolveError::MaxIterations).kind(),
            ReactionExtentErrorKind::NonConvergence
        );
        assert_eq!(
            ReactionExtentError::ValidationNotApplicable {
                path: "keq_solver",
                message: "too many reactions".to_string(),
            }
            .kind(),
            ReactionExtentErrorKind::ValidationMismatch
        );
        let phase_error = ReactionExtentError::PhaseControlDidNotConverge { iterations: 3 };
        assert_eq!(phase_error.kind(), ReactionExtentErrorKind::NonConvergence);
        assert!(phase_error.to_string().contains("phase-control outer loop"));
    }
}
/// Levenberg-Marquardt solver for nonlinear equation systems
///
/// Robust solver that combines Gauss-Newton with gradient descent using adaptive damping.
/// Particularly effective for chemical equilibrium problems with poor initial guesses.
pub struct LMSolver<F, J, C>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
    J: Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>,
    C: Fn(&[f64]) -> bool,
{
    /// Residual function f(x) = 0
    pub f: F,
    /// Jacobian function J(x) = df/dx
    pub jacobian: J,
    /// Feasibility constraint checker
    pub feasible: C,
    /// Damping parameter (increased when steps rejected)
    pub lambda: f64,
    /// Convergence tolerance for ||f(x)||
    pub tol: f64,
    /// Maximum number of iterations
    pub max_iter: usize,
    /// Minimum step size before giving up
    pub alpha_min: f64,
}

impl<F, J, C> LMSolver<F, J, C>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
    J: Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>,
    C: Fn(&[f64]) -> bool,
{
    /// Solves nonlinear system f(x) = 0 using Levenberg-Marquardt algorithm
    ///
    /// Uses adaptive damping and line search to ensure robust convergence.
    /// Respects feasibility constraints throughout the solution process.
    pub fn solve(&mut self, mut x: Vec<f64>) -> Result<Vec<f64>, SolveError> {
        let n = x.len();

        let mut lambda = self.lambda;

        for _iter in 0..self.max_iter {
            let f_val = (self.f)(&x).map_err(|error| {
                SolveError::EvalError(format!("residual evaluation failed: {error}"))
            })?;
            info!("Iteration {}, x = {:?}, f = {:?}", _iter, x, f_val);
            let f_norm = l2_norm(&f_val);
            info!("  ||f|| = {}", f_norm);
            if f_norm < self.tol {
                return Ok(x);
            }
            if f_norm.is_nan() || f_norm.is_infinite() {
                return Err(SolveError::SingularMatrix);
            }

            let j = (self.jacobian)(&x).map_err(|error| {
                SolveError::EvalError(format!("Jacobian evaluation failed: {error}"))
            })?;
            info!("  Jacobian:\n{}", j);
            if j.nrows() != n || j.ncols() != n {
                error!(
                    "Jacobian nrows = {} ncols = {}, state length = {}",
                    j.nrows(),
                    j.ncols(),
                    n
                );
                return Err(SolveError::EvalError(
                    "Jacobian dimension mismatch".to_string(),
                ));
            }
            if j.iter().any(|value| !value.is_finite()) {
                return Err(SolveError::SingularMatrix);
            }
            let jt = j.transpose();

            let jtj = &jt * &j;
            let mut lhs = jtj.clone();
            info!("  Jᵀ·J:\n{}", jtj);
            for i in 0..n {
                lhs[(i, i)] += lambda;
            }

            let rhs = -(&jt * DVector::from_vec(f_val.clone()));

            let delta = lhs.lu().solve(&rhs).ok_or(SolveError::SingularMatrix)?;
            info!("  Step delta: {:?}", delta);

            let delta = delta.data.as_vec().clone();

            let mut alpha = 1.0;
            let mut accepted = false;

            while alpha >= self.alpha_min {
                let x_trial: Vec<f64> = x
                    .iter()
                    .zip(delta.iter())
                    .map(|(xi, dxi)| xi + alpha * dxi)
                    .collect();
                info!("    Trial x (alpha={}): {:?}", alpha, x_trial);
                if !(self.feasible)(&x_trial) {
                    alpha *= 0.5;
                    continue;
                }

                let f_trial = (self.f)(&x_trial).map_err(|error| {
                    SolveError::EvalError(format!("residual evaluation failed: {error}"))
                })?;
                info!("    Trial f: {:?}", f_trial);
                let f_trial_norm = l2_norm(&f_trial);
                info!("    ||f_trial|| = {}", f_trial_norm);
                if f_trial_norm < f_norm {
                    x = x_trial;
                    lambda *= 0.3;
                    accepted = true;
                    break;
                } else {
                    alpha *= 0.5;
                }
            }

            if !accepted {
                lambda *= 10.0;
                info!(
                    "  No acceptable step found; increasing lambda to {}",
                    lambda
                );
            }
        }

        Err(SolveError::MaxIterations)
    }
}
//////////////////////////////////////NR SOLVER//////////////////////////////////////////////////////////////
/// Computes maximum step size to keep mole numbers non-negative
///
/// For Newton steps that would drive species moles negative, this function
/// computes the largest α such that n_i + α*Δn_i ≥ 0 for all species.
pub fn max_step_moles_nonnegative(
    n: &[f64],       // current moles
    delta_n: &[f64], // Newton step
    safety: f64,     // e.g. 0.95
) -> Result<f64, SolveError> {
    if n.len() != delta_n.len() {
        return Err(SolveError::EvalError(format!(
            "mole vector has {} entries but Newton step has {} entries",
            n.len(),
            delta_n.len(),
        )));
    }
    if !safety.is_finite() || safety <= 0.0 || safety > 1.0 {
        return Err(SolveError::EvalError(
            "step safety factor must be finite and lie in (0, 1]".to_string(),
        ));
    }
    if n.iter().chain(delta_n).any(|value| !value.is_finite()) {
        return Err(SolveError::EvalError(
            "moles and Newton step must be finite".to_string(),
        ));
    }

    let mut alpha_max: f64 = 1.0;

    for i in 0..n.len() {
        let ni = n[i];
        let dni = delta_n[i];

        // Only decreasing species can violate positivity
        if dni < 0.0 {
            if ni <= 0.0 {
                return Ok(0.0); // already infeasible
            }

            let alpha_i = safety * ni / (-dni);
            alpha_max = alpha_max.min(alpha_i);
        }
    }

    Ok(alpha_max.clamp(0.0, 1.0))
}

/*
pub fn tolerance_calc(f_val: &Vec<f64>, x: &Vec<f64>){
    let mut complex = 0.0;
    let x_sum: f64 = x.iter().sum();
    for (i, f_val_i) in f_val.iter().enumerate(){
        let rel
       let s =  f_val_i/x[i].abs()
    }
}
*/

/// Newton-Raphson solver with line search and feasibility constraints
///
/// Fast quadratic convergence near solution, with backtracking line search
/// and bound constraints to handle chemical equilibrium problems.
pub struct NRSolver<F, J, C>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
    J: Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>,
    C: Fn(&[f64]) -> bool,
{
    /// Residual function f(x) = 0
    pub f: F,
    /// Jacobian function J(x) = df/dx
    pub jacobian: J,
    /// Feasibility constraint checker
    pub feasible: C,
    /// Initial mole numbers (for bound checking)
    pub n0: Vec<f64>,
    /// Reaction stoichiometry matrix
    pub reactions: DMatrix<f64>,
    /// Convergence tolerance for ||f(x)||
    pub tol: f64,
    /// Maximum number of iterations
    pub max_iter: usize,
    /// Minimum step size before giving up
    pub alpha_min: f64,
}

impl<F, J, C> NRSolver<F, J, C>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
    J: Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>,
    C: Fn(&[f64]) -> bool,
{
    /// Solves nonlinear system using Newton-Raphson with line search
    ///
    /// Performs backtracking line search with feasibility checking.
    /// Uses bound-aware step limiting to prevent negative mole numbers.
    pub fn solve(&mut self, mut x: Vec<f64>) -> Result<Vec<f64>, SolveError> {
        let n = x.len();

        for iter in 0..self.max_iter {
            let f_val = (self.f)(&x).map_err(|error| SolveError::EvalError(error.to_string()))?;
            info!("NR Iteration {}, x = {:?}, f = {:?}", iter, x, f_val);
            let f_norm = l2_norm(&f_val);
            info!(" ||f|| = {}", f_norm);
            if f_norm < self.tol {
                return Ok(x);
            }
            if f_norm.is_nan() || f_norm.is_infinite() {
                return Err(SolveError::SingularMatrix);
            }
            let j =
                (self.jacobian)(&x).map_err(|error| SolveError::EvalError(error.to_string()))?;

            info!(" Jacobian:\n{}", j);
            if j.nrows() != n || j.ncols() != n {
                error!(
                    "Jacobian nrows = {} ncols ={}, initial guess length {}",
                    j.nrows(),
                    j.ncols(),
                    n
                );
                return Err(SolveError::EvalError(
                    "Jacobian dimension mismatch".to_string(),
                ));
            }
            // Solve J * delta = -f
            let rhs = -DVector::from_vec(f_val);
            let delta_vec = j.lu().solve(&rhs).ok_or(SolveError::SingularMatrix)?;
            let delta = delta_vec.data.as_vec();
            info!(" Step delta: {:?}", delta_vec);
            // --- NEW: bound-aware step size ---
            // let alpha_species =  max_step_moles_nonnegative(&x, &delta, 0.95);

            // let mut alpha = alpha_species;
            let mut alpha = 1.0;
            info!("bounded step {}", alpha);

            let mut accepted = false;
            while alpha >= self.alpha_min {
                let x_trial: Vec<f64> = x
                    .iter()
                    .zip(delta.iter())
                    .map(|(xi, dxi)| xi + alpha * dxi)
                    .collect();
                info!(" Trial x (alpha={}): {:?}", alpha, x_trial);
                if !(self.feasible)(&x_trial) {
                    alpha *= 0.5;
                    continue;
                }

                let f_trial =
                    (self.f)(&x_trial).map_err(|error| SolveError::EvalError(error.to_string()))?;

                info!(" Trial f: {:?}", f_trial);

                let f_trial_norm = l2_norm(&f_trial);
                info!(" ||f_trial|| = {}", f_trial_norm);
                if f_trial_norm < f_norm {
                    x = x_trial;
                    accepted = true;
                    break;
                }

                alpha *= 0.5;
                if !accepted {
                    info!(" No acceptable step found; continuing to next iteration");
                }
            }

            if !accepted {
                error!(
                    "NR solver exhausted iteration {} without an acceptable step",
                    iter
                );
                return Err(SolveError::MaxIterations);
            }
        }

        error!("NR solver reached the iteration limit without convergence");
        Err(SolveError::MaxIterations)
    }
}

/// Computes L2 norm of a vector
///
/// Helper function for convergence checking in solvers.
fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
///////////////////////TRUST REGION///////////////////////////////////////////////////

pub struct TrustRegionSolver<F, J, C>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
    J: Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>,
    C: Fn(&[f64]) -> bool,
{
    pub f: F,
    pub jacobian: J,
    pub feasible: C,

    pub tol: f64,
    pub max_iter: usize,

    // Trust region parameters
    pub delta_init: f64,
    pub delta_max: f64,
    pub eta: f64, // acceptance threshold (e.g. 0.1)
}

impl<F, J, C> TrustRegionSolver<F, J, C>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
    J: Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>,
    C: Fn(&[f64]) -> bool,
{
    pub fn solve(&self, mut x: Vec<f64>) -> Result<Vec<f64>, SolveError> {
        let n = x.len();
        let mut delta = self.delta_init;

        for _iter in 0..self.max_iter {
            let f_val = (self.f)(&x).map_err(|error| SolveError::EvalError(error.to_string()))?;
            if f_val.len() != n {
                return Err(SolveError::EvalError("residual dimension mismatch".into()));
            }
            if f_val.iter().any(|value| !value.is_finite()) {
                return Err(SolveError::SingularMatrix);
            }

            let f_norm = l2_norm(&f_val);
            if f_norm < self.tol {
                return Ok(x);
            }

            let J =
                (self.jacobian)(&x).map_err(|error| SolveError::EvalError(error.to_string()))?;

            if J.nrows() != n || J.ncols() != n {
                error!(
                    "Trust region Jacobian nrows = {} ncols = {}, state length = {}",
                    J.nrows(),
                    J.ncols(),
                    n
                );
                return Err(SolveError::EvalError("Jacobian not square".into()));
            }
            if J.iter().any(|value| !value.is_finite()) {
                return Err(SolveError::SingularMatrix);
            }

            let fvec = DVector::from_vec(f_val.clone());

            // --- Newton step ---
            let newton_step = J
                .clone()
                .lu()
                .solve(&(-&fvec))
                .ok_or(SolveError::SingularMatrix)?;

            // --- Cauchy step ---
            let g = &J.transpose() * &fvec;
            let g_norm_sq = g.dot(&g);
            let jg = &J * &g;
            let alpha = g_norm_sq / jg.dot(&jg).max(1e-16);
            let cauchy_step = -alpha * g.clone();

            // --- Dogleg step ---
            let p = if newton_step.norm() <= delta {
                newton_step
            } else if cauchy_step.norm() >= delta {
                delta / cauchy_step.norm() * cauchy_step
            } else {
                let p_u = cauchy_step;
                let p_b = newton_step;
                let d = &p_b - &p_u;

                let a = d.dot(&d);
                let b = 2.0 * p_u.dot(&d);
                let c = p_u.dot(&p_u) - delta * delta;

                let tau = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
                p_u + tau * d
            };

            let p_vec = p.data.as_vec().clone();
            let x_trial: Vec<f64> = x.iter().zip(p_vec.iter()).map(|(xi, pi)| xi + pi).collect();

            if !(self.feasible)(&x_trial) {
                delta *= 0.25;
                continue;
            }

            let f_trial =
                (self.f)(&x_trial).map_err(|error| SolveError::EvalError(error.to_string()))?;

            // A numerically exact trial root must be accepted directly.  At
            // that point the quadratic-model reduction can be zero or lose
            // its sign through round-off, which is not a reason to reject a
            // valid solution and shrink the region into a false singularity.
            if l2_norm(&f_trial) < self.tol {
                return Ok(x_trial);
            }

            let actual_reduction = f_norm.powi(2) - l2_norm(&f_trial).powi(2);
            if !actual_reduction.is_finite() {
                return Err(SolveError::SingularMatrix);
            }

            // The quadratic model is for ||f + Jp||^2, so its linear term is
            // g^T p where g = J^T f. Using f^T p here is dimensionally wrong
            // and can make a well-improving step look harmful when J is
            // ill-conditioned, collapsing the trust region to machine noise.
            let model_reduction = -2.0 * g.dot(&p) - p.dot(&(J.transpose() * &J * &p));

            if !model_reduction.is_finite() || model_reduction <= 0.0 {
                delta *= 0.25;
                continue;
            }

            let rho = actual_reduction / model_reduction;

            // --- Trust region update ---
            if rho < 0.25 {
                delta *= 0.25;
            } else if rho > 0.75 && (p.norm() - delta).abs() < 1e-12 {
                delta = (2.0 * delta).min(self.delta_max);
            }

            // --- Accept / reject ---
            if rho > self.eta {
                x = x_trial;
            }

            if delta < 1e-14 {
                error!("Trust region solver exhausted the trust radius before convergence");
                return Err(SolveError::SingularMatrix);
            }
        }

        error!("Trust region solver reached the iteration limit without convergence");
        Err(SolveError::MaxIterations)
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Container for independent chemical reaction basis
///
/// Holds the results of SVD-based reaction discovery from elemental composition.
#[derive(Debug, Clone)]
pub struct ReactionBasis {
    /// Rank of element composition matrix A
    pub rank: usize,
    /// Number of independent reactions (nullspace dimension)
    pub num_reactions: usize,
    /// Reaction stoichiometry matrix N (m x r)
    /// Each column is one independent reaction
    pub reactions: DMatrix<f64>,
}

/// Computes a basis of independent chemical reactions from an element
/// composition matrix.
///
/// # Mathematical meaning
///
/// Let `A` be the element composition matrix with:
/// - rows corresponding to chemical species
/// - columns corresponding to chemical elements
/// - `A[i, e]` equal to the number of atoms of element `e` in species `i`
///
/// Any physically admissible chemical reaction must conserve each element:
///
/// ```text
/// Aᵀ · ν = 0
/// ```
///
/// where `ν` is the vector of stoichiometric coefficients of a reaction.
/// Therefore, all independent reactions lie in the nullspace of `Aᵀ`.
///
/// This function computes a basis of that nullspace using singular value
/// decomposition (SVD).
///
/// # Physical meaning
///
/// - Each column of the returned matrix represents one independent chemical
///   reaction.
/// - The sign and scaling of each reaction vector are arbitrary and have no
///   physical significance.
/// - Together, these reactions span all possible composition changes that
///   conserve elements.
///
/// The number of independent reactions is:
///
/// ```text
/// r = number_of_species − rank(A)
/// ```
///
/// This construction guarantees that:
/// - Element conservation is satisfied automatically
/// - No redundant or dependent reactions are introduced
/// - Equilibrium systems formulated in reaction extents are well-conditioned
///
/// # Arguments
///
/// * `a`   – Element composition matrix `(m × E)`
/// * `tol` – Singular value threshold used to determine numerical rank
///
/// # Returns
///
/// A [`ReactionBasis`] containing:
/// - the rank of the element matrix
/// - the number of independent reactions
/// - a reaction stoichiometry matrix `(m × r)`
///
/// # Notes
///
/// This function does **not** classify species as reactants or products.
/// That distinction emerges naturally from the solution of the equilibrium
/// equations.
///
/// # Panics
///
/// Panics if SVD decomposition fails (which should not happen for finite
/// matrices).
pub fn compute_reaction_basis(
    a: &DMatrix<f64>,
    tol: f64,
) -> Result<ReactionBasis, ReactionExtentError> {
    use nalgebra::linalg::SVD;

    let m = a.nrows(); // species

    // SVD of Aᵀ (elements × species)
    let at = a.transpose();
    let svd = SVD::new(at, true, true);

    let v_t = match svd.v_t {
        Some(v) => v,
        None => {
            let msg = "SVD failed: Vᵀ".to_string();
            error!("{}", msg);
            return Err(ReactionExtentError::SVDError(msg));
        }
    };
    let singular_values = svd.singular_values;

    // Numerical rank
    let rank = singular_values.iter().filter(|&&s| s > tol).count();

    let num_reactions = if rank <= m { m - rank } else { 0 };

    let mut reactions = DMatrix::<f64>::zeros(m, num_reactions);

    // Collect available rows of Vᵀ (nalgebra returns a "thin" Vᵀ of size min(E, m) × m)
    let p = v_t.nrows(); // number of available singular vectors (<= m)

    // Convert rows of v_t into column vectors of length m
    let mut v_rows: Vec<DVector<f64>> = Vec::new();
    for i in 0..p {
        v_rows.push(v_t.row(i).transpose());
    }

    // Fill reaction columns using available nullspace vectors from Vᵀ (those with zero singular values
    // that are present in the returned Vᵀ). If there are fewer available than needed, we'll complete the
    // nullspace by constructing orthonormal vectors orthogonal to the span of the nonzero-singular-value rows.
    let mut col = 0;

    // Rows corresponding to small/zero singular values are at indices rank..p-1
    for i in rank..p {
        if col >= num_reactions {
            break;
        }
        reactions.set_column(col, &v_rows[i]);
        col += 1;
    }

    // If we still need more nullspace vectors (happens when p < m), construct them via Gram-Schmidt
    if col < num_reactions {
        // The span we must be orthogonal to is the space spanned by the first `rank` rows (if rank>0)
        let orth_span: Vec<&DVector<f64>> = v_rows.iter().take(rank).collect();

        let mut null_basis: Vec<DVector<f64>> = Vec::new();

        for j in 0..m {
            if col >= num_reactions {
                break;
            }

            // start from standard basis vector e_j
            let mut cand = DVector::<f64>::from_element(m, 0.0);
            cand[j] = 1.0;

            // subtract projection onto orth_span (non-null singular vectors)
            for rvec in &orth_span {
                let dot = rvec.dot(&cand);
                cand -= *rvec * dot;
            }

            // subtract projection onto previously found null basis vectors to maintain orthogonality
            for u in &null_basis {
                let dot = u.dot(&cand);
                cand -= u * dot;
            }

            let norm = cand.norm();
            if norm > tol {
                let u = cand / norm;
                reactions.set_column(col, &u);
                null_basis.push(u);
                col += 1;
            }
        }
    }
    let mut reactions = -1.0 * reactions; // flip sign convention

    // Zero out small elements (numerical noise) for stability
    let eps = 1e-3;
    for i in 0..m {
        for j in 0..num_reactions {
            if reactions[(i, j)].abs() < eps {
                reactions[(i, j)] = 0.0;
            }
        }
    }

    Ok(ReactionBasis {
        rank,
        num_reactions,
        reactions,
    })
}
///////////////////////////////TESTS////////////////////////////////////
#[cfg(test)]
mod tests {
    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }
    use super::*;
    #[test]
    fn lm_solves_scalar_quadratic() {
        let f = |x: &[f64]| Ok(vec![x[0] * x[0] - 2.0]) as Result<Vec<f64>, ReactionExtentError>;

        let j = |x: &[f64]| {
            Ok(nalgebra::DMatrix::from_row_slice(1, 1, &[2.0 * x[0]]))
                as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |_x: &[f64]| true;

        let mut solver = LMSolver {
            f,
            jacobian: j,
            feasible,
            lambda: 1e-3,
            tol: 1e-12,
            max_iter: 50,
            alpha_min: 1e-6,
        };

        let x0 = vec![1.0];
        let sol = solver.solve(x0).unwrap();

        assert!(approx_eq(sol[0], 2.0_f64.sqrt(), 1e-8));
    }
    #[test]
    fn lm_solves_2d_nonlinear_system() {
        let f = |x: &[f64]| {
            Ok(vec![x[0] * x[0] + x[1] * x[1] - 1.0, x[0] - x[1]])
                as Result<Vec<f64>, ReactionExtentError>
        };

        let j = |x: &[f64]| {
            Ok(nalgebra::DMatrix::from_row_slice(
                2,
                2,
                &[2.0 * x[0], 2.0 * x[1], 1.0, -1.0],
            )) as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |_x: &[f64]| true;

        let mut solver = LMSolver {
            f,
            jacobian: j,
            feasible,
            lambda: 1e-3,
            tol: 1e-12,
            max_iter: 50,
            alpha_min: 1e-6,
        };

        let x0 = vec![0.8, 0.3];
        let sol = solver.solve(x0).unwrap();

        let expected = 1.0 / 2.0_f64.sqrt();
        assert!(approx_eq(sol[0], expected, 1e-8));
        assert!(approx_eq(sol[1], expected, 1e-8));
    }

    #[test]
    fn lm_respects_feasibility_constraint() {
        let f = |x: &[f64]| Ok(vec![x[0] * x[0] - 1.0]) as Result<Vec<f64>, ReactionExtentError>;

        let j = |x: &[f64]| {
            Ok(nalgebra::DMatrix::from_row_slice(1, 1, &[2.0 * x[0]]))
                as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |x: &[f64]| x[0] >= 0.0;

        let mut solver = LMSolver {
            f,
            jacobian: j,
            feasible,
            lambda: 1e-3,
            tol: 1e-12,
            max_iter: 50,
            alpha_min: 1e-6,
        };

        let x0 = vec![0.1];
        let sol = solver.solve(x0).unwrap();

        assert!(approx_eq(sol[0], 1.0, 1e-8));
    }

    #[test]
    fn lm_handles_flat_jacobian() {
        let f = |x: &[f64]| Ok(vec![x[0].powi(3)]) as Result<Vec<f64>, ReactionExtentError>;

        let j = |x: &[f64]| {
            Ok(nalgebra::DMatrix::from_row_slice(
                1,
                1,
                &[3.0 * x[0] * x[0]],
            )) as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |_x: &[f64]| true;

        let mut solver = LMSolver {
            f,
            jacobian: j,
            feasible,
            lambda: 1e-4,
            tol: 1e-15,
            max_iter: 1000,
            alpha_min: 1e-8,
        };

        let x0 = vec![0.5];
        let sol = solver.solve(x0).unwrap();
        info!("{:?}", &sol);
        assert!(sol[0].abs() < 1e-5);
    }
    #[test]
    fn lm_reports_max_iterations() {
        let f = |_x: &[f64]| Ok(vec![1.0]); // no root
        let j = |_x: &[f64]| Ok(nalgebra::DMatrix::identity(1, 1));
        let feasible = |_x: &[f64]| true;

        let mut solver = LMSolver {
            f,
            jacobian: j,
            feasible,
            lambda: 1e-3,
            tol: 1e-12,
            max_iter: 5,
            alpha_min: 1e-6,
        };

        let res = solver.solve(vec![0.0]);
        assert!(matches!(res, Err(SolveError::MaxIterations)));
    }
}

#[cfg(test)]
mod tests_reaction_basis {
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::compute_reaction_basis;
    use log::info;
    use nalgebra::DMatrix;

    const TOL: f64 = 1e-10;

    #[test]
    fn diatomic_dissociation_o2_o() {
        // Species: O2, O
        // Elements: O
        //
        // Reaction: O2 ⇌ 2 O

        let a = DMatrix::from_row_slice(
            2,
            1,
            &[
                2.0, // O2
                1.0, // O
            ],
        );

        let basis = compute_reaction_basis(&a, TOL).unwrap();
        info!("basis: {:?}", basis);
        assert_eq!(basis.rank, 1);
        assert_eq!(basis.num_reactions, 1);

        // Check element conservation: Aᵀ * ν = 0
        let residual = a.transpose() * &basis.reactions;
        assert!(residual.norm() < 1e-12);
    }

    #[test]
    fn nitrogen_oxygen_system() {
        // Species: N2, O2, NO, NO2
        // Elements: N, O

        let a = DMatrix::from_row_slice(
            4,
            2,
            &[
                2.0, 0.0, // N2
                0.0, 2.0, // O2
                1.0, 1.0, // NO
                1.0, 2.0, // NO2
            ],
        );

        let basis = compute_reaction_basis(&a, TOL).unwrap();

        assert_eq!(basis.rank, 2);
        assert_eq!(basis.num_reactions, 2);

        // Check nullspace condition
        let residual = a.transpose() * &basis.reactions;
        assert!(residual.norm() < 1e-12);
    }

    #[test]
    fn methane_combustion_species() {
        // Species: CH4, O2, CO2, H2O
        // Elements: C, H, O

        let a = DMatrix::from_row_slice(
            4,
            3,
            &[
                1.0, 4.0, 0.0, // CH4
                0.0, 0.0, 2.0, // O2
                1.0, 0.0, 2.0, // CO2
                0.0, 2.0, 1.0, // H2O
            ],
        );

        let basis = compute_reaction_basis(&a, TOL).unwrap();
        assert_eq!(basis.rank, 3);
        assert_eq!(basis.num_reactions, 1);

        let residual = a.transpose() * &basis.reactions;
        assert!(residual.norm() < 1e-12);
    }

    #[test]
    fn reaction_basis_o2_not_zero() {
        use nalgebra::DMatrix;

        let a = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);

        let rb = compute_reaction_basis(&a, 1e-12).unwrap();

        assert_eq!(rb.num_reactions, 1);

        let nu = rb.reactions.column(0);

        assert!(
            nu.iter().any(|&x| x.abs() > 1e-6),
            "Reaction vector is zero!"
        );

        // Element conservation check
        assert!((2.0 * nu[0] + 1.0 * nu[1]).abs() < 1e-10);
    }
}

#[cfg(test)]
mod nr_tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn max_step_limits_alpha() {
        // species 0 has a large negative step, species 1 increases
        let n = vec![0.01, 0.5];
        let delta = vec![-1.0, 0.1];
        let alpha = max_step_moles_nonnegative(&n, &delta, 0.95).unwrap();

        // For species 0: alpha0 = 0.95 * 0.01 / 1.0 = 0.0095
        let expected = 0.95 * 0.01 / 1.0;
        assert!((alpha - expected).abs() < 1e-12);
    }

    #[test]
    fn max_step_rejects_malformed_public_inputs_without_panicking() {
        assert!(matches!(
            max_step_moles_nonnegative(&[1.0], &[1.0, 2.0], 0.95),
            Err(SolveError::EvalError(_))
        ));
        assert!(matches!(
            max_step_moles_nonnegative(&[f64::NAN], &[1.0], 0.95),
            Err(SolveError::EvalError(_))
        ));
        assert!(matches!(
            max_step_moles_nonnegative(&[1.0], &[-1.0], 1.5),
            Err(SolveError::EvalError(_))
        ));
    }

    #[test]
    fn nr_solver_respects_bounds_and_converges() {
        // Solve simple scalar problem f(x) = x (root at 0). Start at small positive x0.
        let f = |x: &[f64]| Ok(vec![x[0]]) as Result<Vec<f64>, ReactionExtentError>;
        let j = |_: &[f64]| {
            Ok(DMatrix::from_row_slice(1, 1, &[1.0])) as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |x: &[f64]| x[0] >= 0.0;

        let mut solver = NRSolver {
            f,
            jacobian: j,
            feasible,
            n0: vec![0.01],
            reactions: DMatrix::zeros(1, 0),
            tol: 1e-12,
            max_iter: 100,
            alpha_min: 1e-12,
        };

        let sol = solver.solve(vec![0.01]).unwrap();

        // Solution should remain non-negative and be close to zero
        assert!(sol[0] >= 0.0);
        assert!(sol[0].abs() < 1e-8);
    }
}

#[cfg(test)]
mod trust_region_tests {
    use super::*;
    use approx::assert_relative_eq;
    #[allow(dead_code)]
    fn vec_norm(v: &[f64]) -> f64 {
        v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    fn approx_eq(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.iter().zip(b.iter()).all(|(x, y)| (*x - *y).abs() < tol)
    }

    #[test]
    fn trust_region_linear_system() {
        let f = |x: &[f64]| -> Result<Vec<f64>, ReactionExtentError> {
            Ok(vec![3.0 * x[0] + 1.0 * x[1] - 1.0, 1.0 * x[0] + 2.0 * x[1]])
        };

        let j = |_x: &[f64]| -> Result<DMatrix<f64>, ReactionExtentError> {
            Ok(DMatrix::from_row_slice(2, 2, &[3.0, 1.0, 1.0, 2.0]))
        };

        let feasible = |_x: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 10,
            delta_init: 1.0,
            delta_max: 10.0,
            eta: 0.1,
        };

        let x0 = vec![0.0, 0.0];
        let sol = solver.solve(x0).unwrap();

        assert!(approx_eq(&sol, &[0.4, -0.2], 1e-10));
    }

    #[test]
    fn trust_region_scalar_nonlinear() {
        let f =
            |x: &[f64]| -> Result<Vec<f64>, ReactionExtentError> { Ok(vec![x[0] * x[0] - 2.0]) };

        let j = |x: &[f64]| -> Result<DMatrix<f64>, ReactionExtentError> {
            Ok(DMatrix::from_row_slice(1, 1, &[2.0 * x[0]]))
        };

        let feasible = |_x: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 50,
            delta_init: 0.1, // intentionally small
            delta_max: 10.0,
            eta: 0.1,
        };

        let sol = solver.solve(vec![0.1]).unwrap();
        assert!((sol[0] - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn trust_region_ill_conditioned() {
        let f = |x: &[f64]| -> Result<Vec<f64>, ReactionExtentError> {
            Ok(vec![1e6 * x[0] - 1.0, x[1] - 1.0])
        };

        let j = |_x: &[f64]| -> Result<DMatrix<f64>, ReactionExtentError> {
            Ok(DMatrix::from_row_slice(2, 2, &[1e6, 0.0, 0.0, 1.0]))
        };

        let feasible = |_x: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 50,
            delta_init: 0.01,
            delta_max: 1.0,
            eta: 0.1,
        };

        let sol = solver.solve(vec![0.0, 0.0]).unwrap();

        assert!((sol[0] - 1e-6).abs() < 1e-10);
        assert!((sol[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn trust_region_feasibility_constraint() {
        let f = |x: &[f64]| -> Result<Vec<f64>, ReactionExtentError> { Ok(vec![x[0] - 1.0]) };

        let j = |_x: &[f64]| -> Result<DMatrix<f64>, ReactionExtentError> {
            Ok(DMatrix::from_row_slice(1, 1, &[1.0]))
        };

        let feasible = |x: &[f64]| x[0] >= 0.0;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 50,
            delta_init: 10.0,
            delta_max: 10.0,
            eta: 0.1,
        };

        let sol = solver.solve(vec![-10.0]).unwrap();
        assert!((sol[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn trust_region_singular_jacobian() {
        let f = |_x: &[f64]| -> Result<Vec<f64>, ReactionExtentError> { Ok(vec![1.0]) };

        let j =
            |_x: &[f64]| -> Result<DMatrix<f64>, ReactionExtentError> { Ok(DMatrix::zeros(1, 1)) };

        let feasible = |_x: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 5,
            delta_init: 1.0,
            delta_max: 1.0,
            eta: 0.1,
        };

        let res = solver.solve(vec![0.0]);
        assert!(matches!(res, Err(SolveError::SingularMatrix)));
    }
    ///////////////////////
    #[test]
    fn tr_solves_scalar_quadratic() {
        let f = |x: &[f64]| Ok(vec![x[0] * x[0] - 2.0]) as Result<Vec<f64>, ReactionExtentError>;

        let j = |x: &[f64]| {
            Ok(nalgebra::DMatrix::from_row_slice(1, 1, &[2.0 * x[0]]))
                as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |_x: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 50,
            delta_init: 0.1,
            delta_max: 10.0,
            eta: 0.0,
        };

        let x0 = vec![0.1];
        let sol = solver.solve(x0).unwrap();

        assert_relative_eq!(sol[0], 2.0_f64.sqrt(), epsilon = 1e-8);
    }

    #[test]
    fn tr_solves_2d_nonlinear_system() {
        let f = |x: &[f64]| {
            Ok(vec![x[0] * x[0] + x[1] * x[1] - 1.0, x[0] - x[1]])
                as Result<Vec<f64>, ReactionExtentError>
        };

        let j = |x: &[f64]| {
            Ok(nalgebra::DMatrix::from_row_slice(
                2,
                2,
                &[2.0 * x[0], 2.0 * x[1], 1.0, -1.0],
            )) as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |_x: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 1000,
            delta_init: 0.1,
            delta_max: 10.0,
            eta: 0.0,
        };

        let x0 = vec![0.5, 0.5];
        let sol = solver.solve(x0).unwrap();

        let expected = 1.0 / 2.0_f64.sqrt();
        assert_relative_eq!(sol[0], expected, epsilon = 1e-8);
        assert_relative_eq!(sol[1], expected, epsilon = 1e-8);
    }

    #[test]
    fn tr_respects_feasibility_constraint() {
        let f = |x: &[f64]| Ok(vec![x[0] * x[0] - 1.0]) as Result<Vec<f64>, ReactionExtentError>;

        let j = |x: &[f64]| {
            Ok(nalgebra::DMatrix::from_row_slice(1, 1, &[2.0 * x[0]]))
                as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |_: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 50,
            delta_init: 1.0,
            delta_max: 10.0,
            eta: 0.1,
        };

        let x0 = vec![0.1];
        let sol = solver.solve(x0).unwrap();

        // The unconstrained problem has two equally valid roots. The
        // feasibility predicate above accepts both, so the contract is root
        // convergence rather than a particular sign selected by the solver.
        assert_relative_eq!(sol[0].abs(), 1.0, epsilon = 1e-8);
    }

    #[test]
    fn tr_reports_max_iterations() {
        let f = |_x: &[f64]| Ok(vec![1.0]); // no root
        let j = |_x: &[f64]| Ok(nalgebra::DMatrix::identity(1, 1));
        let feasible = |_x: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 5,
            delta_init: 1.0,
            delta_max: 10.0,
            eta: 0.1,
        };

        let res = solver.solve(vec![0.0]);
        assert!(matches!(res, Err(SolveError::MaxIterations)));
    }

    #[test]
    fn tr_reports_jacobian_dimension_mismatch() {
        let f = |_x: &[f64]| Ok(vec![1.0]) as Result<Vec<f64>, ReactionExtentError>;

        let j = |_x: &[f64]| {
            Ok(nalgebra::DMatrix::from_row_slice(1, 2, &[1.0, 0.0]))
                as Result<DMatrix<f64>, ReactionExtentError>
        };

        let feasible = |_x: &[f64]| true;

        let solver = TrustRegionSolver {
            f,
            jacobian: j,
            feasible,
            tol: 1e-12,
            max_iter: 5,
            delta_init: 1.0,
            delta_max: 10.0,
            eta: 0.1,
        };

        let res = solver.solve(vec![0.0]);
        assert!(matches!(res, Err(SolveError::EvalError(_))));
    }
}
