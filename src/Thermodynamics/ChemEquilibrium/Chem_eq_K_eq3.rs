//! Numerical solvers and supporting structures for chemical equilibrium calculations.
//!
//! # Purpose
//! This module provides robust numerical solvers for nonlinear equation systems arising
//! in chemical equilibrium problems. It includes Levenberg-Marquardt and Newton-Raphson
//! solvers with line search, feasibility constraints, and reaction basis computation.
//!
//! # Main Structures
//! - [`LMSolver`]: Levenberg-Marquardt solver with damping for robust convergence
//! - [`NRSolver`]: Newton-Raphson solver with line search and bound constraints
//! - [`TrustRegionSolver`] - Trust Region solver
//! - [`ReactionBasis`]: Container for independent reaction stoichiometry from SVD
//! - [`SolveError`]: Error types for solver failures
//! - [`ReactionExtentError`]: Comprehensive error handling for equilibrium calculations
//!
//! # Key Methods
//! - [`LMSolver::solve`]: Solves nonlinear system using Levenberg-Marquardt algorithm
//! - [`NRSolver::solve`]: Solves nonlinear system using Newton-Raphson with line search
//! - [`compute_reaction_basis`]: Computes independent reactions from elemental composition
//! - [`max_step_moles_nonnegative`]: Ensures mole numbers remain non-negative
//!
//! # Examples
//! ```rust
//! use KiThe::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq2::*;
//!
//! // Create Levenberg-Marquardt solver
//! let mut solver = LMSolver {
//!     f: residual_function,
//!     jacobian: jacobian_function,
//!     feasible: |x| x.iter().all(|&xi| xi > 0.0),
//!     lambda: 1e-3,
//!     tol: 1e-12,
//!     max_iter: 100,
//!     alpha_min: 1e-6,
//! };
//!
//! let solution = solver.solve(initial_guess)?;
//! ```
//!
//! # Non-obvious Solutions and Tips
//! - LM solver uses adaptive damping that increases when steps are rejected
//! - NR solver includes bound-aware step limiting to prevent negative mole numbers
//! - SVD-based reaction basis ensures numerical stability and element conservation
//! - Line search with feasibility checking prevents physically invalid solutions
//! - Gram-Schmidt completion handles cases where SVD doesn't return full nullspace
use crate::Thermodynamics::User_substances_error::SubsDataError;
use log::{error, info};
use nalgebra::{DMatrix, DVector};

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

/// Comprehensive error types for chemical equilibrium calculations
#[derive(Debug)]
pub enum ReactionExtentError {
    /// Error from substance database operations
    SubsDataError(SubsDataError),
    /// Error originating from the underlying nonlinear solver
    SolveError(SolveError),
    /// SVD failed when computing reaction basis
    SVDError(String),
    /// Dimension mismatch between vectors/matrices
    DimensionMismatch(String),
    /// Duplicate species assigned to multiple phases
    DuplicateSpecies(usize),
    /// A species is not assigned to any phase
    SpeciesNotAssigned,
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
    /// Invalid DG0 or temperature leading to bad ln_k
    InvalidDG0 { dg0: f64, temperature: f64 },
    /// Other error
    Other(String),
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
            let f_val = (self.f)(&x)
                .map_err(|e| SolveError::EvalError(format!("Residual eval error: {:?}", e)))?;
            info!("Iteration {}, x = {:?}, f = {:?}", _iter, x, f_val);
            let f_norm = l2_norm(&f_val);
            info!("  ||f|| = {}", f_norm);
            if f_norm < self.tol {
                return Ok(x);
            }
            if f_norm.is_nan() || f_norm.is_infinite() {
                return Err(SolveError::SingularMatrix);
            }

            let j = (self.jacobian)(&x)
                .map_err(|e| SolveError::EvalError(format!("Jacobian eval error: {:?}", e)))?;
            info!("  Jacobian:\n{}", j);
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

                let f_trial = (self.f)(&x_trial)
                    .map_err(|e| SolveError::EvalError(format!("Residual eval error: {:?}", e)))?;
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
) -> f64 {
    assert_eq!(n.len(), delta_n.len());

    let mut alpha_max: f64 = 1.0;

    for i in 0..n.len() {
        let ni = n[i];
        let dni = delta_n[i];

        // Only decreasing species can violate positivity
        if dni < 0.0 {
            if ni <= 0.0 {
                return 0.0; // already infeasible
            }

            let alpha_i = safety * ni / (-dni);
            alpha_max = alpha_max.min(alpha_i);
        }
    }

    alpha_max.clamp(0.0, 1.0)
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
            let f_val = (self.f)(&x).map_err(|e| SolveError::EvalError(format!("{:?}", e)))?;
            info!("NR Iteration {}, x = {:?}, f = {:?}", iter, x, f_val);
            let f_norm = l2_norm(&f_val);
            info!(" ||f|| = {}", f_norm);
            if f_norm < self.tol {
                return Ok(x);
            }
            if f_norm.is_nan() || f_norm.is_infinite() {
                return Err(SolveError::SingularMatrix);
            }
            let j = (self.jacobian)(&x).map_err(|e| SolveError::EvalError(format!("{:?}", e)))?;

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
                    (self.f)(&x_trial).map_err(|e| SolveError::EvalError(format!("{:?}", e)))?;

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
                return Err(SolveError::MaxIterations);
            }
        }

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

        for iter in 0..self.max_iter {
            let f_val = (self.f)(&x).map_err(|e| SolveError::EvalError(format!("{:?}", e)))?;

            let f_norm = l2_norm(&f_val);
            if f_norm < self.tol {
                return Ok(x);
            }

            let J = (self.jacobian)(&x).map_err(|e| SolveError::EvalError(format!("{:?}", e)))?;

            if J.nrows() != n || J.ncols() != n {
                return Err(SolveError::EvalError("Jacobian not square".into()));
            }

            let fvec = DVector::from_vec(f_val.clone());

            // --- Newton step ---
            let newton_step = J
                .clone()
                .lu()
                .solve(&(-&fvec))
                .unwrap_or_else(|| DVector::zeros(n));

            // --- Cauchy step ---
            let g = &J.transpose() * &fvec;
            let g_norm_sq = g.dot(&g);
            let jg = &J * &g;
            let alpha = g_norm_sq / jg.dot(&jg).max(1e-16);
            let cauchy_step = -alpha * g;

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
                (self.f)(&x_trial).map_err(|e| SolveError::EvalError(format!("{:?}", e)))?;

            let actual_reduction = f_norm.powi(2) - l2_norm(&f_trial).powi(2);

            let model_reduction = -2.0 * fvec.dot(&p) - p.dot(&(J.transpose() * &J * &p));

            let rho = actual_reduction / model_reduction.max(1e-16);

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
                return Err(SolveError::SingularMatrix);
            }
        }

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
    use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq3::compute_reaction_basis;
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
        let alpha = max_step_moles_nonnegative(&n, &delta, 0.95);

        // For species 0: alpha0 = 0.95 * 0.01 / 1.0 = 0.0095
        let expected = 0.95 * 0.01 / 1.0;
        assert!((alpha - expected).abs() < 1e-12);
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
        assert!(res.is_err());
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

        let feasible = |x: &[f64]| true;

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

        assert_relative_eq!(sol[0], 1.0, epsilon = 1e-8);
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
}
