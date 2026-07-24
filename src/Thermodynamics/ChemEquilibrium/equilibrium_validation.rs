//! Backend-independent acceptance gate for equilibrium solver candidates.
//!
//! # Purpose
//!
//! A nonlinear solver backend may terminate successfully (reporting convergence)
//! while returning a solution that is **non-finite, elementally unconserved, or
//! insufficiently converged**. This module implements the **acceptance gate**
//! that validates candidates independently of which backend produced them.
//!
//! This separation ensures that:
//! - All backends (Legacy LM/NR/TR, RustedSciThe) share the same acceptance criteria.
//! - The acceptance logic is testable without running any solver.
//! - Adding a new backend does not require re-implementing acceptance checks.
//!
//! # Physical and Mathematical Background
//!
//! The acceptance gate checks three conditions:
//!
//! 1. **Residual norm**: `||F(y)||₂ < residual_tolerance` — the nonlinear residual
//!    (element conservation + reaction affinities) must be sufficiently small.
//!
//! 2. **Element balance**: `||A^T·n(y) - b_0||₂ < element_balance_tolerance` —
//!    the solution must conserve elements to within the specified tolerance.
//!    This is an independent check recomputed from physical moles, not from the
//!    residual vector.
//!
//! 3. **Reaction affinity**: `||ν^T·μ(y)/(R·T)||₂ < reaction_affinity_tolerance` —
//!    the reaction affinities must be near zero, indicating thermodynamic equilibrium.
//!
//! Additionally, the gate rejects candidates with:
//! - Non-finite log-mole values (`NaN` or `±Inf`).
//! - Non-finite or negative physical moles after exponentiation.
//! - Dimension mismatches between any of the vectors/matrices.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`EquilibriumCandidateResiduals`] | Typed bundle of log-moles + raw + acceptance residuals |
//! | [`EquilibriumAcceptanceCriteria`] | Three tolerances for the acceptance gate |
//! | [`EquilibriumCandidateReport`] | Full diagnostic report for one candidate |
//!
//! # Key Functions
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`validate_equilibrium_candidate`] | Main acceptance gate — validates one candidate |
//! | [`compare_candidate_reports`] | Compares two candidates to select the better one |
//! | [`select_preferred_candidate_index`] | Selects the best candidate from a list |
//!
//! # Dataflow
//!
//! ```text
//!   Solver backend returns candidate solution (log_moles)
//!     │
//!     v
//!   validate_equilibrium_candidate(
//!     residuals, criteria, elem_comp, n0, stoich, gibbs, T
//!   )
//!     ├── Check log_moles are finite
//!     ├── Compute physical moles: n = exp(log_moles)
//!     ├── Check moles are finite and non-negative
//!     ├── Compute residual norms
//!     ├── Compute element balance: ||A^T·n - b0||
//!     ├── Compute reaction affinities: ||ν^T·μ/(R·T)||
//!     └── Return EquilibriumCandidateReport or error
//!     │
//!     v
//!   EquilibriumLogMoles::solve_backend_cascade()
//!     ├── If accepted: publish candidate
//!     └── If rejected: try next backend or fail
//! ```
//!
//! # Examples
//!
//! ```rust, ignore
//! use equilibrium_validation::*;
//!
//! let criteria = EquilibriumAcceptanceCriteria::new(1e-6, 1e-8, 1e-6).unwrap();
//! let residuals = EquilibriumCandidateResiduals {
//!     log_moles: &log_moles,
//!     raw_residual: &raw_res,
//!     acceptance_residual: &acc_res,
//! };
//! let report = validate_equilibrium_candidate(
//!     &residuals, &criteria, &elem_comp, &n0, &stoich, &gibbs, T
//! ).unwrap();
//! ```
//!
//! # Non-obvious Details
//!
//! - The **raw residual** and **acceptance residual** may differ when row scaling
//!   is applied. The acceptance gate checks both: the raw residual ensures the
//!   unscaled problem is solved, while the acceptance residual ensures the scaled
//!   problem (seen by the backend) is solved.
//! - `compare_candidate_reports` uses a lexicographic ordering: lower residual norm
//!   is preferred, then lower element balance error, then lower reaction affinity.
//! - The `EquilibriumCandidateReport` contains both `L2` and `RMS` norms, plus
//!   `max_abs` values, providing multiple views of convergence quality.
//!
//! # Related Modules
//!
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — orchestrator that calls the acceptance gate
//! - [`equilibrium_solver_policy`](super::equilibrium_solver_policy) — backend selection and cascade
//!

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::compute_species_moles;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use nalgebra::DMatrix;
use std::cmp::Ordering;

/// Typed residual data for one candidate accepted by a nonlinear backend.
///
/// Keeping the log-moles and both residual views together prevents call sites
/// from re-ordering or partially dropping one of the vectors by accident.
#[derive(Debug, Clone, Copy)]
pub struct EquilibriumCandidateResiduals<'a> {
    /// Candidate nonlinear coordinates in solver order.
    pub log_moles: &'a [f64],
    /// Residual before any optional acceptance scaling.
    pub raw_residual: &'a [f64],
    /// Residual used for the backend-independent acceptance gate.
    pub acceptance_residual: &'a [f64],
}

/// Typed acceptance thresholds for one equilibrium candidate.
///
/// The validation gate needs two different tolerances: one for the residual
/// norm seen by the backend, and one for the separately recomputed elemental
/// balance. Keeping them in a single value makes the public contract harder to
/// misuse than passing a pair of anonymous floats around.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumAcceptanceCriteria {
    /// Tolerance for the backend-facing residual norm.
    pub residual_tolerance: f64,
    /// Tolerance for the independently recomputed element balance.
    pub element_balance_tolerance: f64,
    /// Tolerance for the reaction-affinity block.
    pub reaction_affinity_tolerance: f64,
}

impl EquilibriumAcceptanceCriteria {
    /// Builds validated acceptance criteria from two finite, non-negative tolerances.
    pub fn new(
        residual_tolerance: f64,
        element_balance_tolerance: f64,
        reaction_affinity_tolerance: f64,
    ) -> Result<Self, ReactionExtentError> {
        if !residual_tolerance.is_finite()
            || residual_tolerance < 0.0
            || !element_balance_tolerance.is_finite()
            || element_balance_tolerance < 0.0
            || !reaction_affinity_tolerance.is_finite()
            || reaction_affinity_tolerance < 0.0
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "candidate_tolerances",
                message: "candidate tolerances must be finite and non-negative".to_string(),
            });
        }
        Ok(Self {
            residual_tolerance,
            element_balance_tolerance,
            reaction_affinity_tolerance,
        })
    }
}

/// Numerical evidence collected before publishing an equilibrium solution.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumCandidateReport {
    /// Euclidean norm of the residual vector used for acceptance.
    ///
    /// This is the scaled residual norm when the caller enabled residual
    /// scaling, and is therefore the value compared with the tolerance.
    pub residual_l2_norm: f64,
    /// Root-mean-square of the acceptance residual.
    pub residual_rms: f64,
    /// Largest absolute component of the acceptance residual.
    pub max_abs_residual: f64,
    /// Euclidean norm of the unscaled physical formulation residual.
    pub raw_residual_l2_norm: f64,
    /// Root-mean-square of the unscaled physical formulation residual.
    pub raw_residual_rms: f64,
    /// Largest absolute component of the unscaled residual.
    pub raw_max_abs_residual: f64,
    /// Largest absolute elemental-balance error.
    pub max_abs_element_balance_error: f64,
    /// Euclidean norm of the reaction-affinity block of the acceptance residual.
    pub reaction_affinity_l2_norm: f64,
    /// Largest absolute reaction-affinity component.
    pub max_abs_reaction_affinity: f64,
    /// Minimum reconstructed physical mole amount.
    pub min_moles: f64,
}

/// Compares two validated candidates using an explicit numerical criterion.
///
/// The preferred candidate is the one with:
/// 1. smaller backend-facing residual norm,
/// 2. smaller maximum residual component,
/// 3. smaller unscaled residual norm,
/// 4. smaller unscaled maximum residual component,
/// 5. smaller maximum elemental-balance error,
/// 6. larger minimum mole number.
///
/// If all evidence matches, the candidates are considered equivalent and the
/// caller may keep the earlier one for deterministic stability.
pub fn compare_candidate_reports(
    lhs: &EquilibriumCandidateReport,
    rhs: &EquilibriumCandidateReport,
) -> Ordering {
    lhs.residual_l2_norm
        .total_cmp(&rhs.residual_l2_norm)
        .then_with(|| lhs.max_abs_residual.total_cmp(&rhs.max_abs_residual))
        .then_with(|| {
            lhs.raw_residual_l2_norm
                .total_cmp(&rhs.raw_residual_l2_norm)
        })
        .then_with(|| {
            lhs.raw_max_abs_residual
                .total_cmp(&rhs.raw_max_abs_residual)
        })
        .then_with(|| {
            lhs.max_abs_element_balance_error
                .total_cmp(&rhs.max_abs_element_balance_error)
        })
        .then_with(|| rhs.min_moles.total_cmp(&lhs.min_moles))
}

/// Selects the best candidate index from a retained set of validated results.
///
/// The first candidate is kept on ties, which makes the choice deterministic
/// and stable for callers that already retain attempts in execution order.
pub fn select_preferred_candidate_index(
    candidates: &[EquilibriumCandidateReport],
) -> Option<usize> {
    let mut best_index = None;
    for (index, candidate) in candidates.iter().enumerate() {
        match best_index {
            None => best_index = Some(index),
            Some(best) => {
                if compare_candidate_reports(candidate, &candidates[best]) == Ordering::Less {
                    best_index = Some(index);
                }
            }
        }
    }
    best_index
}

/// Checks a log-mole candidate against numerical and conservation invariants.
///
/// `residual` must use the same scaling and tolerance convention as the
/// selected nonlinear backend. Element balances are always recomputed from
/// physical moles, independently of that residual scaling.
pub fn validate_equilibrium_candidate(
    residuals: EquilibriumCandidateResiduals<'_>,
    criteria: EquilibriumAcceptanceCriteria,
    element_composition: &DMatrix<f64>,
    element_totals: &[f64],
) -> Result<EquilibriumCandidateReport, ReactionExtentError> {
    if element_composition.nrows() != residuals.log_moles.len()
        || element_composition.ncols() != element_totals.len()
    {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "candidate has {} species, element matrix is {}x{}, and totals have {} entries",
            residuals.log_moles.len(),
            element_composition.nrows(),
            element_composition.ncols(),
            element_totals.len(),
        )));
    }
    // The canonical log-mole formulation is square. Accepting a truncated or
    // empty residual vector would make its norm look artificially perfect and
    // could publish a candidate that no backend actually validated.
    if residuals.raw_residual.len() != residuals.log_moles.len()
        || residuals.acceptance_residual.len() != residuals.log_moles.len()
    {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "candidate has {} log-moles, raw residual has {}, and acceptance residual has {} entries",
            residuals.log_moles.len(),
            residuals.raw_residual.len(),
            residuals.acceptance_residual.len(),
        )));
    }
    if residuals.log_moles.iter().any(|value| !value.is_finite()) {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "candidate_log_moles",
            message: "candidate contains a non-finite log-mole value".to_string(),
        });
    }
    if residuals
        .raw_residual
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "candidate_raw_residual",
            message: "candidate raw residual contains a non-finite value".to_string(),
        });
    }
    if residuals
        .acceptance_residual
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "candidate_acceptance_residual",
            message: "candidate acceptance residual contains a non-finite value".to_string(),
        });
    }
    if element_totals.iter().any(|value| !value.is_finite()) {
        return Err(ReactionExtentError::InvalidProblem {
            field: "element_totals",
            message: "element totals must be finite".to_string(),
        });
    }

    let moles = compute_species_moles(residuals.log_moles)?;
    let reaction_count = residuals
        .log_moles
        .len()
        .checked_sub(element_composition.ncols())
        .ok_or_else(|| {
            ReactionExtentError::DimensionMismatch(format!(
                "candidate has {} species but {} element columns; the canonical system is not square",
                residuals.log_moles.len(),
                element_composition.ncols(),
            ))
        })?;
    if reaction_count > residuals.acceptance_residual.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "candidate reaction block has {} rows but acceptance residual has {} entries",
            reaction_count,
            residuals.acceptance_residual.len(),
        )));
    }
    let reaction_affinity_block = &residuals.acceptance_residual[..reaction_count];

    let residual_l2_norm = residuals
        .acceptance_residual
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    let residual_rms = residual_l2_norm / (residuals.acceptance_residual.len() as f64).sqrt();
    let max_abs_residual = residuals
        .acceptance_residual
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    if residual_l2_norm > criteria.residual_tolerance {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "candidate_acceptance_residual",
            message: format!(
                "residual L2 norm {residual_l2_norm:e} exceeds tolerance {:e}",
                criteria.residual_tolerance
            ),
        });
    }

    let raw_residual_l2_norm = residuals
        .raw_residual
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    let raw_residual_rms = raw_residual_l2_norm / (residuals.raw_residual.len() as f64).sqrt();
    let raw_max_abs_residual = residuals
        .raw_residual
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);

    let mut max_abs_element_balance_error = 0.0_f64;
    for element in 0..element_composition.ncols() {
        let reconstructed_total = (0..moles.len())
            .map(|species| element_composition[(species, element)] * moles[species])
            .sum::<f64>();
        let error = (reconstructed_total - element_totals[element]).abs();
        if !error.is_finite() || error > criteria.element_balance_tolerance {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "candidate_element_balance",
                message: format!(
                    "element {element} balance error {error:e} exceeds tolerance {:e}",
                    criteria.element_balance_tolerance
                ),
            });
        }
        max_abs_element_balance_error = max_abs_element_balance_error.max(error);
    }

    let reaction_affinity_l2_norm = reaction_affinity_block
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    let max_abs_reaction_affinity = reaction_affinity_block
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    if reaction_affinity_l2_norm > criteria.reaction_affinity_tolerance {
        return Err(ReactionExtentError::InvalidCandidate {
            field: "candidate_reaction_affinity",
            message: format!(
                "reaction-affinity L2 norm {reaction_affinity_l2_norm:e} exceeds tolerance {:e}",
                criteria.reaction_affinity_tolerance
            ),
        });
    }

    Ok(EquilibriumCandidateReport {
        residual_l2_norm,
        residual_rms,
        max_abs_residual,
        raw_residual_l2_norm,
        raw_residual_rms,
        raw_max_abs_residual,
        max_abs_element_balance_error,
        reaction_affinity_l2_norm,
        max_abs_reaction_affinity,
        min_moles: moles.into_iter().fold(f64::INFINITY, f64::min),
    })
}
