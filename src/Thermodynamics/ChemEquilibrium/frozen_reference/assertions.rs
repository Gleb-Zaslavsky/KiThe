//! Narrow shared assertions for frozen-reference accepted-state regressions.
//!
//! The helpers own only contracts that every bounded multiphase result must
//! satisfy. They deliberately do not encode fixture topology, source-specific
//! envelopes, or a particular reference answer.

use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;

/// Numerical limits applied to a published bounded phase-control solution.
#[derive(Debug, Clone, Copy)]
pub(crate) struct AcceptedSolutionContract {
    /// Largest accepted canonical residual L2 norm.
    pub(crate) max_residual_l2_norm: f64,
    /// Largest accepted element-balance residual.
    pub(crate) max_abs_element_balance_error: f64,
}

impl AcceptedSolutionContract {
    /// Constructs a finite, positive validation envelope.
    pub(crate) fn new(max_residual_l2_norm: f64, max_abs_element_balance_error: f64) -> Self {
        assert!(
            max_residual_l2_norm.is_finite() && max_residual_l2_norm > 0.0,
            "accepted-solution residual limit must be finite and positive"
        );
        assert!(
            max_abs_element_balance_error.is_finite() && max_abs_element_balance_error > 0.0,
            "accepted-solution balance limit must be finite and positive"
        );
        Self {
            max_residual_l2_norm,
            max_abs_element_balance_error,
        }
    }
}

/// Verifies the invariant portion of a bounded accepted result.
///
/// The physical view must contain no negative or non-finite amounts. The
/// canonical validation must meet the caller-owned numerical limits, and a
/// bounded lifecycle must have published satisfied complementarity evidence.
pub(crate) fn assert_bounded_solution_accepted(
    solution: &MultiphaseEquilibriumSolution,
    contract: AcceptedSolutionContract,
) {
    assert!(
        solution
            .component_moles()
            .iter()
            .all(|moles| moles.is_finite() && *moles >= 0.0),
        "published physical component moles must be finite and non-negative"
    );
    let validation = solution.accepted_solution().validation();
    assert!(
        validation.residual_l2_norm.is_finite()
            && validation.residual_l2_norm <= contract.max_residual_l2_norm,
        "accepted residual {} exceeds frozen-reference limit {}",
        validation.residual_l2_norm,
        contract.max_residual_l2_norm,
    );
    assert!(
        validation.max_abs_element_balance_error.is_finite()
            && validation.max_abs_element_balance_error <= contract.max_abs_element_balance_error,
        "accepted element balance {} exceeds frozen-reference limit {}",
        validation.max_abs_element_balance_error,
        contract.max_abs_element_balance_error,
    );
    assert!(
        solution
            .acceptance_report()
            .is_some_and(|report| report.complementarity.satisfied),
        "bounded phase-control result must publish satisfied complementarity evidence"
    );
}

/// Confirms that a read-only fixture or local-library byte snapshot did not drift.
pub(crate) fn assert_byte_snapshots_unchanged(
    before: &[(String, Vec<u8>)],
    after: &[(String, Vec<u8>)],
) {
    assert_eq!(
        before, after,
        "frozen/local source bytes changed during a read-only test"
    );
}
