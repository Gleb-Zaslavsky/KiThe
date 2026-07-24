//! Tests for backend-independent equilibrium candidate validation.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::{
    EquilibriumAcceptanceCriteria, EquilibriumCandidateResiduals, compare_candidate_reports,
    select_preferred_candidate_index, validate_equilibrium_candidate,
};
use nalgebra::DMatrix;

#[test]
fn candidate_validation_reports_finite_conserved_solution() {
    let report = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[0.0, 0.0],
            raw_residual: &[1e-8, -1e-8],
            acceptance_residual: &[1e-8, -1e-8],
        },
        EquilibriumAcceptanceCriteria::new(1e-6, 1e-12, 1e-6).unwrap(),
        &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
        &[2.0],
    )
    .unwrap();

    assert!(report.residual_l2_norm < 1e-6);
    assert!((report.residual_rms - 1e-8).abs() < 1e-12);
    assert_eq!(report.raw_residual_l2_norm, report.residual_l2_norm);
    assert_eq!(report.raw_residual_rms, report.residual_rms);
    assert_eq!(report.max_abs_element_balance_error, 0.0);
    assert_eq!(report.reaction_affinity_l2_norm, 1e-8);
    assert_eq!(report.max_abs_reaction_affinity, 1e-8);
    assert_eq!(report.min_moles, 1.0);
}

#[test]
fn candidate_validation_rejects_backend_success_with_large_residual() {
    let result = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[0.0],
            raw_residual: &[1e-2],
            acceptance_residual: &[1e-2],
        },
        EquilibriumAcceptanceCriteria::new(1e-6, 1e-12, 1e-6).unwrap(),
        &DMatrix::from_row_slice(1, 1, &[1.0]),
        &[1.0],
    );

    let err = result.unwrap_err();
    match err {
        ReactionExtentError::InvalidCandidate { field, message } => {
            assert_eq!(field, "candidate_acceptance_residual");
            assert!(message.contains("exceeds tolerance"));
            let rendered = ReactionExtentError::InvalidCandidate {
                field,
                message: message.clone(),
            }
            .to_string();
            assert!(rendered.contains("candidate_acceptance_residual"));
        }
        other => panic!("expected invalid candidate error, got {other:?}"),
    }
}

#[test]
fn candidate_validation_keeps_raw_and_scaled_residuals_distinct() {
    let report = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[0.0, 0.0],
            raw_residual: &[100.0, 0.0],
            acceptance_residual: &[0.5, 0.0],
        },
        EquilibriumAcceptanceCriteria::new(1.0, 1e-12, 1.0).unwrap(),
        &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
        &[2.0],
    )
    .unwrap();

    assert_eq!(report.residual_l2_norm, 0.5);
    assert_eq!(report.residual_rms, 0.5 / 2.0_f64.sqrt());
    assert_eq!(report.raw_residual_l2_norm, 100.0);
    assert_eq!(report.raw_residual_rms, 100.0 / 2.0_f64.sqrt());
    assert_eq!(report.raw_max_abs_residual, 100.0);
    assert_eq!(report.reaction_affinity_l2_norm, 0.5);
    assert_eq!(report.max_abs_reaction_affinity, 0.5);
}

#[test]
fn candidate_validation_rejects_non_finite_log_moles() {
    let result = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[f64::INFINITY],
            raw_residual: &[0.0],
            acceptance_residual: &[0.0],
        },
        EquilibriumAcceptanceCriteria::new(1e-6, 1e-12, 1e-6).unwrap(),
        &DMatrix::from_row_slice(1, 1, &[1.0]),
        &[1.0],
    );

    assert!(matches!(
        result,
        Err(ReactionExtentError::InvalidCandidate { .. })
    ));
}

#[test]
fn candidate_validation_rejects_truncated_or_empty_residual_vectors() {
    let matrix = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);

    let truncated_raw = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[0.0, 0.0],
            raw_residual: &[0.0],
            acceptance_residual: &[0.0, 0.0],
        },
        EquilibriumAcceptanceCriteria::new(1e-6, 1e-12, 1e-6).unwrap(),
        &matrix,
        &[2.0],
    );
    assert!(matches!(
        truncated_raw,
        Err(ReactionExtentError::DimensionMismatch(_))
    ));

    let empty_acceptance = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[0.0, 0.0],
            raw_residual: &[0.0, 0.0],
            acceptance_residual: &[],
        },
        EquilibriumAcceptanceCriteria::new(1e-6, 1e-12, 1e-6).unwrap(),
        &matrix,
        &[2.0],
    );
    assert!(matches!(
        empty_acceptance,
        Err(ReactionExtentError::DimensionMismatch(_))
    ));
}

#[test]
fn acceptance_criteria_rejects_invalid_tolerances_before_validation() {
    let err = EquilibriumAcceptanceCriteria::new(-1.0, 1e-12, 1e-6).unwrap_err();
    assert!(matches!(
        err,
        ReactionExtentError::InvalidProblem {
            field: "candidate_tolerances",
            ..
        }
    ));
}

#[test]
fn candidate_selection_prefers_the_more_accurate_report_and_uses_ties_deterministically() {
    let stronger = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[0.0, 0.0],
            raw_residual: &[1e-9, 1e-9],
            acceptance_residual: &[1e-9, 1e-9],
        },
        EquilibriumAcceptanceCriteria::new(1e-6, 1e-12, 1e-6).unwrap(),
        &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
        &[2.0],
    )
    .unwrap();

    let weaker = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[0.0, 0.0],
            raw_residual: &[1e-6, 1e-6],
            acceptance_residual: &[1e-6, 1e-6],
        },
        EquilibriumAcceptanceCriteria::new(1e-3, 1e-12, 1e-3).unwrap(),
        &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
        &[2.0],
    )
    .unwrap();

    assert_eq!(
        compare_candidate_reports(&stronger, &weaker),
        std::cmp::Ordering::Less
    );
    assert_eq!(
        select_preferred_candidate_index(&[weaker.clone(), stronger.clone()]),
        Some(1)
    );
    assert_eq!(
        select_preferred_candidate_index(&[stronger.clone(), stronger.clone()]),
        Some(0)
    );
}

#[test]
fn candidate_validation_rejects_large_reaction_affinity_even_when_residual_is_small() {
    let result = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &[0.0, 0.0],
            raw_residual: &[1e-4, 1e-4],
            acceptance_residual: &[1.0, 0.0],
        },
        EquilibriumAcceptanceCriteria::new(10.0, 1e-12, 0.1).unwrap(),
        &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
        &[2.0],
    );

    match result.unwrap_err() {
        ReactionExtentError::InvalidCandidate { field, message } => {
            assert_eq!(field, "candidate_reaction_affinity");
            assert!(message.contains("exceeds tolerance"));
        }
        other => panic!("expected reaction-affinity rejection, got {other:?}"),
    }
}
