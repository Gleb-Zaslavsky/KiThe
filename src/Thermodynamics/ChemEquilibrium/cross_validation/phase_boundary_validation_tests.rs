//! Typed failure and root-contract tests for the independent pure-phase
//! boundary validator.
//!
//! These tests intentionally avoid the production active-set runner. Their
//! purpose is to keep the scalar `ln(Q)-ln(K)` reference implementation
//! explicit, deterministic, and independently diagnosable.

use std::rc::Rc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionExtentError, SolveError,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    MOLAR_GAS_CONSTANT, PurePhaseBoundaryElementComposition, PurePhaseBoundaryProblem,
    PurePhaseBoundarySolverSettings, PurePhaseBoundaryStructuralTolerances,
    PurePhaseBoundaryTemperatureSearchSettings, PurePhaseBoundaryTolerances,
    PurePhaseBoundaryValidator, bisect_pure_phase_boundary_temperature,
    evaluate_pure_phase_boundary,
};
use nalgebra::DMatrix;

fn constant_gibbs(value: f64) -> GibbsFn {
    Rc::new(move |_| value)
}

/// The gas-only A/B state has `ln(Q)=ln(2)`.  The candidate closure creates a
/// boundary at `boundary_temperature` with a controlled non-zero slope.
fn analytic_problem(
    temperature: f64,
    boundary_temperature: f64,
    slope_per_kelvin: f64,
) -> PurePhaseBoundaryProblem {
    let pressure = 101_325.0;
    let candidate_gibbs: GibbsFn = Rc::new(move |evaluated_temperature| {
        let ln_k = 2.0_f64.ln() + slope_per_kelvin * (boundary_temperature - evaluated_temperature);
        -MOLAR_GAS_CONSTANT * evaluated_temperature * ln_k
    });
    PurePhaseBoundaryProblem::new(
        vec!["A".to_string(), "B".to_string()],
        vec![1.0, 1.0],
        vec![-2.0, 1.0],
        1.0,
        vec![constant_gibbs(0.0), constant_gibbs(0.0)],
        candidate_gibbs,
        EquilibriumConditions::new(temperature, pressure, pressure).unwrap(),
        "S",
    )
    .unwrap()
}

fn problem_with_gibbs(gas_gibbs: GibbsFn, candidate_gibbs: GibbsFn) -> PurePhaseBoundaryProblem {
    problem_with_gibbs_at_temperature(gas_gibbs, candidate_gibbs, 1_000.0)
}

fn problem_with_gibbs_at_temperature(
    gas_gibbs: GibbsFn,
    candidate_gibbs: GibbsFn,
    temperature: f64,
) -> PurePhaseBoundaryProblem {
    PurePhaseBoundaryProblem::new(
        vec!["A".to_string(), "B".to_string()],
        vec![1.0, 1.0],
        vec![-2.0, 1.0],
        1.0,
        vec![gas_gibbs, constant_gibbs(0.0)],
        candidate_gibbs,
        EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
        "S",
    )
    .unwrap()
}

fn problem_for_target_k(target_k: f64) -> PurePhaseBoundaryProblem {
    problem_with_gibbs(
        constant_gibbs(0.0),
        constant_gibbs(-MOLAR_GAS_CONSTANT * 1_000.0 * target_k.ln()),
    )
}

#[test]
fn temperature_root_accepts_each_bracket_endpoint() {
    let tolerances = PurePhaseBoundaryTolerances::default();
    for (lower, upper, expected) in [(500.0, 700.0, 500.0), (500.0, 700.0, 700.0)] {
        let root = bisect_pure_phase_boundary_temperature(
            lower,
            upper,
            PurePhaseBoundaryTemperatureSearchSettings::default(),
            tolerances,
            |temperature| Ok(analytic_problem(temperature, expected, 0.02)),
        )
        .expect("an exact bracket endpoint is a valid boundary root");
        assert_eq!(root.temperature, expected);
        assert_eq!(root.iterations, 0);
        assert!(root.log_residual.abs() <= tolerances.max_abs_boundary_log_residual);
    }
}

#[test]
fn temperature_root_rejects_unbracketed_and_mismatched_factory_cases() {
    let tolerances = PurePhaseBoundaryTolerances::default();
    let no_sign_change = bisect_pure_phase_boundary_temperature(
        500.0,
        550.0,
        PurePhaseBoundaryTemperatureSearchSettings::default(),
        tolerances,
        |temperature| Ok(analytic_problem(temperature, 700.0, 0.02)),
    )
    .expect_err("a bracket without a sign change must not be bisected");
    assert!(matches!(
        no_sign_change,
        ReactionExtentError::ValidationNotApplicable { .. }
    ));

    let mismatched_factory = bisect_pure_phase_boundary_temperature(
        500.0,
        700.0,
        PurePhaseBoundaryTemperatureSearchSettings::default(),
        tolerances,
        |_temperature| Ok(analytic_problem(600.0, 600.0, 0.02)),
    )
    .expect_err("the factory must preserve the requested temperature");
    assert!(format!("{mismatched_factory}").contains("factory returned"));

    let invalid_bracket = bisect_pure_phase_boundary_temperature(
        700.0,
        500.0,
        PurePhaseBoundaryTemperatureSearchSettings::default(),
        tolerances,
        |temperature| Ok(analytic_problem(temperature, 600.0, 0.02)),
    )
    .expect_err("unordered brackets are invalid input");
    assert!(format!("{invalid_bracket}").contains("ordered"));
}

#[test]
fn temperature_root_fails_explicitly_after_iteration_budget() {
    let result = bisect_pure_phase_boundary_temperature(
        500.0,
        700.0,
        PurePhaseBoundaryTemperatureSearchSettings { max_iterations: 1 },
        PurePhaseBoundaryTolerances {
            max_abs_boundary_log_residual: 1e-14,
            ..PurePhaseBoundaryTolerances::default()
        },
        |temperature| Ok(analytic_problem(temperature, 550.0, 0.02)),
    )
    .expect_err("one bisection step cannot resolve this off-midpoint boundary");
    assert!(matches!(
        result,
        ReactionExtentError::SolveError(SolveError::MaxIterations)
    ));
}

#[test]
fn favorable_boundary_without_interior_root_is_not_applicable() {
    // No gas coefficient is negative, so positive reaction extent has no
    // finite gas-positivity endpoint. The scalar validator deliberately does
    // not claim to solve this open-ended boundary problem.
    let problem = PurePhaseBoundaryProblem::new(
        vec!["A".to_string(), "B".to_string()],
        vec![1.0, 1.0],
        vec![1.0, 0.0],
        1.0,
        vec![constant_gibbs(0.0), constant_gibbs(0.0)],
        constant_gibbs(-MOLAR_GAS_CONSTANT * 1_000.0 * 4.0_f64.ln()),
        EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0).unwrap(),
        "S",
    )
    .unwrap();
    let error = PurePhaseBoundaryValidator::default()
        .validate(&problem)
        .expect_err("an unbounded phase-forming direction is outside scalar-validator scope");
    assert!(matches!(
        error,
        ReactionExtentError::ValidationNotApplicable { .. }
    ));
}

#[test]
fn invalid_settings_and_nonfinite_gibbs_return_typed_errors() {
    let valid_problem = analytic_problem(1_000.0, 1_100.0, 0.02);
    let invalid_boundary_tolerance = evaluate_pure_phase_boundary(
        &valid_problem,
        PurePhaseBoundaryTolerances {
            max_abs_boundary_log_residual: f64::NAN,
            ..PurePhaseBoundaryTolerances::default()
        },
    )
    .expect_err("NaN boundary tolerance must be rejected before evaluation");
    assert!(matches!(
        invalid_boundary_tolerance,
        ReactionExtentError::InvalidProblem { .. }
    ));

    let invalid_solver = PurePhaseBoundaryValidator {
        solver_settings: PurePhaseBoundarySolverSettings {
            feasibility_margin: f64::INFINITY,
            ..PurePhaseBoundarySolverSettings::default()
        },
        ..PurePhaseBoundaryValidator::default()
    }
    .validate(&valid_problem)
    .expect_err("infinite feasibility margin must be rejected before scalar solve");
    assert!(matches!(
        invalid_solver,
        ReactionExtentError::InvalidProblem { .. }
    ));

    let bad_gas = problem_with_gibbs(constant_gibbs(f64::NAN), constant_gibbs(0.0));
    let gas_error = evaluate_pure_phase_boundary(&bad_gas, PurePhaseBoundaryTolerances::default())
        .expect_err("non-finite gas Gibbs closure must preserve its component index");
    assert!(matches!(
        gas_error,
        ReactionExtentError::InvalidDG0 {
            species_index: 0,
            temperature,
            ..
        } if temperature == 1_000.0
    ));

    let bad_candidate = problem_with_gibbs(constant_gibbs(0.0), constant_gibbs(f64::INFINITY));
    let candidate_error =
        evaluate_pure_phase_boundary(&bad_candidate, PurePhaseBoundaryTolerances::default())
            .expect_err("non-finite candidate Gibbs closure must be typed too");
    assert!(matches!(
        candidate_error,
        ReactionExtentError::InvalidDG0 {
            species_index: 2,
            temperature,
            ..
        } if temperature == 1_000.0
    ));
}

#[test]
fn p10_1_pt_boundary_and_structural_tolerance_matrix_returns_typed_errors() {
    let invalid_boundary_tolerances = [
        PurePhaseBoundaryTolerances {
            max_abs_boundary_log_residual: 0.0,
            ..PurePhaseBoundaryTolerances::default()
        },
        PurePhaseBoundaryTolerances {
            max_abs_equilibrium_log_residual: -1.0,
            ..PurePhaseBoundaryTolerances::default()
        },
        PurePhaseBoundaryTolerances {
            max_abs_equilibrium_log_residual: f64::INFINITY,
            ..PurePhaseBoundaryTolerances::default()
        },
    ];
    for tolerances in invalid_boundary_tolerances {
        let error =
            evaluate_pure_phase_boundary(&analytic_problem(1_000.0, 1_100.0, 0.02), tolerances)
                .expect_err("invalid P,T boundary tolerances must never reach the evaluator");
        assert!(matches!(error, ReactionExtentError::InvalidProblem { .. }));
    }

    let composition = PurePhaseBoundaryElementComposition::new(
        vec!["X".to_string(), "Y".to_string()],
        DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]),
        vec![2.0, 1.0],
    )
    .expect("controlled composition must be valid");
    let invalid_structural_tolerances = [
        PurePhaseBoundaryStructuralTolerances {
            max_abs_element_balance: 0.0,
            ..PurePhaseBoundaryStructuralTolerances::default()
        },
        PurePhaseBoundaryStructuralTolerances {
            rank_absolute_tolerance: -1.0,
            ..PurePhaseBoundaryStructuralTolerances::default()
        },
        PurePhaseBoundaryStructuralTolerances {
            rank_absolute_tolerance: f64::NAN,
            ..PurePhaseBoundaryStructuralTolerances::default()
        },
        PurePhaseBoundaryStructuralTolerances {
            rank_relative_tolerance: f64::INFINITY,
            ..PurePhaseBoundaryStructuralTolerances::default()
        },
    ];
    for tolerances in invalid_structural_tolerances {
        let error = match analytic_problem(1_000.0, 1_100.0, 0.02)
            .with_element_composition(composition.clone(), tolerances)
        {
            Ok(_) => {
                panic!("invalid P,T structural tolerances must fail during problem construction")
            }
            Err(error) => error,
        };
        assert!(matches!(error, ReactionExtentError::InvalidProblem { .. }));
    }
}

#[test]
fn finite_extent_root_is_robust_near_each_feasibility_boundary() {
    let validator = PurePhaseBoundaryValidator::default();
    let upper = problem_for_target_k(4.0)
        .positive_extent_upper_bound()
        .unwrap();

    // K only slightly exceeds Q(0)=2, so the physical pure-phase amount is
    // close to zero but still strictly interior.
    let near_zero = validator
        .validate(&problem_for_target_k(2.0 * 1.000_1))
        .unwrap()
        .equilibrium
        .expect("slightly favorable candidate must have an interior root");
    assert!(near_zero.extent > 0.0 && near_zero.extent < 1e-3);
    assert!(near_zero.gas_moles.iter().all(|moles| *moles > 0.0));
    assert!(near_zero.log_residual.abs() <= validator.tolerances.max_abs_equilibrium_log_residual);

    // A very large K moves the root near A(g) exhaustion. The validator must
    // retain strict positivity rather than evaluating an activity at zero.
    let near_upper = validator
        .validate(&problem_for_target_k(1e12))
        .unwrap()
        .equilibrium
        .expect("strongly favorable candidate must still have an interior root");
    assert!(near_upper.extent > 0.0 && near_upper.extent < upper);
    assert!(upper - near_upper.extent < 1e-4);
    assert!(near_upper.gas_moles.iter().all(|moles| *moles > 0.0));
    assert!(near_upper.log_residual.abs() <= validator.tolerances.max_abs_equilibrium_log_residual);
}

#[test]
fn temperature_factory_nonfinite_gibbs_preserves_typed_component_context() {
    let error = bisect_pure_phase_boundary_temperature(
        500.0,
        700.0,
        PurePhaseBoundaryTemperatureSearchSettings::default(),
        PurePhaseBoundaryTolerances::default(),
        |temperature| {
            Ok(problem_with_gibbs_at_temperature(
                constant_gibbs(0.0),
                constant_gibbs(if temperature == 500.0 { f64::NAN } else { 0.0 }),
                temperature,
            ))
        },
    )
    .expect_err("a non-finite factory Gibbs value must stop root evaluation");
    assert!(matches!(
        error,
        ReactionExtentError::InvalidDG0 {
            species_index: 2,
            temperature,
            ..
        } if temperature == 500.0
    ));
}
