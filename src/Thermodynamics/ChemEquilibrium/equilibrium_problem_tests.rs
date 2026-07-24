//! Tests for the typed equilibrium-problem boundary.
//!
//! The tests protect the distinction between physical mole numbers and the
//! log-mole nonlinear iterate, and ensure malformed input is rejected before
//! residual/Jacobian construction can panic.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_component::EquilibriumComponentDescriptor;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    ContinuationSeedPolicy, EquilibriumLogMoles, GibbsFn, Phase, PhaseKind, Solvers,
    compute_species_moles, equilibrium_logmole_jacobian, equilibrium_logmole_residual,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
    ResidualScalingContract, TraceSpeciesSeedPolicy, VariableScalingContract,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverBackend, SolverCascadeBudget, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::{
    EquilibriumAcceptanceCriteria, EquilibriumCandidateResiduals, validate_equilibrium_candidate,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::multiphase_equilibrium_residual_generator_sym;
use crate::Thermodynamics::User_PhaseOrSolution::PhaseModel;
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::PhysicalState;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use std::rc::Rc;

fn assert_jacobian_matches_central_difference(
    prepared: &PreparedEquilibriumProblem,
    log_moles: &[f64],
    step: f64,
    tolerance: f64,
) {
    let analytic = prepared.jacobian(log_moles).unwrap();
    for column in 0..log_moles.len() {
        let mut plus = log_moles.to_vec();
        let mut minus = log_moles.to_vec();
        plus[column] += step;
        minus[column] -= step;

        let residual_plus = prepared.residual(&plus).unwrap();
        let residual_minus = prepared.residual(&minus).unwrap();
        for row in 0..analytic.nrows() {
            let finite_difference = (residual_plus[row] - residual_minus[row]) / (2.0 * step);
            assert!(
                (analytic[(row, column)] - finite_difference).abs() <= tolerance,
                "Jacobian entry ({row}, {column}) differs at log-moles {log_moles:?}: \
                 analytic={}, finite_difference={finite_difference}",
                analytic[(row, column)],
            );
        }
    }
}

fn two_species_problem() -> EquilibriumProblem {
    let initial_moles = vec![1.0, 0.0];
    EquilibriumProblem::new(
        vec!["O2".to_string(), "O".to_string()],
        initial_moles.clone(),
        LogMolesInitialGuess::from_moles(&initial_moles, 1e-20).unwrap(),
        DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
        vec![Rc::new(|_| 0.0) as GibbsFn, Rc::new(|_| 0.0) as GibbsFn],
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }],
        EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
    )
    .unwrap()
}

fn gas_and_pure_condensed_problem() -> EquilibriumProblem {
    let initial_moles = vec![1.0, 1.0, 0.0];
    EquilibriumProblem::new(
        vec!["A_g".to_string(), "B_g".to_string(), "A_s".to_string()],
        initial_moles.clone(),
        LogMolesInitialGuess::from_moles(&initial_moles, 1e-20).unwrap(),
        DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 1.0, 1.0, 0.0]),
        vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| -1_000.0) as GibbsFn,
        ],
        vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![2],
            },
        ],
        EquilibriumConditions::new(1_000.0, 202_650.0, 101_325.0).unwrap(),
    )
    .unwrap()
}

#[test]
fn log_moles_guess_converts_trace_species_without_losing_the_coordinate_system() {
    let guess = LogMolesInitialGuess::from_moles(&[1.0, 0.0], 1e-20).unwrap();

    assert_eq!(guess.as_slice()[0], 0.0);
    assert!((guess.as_slice()[1] - 1e-20_f64.ln()).abs() < 1e-12);
}

#[test]
fn relative_trace_seed_policy_scales_with_the_largest_initial_mole() {
    let guess = LogMolesInitialGuess::from_moles_with_policy(
        &[2.0, 0.0, 0.5],
        TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
            fraction: 1e-6,
            minimum_floor: 1e-30,
        },
    )
    .unwrap();

    assert_eq!(guess.as_slice()[0], 2.0_f64.ln());
    assert!((guess.as_slice()[1] - (2.0e-6_f64).ln()).abs() < 1e-12);
    assert!((guess.as_slice()[2] - 0.5_f64.ln()).abs() < 1e-12);
}

#[test]
fn relative_trace_seed_policy_uses_minimum_floor_for_an_all_zero_system() {
    let guess = LogMolesInitialGuess::from_moles_with_policy(
        &[0.0, 0.0],
        TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
            fraction: 1e-6,
            minimum_floor: 1e-25,
        },
    )
    .unwrap();

    assert!((guess.as_slice()[0] - 1e-25_f64.ln()).abs() < 1e-12);
    assert!((guess.as_slice()[1] - 1e-25_f64.ln()).abs() < 1e-12);
}

#[test]
fn relative_trace_seed_policy_rejects_invalid_parameters_before_seed_construction() {
    let error = LogMolesInitialGuess::from_moles_with_policy(
        &[1.0, 0.0],
        TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
            fraction: 0.0,
            minimum_floor: 1e-25,
        },
    )
    .unwrap_err();

    assert!(matches!(
        error,
        ReactionExtentError::InvalidProblem {
            field: "trace_floor_fraction",
            ..
        }
    ));
}

#[test]
fn canonical_residual_reports_the_specific_bad_gibbs_species() {
    let initial_moles = vec![1.0, 0.0];
    let problem = EquilibriumProblem::new(
        vec!["O2".to_string(), "O".to_string()],
        initial_moles.clone(),
        LogMolesInitialGuess::from_moles(&initial_moles, 1e-20).unwrap(),
        DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
        vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| f64::NAN) as GibbsFn,
        ],
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }],
        EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
    )
    .unwrap();
    let prepared = PreparedEquilibriumProblem::new(problem).unwrap();

    assert!(matches!(
        prepared.residual(&[0.0, -10.0]),
        Err(ReactionExtentError::InvalidDG0 {
            species_index: 1,
            dg0,
            temperature: 1000.0,
        }) if dg0.is_nan()
    ));
}

#[test]
fn default_trace_seed_is_a_named_numerical_coordinate_contract() {
    let guess = LogMolesInitialGuess::from_initial_moles(&[1.0, 0.0]).unwrap();

    assert_eq!(guess.as_slice()[0], 0.0);
    assert_eq!(
        guess.as_slice()[1],
        crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::DEFAULT_TRACE_MOLE_FLOOR.ln()
    );
}

#[test]
fn log_mole_reconstruction_rejects_overflow_and_underflow() {
    assert!(matches!(
        compute_species_moles(&[1000.0]),
        Err(ReactionExtentError::InvalidCandidate { .. })
    ));
    assert!(matches!(
        compute_species_moles(&[-1000.0]),
        Err(ReactionExtentError::InvalidCandidate { .. })
    ));
}

#[test]
fn variable_scaling_contract_applies_to_coordinates_without_touching_residual_rows() {
    let variable_scale = VariableScalingContract::new(vec![2.0, 4.0]).unwrap();
    let residual_scale = ResidualScalingContract::new(vec![3.0, 5.0]).unwrap();
    let iterate = vec![10.0, 20.0];
    let residual = vec![30.0, 50.0];

    assert_eq!(
        variable_scale.apply_iterate(&iterate).unwrap(),
        vec![5.0, 5.0]
    );
    assert_eq!(
        variable_scale.unscale_iterate(&[5.0, 5.0]).unwrap(),
        iterate
    );
    assert_eq!(
        residual_scale.apply_residual(residual.clone()).unwrap(),
        vec![10.0, 10.0]
    );
    assert_eq!(residual, vec![30.0, 50.0]);
    assert_eq!(variable_scale.as_slice(), &[2.0, 4.0]);
}

#[test]
fn variable_scaling_contract_rejects_shape_and_value_errors() {
    let err = VariableScalingContract::new(vec![]).unwrap_err();
    assert!(matches!(
        err,
        ReactionExtentError::InvalidProblem {
            field: "variable_scale",
            ..
        }
    ));

    let err = VariableScalingContract::new(vec![1.0, 0.0]).unwrap_err();
    assert!(matches!(
        err,
        ReactionExtentError::InvalidProblem {
            field: "variable_scale",
            ..
        }
    ));

    let contract = VariableScalingContract::new(vec![1.0, 2.0]).unwrap();
    assert!(matches!(
        contract.apply_iterate(&[1.0]),
        Err(ReactionExtentError::DimensionMismatch(_))
    ));
    assert!(matches!(
        contract.unscale_iterate(&[1.0]),
        Err(ReactionExtentError::DimensionMismatch(_))
    ));
}

#[test]
fn equilibrium_problem_preserves_validated_input_order() {
    let problem = two_species_problem();

    assert_eq!(problem.species(), ["O2".to_string(), "O".to_string()]);
    assert_eq!(problem.initial_moles(), [1.0, 0.0]);
    assert_eq!(problem.phases().len(), 1);
    assert_eq!(problem.conditions().temperature(), 1000.0);
}

#[test]
fn typed_problem_exposes_bounded_identifiers_for_species_phases_elements_and_reactions() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();

    assert_eq!(prepared.species_id(0).unwrap().index(), 0);
    assert_eq!(prepared.species_id(1).unwrap().index(), 1);
    assert!(prepared.species_id(2).is_err());

    assert_eq!(prepared.phase_index(0).unwrap().index(), 0);
    assert!(prepared.phase_index(1).is_err());

    assert_eq!(prepared.element_id(0).unwrap().index(), 0);
    assert!(prepared.element_id(1).is_err());

    assert_eq!(prepared.reaction_id(0).unwrap().index(), 0);
    assert!(prepared.reaction_id(1).is_err());
}

#[test]
fn species_capacity_reports_track_the_limiting_element_for_each_species() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();

    let oxygen_molecule = prepared.species_capacity_report(0).unwrap();
    assert_eq!(oxygen_molecule.species_index, 0);
    assert_eq!(oxygen_molecule.species_name, "O2");
    assert_eq!(oxygen_molecule.limits.len(), 1);
    assert_eq!(oxygen_molecule.maximum_moles(), 1.0);
    assert_eq!(oxygen_molecule.limiting_limit.element_index, 0);
    assert_eq!(oxygen_molecule.limiting_limit.implied_species_capacity, 1.0);

    let oxygen_atom = prepared.species_capacity_report(1).unwrap();
    assert_eq!(oxygen_atom.species_index, 1);
    assert_eq!(oxygen_atom.species_name, "O");
    assert_eq!(oxygen_atom.limits.len(), 1);
    assert_eq!(oxygen_atom.maximum_moles(), 2.0);
    assert_eq!(oxygen_atom.limiting_limit.element_index, 0);
    assert_eq!(oxygen_atom.limiting_limit.implied_species_capacity, 2.0);
}

#[test]
fn equilibrium_problem_preview_exposes_read_only_problem_summary() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let preview = prepared.preview().unwrap();

    assert_eq!(preview.conditions.temperature(), 1000.0);
    assert_eq!(preview.species, vec!["O2".to_string(), "O".to_string()]);
    assert_eq!(preview.initial_moles, vec![1.0, 0.0]);
    assert_eq!(preview.element_totals, vec![2.0]);
    assert_eq!(preview.reaction_basis_rank, 1);
    assert_eq!(preview.reaction_count, 1);
    assert_eq!(preview.species_capacity_reports.len(), 2);
    assert!(preview.formulation_diagnostics.is_none());
}

#[test]
fn equilibrium_problem_preview_with_diagnostics_includes_svd_summary_and_null_directions() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let preview = prepared.preview_with_diagnostics(1e-12).unwrap();
    let diagnostics = preview
        .formulation_diagnostics
        .expect("diagnostics must be present");

    assert_eq!(diagnostics.element_rank, 1);
    assert_eq!(diagnostics.reaction_count, 1);
    assert_eq!(diagnostics.tolerance, 1e-12);
    assert_eq!(diagnostics.singular_values.len(), 1);
    assert!(diagnostics.singular_values[0].is_finite());
    assert!(diagnostics.condition_estimate.is_some());
    assert_eq!(diagnostics.near_null_directions.len(), 1);
    assert_eq!(diagnostics.near_null_directions[0].direction.len(), 2);
}

#[test]
fn equilibrium_problem_preview_summary_rows_are_stable_and_human_readable() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let preview = prepared.preview_with_diagnostics(1e-12).unwrap();
    let rows = preview.summary_rows();

    assert!(
        rows.iter().any(|row| row.section == "problem"
            && row.label == "species_count"
            && row.value == "2")
    );
    assert!(
        rows.iter()
            .any(|row| row.section == "species_capacity" && row.label == "O2")
    );
    assert!(
        rows.iter()
            .any(|row| row.section == "diagnostics" && row.label == "element_rank")
    );

    let rendered = format!("{preview}");
    assert!(rendered.contains("[problem] temperature = 1000.000000"));
    assert!(rendered.contains("[species_capacity] O2 = O2: max 1.000000 mol limited by element 0"));
}

#[test]
fn typed_reports_have_stable_display_output() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let capacity = prepared.species_capacity_report(0).unwrap();
    let diagnostics = prepared.formulation_diagnostics(1e-12).unwrap();

    assert!(format!("{capacity}").contains("O2: max 1.000000 mol limited by element 0"));
    assert!(format!("{diagnostics}").contains("element rank = 1, reaction count = 1"));
    assert!(format!("{diagnostics}").contains("sv[0]"));
}

#[test]
fn formulation_diagnostics_rejects_non_positive_tolerance() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let error = prepared.formulation_diagnostics(0.0).unwrap_err();

    assert!(matches!(
        error,
        ReactionExtentError::InvalidProblem {
            field: "diagnostic_tolerance",
            ..
        }
    ));
}

#[test]
fn equilibrium_problem_rejects_species_rows_without_any_positive_element_support() {
    let initial_moles = vec![1.0];
    let result = EquilibriumProblem::new(
        vec!["X".to_string()],
        initial_moles.clone(),
        LogMolesInitialGuess::from_moles(&initial_moles, 1e-20).unwrap(),
        DMatrix::from_row_slice(1, 1, &[0.0]),
        vec![Rc::new(|_| 0.0) as GibbsFn],
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0],
        }],
        EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
    );

    assert!(matches!(
        result,
        Err(ReactionExtentError::InvalidProblem {
            field: "element_composition",
            ..
        })
    ));
}

#[test]
fn formulation_blocks_remain_aligned_and_dimension_safe_for_the_canonical_two_species_case() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let log_moles = [0.0, -10.0];

    let residual = prepared.residual(&log_moles).unwrap();
    let jacobian = prepared.jacobian(&log_moles).unwrap();

    assert_eq!(residual.len(), 2);
    assert_eq!(jacobian.nrows(), 2);
    assert_eq!(jacobian.ncols(), 2);
    assert!(residual.iter().all(|value| value.is_finite()));
    assert!(jacobian.iter().all(|value| value.is_finite()));

    let moles = compute_species_moles(&log_moles).unwrap();
    assert_eq!(moles.len(), 2);
    assert!(moles.iter().all(|value| value.is_finite() && *value > 0.0));

    let candidate = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &log_moles,
            raw_residual: &residual,
            acceptance_residual: &residual,
        },
        EquilibriumAcceptanceCriteria::new(1e6, 1e6, 1e6).unwrap(),
        prepared.problem().element_composition(),
        prepared.element_totals(),
    )
    .unwrap();

    assert!(candidate.residual_l2_norm.is_finite());
    assert!(candidate.max_abs_element_balance_error.is_finite());
    assert!(candidate.min_moles > 0.0);
}

#[test]
fn typed_problem_builds_a_legacy_solver_without_reinterpreting_moles_as_logs() {
    let solver = EquilibriumLogMoles::from_problem(two_species_problem()).unwrap();

    assert_eq!(solver.n0, vec![1.0, 0.0]);
    assert_eq!(solver.initial_guess, Some(vec![0.0, 1e-20_f64.ln()]));
    assert_eq!(solver.elements_vector, vec![2.0]);
    assert_eq!(solver.species_phase, vec![0, 0]);
}

#[test]
fn invalid_explicit_cascade_budget_stops_before_any_backend_or_publication() {
    let mut solver = EquilibriumLogMoles::from_problem(two_species_problem()).unwrap();
    solver.set_solver_policy(SolverPolicy::Cascade(vec![
        SolverBackend::Legacy(Solvers::LM),
        SolverBackend::Legacy(Solvers::NR),
    ]));
    solver.set_solver_budget(SolverCascadeBudget::new(0, 1, 1));

    assert!(matches!(
        solver.solve(),
        Err(ReactionExtentError::InvalidProblem {
            field: "solver_budget",
            ..
        })
    ));
    assert!(solver.solution.is_empty());
    assert!(solver.moles.is_empty());
    assert!(solver.last_validation_report.is_none());
    assert!(solver.last_solve_report.is_none());
}

#[test]
fn invalid_mutated_input_stops_before_any_backend_attempt_or_state_publication() {
    let mut solver = EquilibriumLogMoles::from_problem(two_species_problem()).unwrap();
    solver.initial_guess = None;
    solver.n0.pop();
    solver.set_solver_policy(SolverPolicy::Cascade(vec![
        SolverBackend::Legacy(Solvers::LM),
        SolverBackend::Legacy(Solvers::NR),
    ]));

    let error = solver.solve().unwrap_err();
    assert!(matches!(error, ReactionExtentError::DimensionMismatch(_)));
    assert!(solver.solution.is_empty());
    assert!(solver.moles.is_empty());
    assert!(solver.last_validation_report.is_none());
    assert!(solver.last_solve_report.is_none());
    assert!(solver.initial_guess.is_none());
}

#[test]
fn failed_solve_after_a_successful_publication_keeps_the_last_accepted_state_intact() {
    let mut solver = EquilibriumLogMoles::from_problem(two_species_problem()).unwrap();
    solver.set_solver_policy(SolverPolicy::Cascade(vec![
        SolverBackend::Legacy(Solvers::LM),
        SolverBackend::Legacy(Solvers::NR),
    ]));

    solver.solve().unwrap();
    let published_solution = solver.solution.clone();
    let published_moles = solver.moles.clone();
    let published_validation = solver.last_validation_report.clone();
    let published_report = solver.last_solve_report.clone();
    let published_initial_guess = solver.initial_guess.clone();

    solver.n0.pop();
    let error = solver.solve().unwrap_err();

    assert!(matches!(error, ReactionExtentError::DimensionMismatch(_)));
    assert_eq!(solver.solution, published_solution);
    assert_eq!(solver.moles, published_moles);
    assert_eq!(solver.last_validation_report, published_validation);
    assert_eq!(solver.last_solve_report, published_report);
    assert_eq!(solver.initial_guess, published_initial_guess);
}

#[test]
fn mutable_initial_guess_setter_rejects_non_finite_and_wrong_sized_seeds() {
    let mut solver = EquilibriumLogMoles::from_problem(two_species_problem()).unwrap();
    let original = solver.initial_guess.clone();

    assert!(matches!(
        solver.set_initial_guess(vec![0.0, f64::NAN]),
        Err(ReactionExtentError::InvalidProblem {
            field: "initial_log_moles",
            ..
        })
    ));
    assert_eq!(solver.initial_guess, original);

    assert!(matches!(
        solver.set_initial_guess(vec![0.0]),
        Err(ReactionExtentError::DimensionMismatch(_))
    ));
    assert_eq!(solver.initial_guess, original);
}

#[test]
fn prepared_problem_owns_deterministic_formulation_data_and_evaluates_it_purely() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let log_moles = prepared.problem().initial_log_moles().as_slice();

    assert_eq!(prepared.element_totals(), [2.0]);
    assert_eq!(prepared.species_phase(), [0, 0]);
    assert_eq!(prepared.reaction_basis().reactions.nrows(), 2);
    assert_eq!(prepared.reaction_basis().reactions.ncols(), 1);
    assert_eq!(prepared.phase_stoichiometry().len(), 1);

    let residual = prepared.residual(log_moles).unwrap();
    let jacobian = prepared.jacobian(log_moles).unwrap();
    assert_eq!(residual.len(), 2);
    assert!(residual.iter().all(|value| value.is_finite()));
    assert_eq!(jacobian.shape(), (2, 2));
    assert!(jacobian.iter().all(|value| value.is_finite()));
}

#[test]
fn prepared_problem_reconstructs_and_packages_accepted_results_in_one_ordering() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let log_moles = vec![0.0, -50.0];
    let moles = prepared.reconstruct_moles(&log_moles).unwrap();
    let residual = prepared.residual(&log_moles).unwrap();
    let validation = validate_equilibrium_candidate(
        EquilibriumCandidateResiduals {
            log_moles: &log_moles,
            raw_residual: &residual,
            acceptance_residual: &residual,
        },
        EquilibriumAcceptanceCriteria::new(100.0, 1e-12, 100.0).unwrap(),
        prepared.problem().element_composition(),
        prepared.element_totals(),
    )
    .unwrap();
    let solution = prepared
        .accepted_solution(log_moles.clone(), validation)
        .unwrap();

    assert_eq!(solution.log_moles(), log_moles);
    assert_eq!(solution.moles(), moles);
    assert_eq!(solution.conditions(), prepared.problem().conditions());
    assert!(matches!(
        prepared.reconstruct_moles(&[0.0]),
        Err(ReactionExtentError::DimensionMismatch(_))
    ));
}

#[test]
fn prepared_problem_scales_residual_and_jacobian_with_one_explicit_contract() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let log_moles = [0.0, -4.0];
    let scaling = prepared.residual_scaling_contract().unwrap();
    let scale = scaling.as_slice().to_vec();
    let raw_residual = prepared.residual(&log_moles).unwrap();
    let raw_jacobian = prepared.jacobian(&log_moles).unwrap();
    let scaled_residual = scaling.apply_residual(raw_residual.clone()).unwrap();
    let scaled_jacobian = scaling.apply_jacobian(raw_jacobian.clone()).unwrap();

    for row in 0..scale.len() {
        assert_eq!(scaled_residual[row], raw_residual[row] / scale[row]);
        for column in 0..raw_jacobian.ncols() {
            assert_eq!(
                scaled_jacobian[(row, column)],
                raw_jacobian[(row, column)] / scale[row]
            );
        }
    }
    assert!(matches!(
        prepared.scaled_residual(&log_moles, &[1.0]),
        Err(ReactionExtentError::DimensionMismatch(_))
    ));
    assert!(matches!(
        ResidualScalingContract::new(vec![1.0, 0.0]),
        Err(ReactionExtentError::InvalidProblem {
            field: "residual_scale",
            ..
        })
    ));
}

#[test]
fn analytic_jacobian_matches_finite_differences_for_normal_and_trace_species() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();

    assert_jacobian_matches_central_difference(&prepared, &[0.0, -2.0], 1e-6, 1e-6);
    assert_jacobian_matches_central_difference(&prepared, &[0.0, -25.0], 1e-6, 1e-6);
}

#[test]
fn legacy_residual_closure_delegates_to_the_pure_prepared_formulation() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let problem = prepared.problem();
    let conditions = problem.conditions();
    let closure = equilibrium_logmole_residual(
        prepared.reaction_basis().reactions.clone(),
        problem.element_composition().clone(),
        prepared.element_totals().to_vec(),
        problem.gibbs().to_vec(),
        problem.phases().to_vec(),
        conditions.temperature(),
        conditions.pressure(),
        conditions.reference_pressure(),
        prepared.species_phase().to_vec(),
        0.0,
        0.0,
    )
    .unwrap();
    let log_moles = [0.0, -4.0];

    assert_eq!(
        closure(&log_moles).unwrap(),
        prepared.residual(&log_moles).unwrap()
    );
}

#[test]
fn numeric_legacy_closure_and_symbolic_residuals_share_one_formulation_contract() {
    // Constant standard Gibbs energies keep this fixture offline and isolate
    // formulation equivalence from database lookup or thermochemical fitting.
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let problem = prepared.problem();
    let conditions = problem.conditions();
    let closure = equilibrium_logmole_residual(
        prepared.reaction_basis().reactions.clone(),
        problem.element_composition().clone(),
        prepared.element_totals().to_vec(),
        problem.gibbs().to_vec(),
        problem.phases().to_vec(),
        conditions.temperature(),
        conditions.pressure(),
        conditions.reference_pressure(),
        prepared.species_phase().to_vec(),
        0.0,
        0.0,
    )
    .unwrap();
    let log_moles = vec![0.0, -4.0];
    let canonical = prepared.residual(&log_moles).unwrap();
    let legacy = closure(&log_moles).unwrap();
    assert_eq!(legacy, canonical);

    let symbolic = multiphase_equilibrium_residual_generator_sym(
        prepared.reaction_basis().reactions.clone(),
        problem.element_composition().clone(),
        prepared.element_totals().to_vec(),
        vec![Expr::Const(0.0); problem.species().len()],
        problem.phases().to_vec(),
        conditions.pressure(),
        conditions.reference_pressure(),
    )
    .unwrap();
    let variables = Expr::IndexedVars(problem.species().len(), "y").0;
    let variable_names: Vec<String> = variables.iter().map(ToString::to_string).collect();
    let variable_refs: Vec<&str> = variable_names.iter().map(String::as_str).collect();

    assert_eq!(symbolic.len(), canonical.len());
    for (index, expression) in symbolic.iter().enumerate() {
        let expression = expression
            .set_variable("T", conditions.temperature())
            .simplify();
        let evaluate = expression.lambdify_borrowed_thread_safe(&variable_refs);
        let symbolic_value = evaluate(&log_moles);
        assert!(
            (symbolic_value - canonical[index]).abs() <= 1e-12,
            "residual {index} diverged: symbolic={symbolic_value}, canonical={}",
            canonical[index],
        );
    }
}

#[test]
fn analytic_jacobian_matches_symbolic_residual_derivatives() {
    // This exercises the same residual expressions through three independent
    // representations: the prepared numeric formulation, symbolic
    // differentiation, and lambdified symbolic evaluation.
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let problem = prepared.problem();
    let conditions = problem.conditions();
    let log_moles = vec![0.0, -4.0];
    let canonical = prepared.jacobian(&log_moles).unwrap();

    let symbolic = multiphase_equilibrium_residual_generator_sym(
        prepared.reaction_basis().reactions.clone(),
        problem.element_composition().clone(),
        prepared.element_totals().to_vec(),
        vec![Expr::Const(0.0); problem.species().len()],
        problem.phases().to_vec(),
        conditions.pressure(),
        conditions.reference_pressure(),
    )
    .unwrap();
    let variables = Expr::IndexedVars(problem.species().len(), "y").0;
    let variable_names: Vec<String> = variables.iter().map(ToString::to_string).collect();
    let variable_refs: Vec<&str> = variable_names.iter().map(String::as_str).collect();

    for (row, residual) in symbolic.iter().enumerate() {
        for (column, variable) in variable_refs.iter().enumerate() {
            let derivative = residual
                .diff(variable)
                .set_variable("T", conditions.temperature())
                .simplify();
            let evaluate = derivative.lambdify_borrowed_thread_safe(&variable_refs);
            let symbolic_value = evaluate(&log_moles);
            assert!(
                (symbolic_value - canonical[(row, column)]).abs() <= 1e-12,
                "Jacobian entry ({row}, {column}) diverged: symbolic={symbolic_value}, canonical={}",
                canonical[(row, column)],
            );
        }
    }
}

#[test]
fn gas_and_pure_condensed_numeric_symbolic_and_jacobian_contracts_agree() {
    let prepared = PreparedEquilibriumProblem::new(gas_and_pure_condensed_problem()).unwrap();
    let problem = prepared.problem();
    let conditions = problem.conditions();
    let log_moles = vec![0.7_f64.ln(), 1.0_f64.ln(), 0.3_f64.ln()];
    let numeric = prepared.residual(&log_moles).unwrap();
    let jacobian = prepared.jacobian(&log_moles).unwrap();

    let symbolic = multiphase_equilibrium_residual_generator_sym(
        prepared.reaction_basis().reactions.clone(),
        problem.element_composition().clone(),
        prepared.element_totals().to_vec(),
        vec![Expr::Const(0.0), Expr::Const(0.0), Expr::Const(-1_000.0)],
        problem.phases().to_vec(),
        conditions.pressure(),
        conditions.reference_pressure(),
    )
    .unwrap();
    let variables = Expr::IndexedVars(problem.species().len(), "y").0;
    let names = variables
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let refs = names.iter().map(String::as_str).collect::<Vec<_>>();

    for (row, residual) in symbolic.iter().enumerate() {
        let residual = residual
            .clone()
            .set_variable("T", conditions.temperature())
            .simplify();
        let evaluate = residual.lambdify_borrowed_thread_safe(&refs);
        assert!((evaluate(&log_moles) - numeric[row]).abs() <= 1e-12);

        for (column, variable) in refs.iter().enumerate() {
            let derivative = residual.diff(variable).simplify();
            let evaluate = derivative.lambdify_borrowed_thread_safe(&refs);
            assert!(
                (evaluate(&log_moles) - jacobian[(row, column)]).abs() <= 1e-12,
                "mixed-phase Jacobian entry ({row}, {column}) diverged: symbolic={}, analytic={}",
                evaluate(&log_moles),
                jacobian[(row, column)]
            );
        }
    }
}

#[test]
fn legacy_jacobian_wrapper_delegates_to_the_pure_prepared_formulation() {
    let prepared = PreparedEquilibriumProblem::new(two_species_problem()).unwrap();
    let problem = prepared.problem();
    let log_moles = [0.0, -4.0];

    let legacy = equilibrium_logmole_jacobian(
        &log_moles,
        &prepared.reaction_basis().reactions,
        problem.element_composition(),
        prepared.species_phase(),
        prepared.phase_stoichiometry(),
        problem.phases().len(),
        0.0,
        0.0,
    )
    .unwrap();

    assert_eq!(legacy, prepared.jacobian(&log_moles).unwrap());
}

#[test]
fn legacy_solver_default_seed_uses_log_moles_when_the_typed_seed_is_not_set() {
    let mut solver = EquilibriumLogMoles::empty();
    solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
    solver.n0 = vec![1.0, 0.0];
    solver.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
    solver.gibbs = vec![Rc::new(|_| 0.0) as GibbsFn, Rc::new(|_| 0.0) as GibbsFn];
    solver.phases = vec![Phase {
        kind: PhaseKind::IdealGas,
        species: vec![0, 1],
    }];
    solver.create_stoich_matrix().unwrap();

    let _ = solver.solve();

    assert_eq!(solver.initial_guess, Some(vec![0.0, 1e-30_f64.ln()]));
}

#[test]
fn solve_problem_with_rejects_invalid_solver_settings_before_publication() {
    let error = EquilibriumLogMoles::solve_problem_with(two_species_problem(), |settings| {
        settings.solver_policy = Some(SolverPolicy::Cascade(Vec::new()));
        settings.continuation_seed_policy = ContinuationSeedPolicy::IndependentPerPoint;
    })
    .unwrap_err();

    assert!(matches!(
        error,
        ReactionExtentError::InvalidProblem {
            field: "solver_policy",
            ..
        }
    ));
}

#[test]
fn canonical_solver_settings_reject_invalid_numerics_before_formulation() {
    let error = EquilibriumLogMoles::solve_problem_with(two_species_problem(), |settings| {
        settings.solver_params.tol = 0.0;
    })
    .unwrap_err();

    assert!(matches!(
        error,
        ReactionExtentError::InvalidProblem {
            field: "solver_params.tol",
            ..
        }
    ));
}

#[test]
fn solve_problem_returns_the_immutable_accepted_snapshot() {
    let solution = EquilibriumLogMoles::solve_problem(two_species_problem()).unwrap();

    assert_eq!(solution.conditions().temperature(), 1000.0);
    assert_eq!(solution.conditions().pressure(), 101325.0);
    assert_eq!(solution.conditions().reference_pressure(), 101325.0);
    assert_eq!(solution.log_moles().len(), 2);
    assert_eq!(solution.moles().len(), 2);
    assert!(solution.moles().iter().all(|value| *value > 0.0));
    assert!(solution.validation().residual_l2_norm.is_finite());
    assert!(
        solution
            .validation()
            .max_abs_element_balance_error
            .is_finite()
    );
}

#[test]
fn solve_publishes_optional_keq_cross_validation_status_when_enabled() {
    let mut solver = EquilibriumLogMoles::from_problem(two_species_problem()).unwrap();
    solver.solver_settings.keq_validation_mode = EquilibriumConstantValidationMode::WhenApplicable;
    solver.solve().unwrap();

    let status = solver
        .last_keq_validation_status
        .as_ref()
        .expect("enabled K_eq validation should publish a status");
    match status {
        EquilibriumConstantCrossValidationStatus::Compared(report) => {
            assert!(report.accepted);
            assert!(report.max_abs_species_mole_delta.is_finite());
        }
        other => panic!("unexpected K_eq validation status: {other:?}"),
    }
}

#[test]
fn solve_keeps_keq_cross_validation_disabled_by_default() {
    let mut solver = EquilibriumLogMoles::from_problem(two_species_problem()).unwrap();
    solver.solve().unwrap();

    assert!(solver.last_keq_validation_status.is_none());
}

#[test]
fn equilibrium_problem_rejects_duplicate_species_before_solver_setup() {
    let result = EquilibriumProblem::new(
        vec!["O2".to_string(), "O2".to_string()],
        vec![1.0, 0.0],
        LogMolesInitialGuess::new(vec![0.0, -20.0]).unwrap(),
        DMatrix::from_row_slice(2, 1, &[2.0, 2.0]),
        vec![Rc::new(|_| 0.0) as GibbsFn, Rc::new(|_| 0.0) as GibbsFn],
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }],
        EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
    );

    assert!(matches!(
        result,
        Err(ReactionExtentError::InvalidProblem {
            field: "components",
            ..
        })
    ));
}

#[test]
fn phase_qualified_components_allow_same_substance_in_distinct_phases() {
    let gas_water = EquilibriumComponentDescriptor::new(
        PhaseComponentId::new(PhaseId::new(Some("gas".to_string())), "H2O"),
        PhysicalState::Gas,
        PhaseModel::IdealGas,
        PhaseActivityModel::IdealGas,
    );
    let liquid_water = EquilibriumComponentDescriptor::new(
        PhaseComponentId::new(PhaseId::new(Some("liquid".to_string())), "H2O"),
        PhysicalState::Liquid,
        PhaseModel::PureCondensed,
        PhaseActivityModel::IdealSolution,
    );

    let problem = EquilibriumProblem::new(
        vec![gas_water.clone(), liquid_water],
        vec![1.0, 0.0],
        LogMolesInitialGuess::new(vec![0.0, -20.0]).unwrap(),
        DMatrix::from_row_slice(2, 2, &[2.0, 1.0, 2.0, 1.0]),
        vec![Rc::new(|_| 0.0) as GibbsFn, Rc::new(|_| 0.0) as GibbsFn],
        vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
        ],
        EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
    )
    .unwrap();

    assert_eq!(problem.species(), ["gas::H2O", "liquid::H2O"]);
    assert_eq!(problem.components()[0], gas_water);

    let duplicate = EquilibriumProblem::new(
        vec![gas_water.clone(), gas_water],
        vec![1.0, 0.0],
        LogMolesInitialGuess::new(vec![0.0, -20.0]).unwrap(),
        DMatrix::from_row_slice(2, 2, &[2.0, 1.0, 2.0, 1.0]),
        vec![Rc::new(|_| 0.0) as GibbsFn, Rc::new(|_| 0.0) as GibbsFn],
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }],
        EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
    );
    assert!(matches!(
        duplicate,
        Err(ReactionExtentError::InvalidProblem {
            field: "components",
            ..
        })
    ));
}

#[test]
fn equilibrium_problem_rejects_invalid_phase_index_before_phase_mapping() {
    let result = EquilibriumProblem::new(
        vec!["O2".to_string(), "O".to_string()],
        vec![1.0, 0.0],
        LogMolesInitialGuess::new(vec![0.0, -20.0]).unwrap(),
        DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
        vec![Rc::new(|_| 0.0) as GibbsFn, Rc::new(|_| 0.0) as GibbsFn],
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 2],
        }],
        EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
    );

    assert!(matches!(
        result,
        Err(ReactionExtentError::InvalidProblem {
            field: "phases",
            ..
        })
    ));
}
