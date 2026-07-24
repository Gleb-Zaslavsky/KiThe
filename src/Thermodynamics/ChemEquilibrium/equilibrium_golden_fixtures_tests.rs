//! Golden regression fixtures for the canonical equilibrium formulation.
//!
//! These tests freeze a couple of representative, currently solvable systems
//! so that future solver migrations do not accidentally change the accepted
//! snapshot shape, fallback trace, or report contract without immediately
//! tripping a regression.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumLogMoles, Phase, PhaseKind, Solvers,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverAttemptOutcome, SolverBackend, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::gas_solver;

fn golden_o2_solver() -> EquilibriumLogMoles {
    let mut solver = gas_solver(
        vec!["O2".to_string(), "O".to_string()],
        500.0,
        101325.0,
        Solvers::LM,
        None,
        false,
    )
    .unwrap();
    solver.n0 = vec![1.0, 1e-5];
    solver.initial_guess = Some(vec![1.0, 1e-5]);
    solver.phases = vec![Phase {
        kind: PhaseKind::IdealGas,
        species: vec![0, 1],
    }];
    solver.create_stoich_matrix().unwrap();
    solver.set_solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
        RustedSciTheSolver::LevenbergMarquardt,
    )));
    solver
}

fn golden_n2_solver() -> EquilibriumLogMoles {
    let mut solver = gas_solver(
        vec!["N2".to_string(), "N".to_string()],
        5500.0,
        101325.0,
        Solvers::LM,
        None,
        false,
    )
    .unwrap();
    solver.n0 = vec![1e-5, 1.0];
    solver.initial_guess = Some(vec![1.0, 1e-5]);
    solver.phases = vec![Phase {
        kind: PhaseKind::IdealGas,
        species: vec![0, 1],
    }];
    solver.create_stoich_matrix().unwrap();
    solver.set_solver_policy(SolverPolicy::Cascade(vec![
        SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
        SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
    ]));
    solver
}

#[test]
fn golden_o2_dissociation_fixture_keeps_a_single_accepted_snapshot() {
    let mut solver = golden_o2_solver();
    solver.solve().unwrap();

    let accepted = solver.accepted_solution().unwrap();
    let report = solver.last_solve_report.as_ref().unwrap();
    let validation = solver.last_validation_report.as_ref().unwrap();

    assert_eq!(report.attempt_count(), 1);
    assert_eq!(report.started_attempt_count(), 1);
    assert_eq!(report.fallback_attempt_count(), 0);
    assert_eq!(report.skipped_attempt_count(), 0);
    assert!(!report.accepted_after_fallback());
    assert_eq!(report.accepted_attempt_index(), Some(0));
    assert_eq!(report.attempt_backends(), vec![report.accepted_backend]);
    assert!(report.summary().contains("attempts=1"));
    assert!(matches!(
        report.accepted_attempt().unwrap().outcome,
        SolverAttemptOutcome::Accepted
    ));

    assert_eq!(accepted.validation(), validation);
    assert_eq!(accepted.moles().len(), 2);
    assert!(accepted.moles().iter().all(|&n| n > 0.0));
    assert!(accepted.moles()[0] > 0.0);
    assert!(accepted.moles()[1] > 0.0);
}

#[test]
fn golden_n2_dissociation_fixture_records_the_fallback_trace() {
    let mut solver = golden_n2_solver();
    let configured_seed = solver.initial_guess.clone();
    let initial_moles = solver.n0.clone();
    solver.solve().unwrap();

    let accepted = solver.accepted_solution().unwrap();
    let report = solver.last_solve_report.as_ref().unwrap();

    assert_eq!(solver.initial_guess, configured_seed);
    assert_eq!(solver.n0, initial_moles);
    assert_eq!(report.attempt_count(), 2);
    assert_eq!(report.started_attempt_count(), 2);
    assert_eq!(report.fallback_attempt_count(), 1);
    assert_eq!(report.skipped_attempt_count(), 0);
    assert!(report.accepted_after_fallback());
    assert_eq!(
        report.attempt_backends(),
        vec![
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
        ]
    );
    assert_eq!(
        report.accepted_backend,
        SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt)
    );
    assert!(report.summary().contains("fallback=1"));
    assert!(matches!(
        report.attempt(0).unwrap().outcome,
        SolverAttemptOutcome::RejectedCandidate { .. }
    ));
    assert!(matches!(
        report.attempt(1).unwrap().outcome,
        SolverAttemptOutcome::Accepted
    ));

    assert_eq!(accepted.log_moles().len(), 2);
    assert_eq!(accepted.moles().len(), 2);
    assert!(accepted.moles()[0] < 0.5);
    assert!(accepted.moles().iter().all(|&n| n > 0.0));
}
