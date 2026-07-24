//! First physical backend matrix for the canonical equilibrium solver.
//!
//! Every selected RST strategy solves the same small dissociation system from
//! the same initial state. The assertions cover the acceptance gate rather
//! than a solver-specific iterate, so this remains a regression test for the
//! physical formulation and the backend boundary together.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{Phase, PhaseKind, Solvers};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverAttemptOutcome, SolverBackend, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::gas_solver;

fn o2_dissociation_solver(
    backend: RustedSciTheSolver,
) -> crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumLogMoles {
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
    solver.set_solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(backend)));
    solver
}

fn n2_dissociation_solver(
    backend: RustedSciTheSolver,
) -> crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumLogMoles {
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
    solver.set_solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(backend)));
    solver
}

fn diluted_o2_solver(
    backend: RustedSciTheSolver,
) -> crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumLogMoles {
    let mut solver = gas_solver(
        vec!["O2".to_string(), "O".to_string(), "N2".to_string()],
        1000.0,
        202_650.0,
        Solvers::LM,
        None,
        false,
    )
    .unwrap();
    solver.n0 = vec![1.0, 1e-5, 5.0];
    solver.initial_guess = Some(vec![1.0, 1e-5, 5.0]);
    solver.phases = vec![Phase {
        kind: PhaseKind::IdealGas,
        species: vec![0, 1, 2],
    }];
    solver.solver_settings.solver_params.tol = 1e-5;
    solver.create_stoich_matrix().unwrap();
    solver.set_solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(backend)));
    solver
}

#[test]
fn every_selected_rst_backend_accepts_the_o2_dissociation_fixture() {
    for backend in RustedSciTheSolver::recommended_cascade() {
        let mut solver = o2_dissociation_solver(backend);
        solver.solve().unwrap_or_else(|error| {
            panic!(
                "{} should solve and pass O2/O acceptance: {error:?}",
                backend.name()
            )
        });

        let validation = solver.last_validation_report.as_ref().unwrap();
        assert!(
            validation.residual_l2_norm.is_finite(),
            "{}",
            backend.name()
        );
        assert!(
            validation.max_abs_element_balance_error <= 1e-5,
            "{}",
            backend.name()
        );
        assert!(validation.min_moles > 0.0, "{}", backend.name());

        let attempt = &solver.last_solve_report.as_ref().unwrap().attempts[0];
        assert_eq!(attempt.outcome, SolverAttemptOutcome::Accepted);
        let metrics = attempt.metrics.as_ref().unwrap();
        assert!(metrics.backend_converged, "{}", backend.name());
        assert!(metrics.residual_evaluations > 0, "{}", backend.name());
        assert!(metrics.jacobian_evaluations > 0, "{}", backend.name());
        assert!(metrics.elapsed_millis < u128::MAX, "{}", backend.name());
    }
}

#[test]
fn every_selected_rst_backend_accepts_the_diluted_fixture_and_reports_metrics() {
    for backend in RustedSciTheSolver::recommended_cascade() {
        let mut solver = diluted_o2_solver(backend);
        solver.solve().unwrap_or_else(|error| {
            panic!(
                "{} should solve and pass the diluted O2/O/N2 acceptance: {error:?}",
                backend.name()
            )
        });

        let validation = solver.last_validation_report.as_ref().unwrap();
        assert!(
            validation.residual_l2_norm.is_finite(),
            "{}",
            backend.name()
        );
        assert!(
            validation.max_abs_element_balance_error <= 1e-5,
            "{}",
            backend.name()
        );
        assert!(validation.min_moles > 0.0, "{}", backend.name());

        let report = solver.last_solve_report.as_ref().unwrap();
        assert_eq!(report.attempts.len(), 1);
        assert_eq!(
            report.accepted_backend,
            SolverBackend::RustedSciThe(backend)
        );
        let attempt = &report.attempts[0];
        assert_eq!(attempt.outcome, SolverAttemptOutcome::Accepted);
        let metrics = attempt.metrics.as_ref().unwrap();
        assert!(metrics.backend_converged, "{}", backend.name());
        assert!(metrics.iterations > 0, "{}", backend.name());
        assert!(metrics.residual_evaluations > 0, "{}", backend.name());
        assert!(metrics.jacobian_evaluations > 0, "{}", backend.name());
        assert!(metrics.elapsed_millis < u128::MAX, "{}", backend.name());

        assert!(solver.moles[2] > 0.0, "{}", backend.name());
    }
}

#[test]
fn high_temperature_n2_dissociation_falls_back_after_rejected_nielsen_candidate() {
    let mut solver = n2_dissociation_solver(RustedSciTheSolver::NielsenLevenbergMarquardt);
    let initial_guess = solver.initial_guess.clone();
    let initial_moles = solver.n0.clone();
    solver.set_solver_policy(SolverPolicy::Cascade(vec![
        SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
        SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
    ]));
    solver.solve().unwrap();

    let validation = solver.last_validation_report.as_ref().unwrap();
    assert!(validation.residual_l2_norm.is_finite());
    assert!(validation.max_abs_element_balance_error <= 1e-5);
    assert!(validation.min_moles > 0.0);
    assert!(solver.moles[0] < 0.5);
    assert_eq!(solver.initial_guess, initial_guess);
    assert_eq!(solver.n0, initial_moles);

    let report = solver.last_solve_report.as_ref().unwrap();
    assert_eq!(report.attempts.len(), 2);
    assert!(matches!(
        report.attempts[0].outcome,
        SolverAttemptOutcome::RejectedCandidate { .. }
    ));
    assert_eq!(
        report.accepted_backend,
        SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt)
    );
    assert_eq!(report.attempts[1].outcome, SolverAttemptOutcome::Accepted);
    assert!(
        report.attempts[1]
            .metrics
            .as_ref()
            .unwrap()
            .backend_converged
    );
}

#[test]
fn rusted_scithe_default_policy_has_a_stable_documented_backend_order() {
    let policy = SolverPolicy::rusted_scithe_default();

    assert_eq!(
        policy.ordered_backends(),
        vec![
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
            SolverBackend::RustedSciThe(RustedSciTheSolver::MinpackLevenbergMarquardt),
            SolverBackend::RustedSciThe(RustedSciTheSolver::NielsenLevenbergMarquardt),
            SolverBackend::RustedSciThe(RustedSciTheSolver::TrustRegionLevenbergMarquardt),
            SolverBackend::RustedSciThe(RustedSciTheSolver::PowellDogleg),
            SolverBackend::RustedSciThe(RustedSciTheSolver::DampedNewton),
        ]
    );
    assert!(
        policy
            .ordered_backends()
            .iter()
            .all(|backend| matches!(backend, SolverBackend::RustedSciThe(_)))
    );
}
