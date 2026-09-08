//! Contract tests for the RustedSciThe equilibrium adapter.
//!
//! These tests deliberately use a tiny symbolic system rather than database
//! data. They verify that every selected RST strategy receives only symbolic
//! residual expressions and lets RST prepare the Jacobian itself.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, RustedSciTheSolveContract, RustedSciTheSolver,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverTermination;
use RustedSciThe::numerical::Nonlinear_systems::prelude::{
    SymbolicNonlinearProblem, SymbolicProblemOptions,
};

fn scalar_symbolic_problem() -> SymbolicNonlinearProblem {
    SymbolicNonlinearProblem::from_strings_with_options(
        vec!["x * x - 9.0".to_string()],
        SymbolicProblemOptions::new().with_variables(vec!["x".to_string()]),
    )
    .unwrap()
}

#[test]
fn every_recommended_rst_strategy_solves_the_same_symbolic_problem() {
    let problem = scalar_symbolic_problem();
    let problem = RstPreparedProblem::new(problem);
    let options = RustedSciTheSolveContract::new(1e-8, 250).unwrap();

    for solver in RustedSciTheSolver::recommended_cascade() {
        let outcome = solver
            .solve(&problem, &[1.0], &options)
            .unwrap_or_else(|error| {
                panic!(
                    "{} should solve the symbolic scalar contract: {error:?}",
                    solver.name()
                )
            });
        assert_eq!(outcome.solution.len(), 1);
        assert!(
            (outcome.solution[0] - 3.0).abs() < 1e-5,
            "{}",
            solver.name()
        );
        assert!(
            outcome.metrics.residual_evaluations > 0,
            "{}",
            solver.name()
        );
        assert!(
            outcome.metrics.jacobian_evaluations > 0,
            "{}",
            solver.name()
        );
        assert!(outcome.metrics.backend_converged, "{}", solver.name());
        assert_eq!(outcome.metrics.termination, SolverTermination::Converged);
        assert!(
            outcome.metrics.evaluation_timing.is_some(),
            "{} must expose callback timing for release diagnosis",
            solver.name()
        );
    }
}

#[test]
fn invalid_rst_contract_rejects_bad_budget_before_solve() {
    let error = RustedSciTheSolveContract::new(1e-8, 0).unwrap_err();

    assert!(matches!(
        error,
        ReactionExtentError::InvalidProblem {
            field: "rst_solve_contract.max_iterations",
            ..
        }
    ));
}

#[test]
fn solve_contract_preserves_finite_log_mole_bounds_in_rst_options() {
    let bounds = [(-700.0, 0.0), (-700.0, 2.0_f64.ln())];
    let options = RustedSciTheSolveContract::new(1e-8, 250)
        .unwrap()
        .with_log_mole_bounds(&bounds)
        .unwrap()
        .to_options()
        .unwrap();

    assert_eq!(
        options
            .bounds
            .as_ref()
            .expect("fixed-P,T contract must pass box bounds to RST")
            .as_slice(),
        bounds
    );
}
