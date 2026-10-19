//! Contract tests for the RustedSciThe equilibrium adapter.
//!
//! These tests deliberately use a tiny symbolic system rather than database
//! data. They verify that every selected RST strategy receives only symbolic
//! residual expressions and lets RST prepare the Jacobian itself.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, RstSymbolicProblem, RustedSciTheSolveContract, RustedSciTheSolver,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverTermination;
use RustedSciThe::numerical::Nonlinear_systems::symbolic::SymbolicLambdifyFrontend;
use RustedSciThe::numerical::SolverTelemetryMode;
use RustedSciThe::symbolic::symbolic_engine::Expr;

fn scalar_symbolic_problem() -> RstPreparedProblem {
    RstPreparedProblem::new(
        RstSymbolicProblem::from_expressions(
            vec![Expr::parse_expression("x * x - 9.0")],
            vec!["x".to_string()],
            None,
            None,
        )
        .unwrap(),
    )
}

#[test]
fn every_recommended_rst_strategy_solves_the_same_symbolic_problem() {
    let problem = scalar_symbolic_problem();
    let options = RustedSciTheSolveContract::new(1e-8, 250)
        .unwrap()
        .with_telemetry_mode(SolverTelemetryMode::Detailed);

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
        assert!(outcome.metrics.counters_collected, "{}", solver.name());
        assert!(outcome.metrics.timings_collected, "{}", solver.name());
        assert_eq!(outcome.metrics.termination, SolverTermination::Converged);
        assert!(
            outcome.metrics.evaluation_timing.is_some(),
            "{} must expose callback timing for release diagnosis",
            solver.name()
        );
        let timing = outcome.metrics.evaluation_timing.as_ref().unwrap();
        assert!(
            timing.linear_system_micros.is_some(),
            "{} must expose native RST linear-system timing",
            solver.name()
        );
        assert!(
            timing.state_residual_evaluations.is_some()
                && timing.trial_residual_evaluations.is_some()
                && timing.state_jacobian_evaluations.is_some()
                && timing.trial_jacobian_evaluations.is_some(),
            "{} must expose state/trial callback counters",
            solver.name()
        );
        assert!(
            timing.linear_factorization_micros.is_some()
                && timing.linear_system_solve_micros.is_some()
                && timing.accepted_steps.is_some()
                && timing.rejected_steps.is_some(),
            "{} must expose native linear-operation and step counters",
            solver.name()
        );
        assert!(
            timing.iteration_telemetry.is_some(),
            "{} must expose per-iteration RST telemetry",
            solver.name()
        );
    }
}

#[test]
fn rst_telemetry_modes_preserve_independent_counter_and_timing_availability() {
    let problem = scalar_symbolic_problem();
    let solver = RustedSciTheSolver::LevenbergMarquardt;

    for mode in [
        SolverTelemetryMode::Off,
        SolverTelemetryMode::Counters,
        SolverTelemetryMode::Timings,
        SolverTelemetryMode::Detailed,
    ] {
        let contract = RustedSciTheSolveContract::new(1e-8, 250)
            .unwrap()
            .with_telemetry_mode(mode);
        let outcome = solver.solve(&problem, &[1.0], &contract).unwrap();
        assert_eq!(
            outcome.metrics.counters_collected,
            mode.collects_counters(),
            "{mode:?}"
        );
        assert_eq!(
            outcome.metrics.timings_collected,
            mode.collects_timings(),
            "{mode:?}"
        );
        assert_eq!(
            outcome.metrics.evaluation_timing.is_some(),
            mode.collects_timings(),
            "{mode:?}"
        );
        if mode.collects_counters() {
            assert!(outcome.metrics.residual_evaluations > 0, "{mode:?}");
        } else {
            assert_eq!(outcome.metrics.residual_evaluations, 0, "{mode:?}");
        }
        if mode == SolverTelemetryMode::Timings {
            let timing = outcome.metrics.evaluation_timing.as_ref().unwrap();
            assert_eq!(timing.state_residual_evaluations, None);
            assert_eq!(timing.iteration_telemetry, None);
        }
        if mode == SolverTelemetryMode::Detailed {
            assert!(
                outcome
                    .metrics
                    .evaluation_timing
                    .as_ref()
                    .unwrap()
                    .iteration_telemetry
                    .is_some()
            );
        }
    }
}

#[test]
fn atom_view_and_expr_legacy_frontends_match_on_parameterized_residual_and_jacobian() {
    let build = |frontend| {
        RstSymbolicProblem::from_expressions_with_timing_and_frontend(
            vec![Expr::parse_expression("(x - target) / scale")],
            vec!["x".to_string()],
            Some(vec!["target".to_string(), "scale".to_string()]),
            Some(nalgebra::DVector::from_vec(vec![2.5, 0.75])),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Disabled,
            frontend,
        )
        .unwrap()
        .0
    };
    let atom_problem = build(SymbolicLambdifyFrontend::AtomViewNative);
    let legacy_problem = build(SymbolicLambdifyFrontend::ExprLegacy);
    let atom = RstPreparedProblem::new(atom_problem);
    let legacy = RstPreparedProblem::new(legacy_problem);
    let point = [1.25];

    let atom_residual = atom.residual_for_test(&point).unwrap();
    let legacy_residual = legacy.residual_for_test(&point).unwrap();
    let max_residual_difference = atom_residual
        .iter()
        .zip(legacy_residual.iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max);
    let atom_jacobian = atom.jacobian_for_test(&point).unwrap();
    let legacy_jacobian = legacy.jacobian_for_test(&point).unwrap();
    let max_jacobian_difference = atom_jacobian
        .iter()
        .zip(legacy_jacobian.iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max);
    assert!(max_residual_difference <= 1e-12);
    assert!(max_jacobian_difference <= 1e-12);
}

#[test]
fn calculator_diagnostics_map_to_the_four_rst_telemetry_levels() {
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_instrumentation::{
        capture, rst_solver_telemetry_mode,
    };

    assert_eq!(rst_solver_telemetry_mode(false), SolverTelemetryMode::Off);
    assert_eq!(
        rst_solver_telemetry_mode(true),
        SolverTelemetryMode::Timings
    );
    let (_, report) = capture(|| {
        assert_eq!(
            rst_solver_telemetry_mode(false),
            SolverTelemetryMode::Counters
        );
        assert_eq!(
            rst_solver_telemetry_mode(true),
            SolverTelemetryMode::Detailed
        );
    });
    assert!(report.enabled());
}

#[test]
fn prepared_rst_problem_rebinds_numeric_parameters_without_rebuilding_the_graph() {
    let mut problem = RstPreparedProblem::new(
        RstSymbolicProblem::from_expressions(
            vec![Expr::parse_expression(
                "(x - target_enthalpy) / enthalpy_scale",
            )],
            vec!["x".to_string()],
            Some(vec![
                "target_enthalpy".to_string(),
                "enthalpy_scale".to_string(),
            ]),
            Some(nalgebra::DVector::from_vec(vec![2.0, 1.0])),
        )
        .unwrap(),
    );
    let first = problem.residual_for_test(&[0.0]).unwrap()[0];
    problem.set_ph_parameters(5.0, 2.0).unwrap();
    let second = problem.residual_for_test(&[0.0]).unwrap()[0];

    assert!((first + 2.0).abs() < 1e-12);
    assert!((second + 2.5).abs() < 1e-12);
}

#[test]
fn rst_preparation_telemetry_reports_detailed_backend_stages() {
    let (_, report) = RstSymbolicProblem::from_expressions_with_timing(
        vec![Expr::parse_expression("x * x - target")],
        vec!["x".to_string()],
        Some(vec!["target".to_string()]),
        Some(nalgebra::DVector::from_vec(vec![9.0])),
        crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Enabled,
    )
    .expect("prepared RST telemetry fixture must be valid");

    assert!(report.enabled());
    assert!(report.backend_stage_breakdown_available());
    assert_eq!(report.backend_kind().as_str(), "lambdify");
    assert_eq!(report.artifact_action().as_str(), "not-applicable");
    assert!(report.backend_build_duration().is_none());
    assert!(report.generated_residual_jobs().is_none());
    assert!(report.generated_jacobian_jobs().is_none());
    assert!(report.total() >= report.problem_materialization());
    assert!(report.symbolic_differentiation() > std::time::Duration::ZERO);
    assert!(report.residual_lambdification() > std::time::Duration::ZERO);
    assert!(report.jacobian_lambdification() > std::time::Duration::ZERO);
    assert!(report.rst_input_validation().is_some());
    assert!(report.rst_problem_assembly().is_some());
    assert!(report.rst_unattributed().is_some());

    let (_, atom_report) = RstSymbolicProblem::from_expressions_with_timing_and_frontend(
        vec![Expr::parse_expression("x * x - target")],
        vec!["x".to_string()],
        Some(vec!["target".to_string()]),
        Some(nalgebra::DVector::from_vec(vec![9.0])),
        crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Enabled,
        SymbolicLambdifyFrontend::AtomViewNative,
    )
    .expect("AtomView RST preparation telemetry fixture must be valid");
    assert!(atom_report.symbolic_differentiation() > std::time::Duration::ZERO);
    assert!(atom_report.residual_lambdification() > std::time::Duration::ZERO);
    assert!(atom_report.jacobian_lambdification() > std::time::Duration::ZERO);
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
