use RustedSciThe::numerical::StiffAPI::{StiffFrontend, StiffMatrixLayout, StiffSolverKind};
use RustedSciThe::numerical::logging::SolverLogLevel;
use RustedSciThe::numerical::{CancellationToken, SolverTelemetryMode};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};

use super::*;

fn species(id: &str, initial_concentration: f64) -> SpeciesSpec {
    SpeciesSpec {
        id: id.to_owned(),
        initial_concentration,
    }
}

fn reaction(
    id: &str,
    reactants: Vec<(&str, f64)>,
    products: Vec<(&str, f64)>,
    rate_law: RateLaw,
) -> ReactionSpec {
    ReactionSpec {
        id: id.to_owned(),
        reactants: reactants
            .into_iter()
            .map(|(s, n)| StoichTerm::new(s, n))
            .collect(),
        products: products
            .into_iter()
            .map(|(s, n)| StoichTerm::new(s, n))
            .collect(),
        rate_law,
    }
}

fn a_to_two_b() -> NetworkProblem {
    NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![reaction(
            "r1",
            vec![("A", 1.0)],
            vec![("B", 2.0)],
            RateLaw::MassAction {
                coefficient: RateCoefficient::Constant(0.5),
            },
        )],
        TemperatureProgram::Constant(700.0),
        0.0,
        2.0,
    )
    .unwrap()
}

#[test]
fn empty_reaction_set_integrates_inert_species_without_changing_state() {
    let problem = NetworkProblem::new(
        vec![species("inert-A", 1.25), species("inert-B", 0.0)],
        Vec::new(),
        TemperatureProgram::Constant(400.0),
        0.0,
        3.0,
    )
    .unwrap();
    let solution = solve_network(&problem, &NetworkIvpOptions::default()).unwrap();
    assert_eq!(solution.species_ids, ["inert-A", "inert-B"]);
    assert!(
        solution
            .concentrations
            .iter()
            .all(|state| state == &[1.25, 0.0])
    );
    assert!(solution.time_seconds.len() >= 2);
}

#[test]
fn mass_action_rhs_uses_reaction_by_species_stoichiometry() {
    let problem = NetworkProblem::new(
        vec![
            species("A", 2.0),
            species("B", 1.0),
            species("C", 3.0),
            species("D", 0.0),
        ],
        vec![
            reaction(
                "r1",
                vec![("A", 1.0)],
                vec![("B", 2.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.5),
                },
            ),
            reaction(
                "r2",
                vec![("B", 1.0), ("C", 1.0)],
                vec![("D", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.1),
                },
            ),
        ],
        TemperatureProgram::Constant(700.0),
        0.0,
        1.0,
    )
    .unwrap();
    let rhs = problem.evaluate_rhs(0.0, &[2.0, 1.0, 3.0, 0.0]).unwrap();
    let expected = [-1.0, 1.7, -0.3, 0.3];
    for (actual, expected) in rhs.iter().zip(expected) {
        assert!((actual - expected).abs() < 1e-12);
    }
    assert_eq!(
        problem.assemble().stoichiometry,
        vec![vec![-1.0, 2.0, 0.0, 0.0], vec![0.0, -1.0, -1.0, 1.0]]
    );
}

#[test]
fn custom_rate_substitutes_prescribed_temperature_once() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![reaction(
            "r",
            vec![("A", 1.0)],
            vec![("B", 1.0)],
            RateLaw::CustomExpr(
                Expr::Var("T".into()) * Expr::Var("c0".into()) / Expr::Const(500.0),
            ),
        )],
        TemperatureProgram::Linear {
            initial_kelvin: 500.0,
            kelvin_per_second: 100.0,
        },
        0.0,
        2.0,
    )
    .unwrap();
    let rhs = problem.evaluate_rhs(1.0, &[2.0, 0.0]).unwrap();
    assert!((rhs[0] + 2.4).abs() < 1e-12);
    assert!((rhs[1] - 2.4).abs() < 1e-12);
}

#[test]
fn empirical_orders_do_not_change_stoichiometric_balance() {
    let problem = NetworkProblem::new(
        vec![species("A", 2.0), species("B", 0.0)],
        vec![reaction(
            "empirical",
            vec![("A", 1.0)],
            vec![("B", 2.0)],
            RateLaw::EmpiricalMassAction {
                coefficient: RateCoefficient::Constant(0.5),
                orders: vec![KineticOrder::new("A", 2.0)],
            },
        )],
        TemperatureProgram::Constant(500.0),
        0.0,
        1.0,
    )
    .unwrap();
    assert_eq!(problem.assemble().stoichiometry, vec![vec![-1.0, 2.0]]);
    let rhs = problem.evaluate_rhs(0.0, &[2.0, 0.0]).unwrap();
    assert!((rhs[0] + 2.0).abs() < 1e-12);
    assert!((rhs[1] - 4.0).abs() < 1e-12);
}

#[test]
fn assembled_state_jacobian_matches_independent_rhs_finite_differences() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.7), species("B", 0.8), species("C", 0.4)],
        vec![
            reaction(
                "arrhenius-mass-action",
                vec![("A", 1.0)],
                vec![("B", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Arrhenius {
                        prefactor: 2.3,
                        activation_energy: 1200.0,
                    },
                },
            ),
            reaction(
                "empirical-order",
                vec![("B", 1.0)],
                vec![("C", 1.0)],
                RateLaw::EmpiricalMassAction {
                    coefficient: RateCoefficient::Constant(0.7),
                    orders: vec![KineticOrder::new("B", 1.5)],
                },
            ),
            reaction(
                "custom-rate",
                vec![("A", 1.0), ("C", 1.0)],
                vec![("B", 1.0)],
                RateLaw::CustomExpr(
                    Expr::Var("c0".into()) * Expr::Var("c2".into()) * Expr::Var("T".into())
                        / Expr::Const(1_000.0),
                ),
            ),
        ],
        TemperatureProgram::Linear {
            initial_kelvin: 650.0,
            kelvin_per_second: 35.0,
        },
        0.0,
        1.0,
    )
    .unwrap();

    let time = 0.37;
    let state = [1.2, 0.9, 0.55];
    let assembled = problem.assemble();
    let mut symbols: Vec<&str> = assembled.variables.iter().map(String::as_str).collect();
    symbols.push("t");
    symbols.extend(assembled.parameter_names.iter().map(String::as_str));
    let mut values = state.to_vec();
    values.push(time);
    values.extend_from_slice(&assembled.parameter_values);

    for (row, equation) in assembled.equations.iter().enumerate() {
        for (column, variable) in assembled.variables.iter().enumerate() {
            let symbolic = equation.diff(variable).eval_expression(&symbols, &values);
            let step = state[column].abs().max(1.0) * 1e-5;
            let mut plus = state;
            let mut minus = state;
            plus[column] += step;
            minus[column] -= step;
            let rhs_plus = problem.evaluate_rhs(time, &plus).unwrap()[row];
            let rhs_minus = problem.evaluate_rhs(time, &minus).unwrap()[row];
            let finite_difference = (rhs_plus - rhs_minus) / (2.0 * step);
            let tolerance = 2e-7 * symbolic.abs().max(finite_difference.abs()).max(1.0);
            assert!(
                (symbolic - finite_difference).abs() <= tolerance,
                "d(rhs[{row}])/d({variable}): symbolic={symbolic:.12e}, finite_difference={finite_difference:.12e}, tolerance={tolerance:.3e}"
            );
        }
    }
}

#[test]
fn explicit_reverse_reactions_preserve_species_order_and_approach_equilibrium() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![
            reaction(
                "forward",
                vec![("A", 1.0)],
                vec![("B", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(2.0),
                },
            ),
            reaction(
                "reverse",
                vec![("B", 1.0)],
                vec![("A", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(1.0),
                },
            ),
        ],
        TemperatureProgram::Constant(600.0),
        0.0,
        5.0,
    )
    .unwrap();
    assert_eq!(
        problem.assemble().stoichiometry,
        vec![vec![-1.0, 1.0], vec![1.0, -1.0]]
    );
    let solution = solve_network(&problem, &NetworkIvpOptions::default()).unwrap();
    let end = solution.concentrations.last().unwrap();
    let transient = (-3.0_f64 * 5.0).exp();
    let expected_a = 1.0 / 3.0 + 2.0 / 3.0 * transient;
    let expected_b = 2.0 / 3.0 * (1.0 - transient);
    assert!((end[0] - expected_a).abs() < 2e-5);
    assert!((end[1] - expected_b).abs() < 2e-5);
    for state in &solution.concentrations {
        assert!((state[0] + state[1] - 1.0).abs() < 2e-6);
    }
}

#[test]
fn rhs_domain_error_identifies_the_nonfinite_state_component() {
    let problem = NetworkProblem::new(
        vec![species("A", 0.0), species("B", 0.0)],
        vec![reaction(
            "singular",
            vec![],
            vec![("B", 1.0)],
            RateLaw::CustomExpr(Expr::Const(1.0) / Expr::Var("c0".into())),
        )],
        TemperatureProgram::Constant(500.0),
        0.0,
        1.0,
    )
    .unwrap();
    let error = problem.evaluate_rhs(0.25, &[0.0, 0.0]).unwrap_err();
    assert!(
        matches!(error, NetworkIvpError::RateDomain { time, component: 1, species_id }
        if time == 0.25 && species_id == "B")
    );
}

#[test]
fn validation_keeps_reaction_context_and_rejects_unknown_symbols() {
    let error = NetworkProblem::new(
        vec![species("A", 1.0)],
        vec![reaction(
            "bad-stage",
            vec![("A", 1.0)],
            vec![],
            RateLaw::CustomExpr(Expr::Var("undefined".into())),
        )],
        TemperatureProgram::Constant(300.0),
        0.0,
        1.0,
    )
    .unwrap_err();
    assert!(
        matches!(error, NetworkIvpError::RateSymbol { reaction_id, symbol }
        if reaction_id == "bad-stage" && symbol == "undefined")
    );

    let error = NetworkProblem::new(
        vec![species("A", 1.0)],
        vec![reaction(
            "bad-ref",
            vec![("missing", 1.0)],
            vec![],
            RateLaw::MassAction {
                coefficient: RateCoefficient::Constant(1.0),
            },
        )],
        TemperatureProgram::Constant(300.0),
        0.0,
        1.0,
    )
    .unwrap_err();
    assert!(
        matches!(error, NetworkIvpError::UnknownSpecies { reaction_id, species_id }
        if reaction_id == "bad-ref" && species_id == "missing")
    );

    let error = NetworkProblem::new(
        vec![species("A", 1.0)],
        vec![reaction(
            "bad-expression",
            vec![("A", 1.0)],
            vec![],
            RateLaw::CustomExpr(Expr::Const(f64::NAN)),
        )],
        TemperatureProgram::Constant(300.0),
        0.0,
        1.0,
    )
    .unwrap_err();
    assert!(matches!(error, NetworkIvpError::InvalidInput { field, .. }
        if field == "reactions[0].rate_law.expression"));

    let failed = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![
            reaction(
                "valid-before",
                vec![("A", 1.0)],
                vec![("B", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(1.0),
                },
            ),
            reaction(
                "malformed-middle",
                vec![("B", 1.0)],
                vec![("A", 1.0)],
                RateLaw::CustomExpr(Expr::Var("not_declared".into())),
            ),
            reaction(
                "valid-after",
                vec![("B", 1.0)],
                vec![],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.1),
                },
            ),
        ],
        TemperatureProgram::Constant(300.0),
        0.0,
        1.0,
    );
    assert!(matches!(
        failed,
        Err(NetworkIvpError::RateSymbol { reaction_id, symbol })
            if reaction_id == "malformed-middle" && symbol == "not_declared"
    ));
}

#[test]
fn bdf_facade_solves_analytical_a_to_two_b_without_telemetry_by_default() {
    let problem = a_to_two_b();
    let options = NetworkIvpOptions::default();
    assert_eq!(options.telemetry, SolverTelemetryMode::Off);
    assert_eq!(options.logging, SolverLogLevel::Off);
    let solution = solve_network(&problem, &options).unwrap();
    assert!(solution.telemetry.is_none());
    assert_eq!(solution.species_ids, ["A", "B"]);
    let final_state = solution.concentrations.last().unwrap();
    let expected_a = (-1.0_f64).exp();
    assert!((final_state[0] - expected_a).abs() < 1e-5);
    assert!((final_state[1] - 2.0 * (1.0 - expected_a)).abs() < 2e-5);
    assert!((2.0 * final_state[0] + final_state[1] - 2.0).abs() < 2e-5);
    for (&time, state) in solution.time_seconds.iter().zip(&solution.concentrations) {
        let analytical_a = (-0.5 * time).exp();
        assert!((state[0] - analytical_a).abs() < 2e-5);
        assert!((state[1] - 2.0 * (1.0 - analytical_a)).abs() < 4e-5);
        assert!((2.0 * state[0] + state[1] - 2.0).abs() < 2e-5);
    }
    assert_eq!(
        solution.temperature_kelvin.len(),
        solution.time_seconds.len()
    );
    assert!(
        solution
            .temperature_kelvin
            .iter()
            .all(|value| *value == 700.0)
    );
}

#[test]
fn native_logging_level_is_forwarded_to_all_stiffapi_families() {
    let problem = a_to_two_b();
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let options = NetworkIvpOptions {
            family,
            logging: SolverLogLevel::Debug,
            ..Default::default()
        };
        let native = super::solver_backend::prepare_native(&problem, problem.assemble(), &options)
            .unwrap_or_else(|error| panic!("{family:?}: {error}"));
        let effective = match family {
            StiffSolverKind::Bdf => native.as_bdf().unwrap().logging_level(),
            StiffSolverKind::Radau => native.as_radau().unwrap().logging_level(),
            StiffSolverKind::BackwardEuler => native.as_backward_euler().unwrap().logging_level(),
            StiffSolverKind::Lsode2 => native.as_lsode2().unwrap().logging_level(),
        };
        assert_eq!(effective, SolverLogLevel::Debug, "{family:?}");
        assert_eq!(
            native.telemetry_mode(),
            SolverTelemetryMode::Off,
            "{family:?}"
        );
    }
}

#[test]
fn bdf_requested_counters_return_a_native_snapshot() {
    let problem = a_to_two_b();
    let options = NetworkIvpOptions {
        telemetry: SolverTelemetryMode::Counters,
        family: StiffSolverKind::Bdf,
        ..Default::default()
    };
    let solution = solve_network(&problem, &options).unwrap();
    assert!(matches!(
        solution.telemetry,
        Some(
            RustedSciThe::numerical::StiffAPI::StiffTelemetrySnapshot::Bdf {
                mode: SolverTelemetryMode::Counters,
                ..
            }
        )
    ));
}

#[test]
fn telemetry_modes_are_forwarded_independently_to_all_stiffapi_families() {
    let problem = a_to_two_b();
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        for telemetry in [
            SolverTelemetryMode::Off,
            SolverTelemetryMode::Counters,
            SolverTelemetryMode::Timings,
            SolverTelemetryMode::Detailed,
        ] {
            let options = NetworkIvpOptions {
                family,
                telemetry,
                ..Default::default()
            };
            let solution = solve_network(&problem, &options)
                .unwrap_or_else(|error| panic!("{family:?}/{telemetry:?}: {error}"));
            assert_eq!(
                solution.telemetry.is_some(),
                telemetry != SolverTelemetryMode::Off,
                "{family:?}/{telemetry:?}"
            );
            match (&solution.stage_diagnostics, telemetry) {
                (None, SolverTelemetryMode::Off) => {}
                (Some(_), SolverTelemetryMode::Off) => {
                    panic!("KiThe stages unexpectedly collected in Off mode for {family:?}")
                }
                (Some(stages), SolverTelemetryMode::Counters) => {
                    assert_eq!(stages.mode, telemetry);
                    assert_eq!(stages.equation_assembly.calls, Some(1));
                    assert_eq!(stages.backend_preparation.calls, Some(1));
                    assert_eq!(stages.integration.calls, Some(1));
                    assert_eq!(stages.output_postprocessing.calls, Some(1));
                    assert!(stages.equation_assembly.elapsed.is_none());
                    assert!(stages.integration.elapsed.is_none());
                }
                (Some(stages), SolverTelemetryMode::Timings) => {
                    assert_eq!(stages.mode, telemetry);
                    assert!(stages.equation_assembly.calls.is_none());
                    assert!(stages.integration.calls.is_none());
                    assert!(stages.equation_assembly.elapsed.is_some());
                    assert!(stages.backend_preparation.elapsed.is_some());
                    assert!(stages.integration.elapsed.is_some());
                    assert!(stages.output_postprocessing.elapsed.is_some());
                }
                (Some(stages), SolverTelemetryMode::Detailed) => {
                    assert_eq!(stages.mode, telemetry);
                    assert_eq!(stages.equation_assembly.calls, Some(1));
                    assert_eq!(stages.integration.calls, Some(1));
                    assert!(stages.equation_assembly.elapsed.is_some());
                    assert!(stages.integration.elapsed.is_some());
                }
                (None, _) => panic!("missing KiThe stages for {family:?}/{telemetry:?}"),
            }
        }
    }
}

#[test]
fn four_stiffapi_families_solve_the_same_analytical_network() {
    let problem = a_to_two_b();
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let options = NetworkIvpOptions {
            family,
            max_step: 0.1,
            ..Default::default()
        };
        let solution = solve_network(&problem, &options)
            .unwrap_or_else(|error| panic!("{family:?} failed on A -> 2B: {error}"));
        let final_state = solution.concentrations.last().unwrap();
        let expected = (-1.0_f64).exp();
        assert!((final_state[0] - expected).abs() < 2e-2, "{family:?}");
        assert!(
            (final_state[1] - 2.0 * (1.0 - expected)).abs() < 4e-2,
            "{family:?}"
        );
        assert!(
            (2.0 * final_state[0] + final_state[1] - 2.0).abs() < 2e-2,
            "{family:?}"
        );
    }
}

#[test]
fn nonlinear_custom_rate_network_matches_logistic_reference_across_frontends_and_families() {
    let total = 1.02_f64;
    let initial_a = 1.0_f64;
    let initial_b = 0.02_f64;
    let rate = 3.0_f64;
    let problem = NetworkProblem::new(
        vec![species("A", initial_a), species("B", initial_b)],
        vec![reaction(
            "autocatalytic-transfer",
            vec![("A", 1.0)],
            vec![("B", 1.0)],
            RateLaw::CustomExpr(Expr::parse_expression("3*c0*c1")),
        )],
        TemperatureProgram::Constant(700.0),
        0.0,
        0.1,
    )
    .unwrap();

    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        for frontend in [StiffFrontend::ExprLegacy, StiffFrontend::AtomViewNative] {
            let mut options = NetworkIvpOptions {
                family,
                frontend,
                max_step: 0.0002,
                ..Default::default()
            };
            if family != StiffSolverKind::BackwardEuler {
                options.rtol = 1e-9;
                options.atol = 1e-11;
            }
            let solution = solve_network(&problem, &options)
                .unwrap_or_else(|error| panic!("{family:?}/{frontend:?}: {error}"));
            for (&time, state) in solution.time_seconds.iter().zip(&solution.concentrations) {
                let reference_b =
                    total / (1.0 + (initial_a / initial_b) * (-rate * total * time).exp());
                let reference_a = total - reference_b;
                assert!(
                    (state[0] - reference_a).abs() < 1e-4,
                    "{family:?}/{frontend:?} A(t={time}): {} != {reference_a}",
                    state[0]
                );
                assert!(
                    (state[1] - reference_b).abs() < 1e-4,
                    "{family:?}/{frontend:?} B(t={time}): {} != {reference_b}",
                    state[1]
                );
                assert!((state[0] + state[1] - total).abs() < 2e-8);
            }
        }
    }
}

#[test]
fn frontend_and_family_matrix_preserves_first_order_solution_across_scales() {
    for scale in [1e-3, 1.0, 1e3] {
        let problem = NetworkProblem::new(
            vec![species("A", scale), species("B", 0.0)],
            vec![reaction(
                "first-order",
                vec![("A", 1.0)],
                vec![("B", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.5),
                },
            )],
            TemperatureProgram::Constant(700.0),
            0.0,
            2.0,
        )
        .unwrap();
        for family in [
            StiffSolverKind::Bdf,
            StiffSolverKind::Radau,
            StiffSolverKind::BackwardEuler,
            StiffSolverKind::Lsode2,
        ] {
            let mut endpoints = Vec::new();
            for frontend in [StiffFrontend::ExprLegacy, StiffFrontend::AtomViewNative] {
                let solution = solve_network(
                    &problem,
                    &NetworkIvpOptions {
                        family,
                        frontend,
                        max_step: 0.02,
                        ..Default::default()
                    },
                )
                .unwrap_or_else(|error| panic!("{family:?}/{frontend:?}/scale={scale}: {error}"));
                let end = solution.concentrations.last().unwrap();
                let expected_a = (-1.0_f64).exp();
                let expected_b = 1.0 - expected_a;
                assert!(
                    (end[0] / scale - expected_a).abs() < 2e-2,
                    "{family:?}/{frontend:?}/scale={scale}: A={}",
                    end[0]
                );
                assert!(
                    (end[1] / scale - expected_b).abs() < 2e-2,
                    "{family:?}/{frontend:?}/scale={scale}: B={}",
                    end[1]
                );
                endpoints.push(end.clone());
            }
            for (index, (expr, atom)) in endpoints[0].iter().zip(&endpoints[1]).enumerate() {
                assert!(
                    (expr - atom).abs() <= 2e-5 * scale,
                    "{family:?}/scale={scale}/component={index}: ExprLegacy={expr:.9e}, AtomView={atom:.9e}"
                );
            }
        }
    }
}

#[test]
fn stiff_consecutive_first_order_network_matches_closed_form_for_all_families() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0), species("C", 0.0)],
        vec![
            reaction(
                "fast-A-to-B",
                vec![("A", 1.0)],
                vec![("B", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(100.0),
                },
            ),
            reaction(
                "slow-B-to-C",
                vec![("B", 1.0)],
                vec![("C", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.2),
                },
            ),
        ],
        TemperatureProgram::Constant(700.0),
        0.0,
        0.1,
    )
    .unwrap();
    let time: f64 = 0.1;
    let fast = 100.0_f64;
    let slow = 0.2_f64;
    let expected_a = (-fast * time).exp();
    let expected_b = fast / (fast - slow) * ((-slow * time).exp() - (-fast * time).exp());
    let expected_c = 1.0 - expected_a - expected_b;

    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let solution = solve_network(
            &problem,
            &NetworkIvpOptions {
                family,
                max_step: 0.005,
                ..Default::default()
            },
        )
        .unwrap_or_else(|error| panic!("{family:?} failed on stiff chain: {error}"));
        let end = solution.concentrations.last().unwrap();
        assert!(
            (end[0] - expected_a).abs() < 2e-3,
            "{family:?}: A={:.9e}, expected={expected_a:.9e}",
            end[0]
        );
        assert!(
            (end[1] - expected_b).abs() < 2e-3,
            "{family:?}: B={:.9e}, expected={expected_b:.9e}",
            end[1]
        );
        assert!(
            (end[2] - expected_c).abs() < 2e-3,
            "{family:?}: C={:.9e}, expected={expected_c:.9e}",
            end[2]
        );
        for state in &solution.concentrations {
            assert!(state.iter().all(|value| value.is_finite() && *value >= 0.0));
            assert!((state.iter().sum::<f64>() - 1.0).abs() < 2e-5);
        }
    }
}

#[test]
fn coupled_two_reaction_network_preserves_invariant_across_stiffapi_families() {
    let problem = NetworkProblem::new(
        vec![
            species("A", 1.0),
            species("B", 0.0),
            species("C", 1.0),
            species("D", 0.0),
        ],
        vec![
            reaction(
                "A-to-two-B",
                vec![("A", 1.0)],
                vec![("B", 2.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.5),
                },
            ),
            reaction(
                "B-plus-C-to-D",
                vec![("B", 1.0), ("C", 1.0)],
                vec![("D", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.1),
                },
            ),
        ],
        TemperatureProgram::Constant(700.0),
        0.0,
        2.0,
    )
    .unwrap();
    let mut endpoints = Vec::new();
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let solution = solve_network(
            &problem,
            &NetworkIvpOptions {
                family,
                max_step: 0.005,
                ..Default::default()
            },
        )
        .unwrap_or_else(|error| panic!("{family:?} failed: {error}"));
        for state in &solution.concentrations {
            assert!(state.iter().all(|value| value.is_finite() && *value >= 0.0));
            assert!(
                (state[2] + state[3] - 1.0).abs() < 2e-6,
                "{family:?} violated C+D=1: {state:?}"
            );
        }
        endpoints.push((family, solution.concentrations.last().unwrap().clone()));
    }

    let reference = &endpoints[0].1;
    for (family, endpoint) in endpoints.iter().skip(1) {
        for (index, (actual, expected)) in endpoint.iter().zip(reference).enumerate() {
            assert!(
                (actual - expected).abs() < 2e-3,
                "{family:?} endpoint component {index}={actual:.9e}, BDF={expected:.9e}"
            );
        }
    }
}

#[test]
fn atom_view_and_expr_legacy_publish_the_same_bdf_solution() {
    let problem = a_to_two_b();
    let expr = solve_network(&problem, &NetworkIvpOptions::default()).unwrap();
    let atom = solve_network(
        &problem,
        &NetworkIvpOptions {
            frontend: StiffFrontend::AtomViewNative,
            ..Default::default()
        },
    )
    .unwrap();
    let expr_final = expr.concentrations.last().unwrap();
    let atom_final = atom.concentrations.last().unwrap();
    for (left, right) in expr_final.iter().zip(atom_final) {
        assert!((left - right).abs() < 1e-8);
    }
}

#[test]
fn unsupported_native_matrix_route_is_typed_preparation_error() {
    let problem = a_to_two_b();
    let options = NetworkIvpOptions {
        family: StiffSolverKind::BackwardEuler,
        matrix_layout: StiffMatrixLayout::Sparse,
        ..Default::default()
    };
    let error = solve_network(&problem, &options).unwrap_err();
    assert!(matches!(error, NetworkIvpError::Preparation(_)));
    assert!(std::error::Error::source(&error).is_some());
}

#[test]
fn backward_euler_does_not_silently_ignore_requested_integration_tolerance() {
    let problem = a_to_two_b();
    let options = NetworkIvpOptions {
        family: StiffSolverKind::BackwardEuler,
        rtol: 1e-5,
        ..Default::default()
    };
    let error = solve_network(&problem, &options).unwrap_err();
    assert!(matches!(
        error,
        NetworkIvpError::UnsupportedCapability {
            family: StiffSolverKind::BackwardEuler,
            ..
        }
    ));
}

#[test]
fn nonautonomous_custom_rate_matches_independent_integral() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![reaction(
            "temperature-driven",
            vec![("A", 1.0)],
            vec![("B", 1.0)],
            RateLaw::CustomExpr(
                Expr::Var("c0".into()) * Expr::Var("T".into()) / Expr::Const(500.0),
            ),
        )],
        TemperatureProgram::Symbolic(
            Expr::Const(500.0) + Expr::Const(100.0) * Expr::Var("t".into()),
        ),
        0.0,
        2.0,
    )
    .unwrap();
    let solution = solve_network(&problem, &NetworkIvpOptions::default()).unwrap();
    for ((&time, state), &temperature) in solution
        .time_seconds
        .iter()
        .zip(&solution.concentrations)
        .zip(&solution.temperature_kelvin)
    {
        let expected = (-(time + 0.1 * time * time)).exp();
        assert!((state[0] - expected).abs() < 2e-5);
        assert!((state[1] - (1.0 - expected)).abs() < 2e-5);
        assert!((temperature - (500.0 + 100.0 * time)).abs() < 1e-12);
    }
}

#[test]
fn nonlinear_temperature_rate_matches_closed_form_across_frontends_and_families() {
    let initial_a = 1.0_f64;
    let rate_scale = 0.4_f64;
    let temperature_scale = 500.0_f64;
    let initial_temperature = 500.0_f64;
    let heating_rate = 100.0_f64;
    let problem = NetworkProblem::new(
        vec![species("A", initial_a), species("B", 0.0)],
        vec![reaction(
            "quadratic-temperature-rate",
            vec![("A", 1.0)],
            vec![("B", 1.0)],
            RateLaw::CustomExpr(
                Expr::Const(rate_scale)
                    * Expr::Var("c0".into())
                    * (Expr::Var("T".into()) / Expr::Const(temperature_scale))
                    * (Expr::Var("T".into()) / Expr::Const(temperature_scale)),
            ),
        )],
        TemperatureProgram::Linear {
            initial_kelvin: initial_temperature,
            kelvin_per_second: heating_rate,
        },
        0.0,
        1.0,
    )
    .unwrap();

    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        for frontend in [StiffFrontend::ExprLegacy, StiffFrontend::AtomViewNative] {
            let mut options = NetworkIvpOptions {
                family,
                frontend,
                max_step: 2e-4,
                ..Default::default()
            };
            if family != StiffSolverKind::BackwardEuler {
                options.rtol = 1e-9;
                options.atol = 1e-11;
            }
            let solution = solve_network(&problem, &options)
                .unwrap_or_else(|error| panic!("{family:?}/{frontend:?}: {error}"));
            for ((&time, state), &temperature) in solution
                .time_seconds
                .iter()
                .zip(&solution.concentrations)
                .zip(&solution.temperature_kelvin)
            {
                let beta = heating_rate / temperature_scale;
                let integrated_rate =
                    rate_scale * (time + beta * time.powi(2) + beta.powi(2) * time.powi(3) / 3.0);
                let reference_a = initial_a * (-integrated_rate).exp();
                let reference_b = initial_a - reference_a;
                assert!(
                    (state[0] - reference_a).abs() < 1e-4,
                    "{family:?}/{frontend:?} A(t={time}): {} != {reference_a}",
                    state[0]
                );
                assert!(
                    (state[1] - reference_b).abs() < 1e-4,
                    "{family:?}/{frontend:?} B(t={time}): {} != {reference_b}",
                    state[1]
                );
                assert!((state[0] + state[1] - initial_a).abs() < 2e-8);
                assert!((temperature - initial_temperature - heating_rate * time).abs() < 1e-10);
            }
        }
    }
}

#[test]
fn stiff_nonlinear_logistic_network_matches_closed_form_across_scales() {
    let total_scaled = 1.02_f64;
    let initial_product_fraction = 0.02_f64;
    let fast_rate = 100.0_f64;
    let slow_rate = 0.1_f64;

    for concentration_scale in [1e-3_f64, 1.0, 1e3] {
        let problem = NetworkProblem::new(
            vec![
                species("A_fast", concentration_scale),
                species("B_fast", initial_product_fraction * concentration_scale),
                species("A_slow", concentration_scale),
                species("B_slow", initial_product_fraction * concentration_scale),
            ],
            vec![
                reaction(
                    "fast-autocatalytic-transfer",
                    vec![("A_fast", 1.0)],
                    vec![("B_fast", 1.0)],
                    RateLaw::CustomExpr(
                        Expr::Const(fast_rate / concentration_scale)
                            * Expr::Var("c0".into())
                            * Expr::Var("c1".into()),
                    ),
                ),
                reaction(
                    "slow-autocatalytic-transfer",
                    vec![("A_slow", 1.0)],
                    vec![("B_slow", 1.0)],
                    RateLaw::CustomExpr(
                        Expr::Const(slow_rate / concentration_scale)
                            * Expr::Var("c2".into())
                            * Expr::Var("c3".into()),
                    ),
                ),
            ],
            TemperatureProgram::Constant(700.0),
            0.0,
            0.1,
        )
        .unwrap();

        for family in [
            StiffSolverKind::Bdf,
            StiffSolverKind::Radau,
            StiffSolverKind::BackwardEuler,
            StiffSolverKind::Lsode2,
        ] {
            for frontend in [StiffFrontend::ExprLegacy, StiffFrontend::AtomViewNative] {
                let mut options = NetworkIvpOptions {
                    family,
                    frontend,
                    max_step: if family == StiffSolverKind::BackwardEuler {
                        2e-6
                    } else {
                        2e-3
                    },
                    ..Default::default()
                };
                if family != StiffSolverKind::BackwardEuler {
                    options.rtol = 1e-9;
                    options.atol = 1e-11;
                }
                let solution = solve_network(&problem, &options).unwrap_or_else(|error| {
                    panic!("{family:?}/{frontend:?}/scale={concentration_scale}: {error}")
                });

                for (&time, state) in solution.time_seconds.iter().zip(&solution.concentrations) {
                    for (pair, reference_rate) in [(0, fast_rate), (2, slow_rate)] {
                        let reference_b = total_scaled
                            / (1.0
                                + (1.0 / initial_product_fraction)
                                    * (-reference_rate * total_scaled * time).exp());
                        let reference_a = total_scaled - reference_b;
                        let normalized_a = state[pair] / concentration_scale;
                        let normalized_b = state[pair + 1] / concentration_scale;
                        assert!(
                            normalized_a.is_finite() && normalized_a >= 0.0,
                            "{family:?}/{frontend:?}/scale={concentration_scale} A[{pair}] at t={time}: {normalized_a}"
                        );
                        assert!(
                            normalized_b.is_finite() && normalized_b >= 0.0,
                            "{family:?}/{frontend:?}/scale={concentration_scale} B[{pair}] at t={time}: {normalized_b}"
                        );
                        assert!(
                            (normalized_a - reference_a).abs() < 2e-4,
                            "{family:?}/{frontend:?}/scale={concentration_scale} A[{pair}] at t={time}: {normalized_a} != {reference_a}"
                        );
                        assert!(
                            (normalized_b - reference_b).abs() < 2e-4,
                            "{family:?}/{frontend:?}/scale={concentration_scale} B[{pair}] at t={time}: {normalized_b} != {reference_b}"
                        );
                        assert!(
                            (normalized_a + normalized_b - total_scaled).abs() < 2e-5,
                            "{family:?}/{frontend:?}/scale={concentration_scale} invariant pair={pair} at t={time}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn prepared_session_reuses_all_four_families_and_retains_last_good_result() {
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let options = NetworkIvpOptions {
            family,
            max_step: 0.1,
            telemetry: SolverTelemetryMode::Counters,
            ..Default::default()
        };
        let problem = a_to_two_b();
        let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
        let first_a = prepared.solve().unwrap().concentrations.last().unwrap()[0];
        assert!((first_a - (-1.0_f64).exp()).abs() < 2e-2, "{family:?}");
        let second = prepared
            .solve_with_initial_concentrations(&[2.0, 0.5])
            .unwrap();
        let second_final = second.concentrations.last().unwrap();
        let expected_a = 2.0 * (-1.0_f64).exp();
        let expected_b = 0.5 + 2.0 * (2.0 - expected_a);
        assert!((second_final[0] - expected_a).abs() < 2e-2, "{family:?}");
        assert!((second_final[1] - expected_b).abs() < 4e-2, "{family:?}");
        assert!(prepared.last_attempt_telemetry().is_some(), "{family:?}");
        let stages = prepared.stage_diagnostics().unwrap();
        assert_eq!(stages.equation_assembly.calls, Some(1), "{family:?}");
        assert_eq!(stages.backend_preparation.calls, Some(1), "{family:?}");
        assert_eq!(stages.integration.calls, Some(2), "{family:?}");
        assert_eq!(stages.output_postprocessing.calls, Some(2), "{family:?}");

        let accepted = prepared
            .last_solution()
            .unwrap()
            .concentrations
            .last()
            .unwrap()
            .clone();
        let error = prepared
            .solve_with_initial_concentrations(&[-1.0, 0.5])
            .unwrap_err();
        assert!(matches!(error, NetworkIvpError::InvalidInput { .. }));
        assert!(prepared.last_attempt_telemetry().is_none());
        assert!(prepared.last_failed_trajectory().is_none());
        assert_eq!(
            prepared
                .last_solution()
                .unwrap()
                .concentrations
                .last()
                .unwrap(),
            &accepted
        );
    }
}

#[test]
fn cancellation_does_not_replace_an_accepted_solution() {
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let problem = a_to_two_b();
        let options = NetworkIvpOptions {
            family,
            max_step: 0.1,
            telemetry: SolverTelemetryMode::Counters,
            ..Default::default()
        };
        let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
        let accepted = prepared
            .solve()
            .unwrap()
            .concentrations
            .last()
            .unwrap()
            .clone();
        let token = CancellationToken::new();
        token.cancel();
        let error = prepared
            .solve_with_cancellation_token(&[2.0, 0.5], token)
            .unwrap_err();
        assert!(
            matches!(error, NetworkIvpError::CancellationAwareIntegration(_)),
            "{family:?}: {error}"
        );
        assert!(std::error::Error::source(&error).is_some());
        assert!(
            prepared.last_attempt_telemetry().is_none(),
            "{family:?} pre-cancel must not inherit previous telemetry"
        );
        assert!(
            prepared.last_failed_trajectory().is_none(),
            "{family:?} pre-cancel must not expose the previous successful trajectory"
        );
        assert_eq!(
            prepared
                .last_solution()
                .unwrap()
                .concentrations
                .last()
                .unwrap(),
            &accepted
        );
    }
}

#[test]
#[ignore = "wall-clock mid-run cancellation characterization"]
fn mid_run_cancellation_reports_native_diagnostic_availability() {
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let problem = NetworkProblem::new(
            vec![species("A", 1.0), species("B", 0.0)],
            vec![reaction(
                "r1",
                vec![("A", 1.0)],
                vec![("B", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.5),
                },
            )],
            TemperatureProgram::Constant(700.0),
            0.0,
            10.0,
        )
        .unwrap();
        let options = NetworkIvpOptions {
            family,
            max_step: 1e-5,
            telemetry: SolverTelemetryMode::Detailed,
            ..Default::default()
        };
        let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
        let cancellation = CancellationToken::new();
        let signal = cancellation.clone();
        let worker = std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(10));
            signal.cancel();
        });
        let result = prepared.solve_with_cancellation_token(&[1.0, 0.0], cancellation);
        worker.join().unwrap();
        let outcome = result.map(|_| "accepted").unwrap_or("typed failure");
        println!(
            "family={family:?} outcome={outcome} telemetry={} partial_points={} accepted={}",
            prepared.last_attempt_telemetry().is_some(),
            prepared
                .last_failed_trajectory()
                .map_or(0, |v| v.time_seconds.len()),
            prepared.last_solution().is_some(),
        );
        assert_eq!(prepared.last_solution().is_some(), outcome == "accepted");
    }
}

#[test]
fn prepared_isothermal_arrhenius_rebinds_temperature_and_initial_state() {
    let gas_constant = 8.314_462_618_153_24;
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![reaction(
            "arrhenius",
            vec![("A", 1.0)],
            vec![("B", 1.0)],
            RateLaw::MassAction {
                coefficient: RateCoefficient::Arrhenius {
                    prefactor: 1.0,
                    activation_energy: gas_constant * 500.0,
                },
            },
        )],
        TemperatureProgram::Constant(500.0),
        0.0,
        2.0,
    )
    .unwrap();
    let assembled = problem.assemble();
    assert_eq!(
        assembled.parameter_names,
        ["p_temperature", "p_rate_0_a", "p_rate_0_e"]
    );
    assert_eq!(assembled.parameter_values[0], 500.0);
    assert_eq!(assembled.parameter_values[1], 1.0);
    assert!((assembled.parameter_values[2] - gas_constant * 500.0).abs() < 1e-10);
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let options = NetworkIvpOptions {
            family,
            max_step: 0.05,
            ..Default::default()
        };
        let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
        let first = prepared.solve().unwrap();
        let expected_first = (-2.0 * (-1.0_f64).exp()).exp();
        assert!(
            (first.concentrations.last().unwrap()[0] - expected_first).abs() < 1e-2,
            "{family:?}"
        );
        let second = prepared
            .solve_with_temperature_and_rate_coefficients(
                &[2.0, 0.5],
                TemperatureBinding::Constant(1000.0),
                &[RateCoefficientBinding::new(
                    "arrhenius",
                    RateCoefficient::Arrhenius {
                        prefactor: 1.5,
                        activation_energy: gas_constant * 500.0,
                    },
                )],
            )
            .unwrap();
        assert_eq!(second.run_metadata.initial_concentrations, [2.0, 0.5]);
        assert_eq!(second.run_metadata.time_interval_seconds, (0.0, 2.0));
        assert_eq!(
            second.run_metadata.temperature_binding,
            Some(TemperatureBinding::Constant(1000.0))
        );
        assert_eq!(second.run_metadata.rate_coefficients.len(), 1);
        assert_eq!(
            second.run_metadata.rate_coefficients[0].reaction_id,
            "arrhenius"
        );
        assert!((second.concentrations[0][0] - 2.0).abs() < 1e-12);
        assert!((second.concentrations[0][1] - 0.5).abs() < 1e-12);
        let expected_second = 2.0 * (-2.0 * 1.5 * (-0.5_f64).exp()).exp();
        let end = second.concentrations.last().unwrap();
        if family == StiffSolverKind::BackwardEuler {
            let rate = 1.5 * (-0.5_f64).exp();
            let discrete_reference = second.time_seconds.windows(2).fold(2.0, |amount, pair| {
                amount / (1.0 + rate * (pair[1] - pair[0]))
            });
            assert!((end[0] - discrete_reference).abs() < 1e-8);
            assert!((end[0] - expected_second).abs() < 2e-2);
        } else {
            assert!((end[0] - expected_second).abs() < 1e-2, "{family:?}");
        }
        assert!(
            (end[1] - (2.5 - expected_second)).abs() < 2e-2,
            "{family:?}"
        );
        assert!(second.temperature_kelvin.iter().all(|&t| t == 1000.0));
        assert_eq!(
            prepared.current_temperature_binding(),
            Some(TemperatureBinding::Constant(1000.0))
        );

        let restored = prepared
            .solve_with_temperature_binding(&[1.0, 0.0], TemperatureBinding::Constant(500.0))
            .unwrap();
        assert_eq!(restored.run_metadata.initial_concentrations, [1.0, 0.0]);
        let expected_restored = (-2.0 * 1.5 * (-1.0_f64).exp()).exp();
        assert!(
            (restored.concentrations.last().unwrap()[0] - expected_restored).abs() < 1e-2,
            "{family:?}"
        );
    }
}

#[test]
fn prepared_rate_coefficients_rebind_by_reaction_id_across_families() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0), species("C", 0.0)],
        vec![
            reaction(
                "first",
                vec![("A", 1.0)],
                vec![("B", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.5),
                },
            ),
            reaction(
                "second",
                vec![("B", 1.0)],
                vec![("C", 1.0)],
                RateLaw::MassAction {
                    coefficient: RateCoefficient::Constant(0.2),
                },
            ),
        ],
        TemperatureProgram::Constant(700.0),
        0.0,
        1.0,
    )
    .unwrap();
    let first_rate = 1.2_f64;
    let second_rate = 0.4_f64;
    let expected_a = (-first_rate).exp();
    let expected_b =
        first_rate / (first_rate - second_rate) * ((-second_rate).exp() - (-first_rate).exp());
    let expected_c = 1.0 - expected_a - expected_b;

    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let options = NetworkIvpOptions {
            family,
            max_step: 0.002,
            ..Default::default()
        };
        let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
        prepared.solve().unwrap();
        let solution = prepared
            .solve_with_rate_coefficients(
                &[1.0, 0.0, 0.0],
                &[
                    RateCoefficientBinding::new("second", RateCoefficient::Constant(second_rate)),
                    RateCoefficientBinding::new("first", RateCoefficient::Constant(first_rate)),
                ],
            )
            .unwrap();
        let end = solution.concentrations.last().unwrap();
        assert!((end[0] - expected_a).abs() < 1e-2, "{family:?}: {end:?}");
        assert!((end[1] - expected_b).abs() < 1e-2, "{family:?}: {end:?}");
        assert!((end[2] - expected_c).abs() < 1e-2, "{family:?}: {end:?}");
    }
}

#[test]
fn invalid_rate_rebind_is_typed_and_preserves_last_solution() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![reaction(
            "r1",
            vec![("A", 1.0)],
            vec![("B", 1.0)],
            RateLaw::MassAction {
                coefficient: RateCoefficient::Constant(0.5),
            },
        )],
        TemperatureProgram::Constant(700.0),
        0.0,
        1.0,
    )
    .unwrap();
    let options = NetworkIvpOptions::default();
    let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
    let accepted = prepared
        .solve()
        .unwrap()
        .concentrations
        .last()
        .unwrap()
        .clone();

    let error = prepared
        .solve_with_rate_coefficients(&[1.0, 0.0], &[])
        .unwrap_err();
    assert!(matches!(
        error,
        NetworkIvpError::RateBindingCountMismatch {
            expected: 1,
            actual: 0
        }
    ));
    let error = prepared
        .solve_with_rate_coefficients(
            &[1.0, 0.0],
            &[
                RateCoefficientBinding::new("r1", RateCoefficient::Constant(0.8)),
                RateCoefficientBinding::new("r1", RateCoefficient::Constant(0.9)),
            ],
        )
        .unwrap_err();
    assert!(matches!(
        error,
        NetworkIvpError::DuplicateRateBinding { .. }
    ));
    let error = prepared
        .solve_with_rate_coefficients(
            &[1.0, 0.0],
            &[RateCoefficientBinding::new(
                "r1",
                RateCoefficient::Arrhenius {
                    prefactor: 1.0,
                    activation_energy: 1.0,
                },
            )],
        )
        .unwrap_err();
    assert!(matches!(
        error,
        NetworkIvpError::RateBindingKindMismatch { .. }
    ));
    let error = prepared
        .solve_with_rate_coefficients(
            &[1.0, 0.0],
            &[RateCoefficientBinding::new(
                "missing",
                RateCoefficient::Constant(0.8),
            )],
        )
        .unwrap_err();
    assert!(matches!(error, NetworkIvpError::UnknownRateBinding { .. }));
    let error = prepared
        .solve_with_rate_coefficients(
            &[1.0, 0.0],
            &[RateCoefficientBinding::new(
                "r1",
                RateCoefficient::Constant(f64::NAN),
            )],
        )
        .unwrap_err();
    assert!(matches!(error, NetworkIvpError::InvalidInput { .. }));
    assert_eq!(
        prepared
            .last_solution()
            .unwrap()
            .concentrations
            .last()
            .unwrap(),
        &accepted
    );
    assert_eq!(
        prepared.current_rate_coefficients()[0].coefficient,
        RateCoefficient::Constant(0.5)
    );
}

#[test]
fn invalid_published_samples_report_component_and_time() {
    let problem = a_to_two_b();
    let axis = DVector::from_vec(vec![0.0, 0.1]);
    let values = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, -1e-12, 1.0]);
    let error =
        super::solver::validate_candidate(&problem, None, &axis, &values, None).unwrap_err();
    assert!(matches!(
        error,
        NetworkIvpError::InvalidConcentration {
            species_id,
            time,
            value,
        } if species_id == "A" && time == 0.1 && value == -1e-12
    ));

    let bad_axis = DVector::from_vec(vec![0.0, 0.0]);
    let valid_values = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 1.0, 0.0]);
    let error = super::solver::validate_candidate(&problem, None, &bad_axis, &valid_values, None)
        .unwrap_err();
    assert!(matches!(
        error,
        NetworkIvpError::InvalidTimeSample {
            index: 1,
            time: 0.0,
            previous_time: Some(0.0),
        }
    ));
}

#[test]
fn prepared_linear_heating_rate_rebinds_without_rebuilding_equations() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![reaction(
            "heating",
            vec![("A", 1.0)],
            vec![("B", 1.0)],
            RateLaw::CustomExpr(
                Expr::Var("c0".into()) * Expr::Var("T".into()) / Expr::Const(500.0),
            ),
        )],
        TemperatureProgram::Linear {
            initial_kelvin: 500.0,
            kelvin_per_second: 100.0,
        },
        0.0,
        2.0,
    )
    .unwrap();
    let assembled = problem.assemble();
    assert_eq!(
        assembled.parameter_names,
        vec!["p_temperature", "p_heating_rate"]
    );
    assert_eq!(assembled.parameter_values, vec![500.0, 100.0]);
    for family in [
        StiffSolverKind::Bdf,
        StiffSolverKind::Radau,
        StiffSolverKind::BackwardEuler,
        StiffSolverKind::Lsode2,
    ] {
        let options = NetworkIvpOptions {
            family,
            max_step: 0.05,
            ..Default::default()
        };
        let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
        let first = prepared.solve().unwrap();
        assert!(
            (first.concentrations.last().unwrap()[0] - (-2.4_f64).exp()).abs() < 1e-2,
            "{family:?}"
        );
        let second = prepared
            .solve_with_temperature_binding(
                &[2.0, 0.5],
                TemperatureBinding::Linear {
                    initial_kelvin: 500.0,
                    kelvin_per_second: 200.0,
                },
            )
            .unwrap();
        let expected_a = 2.0 * (-2.8_f64).exp();
        let end = second.concentrations.last().unwrap();
        assert!((end[0] - expected_a).abs() < 1e-2, "{family:?}");
        assert!((end[1] - (2.5 - expected_a)).abs() < 2e-2, "{family:?}");
        assert!((*second.temperature_kelvin.last().unwrap() - 900.0).abs() < 1e-10);
    }
}

#[test]
fn invalid_temperature_rebind_is_typed_and_does_not_replace_last_solution() {
    let problem = a_to_two_b();
    let options = NetworkIvpOptions::default();
    let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
    let accepted = prepared
        .solve()
        .unwrap()
        .concentrations
        .last()
        .unwrap()
        .clone();
    let error = prepared
        .solve_with_temperature_binding(
            &[1.0, 0.0],
            TemperatureBinding::Linear {
                initial_kelvin: 500.0,
                kelvin_per_second: 1.0,
            },
        )
        .unwrap_err();
    assert!(matches!(
        error,
        NetworkIvpError::TemperatureBindingMismatch { .. }
    ));
    let error = prepared
        .solve_with_temperature_binding(&[1.0, 0.0], TemperatureBinding::Constant(-1.0))
        .unwrap_err();
    assert!(matches!(error, NetworkIvpError::TemperatureDomain { .. }));
    assert_eq!(
        prepared
            .last_solution()
            .unwrap()
            .concentrations
            .last()
            .unwrap(),
        &accepted
    );
    assert_eq!(
        prepared.current_temperature_binding(),
        Some(TemperatureBinding::Constant(700.0))
    );
    let restored = prepared.solve().unwrap();
    assert_eq!(restored.concentrations.last().unwrap(), &accepted);
}

#[test]
fn parameterized_temperature_rebind_has_atom_view_expr_legacy_parity() {
    let problem = NetworkProblem::new(
        vec![species("A", 1.0), species("B", 0.0)],
        vec![reaction(
            "temperature-driven",
            vec![("A", 1.0)],
            vec![("B", 1.0)],
            RateLaw::CustomExpr(
                Expr::Var("c0".into()) * Expr::Var("T".into()) / Expr::Const(500.0),
            ),
        )],
        TemperatureProgram::Linear {
            initial_kelvin: 500.0,
            kelvin_per_second: 100.0,
        },
        0.0,
        2.0,
    )
    .unwrap();
    let expr_options = NetworkIvpOptions {
        max_step: 0.05,
        ..Default::default()
    };
    let atom_options = NetworkIvpOptions {
        frontend: StiffFrontend::AtomViewNative,
        max_step: 0.05,
        ..Default::default()
    };
    let mut expr = PreparedNetworkIvp::new(&problem, &expr_options).unwrap();
    let mut atom = PreparedNetworkIvp::new(&problem, &atom_options).unwrap();
    for binding in [
        TemperatureBinding::Linear {
            initial_kelvin: 500.0,
            kelvin_per_second: 100.0,
        },
        TemperatureBinding::Linear {
            initial_kelvin: 600.0,
            kelvin_per_second: 200.0,
        },
    ] {
        let expr_solution = expr
            .solve_with_temperature_binding(&[1.0, 0.0], binding)
            .unwrap();
        let atom_solution = atom
            .solve_with_temperature_binding(&[1.0, 0.0], binding)
            .unwrap();
        let expr_final = expr_solution.concentrations.last().unwrap();
        let atom_final = atom_solution.concentrations.last().unwrap();
        for (left, right) in expr_final.iter().zip(atom_final) {
            assert!((left - right).abs() < 1e-8);
        }
        assert_eq!(
            expr_solution.temperature_kelvin.last(),
            atom_solution.temperature_kelvin.last()
        );
    }
}
