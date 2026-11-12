//! Story tests for the condensed-reactor IVP path.
//!
//! These tests exercise one canonical condensed-combustion fixture against
//! multiple solver-backend plans. The hypotheses are simple:
//! - the default typed plan should solve a small hydrogen dissociation case;
//! - fixed BDF and fixed Adams should stay usable on a non-stiff story;
//! - the alternate symbolic assembly route should not break a supported matrix
//!   configuration.
//!
//! Expected outcome:
//! - solve status is accepted;
//! - the axis is monotonic and the solution matrix stays finite;
//! - the default backend remains Lambdify + AtomView + Sparse;
//! - the H2 -> 2H conversion trend survives backend changes.

#[cfg(test)]
mod tests {
    use super::super::SimpleReactorIVP::*;
    use super::super::reactor_ivp_test_support::multistep_oxidation_coverage_task;
    use super::super::solver_backend::{
        ReactorIvpCallbackPolicy, ReactorIvpExecutionBackend, ReactorIvpMatrixBackend,
        ReactorIvpMethod, ReactorIvpSolverConfig, ReactorIvpSolverFamily,
        ReactorIvpSymbolicBackend,
    };
    use crate::Kinetics::mechfinder_api::ReactionData;
    use RustedSciThe::numerical::SolverTelemetryMode;
    use approx::assert_abs_diff_eq;
    use std::collections::HashMap;

    fn build_canonical_story_task() -> SimpleReactorTask {
        let mut task = SimpleReactorTask::new();
        task.kindata
            .set_reaction_data_directly(
                vec![ReactionData::new_elementary(
                    "H2=>2H".to_string(),
                    vec![2.5e2, 0.0, 8.0e3],
                    None,
                )],
                None,
            )
            .expect("reaction setup should succeed");
        task.kindata.stecheodata.vec_of_molmasses = Some(vec![2.016, 1.008]);
        task.kindata.substances = vec!["H2".to_string(), "H".to_string()];
        task.thermal_effects = vec![-1.5e5];
        task.ro = 1350.0;
        task.Cp = 1100.0;
        task.Lambda = 0.22;
        task.m = 0.012;
        task.scaling =
            crate::SimpleReactorsBVP::reactor_BVP_utils::ScalingConfig::new(80.0, 0.04, 90.0);
        task.initial_conditions = HashMap::from([
            ("T".to_string(), 520.0),
            ("q".to_string(), 0.05),
            ("H2".to_string(), 0.85),
            ("H".to_string(), 0.15),
        ]);
        task
    }

    fn assert_story_solution(
        task: &SimpleReactorTask,
        snapshot: &IvpSolveSnapshot,
        expected_config: ReactorIvpSolverConfig,
    ) {
        assert_eq!(snapshot.backend_config(), expected_config);
        assert_eq!(snapshot.resolved_backend_plan(), expected_config);
        assert_eq!(snapshot.family(), expected_config.family());
        if let Some(controller) = expected_config.lsode2_controller() {
            assert_eq!(snapshot.lsode2_controller(), Some(controller));
        }
        assert_eq!(snapshot.variable_count(), 4);
        if let Some(summary) = snapshot.summary.as_ref() {
            assert_eq!(summary.time_points, snapshot.time_points());
            assert_eq!(summary.variable_count, snapshot.variable_count());
        }
        assert!(snapshot.status().starts_with("finished"));
        assert!(snapshot.axis.iter().all(|value| value.is_finite()));
        assert!(snapshot.values.iter().all(|value| value.is_finite()));
        assert!(
            snapshot
                .axis
                .as_slice()
                .windows(2)
                .all(|pair| pair[1] >= pair[0])
        );
        assert!(!snapshot.method_label().is_empty());
        assert_eq!(snapshot.values.ncols(), task.canonical_state_names().len());

        let initial_h2 = task.initial_conditions["H2"];
        let initial_h = task.initial_conditions["H"];
        let final_row = snapshot.values.row(snapshot.values.nrows() - 1);

        assert!(
            (final_row[2] - initial_h2).abs() > 1e-9,
            "H2 should change in the canonical story fixture"
        );
        assert!(
            (final_row[3] - initial_h).abs() > 1e-9,
            "H should change in the canonical story fixture"
        );

        assert_abs_diff_eq!(snapshot.axis[0], expected_config.x0, epsilon = 1e-12);
        assert_abs_diff_eq!(
            snapshot.axis[snapshot.axis.len() - 1],
            expected_config.x_bound,
            epsilon = 1e-10
        );
    }

    fn solve_story_with_config(
        config: ReactorIvpSolverConfig,
    ) -> (SimpleReactorTask, IvpSolveSnapshot) {
        let mut task = build_canonical_story_task();
        task.set_solver_backend_config(config);
        task.setup_condensed_ivp()
            .expect("canonical story setup should succeed");
        let snapshot = task.solve().expect("canonical story solve should succeed");
        (task, snapshot)
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_default_lambdify_sparse() {
        let config = ReactorIvpSolverConfig::default();
        let (task, snapshot) = solve_story_with_config(config);

        assert_eq!(
            task.solver_backend_config.execution_backend,
            ReactorIvpExecutionBackend::Lambdify
        );
        assert_eq!(
            task.solver_backend_config.symbolic_backend,
            ReactorIvpSymbolicBackend::AtomView
        );
        assert_eq!(
            task.solver_backend_config.matrix_backend,
            ReactorIvpMatrixBackend::Sparse
        );
        assert_eq!(
            task.solver_backend_config.lsode2_controller(),
            Some(ReactorIvpMethod::Auto)
        );
        assert_story_solution(&task, &snapshot, config);
        assert!(snapshot.stage_telemetry.stages.is_empty());
    }

    #[test]
    fn multistep_oxidation_coverage_fixture_assembles_balanced_network() {
        let mut task = multistep_oxidation_coverage_task();
        task.setup_condensed_ivp()
            .expect("multistep coverage fixture should assemble");

        assert_eq!(task.solver.eq_system.len(), 7);
        assert_eq!(task.solver.unknowns.len(), 7);
        assert_eq!(task.map_eq_rate.len(), 2);
        assert_eq!(task.kindata.vec_of_equations.len(), 2);
    }

    #[test]
    fn callback_execution_policies_preserve_largest_solved_ivp_fixture_and_report_dispatches() {
        let policies = [
            ("sequential", ReactorIvpCallbackPolicy::Sequential),
            ("auto", ReactorIvpCallbackPolicy::Auto { min_work: 0 }),
            (
                "forced_parallel",
                ReactorIvpCallbackPolicy::Parallel { min_work: 0 },
            ),
        ];
        for family in [
            ReactorIvpSolverFamily::Lsode2,
            ReactorIvpSolverFamily::Bdf,
            ReactorIvpSolverFamily::Radau,
        ] {
            let mut reference: Option<Vec<f64>> = None;
            for (label, policy) in policies {
                let mut task = build_canonical_story_task();
                task.kindata
                    .set_reaction_data_directly(
                        vec![ReactionData::new_elementary(
                            "H2=>2H".to_string(),
                            vec![1.0e-3, 0.0, 8.0e3],
                            None,
                        )],
                        None,
                    )
                    .expect("bounded policy fixture reaction should parse");
                task.thermal_effects = vec![0.0];
                task.initial_conditions.insert("T".to_string(), 900.0);
                task.initial_conditions.insert("q".to_string(), 0.0);
                task.initial_conditions.insert("H2".to_string(), 0.70);
                task.initial_conditions.insert("H".to_string(), 0.30);
                let config = ReactorIvpSolverConfig::default()
                    .with_family(family)
                    .with_symbolic_backend(ReactorIvpSymbolicBackend::ExprLegacy)
                    .with_matrix_backend(ReactorIvpMatrixBackend::Dense)
                    .with_integration_domain(0.0, 0.02)
                    .with_max_step(1.0e-3)
                    .with_rtol(1.0e-6)
                    .with_atol(1.0e-8)
                    .with_telemetry(SolverTelemetryMode::Detailed)
                    .with_callback_policy(policy);
                task.set_solver_backend_config(config);
                let snapshot = task.solve().unwrap_or_else(|error| {
                    panic!("largest IVP fixture {family:?}/{label} failed: {error}")
                });
                assert!(snapshot.status().starts_with("finished"));
                assert!(snapshot.time_points() > 1);
                assert!(snapshot.values.iter().all(|value| value.is_finite()));
                let final_state = snapshot
                    .values
                    .row(snapshot.values.nrows() - 1)
                    .iter()
                    .copied()
                    .collect::<Vec<_>>();
                if let Some(reference) = &reference {
                    assert_eq!(reference.len(), final_state.len());
                    for (index, (expected, actual)) in
                        reference.iter().zip(&final_state).enumerate()
                    {
                        assert!(
                            (expected - actual).abs() <= 2e-6 * (1.0 + expected.abs()),
                            "{family:?}/{label} final state[{index}] differs: {expected:e} vs {actual:e}"
                        );
                    }
                } else {
                    reference = Some(final_state);
                }
                let rows = snapshot.diagnostics_report_rows();
                let configured_policy = rows
                    .iter()
                    .find(|row| row.field == "callback_execution_policy")
                    .map(|row| row.value.as_str());
                if let Some(configured_policy) = configured_policy {
                    assert_eq!(
                        configured_policy,
                        format!("{policy:?}"),
                        "{family:?}/{label} did not preserve the requested callback policy"
                    );
                } else {
                    assert_eq!(family, ReactorIvpSolverFamily::Radau);
                }
                let read_counter = |field: &str| {
                    rows.iter()
                        .find(|row| row.field == field)
                        .and_then(|row| row.value.parse::<u64>().ok())
                        .unwrap_or(0)
                };
                let parallel = read_counter("native_counter::parallel_dispatches");
                let sequential = read_counter("native_counter::sequential_dispatches");
                let residual_callbacks = read_counter("native_counter::residual_callbacks");
                let jacobian_callbacks = read_counter("native_counter::jacobian_callbacks");
                let native_residual_calls = read_counter("native_counter::residual_calls");
                let native_jacobian_calls = read_counter("native_counter::jacobian_calls");
                assert!(
                    parallel + sequential > 0,
                    "{family:?}/{label} did not report callback dispatches"
                );
                if family != ReactorIvpSolverFamily::Radau {
                    assert!(
                        residual_callbacks > 0 && jacobian_callbacks > 0,
                        "{family:?}/{label} omitted residual/Jacobian callback counts: residual={residual_callbacks}, jacobian={jacobian_callbacks}"
                    );
                }
                if family == ReactorIvpSolverFamily::Bdf {
                    assert!(
                        native_residual_calls > 0 && native_jacobian_calls > 0,
                        "BDF native solver counters are missing: residual={native_residual_calls}, jacobian={native_jacobian_calls}"
                    );
                }
                if label == "sequential" {
                    assert_eq!(
                        parallel, 0,
                        "{family:?} sequential policy dispatched parallel work"
                    );
                }
                if label == "forced_parallel" {
                    assert!(
                        parallel > 0,
                        "{family:?} forced-parallel policy was not honored"
                    );
                }
                eprintln!(
                    "IVP callback policy family={family:?} policy={label} parallel={parallel} sequential={sequential} residual_callbacks={residual_callbacks} jacobian_callbacks={jacobian_callbacks} native_residual_calls={native_residual_calls} native_jacobian_calls={native_jacobian_calls} telemetry={} ",
                    if parallel + sequential == 0 {
                        "unavailable"
                    } else {
                        "reported"
                    }
                );
            }
        }
    }

    #[test]
    fn symbolic_frontends_match_endpoint_for_each_solver_family() {
        for (scenario, temperature, initial_h2, initial_h) in [
            ("baseline-520K", 520.0, 0.85, 0.15),
            ("warm-inlet-900K", 900.0, 0.70, 0.30),
        ] {
            for family in [
                ReactorIvpSolverFamily::Lsode2,
                ReactorIvpSolverFamily::Bdf,
                ReactorIvpSolverFamily::Radau,
            ] {
                let solve_frontend = |frontend| {
                    let config = ReactorIvpSolverConfig::default()
                        .with_family(family)
                        .with_symbolic_backend(frontend)
                        .with_matrix_backend(ReactorIvpMatrixBackend::Dense)
                        .with_integration_domain(0.0, 0.02)
                        .with_max_step(1.0e-3)
                        .with_rtol(1.0e-6)
                        .with_atol(1.0e-8);
                    let mut task = build_canonical_story_task();
                    task.initial_conditions.insert("T".to_string(), temperature);
                    task.initial_conditions.insert("H2".to_string(), initial_h2);
                    task.initial_conditions.insert("H".to_string(), initial_h);
                    task.kindata
                        .set_reaction_data_directly(
                            vec![ReactionData::new_elementary(
                                "H2=>2H".to_string(),
                                vec![1.0e-3, 0.0, 8.0e3],
                                None,
                            )],
                            None,
                        )
                        .expect("low-rate benchmark reaction should parse");
                    task.thermal_effects = vec![0.0];
                    task.initial_conditions.insert("q".to_string(), 0.0);
                    task.set_solver_backend_config(config);
                    let snapshot = task.solve().unwrap_or_else(|error| {
                        panic!("{scenario}/{family:?}/{frontend:?} failed: {error}")
                    });
                    assert!(snapshot.status().starts_with("finished"));
                    let physical_temperatures = snapshot
                        .values
                        .column(0)
                        .map(|temperature| temperature * task.scaling.T_scale + task.scaling.dT);
                    assert!(physical_temperatures.iter().all(|temperature| {
                        temperature.is_finite() && (0.0..=3000.0).contains(temperature)
                    }));
                    assert!(
                        snapshot
                            .values
                            .columns(2, 2)
                            .iter()
                            .all(|value| (0.0..=1.0).contains(value))
                    );
                    snapshot
                        .values
                        .row(snapshot.values.nrows() - 1)
                        .iter()
                        .copied()
                        .collect::<Vec<_>>()
                };

                let expr = solve_frontend(ReactorIvpSymbolicBackend::ExprLegacy);
                let atom = solve_frontend(ReactorIvpSymbolicBackend::AtomView);
                let max_abs_endpoint_delta = expr
                    .iter()
                    .zip(&atom)
                    .map(|(left, right)| (left - right).abs())
                    .fold(0.0_f64, f64::max);
                println!(
                    "{scenario}/{family:?} ExprLegacy/AtomView max_abs_endpoint_delta={max_abs_endpoint_delta:.6e}"
                );

                assert!(
                    expr.iter()
                        .all(|value| value.is_finite() && value.abs() < 10.0)
                );
                assert!(
                    atom.iter()
                        .all(|value| value.is_finite() && value.abs() < 10.0)
                );
                assert!(
                    (expr[2] - initial_h2).abs() > 1.0e-7,
                    "fixture should evolve H2 measurably"
                );
                assert!(
                    (expr[3] - initial_h).abs() > 1.0e-7,
                    "fixture should evolve H measurably"
                );
                assert!(
                    expr[2] > 0.0 && expr[3] > 0.0,
                    "species values must stay nonnegative"
                );
                assert!(
                    (expr[2] + expr[3] - 1.0).abs() < 1.0e-12,
                    "mass fractions should conserve"
                );

                for (index, (left, right)) in expr.iter().zip(&atom).enumerate() {
                    let tolerance = 1.0e-8 * (1.0 + left.abs());
                    assert!(
                        (left - right).abs() <= tolerance,
                        "{scenario}/{family:?} endpoint state[{index}] differs: ExprLegacy={left:.9e}, AtomView={right:.9e}, tolerance={tolerance:.3e}"
                    );
                }
            }
        }
    }

    #[test]
    fn physical_smoke_scenarios_keep_temperature_and_mass_fractions_bounded() {
        let mut scenarios = vec![("h2_dissociation", build_canonical_story_task())];
        let mut oxidation = multistep_oxidation_coverage_task();
        oxidation
            .kindata
            .set_reaction_data_directly(
                vec![
                    ReactionData::new_elementary(
                        "2CH4+3O2=>2CO+4H2O".to_string(),
                        vec![1.0e-25, 0.0, 8.0e3],
                        None,
                    ),
                    ReactionData::new_elementary(
                        "2CO+O2=>2CO2".to_string(),
                        vec![1.0e-25, 0.0, 5.0e3],
                        None,
                    ),
                ],
                None,
            )
            .expect("low-rate oxidation reactions should parse");
        oxidation.thermal_effects = vec![-2.0e5, -2.83e5];
        oxidation.initial_conditions.insert("CH4".into(), 0.2);
        oxidation.initial_conditions.insert("O2".into(), 0.4);
        oxidation.initial_conditions.insert("CO".into(), 0.1);
        oxidation.initial_conditions.insert("H2O".into(), 0.1);
        oxidation.initial_conditions.insert("CO2".into(), 0.2);
        oxidation.set_solver_backend_config(
            ReactorIvpSolverConfig::default()
                .with_family(ReactorIvpSolverFamily::Lsode2)
                .with_symbolic_backend(ReactorIvpSymbolicBackend::ExprLegacy)
                .with_matrix_backend(ReactorIvpMatrixBackend::Dense)
                .with_integration_domain(0.0, 1.0e-4)
                .with_max_step(1.0e-5)
                .with_rtol(1.0e-6)
                .with_atol(1.0e-9),
        );
        scenarios.push(("two_step_oxidation", oxidation));

        for (name, mut task) in scenarios {
            if name == "h2_dissociation" {
                task.kindata
                    .set_reaction_data_directly(
                        vec![ReactionData::new_elementary(
                            "H2=>2H".to_string(),
                            vec![1.0e-3, 0.0, 8.0e3],
                            None,
                        )],
                        None,
                    )
                    .expect("bounded reaction should parse");
                task.thermal_effects = vec![0.0];
                task.initial_conditions.insert("q".into(), 0.0);
                task.set_solver_backend_config(
                    ReactorIvpSolverConfig::default()
                        .with_family(ReactorIvpSolverFamily::Lsode2)
                        .with_symbolic_backend(ReactorIvpSymbolicBackend::ExprLegacy)
                        .with_matrix_backend(ReactorIvpMatrixBackend::Dense)
                        .with_integration_domain(0.0, 0.02)
                        .with_max_step(1.0e-3)
                        .with_rtol(1.0e-6)
                        .with_atol(1.0e-8),
                );
            }
            let snapshot = task
                .solve()
                .unwrap_or_else(|error| panic!("{name} physical smoke solve failed: {error}"));
            if name == "two_step_oxidation" {
                assert_oxidation_element_balances(&snapshot);
            }
            let temperature_index = snapshot
                .variable_names
                .iter()
                .position(|variable| variable == "Teta")
                .expect("canonical temperature state should exist");
            let physical_temperatures = snapshot
                .values
                .column(temperature_index)
                .map(|temperature| temperature * task.scaling.T_scale + task.scaling.dT);
            assert!(
                physical_temperatures
                    .iter()
                    .all(|temperature| temperature.is_finite()
                        && (0.0..=3000.0).contains(temperature)),
                "{name} produced temperature outside [0, 3000] K: {physical_temperatures:?}"
            );
            for (column, variable) in snapshot.variable_names.iter().enumerate() {
                if variable == "Teta" || variable == "q" {
                    continue;
                }
                let fractions = snapshot.values.column(column);
                assert!(
                    fractions
                        .iter()
                        .all(|fraction| fraction.is_finite() && (0.0..=1.0).contains(fraction)),
                    "{name} produced mass fraction outside [0, 1] for {variable}: min={}, max={}",
                    fractions.iter().copied().fold(f64::INFINITY, f64::min),
                    fractions.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                );
            }
            println!(
                "physical smoke {name}: T={:.3}..{:.3} K, states={}, points={}",
                physical_temperatures
                    .iter()
                    .copied()
                    .fold(f64::INFINITY, f64::min),
                physical_temperatures
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max),
                snapshot.variable_count(),
                snapshot.time_points()
            );
        }
    }

    fn assert_oxidation_element_balances(snapshot: &IvpSolveSnapshot) {
        let species = [
            ("CH4", 16.043, [1.0, 4.0, 0.0]),
            ("O2", 31.998, [0.0, 0.0, 2.0]),
            ("CO", 28.010, [1.0, 0.0, 1.0]),
            ("H2O", 18.015, [0.0, 2.0, 1.0]),
            ("CO2", 44.009, [1.0, 0.0, 2.0]),
        ];
        let columns = species
            .iter()
            .enumerate()
            .map(|(index, (name, _, _))| {
                snapshot
                    .variable_names
                    .iter()
                    .position(|variable| variable == &format!("C{index}"))
                    .unwrap_or_else(|| panic!("missing mass-fraction state for {name}"))
            })
            .collect::<Vec<_>>();
        let elemental_mass_contributions = |row: usize| -> [f64; 3] {
            std::array::from_fn(|element| {
                species
                    .iter()
                    .zip(&columns)
                    .map(|((_, molar_mass, counts), column)| {
                        snapshot.values[(row, *column)]
                            * counts[element]
                            * [12.011, 1.008, 15.999][element]
                            / molar_mass
                    })
                    .sum::<f64>()
            })
        };
        let initial = elemental_mass_contributions(0);
        for row in 0..snapshot.values.nrows() {
            let actual = elemental_mass_contributions(row);
            for (element, (expected, actual)) in initial.iter().zip(actual).enumerate() {
                let tolerance = 2.0e-7_f64.max(expected.abs() * 2.0e-7);
                assert!(
                    (expected - actual).abs() <= tolerance,
                    "oxidation row {row} elemental mass fraction[{element}] drifted: {expected:e} vs {actual:e} (tol {tolerance:e})"
                );
            }
        }
    }

    #[test]
    fn multistep_oxidation_preserves_elemental_inventory_across_solver_families_and_frontends() {
        let mut reference_task = multistep_oxidation_coverage_task();
        reference_task
            .kindata
            .set_reaction_data_directly(
                vec![
                    ReactionData::new_elementary(
                        "2CH4+3O2=>2CO+4H2O".to_string(),
                        vec![1.0e-25, 0.0, 8.0e3],
                        None,
                    ),
                    ReactionData::new_elementary(
                        "2CO+O2=>2CO2".to_string(),
                        vec![1.0e-25, 0.0, 5.0e3],
                        None,
                    ),
                ],
                None,
            )
            .expect("bounded oxidation reactions should parse");
        reference_task.thermal_effects = vec![-2.0e5, -2.83e5];
        reference_task.initial_conditions = HashMap::from([
            ("T".into(), 520.0),
            ("q".into(), 0.0),
            ("CH4".into(), 0.2),
            ("O2".into(), 0.4),
            ("CO".into(), 0.1),
            ("H2O".into(), 0.1),
            ("CO2".into(), 0.2),
        ]);

        for family in [
            ReactorIvpSolverFamily::Lsode2,
            ReactorIvpSolverFamily::Bdf,
            ReactorIvpSolverFamily::Radau,
        ] {
            for frontend in [
                ReactorIvpSymbolicBackend::ExprLegacy,
                ReactorIvpSymbolicBackend::AtomView,
            ] {
                let mut task = reference_task.clone();
                task.set_solver_backend_config(
                    ReactorIvpSolverConfig::default()
                        .with_family(family)
                        .with_symbolic_backend(frontend)
                        .with_matrix_backend(ReactorIvpMatrixBackend::Dense)
                        .with_integration_domain(0.0, 1.0e-4)
                        .with_max_step(1.0e-5)
                        .with_rtol(1.0e-6)
                        .with_atol(1.0e-9),
                );
                let snapshot = task.solve().unwrap_or_else(|error| {
                    panic!("multistep oxidation {family:?}/{frontend:?} failed: {error}")
                });
                assert!(snapshot.status().starts_with("finished"));
                assert_oxidation_element_balances(&snapshot);

                let temperature_index = snapshot
                    .variable_names
                    .iter()
                    .position(|variable| variable == "Teta")
                    .expect("temperature state should be present");
                for value in snapshot.values.column(temperature_index).iter() {
                    let physical_temperature = value * task.scaling.T_scale + task.scaling.dT;
                    assert!((0.0..=3000.0).contains(&physical_temperature));
                }
                for (column, variable) in snapshot.variable_names.iter().enumerate() {
                    if variable == "Teta" || variable == "q" {
                        continue;
                    }
                    assert!(snapshot.values.column(column).iter().all(|fraction| {
                        fraction.is_finite() && (0.0..=1.0).contains(fraction)
                    }));
                }
                println!(
                    "multistep oxidation family={family:?} frontend={frontend:?} points={} T0={:.6} K Tf={:.6} K",
                    snapshot.time_points(),
                    snapshot.values[(0, temperature_index)] * task.scaling.T_scale
                        + task.scaling.dT,
                    snapshot.values[(snapshot.values.nrows() - 1, temperature_index)]
                        * task.scaling.T_scale
                        + task.scaling.dT,
                );
            }
        }
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_fixed_bdf_sparse() {
        let config = ReactorIvpSolverConfig::default()
            .with_lsode2_controller(ReactorIvpMethod::Bdf)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let (task, snapshot) = solve_story_with_config(config);

        assert_eq!(
            task.solver_backend_config.lsode2_controller(),
            Some(ReactorIvpMethod::Bdf)
        );
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_fixed_adams_sparse() {
        let config = ReactorIvpSolverConfig::default()
            .with_lsode2_controller(ReactorIvpMethod::Adams)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let (task, snapshot) = solve_story_with_config(config);

        assert_eq!(
            task.solver_backend_config.lsode2_controller(),
            Some(ReactorIvpMethod::Adams)
        );
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_standalone_bdf() {
        let config = ReactorIvpSolverConfig::default()
            .with_family(ReactorIvpSolverFamily::Bdf)
            .with_telemetry(SolverTelemetryMode::Detailed);
        let (task, snapshot) = solve_story_with_config(config);

        assert_story_solution(&task, &snapshot, config);
        assert_eq!(snapshot.family(), ReactorIvpSolverFamily::Bdf);
        assert_eq!(snapshot.method_label(), "BDF");
        assert!(snapshot.summary.is_none());
        assert!(snapshot.native_telemetry.is_some());
        let rows = snapshot.diagnostics_report_rows();
        assert!(
            rows.iter()
                .any(|row| row.field == "native_counter::residual_calls")
        );
        assert!(
            rows.iter()
                .any(|row| row.field == "native_timing_ms::solve_ms")
        );
        for stage in [
            "contract_validation",
            "scaling_processing",
            "kinetic_processing",
            "peclet_numbers",
            "equation_assembly",
            "pre_solve_validation",
            "solver_request_build",
            "rst_solver_preparation",
            "native_integration",
            "trajectory_materialization",
            "candidate_validation",
            "solution_publication",
            "conservation",
            "postprocessing",
            "total_workflow",
        ] {
            assert!(
                rows.iter()
                    .any(|row| row.field == format!("ivp_stage::{stage}::calls")),
                "missing stage counter for {stage}"
            );
            assert!(
                rows.iter()
                    .any(|row| row.field == format!("ivp_stage::{stage}::elapsed_ms")),
                "missing stage timing for {stage}"
            );
        }
    }

    #[test]
    fn stage_telemetry_modes_keep_counter_and_clock_costs_separate() {
        let mut off = IvpStageTelemetry::new(SolverTelemetryMode::Off);
        off.run("fixture", || Ok::<_, ()>(())).unwrap();
        assert!(off.stages.is_empty());

        let mut counters = IvpStageTelemetry::new(SolverTelemetryMode::Counters);
        counters.run("fixture", || Ok::<_, ()>(())).unwrap();
        assert_eq!(counters.stages[0].calls, Some(1));
        assert_eq!(counters.stages[0].elapsed_ms, None);

        let mut timings = IvpStageTelemetry::new(SolverTelemetryMode::Timings);
        timings.run("fixture", || Ok::<_, ()>(())).unwrap();
        assert_eq!(timings.stages[0].calls, None);
        assert!(timings.stages[0].elapsed_ms.is_some());
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_radau() {
        let config = ReactorIvpSolverConfig::default()
            .with_family(ReactorIvpSolverFamily::Radau)
            .with_telemetry(SolverTelemetryMode::Detailed);
        let (task, snapshot) = solve_story_with_config(config);

        assert_story_solution(&task, &snapshot, config);
        assert_eq!(snapshot.family(), ReactorIvpSolverFamily::Radau);
        assert_eq!(snapshot.method_label(), "Radau5");
        assert!(snapshot.summary.is_none());
        assert!(snapshot.native_telemetry.is_some());
        let rows = snapshot.diagnostics_report_rows();
        assert!(
            rows.iter()
                .any(|row| row.field == "native_counter::residual_calls")
        );
        assert!(
            rows.iter()
                .any(|row| row.field == "native_timing_ms::callback_ms")
        );
    }

    #[test]
    fn hybrid_bdf_document_settings_prepare_the_caller_built_reactor_request() {
        let document = r#"
task
solver: IVP
method: BDF

process_conditions
ro: 1350.0
Cp: 1100.0
Lambda: 0.22
m: 0.012
L: 0.04
thermal_effects: [-150000.0]

scaling
dT: 80.0
T_scale: 90.0
L: 0.04

integration_domain
x0: 0.0
x_bound: 0.02

initial_conditions
T: 520.0
q: 0.05
H2: 0.85
H: 0.15

solver_options
rtol: 1e-6
atol: 1e-8
max_step: 1e-3
max_iterations: 20000
telemetry: detailed
logging: off
"#;
        let spec =
            super::super::task_parser_reactor_IVP::parse_reactor_ivp_document_from_str(document)
                .expect("hybrid BDF settings should parse");
        assert_eq!(
            spec.hybrid_solver_settings
                .as_ref()
                .expect("RST hybrid settings should be retained")
                .solver_kind(),
            RustedSciThe::numerical::StiffAPI::StiffSolverKind::Bdf
        );

        let mut task = build_canonical_story_task();
        spec.apply_to_task(&mut task)
            .expect("KiThe physical document should apply");
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default()
                .with_family(ReactorIvpSolverFamily::Bdf)
                .with_integration_domain(spec.x0, spec.x_bound)
                .with_telemetry(SolverTelemetryMode::Detailed),
        );
        task.setup_condensed_ivp()
            .expect("caller-built reactor equations should assemble");
        let snapshot = task
            .solve()
            .expect("hybrid settings should prepare and solve the KiThe request");
        assert_eq!(snapshot.family, ReactorIvpSolverFamily::Bdf);
        assert!(snapshot.native_telemetry.is_some());
    }

    #[test]
    fn radau_matrix_layouts_match_for_both_symbolic_frontends() {
        let mut reference: Option<Vec<f64>> = None;
        for frontend in [
            ReactorIvpSymbolicBackend::ExprLegacy,
            ReactorIvpSymbolicBackend::AtomView,
        ] {
            for layout in [
                ReactorIvpMatrixBackend::Dense,
                ReactorIvpMatrixBackend::Sparse,
                ReactorIvpMatrixBackend::Banded,
            ] {
                let config = ReactorIvpSolverConfig::default()
                    .with_family(ReactorIvpSolverFamily::Radau)
                    .with_symbolic_backend(frontend)
                    .with_matrix_backend(layout)
                    .with_integration_domain(0.0, 0.02)
                    .with_max_step(1.0e-3)
                    .with_rtol(1.0e-6)
                    .with_atol(1.0e-8);
                let mut task = build_canonical_story_task();
                task.kindata
                    .set_reaction_data_directly(
                        vec![ReactionData::new_elementary(
                            "H2=>2H".to_string(),
                            vec![1.0e-3, 0.0, 8.0e3],
                            None,
                        )],
                        None,
                    )
                    .expect("low-rate Radau fixture should parse");
                task.thermal_effects = vec![0.0];
                task.initial_conditions.insert("q".to_string(), 0.0);
                task.set_solver_backend_config(config);
                let snapshot = task.solve().unwrap_or_else(|error| {
                    panic!("Radau/{frontend:?}/{layout:?} failed: {error}")
                });
                assert!(snapshot.status().starts_with("finished"));
                let physical_temperatures = snapshot
                    .values
                    .column(0)
                    .map(|temperature| temperature * task.scaling.T_scale + task.scaling.dT);
                assert!(physical_temperatures.iter().all(|temperature| {
                    temperature.is_finite() && (0.0..=3000.0).contains(temperature)
                }));
                assert!(
                    snapshot
                        .values
                        .columns(2, 2)
                        .iter()
                        .all(|value| (0.0..=1.0).contains(value))
                );
                let final_state = snapshot.values.row(snapshot.values.nrows() - 1);
                assert!((final_state[2] - 0.85).abs() > 1.0e-7);
                assert!((final_state[3] - 0.15).abs() > 1.0e-7);
                assert!((final_state[2] + final_state[3] - 1.0).abs() <= 1.0e-10);
                let profile = snapshot.values.iter().copied().collect::<Vec<_>>();
                if let Some(expected) = &reference {
                    assert_eq!(profile.len(), expected.len());
                    let max_delta = profile
                        .iter()
                        .zip(expected)
                        .map(|(actual, expected)| (actual - expected).abs())
                        .fold(0.0_f64, f64::max);
                    println!("Radau/{frontend:?}/{layout:?} max_abs_profile_delta={max_delta:.6e}");
                    assert!(max_delta <= 2.0e-6, "profile drift {max_delta:e}");
                } else {
                    reference = Some(profile);
                }
            }
        }
    }

    #[test]
    fn failed_standalone_bdf_solve_retains_partial_native_telemetry() {
        let config = ReactorIvpSolverConfig::default()
            .with_family(ReactorIvpSolverFamily::Bdf)
            .with_bdf_max_steps(1)
            .with_telemetry(SolverTelemetryMode::Detailed);
        let mut task = build_canonical_story_task();
        task.set_solver_backend_config(config);

        let error = task
            .solve_condensed_ivp()
            .expect_err("one accepted step should not reach the requested terminal coordinate");
        match error {
            IvpError::SolverExecution {
                error,
                telemetry:
                    Some(RustedSciThe::numerical::StiffAPI::StiffTelemetrySnapshot::Bdf {
                        mode,
                        statistics,
                        ..
                    }),
            } => {
                assert_eq!(
                    error.solver_kind(),
                    RustedSciThe::numerical::StiffAPI::StiffSolverKind::Bdf
                );
                assert_eq!(mode, SolverTelemetryMode::Detailed);
                assert!(statistics.residual_calls > 0);
                assert!(task.last_solve_snapshot.is_none());
            }
            other => panic!("expected typed BDF failure with partial telemetry, got {other:?}"),
        }
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_exprlegacy_and_banded_matrix() {
        let config = ReactorIvpSolverConfig::default()
            .with_lsode2_controller(ReactorIvpMethod::Auto)
            .with_symbolic_backend(ReactorIvpSymbolicBackend::ExprLegacy)
            .with_matrix_backend(ReactorIvpMatrixBackend::Banded)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let (task, snapshot) = solve_story_with_config(config);

        assert_eq!(
            task.solver_backend_config.symbolic_backend,
            ReactorIvpSymbolicBackend::ExprLegacy
        );
        assert!(task.solver_backend_config.matrix_backend.is_banded());
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_inferred_banded_matrix() {
        let config = ReactorIvpSolverConfig::default()
            .with_lsode2_controller(ReactorIvpMethod::Auto)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_symbolic_backend(ReactorIvpSymbolicBackend::AtomView)
            .with_matrix_backend(ReactorIvpMatrixBackend::Banded)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let (task, snapshot) = solve_story_with_config(config);

        assert!(task.solver_backend_config.matrix_backend.is_banded());
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_sparse_and_banded_routes_match_on_the_canonical_story_fixture() {
        let sparse_config = ReactorIvpSolverConfig::default()
            .with_lsode2_controller(ReactorIvpMethod::Auto)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_symbolic_backend(ReactorIvpSymbolicBackend::AtomView)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let banded_config = ReactorIvpSolverConfig::default()
            .with_lsode2_controller(ReactorIvpMethod::Auto)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_symbolic_backend(ReactorIvpSymbolicBackend::AtomView)
            .with_matrix_backend(ReactorIvpMatrixBackend::Banded)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);

        let (sparse_task, sparse_snapshot) = solve_story_with_config(sparse_config);
        let (banded_task, banded_snapshot) = solve_story_with_config(banded_config);

        assert_eq!(
            sparse_task.solver_backend_config.matrix_backend,
            ReactorIvpMatrixBackend::Sparse
        );
        assert!(banded_task.solver_backend_config.matrix_backend.is_banded());
        assert_eq!(sparse_snapshot.axis.len(), banded_snapshot.axis.len());
        assert_eq!(
            sparse_snapshot.values.ncols(),
            banded_snapshot.values.ncols()
        );
        for (lhs, rhs) in sparse_snapshot.axis.iter().zip(banded_snapshot.axis.iter()) {
            assert_abs_diff_eq!(lhs, rhs, epsilon = 1e-10);
        }
        for (lhs, rhs) in sparse_snapshot
            .values
            .row(sparse_snapshot.values.nrows() - 1)
            .iter()
            .zip(
                banded_snapshot
                    .values
                    .row(banded_snapshot.values.nrows() - 1)
                    .iter(),
            )
        {
            let scale = lhs.abs().max(rhs.abs()).max(1.0);
            let rel = (lhs - rhs).abs() / scale;
            assert!(
                rel <= 1.0e-8,
                "Sparse and banded routes diverged too much on the canonical story fixture: lhs={lhs:e}, rhs={rhs:e}, rel={rel:e}"
            );
        }
    }

    #[test]
    fn test_canonical_story_fixture_reuses_initial_state_shape_after_rerun() {
        let mut task = build_canonical_story_task();
        task.setup_condensed_ivp()
            .expect("canonical story setup should succeed");
        let first = task.solve().expect("first solve should succeed");

        let override_config = ReactorIvpSolverConfig::default()
            .with_lsode2_controller(ReactorIvpMethod::Bdf)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 1.0);
        let second = task
            .solve_with_config(override_config)
            .expect("second solve should succeed");

        assert_eq!(first.variable_count(), second.variable_count());
        assert!(second.time_points() > 0);
        assert_eq!(
            task.last_solve_snapshot.as_ref().unwrap().status(),
            second.status()
        );
    }
}
