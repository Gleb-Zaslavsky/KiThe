#[cfg(test)]
mod tests {
    use super::super::SimpleReactorBVP::NrbvpHandoffConfig;
    use super::super::SimpleReactorBVP::ReactorError;
    use super::super::reactor_bvp_test_support::{
        compact_hmx_reactor, compact_hmx_reactor_with_product_fraction,
    };
    use super::super::solver_backend::ReactorBvpSolverConfig;
    use crate::Kinetics::User_reactions::KinDataState;
    use crate::SimpleReactorsBVP::reactor_BVP_utils::{BoundsConfig, ToleranceConfig};
    use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_frozen::FrozenSolverOptions;
    use RustedSciThe::numerical::BVP_Damp::telemetry::BvpTelemetryMode;
    use RustedSciThe::symbolic::codegen::codegen_provider_api::MatrixBackend;
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend;
    use nalgebra::DMatrix;
    use std::collections::HashMap;
    use std::error::Error;
    use std::path::PathBuf;

    #[test]
    fn checked_in_physical_bvp_documents_load_and_solve() {
        for file in ["hmx_compact.task", "hmx_compact_two_percent_product.task"] {
            let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("task_docs/reactors/bvp")
                .join(file);
            let mut reactor = super::super::SimpleReactorBVP::SimpleReactorTask::new();
            reactor
                .solve_from_file(path)
                .unwrap_or_else(|error| panic!("{file}: file-to-solution story failed: {error}"));

            let data = reactor
                .solution_render_data()
                .unwrap_or_else(|error| panic!("{file}: missing solution data: {error}"));
            assert!(data.x_mesh.iter().all(|value| value.is_finite()), "{file}");
            assert!(
                data.solution.iter().all(|value| value.is_finite()),
                "{file}"
            );
            let temperature_index = data
                .unknowns
                .iter()
                .position(|name| name == "Teta")
                .expect("BVP state should include temperature");
            let temperatures = data.solution.column(temperature_index);
            assert!(
                temperatures
                    .iter()
                    .all(|temperature| (0.0..=3000.0).contains(temperature)),
                "{file}: physical temperature outside [0, 3000] K: {temperatures:?}"
            );
            for (index, name) in data.unknowns.iter().enumerate() {
                if name.starts_with('C') {
                    assert!(
                        data.solution
                            .column(index)
                            .iter()
                            .all(|fraction| (0.0..=1.0).contains(fraction)),
                        "{file}: invalid mass fraction in {name}"
                    );
                }
            }
            println!(
                "BVP task-doc {file}: solved {} species / {} reactions, T={:.2}..{:.2} K, points={}",
                reactor.kindata.substances.len(),
                reactor.kindata.equations().len(),
                temperatures.iter().copied().fold(f64::INFINITY, f64::min),
                temperatures
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max),
                data.x_mesh.len()
            );
        }
    }

    #[test]
    fn correction_newton_physical_task_file_solves_and_reports_telemetry() {
        let source_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("task_docs/reactors/bvp/hmx_compact.task");
        let source = std::fs::read_to_string(&source_path).expect("read compact HMX task");
        let settings_start = source
            .find("solver_settings\n")
            .expect("checked-in task has solver settings");
        let suffix = &source[settings_start..];
        let bounds_start = suffix
            .find("\nbounds\n")
            .map(|relative| settings_start + relative + 1)
            .expect("checked-in task has grouped bounds");
        let correction_settings = "solver_settings\nsolver_family: correction_newton\nassembly: expr_legacy\nexecution: lambdify\nmatrix_layout: sparse\nabsolute_tolerance: 1e-7\ndefault_relative_tolerance: 1e-7\nmax_iterations: 100\nmax_jacobian_age: 5\nmax_damping_trials: 7\ntelemetry: detailed\nloglevel: off\n\n";
        let task_text = format!(
            "{}{}{}",
            &source[..settings_start],
            correction_settings,
            &source[bounds_start..]
        );
        let temp_dir = tempfile::tempdir().expect("temporary task directory");
        let task_path = temp_dir.path().join("hmx_correction_newton.task");
        std::fs::write(&task_path, task_text).expect("write CorrectionNewton task document");

        let mut reactor = super::super::SimpleReactorBVP::SimpleReactorTask::new();
        reactor
            .solve_from_file(task_path)
            .expect("CorrectionNewton should solve the physical HMX task file");

        let data = reactor
            .solution_render_data()
            .expect("CorrectionNewton should publish renderable task data");
        assert!(data.solution.iter().all(|value| value.is_finite()));
        assert!(data.x_mesh.len() >= 2);
        let teta = data
            .unknowns
            .iter()
            .position(|name| name == "Teta")
            .expect("temperature state");
        assert!(data.solution.column(teta).iter().all(|temperature| {
            temperature.is_finite() && *temperature > 0.0 && *temperature < 3000.0
        }));
        let Some(RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::CorrectionNewton(
            Some(telemetry),
        )) = reactor.solver.last_telemetry.as_ref()
        else {
            panic!("task-file CorrectionNewton solve should retain native telemetry");
        };
        assert!(telemetry.residual_evaluations > 0);
        assert!(telemetry.jacobian_evaluations > 0);
        assert!(telemetry.linear_solves > 0);
        let successful_counts = (
            telemetry.residual_evaluations,
            telemetry.jacobian_evaluations,
            telemetry.linear_solves,
        );
        let previous_solution = reactor.solver.solution.clone().expect("published profile");
        let previous_mesh = reactor.solver.x_mesh.clone().expect("published mesh");
        let mut failing_options =
            RustedSciThe::numerical::BVP_CorrectionNewton::CorrectionNewtonOptions {
                absolute_tolerance: 1e-16,
                default_relative_tolerance: 1e-16,
                max_iterations: 1,
                max_damping_trials: 1,
                telemetry_mode: RustedSciThe::numerical::BVP_CorrectionNewton::CorrectionNewtonTelemetryMode::Counters,
                ..Default::default()
            };
        failing_options.relative_tolerances = reactor
            .solver
            .unknowns
            .iter()
            .map(|state| (state.clone(), 1e-16))
            .collect();
        failing_options.bounds = reactor
            .solver
            .unknowns
            .iter()
            .map(|state| {
                let interval = if state == "Teta" {
                    (-100.0, 100.0)
                } else if state == "q" || state.starts_with('J') {
                    (-1e20, 1e20)
                } else {
                    (-0.1, 1.1)
                };
                (state.clone(), interval)
            })
            .collect();
        let state_count = reactor.solver.unknowns.len();
        let failed_retry = reactor.solver.solve_with_hybrid_settings(
            RustedSciThe::numerical::BVP_api::BvpSolverSettings::CorrectionNewton(failing_options),
            DMatrix::from_element(state_count, 32, 0.5),
            0.0,
            1.0,
            32,
        );
        assert!(matches!(failed_retry, Err(ReactorError::BvpApi(_))));
        assert_eq!(reactor.solver.solution.as_ref(), Some(&previous_solution));
        assert_eq!(reactor.solver.x_mesh.as_ref(), Some(&previous_mesh));
        let Some(RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::CorrectionNewton(
            Some(failed_telemetry),
        )) = reactor.solver.last_telemetry.as_ref()
        else {
            panic!("failed CorrectionNewton retry should retain partial telemetry");
        };
        assert!(failed_telemetry.residual_evaluations > 0);
        println!(
            "CorrectionNewton task-file: nodes={} T={:.3}..{:.3} K residuals={} jacobians={} solves={}; failed retry residuals={} retained prior result",
            data.x_mesh.len(),
            data.solution
                .column(teta)
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min),
            data.solution
                .column(teta)
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
            successful_counts.0,
            successful_counts.1,
            successful_counts.2,
            failed_telemetry.residual_evaluations,
        );
    }

    /// End-to-end story test for the modern default reactor BVP setup.
    #[test]
    fn story_compact_hmx_setup_uses_the_default_facade() {
        let mut reactor = compact_hmx_reactor();
        reactor
            .setup_bvp()
            .expect("compact HMX setup should succeed");
        let initial_guess = DMatrix::from_element(reactor.solver.unknowns.len(), 4, 0.5);
        let config = NrbvpHandoffConfig::new(
            initial_guess,
            0.0,
            1.0,
            4,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Sparse".to_string(),
            1e-8,
            None,
            50,
            None,
            Some("info".to_string()),
        )
        .with_solver_backend_config(ReactorBvpSolverConfig::default_lambdify());
        let backend = reactor
            .solver
            .build_nrbvp_backend(config)
            .expect("backend snapshot should build");

        assert_eq!(
            backend.generated_backend_config().symbolic_assembly_backend,
            BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            backend.generated_backend_config().matrix_backend_override,
            Some(MatrixBackend::Banded)
        );
        assert_eq!(
            reactor.solver.unknowns,
            vec!["Teta", "q", "C0", "J0", "C1", "J1"]
        );
        assert_eq!(reactor.kindata.state, KinDataState::Analyzed);
    }

    /// Story-level handoff test that keeps the solver facade and the backend snapshot aligned.
    #[test]
    fn story_build_nrbvp_backend_uses_facade_generated_options() {
        let mut reactor = compact_hmx_reactor();
        reactor
            .setup_bvp()
            .expect("compact HMX setup should succeed");
        let initial_guess = DMatrix::from_element(reactor.solver.unknowns.len(), 4, 0.5);
        let config = NrbvpHandoffConfig::new(
            initial_guess,
            0.0,
            1.0,
            4,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Sparse".to_string(),
            1e-8,
            None,
            50,
            None,
            Some("info".to_string()),
        )
        .with_solver_backend_config(ReactorBvpSolverConfig::default_lambdify());

        let backend = reactor
            .solver
            .build_nrbvp_backend(config)
            .expect("backend snapshot should build");

        assert_eq!(
            backend.generated_backend_config().symbolic_assembly_backend,
            BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            backend.generated_backend_config().matrix_backend_override,
            Some(MatrixBackend::Banded)
        );
        assert_eq!(reactor.kindata.state, KinDataState::Analyzed);
    }

    #[test]
    fn story_compact_reactor_solves_through_damped_bvp_api_with_boundary_conditions() {
        let mut reactor = compact_hmx_reactor();
        reactor
            .setup_bvp()
            .expect("compact reactor should prepare through the production setup");
        let state_count = reactor.solver.unknowns.len();
        let intervals = 32;
        let substances = ["HMX".to_string(), "DecompositionProduct".to_string()];
        let config = NrbvpHandoffConfig::new(
            DMatrix::from_element(state_count, intervals, 0.5),
            0.0,
            1.0,
            intervals,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Banded".to_string(),
            1e-8,
            Some(ToleranceConfig::new(1e-8, 1e-8, 1e-8, 1e-8).to_full_tolerance_map(&substances)),
            80,
            Some(
                BoundsConfig::new((-0.1, 1.1), (-1e20, 1e20), (-100.0, 100.0), (-1e20, 1e20))
                    .to_full_bounds_map(&substances),
            ),
            Some("off".to_string()),
        );

        reactor
            .solver
            .solve_with_bvp_api(config)
            .expect("compact reactor should solve through BVP_api");

        let solution = reactor
            .solver
            .get_solution()
            .expect("successful facade solve should publish a solution");
        let mesh = reactor
            .solver
            .x_mesh
            .as_ref()
            .expect("successful facade solve should publish its mesh");
        assert_eq!(solution.shape(), (mesh.len(), state_count));
        assert_eq!(mesh.len(), intervals + 1);
        assert!(solution.iter().all(|value| value.is_finite()));
        for (state, name) in reactor.solver.unknowns.iter().enumerate() {
            let &(endpoint, expected) = reactor
                .solver
                .BorderConditions
                .get(name)
                .expect("each state should have one endpoint condition");
            let row = if endpoint == 0 { 0 } else { mesh.len() - 1 };
            assert!(
                (solution[(row, state)] - expected).abs() < 1e-5,
                "boundary residual for {name}: actual={}, expected={expected}",
                solution[(row, state)]
            );
        }
    }

    #[test]
    fn story_compact_reactor_frozen_and_naive_publish_full_nodal_profiles() {
        let mut reactor = compact_hmx_reactor();
        reactor
            .setup_bvp()
            .expect("compact reactor should prepare through the production setup");
        let state_count = reactor.solver.unknowns.len();
        let intervals = 32;
        let nodes = intervals + 1;
        let substances = ["HMX".to_string(), "DecompositionProduct".to_string()];
        let damped_config = NrbvpHandoffConfig::new(
            DMatrix::from_element(state_count, intervals, 0.5),
            0.0,
            1.0,
            intervals,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Banded".to_string(),
            1e-8,
            Some(ToleranceConfig::new(1e-8, 1e-8, 1e-8, 1e-8).to_full_tolerance_map(&substances)),
            80,
            Some(
                BoundsConfig::new((-0.1, 1.1), (-1e20, 1e20), (-100.0, 100.0), (-1e20, 1e20))
                    .to_full_bounds_map(&substances),
            ),
            Some("off".to_string()),
        );
        reactor
            .solver
            .solve_with_bvp_api(damped_config)
            .expect("Damped should provide a physical profile to seed Frozen");
        let damped_profile = reactor
            .solver
            .get_solution()
            .expect("Damped solve should publish its profile")
            .clone();
        let initial_guess = DMatrix::from_fn(state_count, nodes, |state, node| {
            damped_profile[(node, state)]
        });
        reactor.solver = reactor.solver.with_bvp_telemetry_mode(
            RustedSciThe::numerical::BVP_api::BvpApiTelemetryMode::Counters,
        );

        reactor
            .solver
            .solve_BVPfrozen(
                initial_guess.clone(),
                FrozenSolverOptions::default()
                    .with_tolerance(1e-6)
                    .with_max_iterations(500),
            )
            .expect("Frozen should publish its converged full nodal profile");

        let frozen_profile = reactor
            .solver
            .get_solution()
            .expect("accepted Frozen profile should be available");
        let frozen_mesh = reactor
            .solver
            .x_mesh
            .as_ref()
            .expect("accepted Frozen mesh should be available");
        assert_eq!(frozen_profile.shape(), (nodes, state_count));
        assert_eq!(frozen_mesh.len(), nodes);
        assert!(frozen_profile.iter().all(|value| value.is_finite()));
        for (state, name) in reactor.solver.unknowns.iter().enumerate() {
            let &(endpoint, expected) = reactor
                .solver
                .BorderConditions
                .get(name)
                .expect("each state should have one endpoint condition");
            let row = if endpoint == 0 {
                0
            } else {
                frozen_mesh.len() - 1
            };
            assert!(
                (frozen_profile[(row, state)] - expected).abs() < 1e-5,
                "Frozen boundary residual for {name}: actual={}, expected={expected}",
                frozen_profile[(row, state)]
            );
        }
        let Some(RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Frozen(snapshot)) =
            reactor.solver.last_telemetry.as_ref()
        else {
            panic!("accepted Frozen attempt should retain its route-tagged telemetry");
        };
        assert_eq!(
            snapshot.telemetry_mode,
            RustedSciThe::numerical::BVP_Damp::telemetry::BvpTelemetryMode::Counters
        );
        assert!(snapshot.counters.iterations > 0);
        assert_eq!(
            snapshot.last_solve_status,
            RustedSciThe::numerical::BVP_Damp::telemetry::BvpSolveStatus::Converged
        );

        reactor
            .solver
            .solve_BVPfrozen(
                initial_guess,
                FrozenSolverOptions::dense_naive()
                    .with_tolerance(1e-6)
                    .with_max_iterations(500),
            )
            .expect("Naive should publish the reconstructed full nodal profile");
        let naive_profile = reactor
            .solver
            .get_solution()
            .expect("accepted Naive profile should be available");
        let naive_mesh = reactor
            .solver
            .x_mesh
            .as_ref()
            .expect("accepted Naive mesh should be available");
        assert_eq!(naive_profile.shape(), (nodes, state_count));
        assert_eq!(naive_mesh.len(), nodes);
        assert!(naive_profile.iter().all(|value| value.is_finite()));
        for (state, name) in reactor.solver.unknowns.iter().enumerate() {
            let &(endpoint, expected) = reactor
                .solver
                .BorderConditions
                .get(name)
                .expect("each state should have one endpoint condition");
            let row = if endpoint == 0 {
                0
            } else {
                naive_mesh.len() - 1
            };
            assert_eq!(naive_profile[(row, state)], expected);
        }
        let Some(RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Frozen(snapshot)) =
            reactor.solver.last_telemetry.as_ref()
        else {
            panic!("accepted Naive attempt should retain its route-tagged telemetry");
        };
        assert_eq!(
            snapshot.telemetry_mode,
            RustedSciThe::numerical::BVP_Damp::telemetry::BvpTelemetryMode::Counters
        );
        assert!(snapshot.counters.iterations > 0);
        assert_eq!(
            snapshot.last_solve_status,
            RustedSciThe::numerical::BVP_Damp::telemetry::BvpSolveStatus::Converged
        );
    }

    #[test]
    fn story_compact_reactor_solves_through_sci_bvp_api_with_boundary_conditions() {
        let mut reactor = compact_hmx_reactor();
        reactor
            .setup_bvp()
            .expect("compact reactor should prepare through the production setup");
        let state_count = reactor.solver.unknowns.len();
        let nodes = 32;
        reactor.solver = reactor.solver.with_bvp_telemetry_mode(
            RustedSciThe::numerical::BVP_api::BvpApiTelemetryMode::Counters,
        );

        reactor
            .solver
            .solve_BVPsci(
                DMatrix::from_element(state_count, nodes, 0.5),
                nodes,
                128,
                1e-5,
            )
            .expect("compact reactor should solve through the BVP_sci route");

        let solution = reactor
            .solver
            .get_solution()
            .expect("successful Sci solve should publish a solution");
        let mesh = reactor
            .solver
            .x_mesh
            .as_ref()
            .expect("successful Sci solve should publish its mesh");
        assert_eq!(solution.shape(), (mesh.len(), state_count));
        assert!(mesh.len() >= nodes);
        assert!(solution.iter().all(|value| value.is_finite()));
        for (state, name) in reactor.solver.unknowns.iter().enumerate() {
            let &(endpoint, expected) = reactor
                .solver
                .BorderConditions
                .get(name)
                .expect("each state should have one endpoint condition");
            let row = if endpoint == 0 { 0 } else { mesh.len() - 1 };
            assert!(
                (solution[(row, state)] - expected).abs() < 1e-5,
                "boundary residual for {name}: actual={}, expected={expected}",
                solution[(row, state)]
            );
        }
    }

    #[test]
    fn story_compact_reactor_inlet_composition_is_accepted_by_damped_and_sci() {
        let mut reactor = compact_hmx_reactor_with_product_fraction(0.02);
        reactor
            .setup_bvp()
            .expect("alternate inlet composition should prepare");
        let state_count = reactor.solver.unknowns.len();
        let intervals = 32;
        let substances = ["HMX".to_string(), "DecompositionProduct".to_string()];
        let damped_config = NrbvpHandoffConfig::new(
            DMatrix::from_element(state_count, intervals, 0.5),
            0.0,
            1.0,
            intervals,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Banded".to_string(),
            1e-8,
            Some(ToleranceConfig::new(1e-8, 1e-8, 1e-8, 1e-8).to_full_tolerance_map(&substances)),
            80,
            Some(
                BoundsConfig::new((-0.1, 1.1), (-1e20, 1e20), (-100.0, 100.0), (-1e20, 1e20))
                    .to_full_bounds_map(&substances),
            ),
            Some("off".to_string()),
        );

        reactor
            .solver
            .solve_with_bvp_api(damped_config)
            .expect("Damped should solve the alternate inlet composition");
        let damped = reactor
            .solver
            .get_solution()
            .expect("Damped should publish its profile")
            .clone();
        let damped_mesh = reactor
            .solver
            .x_mesh
            .as_ref()
            .expect("Damped should publish its mesh")
            .clone();
        assert_eq!(damped.shape(), (damped_mesh.len(), state_count));
        assert!(damped.iter().all(|value| value.is_finite()));
        for (state, name) in reactor.solver.unknowns.iter().enumerate() {
            let &(endpoint, expected) = reactor
                .solver
                .BorderConditions
                .get(name)
                .expect("each state should have an endpoint condition");
            let row = if endpoint == 0 {
                0
            } else {
                damped_mesh.len() - 1
            };
            assert!(
                (damped[(row, state)] - expected).abs() < 1e-5,
                "Damped boundary residual for {name}: actual={}, expected={expected}",
                damped[(row, state)]
            );
        }

        reactor
            .solver
            .solve_BVPsci(
                DMatrix::from_element(state_count, intervals, 0.5),
                intervals,
                128,
                1e-5,
            )
            .expect("BVP_sci should solve the alternate inlet composition");
        let sci = reactor
            .solver
            .get_solution()
            .expect("BVP_sci should publish its profile");
        let sci_mesh = reactor
            .solver
            .x_mesh
            .as_ref()
            .expect("BVP_sci should publish its mesh");
        assert_eq!(sci.ncols(), state_count);
        assert!(sci.nrows() >= intervals);
        assert!(sci.iter().all(|value| value.is_finite()));
        for (state, name) in reactor.solver.unknowns.iter().enumerate() {
            let &(endpoint, expected) = reactor
                .solver
                .BorderConditions
                .get(name)
                .expect("each state should have an endpoint condition");
            let row = if endpoint == 0 { 0 } else { sci_mesh.len() - 1 };
            assert!(
                (sci[(row, state)] - expected).abs() < 1e-5,
                "BVP_sci boundary residual for {name}: actual={}, expected={expected}",
                sci[(row, state)]
            );
        }
    }

    /// The facade publishes a shared API result and retains route-native diagnostics.
    #[test]
    fn story_scalar_bvp_solves_through_bvp_api_and_retains_native_telemetry() {
        let mut solver = super::super::SimpleReactorBVP::BVPSolver::default()
            .with_bvp_telemetry_mode(
                RustedSciThe::numerical::BVP_api::BvpApiTelemetryMode::Counters,
            );
        solver.arg_name = "x".to_string();
        solver.x_range = (0.0, 1.0);
        solver.unknowns = vec!["y".to_string()];
        solver.eq_system = vec![Expr::parse_expression("-y")];
        solver.BorderConditions.insert("y".to_string(), (0, 1.0));
        let state_count = solver.unknowns.len();
        let intervals = 32;
        let config = NrbvpHandoffConfig::new(
            DMatrix::from_element(state_count, intervals, 0.5),
            0.0,
            1.0,
            intervals,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Banded".to_string(),
            1e-8,
            Some(HashMap::from([("y".to_string(), 1e-8)])),
            80,
            Some(HashMap::from([("y".to_string(), (-10.0, 10.0))])),
            Some("off".to_string()),
        );

        solver
            .solve_with_bvp_api(config)
            .expect("scalar BVP should solve through BVP_api");

        let solution = solver
            .get_solution()
            .expect("successful BVP_api solve should publish a solution");
        assert_eq!(solution.nrows(), intervals + 1);
        assert_eq!(solution.ncols(), state_count);
        assert!(solution.iter().all(|value| value.is_finite()));
        let mesh = solver
            .x_mesh
            .as_ref()
            .expect("successful BVP_api solve should publish its mesh");
        assert_eq!(mesh.len(), intervals + 1);
        assert_eq!(mesh[0], 0.0);
        assert_eq!(mesh[intervals], 1.0);
        let Some(RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Damped(telemetry)) =
            solver.last_telemetry.as_ref()
        else {
            panic!("expected the route-tagged Damped telemetry snapshot");
        };
        assert_eq!(telemetry.telemetry_mode, BvpTelemetryMode::Counters);
        assert!(telemetry.counters.residual_calls > 0);

        let previous_solution = solution.clone();
        let previous_mesh = mesh.clone();
        let invalid_retry = NrbvpHandoffConfig::new(
            DMatrix::from_element(state_count, 1, 0.5),
            0.0,
            1.0,
            1,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Banded".to_string(),
            1e-8,
            None,
            20,
            None,
            None,
        );
        assert!(solver.solve_with_bvp_api(invalid_retry).is_err());
        assert_eq!(solver.get_solution(), Some(&previous_solution));
        assert_eq!(solver.x_mesh.as_ref(), Some(&previous_mesh));
        assert!(solver.last_telemetry.is_none());
    }

    #[test]
    fn story_scalar_sci_bvp_uses_shared_api_and_publishes_adaptive_mesh() {
        let mut solver = super::super::SimpleReactorBVP::BVPSolver::default()
            .with_bvp_telemetry_mode(
                RustedSciThe::numerical::BVP_api::BvpApiTelemetryMode::Counters,
            );
        solver.arg_name = "x".to_string();
        solver.x_range = (0.0, 1.0);
        solver.unknowns = vec!["y".to_string()];
        solver.eq_system = vec![Expr::parse_expression("-y")];
        solver.BorderConditions.insert("y".to_string(), (0, 1.0));
        solver
            .solve_BVPsci(DMatrix::from_element(1, 16, 0.8), 16, 128, 1e-5)
            .expect("scalar collocation BVP should solve through BVP_api");

        let solution = solver
            .get_solution()
            .expect("successful BVP_sci route should publish a solution");
        let mesh = solver
            .x_mesh
            .as_ref()
            .expect("successful BVP_sci route should publish its adaptive mesh");
        assert_eq!(solution.nrows(), mesh.len());
        assert_eq!(solution.ncols(), 1);
        assert!(mesh.len() >= 16);
        assert!((mesh[0] - 0.0).abs() < 1e-14);
        assert!((mesh[mesh.len() - 1] - 1.0).abs() < 1e-14);
        assert!(solution.iter().all(|value| value.is_finite()));
        assert!((solution[(solution.nrows() - 1, 0)] - (-1.0f64).exp()).abs() < 2e-3);
        let Some(RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Sci(snapshot)) =
            solver.last_telemetry.as_ref()
        else {
            panic!("expected the route-tagged Sci telemetry snapshot");
        };
        assert_eq!(
            snapshot.mode,
            RustedSciThe::numerical::BVP_sci::new::BvpSciTelemetryMode::Counters
        );
        assert!(snapshot.rhs_calls > 0);

        solver = solver
            .with_bvp_telemetry_mode(RustedSciThe::numerical::BVP_api::BvpApiTelemetryMode::Off);
        solver
            .solve_BVPsci(DMatrix::from_element(1, 16, 0.8), 16, 128, 1e-5)
            .expect("explicit Off mode should preserve the Sci solve");
        let Some(RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Sci(snapshot)) =
            solver.last_telemetry.as_ref()
        else {
            panic!("expected the route-tagged Sci telemetry snapshot");
        };
        assert_eq!(
            snapshot.mode,
            RustedSciThe::numerical::BVP_sci::new::BvpSciTelemetryMode::Off
        );
        assert_eq!(snapshot.rhs_calls, 0);
        assert_eq!(snapshot.jacobian_calls, 0);
    }

    #[test]
    fn story_scalar_frozen_bvp_publishes_boundary_aligned_profile() {
        let mut solver = super::super::SimpleReactorBVP::BVPSolver::default()
            .with_bvp_telemetry_mode(
                RustedSciThe::numerical::BVP_api::BvpApiTelemetryMode::Counters,
            );
        solver.arg_name = "x".to_string();
        solver.x_range = (0.0, 1.0);
        solver.unknowns = vec!["y".to_string()];
        solver.eq_system = vec![Expr::parse_expression("-y")];
        solver.BorderConditions.insert("y".to_string(), (0, 1.0));
        let nodes = 257;
        let options = FrozenSolverOptions::default().with_tolerance(1e-6);

        solver
            .solve_BVPfrozen(DMatrix::from_element(1, nodes, 1.0), options)
            .expect("Frozen should publish the reconstructed full nodal profile");
        let profile = solver
            .get_solution()
            .expect("accepted Frozen profile should be available");
        let mesh = solver
            .x_mesh
            .as_ref()
            .expect("accepted Frozen mesh should be available");
        assert_eq!(profile.shape(), (nodes, 1));
        assert_eq!(mesh.len(), nodes);
        assert_eq!(mesh[0], 0.0);
        assert_eq!(mesh[nodes - 1], 1.0);
        assert_eq!(profile[(0, 0)], 1.0);
        assert!(profile.iter().all(|value| value.is_finite()));
        let Some(RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Frozen(snapshot)) =
            solver.last_telemetry.as_ref()
        else {
            panic!("accepted Frozen attempt should retain route telemetry");
        };
        assert_eq!(
            snapshot.telemetry_mode,
            RustedSciThe::numerical::BVP_Damp::telemetry::BvpTelemetryMode::Counters
        );
        assert!(snapshot.counters.iterations > 0);
        assert_eq!(
            snapshot.last_solve_status,
            RustedSciThe::numerical::BVP_Damp::telemetry::BvpSolveStatus::Converged
        );
    }

    #[test]
    fn frozen_and_naive_routes_match_across_frontends_and_matrix_layouts() {
        let nodes = 33;
        let mut reference: Option<Vec<f64>> = None;

        for (strategy, naive) in [("Frozen", false), ("Naive", true)] {
            for frontend in [
                BvpSymbolicAssemblyBackend::ExprLegacy,
                BvpSymbolicAssemblyBackend::AtomView,
            ] {
                for (layout, options) in [
                    ("Dense", FrozenSolverOptions::dense_frozen()),
                    ("Sparse", FrozenSolverOptions::sparse_frozen()),
                    (
                        "Banded",
                        FrozenSolverOptions::banded_frozen().with_banded_lambdify(),
                    ),
                ] {
                    let mut options = options
                        .with_tolerance(1e-8)
                        .with_max_iterations(500)
                        .with_bvp_telemetry_mode(BvpTelemetryMode::Counters)
                        .with_symbolic_assembly_backend(frontend);
                    if naive {
                        options.strategy = "Naive".to_string();
                        options.strategy_params = None;
                    }
                    options = options.with_matrix_backend(match layout {
                        "Dense" => MatrixBackend::Dense,
                        "Sparse" => MatrixBackend::SparseCol,
                        "Banded" => MatrixBackend::Banded,
                        _ => unreachable!(),
                    });
                    let mut solver = super::super::SimpleReactorBVP::BVPSolver::default();
                    solver.arg_name = "x".to_string();
                    solver.x_range = (0.0, 1.0);
                    solver.unknowns = vec!["y".to_string()];
                    solver.eq_system = vec![Expr::parse_expression("-0.7 * y")];
                    solver.BorderConditions.insert("y".to_string(), (0, 1.0));
                    solver = solver.with_bvp_telemetry_mode(
                        RustedSciThe::numerical::BVP_api::BvpApiTelemetryMode::Counters,
                    );
                    solver
                        .solve_BVPfrozen(DMatrix::from_element(1, nodes, 0.8), options)
                        .unwrap_or_else(|error| {
                            panic!("{strategy}/{frontend:?}/{layout} solve: {error}")
                        });
                    let published = solver
                        .get_solution()
                        .expect("accepted facade route should publish a profile");
                    let published_mesh = solver
                        .x_mesh
                        .as_ref()
                        .expect("accepted facade route should publish its mesh");
                    assert_eq!(published.shape(), (nodes, 1));
                    assert_eq!(published_mesh.len(), nodes);
                    let profile = published_mesh
                        .iter()
                        .enumerate()
                        .map(|(node, x)| {
                            let actual = published[(node, 0)];
                            assert!(actual.is_finite() && (0.0..=1.0).contains(&actual));
                            let error = (actual - (-0.7 * x).exp()).abs();
                            assert!(
                                error <= 4e-2,
                                "{strategy}/{frontend:?}/{layout} at x={x}: error={error:e}"
                            );
                            actual
                        })
                        .collect::<Vec<_>>();
                    if let Some(expected) = &reference {
                        let max_delta = profile
                            .iter()
                            .zip(expected)
                            .map(|(actual, expected)| (actual - expected).abs())
                            .fold(0.0_f64, f64::max);
                        println!(
                            "BVP {strategy}/{frontend:?}/{layout} max_abs_profile_delta={max_delta:.6e}"
                        );
                        assert!(max_delta <= 4e-2, "cross-route profile drift {max_delta:e}");
                    } else {
                        reference = Some(profile);
                    }
                    assert!(matches!(
                        solver.last_telemetry,
                        Some(
                            RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Frozen(_)
                        )
                    ));
                }
            }
        }
    }

    #[test]
    fn story_damped_api_preserves_rst_configuration_error_source() {
        let mut solver = super::super::SimpleReactorBVP::BVPSolver::default();
        solver.unknowns = vec!["y".to_string()];
        solver.eq_system = vec![Expr::parse_expression("-y")];
        solver.BorderConditions.insert("y".to_string(), (0, 1.0));
        let config = NrbvpHandoffConfig::new(
            DMatrix::from_element(1, 2, 0.5),
            0.0,
            1.0,
            2,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Banded".to_string(),
            1e-8,
            None,
            20,
            None,
            None,
        );

        let error = solver
            .solve_with_bvp_api(config)
            .expect_err("RST should reject missing Damped bounds and tolerances");
        assert!(matches!(error, ReactorError::BvpApi(_)));
        let source =
            Error::source(&error).expect("native BVP error should remain in the source chain");
        assert!(source.to_string().contains("Damped"));
        assert!(solver.get_solution().is_none());
    }
}
