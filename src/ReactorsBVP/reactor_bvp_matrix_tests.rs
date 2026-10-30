#[cfg(test)]
mod tests {
    use super::super::SimpleReactorBVP::{
        FastElemReact, NrbvpHandoffConfig, SimpleReactorTask, SymbolicRhsAssemblyPolicy,
    };
    use super::super::reactor_bvp_test_support::compact_hmx_reactor;
    use super::super::solver_backend::{
        ReactorBvpAotCompiler, ReactorBvpMatrixBackend, ReactorBvpSolverConfig,
        ReactorBvpSymbolicBackend,
    };
    use crate::ReactorsBVP::SimpleReactorBVP::ReactorError;
    use crate::ReactorsBVP::reactor_BVP_utils::{BoundsConfig, ScalingConfig, ToleranceConfig};
    use RustedSciThe::numerical::BVP_api::{
        BvpBoundaryCondition, BvpBoundarySide, BvpProblem, BvpSciSparseFactorization, BvpSolver,
        BvpSolverConfig,
    };
    use RustedSciThe::numerical::BVP_sci::new::{
        BvpSciAssembly, BvpSciExecution, BvpSciMatrixLayout, BvpSciOptions, BvpSciTelemetryMode,
    };
    use RustedSciThe::symbolic::bvp::telemetry::BvpLambdifyExecutionPolicy;
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use RustedSciThe::symbolic::symbolic_ivp_generated::SymbolicIvpGeneratedBackendConfig;
    use nalgebra::DMatrix;
    use std::process::Command;

    /// Build a compact initial guess for the tiny HMX combustion fixture.
    fn compact_initial_guess(unknowns: usize, n_steps: usize) -> DMatrix<f64> {
        DMatrix::from_element(unknowns, n_steps, 0.5)
    }

    /// Build the legacy NRBVP handoff wrapper with a selected solver facade.
    fn compact_hmx_handoff(
        initial_guess: DMatrix<f64>,
        solver_backend_config: ReactorBvpSolverConfig,
        n_steps: usize,
    ) -> NrbvpHandoffConfig {
        NrbvpHandoffConfig::new(
            initial_guess,
            0.0,
            1.0,
            n_steps,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Sparse".to_string(),
            1e-8,
            Some(
                ToleranceConfig::new(1e-8, 1e-8, 1e-8, 1e-8)
                    .to_full_tolerance_map(&["HMX".to_string(), "HMXprod".to_string()]),
            ),
            80,
            Some(
                BoundsConfig::new((-0.1, 1.1), (-1e20, 1e20), (-100.0, 100.0), (-1e20, 1e20))
                    .to_full_bounds_map(&["HMX".to_string(), "HMXprod".to_string()]),
            ),
            Some("warn".to_string()),
            false,
        )
        .with_solver_backend_config(solver_backend_config)
    }

    /// Solve the compact combustion fixture through the NRBVP backend snapshot.
    fn solve_compact_backend(
        solver_backend_config: ReactorBvpSolverConfig,
    ) -> Result<(Vec<f64>, Vec<f64>), ReactorError> {
        solve_compact_backend_on_grid(solver_backend_config, 32)
    }

    fn solve_compact_backend_on_grid(
        solver_backend_config: ReactorBvpSolverConfig,
        n_steps: usize,
    ) -> Result<(Vec<f64>, Vec<f64>), ReactorError> {
        let mut reactor = compact_hmx_reactor();
        reactor.setup_bvp()?;

        let initial_guess = compact_initial_guess(reactor.solver.unknowns.len(), n_steps);
        let mut backend = reactor.solver.build_nrbvp_backend(compact_hmx_handoff(
            initial_guess,
            solver_backend_config,
            n_steps,
        ))?;

        backend.before_solve_preprocessing();
        let solve_result = backend
            .try_solve()
            .map_err(|err| ReactorError::CalculationError(format!("{err:?}")))?;
        if solve_result.is_none() {
            return Err(ReactorError::CalculationError(
                "compact backend solve did not converge".to_string(),
            ));
        }

        let solution = backend.get_result().ok_or_else(|| {
            ReactorError::CalculationError(
                "compact backend solve did not return a solution".to_string(),
            )
        })?;
        if solution.nrows() == 0 || solution.ncols() == 0 {
            return Err(ReactorError::CalculationError(
                "compact backend returned an empty solution".to_string(),
            ));
        }
        if solution.iter().any(|value| !value.is_finite()) {
            return Err(ReactorError::InvalidNumericValue(
                "compact backend returned non-finite values".to_string(),
            ));
        }
        if backend.x_mesh.is_empty() {
            return Err(ReactorError::CalculationError(
                "compact backend returned an empty mesh".to_string(),
            ));
        }
        if backend.x_mesh.iter().any(|value| !value.is_finite()) {
            return Err(ReactorError::InvalidNumericValue(
                "compact backend returned non-finite mesh values".to_string(),
            ));
        }

        Ok((
            solution.iter().copied().collect(),
            backend.x_mesh.iter().copied().collect(),
        ))
    }

    fn solve_compact_sci(
        assembly: BvpSciAssembly,
        layout: BvpSciMatrixLayout,
        n_intervals: usize,
    ) -> Vec<f64> {
        solve_compact_sci_with_aot(assembly, layout, n_intervals, None)
    }

    fn solve_compact_sci_with_aot(
        assembly: BvpSciAssembly,
        layout: BvpSciMatrixLayout,
        n_intervals: usize,
        aot_config: Option<SymbolicIvpGeneratedBackendConfig>,
    ) -> Vec<f64> {
        solve_compact_sci_with_factorization(
            assembly,
            layout,
            n_intervals,
            aot_config,
            BvpSciSparseFactorization::StandardFaer,
        )
        .0
    }

    fn solve_compact_sci_with_factorization(
        assembly: BvpSciAssembly,
        layout: BvpSciMatrixLayout,
        n_intervals: usize,
        aot_config: Option<SymbolicIvpGeneratedBackendConfig>,
        sparse_factorization: BvpSciSparseFactorization,
    ) -> (
        Vec<f64>,
        RustedSciThe::numerical::BVP_sci::new::BvpSciTelemetrySnapshot,
    ) {
        let mut reactor = compact_hmx_reactor();
        reactor
            .setup_bvp()
            .expect("compact HMX equations should assemble");
        let state_count = reactor.solver.unknowns.len();
        let mesh = (0..=n_intervals)
            .map(|index| index as f64 / n_intervals as f64)
            .collect::<Vec<_>>();
        let boundary_conditions = reactor
            .solver
            .unknowns
            .iter()
            .enumerate()
            .map(|(state, name)| {
                let (side, value) = reactor.solver.state_boundary_conditions()[name];
                BvpBoundaryCondition {
                    state,
                    side: if side == 0 {
                        BvpBoundarySide::Left
                    } else {
                        BvpBoundarySide::Right
                    },
                    value,
                }
            })
            .collect();
        let problem = BvpProblem::new(
            reactor.solver.eq_system.clone(),
            reactor.solver.unknowns.clone(),
            reactor.solver.arg_name.clone(),
            mesh.clone(),
            boundary_conditions,
        )
        .expect("compact HMX BVP problem should validate");
        let execution = if aot_config.is_some() {
            BvpSciExecution::Aot
        } else {
            BvpSciExecution::Lambdify
        };
        let options = BvpSciOptions::default()
            .with_execution(execution)
            .with_assembly(assembly)
            .with_matrix_layout(layout)
            .with_tolerance(1e-5)
            .with_boundary_tolerance(1e-5)
            .with_limits(512, 100, 20)
            .with_telemetry(BvpSciTelemetryMode::Timings);
        let mut solver = BvpSolver::new(
            problem,
            BvpSolverConfig::Sci {
                options,
                sparse_factorization,
                aot_config,
                initial_guess: DMatrix::from_element(state_count, mesh.len(), 0.5),
            },
        )
        .expect("compact Sci route should prepare");
        solver.solve().expect("compact Sci route should solve");
        let telemetry = solver.telemetry_snapshot();
        let RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot::Sci(telemetry) =
            telemetry
        else {
            panic!("BVP_sci should retain Sci telemetry")
        };
        let view = solver
            .solution_view()
            .expect("Sci result should be retained");
        assert!(view.mesh().len() >= mesh.len());
        let mut values = Vec::with_capacity(view.mesh().len() * state_count);
        for node in 0..view.mesh().len() {
            for state in 0..state_count {
                values.push(view.value(node, state).expect("valid node/state index"));
            }
        }
        for state in 0..state_count {
            let name = &reactor.solver.unknowns[state];
            let (side, expected) = reactor.solver.state_boundary_conditions()[name];
            let node = if side == 0 { 0 } else { view.mesh().len() - 1 };
            let actual = view.value(node, state).expect("valid boundary state");
            assert!(
                (actual - expected).abs() <= 1e-5,
                "{name}: {actual} vs {expected}"
            );
        }
        assert!(values.iter().all(|value| value.is_finite()));
        (values, telemetry)
    }

    #[derive(Clone, Copy, Debug)]
    enum PolicyTestSolver {
        Damped,
        Sci,
    }

    fn solve_large_callback_policy_case(
        solver_family: PolicyTestSolver,
        policy: BvpLambdifyExecutionPolicy,
    ) -> (
        Vec<f64>,
        RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot,
    ) {
        let state_count = 24;
        let intervals = 128;
        let rates = (0..state_count)
            .map(|index| 0.25 + 0.035 * index as f64)
            .collect::<Vec<_>>();
        let names = (0..state_count)
            .map(|index| format!("y{index}"))
            .collect::<Vec<_>>();
        let equations = rates
            .iter()
            .zip(&names)
            .map(|(rate, name)| Expr::parse_expression(&format!("-({rate:.17}) * {name}")))
            .collect();
        let mesh = (0..=intervals)
            .map(|index| index as f64 / intervals as f64)
            .collect::<Vec<_>>();
        let boundaries = (0..state_count)
            .map(|state| BvpBoundaryCondition {
                state,
                side: BvpBoundarySide::Left,
                value: 1.0,
            })
            .collect();
        let problem = BvpProblem::new(equations, names.clone(), "x", mesh.clone(), boundaries)
            .expect("large analytic policy problem should validate");
        let initial_guess = |columns: usize| {
            DMatrix::from_fn(state_count, columns, |state, column| {
                let x = column as f64 / intervals as f64;
                (-rates[state] * x).exp()
            })
        };
        let config = match solver_family {
            PolicyTestSolver::Damped => {
                let tolerances = names
                    .iter()
                    .map(|name| (name.clone(), 1e-8))
                    .collect::<std::collections::HashMap<_, _>>();
                let bounds = names
                    .iter()
                    .map(|name| (name.clone(), (0.0, 1.1)))
                    .collect::<std::collections::HashMap<_, _>>();
                let options = ReactorBvpSolverConfig::sparse_lambdify()
                    .with_symbolic_backend(ReactorBvpSymbolicBackend::AtomView)
                    .to_rusted_options()
                    .expect("Damped sparse options should map")
                    .with_abs_tolerance(1e-8)
                    .with_max_iterations(100)
                    .with_rel_tolerance(tolerances)
                    .with_bounds(bounds)
                    .with_bvp_telemetry_mode(
                        RustedSciThe::numerical::BVP_Damp::telemetry::BvpTelemetryMode::Detailed,
                    );
                let generated = options
                    .generated_backend_config
                    .clone()
                    .with_lambdify_execution_policy(policy);
                let mut generated = generated;
                generated.lambdify_telemetry_mode =
                    RustedSciThe::symbolic::bvp::telemetry::BvpLambdifyTelemetryMode::Detailed;
                BvpSolverConfig::Damped {
                    options: options.with_generated_backend_config(generated),
                    initial_guess: initial_guess(intervals),
                }
            }
            PolicyTestSolver::Sci => {
                let options = BvpSciOptions::default()
                    .with_execution(BvpSciExecution::Lambdify)
                    .with_assembly(BvpSciAssembly::AtomViewNative)
                    .with_matrix_layout(BvpSciMatrixLayout::Sparse)
                    .with_tolerance(1e-6)
                    .with_boundary_tolerance(1e-7)
                    .with_limits(4 * (intervals + 1), 20, 6)
                    .with_execution_policy(policy)
                    .with_telemetry(BvpSciTelemetryMode::Detailed);
                BvpSolverConfig::Sci {
                    options,
                    sparse_factorization: BvpSciSparseFactorization::StandardFaer,
                    aot_config: None,
                    initial_guess: initial_guess(intervals + 1),
                }
            }
        };
        let mut solver = BvpSolver::new(problem, config).expect("policy solver should prepare");
        solver.solve().expect("policy solver should converge");
        let view = solver.solution_view().expect("solution should be retained");
        let values = (0..view.mesh().len())
            .flat_map(|node| (0..state_count).map(move |state| view.value(node, state).unwrap()))
            .collect();
        (values, solver.telemetry_snapshot())
    }

    fn large_kithe_setup_fixture(policy: SymbolicRhsAssemblyPolicy) -> SimpleReactorTask {
        let species = (1..=40).map(|n| format!("C{n}")).collect::<Vec<_>>();
        let reactions = (1..40)
            .map(|n| FastElemReact {
                eq: if n == 1 {
                    "C1+C1=>C2".to_string()
                } else {
                    format!("C{n}+C1=>C{}", n + 1)
                },
                A: 1.0e3,
                n: 0.0,
                E: 1.0e4,
                Q: 0.0,
            })
            .collect();
        let mut reactor = SimpleReactorTask::new().with_symbolic_rhs_assembly_policy(policy);
        reactor
            .fast_react_set(reactions)
            .expect("large KiThe setup reactions should parse");
        let boundary = species
            .iter()
            .enumerate()
            .map(|(i, name)| (name.clone(), if i == 0 { 1.0 } else { 0.0 }))
            .chain([("T".to_string(), 700.0)])
            .collect();
        let diffusion = species.iter().map(|name| (name.clone(), 1.0e-5)).collect();
        reactor.set_parameters(
            vec![0.0; 39],
            101_325.0,
            700.0,
            1_000.0,
            boundary,
            0.1,
            diffusion,
            0.01,
            ScalingConfig::new(600.0, 1.0e-3, 600.0),
        );
        reactor.M = 0.012;
        reactor
    }

    #[test]
    fn story_kithe_pre_rst_symbolic_policies_preserve_large_equation_system() {
        let policies = [
            SymbolicRhsAssemblyPolicy::Sequential,
            SymbolicRhsAssemblyPolicy::Auto,
            SymbolicRhsAssemblyPolicy::Parallel,
        ];
        let mut reference: Option<Vec<String>> = None;
        for policy in policies {
            let mut reactor = large_kithe_setup_fixture(policy);
            reactor
                .setup_bvp()
                .unwrap_or_else(|error| panic!("{policy:?} KiThe setup failed: {error}"));
            assert_eq!(reactor.solver.eq_system.len(), 82);
            let equations = reactor
                .solver
                .eq_system
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>();
            if let Some(reference) = &reference {
                assert_eq!(reference, &equations, "assembly policy changed equations");
            } else {
                reference = Some(equations);
            }
        }
    }

    /// Assert that two compact-task solution profiles are numerically close.
    fn assert_profile_close(reference: &[f64], candidate: &[f64], atol: f64, label: &str) {
        assert_eq!(
            reference.len(),
            candidate.len(),
            "{label} length mismatch: reference={}, candidate={}",
            reference.len(),
            candidate.len()
        );
        for (idx, (expected, actual)) in reference.iter().zip(candidate.iter()).enumerate() {
            assert!(
                (*expected - *actual).abs() <= atol,
                "{label} mismatch at index {idx}: expected {expected}, got {actual}, atol={atol}"
            );
            assert!(
                actual.is_finite(),
                "{label} contains non-finite value at {idx}"
            );
        }
    }

    /// Check whether the gated AOT integration tests are allowed to run on this machine.
    fn aot_tests_enabled(toolchain: &str) -> bool {
        if std::env::var("KITHE_RUN_BVP_AOT_TESTS").is_err() {
            eprintln!("skipping AOT smoke test: set KITHE_RUN_BVP_AOT_TESTS=1 to enable it");
            return false;
        }
        if Command::new(toolchain).arg("--version").output().is_err() {
            eprintln!("skipping AOT smoke test: {toolchain} is not available on this machine");
            return false;
        }
        true
    }

    /// Compare two compact-task solver configurations and optionally gate them on a toolchain.
    fn assert_config_matches_reference(
        reference_config: ReactorBvpSolverConfig,
        candidate_config: ReactorBvpSolverConfig,
        toolchain: Option<&str>,
        atol: f64,
    ) {
        if let Some(toolchain) = toolchain {
            if !aot_tests_enabled(toolchain) {
                return;
            }
        }

        let (reference_solution, reference_mesh) = solve_compact_backend(reference_config)
            .expect("reference backend solve should succeed");
        let (candidate_solution, candidate_mesh) = solve_compact_backend(candidate_config)
            .expect("candidate backend solve should succeed");

        assert_profile_close(&reference_solution, &candidate_solution, atol, "solution");
        assert_profile_close(&reference_mesh, &candidate_mesh, atol, "mesh");
    }

    /// Default and sparse lambdify routes should differ only where the matrix backend changes.
    #[test]
    fn matrix_default_and_sparse_lambdify_configs_are_distinct() {
        assert_eq!(
            ReactorBvpMatrixBackend::default(),
            ReactorBvpMatrixBackend::Banded
        );

        let default_options = ReactorBvpSolverConfig::default_lambdify()
            .to_rusted_options()
            .expect("default config should convert");
        let sparse_options = ReactorBvpSolverConfig::sparse_lambdify()
            .to_rusted_options()
            .expect("sparse config should convert");

        assert_eq!(
            default_options
                .generated_backend_config
                .symbolic_assembly_backend,
            RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            sparse_options
                .generated_backend_config
                .symbolic_assembly_backend,
            RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            default_options
                .generated_backend_config
                .matrix_backend_override,
            Some(RustedSciThe::symbolic::codegen::codegen_provider_api::MatrixBackend::Banded)
        );
        assert_eq!(
            sparse_options
                .generated_backend_config
                .matrix_backend_override,
            None
        );
    }

    /// Compact combustion should solve through the default lambdify facade.
    #[test]
    fn matrix_compact_hmx_default_lambdify_solves() {
        let (solution, mesh) = solve_compact_backend(ReactorBvpSolverConfig::default_lambdify())
            .expect("default lambdify solve should succeed");

        assert!(!solution.is_empty());
        assert!(!mesh.is_empty());
    }

    /// Sparse lambdify should stay usable on the same compact combustion fixture.
    #[test]
    fn matrix_compact_hmx_sparse_lambdify_solves() {
        let (solution, mesh) = solve_compact_backend(ReactorBvpSolverConfig::sparse_lambdify())
            .expect("sparse lambdify solve should succeed");

        assert!(!solution.is_empty());
        assert!(!mesh.is_empty());
    }

    /// ExprLegacy + Sparse should stay compatible with the modern default profile on the same task.
    #[test]
    fn matrix_compact_hmx_legacy_expr_sparse_matches_default_profile() {
        let legacy_sparse = ReactorBvpSolverConfig::sparse_lambdify()
            .with_symbolic_backend(ReactorBvpSymbolicBackend::ExprLegacy);

        assert_config_matches_reference(
            ReactorBvpSolverConfig::default_lambdify(),
            legacy_sparse,
            None,
            1e-6,
        );
    }

    #[test]
    fn matrix_compact_hmx_damped_frontend_and_matrix_matrix_matches_reference() {
        let (reference, _) = solve_compact_backend(ReactorBvpSolverConfig::default_lambdify())
            .expect("default Damped route should solve the compact reactor");
        for matrix in [
            ReactorBvpMatrixBackend::Banded,
            ReactorBvpMatrixBackend::Sparse,
        ] {
            for frontend in [
                ReactorBvpSymbolicBackend::AtomView,
                ReactorBvpSymbolicBackend::ExprLegacy,
            ] {
                let config = ReactorBvpSolverConfig::default_lambdify()
                    .with_matrix_backend(matrix)
                    .with_symbolic_backend(frontend);
                let (candidate, _) = solve_compact_backend(config).unwrap_or_else(|error| {
                    panic!("Damped {matrix:?}/{frontend:?} failed: {error}")
                });
                assert_profile_close(&reference, &candidate, 1e-6, "Damped route parity");
            }
        }
    }

    #[test]
    fn matrix_compact_hmx_large_grid_remains_correct_on_banded_and_sparse_routes() {
        let (reference, mesh) =
            solve_compact_backend_on_grid(ReactorBvpSolverConfig::default_lambdify(), 96)
                .expect("96-interval compact reactor should solve");
        assert_eq!(mesh.len(), 97);
        assert_eq!(reference.len(), 97 * 6);
        for matrix in [
            ReactorBvpMatrixBackend::Banded,
            ReactorBvpMatrixBackend::Sparse,
        ] {
            let candidate_config = ReactorBvpSolverConfig::default_lambdify()
                .with_matrix_backend(matrix)
                .with_symbolic_backend(ReactorBvpSymbolicBackend::ExprLegacy);
            let (candidate, candidate_mesh) = solve_compact_backend_on_grid(candidate_config, 96)
                .unwrap_or_else(|error| panic!("large-grid {matrix:?} failed: {error}"));
            assert_eq!(candidate_mesh, mesh);
            assert_profile_close(&reference, &candidate, 2e-5, "large-grid matrix parity");
        }
    }

    #[test]
    fn matrix_compact_hmx_sci_frontend_and_layout_matrix_matches_boundaries() {
        let mut reference: Option<Vec<f64>> = None;
        for assembly in [BvpSciAssembly::ExprLegacy, BvpSciAssembly::AtomViewNative] {
            for layout in [
                BvpSciMatrixLayout::Dense,
                BvpSciMatrixLayout::Sparse,
                BvpSciMatrixLayout::BandedAuto,
            ] {
                let candidate = solve_compact_sci(assembly, layout, 32);
                if let Some(reference) = &reference {
                    assert_profile_close(reference, &candidate, 2e-5, "Sci frontend/layout parity");
                } else {
                    reference = Some(candidate);
                }
            }
        }
    }

    #[test]
    fn matrix_compact_hmx_sci_superlu_matches_faer_without_fallback() {
        let bridge = ["SUPERLU_LIBRARY", "SUPERLU_EXPERIMENTAL_LIBRARY"]
            .into_iter()
            .filter_map(|name| std::env::var_os(name).map(std::path::PathBuf::from))
            .find(|path| path.is_file());
        let Some(bridge) = bridge else {
            eprintln!(
                "skipping SuperLU correctness route: set SUPERLU_LIBRARY to the native bridge file"
            );
            return;
        };
        assert!(bridge.is_file());
        let assembly = BvpSciAssembly::AtomViewNative;
        let faer = solve_compact_sci_with_factorization(
            assembly,
            BvpSciMatrixLayout::Sparse,
            32,
            None,
            BvpSciSparseFactorization::StandardFaer,
        );
        let superlu = solve_compact_sci_with_factorization(
            assembly,
            BvpSciMatrixLayout::Sparse,
            32,
            None,
            BvpSciSparseFactorization::SuperLu,
        );
        assert_profile_close(&faer.0, &superlu.0, 2e-5, "SuperLU/faer profile parity");
        assert_eq!(superlu.1.superlu_fallbacks, 0);
        assert!(
            superlu.1.superlu.native_factorization_calls > 0,
            "configured SuperLU route did not call native factorization"
        );
    }

    fn callback_dispatch_counts(
        snapshot: &RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot,
    ) -> (u64, u64) {
        use RustedSciThe::numerical::BVP_api::BvpSolverTelemetrySnapshot;
        match snapshot {
            BvpSolverTelemetrySnapshot::Damped(snapshot)
            | BvpSolverTelemetrySnapshot::Frozen(snapshot) => snapshot
                .atom_lambdify
                .as_ref()
                .or(snapshot.legacy_lambdify.as_ref())
                .map(|callback| (callback.parallel_dispatches, callback.sequential_dispatches))
                .unwrap_or_default(),
            BvpSolverTelemetrySnapshot::Sci(snapshot) => {
                (snapshot.parallel_dispatches, snapshot.sequential_dispatches)
            }
        }
    }

    #[test]
    fn story_large_bvp_callback_policies_match_and_report_dispatches() {
        for family in [PolicyTestSolver::Damped, PolicyTestSolver::Sci] {
            let routes = [
                BvpLambdifyExecutionPolicy::Sequential,
                BvpLambdifyExecutionPolicy::Auto { min_work: 0 },
                BvpLambdifyExecutionPolicy::Parallel { min_work: 0 },
            ];
            let mut solutions = Vec::new();
            for policy in routes {
                let (values, telemetry) = solve_large_callback_policy_case(family, policy);
                let dispatches = callback_dispatch_counts(&telemetry);
                eprintln!(
                    "BVP callback policy family={family:?} policy={policy:?} parallel={} sequential={}",
                    dispatches.0, dispatches.1
                );
                if matches!(policy, BvpLambdifyExecutionPolicy::Parallel { .. }) {
                    assert!(dispatches.0 > 0, "forced parallel route did not dispatch");
                }
                solutions.push(values);
            }
            for candidate in &solutions[1..] {
                assert_profile_close(&solutions[0], candidate, 2e-5, "callback policy parity");
            }
        }
    }

    #[test]
    fn matrix_compact_hmx_sci_large_grid_remains_bounded() {
        let values = solve_compact_sci(
            BvpSciAssembly::AtomViewNative,
            BvpSciMatrixLayout::BandedAuto,
            96,
        );
        assert_eq!(values.len() % 6, 0);
        assert!(values.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn matrix_compact_hmx_sci_aot_matches_lambdify_across_frontends_and_layouts() {
        if !aot_tests_enabled("tcc") {
            return;
        }
        for assembly in [BvpSciAssembly::ExprLegacy, BvpSciAssembly::AtomViewNative] {
            for layout in [
                BvpSciMatrixLayout::Dense,
                BvpSciMatrixLayout::Sparse,
                BvpSciMatrixLayout::BandedAuto,
            ] {
                let lambdify = solve_compact_sci(assembly, layout, 32);
                let artifact_dir = tempfile::tempdir().expect("AOT artifact directory");
                let aot = solve_compact_sci_with_aot(
                    assembly,
                    layout,
                    32,
                    Some(
                        SymbolicIvpGeneratedBackendConfig::build_if_missing_release(
                            artifact_dir.path(),
                        )
                        .with_c_tcc(),
                    ),
                );
                assert_profile_close(
                    &lambdify,
                    &aot,
                    2e-5,
                    &format!("Sci AOT parity {assembly:?}/{layout:?}"),
                );
            }
        }
    }

    #[test]
    fn matrix_compact_hmx_damped_aot_matches_lambdify_across_frontends_and_matrices() {
        if !aot_tests_enabled("tcc") {
            return;
        }
        for matrix in [
            ReactorBvpMatrixBackend::Banded,
            ReactorBvpMatrixBackend::Sparse,
        ] {
            for frontend in [
                ReactorBvpSymbolicBackend::AtomView,
                ReactorBvpSymbolicBackend::ExprLegacy,
            ] {
                let lambdify = ReactorBvpSolverConfig::default_lambdify()
                    .with_matrix_backend(matrix)
                    .with_symbolic_backend(frontend);
                let aot = match matrix {
                    ReactorBvpMatrixBackend::Banded => {
                        ReactorBvpSolverConfig::banded_aot(ReactorBvpAotCompiler::CTcc)
                    }
                    ReactorBvpMatrixBackend::Sparse => {
                        ReactorBvpSolverConfig::sparse_aot(ReactorBvpAotCompiler::CTcc)
                    }
                }
                .with_symbolic_backend(frontend);
                assert_config_matches_reference(lambdify, aot, Some("tcc"), 2e-5);
            }
        }
    }

    /// AOT smoke tests are gated because toolchains like tcc are optional on user machines.
    #[test]
    fn matrix_compact_hmx_aot_tcc_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::default_lambdify(),
            ReactorBvpSolverConfig::banded_aot(ReactorBvpAotCompiler::CTcc),
            Some("tcc"),
            1e-6,
        );
    }

    /// GCC-backed AOT smoke test for the same compact combustion fixture.
    #[test]
    fn matrix_compact_hmx_aot_gcc_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::default_lambdify(),
            ReactorBvpSolverConfig::banded_aot(ReactorBvpAotCompiler::CGcc),
            Some("gcc"),
            1e-6,
        );
    }

    /// Zig-backed AOT smoke test for the same compact combustion fixture.
    #[test]
    fn matrix_compact_hmx_aot_zig_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::default_lambdify(),
            ReactorBvpSolverConfig::banded_aot(ReactorBvpAotCompiler::Zig),
            Some("zig"),
            1e-6,
        );
    }

    /// Sparse AOT should remain available as a compatibility matrix route when the toolchain exists.
    #[test]
    fn matrix_compact_hmx_sparse_aot_tcc_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::sparse_lambdify(),
            ReactorBvpSolverConfig::sparse_aot(ReactorBvpAotCompiler::CTcc),
            Some("tcc"),
            1e-6,
        );
    }

    /// Sparse AOT with GCC should keep the same compact-task smoke coverage.
    #[test]
    fn matrix_compact_hmx_sparse_aot_gcc_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::sparse_lambdify(),
            ReactorBvpSolverConfig::sparse_aot(ReactorBvpAotCompiler::CGcc),
            Some("gcc"),
            1e-6,
        );
    }

    /// Sparse AOT with Zig should keep the same compact-task smoke coverage.
    #[test]
    fn matrix_compact_hmx_sparse_aot_zig_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::sparse_lambdify(),
            ReactorBvpSolverConfig::sparse_aot(ReactorBvpAotCompiler::Zig),
            Some("zig"),
            1e-6,
        );
    }
}
