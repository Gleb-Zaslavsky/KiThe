#[cfg(test)]
mod tests {
    use super::super::SimpleReactorBVP::*;
    use crate::Kinetics::User_reactions::KinData;
    use crate::Kinetics::mechfinder_api::{ReactionData, ReactionType};
    use crate::Kinetics::stoichiometry_analyzer::StoichAnalyzer;
    use crate::ReactorsBVP::reactor_BVP_utils::{
        BoundsConfig, ScalingConfig, ToleranceConfig, create_bounds_map, create_tolerance_map,
    };
    use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{
        AdaptiveGridConfig, SolverParams,
    };
    use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;
    use approx::assert_relative_eq;
    use log::info;
    use nalgebra::{DMatrix, DVector};
    use std::collections::HashMap;
    use std::vec;

    fn create_test_reactor() -> SimpleReactorTask {
        let mut reactor = SimpleReactorTask::new();

        // Set basic properties
        reactor.P = 101325.0; // Pa
        reactor.Tm = 500.0; // K
        reactor.Lambda = 0.05; // W/m/K
        reactor.Cp = 1000.0; // J/kg/K
        reactor.m = 0.01; // kg/s

        // Set scaling
        reactor.scaling = ScalingConfig::new(100.0, 0.1, 100.0);

        // Set diffusion coefficients
        let mut diffusion = HashMap::new();
        diffusion.insert("A".to_string(), 1e-5);
        diffusion.insert("B".to_string(), 1.2e-5);
        diffusion.insert("C".to_string(), 1.2e-5);
        reactor.Diffusion = diffusion;

        // Set boundary conditions
        let mut bc = HashMap::new();
        bc.insert("A".to_string(), 0.5);
        bc.insert("B".to_string(), 0.3);
        bc.insert("C".to_string(), 0.2);
        bc.insert("T".to_string(), 450.0);
        reactor.boundary_condition = bc;
        let reactions = vec![
            FastElemReact {
                eq: "A=>B".to_string(),
                A: 1e10,
                n: 0.0,
                E: 50000.0,
                Q: -100000.0,
            },
            FastElemReact {
                eq: "B=>A+C".to_string(),
                A: 1e8,
                n: 0.5,
                E: 30000.0,
                Q: 50000.0,
            },
        ];

        let _ = reactor.fast_react_set(reactions);
        reactor
    }

    /// Create a small reactor snapshot that is ready for BVP equation assembly tests.
    fn create_policy_test_reactor(policy: SymbolicRhsAssemblyPolicy) -> SimpleReactorTask {
        let mut reactor = create_test_reactor().with_symbolic_rhs_assembly_policy(policy);
        reactor.kindata.substances = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        reactor.M = 1.0 / (0.5 / 28.0 + 0.3 / 44.0 + 0.2 / 20.0) / 1000.0;
        reactor
    }

    /// Create a compact HMX-like reactor snapshot that survives the full `setup_bvp` pipeline.
    fn create_setup_bvp_test_reactor() -> SimpleReactorTask {
        let eq = "HMX=>HMXprod".to_string();
        let c_p = 0.35 * 4.184 * 1000.0;
        let lambda_eff = 0.07;
        let n = 0.0;
        let m = 0.077 * (1e6_f64 / 1e5_f64).powf(0.748) / 1e2;
        let a = 1.3e5;
        let e = 5000.0 * 4.184;
        let t0 = 600.0;
        let t_scale = 600.0;
        let p = 1e6;
        let tm = 1500.0;
        let hmx = HashMap::from([
            ("H".to_string(), 4),
            ("N".to_string(), 8),
            ("C".to_string(), 8),
            ("O".to_string(), 8),
        ]);
        let hmxprod = HashMap::from([
            ("H".to_string(), 6),
            ("C".to_string(), 1),
            ("O".to_string(), 1),
        ]);
        let groups = Some(HashMap::from([
            ("HMX".to_string(), hmx.clone()),
            ("HMXprod".to_string(), hmxprod),
        ]));

        let mut reactor = SimpleReactorTask::new();
        let reactions = vec![FastElemReact {
            eq,
            A: a,
            n,
            E: e,
            Q: 3000.0 * 1e3 / 100.0,
        }];

        reactor
            .fast_react_set(reactions)
            .expect("valid elementary reaction should be accepted");
        reactor.kindata.substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        reactor.kindata.groups = groups;

        let diffusion = {
            let ro0 = 34.2e-3 * p / (R_G * t0);
            let d = lambda_eff / (c_p * ro0);
            HashMap::from([("HMX".to_string(), d), ("HMXprod".to_string(), d)])
        };
        let boundary_condition = HashMap::from([
            ("HMX".to_string(), 1.0 - 1e-3),
            ("HMXprod".to_string(), 1e-3),
            ("T".to_string(), t0),
        ]);
        let thermal_effects = vec![3000.0 * 1e3 / 100.0];
        let scaling = ScalingConfig::new(t_scale, 9e-4, t_scale);
        reactor.set_parameters(
            thermal_effects,
            p,
            tm,
            c_p,
            boundary_condition,
            lambda_eff,
            diffusion,
            m,
            scaling,
        );
        reactor.M = 34.2 / 1000.0;
        reactor
    }

    #[test]
    fn test_fast_react_set_valid_input() {
        let mut reactor = SimpleReactorTask::new();

        let reactions = vec![
            FastElemReact {
                eq: "A=>B".to_string(),
                A: 1e10,
                n: 0.0,
                E: 50000.0,
                Q: -100000.0,
            },
            FastElemReact {
                eq: "B=>A+C".to_string(),
                A: 1e8,
                n: 0.5,
                E: 30000.0,
                Q: 50000.0,
            },
        ];

        let result = reactor.fast_react_set(reactions);
        assert!(result.is_ok());
        assert_eq!(reactor.kindata.vec_of_equations.len(), 2);
        assert_eq!(reactor.thermal_effects.len(), 2);
        assert_eq!(reactor.thermal_effects[0], -100000.0);
        assert_eq!(reactor.thermal_effects[1], 50000.0);
    }

    #[test]
    fn test_fast_react_set_empty_equation() {
        let mut reactor = SimpleReactorTask::new();

        let reactions = vec![FastElemReact {
            eq: "".to_string(),
            A: 1e10,
            n: 0.0,
            E: 50000.0,
            Q: -100000.0,
        }];

        let result = reactor.fast_react_set(reactions);
        assert!(result.is_err());
        match result {
            Err(ReactorError::MissingData(msg)) => {
                assert!(msg.contains("No equation in input hashmap"));
            }
            _ => panic!("Expected MissingData error"),
        }
    }

    #[test]
    fn test_fast_react_set_nan_parameters() {
        let mut reactor = SimpleReactorTask::new();
        let _ = reactor.kinetic_processing();
        let reactions = vec![FastElemReact {
            eq: "A=>B".to_string(),
            A: f64::NAN,
            n: 0.0,
            E: 50000.0,
            Q: -100000.0,
        }];

        let result = reactor.fast_react_set(reactions);
        assert!(result.is_err());
        match result {
            Err(ReactorError::MissingData(msg)) => {
                assert!(msg.contains("Missing Arrhenius parameter 'A'"));
            }
            _ => panic!("Expected MissingData error"),
        }
    }

    #[test]
    fn test_mean_molar_mass_calculation() {
        let mut reactor = create_test_reactor();

        // Setup kinetic data with substances and molar masses
        let mut kindata = KinData::new();
        kindata.substances = vec!["A".to_string(), "B".to_string(), "C".to_string()];

        let mut stoich_analyzer = StoichAnalyzer::new();
        stoich_analyzer.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]); // g/mol
        kindata.stecheodata = stoich_analyzer;

        reactor.kindata = kindata;

        let result = reactor.mean_molar_mass();
        assert!(result.is_ok());
        let res = 1.0 / (0.5 / 28.0 + 0.3 / 44.0 + 0.2 / 20.0);
        let res = res / 1000.0;
        // Expected: 0.5 * 28.0 + 0.3 * 44.0 = 14.0 + 13.2 = 27.2
        assert!((reactor.M - res).abs() < 1e-10);
    }

    #[test]
    fn test_mean_molar_mass_missing_data() {
        let mut reactor = create_test_reactor();

        let result = reactor.mean_molar_mass();
        assert!(result.is_err());
        match result {
            Err(ReactorError::MissingData(msg)) => {
                assert!(msg.contains("Molar masses not calculated"));
            }
            _ => panic!("Expected MissingData error"),
        }
    }

    #[test]
    fn test_transport_coefficients() {
        let mut reactor = create_test_reactor();
        let _ = reactor.kinetic_processing();
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        let _ = reactor.mean_molar_mass();
        let _ = reactor.scaling_processing();

        reactor
            .transport_coefficients()
            .expect("Failed to calculate transport coefficients");
        let coeffs = reactor.D_ro_map.clone();
        // Check that coefficients are calculated for substances in diffusion map
        assert!(coeffs.contains_key("A"));
        assert!(coeffs.contains_key("B"));

        // Verify calculation: D*ro = D0*ro0*(T/T0)^0.5
        let M = 1.0 / (0.5 / 28.0 + 0.3 / 44.0 + 0.2 / 20.0);
        let M = M / 1000.0;
        let ro0 = M * reactor.P / (R_G * 298.15);
        let temp_factor = (reactor.Tm / 298.15).powf(0.5);

        let expected_a = 1e-5 * ro0 * temp_factor;
        let expected_b = 1.2e-5 * ro0 * temp_factor;

        assert!((coeffs["A"] - expected_a).abs() < 1e-10);
        assert!((coeffs["B"] - expected_b).abs() < 1e-10);
    }

    #[test]
    fn test_scaling_processing() {
        let mut reactor = create_test_reactor();

        let result = reactor.scaling_processing();
        assert!(result.is_ok());

        // Check that L is set correctly
        assert_eq!(reactor.L, 0.1);
    }

    #[test]
    fn test_scaling_processing_invalid_dt() {
        let mut reactor = SimpleReactorTask::new();
        // Set invalid dT in scaling
        reactor.scaling = ScalingConfig::new(-100.0, 0.1, 100.0); // Negative dT

        let result = reactor.scaling_processing();
        assert!(result.is_err());
        match result {
            Err(ReactorError::InvalidConfiguration(msg)) => {
                assert!(msg.contains("Temperature scaling dT must be positive"));
            }
            _ => panic!("Expected InvalidConfiguration error"),
        }
    }

    #[test]
    fn test_peclet_numbers() {
        let mut reactor = create_test_reactor();
        let _ = reactor.scaling_processing();
        let _ = reactor.kinetic_processing();
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        let _ = reactor.mean_molar_mass();
        assert_eq!(
            reactor.kindata.substances,
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );

        let result = reactor.peclet_numbers();
        assert!(result.is_ok());
        let cached_transport_snapshot = reactor.D_ro_map.clone();
        reactor
            .transport_coefficients()
            .expect("cached transport coefficients should be refreshed in place");
        assert_eq!(cached_transport_snapshot, reactor.D_ro_map);

        // Pe_q = (L * m * Cp) / Lambda
        let expected_pe_q = (0.1 * 0.01 * 1000.0) / 0.05;
        assert!((reactor.Pe_q - expected_pe_q).abs() < 1e-10);

        // Check Pe_D vector length
        assert_eq!(reactor.Pe_D.len(), 3);

        // Verify Pe_D calculations
        let transport_coeffs = reactor.D_ro_map.clone();
        let expected_pe_d_a = (reactor.m * reactor.L) / transport_coeffs["A"];
        let expected_pe_d_b = (reactor.m * reactor.L) / transport_coeffs["B"];

        assert!((reactor.Pe_D[0] - expected_pe_d_a).abs() < 1e-10);
        assert!((reactor.Pe_D[1] - expected_pe_d_b).abs() < 1e-10);
    }

    #[test]
    fn test_peclet_numbers_invalid_lambda() {
        let mut reactor = create_test_reactor();
        reactor.Lambda = 0.0; // Invalid

        let result = reactor.peclet_numbers();
        assert!(result.is_err());
        match result {
            Err(ReactorError::InvalidConfiguration(msg)) => {
                assert!(msg.contains("Lambda must be positive"));
            }
            _ => panic!("Expected InvalidConfiguration error"),
        }
    }

    #[test]
    fn test_peclet_numbers_invalid_mass_flow() {
        let mut reactor = create_test_reactor();
        reactor.m = -0.01; // Invalid

        let result = reactor.peclet_numbers();
        assert!(result.is_err());
        match result {
            Err(ReactorError::InvalidConfiguration(msg)) => {
                assert!(msg.contains("Mass flow rate must be positive"));
            }
            _ => panic!("Expected InvalidConfiguration error"),
        }
    }

    #[test]
    fn test_kindata_symbolic_rhs_policy_builder_sets_policy() {
        let reactor = create_test_reactor()
            .with_symbolic_rhs_assembly_policy(SymbolicRhsAssemblyPolicy::Sequential);
        assert_eq!(
            reactor.symbolic_rhs_assembly_policy,
            SymbolicRhsAssemblyPolicy::Sequential
        );
    }

    #[test]
    fn test_backend_boundary_conditions_preserve_single_value_shape() {
        let mut solver = BVPSolver::default();
        solver
            .BorderConditions
            .insert("Teta".to_string(), (0, 0.25));
        solver.BorderConditions.insert("C0".to_string(), (0, 0.5));
        solver.BorderConditions.insert("J0".to_string(), (1, 1e-10));

        let backend_bc = solver.backend_boundary_conditions();

        assert_eq!(backend_bc.len(), 3);
        assert_eq!(backend_bc.get("Teta"), Some(&vec![(0, 0.25)]));
        assert_eq!(backend_bc.get("C0"), Some(&vec![(0, 0.5)]));
        assert_eq!(backend_bc.get("J0"), Some(&vec![(1, 1e-10)]));
    }

    #[test]
    fn test_build_nrbvp_backend_preserves_solver_snapshot_and_settings() {
        let mut reactor = create_setup_bvp_test_reactor();
        reactor
            .setup_bvp()
            .expect("setup_bvp should succeed before building the backend snapshot");

        let unknown_count = reactor.solver.unknowns.len();
        let initial_guess = DMatrix::from_element(unknown_count, 4, 0.5);
        let rel_tolerance = Some(HashMap::from([
            ("Teta".to_string(), 1e-6),
            ("q".to_string(), 1e-6),
            ("C0".to_string(), 1e-6),
            ("J0".to_string(), 1e-6),
            ("C1".to_string(), 1e-6),
            ("J1".to_string(), 1e-6),
        ]));
        let bounds = Some(HashMap::from([
            ("Teta".to_string(), (-100.0, 100.0)),
            ("q".to_string(), (-1e20, 1e20)),
            ("C0".to_string(), (0.0, 1.0)),
            ("J0".to_string(), (-1e20, 1e20)),
            ("C1".to_string(), (0.0, 1.0)),
            ("J1".to_string(), (-1e20, 1e20)),
        ]));

        let backend = reactor
            .solver
            .build_nrbvp_backend(NrbvpHandoffConfig::new(
                initial_guess.clone(),
                0.0,
                1.0,
                4,
                "forward".to_string(),
                "Damped".to_string(),
                None,
                None,
                "Sparse".to_string(),
                1e-8,
                rel_tolerance.clone(),
                100,
                bounds.clone(),
                Some("info".to_string()),
                false,
            ))
            .expect("default lambdify backend should build NRBVP options");

        assert_eq!(backend.eq_system, reactor.solver.eq_system);
        assert_eq!(backend.values, reactor.solver.unknowns);
        assert_eq!(backend.arg, reactor.solver.arg_name);
        assert_eq!(
            backend.BorderConditions,
            reactor.solver.backend_boundary_conditions()
        );
        assert_eq!(backend.initial_guess, initial_guess);
        assert_eq!(backend.t0, 0.0);
        assert_eq!(backend.t_end, 1.0);
        assert_eq!(backend.rel_tolerance, rel_tolerance);
        assert_eq!(backend.Bounds, bounds);
        assert_eq!(backend.scheme, "forward");
        assert_eq!(backend.strategy, "Damped");
        assert_eq!(backend.method, "Sparse");
        assert_eq!(backend.max_iterations, 100);
        assert_eq!(backend.abs_tolerance, 1e-8);
        assert_eq!(backend.loglevel, Some("info".to_string()));
        assert_eq!(backend.n_steps, 4);
        assert_eq!(
            backend.generated_backend_config().symbolic_assembly_backend,
            RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            backend.generated_backend_config().matrix_backend_override,
            Some(RustedSciThe::symbolic::codegen::codegen_provider_api::MatrixBackend::Banded)
        );
        assert!(backend.x_mesh.len() > 0);
    }

    #[test]
    fn test_reactor_bvp_solver_config_default_is_banded_atomview_lambdify() {
        let options =
            crate::ReactorsBVP::solver_backend::ReactorBvpSolverConfig::default_lambdify()
                .to_rusted_options()
                .expect("default lambdify config should be supported");

        assert_eq!(
            options.generated_backend_config.symbolic_assembly_backend,
            RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            options.generated_backend_config.matrix_backend_override,
            Some(RustedSciThe::symbolic::codegen::codegen_provider_api::MatrixBackend::Banded)
        );
        assert_eq!(options.method, "Sparse");
    }

    #[test]
    fn test_reactor_bvp_solver_config_sparse_lambdify_is_available() {
        let options = crate::ReactorsBVP::solver_backend::ReactorBvpSolverConfig::sparse_lambdify()
            .to_rusted_options()
            .expect("sparse lambdify config should be supported");

        assert_eq!(
            options.generated_backend_config.symbolic_assembly_backend,
            RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            options.generated_backend_config.matrix_backend_override,
            None
        );
    }

    #[test]
    fn test_reactor_bvp_solver_config_banded_aot_tcc_builds_options() {
        let config =
            crate::ReactorsBVP::solver_backend::ReactorBvpSolverConfig::from_generated_backend_name(
                "banded_aot_tcc",
            )
            .expect("AOT aliases should parse");

        let options = config
            .to_rusted_options()
            .expect("AOT config should build RustedSciThe options");
        assert_eq!(
            options.generated_backend_config.symbolic_assembly_backend,
            RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            options.generated_backend_config.matrix_backend_override,
            Some(RustedSciThe::symbolic::codegen::codegen_provider_api::MatrixBackend::Banded)
        );
        assert_eq!(
            options.generated_backend_config.aot_c_compiler,
            Some("tcc".to_string())
        );
        assert!(matches!(
            options.generated_backend_config.aot_build_policy,
            RustedSciThe::numerical::BVP_Damp::generated_solver_handoff::AotBuildPolicy::BuildIfMissing {
                ..
            }
        ));
    }

    #[test]
    fn test_create_bvp_equations_parallel_matches_sequential() {
        let mut sequential = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        let mut parallel = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Parallel);

        sequential
            .scaling_processing()
            .expect("sequential scaling should succeed");
        sequential
            .kinetic_processing()
            .expect("sequential kinetics should succeed");
        sequential.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        sequential
            .mean_molar_mass()
            .expect("sequential mean molar mass should succeed");
        sequential
            .peclet_numbers()
            .expect("sequential Peclet numbers should succeed");
        sequential
            .create_bvp_equations()
            .expect("sequential equation assembly should succeed");
        sequential
            .set_solver_BC()
            .expect("sequential boundary conditions should succeed");
        sequential
            .check_before_solution()
            .expect("sequential invariant check should succeed");

        parallel
            .scaling_processing()
            .expect("parallel scaling should succeed");
        parallel
            .kinetic_processing()
            .expect("parallel kinetics should succeed");
        parallel.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        parallel
            .mean_molar_mass()
            .expect("parallel mean molar mass should succeed");
        parallel
            .peclet_numbers()
            .expect("parallel Peclet numbers should succeed");
        parallel
            .create_bvp_equations()
            .expect("parallel equation assembly should succeed");
        parallel
            .set_solver_BC()
            .expect("parallel boundary conditions should succeed");
        parallel
            .check_before_solution()
            .expect("parallel invariant check should succeed");

        let seq_eqs: Vec<String> = sequential
            .solver
            .eq_system
            .iter()
            .map(|eq| format!("{}", eq))
            .collect();
        let par_eqs: Vec<String> = parallel
            .solver
            .eq_system
            .iter()
            .map(|eq| format!("{}", eq))
            .collect();
        assert_eq!(seq_eqs, par_eqs);
        assert_eq!(sequential.solver.unknowns, parallel.solver.unknowns);
        assert_eq!(
            format!("{}", sequential.heat_release),
            format!("{}", parallel.heat_release)
        );

        let mut seq_map: Vec<(String, String, String)> = sequential
            .map_of_equations
            .iter()
            .map(|(key, (unknown, equation))| {
                (key.clone(), unknown.clone(), format!("{}", equation))
            })
            .collect();
        let mut par_map: Vec<(String, String, String)> = parallel
            .map_of_equations
            .iter()
            .map(|(key, (unknown, equation))| {
                (key.clone(), unknown.clone(), format!("{}", equation))
            })
            .collect();
        seq_map.sort();
        par_map.sort();
        assert_eq!(seq_map, par_map);

        let seq_rates: Vec<String> = sequential
            .kindata
            .vec_of_equations
            .iter()
            .map(|eq| format!("{}", sequential.map_eq_rate.get(eq).unwrap()))
            .collect();
        let par_rates: Vec<String> = parallel
            .kindata
            .vec_of_equations
            .iter()
            .map(|eq| format!("{}", parallel.map_eq_rate.get(eq).unwrap()))
            .collect();
        assert_eq!(seq_rates, par_rates);
    }

    #[test]
    fn test_check_task_accepts_normalized_boundary_fractions() {
        // Use a fully normalized snapshot so this test validates boundary semantics
        // instead of the default placeholder state used by the lower-level helper.
        let reactor = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        assert!(reactor.check_task().is_ok());
    }

    #[test]
    fn test_check_task_rejects_missing_boundary_species() {
        // Start from a valid reactor state so the missing-species branch is the one
        // that gets exercised and reported to the caller.
        let mut reactor = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        reactor.boundary_condition.remove("B");

        let result = reactor.check_task();
        assert!(matches!(result, Err(ReactorError::MissingData(_))));
    }

    #[test]
    fn test_check_task_rejects_zero_sum_boundary_fractions() {
        let mut reactor = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        reactor.boundary_condition.insert("A".to_string(), 0.0);
        reactor.boundary_condition.insert("B".to_string(), 0.0);
        reactor.boundary_condition.insert("C".to_string(), 0.0);

        let result = reactor.check_task();
        assert!(matches!(result, Err(ReactorError::InvalidNumericValue(_))));
    }

    #[test]
    fn test_check_task_rejects_negative_boundary_fraction() {
        let mut reactor = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        reactor.boundary_condition.insert("A".to_string(), -0.1);
        reactor.boundary_condition.insert("B".to_string(), 0.7);
        reactor.boundary_condition.insert("C".to_string(), 0.4);

        let result = reactor.check_task();
        assert!(matches!(result, Err(ReactorError::InvalidNumericValue(_))));
    }

    #[test]
    fn test_setup_bvp_applies_optional_flux_defaults() {
        let mut reactor = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);

        // Assemble only the pieces required by `set_solver_BC()` so this test stays
        // focused on default boundary handling instead of the full BVP pipeline.
        reactor
            .scaling_processing()
            .expect("scaling should succeed");
        reactor
            .kinetic_processing()
            .expect("kinetics should succeed");
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        reactor.M = 1.0 / (0.5 / 28.0 + 0.3 / 44.0 + 0.2 / 20.0) / 1000.0;
        reactor
            .peclet_numbers()
            .expect("Peclet numbers should succeed");
        reactor
            .create_bvp_equations()
            .expect("equations should be assembled");
        reactor
            .set_solver_BC()
            .expect("solver boundary conditions should succeed");

        assert_eq!(reactor.solver.BorderConditions.get("q"), Some(&(1, 1e-10)));
        assert_eq!(reactor.solver.BorderConditions.get("J0"), Some(&(1, 1e-10)));
    }

    #[test]
    fn test_set_solver_bc_rejects_missing_canonical_equation_contract() {
        let mut reactor = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        reactor
            .scaling_processing()
            .expect("scaling should succeed");
        reactor
            .kinetic_processing()
            .expect("kinetics should succeed");
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        reactor
            .mean_molar_mass()
            .expect("mean molar mass should succeed");
        reactor
            .peclet_numbers()
            .expect("Peclet numbers should succeed");
        reactor
            .create_bvp_equations()
            .expect("equations should be assembled");

        reactor.map_of_equations.clear();
        let result = reactor.set_solver_BC();

        assert!(matches!(result, Err(ReactorError::InvalidConfiguration(_))));
    }

    #[test]
    fn test_create_bvp_equations_rejects_missing_transport_snapshot() {
        let mut reactor = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        reactor
            .scaling_processing()
            .expect("scaling should succeed");
        reactor
            .kinetic_processing()
            .expect("kinetics should succeed");
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        reactor
            .mean_molar_mass()
            .expect("mean molar mass should succeed");
        reactor
            .peclet_numbers()
            .expect("Peclet numbers should succeed");
        reactor.D_ro_map.remove("B");

        let result = reactor.create_bvp_equations();
        assert!(matches!(result, Err(ReactorError::MissingData(_))));
    }

    #[test]
    fn test_create_bvp_equations_uses_substance_order_for_transport_coefficients() {
        let mut baseline = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        baseline
            .scaling_processing()
            .expect("scaling should succeed");
        baseline
            .kinetic_processing()
            .expect("kinetics should succeed");
        baseline.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        baseline
            .mean_molar_mass()
            .expect("mean molar mass should succeed");
        baseline
            .peclet_numbers()
            .expect("Peclet numbers should succeed");
        baseline
            .create_bvp_equations()
            .expect("baseline equation assembly should succeed");
        baseline
            .set_solver_BC()
            .expect("baseline boundary conditions should succeed");

        let mut shuffled = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        shuffled
            .scaling_processing()
            .expect("scaling should succeed");
        shuffled
            .kinetic_processing()
            .expect("kinetics should succeed");
        shuffled.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        shuffled
            .mean_molar_mass()
            .expect("mean molar mass should succeed");
        shuffled
            .peclet_numbers()
            .expect("Peclet numbers should succeed");
        shuffled.D_ro_map = HashMap::from([
            ("B".to_string(), shuffled.D_ro_map["B"]),
            ("C".to_string(), shuffled.D_ro_map["C"]),
            ("A".to_string(), shuffled.D_ro_map["A"]),
        ]);
        shuffled
            .create_bvp_equations()
            .expect("shuffled equation assembly should succeed");
        shuffled
            .set_solver_BC()
            .expect("shuffled boundary conditions should succeed");

        assert_eq!(baseline.solver.unknowns, shuffled.solver.unknowns);
        assert_eq!(
            baseline
                .solver
                .eq_system
                .iter()
                .map(|eq| format!("{}", eq))
                .collect::<Vec<_>>(),
            shuffled
                .solver
                .eq_system
                .iter()
                .map(|eq| format!("{}", eq))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            baseline.solver.BorderConditions,
            shuffled.solver.BorderConditions
        );
    }

    #[test]
    fn test_task_report_uses_canonical_order() {
        let mut reactor = create_policy_test_reactor(SymbolicRhsAssemblyPolicy::Sequential);
        reactor
            .scaling_processing()
            .expect("scaling should succeed");
        reactor
            .kinetic_processing()
            .expect("kinetics should succeed");
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        reactor
            .mean_molar_mass()
            .expect("mean molar mass should succeed");
        reactor
            .peclet_numbers()
            .expect("Peclet numbers should succeed");

        let report = reactor.task_report();

        assert_eq!(
            report
                .boundary_conditions
                .iter()
                .map(|(name, _)| name.clone())
                .collect::<Vec<_>>(),
            vec![
                "T".to_string(),
                "A".to_string(),
                "B".to_string(),
                "C".to_string(),
            ]
        );
        assert_eq!(
            report
                .diffusion_coefficients
                .iter()
                .map(|(name, _)| name.clone())
                .collect::<Vec<_>>(),
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
        assert_eq!(
            report
                .transport_coefficients
                .iter()
                .map(|(name, _)| name.clone())
                .collect::<Vec<_>>(),
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
        assert_eq!(
            report
                .peclet_numbers
                .iter()
                .map(|(name, _)| name.clone())
                .collect::<Vec<_>>(),
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
        assert_eq!(report.substances, vec!["A", "B", "C"]);
        assert_eq!(report.reactions.len(), 2);
    }

    #[test]
    fn test_equation_report_rows_follow_solver_order() {
        let mut reactor = create_setup_bvp_test_reactor();
        reactor
            .setup_bvp()
            .expect("setup_bvp should succeed for equation report test");

        let rows = reactor.equation_report_rows();
        let unknowns = reactor.solver.unknowns.clone();

        assert_eq!(rows.len(), unknowns.len());
        assert_eq!(
            rows.iter()
                .map(|(_, unknown, _)| unknown.clone())
                .collect::<Vec<_>>(),
            unknowns
        );
    }

    #[test]
    fn test_balance_report_reflects_cached_quality() {
        let mut reactor = SimpleReactorTask::new();
        reactor.solver.quality.energy_balane_error_abs = 1.25;
        reactor.solver.quality.energy_balane_error_rel = 4.5;
        reactor.solver.quality.sum_of_mass_fractions = vec![(1, 0.97)];
        reactor.solver.quality.atomic_mass_balance_error = vec![(2, 0.02)];

        let report = reactor.balance_report();

        assert_eq!(report.energy_balane_error_abs, 1.25);
        assert_eq!(report.energy_balane_error_rel, 4.5);
        assert_eq!(report.sum_of_mass_fractions, vec![(1, 0.97)]);
        assert_eq!(report.atomic_mass_balance_error, vec![(2, 0.02)]);
    }

    #[test]
    fn test_solution_render_data_returns_owned_snapshot() {
        let mut reactor = create_setup_bvp_test_reactor();
        reactor
            .setup_bvp()
            .expect("setup_bvp should succeed before snapshot extraction");
        reactor.solver.x_mesh = Some(DVector::from_vec(vec![0.0, 0.5, 1.0]));
        reactor.solver.solution = Some(DMatrix::from_row_slice(
            3,
            reactor.solver.unknowns.len(),
            &[
                1.0, 0.2, 0.9, 0.1, 0.8, 0.05, //
                1.0, 0.3, 0.7, 0.2, 0.6, 0.08, //
                1.0, 0.4, 0.5, 0.3, 0.4, 0.1,
            ],
        ));

        let snapshot = reactor
            .solution_render_data()
            .expect("solution snapshot should be available once solution data is present");

        assert_eq!(snapshot.unknowns, reactor.solver.unknowns);
        assert_eq!(snapshot.arg_name, reactor.solver.arg_name);
        assert_eq!(
            snapshot.x_mesh.len(),
            reactor.solver.x_mesh.as_ref().unwrap().len()
        );
        assert_eq!(
            snapshot.solution.ncols(),
            reactor.solver.solution.as_ref().unwrap().ncols()
        );
        assert_eq!(
            snapshot.solution.nrows(),
            reactor.solver.solution.as_ref().unwrap().nrows()
        );
        assert_eq!(snapshot.x_mesh.as_slice(), &[0.0, 0.5, 1.0]);
    }

    #[test]
    fn test_solution_render_data_rejects_missing_solution() {
        let reactor = SimpleReactorTask::new();

        let result = reactor.solution_render_data();

        assert!(matches!(result, Err(ReactorError::MissingData(_))));
    }

    #[test]
    fn test_estimate_values_report_computes_single_reaction_temperature() {
        let mut reactor = SimpleReactorTask::new();
        reactor.kindata.vec_of_equations = vec!["A=>B".to_string()];
        reactor.thermal_effects = vec![1000.0];
        reactor.Cp = 500.0;
        reactor.boundary_condition.insert("T".to_string(), 300.0);

        let report = reactor
            .estimate_values_report()
            .expect("quick estimate should succeed for a valid single-reaction setup");

        assert_eq!(report.reaction_count, 1);
        assert_eq!(report.single_reaction_adiabatic_temperature, Some(302.0));
    }

    #[test]
    fn test_estimate_values_report_rejects_non_positive_heat_capacity() {
        let mut reactor = SimpleReactorTask::new();
        reactor.kindata.vec_of_equations = vec!["A=>B".to_string()];
        reactor.thermal_effects = vec![1000.0];
        reactor.Cp = 0.0;
        reactor.boundary_condition.insert("T".to_string(), 300.0);

        let result = reactor.estimate_values_report();

        assert!(matches!(result, Err(ReactorError::InvalidNumericValue(_))));
    }

    #[test]
    fn test_postprocessing_report_scales_solution_without_mutating_state() {
        let mut reactor = SimpleReactorTask::new();
        reactor.L = 2.0;
        // ScalingConfig::new(dT, L, T_scale).
        reactor.scaling = ScalingConfig::new(10.0, 2.0, 100.0);
        reactor.solver.unknowns = vec![
            "Teta".to_string(),
            "q".to_string(),
            "J0".to_string(),
            "C0".to_string(),
        ];
        reactor.solver.x_mesh = Some(DVector::from_vec(vec![0.0, 0.5]));
        reactor.solver.solution = Some(DMatrix::from_row_slice(
            2,
            4,
            &[
                0.0, 1.0, 2.0, 10.0, //
                1.0, 3.0, 4.0, 20.0,
            ],
        ));

        let report = reactor
            .postprocessing_report()
            .expect("postprocessing snapshot should be available for a complete state");

        assert_eq!(report.x_mesh.as_slice(), &[0.0, 1.0]);
        assert_eq!(report.solution[(0, 0)], 10.0);
        assert_eq!(report.solution[(1, 0)], 110.0);
        assert_eq!(report.solution[(0, 1)], 50.0);
        assert_eq!(report.solution[(1, 1)], 150.0);
        assert_eq!(report.solution[(0, 2)], 1.0);
        assert_eq!(report.solution[(1, 2)], 2.0);
        assert_eq!(report.solution[(0, 3)], 10.0);
        assert_eq!(report.solution[(1, 3)], 20.0);

        // The pure helper must not mutate the solver snapshot.
        assert_eq!(
            reactor.solver.x_mesh.as_ref().unwrap().as_slice(),
            &[0.0, 0.5]
        );
        assert_eq!(reactor.solver.solution.as_ref().unwrap()[(1, 0)], 1.0);
    }

    #[test]
    fn test_postprocessing_mutates_solver_state_from_snapshot() {
        let mut reactor = SimpleReactorTask::new();
        reactor.L = 2.0;
        // ScalingConfig::new(dT, L, T_scale).
        reactor.scaling = ScalingConfig::new(10.0, 2.0, 100.0);
        reactor.solver.unknowns = vec![
            "Teta".to_string(),
            "q".to_string(),
            "J0".to_string(),
            "C0".to_string(),
        ];
        reactor.solver.x_mesh = Some(DVector::from_vec(vec![0.0, 0.5]));
        reactor.solver.solution = Some(DMatrix::from_row_slice(
            2,
            4,
            &[
                0.0, 1.0, 2.0, 10.0, //
                1.0, 3.0, 4.0, 20.0,
            ],
        ));

        reactor
            .postprocessing()
            .expect("postprocessing should scale the stored solver snapshot");

        assert_eq!(
            reactor.solver.x_mesh.as_ref().unwrap().as_slice(),
            &[0.0, 1.0]
        );
        assert_eq!(reactor.solver.solution.as_ref().unwrap()[(0, 0)], 10.0);
        assert_eq!(reactor.solver.solution.as_ref().unwrap()[(1, 0)], 110.0);
        assert_eq!(reactor.solver.solution.as_ref().unwrap()[(0, 1)], 50.0);
        assert_eq!(reactor.solver.solution.as_ref().unwrap()[(1, 1)], 150.0);
    }

    #[test]
    fn test_postprocessing_report_rejects_missing_solution() {
        let reactor = SimpleReactorTask::new();
        let result = reactor.postprocessing_report();
        assert!(matches!(result, Err(ReactorError::MissingData(_))));
    }

    #[test]
    fn test_setup_bvp_builds_consistent_solver_contract() {
        let mut reactor = create_setup_bvp_test_reactor();

        let result = reactor.setup_bvp();
        assert!(
            result.is_ok(),
            "setup_bvp should succeed for a valid reactor snapshot: {result:?}"
        );

        let expected_unknowns = vec![
            "Teta".to_string(),
            "q".to_string(),
            "C0".to_string(),
            "J0".to_string(),
            "C1".to_string(),
            "J1".to_string(),
        ];

        assert_eq!(reactor.solver.unknowns, expected_unknowns);
        assert_eq!(reactor.solver.eq_system.len(), 6);
        assert_eq!(reactor.solver.BorderConditions.len(), 6);
        assert_eq!(reactor.map_of_equations.len(), 6);
        assert_eq!(reactor.Pe_D.len(), 2);

        assert!(reactor.M.is_finite() && reactor.M > 0.0);
        assert!(reactor.Pe_q.is_finite() && reactor.Pe_q > 0.0);
        assert!(
            reactor
                .Pe_D
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
        );

        assert_eq!(reactor.solver.BorderConditions.get("Teta"), Some(&(0, 0.0)));
        assert_eq!(reactor.solver.BorderConditions.get("q"), Some(&(1, 1e-10)));
        assert_eq!(reactor.solver.BorderConditions.get("C0"), Some(&(0, 0.999)));
        assert_eq!(reactor.solver.BorderConditions.get("J0"), Some(&(1, 1e-10)));
        assert_eq!(reactor.solver.BorderConditions.get("C1"), Some(&(0, 1e-3)));
        assert_eq!(reactor.solver.BorderConditions.get("J1"), Some(&(1, 1e-10)));

        for key in ["HMX", "HMXprod", "HMX_flux", "HMXprod_flux"] {
            assert!(
                reactor.map_of_equations.contains_key(key),
                "map_of_equations should contain {key}"
            );
        }
        assert!(reactor.heat_release.to_string().len() > 0);
    }

    #[test]
    fn test_check_balances_succeeds_on_constant_solution_snapshot() {
        let mut reactor = create_setup_bvp_test_reactor();
        reactor
            .setup_bvp()
            .expect("setup_bvp should succeed before balance checks");

        // Build a compact, fully finite solution snapshot so the balance
        // checks can run without a heavy numerical solve.
        let unknowns = reactor.solver.unknowns.clone();
        let mut solution = DMatrix::zeros(3, unknowns.len());
        for (col, name) in unknowns.iter().enumerate() {
            let value = match name.as_str() {
                "Teta" => 0.1,
                "q" => 1e-6,
                value if value.starts_with('C') && value.ends_with('0') => 0.999,
                value if value.starts_with('C') => 0.001,
                value if value.starts_with('J') => 0.0,
                _ => 0.0,
            };
            for row in 0..solution.nrows() {
                solution[(row, col)] = value;
            }
        }
        reactor.solver.solution = Some(solution);
        reactor.solver.x_mesh = Some(DVector::from_vec(vec![0.0, 0.5, 1.0]));

        reactor
            .check_balances()
            .expect("balance checks should accept a finite constant snapshot");

        assert!(reactor.solver.quality.energy_balane_error_abs.is_finite());
        assert!(reactor.solver.quality.energy_balane_error_rel.is_finite());
        assert!(reactor.solver.quality.sum_of_mass_fractions.is_empty());
        assert!(reactor.solver.quality.atomic_mass_balance_error.is_empty());
    }

    #[test]
    fn test_ideal_gas_density() {
        let mut reactor = create_test_reactor();
        let _ = reactor.kinetic_processing();
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        let _ = reactor.mean_molar_mass();

        let M = 1.0 / (0.5 / 28.0 + 0.3 / 44.0 + 0.2 / 20.0);
        let M = M / 1000.0;
        let M_ = reactor.M;
        assert_relative_eq!(M, M_, epsilon = 1e-6);
        let density = reactor
            .ideal_gas_density()
            .expect("Density should be computable");
        info!("Density: {}, M= {}", density, &reactor.M);
        let expected = M * reactor.P / (R_G * reactor.Tm);

        assert!((density - expected).abs() < 1e-10);
    }

    #[test]
    fn test_ideal_gas_density_rejects_invalid_state() {
        let mut reactor = create_test_reactor();
        reactor.M = 0.0;

        let result = reactor.ideal_gas_density();
        assert!(result.is_err());
        match result {
            Err(ReactorError::InvalidNumericValue(msg)) => {
                assert!(msg.contains("Mean molar mass must be finite and positive"));
            }
            _ => panic!("Expected InvalidNumericValue error"),
        }
    }

    #[test]
    fn test_validate_bvp_solution_matrix_rejects_empty_or_non_finite() {
        let empty = DMatrix::zeros(0, 0);
        let empty_result = super::super::SimpleReactorBVP::validate_bvp_solution_matrix(&empty);
        assert!(empty_result.is_err());
        match empty_result {
            Err(ReactorError::CalculationError(msg)) => {
                assert!(msg.contains("empty solution matrix"));
            }
            _ => panic!("Expected CalculationError for empty matrix"),
        }

        let nan_matrix = DMatrix::from_vec(1, 2, vec![1.0, f64::NAN]);
        let nan_result = super::super::SimpleReactorBVP::validate_bvp_solution_matrix(&nan_matrix);
        assert!(nan_result.is_err());
        match nan_result {
            Err(ReactorError::InvalidNumericValue(msg)) => {
                assert!(msg.contains("non-finite value"));
            }
            _ => panic!("Expected InvalidNumericValue for non-finite matrix"),
        }
    }

    #[test]
    fn test_reactor_error_display() {
        let errors = vec![
            ReactorError::MissingData("test data".to_string()),
            ReactorError::InvalidConfiguration("test config".to_string()),
            ReactorError::InvalidNumericValue("test numeric".to_string()),
            ReactorError::CalculationError("test calc".to_string()),
            ReactorError::ParseError("test parse".to_string()),
            ReactorError::IndexOutOfBounds("test index".to_string()),
        ];

        for error in errors {
            let error_string = format!("{}", error);
            assert!(!error_string.is_empty());
        }
    }

    #[test]
    fn test_read_only_accessors_expose_canonical_reactor_state() {
        let reactor = create_test_reactor();

        assert_eq!(reactor.diffusion_coefficients().get("A"), Some(&1e-5));
        assert_eq!(reactor.diffusion_coefficients().get("B"), Some(&1.2e-5));
        assert_eq!(reactor.boundary_conditions().get("T"), Some(&450.0));
        assert_eq!(reactor.transport_cache().get("A"), None);
        assert!(reactor.mass_peclet_numbers().is_empty());
        assert_eq!(reactor.thermal_peclet_number(), 0.0);
        assert_eq!(reactor.temperature_shift(), 100.0);
        assert_eq!(reactor.characteristic_length(), 0.1);
        assert_eq!(reactor.temperature_scale(), 100.0);
        assert_eq!(reactor.scaling_config().dT, 100.0);
        assert_eq!(reactor.scaling_config().L, 0.1);
        assert_eq!(reactor.scaling_config().T_scale, 100.0);
    }

    #[test]
    fn test_scaling_config_accessor_tracks_updates() {
        let mut reactor = SimpleReactorTask::new();

        reactor
            .set_scaling_values(120.0, 0.25, 80.0)
            .expect("scaling should validate");

        assert_eq!(reactor.temperature_shift(), 120.0);
        assert_eq!(reactor.characteristic_length(), 0.25);
        assert_eq!(reactor.temperature_scale(), 80.0);
        let scaling = reactor.scaling_config();
        assert_eq!(scaling.dT, 120.0);
        assert_eq!(scaling.L, 0.25);
        assert_eq!(scaling.T_scale, 80.0);
    }

    #[test]
    fn test_scaling_config_roundtrip_and_validation() {
        let scaling = ScalingConfig::new(120.0, 0.25, 80.0);

        assert!(scaling.validate().is_ok());

        let map = scaling.to_hashmap();
        assert_eq!(map.get("dT"), Some(&120.0));
        assert_eq!(map.get("L"), Some(&0.25));
        assert_eq!(map.get("T_scale"), Some(&80.0));

        let restored =
            ScalingConfig::from_hashmap(&map).expect("scaling should round-trip through a hashmap");
        assert_eq!(restored.dT, 120.0);
        assert_eq!(restored.L, 0.25);
        assert_eq!(restored.T_scale, 80.0);
    }

    #[test]
    fn test_scaling_config_rejects_non_positive_values() {
        let invalid_temperature = ScalingConfig::new(0.0, 0.25, 80.0);
        assert!(invalid_temperature.validate().is_err());

        let invalid_length = ScalingConfig::new(120.0, -0.25, 80.0);
        assert!(invalid_length.validate().is_err());
    }

    #[test]
    fn test_scaling_config_from_hashmap_requires_temperature_scale() {
        let map = HashMap::from([("dT".to_string(), 120.0), ("L".to_string(), 0.25)]);

        let err = ScalingConfig::from_hashmap(&map).unwrap_err();
        assert!(
            err.to_string().contains("T_scale"),
            "missing T_scale should be reported clearly"
        );
    }

    #[test]
    fn test_bvp_solver_default() {
        let solver = BVPSolver::default();

        assert_eq!(solver.arg_name, "x");
        assert_eq!(solver.x_range, (0.0, 1.0));
        assert!(solver.unknowns.is_empty());
        assert!(solver.eq_system.is_empty());
        assert!(solver.BorderConditions.is_empty());
        //  assert!(solver.rhs.is_none());
    }

    #[test]
    fn test_le_number_uses_current_diffusion_snapshot() {
        let mut reactor = SimpleReactorTask::new();
        reactor.Lambda = 0.05;
        reactor.Cp = 1000.0;
        reactor.D_ro_map = HashMap::from([("A".to_string(), 1.0e-5), ("B".to_string(), 2.0e-5)]);

        let le_number = reactor.Le_number().unwrap();
        assert_relative_eq!(le_number, 3.75, max_relative = 1e-12);
    }

    #[test]
    fn test_le_number_rejects_missing_diffusion_map() {
        let mut reactor = SimpleReactorTask::new();
        reactor.Lambda = 0.05;
        reactor.Cp = 1000.0;

        let result = reactor.Le_number();
        assert!(matches!(result, Err(ReactorError::MissingData(_))));
    }

    #[test]
    fn test_with_real_data() {
        let mut kd = KinData::new();
        kd.set_reactions_from_shortcut_range("C1..C3".to_string());
        kd.get_reactions_from_shortcuts();
        kd.reactdata_parsing();
        if let Some(reactdata) = kd.vec_of_reaction_data.as_mut() {
            reactdata.retain(|rd| rd.reaction_type == ReactionType::Elem);
        }
        // we need only elementaty reactions
        let mut kd2 = KinData::new();
        kd2.vec_of_reaction_data = kd.vec_of_reaction_data;
        kd2.equations_from_reactdata().unwrap();
        kd2.analyze_reactions().unwrap();

        let mut reactor = SimpleReactorTask::new();

        reactor.kindata = kd2;

        info!("kindata \n {:?}", reactor.kindata);
        info!("substances {:?}", &reactor.kindata.substances);
        info!("\n  eq {:?}", &reactor.kindata.vec_of_equations);
        // ["H", "O", "OH", "H2", "HO2", "O2"]

        let P = 1e5; // Pa
        let Tm = 1500.0;
        let m = 1e-2; //
        let Cp = 1000.0; // J/kg/K

        let Lambda = 0.027;
        let Diffusion = HashMap::from([
            ("H".to_string(), 1e-3),
            ("O".to_string(), 1e-3),
            ("OH".to_string(), 1e-3),
            ("H2".to_string(), 1e-3),
            ("HO2".to_string(), 1e-3),
            ("O2".to_string(), 1e-3),
        ]);
        let boundary_condition = HashMap::from([
            ("H".to_string(), 0.6),
            ("O".to_string(), 0.4),
            ("OH".to_string(), 1e-3),
            ("H2".to_string(), 1e-3),
            ("HO2".to_string(), 1e-3),
            ("O2".to_string(), 1e-3),
            ("T".to_string(), 450.0),
        ]);
        let thermal_effects = vec![1e4, 1e4];
        let scaling = ScalingConfig::new(100.0, 1e-5, 100.0);
        reactor.set_parameters(
            thermal_effects,
            P,
            Tm,
            Cp,
            boundary_condition,
            Lambda,
            Diffusion,
            m,
            scaling,
        );
        // info!("reactor {:?}", reactor);
        let res = reactor.setup_bvp();
        match res {
            Ok(_) => info!("ok"),
            Err(e) => info!("error {:?}", e),
        }
        let rates = reactor.map_eq_rate;
        for (eq, rate) in rates {
            info!("reaction {} rate {}", eq, rate);
        }
        info!("\n \n");
        let system = reactor.map_of_equations;
        for (subs, (variable, eq)) in system {
            info!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = reactor.solver.BorderConditions;
        info!("bc {:?}", bc);

        // use RustedSciThe::symbolic::symbolic_engine::Expr;
        // let mut A = Expr::parse_expression("-A +B*C");
        //   let newA = A.substitute_variable("A", &Expr::Var("D".to_string()));
        //   info!("newA {:?}", newA);
    }

    #[test]
    fn hmx_test() {
        /////////////////// setting up kinetics
        let mut kd = KinData::new();
        kd.substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        let hmx = HashMap::from([
            ("H".to_string(), 4),
            ("N".to_string(), 8),
            ("C".to_string(), 8),
            ("O".to_string(), 8),
        ]);
        kd.groups = Some(HashMap::from([
            ("HMX".to_string(), hmx.clone()),
            ("HMXprod".to_string(), hmx),
        ]));
        let eq = "HMX=>HMXprod".to_string();
        kd.vec_of_equations = vec![eq.clone()];
        //////////// instance of problem constructor/////////////
        let mut reactor = SimpleReactorTask::new();
        reactor.kindata = kd;
        ///////////// parameters //////////////////
        let Q_g = 3000.0;
        let C_p = 0.35 * 4.184;
        let Lambda_eff = 0.07; // W/m-K 
        let M = 34.2 / 1000.0; //  kg/mol
        let A = 1.3 * 1e5; //
        let E = 5000.0 * 4.184; //  
        let Diffusion = HashMap::from([("HMX".to_string(), 1e-3), ("HMXprod".to_string(), 1e-3)]);
        let boundary_condition = HashMap::from([
            ("HMX".to_string(), 1.0 - 1e-3),
            ("HMXprod".to_string(), 1e-3),
            ("T".to_string(), 500.0),
        ]);
        let thermal_effects = vec![Q_g];
        let P = 1e5; // Pa
        let Tm = 1500.0;
        let m = 1e-3; // м/с
        let scaling = ScalingConfig::new(100.0, 1e-5, 100.0);
        let arrenius = vec![A, 0.0, E];
        let reactdata = ReactionData::new_elementary(eq.clone(), arrenius, None);
        reactor.kindata.vec_of_reaction_data = Some(vec![reactdata]);

        reactor.set_parameters(
            thermal_effects,
            P,
            Tm,
            C_p,
            boundary_condition,
            Lambda_eff,
            Diffusion,
            m,
            scaling,
        );
        let res = reactor.setup_bvp();
        match res {
            Ok(_) => info!("ok"),
            Err(e) => info!("error {:?}", e),
        }
        let rates = reactor.map_eq_rate;
        for (eq, rate) in rates {
            info!("reaction {} rate {}", eq, rate);
        }
        info!("\n \n");
        let system = reactor.map_of_equations;
        for (subs, (variable, eq)) in system {
            info!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = reactor.solver.BorderConditions;
        info!("bc {:?}", bc);
    }

    fn create_hmx(Q_g: f64, L: f64) -> SimpleReactorTask {
        /////////////////// setting up kinetics
        let eq = "HMX=>10HMXprod".to_string();

        let C_p = 0.35 * 4.184 * 1000.0;
        let Lambda_eff = 0.07; // W/m-K 
        let n = 0.0;
        let M = 34.2 / 1000.0; //  kg/mol
        let A = 1.3e5; // 1.3 * 1e5; //
        let E = 5000.0 * 4.184; // 5000.0 * 4.184; //  

        let T0 = 800.0;
        let T_scale = 600.0;
        let P: f64 = 1e6; // Pa
        let Tm = 1500.0;
        let m = 0.077 * (P / 1e5).powf(0.748) / 1e2; // 1000 for sm/s ->м/ы
        let hmx = HashMap::from([
            ("H".to_string(), 4),
            ("N".to_string(), 8),
            ("C".to_string(), 8),
            ("O".to_string(), 8),
        ]);
        let hmxprod = HashMap::from([
            ("H".to_string(), 6),
            ("C".to_string(), 1),
            ("O".to_string(), 1),
        ]);
        let groups = Some(HashMap::from([
            ("HMX".to_string(), hmx.clone()),
            ("HMXprod".to_string(), hmxprod),
        ]));
        //////////// instance of problem constructor/////////////
        let mut reactor = SimpleReactorTask::new();
        let struct_with_params = FastElemReact {
            eq,
            A,
            n,
            E,
            Q: Q_g,
        };
        let vec_of_structs = vec![struct_with_params];

        let _ = reactor.fast_react_set(vec_of_structs);
        reactor.kindata.substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        reactor.kindata.groups = groups;
        // info!("reactor {:?}", reactor.kindata);
        ///////////// parameters //////////////////
        // we assume Le == 1
        let ro0 = M * P / (R_G * T0);
        let D = Lambda_eff / (C_p * ro0);
        info!("D = {}", D);
        // set diffusion coeffisients
        let Diffusion = HashMap::from([("HMX".to_string(), D), ("HMXprod".to_string(), D)]);
        let boundary_condition = HashMap::from([
            ("HMX".to_string(), 1.0 - 1e-3),
            ("HMXprod".to_string(), 1e-3),
            ("T".to_string(), T0),
        ]);
        let thermal_effects = vec![Q_g];
        let scaling = ScalingConfig::new(T_scale, L, T_scale);
        reactor.set_parameters(
            thermal_effects,
            P,
            Tm,
            C_p,
            boundary_condition,
            Lambda_eff,
            Diffusion,
            m,
            scaling,
        );

        reactor.M = M; // calc of M is scipped if it is set manually
        let res = reactor.setup_bvp();
        let _ = reactor.Le_number();
        info!("reactor {:?}", reactor.kindata);
        match res {
            Ok(_) => info!("ok"),
            Err(e) => info!("error {:?}", e),
        }
        let rates = reactor.map_eq_rate.clone();
        for (eq, rate) in rates {
            info!("reaction {} rate {}", eq, rate);
        }
        info!("\n \n");
        let system = reactor.map_of_equations.clone();
        for (subs, (variable, eq)) in system {
            info!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = &reactor.solver.BorderConditions;
        info!("bc {:?}", bc);
        info!(" unknowns{:?}", reactor.solver.unknowns);

        // Override molar masses and recreate equations
        //   reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![340.0, 34.0]);
        //    let _ = reactor.create_bvp_equations();

        reactor
    }

    /// Solve the compact HMX fixture through the new default reactor facade.
    ///
    /// We keep the legacy combustion setup, but the actual backend path now uses
    /// the canonical `Lambdify + AtomView + Banded` solver configuration so the
    /// old smoke tests exercise the same default route as production code.
    fn solve_hmx_with_default_backend(
        reactor: &mut SimpleReactorTask,
        initial_guess: DMatrix<f64>,
        n_steps: usize,
        strategy_params: Option<SolverParams>,
        tolerance_config: ToleranceConfig,
        bounds_config: BoundsConfig,
        max_iterations: usize,
        abs_tolerance: f64,
        loglevel: Option<String>,
    ) -> Result<(), ReactorError> {
        let substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        reactor.solver.solve_NRBVP_with_configs(
            initial_guess,
            n_steps,
            "forward".to_string(),
            "Damped".to_string(),
            strategy_params,
            None,
            "Banded".to_string(),
            abs_tolerance,
            tolerance_config,
            bounds_config,
            &substances,
            max_iterations,
            loglevel,
        )
    }

    #[test]
    fn hmx_test2() {
        let Q_g = 3000.0 * 1e3 / 100.0; // J/kg -> J/mole 
        let L = 9e-4;
        let mut reactor = create_hmx(Q_g, L);
        let n_steps = 100;

        let grid_method = GridRefinementMethod::Pearson(0.0, 3.5);
        // GridRefinementMethod::GrcarSmooke(0.1, 0.1, 3.5);
        //GridRefinementMethod::Pearson(0.09, 3.5);
        // or GridRefinementMethod::Pearson(0.05, 2.5);
        let adaptive = AdaptiveGridConfig {
            version: 1,
            max_refinements: 2,
            grid_method,
        };

        let strategy_params = SolverParams {
            max_jac: Some(3),
            max_damp_iter: Some(3),
            damp_factor: Some(0.5),
            adaptive: None,
        };

        let max_iterations = 100;
        let abs_tolerance = 1e-6;

        let loglevel = Some("info".to_string());
        let scheme = "forward".to_string();
        let method = "Banded".to_string();
        let strategy = "Damped".to_string();
        let linear_sys_method = None;
        // Using the new tolerance helper - much simpler!
        let tolerance_config = ToleranceConfig::new(1e-5, 1e-5, 1e-5, 1e-6);
        // Using the new bounds helper
        let bounds_config = BoundsConfig::new(
            (-0.1, 1.1),     // C bounds
            (-1e20, 1e20),   // J bounds
            (-100.0, 100.0), // Teta bounds
            (-1e20, 1e20),   // q bounds
        );
        let ig = vec![0.99; n_steps * reactor.solver.unknowns.len()];
        let initial_guess = DMatrix::from_vec(reactor.solver.unknowns.len(), n_steps, ig);
        // Using the new convenience method with both configs:
        let _ = reactor.solver.solve_NRBVP_with_configs(
            initial_guess,
            n_steps,
            scheme,
            strategy,
            Some(strategy_params),
            linear_sys_method,
            method,
            abs_tolerance,
            tolerance_config,
            bounds_config,
            &vec!["HMX".to_string(), "HMXprod".to_string()],
            max_iterations,
            loglevel,
        );
        assert!(
            reactor.solver.solution.is_some(),
            "solve_NRBVP_with_configs should store the computed solution"
        );
        assert!(
            reactor.solver.x_mesh.is_some(),
            "solve_NRBVP_with_configs should store the computed mesh"
        );
        //  reactor.postprocessing();
        reactor.gnuplot().unwrap();
        info!("BC {:?}", &reactor.solver.BorderConditions);
        info!("unknowns {:?}", &reactor.solver.unknowns);
        for (i, eq) in reactor.solver.eq_system.clone().iter().enumerate() {
            info!("y = {}, eq {}", reactor.solver.unknowns[i], eq);
        }
        let sol = &reactor.solver.solution.clone().unwrap();
        let T = sol.column(0).clone();
        //  info!("T = {}", T);
        //  info!("q = {}", sol.column(1));
        // info!("C0 = {}", sol.column(2));
    }

    #[test]
    fn test_hmx_molar_mass_bug() {
        use crate::Kinetics::molmass::calculate_molar_mass_of_vector_of_subs;
        // Test the molar mass calculation directly
        let hmx = HashMap::from([
            ("H".to_string(), 4),
            ("N".to_string(), 8),
            ("C".to_string(), 8),
            ("O".to_string(), 8),
        ]);
        let hmxprod = HashMap::from([
            ("H".to_string(), 6),
            ("C".to_string(), 1),
            ("O".to_string(), 1),
        ]);
        let groups = Some(HashMap::from([
            ("HMX".to_string(), hmx.clone()),
            ("HMXprod".to_string(), hmxprod),
        ]));

        let vec_of_formulae = vec!["HMX", "HMXprod"];
        let molar_masses = calculate_molar_mass_of_vector_of_subs(vec_of_formulae, groups).unwrap();
        info!("Test result: {:?}", molar_masses);

        // Expected: HMX should be much larger than HMXprod
        assert!(molar_masses[0] > molar_masses[1]);
        assert!((molar_masses[0] - 340.0).abs() < 1.0); // HMX: 4*1.008 + 8*14.007 + 8*12.011 + 8*15.999 = 344
        assert!((molar_masses[1] - 34.04).abs() < 1.0); // HMXprod: 6*1.008 + 1*12.011 + 1*15.999
    }

    #[test]
    fn hmx_test3() {
        let L = 9e-4;
        let Q_g = 3000.0 * 1e3 / 100.0; // J/kg -> J/mole
        let mut reactor = create_hmx(Q_g, L);
        let n_steps = 300;
        let ig = vec![0.99; n_steps * reactor.solver.unknowns.len()];
        let initial_guess = DMatrix::from_vec(reactor.solver.unknowns.len(), n_steps, ig);
        let tolerance_config = ToleranceConfig::new(1e-5, 1e-5, 1e-5, 1e-6);
        let bounds_config = BoundsConfig::new(
            (-0.1, 1.1),     // C bounds
            (-1e20, 1e20),   // J bounds
            (-100.0, 100.0), // Teta bounds
            (-1e20, 1e20),   // q bounds
        );
        let _ = solve_hmx_with_default_backend(
            &mut reactor,
            initial_guess,
            n_steps,
            None,
            tolerance_config,
            bounds_config,
            n_steps * 2000,
            1e-6,
            Some("info".to_string()),
        );
        assert!(
            reactor.solver.solution.is_some(),
            "solve_NRBVP_with_configs should store the computed solution without extra side effects"
        );
        assert!(
            reactor.solver.x_mesh.is_some(),
            "solve_NRBVP_with_configs should store the computed mesh without extra side effects"
        );
        //  reactor.postprocessing();
        //  reactor.gnuplot();
        info!("BC {:?}", &reactor.solver.BorderConditions);
        info!("unknowns {:?}", &reactor.solver.unknowns);
        for (i, eq) in reactor.solver.eq_system.clone().iter().enumerate() {
            info!("y = {}, eq {}", reactor.solver.unknowns[i], eq);
        }
        // let sol = &reactor.solver.solution.clone().unwrap();
        //let T = sol.column(0).clone();
        //  info!("T = {}", T);
        //  info!("q = {}", sol.column(1));
        // info!("C0 = {}", sol.column(2));
    }

    #[test]
    fn hmx_test4() {
        let Q_g = 3000.0 * 1e3 * 0.034; // J/kg -> J/mole
        let L = 2.7e-4;
        let mut reactor = create_hmx(Q_g, L);

        let n_steps = 50;
        let grid_method = GridRefinementMethod::GrcarSmooke(0.1, 0.1, 2.5);
        // GridRefinementMethod::GrcarSmooke(0.1, 0.1, 3.5);
        //GridRefinementMethod::Pearson(0.09, 3.5);
        // or GridRefinementMethod::Pearson(0.05, 2.5);
        let adaptive = AdaptiveGridConfig {
            version: 1,
            max_refinements: 3,
            grid_method,
        };

        let strategy_params = SolverParams {
            max_jac: Some(3),
            max_damp_iter: Some(10),
            damp_factor: Some(0.5),
            adaptive: Some(adaptive),
        };

        let max_iterations = 100;
        let abs_tolerance = 1e-8;

        let loglevel = Some("info".to_string());
        let scheme = "forward".to_string();
        let method = "Banded".to_string();
        let strategy = "Damped".to_string();
        let linear_sys_method = None;
        // Using the new tolerance helper - much simpler!
        let tolerance_config = ToleranceConfig::new(1e-7, 1e-7, 1e-7, 1e-6);
        let substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        // Using the new bounds helper
        let bounds_config = BoundsConfig::new(
            (-100.0, 100.1), // C bounds
            (-1e20, 1e20),   // J bounds
            (-100.0, 100.0), // Teta bounds
            (-1e20, 1e20),   // q bounds
        );
        let ig = vec![1e-5; n_steps * reactor.solver.unknowns.len()];
        let initial_guess = DMatrix::from_vec(reactor.solver.unknowns.len(), n_steps, ig);
        // Using the new convenience method with both configs:
        let _ = reactor.solver.solve_NRBVP_with_configs(
            initial_guess,
            n_steps,
            scheme,
            strategy,
            Some(strategy_params),
            linear_sys_method,
            method,
            abs_tolerance,
            tolerance_config,
            bounds_config,
            &substances,
            max_iterations,
            loglevel,
        );

        //  reactor.postprocessing();
        reactor.gnuplot().unwrap();
        //  reactor.plot_in_terminal();
        info!("BC {:?}", &reactor.solver.BorderConditions);
        info!("unknowns {:?}", &reactor.solver.unknowns);
        for (i, eq) in reactor.solver.eq_system.clone().iter().enumerate() {
            info!("y = {}, eq {}", reactor.solver.unknowns[i], eq);
        }
        let sol = &reactor.solver.solution.clone().unwrap();
        let T = sol.column(0).clone();
        //  info!("T = {}", T);
        //  info!("q = {}", sol.column(1));
    }

    #[test]
    fn hmx_test5() {
        let L = 5e-4;
        let Q_g = 3000.0 * 1e3 * 0.034; // J/kg -> J/mole
        let mut reactor = create_hmx(Q_g, L);

        let n_steps = 50;
        let ig = vec![1e-3; n_steps * reactor.solver.unknowns.len()];
        let initial_guess = DMatrix::from_vec(reactor.solver.unknowns.len(), n_steps, ig);
        let tolerance_config = ToleranceConfig::new(1e-5, 1e-5, 1e-5, 1e-6);
        let bounds_config = BoundsConfig::new(
            (-0.1, 1.1),     // C bounds
            (-1e20, 1e20),   // J bounds
            (-100.0, 100.0), // Teta bounds
            (-1e20, 1e20),   // q bounds
        );
        let _ = solve_hmx_with_default_backend(
            &mut reactor,
            initial_guess,
            n_steps,
            None,
            tolerance_config,
            bounds_config,
            n_steps * 50,
            1e-6,
            Some("info".to_string()),
        );
        //  reactor.postprocessing();
        //  reactor.gnuplot();
        info!("BC {:?}", &reactor.solver.BorderConditions);
        info!("unknowns {:?}", &reactor.solver.unknowns);
        for (i, eq) in reactor.solver.eq_system.clone().iter().enumerate() {
            info!("y = {}, eq {}", reactor.solver.unknowns[i], eq);
        }
        info!("Pe_D = {:?}, Pe_q ={}", &reactor.Pe_D, &reactor.Pe_q)
        // let sol = &reactor.solver.solution.clone().unwrap();
        //let T = sol.column(0).clone();
        //  info!("T = {}", T);
        //  info!("q = {}", sol.column(1));
        // info!("C0 = {}", sol.column(2));
    }

    #[test]
    fn test_tolerance_helpers() {
        // Test ToleranceConfig struct approach
        let tolerance_config = ToleranceConfig::new(1e-4, 1e-4, 1e-5, 1e-4);
        let substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        let full_tolerance_map = tolerance_config.to_full_tolerance_map(&substances);

        assert_eq!(full_tolerance_map.get("Teta"), Some(&1e-5));
        assert_eq!(full_tolerance_map.get("q"), Some(&1e-4));
        assert_eq!(full_tolerance_map.get("C0"), Some(&1e-4));
        assert_eq!(full_tolerance_map.get("C1"), Some(&1e-4));
        assert_eq!(full_tolerance_map.get("J0"), Some(&1e-4));
        assert_eq!(full_tolerance_map.get("J1"), Some(&1e-4));

        // Test function approach
        let simple_config = HashMap::from([
            ("C".to_string(), 1e-4),
            ("J".to_string(), 1e-4),
            ("Teta".to_string(), 1e-5),
            ("q".to_string(), 1e-4),
        ]);
        let full_tolerance_map2 = create_tolerance_map(simple_config, &substances);

        assert_eq!(full_tolerance_map, full_tolerance_map2);

        // Test default values
        let default_config = ToleranceConfig::default();
        let default_map = default_config.to_full_tolerance_map(&substances);
        assert_eq!(default_map.get("Teta"), Some(&1e-5));
        assert_eq!(default_map.get("q"), Some(&1e-4));
        assert_eq!(default_map.get("C0"), Some(&1e-4));
        assert_eq!(default_map.get("J0"), Some(&1e-4));
    }

    #[test]
    fn test_reactor_tolerance_helper() {
        let mut reactor = create_test_reactor();
        let _ = reactor.kinetic_processing();

        let simple_config = HashMap::from([
            ("C".to_string(), 1e-6),
            ("J".to_string(), 1e-5),
            ("Teta".to_string(), 1e-7),
            ("q".to_string(), 1e-5),
        ]);

        let full_tolerance_map = reactor.create_tolerance_map_for_system(simple_config);

        // Should have entries for all substances
        assert!(full_tolerance_map.contains_key("C0")); // A
        assert!(full_tolerance_map.contains_key("C1")); // B  
        assert!(full_tolerance_map.contains_key("C2")); // C
        assert!(full_tolerance_map.contains_key("J0"));
        assert!(full_tolerance_map.contains_key("J1"));
        assert!(full_tolerance_map.contains_key("J2"));
        assert!(full_tolerance_map.contains_key("Teta"));
        assert!(full_tolerance_map.contains_key("q"));

        assert_eq!(full_tolerance_map.get("Teta"), Some(&1e-7));
        assert_eq!(full_tolerance_map.get("q"), Some(&1e-5));
        assert_eq!(full_tolerance_map.get("C0"), Some(&1e-6));
        assert_eq!(full_tolerance_map.get("J0"), Some(&1e-5));
    }

    #[test]
    fn test_bounds_helpers() {
        // Test BoundsConfig struct approach
        let bounds_config = BoundsConfig::new(
            (0.0, 1.0),    // C bounds
            (-1e20, 1e20), // J bounds
            (-10.0, 10.0), // Teta bounds
            (-1e20, 1e20), // q bounds
        );
        let substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        let full_bounds_map = bounds_config.to_full_bounds_map(&substances);

        assert_eq!(full_bounds_map.get("Teta"), Some(&(-10.0, 10.0)));
        assert_eq!(full_bounds_map.get("q"), Some(&(-1e20, 1e20)));
        assert_eq!(full_bounds_map.get("C0"), Some(&(0.0, 1.0)));
        assert_eq!(full_bounds_map.get("C1"), Some(&(0.0, 1.0)));
        assert_eq!(full_bounds_map.get("J0"), Some(&(-1e20, 1e20)));
        assert_eq!(full_bounds_map.get("J1"), Some(&(-1e20, 1e20)));

        // Test function approach
        let simple_config = HashMap::from([
            ("C".to_string(), (0.0, 1.0)),
            ("J".to_string(), (-1e20, 1e20)),
            ("Teta".to_string(), (-10.0, 10.0)),
            ("q".to_string(), (-1e20, 1e20)),
        ]);
        let full_bounds_map2 = create_bounds_map(simple_config, &substances);

        assert_eq!(full_bounds_map, full_bounds_map2);

        // Test default values
        let default_config = BoundsConfig::default();
        let default_map = default_config.to_full_bounds_map(&substances);
        assert_eq!(default_map.get("Teta"), Some(&(-10.0, 10.0)));
        assert_eq!(default_map.get("q"), Some(&(-1e20, 1e20)));
        assert_eq!(default_map.get("C0"), Some(&(0.0, 1.0)));
        assert_eq!(default_map.get("J0"), Some(&(-1e20, 1e20)));
    }

    #[test]
    fn test_reactor_bounds_helper() {
        let mut reactor = create_test_reactor();
        let _ = reactor.kinetic_processing();

        let simple_config = HashMap::from([
            ("C".to_string(), (0.0, 2.0)),
            ("J".to_string(), (-1e10, 1e10)),
            ("Teta".to_string(), (-5.0, 5.0)),
            ("q".to_string(), (-1e15, 1e15)),
        ]);

        let full_bounds_map = reactor.create_bounds_map_for_system(simple_config);

        // Should have entries for all substances
        assert!(full_bounds_map.contains_key("C0")); // A
        assert!(full_bounds_map.contains_key("C1")); // B  
        assert!(full_bounds_map.contains_key("C2")); // C
        assert!(full_bounds_map.contains_key("J0"));
        assert!(full_bounds_map.contains_key("J1"));
        assert!(full_bounds_map.contains_key("J2"));
        assert!(full_bounds_map.contains_key("Teta"));
        assert!(full_bounds_map.contains_key("q"));

        assert_eq!(full_bounds_map.get("Teta"), Some(&(-5.0, 5.0)));
        assert_eq!(full_bounds_map.get("q"), Some(&(-1e15, 1e15)));
        assert_eq!(full_bounds_map.get("C0"), Some(&(0.0, 2.0)));
        assert_eq!(full_bounds_map.get("J0"), Some(&(-1e10, 1e10)));
    }

    #[test]
    fn test_molar_concentration_alias_matches_legacy_spelling() {
        let mut reactor = create_test_reactor();
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![18.0, 28.0]);

        // Two rows are enough to confirm the alias and the legacy spelling return the same data.
        let matrix_of_mass_fractions = DMatrix::from_row_slice(2, 2, &[0.2, 0.8, 0.6, 0.4]);

        let corrected = reactor
            .from_mass_fractions_to_molar_concentration(matrix_of_mass_fractions.clone())
            .expect("correctly spelled helper should succeed");
        let legacy = reactor
            .from_mass_fractions_to_molar_conentration(matrix_of_mass_fractions)
            .expect("legacy helper should still work");

        assert_eq!(corrected, legacy);
        assert_relative_eq!(corrected[(0, 0)], 0.2 / 0.018, epsilon = 1e-12);
        assert_relative_eq!(corrected[(0, 1)], 0.8 / 0.028, epsilon = 1e-12);
    }
}
