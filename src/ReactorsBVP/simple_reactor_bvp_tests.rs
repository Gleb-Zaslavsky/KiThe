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
    use nalgebra::DMatrix;
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

        let M = reactor.M;
        println!("Molar mass: {}", M);
        let coeffs = reactor.transport_coefficients();
        println!("Coefficients: {:?}", coeffs);
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
        let _ = reactor.transport_coefficients();
        let subs = reactor.kindata.substances.clone();
        assert! {subs == vec!["A".to_string(), "B".to_string(), "C".to_string()]};
        println!("{:?}", subs);
        // Setup substances for transport coefficient calculation
        //let mut kindata = KinData::new();
        //kindata.substances = vec!["A".to_string(), "B".to_string()];
        //reactor.kindata = kindata;

        let result = reactor.peclet_numbers();
        assert!(result.is_ok());

        // Pe_q = (L * m * Cp) / Lambda
        let Cp = reactor.Cp;
        let m = reactor.m;
        let L = reactor.L;
        let Lambda = reactor.Lambda;
        println!("Cp {}, m {}, L {}, Lambda {}", Cp, m, L, Lambda);
        let expected_pe_q = (0.1 * 0.01 * 1000.0) / 0.05;
        assert!((reactor.Pe_q - expected_pe_q).abs() < 1e-10);

        // Check Pe_D vector length
        assert_eq!(reactor.Pe_D.len(), 3);

        // Verify Pe_D calculations
        let transport_coeffs = reactor.transport_coefficients();
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
    fn test_ideal_gas_density() {
        let mut reactor = create_test_reactor();
        let _ = reactor.kinetic_processing();
        reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![28.0, 44.0, 20.0]);
        let _ = reactor.mean_molar_mass();

        let M = 1.0 / (0.5 / 28.0 + 0.3 / 44.0 + 0.2 / 20.0);
        let M = M / 1000.0;
        let M_ = reactor.M;
        assert_relative_eq!(M, M_, epsilon = 1e-6);
        let density = reactor.ideal_gas_density();
        println!("Density: {}, M= {}", density, &reactor.M);
        let expected = M * reactor.P / (R_G * reactor.Tm);

        assert!((density - expected).abs() < 1e-10);
    }

    #[test]
    fn test_reactor_error_display() {
        let errors = vec![
            ReactorError::MissingData("test data".to_string()),
            ReactorError::InvalidConfiguration("test config".to_string()),
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
        kd2.equations_from_reactdata();
        kd2.analyze_reactions();

        let mut reactor = SimpleReactorTask::new();

        reactor.kindata = kd2;

        println!("kindata \n {:?}", reactor.kindata);
        println!("substances {:?}", &reactor.kindata.substances);
        println!("\n  eq {:?}", &reactor.kindata.vec_of_equations);
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
        // println!("reactor {:?}", reactor);
        let res = reactor.setup_bvp();
        match res {
            Ok(_) => println!("ok"),
            Err(e) => println!("error {:?}", e),
        }
        let rates = reactor.map_eq_rate;
        for (eq, rate) in rates {
            println!("reaction {} rate {}", eq, rate);
        }
        println!("\n \n");
        let system = reactor.map_of_equations;
        for (subs, (variable, eq)) in system {
            println!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = reactor.solver.BorderConditions;
        println!("bc {:?}", bc);

        // use RustedSciThe::symbolic::symbolic_engine::Expr;
        // let mut A = Expr::parse_expression("-A +B*C");
        //   let newA = A.substitute_variable("A", &Expr::Var("D".to_string()));
        //   println!("newA {:?}", newA);
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
            Ok(_) => println!("ok"),
            Err(e) => println!("error {:?}", e),
        }
        let rates = reactor.map_eq_rate;
        for (eq, rate) in rates {
            println!("reaction {} rate {}", eq, rate);
        }
        println!("\n \n");
        let system = reactor.map_of_equations;
        for (subs, (variable, eq)) in system {
            println!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = reactor.solver.BorderConditions;
        println!("bc {:?}", bc);
    }

    fn create_hmx(Q_g: f64, L: f64) -> SimpleReactorTask {
        /////////////////// setting up kinetics
        let eq = "HMX=>HMXprod".to_string();

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
        // println!("reactor {:?}", reactor.kindata);
        ///////////// parameters //////////////////
        // we assume Le == 1
        let ro0 = M * P / (R_G * T0);
        let D = Lambda_eff / (C_p * ro0);
        println!("D = {}", D);
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
        let _ = &reactor.Le_number();
        println!("reactor {:?}", reactor.kindata);
        match res {
            Ok(_) => println!("ok"),
            Err(e) => println!("error {:?}", e),
        }
        let rates = reactor.map_eq_rate.clone();
        for (eq, rate) in rates {
            println!("reaction {} rate {}", eq, rate);
        }
        println!("\n \n");
        let system = reactor.map_of_equations.clone();
        for (subs, (variable, eq)) in system {
            println!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = &reactor.solver.BorderConditions;
        println!("bc {:?}", bc);
        println!(" unknowns{:?}", reactor.solver.unknowns);

        // Override molar masses and recreate equations
        //   reactor.kindata.stecheodata.vec_of_molmasses = Some(vec![340.0, 34.0]);
        //    let _ = reactor.create_bvp_equations();

        reactor
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
        let method = "Sparse".to_string();
        let strategy = "Damped".to_string();
        let linear_sys_method = None;
        // Using the new tolerance helper - much simpler!
        let tolerance_config = ToleranceConfig::new(1e-5, 1e-5, 1e-5, 1e-6);
        let substances = vec!["HMX".to_string(), "HMXprod".to_string()];
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
            &substances,
            max_iterations,
            loglevel,
        );
        //  reactor.postprocessing();
        reactor.gnuplot();
        println!("BC {:?}", &reactor.solver.BorderConditions);
        println!("unknowns {:?}", &reactor.solver.unknowns);
        for (i, eq) in reactor.solver.eq_system.clone().iter().enumerate() {
            println!("y = {}, eq {}", reactor.solver.unknowns[i], eq);
        }
        let sol = &reactor.solver.solution.clone().unwrap();
        let T = sol.column(0).clone();
        //  println!("T = {}", T);
        //  println!("q = {}", sol.column(1));
        // println!("C0 = {}", sol.column(2));
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
        let molar_masses = calculate_molar_mass_of_vector_of_subs(vec_of_formulae, groups);
        println!("Test result: {:?}", molar_masses);

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
        let _ = reactor
            .solver
            .solve_BVPsci(initial_guess, n_steps, n_steps * 2000, 1e-6);
        //  reactor.postprocessing();
        //  reactor.gnuplot();
        println!("BC {:?}", &reactor.solver.BorderConditions);
        println!("unknowns {:?}", &reactor.solver.unknowns);
        for (i, eq) in reactor.solver.eq_system.clone().iter().enumerate() {
            println!("y = {}, eq {}", reactor.solver.unknowns[i], eq);
        }
        // let sol = &reactor.solver.solution.clone().unwrap();
        //let T = sol.column(0).clone();
        //  println!("T = {}", T);
        //  println!("q = {}", sol.column(1));
        // println!("C0 = {}", sol.column(2));
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
        let method = "Sparse".to_string();
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
        reactor.gnuplot();
        println!("BC {:?}", &reactor.solver.BorderConditions);
        println!("unknowns {:?}", &reactor.solver.unknowns);
        for (i, eq) in reactor.solver.eq_system.clone().iter().enumerate() {
            println!("y = {}, eq {}", reactor.solver.unknowns[i], eq);
        }
        let sol = &reactor.solver.solution.clone().unwrap();
        let T = sol.column(0).clone();
        //  println!("T = {}", T);
        //  println!("q = {}", sol.column(1));
    }

    #[test]
    fn hmx_test5() {
        let L = 5e-4;
        let Q_g = 3000.0 * 1e3 * 0.034; // J/kg -> J/mole
        let mut reactor = create_hmx(Q_g, L);

        let n_steps = 50;
        let ig = vec![1e-3; n_steps * reactor.solver.unknowns.len()];
        let initial_guess = DMatrix::from_vec(reactor.solver.unknowns.len(), n_steps, ig);
        let _ = reactor
            .solver
            .solve_BVPsci(initial_guess, n_steps, n_steps * 50, 1e-6);
        //  reactor.postprocessing();
        //  reactor.gnuplot();
        println!("BC {:?}", &reactor.solver.BorderConditions);
        println!("unknowns {:?}", &reactor.solver.unknowns);
        for (i, eq) in reactor.solver.eq_system.clone().iter().enumerate() {
            println!("y = {}, eq {}", reactor.solver.unknowns[i], eq);
        }
        println!("Pe_D = {:?}, Pe_q ={}", &reactor.Pe_D, &reactor.Pe_q)
        // let sol = &reactor.solver.solution.clone().unwrap();
        //let T = sol.column(0).clone();
        //  println!("T = {}", T);
        //  println!("q = {}", sol.column(1));
        // println!("C0 = {}", sol.column(2));
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
}
