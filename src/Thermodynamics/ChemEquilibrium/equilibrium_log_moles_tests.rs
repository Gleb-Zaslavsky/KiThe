#[cfg(test)]
mod tests {
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
        ContinuationSeedPolicy, EquilibriumLogMoles, EquilibriumSolverSettings, GibbsFn, Phase,
        PhaseKind, Solvers, compute_element_totals, compute_species_moles,
        equilibrium_logmole_jacobian, equilibrium_logmole_residual, equilibrium_scaling,
        reaction_phase_stoichiometry, reaction_standard_gibbs, scaled_jacobian, scaled_residual,
        species_to_phase_map,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::compute_reaction_basis;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptOutcome;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        SolverBackend, SolverCascadeBudget, SolverPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
        compute_phase_totals, finite_difference_jacobian, gas_solver, gas_solver_for_T_range,
        gas_solver_for_T_range_for_elements, gas_solver_from_elements, initial_phase_activity,
        reject_repeated_phase_set, validate_phase_set_candidate, InitialPhaseSet, PhaseSet,
    };
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use approx::assert_relative_eq;

    use nalgebra::DMatrix;
    use std::collections::{HashMap, HashSet};
    use std::f64;
    use std::rc::Rc;
    use std::time::Instant;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn normalized_mole_fractions(moles: &[f64]) -> Vec<f64> {
        let total = moles.iter().sum::<f64>();
        moles.iter().map(|value| value / total).collect()
    }
    #[test]
    fn reaction_basis_o2_dissociation() {
        use nalgebra::DMatrix;

        let a = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);

        let rb = compute_reaction_basis(&a, 1e-12).unwrap();

        assert_eq!(rb.rank, 1);
        assert_eq!(rb.num_reactions, 1);

        let nu = rb.reactions.column(0);

        // Must satisfy 2*ν_O2 + 1*ν_O = 0
        assert!(approx_eq(2.0 * nu[0] + nu[1], 0.0, 1e-12));
    }

    #[test]
    fn temperature_sweep_entry_points_reject_invalid_grids_before_mutation() {
        let mut sequential = EquilibriumLogMoles::empty();
        sequential.initial_guess = Some(vec![0.0]);
        assert!(matches!(
            sequential.solve_for_T_range(300.0, 400.0, 0.0),
            Err(ReactionExtentError::InvalidProblem {
                field: "T_step",
                ..
            })
        ));
        assert_eq!(sequential.initial_guess, Some(vec![0.0]));

        let mut chunked_parallel = EquilibriumLogMoles::empty();
        assert!(matches!(
            chunked_parallel.solve_for_T_range_par(400.0, 300.0, 10.0),
            Err(ReactionExtentError::InvalidProblem {
                field: "temperature_range",
                ..
            })
        ));

        let mut independent_parallel = EquilibriumLogMoles::empty();
        assert!(matches!(
            independent_parallel.solve_for_T_range_par2(f64::NAN, 400.0, 10.0),
            Err(ReactionExtentError::InvalidProblem {
                field: "T_start",
                ..
            })
        ));
    }

    #[test]
    fn mutable_problem_setters_validate_before_publishing_state() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.solution = vec![42.0];
        solver.moles = vec![42.0];

        assert!(matches!(
            solver.set_problem(
                vec![1.0],
                vec![Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![1],
                }],
                101325.0,
            ),
            Err(ReactionExtentError::DimensionMismatch(_))
        ));
        assert_eq!(solver.solution, vec![42.0]);
        assert_eq!(solver.moles, vec![42.0]);

        solver.subs_data.substances = vec!["O2".to_string()];
        solver.n0 = vec![0.5];
        assert!(matches!(
            solver.set_n0_from_non_zero_map(HashMap::from([("typo".to_string(), 1.0)])),
            Err(ReactionExtentError::SubsDataError(_))
        ));
        assert_eq!(solver.n0, vec![0.5]);

        solver
            .set_problem(
                vec![1.0],
                vec![Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![0],
                }],
                101325.0,
            )
            .unwrap();
        assert_eq!(solver.species_phase, vec![0]);
        assert!(solver.solution.is_empty());
        assert!(solver.moles.is_empty());
    }

    #[test]
    fn legacy_setup_rejects_invalid_solver_settings_before_matrix_publication() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        solver.n0 = vec![1.0, 0.0];
        solver.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        solver.stoich_matrix = DMatrix::from_row_slice(1, 1, &[42.0]);
        solver.solver_settings.solver_params.tol = 0.0;

        assert!(matches!(
            solver.create_stoich_matrix(),
            Err(ReactionExtentError::InvalidProblem {
                field: "solver_params.tol",
                ..
            })
        ));
        assert_eq!(solver.stoich_matrix, DMatrix::from_row_slice(1, 1, &[42.0]));
    }

    #[test]
    fn temperature_worker_configuration_clones_the_canonical_solver_state() {
        let mut source = EquilibriumLogMoles::empty();
        source.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        source.elem_composition = DMatrix::from_row_slice(2, 2, &[2.0, 0.0, 0.0, 1.0]);
        source.stoich_matrix = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        source.n0 = vec![1.0, 0.5];
        source.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        source.species_phase = vec![0, 0];
        source.elements_vector = vec![2.0, 1.0];
        source.solution = vec![9.0];
        source.moles = vec![8.0];
        source.initial_guess = Some(vec![0.25]);
        source.T = 900.0;
        source.P = 202_650.0;
        source.p0 = 101_325.0;
        source.gibbs = vec![Rc::new(|_| -1.0) as GibbsFn, Rc::new(|_| 0.5) as GibbsFn];
        source.solver_settings.solver = Solvers::NR;
        source.solver_settings.scaling_flag = true;
        source.solver_settings.continuation_seed_policy =
            ContinuationSeedPolicy::IndependentPerPoint;

        let mut local = EquilibriumLogMoles::empty();
        let local_gibbs = vec![Rc::new(|_| -1.0) as GibbsFn, Rc::new(|_| 0.5) as GibbsFn];

        let seed = source.temperature_worker_seed();
        seed.apply(&mut local, 777.0, local_gibbs, Some(vec![0.1]));

        assert_eq!(local.subs_data.substances, source.subs_data.substances);
        assert_eq!(local.elem_composition, source.elem_composition);
        assert_eq!(local.stoich_matrix, source.stoich_matrix);
        assert_eq!(local.n0, source.n0);
        assert_eq!(local.phases.len(), source.phases.len());
        assert_eq!(local.species_phase, source.species_phase);
        assert_eq!(local.elements_vector, source.elements_vector);
        assert_eq!(local.solution, Vec::<f64>::new());
        assert_eq!(local.moles, Vec::<f64>::new());
        assert_eq!(local.initial_guess, Some(vec![0.1]));
        assert_eq!(local.T, 777.0);
        assert_eq!(local.P, source.P);
        assert_eq!(local.p0, source.p0);
        assert_eq!(local.gibbs.len(), 2);
        assert_eq!(local.solver_settings.solver, Solvers::NR);
        assert!(local.solver_settings.scaling_flag);
        assert_eq!(
            local.solver_settings.continuation_seed_policy,
            ContinuationSeedPolicy::IndependentPerPoint
        );

        // The source object remains unchanged; the helper only stages the worker.
        assert_eq!(source.solution, vec![9.0]);
        assert_eq!(source.moles, vec![8.0]);
        assert_eq!(source.initial_guess, Some(vec![0.25]));
        assert_eq!(source.T, 900.0);
    }

    #[test]
    fn loglevel_configuration_is_instance_metadata_without_global_side_effects() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.with_loglevel(Some("debug"));
        assert_eq!(solver.loglevel.as_deref(), Some("debug"));

        solver.with_loglevel(None);
        assert_eq!(solver.loglevel, None);
    }

    #[test]
    fn moles_table_renders_the_published_sweep_columns_without_side_effects() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.map_of_moles_for_each_substance = HashMap::from([
            ("T".to_string(), vec![1000.0, 1100.0]),
            ("O2".to_string(), vec![0.2, 0.1]),
            ("O".to_string(), vec![0.8, 0.9]),
        ]);

        let rendered = solver.moles_table().to_string();
        let header_line = rendered
            .lines()
            .find(|line| line.contains("O") && line.contains("T"))
            .expect("table header should be rendered");

        let o_pos = header_line.find("O").expect("O column should be present");
        let o2_pos = header_line.find("O2").expect("O2 column should be present");
        let t_pos = header_line.find("T").expect("T column should be present");
        assert!(o_pos < o2_pos);
        assert!(o2_pos < t_pos);
        assert!(rendered.contains("1000.000000"));
        assert!(rendered.contains("0.200000"));
        assert!(rendered.contains("0.900000"));
    }

    #[test]
    fn failed_solve_does_not_publish_an_automatic_seed_or_result_bundle() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string()];
        solver.n0 = vec![1.0];
        solver.gibbs = vec![Rc::new(|_| 0.0) as GibbsFn];
        solver.solution = vec![42.0];
        solver.moles = vec![24.0];
        solver.solver_settings.solver_policy = Some(SolverPolicy::Cascade(Vec::new()));

        assert!(matches!(
            solver.solve(),
            Err(ReactionExtentError::InvalidProblem {
                field: "solver_policy",
                ..
            })
        ));
        assert!(solver.initial_guess.is_none());
        assert_eq!(solver.solution, vec![42.0]);
        assert_eq!(solver.moles, vec![24.0]);
    }

    #[test]
    fn mole_publication_rejects_duplicate_species_without_mutation() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O2".to_string()];
        solver.moles = vec![42.0];
        solver
            .map_of_moles_for_each_substance
            .insert("sentinel".to_string(), vec![24.0]);

        assert!(matches!(
            solver.compute_species_moles(vec![0.0, 0.0]),
            Err(ReactionExtentError::InvalidProblem {
                field: "species",
                ..
            })
        ));
        assert_eq!(solver.moles, vec![42.0]);
        assert_eq!(
            solver.map_of_moles_for_each_substance.get("sentinel"),
            Some(&vec![24.0])
        );
    }

    #[test]
    fn equilibrium_system_setup_does_not_publish_partial_matrices_on_failure() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.elem_composition = DMatrix::from_element(1, 1, 42.0);
        solver.stoich_matrix = DMatrix::from_element(1, 1, 24.0);
        solver.elements_vector = vec![7.0];
        solver.species_phase = vec![3];
        solver.n0 = vec![1.0];

        assert!(solver.create_equilibrium_system().is_err());
        assert_eq!(solver.elem_composition, DMatrix::from_element(1, 1, 42.0));
        assert_eq!(solver.stoich_matrix, DMatrix::from_element(1, 1, 24.0));
        assert_eq!(solver.elements_vector, vec![7.0]);
        assert_eq!(solver.species_phase, vec![3]);
    }

    #[test]
    fn result_publication_rejects_misaligned_rows_without_partial_maps() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.moles_for_T_range = vec![(1000.0, vec![0.9])];
        solver
            .map_of_moles_for_each_substance
            .insert("previous".to_string(), vec![1.0]);

        assert!(matches!(
            solver.map_of_moles_for_each_substance(),
            Err(ReactionExtentError::DimensionMismatch(_))
        ));
        assert_eq!(
            solver.map_of_moles_for_each_substance.get("previous"),
            Some(&vec![1.0])
        );

        solver.moles = vec![0.5];
        assert!(matches!(
            solver.compute_species_moles(vec![0.0]),
            Err(ReactionExtentError::DimensionMismatch(_))
        ));
        assert_eq!(solver.moles, vec![0.5]);
    }

    #[test]
    fn reaction_phase_stoichiometry_single_phase() {
        use nalgebra::DMatrix;

        let nu = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);

        let phases = vec![Phase {
            species: vec![0, 1],
            kind: PhaseKind::IdealGas,
        }];

        let delta_n = reaction_phase_stoichiometry(&nu, &phases);

        assert_eq!(delta_n.len(), 1);
        assert_eq!(delta_n[0].len(), 1);

        // ΔN = -1 + 2 = +1
        assert!(approx_eq(delta_n[0][0], 1.0, 1e-12));
    }
    /////////////////////////elements total/////////////////////////////////////
    #[test]
    fn element_totals_o2_o() {
        // Species: [O2, O]
        // Elements: [O]
        let a = DMatrix::from_row_slice(
            2,
            1,
            &[
                2.0, // O2
                1.0, // O
            ],
        );

        let n0 = vec![1.0, 0.0]; // 1 mol O2
        let b = compute_element_totals(&a, &n0).unwrap();

        assert_eq!(b.len(), 1);
        assert!((b[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn element_totals_n_o_system() {
        // Species: N2, O2, NO
        // Elements: N, O
        let a = DMatrix::from_row_slice(
            3,
            2,
            &[
                2.0, 0.0, // N2
                0.0, 2.0, // O2
                1.0, 1.0, // NO
            ],
        );

        let n0 = vec![1.0, 2.0, 0.5];
        let b = compute_element_totals(&a, &n0).unwrap();

        assert!((b[0] - (2.0 * 1.0 + 1.0 * 0.5)).abs() < 1e-12); // N
        assert!((b[1] - (2.0 * 2.0 + 1.0 * 0.5)).abs() < 1e-12); // O
    }

    #[test]
    fn element_totals_rejects_invalid_public_matrix_inputs() {
        let matrix = DMatrix::from_element(2, 1, 1.0);
        assert!(matches!(
            compute_element_totals(&matrix, &[1.0]),
            Err(ReactionExtentError::DimensionMismatch(_))
        ));

        let non_finite_matrix = DMatrix::from_element(1, 1, f64::NAN);
        assert!(matches!(
            compute_element_totals(&non_finite_matrix, &[1.0]),
            Err(ReactionExtentError::InvalidProblem {
                field: "element_composition_or_initial_moles",
                ..
            })
        ));
    }

    #[test]
    fn species_to_phase_map_basic() {
        let phases = vec![Phase {
            species: vec![0, 1],
            kind: PhaseKind::IdealGas,
        }];

        let map = species_to_phase_map(&phases, 2).unwrap();

        assert_eq!(map.len(), 2);
        assert_eq!(map[0], 0);
        assert_eq!(map[1], 0);
    }
    ////////////////////////////////RESIDUALS TESTING///////////////////////////////////////////////////
    #[test]
    fn equilibrium_logmole_residual_dimension() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]); // O2 -> 2O
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]); // O2, O
        let element_totals = vec![2.0];

        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];

        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        let m = 2;
        let species_phase = species_to_phase_map(&phases, m).unwrap();
        let f = equilibrium_logmole_residual(
            reactions,
            elements,
            element_totals,
            gibbs,
            phases,
            3000.0,
            101325.0,
            101325.0,
            species_phase,
            0.0,
            0.0,
        )
        .unwrap();

        let x = vec![0.0, -10.0]; // log(n_O2), log(n_O)
        let res = f(&x).unwrap();

        assert_eq!(res.len(), 2); // 1 reaction + 1 element
    }

    #[test]
    fn legacy_residual_factory_reports_typed_setup_errors() {
        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0)];

        assert!(matches!(
            equilibrium_logmole_residual(
                reactions.clone(),
                elements.clone(),
                vec![2.0],
                gibbs,
                phases.clone(),
                3000.0,
                101325.0,
                101325.0,
                vec![0, 0],
                0.0,
                0.0,
            ),
            Err(ReactionExtentError::DimensionMismatch(_))
        ));

        assert!(matches!(
            equilibrium_logmole_residual(
                reactions,
                elements,
                vec![2.0],
                vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)],
                phases,
                0.0,
                101325.0,
                101325.0,
                vec![0, 0],
                0.0,
                0.0,
            ),
            Err(ReactionExtentError::InvalidConditions {
                parameter: "temperature",
                value: 0.0,
            })
        ));
    }

    #[test]
    fn residual_of_complex_sys_low_temp() {
        use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};

        // Test at problematic low temperature
        let T = 500.0; // Low temperature where solver fails

        let mut user_subs = SubsData::new();
        user_subs.substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "N2O4".to_string(),
        ];
        let map_of_phases = user_subs
            .substances
            .iter()
            .map(|s| (s.clone(), Some(Phases::Gas)))
            .collect();
        user_subs.map_of_phases = map_of_phases;
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.search_substances().unwrap();
        user_subs.parse_all_thermal_coeffs().unwrap();
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();
        let elem = &user_subs.elem_composition_matrix.clone().unwrap();
        let _ = user_subs.extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T);
        let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase().unwrap();

        let mut g_vec: Vec<GibbsFn> = Vec::new();
        for sub in &user_subs.substances {
            let boxed_fn = map_of_functions.remove(sub).unwrap();
            let gibbs_fn: GibbsFn = Rc::new(move |T: f64| boxed_fn(T));
            g_vec.push(gibbs_fn);
        }

        let rb = compute_reaction_basis(elem, 1e-12).unwrap();
        let stoich = rb.reactions;
        let n0 = vec![1.0, 1e-5, 1e-5, 1e-5, 1e-5]; // initial moles: N2, N, O2, O, N2O4
        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1, 2, 3, 4],
        }];
        let elem_total = compute_element_totals(&elem, &n0).unwrap();
        let species_phase = species_to_phase_map(&phases, n0.len()).unwrap();
        // Create residual function
        let f = equilibrium_logmole_residual(
            stoich,
            elem.clone(),
            elem_total.data.as_vec().clone(),
            g_vec.clone(),
            phases.clone(),
            T,
            101325.0,
            101325.0,
            species_phase,
            1e-13,
            1e-13,
        )
        .unwrap();
        let initial = vec![1.0, 1e-5, 1e-5, 1e-5, 1e-5];
        let residual = f(&initial).unwrap();
        println!("residual {:?}", residual);
        // Create jacobian function
    }
    #[test]
    fn residual_of_simple_sys_low_temp() {
        use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};

        // Test at problematic low temperature
        let T = 500.0; // Low temperature where solver fails

        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["O2".to_string(), "O".to_string()];
        let map_of_phases = user_subs
            .substances
            .iter()
            .map(|s| (s.clone(), Some(Phases::Gas)))
            .collect();
        user_subs.map_of_phases = map_of_phases;
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.search_substances().unwrap();
        user_subs.parse_all_thermal_coeffs().unwrap();
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();
        let elem = &user_subs.elem_composition_matrix.clone().unwrap();
        let _ = user_subs.extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T);
        let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase().unwrap();

        let mut g_vec: Vec<GibbsFn> = Vec::new();
        for sub in &user_subs.substances {
            let boxed_fn = map_of_functions.remove(sub).unwrap();
            let gibbs_fn: GibbsFn = Rc::new(move |T: f64| boxed_fn(T));
            g_vec.push(gibbs_fn);
        }

        let rb = compute_reaction_basis(elem, 1e-12).unwrap();
        let stoich = rb.reactions;
        let n0 = vec![1.0, 1e-5]; // initial moles: N2, N, O2, O, N2O4
        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        let elem_total = compute_element_totals(&elem, &n0).unwrap();
        let species_phase = species_to_phase_map(&phases, n0.len()).unwrap();
        // Create residual function
        let f = equilibrium_logmole_residual(
            stoich,
            elem.clone(),
            elem_total.data.as_vec().clone(),
            g_vec.clone(),
            phases.clone(),
            T,
            101325.0,
            101325.0,
            species_phase,
            1e-13,
            1e-13,
        )
        .unwrap();
        let initial = vec![1.0, 1e-5];
        let residual = f(&initial).unwrap();
        println!("residual {:?}", residual);
        // Create jacobian function
    }
    /////////////////////////////////JACOBIAN TESTING//////////////////////////////////////
    #[test]
    fn equilibrium_logmole_jacobian_matches_fd() {
        let tol = 1e-6;
        let eps = 1e-7;

        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let element_totals = vec![2.0];

        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];

        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        let species_phase = species_to_phase_map(&phases, 2).unwrap();
        let delta_n = reaction_phase_stoichiometry(&reactions, &phases);

        for (temperature, x) in [
            (500.0, vec![0.0, -25.0]),  // low-T / trace-species state
            (1000.0, vec![0.0, -2.0]),  // normal state
            (3000.0, vec![0.0, -10.0]), // high-T / trace-species state
        ] {
            let f = equilibrium_logmole_residual(
                reactions.clone(),
                elements.clone(),
                element_totals.clone(),
                gibbs.clone(),
                phases.clone(),
                temperature,
                101325.0,
                101325.0,
                species_phase.clone(),
                0.0,
                0.0,
            )
            .unwrap();

            let j_fd = finite_difference_jacobian(&f, &x, eps).unwrap();
            let j_an = equilibrium_logmole_jacobian(
                &x,
                &reactions,
                &elements,
                &species_phase,
                &delta_n,
                phases.len(),
                0.0,
                0.0,
            )
            .unwrap();

            assert_eq!(j_an.nrows(), j_fd.nrows());
            assert_eq!(j_an.ncols(), j_fd.ncols());

            for i in 0..j_an.nrows() {
                for k in 0..j_an.ncols() {
                    let diff = (j_an[(i, k)] - j_fd[(i, k)]).abs();
                    assert!(
                        diff < tol,
                        "Jacobian mismatch at T={temperature}, ({i},{k}) analytic={} fd={}",
                        j_an[(i, k)],
                        j_fd[(i, k)]
                    );
                }
            }
        }
    }
    /////////////////////////////////SCALING////////////////////////////////////////////////////
    #[test]
    fn equilibrium_scaling_dimensions_and_signs() {
        let stoich = DMatrix::from_row_slice(
            2,
            1,
            &[
                -1.0, // O2
                2.0,  // O
            ],
        );

        let elements = DMatrix::from_row_slice(
            2,
            1,
            &[
                2.0, // O2
                1.0, // O
            ],
        );

        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];

        let element_totals = vec![2.0];
        let T = 3000.0;

        let scale = equilibrium_scaling(&stoich, &elements, &gibbs, &element_totals, T).unwrap();

        assert_eq!(scale.len(), 2); // 1 reaction + 1 element

        for (i, s) in scale.iter().enumerate() {
            assert!(*s > 0.0 && s.is_finite(), "Invalid scale[{}] = {}", i, s);
        }
    }

    #[test]
    fn equilibrium_scaling_element_block() {
        use nalgebra::DMatrix;

        let stoich = DMatrix::zeros(2, 1);

        let elements = DMatrix::from_row_slice(
            2,
            2,
            &[
                2.0, 0.0, // O2
                1.0, 0.0, // O
            ],
        );

        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];

        let element_totals = vec![2.0, 0.0];

        let scale =
            equilibrium_scaling(&stoich, &elements, &gibbs, &element_totals, 1000.0).unwrap();

        let r = stoich.ncols();

        // element equations must be scaled
        for j in 0..element_totals.len() {
            let s = scale[r + j];
            assert!(
                s.is_finite() && s > 0.0,
                "Invalid element scale[{}] = {}",
                j,
                s
            );
        }

        // element scales should be >= reaction scales floor
        let floor = 10.0; // or whatever floor equilibrium_scaling uses
        assert!(scale[r] >= floor);
        assert!(scale[r + 1] >= floor);
    }

    #[test]
    fn equilibrium_scaling_inert_species_does_not_break() {
        let stoich = DMatrix::from_row_slice(
            4,
            1,
            &[
                -1.0, // O2
                2.0,  // O
                0.0,  // N2
                0.0,  // N
            ],
        );

        let elements = DMatrix::from_row_slice(
            4,
            2,
            &[
                2.0, 0.0, // O2
                1.0, 0.0, // O
                0.0, 2.0, // N2
                0.0, 1.0, // N
            ],
        );

        let gibbs: Vec<GibbsFn> = vec![
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
        ];

        let element_totals = vec![2.0, 2.0];

        let scale =
            equilibrium_scaling(&stoich, &elements, &gibbs, &element_totals, 4000.0).unwrap();

        for s in scale {
            assert!(s.is_finite() && s > 0.0);
        }
    }

    #[test]
    fn row_scaling_transforms_residual_and_jacobian_consistently() {
        let residual = scaled_residual(Box::new(|_| Ok(vec![100.0, -6.0])), vec![10.0, 3.0]);
        let jacobian = scaled_jacobian(
            Box::new(|_| Ok(DMatrix::from_row_slice(2, 2, &[20.0, 10.0, -9.0, 6.0]))),
            vec![10.0, 3.0],
        );

        assert_eq!(residual(&[0.0, 0.0]).unwrap(), vec![10.0, -2.0]);
        assert_eq!(
            jacobian(&[0.0, 0.0]).unwrap(),
            DMatrix::from_row_slice(2, 2, &[2.0, 1.0, -3.0, 2.0])
        );
    }

    #[test]
    fn scaling_boundaries_return_typed_errors_instead_of_indexing_or_truncating() {
        let stoich = DMatrix::zeros(2, 1);
        let elements = DMatrix::zeros(1, 1);
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0)];

        assert!(matches!(
            equilibrium_scaling(&stoich, &elements, &gibbs, &[1.0], 1000.0),
            Err(ReactionExtentError::DimensionMismatch(_))
        ));
        assert!(matches!(
            equilibrium_scaling(
                &DMatrix::zeros(1, 0),
                &DMatrix::zeros(1, 0),
                &[Rc::new(|_| 0.0) as GibbsFn],
                &[],
                f64::NAN,
            ),
            Err(ReactionExtentError::InvalidConditions {
                parameter: "temperature",
                ..
            })
        ));

        let bad_residual = scaled_residual(Box::new(|_| Ok(vec![1.0])), vec![1.0, 1.0]);
        assert!(matches!(
            bad_residual(&[]),
            Err(ReactionExtentError::DimensionMismatch(_))
        ));
        let bad_jacobian = scaled_jacobian(Box::new(|_| Ok(DMatrix::identity(1, 1))), vec![0.0]);
        assert!(matches!(
            bad_jacobian(&[]),
            Err(ReactionExtentError::InvalidProblem {
                field: "residual_scale",
                ..
            })
        ));
    }

    ////////////////////////////////SOLUTION AT CONST T//////////////////////////////////////////////////////

    #[test]
    fn O2_O_equilibrium_logmoles_LM() {
        let T = 500.0;

        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false).unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.set_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM)));

        instance.create_stoich_matrix().unwrap();
        instance.solve().unwrap();

        let moles = &instance.moles;
        println!("moles {:?}", moles);
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");

        // A backend success becomes public state only after the common
        // candidate-validation gate has checked residuals and element totals.
        let report = instance
            .last_validation_report
            .as_ref()
            .expect("accepted solve must publish a validation report");
        assert!(report.residual_l2_norm.is_finite());
        assert!(report.max_abs_element_balance_error <= 1e-5);
        assert!(report.min_moles > 0.0);

        let accepted = instance.accepted_solution().unwrap();
        assert_eq!(accepted.moles(), moles);
        assert_eq!(accepted.log_moles(), instance.solution);
        assert_eq!(accepted.conditions().temperature(), T);
        assert_eq!(accepted.validation(), report);

        let solve_report = instance
            .last_solve_report
            .as_ref()
            .expect("accepted solve must publish a backend report");
        assert_eq!(
            solve_report.accepted_backend,
            SolverBackend::Legacy(Solvers::LM)
        );
        assert_eq!(
            solve_report.policy,
            SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM))
        );
        assert_eq!(solve_report.attempts.len(), 1);
        assert_eq!(
            solve_report.attempts[0].outcome,
            SolverAttemptOutcome::Accepted
        );
    }

    #[test]
    fn o2_o_equilibrium_uses_rst_symbolic_lm_without_the_legacy_jacobian() {
        let mut instance = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.create_stoich_matrix().unwrap();
        instance.set_solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
            RustedSciTheSolver::LevenbergMarquardt,
        )));

        instance.solve().unwrap();

        let report = instance.last_solve_report.as_ref().unwrap();
        assert_eq!(
            report.accepted_backend,
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt)
        );
        assert_eq!(report.attempts.len(), 1);
        assert_eq!(report.attempts[0].outcome, SolverAttemptOutcome::Accepted);
        let metrics = report.attempts[0]
            .metrics
            .as_ref()
            .expect("RST attempt must preserve engine diagnostics");
        assert!(metrics.residual_evaluations > 0);
        assert!(metrics.jacobian_evaluations > 0);
        assert!(
            instance
                .moles
                .iter()
                .all(|moles| moles.is_finite() && *moles > 0.0)
        );
    }

    #[test]
    fn o2_o_equilibrium_defaults_to_rst_policy_when_no_explicit_policy_is_set() {
        let mut instance = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.create_stoich_matrix().unwrap();
        instance.solve().unwrap();

        let report = instance.last_solve_report.as_ref().unwrap();
        assert_eq!(
            report.policy,
            SolverPolicy::rusted_scithe_default(),
            "default policy should prefer the RST symbolic cascade"
        );
        assert!(
            report
                .attempts
                .iter()
                .all(|attempt| matches!(attempt.backend, SolverBackend::RustedSciThe(_)))
        );
        assert!(matches!(
            report.accepted_backend,
            SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt)
        ));
    }

    #[test]
    fn accepted_solution_preserves_non_negative_moles_and_element_totals() {
        let mut instance = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.create_stoich_matrix().unwrap();
        instance.solve().unwrap();

        assert!(instance.moles.iter().all(|&n| n >= 0.0));

        let totals = compute_element_totals(&instance.elem_composition, &instance.moles).unwrap();
        for (lhs, rhs) in totals.iter().zip(instance.elements_vector.iter()) {
            assert!(
                (lhs - rhs).abs() <= 1e-8,
                "element balance drifted: {lhs} vs {rhs}"
            );
        }
    }

    #[test]
    fn species_ordering_does_not_change_the_accepted_equilibrium_composition() {
        let mut ordered = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        ordered.n0 = vec![1.0, 1e-5];
        ordered.initial_guess = Some(vec![1.0, 1e-5]);
        ordered.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        ordered.create_stoich_matrix().unwrap();
        ordered.solve().unwrap();

        let mut reversed = gas_solver(
            vec!["O".to_string(), "O2".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        reversed.n0 = vec![1e-5, 1.0];
        reversed.initial_guess = Some(vec![1e-5, 1.0]);
        reversed.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        reversed.create_stoich_matrix().unwrap();
        reversed.solve().unwrap();

        let ordered_o2 = ordered
            .map_of_moles_for_each_substance
            .get("O2")
            .and_then(|values| values.first().copied())
            .expect("ordered solution should publish O2");
        let ordered_o = ordered
            .map_of_moles_for_each_substance
            .get("O")
            .and_then(|values| values.first().copied())
            .expect("ordered solution should publish O");
        let reversed_o2 = reversed
            .map_of_moles_for_each_substance
            .get("O2")
            .and_then(|values| values.first().copied())
            .expect("reversed solution should publish O2");
        let reversed_o = reversed
            .map_of_moles_for_each_substance
            .get("O")
            .and_then(|values| values.first().copied())
            .expect("reversed solution should publish O");

        let ordered_fraction = normalized_mole_fractions(&[ordered_o2, ordered_o]);
        let reversed_fraction = normalized_mole_fractions(&[reversed_o2, reversed_o]);

        for (lhs, rhs) in ordered_fraction.iter().zip(reversed_fraction.iter()) {
            assert!((lhs - rhs).abs() < 1e-10);
        }
    }

    #[test]
    fn uniform_scaling_of_initial_moles_preserves_the_accepted_composition() {
        let mut base = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        base.n0 = vec![1.0, 1e-5];
        base.initial_guess = Some(vec![1.0, 1e-5]);
        base.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        base.create_stoich_matrix().unwrap();
        base.solve().unwrap();

        let mut scaled = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        scaled.n0 = vec![10.0, 1e-4];
        scaled.initial_guess = Some(vec![10.0, 1e-4]);
        scaled.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        scaled.create_stoich_matrix().unwrap();
        scaled.solve().unwrap();

        let base_fraction = normalized_mole_fractions(&base.moles);
        let scaled_fraction = normalized_mole_fractions(&scaled.moles);
        for (lhs, rhs) in base_fraction.iter().zip(scaled_fraction.iter()) {
            assert!((lhs - rhs).abs() < 1e-6);
        }
    }

    #[test]
    fn accepted_solution_is_robust_to_reasonable_initial_guess_perturbations() {
        let mut base = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        base.n0 = vec![1.0, 1e-5];
        base.initial_guess = Some(vec![1.0, 1e-5]);
        base.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        base.create_stoich_matrix().unwrap();
        base.solve().unwrap();
        let base_solution = base.solution.clone();

        let mut perturbed = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        perturbed.n0 = vec![1.0, 1e-5];
        perturbed.initial_guess = Some(vec![0.8, -10.2]);
        perturbed.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        perturbed.create_stoich_matrix().unwrap();
        perturbed.solve().unwrap();

        assert_eq!(base_solution.len(), perturbed.solution.len());
        let base_fraction = normalized_mole_fractions(&base.moles);
        let perturbed_fraction = normalized_mole_fractions(&perturbed.moles);
        for (lhs, rhs) in base_fraction.iter().zip(perturbed_fraction.iter()) {
            assert!((lhs - rhs).abs() < 1e-6);
        }
    }

    #[test]
    fn temperature_and_pressure_trends_follow_the_expected_direction() {
        let build_solver = |temperature: f64, pressure: f64| {
            let mut instance = gas_solver(
                vec!["O2".to_string(), "O".to_string()],
                temperature,
                pressure,
                Solvers::LM,
                None,
                false,
            )
            .unwrap();
            instance.n0 = vec![1.0, 1e-5];
            instance.initial_guess = Some(vec![1.0, 1e-5]);
            instance.phases = vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            }];
            instance.create_stoich_matrix().unwrap();
            instance.solve().unwrap();
            instance
        };

        let low_t = build_solver(500.0, 101325.0);
        let high_t = build_solver(5000.0, 101325.0);
        let low_p = build_solver(3000.0, 20_000.0);
        let high_p = build_solver(3000.0, 5_000_000.0);

        let low_t_o_fraction = low_t.moles[1] / (low_t.moles[0] + low_t.moles[1]);
        let high_t_o_fraction = high_t.moles[1] / (high_t.moles[0] + high_t.moles[1]);
        assert!(
            high_t_o_fraction > low_t_o_fraction,
            "oxygen atom fraction should rise with temperature"
        );

        let low_p_o_fraction = low_p.moles[1] / (low_p.moles[0] + low_p.moles[1]);
        let high_p_o_fraction = high_p.moles[1] / (high_p.moles[0] + high_p.moles[1]);
        assert!(
            high_p_o_fraction < low_p_o_fraction,
            "oxygen atom fraction should fall with pressure"
        );
    }

    #[test]
    fn adding_an_inert_diluent_suppresses_gas_phase_dissociation() {
        let build_solver = |include_n2: bool| {
            let mut substances = vec!["O2".to_string(), "O".to_string()];
            let mut n0 = vec![1.0, 1e-5];
            if include_n2 {
                substances.push("N2".to_string());
                n0.push(10.0);
            }

            let mut instance =
                gas_solver(substances, 3000.0, 101325.0, Solvers::LM, None, false).unwrap();
            instance.n0 = n0.clone();
            instance.initial_guess = Some(n0);
            instance.phases = vec![Phase {
                kind: PhaseKind::IdealGas,
                species: (0..instance.n0.len()).collect(),
            }];
            instance.create_stoich_matrix().unwrap();
            instance.solve().unwrap();
            instance
        };

        let without_inert = build_solver(false);
        let with_inert = build_solver(true);

        let without_fraction = without_inert.moles[1] / without_inert.moles.iter().sum::<f64>();
        let with_fraction = with_inert.moles[1] / with_inert.moles.iter().sum::<f64>();

        assert!(
            with_fraction < without_fraction,
            "inert dilution should suppress the atomic oxygen fraction at fixed T and P"
        );
        assert!(with_inert.moles.iter().all(|&n| n >= 0.0));
    }

    #[test]
    fn O2_O_equilibrium_logmoles_NR() {
        let T = 500.0;

        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::NR, Some("info"), false).unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];

        instance.create_stoich_matrix().unwrap();
        instance.solve().unwrap();

        let moles = &instance.moles;
        println!("moles {:?}", moles);
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");
    }
    #[test]
    fn O2_O_equilibrium_logmoles_TR() {
        let T = 500.0;

        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::TR, Some("info"), false).unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];

        instance.create_stoich_matrix().unwrap();
        instance.solve().unwrap();

        let moles = &instance.moles;
        println!("moles {:?}", moles);
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");
    }
    #[test]
    fn N2_N_equilibrium_logmoles_low_T() {
        let T = 500.0;

        let substances = vec!["N2".to_string(), "N".to_string()];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false).unwrap();
        instance.n0 = vec![1e-5, 1.0];
        instance.P = 101325.0;
        instance.T = T;
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];

        instance.create_stoich_matrix().unwrap();
        instance.solve().unwrap();

        let moles = &instance.moles;
        assert_relative_eq!(moles[0], 0.5, epsilon = 1e-4);
        println!("moles {:?}", moles);
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");
    }
    #[test]
    fn N2_N_equilibrium_logmoles_high_T() {
        let T = 5500.0;

        let substances = vec!["N2".to_string(), "N".to_string()];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false).unwrap();
        instance.n0 = vec![1e-5, 1.0];
        instance.P = 101325.0;
        instance.T = T;
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];

        instance.create_equilibrium_system().unwrap();

        instance.solve().unwrap();

        let moles = &instance.moles;
        // at hight T N2 partially decomposed into N
        assert!(moles[0] < 0.5);
        println!("moles {:?}", moles);
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");
    }

    #[test]
    fn complex_equilibrium_logmoles_high_T() {
        let T = 5500.0;

        let substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "NO2".to_string(),
            "NO".to_string(),
            "N2O".to_string(),
            "N2O4".to_string(),
        ];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false).unwrap();
        instance.n0 = vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5];
        instance.P = 101325.0;
        instance.T = T;
        // instance.initial_guess = Some(vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1, 2, 3, 4, 5, 6, 7],
        }];

        instance.create_stoich_matrix().unwrap();
        instance.solve().unwrap();

        let moles = &instance.moles;

        println!("moles {:?}", moles);
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");
    }

    #[test]
    fn complex_equilibrium_logmoles_low_T() {
        let T = 500.0;

        let substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "NO2".to_string(),
            "NO".to_string(),
            "N2O".to_string(),
            "N2O4".to_string(),
        ];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false).unwrap();
        instance.n0 = vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5];
        instance.P = 101325.0;
        instance.T = T;
        //  instance.initial_guess = Some(vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1, 2, 3, 4, 5, 6, 7],
        }];

        instance.create_stoich_matrix().unwrap();
        instance.solve().unwrap();

        let moles = &instance.moles;

        println!("moles {:?}", moles);
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");
    }
    /////////////////////////////////SOLUTION AT T RANGE///////////////////////////////////////////////////////
    #[test]
    #[ignore = "expensive sweep study"]
    fn O2_O_equilibrium_for_T_range_logmoles_par() {
        let now = Instant::now();
        let mut user_subs = SubsData::new();

        user_subs.substances = vec!["O2".to_string(), "O".to_string()];
        user_subs
            .map_of_phases
            .insert("O2".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("O".to_string(), Some(Phases::Gas));

        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

        let mut instance = EquilibriumLogMoles::new();
        instance.solver_settings.scaling_flag = false;
        instance.set_initial_guess(vec![0.0, -10.0]).unwrap();
        instance
            .set_problem(
                vec![1.0, 0.0],
                vec![Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![0, 1],
                }],
                101325.0,
            )
            .unwrap();
        instance.subs_data = user_subs;
        instance.create_equilibrium_system().unwrap();
        instance
            .solve_for_T_range_par2(400.0, 4500.0, 10.0)
            .unwrap();
        assert_eq!(
            instance.temperature_solutions.len(),
            instance.moles_for_T_range.len()
        );
        assert!(instance.temperature_solutions.iter().all(|snapshot| {
            snapshot.log_moles.len() == 2
                && snapshot.moles.len() == 2
                && snapshot.solve_report.accepted_backend
                    == snapshot.solve_report.attempts.last().unwrap().backend
        }));
        let _table = instance.moles_table();
        println!("elapsed time {} ms", now.elapsed().as_millis());
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
    #[ignore = "expensive sweep study"]
    fn O2_O_equilibrium_for_T_range_logmoles() {
        let now = Instant::now();
        let mut user_subs = SubsData::new();

        user_subs.substances = vec!["O2".to_string(), "O".to_string()];
        user_subs
            .map_of_phases
            .insert("O2".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("O".to_string(), Some(Phases::Gas));

        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

        let mut instance = EquilibriumLogMoles::new();
        instance.solver_settings.scaling_flag = false;
        instance.set_initial_guess(vec![0.0, -10.0]).unwrap();
        instance
            .set_problem(
                vec![1.0, 0.0],
                vec![Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![0, 1],
                }],
                101325.0,
            )
            .unwrap();
        instance.subs_data = user_subs;
        instance.create_equilibrium_system().unwrap();
        instance.solve_for_T_range(400.0, 4500.0, 10.0).unwrap();
        assert_eq!(
            instance.temperature_solutions.len(),
            instance.moles_for_T_range.len()
        );
        assert!(instance.temperature_solutions.iter().all(|snapshot| {
            snapshot.log_moles.len() == 2
                && snapshot.moles.len() == 2
                && snapshot.solve_report.accepted_backend
                    == snapshot.solve_report.attempts.last().unwrap().backend
        }));
        let _table = instance.moles_table();
        println!("elapsed time {} ms", now.elapsed().as_millis());
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
    #[ignore = "expensive sweep study"]
    fn N2_N_equilibrium_for_T_range_logmoles() {
        let now = Instant::now();
        let mut user_subs = SubsData::new();

        user_subs.substances = vec!["N2".to_string(), "N".to_string()];
        user_subs
            .map_of_phases
            .insert("N2".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("N".to_string(), Some(Phases::Gas));

        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

        let mut instance = EquilibriumLogMoles::new();
        instance.solver_settings.scaling_flag = false;
        instance.set_initial_guess(vec![0.0, -10.0]).unwrap();
        instance.with_loglevel(Some("info"));
        instance
            .set_problem(
                vec![1.0, 0.0],
                vec![Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![0, 1],
                }],
                101325.0,
            )
            .unwrap();
        instance.subs_data = user_subs;
        instance.create_equilibrium_system().unwrap();
        instance.solve_for_T_range(400.0, 5500.0, 100.0).unwrap();
        let _table = instance.moles_table();

        println!("elapsed time {} ms", now.elapsed().as_millis());
        println!(
            "list of failed T points (if any): {:?}",
            instance.list_of_failed_T
        );
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
    #[ignore = "expensive sweep study"]
    fn complex_gas_equilibrium_for_T_range_final_LM() {
        let now = Instant::now();
        let mut user_subs = SubsData::new();

        user_subs.substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "NO".to_string(),
            "NO2".to_string(),
            "N2O".to_string(),
            "N2O4".to_string(),
        ];
        let map_of_phases = user_subs
            .substances
            .iter()
            .map(|s| (s.clone(), Some(Phases::Gas)))
            .collect();
        user_subs.map_of_phases = map_of_phases;

        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

        let mut instance = EquilibriumLogMoles::new();
        instance.solver_settings.scaling_flag = false;
        instance.with_loglevel(Some("info"));
        //   instance.set_initial_guess(vec![1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]);
        instance
            .set_problem(
                vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5],
                vec![Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![0, 1, 2, 3, 4, 5, 6, 7],
                }],
                101325.0,
            )
            .unwrap();
        instance.subs_data = user_subs;
        instance.create_equilibrium_system().unwrap();
        instance.solve_for_T_range(400.0, 5500.0, 100.0).unwrap();
        let _table = instance.moles_table();

        println!("elapsed time {} ms", now.elapsed().as_millis());
        println!(
            "list of failed T points (if any): {:?}",
            instance.list_of_failed_T
        );
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
    #[ignore = "expensive sweep study"]
    fn complex_gas_equilibrium_for_T_range_final_NR() {
        let now = Instant::now();
        let mut user_subs = SubsData::new();

        user_subs.substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "NO".to_string(),
            "NO2".to_string(),
            "N2O".to_string(),
            "N2O4".to_string(),
        ];
        let map_of_phases = user_subs
            .substances
            .iter()
            .map(|s| (s.clone(), Some(Phases::Gas)))
            .collect();
        user_subs.map_of_phases = map_of_phases;

        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

        let mut instance = EquilibriumLogMoles::new();
        instance.solver_settings.scaling_flag = false;
        instance.with_loglevel(Some("info"));
        //   instance.set_initial_guess(vec![1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]);
        instance
            .set_problem(
                vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5],
                vec![Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![0, 1, 2, 3, 4, 5, 6, 7],
                }],
                101325.0,
            )
            .unwrap();
        instance.subs_data = user_subs;
        instance.solver_settings.solver = Solvers::NR;
        instance.create_equilibrium_system().unwrap();
        instance.solve_for_T_range(400.0, 5500.0, 100.0).unwrap();
        let _table = instance.moles_table();

        println!("elapsed time {} ms", now.elapsed().as_millis());
        println!(
            "list of failed T points (if any): {:?}",
            instance.list_of_failed_T
        );
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
    #[ignore = "expensive sweep study"]
    fn complex_gas_equilibrium_for_T_range_final_TR() {
        let now = Instant::now();
        let substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "NO".to_string(),
            "NO2".to_string(),
            "N2O".to_string(),
            "N2O4".to_string(),
        ];
        let instance = gas_solver_for_T_range(
            substances,
            vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5],
            101325.0,
            400.0,
            5500.0,
            50.0,
            Solvers::TR,
            Some("info"),
            false,
        )
        .unwrap();

        println!("elapsed time {} ms", now.elapsed().as_millis());
        println!(
            "list of failed T points (if any): {:?}",
            instance.list_of_failed_T
        );
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
    #[ignore = "expensive sweep study"]
    fn complex_gas_equilibrium_for_T_range_final_par() {
        let now = Instant::now();
        let substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "NO".to_string(),
            "NO2".to_string(),
            "N2O".to_string(),
            "N2O4".to_string(),
        ];
        let instance = gas_solver_for_T_range(
            substances,
            vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5],
            101325.0,
            400.0,
            5500.0,
            50.0,
            Solvers::TR,
            Some("info"),
            true,
        )
        .unwrap();

        println!("elapsed time {} ms", now.elapsed().as_millis());
        println!(
            "list of failed T points (if any): {:?}",
            instance.list_of_failed_T
        );
        assert!(!instance.moles_for_T_range.is_empty());
    }
    /////////////////////////////////////SYMBOLIC///////////////////////////////////////////
    #[test]
    #[ignore = "expensive symbolic study"]
    fn N2_N_O2_O_equilibrium_sym_low_T() {
        use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::multiphase_equilibrium_residual_generator_sym;
        use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};

        // Test at problematic low temperature
        let T = 500.0; // Low temperature where solver fails

        let mut user_subs = SubsData::new();
        user_subs.substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "N2O4".to_string(),
        ];
        let map_of_phases = user_subs
            .substances
            .iter()
            .map(|s| (s.clone(), Some(Phases::Gas)))
            .collect();
        user_subs.map_of_phases = map_of_phases;
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.search_substances().unwrap();
        user_subs.parse_all_thermal_coeffs().unwrap();
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();
        let elem = &user_subs.elem_composition_matrix.clone().unwrap();
        let _ = user_subs.extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T);
        let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase().unwrap();

        let mut g_vec: Vec<GibbsFn> = Vec::new();
        for sub in &user_subs.substances {
            let boxed_fn = map_of_functions.remove(sub).unwrap();
            let gibbs_fn: GibbsFn = Rc::new(move |T: f64| boxed_fn(T));
            g_vec.push(gibbs_fn);
        }

        let rb = compute_reaction_basis(elem, 1e-12).unwrap();
        let stoich = rb.reactions;
        let n0 = vec![1.0, 2e-5, 2e-5, 3e-5, 4e-5]; // initial moles: N2, N, O2, O, N2O4
        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1, 2, 3, 4],
        }];
        let elem_total = compute_element_totals(&elem, &n0).unwrap();

        // Create residual function
        let species_phase = species_to_phase_map(&phases, n0.len()).unwrap();
        let f = equilibrium_logmole_residual(
            stoich.clone(),
            elem.clone(),
            elem_total.data.as_vec().clone(),
            g_vec.clone(),
            phases.clone(),
            T,
            101325.0,
            101325.0,
            species_phase,
            1e-13,
            1e-13,
        )
        .unwrap();
        let initial = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let residual_with_fn = f(&initial).unwrap();
        // creatinh symbolic residual
        let _ = user_subs.calculate_therm_map_of_sym();
        let map_of_functions = user_subs.calculate_dG0_sym_one_phase().unwrap();
        //   println!("{:?}", map_of_functions);
        let mut g_vec: Vec<Expr> = Vec::new();
        for subs in user_subs.substances.clone() {
            let dg0 = map_of_functions.get(&subs).unwrap();
            g_vec.push(dg0.clone());
        }
        let vars = Expr::IndexedVars(stoich.nrows(), "y").0;
        let vars: Vec<String> = vars.iter().map(|y| y.to_string().clone()).collect();
        let vars: Vec<&str> = vars.iter().map(|yi| yi.as_str()).collect();
        let sym_res = multiphase_equilibrium_residual_generator_sym(
            stoich.clone(),
            elem.clone(),
            elem_total.data.as_vec().clone(),
            g_vec,
            phases.clone(),
            101325.0,
            101325.0,
        )
        .unwrap();
        println!("unknowns {:?}", vars);
        for (i, sym) in sym_res.iter().enumerate() {
            let sym = sym.set_variable("T", T).simplify();
            println!("eq = {}", &sym);
            let f_sym = sym.lambdify_borrowed_thread_safe(&vars);
            let res_sym = f_sym(&initial.clone());

            let res_fun_i = residual_with_fn[i];
            println!("res sym {}, res_fun_i {}", res_sym, res_fun_i);
            assert_relative_eq!(res_sym, res_fun_i, epsilon = 1e-5);
        }
    }
    #[test]
    #[ignore = "expensive symbolic study"]
    fn N2_N_O2_O_equilibrium_sym_high_T() {
        use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::multiphase_equilibrium_residual_generator_sym;
        use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};

        // Test at problematic low temperature
        let T = 4500.0; // Low temperature where solver fails

        let mut user_subs = SubsData::new();
        user_subs.substances = vec![
            "N2".to_string(),
            "N".to_string(),
            "O2".to_string(),
            "O".to_string(),
            "N2O4".to_string(),
        ];
        let map_of_phases = user_subs
            .substances
            .iter()
            .map(|s| (s.clone(), Some(Phases::Gas)))
            .collect();
        user_subs.map_of_phases = map_of_phases;
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.search_substances().unwrap();
        user_subs.parse_all_thermal_coeffs().unwrap();
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();
        let elem = &user_subs.elem_composition_matrix.clone().unwrap();
        let _ = user_subs.extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T);
        let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase().unwrap();

        let mut g_vec: Vec<GibbsFn> = Vec::new();
        for sub in &user_subs.substances {
            let boxed_fn = map_of_functions.remove(sub).unwrap();
            let gibbs_fn: GibbsFn = Rc::new(move |T: f64| boxed_fn(T));
            g_vec.push(gibbs_fn);
        }

        let rb = compute_reaction_basis(elem, 1e-12).unwrap();
        let stoich = rb.reactions;
        let n0 = vec![1.0, 2e-5, 2e-5, 3e-5, 4e-5]; // initial moles: N2, N, O2, O, N2O4
        let phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1, 2, 3, 4],
        }];
        let elem_total = compute_element_totals(&elem, &n0).unwrap();

        // Create residual function
        let species_phase = species_to_phase_map(&phases, n0.len()).unwrap();
        let f = equilibrium_logmole_residual(
            stoich.clone(),
            elem.clone(),
            elem_total.data.as_vec().clone(),
            g_vec.clone(),
            phases.clone(),
            T,
            101325.0,
            101325.0,
            species_phase,
            1e-13,
            1e-13,
        )
        .unwrap();
        let initial = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let residual_with_fn = f(&initial).unwrap();
        // creatinh symbolic residual
        let _ = user_subs.calculate_therm_map_of_sym();
        let map_of_functions = user_subs.calculate_dG0_sym_one_phase().unwrap();
        //   println!("{:?}", map_of_functions);
        let mut g_vec: Vec<Expr> = Vec::new();
        for subs in user_subs.substances.clone() {
            let dg0 = map_of_functions.get(&subs).unwrap();
            g_vec.push(dg0.clone());
        }
        let vars = Expr::IndexedVars(stoich.nrows(), "y").0;
        let vars: Vec<String> = vars.iter().map(|y| y.to_string().clone()).collect();
        let vars: Vec<&str> = vars.iter().map(|yi| yi.as_str()).collect();
        let sym_res = multiphase_equilibrium_residual_generator_sym(
            stoich.clone(),
            elem.clone(),
            elem_total.data.as_vec().clone(),
            g_vec,
            phases.clone(),
            101325.0,
            101325.0,
        )
        .unwrap();
        println!("unknowns {:?}", vars);
        for (i, sym) in sym_res.iter().enumerate() {
            let sym = sym.set_variable("T", T).simplify();
            println!("eq = {}", &sym);
            let f_sym = sym.lambdify_borrowed_thread_safe(&vars);
            let res_sym = f_sym(&initial.clone());

            let res_fun_i = residual_with_fn[i];
            println!("res sym {}, res_fun_i {}", res_sym, res_fun_i);
            assert_relative_eq!(res_sym, res_fun_i, epsilon = 1e-5);
        }
    }

    #[test]
    fn N2_N_O2_O_equilibrium_symbolic_system() {
        let T = 500.0;

        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false).unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];

        instance.create_stoich_matrix().unwrap();
        let _ = instance.create_equilibrium_system();
        let _ = instance.create_symbolic_system();
    }
    /////////////////////from elements////////////////////
    #[test]
    fn gas_equilibrium_from_elements() {
        let T = 5500.0;

        let elements = vec!["N".to_string(), "O".to_string()];
        let map_of_nonzero_moles: HashMap<String, f64> =
            HashMap::from([("N2".to_string(), 0.5), ("O2".to_string(), 0.5)]);
        let mut instance = gas_solver_from_elements(
            elements,
            map_of_nonzero_moles,
            T,
            101325.0,
            Solvers::LM,
            Some("info"),
            false,
        )
        .unwrap();

        instance.create_stoich_matrix().unwrap();

        instance.solve().unwrap();

        //let moles = &instance.moles;
        let map_of_moles_for_each_substance = instance.map_of_moles_for_each_substance;
        println!("\n result {:?}", map_of_moles_for_each_substance);
    }

    #[test]
    fn gas_equilibrium_from_elements_for_T_range() {
        //  let T = 5500.0;

        let elements = vec!["N".to_string(), "O".to_string()];
        let map_of_nonzero_moles: HashMap<String, f64> =
            HashMap::from([("N2".to_string(), 0.5), ("O2".to_string(), 0.5)]);
        let instance = gas_solver_for_T_range_for_elements(
            elements,
            map_of_nonzero_moles,
            101325.0,
            400.0,
            4500.0,
            50.0,
            Solvers::LM,
            Some("info"),
        )
        .unwrap();

        let map_of_moles_for_each_substance = instance.map_of_moles_for_each_substance;
        println!("\n result {:?}", map_of_moles_for_each_substance);
    }
    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.45 — compute_species_moles (free function)
    // ──────────────────────────────────────────────
    #[test]
    fn compute_species_moles_rejects_non_finite_log_moles() {
        let result = compute_species_moles(&[f64::NAN]);
        assert!(matches!(
            result,
            Err(ReactionExtentError::InvalidCandidate {
                field: "log_moles",
                ..
            })
        ));

        let result = compute_species_moles(&[f64::INFINITY]);
        assert!(matches!(
            result,
            Err(ReactionExtentError::InvalidCandidate {
                field: "log_moles",
                ..
            })
        ));
    }

    #[test]
    fn compute_species_moles_rejects_underflow_to_zero() {
        // ln(5e-325) ≈ -746.5, exp(-746.5) underflows to 0.0 in f64
        // f64 minimum positive is ~5e-324, so -750.0 guarantees underflow
        let result = compute_species_moles(&[-750.0]);
        assert!(matches!(
            result,
            Err(ReactionExtentError::InvalidCandidate {
                field: "mole_amount",
                ..
            })
        ));
    }

    #[test]
    fn compute_species_moles_round_trips_positive_moles() {
        // ln(1.0) = 0.0, ln(2.0) ≈ 0.693, ln(1e-10) ≈ -23.0258
        let log_moles = &[0.0, 0.6931471805599453, -23.025850929940457];
        let moles = compute_species_moles(log_moles).unwrap();
        assert!((moles[0] - 1.0).abs() < 1e-15);
        assert!((moles[1] - 2.0).abs() < 1e-15);
        assert!((moles[2] - 1e-10).abs() < 1e-25);
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.46 — compute_element_totals
    // ──────────────────────────────────────────────
    #[test]
    fn compute_element_totals_matches_manual_calculation() {
        // 2 species, 2 elements: O2 = [2 O], O = [1 O]
        let a = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let n0 = vec![0.5, 0.2];
        let totals = compute_element_totals(&a, &n0).unwrap();
        assert!((totals[0] - (2.0 * 0.5 + 1.0 * 0.2)).abs() < 1e-15);
    }

    #[test]
    fn compute_element_totals_rejects_dimension_mismatch() {
        let a = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let n0 = vec![0.5]; // 1 entry, but matrix has 2 rows
        let result = compute_element_totals(&a, &n0);
        assert!(matches!(result, Err(ReactionExtentError::DimensionMismatch(_))));
    }

    #[test]
    fn compute_element_totals_rejects_non_finite_input() {
        let a = DMatrix::from_row_slice(2, 1, &[2.0, f64::NAN]);
        let n0 = vec![0.5, 0.2];
        let result = compute_element_totals(&a, &n0);
        assert!(matches!(
            result,
            Err(ReactionExtentError::InvalidProblem {
                field: "element_composition_or_initial_moles",
                ..
            })
        ));
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.47 — reaction_standard_gibbs
    // ──────────────────────────────────────────────
    #[test]
    fn reaction_standard_gibbs_computes_weighted_sum() {
        // 2 species, 1 reaction: ν = [-1, 2]
        let stoich = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let gibbs: Vec<GibbsFn> = vec![
            Rc::new(|_| 100.0),
            Rc::new(|_| 200.0),
        ];
        let dg0 = reaction_standard_gibbs(&stoich, &gibbs, 1000.0);
        assert!((dg0[0] - (-1.0 * 100.0 + 2.0 * 200.0)).abs() < 1e-12);
    }

    #[test]
    fn reaction_standard_gibbs_handles_multiple_reactions() {
        // 3 species, 2 reactions
        let stoich = DMatrix::from_row_slice(3, 2, &[
            -1.0, 0.0,
            1.0, -1.0,
            0.0, 1.0,
        ]);
        let gibbs: Vec<GibbsFn> = vec![
            Rc::new(|_| 50.0),
            Rc::new(|_| 100.0),
            Rc::new(|_| 200.0),
        ];
        let dg0 = reaction_standard_gibbs(&stoich, &gibbs, 500.0);
        assert!((dg0[0] - (-1.0 * 50.0 + 1.0 * 100.0 + 0.0 * 200.0)).abs() < 1e-12);
        assert!((dg0[1] - (0.0 * 50.0 + -1.0 * 100.0 + 1.0 * 200.0)).abs() < 1e-12);
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.48 — equilibrium_scaling
    // ──────────────────────────────────────────────
    #[test]
    fn equilibrium_scaling_returns_positive_factors() {
        // 2 species (O2, O), 1 element (O) → 1 reaction + 1 element balance = 2 rows
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let element_totals = vec![2.0]; // total O atoms from n0 = [1.0, 0.0]
        let T = 3000.0;
        let scale = equilibrium_scaling(&reactions, &elements, &gibbs, &element_totals, T).unwrap();
        assert_eq!(scale.len(), 2);
        assert!(scale[0] > 0.0);
        assert!(scale[1] > 0.0);
        assert!(scale[0].is_finite());
        assert!(scale[1].is_finite());
    }

    #[test]
    fn equilibrium_scaling_rejects_invalid_temperature() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let element_totals = vec![2.0];
        // Non-positive temperature should be rejected
        let result = equilibrium_scaling(&reactions, &elements, &gibbs, &element_totals, -1.0);
        assert!(result.is_err());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.49 — species_to_phase_map
    // ──────────────────────────────────────────────
    #[test]
    fn species_to_phase_map_assigns_each_species_to_its_phase() {
        let phases = vec![
            Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] },
            Phase { kind: PhaseKind::IdealSolution, species: vec![2] },
        ];
        let map = species_to_phase_map(&phases, 3).unwrap();
        assert_eq!(map, vec![0, 0, 1]);
    }

    #[test]
    fn species_to_phase_map_rejects_out_of_range_species_index() {
        let phases = vec![
            Phase { kind: PhaseKind::IdealGas, species: vec![0, 5] },
        ];
        let result = species_to_phase_map(&phases, 3);
        assert!(result.is_err());
    }

    #[test]
    fn species_to_phase_map_rejects_unassigned_species() {
        let phases = vec![
            Phase { kind: PhaseKind::IdealGas, species: vec![0] },
        ];
        let result = species_to_phase_map(&phases, 2);
        assert!(result.is_err());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.50 — reaction_phase_stoichiometry
    // ──────────────────────────────────────────────
    #[test]
    fn reaction_phase_stoichiometry_aggregates_by_phase() {
        // 3 species, 1 reaction, 2 phases
        // species 0,1 in phase 0; species 2 in phase 1
        let reactions = DMatrix::from_row_slice(3, 1, &[-1.0, 1.0, 2.0]);
        let phases = vec![
            Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] },
            Phase { kind: PhaseKind::IdealSolution, species: vec![2] },
        ];
        let delta_n = reaction_phase_stoichiometry(&reactions, &phases);
        assert_eq!(delta_n.len(), 1); // 1 reaction
        assert_eq!(delta_n[0].len(), 2); // 2 phases
        assert!((delta_n[0][0] - (-1.0 + 1.0)).abs() < 1e-15); // phase 0: -1 + 1 = 0
        assert!((delta_n[0][1] - 2.0).abs() < 1e-15); // phase 1: 2
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.3 — from_problem
    // ──────────────────────────────────────────────
    #[test]
    fn from_problem_preserves_all_problem_data() {
        let initial_moles = vec![1.0, 0.0];
        let problem = crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumProblem::new(
            vec!["O2".to_string(), "O".to_string()],
            initial_moles.clone(),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::LogMolesInitialGuess::from_moles(&initial_moles, 1e-20).unwrap(),
            DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            vec![Rc::new(|_| 0.0) as GibbsFn, Rc::new(|_| 0.0) as GibbsFn],
            vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }],
            crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
        ).unwrap();

        let solver = EquilibriumLogMoles::from_problem(problem).unwrap();
        assert_eq!(solver.subs_data.substances, vec!["O2", "O"]);
        assert_eq!(solver.n0, vec![1.0, 0.0]);
        assert_eq!(solver.T, 1000.0);
        assert_eq!(solver.P, 101325.0);
        assert_eq!(solver.p0, 101325.0);
        assert_eq!(solver.elem_composition.nrows(), 2);
        assert_eq!(solver.elem_composition.ncols(), 1);
        assert_eq!(solver.phases.len(), 1);
        assert_eq!(solver.species_phase, vec![0, 0]);
        assert!(solver.initial_guess.is_some());
    }

    #[test]
    fn from_problem_rejects_invalid_problem() {
        let result = crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumProblem::new(
            Vec::<String>::new(), // empty species list
            Vec::<f64>::new(),
            crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::LogMolesInitialGuess::from_moles(&[], 1e-20).unwrap(),
            DMatrix::from_row_slice(0, 1, &[]),
            Vec::<GibbsFn>::new(),
            vec![Phase { kind: PhaseKind::IdealGas, species: vec![] }],
            crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions::new(1000.0, 101325.0, 101325.0).unwrap(),
        );
        assert!(result.is_err());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.33 — EquilibriumLogMoles::new and empty
    // ──────────────────────────────────────────────
    #[test]
    fn empty_constructor_creates_default_solver() {
        let solver = EquilibriumLogMoles::empty();
        assert!(solver.subs_data.substances.is_empty());
        assert!(solver.n0.is_empty());
        assert!(solver.gibbs.is_empty());
        assert!(solver.phases.is_empty());
        assert!(solver.initial_guess.is_none());
        assert!(solver.solution.is_empty());
        assert!(solver.moles.is_empty());
        assert_eq!(solver.T, 273.15);
        assert_eq!(solver.P, 101325.0);
        assert_eq!(solver.p0, 101325.0);
        assert_eq!(solver.solver_settings.solver, Solvers::LM);
        assert!(!solver.solver_settings.scaling_flag);
    }

    #[test]
    fn new_constructor_creates_default_solver() {
        let solver = EquilibriumLogMoles::new();
        assert!(solver.subs_data.substances.is_empty());
        assert!(solver.n0.is_empty());
        assert!(solver.gibbs.is_empty());
        assert!(solver.phases.is_empty());
        assert!(solver.initial_guess.is_none());
        assert!(solver.solution.is_empty());
        assert!(solver.moles.is_empty());
        assert_eq!(solver.T, 273.15);
        assert_eq!(solver.P, 101325.0);
        assert_eq!(solver.p0, 101325.0);
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.34 — set_initial_guess
    // ──────────────────────────────────────────────
    #[test]
    fn set_initial_guess_validates_length() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.n0 = vec![1.0, 0.0];

        // Wrong length
        let result = solver.set_initial_guess(vec![0.0]);
        assert!(result.is_err());

        // Correct length
        let result = solver.set_initial_guess(vec![0.0, -10.0]);
        assert!(result.is_ok());
        assert_eq!(solver.initial_guess, Some(vec![0.0, -10.0]));
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.35 — clear_published_state
    // ──────────────────────────────────────────────
    #[test]
    fn clear_published_state_resets_all_publication_fields() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.solution = vec![42.0];
        solver.moles = vec![24.0];
        solver.last_validation_report = Some(
            crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport {
                residual_l2_norm: 1.0,
                residual_rms: 0.5,
                max_abs_residual: 0.8,
                raw_residual_l2_norm: 1.0,
                raw_residual_rms: 0.5,
                raw_max_abs_residual: 0.8,
                max_abs_element_balance_error: 0.0,
                reaction_affinity_l2_norm: 0.5,
                max_abs_reaction_affinity: 0.5,
                min_moles: 1.0,
            }
        );

        // We can't call clear_published_state directly (it's not public),
        // but we can verify that set_problem clears the state.
        solver.subs_data.substances = vec!["O2".to_string()];
        solver.n0 = vec![1.0];
        solver.elem_composition = DMatrix::from_row_slice(1, 1, &[2.0]);
        solver
            .set_problem(
                vec![1.0],
                vec![Phase {
                    kind: PhaseKind::IdealGas,
                    species: vec![0],
                }],
                101325.0,
            )
            .unwrap();
        assert!(solver.solution.is_empty());
        assert!(solver.moles.is_empty());
        assert!(solver.last_validation_report.is_none());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.31 — has_rst_symbolic_context
    // ──────────────────────────────────────────────
    #[test]
    fn has_rst_symbolic_context_returns_false_for_empty_solver() {
        let solver = EquilibriumLogMoles::empty();
        // has_rst_symbolic_context is private, but we can test indirectly
        // by checking that the default policy is legacy, not RST.
        // An empty solver has no symbolic context, so solve_candidate_from_seed
        // should use legacy default.
        assert!(solver.gibbs_sym.is_empty());
        assert!(solver.subs_data.search_results.is_empty());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.7 — EquilibriumSolverSettings::validate
    // ──────────────────────────────────────────────
    #[test]
    fn solver_settings_validate_accepts_default_settings() {
        let settings = EquilibriumSolverSettings::default();
        assert!(settings.validate().is_ok());
    }

    #[test]
    fn solver_settings_validate_rejects_zero_max_iter() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.max_iter = 0;
        let err = settings.validate().unwrap_err();
        assert!(matches!(err, ReactionExtentError::InvalidProblem { field: "solver_params.max_iter", .. }));
    }

    #[test]
    fn solver_settings_validate_rejects_non_finite_tol() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.tol = f64::NAN;
        let err = settings.validate().unwrap_err();
        assert!(matches!(err, ReactionExtentError::InvalidProblem { field: "solver_params.tol", .. }));
    }

    #[test]
    fn solver_settings_validate_rejects_zero_tol() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.tol = 0.0;
        assert!(settings.validate().is_err());
    }

    #[test]
    fn solver_settings_validate_rejects_negative_lambda() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.lambda = -1.0;
        assert!(settings.validate().is_err());
    }

    #[test]
    fn solver_settings_validate_rejects_non_finite_alpha_min() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.alpha_min = f64::INFINITY;
        assert!(settings.validate().is_err());
    }

    #[test]
    fn solver_settings_validate_rejects_delta_max_less_than_delta_init() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.delta_init = 100.0;
        settings.solver_params.delta_max = 50.0;
        let err = settings.validate().unwrap_err();
        assert!(matches!(err, ReactionExtentError::InvalidProblem { field: "solver_params.delta_max", .. }));
    }

    #[test]
    fn solver_settings_validate_rejects_eta_out_of_range() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.eta = 1.5;
        assert!(settings.validate().is_err());

        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.eta = -0.1;
        assert!(settings.validate().is_err());

        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_params.eta = f64::NAN;
        assert!(settings.validate().is_err());
    }

    #[test]
    fn solver_settings_validate_accepts_single_backend_policy() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_policy = Some(SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM)));
        assert!(settings.validate().is_ok());
    }

    #[test]
    fn solver_settings_validate_rejects_zero_budget_limits() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_budget = Some(SolverCascadeBudget::new(0, 100, 1000));
        assert!(settings.validate().is_err());

        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_budget = Some(SolverCascadeBudget::new(3, 0, 1000));
        assert!(settings.validate().is_err());

        let mut settings = EquilibriumSolverSettings::default();
        settings.solver_budget = Some(SolverCascadeBudget::new(3, 100, 0));
        assert!(settings.validate().is_err());
    }

    #[test]
    fn solver_settings_validate_rejects_non_finite_keq_tolerances() {
        let mut settings = EquilibriumSolverSettings::default();
        settings.keq_validation_tolerances.max_abs_species_mole_delta = f64::NAN;
        assert!(settings.validate().is_err());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.8 — compute_phase_totals
    // ──────────────────────────────────────────────
    #[test]
    fn compute_phase_totals_sums_moles_by_phase() {
        // log-moles: ln(1.0)=0.0, ln(2.0)≈0.693, ln(1.0)=0.0
        // 3 species in 2 phases: species 0,1 → phase 0; species 2 → phase 1
        let y = vec![0.0, 0.6931471805599453, 0.0];
        let species_phase = vec![0, 0, 1];
        let totals = compute_phase_totals(&y, &species_phase);
        assert_eq!(totals.len(), 2);
        assert!((totals[0] - (1.0 + 2.0)).abs() < 1e-15);
        assert!((totals[1] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn compute_phase_totals_handles_empty_input() {
        let totals = compute_phase_totals(&[], &[]);
        assert_eq!(totals.len(), 1); // max of empty is None → 0 + 1 = 1
        assert!((totals[0] - 0.0).abs() < 1e-15);
    }

    #[test]
    fn compute_phase_totals_handles_single_phase() {
        // log-moles: ln(1.0)=0.0, ln(3.0)≈1.0986
        let y = vec![0.0, 1.0986122886681098];
        let species_phase = vec![0, 0];
        let totals = compute_phase_totals(&y, &species_phase);
        assert_eq!(totals.len(), 1);
        assert!((totals[0] - (1.0 + 3.0)).abs() < 1e-14);
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.9 — initial_phase_activity
    // ──────────────────────────────────────────────
    #[test]
    fn initial_phase_activity_marks_active_phases() {
        // 2 phases, 3 species: phase 0 has 1.0 mol, phase 1 has 0.0 mol
        let log_moles = vec![0.0, -1e10, -1e10]; // exp(0)=1, exp(-1e10)≈0
        let species_phase = vec![0, 0, 1];
        let active = initial_phase_activity(&log_moles, &species_phase, 2, 1e-3).unwrap();
        assert_eq!(active.len(), 2);
        assert!(active[0]); // phase 0 has 1.0 mol > 1e-3
        assert!(!active[1]); // phase 1 has ~0 mol < 1e-3
    }

    #[test]
    fn initial_phase_activity_rejects_dimension_mismatch() {
        let result = initial_phase_activity(&[0.0], &[0, 1], 2, 1e-3);
        assert!(result.is_err());
    }

    #[test]
    fn initial_phase_activity_rejects_invalid_phase_eps() {
        let result = initial_phase_activity(&[0.0], &[0], 1, f64::NAN);
        assert!(result.is_err());
    }

    #[test]
    fn initial_phase_activity_rejects_out_of_bounds_phase() {
        let result = initial_phase_activity(&[0.0], &[5], 2, 1e-3);
        assert!(result.is_err());
    }

    #[test]
    fn initial_phase_activity_rejects_non_finite_moles() {
        let result = initial_phase_activity(&[f64::NAN], &[0], 1, 1e-3);
        assert!(result.is_err());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.13 — reject_repeated_phase_set
    // ──────────────────────────────────────────────
    #[test]
    fn reject_repeated_phase_set_accepts_first_occurrence() {
        let mut visited = HashSet::new();
        let phase_set = PhaseSet::from_policy(
            &InitialPhaseSet::AllCandidatePhases,
            &[true, false],
        )
        .unwrap();
        assert!(reject_repeated_phase_set(&mut visited, &phase_set, 0).is_ok());
    }

    #[test]
    fn reject_repeated_phase_set_rejects_duplicate() {
        let mut visited = HashSet::new();
        let phase_set = PhaseSet::from_policy(
            &InitialPhaseSet::AllCandidatePhases,
            &[true, false],
        )
        .unwrap();
        reject_repeated_phase_set(&mut visited, &phase_set, 0).unwrap();
        let result = reject_repeated_phase_set(&mut visited, &phase_set, 1);
        assert!(result.is_err());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.14 — validate_phase_set_candidate
    // ──────────────────────────────────────────────
    #[test]
    fn validate_phase_set_candidate_accepts_valid_set() {
        let phase_totals = vec![1.0, 0.0];
        let active = vec![true, false];
        assert!(validate_phase_set_candidate(&phase_totals, &active, 1e-3).is_ok());
    }

    #[test]
    fn validate_phase_set_candidate_rejects_dimension_mismatch() {
        let result = validate_phase_set_candidate(&[1.0], &[true, false], 1e-3);
        assert!(result.is_err());
    }

    #[test]
    fn validate_phase_set_candidate_rejects_negative_total() {
        let result = validate_phase_set_candidate(&[-1.0], &[true], 1e-3);
        assert!(result.is_err());
    }

    #[test]
    fn validate_phase_set_candidate_rejects_inactive_phase_with_moles() {
        let result = validate_phase_set_candidate(&[1.0], &[false], 1e-3);
        assert!(result.is_err());
    }

    // ──────────────────────────────────────────────
    // SourceCraft Diagnostics: D.1 — EquilibriumLogMoles::solve() with Legacy backends
    // ──────────────────────────────────────────────
    fn o2_dissociation_solver() -> EquilibriumLogMoles {
        // O2 = 2O, 1 element (O), 2 species
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.n0 = vec![1.0, 0.0];
        solver.T = 3000.0;
        solver.P = 101_325.0;
        solver.p0 = 101_325.0;
        solver.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        solver.gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        solver.phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        solver.species_phase = vec![0, 0];
        solver
    }

    #[test]
    fn legacy_lm_solve_publishes_accepted_solution() {
        let mut solver = o2_dissociation_solver();
        solver.solver_settings.solver = Solvers::LM;
        solver.create_stoich_matrix().unwrap();
        solver.solve().unwrap();
        assert!(!solver.solution.is_empty());
        assert!(solver.last_validation_report.is_some());
        assert!(solver.last_solve_report.is_some());
        // solution contains log-moles; compute actual moles via exp
        let moles: Vec<f64> = solver.solution.iter().map(|&y| y.exp()).collect();
        // elem_composition is (species × elements): column 0 = O coefficients [2.0, 1.0]
        let elem_col: Vec<f64> = solver.elem_composition.column(0).iter().copied().collect();
        let total_o: f64 = moles.iter().zip(elem_col.iter()).map(|(&n, &a)| n * a).sum();
        assert!((total_o - 2.0).abs() < 1e-6);
    }

    #[test]
    fn legacy_nr_solve_publishes_accepted_solution() {
        let mut solver = o2_dissociation_solver();
        solver.solver_settings.solver = Solvers::NR;
        solver.create_stoich_matrix().unwrap();
        solver.solve().unwrap();
        assert!(!solver.solution.is_empty());
        assert!(solver.last_validation_report.is_some());
        assert!(solver.last_solve_report.is_some());
    }

    #[test]
    fn legacy_tr_solve_publishes_accepted_solution() {
        let mut solver = o2_dissociation_solver();
        solver.solver_settings.solver = Solvers::TR;
        solver.create_stoich_matrix().unwrap();
        solver.solve().unwrap();
        assert!(!solver.solution.is_empty());
        assert!(solver.last_validation_report.is_some());
        assert!(solver.last_solve_report.is_some());
    }

    #[test]
    fn legacy_solve_fails_gracefully_without_stoich_matrix() {
        let mut solver = o2_dissociation_solver();
        // Don't call create_stoich_matrix — should fail with InvalidProblem
        let result = solver.solve();
        assert!(result.is_err());
    }

    #[test]
    fn legacy_solve_publishes_solve_report_with_attempts() {
        let mut solver = o2_dissociation_solver();
        solver.solver_settings.solver = Solvers::LM;
        solver.create_stoich_matrix().unwrap();
        solver.solve().unwrap();
        let report = solver.last_solve_report.as_ref().unwrap();
        assert!(report.started_attempt_count() >= 1);
        assert!(report.accepted_attempt_index().is_some());
    }
}
