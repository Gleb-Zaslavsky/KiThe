#[cfg(test)]
mod tests {
    use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq::{
        EquilibriumLogMoles, GibbsFn, Phase, PhaseKind, Solvers, compute_element_totals,
        equilibrium_logmole_jacobian, equilibrium_logmole_residual, equilibrium_scaling,
        reaction_phase_stoichiometry, species_to_phase_map,
    };
    use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq2::{
        finite_difference_jacobian, gas_solver, gas_solver_for_T_range,
        gas_solver_for_T_range_for_elements, gas_solver_from_elements,
    };
    use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq3::compute_reaction_basis;
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use approx::assert_relative_eq;

    use nalgebra::DMatrix;
    use std::collections::HashMap;
    use std::f64;
    use std::rc::Rc;
    use std::time::Instant;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
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
        let b = compute_element_totals(&a, &n0);

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
        let b = compute_element_totals(&a, &n0);

        assert!((b[0] - (2.0 * 1.0 + 1.0 * 0.5)).abs() < 1e-12); // N
        assert!((b[1] - (2.0 * 2.0 + 1.0 * 0.5)).abs() < 1e-12); // O
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
        let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();

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
        let elem_total = compute_element_totals(&elem, &n0);
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
        let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();

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
        let elem_total = compute_element_totals(&elem, &n0);
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
        let f = equilibrium_logmole_residual(
            reactions.clone(),
            elements.clone(),
            element_totals.clone(),
            gibbs.clone(),
            phases.clone(),
            3000.0,
            101325.0,
            101325.0,
            species_phase,
            0.0,
            0.0,
        )
        .unwrap();

        let species_phase = species_to_phase_map(&phases, 2).unwrap();
        let delta_n = reaction_phase_stoichiometry(&reactions, &phases);

        let x = vec![0.0, -10.0];

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
                    "Jacobian mismatch at ({},{}) analytic={} fd={}",
                    i,
                    k,
                    j_an[(i, k)],
                    j_fd[(i, k)]
                );
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

        let scale = equilibrium_scaling(&stoich, &elements, &gibbs, &element_totals, T);

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

        let scale = equilibrium_scaling(&stoich, &elements, &gibbs, &element_totals, 1000.0);

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

        let scale = equilibrium_scaling(&stoich, &elements, &gibbs, &element_totals, 4000.0);

        for s in scale {
            assert!(s.is_finite() && s > 0.0);
        }
    }

    ////////////////////////////////SOLUTION AT CONST T//////////////////////////////////////////////////////

    #[test]
    fn O2_O_equilibrium_logmoles_LM() {
        let T = 500.0;

        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false);
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
    fn O2_O_equilibrium_logmoles_NR() {
        let T = 500.0;

        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::NR, Some("info"), false);
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
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::TR, Some("info"), false);
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
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false);
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
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false);
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
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false);
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
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false);
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
        instance.scaling_flag = false;
        instance.set_initial_guess(vec![0.0, -10.0]);
        instance.set_problem(
            vec![1.0, 0.0],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            }],
            101325.0,
        );
        instance.subs_data = user_subs;
        instance.create_equilibrium_system().unwrap();
        instance
            .solve_for_T_range_par2(400.0, 4500.0, 10.0)
            .unwrap();
        instance.create_moles_table();
        println!("elapsed time {} ms", now.elapsed().as_millis());
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
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
        instance.scaling_flag = false;
        instance.set_initial_guess(vec![0.0, -10.0]);
        instance.set_problem(
            vec![1.0, 0.0],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            }],
            101325.0,
        );
        instance.subs_data = user_subs;
        instance.create_equilibrium_system().unwrap();
        instance.solve_for_T_range(400.0, 4500.0, 10.0).unwrap();
        instance.create_moles_table();
        println!("elapsed time {} ms", now.elapsed().as_millis());
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
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
        instance.scaling_flag = false;
        instance.set_initial_guess(vec![0.0, -10.0]);
        instance.with_loglevel(Some("info"));
        instance.set_problem(
            vec![1.0, 0.0],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            }],
            101325.0,
        );
        instance.subs_data = user_subs;
        instance.create_equilibrium_system().unwrap();
        instance.solve_for_T_range(400.0, 5500.0, 100.0).unwrap();
        instance.create_moles_table();

        println!("elapsed time {} ms", now.elapsed().as_millis());
        println!(
            "list of failed T points (if any): {:?}",
            instance.list_of_failed_T
        );
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
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
        instance.scaling_flag = false;
        instance.with_loglevel(Some("info"));
        //   instance.set_initial_guess(vec![1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]);
        instance.set_problem(
            vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1, 2, 3, 4, 5, 6, 7],
            }],
            101325.0,
        );
        instance.subs_data = user_subs;
        instance.create_equilibrium_system().unwrap();
        instance.solve_for_T_range(400.0, 5500.0, 100.0).unwrap();
        instance.create_moles_table();

        println!("elapsed time {} ms", now.elapsed().as_millis());
        println!(
            "list of failed T points (if any): {:?}",
            instance.list_of_failed_T
        );
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
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
        instance.scaling_flag = false;
        instance.with_loglevel(Some("info"));
        //   instance.set_initial_guess(vec![1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]);
        instance.set_problem(
            vec![1.0, 1e-5, 1.0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1, 2, 3, 4, 5, 6, 7],
            }],
            101325.0,
        );
        instance.subs_data = user_subs;
        instance.solver = Solvers::NR;
        instance.create_equilibrium_system().unwrap();
        instance.solve_for_T_range(400.0, 5500.0, 100.0).unwrap();
        instance.create_moles_table();

        println!("elapsed time {} ms", now.elapsed().as_millis());
        println!(
            "list of failed T points (if any): {:?}",
            instance.list_of_failed_T
        );
        assert!(!instance.moles_for_T_range.is_empty());
    }

    #[test]
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
    fn N2_N_O2_O_equilibrium_sym_low_T() {
        use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq2::multiphase_equilibrium_residual_generator_sym;
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
        let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();

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
        let elem_total = compute_element_totals(&elem, &n0);

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
        let map_of_functions = user_subs.calculate_dG0_sym_one_phase();
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
    fn N2_N_O2_O_equilibrium_sym_high_T() {
        use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq2::multiphase_equilibrium_residual_generator_sym;
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
        let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();

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
        let elem_total = compute_element_totals(&elem, &n0);

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
        let map_of_functions = user_subs.calculate_dG0_sym_one_phase();
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
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false);
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
        );

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
}
