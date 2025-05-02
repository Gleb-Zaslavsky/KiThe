#[cfg(test)]
mod tests {
    use crate::Kinetics::molmass::calculate_molar_mass;
    use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::{
        ThermodynamicCalculations, Thermodynamics,
    };
    use crate::Thermodynamics::User_PhaseOrSolution::{
        CustomSubstance, SubstanceSystemFactory, SubstancesContainer,
    };
    use crate::Thermodynamics::User_substances::{DataType, LibraryPriority, Phases, SubsData};
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use approx::assert_relative_eq;
    use core::panic;
    use nalgebra::{DMatrix, DVector};
    use std::collections::{HashMap, HashSet};
    use std::vec;

    #[test]
    fn test_thermodynamics_new() {
        let thermo = Thermodynamics::new();
        assert_eq!(thermo.T, 298.15);
        assert_eq!(thermo.P, 1e5);
        assert!(thermo.vec_of_subs.is_empty());
        assert!(thermo.unique_elements.is_empty());
        assert!(thermo.dG.is_empty());
        assert!(thermo.dG_sym.is_empty());
        assert!(thermo.dS.is_empty());
        assert!(thermo.dS_sym.is_empty());
    }

    #[test]
    fn test_set_t_and_p() {
        let mut thermo = Thermodynamics::new();
        thermo.set_T(400.0);
        thermo.set_P(101325.0, None);
        assert_eq!(thermo.T, 400.0);
        assert_eq!(thermo.P, 101325.0);
    }

    #[test]
    fn test_calculate_gibbs_free_energy_direct() {
        // Setup test data
        let mut subdata = SubsData::new();
        let subs = vec!["H2".to_string(), "O2".to_string()];

        // Set library priorities
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);

        // Set phases
        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();

        subdata.search_substances();

        // Setup Thermodynamics instance
        let mut thermo = Thermodynamics::new();

        thermo.set_T(400.0);
        thermo.set_P(101325.0, None);
        thermo.subs_container = SubstancesContainer::SinglePhase(subs.clone());
        thermo.vec_of_subs = subs.clone();
        thermo.subdata = CustomSubstance::OnePhase(subdata);
        thermo.create_indexed_variables();
        let mut n = HashMap::new();
        n.insert(None, (Some(1.0), Some(vec![0.5, 0.5])));
        // Calculate Gibbs free energy
        thermo.calculate_Gibbs_free_energy(400.0, n);
        thermo.calculate_Gibbs_fun(400.0);
        thermo.calculate_Gibbs_sym(400.0);
        thermo.set_P_to_sym();
        thermo.set_T_to_sym();

        // println!("symbolic vars: {:?}", vars);
        // Check results
        let g = thermo.dG.get(&None).unwrap();
        let g_sym = thermo.dG_sym.get(&None).unwrap();
        let g_fun = thermo.dG_fun.get(&None).unwrap();
        for gas in subs.clone() {
            assert!(g.contains_key(&gas));
            assert!(g_sym.contains_key(&gas));
            assert!(g_fun.contains_key(&gas));
            let g_sym = g_sym.get(&gas).unwrap().clone();
            println!("I. g_sym for gas: {},\n:  {}", &gas, &g_sym);
            let g_sym_fun = g_sym.lambdify_wrapped();
            let g_sym_value = g_sym_fun(vec![0.5, 1.0]);

            let g = g.get(&gas).unwrap().clone();
            let g_fun = g_fun.get(&gas).unwrap();
            let g_fun_value = g_fun(400.0, Some(vec![0.5, 0.5]), Some(1.0));
            assert_relative_eq!(g_fun_value, g, epsilon = 1e-6);
            assert_relative_eq!(g_sym_value, g_fun_value, epsilon = 1e-6);
        }
        // test with lambdify owned

        for (i, gas) in subs.iter().enumerate() {
            let g_sym = g_sym.get(gas).unwrap().clone();
            println!("II. g_sym for gas: {},\n:  {}", &gas, &g_sym);
            let g_sym_fun = g_sym.lambdify_owned(vec![&format!("N{}", i), "Np"]);
            let g_sym_value = g_sym_fun(vec![0.5, 1.0]);

            let g = g.get(gas).unwrap().clone();
            assert_relative_eq!(g_sym_value, g, epsilon = 1e-6);
        }
        // test with extraction symbolic vars
        for (_, gas) in subs.iter().enumerate() {
            let g_sym = g_sym.get(gas).unwrap().clone();
            let sym_vars = g_sym.all_arguments_are_variables();
            let vars: Vec<&str> = sym_vars.iter().map(|s| s.as_str()).collect::<Vec<&str>>();
            println!("vars {:?}", vars);
            println!("III. g_sym for gas: {},\n:  {}", &gas, &g_sym);
            let g_sym_fun = g_sym.lambdify_owned(vars);
            let g_sym_value = g_sym_fun(vec![0.5, 1.0]);

            let g = g.get(gas).unwrap().clone();
            assert_relative_eq!(g_sym_value, g, epsilon = 1e-6);
        }
        //
        let vars = thermo.symbolic_vars.clone().get(&None).unwrap().clone();
        let Np = vars.0.clone().unwrap().to_string().clone();
        let number_of_moles = thermo
            .symbolic_vars_for_every_subs
            .clone()
            .get(&None)
            .unwrap()
            .clone();
        for (_, gas) in subs.iter().enumerate() {
            let g_sym = g_sym.get(gas).unwrap().clone();
            let moles_num = number_of_moles
                .get(gas)
                .unwrap()
                .clone()
                .to_string()
                .clone();
            println!("IV. g_sym for gas: {},\n:  {}", &gas, &g_sym);
            let g_sym_fun = g_sym.lambdify_owned(vec![&moles_num, &Np]);
            let g_sym_value = g_sym_fun(vec![0.5, 1.0]);

            let g = g.get(gas).unwrap().clone();
            assert_relative_eq!(g_sym_value, g, epsilon = 1e-6);
        }
    }

    /// Test the calculation of symbolic and direct Gibbs free energy for a given
    /// mixture of substances at specified temperature and pressure.
    ///
    /// This test performs the following actions:
    /// 1. Sets up test data with substances, temperature, pressure, and mole fractions.
    /// 2. Calculates Gibbs free energy directly and symbolically.
    /// 3. Asserts that both calculations succeed.
    /// 4. Compares direct and symbolic Gibbs free energy values using lambdified expressions.
    /// 5. Calculates and verifies the element composition matrix.
    /// 6. Sets initial concentrations and verifies element composition.
    /// 7. Compares composition equations between direct and symbolic calculations.
    /// 8. Creates and verifies symbolic variables for Lagrangian multipliers and equilibrium concentrations.
    /// 9. Asserts the correctness of symbolic and functional expressions for equilibrium conditions.
    /// 10. Prints various intermediate results for verification.

    #[test]
    fn test_calculate_gibbs_sym() {
        // Setup test data
        use crate::Thermodynamics::User_PhaseOrSolution::ThermodynamicsCalculatorTrait;
        let subs = vec!["H2".to_string(), "O2".to_string()];
        let T = 400.0;
        let P = 101325.0;
        let mut nv = HashMap::new();
        nv.insert(None, (Some(1.0), (Some(vec![0.5, 0.5]))));

        let mut n = HashMap::new();
        let map_of_moles_num = HashMap::from([("H2".to_string(), 0.5), ("O2".to_string(), 0.5)]);
        n.insert(None, (Some(1.0), (Some(map_of_moles_num))));

        let mut n_sym = HashMap::new();
        n_sym.insert(
            None,
            (
                Some(Expr::Var("Np".to_string())),
                Some(vec![
                    Expr::Var("N0".to_string()),
                    Expr::Var("N1".to_string()),
                ]),
            ),
        );
        let search_in_NIST = false;
        let explicit_search_insructions = None;
        let library_priorities = vec!["NASA_gas".to_string()];
        let permitted_libraries = vec!["NASA_gas".to_string()];
        let container = SubstancesContainer::SinglePhase(subs.clone());
        let mut customsubs = SubstanceSystemFactory::create_system(
            container,
            library_priorities,
            permitted_libraries,
            explicit_search_insructions,
            search_in_NIST,
        )
        .unwrap();

        let dG = customsubs.calcutate_Gibbs_free_energy(T, P, nv.clone());

        let dG_sym = customsubs.calculate_Gibbs_sym(T, n_sym.clone());
        assert!(dG.is_ok());
        assert!(dG_sym.is_ok());
        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(T, P, Some(n.clone()), None);
        assert!(td.is_ok());
        let mut td = td.unwrap();
        td.set_P_to_sym();
        // compare results from direct calculation and from lambdified symbolic expression
        let g_sym = td.dG_sym.get(&None).unwrap();
        let g_sym = g_sym.get(&"H2".to_string()).unwrap().clone();
        let g_sym_fun = g_sym.lambdify_wrapped();
        let g_sym_value = g_sym_fun(vec![400.0, 0.5, 0.5]);
        let g_sym_value_from_sym = g_sym_fun(vec![400.0, 0.5, 0.5]);
        assert_relative_eq!(g_sym_value, g_sym_value_from_sym, epsilon = 1e-6);
        td.calculate_elem_composition_and_molar_mass(None);
        // set initial concentrations
        td.map_of_concentration = HashMap::from([("H2".to_string(), 0.5), ("O2".to_string(), 0.5)]);
        // calculate element composition matrix
        let elem_composition_matrix = td.elem_composition_matrix.clone().unwrap();
        // unique (non duplicated) vector of elements
        let vec_of_elements = &td.unique_elements;
        assert_eq!(vec_of_elements.len(), 2);
        assert_eq!(
            vec_of_elements.clone(),
            vec!["H".to_string(), "O".to_string()]
        );
        let matrix_data = vec![
            2.0, 0.0, // H2
            0.0, 2.0, // O2
        ];
        let matrix_expected = DMatrix::from_row_slice(2, 2, &matrix_data);
        assert_eq!(elem_composition_matrix, matrix_expected);
        td.initial_composition().unwrap();
        let initial_composition = &td.initial_vector_of_elements;
        //sun [2.0, 0.0] *[0.5, 0.5] = sum [1.0, 0.0] = 1.0
        println!("initial_composition: {:?}", initial_composition);
        assert_eq!(initial_composition, &vec![1.0, 1.0]);
        td.composition_equations().unwrap();
        //
        td.composition_equation_sym().unwrap();
        let composition_equations: &Vec<Box<dyn Fn(DVector<f64>) -> f64>> =
            &td.solver.elements_conditions;
        let composition_equation_sym = &td.solver.elements_conditions_sym;
        let n = DVector::from_vec(vec![0.5, 0.5]);
        for (i, ref comp_eq) in composition_equations.iter().enumerate() {
            let eq_for_element = comp_eq(n.clone());
            let eq_for_element_sym = composition_equation_sym[i].clone();
            let eq_for_element_sym_value = eq_for_element_sym.lambdify_wrapped()(vec![0.5, 0.5]);
            assert_relative_eq!(eq_for_element, eq_for_element_sym_value, epsilon = 1e-6);
        }
        println!("composition_equations: {:?}", composition_equation_sym);
        // symbolic variables representing Lagrangian multipliers and equilibrium concentrations
        td.create_indexed_variables();
        let all_mole_nums = n_sym.clone().get(&None).unwrap().1.clone().unwrap();
        assert_eq!(all_mole_nums, td.solver.n.clone());
        let Lamda_sym = td.solver.Lambda.clone();
        let n = td.solver.n.clone();
        println!("Lamda_sym: {:?}", Lamda_sym);
        println!("n: {:?}", n);

        assert_eq!(Lamda_sym[0], Expr::Var("Lambda0".to_string()));
        assert_eq!(Lamda_sym[1], Expr::Var("Lambda1".to_string()));
        assert_eq!(n[0], Expr::Var("N0".to_string()));
        assert_eq!(n[1], Expr::Var("N1".to_string()));
        //
        td.set_T_to_sym();
        td.create_nonlinear_system_sym().unwrap();

        let eq_mu = td.solver.eq_mu.clone();
        assert_eq!(eq_mu.len(), 2);
        println!("eq_mu: {:?}", eq_mu);

        let vars_extracted = td.symbolic_variables_extract().unwrap();

        println!("vars_extracted: {:?}", vars_extracted);
        td.pretty_print_substances_verbose().unwrap();
        td.pretty_print_Lagrange_equations().unwrap();

        td.create_nonlinear_system_fun().unwrap();

        // get symbolic variables
        let mut sym_vars = td.solver.Lambda.clone();
        sym_vars.extend(td.solver.n.clone());
        let sym_vars: Vec<String> = sym_vars.iter().map(|x| x.to_string()).collect();
        let mut sym_vars: Vec<&str> = sym_vars.iter().map(|x| x.as_str()).collect();
        sym_vars.extend(vec!["Np"]);
        println!("list of symbolic variables: {:?}", sym_vars);
        let eq_s = &td.solver.eq_mu.clone();
        let eq_mu_fun = &td.solver.eq_mu_fun;
        let eq_mu_fun_value = eq_mu_fun(T, Some(vec![0.5, 0.5]), None, vec![0.5, 0.5]);
        for (i, eq) in eq_s.clone().iter().enumerate() {
            let subs = &td.vec_of_subs[i];
            println!("eq_{}: {}, {:?}", subs, eq, sym_vars);
            let eq_sym_fun = eq.lambdify_wrapped();
            let eq_sym_fun_value = eq_sym_fun(vec![0.5, 0.5, 1.0]);

            assert_relative_eq!(eq_sym_fun_value, eq_mu_fun_value[i], epsilon = 1e-6);
        }
        for (i, eq) in eq_s.iter().enumerate() {
            let subs = &td.vec_of_subs[i];
            println!("eq_{}: {}, {:?}", subs, eq, &sym_vars);
            let eq_sym_fun = eq.clone().lambdify_owned(sym_vars.clone());
            let eq_sym_fun_value = eq_sym_fun(vec![0.5, 0.5, 0.5, 0.5, 1.0]);

            assert_relative_eq!(eq_sym_fun_value, eq_mu_fun_value[i], epsilon = 1e-6);
        }
        td.create_sum_of_mole_numbers_sym().unwrap();
        td.form_full_system_sym().unwrap();
        let full_system = &td.solver.full_system_sym;
        println!("full_system: {:?}", full_system);
        println!("eq_s: {:?}", eq_s);
        td.pretty_print_full_system();
    }

    #[test]
    fn test_calculate_elem_composition_and_molar_mass() {
        // Setup test data
        let mut subdata = SubsData::new();
        let subs = vec!["H2O".to_string(), "CO2".to_string()];

        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();

        // Calculate molar masses
        let (h2o_mass, _h2o_comp) = calculate_molar_mass("H2O".to_string(), None);
        let (co2_mass, _co2_comp) = calculate_molar_mass("CO2".to_string(), None);

        subdata
            .hasmap_of_molar_mass
            .insert("H2O".to_string(), h2o_mass);
        subdata
            .hasmap_of_molar_mass
            .insert("CO2".to_string(), co2_mass);

        // Setup Thermodynamics instance
        let mut thermo = Thermodynamics::new();
        thermo.subs_container = SubstancesContainer::SinglePhase(subs.clone());
        thermo.vec_of_subs = subs.clone();
        thermo.subdata = CustomSubstance::OnePhase(subdata);

        // Calculate element composition
        thermo.calculate_elem_composition_and_molar_mass(None);

        // Verify results
        assert!(thermo.elem_composition_matrix.is_some());
        assert!(!thermo.unique_elements.is_empty());

        // Check that the unique elements include H, O, and C
        let elements = &thermo.unique_elements;
        assert!(elements.contains(&"H".to_string()));
        assert!(elements.contains(&"O".to_string()));
        assert!(elements.contains(&"C".to_string()));

        // Check composition matrix dimensions
        let matrix = thermo.elem_composition_matrix.as_ref().unwrap();
        assert_eq!(matrix.nrows(), 2); // 2 substances
        assert_eq!(matrix.ncols(), elements.len()); // Number of unique elements
    }

    #[test]
    fn test_initial_composition() {
        // Setup test data with H2O and CO2
        let mut subdata = SubsData::new();
        let subs = vec!["H2O".to_string(), "CO2".to_string()];

        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();

        // Setup Thermodynamics instance
        let mut thermo = Thermodynamics::new();
        thermo.subs_container = SubstancesContainer::SinglePhase(subs.clone());
        thermo.vec_of_subs = subs.clone();
        thermo.subdata = CustomSubstance::OnePhase(subdata);

        // Set concentrations
        thermo.map_of_concentration.insert("H2O".to_string(), 0.7);
        thermo.map_of_concentration.insert("CO2".to_string(), 0.3);

        // Create element composition matrix manually for test
        // H2O has 2 H atoms and 1 O atom
        // CO2 has 1 C atom and 2 O atoms
        // Elements order: H, O, C
        let matrix_data = vec![
            2.0, 1.0, 0.0, // H2O
            0.0, 2.0, 1.0, // CO2
        ];
        let matrix = DMatrix::from_row_slice(2, 3, &matrix_data);
        thermo.elem_composition_matrix = Some(matrix);
        thermo.unique_elements = vec!["H".to_string(), "O".to_string(), "C".to_string()];

        // Calculate initial composition
        let result = thermo.initial_composition();
        assert!(result.is_ok());

        // Verify results
        assert_eq!(thermo.initial_vector_of_elements.len(), 3); // 3 elements

        // Expected values:
        // H: 2.0 * 0.7 + 0.0 * 0.3 = 1.4
        // O: 1.0 * 0.7 + 2.0 * 0.3 = 1.3
        // C: 0.0 * 0.7 + 1.0 * 0.3 = 0.3
        assert_relative_eq!(thermo.initial_vector_of_elements[0], 1.4, epsilon = 1e-10);
        assert_relative_eq!(thermo.initial_vector_of_elements[1], 1.3, epsilon = 1e-10);
        assert_relative_eq!(thermo.initial_vector_of_elements[2], 0.3, epsilon = 1e-10);
    }

    fn setup_test_thermodynamics() -> Thermodynamics {
        // Create a basic SubsData instance
        let mut subdata = SubsData::new();

        // Add some test substances
        subdata.substances = vec!["H2O".to_string(), "H2".to_string(), "O2".to_string()];

        // Create a CustomSubstance from the SubsData
        let custom_substance = CustomSubstance::OnePhase(subdata);

        // Create a Thermodynamics instance
        let mut thermodynamics = Thermodynamics::new();
        thermodynamics.subdata = custom_substance;
        thermodynamics.vec_of_subs = vec!["H2O".to_string(), "H2".to_string(), "O2".to_string()];

        // Create a simple element composition matrix
        // H2O: 2 H, 1 O
        // H2: 2 H, 0 O
        // O2: 0 H, 2 O
        let matrix_data = vec![
            2.0, 1.0, // H2O
            2.0, 0.0, // H2
            0.0, 2.0, // O2
        ];
        let elem_matrix = DMatrix::from_row_slice(3, 2, &matrix_data);
        thermodynamics.elem_composition_matrix = Some(elem_matrix);
        thermodynamics.unique_elements = vec!["H".to_string(), "O".to_string()];

        // Set initial concentrations
        let mut concentrations = HashMap::new();
        concentrations.insert("H2O".to_string(), 0.25);
        concentrations.insert("H2".to_string(), 0.5);
        concentrations.insert("O2".to_string(), 0.25);
        thermodynamics.map_of_concentration = concentrations;
        thermodynamics.create_indexed_variables();
        thermodynamics
    }

    #[test]
    fn test_initial_composition2() {
        let mut thermodynamics = setup_test_thermodynamics();

        // Test the initial_composition method
        let result = thermodynamics.initial_composition();
        assert!(result.is_ok());

        // Check that the initial vector of elements is calculated correctly
        // H: 2*0.25 (H2O) + 2*0.5 (H2) + 0*0.25 (O2) = 1.5
        // O: 1*0.25 (H2O) + 0*0.0 (H2) + 2*0.25 (O2) = 0.75
        println!(
            "Initial vector of elements: {:#?}",
            thermodynamics.initial_vector_of_elements
        );
        assert_eq!(thermodynamics.initial_vector_of_elements.len(), 2);
        assert!((thermodynamics.initial_vector_of_elements[0] - 1.5).abs() < 1e-10);
        assert!((thermodynamics.initial_vector_of_elements[1] - 0.75).abs() < 1e-10);
    }

    #[test]
    fn test_composition_equations() {
        let mut thermodynamics = setup_test_thermodynamics();

        // First set up the initial composition
        let _ = thermodynamics.initial_composition();

        // Test the composition_equations method
        let result = thermodynamics.composition_equations();
        assert!(result.is_ok());

        // Check that the elements_conditions vector has the correct length
        assert_eq!(thermodynamics.solver.elements_conditions.len(), 2);
    }

    #[test]
    fn test_composition_equation_sym() {
        let mut thermodynamics = setup_test_thermodynamics();

        // First set up the initial composition
        let _ = thermodynamics.initial_composition();

        // Create symbolic concentration variables
        let n1 = Expr::Var("N0".to_string());
        let n2 = Expr::Var("N1".to_string());
        let n3 = Expr::Var("N2".to_string());
        let vec_of_concentrations_sym_expect = vec![n1, n2, n3];
        assert_eq!(vec_of_concentrations_sym_expect, thermodynamics.solver.n);
        // Test the composition_equation_sym method
        let result = thermodynamics.composition_equation_sym();
        assert!(result.is_ok());

        // Check that the elements_conditions_sym vector has the correct length
        assert_eq!(thermodynamics.solver.elements_conditions_sym.len(), 2);
    }

    #[test]
    fn test_calculate_equilibrium() {
        // Setup test data

        let subs = vec![
            "H2".to_string(),
            "O2".to_string(),
            "H2O".to_string(),
            "O".to_string(),
            "H".to_string(),
        ];
        let T = 400.0;
        let P = 101325.0;

        let search_in_NIST = false;
        let explicit_search_insructions = None;
        let library_priorities = vec!["NASA_gas".to_string()];
        let permitted_libraries = vec!["NUIG".to_string()];
        let container = SubstancesContainer::SinglePhase(subs.clone());
        let mut customsubs = SubstanceSystemFactory::create_system(
            container,
            library_priorities,
            permitted_libraries,
            explicit_search_insructions,
            search_in_NIST,
        )
        .unwrap();
        let mut n = HashMap::new();
        let map_of_concentration = HashMap::from([
            ("H2".to_string(), 0.5),
            ("O2".to_string(), 0.5),
            ("H2O".to_string(), 0.1),
        ]);
        n.insert(None, (Some(1.0), Some(map_of_concentration)));
        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(T, P, Some(n), None);
        assert!(td.is_ok());
        let mut td = td.unwrap();
        td.set_P_to_sym();
        //  td.initial_composition().unwrap();
        // symbolic variables representing Lagrangian multipliers and equilibrium concentrations
        //   td.create_indexed_variables();
        println!("Lambda: {:#?} \n", td.solver.Lambda);
        println!("n: {:#?} \n", td.solver.n);
        println!("n_sym: {:#?} \n", td.solver.Np);
        // calculate element composition matrix
        td.calculate_elem_composition_and_molar_mass(None);
        // set initial concentrations
        println!(
            "Initial vector of elements: {:#?} \n",
            td.initial_vector_of_elements
        );
        println!(
            "composition: {}, ncols = {}\n",
            &td.clone().elem_composition_matrix.unwrap(),
            &td.clone().elem_composition_matrix.unwrap().ncols()
        );
        td.composition_equations().unwrap();
        //
        td.composition_equation_sym().unwrap();
        //
        td.set_T_to_sym();
        td.create_nonlinear_system_sym().unwrap();
        td.create_nonlinear_system_fun().unwrap();
        td.create_sum_of_mole_numbers_sym().unwrap();
        td.form_full_system_sym().unwrap();
        td.pretty_print_full_system();
    }

    #[test]
    fn test_calculate_equilibrium2() {
        // Setup test data

        let subs = vec![
            "H2".to_string(),
            "O2".to_string(),
            "H2O".to_string(),
            "O".to_string(),
            "H".to_string(),
        ];
        let T = 400.0;
        let P = 101325.0;

        let search_in_NIST = false;
        let explicit_search_insructions = None;
        let library_priorities = vec!["NASA_gas".to_string()];
        let permitted_libraries = vec!["NUIG".to_string()];
        let container = SubstancesContainer::SinglePhase(subs.clone());
        let mut customsubs = SubstanceSystemFactory::create_system(
            container,
            library_priorities,
            permitted_libraries,
            explicit_search_insructions,
            search_in_NIST,
        )
        .unwrap();
        let mut n = HashMap::new();
        let map_of_concentration = HashMap::from([
            ("H2".to_string(), 0.5),
            ("O2".to_string(), 0.5),
            ("H2O".to_string(), 0.1),
        ]);
        n.insert(None, (Some(1.0), Some(map_of_concentration)));
        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(T, P, Some(n), None);
        assert!(td.is_ok());
        let mut td = td.unwrap();

        td.find_composition_for_const_TP().unwrap();
    }

    #[test]
    fn test_calculate_equilibrium3() {
        // Setup test data

        let subs = vec![
            "H2".to_string(),
            "O2".to_string(),
            "H2O".to_string(),
            "O".to_string(),
            "H".to_string(),
        ];
        let T = 400.0;
        let P = 101325.0;

        let search_in_NIST = false;
        let explicit_search_insructions = None;
        let library_priorities = vec!["NASA_gas".to_string()];
        let permitted_libraries = vec!["NUIG".to_string()];
        let container = SubstancesContainer::SinglePhase(subs.clone());
        let mut customsubs = SubstanceSystemFactory::create_system(
            container,
            library_priorities,
            permitted_libraries,
            explicit_search_insructions,
            search_in_NIST,
        )
        .unwrap();
        let mut n = HashMap::new();
        let map_of_concentration = HashMap::from([
            ("H2".to_string(), 0.5),
            ("O2".to_string(), 0.5),
            ("H2O".to_string(), 0.1),
        ]);
        n.insert(None, (Some(1.0), Some(map_of_concentration)));
        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(T, P, None, None);
        assert!(td.is_ok());
        let mut td = td.unwrap();
        td.set_number_of_moles(Some(n));
        td.initial_composition().unwrap();
        td.create_indexed_variables();
        td.find_composition_for_const_TP().unwrap();
    }
    /*
    #[test]
    fn test_indexed_variables() {
        let mut thermodynamics = setup_test_thermodynamics();

        // First set up the initial composition
        let _ = thermodynamics.initial_composition();

        // Set b0 in the solver
        thermodynamics.solver.b0 = vec![3.0, 1.5];

        // Test the indexed_variables method
        let result = thermodynamics.indexed_variables();
        assert!(result.is_ok());

        // Check that Lambda and n vectors have been populated
        assert_eq!(thermodynamics.solver.Lambda.len(), 2);
        assert_eq!(thermodynamics.solver.n.len(), 3);
    }

    #[test]
    fn test_create_nonlinear_system() {
        let mut thermodynamics = setup_test_thermodynamics();

        // Set up the necessary prerequisites
        let _ = thermodynamics.initial_composition();
        let _ = thermodynamics.indexed_variables();

        // Create dummy Gibbs energy expressions
        let mut dG_sym = HashMap::new();
        let mut inner_map = HashMap::new();
        inner_map.insert("H2O".to_string(), Expr::Var("G_H2O".to_string()));
        inner_map.insert("H2".to_string(), Expr::Var("G_H2".to_string()));
        inner_map.insert("O2".to_string(), Expr::Var("G_O2".to_string()));
        dG_sym.insert(None, inner_map);
        thermodynamics.dG_sym = dG_sym;

        // Test the create_nonlinear_system method
        let result = thermodynamics.create_nonlinear_system_sym();

        // This might fail if the subdata doesn't have the necessary methods implemented for testing
        // In a real test environment, you would need to mock these dependencies
        if result.is_ok() {
            assert!(!thermodynamics.solver.eq_mu.is_empty());
        }
    }

    #[test]
    fn test_pretty_print_substances_verbose() {
        let thermodynamics = setup_test_thermodynamics();

        // Test the pretty_print_substances_verbose method
        // This is primarily a visual output test, so we just check it doesn't error
        let result = thermodynamics.pretty_print_substances_verbose();
        assert!(result.is_ok());
    }

    #[test]
    fn test_thermodynamic_calculations_integration() {
        // Create a basic SubsData instance
        let mut subdata = SubsData::new();
        subdata.substances = vec!["H2O".to_string(), "H2".to_string(), "O2".to_string()];

        // Create a CustomSubstance from the SubsData
        let mut custom_substance = CustomSubstance::OnePhase(subdata);

        // Test the create_thermodynamics method
        let temperature = 298.15;
        let pressure = 101325.0;
        let concentrations = Some(vec![0.25, 0.5, 0.25]);

        // This will likely fail in a test environment without proper setup
        // In a real test, you would need to mock the necessary dependencies
        let result = custom_substance.create_thermodynamics(
            temperature,
            pressure,
            concentrations,
            None,
            None,
        );

        // If the test environment is properly set up, this should pass
        if result.is_ok() {
            let thermodynamics = result.unwrap();
            assert_eq!(thermodynamics.T, temperature);
            assert_eq!(thermodynamics.P, pressure);
            assert_eq!(thermodynamics.vec_of_subs.len(), 3);
        }
    }
    */
}
