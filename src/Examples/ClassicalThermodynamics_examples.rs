use std::vec;

use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::ThermodynamicCalculations;
use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics2::Thermodynamics;
use crate::Thermodynamics::User_PhaseOrSolution::SubstanceSystemFactory;
use crate::Thermodynamics::User_PhaseOrSolution::{CustomSubstance, SubstancesContainer};
use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use approx::assert_relative_eq;
use std::collections::HashMap;

pub fn SubsData_examples(thermotask: usize) {
    match thermotask {
        0 => {
            // calculating Gibbs free energy withoun concentration correction RT*ln(P/P°) + RT*ln(w_i)
            let subs = vec!["CO".to_string(), "CO2".to_string()];
            // calling instance of strucutre SubsData created to search substances data in the databases, store
            // search results and calculate thermo properties
            let mut subdata = SubsData::new();
            // Set up library priorities
            subdata.set_multiple_library_priorities(
                vec!["NASA_gas".to_string()],
                LibraryPriority::Priority,
            );
            // Set up phases
            subdata
                .map_of_phases
                .insert(subs[0].clone(), Some(Phases::Gas));
            subdata
                .map_of_phases
                .insert(subs[1].clone(), Some(Phases::Gas));
            subdata.substances = subs.clone();
            // Perform the search
            subdata.search_substances();
            // Calling instance of structure Thermodynamics to calculate thermo dG
            let mut thermo = Thermodynamics::new();
            // Set up basic parameters
            thermo.set_T(400.0);
            thermo.set_P(101325.0, None);
            // Add test substances
            thermo.vec_of_substances = SubstancesContainer::SinglePhase(subs.clone());

            // savng  the search results in the structure Thermodynamics
            thermo.subdata = CustomSubstance::OnePhase(subdata);
            // Calculate Gibbs free energy
            thermo.calculate_Gibbs_free_energy(400.0, None, None);
            // Calculate symbolic Gibbs free energy
            thermo.calculate_Gibbs_sym(400.0, None, None);
            // getting the results of the calculation
            thermo.calculate_Gibbs_fun(400.0);
            let map_of_gibbs = thermo.dG.clone().get(&None).unwrap().clone();
            let map_of_gibbs_sym = thermo.dG_sym.clone().get(&None).unwrap().clone();
            let t = thermo.clone();
            let map_of_gibbs_fun = t.dG_fun.get(&None).unwrap();
            let vec_of_substances =
                if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                    vec.clone()
                } else {
                    panic!("vec_of_substances is not a SinglePhase");
                };
            for substance in &vec_of_substances {
                println!("substance: {:?}", substance);
                println!("map_of_gibbs: {:?}", map_of_gibbs[substance]);
                println!("map_of_gibbs_sym: {:?}", map_of_gibbs_sym[substance]);
            }
            for substance in &vec_of_substances {
                let dG_value = map_of_gibbs[substance];
                let dG_sym = map_of_gibbs_sym[substance].clone();
                let dG_from_sym = dG_sym.lambdify1D()(400.0);
                assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
                let dG_from_fun = map_of_gibbs_fun.get(substance).unwrap()(400.0, None, None);
                assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);
            }

            println!("map_of_gibbs: {:?} \n", map_of_gibbs);
            let _ = thermo.pretty_print();
        }
        1 => {
            // calculating Gibbs free energy withoun concentration correction RT*ln(P/P°) + RT*ln(w_i)
            let subs = vec!["CO".to_string(), "CO2".to_string()];
            // calling instance of strucutre SubsData created to search substances data in the databases, store
            // search results and calculate thermo properties
            let mut subdata = SubsData::new();
            // Set up library priorities
            subdata.set_multiple_library_priorities(
                vec!["NASA_gas".to_string()],
                LibraryPriority::Priority,
            );
            // Set up phases
            subdata
                .map_of_phases
                .insert(subs[0].clone(), Some(Phases::Gas));
            subdata
                .map_of_phases
                .insert(subs[1].clone(), Some(Phases::Gas));
            subdata.substances = subs.clone();
            // Perform the search
            subdata.search_substances();
            // Calling instance of structure Thermodynamics to calculate thermo dG
            let mut thermo = Thermodynamics::new();
            // Set up basic parameters
            thermo.set_T(400.0);
            thermo.set_P(101325.0, None);
            let concentration = Some(vec![0.5, 0.5]);
            // Add test substances
            thermo.vec_of_substances = SubstancesContainer::SinglePhase(subs.clone());

            // savng  the search results in the structure Thermodynamics
            thermo.subdata = CustomSubstance::OnePhase(subdata); //subdata;
            // Calculate Gibbs free energy
            thermo.calculate_Gibbs_free_energy(400.0, concentration.clone(), None);
            // Calculate symbolic Gibbs free energy
            thermo.calculate_Gibbs_sym(
                400.0,
                Some(vec![
                    Expr::Var("w1".to_string()),
                    Expr::Var("w2".to_string()),
                ]),
                None,
            );
            // getting the results of the calculation
            thermo.calculate_Gibbs_fun(400.0);
            let map_of_gibbs = thermo.dG.clone().get(&None).unwrap().clone();
            let map_of_gibbs_sym = thermo.dG_sym.clone().get(&None).unwrap().clone();
            let t = thermo.clone();
            let map_of_gibbs_fun = t.dG_fun.get(&None).unwrap();
            let vec_of_substances =
                if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                    vec.clone()
                } else {
                    panic!("vec_of_substances is not a SinglePhase");
                };
            for substance in &vec_of_substances {
                println!("substance: {:?}", substance);
                println!("map_of_gibbs: {:?}", map_of_gibbs[substance]);
                println!("map_of_gibbs_sym: {:?}", map_of_gibbs_sym[substance]);
            }
            for substance in &vec_of_substances {
                let dG_value = map_of_gibbs[substance];

                let dG_from_fun =
                    map_of_gibbs_fun.get(substance).unwrap()(400.0, concentration.clone(), None);
                assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);

                let dG_sym = map_of_gibbs_sym[substance].clone();
                let dG_sym = dG_sym.set_variable("P", 101325.0);

                let dG_from_sym = dG_sym.lambdify_owned(vec!["T", "w1", "w2"]);
                let dG_from_sym = dG_from_sym(vec![400.0, 0.5, 0.5]);
                assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
            }

            //println!("map_of_gibbs: {:?} \n", map_of_gibbs);
            // let _ = thermo.pretty_print();
        }
        2 => {
            // calculating Gibbs free energy withoun concentration correction RT*ln(P/P°) + RT*ln(w_i)
            let subs = vec!["CO".to_string(), "CO2".to_string()];
            // calling instance of strucutre SubsData created to search substances data in the databases, store
            // search results and calculate thermo properties
            let mut subdata = SubsData::new();
            // Set up library priorities
            subdata.set_multiple_library_priorities(
                vec!["NASA_gas".to_string()],
                LibraryPriority::Priority,
            );
            // Set up phases
            subdata
                .map_of_phases
                .insert(subs[0].clone(), Some(Phases::Gas));
            subdata
                .map_of_phases
                .insert(subs[1].clone(), Some(Phases::Gas));
            subdata.substances = subs.clone();
            // Perform the search
            subdata.search_substances();
            // Calling instance of structure Thermodynamics to calculate thermo dG
            let mut thermo = Thermodynamics::new();
            // Set up basic parameters
            thermo.set_T(400.0);
            thermo.set_P(101325.0, None);
            let concentration = Some(vec![0.5, 0.5]);
            // Add test substances
            thermo.vec_of_substances = SubstancesContainer::SinglePhase(subs.clone());

            // savng  the search results in the structure Thermodynamics
            thermo.subdata = CustomSubstance::OnePhase(subdata); //subdata;
            // Calculate Gibbs free energy
            thermo.calculate_Gibbs_free_energy(400.0, concentration.clone(), None);
            // Calculate symbolic Gibbs free energy
            thermo.calculate_Gibbs_sym(
                400.0,
                Some(vec![
                    Expr::Var("w1".to_string()),
                    Expr::Var("w2".to_string()),
                ]),
                None,
            );
            thermo.set_P_to_sym();
            // getting the results of the calculation
            thermo.calculate_Gibbs_fun(400.0);
            let map_of_gibbs = thermo.dG.clone().get(&None).unwrap().clone();
            let map_of_gibbs_sym = thermo.dG_sym.clone().get(&None).unwrap().clone();
            let t = thermo.clone();
            let map_of_gibbs_fun = t.dG_fun.get(&None).unwrap();
            let vec_of_substances =
                if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                    vec.clone()
                } else {
                    panic!("vec_of_substances is not a SinglePhase");
                };
            for substance in &vec_of_substances {
                println!("substance: {:?}", substance);
                println!("map_of_gibbs: {:?}", map_of_gibbs[substance]);
                println!("map_of_gibbs_sym: {:?}", map_of_gibbs_sym[substance]);
            }
            for substance in &vec_of_substances {
                let dG_value = map_of_gibbs[substance];

                let dG_from_fun =
                    map_of_gibbs_fun.get(substance).unwrap()(400.0, concentration.clone(), None);
                assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);

                let dG_sym = map_of_gibbs_sym[substance].clone();
                let dG_from_sym = dG_sym.lambdify_owned(vec!["T", "w1", "w2"]);
                let dG_from_sym = dG_from_sym(vec![400.0, 0.5, 0.5]);
                assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
            }

            //println!("map_of_gibbs: {:?} \n", map_of_gibbs);
            // let _ = thermo.pretty_print();
        }
        3 => {
            use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::ThermodynamicCalculations;
            use crate::Thermodynamics::User_PhaseOrSolution::{
                SubstanceSystemFactory, SubstancesContainer,
            };

            use crate::Thermodynamics::User_PhaseOrSolution::ThermodynamicsCalculatorTrait;
            use nalgebra::{DMatrix, DVector};
            let subs = vec!["H2".to_string(), "O2".to_string()];
            let T = 400.0;
            let P = 101325.0;
            let mut n = HashMap::new();
            let map_of_moles_num =
                HashMap::from([("H2".to_string(), 0.5), ("O2".to_string(), 0.5)]);
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
                None,
                library_priorities,
                permitted_libraries,
                explicit_search_insructions,
                search_in_NIST,
            )
            .unwrap();
            let mut nv = HashMap::new();
            nv.insert(None, (Some(1.0), (Some(vec![0.5, 0.5]))));
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
            td.map_of_concentration =
                HashMap::from([("H2".to_string(), 0.5), ("O2".to_string(), 0.5)]);
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
            //  let all_mole_nums = n_sym.clone().get(&None).unwrap().1.clone().unwrap();
            td.composition_equation_sym().unwrap();
            let composition_equations: &Vec<Box<dyn Fn(DVector<f64>) -> f64>> =
                &td.solver.elements_conditions;
            let composition_equation_sym = &td.solver.elements_conditions_sym;
            let n = DVector::from_vec(vec![0.5, 0.5]);
            for (i, ref comp_eq) in composition_equations.iter().enumerate() {
                let eq_for_element = comp_eq(n.clone());
                let eq_for_element_sym = composition_equation_sym[i].clone();
                let eq_for_element_sym_value =
                    eq_for_element_sym.lambdify_wrapped()(vec![0.5, 0.5]);
                assert_relative_eq!(eq_for_element, eq_for_element_sym_value, epsilon = 1e-6);
            }
            println!("composition_equations: {:?}", composition_equation_sym);
            // symbolic variables representing Lagrangian multipliers and equilibrium concentrations
            td.create_indexed_variables();
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
        4 => {
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
                None,
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
            td.initial_composition().unwrap();
            // symbolic variables representing Lagrangian multipliers and equilibrium concentrations
            td.create_indexed_variables();
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
        5 => {
            let gas_subs = vec![
                "CO".to_string(),
                "O".to_string(),
                "CO2".to_string(),
                "O2".to_string(),
            ];
            let map_of_subs = HashMap::from([
                ("gas".to_string(), gas_subs.clone()),
                ("solid".to_string(), vec!["C".to_string()]),
            ]);
            let T = 400.0;
            let P = 101325.0;

            let search_in_NIST = false;
            let explicit_search_insructions = None;
            let library_priorities = vec!["NASA_gas".to_string()];
            let permitted_libraries = vec!["NUIG".to_string()];
            let container = SubstancesContainer::MultiPhase(map_of_subs.clone());
            let mut customsubs = SubstanceSystemFactory::create_system(
                container,
                None,
                library_priorities,
                permitted_libraries,
                explicit_search_insructions,
                search_in_NIST,
            )
            .unwrap();
            println!("{:#?} \n", customsubs);

            match customsubs {
                CustomSubstance::PhaseOrSolution(ref phase) => {
                    let subs_data = &phase.subs_data;
                    assert!(subs_data.get(&Some("gas".to_string())).is_some());
                    assert!(subs_data.get(&Some("solid".to_string())).is_some());
                }
                _ => panic!(),
            }

            let mut n = HashMap::new();
            let map_of_gas = HashMap::from([("O2".to_string(), 0.5)]);
            let map_of_solid = HashMap::from([("C".to_string(), 0.5)]);
            n.insert(Some("gas".to_string()), (Some(1.0), Some(map_of_gas)));
            n.insert(Some("solid".to_string()), (Some(1.0), Some(map_of_solid)));
            // create thermodynamics instance
            let td = customsubs.create_thermodynamics(T, P, Some(n), None);
            assert!(td.is_ok());

            let mut td = td.unwrap();
            td.set_P_to_sym();
            //  td.initial_composition().unwrap();
            // symbolic variables representing Lagrangian multipliers and equilibrium concentrations
            //   td.create_indexed_variables();
            println!("\n \n Lambda: {:#?} \n", td.solver.Lambda);
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

            td.solver.solve(None, 1e-4, 100, None, None, None);
        }
        6 => {
            let subs = vec!["CO2".to_string(), "CO".to_string(), "O2".to_string()];
            let T = 273.15;
            let P = 101325.0;

            let search_in_NIST = false;
            let explicit_search_insructions = None;
            let library_priorities = vec!["NASA_gas".to_string()];
            let permitted_libraries = vec!["NUIG".to_string()];
            let container = SubstancesContainer::SinglePhase(subs.clone());
            let mut customsubs = SubstanceSystemFactory::create_system(
                container,
                None,
                library_priorities,
                permitted_libraries,
                explicit_search_insructions,
                search_in_NIST,
            )
            .unwrap();
            let mut n = HashMap::new();
            let map_of_concentration = HashMap::from([
                ("CO".to_string(), 0.4999),
                ("CO2".to_string(), 0.5),
                ("O2".to_string(), 0.0001),
            ]);
            n.insert(None, (Some(1.0), Some(map_of_concentration)));
            // create thermodynamics instance
            let td = customsubs.create_thermodynamics(T, P, Some(n), None);
            assert!(td.is_ok());
            let mut td = td.unwrap();
            td.set_P_to_sym();
            td.initial_composition().unwrap();
            // symbolic variables representing Lagrangian multipliers and equilibrium concentrations
            td.create_indexed_variables();
            println!("Lambda: {:#?} \n", td.solver.Lambda);
            println!("n: {:#?} \n", td.solver.n);
            println!("n_sym: {:#?} \n", td.solver.Np);
            // calculate element composition matrix
            td.calculate_elem_composition_and_molar_mass(None);
            td.find_composition_for_const_TP().unwrap();

            td.solver
                .solve(None, 5.0 * 1e-3, 200, Some(0.005), None, None);
            let result = td.solver.map_of_solutions;

            println!("result: {:#?}", result);
        }
        _ => {
            panic!("Invalid test case");
        }
    }
}
