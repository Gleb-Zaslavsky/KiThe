 #[cfg(test)]
mod tests {
    use crate::Thermodynamics::User_PhaseOrSolution::{
        PhaseOrSolution, SubstanceSystemFactory, SubstancesContainer, SubstancePhaseMapping, CustomSubstance
    };
    use crate::Thermodynamics::User_PhaseOrSolution2::OnePhase;
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use nalgebra::DMatrix;
    use approx::assert_relative_eq;
    use std::collections::HashMap;
 use crate::Thermodynamics::User_PhaseOrSolution::ThermodynamicsCalculatorTrait;
    #[test]
    fn test_substance_phase_mapping_basic() {
        let subs = vec!["A".to_string(), "B".to_string()];
        let phases = vec![Some("p1".to_string()), Some("p2".to_string())];
        let mapping = SubstancePhaseMapping::new(subs.clone(), phases.clone());
        assert_eq!(mapping.all_substances, subs);
        assert_eq!(mapping.substance_to_phase, phases);
        assert_eq!(mapping.get_phase_for_substance(0).unwrap(), &Some("p1".to_string()));
    }

    #[test]
    fn test_create_system_single_and_multi_phase() {
        // Single phase
        let subs = vec!["CO".to_string(), "CO2".to_string()];
        let container = SubstancesContainer::SinglePhase(subs.clone());
        let res = SubstanceSystemFactory::create_system(
            container,
            None,
            vec!["NASA_gas".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );
        assert!(res.is_ok());
        match res.unwrap() {
            crate::Thermodynamics::User_PhaseOrSolution::CustomSubstance::OnePhase(op) => {
                assert_eq!(op.subs_data.substances, subs);
            }
            _ => panic!("expected one phase"),
        }

        // Multi phase
        let mut phases = HashMap::new();
        phases.insert("gas".to_string(), vec!["H2".to_string(), "O2".to_string()]);
        phases.insert("liquid".to_string(), vec!["H2O".to_string()]);
        let container = SubstancesContainer::MultiPhase(phases);
        let res = SubstanceSystemFactory::create_system(
            container,
            None,
            vec!["NASA_gas".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );
        assert!(res.is_ok());
        match res.unwrap() {
            crate::Thermodynamics::User_PhaseOrSolution::CustomSubstance::PhaseOrSolution(p) => {
                assert!(p.subs_data.contains_key(&Some("gas".to_string())));
                assert!(p.subs_data.contains_key(&Some("liquid".to_string())));
            }
            _ => panic!("expected phase or solution"),
        }
    }

    #[test]
    fn test_create_full_map_of_mole_numbers_phase_or_solution() {
        // Build PhaseOrSolution with two phases
        let mut pos = PhaseOrSolution::new();
        let mut s1 = SubsData::new();
        s1.substances = vec!["A".to_string(), "B".to_string()];
        let mut s2 = SubsData::new();
        s2.substances = vec!["B".to_string(), "C".to_string()];
        pos.subs_data.insert(Some("p1".to_string()), s1);
        pos.subs_data.insert(Some("p2".to_string()), s2);

        // Provide non-zero mole numbers only for p1 with only A present
        let mut input: HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)> =
            HashMap::new();
        input.insert(
            Some("p1".to_string()),
            (
                Some(1.0),
                Some(HashMap::from([("A".to_string(), 0.3)])),
            ),
        );

        let (full_map, vec_map, summed) = pos
            .create_full_map_of_mole_numbers(input.clone())
            .unwrap();

        // Check missing keys inserted for p1 (B should be present with 0.0)
        let p1_inner = full_map.get(&Some("p1".to_string())).unwrap();
        assert!(p1_inner.1.as_ref().unwrap().contains_key("B"));
        assert_eq!(p1_inner.1.as_ref().unwrap().get("B").unwrap(), &0.0);

        // Check vec conversion length
        let p1_vec = vec_map.get(&Some("p1".to_string())).unwrap();
        assert!(p1_vec.1.as_ref().unwrap().len() == 2);

        // Check summed map: A=0.3, B=0.0, C absent -> B should be 0.0 present only if present in at least one phase
        assert!(summed.contains_key("A"));
        assert_eq!(summed.get("A").unwrap(), &0.3);
    }

    #[test]
    fn test_indexed_moles_variables_and_symbolic_vars() {
        let mut pos = PhaseOrSolution::new();
        let mut s1 = SubsData::new();
        s1.substances = vec!["X".to_string(), "Y".to_string()];
        let mut s2 = SubsData::new();
        s2.substances = vec!["Z".to_string()];
        pos.subs_data.insert(Some("p1".to_string()), s1);
        pos.subs_data.insert(Some("p2".to_string()), s2);

        let (_map_indexed, vec_of_n_vars, np_vec, map_of_var_each) =
            pos.indexed_moles_variables().unwrap();

        // Np vector length equals number of phases
        assert_eq!(np_vec.len(), 2);
        // total number of n vars equals 3
        assert_eq!(vec_of_n_vars.len(), 3);
        // Check map_of_var_each_substance contains entries
        let p1_map = map_of_var_each.get(&Some("p1".to_string()));
        assert!(p1_map.is_some());
        assert!(p1_map.unwrap().contains_key("X"));
        assert!(p1_map.unwrap().contains_key("Y"));
    }

    #[test]
    fn test_calculate_lagrange_equations_sym_and_fun() {
        let mut pos = PhaseOrSolution::new();
        let mut s = SubsData::new();
        s.substances = vec!["A".to_string(), "B".to_string()];
        pos.subs_data.insert(Some("mix".to_string()), s);

        // set simple dG_sym: A=1.0, B=2.0
        let mut inner = HashMap::new();
        inner.insert("A".to_string(), Expr::Const(1.0));
        inner.insert("B".to_string(), Expr::Const(2.0));
        let mut outer = HashMap::new();
        outer.insert(Some("mix".to_string()), inner);
        pos.dG_sym = outer;

        // A should be nrows = #substances x #elements (we construct original A as (#substances x #elements))
        // For 2 substances and 1 element, create A as 2x1
        let A = DMatrix::from_vec(2, 1, vec![1.0, 1.0]);

        let eqs = pos.calculate_Lagrange_equations_sym(A.clone(), 300.0).unwrap();
        assert_eq!(eqs.len(), 2);

        // Now test functional version
        // build G_fun with simple closures returning constants
        let mut gfun_phase = HashMap::new();
        gfun_phase.insert(
            "A".to_string(),
            Box::new(|_T: f64, _n: Option<Vec<f64>>, _Np: Option<f64>| 1.0f64) as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
        );
        gfun_phase.insert(
            "B".to_string(),
            Box::new(|_T: f64, _n: Option<Vec<f64>>, _Np: Option<f64>| 2.0f64) as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>,
        );
        let mut G_fun = HashMap::new();
        G_fun.insert(Some("mix".to_string()), gfun_phase);

        let f = pos.calculate_Lagrange_equations_fun(A.clone(), G_fun, 300.0).unwrap();

        // Call with arbitrary values: T, n, Np, Lambda (Lambda length equals number of elements = 1)
        let vals = f(300.0, None, None, vec![0.5]);
        assert_eq!(vals.len(), 2);
        // For our simple setup, eq_i = sum_by_elements + G_i/(R*Tm)
        // sum_by_elements = 0.5 (since Lambda=0.5 and element count 1)
        // so eq for A should be 0.5 + 1/(R*300)
        let expected_a = 0.5 + 1.0f64 / (crate::Thermodynamics::User_PhaseOrSolution::R * 300.0);
        assert_relative_eq!(vals[0], expected_a, epsilon = 1e-12);
    }

    #[test]
    fn test_set_T_and_P_in_G_sym() {
        let mut pos = PhaseOrSolution::new();
        let mut s = SubsData::new();
        s.substances = vec!["A".to_string()];
        pos.subs_data.insert(Some("p".to_string()), s);

        let mut inner = HashMap::new();
        // construct simple symbolic expression including T and P
        inner.insert("A".to_string(), Expr::Var("T".to_string()) + Expr::Var("P".to_string()));
        let mut outer = HashMap::new();
        outer.insert(Some("p".to_string()), inner);
        pos.dG_sym = outer;

        // set P and T
        pos.set_P_to_sym_in_G_sym(101325.0);
        pos.set_T_to_sym_in_G_sym(298.15);

        // now dG_sym should contain constants (simplified)
        let v = pos.dG_sym.get(&Some("p".to_string())).unwrap();
        let expr = v.get("A").unwrap();
        // Expect expression to be simplified to a constant: 298.15 + 101325.0
        let expected = Expr::Const(298.15 + 101325.0);
        assert_eq!(format!("{}", expr.simplify()), format!("{}", expected.simplify()));
    }

     #[test]
    fn test_substance_system_factory_multi_phase() {
        let mut phase_substances = HashMap::new();
        phase_substances.insert("gas".to_string(), vec!["N2".to_string(), "O2".to_string()]);
        phase_substances.insert("gas2".to_string(), vec!["H2O".to_string()]);

        let container = SubstancesContainer::MultiPhase(phase_substances);

        let result = SubstanceSystemFactory::create_system(
            container,
            None,
            vec!["NASA_gas".to_string()],
            vec!["NIST".to_string()],
            None,
            false,
        );

        assert!(result.is_ok());

        // Test that we got a PhaseOrSolution
        let mut pos = match result.unwrap() {
            CustomSubstance::PhaseOrSolution(pos) => pos,
            _ => panic!("Expected PhaseOrSolution"),
        };

        // Test extract_all_thermal_coeffs
        let thermal_result = pos.extract_all_thermal_coeffs(400.0);
        assert!(thermal_result.is_ok());

        // Test calculate_therm_map_of_properties
        let therm_props_result = pos.calculate_therm_map_of_properties(400.0);
        assert!(therm_props_result.is_ok());

        // Test calculate_therm_map_of_sym
        let therm_sym_result = pos.calculate_therm_map_of_sym();
        assert!(therm_sym_result.is_ok());

        // Test get_all_substances
        let all_subs = pos.get_all_substances();
        assert_eq!(all_subs.len(), 3);
        assert!(all_subs.contains(&"N2".to_string()));
        assert!(all_subs.contains(&"O2".to_string()));
        assert!(all_subs.contains(&"H2O".to_string()));

        // Test extract_SubstancesContainer
        let container_result = pos.extract_SubstancesContainer();
        assert!(container_result.is_ok());
        match container_result.unwrap() {
            SubstancesContainer::MultiPhase(phases) => {
                assert!(phases.contains_key("gas"));
                assert!(phases.contains_key("gas2"));
                assert_eq!(phases["gas"].len(), 2);
                assert_eq!(phases["gas2"].len(), 1);
            }
            _ => panic!("Expected MultiPhase container"),
        }

        // Test calculate_elem_composition_and_molar_mass
        let elem_comp_result = pos.calculate_elem_composition_and_molar_mass(None);
        assert!(elem_comp_result.is_ok());
        let (matrix, molar_masses, elements) = elem_comp_result.unwrap();
        assert_eq!(matrix.nrows(), 3); // 3 substances
        assert!(matrix.ncols() >= 2); // At least N, O, H elements
        assert!(molar_masses.contains_key("N2"));
        assert!(molar_masses.contains_key("O2"));
        assert!(molar_masses.contains_key("H2O"));
        assert!(elements.contains(&"N".to_string()));
        assert!(elements.contains(&"O".to_string()));
        assert!(elements.contains(&"H".to_string()));

        // Test indexed_moles_variables
        let indexed_vars_result = pos.indexed_moles_variables();
        assert!(indexed_vars_result.is_ok());
        let (map_indexed, vec_of_n_vars, np_vec, map_of_var_each) = indexed_vars_result.unwrap();
        assert_eq!(vec_of_n_vars.len(), 3); // 3 substances total
        assert_eq!(np_vec.len(), 2); // 2 phases
        assert!(map_indexed.contains_key(&Some("gas".to_string())));
        assert!(map_indexed.contains_key(&Some("gas2".to_string())));
        assert!(map_of_var_each.contains_key(&Some("gas".to_string())));
        assert!(map_of_var_each.contains_key(&Some("gas2".to_string())));

        // Test calculate_Gibbs_sym
        let gibbs_sym_result = pos.calculate_Gibbs_sym(400.0);
        assert!(gibbs_sym_result.is_ok());
        assert!(!pos.dG_sym.is_empty());

        // Test calculate_Gibbs_fun
        pos.calculate_Gibbs_fun(400.0, 101325.0);
        assert!(!pos.dG_fun.is_empty());

        // Test calculate_S_sym
        let entropy_sym_result = pos.calculate_S_sym(400.0);
        assert!(entropy_sym_result.is_ok());
        assert!(!pos.dS_sym.is_empty());

        // Test calculate_S_fun
        pos.calculate_S_fun(400.0, 101325.0);
        assert!(!pos.dS_fun.is_empty());

        // Test calcutate_Gibbs_free_energy with mole numbers
        let mut n = HashMap::new();
        let gas_moles = vec![0.4, 0.2];

        let gas2_moles = vec![0.4];
        n.insert(Some("gas".to_string()), (Some(0.6), Some(gas_moles)));
        n.insert(Some("gas2".to_string()), (Some(0.4), Some(gas2_moles)));

        let gibbs_energy_result = pos.calcutate_Gibbs_free_energy(400.0, 101325.0, n.clone());
        assert!(gibbs_energy_result.is_ok());
        let gibbs_map = gibbs_energy_result.unwrap();
        assert!(gibbs_map.contains_key(&Some("gas".to_string())));
        assert!(gibbs_map.contains_key(&Some("gas2".to_string())));

        // Test calculate_S with mole numbers
        let entropy_result = pos.calculate_S(400.0, 101325.0, n);
        assert!(entropy_result.is_ok());
        assert!(!pos.dS.is_empty());

        // Test configure_system_properties
        let molar_masses = HashMap::from([
            ("N2".to_string(), 28.014),
            ("O2".to_string(), 31.998),
            ("H2O".to_string(), 18.015),
        ]);
        let config_result = pos.configure_system_properties(101325.0, Some("Pa".to_string()), molar_masses, Some("g/mol".to_string()));
        assert!(config_result.is_ok());

        // Test if_not_found_go_NIST
        let nist_result = pos.if_not_found_go_NIST();
        assert!(nist_result.is_ok());

        // Test extract_coeffs_if_current_coeffs_not_valid
        let coeffs_result = pos.extract_coeffs_if_current_coeffs_not_valid(400.0);
        assert!(coeffs_result.is_ok());

        // Test create_full_map_of_mole_numbers
        let mut partial_n = HashMap::new();
        let partial_gas_moles = HashMap::from([
            ("N2".to_string(), 0.5),
            // O2 missing - should be filled with 0.0
        ]);
        partial_n.insert(Some("gas".to_string()), (Some(0.6), Some(partial_gas_moles)));
        partial_n.insert(Some("gas2".to_string()), (Some(0.4), Some(HashMap::from([("H2O".to_string(), 0.4)]))));

        let full_map_result = pos.create_full_map_of_mole_numbers(partial_n);
        assert!(full_map_result.is_ok());
        let (full_map, vec_map, summed) = full_map_result.unwrap();

        // Check that missing substances were filled
        let gas_data = full_map.get(&Some("gas".to_string())).unwrap();
        let gas_moles = gas_data.1.as_ref().unwrap();
        assert!(gas_moles.contains_key("O2"));
        assert_relative_eq!(*gas_moles.get("O2").unwrap(), 0.0, epsilon = 1e-10);

        // Check summed map
        assert!(summed.contains_key("N2"));
        assert!(summed.contains_key("O2"));
        assert!(summed.contains_key("H2O"));
        assert_relative_eq!(*summed.get("N2").unwrap(), 0.5, epsilon = 1e-10);
        assert_relative_eq!(*summed.get("O2").unwrap(), 0.0, epsilon = 1e-10);
        assert_relative_eq!(*summed.get("H2O").unwrap(), 0.4, epsilon = 1e-10);

        // Test set_P_to_sym_in_G_sym and set_T_to_sym_in_G_sym
        pos.set_P_to_sym_in_G_sym(200000.0);
        pos.set_T_to_sym_in_G_sym(350.0);
        // Verify that symbolic expressions have been updated
        assert!(!pos.dG_sym.is_empty());
        let gas_sym = pos.dG_sym.get(&Some("gas".to_string())).unwrap();
        let n2_expr = gas_sym.get("N2").unwrap();
        // The expression should now contain constants instead of variables
        assert!(!format!("{}", n2_expr).contains("T"));
        assert!(!format!("{}", n2_expr).contains("P"));

        // Test calculate_Lagrange_equations_sym
        // First need to set up symbolic variables and Gibbs expressions
        let _ = pos.indexed_moles_variables();
        let _ = pos.calculate_Gibbs_sym(400.0);
        let A = DMatrix::from_vec(3, 2, vec![2.0, 0.0, 0.0, 2.0, 1.0, 0.0]); // N2, O2, H2O with N, O elements
        let lagrange_sym_result = pos.calculate_Lagrange_equations_sym(A.clone(), 400.0);
        assert!(lagrange_sym_result.is_ok());
        let lagrange_eqs = lagrange_sym_result.unwrap();
        assert_eq!(lagrange_eqs.len(), 3); // 3 substances
          
        // Test calculate_Lagrange_equations_fun
          /* 
        let lagrange_fun_result = pos.calculate_Lagrange_equations_fun(A.clone(), pos.dG_fun, 400.0);
        assert!(lagrange_fun_result.is_ok());
        let lagrange_fun = lagrange_fun_result.unwrap();
*/
        // Test the function with sample values
     //   let test_lambda = vec![1.0, 2.0]; // Lambda values for N and O elements
      //  let test_vals = lagrange_fun(400.0, None, None, test_lambda);
       // assert_eq!(test_vals.len(), 3); // 3 substances
  
        // Test calculate_Lagrange_equations_fun2
        let lagrange_fun2_result = pos.calculate_Lagrange_equations_fun2(A, 400.0, 101325.0, 400.0);
        assert!(lagrange_fun2_result.is_ok());
        let lagrange_fun2 = lagrange_fun2_result.unwrap();

        // Test the function with sample values
        let test_vals2 = lagrange_fun2(400.0, None, None, vec![0.5, 1.5]);
        assert_eq!(test_vals2.len(), 3); // 3 substances
      
    }
}
