///////////////////////TESTS////////////////////////////////////////////
#[cfg(test)]
mod tests {

    use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
    use crate::Thermodynamics::DBhandlers::transport_api::TransportCalculator;
    use crate::Thermodynamics::User_substances::{
        CalculatorType, DataType, LibraryPriority, Phases, SearchResult, SubsData, WhatIsFound,
    };
    use std::collections::HashMap;

    #[test]
    fn test_user_substances_with_calculators() {
        // Create a new SubsData instance with test substances
        let substances = vec!["CO".to_string(), "H2O".to_string()];
        let mut user_subs = SubsData::new();
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        user_subs.substances = substances.clone();

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Perform the search
        user_subs.search_substances();
        user_subs
            .extract_thermal_coeffs(substances[0].as_str(), 400.0)
            .unwrap();
        user_subs
            .extract_thermal_coeffs(substances[1].as_str(), 400.0)
            .unwrap();

        // Test thermo calculations
        if let Ok((cp, dh, ds)) = user_subs.calculate_thermo_properties("CO", 400.0) {
            println!("CO Thermo properties at 400K:");
            println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
            assert!(cp > 0.0);
        } else {
            panic!("Failed to calculate CO properties");
        }
        let Cp = 33.8;
        user_subs
            .extract_transport_coeffs(substances[0].as_str(), 400.0)
            .unwrap();
        user_subs
            .extract_transport_coeffs(substances[1].as_str(), 400.0)
            .unwrap();
        // Test transport calculations
        user_subs.set_M(
            HashMap::from([("H2O".to_string(), 18.0), ("CO".to_string(), 32.0)]),
            None,
        );
        user_subs.set_P(1e5, None);
        if let Ok((lambda, viscosity)) =
            user_subs.calculate_transport_properties("H2O", 400.0, Some(Cp), None)
        {
            println!("CO Transport properties at 400K:");
            println!("Lambda: {}, Viscosity: {}", lambda, viscosity);
            assert!(lambda > 0.0);
            assert!(viscosity > 0.0);
        } else {
            panic!("Failed to calculate H2O properties");
        }

        // Print full summary
        user_subs.print_search_summary();
    }
    #[test]
    fn test_user_substances_with_calculators_own_functions() {
        // Create a new SubsData instance with test substances
        let substances = vec!["CO".to_string(), "H2O".to_string()];
        let mut user_subs = SubsData::new();
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        user_subs.substances = substances.clone();

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Perform the search
        user_subs.search_substances();
        user_subs
            .extract_thermal_coeffs(substances[0].as_str(), 400.0)
            .unwrap();
        user_subs
            .extract_thermal_coeffs(substances[1].as_str(), 400.0)
            .unwrap();
        // Print full summary
        user_subs.print_search_summary();
        let datamap = user_subs.get_substance_result("CO").unwrap();
        let Thermo = datamap.get(&WhatIsFound::Thermo).unwrap().as_ref().unwrap();
        let Calculator = Thermo.calculator.as_ref().unwrap();
        let Cp;
        match Calculator {
            CalculatorType::Thermo(thermo) => {
                // Test thermo calculations
                let mut thermo = thermo.clone();
                thermo.extract_model_coefficients(400.0).unwrap();
                if let Ok(()) = thermo.calculate_Cp_dH_dS(400.0) {
                    let (cp, dh, ds) = (
                        thermo.get_Cp().unwrap(),
                        thermo.get_dh().unwrap(),
                        thermo.get_ds().unwrap(),
                    ); //thermo.get_Cp().unwrap();
                    println!("CO Thermo properties at 400K:");
                    println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
                    Cp = cp;
                    assert!(cp > 0.0);
                } else {
                    panic!("Failed to calculate CO properties");
                }
            }
            _ => {
                panic!("Failed to calculate CO properties");
            }
        }
        let Transport = datamap
            .get(&WhatIsFound::Transport)
            .unwrap()
            .as_ref()
            .unwrap();
        let Calculator = Transport.calculator.as_ref().unwrap();
        match Calculator {
            CalculatorType::Transport(transport) => {
                // Test transport calculations
                let mut transport = transport.clone();
                let _ = transport.set_M(32.0, None);
                let _ = transport.set_P(1e5, None);
                transport.extract_coefficients(400.0).unwrap();
                if let Ok(L) = transport.calculate_lambda(Some(Cp), None, 400.0) {
                    println!("CO Transport properties at 400K:");
                    println!("Lambda: {}", L);
                    assert!(L > 0.0);
                } else {
                    panic!("Failed to calculate H2O properties");
                }
            }
            _ => {
                panic!("Failed to calculate CO properties");
            }
        }
    }
    #[test]
    fn test_library_priority_setting() {
        let mut user_subs = SubsData::new();

        // Test setting single library priority
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        assert_eq!(
            user_subs.library_priorities.get("NASA_gas"),
            Some(&LibraryPriority::Priority)
        );

        // Test setting multiple library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NIST".to_string(), "CEA".to_string()],
            LibraryPriority::Permitted,
        );
        assert_eq!(
            user_subs.library_priorities.get("NIST"),
            Some(&LibraryPriority::Permitted)
        );
        assert_eq!(
            user_subs.library_priorities.get("CEA"),
            Some(&LibraryPriority::Permitted)
        );

        // Test overriding priority
        user_subs.set_library_priority("NIST".to_string(), LibraryPriority::Priority);
        assert_eq!(
            user_subs.library_priorities.get("NIST"),
            Some(&LibraryPriority::Priority)
        );
    }

    #[test]
    fn test_substance_search() {
        // fails
        let mut user_subs = SubsData::new();
        user_subs.substances = vec![
            "H2O".to_string(),
            "CO2".to_string(),
            "NonExistentSubstance".to_string(),
        ];

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Permitted,
        );

        // Perform search
        user_subs.search_substances();

        // Check results
        assert!(matches!(
            user_subs
                .get_substance_result("H2O")
                .unwrap()
                .get(&WhatIsFound::Thermo)
                .unwrap(),
            Some(SearchResult { .. })
        ));
        let _data_to_expect: HashMap<WhatIsFound, Option<SearchResult>> =
            HashMap::from([(WhatIsFound::NotFound, None)]);
        assert!(matches!(
            user_subs.get_substance_result("NonExistentSubstance"),
            Some(_data_to_expect)
        ));

        // Check substance lists
        let not_found = user_subs.get_not_found_substances();
        assert!(not_found.contains(&"NonExistentSubstance".to_string()));

        let priority_found = user_subs.get_priority_found_substances();
        let permitted_found = user_subs.get_permitted_found_substances();
        println!(
            "Priority found: {:?}, Permitted found: {:?}",
            priority_found, permitted_found
        );
        // At least one of these should contain H2O
        assert!(
            priority_found.contains(&"H2O".to_string())
                || permitted_found.contains(&"H2O".to_string())
        );
    }

    #[test]
    fn test_thermo_property_calculation() {
        let mut user_subs = SubsData::new();
        let substances = vec!["O2".to_string(), "CO".to_string()];
        user_subs.substances = substances.clone();

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("O2".to_string(), Some(Phases::Gas));

        // Perform search
        user_subs.search_substances();
        user_subs
            .extract_thermal_coeffs(substances[0].as_str(), 400.0)
            .unwrap();
        user_subs
            .extract_thermal_coeffs(substances[1].as_str(), 400.0)
            .unwrap();
        // Test for O2
        if let Ok((cp, dh, ds)) = user_subs.calculate_thermo_properties("O2", 400.0) {
            println!(
                "O2 Thermo properties at 400K: Cp={}, dH={}, dS={}",
                cp, dh, ds
            );
            assert!(cp > 0.0);
            assert!(ds > 0.0);
        } else {
            panic!("Failed to calculate O2 properties");
        }

        // Test for CO
        if let Ok((cp, dh, ds)) = user_subs.calculate_thermo_properties("CO", 400.0) {
            println!(
                "CO Thermo properties at 400K: Cp={}, dH={}, dS={}",
                cp, dh, ds
            );
            assert!(cp > 0.0);
            assert!(ds > 0.0);
        } else {
            panic!("Failed to calculate CO properties");
        }
    }

    #[test]
    fn test_transport_property_calculation() {
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["H2O".to_string(), "CO".to_string()];
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();
        user_subs.set_M(
            HashMap::from([("H2O".to_string(), 18.0), ("CO".to_string(), 32.0)]),
            None,
        );
        user_subs.set_P(1e5, None);
        user_subs.extract_transport_coeffs("H2O", 400.0).unwrap();
        user_subs.extract_transport_coeffs("CO", 400.0).unwrap();
        // Test for substances that should have transport data
        for substance in &["H2O", "CO"] {
            if let Ok((lambda, viscosity)) =
                user_subs.calculate_transport_properties(substance, 400.0, Some(33.0), None)
            {
                println!(
                    "{} Transport properties at 400K: Lambda={}, Viscosity={}",
                    substance, lambda, viscosity
                );
                assert!(lambda > 0.0);
                assert!(viscosity > 0.0);
            }
        }
    }

    #[test]
    fn test_property_maps_calculation() {
        let mut user_subs = SubsData::new();
        let substances = vec!["H2O".to_string(), "CO".to_string()];
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        user_subs.substances = substances.clone();

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();
        user_subs
            .extract_thermal_coeffs(substances[0].as_str(), 400.0)
            .unwrap();
        user_subs
            .extract_thermal_coeffs(substances[1].as_str(), 400.0)
            .unwrap();

        // Calculate property values at 400K
        user_subs
            .calculate_therm_map_of_properties(500.0)
            .expect("Failed to calculate property map");

        // Check that values were stored
        for substance in &["H2O", "CO"] {
            let property_map = user_subs.therm_map_of_properties_values.get(*substance);
            assert!(property_map.is_some());

            let property_map = property_map.unwrap();
            assert!(property_map.get(&DataType::Cp).unwrap().is_some());
            assert!(property_map.get(&DataType::dH).unwrap().is_some());
            assert!(property_map.get(&DataType::dS).unwrap().is_some());

            // Check that values are reasonable
            assert!(property_map.get(&DataType::Cp).unwrap().unwrap() > 0.0);
        }
    }

    #[test]
    fn test_function_maps_calculation() {
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["H2O".to_string()];

        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();
        user_subs.extract_thermal_coeffs("H2O", 400.0).unwrap();

        // Calculate function closures
        user_subs
            .calculate_therm_map_of_fun()
            .expect("Failed to calculate function map");

        // Check that functions were stored and work
        let cp_fun = user_subs.get_thermo_function("H2O", DataType::Cp_fun);
        assert!(cp_fun.is_some());

        // Test the function at different temperatures
        let cp_fun = cp_fun.unwrap();
        let cp_300 = cp_fun(300.0);
        let cp_400 = cp_fun(400.0);
        let cp_500 = cp_fun(500.0);

        println!(
            "H2O Cp values: 300K={}, 400K={}, 500K={}",
            cp_300, cp_400, cp_500
        );

        // Cp should be positive and should change with temperature
        assert!(cp_300 > 0.0);
        assert!(cp_400 > 0.0);
        assert!(cp_500 > 0.0);
        assert!(cp_300 != cp_400); // Values should differ with temperature
    }

    #[test]
    fn test_symbolic_maps_calculation() {
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["H2O".to_string()];
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();
        user_subs.extract_thermal_coeffs("H2O", 400.0).unwrap();
        // Calculate symbolic expressions
        user_subs
            .calculate_therm_map_of_sym()
            .expect("Failed to calculate symbolic map");

        // Check that symbolic expressions were stored
        let cp_sym = user_subs.get_thermo_symbolic("H2O", DataType::Cp_sym);
        assert!(cp_sym.is_some());

        // Check that the expression is not empty
        let cp_sym = cp_sym.unwrap();
        println!("H2O Cp symbolic expression: {}", cp_sym.to_string());
        assert!(!cp_sym.to_string().is_empty());
    }

    #[test]
    fn test_integration_all_features() {
        let mut user_subs = SubsData::new();
        let substances = vec!["H2O".to_string(), "CO".to_string(), "CH4".to_string()];
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("CH4".to_string(), Some(Phases::Gas));
        user_subs.substances = substances.clone();

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["NIST".to_string(), "Aramco_transport".to_string()],
            LibraryPriority::Permitted,
        );

        // Perform search
        user_subs.search_substances();
        user_subs
            .extract_thermal_coeffs(substances[0].as_str(), 400.0)
            .unwrap();
        user_subs
            .extract_thermal_coeffs(substances[1].as_str(), 400.0)
            .unwrap();
        user_subs
            .extract_thermal_coeffs(substances[2].as_str(), 400.0)
            .unwrap();
        // Calculate all types of data
        user_subs
            .calculate_therm_map_of_properties(400.0)
            .expect("Failed to calculate property map");
        user_subs
            .calculate_therm_map_of_fun()
            .expect("Failed to calculate function map");
        user_subs
            .calculate_therm_map_of_sym()
            .expect("Failed to calculate symbolic map");

        // Print search summary
        user_subs.print_search_summary();

        // Verify that we have data for all substances
        for substance in &["H2O", "CO", "CH4"] {
            // Check if found in any library
            let result = user_subs.get_substance_result(substance);
            assert!(result.is_some());

            // If found, check that we have property values
            if let Some(_map) = result {
                let property_map = user_subs.therm_map_of_properties_values.get(*substance);
                assert!(property_map.is_some());

                // Check that we have at least one function
                let function_map = user_subs.therm_map_of_fun.get(*substance);
                assert!(function_map.is_some());

                // Check that we have at least one symbolic expression
                let sym_map = user_subs.therm_map_of_sym.get(*substance);
                assert!(sym_map.is_some());
            }
        }
    }
}

#[cfg(test)]
mod tests2 {

    use crate::Thermodynamics::User_substances::{
        DataType, LibraryPriority, Phases, SubsData, WhatIsFound,
    };
    use approx::assert_relative_eq;
    use core::panic;
    use std::collections::HashMap;
    // Helper function to create a test SubsData instance with common setup
    fn create_test_subsdata() -> SubsData {
        let mut user_subs = SubsData::new();

        // Set up test substances
        user_subs.substances = vec!["H2O".to_string(), "CO2".to_string(), "CH4".to_string()];
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("CO2".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("CH4".to_string(), Some(Phases::Gas));

        // Set library priorities
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

        // Set pressure and molar mass
        user_subs.set_P(1.0, Some("atm".to_string()));

        let molar_masses = HashMap::from([
            ("H2O".to_string(), 18.01528),
            ("CO2".to_string(), 44.01),
            ("CH4".to_string(), 16.04),
        ]);
        user_subs.set_M(molar_masses, Some("g/mol".to_string()));

        user_subs
    }

    #[test]
    fn test_new() {
        let user_subs = SubsData::new();
        assert!(user_subs.substances.is_empty());
        assert!(user_subs.library_priorities.is_empty());
        assert!(user_subs.search_results.is_empty());
    }

    #[test]
    fn test_set_library_priority() {
        let mut user_subs = SubsData::new();

        // Set a single library priority
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        assert_eq!(user_subs.library_priorities.len(), 1);
        assert_eq!(
            *user_subs.library_priorities.get("NASA_gas").unwrap(),
            LibraryPriority::Priority
        );

        // Override existing priority
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Permitted);
        assert_eq!(
            *user_subs.library_priorities.get("NASA_gas").unwrap(),
            LibraryPriority::Permitted
        );
    }

    #[test]
    fn test_set_multiple_library_priorities() {
        let mut user_subs = SubsData::new();

        // Set multiple library priorities at once
        let libraries = vec![
            "NASA_gas".to_string(),
            "NUIG_thermo".to_string(),
            "Cantera_nasa_base_gas".to_string(),
        ];
        user_subs.set_multiple_library_priorities(libraries.clone(), LibraryPriority::Priority);

        // Verify all libraries were set
        for lib in libraries {
            assert_eq!(
                *user_subs.library_priorities.get(&lib).unwrap(),
                LibraryPriority::Priority
            );
        }
    }

    #[test]
    fn test_search_substances() {
        let mut user_subs = create_test_subsdata();

        // Perform the search
        user_subs.search_substances();

        // Check that results were populated
        assert!(!user_subs.search_results.is_empty());

        // At least some substances should be found
        let not_found = user_subs.get_not_found_substances();
        assert!(not_found.len() < user_subs.substances.len());
    }

    #[test]
    fn test_get_substance_result() {
        let mut user_subs = create_test_subsdata();
        user_subs.search_substances();

        // Get result for a specific substance
        let result = user_subs.get_substance_result("H2O");
        assert!(result.is_some());
    }

    #[test]
    fn test_get_not_found_substances() {
        let mut user_subs = SubsData::new();

        // Add a substance that likely doesn't exist in any library
        user_subs.substances = vec!["NonExistentSubstance123".to_string()];
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);

        user_subs.search_substances();

        let not_found = user_subs.get_not_found_substances();
        assert_eq!(not_found.len(), 1);
        assert_eq!(not_found[0], "NonExistentSubstance123");
    }

    #[test]
    fn test_extract_thermal_coeffs() {
        let mut user_subs = create_test_subsdata();
        user_subs.search_substances();

        // Find a substance that was found in the search
        let found_substances: Vec<String> = user_subs
            .search_results
            .keys()
            .filter(|&s| {
                if let Some(result_map) = user_subs.get_substance_result(s) {
                    result_map.get(&WhatIsFound::Thermo).is_some()
                        && result_map.get(&WhatIsFound::Thermo).unwrap().is_some()
                } else {
                    false
                }
            })
            .cloned()
            .collect();

        if !found_substances.is_empty() {
            let substance = &found_substances[0];
            let result = user_subs.extract_thermal_coeffs(substance, 298.15);
            assert!(result.is_ok());
        }
    }

    #[test]
    fn test_calculate_thermo_properties() {
        let mut user_subs = create_test_subsdata();
        user_subs.search_substances();

        // Find a substance that was found in the search
        let found_substances: Vec<String> = user_subs
            .search_results
            .keys()
            .filter(|&s| {
                if let Some(result_map) = user_subs.get_substance_result(s) {
                    result_map.get(&WhatIsFound::Thermo).is_some()
                        && result_map.get(&WhatIsFound::Thermo).unwrap().is_some()
                } else {
                    false
                }
            })
            .cloned()
            .collect();

        if !found_substances.is_empty() {
            let substance = &found_substances[0];

            // Extract coefficients first
            let _ = user_subs.extract_thermal_coeffs(substance, 298.15);

            // Calculate properties
            let result = user_subs.calculate_thermo_properties(substance, 298.15);
            assert!(result.is_ok());

            let (cp, _dh, _ds) = result.unwrap();
            assert!(cp > 0.0); // Heat capacity should be positive
        }
    }

    #[test]
    fn test_calculate_therm_map_of_properties() {
        let mut user_subs = create_test_subsdata();
        user_subs.search_substances();

        // Extract coefficients for all substances
        let _ = user_subs.extract_all_thermal_coeffs(298.15);

        // Calculate properties map
        let result = user_subs.calculate_therm_map_of_properties(298.15);
        assert!(result.is_ok());

        // Check that the map was populated
        assert!(!user_subs.therm_map_of_properties_values.is_empty());

        // Check values for found substances
        for substance in user_subs.get_priority_found_substances() {
            if let Some(props) = user_subs.therm_map_of_properties_values.get(&substance) {
                if let Some(Some(cp)) = props.get(&DataType::Cp) {
                    assert!(*cp > 0.0); // Heat capacity should be positive
                }
            }
        }
    }

    #[test]
    fn test_calculate_therm_map_of_fun() {
        let mut user_subs = create_test_subsdata();
        user_subs.search_substances();
        user_subs.print_search_summary();

        // Extract coefficients for all substances
        let _ = user_subs.extract_all_thermal_coeffs(298.15);

        // Calculate function closures
        let result = user_subs.calculate_therm_map_of_fun();
        assert!(result.is_ok());

        // Check that the map was populated
        assert!(!user_subs.therm_map_of_fun.is_empty());

        // Test a function for a found substance
        for substance in user_subs.get_priority_found_substances() {
            if let Some(fun_map) = user_subs.therm_map_of_fun.get(&substance) {
                if let Some(Some(cp_fun)) = fun_map.get(&DataType::Cp_fun) {
                    let cp_300 = cp_fun(300.0);
                    let cp_400 = cp_fun(400.0);

                    assert!(cp_300 > 0.0);
                    assert!(cp_400 > 0.0);
                    // Heat capacity typically increases with temperature
                    // (though not always, so we don't assert this)
                }
            }
        }
    }

    #[test]
    fn test_calculate_therm_map_of_sym() {
        //
        let mut user_subs = create_test_subsdata();
        user_subs.search_substances();

        // Extract coefficients for all substances
        let _ = user_subs.extract_all_thermal_coeffs(298.15);

        // Calculate symbolic expressions
        let result = user_subs.calculate_therm_map_of_sym();
        assert!(result.is_ok());

        // Check that the map was populated
        assert!(!user_subs.therm_map_of_sym.is_empty());
    }

    #[test]
    fn test_transport_properties() {
        //
        let mut user_subs = create_test_subsdata();
        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        user_subs.search_substances();

        // Extract coefficients for all substances
        let _ = user_subs.extract_all_thermal_coeffs(298.15);
        let _ = user_subs.extract_all_transport_coeffs(298.15);

        // Calculate thermal properties first (needed for transport)
        let _ = user_subs.calculate_therm_map_of_properties(298.15);

        // Find a substance that has transport data
        let found_substances: Vec<String> = user_subs
            .search_results
            .keys()
            .filter(|&s| {
                if let Some(result_map) = user_subs.get_substance_result(s) {
                    result_map.get(&WhatIsFound::Transport).is_some()
                        && result_map.get(&WhatIsFound::Transport).unwrap().is_some()
                } else {
                    false
                }
            })
            .cloned()
            .collect();

        if !found_substances.is_empty() {
            let substance = &found_substances[0];

            // Get Cp from thermal properties
            let cp = user_subs
                .therm_map_of_properties_values
                .get(substance)
                .and_then(|props| props.get(&DataType::Cp))
                .and_then(|opt_cp| *opt_cp);

            if let Some(cp) = cp {
                // Calculate transport properties
                let result =
                    user_subs.calculate_transport_properties(substance, 298.15, Some(cp), None);

                if result.is_ok() {
                    let (lambda, visc) = result.unwrap();
                    assert!(lambda > 0.0); // Thermal conductivity should be positive
                    assert!(visc > 0.0); // Viscosity should be positive
                }
            }
        }
    }

    #[test]
    fn test_print_search_summary() {
        let mut user_subs = create_test_subsdata();
        user_subs.search_substances();

        // This just tests that the function runs without errors
        user_subs.print_search_summary();
    }

    #[test]
    fn test_get_priority_and_permitted_found_substances() {
        let mut user_subs = create_test_subsdata();

        // Add more libraries with different priorities
        user_subs.set_library_priority("NASA_cond".to_string(), LibraryPriority::Permitted);
        user_subs.set_library_priority(
            "Cantera_nasa_base_gas".to_string(),
            LibraryPriority::Priority,
        );

        user_subs.search_substances();

        let priority_found = user_subs.get_priority_found_substances();
        let permitted_found = user_subs.get_permitted_found_substances();

        // The sum of priority found, permitted found, and not found should equal total substances
        assert_eq!(
            priority_found.len()
                + permitted_found.len()
                + user_subs.get_not_found_substances().len(),
            user_subs.substances.len()
        );
    }

    #[test]
    fn test_clone() {
        let mut user_subs = create_test_subsdata();
        user_subs.search_substances();

        // Extract coefficients and calculate properties
        let _ = user_subs.extract_all_thermal_coeffs(298.15);
        let _ = user_subs.calculate_therm_map_of_properties(298.15);
        let _ = user_subs.calculate_therm_map_of_fun();

        // Clone the instance
        let cloned_subs = user_subs.clone();

        // Check that basic properties were cloned correctly
        assert_eq!(cloned_subs.substances, user_subs.substances);
        assert_eq!(cloned_subs.library_priorities, user_subs.library_priorities);
        assert_eq!(
            cloned_subs.therm_map_of_properties_values.keys().count(),
            user_subs.therm_map_of_properties_values.keys().count()
        );

        // Check that function maps have the same structure
        assert_eq!(
            cloned_subs.therm_map_of_fun.keys().count(),
            user_subs.therm_map_of_fun.keys().count()
        );
    }
    #[test]
    fn test_if_not_found_go_NIST() {
        let mut user_subs = SubsData::new();

        // Set up test substances
        user_subs.substances = vec![
            "H2O".to_string(),
            "CO2".to_string(),
            "CH4".to_string(),
            "NH4ClO4".to_string(),
        ];
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("CO2".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("CH4".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("NH4ClO4".to_string(), Some(Phases::Solid));
        // Set library priorities
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

        // Set pressure and molar mass
        user_subs.set_P(1.0, Some("atm".to_string()));
        user_subs.search_substances();
        user_subs.if_not_found_go_NIST().unwrap();
        println!("\n\n Search results: {:#?}", user_subs.search_results);
        assert!(user_subs.search_results.contains_key("NH4ClO4"));
        assert!(user_subs.search_results.get("NH4ClO4").is_some());
        assert!(
            user_subs
                .search_results
                .get("NH4ClO4")
                .unwrap()
                .contains_key(&WhatIsFound::Thermo)
        );
        let not_found = user_subs.get_not_found_substances();
        println!("Not found: {:#?}", not_found);
        //panic!("Test failed");
    }
    // Test the calculate_elem_composition_and_molar_mass function
    #[test]
    fn test_calculate_elem_composition_and_molar_mass() {
        // Create a test instance with some substances
        let mut user_subs = create_test_subsdata();

        // Add specific substances with known compositions
        user_subs.substances = vec!["H2O".to_string(), "CO2".to_string(), "CH4".to_string()];

        // Search for substances in the libraries
        user_subs.search_substances();

        // Calculate element composition and molar mass
        let result = user_subs.calculate_elem_composition_and_molar_mass(None);
        assert!(result.is_ok());

        // Check that the element composition matrix was created
        assert!(user_subs.elem_composition_matrix.is_some());

        // Check that molar masses were calculated correctly
        assert!(user_subs.hasmap_of_molar_mass.contains_key("H2O"));
        assert!(user_subs.hasmap_of_molar_mass.contains_key("CO2"));
        assert!(user_subs.hasmap_of_molar_mass.contains_key("CH4"));

        // Check approximate molar mass values
        let h2o_mass = user_subs.hasmap_of_molar_mass.get("H2O").unwrap();
        let co2_mass = user_subs.hasmap_of_molar_mass.get("CO2").unwrap();
        let ch4_mass = user_subs.hasmap_of_molar_mass.get("CH4").unwrap();

        assert_relative_eq!(*h2o_mass, 18.0, epsilon = 0.1);
        assert_relative_eq!(*co2_mass, 44.0, epsilon = 0.1);
        assert_relative_eq!(*ch4_mass, 16.0, epsilon = 0.1);

        // Check the dimensions of the element composition matrix
        let matrix = user_subs.elem_composition_matrix.as_ref().unwrap();
        assert_eq!(matrix.nrows(), 3); // 3 substances
        assert!(matrix.ncols() >= 3); // At least 3 elements (C, H, O)

        // Test with custom chemical groups
        let mut groups = HashMap::new();
        let mut methyl_group = HashMap::new();
        methyl_group.insert("C".to_string(), 1);
        methyl_group.insert("H".to_string(), 3);
        groups.insert("Me".to_string(), methyl_group);

        // Create a new test instance with substances containing the custom group
        let mut user_subs_with_groups = create_test_subsdata();
        user_subs_with_groups.substances = vec!["MeOH".to_string()]; // Methanol

        // Search for substances
        user_subs_with_groups.search_substances();

        // Calculate with custom groups
        let result = user_subs_with_groups.calculate_elem_composition_and_molar_mass(Some(groups));
        assert!(result.is_ok());

        // Check that molar mass for methanol is correct
        if let Some(meoh_mass) = user_subs_with_groups.hasmap_of_molar_mass.get("MeOH") {
            assert_relative_eq!(*meoh_mass, 32.0, epsilon = 0.1); // CH3OH = 32 g/mol
        } else {
            panic!("Molar mass for MeOH not calculated");
        }
    }

    #[test]
    fn test_calculate_elem_composition_matrix_structure() {
        // Create a test instance with substances having known compositions
        let mut user_subs = create_test_subsdata();
        user_subs.substances = vec!["H2".to_string(), "O2".to_string(), "H2O".to_string()];

        // Search for substances
        user_subs.search_substances();

        // Calculate element composition and molar mass
        let result = user_subs.calculate_elem_composition_and_molar_mass(None);
        assert!(result.is_ok());

        // Get the matrix
        let matrix = user_subs.elem_composition_matrix.as_ref().unwrap();
        println!("Matrix: {}", matrix);

        // Check matrix dimensions
        assert_eq!(matrix.nrows(), 3); // 3 substances
        assert_eq!(matrix.ncols(), 2); // 2 elements (H, O)

        // Check specific values in the matrix
        // H2 should have 2 H atoms and 0 O atoms
        assert_eq!(matrix[(0, 0)], 2.0); // H in H2
        assert_eq!(matrix[(0, 1)], 0.0); // O in H2

        // O2 should have 0 H atoms and 2 O atoms
        assert_eq!(matrix[(1, 0)], 0.0); // H in O2
        assert_eq!(matrix[(1, 1)], 2.0); // O in O2

        // H2O should have 2 H atoms and 1 O atom
        assert_eq!(matrix[(2, 0)], 2.0); // H in H2O 
        assert_eq!(matrix[(2, 1)], 1.0); // O in H2O
    }

    #[test]
    fn test_calculate_elem_composition_with_complex_molecules() {
        // Create a test instance with more complex molecules
        let mut user_subs = create_test_subsdata();
        user_subs.substances = vec![
            "C6H12O6".to_string(),
            "C2H5OH".to_string(),
            "NH3".to_string(),
        ];

        // Search for substances
        user_subs.search_substances();

        // Calculate element composition and molar mass
        let result = user_subs.calculate_elem_composition_and_molar_mass(None);
        assert!(result.is_ok());

        // Check molar masses
        let glucose_mass = user_subs.hasmap_of_molar_mass.get("C6H12O6").unwrap();
        let ethanol_mass = user_subs.hasmap_of_molar_mass.get("C2H5OH").unwrap();
        let ammonia_mass = user_subs.hasmap_of_molar_mass.get("NH3").unwrap();

        assert_relative_eq!(*glucose_mass, 180.0, epsilon = 0.5); // C6H12O6 = 180 g/mol
        assert_relative_eq!(*ethanol_mass, 46.0, epsilon = 0.5); // C2H5OH = 46 g/mol
        assert_relative_eq!(*ammonia_mass, 17.0, epsilon = 0.5); // NH3 = 17 g/mol

        // Get the matrix
        let matrix = user_subs.elem_composition_matrix.as_ref().unwrap();

        // Check matrix dimensions
        assert_eq!(matrix.nrows(), 3); // 3 substances
        assert!(matrix.ncols() >= 4); // At least 4 elements (C, H, O, N)
    }
}
