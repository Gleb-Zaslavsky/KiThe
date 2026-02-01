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
        user_subs.search_substances().unwrap();
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
        user_subs.search_substances().unwrap();
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
        let _ = user_subs.search_substances().unwrap();

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
        let _ = user_subs.search_substances().unwrap();
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
        let _ = user_subs.search_substances().unwrap();
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
        let _ = user_subs.search_substances().unwrap();
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
    fn test_search_by_elements() {
        let mut user_subs = SubsData::new();

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Search for substances containing carbon and oxygen
        let elements = vec!["C".to_string(), "O".to_string()];
        let result = user_subs.search_by_elements(elements);

        assert!(result.is_ok());
        let found_substances = result.unwrap();

        // Should find substances like CO, CO2, etc.
        assert!(!found_substances.is_empty());

        // Check that search results were populated
        for substance in &found_substances {
            let search_result = user_subs.get_substance_result(substance);
            assert!(search_result.is_some());

            let result_map = search_result.unwrap();
            assert!(result_map.contains_key(&WhatIsFound::Thermo));
        }
    }

    #[test]
    fn test_search_by_elements_only() {
        let mut user_subs = SubsData::new();

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Search for substances containing only hydrogen and oxygen
        let elements = vec!["H".to_string(), "O".to_string()];
        let result = user_subs.search_by_elements_only(elements);

        assert!(result.is_ok());
        let found_substances = result.unwrap();
        let _ = user_subs.calculate_therm_map_of_properties(400.0);
        // Should find substances like H2O, H2O2 but not CH4, CO2, etc.
        if !found_substances.is_empty() {
            // Check that search results were populated
            for substance in &found_substances {
                let search_result = user_subs.get_substance_result(substance);
                assert!(search_result.is_some());

                let result_map = search_result.unwrap();
                assert!(result_map.contains_key(&WhatIsFound::Thermo));
                let property_map = user_subs.therm_map_of_properties_values.get(substance);
                assert!(property_map.is_some());

                let property_map = property_map.unwrap();
                assert!(property_map.get(&DataType::Cp).unwrap().is_some());
                assert!(property_map.get(&DataType::dH).unwrap().is_some());
                assert!(property_map.get(&DataType::dS).unwrap().is_some());

                // Check that values are reasonable
                // assert!(property_map.get(&DataType::Cp).unwrap().unwrap() > 0.0);
            }
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
        let _ = user_subs.search_substances().unwrap();
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
        let _ = user_subs.search_substances().unwrap();
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
        let _ = user_subs.search_substances().unwrap();
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
        let _ = user_subs.search_substances().unwrap();

        // Check that results were populated
        assert!(!user_subs.search_results.is_empty());

        // At least some substances should be found
        let not_found = user_subs.get_not_found_substances();
        assert!(not_found.len() < user_subs.substances.len());
    }

    #[test]
    fn test_get_substance_result() {
        let mut user_subs = create_test_subsdata();
        let _ = user_subs.search_substances().unwrap();

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

        let _ = user_subs.search_substances().unwrap();

        let not_found = user_subs.get_not_found_substances();
        assert_eq!(not_found.len(), 1);
        assert_eq!(not_found[0], "NonExistentSubstance123");
    }

    #[test]
    fn test_extract_thermal_coeffs() {
        let mut user_subs = create_test_subsdata();
        let _ = user_subs.search_substances().unwrap();

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
        let _ = user_subs.search_substances().unwrap();

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
        let _ = user_subs.search_substances().unwrap();

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
        let _ = user_subs.search_substances().unwrap();
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
        let _ = user_subs.search_substances().unwrap();

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

        let _ = user_subs.search_substances().unwrap();

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
        let _ = user_subs.search_substances().unwrap();

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

        user_subs.search_substances().unwrap();

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
        let _ = user_subs.search_substances().unwrap();

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
        let _ = user_subs.search_substances().unwrap();
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
        let _ = user_subs.search_substances().unwrap();

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
        let _ = user_subs_with_groups.search_substances().unwrap();

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
        user_subs.search_substances().unwrap();

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
        let _ = user_subs.search_substances().unwrap();

        // Calculate element composition and molar mass
        let result = user_subs.calculate_elem_composition_and_molar_mass(None);
        assert!(result.is_ok());
        assert!(user_subs.elem_composition_matrix.is_some());
        assert!(!user_subs.hasmap_of_molar_mass.is_empty());
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

    #[test]
    fn test_calculate_transport_map_of_properties() {
        let mut user_subs = SubsData::new();
        let substances = vec!["H2O".to_string(), "CO".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));

        // Set library priorities for both thermo and transport
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Search substances
        let _ = user_subs.search_substances().unwrap();

        // Calculate molar masses (required prerequisite)
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();
        assert!(user_subs.elem_composition_matrix.is_some());
        assert!(!user_subs.hasmap_of_molar_mass.is_empty());
        println!("Molar masses: {:?}", user_subs.hasmap_of_molar_mass);
        // Set pressure (required for transport calculations)
        user_subs.set_P(1e5, None);

        // Extract thermal coefficients
        user_subs.extract_all_thermal_coeffs(400.0).unwrap();

        // Calculate thermal properties (required prerequisite for Cp values)
        user_subs.calculate_therm_map_of_properties(400.0).unwrap();

        // Extract transport coefficients
        user_subs.extract_all_transport_coeffs(400.0).unwrap();

        // Test the main function
        let result = user_subs.calculate_transport_map_of_properties(400.0);
        assert!(result.is_ok());

        // Verify transport properties were calculated and stored
        for substance in &substances {
            if let Some(transport_props) =
                user_subs.transport_map_of_properties_values.get(substance)
            {
                if let Some(Some(lambda)) = transport_props.get(&DataType::Lambda) {
                    println!("{} Lambda: {}", substance, lambda);
                    assert!(*lambda > 0.0, "Lambda should be positive for {}", substance);
                }
                if let Some(Some(visc)) = transport_props.get(&DataType::Visc) {
                    println!("{} Viscosity: {}", substance, visc);
                    assert!(
                        *visc > 0.0,
                        "Viscosity should be positive for {}",
                        substance
                    );
                }
            }
        }
    }

    #[test]
    fn test_calculate_transport_map_of_properties_CEA() {
        // Create a test instance with substances having known compositions
        let mut user_subs = SubsData::new();
        let substances = vec!["H2O".to_string(), "CO".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));

        // Set library priorities for both thermo and transport
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs
            .set_multiple_library_priorities(vec!["CEA".to_string()], LibraryPriority::Priority);

        // Search substances
        let _ = user_subs.search_substances().unwrap();

        // Calculate molar masses (required prerequisite)
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();
        assert!(user_subs.elem_composition_matrix.is_some());
        assert!(!user_subs.hasmap_of_molar_mass.is_empty());
        println!("Molar masses: {:?}", user_subs.hasmap_of_molar_mass);
        // Set pressure (required for transport calculations)
        user_subs.set_P(1e5, None);
        //   println!("instance before clone: {:#?}", user_subs);
        let search_results = &user_subs.search_results;
        let H2O_results = search_results.get("H2O");
        assert!(H2O_results.is_some());
        let H2O_results = H2O_results.unwrap();
        assert!(H2O_results.get(&WhatIsFound::Thermo).is_some());
        assert!(H2O_results.get(&WhatIsFound::Transport).is_some());
        let H2O_thermo = H2O_results.get(&WhatIsFound::Thermo).unwrap();
        assert!(H2O_thermo.is_some());
        let H2O_transport = H2O_results.get(&WhatIsFound::Transport).unwrap();
        assert!(H2O_transport.is_some());

        // Extract thermal coefficients
        user_subs.extract_all_thermal_coeffs(400.0).unwrap();

        // Calculate thermal properties (required prerequisite for Cp values)
        user_subs.calculate_therm_map_of_properties(400.0).unwrap();

        // Extract transport coefficients
        user_subs.extract_all_transport_coeffs(400.0).unwrap();

        // Test the main function
        let result = user_subs.calculate_transport_map_of_properties(400.0);
        assert!(result.is_ok());

        // Verify transport properties were calculated and stored
        for substance in &substances {
            if let Some(transport_props) =
                user_subs.transport_map_of_properties_values.get(substance)
            {
                if let Some(Some(lambda)) = transport_props.get(&DataType::Lambda) {
                    println!("{} Lambda: {}", substance, lambda);
                    assert!(*lambda > 0.0, "Lambda should be positive for {}", substance);
                }
                if let Some(Some(visc)) = transport_props.get(&DataType::Visc) {
                    println!("{} Viscosity: {}", substance, visc);
                    assert!(
                        *visc > 0.0,
                        "Viscosity should be positive for {}",
                        substance
                    );
                }
            }
        }
    }
}

// Include error handling tests

#[cfg(test)]
mod additional_error_tests {

    use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};
    use crate::Thermodynamics::User_substances_error::SubsDataError;

    #[test]
    fn test_absurdly_high_temperatures() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances().unwrap();

        // Test temperatures beyond physical limits
        let extreme_temps = vec![1_000_000.0, 10_000_000.0, f64::MAX];
        for temp in extreme_temps {
            let result = subs_data.extract_thermal_coeffs("CO", temp);
            // Should either be invalid temperature or coefficient extraction failed
            if let Err(error) = result {
                assert!(matches!(
                    error,
                    SubsDataError::InvalidTemperature(_)
                        | SubsDataError::CoefficientExtractionFailed { .. }
                        | SubsDataError::ThermoError(_)
                ));
            }
        }
    }

    #[test]
    fn test_fictional_substances() {
        let mut subs_data = SubsData::new();
        let fictional_substances = vec![
            "Unobtainium".to_string(),
            "Vibranium".to_string(),
            "Kryptonite".to_string(),
            "Element115".to_string(),
        ];
        subs_data.substances = fictional_substances.clone();
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances().unwrap();

        // All should result in SubstanceNotFound or CalculatorNotAvailable errors
        for substance in &fictional_substances {
            let result = subs_data.extract_thermal_coeffs(substance, 400.0);
            println!("Result for {}: {:?}", substance, result);
            assert!(matches!(
                result,
                Err(SubsDataError::SubstanceNotFound(_))
                    | Err(SubsDataError::CalculatorNotAvailable { .. })
            ));

            let result = subs_data.calculate_thermo_properties(substance, 400.0);
            assert!(matches!(
                result,
                Err(SubsDataError::SubstanceNotFound(_))
                    | Err(SubsDataError::CalculatorNotAvailable { .. })
            ));
        }
    }

    #[test]
    fn test_malformed_substance_names() {
        let mut subs_data = SubsData::new();
        let malformed_names = vec![
            "".to_string(),                                              // Empty string
            "   ".to_string(),                                           // Whitespace only
            "CO@#$%".to_string(),                                        // Special characters
            "123456".to_string(),                                        // Numbers only
            "VeryLongSubstanceNameThatDoesNotExistAnywhere".to_string(), // Very long name
        ];
        subs_data.substances = malformed_names.clone();
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances().unwrap();

        for substance in &malformed_names {
            let result = subs_data.extract_thermal_coeffs(substance, 400.0);
            assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));
        }
    }

    #[test]
    fn test_transport_without_prerequisites() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );
        subs_data.search_substances().unwrap();

        // Try to calculate transport properties without setting pressure
        let result = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);
        assert!(matches!(result, Err(SubsDataError::PressureNotSet)));

        // Set pressure but not molar mass
        subs_data.set_P(101325.0, None);
        let result = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);
        assert!(matches!(result, Err(SubsDataError::MolarMassNotFound(_))));
    }

    #[test]
    fn test_batch_operations_with_mixed_results() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec![
            "CO".to_string(),
            "H2O".to_string(),
            "NonExistent1".to_string(),
            "NonExistent2".to_string(),
        ];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances().unwrap();

        // Batch operation should fail on first non-existent substance
        let result = subs_data.extract_all_thermal_coeffs(400.0);
        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));

        // But individual operations should work for existing substances
        assert!(subs_data.extract_thermal_coeffs("CO", 400.0).is_ok());
        assert!(subs_data.extract_thermal_coeffs("H2O", 400.0).is_ok());
    }

    #[test]
    fn test_error_message_content() {
        // Test that error messages contain useful information
        let error = SubsDataError::CoefficientExtractionFailed {
            substance: "TestSubstance".to_string(),
            temperature: Some(5000.0),
        };
        let message = format!("{}", error);
        assert!(message.contains("TestSubstance"));
        assert!(message.contains("5000"));
        assert!(message.contains("extract coefficients"));

        let error = SubsDataError::FunctionCreationFailed {
            substance: "CO2".to_string(),
            function_type: "test function".to_string(),
        };
        let message = format!("{}", error);
        assert!(message.contains("CO2"));
        assert!(message.contains("test function"));
        assert!(message.contains("Failed to create"));

        let error = SubsDataError::SubstanceNotFound("UnknownSubstance".to_string());
        let message = format!("{}", error);
        assert!(message.contains("UnknownSubstance"));
        assert!(message.contains("not found"));

        let error = SubsDataError::InvalidTemperature(-100.0);
        let message = format!("{}", error);
        assert!(message.contains("-100"));
        assert!(message.contains("Invalid temperature"));
    }
    ///////////////////////////////////////
}

#[cfg(test)]
mod logging_tests {
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
    use crate::Thermodynamics::User_substances_error::ExceptionLogger;
    use std::collections::HashMap;

    #[test]
    fn test_logging_non_existent_substances() {
        use crate::Thermodynamics::User_substances_error::LogErrorType;

        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["FakeSubstance123".to_string(), "Unobtainium".to_string()];
        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        let _ = subs_data.search_substances();

        subs_data.clear_error_logs();

        let _ = subs_data.extract_thermal_coeffs("FakeSubstance123", 400.0);
        let _ = subs_data.calculate_thermo_properties("Unobtainium", 500.0);
        let _ = subs_data.extract_transport_coeffs("FakeSubstance123", 300.0);

        let logs = subs_data.logger.get_logs();
        assert!(!logs.is_empty());

        // Check that SubstanceNotFound errors are logged
        let substance_not_found_count = logs
            .iter()
            .filter(|log| matches!(log.error_type, LogErrorType::SubstanceNotFound(_)))
            .count();
        assert!(
            substance_not_found_count > 0,
            "Expected SubstanceNotFound errors"
        );

        println!("\n=== Non-existent Substances ===");
        subs_data.logger.print_pretty_table();
    }

    #[test]
    fn test_logging_invalid_temperatures() {
        use crate::Thermodynamics::User_substances_error::LogErrorType;

        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        let _ = subs_data.search_substances();

        subs_data.clear_error_logs();

        let invalid_temps = vec![-100.0, 0.0, -273.15];
        for temp in invalid_temps {
            let _ = subs_data.extract_thermal_coeffs("CO", temp);
            let _ = subs_data.extract_transport_coeffs("CO", temp);
        }

        let logs = subs_data.logger.get_logs();
        assert!(!logs.is_empty());

        // Check that InvalidTemperature errors are logged
        let invalid_temp_count = logs
            .iter()
            .filter(|log| matches!(log.error_type, LogErrorType::InvalidTemperature(_)))
            .count();
        assert!(invalid_temp_count > 0, "Expected InvalidTemperature errors");

        println!("\n=== Invalid Temperatures ===");
        subs_data.logger.print_pretty_table();
    }

    #[test]
    fn test_logging_batch_operations() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec![
            "CO".to_string(),
            "NonExistent1".to_string(),
            "CO2".to_string(),
        ];
        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        let _ = subs_data.search_substances();

        subs_data.clear_error_logs();

        let _ = subs_data.extract_all_thermal_coeffs(400.0);
        let _ = subs_data.calculate_therm_map_of_properties(400.0);
        let _ = subs_data.calculate_therm_map_of_fun();

        let logs = subs_data.logger.get_logs();
        assert!(!logs.is_empty());

        println!("\n=== Batch Operations ===");
        subs_data.logger.print_pretty_table();
    }

    #[test]
    fn test_logging_transport_prerequisites() {
        use crate::Thermodynamics::User_substances_error::LogErrorType;

        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_library_priority("Aramco_transport".to_string(), LibraryPriority::Priority);
        let _ = subs_data.search_substances();

        subs_data.clear_error_logs();

        let _ = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);

        subs_data.set_P(101325.0, None);
        let _ = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);

        let logs = subs_data.logger.get_logs();
        assert!(!logs.is_empty());

        // Check for PressureNotSet and MolarMassNotFound errors
        let pressure_not_set_count = logs
            .iter()
            .filter(|log| matches!(log.error_type, LogErrorType::PressureNotSet(_)))
            .count();
        let molar_mass_not_found_count = logs
            .iter()
            .filter(|log| matches!(log.error_type, LogErrorType::MolarMassNotFound(_)))
            .count();

        assert!(
            pressure_not_set_count > 0 || molar_mass_not_found_count > 0,
            "Expected PressureNotSet or MolarMassNotFound errors"
        );

        println!("\n=== Transport Prerequisites ===");
        subs_data.logger.print_pretty_table();
    }

    #[test]
    fn test_logging_transport_batch_operations() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["H2O".to_string(), "FakeGas".to_string(), "CO".to_string()];

        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        subs_data.set_library_priority("Aramco_transport".to_string(), LibraryPriority::Priority);

        for substance in &subs_data.substances {
            subs_data
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        let _ = subs_data.search_substances();
        subs_data.set_P(101325.0, None);
        subs_data.set_M(
            HashMap::from([
                ("H2O".to_string(), 18.0),
                ("CO".to_string(), 28.0),
                ("FakeGas".to_string(), 32.0),
            ]),
            None,
        );

        subs_data.clear_error_logs();

        let _ = subs_data.extract_all_transport_coeffs(400.0);
        let _ = subs_data.extract_all_thermal_coeffs(400.0);
        let _ = subs_data.calculate_therm_map_of_properties(400.0);
        let _ = subs_data.calculate_transport_map_of_properties(400.0);
        let _ = subs_data.calculate_transport_map_of_functions();
        let _ = subs_data.calculate_transport_map_of_sym();

        println!("\n=== Transport Batch Operations ===");
        subs_data.logger.print_pretty_table();
    }

    #[test]
    fn test_logging_extreme_conditions() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        let _ = subs_data.search_substances();

        subs_data.clear_error_logs();

        let extreme_temps = vec![10_000.0, 50_000.0, 100_000.0];
        for temp in extreme_temps {
            let _ = subs_data.extract_thermal_coeffs("CO", temp);
        }

        println!("\n=== Extreme Conditions ===");
        subs_data.logger.print_pretty_table();
    }

    #[test]
    fn test_comprehensive_error_coverage() {
        use crate::Thermodynamics::User_substances_error::LogErrorType;

        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string(), "FakeSubstance".to_string()];

        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        subs_data.set_library_priority("Aramco_transport".to_string(), LibraryPriority::Priority);

        let _ = subs_data.search_substances();
        subs_data.clear_error_logs();

        let _ = subs_data.extract_thermal_coeffs("FakeSubstance", 400.0);
        let _ = subs_data.extract_thermal_coeffs("CO", -100.0);
        let _ = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);

        subs_data.set_P(101325.0, None);
        let _ = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);

        subs_data.set_M(HashMap::from([("CO".to_string(), 28.0)]), None);
        let _ = subs_data.calculate_transport_properties("FakeSubstance", 400.0, Some(30.0), None);

        let logs = subs_data.logger.get_logs();
        println!("\n=== Comprehensive Error Coverage ===");
        println!("Total errors logged: {}", logs.len());
        subs_data.logger.print_pretty_table();

        // Check for specific error types
        let substance_not_found = logs
            .iter()
            .any(|log| matches!(log.error_type, LogErrorType::SubstanceNotFound(_)));
        let invalid_temperature = logs
            .iter()
            .any(|log| matches!(log.error_type, LogErrorType::InvalidTemperature(_)));
        let pressure_not_set = logs
            .iter()
            .any(|log| matches!(log.error_type, LogErrorType::PressureNotSet(_)));
        let molar_mass_not_found = logs
            .iter()
            .any(|log| matches!(log.error_type, LogErrorType::MolarMassNotFound(_)));

        assert!(substance_not_found, "Expected SubstanceNotFound error");
        assert!(invalid_temperature, "Expected InvalidTemperature error");
        assert!(
            pressure_not_set || molar_mass_not_found,
            "Expected PressureNotSet or MolarMassNotFound error"
        );

        let error_types: std::collections::HashSet<_> = logs
            .iter()
            .map(|log| {
                format!("{:?}", log.error_type)
                    .split('(')
                    .next()
                    .unwrap_or("Unknown")
                    .to_string()
            })
            .collect();

        println!("Error types: {:?}", error_types);
        assert!(error_types.len() > 1);
    }

    #[test]
    fn test_logger_functionality() {
        use crate::Thermodynamics::User_substances_error::LogErrorType;

        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["NonExistent".to_string()];
        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        let _ = subs_data.search_substances();

        subs_data.clear_error_logs();
        assert!(subs_data.logger.get_logs().is_empty());

        let _ = subs_data.extract_thermal_coeffs("NonExistent", 400.0);
        let _ = subs_data.extract_thermal_coeffs("NonExistent", -100.0);

        let logs = subs_data.logger.get_logs();
        assert_eq!(logs.len(), 2);

        // Verify specific error types
        assert!(
            logs.iter()
                .any(|log| matches!(log.error_type, LogErrorType::SubstanceNotFound(_)))
        );
        assert!(
            logs.iter()
                .any(|log| matches!(log.error_type, LogErrorType::InvalidTemperature(_)))
        );

        // Verify substance and function fields
        assert!(logs.iter().all(|log| log.substance == "NonExistent"));
        assert!(
            logs.iter()
                .all(|log| log.function == "extract_thermal_coeffs")
        );

        println!("\n=== Logger Functionality Test ===");
        subs_data.logger.print_pretty_table();

        subs_data.clear_error_logs();
        assert!(subs_data.logger.get_logs().is_empty());
    }

    #[test]
    fn test_with_calculator_logging() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        subs_data.set_library_priority("Aramco_transport".to_string(), LibraryPriority::Priority);
        let _ = subs_data.search_substances();

        subs_data.clear_error_logs();

        let _ = subs_data.set_T_range_for_thermo("CO", 10000.0, 5000.0);
        let _ = subs_data.fitting_thermal_coeffs_for_T_interval("CO");
        let _ = subs_data.parse_thermal_coeffs("CO");
        let _ = subs_data.integr_mean("CO");

        println!("\n=== With Calculator Logging ===");
        subs_data.logger.print_pretty_table();
    }
}
