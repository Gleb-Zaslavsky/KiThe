use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
use crate::Thermodynamics::DBhandlers::transport_api::TransportCalculator;
use crate::Thermodynamics::User_substances::{CalculatorType, WhatIsFound};
use crate::Thermodynamics::User_substances::{DataType, LibraryPriority, SubsData};
use std::collections::HashMap;
pub fn collecting_thermo_data(thermotask: usize) {
    match thermotask {
        0 => {
            // Create a new SubsData instance with test substances
            let substances = vec!["CO".to_string(), "H2O".to_string()];
            let mut user_subs = SubsData::new();
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
            user_subs.extract_thermal_coeffs("CO", 400.0).unwrap();
            user_subs.extract_thermal_coeffs("H2O", 400.0).unwrap();

            // Test thermo calculations
            let (cp, dh, ds) = user_subs.calculate_thermo_properties("CO", 400.0).unwrap();
            println!("CO Thermo properties at 400K:");
            println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
            let (cp, dh, ds) = user_subs.calculate_thermo_properties("H2O", 400.0).unwrap();
            println!("H2O Thermo properties at 400K:");
            println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
            let Cp = 33.8;
            // Test transport calculations
            println!("\n \n Transport properties:");
            user_subs.set_M(
                HashMap::from([("H2O".to_string(), 18.0), ("CO".to_string(), 32.0)]),
                None,
            );
            user_subs.set_P(1e5, None);
            user_subs
                .extract_thermal_coeffs(substances[0].as_str(), 400.0)
                .unwrap();
            user_subs
                .extract_thermal_coeffs(substances[1].as_str(), 400.0)
                .unwrap();
            let (lambda, viscosity) = user_subs
                .calculate_transport_properties("H2O", 400.0, Some(Cp), None)
                .unwrap();
            println!("H2O Transport properties at 400K:");
            println!("Lambda: {}, Viscosity: {}", lambda, viscosity);
            // Test transport calculations

            // Print full summary
            user_subs.print_search_summary();
        }
        1 => {
            let mut user_subs = SubsData::new();
            let substances = vec!["H2O".to_string(), "CO".to_string()];
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
            user_subs.print_search_summary();
        }
        2 => {
            // Create a new SubsData instance with test substances
            let substances = vec!["CO".to_string(), "H2O".to_string()];
            let mut user_subs = SubsData::new();
            user_subs.substances = substances;

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
            user_subs.extract_thermal_coeffs("CO", 400.0).unwrap();
            user_subs.extract_thermal_coeffs("H2O", 400.0).unwrap();
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
            user_subs.print_search_summary();
        }
        _ => {
            println!("Invalid task number. Please choose from 1 or 2.");
        }
    }
}
