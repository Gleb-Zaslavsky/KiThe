/// module for classical thermodynamics and chemical equilibrium
pub mod ChemEquilibrium;
#[allow(non_snake_case)]
/// handlers for different formats of thermodynamics and heat-mass transfer data
pub mod DBhandlers;
///
pub mod User_PhaseOrSolution;
/// heat-mass transfer data agregator
/// # Examples 1
/// ```
/// use KiThe::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
/// use KiThe::Thermodynamics::DBhandlers::transport_api::TransportCalculator;
/// use KiThe::Thermodynamics::User_substances::{CalculatorType, WhatIsFound};
/// use KiThe::Thermodynamics::User_substances::{DataType, LibraryPriority, SubsData};
/// use std::collections::HashMap;
///   // Create a new SubsData instance with test substances
///   let substances = vec!["CO".to_string(), "H2O".to_string()];
///   let mut user_subs = SubsData::new();
///   user_subs.substances = substances.clone();
///
///   // Set library priorities
///   user_subs.set_multiple_library_priorities(
///       vec!["NASA_gas".to_string()],
///       LibraryPriority::Priority,
///   );
///   user_subs.set_multiple_library_priorities(
///       vec!["Aramco_transport".to_string()],
///       LibraryPriority::Priority,
///   );
///
///   // Perform the search
///   user_subs.search_substances();
///   user_subs.extract_thermal_coeffs("CO", 400.0).unwrap();
///   user_subs.extract_thermal_coeffs("H2O", 400.0).unwrap();
///
///   // Test thermo calculations
///   let (cp, dh, ds) = user_subs.calculate_thermo_properties("CO", 400.0).unwrap();
///   println!("CO Thermo properties at 400K:");
///   println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
///   let (cp, dh, ds) = user_subs.calculate_thermo_properties("H2O", 400.0).unwrap();
///   println!("H2O Thermo properties at 400K:");
///   println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
///   let Cp = 33.8;
///   // Test transport calculations
///   println!("\n \n Transport properties:");
///   user_subs.set_M(
///       HashMap::from([("H2O".to_string(), 18.0), ("CO".to_string(), 32.0)]),
///       None,
///   );
///   user_subs.set_P(1e5, None);
///   user_subs
///       .extract_thermal_coeffs(substances[0].as_str(), 400.0)
///       .unwrap();
///   user_subs
///       .extract_thermal_coeffs(substances[1].as_str(), 400.0)
///       .unwrap();
///   let (lambda, viscosity) = user_subs
///       .calculate_transport_properties("H2O", 400.0, Some(Cp), None)
///       .unwrap();
///   println!("H2O Transport properties at 400K:");
///   println!("Lambda: {}, Viscosity: {}", lambda, viscosity);
///   // Test transport calculations
///
///   // Print full summary
///   user_subs.print_search_summary();
/// ```
/// heat-mass transfer data agregator
/// # Examples 2
/// ```
/// use KiThe::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
/// use KiThe::Thermodynamics::DBhandlers::transport_api::TransportCalculator;
/// use KiThe::Thermodynamics::User_substances::{CalculatorType, WhatIsFound};
/// use KiThe::Thermodynamics::User_substances::{DataType, LibraryPriority, SubsData};
/// use std::collections::HashMap;
///            // Create a new SubsData instance with test substances
///            let substances = vec!["CO".to_string(), "H2O".to_string()];
///            let mut user_subs = SubsData::new();
///            user_subs.substances = substances;
///
///            // Set library priorities
///            user_subs.set_multiple_library_priorities(
///                vec!["NASA_gas".to_string()],
///                LibraryPriority::Priority,
///            );
///            user_subs.set_multiple_library_priorities(
///                vec!["Aramco_transport".to_string()],
///                LibraryPriority::Priority,
///            );
///
///            // Perform the search
///            user_subs.search_substances();
///            user_subs.extract_thermal_coeffs("CO", 400.0).unwrap();
///            user_subs.extract_thermal_coeffs("H2O", 400.0).unwrap();
///            // Print full summary
///            user_subs.print_search_summary();
///            let datamap = user_subs.get_substance_result("CO").unwrap();
///            let Thermo = datamap.get(&WhatIsFound::Thermo).unwrap().as_ref().unwrap();
///            let Calculator = Thermo.calculator.as_ref().unwrap();
///            let Cp;
///            match Calculator {
///                CalculatorType::Thermo(thermo) => {
///                    // Test thermo calculations
///                    let mut thermo = thermo.clone();
///                    thermo.extract_model_coefficients(400.0).unwrap();
///                    if let Ok(()) = thermo.calculate_Cp_dH_dS(400.0) {
///                        let (cp, dh, ds) = (
///                            thermo.get_Cp().unwrap(),
///                            thermo.get_dh().unwrap(),
///                            thermo.get_ds().unwrap(),
///                        ); //thermo.get_Cp().unwrap();
///                        println!("CO Thermo properties at 400K:");
///                        println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
///                        Cp = cp;
///                        assert!(cp > 0.0);
///                    } else {
///                        panic!("Failed to calculate CO properties");
///                    }
///                }
///                _ => {
///                    panic!("Failed to calculate CO properties");
///                }
///            }
///            let Transport = datamap
///                .get(&WhatIsFound::Transport)
///                .unwrap()
///                .as_ref()
///                .unwrap();
///            let Calculator = Transport.calculator.as_ref().unwrap();
///            match Calculator {
///                CalculatorType::Transport(transport) => {
///                    // Test transport calculations
///                    let mut transport = transport.clone();
///                    let _ = transport.set_M(32.0, None);
///                    let _ = transport.set_P(1e5, None);
///                    transport.extract_coefficients(400.0).unwrap();
///                    if let Ok(L) = transport.calculate_lambda(Some(Cp), None, 400.0) {
///                        println!("CO Transport properties at 400K:");
///                        println!("Lambda: {}", L);
///                        assert!(L > 0.0);
///                    } else {
///                        panic!("Failed to calculate H2O properties");
///                   }
///                }
///                _ => {
///                    panic!("Failed to calculate CO properties");
///                }
///            }
///            user_subs.print_search_summary();
/// ```
pub mod User_substances;
pub mod User_substances2;
/// tests
pub mod User_substances_tests;
/// main functionality to open thermodynamics and heat-mass transfer libraries
pub mod thermo_lib_api;
/// calculations of thermodynamic properties, creating closures and symbolic expressions
pub mod thermo_properties_api;
/// calculations of heat-mass transfer properties, creating closures and symbolic expressions
pub mod transport_properties_api;

mod dG_dS;
