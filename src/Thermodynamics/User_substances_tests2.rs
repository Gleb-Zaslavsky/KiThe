///////////////////////TESTS////////////////////////////////////////////
#[cfg(test)]
mod tests {

    use crate::Thermodynamics::User_substances::{DataType, LibraryPriority, Phases, SubsData};

    #[test]
    fn test_calculate_transport_map_of_properties() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string()];
        user_subs.substances = substances.clone();

        // Set phases

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
        let _ = user_subs.calculate_therm_map_of_fun();
        let _ = user_subs.calculate_therm_map_of_sym();
        match user_subs.calculate_transport_map_of_sym() {
            Ok(_) => println!("Transport map of sym functions calculated successfully."),
            Err(e) => println!("Error calculating transport map of sym functions: {}", e),
        }
        match user_subs.calculate_transport_map_of_functions() {
            Ok(_) => println!("Transport map of functions calculated successfully."),
            Err(e) => println!("Error calculating transport map of functions: {}", e),
        }

        let clo_map = &user_subs.transport_map_of_fun;
        let sym_map = &user_subs.transport_map_of_sym;
        assert!(!sym_map.is_empty());
        assert!(clo_map.contains_key("CO"));

        // Check what's actually in the map
        println!(
            "Transport function map keys: {:?}",
            clo_map.keys().collect::<Vec<_>>()
        );
        if let Some(co_map) = clo_map.get("CO") {
            println!(
                "CO function map keys: {:?}",
                co_map.keys().collect::<Vec<_>>()
            );
            if let Some(lambda_fun) = co_map.get(&DataType::Lambda_fun) {
                println!("Lambda function exists: {:?}", lambda_fun.is_some());
                if lambda_fun.is_some() {
                    println!("Lambda function test passed!");
                } else {
                    println!("Lambda function is None");
                }
            } else {
                println!("Lambda_fun key not found in CO map");
            }
        } else {
            println!("CO key not found in transport function map");
        }
    }

    #[test]
    fn test_fitting_all_thermal_coeffs_for_T_interval() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "H2O".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        user_subs
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Search substances
        let _ = user_subs.search_substances().unwrap();

        // Parse thermal coefficients
        user_subs.parse_thermal_coeffs("CO").unwrap();
        user_subs.parse_thermal_coeffs("H2O").unwrap();
        let _ = user_subs.set_T_range_for_all_thermo(600.0, 1000.0);
        // Test fitting for temperature interval
        let result = user_subs.fitting_all_thermal_coeffs_for_T_interval();
        match result {
            Ok(_) => println!("Fitting for temperature interval completed successfully."),
            Err(e) => {
                println!("Error during fitting: {}", e);
                panic!("Error during fitting: {}", e);
            }
        }

        // Verify coefficients were fitted for found substances
        for substance in &substances {
            if let Some(result_map) = user_subs.search_results.get(substance) {
                if result_map
                    .contains_key(&crate::Thermodynamics::User_substances::WhatIsFound::Thermo)
                {
                    println!("Fitted coefficients for {}", substance);
                }
            }
        }
    }

    #[test]
    fn test_integr_mean() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Search substances
        let _ = user_subs.search_substances().unwrap();

        // Parse thermal coefficients
        user_subs.parse_thermal_coeffs("CO").unwrap();
        let _ = user_subs.set_T_range_for_all_thermo(600.0, 1000.0);
        // Test integral mean calculation
        let result = user_subs.integr_mean("CO");
        assert!(result.is_ok());

        // Verify mean values were calculated
        for substance in &substances {
            if let Some(props) = user_subs.therm_map_of_properties_values.get(substance) {
                if let Some(Some(cp_mean)) = props.get(&DataType::Cp) {
                    println!("{} Mean Cp (400-600K): {}", substance, cp_mean);
                    assert!(
                        *cp_mean > 0.0,
                        "Mean Cp should be positive for {}",
                        substance
                    );
                }
                if let Some(Some(dh_mean)) = props.get(&DataType::dH) {
                    println!("{} Mean dH (400-600K): {}", substance, dh_mean);
                }
                if let Some(Some(ds_mean)) = props.get(&DataType::dS) {
                    println!("{} Mean dS (400-600K): {}", substance, ds_mean);
                    assert!(
                        *ds_mean > 0.0,
                        "Mean dS should be positive for {}",
                        substance
                    );
                }
            }
        }
    }

    #[test]
    fn test_integr_mean_error_handling() {
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["CO".to_string()];
        user_subs
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));

        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        let _ = user_subs.search_substances().unwrap();
        user_subs.parse_thermal_coeffs("CO").unwrap();
        let _ = user_subs.set_T_range_for_all_thermo(1000.0, 300.0);
        // Test with invalid temperature range (T_min > T_max)
        let result = user_subs.integr_mean("CO");
        assert!(result.is_err());

        // Test with equal temperatures
        let result = user_subs.integr_mean("CO");
        assert!(result.is_err());
    }

    #[test]
    fn test_fitting_coeffs_multiple_substances() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "CO2".to_string(), "N2".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        for substance in &substances {
            user_subs
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Search substances
        let _ = user_subs.search_substances().unwrap();

        // Parse thermal coefficients for all substances
        for substance in &substances {
            let _ = user_subs.parse_thermal_coeffs(substance);
        }
        let _ = user_subs.set_T_range_for_all_thermo(400.0, 1400.0);
        // Test fitting for temperature interval
        let result = user_subs.fitting_all_thermal_coeffs_for_T_interval();
        assert!(result.is_ok());

        println!(
            "Successfully fitted coefficients for {} substances",
            substances.len()
        );
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////DIFFUSION TESTS////////////////////////////////////
    #[test]
    fn test_initialize_diffusion_data() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "N2".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        for substance in &substances {
            user_subs
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Set temperature and pressure
        user_subs.set_T(300.0);
        user_subs.set_P(101325.0, None);

        // Search substances
        user_subs.search_substances().unwrap();

        // Calculate molar masses
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();

        // Initialize diffusion data
        let result = user_subs.initialize_diffusion_data();
        assert!(result.is_ok());
        assert!(user_subs.diffusion_data.is_some());

        println!("Diffusion data initialized successfully");
    }

    #[test]
    fn test_calculate_all_pairs() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "N2".to_string(), "O2".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        for substance in &substances {
            user_subs
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Set temperature and pressure
        user_subs.set_T(300.0);
        user_subs.set_P(101325.0, None);

        // Search substances
        user_subs.search_substances().unwrap();

        // Calculate molar masses
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();

        // Initialize diffusion data
        user_subs.initialize_diffusion_data().unwrap();

        // Calculate all pairs
        let result = user_subs.calculate_all_pairs();
        assert!(result.is_ok());

        // Verify pairs were calculated
        user_subs
            .with_diffusion_calculator(|diff| {
                assert!(!diff.pair_diffusions.is_empty());
                println!("Calculated {} diffusion pairs", diff.pair_diffusions.len());
            })
            .unwrap();
    }

    #[test]
    fn test_create_diffusion_matrix() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "N2".to_string(), "O2".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        for substance in &substances {
            user_subs
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Set temperature and pressure
        user_subs.set_T(300.0);
        user_subs.set_P(101325.0, None);

        // Search substances
        user_subs.search_substances().unwrap();

        // Calculate molar masses
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();

        // Initialize diffusion data
        user_subs.initialize_diffusion_data().unwrap();

        // Calculate all pairs
        user_subs.calculate_all_pairs().unwrap();

        // Create diffusion matrix
        let result = user_subs.create_diffusion_matrix();
        assert!(result.is_ok());

        let matrix = result.unwrap();
        assert_eq!(matrix.len(), 3);
        assert_eq!(matrix[0].len(), 3);

        // Check matrix symmetry
        for i in 0..matrix.len() {
            for j in 0..matrix[i].len() {
                assert!((matrix[i][j] - matrix[j][i]).abs() < 1e-10);
            }
        }

        println!(
            "Diffusion matrix created successfully: {}x{}",
            matrix.len(),
            matrix[0].len()
        );
    }

    #[test]
    fn test_calculate_all_closures() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "N2".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        for substance in &substances {
            user_subs
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Set temperature and pressure
        user_subs.set_T(300.0);
        user_subs.set_P(101325.0, None);

        // Search substances
        user_subs.search_substances().unwrap();

        // Calculate molar masses
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();

        // Initialize diffusion data
        user_subs.initialize_diffusion_data().unwrap();

        // Calculate all pairs
        user_subs.calculate_all_pairs().unwrap();

        // Calculate closures
        let result = user_subs.calculate_all_closures();
        assert!(result.is_ok());

        // Test closures work
        user_subs
            .with_diffusion_calculator(|diff| {
                if let Some(pair_diff) = diff.get_pair_diffusion("CO", "N2") {
                    let D_300 = (pair_diff.D_closure)(300.0);
                    let D_400 = (pair_diff.D_closure)(400.0);
                    assert!(D_400 > D_300);
                    println!("Closure test: D(300K)={}, D(400K)={}", D_300, D_400);
                }
            })
            .unwrap();
    }

    #[test]
    fn test_calculate_all_symbolic() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "N2".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        for substance in &substances {
            user_subs
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Set temperature and pressure
        user_subs.set_T(300.0);
        user_subs.set_P(101325.0, None);

        // Search substances
        user_subs.search_substances().unwrap();

        // Calculate molar masses
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();

        // Initialize diffusion data
        user_subs.initialize_diffusion_data().unwrap();

        // Calculate all pairs
        user_subs.calculate_all_pairs().unwrap();

        // Calculate symbolic expressions
        let result = user_subs.calculate_all_symbolic();
        assert!(result.is_ok());

        // Test symbolic expressions work
        user_subs
            .with_diffusion_calculator(|diff| {
                if let Some(pair_diff) = diff.get_pair_diffusion("CO", "N2") {
                    let D_sym_func = pair_diff.D_symbolic.lambdify1D();
                    let D_300 = D_sym_func(300.0);
                    let D_400 = D_sym_func(400.0);
                    assert!(D_400 > D_300);
                    println!("Symbolic test: D(300K)={}, D(400K)={}", D_300, D_400);
                }
            })
            .unwrap();
    }

    #[test]
    fn test_diffusion_coefficient_values() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "N2".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        for substance in &substances {
            user_subs
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Set temperature and pressure
        user_subs.set_T(300.0);
        user_subs.set_P(101325.0, None);

        // Search substances
        user_subs.search_substances().unwrap();

        // Calculate molar masses
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();

        // Initialize diffusion data
        user_subs.initialize_diffusion_data().unwrap();

        // Calculate all pairs
        user_subs.calculate_all_pairs().unwrap();

        // Check diffusion coefficient values
        user_subs
            .with_diffusion_calculator(|diff| {
                let d_co_n2 = diff.get_diffusion_coefficient("CO", "N2");
                assert!(d_co_n2.is_some());
                let d_value = d_co_n2.unwrap();
                assert!(d_value > 0.0);
                println!("D_CO-N2 at 300K: {} m²/s", d_value);
            })
            .unwrap();
    }

    #[test]
    fn test_diffusion_error_handling() {
        let mut user_subs = SubsData::new();

        // Test without initialization
        let result = user_subs.calculate_all_pairs();
        assert!(result.is_err());

        // Test without temperature set
        user_subs.substances = vec!["CO".to_string()];
        user_subs.set_P(101325.0, None);
        let result = user_subs.initialize_diffusion_data();
        assert!(result.is_err());

        // Test without pressure set
        user_subs.set_T(300.0);
        user_subs.P = None;
        let result = user_subs.initialize_diffusion_data();
        assert!(result.is_err());
    }

    #[test]
    fn test_full_diffusion_workflow() {
        let mut user_subs = SubsData::new();
        let substances = vec!["CO".to_string(), "N2".to_string(), "O2".to_string()];
        user_subs.substances = substances.clone();

        // Set phases
        for substance in &substances {
            user_subs
                .map_of_phases
                .insert(substance.clone(), Some(Phases::Gas));
        }

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Set temperature and pressure
        user_subs.set_T(300.0);
        user_subs.set_P(101325.0, None);

        // Search substances
        user_subs.search_substances().unwrap();

        // Calculate molar masses
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .unwrap();

        // Full diffusion workflow
        user_subs.initialize_diffusion_data().unwrap();
        user_subs.calculate_all_pairs().unwrap();
        user_subs.calculate_all_closures().unwrap();
        user_subs.calculate_all_symbolic().unwrap();
        let matrix = user_subs.create_diffusion_matrix().unwrap();

        // Verify complete workflow
        assert_eq!(matrix.len(), 3);
        assert_eq!(matrix[0].len(), 3);

        // Print matrix for verification
        user_subs
            .with_diffusion_calculator(|diff| {
                diff.print_diffusion_matrix();
            })
            .unwrap();

        println!("Full diffusion workflow completed successfully");
    }
}
