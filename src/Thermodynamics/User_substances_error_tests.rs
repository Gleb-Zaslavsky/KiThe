#[cfg(test)]
mod error_handling_tests {
    use crate::Thermodynamics::User_substances::{LibraryPriority, SubsData};
    use crate::Thermodynamics::User_substances_error::SubsDataError;

    #[test]
    fn test_invalid_temperature_errors() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Test negative temperature
        let result = subs_data.extract_thermal_coeffs("CO", -100.0);
        assert!(matches!(
            result,
            Err(SubsDataError::InvalidTemperature(-100.0))
        ));

        // Test zero temperature
        let result = subs_data.extract_thermal_coeffs("CO", 0.0);
        assert!(matches!(
            result,
            Err(SubsDataError::InvalidTemperature(0.0))
        ));

        // Test invalid temperature in calculate_thermo_properties
        let result = subs_data.calculate_thermo_properties("CO", -273.15);
        assert!(matches!(
            result,
            Err(SubsDataError::InvalidTemperature(-273.15))
        ));

        // Test invalid temperature in transport calculations
        let result = subs_data.extract_transport_coeffs("CO", -50.0);
        assert!(matches!(
            result,
            Err(SubsDataError::InvalidTemperature(-50.0))
        ));
    }

    #[test]
    fn test_substance_not_found_errors() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["NonExistentSubstance123".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Test extract_thermal_coeffs with non-existent substance
        let result = subs_data.extract_thermal_coeffs("NonExistentSubstance123", 400.0);
        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));

        // Test calculate_thermo_properties with non-existent substance
        let result = subs_data.calculate_thermo_properties("NonExistentSubstance123", 400.0);
        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));

        // Test with completely unknown substance
        let result = subs_data.extract_thermal_coeffs("CompletelyUnknownSubstance", 400.0);
        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));
    }

    #[test]
    fn test_calculator_not_available_errors() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];

        // Set only transport library, no thermo library
        subs_data.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Try to extract thermal coefficients when only transport calculator is available
        let result = subs_data.extract_thermal_coeffs("CO", 400.0);
        assert!(matches!(
            result,
            Err(SubsDataError::CalculatorNotAvailable { .. })
        ));

        // Try to calculate thermo properties when only transport calculator is available
        let result = subs_data.calculate_thermo_properties("CO", 400.0);
        assert!(matches!(
            result,
            Err(SubsDataError::CalculatorNotAvailable { .. })
        ));
    }

    #[test]
    fn test_pressure_not_set_error() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Don't set pressure - this should cause an error
        let result = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);
        assert!(matches!(result, Err(SubsDataError::PressureNotSet)));
    }

    #[test]
    fn test_molar_mass_not_found_error() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();
        subs_data.set_P(101325.0, None);

        // Don't set molar mass - this should cause an error
        let result = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);
        assert!(matches!(result, Err(SubsDataError::MolarMassNotFound(_))));
    }

    //#[test]
    fn test_heat_capacity_not_available_error() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["meow".to_string(), "H2O".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string(), "Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Don't calculate thermal properties first
        let result = subs_data.calculate_transport_map_of_properties(400.0);
        println!("{:?}", result);
        assert!(matches!(
            result,
            Err(SubsDataError::HeatCapacityNotAvailable(_))
        ));
    }

    #[test]
    fn test_coefficient_extraction_failed_error() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Try to extract coefficients at extremely high temperature that might cause issues
        let result = subs_data.extract_thermal_coeffs("CO", 50000.0);
        // This might succeed or fail depending on the data range, but we test the error type
        if let Err(error) = result {
            assert!(matches!(
                error,
                SubsDataError::CoefficientExtractionFailed { .. } | SubsDataError::ThermoError(_)
            ));
        }
    }

    #[test]
    fn test_error_propagation_in_batch_operations() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec![
            "CO".to_string(),
            "NonExistentSubstance".to_string(),
            "H2O".to_string(),
        ];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Test extract_all_thermal_coeffs with invalid temperature
        let result = subs_data.extract_all_thermal_coeffs(-100.0);
        assert!(matches!(
            result,
            Err(SubsDataError::InvalidTemperature(-100.0))
        ));

        // Test extract_all_thermal_coeffs with non-existent substance
        let result = subs_data.extract_all_thermal_coeffs(400.0);
        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));
    }

    #[test]
    fn test_function_creation_failed_scenarios() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();
        subs_data.extract_thermal_coeffs("CO", 400.0).unwrap();

        // This should normally succeed, but we test the error handling structure
        let result = subs_data.calculate_therm_map_of_fun();
        // If it fails, it should be a FunctionCreationFailed error
        if let Err(error) = result {
            assert!(matches!(
                error,
                SubsDataError::FunctionCreationFailed { .. }
            ));
        }
    }

    #[test]
    fn test_symbolic_creation_failed_scenarios() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();
        subs_data.extract_thermal_coeffs("CO", 400.0).unwrap();

        // This should normally succeed, but we test the error handling structure
        let result = subs_data.calculate_therm_map_of_sym();
        // If it fails, it should be a SymbolicCreationFailed error
        if let Err(error) = result {
            assert!(matches!(
                error,
                SubsDataError::SymbolicCreationFailed { .. }
            ));
        }
    }

    #[test]
    fn test_error_display_messages() {
        // Test that error messages are properly formatted
        let error = SubsDataError::SubstanceNotFound("TestSubstance".to_string());
        let message = format!("{}", error);
        assert!(message.contains("TestSubstance"));
        assert!(message.contains("not found"));

        let error = SubsDataError::InvalidTemperature(-100.0);
        let message = format!("{}", error);
        assert!(message.contains("-100"));
        assert!(message.contains("Invalid temperature"));

        let error = SubsDataError::CalculatorNotAvailable {
            substance: "CO".to_string(),
            calc_type: "Thermo".to_string(),
        };
        let message = format!("{}", error);
        assert!(message.contains("CO"));
        assert!(message.contains("Thermo"));
        assert!(message.contains("not available"));

        let error = SubsDataError::PressureNotSet;
        let message = format!("{}", error);
        assert!(message.contains("Pressure"));
        assert!(message.contains("must be set"));

        let error = SubsDataError::MolarMassNotFound("H2O".to_string());
        let message = format!("{}", error);
        assert!(message.contains("H2O"));
        assert!(message.contains("Molar mass"));

        let error = SubsDataError::HeatCapacityNotAvailable("CO2".to_string());
        let message = format!("{}", error);
        assert!(message.contains("CO2"));
        assert!(message.contains("Heat capacity"));
    }

    #[test]
    fn test_error_chain_with_question_mark_operator() {
        fn test_function(subs_data: &mut SubsData) -> Result<(), SubsDataError> {
            subs_data.extract_thermal_coeffs("NonExistent", 400.0)?;
            subs_data.calculate_thermo_properties("NonExistent", 400.0)?;
            Ok(())
        }

        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["NonExistent".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        let result = test_function(&mut subs_data);
        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));
    }

    #[test]
    fn test_transport_error_scenarios() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()], // Only thermo, no transport
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Try transport operations without transport calculator
        let result = subs_data.extract_transport_coeffs("CO", 400.0);
        assert!(matches!(
            result,
            Err(SubsDataError::CalculatorNotAvailable { .. })
        ));

        let result = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);
        assert!(matches!(
            result,
            Err(SubsDataError::CalculatorNotAvailable { .. })
        ));
    }

    #[test]
    fn test_extreme_temperature_values() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Test extremely low temperatures
        let temperatures = vec![-1000.0, -273.16, -1.0, 0.0];
        for temp in temperatures {
            let result = subs_data.extract_thermal_coeffs("CO", temp);
            assert!(matches!(result, Err(SubsDataError::InvalidTemperature(_))));
        }

        // Test extremely high temperatures (might cause coefficient extraction to fail)
        let result = subs_data.extract_thermal_coeffs("CO", 100000.0);
        if let Err(error) = result {
            // Should be either InvalidTemperature or CoefficientExtractionFailed
            assert!(matches!(
                error,
                SubsDataError::InvalidTemperature(_)
                    | SubsDataError::CoefficientExtractionFailed { .. }
                    | SubsDataError::ThermoError(_)
            ));
        }
    }

    #[test]
    fn test_error_recovery_scenarios() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string(), "H2O".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // First, cause an error with invalid temperature
        let result = subs_data.extract_thermal_coeffs("CO", -100.0);
        assert!(result.is_err());

        // Then, recover with valid temperature
        let result = subs_data.extract_thermal_coeffs("CO", 400.0);
        assert!(result.is_ok());

        // Verify that the system still works after error recovery
        let result = subs_data.calculate_thermo_properties("CO", 400.0);
        assert!(result.is_ok());
    }

    #[test]
    fn test_mixed_success_and_failure_scenarios() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec![
            "CO".to_string(),
            "H2O".to_string(),
            "NonExistentSubstance".to_string(),
        ];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        // Extract coefficients for existing substances should work
        assert!(subs_data.extract_thermal_coeffs("CO", 400.0).is_ok());
        assert!(subs_data.extract_thermal_coeffs("H2O", 400.0).is_ok());

        // But should fail for non-existent substance
        let result = subs_data.extract_thermal_coeffs("NonExistentSubstance", 400.0);
        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));

        // Calculate properties should work for existing substances
        assert!(subs_data.calculate_thermo_properties("CO", 400.0).is_ok());
        assert!(subs_data.calculate_thermo_properties("H2O", 400.0).is_ok());

        // But fail for non-existent substance
        let result = subs_data.calculate_thermo_properties("NonExistentSubstance", 400.0);
        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));
    }
}
