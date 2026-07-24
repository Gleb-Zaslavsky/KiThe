#[cfg(test)]
mod error_handling_tests {
    use crate::Kinetics::error::KineticsError;
    use crate::Thermodynamics::DBhandlers::NIST_parser::NistError as NistParserError;
    use crate::Thermodynamics::DBhandlers::NISTdata::NISTError;
    use crate::Thermodynamics::DBhandlers::thermo_api::ThermoError;
    use crate::Thermodynamics::DBhandlers::transport_api::TransportError;
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
    use crate::Thermodynamics::User_substances_error::{LogErrorType, SubsDataError};
    use std::error::Error;

    #[test]
    fn test_unknown_calculator_library_is_reported() {
        let subs_data = SubsData::new();

        // Unsupported libraries should fail through the typed API instead of panicking.
        let result = subs_data.create_calculator("Definitely_Not_A_Library");
        assert!(matches!(
            result,
            Err(SubsDataError::UnsupportedLibrary(name)) if name == "Definitely_Not_A_Library"
        ));
    }

    #[test]
    fn test_wrapped_thermo_error_preserves_source_chain() {
        let serde_error = serde_json::from_str::<serde_json::Value>("not valid json").unwrap_err();
        let error = SubsDataError::from(ThermoError::SerdeError(serde_error));

        let first_source = error
            .source()
            .expect("thermo wrapper should expose inner source");
        assert!(first_source.source().is_some());
    }

    #[test]
    fn test_wrapped_transport_error_preserves_source_chain() {
        let serde_error = serde_json::from_str::<serde_json::Value>("not valid json").unwrap_err();
        let error = SubsDataError::from(TransportError::SerdeError(serde_error));

        let first_source = error
            .source()
            .expect("transport wrapper should expose inner source");
        assert!(first_source.source().is_some());
    }

    #[test]
    fn test_nist_parser_error_preserves_source_chain_through_wrappers() {
        let parse_error = url::Url::parse("://not-a-url").unwrap_err();
        let nist_error = NISTError::ParserError(NistParserError::UrlError(parse_error));
        let thermo_error = ThermoError::from(nist_error);
        let error = SubsDataError::from(thermo_error);

        let thermo_source = error
            .source()
            .expect("subsdata should expose thermo source");
        let parser_source = thermo_source
            .source()
            .expect("thermo wrapper should expose parser source");
        assert!(
            parser_source.to_string().contains("URL parsing error")
                || parser_source.to_string().contains("url")
        );
    }

    #[test]
    fn test_calculation_failed_constructor_preserves_source() {
        let source = KineticsError::InvalidReactionData("bad chemistry".to_string());
        let error = SubsDataError::calculation_failed_with_source(
            "CO2",
            "molar mass calculation",
            source.to_string(),
            source,
        );

        let wrapped_source = error
            .source()
            .expect("calculation errors should retain their nested source");
        assert!(wrapped_source.to_string().contains("bad chemistry"));
    }

    #[test]
    fn test_nist_retrieval_failed_constructor_preserves_source() {
        let serde_error = serde_json::from_str::<serde_json::Value>("not valid json").unwrap_err();
        let nist_error = NISTError::SerdeError(serde_error);
        let error = SubsDataError::nist_retrieval_failed_with_source(
            "CO2",
            nist_error.to_string(),
            nist_error,
        );

        let wrapped_source = error
            .source()
            .expect("nist retrieval errors should retain their nested source");
        assert!(wrapped_source.to_string().contains("NIST"));
    }

    #[test]
    fn test_element_composition_failure_preserves_calculation_source() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["DefinitelyNotAFormula".to_string()];

        let result =
            SubsData::calculate_elem_composition_and_molar_mass_local(&mut subs_data, None);
        let error = result.expect_err("invalid formula should fail the composition pass");

        let source = error
            .source()
            .expect("calculation failures should preserve the nested source");
        assert!(source.to_string().contains("invalid"));
        assert!(subs_data.elem_composition_matrix.is_none());
        assert!(subs_data.molar_mass_by_substance.is_empty());
        assert!(subs_data.unique_elements.is_empty());
    }

    #[test]
    fn test_nist_fallback_failure_preserves_source() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["bad name".to_string()];
        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        subs_data.insert_not_found_if_absent("bad name");

        let result = subs_data.if_not_found_go_NIST();
        let error = result.expect_err("the invalid NIST query should fail");

        let source = error
            .source()
            .expect("NIST fallback failures should preserve the nested source");
        assert!(
            source.to_string().contains("Failed to retrieve NIST data"),
            "unexpected NIST fallback source: {}",
            source
        );
        assert!(
            source
                .source()
                .expect("the wrapped NIST error should remain accessible")
                .to_string()
                .contains("Substance not found"),
            "unexpected nested NIST source: {}",
            source
        );
    }

    #[test]
    fn test_error_messages_are_domain_specific() {
        let thermo_message = ThermoError::SerdeError(
            serde_json::from_str::<serde_json::Value>("not valid json").unwrap_err(),
        )
        .to_string();
        let nist_message = NISTError::SerdeError(
            serde_json::from_str::<serde_json::Value>("not valid json").unwrap_err(),
        )
        .to_string();

        assert!(thermo_message.contains("thermodynamic data"));
        assert!(nist_message.contains("NIST data"));
        assert!(!thermo_message.contains("NASA data"));
    }

    #[test]
    fn test_parser_failures_remain_nist_retrieval_failures_in_logs() {
        let parse_error = url::Url::parse("://not-a-url").unwrap_err();
        let nist_error = NISTError::ParserError(NistParserError::UrlError(parse_error));
        let log_error: LogErrorType = ThermoError::from(nist_error).into();

        assert!(
            matches!(log_error, LogErrorType::NistRetrievalFailed(msg) if msg.contains("NIST data"))
        );
    }

    #[test]
    fn test_subsdata_error_uses_shared_nist_retrieval_mapping() {
        let parse_error = url::Url::parse("://not-a-url").unwrap_err();
        let nist_error = NISTError::ParserError(NistParserError::UrlError(parse_error));
        let subs_error = SubsDataError::from(ThermoError::from(nist_error));
        let log_error: LogErrorType = (&subs_error).into();

        assert!(
            matches!(log_error, LogErrorType::NistRetrievalFailed(msg) if msg.contains("NIST data"))
        );
    }

    #[test]
    fn test_setters_reject_non_physical_values() {
        let mut subs_data = SubsData::new();

        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0, -1.0] {
            assert!(
                matches!(
                    subs_data.set_T(value),
                    Err(SubsDataError::InvalidPhysicalValue { field, value: got })
                        if field == "temperature" && got.to_bits() == value.to_bits()
                ),
                "temperature {:?}",
                value
            );
            assert!(
                matches!(
                    subs_data.set_P(value, None),
                    Err(SubsDataError::InvalidPhysicalValue { field, value: got })
                        if field == "pressure" && got.to_bits() == value.to_bits()
                ),
                "pressure {:?}",
                value
            );

            let masses = std::collections::HashMap::from([("CO".to_string(), value)]);
            assert!(
                matches!(
                    subs_data.set_M(masses, None),
                    Err(SubsDataError::InvalidPhysicalValue { field, value: got })
                        if field == "molar_mass:CO" && got.to_bits() == value.to_bits()
                ),
                "molar mass {:?}",
                value
            );
        }
    }

    #[test]
    fn test_invalid_temperature_errors() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        let invalid_temperatures = [-100.0, 0.0, -273.15, -50.0];
        for temperature in invalid_temperatures {
            let thermal = subs_data.extract_thermal_coeffs("CO", temperature);
            let calculated = subs_data.calculate_thermo_properties("CO", temperature);
            let _transport = subs_data.extract_transport_coeffs("CO", temperature);

            assert!(matches!(
                thermal,
                Err(SubsDataError::InvalidTemperature(value)) if value == temperature
            ));
            assert!(matches!(
                calculated,
                Err(SubsDataError::InvalidTemperature(value)) if value == temperature
            ));
        }
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
        let _ = subs_data.set_P(101325.0, None);

        // Don't set molar mass - this should cause an error
        let result = subs_data.calculate_transport_properties("CO", 400.0, Some(30.0), None);
        assert!(matches!(result, Err(SubsDataError::MolarMassNotFound(_))));
    }

    //#[test]
    #[allow(dead_code)]
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
    fn test_thermo_batch_map_does_not_publish_partial_state_on_failure() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();
        subs_data.extract_thermal_coeffs("CO", 400.0).unwrap();

        subs_data.calculate_therm_map_of_properties(400.0).unwrap();
        let previous_map = subs_data.therm_map_of_properties_values.clone();

        subs_data.substances = vec!["CO".to_string(), "NonExistentSubstance".to_string()];
        let result = subs_data.calculate_therm_map_of_properties(400.0);

        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));
        assert_eq!(subs_data.therm_map_of_properties_values, previous_map);
    }

    #[test]
    fn test_thermo_batch_map_all_fail_returns_error() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["NonExistentSubstance".to_string()];
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();

        let result = subs_data.calculate_therm_map_of_properties(400.0);

        assert!(matches!(result, Err(SubsDataError::SubstanceNotFound(_))));
        assert!(subs_data.therm_map_of_properties_values.is_empty());
    }

    #[test]
    fn test_transport_batch_map_does_not_publish_partial_state_on_failure() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["CO".to_string()];
        subs_data
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        subs_data.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        subs_data.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );
        let _ = subs_data.search_substances();
        let _ = subs_data.set_P(1.0e5, None);
        let _ = subs_data.set_M(
            std::collections::HashMap::from([("CO".to_string(), 28.0)]),
            None,
        );
        subs_data.extract_thermal_coeffs("CO", 400.0).unwrap();
        subs_data.calculate_therm_map_of_properties(400.0).unwrap();
        subs_data.extract_transport_coeffs("CO", 400.0).unwrap();

        subs_data
            .calculate_transport_map_of_properties(400.0)
            .unwrap();
        let previous_map = subs_data.transport_map_of_properties_values.clone();

        subs_data.substances = vec!["CO".to_string(), "NonExistentSubstance".to_string()];
        let result = subs_data.calculate_transport_map_of_properties(400.0);

        // The batch must fail without publishing a partial cache, but the
        // precise error path depends on which prerequisite is missing first.
        assert!(result.is_err());
        assert_eq!(subs_data.transport_map_of_properties_values, previous_map);
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
