//! Example demonstrating the new error handling system for SubsData operations
//!
//! This example shows how to properly handle errors when working with substance data,
//! including error propagation from thermo_api and transport_api modules.

use KiThe::Thermodynamics::User_substances::{LibraryPriority, SubsData};
use KiThe::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};

fn main() -> SubsDataResult<()> {
    // Create a new SubsData instance
    let mut subs_data = SubsData::new();
    subs_data.set_substances(vec![
        "CO".to_string(),
        "H2O".to_string(),
        "InvalidSubstance".to_string(),
    ]);

    // Set library priorities
    subs_data
        .set_multiple_library_priorities(vec!["NASA_gas".to_string()], LibraryPriority::Priority);

    // Perform the search
    subs_data.search_substances()?;

    // Example 1: Proper error handling for coefficient extraction
    match subs_data.extract_thermal_coeffs("CO", 400.0) {
        Ok(()) => println!("Successfully extracted thermal coefficients for CO"),
        Err(SubsDataError::SubstanceNotFound(substance)) => {
            println!("Substance not found: {}", substance);
        }
        Err(SubsDataError::InvalidTemperature(temp)) => {
            println!("Invalid temperature: {} K", temp);
        }
        Err(e) => println!("Other error: {}", e),
    }

    // Example 2: Error propagation with the ? operator
    subs_data.extract_thermal_coeffs("H2O", 400.0)?;
    println!("Successfully extracted thermal coefficients for H2O");

    // Example 3: Handling calculation errors
    match subs_data.calculate_thermo_properties("CO", 400.0) {
        Ok((cp, dh, ds)) => {
            println!("CO properties at 400K: Cp={}, dH={}, dS={}", cp, dh, ds);
        }
        Err(SubsDataError::CalculatorNotAvailable {
            substance,
            calc_type,
        }) => {
            println!(
                "Calculator '{}' not available for substance '{}'",
                calc_type, substance
            );
        }
        Err(e) => println!("Calculation failed: {}", e),
    }

    // Example 4: Batch operations with proper error handling
    match subs_data.calculate_therm_map_of_properties(400.0) {
        Ok(()) => println!("Successfully calculated thermodynamic properties for all substances"),
        Err(e) => println!("Failed to calculate properties: {}", e),
    }

    // Example 5: Handling invalid temperature
    match subs_data.extract_thermal_coeffs("CO", -100.0) {
        Err(SubsDataError::InvalidTemperature(temp)) => {
            println!("Correctly caught invalid temperature: {} K", temp);
        }
        _ => println!("Should have caught invalid temperature!"),
    }

    // Example 6: Handling missing substance
    match subs_data.calculate_thermo_properties("NonExistentSubstance", 400.0) {
        Err(SubsDataError::SubstanceNotFound(substance)) => {
            println!("Correctly caught missing substance: {}", substance);
        }
        _ => println!("Should have caught missing substance!"),
    }

    println!("Error handling example completed successfully!");
    Ok(())
}

/// Helper function demonstrating error propagation in user code
#[allow(dead_code)]
fn calculate_properties_for_substances(
    subs_data: &mut SubsData,
    substances: &[String],
    temperature: f64,
) -> SubsDataResult<Vec<(f64, f64, f64)>> {
    let mut results = Vec::new();

    for substance in substances {
        // Extract coefficients first
        subs_data.extract_thermal_coeffs(substance, temperature)?;

        // Calculate properties
        let properties = subs_data.calculate_thermo_properties(substance, temperature)?;
        results.push(properties);
    }

    Ok(results)
}

/// Example of custom error handling for specific use cases
#[allow(dead_code)]
fn safe_property_calculation(
    subs_data: &mut SubsData,
    substance: &str,
    temperature: f64,
) -> Option<(f64, f64, f64)> {
    match subs_data.calculate_thermo_properties(substance, temperature) {
        Ok(properties) => Some(properties),
        Err(e) => {
            eprintln!(
                "Warning: Failed to calculate properties for {}: {}",
                substance, e
            );
            None
        }
    }
}
