use crate::Thermodynamics::DBhandlers::thermo_api::{
    ThermoCalculator, ThermoEnum, create_thermal_by_name,
};
use crate::Thermodynamics::DBhandlers::transport_api::{
    TransportCalculator, TransportEnum, create_transport_calculator_by_name,
};
use crate::Thermodynamics::thermo_lib_api::ThermoData;
//use RustedSciThe::symbolic::symbolic_engine::Expr;
use serde_json::Value;
use std::collections::HashMap;
/*
Core Types:
        SearchResult enum: Represents the result of searching for a substance
        LibraryPriority enum: Defines priority levels (Priority or Permitted)
        UserSubstances struct: Main structure managing the substances and search process
Main Functionality:
        new(): Creates a new instance with a list of substances to search for
        set_library_priority(): Sets priority level for a single library
        set_multiple_library_priorities(): Sets priority level for multiple libraries at once
        search_substances(): Performs the search according to priority rules
Query Methods:
        get_substance_result(): Gets result for a specific substance
        get_all_results(): Gets all search results
        get_not_found_substances(): Lists substances not found in any library
        get_priority_found_substances(): Lists substances found in priority libraries
        get_permitted_found_substances(): Lists substances found in permitted libraries
        print_search_summary(): Prints a human-readable summary of search results
Search Logic:
        First searches in priority libraries
        If not found, searches in permitted libraries
        Stores the results including which library the substance was found in
        Maintains priority information in the results
*/
//User Substance module that implements the following functionality 1) A vector of substances is given.
//  Using the thermo_lib_api module, get data from libraries for substances from the vector. Moreover, the
//user somehow divides all libraries into two types "Priority" and "permitted" The substance is first searched
// for in the library assigned to the "priority" group, if it is not there -
//in the "permitted" group, if it is not there we move on to the next substance. Information about which
// library which substance was found was found in and whether it was found at all must be somehow saved
// and available to the user

/// Represents the type of calculator for a substance
#[derive(Debug)]
pub enum CalculatorType {
    Thermo(ThermoEnum),
    Transport(TransportEnum),
}

///SearchResult enum: Represents the result of searching for a substance
#[derive(Debug)]
pub enum SearchResult {
    Found {
        library: String,
        priority_type: LibraryPriority,
        data: Value,
        calculator: Option<CalculatorType>,
    },
    NotFound,
}
///LibraryPriority enum: Defines priority levels (Priority or Permitted)
#[derive(Debug, Clone, PartialEq)]
pub enum LibraryPriority {
    Priority,
    Permitted,
}
///UserSubstances struct: Main structure managing the substances and search process
pub struct UserSubstances {
    /// Vector of substances to search for
    substances: Vec<String>,
    /// Maps libraries to their priority level
    library_priorities: HashMap<String, LibraryPriority>,
    /// Stores search results for each substance
    search_results: HashMap<String, SearchResult>,
    /// Reference to the thermo data library
    thermo_data: ThermoData,
}

impl UserSubstances {
    pub fn new(substances: Vec<String>) -> Self {
        Self {
            substances,
            library_priorities: HashMap::new(),
            search_results: HashMap::new(),
            thermo_data: ThermoData::new(),
        }
    }

    /// Set a library's priority level
    pub fn set_library_priority(&mut self, library: String, priority: LibraryPriority) {
        self.library_priorities.insert(library, priority);
    }

    /// Set multiple libraries' priority levels at once
    pub fn set_multiple_library_priorities(
        &mut self,
        libraries: Vec<String>,
        priority: LibraryPriority,
    ) {
        for library in libraries {
            self.library_priorities.insert(library, priority.clone());
        }
    }

    /// Initialize appropriate calculator for a substance based on its library
    fn create_calculator(&self, library: &str) -> Option<CalculatorType> {
        match ThermoData::what_handler_to_use(library).as_str() {
            "NASA" | "NIST" => Some(CalculatorType::Thermo(create_thermal_by_name(library))),
            "transport" | "CEA" => Some(CalculatorType::Transport(
                create_transport_calculator_by_name(library),
            )),
            _ => None,
        }
    }

    /// Search for all substances in the libraries according to priority
    pub fn search_substances(&mut self) {
        // Get all priority libraries
        let priority_libs: Vec<String> = self
            .library_priorities
            .iter()
            .filter(move |&(_, &ref p)| *p == LibraryPriority::Priority)
            .map(|(lib, _)| lib.clone())
            .collect();

        // Get all permitted libraries
        let permitted_libs: Vec<String> = self
            .library_priorities
            .iter()
            .filter(move |&(_, &ref p)| *p == LibraryPriority::Permitted)
            .map(|(lib, _)| lib.clone())
            .collect();

        // Search for each substance
        for substance in &self.substances {
            let mut found = false;

            // Try priority libraries first
            for lib in &priority_libs {
                if let Some(lib_data) = self.thermo_data.LibThermoData.get(lib) {
                    if let Some(substance_data) = lib_data.get(substance) {
                        let calculator = self.create_calculator(lib);

                        // Initialize the calculator if one was created
                        if let Some(calc_type) = &calculator {
                            match calc_type {
                                CalculatorType::Thermo(thermo) => {
                                    let mut thermo = thermo.clone();
                                    let _ = thermo.from_serde(substance_data.clone());
                                }
                                CalculatorType::Transport(transport) => {
                                    let mut transport = transport.clone();
                                    let _ = transport.from_serde(substance_data.clone());
                                }
                            }
                        }

                        self.search_results.insert(
                            substance.clone(),
                            SearchResult::Found {
                                library: lib.clone(),
                                priority_type: LibraryPriority::Priority,
                                data: substance_data.clone(),
                                calculator,
                            },
                        );
                        found = true;
                        break;
                    }
                }
            }

            // If not found in priority libraries, try permitted libraries
            if !found {
                for lib in &permitted_libs {
                    if let Some(lib_data) = self.thermo_data.LibThermoData.get(lib) {
                        if let Some(substance_data) = lib_data.get(substance) {
                            let calculator = self.create_calculator(lib);

                            // Initialize the calculator if one was created
                            if let Some(calc_type) = &calculator {
                                match calc_type {
                                    CalculatorType::Thermo(thermo) => {
                                        let mut thermo = thermo.clone();
                                        let _ = thermo.from_serde(substance_data.clone());
                                    }
                                    CalculatorType::Transport(transport) => {
                                        let mut transport = transport.clone();
                                        let _ = transport.from_serde(substance_data.clone());
                                    }
                                }
                            }

                            self.search_results.insert(
                                substance.clone(),
                                SearchResult::Found {
                                    library: lib.clone(),
                                    priority_type: LibraryPriority::Permitted,
                                    data: substance_data.clone(),
                                    calculator,
                                },
                            );
                            found = true;
                            break;
                        }
                    }
                }
            }

            // If still not found, mark as NotFound
            if !found {
                self.search_results
                    .insert(substance.clone(), SearchResult::NotFound);
            }
        }
    }

    /// Get the search result for a specific substance
    pub fn get_substance_result(&self, substance: &str) -> Option<&SearchResult> {
        self.search_results.get(substance)
    }

    /// Get all search results
    pub fn get_all_results(&self) -> &HashMap<String, SearchResult> {
        &self.search_results
    }

    /// Get a list of substances that were not found
    pub fn get_not_found_substances(&self) -> Vec<String> {
        self.search_results
            .iter()
            .filter(|(_, result)| matches!(result, SearchResult::NotFound))
            .map(|(substance, _)| substance.clone())
            .collect()
    }

    /// Get a list of substances found in priority libraries
    pub fn get_priority_found_substances(&self) -> Vec<String> {
        self.search_results
            .iter()
            .filter(|(_, result)| matches!(result, SearchResult::Found { priority_type, .. } if *priority_type == LibraryPriority::Priority))
            .map(|(substance, _)| substance.clone())
            .collect()
    }

    /// Get a list of substances found in permitted libraries
    pub fn get_permitted_found_substances(&self) -> Vec<String> {
        self.search_results
            .iter()
            .filter(|(_, result)| matches!(result, SearchResult::Found { priority_type, .. } if *priority_type == LibraryPriority::Permitted))
            .map(|(substance, _)| substance.clone())
            .collect()
    }

    /// Calculate thermodynamic properties for a substance
    pub fn calculate_thermo_properties(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> Result<(f64, f64, f64), String> {
        match self.search_results.get_mut(substance) {
            Some(SearchResult::Found {
                calculator: Some(CalculatorType::Thermo(thermo)),
                ..
            }) => {
                let mut thermo = thermo.clone();

                // Extract coefficients and calculate properties
                if let Err(e) = thermo.extract_model_coefficients(temperature) {
                    return Err(format!("Failed to extract coefficients: {}", e));
                }

                if let Err(e) = thermo.calculate_Cp_dH_dS(temperature) {
                    return Err(format!("Failed to calculate properties: {}", e));
                }

                Ok((
                    thermo.get_Cp().unwrap_or(0.0),
                    thermo.get_dh().unwrap_or(0.0),
                    thermo.get_ds().unwrap_or(0.0),
                ))
            }
            Some(SearchResult::Found {
                calculator: Some(CalculatorType::Transport(_)),
                ..
            }) => Err("Substance found but has transport calculator instead of thermo".to_string()),
            Some(SearchResult::Found {
                calculator: None, ..
            }) => Err("Substance found but has no calculator".to_string()),
            Some(SearchResult::NotFound) => Err("Substance not found".to_string()),
            None => Err("Substance not in search results".to_string()),
        }
    }

    /// Calculate transport properties for a substance
    pub fn calculate_transport_properties(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> Result<(f64, f64), String> {
        match self.search_results.get_mut(substance) {
            Some(SearchResult::Found {
                calculator: Some(CalculatorType::Transport(transport)),
                ..
            }) => {
                let mut transport = transport.clone();

                // Extract coefficients and calculate properties
                if let Err(e) = transport.extract_coefficients(temperature) {
                    return Err(format!("Failed to extract coefficients: {}", e));
                }

                let lambda = transport
                    .calculate_lambda(None, None, temperature)
                    .map_err(|e| format!("Failed to calculate lambda: {}", e))?;

                let viscosity = transport
                    .calculate_viscosity(temperature)
                    .map_err(|e| format!("Failed to calculate viscosity: {}", e))?;

                Ok((lambda, viscosity))
            }
            Some(SearchResult::Found {
                calculator: Some(CalculatorType::Thermo(_)),
                ..
            }) => Err("Substance found but has thermo calculator instead of transport".to_string()),
            Some(SearchResult::Found {
                calculator: None, ..
            }) => Err("Substance found but has no calculator".to_string()),
            Some(SearchResult::NotFound) => Err("Substance not found".to_string()),
            None => Err("Substance not in search results".to_string()),
        }
    }

    /// Print a summary of the search results with calculator types
    pub fn print_search_summary(&self) {
        println!("\nSearch Results Summary:");
        println!("----------------------");

        for (substance, result) in &self.search_results {
            match result {
                SearchResult::Found {
                    library,
                    priority_type,
                    calculator,
                    ..
                } => {
                    let priority_str = match priority_type {
                        LibraryPriority::Priority => "Priority",
                        LibraryPriority::Permitted => "Permitted",
                    };
                    let calc_str = match calculator {
                        Some(CalculatorType::Thermo(_)) => "Thermo",
                        Some(CalculatorType::Transport(_)) => "Transport",
                        None => "None",
                    };
                    println!(
                        "{}: Found in {} library ({}) - Calculator: {}",
                        substance, library, priority_str, calc_str
                    );
                }
                SearchResult::NotFound => {
                    println!("{}: Not found in any library", substance);
                }
            }
        }

        println!("\nStatistics:");
        println!("Total substances: {}", self.substances.len());
        println!(
            "Found in priority libraries: {}",
            self.get_priority_found_substances().len()
        );
        println!(
            "Found in permitted libraries: {}",
            self.get_permitted_found_substances().len()
        );
        println!("Not found: {}", self.get_not_found_substances().len());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_user_substances_with_calculators() {
        // Create a new UserSubstances instance with test substances
        let substances = vec!["CO".to_string(), "H2O".to_string()];
        let mut user_subs = UserSubstances::new(substances);

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Permitted,
        );

        // Perform the search
        user_subs.search_substances();

        // Test thermo calculations
        if let Ok((cp, dh, ds)) = user_subs.calculate_thermo_properties("CO", 400.0) {
            println!("CO Thermo properties at 400K:");
            println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
            assert!(cp > 0.0);
        }

        // Test transport calculations
        if let Ok((lambda, viscosity)) = user_subs.calculate_transport_properties("CO", 400.0) {
            println!("CO Transport properties at 400K:");
            println!("Lambda: {}, Viscosity: {}", lambda, viscosity);
            assert!(lambda > 0.0);
            assert!(viscosity > 0.0);
        }

        // Print full summary
        user_subs.print_search_summary();
    }
}
