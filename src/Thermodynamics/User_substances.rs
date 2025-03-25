use crate::Thermodynamics::DBhandlers::transport_api::TransportCalculator;
use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use RustedSciThe::symbolic::symbolic_engine::Expr;
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

///SearchResult enum: Represents the result of searching for a substance
#[derive(Debug)]
pub enum SearchResult {
    Found {
        library: String,
        priority_type: LibraryPriority,
        data: Value,
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
    pub fn set_multiple_library_priorities(&mut self, libraries: Vec<String>, priority: LibraryPriority) {
        for library in libraries {
            self.library_priorities.insert(library, priority.clone());
        }
    }

    /// Search for all substances in the libraries according to priority
    pub fn search_substances(&mut self) {
        // Get all priority libraries
        let priority_libs: Vec<String> = self.library_priorities
            .iter()
            .filter(move |&(_, &ref p)| *p == LibraryPriority::Priority)
            .map(|(lib, _)| lib.clone())
            .collect();

        // Get all permitted libraries
        let permitted_libs: Vec<String> = self.library_priorities
            .iter()
            .filter(move|&(_, &ref p)| *p == LibraryPriority::Permitted)
            .map(|(lib, _)| lib.clone())
            .collect();

        // Search for each substance
        for substance in &self.substances {
            // First search in priority libraries
            let mut found = false;
            
            // Try priority libraries first
            for lib in &priority_libs {
                if let Some(lib_data) = self.thermo_data.LibThermoData.get(lib) {
                    if let Some(substance_data) = lib_data.get(substance) {
                        self.search_results.insert(
                            substance.clone(),
                            SearchResult::Found {
                                library: lib.clone(),
                                priority_type: LibraryPriority::Priority,
                                data: substance_data.clone(),
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
                            self.search_results.insert(
                                substance.clone(),
                                SearchResult::Found {
                                    library: lib.clone(),
                                    priority_type: LibraryPriority::Permitted,
                                    data: substance_data.clone(),
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
                self.search_results.insert(substance.clone(), SearchResult::NotFound);
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

    /// Print a summary of the search results
    pub fn print_search_summary(&self) {
        println!("\nSearch Results Summary:");
        println!("----------------------");
        
        for (substance, result) in &self.search_results {
            match result {
                SearchResult::Found { library, priority_type, .. } => {
                    let priority_str = match priority_type {
                        LibraryPriority::Priority => "Priority",
                        LibraryPriority::Permitted => "Permitted",
                    };
                    println!("{}: Found in {} library ({})", substance, library, priority_str);
                }
                SearchResult::NotFound => {
                    println!("{}: Not found in any library", substance);
                }
            }
        }
        
        println!("\nStatistics:");
        println!("Total substances: {}", self.substances.len());
        println!("Found in priority libraries: {}", self.get_priority_found_substances().len());
        println!("Found in permitted libraries: {}", self.get_permitted_found_substances().len());
        println!("Not found: {}", self.get_not_found_substances().len());
    } 
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_user_substances() {
        // Create a new UserSubstances instance with some test substances
        let substances = vec!["CO".to_string(), "H2O".to_string(), "N2".to_string()];
        let mut user_subs = UserSubstances::new(substances);

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string(), "NIST".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string(), "CEA".to_string()],
            LibraryPriority::Permitted,
        );

        // Perform the search
        user_subs.search_substances();

        // Print the results
        user_subs.print_search_summary();

        // Test specific substance results
        if let Some(SearchResult::Found { library, priority_type, .. }) = user_subs.get_substance_result("CO") {
            assert!(matches!(priority_type, LibraryPriority::Priority));
            assert!(library == "NASA_gas" || library == "NIST");
        }

        // Verify that we can get lists of substances by their search status
        let priority_found = user_subs.get_priority_found_substances();
        let permitted_found = user_subs.get_permitted_found_substances();
        let not_found = user_subs.get_not_found_substances();

        // The sum of these should equal the total number of substances
        assert_eq!(
            priority_found.len() + permitted_found.len() + not_found.len(),
            3
        );
    }
}