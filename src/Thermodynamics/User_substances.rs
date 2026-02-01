//! # User Substances Module
//!
//! This module provides comprehensive management of thermodynamic and transport properties
//! for chemical substances across multiple databases. It implements a priority-based search
//! system and supports multiple data formats including numerical values, functions, and
//! symbolic expressions.
//!
//! ## Key Features
//!
//! - **Multi-Database Search**: Searches across NASA, NIST, CEA, Aramco and other databases
//! - **Priority System**: Hierarchical search with priority and permitted libraries
//! - **Multiple Data Formats**: Values, functions, and symbolic expressions
//! - **Comprehensive Properties**: Heat capacity, enthalpy, entropy, viscosity, thermal conductivity
//! - **Error Management**: Detailed logging and graceful error handling
//! - **Batch Processing**: Efficient handling of multiple substances
//!
//! ## Usage Example
//!
//! ```rust, ignore
//! use crate::Thermodynamics::User_substances::*;
//!
//! let mut subs_data = SubsData::new();
//! subs_data.set_substances(vec!["CO".to_string(), "CO2".to_string()]);
//! subs_data.set_multiple_library_priorities(
//!     vec!["NASA_gas".to_string()],
//!     LibraryPriority::Priority
//! );
//! subs_data.search_substances().unwrap();
//! subs_data.calculate_therm_map_of_properties(298.15).unwrap();
//! ```
/// Enumeration of all supported thermodynamic and transport property types
///
/// This enum covers three formats for each property:
/// - Direct values (e.g., `Cp`): Numerical results at specific conditions
/// - Functions (e.g., `Cp_fun`): Closures for temperature-dependent calculations  
/// - Symbolic (e.g., `Cp_sym`): Mathematical expressions for analytical work
use crate::Thermodynamics::DBhandlers::Diffusion::MultiSubstanceDiffusion;
/// agregator for all the thermo and transport calculations for user defined substances
use crate::Thermodynamics::DBhandlers::thermo_api::{
    ThermoCalculator, ThermoEnum, ThermoError, create_thermal_by_name,
};
use crate::Thermodynamics::DBhandlers::transport_api::{
    TransportCalculator, TransportEnum, create_transport_calculator_by_name,
};
use crate::Thermodynamics::User_substances_error::SimpleExceptionLogger;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use std::fmt;

//use RustedSciThe::symbolic::symbolic_engine::Expr;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use serde_json::Value;
use std::collections::HashMap;
/*
 * User Substances Module - Comprehensive Thermodynamic and Transport Property Management
 *
 * This module provides a unified interface for managing thermodynamic and transport properties
 * of chemical substances across multiple databases and libraries. It implements a sophisticated
 * search and priority system for substance data retrieval.
 *
 * Core Architecture:
 * ================
 *
 * Data Types:
 * - DataType: Enumeration of all supported property types (Cp, dH, dS, Lambda, Visc, etc.)
 *   with variants for values, functions, and symbolic expressions
 * - CalculatorType: Wrapper for thermodynamic or transport calculators
 * - WhatIsFound: Classification of search results (Thermo, Transport, NotFound)
 * - SearchResult: Complete search result with library info, priority, and calculator
 * - LibraryPriority: Priority classification (Priority, Permitted, Explicit)
 * - Phases: Physical phase enumeration (Liquid, Gas, Solid, Condensed)
 *
 * Main Structure - SubsData:
 * - Manages substance lists and library priorities
 * - Performs hierarchical searches across multiple databases
 * - Stores calculated properties in multiple formats (values, functions, symbolic)
 * - Provides comprehensive error logging and handling
 * - Supports both thermodynamic and transport property calculations
 *
 * Search Strategy:
 * ===============
 * 1. Priority Libraries: Search high-priority databases first
 * 2. Permitted Libraries: Fallback to secondary databases
 * 3. Explicit Instructions: Direct substance-to-library mappings
 * 4. Result Storage: Maintains search provenance and calculator instances
 *
 * Property Management:
 * ===================
 * - Values: Numerical results at specific conditions
 * - Functions: Closures for temperature-dependent calculations
 * - Symbolic: Mathematical expressions for analytical work
 * - Multi-substance: Batch processing capabilities
 *
 * Supported Libraries:
 * ===================
 * Thermodynamic: NASA, NIST, NASA7, NUIG_thermo, Cantera databases
 * Transport: CEA, Aramco_transport for viscosity, thermal conductivity, diffusion
 *
 * Error Handling:
 * ===============
 * - Comprehensive logging system with substance and function tracking
 * - Graceful degradation for missing data
 * - Detailed error propagation and reporting
 */

#[derive(Debug, Eq, PartialEq, Hash, Clone, Copy)]
#[allow(non_camel_case_types)]
pub enum DataType {
    /// Heat capacity at constant pressure [J/mol·K]
    Cp,
    /// Heat capacity function Cp(T)
    Cp_fun,
    /// Heat capacity symbolic expression
    Cp_sym,
    /// Enthalpy change [J/mol]
    dH,
    /// Enthalpy function dH(T)
    dH_fun,
    /// Enthalpy symbolic expression
    dH_sym,
    /// Entropy change [J/mol·K]
    dS,
    /// Entropy function dS(T)
    dS_fun,
    /// Entropy symbolic expression
    dS_sym,
    /// Chemical potential change [J/mol]
    dmu,
    /// Chemical potential function dmu(T)
    dmu_fun,
    /// Chemical potential symbolic expression
    dmu_sym,
    /// Thermal conductivity [W/m·K]
    Lambda,
    /// Thermal conductivity function Lambda(T)
    Lambda_fun,
    /// Thermal conductivity symbolic expression
    Lambda_sym,
    /// Dynamic viscosity [Pa·s]
    Visc,
    /// Viscosity function Visc(T)
    Visc_fun,
    /// Viscosity symbolic expression
    Visc_sym,
}
/// Wrapper enum for different types of property calculators
///
/// Encapsulates either thermodynamic calculators (for Cp, dH, dS) or
/// transport calculators (for viscosity, thermal conductivity, diffusion)
#[derive(Debug, Clone)]
pub enum CalculatorType {
    /// Thermodynamic property calculator (NASA, NIST, etc.)
    Thermo(ThermoEnum),
    /// Transport property calculator (CEA, Aramco, etc.)
    Transport(TransportEnum),
}
/// Classification of search results for substance data
///
/// Indicates what type of data was found during library searches
#[derive(Debug, Copy, Clone, Eq, Hash, PartialEq)]
pub enum WhatIsFound {
    /// Thermodynamic data found (Cp, dH, dS)
    Thermo,
    /// Transport data found (viscosity, thermal conductivity)
    Transport,
    /// No data found in any searched library
    NotFound,
}
/// Complete search result for a substance in a specific library
///
/// Contains all information needed to work with found substance data,
/// including the source library, priority level, raw data, and initialized calculator
#[derive(Debug, Clone)]
pub struct SearchResult {
    /// Name of the library where the substance was found
    pub library: String,
    /// Priority level of the library (Priority, Permitted, Explicit)
    pub priority_type: LibraryPriority,
    /// Raw JSON data from the library
    pub data: Value,
    /// Initialized calculator instance for property calculations
    pub calculator: Option<CalculatorType>,
}
/// Library priority levels for hierarchical searching
///
/// Defines the search order and preference for different databases
#[derive(Debug, Clone, PartialEq)]
pub enum LibraryPriority {
    /// High-priority libraries searched first (most trusted/accurate data)
    Priority,
    /// Secondary libraries used as fallback options
    Permitted,
    /// Direct substance-to-library mapping (bypasses normal search)
    Explicit,
}

/// Physical phases for substance classification
///
/// Used to specify the physical state when multiple phase data exists
#[derive(Debug, Clone, Copy)]
pub enum Phases {
    /// Liquid phase
    Liquid,
    /// Gas/vapor phase
    Gas,
    /// Solid/crystalline phase
    Solid,
    /// Condensed phase (liquid + solid)
    Condensed,
}
/// Main structure for comprehensive substance property management
///
/// `SubsData` is the central hub for managing thermodynamic and transport properties
/// of chemical substances. It provides:
///
/// - **Multi-database integration**: Searches across NASA, NIST, CEA, Aramco databases
/// - **Priority-based search**: Hierarchical search with configurable library priorities
/// - **Multiple data formats**: Stores properties as values, functions, and symbolic expressions
/// - **Comprehensive error handling**: Detailed logging and graceful error recovery
/// - **Batch processing**: Efficient handling of multiple substances simultaneously
///
/// ## Key Components
///
/// ### Search System
/// - `library_priorities`: Maps libraries to priority levels
/// - `search_results`: Stores found data with provenance information
/// - `explicit_search_instructions`: Direct substance-to-library mappings
///
/// ### Property Storage
/// - `therm_map_of_properties_values`: Numerical thermodynamic properties
/// - `therm_map_of_fun`: Function closures for temperature-dependent calculations
/// - `therm_map_of_sym`: Symbolic expressions for analytical work
/// - `transport_map_*`: Equivalent storage for transport properties
///
/// ### Physical Properties
/// - `P`, `T`: Pressure and temperature conditions
/// - `map_of_phases`: Physical phase information for each substance
/// - `hasmap_of_molar_mass`: Molecular weights
/// - `elem_composition_matrix`: Elemental composition data
///
/// ## Usage Pattern
///
/// 1. Create instance and set substances
/// 2. Configure library priorities
/// 3. Perform search across databases
/// 4. Calculate properties in desired formats
/// 5. Access results through getter methods
pub struct SubsData {
    /// System pressure [Pa] for property calculations
    pub P: Option<f64>,
    /// Pressure unit (e.g., "Pa", "atm", "bar")
    pub P_unit: Option<String>,
    /// System temperature [K] for property calculations
    pub T: Option<f64>,
    /// Molar mass unit (e.g., "g/mol", "kg/mol")
    pub Molar_mass_unit: Option<String>,
    /// List of chemical substances to search for and analyze
    pub substances: Vec<String>,
    /// Maps library names to their search priority levels
    pub library_priorities: HashMap<String, LibraryPriority>,
    /// Direct substance-to-library mappings that bypass normal search hierarchy
    pub explicit_search_insructions: HashMap<String, String>,
    /// Density values for each substance [kg/m³]
    pub ro_map: Option<HashMap<String, f64>>,
    /// Symbolic expressions for density calculations
    pub ro_map_sym: Option<HashMap<String, Box<Expr>>>,
    /// Physical phase specification for each substance
    pub map_of_phases: HashMap<String, Option<Phases>>,
    /// Complete search results with library provenance and calculators
    pub search_results: HashMap<String, HashMap<WhatIsFound, Option<SearchResult>>>,
    /// Interface to all thermodynamic and transport databases
    pub thermo_data: ThermoData,
    /// Calculated thermodynamic property values at specific conditions
    pub therm_map_of_properties_values: HashMap<String, HashMap<DataType, Option<f64>>>,
    /// Function closures for temperature-dependent thermodynamic calculations
    pub therm_map_of_fun:
        HashMap<String, HashMap<DataType, Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>>>,
    /// Symbolic expressions for analytical thermodynamic calculations
    pub therm_map_of_sym: HashMap<String, HashMap<DataType, Option<Box<Expr>>>>,
    /// Calculated transport property values at specific conditions
    pub transport_map_of_properties_values: HashMap<String, HashMap<DataType, Option<f64>>>,
    /// Function closures for temperature-dependent transport calculations
    pub transport_map_of_fun:
        HashMap<String, HashMap<DataType, Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>>>,
    /// Symbolic expressions for analytical transport calculations
    pub transport_map_of_sym: HashMap<String, HashMap<DataType, Option<Box<Expr>>>>,
    /// Matrix of elemental composition for all substances
    pub elem_composition_matrix: Option<DMatrix<f64>>,
    /// Molar masses for each substance [g/mol or specified unit]
    pub hasmap_of_molar_mass: HashMap<String, f64>,
    /// Multi-component diffusion coefficient data and calculations
    pub diffusion_data: Option<MultiSubstanceDiffusion>,
    /// List of unique chemical elements present in all substances
    pub unique_elements: Vec<String>,
    /// Error logging system for tracking calculation failures
    pub logger: SimpleExceptionLogger,
}
impl fmt::Debug for SubsData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SubsData")
            .field("substances", &self.substances)
            .field("library_priorities", &self.library_priorities)
            .field("search_results", &self.search_results)
            // Skip thermo_data as it might be large
            .field(
                "therm_map_of_properties_values",
                &self.therm_map_of_properties_values,
            )
            // Skip therm_map_of_fun as it contains closures that don't implement Debug
            .field("therm_map_of_sym", &self.therm_map_of_sym)
            .finish()
    }
}
impl Clone for SubsData {
    fn clone(&self) -> Self {
        // Create a new empty map for the functions
        let mut new_therm_map_of_fun = HashMap::new();

        // We can't directly clone the functions, so we'll need to handle this field specially
        // For now, we'll create an empty map structure that matches the original
        for (substance, type_map) in &self.therm_map_of_fun {
            let mut new_type_map = HashMap::new();
            for (data_type, _) in type_map {
                // We can't clone the functions, so we'll just insert None for each entry
                new_type_map.insert(*data_type, None);
            }
            new_therm_map_of_fun.insert(substance.clone(), new_type_map);
        }
        // Create a new empty map for the functions
        let new_transport_map_of_fun = HashMap::new();

        // We can't directly clone the functions, so we'll need to handle this field specially
        // For now, we'll create an empty map structure that matches the original
        for (substance, type_map) in &self.therm_map_of_fun {
            let mut new_type_map = HashMap::new();
            for (data_type, _) in type_map {
                // We can't clone the functions, so we'll just insert None for each entry
                new_type_map.insert(*data_type, None);
            }
            new_therm_map_of_fun.insert(substance.clone(), new_type_map);
        }
        let mut sd = SubsData {
            P: self.P,
            P_unit: self.P_unit.clone(),
            T: self.T.clone(),
            Molar_mass_unit: self.Molar_mass_unit.clone(),
            substances: self.substances.clone(),
            library_priorities: self.library_priorities.clone(),
            explicit_search_insructions: self.explicit_search_insructions.clone(),
            ro_map: self.ro_map.clone(),
            ro_map_sym: self.ro_map_sym.clone(),
            map_of_phases: self.map_of_phases.clone(),
            search_results: self.search_results.clone(),
            thermo_data: self.thermo_data.clone(),
            therm_map_of_properties_values: self.therm_map_of_properties_values.clone(),
            therm_map_of_fun: new_therm_map_of_fun,
            therm_map_of_sym: self.therm_map_of_sym.clone(),
            transport_map_of_properties_values: self.transport_map_of_properties_values.clone(),
            transport_map_of_fun: new_transport_map_of_fun,
            transport_map_of_sym: self.transport_map_of_sym.clone(),
            elem_composition_matrix: self.elem_composition_matrix.clone(),
            hasmap_of_molar_mass: self.hasmap_of_molar_mass.clone(),
            diffusion_data: self.diffusion_data.clone(),
            unique_elements: self.unique_elements.clone(),
            logger: self.logger.clone(),
        };
        let _ = sd.calculate_transport_map_of_functions();
        let _ = sd.calculate_therm_map_of_fun();
        sd
    }
}
impl SubsData {
    pub fn new() -> Self {
        Self {
            P: None,
            P_unit: None,
            T: None,
            Molar_mass_unit: None,
            substances: Vec::new(),
            library_priorities: HashMap::new(),
            explicit_search_insructions: HashMap::new(),
            search_results: HashMap::new(),
            ro_map: None,
            ro_map_sym: None,
            map_of_phases: HashMap::new(),
            thermo_data: ThermoData::new(),
            therm_map_of_properties_values: HashMap::new(),
            therm_map_of_fun: HashMap::new(),
            therm_map_of_sym: HashMap::new(),
            transport_map_of_properties_values: HashMap::new(),
            transport_map_of_fun: HashMap::new(),
            transport_map_of_sym: HashMap::new(),
            elem_composition_matrix: None,
            hasmap_of_molar_mass: HashMap::new(),
            diffusion_data: None,
            unique_elements: Vec::new(),
            logger: SimpleExceptionLogger::new(),
        }
    }

    pub fn empty() -> Self {
        Self {
            P: None,
            P_unit: None,
            T: None,
            Molar_mass_unit: None,
            substances: Vec::new(),
            library_priorities: HashMap::new(),
            explicit_search_insructions: HashMap::new(),
            search_results: HashMap::new(),
            ro_map: None,
            ro_map_sym: None,
            map_of_phases: HashMap::new(),
            thermo_data: ThermoData::empty(),
            therm_map_of_properties_values: HashMap::new(),
            therm_map_of_fun: HashMap::new(),
            therm_map_of_sym: HashMap::new(),
            transport_map_of_properties_values: HashMap::new(),
            transport_map_of_fun: HashMap::new(),
            transport_map_of_sym: HashMap::new(),
            elem_composition_matrix: None,
            hasmap_of_molar_mass: HashMap::new(),
            diffusion_data: None,
            unique_elements: Vec::new(),
            logger: SimpleExceptionLogger::new(),
        }
    }
    ////////////////////////////SETTERS////////////////////////////////////
    pub fn set_substances(&mut self, substances: Vec<String>) {
        self.substances = substances;
    }

    ////////////////////////SEARCH AND PRIORITY HANDLING///////////////////////////////

    /// Sets the priority level for a single library
    ///
    /// Libraries with `Priority` level are searched first, followed by `Permitted` libraries.
    /// `Explicit` priority is used for direct substance-to-library mappings.
    ///
    /// # Arguments
    ///
    /// * `library` - Name of the library (e.g., "NASA_gas", "CEA", "Aramco_transport")
    /// * `priority` - Priority level for this library
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
    /// subs_data.set_library_priority("NIST".to_string(), LibraryPriority::Permitted);
    /// ```
    pub fn set_library_priority(&mut self, library: String, priority: LibraryPriority) {
        self.library_priorities.insert(library, priority);
    }

    /// Sets the same priority level for multiple libraries simultaneously
    ///
    /// Convenient method for bulk priority assignment, commonly used to set
    /// all preferred libraries to `Priority` level at once.
    ///
    /// # Arguments
    ///
    /// * `libraries` - Vector of library names
    /// * `priority` - Priority level to assign to all libraries
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_multiple_library_priorities(
    ///     vec!["NASA_gas".to_string(), "NASA_cond".to_string()],
    ///     LibraryPriority::Priority
    /// );
    /// ```
    pub fn set_multiple_library_priorities(
        &mut self,
        libraries: Vec<String>,
        priority: LibraryPriority,
    ) {
        for library in libraries {
            self.library_priorities.insert(library, priority.clone());
        }
    }

    /// Sets explicit search instructions for direct substance-to-library mapping
    ///
    /// Bypasses the normal priority-based search for specified substances,
    /// forcing them to be searched only in their designated libraries.
    ///
    /// # Arguments
    ///
    /// * `direct_search` - HashMap mapping substance names to specific library names
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut explicit_map = HashMap::new();
    /// explicit_map.insert("H2O".to_string(), "NIST".to_string());
    /// explicit_map.insert("CO2".to_string(), "NASA_gas".to_string());
    /// subs_data.set_explicis_searh_instructions(explicit_map);
    /// ```
    pub fn set_explicis_searh_instructions(&mut self, direct_search: HashMap<String, String>) {
        self.explicit_search_insructions = direct_search;
    }

    /// Creates the appropriate calculator instance for a given library
    ///
    /// Determines whether a library contains thermodynamic or transport data
    /// and initializes the corresponding calculator type. This is used internally
    /// during the search process to prepare calculators for property calculations.
    ///
    /// # Arguments
    ///
    /// * `library` - Name of the library to create calculator for
    ///
    /// # Returns
    ///
    /// * `Some(CalculatorType)` - Appropriate calculator for the library
    /// * Panics if library type is not recognized
    ///
    /// # Supported Libraries
    ///
    /// **Thermodynamic**: NASA, NIST, NASA7, NUIG_thermo, Cantera databases
    /// **Transport**: CEA, Aramco_transport
    pub fn create_calculator(&self, library: &str) -> Option<CalculatorType> {
        match ThermoData::what_handler_to_use(library).as_str() {
            "NASA"
            | "NIST"
            | "NASA7"
            | "NUIG_thermo"
            | "Cantera_nasa_base_gas"
            | "Cantera_nasa_base_cond"
            | "NASA_gas"
            | "NASA_cond" => Some(CalculatorType::Thermo(create_thermal_by_name(library))),
            "transport" | "CEA" | "Aramco_transport" => Some(CalculatorType::Transport(
                create_transport_calculator_by_name(library),
            )),
            _ => {
                panic!("not found calculator for library: {}", library);
            }
        }
    }

    /// Performs comprehensive search for all substances across configured libraries
    ///
    /// Implements a hierarchical search strategy:
    /// 1. **Priority Libraries**: Searches high-priority databases first
    /// 2. **Permitted Libraries**: Falls back to secondary databases if not found
    /// 3. **Explicit Instructions**: Processes direct substance-to-library mappings
    /// 4. **Result Storage**: Stores all findings with provenance information
    ///
    /// For each found substance, initializes the appropriate calculator and stores
    /// the complete search result including library name, priority type, raw data,
    /// and calculator instance.
    ///
    /// # Returns
    ///
    /// * `Ok(())` - All searches completed successfully
    /// * `Err(SubsDataError)` - Error during search or calculator initialization
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut subs_data = SubsData::new();
    /// subs_data.set_substances(vec!["CO".to_string(), "H2O".to_string()]);
    /// subs_data.set_multiple_library_priorities(
    ///     vec!["NASA_gas".to_string()],
    ///     LibraryPriority::Priority
    /// );
    /// subs_data.search_substances()?;
    /// ```
    pub fn search_substances(&mut self) -> SubsDataResult<()> {
        println!("library priorities: {:?}  \n", self.library_priorities);
        // Get all priority libraries
        let priority_libs: Vec<String> = self
            .library_priorities
            .iter()
            .filter(move |&(_, &ref p)| *p == LibraryPriority::Priority)
            .map(|(lib, _)| lib.clone())
            .collect();
        println!(
            "\n search substances in Priority libraries: {} \n  ",
            priority_libs.join(", ")
        );
        // Get all permitted libraries
        let permitted_libs: Vec<String> = self
            .library_priorities
            .iter()
            .filter(move |&(_, &ref p)| *p == LibraryPriority::Permitted)
            .map(|(lib, _)| lib.clone())
            .collect();
        println!(
            "\n  search substances  in Permitted libraries: {} \n  ",
            permitted_libs.join(", ")
        );
        let subs_to_search_explicit_instruction: Vec<String> =
            self.explicit_search_insructions.keys().cloned().collect();
        let subs_to_search = &mut self.substances.clone();
        // Remove substances that have explicit search instructions
        subs_to_search.retain(|subs_i| !subs_to_search_explicit_instruction.contains(subs_i));
        // Search for each substance
        for substance in subs_to_search {
            println!("\n \n search substance {} \n \n ", substance);
            let mut found = false;
            // Try priority libraries first
            for lib in &priority_libs {
                println!(
                    "\n \n search substance {} in priority library: {} \n \n ",
                    substance, lib
                );
                if let Some(lib_data) = self.thermo_data.LibThermoData.get(lib) {
                    if let Some(substance_data) = lib_data.get(substance) {
                        // Create the calculator for this property
                        let mut calculator = self.create_calculator(lib);

                        // Initialize the calculator if one was created
                        if let Some(calc_type) = calculator.as_mut() {
                            match calc_type {
                                CalculatorType::Thermo(thermo) => {
                                    thermo.from_serde(substance_data.clone())?;
                                    print!("\n \n instance of thermal props struct:  ");
                                    let _ = thermo.print_instance();
                                    self.search_results
                                        .entry(substance.clone()) // insert a new entry if it doesn't exist
                                        .or_insert_with(HashMap::new)
                                        .insert(
                                            WhatIsFound::Thermo,
                                            Some(SearchResult {
                                                library: lib.clone(),
                                                priority_type: LibraryPriority::Priority,
                                                data: substance_data.clone(),
                                                calculator,
                                            }),
                                        );
                                }
                                CalculatorType::Transport(transport) => {
                                    transport.from_serde(substance_data.clone())?;
                                    print!("\n \n instance of transport props struct:  ");
                                    let _ = transport.print_instance();
                                    self.search_results
                                        .entry(substance.clone()) // insert a new entry if it doesn't exist
                                        .or_insert_with(HashMap::new)
                                        .insert(
                                            WhatIsFound::Transport,
                                            Some(SearchResult {
                                                library: lib.clone(),
                                                priority_type: LibraryPriority::Priority,
                                                data: substance_data.clone(),
                                                calculator,
                                            }),
                                        );
                                }
                            }
                        }

                        println!(
                            "\n  substance {} found in priority library {} \n  {:?}  \n",
                            substance,
                            lib,
                            self.search_results.get(substance).unwrap()
                        );
                        found = true;
                        //  break;
                    }
                }
            }

            // If not found in priority libraries, try permitted libraries
            if !found {
                for lib in &permitted_libs {
                    println!(
                        "\n \n search substance {} in Permitted library: {} \n \n ",
                        substance, lib
                    );
                    if let Some(lib_data) = self.thermo_data.LibThermoData.get(lib) {
                        if let Some(substance_data) = lib_data.get(substance) {
                            let mut calculator = self.create_calculator(lib);
                            //  let mut data_to_insert = HashMap::new();
                            // Initialize the calculator if one was created
                            if let Some(calc_type) = calculator.as_mut() {
                                match calc_type {
                                    CalculatorType::Thermo(thermo) => {
                                        thermo.newinstance()?;
                                        thermo.from_serde(substance_data.clone())?;
                                        let _ = thermo.print_instance();
                                        self.search_results
                                            .entry(substance.clone()) // insert a new entry if it doesn't exist
                                            .or_insert_with(HashMap::new)
                                            .insert(
                                                WhatIsFound::Thermo,
                                                Some(SearchResult {
                                                    library: lib.clone(),
                                                    priority_type: LibraryPriority::Permitted,
                                                    data: substance_data.clone(),
                                                    calculator,
                                                }),
                                            );
                                        assert!(
                                            &self
                                                .search_results
                                                .get(substance)
                                                .unwrap()
                                                .get(&WhatIsFound::Thermo)
                                                .unwrap()
                                                .is_some()
                                        );
                                    }
                                    CalculatorType::Transport(transport) => {
                                        //   let _ = transport.newinstance();
                                        transport.from_serde(substance_data.clone())?;
                                        let _ = transport.print_instance();

                                        self.search_results
                                            .entry(substance.clone()) // insert a new entry if it doesn't exist
                                            .or_insert_with(HashMap::new)
                                            .insert(
                                                WhatIsFound::Transport,
                                                Some(SearchResult {
                                                    library: lib.clone(),
                                                    priority_type: LibraryPriority::Permitted,
                                                    data: substance_data.clone(),
                                                    calculator,
                                                }),
                                            );

                                        assert!(
                                            &self
                                                .search_results
                                                .get(substance)
                                                .unwrap()
                                                .get(&WhatIsFound::Transport)
                                                .unwrap()
                                                .is_some()
                                        );
                                    }
                                }
                            }

                            println!(
                                "\n  substance {} found in permitted library {} \n  {:?}  \n",
                                substance,
                                lib,
                                self.search_results.get(substance).unwrap()
                            );
                            found = true;
                            //  break;
                        }
                    }
                }
            }

            // If still not found, mark as NotFound
            if !found {
                println!("\n  substance {} not found \n ", substance);
                let data_to_insert = HashMap::from([(WhatIsFound::NotFound, None)]);
                self.search_results
                    .insert(substance.clone(), data_to_insert);
            }
        } // end of for loop
        // now let's check if we have any explicit search instructions
        for (substance, library) in &self.explicit_search_insructions {
            let mut found = false;
            if let Some(lib_data) = self.thermo_data.LibThermoData.get(library) {
                if let Some(substance_data) = lib_data.get(substance) {
                    // Create the calculator for this property
                    let mut calculator = self.create_calculator(library);
                    // Initialize the calculator if one was created
                    if let Some(calc_type) = calculator.as_mut() {
                        match calc_type {
                            CalculatorType::Thermo(thermo) => {
                                thermo.from_serde(substance_data.clone())?;
                                print!("\n \n instance of thermal props struct:  ");
                                let _ = thermo.print_instance();
                                self.search_results
                                    .entry(substance.clone()) // insert a new entry if it doesn't exist
                                    .or_insert_with(HashMap::new)
                                    .insert(
                                        WhatIsFound::Thermo,
                                        Some(SearchResult {
                                            library: library.clone(),
                                            priority_type: LibraryPriority::Priority,
                                            data: substance_data.clone(),
                                            calculator,
                                        }),
                                    );
                            }
                            CalculatorType::Transport(transport) => {
                                transport.from_serde(substance_data.clone())?;
                                print!("\n \n instance of transport props struct:  ");
                                let _ = transport.print_instance();
                                self.search_results
                                    .entry(substance.clone()) // insert a new entry if it doesn't exist
                                    .or_insert_with(HashMap::new)
                                    .insert(
                                        WhatIsFound::Transport,
                                        Some(SearchResult {
                                            library: library.clone(),
                                            priority_type: LibraryPriority::Priority,
                                            data: substance_data.clone(),
                                            calculator,
                                        }),
                                    );
                            }
                        }
                    }

                    println!(
                        "\n  substance {} found in priority library {} \n  {:?}  \n",
                        substance,
                        library,
                        self.search_results.get(substance).unwrap()
                    );
                    found = true;
                    //  break;
                } // if let Some(substance_data) 
            } // if let Some(lib_data
            // If still not found, mark as NotFound
            if !found {
                println!("\n  substance {} not found \n ", substance);
                let data_to_insert = HashMap::from([(WhatIsFound::NotFound, None)]);
                self.search_results
                    .insert(substance.clone(), data_to_insert);
            }
        } // for loop
        Ok(())
    }
    /////////////////////////////////////GETTERS////////////////////////////////////
    /// Retrieves the complete search results for a specific substance
    ///
    /// Returns all found data types (thermodynamic and/or transport) for the substance,
    /// including library information and calculator instances.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance to query
    ///
    /// # Returns
    ///
    /// * `Some(HashMap)` - Map of data types to search results
    /// * `None` - Substance was not searched or found
    ///
    /// # Example
    ///
    /// ```rust
    /// if let Some(results) = subs_data.get_substance_result("CO2") {
    ///     if let Some(Some(thermo_result)) = results.get(&WhatIsFound::Thermo) {
    ///         println!("Found CO2 thermodynamic data in: {}", thermo_result.library);
    ///     }
    /// }
    /// ```
    pub fn get_substance_result(
        &self,
        substance: &str,
    ) -> Option<&HashMap<WhatIsFound, Option<SearchResult>>> {
        self.search_results.get(substance)
    }

    /// Retrieves mutable search results for a specific substance
    ///
    /// Provides mutable access to search results, allowing modification of
    /// calculator instances or result data. Used internally for calculator operations.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance to query
    ///
    /// # Returns
    ///
    /// * `Some(&mut HashMap)` - Mutable reference to search results
    /// * `None` - Substance was not searched or found
    pub fn get_substance_result_mut(
        &mut self,
        substance: &str,
    ) -> Option<&mut HashMap<WhatIsFound, Option<SearchResult>>> {
        self.search_results.get_mut(substance)
    }

    /// Returns all search results for all substances
    ///
    /// Provides access to the complete search results database, useful for
    /// comprehensive analysis or debugging of the search process.
    ///
    /// # Returns
    ///
    /// Reference to the complete search results HashMap
    pub fn get_all_results(&self) -> &HashMap<String, HashMap<WhatIsFound, Option<SearchResult>>> {
        &self.search_results
    }

    /// Returns a list of substances that were not found in any library
    ///
    /// Useful for identifying missing data and determining which substances
    /// need alternative data sources or manual input.
    ///
    /// # Returns
    ///
    /// Vector of substance names that were not found
    ///
    /// # Example
    ///
    /// ```rust
    /// let not_found = subs_data.get_not_found_substances();
    /// if !not_found.is_empty() {
    ///     println!("Missing data for: {:?}", not_found);
    /// }
    /// ```
    pub fn get_not_found_substances(&self) -> Vec<String> {
        self.search_results
            .iter()
            .filter(|(_, result)| result.get(&WhatIsFound::NotFound).is_some())
            .map(|(substance, _)| substance.clone())
            .collect()
    }

    /// Returns substances found in high-priority libraries
    ///
    /// Lists substances that were successfully found in libraries marked with
    /// `LibraryPriority::Priority`, indicating the most trusted data sources.
    ///
    /// # Returns
    ///
    /// Vector of substance names found in priority libraries
    ///
    /// # Example
    ///
    /// ```rust
    /// let priority_found = subs_data.get_priority_found_substances();
    /// println!("High-quality data available for: {:?}", priority_found);
    /// ```
    pub fn get_priority_found_substances(&self) -> Vec<String> {
        self.search_results
            .iter()
            .filter(|(_, result)| result.values().any(|result| {
                if let Some(result) = result.as_ref() {
                    matches!(result, SearchResult { priority_type, .. } if *priority_type == LibraryPriority::Priority)
                } else {
                    false
                }
                 }))
            .map(|(substance, _)| substance.clone())
            .collect()
    }

    /// Search substances by elements and populate search results
    pub fn search_by_elements(&mut self, elements: Vec<String>) -> SubsDataResult<Vec<String>> {
        let found_substances = self.thermo_data.search_by_elements(elements);
        self.substances = found_substances.clone();
        self.populate_element_search_results(found_substances)
    }

    /// Search substances containing only specified elements and populate search results
    pub fn search_by_elements_only(
        &mut self,
        elements: Vec<String>,
    ) -> SubsDataResult<Vec<String>> {
        let found_substances = self.thermo_data.search_by_elements_only(elements);
        self.substances = found_substances.clone();
        self.populate_element_search_results(found_substances)
    }

    /// Common function to populate search results from element-based searches
    fn populate_element_search_results(
        &mut self,
        found_substances: Vec<String>,
    ) -> SubsDataResult<Vec<String>> {
        for substance in &found_substances {
            if let Some(thermo_data) = self.thermo_data.hashmap_of_thermo_data.get(substance) {
                for (library, data) in thermo_data {
                    let mut calculator = self.create_calculator(library);

                    if let Some(CalculatorType::Thermo(thermo)) = calculator.as_mut() {
                        thermo.from_serde(data.clone())?;
                        self.search_results
                            .entry(substance.clone())
                            .or_insert_with(HashMap::new)
                            .insert(
                                WhatIsFound::Thermo,
                                Some(SearchResult {
                                    library: library.clone(),
                                    priority_type: LibraryPriority::Priority,
                                    data: data.clone(),
                                    calculator,
                                }),
                            );
                    }
                }
            }
        }
        Ok(found_substances)
    }
    ///
    /// Vector of substance names found in permitted libraries
    ///
    /// # Example
    ///
    /// ```rust
    /// let permitted_found = subs_data.get_permitted_found_substances();
    /// if !permitted_found.is_empty() {
    ///     println!("Secondary data sources used for: {:?}", permitted_found);
    /// }
    /// ```
    pub fn get_permitted_found_substances(&self) -> Vec<String> {
        self.search_results
        .iter()
        .filter(|(_, result)| result.values().any(|result| {
            if let Some(result) = result.as_ref() {
                matches!(result, SearchResult { priority_type, .. } if *priority_type == LibraryPriority::Permitted)
            } else {
                false
            }
             }))
        .map(|(substance, _)| substance.clone())
        .collect()
    }

    /// Retrieves a function closure for calculating thermodynamic properties
    ///
    /// Returns a temperature-dependent function for the specified property type.
    /// These functions are created from polynomial fits or other mathematical models.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance
    /// * `data_type` - Type of property function (e.g., `DataType::Cp_fun`)
    ///
    /// # Returns
    ///
    /// * `Some(function)` - Closure that takes temperature and returns property value
    /// * `None` - Function not available for this substance/property combination
    ///
    /// # Example
    ///
    /// ```rust
    /// if let Some(cp_func) = subs_data.get_thermo_function("CO2", DataType::Cp_fun) {
    ///     let cp_at_500k = cp_func(500.0); // Calculate Cp at 500 K
    /// }
    /// ```
    pub fn get_thermo_function(
        &self,
        substance: &str,
        data_type: DataType,
    ) -> Option<&Box<dyn Fn(f64) -> f64 + Send + Sync>> {
        self.therm_map_of_fun
            .get(substance)
            .and_then(|map| map.get(&data_type))
            .and_then(|opt| opt.as_ref())
    }

    /// Retrieves a symbolic expression for thermodynamic property calculations
    ///
    /// Returns a mathematical expression that can be manipulated symbolically,
    /// differentiated, integrated, or converted to other forms for analytical work.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance
    /// * `data_type` - Type of property expression (e.g., `DataType::Cp_sym`)
    ///
    /// # Returns
    ///
    /// * `Some(expression)` - Symbolic mathematical expression
    /// * `None` - Expression not available for this substance/property combination
    ///
    /// # Example
    ///
    /// ```rust
    /// if let Some(cp_expr) = subs_data.get_thermo_symbolic("H2O", DataType::Cp_sym) {
    ///     let cp_func = cp_expr.lambdify1D(); // Convert to evaluable function
    ///     let cp_derivative = cp_expr.diff("T"); // Symbolic differentiation
    /// }
    /// ```
    pub fn get_thermo_symbolic(&self, substance: &str, data_type: DataType) -> Option<&Box<Expr>> {
        self.therm_map_of_sym
            .get(substance)
            .and_then(|map| map.get(&data_type))
            .and_then(|opt| opt.as_ref())
    }
    ///////////////////////////////////////THERMAL PROPERTIES////////////////////////////////////

    /// Generic calculator access pattern - reduces boilerplate for calculator operations
    fn with_thermo_calculator<F, R>(&mut self, substance: &str, f: F) -> SubsDataResult<R>
    where
        F: FnOnce(&mut dyn ThermoCalculator) -> Result<R, ThermoError>,
    {
        let datamap = self
            .get_substance_result_mut(substance)
            .ok_or_else(|| SubsDataError::SubstanceNotFound(substance.to_string()))?;

        if datamap.contains_key(&WhatIsFound::NotFound) {
            return Err(SubsDataError::SubstanceNotFound(substance.to_string()));
        }

        let search_result = datamap
            .get_mut(&WhatIsFound::Thermo)
            .ok_or_else(|| SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Thermo".to_string(),
            })?
            .as_mut()
            .ok_or_else(|| SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Thermo".to_string(),
            })?;

        let calculator = search_result.calculator.as_mut().ok_or_else(|| {
            SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Thermo".to_string(),
            }
        })?;

        match calculator {
            CalculatorType::Thermo(thermo) => {
                let result: SubsDataResult<R> = f(thermo).map_err(Into::into);
                if let Err(ref e) = result {
                    crate::log_error!(self.logger, e, substance, "with_thermo_calculator");
                }
                result
            }
            CalculatorType::Transport(_) => Err(SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Thermo (found Transport instead)".to_string(),
            }),
        }
    }

    /// Extracts polynomial coefficients for thermodynamic calculations at a specific temperature
    ///
    /// Determines which temperature range applies and extracts the appropriate polynomial
    /// coefficients from the substance's thermodynamic data. This is typically required
    /// before performing property calculations.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance
    /// * `temperature` - Temperature [K] to determine coefficient range
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Coefficients extracted successfully
    /// * `Err(SubsDataError)` - Invalid temperature or substance not found
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.extract_thermal_coeffs("CO2", 500.0)?;
    /// ```
    pub fn extract_thermal_coeffs(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<()> {
        if temperature <= 0.0 {
            let error = SubsDataError::InvalidTemperature(temperature);
            return self.log_and_propagate(Err(error), substance, "extract_thermal_coeffs");
        }
        let result = self.with_thermo_calculator(substance, |thermo| {
            thermo.extract_model_coefficients(temperature)
        });
        self.log_and_propagate(result, substance, "extract_thermal_coeffs")
    }

    /// Extracts polynomial coefficients for all substances at a specific temperature
    ///
    /// Batch operation that extracts thermodynamic polynomial coefficients for all
    /// substances in the substance list. Useful for preparing multiple substances
    /// for property calculations at the same temperature.
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature [K] to determine coefficient ranges
    ///
    /// # Returns
    ///
    /// * `Ok(())` - All coefficients extracted successfully
    /// * `Err(SubsDataError)` - Invalid temperature or extraction failure
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.extract_all_thermal_coeffs(298.15)?; // Standard conditions
    /// ```
    pub fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> SubsDataResult<()> {
        if temperature <= 0.0 {
            let error = SubsDataError::InvalidTemperature(temperature);
            crate::log_error!(
                self.logger,
                &error,
                "ALL_SUBSTANCES",
                "extract_all_thermal_coeffs"
            );
            return Err(error);
        }
        for substance in self.substances.clone() {
            let result = self.extract_thermal_coeffs(&substance, temperature);
            self.log_and_propagate(result, &substance, "extract_all_thermal_coeffs")?;
        }
        Ok(())
    }
    /// we have right to calculate thermal properties, construct functions and symbolc expressions
    /// only if the thermal polynomial coefficients corresponds to the given temperature
    /// so this method if the current coefficients for this substance are suitable for the given temperature
    pub fn is_coeffs_valid_for_T(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<bool> {
        if temperature <= 0.0 {
            let error = SubsDataError::InvalidTemperature(temperature);
            return self.log_and_propagate(Err(error), substance, "is_coeffs_valid_for_T");
        }
        let result = self.with_thermo_calculator(substance, |thermo| {
            thermo.is_coeffs_valid_for_T(temperature)
        });
        self.log_and_propagate(result, substance, "is_coeffs_valid_for_T")
    }
    /// we have right to calculate thermal properties, construct functions and symbolc expressions
    /// only if the thermal polynomial coefficients corresponds to the given temperature
    /// so this method checks if the current coefficients for this substance are suitable for the given temperature
    /// and if no the coeddicients are renewed
    pub fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<(bool)> {
        let check_flag = self.is_coeffs_valid_for_T(substance, temperature);
        match check_flag {
            Ok(flag) => match flag {
                true => return Ok(true),
                false => {
                    self.extract_thermal_coeffs(substance, temperature)?;
                    return self.log_and_propagate(
                        check_flag,
                        substance,
                        "extract_coeffs_if_current_coeffs_not_valid",
                    );
                }
            },
            Err(error) => {
                crate::log_error!(
                    self.logger,
                    &error,
                    "ALL_SUBSTANCES",
                    "extract_coeffs_if_current_coeffs_not_valid"
                );
                return Err(error);
            }
        }
    }
    /// we have right to calculate thermal properties, construct functions and symbolc expressions
    /// only if the thermal polynomial coefficients corresponds to the given temperature
    /// so this method checks if the current coefficients for this substance are suitable for the given temperature
    /// and if no the coeddicients are renewed for all substances
    pub fn extract_coeffs_if_current_coeffs_not_valid_for_all_subs(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>> {
        let mut subs_with_changed_coeffs = Vec::new();
        for substance in self.substances.clone() {
            let result = self.extract_coeffs_if_current_coeffs_not_valid(&substance, temperature);
            match result {
                Ok(flag) => {
                    (if flag {
                        subs_with_changed_coeffs.push(substance.clone())
                    } else {
                    })
                }
                Err(_) => {}
            }

            self.log_and_propagate(
                result,
                &substance,
                "extract_coeffs_if_current_coeffs_not_valid_for_all_subs",
            )?;
        }
        Ok(subs_with_changed_coeffs)
    }

    /// Calculates )
    /// Sets the temperature range for thermodynamic calculations for a specific substance
    ///
    /// Defines the valid temperature interval for property calculations. This affects
    /// which polynomial coefficients are used and can enable temperature range validation.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance
    /// * `T_min` - Minimum temperature [K]
    /// * `T_max` - Maximum temperature [K]
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Temperature range set successfully
    /// * `Err(SubsDataError)` - Invalid range or substance not found
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_T_range_for_thermo("H2O", 273.15, 2000.0)?;
    /// ```
    pub fn set_T_range_for_thermo(
        &mut self,
        substance: &str,
        T_min: f64,
        T_max: f64,
    ) -> SubsDataResult<()> {
        let result =
            self.with_thermo_calculator(substance, |thermo| thermo.set_T_interval(T_min, T_max));
        self.log_and_propagate(result, substance, "set_T_range_for_thermo")
    }
    /// Sets the same temperature range for all substances' thermodynamic calculations
    ///
    /// Batch operation that applies the same temperature limits to all substances.
    /// Useful when all substances will be used within the same temperature range.
    ///
    /// # Arguments
    ///
    /// * `T_min` - Minimum temperature [K] for all substances
    /// * `T_max` - Maximum temperature [K] for all substances
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Temperature ranges set for all substances
    /// * `Err(SubsDataError)` - Invalid range or setting failure
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_T_range_for_all_thermo(300.0, 1500.0)?; // Combustion range
    /// ```
    pub fn set_T_range_for_all_thermo(&mut self, T_min: f64, T_max: f64) -> SubsDataResult<()> {
        for substance in self.substances.clone() {
            let result = self.set_T_range_for_thermo(&substance, T_min, T_max);
            self.log_and_propagate(result, &substance, "set_T_range_for_all_thermo")?;
        }
        Ok(())
    }
    /// Parses and validates thermodynamic polynomial coefficients for a substance
    ///
    /// Processes the raw coefficient data from the library, validates the polynomial
    /// structure, and prepares the coefficients for property calculations.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance to parse coefficients for
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Coefficients parsed successfully
    /// * `Err(SubsDataError)` - Parsing failure or substance not found
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.parse_thermal_coeffs("CH4")?;
    /// ```
    pub fn parse_thermal_coeffs(&mut self, substance: &str) -> SubsDataResult<()> {
        let result = self.with_thermo_calculator(substance, |thermo| thermo.parse_coefficients());
        self.log_and_propagate(result, substance, "parse_thermal_coeffs")
    }

    /// Parses thermodynamic coefficients for all substances in the substance list
    ///
    /// Batch operation that processes polynomial coefficients for all substances,
    /// preparing them for property calculations. This is typically done once after
    /// the search process is complete.
    ///
    /// # Returns
    ///
    /// * `Ok(())` - All coefficients parsed successfully
    /// * `Err(SubsDataError)` - Parsing failure for any substance
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.search_substances()?;
    /// subs_data.parse_all_thermal_coeffs()?;
    /// ```
    pub fn parse_all_thermal_coeffs(&mut self) -> SubsDataResult<()> {
        for substance in self.substances.clone() {
            let result = self.parse_thermal_coeffs(&substance);
            self.log_and_propagate(result, &substance, "parse_all_thermal_coeffs")?;
        }
        Ok(())
    }

    /// Fits polynomial coefficients to the specified temperature interval for a substance
    ///
    /// Adjusts or refits the thermodynamic polynomial coefficients to optimize accuracy
    /// within the specified temperature range. This can improve calculation precision
    /// for specific temperature intervals.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance to fit coefficients for
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Coefficients fitted successfully
    /// * `Err(SubsDataError)` - Fitting failure or substance not found
    ///
    /// # Note
    ///
    /// Temperature interval must be set using `set_T_range_for_thermo` before calling this method.
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_T_range_for_thermo("O2", 500.0, 1000.0)?;
    /// subs_data.fitting_thermal_coeffs_for_T_interval("O2")?;
    /// ```
    pub fn fitting_thermal_coeffs_for_T_interval(&mut self, substance: &str) -> SubsDataResult<()> {
        let result =
            self.with_thermo_calculator(substance, |thermo| thermo.fitting_coeffs_for_T_interval());
        self.log_and_propagate(result, substance, "fitting_thermal_coeffs_for_T_interval")
    }

    /// Fits polynomial coefficients for all substances to their specified temperature intervals
    ///
    /// Batch operation that optimizes thermodynamic polynomial coefficients for all substances
    /// within their respective temperature ranges. Improves calculation accuracy across
    /// the entire substance set.
    ///
    /// # Returns
    ///
    /// * `Ok(())` - All coefficients fitted successfully
    /// * `Err(SubsDataError)` - Fitting failure for any substance
    ///
    /// # Note
    ///
    /// Temperature intervals must be set for all substances before calling this method.
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_T_range_for_all_thermo(400.0, 1200.0)?;
    /// subs_data.fitting_all_thermal_coeffs_for_T_interval()?;
    /// ```
    pub fn fitting_all_thermal_coeffs_for_T_interval(&mut self) -> SubsDataResult<()> {
        for substance in self.substances.clone() {
            let result = self.fitting_thermal_coeffs_for_T_interval(&substance);
            self.log_and_propagate(
                result,
                &substance,
                "fitting_all_thermal_coeffs_for_T_interval",
            )?;
        }
        Ok(())
    }

    /// Calculates integral mean values for thermodynamic properties over temperature range
    ///
    /// Computes temperature-averaged values of thermodynamic properties (Cp, dH, dS)
    /// over the specified temperature interval. Useful for mean property calculations
    /// in process design and analysis.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance to calculate integral means for
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Integral means calculated successfully
    /// * `Err(SubsDataError)` - Calculation failure or substance not found
    ///
    /// # Note
    ///
    /// Temperature interval must be set before calling this method.
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_T_range_for_thermo("N2", 300.0, 800.0)?;
    /// subs_data.integr_mean("N2")?;
    /// ```
    pub fn integr_mean(&mut self, substance: &str) -> SubsDataResult<()> {
        let result = self.with_thermo_calculator(substance, |thermo| thermo.integr_mean());
        self.log_and_propagate(result, substance, "integr_mean")
    }

    /// Calculates thermodynamic properties (Cp, dH, dS) for a substance at given temperature
    ///
    /// Computes heat capacity, enthalpy change, and entropy change for the specified
    /// substance at the given temperature using the appropriate thermodynamic model.
    ///
    /// # Arguments
    ///
    /// * `substance` - Name of the substance
    /// * `temperature` - Temperature [K] for property calculation
    ///
    /// # Returns
    ///
    /// * `Ok((cp, dh, ds))` - Tuple of (heat capacity [J/mol·K], enthalpy change [J/mol], entropy change [J/mol·K])
    /// * `Err(SubsDataError)` - Calculation failure, invalid temperature, or substance not found
    ///
    /// # Example
    ///
    /// ```rust
    /// let (cp, dh, ds) = subs_data.calculate_thermo_properties("CO2", 500.0)?;
    /// println!("Cp: {:.2} J/mol·K, dH: {:.2} J/mol, dS: {:.2} J/mol·K", cp, dh, ds);
    /// ```
    pub fn calculate_thermo_properties(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<(f64, f64, f64)> {
        let result = self._calculate_thermo_properties_internal(substance, temperature);
        self.log_and_propagate(result, substance, "calculate_thermo_properties")
    }

    fn _calculate_thermo_properties_internal(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<(f64, f64, f64)> {
        if temperature <= 0.0 {
            return Err(SubsDataError::InvalidTemperature(temperature));
        }

        let datamap = self
            .search_results
            .get(substance)
            .ok_or_else(|| SubsDataError::SubstanceNotFound(substance.to_string()))?;

        // Check if substance was marked as not found
        if datamap.contains_key(&WhatIsFound::NotFound) {
            return Err(SubsDataError::SubstanceNotFound(substance.to_string()));
        }

        let search_result = datamap.get(&WhatIsFound::Thermo).ok_or_else(|| {
            SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Thermo".to_string(),
            }
        })?;

        match search_result {
            Some(SearchResult {
                calculator: Some(CalculatorType::Thermo(thermo)),
                ..
            }) => {
                let mut thermo = thermo.clone();

                thermo.calculate_Cp_dH_dS(temperature)?;

                let cp = thermo.get_Cp()?;
                let dh = thermo.get_dh()?;
                let ds = thermo.get_ds()?;

                Ok((cp, dh, ds))
            }
            Some(SearchResult {
                calculator: Some(CalculatorType::Transport(_)),
                ..
            }) => Err(SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Thermo (found Transport instead)".to_string(),
            }),
            Some(SearchResult {
                calculator: None, ..
            }) => Err(SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Thermo".to_string(),
            }),
            None => Err(SubsDataError::SubstanceNotFound(substance.to_string())),
        }
    }

    /// Sets molar masses for substances with optional unit specification
    ///
    /// Assigns molar mass values to substances, which are used in density calculations,
    /// unit conversions, and other mass-based property calculations.
    ///
    /// # Arguments
    ///
    /// * `M_map` - HashMap mapping substance names to their molar masses
    /// * `M_unit` - Optional unit specification (e.g., "g/mol", "kg/mol")
    ///
    /// # Example
    ///
    /// ```rust
    /// let mut molar_masses = HashMap::new();
    /// molar_masses.insert("CO2".to_string(), 44.01);
    /// molar_masses.insert("H2O".to_string(), 18.015);
    /// subs_data.set_M(molar_masses, Some("g/mol".to_string()));
    /// ```
    pub fn set_M(&mut self, M_map: HashMap<String, f64>, M_unit: Option<String>) {
        self.hasmap_of_molar_mass = M_map;
        self.Molar_mass_unit = M_unit;
    }
    /// Sets system pressure with optional unit specification
    ///
    /// Defines the pressure for property calculations. This affects density calculations,
    /// phase equilibrium, and pressure-dependent transport properties.
    ///
    /// # Arguments
    ///
    /// * `P` - Pressure value
    /// * `P_unit` - Optional unit specification (e.g., "Pa", "atm", "bar")
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_P(101325.0, Some("Pa".to_string())); // Standard atmospheric pressure
    /// subs_data.set_P(1.0, Some("atm".to_string()));      // Alternative specification
    /// ```
    pub fn set_P(&mut self, P: f64, P_unit: Option<String>) {
        self.P = Some(P);
        self.P_unit = P_unit;
    }

    /// Sets system temperature for property calculations
    ///
    /// Defines the default temperature for property calculations. This is used
    /// when temperature is not explicitly specified in calculation methods.
    ///
    /// # Arguments
    ///
    /// * `T` - Temperature [K]
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.set_T(298.15); // Standard temperature (25°C)
    /// ```
    pub fn set_T(&mut self, T: f64) {
        self.T = Some(T);
    }

    /// Generic utility for building property maps with error handling
    ///
    /// Internal utility function that reduces boilerplate code for property calculations
    /// across multiple substances. Handles errors gracefully by inserting empty values
    /// for failed calculations while continuing with successful ones.
    ///
    /// # Type Parameters
    ///
    /// * `T` - Type of property values (f64, function closures, expressions)
    /// * `F` - Calculator function type
    ///
    /// # Arguments
    ///
    /// * `calculator_fn` - Function that calculates properties for a single substance
    /// * `target_map` - Map to store calculated properties
    /// * `property_types` - Types of properties to calculate
    /// * `empty_value_fn` - Function to generate empty values for failed calculations
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Property map built successfully (with possible individual failures)
    /// * `Err(SubsDataError)` - Critical error in map building process
    pub fn build_property_map<T, F>(
        &mut self,
        calculator_fn: F,
        target_map: &mut HashMap<String, HashMap<DataType, Option<T>>>,
        property_types: &[DataType],
        empty_value_fn: fn() -> Option<T>,
    ) -> SubsDataResult<()>
    where
        F: Fn(&mut Self, &str) -> SubsDataResult<HashMap<DataType, Option<T>>>,
    {
        for substance in self.substances.clone() {
            match calculator_fn(self, &substance) {
                Ok(properties) => {
                    target_map.insert(substance, properties);
                }
                Err(e) => {
                    let mut empty_map = HashMap::new();
                    for &prop_type in property_types {
                        empty_map.insert(prop_type, empty_value_fn());
                    }
                    target_map.insert(substance.clone(), empty_map);
                    println!(
                        "Warning: Failed to calculate properties for {}: {}",
                        substance, e
                    );
                }
            }
        }
        Ok(())
    }

    /// Helper function to calculate thermo properties for a single substance
    fn calculate_single_thermo_properties(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<HashMap<DataType, Option<f64>>> {
        let (cp, dh, ds) = self._calculate_thermo_properties_internal(substance, temperature)?;
        let mut property_map = HashMap::new();
        property_map.insert(DataType::Cp, Some(cp));
        property_map.insert(DataType::dH, Some(dh));
        property_map.insert(DataType::dS, Some(ds));
        Ok(property_map)
    }

    /// Calculates and stores thermodynamic property values for all substances at specified temperature
    ///
    /// Batch calculation that computes Cp, dH, and dS for all substances at the given temperature
    /// and stores the results in `therm_map_of_properties_values`. This is the primary method
    /// for obtaining numerical property values.
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature [K] for property calculations
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Properties calculated and stored successfully
    /// * `Err(SubsDataError)` - Invalid temperature or calculation failure
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.search_substances()?;
    /// subs_data.calculate_therm_map_of_properties(500.0)?;
    ///
    /// // Access calculated values
    /// if let Some(co2_props) = subs_data.therm_map_of_properties_values.get("CO2") {
    ///     if let Some(Some(cp)) = co2_props.get(&DataType::Cp) {
    ///         println!("CO2 heat capacity at 500K: {:.2} J/mol·K", cp);
    ///     }
    /// }
    /// ```
    pub fn calculate_therm_map_of_properties(&mut self, temperature: f64) -> SubsDataResult<()> {
        let mut temp_map = HashMap::new();
        for substance in self.substances.clone() {
            let result = self.calculate_single_thermo_properties(&substance, temperature);
            match self.log_and_propagate(result, &substance, "calculate_therm_map_of_properties") {
                Ok(properties) => {
                    temp_map.insert(substance, properties);
                }
                Err(_) => {
                    let mut empty_map = HashMap::new();
                    empty_map.insert(DataType::Cp, None);
                    empty_map.insert(DataType::dH, None);
                    empty_map.insert(DataType::dS, None);
                    temp_map.insert(substance, empty_map);
                }
            }
        }
        self.therm_map_of_properties_values = temp_map;
        Ok(())
    }

    /// Helper function to calculate thermo functions for a single substance
    fn calculate_single_thermo_functions(
        &mut self,
        substance: &str,
    ) -> SubsDataResult<HashMap<DataType, Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>>> {
        self.with_thermo_calculator(substance, |thermo| {
            thermo.create_closures_Cp_dH_dS()?;
            let mut function_map = HashMap::new();

            function_map.insert(DataType::Cp_fun, thermo.get_C_fun().ok());
            function_map.insert(DataType::dH_fun, thermo.get_dh_fun().ok());
            function_map.insert(DataType::dS_fun, thermo.get_ds_fun().ok());

            Ok(function_map)
        })
    }

    /// Creates and stores function closures for thermodynamic property calculations
    ///
    /// Generates temperature-dependent function closures (Cp(T), dH(T), dS(T)) for all substances
    /// and stores them in `therm_map_of_fun`. These functions can be used for efficient
    /// property calculations at multiple temperatures.
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Function closures created and stored successfully
    /// * `Err(SubsDataError)` - Function creation failure
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.search_substances()?;
    /// subs_data.calculate_therm_map_of_fun()?;
    ///
    /// // Use generated functions
    /// if let Some(cp_func) = subs_data.get_thermo_function("H2O", DataType::Cp_fun) {
    ///     let cp_300k = cp_func(300.0);
    ///     let cp_400k = cp_func(400.0);
    ///     println!("H2O Cp: {:.2} at 300K, {:.2} at 400K", cp_300k, cp_400k);
    /// }
    /// ```
    pub fn calculate_therm_map_of_fun(&mut self) -> SubsDataResult<()> {
        let mut temp_map = HashMap::new();
        for substance in self.substances.clone() {
            let result = self.calculate_single_thermo_functions(&substance);
            match self.log_and_propagate(result, &substance, "calculate_therm_map_of_fun") {
                Ok(properties) => {
                    temp_map.insert(substance, properties);
                }
                Err(_) => {
                    let mut empty_map = HashMap::new();
                    empty_map.insert(DataType::Cp_fun, None);
                    empty_map.insert(DataType::dH_fun, None);
                    empty_map.insert(DataType::dS_fun, None);
                    temp_map.insert(substance, empty_map);
                }
            }
        }
        self.therm_map_of_fun = temp_map;
        Ok(())
    }

    /// Helper function to calculate thermo symbolic expressions for a single substance
    fn calculate_single_thermo_symbolic(
        &mut self,
        substance: &str,
    ) -> SubsDataResult<HashMap<DataType, Option<Box<Expr>>>> {
        self.with_thermo_calculator(substance, |thermo| {
            thermo.create_sym_Cp_dH_dS()?;
            let mut sym_map = HashMap::new();

            sym_map.insert(DataType::Cp_sym, thermo.get_Cp_sym().ok().map(Box::new));
            sym_map.insert(DataType::dH_sym, thermo.get_dh_sym().ok().map(Box::new));
            sym_map.insert(DataType::dS_sym, thermo.get_ds_sym().ok().map(Box::new));

            Ok(sym_map)
        })
    }

    /// Creates and stores symbolic expressions for thermodynamic property calculations
    ///
    /// Generates symbolic mathematical expressions for thermodynamic properties (Cp, dH, dS)
    /// and stores them in `therm_map_of_sym`. These expressions can be manipulated symbolically,
    /// differentiated, integrated, or used for analytical work.
    ///
    /// # Returns
    ///
    /// * `Ok(())` - Symbolic expressions created and stored successfully
    /// * `Err(SubsDataError)` - Expression creation failure
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.search_substances()?;
    /// subs_data.calculate_therm_map_of_sym()?;
    ///
    /// // Use symbolic expressions
    /// if let Some(cp_expr) = subs_data.get_thermo_symbolic("CH4", DataType::Cp_sym) {
    ///     let cp_func = cp_expr.lambdify1D();           // Convert to function
    ///     let cp_derivative = cp_expr.diff("T");        // Symbolic differentiation
    ///     let cp_integral = cp_expr.integrate("T");     // Symbolic integration
    /// }
    /// ```
    pub fn calculate_therm_map_of_sym(&mut self) -> SubsDataResult<()> {
        let mut temp_map = HashMap::new();
        for substance in self.substances.clone() {
            let result = self.calculate_single_thermo_symbolic(&substance);
            match self.log_and_propagate(result, &substance, "calculate_therm_map_of_sym") {
                Ok(properties) => {
                    temp_map.insert(substance, properties);
                }
                Err(_) => {
                    let mut empty_map = HashMap::new();
                    empty_map.insert(DataType::Cp_sym, None);
                    empty_map.insert(DataType::dH_sym, None);
                    empty_map.insert(DataType::dS_sym, None);
                    temp_map.insert(substance, empty_map);
                }
            }
        }
        self.therm_map_of_sym = temp_map;
        Ok(())
    }

    /// Internal helper method for error logging and propagation
    ///
    /// Logs errors to the internal logger while preserving the original error for propagation.
    /// This enables comprehensive error tracking while maintaining normal error handling flow.
    ///
    /// # Type Parameters
    ///
    /// * `T` - Type of the result value
    ///
    /// # Arguments
    ///
    /// * `result` - Result to log and propagate
    /// * `substance` - Name of substance associated with the operation
    /// * `function` - Name of function where error occurred
    ///
    /// # Returns
    ///
    /// The original result, unchanged, after logging any errors
    pub fn log_and_propagate<T>(
        &mut self,
        result: SubsDataResult<T>,
        substance: &str,
        function: &str,
    ) -> SubsDataResult<T> {
        if let Err(ref e) = result {
            crate::log_error!(self.logger, e, substance, function);
        }
        result
    }

    /// Prints all logged errors to the console
    ///
    /// Displays a comprehensive list of all errors that occurred during substance
    /// processing, including substance names, function names, and error details.
    /// Useful for debugging and error analysis.
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.search_substances().ok(); // May generate errors
    /// subs_data.calculate_therm_map_of_properties(500.0).ok();
    /// subs_data.print_error_logs(); // Display any errors that occurred
    /// ```
    pub fn print_error_logs(&self) {
        use crate::Thermodynamics::User_substances_error::ExceptionLogger;
        self.logger.print_logs();
    }

    /// Clears all logged errors from the internal logger
    ///
    /// Removes all previously logged errors, providing a clean slate for subsequent
    /// operations. Useful when reusing the same `SubsData` instance for multiple
    /// calculation sessions.
    ///
    /// # Example
    ///
    /// ```rust
    /// subs_data.print_error_logs();  // Show current errors
    /// subs_data.clear_error_logs();  // Clear error history
    /// // Perform new calculations with clean error log
    /// ```
    pub fn clear_error_logs(&mut self) {
        use crate::Thermodynamics::User_substances_error::ExceptionLogger;
        self.logger.clear_logs();
    }
}
