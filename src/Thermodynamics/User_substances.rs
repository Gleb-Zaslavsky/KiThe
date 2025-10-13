/// agregator for all the thermo and transport calculations for user defined substances
use crate::Thermodynamics::DBhandlers::thermo_api::{
    ThermoCalculator, ThermoEnum, create_thermal_by_name,
};
use crate::Thermodynamics::DBhandlers::transport_api::{
    TransportCalculator, TransportEnum, create_transport_calculator_by_name,
};
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use std::fmt;

//use RustedSciThe::symbolic::symbolic_engine::Expr;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use serde_json::Value;
use std::collections::HashMap;
/*
Core Types:
        SearchResult enum: Represents the result of searching for a substance
        LibraryPriority enum: Defines priority levels (Priority or Permitted)
        SubsData struct: Main structure managing the substances and search process
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

#[derive(Debug, Eq, PartialEq, Hash, Clone, Copy)]
#[allow(non_camel_case_types)]
pub enum DataType {
    Cp,

    Cp_fun,

    Cp_sym,

    dH,

    dH_fun,

    dH_sym,

    dS,

    dS_fun,

    dS_sym,

    dmu,

    dmu_fun,

    dmu_sym,

    Lambda,

    Lambda_fun,

    Lambda_sym,

    Visc,

    Visc_fun,

    Visc_sym,
}
/// Represents the type of calculator for a substance
#[derive(Debug, Clone)]
pub enum CalculatorType {
    Thermo(ThermoEnum),
    Transport(TransportEnum),
}
#[derive(Debug, Copy, Clone, Eq, Hash, PartialEq)]
pub enum WhatIsFound {
    Thermo,
    Transport,
    NotFound,
}
///SearchResult enum: Represents the result of searching for a substance
#[derive(Debug, Clone)]
pub struct SearchResult {
    pub library: String,
    pub priority_type: LibraryPriority,
    pub data: Value,
    pub calculator: Option<CalculatorType>,
}
///LibraryPriority enum: Defines priority levels (Priority or Permitted)
#[derive(Debug, Clone, PartialEq)]
pub enum LibraryPriority {
    Priority,
    Permitted,
    Explicit,
}

#[derive(Debug, Clone, Copy)]
pub enum Phases {
    Liquid,
    Gas,
    Solid,
    Condensed,
}
///SubsData struct: Main structure managing the substances and search process

pub struct SubsData {
    pub P: Option<f64>,
    pub P_unit: Option<String>,

    pub Molar_mass_unit: Option<String>,
    /// Vector of substances to search for
    pub substances: Vec<String>,
    /// Maps libraries to their priority level
    pub library_priorities: HashMap<String, LibraryPriority>,
    ///
    pub explicit_search_insructions: HashMap<String, String>,
    pub ro_map: Option<HashMap<String, f64>>,
    ///
    pub ro_map_sym: Option<HashMap<String, Box<Expr>>>,
    /// hashmap of phases of substances
    pub map_of_phases: HashMap<String, Option<Phases>>,
    /// Stores search results for each substance
    pub search_results: HashMap<String, HashMap<WhatIsFound, Option<SearchResult>>>,
    /// Reference to the thermo data library
    pub thermo_data: ThermoData,
    ///  calaculated thermodynamic properties
    pub therm_map_of_properties_values: HashMap<String, HashMap<DataType, Option<f64>>>,
    /// thermodynamic properties as functions
    pub therm_map_of_fun: HashMap<String, HashMap<DataType, Option<Box<dyn Fn(f64) -> f64>>>>,
    /// thermodynamic properties as symbolic expressions
    pub therm_map_of_sym: HashMap<String, HashMap<DataType, Option<Box<Expr>>>>,
    /// calculated transport properties
    pub transport_map_of_properties_values: HashMap<String, HashMap<DataType, Option<f64>>>,
    /// transport properties as functions
    pub transport_map_of_fun: HashMap<String, HashMap<DataType, Option<Box<dyn Fn(f64) -> f64>>>>,
    /// transport properties as symbolic expressions
    pub transport_map_of_sym: HashMap<String, HashMap<DataType, Option<Box<Expr>>>>,
    ///
    pub elem_composition_matrix: Option<DMatrix<f64>>,
    ///
    pub hasmap_of_molar_mass: HashMap<String, f64>,
    ///
    pub unique_elements: Vec<String>,
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
            unique_elements: self.unique_elements.clone(),
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
            unique_elements: Vec::new(),
        }
    }
    ////////////////////////SEARCH AND PRIORITY HANDLING///////////////////////////////
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
    ///
    pub fn set_explicis_searh_instructions(&mut self, direct_search: HashMap<String, String>) {
        self.explicit_search_insructions = direct_search;
    }
    /// Initialize appropriate calculator for a substance based on its library
    /// calculator is ThermoEnum or TransportEnum
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

    /// Search for all substances in the libraries according to priority
    pub fn search_substances(&mut self) {
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
                                    let _ = thermo.from_serde(substance_data.clone());
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
                                    let _ = transport.from_serde(substance_data.clone());
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
                                        let _ = thermo.newinstance();
                                        let _ = thermo.from_serde(substance_data.clone());
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
                                        let _ = transport.from_serde(substance_data.clone());
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
                                let _ = thermo.from_serde(substance_data.clone());
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
                                let _ = transport.from_serde(substance_data.clone());
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
    }
    /////////////////////////////////////GETTERS////////////////////////////////////
    /// Get the search result for a specific substance
    pub fn get_substance_result(
        &self,
        substance: &str,
    ) -> Option<&HashMap<WhatIsFound, Option<SearchResult>>> {
        //Option<&SearchResult> {
        self.search_results.get(substance)
    }

    pub fn get_substance_result_mut(
        &mut self, // Changed to mutable reference
        substance: &str,
    ) -> Option<&mut HashMap<WhatIsFound, Option<SearchResult>>> {
        self.search_results.get_mut(substance) // Use get_mut instead of get
    }
    /// Get all search results
    pub fn get_all_results(&self) -> &HashMap<String, HashMap<WhatIsFound, Option<SearchResult>>> {
        &self.search_results
    }

    /// Get a list of substances that were not found
    pub fn get_not_found_substances(&self) -> Vec<String> {
        self.search_results
            .iter()
            .filter(|(_, result)| result.get(&WhatIsFound::NotFound).is_some())
            .map(|(substance, _)| substance.clone())
            .collect()
    }

    /// Get a list of substances found in priority libraries
    pub fn get_priority_found_substances(&self) -> Vec<String> {
        //search_results = HashMap<String, HashMap<WhatIsFound, Option<SearchResult>>>,
        self.search_results
            .iter()//  HashMap<WhatIsFound, Option<SearchResult> -> SearchResult = { priority_type: LibraryPriority::Priority, .. }
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

    /// Get a list of substances found in permitted libraries
    pub fn get_permitted_found_substances(&self) -> Vec<String> {
        self.search_results
        .iter()//  HashMap<WhatIsFound, Option<SearchResult> -> SearchResult = { priority_type: LibraryPriority::Priority, .. }
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

    ///////////////////////////////////////THERMAL PROPERTIES////////////////////////////////////
    /// Extract coefficients for thermal polynoms for given substance
    pub fn extract_thermal_coeffs(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> Result<(), String> {
        // Get a mutable reference to the substance data
        let datamap = self.get_substance_result_mut(substance).unwrap();

        // Get a mutable reference to the Thermo entry
        let thermo = datamap
            .get_mut(&WhatIsFound::Thermo)
            .unwrap()
            .as_mut()
            .unwrap();

        // Get a mutable reference to the calculator
        let calculator = thermo.calculator.as_mut().unwrap();
        match calculator {
            CalculatorType::Thermo(thermo) => {
                let _ = thermo.extract_model_coefficients(temperature);
                Ok(())
            }
            CalculatorType::Transport(_) => {
                panic!("Substance found but has transport calculator instead of thermo");
            }
        }
    }
    pub fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> Result<(), String> {
        for substance in self.clone().search_results.keys() {
            match self.extract_thermal_coeffs(&substance, temperature) {
                Ok(_) => (),
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
    /// Calculate thermodynamic properties for a substance
    pub fn calculate_thermo_properties(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> Result<(f64, f64, f64), String> {
        //   println!("\n calculating substance {} with data {:?} at {} K \n", substance,  self.search_results.get_mut(substance).unwrap(), temperature);
        // println!( "\n calculating substance {} with data {:?} at {} K \n", substance, self.search_results.get(substance).unwrap(),temperature);
        match self
            .search_results
            .get(substance)
            .unwrap()
            .get(&WhatIsFound::Thermo)
            .unwrap()
        {
            //let's remember that search_results = Found{... calculator: Option<CalculatorType>}
            Some(SearchResult {
                calculator: Some(CalculatorType::Thermo(thermo)),
                ..
            }) => {
                let mut thermo = thermo.clone();
                //thermo.print_instance();

                if let Err(e) = thermo.calculate_Cp_dH_dS(temperature) {
                    return Err(format!("Failed to calculate properties: {}", e));
                }
                let Cp = thermo.get_Cp().unwrap_or(0.0);
                let dh = thermo.get_dh().unwrap_or(0.0);
                let ds = thermo.get_ds().unwrap_or(0.0);
                // dbg!("cp, dh, ds", Cp, dh, ds);
                Ok((Cp, dh, ds))
            }
            Some(SearchResult {
                calculator: Some(CalculatorType::Transport(_)),
                ..
            }) => Err("Substance found but has transport calculator instead of thermo".to_string()),
            Some(SearchResult {
                calculator: None, ..
            }) => Err("Substance found but has no calculator".to_string()),

            None => Err("Substance not in search results".to_string()),
        }
    }
    /// Set Molar mass and molar mass unit
    pub fn set_M(&mut self, M_map: HashMap<String, f64>, M_unit: Option<String>) {
        self.hasmap_of_molar_mass = M_map;
        self.Molar_mass_unit = M_unit;
    }
    /// Set pressure and pressure unit
    pub fn set_P(&mut self, P: f64, P_unit: Option<String>) {
        self.P = Some(P);
        self.P_unit = P_unit;
    }

    /// Calculate and populate therm_map_of_properties_values for all substances at a given temperature
    pub fn calculate_therm_map_of_properties(&mut self, temperature: f64) -> Result<(), String> {
        for substance in &self.substances.clone() {
            match self
                .calculate_thermo_properties(substance, temperature)
                .clone()
            {
                Ok((cp, dh, ds)) => {
                    let mut property_map = HashMap::new();
                    property_map.insert(DataType::Cp, Some(cp));
                    property_map.insert(DataType::dH, Some(dh));
                    property_map.insert(DataType::dS, Some(ds));
                    self.therm_map_of_properties_values
                        .insert(substance.clone(), property_map);
                }
                Err(e) => {
                    // If calculation fails, insert None for all properties
                    let mut property_map = HashMap::new();
                    property_map.insert(DataType::Cp, None);
                    property_map.insert(DataType::dH, None);
                    property_map.insert(DataType::dS, None);
                    self.therm_map_of_properties_values
                        .insert(substance.clone(), property_map);
                    println!(
                        "Warning: Failed to calculate properties for {}: {}",
                        substance, e
                    );
                }
            }
        }
        Ok(())
    }

    /// Calculate and populate therm_map_of_fun with function closures for all substances
    pub fn calculate_therm_map_of_fun(&mut self) -> Result<(), String> {
        for substance in &self.substances.clone() {
            match self
                .search_results
                .get_mut(substance)
                .expect("substance not found among search results")
                .get_mut(&WhatIsFound::Thermo)
                .unwrap()
            {
                Some(SearchResult {
                    calculator: Some(CalculatorType::Thermo(thermo)),
                    ..
                }) => {
                    let mut thermo = thermo.clone();

                    // Create closures for thermodynamic functions
                    if let Err(e) = thermo.create_closures_Cp_dH_dS() {
                        println!(
                            "Warning: Failed to create closures for {}: {}",
                            substance, e
                        );
                        continue;
                    }

                    // Get the closures and store them in the map
                    let mut function_map = HashMap::new();

                    match thermo.get_C_fun() {
                        Ok(cp_fun) => function_map.insert(DataType::Cp_fun, Some(cp_fun)),
                        Err(e) => {
                            println!(
                                "Warning: Failed to get Cp function for {}: {}",
                                substance, e
                            );
                            function_map.insert(DataType::Cp_fun, None)
                        }
                    };

                    match thermo.get_dh_fun() {
                        Ok(dh_fun) => function_map.insert(DataType::dH_fun, Some(dh_fun)),
                        Err(e) => {
                            println!(
                                "Warning: Failed to get dH function for {}: {}",
                                substance, e
                            );
                            function_map.insert(DataType::dH_fun, None)
                        }
                    };

                    match thermo.get_ds_fun() {
                        Ok(ds_fun) => function_map.insert(DataType::dS_fun, Some(ds_fun)),
                        Err(e) => {
                            println!(
                                "Warning: Failed to get dS function for {}: {}",
                                substance, e
                            );
                            function_map.insert(DataType::dS_fun, None)
                        }
                    };

                    self.therm_map_of_fun
                        .insert(substance.clone(), function_map);
                }
                _ => {
                    // If substance not found or doesn't have a thermo calculator, insert None for all functions
                    let mut function_map = HashMap::new();
                    function_map.insert(DataType::Cp_fun, None);
                    function_map.insert(DataType::dH_fun, None);
                    function_map.insert(DataType::dS_fun, None);
                    self.therm_map_of_fun
                        .insert(substance.clone(), function_map);
                }
            }
        }
        Ok(())
    }

    /// Calculate and populate therm_map_of_sym with symbolic expressions for all substances
    pub fn calculate_therm_map_of_sym(&mut self) -> Result<(), String> {
        for substance in &self.substances.clone() {
            match self
                .search_results
                .get_mut(substance)
                .unwrap()
                .get_mut(&WhatIsFound::Thermo)
                .unwrap()
            {
                Some(SearchResult {
                    calculator: Some(CalculatorType::Thermo(thermo)),
                    ..
                }) => {
                    let mut thermo = thermo.clone();

                    // Create symbolic expressions for thermodynamic functions
                    if let Err(e) = thermo.create_sym_Cp_dH_dS() {
                        println!(
                            "Warning: Failed to create symbolic expressions for {}: {}",
                            substance, e
                        );
                        continue;
                    }

                    // Get the symbolic expressions and store them in the map
                    let mut sym_map = HashMap::new();

                    match thermo.get_Cp_sym() {
                        Ok(cp_sym) => sym_map.insert(DataType::Cp_sym, Some(Box::new(cp_sym))),
                        Err(e) => {
                            println!(
                                "Warning: Failed to get Cp symbolic expression for {}: {}",
                                substance, e
                            );
                            sym_map.insert(DataType::Cp_sym, None)
                        }
                    };

                    match thermo.get_dh_sym() {
                        Ok(dh_sym) => sym_map.insert(DataType::dH_sym, Some(Box::new(dh_sym))),
                        Err(e) => {
                            println!(
                                "Warning: Failed to get dH symbolic expression for {}: {}",
                                substance, e
                            );
                            sym_map.insert(DataType::dH_sym, None)
                        }
                    };

                    match thermo.get_ds_sym() {
                        Ok(ds_sym) => sym_map.insert(DataType::dS_sym, Some(Box::new(ds_sym))),
                        Err(e) => {
                            println!(
                                "Warning: Failed to get dS symbolic expression for {}: {}",
                                substance, e
                            );
                            sym_map.insert(DataType::dS_sym, None)
                        }
                    };

                    self.therm_map_of_sym.insert(substance.clone(), sym_map);
                }
                _ => {
                    // If substance not found or doesn't have a thermo calculator, insert None for all symbolic expressions
                    let mut sym_map = HashMap::new();
                    sym_map.insert(DataType::Cp_sym, None);
                    sym_map.insert(DataType::dH_sym, None);
                    sym_map.insert(DataType::dS_sym, None);
                    self.therm_map_of_sym.insert(substance.clone(), sym_map);
                }
            }
        }
        Ok(())
    }

    /// Get a function for calculating a specific thermodynamic property for a substance
    pub fn get_thermo_function(
        &self,
        substance: &str,
        data_type: DataType,
    ) -> Option<&Box<dyn Fn(f64) -> f64>> {
        self.therm_map_of_fun
            .get(substance)
            .and_then(|map| map.get(&data_type))
            .and_then(|opt| opt.as_ref())
    }

    /// Get a symbolic expression for a specific thermodynamic property for a substance
    pub fn get_thermo_symbolic(&self, substance: &str, data_type: DataType) -> Option<&Box<Expr>> {
        self.therm_map_of_sym
            .get(substance)
            .and_then(|map| map.get(&data_type))
            .and_then(|opt| opt.as_ref())
    }
}
