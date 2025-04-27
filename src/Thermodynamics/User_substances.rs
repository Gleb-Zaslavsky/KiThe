/// agregator for all the thermo and transport calculations for user defined substances
use crate::Thermodynamics::DBhandlers::thermo_api::{
    ThermoCalculator, ThermoEnum, create_thermal_by_name,
};
use crate::Thermodynamics::DBhandlers::transport_api::{
    TransportCalculator, TransportEnum, create_transport_calculator_by_name,
};
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use std::{fmt, vec};

//use RustedSciThe::symbolic::symbolic_engine::Expr;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
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
    data: Value,
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
    pub Molar_mass: Option<HashMap<String, f64>>,
    pub Molar_mass_unit: Option<String>,
    /// Vector of substances to search for
    pub substances: Vec<String>,
    /// Maps libraries to their priority level
    pub library_priorities: HashMap<String, LibraryPriority>,
    ///
    pub explicit_search_insructions: HashMap<String, String>,
    ro_map: Option<HashMap<String, f64>>,
    ///
    ro_map_sym: Option<HashMap<String, Box<Expr>>>,
    /// hashmap of phases of substances
    pub map_of_phases: HashMap<String, Option<Phases>>,
    /// Stores search results for each substance
    pub search_results: HashMap<String, HashMap<WhatIsFound, Option<SearchResult>>>,
    /// Reference to the thermo data library
    thermo_data: ThermoData,
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
            Molar_mass: self.Molar_mass.clone(),
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
            Molar_mass: None,
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
    fn create_calculator(&self, library: &str) -> Option<CalculatorType> {
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
        self.Molar_mass = Some(M_map);
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

    /// Calculate 1) matrix of element composition 2) molar mass for all substances

    pub fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> Result<(), String> {
        let (hasmap_of_molar_mass, vec_of_compositions, hashset_of_elems) =
            Self::calculate_elem_composition_and_molar_mass_local(self, groups).unwrap();
        // vector of chemical elements names (unique)
        let mut unique_vec_of_elems = hashset_of_elems.into_iter().collect::<Vec<_>>();
        unique_vec_of_elems.sort();
        let num_rows = unique_vec_of_elems.len();
        let num_cols = vec_of_compositions.len();
        // allocate matrix with num of rows = num of elements and num of cols = num of substances
        let mut matrix: DMatrix<f64> = DMatrix::zeros(num_rows, num_cols);
        for substance_i in 0..self.substances.len() {
            for j in 0..unique_vec_of_elems.len() {
                let element_j = unique_vec_of_elems[j].clone();
                if let Some(count) = vec_of_compositions[substance_i].get(&element_j) {
                    matrix[(j, substance_i)] += *count as f64;
                }
            }
        }
        self.elem_composition_matrix = Some(matrix.transpose());
        self.hasmap_of_molar_mass = hasmap_of_molar_mass;
        self.unique_elements = unique_vec_of_elems;
        Ok(())
    }

    pub fn calculate_elem_composition_and_molar_mass_local(
        sd: &mut SubsData,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> Result<
        (
            HashMap<String, f64>,
            Vec<HashMap<String, f64>>,
            HashSet<String>,
        ),
        String,
    > {
        use crate::Kinetics::molmass::{
            calculate_molar_mass, calculate_molar_mass_for_composition,
        };
        let mut vec_of_compositions: Vec<HashMap<String, f64>> = Vec::new();
        let mut hashset_of_elems: HashSet<String> = HashSet::new();
        let mut hasmap_of_molar_mass: HashMap<String, f64> = HashMap::new();

        for substance in &sd.substances.clone() {
            let (composition, molar_masss) = match sd.search_results.get_mut(substance) {
                Some(result_map) => match result_map.get_mut(&WhatIsFound::Thermo) {
                    Some(Some(SearchResult {
                        calculator: Some(CalculatorType::Thermo(thermo)),
                        ..
                    })) => {
                        let thermo = thermo.clone();

                        if let Some(composition) = thermo.get_composition().unwrap() {
                            // if substance record has composition field then
                            println!(
                                "in the library record for substance {}  found composition: {:#?}",
                                substance, composition
                            );
                            // use it to calculate molar mass and element composition matrix
                            let composition_usize: HashMap<String, usize> = composition
                                .iter()
                                .map(|(k, v)| (k.clone(), *v as usize))
                                .collect();
                            let molar_masss =
                                calculate_molar_mass_for_composition(composition_usize.clone());

                            (composition, molar_masss)
                        } else {
                            // if not then calculate molar mass and element composition matrix from the substance formula
                            let (molar_masss, composition) =
                                calculate_molar_mass(substance.clone(), groups.clone());

                            let composition: HashMap<String, f64> = composition
                                .iter()
                                .map(|(k, v)| (k.clone(), *v as f64))
                                .collect();

                            (composition, molar_masss)
                        }
                    }
                    // Handle the case where Thermo entry exists but doesn't have a calculator
                    _ => {
                        // Calculate from formula since we don't have thermo data
                        let (molar_masss, composition) =
                            calculate_molar_mass(substance.clone(), groups.clone());

                        let composition: HashMap<String, f64> = composition
                            .iter()
                            .map(|(k, v)| (k.clone(), *v as f64))
                            .collect();

                        (composition, molar_masss)
                    }
                },
                // Handle the case where the substance doesn't have any search results
                None => {
                    println!("No search results found for substance: {}", substance);
                    // Calculate from formula since we don't have any data
                    let (molar_masss, composition) =
                        calculate_molar_mass(substance.clone(), groups.clone());

                    let composition: HashMap<String, f64> = composition
                        .iter()
                        .map(|(k, v)| (k.clone(), *v as f64))
                        .collect();

                    (composition, molar_masss)
                }
            }; // let (composition, molar_masss)
            hasmap_of_molar_mass.insert(substance.clone(), molar_masss);
            vec_of_compositions.push(composition.clone());
            let elements = composition.keys().map(|el| el.clone()).collect::<Vec<_>>();
            hashset_of_elems.extend(elements);
        }
        Ok((hasmap_of_molar_mass, vec_of_compositions, hashset_of_elems))
    }
    /////////////////////////////TRANSPORT COEFFICIENTS////////////////////////////////////
    /// Extract transport coefficients for a substance
    pub fn extract_transport_coeffs(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> Result<(), String> {
        // Get a mutable reference to the substance data
        let datamap = self.get_substance_result_mut(substance).unwrap();

        // Get a mutable reference to the Thermo entry
        let transport = datamap
            .get_mut(&WhatIsFound::Transport)
            .unwrap()
            .as_mut()
            .unwrap();

        // Get a mutable reference to the calculator
        let calculator = transport.calculator.as_mut().unwrap();
        match calculator {
            CalculatorType::Thermo(_) => {
                panic!("Substance found but has thermo calculator instead of transport");
            }
            CalculatorType::Transport(transport) => {
                let _ = transport.extract_coefficients(temperature);
                Ok(())
            }
        }
    }
    pub fn extract_all_transport_coeffs(&mut self, temperature: f64) -> Result<(), String> {
        for substance in self.clone().search_results.keys() {
            match self.extract_transport_coeffs(&substance, temperature) {
                Ok(_) => (),
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
    /// Calculate transport properties for a substance
    pub fn calculate_transport_properties(
        &mut self,
        substance: &str,
        temperature: f64,
        Cp: Option<f64>,
        ro: Option<f64>,
    ) -> Result<(f64, f64), String> {
        // println!(  "\n \n found  for substance {} {:?}", substance, &self.search_results.get(substance).unwrap().get(&WhatIsFound::Transport).unwrap());
        match self
            .search_results
            .get_mut(substance)
            .unwrap()
            .get_mut(&WhatIsFound::Transport)
            .unwrap()
        {
            Some(SearchResult {
                calculator: Some(CalculatorType::Transport(transport)),
                ..
            }) => {
                let mut transport = transport.clone();
                println!("transport instance: {:?}", transport);

                //  transport.print_instance();

                let P = self.P.expect("Pressure not set");
                let P_unit = self.P_unit.clone();
                let M_map = self.Molar_mass.clone().expect("Molar mass not set");
                let M = M_map
                    .get(substance)
                    .expect("Molar mass for this substance not found");
                let M_unit = self.Molar_mass_unit.clone();
                let _ = transport.set_M(*M, M_unit);
                let _ = transport.set_P(P, P_unit);
                let lambda = transport
                    .calculate_lambda(Cp, ro, temperature)
                    .map_err(|e| format!("Failed to calculate lambda: {}", e))?;

                let viscosity = transport
                    .calculate_viscosity(temperature)
                    .map_err(|e| format!("Failed to calculate viscosity: {}", e))?;

                Ok((lambda, viscosity))
            }
            Some(SearchResult {
                calculator: Some(CalculatorType::Thermo(_)),
                ..
            }) => Err("Substance found but has thermo calculator instead of transport".to_string()),
            Some(SearchResult {
                calculator: None, ..
            }) => Err("Substance found but has no calculator".to_string()),

            None => Err("Substance not in search results".to_string()),
        }
    }
    /// Calculate and populate transport_map_of_properties_values for all substances at a given temperature
    pub fn calculate_transport_map_of_properties(
        &mut self,
        temperature: f64,
    ) -> Result<(), String> {
        for substance in &self.substances.clone() {
            if let Some(Cp) = self
                .therm_map_of_properties_values
                .get(substance)
                .unwrap()
                .get(&DataType::Cp)
                .unwrap()
                .clone()
            {
                let ro = if let Some(ro_map) = self.ro_map.clone() {
                    ro_map.get(substance).cloned()
                } else {
                    None
                };
                match self
                    .calculate_transport_properties(substance, temperature, Some(Cp), ro)
                    .clone()
                {
                    Ok((Lambda, Visc)) => {
                        let mut property_map = HashMap::new();
                        property_map.insert(DataType::Lambda, Some(Lambda));
                        property_map.insert(DataType::Visc, Some(Visc));

                        self.transport_map_of_properties_values
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
                } // match 
            } else {
                break;
                // return Err("Cp not found".to_string());
            }
        }
        Ok(())
    }

    /// Calculate and populate transport_map_of_sym with symbolic expressions for all substances
    pub fn calculate_transport_map_of_sym(&mut self) -> Result<(), String> {
        for substance in &self.substances.clone() {
            if let Some(Cp) = self
                .therm_map_of_sym
                .get(substance)
                .unwrap()
                .get(&DataType::Cp_sym)
                .unwrap()
                .clone()
            {
                let ro = if let Some(ro_map) = self.ro_map_sym.clone() {
                    Some(*ro_map.get(substance).unwrap().clone())
                } else {
                    None
                };
                match self
                    .search_results
                    .get_mut(substance)
                    .unwrap()
                    .get_mut(&WhatIsFound::Transport)
                    .unwrap()
                {
                    Some(SearchResult {
                        calculator: Some(CalculatorType::Transport(transport)),
                        ..
                    }) => {
                        let mut transport = transport.clone();

                        // Create symbolic expressions for thermodynamic functions
                        if let Err(e) = transport.create_symbolic_lambda(Some(*Cp), ro) {
                            println!(
                                "Warning: Failed to create symbolic expressions for {}: {}",
                                substance, e
                            );
                            continue;
                        }

                        // Get the symbolic expressions and store them in the map
                        let mut sym_map = HashMap::new();

                        match transport.get_lambda_sym() {
                            Ok(Lambda_sym) => {
                                sym_map.insert(DataType::Lambda_sym, Some(Box::new(Lambda_sym)))
                            }
                            Err(e) => {
                                println!(
                                    "Warning: Failed to get Lambda symbolic expression for {}: {}",
                                    substance, e
                                );
                                sym_map.insert(DataType::Lambda_sym, None)
                            }
                        };
                        if let Err(e) = transport.create_symbolic_viscosity() {
                            println!(
                                "Warning: Failed to create symbolic expressions for {}: {}",
                                substance, e
                            );
                            continue;
                        }
                        match transport.get_viscosity_sym() {
                            Ok(visc_sym) => {
                                sym_map.insert(DataType::Visc_sym, Some(Box::new(visc_sym)))
                            }
                            Err(e) => {
                                println!(
                                    "Warning: Failed to get dH symbolic expression for {}: {}",
                                    substance, e
                                );
                                sym_map.insert(DataType::Visc_sym, None)
                            }
                        };

                        self.transport_map_of_sym.insert(substance.clone(), sym_map);
                    }
                    _ => {
                        // If substance not found or doesn't have a thermo calculator, insert None for all symbolic expressions
                        let mut sym_map = HashMap::new();
                        sym_map.insert(DataType::Lambda_sym, None);
                        sym_map.insert(DataType::Visc_sym, None);

                        self.therm_map_of_sym.insert(substance.clone(), sym_map);
                    }
                }
            }
        }
        Ok(())
    }

    pub fn calculate_transport_map_of_functions(&mut self) -> Result<(), String> {
        for substance in &self.substances.clone() {
            if let Some(thermo_info) = self.therm_map_of_properties_values.get(substance) {
                if let Some(Cp) = thermo_info.get(&DataType::Cp).unwrap().clone() {
                    let ro = if let Some(ro_map) = self.ro_map.clone() {
                        Some(ro_map.get(substance).unwrap().clone())
                    } else {
                        None
                    };
                    if self
                        .search_results
                        .get_mut(substance)
                        .unwrap()
                        .get_mut(&WhatIsFound::Transport)
                        .is_none()
                    {
                        continue;
                    }
                    match self
                        .search_results
                        .get_mut(substance)
                        .unwrap()
                        .get_mut(&WhatIsFound::Transport)
                        .unwrap()
                    {
                        Some(SearchResult {
                            calculator: Some(CalculatorType::Transport(transport)),
                            ..
                        }) => {
                            let mut transport = transport.clone();

                            // Create symbolic expressions for thermodynamic functions
                            if let Err(e) = transport.create_lambda_closure(Some(Cp), ro) {
                                println!(
                                    "Warning: Failed to create symbolic expressions for {}: {}",
                                    substance, e
                                );
                                continue;
                            }

                            // Get the symbolic expressions and store them in the map
                            let mut fun_map = HashMap::new();

                            match transport.get_lambda_fun() {
                                Ok(Lambda_fun) => {
                                    fun_map.insert(DataType::Lambda_fun, Some(Lambda_fun))
                                }
                                Err(e) => {
                                    println!(
                                        "Warning: Failed to get Lambda symbolic expression for {}: {}",
                                        substance, e
                                    );
                                    fun_map.insert(DataType::Lambda_fun, None)
                                }
                            };
                            if let Err(e) = transport.create_viscosity_closure() {
                                println!(
                                    "Warning: Failed to create symbolic expressions for {}: {}",
                                    substance, e
                                );
                                continue;
                            }
                            match transport.get_viscosity_fun() {
                                Ok(visc_fun) => fun_map.insert(DataType::Visc_fun, Some(visc_fun)),
                                Err(e) => {
                                    println!(
                                        "Warning: Failed to get dH symbolic expression for {}: {}",
                                        substance, e
                                    );
                                    fun_map.insert(DataType::Visc_fun, None)
                                }
                            };

                            self.transport_map_of_fun.insert(substance.clone(), fun_map);
                        }
                        _ => {
                            // If substance not found or doesn't have a thermo calculator, insert None for all symbolic expressions
                            let mut fun_map = HashMap::new();
                            fun_map.insert(DataType::Lambda_fun, None);
                            fun_map.insert(DataType::Visc_fun, None);

                            self.therm_map_of_sym.insert(substance.clone(), fun_map);
                        }
                    }
                }
            }
        }
        Ok(())
    }
    ///////////////////////////////INPUT/OUTPUT///////////////////////////////////////////
    /// Print a summary of the search results with calculator types
    pub fn print_search_summary(&self) {
        println!("\nSearch Results Summary:");
        println!("----------------------");

        for (substance, map) in &self.search_results.clone() {
            for (what_is_found, result) in map.iter() {
                match what_is_found {
                    WhatIsFound::NotFound => {
                        println!("{}: Not Found", substance);
                    }
                    _ => match result.clone().unwrap() {
                        SearchResult {
                            library,
                            priority_type,
                            calculator,
                            ..
                        } => {
                            let priority_str = match priority_type {
                                LibraryPriority::Priority => "Priority",
                                LibraryPriority::Permitted => "Permitted",
                                LibraryPriority::Explicit => "Explicit",
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
                    },
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
        println!("Not found: {} \n ", self.get_not_found_substances().len());
        println!("--------End of Search Results Summary------------------ \n \n");
    }
    /// if substance is not found in local libraries, go to NIST
    pub fn if_not_found_go_NIST(&mut self) -> Result<(), String> {
        use crate::Thermodynamics::DBhandlers::NIST_parser::{Phase, SearchType};
        let list_of_thermo_libs = vec![
            "NASA",
            "NIST",
            "NASA7",
            "NUIG_thermo",
            "Cantera_nasa_base_gas",
            "Cantera_nasa_base_cond",
            "NASA_gas",
            "NASA_cond",
        ];
        if self.get_not_found_substances().len() > 0 {
            println!(
                "No data found for {} substances",
                self.get_not_found_substances().len()
            );
            println!("\nTrying to get data from NIST...\n");
            // Get a list of substances not found to avoid borrowing issues
            let not_found_substances: Vec<String> = self.get_not_found_substances();
            for substance in not_found_substances {
                println!(
                    "for substance {} which has not been found in local libraries",
                    &substance
                );
                if list_of_thermo_libs // and this data must be thermodynamic
                    .iter()
                    .any(|&x| self.library_priorities.contains_key(x))
                {
                    let mut calculator = self.create_calculator("NIST");
                    if let Some(calc_type) = calculator.as_mut() {
                        if let CalculatorType::Thermo(thermo) = calc_type {
                            let phase = if let Some(Some(phase)) =
                                self.map_of_phases.get(&substance)
                            {
                                match phase {
                                    Phases::Gas => Phase::Gas,
                                    Phases::Solid => Phase::Solid,
                                    Phases::Liquid => Phase::Liquid,
                                    _ => {
                                        println!(
                                            "Unknown phase for substance: {}, defaulting to Gas",
                                            substance
                                        );
                                        Phase::Gas
                                    }
                                }
                            } else {
                                Phase::Gas // Default to Gas if no phase specified
                            };
                            // Try to get data from NIST
                            if let Ok(()) =
                                thermo.renew_base(substance.clone(), SearchType::All, phase)
                            {
                                println!("\nFound data in NIST for {}", substance);
                                // Create a new SearchResult with the NIST data
                                let search_result = SearchResult {
                                    library: "NIST".to_string(),
                                    priority_type: LibraryPriority::Permitted, // or Priority based on your preference
                                    data: serde_json::Value::Null, // This should be the actual data from NIST
                                    calculator: Some(CalculatorType::Thermo(thermo.clone())),
                                }; // search_result
                                // Update the search_results HashMap
                                let mut result_map = HashMap::new();
                                result_map.insert(WhatIsFound::Thermo, Some(search_result));
                                // Remove the NotFound entry and add the new Thermo entry
                                self.search_results.insert(substance.clone(), result_map);
                                println!("Updated search results with NIST data for {}", substance);
                            }
                            // if let Ok(()) = thermo.renew_base(substance.clone(), SearchType::All, phase) {
                            else {
                                println!("Failed to find data in NIST for {}", substance);
                            }
                        } // if let CalculatorType::Thermo(calc
                    } // if let Some(calculator)
                } // list_of_thermo_libs...
            } // for substance..
            Ok(())
        }
        // if substances not found in priority libraries, try to get data from NIST
        else {
            Ok(())
        }
    }
}
