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
    library: String,
    priority_type: LibraryPriority,
    data: Value,
    pub calculator: Option<CalculatorType>,
}
///LibraryPriority enum: Defines priority levels (Priority or Permitted)
#[derive(Debug, Clone, PartialEq)]
pub enum LibraryPriority {
    Priority,
    Permitted,
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
    library_priorities: HashMap<String, LibraryPriority>,
    /// Stores search results for each substance
    pub search_results: HashMap<String, HashMap<WhatIsFound, Option<SearchResult>>>,
    /// Reference to the thermo data library
    thermo_data: ThermoData,
    ///
    pub therm_map_of_properties_values: HashMap<String, HashMap<DataType, Option<f64>>>,
    ///
    pub therm_map_of_fun: HashMap<String, HashMap<DataType, Option<Box<dyn Fn(f64) -> f64>>>>,
    ///
    pub therm_map_of_sym: HashMap<String, HashMap<DataType, Option<Box<Expr>>>>,
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
        Self {
            P: self.P,
            P_unit: self.P_unit.clone(),
            Molar_mass: self.Molar_mass.clone(),
            Molar_mass_unit: self.Molar_mass_unit.clone(),
            substances: self.substances.clone(),
            library_priorities: self.library_priorities.clone(),
            search_results: self.search_results.clone(),
            thermo_data: self.thermo_data.clone(),
            therm_map_of_properties_values: self.therm_map_of_properties_values.clone(),
            therm_map_of_fun: new_therm_map_of_fun,
            therm_map_of_sym: self.therm_map_of_sym.clone(),
        }
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
            search_results: HashMap::new(),
            thermo_data: ThermoData::new(),
            therm_map_of_properties_values: HashMap::new(),
            therm_map_of_fun: HashMap::new(),
            therm_map_of_sym: HashMap::new(),
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
        // Search for each substance
        for substance in &self.substances {
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
        }
    }

    /// Get the search result for a specific substance
    pub fn get_substance_result(
        &self,
        substance: &str,
    ) -> Option<&HashMap<WhatIsFound, Option<SearchResult>>> {
        //Option<&SearchResult> {
        self.search_results.get(substance)
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
        self.search_results
            .iter()
            .filter(|(_, result)| result.values().any(|result| matches!(result.as_ref().unwrap(), SearchResult { priority_type, .. } if *priority_type == LibraryPriority::Priority)))
            .map(|(substance, _)| substance.clone())
            .collect()
    }

    /// Get a list of substances found in permitted libraries
    pub fn get_permitted_found_substances(&self) -> Vec<String> {
        self.search_results
        .iter()
        .filter(|(_, result)| result.values().any(|result| matches!(result.as_ref().unwrap(), SearchResult { priority_type, .. } if *priority_type == LibraryPriority::Permitted)))
        .map(|(substance, _)| substance.clone())
        .collect()
    }
    pub fn extract_thermal_coeffs(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> Result<(), String> {
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
                //thermo.print_instance();
                let mut thermo = thermo.clone();
                // Extract coefficients and calculate properties
                if let Err(e) = thermo.extract_model_coefficients(temperature) {
                    return Err(format!("Failed to extract coefficients: {}", e));
                }
                Ok(())
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

                // Extract coefficients and calculate properties
                if let Err(e) = thermo.extract_model_coefficients(temperature) {
                    return Err(format!("Failed to extract coefficients: {}", e));
                }

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

    pub fn set_M(&mut self, M_map: HashMap<String, f64>, M_unit: Option<String>) {
        self.Molar_mass = Some(M_map);
        self.Molar_mass_unit = M_unit;
    }
    pub fn set_P(&mut self, P: f64, P_unit: Option<String>) {
        self.P = Some(P);
        self.P_unit = P_unit;
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
                // Extract coefficients and calculate properties
                if let Err(e) = transport.extract_coefficients(temperature) {
                    return Err(format!("Failed to extract coefficients: {}", e));
                }
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
    // TODO:test that function
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
                .unwrap()
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_user_substances_with_calculators() {
        // Create a new SubsData instance with test substances
        let substances = vec!["CO".to_string(), "H2O".to_string()];
        let mut user_subs = SubsData::new();
        user_subs.substances = substances;

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Perform the search
        user_subs.search_substances();

        // Test thermo calculations
        if let Ok((cp, dh, ds)) = user_subs.calculate_thermo_properties("CO", 400.0) {
            println!("CO Thermo properties at 400K:");
            println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
            assert!(cp > 0.0);
        } else {
            panic!("Failed to calculate CO properties");
        }
        let Cp = 33.8;
        // Test transport calculations
        user_subs.set_M(
            HashMap::from([("H2O".to_string(), 18.0), ("CO".to_string(), 32.0)]),
            None,
        );
        user_subs.set_P(1e5, None);
        if let Ok((lambda, viscosity)) =
            user_subs.calculate_transport_properties("H2O", 400.0, Some(Cp), None)
        {
            println!("CO Transport properties at 400K:");
            println!("Lambda: {}, Viscosity: {}", lambda, viscosity);
            assert!(lambda > 0.0);
            assert!(viscosity > 0.0);
        } else {
            panic!("Failed to calculate H2O properties");
        }

        // Print full summary
        user_subs.print_search_summary();
    }
    #[test]
    fn test_user_substances_with_calculators_own_functions() {
        // Create a new SubsData instance with test substances
        let substances = vec!["CO".to_string(), "H2O".to_string()];
        let mut user_subs = SubsData::new();
        user_subs.substances = substances;

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Perform the search
        user_subs.search_substances();
        // Print full summary
        user_subs.print_search_summary();
        let datamap = user_subs.get_substance_result("CO").unwrap();
        let Thermo = datamap.get(&WhatIsFound::Thermo).unwrap().as_ref().unwrap();
        let Calculator = Thermo.calculator.as_ref().unwrap();
        let Cp;
        match Calculator {
            CalculatorType::Thermo(thermo) => {
                // Test thermo calculations
                let mut thermo = thermo.clone();
                thermo.extract_model_coefficients(400.0).unwrap();
                if let Ok(()) = thermo.calculate_Cp_dH_dS(400.0) {
                    let (cp, dh, ds) = (
                        thermo.get_Cp().unwrap(),
                        thermo.get_dh().unwrap(),
                        thermo.get_ds().unwrap(),
                    ); //thermo.get_Cp().unwrap();
                    println!("CO Thermo properties at 400K:");
                    println!("Cp: {}, dH: {}, dS: {}", cp, dh, ds);
                    Cp = cp;
                    assert!(cp > 0.0);
                } else {
                    panic!("Failed to calculate CO properties");
                }
            }
            _ => {
                panic!("Failed to calculate CO properties");
            }
        }
        let Transport = datamap
            .get(&WhatIsFound::Transport)
            .unwrap()
            .as_ref()
            .unwrap();
        let Calculator = Transport.calculator.as_ref().unwrap();
        match Calculator {
            CalculatorType::Transport(transport) => {
                // Test transport calculations
                let mut transport = transport.clone();
                let _ = transport.set_M(32.0, None);
                let _ = transport.set_P(1e5, None);
                transport.extract_coefficients(400.0).unwrap();
                if let Ok(L) = transport.calculate_lambda(Some(Cp), None, 400.0) {
                    println!("CO Transport properties at 400K:");
                    println!("Lambda: {}", L);
                    assert!(L > 0.0);
                } else {
                    panic!("Failed to calculate H2O properties");
                }
            }
            _ => {
                panic!("Failed to calculate CO properties");
            }
        }
    }
    #[test]
    fn test_library_priority_setting() {
        let mut user_subs = SubsData::new();

        // Test setting single library priority
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        assert_eq!(
            user_subs.library_priorities.get("NASA_gas"),
            Some(&LibraryPriority::Priority)
        );

        // Test setting multiple library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NIST".to_string(), "CEA".to_string()],
            LibraryPriority::Permitted,
        );
        assert_eq!(
            user_subs.library_priorities.get("NIST"),
            Some(&LibraryPriority::Permitted)
        );
        assert_eq!(
            user_subs.library_priorities.get("CEA"),
            Some(&LibraryPriority::Permitted)
        );

        // Test overriding priority
        user_subs.set_library_priority("NIST".to_string(), LibraryPriority::Priority);
        assert_eq!(
            user_subs.library_priorities.get("NIST"),
            Some(&LibraryPriority::Priority)
        );
    }

    #[test]
    fn test_substance_search() {
        // fails
        let mut user_subs = SubsData::new();
        user_subs.substances = vec![
            "H2O".to_string(),
            "CO2".to_string(),
            "NonExistentSubstance".to_string(),
        ];

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Permitted,
        );

        // Perform search
        user_subs.search_substances();

        // Check results
        assert!(matches!(
            user_subs
                .get_substance_result("H2O")
                .unwrap()
                .get(&WhatIsFound::Thermo)
                .unwrap(),
            Some(SearchResult { .. })
        ));
        let _data_to_expect: HashMap<WhatIsFound, Option<SearchResult>> =
            HashMap::from([(WhatIsFound::NotFound, None)]);
        assert!(matches!(
            user_subs.get_substance_result("NonExistentSubstance"),
            Some(_data_to_expect)
        ));

        // Check substance lists
        let not_found = user_subs.get_not_found_substances();
        assert!(not_found.contains(&"NonExistentSubstance".to_string()));

        let priority_found = user_subs.get_priority_found_substances();
        let permitted_found = user_subs.get_permitted_found_substances();
        println!(
            "Priority found: {:?}, Permitted found: {:?}",
            priority_found, permitted_found
        );
        // At least one of these should contain H2O
        assert!(
            priority_found.contains(&"H2O".to_string())
                || permitted_found.contains(&"H2O".to_string())
        );
    }

    #[test]
    fn test_thermo_property_calculation() {
        // fails
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["O2".to_string(), "CO".to_string()];

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();

        // Test for O2
        if let Ok((cp, dh, ds)) = user_subs.calculate_thermo_properties("O2", 400.0) {
            println!(
                "O2 Thermo properties at 400K: Cp={}, dH={}, dS={}",
                cp, dh, ds
            );
            assert!(cp > 0.0);
            assert!(ds > 0.0);
        } else {
            panic!("Failed to calculate O2 properties");
        }

        // Test for CO
        if let Ok((cp, dh, ds)) = user_subs.calculate_thermo_properties("CO", 400.0) {
            println!(
                "CO Thermo properties at 400K: Cp={}, dH={}, dS={}",
                cp, dh, ds
            );
            assert!(cp > 0.0);
            assert!(ds > 0.0);
        } else {
            panic!("Failed to calculate CO properties");
        }
    }

    #[test]
    fn test_transport_property_calculation() {
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["H2O".to_string(), "CO".to_string()];

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["Aramco_transport".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();
        user_subs.set_M(
            HashMap::from([("H2O".to_string(), 18.0), ("CO".to_string(), 32.0)]),
            None,
        );
        user_subs.set_P(1e5, None);
        // Test for substances that should have transport data
        for substance in &["H2O", "CO"] {
            if let Ok((lambda, viscosity)) =
                user_subs.calculate_transport_properties(substance, 400.0, Some(33.0), None)
            {
                println!(
                    "{} Transport properties at 400K: Lambda={}, Viscosity={}",
                    substance, lambda, viscosity
                );
                assert!(lambda > 0.0);
                assert!(viscosity > 0.0);
            }
        }
    }

    #[test]
    fn test_property_maps_calculation() {
        // fails
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["H2O".to_string(), "CO".to_string()];

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();

        // Calculate property values at 400K
        user_subs
            .calculate_therm_map_of_properties(500.0)
            .expect("Failed to calculate property map");

        // Check that values were stored
        for substance in &["H2O", "CO"] {
            let property_map = user_subs.therm_map_of_properties_values.get(*substance);
            assert!(property_map.is_some());

            let property_map = property_map.unwrap();
            assert!(property_map.get(&DataType::Cp).unwrap().is_some());
            assert!(property_map.get(&DataType::dH).unwrap().is_some());
            assert!(property_map.get(&DataType::dS).unwrap().is_some());

            // Check that values are reasonable
            assert!(property_map.get(&DataType::Cp).unwrap().unwrap() > 0.0);
        }
    }

    #[test]
    fn test_function_maps_calculation() {
        // fails
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["H2O".to_string()];

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();

        // Calculate function closures
        user_subs
            .calculate_therm_map_of_fun()
            .expect("Failed to calculate function map");

        // Check that functions were stored and work
        let cp_fun = user_subs.get_thermo_function("H2O", DataType::Cp_fun);
        assert!(cp_fun.is_some());

        // Test the function at different temperatures
        let cp_fun = cp_fun.unwrap();
        let cp_300 = cp_fun(300.0);
        let cp_400 = cp_fun(400.0);
        let cp_500 = cp_fun(500.0);

        println!(
            "H2O Cp values: 300K={}, 400K={}, 500K={}",
            cp_300, cp_400, cp_500
        );

        // Cp should be positive and should change with temperature
        assert!(cp_300 > 0.0);
        assert!(cp_400 > 0.0);
        assert!(cp_500 > 0.0);
        assert!(cp_300 != cp_400); // Values should differ with temperature
    }

    #[test]
    fn test_symbolic_maps_calculation() {
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["H2O".to_string()];

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );

        // Perform search
        user_subs.search_substances();

        // Calculate symbolic expressions
        user_subs
            .calculate_therm_map_of_sym()
            .expect("Failed to calculate symbolic map");

        // Check that symbolic expressions were stored
        let cp_sym = user_subs.get_thermo_symbolic("H2O", DataType::Cp_sym);
        assert!(cp_sym.is_some());

        // Check that the expression is not empty
        let cp_sym = cp_sym.unwrap();
        println!("H2O Cp symbolic expression: {}", cp_sym.to_string());
        assert!(!cp_sym.to_string().is_empty());
    }

    #[test]
    fn test_integration_all_features() {
        let mut user_subs = SubsData::new();
        user_subs.substances = vec!["H2O".to_string(), "CO".to_string(), "CH4".to_string()];

        // Set library priorities
        user_subs.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        user_subs.set_multiple_library_priorities(
            vec!["NIST".to_string(), "Aramco_transport".to_string()],
            LibraryPriority::Permitted,
        );

        // Perform search
        user_subs.search_substances();

        // Calculate all types of data
        user_subs
            .calculate_therm_map_of_properties(400.0)
            .expect("Failed to calculate property map");
        user_subs
            .calculate_therm_map_of_fun()
            .expect("Failed to calculate function map");
        user_subs
            .calculate_therm_map_of_sym()
            .expect("Failed to calculate symbolic map");

        // Print search summary
        user_subs.print_search_summary();

        // Verify that we have data for all substances
        for substance in &["H2O", "CO", "CH4"] {
            // Check if found in any library
            let result = user_subs.get_substance_result(substance);
            assert!(result.is_some());

            // If found, check that we have property values
            if let Some(_map) = result {
                let property_map = user_subs.therm_map_of_properties_values.get(*substance);
                assert!(property_map.is_some());

                // Check that we have at least one function
                let function_map = user_subs.therm_map_of_fun.get(*substance);
                assert!(function_map.is_some());

                // Check that we have at least one symbolic expression
                let sym_map = user_subs.therm_map_of_sym.get(*substance);
                assert!(sym_map.is_some());
            }
        }
    }
}
