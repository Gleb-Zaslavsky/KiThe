use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
use crate::Thermodynamics::DBhandlers::transport_api::TransportCalculator;
use crate::Thermodynamics::User_substances::{
    CalculatorType, DataType, LibraryPriority, Phases, SearchResult, SubsData, WhatIsFound,
};

use std::vec;

//use RustedSciThe::symbolic::symbolic_engine::Expr;

use nalgebra::DMatrix;

use std::collections::{HashMap, HashSet};
impl SubsData {
    ////////////////////////////////ELEMENT COMPOSITION AND MOLAR MASS////////////////////////////////////
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
    /// Molar mass and Pressure must be set before calling this function
    /// calculate_elem_composition_and_molar_mass must be called before this function
    /// calculate_therm_map_of_properties_values must be called before this function (to get Cp values)
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
                library: lib_name,
                calculator: Some(CalculatorType::Transport(transport)),
                ..
            }) => {
                let mut transport = transport.clone();
                println!("transport instance: {:?}", transport);

                //  transport.print_instance();

                if lib_name.clone() == "CEA" {
                    let lambda = transport
                        .calculate_lambda(Cp, ro, temperature)
                        .map_err(|e| format!("Failed to calculate lambda: {}", e))?;

                    let viscosity = transport
                        .calculate_viscosity(temperature)
                        .map_err(|e| format!("Failed to calculate viscosity: {}", e))?;
                    return Ok((lambda, viscosity));
                } else {
                    let P = self.P.expect("Pressure not set");
                    let P_unit = self.P_unit.clone();
                    let M_map = self.hasmap_of_molar_mass.clone();
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
                    return Ok((lambda, viscosity));
                }
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
    /// Requires Cp values from therm_map_of_properties_values (so to use this function,
    /// first calculate therm_map_of_properties_values, call calculate_therm_map_of_properties_values)
    ///  calculate_elem_composition_and_molar_mass must be called before this function
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
