//! # Extended Substance Data Management Module
//!
//! This module extends the core `SubsData` functionality with advanced capabilities for:
//! - Element composition and molar mass calculations
//! - Transport property calculations (viscosity, thermal conductivity)
//! - Binary diffusion coefficient calculations
//! - NIST database integration for missing substances
//!
//! ## Key Features
//! - **Element Analysis**: Automatic parsing of chemical formulas and matrix generation
//! - **Transport Properties**: Temperature-dependent viscosity and thermal conductivity
//! - **Diffusion Calculations**: Binary diffusion coefficients for all substance pairs
//! - **Symbolic Expressions**: Generate symbolic forms of all properties
//! - **Function Closures**: Create temperature-dependent property functions
//! - **NIST Integration**: Automatic fallback to NIST database for missing data
//!
//! ## Usage Pattern
//! ```rust,ignore
//! let mut subs_data = SubsData::new();
//! subs_data.calculate_elem_composition_and_molar_mass(None)?;
//! subs_data.extract_all_transport_coeffs(300.0)?;
//! subs_data.initialize_diffusion_data()?;
//! subs_data.calculate_all_pairs()?;
//! ```

use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
use crate::Thermodynamics::DBhandlers::transport_api::{
    TransportCalculator, TransportEnum, TransportError,
};
use crate::Thermodynamics::User_substances::{
    CalculatorType, DataType, LibraryPriority, Phases, SearchResult, SubsData, WhatIsFound,
};

use crate::Thermodynamics::DBhandlers::Diffusion::MultiSubstanceDiffusion;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::vec;
//use RustedSciThe::symbolic::symbolic_engine::Expr;

use nalgebra::DMatrix;

use std::collections::{HashMap, HashSet};
impl SubsData {
    ////////////////////////////////ELEMENT COMPOSITION AND MOLAR MASS////////////////////////////////////
    /// Calculate element composition matrix and molar masses for all substances
    ///
    /// Parses chemical formulas to determine elemental composition and molecular weights.
    /// Creates a transposed matrix where rows are substances and columns are elements.
    ///
    /// # Arguments
    /// * `groups` - Optional chemical group definitions (e.g., "Me" -> {"C":1, "H":3})
    ///
    /// # Returns
    /// * `Ok(())` - Successfully calculated composition and masses
    /// * `Err(SubsDataError)` - Formula parsing or calculation error

    pub fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> SubsDataResult<()> {
        let (hasmap_of_molar_mass, vec_of_compositions, hashset_of_elems) =
            Self::calculate_elem_composition_and_molar_mass_local(self, groups)?;
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
    ) -> SubsDataResult<(
        HashMap<String, f64>,
        Vec<HashMap<String, f64>>,
        HashSet<String>,
    )> {
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
    /// Generic transport calculator access pattern
    ///
    /// Provides safe access to transport calculators with automatic error handling
    /// and logging. Used internally by all transport-related methods.
    ///
    /// # Type Parameters
    /// * `F` - Closure that operates on the transport calculator
    /// * `R` - Return type of the closure
    ///
    /// # Arguments
    /// * `substance` - Target substance name
    /// * `f` - Closure to execute with the calculator
    fn with_transport_calculator<F, R>(&mut self, substance: &str, f: F) -> SubsDataResult<R>
    where
        F: FnOnce(&mut dyn TransportCalculator) -> Result<R, TransportError>,
    {
        let datamap = self
            .get_substance_result_mut(substance)
            .ok_or_else(|| SubsDataError::SubstanceNotFound(substance.to_string()))?;

        if datamap.contains_key(&WhatIsFound::NotFound) {
            return Err(SubsDataError::SubstanceNotFound(substance.to_string()));
        }

        let search_result = datamap
            .get_mut(&WhatIsFound::Transport)
            .ok_or_else(|| SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Transport".to_string(),
            })?
            .as_mut()
            .ok_or_else(|| SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Transport".to_string(),
            })?;

        let calculator = search_result.calculator.as_mut().ok_or_else(|| {
            SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Transport".to_string(),
            }
        })?;

        match calculator {
            CalculatorType::Transport(transport) => {
                let result: SubsDataResult<R> = f(transport).map_err(Into::into);
                if let Err(ref e) = result {
                    crate::log_error!(self.logger, e, substance, "with_transport_calculator");
                }
                result
            }
            CalculatorType::Thermo(_) => Err(SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Transport (found Thermo instead)".to_string(),
            }),
        }
    }

    /// Extract transport coefficients for a single substance at given temperature
    ///
    /// Calls the transport calculator to extract Lennard-Jones parameters and
    /// other coefficients needed for property calculations.
    ///
    /// # Arguments
    /// * `substance` - Target substance name
    /// * `temperature` - Temperature [K] for coefficient extraction
    pub fn extract_transport_coeffs(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<()> {
        if temperature <= 0.0 {
            let error = SubsDataError::InvalidTemperature(temperature);
            return self.log_and_propagate(Err(error), substance, "extract_transport_coeffs");
        }
        let result = self.with_transport_calculator(substance, |transport| {
            transport.extract_coefficients(temperature)
        });
        self.log_and_propagate(result, substance, "extract_transport_coeffs")
    }

    /// Extract transport coefficients for all substances at given temperature
    ///
    /// Batch operation that processes all substances in the dataset.
    ///
    /// # Arguments
    /// * `temperature` - Temperature [K] for coefficient extraction
    pub fn extract_all_transport_coeffs(&mut self, temperature: f64) -> SubsDataResult<()> {
        if temperature <= 0.0 {
            return Err(SubsDataError::InvalidTemperature(temperature));
        }
        for substance in self.substances.clone() {
            let result = self.extract_transport_coeffs(&substance, temperature);

            self.log_and_propagate(result, &substance, "extract_all_transport_coeffs")?;
        }
        Ok(())
    }
    /// Calculate viscosity and thermal conductivity for a substance
    ///
    /// Computes temperature-dependent transport properties using kinetic theory.
    /// Requires heat capacity and density for thermal conductivity calculation.
    ///
    /// # Arguments
    /// * `substance` - Target substance name
    /// * `temperature` - Temperature [K]
    /// * `Cp` - Heat capacity [J/mol·K] (optional, uses stored value if None)
    /// * `ro` - Density [kg/m³] (optional, calculates from ideal gas if None)
    ///
    /// # Returns
    /// * `Ok((lambda, viscosity))` - Thermal conductivity and viscosity
    /// * `Err(SubsDataError)` - Calculation or data availability error
    pub fn calculate_transport_properties(
        &mut self,
        substance: &str,
        temperature: f64,
        Cp: Option<f64>,
        ro: Option<f64>,
    ) -> SubsDataResult<(f64, f64)> {
        let result = self._calculate_transport_properties_internal(substance, temperature, Cp, ro);
        self.log_and_propagate(result, substance, "calculate_transport_properties")
    }

    fn _calculate_transport_properties_internal(
        &mut self,
        substance: &str,
        temperature: f64,
        Cp: Option<f64>,
        ro: Option<f64>,
    ) -> SubsDataResult<(f64, f64)> {
        if temperature <= 0.0 {
            return Err(SubsDataError::InvalidTemperature(temperature));
        }

        // Get library name to determine calculation method
        let lib_name = self
            .search_results
            .get(substance)
            .ok_or_else(|| SubsDataError::SubstanceNotFound(substance.to_string()))?
            .get(&WhatIsFound::Transport)
            .ok_or_else(|| SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Transport".to_string(),
            })?
            .as_ref()
            .ok_or_else(|| SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Transport".to_string(),
            })?
            .library
            .clone();

        // Extract values needed for non-CEA libraries
        let (P, M, P_unit, M_unit) = if lib_name != "CEA" {
            let P = self.P.ok_or(SubsDataError::PressureNotSet)?;
            let M = *self
                .hasmap_of_molar_mass
                .get(substance)
                .ok_or_else(|| SubsDataError::MolarMassNotFound(substance.to_string()))?;
            (P, M, self.P_unit.clone(), self.Molar_mass_unit.clone())
        } else {
            (0.0, 0.0, None, None) // Dummy values for CEA
        };

        self.with_transport_calculator(substance, |transport| {
            if lib_name == "CEA" {
                let lambda = transport.calculate_lambda(Cp, ro, temperature)?;
                let viscosity = transport.calculate_viscosity(temperature)?;
                Ok((lambda, viscosity))
            } else {
                transport.set_M(M, M_unit)?;
                transport.set_P(P, P_unit)?;
                let lambda = transport.calculate_lambda(Cp, ro, temperature)?;
                let viscosity = transport.calculate_viscosity(temperature)?;
                Ok((lambda, viscosity))
            }
        })
    }

    /// Helper function to calculate transport properties for a single substance
    fn calculate_single_transport_properties(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<HashMap<DataType, Option<f64>>> {
        let Cp = self
            .therm_map_of_properties_values
            .get(substance)
            .ok_or_else(|| SubsDataError::SubstanceNotFound(substance.to_string()))?
            .get(&DataType::Cp)
            .ok_or_else(|| SubsDataError::HeatCapacityNotAvailable(substance.to_string()))?
            .ok_or_else(|| SubsDataError::HeatCapacityNotAvailable(substance.to_string()))?;

        let ro = self
            .ro_map
            .as_ref()
            .and_then(|ro_map| ro_map.get(substance).cloned());
        let (lambda, visc) =
            self._calculate_transport_properties_internal(substance, temperature, Some(Cp), ro)?;

        let mut property_map = HashMap::new();
        property_map.insert(DataType::Lambda, Some(lambda));
        property_map.insert(DataType::Visc, Some(visc));
        Ok(property_map)
    }

    /// Calculate and populate transport_map_of_properties_values for all substances at a given temperature
    pub fn calculate_transport_map_of_properties(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<()> {
        let mut temp_map = HashMap::new();

        for substance in self.substances.clone() {
            let result = self.calculate_single_transport_properties(&substance, temperature);
            match self.log_and_propagate(
                result,
                &substance,
                "calculate_transport_map_of_properties",
            ) {
                Ok(properties) => {
                    temp_map.insert(substance, properties);
                }
                Err(_) => {} // Error already logged
            }
        }

        self.transport_map_of_properties_values = temp_map;
        Ok(())
    }

    /// Helper function to calculate transport symbolic expressions for a single substance
    fn calculate_single_transport_symbolic(
        &mut self,
        substance: &str,
    ) -> SubsDataResult<HashMap<DataType, Option<Box<Expr>>>> {
        let Cp = self
            .therm_map_of_sym
            .get(substance)
            .ok_or_else(|| SubsDataError::SubstanceNotFound(substance.to_string()))?
            .get(&DataType::Cp_sym)
            .ok_or_else(|| SubsDataError::HeatCapacityNotAvailable(substance.to_string()))?
            .as_ref()
            .ok_or_else(|| SubsDataError::HeatCapacityNotAvailable(substance.to_string()))?
            .clone();
        let Cp = *Cp;
        let P = self.P.ok_or(SubsDataError::PressureNotSet)?;
        let M = *self
            .hasmap_of_molar_mass
            .get(substance)
            .ok_or_else(|| SubsDataError::MolarMassNotFound(substance.to_string()))?;
        let P_unit = self.P_unit.clone();
        let M_unit = self.Molar_mass_unit.clone();

        let ro_expr = if let Some(ro_map) = &self.ro_map_sym {
            ro_map.get(substance).map(|v| *v.clone())
        } else {
            Some(
                RustedSciThe::symbolic::symbolic_engine::Expr::Const(P * M / 8.314)
                    / RustedSciThe::symbolic::symbolic_engine::Expr::Var("T".to_string()),
            )
        };

        self.with_transport_calculator(substance, |transport| {
            transport.set_M(M, M_unit)?;
            transport.set_P(P, P_unit)?;
            transport.create_symbolic_lambda(Some(Cp), ro_expr)?;
            transport.create_symbolic_viscosity()?;

            let mut sym_map = HashMap::new();
            sym_map.insert(
                DataType::Lambda_sym,
                transport.get_lambda_sym().ok().map(Box::new),
            );
            sym_map.insert(
                DataType::Visc_sym,
                transport.get_viscosity_sym().ok().map(Box::new),
            );

            Ok(sym_map)
        })
    }

    /// Calculate and populate transport_map_of_sym with symbolic expressions for all substances
    pub fn calculate_transport_map_of_sym(&mut self) -> SubsDataResult<()> {
        if self.therm_map_of_sym.is_empty() {
            return Err(SubsDataError::HeatCapacityNotAvailable(
                "therm_map_of_sym is empty - call calculate_therm_map_of_sym first".to_string(),
            ));
        }

        let mut temp_map = HashMap::new();
        for substance in self.substances.clone() {
            let result = self.calculate_single_transport_symbolic(&substance);
            match self.log_and_propagate(result, &substance, "calculate_transport_map_of_sym") {
                Ok(properties) => {
                    temp_map.insert(substance, properties);
                }
                Err(_) => {
                    let mut empty_map = HashMap::new();
                    empty_map.insert(DataType::Lambda_sym, None);
                    empty_map.insert(DataType::Visc_sym, None);
                    temp_map.insert(substance, empty_map);
                }
            }
        }

        self.transport_map_of_sym = temp_map;
        Ok(())
    }

    /// Helper function to calculate transport functions for a single substance
    fn calculate_single_transport_functions(
        &mut self,
        substance: &str,
    ) -> SubsDataResult<HashMap<DataType, Option<Box<dyn Fn(f64) -> f64>>>> {
        let Cp = self
            .therm_map_of_properties_values
            .get(substance)
            .ok_or_else(|| SubsDataError::SubstanceNotFound(substance.to_string()))?
            .get(&DataType::Cp)
            .ok_or_else(|| SubsDataError::HeatCapacityNotAvailable(substance.to_string()))?
            .ok_or_else(|| SubsDataError::HeatCapacityNotAvailable(substance.to_string()))?;

        let P = self.P.ok_or(SubsDataError::PressureNotSet)?;
        let M = *self
            .hasmap_of_molar_mass
            .get(substance)
            .ok_or_else(|| SubsDataError::MolarMassNotFound(substance.to_string()))?;
        let P_unit = self.P_unit.clone();
        let M_unit = self.Molar_mass_unit.clone();
        let T = self.T.ok_or(SubsDataError::TemperatureNotSet)?;
        let ro = self
            .ro_map
            .as_ref()
            .and_then(|ro_map| ro_map.get(substance).cloned());
        let ro_value = ro.unwrap_or_else(|| P * M / (8.314 * T));

        self.with_transport_calculator(substance, |transport| {
            transport.set_M(M, M_unit)?;
            transport.set_P(P, P_unit)?;
            transport.create_lambda_closure(Some(Cp), Some(ro_value))?;
            transport.create_viscosity_closure()?;

            let mut function_map = HashMap::new();
            function_map.insert(DataType::Lambda_fun, transport.get_lambda_fun().ok());
            function_map.insert(DataType::Visc_fun, transport.get_viscosity_fun().ok());

            Ok(function_map)
        })
    }

    /// Calculate and populate transport_map_of_fun with function closures for all substances
    pub fn calculate_transport_map_of_functions(&mut self) -> SubsDataResult<()> {
        let mut temp_map = HashMap::new();

        for substance in self.substances.clone() {
            let result = self.calculate_single_transport_functions(&substance);
            match self.log_and_propagate(result, &substance, "calculate_transport_map_of_functions")
            {
                Ok(properties) => {
                    temp_map.insert(substance, properties);
                }
                Err(_) => {
                    let mut empty_map = HashMap::new();
                    empty_map.insert(DataType::Lambda_fun, None);
                    empty_map.insert(DataType::Visc_fun, None);
                    temp_map.insert(substance, empty_map);
                }
            }
        }

        self.transport_map_of_fun = temp_map;
        Ok(())
    }
    ///////////////////////////////////DIFFUSION/////////////////////////////////////
    /// Initialize diffusion calculator with transport data from all substances
    ///
    /// Creates a `MultiSubstanceDiffusion` instance and populates it with
    /// transport parameters and molar masses for all substances in the dataset.
    ///
    /// # Prerequisites
    /// * Temperature and pressure must be set
    /// * Molar masses must be calculated
    /// * Transport data must be available for substances
    pub fn initialize_diffusion_data(&mut self) -> SubsDataResult<()> {
        let T = self.T.ok_or(SubsDataError::TemperatureNotSet)?;
        let P = self.P.ok_or(SubsDataError::PressureNotSet)?;

        // Create MultiSubstanceDiffusion instance
        let mut multi_diff = MultiSubstanceDiffusion::new(T, P);

        // Create HashMap of TransportData with molar masses set
        let mut transport_data_map = HashMap::new();
        let M_unit = self.Molar_mass_unit.clone();
        let P_unit = self.P_unit.clone();

        for substance in &self.substances.clone() {
            // Get molar mass for this substance
            let M = *self
                .hasmap_of_molar_mass
                .get(substance)
                .ok_or_else(|| SubsDataError::MolarMassNotFound(substance.to_string()))?;

            // Get transport data using with_transport_calculator
            let result = self.with_transport_calculator(substance, |transport| {
                transport.set_M(M, M_unit.clone())?;
                transport.set_P(P, P_unit.clone())?;

                Ok(())
            });
            match result {
                Ok(_) => {
                    // Now extract the actual TransportData from the calculator
                    if let Some(search_result) = self
                        .search_results
                        .get(substance)
                        .and_then(|map| map.get(&WhatIsFound::Transport))
                        .and_then(|opt| opt.as_ref())
                    {
                        if let Some(CalculatorType::Transport(TransportEnum::Collision(
                            transport_data,
                        ))) = &search_result.calculator
                        {
                            transport_data_map.insert(substance.clone(), transport_data.clone());
                        }
                    }
                }
                Err(_) => {
                    return Err(SubsDataError::CalculatorNotAvailable {
                        substance: substance.to_string(),
                        calc_type: "Transport".to_string(),
                    });
                }
            }
        }
        multi_diff.substances = transport_data_map;
        multi_diff.substance_names = self.substances.clone();
        self.diffusion_data = Some(multi_diff);
        Ok(())
    }

    /// Generic access pattern for diffusion calculator operations
    ///
    /// Provides safe access to the `MultiSubstanceDiffusion` instance with
    /// automatic error handling for uninitialized diffusion data.
    ///
    /// # Type Parameters
    /// * `F` - Closure that operates on the diffusion calculator
    /// * `R` - Return type of the closure
    pub fn with_diffusion_calculator<F, R>(&mut self, f: F) -> SubsDataResult<R>
    where
        F: FnOnce(&mut MultiSubstanceDiffusion) -> R,
    {
        let diffusion_data =
            self.diffusion_data
                .as_mut()
                .ok_or_else(|| SubsDataError::CalculatorNotAvailable {
                    substance: "diffusion_data".to_string(),
                    calc_type: "Diffusion".to_string(),
                })?;
        Ok(f(diffusion_data))
    }

    /// Calculate binary diffusion coefficients for all substance pairs
    ///
    /// Computes diffusion coefficients for all N(N+1)/2 unique pairs using
    /// Chapman-Enskog kinetic theory with collision integrals.
    pub fn calculate_all_pairs(&mut self) -> SubsDataResult<()> {
        self.with_diffusion_calculator(|diff| diff.calculate_all_pairs())
    }

    /// Create temperature-dependent diffusion coefficient functions
    ///
    /// Generates callable closures for all substance pairs that compute
    /// diffusion coefficients as functions of temperature.
    pub fn calculate_all_closures(&mut self) -> SubsDataResult<()> {
        self.with_diffusion_calculator(|diff| diff.calculate_all_closures())
    }

    /// Generate symbolic expressions for diffusion coefficients
    ///
    /// Creates symbolic expressions D(T) for all substance pairs using
    /// kinetic theory formulations with collision integral expressions.
    pub fn calculate_all_symbolic(&mut self) -> SubsDataResult<()> {
        self.with_diffusion_calculator(|diff| diff.calculate_all_symbolic())
    }

    /// Generate symmetric N×N diffusion coefficient matrix
    ///
    /// Creates a matrix where element [i,j] represents the binary diffusion
    /// coefficient between substances i and j. Matrix is symmetric (D_ij = D_ji).
    ///
    /// # Returns
    /// * `Ok(Vec<Vec<f64>>)` - Symmetric diffusion matrix [m²/s or cm²/s]
    /// * `Err(SubsDataError)` - Diffusion data not initialized
    pub fn create_diffusion_matrix(&mut self) -> SubsDataResult<Vec<Vec<f64>>> {
        self.with_diffusion_calculator(|diff| diff.create_diffusion_matrix())
    }

    ///////////////////////////////INPUT/OUTPUT///////////////////////////////////////////
    /// Display comprehensive search results summary
    ///
    /// Prints formatted output showing which substances were found in which
    /// libraries, their priority levels, and available calculator types.
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
    /// Fallback to NIST database for substances not found in local libraries
    ///
    /// Automatically queries NIST Chemistry WebBook for thermodynamic data
    /// when substances are not available in local database libraries.
    ///
    /// # Behavior
    /// * Only attempts NIST lookup for thermodynamic libraries
    /// * Uses substance phase information if available
    /// * Updates search results with NIST data on success
    pub fn if_not_found_go_NIST(&mut self) -> SubsDataResult<()> {
        use crate::Thermodynamics::DBhandlers::NIST_parser::{Phase, SearchType};
        let list_of_thermo_libs = [
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
                            match thermo.renew_base(substance.clone(), SearchType::All, phase) {
                                Ok(()) => {
                                    println!("\nFound data in NIST for {}", substance);
                                    // Create a new SearchResult with the NIST data
                                    let search_result = SearchResult {
                                        library: "NIST".to_string(),
                                        priority_type: LibraryPriority::Permitted,
                                        data: serde_json::Value::Null,
                                        calculator: Some(CalculatorType::Thermo(thermo.clone())),
                                    };
                                    // Update the search_results HashMap
                                    let mut result_map = HashMap::new();
                                    result_map.insert(WhatIsFound::Thermo, Some(search_result));
                                    // Remove the NotFound entry and add the new Thermo entry
                                    self.search_results.insert(substance.clone(), result_map);
                                    println!(
                                        "Updated search results with NIST data for {}",
                                        substance
                                    );
                                }
                                Err(e) => {
                                    println!(
                                        "Failed to find data in NIST for {}: {}",
                                        substance, e
                                    );
                                }
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
