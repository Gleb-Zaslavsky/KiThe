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

use crate::Thermodynamics::DBhandlers::thermo_api::{ThermoCalculator, ThermoEnum};
use crate::Thermodynamics::DBhandlers::transport_api::{
    TransportCalculator, TransportEnum, TransportError,
};
use crate::Thermodynamics::User_substances::{
    CalculatorType, DataType, LibraryPriority, Phases, PropertySearchState, SearchResult, SubsData,
    TransportCalculationMode, WhatIsFound,
};

use crate::Thermodynamics::DBhandlers::Diffusion::MultiSubstanceDiffusion;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::thermo_lib_api::{LibraryCapability, ThermoData};
use RustedSciThe::symbolic::symbolic_engine::Expr;

use nalgebra::DMatrix;
use tabled::settings::Style;
use tabled::{Table, Tabled};

use std::collections::{HashMap, HashSet};

/// One row in the typed search summary report.
#[derive(Clone, Debug, PartialEq, Eq, Tabled)]
pub struct SearchSummaryRow {
    #[tabled(rename = "Substance")]
    substance: String,
    #[tabled(rename = "Property")]
    property: String,
    #[tabled(rename = "State")]
    state: String,
    #[tabled(rename = "Library")]
    library: String,
    #[tabled(rename = "Source record")]
    record_key: String,
    #[tabled(rename = "Priority")]
    priority: String,
    #[tabled(rename = "Calculator")]
    calculator: String,
}

impl SearchSummaryRow {
    /// Canonical substance name covered by this lookup row.
    pub fn substance(&self) -> &str {
        &self.substance
    }

    /// Requested property family, such as `Thermo` or `Transport`.
    pub fn property(&self) -> &str {
        &self.property
    }

    /// Resolution state rendered for this property.
    pub fn state(&self) -> &str {
        &self.state
    }

    /// Selected source library, or `-` when no record was selected.
    pub fn library(&self) -> &str {
        &self.library
    }

    /// Exact record key selected inside the source library.
    pub fn record_key(&self) -> &str {
        &self.record_key
    }

    /// Lookup tier that selected the record.
    pub fn priority(&self) -> &str {
        &self.priority
    }

    /// Calculator family attached to the resolved record.
    pub fn calculator(&self) -> &str {
        &self.calculator
    }
}

/// Typed summary of the canonical search state.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SearchSummaryReport {
    rows: Vec<SearchSummaryRow>,
    total_substances: usize,
    priority_found: usize,
    permitted_found: usize,
    not_found: usize,
}

impl SearchSummaryReport {
    /// Returns all report rows in deterministic order.
    pub fn rows(&self) -> &[SearchSummaryRow] {
        &self.rows
    }

    /// Returns the number of substances tracked by the report.
    pub fn total_substances(&self) -> usize {
        self.total_substances
    }

    /// Returns the number of priority-found substances.
    pub fn priority_found(&self) -> usize {
        self.priority_found
    }

    /// Returns the number of permitted-found substances.
    pub fn permitted_found(&self) -> usize {
        self.permitted_found
    }

    /// Returns the number of substances that were not found.
    pub fn not_found(&self) -> usize {
        self.not_found
    }

    /// Renders the report as a formatted table without mutating state.
    pub fn render_table(&self) -> String {
        if self.rows.is_empty() {
            return "Search summary: no substances have been resolved yet.".to_string();
        }

        Table::new(self.rows.clone())
            .with(Style::modern_rounded())
            .to_string()
    }
}
impl SubsData {
    ////////////////////////////////ELEMENT COMPOSITION AND MOLAR MASS////////////////////////////////////
    /// Resolves composition and molar mass for one substance.
    ///
    /// The calculator-provided composition is preferred when it exists because
    /// it is the canonical source for already-fitted thermo records. If the
    /// calculator does not provide composition data, we fall back to parsing
    /// the chemical formula string.
    fn resolve_composition_and_molar_mass_for_substance(
        &self,
        substance: &str,
        groups: Option<&HashMap<String, HashMap<String, usize>>>,
        thermo: Option<&ThermoEnum>,
    ) -> SubsDataResult<(HashMap<String, f64>, f64)> {
        use crate::Kinetics::molmass::{
            calculate_molar_mass, calculate_molar_mass_for_composition,
        };

        if let Some(thermo) = thermo {
            if let Ok(Some(composition)) = thermo.get_composition() {
                let composition_usize: HashMap<String, usize> = composition
                    .iter()
                    .map(|(k, v)| (k.clone(), *v as usize))
                    .collect();
                let molar_mass = calculate_molar_mass_for_composition(composition_usize.clone());
                let composition = composition
                    .into_iter()
                    .map(|(k, v)| (k, v as f64))
                    .collect();
                return Ok((composition, molar_mass));
            }
        }

        let (molar_mass, composition) =
            calculate_molar_mass(substance.to_string(), groups.cloned()).map_err(|e| {
                SubsDataError::calculation_failed_with_source(
                    substance.to_string(),
                    "molar mass calculation",
                    e.to_string(),
                    e,
                )
            })?;

        let composition = composition
            .into_iter()
            .map(|(k, v)| (k, v as f64))
            .collect();
        Ok((composition, molar_mass))
    }

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
        let (molar_mass_by_substance, vec_of_compositions, hashset_of_elems) =
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
        self.molar_mass_by_substance = molar_mass_by_substance;
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
        let mut vec_of_compositions: Vec<HashMap<String, f64>> = Vec::new();
        let mut hashset_of_elems: HashSet<String> = HashSet::new();
        let mut molar_mass_by_substance: HashMap<String, f64> = HashMap::new();

        for substance in &sd.substances {
            let thermo =
                sd.get_substance_search_state(substance)
                    .and_then(|state| match &state.thermo {
                        PropertySearchState::Found(SearchResult {
                            calculator: Some(CalculatorType::Thermo(thermo)),
                            ..
                        }) => Some(thermo),
                        _ => None,
                    });

            let (composition, molar_mass) = sd.resolve_composition_and_molar_mass_for_substance(
                substance,
                groups.as_ref(),
                thermo,
            )?;

            molar_mass_by_substance.insert(substance.clone(), molar_mass);
            vec_of_compositions.push(composition.clone());
            let elements = composition.keys().map(|el| el.clone()).collect::<Vec<_>>();
            hashset_of_elems.extend(elements);
        }
        Ok((
            molar_mass_by_substance,
            vec_of_compositions,
            hashset_of_elems,
        ))
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
        let search_result =
            self.get_canonical_search_result_mut(substance, WhatIsFound::Transport)?;
        let calculator = search_result.calculator_mut().ok_or_else(|| {
            SubsDataError::CalculatorNotAvailable {
                substance: substance.to_string(),
                calc_type: "Transport".to_string(),
            }
        })?;

        match calculator {
            CalculatorType::Transport(transport) => {
                let result: SubsDataResult<R> = f(&mut *transport).map_err(Into::into);
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
        let lib_name = self.get_canonical_search_library(substance, WhatIsFound::Transport)?;

        // Extract the exact transport backend inputs instead of routing dummy
        // placeholder numbers through the solver.
        let transport_mode = if lib_name == "CEA" {
            TransportCalculationMode::Cea
        } else {
            let pressure = self.P.ok_or(SubsDataError::PressureNotSet)?;
            let molar_mass = *self
                .molar_mass_by_substance
                .get(substance)
                .ok_or_else(|| SubsDataError::MolarMassNotFound(substance.to_string()))?;
            TransportCalculationMode::Standard {
                pressure,
                molar_mass,
                pressure_unit: self.P_unit.clone(),
                molar_mass_unit: self.Molar_mass_unit.clone(),
            }
        };

        self.with_transport_calculator(substance, |transport| {
            if matches!(&transport_mode, TransportCalculationMode::Cea) {
                let lambda = transport.calculate_lambda(Cp, ro, temperature)?;
                let viscosity = transport.calculate_viscosity(temperature)?;
                Ok((lambda, viscosity))
            } else if let TransportCalculationMode::Standard {
                pressure,
                molar_mass,
                pressure_unit,
                molar_mass_unit,
            } = &transport_mode
            {
                transport.set_M(*molar_mass, molar_mass_unit.clone())?;
                transport.set_P(*pressure, pressure_unit.clone())?;
                let lambda = transport.calculate_lambda(Cp, ro, temperature)?;
                let viscosity = transport.calculate_viscosity(temperature)?;
                Ok((lambda, viscosity))
            } else {
                unreachable!("transport_mode is fully covered by the CEA and standard branches")
            }
        })
    }

    /// Helper function to calculate transport properties for a single substance
    fn calculate_single_transport_properties(
        &mut self,
        substance: &str,
        temperature: f64,
    ) -> SubsDataResult<HashMap<DataType, Option<f64>>> {
        // Read the thermodynamic cache through the explicit readiness state so
        // transport code cannot accidentally treat a missing value as a valid zero.
        let Cp = self
            .therm_value_state(substance, DataType::Cp)
            .value()
            .copied()
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
            let properties = self.log_and_propagate(
                result,
                &substance,
                "calculate_transport_map_of_properties",
            )?;
            temp_map.insert(substance, properties);
        }

        self.transport_map_of_properties_values = temp_map;
        Ok(())
    }

    /// Helper function to calculate transport symbolic expressions for a single substance
    fn calculate_single_transport_symbolic(
        &mut self,
        substance: &str,
    ) -> SubsDataResult<HashMap<DataType, Option<Box<Expr>>>> {
        // Read the symbolic heat capacity through the explicit readiness state.
        let Cp = self
            .therm_symbolic_state(substance, DataType::Cp_sym)
            .value()
            .cloned()
            .ok_or_else(|| SubsDataError::HeatCapacityNotAvailable(substance.to_string()))?;
        let P = self.P.ok_or(SubsDataError::PressureNotSet)?;
        let M = *self
            .molar_mass_by_substance
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
            let properties =
                self.log_and_propagate(result, &substance, "calculate_transport_map_of_sym")?;
            temp_map.insert(substance, properties);
        }

        self.transport_map_of_sym = temp_map;
        Ok(())
    }

    /// Helper function to calculate transport functions for a single substance
    fn calculate_single_transport_functions(
        &mut self,
        substance: &str,
    ) -> SubsDataResult<HashMap<DataType, Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>>> {
        // Reuse cached Cp when available, but fall back to an on-demand thermo
        // evaluation if a temperature change invalidated the property cache.
        let Cp = match self.therm_value_state(substance, DataType::Cp).value() {
            Some(cp) => *cp,
            None => {
                let temperature = self.T.ok_or(SubsDataError::TemperatureNotSet)?;
                let (cp, _, _) = self.calculate_thermo_properties(substance, temperature)?;
                cp
            }
        };

        let P = self.P.ok_or(SubsDataError::PressureNotSet)?;
        let M = *self
            .molar_mass_by_substance
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
            let _ = transport.create_lambda_closure(Some(Cp), Some(ro_value))?;
            let _ = transport.create_viscosity_closure()?;

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
            let properties =
                self.log_and_propagate(result, &substance, "calculate_transport_map_of_functions")?;
            temp_map.insert(substance, properties);
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
        let mut multi_diff = MultiSubstanceDiffusion::new(T, P)?;

        // Create HashMap of TransportData with molar masses set
        let mut transport_data_map = HashMap::new();
        let M_unit = self.Molar_mass_unit.clone();
        let P_unit = self.P_unit.clone();

        for substance in &self.substances.clone() {
            // Get molar mass for this substance
            let M = *self
                .molar_mass_by_substance
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
                    // Extract the transport data from the canonical typed state after
                    // the calculator update has completed.
                    if let Some(search_result) =
                        self.get_canonical_search_result(substance, WhatIsFound::Transport)
                    {
                        if let Some(CalculatorType::Transport(TransportEnum::Collision(
                            transport_data,
                        ))) = search_result.calculator()
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
        self.with_diffusion_calculator(|diff| {
            diff.calculate_all_pairs().map_err(SubsDataError::from)
        })??;
        Ok(())
    }

    /// Create temperature-dependent diffusion coefficient functions
    ///
    /// Generates callable closures for all substance pairs that compute
    /// diffusion coefficients as functions of temperature.
    pub fn calculate_all_closures(&mut self) -> SubsDataResult<()> {
        self.with_diffusion_calculator(|diff| {
            diff.calculate_all_closures().map_err(SubsDataError::from)
        })??;
        Ok(())
    }

    /// Generate symbolic expressions for diffusion coefficients
    ///
    /// Creates symbolic expressions D(T) for all substance pairs using
    /// kinetic theory formulations with collision integral expressions.
    pub fn calculate_all_symbolic(&mut self) -> SubsDataResult<()> {
        self.with_diffusion_calculator(|diff| {
            diff.calculate_all_symbolic().map_err(SubsDataError::from)
        })??;
        Ok(())
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
        self.with_diffusion_calculator(|diff| {
            diff.create_diffusion_matrix().map_err(SubsDataError::from)
        })?
        .map_err(Into::into)
    }

    ///////////////////////////////INPUT/OUTPUT///////////////////////////////////////////
    /// Builds a typed summary of the canonical search state.
    ///
    /// The report is deterministic and can be rendered by CLI, GUI, or tests
    /// without touching the underlying search caches.
    pub fn search_summary_report(&self) -> SearchSummaryReport {
        let mut rows = Vec::new();
        let mut substances = self.search_states.keys().cloned().collect::<Vec<_>>();
        substances.sort();

        for substance in substances {
            if let Some(state) = self.search_states.get(&substance) {
                if state.is_missing() {
                    rows.push(SearchSummaryRow {
                        substance: substance.clone(),
                        property: "Thermo/Transport".to_string(),
                        state: "Missing".to_string(),
                        library: "-".to_string(),
                        record_key: "-".to_string(),
                        priority: "-".to_string(),
                        calculator: "-".to_string(),
                    });
                    continue;
                }

                for (property_name, property_state) in
                    [("Thermo", state.thermo()), ("Transport", state.transport())]
                {
                    let row = match property_state {
                        PropertySearchState::NotSearched => SearchSummaryRow {
                            substance: substance.clone(),
                            property: property_name.to_string(),
                            state: "Pending".to_string(),
                            library: "-".to_string(),
                            record_key: "-".to_string(),
                            priority: "-".to_string(),
                            calculator: "-".to_string(),
                        },
                        PropertySearchState::Missing => SearchSummaryRow {
                            substance: substance.clone(),
                            property: property_name.to_string(),
                            state: "Missing".to_string(),
                            library: "-".to_string(),
                            record_key: "-".to_string(),
                            priority: "-".to_string(),
                            calculator: "-".to_string(),
                        },
                        PropertySearchState::Failed(message) => SearchSummaryRow {
                            substance: substance.clone(),
                            property: property_name.to_string(),
                            state: format!("Failed: {}", message),
                            library: "-".to_string(),
                            record_key: "-".to_string(),
                            priority: "-".to_string(),
                            calculator: "-".to_string(),
                        },
                        PropertySearchState::Found(result) => SearchSummaryRow {
                            substance: substance.clone(),
                            property: property_name.to_string(),
                            state: "Found".to_string(),
                            library: result.library().to_string(),
                            record_key: result.record_key().to_string(),
                            priority: match result.priority_type() {
                                LibraryPriority::Priority => "Priority".to_string(),
                                LibraryPriority::Permitted => "Permitted".to_string(),
                                LibraryPriority::Explicit => "Explicit".to_string(),
                            },
                            calculator: match result.calculator() {
                                Some(CalculatorType::Thermo(_)) => "Thermo".to_string(),
                                Some(CalculatorType::Transport(_)) => "Transport".to_string(),
                                None => "None".to_string(),
                            },
                        },
                    };
                    rows.push(row);
                }
            }
        }

        SearchSummaryReport {
            rows,
            total_substances: self.substances.len(),
            priority_found: self.get_priority_found_substances().len(),
            permitted_found: self.get_permitted_found_substances().len(),
            not_found: self.get_not_found_substances().len(),
        }
    }

    /// Resolves the current query and returns a typed search report.
    ///
    /// This is the narrowest common workflow facade for callers that want to
    /// configure a query, resolve records, and immediately inspect the result
    /// summary without manually chaining several separate calls.
    pub fn resolve_search_and_report(&mut self) -> SubsDataResult<SearchSummaryReport> {
        self.search_substances()?;
        Ok(self.search_summary_report())
    }

    /// Display a comprehensive search results summary.
    ///
    /// This is a formatting adapter over `search_summary_report()`, so the
    /// rendering path stays separate from the canonical data path.
    pub fn print_search_summary(&self) {
        let report = self.search_summary_report();
        println!("\nSearch Results Summary:");
        println!("{}", report.render_table());
        println!(
            "\nStatistics: total={}, priority_found={}, permitted_found={}, not_found={}",
            report.total_substances(),
            report.priority_found(),
            report.permitted_found(),
            report.not_found()
        );
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
        let not_found_substances = self.get_not_found_substances();
        if not_found_substances.is_empty() {
            return Ok(());
        }

        let has_thermo_lookup = self.library_priorities.keys().any(|library| {
            matches!(
                ThermoData::library_capability(library),
                Some(LibraryCapability::Thermo)
            )
        });
        if !has_thermo_lookup {
            let reason = "no thermodynamic libraries are enabled for NIST fallback".to_string();
            let substance = not_found_substances
                .first()
                .cloned()
                .unwrap_or_else(|| "unknown".to_string());
            return Err(SubsDataError::nist_retrieval_failed(substance, reason));
        }

        let mut first_error: Option<SubsDataError> = None;
        let mut found_any = false;

        for substance in not_found_substances {
            let mut calculator = self.create_calculator("NIST")?;
            if let CalculatorType::Thermo(thermo) = &mut calculator {
                let phase = if let Some(Some(phase)) = self.map_of_phases.get(&substance) {
                    match phase {
                        Phases::Gas => Phase::Gas,
                        Phases::Solid => Phase::Solid,
                        Phases::Liquid => Phase::Liquid,
                        _ => Phase::Gas,
                    }
                } else {
                    Phase::Gas
                };

                match thermo.renew_base(substance.clone(), SearchType::All, phase) {
                    Ok(()) => {
                        found_any = true;
                        self.store_search_result(
                            &substance,
                            WhatIsFound::Thermo,
                            "NIST".to_string(),
                            LibraryPriority::Permitted,
                            serde_json::Value::Null,
                            CalculatorType::Thermo(thermo.clone()),
                            true,
                        );
                    }
                    Err(e) => {
                        let error = SubsDataError::nist_retrieval_failed_with_source(
                            substance.clone(),
                            e.to_string(),
                            e,
                        );
                        if first_error.is_none() {
                            first_error = Some(error);
                        }
                    }
                }
            }
        }

        if found_any {
            Ok(())
        } else {
            Err(first_error.unwrap_or_else(|| {
                SubsDataError::nist_retrieval_failed(
                    "unknown",
                    "NIST fallback failed for every unresolved substance",
                )
            }))
        }
    }
}
