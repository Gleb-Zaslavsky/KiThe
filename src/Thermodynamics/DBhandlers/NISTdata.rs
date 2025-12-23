//! # NIST Thermodynamic Data Handler Module
//!
//! ## Aim
//! This module provides a high-level interface to NIST thermodynamic data by wrapping the NIST_parser
//! functionality. It handles NIST Shomate equation format for thermodynamic property calculations
//! and provides integration with the broader thermodynamic calculation framework.
//!
//! ## Main Data Structures and Logic
//! - `NISTdata`: Wrapper structure around `NistInput` providing unified interface
//! - Uses NIST Shomate equation format: Cp = A + BT + CT² + DT³ + E/T²
//! - Integrates with `NistParser` for automatic web scraping of NIST Chemistry WebBook
//! - Supports gas, liquid, and solid phases with appropriate correlations
//!
//! ## Key Methods
//! - `get_data_from_NIST()`: Automatically fetches data from NIST Chemistry WebBook
//! - `extract_coefficients(T)`: Selects appropriate coefficient set for temperature
//! - `calculate_Cp_dH_dS(T)`: Computes thermodynamic properties using Shomate equations
//! - `create_closures_Cp_dH_dS()`: Creates efficient closure functions
//! - `create_sym_Cp_dH_dS()`: Creates symbolic expressions for analytical work
//! - `Taylor_series_Cp_dH_dS()`: Generates Taylor series expansions
//!
//! ## Usage
//! ```rust, ignore
//! let mut nist = NISTdata::new();
//! nist.get_data_from_NIST("CO".to_string(), SearchType::All, Phase::Gas)?;
//! nist.extract_coefficients(400.0)?;
//! nist.calculate_Cp_dH_dS(400.0);
//! let cp = nist.Cp;  // Heat capacity
//! let dh = nist.dh;  // Enthalpy
//! let ds = nist.ds;  // Entropy
//! ```
//!
//! ## Interesting Features
//! - Automatic web scraping integration with NIST Chemistry WebBook
//! - Supports real-time data fetching for any substance in NIST database
//! - Handles multiple temperature ranges with automatic coefficient selection
//! - Implements complex Clone trait with proper function reconstruction
//! - Provides seamless integration between web-scraped data and calculation framework
//! - Uses Shomate equation format optimized for accurate thermodynamic predictions
//! - Supports both numerical and symbolic computation modes

use crate::Thermodynamics::DBhandlers::NIST_parser::{NistInput, NistParser, Phase, SearchType};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use prettytable::{Cell, Row, Table};
use serde_json::Value;
use std::collections::HashMap;
use std::fmt;
use std::{error::Error, fmt::Debug};

#[allow(non_upper_case_globals, non_snake_case)]
const e2: Expr = Expr::Const(2.0);
#[allow(non_upper_case_globals, non_snake_case)]
const e3: Expr = Expr::Const(3.0);
#[allow(non_upper_case_globals, non_snake_case)]
const e4: Expr = Expr::Const(4.0);

#[derive(Debug)]
pub enum NISTError {
    NoCoefficientsFound { temperature: f64, range: String },
    InvalidTemperatureRange,
    SerdeError(serde_json::Error),
    UnsupportedUnit(String),
    SymbolicError(String),
    FittingError(String),
}

impl fmt::Display for NISTError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            NISTError::NoCoefficientsFound { temperature, range } => {
                write!(
                    f,
                    "No coefficients found for temperature {} K. Valid range: {}",
                    temperature, range
                )
            }
            NISTError::InvalidTemperatureRange => {
                write!(f, "Invalid temperature range in coefficient data")
            }
            NISTError::SerdeError(msg) => {
                write!(f, "Failed to deserialize NASA data: {}", msg)
            }
            NISTError::UnsupportedUnit(unit) => {
                write!(
                    f,
                    "Unsupported unit: {}. Only 'J' and 'cal' are supported",
                    unit
                )
            }
            NISTError::SymbolicError(msg) => {
                write!(f, "Symbolic error: {}", msg)
            }
            NISTError::FittingError(msg) => {
                write!(f, "Fitting error: {}", msg)
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct FittingReport {
    pub l2_norm: f64,
    pub max_norm: f64,
    pub number_of_significant_points: usize,
    pub T_range: Option<(f64, f64)>,
}

impl Error for NISTError {}

impl From<serde_json::Error> for NISTError {
    fn from(err: serde_json::Error) -> Self {
        NISTError::SerdeError(err)
    }
}

#[derive(Debug, Clone)]
pub struct Coeffs {
    pub T: (f64, f64),
    pub coeff: (f64, f64, f64, f64, f64, f64, f64, f64),
}

pub struct NISTdata {
    /// data parsed from library
    pub input: NistInput,
    /// Joules or Cal
    pub unit: Option<String>,
    /// unit multiplier for conversions
    pub unit_multiplier: f64,

    pub coeffs_map: HashMap<usize, Coeffs>,
    /// extracted coefficients for current temperature
    pub coeffs: Option<(f64, f64, f64, f64, f64, f64, f64, f64)>,
    /// heat capacity value at T
    pub Cp: f64,
    /// enthalpy value at T
    pub dh: f64,
    /// entropy value at T
    pub ds: f64,
    /// heat capacity function
    #[allow(clippy::type_complexity)]
    pub C_fun: Box<dyn Fn(f64) -> f64 + 'static>,
    /// enthalpy function
    #[allow(clippy::type_complexity)]
    pub dh_fun: Box<dyn Fn(f64) -> f64 + 'static>,
    /// entropy function
    #[allow(clippy::type_complexity)]
    pub ds_fun: Box<dyn Fn(f64) -> f64 + 'static>,
    /// symbolic heat capacity
    pub Cp_sym: Expr,
    /// symbolic enthalpy
    pub dh_sym: Expr,
    /// symbolic entropy
    pub ds_sym: Expr,
    // user defined inteval
    pub T_interval: Option<(f64, f64)>,

    pub current_T_range: Option<(f64, f64)>,
}

impl NISTdata {
    pub fn new() -> Self {
        let input = NistInput::new();
        Self {
            input: input,
            unit: None,
            unit_multiplier: 1.0,
            coeffs_map: HashMap::new(),
            coeffs: None,
            Cp: 0.0,
            dh: 0.0,
            ds: 0.0,
            C_fun: Box::new(|x| x),
            dh_fun: Box::new(|x| x),
            ds_fun: Box::new(|x| x),
            Cp_sym: Expr::Const(0.0),
            dh_sym: Expr::Const(0.0),
            ds_sym: Expr::Const(0.0),
            T_interval: None,
            current_T_range: None,
        }
    }
    /// set energy units: J or calories
    pub fn set_unit(&mut self, unit: &str) -> Result<(), NISTError> {
        match unit {
            "J" => {
                self.unit = Some("J".to_string());
                self.unit_multiplier = 1.0;
            }
            "cal" => {
                self.unit = Some("cal".to_string());
                self.unit_multiplier = 1.0 / 4.184;
            }
            _ => return Err(NISTError::UnsupportedUnit(unit.to_string())),
        }
        Ok(())
    }
    pub fn set_T_interval(&mut self, T_min: f64, T_max: f64) {
        self.T_interval = Some((T_min, T_max));
    }
    /// takes serde Value and parse it into structure
    pub fn from_serde(&mut self, serde: Value) -> Result<(), NISTError> {
        self.input = serde_json::from_value(serde).map_err(|e| NISTError::SerdeError(e))?;
        Ok(())
    }

    pub fn parse_coefficients(&mut self) -> Result<(), NISTError> {
        for (i, T_pairs) in self.input.T.clone().unwrap().iter().enumerate() {
            let coeffs = self.input.cp.clone().unwrap()[i].clone();
            let (a, b, c, d, e, f, g, h) = (
                coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], coeffs[5], coeffs[6],
                coeffs[7],
            );
            let coeffs = Coeffs {
                T: (T_pairs[0], T_pairs[1]),
                coeff: (a, b, c, d, e, f, g, h),
            };
            self.coeffs_map.insert(i, coeffs);
        }
        Ok(())
    }

    fn find_coefficients_for_temperature(
        &mut self,
        t: f64,
    ) -> Result<(f64, f64, f64, f64, f64, f64, f64, f64), NISTError> {
        for coeffs in self.coeffs_map.values() {
            if coeffs.T.0 <= t && t <= coeffs.T.1 {
                self.current_T_range = Some(coeffs.T.clone());
                return Ok(coeffs.coeff);
            }
        }

        let ranges: Vec<String> = self
            .coeffs_map
            .values()
            .map(|c| format!("{}-{}", c.T.0, c.T.1))
            .collect();

        Err(NISTError::NoCoefficientsFound {
            temperature: t,
            range: ranges.join(", "),
        })
    }

    pub fn extract_coefficients(&mut self, T: f64) -> Result<(), NISTError> {
        if self.coeffs_map.is_empty() {
            self.parse_coefficients()?;
        }
        self.coeffs = Some(self.find_coefficients_for_temperature(T)?);
        Ok(())
    }

    pub fn Taylor_series_Cp_dH_dS(
        self,
        T0: f64,
        n: usize,
    ) -> Result<(Expr, Expr, Expr), NISTError> {
        let Cp = self.Cp_sym.clone();

        let Cp_taylor = Cp.taylor_series1D("T", T0, n);

        let dh = self.dh_sym.clone();

        let dh_taylor = dh.taylor_series1D("T", T0, n);

        let ds = self.ds_sym.clone();

        let ds_taylor = ds.taylor_series1D("T", T0, n);
        Ok((Cp_taylor, dh_taylor, ds_taylor))
    }
    pub fn get_data_from_NIST(
        &mut self,
        sub_name: String,
        search_type: SearchType,
        phase: Phase,
    ) -> Result<(), NISTError> {
        let parser_instance = NistParser::new();
        match parser_instance.get_data(&sub_name, search_type, phase) {
            Ok(input) => {
                self.input = input;
                self.parse_coefficients()?;
                Ok(())
            }
            Err(e) => Err(NISTError::NoCoefficientsFound {
                temperature: 0.0, // Default temperature since we don't have a specific one
                range: format!("Failed to get data for {}: {}", sub_name, e),
            }),
        }
    }

    pub fn integr_mean(&mut self) -> Result<(), NISTError> {
        use RustedSciThe::symbolic::symbolic_integration::QuadMethod;
        let (T_min, T_max) = self.T_interval.ok_or(NISTError::InvalidTemperatureRange)?;

        if T_min >= T_max {
            return Err(NISTError::InvalidTemperatureRange);
        }

        let T_range = T_max - T_min;
        let inv_T_range = 1.0 / T_range;

        // Integrate entropy
        let ds_sym = self.ds_sym.clone();

        let ds_integral = ds_sym
            .quad(QuadMethod::GaussLegendre, 3, T_min, T_max, None)
            .map_err(|e| NISTError::SymbolicError(format!("Failed to integrate entropy: {}", e)))?;
        let ds_mean = inv_T_range * ds_integral;

        // Integrate enthalpy
        let dh_sym = self.dh_sym.clone().simplify();
        let dh_integral = dh_sym
            .quad(QuadMethod::GaussLegendre, 3, T_min, T_max, None)
            .map_err(|e| {
                NISTError::SymbolicError(format!("Failed to integrate enthalpy: {}", e))
            })?;
        let dh_mean = inv_T_range * dh_integral;

        // Integrate heat capacity
        let Cp_sym = self.Cp_sym.clone();
        let Cp_integral = Cp_sym
            .quad(QuadMethod::GaussLegendre, 3, T_min, T_max, None)
            .map_err(|e| {
                NISTError::SymbolicError(format!("Failed to integrate heat capacity: {}", e))
            })?;
        let Cp_mean = inv_T_range * Cp_integral;

        // Store mean values
        self.Cp = Cp_mean;
        self.dh = dh_mean;
        self.ds = ds_mean;

        Ok(())
    }

    pub fn calculate_cp_dh_ds(&mut self, T: f64) -> Result<(), NISTError> {
        let um = self.unit_multiplier;
        if let Some(coeffs) = self.coeffs.clone() {
            let T = T / 1000.0;
            let (a, b, c, d, e, f, g, h) = (
                coeffs.0, coeffs.1, coeffs.2, coeffs.3, coeffs.4, coeffs.5, coeffs.6, coeffs.7,
            );
            self.Cp = um * calculate_cp(T, a, b, c, d, e);
            let dh0 = self.input.dh.unwrap_or(0.0); // Use 0.0 for simple substances
            // 1000 is to turn kJ into J
            self.dh = 1000.0 * um * (calculate_dh(T, a, b, c, d, e, f, g, h) + dh0);
            self.ds = um * calculate_s(T, a, b, c, d, e, f, g, h);
            return Ok(());
        }
        Err(NISTError::NoCoefficientsFound {
            temperature: T * 1000.0,
            range: "No coefficients available".to_string(),
        })
    }

    pub fn create_sym_cp_dh_ds(&mut self) -> Result<(), NISTError> {
        let um = Expr::Const(self.unit_multiplier);
        if let Some(coeffs) = self.coeffs.clone() {
            let (a, b, c, d, e, f, g, h) = (
                coeffs.0, coeffs.1, coeffs.2, coeffs.3, coeffs.4, coeffs.5, coeffs.6, coeffs.7,
            );
            self.Cp_sym = um.clone() * calculate_cp_sym(a, b, c, d, e);
            let dh0 = self.input.dh.unwrap_or(0.0); // Use 0.0 for simple substances
            self.dh_sym = Expr::Const(1000.0)
                * um.clone()
                * (calculate_dh_sym(a, b, c, d, e, f, g, h) + Expr::Const(dh0));
            self.ds_sym = um.clone() * calculate_s_sym(a, b, c, d, e, f, g, h);
            return Ok(());
        }
        Err(NISTError::NoCoefficientsFound {
            temperature: 0.0,
            range: "No coefficients available".to_string(),
        })
    }

    pub fn create_closure_cp_dh_ds(&mut self) -> Result<(), NISTError> {
        let um = self.unit_multiplier;
        if let Some(coeffs) = self.coeffs.clone() {
            let (a, b, c, d, e, f, g, h) = (
                coeffs.0, coeffs.1, coeffs.2, coeffs.3, coeffs.4, coeffs.5, coeffs.6, coeffs.7,
            );
            self.C_fun = Box::new(move |t| um * calculate_cp(t / 1000.0, a, b, c, d, e));
            let dh0 = self.input.dh.unwrap_or(0.0); // Use 0.0 for simple substances
            self.dh_fun = Box::new(move |t| {
                1000.0 * um * (calculate_dh(t / 1000.0, a, b, c, d, e, f, g, h) + dh0)
            });
            self.ds_fun = Box::new(move |t| um * calculate_s(t / 1000.0, a, b, c, d, e, f, g, h));
            return Ok(());
        }
        Err(NISTError::NoCoefficientsFound {
            temperature: 0.0,
            range: "No coefficients available".to_string(),
        })
    }

    pub fn pretty_print(&self) -> Result<(), NISTError> {
        let temps = self.input.T.clone().expect("no temperature parsed");
        let coeffs = self.input.cp.clone().expect("no coefficients parsed");
        let mut table = Table::new();
        let mut header_row = vec![Cell::new("Coefficients")];
        for T in temps.clone() {
            let Tstr = format!("{:} - {:}   ", T[0], T[1]);
            header_row.push(Cell::new(&Tstr));
        }
        table.add_row(Row::new(header_row));

        let coeffs_names = ["A", "B", "C", "D", "E", "F", "G", "H"];
        for (i, coeff_name) in coeffs_names.iter().enumerate() {
            let mut row = vec![Cell::new(coeff_name)];
            for coeff_for_every_T in coeffs.iter() {
                let coeff = coeff_for_every_T[i].to_string();
                row.push(Cell::new(&format!("{:}", coeff)));
            }
            table.add_row(Row::new(row));
        }
        table.printstd();

        let mut table = Table::new();
        let M = self.input.molar_mass.clone().expect("no molar mass parsed");
        let dH = self.input.dh.unwrap_or(0.0); // Use 0.0 for simple substances
        let dS = self.input.ds.clone().expect("no dS parsed");
        let header_row = vec![Cell::new("Molar mass"), Cell::new("dH"), Cell::new("dS")];
        table.add_row(Row::new(header_row));
        let row = vec![
            Cell::new(&M.to_string()),
            Cell::new(&dH.to_string()),
            Cell::new(&dS.to_string()),
        ];
        table.add_row(Row::new(row));
        table.printstd();
        Ok(())
    }

    /////////////////////////////////TEMPERATURE RANGE////////////////////////////////////////////////////////////
    fn is_this_T_from_current_T_range(&self, t: f64) -> bool {
        if let Some((T_min, T_max)) = self.current_T_range {
            T_min <= t && t <= T_max
        } else {
            false
        }
    }

    /// Creates closure functions for Cp, dH, dS with automatic temperature range handling
    ///
    /// Automatically selects appropriate NIST Shomate coefficients for the given temperature
    /// and creates closure functions. If temperature is outside current range,
    /// updates coefficients before creating closures.
    ///
    /// # Arguments
    /// * `T` - Temperature in Kelvin
    ///
    /// # Returns
    /// * `Ok(())` - Closures successfully created
    /// * `Err(NISTError)` - If temperature is outside all available ranges
    pub fn create_closures_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), NISTError> {
        let is_this_T_from_current_T_range = self.is_this_T_from_current_T_range(T);
        if !is_this_T_from_current_T_range {
            self.interval_for_this_T(T)?;
            self.extract_coefficients(T)?;
        }
        self.create_closure_cp_dh_ds()?;
        Ok(())
    }

    /// Calculates Cp, dH, dS values with automatic temperature range handling
    ///
    /// Automatically selects appropriate NIST Shomate coefficients for the given temperature
    /// and calculates thermodynamic properties. If temperature is outside current range,
    /// updates coefficients before calculation.
    ///
    /// # Arguments
    /// * `T` - Temperature in Kelvin
    ///
    /// # Returns
    /// * `Ok(())` - Values successfully calculated and stored in Cp, dh, ds fields
    /// * `Err(NISTError)` - If temperature is outside all available ranges
    pub fn calculate_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), NISTError> {
        let is_this_T_from_current_T_range = self.is_this_T_from_current_T_range(T);
        if !is_this_T_from_current_T_range {
            self.interval_for_this_T(T)?;
            self.extract_coefficients(T)?;
        }
        self.calculate_cp_dh_ds(T)?;
        Ok(())
    }

    /// Creates symbolic expressions for Cp, dH, dS with automatic temperature range handling
    ///
    /// Automatically selects appropriate NIST Shomate coefficients for the given temperature
    /// and creates symbolic expressions. If temperature is outside current range,
    /// updates coefficients before creating expressions.
    ///
    /// # Arguments
    /// * `T` - Temperature in Kelvin
    ///
    /// # Returns
    /// * `Ok(())` - Symbolic expressions successfully created
    /// * `Err(NISTError)` - If temperature is outside all available ranges
    pub fn create_sym_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), NISTError> {
        let is_this_T_from_current_T_range = self.is_this_T_from_current_T_range(T);
        if !is_this_T_from_current_T_range {
            self.interval_for_this_T(T)?;
            self.extract_coefficients(T)?;
        }
        self.create_sym_cp_dh_ds()?;
        Ok(())
    }

    /// Finds the coefficient interval that contains the given temperature
    pub fn interval_for_this_T(&mut self, t: f64) -> Result<(usize, Coeffs), NISTError> {
        let intervals = self.coeffs_map.clone();

        for (interval_idx, coeffs) in intervals {
            if coeffs.T.0 <= t && t <= coeffs.T.1 {
                self.current_T_range = Some(coeffs.T.clone());
                return Ok((interval_idx, coeffs.clone()));
            }
        }

        Err(NISTError::InvalidTemperatureRange)
    }

    /// Main interface for temperature interval coefficient fitting with fallback
    pub fn fitting_coeffs_for_T_interval(&mut self) -> Result<(), NISTError> {
        if let Some((T_min, T_max)) = self.T_interval {
            let (i_t_min, coeffs_min) = self.interval_for_this_T(T_min)?;
            let (i_t_max, coeffs_max) = self.interval_for_this_T(T_max)?;

            if i_t_min == i_t_max {
                // Both boundaries are in the same interval, no fitting needed
                self.coeffs = Some(coeffs_min.coeff);
            } else if (i_t_max as isize - i_t_min as isize).abs() == 1 {
                // Boundaries are in adjacent intervals, perform fitting with fallback

                let result =
                    self.fit_cp_dh_ds(coeffs_min.clone(), coeffs_max.clone(), T_min, T_max);

                match result {
                    Ok(_) => {} // fitting_adjacent2 succeeded
                    _ => {
                        // fitting_adjacent2 failed, use weighted fitting as fallback

                        self.fitting_adjacent_weighted(coeffs_min, coeffs_max, T_min, T_max)?;
                    }
                }
            } else if (i_t_max as isize - i_t_min as isize).abs() > 1 {
                // Non-adjacent intervals, collect all coefficients in range
                let mut coeffs_in_range = Vec::new();
                for (_, coeffs) in &self.coeffs_map {
                    if coeffs.T.1 > T_min && coeffs.T.0 < T_max {
                        coeffs_in_range.push(coeffs.clone());
                    }
                }
                let result = self.fitting_non_adjacent(coeffs_in_range.clone(), T_min, T_max);
                match result {
                    Ok(_) => {} // fitting_non_adjacent succeeded
                    _ => {
                        self.fitting_cp_dh_ds_non_adjacent(coeffs_in_range, T_min, T_max)?;
                    }
                }
            } else {
                return Err(NISTError::InvalidTemperatureRange);
            }
            Ok(())
        } else {
            Err(NISTError::InvalidTemperatureRange)
        }
    }
}

////////////////////////////////////////NIST FORMAT FUNCTIONS//////////////////////////////////////////////////////////////////////////
pub fn calculate_cp(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    a + b * t + c * t.powi(2) + d * t.powi(3) + e / t.powi(2)
}

pub fn calculate_dh(
    t: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    _g: f64,
    h: f64,
) -> f64 {
    a * t + (b * t.powi(2)) / 2.0 + (c * t.powi(3)) / 3.0 + (d * t.powi(4)) / 4.0 - e / t + f - h
}

pub fn calculate_s(
    t: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    _f: f64,
    g: f64,
    _h: f64,
) -> f64 {
    a * t.ln() + b * t + (c * t.powi(2)) / 2.0 + (d * t.powi(3)) / 3.0 - e / (2.0 * t.powi(2)) + g
}

pub fn calculate_cp_sym(a: f64, b: f64, c: f64, d: f64, e: f64) -> Expr {
    let T = Expr::Var("T".to_owned());
    let t = T / Expr::Const(1000.0);
    Expr::Const(a)
        + Expr::Const(b) * t.clone()
        + Expr::Const(c) * t.clone().pow(e2)
        + Expr::Const(d) * t.clone().pow(e3)
        + Expr::Const(e) / t.pow(e2)
}

pub fn calculate_dh_sym(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64, _g: f64, h: f64) -> Expr {
    let T = Expr::Var("T".to_owned());
    let t = T / Expr::Const(1000.0);
    Expr::Const(a) * t.clone()
        + (Expr::Const(b) * t.clone().pow(e2)) / e2
        + (Expr::Const(c) * t.clone().pow(e3)) / e3
        + (Expr::Const(d) * t.clone().pow(e4)) / e4
        - Expr::Const(e) / t.clone()
        + Expr::Const(f)
        - Expr::Const(h)
}

pub fn calculate_s_sym(a: f64, b: f64, c: f64, d: f64, e: f64, _f: f64, g: f64, _h: f64) -> Expr {
    let T = Expr::Var("T".to_owned());
    let t = T / Expr::Const(1000.0);
    Expr::Const(a) * t.clone().ln()
        + Expr::Const(b) * t.clone()
        + (Expr::Const(c) * t.clone().pow(e2)) / e2
        + (Expr::Const(d) * t.clone().pow(e3)) / e3
        - Expr::Const(e) / (e2 * t.pow(e2))
        + Expr::Const(g)
}

impl Clone for NISTdata {
    fn clone(&self) -> Self {
        // Clone the necessary data upfront
        // Clone all necessary data upfront

        Self {
            input: self.input.clone(),
            unit: self.unit.clone(),
            unit_multiplier: self.unit_multiplier,
            coeffs_map: self.coeffs_map.clone(),
            coeffs: self.coeffs.clone(),
            Cp: self.Cp,
            dh: self.dh,
            ds: self.ds,
            C_fun: {
                let T_ranges = self.input.T.clone().unwrap();
                let cp_coeffs = self.input.cp.clone().unwrap();
                let unit_multiplier = self.unit_multiplier;

                let coeffs = move |T: f64| -> Vec<f64> {
                    for (i, T_pairs) in T_ranges.iter().enumerate() {
                        if T >= T_pairs[0] && T <= T_pairs[1] {
                            return cp_coeffs[i].clone();
                        }
                    }
                    vec![] // Return an empty vector if no matching range is found
                };

                Box::new(move |T: f64| {
                    let coeff = coeffs(T);
                    if coeff.len() >= 5 {
                        let (a, b, c, d, e) = (coeff[0], coeff[1], coeff[2], coeff[3], coeff[4]);
                        unit_multiplier * calculate_cp(T / 1000.0, a, b, c, d, e)
                    } else {
                        0.0 // Return a default value if coefficients are not available
                    }
                })
            },
            dh_fun: {
                let T_ranges = self.input.T.clone().unwrap();
                let cp_coeffs = self.input.cp.clone().unwrap();
                let unit_multiplier = self.unit_multiplier;
                let dh0 = self.input.dh.unwrap_or(0.0); // Use 0.0 for simple substances

                let coeffs = move |T: f64| -> Vec<f64> {
                    for (i, T_pairs) in T_ranges.iter().enumerate() {
                        if T >= T_pairs[0] && T <= T_pairs[1] {
                            return cp_coeffs[i].clone();
                        }
                    }
                    vec![] // Return an empty vector if no matching range is found
                };

                Box::new(move |T: f64| {
                    let coeff = coeffs(T);
                    if coeff.len() >= 8 {
                        let (a, b, c, d, e, f, g, h) = (
                            coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5], coeff[6],
                            coeff[7],
                        );
                        1000.0
                            * unit_multiplier
                            * (calculate_dh(T / 1000.0, a, b, c, d, e, f, g, h) + dh0)
                    } else {
                        0.0 // Return a default value if coefficients are not available
                    }
                })
            },
            ds_fun: {
                let T_ranges = self.input.T.clone().unwrap();
                let cp_coeffs = self.input.cp.clone().unwrap();
                let unit_multiplier = self.unit_multiplier;

                let coeffs = move |T: f64| -> Vec<f64> {
                    for (i, T_pairs) in T_ranges.iter().enumerate() {
                        if T >= T_pairs[0] && T <= T_pairs[1] {
                            return cp_coeffs[i].clone();
                        }
                    }
                    vec![] // Return an empty vector if no matching range is found
                };

                Box::new(move |T: f64| {
                    let coeff = coeffs(T);
                    if coeff.len() >= 8 {
                        let (a, b, c, d, e, f, g, h) = (
                            coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5], coeff[6],
                            coeff[7],
                        );
                        unit_multiplier * calculate_s(T / 1000.0, a, b, c, d, e, f, g, h)
                    } else {
                        0.0 // Return a default value if coefficients are not available 
                    }
                })
            },
            Cp_sym: self.Cp_sym.clone(),
            dh_sym: self.dh_sym.clone(),
            ds_sym: self.ds_sym.clone(),
            T_interval: self.T_interval.clone(),
            current_T_range: self.current_T_range.clone(),
        }
    }
}
impl Debug for NISTdata {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("NISTdata")
            .field("input", &self.input)
            .field("unit", &self.unit)
            .field("Cp", &self.Cp)
            .field("dh", &self.dh)
            .field("ds", &self.ds)
            .field("Cp_sym", &self.Cp_sym)
            .field("dh_sym", &self.dh_sym)
            .field("ds_sym", &self.ds_sym)
            .finish_non_exhaustive() // Indicates that not all fields are being displayed
    }
}
use super::thermo_api::{EnergyUnit, ThermoCalculator, ThermoError, energy_dimension};
impl ThermoCalculator for NISTdata {
    fn newinstance(&mut self) -> Result<(), ThermoError> {
        *self = NISTdata::new();
        Ok(())
    }
    fn extract_model_coefficients(&mut self, t: f64) -> Result<(), ThermoError> {
        self.extract_coefficients(t)?;
        Ok(())
    }
    fn parse_coefficients(&mut self) -> Result<(), ThermoError> {
        self.parse_coefficients()?;
        Ok(())
    }
    fn set_unit(&mut self, unit: Option<EnergyUnit>) -> Result<(), ThermoError> {
        if let Some(unit) = unit {
            self.set_unit(&energy_dimension(unit))
                .map_err(|e| ThermoError::UnsupportedUnit(e.to_string()))?;
        }
        Ok(())
    }
    fn from_serde(&mut self, serde: Value) -> Result<(), ThermoError> {
        self.from_serde(serde)?;
        Ok(())
    }
    fn calculate_Cp_dH_dS(&mut self, temperature: f64) -> Result<(), ThermoError> {
        self.calculate_cp_dh_ds(temperature)
            .map_err(|e| ThermoError::CalculationError(e.to_string()))
    }
    fn create_closures_Cp_dH_dS(&mut self) -> Result<(), ThermoError> {
        self.create_closure_cp_dh_ds()
            .map_err(|e| ThermoError::CalculationError(e.to_string()))
    }
    fn create_sym_Cp_dH_dS(&mut self) -> Result<(), ThermoError> {
        self.create_sym_cp_dh_ds()
            .map_err(|e| ThermoError::CalculationError(e.to_string()))
    }
    fn Taylor_series_cp_dh_ds(
        &mut self,
        temperature: f64,
        order: usize,
    ) -> Result<(Expr, Expr, Expr), ThermoError> {
        let (Cp, dh, ds) = self.clone().Taylor_series_Cp_dH_dS(temperature, order)?;
        Ok((Cp, dh, ds))
    }
    fn pretty_print_data(&self) -> Result<(), ThermoError> {
        self.pretty_print()?;
        Ok(())
    }
    fn renew_base(
        &mut self,
        sub_name: String,
        search_type: SearchType,
        phase: Phase,
    ) -> Result<(), ThermoError> {
        self.get_data_from_NIST(sub_name, search_type, phase)?;
        println!(
            "Cp found {:?}\n for temperature ranges {:?}\n",
            self.input.cp, self.input.T
        );

        Ok(())
    }
    fn get_coefficients(&self) -> Result<Vec<f64>, ThermoError> {
        let (a, b, c, d, e, f, g, h) = self.coeffs.unwrap();
        Ok(vec![a, b, c, d, e, f, g, h])
    }

    fn print_instance(&self) -> Result<(), ThermoError> {
        println!("{:?}", &self);
        Ok(())
    }
    fn get_Cp(&self) -> Result<f64, ThermoError> {
        Ok(self.Cp)
    }

    fn get_dh(&self) -> Result<f64, ThermoError> {
        Ok(self.dh)
    }

    fn get_ds(&self) -> Result<f64, ThermoError> {
        Ok(self.ds)
    }

    fn get_C_fun(&self) -> Result<Box<dyn Fn(f64) -> f64>, ThermoError> {
        Ok(self.clone().C_fun)
    }

    fn get_dh_fun(&self) -> Result<Box<dyn Fn(f64) -> f64>, ThermoError> {
        Ok(self.clone().dh_fun)
    }

    fn get_ds_fun(&self) -> Result<Box<dyn Fn(f64) -> f64>, ThermoError> {
        Ok(self.clone().ds_fun)
    }

    fn get_Cp_sym(&self) -> Result<Expr, ThermoError> {
        Ok(self.Cp_sym.clone())
    }

    fn get_dh_sym(&self) -> Result<Expr, ThermoError> {
        Ok(self.dh_sym.clone())
    }

    fn get_ds_sym(&self) -> Result<Expr, ThermoError> {
        Ok(self.ds_sym.clone())
    }
    fn get_composition(&self) -> Result<Option<HashMap<String, f64>>, ThermoError> {
        Ok(None)
    }
    fn fitting_coeffs_for_T_interval(&mut self) -> Result<(), ThermoError> {
        self.fitting_coeffs_for_T_interval()
            .map_err(|e| ThermoError::CalculationError(e.to_string()))
    }
    fn integr_mean(&mut self) -> Result<(), ThermoError> {
        self.integr_mean()
            .map_err(|e| ThermoError::CalculationError(e.to_string()))
    }

    fn set_T_interval(&mut self, T_min: f64, T_max: f64) -> Result<(), ThermoError> {
        self.set_T_interval(T_min, T_max);
        Ok(())
    }
    fn calculate_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), ThermoError> {
        self.calculate_Cp_dH_dS_with_T_range(T)?;
        Ok(())
    }
    fn create_closures_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), ThermoError> {
        self.create_closures_Cp_dH_dS_with_T_range(T)?;
        Ok(())
    }
    fn create_sym_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), ThermoError> {
        self.create_sym_Cp_dH_dS_with_T_range(T)?;
        Ok(())
    }
    fn is_coeffs_valid_for_T(&self, T: f64) -> Result<bool, ThermoError> {
        let flag = self.is_this_T_from_current_T_range(T);
        Ok(flag)
    }
}
