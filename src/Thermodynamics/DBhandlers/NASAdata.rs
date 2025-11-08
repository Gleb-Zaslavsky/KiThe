//! # NASA Thermodynamic Data Handler Module
//!
//! ## Aim
//! This module handles NASA polynomial format thermodynamic data for heat capacity (Cp),
//! enthalpy (H), and entropy (S) calculations. It processes NASA-7 coefficient format
//! with temperature-dependent polynomials and provides numerical, functional, and symbolic computations.
//!
//! ## Main Data Structures and Logic
//! - `NASAdata`: Main structure storing NASA-7 coefficients and calculated thermodynamic properties
//! - `NASAinput`: Input structure for parsing NASA database entries with composition and coefficients
//! - `Coeffs`: Structure holding temperature range and 7 NASA coefficients for each interval
//! - `FittingReport`: Quality metrics for curve fitting operations
//! - Uses NASA-7 polynomial format: Cp/R = a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴
//! - Supports multiple temperature ranges (typically 2-3 ranges per substance)
//! - Handles both gas and condensed phase data
//!
//! ## Key Methods
//! ### Basic Operations
//! - `extract_coefficients(T)`: Selects appropriate coefficient set for given temperature
//! - `calculate_Cp_dH_dS(T)`: Computes thermodynamic properties at specified temperature
//! - `create_closures_Cp_dH_dS()`: Creates closure functions for efficient repeated calculations
//! - `create_sym_Cp_dH_dS()`: Creates symbolic expressions for analytical work
//! - `Taylor_series_Cp_dH_dS()`: Generates Taylor series expansions around specified temperature
//! - `pretty_print()`: Displays coefficient data in formatted table
//!
//! ### Temperature Range Fitting
//! - `interval_for_this_T(T)`: Finds coefficient interval containing given temperature
//! - `direct_compare()`: Compares fitted vs original functions and generates quality metrics
//! - `fitting_adjacent()`: Sequential fitting of Cp, dH, dS for adjacent temperature intervals
//! - `fitting_adjacent2()`: Improved fitting starting with entropy, includes quality reports
//! - `fitting_adjacent3()`: Weighted sum fitting (Y = Cp + a*dH + b*dS) for simultaneous optimization
//! - `fitting_coeffs_for_T_interval()`: Main interface for temperature interval fitting with fallback
//!
//! ## Usage
//! ```rust, ignore
//! let mut nasa = NASAdata::new();
//! nasa.from_serde(database_entry)?;
//! nasa.set_unit("J")?;  // or "cal" for calories
//! nasa.set_T_interval(400.0, 800.0);  // Set temperature range for fitting
//! nasa.fitting_coeffs_for_T_interval()?;  // Fit coefficients across temperature range
//! nasa.calculate_Cp_dH_dS(600.0);
//! let cp = nasa.Cp;  // Heat capacity
//! let dh = nasa.dh;  // Enthalpy
//! let ds = nasa.ds;  // Entropy
//! ```
//!
//! ## Advanced Features
//! - Handles complex composition parsing from string format ("{C:1, H:4}")
//! - Supports both Joules and calories with automatic unit conversion
//! - Implements NASA-7 polynomial integration for enthalpy and entropy calculations
//! - Provides Taylor series expansion for local approximations
//! - Uses custom deserialization for composition data parsing
//! - Implements comprehensive error handling for temperature range validation
//! - Creates both numerical functions and symbolic expressions for flexible usage
//! - Panic-safe fitting with automatic fallback mechanisms
//! - Quality assessment with L1, L2, and max norms for fitting validation
//! - Pretty table reporting for fitting quality metrics

use ::RustedSciThe::numerical::optimization::sym_fitting::Fitting;
use RustedSciThe::numerical::optimization::fitting_features::SewTwoFunctions;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use prettytable::{Cell, Row, Table};
use serde::de::Deserializer;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;

use std::error::Error;
use std::fmt;
use std::panic;
use std::str::FromStr;
#[allow(non_upper_case_globals)]
const R: f64 = 1.987; // кал/(K·моль)
#[allow(non_upper_case_globals)]
const Rsym: Expr = Expr::Const(1.987);
use crate::Thermodynamics::DBhandlers::NIST_parser::{Phase, SearchType};
#[derive(Debug)]
pub enum NASAError {
    NoCoefficientsFound { temperature: f64, range: String },
    InvalidTemperatureRange,
    SerdeError(serde_json::Error),
    UnsupportedUnit(String),
}

impl fmt::Display for NASAError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            NASAError::NoCoefficientsFound { temperature, range } => {
                write!(
                    f,
                    "No coefficients found for temperature {} K. Valid range: {}",
                    temperature, range
                )
            }
            NASAError::InvalidTemperatureRange => {
                write!(f, "Invalid temperature range in coefficient data")
            }
            NASAError::SerdeError(msg) => {
                write!(f, "Failed to deserialize NASA data: {}", msg)
            }
            NASAError::UnsupportedUnit(unit) => {
                write!(
                    f,
                    "Unsupported unit: {}. Only 'J' and 'cal' are supported",
                    unit
                )
            }
        }
    }
}

impl Error for NASAError {}

impl From<serde_json::Error> for NASAError {
    fn from(err: serde_json::Error) -> Self {
        NASAError::SerdeError(err)
    }
}

pub fn Cp(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    R * (a + b * t + c * t.powi(2) + d * t.powi(3) + e * t.powi(4))
}
pub fn dh(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> f64 {
    R * t
        * (a + b * t / 2.0
            + c * t.powi(2) / 3.0
            + d * t.powi(3) / 4.0
            + e * t.powi(4) / 5.0
            + f / t)
}
pub fn ds(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64, g: f64) -> f64 {
    R * (a * t.ln() + b * t + c * t.powi(2) / 2.0 + d * t.powi(3) / 3.0 + e * t.powi(4) / 4.0 + g)
}
pub fn Cp_sym(a: f64, b: f64, c: f64, d: f64, e: f64) -> Expr {
    let t = Expr::Var("T".to_owned());
    let (a, b, c, d, e) = (
        Expr::Const(a),
        Expr::Const(b),
        Expr::Const(c),
        Expr::Const(d),
        Expr::Const(e),
    );
    Rsym * (a
        + b * t.clone()
        + c * t.clone().pow(Expr::Const(2.0))
        + d * t.clone().pow(Expr::Const(3.0))
        + e * t.clone().pow(Expr::Const(4.0)))
}
fn dh_sym(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> Expr {
    let t = Expr::Var("T".to_owned());
    let (a, b, c, d, e, f) = (
        Expr::Const(a),
        Expr::Const(b),
        Expr::Const(c),
        Expr::Const(d),
        Expr::Const(e),
        Expr::Const(f),
    );
    Rsym * t.clone()
        * (a + b * t.clone() / Expr::Const(2.0)
            + c * t.clone().pow(Expr::Const(2.0)) / Expr::Const(3.0)
            + d * t.clone().pow(Expr::Const(3.0)) / Expr::Const(4.0)
            + e * t.clone().pow(Expr::Const(4.0)) / Expr::Const(5.0)
            + f / t)
}
fn ds_sym(a: f64, b: f64, c: f64, d: f64, e: f64, g: f64) -> Expr {
    let t = Expr::Var("T".to_owned());
    let (a, b, c, d, e, g) = (
        Expr::Const(a),
        Expr::Const(b),
        Expr::Const(c),
        Expr::Const(d),
        Expr::Const(e),
        Expr::Const(g),
    );
    Rsym * (a * t.clone().ln()
        + b * t.clone()
        + c * t.clone().pow(Expr::Const(2.0)) / Expr::Const(2.0)
        + d * t.clone().pow(Expr::Const(3.0)) / Expr::Const(3.0)
        + e * t.clone().pow(Expr::Const(4.0)) / Expr::Const(4.0)
        + g)
}
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct NASAinput {
    #[serde(deserialize_with = "deserialize_composition")]
    pub composition: Option<HashMap<String, f64>>,
    pub Cp: Vec<f64>,

    pub model: String,
}

fn deserialize_composition<'de, D>(
    deserializer: D,
) -> Result<Option<HashMap<String, f64>>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    if s.is_empty() {
        return Ok(None);
    }

    let s = s.trim_matches(|c| c == '{' || c == '}');
    let pairs: Result<HashMap<String, f64>, _> = s
        .split(',')
        .map(|pair| {
            let mut parts = pair.trim().split(':');
            let key = parts
                .next()
                .ok_or("Missing key")
                .unwrap()
                .trim()
                .to_string();
            let value = f64::from_str(parts.next().ok_or("Missing value").unwrap().trim())
                .map_err(|_| "Invalid float")
                .unwrap();
            Ok((key, value))
        })
        .collect();

    pairs.map(Some)
}
#[derive(Debug, Clone)]
pub struct Coeffs {
    pub T: (f64, f64),
    pub coeff: (f64, f64, f64, f64, f64, f64, f64),
}

pub struct NASAdata {
    /// data parsed from library
    pub input: NASAinput,
    /// Joles or Cal
    pub unit: Option<String>,
    ///
    pub unit_multiplier: f64,
    /// NASA format 7 coefficients
    pub coeffs: (f64, f64, f64, f64, f64, f64, f64),

    pub coeffs_map: HashMap<usize, Coeffs>,
    /// heat capacity value at T
    pub Cp: f64,
    /// enthalpy value at T
    pub dh: f64,
    /// entropy value at T
    pub ds: f64,
    /// heat capacity function
    #[allow(clippy::type_complexity)]
    pub C_fun: Box<dyn Fn(f64) -> f64>,
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
    pub norm_threshold: f64,
    pub T_interval: Option<(f64, f64)>,
    pub fit_point_number: usize,
}

impl NASAdata {
    ///initialize fields
    pub fn new() -> Self {
        let input = NASAinput {
            Cp: Vec::new(),
            composition: None,
            model: String::new(),
        };
        Self {
            input: input,
            unit: None,
            unit_multiplier: 4.184,
            coeffs: (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            coeffs_map: HashMap::new(),
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
            fit_point_number: 1000,
            norm_threshold: 0.1,
        }
    }
    /// set energy unitsЖ J or calories
    pub fn set_unit(&mut self, unit: &str) -> Result<(), NASAError> {
        match unit {
            "J" => {
                self.unit = Some("J".to_string());
                self.unit_multiplier = 4.184;
                Ok(())
            }
            "cal" => {
                self.unit = Some("cal".to_string());
                self.unit_multiplier = 1.0;
                Ok(())
            }
            _ => Err(NASAError::UnsupportedUnit(unit.to_string())),
        }
    }

    pub fn set_T_interval(&mut self, T_min: f64, T_max: f64) {
        self.T_interval = Some((T_min, T_max));
    }
    /// takes serde Value and parse it into structure
    pub fn from_serde(&mut self, serde: Value) -> Result<(), NASAError> {
        self.input = serde_json::from_value(serde).map_err(|e| NASAError::SerdeError(e))?;
        Ok(())
    }

    pub fn parse_coefficients(&mut self) -> Result<HashMap<usize, Coeffs>, NASAError> {
        let c_data = self.input.Cp.as_slice();
        let mut coeffs_map = HashMap::new();

        match c_data.len() {
            25 => {
                let (t1, t2, t3, t4) = (c_data[0], c_data[1], c_data[2], c_data[3]);

                coeffs_map.insert(
                    0,
                    Coeffs {
                        T: (t1, t2),
                        coeff: (
                            c_data[4], c_data[5], c_data[6], c_data[7], c_data[8], c_data[9],
                            c_data[10],
                        ),
                    },
                );

                coeffs_map.insert(
                    1,
                    Coeffs {
                        T: (t2, t3),
                        coeff: (
                            c_data[11], c_data[12], c_data[13], c_data[14], c_data[15], c_data[16],
                            c_data[17],
                        ),
                    },
                );

                coeffs_map.insert(
                    2,
                    Coeffs {
                        T: (t3, t4),
                        coeff: (
                            c_data[18], c_data[19], c_data[20], c_data[21], c_data[22], c_data[23],
                            c_data[24],
                        ),
                    },
                );
            }
            17 => {
                let (t1, t2, t3) = (c_data[0], c_data[1], c_data[2]);

                coeffs_map.insert(
                    0,
                    Coeffs {
                        T: (t1, t2),
                        coeff: (
                            c_data[3], c_data[4], c_data[5], c_data[6], c_data[7], c_data[8],
                            c_data[9],
                        ),
                    },
                );

                coeffs_map.insert(
                    1,
                    Coeffs {
                        T: (t2, t3),
                        coeff: (
                            c_data[10], c_data[11], c_data[12], c_data[13], c_data[14], c_data[15],
                            c_data[16],
                        ),
                    },
                );
            }
            9 => {
                let (t1, t2) = (c_data[0], c_data[1]);

                coeffs_map.insert(
                    0,
                    Coeffs {
                        T: (t1, t2),
                        coeff: (
                            c_data[2], c_data[3], c_data[4], c_data[5], c_data[6], c_data[7],
                            c_data[8],
                        ),
                    },
                );
            }
            _ => return Err(NASAError::InvalidTemperatureRange),
        }
        self.coeffs_map = coeffs_map.clone();
        Ok(coeffs_map)
    }

    fn find_coefficients_for_temperature(
        &self,
        t: f64,
    ) -> Result<(f64, f64, f64, f64, f64, f64, f64), NASAError> {
        for coeffs in self.coeffs_map.values() {
            if coeffs.T.0 <= t && t <= coeffs.T.1 {
                return Ok(coeffs.coeff);
            }
        }

        let ranges: Vec<String> = self
            .coeffs_map
            .values()
            .map(|c| format!("{}-{}", c.T.0, c.T.1))
            .collect();

        Err(NASAError::NoCoefficientsFound {
            temperature: t,
            range: ranges.join(", "),
        })
    }

    /// get the 7 constants of NASA7 format for concrete temperature
    pub fn extract_coefficients(&mut self, t: f64) -> Result<(), NASAError> {
        if self.coeffs_map.is_empty() {
            self.parse_coefficients()?;
        }
        self.coeffs = self.find_coefficients_for_temperature(t)?;
        Ok(())
    }
    /// create functions for heat capacity, enthalpy, and entropy using the NASA7 format
    pub fn create_closures_Cp_dH_dS(&mut self) {
        // calculate Cp, dh, ds using the coefficients and the temperature t
        let (a, b, c, d, e, f, g) = self.coeffs;
        let unit = self.unit_multiplier;
        let Cp =
            move |t: f64| unit * R * (a + b * t + c * t.powi(2) + d * t.powi(3) + e * t.powi(4));
        self.C_fun = Box::new(Cp);
        let dh = move |t: f64| {
            unit * R
                * t
                * (a + b * t / 2.0
                    + c * t.powi(2) / 3.0
                    + d * t.powi(3) / 4.0
                    + e * t.powi(4) / 5.0
                    + f / t)
        };
        self.dh_fun = Box::new(dh);
        let ds = move |t: f64| {
            unit * R
                * (a * t.ln()
                    + b * t
                    + c * t.powi(2) / 2.0
                    + d * t.powi(3) / 3.0
                    + e * t.powi(4) / 4.0
                    + g)
        };
        self.ds_fun = Box::new(ds);
    }
    /// create symbolic expressions for heat capacity, enthalpy, and entropy
    pub fn create_sym_Cp_dH_dS(&mut self) {
        let (a, b, c, d, e, f, g) = self.coeffs;
        let unit = Expr::Const(self.unit_multiplier);
        let C_sym = unit.clone() * Cp_sym(a, b, c, d, e);
        self.Cp_sym = C_sym.simplify();
        let dh_sym = unit.clone() * dh_sym(a, b, c, d, e, f);
        self.dh_sym = dh_sym.simplify();
        let ds_sym = unit * ds_sym(a, b, c, d, e, g);
        self.ds_sym = ds_sym.simplify();
    }
    // Taylor series expension for Cp, dH, dS
    pub fn Taylor_series_Cp_dH_dS(&mut self, T0: f64, n: usize) -> (Expr, Expr, Expr) {
        let Cp = self.Cp_sym.clone();

        let Cp_taylor = Cp.taylor_series1D("T", T0, n);

        let dh = self.dh_sym.clone();

        let dh_taylor = dh.taylor_series1D("T", T0, n);

        let ds = self.ds_sym.clone();

        let ds_taylor = ds.taylor_series1D("T", T0, n);
        (Cp_taylor, dh_taylor, ds_taylor)
    }
    /// calculate heat capacity, enthalpy, and entropy as a function of given temperature using the NASA7 format
    pub fn calculate_Cp_dH_dS(&mut self, t: f64) {
        let (a, b, c, d, e, f, g) = self.coeffs;
        let unit = self.unit_multiplier;
        self.Cp = unit * Cp(t, a, b, c, d, e);
        self.dh = unit * dh(t, a, b, c, d, e, f);
        self.ds = unit * ds(t, a, b, c, d, e, g);
    }

    pub fn pretty_print(&self) {
        let data = self.input.Cp.clone();
        let len = data.len();
        // take temperatrure ranges from the data
        let T: Vec<f64> = if len == 25 {
            data.clone().into_iter().take(4).collect()
        } else if len == 17 {
            data.clone().into_iter().take(3).collect()
        } else if len == 9 {
            data.clone().into_iter().take(2).collect()
        } else {
            vec![]
        };
        // tale coefficients of the polynomial
        let coeffs: Vec<Vec<f64>> = if len == 25 {
            data[4..].chunks(7).map(|chunk| chunk.to_vec()).collect()
        } else if len == 17 {
            data[3..].chunks(7).map(|chunk| chunk.to_vec()).collect()
        } else if len == 9 {
            data[2..].chunks(7).map(|chunk| chunk.to_vec()).collect()
        } else {
            vec![]
        };

        let mut table = Table::new();
        let mut header_row = vec![Cell::new("Coefficients")];
        for pairs in T.windows(2) {
            let Tstr = format!("{:} - {:}   ", pairs[0], pairs[1]);
            header_row.push(Cell::new(&Tstr));
        }
        table.add_row(Row::new(header_row));

        let coeffs_names = vec!["A", "B", "C", "D", "E", "F", "G"];
        for (i, coeff_name) in coeffs_names.iter().enumerate() {
            let mut row = vec![Cell::new(coeff_name)];

            for coeff_for_every_T in coeffs.iter() {
                let coeff = coeff_for_every_T[i].to_string();
                row.push(Cell::new(&format!("{:}", coeff)));
            }
            table.add_row(Row::new(row));
        }

        table.printstd();
    }

    /////////////////////////////////TEMPERATURE RANGE////////////////////////////////////////////////////////////

    /// Finds the coefficient interval that contains the given temperature
    ///
    /// # Arguments
    /// * `t` - Temperature in Kelvin
    ///
    /// # Returns
    /// * `Ok((interval_index, coeffs))` - Index and coefficients for the temperature range
    /// * `Err(NASAError::InvalidTemperatureRange)` - If temperature is outside all ranges
    pub fn interval_for_this_T(&self, t: f64) -> Result<(usize, Coeffs), NASAError> {
        let intervals = self.coeffs_map.clone();

        for (interval_idx, coeffs) in intervals {
            if coeffs.T.0 <= t && t <= coeffs.T.1 {
                return Ok((interval_idx, coeffs.clone()));
            }
        }

        Err(NASAError::InvalidTemperatureRange)
    }

    /// Compares fitted function against original functions from two temperature ranges
    ///
    /// Generates quality metrics by evaluating both original and fitted functions
    /// across the temperature interval and computing various error norms.
    ///
    /// # Arguments
    /// * `func_for_1_range` - Original function for lower temperature range
    /// * `func_for_2_range` - Original function for upper temperature range  
    /// * `T_center` - Temperature boundary between the two ranges
    /// * `fit_function` - Fitted function to compare against originals
    /// * `coeffs_map` - Fitted coefficients to substitute into fit_function
    ///
    /// # Returns
    /// * `FittingReport` - Contains L1/L2/max norms and significant deviation points
    ///
    /// # Panics
    /// * If max_norm exceeds norm_threshold (indicates poor fitting quality)
    pub fn direct_compare(
        &self,
        func_for_1_range: Expr,
        func_for_2_range: Expr,
        T_center: f64,
        fit_function: Expr,
        coeffs_map: HashMap<String, f64>,
    ) -> FittingReport {
        if let Some((T_min, T_max)) = self.T_interval {
            let n = self.fit_point_number / 2;
            let step = (T_max - T_min) / (2.0 * n as f64 - 1.0);
            let T_vec: Vec<f64> = (0..2 * n).map(|i| T_min + i as f64 * step).collect();
            // Generate original data from two ranges
            let mut y_data = func_for_1_range.lambdify1D_from_linspace(T_min, T_center, n);
            let y_data2 = func_for_2_range.lambdify1D_from_linspace(T_center, T_max, n);
            y_data.extend(y_data2);

            // Generate fitted data
            let fitted_data = fit_function.set_variable_from_map(&coeffs_map);
            let y_data_fitted = fitted_data.lambdify1D_from_linspace(T_min, T_max, 2 * n);

            // Calculate norms
            let mut l1_norm = 0.0;
            let mut l2_norm = 0.0;
            let mut max_norm: f64 = 0.0;
            let mut map_of_residuals: HashMap<usize, (f64, f64)> = HashMap::new();
            assert_eq!(
                y_data.len(),
                y_data_fitted.len(),
                "Vectors must have the same length"
            );
            assert_eq!(y_data.len(), 2 * n, "Vectors must have the same length");
            for (i, (y, y_fit)) in y_data.iter().zip(y_data_fitted.iter()).enumerate() {
                let diff = (y - y_fit).abs();

                l1_norm += diff;
                l2_norm += diff * diff;
                let rel_diff = diff / y;
                if rel_diff.abs() > 10e-2 {
                    map_of_residuals.insert(i, (T_vec[i], rel_diff));
                }
                max_norm = max_norm.max((diff / y).abs());
                assert!(
                    !(max_norm > self.norm_threshold),
                    "max norm too higher  {} than therhold {}, diff = {}, y = {}",
                    max_norm,
                    self.norm_threshold,
                    diff,
                    y
                );
            }

            l1_norm = l1_norm / (2.0 * n as f64);
            l2_norm = l2_norm.sqrt() / (2.0 * n as f64);
            let number_of_significant_points = map_of_residuals.len();
            let mut T_range: Option<(f64, f64)> = None;
            if map_of_residuals.len() > 0 {
                let max_T = *map_of_residuals
                    .iter()
                    .map(|(_, (T, _))| T)
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();
                let min_T = *map_of_residuals
                    .iter()
                    .map(|(_, (T, _))| T)
                    .min_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();
                T_range = Some((min_T, max_T));
                println!(
                    "fitted curve significantly deviated from directly calculated in {} points from {} to {} K",
                    number_of_significant_points, min_T, max_T
                );
            }

            println!(
                "L1 norm: {}, L2 norm: {}, Max norm: {}",
                l1_norm, l2_norm, max_norm
            );

            FittingReport {
                l2_norm: l2_norm,
                max_norm: max_norm,
                number_of_significant_points: number_of_significant_points,
                T_range: T_range,
            }
        } else {
            FittingReport {
                l2_norm: 0.0,
                max_norm: 0.0,
                number_of_significant_points: 0,
                T_range: None,
            }
        }
    }

    /// Sequential fitting method for adjacent temperature intervals
    ///
    /// Fits NASA-7 coefficients by sequentially optimizing Cp, then dH, then dS.
    /// Uses the fitted Cp coefficients as constraints for dH and dS fitting.
    /// Includes quality assessment and generates fitting reports.
    ///
    /// # Arguments
    /// * `coeffs_min` - Coefficients for lower temperature range
    /// * `coeffs_max` - Coefficients for upper temperature range
    /// * `T_min` - Lower bound of fitting temperature range
    /// * `T_max` - Upper bound of fitting temperature range
    ///
    /// # Returns
    /// * `Ok(SewTwoFunctions)` - Fitting object with results
    /// * `Err(NASAError)` - If temperature ranges are not adjacent or fitting fails
    pub fn fitting_adjacent(
        &mut self,
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<SewTwoFunctions, NASAError> {
        // Generate fitting reports
        let mut report_map: HashMap<String, FittingReport> = HashMap::new();
        let (_, T_center) = (coeffs_min.T.0, coeffs_min.T.1);
        let (T_center_check, _) = (coeffs_max.T.0, coeffs_max.T.1);
        if T_center != T_center_check {
            return Err(NASAError::InvalidTemperatureRange);
        }
        let func1 = Cp_sym(
            coeffs_min.coeff.0,
            coeffs_min.coeff.1,
            coeffs_min.coeff.2,
            coeffs_min.coeff.3,
            coeffs_min.coeff.4,
        );
        let func2 = Cp_sym(
            coeffs_max.coeff.0,
            coeffs_max.coeff.1,
            coeffs_max.coeff.2,
            coeffs_max.coeff.3,
            coeffs_max.coeff.4,
        );
        let Cp_func_to_fit =
            Expr::parse_expression("1.987* (a + b * t + c * t^2 + d * t^3 + e * t^4)");
        let mut sew = SewTwoFunctions::new(
            func1.clone(),
            func2.clone(),
            T_min,
            T_center,
            T_max,
            self.fit_point_number,
        );
        sew.create_fitting_data();
        let initial_guess = vec![
            (coeffs_min.coeff.0 + coeffs_max.coeff.0) / 2.0,
            (coeffs_min.coeff.1 + coeffs_max.coeff.1) / 2.0,
            (coeffs_min.coeff.2 + coeffs_max.coeff.2) / 2.0,
            (coeffs_min.coeff.3 + coeffs_max.coeff.3) / 2.0,
            (coeffs_min.coeff.4 + coeffs_max.coeff.4) / 2.0,
        ];

        sew.fit(
            Cp_func_to_fit.clone(),
            Some(vec![
                "a".to_string(),
                "b".to_string(),
                "c".to_string(),
                "d".to_string(),
                "e".to_string(),
            ]),
            "t".to_string(),
            initial_guess,
            None,
            None,
            None,
            None,
            None,
        );

        let r_ssquared = sew.get_r_ssquared();
        println!("r_squared for Cp: {}", r_ssquared.unwrap());
        assert!(1.0 - r_ssquared.unwrap() < 1e-2);

        let map_of_solutions_5_coeffs = sew.get_map_of_solutions().unwrap();

        let c_fitting_report = self.direct_compare(
            func1,
            func2,
            T_center,
            Cp_func_to_fit,
            map_of_solutions_5_coeffs.clone(),
        );

        report_map.insert("c fitting report".to_string(), c_fitting_report);
        ///////////////////////////dh fitting///////////////////////////////
        let a = *map_of_solutions_5_coeffs.get("a").unwrap();
        let b = *map_of_solutions_5_coeffs.get("b").unwrap();
        let c = *map_of_solutions_5_coeffs.get("c").unwrap();
        let d = *map_of_solutions_5_coeffs.get("d").unwrap();
        let e = *map_of_solutions_5_coeffs.get("e").unwrap();

        // fitting for coefficient f
        let f = coeffs_min.coeff.5;
        let h_func1 = dh_sym(a, b, c, d, e, f);
        let f = coeffs_max.coeff.5;
        let h_func2 = dh_sym(a, b, c, d, e, f);
        let dh_func_to_fit = Expr::parse_expression(
            "1.987 *t* (a + b * t/2 + (c/3) * t^2 + (d/4) * t^3 + (e/5) * t^4 + f/t)",
        );
        let dh_func_to_fit = dh_func_to_fit
            .clone()
            .set_variable_from_map(&map_of_solutions_5_coeffs);
        let mut sew = SewTwoFunctions::new(
            h_func1.clone(),
            h_func2.clone(),
            T_min,
            T_center,
            T_max,
            self.fit_point_number,
        );
        sew.create_fitting_data();
        sew.fit(
            dh_func_to_fit.clone(),
            Some(vec!["f".to_string()]),
            "t".to_string(),
            vec![1.0],
            None,
            None,
            None,
            None,
            None,
        );
        let r_squared = sew.get_r_ssquared();
        println!("r_squared for h: {}", r_squared.unwrap());
        assert!(1.0 - r_squared.unwrap() < 5e-2);
        let f_map = sew.get_map_of_solutions().unwrap();
        let f = *f_map.get("f").unwrap();
        let mut map_of_solutions_6_coeffs = map_of_solutions_5_coeffs.clone();
        map_of_solutions_6_coeffs.extend(f_map);

        /*
        let h_fitting_report = self.direct_compare(
            h_func1,
            h_func2,
            T_center,
            dh_func_to_fit,
            map_of_solutions_6_coeffs.clone(),
        );
        report_map.insert("h fitting report".to_string(), h_fitting_report);
        */
        //////////////////////////dS fitting///////////////////////////
        let g = coeffs_min.coeff.6;
        let s_func1 = ds_sym(a, b, c, d, e, g);
        let g = coeffs_max.coeff.6;
        let s_func2 = ds_sym(a, b, c, d, e, g);
        let ds_func_to_fit = Expr::parse_expression(
            "1.987 * (a* ln( t ) + b * t + c * t^2 / 2 + d * t^3 / 3.0 + e * t^4 / 4.0 + g)",
        );

        let ds_func_to_fit = ds_func_to_fit
            .clone()
            .set_variable_from_map(&map_of_solutions_5_coeffs);
        let mut sew = SewTwoFunctions::new(
            s_func1.clone(),
            s_func2.clone(),
            T_min,
            T_center,
            T_max,
            self.fit_point_number,
        );
        sew.create_fitting_data();
        sew.fit(
            ds_func_to_fit.clone(),
            Some(vec!["g".to_string()]),
            "t".to_string(),
            vec![(coeffs_min.coeff.6 + coeffs_max.coeff.6) / 2.0],
            None,
            None,
            None,
            None,
            None,
        );
        let r_squared = sew.get_r_ssquared();
        println!("r_squared for s: {}", r_squared.unwrap());

        let g_map = sew.get_map_of_solutions().unwrap();
        let g = *g_map.get("g").unwrap();
        let g = if 1.0 - r_squared.unwrap() < 5e-2 {
            g
        } else {
            (coeffs_min.coeff.6 + coeffs_max.coeff.6) / 2.0
        };
        let mut map_of_solutions_7_coeffs = map_of_solutions_6_coeffs.clone();
        map_of_solutions_7_coeffs.extend(g_map);
        self.coeffs = (a, b, c, d, e, f, g);

        let s_fitting_report = self.direct_compare(
            s_func1,
            s_func2,
            T_center,
            ds_func_to_fit,
            map_of_solutions_7_coeffs,
        );
        report_map.insert("s fitting report".to_string(), s_fitting_report);

        Self::print_quality_table(&report_map);

        Ok(sew)
    }

    /// Weighted sum fitting method for simultaneous optimization
    ///
    /// Fits NASA-7 coefficients by optimizing a weighted combination:
    /// Y = Cp + 1e-3*dH + 1e-1*dS
    ///
    /// This approach optimizes all thermodynamic properties simultaneously rather than
    /// sequentially, potentially providing better overall consistency. The weights
    /// account for the different scales of the properties (dH is typically large).
    ///
    /// # Arguments
    /// * `coeffs_min` - Coefficients for lower temperature range
    /// * `coeffs_max` - Coefficients for upper temperature range
    /// * `T_min` - Lower bound of fitting temperature range
    /// * `T_max` - Upper bound of fitting temperature range
    ///
    /// # Returns
    /// * `Ok(HashMap<String, FittingReport>)` - Quality report for combined fitting
    /// * `Err(NASAError)` - If temperature ranges are not adjacent or fitting fails
    pub fn fitting_adjacent3(
        &mut self,
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NASAError> {
        let mut report_map: HashMap<String, FittingReport> = HashMap::new();
        let (_, T_center) = (coeffs_min.T.0, coeffs_min.T.1);
        let (T_center_check, _) = (coeffs_max.T.0, coeffs_max.T.1);
        if T_center != T_center_check {
            return Err(NASAError::InvalidTemperatureRange);
        }

        // Create combined function Y = Cp + a*dH + b*dS
        let cp_func1 = Cp_sym(
            coeffs_min.coeff.0,
            coeffs_min.coeff.1,
            coeffs_min.coeff.2,
            coeffs_min.coeff.3,
            coeffs_min.coeff.4,
        );
        let dh_func1 = dh_sym(
            coeffs_min.coeff.0,
            coeffs_min.coeff.1,
            coeffs_min.coeff.2,
            coeffs_min.coeff.3,
            coeffs_min.coeff.4,
            coeffs_min.coeff.5,
        );
        let ds_func1 = ds_sym(
            coeffs_min.coeff.0,
            coeffs_min.coeff.1,
            coeffs_min.coeff.2,
            coeffs_min.coeff.3,
            coeffs_min.coeff.4,
            coeffs_min.coeff.6,
        );

        let cp_func2 = Cp_sym(
            coeffs_max.coeff.0,
            coeffs_max.coeff.1,
            coeffs_max.coeff.2,
            coeffs_max.coeff.3,
            coeffs_max.coeff.4,
        );
        let dh_func2 = dh_sym(
            coeffs_max.coeff.0,
            coeffs_max.coeff.1,
            coeffs_max.coeff.2,
            coeffs_max.coeff.3,
            coeffs_max.coeff.4,
            coeffs_max.coeff.5,
        );
        let ds_func2 = ds_sym(
            coeffs_max.coeff.0,
            coeffs_max.coeff.1,
            coeffs_max.coeff.2,
            coeffs_max.coeff.3,
            coeffs_max.coeff.4,
            coeffs_max.coeff.6,
        );

        // Combined weighted function
        let combined_func1 = cp_func1.clone()
            + Expr::Const(1e-3) * dh_func1.clone()
            + Expr::Const(1e-1) * ds_func1.clone();
        let combined_func2 = cp_func2.clone()
            + Expr::Const(1e-3) * dh_func2.clone()
            + Expr::Const(1e-1) * ds_func2.clone();

        // let combined_func_to_fit = Expr::parse_expression(
        //      "1.987* (a + b * t + c * t^2 + d * t^3 + e * t^4) + 1e-3 * 1.987 *t* (a + b * t/2 + (c/3) * t^2 + (d/4) * t^3 + (e/5) * t^4 + f/t) + 1e-1 * 1.987 * (a* ln( t ) + b * t + c * t^2 / 2 + d * t^3 / 3.0 + e * t^4 / 4.0 + g)"
        //   ).simplify_();
        let Cp_func_to_fit =
            Expr::parse_expression("1.987* (a + b * t + c * t^2 + d * t^3 + e * t^4)");
        let dh_func_to_fit = Expr::parse_expression(
            "1.987 *t* (a + b * t/2 + (c/3) * t^2 + (d/4) * t^3 + (e/5) * t^4 + f/t)",
        );
        let ds_func_to_fit = Expr::parse_expression(
            "1.987 * (a* ln( t ) + b * t + c * t^2 / 2 + d * t^3 / 3.0 + e * t^4 / 4.0 + g)",
        );
        let combined_func_to_fit = Cp_func_to_fit.clone()
            + Expr::Const(1e-3) * dh_func_to_fit.clone()
            + Expr::Const(1.0) * ds_func_to_fit.clone();
        let mut sew = SewTwoFunctions::new(
            combined_func1.clone(),
            combined_func2.clone(),
            T_min,
            T_center,
            T_max,
            self.fit_point_number,
        );

        sew.create_fitting_data();
        let initial_guess = vec![
            (coeffs_min.coeff.0 + coeffs_max.coeff.0) / 2.0,
            (coeffs_min.coeff.1 + coeffs_max.coeff.1) / 2.0,
            (coeffs_min.coeff.2 + coeffs_max.coeff.2) / 2.0,
            (coeffs_min.coeff.3 + coeffs_max.coeff.3) / 2.0,
            (coeffs_min.coeff.4 + coeffs_max.coeff.4) / 2.0,
            (coeffs_min.coeff.5 + coeffs_max.coeff.5) / 2.0,
            (coeffs_min.coeff.6 + coeffs_max.coeff.6) / 2.0,
        ];

        sew.fit(
            combined_func_to_fit.clone(),
            Some(vec![
                "a".to_string(),
                "b".to_string(),
                "c".to_string(),
                "d".to_string(),
                "e".to_string(),
                "f".to_string(),
                "g".to_string(),
            ]),
            "t".to_string(),
            initial_guess,
            None,
            None,
            None,
            None,
            None,
        );

        let r_squared = sew.get_r_ssquared();
        println!("r_squared for combined fitting: {}", r_squared.unwrap());

        let map_of_solutions = sew.get_map_of_solutions().unwrap();
        let a = *map_of_solutions.get("a").unwrap();
        let b = *map_of_solutions.get("b").unwrap();
        let c = *map_of_solutions.get("c").unwrap();
        let d = *map_of_solutions.get("d").unwrap();
        let e = *map_of_solutions.get("e").unwrap();
        let f = *map_of_solutions.get("f").unwrap();
        let g = *map_of_solutions.get("g").unwrap();

        self.coeffs = (a, b, c, d, e, f, g);

        // Generate fitting report for combined function
        let combined_fitting_report = self.direct_compare(
            combined_func1,
            combined_func2,
            T_center,
            combined_func_to_fit,
            map_of_solutions,
        );
        report_map.insert(
            "Combined fitting report".to_string(),
            combined_fitting_report,
        );

        Self::print_quality_table(&report_map);

        Ok(report_map)
    }

    /// Improved fitting method starting with entropy optimization
    ///
    /// Fits NASA-7 coefficients by first optimizing all coefficients using entropy (dS),
    /// then fitting enthalpy coefficient f, and finally validating heat capacity.
    /// Generates comprehensive quality reports with pretty table output.
    ///
    /// # Arguments
    /// * `coeffs_min` - Coefficients for lower temperature range
    /// * `coeffs_max` - Coefficients for upper temperature range
    /// * `T_min` - Lower bound of fitting temperature range
    /// * `T_max` - Upper bound of fitting temperature range
    ///
    /// # Returns
    /// * `Ok(HashMap<String, FittingReport>)` - Quality reports for each property
    /// * `Err(NASAError)` - If temperature ranges are not adjacent or fitting fails
    ///
    /// # Panics
    /// * If fitting quality is below acceptable thresholds (R² < 0.95)
    pub fn fitting_adjacent2(
        &mut self,
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NASAError> {
        let mut report_map: HashMap<String, FittingReport> = HashMap::new();
        let (_, T_center) = (coeffs_min.T.0, coeffs_min.T.1);
        let (T_center_check, _) = (coeffs_max.T.0, coeffs_max.T.1);
        if T_center != T_center_check {
            return Err(NASAError::InvalidTemperatureRange);
        }
        /////////////////////////////fitting ds//////////////////////////////////
        let a = coeffs_min.coeff.0;
        let b = coeffs_min.coeff.1;
        let c = coeffs_min.coeff.2;
        let d = coeffs_min.coeff.3;
        let e = coeffs_min.coeff.4;
        let g = coeffs_min.coeff.6;
        let s_func1 = ds_sym(a, b, c, d, e, g);
        let a = coeffs_max.coeff.0;
        let b = coeffs_max.coeff.1;
        let c = coeffs_max.coeff.2;
        let d = coeffs_max.coeff.3;
        let e = coeffs_max.coeff.4;
        let g = coeffs_max.coeff.6;
        let s_func2 = ds_sym(a, b, c, d, e, g);
        let ds_func_to_fit = Expr::parse_expression(
            "1.987 * (a* ln( t ) + b * t + c * t^2 / 2 + d * t^3 / 3.0 + e * t^4 / 4.0 + g)",
        );

        let mut sew = SewTwoFunctions::new(
            s_func1.clone(),
            s_func2.clone(),
            T_min,
            T_center,
            T_max,
            self.fit_point_number,
        );

        sew.create_fitting_data();
        let initial_guess = vec![
            (coeffs_min.coeff.0 + coeffs_max.coeff.0) / 2.0,
            (coeffs_min.coeff.1 + coeffs_max.coeff.1) / 2.0,
            (coeffs_min.coeff.2 + coeffs_max.coeff.2) / 2.0,
            (coeffs_min.coeff.3 + coeffs_max.coeff.3) / 2.0,
            (coeffs_min.coeff.4 + coeffs_max.coeff.4) / 2.0,
            (coeffs_min.coeff.6 + coeffs_max.coeff.6) / 2.0,
        ];
        sew.fit(
            ds_func_to_fit.clone(),
            Some(vec![
                "a".to_string(),
                "b".to_string(),
                "c".to_string(),
                "d".to_string(),
                "e".to_string(),
                "g".to_string(),
            ]),
            "t".to_string(),
            initial_guess,
            None,
            None,
            None,
            None,
            None,
        );

        let r_squared = sew.get_r_ssquared();
        println!("r_squared for s: {}", r_squared.unwrap());
        assert!(1.0 - r_squared.unwrap() < 5.0e-2);
        let map_of_solutions_6_coeffs = sew.get_map_of_solutions().unwrap();
        let a = *map_of_solutions_6_coeffs.get("a").unwrap();
        let b = *map_of_solutions_6_coeffs.get("b").unwrap();
        let c = *map_of_solutions_6_coeffs.get("c").unwrap();
        let d = *map_of_solutions_6_coeffs.get("d").unwrap();
        let e = *map_of_solutions_6_coeffs.get("e").unwrap();
        let g = *map_of_solutions_6_coeffs.get("g").unwrap();
        let s_fitting_report = self.direct_compare(
            s_func1,
            s_func2,
            T_center,
            ds_func_to_fit,
            map_of_solutions_6_coeffs.clone(),
        );
        report_map.insert("S fitting report".to_string(), s_fitting_report);
        /////////////////// fitting h ///////////////////////
        let f = coeffs_min.coeff.5;
        let h_func1 = dh_sym(a, b, c, d, e, f);
        let f = coeffs_max.coeff.5;
        let h_func2 = dh_sym(a, b, c, d, e, f);
        let dh_func_to_fit = Expr::parse_expression(
            "1.987 *t* (a + b * t/2 + (c/3) * t^2 + (d/4) * t^3 + (e/5) * t^4 + f/t)",
        );
        let dh_func_to_fit = dh_func_to_fit
            .clone()
            .set_variable_from_map(&map_of_solutions_6_coeffs);
        let mut sew = SewTwoFunctions::new(
            h_func1.clone(),
            h_func2.clone(),
            T_min,
            T_center,
            T_max,
            self.fit_point_number,
        );
        sew.create_fitting_data();
        sew.fit(
            dh_func_to_fit.clone(),
            Some(vec!["f".to_string()]),
            "t".to_string(),
            vec![1.0],
            None,
            None,
            None,
            None,
            None,
        );
        let r_squared = sew.get_r_ssquared();
        println!("r_squared for h: {}", r_squared.unwrap());
        assert!(1.0 - r_squared.unwrap() < 5.0e-2);
        let f_map = sew.get_map_of_solutions().unwrap();
        let mut map_of_7_coeffs = map_of_solutions_6_coeffs.clone();
        map_of_7_coeffs.extend(f_map.clone());
        let f = *f_map.get("f").unwrap();
        self.coeffs = (a, b, c, d, e, f, g);

        let h_fitting_report = self.direct_compare(
            h_func1,
            h_func2,
            T_center,
            dh_func_to_fit,
            map_of_7_coeffs.clone(),
        );
        report_map.insert("h fitting report".to_string(), h_fitting_report);
        ////////////////////////////ensure that fitting is correct///////////////
        let Cp_range1 = Cp_sym(
            coeffs_min.coeff.0,
            coeffs_min.coeff.1,
            coeffs_min.coeff.2,
            coeffs_min.coeff.3,
            coeffs_min.coeff.4,
        );
        let Cp_range2 = Cp_sym(
            coeffs_max.coeff.0,
            coeffs_max.coeff.1,
            coeffs_max.coeff.2,
            coeffs_max.coeff.3,
            coeffs_max.coeff.4,
        );
        let Cp_func_to_fit =
            Expr::parse_expression("1.987* (a + b * t + c * t^2 + d * t^3 + e * t^4)");
        let c_fitting_report = self.direct_compare(
            Cp_range1,
            Cp_range2,
            T_center,
            Cp_func_to_fit,
            map_of_solutions_6_coeffs,
        );
        report_map.insert("c fitting report".to_string(), c_fitting_report);

        // Print fitting quality report table
        Self::print_quality_table(&report_map);

        Ok(report_map)
    }
    /// Prints a formatted table of fitting quality metrics
    ///
    /// Displays L1 norm, L2 norm, max norm, number of significant deviation points,
    /// and temperature range where significant deviations occur for each fitted property.
    ///
    /// # Arguments
    /// * `report_map` - HashMap containing fitting reports for different properties
    pub fn print_quality_table(report_map: &HashMap<String, FittingReport>) {
        // Print fitting quality report table
        let mut table = Table::new();
        table.add_row(Row::new(vec![
            Cell::new("Property"),
            Cell::new("L2 Norm"),
            Cell::new("Max Norm"),
            Cell::new("Points"),
            Cell::new("T Range (K)"),
        ]));

        for (property, report) in report_map {
            let t_range_str = match report.T_range {
                Some((min, max)) => format!("{:.1}-{:.1}", min, max),
                None => "N/A".to_string(),
            };

            table.add_row(Row::new(vec![
                Cell::new(property),
                Cell::new(&format!("{:.2e}", report.l2_norm)),
                Cell::new(&format!("{:.2e}", report.max_norm)),
                Cell::new(&report.number_of_significant_points.to_string()),
                Cell::new(&t_range_str),
            ]));
        }

        println!("\nFitting Quality Report:");
        table.printstd();
    }
    /// Main interface for temperature interval coefficient fitting with fallback
    ///
    /// Determines the appropriate fitting strategy based on temperature interval:
    /// - Same interval: Uses existing coefficients
    /// - Adjacent intervals: Attempts fitting_adjacent2, falls back to fitting_adjacent on panic
    /// - Non-adjacent intervals: Returns error
    ///
    /// The fallback mechanism catches panics from assertion failures in fitting_adjacent2
    /// and automatically retries with the more robust fitting_adjacent method.
    ///
    /// # Returns
    /// * `Ok(())` - Coefficients successfully fitted or selected
    /// * `Err(NASAError::InvalidTemperatureRange)` - If no T_interval set or intervals invalid
    ///
    /// # Requires
    /// * `T_interval` must be set via `set_T_interval()` before calling
    /// * Coefficient map must be populated via `parse_coefficients()`
    pub fn fitting_coeffs_for_T_interval(&mut self) -> Result<(), NASAError> {
        if let Some((T_min, T_max)) = self.T_interval {
            let (i_t_min, coeffs_min) = self.interval_for_this_T(T_min)?;
            let (i_t_max, coeffs_max) = self.interval_for_this_T(T_max)?;
            if i_t_min == i_t_max {
                // Both boundaries are in the same interval, no fitting needed
                self.coeffs = coeffs_min.coeff.clone();
            } else if (i_t_max as isize - i_t_min as isize).abs() == 1 {
                // Boundaries are in adjacent intervals, perform fitting
                let result = panic::catch_unwind(panic::AssertUnwindSafe(|| {
                    self.fitting_adjacent2(coeffs_min.clone(), coeffs_max.clone(), T_min, T_max)
                }));

                match result {
                    Ok(Ok(_)) => {} // fitting_adjacent2 succeeded
                    _ => {
                        // fitting_adjacent2 panicked or failed, use fitting_adjacent as fallback
                        println!("fitting_adjacent2 failed, falling back to fitting_adjacent");
                        self.fitting_adjacent3(coeffs_min, coeffs_max, T_min, T_max)?;
                    }
                }
            } else {
                return Err(NASAError::InvalidTemperatureRange);
            }

            Ok(())
        } else {
            Err(NASAError::InvalidTemperatureRange)
        }
    }
}
#[derive(Debug, Clone)]
pub struct FittingReport {
    l2_norm: f64,
    max_norm: f64,
    number_of_significant_points: usize,
    T_range: Option<(f64, f64)>,
}
impl fmt::Debug for NASAdata {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("NASAdata")
            .field("input", &self.input)
            .field("unit", &self.unit)
            .field("unit_multiplier", &self.unit_multiplier)
            .field("coeffs", &self.coeffs)
            .field("Cp", &self.Cp)
            .field("dh", &self.dh)
            .field("ds", &self.ds)
            .field("Cp_sym", &self.Cp_sym)
            .field("dh_sym", &self.dh_sym)
            .field("ds_sym", &self.ds_sym)
            .finish_non_exhaustive() //  The finish_non_exhaustive() method is used to indicate that not all fields are being displayed.
    }
}

impl Clone for NASAdata {
    fn clone(&self) -> Self {
        let (a, b, c, d, e, f, g) = self.coeffs;
        let unit = self.unit_multiplier;
        NASAdata {
            input: self.input.clone(),
            unit: self.unit.clone(),
            unit_multiplier: self.unit_multiplier,
            coeffs: self.coeffs,
            coeffs_map: self.coeffs_map.clone(),
            Cp: self.Cp,
            dh: self.dh,
            ds: self.ds,
            C_fun: Box::new(move |t| {
                unit * R * (a + b * t + c * t.powi(2) + d * t.powi(3) + e * t.powi(4))
            }),
            dh_fun: Box::new(move |t| {
                unit * R
                    * t
                    * (a + b * t / 2.0
                        + c * t.powi(2) / 3.0
                        + d * t.powi(3) / 4.0
                        + e * t.powi(4) / 5.0
                        + f / t)
            }),
            ds_fun: Box::new(move |t| {
                unit * R
                    * (a * t.ln()
                        + b * t
                        + c * t.powi(2) / 2.0
                        + d * t.powi(3) / 3.0
                        + e * t.powi(4) / 4.0
                        + g)
            }),
            Cp_sym: self.Cp_sym.clone(),
            dh_sym: self.dh_sym.clone(),
            ds_sym: self.ds_sym.clone(),
            T_interval: self.T_interval.clone(),
            fit_point_number: self.fit_point_number,
            norm_threshold: self.norm_threshold.clone(),
        }
    }
}
use super::thermo_api::{EnergyUnit, ThermoCalculator, ThermoError, energy_dimension};
impl ThermoCalculator for NASAdata {
    fn newinstance(&mut self) -> Result<(), ThermoError> {
        *self = NASAdata::new();
        Ok(())
    }
    fn extract_model_coefficients(&mut self, t: f64) -> Result<(), ThermoError> {
        self.extract_coefficients(t)?;
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
        self.calculate_Cp_dH_dS(temperature);
        Ok(())
    }
    fn create_closures_Cp_dH_dS(&mut self) -> Result<(), ThermoError> {
        self.create_closures_Cp_dH_dS();
        Ok(())
    }
    fn create_sym_Cp_dH_dS(&mut self) -> Result<(), ThermoError> {
        self.create_sym_Cp_dH_dS();
        Ok(())
    }
    fn Taylor_series_cp_dh_ds(
        &mut self,
        temperature: f64,
        order: usize,
    ) -> Result<(Expr, Expr, Expr), ThermoError> {
        let (Cp, dh, ds) = self.Taylor_series_Cp_dH_dS(temperature, order);
        Ok((Cp, dh, ds))
    }
    fn pretty_print_data(&self) -> Result<(), ThermoError> {
        self.pretty_print();
        Ok(())
    }
    fn renew_base(
        &mut self,
        _sub_name: String,
        _search_type: SearchType,
        _phase: Phase,
    ) -> Result<(), ThermoError> {
        Ok(())
    }
    fn get_coefficients(&self) -> Result<Vec<f64>, ThermoError> {
        let (a, b, c, d, e, f, g) = self.coeffs;
        Ok(vec![a, b, c, d, e, f, g])
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
        Ok(self.input.composition.clone())
    }
}
