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

use RustedSciThe::symbolic::symbolic_engine::Expr;

use prettytable::{Cell, Row, Table};
use serde::de::Deserializer;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;

use std::error::Error;
use std::fmt;

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
    SymbolicError(String),
    FittingError(String),
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
            NASAError::SymbolicError(msg) => {
                write!(f, "Symbolic computation error: {}", msg)
            }
            NASAError::FittingError(msg) => {
                write!(f, "Fitting error: {}", msg)
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
pub fn dh_sym(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> Expr {
    let t = Expr::Var("T".to_owned());
    let (a, b, c, d, e, f) = (
        Expr::Const(a),
        Expr::Const(b),
        Expr::Const(c),
        Expr::Const(d),
        Expr::Const(e),
        Expr::Const(f),
    );
    Rsym * (a * t.clone()
        + b * t.clone().pow(Expr::Const(2.0)) / Expr::Const(2.0)
        + c * t.clone().pow(Expr::Const(3.0)) / Expr::Const(3.0)
        + d * t.clone().pow(Expr::Const(4.0)) / Expr::Const(4.0)
        + e * t.clone().pow(Expr::Const(5.0)) / Expr::Const(5.0)
        + f)
}
pub fn ds_sym(a: f64, b: f64, c: f64, d: f64, e: f64, g: f64) -> Expr {
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
    // check what coefficients inerval belongs current temperature
    pub current_T_range: Option<(f64, f64)>,
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
            current_T_range: None,
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

    pub fn parse_coefficients(&mut self) -> Result<(), NASAError> {
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
        Ok(())
    }

    fn find_coefficients_for_temperature(
        &mut self,
        t: f64,
    ) -> Result<(f64, f64, f64, f64, f64, f64, f64), NASAError> {
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

        let coeffs_names = ["A", "B", "C", "D", "E", "F", "G"];
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
    fn is_this_T_from_current_T_range(&self, t: f64) -> bool {
        if let Some((T_min, T_max)) = self.current_T_range {
            T_min <= t && t <= T_max
        } else {
            false
        }
    }

    pub fn with_T_range<F>(&mut self, T: f64, function: F) -> Result<(), NASAError>
    where
        F: FnOnce(),
    {
        let is_this_T_from_current_T_range = self.is_this_T_from_current_T_range(T);
        if is_this_T_from_current_T_range {
            function();
            Ok(())
        } else {
            self.interval_for_this_T(T)?;
            // if no - change coefficients
            self.extract_coefficients(T)?;
            function();
            Ok(())
        }
    }

    pub fn with_T_range_mut<F, R>(&mut self, T: f64, function: F) -> Result<R, NASAError>
    where
        F: FnOnce() -> R,
    {
        let is_this_T_from_current_T_range = self.is_this_T_from_current_T_range(T);
        if is_this_T_from_current_T_range {
            Ok(function())
        } else {
            self.interval_for_this_T(T)?;
            self.extract_coefficients(T)?;
            Ok(function())
        }
    }

    /// Creates closure functions for Cp, dH, dS with automatic temperature range handling
    ///
    /// Automatically selects appropriate NASA coefficients for the given temperature
    /// and creates closure functions. If temperature is outside current range,
    /// updates coefficients before creating closures.
    ///
    /// # Arguments
    /// * `T` - Temperature in Kelvin
    ///
    /// # Returns
    /// * `Ok(())` - Closures successfully created
    /// * `Err(NASAError)` - If temperature is outside all available ranges
    pub fn create_closures_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), NASAError> {
        let is_this_T_from_current_T_range = self.is_this_T_from_current_T_range(T);
        if !is_this_T_from_current_T_range {
            self.interval_for_this_T(T)?;
            self.extract_coefficients(T)?;
        }
        self.create_closures_Cp_dH_dS();
        Ok(())
    }

    /// Calculates Cp, dH, dS values with automatic temperature range handling
    ///
    /// Automatically selects appropriate NASA coefficients for the given temperature
    /// and calculates thermodynamic properties. If temperature is outside current range,
    /// updates coefficients before calculation.
    ///
    /// # Arguments
    /// * `T` - Temperature in Kelvin
    ///
    /// # Returns
    /// * `Ok(())` - Values successfully calculated and stored in Cp, dh, ds fields
    /// * `Err(NASAError)` - If temperature is outside all available ranges
    pub fn calculate_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), NASAError> {
        let is_this_T_from_current_T_range = self.is_this_T_from_current_T_range(T);
        if !is_this_T_from_current_T_range {
            self.interval_for_this_T(T)?;
            self.extract_coefficients(T)?;
        }
        self.calculate_Cp_dH_dS(T);
        Ok(())
    }

    /// Creates symbolic expressions for Cp, dH, dS with automatic temperature range handling
    ///
    /// Automatically selects appropriate NASA coefficients for the given temperature
    /// and creates symbolic expressions. If temperature is outside current range,
    /// updates coefficients before creating expressions.
    ///
    /// # Arguments
    /// * `T` - Temperature in Kelvin
    ///
    /// # Returns
    /// * `Ok(())` - Symbolic expressions successfully created
    /// * `Err(NASAError)` - If temperature is outside all available ranges
    pub fn create_sym_Cp_dH_dS_with_T_range(&mut self, T: f64) -> Result<(), NASAError> {
        let is_this_T_from_current_T_range = self.is_this_T_from_current_T_range(T);
        if !is_this_T_from_current_T_range {
            self.interval_for_this_T(T)?;
            self.extract_coefficients(T)?;
        }
        self.create_sym_Cp_dH_dS();
        Ok(())
    }
    /// Finds the coefficient interval that contains the given temperature
    ///
    /// # Arguments
    /// * `t` - Temperature in Kelvin
    ///
    /// # Returns
    /// * `Ok((interval_index, coeffs))` - Index and coefficients for the temperature range
    /// * `Err(NASAError::InvalidTemperatureRange)` - If temperature is outside all ranges
    pub fn interval_for_this_T(&mut self, t: f64) -> Result<(usize, Coeffs), NASAError> {
        let intervals = self.coeffs_map.clone();

        for (interval_idx, coeffs) in intervals {
            if coeffs.T.0 <= t && t <= coeffs.T.1 {
                self.current_T_range = Some(coeffs.T);
                return Ok((interval_idx, coeffs.clone()));
            }
        }

        Err(NASAError::InvalidTemperatureRange)
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

                let result =
                    self.fitting_adjacent2(coeffs_min.clone(), coeffs_max.clone(), T_min, T_max);

                match result {
                    Ok(_) => {
                        self.current_T_range = self.T_interval;
                    } // fitting_adjacent2 succeeded
                    _ => {
                        // fitting_adjacent2 panicked or failed, use fitting_adjacent as fallback
                        println!("fitting_adjacent2 failed, falling back to fitting_adjacent");
                        self.fit_cp_dh_ds(coeffs_min.clone(), coeffs_max.clone(), T_min, T_max)?;
                        self.current_T_range = self.T_interval;
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
                    Ok(_) => {
                        self.current_T_range = self.T_interval;
                    } // fitting_non_adjacent succeeded
                    _ => {
                        self.fitting_cp_dh_ds_non_adjacent(coeffs_in_range, T_min, T_max)?;
                        self.current_T_range = self.T_interval;
                    }
                }
            } else {
                // this is unreachable
                return Err(NASAError::InvalidTemperatureRange);
            }
            Ok(())
        } else {
            Err(NASAError::InvalidTemperatureRange)
        }
    }
    /// Calculates temperature-averaged mean values of thermodynamic properties
    ///
    /// Computes the mean values of heat capacity (Cp), enthalpy (dH), and entropy (dS)
    /// over the specified temperature interval using symbolic integration:
    ///
    /// mean_value = (1/(T_max - T_min)) * ∫[T_min to T_max] property(T) dT
    ///
    /// The method uses the symbolic expressions stored in `Cp_sym`, `dh_sym`, and `ds_sym`
    /// to perform definite integration over the temperature range, then divides by the
    /// temperature interval to obtain mean values.
    ///
    /// # Returns
    /// * `Ok(())` - Mean values successfully calculated and stored in Cp, dh, ds fields
    /// * `Err(NASAError::InvalidTemperatureRange)` - If T_interval is not set or invalid
    /// * `Err(NASAError::SymbolicError)` - If symbolic integration fails
    ///
    /// # Requires
    /// * `T_interval` must be set via `set_T_interval()` before calling
    /// * Symbolic expressions must be created via `create_sym_Cp_dH_dS()` before calling
    /// * T_min < T_max for valid temperature range
    ///
    /// # Example
    /// ```rust, ignore
    /// nasa.set_T_interval(400.0, 800.0);
    /// nasa.create_sym_Cp_dH_dS();
    /// nasa.integr_mean()?;
    /// let mean_cp = nasa.Cp;  // Average heat capacity over 400-800K
    /// ```
    pub fn integr_mean(&mut self) -> Result<(), NASAError> {
        let (T_min, T_max) = self.T_interval.ok_or(NASAError::InvalidTemperatureRange)?;

        if T_min >= T_max {
            return Err(NASAError::InvalidTemperatureRange);
        }

        let T_range = T_max - T_min;
        let inv_T_range = 1.0 / T_range;

        // Integrate entropy
        let ds_sym = self.ds_sym.clone();
        let ds_integral = ds_sym
            .definite_integrate("T", T_min, T_max)
            .map_err(|e| NASAError::SymbolicError(format!("Failed to integrate entropy: {}", e)))?;
        let ds_mean = inv_T_range * ds_integral;

        // Integrate enthalpy
        let dh_sym = self.dh_sym.clone().simplify();
        let dh_integral = dh_sym.definite_integrate("T", T_min, T_max).map_err(|e| {
            NASAError::SymbolicError(format!("Failed to integrate enthalpy: {}", e))
        })?;
        let dh_mean = inv_T_range * dh_integral;

        // Integrate heat capacity
        let Cp_sym = self.Cp_sym.clone();
        let Cp_integral = Cp_sym.definite_integrate("T", T_min, T_max).map_err(|e| {
            NASAError::SymbolicError(format!("Failed to integrate heat capacity: {}", e))
        })?;
        let Cp_mean = inv_T_range * Cp_integral;

        // Store mean values
        self.Cp = Cp_mean;
        self.dh = dh_mean;
        self.ds = ds_mean;

        Ok(())
    }
}
//////////////////////////////////////////////////////////////////////////////////////////

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
            current_T_range: self.current_T_range.clone(),
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
    fn fitting_coeffs_for_T_interval(&mut self) -> Result<(), ThermoError> {
        self.fitting_coeffs_for_T_interval()?;
        Ok(())
    }
    fn integr_mean(&mut self) -> Result<(), ThermoError> {
        self.integr_mean()?;
        Ok(())
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
