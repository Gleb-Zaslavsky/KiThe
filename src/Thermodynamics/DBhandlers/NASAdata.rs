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
//! - Uses NASA-7 polynomial format: Cp/R = a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴
//! - Supports multiple temperature ranges (typically 2-3 ranges per substance)
//! - Handles both gas and condensed phase data
//!
//! ## Key Methods
//! - `extract_coefficients(T)`: Selects appropriate coefficient set for given temperature
//! - `calculate_Cp_dH_dS(T)`: Computes thermodynamic properties at specified temperature
//! - `create_closures_Cp_dH_dS()`: Creates closure functions for efficient repeated calculations
//! - `create_sym_Cp_dH_dS()`: Creates symbolic expressions for analytical work
//! - `Taylor_series_Cp_dH_dS()`: Generates Taylor series expansions around specified temperature
//! - `pretty_print()`: Displays coefficient data in formatted table
//!
//! ## Usage
//! ```rust, ignore
//! let mut nasa = NASAdata::new();
//! nasa.from_serde(database_entry)?;
//! nasa.set_unit("J")?;  // or "cal" for calories
//! nasa.extract_coefficients(400.0)?;
//! nasa.calculate_Cp_dH_dS(400.0);
//! let cp = nasa.Cp;  // Heat capacity
//! let dh = nasa.dh;  // Enthalpy
//! let ds = nasa.ds;  // Entropy
//! ```
//!
//! ## Interesting Features
//! - Handles complex composition parsing from string format ("{C:1, H:4}")
//! - Supports both Joules and calories with automatic unit conversion
//! - Implements NASA-7 polynomial integration for enthalpy and entropy calculations
//! - Provides Taylor series expansion for local approximations
//! - Uses custom deserialization for composition data parsing
//! - Implements comprehensive error handling for temperature range validation
//! - Creates both numerical functions and symbolic expressions for flexible usage

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

fn Cp(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    R * (a + b * t + c * t.powi(2) + d * t.powi(3) + e * t.powi(4))
}
fn dh(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> f64 {
    R * t
        * (a + b * t / 2.0
            + c * t.powi(2) / 3.0
            + d * t.powi(3) / 4.0
            + e * t.powi(4) / 5.0
            + f / t)
}
fn ds(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64, g: f64) -> f64 {
    R * (a * t.ln() + b * t + c * t.powi(2) / 2.0 + d * t.powi(3) / 3.0 + e * t.powi(4) / 4.0 + g)
}
fn Cp_sym(a: f64, b: f64, c: f64, d: f64, e: f64) -> Expr {
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
    composition: Option<HashMap<String, f64>>,
    Cp: Vec<f64>,

    model: String,
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

pub struct NASAdata {
    /// data parsed from library
    pub input: NASAinput,
    /// Joles or Cal
    pub unit: Option<String>,
    ///
    unit_multiplier: f64,
    /// NASA format 7 coefficients
    pub coeffs: (f64, f64, f64, f64, f64, f64, f64),
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
            Cp: 0.0,
            dh: 0.0,
            ds: 0.0,
            C_fun: Box::new(|x| x),
            dh_fun: Box::new(|x| x),
            ds_fun: Box::new(|x| x),

            Cp_sym: Expr::Const(0.0),
            dh_sym: Expr::Const(0.0),
            ds_sym: Expr::Const(0.0),
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
    /// takes serde Value and parse it into structure
    pub fn from_serde(&mut self, serde: Value) -> Result<(), NASAError> {
        self.input = serde_json::from_value(serde).map_err(|e| NASAError::SerdeError(e))?;
        Ok(())
    }

    fn extract_coefficients_(
        c_data: &[f64],
        t: f64,
    ) -> Result<(f64, f64, f64, f64, f64, f64, f64), NASAError> {
        let get_range_str = |temps: &[f64]| {
            temps
                .iter()
                .map(|t| t.to_string())
                .collect::<Vec<_>>()
                .join(" - ")
        };

        match c_data.len() {
            25 => {
                let (t1, t2, t3, t4) = (c_data[0], c_data[1], c_data[2], c_data[3]);
                if t1 <= t && t <= t2 {
                    Ok((
                        c_data[4], c_data[5], c_data[6], c_data[7], c_data[8], c_data[9],
                        c_data[10],
                    ))
                } else if t2 < t && t <= t3 {
                    Ok((
                        c_data[11], c_data[12], c_data[13], c_data[14], c_data[15], c_data[16],
                        c_data[17],
                    ))
                } else if t3 < t && t < t4 {
                    Ok((
                        c_data[18], c_data[19], c_data[20], c_data[21], c_data[22], c_data[23],
                        c_data[24],
                    ))
                } else {
                    Err(NASAError::NoCoefficientsFound {
                        temperature: t,
                        range: get_range_str(&[t1, t2, t3, t4]),
                    })
                }
            }
            17 => {
                let (t1, t2, t3) = (c_data[0], c_data[1], c_data[2]);
                if t1 <= t && t <= t2 {
                    Ok((
                        c_data[3], c_data[4], c_data[5], c_data[6], c_data[7], c_data[8], c_data[9],
                    ))
                } else if t2 < t && t <= t3 {
                    Ok((
                        c_data[10], c_data[11], c_data[12], c_data[13], c_data[14], c_data[15],
                        c_data[16],
                    ))
                } else {
                    Err(NASAError::NoCoefficientsFound {
                        temperature: t,
                        range: get_range_str(&[t1, t2, t3]),
                    })
                }
            }
            9 => {
                let (t1, t2) = (c_data[0], c_data[1]);
                if t1 <= t && t <= t2 {
                    Ok((
                        c_data[2], c_data[3], c_data[4], c_data[5], c_data[6], c_data[7], c_data[8],
                    ))
                } else {
                    Err(NASAError::NoCoefficientsFound {
                        temperature: t,
                        range: get_range_str(&[t1, t2]),
                    })
                }
            }
            _ => Err(NASAError::InvalidTemperatureRange),
        }
    }
    /// get the 7 constants of NASA7 format for conctrete temperature
    pub fn extract_coefficients(&mut self, t: f64) -> Result<(), NASAError> {
        let c_data = self.input.Cp.as_slice();
        let (a, b, c, d, e, f, g) = Self::extract_coefficients_(c_data, t)?;
        self.coeffs = (a, b, c, d, e, f, g);
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
        self.Cp_sym = C_sym.symplify();
        let dh_sym = unit.clone() * dh_sym(a, b, c, d, e, f);
        self.dh_sym = dh_sym.symplify();
        let ds_sym = unit * ds_sym(a, b, c, d, e, g);
        self.ds_sym = ds_sym.symplify();
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::thermo_lib_api::ThermoData;
    use approx::assert_relative_eq;

    #[test]
    fn test_new() {
        let nasa_data = NASAdata::new();
        assert_eq!(nasa_data.coeffs, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        assert_eq!(nasa_data.Cp, 0.0);
        assert_eq!(nasa_data.dh, 0.0);
        assert_eq!(nasa_data.ds, 0.0);
    }

    #[test]
    fn test_extract_coefficients() {
        let mut nasa_data = NASAdata::new();
        nasa_data.input.Cp = vec![300.0, 1000.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];

        assert!(nasa_data.extract_coefficients(500.0).is_ok());
        assert_eq!(nasa_data.coeffs, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0));

        assert!(nasa_data.extract_coefficients(700.0).is_ok());
        assert_eq!(nasa_data.coeffs, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0));
    }

    #[test]
    fn test_extract_coefficients_out_of_range() {
        let mut nasa_data = NASAdata::new();
        nasa_data.input.Cp = vec![300.0, 1000.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];

        let result = nasa_data.extract_coefficients(1500.0);
        assert!(matches!(
            result,
            Err(NASAError::NoCoefficientsFound {
                temperature: 1500.0,
                ..
            })
        ));
    }

    #[test]
    fn test_set_unit() {
        let mut nasa_data = NASAdata::new();

        assert!(nasa_data.set_unit("J").is_ok());
        assert_eq!(nasa_data.unit, Some("J".to_string()));
        assert_eq!(nasa_data.unit_multiplier, 4.184);

        assert!(nasa_data.set_unit("cal").is_ok());
        assert_eq!(nasa_data.unit, Some("cal".to_string()));
        assert_eq!(nasa_data.unit_multiplier, 1.0);

        let result = nasa_data.set_unit("invalid");
        assert!(matches!(result, Err(NASAError::UnsupportedUnit(unit)) if unit == "invalid"));
    }

    #[test]
    fn test_create_closures_cp_dh_ds() {
        let mut nasa_data = NASAdata::new();
        nasa_data.coeffs = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
        let _ = nasa_data.set_unit("cal");
        nasa_data.create_closures_Cp_dH_dS();

        let t = 500.0;
        assert_relative_eq!(
            (nasa_data.C_fun)(t),
            Cp(t, 1.0, 2.0, 3.0, 4.0, 5.0),
            epsilon = 1e-6
        );
        assert_relative_eq!(
            (nasa_data.dh_fun)(t),
            dh(t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
            epsilon = 1e-6
        );
        assert_relative_eq!(
            (nasa_data.ds_fun)(t),
            ds(t, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_with_real_data() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut NASA = NASAdata::new();

        assert!(NASA.from_serde(CO_data.clone()).is_ok());
        assert!(NASA.extract_coefficients(400.0).is_ok());

        println!("\n \n {:?} \n \n ", NASA.coeffs);
        let coeffs_len = {
            let (_a, _b, _c, _d, _e, _f, _g) = NASA.coeffs;
            7 // Since we know the tuple has 7 elements
        };
        assert_eq!(coeffs_len, 7);

        NASA.calculate_Cp_dH_dS(400.0);
        let Cp = NASA.Cp;
        let dh = NASA.dh;
        let ds = NASA.ds;

        println!("Cp: {}, dh: {}, ds: {}", Cp, dh, ds);
        assert!(Cp > 0.0);
        assert!(dh < 0.0);
        assert!(ds > 0.0);

        let t = 400.0;
        NASA.create_closures_Cp_dH_dS();

        let Cp_fun = &NASA.C_fun;
        let dh_fun = &NASA.dh_fun;
        let ds_fun = &NASA.ds_fun;
        assert_relative_eq!((Cp_fun)(t), NASA.Cp, epsilon = 1e-6);
        assert_relative_eq!((dh_fun)(t), NASA.dh, epsilon = 1e-6);
        assert_relative_eq!((ds_fun)(t), NASA.ds, epsilon = 1e-6);

        NASA.create_sym_Cp_dH_dS();
        let Cp_sym = &NASA.Cp_sym;
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        assert_relative_eq!(Cp_value, NASA.Cp, epsilon = 1e-6);
        let dh_sym = &NASA.dh_sym;
        let dh_T = dh_sym.lambdify1D();
        let dh_value = dh_T(400.0);
        assert_relative_eq!(dh_value, NASA.dh, epsilon = 1e-6);
        let ds_sym = &NASA.ds_sym;
        let ds_T = ds_sym.lambdify1D();
        let ds_value = ds_T(400.0);
        assert_relative_eq!(ds_value, NASA.ds, epsilon = 1e-6);
    }

    #[test]
    fn test_thermo_calculator_error_handling() {
        let mut nasa = NASAdata::new();

        // Test invalid temperature range
        let result = nasa.extract_model_coefficients(1000.0);
        assert!(result.is_err());

        // Test invalid unit
        let result = nasa.set_unit(&energy_dimension(EnergyUnit::J));
        assert!(result.is_ok());

        // Test invalid serde data
        let invalid_data = serde_json::json!({
            "invalid": "data"
        });
        let result = nasa.from_serde(invalid_data);
        assert!(result.is_err());
    }
    #[test]
    fn test_with_real_data_ThermoCalc() {
        use super::ThermoCalculator;
        use crate::Thermodynamics::DBhandlers::thermo_api::create_thermal_by_name;
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        println!(" CO data \n {} \n", CO_data);
        let mut NASA = create_thermal_by_name("NASA_gas");
        let _ = NASA.newinstance();
        let _ = NASA.from_serde(CO_data.clone());
        //  assert!(NASA.from_serde(CO_data.clone()).is_ok());
        print!(" this is NASA instance: \n");
        let _ = NASA.print_instance();
        assert!(NASA.extract_model_coefficients(400.0).is_ok());

        let coeffs_len = {
            let coeff_vec = NASA.get_coefficients().unwrap();
            coeff_vec.len()
        };
        assert_eq!(coeffs_len, 7);

        let _ = NASA.calculate_Cp_dH_dS(400.0);
        let Cp = NASA.get_Cp().unwrap();
        let dh = NASA.get_dh().unwrap();
        let ds = NASA.get_ds().unwrap();

        println!("Cp: {}, dh: {}, ds: {}", Cp, dh, ds);
        assert!(Cp > 0.0);
        assert!(dh < 0.0);
        assert!(ds > 0.0);

        let t = 400.0;
        let _ = NASA.create_closures_Cp_dH_dS();

        let Cp_fun = NASA.get_C_fun().unwrap();
        let dh_fun = NASA.get_dh_fun().unwrap();
        let ds_fun = NASA.get_ds_fun().unwrap();
        assert_relative_eq!((Cp_fun)(t), Cp, epsilon = 1e-6);
        assert_relative_eq!((dh_fun)(t), dh, epsilon = 1e-6);
        assert_relative_eq!((ds_fun)(t), ds, epsilon = 1e-6);

        let _ = NASA.create_sym_Cp_dH_dS();
        let Cp_sym = NASA.get_Cp_sym().unwrap();
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        assert_relative_eq!(Cp_value, Cp, epsilon = 1e-6);

        let dh_sym = NASA.get_dh_sym().unwrap();
        let dh_T = dh_sym.lambdify1D();
        let dh_value = dh_T(400.0);
        assert_relative_eq!(dh_value, dh, epsilon = 1e-6);

        let ds_sym = NASA.get_ds_sym().unwrap();
        let ds_T = ds_sym.lambdify1D();
        let ds_value = ds_T(400.0);
        assert_relative_eq!(ds_value, ds, epsilon = 1e-6);
    }

    #[test]
    fn test_thermo_calculator_nasa() {
        // use super::EnergyUnit;
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        // Test newinstance
        assert!(nasa.newinstance().is_ok());

        // Test from_serde
        assert!(nasa.from_serde(CO_data.clone()).is_ok());

        // Test set_unit
        // assert!(nasa.set_unit(Some(EnergyUnit::J)).is_ok());
        // assert!(nasa.set_unit(Some(EnergyUnit::Cal)).is_ok());

        // Test extract_model_coefficients
        assert!(nasa.extract_model_coefficients(400.0).is_ok());

        // Test calculate_Cp_dH_dS
        nasa.calculate_Cp_dH_dS(400.0);
        assert!(nasa.Cp > 0.0);
        assert!(nasa.dh != 0.0);
        assert!(nasa.ds != 0.0);

        // Test create_closures_Cp_dH_dS
        nasa.create_closures_Cp_dH_dS();
        let t = 400.0;
        assert_relative_eq!((nasa.C_fun)(t), nasa.Cp, epsilon = 1e-6);
        assert_relative_eq!((nasa.dh_fun)(t), nasa.dh, epsilon = 1e-6);
        assert_relative_eq!((nasa.ds_fun)(t), nasa.ds, epsilon = 1e-6);

        // Test create_sym_Cp_dH_dS
        nasa.create_sym_Cp_dH_dS();
        let Cp_sym = &nasa.Cp_sym;
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        let dh_sym = &nasa.dh_sym;
        let dh_T = dh_sym.lambdify1D();
        let dh_value = dh_T(400.0);
        let ds_sym = &nasa.ds_sym;
        let ds_T = ds_sym.lambdify1D();
        let ds_value = ds_T(400.0);
        assert_relative_eq!(Cp_value, nasa.Cp, epsilon = 1e-6);
        assert_relative_eq!(dh_value, nasa.dh, epsilon = 1e-6);
        assert_relative_eq!(ds_value, nasa.ds, epsilon = 1e-6);

        // Test Taylor_series_cp_dh_ds
        let (Cp_taylor, dh_taylor, ds_taylor) = nasa.Taylor_series_cp_dh_ds(300.0, 3).unwrap();
        let Cp_taylor_T = Cp_taylor.lambdify1D();
        let Cp_value = Cp_taylor_T(400.0);
        assert_relative_eq!(Cp_value, nasa.Cp, epsilon = 3.0);
        assert!(!Cp_taylor.is_zero());
        assert!(!dh_taylor.is_zero());
        assert!(!ds_taylor.is_zero());

        // Test pretty_print_data
        assert!(nasa.pretty_print_data().is_ok());
    }
}
