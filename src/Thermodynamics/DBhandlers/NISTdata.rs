use crate::Thermodynamics::DBhandlers::NIST_parser::{NistInput, NistParser};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::error::Error;
use std::fmt;
use std::str::FromStr;

#[derive(Debug)]
pub enum NISTError {
    NoCoefficientsFound { temperature: f64, range: String },
    InvalidTemperatureRange,
    DeserializationError(String),
    UnsupportedUnit(String),
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
            NISTError::DeserializationError(msg) => {
                write!(f, "Failed to deserialize NASA data: {}", msg)
            }
            NISTError::UnsupportedUnit(unit) => {
                write!(
                    f,
                    "Unsupported unit: {}. Only 'J' and 'cal' are supported",
                    unit
                )
            }
        }
    }
}

impl Error for NISTError {}

pub struct NISTdata {
    /// data parsed from library
    pub input: NistInput,
    /// Joles or Cal
    pub unit: Option<String>,
    ///

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

impl NISTdata {
    pub fn new() -> Self {
        let input = NistInput::new();
        Self {
            input: input,
            unit: None,

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
    pub fn set_unit(&mut self, unit: &str) -> Result<(), NISTError> {
        self.input.set_unit(unit);
        Ok(())
    }
    /// takes serde Value and parse it into structure
    pub fn from_serde(&mut self, serde: Value) -> Result<(), NISTError> {
        self.input = serde_json::from_value(serde)
            .map_err(|e| NISTError::DeserializationError(e.to_string()))?;
        Ok(())
    }
    pub fn create_closures_Cp_dH_dS(&mut self, T: f64) {
        let (C_fun, dh_fun, ds_fun) = self
            .input
            .create_closure_cp_dh_ds(T)
            .expect("Error calculating cp, dh, ds");
        (self.C_fun, self.dh_fun, self.ds_fun) = (C_fun, dh_fun, ds_fun);
    }
    pub fn create_sym_Cp_dH_dS(&mut self, T: f64) {
        (self.Cp_sym, self.dh_sym, self.ds_sym) = self
            .input
            .create_sym_cp_dh_ds(T)
            .expect("Error calculating cp, dh, ds");
    }
    pub fn calculate_Cp_dH_dS(&mut self, T: f64) {
        (self.Cp, self.dh, self.ds) = self
            .input
            .caclc_cp_dh_ds(T)
            .expect("Error calculating cp, dh, ds");
    }
    pub fn pretty_print(&self) -> Result<(), NISTError> {
        self.input.pretty_print();
        Ok(())
    }
}
