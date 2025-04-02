use crate::Thermodynamics::DBhandlers::NIST_parser::{Phase, SearchType};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use enum_dispatch::enum_dispatch;
use serde_json::Value;
use std::error::Error;
use std::fmt;
#[derive(Debug)]
pub enum ThermoError {
    NoCoefficientsFound { temperature: f64, range: String },
    InvalidTemperatureRange,

    UnsupportedUnit(String),
    SerdeError(serde_json::Error),
}

impl fmt::Display for ThermoError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ThermoError::NoCoefficientsFound { temperature, range } => {
                write!(
                    f,
                    "No coefficients found for temperature {} K. Valid range: {}",
                    temperature, range
                )
            }
            ThermoError::InvalidTemperatureRange => {
                write!(f, "Invalid temperature range in coefficient data")
            }
            ThermoError::SerdeError(msg) => {
                write!(f, "Failed to deserialize NASA data: {}", msg)
            }
            ThermoError::UnsupportedUnit(unit) => {
                write!(
                    f,
                    "Unsupported unit: {}. Only 'J' and 'cal' are supported",
                    unit
                )
            }
        }
    }
}

impl Error for ThermoError {}

impl From<super::NISTdata::NISTError> for ThermoError {
    fn from(err: super::NISTdata::NISTError) -> Self {
        match err {
            super::NISTdata::NISTError::NoCoefficientsFound { temperature, range } => {
                ThermoError::NoCoefficientsFound { temperature, range }
            }
            super::NISTdata::NISTError::InvalidTemperatureRange => {
                ThermoError::InvalidTemperatureRange
            }
            super::NISTdata::NISTError::SerdeError(msg) => ThermoError::SerdeError(msg),
            super::NISTdata::NISTError::UnsupportedUnit(unit) => ThermoError::UnsupportedUnit(unit),
            //    super::NISTdata::NISTError::SerdeError(err) => ThermoError::SerdeError(err),
        }
    }
}

impl From<super::NASAdata::NASAError> for ThermoError {
    fn from(err: super::NASAdata::NASAError) -> Self {
        match err {
            super::NASAdata::NASAError::NoCoefficientsFound { temperature, range } => {
                ThermoError::NoCoefficientsFound { temperature, range }
            }
            super::NASAdata::NASAError::InvalidTemperatureRange => {
                ThermoError::InvalidTemperatureRange
            }
            super::NASAdata::NASAError::SerdeError(msg) => ThermoError::SerdeError(msg),
            super::NASAdata::NASAError::UnsupportedUnit(unit) => ThermoError::UnsupportedUnit(unit),
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum EnergyUnit {
    J,
    Cal,
}
#[enum_dispatch]
pub trait ThermoCalculator {
    fn newinstance(&mut self) -> Result<(), ThermoError>;
    fn set_unit(&mut self, unit: Option<EnergyUnit>) -> Result<(), ThermoError>;
    fn from_serde(&mut self, serde: Value) -> Result<(), ThermoError>;
    fn extract_model_coefficients(&mut self, t: f64) -> Result<(), ThermoError>;
    fn calculate_Cp_dH_dS(&mut self, t: f64) -> Result<(), ThermoError>;
    fn create_closures_Cp_dH_dS(&mut self) -> Result<(), ThermoError>;
    fn create_sym_Cp_dH_dS(&mut self) -> Result<(), ThermoError>;
    fn Taylor_series_cp_dh_ds(
        &mut self,
        temperature: f64,
        order: usize,
    ) -> Result<(Expr, Expr, Expr), ThermoError>;
    fn pretty_print_data(&self) -> Result<(), ThermoError>;
    fn renew_base(
        &mut self,
        sub_name: String,
        search_type: SearchType,
        phase: Phase,
    ) -> Result<(), ThermoError>;
    fn get_coefficients(&self) -> Result<Vec<f64>, ThermoError>;

    // Add getter methods for thermodynamic properties
    fn print_instance(&self) -> Result<(), ThermoError>;
    fn get_Cp(&self) -> Result<f64, ThermoError>;
    fn get_dh(&self) -> Result<f64, ThermoError>;
    fn get_ds(&self) -> Result<f64, ThermoError>;
    fn get_C_fun(&self) -> Result<Box<dyn Fn(f64) -> f64>, ThermoError>;
    fn get_dh_fun(&self) -> Result<Box<dyn Fn(f64) -> f64>, ThermoError>;
    fn get_ds_fun(&self) -> Result<Box<dyn Fn(f64) -> f64>, ThermoError>;
    fn get_Cp_sym(&self) -> Result<Expr, ThermoError>;
    fn get_dh_sym(&self) -> Result<Expr, ThermoError>;
    fn get_ds_sym(&self) -> Result<Expr, ThermoError>;
}
#[derive(Clone, Debug)]
#[enum_dispatch(ThermoCalculator)]
pub enum ThermoEnum {
    NIST(super::NISTdata::NISTdata),
    NASA(super::NASAdata::NASAdata),
}

pub enum ThermoType {
    NIST,
    NASA,
}

pub fn create_thermal(calc_type: ThermoType) -> ThermoEnum {
    match calc_type {
        ThermoType::NIST => ThermoEnum::NIST(super::NISTdata::NISTdata::new()),
        ThermoType::NASA => ThermoEnum::NASA(super::NASAdata::NASAdata::new()),
    }
}

pub fn create_thermal_by_name(calc_name: &str) -> ThermoEnum {
    match calc_name {
        "Cantera_nasa_base_gas"
        | "NASA_gas"
        | "NASA_cond"
        | "Cantera_nasa_base_cond"
        | "nuig_thermo"
        | "NASA"
        | "NASA7" => ThermoEnum::NASA(super::NASAdata::NASAdata::new()),
        "NIST" | "NIST9" => ThermoEnum::NIST(super::NISTdata::NISTdata::new()),
        _ => panic!("no such library!"),
    }
}
pub fn energy_dimension(unit: EnergyUnit) -> String {
    match unit {
        EnergyUnit::J => "J".to_owned(),
        EnergyUnit::Cal => "cal".to_owned(),
    }
}
