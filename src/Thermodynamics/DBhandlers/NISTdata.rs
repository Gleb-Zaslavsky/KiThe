use crate::Thermodynamics::DBhandlers::NIST_parser::{
    NistInput, calculate_cp, calculate_dh, calculate_s,
};
use crate::Thermodynamics::DBhandlers::NIST_parser::{NistParser, Phase, SearchType};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use serde_json::Value;
use std::collections::HashMap;
use std::fmt;
use std::{error::Error, fmt::Debug};

#[derive(Debug)]
pub enum NISTError {
    NoCoefficientsFound { temperature: f64, range: String },
    InvalidTemperatureRange,
    SerdeError(serde_json::Error),
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
        }
    }
}

impl Error for NISTError {}

impl From<serde_json::Error> for NISTError {
    fn from(err: serde_json::Error) -> Self {
        NISTError::SerdeError(err)
    }
}

pub struct NISTdata {
    /// data parsed from library
    pub input: NistInput,
    /// Joules or Cal
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
        self.input = serde_json::from_value(serde).map_err(|e| NISTError::SerdeError(e))?;
        Ok(())
    }
    pub fn extract_coefficients(&mut self, T: f64) -> Result<(), NISTError> {
        let _ = self.input.extract_coefficients(T);
        Ok(())
    }
    pub fn create_closures_Cp_dH_dS(&mut self) {
        let (C_fun, dh_fun, ds_fun) = self
            .input
            .create_closure_cp_dh_ds()
            .expect("Error calculating cp, dh, ds");
        (self.C_fun, self.dh_fun, self.ds_fun) = (C_fun, dh_fun, ds_fun);
    }
    pub fn create_sym_Cp_dH_dS(&mut self) {
        (self.Cp_sym, self.dh_sym, self.ds_sym) = self
            .input
            .create_sym_cp_dh_ds()
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
                Ok(())
            }
            Err(e) => Err(NISTError::NoCoefficientsFound {
                temperature: 0.0, // Default temperature since we don't have a specific one
                range: format!("Failed to get data for {}: {}", sub_name, e),
            }),
        }
    }
}

impl Clone for NISTdata {
    fn clone(&self) -> Self {
        // Clone the necessary data upfront
        // Clone all necessary data upfront

        Self {
            input: self.input.clone(),
            unit: self.unit.clone(),
            Cp: self.Cp,
            dh: self.dh,
            ds: self.ds,
            C_fun: {
                let T_ranges = self.input.T.clone().unwrap();
                let cp_coeffs = self.input.cp.clone().unwrap();
                let unit_multiplier = self.input.unit_multiplier;

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
                let unit_multiplier = self.input.unit_multiplier;
                let dh0 = self.input.dh.clone().unwrap();

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
                        unit_multiplier * (calculate_dh(T / 1000.0, a, b, c, d, e, f, g, h) + dh0)
                    } else {
                        0.0 // Return a default value if coefficients are not available
                    }
                })
            },
            ds_fun: {
                let T_ranges = self.input.T.clone().unwrap();
                let cp_coeffs = self.input.cp.clone().unwrap();
                let unit_multiplier = self.input.unit_multiplier;

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
        let (a, b, c, d, e, f, g, h) = self.input.coeffs.unwrap();
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
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::ThermoCalculator;
    use crate::Thermodynamics::DBhandlers::thermo_api::create_thermal_by_name;
    use approx::assert_relative_eq;

    #[test]
    fn test_thermo_calculator_nist() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);

        // Test newinstance
        //   assert!(nist.newinstance().is_ok());

        // Test from_serde

        // Test set_unit
        //   assert!(nist.set_unit(EnergyUnit::J).is_ok()) ;
        //   assert!(nist.set_unit(Some(EnergyUnit::Cal)).is_ok());

        // Test extract_model_coefficients
        // assert!(nist.extract_model_coefficients(400.0).is_ok());

        // Test calculate_Cp_dH_dS
        let _ = nist.extract_coefficients(400.0);
        nist.calculate_Cp_dH_dS(400.0);
        assert!(nist.Cp > 0.0);
        assert!(nist.dh != 0.0);
        assert!(nist.ds != 0.0);

        // Test create_closures_Cp_dH_dS
        nist.create_closures_Cp_dH_dS();
        let t = 400.0;
        assert_relative_eq!((nist.C_fun)(t), nist.Cp, epsilon = 1e-6);
        assert_relative_eq!((nist.dh_fun)(t), nist.dh, epsilon = 1e-6);
        assert_relative_eq!((nist.ds_fun)(t), nist.ds, epsilon = 1e-6);

        // Test create_sym_Cp_dH_dS
        nist.create_sym_Cp_dH_dS();
        let Cp_sym = &nist.Cp_sym;
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        assert_relative_eq!(Cp_value, nist.Cp, epsilon = 1e-6);

        // Test Taylor_series_cp_dh_ds
        let (Cp_taylor, dh_taylor, ds_taylor) = nist.Taylor_series_cp_dh_ds(400.0, 3).unwrap();
        assert!(!Cp_taylor.is_zero());
        assert!(!dh_taylor.is_zero());
        assert!(!ds_taylor.is_zero());
        nist.pretty_print_data().unwrap();
        // Test pretty_print_data
        assert!(nist.pretty_print_data().is_ok());
    }

    #[test]
    fn test_thermo_calculator_nist_error_handling() {
        let mut nist = NISTdata::new();

        // Test invalid serde data
        let invalid_data = serde_json::json!({
            "invalid": "data"
        });
        let result = nist.from_serde(invalid_data);
        assert!(result.is_err());
    }

    #[test]
    fn test_nist_clone() {
        let mut nist = NISTdata::new();

        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);

        let _ = nist.extract_coefficients(400.0);
        nist.calculate_Cp_dH_dS(400.0);
        nist.create_closures_Cp_dH_dS();
        nist.create_sym_Cp_dH_dS();

        // Clone the instance
        let mut nist_clone = nist.clone();
        let _ = nist_clone.extract_coefficients(400.0);
        // Test that the clone has the same values
        assert_relative_eq!(nist_clone.Cp, nist.Cp, epsilon = 1e-6);
        assert_relative_eq!(nist_clone.dh, nist.dh, epsilon = 1e-6);
        assert_relative_eq!(nist_clone.ds, nist.ds, epsilon = 1e-6);

        // Test that the cloned functions work
        let t = 400.0;
        assert_relative_eq!((nist_clone.C_fun)(t), (nist.C_fun)(t), epsilon = 1e-6);
        assert_relative_eq!((nist_clone.dh_fun)(t), (nist.dh_fun)(t), epsilon = 1e-6);
        assert_relative_eq!((nist_clone.ds_fun)(t), (nist.ds_fun)(t), epsilon = 1e-6);

        // Test that the symbolic expressions are the same
        assert_eq!(nist_clone.Cp_sym.to_string(), nist.Cp_sym.to_string());
        assert_eq!(nist_clone.dh_sym.to_string(), nist.dh_sym.to_string());
        assert_eq!(nist_clone.ds_sym.to_string(), nist.ds_sym.to_string());
    }
    #[test]
    fn ThermoCalculator_nist() {
        let mut nist = create_thermal_by_name("NIST");
        let _ = nist.newinstance();
        let _ = nist.renew_base("CO".to_owned(), SearchType::All, Phase::Gas);
        let T = 400.0;
        let _ = nist.extract_model_coefficients(T);
        let _ = nist.calculate_Cp_dH_dS(400.0);
        let Cp = nist.get_Cp().unwrap();
        let dh = nist.get_dh().unwrap();
        let ds = nist.get_ds().unwrap();
        assert!(Cp > 0.0);
        assert!(dh != 0.0);
        assert!(ds != 0.0);

        // Test create_closures_Cp_dH_dS
        let _ = nist.create_closures_Cp_dH_dS();
        let t = 400.0;
        assert_relative_eq!((nist.get_C_fun().unwrap())(t), Cp, epsilon = 1e-6);
        assert_relative_eq!((nist.get_dh_fun().unwrap())(t), dh, epsilon = 1e-6);
        assert_relative_eq!((nist.get_ds_fun().unwrap())(t), ds, epsilon = 1e-6);

        // Test create_sym_Cp_dH_dS
        let _ = nist.create_sym_Cp_dH_dS();
        let Cp_sym = &nist.get_Cp_sym().unwrap();
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        assert_relative_eq!(Cp_value, Cp, epsilon = 1e-6);

        // Test Taylor_series_cp_dh_ds
        let (Cp_taylor, dh_taylor, ds_taylor) = nist.Taylor_series_cp_dh_ds(400.0, 3).unwrap();
        assert!(!Cp_taylor.is_zero());
        assert!(!dh_taylor.is_zero());
        assert!(!ds_taylor.is_zero());

        // Test pretty_print_data
        assert!(nist.pretty_print_data().is_ok());
    }
}
