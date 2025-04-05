use RustedSciThe::symbolic::symbolic_engine::Expr;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fmt;

use serde_json::Value;
//use super::transport_api::{TransportCalculator, TransportError, LambdaUnit, ViscosityUnit};
//use super::transport_api::{validate_temperature, validate_pressure, validate_molar_mass, validate_density};
//use super::transport_api::{lambda_unit_to_multiplier, viscosity_unit_to_multiplier};

#[derive(Debug)]
pub enum CEAError {
    InvalidUnit(String),
    InvalidTemperature(f64),
    MissingCoefficients(String),
    ParseError(String),
    SerdeError(serde_json::Error),
    MissingData(&'static str),
}

impl fmt::Display for CEAError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CEAError::InvalidUnit(unit) => write!(f, "Invalid unit: {}", unit),
            CEAError::InvalidTemperature(temp) => write!(f, "Temperature {} is out of range", temp),
            CEAError::MissingCoefficients(msg) => write!(f, "Missing coefficients: {}", msg),
            CEAError::ParseError(msg) => write!(f, "Failed to parse data: {}", msg),
            CEAError::SerdeError(e) => write!(f, "Serde error: {}", e),
            CEAError::MissingData(field) => write!(f, "Missing required data: {}", field),
        }
    }
}

impl Error for CEAError {}

impl From<serde_json::Error> for CEAError {
    fn from(err: serde_json::Error) -> Self {
        CEAError::SerdeError(err)
    }
}

fn calculate_L(t: f64, e: f64, f: f64, k: f64, g: f64) -> f64 {
    //Вт/м*К
    (1E-4) * (e * t.ln() + f / t + k / t.powi(2) + g).exp()
}

fn calculate_V(t: f64, e: f64, f: f64, k: f64, g: f64) -> f64 {
    // кг/м*с
    (1E-7) * (e * t.ln() + f / t + k / t.powi(2) + g).exp()
}

fn calculate_L_sym(e: f64, f: f64, k: f64, g: f64) -> Expr {
    let (e, f, k, g) = (
        Expr::Const(e),
        Expr::Const(f),
        Expr::Const(k),
        Expr::Const(g),
    );
    let t = Expr::Var("T".to_owned());
    Expr::Const(1E-4)
        * (e * t.clone().ln() + f / t.clone() + k / t.clone().pow(Expr::Const(2.0)) + g).exp()
}
fn calculate_V_sym(e: f64, f: f64, k: f64, g: f64) -> Expr {
    let (e, f, k, g) = (
        Expr::Const(e),
        Expr::Const(f),
        Expr::Const(k),
        Expr::Const(g),
    );
    let t = Expr::Var("T".to_owned());
    Expr::Const(1E-7)
        * (e * t.clone().ln() + f / t.clone() + k / t.clone().pow(Expr::Const(2.0)) + g).exp()
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct CEAinput {
    pub LC: Vec<String>,
    pub model: Vec<f64>,
}

pub struct CEAdata {
    input: CEAinput,
    pub L_unit: Option<String>,
    ///
    L_unit_multiplier: f64,
    pub V_unit: Option<String>,
    V_unit_multiplier: f64,

    pub coeffs: HashMap<String, HashMap<String, Vec<f64>>>,
    pub coeff_Visc: Vec<f64>,
    pub coeff_Lambda: Vec<f64>,
    pub Lambda: f64,
    pub V: f64,

    pub Lambda_fun: Box<dyn Fn(f64) -> f64>,
    pub V_fun: Box<dyn Fn(f64) -> f64>,

    pub Lambda_sym: Option<Expr>,
    pub V_sym: Option<Expr>,
}

impl CEAdata {
    pub fn new() -> Self {
        let input = CEAinput {
            LC: vec!["".to_string()],
            model: Vec::new(),
        };
        Self {
            input: input,
            L_unit: None,
            L_unit_multiplier: 1.0,
            V_unit: None,
            V_unit_multiplier: 1.0,
            coeffs: HashMap::new(),
            coeff_Visc: Vec::new(),
            coeff_Lambda: Vec::new(),

            Lambda: 0.0,
            V: 0.0,
            Lambda_fun: Box::new(|x| x),
            V_fun: Box::new(|x| x),
            Lambda_sym: None,
            V_sym: None,
        }
    }

    pub fn set_input(&mut self, input: CEAinput) {
        self.input = input;
        let _ = self.parse_coefficients();
    }

    pub fn set_lambda_unit(&mut self, unit: Option<String>) -> Result<(), CEAError> {
        if let Some(unit) = unit {
            self.L_unit_multiplier = match unit.as_str() {
                "W/m/K" => 1.0,
                "mW/m/K" => 1E3,
                "mkW/m/K" => 1E6,
                "mkW/sm/K" => 1E+4,
                _ => return Err(CEAError::InvalidUnit(unit)),
            };
            self.L_unit = Some(unit);
        }
        Ok(())
    }

    pub fn set_V_unit(&mut self, unit: Option<String>) -> Result<(), CEAError> {
        if let Some(unit) = unit {
            self.V_unit_multiplier = match unit.as_str() {
                "kg/m/s" => 1.0,
                "Pa*s" => 1.0,
                "mkPa*s" => 1E6,
                _ => return Err(CEAError::InvalidUnit(unit)),
            };
            self.V_unit = Some(unit);
        }
        Ok(())
    }

    pub fn from_serde(&mut self, serde: Value) -> Result<(), CEAError> {
        self.input = serde_json::from_value(serde).map_err(|e| CEAError::SerdeError(e))?;
        let _ = self.parse_coefficients();
        Ok(())
    }

    pub fn parse_coefficients(&mut self) -> Result<(), CEAError> {
        let LC = &self.input.LC;
        let model = &self.input.model;

        if LC.is_empty() {
            return Err(CEAError::MissingData("LC is empty"));
        }

        let mut V: HashMap<String, HashMap<String, Vec<f64>>> = HashMap::new();
        let mut i = 0;
        let mut k: usize = 0;

        while i < LC.len() {
            let lc = &LC[i];
            let elements_to_take = 6;

            if k + elements_to_take > model.len() {
                return Err(CEAError::ParseError(format!(
                    "Not enough data in model for LC entry {}",
                    lc
                )));
            }

            let start: usize = elements_to_take - 4;
            let (T0, T1) = (model[k], model[k + 1]);
            let coefficients = &model[k + start..k + elements_to_take];

            V.entry(lc.to_string())
                .or_insert_with(HashMap::new)
                .entry("T".to_string())
                .or_insert_with(Vec::new)
                .extend_from_slice(&[T0, T1]);

            V.entry(lc.to_string())
                .or_insert_with(HashMap::new)
                .entry("coeff".to_string())
                .or_insert_with(Vec::new)
                .extend_from_slice(coefficients);

            k += elements_to_take;
            i += 1;
        }

        self.coeffs = V;
        Ok(())
    }

    pub fn extract_coefficients(&mut self, t: f64) -> Result<(), CEAError> {
        // Extract viscosity coefficients
        let V = self
            .coeffs
            .get("V")
            .ok_or_else(|| CEAError::MissingData("Viscosity data not found"))?;
        let T = V
            .get("T")
            .ok_or_else(|| CEAError::MissingData("Temperature range for viscosity not found"))?;
        let coeffs = V
            .get("coeff")
            .ok_or_else(|| CEAError::MissingData("Viscosity coefficients not found"))?;

        let i = T
            .windows(2)
            .position(|pair| pair[0] <= t && pair[1] >= t)
            .ok_or_else(|| CEAError::InvalidTemperature(t))?
            / 2;

        if 4 * i + 4 > coeffs.len() {
            return Err(CEAError::MissingCoefficients(
                "Insufficient viscosity coefficients".to_string(),
            ));
        }

        self.coeff_Visc = coeffs[4 * i..4 * i + 4].to_vec();

        // Extract thermal conductivity coefficients
        let L = self
            .coeffs
            .get("C")
            .ok_or_else(|| CEAError::MissingData("Thermal conductivity data not found"))?;
        let T = L
            .get("T")
            .ok_or_else(|| CEAError::MissingData("Temperature range for conductivity not found"))?;
        let coeffs = L
            .get("coeff")
            .ok_or_else(|| CEAError::MissingData("Conductivity coefficients not found"))?;

        let i = T
            .windows(2)
            .position(|pair| pair[0] <= t && pair[1] >= t)
            .ok_or_else(|| CEAError::InvalidTemperature(t))?
            / 2;

        if 4 * i + 4 > coeffs.len() {
            return Err(CEAError::MissingCoefficients(
                "Insufficient conductivity coefficients".to_string(),
            ));
        }

        self.coeff_Lambda = coeffs[4 * i..4 * i + 4].to_vec();
        Ok(())
    }

    pub fn calculate_Lambda(&mut self, t: f64) -> Result<f64, CEAError> {
        if self.coeff_Lambda.len() != 4 {
            return Err(CEAError::MissingCoefficients(
                "Lambda coefficients not properly initialized".to_string(),
            ));
        }
        let c = self.coeff_Lambda.clone();
        let (e, f, k, g) = (c[0], c[1], c[2], c[3]);
        let Lambda = self.L_unit_multiplier * calculate_L(t, e, f, k, g);
        self.Lambda = Lambda;
        Ok(Lambda)
    }

    pub fn calculate_Visc(&mut self, t: f64) -> Result<f64, CEAError> {
        if self.coeff_Visc.len() != 4 {
            return Err(CEAError::MissingCoefficients(
                "Viscosity coefficients not properly initialized".to_string(),
            ));
        }
        let c = self.coeff_Visc.clone();
        let (e, f, k, g) = (c[0], c[1], c[2], c[3]);
        let Visc = self.V_unit_multiplier * calculate_V(t, e, f, k, g);
        self.V = Visc;
        Ok(Visc)
    }

    pub fn create_closure_Lambda(&mut self) -> Result<Box<dyn Fn(f64) -> f64>, CEAError> {
        if self.coeff_Lambda.len() != 4 {
            return Err(CEAError::MissingCoefficients(
                "Lambda coefficients not properly initialized".to_string(),
            ));
        }
        let c = self.coeff_Lambda.clone();
        let um = self.L_unit_multiplier;
        let (e, f, k, g) = (c[0], c[1], c[2], c[3]);
        let Lambda = move |t: f64| um * calculate_L(t, e, f, k, g);
        self.Lambda_fun = Box::new(Lambda.clone());
        Ok(Box::new(Lambda))
    }

    pub fn create_closure_Visc(&mut self) -> Result<Box<dyn Fn(f64) -> f64>, CEAError> {
        if self.coeff_Visc.len() != 4 {
            return Err(CEAError::MissingCoefficients(
                "Viscosity coefficients not properly initialized".to_string(),
            ));
        }
        let c = self.coeff_Visc.clone();
        let um = self.V_unit_multiplier;
        let (e, f, k, g) = (c[0], c[1], c[2], c[3]);
        let V = move |t: f64| um * calculate_V(t, e, f, k, g);
        self.V_fun = Box::new(V.clone());
        Ok(Box::new(V))
    }

    pub fn create_sym_Lambda(&mut self) -> Result<(), CEAError> {
        if self.coeff_Lambda.len() != 4 {
            return Err(CEAError::MissingCoefficients(
                "Lambda coefficients not properly initialized".to_string(),
            ));
        }
        let c = self.coeff_Lambda.clone();
        let um = self.L_unit_multiplier;
        let (e, f, k, g) = (c[0], c[1], c[2], c[3]);
        let L_sym = (Expr::Const(um) * calculate_L_sym(e, f, k, g)).symplify();
        self.Lambda_sym = Some(L_sym);
        Ok(())
    }

    pub fn create_sym_Visc(&mut self) -> Result<(), CEAError> {
        if self.coeff_Visc.len() != 4 {
            return Err(CEAError::MissingCoefficients(
                "Viscosity coefficients not properly initialized".to_string(),
            ));
        }
        let c = self.coeff_Visc.clone();
        let um = self.V_unit_multiplier;
        let (e, f, k, g) = (c[0], c[1], c[2], c[3]);
        let V_sym = (Expr::Const(um) * calculate_V_sym(e, f, k, g)).symplify();
        self.V_sym = Some(V_sym);
        Ok(())
    }
    pub fn Taylor_series_Lambda(&mut self, T0: f64, n: usize) -> Result<Expr, CEAError> {
        if T0 <= 0.0 {
            return Err(CEAError::InvalidTemperature(T0));
        }
        self.create_sym_Lambda()?;
        let Lambda_series = self
            .Lambda_sym
            .clone()
            .ok_or_else(|| CEAError::MissingCoefficients("Lambda_sym not calculated".to_string()))?
            .taylor_series1D_("T", T0, n);
        Ok(Lambda_series)
    }
}

impl Clone for CEAdata {
    fn clone(&self) -> Self {
        Self {
            input: self.input.clone(),
            L_unit: self.L_unit.clone(),
            L_unit_multiplier: self.L_unit_multiplier,
            V_unit: self.V_unit.clone(),
            V_unit_multiplier: self.V_unit_multiplier,
            coeffs: self.coeffs.clone(),
            coeff_Visc: self.coeff_Visc.clone(),
            coeff_Lambda: self.coeff_Lambda.clone(),
            Lambda: self.Lambda,
            V: self.V,
            Lambda_fun: {
                let c = self.coeff_Lambda.clone();
                let um = self.L_unit_multiplier;
                let (e, f, k, g) = (c[0], c[1], c[2], c[3]);
                let Lambda = move |t: f64| um * calculate_L(t, e, f, k, g);

                Box::new(Lambda)
            },
            V_fun: {
                let c = self.coeff_Visc.clone();
                let um = self.V_unit_multiplier;
                let (e, f, k, g) = (c[0], c[1], c[2], c[3]);
                let V = move |t: f64| um * calculate_V(t, e, f, k, g);
                Box::new(V)
            },
            Lambda_sym: self.Lambda_sym.clone(),
            V_sym: self.V_sym.clone(),
        }
    }
}

impl fmt::Debug for CEAdata {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CEAdata")
            .field("input", &self.input)
            .field("L_unit", &self.L_unit)
            .field("L_unit_multiplier", &self.L_unit_multiplier)
            .field("V_unit", &self.V_unit)
            .field("V_unit_multiplier", &self.V_unit_multiplier)
            .field("coeffs", &self.coeffs)
            .field("coeff_Visc", &self.coeff_Visc)
            .field("coeff_Lambda", &self.coeff_Lambda)
            .field("Lambda", &self.Lambda)
            .field("V", &self.V)
            .field("Lambda_fun", &"<closure>")
            .field("V_fun", &"<closure>")
            .field("Lambda_sym", &self.Lambda_sym)
            .field("V_sym", &self.V_sym)
            .finish_non_exhaustive() //  The finish_non_exhaustive() method is used to indicate that not all fields are being displayed.
    }
}

use super::transport_api::{LambdaUnit, TransportCalculator, ViscosityUnit};
use super::transport_api::{lambda_dimension, validate_temperature, viscosity_dimension};
impl TransportCalculator for CEAdata {
    fn extract_coefficients(&mut self, t: f64) -> Result<(), super::transport_api::TransportError> {
        self.parse_coefficients()?;
        self.extract_coefficients(t)?;
        Ok(())
    }
    fn calculate_lambda(
        &mut self,
        _C: Option<f64>,
        _ro: Option<f64>,
        T: f64,
    ) -> Result<f64, super::transport_api::TransportError> {
        validate_temperature(T)?;
        let Lambda = self.calculate_Lambda(T)?;
        self.Lambda = Lambda;
        Ok(Lambda)
    }

    fn calculate_viscosity(&mut self, T: f64) -> Result<f64, super::transport_api::TransportError> {
        validate_temperature(T)?;
        let eta = self.calculate_Visc(T)?;
        self.V = eta;
        Ok(eta)
    }

    fn set_lambda_unit(
        &mut self,
        unit: Option<LambdaUnit>,
    ) -> Result<(), super::transport_api::TransportError> {
        if let Some(unit) = unit {
            self.L_unit = Some(lambda_dimension(unit));
        }
        Ok(())
    }

    fn set_viscosity_unit(
        &mut self,
        unit: Option<ViscosityUnit>,
    ) -> Result<(), super::transport_api::TransportError> {
        if let Some(unit) = unit {
            self.V_unit = Some(viscosity_dimension(unit));
        }
        Ok(())
    }

    fn create_lambda_closure(
        &mut self,
        _C: Option<f64>,
        _ro: Option<f64>,
    ) -> Result<Box<dyn Fn(f64) -> f64>, super::transport_api::TransportError> {
        let Lambda_fun = self.create_closure_Lambda()?;
        Ok(Lambda_fun)
    }

    fn create_viscosity_closure(
        &mut self,
    ) -> Result<Box<dyn Fn(f64) -> f64>, super::transport_api::TransportError> {
        let visc = self.create_closure_Visc()?;

        Ok(visc)
    }

    fn create_symbolic_lambda(
        &mut self,
        _C: Option<Expr>,
        _ro: Option<Expr>,
    ) -> Result<(), super::transport_api::TransportError> {
        self.create_sym_Lambda()?;
        Ok(())
    }

    fn create_symbolic_viscosity(&mut self) -> Result<(), super::transport_api::TransportError> {
        self.create_sym_Visc()?;
        Ok(())
    }

    fn taylor_series_lambda(
        &mut self,
        _C: Option<Expr>,
        _ro: Option<Expr>,
        t0: f64,
        n: usize,
    ) -> Result<Expr, super::transport_api::TransportError> {
        validate_temperature(t0)?;
        let Lambda_series = self.Taylor_series_Lambda(t0, n)?;
        Ok(Lambda_series)
    }

    fn from_serde(
        &mut self,
        data: serde_json::Value,
    ) -> Result<(), super::transport_api::TransportError> {
        self.input = serde_json::from_value(data).map_err(|e| CEAError::SerdeError(e))?;
        Ok(())
    }

    fn set_M(
        &mut self,
        _M: f64,
        _M_unit: Option<String>,
    ) -> Result<(), super::transport_api::TransportError> {
        Ok(())
    }
    fn set_P(
        &mut self,
        _P: f64,
        _P_unit: Option<String>,
    ) -> Result<(), super::transport_api::TransportError> {
        Ok(())
    }
    fn print_instance(&self) -> Result<(), super::transport_api::TransportError> {
        println!("{:?}", &self);
        Ok(())
    }
    fn get_lambda_sym(&self) -> Result<Expr, super::transport_api::TransportError> {
        match &self.Lambda_sym {
            Some(lambda_sym) => Ok(lambda_sym.clone()),
            None => Err(super::transport_api::TransportError::CalculationError(
                "Lambda_sym not calculated".to_string(),
            )),
        }
    }
    fn get_viscosity_sym(&self) -> Result<Expr, super::transport_api::TransportError> {
        match &self.V_sym {
            Some(viscosity_sym) => Ok(viscosity_sym.clone()),
            None => Err(super::transport_api::TransportError::CalculationError(
                "V_sym not calculated".to_string(),
            )),
        }
    }
    fn get_lambda_fun(
        &self,
    ) -> Result<Box<dyn Fn(f64) -> f64>, super::transport_api::TransportError> {
        Ok(self.clone().Lambda_fun)
    }
    fn get_viscosity_fun(
        &self,
    ) -> Result<Box<dyn Fn(f64) -> f64>, super::transport_api::TransportError> {
        Ok(self.clone().V_fun)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::thermo_lib_api::ThermoData;
    use approx::assert_relative_eq;

    fn setup_cea_data() -> CEAdata {
        let mut CEA = CEAdata::new();
        let input = CEAinput {
            LC: vec![
                "V".to_string(),
                "V".to_string(),
                "C".to_string(),
                "C".to_string(),
            ],
            model: vec![
                200.0,
                1000.0,
                0.62526577,
                -31.779652,
                -1640.7983,
                1.7454992,
                1000.0,
                5000.0,
                0.87395209,
                561.52222,
                -173948.09,
                -0.39335958,
                200.0,
                1000.0,
                0.85439436,
                105.73224,
                -12347.848,
                0.47793128,
                1000.0,
                5000.0,
                0.88407146,
                133.57293,
                -11429.64,
                0.24417019,
            ],
        };
        CEA.set_input(input);
        CEA
    }

    #[test]
    fn test_parse_coefficients() {
        let mut CEA = setup_cea_data();
        assert!(CEA.parse_coefficients().is_ok());

        let v_coeffs = CEA.coeffs.get("V").unwrap();
        assert_eq!(
            v_coeffs.get("T").unwrap(),
            &vec![200.0, 1000.0, 1000.0, 5000.0]
        );
        assert_eq!(
            v_coeffs.get("coeff").unwrap(),
            &vec![
                0.62526577,
                -31.779652,
                -1640.7983,
                1.7454992,
                0.87395209,
                561.52222,
                -173948.09,
                -0.39335958
            ]
        );

        // Test error case
        let mut CEA = CEAdata::new();
        let input = CEAinput {
            LC: vec!["V".to_string()],
            model: vec![1.0], // Not enough data
        };
        CEA.set_input(input);
        assert!(matches!(
            CEA.parse_coefficients(),
            Err(CEAError::ParseError(_))
        ));
    }

    #[test]
    fn test_extract_coefficients() {
        let mut CEA = setup_cea_data();

        // Test valid temperature
        assert!(CEA.extract_coefficients(500.0).is_ok());
        assert_eq!(
            CEA.coeff_Visc,
            vec![0.62526577, -31.779652, -1640.7983, 1.7454992]
        );
        assert_eq!(
            CEA.coeff_Lambda,
            vec![0.85439436, 105.73224, -12347.848, 0.47793128]
        );

        // Test temperature out of range
        assert!(matches!(
            CEA.extract_coefficients(100.0),
            Err(CEAError::InvalidTemperature(_))
        ));
    }

    #[test]
    fn test_calculate_lambda_and_visc() {
        let mut CEA = setup_cea_data();
        CEA.extract_coefficients(500.0).unwrap();
        CEA.set_lambda_unit(Some("mW/m/K".to_string())).unwrap();
        CEA.set_V_unit(Some("mkPa*s".to_string())).unwrap();

        let lambda = CEA.calculate_Lambda(500.0).unwrap();
        let visc = CEA.calculate_Visc(500.0).unwrap();

        assert_relative_eq!(lambda, 39.2, epsilon = 1.0);
        assert_relative_eq!(visc, 25.8, epsilon = 1.0);

        // Test invalid units
        assert!(matches!(
            CEA.set_lambda_unit(Some("invalid".to_string())),
            Err(CEAError::InvalidUnit(_))
        ));
        assert!(matches!(
            CEA.set_V_unit(Some("invalid".to_string())),
            Err(CEAError::InvalidUnit(_))
        ));
    }

    #[test]
    fn test_create_closures() {
        let mut CEA = setup_cea_data();
        CEA.extract_coefficients(500.0).unwrap();
        CEA.set_lambda_unit(Some("mW/m/K".to_string())).unwrap();
        CEA.set_V_unit(Some("mkPa*s".to_string())).unwrap();

        let lambda_closure = CEA.create_closure_Lambda().unwrap();
        let visc_closure = CEA.create_closure_Visc().unwrap();

        assert_relative_eq!(lambda_closure(500.0), 39.2, epsilon = 1.);
        assert_relative_eq!(visc_closure(500.0), 25.8, epsilon = 1.0);

        // Test error case when coefficients are not initialized
        let mut CEA = CEAdata::new();
        assert!(matches!(
            CEA.create_closure_Lambda(),
            Err(CEAError::MissingCoefficients(_))
        ));
        assert!(matches!(
            CEA.create_closure_Visc(),
            Err(CEAError::MissingCoefficients(_))
        ));
    }

    #[test]
    fn test_create_symbolic_expressions() {
        let mut CEA = setup_cea_data();
        CEA.extract_coefficients(500.0).unwrap();
        CEA.set_lambda_unit(Some("mW/m/K".to_string())).unwrap();
        CEA.set_V_unit(Some("mkPa*s".to_string())).unwrap();

        assert!(CEA.create_sym_Lambda().is_ok());
        assert!(CEA.create_sym_Visc().is_ok());

        let lambda_sym = CEA.Lambda_sym.as_ref().unwrap();
        let lambda_val = lambda_sym.lambdify1D()(500.0);
        let visc_sym = CEA.V_sym.as_ref().unwrap();
        let visc_val = visc_sym.lambdify1D()(500.0);

        assert_relative_eq!(lambda_val, 39.2, epsilon = 1.);
        assert_relative_eq!(visc_val, 25.8, epsilon = 1.0);

        // Test error case when coefficients are not initialized
        let mut CEA = CEAdata::new();
        assert!(matches!(
            CEA.create_sym_Lambda(),
            Err(CEAError::MissingCoefficients(_))
        ));
        assert!(matches!(
            CEA.create_sym_Visc(),
            Err(CEAError::MissingCoefficients(_))
        ));
    }

    #[test]
    fn test_with_real_data() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("CEA").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut CEA = CEAdata::new();
        CEA.from_serde(CO_data.clone()).unwrap();
        CEA.set_lambda_unit(Some("mW/m/K".to_string())).unwrap();
        CEA.set_V_unit(Some("mkPa*s".to_string())).unwrap();
        CEA.parse_coefficients().unwrap();
        CEA.extract_coefficients(500.0).unwrap();
        let lambda = CEA.calculate_Lambda(500.0).unwrap();
        let visc = CEA.calculate_Visc(500.0).unwrap();
        println!("Lambda, mW/m/K: {:?}, Visc: {:?}", lambda, visc);
    }
}
