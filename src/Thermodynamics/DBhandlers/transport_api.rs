use RustedSciThe::symbolic::symbolic_engine::Expr;

use enum_dispatch::enum_dispatch;
use std::error::Error;
use std::fmt;
// Common error type for transport calculations
#[derive(Debug)]
pub enum TransportError {
    SerdeError(serde_json::Error),
    InvalidUnit(String),
    InvalidTemperature(f64),
    InvalidPressure(f64),
    InvalidMolarMass(f64),
    InvalidDensity(f64),
    CalculationError(String),
    MissingCoefficients(String),
    MissingData(&'static str),
    ParseError(String),
}

impl fmt::Display for TransportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TransportError::SerdeError(e) => write!(f, "Serde error: {}", e),
            TransportError::InvalidUnit(unit) => write!(f, "Invalid unit: {}", unit),
            TransportError::InvalidTemperature(t) => write!(f, "Invalid temperature: {}", t),
            TransportError::InvalidPressure(p) => write!(f, "Invalid pressure: {}", p),
            TransportError::InvalidMolarMass(m) => write!(f, "Invalid molar mass: {}", m),
            TransportError::InvalidDensity(d) => write!(f, "Invalid density: {}", d),
            TransportError::CalculationError(msg) => write!(f, "Calculation error: {}", msg),
            TransportError::MissingCoefficients(msg) => write!(f, "Missing coefficients: {}", msg),
            TransportError::MissingData(field) => write!(f, "Missing required data: {}", field),
            TransportError::ParseError(field) => write!(f, "Missing required data: {}", field),
        }
    }
}

impl Error for TransportError {}

impl From<super::TRANSPORTdata::TransportError> for TransportError {
    fn from(err: super::TRANSPORTdata::TransportError) -> Self {
        match err {
            super::TRANSPORTdata::TransportError::InvalidUnit(unit) => {
                TransportError::InvalidUnit(unit)
            }
            super::TRANSPORTdata::TransportError::SerdeError(e) => TransportError::SerdeError(e),
            super::TRANSPORTdata::TransportError::InvalidTemperature(t) => {
                TransportError::InvalidTemperature(t)
            }
            super::TRANSPORTdata::TransportError::InvalidPressure(p) => {
                TransportError::InvalidPressure(p)
            }
            super::TRANSPORTdata::TransportError::InvalidMolarMass(m) => {
                TransportError::InvalidMolarMass(m)
            }
            super::TRANSPORTdata::TransportError::InvalidDensity(d) => {
                TransportError::InvalidDensity(d)
            }
            super::TRANSPORTdata::TransportError::CalculationError(msg) => {
                TransportError::CalculationError(msg)
            }
        }
    }
}

impl From<super::CEAdata::CEAError> for TransportError {
    fn from(err: super::CEAdata::CEAError) -> Self {
        match err {
            super::CEAdata::CEAError::InvalidUnit(unit) => TransportError::InvalidUnit(unit),
            super::CEAdata::CEAError::SerdeError(e) => TransportError::SerdeError(e),
            super::CEAdata::CEAError::InvalidTemperature(t) => {
                TransportError::InvalidTemperature(t)
            }
            super::CEAdata::CEAError::MissingCoefficients(p) => {
                TransportError::MissingCoefficients(p)
            }
            super::CEAdata::CEAError::ParseError(m) => TransportError::ParseError(m),
            super::CEAdata::CEAError::MissingData(d) => TransportError::MissingData(d),
        }
    }
}
// Common unit types
#[derive(Debug, Clone, Copy)]
pub enum LambdaUnit {
    WPerMK,
    MWPerMK,
    MKWPerMK,
    MKWPerSMK,
}

#[derive(Debug, Clone, Copy)]
pub enum ViscosityUnit {
    KgPerMS,
    PaS,
    MKPaS,
}

// Common trait for transport calculations
#[enum_dispatch]
pub trait TransportCalculator {
    fn extract_coefficients(&mut self, t: f64) -> Result<(), TransportError>;
    // Basic calculations
    fn calculate_lambda(
        &mut self,
        C: Option<f64>,
        ro: Option<f64>,
        t: f64,
    ) -> Result<f64, TransportError>;
    fn calculate_viscosity(&mut self, t: f64) -> Result<f64, TransportError>;

    // Unit settings
    fn set_lambda_unit(&mut self, unit: Option<LambdaUnit>) -> Result<(), TransportError>;
    fn set_viscosity_unit(&mut self, unit: Option<ViscosityUnit>) -> Result<(), TransportError>;

    // Closure creation
    fn create_lambda_closure(
        &mut self,
        C: Option<f64>,
        ro: Option<f64>,
    ) -> Result<Box<dyn Fn(f64) -> f64>, TransportError>;
    fn create_viscosity_closure(&mut self) -> Result<Box<dyn Fn(f64) -> f64>, TransportError>;

    // Symbolic calculations
    fn create_symbolic_lambda(
        &mut self,
        C: Option<Expr>,
        ro: Option<Expr>,
    ) -> Result<(), TransportError>;
    fn create_symbolic_viscosity(&mut self) -> Result<(), TransportError>;
    fn taylor_series_lambda(
        &mut self,
        C: Option<Expr>,
        ro: Option<Expr>,
        t0: f64,
        n: usize,
    ) -> Result<Expr, TransportError>;

    // Data loading
    fn from_serde(&mut self, data: serde_json::Value) -> Result<(), TransportError>;
    fn set_M(&mut self, M: f64, M_unit: Option<String>) -> Result<(), TransportError>;
    fn set_P(&mut self, P: f64, P_unit: Option<String>) -> Result<(), TransportError>;
}

// Helper functions for unit conversion

pub fn lambda_dimension(unit: LambdaUnit) -> String {
    match unit {
        LambdaUnit::WPerMK => "W/m/K".to_owned(),
        LambdaUnit::MWPerMK => "mW/m/K".to_owned(),
        LambdaUnit::MKWPerMK => "mkW/m/K".to_owned(),
        LambdaUnit::MKWPerSMK => "mkW/sm/K".to_owned(),
    }
}
pub fn viscosity_dimension(unit: ViscosityUnit) -> String {
    match unit {
        ViscosityUnit::KgPerMS => "kg/m/s".to_owned(),
        ViscosityUnit::PaS => "Pa*s".to_owned(),
        ViscosityUnit::MKPaS => "mkPa*s".to_owned(),
    }
}

// Helper functions for validation
pub fn validate_temperature(t: f64) -> Result<(), TransportError> {
    if t <= 0.0 {
        Err(TransportError::InvalidTemperature(t))
    } else {
        Ok(())
    }
}

pub fn validate_pressure(p: f64) -> Result<(), TransportError> {
    if p <= 0.0 {
        Err(TransportError::InvalidPressure(p))
    } else {
        Ok(())
    }
}

pub fn validate_molar_mass(m: f64) -> Result<(), TransportError> {
    if m <= 0.0 {
        Err(TransportError::InvalidMolarMass(m))
    } else {
        Ok(())
    }
}

pub fn validate_density(d: f64) -> Result<(), TransportError> {
    if d <= 0.0 {
        Err(TransportError::InvalidDensity(d))
    } else {
        Ok(())
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
// FACTORY METHODS  ////////////////////////////////////////////////////////////////////
#[derive(Debug, Clone)]
#[enum_dispatch(TransportCalculator)]
pub enum TransportEnum {
    CEA(super::CEAdata::CEAdata),
    Collision(super::TRANSPORTdata::TransportData),
}

pub enum TransportType {
    CEA,
    Collision,
}

pub fn create_transport_calculator(calc_type: TransportType) -> TransportEnum {
    match calc_type {
        TransportType::CEA => TransportEnum::CEA(super::CEAdata::CEAdata::new()),
        TransportType::Collision => {
            TransportEnum::Collision(super::TRANSPORTdata::TransportData::new())
        }
    }
}

pub fn create_transport_calculator_by_name(calc_name: &str) -> TransportEnum {
    match calc_name {
        "CEA" => TransportEnum::CEA(super::CEAdata::CEAdata::new()),
        "Aramco_transpot" => TransportEnum::Collision(super::TRANSPORTdata::TransportData::new()),
        _ => panic!("No such library!"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::thermo_lib_api::ThermoData;
    // use approx::assert_relative_eq;

    #[test]
    fn test_unit_conversion() {}

    #[test]
    fn test_validation() {
        assert!(validate_temperature(300.0).is_ok());
        assert!(validate_temperature(-1.0).is_err());
        assert!(validate_temperature(0.0).is_err());

        assert!(validate_pressure(101325.0).is_ok());
        assert!(validate_pressure(-1.0).is_err());
        assert!(validate_pressure(0.0).is_err());

        assert!(validate_molar_mass(28.0).is_ok());
        assert!(validate_molar_mass(-1.0).is_err());
        assert!(validate_molar_mass(0.0).is_err());

        assert!(validate_density(1.225).is_ok());
        assert!(validate_density(-1.0).is_err());
        assert!(validate_density(0.0).is_err());
    }

    #[test]
    fn test_calculator_creation() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let _CO_data = sublib.get("CO").unwrap();
        /*
        let mut calculator = create_transport_calculator("TRANSPORT");
        calculator.from_serde(CO_data.clone()).unwrap();
        calculator.set_lambda_unit(Some(LambdaUnit::MWPerMK)).unwrap();
        calculator.set_viscosity_unit(Some(ViscosityUnit::MKPaS)).unwrap();

        let lambda = calculator.calculate_lambda(500.0).unwrap();
        let visc = calculator.calculate_viscosity(500.0).unwrap();

        assert!(lambda > 0.0);
        assert!(visc > 0.0);
        */
    }
    #[test]
    fn test_factory_method_CEA() {
        let thermo_data = ThermoData::new();

        // Test CEA calculator creation
        let mut cea_calc = create_transport_calculator(TransportType::CEA);
        let sublib = thermo_data.LibThermoData.get("CEA").unwrap();
        let co_data = sublib.get("CO").unwrap();
        cea_calc.from_serde(co_data.clone()).unwrap();
        cea_calc.set_lambda_unit(Some(LambdaUnit::MWPerMK)).unwrap();
        cea_calc
            .set_viscosity_unit(Some(ViscosityUnit::MKPaS))
            .unwrap();

        // Test calling methods on the enum
        cea_calc.extract_coefficients(500.0).unwrap();
        let lambda_cea = cea_calc.calculate_lambda(None, None, 500.0).unwrap();
        let visc_cea = cea_calc.calculate_viscosity(500.0).unwrap();
        assert!(lambda_cea > 0.0);
        assert!(visc_cea > 0.0);
    }
    #[test]
    fn test_factory_method_Collision() {
        // Test Collision calculator creation
        let thermo_data = ThermoData::new();
        let mut collision_calc = create_transport_calculator(TransportType::Collision);

        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let co_data = sublib.get("CO").unwrap();
        let T = 500.0;
        let _ = collision_calc.set_M(0.028, None);
        let _ = collision_calc.set_P(101325.0, Some("Pa".to_owned()));
        collision_calc.from_serde(co_data.clone()).unwrap();
        collision_calc
            .set_lambda_unit(Some(LambdaUnit::MWPerMK))
            .unwrap();
        collision_calc
            .set_viscosity_unit(Some(ViscosityUnit::MKPaS))
            .unwrap();
        let _ = collision_calc.extract_coefficients(T);
        // we need Cp info lets take it from library___________________________________________
        let sublib = thermo_data.LibThermoData.get("nuig_thermo").unwrap();
        let co_data = sublib.get("CO").unwrap();
        use crate::Thermodynamics::DBhandlers::NASAdata::NASAdata;
        let mut nasa = NASAdata::new();
        let _ = nasa.from_serde(co_data.clone());
        let _ = nasa.calculate_Cp_dH_dS(T);
        //  let ro = P;

        let Cp = nasa.Cp;
        //______________________________________________________________________
        let lambda_collision = collision_calc
            .calculate_lambda(Some(Cp), None, 500.0)
            .unwrap();
        let visc_collision = collision_calc.calculate_viscosity(T).unwrap();
        assert!(lambda_collision > 0.0);
        assert!(visc_collision > 0.0);
    }
}
