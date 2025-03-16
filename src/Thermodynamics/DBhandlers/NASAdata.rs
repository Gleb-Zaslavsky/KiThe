use RustedSciThe::symbolic::symbolic_engine::Expr;
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

#[derive(Debug)]
pub enum NASAError {
    NoCoefficientsFound { temperature: f64, range: String },
    InvalidTemperatureRange,
    DeserializationError(String),
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
            NASAError::DeserializationError(msg) => {
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
        self.input = serde_json::from_value(serde)
            .map_err(|e| NASAError::DeserializationError(e.to_string()))?;
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
    pub fn Taylor_series_Cp_dH_dS(&mut self, T0: f64) -> (Expr, Expr, Expr) {
        let T = Expr::Var("T".to_owned());
        let T0exp = Expr::Const(T0);
        let Cp = self.Cp_sym.clone();
        let Cp_at_T0 = Cp.lambdify1D()(T0);
        let Cp_at_T0_sym = Expr::Const(Cp_at_T0);
        let dC_dt = Cp.diff("T");
        let dC_dt_at_T0 = dC_dt.lambdify1D()(T0);
        let dC_dt_at_T0_sym = Expr::Const(dC_dt_at_T0);
        let Cp_taylor = Cp_at_T0_sym + dC_dt_at_T0_sym * (T0exp.clone() - T.clone());

        let dh = self.dh_sym.clone();
        let dh_at_T0 = dh.lambdify1D()(T0);
        let dh_at_T0_sym = Expr::Const(dh_at_T0);
        let dh_dt = dh.diff("T");
        let dh_dt_at_T0 = dh_dt.lambdify1D()(T0);
        let dh_dt_at_T0_sym = Expr::Const(dh_dt_at_T0);
        let dh_taylor = dh_at_T0_sym + dh_dt_at_T0_sym * (T0exp.clone() - T.clone());

        let ds = self.ds_sym.clone();
        let ds_at_T0 = ds.lambdify1D()(T0);
        let ds_at_T0_sym = Expr::Const(ds_at_T0);
        let ds_dt = ds.diff("T");
        let ds_dt_at_T0 = ds_dt.lambdify1D()(T0);
        let ds_dt_at_T0_sym = Expr::Const(ds_dt_at_T0);
        let ds_taylor = ds_at_T0_sym + ds_dt_at_T0_sym * (T0exp.clone() - T.clone());
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
}
