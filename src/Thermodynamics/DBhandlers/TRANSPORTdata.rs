//! # Transport Properties Data Handler Module
//!
//! ## Aim
//! This module handles transport properties (viscosity, thermal conductivity, diffusion) calculations
//! using collision theory and molecular kinetic theory. It implements the Aramco transport database
//! format and provides both numerical and symbolic calculations for gas transport properties.
//!
//! ## Main Data Structures and Logic
//! - `TransportData`: Main structure containing molecular parameters (sigma, epsilon, dipole moment)
//!   and calculated transport properties. Implements collision integrals and kinetic theory equations.
//! - `TransportInput`: Input structure for molecular parameters from database
//! - `TransportError`: Comprehensive error handling for invalid parameters and calculations
//!
//! ## Key Methods
//! - `calculate_Visc()`: Calculates dynamic viscosity using Chapman-Enskog theory
//! - `calculate_Lambda()`: Calculates thermal conductivity using kinetic theory with rotational/vibrational contributions
//! - `create_sym_*()`: Creates symbolic expressions for temperature-dependent properties
//! - `from_serde()`: Parses transport data from JSON database format
//!
//! ## Usage
//! ```rust, ignore
//! let mut transport = TransportData::new();
//! transport.from_serde(database_entry)?;
//! transport.set_M(28.0, Some("g/mol".to_string()))?;
//! transport.set_P(101325.0, Some("Pa".to_string()))?;
//! let viscosity = transport.calculate_Visc(500.0)?;
//! let conductivity = transport.calculate_Lambda(heat_capacity, density, 500.0)?;
//! ```
//!
//! ## Interesting Features
//! - Implements collision integrals (Ω₁₁, Ω₂₂) with dipole moment corrections
//! - Handles molecular geometry effects (linear vs non-linear molecules)
//! - Provides both numerical functions and symbolic expressions for integration with symbolic math
//! - Uses Hirschfelder molecular theory for accurate gas transport predictions
//! - Supports multiple unit systems with automatic conversion

use RustedSciThe::symbolic::symbolic_engine::Expr;

use serde::{Deserialize, Serialize};

use std::error::Error;
use std::fmt;

#[derive(Debug)]
pub enum TransportError {
    InvalidUnit(String),
    SerdeError(serde_json::Error),
    InvalidTemperature(f64),
    InvalidPressure(f64),
    InvalidMolarMass(f64),
    InvalidDensity(f64),
    CalculationError(String),
}

impl fmt::Display for TransportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TransportError::InvalidUnit(unit) => write!(f, "Invalid unit: {}", unit),
            TransportError::SerdeError(e) => write!(f, "Serde error: {}", e),
            TransportError::InvalidTemperature(t) => write!(f, "Invalid temperature value: {}", t),
            TransportError::InvalidPressure(p) => write!(f, "Invalid pressure value: {}", p),
            TransportError::InvalidMolarMass(m) => write!(f, "Invalid molar mass value: {}", m),
            TransportError::InvalidDensity(d) => write!(f, "Invalid density value: {}", d),
            TransportError::CalculationError(msg) => write!(f, "Calculation error: {}", msg),
        }
    }
}

impl Error for TransportError {}

impl From<serde_json::Error> for TransportError {
    fn from(err: serde_json::Error) -> Self {
        TransportError::SerdeError(err)
    }
}

use std::f64::consts::PI;
const K_B: f64 = 1.38e-23;
const R: f64 = 8.314;
#[allow(non_upper_case_globals)]
const R_sym: Expr = Expr::Const(R);
// Collision integral
// Гиршфельдер. Молекулярная теория газов и жидкостей
//Hirschfelder. The Molecular Theory of Gases and Liquids
fn omega_22_calc(T: f64, e_k: f64, mu: f64, sigma_k: f64) -> f64 {
    let T1 = T / e_k;
    let a1 = 1.16145;
    let b1 = 0.14874;
    let c1 = 0.52487;
    let d1 = 0.77320;
    let e1 = 2.16178;
    let f1 = 2.43787;

    let delta = 1e-19 * mu.powi(2) / (2.0 * e_k * K_B * sigma_k.powi(3));

    let Col =
        a1 / T1.powf(b1) + c1 / (d1 * T1).exp() + e1 / (f1 * T1).exp() + 0.2 * delta.powi(2) / T1;
    //  println!("T1: {}, Col: {}", T1, Col);
    Col
}

// collision integral for diffusion
// Гиршфельдер. Молекулярная теория газов и жидкостей
//Hirschfelder. The Molecular Theory of Gases and Liquids
fn omega_11_calc(T: f64, e_k: f64, mu: f64, sigma_k: f64) -> f64 {
    let a = 1.06036;
    let b = 0.15610;
    let c = 0.19300;
    let d = 0.47635;
    let e = 1.03587;
    let f = 1.52996;
    let g = 1.76474;
    let h = 3.89411;
    let T1 = T / e_k;
    let delta = 1e-19 * mu.powi(2) / (2.0 * e_k * K_B * sigma_k.powi(3));
    a / T1.powf(b)
        + c / (d * T1).exp()
        + e / (f * T1).exp()
        + g / (h * T1).exp()
        + 0.19 * delta.powi(2) / T1
}
// The properties of gases and liquids (McGraw-Hill chemical engineering series) 3rd Edition
// by Robert C. Reid, John M. Prausnitz, Thomas K. Sherwood (p. 346)
//  Рид, Праусниц, Шервуд. Свойства газов и жидкостей. 1982. стр 346
fn visc(M: f64, T: f64, e_k: f64, sigma_k: f64, mu: f64) -> f64 {
    // N*s/m2 to transfer to mkPoise *10^7
    let M = M * 1000.0; // from kg/mol to g/mol
    let omega_22 = omega_22_calc(T, e_k, mu, sigma_k);
    1e-7 * 26.69 * (M * T).sqrt() / (sigma_k.powi(2) * omega_22)
}

fn calculate_Lambda_(p: TransportInput, um: f64, M: f64, P: f64, C: f64, ro: f64, T: f64) -> f64 {
    let (form, sigma_k, e_k, rot_relax, mu) = (p.Form, p.diam, p.well_depth, p.rot_relax, p.dipole);
    let C = C - R; // from heat capacity at constant pressure to heat capacity at constant volume
    //   println!("form: {}, sigma_k: {}, e_k: {}, rot_relax: {}, mu: {}", form, sigma_k, e_k, rot_relax, mu);
    let (C_trans, C_rot, C_vib) = match form as i32 {
        1 => ((3.0 / 2.0) * R, R, C - (5.0 / 2.0) * R),
        2 => ((3.0 / 2.0) * R, (3.0 / 2.0) * R, C - 3.0 * R),
        _ => ((3.0 / 2.0) * R, 0.0, 0.0),
    };

    let omega_11 = omega_11_calc(T, e_k, mu, sigma_k);
    // println!("Omega_11 {}, ro {}, M1, {}, {}, {}, {}, {}", omega_11 , ro, M, sigma_k, e_k, rot_relax, mu);
    let f_t = 1.0
        + (PI.powf(1.5) / 2.0) * (e_k / T).sqrt()
        + (2.0 + PI.powi(2) / 4.0) * (e_k / T)
        + PI.powf(1.5) * (e_k / T).powf(3.2);
    let f_298 = 1.0
        + (PI.powf(1.5) / 2.0) * (e_k / 298.0).sqrt()
        + (2.0 + PI.powi(2) / 4.0) * (e_k / 298.0)
        + PI.powf(1.5) * (e_k / 298.0).powf(3.2);
    let z_rot = rot_relax * f_t / f_298;
    // self - diffusion coefficient
    let D_kk = 1e+20 * (3.0 / 16.0) * T.powf(1.5) * K_B * (2.0 * R / (M * PI)).sqrt()
        / (P * sigma_k.powi(2) * omega_11);
    let eta = visc(M, T, e_k, sigma_k, mu);

    let a = (5.0 / 2.0) - ro * D_kk / eta;
    let b = z_rot + (2.0 / PI) * ((5.0 / 3.0) * (C_rot / R) + ro * D_kk / eta);
    let f_trans = (5.0 / 2.0) * (1.0 - (2.0 / PI) * (C_rot / C_trans) * (a / b));
    let f_rot = ro * D_kk * (1.0 + (2.0 / PI) * (a / b)) / eta;
    let f_vib = ro * D_kk / eta;

    if form == 1.0 || form == 2.0 {
        um * (eta / M) * (f_trans * C_trans + f_rot * C_rot + f_vib * C_vib)
    } else {
        um * (eta / M) * ((5.0 / 2.0) * (3.0 / 2.0) * R)
    }
}

fn omega_11_calc_sym(e_k: f64, mu: f64, sigma_k: f64) -> Expr {
    let a = Expr::Const(1.06036);
    let b = Expr::Const(0.15610);
    let c = Expr::Const(0.19300);
    let d = Expr::Const(0.47635);
    let e = Expr::Const(1.03587);
    let f = Expr::Const(1.52996);
    let g = Expr::Const(1.76474);
    let h = Expr::Const(3.89411);
    let T = Expr::Var("T".to_owned());
    let T1 = T / Expr::Const(e_k);
    let delta = Expr::Const(1e-19 * mu.powi(2) / (2.0 * e_k * K_B * sigma_k.powi(3)));
    a / (T1.clone()).pow(b)
        + c / (d * T1.clone()).exp()
        + e / (f * T1.clone()).exp()
        + g / (h * T1.clone()).exp()
        + Expr::Const(0.19) * delta.pow(Expr::Const(2.0)) / T1
}
fn omega_22_calc_sym(e_k: f64, mu: f64, sigma_k: f64) -> Expr {
    let a1 = Expr::Const(1.16145);
    let b1 = Expr::Const(0.14874);
    let c1 = Expr::Const(0.52487);
    let d1 = Expr::Const(0.77320);
    let e1 = Expr::Const(2.16178);
    let f1 = Expr::Const(2.43787);
    let T = Expr::Var("T".to_owned());
    let T1 = T / Expr::Const(e_k);
    let delta = Expr::Const(1e-19 * mu.powi(2) / (2.0 * e_k * K_B * sigma_k.powi(3)));

    let Col = a1 / (T1.clone()).pow(b1)
        + c1 / (d1 * T1.clone()).exp()
        + e1 / (f1 * T1.clone()).exp()
        + Expr::Const(0.2) * delta.pow(Expr::Const(2.0)) / T1;
    //  println!("T1: {}, Col: {}", T1, Col);
    Col
}

fn visc_sym(M: f64, e_k: f64, sigma_k: f64, mu: f64) -> Expr {
    let T = Expr::Var("T".to_owned());
    // N*s/m2 to transfer to mkPoise *10^7
    let M = Expr::Const(M * 1000.0); // from kg/mol to g/mol
    let omega_22 = omega_22_calc_sym(e_k, mu, sigma_k);
    Expr::Const(1e-7 * 26.69) * (M * T).pow(Expr::Const(0.5))
        / (Expr::Const(sigma_k.powi(2)) * omega_22)
}

fn calculate_Lambda_sym(p: TransportInput, um: f64, M: f64, P: f64, C: Expr, ro: Expr) -> Expr {
    let (form, sigma_k, e_k, rot_relax, mu) = (p.Form, p.diam, p.well_depth, p.rot_relax, p.dipole);

    let C = C - R_sym; // from heat capacity at constant pressure to heat capacity at constant volume
    //   println!("form: {}, sigma_k: {}, e_k: {}, rot_relax: {}, mu: {}", form, sigma_k, e_k, rot_relax, mu);
    let (C_trans, C_rot, C_vib) = match form as i32 {
        1 => (
            Expr::Const(3.0 / 2.0) * R_sym,
            R_sym,
            C - Expr::Const(5.0 / 2.0) * R_sym,
        ),
        2 => (
            Expr::Const(3.0 / 2.0) * R_sym,
            Expr::Const(3.0 / 2.0) * R_sym,
            C - Expr::Const(3.0) * R_sym,
        ),
        _ => (
            Expr::Const(3.0 / 2.0) * R_sym,
            Expr::Const(0.0),
            Expr::Const(0.0),
        ),
    };
    let T = Expr::Var("T".to_owned());
    let omega_11 = omega_11_calc_sym(e_k, mu, sigma_k);
    // println!("Omega_11 {}, ro {}, M1, {}, {}, {}, {}, {}", omega_11 , ro, M, sigma_k, e_k, rot_relax, mu);
    let T1 = Expr::Const(e_k) / T.clone();
    let f_t = Expr::Const(1.0)
        + Expr::Const(PI.powf(1.5) / 2.0) * (T1.clone()).pow(Expr::Const(0.5))
        + Expr::Const(2.0 + PI.powi(2) / 4.0) * T1.clone()
        + Expr::Const(PI.powf(1.5)) * T1.clone().pow(Expr::Const(3.2));
    let f_298 = Expr::Const(1.0)
        + Expr::Const(
            (PI.powf(1.5) / 2.0) * (e_k / 298.0).sqrt()
                + (2.0 + PI.powi(2) / 4.0) * (e_k / 298.0)
                + PI.powf(1.5) * (e_k / 298.0).powf(3.2),
        );
    let z_rot = Expr::Const(rot_relax) * f_t / f_298;
    // self - diffusion coefficient
    let D_kk = Expr::Const(1e+20 * (3.0 / 16.0))
        * T.pow(Expr::Const(1.5))
        * Expr::Const(K_B * (2.0 * R / (M * PI)).sqrt())
        / (Expr::Const(P * sigma_k.powi(2)) * omega_11);
    let eta = visc_sym(M, e_k, sigma_k, mu);

    let a = Expr::Const(5.0 / 2.0) - ro.clone() * D_kk.clone() / eta.clone();
    let b = z_rot
        + Expr::Const(2.0 / PI)
            * (Expr::Const((5.0 / 3.0) / R) * C_rot.clone()
                + ro.clone() * D_kk.clone() / eta.clone());
    let f_trans = Expr::Const(5.0 / 2.0)
        * (Expr::Const(1.0)
            - Expr::Const(2.0 / PI) * (C_rot.clone() / C_trans.clone()) * (a.clone() / b.clone()));
    let f_rot = ro.clone() * D_kk.clone() * (Expr::Const(1.0) + Expr::Const(2.0 / PI) * (a / b))
        / eta.clone();
    let f_vib = ro.clone() * D_kk.clone() / eta.clone();

    if form == 1.0 || form == 2.0 {
        Expr::Const(um)
            * (eta / Expr::Const(M))
            * (f_trans * C_trans + f_rot * C_rot + f_vib * C_vib)
    } else {
        Expr::Const(um) * (eta / Expr::Const(M)) * Expr::Const((5.0 / 2.0) * (3.0 / 2.0) * R)
    }
}
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct TransportInput {
    pub Altname: Option<String>,
    /// Form of the molecule
    ///   Atom=0 Linear=1 Nonlinear=2
    pub Form: f64,
    ///  Lennard-Jones collision diameter in angstroms
    pub diam: f64,
    ///  dipole moment in Debye. Default: 0.0
    pub dipole: f64,
    ///  Polarizability in A^3. Default: 0.0
    pub polar: f64,
    ///  Number of rotational relaxation collisions at 298 K.  Dimensionless. Default: 0.0
    pub rot_relax: f64,
    /// Lennard-Jones well depth in Kelvin
    pub well_depth: f64,
}
pub struct TransportData {
    pub input: TransportInput,
    pub M: f64,
    pub P: f64,
    pub ro: Option<f64>,
    pub ro_sym: Option<Expr>,
    /// molar mass units choice
    pub M_unit: Option<String>,
    M_unit_multiplier: f64,
    /// pressure units choice
    pub P_unit: Option<String>,
    P_unit_multiplier: f64,
    /// density units choice
    pub ro_unit: Option<String>,
    ro_unit_multiplier: f64,
    /// thermal conductivity units choice
    pub L_unit: Option<String>,
    L_unit_multiplier: f64,
    /// viscosity units choice
    pub V_unit: Option<String>,
    V_unit_multiplier: f64,

    pub Lambda: f64,
    pub V: f64,

    pub Lambda_fun: Box<dyn Fn(f64) -> f64 + 'static>,
    pub V_fun: Box<dyn Fn(f64) -> f64>,

    pub Lambda_sym: Option<Expr>,
    pub V_sym: Option<Expr>,
}

impl TransportData {
    pub fn new() -> Self {
        let input = TransportInput {
            Altname: None,
            Form: 0.0,
            diam: 0.0,
            dipole: 0.0,
            polar: 0.0,
            rot_relax: 0.0,
            well_depth: 0.0,
        };
        Self {
            input: input,
            M: 0.0,
            P: 0.0,
            ro: None,
            ro_sym: None,
            M_unit: None,
            M_unit_multiplier: 1.0,
            P_unit: None,
            P_unit_multiplier: 1.0,
            ro_unit: None,
            ro_unit_multiplier: 1.0,
            L_unit: None,
            L_unit_multiplier: 1.0,
            V_unit: None,
            V_unit_multiplier: 1.0,

            Lambda: 0.0,
            V: 0.0,
            Lambda_fun: Box::new(|x| x),
            V_fun: Box::new(|x| x),
            Lambda_sym: None,
            V_sym: None,
        }
    } // end of new
    pub fn set_input(&mut self, input: TransportInput) {
        self.input = input.clone();
    }

    pub fn set_lambda_unit(&mut self, unit: Option<String>) -> Result<(), TransportError> {
        if let Some(unit) = unit {
            self.L_unit = Some(unit.clone());
            match unit.as_str() {
                "W/m/K" => self.L_unit_multiplier = 1.0,
                "mW/m/K" => self.L_unit_multiplier = 1E3,
                "mkW/m/K" => self.L_unit_multiplier = 1E6,
                "mkW/sm/K" => self.L_unit_multiplier = 1E+4,
                _ => return Err(TransportError::InvalidUnit(unit)),
            }
        }
        Ok(())
    }

    pub fn set_V_unit(&mut self, unit: Option<String>) -> Result<(), TransportError> {
        if let Some(unit) = unit {
            self.V_unit = Some(unit.clone());
            match unit.as_str() {
                "kg/m/s" | "Pa*s" => self.V_unit_multiplier = 1.0,
                "mkPa*s" => self.V_unit_multiplier = 1E6,
                _ => return Err(TransportError::InvalidUnit(unit)),
            }
        }
        Ok(())
    }

    pub fn set_M_unit(&mut self, unit: Option<String>) -> Result<(), TransportError> {
        if let Some(unit) = unit {
            self.M_unit = Some(unit.clone());
            match unit.as_str() {
                "kg/mol" => self.M_unit_multiplier = 1.0,
                "g/mol" => self.M_unit_multiplier = 1E-3,
                _ => return Err(TransportError::InvalidUnit(unit)),
            }
        }
        Ok(())
    }

    pub fn set_ro_unit(&mut self, unit: Option<String>) -> Result<(), TransportError> {
        if let Some(unit) = unit {
            self.ro_unit = Some(unit.clone());
            match unit.as_str() {
                "kg/m3" => self.ro_unit_multiplier = 1.0,
                "g/cm3" => self.ro_unit_multiplier = 1E3,
                _ => return Err(TransportError::InvalidUnit(unit)),
            }
        }
        Ok(())
    }

    pub fn set_P_unit(&mut self, unit: Option<String>) -> Result<(), TransportError> {
        if let Some(unit) = unit {
            self.P_unit = Some(unit.clone());
            match unit.as_str() {
                "Pa" => self.P_unit_multiplier = 1.0,
                "atm" => self.P_unit_multiplier = 101325.0,
                "bar" => self.P_unit_multiplier = 101325.0 / 1000.0,
                _ => return Err(TransportError::InvalidUnit(unit)),
            }
        }
        Ok(())
    }

    pub fn from_serde(
        &mut self,
        data: serde_json::Value,
    ) -> Result<(), super::transport_api::TransportError> {
        self.input = serde_json::from_value(data)
            .map_err(|e| super::transport_api::TransportError::SerdeError(e))?;
        Ok(())
    }
    /*
    e_k 'well_depth'
    sigma_k 'diam'
    mu - 'dipole'
    */
    pub fn calculate_Visc(&mut self, T: f64) -> Result<f64, TransportError> {
        if T <= 0.0 {
            return Err(TransportError::InvalidTemperature(T));
        }
        if self.M <= 0.0 {
            return Err(TransportError::InvalidMolarMass(self.M));
        }

        let p = self.input.clone();
        let M = self.M * self.M_unit_multiplier;
        let (_, sigma_k, e_k, _, mu) = (p.Form, p.diam, p.well_depth, p.rot_relax, p.dipole);
        let eta = visc(M, T, e_k, sigma_k, mu) * self.V_unit_multiplier;
        self.V = eta;
        Ok(eta)
    }
    pub fn ideal_gas(&mut self, T: f64) {
        let ro = (self.P * self.P_unit_multiplier) * (self.M * self.M_unit_multiplier) / (R * T);
        self.ro = Some(ro);
    }

    pub fn ideal_gas_sym(&mut self) {
        let P = self.P * self.P_unit_multiplier; // from other unit to Pa = N/m2
        let M = self.M * self.M_unit_multiplier; // from other unit to kg/mol
        let um = self.ro_unit_multiplier;
        let ro_sym = Expr::Const(P * M * um / R) / Expr::Var("T".to_owned());
        self.ro_unit = Some("kg/m3".to_string());
        self.ro_unit_multiplier = 1.0;
        self.ro_sym = Some(ro_sym);
    }

    pub fn calculate_Lambda(
        &mut self,
        C: f64,
        ro: Option<f64>,
        T: f64,
    ) -> Result<f64, TransportError> {
        if T <= 0.0 {
            return Err(TransportError::InvalidTemperature(T));
        }
        if self.M <= 0.0 {
            return Err(TransportError::InvalidMolarMass(self.M));
        }
        if self.P <= 0.0 {
            return Err(TransportError::InvalidPressure(self.P));
        }
        let ro = if let Some(ro) = ro {
            ro
        } else {
            self.ideal_gas(T);
            self.ro.unwrap()
        };
        if ro <= 0.0 {
            return Err(TransportError::InvalidDensity(ro));
        }

        let p = self.input.clone();
        let um = self.L_unit_multiplier.clone();
        let M = self.M * self.M_unit_multiplier; // from other unit to kg/mol
        let P = self.P * self.P_unit_multiplier; // from other unit to Pa = N/m2
        let ro = ro * self.ro_unit_multiplier; // from other unit to kg/m3
        let Lambda = calculate_Lambda_(p, um, M, P, C, ro, T);
        self.Lambda = Lambda;
        Ok(Lambda)
    }
    pub fn create_closure_Lambda(&mut self, C: f64, ro: Option<f64>) -> Result<(), TransportError> {
        let um = self.ro_unit_multiplier;
        let M = self.M * self.M_unit_multiplier; // from other unit to kg/mol
        let P = self.P * self.P_unit_multiplier; // from other unit to Pa = N/m2
        if self.M <= 0.0 {
            return Err(TransportError::InvalidMolarMass(self.M));
        }
        if self.P <= 0.0 {
            return Err(TransportError::InvalidPressure(self.P));
        }
        let ro_fn: Box<dyn Fn(f64) -> f64> = if let Some(ro) = ro {
            Box::new(move |_| ro * um)
        } else {
            Box::new(move |T| um * P * M / (R * T))
        };
        if ro_fn(300.0) <= 0.0 {
            return Err(TransportError::InvalidDensity(ro_fn(300.0)));
        }

        let p = self.input.clone();
        let um = self.L_unit_multiplier.clone();
        let Lambda = move |t: f64| calculate_Lambda_(p.clone(), um, M, P, C, ro_fn(t), t);
        self.Lambda_fun = Box::new(Lambda);
        Ok(())
    }
    pub fn create_closure_visc(&mut self) -> Result<(), TransportError> {
        if self.M <= 0.0 {
            return Err(TransportError::InvalidMolarMass(self.M));
        }

        let p = self.input.clone();
        let M = self.M * self.M_unit_multiplier;
        let vu = self.V_unit_multiplier;
        let (_, sigma_k, e_k, _, mu) = (p.Form, p.diam, p.well_depth, p.rot_relax, p.dipole);
        let visc = move |t: f64| visc(M, t, e_k, sigma_k, mu) * vu;
        self.V_fun = Box::new(visc);
        Ok(())
    }
    pub fn calculate_Lambda_sym(&mut self, C: Expr, ro: Expr) -> Result<(), TransportError> {
        if self.M <= 0.0 {
            return Err(TransportError::InvalidMolarMass(self.M));
        }
        if self.P <= 0.0 {
            return Err(TransportError::InvalidPressure(self.P));
        }

        let p = self.input.clone();
        let um = self.L_unit_multiplier.clone();
        let M = self.M * self.M_unit_multiplier; // from other unit to kg/mol
        let P = self.P * self.P_unit_multiplier; // from other unit to Pa = N/m2
        let ro = ro * Expr::Const(self.ro_unit_multiplier); // from other unit to kg/m3
        let Lambda = calculate_Lambda_sym(p, um, M, P, C, ro);
        self.Lambda_sym = Some(Lambda);
        Ok(())
    }
    pub fn Taylor_series_Lambda(
        &mut self,
        C: Expr,
        ro: Expr,
        T0: f64,
        n: usize,
    ) -> Result<Expr, TransportError> {
        if T0 <= 0.0 {
            return Err(TransportError::InvalidTemperature(T0));
        }
        self.calculate_Lambda_sym(C, ro)?;
        let Lambda_series = self
            .Lambda_sym
            .clone()
            .ok_or_else(|| {
                TransportError::CalculationError("Lambda_sym not calculated".to_string())
            })?
            .taylor_series1D_("T", T0, n);
        Ok(Lambda_series)
    }
}

impl fmt::Debug for TransportData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TransportData")
            .field("input", &self.input)
            .field("M", &self.M)
            .field("P", &self.P)
            .field("M_unit", &self.M_unit)
            .field("M_unit_multiplier", &self.M_unit_multiplier)
            .field("P_unit", &self.P_unit)
            .field("P_unit_multiplier", &self.P_unit_multiplier)
            .field("ro_unit", &self.ro_unit)
            .field("ro_unit_multiplier", &self.ro_unit_multiplier)
            .field("L_unit", &self.L_unit)
            .field("L_unit_multiplier", &self.L_unit_multiplier)
            .field("V_unit", &self.V_unit)
            .field("V_unit_multiplier", &self.V_unit_multiplier)
            .field("Lambda", &self.Lambda)
            .field("V", &self.V)
            .field("Lambda_fun", &"<closure>")
            .field("V_fun", &"<closure>")
            .field("Lambda_sym", &self.Lambda_sym)
            .field("V_sym", &self.V_sym)
            .finish_non_exhaustive() //  The finish_non_exhaustive() method is used to indicate that not all fields are being displayed.
    }
}
use super::transport_api::{LambdaUnit, ViscosityUnit};
use super::transport_api::{
    lambda_dimension, validate_molar_mass, validate_pressure, validate_temperature,
    viscosity_dimension,
};
impl super::transport_api::TransportCalculator for TransportData {
    fn extract_coefficients(
        &mut self,
        _t: f64,
    ) -> Result<(), super::transport_api::TransportError> {
        Ok(())
    }
    fn calculate_lambda(
        &mut self,
        C: Option<f64>,
        ro: Option<f64>,
        T: f64,
    ) -> Result<f64, super::transport_api::TransportError> {
        let C = C.unwrap();
        // let ro = ro.unwrap();
        validate_temperature(T)?;
        validate_molar_mass(self.M)?;

        let Lambda = self.calculate_Lambda(C, ro, T)?;
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
        C: Option<f64>,
        ro: Option<f64>,
    ) -> Result<Box<dyn Fn(f64) -> f64>, super::transport_api::TransportError> {
        validate_molar_mass(self.M)?;
        validate_pressure(self.P)?;
        let C = C.unwrap();
        let ro = ro.unwrap();
        let p = self.input.clone();
        let um = self.L_unit_multiplier.clone();
        let M = self.M * self.M_unit_multiplier; // from other unit to kg/mol
        let P = self.P * self.P_unit_multiplier; // from other unit to Pa = N/m2
        let ro = ro * self.ro_unit_multiplier; // from other unit to kg/m3
        let Lambda = move |t: f64| calculate_Lambda_(p.clone(), um, M, P, C, ro, t);
        self.Lambda_fun = Box::new(Lambda.clone());
        Ok(Box::new(Lambda))
    }

    fn create_viscosity_closure(
        &mut self,
    ) -> Result<Box<dyn Fn(f64) -> f64>, super::transport_api::TransportError> {
        validate_molar_mass(self.M)?;

        let p = self.input.clone();
        let M = self.M * self.M_unit_multiplier;
        let vu = self.V_unit_multiplier;
        let (_, sigma_k, e_k, _, mu) = (p.Form, p.diam, p.well_depth, p.rot_relax, p.dipole);
        let visc = move |t: f64| visc(M, t, e_k, sigma_k, mu) * vu;
        self.V_fun = Box::new(visc.clone());
        Ok(Box::new(visc))
    }

    fn create_symbolic_lambda(
        &mut self,
        C: Option<Expr>,
        ro: Option<Expr>,
    ) -> Result<(), super::transport_api::TransportError> {
        validate_molar_mass(self.M)?;
        validate_pressure(self.P)?;

        let p = self.input.clone();
        let um = self.L_unit_multiplier;
        let M = self.M * self.M_unit_multiplier;
        let P = self.P * self.P_unit_multiplier;
        let ro = ro.unwrap();
        let C = C.unwrap();

        let Lambda = calculate_Lambda_sym(p, um, M, P, C, ro);
        self.Lambda_sym = Some(Lambda);
        Ok(())
    }

    fn create_symbolic_viscosity(&mut self) -> Result<(), super::transport_api::TransportError> {
        validate_molar_mass(self.M)?;

        let p = self.input.clone();
        let M = self.M * self.M_unit_multiplier;
        let vu = self.V_unit_multiplier;
        let (_, sigma_k, e_k, _, mu) = (p.Form, p.diam, p.well_depth, p.rot_relax, p.dipole);
        let V_sym = visc_sym(M, e_k, sigma_k, mu) * Expr::Const(vu);
        self.V_sym = Some(V_sym);
        Ok(())
    }

    fn taylor_series_lambda(
        &mut self,
        C: Option<Expr>,
        ro: Option<Expr>,
        t0: f64,
        n: usize,
    ) -> Result<Expr, super::transport_api::TransportError> {
        validate_temperature(t0)?;
        self.create_symbolic_lambda(C, ro)?;
        let Lambda_series = self
            .Lambda_sym
            .clone()
            .ok_or_else(|| {
                super::transport_api::TransportError::CalculationError(
                    "Lambda_sym not calculated".to_string(),
                )
            })?
            .taylor_series1D_("T", t0, n);
        Ok(Lambda_series)
    }

    fn from_serde(
        &mut self,
        data: serde_json::Value,
    ) -> Result<(), super::transport_api::TransportError> {
        self.input = serde_json::from_value(data)
            .map_err(|e| super::transport_api::TransportError::SerdeError(e))?;
        Ok(())
    }
    fn set_M(
        &mut self,
        M: f64,
        M_unit: Option<String>,
    ) -> Result<(), super::transport_api::TransportError> {
        self.M = M;
        self.set_M_unit(M_unit)?;
        Ok(())
    }
    fn set_P(
        &mut self,
        P: f64,
        P_unit: Option<String>,
    ) -> Result<(), super::transport_api::TransportError> {
        self.P = P;
        self.set_P_unit(P_unit)?;
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

impl Clone for TransportData {
    fn clone(&self) -> Self {
        // Create a new instance with default closures
        let mut cloned = TransportData {
            input: self.input.clone(),
            M: self.M,
            P: self.P,
            ro: self.ro,
            ro_sym: self.ro_sym.clone(),
            M_unit: self.M_unit.clone(),
            M_unit_multiplier: self.M_unit_multiplier,
            P_unit: self.P_unit.clone(),
            P_unit_multiplier: self.P_unit_multiplier,
            ro_unit: self.ro_unit.clone(),
            ro_unit_multiplier: self.ro_unit_multiplier,
            L_unit: self.L_unit.clone(),
            L_unit_multiplier: self.L_unit_multiplier,
            V_unit: self.V_unit.clone(),
            V_unit_multiplier: self.V_unit_multiplier,
            Lambda: self.Lambda,
            V: self.V,
            Lambda_fun: Box::new(|x| x), // Default closure
            V_fun: Box::new(|x| x),      // Default closure
            Lambda_sym: self.Lambda_sym.clone(),
            V_sym: self.V_sym.clone(),
        };

        // Recreate the closures if necessary

        let _ = cloned.create_closure_Lambda(self.Lambda, self.ro);

        let _ = cloned.create_closure_visc();

        cloned
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::Thermodynamics::DBhandlers::NASAdata::NASAdata;
    use crate::Thermodynamics::thermo_lib_api::ThermoData;
    use approx::assert_relative_eq;

    #[test]
    fn lambda_calc() {
        /*     N2   */
        let mut tr = TransportData::new();
        let input = TransportInput {
            Altname: None,
            Form: 1.0,
            diam: 3.3,   // sigma_k
            dipole: 0.0, // mu
            polar: 0.0,
            rot_relax: 4.0,
            well_depth: 97.0, // e_k
        };
        tr.input = input;
        tr.set_M_unit(Some("g/mol".to_owned())).unwrap();
        tr.set_P_unit(Some("atm".to_owned())).unwrap();
        tr.set_V_unit(Some("mkPa*s".to_owned())).unwrap();
        tr.set_lambda_unit(Some("mW/m/K".to_owned())).unwrap();
        tr.M = 28.0;
        tr.P = 1.0;
        let T = 300.0;
        let C = 29.15;
        let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
        let L = tr.calculate_Lambda(C, Some(ro), T).unwrap();
        assert_relative_eq!(L, 26.0, epsilon = 5.0);
        println!("Lambda: {}", L);
    }

    #[test]
    fn vis_calc() {
        // Пример из книги Рид, Праусниц стр 353
        // at T = 220C => 293 K has viscosity 169 mcPoise
        let NH3visc = visc(17.0 / 1000.0, 493.0, 358.0, 3.15, 1.47);
        // at T = 298K => 298 K has viscosity 10.0 mPa*s
        //   let N2visc =  visc(28.0/1000.0,  473.0, 97.0, 3.3, 0.0);
        // to mkPa*s
        //  println!("NH3: {}, N2: {}", NH3visc*1e7, N2visc*1e7);
        //  assert_relative_eq!( N2visc*1e7, 251.0, epsilon = 5.0);
        assert_relative_eq!(NH3visc * 1e7, 169.0, epsilon = 5.0);
    }

    #[test]
    fn test_with_real_data() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        println!("CO_data: {}", CO_data);
        let mut tr = TransportData::new();
        tr.from_serde(CO_data.clone()).unwrap();
        tr.set_M_unit(Some("g/mol".to_owned())).unwrap();
        tr.set_P_unit(Some("atm".to_owned())).unwrap();
        tr.set_V_unit(Some("mkPa*s".to_owned())).unwrap();
        tr.set_lambda_unit(Some("mW/m/K".to_owned())).unwrap();
        let T = 473.15; // K 
        tr.P = 1.0;
        tr.M = 28.0; // g/mol
        let _ = tr.calculate_Visc(T);
        assert_relative_eq!(tr.V, 25.2, epsilon = 5.0);
        println!("Viscosity: {:?} mkPa*s", tr.V);

        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut NASA = NASAdata::new();
        let _ = NASA.from_serde(CO_data.clone());
        let _ = NASA.extract_coefficients(T);
        let _ = NASA.calculate_Cp_dH_dS(T);
        let Cp = NASA.Cp;
        println!("Cp: {}", Cp);
        let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
        let L = tr.calculate_Lambda(Cp, Some(ro), T).unwrap();
        println!("Lambda: {}", L);
    }

    #[test]
    fn test_with_real_data_sym() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        println!("CO_data: {}", CO_data);
        let mut tr = TransportData::new();
        tr.from_serde(CO_data.clone()).unwrap();
        tr.set_M_unit(Some("g/mol".to_owned())).unwrap();
        tr.set_P_unit(Some("atm".to_owned())).unwrap();
        tr.set_V_unit(Some("mkPa*s".to_owned())).unwrap();
        tr.set_lambda_unit(Some("mW/m/K".to_owned())).unwrap();
        let T = 473.15; // K 
        tr.P = 1.0;
        tr.M = 28.0; // g/mol
        let _ = tr.calculate_Visc(T);
        assert_relative_eq!(tr.V, 25.2, epsilon = 5.0);
        println!("Viscosity: {:?} mkPa*s", tr.V);

        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut NASA = NASAdata::new();
        let _ = NASA.from_serde(CO_data.clone());
        let _ = NASA.extract_coefficients(T);
        let _ = NASA.calculate_Cp_dH_dS(T);
        let Cp = NASA.Cp;
        println!("Cp: {}", Cp);
        let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
        let L = tr.calculate_Lambda(Cp, Some(ro), T).unwrap();
        tr.create_closure_Lambda(Cp, Some(ro)).unwrap();
        let Lambda_closure = &mut tr.Lambda_fun;
        let Lambda_from_closure = Lambda_closure(T);

        assert_eq!(Lambda_from_closure, L);

        tr.calculate_Lambda_sym(Expr::Const(Cp), Expr::Const(ro))
            .unwrap();
        let L_sym = tr.Lambda_sym.unwrap();
        let L_from_sym = L_sym.lambdify1D()(T);

        assert_relative_eq!(L_from_sym, L.clone(), epsilon = 1e-5);

        // let lambda = tr.calculate_Lambda(C, ro, T, M, P);
        //  println!("Lamnda {}", lambda);
        //  let visc = CEA.calculate_Visc(500.0);
        //  println!("Lambda, mW/m/K: {:?}, Visc: {:?}", lambda, visc);
    }
    #[test]
    fn test_with_real_data_sym_ideal_gas() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        println!("CO_data: {}", CO_data);
        let mut tr = TransportData::new();
        tr.from_serde(CO_data.clone()).unwrap();
        tr.set_M_unit(Some("g/mol".to_owned())).unwrap();
        tr.set_P_unit(Some("atm".to_owned())).unwrap();
        tr.set_V_unit(Some("mkPa*s".to_owned())).unwrap();
        tr.set_lambda_unit(Some("mW/m/K".to_owned())).unwrap();
        let T = 473.15; // K 
        tr.P = 1.0;
        tr.M = 28.0; // g/mol
        let _ = tr.calculate_Visc(T);
        assert_relative_eq!(tr.V, 25.2, epsilon = 5.0);
        println!("Viscosity: {:?} mkPa*s", tr.V);

        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut NASA = NASAdata::new();
        let _ = NASA.from_serde(CO_data.clone());
        let _ = NASA.extract_coefficients(T);
        let _ = NASA.calculate_Cp_dH_dS(T);
        let Cp = NASA.Cp;
        println!("Cp: {}", Cp);
        //  let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
        let L = tr.calculate_Lambda(Cp, None, T).unwrap();
        tr.create_closure_Lambda(Cp, None).unwrap();
        let Lambda_closure = &mut tr.Lambda_fun;
        let Lambda_from_closure = Lambda_closure(T);

        assert_eq!(Lambda_from_closure, L);
        tr.ideal_gas_sym();
        tr.calculate_Lambda_sym(Expr::Const(Cp), tr.ro_sym.clone().unwrap())
            .unwrap();
        let L_sym = tr.Lambda_sym.unwrap();
        let L_from_sym = L_sym.lambdify1D()(T);

        assert_relative_eq!(L_from_sym, L.clone(), epsilon = 1e-5);

        // let lambda = tr.calculate_Lambda(C, ro, T, M, P);
        //  println!("Lamnda {}", lambda);
        //  let visc = CEA.calculate_Visc(500.0);
        //  println!("Lambda, mW/m/K: {:?}, Visc: {:?}", lambda, visc);
    }
    #[test]
    fn test_with_real_data_taylor_sym() {
        use std::time::Instant;
        let now = Instant::now();
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        println!("CO_data: {}", CO_data);
        let mut tr = TransportData::new();
        tr.from_serde(CO_data.clone()).unwrap();
        tr.set_M_unit(Some("g/mol".to_owned())).unwrap();
        tr.set_P_unit(Some("atm".to_owned())).unwrap();
        tr.set_V_unit(Some("mkPa*s".to_owned())).unwrap();
        tr.set_lambda_unit(Some("mW/m/K".to_owned())).unwrap();
        let T = 473.15; // K 
        tr.P = 1.0;
        tr.M = 28.0; // g/mol
        let _ = tr.calculate_Visc(T);
        assert_relative_eq!(tr.V, 25.2, epsilon = 5.0);
        println!("Viscosity: {:?} mkPa*s", tr.V);

        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut NASA = NASAdata::new();
        let _ = NASA.from_serde(CO_data.clone());
        let _ = NASA.extract_coefficients(T);
        let _ = NASA.create_sym_Cp_dH_dS();
        let Cp_sym = NASA.clone().Cp_sym;
        NASA.calculate_Cp_dH_dS(T);
        let Cp = NASA.Cp;
        println!("Cp: {}", Cp);
        let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
        let L = tr.calculate_Lambda(Cp, Some(ro), T).unwrap();
        tr.create_closure_Lambda(Cp, Some(ro)).unwrap();
        let Lambda_closure = &mut tr.Lambda_fun;
        let Lambda_from_closure = Lambda_closure(T);

        assert_eq!(Lambda_from_closure, L);

        let taylor_series_Lambda = tr
            .Taylor_series_Lambda(Cp_sym, Expr::Const(ro), 400.0, 4)
            .unwrap();
        let taylor_series_Lambda = taylor_series_Lambda.lambdify1D()(T);
        let elapsed = now.elapsed().as_secs_f64();
        println!("Elapsed: {:.2?}", elapsed);
        assert_relative_eq!(taylor_series_Lambda, L, epsilon = 1.0);
    }

    #[test]
    fn test_invalid_units() {
        let mut tr = TransportData::new();
        assert!(tr.set_lambda_unit(Some("invalid_unit".to_owned())).is_err());
        assert!(tr.set_V_unit(Some("invalid_unit".to_owned())).is_err());
        assert!(tr.set_M_unit(Some("invalid_unit".to_owned())).is_err());
        assert!(tr.set_ro_unit(Some("invalid_unit".to_owned())).is_err());
        assert!(tr.set_P_unit(Some("invalid_unit".to_owned())).is_err());
    }

    #[test]
    fn test_invalid_parameters() {
        let mut tr = TransportData::new();
        let input = TransportInput {
            Altname: None,
            Form: 1.0,
            diam: 3.3,
            dipole: 0.0,
            polar: 0.0,
            rot_relax: 4.0,
            well_depth: 97.0,
        };
        tr.input = input;
        tr.M = -1.0;
        assert!(tr.calculate_Visc(300.0).is_err());
        assert!(tr.calculate_Lambda(29.15, Some(1.0), -300.0).is_err());
    }
}
