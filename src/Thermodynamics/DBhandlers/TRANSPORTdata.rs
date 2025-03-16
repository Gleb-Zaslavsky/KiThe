use RustedSciThe::symbolic::symbolic_engine::Expr;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use std::fmt;

use std::f64::consts::PI;
const K_B: f64 = 1.38e-23;
const R: f64 = 8.314;
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
    let T1 = (Expr::Const(e_k) / T.clone());
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
    input: TransportInput,
    pub M: f64,
    pub P: f64,
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

    pub Lambda_fun: Box<dyn FnMut(f64) -> f64 + 'static>,
    pub V_fun: Box<dyn FnMut(f64) -> f64>,

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

    pub fn set_lambda_unit(&mut self, unit: Option<String>) {
        if let Some(unit) = unit {
            match unit.as_str() {
                "W/m/K" => self.L_unit_multiplier = 1.0,
                "mW/m/K" => self.L_unit_multiplier = 1E3,
                "mkW/m/K" => self.L_unit_multiplier = 1E6,

                "mkW/sm/K" => self.L_unit_multiplier = 1E+4,
                _ => panic!("Invalid unit: {}", unit),
            }
        }
    }

    pub fn set_V_unit(&mut self, unit: Option<String>) {
        if let Some(unit) = unit {
            match unit.as_str() {
                "kg/m/s" => self.V_unit_multiplier = 1.0,
                "Pa*s" => self.V_unit_multiplier = 1.0,
                "mkPa*s" => self.V_unit_multiplier = 1E6,

                _ => panic!("Invalid unit: {}", unit),
            }
        }
    }
    pub fn set_M_unit(&mut self, unit: Option<String>) {
        if let Some(unit) = unit {
            match unit.as_str() {
                "kg/mol" => self.M_unit_multiplier = 1.0,
                "g/mol" => self.M_unit_multiplier = 1E-3,
                _ => panic!("Invalid unit: {}", unit),
            }
        }
    }
    pub fn set_ro_unit(&mut self, unit: Option<String>) {
        if let Some(unit) = unit {
            match unit.as_str() {
                "kg/m3" => self.ro_unit_multiplier = 1.0,
                "g/cm3" => self.ro_unit_multiplier = 1E3,

                _ => panic!("Invalid unit: {}", unit),
            }
        }
    }
    pub fn set_P_unit(&mut self, unit: Option<String>) {
        if let Some(unit) = unit {
            match unit.as_str() {
                "Pa" => self.P_unit_multiplier = 1.0,
                "atm" => self.P_unit_multiplier = 101325.0,
                "bar" => self.P_unit_multiplier = 101325.0 / 1000.0,
                _ => panic!("Invalid unit: {}", unit),
            }
        }
    }
    pub fn from_serde(&mut self, serde: Value) {
        self.input = serde_json::from_value(serde).unwrap();
    }
    /*
    e_k 'well_depth'
    sigma_k 'diam'
    mu - 'dipole'
    */
    pub fn calculate_Visc(&mut self, T: f64) -> f64 {
        let p = self.input.clone();
        let M = self.M * self.M_unit_multiplier;
        let (_, sigma_k, e_k, _, mu) = (p.Form, p.diam, p.well_depth, p.rot_relax, p.dipole);
        let eta = visc(M, T, e_k, sigma_k, mu) * self.V_unit_multiplier;
        self.V = eta;
        eta
    }

    pub fn calculate_Lambda(&mut self, C: f64, ro: f64, T: f64) -> f64 {
        let p = self.input.clone();
        let um = self.L_unit_multiplier.clone();

        let M = self.M * self.M_unit_multiplier; // from other unit to kg/mol
        let P = self.P * self.P_unit_multiplier; // from other unit to Pa = N/m2
        let ro = ro * self.ro_unit_multiplier; // from other unit to kg/m3
        let Lambda = calculate_Lambda_(p, um, M, P, C, ro, T);
        self.Lambda = Lambda.clone();
        Lambda
    }
    pub fn create_closure_Lambda(&mut self, C: f64, ro: f64) {
        let p = self.input.clone();
        let um = self.L_unit_multiplier.clone();

        let M = self.M * self.M_unit_multiplier; // from other unit to kg/mol
        let P = self.P * self.P_unit_multiplier; // from other unit to Pa = N/m2
        let ro = ro * self.ro_unit_multiplier; // from other unit to kg/m3
        let Lambda = move |t: f64| calculate_Lambda_(p.clone(), um, M, P, C, ro, t);
        self.Lambda_fun = Box::new(Lambda);
        //    let Lambda = move |t: f64| self.calculate_Lambda_(C, ro, t);
        //  Box::new(Lambda);
    }
    pub fn create_closure_visc(&mut self) {
        let p = self.input.clone();
        let M = self.M * self.M_unit_multiplier;
        let vu = self.V_unit_multiplier;
        let (_, sigma_k, e_k, _, mu) = (p.Form, p.diam, p.well_depth, p.rot_relax, p.dipole);
        let visc = move |t: f64| visc(M, t, e_k, sigma_k, mu) * vu;
        self.V_fun = Box::new(visc);
    }
    pub fn calculate_Lambda_sym(&mut self, C: Expr, ro: Expr, T: Expr) {
        let p = self.input.clone();
        let um = self.L_unit_multiplier.clone();

        let M = self.M * self.M_unit_multiplier; // from other unit to kg/mol
        let P = self.P * self.P_unit_multiplier; // from other unit to Pa = N/m2
        let ro = ro * Expr::Const(self.ro_unit_multiplier); // from other unit to kg/m3
        let Lambda = calculate_Lambda_sym(p, um, M, P, C, ro);
        self.Lambda_sym = Some(Lambda.clone());
    }
}

impl fmt::Debug for TransportData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TransportData")
            .field("input", &self.input)
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
        tr.set_M_unit(Some("g/mol".to_owned()));
        tr.set_P_unit(Some("atm".to_owned()));
        tr.set_V_unit(Some("mkPa*s".to_owned()));
        tr.set_lambda_unit(Some("mW/m/K".to_owned()));
        //  tr.set_ro_unit(Some("g/cm3".to_owned())  );
        tr.M = 28.0;
        tr.P = 1.0;
        let T = 300.0;
        let C = 29.15;
        let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
        let L = tr.calculate_Lambda(C, ro, T);
        assert_relative_eq!(L, 26.0, epsilon = 5.0);
        println!("Lambda: {}", L);
        /*
        e_k 'well_depth'
        sigma_k 'diam'
        mu - 'dipole'
        */ //                                 sigma_k, e_k, rot_relax, mu
        //  calculate_Lambda(1, 29.15,  1.6, 300,  28, 1,  3.3, 97, 4, 0)
        // print('H2O',Lambda_calc(2, 35.3,  1.6, 400,  16, 1,   2.605, 572.400, 4, 1.844) )
        // print(Lambda_calc(  2,  2,   1, 600,   115, 10, 6.78, 521.551, 0.0, 0.0))
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
        tr.from_serde(CO_data.clone());
        tr.set_M_unit(Some("g/mol".to_owned()));
        tr.set_P_unit(Some("atm".to_owned()));
        tr.set_V_unit(Some("mkPa*s".to_owned()));
        tr.set_lambda_unit(Some("mW/m/K".to_owned()));
        let T = 473.15; // K 
        tr.P = 1.0;
        tr.M = 28.0; // g/mol
        tr.calculate_Visc(T);
        assert_relative_eq!(tr.V, 25.2, epsilon = 5.0);
        println!("Viscosity: {:?} mkPa*s", tr.V);
        let P = 1.0; // atm

        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut NASA = NASAdata::new();
        NASA.from_serde(CO_data.clone());
        NASA.extract_coefficients(T);
        NASA.calculate_Cp_dH_dS(T);
        let Cp = NASA.Cp;
        println!("Cp: {}", Cp,);
        let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
        let L = tr.calculate_Lambda(Cp, ro, T);
        println!("Lambda: {}", L);

        // let lambda = tr.calculate_Lambda(C, ro, T, M, P);
        //  println!("Lamnda {}", lambda);
        //  let visc = CEA.calculate_Visc(500.0);
        //  println!("Lambda, mW/m/K: {:?}, Visc: {:?}", lambda, visc);
    }
}
