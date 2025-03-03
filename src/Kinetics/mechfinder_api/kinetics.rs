//v 0.1.1
#![allow(warnings)]
use RustedSciThe::symbolic::symbolic_engine::Expr;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::f64;

const R: f64 = 8.314;
const Rsym: Expr = Expr::Const(8.314);
// Different types of reactions proceeding here

// Struct for reaction type "elementary" with simplest form of
// kinetic constant - easy Arrhenius form  A*Temp.powf(*n)exp(-E/(Temp*R) )
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct ElementaryStruct {
    pub Arrenius: Vec<f64>,
}

impl ElementaryStruct {
    pub fn new(Arrenius: Vec<f64>) -> Self {
        Self { Arrenius }
    }
    pub fn K_const(&self, Temp: f64) -> f64 {
        let A: &f64 = &self.Arrenius[0];
        let n: &f64 = &self.Arrenius[1];
        let E: &f64 = &self.Arrenius[2];
        let K_const_: f64 = A * Temp.powf(*n) * f64::exp(-E / (Temp * R));
        return K_const_;
    }
    pub fn K_expr(&self, T: Expr) -> Expr {
        let A: f64 = self.Arrenius[0];
        let n: f64 = self.Arrenius[1];
        let E: f64 = self.Arrenius[2];
        let A = Expr::Const(A);
        let n = Expr::Const(n);
        let E = Expr::Const(E);
        let k0 = A * (T.clone()).pow(n);
        let k = k0 * (E / (Rsym * T)).exp();
        return k;
    }
}

// Struct for reaction type "falloff" - form of kinetic constant is more compicated
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct FalloffStruct {
    pub low_rate: Vec<f64>,
    pub high_rate: Vec<f64>,
    pub eff: HashMap<String, f64>,
    pub troe: Option<Vec<f64>>,
}
impl FalloffStruct {
    pub fn new(
        low_rate: Vec<f64>,
        high_rate: Vec<f64>,
        eff: HashMap<String, f64>,
        troe: Option<Vec<f64>>,
    ) -> Self {
        Self {
            low_rate,
            high_rate,
            eff,
            troe: None,
        }
    }
    pub fn K_const(&self, Temp: f64, Concentrations: HashMap<String, f64>) -> f64 {
        //
        let k_l = &self.low_rate[0];
        let b_l = &self.low_rate[1];
        let E_l = &self.low_rate[2];
        //
        let k_h = &self.high_rate[0];
        let b_h = &self.high_rate[1];
        let E_h = &self.high_rate[2];

        let K0 = k_l * f64::exp(-E_l / (R * Temp)) * Temp.powf(*b_l);
        let K_inf = k_h * f64::exp(-E_h / (R * Temp)) * Temp.powf(*b_h);
        let P_r = K0 / K_inf;
        // Calculate effective concentrations, e.g., by multiplying concentrations with coefficients from self.eff
        // Hashmap {substance: concentration}
        let mut Eff: f64 = 0.0;
        for (subs_name, C_i) in Concentrations.iter() {
            if self.eff.get(subs_name).is_some() {
                self.eff.get(subs_name).map(|&eff_i| Eff += eff_i * C_i);
            } else {
                Eff += C_i;
            };
            return Eff;
        }
        let k: f64 = {
            if let Some(troe) = &self.troe {
                let F_c = if troe.len() == 3 {
                    let A: f64 = troe[0];
                    let T_3 = troe[1];
                    let T_1 = troe[2];
                    let F_c = (1.0 - A) * f64::exp(-Temp / T_3) + A * f64::exp(-Temp / T_1);
                    F_c
                }
                //troe.len()==3
                else if troe.len() == 4 {
                    let A: f64 = troe[0];
                    let T_3 = troe[1];
                    let T_1 = troe[2];
                    let T_2 = troe[3];
                    let F_c = (1.0 - A) * f64::exp(-Temp / T_3)
                        + A * (f64::exp(-Temp / T_1) + f64::exp(-Temp / T_2));
                    F_c
                }
                //troe.len()==4
                else {
                    println!("Error in Troe parameters");
                    return 0.0;
                }; //troe.len()!=3,4
                let C: f64 = -0.4 - 0.67 * f64::log(F_c, 10.0);
                let N: f64 = 0.75 - 1.27 * f64::log(F_c, 10.0);
                let f_1: f64 = (f64::log(P_r, 10.0) + C) / (N - 0.14 * (f64::log(P_r, 10.0) + C));
                let F = 10.0_f64.powf(f64::log(F_c, 10.0) / (1.0 + f_1.powf(2.0)));
                let k = K_inf * (P_r / (1.0 + P_r)) * F;
                return k;
            } else {
                // there is no troe field
                let k = K_inf * (P_r / (1.0 + P_r));
                return k;
            }; //troe
        }; //k

        let K_const_: f64 = Eff * k;
        return K_const_;
    }
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct PressureStruct {
    pub Arrenius: HashMap<i32, Vec<f64>>,
}

impl PressureStruct {
    pub fn new(r#type: String, Arrenius: HashMap<i32, Vec<f64>>, eq: String) -> Self {
        Self { Arrenius }
    }
    pub fn K_const(&self, Temp: &f64, P: f64) -> f64 {
        let pressures: Vec<f64> = self.Arrenius.keys().map(|p| *p as f64).collect::<Vec<_>>();

        let location =
            pressures.binary_search_by(|v| v.partial_cmp(&Temp).expect("Couldn't compare values"));
        let k: f64 = {
            match location {
                //
                Ok(i) => println!("Found at {}", i),
                Err(i) => println!("Not found, could be inserted at {}", i),
            }
            return 0.0;
        };
    }
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct ThreeBodyStruct {
    pub Arrenius: Vec<f64>,
    pub eff: HashMap<String, f64>,
}

//
impl ThreeBodyStruct {
    pub fn new(Arrenius: Vec<f64>, eff: HashMap<String, f64>) -> Self {
        Self { Arrenius, eff }
    }
    pub fn K_const(&self, Temp: f64, Concentrations: HashMap<String, f64>) -> f64 {
        let A: &f64 = &self.Arrenius[0];
        let n: &f64 = &self.Arrenius[1];
        let E: &f64 = &self.Arrenius[2];
        // Hashmap {substance: concentration}
        let mut Eff: f64 = 0.0;
        for (subs_name, C_i) in Concentrations.iter() {
            if self.eff.get(subs_name).is_some() {
                self.eff.get(subs_name).map(|&eff_i| Eff += eff_i * C_i);
            } else {
                Eff += C_i;
            };
            return Eff;
        }
        let K_const_: f64 = Eff * A * Temp.powf(*n) * f64::exp(-E / (Temp * R));
        return K_const_;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_elementary_reaction() {
        let elementary_reaction = ElementaryStruct::new(vec![1.0, 2.0, 300.0]);

        let temp = 298.0; // 298K
        let expected_k_const = 1.0 * (298.0_f64).powf(2.0) * f64::exp(-300.0 / (298.0 * 8.314));
        assert!((elementary_reaction.K_const(temp) - expected_k_const).abs() < 1e-6);
    }

    #[test]
    fn test_falloff_reaction() {
        let falloff_reaction = FalloffStruct::new(
            vec![1.0, 2.0, 300.0],
            vec![10.0, 1.5, 400.0],
            HashMap::from([("H2".to_string(), 2.0), ("M".to_string(), 1.0)]),
            Some(vec![0.5, 300.0, 1000.0]),
        );

        let temp = 298.0; // 298K
        let concentrations = HashMap::from([("H2".to_string(), 2.0), ("M".to_string(), 1.0)]);
        let expected_k_const = falloff_reaction.K_const(temp, concentrations);
        assert!(expected_k_const > 0.0);
    }
    /*
    #[test]
    fn test_pressure_reaction() {
        let pressure_reaction = PressureStruct::new(
            "pressure".to_string(),
            HashMap::from([(1.0, vec![1.0, 2.0, 300.0])]),
            "H2+M<=>H+H+M".to_string(),
        );

        let temp = 298.0; // 298K
        let pressure = 1.0; // 1 atm
        // The K_const method for pressure reactions is not implemented in the provided code.
        // You would need to implement it based on the specific requirements of pressure-dependent reactions.
        // For now, we'll just assert that the method does not panic.
        assert!(!pressure_reaction.K_const(&temp, pressure).is_nan());
    }
    */
    #[test]
    fn test_threebody_reaction() {
        let threebody_reaction = ThreeBodyStruct::new(
            vec![1.0, 2.0, 300.0],
            HashMap::from([("H2".to_string(), 2.0), ("M".to_string(), 1.0)]),
        );

        let temp = 298.0; // 298K
        let concentrations = HashMap::from([("H2".to_string(), 2.0), ("M".to_string(), 1.0)]);
        let expected_k_const = threebody_reaction.K_const(temp, concentrations);
        assert!(expected_k_const > 0.0);
    }
}
