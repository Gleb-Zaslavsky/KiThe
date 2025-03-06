//v 0.1.1
#![allow(warnings)]
use RustedSciThe::symbolic::symbolic_engine::Expr;
use serde::{Deserialize, Serialize};
use serde::de::{self, Deserializer, MapAccess, Visitor};
use std::collections::{HashMap, HashSet};
use std::f64;
use std::fmt;
use std::str::FromStr;
const R: f64 = 8.314;
const Rsym: Expr = Expr::Const(8.314);
// Different types of reactions proceeding here
/////////////////////////ELEMENTARTY KINETICS///////////////////////////////////////////////////////////////
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
/////////////////////////FALLOFF KINETICS///////////////////////////////////////////////////////////////
// Struct for reaction type "falloff" - form of kinetic constant is more compicated
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct FalloffStruct {
    pub low_rate: Vec<f64>,
    pub high_rate: Vec<f64>,
    pub eff: Option<HashMap<String, f64>>,
    pub troe: Option<Vec<f64>>,
}
impl FalloffStruct {
    pub fn new(
        low_rate: Vec<f64>,
        high_rate: Vec<f64>,
        eff: Option<HashMap<String, f64>>,
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
        let mut Eff: f64 = 1.0;
        if let Some(eff) = self.eff.clone(){
        let mut Eff: f64 = 0.0;
        for (subs_name, C_i) in Concentrations.iter() {
            if eff.get(subs_name).is_some() {
                eff.get(subs_name).map(|&eff_i| Eff += eff_i * C_i);
            } else {
                Eff += C_i;
            };
            return Eff;
        }
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
    pub fn K_expr(&self, Temp: Expr, Concentrations: HashMap<String, Expr>) -> Expr {
        let k_l = Expr::Const(self.low_rate[0]);
        let b_l = Expr::Const( self.low_rate[1]);
        let E_l = Expr::Const( self.low_rate[2]);
        //
        let k_h = Expr::Const(self.high_rate[0]);
        let b_h = Expr::Const(self.high_rate[1]);
        let E_h = Expr::Const(self.high_rate[2]);
    
        let K0 = k_l * (-E_l / (Rsym * Temp.clone())).exp() * Temp.clone().pow(b_l);
        let K_inf = k_h * (-E_h / (Rsym * Temp.clone())).exp() * (Temp.clone()).pow(b_h);
        let P_r = K0 / K_inf.clone();
        // Calculate effective concentrations, e.g., by multiplying concentrations with coefficients from self.eff
        // Hashmap {substance: concentration}
        let mut Eff  = Expr::Const(1.0);
        if let Some(eff) = self.eff.clone(){
        let mut Eff: Expr = Expr::Const(0.0);
            for (subs_name, C_i) in Concentrations.iter() {
                if eff.get(subs_name).is_some() {
                    eff.get(subs_name).map(|&eff_i| Eff += (Expr::Const(eff_i).clone() * C_i.clone()));
                } else {
                    Eff += C_i.clone();
                };
                return Eff;
                }
        }
        let k= {
            if let Some(troe) = &self.troe {
                let F_c = if troe.len() == 3 {
                    let A  = Expr::Const(troe[0]);
                    let T_3 = Expr::Const(troe[1]);
                    let T_1 =Expr::Const( troe[2]);
                    let F_c = (Expr::Const(1.0) - A.clone()) * (-Temp.clone() / T_3).exp() + A * (-Temp.clone() / T_1).exp();
                    F_c.symplify()
                }
                //troe.len()==3
                else if troe.len() == 4 {
                    let A = Expr::Const(troe[0]);
                    let T_3 = Expr::Const(troe[1]);
                    let T_1 = Expr::Const(troe[2]);
                    let T_2 = Expr::Const(troe[3]);
                    let F_c = (Expr::Const(1.0) - A) * (-Temp.clone() / T_3).exp()
                        + Temp.clone() * ((-Temp.clone() / T_1).exp() + (-Temp.clone() / T_2).exp());
                    F_c.symplify()
                }
                //troe.len()==4
                else {
                    println!("Error in Troe parameters");
                    return Expr::Const(0.0);
                }; //troe.len()!=3,4
                
                let C =  Expr::Const(-0.4) -  Expr::Const(0.67) * (F_c.clone()).log10(); 
                let N = Expr::Const(0.75) - Expr::Const(1.27) * (F_c.clone()).log10();
                let f_1 = ((P_r.clone()).log10() + C.clone()) / (N -  Expr::Const(0.14) * ((P_r.clone()).log10() + C));
                let F = Expr::Const(10.0).pow( (F_c.clone()).log10() / (Expr::Const(1.0) + f_1.pow(Expr::Const(2.0)  )));
                let k = K_inf * (P_r.clone() / (Expr::Const(1.0) + P_r.clone())) * F;
                return k.symplify();
            } else {
                // there is no troe field
                let k = K_inf * (P_r.clone() / (Expr::Const(1.0) + P_r));
                return k.symplify();
            }; //troe
        }; //k

        let K_const_ = Eff * k;
        return K_const_;
    }
}
/////////////////////////////THREE-BODY KINETICS////////////////////////////////
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
    pub fn  K_expr(&self, Temp: Expr, Concentrations: HashMap<String, Expr>) -> Expr {
        let A: f64 = self.Arrenius[0];
        let n: f64 = self.Arrenius[1];
        let E: f64 = self.Arrenius[2];
        let A = Expr::Const(A);
        let n = Expr::Const(n);
        let E = Expr::Const(E);
        let k0 = A * (Temp.clone()).pow(n);
        let k = k0 * (E / (Rsym * Temp)).exp();
          // Calculate effective concentrations, e.g., by multiplying concentrations with coefficients from self.eff
        // Hashmap {substance: concentration}
        let eff = self.eff.clone();
        let mut Eff: Expr = Expr::Const(0.0);
            for (subs_name, C_i) in Concentrations.iter() {
                if eff.get(subs_name).is_some() {
                    eff.get(subs_name).map(|&eff_i| Eff += (Expr::Const(eff_i).clone() * C_i.clone()));
                } else {
                    Eff += C_i.clone();
                };
                return Eff;
                }
        
        let K_sym = (Eff*k).symplify();
        return K_sym;
    }
}

/////////////////////////PRESSURE DEPENDENT KINETICS///////////////////////////////////////////////////////////////
#[derive(Debug, Deserialize,  Serialize, Clone)]
pub struct PressureStruct {
    pub Arrenius: HashMap<String, Vec<f64>>,
}
impl PressureStruct {
    pub fn new( Arrenius: HashMap<String, Vec<f64>>, ) -> Self {
        Self { Arrenius }
    }

    pub fn K_const(&self, Temp: &f64, P: f64) -> f64 {
        let pressures_str: Vec<String> = self.Arrenius.keys().cloned().collect();
        let pressures = match convert_strings_to_f64(pressures_str.clone()) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("Error converting pressures to f64: {}", e);
                return 0.0; // or handle the error in a way that makes sense for your application
            }
        };
        let location:Result<usize, usize> = pressures.binary_search_by(|v| v.partial_cmp(&P).expect("Couldn't compare values"));
        
        match location {
            Ok(i) => {
                // Exact pressure found
                let arr_params = &self.Arrenius[pressures_str[i].as_str()];
                calculate_k(arr_params, Temp)
            },
            Err(i) if i > 0 && i < pressures.len() => {
                // Interpolate between two pressures
                let p_low = pressures[i-1];
                let p_high = pressures[i];
                let arr_low = &self.Arrenius[pressures_str[i].as_str()];
                let arr_high = &self.Arrenius[pressures_str[i].as_str()];
                let k_low = calculate_k(arr_low, Temp);
                let k_high = calculate_k(arr_high, Temp);
                interpolate(P, p_low, p_high, k_low, k_high)
            },
            _ => {
                // Pressure out of range, use closest available pressure
                let i = location.unwrap();
                let closest_p = if i == 0 { pressures[0] } else { pressures[pressures.len() - 1] };
                let arr_params = &self.Arrenius[closest_p.to_string().as_str()];
                calculate_k(arr_params, Temp)
            }
        }
    }


    pub fn K_expr(&self, Temp: Expr, P: f64) -> Expr {
        let pressures_str: Vec<String> = self.Arrenius.keys().cloned().collect();
        let pressures = match convert_strings_to_f64(pressures_str.clone()) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("Error converting pressures to f64: {}", e);
                return Expr::Const(0.0);
            }
        };
        let Ps = Expr::Const(P);
        let interpolate_expr = |p_low: f64, p_high: f64, k_low: Expr, k_high: Expr| -> Expr {
            let p_low = Expr::Const(p_low);
            let p_high = Expr::Const(p_high);
            k_low.clone() + (Ps.clone() - p_low.clone()) * (k_high - k_low) / (p_high - p_low)
        };

        let calculate_k_expr = |arr_params: &[f64], temp: Expr| -> Expr {
            let a = Expr::Const(arr_params[0]);
            let b = Expr::Const(arr_params[1]);
            let e = Expr::Const(arr_params[2]);
            a * temp.clone().pow(b) * (-e / (Rsym * temp)).exp()
        };

        let result:Result<usize, usize> = pressures.binary_search_by(|v| v.partial_cmp(&P).expect("Couldn't compare values"));

        match result {
            Ok(i) => {
                // Exact pressure found
                let arr_params = &self.Arrenius[&pressures_str[i]];
                calculate_k_expr(arr_params, Temp)
            },
            Err(i) if i > 0 && i < pressures.len() => {
                // Interpolate between two pressures
                let p_low = pressures[i-1];
                let p_high = pressures[i];
                let arr_low = &self.Arrenius[&pressures_str[i-1]];
                let arr_high = &self.Arrenius[&pressures_str[i]];
                let k_low = calculate_k_expr(arr_low, Temp.clone());
                let k_high = calculate_k_expr(arr_high, Temp);
                interpolate_expr(p_low, p_high, k_low, k_high)
            },
            _ => {
                // Pressure out of range, use closest available pressure
                let i = if result.is_ok() { result.unwrap() } else { result.unwrap_err() };
                let closest_p = if i == 0 { pressures[0] } else { pressures[pressures.len() - 1] };
                let arr_params = &self.Arrenius[&closest_p.to_string()];
                calculate_k_expr(arr_params, Temp)
            }
        }.symplify()
    }


}
fn calculate_k(arr_params: &[f64], temp: &f64) -> f64 {
    let a = arr_params[0];
    let b = arr_params[1];
    let e = arr_params[2];
    a * temp.powf(b) * (-e / (R * temp)).exp()
}

fn interpolate(x: f64, x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    y1 + (x - x1) * (y2 - y1) / (x2 - x1)
}

fn convert_strings_to_f64(strings: Vec<String>) -> Result<Vec<f64>, std::num::ParseFloatError> {
    strings.into_iter().map(|s| s.parse::<f64>()).collect()
}
//________________________________________________________________
use std::hash::{Hash, Hasher};
use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct F64Wrapper(pub f64);

impl PartialOrd for F64Wrapper {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Ord for F64Wrapper {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

impl Eq for F64Wrapper {}

impl Hash for F64Wrapper {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let bits = self.0.to_bits();
        bits.hash(state);
    }
}


fn deserialize_f64_wrapper_map<'de, D>(deserializer: D) -> Result<HashMap<F64Wrapper, Vec<f64>>, D::Error>
where
    D: Deserializer<'de>,
{
    struct F64WrapperMapVisitor;

    impl<'de> Visitor<'de> for F64WrapperMapVisitor {
        type Value = HashMap<F64Wrapper, Vec<f64>>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("a map with f64 keys and vec<f64> values")
        }

        fn visit_map<M>(self, mut access: M) -> Result<Self::Value, M::Error>
        where
            M: MapAccess<'de>,
        {
            let mut map = HashMap::with_capacity(access.size_hint().unwrap_or(0));

            while let Some((key, value)) = access.next_entry::<f64, Vec<f64>>()? {
                map.insert(F64Wrapper(key), value);
            }
            println!("map: {:#?} \n", map);
            Ok(map)
        }
    }

    deserializer.deserialize_map(F64WrapperMapVisitor)
}

fn deserialize_f64_wrapper_map2<'de, D>(deserializer: D) -> Result<HashMap<F64Wrapper, Vec<f64>>, D::Error>
where
    D: Deserializer<'de>,
{
    let map = HashMap::<String, Vec<f64>>::deserialize(deserializer)?;
    Ok(map.into_iter().map(|(k, v)| (F64Wrapper(k.parse::<f64>().unwrap()), v)).collect())
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
            Some(HashMap::from([("H2".to_string(), 2.0), ("M".to_string(), 1.0)])  ),
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
