//v 0.1.1
use RustedSciThe::symbolic::symbolic_engine::Expr;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::f64;
const R: f64 = 8.314;
const R_SYM: Expr = Expr::Const(8.314);

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
    pub fn K_expr(&self) -> Expr {
        let T = Expr::Var("T".to_owned());
        let A: f64 = self.Arrenius[0];
        let n: f64 = self.Arrenius[1];
        let E: f64 = self.Arrenius[2];
        let A = Expr::Const(A);
        let n = Expr::Const(n);
        let E = Expr::Const(E);
        let k0 = A * (T.clone()).pow(n);
        let k = k0 * (-E / (R_SYM * T)).exp();
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
            troe,
        }
    }
    /// Evaluate the falloff rate coefficient at temperature in kelvin.
    ///
    /// Collider concentrations use the same concentration-unit convention as
    /// the source mechanism's Arrhenius parameters; no implicit conversion is
    /// performed. CHEMKIN-style mechanism data commonly use mol/cm^3. Third-body
    /// efficiencies are dimensionless.
    pub fn K_const(&self, Temp: f64, Concentrations: HashMap<String, f64>) -> f64 {
        let k0 = self.low_rate[0]
            * f64::exp(-self.low_rate[2] / (R * Temp))
            * Temp.powf(self.low_rate[1]);
        let k_inf = self.high_rate[0]
            * f64::exp(-self.high_rate[2] / (R * Temp))
            * Temp.powf(self.high_rate[1]);
        let effective_concentration = Concentrations
            .iter()
            .map(|(species, concentration)| {
                self.eff
                    .as_ref()
                    .and_then(|efficiencies| efficiencies.get(species))
                    .copied()
                    .unwrap_or(1.0)
                    * concentration
            })
            .sum::<f64>();

        if effective_concentration == 0.0 {
            return 0.0;
        }

        let P_r = k0 * effective_concentration / k_inf;
        let k = if let Some(troe) = &self.troe {
            let F_c = match troe.len() {
                3 => {
                    let A: f64 = troe[0];
                    let T_3 = troe[1];
                    let T_1 = troe[2];
                    (1.0 - A) * f64::exp(-Temp / T_3) + A * f64::exp(-Temp / T_1)
                }
                4 => {
                    let A: f64 = troe[0];
                    let T_3 = troe[1];
                    let T_1 = troe[2];
                    let T_2 = troe[3];
                    (1.0 - A) * f64::exp(-Temp / T_3)
                        + A * f64::exp(-Temp / T_1)
                        + f64::exp(-T_2 / Temp)
                }
                _ => return 0.0,
            }; //troe.len()!=3,4
            let C: f64 = -0.4 - 0.67 * f64::log(F_c, 10.0);
            let N: f64 = 0.75 - 1.27 * f64::log(F_c, 10.0);
            let f_1: f64 = (f64::log(P_r, 10.0) + C) / (N - 0.14 * (f64::log(P_r, 10.0) + C));
            let F = 10.0_f64.powf(f64::log(F_c, 10.0) / (1.0 + f_1.powf(2.0)));
            k_inf * (P_r / (1.0 + P_r)) * F
        } else {
            k_inf * (P_r / (1.0 + P_r))
        };
        k
    }
    pub fn K_expr(&self, Concentrations: HashMap<String, Expr>) -> Expr {
        let T = Expr::Var("T".to_owned());
        let k_l = Expr::Const(self.low_rate[0]);
        let b_l = Expr::Const(self.low_rate[1]);
        let E_l = Expr::Const(self.low_rate[2]);
        //
        let k_h = Expr::Const(self.high_rate[0]);
        let b_h = Expr::Const(self.high_rate[1]);
        let E_h = Expr::Const(self.high_rate[2]);

        let K0 = k_l * (-E_l / (R_SYM * T.clone())).exp() * T.clone().pow(b_l);
        let K_inf = k_h * (-E_h / (R_SYM * T.clone())).exp() * (T.clone()).pow(b_h);
        let mut effective_concentration = Expr::Const(0.0);
        for (species, concentration) in &Concentrations {
            let efficiency = self
                .eff
                .as_ref()
                .and_then(|efficiencies| efficiencies.get(species))
                .copied()
                .unwrap_or(1.0);
            effective_concentration += Expr::Const(efficiency) * concentration.clone();
        }
        let effective_concentration = effective_concentration.simplify();
        if effective_concentration.is_zero() {
            return Expr::Const(0.0);
        }
        let P_r = K0 * effective_concentration / K_inf.clone();
        let k = if let Some(troe) = &self.troe {
            let F_c = match troe.len() {
                3 => {
                    let A = Expr::Const(troe[0]);
                    let T_3 = Expr::Const(troe[1]);
                    let T_1 = Expr::Const(troe[2]);
                    ((Expr::Const(1.0) - A.clone()) * (-T.clone() / T_3).exp()
                        + A * (-T.clone() / T_1).exp())
                    .simplify()
                }
                4 => {
                    let A = Expr::Const(troe[0]);
                    let T_3 = Expr::Const(troe[1]);
                    let T_1 = Expr::Const(troe[2]);
                    let T_2 = Expr::Const(troe[3]);
                    ((Expr::Const(1.0) - A) * (-T.clone() / T_3).exp()
                        + Expr::Const(troe[0]) * (-T.clone() / T_1).exp()
                        + (-T_2 / T.clone()).exp())
                    .simplify()
                }
                _ => return Expr::Const(0.0),
            }; //troe.len()!=3,4

            let C = Expr::Const(-0.4) - Expr::Const(0.67) * (F_c.clone()).log10();
            let N = Expr::Const(0.75) - Expr::Const(1.27) * (F_c.clone()).log10();
            let f_1 = ((P_r.clone()).log10() + C.clone())
                / (N - Expr::Const(0.14) * ((P_r.clone()).log10() + C));
            let F = Expr::Const(10.0)
                .pow((F_c.clone()).log10() / (Expr::Const(1.0) + f_1.pow(Expr::Const(2.0))));
            (K_inf * (P_r.clone() / (Expr::Const(1.0) + P_r.clone())) * F).simplify()
        } else {
            (K_inf * (P_r.clone() / (Expr::Const(1.0) + P_r))).simplify()
        };
        k
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
            //  return Eff;
            // TODODO: handle case when substance is not in eff map
        }
        let K_const_: f64 = Eff * A * Temp.powf(*n) * f64::exp(-E / (Temp * R));
        return K_const_;
    }
    pub fn K_expr(&self, Concentrations: HashMap<String, Expr>) -> Expr {
        let T = Expr::Var("T".to_owned());
        let A: f64 = self.Arrenius[0];
        let n: f64 = self.Arrenius[1];
        let E: f64 = self.Arrenius[2];
        let A = Expr::Const(A);
        let n = Expr::Const(n);
        let E = Expr::Const(E);
        let k0 = A * (T.clone()).pow(n);
        let k = k0 * (-E / (R_SYM * T)).exp();
        // Calculate effective concentrations, e.g., by multiplying concentrations with coefficients from self.eff
        // Hashmap {substance: concentration}
        let eff = self.eff.clone();
        let mut Eff: Expr = Expr::Const(0.0);
        for (subs_name, C_i) in Concentrations.iter() {
            if let Some(&eff_i) = eff.get(subs_name) {
                Eff += Expr::Const(eff_i) * C_i.clone();
            } else {
                Eff += C_i.clone();
            };
            //   return Eff;
        }

        let K_sym = (Eff * k).simplify();
        return K_sym;
    }
}

/////////////////////////PRESSURE DEPENDENT KINETICS///////////////////////////////////////////////////////////////
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct PressureStruct {
    pub Arrenius: HashMap<String, Vec<f64>>,
}
impl PressureStruct {
    pub fn new(Arrenius: HashMap<String, Vec<f64>>) -> Self {
        Self { Arrenius }
    }

    fn sorted_pressure_points(&self) -> Vec<(f64, &Vec<f64>)> {
        let mut points: Vec<(f64, &Vec<f64>)> = self
            .Arrenius
            .iter()
            .filter_map(|(pressure, params)| pressure.parse::<f64>().ok().map(|p| (p, params)))
            .collect();
        points.sort_by(|a, b| a.0.total_cmp(&b.0));
        points
    }

    fn pressure_interval(pressures: &[f64], P: f64) -> (usize, usize) {
        match pressures.binary_search_by(|v| v.total_cmp(&P)) {
            Ok(i) => (i, i),
            Err(0) => (0, 0),
            Err(i) if i >= pressures.len() => (pressures.len() - 1, pressures.len() - 1),
            Err(i) => (i - 1, i),
        }
    }

    pub fn K_const(&self, Temp: f64, P: f64) -> f64 {
        let points = self.sorted_pressure_points();
        if points.is_empty() {
            return 0.0;
        }
        let pressures: Vec<f64> = points.iter().map(|(pressure, _)| *pressure).collect();
        let (low_index, high_index) = Self::pressure_interval(&pressures, P);

        if low_index == high_index {
            let arr_params = points[low_index].1;
            return calculate_k(arr_params, &Temp);
        }

        let p_low = pressures[low_index];
        let p_high = pressures[high_index];
        let arr_low = points[low_index].1;
        let arr_high = points[high_index].1;
        let k_low = calculate_k(arr_low, &Temp);
        let k_high = calculate_k(arr_high, &Temp);
        interpolate(P, p_low, p_high, k_low, k_high)
    }

    pub fn K_expr(&self, P: f64) -> Expr {
        let T = Expr::Var("T".to_owned());
        let points = self.sorted_pressure_points();
        if points.is_empty() {
            return Expr::Const(0.0);
        }
        let pressures: Vec<f64> = points.iter().map(|(pressure, _)| *pressure).collect();
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
            a * temp.clone().pow(b) * (-e / (R_SYM * temp)).exp()
        };

        let (low_index, high_index) = Self::pressure_interval(&pressures, P);

        if low_index == high_index {
            let arr_params = points[low_index].1;
            return calculate_k_expr(arr_params, T).simplify();
        }

        let p_low = pressures[low_index];
        let p_high = pressures[high_index];
        let arr_low = points[low_index].1;
        let arr_high = points[high_index].1;
        let k_low = calculate_k_expr(arr_low, T.clone());
        let k_high = calculate_k_expr(arr_high, T);
        interpolate_expr(p_low, p_high, k_low, k_high).simplify()
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

use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::mechfinder_api::{
        KineticsError, ReactionData, ReactionType, parse_kinetic_data_vec,
    };
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use approx::assert_relative_eq;

    fn eval_expr(expr: &Expr) -> f64 {
        match expr {
            Expr::Const(value) => *value,
            Expr::Add(lhs, rhs) => eval_expr(lhs) + eval_expr(rhs),
            Expr::Sub(lhs, rhs) => eval_expr(lhs) - eval_expr(rhs),
            Expr::Mul(lhs, rhs) => eval_expr(lhs) * eval_expr(rhs),
            Expr::Div(lhs, rhs) => eval_expr(lhs) / eval_expr(rhs),
            Expr::Pow(base, exp) => eval_expr(base).powf(eval_expr(exp)),
            Expr::Exp(inner) => eval_expr(inner).exp(),
            Expr::Ln(inner) => eval_expr(inner).ln(),
            other => panic!("unexpected symbolic node: {:?}", other),
        }
    }

    #[test]
    fn test_elementary_reaction() {
        let elementary_reaction = ElementaryStruct::new(vec![1.0, 2.0, 300.0]);

        let temp = 298.0; // 298K
        let expected_k_const = 1.0 * (298.0_f64).powf(2.0) * f64::exp(-300.0 / (298.0 * 8.314));
        assert!((elementary_reaction.K_const(temp) - expected_k_const).abs() < 1e-6);
    }

    #[test]
    fn test_falloff_reaction() {
        let troe = Some(vec![0.5, 300.0, 1000.0]);
        let falloff_reaction = FalloffStruct::new(
            vec![1.0, 2.0, 300.0],
            vec![10.0, 1.5, 400.0],
            Some(HashMap::from([
                ("H2".to_string(), 2.0),
                ("M".to_string(), 1.0),
            ])),
            troe.clone(),
        );

        let temp = 298.0; // 298K
        let concentrations = HashMap::from([("H2".to_string(), 2.0), ("M".to_string(), 1.0)]);
        let expected_k_const = falloff_reaction.K_const(temp, concentrations);
        assert!(expected_k_const > 0.0);
        assert_eq!(falloff_reaction.troe, troe);
    }

    #[test]
    fn falloff_uses_effective_collider_concentration_and_has_physical_limits() {
        let falloff = FalloffStruct::new(
            vec![1.0e12, 0.0, 5000.0],
            vec![1.0e10, 0.0, 10000.0],
            Some(HashMap::from([("H2".to_string(), 2.0)])),
            None,
        );

        let no_colliders = falloff.K_const(1000.0, HashMap::from([("H2".to_string(), 0.0)]));
        assert_eq!(no_colliders, 0.0);

        // [H2]=0.5 with efficiency 2 and [Ar]=1 with default efficiency 1 gives [M]eff=2.
        let concentrations = HashMap::from([("H2".to_string(), 0.5), ("Ar".to_string(), 1.0)]);
        let at_reference_concentration = falloff.K_const(1000.0, concentrations);
        assert!((at_reference_concentration - 2.9953411585551653e9).abs() < 1.0e-3);

        let without_efficiency_map = FalloffStruct::new(
            vec![1.0e12, 0.0, 5000.0],
            vec![1.0e10, 0.0, 10000.0],
            None,
            None,
        );
        let with_empty_efficiency_map = FalloffStruct::new(
            vec![1.0e12, 0.0, 5000.0],
            vec![1.0e10, 0.0, 10000.0],
            Some(HashMap::new()),
            None,
        );
        let default_efficiency_concentrations =
            HashMap::from([("H2".to_string(), 0.5), ("Ar".to_string(), 1.0)]);
        let no_map_k =
            without_efficiency_map.K_const(1000.0, default_efficiency_concentrations.clone());
        let empty_map_k =
            with_empty_efficiency_map.K_const(1000.0, default_efficiency_concentrations);
        assert!((no_map_k - 2.9926151416447945e9).abs() < 1.0e-3);
        assert_eq!(empty_map_k, no_map_k);

        let high_pressure_limit =
            falloff.K_const(1000.0, HashMap::from([("H2".to_string(), 1.0e8)]));
        let k_inf = 1.0e10 * (-10000.0_f64 / (R * 1000.0)).exp();
        assert!((high_pressure_limit / k_inf - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn troe_three_and_four_parameter_forms_match_independent_references_and_symbolic_route() {
        // Frozen evaluations of Cantera's published Troe equations
        // (rate-constants.html#the-troe-falloff-function; R=8.314) for these
        // coefficients, temperatures, and effective collider concentrations.
        let expected = [
            (
                vec![0.5, 1000.0, 5000.0],
                [
                    (700.0, 0.2, 1.60637874974864e9),
                    (1000.0, 2.0, 2.8624702685039926e9),
                    (1500.0, 20.0, 4.346615675659192e9),
                ],
            ),
            (
                vec![0.5, 1000.0, 5000.0, 10000.0],
                [
                    (700.0, 0.2, 1.606379186053633e9),
                    (1000.0, 2.0, 2.8625000064569616e9),
                    (1500.0, 20.0, 4.347542424109653e9),
                ],
            ),
        ];

        for (troe, reference_cases) in expected {
            let falloff = FalloffStruct::new(
                vec![1.0e12, 0.0, 5000.0],
                vec![1.0e10, 0.0, 10000.0],
                Some(HashMap::from([("H2".to_string(), 2.0)])),
                Some(troe),
            );

            for (temperature, effective_concentration, reference_k) in reference_cases {
                let concentrations =
                    HashMap::from([("H2".to_string(), effective_concentration / 2.0)]);
                let symbolic_concentrations: HashMap<String, Expr> = concentrations
                    .iter()
                    .map(|(species, concentration)| (species.clone(), Expr::Const(*concentration)))
                    .collect();

                let numeric_k = falloff.K_const(temperature, concentrations);
                assert!(
                    (numeric_k - reference_k).abs() < 1.0e-3,
                    "numeric rate {numeric_k:e} differs from reference {reference_k:e}"
                );

                let symbolic_k = falloff
                    .K_expr(symbolic_concentrations)
                    .set_variable("T", temperature);
                let symbolic_value = eval_expr(&symbolic_k);
                assert!(
                    (symbolic_value - numeric_k).abs() < numeric_k.abs() * 1.0e-9,
                    "symbolic rate {symbolic_value:e} differs from numeric rate {numeric_k:e}"
                );
            }
        }
    }

    #[test]
    fn troe_symbolic_temperature_derivative_matches_independent_central_difference() {
        for troe in [
            vec![0.5, 1000.0, 5000.0],
            vec![0.5, 1000.0, 5000.0, 10000.0],
        ] {
            let falloff = FalloffStruct::new(
                vec![1.0e12, 0.3, 5000.0],
                vec![1.0e10, -0.2, 10000.0],
                Some(HashMap::from([("H2".to_string(), 2.0)])),
                Some(troe),
            );
            let concentrations = HashMap::from([("H2".to_string(), 1.0)]);
            let symbolic_concentrations = HashMap::from([("H2".to_string(), Expr::Const(1.0))]);
            let rate_expression = falloff.K_expr(symbolic_concentrations);
            let derivative_expression = rate_expression.diff("T");

            for temperature in [700.0, 1000.0, 1500.0] {
                let step = temperature * 1.0e-5;
                let finite_difference = (falloff
                    .K_const(temperature + step, concentrations.clone())
                    - falloff.K_const(temperature - step, concentrations.clone()))
                    / (2.0 * step);
                let symbolic_derivative =
                    derivative_expression.eval_expression(&["T"], &[temperature]);
                let relative_error = (symbolic_derivative - finite_difference).abs()
                    / finite_difference.abs().max(1.0);
                assert!(
                    symbolic_derivative.is_finite() && relative_error < 2.0e-5,
                    "Troe {:?} at T={temperature}: symbolic dK/dT={symbolic_derivative:e}, finite difference={finite_difference:e}, relative error={relative_error:e}",
                    falloff.troe
                );
            }
        }
    }

    #[test]
    fn troe_symbolic_rate_is_zero_when_effective_collider_concentration_is_zero() {
        let falloff = FalloffStruct::new(
            vec![1.0e12, 0.0, 5000.0],
            vec![1.0e10, 0.0, 10000.0],
            None,
            Some(vec![0.5, 1000.0, 5000.0, 10000.0]),
        );
        let symbolic = falloff.K_expr(HashMap::from([("M".to_string(), Expr::Const(0.0))]));

        assert_eq!(symbolic, Expr::Const(0.0));
        assert_eq!(
            falloff.K_const(1000.0, HashMap::from([("M".to_string(), 0.0)])),
            0.0
        );
    }

    #[test]
    fn reaction_data_falloff_dispatch_preserves_numeric_and_symbolic_reference() {
        let reaction = crate::Kinetics::mechfinder_api::ReactionData::new_falloff(
            "A (+M) => B (+M)".to_string(),
            vec![1.0e12, 0.0, 5000.0],
            vec![1.0e10, 0.0, 10000.0],
            Some(HashMap::from([("H2".to_string(), 2.0)])),
            Some(vec![0.5, 1000.0, 5000.0, 10000.0]),
            None,
        );
        let concentrations = HashMap::from([("H2".to_string(), 1.0)]);
        let numeric_k = reaction
            .K_const(1000.0, None, Some(concentrations.clone()))
            .expect("falloff reaction accepts collider concentrations");
        assert!((numeric_k - 2.8625000064569616e9).abs() < 1.0e-3);

        let symbolic_concentrations = HashMap::from([("H2".to_string(), Expr::Const(1.0))]);
        let symbolic_k = reaction
            .K_sym(None, Some(symbolic_concentrations))
            .expect("falloff symbolic route accepts collider concentrations")
            .set_variable("T", 1000.0);
        assert!((eval_expr(&symbolic_k) - numeric_k).abs() < numeric_k * 1.0e-9);
    }

    #[test]
    fn public_numeric_dispatch_returns_typed_errors_for_invalid_payloads_and_domains() {
        let valid = ReactionData::new_elementary("A -> B".into(), vec![1.0, 0.0, 0.0], None);
        assert!(valid.K_const(0.0, None, None).is_err());
        assert!(valid.K_const(f64::NAN, None, None).is_err());

        let short_arrhenius = ReactionData::new_elementary("A -> B".into(), vec![1.0, 0.0], None);
        assert!(short_arrhenius.K_const(300.0, None, None).is_err());
        assert!(short_arrhenius.K_sym(None, None).is_err());
        let non_finite_parameters =
            ReactionData::new_elementary("A -> B".into(), vec![1.0, f64::NAN, 0.0], None);
        assert!(non_finite_parameters.K_const(300.0, None, None).is_err());
        let negative_pre_exponential =
            ReactionData::new_elementary("A -> B".into(), vec![-1.0, 0.0, 0.0], None);
        assert!(matches!(
            negative_pre_exponential.K_const(300.0, None, None),
            Err(KineticsError::InvalidReactionPayload { field, .. }) if field == "Arrenius"
        ));
        assert!(matches!(
            valid.K_const(-1.0, None, None),
            Err(KineticsError::InvalidKineticInput {
                field: "temperature",
                ..
            })
        ));
        assert!(matches!(
            valid.K_const(f64::NAN, None, None),
            Err(KineticsError::InvalidKineticInput {
                field: "temperature",
                ..
            })
        ));

        let falloff = ReactionData::new_falloff(
            "A (+M) -> B (+M)".into(),
            vec![1.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            None,
            None,
            None,
        );
        assert!(matches!(
            falloff.K_const(300.0, None, None),
            Err(KineticsError::MissingKineticInput {
                field: "concentrations",
                ..
            })
        ));
        assert!(matches!(
            falloff.K_const(300.0, None, Some(HashMap::from([("M".into(), -1.0)]))),
            Err(KineticsError::InvalidKineticInput {
                field: "concentrations",
                ..
            })
        ));
        let pressure = ReactionData::new_pressure(
            "A -> B".into(),
            HashMap::from([("1.0".into(), vec![1.0, 0.0, 0.0])]),
            None,
        );
        assert!(matches!(
            pressure.K_const(300.0, None, None),
            Err(KineticsError::MissingKineticInput {
                field: "pressure",
                ..
            })
        ));
        assert!(matches!(
            pressure.K_const(300.0, Some(0.0), None),
            Err(KineticsError::InvalidKineticInput {
                field: "pressure",
                ..
            })
        ));
        let negative_temperature_exponent_and_activation_energy =
            ReactionData::new_elementary("A -> B".into(), vec![1.0, -1.4, -1000.0], None);
        assert!(
            negative_temperature_exponent_and_activation_energy
                .K_const(300.0, None, None)
                .unwrap()
                .is_finite()
        );
        let zero_pre_exponential =
            ReactionData::new_elementary("A -> B".into(), vec![0.0, 0.0, 0.0], None);
        assert_eq!(
            zero_pre_exponential.K_const(300.0, None, None).unwrap(),
            0.0
        );
        let underflowing_rate =
            ReactionData::new_elementary("A -> B".into(), vec![1.0, 0.0, 1.0e308], None);
        assert_eq!(underflowing_rate.K_const(300.0, None, None).unwrap(), 0.0);

        let bad_troe = ReactionData::new_falloff(
            "A (+M) -> B (+M)".into(),
            vec![1.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            None,
            Some(vec![0.5, 1000.0]),
            None,
        );
        let concentrations = HashMap::from([("M".to_string(), 1.0)]);
        assert!(
            bad_troe
                .K_const(1000.0, None, Some(concentrations))
                .is_err()
        );
        let invalid_troe_temperature = ReactionData::new_falloff(
            "A (+M) -> B (+M)".into(),
            vec![1.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            None,
            Some(vec![1.1, 1000.0, 5000.0]),
            None,
        );
        assert!(
            invalid_troe_temperature
                .K_const(1000.0, None, Some(HashMap::from([("M".into(), 1.0)])))
                .is_err()
        );
        assert!(
            bad_troe
                .K_sym(
                    None,
                    Some(HashMap::from([("M".to_string(), Expr::Const(1.0))]))
                )
                .is_err()
        );

        let three_body = ReactionData::new_three_body(
            "A + M -> B + M".into(),
            vec![1.0, 0.0, 0.0],
            HashMap::new(),
            None,
        );
        assert!(
            three_body
                .K_const(1000.0, None, Some(HashMap::from([("M".to_string(), -1.0)])))
                .is_err()
        );
        assert!(
            three_body
                .K_const(
                    1000.0,
                    None,
                    Some(HashMap::from([("M".to_string(), f64::INFINITY)]))
                )
                .is_err()
        );
        let negative_efficiency = ReactionData::new_three_body(
            "A + M -> B + M".into(),
            vec![1.0, 0.0, 0.0],
            HashMap::from([("M".into(), -0.1)]),
            None,
        );
        assert!(
            negative_efficiency
                .K_const(1000.0, None, Some(HashMap::from([("M".into(), 1.0)])))
                .is_err()
        );

        let pressure = ReactionData::new_pressure(
            "A -> B".into(),
            HashMap::from([("1.0".into(), vec![1.0, 0.0, 0.0])]),
            None,
        );
        assert!(pressure.K_const(300.0, Some(0.0), None).is_err());
        assert!(pressure.K_sym(Some(f64::INFINITY), None).is_err());
        let malformed_pressure_node = ReactionData::new_pressure(
            "A -> B".into(),
            HashMap::from([("1.0".into(), vec![1.0, 0.0])]),
            None,
        );
        assert!(
            malformed_pressure_node
                .K_const(300.0, Some(1.0), None)
                .is_err()
        );
        let empty_pressure_table =
            ReactionData::new_pressure("A -> B".into(), HashMap::new(), None);
        assert!(
            empty_pressure_table
                .K_const(300.0, Some(1.0), None)
                .is_err()
        );
        let invalid_pressure_key = ReactionData::new_pressure(
            "A -> B".into(),
            HashMap::from([("not-a-pressure".into(), vec![1.0, 0.0, 0.0])]),
            None,
        );
        assert!(
            invalid_pressure_key
                .K_const(300.0, Some(1.0), None)
                .is_err()
        );
        let duplicate_pressure_keys = ReactionData::new_pressure(
            "A -> B".into(),
            HashMap::from([
                ("1".into(), vec![1.0, 0.0, 0.0]),
                ("1.0".into(), vec![2.0, 0.0, 0.0]),
            ]),
            None,
        );
        assert!(
            duplicate_pressure_keys
                .K_const(300.0, Some(1.0), None)
                .is_err()
        );

        let falloff = ReactionData::new_falloff(
            "A (+M) -> B (+M)".into(),
            vec![1.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            None,
            None,
            None,
        );
        assert!(
            falloff
                .K_sym(
                    None,
                    Some(HashMap::from([("M".into(), Expr::Const(f64::NAN))]))
                )
                .is_err()
        );
        let nested_non_finite = Expr::Add(
            Box::new(Expr::Var("M".into())),
            Box::new(Expr::Mul(
                Box::new(Expr::Const(1.0)),
                Box::new(Expr::Const(f64::NAN)),
            )),
        );
        assert!(
            falloff
                .K_sym(None, Some(HashMap::from([("M".into(), nested_non_finite)])))
                .is_err()
        );
        let finite_symbolic =
            Expr::Add(Box::new(Expr::Var("M".into())), Box::new(Expr::Const(1.0)));
        assert!(
            falloff
                .K_sym(None, Some(HashMap::from([("M".into(), finite_symbolic)])))
                .is_ok()
        );

        let overflow = ReactionData::new_elementary("A -> B".into(), vec![1.0e308, 2.0, 0.0], None);
        assert!(matches!(
            overflow.K_const(1.0e308, None, None),
            Err(KineticsError::NonFiniteRate { .. })
        ));
        assert!(ReactionData::create_uniform_vector(0.0, 300.0, 10).is_err());
        assert!(ReactionData::create_uniform_vector(300.0, f64::INFINITY, 10).is_err());

        let mut mismatched = valid;
        mismatched.reaction_type = ReactionType::Pressure;
        assert!(mismatched.K_const(300.0, None, None).is_err());
        mismatched.reaction_type = ReactionType::Empirical;
        assert!(matches!(
            mismatched.K_const(300.0, None, None),
            Err(KineticsError::UnsupportedReactionModel(_))
        ));
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

    #[test]
    fn test_threebody_symbolic_matches_numeric() {
        let threebody_reaction = ThreeBodyStruct::new(
            vec![1.0, 2.0, 300.0],
            HashMap::from([("H2".to_string(), 2.0), ("M".to_string(), 1.0)]),
        );

        let temp = 298.0;
        let concentrations = HashMap::from([
            ("H2".to_string(), Expr::Const(2.0)),
            ("M".to_string(), Expr::Const(1.0)),
        ]);

        let expected_k_const = threebody_reaction.K_const(
            temp,
            HashMap::from([("H2".to_string(), 2.0), ("M".to_string(), 1.0)]),
        );
        let symbolic = threebody_reaction
            .K_expr(concentrations)
            .set_variable("T", temp);
        let symbolic_value = eval_expr(&symbolic);

        assert_relative_eq!(symbolic_value, expected_k_const, epsilon = 1e-10);
    }

    #[test]
    fn test_pressure_reaction_is_order_independent_and_interpolates() {
        let mut arrenius = HashMap::new();
        arrenius.insert("100.0".to_string(), vec![3.0, 1.0, 200.0]);
        arrenius.insert("1.0".to_string(), vec![1.0, 2.0, 100.0]);
        arrenius.insert("10.0".to_string(), vec![2.0, 0.5, 150.0]);

        let pressure_reaction = PressureStruct::new(arrenius);
        let temp = 500.0;

        let k_exact_low = pressure_reaction.K_const(temp, 1.0);
        let k_exact_mid = pressure_reaction.K_const(temp, 10.0);
        let k_exact_high = pressure_reaction.K_const(temp, 100.0);
        assert_relative_eq!(
            k_exact_low,
            calculate_k(&[1.0, 2.0, 100.0], &temp),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            k_exact_mid,
            calculate_k(&[2.0, 0.5, 150.0], &temp),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            k_exact_high,
            calculate_k(&[3.0, 1.0, 200.0], &temp),
            epsilon = 1e-12
        );

        let expected_interp = interpolate(
            5.0,
            1.0,
            10.0,
            calculate_k(&[1.0, 2.0, 100.0], &temp),
            calculate_k(&[2.0, 0.5, 150.0], &temp),
        );
        let k_interp = pressure_reaction.K_const(temp, 5.0);
        assert_relative_eq!(k_interp, expected_interp, epsilon = 1e-12);

        let low_clamped = pressure_reaction.K_const(temp, 0.1);
        assert_relative_eq!(low_clamped, k_exact_low, epsilon = 1e-12);
        let high_clamped = pressure_reaction.K_const(temp, 1000.0);
        assert_relative_eq!(high_clamped, k_exact_high, epsilon = 1e-12);
    }

    #[test]
    fn test_pressure_symbolic_matches_numeric() {
        let mut arrenius = HashMap::new();
        arrenius.insert("10.0".to_string(), vec![2.0, 0.5, 150.0]);
        arrenius.insert("1.0".to_string(), vec![1.0, 2.0, 100.0]);

        let pressure_reaction = PressureStruct::new(arrenius);
        let temp = 500.0;
        let symbolic = pressure_reaction.K_expr(5.0).set_variable("T", temp);
        let symbolic_value = eval_expr(&symbolic);
        let numeric_value = pressure_reaction.K_const(temp, 5.0);

        assert_relative_eq!(symbolic_value, numeric_value, epsilon = 1e-12);
    }

    #[test]
    fn test_create_uniform_vector_rejects_single_point() {
        let result = ReactionData::create_uniform_vector(300.0, 400.0, 1);

        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            KineticsError::InvalidReactionData(message)
            if message.contains("n must be at least 2")
        ));
    }

    #[test]
    fn temperature_sweep_rejects_requests_above_configured_point_limit() {
        let limit = crate::Kinetics::mechfinder_api::MAX_TEMPERATURE_SWEEP_POINTS;
        let requested = limit + 1;
        assert!(matches!(
            ReactionData::create_uniform_vector(300.0, 400.0, requested),
            Err(KineticsError::ResourceLimitExceeded {
                resource: "temperature sweep",
                requested: count,
                limit: actual_limit
            }) if count == requested && actual_limit == limit
        ));
    }

    #[test]
    fn test_k_const_for_t_range_rejects_single_point() {
        let reaction =
            ReactionData::new_elementary("A -> B".to_string(), vec![1.0, 0.0, 100.0], None);

        let result = reaction.K_const_for_T_range(300.0, 400.0, 1, None, None);

        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            KineticsError::InvalidReactionData(message)
            if message.contains("n must be at least 2")
        ));
    }

    #[test]
    fn test_validate_reaction_type_rejects_mismatch() {
        // Keep the reaction payload structurally valid but pair it with the wrong type tag.
        let reaction = ReactionData {
            reaction_type: crate::Kinetics::mechfinder_api::ReactionType::Falloff,
            eq: "A -> B".to_string(),
            react: None,
            data: crate::Kinetics::mechfinder_api::ReactionKinetics::Elementary(
                crate::Kinetics::mechfinder_api::ElementaryStruct::new(vec![1.0, 0.0, 100.0]),
            ),
        };

        let err = reaction.validate_reaction_type().unwrap_err();
        assert!(matches!(
            err,
            KineticsError::InvalidReactionData(message)
            if message.contains("Mismatch between reaction_type")
        ));
    }

    #[test]
    fn test_parse_kinetic_data_vec_rejects_malformed_payload() {
        // Missing the kinetic payload should fail before any state is cached.
        let malformed = vec![serde_json::json!({
            "type": "elem",
            "eq": "A -> B"
        })];

        let err = parse_kinetic_data_vec(malformed).unwrap_err();
        assert!(matches!(
            &err,
            KineticsError::Context { context, .. } if context.contains("index 0")
        ));
        let KineticsError::Context { source, .. } = &err else {
            panic!("missing reaction index context: {err}");
        };
        assert!(matches!(**source, KineticsError::Json(_)));
    }
}
