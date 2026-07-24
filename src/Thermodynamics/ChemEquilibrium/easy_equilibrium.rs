//! Simplified equilibrium interface using the law of mass action.
//!
//! # Purpose
//!
//! This module provides a lightweight, single-reaction equilibrium calculator
//! that works directly with the **law of mass action** rather than the full
//! Gibbs free energy minimization used elsewhere in the ChemEquilibrium module.
//! It is suitable for quick estimates, educational examples, and problems where
//! a single independent reaction dominates the chemistry.
//!
//! # Physical and Mathematical Background
//!
//! For a single reaction `Σ ν_i A_i = 0`, the equilibrium constant `K(T)` is
//! related to the standard Gibbs free energy change:
//!
//! ```text
//! ΔG°(T) = Σ ν_i · ΔG_f°(T, A_i)
//! K(T) = exp(-ΔG°(T) / (R·T))
//! ```
//!
//! The extent `ξ` satisfies the mass-action equation:
//!
//! ```text
//! K(T) = Π (n_i^0 + ν_i·ξ)^(ν_i) · (P / (P0 · n_total))^(Σ ν_i)
//! ```
//!
//! where `n_i^0` are initial moles, `ν_i` are stoichiometric coefficients,
//! `P` is pressure, `P0` is reference pressure, and `n_total` is the total
//! mole number. This module solves for `ξ` using a scalar root finder.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`EasyEquilibrium`] | Single-reaction equilibrium container |
//!
//! # Dataflow
//!
//! ```text
//!   User provides:
//!     P, subs_data, subs_coeffs (HashMap<species, coeff>), initial_moles
//!     │
//!     v
//!   EasyEquilibrium::new()
//!     ├── stores P/101325 (dimensionless pressure)
//!     ├── computes n_total = Σ initial_moles
//!     └── computes dnu = Σ ν_i
//!     │
//!     v
//!   create_equilibrium_const()
//!     ├── builds symbolic expression for ΔG°(T) from NASA polynomials
//!     ├── builds symbolic K(T) = exp(-ΔG°/(R·T))
//!     └── returns closure: f(ξ) = K - Π (n_i^0 + ν_i·ξ)^(ν_i) · (P/P0)^(dnu) / n_total^(dnu)
//!     │
//!     v
//!   ScalarRootFinder solves f(ξ) = 0 for extent ξ
//! ```
//!
//! # Limitations
//!
//! - Only handles **one independent reaction** at a time.
//! - Uses symbolic expressions via `RustedSciThe` for ΔG°(T) computation.
//! - Does not support multiphase or activity coefficient models.
//! - For complex multi-reaction systems, use [`EquilibriumLogMoles`](super::equilibrium_log_moles) instead.
//!
//! # Examples
//!
//! ```rust, ignore
//! use KiThe::Thermodynamics::ChemEquilibrium::easy_equilibrium::EasyEquilibrium;
//!
//! let eq = EasyEquilibrium::new(
//!     101325.0,                          // pressure in Pa
//!     subs_data,                         // substance database
//!     HashMap::from([("CO2".into(), 1.0), ("CO".into(), -1.0), ("O2".into(), -0.5)]),
//!     HashMap::from([("CO2".into(), 1.0), ("CO".into(), 0.0), ("O2".into(), 0.0)]),
//! );
//! let k_eq = eq.create_equilibrium_const();
//! let xi = ScalarRootFinder::new(k_eq, 0.0, 1.0).solve();
//! ```
//!
use crate::Thermodynamics::User_substances::{DataType, SubsData};
use RustedSciThe::numerical::optimization::minimize_scalar::{ClosureFunction, ScalarRootFinder};
use RustedSciThe::symbolic::symbolic_engine::Expr;

use std::collections::HashMap;

/// Simplified single-reaction equilibrium calculator using the law of mass action.
///
/// Fields are public for legacy callers that construct this struct directly.
/// New code should prefer [`EasyEquilibrium::new()`].
#[derive(Clone, Debug)]
pub struct EasyEquilibrium {
    /// System pressure in Pascals (Pa).
    pub P: f64,
    /// Substance database with thermochemical data for all species.
    pub subs_data: SubsData,
    /// Stoichiometric coefficients keyed by substance name.
    /// Positive for products, negative for reactants.
    pub subs_coeffs: HashMap<String, f64>,
    /// Initial physical mole numbers keyed by substance name.
    pub initial_moles: HashMap<String, f64>,
    /// Sum of stoichiometric coefficients `Σ ν_i` (used for mole-fraction corrections).
    dnu: f64,
    /// Total initial moles `Σ n_i^0` (used for mole-fraction normalization).
    n_total: f64,
}

impl EasyEquilibrium {
    pub fn new(
        P: f64,
        subs_data: SubsData,
        subs_coeffs: HashMap<String, f64>,
        initial_moles: HashMap<String, f64>,
    ) -> Self {
        let P = P / 101325.0;
        let mut n_total = 0.0;
        for (_, n0i) in initial_moles.clone() {
            n_total += n0i;
        }
        let dnu = subs_coeffs.values().sum();
        EasyEquilibrium {
            P,
            subs_data,
            subs_coeffs,
            initial_moles,
            dnu,
            n_total,
        }
    }
    pub fn create_equilibrium_const(&self) -> Box<dyn Fn(f64) -> f64> {
        let subs_data = self.subs_data.clone();
        let subs_coeffs = self.subs_coeffs.clone();

        let T = Expr::Var("T".to_string());
        let mut G_react = Expr::Const(0.0);
        for (subs_i, coeff_i) in subs_coeffs {
            let sym_data = subs_data.therm_map_of_sym.get(&subs_i).unwrap();
            let dS_i = *sym_data
                .get(&DataType::dS_sym)
                .clone()
                .unwrap()
                .clone()
                .unwrap();
            let dH_i = *sym_data
                .get(&DataType::dH_sym)
                .clone()
                .unwrap()
                .clone()
                .unwrap();
            let dG_i = dH_i - T.clone() * dS_i;
            G_react += dG_i * Expr::Const(coeff_i);
        }
        let G_react_simplified = G_react.simplify();
        println!("G_react_simplified: {}", G_react_simplified.pretty_print());
        let dG_fun = G_react_simplified.lambdify1D();
        let K = Box::new(move |T| f64::exp(-dG_fun(T) / (8.31446261815324 * T)));
        K
    }

    pub fn eta_equauion(&mut self) -> Expr {
        let mut LHS = Expr::Const(1.0);
        let subs_coeffs = self.subs_coeffs.clone();
        let eta = Expr::Var("eta".to_string());
        let _n_total = Expr::Const(self.n_total.clone());
        let _dnu = Expr::Const(self.dnu.clone());
        let P_pow_dnu = Expr::Const(self.P.clone().powf(self.dnu));
        for (subs_i, n0i) in self.initial_moles.clone() {
            let coeff_i = subs_coeffs.get(&subs_i).unwrap();
            LHS *= (Expr::Const(n0i) + Expr::Const(coeff_i.clone()) * eta.clone())
                .pow(Expr::Const(*coeff_i));
            //println!("LHS: {}", LHS.pretty_print());
        }
        let LHS = P_pow_dnu * LHS.simplify(); // / (n_total + dnu.clone() * eta).pow(dnu);
        LHS
    }

    pub fn solve(&mut self, T0: f64, T_step: f64, T_end: f64) -> HashMap<String, Vec<f64>> {
        let K = self.create_equilibrium_const();
        let LHS = self.eta_equauion();
        let LHS_fun = LHS.lambdify1D();
        let solver = ScalarRootFinder::new();
        let mut map_of_solution: HashMap<String, Vec<f64>> = self
            .initial_moles
            .clone()
            .into_iter()
            .map(|(k, _v)| (k, Vec::new()))
            .collect();

        let mut eta_vector = Vec::new();
        let mut T = T0;
        while T < T_end {
            let K = K(T);
            let eq = |eta| LHS_fun(eta) - K;
            let func = ClosureFunction::new(eq, "".to_string());

            let eta = solver.secant(&func, 0.2, 0.7).unwrap();
            let eta = eta.root;
            dbg!(eta);
            for subs_i in self.initial_moles.keys() {
                let n0i = self.initial_moles.get(subs_i).unwrap();
                let coeff_i = self.subs_coeffs.get(subs_i).unwrap();
                let n_i = n0i + coeff_i * eta;
                map_of_solution.get_mut(subs_i).unwrap().push(n_i);
            }
            eta_vector.push(eta);
            T += T_step;
        }
        map_of_solution
    }
}

//////////////////////////////////////////TESTS/////////////////////////////////////

#[cfg(test)]
mod easy_equilibrium_tests {
    use super::*;
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases};

    #[test]
    fn test_easy_equilibrium_basic() {
        let subs = vec!["H2".to_string(), "H".to_string()];
        let mut subdata = SubsData::new();
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);
        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();
        let _ = subdata.search_substances();
        let _ = subdata.extract_thermal_coeffs(subs[0].as_str(), 3000.0);
        let _ = subdata.extract_thermal_coeffs(subs[1].as_str(), 3000.0);
        let _ = subdata.calculate_therm_map_of_sym();

        let mut subs_coeffs = HashMap::new();
        subs_coeffs.insert("H2".to_string(), -1.0);
        subs_coeffs.insert("H".to_string(), 2.0);

        let mut initial_moles = HashMap::new();
        initial_moles.insert("H2".to_string(), 1.0);
        initial_moles.insert("H".to_string(), 0.0);

        let eq = EasyEquilibrium::new(101325.0, subdata, subs_coeffs, initial_moles);
        assert_eq!(eq.P, 1.0);
        assert_eq!(eq.dnu, 1.0);
        assert_eq!(eq.n_total, 1.0);
    }
    #[test]
    fn test_easy_equilibrium_basic2() {
        let subs = vec!["N2".to_string(), "N".to_string()];
        let mut subdata = SubsData::new();
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);
        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();
        let _ = subdata.search_substances();
        let _ = subdata.extract_thermal_coeffs(subs[0].as_str(), 3000.0);
        let _ = subdata.extract_thermal_coeffs(subs[1].as_str(), 3000.0);
        let _ = subdata.calculate_therm_map_of_sym();

        let mut subs_coeffs = HashMap::new();
        subs_coeffs.insert("N2".to_string(), -1.0);
        subs_coeffs.insert("N".to_string(), 2.0);

        let mut initial_moles = HashMap::new();
        initial_moles.insert("N2".to_string(), 1.0);
        initial_moles.insert("N".to_string(), 0.0);

        let mut eq = EasyEquilibrium::new(101325.0, subdata, subs_coeffs, initial_moles);
        let solutions = eq.solve(400.0, 100.0, 5500.0);
        let n2_sol = solutions.get("N2").unwrap();
        let n_sol = solutions.get("N").unwrap();
        println!("T\t N2\t N");
        let mut T = 400.0;
        for i in 0..n2_sol.len() {
            println!("{}\t {:.6}\t {:.6}", T, n2_sol[i], n_sol[i]);
            T += 100.0;
        }

        assert_eq!(eq.P, 1.0);
        assert_eq!(eq.dnu, 1.0);
        assert_eq!(eq.n_total, 1.0);
    }
}
