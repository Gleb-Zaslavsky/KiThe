//! # Single-Phase Thermodynamics Module
//!
//! This module provides optimized thermodynamic calculations for single-phase systems.
//! It offers better performance and simpler API compared to multi-phase systems when
//! only one phase is involved.
//!
//! ## Key Structure
//!
//! - [`OnePhase`]: Manages a single phase with direct access to substance data
//!
//! ## None Key Convention
//!
//! To maintain API compatibility with multi-phase systems, this module uses `None`
//! as the phase key in all HashMap returns. This allows seamless integration with
//! the unified [`CustomSubstance`] interface.
//!
//! ## Example Usage
//!
//! ```rust
//! use std::collections::HashMap;
//!
//! // Single-phase gas system
//! let mut system = OnePhase::new();
//! // system.subs_data.substances = vec!["CO2".to_string(), "H2O".to_string()];
//!
//! // Calculate properties - results use None key for compatibility
//! // let result = system.calculate_Gibbs_sym(298.15)?;
//! // assert!(result.contains_key(&None));
//! ```
//!
//! ## Performance Benefits
//!
//! - Direct field access instead of HashMap lookups
//! - No phase iteration overhead
//! - Optimized data structures for single-phase operations

use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use std::fmt;

use crate::Thermodynamics::User_PhaseOrSolution::{
    R, R_sym, SubstancesContainer, ThermodynamicsCalculatorTrait,
};
use std::collections::HashMap;
use std::f64;

/// Single-phase thermodynamic system with optimized data access.
///
/// Provides direct access to substance data without phase-level indirection.
/// All HashMap returns use `None` as the key for API compatibility with multi-phase systems.
///
/// # Example
/// ```rust
/// let mut system = OnePhase::new();
/// // Direct access to substance data
/// // system.subs_data.substances.push("CO2".to_string());
/// ```
pub struct OnePhase {
    pub subs_data: SubsData,
    pub symbolic_vars: (Option<Expr>, Option<Vec<Expr>>),
    pub vec_of_n_vars: Vec<Expr>,
    pub Np: Expr,
    pub map_of_var_each_substance: HashMap<String, Expr>,
    /// hashmap of Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    pub dG: HashMap<String, f64>,
    /// hashmap of Gibbs free energy of a given mixure of substances (function at given T, P, concentration)
    pub dG_fun: HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    /// hashmap of Gibbs free energy of a given mixure of substances (symbolic result at given T, P, concentration)
    pub dG_sym: HashMap<String, Expr>,

    /// hashmap of entropy of a given mixure of substances (numerical result at given T, P, concentration)
    pub dS: HashMap<String, f64>,
    /// hashmap of entropy functions
    pub dS_fun: HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    /// hashmap of entropy of a given mixure of substances (symbolic result at given T, P, concentration)
    pub dS_sym: HashMap<String, Expr>,
}
impl OnePhase {
    /// Creates a new single-phase system with empty data.
    pub fn new() -> Self {
        OnePhase {
            subs_data: SubsData::new(),
            symbolic_vars: (None, None),
            vec_of_n_vars: Vec::new(),
            Np: Expr::Var("Np".to_string()),
            map_of_var_each_substance: HashMap::new(),
            dG: HashMap::new(),
            dG_fun: HashMap::new(),
            dG_sym: HashMap::new(),
            dS: HashMap::new(),
            dS_fun: HashMap::new(),
            dS_sym: HashMap::new(),
        }
    }

    /// Wraps any result with None key for trait compatibility.
    /// Uses optimized HashMap capacity for single-item storage.
    fn wrap_result<T>(&self, result: T) -> HashMap<Option<String>, T> {
        let mut outer = HashMap::with_capacity(1);
        outer.insert(None, result);
        outer
    }

    /// Extracts data from trait's expected format with proper error handling.
    /// Looks for None key which represents the single phase.
    fn extract_phase_data<'a, T>(
        &self,
        n: &'a HashMap<Option<String>, T>,
    ) -> SubsDataResult<&'a T> {
        n.get(&None).ok_or_else(move || SubsDataError::MissingData {
            field: "phase data".to_string(),
            substance: "None".to_string(),
        })
    }

    fn indexed_moles_variables_local(
        &mut self,
        Np: Option<Expr>,
    ) -> Result<
        (
            (Option<Expr>, Option<Vec<Expr>>),
            Vec<Expr>,
            Vec<Expr>,
            HashMap<String, Expr>,
        ),
        SubsDataError,
    > {
        let n = Expr::IndexedVars(self.subs_data.substances.len(), "N").0;
        let Np = Np.unwrap_or_else(|| Expr::Var("Np".to_string()));

        let mut map_of_var_each_substance = HashMap::with_capacity(self.subs_data.substances.len());
        for (i, substance) in self.subs_data.substances.iter().enumerate() {
            let var = Expr::Var(format!("N{}", i));
            map_of_var_each_substance.insert(substance.clone(), var);
        }

        self.symbolic_vars = (Some(Np.clone()), Some(n.clone()));
        self.vec_of_n_vars = n.clone();
        self.Np = Np.clone();
        self.map_of_var_each_substance = map_of_var_each_substance.clone();

        Ok((
            (Some(Np.clone()), Some(n.clone())),
            n,
            vec![Np],
            map_of_var_each_substance.clone(),
        ))
    }

    /// Calculates Gibbs free energy for given conditions and mole numbers.
    pub fn calcutate_Gibbs_free_energy_local(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> SubsDataResult<()> {
        let n_Np = self.extract_phase_data(&n)?.clone();
        let Np = n_Np.0;
        let n = n_Np.1;
        let dG_phase = self.subs_data.calc_dG_for_one_phase(P, T, n, Np);
        self.dG = dG_phase.clone();
        Ok(())
    }

    pub fn calculate_Gibbs_fun_unmut(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>> {
        self.subs_data.calculate_Gibbs_fun_one_phase(P, T)
    }
    /// Returns entropy functions wrapped for trait compatibility.
    pub fn calculate_S_fun_unmut(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>> {
        self.subs_data.calculate_S_fun_for_one_phase(P, T)
    }

    pub fn calculate_Lagrange_equations_fun(
        &mut self,
        A: DMatrix<f64>,
        G_fun: HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,

        Tm: f64,
    ) -> Result<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
        SubsDataError,
    > {
        let G = G_fun;
        let subs = self.subs_data.substances.clone();
        let fun =
            move |T: f64, n: Option<Vec<f64>>, Np: Option<f64>, Lambda: Vec<f64>| -> Vec<f64> {
                let A = A.transpose();
                let n_substances = A.ncols();

                let mut vec_of_eqs: Vec<f64> = Vec::new();

                for i in 0..n_substances {
                    let subst_i = subs[i].clone();
                    let G_i = G.get(&subst_i).unwrap();
                    let col_of_subs_i = A.column(i).iter().map(|&v| v).collect::<Vec<f64>>();
                    let sum_by_elemnts: f64 = col_of_subs_i
                        .iter()
                        .zip(Lambda.iter())
                        .map(|(aij, Lambda_j)| aij * Lambda_j)
                        .fold(0.0, |acc, x| acc + x);
                    let eq_i = sum_by_elemnts + G_i(T, n.clone(), Np.clone()) / (R * Tm);
                    vec_of_eqs.push(eq_i);
                }
                vec_of_eqs
            };
        Ok(Box::new(fun))
    }
}

impl Clone for OnePhase {
    fn clone(&self) -> Self {
        let mut new_dG_fun = HashMap::new();
        for (substance, _) in &self.dG_fun {
            let placeholder: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static> =
                Box::new(|_: f64, _: Option<Vec<f64>>, _: Option<f64>| 0.0);
            new_dG_fun.insert(substance.clone(), placeholder);
        }

        let mut new_dS_fun = HashMap::new();
        for (substance, _) in &self.dS_fun {
            let placeholder: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static> =
                Box::new(|_: f64, _: Option<Vec<f64>>, _: Option<f64>| 0.0);
            new_dS_fun.insert(substance.clone(), placeholder);
        }

        Self {
            subs_data: self.subs_data.clone(),
            symbolic_vars: self.symbolic_vars.clone(),
            vec_of_n_vars: self.vec_of_n_vars.clone(),
            Np: self.Np.clone(),
            map_of_var_each_substance: self.map_of_var_each_substance.clone(),
            dG: self.dG.clone(),
            dG_fun: new_dG_fun,
            dG_sym: self.dG_sym.clone(),
            dS: self.dS.clone(),
            dS_fun: new_dS_fun,
            dS_sym: self.dS_sym.clone(),
        }
    }
}

impl fmt::Debug for OnePhase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("OnePhase")
            .field("subs_data", &self.subs_data)
            .field("symbolic_vars", &self.symbolic_vars)
            .field("vec_of_n_vars", &self.vec_of_n_vars)
            .field("Np", &self.Np)
            .field("map_of_var_each_substance", &self.map_of_var_each_substance)
            .field("dG", &self.dG)
            .field("dG_sym", &self.dG_sym)
            .field("dS", &self.dS)
            .field("dS_sym", &self.dS_sym)
            .finish()
    }
}
/// Implementation of ThermodynamicsCalculatorTrait for single-phase systems.
/// All methods maintain compatibility with multi-phase interface using None keys.
impl ThermodynamicsCalculatorTrait for OnePhase {
    /// Extracts thermal coefficients for the given temperature.
    fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.subs_data.extract_all_thermal_coeffs(temperature)?;
        Ok(())
    }

    /// Calculates thermodynamic property maps at given temperature.
    fn calculate_therm_map_of_properties(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.subs_data
            .calculate_therm_map_of_properties(temperature)?;
        Ok(())
    }

    /// Creates symbolic expressions for thermodynamic properties.
    fn calculate_therm_map_of_sym(&mut self) -> SubsDataResult<()> {
        self.subs_data.calculate_therm_map_of_sym()?;
        Ok(())
    }

    /// Validates and extracts coefficients for temperature range.
    fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>> {
        let v = self
            .subs_data
            .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(temperature)?;
        Ok(v)
    }

    /// Calculates Gibbs free energy for given conditions and mole numbers.
    fn calcutate_Gibbs_free_energy(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        self.calcutate_Gibbs_free_energy_local(T, P, n)?;
        let dG_phase = self.dG.clone();
        Ok(self.wrap_result(dG_phase))
    }

    /// Creates symbolic expressions for Gibbs free energy.
    fn calculate_Gibbs_sym(&mut self, T: f64) -> SubsDataResult<()> {
        let Np = self.symbolic_vars.0.clone();
        let n = self.symbolic_vars.1.clone();
        self.dG_sym = self.subs_data.calculate_Gibbs_sym_one_phase(T, n, Np);
        Ok(())
    }

    /// Sets pressure value in symbolic Gibbs expressions.
    fn set_P_to_sym_in_G_sym(&mut self, P: f64) {
        for (_, sym_fun) in self.dG_sym.iter_mut() {
            *sym_fun = sym_fun.set_variable("P", P).simplify()
        }
    }

    /// Sets temperature value in symbolic Gibbs expressions.
    fn set_T_to_sym_in_G_sym(&mut self, T: f64) {
        for (_, sym_fun) in self.dG_sym.iter_mut() {
            *sym_fun = sym_fun.set_variable("T", T).simplify()
        }
    }

    /// Creates Gibbs free energy functions for given conditions.
    fn calculate_Gibbs_fun(&mut self, T: f64, P: f64) {
        self.dG_fun = self.subs_data.calculate_Gibbs_fun_one_phase(P, T);
    }

    /// Calculates entropy for given conditions and mole numbers.
    fn calculate_S(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> SubsDataResult<()> {
        let n_Np = self.extract_phase_data(&n)?.clone();
        let Np = n_Np.0;
        let n = n_Np.1;
        self.dS = self.subs_data.calculate_S_for_one_phase(P, T, n, Np);
        Ok(())
    }

    /// Creates symbolic expressions for entropy.
    fn calculate_S_sym(&mut self, T: f64) -> SubsDataResult<()> {
        let Np = self.symbolic_vars.0.clone();
        let n = self.symbolic_vars.1.clone();
        self.dS_sym = self.subs_data.calculate_S_sym_for_one_phase(T, n, Np);
        Ok(())
    }

    /// Creates entropy functions for given conditions.
    fn calculate_S_fun(&mut self, T: f64, P: f64) {
        self.dS_fun = self.subs_data.calculate_S_fun_for_one_phase(P, T);
    }

    fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()> {
        // Set pressure and molar masses for each phase
        self.subs_data.set_P(pressure, pressure_unit.clone());
        self.subs_data
            .set_M(molar_masses.clone(), mass_unit.clone());
        Ok(())
    }
    fn if_not_found_go_NIST(&mut self) -> SubsDataResult<()> {
        self.subs_data.if_not_found_go_NIST()?;
        Ok(())
    }
    fn extract_SubstancesContainer(&mut self) -> SubsDataResult<SubstancesContainer> {
        let substances = self.subs_data.substances.clone();
        Ok(SubstancesContainer::SinglePhase(substances))
    }
    fn get_all_substances(&mut self) -> Vec<String> {
        self.subs_data.substances.clone()
    }
    fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> SubsDataResult<(DMatrix<f64>, HashMap<String, f64>, Vec<String>)> {
        self.subs_data
            .calculate_elem_composition_and_molar_mass(groups)?;
        let elem_composition_matrix =
            self.subs_data
                .elem_composition_matrix
                .clone()
                .ok_or_else(|| {
                    SubsDataError::MatrixOperationFailed(
                        "Element composition matrix not calculated".to_string(),
                    )
                })?;

        let unique_elements = self.subs_data.unique_elements.clone();
        let hashmap_of_molar_masses = self.subs_data.hasmap_of_molar_mass.clone();
        Ok((
            elem_composition_matrix,
            hashmap_of_molar_masses,
            unique_elements,
        ))
    }
    fn indexed_moles_variables(
        &mut self,
    ) -> Result<
        (
            HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
            Vec<Expr>,
            Vec<Expr>,
            HashMap<Option<String>, HashMap<String, Expr>>,
        ),
        SubsDataError,
    > {
        let n = Expr::IndexedVars(self.subs_data.substances.len(), "N").0;
        let Np = Expr::Var("Np".to_string());

        let mut map_of_var_each_substance = HashMap::with_capacity(self.subs_data.substances.len());
        for (i, substance) in self.subs_data.substances.iter().enumerate() {
            let var = Expr::Var(format!("N{}", i));
            map_of_var_each_substance.insert(substance.clone(), var);
        }

        self.symbolic_vars = (Some(Np.clone()), Some(n.clone()));
        self.vec_of_n_vars = n.clone();
        self.Np = Np.clone();
        self.map_of_var_each_substance = map_of_var_each_substance.clone();

        Ok((
            self.wrap_result((Some(Np.clone()), Some(n.clone()))),
            n,
            vec![Np],
            self.wrap_result(map_of_var_each_substance),
        ))
    }
    /////////////////////////////equations for G->min///////////////////////////////////////
    fn calculate_Lagrange_equations_sym(
        &mut self,
        A: DMatrix<f64>,
        Tm: f64,
    ) -> Result<Vec<Expr>, SubsDataError> {
        let A = A.clone().transpose();
        let n_elements = A.nrows();
        let Lambda = Expr::IndexedVars(n_elements, "Lambda").0;
        let n_substances = A.ncols();
        let subs = self.subs_data.substances.clone();
        let mut vec_of_eqs: Vec<Expr> = Vec::new();

        for i in 0..n_substances {
            let Tm = Expr::Const(Tm);
            let subst_i = subs[i].clone();
            let G_i = self.dG_sym.get(&subst_i).unwrap().clone();
            let col_of_subs_i = A
                .column(i)
                .iter()
                .map(|&v| Expr::Const(v))
                .collect::<Vec<Expr>>();
            let sum_by_elemnts: Expr = col_of_subs_i
                .iter()
                .zip(Lambda.iter())
                .map(|(aij, Lambda_j)| aij.clone() * Lambda_j.clone())
                .fold(Expr::Const(0.0), |acc, x| acc + x)
                .simplify();
            let eq_i = sum_by_elemnts + (G_i.clone() / (R_sym * Tm)).simplify();
            vec_of_eqs.push(eq_i.simplify());
        }
        Ok(vec_of_eqs)
    }

    fn calculate_Lagrange_equations_fun2(
        &mut self,
        A: DMatrix<f64>,
        T: f64,
        P: f64,
        Tm: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    > {
        let dG_fun = self.calculate_Gibbs_fun_unmut(T, P);
        let Lagrange_equations = self.calculate_Lagrange_equations_fun(A, dG_fun, Tm)?;
        Ok(Lagrange_equations)
    }
    //////////////////////////////equations for Keq///////////////////////////////////////
    /*
        fn calculate_K_eq_sym(
            &mut self,
            stolich_matrix: DMatrix<f64>,

            map_of_var_each_substance: HashMap<Option<String>, HashMap<String, Expr>>,
        ) {
            let r = stolich_matrix.ncols();
            let substances = self.substances.clone();
            assert_eq!(r, substances.len());
            let T = Expr::Var("T".to_string());

            let subs_data = self.calculate_Gibbs_sym_one_phase(298.0, None, None);
            let mut K_vec: Vec<Expr> = Vec::with_capacity(r);
            for j in 0..r {
                let mut dG_reaction_j = Expr::Const(0.0);
                for (i, subs_i) in substances.iter().enumerate() {
                    let dG_i = subs_data.get(subs_i).unwrap().clone();
                    let nu_ij = stolich_matrix[(i, j)];
                    dG_reaction_j += Expr::Const(nu_ij) * dG_i;
                }
                let ln_K = -dG_reaction_j / (R_sym * T.clone());
                K_vec.push(ln_K);
            }
        }

        fn calculate_K_eq_fun(
            &mut self,
            stolich_matrix: DMatrix<f64>,
            G_fun: HashMap<
                Option<String>,
                HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
            >,
        ) {
            let r = stolich_matrix.ncols();
            let substances = self.substances.clone();
            assert_eq!(r, substances.len());
            let P = self.P.unwrap_or(101325.0);
            let dG = self.calculate_Gibbs_fun_one_phase(P, 298.0);

            let stol = stolich_matrix; // move into closure
            let K_vec: Box<dyn Fn(f64) -> Vec<f64> + 'static> = Box::new(move |T: f64| -> Vec<f64> {
                let mut K_vec: Vec<f64> = Vec::with_capacity(r);
                for j in 0..r {
                    let mut dG_reaction_j = 0.0;
                    for (i, subs_i) in substances.iter().enumerate() {
                        let dG_i = dG
                            .get(subs_i)
                            .expect("Missing Gibbs function for substance");
                        let nu_ij = stol[(i, j)];
                        dG_reaction_j += nu_ij * dG_i(T, None, None);
                    }
                    let ln_K = -dG_reaction_j / (R * T);
                    K_vec.push(ln_K);
                }
                K_vec
            });
        }
        fn calculate_K_eq_equation_sym(
            &mut self,
            stolich_matrix: DMatrix<f64>,

            map_of_var_each_substance: HashMap<Option<String>, HashMap<String, Expr>>,
        ) {
            let r = stolich_matrix.ncols();
            let substances = self.substances.clone();
            assert_eq!(r, substances.len());
            let T = Expr::Var("T".to_string());
            let map_of_vars = map_of_var_each_substance.get(&None).unwrap().clone();
            let subs_data = self.calculate_Gibbs_sym_one_phase(298.0, None, None);
            let mut residual: Vec<Expr> = Vec::new();
            for j in 0..r {
                let mut dG_reaction_j = Expr::Const(0.0);
                let mut sum_ln_n = Expr::Const(0.0);
                for (i, subs_i) in substances.iter().enumerate() {
                    let dG_i = subs_data.get(subs_i).unwrap().clone();
                    let nu_ij = Expr::Const(stolich_matrix[(i, j)]);
                    dG_reaction_j += nu_ij.clone() * dG_i;

                    let n_i = map_of_vars.get(subs_i).unwrap();
                    sum_ln_n += nu_ij * n_i.ln();
                }
                let ln_K = -dG_reaction_j / (R_sym * T.clone());
                let residual_j = sum_ln_n - sum_ln_n_phase - ln_K;
                residual.push(ln_K);
            }
        }
        fn calculate_K_eq_equation_fun(
            &mut self,
            stolich_matrix: DMatrix<f64>,
            G_fun: HashMap<
                Option<String>,
                HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
            >,
        ) {
            let r = stolich_matrix.ncols();
            let substances = self.substances.clone();
            assert_eq!(r, substances.len());
            let P = self.P.unwrap_or(101325.0);
            let dG = self.calculate_Gibbs_fun_one_phase(P, 298.0);

            let stol = stolich_matrix; // move into closure
            let K_vec: Box<dyn Fn(f64, Vec<f64>, Option<f64>) -> Vec<f64> + 'static> =
                Box::new(move |T: f64, n: Vec<f64>, Np: Option<f64>| -> Vec<f64> {
                    let mut K_vec: Vec<f64> = Vec::with_capacity(r);
                    let mut sum_ln_n = 0.0;
                    for j in 0..r {
                        let mut dG_reaction_j = 0.0;
                        for (i, subs_i) in substances.iter().enumerate() {
                            let dG_i = dG
                                .get(subs_i)
                                .expect("Missing Gibbs function for substance");
                            let nu_ij = stol[(i, j)];
                            dG_reaction_j += nu_ij * dG_i(T, None, None);

                            sum_ln_n += nu_ij * n[i].ln();
                        }
                        let ln_K = -dG_reaction_j / (R * T.clone());
                        let residual_j = sum_ln_n - sum_ln_n_phase - ln_K;
                        residual.push(ln_K);
                    }
                    K_vec
                });
        }
    */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    fn create_full_map_of_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
        // physiscal_nature_of_phase_or_solution: Option<HashMap<String, Phases>>,
    ) -> Result<
        (
            HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
            HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
            HashMap<String, f64>,
        ),
        SubsDataError,
    > {
        let mut full_map_of_mole_numbers = non_zero_number_of_moles.clone();
        let mut map_of_mole_num_vecs: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)> =
            HashMap::new();

        let mut initial_map_of_mole_numbers: HashMap<String, f64> = HashMap::new();
        let subs = self.subs_data.substances.clone();

        for (_, value) in full_map_of_mole_numbers.iter_mut() {
            if let (Np, Some(inner_map)) = (value.0.clone(), &mut value.1) {
                // For each required key, insert with 0.0 if not present
                for key in &subs {
                    inner_map.entry(key.clone()).or_insert(0.0);
                }
                let inner_map = inner_map.clone();
                let vec_of_mole_nums = inner_map.values().cloned().collect::<Vec<f64>>();
                map_of_mole_num_vecs.insert(None, (Np, Some(vec_of_mole_nums)));

                // making initial map of mole numbers
                for (keys, values) in inner_map.iter().clone() {
                    *initial_map_of_mole_numbers
                        .entry(keys.clone())
                        .or_insert(0.0) += values;
                }
            }
        }
        //  println!("full_map_of_mole_numbers,: {:#?},\n map_of_mole_num_vecs {:?}, \n initial_map_of_mole_numbers {:?} ",full_map_of_mole_numbers, map_of_mole_num_vecs, initial_map_of_mole_numbers); panic!();

        Ok((
            full_map_of_mole_numbers,
            map_of_mole_num_vecs,
            initial_map_of_mole_numbers,
        ))
    }
}
