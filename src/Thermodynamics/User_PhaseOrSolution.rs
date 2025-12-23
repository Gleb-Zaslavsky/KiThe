//! # Multi-Phase Thermodynamics Module
//!
//! This module provides thermodynamic calculations for systems with multiple phases or solutions.
//! It handles complex systems where substances can exist in different phases (gas, liquid, solid)
//! or different solutions simultaneously.
//!
//! ## Key Structures
//!
//! - [`PhaseOrSolution`]: Manages multiple phases, each with their own substance data
//! - [`CustomSubstance`]: Enum that unifies single-phase and multi-phase systems
//!
//! ## Phase Key Convention
//!
//! - **Multi-phase systems**: Use `Some("phase_name")` as keys (e.g., `Some("gas")`, `Some("liquid")`)
//! - **Single-phase systems**: Use `None` as the key to represent the single phase
//!
//! ## Example Usage
//!
//! ```rust
//! use std::collections::HashMap;
//! 
//! // Multi-phase system with gas and liquid phases
//! let mut system = PhaseOrSolution::new();
//! // system.subs_data.insert(Some("gas".to_string()), gas_data);
//! // system.subs_data.insert(Some("liquid".to_string()), liquid_data);
//! 
//! // Calculate Gibbs energy for all phases
//! // system.calculate_Gibbs_sym(298.15)?;
//! ```

use std::fmt;

use crate::Thermodynamics::User_PhaseOrSolution2::OnePhase;
use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use enum_dispatch::enum_dispatch;
use nalgebra::{DMatrix, DVector};

use std::collections::{HashMap, HashSet};
use std::f64;
pub const R: f64 = 8.314;
#[allow(non_upper_case_globals)]
pub const R_sym: Expr = Expr::Const(R);
/// Unified interface for single-phase and multi-phase thermodynamic systems.
/// 
/// Uses `None` key for single-phase systems and `Some(phase_name)` for multi-phase.
/// Provides consistent API regardless of system complexity.
#[derive(Debug, Clone)]
#[enum_dispatch(ThermodynamicsCalculatorTrait)]
pub enum CustomSubstance {
    OnePhase(OnePhase),
    PhaseOrSolution(PhaseOrSolution),
}

/// Maps substances to their respective phases while maintaining order.
/// Essential for tracking which substances belong to which phases in multi-phase systems.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubstancePhaseMapping {
    pub all_substances: Vec<String>,
    pub substance_to_phase: Vec<Option<String>>,
}

impl SubstancePhaseMapping {
    /// Creates a new substance-phase mapping.
    /// Substances and phases vectors must have the same length.
    pub fn new(substances: Vec<String>, phases: Vec<Option<String>>) -> Self {
        assert_eq!(substances.len(), phases.len());
        Self {
            all_substances: substances,
            substance_to_phase: phases,
        }
    }

    /// Gets the phase for a substance at the given index.
    pub fn get_phase_for_substance(&self, index: usize) -> Option<&Option<String>> {
        self.substance_to_phase.get(index)
    }
}

/// Container for substance names in single or multiple phases.
/// Provides unified interface for accessing substances regardless of phase structure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubstancesContainer {
    /// Names of substances in a single phase system
    SinglePhase(Vec<String>),
    /// Map of phase names to their substances in multi-phase systems
    MultiPhase(HashMap<String, Vec<String>>),
}

impl SubstancesContainer {
    /// Returns all substances as a flat, sorted, deduplicated vector.
    /// Useful for getting complete substance list regardless of phase structure.
    pub fn get_all_substances(&self) -> Vec<String> {
        let subs = match self {
            SubstancesContainer::SinglePhase(substances) => substances.clone(),
            SubstancesContainer::MultiPhase(phase_substances) => {
                let mut all_substances = Vec::new();
                for substances in phase_substances.values() {
                    all_substances.extend(substances.clone());
                }
                all_substances
            }
        }; //getting rid of duplicates
        let subs: HashSet<String> = HashSet::from_iter(subs.clone());
        let mut subs: Vec<String> = subs.into_iter().collect();
        subs.sort();
        subs
    } // get all substances in all phases
}

impl CustomSubstance {
    /// Extracts substances container from the system.
    /// Returns SinglePhase for both OnePhase and flattened PhaseOrSolution.
    pub fn extract_SubstancesContainer(&mut self) -> SubsDataResult<SubstancesContainer> {
        match self {
            CustomSubstance::OnePhase(subs_data) => Ok(SubstancesContainer::SinglePhase(
                subs_data.subs_data.substances.clone(),
            )),
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let mut all_substances = Vec::new();
                for (_, subs_data) in &phase_or_solution.subs_data {
                    all_substances.extend(subs_data.substances.clone());
                }
                let subs: HashSet<String> = HashSet::from_iter(all_substances.clone());
                let mut subs: Vec<String> = subs.into_iter().collect();
                subs.sort();
                Ok(SubstancesContainer::SinglePhase(subs))
            }
        }
    } // get all substances in all phases

    /// Creates complete mole number maps from sparse input.
    /// Fills missing substances with 0.0 and provides multiple formats for calculations.
    pub fn create_full_map_of_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
        //  physiscal_nature_of_phase_or_solution: Option<HashMap<String, Phases>>,
    ) -> SubsDataResult<(
        HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
        HashMap<String, f64>,
    )> {
        match self {
            CustomSubstance::OnePhase(one_phase) => {
                one_phase.create_full_map_of_mole_numbers(non_zero_number_of_moles)
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                phase_or_solution.create_full_map_of_mole_numbers(non_zero_number_of_moles)
            }
        }
    }
    
    /// Gets mutable reference to Gibbs symbolic expressions.
    /// Returns normalized format with None key for single-phase systems.
    pub fn get_dG_sym_mut(&mut self) -> HashMap<Option<String>, HashMap<String, Expr>> {
        match self {
            CustomSubstance::OnePhase(one_phase) => {
                let mut result = HashMap::new();
                result.insert(None, one_phase.dG_sym.clone());
                result
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => phase_or_solution.dG_sym.clone(),
        }
    }

    /// Gets reference to Gibbs symbolic expressions.
    /// Returns normalized format with None key for single-phase systems.
    pub fn get_dG_sym(&self) -> HashMap<Option<String>, HashMap<String, Expr>> {
        match self {
            CustomSubstance::OnePhase(one_phase) => {
                let mut result = HashMap::new();
                result.insert(None, one_phase.dG_sym.clone());
                result
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => phase_or_solution.dG_sym.clone(),
        }
    }
    

 
    /// Gets Gibbs free energy values.
    /// Returns normalized format with None key for single-phase systems.
    pub fn get_dG(&self) -> HashMap<Option<String>, HashMap<String, f64>> {
        match self {
            CustomSubstance::OnePhase(one_phase) => {
                let mut result = HashMap::new();
                result.insert(None, one_phase.dG.clone());
                result
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => phase_or_solution.dG.clone(),
        }
    }
}

/// Multi-phase thermodynamic system managing multiple phases or solutions.
/// 
/// Each phase is identified by an `Option<String>` key where `Some(name)` represents
/// named phases like "gas", "liquid", "solid", etc.
/// 
/// # Example
/// ```rust
/// let mut system = PhaseOrSolution::new();
/// // Add gas phase data
/// // system.subs_data.insert(Some("gas".to_string()), gas_substances);
/// ```
pub struct PhaseOrSolution {
    pub subs_data: HashMap<Option<String>, SubsData>,
    pub substance_phase_mapping: Option<SubstancePhaseMapping>,
    pub symbolic_vars: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    pub vec_of_n_vars: Vec<Expr>,
    pub Np_vec: Vec<Expr>,
    pub map_of_var_each_substance: HashMap<Option<String>, HashMap<String, Expr>>,
    /// hashmap of Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    pub dG: HashMap<Option<String>, HashMap<String, f64>>,
    /// hashmap of Gibbs free energy of a given mixure of substances (function at given T, P, concentration)
    pub dG_fun: HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,
    >,
    /// hashmap of Gibbs free energy of a given mixure of substances (symbolic result at given T, P, concentration)
    pub dG_sym: HashMap<Option<String>, HashMap<String, Expr>>,

    /// hashmap of entropy of a given mixure of substances (numerical result at given T, P, concentration)
    pub dS: HashMap<Option<String>, HashMap<String, f64>>,
    /// hashmap {phase or solution name:{ substance name: entropy function(T, vector moles of substances, totsl number of moles in phase or solution)} }
    pub dS_fun: HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    >,
    /// hashmap of entropy of a given mixure of substances (symbolic result at given T, P, concentration)
    pub dS_sym: HashMap<Option<String>, HashMap<String, Expr>>,
}

impl PhaseOrSolution {
    /// Creates a new empty multi-phase system.
    pub fn new() -> Self {
        Self {
            subs_data: HashMap::new(),
            substance_phase_mapping: None,
            symbolic_vars: HashMap::new(),
            vec_of_n_vars: Vec::new(),
            Np_vec: Vec::new(),
            map_of_var_each_substance: HashMap::new(),
            dG: HashMap::new(),
            dG_fun: HashMap::new(),
            dG_sym: HashMap::new(),
            dS: HashMap::new(),
            dS_fun: HashMap::new(),
            dS_sym: HashMap::new(),
        }
    }

    /// Sets the substance-phase mapping for the system.
    /// Used to track which substances belong to which phases.
    pub fn set_substance_phase_mapping(&mut self, mapping: SubstancePhaseMapping) {
        self.substance_phase_mapping = Some(mapping);
    }

    /// Applies a fallible function to all phases and collects results.
    /// Stops on first error and returns it.
    fn apply_to_all_phases<R>(&mut self, f: impl Fn(&mut SubsData) -> SubsDataResult<R>) -> SubsDataResult<HashMap<Option<String>, R>> {
        let mut result = HashMap::with_capacity(self.subs_data.len());
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let phase_result = f(subsdata)?;
            result.insert(phase_name.clone(), phase_result);
        }
        Ok(result)
    }

    /// Applies an infallible function to all phases and collects results.
    /// More efficient than the fallible version for operations that cannot fail.
    fn apply_to_all_phases_infallible<R>(&mut self, f: impl Fn(&mut SubsData) -> R) -> HashMap<Option<String>, R> {
        self.subs_data.iter_mut()
            .map(|(phase_name, subsdata)| (phase_name.clone(), f(subsdata)))
            .collect()
    }
}

impl Clone for PhaseOrSolution {
    fn clone(&self) -> Self {
        let mut new_dG_fun = HashMap::new();
        for (phase_or_solution, map_dG_fun) in &self.dG_fun {
            let mut inner_map = HashMap::new();
            for (substance, _) in map_dG_fun {
                let placeholder: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static> =
                    Box::new(|_: f64, _: Option<Vec<f64>>, _: Option<f64>| 0.0);
                inner_map.insert(substance.clone(), placeholder);
            }
            new_dG_fun.insert(phase_or_solution.clone(), inner_map);
        }

        let mut new_dS_fun = HashMap::new();
        for (phase_or_solution, map_dS_fun) in &self.dS_fun {
            let mut inner_map = HashMap::new();
            for (substance, _) in map_dS_fun {
                let placeholder: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static> =
                    Box::new(|_: f64, _: Option<Vec<f64>>, _: Option<f64>| 0.0);
                inner_map.insert(substance.clone(), placeholder);
            }
            new_dS_fun.insert(phase_or_solution.clone(), inner_map);
        }

        Self {
            subs_data: self.subs_data.clone(),
            substance_phase_mapping: self.substance_phase_mapping.clone(),
            symbolic_vars: self.symbolic_vars.clone(),
            vec_of_n_vars: self.vec_of_n_vars.clone(),
            Np_vec: self.Np_vec.clone(),
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

impl fmt::Debug for PhaseOrSolution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PhaseOrSolution")
            .field("subs_data", &self.subs_data)
            .field("substance_phase_mapping", &self.substance_phase_mapping)
            .field("symbolic_vars", &self.symbolic_vars)
            .field("vec_of_n_vars", &self.vec_of_n_vars)
            .field("Np_vec", &self.Np_vec)
            .field("map_of_var_each_substance", &self.map_of_var_each_substance)
            .field("dG", &self.dG)
            .field("dG_sym", &self.dG_sym)
            .field("dS", &self.dS)
            .field("dS_sym", &self.dS_sym)
            .finish()
    }
}
#[enum_dispatch]
pub trait ThermodynamicsCalculatorTrait {
    fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> SubsDataResult<()>;
    /// extract from libraries all thermal polynomials coefficents
    fn calculate_therm_map_of_properties(&mut self, temperature: f64) -> SubsDataResult<()>;
    fn calculate_therm_map_of_sym(&mut self) -> SubsDataResult<()>;

    fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>>;
    /// calculating Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    fn calcutate_Gibbs_free_energy(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>>;
    /// calculating Gibbs free energy of a given mixure of substances (symbolic result at given T, P, concentration)
    fn calculate_Gibbs_sym(&mut self, T: f64) -> SubsDataResult<()>;
    ///  calculating Gibbs free energy of a given mixure of substances (returns a function that calculates Gibbs free energy at given T, P, concentration)
    fn calculate_Gibbs_fun(&mut self, T: f64, P: f64);

    fn calculate_Gibbs_fun_unmut(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,
    >;
     fn set_P_to_sym_in_G_sym(&mut self, P:f64) ;

     fn set_T_to_sym_in_G_sym(&mut self, T:f64) ;
    //////////////////////////////////////////S - entropy///////////////////////////////////////
    ///  calculating entropy of a given mixure of substances (numerical result at given T, P, concentration)
    fn calculate_S(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> SubsDataResult<()>;
    ///  calculating entropy of a given mixure of substances (symbolic result at given T, P, concentration)
    fn calculate_S_sym(&mut self, T: f64) -> SubsDataResult<()>;
    ///  calculating entropy of a given mixure of substances (returns a function that calculates entropy at given T, P, concentration)
    fn calculate_S_fun(&mut self, T: f64, P: f64);
    fn calculate_S_fun_unmut(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    >;
    ///  set the system pressure, pressure unit, molar masses and mass unit
    fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()>;
    /// if not found in the library, go to NIST
    fn if_not_found_go_NIST(&mut self) -> SubsDataResult<()>;
    /// returns a SubstancesContainer
    fn extract_SubstancesContainer(&mut self) -> SubsDataResult<SubstancesContainer>;
    /// returns all substances in all phases
    fn get_all_substances(&mut self) -> Vec<String>;
    /// calculating element composition matrix and hashmap of molar masses
    fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> SubsDataResult<(DMatrix<f64>, HashMap<String, f64>, Vec<String>)>;
    fn indexed_moles_variables(
        &mut self,
    ) -> SubsDataResult<(
        HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
        Vec<Expr>,
        Vec<Expr>,
        HashMap<Option<String>, HashMap<String, Expr>>,
    )>;

    /////////////////////////////equations for G->min///////////////////////////////////////
    fn calculate_Lagrange_equations_sym(
        &mut self,
        A: DMatrix<f64>,

        Tm: f64,
    ) -> SubsDataResult<Vec<Expr>>;
    fn calculate_Lagrange_equations_fun(
        &mut self,
        A: DMatrix<f64>,
        dG_fun: HashMap<
            Option<String>,
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,
        >,
        Tm: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    >;
    fn calculate_Lagrange_equations_fun2(
        &mut self,
        A: DMatrix<f64>,
              T: f64,
        P: f64,
        Tm: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    >;
    //////////////////////////////equations for Keq/////////////////////////////////////////
    /*
    fn calculate_K_eq_sym(
        &mut self,
        stolich_matrix: DMatrix<f64>,

        map_of_indexed_vars: HashMap<String, (Expr, Vec<Expr>)>,
        G: HashMap<Option<String>, HashMap<String, Expr>>,
    );
    fn calculate_K_eq_fun(
        &mut self,
        stolich_matrix: DMatrix<f64>,
        G_fun: HashMap<
            Option<String>,
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
        >,
    );
    fn calculate_K_eq_equation_sym(
        &mut self,
        stolich_matrix: DMatrix<f64>,
        subs: Vec<String>,
        map_of_var_each_substance: HashMap<Option<(String, Expr)>, HashMap<String, Expr>>,
        G: HashMap<Option<String>, HashMap<String, Expr>>,
    );
    fn calculate_K_eq_equation_fun(
        &mut self,
        stolich_matrix: DMatrix<f64>,
        G_fun: HashMap<
            Option<String>,
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
        >,
    );
    */
    ////////////////////////////////////////////////////////////////////////////////////////
    fn create_full_map_of_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
        //  physiscal_nature_of_phase_or_solution: Option<HashMap<String, Phases>>,
    ) -> SubsDataResult<(
        HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
        HashMap<String, f64>,
    )>;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Implementation of the ThermodynamicsCalculatorTrait for PhaseOrSolution (multiole phases case)
impl ThermodynamicsCalculatorTrait for PhaseOrSolution {
    fn extract_all_thermal_coeffs(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.apply_to_all_phases(|subsdata| subsdata.extract_all_thermal_coeffs(temperature))?;
        Ok(())
    }
    
    fn calculate_therm_map_of_properties(&mut self, temperature: f64) -> SubsDataResult<()> {
        self.apply_to_all_phases(|subsdata| subsdata.calculate_therm_map_of_properties(temperature))?;
        Ok(())
    }
    
    fn calculate_therm_map_of_sym(&mut self) -> SubsDataResult<()> {
        self.apply_to_all_phases(|subsdata| subsdata.calculate_therm_map_of_sym())?;
        Ok(())
    }
    fn extract_coeffs_if_current_coeffs_not_valid(
        &mut self,
        temperature: f64,
    ) -> SubsDataResult<Vec<String>> {
        let mut results = Vec::new();
        for (_, subsdata) in self.subs_data.iter_mut() {
            let result =
                subsdata.extract_coeffs_if_current_coeffs_not_valid_for_all_subs(temperature)?;
            results.extend(result);
        }
        Ok(results)
    }
    fn calcutate_Gibbs_free_energy(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let n_Np = n
                .get(&phase_name)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "phase data".to_string(),
                    substance: phase_name.clone().unwrap_or("None".to_string()),
                })?
                .clone();
            let Np = n_Np.0;
            let n = n_Np.1;
            let map_to_insert_phase = subsdata.calc_dG_for_one_phase(P, T, n.clone(), Np);
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        self.dG = map_to_insert.clone();
        Ok(map_to_insert)
    }
    fn calculate_Gibbs_sym(&mut self, T: f64) -> SubsDataResult<()> {
        let n = self.symbolic_vars.clone();
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let n_Np = n
                .get(&phase_name)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "phase data".to_string(),
                    substance: phase_name.clone().unwrap_or("None".to_string()),
                })?
                .clone();
            let Np = n_Np.0;
            let n = n_Np.1;
            let map_to_insert_phase =
                subsdata.calculate_Gibbs_sym_one_phase(T, n.clone(), Np.clone());
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        self.dG_sym = map_to_insert;
        Ok(())
    }

     fn set_P_to_sym_in_G_sym(&mut self, P:f64) {
         let dG_sym = &mut self.dG_sym;
        for (_phase_or_solution_name, sym_fun) in dG_sym.iter_mut() {
            for (_, sym_fun) in sym_fun.iter_mut() {
                *sym_fun = sym_fun.set_variable("P", P).simplify()
            }
        }
    }

     fn set_T_to_sym_in_G_sym(&mut self, T:f64) {
              let dG_sym = &mut self.dG_sym;
        for (_phase_or_solution_name, sym_fun) in dG_sym.iter_mut() {
            for (_, sym_fun) in sym_fun.iter_mut() {
                *sym_fun = sym_fun.set_variable("T", T).simplify()
            }
        }
    }

    fn calculate_Gibbs_fun(&mut self, T: f64, P: f64) {
        self.dG_fun = self.apply_to_all_phases_infallible(|subsdata| {
            subsdata.calculate_Gibbs_fun_one_phase(P, T)
        });
    }

    fn calculate_Gibbs_fun_unmut(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,
    > {
        self.apply_to_all_phases_infallible(|subsdata| {
            subsdata.calculate_Gibbs_fun_one_phase(P, T)
        })
    }
    fn calculate_S(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> SubsDataResult<()> {
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let n_Np = n
                .get(&phase_name)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "phase data".to_string(),
                    substance: phase_name.clone().unwrap_or("None".to_string()),
                })?
                .clone();
            let Np = n_Np.0;
            let n = n_Np.1;
            let map_to_insert_phase = subsdata.calculate_S_for_one_phase(P, T, n.clone(), Np);
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        self.dS = map_to_insert.clone();
        Ok(())
    }
    fn calculate_S_sym(&mut self, T: f64) -> SubsDataResult<()> {
        let n = self.symbolic_vars.clone();
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let n_Np = n
                .get(&phase_name)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "phase data".to_string(),
                    substance: phase_name.clone().unwrap_or("None".to_string()),
                })?
                .clone();
            let Np = n_Np.0;
            let n = n_Np.1;
            let map_to_insert_phase =
                subsdata.calculate_S_sym_for_one_phase(T, n.clone(), Np.clone());
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        self.dS_sym = map_to_insert.clone();
        Ok(())
    }

    fn calculate_S_fun(&mut self, T: f64, P: f64) {
        self.dS_fun = self.apply_to_all_phases_infallible(|subsdata| {
            subsdata.calculate_S_fun_for_one_phase(P, T)
        });
    }
    
    fn calculate_S_fun_unmut(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    > {
        self.apply_to_all_phases_infallible(|subsdata| {
            subsdata.calculate_S_fun_for_one_phase(P, T)
        })
    }
    fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> SubsDataResult<()> {
        self.apply_to_all_phases(|subs_data| {
            subs_data.set_P(pressure, pressure_unit.clone());
            subs_data.set_M(molar_masses.clone(), mass_unit.clone());
            Ok(())
        })?;
        Ok(())
    }
    
    fn if_not_found_go_NIST(&mut self) -> SubsDataResult<()> {
        self.apply_to_all_phases(|subs_data| subs_data.if_not_found_go_NIST())?;
        Ok(())
    }
    fn extract_SubstancesContainer(&mut self) -> SubsDataResult<SubstancesContainer> {
        let mut map_of_subs = HashMap::new();
        for (phase_name, subs_data) in &self.subs_data {
            let substances = subs_data.substances.clone();
            let phase_key = phase_name
                .clone()
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "phase name".to_string(),
                    substance: "unknown".to_string(),
                })?;
            map_of_subs.insert(phase_key, substances);
        }
        Ok(SubstancesContainer::MultiPhase(map_of_subs))
    }
    fn get_all_substances(&mut self) -> Vec<String> {
        let mut all_substances = HashSet::new();
        for (_, subs_data) in &mut self.subs_data {
            let substances = subs_data.substances.clone();
            all_substances.extend(substances);
        }
        let mut all_substances: Vec<String> = all_substances.into_iter().collect();
        all_substances.sort();
        all_substances
    }
    fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> SubsDataResult<(DMatrix<f64>, HashMap<String, f64>, Vec<String>)> {
        let mut hashmap_of_compositions: HashMap<String, HashMap<String, f64>> = HashMap::new();
        let mut hashset_of_elems_: HashSet<String> = HashSet::new();
        let mut hasmap_of_molar_mass_: HashMap<String, f64> = HashMap::new();
        for (_phase_name, subsdata) in self.subs_data.iter_mut() {
            let substances = &subsdata.substances.to_owned();
            let (hasmap_of_molar_mass, vec_of_compositions, hashset_of_elems) =
                SubsData::calculate_elem_composition_and_molar_mass_local(subsdata, groups.clone())
                    .map_err(|e| SubsDataError::MatrixOperationFailed(e.to_string()))?;
            // forming hashmap of compositions of substancess
            let map: HashMap<String, HashMap<String, f64>> = substances
                .iter()
                .zip(vec_of_compositions.iter())
                .map(|(substance, composition)| (substance.clone(), composition.clone()))
                .collect();
            hashmap_of_compositions.extend(map);
            hashset_of_elems_.extend(hashset_of_elems);
            hasmap_of_molar_mass_.extend(hasmap_of_molar_mass);
        }
        let mut unique_vec_of_elems = hashset_of_elems_.into_iter().collect::<Vec<_>>();
        unique_vec_of_elems.sort();
        let num_rows = unique_vec_of_elems.len();
        let substances = self.get_all_substances();
        let num_cols = substances.len();
        // allocate matrix with num of rows = num of elements and num of cols = num of substances

        let mut matrix: DMatrix<f64> = DMatrix::zeros(num_rows, num_cols);
        for (i, substance_i) in substances.iter().enumerate() {
            for j in 0..unique_vec_of_elems.len() {
                let element_j = unique_vec_of_elems[j].clone();
                if let Some(count) = hashmap_of_compositions
                    .get(&substance_i.clone())
                    .unwrap()
                    .get(&element_j)
                {
                    matrix[(j, i)] += *count as f64;
                }
            }
        }
        Ok((
            matrix.transpose(),
            hasmap_of_molar_mass_,
            unique_vec_of_elems,
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
        let mut vec_of_n_vars: Vec<Expr> = Vec::new();
        let mut Np_vec = Vec::new();
        let mut map_of_indexed_vars = HashMap::new();
        let mut map_of_var_each_substance: HashMap<Option<String>, HashMap<String, Expr>> =
            HashMap::new();
        for (i, (phase_name, subdata)) in self.subs_data.iter().enumerate() {
            let Np = Expr::Var(format!("Np{}", i));
            Np_vec.push(Np.clone());
            map_of_var_each_substance.insert(phase_name.clone(), HashMap::new());
            let mut v = Vec::new();
            for (j, _subdata) in subdata.substances.iter().enumerate() {
                let n = Expr::Var(format!("n{}", i + j));
                vec_of_n_vars.push(n.clone());
                v.push(n.clone());
                map_of_var_each_substance
                    .get_mut(phase_name)
                    .unwrap()
                    .insert(_subdata.clone(), n);
            }
            map_of_indexed_vars.insert(phase_name.clone(), (Some(Np), Some(v)));
        }
        self.symbolic_vars = map_of_indexed_vars.clone();
        self.vec_of_n_vars = self.vec_of_n_vars.clone();
        self.Np_vec = Np_vec.clone();
        self.map_of_var_each_substance = map_of_var_each_substance.clone();
        Ok((
            map_of_indexed_vars,
            vec_of_n_vars,
            Np_vec,
            map_of_var_each_substance,
        ))
    }
    /////////////////////////////equations for G->min///////////////////////////////////////
    fn calculate_Lagrange_equations_sym(
        &mut self,
        A: DMatrix<f64>,

        Tm: f64,
    ) -> SubsDataResult<Vec<Expr>> {
        let G = self.dG_sym.clone();
        let A = A.clone().transpose(); // now rows are elements and columns are substances
        let n_elements = A.nrows();

        // create a vector of Lagrange multipliers for each element
        let Lambda = Expr::IndexedVars(n_elements, "Lambda").0;
        let n_substances = A.ncols();
        // number of key values pairs in the nested hashmaps = number of substances
        let G_len: usize = G.values().map(|inner_map| inner_map.len()).sum();
        if n_substances != G_len {
            println!("A {} \n  G = {:?} \n, Lambda = {:?}", A, G, Lambda);
            return Err(SubsDataError::CalculationFailed {
                substance: "multiple".to_string(),
                operation: "Lagrange equations".to_string(),
                reason: format!(
                    "The number of substances ({}) does not match the number of mole variables ({})",
                    n_substances, G_len
                ),
            });
        }
        let substances = self.get_all_substances();
        let mut vec_of_eqs: Vec<Expr> = Vec::new();
        for (_phasse_or_solution, G_phase) in G.iter() {
            for i in 0..n_substances {
                let Tm = Expr::Const(Tm);
                let substance_i = &substances[i];
                if let Some(G_i) = G_phase.get(substance_i) {
                    // get vector of numbers of elements corresponding to the i-th substance
                    let col_of_subs_i = A
                        .column(i)
                        .iter()
                        .map(|&v| Expr::Const(v))
                        .collect::<Vec<Expr>>();
                    // sum over all elements of the i-th substance multiplied by the corresponding Lagrange multiplier
                    let sum_by_elemnts: Expr = col_of_subs_i
                        .iter()
                        .zip(Lambda.iter())
                        .map(|(aij, Lambda_j)| aij.clone() * Lambda_j.clone())
                        .fold(Expr::Const(0.0), |acc, x| acc + x)
                        .simplify();
                    let eq_i = sum_by_elemnts + (G_i.clone() / (R_sym * Tm)).simplify();
                    vec_of_eqs.push(eq_i.simplify());
                }
            }
        }
        Ok(vec_of_eqs)
    }

    fn calculate_Lagrange_equations_fun(
        &mut self,
        A: DMatrix<f64>,
        G_fun: HashMap<
            Option<String>,
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,
        >,

        Tm: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    > {
        let substances = self.clone().get_all_substances();
        let f = move |T: f64, n: Option<Vec<f64>>, Np: Option<f64>, Lambda: Vec<f64>| -> Vec<f64> {
            let A = A.transpose(); // now rows are elements and columns are substances

            // create a vector of Lagrange multipliers for each element
            let n_substances = A.ncols();
            let Lambda: DVector<f64> = DVector::from_vec(Lambda);

            let mut vec_of_eqs: Vec<f64> = Vec::new();
            for (_phasse_or_solution, G_phase) in G_fun.iter() {
                for i in 0..n_substances {
                    let substance_i = &substances[i];
                    if let Some(G_i) = G_phase.get(substance_i) {
                        let col_i: Vec<f64> = A.column(i).iter().map(|&v| v).collect();
                        // get vector of numbers of elements corresponding to the i-th substance
                        let col_of_subs_i: DVector<f64> = DVector::from(col_i);
                        // sum over all elements of the i-th substance multiplied by the corresponding Lagrange multiplier
                        let sum_by_elemnts: f64 = col_of_subs_i.dot(&Lambda);
                        let eq_i = sum_by_elemnts + G_i(T, n.clone(), Np) / (R * Tm);
                        vec_of_eqs.push(eq_i);
                    } // if let Some(G_i)
                } // for i
            } // for 
            vec_of_eqs
        };
        Ok(Box::new(f))
    }

        fn calculate_Lagrange_equations_fun2(
        &mut self,
        A: DMatrix<f64>,
        T: f64,
        P: f64,
        Tm: f64,
    ) -> SubsDataResult<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    >{

            let dG_fun = self.calculate_Gibbs_fun_unmut(T, P);
             let Lagrange_equations: Box< dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static> = self
            .calculate_Lagrange_equations_fun(A, dG_fun, Tm)
            .unwrap();
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
        let s = stolich_matrix.nrows();
        let all_substances = self.get_all_substances();
        assert_eq!(all_substances.len(), s);
        let T = Expr::Var("T".to_string());

        for (phase, subdata) in self.subs_data{
            let substance = subdata.substances;
            let subs_data = subdata.calculate_Gibbs_sym_one_phase(298.0, None, None);

        }

        let mut K_vec: Vec<Expr> = Vec::with_capacity(r);
        for j in 0..r {
            let mut dG_reaction_j = Expr::Const(0.0);
            for (i, subs_i) in all_substances.iter().enumerate() {

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
    }
    fn calculate_K_eq_equation_sym(
        &mut self,
        stolich_matrix: DMatrix<f64>,
        subs: Vec<String>,
        map_of_var_each_substance: HashMap<Option<String>, HashMap<String, Expr>>,
        G: HashMap<Option<String>, HashMap<String, Expr>>,
    ) {
        let m = subs.len();
        let r = stolich_matrix.ncols();
        let mut f = vec![0.0; r];

        // for each reaction k
        for k in 0..r {
            let mut sum_ln_n = 0.0;
            let mut dg0 = 0.0;

            for i in 0..m {
                sum_ln_n += stolich_matrix[(i, k)] * n[i].ln();
                dg0 += reactions[(i, k)] * gibbs[i](temperature);
            }

            let mut sum_ln_n_phase = Expr::Const(0.0);
            let mut sum_phi = 0.0;
            for (phase_name, subs_data) in &self.subs_data {
                let phase_name = phase_name.unwrap();
                let (Npj, vec_of_nij) = map_of_indexed_vars.get(&phase_name).unwrap();
                let dNk = delta_n.get(&phase_name)[r];
                if dNk != 0.0 {
                    let G_vec = subs_data.calculate_Gibbs_sym_one_phase(T, None, None);
                    sum_ln_n_phase += Expr::Const(dNk) * Np.ln();
                }
            }
            for j in 0..phases.len() {
                let dnk = delta_n[k][j];
                if dnk != 0.0 {
                    sum_ln_n_phase += dnk * n_phase[j].ln();
                }
            }

            let ln_k = -dg0 / (r_gas * temperature);

            f[k] = sum_ln_n - sum_ln_n_phase + sum_phi - ln_k;
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
    }
    */
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// This function, create_full_map_of_mole_numbers, takes a hashmap of non-zero mole numbers and returns a tuple of three hashmaps.
    /// The first hashmap is the original input with missing substances inserted with a value of 0.0.
    /// The second hashmap contains the same information as the first, but with the inner hashmap values converted to a vector of floats.
    /// The third hashmap contains the sum of the values of the inner hashmap for each substance across all phases.
    fn create_full_map_of_mole_numbers(
        &self,
        non_zero_number_of_moles: HashMap<
            Option<String>,
            (Option<f64>, Option<HashMap<String, f64>>),
        >,
        //  physiscal_nature_of_phase_or_solution: Option<HashMap<String, Phases>>,
    ) -> SubsDataResult<(
        HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
        HashMap<String, f64>,
    )> {
        //  hashmap  HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>, we iterate through it and when some string from vector Vec<String> is not
        // found as the key in the inner map it is inserted in the inner map with value 0.0
        let mut full_map_of_mole_numbers = non_zero_number_of_moles.clone();

        let mut map_of_mole_num_vecs: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)> =
            HashMap::new();

        //hashmap  HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)> keys of the inner may can be the same.
        //let's generate map HashMap<String, f64> where f64 is the sum of the values of inner map with same keys
        let mut initial_map_of_mole_numbers: HashMap<String, f64> = HashMap::new();
        let subs_cloned: HashMap<Option<String>, SubsData> = self.subs_data.clone();
        for (phase, value) in full_map_of_mole_numbers.iter_mut() {
            if let (Np, Some(inner_map)) = (value.0.clone(), &mut value.1) {
                //  println!("phase: {}", phase.clone().is_none());
                //  println!("subdata: {:?}", subs_cloned);
                let subsdata = subs_cloned.get(phase).unwrap();
                let subs = subsdata.clone().substances;
                // For each required key, insert with 0.0 if not present
                for key in &subs {
                    inner_map.entry(key.clone()).or_insert(0.0);
                }
                // making
                let inner_map = inner_map.clone();
                let vec_of_mole_nums = inner_map.values().cloned().collect::<Vec<f64>>();
                map_of_mole_num_vecs.insert(phase.clone(), (Np, Some(vec_of_mole_nums)));
                // making initial map of mole numbers
                for (keys, values) in inner_map.iter().clone() {
                    *initial_map_of_mole_numbers
                        .entry(keys.clone())
                        .or_insert(0.0) += values;
                }
            }
        }

        Ok((
            full_map_of_mole_numbers,
            map_of_mole_num_vecs,
            initial_map_of_mole_numbers,
        ))
    }
}

////////////////////FACTORY METHODS////////////////////
pub struct SubstanceSystemFactory;
//use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::Thermodynamics;
impl SubstanceSystemFactory {
    /// Creates a substance system based on the provided container type
    ///
    /// # Arguments
    /// * `container` - The container specifying substances and phases
    /// * `library_priorities` - Libraries to prioritize for data lookup
    /// * `permitted_libraries` - Libraries that are permitted but not prioritized
    ///
    /// # Returns
    /// A configured CustomSubstance ready for thermodynamic calculations
    pub fn create_system(
        container: SubstancesContainer,
        physiscal_nature_of_phase_or_solution: Option<HashMap<String, Phases>>,
        library_priorities: Vec<String>,
        permitted_libraries: Vec<String>,
        explicit_search_insructions: Option<HashMap<String, String>>,
        search_in_NIST: bool,
    ) -> Result<CustomSubstance, String> {
        match container {
            SubstancesContainer::SinglePhase(substances) => {
                // Create a single phase system
                let mut subs_data = SubsData::new();
                subs_data.substances = substances;

                // Configure library priorities
                subs_data
                    .set_multiple_library_priorities(library_priorities, LibraryPriority::Priority);
                subs_data.set_multiple_library_priorities(
                    permitted_libraries,
                    LibraryPriority::Permitted,
                );
                if let Some(explicit_search_instructions) = explicit_search_insructions.clone() {
                    subs_data.set_explicis_searh_instructions(explicit_search_instructions);
                };
                // Set default phases to Gas
                for substance in &subs_data.substances {
                    subs_data
                        .map_of_phases
                        .insert(substance.clone(), Some(Phases::Gas));
                }

                // Search for substance data
                let _ = subs_data.search_substances();
                // if not found go to NIST
                if search_in_NIST {
                    let _ = subs_data.if_not_found_go_NIST();
                }

                // Extract thermal coefficients at standard temperature
                //  let _ = subs_data.extract_all_thermal_coeffs(298.15);
                let _ = subs_data.parse_all_thermal_coeffs();
                let mut op = OnePhase::new();
                op.subs_data = subs_data;
                Ok(CustomSubstance::OnePhase(op))
            }
            SubstancesContainer::MultiPhase(phase_substances) => {
                // Create a multi-phase system
                let mut phase_or_solution = PhaseOrSolution::new();

                // Create and configure SubsData for each phase
                for (phase_name, substances) in phase_substances {
                    let mut phase_data = SubsData::new();
                    phase_data.substances = substances;

                    // Configure library priorities
                    phase_data.set_multiple_library_priorities(
                        library_priorities.clone(),
                        LibraryPriority::Priority,
                    );
                    phase_data.set_multiple_library_priorities(
                        permitted_libraries.clone(),
                        LibraryPriority::Permitted,
                    );
                    /*
                    // Set default phases based on phase name or default to Gas
                    let default_phase = match phase_name.to_lowercase().as_str() {
                        "gas" | "vapor" => Phases::Gas,
                        "liquid" => Phases::Liquid,
                        "solid" => Phases::Solid,
                        _ => Phases::Gas,
                    };
                    */
                    let default_phase = match physiscal_nature_of_phase_or_solution {
                        Some(ref map_of_phases) => {
                            let phase = map_of_phases.get(&phase_name);
                            phase.unwrap_or(&Phases::Gas)
                        }
                        None => &Phases::Gas,
                    };

                    for substance in &phase_data.substances {
                        phase_data
                            .map_of_phases
                            .insert(substance.clone(), Some(default_phase.clone()));
                    }
                    if let Some(ref explicit_search_instructions) = explicit_search_insructions {
                        phase_data
                            .set_explicis_searh_instructions(explicit_search_instructions.clone());
                    };
                    // Search for substance data
                    let _ = phase_data.search_substances();
                    // if not found go to NIST
                    if search_in_NIST {
                        let _ = phase_data.if_not_found_go_NIST();
                    }
                    // Extract thermal coefficients at standard temperature
                    //  let _ = phase_data.extract_all_thermal_coeffs(298.15);

                    // Add to phase collection

                    //  println!("creating phase or solution with phase_name: {} and phase_data {:?}", phase_name, phase_data.clone());
                    let _ = phase_data.parse_all_thermal_coeffs();
                    phase_or_solution
                        .subs_data
                        .insert(Some(phase_name), phase_data);
                }

                Ok(CustomSubstance::PhaseOrSolution(phase_or_solution))
            }
        } //end match
    } //end fn

    /*
    fn create_system2(
        container: SubstancesContainer,
        library_priorities: Vec<String>,
        permitted_libraries: Vec<String>,
        explicit_search_insructions: Option<HashMap<String, String>>,
        search_in_NIST: bool,
    ) -> Result<Thermodynamics, String> {
        let mut system = Thermodynamics::new();
        let CustomSubstance_instance = Self::create_system(
            container,
            library_priorities,
            permitted_libraries,
            explicit_search_insructions,
            search_in_NIST,
        )?;

    }
     */
}
