use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
use crate::Thermodynamics::dG_dS::{
    calc_dG_for_one_phase, calculate_Gibbs_fun_one_phase, calculate_Gibbs_sym_one_phase,
    calculate_S_for_one_phase, calculate_S_fun_for_one_phase, calculate_S_sym_for_one_phase,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use enum_dispatch::enum_dispatch;
use nalgebra::{DMatrix, DVector};

use std::collections::{HashMap, HashSet};
use std::f64;

/// One phase or multiple phases or solutions
#[derive(Debug, Clone)]
#[enum_dispatch(ThermodynamicsCalculatorTrait)]
pub enum CustomSubstance {
    OnePhase(SubsData),
    PhaseOrSolution(PhaseOrSolution),
}
/// container of substances names in one phase or in multiple phases
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubstancesContainer {
    ///names of substances in one phase
    SinglePhase(Vec<String>),
    /// map of phases and their substances
    MultiPhase(HashMap<String, Vec<String>>),
}

impl SubstancesContainer {
    // Helper method to get all substances as a flat vector
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
///  structure contains map: phase name - substances data in this phase
#[derive(Debug, Clone)]
pub struct PhaseOrSolution {
    pub subs_data: HashMap<Option<String>, SubsData>,
}

impl PhaseOrSolution {
    pub fn new() -> Self {
        Self {
            subs_data: HashMap::new(),
        }
    }
}
#[enum_dispatch]
pub trait ThermodynamicsCalculatorTrait {
    /// calculating Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    fn calcutate_Gibbs_free_energy(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, f64>>, String>;
    /// calculating Gibbs free energy of a given mixure of substances (symbolic result at given T, P, concentration)
    fn calculate_Gibbs_sym(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, Expr>>, String>;
    ///  calculating Gibbs free energy of a given mixure of substances (returns a function that calculates Gibbs free energy at given T, P, concentration)
    fn calculate_Gibbs_fun(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    >;
    ///  calculating entropy of a given mixure of substances (numerical result at given T, P, concentration)
    fn calculate_S(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, f64>>, String>;
    ///  calculating entropy of a given mixure of substances (symbolic result at given T, P, concentration)
    fn calculate_S_sym(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, Expr>>, String>;
    ///  calculating entropy of a given mixure of substances (returns a function that calculates entropy at given T, P, concentration)
    fn calculate_S_fun(
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
    ) -> Result<(), String>;
    /// if not found in the library, go to NIST
    fn if_not_found_go_NIST(&mut self) -> Result<(), String>;
    /// returns a SubstancesContainer
    fn extract_SubstancesContainer(&mut self) -> Result<SubstancesContainer, String>;
    /// returns all substances in all phases
    fn get_all_substances(&mut self) -> Vec<String>;
    /// calculating element composition matrix and hashmap of molar masses
    fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> Result<(DMatrix<f64>, HashMap<String, f64>, Vec<String>), String>;
    fn indexed_moles_variables(
        &mut self,
    ) -> Result<
        (
            HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
            Vec<Expr>,
            Vec<Expr>,
            HashMap<Option<String>, HashMap<String, Expr>>,
        ),
        String,
    >;
    fn calculate_Lagrange_equations_sym(
        &mut self,
        A: DMatrix<f64>,
        G: HashMap<Option<String>, HashMap<String, Expr>>,
    ) -> Result<Vec<Expr>, String>;
    fn calculate_Lagrange_equations_fun(
        &mut self,
        A: DMatrix<f64>,
        G_fun: HashMap<
            Option<String>,
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
        >,
    ) -> Result<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
        String,
    >;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Implementation of the ThermodynamicsCalculatorTrait for PhaseOrSolution (multiole phases case)
impl ThermodynamicsCalculatorTrait for PhaseOrSolution {
    fn calcutate_Gibbs_free_energy(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, f64>>, String> {
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let n_Np = n.get(&phase_name).unwrap().clone();
            let Np = n_Np.0;
            let n = n_Np.1;
            let map_to_insert_phase = calc_dG_for_one_phase(P, subsdata, T, n.clone(), Np);
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        let dG = map_to_insert;
        Ok(dG)
    }
    fn calculate_Gibbs_sym(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, Expr>>, String> {
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let n_Np = n.get(&phase_name).unwrap().clone();
            let Np = n_Np.0;
            let n = n_Np.1;
            let map_to_insert_phase =
                calculate_Gibbs_sym_one_phase(subsdata, T, n.clone(), Np.clone());
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        Ok(map_to_insert)
    }
    fn calculate_Gibbs_fun(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    > {
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let map_to_insert_phase = calculate_Gibbs_fun_one_phase(subsdata, P, T);
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        map_to_insert
    }
    fn calculate_S(
        &mut self,
        T: f64,
        P: f64,

        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, f64>>, String> {
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let n_Np = n.get(&phase_name).unwrap().clone();
            let Np = n_Np.0;
            let n = n_Np.1;
            let map_to_insert_phase = calculate_S_for_one_phase(P, subsdata, T, n.clone(), Np);
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        let dG = map_to_insert;
        Ok(dG)
    }
    fn calculate_S_sym(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, Expr>>, String> {
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let n_Np = n.get(&phase_name).unwrap().clone();
            let Np = n_Np.0;
            let n = n_Np.1;
            let map_to_insert_phase =
                calculate_S_sym_for_one_phase(subsdata, T, n.clone(), Np.clone());
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        Ok(map_to_insert)
    }
    fn calculate_S_fun(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    > {
        let mut map_to_insert = HashMap::new();
        for (phase_name, subsdata) in self.subs_data.iter_mut() {
            let map_to_insert_phase = calculate_S_fun_for_one_phase(subsdata, P, T);
            map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
        }
        map_to_insert
    }
    fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> Result<(), String> {
        // Set pressure and molar masses for each phase
        for (_, subs_data) in &mut self.subs_data {
            subs_data.set_P(pressure, pressure_unit.clone());
            subs_data.set_M(molar_masses.clone(), mass_unit.clone());
        }
        Ok(())
    }
    fn if_not_found_go_NIST(&mut self) -> Result<(), String> {
        for (_, subs_data) in &mut self.subs_data {
            subs_data.if_not_found_go_NIST()?;
        }
        Ok(())
    }
    fn extract_SubstancesContainer(&mut self) -> Result<SubstancesContainer, String> {
        let mut map_of_subs = HashMap::new();
        for (phase_name, subs_data) in &self.subs_data {
            let substances = subs_data.substances.clone();
            map_of_subs.insert(phase_name.clone().unwrap(), substances);
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
    ) -> Result<(DMatrix<f64>, HashMap<String, f64>, Vec<String>), String> {
        let mut hashmap_of_compositions: HashMap<String, HashMap<String, f64>> = HashMap::new();
        let mut hashset_of_elems_: HashSet<String> = HashSet::new();
        let mut hasmap_of_molar_mass_: HashMap<String, f64> = HashMap::new();
        for (_phase_name, subsdata) in self.subs_data.iter_mut() {
            let substances = &subsdata.substances.to_owned();
            let (hasmap_of_molar_mass, vec_of_compositions, hashset_of_elems) =
                SubsData::calculate_elem_composition_and_molar_mass_local(subsdata, groups.clone())
                    .unwrap();
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
        String,
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
        Ok((
            map_of_indexed_vars,
            vec_of_n_vars,
            Np_vec,
            map_of_var_each_substance,
        ))
    }
    fn calculate_Lagrange_equations_sym(
        &mut self,
        A: DMatrix<f64>, // hasmap {phase or solution, {substance, expression}}
        G: HashMap<Option<String>, HashMap<String, Expr>>,
    ) -> Result<Vec<Expr>, String> {
        let A = A.transpose(); // now rows are elements and columns are substances
        let n_elements = A.nrows();
        // create a vector of Lagrange multipliers for each element
        let Lambda = Expr::IndexedVars(n_elements, "Lambda").0;
        let n_substances = A.ncols();

        if n_substances != G.len() {
            return Err(format!(
                "The number of substances ({}) does not match the number of mole variables ({})",
                n_substances,
                G.len()
            ));
        }
        let substances = self.get_all_substances();
        let mut vec_of_eqs: Vec<Expr> = Vec::new();
        for (_phasse_or_solution, G_phase) in G.iter() {
            for i in 0..n_substances {
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
                        .symplify();
                    let eq_i = sum_by_elemnts - G_i.clone();
                    vec_of_eqs.push(eq_i);
                }
            }
        }
        Ok(vec_of_eqs)
    }

    fn calculate_Lagrange_equations_fun(
        &mut self,
        A: DMatrix<f64>,
        G_fun: HashMap<
            Option<String>, //T, number of moles of each substance, total number of moles in phase
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
        >,
    ) -> Result<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
        String,
    > {
        let substances = self.get_all_substances();
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
                        let eq_i = sum_by_elemnts - G_i(T, n.clone(), Np);
                        vec_of_eqs.push(eq_i);
                    } // if let Some(G_i)
                } // for i
            } // for 
            vec_of_eqs
        };
        Ok(Box::new(f))
    }
}
/// Implementation of the ThermodynamicsCalculatorTrait for SubsData (single phase case)
impl ThermodynamicsCalculatorTrait for SubsData {
    fn calcutate_Gibbs_free_energy(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, f64>>, String> {
        let mut map_to_insert = HashMap::new();
        let n_Np = n.get(&None).unwrap().clone();
        let Np = n_Np.0;
        let n = n_Np.1;
        let map_to_insert_phase = calc_dG_for_one_phase(P, self, T, n, Np);
        map_to_insert.insert(None, map_to_insert_phase);
        let dG = map_to_insert;
        Ok(dG)
    }
    fn calculate_Gibbs_sym(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, Expr>>, String> {
        let mut map_to_insert = HashMap::new();
        let n_Np = n.get(&None).unwrap().clone();
        let Np = n_Np.0;
        let n = n_Np.1;
        let map_to_insert_phase = calculate_Gibbs_sym_one_phase(self, T, n, Np);
        map_to_insert.insert(None, map_to_insert_phase);
        Ok(map_to_insert)
    }
    fn calculate_Gibbs_fun(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    > {
        let mut map_to_insert = HashMap::new();
        let map_to_insert_phase = calculate_Gibbs_fun_one_phase(self, P, T);
        map_to_insert.insert(None, map_to_insert_phase);
        map_to_insert
    }
    fn calculate_S(
        &mut self,
        T: f64,
        P: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, f64>>, String> {
        let n_Np = n.get(&None).unwrap().clone();
        let Np = n_Np.0;
        let n = n_Np.1;
        let mut map_to_insert = HashMap::new();
        let map_to_insert_phase = calculate_S_for_one_phase(P, self, T, n, Np);
        map_to_insert.insert(None, map_to_insert_phase);
        let dG = map_to_insert;
        Ok(dG)
    }
    fn calculate_S_sym(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    ) -> Result<HashMap<Option<String>, HashMap<String, Expr>>, String> {
        let mut map_to_insert = HashMap::new();
        let n_Np = n.get(&None).unwrap().clone();
        let Np = n_Np.0;
        let n = n_Np.1;
        let map_to_insert_phase = calculate_S_sym_for_one_phase(self, T, n, Np);
        map_to_insert.insert(None, map_to_insert_phase);
        Ok(map_to_insert)
    }
    fn calculate_S_fun(
        &mut self,
        T: f64,
        P: f64,
    ) -> HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    > {
        let mut map_to_insert = HashMap::new();
        let map_to_insert_phase = calculate_S_fun_for_one_phase(self, P, T);
        map_to_insert.insert(None, map_to_insert_phase);
        map_to_insert
    }
    fn configure_system_properties(
        &mut self,
        pressure: f64,
        pressure_unit: Option<String>,
        molar_masses: HashMap<String, f64>,
        mass_unit: Option<String>,
    ) -> Result<(), String> {
        // Set pressure and molar masses for each phase
        self.set_P(pressure, pressure_unit.clone());
        self.set_M(molar_masses.clone(), mass_unit.clone());
        Ok(())
    }
    fn if_not_found_go_NIST(&mut self) -> Result<(), String> {
        self.if_not_found_go_NIST()?;
        Ok(())
    }
    fn extract_SubstancesContainer(&mut self) -> Result<SubstancesContainer, String> {
        let substances = self.substances.clone();
        Ok(SubstancesContainer::SinglePhase(substances))
    }
    fn get_all_substances(&mut self) -> Vec<String> {
        self.substances.clone()
    }
    fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> Result<(DMatrix<f64>, HashMap<String, f64>, Vec<String>), String> {
        self.calculate_elem_composition_and_molar_mass(groups)?;
        let elem_composition_matrix = self.elem_composition_matrix.clone().unwrap();

        let unique_elements = self.unique_elements.clone();
        let hashmap_of_molar_masses = self.hasmap_of_molar_mass.clone();
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
        String,
    > {
        // create a vector of indexed variables for the number of substances
        let n = Expr::IndexedVars(self.substances.len(), "N").0;
        // only one variable for the total number of moles
        let Np = Expr::Var("Np".to_string());
        let mut map_of_var_each_substance: HashMap<Option<String>, HashMap<String, Expr>> =
            HashMap::new();
        map_of_var_each_substance.insert(None, HashMap::new());
        for i in 0..self.substances.len() {
            let var = Expr::Var(format!("N{}", i));
            map_of_var_each_substance
                .get_mut(&None)
                .unwrap()
                .insert(self.substances[i].clone(), var);
        }
        let mut map_of_indexed_vars: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)> =
            HashMap::new();
        map_of_indexed_vars.insert(None, (Some(Np.clone()), Some(n.clone())));
        Ok((map_of_indexed_vars, n, vec![Np], map_of_var_each_substance))
    }
    fn calculate_Lagrange_equations_sym(
        &mut self,
        A: DMatrix<f64>,
        G: HashMap<Option<String>, HashMap<String, Expr>>,
    ) -> Result<Vec<Expr>, String> {
        let A = A.transpose(); // now rows are elements and columns are substances
        let n_elements = A.nrows();
        // create a vector of Lagrange multipliers for each element
        let Lambda = Expr::IndexedVars(n_elements, "Lambda").0;
        let n_substances = A.ncols();
        let G = G.get(&None).unwrap().clone();
        if n_substances != G.len() {
            return Err(format!(
                "The number of substances ({}) does not match the number of mole variables ({})",
                n_substances,
                G.len()
            ));
        }
        let subs = self.substances.clone();
        let mut vec_of_eqs: Vec<Expr> = Vec::new();
        for i in 0..n_substances {
            let subst_i = subs[i].clone();
            let G_i = G.get(&subst_i).unwrap().clone();
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
                .symplify();
            let eq_i = sum_by_elemnts - G_i.clone();
            vec_of_eqs.push(eq_i.symplify());
        }
        Ok(vec_of_eqs)
    }

    fn calculate_Lagrange_equations_fun(
        &mut self,
        A: DMatrix<f64>,
        G_fun: HashMap<
            Option<String>,
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
        >,
    ) -> Result<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
        String,
    > {
        let subs = self.substances.clone();
        let fun =
            move |T: f64, n: Option<Vec<f64>>, Np: Option<f64>, Lambda: Vec<f64>| -> Vec<f64> {
                let A = A.transpose(); // now rows are elements and columns are substances

                let n_substances = A.ncols();
                let G = G_fun.get(&None).unwrap();
                if n_substances != G.len() {
                    panic!("The number of substances does not match the number of mole variables");
                }

                let mut vec_of_eqs: Vec<f64> = Vec::new();
                for i in 0..n_substances {
                    let subst_i = subs[i].clone();
                    let G_i = G.get(&subst_i).unwrap();
                    // get vector of numbers of elements corresponding to the i-th substance
                    let col_of_subs_i = A.column(i).iter().map(|&v| v).collect::<Vec<f64>>();
                    // sum over all elements of the i-th substance multiplied by the corresponding Lagrange multiplier
                    let sum_by_elemnts: f64 = col_of_subs_i
                        .iter()
                        .zip(Lambda.iter())
                        .map(|(aij, Lambda_j)| aij * Lambda_j)
                        .fold(0.0, |acc, x| acc + x);
                    let eq_i = sum_by_elemnts - G_i(T, n.clone(), Np.clone());
                    vec_of_eqs.push(eq_i);
                }
                vec_of_eqs
            };
        Ok(Box::new(fun))
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
                subs_data.search_substances();
                // if not found go to NIST
                if search_in_NIST {
                    subs_data.if_not_found_go_NIST()?;
                }

                // Extract thermal coefficients at standard temperature
                //  let _ = subs_data.extract_all_thermal_coeffs(298.15);

                Ok(CustomSubstance::OnePhase(subs_data))
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

                    // Set default phases based on phase name or default to Gas
                    let default_phase = match phase_name.to_lowercase().as_str() {
                        "gas" | "vapor" => Phases::Gas,
                        "liquid" => Phases::Liquid,
                        "solid" => Phases::Solid,
                        _ => Phases::Gas,
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
                    phase_data.search_substances();
                    // if not found go to NIST
                    if search_in_NIST {
                        phase_data.if_not_found_go_NIST()?;
                    }
                    // Extract thermal coefficients at standard temperature
                    //  let _ = phase_data.extract_all_thermal_coeffs(298.15);

                    // Add to phase collection
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
