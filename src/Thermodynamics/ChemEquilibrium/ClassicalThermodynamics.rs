use crate::Thermodynamics::User_substances::{DataType, SubsData};

use crate::Thermodynamics::User_PhaseOrSolution::{
    CustomSubstance, SubstancesContainer, ThermodynamicsCalculatorTrait,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};
use prettytable::{Cell, Row, Table, row};

use std::collections::{HashMap, HashSet};
use std::fmt;

use std::{f64, vec};
////////////////////////////////////////////////////////////////////////////////////////////
/// Extension trait for CustomSubstance to provide a unified calculation interface
pub trait ThermodynamicCalculations {
    /// Performs a complete thermodynamic analysis at the specified conditions
    ///
    /// # Arguments
    /// * `temperature` - Temperature in K
    /// * `pressure` - Pressure in Pa
    /// * `concentrations` - Optional vector of mole/mass fractions
    ///
    /// # Returns
    /// A comprehensive thermodynamic analysis including Gibbs energy and entropy
    fn create_thermodynamics(
        &mut self,
        temperature: f64,
        pressure: f64,
        non_zero_moles_number: Option<
            HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        >,

        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> Result<Thermodynamics, String>;
}

impl ThermodynamicCalculations for CustomSubstance {
    fn create_thermodynamics(
        &mut self,
        temperature: f64,
        pressure: f64,
        // non-zero concentrations of substances in the mixture
        non_zero_moles_number: Option<
            HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        >,

        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> Result<Thermodynamics, String> {
        let mut thermodynamics = Thermodynamics::new();
        thermodynamics.subdata = self.clone();
        let subs_container = self.extract_SubstancesContainer()?;
        thermodynamics.subs_container = subs_container;
        if let Some(non_zero_moles_number) = non_zero_moles_number {
            let full_map_of_moles_number = &self
                .create_full_map_of_mole_numbers(non_zero_moles_number)
                .unwrap();
            thermodynamics.map_of_map_mole_numbers_and_phases = full_map_of_moles_number.0.clone();
            thermodynamics.map_of_vec_mole_numbers_and_phases = full_map_of_moles_number.1.clone();
            thermodynamics.map_of_concentration = full_map_of_moles_number.2.clone();
            thermodynamics
                .calculate_Gibbs_free_energy(temperature, full_map_of_moles_number.1.clone());
            thermodynamics.calculate_S(temperature, full_map_of_moles_number.1.clone());
            let vec_of_substances = self.get_all_substances();
            thermodynamics.vec_of_subs = vec_of_substances;
            thermodynamics.calculate_elem_composition_and_molar_mass(groups);
            thermodynamics.initial_composition().unwrap();
            thermodynamics.create_indexed_variables();
        } else {
            let vec_of_substances = self.get_all_substances();
            thermodynamics.vec_of_subs = vec_of_substances;
            thermodynamics.calculate_elem_composition_and_molar_mass(groups);
            // with no concentrations specified, no initiial composition is calculated, self.solver.b0 is not calculated, and Lambda vector is created empty
            thermodynamics.create_indexed_variables();
        }
        //   let vec_of_substances = subs_container.get_all_substances();

        thermodynamics.T = temperature;
        thermodynamics.P = pressure;

        thermodynamics.calculate_Gibbs_sym(temperature);
        thermodynamics.calculate_Gibbs_fun(temperature);

        thermodynamics.calculate_S_sym(temperature);
        thermodynamics.calculate_S_fun(temperature);

        thermodynamics.set_P_to_sym();

        Ok(thermodynamics)
    }
}
///////////////////////////////////////////////////////////////////////////////////
/// strucutre for chemical thermodynamics calculations. Calculaton of Gibbs free energy of a given mixure of substances (numerical result at given T, P,
///  concentration) and symbolic Gibbs free energy of a given mixure of substances (symbolic result at given T, P, concentration).
pub struct Thermodynamics {
    pub subs_container: SubstancesContainer,
    /// map of concentration of substances
    pub vec_of_subs: Vec<String>,
    /// map of values {Option<phase or solution>, Option{total number of moles in phase or solution, map of moles of substances in phase or solution}}
    pub map_of_map_mole_numbers_and_phases:
        HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
    /// map of values {Option<phase or solution>, {total number of moles in phase or solution, vec of moles of substances in phase or solution}}  
    pub map_of_vec_mole_numbers_and_phases:
        HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    /// initial composition of the mixture
    pub map_of_concentration: HashMap<String, f64>,
    /// map of symbolic variables {phase or solutuion, {total number of moles in phase or solution, symbolic vector of number of moles}}
    pub symbolic_vars: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,

    pub symbolic_vars_for_every_subs: HashMap<Option<String>, HashMap<String, Expr>>,
    /// stoichiometric matrix of the reactions, where columns are substances and rows are reactions
    pub stoichio_matrix: Option<Vec<Vec<f64>>>,
    /// matrix of elements composition of substances, where columns are elements and rows are substances
    pub elem_composition_matrix: Option<DMatrix<f64>>,
    /// vector of elements names of the substances
    pub unique_elements: Vec<String>,
    ///
    pub initial_vector_of_elements: Vec<f64>,
    pub T: f64,
    pub P: f64,
    /// INSTANCE OF STRUCTURE SUBSDATA CREATED TO SEARCH SUBSTANCES DATA IN THE DATABASES, STORE SEARCH RESULTS AND CALCULATE THERMO PROPERTIES
    pub subdata: CustomSubstance,
    /// hashmap of entropy of a given mixure of substances (numerical result at given T, P, concentration)
    pub dS: HashMap<Option<String>, HashMap<String, f64>>,
    /// hashmap {phase or solution name:{ substance name: entropy function(T, vector moles of substances, totsl number of moles in phase or solution)} }
    pub dS_fun: HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    >,
    /// hashmap of entropy of a given mixure of substances (symbolic result at given T, P, concentration)
    pub dS_sym: HashMap<Option<String>, HashMap<String, Expr>>,
    /// hashmap of Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    pub dG: HashMap<Option<String>, HashMap<String, f64>>,
    /// hashmap of Gibbs free energy of a given mixure of substances (function at given T, P, concentration)
    pub dG_fun: HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    >,
    /// hashmap of Gibbs free energy of a given mixure of substances (symbolic result at given T, P, concentration)
    pub dG_sym: HashMap<Option<String>, HashMap<String, Expr>>,
    /// Structure for solver of nonlinear equations
    pub solver: Solver,
}

impl Clone for Thermodynamics {
    fn clone(&self) -> Self {
        let mut new_dG = HashMap::new();

        // We can't directly clone the functions, so we'll need to handle this field specially
        // For now, we'll create an empty map structure that matches the original
        for (phase_or_solution, map_dG_fun) in &self.dG_fun {
            let mut inner_map = HashMap::new();
            for (substance, _) in map_dG_fun {
                let placeholder: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static> =
                    Box::new(|_: f64, _: Option<Vec<f64>>, _: Option<f64>| {
                        // your function implementation here
                        0.0
                    });
                inner_map.insert(substance.clone(), placeholder);
            }
            new_dG.insert(phase_or_solution.clone(), inner_map);
        }
        let mut new_dS = HashMap::new();

        // We can't directly clone the functions, so we'll need to handle this field specially
        // For now, we'll create an empty map structure that matches the original
        for (phase_or_solution, map_dS_fun) in &self.dS_fun {
            let mut inner_map = HashMap::new();
            for (substance, _) in map_dS_fun {
                let placeholder: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static> =
                    Box::new(|_: f64, _: Option<Vec<f64>>, _: Option<f64>| {
                        // your function implementation here
                        0.0
                    });
                inner_map.insert(substance.clone(), placeholder);
            }
            new_dS.insert(phase_or_solution.clone(), inner_map);
        }
        let placeholder: Vec<Box<dyn Fn(DVector<f64>) -> f64>> = vec![
            Box::new(|x: DVector<f64>| x[0]),
            Box::new(|x: DVector<f64>| x.sum()),
        ];
        let solver_cloned = Solver {
            b0: self.solver.b0.clone(),
            elements_conditions: placeholder,
            elements_conditions_sym: self.solver.elements_conditions_sym.clone(),
            Lambda: self.solver.Lambda.clone(),
            n: self.solver.n.clone(),
            Np: self.solver.Np.clone(),
            eq_mu: self.solver.eq_mu.clone(),
            eq_mu_fun: Box::new(|_, _, _, _| vec![]),
            eq_sum_mole_numbers: self.solver.eq_sum_mole_numbers.clone(),
            full_system_sym: self.solver.full_system_sym.clone(),
            all_unknowns: self.solver.all_unknowns.clone(),
        };
        let mut cloned_instance = Thermodynamics {
            subs_container: self.subs_container.clone(),
            vec_of_subs: self.vec_of_subs.clone(),
            map_of_map_mole_numbers_and_phases: self.map_of_map_mole_numbers_and_phases.clone(),
            map_of_vec_mole_numbers_and_phases: self.map_of_vec_mole_numbers_and_phases.clone(),
            map_of_concentration: self.map_of_concentration.clone(),
            symbolic_vars: self.symbolic_vars.clone(),
            symbolic_vars_for_every_subs: self.symbolic_vars_for_every_subs.clone(),

            T: self.T,
            P: self.P,
            stoichio_matrix: self.stoichio_matrix.clone(),
            elem_composition_matrix: self.elem_composition_matrix.clone(),
            unique_elements: self.unique_elements.clone(),
            initial_vector_of_elements: self.initial_vector_of_elements.clone(),
            subdata: self.subdata.clone(),
            dS: self.dS.clone(),
            dS_fun: new_dS,
            dS_sym: self.dS_sym.clone(),
            dG: self.dG.clone(),
            dG_fun: new_dG,
            dG_sym: self.dG_sym.clone(),
            solver: solver_cloned,
        };
        // Clone the functions
        let _ = cloned_instance.calculate_Gibbs_fun(self.T);
        let _ = cloned_instance.calculate_S_fun(self.T);
        let _ = cloned_instance.composition_equations();
        let _ = cloned_instance.create_nonlinear_system_fun();
        cloned_instance
    }
}

impl fmt::Debug for Thermodynamics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Thermodynamics")
            .field("subs_container", &self.subs_container)
            .field("vec_of_subs", &self.vec_of_subs)
            .field("T", &self.T)
            .field("P", &self.P)
            .field("subdata", &self.subdata)
            .field("dmu", &self.dG)
            //.field("dmu_fun", &self.dmu_fun.is_some())
            .field("dmu_sym", &self.dG_sym)
            .finish()
    }
}
impl Thermodynamics {
    pub fn new() -> Self {
        let solver = Solver::new();
        Self {
            subs_container: SubstancesContainer::SinglePhase(Vec::new()),
            vec_of_subs: Vec::new(),
            map_of_map_mole_numbers_and_phases: HashMap::new(),
            map_of_vec_mole_numbers_and_phases: HashMap::new(),
            map_of_concentration: HashMap::new(),
            symbolic_vars: HashMap::new(),
            symbolic_vars_for_every_subs: HashMap::new(),
            T: 298.15,
            P: 1e5,
            stoichio_matrix: None,
            elem_composition_matrix: None,
            unique_elements: Vec::new(),
            initial_vector_of_elements: Vec::new(),
            subdata: CustomSubstance::OnePhase(SubsData::new()),
            dG: HashMap::new(),
            dG_fun: HashMap::new(),
            dG_sym: HashMap::new(),
            dS: HashMap::new(),
            dS_fun: HashMap::new(),
            dS_sym: HashMap::new(),
            solver: solver,
        }
    }
    /////////////////////////////////////////////////////////////SETTERS//////////////////////////////////////////////
    pub fn set_T(&mut self, T: f64) {
        self.T = T;
    }
    pub fn set_P(&mut self, P: f64, _P_unit: Option<String>) {
        self.P = P;
    }
    ///  a method to update substances based on subdata
    pub fn update_substances_from_subdata(&mut self) {
        match &self.subdata {
            CustomSubstance::OnePhase(subdata) => {
                self.subs_container = SubstancesContainer::SinglePhase(subdata.substances.clone());
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let mut phase_substances = HashMap::new();
                for (phase_name, subsdata) in &phase_or_solution.subs_data {
                    if let Some(phase_name) = phase_name {
                        phase_substances.insert(phase_name.clone(), subsdata.substances.clone());
                    }
                }
                self.subs_container = SubstancesContainer::MultiPhase(phase_substances);
            }
        }
    }
    /// function to calculate chemical elements composition and molar mass of all substances
    pub fn calculate_elem_composition_and_molar_mass(
        &mut self,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) {
        let (composition_matrix, _hasmap_of_molar_mass, vec_of_elements) = self
            .subdata
            .calculate_elem_composition_and_molar_mass(groups)
            .unwrap();
        //  println!("composition_matrix = {:#?} \n vec_of_elements = {:#?}", composition_matrix, vec_of_elements);
        self.elem_composition_matrix = Some(composition_matrix);
        self.unique_elements = vec_of_elements;
    }
    ///create symbolic variables for equations:
    ///  - Lambda - lafrangian multipliers ( dG per chemical element)
    ///  - Np - mole numbers for each phase or solution
    ///  - n - number of moles of each substance
    pub fn create_indexed_variables(&mut self) {
        let symbolic_vars = self.subdata.indexed_moles_variables().unwrap();
        self.symbolic_vars = symbolic_vars.0;
        let number_of_elements = self.solver.b0.len();
        //  if number_of_elements == 0 {panic!("number_of_elements == 0");}
        let Lambda = Expr::IndexedVars(number_of_elements, "Lambda").0;
        self.symbolic_vars_for_every_subs = symbolic_vars.3;
        self.solver.Lambda = Lambda;
        self.solver.n = symbolic_vars.1;
        self.solver.Np = symbolic_vars.2;
    }
    /// function takes map {Option<Phase_name> : (Option<Number of moles in phase>, Option<{substance_name : Number of moles}>}> }
    /// where Number of moles are non zero and creates full map of mole numbers
    ///
    pub fn set_number_of_moles(
        &mut self,
        non_zero_moles_number: Option<
            HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        >,
    ) {
        if let Some(non_zero_moles_number) = non_zero_moles_number {
            let full_map_of_moles_number = &self
                .subdata
                .create_full_map_of_mole_numbers(non_zero_moles_number)
                .unwrap();
            self.map_of_map_mole_numbers_and_phases = full_map_of_moles_number.0.clone();
            self.map_of_vec_mole_numbers_and_phases = full_map_of_moles_number.1.clone();
            self.map_of_concentration = full_map_of_moles_number.2.clone();
        }
    }

    //////////////////////////////////////////dG//////////////////////////////////////////////

    /// Function for calculating Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    pub fn calculate_Gibbs_free_energy(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) {
        self.T = T;
        let P = self.P;
        let sd = &mut self.subdata;
        self.dG = sd.calcutate_Gibbs_free_energy(T, P, n).unwrap();
    }

    /// Symbolic function for calculating Gibbs free energy of a given mixure of substances (symbolic function at given T, P, concentration)
    pub fn calculate_Gibbs_sym(&mut self, T: f64) {
        self.T = T;
        let n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)> =
            self.symbolic_vars.clone();
        let sd = &mut self.subdata;
        self.dG_sym = sd.calculate_Gibbs_sym(T, n).unwrap();
    }
    /// Function for making Gibbs free energy closure of a given mixure of substances (as function of given T, P, concentration)
    pub fn calculate_Gibbs_fun(&mut self, T: f64) {
        self.T = T;
        let P = self.P;
        let sd = &mut self.subdata;
        self.dG_fun = sd.calculate_Gibbs_fun(T, P);
    }
    ////////////////////////////////////////////S - ENTHROPY///////////////////////////////////////////////////////////
    /// Function for calculating enthropy of a given mixure of substances (numerical result at given T, P, concentration)
    pub fn calculate_S(
        &mut self,
        T: f64,
        n: HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    ) {
        self.T = T;
        let P = self.P;
        let sd = &mut self.subdata;
        self.dS = sd.calculate_S(T, P, n).unwrap();
    }

    ///Symbolic  function for calculating enthropy of a given mixure of substances ( symbolic result at given T, P, concentration)
    pub fn calculate_S_sym(&mut self, T: f64) {
        self.T = T;
        let n: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)> =
            self.symbolic_vars.clone();
        let sd = &mut self.subdata;
        self.dS_sym = sd.calculate_S_sym(T, n).unwrap();
    }

    /// Function for making enthropy closure of a given mixure of substances (as function of given T, P, concentration)
    pub fn calculate_S_fun(&mut self, T: f64) {
        self.T = T;
        let sd = &mut self.subdata;
        self.dS_fun = sd.calculate_S_fun(T, self.P);
    }
    /// convert synbolic variable P in G symbolic function to numerical value
    pub fn set_P_to_sym(&mut self) {
        let P = self.P;
        for (_phase_or_solution_name, sym_fun) in self.dG_sym.iter_mut() {
            for (_, sym_fun) in sym_fun.iter_mut() {
                *sym_fun = sym_fun.set_variable("P", P).symplify()
            }
        }
    }
    /// convert synbolic variable T in G symbolic function to numerical value
    pub fn set_T_to_sym(&mut self) {
        let T = self.T;
        for (_phase_or_solution_name, sym_fun) in self.dG_sym.iter_mut() {
            for (_, sym_fun) in sym_fun.iter_mut() {
                *sym_fun = sym_fun.set_variable("T", T).symplify()
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Functions for initial composition
    /// This function calculates the initial composition of elements in a system given the composition matrix and initial concentrations
    /// of substances. It checks if the number of columns matches the number of elements, then performs a matrix multiplication to
    /// calculate the initial amount of each element. The result is stored in self.initial_vector_of_elements and self.solver.b0.
    pub fn initial_composition(&mut self) -> Result<(), std::io::Error> {
        let vec_of_substances = self.vec_of_subs.clone();

        if let Some(A) = self.elem_composition_matrix.clone() {
            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
                println!(
                    "A.ncols() {} != vec_of_substances.len() {}",
                    A.ncols(),
                    vec_of_substances.len()
                );
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "The number of columns is not equal to the number of elements",
                ));
            }
            let vec_of_initional_concentrations: Vec<f64> = self
                .vec_of_subs
                .iter()
                .map(|v| {
                    let conc = self.map_of_concentration.get(v).unwrap();
                    *conc
                })
                .collect();

            //let vec_of_initional_concentrations = DVector::from_vec(vec_of_initional_concentrations);
            let mut initial_vector_of_elements: Vec<f64> = Vec::new();
            for i in 0..n_elements {
                let row_of_element_i = A.row(i).iter().map(|&v| v).collect::<Vec<f64>>();
                let row_of_element_i: DVector<_> = DVector::from_vec(row_of_element_i);

                let vec_of_initional_concentrations =
                    DVector::from_vec(vec_of_initional_concentrations.clone());
                if vec_of_initional_concentrations.len() != row_of_element_i.len() {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        "The number of substances is not equal to the number of elements",
                    ));
                }

                let b_i = row_of_element_i.dot(&vec_of_initional_concentrations);
                //  let b_i = b_i_vec.sum();
                initial_vector_of_elements.push(b_i);
            }
            self.initial_vector_of_elements = initial_vector_of_elements;
            self.solver.b0 = self.initial_vector_of_elements.clone();
        }
        Ok(())
    }
    /// Function to create closure for conservation of  the elements  condition
    ///  SUM_j (a_ij * x_j) = b_i
    pub fn composition_equations(&mut self) -> Result<(), std::io::Error> {
        if let Some(A) = self.elem_composition_matrix.clone() {
            let vec_of_substances = self.vec_of_subs.clone();

            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "The number of substances is not equal to the number of elements",
                ));
            }

            let mut vector_of_conditions: Vec<Box<dyn Fn(DVector<f64>) -> f64>> = Vec::new();
            for i in 0..n_elements {
                let initial_vector_of_elements = self.initial_vector_of_elements.clone();
                let row_of_element_i = A.row(i).iter().map(|&v| v).collect::<Vec<f64>>();
                let row_of_element_i: DVector<_> = DVector::from_vec(row_of_element_i);
                let condition_closure = move |vec_of_concentrations: DVector<f64>| {
                    let b_i = row_of_element_i.dot(&vec_of_concentrations);

                    let condition_i = b_i - initial_vector_of_elements[i].clone();
                    condition_i
                };
                vector_of_conditions.push(Box::new(condition_closure));
            }
            self.solver.elements_conditions = vector_of_conditions;
        }
        Ok(())
    }
    /// This function generates symbolic equations for conservation of elements in a chemical reaction. It takes a vector of symbolic
    ///  concentrations as input and returns a vector of symbolic equations representing the conservation of each element.
    /// SUM_j (a_ij * x_j) = b_i
    /// The equations are constructed by multiplying the element composition matrix (transposed) with the symbolic concentrations and
    ///  equating the result to the initial amount of each element. The resulting equations are stored in self.solver.elements_conditions_sym.
    pub fn composition_equation_sym(&mut self) -> Result<(), std::io::Error> {
        let vec_of_concentrations_sym: Vec<Expr> = self.solver.n.clone();
        if let Some(A) = self.elem_composition_matrix.clone() {
            let vec_of_substances = self.vec_of_subs.clone();

            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
                println!(
                    "The number of columns {} is not equal to the number of elements {}",
                    A.ncols(),
                    vec_of_substances.len()
                );
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "The number of substances is not equal to the number of elements",
                ));
            }

            let initial_vector_of_elements = self.initial_vector_of_elements.clone();
            let initial_vector_of_element_sym = initial_vector_of_elements
                .iter()
                .map(|&v| Expr::Const(v))
                .collect::<Vec<Expr>>();

            let mut vector_of_conditions: Vec<Expr> = Vec::new();
            for i in 0..n_elements {
                let row_of_element_i = A
                    .row(i)
                    .iter()
                    .map(|&v| Expr::Const(v))
                    .collect::<Vec<Expr>>();

                let b_i_sum: Expr = row_of_element_i
                    .iter()
                    .zip(vec_of_concentrations_sym.iter())
                    .map(|(row_of_element_i, vec_of_concentrations_sym)| {
                        row_of_element_i.clone() * vec_of_concentrations_sym.clone()
                    })
                    .fold(Expr::Const(0.0), |acc, x| acc + x);

                let condition_i =
                    b_i_sum.symplify() - initial_vector_of_element_sym[i].clone().symplify();
                vector_of_conditions.push(condition_i.symplify());
            }
            self.solver.elements_conditions_sym = vector_of_conditions;
        }
        Ok(())
    }

    /// Creates the nonlinear system of equations for the classical thermodynamics solver.
    ///
    /// This function creates the system of equations needed to solve the problem of finding the
    /// chemical equilibrium composition of a system. The system of equations is built using the
    /// Lagrange multipliers method, which is used to minimize the total free energy of the system
    /// subject to the constraints of element conservation.
    ///
    /// The function takes the transpose of the element composition matrix, and then uses the
    /// `calculate_Lagrange_equations_sym` function of the `SubsData` struct to calculate the
    /// Lagrange equations for the system. Finally, the function stores the equations in the
    /// `solver` struct.
    ///
    /// # Errors
    ///
    /// The function returns an `Err` if the `elem_composition_matrix` is not set in the `SubsData`
    /// struct, or if the `calculate_Lagrange_equations_sym` function fails.
    pub fn create_nonlinear_system_sym(&mut self) -> Result<(), std::io::Error> {
        let A = self.elem_composition_matrix.clone().unwrap();
        // println!("A: {}, ncols: {}, nrows: {}", A, A.ncols(), A.nrows()) ;
        let eq = self
            .subdata
            .calculate_Lagrange_equations_sym(A, self.dG_sym.clone())
            .unwrap();
        self.solver.eq_mu = eq;
        Ok(())
    }

    pub fn create_nonlinear_system_fun(&mut self) -> Result<(), std::io::Error> {
        let T = self.T;
        self.calculate_Gibbs_fun(T);
        let A = self.elem_composition_matrix.clone().unwrap();

        let dG_fun = self.subdata.calculate_Gibbs_fun(T, self.P);
        let eq = self.create_nonlinear_system_fun_local(A, dG_fun).unwrap();
        self.solver.eq_mu_fun = eq;
        Ok(())
    }

    fn create_nonlinear_system_fun_local(
        &mut self,
        A: DMatrix<f64>,
        dG_fun: HashMap<
            Option<String>,
            HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
        >,
    ) -> Result<
        Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
        std::io::Error,
    > {
        let Lagrange_equations: Box<
            dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static,
        > = self
            .subdata
            .calculate_Lagrange_equations_fun(A, dG_fun)
            .unwrap();
        let eq = move |T: f64,
                       vec_of_concentrations: Option<Vec<f64>>,
                       Np: Option<f64>,
                       Lambda: Vec<f64>| {
            Lagrange_equations(T, vec_of_concentrations, Np, Lambda)
        };
        let f = Box::new(eq)
            as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>;
        Ok(f)
    }
    /// creates symbolic equations representing the idea that sum  of mole numbers in this phase is equal to the total number of moles in this phase
    pub fn create_sum_of_mole_numbers_sym(&mut self) -> Result<Vec<Expr>, std::io::Error> {
        let sym_variables = self.symbolic_vars.clone();
        let mut vec_of_sum_of_mole_numbers: Vec<Expr> = Vec::new();
        for (_phase, phase_data) in sym_variables.iter() {
            let phase_data = phase_data.clone();
            if let (Some(Np), Some(vec_of_mole_numbers)) = (phase_data.0, phase_data.1) {
                let sum_of_mole_numbers = vec_of_mole_numbers
                    .iter()
                    .fold(Expr::Const(0.0), |acc, x| acc + x.clone());
                let eq = sum_of_mole_numbers.symplify() - Np;
                vec_of_sum_of_mole_numbers.push(eq);
            }
        }
        //  println!("vec_of_sum_of_mole_numbers {:?}", vec_of_sum_of_mole_numbers);
        self.solver.eq_sum_mole_numbers = vec_of_sum_of_mole_numbers.clone();
        Ok(vec_of_sum_of_mole_numbers)
    }
    /// construct the full system of equations for the nonlinear solver.
    /// It takes the Lagrange equations, the element conservation equations, and the sum of mole numbers equations and creates a single system of equations.
    /// Check if the number of equations is equal to the number of variables
    pub fn form_full_system_sym(&mut self) -> Result<(), std::io::Error> {
        let mut eq = self.solver.eq_mu.clone();
        //adding equations SUM j (a_ij * x_j) - b_i = 0
        eq.extend(self.solver.elements_conditions_sym.clone());
        eq.extend(self.solver.eq_sum_mole_numbers.clone());
        //  println!("equation system {:?}", eq);
        let n_eq = eq.len(); // number of equations
        let mut x = self.solver.n.clone();
        x.extend(self.solver.Lambda.clone());
        x.extend(self.solver.Np.clone());
        let n_vars = x.len(); // number of variables

        if n_eq != n_vars {
            println!(
                "The number of equations {}:{:?} the number of variables {}:{:?}",
                n_eq, eq, n_vars, x
            );
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "The number of equations is not equal to the number of variables",
            ));
        }
        self.solver.full_system_sym = eq;
        self.solver.all_unknowns = x;
        Ok(())
    }
    /// create symbolic variables for the nonlinear solver (function for checking and debugging)
    pub fn symbolic_variables_extract(&mut self) -> Result<Vec<Expr>, std::io::Error> {
        let mut vec_of_subs: Vec<String> = Vec::new();
        for eq in &self.solver.eq_mu {
            let vars = eq.all_arguments_are_variables();
            let mut vars = vars
                .iter()
                .map(|v| v.trim().to_string())
                .collect::<Vec<String>>();
            vars.dedup();

            vec_of_subs.extend(vars);
        }
        let vars_unique: HashSet<String> = HashSet::from_iter(vec_of_subs.clone());
        let mut vars_unique: Vec<String> = vars_unique.into_iter().collect();
        vars_unique.sort();
        //  println!("vec_of_subs: {:#?}", vars_unique);
        let vec_of_subs =
            Expr::parse_vector_expression(vars_unique.iter().map(|v| v.as_str()).collect());
        Ok(vec_of_subs)
    }
    /// shortcut function to create the full system equations for a given T and P
    /// runs the following functions:
    ///  set_P_to_sym();
    ///  composition_equations().unwrap();
    ///  composition_equation_sym().unwrap();
    ///  set_T_to_sym();
    ///  create_nonlinear_system_sym().unwrap();
    ///  create_nonlinear_system_fun().unwrap();
    ///  create_sum_of_mole_numbers_sym().unwrap();
    ///  form_full_system_sym().unwrap();
    ///  pretty_print_full_system();
    pub fn find_composition_for_const_TP(&mut self) -> Result<(), std::io::Error> {
        self.set_P_to_sym();
        self.composition_equations().unwrap();
        self.composition_equation_sym().unwrap();
        self.set_T_to_sym();
        self.create_nonlinear_system_sym().unwrap();
        self.create_nonlinear_system_fun().unwrap();
        self.create_sum_of_mole_numbers_sym().unwrap();
        self.form_full_system_sym().unwrap();
        self.pretty_print_full_system();

        Ok(())
    }
    // pub fn construct_nonlinear_syste(&mut self) -> Result<(), std::io::Error> {}
    ////////////////////////INPUT/OUTPUT////////////////////////////////////////////////////////

    /// Prints the data of the substance to the console
    pub fn pretty_print_thermo(&self) -> Result<(), std::io::Error> {
        println!("__________subs properties at {} K__________", self.T);
        let mut table = Table::new();

        let subs_container = match &self.subs_container {
            SubstancesContainer::SinglePhase(vec) => vec.clone(),
            SubstancesContainer::MultiPhase(map) => {
                map.iter().map(|(substance, _)| substance.clone()).collect()
            }
        };

        let sd = &self.subdata;
        match sd {
            CustomSubstance::OnePhase(sd) => {
                // Add the header row
                table.add_row(row!["substance", "Cp", "dH", "dS", "dG",]);
                // Add the data rows
                for (substance, data) in self.dG.get(&None).unwrap() {
                    let map_property_values =
                        sd.therm_map_of_properties_values.get(substance).unwrap();
                    let Cp = map_property_values.get(&DataType::Cp).unwrap().unwrap();
                    let dh = map_property_values.get(&DataType::dH).unwrap().unwrap();
                    let ds = map_property_values.get(&DataType::dS).unwrap().unwrap();
                    let dG = data.clone();

                    table.add_row(row![substance, Cp, dh, ds, dG,]);
                }
                table.printstd();
                println!("_____________________________________________________________");
                println!(
                    "\n_________ANALYTIC GIBBS FREE ENERGY AT {} K__________________________________",
                    self.T
                );
                let mut table2 = Table::new();
                table2.add_row(row!["substance", "dG_sym",]);
                for substance in &subs_container {
                    let dG_sym = self
                        .dG_sym
                        .get(&None)
                        .unwrap()
                        .get(substance)
                        .unwrap()
                        .clone();
                    table2.add_row(row![substance, dG_sym,]);
                }
                table2.printstd();
                println!("_____________________________________________________________");
            } // CustomSubstance::OnePhase(sd)
            CustomSubstance::PhaseOrSolution(sd) => {
                for (phase_name, subsdata) in sd.subs_data.iter() {
                    let sd = subsdata;
                    if let Some(phase_name) = phase_name {
                        println!("__________subs properties at {} K__________", self.T);
                        let mut table = Table::new();
                        // Add the header row
                        table.add_row(row![phase_name]);
                        table.add_row(row!["substance", "Cp", "dH", "dS", "dG",]);
                        for (substance, data) in self.dG.get(&Some(phase_name.clone())).unwrap() {
                            let map_property_values =
                                sd.therm_map_of_properties_values.get(substance).unwrap();
                            let Cp = map_property_values.get(&DataType::Cp).unwrap().unwrap();
                            let dh = map_property_values.get(&DataType::dH).unwrap().unwrap();
                            let ds = map_property_values.get(&DataType::dS).unwrap().unwrap();
                            let dG = data.clone();
                            table.add_row(row![substance, Cp, dh, ds, dG,]);
                        }
                        table.printstd();
                        println!("_____________________________________________________________");
                        println!(
                            "\n_________ANALYTIC GIBBS FREE ENERGY AT {} K__________________________________",
                            self.T
                        );
                        let mut table2 = Table::new();
                        table2.add_row(row!["substance", "dG_sym",]);
                        for substance in &subs_container {
                            let dG_sym = self
                                .dG_sym
                                .get(&Some(phase_name.clone()))
                                .unwrap()
                                .get(substance)
                                .unwrap()
                                .clone();
                            table2.add_row(row![substance, dG_sym,]);
                        }
                        table2.printstd();
                        println!("_____________________________________________________________");
                    } // if let Some(phase_name) = phase_name
                } // for (phase_name, subsdata) in sd.subs_data.iter()
            } // CustomSubstance::PhaseOrSolution(sd)
        } // match
        Ok(())
    }

    pub fn pretty_print_substances_verbose(&self) -> Result<(), std::io::Error> {
        let mut elem_table = Table::new();
        println!("___________________ELEMENT COMPOSITION MATRIX________________________");
        let mut header_row = vec![Cell::new("Substances/Elements")];
        let unique_vec_of_elems = self.unique_elements.clone();
        let elem_matrix = self.elem_composition_matrix.clone().unwrap();
        let subs = self.vec_of_subs.clone();
        if subs.len() == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "No substances in the system",
            ));
        }
        for elems in unique_vec_of_elems.clone() {
            header_row.push(Cell::new(&elems));
        }
        elem_table.add_row(Row::new(header_row));
        for (i, sub) in subs.iter().enumerate() {
            let mut row = vec![Cell::new(sub)];
            for j in 0..unique_vec_of_elems.len() {
                row.push(Cell::new(&format!(
                    "{:.4}",
                    elem_matrix.get((i, j)).unwrap()
                )));
            }
            elem_table.add_row(Row::new(row));
        }
        elem_table.printstd();
        println!("_____________________________________________________________");
        Ok(())
    } //pretty_print_substances_verbose
    /// Print Lagrange equations as a table
    pub fn pretty_print_Lagrange_equations(&self) -> Result<(), std::io::Error> {
        println!("___________________NONLINEAR EQUATIONS________________________");
        let mut elem_table = Table::new();
        let subs = self.vec_of_subs.clone();
        if subs.len() == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "No substances in the system",
            ));
        }
        if self.solver.eq_mu.len() == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "No Lagrange equations in the system",
            ));
        }
        for (i, sub) in subs.iter().enumerate() {
            let mut row = vec![Cell::new(sub)];
            let eq_i = self.solver.eq_mu[i].clone();
            row.push(Cell::new(&format!("{}", eq_i)));
            elem_table.add_row(Row::new(row));
        }
        elem_table.printstd();
        println!("_____________________________________________________________");
        Ok(())
    }
    pub fn pretty_print_composition_equations(&self) -> Result<(), std::io::Error> {
        println!("___________________COMPOSITION EQUATIONS________________________");
        let mut elem_table = Table::new();
        let Lambda = self.solver.Lambda.clone();
        if Lambda.len() == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "No substances in the system",
            ));
        }
        for (i, eq) in self.solver.elements_conditions_sym.iter().enumerate() {
            let mut row = vec![Cell::new(&Lambda[i].to_string().clone())];

            row.push(Cell::new(&eq.to_string().clone()));
            elem_table.add_row(Row::new(row));
        }
        elem_table.printstd();
        println!("_____________________________________________________________");
        Ok(())
    }
    /// Print sum of mole numbers as a table
    pub fn pretty_print_sum_mole_numbers(&self) -> Result<(), std::io::Error> {
        println!("___________________SUM OF MOLE NUMBERS________________________");
        let mut elem_table = Table::new();
        let Np = self.solver.Np.clone();
        if Np.len() == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "No substances in the system",
            ));
        }
        for (i, eq) in self.solver.eq_sum_mole_numbers.iter().enumerate() {
            let mut row = vec![Cell::new(&Np[i].to_string().clone())];
            row.push(Cell::new(&eq.to_string().clone()));
            elem_table.add_row(Row::new(row));
        }
        elem_table.printstd();
        println!("_____________________________________________________________");
        Ok(())
    }
    /// Print all Lagrange equations, composition equations and sum of mole numbers as a table
    pub fn pretty_print_full_system(&self) {
        self.pretty_print_Lagrange_equations().unwrap();
        self.pretty_print_composition_equations().unwrap();
        self.pretty_print_sum_mole_numbers().unwrap();
    }
}
/// structure for solving chemical equilibrium problem
/// its fields include:
/// 1) initial elements composition
/// 2) Lagrange equations
/// 3) equations representing the idea that sum  of mole numbers in this phase is equal to the total number of moles in this phase  
/// 4) unknowns: mole numbers and Lagrange multipliers
/// 5) instance of the solver of the system of nonlinear equations
pub struct Solver {
    /// Initial elements composition
    pub b0: Vec<f64>,
    /// closure for elements conditions: SUM_j(a_ij*n_j) = b_i =0  , i = 1..m
    pub elements_conditions: Vec<Box<dyn Fn(DVector<f64>) -> f64>>,
    /// symbolic elements conditions: SUM_j(a_ij*n_j) = b_i =0  , i = 1..m
    pub elements_conditions_sym: Vec<Expr>,
    // Lagrange multipliers
    pub Lambda: Vec<Expr>,
    /// mole numbers
    pub n: Vec<Expr>,
    ///  mole numbers of each phase
    pub Np: Vec<Expr>,
    /// equations
    pub eq_mu: Vec<Expr>, // T, n,
    ///
    pub eq_mu_fun: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    ///  symbolic equations representing the idea that sum  of mole numbers in this phase is equal to the total number of moles in this phase
    pub eq_sum_mole_numbers: Vec<Expr>,
    /// full system of equations
    pub full_system_sym: Vec<Expr>,
    /// Lambda, n, Np,
    pub all_unknowns: Vec<Expr>,
}
impl Solver {
    pub fn new() -> Self {
        Self {
            b0: vec![],
            elements_conditions: vec![],
            elements_conditions_sym: vec![],
            Lambda: vec![],
            n: vec![],
            Np: vec![],
            eq_mu: vec![],
            eq_mu_fun: Box::new(|_, _, _, _| vec![]),
            eq_sum_mole_numbers: vec![],
            full_system_sym: vec![],
            all_unknowns: vec![],
        }
    }
}
