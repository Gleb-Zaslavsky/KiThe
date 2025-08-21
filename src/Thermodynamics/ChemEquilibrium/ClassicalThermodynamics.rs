use crate::Thermodynamics::User_PhaseOrSolution::{
    CustomSubstance, SubstancesContainer, ThermodynamicsCalculatorTrait,
};

use super::ClassicalThermodynamicsSolver::Solver;
use crate::Thermodynamics::User_substances::SubsData;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::{f64, vec};
/////////////////////ERROR HANDLING////////////////////////////////////////////////////////

use std::io;
#[derive(Debug)]
pub enum ThermodynamicsError {
    MatrixDimensionMismatch(String),
    MissingData(String),
    CalculationError(String),
    SubstanceNotFound(String),
    InvalidComposition(String),
    SolverError(String),
    IoError(io::Error),
    ParseError(String),
}

impl fmt::Display for ThermodynamicsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ThermodynamicsError::MatrixDimensionMismatch(msg) => {
                write!(f, "Matrix dimension mismatch: {}", msg)
            }
            ThermodynamicsError::MissingData(msg) => write!(f, "Missing data: {}", msg),
            ThermodynamicsError::CalculationError(msg) => write!(f, "Calculation error: {}", msg),
            ThermodynamicsError::SubstanceNotFound(msg) => {
                write!(f, "Substance not found: {}", msg)
            }
            ThermodynamicsError::InvalidComposition(msg) => {
                write!(f, "Invalid composition: {}", msg)
            }
            ThermodynamicsError::SolverError(msg) => write!(f, "Solver error: {}", msg),
            ThermodynamicsError::IoError(err) => write!(f, "IO error: {}", err),
            ThermodynamicsError::ParseError(msg) => write!(f, "Parse error: {}", msg),
        }
    }
}

impl std::error::Error for ThermodynamicsError {}

// Conversion from io::Error to ThermodynamicsError
impl From<io::Error> for ThermodynamicsError {
    fn from(error: io::Error) -> Self {
        ThermodynamicsError::IoError(error)
    }
}

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
    ) -> Result<Thermodynamics, ThermodynamicsError>;
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
    ) -> Result<Thermodynamics, ThermodynamicsError> {
        let mut thermodynamics = Thermodynamics::new();
        thermodynamics.subdata = self.clone();
        let subs_container = self
            .extract_SubstancesContainer()
            .map_err(|e| ThermodynamicsError::MissingData(e))?;
        thermodynamics.subs_container = subs_container;
        if let Some(non_zero_moles_number) = non_zero_moles_number {
            let full_map_of_moles_number = &self
                .create_full_map_of_mole_numbers(non_zero_moles_number)
                .map_err(|e| ThermodynamicsError::CalculationError(e))?;
            thermodynamics.map_of_map_mole_numbers_and_phases = full_map_of_moles_number.0.clone();
            thermodynamics.map_of_vec_mole_numbers_and_phases = full_map_of_moles_number.1.clone();
            thermodynamics.map_of_concentration = full_map_of_moles_number.2.clone();

            thermodynamics
                .calculate_Gibbs_free_energy(temperature, full_map_of_moles_number.1.clone());
            thermodynamics.calculate_S(temperature, full_map_of_moles_number.1.clone());
            let vec_of_substances = self.get_all_substances();
            check_task(
                &thermodynamics.map_of_map_mole_numbers_and_phases,
                &vec_of_substances,
            )
            .map_err(|e| ThermodynamicsError::InvalidComposition(format!("{:?}", e)))?;
            thermodynamics.vec_of_subs = vec_of_substances;
            thermodynamics.calculate_elem_composition_and_molar_mass(groups);
            thermodynamics
                .initial_composition()
                .map_err(|e| ThermodynamicsError::CalculationError(format!("{}", e)))?;
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
///
fn check_task(
    full_map_of_moles_number: &HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
    vec_of_substances: &Vec<String>,
) -> Result<(), ThermodynamicsError> {
    let mut vec_of_substances_in_map = Vec::new();
    // println!("full_map_of_moles_number: {:#?}, substances: {:#?}", full_map_of_moles_number, vec_of_substances);

    for value in full_map_of_moles_number.values() {
        let map_of_subs = value.1.clone();
        if let Some(map_of_subs) = map_of_subs {
            for (substance, _) in map_of_subs.iter() {
                vec_of_substances_in_map.push(substance.clone());
            }
        }
    }
    let set_of_substances_in_map: HashSet<String> =
        HashSet::from_iter(vec_of_substances_in_map.clone());
    let set_of_substances_in_vec: HashSet<String> = HashSet::from_iter(vec_of_substances.clone());
    if set_of_substances_in_map != set_of_substances_in_vec {
        return Err(ThermodynamicsError::InvalidComposition(
            "Number of substances in map != number of substances in vector".to_string(),
        ));
    }

    Ok(())
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
    pub Tm: f64,
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
            solution: self.solver.solution.clone(),
            map_of_solutions: self.solver.map_of_solutions.clone(),
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
            Tm: self.Tm,
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
            Tm: 1e5,
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
    ///  - Lambda - laпrangian multipliers ( dG per chemical element)
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
        //  physiscal_nature_of_phase_or_solution: Option<HashMap<String, Phases>>,
    ) -> Result<(), ThermodynamicsError> {
        if let Some(non_zero_moles_number) = non_zero_moles_number {
            let full_map_of_moles_number = &self
                .subdata
                .create_full_map_of_mole_numbers(non_zero_moles_number)
                .map_err(|e| ThermodynamicsError::CalculationError(e))?;
            self.map_of_map_mole_numbers_and_phases = full_map_of_moles_number.0.clone();
            self.map_of_vec_mole_numbers_and_phases = full_map_of_moles_number.1.clone();
            self.map_of_concentration = full_map_of_moles_number.2.clone();
            check_task(&self.map_of_map_mole_numbers_and_phases, &self.vec_of_subs)?;
        }
        Ok(())
    }

    /// convert symbolic variable P in G symbolic function to numerical value
    pub fn set_P_to_sym(&mut self) {
        let P = self.P;
        for (_phase_or_solution_name, sym_fun) in self.dG_sym.iter_mut() {
            for (_, sym_fun) in sym_fun.iter_mut() {
                *sym_fun = sym_fun.set_variable("P", P).symplify()
            }
        }
    }
    /// convert symbolic variable T in G symbolic function to numerical value
    pub fn set_T_to_sym(&mut self) {
        let T = self.T;
        for (_phase_or_solution_name, sym_fun) in self.dG_sym.iter_mut() {
            for (_, sym_fun) in sym_fun.iter_mut() {
                *sym_fun = sym_fun.set_variable("T", T).symplify()
            }
        }
        if self.solver.eq_mu.len() > 0 {
            for mu in self.solver.eq_mu.iter_mut() {
                *mu = mu.set_variable("T", T).symplify()
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Functions for initial composition
    /// This function calculates the initial composition of elements in a system given the composition matrix and initial concentrations
    /// of substances. It checks if the number of columns matches the number of elements, then performs a matrix multiplication to
    /// calculate the initial amount of each element. The result is stored in self.initial_vector_of_elements and self.solver.b0.
    pub fn initial_composition(&mut self) -> Result<(), ThermodynamicsError> {
        let vec_of_substances = self.vec_of_subs.clone();

        if let Some(A) = self.elem_composition_matrix.clone() {
            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
                println!(
                    "A.ncols() {} != vec_of_substances.len() {}",
                    A.ncols(),
                    vec_of_substances.len()
                ); // "The number of columns is not equal to the number of elements",
                return Err(ThermodynamicsError::MatrixDimensionMismatch(format!(
                    "The number of columns ({}) is not equal to the number of elements ({})",
                    A.ncols(),
                    vec_of_substances.len()
                )));
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
                    return Err(ThermodynamicsError::MatrixDimensionMismatch(format!(
                        "The length of vec of initional concentrations ({}) is not equal to the length of row of element i ({})",
                        vec_of_initional_concentrations.len(),
                        row_of_element_i.len()
                    )));
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

        self.create_nonlinear_system_sym().unwrap();
        self.set_T_to_sym();
        self.create_nonlinear_system_fun().unwrap();
        self.create_sum_of_mole_numbers_sym().unwrap();
        self.form_full_system_sym().unwrap();
        self.pretty_print_full_system();

        Ok(())
    }
    // pub fn construct_nonlinear_syste(&mut self) -> Result<(), std::io::Error> {}
}
