// module containing basic formulae of chemical thermodynamics
//use crate::Kinetics::molmass::create_elem_composition_matrix;

use crate::Thermodynamics::User_PhaseOrSolution::{CustomSubstance, SubstancesContainer};
use crate::Thermodynamics::User_substances::{DataType, SubsData};
use crate::Thermodynamics::dG_dS::{
    calc_dG_for_one_phase, calculate_Gibbs_fun_one_phase, calculate_Gibbs_sym_one_phase,
    calculate_S_for_one_phase, calculate_S_fun_for_one_phase, calculate_S_sym_for_one_phase,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;

use nalgebra::{DMatrix, DVector};
use prettytable::{Table, row};
use std::collections::HashMap;
use std::fmt;
use std::{f64, vec};

/// strucutre for chemical thermodynamics calculations. Calculaton of Gibbs free energy of a given mixure of substances (numerical result at given T, P,
///  concentration) and symbolic Gibbs free energy of a given mixure of substances (symbolic result at given T, P, concentration).
pub struct Thermodynamics {
    pub vec_of_substances: SubstancesContainer,
    /// map of initial concentration of initial substances
    pub map_of_concentration: HashMap<String, f64>,

    /// stoichiometric matrix of the reactions, where columns are substances and rows are reactions
    pub stoichio_matrix: Option<Vec<Vec<f64>>>,
    /// matrix of elements composition of substances, where columns are elements and rows are substances
    pub elem_composition_matrix: Option<DMatrix<f64>>,
    //
    pub initial_vector_of_elements: Vec<f64>,
    pub T: f64,
    pub P: f64,
    /// INSTANCE OF STRUCTURE SUBSDATA CREATED TO SEARCH SUBSTANCES DATA IN THE DATABASES, STORE SEARCH RESULTS AND CALCULATE THERMO PROPERTIES
    pub subdata: CustomSubstance,
    /// hashmap of entropy of a given mixure of substances (numerical result at given T, P, concentration)
    pub dS: HashMap<Option<String>, HashMap<String, f64>>,
    /// hashmap of entropy of a given mixure of substances (function at given T, P, concentration)
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
        let mut cloned_instance = Thermodynamics {
            vec_of_substances: self.vec_of_substances.clone(),
            map_of_concentration: self.map_of_concentration.clone(),

            T: self.T,
            P: self.P,
            stoichio_matrix: self.stoichio_matrix.clone(),
            elem_composition_matrix: self.elem_composition_matrix.clone(),
            initial_vector_of_elements: self.initial_vector_of_elements.clone(),
            subdata: self.subdata.clone(),
            dS: self.dS.clone(),
            dS_fun: new_dS,
            dS_sym: self.dS_sym.clone(),
            dG: self.dG.clone(),
            dG_fun: new_dG,
            dG_sym: self.dG_sym.clone(),
        };
        // Clone the functions
        let _ = cloned_instance.calculate_Gibbs_fun(self.T);
        let _ = cloned_instance.calculate_S_fun(self.T);
        cloned_instance
    }
}

impl fmt::Debug for Thermodynamics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Thermodynamics")
            .field("vec_of_substances", &self.vec_of_substances)
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
        Self {
            vec_of_substances: SubstancesContainer::SinglePhase(Vec::new()),
            map_of_concentration: HashMap::new(),
            T: 298.15,
            P: 1e5,
            stoichio_matrix: None,
            elem_composition_matrix: None,
            initial_vector_of_elements: Vec::new(),
            subdata: CustomSubstance::OnePhase(SubsData::new()),
            dG: HashMap::new(),
            dG_fun: HashMap::new(),
            dG_sym: HashMap::new(),
            dS: HashMap::new(),
            dS_fun: HashMap::new(),
            dS_sym: HashMap::new(),
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
                self.vec_of_substances =
                    SubstancesContainer::SinglePhase(subdata.substances.clone());
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let mut phase_substances = HashMap::new();
                for (phase_name, subsdata) in &phase_or_solution.subs_data {
                    if let Some(phase_name) = phase_name {
                        phase_substances.insert(phase_name.clone(), subsdata.substances.clone());
                    }
                }
                self.vec_of_substances = SubstancesContainer::MultiPhase(phase_substances);
            }
        }
    }

    //////////////////////////////////////////dG//////////////////////////////////////////////

    /// Function for calculating Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration)
    pub fn calculate_Gibbs_free_energy(&mut self, T: f64, n: Option<Vec<f64>>, Np: Option<f64>) {
        self.T = T;
        let sd = &mut self.subdata;
        match sd {
            CustomSubstance::OnePhase(subdata) => {
                let P = self.P;
                let map_to_insert = calc_dG_for_one_phase(P, subdata, T, n, Np);
                self.dG.insert(None, map_to_insert);
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let P = self.P;
                let mut map_to_insert = HashMap::new();
                for (phase_name, subsdata) in phase_or_solution.subs_data.iter_mut() {
                    let map_to_insert_phase = calc_dG_for_one_phase(P, subsdata, T, n.clone(), Np);
                    map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
                }
                self.dG = map_to_insert;
            } // CustomSubstance::PhaseOrSolution
        } // match sd
    }

    /// Symbolic function for calculating Gibbs free energy of a given mixure of substances (symbolic function at given T, P, concentration)
    pub fn calculate_Gibbs_sym(&mut self, T: f64, n: Option<Vec<Expr>>, Np: Option<Expr>) {
        self.T = T;
        let sd = &mut self.subdata;
        match sd {
            CustomSubstance::OnePhase(subdata) => {
                let map_to_insert = calculate_Gibbs_sym_one_phase(subdata, T, n, Np);
                self.dG_sym.insert(None, map_to_insert);
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let mut map_to_insert = HashMap::new();
                for (phase_name, subsdata) in phase_or_solution.subs_data.iter_mut() {
                    let map_to_insert_phase =
                        calculate_Gibbs_sym_one_phase(subsdata, T, n.clone(), Np.clone());
                    map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
                }
                self.dG_sym = map_to_insert;
            } // CustomSubstance::PhaseOrSolution
        } // match sd
    }

    /// set P to all symbolic dG functions
    pub fn set_P_to_sym(&mut self) {
        let P = self.P;
        for (_phase_or_solution_name, sym_fun) in self.dG_sym.iter_mut() {
            for (_, sym_fun) in sym_fun.iter_mut() {
                *sym_fun = sym_fun.set_variable("P", P).simplify()
            }
        }
    }

    pub fn calculate_Gibbs_fun(&mut self, T: f64) {
        self.T = T;
        let sd = &mut self.subdata;
        match sd {
            CustomSubstance::OnePhase(subdata) => {
                let P = self.P;
                let map_to_insert = calculate_Gibbs_fun_one_phase(subdata, P, T);
                self.dG_fun.insert(None, map_to_insert);
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let P = self.P;
                let mut map_to_insert = HashMap::new();
                for (phase_name, subsdata) in phase_or_solution.subs_data.iter_mut() {
                    let map_to_insert_phase = calculate_Gibbs_fun_one_phase(subsdata, P, T);
                    map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
                }
                self.dG_fun = map_to_insert;
            } // CustomSubstance::PhaseOrSolution
        } // match sd
    }
    ////////////////////////////////////////////S - ENTHROPY///////////////////////////////////////////////////////////

    pub fn calculate_S(&mut self, T: f64, n: Option<Vec<f64>>, Np: Option<f64>) {
        self.T = T;
        let sd = &mut self.subdata;
        match sd {
            CustomSubstance::OnePhase(subdata) => {
                let P = self.P;
                let map_to_insert = calculate_S_for_one_phase(P, subdata, T, n, Np);
                self.dS.insert(None, map_to_insert);
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let P = self.P;
                let mut map_to_insert = HashMap::new();
                for (phase_name, subsdata) in phase_or_solution.subs_data.iter_mut() {
                    let map_to_insert_phase =
                        calculate_S_for_one_phase(P, subsdata, T, n.clone(), Np);
                    map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
                }
                self.dS = map_to_insert;
            } // CustomSubstance::PhaseOrSolution
        } // match sd
    }

    /// Function for calculating enthropy of a given mixure of substances ( symbolic result at given T, P, concentration)
    pub fn calculate_S_sym(&mut self, T: f64, n: Option<Vec<Expr>>, Np: Option<Expr>) {
        self.T = T;
        let sd = &mut self.subdata;
        match sd {
            CustomSubstance::OnePhase(subdata) => {
                let map_to_insert = calculate_S_sym_for_one_phase(subdata, T, n, Np);
                self.dS_sym.insert(None, map_to_insert);
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let mut map_to_insert = HashMap::new();
                for (phase_name, subsdata) in phase_or_solution.subs_data.iter_mut() {
                    let map_to_insert_phase =
                        calculate_S_sym_for_one_phase(subsdata, T, n.clone(), Np.clone());
                    map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
                }
                self.dS_sym = map_to_insert;
            } // CustomSubstance::PhaseOrSolution
        } // match sd
    }

    /// Function for calculating enthropy of a given mixure of substances (numerical result at given T, P, concentration)
    pub fn calculate_S_fun(&mut self, T: f64) {
        self.T = T;
        let sd = &mut self.subdata;
        match sd {
            CustomSubstance::OnePhase(subdata) => {
                let P = self.P;
                let map_to_insert = calculate_S_fun_for_one_phase(subdata, P, T);
                self.dS_fun.insert(None, map_to_insert);
            }
            CustomSubstance::PhaseOrSolution(phase_or_solution) => {
                let P = self.P;
                let mut map_to_insert = HashMap::new();
                for (phase_name, subsdata) in phase_or_solution.subs_data.iter_mut() {
                    let map_to_insert_phase = calculate_S_fun_for_one_phase(subsdata, P, T);
                    map_to_insert.insert(phase_name.clone(), map_to_insert_phase);
                }
                self.dS_fun = map_to_insert;
            } // CustomSubstance::PhaseOrSolution
        } // match sd
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Functions for initial composition
    pub fn initial_composition(&mut self) -> Result<(), std::io::Error> {
        let vec_of_substances = match &self.vec_of_substances {
            SubstancesContainer::SinglePhase(vec) => vec.clone(),
            SubstancesContainer::MultiPhase(map) => {
                map.iter().map(|(substance, _)| substance.clone()).collect()
            }
        };

        if let Some(A) = self.elem_composition_matrix.clone() {
            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "The number of substances is not equal to the number of elements",
                ));
            }
            let vec_of_initional_concentrations = self
                .map_of_concentration
                .clone()
                .iter()
                .map(|(_, v)| *v)
                .collect::<Vec<f64>>();
            //let vec_of_initional_concentrations = DVector::from_vec(vec_of_initional_concentrations);
            let mut initial_vector_of_elements: Vec<f64> = Vec::new();
            for i in 0..n_elements {
                let row_of_element_i = A.row(i).iter().map(|&v| v).collect::<Vec<f64>>();
                let row_of_element_i: DVector<_> = DVector::from_vec(row_of_element_i);

                let vec_of_initional_concentrations =
                    DVector::from_vec(vec_of_initional_concentrations.clone());
                let b_i_vec = row_of_element_i * vec_of_initional_concentrations.clone();
                let b_i = b_i_vec.sum();
                initial_vector_of_elements.push(b_i);
            }
            self.initial_vector_of_elements = initial_vector_of_elements;
        }
        Ok(())
    }
    /// Function for conservation of  the elements  condition
    pub fn composition_equations(&mut self) -> Result<(), std::io::Error> {
        if let Some(A) = self.elem_composition_matrix.clone() {
            let vec_of_substances = match &self.vec_of_substances {
                SubstancesContainer::SinglePhase(vec) => vec.clone(),
                SubstancesContainer::MultiPhase(map) => {
                    map.iter().map(|(substance, _)| substance.clone()).collect()
                }
            };

            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "The number of substances is not equal to the number of elements",
                ));
            }

            let mut vector_of_conditions: Vec<Box<dyn FnOnce(DVector<f64>) -> f64>> = Vec::new();
            for i in 0..n_elements {
                let initial_vector_of_elements = self.initial_vector_of_elements.clone();
                let row_of_element_i = A.row(i).iter().map(|&v| v).collect::<Vec<f64>>();
                let row_of_element_i: DVector<_> = DVector::from_vec(row_of_element_i);
                let condition_closure = move |vec_of_concentrations: DVector<f64>| {
                    let b_i_vec = row_of_element_i * vec_of_concentrations.clone();
                    let b_i = b_i_vec.sum();
                    let condition_i = b_i - initial_vector_of_elements[i].clone();
                    condition_i
                };
                vector_of_conditions.push(Box::new(condition_closure));
            }
        }
        Ok(())
    }
    pub fn composition_equation_sym(
        &mut self,
        vec_of_concentrations_sym: Vec<Expr>,
    ) -> Result<(), std::io::Error> {
        if let Some(A) = self.elem_composition_matrix.clone() {
            let vec_of_substances = match &self.vec_of_substances {
                SubstancesContainer::SinglePhase(vec) => vec.clone(),
                SubstancesContainer::MultiPhase(map) => {
                    map.iter().map(|(substance, _)| substance.clone()).collect()
                }
            };

            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
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

                let condition_i = b_i_sum - initial_vector_of_element_sym[i].clone();
                vector_of_conditions.push(condition_i);
            }
        }
        Ok(())
    }
    /// Lagrange multipliers

    ////////////////////////INPUT/OUTPUT////////////////////////////////////////////////////////

    /// Prints the data of the substance to the console
    pub fn pretty_print(&self) -> Result<(), std::io::Error> {
        println!("__________subs properties at {} K__________", self.T);
        let mut table = Table::new();

        let vec_of_substances = match &self.vec_of_substances {
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
                for substance in &vec_of_substances {
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
                        for substance in &vec_of_substances {
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases};

    use approx::assert_relative_eq;
    fn setup_test_data() -> Thermodynamics {
        let subs = vec!["CO".to_string(), "CO2".to_string()];
        let mut subdata = SubsData::new();
        subdata.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        subdata
            .map_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));

        subdata
            .map_of_phases
            .insert("CO2".to_string(), Some(Phases::Gas));
        subdata.substances = subs.clone();
        subdata.search_substances();

        let mut thermo = Thermodynamics::new();
        // Set up basic parameters
        thermo.set_T(400.0);
        thermo.set_P(101325.0, None);
        // Add test substances
        thermo.vec_of_substances = SubstancesContainer::SinglePhase(subs.clone());
        // Set up phases

        // Set up SubsData with mock thermodynamic data
        thermo.subdata = CustomSubstance::OnePhase(subdata);
        thermo
    }

    #[test]
    fn test_new() {
        let thermo = Thermodynamics::new();
        let vec_of_substances =
            if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                vec.clone()
            } else {
                panic!("vec_of_substances is not a SinglePhase");
            };

        assert!(vec_of_substances.is_empty());
        let subdata = if let CustomSubstance::OnePhase(subdata) = &thermo.subdata {
            subdata
        } else {
            panic!("subdata is not a OnePhase");
        };
        assert!(subdata.map_of_phases.is_empty());
        assert_eq!(thermo.T, 298.15);
        assert_eq!(thermo.P, 1e5);
        assert!(thermo.dG.is_empty());
        assert!(thermo.dG_fun.is_empty());
        assert!(thermo.dG_sym.is_empty());
    }

    #[test]
    fn test_set_t_and_p() {
        let mut thermo = Thermodynamics::new();

        // Test setting temperature
        thermo.set_T(500.0);
        assert_eq!(thermo.T, 500.0);

        // Test setting pressure
        thermo.set_P(2e5, None);
        assert_eq!(thermo.P, 2e5);
    }

    #[test]
    fn test_calculate_gibbs_free_energy() {
        let mut thermo = setup_test_data();

        // Calculate Gibbs free energy without concentration correction
        thermo.calculate_Gibbs_free_energy(400.0, None, None);
        let vec_of_substances =
            if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                vec.clone()
            } else {
                panic!("vec_of_substances is not a SinglePhase");
            };
        // Check results
        for substance in &vec_of_substances {
            assert!(thermo.dG.get(&None).unwrap().contains_key(substance));
        }

        // Test with concentration correction
        let concentrations = vec![0.8, 0.2]; // Mock mole fractions
        thermo.calculate_Gibbs_free_energy(400.0, Some(concentrations.clone()), None);
        println!("thermo.dG: {:?} \n", thermo.dG);
        // Check results with gas correction
        for (_i, substance) in vec_of_substances.iter().enumerate() {
            // dG = dH - T*dS + RT*ln(P/P°) + RT*ln(w_i)
            // let gas_correction = R * 400.0 * f64::ln(thermo.P / 101325.0) + R * 400.0 * f64::ln(concentrations[i]);
            // let expected_dg = -110.0 - 400.0 * 200.0 + gas_correction;
            //assert!((thermo.dG[substance] - expected_dg).abs() < 1e-10);
            assert!(thermo.dG.get(&None).unwrap().contains_key(substance));
        }
    }

    #[test]
    fn test_calculate_gibbs_sym() {
        let mut thermo = setup_test_data();
        let vec_of_substances =
            if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                vec.clone()
            } else {
                panic!("vec_of_substances is not a SinglePhase");
            };
        // Calculate symbolic Gibbs free energy
        thermo.calculate_Gibbs_sym(400.0, None, None);
        // Calculate Gibbs free energy without concentration correction
        thermo.calculate_Gibbs_free_energy(400.0, None, None);
        // Check results
        for substance in &vec_of_substances {
            assert!(thermo.dG_sym.get(&None).unwrap().contains_key(substance));

            // Verify the structure of the symbolic expression
            // For our test data, it should be a simple expression: -110 - T*200
            let expr = &thermo.dG_sym.get(&None).unwrap()[substance];
            let dG_function = expr.lambdify1D();
            let dG_i_from_sym = dG_function(400.0);
            let dG_direct = thermo.dG.get(&None).unwrap()[substance];
            assert_relative_eq!(dG_i_from_sym, dG_direct, epsilon = 1e-10);

            // This is a simplified check - in a real test you might want to evaluate the expression
            // or check its structure more thoroughly
            //   assert!(expr.to_string().contains("-110"));
        }

        // Test with concentration correction
        let concentrations = vec![Expr::Var("w1".to_owned()), Expr::Var("w2".to_owned())];
        thermo.calculate_Gibbs_sym(400.0, Some(concentrations), None);

        // Check that the expressions now include concentration terms
        for substance in &vec_of_substances {
            let expr = &thermo.dG_sym.get(&None).unwrap()[substance];
            let expr_str = expr.to_string();

            // The expression should now include ln terms for the gas correction
            assert!(expr_str.contains("ln"));
        }
    }

    #[test]
    fn test_calculate_gibbs_fun() {
        let mut thermo = setup_test_data();
        let vec_of_substances =
            if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                vec.clone()
            } else {
                panic!("vec_of_substances is not a SinglePhase");
            };
        // Calculate Gibbs free energy functions
        thermo.calculate_Gibbs_fun(500.0);
        thermo.calculate_Gibbs_free_energy(500.0, None, None);
        // Check that functions were created for each substance
        for substance in &vec_of_substances {
            assert!(thermo.dG_fun.get(&None).unwrap().contains_key(substance));

            // Test the function at a specific temperature
            let dg_fun = &thermo.dG_fun.get(&None).unwrap()[substance];
            let dg_value = dg_fun(500.0, None, None);
            let dG_direct = thermo.dG.get(&None).unwrap()[substance];
            assert_relative_eq!(dg_value, dG_direct, epsilon = 1e-10);
            // For our test data with constant values:

            // Test with concentration correction
            let dg_value_with_conc = dg_fun(500.0, Some(vec![0.5, 0.5]), None);

            // The value should now include the gas correction term
            assert!(dg_value_with_conc != dg_value);
        }
    }

    #[test]
    fn test_calculate_therm_map_of_fun_local() {
        use crate::Thermodynamics::dG_dS::calculate_therm_map_of_fun_local;
        let thermo = setup_test_data();
        let vec_of_substances =
            if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                vec.clone()
            } else {
                panic!("vec_of_substances is not a SinglePhase");
            };
        let sd = thermo.subdata;
        let mut subdata = if let CustomSubstance::OnePhase(subdata) = sd {
            subdata
        } else {
            panic!()
        };
        let _ = subdata.extract_all_thermal_coeffs(400.0);
        subdata.print_search_summary();

        let result = calculate_therm_map_of_fun_local(&mut subdata);

        assert!(result.is_ok());
        let (cp_funs, dh_funs, ds_funs) = result.unwrap();

        // Check that we got the expected number of functions
        assert_eq!(cp_funs.len(), vec_of_substances.len());
        assert_eq!(dh_funs.len(), vec_of_substances.len());
        assert_eq!(ds_funs.len(), vec_of_substances.len());

        // Test the functions
        for i in 0..vec_of_substances.len() {
            let cp_fun = cp_funs[i].as_ref().unwrap();
            let dh_fun = dh_funs[i].as_ref().unwrap();
            let ds_fun = ds_funs[i].as_ref().unwrap();

            // For our test data, these should return constant values
            assert!(cp_fun(400.0) > 0.0);
            assert!(dh_fun(400.0) < 0.0);
            assert!(ds_fun(400.0) > 0.0);
        }
    }
    #[test]
    fn test_all_features() {
        let subs = vec!["CO".to_string(), "CO2".to_string()];
        // calling instance of strucutre SubsData created to search substances data in the databases, store
        // search results and calculate thermo properties
        let mut subdata = SubsData::new();
        // Set up library priorities
        subdata.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        subdata.substances = subs.clone();
        // Perform the search
        subdata.search_substances();

        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));

        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        // Calling instance of structure Thermodynamics to calculate thermo dG
        let mut thermo = Thermodynamics::new();
        // Set up basic parameters
        thermo.set_T(400.0);
        thermo.set_P(101325.0, None);
        // Add test substances
        thermo.vec_of_substances = SubstancesContainer::SinglePhase(subs.clone());
        // Set up phases

        // savng  the search results in the structure Thermodynamics
        thermo.subdata = CustomSubstance::OnePhase(subdata.clone());
        // Calculate Gibbs free energy
        thermo.calculate_Gibbs_free_energy(400.0, None, None);
        // Calculate symbolic Gibbs free energy
        thermo.calculate_Gibbs_sym(400.0, None, None);
        // getting the results of the calculation
        thermo.calculate_Gibbs_fun(400.0);
        let vec_of_substances =
            if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                vec.clone()
            } else {
                panic!("vec_of_substances is not a SinglePhase");
            };
        let g = thermo.dG.clone();
        let map_of_gibbs = g.get(&None).unwrap();
        let g_sym = thermo.dG_sym.clone();
        let map_of_gibbs_sym = g_sym.get(&None).unwrap();
        let g_fun = thermo.clone().dG_fun;
        let map_of_gibbs_fun = g_fun.get(&None).unwrap();
        for substance in &vec_of_substances {
            println!("substance: {:?}", substance);
            println!("map_of_gibbs: {:?}", map_of_gibbs[substance]);
            println!("map_of_gibbs_sym: {:?}", map_of_gibbs_sym[substance]);
        }
        for substance in &vec_of_substances {
            let dG_value = map_of_gibbs[substance];
            let dG_sym = map_of_gibbs_sym[substance].clone();
            let dG_from_sym = dG_sym.lambdify1D()(400.0);
            assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
            let dG_from_fun = map_of_gibbs_fun.get(substance).unwrap()(400.0, None, None);
            assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);
        }

        println!("map_of_gibbs: {:?} \n", map_of_gibbs);
        let _ = thermo.pretty_print();
    }

    #[test]
    fn test_all_features_witn_conc() {
        // calculating Gibbs free energy withoun concentration correction RT*ln(P/P°) + RT*ln(w_i)
        let subs = vec!["CO".to_string(), "CO2".to_string()];
        // calling instance of strucutre SubsData created to search substances data in the databases, store
        // search results and calculate thermo properties
        let mut subdata = SubsData::new();
        // Set up library priorities
        subdata.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        // Set up phases
        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();
        // Perform the search
        subdata.search_substances();
        // Calling instance of structure Thermodynamics to calculate thermo dG
        let mut thermo = Thermodynamics::new();
        // Set up basic parameters
        thermo.set_T(400.0);
        thermo.set_P(101325.0, None);
        let concentration = Some(vec![0.5, 0.5]);
        // Add test substances
        thermo.vec_of_substances = SubstancesContainer::SinglePhase(subs.clone());

        // savng  the search results in the structure Thermodynamics
        thermo.subdata = CustomSubstance::OnePhase(subdata);
        // Calculate Gibbs free energy
        thermo.calculate_Gibbs_free_energy(400.0, concentration.clone(), None);
        // Calculate symbolic Gibbs free energy
        thermo.calculate_Gibbs_sym(
            400.0,
            Some(vec![
                Expr::Var("w1".to_string()),
                Expr::Var("w2".to_string()),
            ]),
            None,
        );
        thermo.set_P_to_sym();
        // getting the results of the calculation
        thermo.calculate_Gibbs_fun(400.0);
        let map_of_gibbs = thermo.dG.clone().get(&None).unwrap().clone();
        let map_of_gibbs_sym = thermo.dG_sym.get(&None).unwrap().clone();
        let binding = thermo.clone();
        let map_of_gibbs_fun = binding.dG_fun.get(&None).unwrap();
        let vec_of_substances =
            if let SubstancesContainer::SinglePhase(vec) = &thermo.vec_of_substances {
                vec.clone()
            } else {
                panic!("vec_of_substances is not a SinglePhase");
            };
        for substance in &vec_of_substances {
            println!("substance: {:?}", substance);
            println!("map_of_gibbs: {:?}", map_of_gibbs[substance]);
            println!("map_of_gibbs_sym: {:?}", map_of_gibbs_sym[substance]);
        }
        for substance in &vec_of_substances {
            let dG_value = map_of_gibbs[substance];

            let dG_from_fun =
                map_of_gibbs_fun.get(substance).unwrap()(400.0, concentration.clone(), None);
            assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);

            let dG_sym = map_of_gibbs_sym[substance].clone();
            let dG_from_sym = dG_sym.lambdify_borrowed_thread_safe(&["T", "w1", "w2"]);
            let dG_from_sym = dG_from_sym(&[400.0, 0.5, 0.5]);
            assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
        }

        //println!("map_of_gibbs: {:?} \n", map_of_gibbs);
        // let _ = thermo.pretty_print();
    }
}
