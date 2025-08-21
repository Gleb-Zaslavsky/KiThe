use super::ClassicalThermodynamics::Thermodynamics;
use crate::Thermodynamics::User_PhaseOrSolution::{CustomSubstance, SubstancesContainer};
use crate::Thermodynamics::User_substances::DataType;
use prettytable::{Cell, Row, Table, row};

use std::vec;

impl Thermodynamics {
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
    pub fn pretty_print_unknowns(&self) {
        println!("___________________UNKNOWN VARIABLES________________________");
        let mut elem_table = Table::new();
        let unknowns = self.solver.all_unknowns.clone();
        if unknowns.len() == 0 {
            return;
        }
        for (i, eq) in unknowns.iter().enumerate() {
            let row = vec![
                Cell::new(&i.to_string()),
                Cell::new(&eq.to_string().clone()),
            ];
            elem_table.add_row(Row::new(row));
        }
        elem_table.printstd();
        println!("_____________________________________________________________");
    }
    /// Print all Lagrange equations, composition equations and sum of mole numbers as a table
    pub fn pretty_print_full_system(&self) {
        self.pretty_print_Lagrange_equations().unwrap();
        self.pretty_print_composition_equations().unwrap();
        self.pretty_print_sum_mole_numbers().unwrap();
        self.pretty_print_unknowns();
    }
}
