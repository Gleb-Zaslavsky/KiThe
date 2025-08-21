use super::SimpleReactorBVP::SimpleReactorTask;
use std::collections::HashMap;

impl SimpleReactorTask {
    ////////////////////////////////////////////////I/O/////////////////////////////////////////////////////
    pub fn parse_from_file(&mut self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(file_path)?;
        let sections: Vec<&str> = content.split("\n\n").collect();

        for section in sections {
            let lines: Vec<&str> = section.lines().collect();
            if lines.is_empty() {
                continue;
            }

            match lines[0].trim() {
                "DESCRIPTION" => {
                    let description = lines[1..].join("\n");
                    self.set_problem_description(&description);
                }
                "REACTIONS" => {
                    let mut reaction_data = Vec::new();
                    let mut vector_of_reactions = Vec::new();
                    for line in &lines[1..] {
                        let parts: Vec<&str> = line.split_whitespace().collect();
                        if parts.len() < 5 {
                            continue;
                        }

                        // Parse the line: reaction Q A n Ea
                        let reaction_eq = parts[0].to_string();
                        let Q: f64 = parts[1].parse()?;
                        let A: f64 = parts[2].parse()?;
                        let n: f64 = parts[3].parse()?;
                        let Ea: f64 = parts[4].parse()?;

                        let mut reaction_map = HashMap::new();

                        reaction_map.insert("Q".to_string(), vec![Q]);
                        reaction_map.insert("kin".to_string(), vec![A, n, Ea]);

                        reaction_data.push(reaction_map);
                        vector_of_reactions.push(reaction_eq)
                    }
                }
                _ => continue,
            }
        }
        Ok(())
    }
    ////////////////////////PRETTY PRINTING/////////////////////////////

    pub fn pretty_print_task(&self) {
        use prettytable::{Table, row};

        println!("\n=== REACTOR TASK SUMMARY ===");

        // Basic info
        println!("Problem Name: {:?}", self.problem_name);
        println!("Problem Description: {:?}", self.problem_description);

        // Operating conditions
        let mut table = Table::new();
        table.add_row(row!["Parameter", "Value", "Units"]);
        table.add_row(row!["Pressure (P)", format!("{:.2}", self.P), "Pa"]);
        table.add_row(row!["Temperature (Tm)", format!("{:.2}", self.Tm), "K"]);
        table.add_row(row![
            "Heat Capacity (Cp)",
            format!("{:.2}", self.Cp),
            "J/kg/K"
        ]);
        table.add_row(row![
            "Thermal Conductivity (Lambda)",
            format!("{:.6}", self.Lambda),
            "W/m/K"
        ]);
        table.add_row(row!["Mass Flow (m)", format!("{:.6}", self.m), "kg/s"]);
        table.add_row(row!["Length Scale (L)", format!("{:.3}", self.L), "m"]);
        table.add_row(row![
            "Mean Molar Mass (M)",
            format!("{:.3}", self.M),
            "g/mol"
        ]);
        table.add_row(row!["Peclet Heat (Pe_q)", format!("{:.2}", self.Pe_q), "-"]);
        println!("\nOperating Conditions:");
        table.printstd();

        // Scaling
        println!("\nScaling Parameters:");
        println!("  dT: {:.3}", self.scaling.dT);
        println!("  L: {:.3}", self.scaling.L);
        println!("T_scaling: {}", self.T_scaling);

        // Boundary conditions
        println!("\nBoundary Conditions:");
        for (key, value) in &self.boundary_condition {
            println!("  {}: {:.6}", key, value);
        }

        // Diffusion coefficients
        println!("\nDiffusion Coefficients:");
        for (substance, coeff) in &self.Diffusion {
            println!("  {}: {:.2e}", substance, coeff);
        }

        // Transport coefficients
        println!("\nTransport Coefficients (D*ro):");
        for (substance, coeff) in &self.D_ro_map {
            println!("  {}: {:.2e}", substance, coeff);
        }

        // Peclet numbers for diffusion
        println!("\nPeclet Numbers (Diffusion):");
        for (i, pe_d) in self.Pe_D.iter().enumerate() {
            if let Some(substance) = self.kindata.substances.get(i) {
                println!("  {}: {:.7}", substance, pe_d);
            }
        }

        // Thermal effects
        println!("\nThermal Effects:");
        for (i, q) in self.thermal_effects.iter().enumerate() {
            println!("  Reaction {}: {:.2e} J/mol", i + 1, q);
        }

        // Kinetic data
        println!("\nSubstances ({}):", self.kindata.substances.len());
        for (i, substance) in self.kindata.substances.iter().enumerate() {
            println!("  {}: {}", i + 1, substance);
        }

        println!("\nReactions ({}):", self.kindata.vec_of_equations.len());
        for (i, equation) in self.kindata.vec_of_equations.iter().enumerate() {
            println!("  {}: {}", i + 1, equation);
        }

        // Stoichiometric data
        if let Some(mol_masses) = &self.kindata.stecheodata.vec_of_molmasses {
            println!("\nMolar Masses:");
            for (i, mass) in mol_masses.iter().enumerate() {
                if let Some(substance) = self.kindata.substances.get(i) {
                    println!("  {}: {:.3} g/mol", substance, mass);
                }
            }
        }

        // Symbolic rate constants
        if let Some(k_sym_vec) = &self.kindata.K_sym_vec {
            println!("\nSymbolic Rate Constants ({}):", k_sym_vec.len());
            for (i, k_sym) in k_sym_vec.iter().enumerate() {
                if let Some(equation) = self.kindata.vec_of_equations.get(i) {
                    println!("  {}: {}", equation, k_sym);
                }
            }
        }

        println!("\n=== END TASK SUMMARY ===\n");
    }

    /// Display map_of_equations as a table with columns: "balance", "unknown var", "equation"
    pub fn pretty_print_equations(&self) {
        println!("____________________EQUATIONS_________________________");
        use prettytable::{Cell, Row, Table, row};

        let mut table = Table::new();
        table.add_row(row!["Balance", "Unknown Var", "Equation"]);

        for (balance, (unknown_var, equation)) in &self.map_of_equations {
            table.add_row(Row::new(vec![
                Cell::new(balance),
                Cell::new(unknown_var),
                Cell::new(&format!("{}", equation)),
            ]));
        }

        table.printstd();
    }

    /// Display map_eq_rate as a table with columns: "reaction", "rate"
    pub fn pretty_print_reaction_rates(&self) {
        println!("____________________REACTION RATES_________________________");
        use prettytable::{Cell, Row, Table, row};

        let mut table = Table::new();
        table.add_row(row!["Reaction", "Rate"]);

        for (reaction, rate) in &self.map_eq_rate {
            table.add_row(Row::new(vec![
                Cell::new(reaction),
                Cell::new(&format!("{}", rate)),
            ]));
        }

        table.printstd();
    }
    ////////////////////////////////////////////////////////////////////////
}
