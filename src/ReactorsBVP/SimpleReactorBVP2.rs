//! # Pretty Printing Module for Reactor BVP Tasks
//!
//! This module provides formatted output methods for displaying reactor task data,
//! including operating conditions, boundary conditions, kinetic data, and equations.
//! All output is formatted using prettytable for clear tabular presentation.

use super::SimpleReactorBVP::SimpleReactorTask;

impl SimpleReactorTask {
    /// Displays a comprehensive summary of the reactor task including all parameters,
    /// conditions, and kinetic data in formatted tables.
    ///
    /// # Output Sections
    /// - Basic problem information
    /// - Operating conditions (P, T, Cp, λ, m, L, M, Pe_q)
    /// - Scaling parameters
    /// - Boundary conditions
    /// - Transport properties (diffusion coefficients, Peclet numbers)
    /// - Thermal effects
    /// - Substances and reactions
    /// - Molar masses
    /// - Symbolic rate constants
    pub fn pretty_print_task(&self) {
        use prettytable::{Table, row};

        println!("\n=== REACTOR TASK SUMMARY ===");

        // Basic info
        println!("Problem Name: {:?}", self.problem_name);
        println!("Problem Description: {:?}", self.problem_description);

        // Operating conditions
        let mut table = Table::new();
        table.add_row(row!["Parameter", "Value", "Units"]);
        table.add_row(row!["Pressure (P)", format!("{:.2e}", self.P), "Pa"]);
        table.add_row(row!["Temperature (Tm)", format!("{:.2}", self.Tm), "K"]);
        table.add_row(row![
            "Heat Capacity (Cp)",
            format!("{:.2e}", self.Cp),
            "J/kg/K"
        ]);
        table.add_row(row![
            "Thermal Conductivity (Lambda)",
            format!("{:.6}", self.Lambda),
            "W/m/K"
        ]);
        table.add_row(row!["Mass Flow (m)", format!("{:.6}", self.m), "kg/s"]);
        table.add_row(row!["Length Scale (L)", format!("{:.6}", self.L), "m"]);
        table.add_row(row![
            "Mean Molar Mass (M)",
            format!("{:.3}", self.M),
            "g/mol"
        ]);
        table.add_row(row!["Peclet Heat (Pe_q)", format!("{:.7}", self.Pe_q), "-"]);

        // Add vector data
        if !self.Pe_D.is_empty() {
            let pe_d_str = self
                .Pe_D
                .iter()
                .map(|x| format!("{:.2e}", x))
                .collect::<Vec<_>>()
                .join(", ");
            table.add_row(row!["Peclet Mass (Pe_D)", format!("[{}]", pe_d_str), "-"]);
        }

        if !self.thermal_effects.is_empty() {
            let thermal_str = self
                .thermal_effects
                .iter()
                .map(|x| format!("{:.2e}", x))
                .collect::<Vec<_>>()
                .join(", ");
            table.add_row(row![
                "Thermal Effects",
                format!("[{}]", thermal_str),
                "J/mol"
            ]);
        }

        println!("\nOperating Conditions:");
        table.printstd();
        // Scaling
        println!("\nScaling Parameters:");
        println!("  dT: {:.3}", self.scaling.dT);
        println!("  L: {:.7}", self.scaling.L);
        println!("T_scaling: {}", self.T_scaling);

        // Boundary conditions table
        if !self.boundary_condition.is_empty() {
            println!("\nBoundary Conditions:");
            let mut bc_table = Table::new();
            bc_table.add_row(row!["Variable", "Value"]);
            for (key, value) in &self.boundary_condition {
                bc_table.add_row(row![key, format!("{:.6}", value)]);
            }
            bc_table.printstd();
        }

        // Diffusion coefficients table
        if !self.Diffusion.is_empty() {
            println!("\nDiffusion Coefficients:");
            let mut diff_table = Table::new();
            diff_table.add_row(row!["Substance", "D (m²/s)"]);
            for (substance, coeff) in &self.Diffusion {
                diff_table.add_row(row![substance, format!("{:.2e}", coeff)]);
            }
            diff_table.printstd();
        }

        // Transport coefficients table
        if !self.D_ro_map.is_empty() {
            println!("\nTransport Coefficients (D*ρ):");
            let mut transport_table = Table::new();
            transport_table.add_row(row!["Substance", "D*ρ (kg/m/s)"]);
            for (substance, coeff) in &self.D_ro_map {
                transport_table.add_row(row![substance, format!("{:.2e}", coeff)]);
            }
            transport_table.printstd();
        }

        // Peclet numbers table
        if !self.Pe_D.is_empty() {
            println!("\nPeclet Numbers (Diffusion):");
            let mut pe_table = Table::new();
            pe_table.add_row(row!["Substance", "Pe_D"]);
            for (i, pe_d) in self.Pe_D.iter().enumerate() {
                if let Some(substance) = self.kindata.substances.get(i) {
                    pe_table.add_row(row![substance, format!("{:.2e}", pe_d)]);
                }
            }
            pe_table.printstd();
        }

        // Thermal effects table
        if !self.thermal_effects.is_empty() {
            println!("\nThermal Effects:");
            let mut thermal_table = Table::new();
            thermal_table.add_row(row!["Reaction", "Q (J/mol)"]);
            for (i, q) in self.thermal_effects.iter().enumerate() {
                thermal_table.add_row(row![format!("Reaction {}", i + 1), format!("{:.2e}", q)]);
            }
            thermal_table.printstd();
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

        // Molar masses table
        if let Some(mol_masses) = &self.kindata.stecheodata.vec_of_molmasses {
            println!("\nMolar Masses:");
            let mut mass_table = Table::new();
            mass_table.add_row(row!["Substance", "Molar Mass (g/mol)"]);
            for (i, mass) in mol_masses.iter().enumerate() {
                if let Some(substance) = self.kindata.substances.get(i) {
                    mass_table.add_row(row![substance, format!("{:.3}", mass)]);
                }
            }
            mass_table.printstd();
        }

        // Symbolic rate constants table
        if let Some(k_sym_vec) = &self.kindata.K_sym_vec {
            println!("\nSymbolic Rate Constants ({}):", k_sym_vec.len());
            let mut rate_table = Table::new();
            rate_table.add_row(row!["Reaction", "Rate Constant"]);
            for (i, k_sym) in k_sym_vec.iter().enumerate() {
                if let Some(equation) = self.kindata.vec_of_equations.get(i) {
                    rate_table.add_row(row![equation, format!("{}", k_sym)]);
                }
            }
            rate_table.printstd();
        }

        println!("\n=== END TASK SUMMARY ===\n");
    }

    /// Displays the system of differential equations in tabular format.
    ///
    /// Shows the balance equations with their corresponding unknown variables
    /// and symbolic expressions. Each row represents one equation in the BVP system.
    ///
    /// # Table Columns
    /// - **Balance**: Type of balance equation (mass, energy)
    /// - **Unknown Var**: Variable being solved for
    /// - **Equation**: Symbolic mathematical expression
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

    /// Displays reaction rate expressions in tabular format.
    ///
    /// Shows each chemical reaction with its corresponding rate expression,
    /// including kinetic constants and concentration dependencies.
    ///
    /// # Table Columns
    /// - **Reaction**: Chemical reaction equation
    /// - **Rate**: Symbolic rate expression with kinetic parameters
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
