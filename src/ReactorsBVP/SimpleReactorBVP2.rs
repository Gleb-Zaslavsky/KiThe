//! # Pretty Printing Module for Reactor BVP Tasks
//!
//! This module provides formatted output methods for displaying reactor task data,
//! including operating conditions, boundary conditions, kinetic data, and equations.
//! All output is formatted using prettytable for clear tabular presentation.

use super::SimpleReactorBVP::SimpleReactorTask;
use log::info;

/// Canonical snapshot for reactor summary output.
///
/// The reporting layer uses this structure so tests can assert on stable data
/// instead of depending on formatted console output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactorTaskReport {
    pub problem_name: Option<String>,
    pub problem_description: Option<String>,
    pub operating_conditions: Vec<(String, String, String)>,
    pub scaling_parameters: Vec<(String, String)>,
    pub boundary_conditions: Vec<(String, String)>,
    pub diffusion_coefficients: Vec<(String, String)>,
    pub transport_coefficients: Vec<(String, String)>,
    pub peclet_numbers: Vec<(String, String)>,
    pub thermal_effects: Vec<(String, String)>,
    pub substances: Vec<String>,
    pub reactions: Vec<String>,
    pub molar_masses: Vec<(String, String)>,
    pub symbolic_rate_constants: Vec<(String, String)>,
}

impl SimpleReactorTask {
    /// Build a stable summary snapshot for display and testing.
    pub fn task_report(&self) -> ReactorTaskReport {
        let mut boundary_conditions = Vec::with_capacity(self.kindata.substances.len() * 2 + 2);
        if let Some(t) = self.boundary_condition.get("T") {
            boundary_conditions.push(("T".to_string(), format!("{:.6}", t)));
        }
        if let Some(q) = self.boundary_condition.get("q") {
            boundary_conditions.push(("q".to_string(), format!("{:.6}", q)));
        }
        for substance in &self.kindata.substances {
            if let Some(value) = self.boundary_condition.get(substance) {
                boundary_conditions.push((substance.clone(), format!("{:.6}", value)));
            }
            let mut flux_key = String::with_capacity(substance.len() + 5);
            flux_key.push_str(substance);
            flux_key.push_str("_flux");
            if let Some(value) = self.boundary_condition.get(&flux_key) {
                boundary_conditions.push((flux_key, format!("{:.6}", value)));
            }
        }

        let diffusion_coefficients = self
            .kindata
            .substances
            .iter()
            .filter_map(|substance| {
                self.Diffusion
                    .get(substance)
                    .map(|value| (substance.clone(), format!("{:.2e}", value)))
            })
            .collect::<Vec<_>>();

        let transport_coefficients = self
            .kindata
            .substances
            .iter()
            .filter_map(|substance| {
                self.D_ro_map
                    .get(substance)
                    .map(|value| (substance.clone(), format!("{:.2e}", value)))
            })
            .collect::<Vec<_>>();

        let peclet_numbers = self
            .kindata
            .substances
            .iter()
            .enumerate()
            .filter_map(|(index, substance)| {
                self.Pe_D
                    .get(index)
                    .map(|value| (substance.clone(), format!("{:.2e}", value)))
            })
            .collect::<Vec<_>>();

        let thermal_effects = self
            .thermal_effects
            .iter()
            .enumerate()
            .map(|(index, q)| (format!("Reaction {}", index + 1), format!("{:.2e}", q)))
            .collect::<Vec<_>>();

        let molar_masses = self
            .kindata
            .stecheodata
            .vec_of_molmasses
            .as_ref()
            .map(|molar_masses| {
                self.kindata
                    .substances
                    .iter()
                    .enumerate()
                    .filter_map(|(index, substance)| {
                        molar_masses
                            .get(index)
                            .map(|mass| (substance.clone(), format!("{:.3}", mass)))
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();

        let symbolic_rate_constants = self
            .kindata
            .vec_of_equations
            .iter()
            .enumerate()
            .filter_map(|(index, equation)| {
                self.kindata
                    .K_sym_vec
                    .as_ref()
                    .and_then(|rates| rates.get(index))
                    .map(|rate| (equation.clone(), format!("{}", rate)))
            })
            .collect::<Vec<_>>();

        ReactorTaskReport {
            problem_name: self.problem_name.clone(),
            problem_description: self.problem_description.clone(),
            operating_conditions: vec![
                (
                    "Pressure (P)".to_string(),
                    format!("{:.2e}", self.P),
                    "Pa".to_string(),
                ),
                (
                    "Temperature (Tm)".to_string(),
                    format!("{:.2}", self.Tm),
                    "K".to_string(),
                ),
                (
                    "Heat Capacity (Cp)".to_string(),
                    format!("{:.2e}", self.Cp),
                    "J/kg/K".to_string(),
                ),
                (
                    "Thermal Conductivity (Lambda)".to_string(),
                    format!("{:.6}", self.Lambda),
                    "W/m/K".to_string(),
                ),
                (
                    "Mass Flow (m)".to_string(),
                    format!("{:.6}", self.m),
                    "kg/s".to_string(),
                ),
                (
                    "Length Scale (L)".to_string(),
                    format!("{:.6}", self.L),
                    "m".to_string(),
                ),
                (
                    "Mean Molar Mass (M)".to_string(),
                    format!("{:.3}", self.M),
                    "g/mol".to_string(),
                ),
                (
                    "Peclet Heat (Pe_q)".to_string(),
                    format!("{:.7}", self.Pe_q),
                    "-".to_string(),
                ),
            ],
            scaling_parameters: vec![
                ("dT".to_string(), format!("{:.3}", self.scaling.dT)),
                ("L".to_string(), format!("{:.7}", self.scaling.L)),
                ("T_scaling".to_string(), format!("{}", self.T_scaling)),
            ],
            boundary_conditions,
            diffusion_coefficients,
            transport_coefficients,
            peclet_numbers,
            thermal_effects,
            substances: self.kindata.substances.clone(),
            reactions: self.kindata.vec_of_equations.clone(),
            molar_masses,
            symbolic_rate_constants,
        }
    }

    /// Return the equations in solver order instead of hash-map order.
    pub fn equation_report_rows(&self) -> Vec<(String, String, String)> {
        let mut rows = Vec::with_capacity(self.solver.unknowns.len());
        for unknown in &self.solver.unknowns {
            if let Some((balance, (unknown_var, equation))) = self
                .map_of_equations
                .iter()
                .find(|(_, (stored_unknown, _))| stored_unknown == unknown)
            {
                rows.push((
                    balance.clone(),
                    unknown_var.clone(),
                    format!("{}", equation),
                ));
            }
        }
        rows
    }

    /// Return the reaction rate rows in reaction order.
    pub fn reaction_rate_report_rows(&self) -> Vec<(String, String)> {
        self.kindata
            .vec_of_equations
            .iter()
            .enumerate()
            .filter_map(|(index, equation)| {
                self.kindata
                    .K_sym_vec
                    .as_ref()
                    .and_then(|rates| rates.get(index))
                    .map(|rate| (equation.clone(), format!("{}", rate)))
            })
            .collect::<Vec<_>>()
    }

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
        let report = self.task_report();

        info!("\n=== REACTOR TASK SUMMARY ===");

        // Basic info
        info!("Problem Name: {:?}", report.problem_name);
        info!("Problem Description: {:?}", report.problem_description);

        // Operating conditions
        let mut table = Table::new();
        table.add_row(row!["Parameter", "Value", "Units"]);
        for (label, value, units) in &report.operating_conditions {
            table.add_row(row![label, value, units]);
        }

        // Add vector data
        if !report.peclet_numbers.is_empty() {
            let pe_d_str = report
                .peclet_numbers
                .iter()
                .map(|(_, value)| value.clone())
                .collect::<Vec<_>>()
                .join(", ");
            table.add_row(row!["Peclet Mass (Pe_D)", format!("[{}]", pe_d_str), "-"]);
        }

        if !report.thermal_effects.is_empty() {
            let thermal_str = report
                .thermal_effects
                .iter()
                .map(|(_, value)| value.clone())
                .collect::<Vec<_>>()
                .join(", ");
            table.add_row(row![
                "Thermal Effects",
                format!("[{}]", thermal_str),
                "J/mol"
            ]);
        }

        info!("\nOperating Conditions:\n{}", table);
        // Scaling
        info!("\nScaling Parameters:");
        for (label, value) in &report.scaling_parameters {
            info!("  {}: {}", label, value);
        }

        // Boundary conditions table
        if !report.boundary_conditions.is_empty() {
            info!("\nBoundary Conditions:");
            let mut bc_table = Table::new();
            bc_table.add_row(row!["Variable", "Value"]);
            for (key, value) in &report.boundary_conditions {
                bc_table.add_row(row![key, value]);
            }
            info!("{}", bc_table);
        }

        // Diffusion coefficients table
        if !report.diffusion_coefficients.is_empty() {
            info!("\nDiffusion Coefficients:");
            let mut diff_table = Table::new();
            diff_table.add_row(row!["Substance", "D (m²/s)"]);
            for (substance, coeff) in &report.diffusion_coefficients {
                diff_table.add_row(row![substance, coeff]);
            }
            info!("{}", diff_table);
        }

        // Transport coefficients table
        if !report.transport_coefficients.is_empty() {
            info!("\nTransport Coefficients (D*ρ):");
            let mut transport_table = Table::new();
            transport_table.add_row(row!["Substance", "D*ρ (kg/m/s)"]);
            for (substance, coeff) in &report.transport_coefficients {
                transport_table.add_row(row![substance, coeff]);
            }
            info!("{}", transport_table);
        }

        // Peclet numbers table
        if !report.peclet_numbers.is_empty() {
            info!("\nPeclet Numbers (Diffusion):");
            let mut pe_table = Table::new();
            pe_table.add_row(row!["Substance", "Pe_D"]);
            for (substance, pe_d) in &report.peclet_numbers {
                pe_table.add_row(row![substance, pe_d]);
            }
            info!("{}", pe_table);
        }

        // Thermal effects table
        if !report.thermal_effects.is_empty() {
            info!("\nThermal Effects:");
            let mut thermal_table = Table::new();
            thermal_table.add_row(row!["Reaction", "Q (J/mol)"]);
            for (reaction, q) in &report.thermal_effects {
                thermal_table.add_row(row![reaction, q]);
            }
            info!("{}", thermal_table);
        }

        // Kinetic data
        info!("\nSubstances ({}):", report.substances.len());
        for (i, substance) in report.substances.iter().enumerate() {
            info!("  {}: {}", i + 1, substance);
        }

        info!("\nReactions ({}):", report.reactions.len());
        for (i, equation) in report.reactions.iter().enumerate() {
            info!("  {}: {}", i + 1, equation);
        }

        // Molar masses table
        if !report.molar_masses.is_empty() {
            info!("\nMolar Masses:");
            let mut mass_table = Table::new();
            mass_table.add_row(row!["Substance", "Molar Mass (g/mol)"]);
            for (substance, mass) in &report.molar_masses {
                mass_table.add_row(row![substance, mass]);
            }
            info!("{}", mass_table);
        }

        // Symbolic rate constants table
        if !report.symbolic_rate_constants.is_empty() {
            info!(
                "\nSymbolic Rate Constants ({}):",
                report.symbolic_rate_constants.len()
            );
            let mut rate_table = Table::new();
            rate_table.add_row(row!["Reaction", "Rate Constant"]);
            for (equation, k_sym) in &report.symbolic_rate_constants {
                rate_table.add_row(row![equation, k_sym]);
            }
            info!("{}", rate_table);
        }

        info!("\n=== END TASK SUMMARY ===\n");
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
        info!("____________________EQUATIONS_________________________");
        use prettytable::{Cell, Row, Table, row};

        let mut table = Table::new();
        table.add_row(row!["Balance", "Unknown Var", "Equation"]);

        for (balance, unknown_var, equation) in self.equation_report_rows() {
            table.add_row(Row::new(vec![
                Cell::new(&balance),
                Cell::new(&unknown_var),
                Cell::new(&equation),
            ]));
        }

        info!("{}", table);
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
        info!("____________________REACTION RATES_________________________");
        use prettytable::{Cell, Row, Table, row};

        let mut table = Table::new();
        table.add_row(row!["Reaction", "Rate"]);

        for (reaction, rate) in self.reaction_rate_report_rows() {
            table.add_row(Row::new(vec![Cell::new(&reaction), Cell::new(&rate)]));
        }

        info!("{}", table);
    }
    ////////////////////////////////////////////////////////////////////////
}
