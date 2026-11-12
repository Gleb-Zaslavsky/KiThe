//! # Reactor BVP Post-Processing and Analysis Module
//!
//! ## Aim and General Description
//! This module provides comprehensive post-processing, validation, and analysis capabilities for
//! 1D steady-state plug-flow/combustion reactor boundary value problems (BVP). It handles the
//! conversion of dimensionless solver results back to physical units, performs rigorous energy
//! and mass balance checks, and provides various output formats for visualization and data export.
//!
//! The module is designed to work with chemical reactor simulations involving:
//! - Temperature profiles along reactor length
//! - Species concentration distributions  
//! - Heat flux calculations
//! - Chemical reaction kinetics validation
//!
//! ## Main Structures, Methods and Functions
//!
//! ### Core Implementation (`SimpleReactorTask` impl block)
//! - **`check_balances()`** - Master validation method calling energy and material balance checks
//! - **`check_energy_balance()`** - Validates energy conservation using integral heat balance equations
//! - **`check_material_balance()`** - Verifies atomic mass conservation throughout the reactor
//! - **`postprocessing()`** - Converts dimensionless solver variables back to physical units
//! - **`estimate_values()`** - Provides quick estimates for adiabatic reactor temperature rise
//!
//! ### Conversion Methods
//! - **`from_mass_to_molar_fractions()`** - Converts mass fractions to molar fractions using molecular weights
//! - **`from_mass_fractions_to_molar_concentration()`** - Converts mass fractions to molar concentrations
//! - **`get_only_concentrations()`** - Extracts concentration data from full solution matrix
//!
//! ### I/O and Visualization
//! - **`plot()`, `gnuplot()`, `plot_in_terminal()`** - Various plotting backends for solution visualization
//! - **`save_to_file()`, `save_to_csv()`** - Data export in different formats
//! - **`pretty_print_balances()`** - Formatted output of balance validation results
//!
//! ### Numerical Integration Functions
//! - **`trapezoidal()`** - Trapezoidal rule integration for non-uniform grids
//! - **`simpsons()`** - Adaptive Simpson's rule with fallback to trapezoidal for non-uniform spacing
//! - **`estimate_error_richardson()`** - Richardson extrapolation for integration error estimation
//! - **`estimate_error_simpsons_richardson()`** - Error estimation for Simpson's method
//!
//! ### Supporting Structures
//! - **`AfterSolution`** - Marker struct for post-solution operations
//!
//! ## Interesting Non-Obvious Code Features and Tips
//!
//! ### 1. Dimensionless Variable Scaling
//! The solver works with dimensionless variables for numerical stability. Key scaling relationships:
//! - Temperature: `Teta = (T - T_ref) / T_scale`
//! - Heat flux: `q_dimensionless = q * L / T_scale`  
//! - Mass flux: `J_dimensionless = J * L`
//!
//! ### 2. Adaptive Integration Strategy
//! The `simpsons()` function intelligently switches between Simpson's rule and trapezoidal rule:
//! ```rust, ignore
//! if (h1 - h2).abs() < 1e-10 {
//!     // Uniform spacing - use Simpson's rule
//!     integral += h / 3.0 * (y[i] + 4.0 * y[i + 1] + y[i + 2]);
//! } else {
//!     // Non-uniform spacing - fallback to trapezoidal
//!     integral += 0.5 * (y[i] + y[i + 1]) * h1;
//! }
//! ```
//!
//! ### 3. Richardson Extrapolation Error Estimation
//! Uses coarse grid (every 2nd point) vs fine grid comparison:
//! - Trapezoidal: `E_fine = (I_coarse - I_fine) / 3` (O(h²) method)
//! - Simpson's: `E_fine = (I_coarse - I_fine) / 15` (O(h⁴) method)
//!
//! ### 4. Robust Balance Checking
//! Energy balance uses integral form: `-∇q - ∫F dx + ṁCₚΔT = 0`
//! where F is the heat release function, avoiding numerical differentiation errors.
//!
//! ### 5. Mass Balance via Elemental Conservation
//! Uses elemental composition matrix multiplication to check atomic conservation:
//! ```rust, ignore
//! let vector_of_elements = matrix_of_elements * concentrations_at_step;
//! ```
//! This is more robust than checking individual species balances.
//!
//! ### 6. Iterator Chaining for Efficient Processing
//! Extensive use of iterator patterns for performance:
//! ```rust, ignore
//! x.windows(2).zip(y.windows(2)).map(|(x_pair, y_pair)| { ... }).sum()
//! ```
//!
//! ### 7. Comprehensive Test Coverage
//! Module includes extensive tests covering edge cases, error conditions, and numerical accuracy
//! validation against analytical solutions.

use super::SimpleReactorBVP::{ReactorError, SimpleReactorTask};
use RustedSciThe::Utils::logger::{save_matrix_to_csv, save_matrix_to_file};
use RustedSciThe::Utils::plots::{plots, plots_gnulot, plots_terminal};
use log::{info, warn};
use nalgebra::{DMatrix, DVector};
pub struct AfterSolution {}

/// Snapshot of balance-check results for display and testing.
#[derive(Debug, Clone, PartialEq)]
pub struct BalanceReport {
    pub energy_balane_error_abs: f64,
    pub energy_balane_error_rel: f64,
    pub sum_of_mass_fractions: Vec<(usize, f64)>,
    pub atomic_mass_balance_error: Vec<(usize, f64)>,
}

/// Snapshot of the quick estimation helper.
///
/// Keeping the arithmetic in a data-returning helper makes the UX easier to test
/// and avoids coupling the estimate logic to logging side effects.
#[derive(Debug, Clone, PartialEq)]
pub struct EstimateValuesReport {
    pub reaction_count: usize,
    pub single_reaction_adiabatic_temperature: Option<f64>,
}

/// Owned snapshot of the postprocessed solver state.
///
/// The wrapper method can apply this snapshot back to the solver, while tests can
/// inspect the transformed values without mutating the reactor.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PostprocessingReport {
    pub x_mesh: DVector<f64>,
    pub solution: DMatrix<f64>,
}

/// Owned snapshot of solved reactor data for plotting and export.
///
/// Keeping this in one helper makes the side-effecting methods thin wrappers.
#[derive(Debug, Clone)]
pub(crate) struct SolutionRenderData {
    pub x_mesh: DVector<f64>,
    pub solution: DMatrix<f64>,
    pub unknowns: Vec<String>,
    pub arg_name: String,
}

impl SimpleReactorTask {
    /// Returns the solved matrix if the reactor has already been solved.
    ///
    /// A dedicated helper keeps all runtime checks for solution availability in one place.
    fn solution_ref(&self) -> Result<&DMatrix<f64>, ReactorError> {
        self.solver.solution.as_ref().ok_or_else(|| {
            ReactorError::MissingData(
                "Solver solution is not available; run the BVP solve first".to_string(),
            )
        })
    }

    /// Returns the x-mesh if the solver has already produced one.
    ///
    /// This avoids repeating the same `Option` handling in every post-processing helper.
    fn x_mesh_ref(&self) -> Result<&DVector<f64>, ReactorError> {
        self.solver.x_mesh.as_ref().ok_or_else(|| {
            ReactorError::MissingData(
                "Solver x mesh is not available; run the BVP solve first".to_string(),
            )
        })
    }

    /// Returns the molecular weights used by the concentration conversion helpers.
    ///
    /// The conversion code needs a consistent molar-mass vector and should fail explicitly
    /// instead of assuming the kinetic preprocessing already populated it.
    fn molar_masses_ref(&self) -> Result<&[f64], ReactorError> {
        self.kindata
            .stecheodata
            .vec_of_molmasses
            .as_deref()
            .ok_or_else(|| {
                ReactorError::MissingData(
                    "Molar masses are not available in stecheodata".to_string(),
                )
            })
    }

    /// Returns the index of a solver variable by name.
    ///
    /// Using a typed error here keeps downstream balance code from panicking if a variable name
    /// is missing after a future model change.
    fn unknown_index(&self, name: &str) -> Result<usize, ReactorError> {
        self.solver
            .unknowns
            .iter()
            .position(|value| value == name)
            .ok_or_else(|| {
                ReactorError::MissingData(format!(
                    "Solver variable `{}` is not present in the current state",
                    name
                ))
            })
    }

    /// Performs comprehensive balance validation for the reactor solution.
    ///
    /// Calls both energy and material balance checking methods to validate
    /// the numerical solution against conservation laws.
    pub fn check_balances(&mut self) -> Result<(), ReactorError> {
        self.check_energy_balance()?;
        self.check_material_balance()?;
        Ok(())
    }

    /// Return the current balance metrics as a plain data snapshot.
    ///
    /// This keeps the reporting layer testable without parsing log output.
    pub fn balance_report(&self) -> BalanceReport {
        let quality = &self.solver.quality;
        BalanceReport {
            energy_balane_error_abs: quality.energy_balane_error_abs,
            energy_balane_error_rel: quality.energy_balane_error_rel,
            sum_of_mass_fractions: quality.sum_of_mass_fractions.clone(),
            atomic_mass_balance_error: quality.atomic_mass_balance_error.clone(),
        }
    }

    /// Compute the quick estimate data without printing anything.
    ///
    /// The helper validates the thermal inputs and returns the adiabatic
    /// temperature estimate only for the single-reaction case.
    pub fn estimate_values_report(&self) -> Result<EstimateValuesReport, ReactorError> {
        let reaction_count = self.kindata.vec_of_equations.len();
        let single_reaction_adiabatic_temperature = if reaction_count == 1 {
            let q = *self.thermal_effects.first().ok_or_else(|| {
                ReactorError::MissingData(
                    "Single-reaction estimate requires one thermal effect value".to_string(),
                )
            })?;
            if !q.is_finite() {
                return Err(ReactorError::InvalidNumericValue(
                    "Thermal effect must be finite for quick estimates".to_string(),
                ));
            }

            let cp = self.Cp;
            if !cp.is_finite() || cp <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(
                    "Heat capacity must be finite and positive for quick estimates".to_string(),
                ));
            }

            let t0 = *self.boundary_condition.get("T").ok_or_else(|| {
                ReactorError::MissingData("Boundary condition does not contain `T`".to_string())
            })?;
            if !t0.is_finite() {
                return Err(ReactorError::InvalidNumericValue(
                    "Boundary temperature must be finite for quick estimates".to_string(),
                ));
            }

            Some(t0 + q / cp)
        } else {
            None
        };

        Ok(EstimateValuesReport {
            reaction_count,
            single_reaction_adiabatic_temperature,
        })
    }

    /// Build an owned snapshot of the solved state for plotting and export.
    ///
    /// This keeps all solution availability checks in one place and lets tests
    /// validate the snapshot without invoking file or plotting backends.
    pub(crate) fn solution_render_data(&self) -> Result<SolutionRenderData, ReactorError> {
        Ok(SolutionRenderData {
            x_mesh: self.x_mesh_ref()?.clone(),
            solution: self.solution_ref()?.clone(),
            unknowns: self.solver.unknowns.clone(),
            arg_name: self.solver.arg_name.clone(),
        })
    }

    /// Build a postprocessed snapshot without mutating the solver state.
    ///
    /// This helper keeps the scaling rules in one place and lets tests check the
    /// transformed values before they are written back to the solver snapshot.
    pub(crate) fn postprocessing_report(&self) -> Result<PostprocessingReport, ReactorError> {
        let mut x_mesh = self.x_mesh_ref()?.clone();
        let mut solution = self.solution_ref()?.clone();
        let unknowns = &self.solver.unknowns;

        if unknowns.len() != solution.ncols() {
            return Err(ReactorError::CalculationError(format!(
                "Solution has {} columns but {} unknowns are tracked",
                solution.ncols(),
                unknowns.len()
            )));
        }

        let dT = self.scaling.dT;
        let T_scale = self.scaling.T_scale;
        let L = self.L;

        x_mesh.iter_mut().for_each(|xi| *xi *= L);
        for (i, mut sol_for_var) in solution.column_iter_mut().enumerate() {
            if unknowns[i] == "Teta" {
                sol_for_var
                    .iter_mut()
                    .for_each(|Teta_i| *Teta_i = *Teta_i * T_scale + dT);
            }
            if unknowns[i] == "q" {
                sol_for_var
                    .iter_mut()
                    .for_each(|q_i| *q_i = *q_i * T_scale / L);
            }
            if unknowns[i].starts_with("J") {
                sol_for_var.iter_mut().for_each(|J_i| *J_i = *J_i / L);
            }
        }

        Ok(PostprocessingReport { x_mesh, solution })
    }
    /// Validates energy conservation using integral heat balance equations.
    ///
    /// Checks the energy balance equation: -∇q - ∫F dx + ṁCₚΔT = 0
    /// where F is the heat release function. Calculates both absolute and
    /// relative errors and stores them in solver quality metrics.
    pub fn check_energy_balance(&mut self) -> Result<(), ReactorError> {
        let L = self.L;
        let T_scale = self.scaling.T_scale;

        let solution = self.solution_ref()?;
        let q_index = self.unknown_index("q")?;
        let q_profile: Vec<f64> = solution.column(q_index).iter().copied().collect();
        // return from dimensionless heat flow to dimension
        let q_f = q_profile.last().ok_or_else(|| {
            ReactorError::CalculationError("Heat flux column is empty".to_string())
        })? * T_scale
            / L;
        let q_0 = q_profile.first().ok_or_else(|| {
            ReactorError::CalculationError("Heat flux column is empty".to_string())
        })? * T_scale
            / L;
        // Lambda (dT/dx)_last_ploint - Lambda (dT/dx)_first_ploint
        let dq = q_f - q_0;
        let t_index = self.unknown_index("Teta")?;
        let T_profile: Vec<f64> = solution.column(t_index).iter().copied().collect();
        let dT = self.scaling.dT;
        // return from dimensionless T to dimension
        let T_f = T_profile.last().ok_or_else(|| {
            ReactorError::CalculationError("Temperature column is empty".to_string())
        })? * T_scale
            + dT;
        let T_0 = T_profile.first().ok_or_else(|| {
            ReactorError::CalculationError("Temperature column is empty".to_string())
        })? * T_scale
            + dT;
        let m = self.m;
        let Cp = self.Cp;
        let dT = m * Cp * (T_f - T_0);
        // Keep the solver snapshot borrowed so the balance check stays cheap.
        let unknowns: Vec<&str> = self
            .solver
            .unknowns
            .iter()
            .map(|unknown| unknown.as_str())
            .collect();
        let heat_release = &self.heat_release;
        // compare heat release profiles calculated by 2 methods
        // via lambdify
        let heat_release_fun = heat_release.lambdify_borrowed_thread_safe(unknowns.as_slice());
        let mut heat_releas_val_via_lambdify = Vec::new();
        for solution_for_timestep in solution.row_iter() {
            let solution_for_timestep = solution_for_timestep.iter().cloned().collect::<Vec<f64>>();
            heat_releas_val_via_lambdify.push(heat_release_fun(solution_for_timestep.as_slice()));
        }
        // via evaluation
        let mut heat_releas_val = Vec::new();
        for solution_for_timestep in solution.row_iter() {
            let solution_for_timestep = solution_for_timestep.iter().cloned().collect::<Vec<f64>>();
            let heat_release_for_timestep =
                heat_release.eval_expression(unknowns.as_slice(), &solution_for_timestep);
            heat_releas_val.push(heat_release_for_timestep);
        }
        let heat_rel_2_methods_difference = DVector::from_vec(heat_releas_val_via_lambdify.clone())
            - DVector::from_vec(heat_releas_val.clone());
        let heat_rel_2_methods_difference_error = heat_rel_2_methods_difference.norm();
        // return from dimensionless x to dimension
        let x_mesh: Vec<f64> = self.x_mesh_ref()?.iter().map(|x| x * L).collect();
        // Q = Integral_0_f (F)dx where F - is heat release function - source of the energy balance eq
        let (total_heat_release, heat_release_error) =
            estimate_error_simpsons_richardson(&x_mesh, &heat_releas_val)?;
        let heat_release_integration_rel_error =
            100.0 * (heat_release_error / total_heat_release).abs();
        let (integrated_q, error_integrated_q) =
            estimate_error_simpsons_richardson(&x_mesh, &q_profile)?;
        let integrated_q = integrated_q * T_scale / L;
        let error_integrated_q = error_integrated_q * T_scale / L;
        let integrated_q_rel_error = 100.0 * (error_integrated_q / integrated_q).abs();
        let T_discrete_error = integrated_q / self.Lambda - (T_f - T_0);
        // integral form of heat balance equation
        // - (Lambda (dT/dx)_f - Lambda (dT/dx)_0) -  Integral_0_f (F)dx  +  m * Cp * (T_f - T_0);
        let absolute_energy_error = -dq - total_heat_release + dT;
        let rel_energy_error = 100.0 * absolute_energy_error / total_heat_release;
        self.solver.quality.energy_balane_error_abs = absolute_energy_error;
        self.solver.quality.energy_balane_error_rel = rel_energy_error;
        info!(
            "\n \n relative error in heat balance is {:.3} %",
            rel_energy_error
        );
        if rel_energy_error > 5.0 {
            warn!(
                "ATTENTION! Relative error in heat balance is {:.3} %",
                rel_energy_error
            )
        }
        info!(
            "heat release integration error {:.2} %",
            heat_release_integration_rel_error
        );
        info!(
            "heat_rel_2_methods_difference_error {:}",
            heat_rel_2_methods_difference_error
        );
        info!(
            "T_discrete_error {:.2} % with integration error {:.3}",
            100.0 * T_discrete_error / (T_f - T_0),
            integrated_q_rel_error
        );
        Ok(())
    }

    /// Extracts concentration variables from the full solution matrix.
    ///
    /// Filters the solution matrix to return only columns corresponding to
    /// concentration variables (those starting with "C").
    ///
    /// # Returns
    /// Matrix containing only concentration data with dimensions (n_points, n_species)
    fn get_only_concentrations(&self) -> Result<DMatrix<f64>, ReactorError> {
        let solution = self.solution_ref()?;
        let mut concentration_indices = Vec::new();
        for (i, unknown) in self.solver.unknowns.iter().enumerate() {
            if unknown.starts_with("C") {
                concentration_indices.push(i);
            }
        }
        let nrows = solution.nrows();
        let ncols = concentration_indices.len();
        let mut only_concentrations = DMatrix::zeros(nrows, ncols);
        for (j, &col_idx) in concentration_indices.iter().enumerate() {
            only_concentrations.set_column(j, &solution.column(col_idx));
        }
        Ok(only_concentrations)
    }
    /// Converts mass fractions to molar fractions using molecular weights.
    ///
    /// Performs the conversion: x_i = (w_i/M_i) / Σ(w_j/M_j)
    /// where w_i is mass fraction, M_i is molecular weight, x_i is molar fraction.
    ///
    /// # Arguments
    /// * `matrix_of_mass_fractions` - Matrix of mass fractions (n_points × n_species)
    ///
    /// # Returns
    /// Matrix of molar fractions with same dimensions
    pub fn from_mass_to_molar_fractions(
        &self,
        matrix_of_mass_fractions: DMatrix<f64>,
    ) -> Result<DMatrix<f64>, ReactorError> {
        let Mi = self.molar_masses_ref()?;
        if Mi.iter().any(|m| !m.is_finite() || *m <= 0.0) {
            return Err(ReactorError::InvalidNumericValue(
                "Molar masses must be finite and positive for fraction conversion".to_string(),
            ));
        }
        let Mi: Vec<f64> = Mi.iter().map(|m| *m / 1000.0).collect();
        let Mi = DVector::from_vec(Mi).transpose();
        let mut matrix_of_molar_fractions = DMatrix::zeros(
            matrix_of_mass_fractions.nrows(),
            matrix_of_mass_fractions.ncols(),
        );
        for (i, row) in matrix_of_mass_fractions.row_iter().enumerate() {
            let Minv = Mi.map(|x| 1.0 / x);
            let M_dot_x = row.dot(&Minv);
            let row_new = row
                .iter()
                .enumerate()
                .map(|(i, x)| (x / Mi[i]) / M_dot_x)
                .collect::<Vec<f64>>();
            let row_new = DVector::from_vec(row_new).transpose();
            matrix_of_molar_fractions.set_row(i, &row_new);
        }
        Ok(matrix_of_molar_fractions)
    }

    /// Converts mass fractions to molar concentrations.
    ///
    /// Performs the conversion: C_i = w_i / M_i
    /// where w_i is mass fraction, M_i is molecular weight, C_i is molar concentration.
    ///
    /// # Arguments
    /// * `matrix_of_mass_fractions` - Matrix of mass fractions (n_points × n_species)
    ///
    /// # Returns
    /// Matrix of molar concentrations with same dimensions
    pub fn from_mass_fractions_to_molar_concentration(
        &self,
        matrix_of_mass_fractions: DMatrix<f64>,
    ) -> Result<DMatrix<f64>, ReactorError> {
        let Mi = self.molar_masses_ref()?;
        if Mi.iter().any(|m| !m.is_finite() || *m <= 0.0) {
            return Err(ReactorError::InvalidNumericValue(
                "Molar masses must be finite and positive for concentration conversion".to_string(),
            ));
        }
        let Mi: Vec<f64> = Mi.iter().map(|m| *m / 1000.0).collect();
        let Mi = DVector::from_vec(Mi).transpose();
        let mut matrix_of_molar_fractions = DMatrix::zeros(
            matrix_of_mass_fractions.nrows(),
            matrix_of_mass_fractions.ncols(),
        );
        for (i, row) in matrix_of_mass_fractions.row_iter().enumerate() {
            let row_new = row
                .iter()
                .enumerate()
                .map(|(i, x)| x / Mi[i])
                .collect::<Vec<f64>>();
            let row_new = DVector::from_vec(row_new).transpose();
            matrix_of_molar_fractions.set_row(i, &row_new);
        }
        Ok(matrix_of_molar_fractions)
    }

    /// Backward-compatible wrapper for the legacy misspelled API name.
    ///
    /// The corrected method name should be preferred by new code, but this
    /// wrapper keeps existing callers working while the public API is cleaned up.
    pub fn from_mass_fractions_to_molar_conentration(
        &self,
        matrix_of_mass_fractions: DMatrix<f64>,
    ) -> Result<DMatrix<f64>, ReactorError> {
        self.from_mass_fractions_to_molar_concentration(matrix_of_mass_fractions)
    }
    /// Verifies atomic mass conservation throughout the reactor.
    ///
    /// Checks that elemental composition is conserved by multiplying
    /// the elemental composition matrix with molar concentrations.
    /// Validates that sum of mass fractions equals 1.0 at each point.
    fn check_material_balance(&mut self) -> Result<(), ReactorError> {
        let matrix_of_elements = self
            .kindata
            .stecheodata
            .matrix_of_elements
            .as_ref()
            .ok_or_else(|| {
                ReactorError::MissingData(
                    "Element matrix is not available for material balance checks".to_string(),
                )
            })?
            .transpose();
        // let us check if sum of mass fractions is close to 1 in each step
        let concentrations = self.get_only_concentrations()?;
        if concentrations.ncols() != self.kindata.substances.len() {
            return Err(ReactorError::CalculationError(format!(
                "Concentration matrix has {} columns but {} substances are tracked",
                concentrations.ncols(),
                self.kindata.substances.len()
            )));
        }
        let mut vec_of_sums_of_mass_fractions_at_each_step: Vec<(usize, f64)> = Vec::new();
        for (i, row) in concentrations.row_iter().enumerate() {
            let mut sum = 0.0;
            for (_j, &element) in row.iter().enumerate() {
                sum += element;
            }
            if (sum - 1.0).abs() > 0.01 {
                warn!(
                    "ATTENTION! Sum of concentrations in row {} is {:.3} not 1.0",
                    i, sum
                );
                vec_of_sums_of_mass_fractions_at_each_step.push((i, sum));
            }
        }
        self.solver.quality.sum_of_mass_fractions = vec_of_sums_of_mass_fractions_at_each_step;
        let molar_concentrations = self.from_mass_to_molar_fractions(concentrations)?;
        //    dbg!(&molar_concentrations.column(1).transpose().data.as_slice());

        let initial_concentrations: DVector<f64> = molar_concentrations.row(0).transpose().into();
        let initial_vector_of_elemnts = &matrix_of_elements * initial_concentrations;
        let final_concentrations: DVector<f64> = molar_concentrations
            .row(molar_concentrations.nrows() - 1)
            .transpose()
            .into();
        let final_vector_of_elemnts = &matrix_of_elements * final_concentrations;
        if (&initial_vector_of_elemnts - &final_vector_of_elemnts).norm() > 0.01 {
            warn!("ATTENTION! Initial and final vectors of elements are not the same");
        }
        //(initial_vector_of_elemnts.clone().data.as_slice(), &final_vector_of_elemnts.as_slice());

        let mut atomic_mass_balance_error_vec: Vec<(usize, f64)> = Vec::new();
        for i in 0..molar_concentrations.nrows() {
            let concentrations_at_step: DVector<f64> =
                molar_concentrations.row(i).transpose().into();
            let vector_of_elemnts = &matrix_of_elements * concentrations_at_step;

            let mass_balance_error_for_step =
                (&initial_vector_of_elemnts - &vector_of_elemnts).norm();
            if mass_balance_error_for_step > 0.01 {
                atomic_mass_balance_error_vec.push((i, mass_balance_error_for_step));
                warn!(
                    "ATTENTION! Mass balance error in step {} is {:.3}",
                    i, mass_balance_error_for_step
                );
            }
        }
        self.solver.quality.atomic_mass_balance_error = atomic_mass_balance_error_vec;
        Ok(())
    }
    /////////////////////////////////////////POSTPROCESSING///////////////////////////////////////////
    // solver returns scaled variables - now we must return them from dimensionless to dimension form
    /// Converts dimensionless solver variables back to physical units.
    ///
    /// Applies scaling transformations to return variables from dimensionless
    /// form to their physical dimensions:
    /// - Temperature: T = Teta * T_scale + dT
    /// - Heat flux: q = q_dimensionless * T_scale / L
    /// - Mass flux: J = J_dimensionless / L
    pub fn postprocessing(&mut self) -> Result<(), ReactorError> {
        let report = self.postprocessing_report()?;
        self.solver.x_mesh = Some(report.x_mesh);
        self.solver.solution = Some(report.solution);
        Ok(())
    }
    ////////////////////////////////////////////////I/O/////////////////////////////////////////////////////
    /// Generates plots of the solution using the default plotting backend.
    ///
    /// Creates visualization of all solution variables vs spatial coordinate.
    pub fn plot(&self) -> Result<(), ReactorError> {
        let render_data = self.solution_render_data()?;
        plots(
            render_data.arg_name,
            render_data.unknowns,
            render_data.x_mesh,
            render_data.solution,
        );
        Ok(())
    }
    /// Generates plots using gnuplot backend.
    ///
    /// Creates gnuplot-compatible output files for solution visualization.
    pub fn gnuplot(&self) -> Result<(), ReactorError> {
        let render_data = self.solution_render_data()?;
        plots_gnulot(
            render_data.arg_name,
            render_data.unknowns,
            render_data.x_mesh,
            render_data.solution,
        );
        Ok(())
    }
    /// Displays plots directly in the terminal.
    ///
    /// Creates ASCII-based plots for quick visualization without external tools.
    pub fn plot_in_terminal(&self) -> Result<(), ReactorError> {
        let render_data = self.solution_render_data()?;
        plots_terminal(
            render_data.arg_name,
            render_data.unknowns,
            render_data.x_mesh,
            render_data.solution,
        );
        Ok(())
    }
    /// Saves solution data to a text file.
    ///
    /// # Arguments
    /// * `filename` - Optional filename prefix. Defaults to "result" if None.
    pub fn save_to_file(&self, filename: Option<String>) -> Result<(), ReactorError> {
        let name = if let Some(name) = filename {
            format!("{}.txt", name)
        } else {
            "result.txt".to_string()
        };

        let render_data = self.solution_render_data()?;
        save_matrix_to_file(
            &render_data.solution,
            &render_data.unknowns,
            &name,
            &render_data.x_mesh,
            &render_data.arg_name,
        )
        .map_err(|err| ReactorError::CalculationError(format!("{}", err)))?;
        Ok(())
    }
    /// Saves solution data to CSV format.
    ///
    /// # Arguments
    /// * `filename` - Optional filename prefix. Defaults to "result_table" if None.
    pub fn save_to_csv(&self, filename: Option<String>) -> Result<(), ReactorError> {
        let name = if let Some(name) = filename {
            name
        } else {
            "result_table".to_string()
        };

        let render_data = self.solution_render_data()?;
        save_matrix_to_csv(
            &render_data.solution,
            &render_data.unknowns,
            &name,
            &render_data.x_mesh,
            &render_data.arg_name,
        )
        .map_err(|err| ReactorError::CalculationError(format!("{}", err)))?;
        Ok(())
    }

    /// Provides quick estimates and prints balance validation results.
    ///
    /// For single-reaction systems, estimates adiabatic temperature rise.
    /// Always calls pretty_print_balances() to display validation results.
    pub fn estimate_values(&self) -> Result<(), ReactorError> {
        let report = self.estimate_values_report()?;

        if let Some(t_fin) = report.single_reaction_adiabatic_temperature {
            info!(
                "If the reactor behaves like an ideal mixed adiabatic reactor, the temperature will be {:?}",
                t_fin
            );
        }
        info!(
            "Quick estimate prepared for {} reaction(s)",
            report.reaction_count
        );

        self.pretty_print_balances();
        Ok(())
    }

    /// Displays balance validation results in a formatted table.
    ///
    /// Shows energy balance errors, mass fraction deviations, and
    /// atomic balance violations in a readable tabular format.
    fn pretty_print_balances(&self) {
        use prettytable::{Table, row};
        let report = self.balance_report();

        let mut table = Table::new();
        table.add_row(row!["Parameter", "Value"]);
        table.add_row(row![
            "Absolute energy balance error",
            report.energy_balane_error_abs
        ]);
        table.add_row(row![
            "Relative energy balance error, %",
            report.energy_balane_error_rel
        ]);
        table.add_row(row![
            "Points where mass fractions deviated from 1",
            report.sum_of_mass_fractions.len()
        ]);
        table.add_row(row![
            "Atomic balance violated in points",
            report.atomic_mass_balance_error.len()
        ]);
        info!("{}", table);
    }
}

/// Performs numerical integration using the trapezoidal rule.
///
/// Integrates function values y over non-uniform grid x using the composite
/// trapezoidal rule: ∫f(x)dx ≈ Σ[(x_{i+1} - x_i) * (y_i + y_{i+1})/2]
///
/// # Arguments
/// * `y` - Function values at grid points
/// * `x` - Grid points (must be same length as y)
///
/// # Returns
/// Result containing the integral value or error message
fn trapezoidal(y: &Vec<f64>, x: &Vec<f64>) -> Result<f64, String> {
    if x.len() != y.len() {
        return Err("x and y must have the same length".to_string());
    };
    if x.len() <= 2 {
        return Err("At least 2 points required for integration".to_string());
    };
    if y.len() <= 2 {
        return Err("At least 2 points required for integration".to_string());
    };
    // Zip x and y into pairs, then create an iterator over the segments.
    // For each segment [i, i+1], calculate its area: (x_{i+1} - x_i) * (y_i + y_{i+1}) / 2.0
    let sum: f64 = x
        .windows(2)
        .zip(y.windows(2))
        .map(|(x_pair, y_pair)| {
            let dx = x_pair[1] - x_pair[0];
            let avg_height = (y_pair[0] + y_pair[1]) / 2.0;
            dx * avg_height
        })
        .sum();

    Ok(sum)
}
#[allow(dead_code)]
/// Alternative implementation of trapezoidal integration.
///
/// Simpler loop-based implementation for comparison and testing.
///
/// # Arguments
/// * `heat_values` - Function values at mesh points
/// * `x_mesh` - Spatial mesh points
///
/// # Returns
/// Result containing the integral value or error message
fn trapezoidal2(heat_values: &Vec<f64>, x_mesh: &Vec<f64>) -> Result<f64, String> {
    assert_eq!(
        heat_values.len(),
        x_mesh.len(),
        "Heat values and x_mesh must have same length"
    );

    if heat_values.len() < 2 {
        return Err("At least 2 heat values required for integration".to_string());
    }

    let mut integral = 0.0;
    for i in 0..heat_values.len() - 1 {
        let dx = x_mesh[i + 1] - x_mesh[i];
        integral += 0.5 * (heat_values[i] + heat_values[i + 1]) * dx;
    }
    Ok(integral) // Return the integral
}

/// Performs numerical integration using adaptive Simpson's rule.
///
/// Uses Simpson's rule for uniform spacing, falls back to trapezoidal
/// for non-uniform grids. More accurate than trapezoidal for smooth functions.
///
/// # Arguments
/// * `y` - Function values at grid points
/// * `x` - Grid points (must be same length as y)
///
/// # Returns
/// Result containing the integral value or error message
pub fn simpsons(y: &Vec<f64>, x: &Vec<f64>) -> Result<f64, String> {
    if x.len() != y.len() {
        return Err("x and y must have the same length".to_string());
    }
    if x.len() < 3 {
        return Err("At least 3 points required for Simpson's rule".to_string());
    }

    let mut integral = 0.0;

    // For non-uniform grids, use composite Simpson's rule on pairs of intervals
    let mut i = 0;
    while i + 2 < x.len() {
        let h1 = x[i + 1] - x[i];
        let h2 = x[i + 2] - x[i + 1];

        if (h1 - h2).abs() < 1e-10 {
            // Uniform spacing - use standard Simpson's rule
            let h = h1;
            integral += h / 3.0 * (y[i] + 4.0 * y[i + 1] + y[i + 2]);
            i += 2;
        } else {
            // Non-uniform spacing - use trapezoidal rule for this interval
            integral += 0.5 * (y[i] + y[i + 1]) * h1;
            i += 1;
        }
    }

    // Handle remaining interval if any
    if i + 1 < x.len() {
        let h = x[i + 1] - x[i];
        integral += 0.5 * (y[i] + y[i + 1]) * h;
    }

    Ok(integral)
}
/// Estimates the error of the trapezoidal integration using Richardson Extrapolation.
/// This is done by comparing the result from the full dataset (fine) with the result
/// from a dataset with every second point removed (coarse).
/// Returns a tuple (integral_fine, error_estimate)
#[allow(dead_code)]
fn estimate_error_richardson(x: &Vec<f64>, y: &Vec<f64>) -> Result<(f64, f64), ReactorError> {
    // Calculate the integral with the full (fine) dataset
    let i_fine = trapezoidal(y, x).map_err(ReactorError::CalculationError)?;

    // Create coarse datasets by taking every second point (starting from index 0).
    // This is valid for both even and odd-length vectors.
    let x_coarse: Vec<f64> = x.iter().step_by(2).cloned().collect();
    let y_coarse: Vec<f64> = y.iter().step_by(2).cloned().collect();

    // If the coarse dataset has less than 2 points, Richardson is not possible.
    if x_coarse.len() < 2 {
        return Ok((i_fine, 0.0)); // Can't estimate error, return 0.
    }

    let i_coarse = trapezoidal(&y_coarse, &x_coarse).map_err(ReactorError::CalculationError)?;

    // The error estimate for the FINE solution is (I_coarse - I_fine) / 3.
    // The theory: I_true = I_fine + E_fine, I_true = I_coarse + E_coarse
    // where E_coarse ≈ 4*E_fine (since error is O(h^2) and h_coarse = 2*h_fine)
    // So: I_fine + E_fine = I_coarse + 4*E_fine -> E_fine = (I_coarse - I_fine)/3
    let error_estimate = (i_coarse - i_fine) / 3.0;

    Ok((i_fine, error_estimate))
}

/// Estimates the error of the trapezoidal integration by approximating the second derivative.
/// This method assumes the function is smooth and the data is not too noisy.
/// Returns a tuple (integral, error_estimate)
#[allow(dead_code)]
fn estimate_error_second_derivative(
    x: &Vec<f64>,
    y: &Vec<f64>,
) -> Result<(f64, f64), ReactorError> {
    let integral = trapezoidal(y, x).map_err(ReactorError::CalculationError)?;
    let n = x.len();

    if n < 3 {
        return Ok((integral, 0.0)); // Need at least 3 points to calculate a 2nd derivative.
    }

    // Calculate an average spacing `h_avg` for the entire dataset.
    // This is a heuristic since the mesh is non-uniform.
    let first = *x.first().ok_or_else(|| {
        ReactorError::CalculationError("At least 2 points required for integration".to_string())
    })?;
    let last = *x.last().ok_or_else(|| {
        ReactorError::CalculationError("At least 2 points required for integration".to_string())
    })?;
    let total_interval = last - first;
    let h_avg = total_interval / (n - 1) as f64;

    // Calculate the second finite difference for each interior point i.
    // D2_i = (y_{i-1} - 2y_i + y_{i+1}) / (h_actual^2)
    // Since h is not constant, we use the local spacing around point i.
    // We approximate the local h^2 as ( (x_i - x_{i-1}) * (x_{i+1} - x_i) )
    let mut second_derivatives = Vec::new();

    for i in 1..n - 1 {
        let h_left = x[i] - x[i - 1];
        let h_right = x[i + 1] - x[i];
        // A better heuristic is to use a local average h for the denominator.
        // Using the geometric mean of the adjacent intervals for h^2.
        let local_h_sq = h_left * h_right;
        // Avoid division by zero in case of duplicate nodes (shouldn't happen).
        if local_h_sq > 1e-12 {
            let d2_i = (y[i - 1] - 2.0 * y[i] + y[i + 1]) / local_h_sq;
            second_derivatives.push(d2_i.abs());
        }
    }

    // Find the maximum value of the absolute second derivative approximation.
    let max_f2 = second_derivatives
        .into_iter()
        .fold(0.0_f64, |acc, value| acc.max(value));

    // The error formula for the composite trapezoidal rule is:
    // E = - ( (b-a) * h^2 / 12 ) * f''(ξ)
    // We use h_avg for h and max_f2 as an estimate for |f''(ξ)|.
    // This gives us an UPPER BOUND estimate for the error magnitude.
    let error_estimate = (total_interval * h_avg.powi(2) / 12.0) * max_f2;

    Ok((integral, error_estimate))
}

use std::f64;

/// Estimates the error of Simpson's integration using Richardson Extrapolation
/// Returns a tuple (integral_fine, error_estimate)
/// For Simpson's method with O(h⁴) error: E ≈ (I_fine - I_coarse) / 15
#[allow(dead_code)]
fn estimate_error_simpsons_richardson(
    x: &Vec<f64>,
    y: &Vec<f64>,
) -> Result<(f64, f64), ReactorError> {
    // Calculate the integral with the full (fine) dataset
    let i_fine = simpsons(y, x).map_err(ReactorError::CalculationError)?;

    // Create coarse datasets by taking every second point
    let x_coarse: Vec<f64> = x.iter().step_by(2).cloned().collect();
    let y_coarse: Vec<f64> = y.iter().step_by(2).cloned().collect();

    // Simpson's rule needs at least 3 points
    if x_coarse.len() < 3 {
        return Ok((i_fine, 0.0)); // Can't estimate error meaningfully
    }

    let i_coarse = simpsons(&y_coarse, &x_coarse).map_err(ReactorError::CalculationError)?;

    // Richardson extrapolation for O(h⁴) methods: E = (I_fine - I_coarse) / 15
    let error_estimate = (i_coarse - i_fine) / 15.0;

    Ok((i_fine, error_estimate))
}

#[cfg(test)] // This attribute tells Rust to compile this code only when running tests
mod tests_trapezoide {
    use super::*; // Import all functions from the parent module
    use approx::assert_abs_diff_eq; // For floating-point comparisons with tolerance

    /// Helper function to generate a uniform mesh for testing
    fn generate_uniform_mesh(start: f64, end: f64, n_points: usize) -> Vec<f64> {
        let step = (end - start) / (n_points - 1) as f64;
        (0..n_points).map(|i| start + i as f64 * step).collect()
    }

    /// Helper function to generate a non-uniform mesh for testing
    fn generate_nonuniform_mesh(start: f64, end: f64, n_points: usize) -> Vec<f64> {
        // Create a mesh that is denser near the start and sparser near the end
        let mut x = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let t = i as f64 / (n_points - 1) as f64;
            // Use a quadratic transformation to make it non-uniform
            x.push(start + (end - start) * t * t);
        }
        x
    }

    #[test]
    fn test_trapezoidal_basic() {
        // Test a simple linear function: f(x) = 2x
        // Integral from 0 to 3 is x^2 |_0^3 = 9
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 2.0, 4.0, 6.0];
        let result = trapezoidal(&y, &x).unwrap();
        assert_abs_diff_eq!(result, 9.0, epsilon = 1e-10);
    }

    #[test]
    fn test_trapezoidal_constant() {
        // Test a constant function: f(x) = 5
        // Integral from 1 to 4 is 5 * (4-1) = 15
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = vec![5.0, 5.0, 5.0, 5.0];
        let result = trapezoidal(&y, &x).unwrap();
        assert_abs_diff_eq!(result, 15.0, epsilon = 1e-10);
    }

    #[test]
    fn test_richardson_on_linear_function() {
        // For a linear function, the trapezoidal rule is exact.
        // Therefore, the error estimate should be very close to zero.
        let x = generate_uniform_mesh(0.0, 5.0, 11); // 11 points
        let y = x.iter().map(|x_val| 2.0 * x_val + 3.0).collect(); // f(x)=2x+3

        let (integral, error_estimate) = estimate_error_richardson(&x, &y).unwrap();

        // The integral should be exact: ∫(2x+3)dx from 0 to 5 = (x^2 + 3x)|_0^5 = 25 + 15 = 40
        assert_abs_diff_eq!(integral, 40.0, epsilon = 1e-10);
        // The error estimate should be very small
        assert!(
            error_estimate.abs() < 1e-10,
            "Error estimate for linear function should be near zero, got {}",
            error_estimate
        );
    }

    #[test]
    fn test_second_derivative_on_linear_function() {
        // For a linear function, the second derivative is zero.
        // Therefore, the error estimate should be zero.
        let x = generate_uniform_mesh(0.0, 5.0, 11);
        let y = x.iter().map(|x_val| 2.0 * x_val + 3.0).collect(); // f(x)=2x+3

        let (integral, error_estimate) = estimate_error_second_derivative(&x, &y).unwrap();

        assert_abs_diff_eq!(integral, 40.0, epsilon = 1e-10);
        assert_abs_diff_eq!(error_estimate, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_richardson_on_quadratic_function() {
        // Test with a quadratic function f(x) = x^2.
        // The trapezoidal rule has a known, non-zero error for quadratics.
        let n_points = 21;
        let start = 0.0;
        let end = 4.0;
        let x = generate_uniform_mesh(start, end, n_points);
        let y = x.iter().map(|x_val| x_val * x_val).collect();

        let true_integral = (end.powi(3) - start.powi(3)) / 3.0; // 64/3 ≈ 21.333...

        let (integral, error_estimate) = estimate_error_richardson(&x, &y).unwrap();

        // Check that the actual error (integral - true_integral) and the
        // Richardson estimate are reasonably close. For a smooth function like
        // x^2, Richardson should be very accurate.
        let actual_error = integral - true_integral;
        assert_abs_diff_eq!(error_estimate, actual_error, epsilon = 1e-3);
    }

    #[test]
    fn test_second_derivative_on_quadratic_function() {
        // For f(x) = x^2, f''(x) = 2.0 everywhere.
        let n_points = 21;
        let start = 0.0;
        let end = 4.0;
        let x = generate_uniform_mesh(start, end, n_points);
        let y = x.iter().map(|x_val| x_val * x_val).collect();

        let (integral, error_estimate) = estimate_error_second_derivative(&x, &y).unwrap();
        let true_integral = (end.powi(3) - start.powi(3)) / 3.0; // 64/3

        // The theoretical error bound for uniform mesh is:
        // E = - ( (b-a) * h^2 / 12 ) * f''(ξ)
        let h = (end - start) / (n_points - 1) as f64;
        let theoretical_error_bound = ((end - start) * h.powi(2) / 12.0) * 2.0;

        // The calculated error estimate should be close to the theoretical bound.
        // It might not be exact due to the averaging heuristics for non-uniform meshes,
        // but for a uniform mesh it should be very close.
        assert_abs_diff_eq!(error_estimate, theoretical_error_bound, epsilon = 1e-3);
        // The actual error should be less than this bound
        info!("{}", (integral - true_integral).abs());
        assert!((integral - true_integral).abs() < 10.0 * error_estimate);
    }

    #[test]
    fn test_non_uniform_mesh() {
        // Test that the functions work on a non-uniform mesh
        let x = generate_nonuniform_mesh(0.0, 3.0, 15);
        let y = x.iter().map(|x_val| x_val.sin()).collect();

        let (integral_rich, error_rich) = estimate_error_richardson(&x, &y).unwrap();
        let (integral_sec, error_sec) = estimate_error_second_derivative(&x, &y).unwrap();

        // Both methods should complete without panicking and return values
        assert!(integral_rich.is_finite());
        assert!(error_rich.is_finite());
        assert!(integral_sec.is_finite());
        assert!(error_sec.is_finite());

        // For a smooth function like sin(x), the two integral calculations
        // should be very close to each other.
        assert_abs_diff_eq!(integral_rich, integral_sec, epsilon = 1e-10);
    }

    #[test]
    fn test_richardson_with_too_few_points() {
        // Test behavior with a very small dataset
        let x = vec![0.0, 1.0]; // Only 2 points -> coarse mesh will have 1 point
        let y = vec![1.0, 2.0];

        let result = estimate_error_richardson(&x, &y);
        assert!(result.is_err());
    }

    #[test]
    fn test_second_derivative_with_too_few_points() {
        // Test behavior with a very small dataset
        let x = vec![0.0, 1.0]; // Only 2 points -> can't calculate 2nd derivative
        let y = vec![1.0, 2.0];

        let result = estimate_error_second_derivative(&x, &y);
        assert!(result.is_err());
    }

    #[test]
    fn test_trapezoidal_panic_with_one_point() {
        let x = vec![1.0];
        let y = vec![5.0];
        assert!(trapezoidal(&y, &x).is_err());
    }

    #[test]
    fn test_trapezoidal_panic_with_mismatched_lengths() {
        let x = vec![1.0, 2.0];
        let y = vec![5.0];
        assert!(trapezoidal(&x, &y).is_err());
    }
}
