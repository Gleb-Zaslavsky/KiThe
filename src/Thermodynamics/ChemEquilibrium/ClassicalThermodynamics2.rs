use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::Thermodynamics;
use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::ThermodynamicsError;
use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamicsSolver::SolverType;
use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamicsSolver::{LMParams, NRParams};
use RustedSciThe::Utils::plots::plots;
use RustedSciThe::numerical::optimization::inter_n_extrapolate::{
    lagrange_interpolate, newton_interpolate,
};
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;
use std::time::Instant;

#[derive(Clone, Copy)]
pub enum InterpolationMethod {
    Lagrange,
    Newton,
}

impl InterpolationMethod {
    fn interpolate(&self, x: f64, x_vals: &DVector<f64>, y_vals: &DVector<f64>) -> f64 {
        match self {
            InterpolationMethod::Lagrange => lagrange_interpolate(x, x_vals, y_vals),
            InterpolationMethod::Newton => newton_interpolate(x, x_vals, y_vals),
        }
    }
}
//use rayon::prelude::*;
//use std::sync::Mutex;
impl Thermodynamics {
    pub fn calculate_equilibrium_for_T_range(
        &mut self,
        nr_params: LMParams,
        non_zero_moles_number: Option<
            HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        >,
        T0: f64,
        Tend: f64,
        step: f64,
    ) -> (HashMap<String, Vec<f64>>, Vec<f64>) {
        let now = Instant::now();
        self.set_number_of_moles(non_zero_moles_number).unwrap();
        self.initial_composition().unwrap();
        self.create_indexed_variables();
        self.create_full_system_of_equations().unwrap();

        self.set_solver_type(SolverType::LevenbergMarquardt(nr_params));
        //  self.set_solver_type(SolverType::NewtonRaphson(nr_params));
        self.generate_eqs();
        let mut empty_map_of_subs = self
            .vec_of_subs
            .clone()
            .into_iter()
            .map(|x| (x, vec![]))
            .collect::<HashMap<String, Vec<f64>>>();
        let mut T = T0;
        let mut T_vec = Vec::new();
        while T <= Tend {
            let map_of_substances: HashMap<String, f64> = self.solve_for_T(T);
            let vec_of_sol = self.solver.get_vector_of_solution();
            self.solver.initial_guess = Some(vec_of_sol);
            // Push values from map_of_substances into empty_map_of_subs
            for (key, value) in &map_of_substances {
                if let Some(vec) = empty_map_of_subs.get_mut(key) {
                    vec.push(*value);
                }
            }
            T_vec.push(T);
            T += step;
        }

        let elapsed = now.elapsed().as_secs();
        println!("Elapsed: {:.2} seconds", elapsed);
        (empty_map_of_subs, T_vec)
    }

    /// Updates the nonlinear system equations when thermodynamic coefficients change
    ///
    /// This method regenerates the system equations after coefficient updates,
    /// typically called when temperature-dependent coefficients become invalid
    pub fn update_system(&mut self) -> Result<(), ThermodynamicsError> {
        // Create symbolic nonlinear system equations
        self.create_nonlinear_system_sym()?;
        // Create function nonlinear system equations
        self.create_nonlinear_system_fun()?;
        // Form the complete symbolic system
        self.form_full_system_sym()?;
        // Generate equations for the solver
        self.generate_eqs();
        Ok(())
    }

    /// Calculates chemical equilibrium over a temperature range with coefficient validation
    ///
    /// This improved version checks for outdated thermodynamic coefficients at each temperature
    /// and updates the system equations when necessary, ensuring accuracy across wide temperature ranges.
    ///
    /// # Arguments
    /// * `nr_params` - Levenberg-Marquardt solver parameters
    /// * `non_zero_moles_number` - Initial mole numbers for substances
    /// * `T0` - Starting temperature [K]
    /// * `Tend` - Ending temperature [K]
    /// * `step` - Temperature step size [K]
    ///
    /// # Returns
    /// Tuple of (substance concentrations vs T, temperature vector)
    pub fn calculate_equilibrium_for_T_range2(
        &mut self,
        nr_params: LMParams,
        non_zero_moles_number: Option<
            HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
        >,
        T0: f64,
        Tend: f64,
        step: f64,
    ) -> Result<(HashMap<String, Vec<f64>>, Vec<f64>), ThermodynamicsError> {
        let now = Instant::now();
        // Set initial mole numbers for the system
        self.set_number_of_moles(non_zero_moles_number)?;
        // Calculate initial elemental composition
        self.initial_composition()?;
        // Create symbolic variables for equations
        self.create_indexed_variables();
        // Convert pressure to symbolic form
        self.set_P_to_sym();
        // Create composition constraint equations
        self.composition_equations()?;
        // Create symbolic composition equations
        self.composition_equation_sym()?;
        // Create sum of mole numbers equations
        self.create_sum_of_mole_numbers_sym()?;
        // Set variable bounds for solver
        self.create_variable_bounds();
        // Configure Levenberg-Marquardt solver
        self.set_solver_type(SolverType::LevenbergMarquardt(nr_params));
        let mut empty_map_of_subs = self
            .vec_of_subs
            .clone()
            .into_iter()
            .map(|x| (x, vec![]))
            .collect::<HashMap<String, Vec<f64>>>();
        let mut T = T0;
        let mut T_vec = Vec::new();
        while T <= Tend {
            // Check if thermodynamic polinomials coefficients are valid at current temperature
            let subs_with_outdated_coeffs = self.extract_coeffs_if_current_coeffs_not_valid(T);

            let subs_with_outdated_coeffs = subs_with_outdated_coeffs?;
            // Update system equations if coefficients changed
            if subs_with_outdated_coeffs.len() > 0 {
                self.update_system()?;
            }

            // Solve equilibrium at current temperature
            let map_of_substances: HashMap<String, f64> = self.solve_for_T(T);
            // Use current solution as initial guess for next temperature
            let vec_of_sol = self.solver.get_vector_of_solution();
            self.solver.initial_guess = Some(vec_of_sol);
            // Store results for each substance
            for (key, value) in &map_of_substances {
                if let Some(vec) = empty_map_of_subs.get_mut(key) {
                    vec.push(*value);
                }
            }
            T_vec.push(T);
            T += step;
        }

        let elapsed = now.elapsed().as_secs();
        println!("Elapsed: {:.2} seconds", elapsed);
        Ok((empty_map_of_subs, T_vec))
    }

    pub fn interpolate_data(
        &self,
        x_vals: Vec<f64>,
        y_map: HashMap<String, Vec<f64>>,
        num_points: usize,
        method: InterpolationMethod,
    ) -> (Vec<f64>, HashMap<String, Vec<f64>>) {
        let x_min = x_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = x_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let step = (x_max - x_min) / (num_points - 1) as f64;

        let new_x: Vec<f64> = (0..num_points).map(|i| x_min + i as f64 * step).collect();

        let x_dv = DVector::from_vec(x_vals);
        let mut new_y_map = HashMap::new();

        for (key, y_vals) in y_map {
            let y_min = y_vals.iter().cloned().fold(f64::INFINITY, f64::min);
            let y_max = y_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let y_dv = DVector::from_vec(y_vals);

            let interpolated: Vec<f64> = new_x
                .iter()
                .map(|&x| {
                    let val = method.interpolate(x, &x_dv, &y_dv);
                    // Clamp values to reasonable bounds for mole fractions
                    if y_max <= 1.0 && y_min >= 0.0 {
                        val.max(0.0).min(1.0)
                    } else {
                        val.max(y_min * 0.1).min(y_max * 10.0)
                    }
                })
                .collect();
            new_y_map.insert(key, interpolated);
        }

        (new_x, new_y_map)
    }

    pub fn plot(&self, x: Vec<f64>, y_map: HashMap<String, Vec<f64>>) {
        let subs = self.vec_of_subs.clone();
        let y: DMatrix<f64> = DMatrix::from_columns(
            &subs
                .iter()
                .map(|sub| {
                    DVector::from_vec(
                        y_map
                            .get(sub)
                            .expect(&format!("No data for substance {}", sub))
                            .clone(),
                    )
                })
                .collect::<Vec<DVector<f64>>>(),
        );
        let x_dv = DVector::from_vec(x);
        plots("T".to_string(), subs, x_dv, y);
    }

    pub fn print_table(&self, map_of_substances: HashMap<String, Vec<f64>>, T_vec: Vec<f64>) {
        use prettytable::{Cell, Row, Table};

        let mut table = Table::new();
        let mut headers = vec![Cell::new("T")];
        let mut substance_keys: Vec<String> = map_of_substances.keys().cloned().collect();
        substance_keys.sort();

        for key in &substance_keys {
            headers.push(Cell::new(key));
        }
        table.add_row(Row::new(headers));

        for (i, &temp) in T_vec.iter().enumerate() {
            let mut row = vec![Cell::new(&temp.to_string())];
            for key in &substance_keys {
                let value = map_of_substances.get(key).unwrap()[i];
                row.push(Cell::new(&format!("{:.6e}", value)));
            }
            table.add_row(Row::new(row));
        }

        table.printstd();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::ThermodynamicCalculations;
    use crate::Thermodynamics::User_PhaseOrSolution::{
        CustomSubstance, SubstanceSystemFactory, SubstancesContainer,
    };
    use RustedSciThe::numerical::Nonlinear_systems::NR::Method;
    fn create_N_plus_O_equilibrium(T: f64) -> Thermodynamics {
        let P = 101325.0;
        let gas_subs = vec![
            "N2".to_string(),
            "O2".to_string(),
            "NO2".to_string(),
            "N2O".to_string(),
            "NO".to_string(),
            "O".to_string(),
            "N2O4".to_string(),
            "N".to_string(),
        ];
        let search_in_NIST = false;
        let explicit_search_insructions = None;
        let library_priorities = vec!["NASA_gas".to_string()];
        let permitted_libraries = vec!["NUIG".to_string()];
        let container = SubstancesContainer::SinglePhase(gas_subs.clone());
        let mut customsubs = SubstanceSystemFactory::create_system(
            container,
            None,
            library_priorities,
            permitted_libraries,
            explicit_search_insructions,
            search_in_NIST,
        )
        .unwrap();
        println!("{:#?} \n", customsubs);

        match &customsubs {
            CustomSubstance::OnePhase(one_phase) => {
                let subs = &one_phase.subs_data.substances;
                assert_eq!(subs.len(), 8);
            }
            _ => panic!(),
        }

        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(Some(T), P, None, None);
        assert!(td.is_ok());

        let td = td.unwrap();

        td
    }
    //  //cargo test tests::test_calculate_equilibrium_for_T_range --release -- --nocapture
    #[test]
    fn test_calculate_equilibrium_for_T_range() {
        let mut td = create_N_plus_O_equilibrium(600.0);
        td.Tm = 1e5;
        /*
                let params = NRParams {
            initial_guess: None,
            tolerance: 1e-6,
            method: Some(Method::damped),
            loglevel: None //Some("info".to_string()),
        };
         */
        let params = LMParams {
            initial_guess: None,
            tolerance: Some(1e-6),
            max_iterations: Some(10_000),
            loglevel: Some("none".to_string()),
        };
        let mut n = HashMap::new();
        let map_of_gas = HashMap::from([("N2".to_string(), 1.0), ("O2".to_string(), 1.0)]);

        n.insert(Some("gas".to_string()), (Some(1.0), Some(map_of_gas)));
        let T0 = 500.0;
        let Tend = 4500.0;
        let step = 20.0;
        let (map_of_substances, T_vec) = td
            .calculate_equilibrium_for_T_range2(params, Some(n), T0, Tend, step)
            .unwrap();

        assert!(map_of_substances.len() > 0);
        //  let ( T_vec, map_of_substances)  =  td.interpolate_data(T_vec, map_of_substances, 1000, InterpolationMethod::Lagrange);
        //  println!("{:#?}", map_of_substances);
        td.print_table(map_of_substances.clone(), T_vec.clone());
        //    td.plot(T_vec, map_of_substances);
    }
}
