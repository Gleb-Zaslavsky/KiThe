//! Chemical equilibrium calculator using log-mole formulation and reaction extent approach.
//!
//! # Purpose
//! This module provides a robust equilibrium calculator that solves chemical equilibrium problems
//! by minimizing Gibbs free energy using reaction extents. It handles both single-temperature
//! and temperature-range calculations for gas-phase and multi-phase systems.
//!
//! # Main Structures
//! - [`EquilibriumLogMoles`]: Main solver structure containing all equilibrium data and methods
//! - [`SolverParams`]: Configuration parameters for numerical solvers
//! - [`Phase`]: Represents different thermodynamic phases (ideal gas, ideal solution)
//! - [`Solvers`]: Enum for choosing between Levenberg-Marquardt and Newton-Raphson solvers
//!
//! # Key Methods
//! - [`EquilibriumLogMoles::solve`]: Solves equilibrium for single temperature
//! - [`EquilibriumLogMoles::solve_for_T_range`]: Solves equilibrium over temperature range
//! - [`EquilibriumLogMoles::create_equilibrium_system`]: Sets up the equilibrium system
//! - [`gas_solver`]: Convenience function for gas-phase equilibrium calculations
//!
//! # Examples
//! ```rust, ignore
//! use KiThe::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq::*;
//!
//! // Create equilibrium solver for gas-phase system
//! let mut solver = gas_solver(
//!     vec!["CO".to_string(), "CO2".to_string(), "O2".to_string()],
//!     1000.0, // Temperature in K
//!     101325.0, // Pressure in Pa
//!     Solvers::LM,
//!     Some("info"),
//!     true // enable scaling
//! );
//!
//! // Set initial moles and solve
//! solver.set_problem(vec![1.0, 0.1, 0.5], phases, 101325.0);
//! solver.solve().unwrap();
//! ```
//!
//! # Non-obvious Solutions and Tips
//! - Uses log-mole variables (y_i = ln(n_i)) to handle species with very small concentrations
//! - Automatic reaction basis generation from elemental composition matrix using SVD
//! - Scaling improves numerical conditioning for systems with vastly different magnitudes
//! - Line search with feasibility constraints prevents negative mole numbers
//! - Temperature-range solving reuses previous solution as initial guess for better convergence
use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq2::{
    PhaseManager, multiphase_equilibrium_residual_generator_sym,
};
use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq3::{
    LMSolver, NRSolver, ReactionBasis, ReactionExtentError, TrustRegionSolver,
    compute_reaction_basis,
};
use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::SubsDataError;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use RustedSciThe::symbolic::symbolic_functions::Jacobian;
use log::{error, info};
use nalgebra::{DMatrix, DVector};
use prettytable::{Cell, Row, Table};
use simplelog::{ColorChoice, Config, LevelFilter, SimpleLogger, TermLogger, TerminalMode};
use std::collections::HashMap;
use std::default::Default;
use std::f64;
use std::rc::Rc;
use std::time::Instant;
//TODO!
// fix scaling
// fix subs disappear
/// Universal gas constant in J/(mol·K)
pub const R: f64 = 8.314;

/// Configuration parameters for numerical solvers
#[derive(Clone)]
pub struct SolverParams {
    /// Maximum number of iterations before giving up
    pub max_iter: usize,
    /// Convergence tolerance for residual norm
    pub tol: f64,
    /// Initial damping parameter for Levenberg-Marquardt
    pub lambda: f64,
    /// Minimum step size in line search before failure
    pub alpha_min: f64,
    // Trust region parameters
    pub delta_init: f64,
    pub delta_max: f64,
    pub eta: f64, // acceptance threshold (e.g. 0.1)
}

/// Available numerical solvers for equilibrium calculations
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Solvers {
    /// Levenberg-Marquardt solver with damping
    LM,
    /// Newton-Raphson solver with line search
    NR,
    //Trust region
    TR,
}
impl Default for SolverParams {
    fn default() -> Self {
        Self {
            max_iter: 50,
            tol: 1e-6,
            lambda: 1e-3,
            alpha_min: 1e-6,
            delta_init: 1.0,
            delta_max: 100.0,
            eta: 0.1,
        }
    }
}
pub struct EquilibriumLogMoles {
    pub subs_data: SubsData,
    pub elem_composition: DMatrix<f64>,
    pub reaction_basis: ReactionBasis,
    pub initial_guess: Option<Vec<f64>>,
    pub P: f64,
    pub T: f64,
    pub stoich_matrix: DMatrix<f64>,
    pub n0: Vec<f64>,
    pub gibbs: Vec<GibbsFn>, // m
    pub gibbs_sym: Vec<Expr>,
    pub phases: Vec<Phase>,
    pub species_phase: Vec<usize>,
    pub list_of_failed_T: Vec<f64>,
    pub solution: Vec<f64>,
    pub moles: Vec<f64>,
    pub elements_vector: Vec<f64>,
    pub loglevel: Option<String>,
    pub scaling_flag: bool,
    pub species_eps: f64,
    pub substate_eps: f64,
    pub map_of_moles_for_each_substance: HashMap<String, Vec<f64>>,
    pub moles_for_T_range: Vec<(f64, Vec<f64>)>,
    pub solver_params: SolverParams,
    pub solver: Solvers,
    pub phase_manager: PhaseManager,
    pub p0: f64,
    tol: f64,
}
impl EquilibriumLogMoles {
    pub fn new() -> Self {
        let ReactionBasis0 = ReactionBasis {
            rank: 0,
            num_reactions: 0,
            reactions: DMatrix::zeros(0, 0),
        };
        Self {
            subs_data: SubsData::new(),
            elem_composition: DMatrix::zeros(0, 0),
            reaction_basis: ReactionBasis0,
            initial_guess: None,
            P: 101325.0,
            T: 273.15,
            stoich_matrix: DMatrix::zeros(0, 0),
            n0: Vec::new(),
            gibbs: Vec::new(),
            gibbs_sym: Vec::new(),
            phases: Vec::new(),
            species_phase: Vec::new(),
            solution: Vec::new(),
            list_of_failed_T: Vec::new(),
            map_of_moles_for_each_substance: HashMap::new(),
            moles_for_T_range: Vec::new(),
            moles: Vec::new(),
            elements_vector: Vec::new(),
            loglevel: None,
            scaling_flag: false,
            species_eps: 1e-30,
            substate_eps: 1e-30,
            solver_params: SolverParams::default(),
            solver: Solvers::LM,
            phase_manager: PhaseManager::default(),
            p0: 101325.0,

            tol: 1e-6,
        }
    }

    pub fn empty() -> Self {
        let ReactionBasis0 = ReactionBasis {
            rank: 0,
            num_reactions: 0,
            reactions: DMatrix::zeros(0, 0),
        };
        Self {
            subs_data: SubsData::empty(),
            elem_composition: DMatrix::zeros(0, 0),
            reaction_basis: ReactionBasis0,
            initial_guess: None,
            P: 101325.0,
            T: 273.15,
            stoich_matrix: DMatrix::zeros(0, 0),
            n0: Vec::new(),
            gibbs: Vec::new(),
            gibbs_sym: Vec::new(),
            phases: Vec::new(),
            species_phase: Vec::new(),
            solution: Vec::new(),
            list_of_failed_T: Vec::new(),
            map_of_moles_for_each_substance: HashMap::new(),
            moles_for_T_range: Vec::new(),
            moles: Vec::new(),
            elements_vector: Vec::new(),
            loglevel: None,
            scaling_flag: false,
            species_eps: 1e-30,
            substate_eps: 1e-30,
            solver_params: SolverParams::default(),
            solver: Solvers::LM,
            phase_manager: PhaseManager::default(),
            p0: 101325.0,

            tol: 1e-6,
        }
    }
    pub fn create_stoich_matrix(&mut self) -> Result<(), ReactionExtentError> {
        let rb = compute_reaction_basis(&self.elem_composition, self.tol)?;
        self.reaction_basis = rb;
        self.stoich_matrix = self.reaction_basis.reactions.clone();
        let elements_vector = compute_element_totals(&self.elem_composition, &self.n0);
        self.elements_vector = elements_vector.data.as_vec().clone();
        let m = self.n0.len();
        let species_phase = species_to_phase_map(&self.phases, m)?;
        self.species_phase = species_phase;
        Ok(())
    }

    /// Construct with optional log level. If `loglevel` is `None` logging is left disabled.
    pub fn with_loglevel(&mut self, loglevel: Option<&str>) -> Self {
        let mut s = Self::new();
        s.loglevel = loglevel.map(|s| s.to_string());
        if let Some(ref lv) = s.loglevel {
            Self::init_logger(lv);
        }
        s
    }

    fn init_logger(level: &str) {
        // parse common level names
        let lf = match level.to_lowercase().as_str() {
            "off" | "none" => LevelFilter::Off,
            "error" => LevelFilter::Error,
            "warn" | "warning" => LevelFilter::Warn,
            "info" => LevelFilter::Info,
            "debug" => LevelFilter::Debug,
            "trace" => LevelFilter::Trace,
            _ => LevelFilter::Info,
        };

        if lf == LevelFilter::Off {
            return;
        }

        // Try SimpleLogger first, fall back to TermLogger if necessary
        let _ = SimpleLogger::init(lf, Config::default()).or_else(|_| {
            TermLogger::init(
                lf,
                Config::default(),
                TerminalMode::Mixed,
                ColorChoice::Auto,
            )
        });
    }
    pub fn set_initial_guess(&mut self, guess: Vec<f64>) {
        self.initial_guess = Some(guess.clone());
    }
    pub fn set_problem(&mut self, n0: Vec<f64>, phases: Vec<Phase>, P: f64) {
        self.n0 = n0;
        self.phases = phases;
        self.P = P;
    }
    pub fn set_n0_from_non_zero_map(&mut self, map_of_nonzero_moles: HashMap<String, f64>) {
        let subs = self.subs_data.substances.clone();
        let mut n0 = vec![1e-8; subs.len()];
        for (i, subs_i) in subs.iter().enumerate() {
            if let Some(moles_i) = map_of_nonzero_moles.get(subs_i) {
                n0[i] = *moles_i;
            }
        }
        self.n0 = n0;
    }
    /// Create the equilibrium system: search substances, parse coeffs, calculate elem comp and molar mass, create stoich matrix
    pub fn create_equilibrium_system(&mut self) -> Result<(), ReactionExtentError> {
        let user_subs = &mut self.subs_data;
        user_subs
            .search_substances()
            .map_err(|e| ReactionExtentError::SubsDataError(e))?;
        user_subs
            .parse_all_thermal_coeffs()
            .map_err(|e| ReactionExtentError::SubsDataError(e))?;
        user_subs
            .calculate_elem_composition_and_molar_mass(None)
            .map_err(|e| ReactionExtentError::SubsDataError(e))?;
        self.elem_composition = user_subs.clone().elem_composition_matrix.unwrap().clone();
        self.create_stoich_matrix()?;

        Ok(())
    }
    pub fn create_symbolic_system(&mut self) -> Result<(), ReactionExtentError> {
        // creating variables
        let stoich = self.stoich_matrix.clone();
        let vars = Expr::IndexedVars(stoich.nrows(), "y").0;
        let vars_Str: Vec<String> = vars.iter().map(|y| y.to_string().clone()).collect();
        //let vars_str: Vec<&str> = vars_Str.iter().map(|yi| yi.as_str()).collect();
        let mut jac = Jacobian::new();
        // creating vector of symbolic dG0
        let user_subs = &mut self.subs_data;
        user_subs
            .calculate_therm_map_of_sym()
            .map_err(|e| ReactionExtentError::SubsDataError(e))?;
        let map_of_functions = user_subs.calculate_dG0_sym_one_phase();
        let mut g_vec: Vec<Expr> = Vec::new();
        for subs in user_subs.substances.clone() {
            let dg0 = map_of_functions.get(&subs).unwrap();
            g_vec.push(dg0.clone());
        }
        // vector of elements
        // let elements_vector = compute_element_totals(&self.elem_composition, &self.n0);
        // self.elements_vector = elements_vector.data.as_vec().clone();
        let elements_vector = self.elements_vector.clone();
        // creating vector of symbolic equations for equilibrium
        let vec_of_functions = multiphase_equilibrium_residual_generator_sym(
            self.stoich_matrix.clone(),
            self.elem_composition.clone(),
            elements_vector.clone(),
            g_vec,
            self.phases.clone(),
            self.P,
            self.p0,
        )?;
        jac.set_vector_of_functions(vec_of_functions);
        jac.vector_of_variables = vars;
        jac.variable_string = vars_Str;
        jac.parameters_string = vec!["T".to_string()];
        jac.calc_jacobian();

        Ok(())
    }
    pub fn solve_for_T_range(
        &mut self,
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Result<Vec<(f64, Vec<f64>)>, ReactionExtentError> {
        let now = Instant::now();
        let mut T = T_start;
        let mut moles_for_T_range = Vec::new();

        while T < T_end {
            info!("\n solving equilibrium for T = {} \n", T);
            // Build Gibbs functions while holding a short-lived mutable borrow to subs_data
            let user_subs = &mut self.subs_data;
            let subs = user_subs.substances.clone();
            let vec_of_coeffs = user_subs
                .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T)
                .map_err(|_| {
                    ReactionExtentError::SubsDataError(SubsDataError::CoefficientExtractionFailed {
                        substance: "all_subs".to_string(),
                        temperature: Some(T),
                    })
                })?;
            if vec_of_coeffs.len() != 0 || self.gibbs.len() == 0 {
                info!("\n coefficients outdated ... \n renewing pipeline ...");
                let g_vec: Vec<GibbsFn> = {
                    let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();

                    let mut local_g_vec: Vec<GibbsFn> = Vec::new();
                    for key in &subs {
                        let boxed_fn = map_of_functions.remove(key).unwrap();
                        // take ownership of the boxed function
                        let boxed_fn_owned = boxed_fn;
                        // move the owned boxed fn into the closure and wrap in Rc
                        let gibbs_fn: GibbsFn = Rc::new(move |T: f64| boxed_fn_owned(T));
                        local_g_vec.push(gibbs_fn);
                    }

                    local_g_vec
                }; // user_subs borrow ends here

                self.gibbs = g_vec;
            }
            self.T = T; // Temperature in K
            match self.solve() {
                Ok(()) => {}
                Err(e) => {
                    error!("Equilibrium calculation failed at T = {} K: {:?}", T, e);
                    self.list_of_failed_T.push(T);
                    T += T_step;
                    continue;
                }
            }

            let xi_sol = &self.solution;
            self.initial_guess = Some(xi_sol.clone());
            let moles = self.moles.clone();

            moles_for_T_range.push((T, moles.clone()));
            //  println!("moles {:?}", moles);
            // println!("T = {}, eta = {:?}", T, xi_sol);
            // info!("Reaction extents solution: {:?}", xi_sol);

            T += T_step;
        }
        self.moles_for_T_range = moles_for_T_range.clone();
        self.map_of_moles_for_each_substance();
        println!(
            "total elapsed time non-parallel solver{} ms",
            now.elapsed().as_millis()
        );
        Ok(Vec::new())
    }

    /// parallel solver
    pub fn solve_for_T_range_par(
        &mut self,
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Result<Vec<(f64, Vec<f64>)>, ReactionExtentError> {
        use rayon::prelude::*;
        use std::sync::Arc;
        let now = Instant::now();

        let subs = self.subs_data.substances.clone();

        // Step 1: Identify temperature ranges where coefficients change (single O(n) pass)
        let mut temp_ranges: Vec<(f64, f64)> = Vec::new();
        let mut T = T_start;

        while T < T_end {
            let user_subs = &mut self.subs_data;
            let vec_of_coeffs = user_subs
                .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T)
                .map_err(|_| {
                    ReactionExtentError::SubsDataError(SubsDataError::CoefficientExtractionFailed {
                        substance: "all_subs".to_string(),
                        temperature: Some(T),
                    })
                })?;

            if vec_of_coeffs.len() != 0 || self.gibbs.len() == 0 {
                let range_start = T;
                let mut T_next = T + T_step;

                // Find next temperature where coefficients change
                while T_next < T_end {
                    if user_subs
                        .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T_next)
                        .unwrap_or_default()
                        .len()
                        != 0
                    {
                        break;
                    }
                    T_next += T_step;
                }

                temp_ranges.push((range_start, T_next));
                T = T_next;
            } else {
                T += T_step;
            }
        }
        let step2 = now.elapsed().as_millis();
        // Step 2: Pre-compute Gibbs functions for each identified range
        let mut gibbs_cache: HashMap<
            usize,
            (Vec<Arc<dyn Fn(f64) -> f64 + Send + Sync>>, f64, f64),
        > = HashMap::new();

        for (range_id, (T_start_range, T_end_range)) in temp_ranges.iter().enumerate() {
            let user_subs = &mut self.subs_data;
            user_subs
                .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(*T_start_range)
                .map_err(|_| {
                    ReactionExtentError::SubsDataError(SubsDataError::CoefficientExtractionFailed {
                        substance: "all_subs".to_string(),
                        temperature: Some(*T_start_range),
                    })
                })?;

            // Create Gibbs functions for this range as Arc (thread-safe)
            let g_vec: Vec<Arc<dyn Fn(f64) -> f64 + Send + Sync>> = {
                let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();
                let mut local_g_vec: Vec<Arc<dyn Fn(f64) -> f64 + Send + Sync>> = Vec::new();
                for key in &subs {
                    let boxed_fn = map_of_functions.remove(key).unwrap();
                    let gibbs_fn: Arc<dyn Fn(f64) -> f64 + Send + Sync> =
                        Arc::new(move |T: f64| boxed_fn(T));
                    local_g_vec.push(gibbs_fn);
                }
                local_g_vec
            };

            gibbs_cache.insert(range_id, (g_vec, *T_start_range, *T_end_range));
        }
        let step3 = now.elapsed().as_millis();
        // Step 2: Generate all temperature points with their range indices
        let mut temp_range_pairs: Vec<(f64, usize)> = Vec::new();
        T = T_start;

        while T < T_end {
            // Find which range this temperature belongs to
            let mut found_range_id = None;
            for (id, (_gibbs_vec, t_start, t_end)) in &gibbs_cache {
                if T >= *t_start && T < *t_end {
                    found_range_id = Some(*id);
                    break;
                }
            }

            if let Some(range_id) = found_range_id {
                temp_range_pairs.push((T, range_id));
            }
            T += T_step;
        }
        let step4 = now.elapsed().as_millis();
        // Step 3: Parallel computation
        let elem_composition = self.elem_composition.clone();
        let stoich_matrix = self.stoich_matrix.clone();
        let n0 = self.n0.clone();
        let phases = self.phases.clone();
        let elements_vector = self.elements_vector.clone();
        let P = self.P;
        let p0 = self.p0;
        let solver_params = self.solver_params.clone();
        let species_phase = self.species_phase.clone();
        let solver = self.solver;
        let scaling_flag = self.scaling_flag;
        let species_eps = self.species_eps;
        let substate_eps = self.substate_eps;
        let gibbs_cache = Arc::new(gibbs_cache);
        println!("time for creation closure {:?}", now.elapsed());
        println!(" {}, {}, {}", step2, step3, step4);
        let now = Instant::now();
        // Chunk the temperature tasks to reduce per-task scheduling overhead.
        // Compute chunk size based on Rayon thread count so work is balanced.
        let num_threads = rayon::current_num_threads();
        let chunk_size = std::cmp::max(1, (temp_range_pairs.len() + num_threads - 1) / num_threads);

        let chunks: Vec<Vec<(f64, usize)>> = temp_range_pairs
            .chunks(chunk_size)
            .map(|c| c.to_vec())
            .collect();

        let mut chunk_results: Vec<
            Vec<Result<(f64, Vec<f64>, Vec<f64>, f64), (f64, ReactionExtentError)>>,
        > = Vec::new();

        for chunk in chunks {
            let initial_guess = self.initial_guess.clone();
            println!("new initial guess {:?}", initial_guess);
            let local_results: Vec<
                Result<(f64, Vec<f64>, Vec<f64>, f64), (f64, ReactionExtentError)>,
            > = chunk
                .par_iter()
                .map(|(T, range_id)| {
                    // Retrieve the Gibbs functions for this range
                    let gibbs_vec = if let Some((g_vec, _, _)) = gibbs_cache.get(range_id) {
                        g_vec.clone()
                    } else {
                        return Err((
                            *T,
                            ReactionExtentError::Other("Gibbs cache not found".to_string()),
                        ));
                    };
                    let now = Instant::now();
                    let mut local_solver = EquilibriumLogMoles::empty();
                    let elapsed = now.elapsed();
                    local_solver.subs_data.substances = subs.clone();
                    local_solver.elem_composition = elem_composition.clone();
                    local_solver.stoich_matrix = stoich_matrix.clone();
                    local_solver.species_phase = species_phase.clone();
                    local_solver.n0 = n0.clone();
                    local_solver.phases = phases.clone();
                    local_solver.elements_vector = elements_vector.clone();
                    local_solver.P = P;
                    local_solver.p0 = p0;
                    local_solver.T = *T;
                    local_solver.solver_params = solver_params.clone();
                    local_solver.solver = solver;
                    local_solver.scaling_flag = scaling_flag;
                    local_solver.species_eps = species_eps;
                    local_solver.substate_eps = substate_eps;
                    local_solver.initial_guess = initial_guess.clone();

                    // Convert Arc-based gibbs functions back to Rc for internal use
                    let gibbs_as_rc: Vec<GibbsFn> = gibbs_vec
                        .iter()
                        .map(|arc_fn| {
                            let fn_clone = arc_fn.clone();
                            Rc::new(move |T: f64| fn_clone(T)) as GibbsFn
                        })
                        .collect();

                    local_solver.gibbs = gibbs_as_rc;

                    match local_solver.solve() {
                        Ok(()) => Ok((
                            *T,
                            local_solver.solution.clone(),
                            local_solver.moles.clone(),
                            elapsed.as_millis() as f64,
                        )),
                        Err(e) => Err((*T, e)),
                    }
                })
                .collect();
            let last_solution = local_results
                .last()
                .unwrap()
                .as_ref()
                .ok()
                .map(|(_, sol, _, _)| sol.clone())
                .unwrap();
            let time = local_results
                .iter()
                .map(|res| match res {
                    Ok((_, _, _, t)) => *t,
                    Err(_) => 0.0,
                })
                .sum::<f64>();
            println!("time for chunk {:?} ms", time);
            self.initial_guess = Some(last_solution);
            // Use last solution of this chunk as initial guess for next chunk
            chunk_results.push(local_results);
        }

        // Flatten chunk results into a single results Vec
        let mut results: Vec<Result<(f64, Vec<f64>, Vec<f64>, f64), (f64, ReactionExtentError)>> =
            chunk_results.into_iter().flatten().collect();

        // Step 4: Process results
        let mut moles_for_T_range = Vec::new();
        let mut failed_temps = Vec::new();

        for result in results {
            match result {
                Ok((T, _solution, moles, time)) => moles_for_T_range.push((T, moles)),
                Err((T, e)) => {
                    error!("Equilibrium calculation failed at T = {} K: {:?}", T, e);
                    failed_temps.push(T);
                }
            }
        }

        moles_for_T_range.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        self.moles_for_T_range = moles_for_T_range.clone();
        self.list_of_failed_T = failed_temps;
        self.map_of_moles_for_each_substance();
        println!("time for parallel solving {:?} \n", now.elapsed());
        Ok(moles_for_T_range)
    }

    pub fn solve_for_T_range_par2(
        &mut self,
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Result<Vec<(f64, Vec<f64>)>, ReactionExtentError> {
        use rayon::prelude::*;
        use std::sync::Arc;
        let now = Instant::now();

        let subs = self.subs_data.substances.clone();

        // Step 1: Identify temperature ranges where coefficients change (single O(n) pass)
        let mut temp_ranges: Vec<(f64, f64)> = Vec::new();
        let mut T = T_start;

        while T < T_end {
            let user_subs = &mut self.subs_data;
            let vec_of_coeffs = user_subs
                .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T)
                .map_err(|_| {
                    ReactionExtentError::SubsDataError(SubsDataError::CoefficientExtractionFailed {
                        substance: "all_subs".to_string(),
                        temperature: Some(T),
                    })
                })?;

            if vec_of_coeffs.len() != 0 || self.gibbs.len() == 0 {
                let range_start = T;
                let mut T_next = T + T_step;

                // Find next temperature where coefficients change
                while T_next < T_end {
                    if user_subs
                        .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T_next)
                        .unwrap_or_default()
                        .len()
                        != 0
                    {
                        break;
                    }
                    T_next += T_step;
                }

                temp_ranges.push((range_start, T_next));
                T = T_next;
            } else {
                T += T_step;
            }
        }
        let step2 = now.elapsed().as_millis();
        // Step 2: Pre-compute Gibbs functions for each identified range
        let mut gibbs_cache: HashMap<
            usize,
            (Vec<Arc<dyn Fn(f64) -> f64 + Send + Sync>>, f64, f64),
        > = HashMap::new();

        for (range_id, (T_start_range, T_end_range)) in temp_ranges.iter().enumerate() {
            let user_subs = &mut self.subs_data;
            user_subs
                .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(*T_start_range)
                .map_err(|_| {
                    ReactionExtentError::SubsDataError(SubsDataError::CoefficientExtractionFailed {
                        substance: "all_subs".to_string(),
                        temperature: Some(*T_start_range),
                    })
                })?;

            // Create Gibbs functions for this range as Arc (thread-safe)
            let g_vec: Vec<Arc<dyn Fn(f64) -> f64 + Send + Sync>> = {
                let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();
                let mut local_g_vec: Vec<Arc<dyn Fn(f64) -> f64 + Send + Sync>> = Vec::new();
                for key in &subs {
                    let boxed_fn = map_of_functions.remove(key).unwrap();
                    let gibbs_fn: Arc<dyn Fn(f64) -> f64 + Send + Sync> =
                        Arc::new(move |T: f64| boxed_fn(T));
                    local_g_vec.push(gibbs_fn);
                }
                local_g_vec
            };

            gibbs_cache.insert(range_id, (g_vec, *T_start_range, *T_end_range));
        }
        let step3 = now.elapsed().as_millis();
        // Step 2: Generate all temperature points with their range indices
        let mut temp_range_pairs: Vec<(f64, usize)> = Vec::new();
        T = T_start;

        while T < T_end {
            // Find which range this temperature belongs to
            let mut found_range_id = None;
            for (id, (_gibbs_vec, t_start, t_end)) in &gibbs_cache {
                if T >= *t_start && T < *t_end {
                    found_range_id = Some(*id);
                    break;
                }
            }

            if let Some(range_id) = found_range_id {
                temp_range_pairs.push((T, range_id));
            }
            T += T_step;
        }
        let step4 = now.elapsed().as_millis();
        // Step 3: Parallel computation
        let elem_composition = self.elem_composition.clone();
        let stoich_matrix = self.stoich_matrix.clone();
        let n0 = self.n0.clone();
        let phases = self.phases.clone();
        let elements_vector = self.elements_vector.clone();
        let P = self.P;
        let p0 = self.p0;
        let solver_params = self.solver_params.clone();
        let species_phase = self.species_phase.clone();
        let solver = self.solver;
        let scaling_flag = self.scaling_flag;
        let species_eps = self.species_eps;
        let substate_eps = self.substate_eps;
        let gibbs_cache = Arc::new(gibbs_cache);
        println!("time for creation closure {:?}", now.elapsed());
        println!(" {}, {}, {}", step2, step3, step4);
        let results: Vec<Result<(f64, Vec<f64>), (f64, ReactionExtentError)>> = temp_range_pairs
            .par_iter()
            .map(|(T, range_id)| {
                // Retrieve the Gibbs functions for this range
                let gibbs_vec = if let Some((g_vec, _, _)) = gibbs_cache.get(range_id) {
                    g_vec.clone()
                } else {
                    return Err((
                        *T,
                        ReactionExtentError::Other("Gibbs cache not found".to_string()),
                    ));
                };

                let mut local_solver = EquilibriumLogMoles::empty();
                local_solver.subs_data.substances = subs.clone();
                local_solver.elem_composition = elem_composition.clone();
                local_solver.stoich_matrix = stoich_matrix.clone();
                local_solver.species_phase = species_phase.clone();
                local_solver.n0 = n0.clone();
                local_solver.phases = phases.clone();
                local_solver.elements_vector = elements_vector.clone();
                local_solver.P = P;
                local_solver.p0 = p0;
                local_solver.T = *T;
                local_solver.solver_params = solver_params.clone();
                local_solver.solver = solver;
                local_solver.scaling_flag = scaling_flag;
                local_solver.species_eps = species_eps;
                local_solver.substate_eps = substate_eps;

                // Convert Arc-based gibbs functions back to Rc for internal use
                let gibbs_as_rc: Vec<GibbsFn> = gibbs_vec
                    .iter()
                    .map(|arc_fn| {
                        let fn_clone = arc_fn.clone();
                        Rc::new(move |T: f64| fn_clone(T)) as GibbsFn
                    })
                    .collect();

                local_solver.gibbs = gibbs_as_rc;

                match local_solver.solve() {
                    Ok(()) => Ok((*T, local_solver.moles.clone())),
                    Err(e) => Err((*T, e)),
                }
            })
            .collect();

        // Step 4: Process results
        let mut moles_for_T_range = Vec::new();
        let mut failed_temps = Vec::new();

        for result in results {
            match result {
                Ok((T, moles)) => moles_for_T_range.push((T, moles)),
                Err((T, e)) => {
                    error!("Equilibrium calculation failed at T = {} K: {:?}", T, e);
                    failed_temps.push(T);
                }
            }
        }

        moles_for_T_range.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        self.moles_for_T_range = moles_for_T_range.clone();
        self.list_of_failed_T = failed_temps;
        self.map_of_moles_for_each_substance();
        println!("time for parallel solving {:?} \n", now.elapsed());
        Ok(moles_for_T_range)
    }

    pub fn map_of_moles_for_each_substance(&mut self) {
        let user_subs = &mut self.subs_data.substances;
        let mut vec_of_T = Vec::new();
        for (T, _moles) in &self.moles_for_T_range {
            vec_of_T.push(*T);
        }
        for (i, sub_i) in user_subs.iter().enumerate() {
            let name = sub_i.clone();
            let mut vec_of_moles: Vec<f64> = Vec::new();
            for (_T, moles) in &self.moles_for_T_range {
                vec_of_moles.push(moles[i]);
            }
            self.map_of_moles_for_each_substance
                .insert(name, vec_of_moles);
        }
        self.map_of_moles_for_each_substance
            .insert("T".to_string(), vec_of_T);
    }

    pub fn create_moles_table(&self) {
        let mut table = Table::new();

        // Get the keys (headers) and sort them for consistent order
        let mut keys: Vec<&String> = self.map_of_moles_for_each_substance.keys().collect();
        keys.sort();

        // Create header row
        let header_row: Vec<Cell> = keys.iter().map(|k| Cell::new(k)).collect();
        table.add_row(Row::new(header_row));

        // Determine the number of rows (length of the vectors)
        let num_rows = if let Some(first_vec) = self.map_of_moles_for_each_substance.values().next()
        {
            first_vec.len()
        } else {
            0
        };

        // Add data rows
        for row_idx in 0..num_rows {
            let mut row: Vec<Cell> = Vec::new();
            for key in &keys {
                if let Some(vec) = self.map_of_moles_for_each_substance.get(*key) {
                    if row_idx < vec.len() {
                        row.push(Cell::new(&format!("{:.6}", vec[row_idx])));
                    } else {
                        row.push(Cell::new(""));
                    }
                } else {
                    row.push(Cell::new(""));
                }
            }
            table.add_row(Row::new(row));
        }

        table.printstd();
    }

    /// Solve the equilibrium and return an error enum on failure
    pub fn solve(&mut self) -> Result<(), ReactionExtentError> {
        // 1. Stoichiometric matrix ν_{ik}

        let m = self.n0.len();
        assert_eq!(self.gibbs.len(), m);
        let subs_eps = self.substate_eps;
        let phase_eps = self.species_eps;
        let initial_guess = match &self.initial_guess {
            Some(g) => g.clone(),
            // safer default: no reaction progress
            None => {
                self.initial_guess = Some(self.n0.clone());
                self.n0.clone()
            }
        };
        self.check_task()?;
        //

        // 3. Structural (ξ-independent) data
        // Compute species-to-phase map and related data (fail early if mapping invalid)
        let species_phase_ = self.species_phase.clone();
        info!("Species to phase map: {:?}", species_phase_);
        info!("Stoichiometric matrix:\n{}", self.stoich_matrix);
        let delta_n = reaction_phase_stoichiometry(&self.stoich_matrix, &self.phases);
        info!("Δn matrix: {:?}", delta_n);
        let n_phases = self.phases.len();
        info!("Number of phases: {}", n_phases);

        // Shared immutable data
        let stoich = Rc::new(self.stoich_matrix.clone());
        let n0 = Rc::new(self.n0.clone());
        let species_phase = Rc::new(species_phase_.clone());
        let delta_n = Rc::new(delta_n);

        // 4. Feasibility: species disappearance
        let feasible = {
            move |xi: &[f64]| {
                let n = compute_species_moles(xi);
                n.iter().all(|&ni| ni >= -1e-12)
            }
        };

        // 5. Residuals
        let elem_composition = self.elem_composition.clone();

        let elements_vector = self.elements_vector.clone();
        let f_raw = equilibrium_logmole_residual(
            (*stoich).clone(),
            elem_composition.clone(),
            elements_vector.clone(),
            self.gibbs.clone(),
            self.phases.clone(),
            self.T,
            self.P,
            self.p0,
            species_phase_,
            subs_eps,
            phase_eps,
        )?;
        info!("Initial guess: {:?}", &initial_guess);
        info!("n0 (species initial moles): {:?}", &self.n0);
        info!("Stoichiometric matrix:\n{}", &self.stoich_matrix);

        let try_step = f_raw(&initial_guess)?;
        if try_step.iter().any(|v| v.is_nan() || v.is_infinite()) {
            error!(
                "Initial guess produces invalid residuals: {:?}",
                &try_step.clone()
            );
            println!("Initial guess produces invalid residuals: {:?}", try_step);
            return Err(ReactionExtentError::InvalidInitialResiduals(try_step));
        }
        // 6. Jacobian
        let j_raw = {
            let stoich = stoich.clone();
            let _n0 = n0.clone();
            let species_phase = species_phase.clone();
            let delta_n = delta_n.clone();

            Box::new(move |y: &[f64]| {
                equilibrium_logmole_jacobian(
                    y,
                    &stoich,
                    &elem_composition,
                    &species_phase,
                    delta_n.as_ref(),
                    n_phases,
                    subs_eps,
                    phase_eps,
                )
            }) as Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>
        };
        // Scaled closures
        let scale = equilibrium_scaling(
            &stoich,
            &self.elem_composition,
            &self.gibbs,
            &elements_vector,
            self.T,
        );
        let f: Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>> = if self.scaling_flag {
            scaled_residual(Box::new(f_raw), scale.clone())
        } else {
            Box::new(f_raw)
        };
        let j: Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>> =
            if self.scaling_flag {
                scaled_jacobian(Box::new(j_raw), scale)
            } else {
                Box::new(j_raw)
            };
        // 7. Solve
        let sol = self.solver_impl(initial_guess, f, j, Box::new(feasible))?;
        // from log moles to moles
        self.compute_species_moles(sol.clone());
        self.solution = sol;
        Ok(())
    }
    /// Automatic Fallback : If the user's chosen solver fails, it automatically tries the other two solvers
    /// No Duplicates : Creates a unique list of solvers starting with the user's choice
    ///Ordered Attempts : Tries solvers in order: [User's Choice, LM, NR, TR] (removing duplicates)
    /// Logging : Logs when fallback occurs and which solver eventually succeeds
    fn solver_impl(
        &self,
        initial_guess: Vec<f64>,
        f: Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>,
        j: Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>,
        feasible: Box<impl Fn(&[f64]) -> bool>,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        let solvers = [self.solver, Solvers::LM, Solvers::NR, Solvers::TR];
        let mut unique_solvers = Vec::new();

        // Remove duplicates while preserving order
        for solver in solvers {
            if !unique_solvers.contains(&solver) {
                unique_solvers.push(solver);
            }
        }

        for (attempt, &solver) in unique_solvers.iter().enumerate() {
            let result = match solver {
                Solvers::LM => {
                    let mut lm_solver = LMSolver {
                        f: f.as_ref(),
                        jacobian: j.as_ref(),
                        feasible: feasible.as_ref(),
                        lambda: self.solver_params.lambda,
                        tol: 1e-12,
                        max_iter: self.solver_params.max_iter,
                        alpha_min: self.solver_params.alpha_min,
                    };
                    lm_solver.solve(initial_guess.clone())
                }
                Solvers::NR => {
                    let mut nr_solver = NRSolver {
                        f: f.as_ref(),
                        jacobian: j.as_ref(),
                        feasible: feasible.as_ref(),
                        n0: self.n0.clone(),
                        reactions: self.stoich_matrix.clone(),
                        tol: self.solver_params.tol,
                        max_iter: self.solver_params.max_iter,
                        alpha_min: self.solver_params.alpha_min,
                    };
                    nr_solver.solve(initial_guess.clone())
                }
                Solvers::TR => {
                    let tr_solver = TrustRegionSolver {
                        f: f.as_ref(),
                        jacobian: j.as_ref(),
                        feasible: feasible.as_ref(),
                        tol: self.solver_params.tol,
                        max_iter: self.solver_params.max_iter,
                        delta_init: self.solver_params.delta_init,
                        delta_max: self.solver_params.delta_max,
                        eta: self.solver_params.eta,
                    };
                    tr_solver.solve(initial_guess.clone())
                }
            };

            match result {
                Ok(solution) => {
                    if attempt > 0 {
                        info!(
                            "Solver {:?} succeeded after {:?} failed",
                            solver, self.solver
                        );
                    }
                    return Ok(solution);
                }
                Err(e) => {
                    if attempt == unique_solvers.len() - 1 {
                        return Err(ReactionExtentError::SolveError(e));
                    }
                    info!("Solver {:?} failed, trying next solver", solver);
                }
            }
        }

        Err(ReactionExtentError::Other("All solvers failed".to_string()))
    }

    /// Converts log-mole solution to actual mole numbers
    ///
    /// Transforms the solver output (ln(n_i)) back to mole numbers (n_i = exp(ln(n_i))).
    pub fn compute_species_moles(&mut self, sol: Vec<f64>) {
        let moles = compute_species_moles(&sol);
        let map_of_moles_for_each_substance: HashMap<String, Vec<f64>> = self
            .subs_data
            .substances
            .clone()
            .iter()
            .zip(moles.clone())
            .map(|(subs, n)| (subs.clone(), vec![n]))
            .collect();
        self.moles = moles;
        self.map_of_moles_for_each_substance = map_of_moles_for_each_substance;
    }

    /// Validates problem setup for dimensional consistency
    ///
    /// Ensures that initial moles, substances, and initial guess vectors
    /// all have consistent dimensions.
    fn check_task(&self) -> Result<(), ReactionExtentError> {
        let n_subs = self.subs_data.substances.len();
        if self.n0.len() != n_subs {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial moles number of substances {} doesn't
            match number of substances {}",
                self.n0.len(),
                n_subs
            )));
        }
        let ig = self
            .initial_guess
            .clone()
            .expect("initial guess not defined")
            .len();
        if ig != n_subs {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial guess for solution {} doesn't
            match number of substances {}",
                ig, n_subs
            )));
        }
        Ok(())
    }
}
/// Computes species moles from log-mole variables
///
/// Converts log-mole solution vector back to actual mole numbers by exponentiating.
pub fn compute_species_moles(n: &[f64]) -> Vec<f64> {
    let n: Vec<f64> = n.to_vec().iter().map(|yi| yi.exp()).collect();
    n
}

/// Computes total element amounts from species composition
///
/// Multiplies elemental composition matrix by initial mole vector to get
/// total amount of each element in the system.
pub fn compute_element_totals(a: &DMatrix<f64>, n0: &Vec<f64>) -> DVector<f64> {
    let n0 = DVector::from_vec(n0.clone());
    &a.transpose() * n0
}

/// Calculates standard Gibbs free energy change for each reaction
///
/// Computes ΔG°_rxn = Σ ν_i * G°_i for each independent reaction.
pub fn reaction_standard_gibbs(stoich: &DMatrix<f64>, gibbs: &[GibbsFn], T: f64) -> Vec<f64> {
    let m = stoich.nrows();
    let r = stoich.ncols();

    let mut dg0 = vec![0.0; r];

    for k in 0..r {
        let mut sum = 0.0;
        for i in 0..m {
            sum += stoich[(i, k)] * gibbs[i](T);
        }
        dg0[k] = sum;
    }
    dg0
}

/// Computes scaling factors for equilibrium equations
///
/// Generates appropriate scaling to improve numerical conditioning
/// by normalizing reaction and element balance equations.
pub fn equilibrium_scaling(
    stoich: &DMatrix<f64>,   // m × r
    elements: &DMatrix<f64>, // m × E
    gibbs: &[GibbsFn],       // m
    element_totals: &[f64],  // E
    T: f64,
) -> Vec<f64> {
    let m = stoich.nrows();
    let r = stoich.ncols();
    let e = elements.ncols();

    let rt = 8.314462618 * T;

    // --- reaction scaling ---
    let dg0 = reaction_standard_gibbs(stoich, gibbs, T);
    let mut scale = vec![0.0; r + e];

    for k in 0..r {
        let mut nu_norm = 0.0;
        for i in 0..m {
            let nu = stoich[(i, k)];
            nu_norm += nu * nu;
        }
        nu_norm = nu_norm.sqrt();

        scale[k] = dg0[k].abs().max(rt * nu_norm).max(10.0); // dimensionless, log-scale
    }

    // --- element balance scaling ---
    for el in 0..e {
        scale[r + el] = element_totals[el].abs().max(10.0);
    }

    scale
}

/// Types of thermodynamic phases supported
#[derive(Debug, Clone, Copy)]
pub enum PhaseKind {
    /// Ideal gas phase with pressure dependence
    IdealGas,
    /// Ideal solution phase (condensed phases)
    IdealSolution,
}

/// Represents a thermodynamic phase containing multiple species
#[derive(Debug, Clone)]
pub struct Phase {
    /// Type of phase (ideal gas or ideal solution)
    pub kind: PhaseKind,
    /// Indices of species belonging to this phase
    pub species: Vec<usize>,
}

/// Function type for Gibbs free energy evaluation
///
/// Takes temperature in K and returns standard Gibbs free energy in J/mol.
pub type GibbsFn = Rc<dyn Fn(f64) -> f64>;

/// Maps species indices to their corresponding phase indices
///
/// Validates that each species belongs to exactly one phase and returns
/// a vector where element i gives the phase index for species i.
pub fn species_to_phase_map(
    phases: &[Phase],
    num_species: usize,
) -> Result<Vec<usize>, ReactionExtentError> {
    let mut map = vec![None; num_species];

    for (j, phase) in phases.iter().enumerate() {
        for &i in &phase.species {
            if map[i].is_some() {
                error!(
                    "Species {} appears in more than one phase. Duplicate species per phase instead (e.g. O2_g, O2_l).",
                    i
                );
                return Err(ReactionExtentError::DuplicateSpecies(i));
            }
            // Assign species i to phase j
            map[i] = Some(j);
        }
    }

    let mut out = Vec::with_capacity(num_species);
    for x in map.into_iter() {
        match x {
            Some(v) => out.push(v),
            None => {
                error!("Species not assigned to any phase");
                return Err(ReactionExtentError::SpeciesNotAssigned);
            }
        }
    }
    Ok(out)
}

/// Computes change in mole numbers per phase for each reaction
///
/// For each reaction k and phase j, calculates Δn[k][j] = sum of stoichiometric
/// coefficients for all species in phase j participating in reaction k.
pub fn reaction_phase_stoichiometry(
    reactions: &DMatrix<f64>, // m × r
    phases: &[Phase],
) -> Vec<Vec<f64>> {
    let r = reactions.ncols();

    let mut delta_n = vec![vec![0.0; phases.len()]; r];

    for (j, phase) in phases.iter().enumerate() {
        for &i in &phase.species {
            for k in 0..r {
                delta_n[k][j] += reactions[(i, k)];
            }
        }
    }

    delta_n
}

#[inline]
fn species_active(n: f64, eps: f64) -> bool {
    n > eps
}

#[inline]
fn phase_active(n: f64, eps: f64) -> bool {
    n > eps
}
/// Generates residual function for log-mole equilibrium equations
///
/// Creates closure that evaluates equilibrium residuals for given log-mole variables.
/// Handles both reaction equilibrium and element balance constraints.
pub fn equilibrium_logmole_residual(
    reactions: DMatrix<f64>,  // m × r
    elements: DMatrix<f64>,   // m × E
    element_totals: Vec<f64>, // E
    gibbs: Vec<GibbsFn>,      // m
    phases: Vec<Phase>,
    temperature: f64,
    pressure: f64,
    p0: f64,
    species_phase: Vec<usize>,
    subs_eps: f64,
    phase_eps: f64,
) -> Result<Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    if gibbs.len() != m || element_totals.len() != e {
        return Err(ReactionExtentError::Other(
            "Dimension mismatch in log-mole residual".to_string(),
        ));
    }

    let delta_n = reaction_phase_stoichiometry(&reactions, &phases);
    let rt = R * temperature;

    Ok(Box::new(move |y: &[f64]| {
        if y.len() != m {
            return Err(ReactionExtentError::ResidualEvaluation(
                "log-mole vector length mismatch".to_string(),
            ));
        }

        // --- species moles ---
        let mut n = vec![0.0; m];
        for i in 0..m {
            n[i] = y[i].exp();
        }

        // --- phase totals ---
        let mut n_phase = vec![0.0; phases.len()];
        for i in 0..m {
            n_phase[species_phase[i]] += n[i];
        }

        // --- phase offsets ---
        let mut phi = vec![0.0; phases.len()];
        for (j, phase) in phases.iter().enumerate() {
            phi[j] = match phase.kind {
                PhaseKind::IdealGas => (pressure / p0).ln(),
                PhaseKind::IdealSolution => 0.0,
            };
        }

        let mut f = vec![0.0; m];

        // -------------------------
        // Reaction equilibrium eqs
        // -------------------------
        for k in 0..r {
            let mut sum_ln_n = 0.0;
            let mut sum_ln_n_phase = 0.0;
            let mut sum_phi = 0.0;
            let mut dg0 = 0.0;

            for i in 0..m {
                let nu = reactions[(i, k)];
                if nu != 0.0 {
                    sum_ln_n += nu * y[i];
                    dg0 += nu * gibbs[i](temperature);
                }
            }

            for j in 0..phases.len() {
                let dnk = delta_n[k][j];
                if dnk != 0.0 {
                    sum_ln_n_phase += dnk * n_phase[j].ln();
                    sum_phi += dnk * phi[j];
                }
            }

            let ln_k = -dg0 / rt;
            f[k] = sum_ln_n - sum_ln_n_phase + sum_phi - ln_k;
        }

        // -------------------------
        // Element balance eqs
        // -------------------------
        for el in 0..elements.ncols() {
            let mut sum = 0.0;
            for i in 0..m {
                sum += elements[(i, el)] * n[i];
            }
            f[r + el] = sum - element_totals[el];
        }

        Ok(f)
    }))
}

pub fn equilibrium_logmole_residual2(
    reactions: DMatrix<f64>,  // m × r
    elements: DMatrix<f64>,   // m × E
    element_totals: Vec<f64>, // E
    gibbs: Vec<GibbsFn>,      // m
    phases: Vec<Phase>,
    temperature: f64,
    pressure: f64,
    p0: f64,
    subs_eps: f64,
    phase_eps: f64,
) -> Result<Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    let species_phase = species_to_phase_map(&phases, m)?;
    let delta_n = reaction_phase_stoichiometry(&reactions, &phases);
    let rt = R * temperature;

    Ok(Box::new(move |y: &[f64]| {
        if y.len() != m {
            return Err(ReactionExtentError::ResidualEvaluation(
                "log-mole vector length mismatch".to_string(),
            ));
        }

        // --- species moles ---
        let mut n = vec![0.0; m];
        let mut active = vec![false; m];
        for i in 0..m {
            n[i] = y[i].exp();
            active[i] = species_active(n[i], subs_eps);
        }

        // --- phase totals ---
        let mut n_phase = vec![0.0; phases.len()];
        let mut phase_active_mask = vec![false; phases.len()];
        for i in 0..m {
            if active[i] {
                n_phase[species_phase[i]] += n[i];
            }
        }
        for j in 0..phases.len() {
            phase_active_mask[j] = phase_active(n_phase[j], phase_eps);
        }

        // --- phase offsets ---
        let mut phi = vec![0.0; phases.len()];
        for (j, phase) in phases.iter().enumerate() {
            if phase_active_mask[j] {
                phi[j] = match phase.kind {
                    PhaseKind::IdealGas => (pressure / p0).ln(),
                    PhaseKind::IdealSolution => 0.0,
                };
            }
        }

        let mut f = vec![0.0; m];

        // -------------------------
        // Reaction equilibrium eqs
        // -------------------------
        for k in 0..r {
            let mut sum_ln_n = 0.0;
            let mut sum_ln_n_phase = 0.0;
            let mut sum_phi = 0.0;
            let mut dg0 = 0.0;

            for i in 0..m {
                if !active[i] {
                    continue;
                }
                let nu = reactions[(i, k)];
                if nu != 0.0 {
                    sum_ln_n += nu * y[i];
                    dg0 += nu * gibbs[i](temperature);
                }
            }

            for j in 0..phases.len() {
                if !phase_active_mask[j] {
                    continue;
                }
                let dnk = delta_n[k][j];
                if dnk != 0.0 {
                    sum_ln_n_phase += dnk * n_phase[j].ln();
                    sum_phi += dnk * phi[j];
                }
            }

            let ln_k = -dg0 / rt;
            f[k] = sum_ln_n - sum_ln_n_phase + sum_phi - ln_k;
        }

        // -------------------------
        // Element balance eqs
        // -------------------------
        for el in 0..e {
            let mut sum = 0.0;
            for i in 0..m {
                if active[i] {
                    sum += elements[(i, el)] * n[i];
                }
            }
            f[r + el] = sum - element_totals[el];
        }

        // -------------------------
        // Inactive species constraints
        // -------------------------
        for i in 0..m {
            if !active[i] {
                f[i] = y[i]; // push to -∞
            }
        }

        Ok(f)
    }))
}

/// Creates scaled residual function for improved numerical conditioning
///
/// Wraps original residual function with scaling factors to normalize
/// equations with vastly different magnitudes.
pub fn scaled_residual(
    f: Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>,
    scale: Vec<f64>,
) -> Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>> {
    Box::new(move |xi: &[f64]| {
        let mut r = f(xi)?;

        if r.len() != scale.len() {
            return Err(ReactionExtentError::Other(
                "Residual/scale dimension mismatch".to_string(),
            ));
        }

        for k in 0..r.len() {
            let s = scale[k];
            if s <= 0.0 || !s.is_finite() {
                return Err(ReactionExtentError::Other(format!(
                    "Invalid scale[{}] = {}",
                    k, s
                )));
            }
            r[k] /= s;
        }
        Ok(r)
    })
}

/// Computes Jacobian matrix for log-mole equilibrium equations
///
/// Evaluates analytical Jacobian of equilibrium residuals with respect
/// to log-mole variables for Newton-type solvers.
pub fn equilibrium_logmole_jacobian(
    y: &[f64],
    reactions: &DMatrix<f64>, // m × r
    elements: &DMatrix<f64>,  // m × E
    species_phase: &[usize],
    delta_n: &[Vec<f64>],
    n_phases: usize,

    subs_eps: f64,
    phase_eps: f64,
) -> Result<DMatrix<f64>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    let mut n = vec![0.0; m];
    for i in 0..m {
        n[i] = y[i].exp();
    }

    let mut n_phase = vec![0.0; n_phases];
    for i in 0..m {
        n_phase[species_phase[i]] += n[i];
    }

    let mut jmat = DMatrix::<f64>::zeros(m, m);

    // --- equilibrium equations ---
    for k in 0..r {
        for i in 0..m {
            let nu = reactions[(i, k)];
            if nu == 0.0 {
                continue;
            }
            let phase = species_phase[i];
            let term = delta_n[k][phase] * n[i] / n_phase[phase];
            jmat[(k, i)] = nu - term;
        }
    }

    // --- element balances ---
    for el in 0..e {
        for i in 0..m {
            jmat[(r + el, i)] = elements[(i, el)] * n[i];
        }
    }

    Ok(jmat)
}

pub fn equilibrium_logmole_jacobian2(
    y: &[f64],
    reactions: &DMatrix<f64>, // m × r
    elements: &DMatrix<f64>,  // m × E
    species_phase: &[usize],
    delta_n: &[Vec<f64>],
    n_phases: usize,
    subs_eps: f64,
    phase_eps: f64,
) -> Result<DMatrix<f64>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    let mut n = vec![0.0; m];
    let mut active = vec![false; m];
    for i in 0..m {
        n[i] = y[i].exp();
        active[i] = species_active(n[i], subs_eps);
    }

    let mut n_phase = vec![0.0; n_phases];
    let mut phase_active_mask = vec![false; n_phases];
    for i in 0..m {
        if active[i] {
            n_phase[species_phase[i]] += n[i];
        }
    }
    for j in 0..n_phases {
        phase_active_mask[j] = phase_active(n_phase[j], phase_eps);
    }

    let mut jmat = DMatrix::<f64>::zeros(m, m);

    // --- reaction equilibrium rows ---
    for k in 0..r {
        for i in 0..m {
            if !active[i] {
                continue;
            }
            let nu = reactions[(i, k)];
            if nu == 0.0 {
                continue;
            }

            let phase = species_phase[i];
            if phase_active_mask[phase] {
                let term = delta_n[k][phase] * n[i] / n_phase[phase];
                jmat[(k, i)] = nu - term;
            } else {
                jmat[(k, i)] = nu;
            }
        }
    }

    // --- element balances ---
    for el in 0..e {
        for i in 0..m {
            if active[i] {
                jmat[(r + el, i)] = elements[(i, el)] * n[i];
            }
        }
    }

    // --- inactive species diagonal ---
    for i in 0..m {
        if !active[i] {
            jmat[(i, i)] = 1.0;
        }
    }

    Ok(jmat)
}

/// Creates scaled Jacobian function for improved conditioning
///
/// Applies same scaling factors used for residuals to the Jacobian matrix.
pub fn scaled_jacobian(
    j: Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>,
    scale: Vec<f64>,
) -> Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>> {
    Box::new(move |xi: &[f64]| {
        let mut J = j(xi)?;

        if J.nrows() != scale.len() {
            return Err(ReactionExtentError::Other(
                "Jacobian/scale dimension mismatch".to_string(),
            ));
        }

        for k in 0..J.nrows() {
            let s = scale[k];
            if s <= 0.0 || !s.is_finite() {
                return Err(ReactionExtentError::Other(format!(
                    "Invalid scale[{}] = {}",
                    k, s
                )));
            }
            for jcol in 0..J.ncols() {
                J[(k, jcol)] /= s;
            }
        }
        Ok(J)
    })
}
