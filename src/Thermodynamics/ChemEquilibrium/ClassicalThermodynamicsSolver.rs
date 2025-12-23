//! Classical Thermodynamics Solver Module
//!
//! This module provides solvers for chemical equilibrium problems using classical thermodynamics.
//! It implements Newton-Raphson and Levenberg-Marquardt methods to solve systems of nonlinear
//! equations arising from Gibbs free energy minimization with element conservation constraints.
//!
//! The main components include:
//! - Parameter structures for different solver types
//! - Unified solver interface with enum-based type selection
//! - Sanity checking for numerical stability
//! - Support for both temperature-dependent and temperature-independent solving

use RustedSciThe::numerical::Nonlinear_systems::NR::{Method, NR};
use RustedSciThe::numerical::optimization::sym_wrapper::LM;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use RustedSciThe::symbolic::symbolic_functions::Jacobian;
use nalgebra::{DMatrix, DVector, SVD};
use std::collections::HashMap;
use std::hash::Hash;

/// Parameters for Newton-Raphson solver
#[derive(Debug, Clone)]
pub struct NRParams {
    /// Initial guess for the solution vector
    pub initial_guess: Option<Vec<f64>>,
    /// Convergence tolerance for the solver
    pub tolerance: f64,
    /// Maximum number of iterations allowed
    pub max_iterations: usize,
    /// Damping factor for step size control
    pub damping_factor: Option<f64>,
    /// Logging level for solver output
    pub loglevel: Option<String>,
    /// Newton-Raphson method variant to use
    pub method: Option<Method>,
}

impl Default for NRParams {
    fn default() -> Self {
        Self {
            initial_guess: None,
            tolerance: 1e-6,
            max_iterations: 100,
            damping_factor: None,
            loglevel: None,
            method: None,
        }
    }
}

/// Parameters for Levenberg-Marquardt solver
#[derive(Debug, Clone)]
pub struct LMParams {
    /// Initial guess for the solution vector
    pub initial_guess: Option<Vec<f64>>,
    /// Convergence tolerance for the solver
    pub tolerance: Option<f64>,
    /// Maximum number of iterations allowed
    pub max_iterations: Option<usize>,
    /// Logging level for solver output
    pub loglevel: Option<String>,
}

impl Default for LMParams {
    fn default() -> Self {
        Self {
            initial_guess: None,
            tolerance: Some(1e-6),
            max_iterations: Some(1000),
            loglevel: Some("none".to_string()),
        }
    }
}

/// Enum for solver type selection with associated parameters
#[derive(Debug, Clone)]
pub enum SolverType {
    /// Newton-Raphson solver with its parameters
    NewtonRaphson(NRParams),
    /// Levenberg-Marquardt solver with its parameters
    LevenbergMarquardt(LMParams),
}

impl Default for SolverType {
    fn default() -> Self {
        SolverType::NewtonRaphson(NRParams::default())
    }
}
/// Enum holding actual solver instances
pub enum SolverInstance {
    /// Newton-Raphson solver instance
    NR(NR),
    /// Levenberg-Marquardt solver instance
    LM(LM),
}

impl Default for SolverInstance {
    fn default() -> Self {
        SolverInstance::NR(NR::new())
    }
}
/// Main solver structure for chemical equilibrium problems
///
/// Solves systems of nonlinear equations arising from Gibbs free energy minimization
/// with element conservation constraints using Lagrange multipliers method.
/// its fields include:
/// 1) initial elements composition
/// 2) Lagrange equations
/// 3) equations representing the idea that sum  of mole numbers in this phase is equal to the total number of moles in this phase  
/// 4) unknowns: mole numbers and Lagrange multipliers
/// 5) instance of the solver of the system of nonlinear equations
pub struct Solver {
    /// Initial elements composition vector (b₀)
    pub b0: Vec<f64>,
    /// Closure functions for element conservation conditions: ∑ⱼ(aᵢⱼ*nⱼ) = bᵢ = 0
    pub elements_conditions: Vec<Box<dyn Fn(DVector<f64>) -> f64>>,
    /// Symbolic element conservation conditions: ∑ⱼ(aᵢⱼ*nⱼ) = bᵢ = 0
    pub elements_conditions_sym: Vec<Expr>,
    /// Lagrange multipliers (λ) for element conservation constraints
    pub Lambda: Vec<Expr>,
    /// Mole numbers of individual substances (n)
    pub n: Vec<Expr>,
    /// Total mole numbers for each phase (Nₚ)
    pub Np: Vec<Expr>,
    /// Chemical potential equilibrium equations (μ equations)
    pub eq_mu: Vec<Expr>,
    /// Function version of chemical potential equations
    pub eq_mu_fun: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    /// Equations ensuring sum of mole numbers equals total phase moles
    pub eq_sum_mole_numbers: Vec<Expr>,
    /// Complete system of symbolic equations
    pub full_system_sym: Vec<Expr>,
    /// All unknown variables (λ, n, Nₚ)
    pub all_unknowns: Vec<Expr>,
    /// Initial guess for solution vector
    pub initial_guess: Option<Vec<f64>>,
    /// Final solution vector
    pub solution: Vec<f64>,
    /// Solution mapped by variable names
    pub map_of_solutions: HashMap<String, f64>,
    /// Variable bounds for constrained optimization
    pub bounds: Option<HashMap<String, (f64, f64)>>,
    /// Selected solver type and parameters
    pub solver_type: SolverType,
    /// Active solver instance
    pub solver_instance: SolverInstance,
    /// Logging level for solver output
    pub loglevel: Option<String>,
}
impl Solver {
    /// Creates a new solver instance with default values
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
            initial_guess: None,
            solution: vec![],
            map_of_solutions: HashMap::new(),
            bounds: None,
            solver_type: SolverType::default(),
            solver_instance: SolverInstance::default(),
            loglevel: Some("none".to_string()),
        }
    } // new

    /// Sets the solver type and parameters
    pub fn set_solver_type(&mut self, solver_type: SolverType) {
        self.solver_type = solver_type;
    }

    /// Solves the system using the configured solver type
    pub fn solve(&mut self) {
        match &self.solver_type {
            SolverType::NewtonRaphson(params) => self.solve_NR(params.clone()),
            SolverType::LevenbergMarquardt(params) => self.solve_lm(params.clone()),
        }
    }

    /// Generates equations for the configured solver type
    pub fn generate_eqs(&mut self) {
        match &self.solver_type {
            SolverType::NewtonRaphson(params) => self.set_NR_solver(params.clone()),
            SolverType::LevenbergMarquardt(params) => self.set_lm_solver(params.clone()),
        }
    }

    /// Solves the system for a specific temperature value
    pub fn solve_for_T(&mut self, T: f64) {
        match &self.solver_type {
            SolverType::NewtonRaphson(_) => self.solve_NR_with_T(T),
            SolverType::LevenbergMarquardt(_) => self.solve_LM_with_T(T),
        }
    }

    /// Returns the solution as a HashMap of variable names to values
    pub fn get_result(&self) -> HashMap<String, f64> {
        self.map_of_solutions.clone()
    }

    ////////////////////////////////////////////NEWTON RAPHSON//////////////////////////////////////
    /// Sets up Newton-Raphson solver with temperature parameter support
    fn set_NR_solver(&mut self, params: NRParams) {
        let bounds = self.bounds.clone();
        let loglevel = if params.loglevel.is_some() {
            params.loglevel
        } else {
            self.loglevel.clone()
        };
        let initial_guess = params
            .initial_guess
            .unwrap_or(self.initial_guess.clone().unwrap());
        //   self.sanity_check(initial_guess.clone());
        let unknowns: Vec<String> = self.all_unknowns.iter().map(|x| x.to_string()).collect();
        assert_eq!(
            self.all_unknowns.len(),
            initial_guess.len(),
            "the number of unknowns should be equal to the number of initial guesses"
        );

        let mut solver = NR::new();
        solver.set_equation_system(
            self.full_system_sym.clone(),
            Some(unknowns),
            initial_guess,
            params.tolerance,
            params.max_iterations,
        );
        solver.set_solver_params(
            loglevel,
            None,
            params.damping_factor,
            bounds,
            params.method,
            None,
        );
        solver.set_eq_params(vec!["T".to_string()]);
        solver.eq_generate();
        self.solver_instance = SolverInstance::NR(solver)
    }

    /// Solves using Newton-Raphson method with specific temperature
    fn solve_NR_with_T(&mut self, T: f64) {
        if let SolverInstance::NR(ref mut nr) = self.solver_instance {
            nr.set_eq_params_values(DVector::from_vec(vec![T]));
            nr.solve();
            let solution = nr.get_result().clone().unwrap();
            self.solution = solution.data.into();
            self.map_of_solutions = self
                .all_unknowns
                .iter()
                .zip(self.solution.iter())
                .map(|(k, v)| (k.to_string(), *v))
                .collect();
        }
    }
    /// Solves using Newton-Raphson method (requires T to be substituted in equations)
    fn solve_NR(&mut self, params: NRParams) {
        let bounds = self.bounds.clone();
        let loglevel = if params.loglevel.is_some() {
            params.loglevel
        } else {
            self.loglevel.clone()
        };
        let initial_guess = params
            .initial_guess
            .unwrap_or(self.initial_guess.clone().unwrap());
        self.sanity_check(initial_guess.clone());
        let unknowns: Vec<String> = self.all_unknowns.iter().map(|x| x.to_string()).collect();
        assert_eq!(
            self.all_unknowns.len(),
            initial_guess.len(),
            "the number of unknowns should be equal to the number of initial guesses"
        );

        let mut solver = NR::new();
        solver.set_equation_system(
            self.full_system_sym.clone(),
            Some(unknowns),
            initial_guess,
            params.tolerance,
            params.max_iterations,
        );
        solver.set_solver_params(
            loglevel,
            None,
            params.damping_factor,
            bounds,
            params.method,
            None,
        );
        solver.eq_generate();
        solver.solve();
        let solution = solver.get_result().expect("Failed to get result");
        self.solution = solution.data.into();
        self.map_of_solutions = self
            .all_unknowns
            .iter()
            .zip(self.solution.iter())
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
    }
    ///////////////////////////////LEVENBERG/////////////////////////////////////////////
    /// Solves using Levenberg-Marquardt method (requires T to be substituted in equations)
    fn solve_lm(&mut self, params: LMParams) {
        let unknowns: Vec<String> = self.all_unknowns.iter().map(|x| x.to_string()).collect();
        let initial_guess = params
            .initial_guess
            .unwrap_or(vec![0.5; self.all_unknowns.len()]);
        let mut LM = LM::new();
        LM.set_equation_system(
            self.full_system_sym.clone(),
            Some(unknowns),
            None,
            initial_guess,
            params.tolerance,
            Some(1e-6),
            Some(1e-6),
            Some(true),
            params.max_iterations,
        );
        let loglevel = if params.loglevel.is_some() {
            params.loglevel
        } else {
            self.loglevel.clone()
        };
        LM.loglevel = loglevel;
        LM.eq_generate();
        LM.solve();
        self.map_of_solutions = LM.map_of_solutions.unwrap();
    }

    /// Sets up Levenberg-Marquardt solver with temperature parameter support
    fn set_lm_solver(&mut self, params: LMParams) {
        let unknowns: Vec<String> = self.all_unknowns.iter().map(|x| x.to_string()).collect();
        let initial_guess = params
            .initial_guess
            .unwrap_or(vec![0.5; self.all_unknowns.len()]);
        let mut LM = LM::new();
        LM.set_equation_system(
            self.full_system_sym.clone(),
            Some(unknowns),
            None,
            initial_guess,
            params.tolerance,
            Some(1e-6),
            Some(1e-6),
            Some(true),
            params.max_iterations,
        );
        LM.parameters = Some(vec!["T".to_string()]);
        let loglevel = if params.loglevel.is_some() {
            params.loglevel
        } else {
            self.loglevel.clone()
        };
        LM.loglevel = loglevel;
        LM.eq_generate_with_params();

        self.solver_instance = SolverInstance::LM(LM)
    }

    /// Solves using Levenberg-Marquardt method with specific temperature
    fn solve_LM_with_T(&mut self, T: f64) {
        if let SolverInstance::LM(ref mut lm) = self.solver_instance {
            lm.solve_with_params(vec![T]);
            self.map_of_solutions = lm.map_of_solutions.clone().unwrap().clone();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    /// Returns the current solution mapping
    pub fn get_solution(&self) -> HashMap<String, f64> {
        self.map_of_solutions.clone()
    }

    pub fn get_vector_of_solution(&self) -> Vec<f64> {
        let vec_of_res = self
            .all_unknowns
            .iter()
            .map(|x| self.map_of_solutions.get(&x.to_string()).unwrap().clone())
            .collect();
        vec_of_res
    }
    /// Performs numerical stability analysis of the Jacobian matrix
    ///
    /// Checks for singularity, poor conditioning, and computes SVD diagnostics
    pub fn sanity_check(&self, initial_guess: Vec<f64>) {
        use RustedSciThe::somelinalg::linear_sys_diagnostics::{
            SVD_diagnostics, is_singular, poorly_conditioned,
        };
        use nalgebra::DMatrix;
        let mut jac = Jacobian::new();
        let unknowns: Vec<String> = self.all_unknowns.iter().map(|x| x.to_string()).collect();
        jac.set_vector_of_functions(self.full_system_sym.clone());
        jac.set_vector_of_variables(self.all_unknowns.clone());
        jac.set_variables(unknowns.clone().iter().map(|x| x.as_str()).collect());
        jac.calc_jacobian();

        jac.jacobian_generate(unknowns.iter().map(|x| x.as_str()).collect());
        let jac_func = jac.function_jacobian;
        let mut DMatrix: DMatrix<f64> = DMatrix::zeros(unknowns.len(), unknowns.len());
        for i in 0..unknowns.len() {
            for j in 0..unknowns.len() {
                DMatrix[(i, j)] = jac_func[i][j](initial_guess.clone());
            }
        }
        let singularity = match is_singular(DMatrix.clone(), 1e-4) {
            true => "Matrix is singular",
            false => "Matrix is not singular",
        };
        println!("{}", singularity);
        let poorly_conditioned = match poorly_conditioned(DMatrix.clone(), 1e5) {
            true => "Matrix is poorly conditioned",
            false => "Matrix is not poorly conditioned",
        };
        println!("{}", poorly_conditioned);
        let SVD = match SVD_diagnostics(DMatrix.clone()) {
            true => "Matrix has SVD",
            false => "Matrix does not have SVD",
        };
        println!("{}", SVD);

        let svd = SVD::new(DMatrix.clone(), true, true);
        let singular_values = svd.singular_values;
        println!("singular values: {:?}", singular_values);
        let tol = 1e-12;
        let rank = singular_values.iter().filter(|&&v| v > tol).count();
        println!("rank = {}", rank);
        if let Some(vt) = svd.v_t {
            // columns of V corresponding to singular_values <= tol are the nullspace basis
            for (i, &sv) in singular_values.iter().enumerate() {
                if sv <= tol {
                    // i-th singular value is near zero => corresponding column in V is in nullspace
                    // vt.row(i) is the transpose: vt.row(i) is one nullspace vector (length n)
                    println!("nullspace row (vt.row({})): {:?}", i, vt.row(i));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use std::vec;
    #[test]
    fn test_solver_struct() {
        let symbolic = Expr::Symbols("N0, N1, N2, Np, Lambda0, Lambda1");
        let dG0 = Expr::Const(-10.0e3);
        let dG1 = Expr::Const(-10.0e3);
        let dG2 = Expr::Const(-5e3);
        let dGm = Expr::Const(8.314 * 1e4);
        let mut solver = Solver::new();
        solver.all_unknowns = symbolic.clone();
        let N0 = symbolic[0].clone();
        let N1 = symbolic[1].clone();
        let N2 = symbolic[2].clone();
        let Np = symbolic[3].clone();
        let Lambda0 = symbolic[4].clone();
        let Lambda1 = symbolic[5].clone();

        let RT = Expr::Const(8.314) * Expr::Const(273.15);
        let eq_mu = vec![
            Lambda0.clone()
                + Expr::Const(2.0) * Lambda1.clone()
                + (dG0.clone() + RT.clone() * Expr::ln(N0.clone() / Np.clone())) / dGm.clone(),
            Lambda0
                + Lambda1.clone()
                + (dG1 + RT.clone() * Expr::ln(N1.clone() / Np.clone())) / dGm.clone(),
            Expr::Const(2.0) * Lambda1
                + (dG2 + RT * Expr::ln(N2.clone() / Np.clone())) / dGm.clone(),
        ];
        let eq_sum_mole_numbers = vec![N0.clone() + N1.clone() + N2.clone() - Np.clone()];
        let composition_eq = vec![
            N0.clone() + N1.clone() - Expr::Const(0.999),
            Expr::Const(2.0) * N0.clone() + N1.clone() + Expr::Const(2.0) * N2 - Expr::Const(1.501),
        ];

        let mut full_system_sym = Vec::new();
        full_system_sym.extend(eq_mu.clone());
        full_system_sym.extend(eq_sum_mole_numbers.clone());
        full_system_sym.extend(composition_eq.clone());
        solver.full_system_sym = full_system_sym;

        let initial_guess = vec![0.5, 0.5, 0.5, 1.0, 2.0, 2.0];

        solver.initial_guess = Some(initial_guess.clone());
        //  solver.sanity_check(initial_guess.clone());
        //panic!("stop");
        let params = NRParams {
            initial_guess: None,
            tolerance: 5.0 * 1e-3,
            max_iterations: 200,
            damping_factor: Some(1.0),
            ..Default::default()
        };
        solver.set_solver_type(SolverType::NewtonRaphson(params));
        solver.solve();
        let map_of_solutions = solver.map_of_solutions;
        let N0 = map_of_solutions.get("N0").unwrap();
        let N1 = map_of_solutions.get("N1").unwrap();
        let N2 = map_of_solutions.get("N2").unwrap();
        let Np = map_of_solutions.get("Np").unwrap();
        let _Lambda0 = map_of_solutions.get("Lambda0").unwrap();
        let _Lambda1 = map_of_solutions.get("Lambda1").unwrap();
        let d1 = *N0 + *N1 - 0.999;
        let d2 = N0 + N1 + N2 - Np;
        let d3 = 2.0 * N0 + N1 + 2.0 * N2 - 1.501;
        println!("d1: {}", d1);
        println!("d2: {}", d2);
        println!("d3: {}", d3);
        assert!(d1.abs() < 1e-3 && d2.abs() < 1e-3 && d3.abs() < 1e-3);
        println!("map_of_solutions: {:?}", map_of_solutions);
        for eq in &solver.full_system_sym {
            println!("eq: {}", eq);
        }
    }

    #[test]
    fn test_solver_struct2() {
        let symbolic = Expr::Symbols("N0, N1, N2, Np, Lambda0, Lambda1");
        let dG0 = Expr::Const(-450.0e3);
        let dG1 = Expr::Const(-150.0e3);
        let dG2 = Expr::Const(-50e3);
        let dGm0 = Expr::Const(8.314 * 450e5);
        let dGm1 = Expr::Const(8.314 * 150e5);
        let dGm2 = Expr::Const(8.314 * 50e5);
        let mut solver = Solver::new();
        solver.all_unknowns = symbolic.clone();
        let N0 = symbolic[0].clone();
        let N1 = symbolic[1].clone();
        let N2 = symbolic[2].clone();
        let Np = symbolic[3].clone();
        let Lambda0 = symbolic[4].clone();
        let Lambda1 = symbolic[5].clone();

        let RT = Expr::Const(8.314) * Expr::Const(273.15);
        let eq_mu = vec![
            Lambda0.clone()
                + Expr::Const(2.0) * Lambda1.clone()
                + (dG0.clone() + RT.clone() * Expr::ln(N0.clone() / Np.clone())) / dGm0.clone(),
            Lambda0
                + Lambda1.clone()
                + (dG1 + RT.clone() * Expr::ln(N1.clone() / Np.clone())) / dGm1.clone(),
            Expr::Const(2.0) * Lambda1
                + (dG2 + RT * Expr::ln(N2.clone() / Np.clone())) / dGm2.clone(),
        ];
        let eq_sum_mole_numbers = vec![N0.clone() + N1.clone() + N2.clone() - Np.clone()];
        let composition_eq = vec![
            N0.clone() + N1.clone() - Expr::Const(0.999),
            Expr::Const(2.0) * N0.clone() + N1.clone() + Expr::Const(2.0) * N2 - Expr::Const(1.501),
        ];

        let mut full_system_sym = Vec::new();
        full_system_sym.extend(eq_mu.clone());
        full_system_sym.extend(eq_sum_mole_numbers.clone());
        full_system_sym.extend(composition_eq.clone());
        solver.full_system_sym = full_system_sym;

        let initial_guess = vec![0.5, 0.5, 0.5, 1.0, 2.0, 2.0];
        solver.initial_guess = Some(initial_guess.clone());
        let params = NRParams {
            initial_guess: None,
            tolerance: 5.0 * 1e-3,
            max_iterations: 1000,
            damping_factor: Some(0.009),
            ..Default::default()
        };
        solver.set_solver_type(SolverType::NewtonRaphson(params));
        solver.solve();
        solver.sanity_check(initial_guess.clone());
        //  panic!("stop");
        let map_of_solutions = solver.map_of_solutions;
        let N0 = map_of_solutions.get("N0").unwrap();
        let N1 = map_of_solutions.get("N1").unwrap();
        let N2 = map_of_solutions.get("N2").unwrap();
        let Np = map_of_solutions.get("Np").unwrap();
        let _Lambda0 = map_of_solutions.get("Lambda0").unwrap();
        let _Lambda1 = map_of_solutions.get("Lambda1").unwrap();
        let d1 = *N0 + *N1 - 0.999;
        let d2 = N0 + N1 + N2 - Np;
        let d3 = 2.0 * N0 + N1 + 2.0 * N2 - 1.501;
        println!("d1: {}", d1);
        println!("d2: {}", d2);
        println!("d3: {}", d3);
        assert!(d1.abs() < 1e-3);
        assert!(d2.abs() < 1e-2);
        assert!(d3.abs() < 1e-2);
        println!("map_of_solutions: {:?}", map_of_solutions);
        for eq in &solver.full_system_sym {
            println!("eq: {}", eq);
        }
    }
}
