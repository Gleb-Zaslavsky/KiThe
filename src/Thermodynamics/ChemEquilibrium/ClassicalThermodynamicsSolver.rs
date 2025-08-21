use RustedSciThe::numerical::Nonlinear_systems::NR::{Method, NR};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DVector;
use std::collections::HashMap;
/// structure for solving chemical equilibrium problem
/// its fields include:
/// 1) initial elements composition
/// 2) Lagrange equations
/// 3) equations representing the idea that sum  of mole numbers in this phase is equal to the total number of moles in this phase  
/// 4) unknowns: mole numbers and Lagrange multipliers
/// 5) instance of the solver of the system of nonlinear equations
pub struct Solver {
    /// Initial elements composition
    pub b0: Vec<f64>,
    /// closure for elements conditions: SUM_j(a_ij*n_j) = b_i =0  , i = 1..m
    pub elements_conditions: Vec<Box<dyn Fn(DVector<f64>) -> f64>>,
    /// symbolic elements conditions: SUM_j(a_ij*n_j) = b_i =0  , i = 1..m
    pub elements_conditions_sym: Vec<Expr>,
    // Lagrange multipliers
    pub Lambda: Vec<Expr>,
    /// mole numbers
    pub n: Vec<Expr>,
    ///  mole numbers of each phase
    pub Np: Vec<Expr>,
    /// equations
    pub eq_mu: Vec<Expr>, // T, n,
    ///
    pub eq_mu_fun: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>,
    ///  symbolic equations representing the idea that sum  of mole numbers in this phase is equal to the total number of moles in this phase
    pub eq_sum_mole_numbers: Vec<Expr>,
    /// full system of equations
    pub full_system_sym: Vec<Expr>,
    /// Lambda, n, Np,
    pub all_unknowns: Vec<Expr>,

    pub solution: Vec<f64>,

    pub map_of_solutions: HashMap<String, f64>,
}
impl Solver {
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
            solution: vec![],
            map_of_solutions: HashMap::new(),
        }
    } // new
    pub fn solve(
        &mut self,
        initial_guess: Option<Vec<f64>>,
        tolerance: f64,
        max_iterations: usize,
        damping_factor: Option<f64>,
        Bounds: Option<HashMap<String, (f64, f64)>>,
        method: Option<Method>,
    ) {
        let initial_guess = initial_guess.unwrap_or(vec![0.5; self.all_unknowns.len()]);

        let unknowns = self.all_unknowns.clone();
        assert_eq!(
            self.all_unknowns.len(),
            initial_guess.len(),
            "the number of unknowns should be equal to the number of initial guesses"
        );
        let unknowns: Vec<String> = unknowns.iter().map(|x| x.to_string()).collect();
        let mut solver = NR::new();
        solver.set_equation_system(
            self.full_system_sym.clone(),
            Some(unknowns),
            initial_guess,
            tolerance,
            max_iterations,
        );
        solver.set_solver_params(
            Some("info".to_string()),
            None,
            damping_factor,
            Bounds,
            method,
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
}

#[cfg(test)]
mod tests {
    use super::Solver;
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
        for eq in &solver.full_system_sym {
            println!("eq: {}", eq);
        }
        solver.solve(
            Some(vec![0.5, 0.5, 0.5, 1.0, 2.0, 2.0]),
            5.0 * 1e-3,
            200,
            Some(1.0),
            None,
            None,
        );
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
        for eq in &solver.full_system_sym {
            println!("eq: {}", eq);
        }
        solver.solve(
            Some(vec![0.5, 0.5, 0.5, 1.0, 2.0, 2.0]),
            5.0 * 1e-3,
            1000,
            Some(0.009),
            None,
            None,
        );
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
    }
}
