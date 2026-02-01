use super::ClassicalThermodynamics::{Thermodynamics, ThermodynamicsError};
use crate::Thermodynamics::User_PhaseOrSolution::ThermodynamicsCalculatorTrait;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;

impl Thermodynamics {
    /// Function to create closure for conservation of  the elements  condition
    ///  SUM_j (a_ij * x_j) = b_i
    pub fn composition_equations(&mut self) -> Result<(), ThermodynamicsError> {
        if let Some(A) = self.elem_composition_matrix.clone() {
            let vec_of_substances = self.vec_of_subs.clone();

            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
                return Err(ThermodynamicsError::MatrixDimensionMismatch(format!(
                    "Number of reactions from matrix ({}) doesn't match substance vector length ({})",
                    A.ncols(),
                    vec_of_substances.len()
                )));
            }

            let mut vector_of_conditions: Vec<Box<dyn Fn(DVector<f64>) -> f64>> = Vec::new();
            for i in 0..n_elements {
                let initial_vector_of_elements = self.initial_vector_of_elements.clone();
                let row_of_element_i = A.row(i).iter().map(|&v| v).collect::<Vec<f64>>();
                let row_of_element_i: DVector<_> = DVector::from_vec(row_of_element_i);
                let condition_closure = move |vec_of_concentrations: DVector<f64>| {
                    let b_i = row_of_element_i.dot(&vec_of_concentrations);

                    let condition_i = b_i - initial_vector_of_elements[i].clone();
                    condition_i
                };
                vector_of_conditions.push(Box::new(condition_closure));
            }
            self.solver.elements_conditions = vector_of_conditions;
        }
        Ok(())
    }
    /// This function generates symbolic equations for conservation of elements in a chemical reaction. It takes a vector of symbolic
    ///  concentrations as input and returns a vector of symbolic equations representing the conservation of each element.
    /// SUM_j (a_ij * x_j) = b_i
    /// The equations are constructed by multiplying the element composition matrix (transposed) with the symbolic concentrations and
    ///  equating the result to the initial amount of each element. The resulting equations are stored in self.solver.elements_conditions_sym.
    pub fn composition_equation_sym(&mut self) -> Result<(), ThermodynamicsError> {
        let vec_of_concentrations_sym: Vec<Expr> = self.solver.n.clone();
        if let Some(A) = self.elem_composition_matrix.clone() {
            let vec_of_substances = self.vec_of_subs.clone();

            let A = A.transpose(); // now rows are elements and columns are substances
            let n_elements = A.nrows();
            if A.ncols() != vec_of_substances.len() {
                println!(
                    "The number of columns {} is not equal to the number of elements {}",
                    A.ncols(),
                    vec_of_substances.len()
                );
                return Err(ThermodynamicsError::MatrixDimensionMismatch(format!(
                    "Number of reactions from matrix ({}) doesn't match substance vector length ({})",
                    A.ncols(),
                    vec_of_substances.len()
                )));
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

                let condition_i =
                    b_i_sum.simplify() - initial_vector_of_element_sym[i].clone().simplify();
                vector_of_conditions.push(condition_i.simplify());
            }
            self.solver.elements_conditions_sym = vector_of_conditions;
        }
        Ok(())
    }

    /// Creates the nonlinear system of equations for the classical thermodynamics solver.
    ///
    /// This function creates the system of equations needed to solve the problem of finding the
    /// chemical equilibrium composition of a system. The system of equations is built using the
    /// Lagrange multipliers method, which is used to minimize the total free energy of the system
    /// subject to the constraints of element conservation.
    ///
    /// The function takes the transpose of the element composition matrix, and then uses the
    /// `calculate_Lagrange_equations_sym` function of the `SubsData` struct to calculate the
    /// Lagrange equations for the system. Finally, the function stores the equations in the
    /// `solver` struct.
    ///
    /// # Errors
    ///
    /// The function returns an `Err` if the `elem_composition_matrix` is not set in the `SubsData`
    /// struct, or if the `calculate_Lagrange_equations_sym` function fails.
    pub fn create_nonlinear_system_sym(&mut self) -> Result<(), ThermodynamicsError> {
        let A = self.elem_composition_matrix.clone().unwrap();
        // println!("A: {}, ncols: {}, nrows: {}", A, A.ncols(), A.nrows()) ;
        let Tm = self.Tm;

        let eq = self
            .subdata
            .calculate_Lagrange_equations_sym(A, Tm)
            .or_else(|_| {
                Err(ThermodynamicsError::CalculationError(format!(
                    "Failed to calculate Lagrange equations"
                )))
            })?;
        self.solver.eq_mu = eq;
        Ok(())
    }

    pub fn create_nonlinear_system_fun(&mut self) -> Result<(), ThermodynamicsError> {
        let T = self.T;
        let Tm = self.Tm;
        self.calculate_Gibbs_fun(T);
        let A = self.elem_composition_matrix.clone().unwrap();

        self.subdata.calculate_Gibbs_fun(T, self.P);

        let Lagrange_equations: Box<
            dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static,
        > = self
            .subdata
            .calculate_Lagrange_equations_fun2(A, T, self.P, Tm)
            .unwrap();
        let eq = move |T: f64,
                       vec_of_concentrations: Option<Vec<f64>>,
                       Np: Option<f64>,
                       Lambda: Vec<f64>| {
            Lagrange_equations(T, vec_of_concentrations, Np, Lambda)
        };
        let f = Box::new(eq)
            as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>;
        self.solver.eq_mu_fun = f;
        Ok(())
    }

    /// creates symbolic equations representing the idea that sum  of mole numbers in this phase is equal to the total number of moles in this phase
    pub fn create_sum_of_mole_numbers_sym(&mut self) -> Result<Vec<Expr>, std::io::Error> {
        let sym_variables = self.symbolic_vars.clone();
        let mut vec_of_sum_of_mole_numbers: Vec<Expr> = Vec::new();
        for (_phase, phase_data) in sym_variables.iter() {
            let phase_data = phase_data.clone();
            if let (Some(Np), Some(vec_of_mole_numbers)) = (phase_data.0, phase_data.1) {
                let sum_of_mole_numbers = vec_of_mole_numbers
                    .iter()
                    .fold(Expr::Const(0.0), |acc, x| acc + x.clone());
                let eq = sum_of_mole_numbers.simplify() - Np;
                vec_of_sum_of_mole_numbers.push(eq);
            }
        }
        //  println!("vec_of_sum_of_mole_numbers {:?}", vec_of_sum_of_mole_numbers);
        self.solver.eq_sum_mole_numbers = vec_of_sum_of_mole_numbers.clone();
        Ok(vec_of_sum_of_mole_numbers)
    }
    /// construct the full system of equations for the nonlinear solver.
    /// It takes the Lagrange equations, the element conservation equations, and the sum of mole numbers equations and creates a single system of equations.
    /// Check if the number of equations is equal to the number of variables
    pub fn form_full_system_sym(&mut self) -> Result<(), std::io::Error> {
        let mut eq = self.solver.eq_mu.clone();
        //adding equations SUM j (a_ij * x_j) - b_i = 0
        eq.extend(self.solver.elements_conditions_sym.clone());
        eq.extend(self.solver.eq_sum_mole_numbers.clone());
        //  println!("equation system {:?}", eq);
        let n_eq = eq.len(); // number of equations
        let mut x = self.solver.n.clone();
        x.extend(self.solver.Lambda.clone());
        x.extend(self.solver.Np.clone());
        let n_vars = x.len(); // number of variables

        if n_eq != n_vars {
            println!(
                "The number of equations {}:{:?} the number of variables {}:{:?}",
                n_eq, eq, n_vars, x
            );
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "The number of equations is not equal to the number of variables",
            ));
        }
        let full_system_sym: Vec<Expr> = eq.iter().map(|x| x.clone().simplify()).collect();
        self.solver.full_system_sym = full_system_sym;
        self.solver.all_unknowns = x;
        Ok(())
    }

    /// Creates bounds for all variables in the solver.
    /// For mole numbers (n), the lower bound is 0.0.
    /// Returns HashMap {variable_name: (lower_bound, upper_bound)}
    pub fn create_variable_bounds(&mut self) {
        let mut bounds: HashMap<String, (f64, f64)> = HashMap::new();
        // rows are substances and columns are elements
        let A = self.elem_composition_matrix.clone().unwrap();
        // vector of elements in the initial composition
        let b = self.initial_vector_of_elements.clone();
        let mut initial_guess: Vec<f64> = Vec::new();
        // Bounds for mole numbers (n) - must be non-negative
        for (i, n_var) in self.solver.n.iter().enumerate() {
            let estimation = Self::upper_estimation_of_subs_mole_number(i, &A, &b);
            initial_guess.push(estimation / 2.0);
            bounds.insert(n_var.to_string(), (0.0, estimation));
        }

        // Bounds for Lagrange multipliers (Lambda) - unbounded
        for lambda_var in &self.solver.Lambda {
            initial_guess.push(10.0);
            bounds.insert(lambda_var.to_string(), (f64::NEG_INFINITY, f64::INFINITY));
        }

        // Bounds for total mole numbers per phase (Np) - must be non-negative
        for np_var in &self.solver.Np {
            initial_guess.push(1.0);
            bounds.insert(np_var.to_string(), (0.0, f64::INFINITY));
        }

        self.solver.bounds = Some(bounds.clone());
        self.solver.initial_guess = Some(initial_guess.clone());
    }

    fn upper_estimation_of_subs_mole_number(i: usize, A: &DMatrix<f64>, b: &Vec<f64>) -> f64 {
        let elem_composition_of_subs_i = A.row(i).transpose();
        assert_eq!(elem_composition_of_subs_i.len(), b.len());
        let mut estimation_vec = Vec::new();
        for (j, sum_of_elements_j) in b.iter().enumerate() {
            //quantity element j in the substance
            let E_ij = elem_composition_of_subs_i[j];
            if E_ij > 0.0 {
                estimation_vec.push(sum_of_elements_j / E_ij);
            }
        }
        // 1.2 is safety parameter
        let estimation = 1.2 * estimation_vec.iter().min_by(|a, b| a.total_cmp(b)).unwrap();
        estimation.clone()
    }
}
