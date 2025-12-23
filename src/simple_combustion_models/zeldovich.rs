#![allow(non_snake_case)]
use RustedSciThe::symbolic::symbolic_engine::Expr;

use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_frozen::NRBVP;
use RustedSciThe::numerical::BVP_sci::BVP_sci_symb::BVPwrap;
use std::collections::HashMap;

use nalgebra::{DMatrix, DVector};
/*    [J/(s*m*K)]*[K/m2] = [g/s*m2]*[J/g*K]*[K/m] - [J/g]*[g/m3]^n*[ (*[g/m3]^1-n)/s ]
        [J/s*m3] =           [J/s*m3] -                 [J/s*m3]
        Lambda*d^2T/dx^2 =  m*C_p*dT/dx -  Q*ro^r*a^n*k*exp(-E/T*R)
        Lambda/l^2*d2T/d(x/l)^2 =  m*C_p/l*dT/d(x/l) -   Q*ro^r*a^n*k*k*exp(-E/T*R)
         x1 = x/l
        Lambda/l^2*d^2T/dx1^2 =  m*C_p/l*dT/dx1 -   Q*ro^r*a^n*k*exp(-E/T*R)
        d^2T/dx1^2 =  (m*C_p*l/Lambda)*dT/dx1 -  l^2/Lambda* Q*ro^r*a^n*k*exp(-E/T*R)
        Pe = m*C_p*l/Lambda
         d^2T/dx1^2 =  Pe*dT/dx1 -  l^2/Lambda*Q*ro^r*a^n*k*exp(-E/T*R)
        dT = T_f - T0
        T1 = T/dT
         d^2T1/dx1^2 =  Pe*dT1/dx1 -  l^2/Lambda* Q*ro^r*a^n*k*exp(-E/dT*T1*R)/dT
        [J/g] = [K]*[J/g*K]
        Q*a         = (T_f - T)*C_p =>    a = (T_f - T)/(T_f - T0) = (T_f - T)/dT = Tfm - T1
        ro = P/(R*T) = P/(R*dT*T1)
        ro_s = P/(R*dT)
         ro =  ro_s/T1
        E1 = E/(dT*R)
        Q/dT = C_p
        l^2/Lambda* Q*(ro*a)^n*k*exp(-E/dT*T1*R)/dT => l^2*C_p*ro_s^r*((Tfm - T1)/T1 )^n*k*exp(-E1/T1)/Lambda
         d^2T1/dx1^2 =  Pe*dT1/dx1 -  Da*((Tfm - T1)/T1 )^n*exp(-E1/T1)
         qs =  Pe*dT1/dx1
        qs' = qs -  Da*((Tfm - T1)/T1 )^n*exp(-E1/T1)


*/
#[derive(Debug, Clone)]
pub struct Par {
    pub Q_g: f64,
    pub C_p: f64,
    pub Lambda_eff: f64,
    pub M: f64,
    pub B_g: f64,
    pub E: f64,
    pub n: f64,
    pub r: f64,
}
impl Default for Par {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone)]
pub struct SolverParams {
    pub n_steps: usize,
    pub method: String,
    pub strategy: String,
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl SolverParams {
    pub fn new() -> Self {
        Self {
            n_steps: 100,
            method: "Dense".to_string(),
            strategy: "Naive".to_string(),
            tolerance: 1e-5,
            max_iterations: 100,
        }
    }
}
impl Default for SolverParams {
    fn default() -> Self {
        Self::new()
    }
}
impl Par {
    pub fn new() -> Self {
        Self {
            Q_g: 0.0,
            C_p: 0.0,
            Lambda_eff: 0.0,
            M: 0.0,
            B_g: 0.0,
            E: 0.0,
            n: 1.0,
            r: 1.0,
        }
    }

    pub fn set_params(
        Q_g: f64,
        C_p: f64,
        Lambda_eff: f64,
        M: f64,
        B_g: f64,
        E: f64,
        n: f64,
        r: f64,
    ) -> Self {
        Self {
            Q_g,
            C_p,
            Lambda_eff,
            M,
            B_g,
            E,
            n,
            r,
        }
    }
}

impl Default for Zeldovich {
    fn default() -> Self {
        Self::new()
    }
}

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct Zeldovich {
    params: Par,
    solver_params: SolverParams,
    P: f64,
    m: f64,
    L: f64,
    T0: f64,
    T0s: f64,
    Tfs: f64,
    vec_of_eq: Vec<Expr>,
    vec_of_vars: Vec<String>,
    dT: f64,
    Pe: f64,
    x_mesh: DVector<f64>,
    T_sol: DVector<f64>,
    qs_sol: DVector<f64>,
}

impl Zeldovich {
    pub fn new() -> Self {
        Self {
            params: Par::new(),
            solver_params: SolverParams::new(),
            P: 0.0,
            m: 0.0,
            L: 0.0,
            T0: 0.0,
            T0s: 0.0,
            Tfs: 0.0,
            dT: 0.0,
            Pe: 0.0,
            vec_of_eq: Vec::new(),
            vec_of_vars: Vec::new(),
            x_mesh: DVector::zeros(0),
            T_sol: DVector::zeros(0),
            qs_sol: DVector::zeros(0),
        }
    }

    /// Sets the parameters of the Zeldovich model.
    ///
    /// # Arguments
    /// * `Q_g` - Heat of gasification [J/g]
    /// * `C_p` - Specific heat capacity of the gas [J/(g·K)]
    /// * `Lambda_eff` - Effective thermal conductivity of the gas [W/(m·K)]
    /// * `M` - Molecular weight of the gas [g/mol]
    /// * `B_g` - Thermal expansion coefficient of the gas [1/K]
    /// * `E` - Activation energy of the gasification reaction [J/mol]
    pub fn set_params(
        &mut self,
        Q_g: f64,
        C_p: f64,
        Lambda_eff: f64,
        M: f64,
        B_g: f64,
        E: f64,
        n: f64,
        r: f64,
    ) {
        self.params = Par::set_params(Q_g, C_p, Lambda_eff, M, B_g, E, n, r);
    }

    pub fn set_solver_params(&mut self, solver_params: SolverParams) {
        self.solver_params = solver_params;
    }
    pub fn set_problem(&mut self, P: f64, m: f64, L: f64, T0: f64) {
        self.P = P;
        self.m = m;
        self.L = L;
        self.T0 = T0;
    }
    pub fn calc_constants(&mut self) {
        let Par {
            Q_g,
            C_p,
            Lambda_eff,
            M,
            B_g,
            E,
            n,
            r,
        } = self.params;
        let P = self.P;
        let m = self.m;
        let L = self.L;
        let T0 = self.T0;

        let Pe = m * C_p * L / Lambda_eff;
        self.Pe = Pe;

        let T_f = T0 + Q_g / C_p;
        let dT = T_f - T0;
        self.dT = dT;
        let T0s = T0 / dT;
        self.T0s = T0s;
        let Tfs = T_f / dT;
        self.Tfs = Tfs;
        let R_g = 8.314;
        let Em = E / (R_g * dT);
        let ro_s = P / (R_g * dT);

        let Da = ro_s.powf(r) * L.powf(2.0) * C_p * B_g / Lambda_eff;
        println!("Da = {} \n Pe = {}, Tf = {}", Da, Pe, T_f);
        let T = Expr::Var("T".to_string());

        let n_expr = Expr::Const(n);
        let r_expr = Expr::Const(r);

        let Source = Expr::Pow(Box::new(Expr::Const(Tfs) - T.clone()), Box::new(n_expr))
            * T.clone().pow(-r_expr)
            * Expr::Const(Da)
            * Expr::Exp(Box::new(-Expr::Const(Em) / T));
        println!("RHS = {}", Source.clone());
        let qs = Expr::Var("qs".to_string());

        let RHS_T = qs.clone() / Expr::Const(Pe);
        let RHS_q = qs - Source;
        let vec_of_eq = vec![RHS_T, RHS_q];
        self.vec_of_eq = vec_of_eq;
    }
    /*
       pub fn solve(&mut self) {
           let vec_of_eq = self.vec_of_eq.clone();
           let vec_of_vars = vec!["T".to_string(), "qs".to_string()];
           self.vec_of_vars = vec_of_vars.clone();
           let T0s = self.T0s;
           let Tfs = self.Tfs;


           let arg = "x".to_string();

           let mut boundary_conditions = HashMap::new();
          // boundary_conditions.insert("T".to_string(), vec![(0, T0s), (1, Tfs)]);
            boundary_conditions.insert("T".to_string(), vec![(0, T0s)]);
            boundary_conditions.insert("qs".to_string(), vec![ (1, 1e-7)]);
           let mut bvp = BVPShooting::new( vec_of_eq, vec_of_vars, arg, boundary_conditions, (0.0, 1.0));

           let params = HashMap::from([
               ("rtol".to_string(), SolverParam::Float(1e-6)),
               ("atol".to_string(), SolverParam::Float(1e-8)),
           ]);

           bvp.solve_with_certain_ivp(
               1.0,
               1e-6,
               100,
               0.01,
               SolverType::BDF,
               params,
           );
           let sol = bvp.get_solution();

           let boundpoint = sol.bound_values;
           let x_mesh = bvp.get_x();
           let y_mesh = bvp.get_y();
           let T_sol: DVector<f64> = y_mesh.row(0).transpose().into_owned();
           let qs_sol: DVector<f64> = y_mesh.row(1).transpose().into_owned();
           self.x_mesh = x_mesh;
           self.T_sol = T_sol;
           self.qs_sol = qs_sol;



       }
    */
    fn solve_bvp(&mut self) {
        let eq_system = self.vec_of_eq.clone();
        let values = vec!["T".to_string(), "qs".to_string()];
        self.vec_of_vars = values.clone();
        let T0s = self.T0s;
        let arg = "x".to_string();

        let mut BorderConditions = HashMap::new();
        BorderConditions.insert("T".to_string(), vec![(0, T0s)]);
        BorderConditions.insert(
            "qs".to_string(),
            vec![(1, 1e-10 * self.Pe * (self.L / self.dT))],
        );

        let t0 = 0.0;
        let t_end = 1.0;
        let SolverParams {
            n_steps,
            method,
            strategy,
            tolerance,
            max_iterations,
        } = &self.solver_params;
        let ones = vec![1.0; values.len() * n_steps];
        let initial_guess: DMatrix<f64> =
            DMatrix::from_column_slice(values.len(), *n_steps, DVector::from_vec(ones).as_slice());
        let mut nr = NRBVP::new(
            eq_system,
            initial_guess,
            values,
            arg,
            BorderConditions,
            t0,
            t_end,
            *n_steps,
            strategy.clone(),
            None,
            None,
            method.clone(),
            *tolerance,
            *max_iterations,
        );

        nr.eq_generate();

        nr.dont_save_log(true);

        let result = nr.solve();
        match result {
            Some(sol) => {
                println!("{:?}", sol);
                let dT = self.dT;
                let x_mesh = self.L * nr.x_mesh.clone();
                let y_mesh = nr.get_result().unwrap();
                let T_sol: DVector<f64> = dT * y_mesh.column(1).into_owned();
                let qs_sol: DVector<f64> =
                    (self.dT / self.L) * y_mesh.column(0).into_owned() / self.Pe;
                self.x_mesh = x_mesh;
                self.T_sol = T_sol;
                self.qs_sol = qs_sol;
                println!("BVP solved successfully");
                //    println!("x_mesh = {:?}", self.x_mesh);
                //   println!("T_sol = {:?}", self.T_sol);
                //   println!("qs_sol = {:?}", self.qs_sol);
            }
            None => {
                println!("Error during BVP solving");
            }
        }
    }

    pub fn solve_bvp_sci(&mut self) {
        let eq_system = self.vec_of_eq.clone();
        let values = vec!["T".to_string(), "qs".to_string()];
        self.vec_of_vars = values.clone();
        let T0s = self.T0s;
        let arg = "x".to_string();

        let mut BorderConditions = HashMap::new();
        BorderConditions.insert("T".to_string(), vec![(0, T0s)]);
        BorderConditions.insert(
            "qs".to_string(),
            vec![(1, self.Pe * (self.L / self.dT) * 1e-10)],
        );
        let n = 100;
        let t0 = 0.0;
        let t_end = 1.0;
        let ones = vec![1.0; values.len() * n];
        let initial_guess: DMatrix<f64> =
            DMatrix::from_column_slice(values.len(), n, DVector::from_vec(ones).as_slice());

        let mut bvp_solver = BVPwrap::new(
            None,
            Some(t0),
            Some(t_end),
            Some(n),
            eq_system,
            values,
            vec![],
            None,
            BorderConditions,
            arg,
            1e-8,
            100000,
            initial_guess,
        );
        bvp_solver.solve();
        let sol = bvp_solver.get_result();
        match sol {
            Some(y_mesh) => {
                let dT = self.dT;
                let x_mesh = self.L * bvp_solver.x_mesh.clone();
                let T_sol: DVector<f64> = dT * y_mesh.column(0).into_owned();
                let qs_sol: DVector<f64> =
                    (self.dT / self.L) * y_mesh.column(1).into_owned() / self.Pe;
                self.x_mesh = x_mesh;
                self.T_sol = T_sol;
                self.qs_sol = qs_sol;
                println!("BVP solved successfully with Sci method");
                //println!("x_mesh = {:?}", self.x_mesh);
                println!("T_sol = {:?}", self.T_sol);
                println!("qs_sol = {:?}", self.qs_sol);
            }
            None => {
                println!("Error during BVP solving with Sci method");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zeldovich() {
        let Q_g = 3000.0; //#350*4.184#?
        let C_p = 0.35 * 4.184;
        let Lambda_eff = 0.07; // W/m-K 
        let M = 34.2 / 1000.0; // # kg/mol
        // let B_g  =  10.0_f64.powf(12.5); // c^-1
        // let E =  39_500.0*4.184;

        let B_g = 1.3 * 1e5;
        let E = 5000.0 * 4.184;
        ///////////////////////////////////////
        let P = 40.0;
        let m = 0.7 * 10.0;
        let L = 1e-3;
        let T0 = 450.0;
        let mut z = Zeldovich::new();
        z.set_params(Q_g, C_p, Lambda_eff, M, B_g, E, 1.0, 1.0);
        let P = P * 101325.0;
        z.set_problem(P, m, L, T0);
        z.calc_constants();
        for eq_i in z.vec_of_eq.iter() {
            println!("eq_i = {}", eq_i);
        }
        z.solve_bvp();
    }
}
