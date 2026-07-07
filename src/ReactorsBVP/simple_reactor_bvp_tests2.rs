#[cfg(test)]
mod tests {
    use super::super::SimpleReactorBVP::*;
    use crate::ReactorsBVP::solver_backend::ReactorBvpSolverConfig;
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use nalgebra::DMatrix;

    use crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig;
    use crate::ReactorsBVP::reactor_BVP_utils::{BoundsConfig, ToleranceConfig};
    use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{
        AdaptiveGridConfig, SolverParams,
    };
    use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;
    use log::info;
    use std::{collections::HashMap, vec};

    /// Build a damped-solver option block from the new KiThe facade while preserving
    /// the legacy test knobs that still matter for the direct NRBVP smoke cases.
    fn direct_default_solver_options(
        scheme: &str,
        strategy: &str,
        strategy_params: Option<SolverParams>,
        linear_sys_method: Option<String>,
        method: &str,
        abs_tolerance: f64,
        rel_tolerance: Option<HashMap<String, f64>>,
        max_iterations: usize,
        bounds: Option<HashMap<String, (f64, f64)>>,
        loglevel: Option<String>,
    ) -> RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::DampedSolverOptions {
        let mut options = ReactorBvpSolverConfig::default_lambdify()
            .to_rusted_options()
            .expect("default reactor BVP options should build");
        options = options
            .with_scheme_name(scheme.to_string())
            .with_strategy_params(strategy_params)
            .with_abs_tolerance(abs_tolerance)
            .with_max_iterations(max_iterations)
            .with_loglevel(loglevel);
        options.strategy = strategy.to_string();
        options.linear_sys_method = linear_sys_method;
        options.method = method.to_string();
        if let Some(rel_tolerance) = rel_tolerance {
            options = options.with_rel_tolerance(rel_tolerance);
        }
        if let Some(bounds) = bounds {
            options = options.with_bounds(bounds);
        }
        options
    }

    /// Solve the compact HMX fixture through the new default reactor facade.
    ///
    /// This helper keeps the old smoke tests focused on the combustion setup
    /// while making the backend path canonical: lambdify execution, AtomView
    /// symbolic assembly, and the banded matrix route.
    fn solve_hmx_with_default_backend(
        reactor: &mut SimpleReactorTask,
        initial_guess: DMatrix<f64>,
        n_steps: usize,
        strategy_params: Option<SolverParams>,
        tolerance_config: ToleranceConfig,
        bounds_config: BoundsConfig,
        max_iterations: usize,
        abs_tolerance: f64,
        loglevel: Option<String>,
    ) -> Result<(), ReactorError> {
        let substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        reactor.solver.solve_NRBVP_with_configs(
            initial_guess,
            n_steps,
            "forward".to_string(),
            "Damped".to_string(),
            strategy_params,
            None,
            "Banded".to_string(),
            abs_tolerance,
            tolerance_config,
            bounds_config,
            &substances,
            max_iterations,
            loglevel,
        )
    }
    #[test]
    fn hmx_test3() {
        /////////////////// setting up kinetics
        let eq = "HMX=>HMXprod".to_string();
        let Q_g = 0.0; // Дж/кг
        let C_p = 0.35 * 4.184 * 1000.0;
        let Lambda_eff = 0.07; // W/m-K 
        let n = 0.0;
        let M = 34.2 / 1000.0; //  kg/mol
        let A = 1.3e5; // 1.3 * 1e5; //
        let E = 5000.0 * 4.184; // 5000.0 * 4.184; //  
        let L = 3e-4;
        let T0 = 1000.0;
        let P: f64 = 2e6; // Pa
        let Tm = 1500.0;
        let hmx = HashMap::from([
            ("H".to_string(), 4),
            ("N".to_string(), 8),
            ("C".to_string(), 8),
            ("O".to_string(), 8),
        ]);
        let hmxprod = HashMap::from([
            ("H".to_string(), 6),
            ("C".to_string(), 1),
            ("O".to_string(), 1),
        ]);
        let groups = Some(HashMap::from([
            ("HMX".to_string(), hmx.clone()),
            ("HMXprod".to_string(), hmxprod),
        ]));
        //////////// instance of problem constructor/////////////
        let mut reactor = SimpleReactorTask::new();
        let struct_with_params = FastElemReact {
            eq,
            A,
            n,
            E,
            Q: Q_g,
        };
        let vec_of_structs = vec![struct_with_params];

        let _ = reactor.fast_react_set(vec_of_structs);
        reactor.kindata.substances = vec!["HMX".to_string(), "HMXprod".to_string()];
        reactor.kindata.groups = groups;
        // info!("reactor {:?}", reactor.kindata);
        ///////////// parameters //////////////////
        let ro0 = M * P / (R_G * T0);
        let D = Lambda_eff / (C_p * ro0);
        info!("D = {}", D);
        let Diffusion = HashMap::from([("HMX".to_string(), D), ("HMXprod".to_string(), D)]);
        let boundary_condition = HashMap::from([
            ("HMX".to_string(), 1.0 - 1e-3),
            ("HMXprod".to_string(), 1e-3),
            ("T".to_string(), T0),
        ]);
        let thermal_effects = vec![Q_g];

        let m = 0.077 * (P / 1e5).powf(0.748); //
        let scaling = ScalingConfig::new(800.0, L, 800.0);

        reactor.set_parameters(
            thermal_effects,
            P,
            Tm,
            C_p,
            boundary_condition,
            Lambda_eff,
            Diffusion,
            m,
            scaling,
        );

        reactor.M = M;
        let res = reactor.setup_bvp();
        let _ = reactor.Le_number();
        info!("reactor {:?}", reactor.kindata);
        match res {
            Ok(_) => info!("ok"),
            Err(e) => info!("error {:?}", e),
        }

        let n_steps = 200;
        let ig = vec![0.99; n_steps * reactor.solver.unknowns.len()];
        let initial_guess = DMatrix::from_vec(reactor.solver.unknowns.len(), n_steps, ig);
        let tolerance_config = ToleranceConfig::new(1e-5, 1e-5, 1e-5, 1e-6);
        let bounds_config = BoundsConfig::new(
            (-0.1, 1.1),     // C bounds
            (-1e20, 1e20),   // J bounds
            (-100.0, 100.0), // Teta bounds
            (-1e20, 1e20),   // q bounds
        );
        let _ = solve_hmx_with_default_backend(
            &mut reactor,
            initial_guess,
            n_steps,
            None,
            tolerance_config,
            bounds_config,
            n_steps * 2000,
            1e-4,
            Some("info".to_string()),
        );
        reactor.postprocessing().unwrap();
        //  reactor.gnuplot();
        info!("BC {:?}", &reactor.solver.BorderConditions);
        info!("unknowns {:?}", &reactor.solver.unknowns);
        let rates = reactor.map_eq_rate.clone();
        for (eq, rate) in rates {
            info!("reaction {} rate {}", eq, rate);
        }
        info!("\n \n");
        let system = reactor.map_of_equations.clone();
        for (subs, (variable, eq)) in system {
            info!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = &reactor.solver.BorderConditions;
        info!("bc {:?}", bc);
        info!(" unknowns{:?}", reactor.solver.unknowns);
        // let sol = &reactor.solver.solution.clone().unwrap();
        //let T = sol.column(0).clone();
        //  info!("T = {}", T);
        //  info!("q = {}", sol.column(1));
        // info!("C0 = {}", sol.column(2));
    }
    #[test]
    fn test_direct_eq() {
        use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::NRBVP;
        use prettytable::{Table, row};
        // variables
        let unknowns_str: Vec<&str> = vec!["Teta", "q", "C0", "J0", "C1", "J1"];
        let unhnowns_Str: Vec<String> = unknowns_str.iter().map(|s| s.to_string()).collect();
        let unknowns: Vec<Expr> = Expr::parse_vector_expression(unknowns_str);
        let Teta = unknowns[0].clone();
        let q = unknowns[1].clone();
        let C0 = unknowns[2].clone();
        let J0 = unknowns[3].clone();

        let J1 = unknowns[5].clone();
        // Parameters
        let Q = 3000.0 * 1e3 * 0.034;
        let dT = 600.0;
        let T_scale = 600.0;
        let L = 1e-4;
        let M0 = 34.2 / 1000.0;
        let Lambda = 0.07;
        let Cp = 0.35 * 4.184 * 1000.0;
        let P = 2e6;
        let T0 = 800.0;
        let Tm = 1500.0;
        let C1_0 = 1.0;
        let T_initial = 1000.0;
        // problem settings

        // coefficients
        let m = 0.077 * (P / 1e5 as f64).powf(0.748) / 100.0;
        let ro0 = M0 * P / (R_G * T0);
        let D = Lambda / (Cp * ro0);
        let Pe_q = (L * Cp * m) / Lambda;

        let D_ro = D * ro0 * (Tm / T0).powf(0.5);
        let Pe_D = m * L / D_ro;
        // conversion to sym
        let dT_sym = Expr::Const(dT);
        let T_scale_sym = Expr::Const(T_scale);

        let Lambda_sym = Expr::Const(Lambda);
        let Q = Expr::Const(Q);
        let A = Expr::Const(1.3e5);
        let E = Expr::Const(5000.0 * 4.184);
        let M = Expr::Const(M0);
        let R_g = Expr::Const(8.314);
        let ro_m = Expr::Const(M0 * P / (8.314 * Tm));
        let qm = Expr::Const(L.powf(2.0) / T_scale);
        let qs = Expr::Const(L.powf(2.0));
        let Pe_q_sym = Expr::Const(Pe_q);
        let ro_D = Expr::Const(D_ro);
        let ro_D = vec![ro_D.clone(), ro_D.clone()];
        let Pe_D = vec![Expr::Const(Pe_D), Expr::Const(Pe_D)];
        let minus = Expr::Const(-1.0);
        // EQ SYSTEM

        let Rate = A
            * Expr::exp(-E / (R_g * (Teta * T_scale_sym + dT_sym)))
            * C0
            * (ro_m.clone() / M.clone());
        let eq_T = q.clone() / Lambda_sym;
        let eq_q = q * Pe_q_sym - Q * Rate.clone() * qm;
        let eq_C0 = J0.clone() / ro_D[0].clone();
        let eq_J0 = J0 * Pe_D[0].clone()
            - (M.clone() * minus * Rate.clone() * ro_m.clone() / M.clone()) * qs.clone();
        let eq_C1 = J1.clone() / ro_D[1].clone();
        let eq_J1 = J1 * Pe_D[1].clone() - (M.clone() * Rate * ro_m / M) * qs;

        let eqs = vec![eq_T, eq_q, eq_C0, eq_J0, eq_C1, eq_J1];
        let eq_and_unknowns = unknowns.clone().into_iter().zip(eqs.clone());

        // Pretty print coefficients table
        let mut coeff_table = Table::new();
        coeff_table.add_row(row!["Parameter", "Value"]);
        coeff_table.add_row(row!["Q", format!("{:.2e}", 102000.0)]);
        coeff_table.add_row(row!["dT", format!("{:.1}", 600.0)]);
        coeff_table.add_row(row!["L", format!("{:.2e}", L)]);
        coeff_table.add_row(row!["M0", format!("{:.4}", M0)]);
        coeff_table.add_row(row!["Lambda", format!("{:.3}", Lambda)]);
        coeff_table.add_row(row!["Pe_q", format!("{:.4e}", Pe_q)]);
        coeff_table.add_row(row!["Pe_D", format!("{:?}", Pe_D)]);
        coeff_table.add_row(row!["D_ro", format!("{:.2e}", D_ro)]);

        info!("\n=== COEFFICIENTS ===");
        info!("{}", coeff_table);

        // Pretty print equations table
        let mut eq_table = Table::new();
        eq_table.add_row(row!["Unknown", "Equation"]);

        for (unknown, equation) in eq_and_unknowns {
            eq_table.add_row(row![format!("{}", unknown), format!("{}", equation)]);
        }
        info!("{}", eq_table);

        ////////////////////////////////////////////
        //solver

        let Teta_initial = (T_initial - dT) / T_scale;
        let BoundaryConditions = HashMap::from([
            ("Teta".to_string(), vec![(0, Teta_initial)]),
            ("q".to_string(), vec![(1, 1e-10)]),
            ("C0".to_string(), vec![(0, C1_0)]),
            ("J0".to_string(), vec![(1, 1e-7)]),
            ("C1".to_string(), vec![(0, 1e-3)]),
            ("J1".to_string(), vec![(1, 1e-7)]),
        ]);

        let Bounds = HashMap::from([
            ("Teta".to_string(), (-100.0, 100.0)),
            ("q".to_string(), (-1e20, 1e20)),
            ("C0".to_string(), (-1.0, 1.5)),
            ("J0".to_string(), (-1e20, 1e20)),
            ("C1".to_string(), (-1.0, 1.5)),
            ("J1".to_string(), (-1e20, 1e20)),
        ]);
        let n_steps = 30;
        let grid_method = GridRefinementMethod::GrcarSmooke(0.01, 0.01, 1.5);
        // or GridRefinementMethod::Pearson(0.05, 2.5);
        let adaptive = AdaptiveGridConfig {
            version: 1,
            max_refinements: 3,
            grid_method,
        };
        let strategy_params = SolverParams {
            max_jac: Some(5),
            max_damp_iter: Some(5),
            damp_factor: Some(0.5),
            adaptive: Some(adaptive),
        };

        let rel_tolerance = HashMap::from([
            ("Teta".to_string(), 1e-5),
            ("q".to_string(), 1e-5),
            ("C0".to_string(), 1e-5),
            ("J0".to_string(), 1e-5),
            ("C1".to_string(), 1e-5),
            ("J1".to_string(), 1e-5),
        ]);
        let ig = vec![0.99; n_steps * unknowns.len()];
        let initial_guess = DMatrix::from_vec(unknowns.len(), n_steps, ig);
        let max_iterations = 100;
        let abs_tolerance = 1e-6;
        let loglevel = Some("info".to_string());
        let scheme = "forward".to_string();
        let method = "Banded".to_string();
        let strategy = "Damped".to_string();
        let linear_sys_method = None;

        // Using the new tolerance helper - much simpler!
        let options = direct_default_solver_options(
            &scheme,
            &strategy,
            Some(strategy_params),
            linear_sys_method,
            &method,
            abs_tolerance,
            Some(rel_tolerance),
            max_iterations,
            Some(Bounds),
            loglevel,
        );

        let mut bvp = NRBVP::new_with_options(
            eqs.clone(),
            initial_guess,
            unhnowns_Str,
            "x".to_string(),
            BoundaryConditions,
            0.0,
            1.0,
            n_steps,
            options,
        );
        bvp.dont_save_log(false);
        bvp.solve();
        bvp.gnuplot_result();
        let eq_and_unknowns = unknowns.clone().into_iter().zip(eqs.clone());
        for (unknown, equation) in eq_and_unknowns {
            info!("unknown: {} | equation: {}", unknown, equation);
        }
        info!("\n=== EQUATIONS SYSTEM ===");
        info!("{}", coeff_table);
    }

    #[test]
    fn test_direct_eq2() {
        use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::NRBVP;
        use prettytable::{Table, row};
        // variables
        let unknowns_str: Vec<&str> = vec!["Teta", "q", "C0", "J0", "C1", "J1"];
        let unhnowns_Str: Vec<String> = unknowns_str.iter().map(|s| s.to_string()).collect();
        let unknowns: Vec<Expr> = Expr::parse_vector_expression(unknowns_str);
        let Teta = unknowns[0].clone();
        let q = unknowns[1].clone();
        let C0 = unknowns[2].clone();
        let J0 = unknowns[3].clone();

        let J1 = unknowns[5].clone();
        // Parameters
        let Q = 3000.0 * 1e3 * 0.034;
        let dT = 600.0;
        let T_scale = 600.0;
        let L: f64 = 1e-4;
        let M0 = 34.2 / 1000.0;
        let Lambda = 0.07;

        let P = 2e6;

        let Tm = 1500.0;
        let C1_0 = 1.0;
        let T_initial = 1000.0;
        // problem settings

        // coefficients

        let Pe_q = 0.0090168;

        let D_ro = 2.88e-4;
        let Pe_D = 1.50e-3;
        let ro_m_ = M0 * P / (8.314 * Tm);
        // conversion to sym
        let dT_sym = Expr::Const(dT);
        let T_scale_sym = Expr::Const(T_scale);

        let Lambda_sym = Expr::Const(Lambda);
        let Q = Expr::Const(Q);
        let A = Expr::Const(1.3e5);
        let E = Expr::Const(5000.0 * 4.184);
        let M = Expr::Const(M0);
        let R_g = Expr::Const(8.314);
        let ro_m = Expr::Const(ro_m_);
        let qm = Expr::Const(L.powf(2.0) / T_scale);
        let qs = Expr::Const(L.powf(2.0));
        let Pe_q_sym = Expr::Const(Pe_q);
        let ro_D = Expr::Const(D_ro);
        let ro_D = vec![ro_D.clone(), ro_D.clone()];
        let Pe_D = vec![Expr::Const(Pe_D), Expr::Const(Pe_D)];
        let minus = Expr::Const(-1.0);
        let M_reag = Expr::Const(0.342);
        // EQ SYSTEM

        let Rate = A
            * Expr::exp(-E / (R_g * (Teta * T_scale_sym + dT_sym)))
            * C0
            * (ro_m.clone() / M_reag.clone());
        let eq_T = q.clone() / Lambda_sym;
        let eq_q = q * Pe_q_sym - Q * Rate.clone() * qm;
        let eq_C0 = J0.clone() / ro_D[0].clone();
        let eq_J0 = J0 * Pe_D[0].clone()
            - (M.clone() * minus * Rate.clone() * ro_m.clone() / M.clone()) * qs.clone();
        let eq_C1 = J1.clone() / ro_D[1].clone();
        let eq_J1 = J1 * Pe_D[1].clone() - (M.clone() * Rate * ro_m / M) * qs;

        let eqs = vec![eq_T, eq_q, eq_C0, eq_J0, eq_C1, eq_J1];
        let eq_and_unknowns = unknowns.clone().into_iter().zip(eqs.clone());

        // Pretty print coefficients table
        let mut coeff_table = Table::new();
        coeff_table.add_row(row!["Parameter", "Value"]);
        coeff_table.add_row(row!["Q", format!("{:.2e}", 102000.0)]);
        coeff_table.add_row(row!["dT", format!("{:.1}", 600.0)]);
        coeff_table.add_row(row!["L", format!("{:.2e}", L)]);
        coeff_table.add_row(row!["M0", format!("{:.4}", M0)]);
        coeff_table.add_row(row!["Lambda", format!("{:.3}", Lambda)]);
        coeff_table.add_row(row!["Pe_q", format!("{:.4e}", Pe_q)]);
        coeff_table.add_row(row!["Pe_D", format!("{:?}", Pe_D)]);
        coeff_table.add_row(row!["D_ro", format!("{:.2e}", D_ro)]);
        coeff_table.add_row(row!["ro", format!("{:.2e}", ro_m_)]);

        info!("\n=== COEFFICIENTS ===");
        info!("{}", coeff_table);

        // Pretty print equations table
        let mut eq_table = Table::new();
        eq_table.add_row(row!["Unknown", "Equation"]);

        for (unknown, equation) in eq_and_unknowns {
            eq_table.add_row(row![format!("{}", unknown), format!("{}", equation)]);
        }
        info!("{}", eq_table);

        ////////////////////////////////////////////
        //solver

        let Teta_initial = (T_initial - dT) / T_scale;
        let BoundaryConditions = HashMap::from([
            ("Teta".to_string(), vec![(0, Teta_initial)]),
            ("q".to_string(), vec![(1, 1e-10)]),
            ("C0".to_string(), vec![(0, C1_0)]),
            ("J0".to_string(), vec![(1, 1e-7)]),
            ("C1".to_string(), vec![(0, 1e-3)]),
            ("J1".to_string(), vec![(1, 1e-7)]),
        ]);

        let Bounds = HashMap::from([
            ("Teta".to_string(), (-100.0, 100.0)),
            ("q".to_string(), (-1e20, 1e20)),
            ("C0".to_string(), (-1.0, 1.5)),
            ("J0".to_string(), (-1e20, 1e20)),
            ("C1".to_string(), (-1.0, 1.5)),
            ("J1".to_string(), (-1e20, 1e20)),
        ]);
        let n_steps = 50;
        let grid_method = GridRefinementMethod::GrcarSmooke(0.1, 0.1, 2.5);
        // or GridRefinementMethod::Pearson(0.05, 2.5);
        let adaptive = AdaptiveGridConfig {
            version: 1,
            max_refinements: 3,
            grid_method,
        };
        let strategy_params = SolverParams {
            max_jac: Some(5),
            max_damp_iter: Some(5),
            damp_factor: Some(0.5),
            adaptive: Some(adaptive),
        };

        let rel_tolerance = HashMap::from([
            ("Teta".to_string(), 1e-5),
            ("q".to_string(), 1e-5),
            ("C0".to_string(), 1e-5),
            ("J0".to_string(), 1e-5),
            ("C1".to_string(), 1e-5),
            ("J1".to_string(), 1e-5),
        ]);
        let ig = vec![0.99; n_steps * unknowns.len()];
        let initial_guess = DMatrix::from_vec(unknowns.len(), n_steps, ig);
        let max_iterations = 100;
        let abs_tolerance = 1e-6;
        let loglevel = Some("info".to_string());
        let scheme = "forward".to_string();
        let method = "Banded".to_string();
        let strategy = "Damped".to_string();
        let linear_sys_method = None;

        // Using the new tolerance helper - much simpler!
        let options = direct_default_solver_options(
            &scheme,
            &strategy,
            Some(strategy_params),
            linear_sys_method,
            &method,
            abs_tolerance,
            Some(rel_tolerance),
            max_iterations,
            Some(Bounds),
            loglevel,
        );

        let mut bvp = NRBVP::new_with_options(
            eqs.clone(),
            initial_guess,
            unhnowns_Str,
            "x".to_string(),
            BoundaryConditions,
            0.0,
            1.0,
            n_steps,
            options,
        );
        bvp.dont_save_log(false);
        bvp.solve();
        bvp.gnuplot_result();
        let eq_and_unknowns = unknowns.clone().into_iter().zip(eqs.clone());
        for (unknown, equation) in eq_and_unknowns {
            info!("unknown: {} | equation: {}", unknown, equation);
        }
        info!("\n=== EQUATIONS SYSTEM ===");
        info!("{}", coeff_table);
    }
}
