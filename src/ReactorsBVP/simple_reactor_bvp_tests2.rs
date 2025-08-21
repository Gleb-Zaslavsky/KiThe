#[cfg(test)]
mod tests {
    use super::super::SimpleReactorBVP::*;
    use RustedSciThe::{
         symbolic::symbolic_engine::Expr,
    };
    use nalgebra::DMatrix;
    use std::{collections::HashMap, vec};
use crate::ReactorsBVP::reactor_BVP_utils::{
   
    ScalingConfig,
   
};
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
        // println!("reactor {:?}", reactor.kindata);
        ///////////// parameters //////////////////
        let ro0 = M * P / (R_G * T0);
        let D = Lambda_eff / (C_p * ro0);
        println!("D = {}", D);
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
        let _ = &reactor.Le_number();
        println!("reactor {:?}", reactor.kindata);
        match res {
            Ok(_) => println!("ok"),
            Err(e) => println!("error {:?}", e),
        }

        let n_steps = 200;
        let ig = vec![0.99; n_steps * reactor.solver.unknowns.len()];
        let initial_guess = DMatrix::from_vec(reactor.solver.unknowns.len(), n_steps, ig);
        let _ = reactor
            .solver
            .solve_BVPsci(initial_guess, n_steps, n_steps * 2000, 1e-4);
        reactor.postprocessing();
        //  reactor.gnuplot();
        println!("BC {:?}", &reactor.solver.BorderConditions);
        println!("unknowns {:?}", &reactor.solver.unknowns);
        let rates = reactor.map_eq_rate.clone();
        for (eq, rate) in rates {
            println!("reaction {} rate {}", eq, rate);
        }
        println!("\n \n");
        let system = reactor.map_of_equations.clone();
        for (subs, (variable, eq)) in system {
            println!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = &reactor.solver.BorderConditions;
        println!("bc {:?}", bc);
        println!(" unknowns{:?}", reactor.solver.unknowns);
        // let sol = &reactor.solver.solution.clone().unwrap();
        //let T = sol.column(0).clone();
        //  println!("T = {}", T);
        //  println!("q = {}", sol.column(1));
        // println!("C0 = {}", sol.column(2));
    }
    #[test]
    fn test_direct_eq() {
        use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::NRBVP;
        use prettytable::{Table, row};
        // variables
        let unknowns_str: Vec<&str> = vec!["Teta", "q", "C0", "J0", "C1", "J1"];
        let unhnowns_Str:Vec<String> = unknowns_str.iter().map(|s| s.to_string()).collect();
        let unknowns: Vec<Expr> = Expr::parse_vector_expression(unknowns_str);
        let Teta = unknowns[0].clone();
        let q = unknowns[1].clone();
        let C0 = unknowns[2].clone();
        let J0 = unknowns[3].clone();

        let J1 = unknowns[5].clone();
        // Parameters
        let Q = 3000.0 * 1e3 * 0.034;;
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

        println!("\n=== COEFFICIENTS ===");
        coeff_table.printstd();

        // Pretty print equations table
        let mut eq_table = Table::new();
        eq_table.add_row(row!["Unknown", "Equation"]);

        for (unknown, equation) in eq_and_unknowns {
            eq_table.add_row(row![format!("{}", unknown), format!("{}", equation)]);
        }

        println!("\n=== EQUATIONS SYSTEM ===");
        eq_table.printstd();
    

        ////////////////////////////////////////////
        //solver

            let Teta_initial = (T_initial - dT) / T_scale;
        let BoundaryConditions = HashMap::from([
            ("Teta".to_string(), (0, Teta_initial)),
            ("q".to_string(), (1, 1e-10)),
            ("C0".to_string(), (0, C1_0)),
            ("J0".to_string(), (1, 1e-7)),
            ("C1".to_string(), (0, 1e-3)),
            ("J1".to_string(), (1, 1e-7)),
        ]);

        let Bounds = HashMap::from([
            ("Teta".to_string(), (-100.0, 100.0)),
            ("q".to_string(), (-1e20, 1e20)),
            ("C0".to_string(), (0.0, 1.0)),
            ("J0".to_string(), (-1e20, 1e20)),
            ("C1".to_string(), (0.0, 1.0)),
            ("J1".to_string(), (-1e20, 1e20)),
        ]);
        let n_steps = 400;
        let strategy_params = Some(HashMap::from([
            ("max_jac".to_string(), None),
            ("maxDampIter".to_string(), None),
            ("DampFacor".to_string(), None),
            ("adaptive".to_string(), None),
        ]));
        let rel_tolerance =HashMap::from([
            ("Teta".to_string(), 1e-3),
            ("q".to_string(), 1e-3),
            ("C0".to_string(), 1e-3),
            ("J0".to_string(), 1e-3),
            ("C1".to_string(), 1e-3),
            ("J1".to_string(), 1e-3),
        ]);
                let ig = vec![0.99; n_steps * unknowns.len()];
        let initial_guess = DMatrix::from_vec(unknowns.len(), n_steps, ig);
        let max_iterations = 100;
        let abs_tolerance = 1e-6;
        let loglevel = Some("info".to_string());
        let scheme = "forward".to_string();
        let method = "Sparse".to_string();
        let strategy = "Damped".to_string();
        let linear_sys_method = None;

        // Using the new tolerance helper - much simpler!
        let mut bvp = NRBVP::new(
            eqs.clone(),
            initial_guess,
             unhnowns_Str,
            "x".to_string(),
            BoundaryConditions,
            0.0,
            1.0,
            n_steps,
            scheme,
            strategy,
            strategy_params,
            linear_sys_method,
            method,
            abs_tolerance,
            Some(rel_tolerance),
            max_iterations,
            Some(Bounds),
            loglevel,
        );
        bvp.solve();
        bvp.gnuplot_result();
        let eq_and_unknowns = unknowns.clone().into_iter().zip(eqs.clone());
        for (unknown, equation) in  eq_and_unknowns{
            println!("unknown: {} | equation: {}", unknown, equation);
        }
  
    }
}
