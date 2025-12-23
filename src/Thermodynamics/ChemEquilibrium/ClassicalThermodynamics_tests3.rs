#[cfg(test)]
mod tests {
    use RustedSciThe::numerical::Nonlinear_systems::NR::Method;

    use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::{
        ThermodynamicCalculations, Thermodynamics,
    };
    use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamicsSolver::{
        LMParams, NRParams, SolverType,
    };
    use crate::Thermodynamics::User_PhaseOrSolution::{
        CustomSubstance, SubstanceSystemFactory, SubstancesContainer,
    };
    use crate::Thermodynamics::User_substances::{DataType, LibraryPriority, Phases, SubsData};
    use RustedSciThe::symbolic::symbolic_engine::Expr;
    use approx::assert_relative_eq;
    use core::panic;
    use std::time::Instant;

    use crate::Thermodynamics::ChemEquilibrium::easy_equilibrium::EasyEquilibrium;
    use std::collections::HashMap;
    use std::vec;
    /////////////////DISSOCIATION TESTS ///////////////////////
    pub fn set_equilibrium_system_dissociation(
        subs: Vec<String>,
        T: f64,
    ) -> Box<dyn Fn(f64) -> f64> {
        assert!(subs.len() == 2, "This test is designed for two substances.");
        let mut subdata = SubsData::new();

        // Set library priorities
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);

        // Set phases
        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();

        let _ = subdata.search_substances();
        let _ = subdata.extract_thermal_coeffs(subs[0].as_str(), T);
        let _ = subdata.extract_thermal_coeffs(subs[1].as_str(), T);
        let _ = subdata.calculate_therm_map_of_sym();

        let A2 = subdata.therm_map_of_sym.get(&subs[0].clone()).unwrap();
        let S_A2 = *A2.get(&DataType::dS_sym).clone().unwrap().clone().unwrap();
        let H_A2 = *A2.get(&DataType::dH_sym).clone().unwrap().clone().unwrap();
        let T = Expr::Var("T".to_string());
        let dG_A2 = H_A2 - T.clone() * S_A2;
        println!("dG_A2: {}", dG_A2);
        let A = subdata.therm_map_of_sym.get(&subs[1].clone()).unwrap();
        let S_A = *A.get(&DataType::dS_sym).clone().unwrap().clone().unwrap();
        let H_A = *A.get(&DataType::dH_sym).clone().unwrap().clone().unwrap();
        let dG_A = H_A - T * S_A;
        println!("dG_A2: {}", dG_A2);
        let dG_mix = Expr::Const(2.0) * dG_A - dG_A2;
        let dG_mix_simplified = dG_mix.simplify();
        let dG_mix = dG_mix_simplified.lambdify1D();
        let K = Box::new(move |T| f64::exp(-dG_mix(T) / (8.31446261815324 * T)));
        K
    }

    /// Compute equilibrium composition for A2 <=> 2 A (ideal gases)
    /// Kp: equilibrium constant in PA (Pa).  -- NOTE: Barklem's table gives log10(pK) where pK is in Pa.
    /// P:  total system pressure in PA (Pa) -- be consistent with Kp units.
    ///
    /// Returns (x_A, x_A2) mole fractions (sum ≈ 1).
    pub fn dissociation_a2_to_2a_from_kp_pa(kp_pa: f64, p_pa: f64) -> (f64, f64) {
        // From algebra in chat: alpha^2 = Kp / (Kp + 4 P)
        // where alpha = degree of dissociation of initial 1 mol A2.
        let denom = kp_pa + 4.0 * p_pa;
        // protect against tiny negative noise
        let ratio = if denom > 0.0 { kp_pa / denom } else { 0.0 };
        let alpha = if ratio > 0.0 { ratio.sqrt() } else { 0.0 };
        let x_a = 2.0 * alpha / (1.0 + alpha);
        let x_a2 = (1.0 - alpha) / (1.0 + alpha);
        (x_a, x_a2)
    }
    #[test]
    fn test_dissociation_a2_to_2a() {
        let subs = vec!["O2".to_string(), "O".to_string()];

        let temperatures = vec![1000.0, 2000.0, 3000.0, 4000.0, 5000.0];
        let pressure_pa = 101325.0; // 1 atm in Pa

        for &T in &temperatures {
            let kp_func = set_equilibrium_system_dissociation(subs.clone(), T);
            let kp_pa = kp_func(T) * pressure_pa; // Convert K to Kp in Pa
            let (x_o, x_o2) = dissociation_a2_to_2a_from_kp_pa(kp_pa, pressure_pa);
            println!(
                "   \n \n T: {} K, Kp: {:.3e} Pa => x_O: {:.5}, x_O2: {:.5}",
                T, kp_pa, x_o, x_o2
            );
            assert!(
                (x_o + x_o2 - 1.0).abs() < 1e-8,
                "Mole fractions do not sum to 1"
            );
        }
    }

    fn calculate_A2_2A_decomposition(T: f64, gas_subs: Vec<String>) -> (f64, f64) {
        let P = 101325.0;

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
            CustomSubstance::OnePhase(onephase) => {
                let subs = &onephase.subs_data.substances;
                assert_eq!(subs.len(), 2);
            }
            _ => panic!(),
        }

        let mut n = HashMap::new();
        let map_of_gas = HashMap::from([(gas_subs.clone()[0].clone(), 1.0)]);

        n.insert(Some("gas".to_string()), (Some(1.0), Some(map_of_gas)));

        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(Some(T), P, None, None);
        assert!(td.is_ok());

        let mut td = td.unwrap();
        td.Tm = 1e3;
        td.set_number_of_moles(Some(n)).unwrap();
        td.initial_composition().unwrap();
        td.create_indexed_variables();
        td.create_full_system_of_equations_with_const_T().unwrap();

        let params = LMParams {
            initial_guess: Some(vec![0.1, 0.1, 1.0, 0.1]),
            tolerance: Some(1e-6),
            max_iterations: Some(1000),
            loglevel: Some("info".to_string()),
        };
        td.set_solver_type(SolverType::LevenbergMarquardt(params));
        td.generate_eqs();
        td.solve_for_T(T);
        let sol = td.solver.get_solution();
        let n_O2 = *sol.get("N0").unwrap();
        let n_O = *sol.get("N1").unwrap();
        let x_O2 = n_O2 / (n_O2 + n_O);
        let x_O = n_O / (n_O2 + n_O);

        (x_O2, x_O)
    }

    #[test]
    fn test_calculate_O2_2O_decomposition() {
        let gas_subs = vec!["O2".to_string(), "O".to_string()];
        let temperatures = vec![3000.0, 4000.0, 5000.0];
        let pressure_pa = 101325.0; // 1 atm in Pa
        for &T in &temperatures {
            println!("\n Calculating at T = {} K", T);
            let (x_O2, x_O) = calculate_A2_2A_decomposition(T, gas_subs.clone());

            let kp_func = set_equilibrium_system_dissociation(gas_subs.clone(), T);
            let kp_pa = kp_func(T) * pressure_pa; // Convert K to Kp in Pa
            let (x_o_direct, x_o2_direct) = dissociation_a2_to_2a_from_kp_pa(kp_pa, pressure_pa);
            println!(
                "   \n T {} K, from Kp => x_O: {:.5}, x_O2: {:.5}  (from Kp) ",
                T, x_o_direct, x_o2_direct
            );
            println!(
                "   \n T {} K, from system => x_O: {:.5}, x_O2: {:.5}  (from Kp) ",
                T, x_O, x_O2
            );
            assert_relative_eq!(x_O, x_o_direct, epsilon = 1e-3);
        }
    }
    #[test]
    fn test_calculate_N2_2N_decomposition() {
        let gas_subs = vec!["N2".to_string(), "N".to_string()];
        let temperatures = vec![3000.0, 4000.0, 5000.0];
        let pressure_pa = 101325.0; // 1 atm in Pa
        for &T in &temperatures {
            println!("\n Calculating at T = {} K", T);
            let (x_O2, x_O) = calculate_A2_2A_decomposition(T, gas_subs.clone());

            let kp_func = set_equilibrium_system_dissociation(gas_subs.clone(), T);
            let kp_pa = kp_func(T) * pressure_pa; // Convert K to Kp in Pa
            let (x_o_direct, x_o2_direct) = dissociation_a2_to_2a_from_kp_pa(kp_pa, pressure_pa);
            println!(
                "   \n T {} K, from Kp => x_O: {:.5}, x_O2: {:.5}  (from Kp) ",
                T, x_o_direct, x_o2_direct
            );
            println!(
                "   \n T {} K, from system => x_O: {:.5}, x_O2: {:.5}  (from Kp) ",
                T, x_O, x_O2
            );
            assert_relative_eq!(x_O, x_o_direct, epsilon = 1e-3);
        }
    }

    #[test]
    fn test_calculate_H2_2H_decomposition() {
        let gas_subs = vec!["H2".to_string(), "H".to_string()];
        let temperatures = vec![3000.0, 4000.0, 5000.0];
        let pressure_pa = 101325.0; // 1 atm in Pa
        for &T in &temperatures {
            println!("\n Calculating at T = {} K", T);
            let (x_O2, x_O) = calculate_A2_2A_decomposition(T, gas_subs.clone());

            let kp_func = set_equilibrium_system_dissociation(gas_subs.clone(), T);
            let kp_pa = kp_func(T) * pressure_pa; // Convert K to Kp in Pa
            let (x_o_direct, x_o2_direct) = dissociation_a2_to_2a_from_kp_pa(kp_pa, pressure_pa);
            println!(
                "   \n T {} K, from Kp => x_O: {:.5}, x_O2: {:.5}  (from Kp) ",
                T, x_o_direct, x_o2_direct
            );
            println!(
                "   \n T {} K, from system => x_O: {:.5}, x_O2: {:.5}  (from Kp) ",
                T, x_O, x_O2
            );
            assert_relative_eq!(x_O, x_o_direct, epsilon = 1e-3);
        }
    }

    #[test]
    fn test_calculate_Cl2_2Cl_decomposition() {
        let gas_subs = vec!["CL2".to_string(), "CL".to_string()];
        let temperatures = vec![3000.0, 4000.0, 5000.0];
        let pressure_pa = 101325.0; // 1 atm in Pa
        for &T in &temperatures {
            println!("\n Calculating at T = {} K", T);
            let (x_O2, x_O) = calculate_A2_2A_decomposition(T, gas_subs.clone());

            let kp_func = set_equilibrium_system_dissociation(gas_subs.clone(), T);
            let kp_pa = kp_func(T) * pressure_pa; // Convert K to Kp in Pa
            let (x_o_direct, x_o2_direct) = dissociation_a2_to_2a_from_kp_pa(kp_pa, pressure_pa);
            println!(
                "   \n T {} K, from Kp => x_O: {:.5}, x_O2: {:.5}  (from Kp) ",
                T, x_o_direct, x_o2_direct
            );
            println!(
                "   \n T {} K, from system => x_O: {:.5}, x_O2: {:.5}  (from Kp) ",
                T, x_O, x_O2
            );
            assert_relative_eq!(x_O, x_o_direct, epsilon = 1e-3);
        }
    }

    //////////////////////////SIMPLE GAS MIXTURE ( a A = b B + c C) TESTS ///////////////////////

    /// Create equilibrium constant function K(T) for reaction: a A = b B + c C
    ///
    /// `subs` expected as [A, B, C]
    ///
    ///
    pub fn Const_equilibrium_system_abc(
        subs: Vec<String>,
        _T: f64,
        a: f64,
        b: f64,
        c: f64,
    ) -> Box<dyn Fn(f64) -> f64> {
        assert!(
            subs.len() == 3,
            "This test is designed for three substances (A, B, C)."
        );

        let mut subdata = SubsData::new();

        // Set library priorities
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);

        // Set phases (assume gas for all three)
        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[2].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();

        let _ = subdata.search_substances();
        let _ = subdata.extract_thermal_coeffs(subs[0].as_str(), _T);
        let _ = subdata.extract_thermal_coeffs(subs[1].as_str(), _T);
        let _ = subdata.extract_thermal_coeffs(subs[2].as_str(), _T);
        let _ = subdata.calculate_therm_map_of_sym();

        // Get symbolic dH and dS for each species
        let A_sym = subdata.therm_map_of_sym.get(&subs[0].clone()).unwrap();
        let S_A = *A_sym
            .get(&DataType::dS_sym)
            .clone()
            .unwrap()
            .clone()
            .unwrap();
        let H_A = *A_sym
            .get(&DataType::dH_sym)
            .clone()
            .unwrap()
            .clone()
            .unwrap();

        let B_sym = subdata.therm_map_of_sym.get(&subs[1].clone()).unwrap();
        let S_B = *B_sym
            .get(&DataType::dS_sym)
            .clone()
            .unwrap()
            .clone()
            .unwrap();
        let H_B = *B_sym
            .get(&DataType::dH_sym)
            .clone()
            .unwrap()
            .clone()
            .unwrap();

        let C_sym = subdata.therm_map_of_sym.get(&subs[2].clone()).unwrap();
        let S_C = *C_sym
            .get(&DataType::dS_sym)
            .clone()
            .unwrap()
            .clone()
            .unwrap();
        let H_C = *C_sym
            .get(&DataType::dH_sym)
            .clone()
            .unwrap()
            .clone()
            .unwrap();

        let T_var = Expr::Var("T".to_string());
        let dG_A = H_A - T_var.clone() * S_A;
        let dG_B = H_B - T_var.clone() * S_B;
        let dG_C = H_C - T_var * S_C;

        // ΔG_mix = b*G_B + c*G_C - a*G_A
        let dG_mix = Expr::Const(b) * dG_B + Expr::Const(c) * dG_C - Expr::Const(a) * dG_A;
        let dG_mix_simplified = dG_mix.simplify();
        let dG_fun = dG_mix_simplified.lambdify1D();

        // Return K(T) = exp( -ΔG_mix / (R T) ), R in J/mol/K
        let R = 8.31446261815324_f64;
        Box::new(move |T: f64| f64::exp(-dG_fun(T) / (R * T)))
    }
    /// wrapper function for calculation aA<=>bB+cC equilibrium using crate's native methods
    fn calculate_aA_bB_cC_decomposition(T: f64, gas_subs: Vec<String>) -> (f64, f64, f64) {
        let P = 101325.0;

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
            CustomSubstance::OnePhase(onephase) => {
                let subs = &onephase.subs_data.substances;
                assert_eq!(subs.len(), 3);
            }
            _ => panic!(),
        }

        let mut n = HashMap::new();
        let map_of_gas = HashMap::from([(gas_subs.clone()[0].clone(), 1.0)]);

        n.insert(Some("gas".to_string()), (Some(1.0), Some(map_of_gas)));

        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(Some(T), P, None, None);
        assert!(td.is_ok());

        let mut td = td.unwrap();
        td.Tm = 1e3;
        td.set_number_of_moles(Some(n)).unwrap();
        td.initial_composition().unwrap();
        td.create_indexed_variables();
        td.create_full_system_of_equations().unwrap();

        use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamicsSolver::{
            NRParams, SolverType,
        };
        let params = NRParams {
            initial_guess: None,
            tolerance: 1e-6,
            max_iterations: 10000,
            damping_factor: Some(1.0),
            method: Some(Method::damped),
            loglevel: Some("info".to_string()),
            ..Default::default()
        };
        td.set_solver_type(SolverType::NewtonRaphson(params));
        td.generate_eqs();
        let sol = td.solve_for_T(T);
        //  let sol = td.solver.get_solution();
        println!("{:#?}", sol);
        let A = gas_subs[0].clone();
        let B = gas_subs[1].clone();
        let C = gas_subs[2].clone();
        let n_A = *sol.get(&A).unwrap();
        let n_B = *sol.get(&B).unwrap();
        let n_C = *sol.get(&C).unwrap();
        let N = n_A + n_B + n_C;
        let x_A = n_A / N;
        let x_B = n_B / N;
        let x_C = n_C / N;

        (x_A, x_B, x_C)
    }
    //////////////////////2NO<=>N2 + O2/////////////////////
    /*
    nNO,0​=n0​
    nN2,0​= nO2,0​=0
    Consider the reaction: 2 NO  <=>  N2 + O2
    nN2=ξ
    nO2=ξ
    nNO2=n0​−2ξ

    Kp​=​pN2​​pO2​​ /pNO^2​=(x(N2)​​P)(x(O2)​​P)/(x(NO)​P)^2​=ξ*ξ/(n0​−2ξ)^2
    Kp*n0^2 -4Kp*n0*ξ + 4Kpξ^2 = ξ^2
    (1−4Kp​)ξ^2 = Kp*n0^2 -4Kp*n0*ξ
    Rearranging gives a quadratic in ξ:
    (1−4Kp​)ξ^2+(4Kp​n0​)ξ+(−Kp​n02​)=0.
    set n0​=1 for simplicity
    (1−4Kp​)ξ^2+(4Kp​)ξ+(−Kp​)=0.
    ξ = K^0.5 / (1 + 2 K^0.5)
       */
    fn NO_dissociation_from_kp_pa(kp_pa: f64) -> (f64, f64, f64) {
        // Solve quadratic for ξ
        let eta = (kp_pa).sqrt() / (1.0 + 2.0 * (kp_pa).sqrt());

        let n_NO = 1.0 - 2.0 * eta;
        let n_N2 = eta;
        let n_O2 = eta;
        let N = n_NO + n_N2 + n_O2;
        let x_NO = n_NO / N;
        let x_N2 = n_N2 / N;
        let x_O2 = n_O2 / N;
        (x_NO, x_N2, x_O2)
    }
    #[test]
    fn test_Keq_NO_decomposition() {
        let subs = vec!["NO".to_string(), "N2".to_string(), "O2".to_string()];
        let T = 4500.0;
        let K = Const_equilibrium_system_abc(subs, T, 2.0, 1.0, 1.0);
        let kp_pa = K(T);

        let (x_NO, x_N2, x_O2) = NO_dissociation_from_kp_pa(kp_pa);
        println!(
            "T: {} K, Kp: {:.3e} Pa => x_NO: {:.5}, x_N2: {:.5}, x_O2: {:.5}",
            T, kp_pa, x_NO, x_N2, x_O2
        );
    }
    #[test]
    fn test_simple_gas_mixture_no_n2_o2_system1() {
        let T = 4500.0;
        let gas_subs = vec!["NO".to_string(), "N2".to_string(), "O2".to_string()];
        let (x_NO, x_N2, x_O2) = calculate_aA_bB_cC_decomposition(T, gas_subs.clone());
        println!(
            "T: {} K, x_NO: {:.5}, x_N2: {:.5}, x_O2: {:.5}",
            T, x_NO, x_N2, x_O2
        );
    }

    #[test]
    fn test_NO_compare_methods() {
        let subs = vec!["NO".to_string(), "N2".to_string(), "O2".to_string()];
        let temperatures = vec![3000.0, 4000.0, 5000.0];

        for &T in &temperatures {
            let K = Const_equilibrium_system_abc(subs.clone(), T, 2.0, 1.0, 1.0);
            let kp_pa = K(T);

            let (x_NO_direct, x_N2_direct, x_O2_direct) = NO_dissociation_from_kp_pa(kp_pa);
            let (x_NO, x_N2, x_O2) = calculate_aA_bB_cC_decomposition(T, subs.clone());
            assert_relative_eq!(x_NO_direct, x_NO, epsilon = 1e-2);
            assert_relative_eq!(x_N2_direct, x_N2, epsilon = 1e-2);
            assert_relative_eq!(x_O2_direct, x_O2, epsilon = 1e-2);
        }
    }
    //////////////////////////////debug//////////////////////////////////////////////////
    #[test]
    fn debug_NO_equilibrium_methods_comparison() {
        let T = 4500.0;
        let gas_subs = vec!["NO".to_string(), "N2".to_string(), "O2".to_string()];

        println!(
            "\n=== DEBUGGING NO EQUILIBRIUM COMPARISON AT T = {} K ===",
            T
        );

        // Method 1: Equilibrium constant approach
        println!("\n--- Method 1: Equilibrium Constant ---");
        let K_func = Const_equilibrium_system_abc(gas_subs.clone(), T, 2.0, 1.0, 1.0);
        let kp = K_func(T);
        println!("K(T): {:.6e}", kp);

        let (x_NO_kp, x_N2_kp, x_O2_kp) = NO_dissociation_from_kp_pa(kp);
        println!(
            "From Kp: x_NO={:.6}, x_N2={:.6}, x_O2={:.6}",
            x_NO_kp, x_N2_kp, x_O2_kp
        );

        // Method 2: System approach
        println!("\n--- Method 2: Thermodynamic System ---");
        let (x_NO_sys, x_N2_sys, x_O2_sys) = calculate_aA_bB_cC_decomposition(T, gas_subs.clone());
        println!(
            "From System: x_NO={:.6}, x_N2={:.6}, x_O2={:.6}",
            x_NO_sys, x_N2_sys, x_O2_sys
        );

        // Compare
        println!("\n--- Differences ---");
        println!("Δx_NO = {:.6e}", (x_NO_sys - x_NO_kp).abs());
        println!("Δx_N2 = {:.6e}", (x_N2_sys - x_N2_kp).abs());
        println!("Δx_O2 = {:.6e}", (x_O2_sys - x_O2_kp).abs());

        if (x_NO_sys - x_NO_kp).abs() > 1e-3 {
            println!("\n*** SIGNIFICANT MISMATCH DETECTED ***");
            debug_equilibrium_constant_calculation(&gas_subs, T);
        }
    }

    fn debug_equilibrium_constant_calculation(subs: &[String], T: f64) {
        println!("\n--- Debugging Equilibrium Constant Calculation ---");

        let mut subdata = SubsData::new();
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);

        for sub in subs {
            subdata.map_of_phases.insert(sub.clone(), Some(Phases::Gas));
        }
        subdata.substances = subs.to_vec();

        let _ = subdata.search_substances();
        for sub in subs {
            let _ = subdata.extract_thermal_coeffs(sub.as_str(), T);
        }
        let _ = subdata.calculate_therm_map_of_sym();

        // Calculate individual Gibbs energies
        let R = 8.31446261815324_f64;
        for (i, sub) in subs.iter().enumerate() {
            if let Some(therm_data) = subdata.therm_map_of_sym.get(sub) {
                if let Some(Some(h_sym)) = therm_data.get(&DataType::dH_sym) {
                    if let Some(Some(s_sym)) = therm_data.get(&DataType::dS_sym) {
                        let h_func = h_sym.lambdify1D();
                        let s_func = s_sym.lambdify1D();
                        let h_val = h_func(T);
                        let s_val = s_func(T);
                        let g_val = h_val - T * s_val;
                        println!(
                            "{}: H={:.3e}, S={:.3e}, G={:.3e} J/mol",
                            sub, h_val, s_val, g_val
                        );
                    }
                }
            }
        }

        // Calculate reaction Gibbs energy: 2NO -> N2 + O2
        // ΔG = G_N2 + G_O2 - 2*G_NO
        let get_gibbs = |sub: &str| -> f64 {
            if let Some(therm_data) = subdata.therm_map_of_sym.get(sub) {
                if let Some(Some(h_sym)) = therm_data.get(&DataType::dH_sym) {
                    if let Some(Some(s_sym)) = therm_data.get(&DataType::dS_sym) {
                        let h_func = h_sym.lambdify1D();
                        let s_func = s_sym.lambdify1D();
                        return h_func(T) - T * s_func(T);
                    }
                }
            }
            0.0
        };

        let g_NO = get_gibbs("NO");
        let g_N2 = get_gibbs("N2");
        let g_O2 = get_gibbs("O2");

        let delta_G = g_N2 + g_O2 - 2.0 * g_NO;
        let K_calculated = f64::exp(-delta_G / (R * T));

        println!("\nReaction: 2NO -> N2 + O2");
        println!("ΔG_reaction = {:.3e} J/mol", delta_G);
        println!("K_calculated = exp(-ΔG/RT) = {:.6e}", K_calculated);
    }

    #[test]
    fn verify_NO_dissociation_formula() {
        println!("\n=== Verifying NO Dissociation Formula ===");

        // Test with known values
        let test_kp = 1e-6;
        let (x_NO, x_N2, x_O2) = NO_dissociation_from_kp_pa(test_kp);

        println!("Test Kp = {:.1e}", test_kp);
        println!(
            "Results: x_NO={:.6}, x_N2={:.6}, x_O2={:.6}",
            x_NO, x_N2, x_O2
        );
        println!("Sum = {:.8}", x_NO + x_N2 + x_O2);

        // Verify equilibrium condition: Kp = (x_N2 * x_O2) / x_NO^2
        let kp_check = (x_N2 * x_O2) / (x_NO * x_NO);
        println!("Kp check = {:.6e} (should equal {:.6e})", kp_check, test_kp);
        println!(
            "Relative error = {:.2e}",
            (kp_check - test_kp).abs() / test_kp
        );

        assert!(
            (kp_check - test_kp).abs() / test_kp < 1e-10,
            "Equilibrium formula verification failed"
        );
    }

    #[test]
    fn debug_equilibrium_constant_units_and_calculation() {
        println!("\n=== Debugging Equilibrium Constant Calculation ===");

        let T = 4500.0;
        let gas_subs = vec!["NO".to_string(), "N2".to_string(), "O2".to_string()];

        // Method 1: Manual calculation
        println!("\n--- Manual Calculation ---");
        let mut subdata = SubsData::new();
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);

        for sub in &gas_subs {
            subdata.map_of_phases.insert(sub.clone(), Some(Phases::Gas));
        }
        subdata.substances = gas_subs.clone();

        let _ = subdata.search_substances();
        for sub in &gas_subs {
            let _ = subdata.extract_thermal_coeffs(sub.as_str(), T);
        }
        let _ = subdata.calculate_therm_map_of_sym();

        let get_gibbs = |sub: &str| -> f64 {
            if let Some(therm_data) = subdata.therm_map_of_sym.get(sub) {
                if let Some(Some(h_sym)) = therm_data.get(&DataType::dH_sym) {
                    if let Some(Some(s_sym)) = therm_data.get(&DataType::dS_sym) {
                        let h_func = h_sym.lambdify1D();
                        let s_func = s_sym.lambdify1D();
                        return h_func(T) - T * s_func(T);
                    }
                }
            }
            0.0
        };

        let g_NO = get_gibbs("NO");
        let g_N2 = get_gibbs("N2");
        let g_O2 = get_gibbs("O2");

        println!("G_NO = {:.3e} J/mol", g_NO);
        println!("G_N2 = {:.3e} J/mol", g_N2);
        println!("G_O2 = {:.3e} J/mol", g_O2);

        // For reaction: 2NO -> N2 + O2
        let delta_G = g_N2 + g_O2 - 2.0 * g_NO;
        println!("ΔG_reaction = {:.3e} J/mol", delta_G);

        let R = 8.31446261815324_f64;
        let K_manual = f64::exp(-delta_G / (R * T));
        println!("K_manual = {:.6e}", K_manual);

        // Method 2: Using the function
        println!("\n--- Function Calculation ---");
        let K_func = Const_equilibrium_system_abc(gas_subs.clone(), T, 2.0, 1.0, 1.0);
        let K_function = K_func(T);
        println!("K_function = {:.6e}", K_function);

        println!("\n--- Comparison ---");
        println!("Difference = {:.6e}", (K_manual - K_function).abs());
        println!(
            "Relative error = {:.6e}",
            (K_manual - K_function).abs() / K_manual
        );

        // Test both with the dissociation formula
        println!("\n--- Testing with Dissociation Formula ---");
        let (x_NO_manual, x_N2_manual, x_O2_manual) = NO_dissociation_from_kp_pa(K_manual);
        let (x_NO_func, x_N2_func, x_O2_func) = NO_dissociation_from_kp_pa(K_function);

        println!(
            "Manual K: x_NO={:.6}, x_N2={:.6}, x_O2={:.6}",
            x_NO_manual, x_N2_manual, x_O2_manual
        );
        println!(
            "Function K: x_NO={:.6}, x_N2={:.6}, x_O2={:.6}",
            x_NO_func, x_N2_func, x_O2_func
        );
    }

    #[test]
    fn debug_NO2_equilibrium_comparison() {
        let T = 500.0;
        let gas_subs = vec!["NO2".to_string(), "N2".to_string(), "O2".to_string()];

        println!("\n=== DEBUGGING NO2 EQUILIBRIUM AT T = {} K ===", T);

        // Method 1: Equilibrium constant approach
        println!("\n--- Method 1: Equilibrium Constant (2NO2 -> N2 + 2O2) ---");
        let K_func = Const_equilibrium_system_abc(gas_subs.clone(), T, 2.0, 1.0, 2.0);
        let kp = K_func(T);
        println!("K(T): {:.6e}", kp);

        let P = 101325.0;
        let (x_NO2_kp, x_N2_kp, x_O2_kp) = NO2_dissociation_from_kp_pa(kp, P);
        println!(
            "From Kp: x_NO2={:.6}, x_N2={:.6}, x_O2={:.6}",
            x_NO2_kp, x_N2_kp, x_O2_kp
        );

        // Method 2: System approach
        println!("\n--- Method 2: Thermodynamic System ---");
        let (x_NO2_sys, x_N2_sys, x_O2_sys) = calculate_aA_bB_cC_decomposition(T, gas_subs.clone());
        println!(
            "From System: x_NO2={:.6}, x_N2={:.6}, x_O2={:.6}",
            x_NO2_sys, x_N2_sys, x_O2_sys
        );

        // Compare
        println!("\n--- Differences ---");
        println!("Δx_NO2 = {:.6e}", (x_NO2_sys - x_NO2_kp).abs());
        println!("Δx_N2 = {:.6e}", (x_N2_sys - x_N2_kp).abs());
        println!("Δx_O2 = {:.6e}", (x_O2_sys - x_O2_kp).abs());

        // Debug the NO2 dissociation formula
        println!("\n--- Debugging NO2 Formula ---");
        println!("Initial assumption: 1 mol NO2, 0 mol N2, 0 mol O2");
        println!("Reaction: 2NO2 -> N2 + 2O2");
        println!("If ξ is extent of reaction:");
        println!("n_NO2 = 1 - 2ξ, n_N2 = ξ, n_O2 = 2ξ");
        println!("Total moles = 1 - 2ξ + ξ + 2ξ = 1 + ξ");
        println!(
            "Kp =   (x_NO2 * P)^2/(x_N2 * P) * (x_O2 * P)^2 = n_total* (n_NO2 * P)^2/(n_N2 * P) * (n_O2 * P)^2"
        );
        println!("Kp*P =  (1 - 2ξ )^2/(ξ ) * (2ξ )^2*(1 + ξ)");
        println!("Kp*P*4*ξ ^3  =  (1 + ξ)*(1 - 2ξ )^2");

        if (x_NO2_sys - x_NO2_kp).abs() > 1e-3 {
            println!("\n*** SIGNIFICANT MISMATCH DETECTED ***");
            debug_NO2_formula_verification(kp, P);
        }
    }

    fn debug_NO2_formula_verification(kp: f64, P: f64) {
        println!("\n--- Verifying NO2 Formula Implementation ---");

        // Check if the current formula is correct
        let (x_NO2, x_N2, x_O2) = NO2_dissociation_from_kp_pa(kp, P);

        // Back-calculate Kp from the results
        let kp_check = (x_N2 * P) * (x_O2 * P).powi(2) / (x_NO2 * P).powi(2);
        println!("Original Kp: {:.6e}", kp);
        println!("Back-calculated Kp: {:.6e}", kp_check);
        println!("Relative error: {:.6e}", (kp_check - kp).abs() / kp);

        if (kp_check - kp).abs() / kp > 1e-6 {
            println!("*** ERROR: NO2 dissociation formula is incorrect! ***");
        } else {
            println!("NO2 dissociation formula is mathematically correct.");
        }
    }

    ////////////////////2N2O<=>2N2 +2O///////////////////////////

    /*
    nN2O,0​=n0=1.0​
    nN2,0​= nO2,0​=0
    4(Kp​−P)ξ^3−3Kp​ξ+Kp​=0.

         */
    fn N2O_dissociation_from_kp_pa(kp: f64, P: f64) -> (f64, f64, f64) {
        use RustedSciThe::numerical::optimization::minimize_scalar::{
            ClosureFunction, ScalarRootFinder,
        };
        let solver = ScalarRootFinder::new();

        let func = ClosureFunction::new(
            |x| 4.0 * (1.0 - P / kp) * x.powi(3) - 3.0 * x + 1.0,
            "".to_string(),
        );

        let result = solver.secant(&func, 0.2, 0.7).unwrap();
        let eta = result.root;
        let n_N2O = 1.0 - 2.0 * eta;
        let n_N2 = 2.0 * eta;
        let n_O2 = eta;
        let N = n_N2O + n_N2 + n_O2;
        let x_N2O = n_N2O / N;
        let x_N2 = n_N2 / N;
        let x_O2 = n_O2 / N;
        (x_N2O, x_N2, x_O2)
    }
    #[test]
    fn test_Keq_N2O_decomposition() {
        let subs = vec!["N2O".to_string(), "N2".to_string(), "O2".to_string()];
        let T = 4500.0;
        let K = Const_equilibrium_system_abc(subs, T, 2.0, 2.0, 1.0);
        let kp_pa = K(T);
        let P = 101325.0;
        let (x_N2O, x_N2, x_O2) = N2O_dissociation_from_kp_pa(kp_pa, P);
        println!(
            "T: {} K, Kp: {:.3e} Pa => x_N2O: {:.5}, x_N2: {:.5}, x_O2: {:.5}",
            T, kp_pa, x_N2O, x_N2, x_O2
        );
    }
    #[test]
    fn test_simple_gas_mixture_n2o_n2_o2_system1() {
        let T = 4500.0;
        let gas_subs = vec!["N2O".to_string(), "N2".to_string(), "O2".to_string()];
        let (x_N2O, x_N2, x_O2) = calculate_aA_bB_cC_decomposition(T, gas_subs.clone());
        println!(
            "T: {} K, x_N2O: {:.5}, x_N2: {:.5}, x_O2: {:.5}",
            T, x_N2O, x_N2, x_O2
        );
    }

    #[test]
    fn test_N2O_compare_methods() {
        let subs = vec!["N2O".to_string(), "N2".to_string(), "O2".to_string()];
        let temperatures = vec![3000.0, 4000.0, 5000.0];
        let P = 101325.0;
        for &T in &temperatures {
            let K = Const_equilibrium_system_abc(subs.clone(), T, 2.0, 2.0, 1.0);
            let kp_pa = K(T);

            let (x_N2O_direct, x_N2_direct, x_O2_direct) = N2O_dissociation_from_kp_pa(kp_pa, P);
            let (x_N2O, x_N2, x_O2) = calculate_aA_bB_cC_decomposition(T, subs.clone());
            assert_relative_eq!(x_N2O_direct, x_N2O, epsilon = 1e-2);
            assert_relative_eq!(x_N2_direct, x_N2, epsilon = 1e-2);
            assert_relative_eq!(x_O2_direct, x_O2, epsilon = 1e-2);
        }
    }
    ////////////////2NO2 <=> N2 + 2O2////////////////////////
    /*
    nNO2,0​=n0​
    nN2,0​= nO2,0​=0
    (n0 + 2ξ)^2(n0 - ξ)+ 4(KpP)ξ^3 = 0
    */
    fn NO2_dissociation_from_kp_pa(kp: f64, P: f64) -> (f64, f64, f64) {
        use RustedSciThe::numerical::optimization::minimize_scalar::{
            ClosureFunction, ScalarRootFinder,
        };
        let solver = ScalarRootFinder::new();

        let func = ClosureFunction::new(
            |x| (1.0 + 2.0 * x).powi(2) * (1.0 - x) + 4.0 * P * kp * x.powi(3),
            "".to_string(),
        );

        let result = solver.secant(&func, 0.2, 0.7).unwrap();
        let eta = result.root;
        let n_NO2 = 1.0 + 2.0 * eta;
        let n_N2 = -eta;
        let n_O2 = -2.0 * eta;
        let N = n_NO2 + n_N2 + n_O2;
        let x_NO2 = n_NO2 / N;
        let x_N2 = n_N2 / N;
        let x_O2 = n_O2 / N;
        (x_NO2, x_N2, x_O2)
    }
    #[test]
    fn test_Keq_NO2_decomposition() {
        let subs = vec!["NO2".to_string(), "N2".to_string(), "O2".to_string()];
        let T = 1500.0;
        let K = Const_equilibrium_system_abc(subs, T, 2.0, 1.0, 2.0);
        let kp_pa = K(T);
        let P = 101325.0;
        let (x_NO2, x_N2, x_O2) = NO2_dissociation_from_kp_pa(kp_pa, P);
        println!(
            "T: {} K, Kp: {:.3e} Pa => x_NO2: {:.5}, x_N2: {:.5}, x_O2: {:.5}",
            T, kp_pa, x_NO2, x_N2, x_O2
        );
    }

    #[test]
    fn test_simple_gas_mixture_no2_n2_o2_system1() {
        let gas_subs = vec!["NO2".to_string(), "N2".to_string(), "O2".to_string()];

        let T = 500.0;

        let (x_A, x_B, x_C) = calculate_aA_bB_cC_decomposition(T, gas_subs.clone());
        println!(
            "\n T: {} K,  => x_NO2: {:.5}, x_N2: {:.5}, x_O2: {:.5}",
            T, x_A, x_B, x_C
        );
    }

    #[test]
    fn test_simple_gas_mixture_no2_n2_o2_system2() {
        let gas_subs = vec!["NO2".to_string(), "N2".to_string(), "O2".to_string()];

        let temperatures = vec![500.0];

        for &T in &temperatures {
            let (x_NO2, x_N2, x_O2) = calculate_aA_bB_cC_decomposition(T, gas_subs.clone());

            assert!(
                (x_NO2 + x_N2 + x_O2 - 1.0).abs() < 1e-8,
                "Mole fractions do not sum to 1"
            );

            let K = Const_equilibrium_system_abc(gas_subs.clone(), T, 2.0, 1.0, 2.0);
            let kp_pa = K(T);
            let P = 101325.0;
            let (x_NO2_direct, x_N2_direct, x_O2_direct) = NO2_dissociation_from_kp_pa(kp_pa, P);

            println!("\nT: {} K, Kp: {:.3e}", T, kp_pa);
            println!(
                "From Kp: x_NO2={:.6}, x_N2={:.6}, x_O2={:.6}",
                x_NO2_direct, x_N2_direct, x_O2_direct
            );
            println!(
                "From System: x_NO2={:.6}, x_N2={:.6}, x_O2={:.6}",
                x_NO2, x_N2, x_O2
            );
            println!(
                "Differences: Δx_NO2={:.6e}, Δx_N2={:.6e}, Δx_O2={:.6e}",
                (x_NO2 - x_NO2_direct).abs(),
                (x_N2 - x_N2_direct).abs(),
                (x_O2 - x_O2_direct).abs()
            );

            assert_relative_eq!(x_NO2_direct, x_NO2, epsilon = 1e-2);
            assert_relative_eq!(x_N2_direct, x_N2, epsilon = 1e-2);
            assert_relative_eq!(x_O2_direct, x_O2, epsilon = 1e-2);
        }
    }

    /////////////////////////////COMPLICATED GAS SYSTEMS//////////////////////
    ///

    fn create_N_plus_O_equilibrium(T: f64) -> Thermodynamics {
        let P = 101325.0;
        let gas_subs = vec![
            "N2".to_string(),
            "O2".to_string(),
            "NO2".to_string(),
            "N2O".to_string(),
            "NO".to_string(),
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
            CustomSubstance::OnePhase(onephase) => {
                let subs = &onephase.subs_data.substances;
                assert_eq!(subs.len(), 5);
            }
            _ => panic!(),
        }

        let mut n = HashMap::new();
        let map_of_gas = HashMap::from([("N2".to_string(), 1.0), ("O2".to_string(), 1.0)]);

        n.insert(Some("gas".to_string()), (Some(1.0), Some(map_of_gas)));

        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(Some(T), P, None, None);
        assert!(td.is_ok());

        let mut td = td.unwrap();
        td.Tm = 1e3;
        td.set_number_of_moles(Some(n)).unwrap();
        td.initial_composition().unwrap();
        td.create_indexed_variables();
        td.create_full_system_of_equations().unwrap();
        let params = NRParams {
            initial_guess: None,
            tolerance: 1e-6,
            max_iterations: 10_000,
            damping_factor: Some(1.0),
            method: Some(Method::damped),
            loglevel: None, // Some("info".to_string()),
        };
        td.set_solver_type(SolverType::NewtonRaphson(params));
        td.generate_eqs();

        td
        //  HashMap::new()
    }

    fn create_N_plus_O_equilibrium_LM(T: f64) -> Thermodynamics {
        let P = 101325.0;
        let gas_subs = vec![
            "N2".to_string(),
            "O2".to_string(),
            "NO2".to_string(),
            "N2O".to_string(),
            "NO".to_string(),
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
            CustomSubstance::OnePhase(onephase) => {
                let subs = &onephase.subs_data.substances;
                assert_eq!(subs.len(), 5);
            }
            _ => panic!(),
        }

        let mut n = HashMap::new();
        let map_of_gas = HashMap::from([("N2".to_string(), 1.0), ("O2".to_string(), 1.0)]);

        n.insert(Some("gas".to_string()), (Some(1.0), Some(map_of_gas)));

        // create thermodynamics instance
        let td = customsubs.create_thermodynamics(Some(T), P, None, None);
        assert!(td.is_ok());

        let mut td = td.unwrap();
        td.Tm = 1e5;
        td.set_number_of_moles(Some(n)).unwrap();
        td.initial_composition().unwrap();
        td.create_indexed_variables();
        td.create_full_system_of_equations().unwrap();
        let params = LMParams::default();
        td.set_solver_type(SolverType::LevenbergMarquardt(params));
        td.generate_eqs();

        td
        //  HashMap::new()
    }

    pub fn solve_N_plus_O_equilibrium(T: f64, td: &mut Thermodynamics) -> HashMap<String, f64> {
        let map_of_substances = td.solve_for_T(T);
        // let map_of_substances = td.solve_NR(None, 1e-6, 10000, Some(0.01), None, None);
        let sol = td.solver.get_solution();

        println!("{:#?}", sol);
        println!("{:#?}", map_of_substances);
        map_of_substances
    }
    #[test]
    fn N_plus_O_equilibrium_500K() {
        let mut td = create_N_plus_O_equilibrium(500.0);
        let map_of_substances = solve_N_plus_O_equilibrium(500.0, &mut td);

        let n_O2 = *map_of_substances.get("O2").unwrap();
        let n_N2 = *map_of_substances.get("N2").unwrap();
        assert_relative_eq!(n_N2, 1.0, epsilon = 1e-4);
        assert_relative_eq!(n_O2, 1.0, epsilon = 1e-4);
    }

    #[test]
    fn N_plus_O_equilibrium_2500K_LM() {
        let mut td = create_N_plus_O_equilibrium_LM(2500.0);
        let map_of_substances = solve_N_plus_O_equilibrium(2500.0, &mut td);

        let _n_O2 = *map_of_substances.get("O2").unwrap();
        let _n_N2 = *map_of_substances.get("N2").unwrap();
    }
    #[test]
    fn N_plus_O_equilibrium_4500K() {
        let mut td = create_N_plus_O_equilibrium(2500.0);
        let map_of_substances = solve_N_plus_O_equilibrium(2500.0, &mut td);
        let _n_O2 = *map_of_substances.get("O2").unwrap();
        let _n_N2 = *map_of_substances.get("N2").unwrap();
    }
    //cargo test tests::N_plus_O_equilibrium_load_test --release -- --nocapture
    #[test]
    fn N_plus_O_equilibrium_load_test() {
        let start = 500.0_f64;
        let end = 2500.0_f64;
        let step = 10.0_f64;

        let mut T = start;
        let now = Instant::now();
        let mut td = create_N_plus_O_equilibrium(500.0);

        while T <= end {
            solve_N_plus_O_equilibrium(T, &mut td);
            T += step;
        }
        let elapsed = now.elapsed().as_secs();
        println!("Elapsed: {:.2} seconds", elapsed);
    }

    ///////////////////////// EASY EQUILIBRIUM TESTS /////////////////////////

    #[test]
    fn test_easy_equilibrium_o2_dissociation() {
        let subs = vec!["O2".to_string(), "O".to_string()];
        let mut subdata = SubsData::new();
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);
        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();
        let _ = subdata.search_substances();
        let _ = subdata.extract_thermal_coeffs(subs[0].as_str(), 3000.0);
        let _ = subdata.extract_thermal_coeffs(subs[1].as_str(), 3000.0);
        let _ = subdata.calculate_therm_map_of_sym();

        let mut subs_coeffs = HashMap::new();
        subs_coeffs.insert("O2".to_string(), -1.0);
        subs_coeffs.insert("O".to_string(), 2.0);

        let mut initial_moles = HashMap::new();
        initial_moles.insert("O2".to_string(), 1.0);
        initial_moles.insert("O".to_string(), 0.0);

        let mut eq = EasyEquilibrium::new(101325.0, subdata, subs_coeffs, initial_moles);
        let K = eq.create_equilibrium_const();
        let T = 3000.0;
        let kp = K(T);
        println!(
            "\n \n eq {:?} \n \n LHS: {}",
            &eq.clone(),
            &eq.eta_equauion()
        );
        let kp_func = set_equilibrium_system_dissociation(subs.clone(), T);
        let kp_expected = kp_func(T);

        assert_relative_eq!(kp, kp_expected, epsilon = 1e-6);
    }

    #[test]
    fn test_easy_equilibrium_no_decomposition() {
        let subs = vec!["NO".to_string(), "N2".to_string(), "O2".to_string()];
        let mut subdata = SubsData::new();
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);
        for sub in &subs {
            subdata.map_of_phases.insert(sub.clone(), Some(Phases::Gas));
        }
        subdata.substances = subs.clone();
        let _ = subdata.search_substances();
        for sub in &subs {
            let _ = subdata.extract_thermal_coeffs(sub.as_str(), 4500.0);
        }
        let _ = subdata.calculate_therm_map_of_sym();

        let mut subs_coeffs = HashMap::new();
        subs_coeffs.insert("NO".to_string(), -2.0);
        subs_coeffs.insert("N2".to_string(), 1.0);
        subs_coeffs.insert("O2".to_string(), 1.0);

        let mut initial_moles = HashMap::new();
        initial_moles.insert("NO".to_string(), 1.0);
        initial_moles.insert("N2".to_string(), 0.0);
        initial_moles.insert("O2".to_string(), 0.0);

        let mut eq = EasyEquilibrium::new(101325.0, subdata, subs_coeffs, initial_moles);
        let K = eq.create_equilibrium_const();
        let T = 4500.0;
        let kp = K(T);
        println!(
            "\n \n eq {:?} \n \n LHS: {}",
            &eq.clone(),
            &eq.eta_equauion()
        );
        let K_func = Const_equilibrium_system_abc(subs.clone(), T, 2.0, 1.0, 1.0);
        let kp_expected = K_func(T);

        assert_relative_eq!(kp, kp_expected, epsilon = 1e-6);
    }

    #[test]
    fn test_easy_equilibrium_solve() {
        let subs = vec!["O2".to_string(), "O".to_string()];
        let mut subdata = SubsData::new();
        subdata
            .library_priorities
            .insert("NASA_gas".to_string(), LibraryPriority::Priority);
        subdata
            .map_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        subdata
            .map_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        subdata.substances = subs.clone();
        let _ = subdata.search_substances();
        let _ = subdata.extract_thermal_coeffs(subs[0].as_str(), 500.0);
        let _ = subdata.extract_thermal_coeffs(subs[1].as_str(), 500.0);
        let _ = subdata.calculate_therm_map_of_sym();

        let mut subs_coeffs = HashMap::new();
        subs_coeffs.insert("O2".to_string(), -1.0);
        subs_coeffs.insert("O".to_string(), 2.0);

        let mut initial_moles = HashMap::new();
        initial_moles.insert("O2".to_string(), 1.0);
        initial_moles.insert("O".to_string(), 0.0);

        let mut eq = EasyEquilibrium::new(101325.0, subdata, subs_coeffs, initial_moles);
        println!(
            "\n \n eq {:?} \n \n LHS: {}",
            &eq.clone(),
            &eq.eta_equauion()
        );
        let res = eq.solve(500.0, 50.0, 2000.0);
        println!("{:#?}", res);
        println!("EasyEquilibrium solve test completed");
    }
}
