#[cfg(test)]
mod debug_tests {
    use super::super::ClassicalThermodynamics_tests3::tests::*;
    use crate::Thermodynamics::User_substances::{DataType, LibraryPriority, Phases, SubsData};
    use RustedSciThe::symbolic::symbolic_engine::Expr;

    #[test]
    fn debug_NO_equilibrium_comparison() {
        let T = 4500.0;
        let P = 101325.0;
        let gas_subs = vec!["NO".to_string(), "N2".to_string(), "O2".to_string()];
        
        println!("=== DEBUGGING NO EQUILIBRIUM AT T = {} K ===", T);
        
        // Method 1: Naive approach using equilibrium constant
        println!("\n--- Method 1: Equilibrium Constant Approach ---");
        let K_func = Const_equilibrium_system_abc(gas_subs.clone(), T, 2.0, 1.0, 1.0);
        let kp = K_func(T);
        println!("Equilibrium constant K: {:.6e}", kp);
        
        let (x_NO_naive, x_N2_naive, x_O2_naive) = NO_dissociation_from_kp_pa(kp);
        println!("Naive results: x_NO = {:.6}, x_N2 = {:.6}, x_O2 = {:.6}", 
                x_NO_naive, x_N2_naive, x_O2_naive);
        println!("Sum check: {:.8}", x_NO_naive + x_N2_naive + x_O2_naive);
        
        // Method 2: Sophisticated thermodynamic system
        println!("\n--- Method 2: Thermodynamic System Approach ---");
        let (x_NO_system, x_N2_system, x_O2_system) = calculate_aA_bB_cC_decomposition(T, gas_subs.clone());
        println!("System results: x_NO = {:.6}, x_N2 = {:.6}, x_O2 = {:.6}", 
                x_NO_system, x_N2_system, x_O2_system);
        println!("Sum check: {:.8}", x_NO_system + x_N2_system + x_O2_system);
        
        // Compare results
        println!("\n--- Comparison ---");
        println!("Δx_NO = {:.6e}", (x_NO_system - x_NO_naive).abs());
        println!("Δx_N2 = {:.6e}", (x_N2_system - x_N2_naive).abs());
        println!("Δx_O2 = {:.6e}", (x_O2_system - x_O2_naive).abs());
        
        // Debug thermodynamic data consistency
        println!("\n--- Thermodynamic Data Debug ---");
        debug_thermo_data_consistency(&gas_subs, T);
    }
    
    fn debug_thermo_data_consistency(subs: &[String], T: f64) {
        let mut subdata = SubsData::new();
        subdata.library_priorities.insert("NASA_gas".to_string(), LibraryPriority::Priority);
        
        for sub in subs {
            subdata.map_of_phases.insert(sub.clone(), Some(Phases::Gas));
        }
        subdata.substances = subs.to_vec();
        
        let _ = subdata.search_substances();
        for sub in subs {
            let _ = subdata.extract_thermal_coeffs(sub.as_str(), T);
        }
        let _ = subdata.calculate_therm_map_of_sym();
        
        for (i, sub) in subs.iter().enumerate() {
            println!("Substance {}: {}", i, sub);
            if let Some(therm_data) = subdata.therm_map_of_sym.get(sub) {
                if let Some(Some(h_sym)) = therm_data.get(&DataType::dH_sym) {
                    if let Some(Some(s_sym)) = therm_data.get(&DataType::dS_sym) {
                        let h_func = h_sym.lambdify1D();
                        let s_func = s_sym.lambdify1D();
                        let h_val = h_func(T);
                        let s_val = s_func(T);
                        let g_val = h_val - T * s_val;
                        println!("  H({}) = {:.3e} J/mol", T, h_val);
                        println!("  S({}) = {:.3e} J/mol/K", T, s_val);
                        println!("  G({}) = {:.3e} J/mol", T, g_val);
                    }
                }
            }
        }
    }
}