// module containing basic formulae of chemical thermodynamics
//use crate::Kinetics::molmass::create_elem_composition_matrix;
use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
use crate::Thermodynamics::User_substances::{
    CalculatorType, DataType, Phases, SearchResult, SubsData, WhatIsFound,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use prettytable::{Table, row};
use std::collections::HashMap;
use std::f64;
use std::fmt;
const R: f64 = 8.314;
#[allow(non_upper_case_globals)]
const R_sym: Expr = Expr::Const(R);
pub struct Thermodynamics {
    pub vec_of_substances: Vec<String>,
    pub vec_of_phases: HashMap<String, Option<Phases>>,
    pub T: f64,
    pub P: f64,

    pub subdata: SubsData,
    pub dG: HashMap<String, f64>,
    pub dG_fun: HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>) -> f64 + 'static>>,
    pub dG_sym: HashMap<String, Expr>,
}

impl Clone for Thermodynamics {
    fn clone(&self) -> Self {
        let mut new_dG = HashMap::new();

        // We can't directly clone the functions, so we'll need to handle this field specially
        // For now, we'll create an empty map structure that matches the original
        for (substance, _) in &self.dG_fun {
            let placeholder: Box<dyn Fn(f64, Option<Vec<f64>>) -> f64 + 'static> =
                Box::new(|_: f64, _: Option<Vec<f64>>| {
                    // your function implementation here
                    0.0
                });
            new_dG.insert(substance.clone(), placeholder);
        }
        let mut cloned_instance = Thermodynamics {
            vec_of_substances: self.vec_of_substances.clone(),
            vec_of_phases: self.vec_of_phases.clone(),
            T: self.T,
            P: self.P,
            subdata: self.subdata.clone(),
            dG: self.dG.clone(),
            dG_fun: new_dG,
            dG_sym: self.dG_sym.clone(),
        };
        // Clone the functions
        let _ = cloned_instance.calculate_Gibbs_fun(self.T);
        cloned_instance
    }
}

impl fmt::Debug for Thermodynamics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Thermodynamics")
            .field("vec_of_substances", &self.vec_of_substances)
            .field("T", &self.T)
            .field("P", &self.P)
            .field("subdata", &self.subdata)
            .field("dmu", &self.dG)
            //.field("dmu_fun", &self.dmu_fun.is_some())
            .field("dmu_sym", &self.dG_sym)
            .finish()
    }
}

impl Thermodynamics {
    pub fn new() -> Self {
        Self {
            vec_of_substances: Vec::new(),
            vec_of_phases: HashMap::new(),
            T: 298.15,
            P: 1e5,
            subdata: SubsData::new(),
            dG: HashMap::new(),
            dG_fun: HashMap::new(),
            dG_sym: HashMap::new(),
        }
    }
    pub fn set_T(&mut self, T: f64) {
        self.T = T;
    }
    pub fn set_P(&mut self, P: f64, _P_unit: Option<String>) {
        self.P = P;
    }
    pub fn calculate_Gibbs_free_energy(&mut self, T: f64, w: Option<Vec<f64>>) {
        self.T = T;
        let sd = &mut self.subdata;
        let _ = sd.extract_all_thermal_coeffs(T);
        let _ = sd.calculate_therm_map_of_properties(T);

        for (i, substance) in sd.substances.iter().enumerate() {
            let map_property_values = sd.therm_map_of_properties_values.get(substance).unwrap();

            let dh = map_property_values.get(&DataType::dH).unwrap().unwrap();
            let ds = map_property_values.get(&DataType::dS).unwrap().unwrap();
            let dG = dh - T * ds;
            // gas_correrction = 0 if w = None or if Phase is not Gas, else gas_correrction = R*T*ln(P/101325) + R*T*ln(w[i])
            // if Phase is Gas or Phase is None
            let gas_correrction: f64 = if let Some(ref w) = w {
                let correction: f64 = R * T * f64::ln(self.P / 101325.0) + R * T * f64::ln(w[i]);
                if let Some(phase) = self.vec_of_phases.get(substance).unwrap() {
                    match phase {
                        Phases::Gas => correction,
                        _ => 0.0,
                    }
                } else {
                    correction
                }
            } else {
                0.0
            };
            let dG = dG + gas_correrction;
            self.dG.insert(substance.clone(), dG);
        }
    }
    pub fn calculate_Gibbs_sym(&mut self, T: f64, w: Option<Vec<Expr>>) {
        let sd = &mut self.subdata;
        let _ = sd.extract_all_thermal_coeffs(T);
        let _ = sd.calculate_therm_map_of_sym();
        let T = Expr::Var("T".to_owned());
        let P = Expr::Var("P".to_owned());
        for (i, substance) in sd.substances.iter().enumerate() {
            let map_sym = sd
                .clone()
                .therm_map_of_sym
                .get(&substance.clone())
                .unwrap()
                .clone();
            let dh_sym = *map_sym.get(&DataType::dH_sym).unwrap().clone().unwrap();
            let ds_sym = *map_sym.get(&DataType::dS_sym).unwrap().clone().unwrap();
            let dG_sym = dh_sym - T.clone() * ds_sym;
            let gas_correrction: Expr = if let Some(ref w) = w {
                let correction: Expr =
                    R_sym * T.clone() * Expr::ln(P.clone() / Expr::Const(101325.0))
                        + R_sym * T.clone() * Expr::ln(w[i].clone());
                let gas_correrction: Expr =
                    if let Some(phase) = self.vec_of_phases.get(substance).unwrap() {
                        match phase {
                            Phases::Gas => correction,
                            _ => Expr::Const(0.0),
                        }
                    } else {
                        correction
                    };
                gas_correrction
            } else {
                Expr::Const(0.0)
            };

            let dG_sym = dG_sym.symplify() + gas_correrction.symplify();
            self.dG_sym.insert(substance.clone(), dG_sym.clone());
        }
    }
    pub fn set_P_to_sym(&mut self) {
        let P =self.P;
        for (_, sym_fun) in self.dG_sym.iter_mut() {
            *sym_fun = sym_fun.set_variable("P", P).symplify()
        }
    }
    pub fn calculate_Gibbs_fun(&mut self, T: f64) {
        let sd = &mut self.subdata;
        let _ = sd.extract_all_thermal_coeffs(T);

        let (_Cp, mut dh_vec, mut ds_vec) = calculate_therm_map_of_fun_local(sd).unwrap();
        for (i, substance) in sd.substances.iter().enumerate() {
            let dh = dh_vec[i].take().unwrap();
            let ds = ds_vec[i].take().unwrap();
            let phases = self.vec_of_phases.clone();
            let gas_correction = {
                let substance = substance.clone();

                let P = self.P;
                move |T: f64, w: Option<Vec<f64>>| -> f64 {
                    if let Some(w) = w {
                        let correction = R * T * f64::ln(P / 101325.0) + R * T * f64::ln(w[i]);
                        match phases.get(&substance) {
                            Some(Some(Phases::Gas)) => correction,
                            _ => 0.0,
                        }
                    } else {
                        0.0
                    }
                }
            };
            let dG = Box::new(move |t: f64, w: Option<Vec<f64>>| {
                dh(t) - t * ds(t) + gas_correction(t, w)
            });
            self.dG_fun.insert(substance.clone(), dG);
        }
    }
    ////////////////////////INPUT/OUTPUT////////////////////////////////////////////////////////
    /// Prints the data of the substance to the console
    pub fn pretty_print(&self) -> Result<(), std::io::Error> {
        let sd = &mut self.subdata.clone();
        println!("__________subs properties at {} K__________", self.T);
        let mut table = Table::new();
        // Add the header row
        table.add_row(row!["substance", "Cp", "dH", "dS", "dG",]);
        // Add the data rows
        for (substance, data) in &self.dG {
            let map_property_values = sd.therm_map_of_properties_values.get(substance).unwrap();
            let Cp = map_property_values.get(&DataType::Cp).unwrap().unwrap();
            let dh = map_property_values.get(&DataType::dH).unwrap().unwrap();
            let ds = map_property_values.get(&DataType::dS).unwrap().unwrap();
            let dG = data.clone();

            table.add_row(row![substance, Cp, dh, ds, dG,]);
        }
        table.printstd();
        println!("_____________________________________________________________");
        println!(
            "\n_________ANALYTIC GIBBS FREE ENERGY AT {} K__________________________________",
            self.T
        );
        let mut table2 = Table::new();
        table2.add_row(row!["substance", "dG_sym",]);
        for substance in &self.vec_of_substances {
            let dG_sym = self.dG_sym.get(substance).unwrap().clone();
            table2.add_row(row![substance, dG_sym,]);
        }
        table2.printstd();
        println!("_____________________________________________________________");
        Ok(())
    }
}
/// Creates closures for the thermodynamic properties of the substances
pub fn calculate_therm_map_of_fun_local(
    sd: &mut SubsData,
) -> Result<
    (
        Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>>,
        Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>>,
        Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>>,
    ),
    String,
> {
    let mut vec_of_cp: Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>> = Vec::new();
    let mut vec_of_dh: Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>> = Vec::new();
    let mut vec_of_ds: Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>> = Vec::new();
    for substance in &sd.substances.clone() {
        match sd
            .search_results
            .get_mut(substance)
            .unwrap()
            .get_mut(&WhatIsFound::Thermo)
            .unwrap()
        {
            Some(SearchResult {
                calculator: Some(CalculatorType::Thermo(thermo)),
                ..
            }) => {
                let mut thermo = thermo.clone();

                // Create closures for thermodynamic functions
                if let Err(e) = thermo.create_closures_Cp_dH_dS() {
                    println!(
                        "Warning: Failed to create closures for {}: {}",
                        substance, e
                    );
                    continue;
                }

                match thermo.get_C_fun() {
                    Ok(cp_fun) => {
                        vec_of_cp.push(Some(cp_fun));
                    }
                    Err(e) => {
                        println!(
                            "Warning: Failed to get Cp function for {}: {}",
                            substance, e
                        );
                    }
                };

                match thermo.get_dh_fun() {
                    Ok(dh_fun) => {
                        vec_of_dh.push(Some(dh_fun));
                    }
                    Err(e) => {
                        println!(
                            "Warning: Failed to get dH function for {}: {}",
                            substance, e
                        );
                    }
                };

                match thermo.get_ds_fun() {
                    Ok(ds_fun) => {
                        vec_of_ds.push(Some(ds_fun));
                    }
                    Err(e) => {
                        println!(
                            "Warning: Failed to get dS function for {}: {}",
                            substance, e
                        );
                    }
                };
            }
            _ => {
                continue;
            }
        }
    }
    Ok((vec_of_cp, vec_of_dh, vec_of_ds))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::User_substances::LibraryPriority;
    use approx::assert_relative_eq;
    fn setup_test_data() -> Thermodynamics {
        let subs = vec!["CO".to_string(), "CO2".to_string()];
        let mut subdata = SubsData::new();
        subdata.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        subdata.substances = subs.clone();
        subdata.search_substances();

        let mut thermo = Thermodynamics::new();
        // Set up basic parameters
        thermo.set_T(400.0);
        thermo.set_P(101325.0, None);
        // Add test substances
        thermo.vec_of_substances = subs.clone();
        // Set up phases
        thermo
            .vec_of_phases
            .insert("CO".to_string(), Some(Phases::Gas));
        thermo
            .vec_of_phases
            .insert("CO2".to_string(), Some(Phases::Gas));
        // Set up SubsData with mock thermodynamic data
        thermo.subdata = subdata;
        thermo
    }

    #[test]
    fn test_new() {
        let thermo = Thermodynamics::new();
        assert!(thermo.vec_of_substances.is_empty());
        assert!(thermo.vec_of_phases.is_empty());
        assert_eq!(thermo.T, 298.15);
        assert_eq!(thermo.P, 1e5);
        assert!(thermo.dG.is_empty());
        assert!(thermo.dG_fun.is_empty());
        assert!(thermo.dG_sym.is_empty());
    }

    #[test]
    fn test_set_t_and_p() {
        let mut thermo = Thermodynamics::new();

        // Test setting temperature
        thermo.set_T(500.0);
        assert_eq!(thermo.T, 500.0);

        // Test setting pressure
        thermo.set_P(2e5, None);
        assert_eq!(thermo.P, 2e5);
    }

    #[test]
    fn test_calculate_gibbs_free_energy() {
        let mut thermo = setup_test_data();

        // Calculate Gibbs free energy without concentration correction
        thermo.calculate_Gibbs_free_energy(400.0, None);

        // Check results
        for substance in &thermo.vec_of_substances {
            assert!(thermo.dG.contains_key(substance));
        }

        // Test with concentration correction
        let concentrations = vec![0.8, 0.2]; // Mock mole fractions
        thermo.calculate_Gibbs_free_energy(400.0, Some(concentrations.clone()));
        println!("thermo.dG: {:?} \n", thermo.dG);
        // Check results with gas correction
        for (_i, substance) in thermo.vec_of_substances.iter().enumerate() {
            // dG = dH - T*dS + RT*ln(P/P°) + RT*ln(w_i)
            // let gas_correction = R * 400.0 * f64::ln(thermo.P / 101325.0) + R * 400.0 * f64::ln(concentrations[i]);
            // let expected_dg = -110.0 - 400.0 * 200.0 + gas_correction;
            //assert!((thermo.dG[substance] - expected_dg).abs() < 1e-10);
            assert!(thermo.dG.contains_key(substance));
        }
    }

    #[test]
    fn test_calculate_gibbs_sym() {
        let mut thermo = setup_test_data();

        // Calculate symbolic Gibbs free energy
        thermo.calculate_Gibbs_sym(400.0, None);
        // Calculate Gibbs free energy without concentration correction
        thermo.calculate_Gibbs_free_energy(400.0, None);
        // Check results
        for substance in &thermo.vec_of_substances {
            assert!(thermo.dG_sym.contains_key(substance));

            // Verify the structure of the symbolic expression
            // For our test data, it should be a simple expression: -110 - T*200
            let expr = &thermo.dG_sym[substance];
            let dG_function = expr.lambdify1D();
            let dG_i_from_sym = dG_function(400.0);
            let dG_direct = thermo.dG[substance];
            assert_relative_eq!(dG_i_from_sym, dG_direct, epsilon = 1e-10);

            // This is a simplified check - in a real test you might want to evaluate the expression
            // or check its structure more thoroughly
            //   assert!(expr.to_string().contains("-110"));
        }

        // Test with concentration correction
        let concentrations = vec![Expr::Var("w1".to_owned()), Expr::Var("w2".to_owned())];
        thermo.calculate_Gibbs_sym(400.0, Some(concentrations));

        // Check that the expressions now include concentration terms
        for substance in &thermo.vec_of_substances {
            let expr = &thermo.dG_sym[substance];
            let expr_str = expr.to_string();

            // The expression should now include ln terms for the gas correction
            assert!(expr_str.contains("ln"));
        }
    }

    #[test]
    fn test_calculate_gibbs_fun() {
        let mut thermo = setup_test_data();

        // Calculate Gibbs free energy functions
        thermo.calculate_Gibbs_fun(500.0);
        thermo.calculate_Gibbs_free_energy(500.0, None);
        // Check that functions were created for each substance
        for substance in &thermo.vec_of_substances {
            assert!(thermo.dG_fun.contains_key(substance));

            // Test the function at a specific temperature
            let dg_fun = &thermo.dG_fun[substance];
            let dg_value = dg_fun(500.0, None);
            let dG_direct = thermo.dG[substance];
            assert_relative_eq!(dg_value, dG_direct, epsilon = 1e-10);
            // For our test data with constant values:

            // Test with concentration correction
            let dg_value_with_conc = dg_fun(500.0, Some(vec![0.5, 0.5]));

            // The value should now include the gas correction term
            assert!(dg_value_with_conc != dg_value);
        }
    }

    #[test]
    fn test_calculate_therm_map_of_fun_local() {
        let mut thermo = setup_test_data();
        let _ = thermo.subdata.extract_all_thermal_coeffs(400.0);
        thermo.subdata.print_search_summary();

        let result = calculate_therm_map_of_fun_local(&mut thermo.subdata);

        assert!(result.is_ok());
        let (cp_funs, dh_funs, ds_funs) = result.unwrap();

        // Check that we got the expected number of functions
        assert_eq!(cp_funs.len(), thermo.vec_of_substances.len());
        assert_eq!(dh_funs.len(), thermo.vec_of_substances.len());
        assert_eq!(ds_funs.len(), thermo.vec_of_substances.len());

        // Test the functions
        for i in 0..thermo.vec_of_substances.len() {
            let cp_fun = cp_funs[i].as_ref().unwrap();
            let dh_fun = dh_funs[i].as_ref().unwrap();
            let ds_fun = ds_funs[i].as_ref().unwrap();

            // For our test data, these should return constant values
            assert!(cp_fun(400.0) > 0.0);
            assert!(dh_fun(400.0) < 0.0);
            assert!(ds_fun(400.0) > 0.0);
        }
    }
    #[test]
    fn test_all_features() {
        let subs = vec!["CO".to_string(), "CO2".to_string()];
        // calling instance of strucutre SubsData created to search substances data in the databases, store
        // search results and calculate thermo properties
        let mut subdata = SubsData::new();
        // Set up library priorities
        subdata.set_multiple_library_priorities(
            vec!["NASA_gas".to_string()],
            LibraryPriority::Priority,
        );
        subdata.substances = subs.clone();
        // Perform the search
        subdata.search_substances();
        // Calling instance of structure Thermodynamics to calculate thermo dG
        let mut thermo = Thermodynamics::new();
        // Set up basic parameters
        thermo.set_T(400.0);
        thermo.set_P(101325.0, None);
        // Add test substances
        thermo.vec_of_substances = subs.clone();
        // Set up phases
        thermo
            .vec_of_phases
            .insert(subs[0].clone(), Some(Phases::Gas));
        thermo
            .vec_of_phases
            .insert(subs[1].clone(), Some(Phases::Gas));
        // savng  the search results in the structure Thermodynamics
        thermo.subdata = subdata;
        // Calculate Gibbs free energy
        thermo.calculate_Gibbs_free_energy(400.0, None);
        // Calculate symbolic Gibbs free energy
        thermo.calculate_Gibbs_sym(400.0, None);
        // getting the results of the calculation
        thermo.calculate_Gibbs_fun(400.0);
        let map_of_gibbs = thermo.dG.clone();
        let map_of_gibbs_sym = thermo.dG_sym.clone();
        let map_of_gibbs_fun = thermo.clone().dG_fun;
        for substance in &thermo.vec_of_substances {
            println!("substance: {:?}", substance);
            println!("map_of_gibbs: {:?}", map_of_gibbs[substance]);
            println!("map_of_gibbs_sym: {:?}", map_of_gibbs_sym[substance]);
        }
        for substance in &thermo.vec_of_substances {
            let dG_value = map_of_gibbs[substance];
            let dG_sym = map_of_gibbs_sym[substance].clone();
            let dG_from_sym = dG_sym.lambdify1D()(400.0);
            assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
            let dG_from_fun = map_of_gibbs_fun.get(substance).unwrap()(400.0, None);
            assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);
        }

        println!("map_of_gibbs: {:?} \n", map_of_gibbs);
        let _ = thermo.pretty_print();
    }

    #[test]
    fn test_all_features_witn_conc() {
         // calculating Gibbs free energy withoun concentration correction RT*ln(P/P°) + RT*ln(w_i)
         let subs =  vec!["CO".to_string(), "CO2".to_string()];
         // calling instance of strucutre SubsData created to search substances data in the databases, store 
         // search results and calculate thermo properties
         let mut subdata = SubsData::new();
         // Set up library priorities
         subdata.set_multiple_library_priorities(
             vec!["NASA_gas".to_string()],
             LibraryPriority::Priority,
         );
         subdata.substances = subs.clone();
         // Perform the search
         subdata.search_substances();
         // Calling instance of structure Thermodynamics to calculate thermo dG 
         let mut thermo = Thermodynamics::new();  
         // Set up basic parameters
         thermo.set_T(400.0);
         thermo.set_P(101325.0, None); 
         let concentration = Some(vec![0.5, 0.5]);
         // Add test substances
         thermo.vec_of_substances = subs.clone();   
         // Set up phases
         thermo.vec_of_phases.insert(subs[0].clone(), Some(Phases::Gas));
         thermo.vec_of_phases.insert(subs[1].clone(), Some(Phases::Gas));
         // savng  the search results in the structure Thermodynamics
         thermo.subdata = subdata;
         // Calculate Gibbs free energy
         thermo.calculate_Gibbs_free_energy(400.0,  concentration.clone());
         // Calculate symbolic Gibbs free energy
         thermo.calculate_Gibbs_sym(400.0, Some(vec![Expr::Var("w1".to_string()), Expr::Var("w2".to_string())]));
         thermo.set_P_to_sym();
         // getting the results of the calculation
         thermo.calculate_Gibbs_fun( 400.0);
         let map_of_gibbs = thermo.dG.clone();
         let map_of_gibbs_sym = thermo.dG_sym.clone();
         let map_of_gibbs_fun = thermo.clone().dG_fun;
        for substance in &thermo.vec_of_substances {
                 println!("substance: {:?}", substance);
                 println!("map_of_gibbs: {:?}", map_of_gibbs[substance]);
                 println!("map_of_gibbs_sym: {:?}", map_of_gibbs_sym[substance]);
             
             }
        for substance in &thermo.vec_of_substances {
                 let dG_value = map_of_gibbs[substance];
     
                 let dG_from_fun = map_of_gibbs_fun.get(substance).unwrap()(400.0, concentration.clone());
                 assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);
     
                 let dG_sym = map_of_gibbs_sym[substance].clone();
                 let dG_from_sym  = dG_sym.lambdify_owned(vec!["T", "w1", "w2"]);
                 let dG_from_sym  = dG_from_sym(vec![400.0, 0.5, 0.5]);
                 assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
             
             }
             
         //println!("map_of_gibbs: {:?} \n", map_of_gibbs);
        // let _ = thermo.pretty_print();
     
    }
}
