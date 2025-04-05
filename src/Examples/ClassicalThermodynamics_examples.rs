use std::vec;

use crate::Thermodynamics::ChemEquilibrium::ClassicalThermodynamics::Thermodynamics;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use crate::Thermodynamics::User_substances::{
 LibraryPriority, SubsData,  Phases
};
use approx::assert_relative_eq;
pub fn SubsData_examples(thermotask: usize) {
    match thermotask {
        0 => { // calculating Gibbs free energy withoun concentration correction RT*ln(P/P°) + RT*ln(w_i)
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
        // Add test substances
        thermo.vec_of_substances = subs.clone();   
        // Set up phases
        thermo.vec_of_phases.insert(subs[0].clone(), Some(Phases::Gas));
        thermo.vec_of_phases.insert(subs[1].clone(), Some(Phases::Gas));
        // savng  the search results in the structure Thermodynamics
        thermo.subdata = subdata;
        // Calculate Gibbs free energy
        thermo.calculate_Gibbs_free_energy(400.0, None);
        // Calculate symbolic Gibbs free energy
       thermo.calculate_Gibbs_sym(400.0, None);
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
            let dG_sym = map_of_gibbs_sym[substance].clone();
            let dG_from_sym  = dG_sym.lambdify1D()(400.0);
            assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
            let dG_from_fun = map_of_gibbs_fun.get(substance).unwrap()(400.0, None);
            assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);
         
        }
        
       println!("map_of_gibbs: {:?} \n", map_of_gibbs);
       let _ = thermo.pretty_print();
    }
    1 =>{ // calculating Gibbs free energy withoun concentration correction RT*ln(P/P°) + RT*ln(w_i)
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
                let dG_sym = dG_sym.set_variable("P", 101325.0);
          
                let dG_from_sym  = dG_sym.lambdify_owned(vec!["T", "w1", "w2"]);
                let dG_from_sym  = dG_from_sym(vec![400.0, 0.5, 0.5]);
                assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
            
            }
            
        //println!("map_of_gibbs: {:?} \n", map_of_gibbs);
       // let _ = thermo.pretty_print();


    }
    2 =>{ 
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
    _ => {}
    }
}