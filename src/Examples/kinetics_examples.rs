

pub fn kin_examples(kintask:usize) {
    //

    match kintask {
        0 => {
            // STOICHIOMETRIC ANALYSIS 
            use crate::Kinetics::stoichiometry_analyzer::StoichAnalyzer;
            let mut  StoichAnalyzer_instance = StoichAnalyzer::new();
            let reactions_: Vec<&str> = vec!["A=2BM)", "B=>A + 3C_DUP", "2B+A=D**0.5"];
            let reaction = reactions_.iter().map(|s| s.to_string()).collect();
            StoichAnalyzer_instance.reactions = reaction;
            StoichAnalyzer_instance.search_substances();
            StoichAnalyzer_instance.analyse_reactions();
          
            let stecheo_matrx =  StoichAnalyzer_instance.stecheo_matrx;

            let result = [
                [-1.0, 2.0, 0.0, 0.0],
                [1.0, -1.0, 3.0, 0.0],
                [-1.0, -2.0, 0.0, 1.0],
            ];
            let result: Vec<Vec<f64>> = result.iter().map(|row| row.to_vec()).collect();
            assert_eq!(stecheo_matrx, result);
            println!("substances: {:?}", StoichAnalyzer_instance.substances);
            println!("stecheo_matrx {:?}", stecheo_matrx);

        }
        1 => { //  Calculation of atomic composition, molar masses and matrix of atomic composition
            use crate::Kinetics::molmass::{calculate_molar_mass, parse_formula, calculate_molar_mass_of_vector_of_subs, create_elem_composition_matrix};
            let formula = "C6H8O6";
            let (molar_mass, element_composition) = calculate_molar_mass(formula.to_string(), None); 
            println!("Element counts: {:?}", element_composition);
            println!("Molar mass: {:?} g/mol", molar_mass);
            
            let formula = "Na(NO3)2".to_string();
            let atomic_composition = parse_formula(formula, None);
            println!("{:?}", atomic_composition);


            let vec_of_formulae = vec!["H2O", "NaCl", "C6H8O6", "Ca(NO3)2"];
            let expected_molar_masses = vec![18.01528, 58.44316, 176.12, 164.093];
            let calculated_molar_masses = calculate_molar_mass_of_vector_of_subs(vec_of_formulae, None);
    
            for (i, &expected_molar_mass) in expected_molar_masses.iter().enumerate() {
                println!("molar mass: {:?} g/mol", calculated_molar_masses[i]);
                assert!((calculated_molar_masses[i] - expected_molar_mass).abs() < 1e-2);
            }

            let vec_of_formulae = vec!["H2O", "NaCl", "C3H8", "CH4"]; // 5 elements
            let (matrix, vec_of_formulae) = create_elem_composition_matrix(vec_of_formulae, None);
            println!("{}", matrix);
            println!("{:?}", vec_of_formulae);
        }
        2 => {// kinetics library api
            use crate::Kinetics::kinetics_lib_api::KineticData;
            let mut kin_instance = KineticData::new();
            // collecting reaction data for library name lib
            let lib = "NUIG";

            kin_instance.open_json_files(lib);
            // veiew all reactions in library
            kin_instance.print_all_reactions();
      
            let reaction1 = kin_instance.search_reactdata_by_reaction_id("1");
            println!("reaction1: {:?}", reaction1);
            // search reactions by substances 
            kin_instance.search_reaction_by_reagents_and_products(vec!["CO".to_string()]  );
            println!("reactions where CO is product: {:?}", kin_instance.FoundReactionsByProducts);
            println!("reactions where CO is reagent: {:?}", kin_instance.FoundReactionsByReagents);

        }
        3 => { // mechanism construction
            use crate::Kinetics::mechfinder_api::Mechanism_search;
            let mut mech_search = Mechanism_search::new(
                vec!["O".to_string(), "NH3".to_string(), "NO".to_string()],
                "NUIG".to_string(),
            );
    
            let (mechanism, reactants, vec_of_reactions) = mech_search.mechfinder_api();
    
            println!("mechanism (reaction ID's) : {:?}", mechanism);
            println!("reactants: {:?}", reactants);
            println!("reaction data: {:?}", vec_of_reactions);
            println!("vector of ReactionData structs with parsed data: {:#?}", mech_search.reactdata);
        }
        4 => { // struct KinData aggregate most of functionality of the Kinetics module 
            use crate::Kinetics::User_reactions::KinData;
            // let our journey begin with a new instance of KinData
            let mut kd = KinData::new();
            // set the shortcut reactions for our KineticData instance
            // it means we want reactions from Cantera sub-librarie from number 1 to number 10
            kd.set_reactions_from_shortcut_range("C1..C10".to_string());
            // searching for reactions in data base
            kd.get_reactions_from_shortcuts();
            kd.kinetic_main(); // parsing reaction data into structures and stoichometric calculations under one hood
            kd.pretty_print_kindata();
           // kd.pretty_print_kindata_verbose();

        }
        5=> {
            use crate::Kinetics::User_reactions::KinData;
            // let our journey begin with a new instance of KinData
            let mut kd = KinData::new();
            let task_substances = vec!["O".to_string(), "NH3".to_string(), "NO".to_string()];
            let task_library = "NUIG".to_string();
            kd.construct_mechanism(task_substances, task_library);
            kd.kinetic_main(); // parsing reaction data into structures and stoichometric calculations under one hood
            kd.pretty_print_kindata();
            println!("vector of reactions \n\n {:#?}", kd.vec_of_equations);
            println!("vector of substances \n\n {:#?}", kd.substances);

        }
        6=> { // parsing kinetcis data from json
            use crate::Kinetics::mechfinder_api::kinetics::{PressureStruct, ThreeBodyStruct};
            let ThreeBodyStruct_test_data: &str = r#"{"Arrenius": [4.577e+19, -1.4, 436705.19999999995],
            "eff": {"H2": 2.5, "H2O": 12.0, "CO": 1.9, "CO2": 3.8, "HE": 0.83, "CH4": 2.0, "C2H6": 3.0} }"#;
           
            let tdata =
            serde_json::from_value::<ThreeBodyStruct>(
                serde_json::from_str(ThreeBodyStruct_test_data).unwrap(),
            );

            if let Ok(tdata)= tdata {
                println!("tdata: {:#?}", tdata);}
                else { panic!("Expected ThreeBodyStruct variant");}
              
              const PRES_TESTING_JSON: &str = r#"{
                        "Arrenius":{"0.01": [2.89e+40, -9.76, 140552.983],
                                    "0.1": [1.8e+44, -10.5, 154800.281]
                                    }}"#;
              let pres_val  = serde_json::from_str(PRES_TESTING_JSON).expect("Error parsing JSON: {err:?}");                      
              let pres_data = serde_json::from_value::<PressureStruct>(pres_val);
            
                if let Ok(pres)= pres_data {
                    println!("pres_data: {:#?}", pres);
                    // Success
                } else { panic!("Expected Pressure variant");}

            }


    
        _ => {
            println!("Wrong task number");
        }
    }
}