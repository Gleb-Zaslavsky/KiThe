use crate::Thermodynamics::DBhandlers::NIST_parser::{ Phase, SearchType, NistParser};


pub fn NIST_examples(nist_examples:usize) {
  match nist_examples {
      
  
        0 =>  {
                let parser = NistParser::new();

                // Example usage
                let substance = "CH4";
                match parser.get_data(substance, SearchType::All, Phase::Gas) {
                    Ok(data) =>{ println!("Data for {}: {:?}", substance, data);
                
                    data.pretty_print();
                    #[allow(non_snake_case)]
                    let (Cp, dh, ds) = data.caclc_cp_dh_ds(298.15).expect("Error calculating cp, dh, ds");
                    println!("Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}", Cp, dh, ds);
                },
                    Err(e) => eprintln!("Error: {}", e),
                    
                }
        }

   1=> {
                let parser = NistParser::new();
            let substance = "NaCl";
            match parser.get_data(substance, SearchType::All, Phase::Solid) {
                Ok(data) =>{ println!("Data for {}: {:?}", substance, data);
            
                data.pretty_print();
                #[allow(non_snake_case)]
                let (Cp, dh, ds) = data.caclc_cp_dh_ds(298.15).expect("Error calculating cp, dh, ds");
                println!("Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}", Cp, dh, ds);
            },
                Err(e) => eprintln!("Error: {}", e),
                
            }

            match parser.get_data(substance, SearchType::All, Phase::Liquid) {
                Ok(data) =>{ println!("Data for {}: {:?}", substance, data);
            
                data.pretty_print();
                #[allow(non_snake_case)]
                let (Cp, dh, ds) = data.caclc_cp_dh_ds(1200.15).expect("Error calculating cp, dh, ds");
                println!("Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}", Cp, dh, ds);
            },
                Err(e) => eprintln!("Error: {}", e),
                
            }
        }
    _=> {println!("non existing examples");}
    }
}

