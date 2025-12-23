use crate::Thermodynamics::DBhandlers::NIST_parser::{NistParser, Phase, SearchType};
use crate::Thermodynamics::DBhandlers::NISTdata::NISTdata;

pub fn NIST_examples(nist_examples: usize) {
    match nist_examples {
        0 => {
            let parser = NistParser::new();
            let mut nist_data = NISTdata::new();

            // Example usage
            let substance = "CH4";
            match parser.get_data(substance, SearchType::All, Phase::Gas) {
                Ok(data) => {
                    println!("Data for {}: {:?}", substance, data);
                    nist_data.input = data;
                    let _ = nist_data.pretty_print();
                    let _ = nist_data.extract_coefficients(298.15);
                    let _ = nist_data.calculate_cp_dh_ds(298.15);
                    println!(
                        "Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}",
                        nist_data.Cp, nist_data.dh, nist_data.ds
                    );
                }
                Err(e) => eprintln!("Error: {}", e),
            }
        }

        1 => {
            let parser = NistParser::new();
            let mut nist_data = NISTdata::new();
            let substance = "NaCl";

            match parser.get_data(substance, SearchType::All, Phase::Solid) {
                Ok(data) => {
                    println!("Data for {}: {:?}", substance, data);
                    nist_data.input = data;
                    let _ = nist_data.pretty_print();
                    let _ = nist_data.extract_coefficients(298.15);
                    let _ = nist_data.calculate_cp_dh_ds(298.15);
                    println!(
                        "Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}",
                        nist_data.Cp, nist_data.dh, nist_data.ds
                    );
                }
                Err(e) => eprintln!("Error: {}", e),
            }

            match parser.get_data(substance, SearchType::All, Phase::Liquid) {
                Ok(data) => {
                    println!("Data for {}: {:?}", substance, data);
                    nist_data.input = data;
                    let _ = nist_data.pretty_print();
                    let _ = nist_data.extract_coefficients(1200.15);
                    let _ = nist_data.calculate_cp_dh_ds(1200.15);
                    println!(
                        "Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}",
                        nist_data.Cp, nist_data.dh, nist_data.ds
                    );
                }
                Err(e) => eprintln!("Error: {}", e),
            }
        }
        _ => {
            println!("non existing examples");
        }
    }
}
