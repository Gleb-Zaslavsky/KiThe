use crate::Examples::ClassicalThermodynamics_examples::SubsData_examples;
use crate::Examples::NIST_examples::NIST_examples;
use crate::Examples::SubstanceDataCollecting::collecting_thermo_data;
use crate::Examples::kinetics_examples::kin_examples;
use crate::Examples::thermo_examples::thermo_examples;
use std::io::{self, Write};

pub fn examples_menu() {
    loop {
        println!("\n=== Examples ===");
        println!("1. Kinetics Examples");
        println!("2. Thermodynamics Examples");
        println!("3. NIST Examples");
        println!("4. Substance Data Collection");
        println!("5. Classical Thermodynamics");
        println!("0. Back to main menu");
        print!("Enter your choice: ");
        io::stdout().flush().unwrap();

        let choice = get_user_input();
        match choice.trim() {
            "1" => kin_examples(3),
            "2" => thermo_examples(6),
            "3" => NIST_examples(6),
            "4" => collecting_thermo_data(6),
            "5" => SubsData_examples(6),
            "0" => break,
            _ => println!("Invalid choice. Please try again."),
        }
    }
}

fn get_user_input() -> String {
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed to read input");
    input
}
