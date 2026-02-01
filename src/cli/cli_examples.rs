use crate::Examples::ClassicalThermodynamics_examples::SubsData_examples;
use crate::Examples::NIST_examples::NIST_examples;
use crate::Examples::SubstanceDataCollecting::collecting_thermo_data;
use crate::Examples::kinetics_examples::kin_examples;
use crate::Examples::thermo_examples::thermo_examples;
use dialoguer::{Select, theme::ColorfulTheme};

pub fn examples_menu() {
    let theme = ColorfulTheme::default();
    loop {
        let items = vec![
            "Kinetics Examples",
            "Thermodynamics Examples",
            "NIST Examples",
            "Substance Data Collection",
            "Classical Thermodynamics",
            "Back to main menu",
        ];

        let selection = Select::with_theme(&theme)
            .with_prompt("=== Examples ===")
            .items(&items)
            .default(0)
            .interact()
            .unwrap();

        match selection {
            0 => kin_examples(3),
            1 => thermo_examples(6),
            2 => NIST_examples(6),
            3 => collecting_thermo_data(6),
            4 => SubsData_examples(6),
            5 => break,
            _ => {}
        }
    }
}
