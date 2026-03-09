use super::cli_bvp_reactor::reactor_menu;
use super::cli_examples::examples_menu;
use super::cli_solid_state_ivp::solid_state_ivp_menu;
use crate::gui::gui_main::gui_main;
use dialoguer::{Select, theme::ColorfulTheme};
pub fn run_interactive_menu() {
    let theme = ColorfulTheme::default();
    loop {
        println!(
            "\x1b[34m\n Wellcome to KiThe: Toolkit for chemical engineering, combustion,\n
    chemical kinetics, chemical thermodynamics and more \n
    (c) Gleb E. Zaslavsky, 2024-2026 \n \x1b[0m"
        );

        let items = vec![
            "Reactor BVP Problems",
            "Solid State IVP Problems",
            "Examples",
            "GUI",
            "Exit",
        ];

        let selection = Select::with_theme(&theme)
            .with_prompt("Select an option")
            .items(&items)
            .default(0)
            .interact()
            .unwrap();

        match selection {
            0 => reactor_menu(),
            1 => solid_state_ivp_menu(),
            2 => examples_menu(),
            3 => {
                gui_main().unwrap();
            }
            4 => {
                println!("Goodbye!");
                break;
            }
            _ => {}
        }
    }
}
/* colors
Blue (\x1b[34m) - Welcome header text

Yellow (\x1b[33m) - Menu options (1, 2, 0)

Cyan (\x1b[36m) - "Enter your choice:" prompt

Reset (\x1b[0m) - Returns to normal color after each colored section
*/
