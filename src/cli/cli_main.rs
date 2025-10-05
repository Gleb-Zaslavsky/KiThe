use super::cli_bvp_reactor::reactor_menu;
use super::cli_examples::examples_menu;
use super::cli_solid_state_ivp::solid_state_ivp_menu;
use std::io::{self, Write};
use crate::gui::gui_main::gui_main;
pub fn run_interactive_menu() {
    loop {
        show_main_menu();
        let choice = get_user_input();

        match choice.trim() {
            "1" => reactor_menu(),
            "2" => solid_state_ivp_menu(),
            "3" => examples_menu(),
            "4"=> gui_main().unwrap(),
            "0" => {
                println!("Goodbye!");
                break;
            }
            _ => println!("Invalid choice. Please try again."),
        }
    }
}
/* colors
Blue (\x1b[34m) - Welcome header text

Yellow (\x1b[33m) - Menu options (1, 2, 0)

Cyan (\x1b[36m) - "Enter your choice:" prompt

Reset (\x1b[0m) - Returns to normal color after each colored section
*/
fn show_main_menu() {
    println!(
        "\x1b[34m\n Wellcome to KiThe: Toolkit for chemical engineering, combustion,\n
    chemical kinetics, chemical thermodynamics and more \n
    (c) Gleb E. Zaslavsky, 2024 \n \x1b[0m"
    );
    println!("\x1b[33m1. Reactor BVP Problems\x1b[0m");
    println!("\x1b[33m2. Solid State IVP Problems\x1b[0m");
    println!("\x1b[33m3. Examples\x1b[0m");
    println!("\x1b[33m4. GUI\x1b[0m");
    println!("\x1b[33m0. Exit\x1b[0m");
    print!("\x1b[36mEnter your choice: \x1b[0m");
    io::stdout().flush().unwrap();
}

fn get_user_input() -> String {
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed to read input");
    input
}
