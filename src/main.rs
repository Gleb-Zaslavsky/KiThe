#[allow(non_snake_case)]
pub mod Examples;
#[allow(non_snake_case)]
pub mod Kinetics;
#[allow(non_snake_case)]
pub mod ReactorsBVP;
#[allow(non_snake_case)]
pub mod Thermodynamics;
#[allow(non_snake_case)]
pub mod Utils;
#[allow(non_snake_case)]
pub mod cli;
pub mod gui;

use cli::cli_main::run_interactive_menu;
use gui::gui_main::gui_main;

pub mod simple_combustion_models;
pub fn main() {
    run_interactive_menu();
    // gui_main().unwrap();
}
