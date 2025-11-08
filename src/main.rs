#[allow(non_snake_case)]
pub mod Examples;
#[allow(non_snake_case)]
pub mod Kinetics;
#[allow(non_snake_case)]
pub mod ReactorsBVP;
#[allow(non_snake_case)]
pub mod ReactorsIVP;
#[allow(non_snake_case)]
pub mod Thermodynamics;
#[allow(non_snake_case)]
pub mod Utils;
#[allow(non_snake_case)]
pub mod cli;
pub mod gui;
use cli::cli_main::run_interactive_menu;
use gui::gui_main::gui_main;
pub mod library_manager;
pub mod settings;
pub mod simple_combustion_models;
pub fn main() {
    run_interactive_menu();
    // gui_main().unwrap();
}

/*
TODO list
 - transport gui fix for Aramco transport
 - NISTparser - opportunity to load from library
 - gui for NIST parser
 - All substances comparing gui instrument
 - fix chemical equilibrium 
 - add simple combustion models
 - add simple combustion models gui
 - interfaces for temperature range Cp, dH, dS
 - integral mean properties
 - non adjacent temperature intervals 




*/