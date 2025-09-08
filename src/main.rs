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

use ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;
use cli::cli_main::run_interactive_menu;
use std::env;
use std::path::PathBuf;

pub fn main() {
    run_interactive_menu();
    /*
    let args: Vec<String> = env::args().collect();

    if args.len() > 1 {
        // Command line argument provided - direct file solving
        let file_path = PathBuf::from(&args[1]);
        println!("Solving reactor problem from file: {:?}", file_path);
        let mut reactor = SimpleReactorTask::new();
        reactor.solve_from_file(file_path);
    } else {
        // No arguments - show interactive menu
        run_interactive_menu();
    }
    */
}
