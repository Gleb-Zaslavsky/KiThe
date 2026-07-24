//! A small condensed-reactor IVP preview example.
//!
//! This example builds a minimal condensed-phase task, folds it into the typed
//! IVP editor model, and prints the preview tables that the GUI would show
//! before a real solve.

use KiThe::Kinetics::User_reactions::KinData;
use KiThe::ReactorsIVP::SimpleReactorIVP::{IvpError, SimpleReactorTask};
use KiThe::gui::reactor_ivp_gui::{IvpGuiConfig, build_reactor_ivp_task_preview_snapshot};
use std::collections::HashMap;

fn main() -> Result<(), IvpError> {
    let mut task = SimpleReactorTask::new();
    task.kindata = KinData::from_direct_reactions(vec!["A=>B".to_string()])
        .map_err(|err| IvpError::InvalidConfiguration(err.to_string()))?;
    task.set_problem_name("Condensed burn");
    task.set_problem_description("Preview-only IVP example");
    task.set_density(1200.0)?;
    task.set_transport_properties(0.25, 1000.0);
    task.m = 0.015;
    task.L = 0.05;
    task.set_initial_conditions(HashMap::from([
        ("T".to_string(), 450.0),
        ("q".to_string(), 0.15),
        ("A".to_string(), 0.7),
        ("B".to_string(), 0.3),
    ]));
    task.set_thermal_effects(vec![-2.0e5]);

    let config = IvpGuiConfig::from_task(&task);
    let preview = build_reactor_ivp_task_preview_snapshot(&task, &config)?;

    println!("=== Reactor IVP preview ===");
    println!("{}", tabled::Table::new(preview.summary_rows.clone()));
    println!("=== Reactor IVP equations ===");
    println!("{}", tabled::Table::new(preview.equation_rows.clone()));

    Ok(())
}
