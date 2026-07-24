//! Guide: canonical gas-phase equilibrium solve.
//!
//! This is the production equilibrium path: build a gas-phase solver, run the
//! solve, and inspect the accepted solution snapshot plus the reconstructed
//! mole table.

use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_workflows::gas_solver;

fn main() {
    let mut solver = gas_solver(
        vec!["CO".to_string(), "CO2".to_string(), "O2".to_string()],
        1500.0,
        101_325.0,
        Solvers::LM,
        Some("info"),
        true,
    )
    .expect("failed to prepare equilibrium solver");

    solver.solve().expect("equilibrium solve failed");
    let accepted = solver
        .accepted_solution()
        .expect("equilibrium solve did not publish an accepted solution");

    println!("accepted solution: {accepted:?}");
    println!("{}", solver.moles_table());
}
