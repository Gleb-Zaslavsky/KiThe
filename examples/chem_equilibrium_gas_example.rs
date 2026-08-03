//! Guide: canonical gas-phase equilibrium solve.
//!
//! The example uses the same production facade as resolved multiphase users:
//! declare a phase specification, let the repository-backed pipeline resolve
//! its thermochemistry, and inspect the immutable outcome. The handwritten
//! nonlinear solvers remain fallback implementations selected by the typed
//! cascade; they are not the orchestration API shown to users.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumConditions, EquilibriumSolveOptions, PhaseEquilibriumPipelineRequest,
    SubstanceSystemSpecBuilder, SubstancesContainer,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
        "CO".to_string(),
        "CO2".to_string(),
        "O2".to_string(),
    ]))
    .with_library_priorities(vec!["NASA_gas".to_string()])
    .with_search_in_nist(false)
    .build()?;

    let outcome = PhaseEquilibriumPipelineRequest::new(
        spec,
        vec![0.25, 0.25, 0.5],
        EquilibriumConditions::new(1_500.0, 101_325.0, 101_325.0)?,
    )
    .with_solve_options(EquilibriumSolveOptions::new().with_production_cascade())
    .solve()?;

    println!("lookup report: {:?}", outcome.lookup_report());
    println!("accepted solution:\n{}", outcome.solution());
    Ok(())
}
