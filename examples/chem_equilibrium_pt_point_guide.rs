//! Guide: one fixed-pressure, fixed-temperature equilibrium point.
//!
//! The pipeline resolves local thermochemistry first, then solves the declared
//! phase system transactionally. No mutable solver state or database payload
//! is exposed to the caller.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumConditions, EquilibriumSolveOptions, PhaseEquilibriumPipelineRequest,
    SubstanceSystemSpecBuilder, SubstancesContainer,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
        "H2".to_string(),
        "O2".to_string(),
        "H2O".to_string(),
    ]))
    .with_library_priorities(vec!["NASA_gas".to_string()])
    .with_search_in_nist(false)
    .build()?;

    let outcome = PhaseEquilibriumPipelineRequest::new(
        spec,
        vec![0.1, 0.05, 1.9],
        EquilibriumConditions::new(2_500.0, 101_325.0, 101_325.0)?,
    )
    .with_solve_options(EquilibriumSolveOptions::new().with_production_cascade())
    .solve()?;

    println!("lookup provenance: {:?}", outcome.lookup_report());
    println!("equilibrium point:\n{}", outcome.solution());
    Ok(())
}
