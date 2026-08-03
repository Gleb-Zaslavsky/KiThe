//! Guide: one repository-resolved phase system through the fixed-`P,T` facade.
//!
//! The example deliberately uses only local NASA data. The typed pipeline
//! resolves thermochemistry first and then consumes the immutable resolved
//! system without exposing the lookup maps to the caller or reopening mutable
//! library state while solving.

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
        vec![2.0, 1.0, 0.0],
        EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0)?,
    )
    .with_solve_options(EquilibriumSolveOptions::new().with_production_cascade())
    .solve()?;

    println!("lookup report: {:?}", outcome.lookup_report());
    println!("accepted solution:\n{}", outcome.solution());
    Ok(())
}
