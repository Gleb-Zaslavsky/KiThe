//! Guide: a fixed-pressure temperature range with accepted-state continuation.
//!
//! The immutable lookup result, layout, and prepared formulation are reused
//! across points. The range publishes only if every requested temperature is
//! accepted.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumConditions, EquilibriumSolveOptions, PhaseEquilibriumPipelineRequest,
    SubstanceSystemSpecBuilder, SubstancesContainer, TemperatureGrid,
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

    let range = PhaseEquilibriumPipelineRequest::new(
        spec,
        vec![0.1, 0.05, 1.9],
        EquilibriumConditions::new(2_300.0, 101_325.0, 101_325.0)?,
    )
    .with_solve_options(EquilibriumSolveOptions::new().with_production_cascade())
    .solve_temperature_range(TemperatureGrid::new(vec![2_300.0, 2_500.0, 2_700.0])?)?;

    println!(
        "range: points={}, builds={}, reuses={}",
        range.report().point_count(),
        range.report().formulation_builds(),
        range.report().formulation_reuses(),
    );
    for point in range.points() {
        println!(
            "T = {:.0} K, moles = {:?}",
            point.report().temperature(),
            point.solution().component_moles()
        );
    }
    Ok(())
}
