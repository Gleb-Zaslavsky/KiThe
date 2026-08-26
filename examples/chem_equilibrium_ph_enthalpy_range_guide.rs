//! Guide: a fixed-pressure sweep over total-enthalpy targets.
//!
//! Each target is solved transactionally. Continuation carries only the prior
//! accepted composition and temperature, never an iterate from a failed solve.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumConditions, EquilibriumSolveOptions, MultiphaseEquilibriumLayout,
    MultiphaseInitialComposition, PhEnthalpyGrid, PhRangeRequest, ResolvedPhaseEquilibriumRequest,
    ResolvedThermochemistry, SubstanceSystemFactory, SubstanceSystemSpecBuilder,
    SubstancesContainer, TemperatureBounds, solve_resolved_pt,
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
    let resolved = SubstanceSystemFactory::resolve_phase_system(spec)?;
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
    let initial = MultiphaseInitialComposition::from_dense(&layout, vec![0.1, 0.05, 1.9])?;
    let options = EquilibriumSolveOptions::new().with_production_cascade();
    let pressure = 101_325.0;
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;

    // Use P,T reference points only to construct deterministic targets for
    // this guide. Real callers may provide measured total enthalpies instead.
    let mut targets = [2_300.0, 2_500.0, 2_700.0]
        .into_iter()
        .map(|temperature| {
            let point = solve_resolved_pt(
                ResolvedPhaseEquilibriumRequest::new(
                    &resolved,
                    EquilibriumConditions::new(temperature, pressure, pressure)?,
                    initial.clone(),
                )
                .with_solve_options(options.clone()),
            )?;
            Ok::<_, Box<dyn std::error::Error>>((
                thermochemistry
                    .enthalpy_model()
                    .evaluate_total(point.component_moles(), temperature)?,
                temperature,
            ))
        })
        .collect::<Result<Vec<_>, _>>()?;
    targets.sort_by(|left, right| left.0.total_cmp(&right.0));

    let range = PhRangeRequest::from_resolved_thermochemistry(
        &resolved,
        initial,
        pressure,
        pressure,
        PhEnthalpyGrid::new(targets.iter().map(|(enthalpy, _)| *enthalpy).collect())?,
        TemperatureBounds::new(2_100.0, 2_900.0)?,
        targets[0].1,
        thermochemistry,
    )?
    .with_solve_options(options)
    .solve()?;

    println!(
        "P,H range: points={}, continuation={}, builds={}, reuses={}",
        range.report().point_count(),
        range.report().continuation_points(),
        range.report().formulation_builds(),
        range.report().formulation_reuses(),
    );
    for point in range.points() {
        println!(
            "H = {:.6e} J, T = {:.6} K",
            point.report().target_enthalpy_joules(),
            point.report().solved_temperature()
        );
    }
    Ok(())
}
