//! Guide: one fixed-pressure, fixed-total-enthalpy equilibrium point.
//!
//! A known `P,T` state supplies a reproducible enthalpy target. The subsequent
//! `P,H` solve recovers temperature and equilibrium composition through the
//! same resolved local thermochemistry.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumConditions, EquilibriumConstraint, EquilibriumSolveOptions,
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition, PhSolveMode,
    ResolvedPhaseEnthalpyRequest, ResolvedPhaseEquilibriumRequest, ResolvedThermochemistry,
    SubstanceSystemFactory, SubstanceSystemSpecBuilder, SubstancesContainer, TemperatureBounds,
    TotalEnthalpyJoules, solve_resolved_ph, solve_resolved_pt,
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
    let reference_temperature = 2_500.0;

    let reference = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(reference_temperature, pressure, pressure)?,
            initial.clone(),
        )
        .with_solve_options(options.clone()),
    )?;
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(reference.component_moles(), reference_temperature)?;

    let solution = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial,
            EquilibriumConstraint::ph_joules(
                pressure,
                pressure,
                TotalEnthalpyJoules::new(target_enthalpy)?,
                reference_temperature,
            )?,
            TemperatureBounds::new(2_100.0, 2_900.0)?,
            thermochemistry,
        )?
        .with_ph_solve_mode(PhSolveMode::Auto)
        .with_solve_options(options),
    )?;

    println!("P,H temperature = {:.6} K", solution.temperature());
    println!(
        "P,H equilibrium moles = {:?}",
        solution.equilibrium().component_moles()
    );
    println!(
        "P,H path = {:?}, enthalpy error = {:.3e} J",
        solution.report().solve_path(),
        solution.enthalpy_error(),
    );
    Ok(())
}
