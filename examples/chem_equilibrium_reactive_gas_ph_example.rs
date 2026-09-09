//! Guide: a realistic fixed-pressure, fixed-enthalpy reactive-gas calculation.
//!
//! The feed is a simple CO/H2 synthesis-gas stream diluted with air. The
//! candidate universe contains stable products, radicals, and nitrogen oxides;
//! it is intentionally wider than the feed. All thermochemistry is resolved
//! from the bundled local NASA gas library before the P,H solve begins.
//!
//! For a runnable and deterministic guide, the target enthalpy is generated
//! from a known P,T reference state. A real caller normally replaces that one
//! reference block with a measured or process-model total enthalpy in J.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumConditions, EquilibriumConstraint, EquilibriumPresentationReport,
    EquilibriumSolveOptions, MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
    PhSolveMode, ResolvedPhaseEnthalpyRequest, ResolvedPhaseEquilibriumRequest,
    ResolvedThermochemistry, SubstanceSystemFactory, SubstanceSystemSpecBuilder,
    SubstancesContainer, TemperatureBounds, TotalEnthalpyJoules,
    format_ph_solution_execution_summary, solve_resolved_ph, solve_resolved_pt,
};

const PRESSURE_PA: f64 = 1_000_000.0;
const GAS_STANDARD_STATE_PRESSURE_PA: f64 = 100_000.0;
const REFERENCE_TEMPERATURE_K: f64 = 2_400.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // The order here is also the dense-inventory order below. Species omitted
    // from the feed start at zero and receive only the configured trace seed.
    let candidates = [
        "CO", "H2", "O2", "N2", "Ar", "CO2", "H2O", "H", "O", "OH", "NO", "NO2", "NH3", "CH4",
    ];
    let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(
        candidates.iter().map(|name| (*name).to_owned()).collect(),
    ))
    .with_library_priorities(vec!["NASA_gas".to_owned()])
    .with_search_in_nist(false)
    .build()?;

    let resolved = SubstanceSystemFactory::resolve_phase_system(spec)?;
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
    let initial = MultiphaseInitialComposition::from_dense(
        &layout,
        vec![
            1.0,   // CO
            2.0,   // H2
            1.0,   // O2
            3.76,  // N2 from air
            0.044, // Ar from air
            0.001, // CO2 impurity in air/feed
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ],
    )?;
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
    let bounds = thermochemistry.temperature_bounds();
    let options = EquilibriumSolveOptions::new().with_production_cascade();

    // In production, this is a measured total enthalpy of the supplied feed.
    // The reference solve only makes this example self-contained and repeatable.
    let reference = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(
                REFERENCE_TEMPERATURE_K,
                PRESSURE_PA,
                GAS_STANDARD_STATE_PRESSURE_PA,
            )?,
            initial.clone(),
        )
        .with_solve_options(options.clone()),
    )?;
    let target_enthalpy_j = thermochemistry
        .enthalpy_model()
        .evaluate_total(reference.component_moles(), REFERENCE_TEMPERATURE_K)?;

    let solution = solve_resolved_ph(
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            &resolved,
            initial,
            EquilibriumConstraint::ph_joules(
                PRESSURE_PA,
                GAS_STANDARD_STATE_PRESSURE_PA,
                TotalEnthalpyJoules::new(target_enthalpy_j)?,
                2_000.0,
            )?,
            TemperatureBounds::new(bounds.lower(), bounds.upper())?,
            thermochemistry,
        )?
        // This wide NASA candidate universe crosses native polynomial
        // intervals. Nested P,H remains the canonical compatible route;
        // monolithic symbolic P,H is intentionally reserved for a common
        // single-coefficient interval.
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
        .with_solve_options(options),
    )?;

    println!(
        "reactive-gas P,H: target_H={target_enthalpy_j:.6e} J, solved_T={:.6} K",
        solution.temperature(),
    );
    println!("{}", format_ph_solution_execution_summary(&solution));
    println!(
        "{}",
        EquilibriumPresentationReport::from_solution(solution.equilibrium()).render_compact()
    );
    Ok(())
}
