//! Guide: concise application-facing P,H equilibrium request.
//!
//! This is the public calculator facade. For the explicit resolved-data and
//! request assembly underneath it, see `chem_equilibrium_reactive_gas_ph_example`.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumCalculator, EquilibriumCalculatorOutcome, EquilibriumTimingMode, PhSolveMode,
    TemperatureBounds, TotalEnthalpyJoules, format_ph_solution_execution_summary,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // In an application this extensive value normally comes from a measured
    // feed or an upstream energy model, in J on the same inventory basis.
    let measured_total_enthalpy_j = TotalEnthalpyJoules::new(-7.267_918e4)?;

    let outcome = EquilibriumCalculator::builder()
        .single_ideal_gas([
            "CO", "H2", "O2", "N2", "Ar", "CO2", "H2O", "H", "O", "OH", "NO", "NO2", "NH3", "CH4",
        ])?
        .initial_moles(vec![
            1.0, 2.0, 1.0, 3.76, 0.044, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ])
        .prefer_libraries(["NASA_gas"])
        .offline_only()
        .pressure_pa(1_000_000.0)
        .reference_pressure_pa(100_000.0)
        .at_total_enthalpy(
            measured_total_enthalpy_j,
            2_000.0,
            TemperatureBounds::new(200.0, 6_000.0)?,
        )
        // The selected NASA records span native polynomial intervals, so the
        // nested P,H route is the compatible formulation for this wide range.
        .ph_solve_mode(PhSolveMode::NestedTemperature)
        .production_cascade()
        .timing(EquilibriumTimingMode::Enabled)
        .solve()?;

    let EquilibriumCalculatorOutcome::PhPoint(point) = outcome else {
        unreachable!("builder mode determines the matching typed outcome");
    };
    println!("{}", format_ph_solution_execution_summary(point.solution()));
    println!(
        "lookup: phases={}, NIST fallback enabled={}",
        point.resolved().phase_specs().len(),
        point.resolved().report().nist_fallback_enabled(),
    );
    Ok(())
}
