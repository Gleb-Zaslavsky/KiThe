//! Public-facade example: solve a fixed-pressure, fixed-total-enthalpy
//! problem from an elemental inventory.
//!
//! Run with:
//! `cargo run --example chem_equilibrium_element_inventory_ph_example --no-default-features`
//!
//! The P,H facade treats total enthalpy as an extensive target. Scaling the
//! inventory and the target together should preserve the recovered temperature
//! and normalized composition; the frozen TP-1906 stories exercise this same
//! contract on the larger CHON + graphite fixture.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    ElementInventory, EquilibriumCalculator, EquilibriumCalculatorOutcome, PhSolveMode,
    PhaseEquilibriumInputKind, TemperatureBounds, TotalEnthalpyJoules,
};

fn solve(scale: f64) -> Result<f64, Box<dyn std::error::Error>> {
    // H and O are physical mol-atom totals. Both the inventory and enthalpy
    // target use the same scale so that the thermodynamic state is comparable.
    let inventory = ElementInventory::from_amounts([("H", 4.0 * scale), ("O", 2.0 * scale)])?;
    let target_enthalpy = TotalEnthalpyJoules::new(-72_679.18 * scale)?;

    let outcome = EquilibriumCalculator::builder()
        .single_ideal_gas(["H2", "O2", "H2O"])?
        .element_inventory(inventory)
        .prefer_libraries(["NASA_gas"])
        .offline_only()
        .pressure_pa(1_000_000.0)
        .reference_pressure_pa(100_000.0)
        .at_total_enthalpy(
            target_enthalpy,
            2_000.0,
            TemperatureBounds::new(200.0, 6_000.0)?,
        )
        // The broad interval crosses native NASA polynomial boundaries. The
        // nested temperature route handles that piecewise thermochemistry;
        // symbolic P,H is intended for a single native interval.
        .ph_solve_mode(PhSolveMode::NestedTemperature)
        .production_cascade()
        .solve()?;

    let EquilibriumCalculatorOutcome::PhPoint(point) = outcome else {
        unreachable!("at_total_enthalpy selects the typed P,H point outcome");
    };
    let report = point.solution().equilibrium().build_report();
    assert_eq!(
        report.input_kind(),
        PhaseEquilibriumInputKind::ElementInventory
    );
    assert!(point.solution().enthalpy_error().is_finite());
    println!(
        "scale={scale:.1e} temperature={:.6} K enthalpy_error={:.3e} J",
        point.solution().temperature(),
        point.solution().enthalpy_error(),
    );
    Ok(point.solution().temperature())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let base_temperature = solve(1.0)?;
    let scaled_temperature = solve(10.0)?;

    // Temperature is intensive, so the two facade calls should recover the
    // same value even though the physical inventory is ten times larger.
    // The two solves may choose slightly different nonlinear paths, so retain
    // a tight numerical envelope instead of requiring bit-for-bit equality.
    assert!((scaled_temperature - base_temperature).abs() < 1.0e-4);
    println!("intensive temperature invariant: OK");
    Ok(())
}
