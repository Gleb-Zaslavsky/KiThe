//! Public-facade example: solve a fixed-pressure, fixed-temperature problem
//! from a physical elemental inventory.
//!
//! Run with:
//! `cargo run --example chem_equilibrium_element_inventory_pt_example --no-default-features`
//!
//! The important distinction is that the inventory is not a synthetic species
//! and is not passed to a hand-built reaction solver. `ElementInventory` is the
//! physical conservation vector `b`; `EquilibriumCalculator` resolves the real
//! H2/O2/H2O universe and runs the canonical production workflow.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    ElementInventory, EquilibriumCalculator, EquilibriumCalculatorOutcome,
    PhaseEquilibriumInputKind,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // These are mol-atom amounts, not molecule amounts. They describe the
    // closed system H=4 mol-atoms and O=2 mol-atoms. The solver is free to
    // redistribute those atoms among the declared molecular components.
    let inventory = ElementInventory::from_amounts([("H", 4.0), ("O", 2.0)])?;

    // For a difficult multiphase basin, the same builder also accepts an
    // optional `element_numerical_seed(...)`. That value only changes where
    // nonlinear solving starts; `inventory` remains the physical `b` and the
    // published report remains `ElementInventory`.
    //
    // The facade owns phase/layout construction, thermochemistry lookup,
    // input validation, nonlinear fallback, and accepted-state publication.
    // Applications should stop at this boundary instead of assembling those
    // internal objects themselves.
    let outcome = EquilibriumCalculator::builder()
        .single_ideal_gas(["H2", "O2", "H2O"])?
        .element_inventory(inventory.clone())
        .prefer_libraries(["NASA_gas"])
        .offline_only()
        .pressure_pa(101_325.0)
        .reference_pressure_pa(101_325.0)
        .at_temperature(2_500.0)
        .production_cascade()
        .solve()?;

    let EquilibriumCalculatorOutcome::PtPoint(point) = outcome else {
        unreachable!("at_temperature selects the typed P,T point outcome");
    };

    let report = point.solution().build_report();
    assert_eq!(
        report.input_kind(),
        PhaseEquilibriumInputKind::ElementInventory
    );
    assert_eq!(report.element_labels(), ["H", "O"]);
    assert_eq!(report.element_totals(), [4.0, 2.0]);

    println!("inventory: {:?}", inventory.entries().collect::<Vec<_>>());
    println!("input kind: {:?}", report.input_kind());
    println!("components: {:?}", point.solution().metadata().components());
    println!(
        "equilibrium moles: {:?}",
        point.solution().component_moles()
    );
    println!("element totals: {:?}", report.element_totals());
    println!(
        "validation: {:?}",
        point.solution().accepted_solution().validation()
    );

    Ok(())
}
