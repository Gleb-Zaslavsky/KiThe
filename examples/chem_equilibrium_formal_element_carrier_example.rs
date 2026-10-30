//! Public-facade example: use formal carriers as input syntax only.
//!
//! Run with:
//! `cargo run --example chem_equilibrium_formal_element_carrier_example --no-default-features`
//!
//! A formal carrier such as `H2` is parsed into elemental amounts and then
//! disappears from the physical solve. It is not a thermochemical component,
//! a phase, or a numerical seed. This is useful when an application receives
//! compact formula-shaped feed data but must still solve a closed elemental
//! problem through the normal facade.

use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    ElementInventory, EquilibriumCalculator, EquilibriumCalculatorOutcome, FormalElementCarrier,
    PhaseEquilibriumInputKind,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Two formal carriers contribute H=2*2=4 and O=1*2=2 mol-atoms.
    // The formulas are input notation only; the resulting inventory retains
    // canonical element labels and stable numeric totals.
    let carriers = [
        FormalElementCarrier::new("H2", 2.0)?,
        FormalElementCarrier::new("O2", 1.0)?,
    ];
    let inventory = ElementInventory::from_formal_carriers(carriers)?;

    // This is the same application-facing facade used by the direct elemental
    // example. No formal carrier is added to the species list.
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
    assert_eq!(
        inventory.entries().collect::<Vec<_>>(),
        vec![("H", 4.0), ("O", 2.0)]
    );
    assert_eq!(report.element_totals(), [4.0, 2.0]);

    println!("formal carriers: H2(2.0), O2(1.0)");
    println!(
        "canonical inventory: {:?}",
        inventory.entries().collect::<Vec<_>>()
    );
    println!("published input kind: {:?}", report.input_kind());
    println!(
        "equilibrium moles: {:?}",
        point.solution().component_moles()
    );

    Ok(())
}
