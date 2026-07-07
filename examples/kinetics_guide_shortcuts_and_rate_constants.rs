//! Guide: select database reactions by shortcuts and calculate rate constants.
//!
//! Run with:
//! `cargo run --example kinetics_guide_shortcuts_and_rate_constants`

use KiThe::Kinetics::User_reactions::KinData;
use KiThe::Kinetics::error::KineticsResult;

fn main() -> KineticsResult<()> {
    // C1..C3 means Cantera reactions 1, 2, and 3.
    // At this stage the shortcuts are only a user assignment, not materialized reaction data.
    let mut kinetics = KinData::from_shortcut_range("C1..C3".to_string())?;
    println!("Selected shortcuts: {:?}", kinetics.shortcut_names());

    // Resolve shortcut ids into local database records, then parse and analyze the equations.
    kinetics.get_reactions_from_shortcuts()?;
    kinetics.kinetic_main()?;

    println!("Resolved equations:");
    for equation in kinetics.equations() {
        println!("  {equation}");
    }

    // Calculate all rate constants at one temperature.
    // `save_rearranged = Some(false)` keeps the KinData order unchanged.
    let temperature_k = 1000.0;
    let rate_constants =
        kinetics.calc_K_const_for_all_reactions(temperature_k, None, None, Some(false))?;
    for (equation, k) in kinetics.equations().iter().zip(rate_constants.iter()) {
        println!("k({temperature_k:.0} K) = {k:.6e} for {equation}");
    }

    // A temperature sweep returns the maximum rate constant found on the grid for each reaction.
    let max_rate_constants = kinetics.calc_K_const_for_all_reactions_forTrange(
        800.0,
        1200.0,
        8,
        None,
        None,
        Some(false),
    )?;
    println!("Max k values on the 800..1200 K grid:");
    for (index, k) in max_rate_constants.iter().enumerate() {
        println!("  reaction #{index}: max k = {k:.6e}");
    }

    Ok(())
}
