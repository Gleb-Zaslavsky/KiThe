//! Guide: calculate numeric rate constants from user-supplied reaction payloads.
//!
//! This path does not require a local kinetics database. It is useful when the user already
//! has kinetic parameters in the same JSON shape as the bundled libraries.
//!
//! Run with:
//! `cargo run --example kinetics_guide_custom_reaction_payloads`

use KiThe::Kinetics::User_reactions::KinData;
use KiThe::Kinetics::error::KineticsResult;
use serde_json::json;

fn main() -> KineticsResult<()> {
    let mut kinetics = KinData::new();

    // Minimal elementary reactions: type, equation, and Arrhenius parameters [A, n, E].
    // The activation energy convention is the same as in the bundled database records.
    kinetics.append_reaction(vec![
        json!({
            "type": "elem",
            "eq": "A -> B",
            "Arrenius": [1.0e12, 0.0, 20_000.0]
        }),
        json!({
            "type": "elem",
            "eq": "B -> C",
            "Arrenius": [1.0e13, 0.0, 25_000.0]
        }),
    ])?;

    // Parse JSON payloads into typed ReactionData. Full stoichiometric analysis is optional
    // for this numeric rate-constant calculation.
    kinetics.reactdata_parsing()?;

    let temperature_k = 1000.0;
    let rate_constants =
        kinetics.calc_K_const_for_all_reactions(temperature_k, None, None, Some(false))?;

    for (equation, k) in kinetics.equations().iter().zip(rate_constants.iter()) {
        println!("k({temperature_k:.0} K) = {k:.6e} for {equation}");
    }

    Ok(())
}
