//! Guide: run direct kinetic payloads through the KinData facade.
//!
//! This example needs no local mechanism database. It parses user-supplied
//! Arrhenius records, analyzes stoichiometry, evaluates numeric rate constants,
//! and prepares symbolic constants with a scaled temperature variable.
//!
//! Run with:
//! `cargo run --example kinetics_guide_direct_facade_and_symbolic`

use KiThe::Kinetics::User_reactions::KinData;
use KiThe::Kinetics::error::KineticsResult;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use serde_json::json;

fn main() -> KineticsResult<()> {
    let mut kinetics = KinData::new();
    kinetics.append_reaction(vec![
        json!({
            "type": "elem",
            "eq": "2H2 + O2 -> 2H2O",
            "Arrenius": [2.0e6, 0.5, 1200.0]
        }),
        json!({
            "type": "elem",
            "eq": "H2O -> H2 + 0.5O2",
            "Arrenius": [7.0e5, 0.0, 800.0]
        }),
    ])?;

    // The facade parses payloads and builds stoichiometric and elemental views.
    kinetics.kinetic_main()?;
    println!("Species: {:?}", kinetics.substances());
    println!(
        "Stoichiometry: {:?}",
        kinetics.stoichiometric_analyzer().stecheo_matrx
    );

    let temperature_k = 900.0;
    for (index, equation) in kinetics.equations().iter().enumerate() {
        let rate = kinetics.calc_K_const_for_1_react(index, temperature_k, None, None)?;
        println!("k({temperature_k:.0} K) = {rate:.6e} for {equation}");
    }

    // Elementary reactions do not need pressure or collider concentrations.
    // T_scaled is substituted for the symbolic temperature T in every expression.
    kinetics.calc_sym_constants(None, None, Some(Expr::Var("T_scaled".to_string())))?;
    for (equation, expression) in kinetics
        .equations()
        .iter()
        .zip(kinetics.symbolic_constants().unwrap())
    {
        println!("k(T_scaled) for {equation}: {expression:?}");
    }

    Ok(())
}
