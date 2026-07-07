//! Guide: construct a reaction mechanism from seed substances.
//!
//! Run with:
//! `cargo run --example kinetics_guide_construct_mechanism`

use KiThe::Kinetics::User_reactions::KinData;
use KiThe::Kinetics::error::KineticsResult;

fn main() -> KineticsResult<()> {
    // The mechanism builder searches the chosen local library for reactions connected to
    // the seed substances and their discovered products.
    let seed_substances = vec!["O".to_string(), "NH3".to_string(), "NO".to_string()];
    let library = "NUIG".to_string();

    let mut kinetics = KinData::new();
    kinetics.construct_mechanism(seed_substances, library)?;

    // `construct_mechanism` already materializes parsed ReactionData; `kinetic_main` completes
    // stoichiometric analysis for the generated mechanism.
    kinetics.kinetic_main()?;

    println!("Mechanism contains {} reactions", kinetics.reaction_count());
    println!("Discovered substances: {:?}", kinetics.substances());
    println!("First few mechanism equations:");
    for equation in kinetics.equations().iter().take(10) {
        println!("  {equation}");
    }

    Ok(())
}
