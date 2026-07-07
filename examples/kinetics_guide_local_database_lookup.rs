//! Guide: look up reactions in a local kinetics database.
//!
//! Run with:
//! `cargo run --example kinetics_guide_local_database_lookup`

use KiThe::Kinetics::error::KineticsResult;
use KiThe::Kinetics::kinetics_lib_api::KineticData;

fn main() -> KineticsResult<()> {
    // Load one of the bundled local reaction libraries.
    // The JSON files are cached by the library API, so repeated loads are cheap.
    let mut database = KineticData::new();
    database.open_json_files("Cantera")?;

    // Build equation/id indexes for the loaded library.
    database.print_all_reactions()?;
    println!("Loaded {} Cantera reactions", database.AllEquations.len());

    // Direct ID lookup returns the raw serde_json payload for that reaction.
    let reaction_id = "1";
    let reaction = database.search_reactdata_by_reaction_id(reaction_id)?;
    println!("Reaction {reaction_id}: {}", reaction["eq"]);

    // Product/reagent search works against the precomputed local relationship map.
    database.search_reaction_by_reagents_and_products(vec!["H".to_string(), "O".to_string()])?;
    println!(
        "Reactions where the task substances can be reagents: {:?}",
        database.FoundReactionsByReagents
    );
    println!(
        "Reactions where the task substances can be products: {:?}",
        database.FoundReactionsByProducts
    );

    Ok(())
}
