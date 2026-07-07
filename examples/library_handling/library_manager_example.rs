use KiThe::library_manager::{ with_library_manager, with_library_manager_mut};
use KiThe::settings::Settings;
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Library Manager Example ===\n");

    // 1. Basic usage - get current library paths
    println!("1. Current library configuration:");
    with_library_manager(|manager| {
        println!("  Substance Base: {}", manager.substance_base_path());
        println!(
            "  All Keys Substance: {}",
            manager.all_keys_substance_path()
        );
        println!("  Reactbase: {}", manager.reactbase_path());
        println!("  Dict Reaction: {}\n", manager.dict_reaction_path());
    });

    // 2. Using Settings for user-friendly interface
    println!("2. Using Settings interface:");
    let mut settings = Settings::new();

    println!("Available libraries:");
    for lib in settings.get_available_libraries() {
        println!(
            "  - {}: {}",
            lib,
            settings.get_library_version(lib).unwrap()
        );
    }
    println!();

    // 3. Programmatic library version change (for library use)
    println!("3. Programmatic usage example:");
    with_library_manager_mut(|manager_mut| {
        // Example: Change to version 3 files (if they exist)
        if std::path::Path::new("substance_base_v3.json").exists() {
            match manager_mut.set_substance_base("substance_base_v3.json") {
                Ok(_) => println!("  ✓ Successfully updated substance base to v3"),
                Err(e) => println!("  ✗ Failed to update: {}", e),
            }
        } else {
            println!("  ℹ substance_base_v3.json not found, skipping update");
        }
    });

    // 4. Batch update example
    println!("\n4. Batch update example:");
    let mut updates = HashMap::new();

    // Only update if files exist
    if std::path::Path::new("substance_base_v2.json").exists() {
        updates.insert("substance_base", "substance_base_v2.json");
    }
    if std::path::Path::new("all_keys_substance.json").exists() {
        updates.insert("all_keys_substance", "all_keys_substance.json");
    }

    if !updates.is_empty() {
        with_library_manager_mut(|manager_mut| match manager_mut.update_libraries(updates) {
            Ok(_) => println!("  ✓ Batch update successful"),
            Err(e) => println!("  ✗ Batch update failed: {}", e),
        });
    }

    // 5. Settings-based update (for GUI use)
    println!("\n5. Settings-based update (GUI-friendly):");
    if std::path::Path::new("substance_base_v2.json").exists() {
        match settings.set_library_version("Substance Base", "substance_base_v2.json") {
            Ok(_) => println!("  ✓ Updated Substance Base via Settings"),
            Err(e) => println!("  ✗ Failed to update via Settings: {}", e),
        }
    }

    // 6. Reset to defaults
    println!("\n6. Reset to defaults:");
    match settings.reset_to_defaults() {
        Ok(_) => println!("  ✓ Reset to default configuration"),
        Err(e) => println!("  ✗ Failed to reset: {}", e),
    }

    println!("\n=== Example completed ===");
    Ok(())
}
