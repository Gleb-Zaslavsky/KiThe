//! Guide: configure thermo lookup and use the canonical read-only accessors.
//!
//! This example focuses on the stable, user-facing entry points:
//! - typed library identifiers
//! - explicit search instructions
//! - canonical read-only views over the configured query
//!
//! Run with:
//! `cargo run --example thermodynamics_guide_canonical_accessors`

use KiThe::Thermodynamics::User_substances::{LibraryPriority, SubsData};
use KiThe::Thermodynamics::thermo_lib_api::LibraryId;

fn main() {
    let mut subs_data = SubsData::new();

    // Configure a small query set.
    subs_data.set_substances(vec!["H2O".to_string(), "CO2".to_string()]);

    // Prefer canonical library identifiers instead of ad-hoc strings.
    subs_data.set_library_priority_id(LibraryId::NasaGas, LibraryPriority::Priority);
    subs_data.set_explicit_search_instruction("H2O".to_string(), LibraryId::NasaGas);

    // Inspect the normalized, read-only configuration.
    println!("Substances: {:?}", subs_data.substances());
    println!("Library priorities: {:?}", subs_data.library_priorities());
    println!("Explicit search map: {:?}", subs_data.explicit_search_map());

    // Molar masses are exposed through a clean accessor instead of a public typo field.
    println!("H2O molar mass: {:?}", subs_data.molar_mass_of("H2O"));
}
