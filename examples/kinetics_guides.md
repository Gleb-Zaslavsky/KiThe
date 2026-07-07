# Kinetics Example Guides

These examples cover the ordinary `src/Kinetics/` workflows and intentionally do not use
`src/Kinetics/experimental_kinetics/`.

Run them with:

```powershell
cargo run --example kinetics_guide_local_database_lookup
cargo run --example kinetics_guide_shortcuts_and_rate_constants
cargo run --example kinetics_guide_construct_mechanism
cargo run --example kinetics_guide_custom_reaction_payloads
```

The examples demonstrate:

- local database lookup by library, reaction id, reagents, and products;
- shortcut selection such as `C1..C3`;
- mechanism construction from seed substances;
- numeric rate-constant calculations from database records or user-provided JSON payloads.
