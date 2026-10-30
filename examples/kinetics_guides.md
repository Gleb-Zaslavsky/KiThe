# Kinetics Example Guides

These examples cover the ordinary `src/Kinetics/` workflows and intentionally do not use
`src/Kinetics/experimental_kinetics/`.

Run them with:

```powershell
cargo run --example kinetics_guide_local_database_lookup
cargo run --example kinetics_guide_shortcuts_and_rate_constants
cargo run --example kinetics_guide_construct_mechanism
cargo run --example kinetics_guide_custom_reaction_payloads
cargo run --example kinetics_guide_direct_facade_and_symbolic
```

The examples demonstrate:

- local database lookup by library, reaction id, reagents, and products;
- shortcut selection such as `C1..C3`;
- mechanism construction from seed substances;
- numeric rate-constant calculations from database records or user-provided JSON payloads.
- direct facade parsing, stoichiometric analysis, and numeric plus symbolic constants without a local database.

## Audit before production use

Catalog audits are diagnostic and non-transactional: each requested reaction
gets its own result, and all failures are retained in the report. Audit the
payloads and equations first; only a clean selection should proceed through the
ordinary transactional `KinData` workflow. The audit report does not mutate or
populate `KinData`.

```rust,no_run
use KiThe::Kinetics::error::{KineticsError, KineticsResult};
use KiThe::Kinetics::library_audit::{
    audit_library_payload_range, audit_parsed_equations,
};
use KiThe::Kinetics::User_reactions::KinData;

fn audited_cantera_subset() -> KineticsResult<KinData> {
    let payload_report = audit_library_payload_range("Cantera", 1..=3)?;
    let equation_report = audit_parsed_equations(&payload_report);
    if !payload_report.is_clean() || !equation_report.is_clean() {
        payload_report.write_json("cantera-payload-audit.json")?;
        equation_report.write_json("cantera-equation-audit.json")?;
        return Err(KineticsError::InvalidReactionData(
            "selected reaction set failed audit".into(),
        ));
    }

    let mut kinetics = KinData::new();
    kinetics.set_reactions_from_shortcut_range("C1..3".into())?;
    kinetics.get_reactions_from_shortcuts()?;
    kinetics.reactdata_parsing()?;
    kinetics.analyze_reactions()?;
    Ok(kinetics)
}
```
