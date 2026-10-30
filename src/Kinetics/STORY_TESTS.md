# Kinetics Story Tests

These are correctness stories, not performance baselines. Commands use the
debug test profile with one Cargo build job and one test thread. The Kinetics
campaign includes ignored stories and excludes `experimental_kinetics`.

## Story KIN-12-01: Direct facade prepares numeric and symbolic kinetics

**Name:** Direct reaction payloads retain coherent stoichiometry and rate views

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features kinetics_facade_direct_payload_symbolic_and_stoichiometry_story -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can a caller provide reaction payloads without a
database, process them through `KinData`, and obtain consistent stoichiometric,
numeric-rate, and symbolic-rate views?

**Fixture and units:** Two elementary reactions, `2H2 + O2 -> 2H2O` and
`H2O -> H2 + 0.5O2`, with Arrhenius vectors `[A, n, E]`. Temperature is `900 K`;
activation energies use the crate's existing Arrhenius convention and `R =
8.314` in the matching units.

**Acceptance:** The facade discovers three species; the stoichiometric matrix
is exactly `[[-2, -1, 2], [1, 0.5, -1]]`; the first numeric rate agrees with the
independent Arrhenius expression to relative `1e-12`; and both reactions have
published symbolic constants in the normalized view after temperature scaling.

**Observed result (2026-10-26):** 1 passed, 0 failed.

**Conclusion and limits:** This verifies the direct payload -> parse -> analyze
-> numeric/symbolic facade path for elementary reactions. It does not validate
database provenance, falloff/pressure models, or external experimental rates.

## Story KIN-12-02: Batch failures identify the failing reaction without partial publication

**Name:** Malformed reaction batches fail at the exact row and retain a coherent stage

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features kinetics_batch_parse_failures_are_located_and_do_not_publish_partial_results_story -- --nocapture --test-threads=1
```

**Question / hypothesis:** If a payload or equation in the middle of a reaction
batch is malformed, does the error identify its position and does `KinData`
avoid publishing a partial parsed or analyzed result?

**Fixture:** Three elementary payloads keyed to shortcuts `C_10`, `C_11`, and
`C_12`, with the invalid record at zero-based index 1. The first case has an
invalid Arrhenius vector; the second has a
malformed equation (`O2 + H2`, no reaction arrow) with otherwise valid kinetic data.

**Acceptance:** Payload parsing reports shortcut `C_11` and index 1, leaving the
full raw batch in `ShortcutsSelected`, with no typed equations or matrix.
Equation analysis reports shortcut `C_11`, row 1, and the offending equation;
the complete three-record typed batch remains at `ReactionDataParsed`, while no stoichiometric matrix is published.
The error remains a typed `KineticsError` with its original cause.

**Observed result (2026-10-26):** 1 passed, 0 failed.

**Conclusion and limits:** Failure is fail-fast and stage-transactional, not a
best-effort batch. A complete parsed stage is retained when a later analysis
stage fails; it is not a partially parsed list. The test supplies shortcut
identities directly so diagnostics can be checked without depending on mutable
bundled catalog contents.

## Story KIN-11-01: Parse the complete Cantera library

**Name:** Every bundled Cantera payload deserializes and validates as typed reaction data

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features story_cantera_full_library_deserializes_into_typed_reactions -- --include-ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Does the complete bundled Cantera reaction library,
not just a hand-picked sample, satisfy the current `ReactionData` schema and
validation contract?

**Fixture:** The checked-in `Cantera` library, all records in deterministic
reaction-ID order. Parsing uses the same `parse_kinetic_data_vec` function as
the facade, without stoichiometric interpretation or rate evaluation.

**Acceptance:** Catalog loading succeeds; every library record produces one
validated `ReactionData` and one equation; counts match exactly. The test prints
the number of records parsed.

**Observed result (2026-10-26, debug):** 324 payloads parsed into 324 typed
reactions and 324 equations; 1 passed, 0 failed.

**Conclusion and limits:** This is a full-catalog schema/parameter stress gate.
It does not prove that every equation is stoichiometrically valid or that every
rate model has appropriate experimental provenance.

## Story KIN-09-01: Recursive Cantera mechanism closure

**Name:** Seeded Cantera mechanism construction returns a coherent recursive closure

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features story_cantera_recursive_mechanism_construction_reaches_a_consistent_closure -- --include-ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Does recursive directed closure from a broad set of
common seed species in the Cantera catalog terminate and return aligned
reaction IDs, equations, and typed records while retaining the seeds?

**Fixture:** The bundled Cantera catalog and seeds `H`, `H2`, `O`, `O2`, `OH`,
`H2O`, `CO`, `CO2`, and `CH4`.

**Acceptance:** The closure is nonempty, retains every seed, and has identical
counts for reaction IDs, equations, and typed `ReactionData`; each typed record
equation matches its aligned equation entry. The test prints closure sizes.

**Observed result (2026-10-26, debug):** 212 reactions, 34 discovered species,
9 seeds; 1 passed, 0 failed (closure call 0.15 s in this run).

**Conclusion and limits:** This stresses recursive mechanism construction
against the checked-in Cantera catalog. It is correctness evidence for this seed
set, not a guarantee for every catalog or arbitrary mechanism size.

## Story KIN-13-01: Non-transactional catalog and equation audit

**Name:** Bounded library audit retains successful records and localizes failures

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features story_cantera_range_audit_writes_payload_and_equation_reports -- --include-ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Can a selected library range be audited without
fail-fast behavior, with payload deserialization/validation and equation
analysis reported as independent stages?

**Fixture:** Cantera reaction IDs `0..=100`. Each ID receives an individual
payload outcome. Equation analysis runs only for successful typed payloads and
also records each result independently.

**Acceptance:** The payload report contains exactly 101 entries, including
missing IDs; the equation report contains one entry for every successful
payload. Both reports are written as JSON and read back successfully. Failures
retain reaction ID, stage, typed category, message, and raw payload when one
was present.

**Observed result (2026-10-26, debug):** 101 requested payloads parsed and
validated; 101 equations analyzed; both JSON reports round-tripped; 1 passed,
0 failed.

**Conclusion and limits:** This diagnostic lifecycle is deliberately separate
from transactional `KinData` production methods. Its partial successes are
audit findings and must not be published as a production mechanism.

**Related catalog finding:** A broader NUIG seed probe reaches a catalog record
whose Troe payload fails `ReactionData::validate` (`C2H4+H(+M)<=>C2H5(+M)`,
reported at selected reaction index 223). That failure is not suppressed or
included in this Cantera closure gate; repairing the bundled NUIG data is a
separate correctness task.

## Story KIN-04-01: Stoichiometry parser and composition cache contract

**Name:** Supported notation is preserved and species changes invalidate derived composition

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features Kinetics::stoichiometry_analyzer::tests -- --nocapture --test-threads=1
```

**Question / hypothesis:** Do supported phase/charge/formula-parenthesis and
fractional/scientific coefficient forms parse without losing species identity,
and do elemental/molar-mass views follow species and group revisions?

**Fixture and units:** Parser fixtures include `H+`, `OH-(aq)`, `Fe^3+`,
`SO4^2-`, `Ca(OH)2(s)`, `0.5` and `1e-3` coefficients. Formula composition
projection removes supported terminal phase/charge notation while preserving
the kinetic species names. Rational text `1/2` is intentionally unsupported.

**Acceptance:** All 16 analyzer tests pass, including malformed-input rejection,
cache invalidation after species/group changes, and preservation of a matching
explicit molar-mass override.

**Observed result (2026-10-26):** 16 passed, 0 failed.

**Conclusion and limits:** The story locks down the current parser contract and
invalidation behavior. It does not claim arbitrary chemical notation support
or universal phase-tag recognition.
