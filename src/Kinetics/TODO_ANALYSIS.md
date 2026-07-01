# Kinetics TODO analysis

Scope: `src/Kinetics/`, excluding `src/Kinetics/experimental_kinetics/`.

This file is an analytical TODO, not an implementation plan locked in stone. It records current risks in correctness, API ergonomics, test coverage, and solver integration.

## P0: correctness and panic-safety

- Replace user-facing `panic!`, `assert!`, `unwrap`, and `expect` paths with typed errors.
  Main hot spots: `User_reactions.rs`, `mechfinder_api.rs`, `kinetics_lib_api.rs`, `stoichiometry_analyzer.rs`, `molmass.rs`, `solid_state_kinetics_IVP.rs`.
  Public API should prefer `Result<_, KineticsError>` for invalid shortcuts, missing libraries/reactions, malformed equations, missing pressure/concentration data, failed file IO, JSON parse errors, and solver failures.

- Normalize `KinData` state.
  `KinData` is the top-level facade, but it currently stores reaction state in many parallel optional vectors: shortcuts, maps, equation strings, parsed values, reaction data, symbolic constants, stoichiometry, and `every_reaction`. This makes it easy to create inconsistent states. `every_reaction` looks like the right future single source of truth, but it is not yet the real backbone.

- Fix direct-reaction workflow.
  `set_reactions_directly` stores raw equations in `shortcut_reactions`, then `kinetic_main` routes through shortcut parsing. That contradicts the documented idea of directly supplied reactions and can leave `vec_of_reaction_Values` empty while later steps still continue.

- [x] Fix document generation with library maps.
  `create_kinetics_document` enumerates `map_of_reactions` by library and indexes `vec_of_reaction_data` with that library index. If a library contains several reactions, the same `ReactionData` can be paired with multiple equations or indexing can drift. This is a real data corruption/panic risk.

- [x] Fix falloff metadata loss.
  In `mechfinder_api/kinetics.rs`, `FalloffStruct::new` accepts `troe`, but initializes `troe: None`. Troe/SRI metadata is silently discarded.

- [x] Fix symbolic/numeric mismatch for three-body reactions.
  `ThreeBodyStruct::K_const` uses `exp(-E / RT)`, while symbolic `K_expr` builds `Exp(E / RT)`. Symbolic output can therefore disagree with numeric kinetics.

- Fix pressure-dependent interpolation.
  `PressureStruct::K_const` builds pressure arrays from `HashMap` key iteration, so ordering is not guaranteed before binary search. It also uses suspect low/high indexing in the interpolation branch and can call `location.unwrap()` after an `Err` result for out-of-range pressure.

- Revisit stoichiometric power parsing.
  In `stoichiometry_analyzer.rs`, the branch that should probably copy a stoichiometric coefficient into the kinetic power currently uses comparison (`power_coeff == stec_coeff`) instead of assignment. Tests currently tolerate powers of `1.0` for explicit stoichiometric coefficients, so the intended chemistry contract needs to be clarified and then tested.

- Make compute methods pure by default.
  Methods such as `analyze_reactions`, `calc_K_const_for_all_reactions`, `calc_sym_constants`, molar-mass parsing, and stoichiometry parsing print during normal computation. Computation should return data; formatting/logging should be opt-in.

## P1: LSODE2 support in solid-state IVP

- Add an idiomatic LSODE2 entry point.
  `solid_state_kinetics_IVP.rs` accepts arbitrary `SolverType` in `KineticModelIVP::new`, but there is no convenient constructor or shortcut for the now-standard `SolverType::LSODE2`. Candidate API:
  - `KineticModelIVP::new_lsode2(...)`
  - `KineticModelIVP::with_lsode2_defaults(...)`
  - or a solver-profile builder that makes LSODE2 the ergonomic default while still allowing RK/BDF/Radau selection.

- Decide whether LSODE2 should become the default solver.
  If LSODE2 is the de facto default for RustedSciThe IVP work, `KineticModelIVP::new` should either default to it through a high-level constructor or the docs should explain when to choose LSODE2 versus the existing methods.

- Add focused LSODE2 tests.
  Required coverage:
  - smoke solve for a simple solid-state model, e.g. `F1`, using `SolverType::LSODE2`;
  - comparison against an existing trusted solver on the same problem within domain-appropriate tolerance;
  - a stiffer/high-activation-energy case where LSODE2 is expected to be useful;
  - result shape invariants: `time`, `temperature`, `conversion`, `conversion_rate`, `rate_constant` lengths match;
  - no-solve and invalid-parameter cases return errors instead of panicking;
  - solver status/statistics are accessible or at least solver failure is not hidden.

- Check solver error propagation.
  `KineticModelIVP::solve` currently calls the universal solver and then returns `Ok(())`. Confirm whether `ode.solve()` exposes backend failure details and propagate them instead of treating every call as success.

- Keep plot/save tests out of ordinary unit tests.
  Existing tests call plotting and saving paths directly. LSODE2 tests should avoid filesystem/graphics side effects unless they use temp directories and are explicitly marked as integration tests.

## P1: API and UX

- Introduce typed state transitions for `KinData`.
  A safer user story would be: load/select reactions -> parse -> analyze stoichiometry -> compute constants -> export/format. Today many methods assume hidden prior calls and panic when state is missing.

- Split library loading from reaction lookup.
  `get_reactions_from_shortcuts` calls `open_json_files` per reaction. Load each library once, cache it, then resolve all reaction ids.

- Avoid public mutable database internals.
  `KineticData` exposes mutable fields like `Reactbase` and `dict_reactions`. Use snake_case and controlled methods, or make raw structures clearly low-level.

- Make output methods return structured data or strings.
  `pretty_print_*` can remain as convenience wrappers, but core code should return tables/records that GUI, CLI, notebooks, and tests can reuse.

- Clarify naming.
  Mixed names like `K_sym_vec`, `vec_of_reaction_Values`, `stecheodata`, `every_reaction`, `shortcut_reactions`, and `map_of_reactions` obscure which level of representation each field stores. Align naming with a small data model: shortcut, equation, library metadata, parsed reaction, stoichiometry, kinetic expression.

## P1: tests missing or too weak

- Add negative tests for user-facing failures:
  invalid shortcut format, unknown library prefix, missing reaction id, malformed equation, missing arrow, missing pressure/concentration data, empty temperature range, `n < 2` in temperature sweep, and calling result methods before solve.

- Add regression tests for the P0 correctness bugs:
  falloff `troe` preservation, three-body symbolic sign, pressure interpolation order/out-of-range behavior, direct-reaction workflow, and `create_kinetics_document` with several reactions from one library.

- Add tests for `KinData` state consistency.
  After append/remove/manual/direct/shortcut workflows, all per-reaction collections should have matching lengths or the normalized reaction record should make this invariant trivial.

- Add formula and stoichiometry edge-case tests.
  Cover parentheses, multipliers after brackets, fractional stoichiometry, reversible arrows, species names with charges or phase markers if supported, and malformed formulas.

- Fix test discovery gaps.
  `mechfinder_api.rs` contains `test_pres_deserialization` without `#[test]`, so it is not executed as an ordinary test.

## P2: performance and maintainability

- Remove repeated full-database clones and repeated JSON reads.
  `kinetics_lib_api.rs` and `mechfinder_api/mechfinder.rs` clone or reload large maps in hot lookup flows.

- Replace stdout debugging with logging or caller-controlled reporting.
  The current output volume makes CLI use, tests, and GUI integration noisy.

- Split large modules by responsibility.
  Good candidates: errors, reaction identifiers/shortcuts, library IO, reaction metadata, stoichiometry, symbolic constants, document/export, and IVP models.

- Add examples as API contracts.
  Short examples for `KinData` shortcut workflow, direct equation workflow, symbolic constants, kinetic document generation, and solid-state IVP with LSODE2 would double as usability checks.

## Suggested verification commands after implementation

```powershell
cargo test Kinetics::User_reactions_tests --lib
cargo test Kinetics::mechfinder_api --lib
cargo test Kinetics::solid_state_kinetics_IVP --lib
cargo test lsode2 --lib
git diff --check
```

## Open questions before coding

- Should `SolverType::LSODE2` become the public default for `KineticModelIVP`, or only the recommended constructor?
- Is the kinetic power in `stoichiometry_analyzer.rs` supposed to default to the stoichiometric coefficient or to `1.0` unless explicitly specified?
- How much backward compatibility is required for current public mutable fields and printing behavior?
- Should `KinData` preserve legacy parallel vectors temporarily, or should the next refactor move directly to a normalized per-reaction record?
