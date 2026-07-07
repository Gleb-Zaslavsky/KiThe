# Kinetics TODO analysis

Scope: `src/Kinetics/`, excluding `src/Kinetics/experimental_kinetics/`.

This file is an analytical TODO, not an implementation plan locked in stone. It records current risks in correctness, API ergonomics, test coverage, and solver integration.

## P0: correctness and panic-safety

### Done so far

- [x] Replace user-facing panic and unwrap paths in the main runtime flows with typed errors.
- [x] Fix direct-reaction workflow so direct equations can go straight through analysis.
- [x] Fix document generation with library maps.
- [x] Fix falloff metadata loss.
- [x] Fix symbolic and numeric mismatch for three-body reactions.
- [x] Fix pressure-dependent interpolation.
- [x] Revisit stoichiometric power parsing.
- [x] Split library loading from reaction lookup.
- [x] Tighten mutation invalidation for `KinData`.
- [x] Make `KinData` parallel/cache fields crate-private and expose read-only facade accessors.

- Replace user-facing `panic!`, `assert!`, `unwrap`, and `expect` paths with typed errors.
  Main hot spots: `User_reactions.rs`, `mechfinder_api.rs`, `kinetics_lib_api.rs`, `stoichiometry_analyzer.rs`, `molmass.rs`, `solid_state_kinetics_IVP.rs`.
  Public API should prefer `Result<_, KineticsError>` for invalid shortcuts, missing libraries/reactions, malformed equations, missing pressure/concentration data, failed file IO, JSON parse errors, and solver failures.
  Progress: `mechfinder_api.rs`, `mechfinder_api/kinetics.rs`, `mechfinder_api/mechfinder.rs`, `User_reactions.rs`, `stoichiometry_analyzer.rs`, `molmass.rs`, and the main `User_substances2.rs` fallback path now use typed errors instead of user-facing panics in their primary runtime flows. Remaining `unwrap`/`expect` sites are concentrated in tests, examples, and a few older compatibility paths.

- Normalize `KinData` state.
  `KinData` is the top-level facade, but it currently stores reaction state in many parallel optional vectors: shortcuts, maps, equation strings, parsed values, reaction data, symbolic constants, stoichiometry, and `every_reaction`. This makes it easy to create inconsistent states. `every_reaction` looks like the right future single source of truth, but it is not yet the real backbone.
  Status: partially addressed in `User_reactions.rs` with `KinDataState`, `ReactionIdentity`, `EveryReaction`, `ReactionBranchKind`, builder/helpers, append/remove invalidation, direct-reaction normalization, and a centralized checked state-transition helper. Legacy parallel fields still exist, but more read/write paths now invalidate stale derived views instead of silently drifting.
  Closed for now: the public shortcut/library/reaction-map accessors are canonical-only, while `get_reactions_from_shortcuts()` still resolves the active raw shortcut branch directly so the workflow keeps working without reintroducing legacy equality between caches and truth.
  The normalized export path now also preserves an already-coherent `every_reaction` cache when there is no parsed reaction data to rebuild from, so document export can keep using the prepared normalized view instead of discarding it.
  Progress: `reaction_branch_kind()` now prefers shortcut and library branches ahead of the direct raw-payload fallback, so `refresh_state_from_available_data()` no longer collapses shortcut/library workflows into direct just because raw values are present. New regression tests cover both shortcut and library raw-payload cases.
  Progress: branch snapshots no longer carry their own cached state tag; `refresh_state_from_available_data()` now derives the visible state from payload shape, with mechanism branches kept distinct from normalized sorted views. That removes one more duplicated source of truth and one more dead helper.
  Progress: the mutation tail now carries an explicit optional state override, so parsing and analysis steps no longer hand-roll separate `mark_*` calls after the shared refresh helper. That keeps the state publication path centralized while still letting workflow stages publish `ReactionDataParsed` and `Analyzed` intentionally.
  Progress: the dead `map_of_reactions` storage layer has been removed from `KinData`; grouped reaction views now live only as derived compatibility output from canonical reaction ids.
  Progress: `rebuild_reaction_ids_from_available_metadata()` now only reconstructs ids when they are missing, so stale legacy metadata can no longer overwrite already canonical ids during refresh.
  Progress: the normalized snapshot path is now split into a pure capture helper plus a single id-refresh step, so `sync_reaction_views_after_mutation()` no longer rebuilds canonical ids twice in one mutation flow.
  Progress: the rate-constant sorting path now builds values and reordered views from the same normalized snapshot, so `calc_K_const_for_all_reactions*()` no longer rereads a separate live reaction branch after the initial normalization.
  Progress: the normalized snapshot used for rebuild and reorder no longer carries shortcut/pair compatibility fields as inputs; those labels are now derived from canonical ids, which further shrinks the number of quasi-sources inside the internal snapshot model.
  Progress: `rebuild_reaction_ids_from_available_metadata()` now expands to the widest available reaction row count before filling missing ids, so export and refresh paths can recover missing canonical ids from a fuller parsed or raw branch without weakening the canonical-first state selection.
  Progress: `get_reactions_from_shortcuts()` now prefers canonical shortcut ids over the legacy `shortcut_reactions` cache, so a stale shortcut list can no longer redirect library resolution when the canonical branch is already present.
  Progress: `get_reactions_from_shortcuts()` now rejects non-shortcut canonical ids instead of falling back to stale legacy shortcut names, so shortcut resolution cannot silently resurrect an incompatible branch.
  Progress: the public derived compatibility views now reject mixed canonical-id sets instead of partially projecting them, so shortcut, library-pair, and reaction-map outputs must describe one coherent branch or return `None`.
  Progress: raw append branch selection now classifies `ReactionIdentity::Document` through a single internal helper and a typed `DocumentSourceKind`, so `append_reaction_document_source()` no longer reasons about document-source strings directly.
  Progress: visible-state publication now goes through one internal helper for both ordinary refreshes and mutation tails, so the refresh path and the post-mutation path share the same state publication contract.
  Progress: branch classification helpers now hang off `ReactionIdentity` itself, so shortcut, library, and document detection share one internal vocabulary instead of repeating pattern matches in several `has_*_branch()` checks.
  Progress: workflow-state selection now snapshots branch kind, normalized view readiness, parsing, and analysis artifacts once, so `state_from_available_data()` and `validate_state_contract()` no longer each assemble the same branch picture from scratch.
  Progress: the `Analyzed` state now explicitly accepts existing analyzed artifacts as well as parsed data/equations, so the visible state and the validation contract no longer disagree about already-finished workspaces.
  Progress: metadata validation now checks the actual stored compatibility caches as well as the canonical rows, so stale `shortcut_reactions`, `vec_of_pairs`, and `every_reaction` lengths can no longer hide behind the derived accessors during normalization.
  Progress: canonical reaction counts now prefer canonical ids first, which keeps rebuild and refresh flows anchored to the canonical branch instead of letting raw payload counts outrank it.
  Progress: `KinData`'s parallel/cache fields are now `pub(crate)` instead of externally mutable `pub` fields. External users should go through facade accessors such as `canonical_reaction_ids()`, `shortcut_names()`, `library_pairs()`, `reaction_values()`, `reaction_data()`, `equations()`, `substances()`, `groups()`, `symbolic_constants()`, `normalized_reactions()`, and `stoichiometric_analyzer()`. This preserves the domain caches internally while removing the most dangerous external desynchronization path.
  Progress: the assignment/result distinction is now explicit in the normalized-read path. `reaction_ids` describes the requested/selected identities, while `every_reaction` is the materialized result used by export and formatting. A materialized `every_reaction` can be exported without a separate assignment cache, failed library resolution leaves the shortcut assignment intact without creating a result, and failed direct-equation analysis is transactional instead of leaving partial substance/stoichiometry state.
  Concrete next pass, in order:
  1. `refresh_state_from_available_data()` should become the only place that decides the visible state from the current payload shape.
  2. `rebuild_reaction_ids_from_available_metadata()` should be the only place that recreates canonical ids from legacy caches.
  3. `rebuild_every_reaction()` and `build_every_reaction_view()` should read from canonical ids only and stop consulting parallel branches directly.
  4. `sync_reaction_views_after_mutation()` and `build_reordered_reaction_views()` should keep derived views aligned from one shared normalized snapshot, not from several ad hoc sources.
  5. `set_reactions_directly()`, `set_reactions_from_shortcut_range()`, `append_reaction()`, `append_reaction_with_shortcut()`, `append_reaction_from_map()`, `remove_by_index()`, and `load_reactions_from_json()` should route through the same state-refresh helper instead of each one reasoning about consistency independently.
  6. `derived_reaction_map()`, `derived_library_pairs()`, and `derived_shortcut_names()` should keep shrinking their fallback behavior until canonical ids are present everywhere that matters.
  Finite implementation batch for the next pass:
  1. `refresh_state_from_available_data()` as the single visible-state publisher.
  2. `sync_reaction_views_after_mutation()` as the shared derived-view invalidation/rebuild helper.
  3. `install_reaction_branch()` as the transactional branch switch path.
  4. `write_reaction_payload()` and `append_reaction_with_shortcut()` as the two raw-mutation entry points.
  5. `remove_by_index()` as the structural mutation path that must keep every view aligned.
  Done means:
  - one canonical write path for ids and per-reaction views;
  - one canonical refresh path for visible state;
  - legacy fields surviving only as compatibility storage, not as equal readers of truth;
  - new tests that fail if stale shortcut/pair caches can override canonical ids again.
  Progress: the rebuild path now snapshots equations, shortcuts, library pairs, and symbolic constants once per normalized rebuild instead of rereading the same accessors inside the per-reaction loop. That makes the `every_reaction` rebuild more clearly one-pass and keeps the canonical snapshot obvious. The mutation tail for branch install, raw appends, shortcut appends, removal, parsing, equation-cache rebuilds, and symbolic-constant updates now also goes through a shared `finalize_reaction_mutation()` helper, so we are closing the last obvious sync+refresh duplication without widening the API surface. JSON shortcut loading now also goes through a dedicated `ReactionBranchSnapshot::shortcut_values()` constructor instead of hand-building the branch inline, so one more input path now shares the same branch-install contract. Raw appends now use a branch snapshot constructor too, and the append helper now preserves append-style accumulation for document branches while still restarting from zero when a non-document branch is replaced, so the branch semantics are explicit instead of accidental. Shortcut appends now remain shortcut-selected rather than masquerading as direct input, so the branch semantics are clearer and the stale-state replacement contract is explicit. The shortcut, library, and reaction-map derived views now share common helpers for reading canonical ids, which trims one more layer of duplicate fallback logic. Reordered views now also rebuild from the same normalized snapshot, so sorting no longer re-reads a separate pile of live fields.

- [x] Fix direct-reaction workflow.
  `set_reactions_directly` now feeds `kinetic_main` through the equation-only path when reaction values are absent, so direct reactions can go straight to analysis without pretending to be shortcut-driven input.

- [x] Fix document generation with library maps.
  `create_kinetics_document` enumerates `map_of_reactions` by library and indexes `vec_of_reaction_data` with that library index. If a library contains several reactions, the same `ReactionData` can be paired with multiple equations or indexing can drift. This is a real data corruption/panic risk.

- [x] Fix falloff metadata loss.
  In `mechfinder_api/kinetics.rs`, `FalloffStruct::new` accepts `troe`, but initializes `troe: None`. Troe/SRI metadata is silently discarded.

- [x] Fix symbolic/numeric mismatch for three-body reactions.
  `ThreeBodyStruct::K_const` uses `exp(-E / RT)`, while symbolic `K_expr` builds `Exp(E / RT)`. Symbolic output can therefore disagree with numeric kinetics.

- [x] Fix pressure-dependent interpolation.
  `PressureStruct::K_const` builds pressure arrays from `HashMap` key iteration, so ordering is not guaranteed before binary search. It also uses suspect low/high indexing in the interpolation branch and can call `location.unwrap()` after an `Err` result for out-of-range pressure.

- [x] Revisit stoichiometric power parsing.
  In `stoichiometry_analyzer.rs`, the branch that should probably copy a stoichiometric coefficient into the kinetic power currently uses comparison (`power_coeff == stec_coeff`) instead of assignment. Tests currently tolerate powers of `1.0` for explicit stoichiometric coefficients, so the intended chemistry contract needs to be clarified and then tested.

- Make compute methods pure by default.
  Methods such as `analyze_reactions`, `calc_K_const_for_all_reactions`, `calc_sym_constants`, molar-mass parsing, and stoichiometry parsing print during normal computation. Computation should return data; formatting/logging should be opt-in.
  Progress: the main `KinData` compute path no longer prints substance discovery or reaction-rate tables, so the numeric and parsing steps are much closer to pure functions. The pretty-print methods remain as explicit output helpers.

## P1: LSODE2 support in solid-state IVP

- Add an idiomatic LSODE2 entry point.
  `solid_state_kinetics_IVP.rs` accepts arbitrary `SolverType` in `KineticModelIVP::new`, but there is no convenient constructor or shortcut for the now-standard `SolverType::LSODE2`. Candidate API:
  - `KineticModelIVP::new_lsode2(...)`
  - `KineticModelIVP::with_lsode2_defaults(...)`
  - or a solver-profile builder that makes LSODE2 the ergonomic default while still allowing RK/BDF/Radau selection.
  Progress: `KineticModelIVP::new_lsode2()` now exists and a smoke test covers `F1` with `SolverType::LSODE2`. The ergonomic surface is still minimal and can be expanded later if we want a richer profile builder.

- Decide whether LSODE2 should become the default solver.
  If LSODE2 is the de facto default for RustedSciThe IVP work, `KineticModelIVP::new` should either default to it through a high-level constructor or the docs should explain when to choose LSODE2 versus the existing methods.
  Progress: `KineticModelIVP::default()` now routes to `LSODE2`, so the ergonomic default is explicit without changing the older `new(solvertype)` constructor.

- Add focused LSODE2 tests.
  Required coverage:
  - smoke solve for a simple solid-state model, e.g. `F1`, using `SolverType::LSODE2`;
  - comparison against an existing trusted solver on the same problem within domain-appropriate tolerance;
  - a stiffer/high-activation-energy case where LSODE2 is expected to be useful;
  - result shape invariants: `time`, `temperature`, `conversion`, `conversion_rate`, `rate_constant` lengths match;
  - no-solve and invalid-parameter cases return errors instead of panicking;
  - solver status/statistics are accessible or at least solver failure is not hidden.
  Note: `IVPSolution` currently exposes `time`, `temperature`, `conversion`, and `conversion_rate`; `rate_constant` is not yet part of the public solution record, so the shape contract still needs a decision.
  Progress: smoke coverage now exists for `LSODE2`, a direct `LSODE2` vs `BDF` regression test now compares interpolated conversion values on a shared time window for `F1`, and a stiffer `D3` smoke test now covers the same solver path.

- Check solver error propagation.
  `KineticModelIVP::solve` currently calls the universal solver and then returns `Ok(())`. Confirm whether `ode.solve()` exposes backend failure details and propagate them instead of treating every call as success.
  Progress: `solve()` now uses `UniversalODESolver::try_solve()` so backend solver failures can bubble up as `Err(String)` instead of being hidden behind a blanket success.

- Keep plot/save tests out of ordinary unit tests.
  Existing tests call plotting and saving paths directly. LSODE2 tests should avoid filesystem/graphics side effects unless they use temp directories and are explicitly marked as integration tests.

## P1: API and UX

- Introduce typed state transitions for `KinData`.
  A safer user story would be: load/select reactions -> parse -> analyze stoichiometry -> compute constants -> export/format. Today many methods assume hidden prior calls and panic when state is missing.
  Progress: `KinDataState` now acts as a real gate for the analysis pipeline. `reactdata_parsing()`, `analyze_reactions()`, and `calc_sym_constants()` check the current workflow state before proceeding, so invalid stage jumps now return `InvalidState` instead of quietly running in the wrong phase. State publication itself now also goes through a checked helper, so internal transitions cannot silently publish an invalid phase.

- [x] Split library loading from reaction lookup.
  `get_reactions_from_shortcuts` now caches `KineticData` per library so each JSON library is parsed once per workflow, then all reaction ids are resolved from that cache.

- [x] Tighten mutation invalidation for `KinData`.
  `append_reaction*`, `set_reactions_directly`, and the builder now clear stale derived views so a new input branch does not retain old parsed/analyzed state.

- [x] Make `library_pairs()` a derived compatibility view.
  The facade now rebuilds library/reaction pairs from canonical ids instead of treating `vec_of_pairs` as the main read source, so the pair cache is no longer the primary contract.

- [x] Make `KinData` caches read-only for external callers.
  The reaction ids, shortcut cache, library-pair cache, raw payloads, parsed data, equation cache, stoichiometry, symbolic constants, workflow state, and normalized reaction view are now crate-private fields. Public callers can inspect them through accessors, but cannot mutate them directly and bypass the normalization helpers.

- Avoid public mutable database internals.
  `KineticData` exposes mutable fields like `Reactbase` and `dict_reactions`. Use snake_case and controlled methods, or make raw structures clearly low-level.

- Make output methods return structured data or strings.
  `pretty_print_*` can remain as convenience wrappers, but core code should return tables/records that GUI, CLI, notebooks, and tests can reuse.
  Progress: the main pretty-print flows now have structured builders behind them. `substances_verbose_tables()` returns the stoichiometric and elemental tables, `normalized_reaction_data_table()` exposes the already-normalized reaction table without mutating state, `reaction_data_table()` normalizes on demand, and `kinetics_document_map()` / `kinetics_document_map_from_normalized()` build the export payload without writing to disk. The normalized export path now relies on `every_reaction` first and only checks raw reaction values when they still exist as a consistency hint. Empty normalized views now fail with a consistent typed error and have regression coverage. The human-facing pretty-print methods now operate on a cloned normalized workspace, so output no longer mutates the caller's live `KinData` state.
  Progress: pure normalized print helpers now have regression coverage that confirms they read an already prepared normalized cache without mutating workflow state.
  Progress: repeated print/export reads now short-circuit when the normalized view is already coherent, so already-normalized workspaces do not pay for a full rebuild or a second clone just to render the same data again.

- Clarify naming.
  Mixed names like `K_sym_vec`, `vec_of_reaction_Values`, `stecheodata`, `every_reaction`, `shortcut_reactions`, and `map_of_reactions` obscure which level of representation each field stores. Align naming with a small data model: shortcut, equation, library metadata, parsed reaction, stoichiometry, kinetic expression.
  The builder surface is also being trimmed so the public API has fewer overlapping entry points for the same workflow shape. Redundant builder branches for raw values and parsed reaction data were removed in favor of the clearer direct/shortcut/library paths. There is now a cleaner facade layer with `workflow_state()`, `canonical_reaction_ids()`, `shortcut_names()`, `library_pairs()`, `reaction_values()`, `reaction_data()`, `equations()`, `symbolic_constants()`, and `normalized_reactions()`, plus `with_*` builder aliases and `from_*` constructors that are more idiomatic than the older method names. Builder setters for shortcut and library workflows now also return typed results instead of relying on internal `expect`, so the remaining cleanup is mostly about tightening the legacy field names and shrinking the internal model further.

## P1: tests missing or too weak

- Add negative tests for user-facing failures:
  invalid shortcut format, unknown library prefix, missing reaction id, malformed equation, missing arrow, missing pressure/concentration data, empty temperature range, `n < 2` in temperature sweep, and calling result methods before solve.
  Progress: `get_result()` and `get_solution()` now return a structured error when `solve()` has not been called, and `plot_in_terminal()` no longer panics on an empty solver state. Shortcut-range validation now also rejects missing separators like `C1-C5`. Library lookup and missing-reaction paths now have regression tests in both `kinetics_lib_api.rs` and `User_reactions_tests.rs`, and the reaction parser now has negative coverage for malformed payloads and unknown reaction types. `parse_kinetic_data_vec` now also rejects missing reaction equations and missing kinetic parameters through explicit regression tests.

- Add regression tests for the P0 correctness bugs:
  falloff `troe` preservation, three-body symbolic sign, pressure interpolation order/out-of-range behavior, direct-reaction workflow, and `create_kinetics_document` with several reactions from one library.

- Add tests for `KinData` state consistency.
  After append/remove/manual/direct/shortcut workflows, all per-reaction collections should have matching lengths or the normalized reaction record should make this invariant trivial.
  Progress: append/remove/direct/shortcut paths now have regression tests for invalidation and failure cases, plus a consolidated invariant test covering raw and normalized workflows. Repeated raw appends are also covered so canonical ids cannot drift during incremental input, library-only removal now preserves `LibraryResolved`, stale-state reset before `set_reactions_directly` is now explicitly covered too, and `load_reactions_from_json` now has a stale-state replacement regression test. Export normalization now also has coverage for missing equations in a normalized view and for rebuilding raw payload ids from the current data branch.

- Add stronger negative tests around reaction parsing and numeric sweeps.
  The `mechfinder_api` layer now has direct coverage for mismatched reaction-type tags and malformed payloads, and the numeric sweep path rejects `n < 2` at the reaction-data layer. This still leaves room for formula-parser edge cases and a few specialized malformed-equation cases.

- Add formula and stoichiometry edge-case tests.
  Cover parentheses, multipliers after brackets, fractional stoichiometry, reversible arrows, species names with charges or phase markers if supported, and malformed formulas.
  Progress: `analyse_substances` now defaults kinetic powers to stoichiometric coefficients unless an explicit `^` override is present, and formula parsing now has regression coverage for phase markers combined with groups plus a molar-mass smoke test on phased input.

- Fix test discovery gaps.
  `mechfinder_api.rs` contains `test_pres_deserialization` without `#[test]`, so it is not executed as an ordinary test.

## P2: performance and maintainability

- Remove repeated full-database clones and repeated JSON reads.
  `kinetics_lib_api.rs` and `mechfinder_api/mechfinder.rs` clone or reload large maps in hot lookup flows.
  Progress: the shared library loaders now cache the parsed `Reactbase` and reaction-relationship JSON in memory, so repeated library loads and mechanism searches no longer reread and reparse the same files on every call. The substance-matching search path now also checks subsets without cloning each relationship set first, so one more hot allocation path is gone. The next cleanup pass can focus on shrinking the remaining per-query clones instead of paying the disk/parsing tax repeatedly.
  Progress: `mechfinder_api/mechfinder.rs` now walks the cached owned reaction database directly instead of building a borrowed mirror of the whole map first, so the mechanism search no longer pays for a second full database conversion on each call.
  Progress: `mechfinder_api/mechfinder.rs` now returns a typed `Result`, reuses the cached libraries without `unwrap()`, and `parse_kinetic_data_vec` consumes reaction payloads by value instead of cloning each JSON `Value` before parsing.
  Progress: reaction cleanup inside `mechfinder.rs` now keeps the owned `String` instead of cloning after `clean_off_DUP`, and the parsing error path keeps the old user-facing wording while avoiding a second payload clone.

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
