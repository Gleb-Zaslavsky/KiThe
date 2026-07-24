# Thermodynamics DB and Aggregator Audit

## Scope

This checklist covers the thermodynamic/transport database layer and the
`SubsData` property aggregator:

- `DBhandlers.rs` and `DBhandlers/**`
- `thermo_lib_api.rs` and `prelude.rs`
- `User_substances.rs`, `User_substances2.rs`, and
  `User_substances_presets.rs`
- `User_substances_error.rs` and the corresponding tests
- `dG_dS.rs` only where it consumes or mutates `SubsData` caches
- `library_manager2.rs` only where it implements the internal thermodynamics
  write/drop workflow tracked by P2.5

The following areas were deliberately out of scope for the original database
and aggregator audit:

- `ChemEquilibrium/**` and both chemical-equilibrium calculation branches
- phase/solution thermodynamic modelling, which is now covered by the dedicated
  audit and execution plan at the end of this document

## Current Data Flow

1. `ThermoRepository` loads and validates local JSON libraries once per
   explicit repository lifecycle; `ThermoData` is a lightweight query handle
   over the shared immutable catalogs.
2. `SubsData` accepts substance names, library preferences, explicit lookup
   instructions, temperature, pressure, and units.
3. `search_substances()` resolves raw records deterministically into canonical
   per-property `SubstanceSearchState` values with provenance and calculators.
4. Calculation methods publish numeric values, closures, symbolic expressions,
   composition, molar mass, diffusion data, and derived-property caches behind
   controlled mutation and read-only readiness accessors.

The original public-mutation, non-deterministic lookup, partial-publication,
and repeated-library-loading problems are resolved. Remaining closeout work is
concentrated in legacy Gibbs/entropy consumers, NIST/CEA fitting failure paths,
property-report consistency, and deterministic offline test coverage.

## Priority Definitions

- **P0**: may panic on user/library data, return physically invalid values,
  silently lose failures, or select the wrong database record.
- **P1**: architectural changes required to make state, errors, I/O, and lookup
  behaviour explicit and reliably maintainable.
- **P2**: performance, API clarity, deterministic output, and separation of
  computation from presentation.
- **P3**: long-term test infrastructure, documentation, benchmarks, and module
  cleanup.

# P0: Correctness and Failure Containment

P0 must be completed before changing the public API or optimizing the hot read
paths. Every item requires regression tests that fail against the old
behaviour.

## P0.1 Remove panics from library and record parsing

- [x] Introduce a typed `ThermoLibraryError` (or equivalent) covering file
  open/read failures, JSON decoding, malformed records, unsupported libraries,
  and schema validation.
  - Replace the infallible `ThermoData::new()` I/O path with `try_new()` or an
    injected repository constructor.
  - Do not silently replace a failed elements-library load with an empty map.
  - Keep a convenience constructor only if its failure policy is explicit.
- [x] Make NASA composition deserialization fully fallible.
  - Remove `unwrap()` calls from `NASAdata::deserialize_composition`.
  - Reject missing element names, missing coefficients, malformed numbers, and
    duplicate malformed entries through `serde::de::Error::custom`.
- [x] Validate the complete NIST interval schema before indexing it.
  - Cover `NISTdata::extract_coefficients`, `NistInput::extract_coefficients`,
    `NistInput::pretty_print`, and `NISTdata::get_coefficients`.
  - Validate required `T` and `Cp` fields, equal interval counts, two bounds per
    interval, minimum coefficient count, ordered finite bounds, and finite
    coefficients.
- [x] Replace user-facing panics in calculator factories.
  - `thermo_api::create_thermal_by_name`
  - `transport_api::create_transport_by_name`
  - `SubsData::what_handler_to_use`
  - Return an `UnsupportedLibrary` error containing the requested name.
- [x] Remove panic paths from `ThermoData::remove_by_name()` and
  `ThermoData::create_substance_document()`.
  - Missing substances and line-read errors must be typed results.

Required tests:

- [x] Corrupt and missing library files return typed errors without panicking.
- [x] Malformed NASA composition strings return serde errors.
- [x] Missing, uneven, short, unordered, and non-finite NIST intervals are
  rejected before calculator construction.
- [x] Unknown calculator/library names return `UnsupportedLibrary`.
- [x] Missing remove targets and malformed document input remain recoverable.

Completion criterion: no production parsing, library-loading, or calculator
factory path in scope contains `panic!`, `assert!`, `unwrap()`, or `expect()`.

## P0.2 Enforce finite physical inputs and outputs

- [x] Add shared validation helpers for finite values, finite positive values,
  ranges, and compatible vector lengths.
  - Transport density/temperature/pressure/molar-mass helpers now share a
    local finite-value validation layer, and transport prerequisites validate
    collision diameter plus the other scalar inputs before derived formulas run.
- [x] Transport scalar inputs reject `NaN`, infinities, and non-positive
  values before derived calculations.
- [x] Pair diffusion calculations reject non-finite temperature, pressure, and
  mass state before numerical evaluation.
- [x] Apply finite-positive validation to temperature, pressure, density, molar
  mass, collision diameter, heat capacity, viscosity inputs, and every divisor
  used by transport/diffusion formulas.
- [x] Reject `NaN`, positive/negative infinity, zero where physically invalid,
  and negative values before derived calculations.
- [x] Validate calculated scalar outputs and report `NonFiniteOutput` with the
  property and substance name instead of caching `NaN`/`Inf`.
- [x] Make `PairDiffusion` and `MultiSubstanceDiffusion` construction and
  calculation fallible.
  - Validate `T`, `P`, `M`, `sigma`, and `e_k`.
  - Validate consistency between `substance_names` and the substance map.
  - Constructor-level rejection now happens for invalid thermodynamic state;
    derived coefficients are still checked again before use.
- [x] Align TRANSPORT trait and inherent APIs.
  - The density fallback is now centralized and consistent across numeric,
    closure, and symbolic transport paths.
  - `ideal_gas()` and `ideal_gas_sym()` now validate their thermodynamic state
    before publishing density results.

Required tests:

- [x] Table-driven `NaN`, `+Inf`, `-Inf`, zero, and negative input tests for
  thermo, transport, and diffusion APIs.
- [x] Table-driven finite-input validation covers transport scalar helpers.
- [x] Pair diffusion rejects non-finite state before computation.
- [x] No public calculation path returns or caches a non-finite value as a
  successful result.
- [x] Numeric, closure, and symbolic transport paths apply the same prerequisite
  and density rules.
- [x] Diffusion rejects inconsistent maps and invalid pair parameters without
  indexing panics.

Completion criterion: successful numeric results are finite by contract, and
invalid physical inputs cannot enter derived caches.

## P0.3 Make database lookup deterministic and property-specific

- [x] Replace `HashMap` iteration as the effective library-priority order in
  `SubsData::search_substances()` with an ordered lookup plan.
  - Priority, permitted, and element-based lookup now traverse libraries in a
    stable sorted order instead of relying on hash iteration order.
- [x] Track thermo and transport lookup independently.
  - Finding thermodynamic data must not suppress permitted transport fallback.
  - Finding transport data must not suppress permitted thermodynamic fallback.
- [x] Stop overwriting successful same-priority results according to randomized
  map iteration order.
- [x] Preserve explicit lookup provenance as `LibraryPriority::Explicit` rather
  than converting it to ordinary priority lookup.
- [x] Ensure a failed explicit lookup updates only the requested property and
  does not erase an independently found property.
- [x] Apply the same deterministic rule to element-library lookup in
  `populate_element_search_results()`.
- [x] Replace production `assert!`/`unwrap()` in permitted-library traversal
  with typed internal-invariant errors.

Required tests:

- [x] Repeated searches produce identical sources and results across runs.
- [x] Thermo-found/transport-fallback and transport-found/thermo-fallback
  scenarios both work.
- [x] Two matching libraries at the same priority resolve by documented order.
- [x] Explicit lookup success/failure preserves the other property category.
- [x] Element lookup has deterministic provenance.

Completion criterion: lookup results depend only on the declared ordered
policy, never on hash iteration or the success of an unrelated property search.

## P0.4 Stop reporting partial failure as unconditional success

- [x] Replace batch methods that log per-substance errors and then return
  `Ok(())` with a fail-fast policy that only publishes a rebuilt cache after
  the full loop succeeds.
- [x] Apply the policy consistently to numeric, closure, and symbolic thermo
  maps in `User_substances.rs`.
- [x] Apply the same policy to transport maps in `User_substances2.rs`.
- [x] Apply it to NIST fallback/network lookup: an all-failed request must not
  be reported as success.
- [x] Keep the error logger as presentation/diagnostic support, not as a
  replacement for the return contract.

Required tests:

- [x] One success plus one failure returns the first error and leaves the old
  cache untouched.
- [x] All-failed batches cannot return an empty successful result.
- [x] Fail-fast mode publishes no partially rebuilt cache.
- [ ] Best-effort mode identifies every failed substance and property.

Completion criterion: callers can distinguish success from failure without
inspecting logs, and a failed batch never replaces a previously valid cache
with a partially rebuilt one.

## P0.5 Invalidate stale aggregator state on input mutation

- [x] Define the dependency graph for substances, lookup instructions,
  temperature, pressure, units, phases, and derived caches.
  - The invalidation helpers now document the cache dependency graph directly
    in `User_substances.rs`, so the policy lives next to the implementation.
- [x] Make `set_substances()` invalidate search results and every dependent
  thermo, transport, composition, molar-mass, diffusion, and derived cache.
- [x] Add controlled setters for all other inputs and apply the same targeted
  invalidation policy.
  - `set_library_priority`, `set_multiple_library_priorities`,
    `set_explicis_searh_instructions`, `set_T`, `set_P`, and `set_M` now
    invalidate the caches that depend on them.
- [x] Prevent public field mutation from bypassing invalidation (completed
  structurally in P1.2).
  - `SubsData` mutation-critical fields are now crate-private and can only be
    updated through controlled setters inside the module.
- [x] Build replacement maps locally and publish them only after the selected
  batch policy completes.

Required tests:

- [x] Mutation-path invariant test covering substances, T, P, units, library
  policy, and explicit instructions.
- [x] No cache can retain entries for removed substances.
- [x] Failed recomputation does not mix old and new values.

Completion criterion: after every supported mutation, all readable caches are
either valid for the current inputs or explicitly absent.

## P0.6 Remove residual production panic paths

- [x] Refactor public calculation paths in `dG_dS.rs` to consume typed
  readiness accessors and return `SubsDataResult` instead of unwrapping phase,
  numeric, symbolic, and closure caches.
  - Gibbs/entropy/Helmholtz read paths now degrade to empty maps for incomplete
    calculator state and are covered by a regression test; the full typed error
    propagation still remains to be finished. The error boundary around
    aggregate calculation failures now preserves nested sources for the main
    thermo/NIST paths, so the remaining work is narrower than before.
  - Gas-correction calculations no longer index concentration vectors directly;
    short concentration vectors now degrade safely with a warning instead of a
    panic.
  - The `PhaseOrSolution` and `OnePhase` facades now propagate typed
    calculation failures instead of silently discarding them.
- [x] Replace `thermo.get_composition().unwrap()` in
  `SubsData::calculate_elem_composition_and_molar_mass_local()` with propagated
  calculator context and preserve the fallback formula parser only for a valid
  `Ok(None)` result.
  - The composition/molar-mass path now uses a dedicated helper that prefers
    calculator-provided composition, falls back to formula parsing only when
    needed, and leaves no partial state behind on failure.
- [x] Make NIST LM fitting treat missing/non-finite solver results and poor
  termination as `NISTError::FittingError`; remove production `assert!` and
  `unwrap()` from `fitting_adjacent_weighted()` and `fitting_non_adjacent()`.
  The LM branches now validate published `R^2`, coefficient count, and
  finiteness before any coefficients are published. A dedicated regression test
  covers the new typed error gate, and the empty-input/empty-interval
  guardrails remain in place for the local fitting helpers.
- [x] Make CEA adjacent fitting validate symbolic expressions, the complete LM
  coefficient map, and `R²` before publishing `coeff_Visc`/`coeff_Lambda`.
  The adjacent fit now works on temporary clones and only publishes the
  fitted coefficients after all validation gates pass.
- [ ] Add regression tests proving each path returns a typed error and leaves
  the previously published state unchanged when its prerequisite or solver
  result is absent. The NIST empty-vector and empty-selection regressions are
  now covered, and the adjacent/non-adjacent NIST fit failure paths now verify
  that previously published coefficients survive an error return. The remaining
  work is the broader typed-error sweep.

Completion criterion: the scoped production paths contain no panic triggered
by user, library, calculator, or optimizer state.

# P1: Canonical State, Repository, and Error Architecture

## P1.1 Introduce a shared immutable library repository

- [x] Separate library loading from `ThermoData` query operations.
- [x] Represent loaded catalogs with an immutable repository that can be shared
  through `Arc` or an explicit application-owned handle.
- [x] Avoid opening and deserializing the same JSON files for every
  `SubsData::new()`.
- [x] Preserve an explicit reload path for development/custom-library workflows;
  do not hide mutable reload behaviour behind a process-global cache.
- [x] Inject repository paths/data in tests instead of relying on the working
  directory and installed local files.
- [x] Store library capabilities and aliases in one registry instead of
  inferring the handler from the first JSON record (for example, by looking for
  a `Cp` field).

Required tests:

- [x] Two aggregators share one immutable repository without sharing mutable
  calculation state.
- [x] Custom temporary repositories load independently.
- [x] Reload is explicit and old handles remain internally consistent.
- [x] Empty libraries and mixed/invalid schemas return typed errors.

## P1.2 Give `SubsData` one controlled state boundary

- [x] Make assignment and derived fields private where external mutation can
  break invariants.
- [x] Split the model conceptually into:
  - user query/configuration,
  - immutable catalog handle,
  - per-property lookup results with provenance,
  - numeric/symbolic/closure derived caches,
  - diagnostics.
  - The canonical typed boundary is now in place: `search_states` carries the lookup contract, the derived caches have explicit readiness states, and the public API is read-only by default.
  - Several long-form tests now query `SubstanceSearchState` and derived readiness helpers directly instead of reading compatibility snapshots first.
  - Remaining work is narrower now: keep shrinking compatibility reads and remove the last direct legacy-map dependencies from internal helpers and tests.
- [x] Replace `HashMap<WhatIsFound, Option<SearchResult>>` with a typed
  per-property state such as `NotSearched`, `Found`, `Missing`, or `Failed`.
  - `search_states` has been introduced as the canonical layer; the legacy
    compatibility view still remains for now.
  - The main lookup queries now read from `search_states`, so `search_results`
    is increasingly treated as a compatibility/export view instead of the
    authoritative state.
  - The element-search path and typed lookup summary helpers now also read
    from `search_states`, so fewer code paths depend on the legacy map.
  - `get_substance_result()` and `get_all_results()` now reconstruct their
    compatibility snapshots from `search_states` instead of returning the
    stored legacy map directly.
  - Calculator access helpers now resolve mutable `SearchResult` instances from
    `search_states`, so lookup operations no longer depend on the legacy map.
- [x] Decide how canonical calculator data is represented inside `SubsData`.
  - `search_states` is the authoritative lookup layer.
  - `search_results` remains a compatibility/export snapshot only.
  - `SearchResult` now behaves as a read-only payload: provenance JSON and the
    attached calculator live together inside the canonical found record instead
    of as parallel public mutable fields.
- [x] Expose read-only getters and narrow mutation methods that centralize
  invalidation.
  - Public read access now goes through typed getters on `SearchResult`,
    `SubstanceSearchState`, and `SubsData` instead of direct field mutation.
  - The read-only surface now covers config, repository, lookup state, and the
    major derived caches.
- [x] Represent "not calculated" explicitly instead of initializing values and
  closures to valid-looking zero/identity placeholders.
  - The lookup layer now exposes explicit pending/resolved/missing/failed
    helpers, and the numeric caches now have explicit `NotCalculated`/`Ready`
    read-only accessors instead of forcing callers to infer readiness from map
    shape.
  - Transport calculation no longer routes dummy CEA placeholder values through
    the solver; the branch selection is now explicit.
  - The first production consumer has already switched to the explicit
    readiness accessor, so the migration is now real rather than purely
    theoretical.
  - Gibbs/entropy consumers now also read `dH` and `dS` through the explicit
    readiness state instead of reaching into the nested cache directly.
  - Closure caches now also have typed readiness accessors, so tests and
    consumers can inspect function availability without pattern matching on the
    raw nested map shape.
  - Thermo and transport regression tests now assert readiness explicitly
    instead of only peeking into nested `Option` maps.
  - Temperature invalidation now has a dedicated regression test that proves
    numeric and transport caches return to `NotCalculated` while reusable thermo
    closures remain explicitly `Ready`.

Required tests:

- [ ] Compile-time/API tests demonstrate that callers cannot mutate caches
  directly.
- [ ] State-transition tests cover new, searched, partially found, calculated,
  invalidated, and failed states.
- [ ] Search provenance survives calculation and export.

## P1.3 Make clone semantics truthful and cheap

- [x] Audit `SubsData::clone()` and remove silent cache reconstruction.
  - The current implementation cannot clone boxed closures faithfully.
  - It silently ignores closure recreation failures.
  - Transport closure state can be lost or rebuilt into the wrong map.
- [x] Prefer `Arc<dyn Fn... + Send + Sync>` for shareable compiled closures.
  - CEA transport, transport property, and pair-diffusion closures now clone by
    sharing the compiled callable instead of recreating a placeholder.
- [x] Apply the same rule to NASA/NIST/CEA/TRANSPORT calculators and diffusion
  objects; a clone must not replace a working closure with a zero placeholder.
  - NASA/NIST clones still rebuild from coefficients, while the transport and
    diffusion calculators now preserve active callable state directly.
- [x] If semantic cloning remains impossible, remove `Clone` and provide an
  explicit fallible rebuild operation.
  - The remaining scoped calculators can now clone truthfully, so the fallback
    path is not needed here.

Required tests:

- [x] Original and clone produce equivalent numeric and closure results.
- [x] Clone does not perform file I/O or silently recalculate state.
- [x] A clone failure is explicit rather than swallowed.

## P1.4 Consolidate typed errors without losing source context

- [ ] Define error boundaries for repository I/O/schema, lookup, calculator
  construction, physical validation, network access, and batch aggregation.
- [x] Preserve `source()` chains instead of converting network/serde failures
  into generic `NoCoefficientsFound` values.
  - Thermo, NIST, and transport wrappers now expose the underlying source
    chain instead of flattening it into an ambiguous lookup-style outcome.
- [x] Preserve those concrete sources through the aggregator boundary.
  `SubsDataError::NistRetrievalFailed` and `SubsDataError::CalculationFailed`
  now retain an optional nested source so callers can inspect the original
  calculator or NIST error when one is available.
- [x] Remove duplicated error-message mapping between format errors,
  `ThermoError`/`TransportError`, `SubsDataError`, and `LogErrorType` where a
  structured source can be retained.
  - The owned and borrowed conversions for thermo, transport, and `SubsDataError`
    now share local helpers instead of carrying duplicated `match` blocks.
- [x] Fix the most misleading user-facing message paths, especially NIST
  retrieval failures that were reported as missing coefficients.
- [x] Replace `Box<dyn Error>` in public file APIs with stable typed errors.
  - No in-scope `DBhandlers`/`User_substances` public file API currently uses
    `Box<dyn Error>`; the remaining occurrences are in `ChemEquilibrium`, which
    was intentionally excluded from this thermodynamics audit.

Required tests:

- [x] Error variants retain path, library, substance, property, and underlying
  source where applicable.
- [x] Display strings are format-specific and actionable.
- [x] Error conversion does not turn infrastructure failure into data absence.
- [x] Aggregator-level NIST fallback and calculation failures expose their
  original typed error through `source()` instead of only formatted text.
  Constructor-level coverage is in place; add one end-to-end regression for the
  full fallback path before closing this fully.

## P1.5 Harden the NIST network boundary

  - [x] Inject/configure the HTTP client with finite connect/request timeouts and
    a stable user agent.
  - [x] Define retry behaviour explicitly; do not retry schema/4xx failures.
  - [x] Bound accepted response sizes and validate content before parsing.
- [x] Separate offline deterministic tests from live NIST integration tests.
    Parser integration tests are isolated, and the live NIST suites are now
    marked `#[ignore]` instead of running by default.
  - [x] Return a batch report for multi-substance fetching and fail truthfully if
    every request fails.
  - [x] Remove direct URL/data printing from the network and parser core.

Required tests:

- [x] Mocked timeout, HTTP error, malformed HTML, missing coefficients, partial
  success, and complete failure.
- [x] Ordinary CI does not require network access.
- [x] Optional live tests are explicitly ignored.

# P2: Performance, API, and Maintainability

## P2.1 Remove repeated large clones and redundant scans

- [ ] Profile before changing algorithms and record representative catalog sizes.
- [x] Pass borrowed records through lookup instead of cloning full JSON values,
  search-result maps, and allowed-library vectors.
- [x] Replace explicit-instruction membership scans with a set-based lookup in
  `search_substances()`.
- [x] Make `ThermoData::search_libs()` return matching entries rather than
  allocating empty-string placeholders.
- [x] Rework `search_libs_for_subs()` to avoid cloning whole maps and duplicate
  index vectors.
- [x] Use sets for library comparison instead of repeated `Vec::contains()`
  scans.
- [x] Iterate directly over unique diffusion pairs rather than entering all
  `n * n` combinations and skipping half.
- [x] Capture immutable transport parameters once in closures instead of cloning
  complete parameter structs on each evaluation.
- [ ] Keep symbolic expressions and closures shared where their ownership model
  permits it.

Required measurements/tests:

- [ ] Criterion benchmarks for cold/warm repository load, 100/1000-substance
  lookup, numeric/closure/symbolic map construction, and diffusion-pair setup.
- [ ] Allocation counts or profiler evidence accompany clone-removal changes.
- [ ] Behavioural equivalence tests protect every optimized path.

## P2.2 Make computation pure by default

- [ ] Remove unconditional `println!`, `print!`, and `Table::printstd()` calls
  from search, parsing, calculation, and network methods.
  Remaining core output includes NASA/NIST fitting quality tables, CEA fitting
  coefficients/`R²`, and a Gibbs calculation diagnostic in `dG_dS.rs`.
- [x] Return structured reports/tables/data from core methods.
- [x] Keep explicit `print_*` adapters or caller-owned logging for interactive
  workflows.
- [x] Ensure pretty-print methods never mutate or rebuild calculation state.

Required tests:

- [x] Pure report methods return deterministic data suitable for CLI, GUI, and
  library callers.
- [x] Printing adapters format the same report without changing state.

## P2.3 Normalize the public API and naming

- [x] Introduce typed library identifiers/capabilities and centralize aliases;
  remove divergent string matching and the `Aramco_transpot` typo.
- [x] Separate property identity from representation.
  - Prefer `PropertyKind::Cp` plus numeric/closure/symbolic accessors over
    variants such as `Cp`, `Cp_fun`, and `Cp_sym` in one enum.
- [x] Return `Option<&T>` rather than `Option<&Box<T>>`.
- [x] Accept `AsRef<Path>` for file APIs and typed units internally.
- [x] Correct API typos such as `hasmap`, `explicit_search_insructions`, and
  `set_explicis_searh_instructions`, with temporary deprecation aliases only if
  downstream compatibility requires them.
- [x] Preserve conventional scientific notation (`T`, `P`, `Cp`, `Lambda`)
  where it improves domain readability; naming cleanup must not erase accepted
  scientific symbols.
- [x] Define a narrow facade for common workflows: configure query, resolve
  records, calculate requested properties, inspect a typed report.

Required tests:

- [x] Compile-oriented API examples cover the primary workflows.
- [x] Deprecated aliases, if retained, forward exactly to canonical methods.

## P2.4 Validate fitting and representation equivalence

- [ ] Validate fitting input lengths, finiteness, interval ordering, rank, and
  solver output before publishing fitted coefficients.
  NASA now follows this contract; NIST and CEA LM branches still need equivalent
  typed solver-result validation.
- [x] Test continuity and boundary ownership between adjacent temperature
  intervals.
- [x] Decide whether adjacent interval fitting is allowed to drift to a
  different coefficient vector when both intervals start from the same
  coefficients, or whether exact coefficient preservation is part of the
  contract. Identical adjacent NASA-7 vectors are now an idempotent case:
  all adjacent fitting paths preserve the vector exactly and only calculate
  quality diagnostics. For genuinely different ranges, coefficient drift is
  allowed and the contract is expressed in property-space fitting errors.
- [x] Strengthen NASA adjacent LM fitting and regression coverage.
  Sequential LM and stacked simultaneous LM are checked against real NASA
  records with independent `Cp`, `dH`, and `dS` error bounds. Stacked residual
  blocks make the simultaneous LM formulation more direct and idiomatic; the
  pre-existing property-wise reports remain the final acceptance check.
  Unsuccessful RustedSciThe fits become typed errors.
  Fitting acceptance is based on `reports_for_coefficients` comparing raw
  property samples with the published coefficient vector; finite coefficients
  or a good aggregate LM `R²` are necessary diagnostics, not quality proof.
- [x] Define one shared property-space report contract for NASA and NIST.
  NASA and NIST now use the same shared helper for property-wise fitting
  diagnostics and finite/length validation. CEA still needs the same treatment
  before the entire thermodynamic fitting surface is fully unified.
  - Reject mismatched/empty sample vectors without `assert_eq!`.
  - Use one documented relative-error formula, including a near-zero reference
    policy, for maximum error and significant-point classification.
  - Make the documented significant residual threshold match the implementation.
  - Correct the NASA `diff / y^2` significant-residual calculation and clarify
    whether `l2_norm` is absolute, relative, or RMS.
- [ ] Route NIST and CEA fitted outputs through the same property-wise
  acceptance principle already used by NASA before publishing coefficients.
  - NIST now has a shared `validate_property_reports` helper plus regression
    tests for the acceptance boundary, but the helper is still not wired as a
    hard publish gate because several legacy fitting fixtures intentionally
    expect the older permissive behavior.
  - CEA still needs the same property-wise report layer.
- [ ] Cross-check numeric, closure, and symbolic implementations over the same
  valid domain for NASA, NIST, CEA, TRANSPORT, and diffusion formulas.
- [ ] Define tolerances from the formulas/data precision rather than weakening
  assertions after failures.

## P2.5 Scenario presets and substance write workflows

- [x] Add a thermodynamics `prelude` module for the core read-only workflow
  imports that preset code and application code reuse most often.
- [x] Add typed presets over `SubsData` for the common application scenarios:
  - `StirredTank` should request only thermochemistry (`Cp`, `dH`, `dS`, `M`).
  - `PlugFlowGas` should request thermochemistry plus diffusion and thermal
    conductivity.
  - `PlugFlowGas` variants that do not need diffusion should request
    thermochemistry plus thermal conductivity only.
  - `PlugFlowCond` should request thermochemistry plus thermal conductivity
    without diffusion or viscosity.
  - Each preset should also carry library preference, temperature interval, and
    any other required lookup policy in one structured object.
  - The first typed preset scaffold now lives in `User_substances_presets.rs`
    and reuses the existing `SubsData` setters plus the read-only prelude.
  - The preset layer now has a typed `ThermoLookupPolicy` and
    `ThermoOutputPolicy`, plus preview/summary helpers for GUI and examples.
- [x] Define the contract for fitting drift between adjacent ranges before more
  fixtures are rewritten around the current implementation. Exact preservation
  applies to identical input vectors; approximation of different vectors is
  validated by `Cp`, `dH`, and `dS` errors rather than coefficient proximity.
- [x] Design a typed JSON write workflow for new substances.
  - Write the substance in the correct format-specific library representation
    (`NASA`, `NIST`, `Aramco`, etc.).
  - Update the address library that stores library-substance pairs.
  - Update the element-composition library.
  - [x] Append a machine-readable journal entry that records which canonical
    names were written and to which files.
  - [x] Add an internal crate-only upsert/drop path that updates one substance
    in one library and syncs the derived address and element catalogs.
- [x] Add roundtrip tests for the write path and for the preset constructors.
- [x] Complete the thermodynamics prelude for preset-oriented application code
  by re-exporting `ThermoLookupPolicy`, `ThermoOutputPolicy`,
  `ThermoPresetPreviewReport`, and `ThermoRequestedData`.

# P3: Test Infrastructure, Documentation, and Cleanup

## P3.1 Build a deterministic test matrix

- [ ] Add property-based/fuzz tests for malformed JSON, NASA composition text,
  NIST interval shapes, units, and finite-domain inputs.
- [x] Add repository story tests using temporary fixture directories.
- [x] Add an end-to-end report regression covering explicit lookup and
  priority fallback.
- [ ] Add end-to-end aggregator tests for explicit lookup, priority fallback,
  mixed thermo/transport sources, invalidation, partial batches, and export.
- [ ] Replace tests that only print values or discard errors with assertions on
  contracts and invariants.
- [x] Keep live network tests outside the default test suite.
  Replace live-data fitting fixtures in `NISTdata_tests.rs` and
  `NISTdata_fitting.rs` with checked-in/local JSON records; retain only a small
  explicitly ignored live smoke-test group.

## P3.2 Document supported contracts and units

- [x] Document every supported library format, required fields, unit system,
  valid temperature range, and transport prerequisites.
- [x] Add examples for local-library lookup, explicit mixed-source lookup,
  numeric values, closures, symbolic expressions, and diffusion coefficients.
- [x] Document best-effort versus fail-fast batch semantics.
- [x] Document repository loading/reload and offline NIST behaviour.

## P3.3 Remove dead or misleading module surface

- [x] Either implement real facades in `thermo_properties_api.rs` and
  `transport_properties_api.rs` or remove the empty public modules.
- [x] Remove `User_substances3.rs` if it remains an empty implementation block.
- [x] Remove unused imports and stale commented code from the thermo surface.
- [x] Remove duplicated aliases from the thermodynamics surface.
- [x] Remove misleading zero/default calculator state.
- [x] Move filenames and module names toward Rust conventions where this can be
  done without obscuring domain notation.
  - Conventional snake_case module aliases are now available for new code:
    `chem_equilibrium`, `db_handlers`, `user_substances`,
    `user_substances2`, `user_phase_or_solution`, and
    `user_phase_or_solution2`.
  - Legacy module names remain as compatibility shims so the broad public
    surface does not break in one shot.

# Linear Execution Plan

To avoid repeatedly rewriting the same lookup and cache functions, execute the
work in this order:

1. **Pass 1: parsing and numeric safety**
   Complete P0.1 and P0.2. Do not redesign `SubsData` yet.
2. **Pass 2: lookup correctness and batch truthfulness**
   Complete P0.3 and P0.4 with deterministic per-property reports.
3. **Pass 3: mutation invariants**
   Complete P0.5 and its full mutation-path test matrix.
4. **Pass 4: repository and canonical state**
   Complete P1.1 and P1.2 together; this is the only planned structural rewrite
   of `SubsData`.
5. **Pass 5: clone, errors, and network boundary**
   Complete P1.3-P1.5 after the canonical ownership model exists.
6. **Pass 6: measured optimization and API cleanup**
   Complete P2 only against the stabilized contracts and benchmark baseline.
7. **Pass 7: documentation and final debt removal**
   Complete P3 and remove compatibility scaffolding that is no longer needed.

## Residual Closeout Plan

This plan supersedes repeated broad re-audits. Work linearly and do not reopen
the completed repository or canonical `SubsData` redesign unless a regression
demonstrates a concrete contract violation.

1. **Closeout pass 1: panic containment**
   Complete P0.6 in `dG_dS.rs`, `User_substances2.rs`,
   `NISTdata_fitting.rs`, and `CEAdata.rs`, with state-preservation error tests.
2. **Closeout pass 2: fitting and report contract**
   Normalize NASA/NIST property reports, then apply typed LM result validation
   and property-wise acceptance to NIST and CEA.
3. **Closeout pass 3: deterministic tests and pure core**
   Replace default live NIST tests with fixtures and move all fitting/Gibbs
   output behind explicit `print_*` adapters.
4. **Closeout pass 4: API/test completeness**
   Finish typed error boundaries, best-effort batch reporting, compile-time
   cache privacy tests, provenance-through-export coverage, and prelude exports.
5. **Optional measured pass: performance only**
   Add benchmarks and make further clone/allocation changes only when profiling
   demonstrates a material cost.

## Final Completion Criteria

- [ ] No scoped production path panics on user, file, network, or library data.
- [ ] Every successful numeric result is finite and valid for its declared
  physical contract.
- [ ] Lookup is deterministic, property-specific, and preserves provenance.
- [ ] Partial failures are visible in returned typed reports.
- [ ] Input mutation cannot expose stale derived data.
- [ ] Library files are loaded once per explicit repository lifecycle, not once
  per aggregator.
- [ ] Clone behaviour is semantically faithful or deliberately unavailable.
- [ ] Core computations do not print or mutate unrelated state.
- [ ] Default tests are deterministic and offline.
- [ ] Benchmarks demonstrate that performance changes address measured costs.

# Dedicated Phase/Solution Architecture Audit

## Scope and current assessment

This audit covers `User_PhaseOrSolution.rs`, `User_PhaseOrSolution2.rs`, and
their direct equilibrium-facing contract. It does not propose a rewrite of the
equilibrium solvers themselves, but the adapter between the phase model and the
solvers must migrate with the data model.

The existing seven focused tests pass, but they mainly prove that the current
surface can be called. They do not protect component identity, deterministic
matrix/vector alignment, clone equivalence, failure atomicity, or a substance
appearing in more than one phase. Those are the central correctness contracts
for this module.

The intended canonical model is:

- `PhaseSystemSpec`: validated user input and lookup policy.
- `PhaseSpec`: typed phase id, explicit thermodynamic/activity model, and an
  ordered component list.
- `ResolvedPhaseSystem`: resolved `SubsData` records plus one immutable
  `SystemLayout` shared by every numeric, symbolic, and solver-facing path.
- `SystemLayout`: the sole ordered list of `(phase_id, substance_id)` component
  identities, phase ranges, and reverse indices used by vectors and matrices.
- Typed composition and result snapshots instead of nested
  `HashMap<Option<String>, (Option<f64>, Option<Vec<_>>)>` values.

`OnePhase` should become a convenience constructor or thin newtype over the
same canonical phase-system engine. It should not maintain a second copy of all
thermodynamic and equation-building algorithms.

## Phase/Solution P0: correctness before redesign

### P0.1 Preserve phase-component identity and deterministic order

- [x] Introduce typed `PhaseId` and `PhaseComponentId` (phase plus substance).
  A substance present in gas and liquid phases is two solver components even
  though both components share the same molecular composition.
- [x] Replace the canonical `HashMap<Option<String>, SubsData>` iteration order
  with an ordered phase/component layout. A `Vec` plus lookup indices is enough;
  no new container dependency is required.
- [x] Make `SubstancesContainer::get_all_substances()` distinguish between a
  unique catalog view and an ordered solver-component view. The current
  deduplication is invalid for multiphase equilibrium.
  - `get_all_substances()` now stays a deduplicated catalog view, while
    `get_ordered_component_labels()` and `get_ordered_substances()` expose the
    solver-facing ordered layout without flattening phase identity.
- [x] Rewrite both `indexed_moles_variables()` implementations from the shared
  layout.
- [x] Fix the current multiphase name collision produced by `format!("n{}", i + j)`.
  - [x] Fix the assignment that currently writes
    `self.vec_of_n_vars = self.vec_of_n_vars.clone()` instead of publishing the
    newly constructed vector.
  - [x] Make variable names stable under different `HashMap` insertion orders.
- [x] Rewrite both `create_full_map_of_mole_numbers()` implementations so mole
  vectors follow the declared component order, never `HashMap::values()` order.
- [x] Build the element-composition matrix, Gibbs/entropy vectors, symbolic
  variables, initial composition, and result labels from the same
  `SystemLayout`.
  - Matrix assembly and result-label generation now follow the canonical
    ordered layout instead of separate ad hoc phase iteration paths.

Required tests:

- [x] The same substance in two phases produces two distinct component ids,
  variables, matrix rows, and phase amounts.
- [x] Shuffling phase-map insertion order produces an identical layout and
  identical equations.
- [x] Mole vectors have the exact declared component order, not merely the
  expected length.
- [x] All generated variable names are globally unique for multiple uneven
  phase sizes.

### P0.2 Remove false-success clone and panic paths

- [x] Remove the manual `Clone` implementations that replace working `dG_fun`
  and `dS_fun` closures with functions returning `0.0`.
- [x] Store shareable compiled functions as `Arc<dyn Fn...>`, or
  remove `Clone` and provide an explicit fallible rebuild. A clone must never
  be numerically different while looking ready.
- [x] Make `SubstancePhaseMapping::new()` fallible instead of asserting on
  user-provided vector lengths.
- [x] Remove production `unwrap()` calls from:
  - `calculate_elem_composition_and_molar_mass()`;
  - `indexed_moles_variables()`;
  - `calculate_Lagrange_equations_sym()` and
    `calculate_Lagrange_equations_fun2()`;
  - both `calculate_Lagrange_equations_fun()` closure constructors;
  - both `create_full_map_of_mole_numbers()` implementations.
- [x] Validate all matrix dimensions, component lookups, and closure inputs at
  construction time. Solver closures cannot return typed errors, so an invalid
  closure must never be published.
  - `calculate_Lagrange_equations_sym()` and `calculate_Lagrange_equations_fun()`
    now reject mismatched component counts and missing Gibbs data before a
    closure is published; the runtime closure guards against wrong lambda/vector
    lengths by returning an empty result instead of panicking.
- [x] Remove direct `println!` diagnostics from equation construction and
  retain the complete context in typed errors.

Required tests:

- [x] Original and clone return equivalent numeric, symbolic, and compiled
  Gibbs/entropy values.
- [ ] Missing Gibbs data, wrong matrix dimensions, missing phases, and wrong
  lambda/composition lengths return typed errors without partial publication.
  - [x] The pure Lagrange adapter rejects missing Gibbs functions and invalid
    matrix dimensions before it returns a closure; malformed runtime solver
    vectors return the established empty-result signal rather than panicking.
  - [x] Symbolic wrong-matrix and numeric missing-Gibbs regressions now also
    prove that the existing cache family and revision remain unchanged after
    rejection.
- [ ] No public constructor or mutation path panics on malformed input.
  - [x] Public cache views no longer expose an indexing operator that panics for
    an absent phase; callers use the explicit optional `get()` contract.
  - [x] Public typed constructors and pure Gibbs/entropy evaluation facades are
    covered by an unwind regression test for malformed input. Internal legacy
    mutation adapters remain a separate migration concern.

### P0.3 Validate compositions and phase definitions

- [x] Replace nested `Option` tuples with a typed `PhaseComposition` containing
  an ordered amount vector and an explicit total amount.
  - [x] Numeric one-phase and multi-phase read paths now use a typed
    `PhaseComposition` wrapper instead of passing the raw tuple around as the
    primary internal representation.
  - [x] The typed composition now has a fallible numeric constructor, so
    invalid phase totals and component vectors fail before entering the
    calculation paths.
- [x] Require finite non-negative component amounts and a finite positive phase
  total where logarithmic activity terms are evaluated.
  - [x] Numeric phase-composition paths now validate finite, non-negative,
    summed phase totals before Gibbs/entropy evaluation.
- [x] Validate finite positive temperature and pressure together with an exact
  phase/component composition request before evaluating any phase.
  - `ThermoEvaluationConditions` now makes the physical inputs explicit, and
    `evaluate_gibbs`/`evaluate_entropy` reject missing, extra, or wrongly sized
    phase compositions before a lower-level calculator is invoked.
  - The read-only evaluation paths now run against local `SubsData` snapshots,
    so both successful and failed evaluations leave the wrapper and its
    published result caches unchanged.
- [x] Define and validate whether `Np` must equal the sum of component amounts;
  do not silently accept contradictory values.
  - [x] The numeric composition constructor now enforces total/component sum
    agreement instead of letting contradictory values reach the solver.
- [x] Reject unknown, missing, and extra phases in numeric requests; reject
  empty/duplicate typed phase specifications and non-finite `T`/`P` before a
  phase calculator runs.
- [x] Validate finite positive reference temperature `Tm` at the symbolic and
  numeric Lagrange adapter boundary.
- [x] Stop silently assigning `Phases::Gas` when a named multiphase entry is
  missing from `physiscal_nature_of_phase_or_solution`.
  - `SubstancesContainer::MultiPhase` now requires a complete physical-state
    map. Only the explicit `SinglePhase` compatibility form defaults to an
    ideal-gas `PhaseSpec`.
- [x] Make all multi-phase operations transactional: failure in a later phase
  must not leave earlier `SubsData` instances or derived caches updated.
  - [x] Phase-local preparation and symbolic/closure construction now use one
    commit-on-success payload barrier. A failure after any earlier temporary
    phase update leaves the visible payload map, derived property caches, and
    cache shape unchanged.
  - [x] Numeric Gibbs/entropy publication now validates the request, canonical
    phase/component cache shape, and finite values before replacing any cache
    family or context; rejected publication retains the previous revision and
    cache state.
  - [x] Symbolic and closure builders now validate exact canonical
    phase/component coverage before commit; raw cache-map insertion is
    restricted to test fixtures and cannot create a production bundle for an
    unknown phase.
  - [x] Typed multi-phase resolution now installs payloads, canonical layout,
    empty cache bundles, and lookup provenance through one consuming
    `ResolvedPhaseSystem -> PhaseSystem` transition rather than exposing a
    payload-only intermediate state.
  - [x] The single-phase factory now consumes the same resolved payload through
    one validated install transition; it no longer publishes data, phase spec,
    and provenance through separate mutations.
  - [x] `OnePhase` rejects a multi-phase resolved payload before any state is
    published, preventing accidental loss of phase-qualified components.
  - [x] Symbolic Gibbs substitutions now publish a new cache revision, so a
    snapshot cannot look current after its expressions have changed.
  - [x] `calculate_lagrange_equations_fun2()` now stages Gibbs closure
    construction, validates equation assembly, and only then publishes the
    refreshed payload/cache pair. A matrix-layout failure preserves the prior
    cache shape and revision.

Required tests:

- [x] Negative, NaN, infinite, inconsistent-total, missing-phase, extra-phase,
  and duplicate-component cases fail with actionable typed errors.
  - [x] `Tm` is covered directly for zero, `NaN`, and both infinities in the
    symbolic and numeric Lagrange adapters.
- [x] A failure in the second phase leaves the complete system state unchanged.
- [ ] Zero mole fractions follow one explicitly documented boundary policy and
  never accidentally publish NaN.
  - [x] Typed numeric Gibbs/entropy requests reject a zero component amount
    whenever the raw calculator would apply its ideal-gas `ln(n_i / Np)`
    correction. The check runs before calculator preparation; zero condensed
    components remain valid, and absent phase metadata continues to mean gas
    to match the legacy calculator contract.
  - [ ] Defer migration of the remaining raw Gibbs/entropy closures to a
    fallible or guarded evaluation contract until the ChemEquilibrium facade
    is migrated. A bare `Fn(...) -> f64` cannot surface the same typed domain
    error at invocation time, so it remains a compatibility-only API and must
    not become the primary public numeric evaluation path.

### Deferred until the ChemEquilibrium migration

These items depend on changing equilibrium consumers, not on further changes
to the phase-data engine. Keep them out of independent phase-system passes.

- [ ] Migrate equilibrium consumers from `ThermodynamicsCalculatorTrait` to
  `PhaseLayoutAccess`, `PhasePropertyEvaluator`,
  `PhaseSymbolicPropertyBuilder`, and `PhaseEquilibriumAssembly`.
- [ ] Replace equilibrium closure input `Option<Vec<f64>>` with a borrowed
  typed composition view aligned to `SystemLayout`.
- [ ] Migrate legacy raw Gibbs/entropy `Fn(...) -> f64` closures together with
  the equilibrium facade. Do not add another parallel phase-cache format
  before this migration: equilibrium consumers still require the infallible
  callback signature.
  - Return to this item only after the equilibrium facade accepts the narrow
    typed composition API and can propagate `SubsDataResult` at callback
    invocation time.
- [ ] Remove `CustomSubstance`/`enum_dispatch` once no equilibrium-facing
  consumer requires its broad legacy surface.
- [ ] Add equilibrium integration story-tests for `PhaseEvaluationRequest`,
  cache revision invalidation, and the same substance represented in two
  phases.

## Phase/Solution P1: canonical domain model

### P1.1 Separate phase identity, physical state, and activity model

- [x] Add a typed `PhaseSpec` with `PhaseId`, ordered components, and a
  `PhaseModel`.
  - `SubstanceSystemSpec` now stores ordered `PhaseSpec` values as its canonical
    input. `SubstancesContainer` and `with_phase_natures()` are retained only as
    a fallible adapter for existing callers.
- [x] Keep physical state (`Gas`, `Liquid`, `Solid`, `Condensed`) separate from
  the thermodynamic/activity model.
  - [x] Route physical state into local record selection through one typed
    `ThermoRecordQuery`, rather than requiring callers to encode a library key
    such as `Fe(L)` by hand.
  - [x] Make suffix handling library-aware: `NASA_cond`/NIST interpret their
    documented `(g)/(l)/(s)` and NASA-condensed allotrope conventions, while
    `nuig_thermo` parenthesized chemistry labels remain literal identifiers.
  - [x] Preserve the selected source record key in `SearchResult` provenance
    and reject ambiguous state requests (for example solid `Fe` with several
    NASA-condensed polymorphs) rather than selecting an arbitrary record.
  - [x] Keep `PhaseSpec` physical state as phase/NIST context, not an implicit
    local-library filter. State-aware lookup is an explicit `SubsData` opt-in
    and its bulk setter remains transactional when a caller chooses it.
  - [ ] Extend the catalog format with explicit state metadata when upstream
    libraries provide it. Key conventions remain a compatibility bridge, not a
    universal chemical-identity parser.
- [x] Define the supported models explicitly, starting with the models the code
  really implements: `IdealGas` and `PureCondensed`.
- [x] Do not call a container a "solution" until an explicit solution model is
  selected. The current `PhaseModel` deliberately has no solution variant:
  `Phases::Liquid` or `Phases::Solid` resolve only as `PureCondensed`, not as an
  implicit ideal or non-ideal solution.
- [x] Require legacy multi-phase state maps to match declared phase names
  exactly: both missing states and states for unknown phases are typed errors.
- [ ] Decide the first real solution contract (`IdealSolution` and its standard
  state, or a deliberately unsupported typed variant) before exposing it in
  the public builder.

### P1.2 Use one engine for one or many phases

- [x] Replace duplicated `OnePhase` and `PhaseOrSolution` algorithm/state fields
  with one `PhaseSystem` implementation.
  - Both public facades are now thin views over the same crate-private engine.
    The only remaining duplicated surface is deliberately compatibility-facing
    and is tracked separately under the deferred `CustomSubstance` and
    equilibrium migration items below.
- [x] Move the canonical typed metadata (`PhaseSpec` and resolved
  `SystemLayout`) into `PhaseSystem`, so both wrappers now read one common
  declaration/layout owner after typed resolution.
  - [x] Move canonical phase-keyed `SubsData` storage into `PhaseSystem`.
    `OnePhase` no longer owns a separate `SubsData` field and `PhaseOrSolution`
    no longer owns a separate raw phase map; both facades now read the same
    engine-owned payload through read-only views.
  - [x] Replace direct consumer/test access with `subs_data_view()` and
    `phase_data_view()`. Crate-private replacement and transition helpers
    invalidate metadata and derived caches after the final payload change.
  - [x] Keep `PhaseSystem` crate-private so external code cannot mutate cache
    maps or raw payloads behind the facade invariants.
  - [x] Route pure Gibbs and entropy evaluation through `PhaseSystem`, so the
    one-phase facade no longer carries a second property-evaluation path.
  - [x] Move shared payload preparation, NIST fallback, numeric result
    publication, symbolic/closure cache construction, and symbolic `T`/`P`
    substitution into `PhaseSystem`. `OnePhase` and `PhaseOrSolution` now
    delegate these operations instead of maintaining competing implementations.
  - [x] Centralize layout-derived element matrices, symbolic mole variables,
    and sparse mole-number normalization in `PhaseSystem`.
    - Element/molar-mass assembly now evaluates cloned phase payloads for both
      facades, instead of mutating only the one-phase payload cache.
    - Both facades now use the same `Np{phase}` / `n{phase}_{component}` names
      and generate solver vectors exclusively in `SystemLayout` order.
    - The one-phase path no longer derives vectors from `HashMap::values()`.
  - [x] Introduce `MoleNumberSnapshot` as the typed result of sparse mole
    normalization. It carries the `SystemLayout`, normalized named amounts,
    ordered phase vectors, and cross-phase substance totals; the former tuple
    of maps now exists only as an explicit compatibility projection.
  - [x] Introduce `PhaseEvaluationRequest` as the atomic numeric-evaluation
    boundary. Temperature, pressure, and phase compositions now travel
    together through the pure evaluator and numeric cache provenance, instead
    of being independently threaded through parallel method parameters.
  - [x] Replace anonymous symbolic mole-variable bookkeeping with the named
    `PhaseSymbolicLayout` snapshot. The former four-part tuple is now projected
    only for the broad compatibility trait; engine code and tests consume named
    vectors/maps generated from the canonical `SystemLayout`.
  - [x] Keep `PhaseSystem` as the sole owner of `SystemLayout`. The symbolic
    snapshot stores only derived variable data; it no longer duplicates a
    second layout/provenance copy that could diverge from `resolved_layout`.
- [x] Introduce a shared read-only `PhaseDataView` and route pure numeric
  Gibbs/entropy evaluation through it for both wrappers.
  - The bridge preserves the same ordered layout for one and many phases and
    evaluates only cloned `SubsData` payloads, so a read-only evaluation cannot
    publish a partial result.
  - Crate-private replacement/mutation compensators now invalidate phase specs,
    layouts, and all caches before exposing legacy raw payload edits. Their
    regression tests cover both one- and multi-phase paths.
- [x] Retain `OnePhase` as a one-entry facade over `PhaseSystem` storage and
  cache state. Remaining method-level duplication is tracked by the parent
  engine-consolidation item above.
- [ ] Reduce or remove `CustomSubstance`/`enum_dispatch` once both variants use
  the same engine.
  - [x] `CustomSubstance` now implements the narrow preparation, symbolic, and
    equilibrium roles by delegation. New consumers no longer need the broad
    enum-dispatch trait; removing the enum itself remains deferred until the
    equilibrium migration.
- [ ] Split the oversized `ThermodynamicsCalculatorTrait` into narrow roles:
  resolved data access, phase-property evaluation, layout/composition access,
  and equilibrium equation assembly.
  - [x] Introduce `PhaseLayoutAccess` and `PhasePropertyEvaluator` for both
    facades. `CustomSubstance::system_layout()` now consumes the layout role,
    while the old enum-dispatch trait remains a temporary compatibility adapter
    for legacy equilibrium code.
  - [x] Introduce `PhaseDataPreparation`, `PhaseSymbolicPropertyBuilder`, and
    `PhaseEquilibriumAssembly` for both facades. New solver-facing code can
    depend on one capability at a time; `ThermodynamicsCalculatorTrait` stays
    intact only as a compatibility adapter until its legacy consumers migrate.
    Its preparation, symbolic, and Lagrange methods now delegate through those
    narrow roles instead of maintaining a second direct `PhaseSystem` route.
- [ ] Move equilibrium-only Lagrange equation builders out of the data
  acquisition facade and into an equilibrium adapter module.
  - [x] Symbolic Gibbs-minimisation assembly now lives in
    `phase_equilibrium_adapter.rs`. It is a pure operation over
    `SystemLayout`, a borrowed Gibbs view, and the component-by-element matrix.
    The adapter uses global component indices, preventing the former second
    phase matrix-column reuse bug.
  - [x] Numeric closure assembly now also lives in
    `phase_equilibrium_adapter.rs` and uses the global `SystemLayout` index for
    each component. This fixes the former multi-phase bug where every phase
    reused matrix columns from zero.
  - [x] Move the remaining facade orchestration into `PhaseSystem`: symbolic
    builder dispatch, closure-map conversion, and the `calculate_*_fun2`
    workflow now share one engine path for one and many phases. The adapter
    remains pure and receives only layout, property views, matrices, and
    explicit legacy closure maps.
  - [ ] Replace the compatibility closure input `Option<Vec<f64>>` with a
    borrowed typed composition view once the equilibrium solver facade itself
    can accept that narrower contract.
- [x] Make read-only methods take `&self`; methods such as
  `get_all_substances()` and `extract_SubstancesContainer()` do not need
  `&mut self`.

### P1.3 Replace parallel public caches with typed snapshots

- [x] Make phase-system fields private and expose invariant-preserving builders,
  mutation methods, and borrowed views.
  - `OnePhase` and `PhaseOrSolution` cache fields are now private or
    crate-private and are read through explicit accessor/view methods.
    Layout revision is also hidden behind accessors, so callers no longer need
    direct field access for the visible contract. Test code now reaches cache
    state only through debug-only mutation hooks or read-only views.
  - Single-phase raw mutable cache accessors are now crate-private and are
    reachable only through explicitly named debug hooks; external callers must
    use validated publication APIs.
  - Bulk raw `replace_dG*`/`replace_dS*` fixture helpers are now compiled only
    for tests. Production builds retain no unvalidated cache-replacement path.
- [x] Replace parallel `dG`, `dG_fun`, `dG_sym`, `dS`, `dS_fun`, and `dS_sym`
  maps with typed result snapshots aligned to `SystemLayout` and carrying their
  thermodynamic conditions.
  - [x] PhaseSystem now stores per-phase Gibbs and entropy property families,
    each grouping numeric, closure, and symbolic representations. The result
    families stay aligned inside one phase-scoped record rather than six
    independently named maps.
  - [x] Snapshot construction is now centralized through shared typed helpers
    for single-phase and multi-phase views, and the public facade now exposes a
    typed `ThermoResultSnapshot` aligned to `SystemLayout` and explicit
    temperature/pressure metadata. The remaining work is to reduce the
    underlying parallel caches themselves, not just their public packaging.
  - [x] Cache write paths now also go through shared revision-aware helpers
    instead of repeating raw `field = value; bump_revision()` sequences.
  - [x] Numeric Gibbs and entropy caches now retain a typed request context
    (temperature, pressure, phase compositions, and configuration revision).
    Consumers can ask whether a cache is valid for an exact new request instead
    of trusting the mere presence of values.
  - [x] The symbolic Lagrange read path now borrows the cache view instead of
    cloning the full `dG_sym` map before evaluation.
  - [x] The `calculate_Lagrange_equations_fun2()` helpers now consume the
    canonical cached Gibbs-function view instead of rebuilding function maps
    from `SubsData` again.
- [x] Prefer pure `evaluate_*`/`build_*` methods that return results, and tag
  published cache products with the inputs/layout revision that makes them
  valid.
  - [x] `OnePhase` and `PhaseOrSolution` now expose `evaluate_gibbs()` and
    `evaluate_entropy()` paths which return values without publishing a
    `PhaseSystem` cache. The legacy calculation methods explicitly publish the
    returned values afterwards.
  - [x] `SubsData` derived builders now run on a cloned working copy, so the
    public one-phase Gibbs/entropy paths no longer publish cache state while
    computing a result. The remaining `&mut self` cache-building methods are
    now explicit publication paths rather than hidden side effects.
- [x] Remove or rename `get_dG_sym_mut()`: it currently returns an owned clone,
  not mutable access.
- [x] Return borrowed views from `get_dG_sym()`/`get_dG()` instead of cloning
  complete nested maps on every read.
  - Key read paths now use borrowed nested-map views, and the old cloning-based
    accessors remain only as compatibility snapshots.
  - `PhaseThermoPropertyView` now also exposes `iter_ref()` plus explicit
    `phase_count` and `populated_phase_count`, so borrowed consumers do not
    need a boxed owned-key iterator merely to inspect cache state.
  - `NestedPhaseCacheView` now mirrors the same borrowed traversal contract
    through `NestedPhaseCacheIter`, including the single-phase `None` key.
    Both cache views now distinguish canonical phase count from populated
    property-cache count explicitly.
- [x] Invalidate layout, variables, property snapshots, and equations through
  one revision boundary whenever phases, components, lookup policy, or resolved
  data change.
  - `PhaseOrSolution` and `OnePhase` now carry `layout_revision` counters, and
    the main layout/property-building paths bump them on mutation while
    read-only Lagrange builders stay stable.
  - The revision boundary now also recreates cache bundles from canonical
    payload keys, so borrowed cache views cannot reveal values from the prior
    configuration merely because their request context was cleared.

### P1.4 Separate specification, resolution, and computation

 - [x] Replace `SubstanceSystemFactory::create_system(...many primitive
  arguments...) -> Result<_, String>` with a typed builder and typed error.
  - `SubstanceSystemSpec` and its fluent builder now carry the construction
    inputs, and `resolve()` performs the typed lookup step. The legacy wrapper
    remains only as a temporary compatibility layer.
 - [x] Make system construction validate only the specification; make database
  lookup/NIST access an explicit fallible `resolve()` step.
  - The new builder validates empty/degenerate containers before lookup, and
    resolution is now an explicit typed step rather than implicit inside the
    primitive argument list.
  - Lookup-policy application is centralized so single-phase and multi-phase
    resolution share the same configuration path.
  - [x] `ResolvedPhaseSystem` now atomically carries resolved `SubsData`, the
    corresponding ordered `PhaseSpec` values, and the derived `SystemLayout`
    before a single- or multi-phase facade is selected.
- [x] Reuse a shared `ThermoRepository`/lookup policy across phases instead of
  treating every phase as an independent database lifecycle.
  - Phase resolution can now receive one explicit `Arc<ThermoRepository>`;
    every resolved `SubsData` receives that same immutable catalog while
    keeping its own query state and calculator instances. The default path
    obtains the same shared catalog through `ThermoData::try_default_repository`.
- [ ] Store phase kind once on the phase, then derive the per-substance phase
  view needed by `SubsData`; do not duplicate the same value manually into
  every component as canonical state.
  - [x] `PhaseSpec` now derives the legacy `SubsData::map_of_phases` projection
    in one place. The raw per-substance map remains only because existing
    lower-level property calculations still consume that compatibility shape.
- [x] Keep explicit NIST fallback policy visible in the resolved-system report,
  including per-component provenance.
  - `ResolvedPhaseSystemReport` now records the NIST fallback policy and an
    ordered per-phase `SearchSummaryReport`; each row exposes the selected
    library, priority tier, state, and calculator family without raw cache
    access.

## Phase/Solution P2: API, performance, and module boundaries

- [ ] Replace closure arguments `Option<Vec<f64>>` and `Option<f64>` with
  borrowed typed composition/state views. Avoid cloning the full composition
  vector once per species evaluation.
  - Deferred until ChemEquilibrium accepts a fallible typed composition
    callback. The phase engine already exposes `PhaseComposition`; changing
    only this side would create a second adapter around the same legacy
    closure signature.
- [ ] Transpose and validate matrices once when constructing a numeric closure,
  not on every solver call.
  - Deferred until the equilibrium adapter owns the numeric closure boundary.
    Matrix orientation belongs to Lagrange assembly, not to phase lookup or
    property evaluation.
- [ ] Use dense ordered vectors for solver hot paths and maps only as lookup or
  presentation indices.
  - [x] `ThermoResultSnapshot` now projects numeric Gibbs and entropy caches
    into validated `SystemLayout`-ordered vectors. The remaining migration is
    to make equilibrium consumers accept these vectors directly instead of
    rebuilding map lookups inside their legacy closure path.
- [ ] Share immutable symbolic expressions and compiled functions where their
  ownership permits it.
- [ ] Normalize spelling and names (`calculate`, `mixture`, `coefficients`,
  `phase`, `Np`, and `dG`/`dS` domain notation) without erasing conventional
  scientific symbols.
- [x] Delete the large commented-out equilibrium-constant implementation.
  It was an uncompiled duplicate with no validated contract; any future
  equilibrium-constant behavior must enter through the dedicated equilibrium
  workstream with tests rather than being revived from a stale comment block.
- [x] Split the current files by responsibility: domain model/layout, builder
  and resolution, property evaluation, and equilibrium adapter. Keep public
  re-exports as the stable facade.
  - [x] Extract `phase_domain.rs`: typed phase state/model/specification and
    resolved-system payloads now have one focused home.
  - [x] Extract `phase_evaluation.rs`: numeric evaluation requests, condition
    validation, composition validation, and the pure `SubsData` bridge are no
    longer mixed with facade/cache code.
  - [x] Extract `phase_cache.rs`: symbolic layouts, per-phase cache bundles,
    borrowed cache views, and revision-aligned snapshots are separated from
    `PhaseSystem` orchestration.
  - [x] Extract `phase_factory.rs`: compatibility container adapters, typed
    specification/builder validation, lookup policy, repository resolution,
    and factory errors are separated from `PhaseSystem` orchestration.
  - [x] Extract `phase_facade.rs`: `CustomSubstance` is now a read-only
    compatibility facade over `OnePhase` and `PhaseOrSolution`, with no direct
    ownership of raw payloads or cache state.
  - [x] Extract `phase_interfaces.rs`: narrow consumer capabilities and their
    facade implementations are separated from the engine and from the broad
    legacy calculator trait.
  - [x] Extract `phase_system_facade.rs`: the public `PhaseOrSolution` API,
    pure evaluation entry points, read-only cache views, snapshots, and
    crate-private transition helpers no longer obscure `PhaseSystem` state
    orchestration.
  - [x] Extract `phase_legacy_api.rs`: the deprecated broad
    `ThermodynamicsCalculatorTrait` and its multi-phase compatibility
    implementation now live outside the engine; the disabled duplicate and
    unmaintained equilibrium stubs were removed from the facade module.
  - [x] Extract `phase_moles.rs`: canonical `n_i`/`Np` variable construction
    and sparse-to-ordered mole normalization now live outside lifecycle and
    cache orchestration, while `PhaseSystem` remains their single owner.
  - [x] Extract `phase_elements.rs`: pure component-by-element matrix and
    molar-mass assembly are now isolated from cache lifecycle and operate on
    cloned phase payloads aligned to `SystemLayout`.
  - [x] Extract `phase_operations.rs`: fallible phase-local preparation now
    stages cloned payloads transactionally outside `PhaseSystem`; the engine
    remains responsible only for committing a successful stage and invalidating
    derived state.
  - [x] Extract `phase_property_builders.rs`: symbolic Gibbs/entropy builders
    and compiled closure construction now return staged payloads and cache
    values without owning publication or revision policy.
  - [x] Document `phase_layout.rs` as the canonical solver-order contract,
    including the distinction between a substance and a phase-qualified
    component.
  - [x] Keep the remaining `PhaseSystem` orchestration together. It is the
    single owner of transactional state transitions and revision-aware cache
    publication; extracting it further would recreate the second orchestration
    facade this refactor removed.

## Phase/Solution P3: contract and story-test matrix

- [x] Add numeric/closure/symbolic equivalence tests at the same `T`, `P`, and
  composition for every supported `PhaseModel`.
  - [x] The `IdealGas` Gibbs path is now checked end-to-end against one resolved
    local NASA fixture: pure numeric evaluation, the cached closure, and the
    symbolic expression after `T`, `P`, and composition substitution agree.
  - [x] The `PureCondensed` Gibbs path is checked through a named liquid phase.
    It reuses a NASA source record while proving that numeric, closure, and
    symbolic paths all suppress the ideal-gas concentration correction.
  - [x] Entropy now has the same three-way equivalence coverage for `IdealGas`
    and `PureCondensed`; the latter proves the condensed path remains free of
    ideal-gas concentration terms.
  - [x] The standalone ChemEquilibrium log-moles formulation now has an
    offline fixture proving that the canonical numeric residual, its legacy
    closure adapter, and the symbolic residual generator evaluate identically
    for the same prepared problem and log-mole vector.
- [ ] Add a two-phase story with a shared substance across phases and verify
  element conservation, equation count, variable labels, and phase totals.
  - [x] The typed `spec -> resolve -> layout -> normalize` path now covers H2
    in gas and liquid phases. The story verifies distinct phase-qualified
    component labels, `IdealGas` versus `PureCondensed` models, independent
    mole variables/totals, and the expected cross-phase H2 total. Element and
    equilibrium-equation assertions remain part of the deferred equilibrium
    workstream.
- [ ] Add pure-condensed, ideal-gas, and explicitly unsupported solution-model
  fixtures.
- [x] Add state-transition tests for spec -> resolved -> layout -> evaluated ->
  invalidated/rebuilt.
  - [x] The one-phase lifecycle now has one story-test covering resolution,
    pure evaluation, explicit cache publication, payload mutation that clears
    layout/cache/report state, and clean reconstruction from the original spec.
  - [x] A matching multi-phase story now publishes Gibbs values for gas and
    condensed phases, mutates only the gas payload, and proves that all
    system-wide cache/layout/report state is invalidated before rebuilding.
- [x] Add compile-time/API tests proving callers cannot mutate phase/component
  collections or derived snapshots directly.
  - [x] The public `PhaseOrSolution` documentation now contains a
    `compile_fail` boundary test: an external caller cannot access the private
    `PhaseSystem` engine and therefore cannot bypass revision-aware transition
    methods or cache invalidation.
- [x] Add local-repository story tests for mixed library provenance without
  network access.
  - [x] A temporary two-library NASA fixture now resolves gas H2 from
    `NASA_gas` and liquid H2O from `NASA_cond` through one injected repository.
    The phase facade report proves both selected libraries while NIST fallback
    remains disabled.
- [ ] Add property tests for phase insertion order and composition-map insertion
  order invariance.
  - [x] Deterministic regression coverage now proves that permuting both raw
    phase records and `PhaseSpec` input leaves the resolved solver layout and
    phase order unchanged. Randomized/property coverage remains open.
  - [x] An exhaustive three-phase matrix now covers all six `PhaseOrSolution`
    insertion orders together with reverse component-map insertion, and proves
    stable component labels, symbolic mole variables, and normalized vectors.
- [ ] Replace length/contains-only assertions with exact component order,
  equation identity, numeric value, and typed-error assertions.

### Phase/Solution technical stage close-out

The independent engine work is complete for this stage: canonical layouts,
transactional payload/cache publication, pure numeric evaluation, provenance,
repository sharing, module boundaries, and deterministic lifecycle stories are
covered. Do not reopen this layer for cosmetic refactors alone.

### Deferred phase backlog

These are explicit future tasks, not unfinished work hidden inside the current
phase-engine implementation. Do not start them without entering the named
workstream and agreeing on its contract.

#### P1: physical-model decisions

- [ ] Define the first supported solution model before adding a `Solution`
  `PhaseModel` variant: ideal/non-ideal activity law, standard state,
  composition domain, and the meaning of pressure for that phase.
- [ ] Add fixtures for the chosen solution model and a typed rejection fixture
  for every still-unsupported model. The existing `IdealGas` and
  `PureCondensed` matrix is complete; this item is about new physics only.

#### P1: ChemEquilibrium migration

- [x] Give `ReactionExtentError` and `SolveError` stable `Display` and
  `std::error::Error` contracts. Nested `SubsDataError` and `SolveError`
  sources now remain discoverable through `source()` instead of forcing GUI or
  workflow code to parse `Debug` output.
- [ ] Remove the deprecated threshold-based `equilibrium_logmole_residual2` /
  `equilibrium_logmole_jacobian2` pair after any external experimental users
  migrate to the canonical validated log-moles formulation. The core solver
  and all new tests already avoid this alternate activity-mask model.
- [ ] Migrate every equilibrium consumer from `CustomSubstance` and
  `ThermodynamicsCalculatorTrait` to the narrow phase interfaces and the
  `phase_equilibrium_adapter`.
- [ ] Replace legacy `Fn(..., Option<Vec<f64>>, Option<f64>) -> f64` callbacks
  with a borrowed, fallible composition/state contract at the equilibrium
  solver boundary. Only then remove the legacy infallible closure adapters.
- [ ] Move one-time matrix orientation/validation to that new numeric closure
  boundary, then remove per-call matrix work from legacy Lagrange paths.
- [ ] Complete equilibrium stories: phase-qualified element conservation,
  equation count/identity, and two occurrences of the same substance in
  different phases.

#### P2: lower-level data and performance follow-up

- [ ] Remove the `SubsData::map_of_phases` compatibility projection after its
  property calculators accept phase-level state from `PhaseSpec` directly.
- [ ] Benchmark the phase/equilibrium hot path before sharing more `Expr`
  values or compiled functions. `PhaseFunction` is already shared through
  `Arc`; do not add ownership machinery without a measured bottleneck.
- [ ] Add randomized/property permutations beyond the existing exhaustive
  three-phase ordering matrix when a property-test dependency and stable
  generator policy are selected.

#### P3: repository-wide documentation debt

- [ ] Repair or explicitly mark as `text`/`ignore` the unrelated broken
  doctests in `User_substances`, `ChemEquilibrium`, `settings`, and
  `library_manager`. The phase-layout diagram is already `text`; this task
  must not be mixed with phase-engine behavior changes.
- [ ] Triage the remaining legacy length/contains-only assertions. New phase
  engine stories already use exact order, numeric values, and typed errors;
  broaden old assertions only together with the contract they are intended to
  protect, especially for equilibrium behavior.

## Phase/Solution Linear Execution Plan

1. **Pass 1: correctness scaffold**
   Add typed ids, `SystemLayout`, and failing regression tests for repeated
   species, ordering, variable collisions, and clone equivalence. Adapt the old
   structures internally without changing the public facade yet.
2. **Pass 2: composition and equation paths**
   Migrate `indexed_moles_variables()`, `create_full_map_of_mole_numbers()`, the
   element matrix, and Lagrange builders to the shared layout. Remove the P0
   panic paths and validate all closure inputs before publication.
3. **Pass 3: canonical phase engine**
   Introduce `PhaseSpec`/`PhaseModel`/`ResolvedPhaseSystem`, make construction
   and resolution separate typed steps, and turn `OnePhase` into a thin facade.
4. **Pass 4: state and API compression**
   Replace parallel public caches with typed snapshots, split the broad trait,
   privatize mutation, and migrate equilibrium consumers through the adapter.
5. **Pass 5: optimization and final matrix**
   Borrow hot-path inputs, share closures/expressions, remove duplicated/dead
   code, and complete the deterministic story/property test matrix.

Completion means the model can represent one substance in multiple phases
without identity loss, every matrix/vector/expression uses one deterministic
layout, cloning is truthful, invalid inputs are typed failures, and "solution"
has an explicit thermodynamic model rather than being only a name.
