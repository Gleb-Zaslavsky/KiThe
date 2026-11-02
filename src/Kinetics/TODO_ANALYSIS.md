# Kinetics TODO analysis

Scope: `src/Kinetics/`, including `mechfinder_api.rs` and `mechfinder_api/`.
Exclude `experimental_kinetics/` and its module entry file
`experimental_kinetics.rs` from reading, searching, changes, and test execution
in this workstream. Do not follow callers into the excluded subsystem.

This file is an analytical TODO, not an implementation plan locked in stone. It records current risks in correctness, API ergonomics, test coverage, and solver integration.

## Current audit and implementation queue (2026-10-26)

This section supersedes the priorities and broad completion claims in the
historical notes below. The initial audit was read-only; implementation has
since begun with KIN-01 and KIN-02, including new regression tests. Findings
marked **confirmed by inspection** were established from the code; they are not
claims that every finding has a reproducer yet. Passing tests establish only
their tested contracts.

Follow [KiThe Architectural Guideline](../../KiThe_ARCHITECTURAL_GUIDELINE.md).
Correct physical outputs, input validation, and consistent publication take
precedence over API cosmetics and performance work. Never repair an erroneous
formula by updating expected values to match production, or make a test pass
by weakening its physical acceptance criteria.

### Priority map

| ID | Priority | Work package | Evidence/status |
| --- | --- | --- | --- |
| KIN-01 | P0 | Restore concentration-dependent falloff | Concentration basis documented as source-mechanism units; numeric/symbolic limits pass |
| KIN-02 | P0 | Correct four-parameter Troe and establish independent references | Independent values plus analytic-vs-finite-difference temperature derivative pass |
| KIN-03 | P0 decision | Establish pressure-rate model and units before changing interpolation | Linear implementation confirmed; PLOG provenance unresolved |
| KIN-04 | P0 | Make stoichiometry correct, repeatable, and panic-safe | Grammar edge policy and species-dependent composition invalidation implemented |
| KIN-05 | P1 | Validate rate payloads and numerical domains | Complete at `ReactionData` ingestion/evaluation boundary; unchecked low-level compatibility APIs documented |
| KIN-06 | P1 | Fix solid-state IVP preparation/result lifecycle | Confirmed by inspection; open |
| KIN-07 | P1 | Reject invalid formula/mass inputs safely | Empty formulas and zero atom counts are rejected; focused tests pass |
| KIN-08 | P1 | Make catalog loading and cache identity coherent | Multimap lookup, content provenance hash, transactions and instance close implemented |
| KIN-09 | P1 | Preserve mechanism, shortcut, and JSON input meaning | Directed closure semantics and local fixture implemented |
| KIN-10 | P1 | Complete typed error boundaries | Partial typing exists; open |
| KIN-11 | P2 | Add opt-in catalog I/O diagnostics for caller telemetry | Implemented; no performance gates or timing guarantees |
| KIN-12 | P2 | Consolidate ownership, remove proven dead helpers, and finish reporting | Representation/helper audit and facade story/example complete; campaign scripts remain last |

Recommended sequence: KIN-01/02 with their independent tests; decide KIN-03;
KIN-04; KIN-06; KIN-07/08/09. Apply KIN-10 incrementally at each touched
boundary, not as a second rewrite of all signatures. KIN-11 is observability
for catalog loading, not a performance-optimization campaign. Keep KIN-12
cleanup scoped to demonstrated complexity and finish campaign infrastructure
only after correctness contracts stabilize.

### KIN-01: concentration-dependent falloff

Owner: [FalloffStruct](mechfinder_api/kinetics.rs), numeric `K_const` and symbolic
`K_expr`; public dispatch in [ReactionData](mechfinder_api.rs).

**Confirmed by inspection:** reduced pressure is `K0 / K_inf`, without the
effective collider concentration. The computed `Eff` shadows its outer binding,
and early returns bypass the final multiplication anyway. Both frontends can
therefore ignore changes to the supplied composition. The current test checks
only a positive constant and preservation of Troe metadata.

- [x] Define concentration and rate-coefficient units at the public boundary:
  concentrations must match the source mechanism's Arrhenius convention; no
  implicit conversion occurs. CHEMKIN-style data commonly use mol/cm^3.
- [x] Implement `M_eff = sum(efficiency_i * concentration_i)`, with default
  efficiency one for unlisted species, including the no-efficiency-map case.
- [x] Use `Pr = k0 * M_eff / k_inf` in the falloff law. Do not repair the bug by
  merely multiplying the final rate by `M_eff`; that changes the high-pressure
  limit incorrectly.
- [x] Define the zero-collider limit without evaluating `log10(0)` in Troe.
- [x] Test zero composition, explicit efficiencies, default efficiency for an
  unlisted collider, concentration dependence, and low/high-pressure limits.
- [x] Compare numeric and symbolic outputs on identical inputs and
  also compare against independently calculated reference values.

**Done when:** the public facade obeys the concentration limits and agrees with
independent fixtures. Constructor metadata preservation alone is not closure.

**Progress (2026-10-26):** numeric and symbolic falloff routes now include
effective collider concentration in reduced pressure, and exact zero collider
concentration returns zero before Troe logarithms. Added Lindemann limiting and
reference checks plus Troe reference/parity checks. Parameter validation is
covered under KIN-05. Concentration values must use the source mechanism's
coefficient convention; CHEMKIN-style data commonly use mol/cm^3. Pressure
model units remain open under KIN-03.

Reference: [Cantera falloff definitions](https://cantera.org/stable/reference/kinetics/rate-constants.html#falloff-reactions).

### KIN-02: four-parameter Troe correctness

Owner: [FalloffStruct](mechfinder_api/kinetics.rs).

**Confirmed by inspection:** the numeric four-parameter branch uses a term
proportional to `A * exp(-T/T2)`, while the conventional fourth term is
`exp(-T2/T)`. The symbolic branch additionally substitutes a multiplier `T`
where the numeric branch uses `A`. The numeric/symbolic mismatch and the
independent formula error are separate defects.

- [x] Verify imported parameter order and implement the documented three- and
  four-parameter forms consistently.
- [x] Freeze independent values for `F_cent`, broadening `F`, and final `k` at
  several temperatures and reduced pressures; retain source and units.
- [x] Add explicit three-parameter, four-parameter, and absent-Troe tests.
- [x] Validate parameter lengths and the domains required by divisions and logs.
  Invalid parameters must produce typed failures, not printed messages and zero.
- [x] Check temperature derivatives of the symbolic rate at regular interior
  points against an independent finite-difference check where appropriate.

**Done when:** both routes satisfy the reference and parity checks, including
the fourth term. Do not bless agreement between two equally wrong formulas.

**Progress (2026-10-26):** numeric and symbolic `F_cent` now use
`(1-A) exp(-T/T3) + A exp(-T/T1) + exp(-T2/T)` for four parameters; the
three-parameter case omits the final term. Three temperatures and effective
collider levels are frozen for both forms, and public `ReactionData` dispatch
is covered. Symbolic parity uses relative tolerance `1e-9` because the current
RST `Expr::log10` representation divides by the rounded constant `2.30258509`;
the numeric route separately retains a tight absolute reference check.
Malformed parameter lengths/domains are covered under KIN-05. A central
finite-difference regression now checks `dK/dT` for three- and four-parameter
Troe forms at three interior temperatures, using the numeric route as an
independent finite-difference oracle.

Reference: [Cantera Troe definition](https://cantera.org/stable/reference/kinetics/rate-constants.html#the-troe-falloff-function).

### KIN-03: pressure model, interpolation, and units

Owner: [PressureStruct](mechfinder_api/kinetics.rs), source reaction catalogs,
and [ReactionData pressure entry points](mechfinder_api.rs).

**Confirmed:** interpolation is linear in pressure and rate. Sorting and
endpoint lookup have already been improved. **Unresolved:** whether every
imported `pres` record represents PLOG or an intentionally different model.
Do not label the physical model fixed merely because interpolation is ordered.

- [ ] Trace catalog provenance and document pressure, activation-energy,
  concentration, and pre-exponential-factor units. Do not infer units from
  coefficient magnitude or silently reinterpret existing payloads.
- [ ] If records are PLOG, implement interpolation in `log(P)` and `log(k)`;
  if linear interpolation is intentional, expose/document it as a distinct
  model rather than silently calling it PLOG.
- [ ] Decide how multiple Arrhenius terms at one pressure are represented and
  summed. A single `Vec<f64>` per textual pressure key does not express that
  general PLOG contract.
- [x] Reject empty pressure tables, malformed/non-finite/non-positive pressure
  keys, invalid node triplets, and ambiguous numeric duplicates such as `"1"`
  and `"1.0"` at the `ReactionData` boundary. The interpolation model itself is
  still unresolved.
- [ ] Test exact knots, interior points, the documented outside-range policy,
  single-node data, empty data, and invalid pressure inputs.
- [ ] Replace the self-referential interpolation oracle with independently
  calculated values; retain numeric/symbolic parity as a separate gate.

**Done when:** data provenance determines the model, units are explicit, and
tests protect that physical contract rather than only the current helper.

Reference: [Cantera PLOG definition](https://cantera.org/stable/reference/kinetics/rate-constants.html#pressure-dependent-arrhenius-rate-expressions-p-log).

### KIN-04: stoichiometry and reaction parsing

Owner: [StoichAnalyzer and analyse_substances](stoichiometry_analyzer.rs);
facade integration in [KinData::analyze_reactions](User_reactions.rs).

**Original audit findings (now regression-covered):** `resize` preserved old
matrix values, repeated analysis accumulated coefficients, `.position()` sent
identical reactions to one row, repeated participants overwrote entries,
unchecked side indexing could panic, and discovery could duplicate species.

- [x] Parse an equation into a checked intermediate representation once, with
  exactly one supported arrow and valid nonempty sides/tokens.
- [x] Define support for parentheses, phase/state markers, charges, explicit
  powers, and fractional coefficients. Reject unsupported syntax explicitly;
  do not remove meaningful parentheses as generic database artifacts.
- [x] Build fresh matrices and enumerate reaction rows; never infer row identity
  from equation-string equality. Matrices are built off to the side and only
  published after every reaction succeeds.
- [x] Aggregate repeated species consistently. Repeated terms sum both
  stoichiometric coefficients and their per-term kinetic orders, including
  explicit `^` overrides; this is now the documented parser behavior.
- [x] Deduplicate species during discovery while preserving first-seen order.
- [x] Reject a supplied species universe that omits a participant rather than
  silently dropping its coefficient.
- [x] Commit reaction and matrix analysis only after all equations pass.
- [x] Invalidate or rebuild dependent elemental and mass data when the species
  universe changes, respecting explicit overrides.
- [x] Add tests for repeated analysis, duplicate equations, repeated
  participants, missing/multiple arrows, malformed tokens, unknown species, and
-  a failure after a valid first reaction, matrix dimension changes, and a
  balanced elemental reference reaction.
- [x] For a balanced reference reaction, assert elemental conservation and
  agreement of stoichiometric and kinetic-order representations where expected.

**Done when:** direct analyzer and facade routes produce consistent matrices,
are idempotent for unchanged inputs, and fail without partial publication.

**Progress (2026-10-26):** the analyzer now parses each equation into a checked
internal representation shared by discovery and matrix generation; it rejects
absent/multiple arrows, empty sides/participants, malformed powers, and
participants missing from an explicit species list. Each reaction is assigned
its own row, all matrices are rebuilt from zero and published transactionally,
repeated terms sum coefficients/orders, and discovery preserves first-seen
unique species. Balanced formula parentheses are preserved; third-body markers
are excluded; recognized phase tags remain in kinetic species identity but
are stripped for elemental/mass projection. Terminal ionic `+`/`-` and caret
magnitudes such as `Fe^3+` are supported. Stoichiometric coefficients accept
positive finite decimal/scientific notation; rational text such as `1/2` is
rejected. Element/mass data are invalidated when species or group definitions
change; `set_substances` invalidates eagerly, and legacy direct field mutation
is detected at the next analyzer operation. Explicit molar-mass overrides are
retained only for an unchanged species/group snapshot and matching vector
length. Verification: 16 analyzer tests pass in debug mode; no performance
claim is made.

### KIN-05: numerical input and result validation

Owners: [reaction payloads](mechfinder_api/kinetics.rs),
[dispatch/deserialization](mechfinder_api.rs), and IVP model constructors.

- [x] Validate payloads immediately after deserialization in both mechanism
  ingestion APIs. `ReactionData::validate()` is also public for standalone
  `serde` users.
- [x] Add `try_new_*` constructors that validate each reaction payload
  immediately.
- [x] Preserve the serialized/public vector payload shape for compatibility;
  provide validating `try_new_*` constructors and validate before every
  `ReactionData` evaluation. Raw `new_*` constructors and low-level kernel
  structs remain unchecked compatibility APIs; a fixed-shape representation
  is deferred to a major API revision rather than silently breaking callers.
- [x] Validate finite payload parameters, positive finite temperatures and
  pressure inputs, and finite non-negative numeric concentrations/efficiencies
  at the `ReactionData` evaluation boundary. Symbolic concentration trees are
  recursively checked for non-finite constants.
- [x] Define Arrhenius-domain rules: `A >= 0`; temperature exponent and
  activation energy may be any finite values, including negatives. Pressure
  interpolation and units remain a separate KIN-03 contract.
- [x] Validate temperature sweep boundaries/counts before allocation. Boundaries
  must be finite/positive, and requests are limited to 100,000 points with a
  typed `ResourceLimitExceeded`; consider streaming only if larger sweeps become
  a real use case.
- [x] Reject non-finite computed rates with a dedicated error. Document finite
  `0.0` as valid, including exact zero parameters and IEEE-754 underflow.
- [x] Ensure `ReactionType` and payload agree at public evaluation boundaries;
  `Empirical` currently returns a typed not-implemented error at evaluation.
  This is a planned symbolic arbitrary-form rate model, not invalid input; no
  production task currently exercises it.
- [x] Recursively reject non-finite constants anywhere in symbolic concentration
  trees; regression covers nested NaN and a valid composed expression. Existing
  coverage includes malformed/empty pressure tables and duplicate
  pressure keys, direct non-finite symbolic constants, short/non-finite arrays,
  invalid T/P/concentrations/efficiencies, malformed Troe, overflow, and direct
  construction that bypasses constructors.

**Done when:** at supported `ReactionData` ingestion/evaluation boundaries,
malformed input yields a contextual typed error rather than panic, `Ok(NaN)`,
or an invented zero coefficient.

**Progress (2026-10-26):** `ReactionData::K_const` and `K_sym` now validate the
payload/tag pair before dispatch. Numeric evaluation checks Arrhenius triplets,
non-negative pre-exponential factors while allowing negative exponents and
activation energies, Troe shape/parameter domain and characteristic
temperatures, efficiency tables, positive finite temperature/pressure, and
finite non-negative concentrations. Rate overflow returns `NonFiniteRate`;
finite zero, including floating-point underflow, is documented and tested.
`Empirical` is recognized as a planned symbolic arbitrary-form rate but has no
implementation yet. Temperature sweep boundaries are checked before building
the grid, and both mechanism parser entry points validate payloads on ingestion.
Empty/invalid/ambiguous pressure nodes are rejected before interpolation.
Validating `try_new_*` constructors coexist with explicitly unchecked legacy
constructors. Materialized sweeps reject more than 100,000 samples before
allocation via `ResourceLimitExceeded`. Verification: 16
`mechfinder_api::kinetics` tests, 14 `mechfinder_api` tests, and 107
`User_reactions` tests passed in debug mode. The pressure interpolation law and
units remain KIN-03 work; arbitrary symbolic operations are validated for
finite constants, not for runtime domains at every possible variable value.

**Status:** KIN-05 is complete for the supported `ReactionData` ingestion and
evaluation boundary. The legacy low-level structs and infallible constructors
remain unchecked for compatibility; callers needing input validation should
use `try_new_*` or call `ReactionData::validate()` before evaluation. Moving
those raw APIs to fixed-shape validated types is a major-version follow-up.

### KIN-06: solid-state IVP lifecycle and publication

Owner: [KineticModelIVP](solid_state_kinetics_IVP.rs).

**Confirmed by inspection:** `set_model` embeds current thermal/Arrhenius values
in the expression. Later setters do not rebuild it or invalidate old solver
results. `get_solution` uses current fields for temperatures and conversion
rates, potentially combining a previous solution with a new problem.

- [ ] Separate problem configuration, prepared equations, and accepted result.
  Choose lazy rebuild at solve or consistent setter invalidation, with one owner.
- [ ] Invalidate results after a successful relevant mutation. Invalid setter
  input should not partially modify the previous valid configuration.
- [ ] Make setter ordering irrelevant to the final configured problem.
- [ ] Publish result metadata from the solved problem snapshot, not mutable
  current fields. Define behavior after a failed solve following a successful one.
- [ ] Validate finiteness as well as sign; audit model parameter domains and
  bounds near singular conversion endpoints.
- [ ] Audit solver initialization and failure propagation through the current
  RST API; preserve typed causes, status, and available diagnostics.
- [ ] Decide the purpose of `jac`: it is stored/checked but not passed to the
  universal solver here. A zero derivative is not automatically an invalid model.
- [ ] Test setter permutations, parameter changes after preparation/solve,
  successful solve followed by failure, and output shape/finite/domain checks.
- [ ] Check analytic derivatives numerically, rather than testing only that a
  symbolic string is nontrivial. Keep backend parity and independent model
  references separate.

**Done when:** configuration, integrated equations, and every published output
describe the same problem; stale results cannot masquerade as a new solution.

### KIN-07: formula parsing and molecular mass

Owner: [molmass](molmass.rs), including calls from stoichiometric analysis.

**Progress (2026-10-26):** physical mass/composition now reject unsupported
symbols with typed errors; symbolic species parsing remains separate. Formula
syntax is ASCII-validated before byte indexing, bounded by byte length and
nesting, and count arithmetic is checked. Custom group definitions are
validated and nested/cyclic expansion is explicitly unsupported. Only a
recognized terminal phase suffix is stripped. Element-matrix columns have a
deterministic order and the species-by-element orientation is tested. The
public composition-mass helper now returns `KineticsResult<f64>`; its sole
in-workspace production caller preserves the typed source.

- [x] Reject unknown elements when computing physical mass/composition. Preserve
  explicit custom groups and distinguish symbolic species names from formulas.
- [x] Decide empty-formula and zero-count semantics; do not silently publish an
  empty physical composition or zero mass for invalid input.
- [x] Use a consistent byte/character parsing strategy. If only ASCII formula
  syntax is supported, reject other input safely before slicing.
- [x] Use checked atom-count arithmetic and bound nesting/input size to avoid
  overflow, pathological recursion, or excessive parser work.
- [x] Validate custom group expansions, including unknown atoms and any supported
  nesting; reject unsupported cyclic/ambiguous definitions explicitly.
- [x] Restrict phase-marker removal to recognized positions so it cannot erase
  genuine parenthesized chemical groups.
- [x] Give elemental matrix columns a deterministic order and return that order
  with the matrix. Preserve documented orientation and downstream contracts.
- [x] Cover unknown symbols, UTF-8, huge counts, deep/unbalanced parentheses,
  custom groups, phase suffixes, and stable repeated matrix construction.

**Done when:** invalid formulas cannot silently lose mass and all public parser
failures are deterministic typed errors rather than panic or wrapped arithmetic.
Empty formulas, empty physical compositions, and zero atom counts are rejected
as invalid physical input. The regression covers formula syntax, composition
mass, and zero-count custom-group definitions.

Focused verification: 14 `Kinetics::molmass::tests`, 12
`Kinetics::stoichiometry_analyzer::tests`, the relevant `User_substances`
composition/molar-mass integration test, and the calculation-source error test
pass in debug mode. The full Kinetics test suite was not used because this
workstream explicitly excludes the separate IVP module.

### KIN-08: catalog transactions and cache identity

Owner: [KineticData and shared caches](kinetics_lib_api.rs); configuration
boundary in the repository's library manager.

**Progress (2026-10-26):** catalog pairs are resolved from one configuration
snapshot, each process cache is keyed by configured path plus file length and
modification time, and a file changing during load is retried/rejected.
`open_json_files` resolves both library halves before mutating the receiver;
successful switching clears dependent equation/search state, while lookup
failure preserves the previous snapshot. Initialization is concurrency-tested.
The cache revision remains metadata-based (configured path, size, and mtime) to
avoid rereading large catalogs on every cache hit. Cold loads compute and expose
SHA-256 of the exact loaded bytes as content provenance; this is not used as a
hot-path invalidation key, so same-size/same-mtime external rewrites are not
detected until metadata changes or the cache is otherwise refreshed.

- [x] Resolve both catalog halves before committing a library switch.
- [x] On successful switching, rebuild or invalidate equation indexes and search
  results; on failure preserve the previous coherent library snapshot.
- [x] Establish cache identity and invalidation for configured paths/data
  revisions. Do not mix parameter and relationship catalogs from different
  configurations; record their provenance together.
- [x] Test first-load failure/retry, A-to-B switching, failed switching,
  configuration changes, and concurrent initialization.
- [x] Decide how repeated equations map to multiple reaction IDs; a single
  `equation -> id` mapping must not silently discard valid duplicate reactions.
- [x] Preserve deterministic ordering in public query results where reproducible
  identity/order is part of the contract.

**Done when:** catalog identity, payloads, indexes, and selected results remain
coherent across success, failure, configuration updates, and reuse. Duplicate
equation lookup preserves every matching ID; plural lookup returns all matches,
while singular lookup returns a typed ambiguity error.

Focused verification: all 20 `Kinetics::kinetics_lib_api::tests` pass, including
cache revision/retry, concurrent initialization, failed-switch snapshot
preservation, and library switching. Cache freshness uses path, file length,
and modification time; the exposed SHA-256 identifies loaded bytes but does
not invalidate a warm entry when metadata is unchanged.
Lifecycle contract: `KineticData::close` releases the instance-owned library
and derived query state. It intentionally does not evict the process-wide parsed
catalog cache, which may be reused by other callers or later opens. The cache
retains its last parsed revision per catalog for process lifetime.

### KIN-09: mechanism and input fidelity

- [ ] Audit and repair or quarantine malformed Troe payload(s) in the default
  GUI mechanism. Bulk parsing reports a contextual reaction index, but that
  index is the iteration position of an unordered map and is not a stable
  source-record identifier. The payload is correctly rejected by
  `ReactionData::validate` for invalid Troe characteristic temperatures.
  The GUI bulk-add path must remain all-or-nothing; restore the success story
  after the source catalog record is corrected or explicitly excluded by the
  catalog's supported-mechanism policy.

Owners: [mechfinder](mechfinder_api/mechfinder.rs),
[public parsing](mechfinder_api.rs), [shortcuts](parsetask.rs), and `KinData`.

- [x] Check ID/payload lengths before `parse_kinetic_data` uses `zip`; reject
  mismatches rather than silently truncating the longer input.
- [x] Define duplicate-ID handling instead of silently overwriting an earlier
  reaction payload.
- [x] Extract equation strings as strings, not JSON-encoded strings.
  `mechfinder` currently uses `serde_json::to_string` for `eq`; trace and test
  `construct_mechanism -> kinetic_main` to catch quoted participant names or
  disagreements with `ReactionData.eq`.
- [x] Define mechanism reachability semantics, including reversible reactions,
  third-body markers, missing relationship fields, and empty initial species.
  Audit whether forward-only expansion matches the intended public contract.
- [x] Test graph closure and determinism using a small local fixture before
  relying on a large database smoke test.
- [x] Validate shortcut grammar and preserve meaningful reaction IDs. Audit
  library aliases against actual catalog keys: lowercasing makes uppercase
  `B` match arms unreachable, and `ARAMCO` versus `Aramco` needs reconciliation.
- [x] Replace `Unknown` sentinel strings with explicit parse/lookup outcomes and
  remove unconditional parsing debug prints.

**Done when:** selection, parsing, analysis, and export preserve reaction
identity and equation content end to end without silent omission or rewriting.
Shortcut parsing now accepts only an ASCII library prefix plus optional
separator and an integer ID, returns typed errors, and maps aliases to catalog
keys (`Aramco`, not `ARAMCO`). Mechanism parsing rejects repeated IDs before
insertion. Mechanism IDs are sorted numerically as strings, preserving leading
zeros and values larger than machine integer types; nonnumeric IDs produce a
contextual typed error. The public grouped-map convenience view remains `Option`-based and
therefore reports invalid canonical shortcut metadata as unavailable; the
resolution operation itself preserves the typed parse error. Reachability is
directed: a row is selected when every non-inert reagent is known, then its
products are added. Reversible-arrow interpretation is not inferred from
equation text; the relationship catalog is authoritative. Empty seed sets
select nothing. Third-body marker tokens are ignored on either side, and
missing `reagents` or `products` fields produce a typed payload error with
reaction ID. A local graph fixture covers multi-step closure, multi-reagent
gating, unreachable rows, markers, empty seeds, and malformed rows.
Focused verification: 3 shortcut parser tests, 6 mechanism ID/equation/graph
helper tests, 17 mechanism API tests, and all 107 `User_reactions` tests pass in debug
mode. Deterministic local assertions now prove equation extraction returns the
string content rather than JSON-encoded text; the existing live mechanism
smoke test covers integration, while local graph assertions establish the
catalog-independent closure contract.

### KIN-10: typed errors and idiomatic boundaries

Owner: [KineticsError](error.rs) and each affected public API.

- [ ] Retain `KineticsResult` as the common boundary and introduce contextual
  variants for input/domain, equation/formula parse, payload shape, lookup,
  unsupported capability, lifecycle, numerical failure, and output rejection.
- [x] Preserve serde causes from both reaction parser entry points and attach
  reaction index, ID, and mechanism context where available. Add source-chain
  regression coverage.
- [x] Give stoichiometric reaction parsing a distinct syntax-error category,
  preserve regex-construction sources, and report unknown reaction species as
  a structured lookup error containing both the species and equation.
- [x] Distinguish malformed serialized rate payloads from invalid runtime
  values and missing runtime inputs; include reaction identity and field name.
- [x] Return structured operation/current-state fields when a `KinData`
  workflow operation is not allowed by its lifecycle state.
- [x] Report invalid lifecycle payloads as a state plus a typed unmet
  `KinDataStateRequirement`, not a formatted generic `InvalidState` string.
- [x] For the in-scope non-IVP APIs, include property/field, catalog/path, and
  stage context wherever it is useful and stable.
- [ ] Audit public constructors and mutable fields that bypass validation.
  Make validated ownership enforce invariants without requiring repeated full
  payload validation on every hot callback.
- [ ] Audit remaining in-scope modules for blanket warning suppression and
  replace it only where the underlying code can be safely cleaned or justified.
- [x] Remove blanket warning suppression from `mechfinder_api` and its rate
  model module; retain only a narrow allowance for the historical exported type
  spelling. Remove the unreachable falloff tails.
- [x] Test parser error categories, context, and serde source chains separately
  from formatted error wording. IVP error-category coverage remains open.

**Progress (2026-10-26):** preparatory cleanup removed blanket warning
suppression from the rate-model modules, deleted unused private conversion and
deserialization helpers, removed the unreachable falloff result tails and
debug printing, and replaced ignored test `Result`s with explicit expectations.
The legacy `Mechanism_search` name is retained with a narrow compatibility
allowance. Focused rate-model and mechanism-parser tests pass. This cleanup does
not yet implement the structured error taxonomy or migrate IVP string errors.

The first typed boundary is now in place: parser errors are wrapped in a
source-preserving `Context` variant, and parser tests assert both the typed
serde cause and record/mechanism context. `parse_kinetic_data` also rejects
ID/payload length mismatches using `LengthMismatch` instead of truncating via
`zip`. Stoichiometry now returns `ReactionSyntax`, `RegexCompile`, and
`UnknownReactionSpecies` rather than encoding these cases as generic payload
strings. Tests assert the structured fields, and the facade suite passes through
the new errors unchanged. Broader structured domain/lifecycle errors remain.

The rate-evaluation boundary now has separate `InvalidReactionPayload`,
`InvalidKineticInput`, and `MissingKineticInput` categories. Payload failures
carry reaction and field; runtime-domain failures carry reaction and input
field; absent pressure/concentration arguments are no longer confused with
malformed reaction data. Regression tests assert these distinctions through
numeric dispatch and both mechanism parser entry points. Current focused
verification: `mechfinder_api::kinetics::tests` (16 passed) and
`mechfinder_api::tests::` (16 passed). Generic errors remain in the standalone
temperature-sweep utility and reaction-tag/payload mismatch path for a later
focused taxonomy pass. Mechanism lookup now distinguishes a missing parameter
catalog from a missing relationship catalog; catalog acquisition and graph
closure failures carry operation context while retaining typed causes.

Unchecked `ReactionData::new_*` builders remain public because GUI and reactor
clients use them, including the separately scoped IVP subsystem; this audit
does not migrate those callers. The validated `try_new_*` methods are the
recommended construction path. Public fields and serde deserialization also
permit invalid states, so rate evaluation continues to validate at its
boundary. Fully validated ownership and avoiding repeated callback validation
remain open until the public-field contract is migrated to accessors or a
validated wrapper.

The kinetic-library cache loaders now attach catalog and configured path
context at file open, read, and JSON parse stages while preserving the typed
`Io` or `Json` source. Failure-path regression tests verify both the context
and source category; the complete `kinetics_lib_api::tests` suite passes (12).
Poisoned Reactbase and reaction-database locks now return the structured
`CacheLockPoisoned` error with the cache identity. Remaining broad taxonomy
work includes internal state-contract mismatches and legacy validation strings.
Public lifecycle guards now return `OperationNotAllowedInState` with structured
operation and current-state fields; state-machine regressions assert those
fields directly. `validate_state_contract` now reports
`StateContractMismatch { state, requirement }`, with a typed requirement enum
covering shortcut, library, direct, parsed, mechanism, analyzed, and normalized
state prerequisites. The existing state-contract tests assert variants and
fields instead of diagnostic prose.

**Deferred separate workstream:** `solid_state_kinetics_IVP.rs` and its solver
API are intentionally excluded from this KIN-10 pass. They form an independent
subsystem and need a coordinated lifecycle/API refactor; do not migrate their
`Result<_, String>` signatures piecemeal here. Track that work with KIN-06 and
return to its typed-error policy as part of that dedicated refactor.

**Done when:** callers can decide how to handle a failure from structured data,
and a public `Result` method does not hide reachable indexing panics below it.

### KIN-11: opt-in catalog I/O diagnostics

Kinetics is a domain library commonly called from higher-level solvers that
already own workflow telemetry. Do not build a second solver-wide telemetry
system here. Formula construction and rate-expression preparation are not
current performance targets; the material observability gap is loading and
parsing the JSON catalogs.

- [x] Expose opt-in diagnostics for catalog configuration snapshot, cache hit or
  miss/revision invalidation, file open/read, JSON deserialization, and selected
  library materialization. Include catalog identity/path, bytes where available,
  and elapsed duration for each measured operation.
- [x] Make diagnostics easy for a caller's existing telemetry collector to
  receive and aggregate. Keep the default path quiet and avoid global logger
  setup.
- [x] Distinguish a cache hit from a disk read/parse and report failed stages;
  an unavailable duration/counter is not zero.
- [x] Document the load-use-release lifecycle: open/parse a catalog once, make
  repeated queries against that loaded view, then release the view. Clarify that
  the current shared process cache survives dropping a `KineticData` instance;
  decide whether explicit cache eviction is actually needed.
- [x] Do not add performance gates, benchmark pass/fail thresholds, or speedup
  claims for OS-sensitive file I/O and JSON parsing. Timings are diagnostic
  observations for the calling application's telemetry, not a performance
  contract of this crate.

**Done when:** higher-level callers can attribute catalog load and parse work
through their own telemetry, reuse a loaded library for a batch of queries, and
understand exactly which memory is released at the end of that batch.

**Progress (2026-10-26):** `open_json_files_with_diagnostics` returns opt-in
timings for configuration snapshot, cache revision/lookup, file open/read, JSON
deserialization, library selection, and instance materialization. Failures
carry a partial report and the original typed error. The ordinary
`open_json_files` path does not call clocks or construct a report. `close()`
releases the selected instance view and query/index state while preserving the
shared process cache. Regression tests cover cold parse versus cache hit,
parse-failure attribution, close/reopen reuse, and selection-failure reports.
Durations are diagnostic only; no timing thresholds or performance gates were
added.

### KIN-12: maintainability and final evidence campaign

- [x] Audit current reaction representation ownership. `ReactionIdentity` is
  canonical identity; raw JSON and parsed `ReactionData` represent distinct
  input/lifecycle stages; equations, stoichiometry, and symbolic constants are
  derived views; `EveryReaction` is the materialized normalized/export view.
  These are not all interchangeable duplicates. Preserve current transactional
  synchronization and invariant tests; defer removing compatibility fields or
  collapsing storage until a separate migration has a concrete invariant/API
  plan. Do not do a broad naming-only rewrite.
- [x] Audit apparently unused helpers before deprecation/removal:
  `convert_strings_to_f64` and `after_bracket_stoichio` are absent from the
  current source tree. Removed `clean_off_artifacts`, an obsolete generic
  string cleaner used only by its own unit test; production parsing now has
  targeted third-body handling and preserves meaningful formula parentheses.
  Limited the JSON loader and cache accessor wrappers plus sample JSON payloads
  to test builds, since their only call sites are tests.
  `F64Wrapper` remains public but has no in-repo callers; retain it for now
  because deleting a public type may break external users. This audit excluded
  `experimental_kinetics` and IVP by design.
- [ ] Separate parsing/evaluation from opt-in logs, tables, plotting, and export.
  Keep graphical and file-writing demonstrations outside ordinary unit gates.
- [x] Add compact facade examples for direct/library/mechanism workflows and
  numeric/symbolic constants. IVP examples remain explicitly out of scope.
- [x] Add curated [STORY_TESTS.md](STORY_TESTS.md) entries with
  name, exact command/profile, hypothesis, fixture/units, acceptance criteria,
  measured output, conclusion, and limits. Do not invent historical measurements.
- [ ] As the final stage, add a non-fail-fast PowerShell campaign in `scripts/`
  using the existing repository reporting conventions. Use timestamped run
  directories with `reports/` for structured tables and `technical/` for raw
  output; record every step/status and continue after individual failures.
- [ ] Keep campaign artifacts, curated stories, and Criterion baselines separate.
  Prefer `tabled` summaries; use RST reporting utilities only if their API fits
  the local campaign without unnecessary coupling.

**KIN-12 progress (2026-10-26):** added a database-independent direct-facade
example that demonstrates stoichiometric analysis and numeric/symbolic rate
constants, plus matching correctness stories for facade processing and located
batch failures with shortcut/library identity, batch position, and equation
context. Existing guides already cover
local lookup, shortcut selection, and mechanism construction. The helper audit
also removed an obsolete test-only cleaner and constrained test fixtures/cache
wrappers to test builds. The in-repo reaction representations were reviewed;
they encode separate lifecycle/derived views, so no broad storage rewrite was
justified in this cleanup. The remaining KIN-12 task is the repository-wide
campaign script/reporting system; it stays last as planned.

### Observed test baseline and missing gates

The initial audit ran these existing debug tests with `--lib --no-default-features`,
one Cargo build job, and one test thread. All 79 selected tests passed:

| Filter under `Kinetics::` | Passed | Scope |
| --- | ---: | --- |
| `mechfinder_api::kinetics::tests` | 10 | Existing rate/parity and validation checks |
| `mechfinder_api::tests::` | 12 | Existing parsing and mechanism API checks |
| `stoichiometry_analyzer::tests` | 6 | Existing parser/matrix examples |
| `molmass::tests` | 8 | Existing formula/mass examples |
| `solid_state_kinetics_IVP::tests::` | 9 | Model construction/expression checks |
| `solid_state_kinetics_IVP::kinetic_model_ivp_tests::` | 23 | Existing IVP checks, excluding `solve_with_` plot demos |
| `kinetics_lib_api::tests::` | 11 | Existing catalog/cache/query checks |

This was not a full Kinetics test campaign. `User_reactions_tests` and the
remaining IVP test group were not rerun in this audit. No tests from
`experimental_kinetics` were selected. The falloff positivity check and the
symbolic-shape-only Jacobian check do not prove physical formula correctness.

Commands reproducing the selected baseline:

```powershell
cargo test -j 1 --lib --no-default-features Kinetics::mechfinder_api::kinetics::tests -- --test-threads=1
cargo test -j 1 --lib --no-default-features Kinetics::mechfinder_api::tests:: -- --test-threads=1
cargo test -j 1 --lib --no-default-features Kinetics::stoichiometry_analyzer::tests -- --test-threads=1
cargo test -j 1 --lib --no-default-features Kinetics::molmass::tests -- --test-threads=1
cargo test -j 1 --lib --no-default-features Kinetics::solid_state_kinetics_IVP::tests:: -- --test-threads=1
cargo test -j 1 --lib --no-default-features Kinetics::solid_state_kinetics_IVP::kinetic_model_ivp_tests:: -- --test-threads=1 --skip solve_with_
cargo test -j 1 --lib --no-default-features Kinetics::kinetics_lib_api::tests:: -- --test-threads=1
```

For remaining findings, add focused failing regressions before implementation,
then run both low-level and facade tests. Expand to the applicable
non-experimental corpus once side-effecting demos are isolated. Use explicit
filters rather than an unqualified `Kinetics::` run in this scoped workstream.
KIN-11 intentionally has no release benchmark or performance gate: its timings
are operational diagnostics for the caller, not speedup evidence or a crate
performance contract. KIN-01/02 evidence is debug correctness coverage.

**KIN-01/02 verification (2026-10-26):** after correcting falloff and Troe, the
following focused debug tests passed (34 tests total):

| Test filter | Passed | Coverage added or rechecked |
| --- | ---: | --- |
| `Kinetics::mechfinder_api::kinetics::tests` | 17 | Effective collider and concentration limits; 3/4 parameter Troe references; symbolic parity; central finite-difference temperature derivative at three interior temperatures for both Troe forms; `ReactionData` dispatch |
| `Kinetics::mechfinder_api::tests::` | 17 | Public mechanism parsing/deserialization and related facade regression |

The current `Expr::log10` implementation uses `ln(x) / 2.30258509`; observed
symbolic/numeric parity is therefore asserted at relative `1e-9`, while the
numeric implementation is checked against frozen equation evaluations at
absolute `1e-3`. The derivative check uses a central difference of the numeric
rate path, independent of symbolic differentiation. These fixtures establish
formula regression, not external experimental rate-data validation. Malformed
rate parameters/domains are covered separately under KIN-05, not by these 34
KIN-01/02 tests.

## Historical notes: correctness and panic-safety

The sections below retain earlier findings and progress. Descriptions beneath
completed items may describe the old defect, not current code. Read their
guarantees narrowly; the current KIN queue above governs new work.

### Done so far

- [x] Introduce typed errors in selected primary runtime flows. This is not complete panic-safety; see KIN-04/05/07/10.
- [x] Fix direct-reaction workflow so direct equations can go straight through analysis.
- [x] Fix document generation with library maps.
- [x] Fix falloff metadata loss.
- [x] Fix symbolic and numeric mismatch for three-body reactions.
- [x] Fix pressure-node ordering and endpoint lookup. Physical interpolation law remains under KIN-03.
- [x] Revisit stoichiometric power parsing.
- [x] Split library loading from reaction lookup.
- [x] Tighten mutation invalidation for `KinData`.
- [x] Make `KinData` parallel/cache fields crate-private and expose read-only facade accessors.

- Replace user-facing `panic!`, `assert!`, `unwrap`, and `expect` paths with typed errors.
  Main hot spots: `User_reactions.rs`, `mechfinder_api.rs`, `kinetics_lib_api.rs`, `stoichiometry_analyzer.rs`, `molmass.rs`, `solid_state_kinetics_IVP.rs`.
  Public API should prefer `Result<_, KineticsError>` for invalid shortcuts, missing libraries/reactions, malformed equations, missing pressure/concentration data, failed file IO, JSON parse errors, and solver failures.
  Historical progress: selected facade and lookup paths were converted to typed errors. The 2026-10-26 audit found reachable unchecked indexing, invalid-domain outputs, and parser arithmetic/slicing hazards; absence of explicit `unwrap` is not proof of panic-safety. See KIN-04/05/07/10.

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

## Historical notes: LSODE2 support in solid-state IVP

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

## Historical notes: API and UX

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

## Historical notes: tests missing or too weak

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

- [x] Fix the `test_pres_deserialization` discovery gap.
  Rechecked on 2026-10-26: it has `#[test]` and passed in the selected baseline.

## Historical notes: performance and maintainability

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

## Historical verification and design questions

Use the scoped commands in the current baseline above. Earlier broad filters
such as `lsode2` do not enforce the excluded-subsystem boundary and are not the
default campaign commands for this audit.

- Resolved: `KineticModelIVP::default()` currently selects LSODE2.
- Resolved for the existing default: kinetic powers follow stoichiometric coefficients unless explicitly overridden. Repeated-participant aggregation and override grammar remain in KIN-04.
- How much backward compatibility is required for current public mutable fields and printing behavior?
- Should `KinData` preserve legacy parallel vectors temporarily, or should the next refactor move directly to a normalized per-reaction record?
