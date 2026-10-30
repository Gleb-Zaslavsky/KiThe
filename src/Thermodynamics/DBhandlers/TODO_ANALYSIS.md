# DBhandlers Technical Plan

This checklist covers correctness, error architecture, and maintainability of
the thermodynamic and transport data handlers. It complements the repository
guideline: a numerically converged fit is publishable only when it reproduces
the requested physical properties in the declared units.

## P0: NIST Shomate stitching correctness

- [x] Define one canonical Shomate design matrix shared by fitting,
  reconstruction, and validation.
  - `Cp` uses `A..E` in `J/(mol K)`.
  - `H-H298` uses `A..E` and the identifiable offset `F-H` in `kJ/mol`.
  - `S` uses `A..E,G` in `J/(mol K)`.
  - Do not solve independently for `F` and `H`: property samples identify only
    their difference. Preserve an explicit source `H` gauge and reconstruct
    `F` from the fitted offset.
- [x] Establish a deterministic linear SVD/QR route as the reference and
  preferred production implementation. The Shomate model is linear in its
  coefficients; symbolic nonlinear fitting is not required for the normal
  path.
- [x] Add offline exact-polynomial and two-range fixtures. Check `Cp`, `H`, and
  `S` independently over the complete range, including the join and endpoints.
  - Exact single-polynomial round trip and non-identical adjacent two-range
    fixture are active offline tests; the latter checks a 1 K grid with
    absolute-plus-relative tolerances and both source branches at the join.
- [x] Replace the one-size-fits-all relative error near zero with explicit
  property-specific absolute-plus-relative tolerances. In particular,
  `H-H298` legitimately crosses zero and must not be divided by an arbitrary
  tiny floor as its sole acceptance metric.
  - Production NIST SVD uses `abs_error <= 0.05 + 0.02 * |reference|` for Cp
    and S in their native units, and H in kJ/mol. The existing comparison API
    used by NASA and symbolic compatibility paths retains its prior semantics.
    Tests cover both accepted near-zero error and rejection of a larger error.
- [x] Reject malformed ranges, empty selected samples, non-finite inputs,
  rank-deficient systems, and non-finite coefficients with typed errors.
  - Production NIST SVD exposes categorized input, preparation, numerical, and
    quality errors; the facade preserves the typed source. Coefficients are
    committed only after all reports pass, and that commit is an infallible
    state replacement rather than a separate fallible operation.
- [x] Validate all three property reports before publishing coefficients.
  Failed fitting or validation must leave previously published state intact.
- [x] Replace tests that merely assert different algorithms produce different
  coefficients with output-space acceptance checks. Coefficient vectors are
  not unique under the `F/H` gauge.

## P0: Units and sampling contract

- [x] Centralize the NIST units contract and remove accidental conversions in
  fitting reports. Raw and fitted `H` must both be compared in `kJ/mol`; public
  calculator output may convert to `J/mol` exactly once at publication.
- [x] Replace `f64 -> usize` range truncation with a deterministic finite-grid
  builder that supports non-integral bounds, includes intended endpoints, and
  never wraps negative values. Reject non-positive sample temperatures and
  excessive grids before allocation.
- [x] Sort intervals and samples deterministically, define overlap and boundary
  ownership, and reject contradictory overlapping records.
  - Active regressions cover shuffled adjacent intervals, fractional requested
    endpoints, single ownership of the shared boundary, and overlap rejection.

## P1: Typed errors and transactional publication

- [ ] Split catch-all fitting strings into typed input, preparation, backend,
  convergence, quality-rejection, and publication errors. Preserve upstream
  RST errors as sources when compatibility fitting is invoked.
- [x] Make fallback evidence observable: primary failure, fallback identity,
  and final property acceptance must remain available to callers or optional
  diagnostics.
  - NASA's report API identifies the selected adjacent/non-adjacent symbolic
    or linear-SVD route, retains primary failure text when fallback succeeds,
    and returns per-property quality reports. If both attempts fail, the
    aggregate typed error preserves both causes through the thermo facade and
    user-facing error logging.
- [ ] Audit `NISTError`, `NASAError`, `CEAError`, `TransportError`, and their
  facade conversions. Keep domain-specific context instead of flattening it to
  `CalculationError(String)` or `FittingError(String)`.
  - NIST schema errors now remain distinct at the facade, and parser setup
    preserves the concrete parser error as an error source. NASA's two-route
    failure preserves both typed attempts; its symbolic RST failures now carry
    stage and source. CEA and transport conversions preserve the RST fitting
    error and stage. Remaining work is to audit CEA parsing/symbolic categories
    and TransportData's internal catch-all calculations, then add end-to-end
    source-chain tests for those paths.
- [ ] Ensure every mutating calculation builds temporary state and commits it
  only after full validation.
  - NIST interval routing now resolves both endpoints without mutating state,
    then commits the selected interval only after the fit/selection succeeds;
    a failing adjacent fit preserves both coefficients and current range.
  - CEA interval routing now stages on a clone, and adjacent-fit coefficients
    are committed only after R² acceptance. A failing interval fit preserves
    coefficients, flags, and requested range.

  - CEA deserialization now propagates coefficient-parse failures and commits
    the new input/coefficient map only after successful parsing.
  - Transport unit setters validate before publication. Facade `M`/`P` setters
    validate the numeric input and unit before replacing the stored value.
  - NASA adjacent and non-adjacent fitting now computes property reports before
    publishing the candidate coefficient tuple.

## P1: Remove obsolete symbolic fitting from production paths

- [x] Compare the corrected linear reference with RST LM on identical samples,
  weights, units, and acceptance criteria.
  - Offline parity fixture feeds the same sampled Cp/H/S rows and inverse-RMS
    weights to the 7-coordinate SVD solve and RST's three-equation LM solve.
    Both use the same physical property acceptance; weighted objectives and
    componentwise published Cp/H/S predictions are compared.
- [ ] Extend the identical-objective parity fixture to RST Gavin after the
  repository's fitting adapter exposes that backend with the same report data.
- [x] If the linear route covers the supported model, remove symbolic NIST
  stitching from the production fallback cascade. Keep only a
  focused diagnostic parity test if it still provides value.
  - Production adjacent and non-adjacent stitching are now SVD-only. The old
    public symbolic methods remain compatibility-only pending API deprecation.
- [ ] Apply the same analysis to NASA and CEA: use linear algebra whenever the
  model is linear in fitted coefficients; retain nonlinear fitting only where
  bounds or genuinely nonlinear parameters require it.

## P2: Handler idiomaticity

- [x] Make cloning total and semantics-preserving. `NISTdata::clone()` must not
  unwrap absent source arrays or silently rebuild closures from a different
  coefficient source than the published fitted state.
  - NIST closures are shared as `Arc<dyn Fn>` and cloning preserves the exact
    callable state; an uninitialized instance clones without panic.
- [ ] Replace public mutation-heavy workflows with prepared immutable data plus
  explicit result objects where practical. Keep compatibility methods thin.
- [ ] Remove unconditional tables, sleeps, and network access from ordinary
  unit tests. Live NIST tests must be isolated and labeled; offline fixtures are
  the release correctness gate.
  - NASA fitting no longer prints tables or R-squared lines from normal fit
    calls; callers can explicitly invoke the retained table formatter. The
    remaining ignored live tests and test-local sleeps still need inventory.
- [ ] Add a stable, local NIST/Shomate fixture that exercises record parsing,
  coefficient stitching, and independent Cp/H/S evaluation without network
  access. Keep live WebBook retrieval as a separately labeled integration test.
- [ ] Add a DBhandlers campaign script and compact report only after the
  correctness matrix is stable.

## P0: Public catalog lifecycle, audit, and telemetry

- [x] Reuse the existing immutable `ThermoRepository` lifecycle instead of
  adding a competing database cache. `ThermoData::try_new()` shares the default
  repository through `Arc`; `try_new_from_paths()` creates an isolated
  repository snapshot; `try_new_fresh()` explicitly reloads a fresh snapshot.
- [x] Expose the repository lifecycle as an ergonomic public open/use/close
  session. `close()` must release only the caller's handle and local query
  state; it must not invalidate snapshots still held by other callers or the
  process-wide default repository.
- [x] Add a separate non-transactional, per-record library audit API. It must
  retain each source record and report library, record key, stage, typed error
  category, and message while continuing to inspect later records.
  - Initial API is in `DBhandlers::library_audit`: explicit record selection or
    full-library scan; JSON report persistence; thermo and transport payload
    deserialization/parsing. Audit results are diagnostic only and never
    populate or mutate `SubsData`/production calculators.
- [x] Add full-catalog stress stories for every bundled record family. The
  audit accounts for all rows in NASA gas, NUIG thermo, condensed NASA, CEA,
  and Aramco transport; NIST is an external/dynamic source, not a bundled
  catalog family, and needs a separate offline fixture.
- [x] Add opt-in catalog-load telemetry at the repository boundary: shared
  cache hit/miss, per-file bytes/open/read/deserialize durations, schema
  validation, index materialization, total duration, and partial per-file
  failure observations. Ordinary open/lookup paths remain uninstrumented.
- [ ] Extend transactionality coverage across every production property and
  coefficient-publication route. The primary `search_substances` parse path
  now has a failure-in-the-middle rollback gate; inject failures during
  derived-value construction/publication and verify caches, closures,
  coefficients, and prior search results remain unchanged.
- [x] Preserve canonical library, record key, parsing stage, and typed source
  error when a selected production thermo/transport payload fails to load.
- [ ] Extend typed context to record resolution and NIST retrieval failures:
  report the canonical library/query or requested substance, failing stage,
  and preserve the typed lookup/parser/network source through public error and
  logging boundaries. Do not mislabel a missing record as a payload parse
  failure.
- [x] Add telemetry and audit examples showing the intended lifecycle:
  open one repository/session, perform multiple selected lookups, close the
  caller-owned session; use audit reports to inspect defects, then run clean
  selections through the ordinary transactional production APIs.
  - Runnable example: `examples/dbhandlers_catalog_audit.rs`.

## Completion criteria

- Offline fixtures prove property-level accuracy and unit consistency.
- No failed fit publishes partial coefficients or loses its typed root cause.
- The production NIST stitching path does not require symbolic nonlinear
  optimization for a coefficient-linear model.
- Live-data tests are advisory; release correctness is reproducible offline.
