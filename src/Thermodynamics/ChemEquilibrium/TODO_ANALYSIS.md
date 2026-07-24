# Chemical Equilibrium Refactoring Plan

## Scope and decisions

This checklist covers `Thermodynamics/ChemEquilibrium` only. Its goal is one
maintainable equilibrium engine with explicit numerical policies, typed failure
reporting, independent validation, and deterministic tests.

- [x] Treat the logarithmic `Chem_eq_K_eq*` formulation as the candidate
  source of truth for equilibrium calculations.
- [x] Treat the legacy classical equilibrium stack as a migration/reference path, not as a
  second production architecture where new equilibrium logic should be added.
- [x] Keep hand-written nonlinear solvers only as explicit temporary fallback
  backends. They must not silently define the public equilibrium contract.
- [x] Use the reaction-equilibrium-constant approach as an independent
  validation solver for systems within its documented applicability range.
- [x] Defer the final `phase_*` bridge and GUI until the core problem, result,
  validation, and solver-policy contracts are stable.
- [x] Replace the misleading `Chem_eq_K_eq*` module family with names that
  reflect its actual responsibilities: log-moles formulation, workflows, and
  nonlinear support.

## Non-negotiable contracts

A nonlinear backend reporting convergence is not sufficient for accepting an
equilibrium result. Every accepted result must pass the same backend-independent
validation gate.

- [x] Require finite positive temperature and pressure and finite input data.
- [x] Require a dimensionally consistent species list, initial composition,
  elemental-composition matrix, phase assignment, and thermochemical vector.
- [x] Define the contract for zero and trace species in log space. Do not hide
  clipping, floors, or exponent saturation inside residual evaluation.
  The phase-control helpers now use a named trace-floor constant instead of a
  hidden numeric literal, and the log-mole seed path stays centralized in the
  typed initial-guess policy.
- [x] Reject NaN/Inf residuals, Jacobians, iterates, and reconstructed mole
  numbers with typed errors.
- [x] Check non-negative mole numbers, elemental conservation, normalized
  nonlinear residuals, and equilibrium optimality before publishing a result.
  - [x] The candidate report now exposes both L2 and RMS residual norms, so
    callers can reason about normalized residual evidence without recomputing
    it from scratch.
  - [x] Reaction-affinity optimality now has its own tolerance in the
    backend-independent acceptance gate, instead of being only reported
    passively after the solve.
  - [x] The validation gate now rejects candidates whose reaction-affinity
    block exceeds its tolerance, even when the backend-facing residual is
    small.
- [x] Keep species and reaction ordering explicit and deterministic. The
  validated problem and prepared snapshot now preserve canonical species,
  phase, and reaction ordering, and the golden fixtures plus solver-report
  tests lock that behavior down.
- [x] Never publish partial state from a failed setup, solve, fallback attempt,
  temperature sweep, or phase-control restart.
  - [x] Legacy matrix/system setup now validates the numerical settings before
    staging, so malformed tolerances or budgets cannot publish a partly rebuilt
    reaction basis.
  - [x] Single-temperature solve now stages automatic seeds and reconstructed
    mole/result bundles locally, publishing them only after an accepted solve.
  - [x] Sequential and parallel sweep paths now validate their complete
    column-oriented mole table before replacing any published sweep results.
  - [x] A failed solve after a previously accepted solve keeps the published
    solution, mole state, validation report, and solve report intact.

## P0 - Establish a safe canonical core

### P0.1 Characterize the current behavior

- [x] Inventory the public entry points in `Chem_eq_K_eq*`,
  `easy_equilibrium.rs`, and the retired classical equilibrium stack.
  A current behavior snapshot now lives in `behavior_inventory.md` and
  separates canonical modules, transitional modules, and their test families.
- [x] Record which current tests are physical regression fixtures, solver unit
  tests, debug experiments, or duplicates.
  Canonical regression and solver-policy tests live in the `equilibrium_*`
  modules, low-level numerical characterization lives in
  `equilibrium_nonlinear.rs`, legacy/classical coverage stays in the
  the retired classical stack and `NR_Legacy.rs` suites, while retained
  `debug_*` studies remain explicitly experimental. `Untitled-1.rs` has been
  removed.
  - [x] Audit the retired classical test family before deletion.
    The migration notes now separate adopted canonical contracts, P4
    phase-bridge scenarios, and tests tied only to obsolete state.
- [x] Preserve representative currently solvable systems as golden regression
  fixtures before changing numerical backends. The representative O2/O and
  N2/N fixtures now live in `equilibrium_golden_fixtures_tests.rs` and the
  RST matrix tests, so backend changes immediately trip a stable regression.
- [x] Record difficult and currently failing systems separately. A known
  failure must not be converted into a weak passing assertion.
  - [x] Hard fallback-prone fixtures are listed in `behavior_inventory.md` so
    future refactors can keep them distinct from easy acceptance regressions.
  - [x] No canonical equilibrium failure is being softened into a passing
    assertion; retained hard cases remain explicit regression fixtures.
- [x] Identify all `panic!`, `assert!`, `expect`, unchecked indexing, and
  production `unwrap` paths in setup, solving, phase control, and output.
  The next audit pass should focus on this boundary now that the public
  surface inventory is explicit.
  Current hotspots are concentrated in the retained legacy/classical modules;
  the canonical `equilibrium_*` path is already much more typed and its
  remaining `unwrap` / `assert` hits are in tests or legacy-only code.

### P0.2 Typed domain boundary

- [x] Introduce a validated `EquilibriumProblem` containing species identity,
  conditions, initial composition, elemental constraints, thermochemical data,
  phase/activity information, and numerical formulation options.
- [x] Introduce an immutable `PreparedEquilibriumProblem` for matrices,
  reaction bases, element totals, and phase-derived data.
- [x] Add residual scaling to the prepared snapshot as an explicit typed
  contract once its units and backend contract are defined. The prepared
  boundary now owns `ResidualScalingContract`; variable scaling remains a
  separate future design because it is a different coordinate transform.
- [x] Introduce `EquilibriumSolution`; accepted legacy solves expose an
  immutable snapshot instead of requiring callers to combine mutable fields.
- [x] Introduce `EquilibriumSolveReport` with policy, ordered attempts, and
  accepted backend details. RST attempts additionally expose termination,
  iteration/evaluation counters, linear solves, and elapsed time; comparable
  legacy counters remain unavailable until that temporary adapter is replaced.
  The report now also exposes summary helpers so UI/tests can reason about
  fallbacks and skipped attempts without hand-walking the raw attempt vector.
- [x] Introduce typed identifiers for species, elements, phases, and reactions
  where raw indices can be confused. Typed boundary accessors now exist on the
  validated problem and prepared snapshot; the internal ordered arrays are
  still the storage format for now.
- [x] Separate domain settings from backend-specific solver settings.
  `EquilibriumLogMoles` now owns a dedicated `EquilibriumSolverSettings`
  bundle, while the domain/problem boundary stays focused on chemistry and
  composition data.
- [x] Replace public mutation sequences such as `set_problem -> create_* ->
  solve` with validated construction plus one explicit solve operation.
  - [x] The canonical `solve_problem_with` configuration hook now receives
    only `EquilibriumSolverSettings`, so post-validation tuning cannot mutate
    chemistry, phase layout, conditions, or published state.
  - [x] Public workflow constructors that still need a mutable sweep engine now
    route through `EquilibriumProblem` + `EquilibriumLogMoles::from_problem`
    instead of rebuilding the solver by staging mutable fields step by step.

### P0.3 Typed errors

- [x] Split `ReactionExtentError` into meaningful categories: invalid input,
  thermochemical lookup, formulation, residual/Jacobian evaluation, backend
  failure, non-convergence, invalid candidate solution, validation mismatch,
  unsupported model, and all-backends-failed. RST linear solve, singular
  Jacobian, and numerical-breakdown errors now have a typed backend-failure
  category instead of being mislabelled as legacy evaluation errors. Canonical
  Gibbs-cache, residual-generator dimension, and condition failures no longer
  use the catch-all `Other` variant. The `Other` variant has been removed from
  the public error enum; the remaining work is to make every category carry
  enough domain identity for a user to locate the bad input.
  - [x] `ReactionExtentErrorKind` now exposes a machine-readable top-level
    classification for the requested failure families, so tests and policy code
    can branch on categories instead of parsing display text.
- [x] Preserve error sources instead of converting them to debug strings in
  every backend adapter and legacy compatibility path.
  - [x] `ReactionExtentError` and `SolveError` now implement stable `Display`
    and `std::error::Error` contracts. Nested `SubsDataError`/`SolveError`
    values remain available through `source()`.
  - [x] Solver-cascade attempt reports now retain stable `Display` diagnostics
    rather than private `Debug` formatting.
  - [x] Recoverable attempt failures carry a typed report category so callers
    do not need to infer retry diagnostics by parsing a message.
- [x] Include species/reaction identifiers and offending values in diagnostics.
  - [x] Canonical residual/Jacobian condition checks now identify the exact
    invalid scalar; non-finite Gibbs values identify the species index and
    temperature; invalid phase totals identify the phase index and value.
  - [x] Invalid equilibrium candidates now carry the offending field name in
    the error itself, so setup/acceptance failures are no longer just opaque
    "candidate rejected" messages.
- [x] Make invalid input and unsupported physics non-retryable. The retry
  boundary now has typed helpers on `ReactionExtentError`, and the cascade
  helper explicitly classifies only numerical backend failures as retryable.
- [x] Return every backend attempt in `AllBackendsFailed`, rather than only the
  final hand-written solver error.
  - [x] A cascade that reaches a non-retryable error after earlier attempts now
    returns `CascadeAborted { attempts, cause }`; ordinary exhaustion remains
    `AllBackendsFailed`, so neither case loses its diagnostic trace.

### P0.4 One logarithmic formulation

- [x] Reject non-square systems before closure construction: independent
  reaction equations plus elemental-balance equations must equal the number of
  log-mole unknowns.
- [x] Move residual, Jacobian, scaling, and mole reconstruction into pure
  formulation functions over `PreparedEquilibriumProblem`. The canonical
  residual and Jacobian now have pure `evaluate_*` entry points, while legacy
  closure/wrapper APIs delegate to the same formulas. Mole reconstruction,
  immutable accepted-result packaging, and explicit row scaling now also live
  on the prepared boundary. Variable scaling remains a separate coordinate
  transformation and is intentionally not implied by row scaling.
- [x] Ensure every numerical backend receives exactly the same residual and
  Jacobian contract.
  - [x] The offline formulation fixture now checks every analytic Jacobian
    entry against the derivative of the corresponding symbolic residual.
- [x] Remove duplicate closure construction from single solve, parallel solve,
  temperature sweep, and phase-control paths.
  - [x] Single-solve and parallel-solve Gibbs closure construction now share
    one ordered lookup helper, so the closure ordering contract is no longer
    duplicated in two separate code paths.
  - [x] Temperature-sweep and workflow constructors now reuse the same shared
    Gibbs-function assembly path instead of hand-building another closure
    ordering layer.
- [x] Validate the analytic Jacobian against finite differences across normal,
  trace-species, high-temperature, and low-temperature states. The regression
  matrix now covers several fixed temperature/state pairs on the canonical
  typed problem; richer real thermochemistry temperature fixtures remain a
  separate future layer.
- [x] Document the meaning and units of every residual block and scale factor.
- [x] Guard `exp(log_n)` reconstruction against overflow and underflow without
  silently accepting a distorted solution.

### P0.5 Backend-independent acceptance gate

- [x] Add a pure `validate_equilibrium_candidate` function.
  The acceptance gate now takes typed residual data plus typed acceptance
  criteria, rather than a bundle of unstructured vectors and bare tolerances.
- [x] Report absolute and scaled residual norms separately.
- [x] Report maximum elemental-balance error by element.
- [x] Report minimum mole number and all non-finite values.
- [x] Reject truncated or empty raw/acceptance residual vectors before their
  norms can be interpreted as convergence evidence.
- [x] Check reaction affinity or chemical-potential optimality appropriate to
  the implemented model.
  The acceptance gate now reports the reaction-affinity block separately from
  the elemental-balance block, so the physical optimality evidence is explicit
  and can be compared without reusing the wrong residual rows.
- [x] Distinguish `BackendConverged` from `SolutionAccepted` in reports.
- [x] Publish solution/state only after this gate succeeds.

### P0.6 Isolate legacy nonlinear solvers

- [x] Move `LMSolver`, `NRSolver`, and `TrustRegionSolver` behind one temporary
  legacy adapter. `equilibrium_legacy_backend` now owns their mechanical
  dispatch; policy, budgets, validation, retry, and publication remain in the
  canonical orchestration path.
- [x] Remove panics, assertions on user data, global logger initialization, and
  direct console output from legacy production paths. The canonical log-moles
  engine no longer initializes a process-global logger, and the legacy trust
  region solver now exposes a singular LU step as `SolveError::SingularMatrix`
  instead of silently substituting a zero step. The public Newton step-limit
  helper now returns typed input errors instead of asserting vector lengths.
  Remaining console output and legacy solver diagnostics still require
  isolation.
  - [x] Canonical sweep-table construction is now pure (`moles_table`); the
    solver and workflow APIs no longer print directly to stdout.
  - [x] Correct the Trust Region quadratic-model reduction to use `J^T f` for
    its linear term. The previous `f^T p` expression was mathematically wrong
    and caused ill-conditioned systems to shrink the trust radius to machine
    noise; the regression fixture now converges deterministically.
  - [x] Canonical `equilibrium_*` runtime paths are now free of production
    `panic!/unwrap/expect` sites; the remaining occurrences are confined to
    tests, doc examples, or the legacy/classical layer that is explicitly out
    of the canonical path.
- [x] Do not select legacy solvers implicitly outside an explicit solver policy.
  The default no longer falls back to the handwritten cascade; legacy
  backends remain available only through an explicit policy.
- [x] Add characterization tests before changing their behavior.
  Canonical regression, backend-matrix, and fallback-story fixtures now pin
  the current behavior before more structural changes land.
- [x] Define the legacy-backend retention rule: the handwritten backend stays
  as an explicit compatibility/fallback path, even after the RST matrix covers
  the main regression set. The remaining work is to keep its policy boundary
  strict and its diagnostics typed.
  - [x] Review the retained handwritten solvers as production-quality
    fallbacks, not as disposable demo code. `LMSolver`, `NRSolver`, and
    `TrustRegionSolver` should keep explicit step control, bounded failure
    modes, and deterministic test coverage even though they are legacy.
  - [x] Remove `panic!` / `assert!` / silent `None` publication paths from
    `NR_Legacy.rs` where a typed failure or an explicit best-effort step can
    be returned instead.
  - [x] Decide whether the two currently failing legacy regression tests are
    strict contract tests or characterization tests. If the solver is allowed
    to stop without convergence for a given fixture, encode that explicitly in
    the test instead of relying on an implicit panic or unwrap.
  - [x] Audit the inner step control in `LMSolver`, `NRSolver`, and
    `TrustRegionSolver` for finite-value checks, dimension checks, and clear
    iteration-exhaustion reports.

## P1 - Reliable solver architecture

### P1.1 Common backend adapter

- [x] Define an internal `EquilibriumNonlinearBackend` interface receiving a
  prepared problem, initial iterate, residual, Jacobian, and attempt budget.
  The adapter now owns one mechanical request/result boundary for both legacy
  and RST backends.
- [x] Return a typed `SolverAttemptReport` with backend name, termination
  reason, iteration/evaluation counts, norms, elapsed time, and candidate.
- [x] Keep RustedSciThe types inside adapter modules so domain types do not
  depend on one numerical library's API. The prepared symbolic problem is now
  carried as an opaque adapter-level wrapper instead of leaking the concrete
  RST problem type into the domain layer.
- [x] Extend the typed problem boundary with immutable symbolic thermochemical
  expressions (or a dedicated symbolic formulation snapshot) before making an
  RST-symbolic policy the default for every `EquilibriumProblem`. The current
  typed input intentionally stores numeric Gibbs closures only. The adapter
  now builds a dedicated immutable symbolic thermochemistry snapshot before it
  assembles the RST problem, so the lookup pipeline is no longer intertwined
  with the symbolic bridge itself.
- [x] Support deterministic fake backends for cascade tests.

### P1.2 RustedSciThe migration

- [x] Audit the current RustedSciThe 0.4.12 `Nonlinear_systems` API through its
  `NonlinearProblem`, `JacobianProvider`, `SolverEngine`, and typed result APIs.
- [x] Implement symbolic RST adapters for ordinary LM, MINPACK LM, Nielsen LM,
  trust-region LM, Powell dogleg, and damped Newton. The canonical adapter
  passes residual `Expr` plus temperature/settings to
  `SymbolicNonlinearProblem`; RST owns Lambdify/Jacobian preparation.
- [x] Benchmark and characterize candidates before fixing the default order.
  Do not choose a default solely from method names or one easy fixture.
  - [x] Current RST default order is now characterized by a dedicated matrix
    test, so the chosen order is explicit and regression-protected even before
    any future performance-driven reordering.
- [x] Map RST termination reasons without losing details.
- [x] Verify scaling, bounds, and stopping tolerances have the same meaning at
  the domain and backend boundaries.
  - [x] RST solve-contract construction is now centralized as a typed
    `RustedSciTheSolveContract`, so tolerance and iteration budget are no
    longer scattered across call sites. Feasibility bounds remain a
    formulation-side concern and are not silently implied by the backend
    contract.
- [x] Use the RST symbolic Lambdify backend for the first migration; keep AOT
  out until a measurement shows a material evaluation bottleneck.

### P1.3 Explicit solver cascade

- [x] Replace the loop inside `solver_impl` with a typed `SolverPolicy`.
- [x] Support at least `Single(backend)` and `Cascade(Vec<backend>)`; provide a
  documented `Auto` policy only after the backend matrix is measured.
- [x] Give every attempt an explicit iteration budget and enforce a global
  cascade budget. `SolverCascadeBudget` caps started backends, per-backend
  iterations, and total allocated iterations; unstarted backends are retained
  as ordered `Skipped` entries in the report. Evaluation budgets still depend
  on support from the numerical backend API.
- [x] Restart each backend from the original validated initial guess unless a
  named warm-start policy explicitly permits reuse of a previous candidate.
  The cascade now documents and tests this invariant; no implicit candidate
  reuse is allowed.
- [x] Retry only recoverable numerical failures such as non-convergence,
  singular steps, or rejected candidate solutions.
- [x] Never retry invalid input, missing thermochemical data, dimension errors,
  or unsupported activity/phase models.
- [x] Run the common acceptance gate after every backend success. If the
  candidate is invalid, record why and continue the cascade.
- [x] Preserve deterministic backend order and expose all attempts in the final
  solve report.
- [x] If several valid candidates are retained, choose between them only by a
  documented physical/numerical criterion, never by "first finite vector".
  A dedicated comparator now ranks accepted candidates by residual norm,
  raw residual norm, balance error, and minimum mole evidence. The current
  cascade still stops on the first accepted backend, so this criterion is now
  explicit and testable for any future multi-candidate retention path.
- [x] Make logging a view of the report, not the only record of fallback.
  - [x] Solve reports now expose typed summary accessors for accepted,
    fallback, skipped, and ordered attempt views.
  - [x] The main solve path now logs the typed solve report summary instead of
    a separate ad hoc fallback message.

### P1.4 Initialization, scaling, and continuation

- [x] Centralize construction and validation of log-mole initial guesses. The
  typed constructor, solver settings, and legacy mutable setter now share the
  same finite-value and dimension validation; temperature-sweep/phase-control
  seed paths now reuse the same trace-seed policy instead of re-encoding the
  floor contract in multiple places.
- [x] Define a deterministic trace-species floor contract. A named absolute
  coordinate default now replaces hidden duplicated literals, and a typed
  system-scale-relative policy is available when callers want scale-aware
  seeding explicitly. The production default remains the absolute floor until
  we finish the numerical characterization for a relative default.
- [x] Test row scaling and variable scaling independently. Row scaling now has
  typed dimension/scale validation and direct residual/Jacobian tests.
  Variable scaling now has its own typed coordinate-scaling contract and
  tests, but it remains intentionally separate from the solver path; the
  remaining work is integration, not conceptual separation.
- [x] Treat temperature continuation as an explicit policy with transactional
  per-temperature results and a clear warm-start contract. Sequential and
  chunked sweeps now use `ContinuationSeedPolicy`: either reuse the previous
  accepted solution or restart every point from the configured seed. The fully
  parallel `par2` path is explicitly independent because no ordered prior
  point exists while tasks run concurrently. All sweep paths publish
  `TemperatureSolveSnapshot` evidence for accepted points.
- [x] Make legacy mutable `set_problem` and explicit-mole-map updates validate
  before publishing state. They now reject invalid pressure, non-finite or
  negative amounts, out-of-range/duplicate phase indices, and unknown mapped
  substances without clearing a previously published solution.
- [x] On continuation failure, record the failed temperature and attempts;
  never silently leave vectors with different lengths or stale state. Sweeps
  publish `TemperatureSolveFailure` records with the temperature, diagnostic
  message, and cascade trace when available, alongside separate accepted-point
  snapshots.
- [ ] Add multi-start only if regression evidence justifies its complexity.

## P2 - Independent equilibrium-constant validation

The K-equilibrium validator is valuable precisely because it is a second
formulation. It may share immutable thermochemical input and common final
invariant checks, but it must not reuse the main residual/Jacobian assembly.

### P2.1 Reaction-basis contract

- [x] Make `ReactionBasis` typed and validate its dimensions, rank, nullspace
  accuracy, species ordering, and elemental conservation.
  - The independent validator now owns `ValidatedReactionBasis`, which binds
    the species order to a finite species-by-reaction matrix and validates
    dimensions, declared rank, reaction count, and `A^T * N`. Replacing the
    canonical solver's older public `ReactionBasis` remains separate work.
- [x] Detect underdetermined, overdetermined, and numerically rank-ambiguous
  bases with typed errors.
- [x] Make basis normalization/sign conventions deterministic so comparisons
  and fixtures are stable.
  - Explicit validator bases now normalize each reaction by its first
    significant species coefficient and orient that coefficient as a
    reactant. Multi-dimensional SVD basis canonicalization is still open.
- [x] Test basis invariance under species and element permutations.

### P2.2 Independent K_eq solver

- [x] Implement a distinct reaction-extent/K_eq problem and solver path for
  small systems with a known independent reaction basis.
  The first pass now solves one-reaction systems in extent space with a
  safeguarded Newton/bracketing loop and returns its own typed solve report
  plus an independent validation report.
- [x] Compute `ln(K)` and reaction quotients through code independent from the
  main log-moles residual builder.
  `EquilibriumConstantProblem` evaluates standard reaction Gibbs energies,
  dimensionless ideal-gas activities, and `ln(Q) - ln(K)` without calling the
  canonical residual/Jacobian implementation.
- [x] Define supported activity models and phase combinations explicitly.
  The first contract supports one ideal-gas phase only; condensed and mixed
  phases must become explicit variants rather than implicit corrections.
- [x] Return `ValidationNotApplicable` for unsupported or excessively large
  systems; do not report a false solver failure.
  The first-pass solver refuses multi-reaction systems explicitly instead of
  pretending that they failed numerically.
- [x] Give the validation solver its own typed report and numerical tolerances.
  Candidate assessment now returns ordered per-reaction `ln(Q)`, `ln(K)`, and
  residual evidence plus the maximum residual and acceptance decision. The
  future extent solver will add iteration/backend evidence to this report.
- [x] Keep it optional in production (`ValidationMode::Off/WhenApplicable/
  Required`) and enabled broadly in tests.
  The new solver already exposes an explicit on/off/applicable policy, and
  `solve_if_applicable` now returns `None` immediately in `Off` mode while the
  solver tests cover applicable, `Off`, and `Required`-mode paths. The
  top-level solve boundary now threads the validation mode and publishes an
  optional independent K_eq status alongside the main accepted solution.

### P2.3 Cross-validation report

- [x] Compare main and K_eq solutions by species moles/mole fractions with
  absolute, relative, and trace-species-aware tolerances.
  The new cross-validation report compares deterministic species moles and
  mole fractions and records the largest disagreement.
- [x] Compare elemental balances independently for both solutions.
  Canonical balance evidence remains in the canonical candidate report, while
  the independent validator remains logically separate.
- [x] Evaluate `ln(Q) - ln(K)` for every independent reaction.
  The independent validation report carries the per-reaction residuals, and
  the comparison layer exposes the largest one as part of the bridge report.
- [x] Compare total Gibbs energy or another documented equilibrium objective.
  The bridge report now computes a total ideal-gas Gibbs comparison directly
  from the thermochemical inputs and mole numbers.
- [x] Report the largest disagreement with species/reaction identity.
  Each comparison row carries the species name, index, canonical value, K_eq
  value, and both absolute deltas.
- [x] Distinguish main-solver failure, validator failure, not-applicable, and
  genuine cross-validation mismatch.
  The report layer is in place, and a typed status helper now classifies
  canonical failure, validator failure, not-applicable, and compared-result
  branches separately. The top-level solve path now routes and publishes the
  resulting status instead of collapsing those cases into one branch.

### P2.4 Validation fixture matrix

- [x] Add small analytically tractable dissociation/association systems.
  - Offline contracts now cover `A2 <=> 2A` with pressure dependence and
    `2NO <=> N2 + O2` with the closed-form extent inherited from the old test
    module, now exercised through the independent extent solver.
    The inherited `2N2O <=> 2N2 + O2` and `2NO2 <=> N2 + 2O2` cases now also
    cover three-species reconstruction and the nonzero pressure exponent.
    Thermochemical-library association fixtures remain open.
- [x] Add water-gas and hydrogen/oxygen equilibrium fixtures.
  - The cross-validation matrix now includes a water-gas-shift fixture in
    addition to the earlier diatomic and nitric-oxide contracts.
- [x] Add methane/air lean, stoichiometric, and rich fixtures.
  - The methane-combustion baseline, lean, and rich cases are now covered in
    the independent validation matrix.
- [x] Add inert dilution and pressure-variation fixtures.
  - A hydrogen/oxygen fixture with inert `N2` now exercises both pressure
    dependence and diluent handling in the independent validator bridge.
- [x] Add low-, medium-, and high-temperature sweeps.
  - The methane-combustion cross-validation matrix now exercises a
    low/medium/high temperature sweep at fixed composition.
- [x] Compare the canonical solver and K_eq validator for every applicable
    fixture, not only for one final scalar.
- [x] Keep all validation fixtures offline and deterministic.

## P3 - Test and reliability matrix

### P3.1 Formulation unit tests

- [x] Test residual blocks, Jacobian blocks, scaling, log/mole conversion, and
  phase totals independently.
- [x] Test analytic versus finite-difference Jacobians over a state matrix.
- [x] Test NaN/Inf, overflow, underflow, zero totals, duplicate species,
  unassigned species, singular composition matrices, and malformed phases.
  The public phase-index boundary now returns a dimension error instead of
  indexing past the phase map; broaden the malformed-phase matrix further.
- [x] Test exact dimensions and ordering at every boundary.

### P3.2 Backend matrix

- [x] Run every supported RST backend against the same curated fixture set.
  The first `O2/O` case is covered, high-temperature `N2/N` provides a real
  rejected-candidate -> RST fallback story, and the dilute inert-diluent
  matrix now exercises the same accepted-backend reporting contract.
- [x] Record convergence, acceptance, iterations, residuals, balance errors,
  and runtime without making runtime assertions flaky.
- [ ] Keep legacy backend parity tests only until its deletion criteria pass.
- [x] Include difficult cases that require fallback, not only easy systems on
  which every solver succeeds immediately.

### P3.3 Cascade story tests

- [x] First backend fails recoverably; second succeeds and is accepted.
- [x] First backend returns a nominal success with an invalid candidate; the
  acceptance gate rejects it and the next backend succeeds.
- [x] Invalid input produces zero backend attempts.
- [x] Every backend fails; the error contains every attempt in exact order.
- [x] Global budget exhaustion stops the cascade deterministically.
- [x] A failed attempt cannot mutate the original problem, initial guess,
  published solution, or temperature-sweep results. The RST rejected-candidate
  fallback story now verifies that the original seed and initial composition
  remain unchanged; failed numerical attempts and temperature sweeps remain.
- [x] Warm-start behavior occurs only when explicitly selected.

### P3.4 Physical and metamorphic tests

- [x] Check non-negative concentrations and elemental conservation for all
  accepted fixtures.
- [x] Check invariance under species, element, and reaction reordering.
- [x] Check consistent results when all initial mole totals are scaled.
- [x] Check robustness to reasonable initial-guess perturbations.
- [x] Check inert-species addition and removal where the model predicts it.
- [x] Check temperature/pressure trends against known qualitative behavior.

### P3.5 Integration and regression tests

  - [x] Build problems from real offline `SubsData` thermochemistry and retain
    provenance in reports.
    Local `SubsData` thermochemistry is the normal path; NIST fallback should
    remain a rare and explicitly reported fallback path, not the default.
  - [x] Cover mixed-library thermochemical inputs without network access.
- [x] Add transactional tests for setup, solve, phase-control restarts, and
  serial/parallel temperature sweeps.
- [x] Replace debug binaries/scripts with assertions in named test modules or
  documented examples.
  The standalone debug comparison file has been removed, and no separate
  debug study source files remain in this subsystem.
- [x] Separate fast unit tests, deterministic integration tests, and explicitly
  ignored expensive studies.

### P3.6 Retained ideas from the classical implementation

The retired classical equilibrium modules are not a compatibility target. The
following list records the useful contracts discovered during their deletion
audit so the implementation can be removed without losing sound ideas.

Technical work that is independent from the future phase bridge:

- [x] Add a pure typed `SpeciesCapacityReport` derived from elemental totals:
  `n_i_max = min(b_e / a_ie)` over elements present in species `i`.
  - Do not copy the legacy `1.2` safety multiplier or string-keyed solver
    bounds.
  - Use the exact capacity for candidate diagnostics, initial-guess checks,
    impossible-composition detection, and composition-ordering tests.
  - Keep log-moles positivity as the production feasibility mechanism; this
    report is evidence and validation, not a second hidden solver policy.
- [x] `EquilibriumProblem::validate` now rejects species rows without any
  positive elemental support, so impossible compositions are caught before
  preparation.
- [x] Add optional, pure `FormulationDiagnostics` for prepared problems.
  Include numerical rank, singular values, a documented condition estimate,
  and suspicious null directions where available.
  - Compute expensive SVD diagnostics only during preparation or by explicit
    request, never on every nonlinear iteration.
  - Return typed data for reports/tests; do not print directly or let a
    diagnostic heuristic silently reject an otherwise valid solution.
  - The prepared boundary now exposes `formulation_diagnostics(tolerance)` and
    `preview_with_diagnostics(tolerance)` without changing solve behavior.
- [x] Add a typed `EquilibriumProblemPreview` (or equivalent report rows) for
  CLI, GUI, examples, and snapshot tests. It should expose conditions, species
  and element ordering, element matrix, reaction basis, initial inventory,
  scaling contract, solver policy, thermochemical provenance, and optional
  formulation diagnostics without owning any solving behavior.
  - The current preview exposes the validated problem state, reaction basis
    size, exact element totals, and the species-capacity evidence that can be
    derived without the future phase bridge.
  - The preview also now exposes stable summary rows and `Display` output for
    CLI/snapshot consumers without adding any solving behavior.
  - `SpeciesCapacityReport` and `FormulationDiagnostics` now also have stable
    `Display` output for human-readable reports and regression snapshots.
- [x] Give the independent K_eq cross-validation report the same typed
  summary/`Display` treatment so comparisons can be inspected without ad hoc
  string assembly.

Ideas intentionally deferred until P4:

- [ ] Validate named initial compositions transactionally at the
  `ResolvedPhaseSystem -> EquilibriumProblem` adapter: every component must be
  known, represented exactly once by its typed component id, and aligned with
  deterministic solver ordering. Do not restore the legacy
  `HashMap<Option<String>, ...>` contract.
- [ ] Reconsider a small independent Lagrange-stationarity validator only if
  multi-reaction regression evidence shows a real validation gap after the
  phase bridge exists.
  - Its possible value is mathematical independence from the reaction-basis
    residual and applicability beyond the one-reaction extent validator.
  - Do not migrate the legacy mutable `Solver`, auxiliary `Np` unknowns, raw
    mole variables, or unscaled equations. Any future experiment must use
    log-moles, typed phase/component ids, scaled residuals, and a typed report.

Classical ideas already superseded by the canonical engine:

- [x] SVD reaction-basis discovery and elemental-conservation validation.
- [x] Temperature-band discovery and thermochemical coefficient refresh.
- [x] Explicit sequential warm-start/continuation policy.
- [x] Ordered nonlinear backend cascade, acceptance gate, budgets, and attempt
  reports.
- [x] Independent reaction `delta G`, `ln(K)`, quotient, and extent
  reconstruction through the K_eq validation subsystem.

Classical behavior that must not be migrated:

- [x] Global Lagrange/Newton interpolation followed by value clamping. Global
  polynomial oscillation and clamping can hide invalid composition and break
  elemental balances; future plotting interpolation belongs to a separate
  postprocessing layer with explicit error bounds.
- [x] Parallel mutable maps, placeholder zero closures during clone, string
  variable names such as `N0`/`Lambda0`, direct `println!`, and unchecked
  `unwrap`-driven output.
- [x] Auxiliary phase-total unknowns when phase totals can be derived from the
  ordered species mole vector.
- [x] The classical solver cascade and residual-norm-only acceptance rule.

## P4 - Unify the two equilibrium worlds

This phase begins only after P0-P3 contracts are stable.

The first production target is a closed reacting system at fixed pressure and
temperature. The supported phase-model slice is deliberately narrow:

- one ideal-gas phase containing any number of gas components;
- any number of one-component pure liquid/solid/condensed phases;
- no ideal/non-ideal solution phase and no multiple gas phases until those
  activity models have explicit equations and validation fixtures.

This is not merely an adapter task. `ResolvedPhaseSystem` already has the
correct phase-qualified component identity, while `EquilibriumProblem` still
uses unique bare strings and integer-only phase records. The domain boundary
must be corrected before the two systems are connected.

### P4.0 Freeze the fixed-P,T multiphase contract

- [x] Document the thermodynamic ensemble as closed-system Gibbs minimization
  at fixed `P`, `T`, and conserved elemental totals.
- [x] Define the supported-model matrix in code and documentation:
  `IdealGas` maps to the ideal-gas activity law and `PureCondensed` maps to
  unit activity.
- [x] Require a `PureCondensed` phase to contain exactly one component. A
  multi-component condensed phase is a solution and must return
  `UnsupportedModel`, not silently reuse the ideal-solution equation.
- [x] Reject more than one ideal-gas phase until a physically meaningful
  immiscible/multiple-gas model exists.
- [x] Keep physical state separate from activity model. `Liquid`, `Solid`, and
  `Condensed` select records and describe output; `PhaseModel` selects the
  chemical-potential equation.
- [x] State explicitly that zero initial amount means "candidate phase absent
  initially", not "phase excluded from equilibrium". Exclusion must be a
  separate typed input choice.
  `equilibrium_multiphase_domain` now makes this distinction at the physical
  input boundary; a later `PhaseSet` bridge owns explicit exclusion.

### P4.1 Make phase-qualified identity canonical in the equilibrium problem

- [x] Introduce an equilibrium component descriptor that retains
  `PhaseComponentId`, bare substance name, physical state, phase model, and a
  stable display label. The bare name may repeat across phases; the qualified
  component id may not.
  `equilibrium_component::EquilibriumComponentDescriptor` now derives the
  bare name and stable label from the qualified id and carries the exact
  physical-state, phase-model, and solver activity-model mapping.
- [x] Replace `EquilibriumProblem::species: Vec<String>` as the canonical
  identity with the ordered component descriptors. Keep labels as a derived
  view only.
  `EquilibriumProblem` now owns ordered `EquilibriumComponentDescriptor`
  values; `species()` is a derived label view retained only for numerical
  compatibility while the dense legacy formulation is migrated.
- [x] Replace or rename the integer-only equilibrium `PhaseId` so it cannot be
  confused with the phase subsystem's semantic `PhaseId`. Use a typed phase
  index internally and preserve the semantic id at the boundary.
  The dense solver id is now `equilibrium_ids::PhaseIndex`; bridge descriptors
  retain `phase_layout::PhaseId` as the semantic identity.
- [ ] Store an ordered equilibrium phase descriptor instead of the current
  bare `Phase { kind, species: Vec<usize> }`; derive index ranges/maps once
  from `SystemLayout`.
  - [x] `EquilibriumPhaseDescriptor` is now the bridge-level source of truth.
  - [ ] Remove the legacy numeric `Phase` projection from
    `EquilibriumProblem` after P4.3 builds the complete numerical payload.
- [x] Permit the same chemical substance in multiple phases while continuing
  to reject duplicate `PhaseComponentId` values.
- [x] Carry a layout fingerprint/revision in the bridge result so a solution
  cannot be applied to a different resolved component order.
  `MultiphaseInitialComposition` is now fingerprint-bound to its layout; the
  future solver bridge/result uses the same fingerprint.

Required tests:

- [x] `gas::H2O` and `liquid::H2O` become two solver unknowns with one shared
  molecular formula and two distinct thermochemical records.
- [x] Permuting input map insertion order does not change component, phase,
  matrix, residual, or result order.
  The bridge regression builds the same local NASA gas/condensed system from
  opposite `HashMap` insertion orders and compares metadata, descriptors,
  labels, element matrix, initial coordinates, and `G0(T)` ordering.
- [ ] Duplicate qualified ids, duplicate phase ids, and mismatched phase ranges
  fail before thermochemistry or a nonlinear backend is invoked.

### P4.2 Add a typed initial-composition boundary

- [x] Introduce `MultiphaseInitialComposition` aligned to `SystemLayout`, with
  finite non-negative physical mole numbers and at least one positive amount.
- [x] Provide a strict sparse constructor keyed by `PhaseComponentId`. Unknown
  keys are errors; omitted known components become physical zero only through
  an explicit sparse-input policy.
- [x] Reject bare-name maps when one substance occurs in more than one phase.
  Never guess whether `H2O` means gas or liquid.
  The typed boundary accepts only dense layout-order input or sparse
  `PhaseComponentId` input; there is deliberately no bare-name constructor.
- [x] Compute conserved element totals from physical initial amounts before
  trace floors are introduced. Numerical log-mole seeds must not create mass.
- [x] Convert physical zeroes to log-coordinate seeds through the existing
  typed trace policy only after the physical inventory has been validated.
  The bridge computes its element totals from physical input before creating
  `LogMolesInitialGuess`; the local NASA-gas solve verifies that a trace-seeded
  zero `H2O` leaves the `H=4`, `O=2` inventory exactly unchanged.
- [x] Make request construction transactional: failed composition validation
  must not publish a partially prepared problem, cache, or lookup report.
  `PhaseEquilibriumBuildRequest::new` validates the resolved layout,
  composition fingerprint, supported-model policy, and immutable metadata
  before returning the request; it mutates neither subsystem.

Required tests:

- [x] Dense and sparse constructors produce the same ordered vector.
- [x] Ambiguous bare names, negative/NaN/Inf amounts, unknown components, and
  an all-zero inventory return typed errors.
- [x] Trace seeding leaves computed elemental totals bit-for-bit unchanged.

### P4.3 Build one `ResolvedPhaseSystem -> EquilibriumProblem` adapter

- [x] Add the adapter inside `ChemEquilibrium` (for example
  `phase_equilibrium_problem.rs`). The phase subsystem supplies data and typed
  layout; it must not know solver policies, residuals, or backend types.
  - [x] The foundational `phase_equilibrium_problem` module now owns the typed
    request and immutable structural metadata.
  - [x] Add standard-state extraction and final `EquilibriumProblem`
    construction without leaking solver policy into `phase_*`.
- [x] Define a single request containing `&ResolvedPhaseSystem`, fixed
  `EquilibriumConditions`, typed initial composition, trace-seed policy, and
  supported-model policy.
- [x] Return a typed bridge bundle containing the canonical
  `EquilibriumProblem`, `SystemLayout`, component/phase index maps, immutable
  lookup provenance, and a build report. Do not return parallel unrelated
  vectors that callers must keep aligned manually.
  - [x] Introduce immutable `PhaseEquilibriumMetadata` with the canonical
    layout, fingerprint, component/phase descriptors, index maps, and retained
    `ResolvedPhaseSystemReport`.
  - [x] Add the numerical `EquilibriumProblem` and standard-state build report
    only after all thermochemical extraction has succeeded.
- [x] Build the element-composition matrix and element labels in exact
  `SystemLayout` component order. Validate that the same substance resolved in
  two phases has identical molecular composition even when its thermochemical
  records differ.
- [x] Extract one standard-state `G0(T)` function/expression per qualified
  component from its own phase-local `SubsData` record. Do not pass the
  composition-corrected Gibbs value from `evaluate_gibbs`: the canonical
  equilibrium residual already applies mixing and pressure activity terms.
- [x] Evaluate every standard-state model once at the requested temperature
  during preparation and reject missing, out-of-range, NaN, or infinite data
  with the offending `PhaseComponentId` and source provenance.
- [x] Map `IdealGas` and one-component `PureCondensed` to canonical activity
  models; return a typed unsupported-model error for every other combination.
- [x] Preserve `ResolvedPhaseSystemReport` provenance by component, including
  library, record key, physical-state match, and NIST fallback evidence.
- [x] Add a stable preview/report surface showing conditions, ordered
  components, phase models, initial amounts, element totals, thermochemical
  sources, and solver-facing labels without running a solver.

Required tests:

- [x] Adapter output agrees exactly with direct `SubsData` standard-state Gibbs
  and elemental-composition calculations for each component.
- [x] A regression proves that ideal-gas mixing/pressure terms are applied once,
  not once in `SubsData` and again in the equilibrium residual.
- [x] Missing data in the last phase leaves no partially built bridge result.
- [x] Local NASA gas plus NASA condensed records build fully offline and keep
  per-component provenance.

### P4.4 Solve a fixed active phase set through the canonical engine

- [x] Make the pure numeric and symbolic formulations consume the same typed
  phase/activity descriptors; remove `PhaseKind` switches duplicated across
  residual, Jacobian, symbolic generation, validation, and K_eq support.
  `Phase` now stores the canonical `PhaseActivityModel` directly. The numeric
  residual, RST symbolic generator, phase-stability logic, and K_eq guard read
  that same value. `PhaseKind` is only a temporary type alias for older test
  fixture literals, not a second runtime representation.
- [x] Implement and test the chemical-potential/activity term for each
  supported model:
  - ideal gas: `ln(x_i * P / P0)`;
  - one-component pure condensed phase: zero activity correction.
  - `equilibrium_activity` now owns the common `ln(a_i)` contract and the
    pressure offset used by both numeric and symbolic residual construction.
- [x] Build an immutable active-set projection that maps global
  `PhaseComponentId` values to local nonlinear coordinates and can scatter an
  accepted local solution back to the full `SystemLayout`.
  - The ChemEquilibrium-local `ActiveSetProjection` currently uses typed
    `SpeciesId`/`PhaseId`; upgrading its global identity to `PhaseComponentId`
    belongs to the later phase-subsystem bridge.
- [x] Solve the all-declared-phases-active case first through the existing
  backend cascade and backend-independent acceptance gate.
  `PhaseEquilibriumProblemBundle::solve_with` now consumes the prepared
  bridge bundle, configures only `EquilibriumSolverSettings`, and returns one
  immutable solution bundle with the accepted snapshot, build provenance, and
  backend trace. The offline NASA-gas contract covers this end-to-end path.
- [x] Extend the RST symbolic thermochemistry snapshot to build from the typed
  bridge bundle rather than from the legacy single `solver.subs_data` field.
  The bridge now owns phase-qualified symbolic `G0(T)` expressions alongside
  numeric closures. It injects them into the canonical solver, and the RST
  adapter prefers that immutable ordered snapshot before falling back to the
  historical mutable `SubsData` facade. The local NASA-gas bridge regression
  asserts that the accepted backend is RustedSciThe.
- [x] Ensure legacy nonlinear backends can consume the same prepared active-set
  problem as an explicit fallback, without creating a second multiphase
  formulation.
- [x] Keep fixed-active-set solve publication transactional: only an accepted,
  globally re-expanded solution and its complete reports become visible.
  The non-phase-controlled bridge solve is also consuming and transactional:
  invalid backend settings return no accepted result and leave the resolved
  phase data untouched.

Required tests:

- [x] Numeric residual, analytic Jacobian, and RST symbolic residual agree for
  ideal-gas plus pure-condensed fixtures.
  - The mixed fixture runs at `P != P0` and also caught a corrected Jacobian
    defect: phase-total derivatives apply to non-reacting species in the same
    phase as a reaction participant.
- [x] Every configured RST backend sees identical equations and component
  order; fallback attempts preserve the same active-set projection.
  The bridge matrix configures every RST backend independently against one
  local NASA-gas bundle. Each reaches the same bridge-owned symbolic contract;
  methods that do not accept the fixture report a normal one-attempt backend
  failure rather than a setup error. A Nielsen -> LM cascade preserves the
  exact qualified component order and accepts the second attempt.
- [x] A failed backend cascade cannot overwrite a previous accepted
  multiphase result.
  Bridge solves publish only immutable `PhaseEquilibriumSolutionBundle` values.
  A regression retains an accepted NASA-gas snapshot, forces a later Nielsen
  failure, and verifies both accepted moles and retained source provenance are
  unchanged.

### P4.5 Replace legacy phase control with a bounded active-set algorithm

This subsection is about the current helper-level API surface in
`equilibrium_workflows.rs`: `solve_with_phase_control`,
`compute_phase_creation_dg`, `activate_phase`,
`deactivate_phases`, and `deactivate_phases_mass_conserving`.
The complaint is not the phase subsystem interface itself; it is the behavior
and contract of these old helpers.

- [x] Preserve the updated log-mole seed before every outer-loop `continue`
  so the next nonlinear attempt starts from the new phase-adjusted state.
- [x] Add a hard outer-loop iteration cap and a typed
  `PhaseControlDidNotConverge { iterations }` error for phase-control
  non-convergence.
- [x] Retire the current `deactivate_phases_mass_conserving` behavior from
  production paths. Phase removal should seed the deactivated phase to trace
  floor only, not transfer its moles to an arbitrary receiver species.
- [x] Replace the current phase-creation criterion with an explicit phase
  stability contract. Production should support only the physically justified
  one-component condensed-phase case first; multicomponent ideal-solution
  stability should return `ValidationNotApplicable` until tangent-plane
  minimization exists.
- [x] Make `activate_phase` a seed-only helper: mark the phase active, seed it
  with a small positive mole number, and let the next full nonlinear solve
  restore the element balance.
- [x] Keep the active-phase set as persistent state across restarts instead of
  recreating it opportunistically inside one solve call.
- [x] Encode hysteresis explicitly with separate create/keep thresholds and a
  typed transition plan so phase flip-flopping cannot consume the outer loop.
- [x] Do not promote the current `PhaseManager`, `compute_phase_creation_dg`,
  or in-place `activate_phase`/`deactivate_phases` helpers to the canonical
  path. Characterize any useful behavior, then replace their unbounded mutable
  loop with typed orchestration.
  - The raw creation-score and compatibility activation/deactivation helpers
    have been removed. The canonical boundary now exposes only physical
    stability reports and typed seed-only transitions.
- [x] Introduce `PhaseSet`/`PhaseStatus` for declared, active, inactive,
  excluded, appeared, and disappeared phases. Every transition must retain the
  semantic `PhaseId`.
  - `InitialPhaseSet` supports inventory-derived, all-candidate, and explicit
    active/excluded policies. Excluded phases do not enter stability
    classification, while the nonlinear backend sees only the immutable active
    mask for one fixed-set attempt.
- [x] Solve only active components, then compute a typed phase-stability report
  for every inactive candidate. For the first supported slice, implement the
  mathematically justified stability criterion for one-component pure
  condensed phases against the accepted gas/condensed state.
  - Each outer iteration now builds a reduced species/phase projection,
    recomputes its reaction basis, solves it against the unchanged original
    element totals, and expands only the accepted candidate back into declared
    ordering.
- [x] Derive any elemental potentials/reduced chemical potentials required by
  the stability test from the accepted state with a checked linear solve and a
  reported residual; do not infer phase stability from raw `sum(G0)`.
- [x] Remove an active phase only when its amount is below the destruction
  threshold and the resulting inactive phase satisfies the stability
  inequality. Never transfer its moles to an arbitrary first species.
- [x] Add one most-unstable inactive phase per restart, apply only a numerical
  trace seed, and recompute the reduced active problem against the original
  physical element totals. The seed is not claimed to conserve elements; every
  solved candidate must pass the common element-balance gate before it can
  trigger another transition or be published.
- [x] Use separate create/keep thresholds (hysteresis), a maximum restart
  count, and visited-phase-set cycle detection. Report oscillation and budget
  exhaustion as typed non-convergence errors.
- [x] Stage the whole outer solve locally. If any restart, stability test, or
  backend cascade fails, retain the previously published solution unchanged.
  - `solve_candidate_from_seed` now returns a validated unpublished bundle;
    `solve_with_phase_control` publishes it only after the transition
    classifier reaches a fixed point.
- [x] Accept the final result only when the nonlinear candidate gate passes,
  every active phase is internally valid, and every inactive candidate meets
  the phase-stability tolerance.
  - The outer loop refuses to publish material in an inactive phase, rejects
    unsupported multicomponent solution models, rejects missing/rank-deficient
    elemental-potential references, and checks the elemental-potential fit
    residual before using a phase driving force.

Required tests:

- [x] A stable absent condensed phase stays absent.
- [x] An unstable absent condensed phase appears and converges after a bounded
  restart.
- [x] A vanishing unstable/stable boundary fixture exercises hysteresis without
  cycling.
- [x] Deliberately alternating phase sets terminate with a typed cycle/budget
  report instead of looping forever.
- [x] A failure after an accepted fixed-set candidate but before outer-loop
  convergence leaves the previously published solution and reports unchanged.
- [x] Every intermediate and final state conserves each element within the
  common acceptance tolerance.
  - Transition records retain the candidate validation evidence that admitted
    each fixed-set solve; trace seeds are restart coordinates, not published
    thermodynamic states.

### P4.5a Production hardening and regression evidence

The bounded active-set architecture is now canonical, but production readiness
requires proof that every reduced solve preserves the original closed-system
inventory and that near-boundary numerical behavior cannot manufacture phase
transitions.

- [x] Assert in `solve_fixed_active_set_candidate` and its regression tests
  that element totals always originate from the full physical problem, never
  from a reduced seed or a trace-seeded inactive component vector.
  The local solver now receives totals recomputed from full `n0` and the full
  element matrix; the regression compares the accepted phase-control moles
  with that original inventory.
- [x] Validate active-set feasibility before invoking a nonlinear backend:
  the original element-total vector must belong to the column space of the
  active species element matrix. Reject an impossible active set with typed
  `InvalidProblem { field: "phase_active_set", .. }`.
  `ActiveSetProjection` now checks the reduced SVD feasibility residual before
  any backend call. Its pre-existing square/rank gate is now classified as the
  same typed invalid active-set error instead of `ValidationNotApplicable`.
- [x] Document and test the strict reduced-system contract: the active element
  matrix must retain every conserved-element direction required by the
  physical inventory, and the reduced reaction/elements equation count must
  be square in the present log-moles formulation. Explain why a phase set
  without a carbon-bearing component cannot solve a carbon-containing system.
  The projection documentation and its impossible H/C active-set regression
  now make this boundary explicit.
- [x] Define the inactive ideal-gas contract. Until a true gas-phase split or
  tangent-plane criterion exists, trace-floor gas species must not generate
  an artificial stability driving force or activate a second gas phase.
  `compute_phase_stability_reports` classifies every gas phase as
  `FixedIdealGas` with no driving force; a regression uses an extreme Gibbs
  value on an inactive trace gas and proves it cannot enter the transition
  plan.
- [x] Complete the activity-law audit with a source-level regression: numeric
  residual, analytic Jacobian, symbolic residual, and phase-stability code
  must delegate to `equilibrium_activity`; no duplicate `ln(a_i)` formula may
  be added outside that module.
  The deprecated residual path now also routes phase offsets through the same
  helper, and a regression test locks the canonical and deprecated residual
  behavior together.
 - [x] Expand `ActiveSetProjection` tests for non-consecutive global indices,
   multiple active phases, local phase remapping, scatter/project round trips,
   preserved ordering, inactive trace floors, reduced basis dimensions,
   dependent element columns, and impossible active-set rejection.
  - [x] Sparse non-consecutive active species now round-trip through
    `project_log_moles` and `scatter_log_moles` while preserving global
    ordering and inactive floors.
  - [x] Projected local values now round-trip through `scatter_log_moles` and
    `project_log_moles` without reordering.
  - [x] The reduced basis reports the expected reaction dimension for a
    nontrivial active set, so the projection contract is pinned instead of
    inferred.
  - [x] Impossible sparse active sets now fail before backend publication
    instead of being silently coerced into a reducible phase inventory.
- [x] Expand phase-control stories: appearance, disappearance, appearance
  followed by disappearance, disappearance followed by reappearance, several
  sequential transitions, maximum-restart termination, A -> B -> A cycle,
  and transactional preservation after a failed intermediate restart.
  - [x] A first transactional story now covers appearance followed by a
    failed restaging attempt, and proves that the previously published
    solution and moles remain intact.
  - [x] The intermediate publication also survives a second solve attempt
    that fails before changing the published active set or mole state.
  - [x] Maximum-restart termination is covered by a typed story that confirms
    the solver stops cleanly after one transition budget is exhausted.
  - [x] The activation/deactivation seed helpers now cover a full
    appearance -> disappearance -> reappearance cycle without relying on a
    brittle full-solve stress case.
  - [x] Several sequential transition scenarios are now covered through
    fresh fixed-fixture solves so the publication contract remains stable
    across appearance, disappearance, and reappearance cases.
- [x] Add deterministic numerical stress fixtures covering trace species,
  large Gibbs-energy scales, nearly rank-deficient reaction bases, several
  pure condensed candidates, and low/high temperatures. Assert finite values,
  accepted element conservation, and typed failure rather than NaN/Inf.
  - [x] A first stress fixture now exercises trace species against extreme
    Gibbs-energy scales across low, nominal, and high temperatures, while
    checking finite accepted moles and exact element conservation.
  - [x] A second stress fixture now covers multiple pure condensed
    candidates and confirms that the phase-control solve remains finite and
    element-conserving.
- [x] Run at least LM and Trust-Region policies through the same phase-control
  fixtures and compare accepted phase sets, mole vectors, and validation
  reports within documented tolerances.
  - [x] LM and TR now match on both a pure-condensed appearance fixture and
    a stable gas-only fixture, with checked phase sets, mole vectors, and
    element-balance tolerances.
- [x] Make hysteresis temperature-aware. Store dimensionless create/keep
  coefficients or derive thresholds from the solve temperature immediately
  before phase control; remove the current accidental 298 K default from the
  physical policy contract. Phase control now resolves the hysteresis band
  from the current solve temperature instead of baking in a fixed 298 K
  constant, and the policy itself can be expressed either explicitly or as
  dimensionless RT factors.
- [x] Add ignored, opt-in performance characterizations for 50--100 species
  and several phases. Record outer iterations, nonlinear iterations,
  transitions, and elapsed time, but do not impose flaky wall-clock limits in
  the default test suite.
  - [x] A first ignored characterization scaffold now runs a larger
    multi-phase synthetic system and prints iterations, transitions, and
    elapsed time without affecting the default suite.
- [x] Establish an always-offline long-term regression matrix for O2 <-> 2O,
  N2 <-> 2N, N/O mixtures, condensed appearance/disappearance, trace species,
  and deterministic generated small systems. Each fixture must state a
  physical invariant, not merely retain a historical iterate.
  - [x] A first offline matrix now covers O2/O dissociation, N2/N
    dissociation, a diluted N/O gas mixture, and synthetic condensed
    appearance/disappearance, with finite accepted states and bounded
    element-balance drift.
  - [x] A second offline matrix now covers additional deterministic small gas
    systems, including NO/N2/O2 and diluted O2/O/N2 cases, without network
    access.

### P4.6 Publish a phase-aware immutable solution and reports

- [ ] Add `MultiphaseEquilibriumSolution` containing fixed conditions, layout
  revision/fingerprint, ordered component moles, phase totals, phase-local mole
  fractions, active/inactive status, and the accepted canonical solution.
  - [x] The fixed-active bridge result is now published as an immutable
    `MultiphaseEquilibriumSolution` with fingerprint checks, phase totals,
    local mole fractions, qualified lookups, provenance, and backend evidence.
  - [ ] Wire accepted `PhaseControlledSolveReport`/`PhaseSet` states into this
    same result once the bounded outer loop consumes the bridge rather than the
    historical mutable solver.
- [x] Provide lookups by `PhaseComponentId` and `PhaseId`; expose aggregate
  totals by bare substance only as an explicit derived view.
- [x] Retain build/lookup provenance and the complete nested backend solve
  report in the result bundle.
- [x] Add `PhaseTransitionReport` entries with previous/new phase sets, reason,
  seed, stability metric, backend outcome, and elemental-balance evidence.
  - The immutable `PhaseControlledSolveReport` now retains typed initial/final
    phase sets, ordered transitions, explicit physical reasons, exact restart
    seeds, phase totals, driving forces, per-candidate validation evidence, the
    final validation report, and every nonlinear backend report. The remaining
    report work is integration into the final multiphase solution bundle.
- [ ] Add a final `MultiphaseAcceptanceReport` combining canonical residual and
  element checks with phase stability/complementarity checks.
  - [x] The final bundle now exists and combines phase-control evidence,
    canonical validation, phase-stability reports, and a complementarity
    summary with stable rows and `Display` output.
- [x] Provide stable summary rows and `Display` output for CLI, future GUI, and
  snapshot tests; core code must not print directly.
  - [x] `PhaseControlledSolveReport` now exposes stable summary rows and a
    `Display` implementation for CLI and snapshot-friendly output.
- [ ] Reject result reconstruction when the result layout does not match the
  resolved-system fingerprint.
  - [x] `MultiphaseInitialComposition` already refuses reconstruction against
    a foreign layout fingerprint, and the regression tests now pin that
    contract down explicitly.

### P4.7 Complete the multiphase test matrix

- [ ] Add a dedicated `equilibrium_phase_bridge_tests.rs` for identity,
  ordering, adapter validation, provenance, and transactional failures.
- [x] Add a dedicated `equilibrium_multiphase_story_tests.rs` with module-level
  documentation describing each physical hypothesis and expected result.
- [ ] Cover at minimum:
  - ideal-gas-only parity with the current canonical solver;
  - gas plus one stable pure solid;
  - gas plus one stable pure liquid;
  - the same molecule represented in gas and condensed phases;
  - two independent pure condensed candidate phases;
  - unsupported multi-component condensed solution;
  - phase appearance, disappearance, hysteresis, and cycle detection;
  - offline mixed NASA gas/NASA condensed lookup provenance.
- [ ] For every applicable small fixed-phase set, compare the accepted result
  with the independent equilibrium-constant solver. Report non-applicability
  explicitly for phase-appearance decisions that the K_eq problem does not
  model.
- [ ] Add exact component-order and numeric/closure/symbolic equivalence tests,
  then run the existing RST backend/fallback matrix over at least one physical
  multiphase fixture.
- [ ] Adopt physically meaningful fixtures from the retired classical stack
  only after restating their expected invariants. Do not preserve tests whose
  only contract is legacy mutable state or one historical numeric iterate.
- [ ] Keep all default tests offline and leave thermochemical library files
  byte-for-byte unchanged.

### P4.8 Migrate workflows and retire the duplicate engine

- [ ] Add one public one-shot facade such as `solve_resolved_pt` accepting the
  resolved phase system, typed initial composition, conditions, solver policy,
  and phase-control policy.
  - [x] `ResolvedPhaseEquilibriumRequest` and `solve_resolved_pt` now own the
    fixed-declared-phase bridge transaction, including typed numerical
    settings and immutable result publication.
  - [ ] Extend the same request with the bridge-backed bounded phase-control
    policy only after phase-control no longer relies on the mutable legacy
    workflow. Do not claim the existing helper is a production-equivalent
    public mode.
- [ ] Migrate useful legacy equilibrium workflows and examples one
  physical scenario at a time onto that facade.
- [ ] During migration, compare both implementations only in characterization
  tests; do not expose two production APIs as equivalent long-term choices.
- [ ] Migrate equilibrium consumers away from the broad
  `ThermodynamicsCalculatorTrait`, raw Gibbs closures, and nested legacy phase
  maps to the narrow typed bridge.
- [ ] Deprecate `solve_with_phase_control` and the old mutable phase helpers
  after the active-set story matrix passes.
- [ ] Deprecate and then delete the remaining legacy equilibrium solving
  after all retained physical scenarios and public examples have migrated.
- [ ] Update the phase-subsystem architecture documents with the final
  dependency direction and add a fixed-P,T multiphase usage example.

### Recommended P4 implementation passes

1. [ ] Component/phase identity refactor in `EquilibriumProblem` plus tests.
2. [ ] Typed initial composition and pure `ResolvedPhaseSystem` adapter.
3. [ ] Standard-state thermochemistry/provenance bridge and preview report.
4. [ ] Fixed-active-set numeric, symbolic, RST, and legacy-fallback parity.
5. [ ] Bounded phase-stability active-set orchestration.
6. [ ] Phase-aware solution/report boundary and full acceptance gate.
7. [ ] Offline physical story matrix and independent K_eq cross-validation.
8. [ ] Classical workflow migration, deprecation, and duplicate-engine removal.

### P4 definition of done

- [ ] One typed request can solve a resolved closed multiphase system at fixed
  `P,T` without flattening phase-qualified identity.
- [ ] Standard-state thermochemistry is evaluated once per component and every
  record retains lookup provenance.
- [ ] All supported activity terms are explicit and tested; unsupported phase
  models fail before a backend starts.
- [ ] Phase appearance/disappearance is bounded, transactional, conservative,
  and included in the acceptance report.
- [ ] The result can be queried unambiguously by phase and component and is
  tied to the exact resolved layout that produced it.
- [ ] Numeric, analytic-Jacobian, and symbolic formulations agree; the RST
  backend matrix and explicit legacy fallback pass the same physical fixture.
- [ ] Offline gas/condensed stories pass elemental, residual, stability, and
  applicable independent K_eq validation gates.
- [ ] No retained production workflow requires the duplicate classical
  equilibrium engine.

### Post-P4 optimization and review notes

These items are useful, but they are not blocking the core correctness work
above. Keep them in mind once the phase-control and multiphase acceptance
paths are stable.

- [x] Verify the phase-driving-force contract directly from chemical
  potentials, and add one regression that makes the formula visible at the
  phase-boundary level rather than only through downstream transition plans.
  The new regression checks a pure-condensed candidate against an active gas
  reference and pins the `g0 + RT ln(a)` activity contribution in the
  reported driving force.
- [x] Add a dedicated hysteresis story that exercises
  `appear -> keep -> disappear -> reappear` behavior explicitly and pins the
  transition thresholds against accidental regression.
  The regression now exercises the public phase-classification path with an
  explicit hysteresis band and checks appearance, keep, disappearance, and
  reappearance decisions in order.
- [ ] Review whether `ActiveSetProjection` and the surrounding fixed-layout
  objects should cache reusable projections for the `A -> B -> A` case once
  correctness is fully stable.
- [x] Add a separate benchmark suite for larger synthetic systems
  (20, 50, 100, 200 species) that reports outer iterations, nonlinear
  iterations, projection build time, Jacobian time, and total solve time.
  The ignored benchmark sweep now covers 20, 50, 100, and 200 species
  synthetic inventories and prints projection-build and solve timing together
  with iteration counts.
- [ ] Revisit allocation-heavy paths only after the benchmark data shows a
  meaningful regression; do not prematurely optimize them before the
  correctness and acceptance story is complete.
- [ ] After the phase* bridge is connected, add live resolved-phase-system
  regression fixtures and gradually migrate selected synthetic stories to
  those real inputs. Keep the synthetic cases until the live coverage reaches
  the same physical invariants and diagnostic strength.

## P5 - Public API, GUI, and cleanup

- [ ] Provide a small public facade: validated problem builder, solver policy,
  solve method, solution, solve report, and optional validation report.
- [ ] Expose backend selection and cascade diagnostics in GUI only through typed
  controls; never require users to type internal enum names.
- [ ] Show conservation, residual, fallback-attempt, and K_eq validation status
  as first-class result sections.
- [ ] Add GUI story tests for success, fallback success, all-backends-failed,
  invalid input, validation mismatch, and save/load roundtrip.
- [ ] Split the implementation by responsibility: `domain`, `formulation`,
  `validation`, `backend`, `solver_policy`, `keq_validation`, and `report`.
- [ ] Review `easy_equilibrium.rs` as the future facade or replace it cleanly.
- [x] Rename the canonical modules and their tests so public paths no longer
  imply that the main formulation is only an equilibrium-constant solver.
- [ ] Remove or convert leftover debug modules, stale commented code,
  direct `println!`, and duplicated temperature-sweep implementations.
  `Untitled-1.rs` has been deleted; keep this item focused on the remaining
  experimental and legacy cleanup.

## Recommended implementation passes

1. [ ] Characterize current fixtures and introduce typed problem/result/error
   boundaries without changing numerical behavior.
2. [ ] Extract the pure log-moles formulation and common acceptance gate.
3. [ ] Characterize the RST backend matrix, then complete attempt counters,
   termination mapping, and a measured production default.
4. [ ] Isolate legacy solvers and keep them behind an explicit fallback
   policy with typed diagnostics.
5. [ ] Build the independent K_eq validator and cross-validation fixtures.
6. [ ] Bridge `ResolvedPhaseSystem`, migrate retained classical workflows, and
   remove the duplicate equilibrium engine.
7. [ ] Build the typed GUI on the stable facade and reports.

## Definition of done

- [ ] There is one canonical equilibrium problem and solution model.
- [ ] Every accepted solution has passed backend-independent numerical and
  physical validation.
- [ ] Fallback order, retry rules, budgets, and attempts are explicit and fully
  reported.
- [ ] RustedSciThe provides the production nonlinear backends.
- [ ] The K_eq path independently validates every applicable small-system
  fixture and clearly reports when it is not applicable.
- [ ] Tests cover formulation, every backend, cascade behavior, physical
  invariants, cross-validation, transactionality, and offline integration.
- [ ] Legacy solvers and the duplicate classical equilibrium engine are removed
  after their migration gates pass.

---

## SourceCraft Diagnostics

Результаты ревизии кода модуля `ChemEquilibrium`, проведённой 23.07.2026.

### A. Потенциальные ошибки (bugs)

#### A.1 `R = 8.314` — неверное значение газовой постоянной

**Где:** [`equilibrium_log_moles.rs:100`](equilibrium_log_moles.rs:100)

**Проблема:** `pub const R: f64 = 8.314;` — это значение не является стандартной молярной газовой постоянной. Правильное значение: `8.314462618`. Ошибка в 0.00046 J/(mol·K) даёт систематическое смещение в логарифме константы равновесия: `ln(K) = -ΔG/(RT)`. При T=1000K, ΔG=-100 kJ/mol: `ln(K)_true = 100000/(8.31446*1000) = 12.027`, `ln(K)_current = 100000/(8.314*1000) = 12.028`. Ошибка мала (~0.01%), но накапливается в температурных сериях.

**Важность:** Средняя. Для инженерных расчётов ошибка незначительна, но для научных публикаций — недопустима.

#### A.2 `has_rst_symbolic_context()` — некорректная эвристика

**Где:** [`equilibrium_log_moles.rs:490-496`](equilibrium_log_moles.rs:490)

**Проблема:** Метод проверяет наличие символьного контекста через `gibbs_sym.len() == substances.len() && !gibbs_sym.is_empty()`, но также возвращает `true`, если `search_results`, `search_states`, `therm_map_of_sym` или `therm_map_of_fun` непусты. Последние четыре условия не гарантируют, что символьные Gibbs-функции действительно построены для всех веществ. Это может привести к тому, что `SolverPolicy::rusted_scithe_default()` будет выбран, но RST-адаптер не сможет построить SymbolicNonlinearProblem.

**Важность:** Средняя. Может проявляться как трудноотлавливаемая ошибка "RST backend failed" при определённых последовательностях вызовов.

#### A.3 `compute_species_moles` — публичная функция с побочным эффектом

**Где:** [`equilibrium_log_moles.rs:2173-2178`](equilibrium_log_moles.rs:2173)

**Проблема:** Публичный метод `compute_species_moles(&mut self, sol: Vec<f64>)` мутирует `self.moles` и `self.map_of_moles_for_each_substance`. Название предполагает чистую функцию (как у свободной функции `compute_species_moles(log_moles: &[f64])`), но на самом деле это метод с побочным эффектом. Это нарушает принцип наименьшего удивления.

**Важность:** Низкая. Косметическая проблема именования.

#### A.4 `EquilibriumLogMoles` — 22 публичных поля

**Где:** [`equilibrium_log_moles.rs:275-320`](equilibrium_log_moles.rs:275)

**Проблема:** Структура имеет 22 публичных поля, которые можно изменять извне в любом порядке. Нет гарантии, что после изменения одного поля (например, `T`) остальные поля (например, `gibbs`) остаются консистентными. Это прямой путь к багам "забыл обновить gibbs после смены T".

**Важность:** Высокая. Основной источник производственных ошибок.

### B. Overengineering

#### B.1 `SolverCascadeBudget` — избыточная сложность для текущих потребностей

**Где:** [`equilibrium_solver_policy.rs:38-46`](equilibrium_solver_policy.rs:38)

**Проблема:** Бюджет каскада содержит три независимых лимита (`max_attempts`, `max_iterations_per_attempt`, `max_total_iterations`), которые проверяются в [`solve_backend_cascade`](equilibrium_log_moles.rs:2064-2073). При этом `max_total_iterations` вычисляется как `max_iter * backends.len()` — то есть всегда пропорционален двум другим. Три лимита вместо одного — это overengineering. Достаточно двух: `max_attempts` и `max_iterations_per_attempt`.

**Важность:** Средняя. Усложняет понимание без реальной выгоды.

#### B.2 `VariableScalingContract` — объявлен, но нигде не используется

**Где:** [`equilibrium_problem.rs:313-384`](equilibrium_problem.rs:313)

**Проблема:** `VariableScalingContract` полностью реализован (с `apply_iterate`, `unscale_iterate`, валидацией), но не используется ни в одном бэкенде. Row scaling (`ResidualScalingContract`) используется, а variable scaling — нет. Это мёртвый код.

**Важность:** Средняя. Увеличивает когнитивную нагрузку при чтении.

#### B.3 `EquilibriumConstantSolverMode::Required` — нереализуемый контракт

**Где:** [`equilibrium_constant_solver.rs:74-75`](equilibrium_constant_solver.rs:74)

**Проблема:** Режим `Required` требует, чтобы K_eq валидация была выполнена для каждого решения. Но K_eq решатель работает только для однореакционных систем. Для многореакционных систем `Required` гарантированно упадёт с ошибкой. Этот режим нельзя использовать в production, только в тестах.

**Важность:** Низкая. Но вводит в заблуждение.

#### B.4 `TemperaturePostprocessingPolicy` — избыточная гибкость

**Где:** [`equilibrium_temperature_postprocessing.rs`](equilibrium_temperature_postprocessing.rs)

**Проблема:** Политика постобработки поддерживает три типа сеток (RawOnly, Uniform, Explicit) и два типа интерполяции (Linear, Log). При этом в коде нет ни одного вызова с Explicit или Log. Вся гибкость существует только для тестов.

**Важность:** Низкая. YAGNI-нарушение.

#### B.5 `PhaseHysteresisPolicy::Explicit` — дублирование TemperatureScaled

**Где:** [`equilibrium_workflows.rs:1017-1025`](equilibrium_workflows.rs:1017)

**Проблема:** `Explicit { dg_create, dg_keep }` — это фактически `TemperatureScaled` с замороженной температурой. При T=const они эквивалентны. Можно было бы обойтись одним вариантом.

**Важность:** Низкая.

### C. Недостатки дизайна

#### C.1 `EquilibriumLogMoles` — God Object

**Где:** [`equilibrium_log_moles.rs:275-320`](equilibrium_log_moles.rs:275)

**Проблема:** Структура содержит 22 поля, объединяющие термохимические данные, настройки решателя, состояние температурной серии, опубликованные решения, отчёты валидации и фазы. Это классический God Object. Новый API (`EquilibriumProblem` + `EquilibriumSolution`) уже решает эту проблему, но старый фасад остаётся в production.

**Важность:** Высокая. Затрудняет тестирование и поддержку.

#### C.2 Две параллельные иерархии ошибок

**Где:** [`equilibrium_nonlinear.rs:55-72`](equilibrium_nonlinear.rs:55) и [`equilibrium_nonlinear.rs:89-188`](equilibrium_nonlinear.rs:89)

**Проблема:** `SolveError` (для численных решателей) и `ReactionExtentError` (для всего остального) — две пересекающиеся иерархии. `SolveError` мог бы быть вариантом `ReactionExtentError`, но они разделены. Это приводит к тому, что в некоторых местах ошибка оборачивается (`ReactionExtentError::SolveError`), а в некоторых — нет.

**Важность:** Средняя.

#### C.3 `from_problem` — деструктуризация типобезопасности

**Где:** [`equilibrium_log_moles.rs:615-644`](equilibrium_log_moles.rs:615)

**Проблема:** Метод `from_problem` принимает типобезопасный `EquilibriumProblem`, но раскладывает его в 22 публичных поля `EquilibriumLogMoles`. Вся типобезопасность теряется. Это мост для миграции, но он должен быть временным.

**Важность:** Средняя.

#### C.4 `solve_fixed_active_set_candidate` — клонирование всего solver'а

**Где:** [`equilibrium_workflows.rs:1439-1476`](equilibrium_workflows.rs:1439)

**Проблема:** Для решения с проекцией активного набора создаётся полная копия `EquilibriumLogMoles` (`let mut local = EquilibriumLogMoles::empty()`), в которую копируются поля по одному. Это 30+ строк ручного копирования. Любое новое поле в `EquilibriumLogMoles` нужно не забыть добавить и сюда.

**Важность:** Средняя. Хрупкий код.

### D. Пробелы в тестовом покрытии — статус

> **Статус:** 69 из 70 пунктов закрыты. Остаётся D.2 (`solve_for_T_range` + phase control).

#### D.1 ✅ `EquilibriumLogMoles::solve()` с Legacy бэкендами

**Статус:** 5 тестов в `equilibrium_log_moles_tests.rs` (`legacy_lm_solve_publishes_accepted_solution`, `legacy_nr_solve_publishes_accepted_solution`, `legacy_tr_solve_publishes_accepted_solution`, `legacy_solve_fails_gracefully_without_stoich_matrix`, `legacy_solve_publishes_solve_report_with_attempts`).

#### D.2 ❌ `solve_for_T_range` с phase control

**Проблема:** Температурные серии и управление фазами не тестируются вместе. Нет теста, который бы запустил `solve_with_phase_control()` внутри `solve_for_T_range()`.

#### D.3 ✅ `EquilibriumLogMoles::from_problem()`

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`from_problem_preserves_all_problem_data`, `from_problem_rejects_invalid_problem`).

#### D.4 ✅ `accepted_solution()` после `solve_with_phase_control()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`accepted_solution_after_phase_control_returns_ok_for_single_phase`, `accepted_solution_after_phase_control_preserves_element_totals`, `accepted_solution_after_phase_control_fails_without_solve`).

#### D.5 ✅ `current_problem_snapshot()`

**Статус:** Косвенно через `from_problem` тесты (D.3). Прямой вызов невозможен — метод приватный.

#### D.6 ✅ `run_keq_cross_validation()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`solve_with_keq_validation_enabled_does_not_crash`).

#### D.7 ✅ `EquilibriumSolverSettings::validate()`

**Статус:** 11 тестов в `equilibrium_log_moles_tests.rs` (`solver_settings_validate_accepts_default_settings`, `solver_settings_validate_rejects_zero_max_iter`, `solver_settings_validate_rejects_non_finite_tol`, `solver_settings_validate_rejects_zero_tol`, `solver_settings_validate_rejects_negative_lambda`, `solver_settings_validate_rejects_non_finite_alpha_min`, `solver_settings_validate_rejects_delta_max_less_than_delta_init`, `solver_settings_validate_rejects_eta_out_of_range`, `solver_settings_validate_accepts_single_backend_policy`, `solver_settings_validate_rejects_zero_budget_limits`, `solver_settings_validate_rejects_non_finite_keq_tolerances`).

#### D.8 ✅ `compute_phase_totals()`

**Статус:** 3 теста в `equilibrium_log_moles_tests.rs` (`compute_phase_totals_sums_moles_by_phase`, `compute_phase_totals_handles_empty_input`, `compute_phase_totals_handles_single_phase`).

#### D.9 ✅ `initial_phase_activity()`

**Статус:** 5 тестов в `equilibrium_log_moles_tests.rs` (`initial_phase_activity_marks_active_phases`, `initial_phase_activity_rejects_dimension_mismatch`, `initial_phase_activity_rejects_invalid_phase_eps`, `initial_phase_activity_rejects_out_of_bounds_phase`, `initial_phase_activity_rejects_non_finite_moles`).

#### D.10 ✅ `build_multiphase_acceptance_report()`

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`build_multiphase_acceptance_report_rejects_dimension_mismatch`, `build_multiphase_acceptance_report_accepts_matching_dimensions`).

#### D.11 ✅ `PhaseManager::detect_phase_destruction()`

**Статус:** 4 теста в `equilibrium_log_moles_tests2.rs` (`detect_phase_destruction_returns_empty_when_no_phase_below_threshold`, `detect_phase_destruction_ignores_inactive_phases_below_threshold`, `detect_phase_destruction_finds_active_phase_below_threshold`, `detect_phase_destruction_handles_empty_input`).

#### D.12 ✅ `PhaseManager::classify_phases()` (без температуры)

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`classify_phases_returns_no_transition_when_all_phases_stable`, `classify_phases_rejects_temperature_scaled_hysteresis`, `classify_phases_rejects_dimension_mismatch`).

#### D.13 ✅ `reject_repeated_phase_set()` с пустым множеством

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`reject_repeated_phase_set_accepts_first_occurrence`, `reject_repeated_phase_set_rejects_duplicate`).

#### D.14 ✅ `validate_phase_set_candidate()` с корректным набором

**Статус:** 4 теста в `equilibrium_log_moles_tests.rs` (`validate_phase_set_candidate_accepts_valid_set`, `validate_phase_set_candidate_rejects_dimension_mismatch`, `validate_phase_set_candidate_rejects_negative_total`, `validate_phase_set_candidate_rejects_inactive_phase_with_moles`).

#### D.15 ✅ `seed_activated_phase()` с разными `PhaseSeedPolicy`

**Статус:** 7 тестов в `equilibrium_log_moles_tests2.rs` (`seed_activated_phase_trace_floor_sets_log_moles_to_floor`, `seed_activated_phase_absolute_per_species_sets_positive_moles`, `seed_activated_phase_relative_to_system_total_scales_with_inventory`, `seed_activated_phase_rejects_phase_with_no_species`, `seed_activated_phase_rejects_dimension_mismatch`, `seed_activated_phase_absolute_rejects_non_positive_moles`, `seed_activated_phase_relative_rejects_invalid_fraction`).

#### D.16 ✅ `deactivate_phases_seed_only()` с несколькими фазами

**Статус:** 5 тестов в `equilibrium_log_moles_tests2.rs` (`deactivate_phases_seed_only_sets_species_to_trace_floor`, `deactivate_phases_seed_only_handles_multiple_phases`, `deactivate_phases_seed_only_rejects_non_positive_trace_floor`, `deactivate_phases_seed_only_rejects_non_finite_trace_floor`, `deactivate_phases_seed_only_handles_empty_deactivation_list`).

#### D.17 ✅ `PhaseEquilibriumProblemBundle::solve_with()`

**Статус:** Тесты существуют в `phase_equilibrium_problem_tests.rs` (например, `local_phase_data_solves_through_one_accepted_bridge_bundle`).

#### D.18 ✅ `PhaseEquilibriumSolutionBundle::into_multiphase_solution()`

**Статус:** Тесты существуют в `phase_equilibrium_problem_tests.rs` и `equilibrium_multiphase_story_tests.rs`.

#### D.19 ✅ `MultiphaseEquilibriumSolution::moles_for()` и `mole_fraction_for()`

**Статус:** Тесты существуют в `equilibrium_multiphase_story_tests.rs` (`accepted_fixed_phase_solution_exposes_qualified_amounts_totals_and_summary`).

#### D.20 ✅ `MultiphaseEquilibriumSolution::aggregate_moles_by_substance()`

**Статус:** Тесты существуют в `equilibrium_multiphase_story_tests.rs`.

#### D.21 ✅ `TemperatureSolveSnapshot` и `TemperatureSolveFailure`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`temperature_solve_snapshot_holds_expected_fields`, `temperature_solve_failure_holds_error_message`, `temperature_solve_failure_carries_backend_attempts`).

#### D.22 ✅ `TemperatureWorkerSeed::apply()`

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`temperature_worker_seed_apply_copies_state_to_local_solver`, `temperature_worker_seed_apply_clears_previous_solution`).

#### D.23 ✅ `continuation_seed_for_point()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`continuation_seed_for_point_uses_previous_accepted_when_policy_says_so`, `continuation_seed_for_point_uses_configured_seed_when_independent`, `continuation_seed_for_point_returns_none_when_no_seed_available`).

#### D.24 ✅ `build_temperature_ranges()` и `build_temperature_gibbs_cache()`

**Статус:** 4 теста в `equilibrium_log_moles_tests2.rs` (`build_temperature_ranges_returns_single_range_on_empty_solver`, `build_temperature_gibbs_cache_returns_empty_on_empty_solver`, `build_temperature_point_index_returns_empty_for_empty_ranges`, `build_temperature_point_index_maps_temperatures_to_range_ids`).

#### D.25 ✅ `solve_temperature_point_from_seed()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`solve_temperature_point_from_seed_fails_on_empty_seed`).

#### D.26 ✅ `snapshot_published_solution()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`snapshot_published_solution_fails_without_publication`).

#### D.27 ✅ `ordered_gibbs_functions()` и `build_gibbs_functions()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`ordered_gibbs_functions_returns_error_for_missing_substance`, `ordered_gibbs_functions_returns_empty_for_empty_input`, `build_gibbs_functions_returns_empty_on_empty_solver`).

#### D.28 ✅ `build_parallel_gibbs_functions()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`build_parallel_gibbs_functions_returns_empty_on_empty_solver`).

#### D.29 ✅ `reconstructed_mole_state()`

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`reconstructed_mole_state_rejects_dimension_mismatch`, `reconstructed_mole_state_rejects_non_finite_log_moles`).

#### D.30 ✅ `check_task()` с некорректными данными

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`check_task_fails_on_empty_solver_without_settings`).

#### D.31 ✅ `has_rst_symbolic_context()`

**Статус:** 2 теста: 1 в `equilibrium_log_moles_tests.rs` (`has_rst_symbolic_context_returns_false_for_empty_solver`), 1 в `equilibrium_log_moles_tests2.rs` (`has_rst_symbolic_context_returns_false_for_empty_solver`).

#### D.32 ✅ `validate_temperature_range()` с граничными значениями

**Статус:** 4 прямых теста в `equilibrium_log_moles_tests2.rs` + 7 косвенных через `solve_for_T_range`/`par`/`par2` rejection тесты.

#### D.33 ✅ `EquilibriumLogMoles::new()` и `EquilibriumLogMoles::empty()`

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`empty_constructor_creates_default_solver`, `new_constructor_creates_default_solver`).

#### D.34 ✅ `set_initial_guess()`

**Статус:** 1 тест в `equilibrium_log_moles_tests.rs` (`set_initial_guess_validates_length`).

#### D.35 ✅ `clear_published_state()`

**Статус:** 2 теста: 1 в `equilibrium_log_moles_tests.rs` (`clear_published_state_resets_all_publication_fields`), 1 в `equilibrium_log_moles_tests2.rs` (`clear_published_state_clears_all_fields`).

#### D.36 ✅ `validate_mutable_problem_shape()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`validate_mutable_problem_shape_accepts_valid_input`, `validate_mutable_problem_shape_rejects_non_finite_moles`, `validate_mutable_problem_shape_rejects_negative_moles`).

#### D.37 ✅ `stage_equilibrium_system()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`stage_equilibrium_system_panics_on_empty_solver`, `#[should_panic]`).

#### D.38 ✅ `resolved_initial_guess()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`resolved_initial_guess_derives_seed_from_n0_when_no_explicit_guess`, `resolved_initial_guess_uses_explicit_seed_when_provided`, `resolved_initial_guess_rejects_dimension_mismatch`).

#### D.39 ✅ `build_solve_contract()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`build_solve_contract_fails_on_empty_solver`).

#### D.40 ✅ `temperature_worker_seed()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`temperature_worker_seed_getters_return_expected_values`).

#### D.41 ✅ `publish_solve_candidate()` с частичными данными

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`publish_solve_candidate_updates_all_publication_fields`, `publish_solve_candidate_accepts_missing_keq_validation`).

#### D.42 ✅ `solver_impl()` с разными политиками

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`solver_impl_fails_on_empty_solver_with_empty_policy`).

#### D.43 ✅ `solve_backend_cascade()` с пустым списком бэкендов

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`solve_backend_cascade_rejects_empty_backends`, `solve_backend_cascade_rejects_zero_budget`).

#### D.44 ✅ `recoverable_backend_failure_kind()` со всеми вариантами ошибок

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`recoverable_backend_failure_kind_returns_none_for_invalid_problem`, `recoverable_backend_failure_kind_returns_solver_for_solve_error`).

#### D.45 ✅ `compute_species_moles()` (свободная функция)

**Статус:** 3 теста в `equilibrium_log_moles_tests.rs` (`compute_species_moles_rejects_non_finite_log_moles`, `compute_species_moles_rejects_underflow_to_zero`, `compute_species_moles_round_trips_positive_moles`).

#### D.46 ✅ `compute_element_totals()`

**Статус:** 3 теста в `equilibrium_log_moles_tests.rs` (`compute_element_totals_matches_manual_calculation`, `compute_element_totals_rejects_dimension_mismatch`, `compute_element_totals_rejects_non_finite_input`).

#### D.47 ✅ `reaction_standard_gibbs()`

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`reaction_standard_gibbs_computes_weighted_sum`, `reaction_standard_gibbs_handles_multiple_reactions`).

#### D.48 ✅ `equilibrium_scaling()`

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`equilibrium_scaling_returns_positive_factors`, `equilibrium_scaling_rejects_invalid_temperature`).

#### D.49 ✅ `species_to_phase_map()`

**Статус:** 3 теста в `equilibrium_log_moles_tests.rs` (`species_to_phase_map_assigns_each_species_to_its_phase`, `species_to_phase_map_rejects_out_of_range_species_index`, `species_to_phase_map_rejects_unassigned_species`).

#### D.50 ✅ `reaction_phase_stoichiometry()`

**Статус:** 1 тест в `equilibrium_log_moles_tests.rs` (`reaction_phase_stoichiometry_aggregates_by_phase`).

#### D.51 ✅ `validate_logmole_system_dimensions()`

**Статус:** 5 тестов в `equilibrium_log_moles_tests2.rs` (`validate_logmole_system_dimensions_rejects_element_matrix_mismatch`, `validate_logmole_system_dimensions_rejects_species_phase_mismatch`, `validate_logmole_system_dimensions_rejects_out_of_bounds_phase`, `validate_logmole_system_dimensions_rejects_delta_n_mismatch`, `validate_logmole_system_dimensions_accepts_valid_dimensions`).

#### D.52 ✅ `validate_residual_conditions()`

**Статус:** 4 теста в `equilibrium_log_moles_tests2.rs` (`validate_residual_conditions_accepts_valid_conditions`, `validate_residual_conditions_rejects_non_finite_temperature`, `validate_residual_conditions_rejects_zero_pressure`, `validate_residual_conditions_rejects_negative_reference_pressure`).

#### D.53 ✅ `scaled_residual()` и `scaled_jacobian()`

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`scaled_residual_applies_row_scaling`, `scaled_jacobian_applies_row_scaling`).

#### D.54 ✅ `scale_residual_rows()` и `scale_jacobian_rows()` с некорректными scale

**Статус:** 7 тестов в `equilibrium_log_moles_tests2.rs` (`scale_residual_rows_rejects_dimension_mismatch`, `scale_residual_rows_rejects_non_positive_scale`, `scale_residual_rows_rejects_non_finite_scale`, `scale_residual_rows_applies_correct_division`, `scale_jacobian_rows_rejects_dimension_mismatch`, `scale_jacobian_rows_rejects_non_positive_scale`, `scale_jacobian_rows_applies_correct_division`).

#### D.55 ✅ `evaluate_equilibrium_logmole_residual()` с вырожденными входными данными

**Статус:** 4 теста в `equilibrium_log_moles_tests2.rs` (`evaluate_equilibrium_logmole_residual_rejects_dimension_mismatch`, `evaluate_equilibrium_logmole_residual_rejects_invalid_temperature`, `evaluate_equilibrium_logmole_residual_rejects_non_finite_log_moles`, `evaluate_equilibrium_logmole_residual_returns_correct_length`).

#### D.56 ✅ `evaluate_equilibrium_logmole_jacobian()` с вырожденными входными данными

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`evaluate_equilibrium_logmole_jacobian_rejects_dimension_mismatch`, `evaluate_equilibrium_logmole_jacobian_returns_correct_shape`, `evaluate_equilibrium_logmole_jacobian_rejects_non_finite_log_moles`).

#### D.57 ✅ `equilibrium_logmole_residual2()` и `equilibrium_logmole_jacobian2()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`equilibrium_logmole_residual2_creates_closure`, `equilibrium_logmole_residual2_rejects_wrong_log_moles_length`, `equilibrium_logmole_jacobian2_returns_correct_shape`).

#### D.58 ✅ `PhaseActivityModel::log_activity()` с разными моделями

**Статус:** Тесты существуют в `equilibrium_activity.rs` (`gas_and_solution_activity_offsets_follow_their_contracts`, `pure_condensed_phase_has_unit_activity`).

#### D.59 ✅ `EquilibriumComponentDescriptor` с пустым веществом

**Статус:** Тесты существуют в `equilibrium_problem_tests.rs`.

#### D.60 ✅ `EquilibriumConstantProblem` с некорректными данными

**Статус:** Тесты существуют в `equilibrium_constant_tests.rs` и `equilibrium_constant_solver_tests.rs`.

#### D.61 ✅ `EquilibriumConstantSolver::solve()` с некорректными границами

**Статус:** Тесты существуют в `equilibrium_constant_solver_tests.rs`.

#### D.62 ✅ `EquilibriumConstantSolver::solve()` с неограниченной задачей

**Статус:** Тесты существуют в `equilibrium_constant_solver_tests.rs`.

#### D.63 ✅ `compare_equilibrium_constant_solutions()` с идентичными решениями

**Статус:** Тесты существуют в `equilibrium_constant_cross_validation.rs`.

#### D.64 ✅ `classify_equilibrium_constant_cross_validation()` со всеми статусами

**Статус:** Тесты существуют в `equilibrium_constant_cross_validation.rs`.

#### D.65 ✅ `MultiphaseEquilibriumLayout` с некорректными PhaseSpec

**Статус:** Тесты существуют в `equilibrium_multiphase_domain_tests.rs`.

#### D.66 ✅ `MultiphaseInitialComposition::from_sparse()`

**Статус:** Тесты существуют в `equilibrium_multiphase_domain_tests.rs`.

#### D.67 ✅ `MultiphaseInitialComposition::element_totals()`

**Статус:** Тесты существуют в `equilibrium_multiphase_domain_tests.rs`.

#### D.68 ✅ `PhaseEquilibriumMetadata::from_resolved()` с пустой фазовой системой

**Статус:** Тесты существуют в `phase_equilibrium_problem_tests.rs`.

#### D.69 ✅ `build_phase_equilibrium_problem()` с некорректными условиями

**Статус:** Тесты существуют в `phase_equilibrium_problem_tests.rs`.

#### D.70 ✅ `solve_resolved_pt()` с разными `PhaseEquilibriumSolveMode`

**Статус:** Тесты существуют в `equilibrium_multiphase_story_tests.rs`.

### E. Прочие замечания

#### E.1 Магические числа

- `1e-8` в [`solve_with_phase_control`](equilibrium_workflows.rs:1617) — seed для активируемой фазы. Не вынесено в константу.
- `1e-12` в [`solve_candidate_from_seed`](equilibrium_log_moles.rs:1921) — допуск для отрицательных молей в feasibility check. Не вынесено в константу.
- `10.0` в [`solve_candidate_from_seed`](equilibrium_log_moles.rs:1937) — множитель для element_balance_tolerance. Не вынесено в константу.

#### E.2 Избыточное логирование

В [`solve_candidate_from_seed`](equilibrium_log_moles.rs:1908-1910) три `info!` подряд с полным дампом initial_guess, n0 и stoich_matrix. Для production это слишком подробно. Должно быть `debug!`.

#### E.3 Дублирование кода в `solve_for_T_range`, `solve_for_T_range_par`, `solve_for_T_range_par2`

Три реализации температурной серии содержат значительное дублирование логики обновления Gibbs, валидации и сохранения результатов. Можно было бы выделить общий шаг итерации.

#### E.4 `PHASE_CONTROL_TRACE_MOLE_FLOOR = 1e-300`

**Где:** [`equilibrium_workflows.rs:31`](equilibrium_workflows.rs:31)

Значение `1e-300` находится на грани представимости в f64 (min positive normal ≈ 2.2e-308). `ln(1e-300) ≈ -690.8`, что далеко от -Inf, но при дальнейшем делении на phase totals может привести к underflow. Рекомендуется `1e-100` или документировать риск.

#### E.5 `DEFAULT_TRACE_MOLE_FLOOR = 1e-30`

**Где:** [`equilibrium_problem.rs:28`](equilibrium_problem.rs:28)

`ln(1e-30) ≈ -69.1`. Это безопасно для f64. Замечаний нет.

#### E.6 Неконсистентность: `scaling_flag` в `EquilibriumSolverSettings` по умолчанию `false`

**Где:** [`equilibrium_log_moles.rs:153`](equilibrium_log_moles.rs:153)

Масштабирование выключено по умолчанию, хотя в архитектурном обзоре оно описано как важная фича. Рекомендуется включить по умолчанию.

---

### Приоритеты для исправления

| Приоритет | Категория | Пункт | Описание |
|-----------|-----------|-------|----------|
| P0 | Bug | A.4 | 22 публичных поля — прямой путь к багам |
| P0 | Design | C.1 | God Object — затрудняет поддержку |
| P1 | Bug | A.2 | Некорректная эвристика выбора бэкенда |
| P1 | Bug | A.1 | Неверное значение R |
| P1 | Overengineering | B.1 | Избыточный бюджет каскада |
| P1 | Overengineering | B.2 | Мёртвый код VariableScalingContract |
| P1 | Test | D.1-D.70 | 70 пробелов в тестовом покрытии — 69 закрыто, остаётся D.2 |
| P2 | Design | C.2 | Две иерархии ошибок |
| P2 | Design | C.4 | Ручное копирование полей |
| P2 | Code | E.1 | Магические числа |
| P2 | Code | E.4 | Риск underflow в PHASE_CONTROL_TRACE_MOLE_FLOOR |
| P3 | Code | E.2 | Избыточное логирование |
| P3 | Code | E.3 | Дублирование кода температурных серий |
