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
- [x] Do not add multi-start without regression evidence justifying its
  complexity. The deterministic backend cascade remains the production policy;
  multi-start is an explicitly deferred extension, not a release blocker.

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
- [x] Keep legacy backend parity tests as characterization while the handwritten
  LM/NR/TR implementations remain supported fallback backends. Their deletion
  criteria are intentionally deferred; this does not make them a second
  production orchestration path.
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

- [x] Validate named initial compositions transactionally at the
  `ResolvedPhaseSystem -> EquilibriumProblem` adapter: every component must be
  known, represented exactly once by its typed component id, and aligned with
  deterministic solver ordering. `new_with_sparse_initial_composition` now
  enforces this contract. Do not restore the legacy
  `HashMap<Option<String>, ...>` contract.
- [x] Defer a small independent Lagrange-stationarity validator unless
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

### Current production top-level priorities

The canonical fixed-`P,T` equations, backend cascade, independent K_eq
validation, bounded phase control, immutable result, and resolved-phase bridge
are already implemented. The remaining work must now be judged by whether it
makes `solve_resolved_pt` the single reliable production entry point.

#### Production-readiness audit (2026-07-26)

**Verdict:** the supported fixed-`P,T` slice is a strong production candidate,
but it is not yet the only internally canonical and externally unambiguous
workflow. The first production milestone remains deliberately limited to one
ideal-gas phase plus any number of one-component pure condensed phases.
Multiple gas phases and condensed solution models are separate physical-model
features, not defects in this milestone.

Release blockers, in execution order:

1. [x] Remove `EquilibriumLogMoles` as the mutable orchestration host from the
   canonical `solve_resolved_pt` path. Both fixed and bounded modes now use
   immutable prepared inner problems. `PreparedPhaseControlRunner` owns only
   active-set transitions, restart seeds, and reports; the historical
   `solve_with_phase_control` remains available solely for compatibility
   callers and is no longer reached by the resolved-data facade.
2. [x] Honor independent-validation policy in bounded phase control.
   `solve_fixed_active_set_candidate` currently forces
   `keq_validation_mode = Off`, and `from_phase_control_parts` publishes
   `keq_validation_status = None`. In particular, `Required` must never be
   silently downgraded. Either validate the final applicable active set or
   return an explicit typed `NotApplicable`/error according to the requested
   mode.
3. [x] Separate physical published composition from internal positive
   log-coordinate trace floors. An inactive phase currently remains visible
   through `component_moles()` and `phase_total()` as a tiny positive amount
   while its status is `Inactive`. Define and test one unambiguous result
   contract: physical public amounts are zero for inactive phases, while
   numerical trace coordinates remain available only through explicit
   `numerical_*` diagnostic accessors. Bounded mixed-phase stories now cover
   excluded phase publication and trace-floor separation.
4. [x] Finish the independent typed policy boundary. `EquilibriumSolveOptions`
   and `PhaseControlPolicy` now expose the supported configuration through
   validated typed builders. There is one source of truth for trace policy;
   `with_solver_backend` clears an explicit conflicting policy; scalar phase
   limits, explicit-set duplicates, and resolved phase-index bounds are
   rejected before the bounded loop. Raw settings/manager constructors remain
   crate-private migration ingress only, so the public prelude does not leak
   mutable backend orchestration objects.
5. [x] Add a phase-qualified sparse/named initial-composition constructor to
   `PhaseEquilibriumPipelineRequest`. The lower layer already provides
   `MultiphaseInitialComposition::from_sparse`; the top-level pipeline should
   not require users to supply an anonymous `Vec<f64>` in an ordering that is
   only known after resolution. `new_with_sparse_initial_composition` now
   provides this contract and has dense-equivalence and duplicate-entry tests.

Mandatory reliability evidence before the production label:

- [x] Record live offline continued-lifecycle cycle detection as a deferred
  evidence gap. The
  supported fixed-P,T phase models do not naturally generate a physical cycle;
  add this only when a real cyclic fixture exists rather than manufacturing
  one by corrupting the phase policy. The existing synthetic cycle tests remain
  the sharper unit layer.
  - [x] Real disappearance, hysteresis retention, budget termination, and
    failed-transition rollback are covered by `equilibrium_live_data_tests`.
  The zero-inventory `AllCandidatePhases` boundary is now normalized before
  the positive log-moles solve and covered by a live water regression. A
  separate positive-inventory boundary-recovery scenario is also covered: the
  live high-temperature water fixture deactivates liquid water through the
  typed `BoundaryUnstableActivePhase` transition instead of publishing a
  backend NaN/`AllBackendsFailed` error.
- [x] Add inventory-scale metamorphic tests over several orders of magnitude
  on the public resolved-phase facade. The acceptance gate now applies the
  explicit contract `error <= absolute_tolerance +
  relative_tolerance * abs(original_element_total)` and the live H/O fixture
  exercises both small and large inventory scales.
- [x] Add explicit tests for bounded `K_eq` modes (`Off`, `WhenApplicable`,
  `Required`) and for the typed solver/phase-policy precedence rules covered
  by the public facade. The default `Off` path is now asserted not to publish
  a validator report; `WhenApplicable` publishes either comparison evidence or
  typed `ValidatorNotApplicable`, while `Required` rejects an unsupported
  mixed-phase request.
- [x] Make retained fallback backends library-quiet by default. The reachable
  legacy NR and shared LM/NR/TR fallback internals now route iteration-level
  diagnostics through `debug`; direct production `println!` calls were removed.
  Test-only demonstrations may still print their intermediate values.
- [x] Run the ignored 20/50/100/200-species benchmark after the immutable
  orchestration path is final. The baseline is recorded in
  `PHASE_CONTROL_BENCHMARK.md`. The original debug profile grew from `2.14 ms`
  at 20 species to `8.81 s` at 200 species; a later release-build run measured
  `67.2 us` and `70.5061 ms` respectively. Both profiles are retained as
  machine/build-specific characterization, and the release run is the useful
  deployment baseline. Cache/allocation work still requires a realistic
  transition benchmark.

Required migration and documentation closure:

- [x] Migrate `chem_equilibrium_gas_example.rs` and the legacy gas section of
  `chem_equilibrium_guides.md` to the typed resolved-phase pipeline facade.
  Retained internal consumers still need a separate migration audit.
- [x] Mark `easy_equilibrium` and broad mutable workflow entry points as
  compatibility/experimental APIs. `EasyEquilibrium`, the `gas_solver*`
  helpers, and `solve_with_phase_control` now carry Rust deprecation messages. Their
  unchecked indexing, `unwrap`, and direct console output remain outside the
  production prelude; the APIs stay available for characterization and
  migration.
- [x] Update architecture documents after orchestration migration so the
  dependency graph names exactly one production entry point and clearly
  distinguishes retained numerical fallback implementations from legacy
  orchestration. The fixed-P,T boundary and compatibility status are recorded
  in `ARCHITECTURE_RU.md`.

The following do **not** block the first production milestone: GUI work,
temperature-range orchestration, multiple ideal-gas phases, condensed solution
activity models, projection caching, speculative multi-start, and deletion of
the handwritten LM/NR/TR fallback implementations.

#### P0 - Remove the remaining canonical-path correctness risks

- [x] Replace the internal mutable `EquilibriumLogMoles` host used by
  `PhaseEquilibriumProblemBundle::solve_with_bounded_phase_control` with
  canonical orchestration over immutable `EquilibriumProblem` /
  `PreparedEquilibriumProblem`, active-set projections, backend policies, and
  acceptance reports. Retain the handwritten legacy nonlinear backends only as
  explicit fallback implementations behind the common backend contract.
- [x] Build an independent element-constraint basis for every active-set
  projection. A physically admissible reduced phase set must not be rejected
  merely because the full declared element matrix contains columns that become
  dependent in that active set. Preserve the mapping back to full element
  labels and validate conservation against the original full totals. The active
  set now carries the retained independent element rank explicitly, and the
  projection accepts physically valid rank-deficient cases such as single
  `H2O(g)` H/O inventories and permuted element columns.
- [x] Add regression tests for rank-deficient but physically valid active sets,
  including the one-component `H2O(g)` H/O case, element permutations, and
  transactional failure when the reduced basis is genuinely infeasible. The
  active-set test matrix now covers acceptance of valid rank-deficient
  projections, element-order permutations, and infeasible total vectors.

#### P1 - Finish the production request and data boundary

- [x] Replace the broad legacy-facing `EquilibriumSolverSettings` and mutable
  public-field `PhaseManager` at the facade boundary with small validated
  `EquilibriumSolveOptions` and `PhaseControlPolicy` types. Backend order,
  global budget, acceptance tolerances, K_eq validation, trace policy, and
  phase hysteresis must each have one source of truth.
  - The typed wrappers now exist and are wired into the high-level
    pipeline facade. Legacy setters remain as a migration shim, but the public
    request objects no longer have to expose the raw backend controller types
    as their only shape. The old raw-settings setters are now marked
    deprecated so the typed path is visually and semantically primary, and the
    main story consumers have already moved to `EquilibriumSolveOptions` and
    `PhaseControlPolicy`, including the bounded phase-control story paths. The
    request objects now also expose an explicit reset helper for the safe fixed
    path, while the enum-based mode setter is transition-only.
- [x] Add a narrow public re-export/prelude for the supported production path:
  phase-system specification/resolution, conditions, initial composition,
  solve options, `solve_resolved_pt`, immutable solution, and typed reports.
  Internal formulation and compatibility modules must not be required imports
  for normal users. A typed `ChemEquilibrium::prelude` now re-exports the
  production path while leaving legacy modules opt-in.
- [x] Connect the existing `SubstanceSystemFactory` /
  `ThermoRepository` resolution workflow to a high-level equilibrium builder
  so callers can request `PhaseSpec -> resolve -> solve` without manually
  constructing per-phase `SubsData` maps. Keep resolution separately callable
  for users who need to inspect or reuse the immutable resolved system. A
  typed `PhaseEquilibriumPipelineRequest` now resolves the spec and can return
  both the immutable resolved system and the accepted solution as one
  transaction.
- [x] Preserve lookup provenance and NIST fallback decisions from the
  repository through build, solve, and final result reports. The resolved
  report now survives into build and solution bundles, and the public pipeline
  outcome forwards the original lookup report directly.
- [x] Decide and test the safe multiphase default explicitly:
  `FixedDeclaredPhases` and `BoundedPhaseControl` must never be selected by an
  accidental or undocumented default. The facade now exposes explicit helper
  constructors, and the default is tested to remain the fixed-declared-phase
  path.

#### P2 - Prove the top-level workflow and retire duplicate orchestration

- [x] Defer a live offline resolved-system story for cycle detection until a
  physically credible fixed-P,T fixture exists. The
  real disappearance, hysteresis, budget-termination, and failed-transition
  rollback stories are now covered; a physical fixed-P,T cycle fixture is
  intentionally not manufactured for a checklist.
- [x] Assert that the live integration layer leaves all thermochemical library
  files byte-for-byte unchanged.
- [x] Migrate the retained user-facing examples and production entry point to
  `solve_resolved_pt`; old/new parity remains confined to characterization
  tests and compatibility modules.
- [x] Deprecate the broad mutable equilibrium facade and old phase-control
  helpers after the active-set matrix passed. The explicit legacy nonlinear
  fallback backends remain supported.
- [x] Update the architecture documents and examples so only the typed
  resolved-phase facade is described as the production path.

#### P3 - Optimize only from measured evidence

- [x] Provide final canonical characterization at the supported scales and
  record projection, Jacobian, nonlinear, outer-loop, and total timings.
  Synthetic 20/50/100/200 active-set timings are recorded in
  `PHASE_CONTROL_BENCHMARK.md`; real exact-element NASA coverage exercises
  20/50/100 candidates from the local NASA gas catalog, using the allowed
  element alphabet `{C,H,O}`. The strict all-three-element subset contains
  only about 20 records and is covered separately.
  - [x] Add an ignored release story over 20/50/100 real C/H/O-limited NASA
    candidates using the same legacy-NR contract and a byte-for-byte library
    immutability check. Release evidence was captured on 2026-07-28 and is
    recorded in `PHASE_CONTROL_BENCHMARK.md`.
    The first run exposed a fixture-contract error rather than a solver limit:
    strict all-three-element matching has only about 20 local NASA gas records,
    while the allowed-element (`SubsetOf`) search supplies 100+ candidates.
    The live regression now pins both semantics and the low-level search keeps
    element sets isolated per `(library, substance)`.
- [x] Cache `ActiveSetProjection` and reduced numeric formulations for repeated
  `A -> B -> A` phase sets. The bounded range runner now retains one numeric
  and, for RST policies, one symbolic preparation per active mask; cache
  cardinality is published in the typed report and covered by live stories.

Explicitly deferred beyond this production engine milestone: GUI, multiple
ideal-gas phases, ideal/non-ideal condensed solutions, and speculative
multi-start. The typed fixed-`P,T` and temperature-range facades are in scope;
physical extensions of the supported phase models remain separate work.

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
- [x] Store an ordered equilibrium phase descriptor instead of the current
  bare `Phase { kind, species: Vec<usize> }`; derive index ranges/maps once
  from `SystemLayout`.
  - [x] `EquilibriumPhaseDescriptor` now lives beside the canonical component
    descriptor and is the source of truth in both bridge metadata and
    `EquilibriumProblem`. The legacy numeric `Phase` vector is derived once
    from descriptor ranges for existing residual/Jacobian backends.
  - [x] Retain the derived legacy numeric `Phase` projection only inside the
    numerical compatibility payload. It is not canonical identity and is not
    exposed by the production facade; removing it is tied to a future backend
    payload rewrite and is not useful release work by itself.
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
- [x] Duplicate qualified ids, duplicate phase ids, and mismatched phase ranges
  fail before thermochemistry or a nonlinear backend is invoked.
  `EquilibriumProblem::new_with_phase_descriptors` now rejects duplicate
  semantic phase ids and non-contiguous/overlapping descriptor ranges at its
  typed input boundary; tests assert that these failures occur before solver
  setup.

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

- [x] Add `MultiphaseEquilibriumSolution` containing fixed conditions, layout
  revision/fingerprint, ordered component moles, phase totals, phase-local mole
  fractions, active/inactive status, and the accepted canonical solution.
  - [x] The fixed-active bridge result is now published as an immutable
    `MultiphaseEquilibriumSolution` with fingerprint checks, phase totals,
    local mole fractions, qualified lookups, provenance, and backend evidence.
  - [x] Bridge-backed bounded phase control now publishes its accepted
    `PhaseControlledSolveReport`, final `PhaseSet` statuses, and complementarity
    acceptance report through the same immutable solution type. The numerical
    inner loop is now owned by `PreparedPhaseControlRunner`, which constructs
    immutable reduced problems from the bridge-owned data and symbolic `G0`
    snapshot; `EquilibriumLogMoles` is retained only for compatibility paths.
- [x] Provide lookups by `PhaseComponentId` and `PhaseId`; expose aggregate
  totals by bare substance only as an explicit derived view.
- [x] Retain build/lookup provenance and the complete nested backend solve
  report in the result bundle.
  - [x] When fixed-phase K_eq validation is enabled, its typed status is also
    retained and emitted as a stable result-summary row instead of being lost
    inside the mutable solver host.
- [x] Add `PhaseTransitionReport` entries with previous/new phase sets, reason,
  seed, stability metric, backend outcome, and elemental-balance evidence.
  - The immutable `PhaseControlledSolveReport` now retains typed initial/final
    phase sets, ordered transitions, explicit physical reasons, exact restart
    seeds, phase totals, driving forces, per-candidate validation evidence, the
    final validation report, and every nonlinear backend report. The remaining
    report work is integration into the final multiphase solution bundle.
- [x] Add a final `MultiphaseAcceptanceReport` combining canonical residual and
  element checks with phase stability/complementarity checks.
  - [x] The final bundle now exists and combines phase-control evidence,
    canonical validation, phase-stability reports, and a complementarity
    summary with stable rows and `Display` output.
- [x] Provide stable summary rows and `Display` output for CLI, future GUI, and
  snapshot tests; core code must not print directly.
  - [x] `PhaseControlledSolveReport` now exposes stable summary rows and a
    `Display` implementation for CLI and snapshot-friendly output.
- [x] Reject result reconstruction when the result layout does not match the
  resolved-system fingerprint.
  - [x] `MultiphaseInitialComposition` already refuses reconstruction against
    a foreign layout fingerprint, and the regression tests now pin that
    contract down explicitly.
  - [x] `MultiphaseEquilibriumSolution` now independently rejects metadata and
    build provenance originating from distinct valid resolved layouts before
    publishing a queryable result.

### P4.7 Complete the multiphase test matrix

- [x] Add a dedicated `equilibrium_phase_bridge_tests.rs` for identity,
  ordering, adapter validation, provenance, and transactional failures.
- [x] Add a dedicated `equilibrium_multiphase_story_tests.rs` with module-level
  documentation describing each physical hypothesis and expected result.
  - [x] Cover the remaining live transition scenarios supported by the current
    fixed-P,T model contract:
  - [x] ideal-gas-only parity with the current canonical solver;
  - [x] gas plus one stable pure solid;
  - [x] gas plus one stable pure liquid;
  - [x] the same molecule represented in gas and condensed phases;
  - [x] two independent pure condensed candidate phases;
  - [x] unsupported multi-component condensed solution;
  - [x] phase appearance on live resolved thermochemical data;
  - [x] disappearance within one continued solve on live resolved
    thermochemical data, including positive-inventory liquid boundary recovery;
  - [x] hysteresis retention on live resolved thermochemical data;
  - [x] cycle detection remains covered by deterministic synthetic state-machine
    regressions; a live fixture is explicitly deferred until credible physics
    produces one;
  - [x] offline mixed NASA gas/NASA condensed lookup provenance.
  - [x] The dedicated equilibrium_phase_bridge_tests.rs module now covers real offline NASA gas/condensed provenance, a real NASA gas plus solid H2O(s) fixture, the same molecule as two phase-qualified components, bounded inactive-liquid startup, and transactional no-mutation-on-failure publication.
  - [x] The new `equilibrium_live_data_tests.rs` module adds a separate live-data layer over the bundled local thermochemistry repository: explicit offline lookup policy, real NASA gas resolution, one stable real gas solve, one bounded real multiphase solve, and transactional non-mutation of the resolved source view.
  - [x] The live-data layer now also includes a real water gas/liquid
    temperature-shift regression: 350 K carries a non-trace liquid inventory,
    550 K shifts the inventory toward vapor, both solves stay inside the
    `H2O(L)` record range, and layout identity remains unchanged.
  - [x] The live bounded pipeline now snapshots the canonical substance-base,
    address-catalog, and element-composition JSON files before solve and proves
    that all three remain byte-for-byte unchanged afterward.
  - [x] A live `NASA_gas` H2/O2/H2O reaction fixture now solves through the
    public pipeline with element conservation and an accepted independent
    `K_eq` cross-validation report. `EquilibriumSolveOptions` exposes this as
    an explicit typed opt-in rather than requiring callers to mutate backend
    settings directly.
  - [x] The same reactive fixture pins stable fingerprints of its three source
    NASA7 records. Updating local thermochemical data now makes the fixture
    fail loudly until its physical contract is deliberately reviewed.
  - [x] A live explicit-state lookup now resolves `gas::H2O` and
    `solid::H2O(s)` as distinct phase-qualified records from the local
    NASA_gas/NASA_cond catalog, with NIST fallback disabled.
  - [x] Fix phase-qualified live lookup and thermochemistry preparation.
    `PhaseSpec::physical_state` is now installed as a typed `SubsData` lookup
    constraint, and bridge construction selects the NASA coefficient interval
    at the requested temperature before snapshotting numeric/symbolic `G0(T)`.
    Previously liquid water could resolve as `NASA_gas:H2O`, while freshly
    parsed NASA calculators silently published zero-coefficient `G0(T)=0`.
  - [x] Replace the invalid ice diagnostic built on zero `G0`. The real
    `H2O(g)/H2O(s)` fixture at 250 K resolves
    `NASA_gas:H2O` plus `NASA_cond:H2O(s)`, activates the initially absent ice
    phase, and transfers almost all water inventory into the solid while
    preserving the gas oxygen inventory.
  - [x] Defer Fe phase-lifecycle fixtures until defining an exact-record
    selection policy for condensed polymorphs (`Fe(a)`, `Fe(c)`, `Fe(d)`).
    The earlier zero-driving-force observation is invalid because it was made
    before the zero-`G0` bridge defect was fixed; do not retain it as physical
    evidence or let a state-only resolver guess a polymorph.
  - [x] Promote a physically scaled Boudouard fixture
    `2 CO(g) <=> CO2(g) + C(gr)` into the live matrix. A finite CO/CO2 initial
    mixture avoids trace-floor chemical potentials, activates real
    `NASA_cond:C(gr)` at 700 K, and shows lower graphite inventory at 1400 K.
    The deliberately harsher all-CO boundary case remains a separate
    initialization benchmark rather than the production lifecycle fixture.
  - [x] The local NASA_gas H2O/O2 plus NASA_cond H2O fixture now reaches
    the bounded bridge with the zero-inventory liquid phase initially inactive.
    The outer loop derives its first phase set from physical `n0`, not the
    positive numerical trace coordinate, and solves the square gas-only
    projection without relaxing acceptance tolerances.
  - [x] `AllCandidatePhases` now has the same physical boundary contract:
    inventory-free candidates are normalized to inactive before a positive
    log-moles solve, while remaining eligible for stability-driven activation.
    The live high-temperature water story exercises this path without a
    synthetic trace inventory.
  - [x] The public water phase-pair story now shows temperature-driven vapor /
    condensed dominance on the real `solve_resolved_pt` facade rather than only
    in synthetic workflow tests.
  - [x] Accept that cycle detection during one continued live run is pinned
    primarily by synthetic workflow regressions. Positive-inventory
    disappearance, hysteresis retention, budget termination, and rollback are
    now covered by real local thermochemistry; these must not be conflated
    with zero-inventory normalization.
    Ice, liquid water, and graphite now provide stable offline appearance and
    temperature-shift fixtures for extending that matrix without inventing
    thermochemistry.
  - [x] Keep rank-deficient-but-representable active sets valid. The active-set
    projection derives the independent reaction/element rank and separately
    verifies that the closed-system element totals lie in the range of
    `A_active^T`; it no longer equates the number of declared element columns
    with the number of independent constraints.
- [x] Define an explicit element-candidate search contract for real-data
  equilibrium fixtures. `ThermoData` and `SubsData` now expose
  `ElementSearchMode::{AnyRequested, SubsetOf, ExactSet}` and a dedicated
  `search_by_exact_elements` method. The historical `search_by_elements_only`
  remains the broader subset-of-elements mode; it is not silently redefined.
- [x] Build a production candidate-selection layer on top of element search.
  It must intersect the chosen element mode with ordered thermochemical
  library preference, physical-state/phase policy, temperature-interval
  support, and per-record provenance before constructing an equilibrium
  problem. `EquilibriumCandidateSelector` now performs this as a read-only
  deterministic transaction and returns selected records plus explicit
  rejection reasons. It chooses one acceptable record per substance according
  to library preference, preserves physical-state evidence and record keys,
  screens recognizable coefficient intervals, and leaves unknown interval
  schemas visible as `Unknown` rather than silently rejecting them.
  - [x] Add `EquilibriumCandidatePhasePlan` and the
    `PhaseEquilibriumPipelineRequest::from_candidate_selection` boundary.
    Every selected exact `record_key` must be assigned once to an explicit
    phase/model; the selected library is pinned through explicit lookup
    instructions. The builder rejects omitted, duplicate, or unknown records
    transactionally and never guesses a phase model from element data.
- [x] For every applicable small fixed-phase set, compare the accepted result
  with the independent equilibrium-constant solver. Report non-applicability
  explicitly for phase-appearance decisions that the K_eq problem does not
  model. The gas-only fixed-phase bridge now retains a `Compared(report)` K_eq
  validation status in the immutable summary, and the acceptance evidence is
  checked in the dedicated story test.
  - [x] Keep `WhenApplicable` validation observational: failure of the
    independent K_eq solver is retained as `ValidatorFailed` evidence and does
    not reject an otherwise accepted canonical solution. `Required` remains
    the explicit fail-closed mode.
  - [x] Scale the default total-Gibbs comparison tolerance for live problems.
    The validator still requires close species amounts and fractions, while
    avoiding false rejection from a few microjoules of absolute objective
    drift on megajoule-scale solutions.
- [x] Add exact component-order and numeric/closure/symbolic equivalence
  tests, then run the existing RST backend/fallback matrix over at least one
  physical multiphase fixture. The real NASA gas + condensed-water bridge test
  now compares prepared residual/Jacobian, legacy closure/Jacobian, and
  symbolic residual/Jacobian on one real multiphase fixture, while the
  bounded phase-control matrix still pins the explicit legacy-vs-RST policy
  parity on the same physical system.
  - [x] A physical bounded phase-control fixture now compares explicit legacy
    and explicit RST solver policies on the same resolved multiphase system and
    confirms that the accepted component order and moles remain aligned.
  - [x] Make the implicit production solver policy RST-first while retaining
    LM/NR/TR legacy fallbacks. Explicit single-backend policies remain strict
    and never acquire an undeclared fallback.
- [x] Adopt physically meaningful fixtures from the retired classical stack
  only after restating their expected invariants. Oxygen dissociation,
  reactive H/O, water phase pairs, ice, and graphite now test conservation,
  residual, provenance, and lifecycle contracts; fixtures tied only to legacy
  mutable state or one historical iterate were not retained.
- [x] Keep all default tests offline and leave thermochemical library files
  byte-for-byte unchanged. The dedicated live-data regression snapshots every
  canonical JSON file used by normal thermo lookup around a real bounded solve.

### P4.8 Migrate workflows and retire the duplicate engine

- [x] Add one public one-shot facade such as `solve_resolved_pt` accepting the
  resolved phase system, typed initial composition, conditions, solver policy,
  and phase-control policy.
  - [x] `ResolvedPhaseEquilibriumRequest` and `solve_resolved_pt` now own the
    fixed-declared-phase bridge transaction, including typed numerical
    settings and immutable result publication.
  - [x] `PhaseEquilibriumSolveMode::BoundedPhaseControl(PhaseManager)` now
    routes through the same bridge-owned thermochemistry and publishes one
    `MultiphaseEquilibriumSolution` with phase-control and complementarity
    evidence. The mutable solver remains an internal compatibility host, not
    an alternative public workflow.
- [x] Migrate useful legacy equilibrium workflows and examples one
  physical scenario at a time onto that facade.
  - [x] The local NASA ideal-gas guide now uses `ResolvedPhaseSystem`, typed
    initial composition, and `solve_resolved_pt`; it is the first retained
    user-facing scenario that does not call legacy `gas_solver` directly.
  - [x] The resolved-phase example now uses the repository-backed typed
    pipeline as well; callers no longer need to assemble `SubsData` maps or
    `ResolvedPhaseSystem` manually for the standard guide scenario.
  - [x] A source audit found no non-test consumer outside the compatibility
    module that still calls `gas_solver*`, `solve`, or
    `solve_with_phase_control`; remaining calls are characterization coverage.
    The mutable temperature-range methods remain deprecated now that the typed
    range facade is available.
- [x] During migration, compare both implementations only in characterization
  tests; do not expose two production APIs as equivalent long-term choices.
- [x] Add a typed temperature-range facade above the fixed-`P,T` production
  solve. Resolve the real-data candidate set once, reuse the element matrix,
  reaction basis, and active-set projection while the layout is unchanged,
  refresh temperature-dependent standard-state data for each point, use the
  previous accepted physical solution as the continuation seed, and publish
  each point transactionally. Rebuild only after an accepted phase transition.
  - [x] Make the range request own an immutable resolved layout and a typed
    point grid; do not route the canonical range path through the mutable
    `EquilibriumLogMoles` sweep. `TemperatureRangeRequest` and
    `PhaseEquilibriumPipelineRequest::solve_temperature_range` now provide
    this fixed-declared-phase path.
  - [x] Add per-point timing records for repository lookup, coefficient
    refresh, numeric closure refresh/construction, symbolic refresh or reuse,
    equation/Jacobian preparation, nonlinear solve, phase control, validation,
    and postprocessing. Keep a sweep summary with total, mean, median, and
    worst-point durations. The fixed-declared range now publishes the full
    typed point timing report produced by the solve and a
    `TemperatureRangeDurationSummary`; bounded points additionally publish
    transition counts, phase-set reuse, and projection-cache cardinality.
  - [x] Distinguish one-time setup from per-point work in the report. The
    initial baseline of approximately 445.9 ms (NR), 457.9 ms (LM), and
    466.2 ms (TR) for five points was diagnosed as an accidental RST
    promotion: the compatibility wrapper left `solver_policy` unset even
    after receiving a legacy `Solvers` selector. After the boundary fix,
    the same debug run measured approximately 5.6 ms (NR), 9.2 ms (TR), and
    16.2 ms (LM) for five points. The typed report now separates
    `initial_formulation_timing` from accepted-point timing and exposes
    total/mean/median/worst point durations. The release baseline is recorded
    below for the retained real-data characterization story.
    - [x] Release characterization is now recorded for the real 20-species,
      three-point story: legacy NR took 494 us ascending / 409.4 us descending
      across accepted points; default RST took 106.5217 ms ascending. The RST
      run reused symbolic preparation for 2 of 3 points; the first point
      crosses a real NASA coefficient interval boundary. These are measured
      point totals, not a claim that the two backend families have equivalent
      setup costs.
  - [x] Add a regression comparing one fixed-`P,T` solve with the equivalent
    one-point range solve. The real local fixture compares accepted moles under
    the same legacy backend and keeps any setup cost visible in the timing
    report.
  - [x] Reuse the element matrix, phase/component layout, reaction basis, and
    allocation buffers when the active phase set is unchanged. Add counters
    for reused versus rebuilt objects so a faster result cannot be caused by
    accidentally dropping required work. The fixed range reports one
    formulation build and point-level formulation reuse counters.
  - [x] Refresh only temperature-dependent standard-state data at each point;
    measure numeric closure rebuild separately from symbolic expression reuse.
    Do not assume symbolic expressions are reusable if their captured values
    are temperature-specific; the current reusable RST path updates the
    shared `T` parameter and records that decision explicitly. Coefficient
    interval changes still require the symbolic-refresh branch below.
  - [x] Use the previous accepted physical solution as the continuation seed
    and test both ascending and descending grids. Record seed source and
    rejected/failed points without publishing a partial sweep as complete.
    The ignored real-data release story covers both directions and verifies
    that a failed point cannot publish a partial `TemperatureRangeSolution`.
  - [x] Add a real-data backend matrix over the same temperature grid for RST
    LM/Minpack/Nielsen/trust-region/Powell/Newton and legacy LM/NR/TR. The
    ignored release story records success/failure, per-point residual,
    conservation error, and timing; every run uses `SolverPolicy::Single`, so
    a cascade cannot hide a backend-specific failure. Target-machine release
    measurements remain an operator-run evidence step.
  - [x] Add phase-control range stories with an accepted phase transition and
    verify that projection/layout rebuild counters increase only after that
    transition. A temperature point without a transition must not rebuild the
    full active-set problem. The real offline water/ice fixture now covers the
    transition path, while the bounded live story covers the no-transition
    reuse path and reports both projection and reduced-formulation caches.
  - [x] Extend bounded-range reuse to the RST symbolic backend itself. Each
    active-set cache entry retains one `RstPreparedProblem`; unchanged Gibbs
    expressions update only its shared `T` parameter, while a new active mask
    or symbolic snapshot creates a new entry. The range report exposes the
    retained symbolic-entry count and live stories assert it.
  - [x] Add a scale-aware conservation and residual contract for every bounded
    range point, plus a byte-for-byte no-file-mutation check for the live local
    repository.
  - [x] Add the release characterization in release mode over at least 20,
  50, and 100 real species where exact element search supplies enough local
  records; retain the smaller five-point story as the fast regression
  fixture. The printed timing baseline still has to be collected on the
  target machine.
- [x] Make the retained `gas_solver_for_T_range` compatibility boundary honor
  its explicit legacy `Solvers` argument. It now installs
  `SolverPolicy::legacy_default` instead of allowing symbolic context to
  select the RST production cascade implicitly. The real T-range story pins
  that every accepted point reports a legacy backend.
- [x] Add an opt-in typed timing report to the canonical fixed-`P,T` solve.
  `EquilibriumTimingMode::Enabled` publishes repository lookup,
  thermochemistry, numeric-closure, symbolic, equation, numerical-preparation,
  projection, nonlinear-solve, phase-control, validation, postprocessing, and
  total durations without logging or clock reads on the default disabled path.
   - [x] Extend the timing report to the typed temperature-range facade. Record
  per-point durations plus cache/layout/projection reuse and rebuild counters;
  distinguish coefficient refresh from a genuine active-phase transition so
     range benchmarks expose the cost of changing temperature rather than only
     the cost of one fixed-`P,T` solve. Fixed-declared and bounded-range timing,
     formulation-reuse counters, active-set cache cardinality, symbolic RST
     cache cardinality, and total/mean/median/worst point summaries now exist.
- [x] Add a typed postprocessing adapter from `TemperatureRangeSolution` to the
  immutable interpolation/report layer. It preserves phase-qualified component
  labels and physical component amounts, supports descending solver grids by
  sorting a copied view, and keeps raw-row helpers as a lower-level API.
- [x] Move the direct mutable orchestration modules behind crate-private
  implementation boundaries and expose their retained compatibility surface
  only through `ChemEquilibrium::legacy`; LM/NR/TR numerical fallback code is
  intentionally retained there. A compatibility test now resolves
  `legacy::gas_solver` through that namespace, and the guide no longer presents
  the historical module paths as peer production APIs.
- [x] Migrate production equilibrium consumers away from the broad
  `ThermodynamicsCalculatorTrait`, raw Gibbs closures, and nested legacy phase
  maps to the narrow typed bridge.
  A source audit found no non-test production consumer of those representations;
  they remain only in the explicitly namespaced compatibility subsystem.
- [x] Deprecate `solve_with_phase_control` and the old mutable phase helpers
  after the active-set story matrix passes. The deprecation is advisory: the
  legacy implementation remains available as a compatibility host, while the
  resolved facade never calls it.
- [x] Remove duplicate mutable orchestration from the production public surface.
  The old workflow remains available only under `ChemEquilibrium::legacy` by
  explicit compatibility policy. Keep the handwritten legacy nonlinear
  backends as explicit policy-selected fallback implementations.
- [x] Update the phase-subsystem architecture documents with the final
  dependency direction and add a fixed-P,T multiphase usage example. The
  retained legacy dataflow is labeled historical compatibility rather than
  presented as a second production architecture.

### Recommended P4 implementation passes

1. [x] Component/phase identity refactor in `EquilibriumProblem` plus tests.
2. [x] Typed initial composition and pure `ResolvedPhaseSystem` adapter.
3. [x] Standard-state thermochemistry/provenance bridge and preview report.
4. [x] Fixed-active-set numeric, symbolic, RST, and legacy-fallback parity.
5. [x] Bounded phase-stability active-set orchestration.
6. [x] Phase-aware solution/report boundary and full acceptance gate.
7. [x] Initial offline physical story matrix and independent K_eq
   cross-validation. The remaining live transition matrix is tracked in the
   production P2 queue above.
8. [x] Classical workflow migration and production-boundary cleanup. Retained
   compatibility workflows are isolated under `ChemEquilibrium::legacy`.

### P4 definition of done

- [x] One typed request can solve a resolved closed multiphase system at fixed
  `P,T` without flattening phase-qualified identity.
- [x] Standard-state thermochemistry is evaluated once per component and every
  record retains lookup provenance.
- [x] All supported activity terms are explicit and tested; unsupported phase
  models fail before a backend starts.
- [x] Phase appearance/disappearance is bounded, transactional, conservative,
  and included in the acceptance report.
- [x] The result can be queried unambiguously by phase and component and is
  tied to the exact resolved layout that produced it.
- [x] Numeric, analytic-Jacobian, and symbolic formulations agree; the RST
  backend matrix and explicit legacy fallback pass the same physical fixture.
- [x] Offline gas/condensed stories pass elemental, residual, stability, and
  applicable independent K_eq validation gates.
- [x] No retained production workflow requires the duplicate classical
  or mutable equilibrium orchestration. The handwritten legacy nonlinear
  backends remain supported only as explicit fallbacks.

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
- [x] Cache `ActiveSetProjection` and surrounding fixed-layout preparations for
  repeated active sets, including `A -> B -> A`. Bounded temperature sweeps
  expose projection, prepared-formulation, and symbolic cache cardinalities.
- [x] Add a separate benchmark suite for larger synthetic systems
  (20, 50, 100, 200 species) that reports outer iterations, nonlinear
  iterations, projection build time, Jacobian time, and total solve time.
  The ignored benchmark sweep now covers 20, 50, 100, and 200 species
  synthetic inventories and prints projection-build and solve timing together
  with iteration counts.
- [x] Do not rewrite allocation-heavy paths without measured evidence. Current
  characterization identifies nonlinear solve/backend setup as the dominant
  real-data cost; further allocation work is explicitly performance-driven.
- [x] After the phase* bridge was connected, add live resolved-phase-system
  regression fixtures and gradually migrate selected synthetic stories to
  those real inputs. Water gas/liquid/ice, graphite, reactive H/O, lookup
  provenance, conservation, lifecycle, and transactionality are now covered;
  synthetic cases remain as the fast unit layer.
- [x] Add an opt-in real large-system characterization: select 20
  C/H/O-limited candidates from the local element catalog, resolve them from
  `NASA_gas`,
  solve the fixed-`P,T` system, and print the stage timing report. This is a
  performance/diagnostic fixture, not a replacement for the smaller stable
  physical regressions.
- [x] Run the same real 20-species problem through each concrete RST and
  legacy backend as isolated `Single` policies. Keep backend failures visible
  instead of hiding them behind the production cascade, and require every
  accepted backend result to satisfy finite-mole, residual, and scale-aware
  elemental-balance checks.
- [x] Add a real 20-species, five-point temperature-range characterization for
  legacy LM/NR/TR. It verifies ordered accepted points, no failed temperatures,
  finite positive moles, preserves per-point backend reports, rejects hidden
  RST promotion, and records one sweep duration per backend. This is
  deliberately labeled compatibility coverage until the typed range facade
  replaces the mutable sweep.
- [x] Add the corresponding real-data temperature-range characterization after
  the typed range facade is finalized; report per-point coefficient refresh,
  continuation, rebuild, and solve timings. The ignored live story
  `live_large_element_limited_typed_temperature_range_story` covers real local
  NASA data in both directions and checks one formulation build, continuation,
  symbolic reuse, and point reports. The release run now passes and prints the
  measured point timing baseline; a broader 20/50/100-species release matrix
  remains separate characterization work.
- [x] Add an ignored release characterization over 100 real local NASA
  C/H/O-limited species and 50 temperature points, running every supported RST
  and retained legacy backend as an isolated `Single` policy. The test reports
  per-backend total/mean/median/worst timings, formulation builds and reuses,
  symbolic updates, residuals, conservation error, continuation usage, and
  preserves the JSON-library snapshot. Rebuilds are allowed at genuine NASA
  coefficient-interval boundaries; the contract requires complete build/reuse
  accounting rather than assuming one formulation for the whole grid. The
  target-machine release run is now recorded in `PHASE_CONTROL_BENCHMARK.md`:
  six backends succeeded and three reported explicit first-point failures.
- [x] Diagnose backend-specific failures and rare slow points from the 100 x 50
  release matrix before changing the production default. The range error
  boundary now preserves the boxed `ReactionExtentError` instead of flattening
  the backend trace into a string; the ignored 100 x 50 story prints the three
  slowest accepted points for each successful backend with seed-to-result
  log-coordinate movement, detailed stage timings, and RST counters, plus a
  separate typed attempt table for every failure. The strict monolithic P,H
  matrix additionally records failure point and `ReactionExtentErrorKind`.

  A concrete defect has now been corrected in the canonical fixed-`P,T` RST
  path: symbolic methods previously received unbounded log-mole coordinates,
  so a trial step could make `exp(y)` overflow before the common candidate
  validator ran. `PreparedEquilibriumProblem` now derives finite per-species
  box bounds from exact elemental capacities and the caller's trace seed, and
  passes them through `RustedSciTheSolveContract` into RST. The focused live
  20-species C/H/O backend matrix now accepts Nielsen LM, Powell Dogleg, and
  legacy TR; the old Nielsen non-finite residual no longer reproduces there.
  Damped Newton instead exposes an independent stagnation/high-residual
  rejection, which remains a valid method-specific outcome.

  The enhanced target-machine report isolated the common 1010 K anomaly:
  RST residual/Jacobian callbacks and the nonlinear solve took only about
  10--20 ms, while KiThe rebuilt `SymbolicNonlinearProblem` after a NASA
  coefficient-interval change. The canonical fixed-`P,T` graph now carries
  `G0_i` as checked equation parameters. A range point updates `[T, G0...]`
  in place, so an unchanged active set may never rebuild its symbolic graph;
  the same rule now applies to the cached bounded phase-control runner.
  Initial graph construction is recorded in `initial_formulation_timing`, and
  live 20-species stories assert zero per-point formulation build across an
  unchanged layout / active set. The 100 x 50 operator table now prints
  initial setup, initial symbolic construction, and initial numerical-problem
  preparation separately from later per-point formulation builds. The bridge
  builder now also records its outer `total`; a focused local-NASA test guards
  the invariant that setup total encloses every recorded setup stage. A second
  local-NASA contract compares residual and Jacobian entries of the
  parameterized RST graph directly against `PreparedEquilibriumProblem`, so a
  performance change cannot silently alter the fixed-`P,T` equations.

  - [x] Re-run the 100 x 50 release matrix after parameterizing `G0_i`.
    The recorded release table confirms 49 reuses, 49 parameter updates, and
    zero later `formulation_build` for every successful RST backend. The old
    9--11 second coefficient-boundary point is gone; the remaining roughly
    4.2--4.3 second cost is the explicitly reported one-time symbolic setup.
    A test-only baked-vs-parameterized 20-species Damped-Newton comparison
    rejects both graphs under the same strict gate, while an entrywise
    residual/Jacobian contract proves the reusable graph matches the canonical
    real NASA formulation. Rebuilding the historical baked 100-species graph
    is itself too expensive for a practical regression (it exceeded two
    release minutes), which is precisely why it is not retained as a fallback.
  - [ ] Classify the remaining strict-single-backend outcomes without weakening
    acceptance: the latest 100 x 50 release matrix shows Nielsen LM stopping
    at the first point with a finite residual after `MaxIterations`, Damped
    Newton stagnating after two iterations with a high residual, and legacy TR
    exhausting its iteration budget. The Damped Newton 20-species real-data
    baked-vs-parameterized regression already rules out a reusable-graph
    regression for that method. Add minimal regression only for behaviour that
    changes unexpectedly (especially renewed non-finite residual), rather
    than forcing every method to be a production default.

  RST attempt metrics split residual-callback time, Jacobian-callback time,
  and remaining engine overhead; the enhanced report prints that split for
  each slow accepted point and typed failure attempt. Legacy methods
  deliberately leave this optional evidence absent rather than inventing
  incomparable counters. Do not weaken the common acceptance gate or select a
  universal default from one inventory.

- [x] Add an explicit fixed-`P,T` multi-start recovery policy. The typed
  `ResolvedPhaseEquilibriumRequest` and high-level pipeline now accept an
  ordered list of `LogMolesInitialGuess` values, reuse the prepared formulation
  and one symbolic RST problem across seeds, compare all accepted candidates
  through the common validation ordering, and publish selected-seed evidence in
  both solution bundles. Multi-start is intentionally rejected for bounded
  phase control and temperature ranges until those workflows have a separate
  lifecycle-aware policy; it is not silently applied to an incompatible outer
  loop.
- [x] Keep the live-cycle evidence boundary explicit. Real local water/ice and
  liquid fixtures cover appearance, disappearance, hysteresis retention,
  rollback, and budget termination. A physically credible continued run that
  actually revisits an active phase set has not been observed; synthetic cycle
  state-machine tests remain the mandatory correctness layer, while a live
  cycle remains deferred rather than being manufactured with an invalid policy.

## P5 - Public API, GUI, and cleanup

The original P5 checklist below is retained as historical planning context;
the fixed-P,T public API items that were completed during P4 are synchronized
with their current status here. GUI and deferred physical models remain open.

- [x] Provide a small public facade: validated problem builder, solver policy,
  solve method, solution, solve report, and optional validation report.
  A dedicated `equilibrium_public_api_tests` module imports only
  `ChemEquilibrium::prelude` and covers fixed-`P,T`, one-point range parity,
  typed phase-qualified lookup, duplicate sparse-inventory rejection, backend
  policy selection, element-selection policy, and physical phase assignment.
- [x] Expose backend selection and cascade diagnostics in GUI only through typed
  controls; never require users to type internal enum names. The editor uses
  `GuiSolverBackend` and ordered typed cascade controls; egui tests cover every
  concrete backend and duplicate rejection.
- [x] Show conservation, residual, fallback-attempt, and K_eq validation status
  as first-class result sections. The result view reads the immutable accepted
  report and does not reconstruct solver diagnostics in the GUI.
- [x] Complete GUI story tests for success, fallback success, all-backends-
  failed, invalid input, validation mismatch, and save/load roundtrip.
  - [x] Success, invalid editor input, worker failure/rollback, lifecycle, and
    document roundtrip stories are covered.
  - [x] The ignored local fallback story covers accepted backend attempts and
    renders the diagnostic sections.
  - [x] The ignored local H2/O2/H2O P,H stories cover inner fallback and a
    true engine `AllBackendsFailed` publication failure; a generic worker
    error is not used as a substitute.
  - [x] Add a deterministic accepted-result layout-mismatch fixture. The
    ignored GUI result-layer story solves two real NASA systems, combines their
    accepted snapshots as one range payload, and verifies transactional
    rejection before any partial table can be published.
- [x] Split the implementation by responsibility: `domain`, `formulation`,
  `validation`, `backend`, `solver_policy`, `keq_validation`, and `report`.
- [x] Review `easy_equilibrium.rs` as the future facade: it is rejected as the
  production facade because it models only one reaction and has no phase-aware
  acceptance contract; it is retained as a deprecated compatibility helper.
- [x] Rename the canonical modules and their tests so public paths no longer
  imply that the main formulation is only an equilibrium-constant solver.
- [x] Keep the production facade free of direct console output and duplicate
  temperature-sweep orchestration. `Untitled-1.rs` is deleted; timing output is
  confined to opt-in ignored characterization tests. The retained
  single-reaction legacy helper no longer emits `println!`/`dbg!` output.

## P6 - Production fixed-pressure, fixed-enthalpy equilibrium

This stage adds the closed-system `P,H = const` formulation after the fixed
`P,T` engine has become the canonical production core. Temperature is an
unknown result, while pressure, total elemental inventory, and total enthalpy
are prescribed.

The first implemented path is a bounded scalar temperature solve around the
existing accepted `P,T` solve:

```text
target H, pressure, temperature bracket
                    |
                    v
        evaluate F(T) = H_eq(P, T) - H_target
                    |
                    v
        canonical solve_resolved_pt at trial T
                    |
                    v
      accepted composition + phase-control report
```

This path remains valuable as an independent reference and safeguarded
fallback. It must not, however, be presented as the final production
architecture: the canonical target is a coupled composition-temperature
system with an analytic Jacobian.

### P6.A Architectural correction: monolithic P,H is the production target

Current dependency map (30.07.2026):

```text
ResolvedPhaseEnthalpyRequest
  |
  +-- ResolvedThermochemistry / EnthalpyModel
  |
  +-- solve_resolved_ph
        |
        +-- fixed declared phases:
        |     PreparedMonolithicPhRunner          [enabled, RST symbolic + analytic legacy backends]
        |
        +-- nested/reference path:
        |     solve_bracketed_temperature_from
        |       |
        |       +-- PhTrialEvaluator::evaluate(T)
        |             |
        |             +-- PreparedPhaseEquilibriumTemplate::solve_at(T)
        |
        +-- bounded monolithic path:
              PreparedPhaseControlTemplate
                |
                +-- PreparedPhaseControlRunner
                      |
                      +-- solve_monolithic_active_set_candidate
```

The fixed-`P,T` equations that must be reused rather than copied are:

- `PreparedEquilibriumProblem::residual` and
  `evaluate_equilibrium_logmole_residual` for the canonical reaction and
  element rows;
- `PreparedEquilibriumProblem::jacobian` and
  `evaluate_equilibrium_logmole_jacobian` for the analytic log-mole block;
- `ResidualScalingContract` for matching residual/Jacobian row scaling;
- `PreparedEquilibriumRunner` and `EquilibriumNonlinearBackend` as the current
  backend-cascade boundary. The dedicated P,H RST payload now carries the
  coupled `(ln(n), theta_T)` residual rather than reusing a fixed-P,T problem;
  its exact symbolic route is deliberately limited to one native coefficient
  interval per component. A range crossing stays on the analytic P,H path
  until a piecewise-symbolic contract is designed explicitly.

The active-set boundary is `PreparedPhaseControlRunner::solve`: phase
activation/deactivation, hysteresis, projection caches, and transition reports
remain outside the nonlinear fixed-active-set solve. P,H must provide one
monolithic fixed-set candidate solver to this outer loop; phase-control state
must never enter the P,H residual.

The existing P,H code is split by ownership as follows:

- `equilibrium_ph_thermochemistry.rs`:
  `ThermochemistryProvenance`, `ResolvedThermochemistry`,
  `MolarThermoFunction`, `EnthalpyModel`, `EnthalpyEvaluation`, and private
  property/capability construction;
- `equilibrium_ph_formulation.rs`:
  monolithic unknown/row layouts, bounded temperature transform, P,H
  residual/Jacobian evaluation, row scaling, and dimension validation;
- `equilibrium_ph_monolithic.rs`:
  prepared fixed-active-set P,H problem, backend cascade, candidate
  reconstruction, and monolithic diagnostics;
- `equilibrium_ph_nested.rs`:
  the current safeguarded bracket, trial evaluator, nested budgets,
  monotonicity policy, and nested evidence;
- `equilibrium_ph_workflow.rs`:
  public request/result facade, `PhSolveMode`, classified fallback policy,
  phase-control integration, and final publication validation.

Reference tests retained during migration:

- `nonreacting_sensible_heat_fixture_recovers_the_analytic_root`;
- `inventory_scaling_preserves_temperature_and_mole_fractions`;
- `sampled_non_monotone_branch_is_rejected_by_default_but_explicit_compatibility_allows_it`;
- `ph_temperature_seed_is_a_real_trial_and_can_accept_the_root`;
- `bracketed_solver_rejects_unreachable_target`;
- live `live_reactive_pt_to_h_to_ph_recovers_temperature_and_composition`;
- live `live_bounded_water_pt_to_h_to_ph_preserves_phase_evidence_and_json`;
- live `live_reactive_gas_ph_backend_matrix`.

Migration contract:

- [x] Introduce `PhSolveMode::{Monolithic, NestedTemperature, Auto}`. The
  resolved-thermochemistry facade now defaults to `Monolithic`, the coupled
  production target. The generic closure constructor remains explicitly
  nested because it cannot provide the Gibbs/Cp bundle required by the
  coupled formulation; `NestedTemperature` remains available as an
  independent reference route and `Auto` as the classified recovery mode.
- [x] Make `Auto` run monolithic first and fall back solely for classified
  numerical failures. Input, thermochemistry, dimension, and
  unsupported-physics errors never trigger fallback; the immutable report
  retains the typed fallback reason. The default does not hide a formulation
  failure behind nested solving; callers select `Auto` when recovery is
  desired.
- [x] Add unknown vector `[ln(n_0), ..., ln(n_{m-1}), theta_T]` with a smooth
  bounded temperature transform. Do not clip temperature inside residual
  evaluation. `equilibrium_ph_formulation` now owns deterministic unknown and
  row layouts plus a tested logistic transform that rejects floating-point
  saturation instead of clipping.
- [x] Add residual rows `[reaction equilibrium, element balances, enthalpy]`
  while preserving the exact P,T row ordering and one explicit energy scale.
  The prepared formulation is solver-agnostic and is published through the
  fixed-declared-phase facade only after the common candidate gate accepts it.
- [x] Reuse the analytic P,T log-mole Jacobian as the upper-left block and add
  the temperature column plus enthalpy row analytically.
- [x] Require component-aligned `Cp(T)` for the production monolithic path.
  Any finite-difference fallback must be a typed thermochemistry capability
  with provenance, not an implicit numerical trick in the solver. The prepared
  monolithic formulation rejects a missing component capability explicitly.
- [x] Keep resolved `G0(T)` synchronized with the coefficient interval used by
  the canonical `P,T` bridge. `ResolvedThermochemistry` now owns a private
  per-phase, temperature-keyed Gibbs snapshot: it selects all phase
  coefficients and constructs the phase `G0` functions once per temperature,
  then shares them across component residual rows. This fixed the live
  `P,T -> H -> monolithic P,H` mismatch without mutating resolved data or JSON
  libraries.
- [x] Extend the activity contract with `d ln(a) / dT`. The current ideal-gas
  and ideal-solution models have zero derivative at fixed pressure; future
  non-ideal models must provide or explicitly decline this capability.
- [x] Compare the full analytic P,H Jacobian against central finite
  differences block by block, including the temperature chain rule.
  The initial test covers reaction, element, and enthalpy rows against all
  log-mole and temperature columns; live thermochemistry fixtures remain part
  of backend integration.
- [x] Run one monolithic solve per fixed active set and carry its accepted
  composition and temperature into the next phase-control transition.
  - [x] `equilibrium_ph_monolithic::PreparedMonolithicPhRunner` now performs
    the coupled fixed-active-set solve through the common legacy backend
    cascade, reuses the shared P,T candidate gate, and rejects the existing
    P,T-only RST symbolic payload explicitly instead of faking compatibility.
  - [x] The fixed declared public facade bridges the accepted monolithic
    snapshot into the ordinary immutable phase-aware result with the solved
    temperature, matching lookup report, and backend evidence. No mutable
    `EquilibriumLogMoles` compatibility object participates in that path.
  - [x] Let `PreparedPhaseControlRunner` invoke this fixed-set runner after
    active-set transitions. This is implemented through the existing injected
    candidate callback, not a second orchestration facade. A monolithic
    candidate may probe a trace-seeded wider mask when the current reduced
    branch has no enthalpy root; the runner records that mask expansion as an
    ordinary activation transition before final publication.
    - [x] Extract the shared active-set lifecycle behind an injected fixed-set
      candidate callback. The existing P,T runner uses that path now, while
      hysteresis, boundary recovery, cycle detection, transition reports, and
      publication remain owned by `PreparedPhaseControlRunner`.
    - [x] Make phase-stability evaluation candidate-local: it now consumes the
      accepted candidate's Gibbs capabilities and conditions rather than the
      runner construction temperature. This is required before a P,H candidate
      with solved temperature can enter the same loop.
    - [x] Add the reduced monolithic P,H candidate adapter. It projects
      both `PreparedEquilibriumProblem` and `ResolvedThermochemistry` by the
      active mask, solve `[ln(n_active), theta_T]`, scatter back to the full
      layout, and provides full-layout `G0(T_solution)` to stability checks.
- [x] Return one result type for monolithic, nested, and fallback paths, with
  the selected path, fallback reason, backend attempts, residual blocks,
  energy error, conservation evidence, and timing stated explicitly.
  - [x] The existing `FixedPressureEnthalpySolution` now exposes the actual
    `PhSolvePath` through its common report; nested and monolithic live paths
    are both asserted, including `MonolithicPhaseControl`. A classified
    fallback reason and `Auto` policy remain open.

### P6.0 Freeze the physical and dimensional contract

- [x] Introduce a typed top-level equilibrium constraint instead of boolean
  flags. Keep fixed pressure/reference pressure in the existing condition
  types and distinguish:
  - prescribed `P,T`;
  - prescribed `P,H`;
  - the `P,H` temperature seed or bracket hint;
  - the accepted equilibrium temperature;
  - the prescribed total extensive enthalpy.
  `equilibrium_constraints::EquilibriumConstraint` now provides this boundary
  and is consumed by the typed P,H facade.
- [x] Make total enthalpy in joules the canonical engine input. The total is
  tied to the supplied closed-system inventory; any GUI molar, mass-specific,
  or mixture-normalized input is an explicit conversion before entering the
  engine. `TotalEnthalpyJoules` is now the typed `PH` boundary, while the raw
  constructor remains as a transition convenience.
- [x] Document the initial supported physical scope: ideal-gas phases and pure
  condensed phases at fixed pressure, using the current standard-state
  thermochemistry. Do not imply support for non-ideal solution enthalpy,
  pressure-dependent real-fluid enthalpy, excess enthalpy, kinetic/potential
  energy, or unmodelled heat/work terms. The public module documentation now
  makes this boundary explicit before any GUI exposes the workflow.
- [x] Define the enthalpy reference convention. All component enthalpies in one
  solve must come from thermochemically compatible records and reference
  states; `H_target` must use the same convention. The public contract states
  this explicitly and `ResolvedThermochemistry` already pins every component
  to its selected record and provenance.
- [x] State precisely when `P,H` represents an adiabatic isobaric calculation:
  no unaccounted heat transfer, no shaft/electrical work, and no omitted
  kinetic or potential energy. This is now part of the public module contract,
  not a hidden assumption of the scalar solver.
- [x] Do not require source compatibility at any cost. Preserve `P,T`
  numerical behavior and provide a clear migration path, but remove or
  deprecate an old constructor if retaining it would duplicate the canonical
  condition model. Compatibility constructors that duplicate the canonical
  resolved-thermochemistry boundary are now explicitly deprecated; the narrow
  `from_resolved_thermochemistry(...)` builder remains the production entry.

### P6.1 Extend the resolved thermochemistry capability

- [x] Add the first resolved-system `dH(T)` bridge. `EnthalpyModel::from_resolved_system`
  aligns capabilities to `SystemLayout` and rebuilds temperature-dependent
  property readers in private `SubsData` copies.
- [x] Add a typed per-component thermochemistry capability aligned exactly to
  `SystemLayout`, carrying at least:
  - standard Gibbs free energy `g0(T)`;
  - molar enthalpy `h(T)`;
  - optional heat capacity `cp(T)`;
  - valid temperature interval;
  - source library, record key, phase/state evidence, and lookup provenance.
  `ResolvedThermochemistry` carries this bundle and the accepted P,H result
  retains it for diagnostics.
- [x] Reuse the existing `SubsData` `dH(T)` and `Cp(T)` calculators and
  closures. Do not reconstruct enthalpy by differentiating Gibbs free energy,
  and do not copy NASA/NIST formula implementations into `ChemEquilibrium`.
  The resolved P,H bundle is format-agnostic: it consumes the common
  `ThermoCalculator` capability API after lookup and never branches on a
  library name or polynomial representation.
- [x] Build Gibbs, enthalpy, and optional heat-capacity capabilities from the
  same selected thermochemical record. Reject mixed-record property bundles
  unless an explicit, validated reconciliation policy is introduced later.
- [x] Preserve the current read-only repository transaction: resolving or
  evaluating `P,H` thermochemistry must never mutate the JSON libraries.
- [x] Compute the admissible temperature domain as the intersection of all
  selected records' valid intervals. Reject an empty intersection before any
  nonlinear or scalar solver starts.
- [ ] Replace the current single interval/envelope representation with an
  exact ordered set of disjoint valid segments if a supported NASA/NIST
  record can contain gaps. Until then, runtime property evaluation must still
  reject a temperature that falls inside an envelope gap; the envelope is not
  permission to extrapolate across missing coefficients.
- [x] Treat enthalpy as mandatory for the first production outer-temperature
  workflow. Keep `Cp` optional there; it becomes mandatory only for an
  analytic derivative/Newton acceleration or the monolithic formulation.
- [x] Expose a pure additive derivative-ready primitive for
  `dH/dT|n = sum_i n_i Cp_i`. Its result is explicitly named a partial
  derivative and does not pretend to include the implicit equilibrium term
  `sum_i h_i d(n_i)/dT`.
- [ ] Define a future extension point for phase-model enthalpy contributions.
  The current additive contract
  `H = sum_i n_i h_i(T)` must fail explicitly for any model requiring an
  unimplemented excess/mixing enthalpy term.

### P6.2 Add the canonical outer-temperature workflow

- [x] Add the first closure-backed bracketed workflow. The
  `equilibrium_ph_workflow` module provides `EnthalpyModel`, the
  `ResolvedThermochemistry` bundle, validated scalar controls, trial/report
  types, and `solve_resolved_ph` over the canonical `solve_resolved_pt`
  facade. The bundle-backed constructor is now the resolved-data boundary;
  continuation, budgets, and phase-transition policy remain open below.
- [x] Promote the prototype into the final narrow typed `P,H` request above
  the resolved-phase facade. The production request must own:
  - immutable `ResolvedPhaseSystem`;
  - typed initial composition;
  - pressure and reference pressure;
  - target total enthalpy;
  - validated temperature bounds and optional seed;
  - fixed-`P,T` solve options and phase-control policy;
  - scalar temperature-solver policy, budgets, cancellation, and timing mode.
  `ResolvedPhaseEnthalpyRequest` now captures an owned immutable
  `ResolvedPhaseSystem` snapshot, typed composition/constraint/bounds, inner
  solver and phase policy, scalar policy, budgets, cancellation, and timing.
- [x] Define a pure trial evaluation:
  `F(T) = H(solution_of_P_T(T), T) - H_target`.
  A trial is usable only after the inner `P,T` candidate passes the existing
  backend-independent acceptance gate. `PhTrialEvaluator` now creates the
  local immutable outcome while the outer workflow alone owns budgets,
  progress, continuation, and final publication.
- [x] Use a proven safeguarded bracketed scalar method (Brent-style or
  equivalent bisection/interpolation hybrid) as the production default.
  The current default accepts a secant proposal only inside a guarded interior
  portion of the valid sign bracket and otherwise uses bisection. Trial
  reports retain `LowerBound`, `UpperBound`, `Seed`, `Interpolation`, or
  `Bisection`, so the scalar path is auditable. No unguarded Newton step is
  enabled.
- [x] Add deterministic bracket construction and validation:
  - explicit user bracket takes precedence;
  - the required `P,H` temperature seed is evaluated as an additional trial
    strictly inside, and never instead of, user bounds; it deterministically
    narrows an already valid bracket when possible;
  - [x] endpoint and midpoint trial failures retain the zero-based trial
    index, temperature, typed inner cause, source chain, and retryability
    classification instead of becoming an unqualified scalar-solver error;
  - an unbracketed or unreachable target returns a typed error; a seed that
    exposes two sign-change intervals returns `enthalpy_multiple_brackets`
    instead of selecting a root by incidental evaluation order.
- [x] Reuse expensive immutable preparation across fixed-declared-phase trial
  temperatures: `PreparedPhaseEquilibriumTemplate` now retains the resolved
  layout, element matrix, reaction basis, provenance, numeric closures, and
  optional RST symbolic problem. The P,H workflow creates it lazily at the
  first real trial, then retargets only temperature-dependent
  thermochemistry while preserving deterministic multi-start recovery. The
  report publishes one build and per-trial reuse evidence. Bounded phase
  control remains intentionally separate because its active-set lifecycle
  cannot be reused across the non-monotone temperature order of a scalar
  bracket without changing the physical contract.
- [x] Use the nearest accepted fixed-declared-phase composition as an
  additional typed log-mole continuation seed without publishing trial state.
  Bounded phase control intentionally keeps its ordinary inventory seed until
  active-set-aware continuation is implemented; no active set is assumed to
  remain unchanged across a phase transition.
- [x] Add one global scalar temperature-evaluation budget and enforce it in
  the bracket helper itself, so endpoint and interior evaluations share the
  same limit rather than receiving a fresh budget per trial.
- [x] Forward `EquilibriumExecutionControl` into the outer P,H workflow and
  every inner P,T request, with explicit temperature-trial progress stages.
- [x] Enforce an optional wall-time budget at the scalar evaluation gate. The
  budget is global to the outer operation; inner cooperative cancellation is
  still required for interruption during a long backend call. The deadline is
  now checked after evaluator return as well, so a trial that finishes after
  the deadline cannot be accepted retroactively.
- [x] Make the complete `P,H` solve transactional. Failed trial solves,
  rejected phase transitions, cancellation, exhausted budgets, or a failed
  final acceptance check must not publish a partial result or mutate the
  caller's resolved system.
- [x] Define nested resource budgets explicitly: scalar evaluations, total
  inner backend attempts, total nonlinear iterations/evaluations, wall time,
  and phase-control transitions. A per-trial budget must not accidentally
  multiply into an unbounded outer budget. The scalar evaluation budget is
  global, and `PhTemperatureSolveOptions` now adds optional global
  `max_inner_backend_attempts`, `max_inner_nonlinear_iterations`, and
  `max_phase_control_transitions` limits.
  The outer transaction accounts for every started inner backend attempt and
  every reported nonlinear iteration before publishing its trial result,
  including all work belonging to every continuation multi-start seed. Phase
  transitions are counted from the immutable solution report and are subject
  to the same global outer budget.
  The request owns an immutable resolved snapshot, while accepted trials and
  the final result remain local until the scalar bracket and enthalpy contract
  succeed.
- [x] Integrate cooperative cancellation and progress events at both levels:
  bracket preparation, temperature trial start/accept/reject, inner backend
  attempt, phase transition, and final publication. Trial-level cancellation
  and progress are wired, including a typed rejected-trial event for inner
  P,T or enthalpy-evaluation failures. The P,H report now also emits inner
  backend start/finish, accepted transition, and publication start/finish
  events; cancellation is checked at each publication boundary.

### P6.3 Handle phase transitions and non-smooth enthalpy honestly

- [ ] Do not assume `H_eq(P,T)` is globally smooth or monotone. Phase
  appearance/disappearance can create derivative discontinuities, and
  coexistence can create very flat or discontinuous-looking numerical
  branches.
  - [x] The safeguarded outer solver now has an explicit
  `PhMonotonicityPolicy`. Its default rejects an observed reversal in the
  sampled enthalpy branch with a typed error; the legacy sign-bracket
  behavior is available only through an explicit compatibility policy and
  is recorded in the immutable solve report. This is sampled evidence, not
  a claim of global monotonicity.
  - [x] Every accepted temperature trial now retains the exact phase-id to
    lifecycle-status snapshot accepted by its inner P,T solve. Consumers can
    therefore compare adjacent active sets instead of inferring a transition
    only from aggregate counters.
- [ ] Define how the scalar solver treats active-set changes inside a bracket.
  Retain transition evidence per trial and prevent interpolation steps from
  treating values from incompatible failed branches as a smooth derivative.
  - [x] The bounded phase-control P,H path now forces bisection. Safeguarded
    secant interpolation remains available only for fixed declared phases,
    where the phase layout is invariant. This prevents a numerical slope from
    being inferred across two potentially different active sets.
  - [x] Every published trial now exposes a typed preparation reason:
    `FixedFormulationInitial`, `FixedFormulationReused`, or
    `BoundedPhaseControlIsolated`. The latter makes the absence of active-set
    reuse explicit rather than looking like a missed performance metric.
- [ ] Detect and report multiple brackets/multiple roots when sampled evidence
  reveals them. The first production policy must be deterministic (for
  example, nearest valid root to the requested seed), not dependent on hash or
  backend iteration order.
  - [x] The mandatory interior seed already detects the first concrete
    multiple-bracket evidence: equal-sign endpoints with an opposite-sign seed
    are rejected as `enthalpy_multiple_brackets`. Broader scan-based root
    enumeration remains intentionally deferred until phase-boundary semantics
    are finalized.
- [ ] Define boundary behavior for latent-heat/coexistence cases. If the
  requested enthalpy is achieved by changing phase fractions at nearly fixed
  temperature, the accepted result still must satisfy composition,
  complementarity, conservation, and enthalpy tolerances.
  - [x] Add a release-only real gas/ice `P,T -> H -> P,H` boundary fixture.
    The monolithic active-set candidate exhausts its backend cascade, but
    explicit `Auto` deterministically recovers through bounded nested `P,T`
    phase control, retains `AllBackendsFailed` as fallback evidence, activates
    the solid phase, preserves balances/enthalpy, and leaves local libraries
    unchanged.
  - [ ] Define a true phase-fraction/coexistence formulation only after the
    physical model and complementarity contract are specified. The successful
    ice recovery is evidence for one point solution, not proof that every
    latent-heat plateau is represented correctly.
- [x] Reuse the bounded phase-control hysteresis and cycle/budget protections.
  The bounded monolithic P,H adapter now runs through
  `PreparedPhaseControlRunner`, so hysteresis, transition limits, cycle
  detection, rollback, and typed lifecycle reports are shared with the P,T
  path. A scalar P,H trial cannot hide an inner phase-control cycle as a
  generic scalar-function failure; the accepted water fixture also verifies
  trace-seeded activation of a zero-inventory liquid phase.

### P6.4 Publish an immutable `P,H` solution and diagnostics

- [x] Return an immutable result that contains the accepted temperature,
  pressure, component moles, phase states, target/calculated total enthalpy,
  raw enthalpy error in joules, relative/scaled error, and the normal
  fixed-`P,T` acceptance reports. Bundle-backed results also retain the
  component-aligned thermochemistry provenance and common domain.
- [x] Add a validated enthalpy acceptance contract with both absolute and
  scale-aware relative tolerances. The scalar solver accepts only when
  `|H-H_target| <= max(abs_tol, relative_tol * scale)`; near-zero targets
  therefore retain a meaningful absolute floor.
- [x] Define one positive finite enthalpy scale from the immutable request and
  initial state, not from arbitrary solver iterates. Record the scale and both
  tolerance components in the report so acceptance is reproducible.
- [x] Retain complete nested evidence:
  - bracket endpoints and accepted root;
  - every temperature trial and its status;
  - inner backend attempts and solver metrics;
  - phase transitions and active-set changes;
  - rebuild/reuse reasons;
  - thermochemistry provenance;
  - timing by repository, property refresh, preparation, nonlinear solve,
    phase control, enthalpy evaluation, scalar orchestration, and
    postprocessing.
  Every accepted trial now retains its complete inner backend cascade,
  optional multi-start comparison, phase-control lifecycle, and final
  acceptance snapshots alongside compact counters and timing. Rejected scalar
  trials remain typed errors rather than partial result rows. Each published
  trial now also carries an immutable `PhTrialTimingReport` separating trial
  wall time, nested `P,T` time, and additive enthalpy evaluation. The timing
  report is opt-in and covered by the P,H unit evidence test. Each inner
  `PhaseTransitionRecord` now also carries a measured control-pass duration;
  live activation/deactivation stories verify that this timing is published.
  - [x] Publish a point-level `formulation_build` duration for typed T-range
    points. It accounts for reduced formulation and RST symbolic construction
    performed during that point and is covered by the live ice-transition
    range story.
   - [x] Preserve separate stopwatch intervals for every individual formulation
     cache entry when one point creates multiple active-set formulations. The
     public point duration remains an aggregate, while each accepted range
     point now carries a deterministic active-mask/cache-entry timing snapshot.
  - [x] Per-trial and aggregate backend-attempt, nonlinear-iteration, and
    phase-transition counters are now published together with nested timing;
    trial/evidence cardinality is validated before report publication.
  - [x] Fixed-declared-phase P,H reports now state the number of immutable
  formulation builds and retarget reuses. Bounded phase control deliberately
  reports no such reuse because rejected scalar trials cannot safely mutate
  or seed its active-set lifecycle.
  - [x] `Auto` now preserves both typed failure trees when its monolithic
    attempt and nested recovery both fail. Successful recovery retains the
    classified fallback reason in the immutable P,H report; failed recovery
    returns `PhAutoFallbackFailed { monolithic, nested }` rather than erasing
    the primary active-set evidence.
  - [x] Keep the default production cascade free of unsolicited console
    output. RustedSciThe Powell Dogleg remains an explicit benchmark/backend
    option, but is excluded from the default cascade because version 0.4.12
    prints internal `beta` diagnostics and rejected the large real-data case.
- [x] Keep independent `K_eq` validation scoped to applicable inner
  fixed-phase chemical-equilibrium candidates. It does not independently
  validate the outer enthalpy root and must not be presented as doing so. The
  P,H module documentation now explicitly separates this observational inner
  evidence from enthalpy-root acceptance.
- [x] Add a compact stable public facade and prelude exports only after the
  request/result/error contracts are validated. Keep raw prepared state and
  scalar orchestration internals crate-private. `solve_resolved_ph`, its typed
  request/result/report types, and constrained solver controls are exported
  from `ChemEquilibrium::prelude`; prepared runner/template internals remain
  crate-private.

### P6.5 Verification and release evidence

- [x] Add a typed fixed-pressure, fixed-enthalpy target-range facade with a
  strictly monotone ascending/descending enthalpy grid, continuation from the
  previous accepted physical composition and temperature, per-point timing,
  phase-transition evidence, and transactional publication. `PhRangeRequest`
  is the canonical engine boundary; it does not publish a partially solved
  batch after a point failure.
- [x] Cover the continuation contract with unit tests for grid/error/report
  semantics and an ignored live NASA H/O story that checks ascending and
  descending targets, accepted temperatures, provenance, conservation, and
  the seed hand-off between adjacent points.
- [x] Reuse prepared monolithic P,H formulation state across fixed-phase target
  points. Reaction basis, element totals, phase projection, analytic row
  structure, and RST symbolic expressions are built once; accepted
  `[log-moles, T]` state plus target/scale parameters are retargeted between
  points. The range report exposes formulation build/reuse counters.
- [x] Extend the same prepared-state contract to nested P,H batches without
  conflating it with the monolithic formulation. The nested route now shares
  only its fixed-P,T prepared template across target batches; every target
  retains an independent scalar bracket, while bounded phase-control remains
  isolated whenever active-set lifecycle can change. The live continuation
  story asserts one template build across the nested target range.
- [x] Add a live backend matrix for the fixed-phase monolithic P,H target
  range with strict single-backend policies and per-point timing. The ignored
  matrix keeps backend failure visible instead of hiding it behind a cascade.
- [x] Extend the live route matrix with nested/Auto batches, a real
  phase-transition range, and rollback after an unreachable target; these are
  separate lifecycle stories and must not be implied by the fixed-phase
  matrix. The ignored real-data route matrix now prints all three rows and
  keeps the rollback point index typed. Its release evidence is recorded in
  `STORY_TESTS.md`.
- [x] Add a large real-data nested P,H continuation story: twenty exact
  element-limited C/H/O species, nine interior target enthalpies, accepted
  temperature/composition hand-off, one fixed-P,T template build, conservation
  and enthalpy acceptance at every point.
- [x] Add a denser real water/ice Auto P,H story: five target enthalpies from
  250 to 270 K, route-dependent monolithic/nested acceptance, retained
  fallback reason, phase-transition evidence, and transactional output.
- [x] Record release evidence for the nested/Auto route matrix, the large
  twenty-species nested story, and the dense water/ice Auto story in
  `STORY_TESTS.md`. Keep these large live tests out of the default suite.
- [x] Implement and document the strict fixed-phase monolithic P,H target-range
  backend matrix. It remains separate because backend failures are deliberate
  characterization evidence rather than a route-lifecycle failure; the release
  command and aligned backend table are recorded in `STORY_TESTS.md` for the
  target-machine evidence refresh.
  - [x] The matrix now emits one aligned row per backend with point count,
    formulation builds/reuses, solve and wall timing, maximum residual,
    maximum element-balance error, and the complete typed failure reason.
    A debug characterization currently shows \`legacy_tr\` completing all
    three points while the other rows expose backend-specific acceptance or
    iteration-limit failures; this is evidence, not a production-default
    decision.

- [x] Keep every existing fixed-`P,T` regression green and add direct parity
  between the resolve/build pipeline request and the resolved `P,T` facade.
  The parity story fixes one backend policy, compares component amounts,
  residual quality, and layout fingerprint, and therefore checks orchestration
  equivalence without conflating it with backend-cascade selection.
- [x] Add the primary inverse story:
  solve a stable `P,T` problem at `T*`, calculate `H*`, then solve the same
  inventory at `P,H*` from a different seed/bracket and recover temperature,
  composition, conservation, and enthalpy within explicit tolerances.
  The ignored live NASA H2/O2/H2O story now exercises this contract through
  `ResolvedThermochemistry::from_resolved_system` and `solve_resolved_ph`.
- [x] Cover synthetic analytic fixtures where `h(T)` and the expected
  temperature/root are known exactly. Include a non-reacting sensible-heat
  case so failures in chemistry cannot mask errors in the energy equation;
  the same unit layer also checks inventory scaling invariance.
- [ ] Add real offline cross-format gas and gas-plus-pure-condensed stories
  with pinned lookup policy and provenance. Snapshot all canonical JSON files
  before and after the tests. The solver contract is library/polynomial-format
  agnostic, so this must include at least one stable non-NASA fixture rather
  than treating NASA as the production format.
  - [x] Ignored release stories now cover NASA-gas reactive P,T -> H -> P,H
    inversion and NASA-gas plus NASA-condensed bounded water inversion. Both
    retain provenance, verify the P,H reuse/phase evidence appropriate to
    their solve mode, and compare the canonical library files byte-for-byte.
  - [ ] Add a pinned local NIST P,T -> H -> P,H inverse fixture once the
    read-only repository actually contains a stable NIST record set. The
    current `all_keys_substance.json` advertises historical NIST addresses,
    but the canonical local repository does not contain the corresponding
    H2/O2/H2O records, so an offline test would only pin stale-index failure.
    Repair the library/index consistency first; do not turn this into a
    network-dependent regression. This fixture is evidence for the
    format-agnostic capability boundary, not an alternate solver path. See
    **F1** for the required atomic data-release and full offline matrix.
  - [x] Make the fallback contract explicit before adding that fixture:
    local data is authoritative; an enabled canonical fallback queries NIST
    only for the requested gas/liquid/solid state; an unspecified state fails
    closed; and parser/network failures are preserved as failures. The old
    unconstrained gas fallback remains compatibility-only and deprecated.
  - [x] Add an ignored online NIST state-matrix smoke test for the parser
    boundary. It may validate gas/liquid/solid payload availability and finite
    Shomate data, but it must not be used as equilibrium evidence or mutate
    local JSON. Offline NIST parity remains blocked until real payloads are
    checked into the repository.
    - [x] The current H2O run reports complete gas/liquid payloads and an
      incomplete solid navigation payload without Cp intervals. The latter is
      recorded as unavailable solid thermochemistry, never as a phase fallback.
  - [x] Add a read-only catalog consistency report before repairing that
    fixture. It canonicalizes library aliases and reports missing payload,
    orphan payload, and duplicate index pairs without modifying JSON. The
    live release diagnostic currently reports `5519` indexed pairs, `5472`
    payload pairs, and `47` indexed records without payload; this is now an
    explicit data-release blocker rather than a hidden solver failure.
- [x] Add phase-transition stories for water vapor/liquid/ice and another
  physically credible condensed system. Check latent-heat/coexistence
  behavior, hysteresis, rollback, and final complementarity.
   - [x] Real P,T water liquid/ice appearance and disappearance stories are
     covered; the ignored P,H gas/ice inverse story verifies bounded nested
     recovery after a rejected monolithic candidate. The release rerun
     recovered `T=250 K`, `solid_moles=4.998093e-1`, one transition, and zero
     scaled enthalpy error through `Auto`.
   - [x] The ignored real phase-transition release matrix now reports water
     gas/ice, water gas/liquid, hot-water disappearance, and graphite
     appearance/disappearance in one table. Every row checks finite
     non-negative amounts, conservation, complementarity, transition evidence,
     and byte-for-byte JSON immutability. The release run passed for all five
     scenarios, with per-case solve times from roughly 1.1 to 5.0 ms.
   - [x] The real ice temperature-range story now retains one deterministic
     build-duration snapshot per prepared active-set cache entry, instead of
     exposing only the aggregate point duration.
- [x] Test invalid and unreachable inputs: non-finite enthalpy, invalid
  pressure/temperature bounds, missing enthalpy capability, incompatible
  record intervals, unsupported phase enthalpy model, unbracketed target,
  inner all-backends-failed, cancellation, and exhausted global budget.
  - [x] The P,H unit layer already covers non-finite targets, invalid scalar
    budgets, missing/misaligned `Cp`, wrong enthalpy component counts,
    incompatible bundle domains, unbracketed targets, observed multiple
    brackets, endpoint/midpoint inner failures, cancellation before and after
    inner start, and global evaluation/wall-time budgets. These paths return
    typed errors and never publish a partial `FixedPressureEnthalpySolution`.
  - [ ] Add the remaining resolved-data cases only when their physical source
    exists: an actual record without `dH`, an unsupported non-additive phase
    model, and a repaired local NIST fixture. Do not fabricate malformed
    production records merely to make this checklist look complete.
  - [x] The ignored real-data validation matrix covers non-finite target,
    reversed bounds, non-positive pressure, unreachable target rollback, and
    JSON immutability.
- [x] Add metamorphic inventory-scaling tests. Multiplying every initial mole
  and total target enthalpy by the same factor should preserve equilibrium
  temperature and mole fractions while scaling extensive amounts. The
  analytic outer-solver fixture now pins this contract independently of
  chemistry.
- [x] Add a backend/cascade matrix for the inner `P,T` solves and a release
  characterization over small, medium, and large real systems. Report outer
  evaluations separately from measured inner nonlinear solve time; do not
  select a production default from wall time alone.
  - [x] The ignored live H2/O2/H2O `P,H` backend matrix runs each RST and
    legacy backend under a strict `Single` policy against one independently
    constructed `H*`. Its rounded table reports outer trials, total/wall and
    inner nonlinear timings, formulation build/reuse counts, inner work,
    final temperature, energy error, residual, balance, and an explicit error
    row for failed methods. It snapshots canonical JSON libraries before and
  after the complete matrix.
   - [x] The strict monolithic P,H target-range matrix now emits aligned
     per-backend timing, attempts, accepted-backend, enthalpy, residual,
     conservation, failure-point, failure-kind, and full typed-error evidence.
     Release execution remains the final characterization step.
   - [x] `live_nested_ph_inner_pt_cascade_release_matrix` now covers 5, 20,
     and 100 real local NASA candidates over three target enthalpies. It reports
     outer evaluations, inner backend attempts, inner solve time, formulation
     reuse, conservation, and immutable-library evidence. The debug baseline
     passes; the release command is recorded in `STORY_TESTS.md`.
- [x] Verify deterministic ascending/descending enthalpy sweeps in the typed
  P,H batch facade. Reuse accepted neighboring solutions transactionally,
  while keeping the semantic distinction explicit: one P,H target has one
  solved temperature.

### P6.6 Evaluate and promote the monolithic formulation

- [x] Only after the outer-temperature workflow was accepted, prototype a
  coupled unknown layout with log-moles and bounded/logarithmic temperature.
  The canonical reaction and element rows remain embedded in the prepared
  problem; no second multiplier-based formulation was introduced. Keep the
  layout typed; do not spread arithmetic such as `n_species + n_elements`
  through residual code.
- [x] Reuse the existing chemical and balance blocks and append exactly one
  scaled enthalpy equation. Do not create a second copy of the `P,T`
  formulation.
- [x] Derive and test the full temperature column of every chemical residual,
  including `g0(T)/(RT)`, pressure/activity terms, and any phase-model
  temperature dependency. A zero temperature column is invalid.
- [x] For additive ideal enthalpy, verify the analytic energy derivatives:
  `dH/dln(n_i) = n_i h_i(T)` and
  `dH/dln(T) = T sum_i n_i cp_i(T)`.
  Extend these formulas before enabling any excess-enthalpy model.
- [ ] If analytic `Cp` or activity derivatives are unavailable, isolate a
  scale-aware finite-difference fallback and compare the complete Jacobian
  against central differences away from phase boundaries. Do not replace the
  established analytic composition Jacobian wholesale.
- [x] Bound temperature through a validated transform or safeguarded step
  policy over the common thermochemistry interval; `ln(T)` alone enforces
  positivity but does not enforce the upper/lower data bounds.
- [ ] Compare monolithic and outer-temperature formulations on convergence
  basin, phase transitions, backend fallback behavior, residual quality,
  thermochemistry evaluations, and release timing. The real reactive inverse
  story now compares both paths and establishes monolithic as the default for
  resolved requests. A real water/liquid phase-activation story now directly
  compares explicit monolithic and nested routes with the same `P,H` target,
  phase lifecycle evidence, accepted temperature, composition, residual, and
  conservation. The real gas/ice story now additionally proves the hard-route
  contract: direct monolithic returns `AllBackendsFailed`, explicit nested
  accepts, and `Auto` publishes the nested-equivalent state while retaining the
  same fallback reason. The broader phase-transition/release matrix remains
  before declaring promotion complete.

### P6.7 GUI integration gate

- [x] The fixed-phase ideal/pure-condensed engine facade, immutable reports,
  typed errors, cancellation, budgets, and release backend evidence are now
  sufficient to begin GUI integration within the documented P6.0 scope. GUI
  must not advertise non-ideal phases, latent-heat coexistence, or an
  enthalpy-temperature range sweep as implemented capabilities.
- [x] Expose total enthalpy in explicit joules, temperature seed and bracket,
  inner solver policy, progress, cancellation, immutable energy evidence, and
  route-dependent diagnostics through typed GUI controls. Nested routes expose
  scalar trial rows; monolithic routes expose backend attempts, accepted
  backend, acceptance, and phase-control evidence without fabricated trials.
  The GUI selects the canonical monolithic route when the resolved
  thermochemistry supports one exact symbolic interval; a bracket crossing a
  native coefficient switch remains on the nested numeric reference route.
- [x] Do not offer a temperature-range control in `P,H` mode. A future batch
  operation should sweep target enthalpy and report the solved temperature at
  each point.
- [x] Add the basic GUI stories: valid local P,H solve, immutable energy and
  outer-trial diagnostics, lifecycle staleness, cancellation boundary, and
  document roundtrip.
- [x] Add GUI stories for unreachable target, inner fallback success,
  all-backends-failed, cancellation/rollback during an active P,H worker, and
  a real phase-transition publication using stable local fixtures.
  - [x] The ignored local GUI story now covers an unreachable finite P,H target:
    the worker fails visibly, publishes no partial snapshot, and leaves the
    canonical JSON files unchanged.
  - [x] Existing local GUI stories cover fallback acceptance, transactional
    range failure, cancellation/stale-result rejection, and real water/ice
    phase publication.
  - [x] All-backends-failed is a separate ignored local fixture and asserts
    the typed engine error plus absence of a partial snapshot; it does not
    manufacture a generic worker failure solely to satisfy the checklist.

## P6.8 Release-oriented monolithic P,H architecture hardening

The external architecture review is accepted as a follow-up plan, with the
following corrections to its scope. The coupled P,H mathematics and the
nested reference route already exist; this work must preserve the P,T facade,
the phase-control boundary, and the existing typed result contract. The GUI
now covers the basic P,H editor, worker, immutable energy/outer diagnostics,
trial table, cancellation boundary, and document lifecycle. The remaining
items below are engine architecture and release evidence, not a reason to
rebuild the GUI.

### P6.8.1 Dependency map and module boundaries

- [x] Produce a checked dependency map for `equilibrium_ph_workflow.rs`,
  `equilibrium_ph_formulation.rs`, `equilibrium_ph_monolithic.rs`, and the
  nested implementation. Record every shared type, thermochemistry adapter,
  `GibbsFn` construction site, `Cp` capability check, option/report type, and
  fallback decision before moving code. The map above is now synchronized
  with the extracted nested/thermochemistry modules and the shared typed
  nonlinear-system adapter.
- [x] Move `ResolvedThermochemistry`, its provenance, molar property
  function types, interval intersection, and Gibbs/enthalpy/Cp evaluation
  into a lower-level thermochemistry module. Formulation and runners may
  depend on this module; it must not depend on workflow orchestration,
  publication, phase lifecycle, GUI progress, or scalar bracket policy. The
  workflow retains only compatibility re-exports and the nested enthalpy
  adapter.
- [x] Extract the nested scalar-temperature algorithm and its trial types,
  bracket policy, monotonicity policy, continuation, budgets, and
  nested-specific report into `equilibrium_ph_nested.rs`. Keep the workflow
  as a facade, mode selector, fallback coordinator, and common publication
  boundary.
  - [x] Extract the stateless safeguarded interpolation, wall-time guard,
    sampled-monotonicity helpers, and the scalar bracket loop with unit tests.
    The module also owns `PhTemperatureTrial`, phase-state snapshots, timing,
    and inner-solve evidence. The workflow re-exports those types, adapts
    scalar records, attaches phase evidence, and publishes the final common
    result.
- [x] Do not introduce a second thermochemistry bundle or a second public
  result merely to perform this extraction. Reuse the existing immutable
  bundle and expose accessors where the current result already owns the data.
  The nested module owns scalar trial/evidence snapshots only; the common
  `FixedPressureEnthalpySolution` remains the single published result.

### P6.8.2 Typed thermochemistry evaluation

- [x] Remove every `unwrap_or(f64::NAN)` conversion from the new monolithic
  P,H path. An out-of-domain evaluation, missing coefficient, database error,
  or property evaluation failure must retain its typed cause through
  preparation/residual evaluation and must not be reported as an anonymous
  numerical backend breakdown. The remaining infallible `GibbsFn` crossing is
  an explicit, finite snapshot adapter at the accepted temperature.
- [x] Choose one explicit fallible boundary for P,H: either a fallible
  prepared-property/evaluation context or residual helpers that consume
  already evaluated Gibbs/enthalpy/Cp vectors. Do not change the established
  P,T closure API wholesale solely to introduce a new closure alias. The
  monolithic residual already consumes the typed evaluated context; phase
  control remains a documented compatibility boundary until its callback is
  made fallible.
- [ ] Replace the snapshot compatibility adapter with a fallible stability
  evaluator once phase-control no longer requires `GibbsFn = Fn(T) -> f64`.
- [x] Move mandatory `Cp(T)` capability validation to formulation preparation.
  Validate component count/order, common temperature domain, finite initial
  Gibbs/enthalpy/Cp values, and the selected capability policy before a
  nonlinear backend starts. Missing Cp is now an `InvalidProblem` before
  backend execution; a deliberately enabled finite-difference policy remains
  future work.
- [ ] If finite-difference Cp is supported, add explicit analytic/central,
  forward, and backward source diagnostics, boundary-aware step selection,
  and a Jacobian comparison test. Do not enforce an unconditional positive-Cp
  rule without first defining its physical scope.

### P6.8.3 Separate solve policies and diagnostics

- [x] Separate common physical/acceptance/execution options from monolithic
  nonlinear options and nested scalar/bracket options. Nested-only controls
  must not affect monolithic solves, and monolithic-only controls must not be
  silently applied to nested recovery.
  - [x] `PhAcceptanceOptions`, `PhMonolithicOptions`, and `PhNestedOptions`
    are typed route contracts; `PhTemperatureSolveOptions` exposes validated
    projections into them while retaining the compatibility builder surface.
  - [x] The fixed-active-set monolithic runner now receives only an immutable
    enthalpy acceptance contract; scalar bracket limits, phase-transition
    budgets, and progress controls remain in the workflow. The nested scalar
    engine now receives its own smaller `NestedBracketOptions` contract.
- [x] Separate monolithic diagnostics from nested diagnostics. A monolithic
  Newton iteration must not be represented as an outer temperature trial;
  nested reports must retain trials, bracket evidence, inner P,T work, and
  scalar residuals. Shared counters may remain in a common immutable result
  only when their semantics are identical. Monolithic backend and phase
  lifecycle evidence now lives in `PhMonolithicEvidence`, while `trials` is
  populated only by the nested scalar route.
- [x] Isolate the compatibility constructor
  `ResolvedPhaseEnthalpyRequest::new(...)` from the canonical API. It is now
  deprecated, explicitly documented as nested/reference-only, and the
  production `from_resolved_thermochemistry(...)` builder no longer calls it.
- [x] Audit `Auto`: retry nested only for classified retryable numerical
  failures; never retry for invalid input, missing capability, layout/data
  errors, cancellation, unsupported phase models, or invariant failures.
  Preserve the original monolithic error tree and test the classification
  table explicitly. `AllBackendsFailed` is retryable only when its non-empty
  trace consists solely of started numerical failures; rejected candidates,
  skipped attempts, and empty traces fail fast.
- [x] Keep phase-control lifecycle outside the monolithic residual and keep
  the nested route as a reference/fallback. Do not move hysteresis or phase
  transitions into the coupled equation itself without a new physical
  complementarity contract.

### P6.8.4 Backend and publication contracts

- [x] Remove the monolithic runner's dependence on `EquilibriumLogMoles` as a
  mutable orchestration host. Reuse neutral residual/backend helpers only;
  preparation, iteration, phase lifecycle, and result publication must use
  immutable prepared state and local candidates. The backend cascade is now a
  stateless adapter function; the old associated method is only a legacy
  wrapper.
- [x] Define a domain-independent prepared nonlinear-system contract that can
  represent both fixed P,T and P,H residual/Jacobian dimensions. RST must
  receive residual/Jacobian capabilities without knowing which coordinate is
  temperature or which row is enthalpy. The adapter now carries
  `PreparedNonlinearSystem` with dimension, residual, Jacobian, feasibility,
  and optional backend payload without attaching thermodynamic meaning to any
  coordinate.
  - [x] Add a dedicated symbolic monolithic `P,H` payload for RST. It uses the
    exact resolved component order, substitutes the bounded `T(theta)`
    expression into phase-qualified `G0(T)` and `H(T)`, appends the scaled
    enthalpy row, and gives RST the complete `(N + 1)` system. The analytic
    formulation remains an independent legacy/fallback route; no fixed-P,T
    payload or fake zero temperature column is reused.
  - [x] Keep symbolic thermochemistry optional in `ResolvedThermochemistry`.
    Numeric resolved data remain valid for legacy and nested P,H solving.
    Real symbolic capabilities are materialised lazily from private
    phase-local `SubsData` copies for the requested P,H bounds and are
    accepted only when every component remains inside one native NASA/NIST
    coefficient interval. Selecting RST without `G0/H` expressions, or across
    a native polynomial boundary, fails with a typed unsupported-capability
    error instead of fitting/extrapolating silently or changing lookup
    semantics.
- [x] Retain one physical P,H publication gate for both paths: finite and
  in-range temperature, finite log-moles/reconstructed moles, reaction and
  element residuals, enthalpy error, conservation, phase stability, and
  active-set consistency. Diagnostics may differ, acceptance meaning may not.
- [x] Refresh the immutable build report at the accepted monolithic
  temperature before publication. A fixed-active P,H solve is prepared at its
  seed temperature but may converge elsewhere; the report and accepted
  solution must carry the same `(P,T)` conditions. The live GUI monolithic
  story covers this publication contract.
- [x] Add a typed capability error that names the requested incompatible
  backend and the supported alternatives. The monolithic RST guard now uses
  `UnsupportedBackendCapability` both for absent symbolic capabilities and
  for a requested symbolic range crossing a native coefficient boundary. The
  error is non-retryable and therefore cannot be hidden by fallback.

### P6.8.5 Architecture tests and release characterization

- [x] Add architecture tests proving that monolithic mode performs no outer
  nested temperature trials, nested mode reports real scalar trials, and
  `Auto` records monolithic-then-nested fallback only after a retryable
  failure. Include an invalid-input case proving nested is not started.
  - [x] The real-data P,H story now asserts that monolithic reports no scalar
    trials and exposes `PhMonolithicEvidence`; the unit contract also records
    that invalid monolithic input emits no `TemperatureTrialStarted` event.
- [x] Add a route-dependent GUI diagnostics snapshot. The monolithic GUI view
  now presents backend-attempt summaries, accepted backend, acceptance rows,
  phase-control rows, and inner timing; it does not render a zero-valued
  temperature-trial section. Nested routes retain their real trial table.
- [ ] Extend the analytic Jacobian matrix with missing-Cp, typed property
  error, analytic/finite-difference Cp, zero activity-temperature derivative,
  non-zero activity-temperature derivative, lower/upper interval boundary,
  scaling, and every P,H block-column/row comparison.
  - [x] Add a synthetic RST P,H regression that verifies the symbolic
    `(ln(n), theta_T)` system converges through the normal candidate gate and
    never reuses the fixed-P,T symbolic payload.
  - [x] Add an ignored real-data NASA-gas parity case for RST-monolithic P,H
    against the independent analytic/legacy route. The 1900..2900 K fixture
    stays within the common native interval and compares accepted temperature,
    composition, conservation, and enthalpy rather than backend iteration
    counts. A companion real-data guard proves that 900..1100 K is rejected
    for the RST symbolic path when it crosses a coefficient switch, while the
    JSON repository remains unchanged.
   - [x] Add an ignored real-data NASA-gas Jacobian central-difference matrix.
     The fixture uses five locally resolved C/H/O records and points immediately
     inside both common temperature boundaries plus an interior point. It
     compares every analytic P,H block entry and verifies that the JSON source
     files remain unchanged. The test intentionally keeps ideal-model activity
     temperature derivatives at their physical zero value.
    - [x] Extend that real matrix across four inventory scales from `1e-6` to
      `1e3`, while retaining both native-interval boundaries and the complete
      residual/Jacobian block comparison.
    - [x] Release rerun confirms the same real Jacobian matrix: 432 entries,
      four inventory scales, three boundary/interior temperatures, and a
      six-dimensional coupled formulation.
  - [ ] Add the corresponding ignored NIST P,H parity fixture after a stable
    locally bundled NIST record set is selected. Do not weaken the exact
    single-native-interval contract to make this fixture pass.
- [x] Add one release characterization over the same real problem for
  monolithic and nested paths. Report backend attempts, residual/Jacobian and
  thermochemistry evaluations, active-set work, temperature trials, total
  inner P,T solves, reuse, and wall time. Keep performance ratios out of unit
  test assertions; use an ignored release benchmark instead. The ignored live
  P,H story now prints a table for nested, analytic monolithic, symbolic
  RST-LM monolithic, and Auto with backend metrics plus
  thermochemistry-preparation timing.
- [x] Add the remaining GUI stories when their engine fixtures exist:
  unreachable target, P,H inner fallback, all-backends-failed,
  cancellation/rollback, and a real phase-transition publication. The stories
  now exist in `equilibrium_gui_tests`; catalog-dependent cases remain
  explicitly ignored and preserve JSON immutability checks. The mixed
  gas/condensed phase-publication story selects Legacy NR explicitly because
  the symbolic monolithic route requires one native coefficient interval per
  component. The basic valid P,H, diagnostics, lifecycle, and
  document-roundtrip stories are covered.

### P6.8.6 Documentation and explicit non-goals

- [x] Update module rustdoc after extraction: workflow as facade, monolithic
  as fixed-active-set coupled solve, nested as reference/recovery algorithm,
  thermochemistry as capability/provenance layer, and formulation as the
  complete `[F_rxn, F_elem, F_H]` system with its block Jacobian. The module
  headers now describe these ownership boundaries and explicitly distinguish
  outer scalar trials from monolithic coupled evidence.
- [x] Define and implement typed P,H temperature-range batch solving with
  accepted-state continuation, transactional publication, and per-point
  evidence. Export/import, latent-heat coexistence modeling, and non-ideal
  phase physics remain outside this technical hardening pass until their
  separate contracts are defined.
- [x] Re-run the full P,T regression suite after every logical extraction;
  the current full `ChemEquilibrium` library run passes after the option and
  prepared-system extractions, and no P,T compatibility workaround became a
  hidden dependency of the monolithic P,H path.

## P7 - Replace phase stability with reacting-system TPD

This is a correctness replacement for the current ideal phase-control path,
not a future non-ideal extension. The contracts in `task.md` and
`arch_contract.md` were compared with the implementation on 24.08.2026.

### P7 review verdict

- [x] **Confirmed:** the present implementation already evaluates canonical
  chemical potentials as `mu_i = g_i^0(T) + R*T*ln(a_i)` through
  `PhaseActivityModel::log_activity`. Preserve and reuse that boundary; do not
  introduce a second activity implementation inside phase stability.
- [x] **Confirmed:** the current stability criterion is not general. It
  special-cases one-component condensed phases, rejects a multicomponent
  `IdealSolution`, and rejects a rank-deficient active elemental assemblage.
  Those restrictions are implementation artifacts rather than the required
  reacting-system stability contract.
- [x] **Confirmed:** `PhaseStabilityModel`, `driving_force`, and the current
  per-species `PhaseSeedPolicy` encode the old `q = 1` workflow. In particular,
  phase activation does not seed the candidate with the composition that
  minimizes its tangent-plane distance.
- [x] **Confirmed:** stability mathematics currently lives in
  `equilibrium_workflows.rs`, while both the prepared `P,T` runner and the
  monolithic `P,H` route consume its limited report and seed semantics.
- [x] **Retain:** active-set projection, bounded outer iterations, hysteresis,
  cycle detection, rollback, and transactional publication remain valid
  orchestration mechanisms. Rewire them to typed TPD evidence rather than
  reimplementing them.
- [x] **Qualified:** a one-component phase may retain an analytical fast path,
  but only as an internal optimization proven equivalent to the general
  `q = 1` TPD problem. It must not remain a separate public physical criterion
  or a hidden fallback.
- [x] **Qualified:** the solver core has `PhaseActivityModel::IdealSolution`,
  while the resolved phase domain currently exposes only `IdealGas` and
  `PureCondensed`. The migration must introduce an explicit ideal-solution
  capability at the phase-domain boundary before the top-level facade claims
  support for a multicomponent condensed solution. This does not imply any
  non-ideal excess-Gibbs model.

### P7.0 Freeze the contract and dependency graph (P0)

- [x] Record every producer and consumer of `compute_phase_stability_reports`,
  `PhaseStabilityReport`, `PhaseStabilityModel`, `ElementPotentialReport`,
  `driving_force`, `dg_create`, `dg_keep`, `seed_activated_phase`,
  `PhaseSeedPolicy`, `PhaseManager`, `PreparedActiveSetCandidate`,
  `PreparedPhaseControlRunner`, monolithic `P,H` active-set restarts,
  complementarity reports, GUI snapshots, presentation/export helpers, and
  reproducibility capsules.
- [x] Freeze the mathematical units and signs: `minimum_tpd`, `dg_create`, and
  `dg_keep` are molar Gibbs quantities in `J/mol`; a sufficiently negative
  minimum requests activation, zero is coexistence within tolerance, and a
  positive minimum is stable against creation of that candidate phase.
- [x] Define phase capability explicitly. `Excluded` means intentionally not
  evaluated and must remain distinct from an unsupported activity model, an
  elementally infeasible candidate, and a failed TPD minimization.
- [x] Add `IdealSolution` to the resolved phase-model contract with validation
  of its allowed physical states and ordered components. Keep
  `PureCondensed` as a one-component convenience constructor or normalize it
  to the same ideal-solution activity implementation; do not infer ideal
  mixing merely from a liquid/solid phase name.
  - [x] **Corrective audit closed:** `PhaseModel::IdealSolution` is now an
    explicit condensed ideal-mixing declaration, while `PureCondensed` rejects
    every component count other than one. `MultiphaseEquilibriumLayout` accepts
    the new semantic model and the bridge maps it explicitly to
    `PhaseActivityModel::IdealSolution`.
- [x] Treat this as an intentional API replacement. Do not add adapters that
  fabricate old `driving_force` or `PureCondensedSpecies` semantics from a new
  TPD result. Migrate all retained consumers in one dependency-directed pass.
  - [x] The audit leaves exactly one mutable owner: `PreparedPhaseControlRunner`.
    The pure `equilibrium_phase_stability` service receives immutable state and
    returns evidence only; workflows classify that evidence and no production
    caller retains obsolete driving-force or per-species seed semantics.
  - [x] `PhaseStabilityStatus` distinguishes policy exclusion,
    `FixedGasAssemblage`, absence of an independent reference assemblage, and
    a genuine evaluated TPD result. Typed errors cover failed construction or
    infeasible minimization rather than being encoded as a numeric status.

### P7.0a Corrective semantic-boundary audit (P0)

The post-refactor architecture review confirms that the constrained TPD kernel
is the right foundation and must not be replaced by the historical pure-phase
criterion. The unfinished work is at the semantic construction boundary and
in lifecycle evidence around it.

- [x] **Accept the diagnostic's central finding.** The current chain is
  inconsistent: `PhaseModel` has only `IdealGas` and `PureCondensed`, while
  `PhaseActivityModel` already has `IdealSolution`; `activity_model_for`
  currently maps `PureCondensed` to that numerical law, and the production
  layout enforces one component for the semantic condensed model. Unit-level
  multicomponent TPD tests therefore do not prove production support.
- [x] Introduce an explicit semantic `PhaseModel::IdealSolution` and validate
  the physical states for which the present ideal-mixing contract is defined.
  Keep `PureCondensed` semantically distinct and one-component, even though its
  numerical activity law is the `q = 1` special case of `IdealSolution`.
- [x] Migrate the complete definition chain without a compatibility adapter:
  `PhaseSpec -> ResolvedPhaseSystem -> MultiphaseEquilibriumLayout ->
  PhaseEquilibriumMetadata -> EquilibriumPhaseDescriptor ->
  PhaseActivityModel`. Update canonical ordering/fingerprints, factories,
  candidate selection, public prelude types, and all exhaustive matches.
- [x] Add a production bridge regression that constructs a multicomponent
  `IdealSolution` through `SubstanceSystemFactory` and the local repository,
  then through `ResolvedPhaseSystem` and `PhaseEquilibriumBuildRequest`, and
  proves that every ordered component,
  thermochemistry capability, element row, phase range, and
  `PhaseActivityModel::IdealSolution` reaches the resulting problem bundle.
  No branch on this path may require `components.len() == 1` for
  `IdealSolution`. The local NASA-condensed fixture proves two thermochemistry
  records, lookup provenance, element rows, descriptors, activity laws, and
  the full component range reach `PhaseEquilibriumProblemBundle` unchanged.
- [ ] Add an end-to-end fixed-`P,T` lifecycle story in which an initially
  inactive multicomponent ideal solution has negative accepted `minimum_tpd`,
  is seeded with the reported `incipient_composition`, is re-solved, and is
  published only after residual, conservation, feasibility, KKT, and
  complementarity validation pass.
- [x] Redesign `SupportedPhaseModelPolicy` after its dependency audit. The
  retained versioned capability contract is now `IdealPhaseModelsV1`; its
  documentation names ideal gas, pure condensed, multicomponent ideal
  solution, and canonical constrained TPD. The obsolete
  `FixedPressureTemperatureV1` spelling has no compatibility alias.
- [x] Extend the GUI phase-model enum and validation after the engine
  semantic bridge is canonical. The GUI must serialize/construct
  `IdealSolution` explicitly rather than disguising it as `PureCondensed`.
- [x] Audit reference-assemblage construction for inactive candidates, active
  retention, a single active phase, multiple active condensed phases, and a
  rank-deficient active elemental space. The audit retained the canonical rule:
  every candidate is compared only against other active phases, a single active
  condensed phase reports `ActiveWithoutReferenceAssemblage`, policy-excluded
  phases are never assigned a numeric TPD, and rank-deficient reference rows
  use the SVD geometry rather than an accidental full-rank assumption. Focused
  orchestration regressions cover exclusion, the missing-reference state,
  multiple active condensed phases, and a rank-one H/O reference space.
- [ ] Add the remaining monolithic `P,H` recovery matrix with typed evidence:
  reduced active-set success without a probe; reduced failure followed by a
  successful bounded all-active probe; a probed but TPD-stable phase remaining
  inactive; a TPD-unstable phase activated from `x*` rather than the neutral
  probe composition; stability evaluated at accepted `T*`; and complete
  rollback when both reduced and recovery solves fail. The recovery probe is
  numerical branch discovery, never an activation criterion.
  - [x] The runner now treats all-active probe occupancy as non-publishable.
    It restores inactive probe phases to trace and permits one deterministic
    transition only when `minimum_tpd < dg_create`, seeded from `x*`; after an
    activation it re-solves the new fixed set before acceptance. Real water
    stories prove that strict `Monolithic` and `Auto` reject a TPD-stable
    liquid probe rather than manufacturing appearance. Deterministic
    reduced-success/recovery/restart/rollback fixtures now protect the
    lifecycle itself; the remaining observability gap is typed probe counters
    in report evidence.
  - [x] **Regression review, retained evidence:** real strict/`Auto` water
    stories prove that a wider numerical probe with `minimum_tpd >= dg_create`
    is rejected rather than published; the dense water/ice `P,H` range proves
    every published TPD report carries the accepted point temperature; generic
    runner rollback preserves accepted continuation; direct multicomponent
    deactivation uses total phase amount; and capsules preserve
    `CanonicalTpdV1` plus named `IdealSolution` semantics.
  - [x] **Deterministic lifecycle harness:** crate-private injected-candidate
    stories now force, independently: (a) a reduced fixed-set success with no
    wider candidate; (b) a numerically wider TPD-stable candidate rejection;
    (c) a TPD-unstable multicomponent `IdealSolution` activation from a
    visibly non-uniform `x*`, followed by a new fixed-set solve; and (d)
    rollback when that restart fails. They assert masks, transition reason,
    `minimum_tpd`, `incipient_composition`, seed ratios, accepted candidate
    temperature, validation evidence, and unchanged accepted continuation.
    The real P,H stories remain responsible for enthalpy residuals; neither
    suite fixes optimizer iteration counts or SVD-call order.
  - [ ] Add compact typed numerical-probe evidence to the P,H lifecycle report
    for the preceding stories: attempted wider masks, outcome class
    (unneeded/failed/TPD-rejected/TPD-activated), and restart outcome. It must
    describe numerical branch discovery without presenting a probe candidate
    as a physical phase transition or adding a public mutable API.
- [x] Strengthen multicomponent deactivation regressions. A phase with one
  trace-small component but substantial `sum_i n_i` must remain active; a
  phase may be deactivated only from its total amount plus the accepted
  retention TPD criterion. The direct regression covers both a trace-small
  component with substantial phase total and a truly vanishing total phase.
- [x] Replace stale production wording such as `inactive pure phase` and
  `Detect phases that must be created (Delta G < 0)` with candidate-phase,
  minimum-TPD, and phase-stability terminology. Retain legitimate historical
  `q = 1` explanations and schema tests that deliberately reject the old
  `driving_force` representation.
- [x] Reassess reproducibility schema compatibility after the serialized
  semantic `PhaseModel` change. `EquilibriumPhaseSpecSnapshot` stores explicit
  model names rather than enum ordinals, so prior `IdealGas` and
  `PureCondensed` capsules retain their exact meaning while current capsules
  can record `IdealSolution`. `CanonicalTpdV1` and capsule schema v2 remain
  correct; the regression round-trips the new symbolic model value. GUI
  documents have the same named-enum property, so their schema remains v1 and
  older documents need no migration.
- [x] **Already satisfied:** canonical state construction, element-potential
  reconstruction, exact elemental feasibility, constrained ideal TPD,
  TPD-derived activation composition, transactional publication, typed
  stability evidence, presentation-only diagnostics, and rejection of stale
  `driving_force` fields are retained. No parallel pure-phase fast path or
  replacement minimizer is required by this audit.
- [x] **Deferred physical decision:** `FixedGasAssemblage` remains explicit
  until the project decides whether multiple gas declarations represent one
  shared mixture, separate compartments, or competing phases. Generic gas
  TPD cannot be chosen as an engineering cleanup because activity
  normalization depends on that physical contract.
- [x] **Out of scope:** non-ideal activity coefficients, fugacity/EOS phase
  split, and global non-convex multi-start TPD remain in the future-physics
  section. The semantic ideal-solution bridge must not imply these models.

### P7.1 Extract the stability domain and canonical inputs (P0)

- [x] Create `equilibrium_phase_stability.rs`. It owns chemical-potential
  snapshots, element-potential reconstruction, elemental-direction
  feasibility, candidate TPD problems, minimization, and typed stability
  reports. It must not mutate `PhaseSet`, choose transitions, publish a solve,
  or own cycle/hysteresis policy.
- [x] Keep `equilibrium_workflows.rs` and `PreparedPhaseControlRunner`
  responsible only for orchestration: fixed-active-set solve, request a pure
  stability analysis, classify the returned evidence, seed/reproject, detect
  cycles/budgets, and transactionally publish an accepted result.
- [x] Build one immutable accepted-state input containing solver-order moles,
  phase totals, canonical `g_i^0`, canonical activities/chemical potentials,
  conditions, element matrix, active/candidate masks, and layout identity.
  Validate dimensions, finiteness, positivity requirements, and snapshot
  alignment once at construction.
- [x] Reuse `PhaseActivityModel` for `ln(a_i)` and expose only the additional
  composition-level operation needed by TPD. The residual, Jacobian, and
  stability paths must share the same standard-state and pressure convention.
  - [x] First extraction: `equilibrium_phase_stability.rs` now owns an
    immutable canonical-state snapshot. The retained pure-phase workflow
    obtains `mu_i` from that snapshot rather than recomputing its own activity
    convention.

### P7.2 Reconstruct element potentials without a full-rank shortcut (P0)

- [x] Solve `A_active * lambda ~= mu_active` with a documented SVD tolerance.
  Publish `lambda`, numerical rank, singular-value/tolerance evidence, maximum
  absolute residual, scaled residual, and the solver species used in the fit.
- [x] Validate representability of accepted active chemical potentials with
  an absolute-plus-relative residual contract. A poor fit is a typed physical
  validation failure; it must not silently produce phase decisions.
- [x] Do not reject `rank < element_count`. Element potentials may be
  non-unique; stability must remain invariant to null-space-equivalent choices
  of `lambda` for an elementally feasible candidate direction.
- [x] Construct the active elemental range/null space once per unchanged
  active set and cache it in the prepared phase-control state. Rebuild it only
  after an active-set/layout change and report the rebuild/reuse decision.
  - [x] SVD fit now records rank, tolerance, residual, and reference species,
    accepts a rank-deficient valid fit, and the phase-control workflow no
    longer reintroduces a `rank < element_count` veto after reconstruction.
    A null-space-equivalence regression protects the resulting TPD value.
  - [x] `PreparedPhaseControlRunner` now owns a runner-scoped cache of immutable
    reference-assemblage geometry. One cache entry retains the SVD fit factors
    and full `A_active^T A_active` null space; unchanged reference species
    reuse it for elemental-potential fitting, feasibility, and constrained TPD
    constraints. `TemperatureRangeSolveReport` publishes entry/build/reuse
    counters, while a unit regression proves one build plus one reuse leaves
    the fit unchanged. Temperature, activities, chemical potentials, and TPD
    values remain uncached candidate-state data.

### P7.3 Enforce candidate elemental-direction feasibility (P0)

- [x] For candidate composition `x`, evaluate
  `c_beta(x) = A_beta^T * x` and require it to lie in
  `Range(A_active^T)`. Implement the check through an SVD/null-space basis with
  explicit absolute-plus-relative tolerance, not through a full-rank guard.
- [x] Include simplex constraints `x_i >= 0` and `sum(x_i) = 1` in the
  candidate problem itself. Never normalize an unconstrained answer after the
  optimizer and call it a constrained minimum.
- [x] Return a typed infeasible-candidate result when the admissible set is
  empty. Do not evaluate or classify a TPD value for a physically infeasible
  trial composition.
- [x] Report the candidate elemental composition and feasibility residual for
  every evaluated phase so complementarity and release evidence can audit the
  decision.
  - [x] The canonical ideal TPD path constructs elemental null-space rows
    `B*x = 0`, solves them as part of the candidate simplex problem, and then
    independently verifies `A_active^T*y ~= A_candidate^T*x`. No full-rank
    shortcut remains in the connected phase-control path.
  - [x] `CanonicalPhaseState` is the only accepted-state ingress for the
    connected TPD workflow and delegates activity evaluation to the shared
    `PhaseActivityModel`. `ElementalFeasibilityReport` retains candidate
    elemental totals, residual, and tolerance in every evaluated
    `PhaseStabilityReport`; an empty feasible simplex returns a typed error
    before any TPD classification or transition publication.

### P7.4 Implement the general ideal TPD minimum (P0)

- [x] Define the single canonical objective
  `TPD_beta(x) = sum_i x_i * (mu_i^beta(x,T,P) - a_i*lambda)` over the
  admissible simplex. Store the minimum and all thresholds in `J/mol`; no
  arbitrary trace amount belongs to the TPD definition.
- [x] Cover `IdealSolution` candidates with the convex analytical/specialized
  structure: both the full-rank closed-form case and the rank-deficient
  constrained case use deterministic minimization with independently checked
  KKT/feasibility evidence.
- [x] Return both `minimum_tpd` and normalized `incipient_composition = argmin
  TPD` in declared phase-component order. Validate finite objective values,
  finite/nonnegative composition, unit sum, elemental feasibility, and
  deterministic tie handling.
- [x] Handle boundary optima without evaluating `ln(0)`. Zero composition is
  valid in the mathematical minimizer; a positive numerical floor is allowed
  only when constructing log-mole restart coordinates.
- [x] If minimization, feasibility, or validation fails, return a typed error
  with phase/context evidence and abort the transition transaction. Do not
  fall back to uniform composition, a random start, the old pure-phase
  formula, or a penalty-only answer.
- [x] Keep the TPD problem interface ready for a future non-convex activity
  model and multi-start/global policy, but do not implement or imply a fake
  non-ideal model in this stage.
  - [x] Interior pass: `IdealTpdProblem` now supplies an analytical softmax
    minimum for full-rank inputs and a deterministic dual-Newton minimizer for
    rank-deficient interior simplexes. It validates independent elemental
    constraints, feasibility, and KKT residuals; the connected workflow uses
    it for both pure and multicomponent `IdealSolution` phases. Unit tests
    cover `q = 1`, the analytic multicomponent result, rank-deficient interior
    feasibility, and invariance to non-unique lambda representatives.
  - [x] Boundary pass: a deterministic phase-I linear program identifies the
    relative-interior feasible support, then dual Newton runs only on that
    face. Exact zero components remain in `incipient_composition`; the
    positive floor is introduced only by the separate log-mole seed builder.
    Regression tests cover a boundary-only candidate and seed-total recovery.
  - [ ] **Physical decision gate: separately declared `IdealGas` phases.** The
    present lifecycle intentionally treats them as one fixed gas assemblage,
    reported as `FixedGasAssemblage`, rather than as competing candidate
    phases. Before implementing generic inactive-`IdealGas` TPD coverage,
    decide whether multiple gas declarations mean one shared mixture, separate
    compartments, or genuinely competing gas phases. The answer determines
    activity normalization and is therefore not an engineering cleanup.

### P7.5 Replace reports, hysteresis, and activation seeding (P1)

- [x] Replace the old physical-model/driving-force report with a typed result
  containing phase identity, active state, evaluation status,
  `minimum_tpd`, `incipient_composition`, element-potential evidence,
  feasibility evidence, minimizer diagnostics, conditions, and layout
  identity. Excluded/not-evaluated phases carry no fabricated numeric value.
  - [x] `PhaseStabilityReport` now has `PhaseStabilityStatus` plus
    `minimum_tpd`; obsolete `PhaseStabilityModel` and the misleading
    `driving_force` field are removed. Transition and complementarity evidence
    now use the same TPD terminology. Skipped phases report an explicit status
    and no numeric minimum.
  - [x] Evaluated reports retain conditions, canonical ordered layout evidence,
    elemental-feasibility residuals, and constrained-minimizer/KKT diagnostics.
    Skipped reports retain the same identity and conditions but no fabricated
    numerical stability evidence.
- [x] Apply `dg_create` and `dg_keep` only after the mathematical minimum has
  been accepted. Hysteresis is a transition policy over `minimum_tpd`; it must
  not alter the objective or define phase stability.
  - [x] `PhaseManager` classifies only accepted `minimum_tpd` values;
    thresholding is downstream of TPD construction and minimization.
- [x] Replace per-species seed semantics with a total phase-seed amount policy.
  For an activated phase construct `n_i = n_phase_seed*x_i_star`, apply a
  strictly positive floor only for log coordinates, and renormalize so the
  seeded component sum equals the requested total seed.
  - [x] `PhaseSeedPolicy` and the uniform `seed_activated_phase` helper have
    been removed after their last production caller migrated. The monolithic
    all-active recovery uses the same total-seed API with a documented neutral
    numerical probe; physical phase appearance always uses the TPD minimizer.
- [x] Preserve transactional mass handling: failed seed construction,
  projection, reduced solve, candidate validation, or complementarity must
  leave the last accepted composition, temperature, active set, caches, and
  published report unchanged.
  - [x] Phase transition records now retain the chosen `incipient_composition`
    independently of the floored `restart_seed`; both prepared and mutable
    orchestration paths seed an activated phase through
    `seed_activated_phase_with_composition`.
  - [x] A failed prepared-runner lifecycle now restores the earlier accepted
    continuation seed and phase set. The injected-failure regression proves a
    rejected solve cannot consume continuation intended for the next point.
    Projection/formulation/TPD geometry caches are immutable-layout
    memoization, never accepted-solution publication; every numerical use is
    retargeted before solving.
- [x] Update boundary recovery, transition records, cycle fingerprints, and
  complementarity checks to consume `minimum_tpd` and the minimizer-derived
  seed. Record the chosen seed composition in transition evidence.
  - [x] Transition and complementarity paths use `minimum_tpd`; activated
    transitions retain the exact TPD minimizer independently of the floored
    restart seed. Cycle fingerprints deliberately remain phase-set based until
    a composition-sensitive cycle policy has a physical contract.

### P7.6 Integrate the same criterion into P,T and P,H (P1)

- [x] Wire the prepared fixed-`P,T` phase-control runner to the extracted TPD
  service without rebuilding invariant matrices or activity metadata on every
  outer iteration when the active set is unchanged.
- [x] Run monolithic `P,H` stability analysis at the accepted solution
  temperature `T*`, never at the initial temperature seed. Do not create a
  separate enthalpy-specific phase criterion.
- [x] After a `P,H` transition, seed the phase with `x*`, retain accepted `T*`
  as the continuation temperature, rebuild the active-set formulation, and
  solve again against the same target enthalpy transactionally.
- [x] Propagate typed TPD failure and transition evidence through point,
  temperature-range, and enthalpy-range reports. A failed point must not seed
  the next continuation point.
  - [x] The prepared runner is the sole lifecycle owner for fixed `P,T` and
    bounded monolithic `P,H`: it evaluates TPD from the candidate's accepted
    conditions, seeds transitions from `x*`, and publishes the same immutable
    acceptance report through each point/range solution. `PhRangeRequest`
    advances composition and temperature only after a point is accepted; its
    real water/ice release story retains transition and phase-control evidence
    on each accepted point.

### P7.7 Migrate diagnostics and public consumers (P1)

- [x] Update `MultiphaseAcceptanceReport`, phase-control traces, presentation
  rows, GUI route-specific snapshots, tables, plotting boundaries,
  reproducibility capsules, and story-test output to use TPD terminology and
  retain minimizer/feasibility evidence where appropriate.
  - [x] Acceptance summary rows and the GUI phase-evidence table now retain
    per-phase status, accepted `minimum_tpd`, elemental-feasibility residual,
    and KKT residual when a phase was evaluated. Aggregate complementarity is
    displayed separately, so it cannot hide a weak individual minimization.
- [x] Update stable facade/prelude exports to expose the new immutable report
  types while keeping minimizer implementation details crate-private.
  - [x] `ChemEquilibrium::prelude` now exposes the immutable phase-control,
    transition, TPD, feasibility, minimizer, validation, and typed-index
    evidence. The canonical-state and minimizer implementation remains
    crate-private; external consumers can inspect accepted reports but cannot
    invoke a second mutable orchestration path.
- [x] Add a schema/version decision for serialized reproducibility evidence;
  reject stale old stability fields explicitly rather than interpreting them
  under new semantics.
  - [x] Reproducibility capsule schema v2 records `CanonicalTpdV1` explicitly.
    Its typed JSON loader rejects every older schema and recursively rejects
    obsolete `driving_force`/`PhaseStabilityModel` fields rather than silently
    reinterpreting them as TPD evidence.
- [x] Remove console-only diagnostics from minimization. Timings, iterations,
  termination, KKT/feasibility residuals, and cache reuse belong in typed
  reports and deterministic story tables.
  - [x] The canonical TPD/minimizer, prepared phase-control runner, and
    fixed-`P,T` range route contain no direct `print!`/`println!` calls.
    Runtime evidence is carried by immutable reports and story-table formatters;
    legacy backend debug logging is outside this minimization contract.

### P7.8 Required test matrix (P0/P1)

- [x] Pure `q = 1` regression: general TPD equals
  `g_beta^0(T) - a_beta*lambda`, including sign, units, and minimizer `[1]`.
- [x] Stable and unstable multicomponent ideal-solution fixtures: verify the
  sign and value of `minimum_tpd`, normalized `x*`, and transition decision.
- [x] Full-rank analytical ideal minimum: compare both `minimum_tpd` and `x*`
  against the closed-form reference, not only the final active set.
- [x] Seed regression: component seed ratios match `x*`, every log-mole seed is
  positive, and physical seeded moles sum to the requested phase amount after
  floor-and-renormalize handling.
- [x] Boundary optimum: permit exact zero components in `x*`, avoid `ln(0)`,
  and prove that the log-coordinate floor does not change the reported TPD
  minimizer.
- [x] Rank-deficient active assemblage: an elementally feasible candidate is
  accepted, and TPD is invariant under admissible null-space changes to
  `lambda` within tolerance.
- [x] Elementally infeasible candidate: return explicit infeasibility with no
  TPD classification, activation, or published partial state.
- [x] Hysteresis story: activation below `dg_create` is retained when the next
  accepted minimum lies between `dg_create` and `dg_keep`; test the reverse
  direction separately.
- [x] P,H regression: stability uses accepted `T*`; activation retains `T*` as
  the next continuation seed and preserves the target enthalpy contract.
- [x] Property/metamorphic tests: normalization, non-negativity, finite TPD,
  feasibility residual, element-potential residual, permutation invariance,
  inventory scaling, mass conservation, rollback, cycle/budget termination,
  and cache reuse across repeated active sets. The mathematical unit matrix
  includes direct component-permutation invariance; live reactive-gas stories
  cover inventory scaling, while phase-control regressions cover rollback,
  cycle/budget termination, and runner-scoped geometry-cache reuse.
- [ ] Re-run all retained `P,T`, `P,H`, fixed-phase, phase-transition, GUI, and
  release story matrices. Add offline real-data water/ice, water/liquid, and
  carbon/graphite evidence without mutating JSON libraries; synthetic tests
  remain as precise mathematical fixtures.
  - [x] The current release ledger already records the real TPD
    water/ice, water/liquid, hot-water, graphite, and hot-carbon transition
    matrix in `STORY_TESTS.md`, including conservation, complementarity,
    transition count, timing, and JSON immutability. This is retained evidence,
    not a substitute for the final whole-suite release campaign.
  - [x] The paired release ice T-range story records `tpd_geometry=1/1/3`:
    one runner-scoped active-element geometry build and three reuses while the
    real range crosses a phase transition. It simultaneously verifies the
    projection/prepared/RST cache layout, continuation, conservation, and JSON
    immutability.

### P7.9 Cleanup and definition of done (P1)

- [x] Delete obsolete `PhaseStabilityModel` variants/enum, old
  `driving_force` data paths, the multicomponent capability veto,
  full-element-rank rejection, uniform per-species activation seeding, and
  tests asserting the old limitation. The remaining `ValidationNotApplicable`
  uses are generic typed applicability errors, not a solution-model veto.
  `PureCondensed` retains its intentional one-component semantic constraint;
  the explicit `IdealSolution` model carries multicomponent mixtures. The
  legacy indexed numerical API now labels an undifferentiated ideal-solution
  activity law as `IdealSolution`, never falsely as `PureCondensed`.
- [ ] Repair the stale/corrupted phase-stability subsection in
  `ARCHITECTURE_RU.md`. It still describes scalar `Delta G` driving force and
  says multicomponent solutions are unsupported; rewrite it around constrained
  `minimum_tpd`, `incipient_composition`, and `FixedGasAssemblage` without
  altering the executable stability kernel.
- [x] Remove duplicated pure-phase stability helpers and any internal
  compatibility translation. An optional `q = 1` fast path must be tested
  against the general formulation and remain invisible to consumers.
- [ ] Mark P7 complete only when the canonical dataflow is:
  `accepted fixed-set state -> canonical mu -> element potentials -> feasible
  TPD minimization -> minimum_tpd + x* -> hysteresis + seed -> rebuilt
  active-set solve -> transactional acceptance` for both `P,T` and `P,H`.
  The kernel and current pure-condensed production path satisfy this flow, but
  P7 remains open until a multicomponent semantic `IdealSolution` traverses
  the same production facade and the bounded `P,H` recovery matrix proves that
  a neutral all-active probe cannot leak into physical activation semantics.
- [x] Document the remaining boundary honestly: P7 completes general phase
  stability for the currently supported ideal activity models. Non-ideal
  excess-Gibbs/fugacity models and global multi-start policy remain in F2.

## Deferred Foundation: Local NIST Fixture and New Physics

The current production path is deliberately scoped to ideal, fixed-pressure
equilibrium with local data. The items in this section are not small cleanup
tasks. They require an explicit data release or a new physical contract, and
must not be quietly approximated by the existing NASA-gas fixtures.

### F1. Local NIST fixture is a data-release blocker

- [ ] **Create a stable, read-only local NIST thermochemistry fixture.**
  The NIST parser is an online acquisition boundary, while production and
  regression equilibrium calculations must remain offline and deterministic.
  The current local index advertises historical NIST addresses but has no
  matching payload records; that is not an acceptable cross-format fixture.

  Required deliverables:
  - choose a small physically meaningful record set with pinned canonical
    names, library keys, physical states, temperature intervals, elemental
    composition, and NIST/Shomate payloads;
  - write the payload, address/index, and elemental-composition entries as one
    reviewed atomic data update, plus a manifest/fingerprint that identifies
    the data release used by the tests;
  - make the repository consistency report show every selected NIST index
    entry paired with an exact payload, without orphan aliases;
  - add offline `P,T` point/range and `P,H` inverse/range stories that retain
    NIST provenance, verify conservation and acceptance, and prove that no
    network call or JSON mutation occurred;
  - add at least one mixed local NASA/NIST story only after the individual
    NIST records are stable, with explicit per-component provenance and a
    documented temperature-domain intersection.

  Do not use live NIST access, a parser mock, or an incomplete historical
  index as a substitute. This fixture proves the format-agnostic contract; it
  is not a second solver implementation.

### F2. Future physical extensions: explicitly outside the current model

- [ ] **Non-ideal phase thermodynamics.** Introduce activity/fugacity models
  with typed parameters, temperature/pressure validity ranges, and analytic or
  independently validated derivative contracts. Ideal-gas and pure-condensed
  activity semantics must remain selectable baselines, not implicit fallbacks.
- [ ] **Liquid/solid solution models.** Support real solution excess Gibbs
  models, composition-dependent chemical potentials, and phase-qualified
  parameter provenance. Reject unsupported solution models before a nonlinear
  backend starts; never silently treat them as ideal.
- [ ] **Phase coexistence and latent heat.** Define the physical contract for
  liquid-solid, liquid-vapor, and solid-solid coexistence, phase fractions,
  and latent contributions in `P,H` solves. Existing phase-control appearance
  evidence is not a general coexistence model.
- [ ] **Non-ideal/global TPD extension.** After P7 establishes the general TPD
  contract for the current ideal activity models, add model-specific
  non-convex minimization, deterministic multi-start/global policy, and
  independently validated chemical-potential derivatives. Do not weaken P7
  or reuse an ideal closed-form minimum for a non-ideal phase.
- [ ] **Additional constraints and state variables.** Design separate typed
  workflows for `P,V`, `U,V`, pressure sweeps, and reactive-flash problems.
  They must not be encoded as ad-hoc switches inside the fixed-pressure
  `P,T`/`P,H` request types.
- [ ] **Electrolyte/charged-species support.** Add electroneutrality,
  reference-state conventions, ionic-strength/activity models, and a clear
  contract for aqueous phases only when suitable local data and validation
  fixtures exist.

### F3. Quality-of-life and observability backlog

These items are safe to implement independently of new physics. They must
consume immutable solution/report snapshots and must not reopen `SubsData` or
change accepted numerical results.

- [x] **Human-readable report views.** `equilibrium_presentation` now projects
  an accepted immutable solution into deterministic phase, component/provenance,
  backend-attempt, and timing rows for GUI, CLI, or export consumers. Its
  compact ASCII renderer is a convenience view only; structured rows remain the
  canonical presentation contract and retain failure evidence.
- [x] **Range presentation and plotting series.** `TemperatureRangePresentationReport`
  now exposes stable ordered `P,T` component-mole and phase-total series,
  point-level continuation/timing/validation rows, and explicit phase-transition
  boundaries. `PhRangePresentationReport` now mirrors this with a correctly
  typed target-enthalpy axis plus solved-temperature/component/phase raw series
  and route/fallback/continuation evidence. Both range views also publish raw
  residual/balance/per-point-time columns and phase-local component mole
  fractions. `EquilibriumDisplayPolicy::visible_amount_columns()` returns
  display-column indices without rebuilding or dropping the canonical raw
  series, so plotting/UI code can apply trace filtering safely.
- [x] **Phase-aware postprocessing.** PCHIP remains an optional presentation
  layer, but typed `P,T` range resampling is now forbidden when accepted
  phase-control transitions occurred; callers keep raw points or explicitly
  split phase-stable ranges. Raw and resampled rows remain distinct, and the
  existing policy still makes the linear/log interpolation space explicit.
  Future range presentation must extend the same guard to failed gaps and any
  future layout-changing workflow.
- [x] **Run comparison report.** `EquilibriumComparisonReport` compares two
  accepted single-point solutions only after exact phase-qualified layout
  validation. It reports component/phase mole deltas, conditions, backend,
  balance/residual evidence, and lookup provenance changes. A future range
  comparison may aggregate these immutable point comparisons without altering
  the strict layout contract.
- [x] **Reproducibility capsule.** `equilibrium_reproducibility` now exports
  a JSON-compatible immutable snapshot of solved conditions, phase specs,
  selected per-component record identities/provenance, effective backend
  cascade and numerical options, acceptance evidence, optional candidate
  selection, catalog consistency, and a caller-supplied data-release label.
  The record fingerprint is deliberately an identity fingerprint and the
  catalog value is structural evidence, not a false claim of a payload-content
  hash; **F1** remains responsible for a versioned local-data manifest. This
  is provenance/report export, not mutable project import/export.
- [x] **Failure-focused diagnostics.** `backend_attempt_rows_from_error()`
  exposes the complete ordered trace retained by `AllBackendsFailed` and
  `CascadeAborted`, including termination, iterations, callback timing, and
  typed cause. Range-point/seed context remains owned by the range reports and
  can be paired with these rows without fabricating a successful solution.
- [x] **Display thresholds and units policy.** `equilibrium_display` now
  provides validated trace-row filtering, scientific/fixed/engineering-SI
  number styles, and explicit fraction-versus-percent rendering. It only
  borrows or projects immutable presentation rows: thresholding never alters
  solver inputs, conservation checks, or raw data retained for export.
- [x] **Opt-in structured phase-lifecycle diagnostics.**
  `equilibrium_diagnostics` now records bounded typed events for the canonical
  `PreparedPhaseControlRunner`: solve/outer-iteration boundaries, accepted
  fixed-set candidates, constrained TPD evidence, hysteresis holds, accepted
  activation/deactivation, boundary recovery, and final publication. The
  default is disabled; enabled reports are attached immutably to the accepted
  solution and an optional live sink can observe the same events. The separate
  `equilibrium_diagnostics_display` adapter resolves phase/component labels and
  can emit a human-readable tree through the existing `log` facade without
  putting formatting or logging side effects in the physical solver.
- [x] **Complete engine diagnostics coverage for rejected paths and long
  ranges.** Typed rejected-candidate, boundary-probe, rollback, outer-budget,
  cycle, and `P,H` `Auto` route-fallback events stream through the live sink.
  Bounded P,T and P,H ranges support endpoints, every-N, every-point, and
  `TransitionsOnly` policies. The last policy defers the sink until an
  accepted point actually changes the phase set, then retains/replays that
  point trace only. Accepted `P,H Auto` results retain an immutable
  `PhRouteDecision::AutoFallback` independently of scalar trials.
- [x] **Expose diagnostics in the equilibrium GUI.** The document now maps an
  explicit lifecycle detail level and range retention policy to the canonical
  engine options. Accepted point snapshots project the immutable typed report
  into a phase-qualified collapsible decision tree; P,H also exposes its
  route decisions independently from scalar trials. While a worker runs, the
  same opt-in engine sink streams bounded transient events through the
  ticket/fingerprint gate. GUI formatting remains a consumer of engine
  evidence, never a second source of phase-control decisions.

## Recommended implementation passes

1. [x] Characterize current fixtures and introduce typed problem/result/error
   boundaries without changing numerical behavior.
2. [x] Extract the pure log-moles formulation and common acceptance gate.
3. [x] Characterize the RST backend matrix, then complete attempt counters,
   termination mapping, and a measured production default.
4. [x] Isolate legacy solvers and keep them behind an explicit fallback
   policy with typed diagnostics.
5. [x] Build the independent K_eq validator and cross-validation fixtures.
6. [x] Bridge `ResolvedPhaseSystem`, migrate retained classical workflows, and
   remove duplicate orchestration from the production surface.
7. [ ] Build the typed GUI on the stable facade and reports.

## Historical definition of done

The original checklist below is retained for traceability and synchronized
with the current fixed-`P,T` production contract. Legacy orchestration removal
is intentionally deferred until compatibility consumers have migrated; the
handwritten LM/NR/TR implementations remain supported numerical fallbacks.

- [x] There is one canonical equilibrium problem and solution model.
- [x] Every accepted solution has passed backend-independent numerical and
  physical validation.
- [x] Fallback order, retry rules, budgets, and attempts are explicit and fully
  reported.
- [x] RustedSciThe provides the production nonlinear backends.
- [x] The K_eq path independently validates every applicable small-system
  fixture and clearly reports when it is not applicable.
- [x] Tests cover formulation, every backend, cascade behavior, physical
  invariants, cross-validation, transactionality, and offline integration.
- [x] Legacy mutable orchestration is absent from the production facade. Its
  deprecated compatibility entry points live under `ChemEquilibrium::legacy`,
  while legacy numerical fallback backends remain intentionally supported.

---

## SourceCraft Diagnostics

Результаты ревизии кода модуля `ChemEquilibrium`, проведённой 23.07.2026.

### Проверка диагностики (26.07.2026)

Каждый пункт ниже имеет итоговый бинарный статус: `✅ выполнен` либо
`❌ отклонён`. Формулировка «отклонён как локальное исправление» означает, что
сама архитектурная тема признана реальной, но предложенный SourceCraft патч
небезопасен или преждевременен; связанная миграция остаётся в основном плане.

- [x] **A.1 подтверждён и исправлен.** `R` теперь использует полное значение
  CODATA 2018 `8.314_462_618_153_24 J/(mol K)`; есть прямой regression-test.
- [x] **A.2 подтверждён и исправлен.** Автоматический RST путь требует либо
  полный `gibbs_sym` snapshot в solver order, либо успешный thermochemistry
  lookup для каждого requested substance. Частично заполненные кэши больше не
  считаются символическим контекстом.
- [x] **A.3 подтверждён и исправлен.** Мутирующий compatibility method переименован
  в crate-private `publish_reconstructed_moles`; чистое преобразование остаётся
  свободной функцией `compute_species_moles`.
- [x] **B.2 подтверждён и исправлен.** Неиспользуемый
  `VariableScalingContract` и тесты его изолированной арифметики удалены. Если
  variable scaling понадобится, его следует вводить только вместе с реальной
  передачей масштаба во все backend contracts.
- [x] **E.1/E.2 подтверждены и исправлены.** Численные допуски получили имена,
  а подробные дампы initial guess, inventory и stoichiometric matrix понижены
  с `info!` до `debug!`.
- [x] **A.4/C.1/C.2/C.3/C.4 отклонены как локальные SourceCraft-патчи.** Полная
  миграция mutable legacy facade к `EquilibriumProblem`/`PreparedEquilibriumProblem`
  и устранение ручной active-set reconstruction остаются отдельным этапом, а
  не безопасным локальным патчем.
- [x] **D.2 подтверждён и исправлен.** Добавлен serial-only
  `solve_for_T_range_with_phase_control`: он разделяет canonical обновление
  Gibbs/publish sweep с ordinary serial path, переносит accepted seed и phase
  mask между точками и покрыт offline NASA-gas regression test. Parallel
  phase-control sweep остаётся намеренно неподдержанным, поскольку ему нужна
  отдельная independent semantics.
- [x] **B.1, B.3, B.4, B.5, E.4-E.6 отклонены как дефекты.** Это осознанные
  контракты или policy choices: глобальный budget независим от per-attempt
  лимита, `Required` честно сообщает о неприменимости K_eq validation,
  explicit hysteresis имеет размерную семантику, а trace floors/default scaling
  нельзя менять без отдельных stability benchmarks.

### A. Потенциальные ошибки (bugs)

#### A.1 `R = 8.314` — неверное значение газовой постоянной

**Где:** [`equilibrium_log_moles.rs:100`](equilibrium_log_moles.rs:100)

**Проблема:** `pub const R: f64 = 8.314;` — это значение не является стандартной молярной газовой постоянной. Правильное значение: `8.314462618`. Ошибка в 0.00046 J/(mol·K) даёт систематическое смещение в логарифме константы равновесия: `ln(K) = -ΔG/(RT)`. При T=1000K, ΔG=-100 kJ/mol: `ln(K)_true = 100000/(8.31446*1000) = 12.027`, `ln(K)_current = 100000/(8.314*1000) = 12.028`. Ошибка мала (~0.01%), но накапливается в температурных сериях.

**Важность:** Средняя. Для инженерных расчётов ошибка незначительна, но для научных публикаций — недопустима.

**Статус (26.07.2026):** ✅ **Выполнен.** Константа заменена на полное
значение CODATA 2018; прямой regression-test защищает значение.

#### A.2 `has_rst_symbolic_context()` — некорректная эвристика

**Где:** [`equilibrium_log_moles.rs:490-496`](equilibrium_log_moles.rs:490)

**Проблема:** Метод проверяет наличие символьного контекста через `gibbs_sym.len() == substances.len() && !gibbs_sym.is_empty()`, но также возвращает `true`, если `search_results`, `search_states`, `therm_map_of_sym` или `therm_map_of_fun` непусты. Последние четыре условия не гарантируют, что символьные Gibbs-функции действительно построены для всех веществ. Это может привести к тому, что `SolverPolicy::rusted_scithe_default()` будет выбран, но RST-адаптер не сможет построить SymbolicNonlinearProblem.

**Важность:** Средняя. Может проявляться как трудноотлавливаемая ошибка "RST backend failed" при определённых последовательностях вызовов.

**Статус (26.07.2026):** ✅ **Выполнен.** Эвристика требует либо полного
`gibbs_sym` snapshot, либо успешного thermo lookup для каждого вещества;
частичные caches больше не включают RST автоматически.

#### A.3 `compute_species_moles` — публичная функция с побочным эффектом

**Где:** [`equilibrium_log_moles.rs:2173-2178`](equilibrium_log_moles.rs:2173)

**Проблема:** Публичный метод `compute_species_moles(&mut self, sol: Vec<f64>)` мутирует `self.moles` и `self.map_of_moles_for_each_substance`. Название предполагает чистую функцию (как у свободной функции `compute_species_moles(log_moles: &[f64])`), но на самом деле это метод с побочным эффектом. Это нарушает принцип наименьшего удивления.

**Важность:** Низкая. Косметическая проблема именования.

**Статус (26.07.2026):** ✅ **Выполнен.** Мутирующий метод стал
crate-private `publish_reconstructed_moles`; чистая свободная функция сохранила
имя `compute_species_moles`.

#### A.4 `EquilibriumLogMoles` — 22 публичных поля

**Где:** [`equilibrium_log_moles.rs:275-320`](equilibrium_log_moles.rs:275)

**Проблема:** Структура имеет 22 публичных поля, которые можно изменять извне в любом порядке. Нет гарантии, что после изменения одного поля (например, `T`) остальные поля (например, `gibbs`) остаются консистентными. Это прямой путь к багам "забыл обновить gibbs после смены T".

**Важность:** Высокая. Основной источник производственных ошибок.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.** Это реальный
долг legacy facade, но безопасный путь — дальнейшая миграция callers на
`EquilibriumProblem`/`PreparedEquilibriumProblem`, а не приватизация полей
одним разрушающим патчем.

### B. Overengineering

#### B.1 `SolverCascadeBudget` — избыточная сложность для текущих потребностей

**Где:** [`equilibrium_solver_policy.rs:38-46`](equilibrium_solver_policy.rs:38)

**Проблема:** Бюджет каскада содержит три независимых лимита (`max_attempts`, `max_iterations_per_attempt`, `max_total_iterations`), которые проверяются в [`solve_backend_cascade`](equilibrium_log_moles.rs:2064-2073). При этом `max_total_iterations` вычисляется как `max_iter * backends.len()` — то есть всегда пропорционален двум другим. Три лимита вместо одного — это overengineering. Достаточно двух: `max_attempts` и `max_iterations_per_attempt`.

**Важность:** Средняя. Усложняет понимание без реальной выгоды.

**Статус (26.07.2026):** ❌ **Отклонён.** `max_total_iterations` не обязан
равняться произведению двух остальных лимитов: он задаёт независимый глобальный
resource cap и уже проверяется regression-тестами.

#### B.2 `VariableScalingContract` — объявлен, но нигде не используется

**Где:** [`equilibrium_problem.rs:313-384`](equilibrium_problem.rs:313)

**Проблема:** `VariableScalingContract` полностью реализован (с `apply_iterate`, `unscale_iterate`, валидацией), но не используется ни в одном бэкенде. Row scaling (`ResidualScalingContract`) используется, а variable scaling — нет. Это мёртвый код.

**Важность:** Средняя. Увеличивает когнитивную нагрузку при чтении.

**Статус (26.07.2026):** ✅ **Выполнен.** Неиспользуемый
`VariableScalingContract` удалён вместе с изолированными тестами; при реальной
потребности его следует вводить сразу через backend contracts.

#### B.3 `EquilibriumConstantSolverMode::Required` — нереализуемый контракт

**Где:** [`equilibrium_constant_solver.rs:74-75`](equilibrium_constant_solver.rs:74)

**Проблема:** Режим `Required` требует, чтобы K_eq валидация была выполнена для каждого решения. Но K_eq решатель работает только для однореакционных систем. Для многореакционных систем `Required` гарантированно упадёт с ошибкой. Этот режим нельзя использовать в production, только в тестах.

**Важность:** Низкая. Но вводит в заблуждение.

**Статус (26.07.2026):** ❌ **Отклонён.** `Required` — намеренный строгий
контракт для workflows, где независимая validation обязательна. Неприменимость
к многореакционной задаче возвращается как typed error, а не маскируется.

#### B.4 `TemperaturePostprocessingPolicy` — избыточная гибкость

**Где:** [`equilibrium_temperature_postprocessing.rs`](equilibrium_temperature_postprocessing.rs)

**Проблема:** Политика постобработки поддерживает три типа сеток (RawOnly, Uniform, Explicit) и два типа интерполяции (Linear, Log). При этом в коде нет ни одного вызова с Explicit или Log. Вся гибкость существует только для тестов.

**Важность:** Низкая. YAGNI-нарушение.

**Статус (26.07.2026):** ❌ **Отклонён.** Это публичная policy граница для
postprocessing и GUI/CLI consumers; отсутствие текущих production callers не
делает корректные modes мёртвым кодом.

#### B.5 `PhaseHysteresisPolicy::Explicit` — дублирование TemperatureScaled

**Где:** [`equilibrium_workflows.rs:1017-1025`](equilibrium_workflows.rs:1017)

**Проблема:** `Explicit { dg_create, dg_keep }` — это фактически `TemperatureScaled` с замороженной температурой. При T=const они эквивалентны. Можно было бы обойтись одним вариантом.

**Важность:** Низкая.

**Статус (26.07.2026):** ❌ **Отклонён.** `TemperatureScaled` хранит
безразмерные множители `R*T`, а `Explicit` — физические пороги `ΔG` в J/mol.
Они намеренно различаются при температурных сериях.

### C. Недостатки дизайна

#### C.1 `EquilibriumLogMoles` — God Object

**Где:** [`equilibrium_log_moles.rs:275-320`](equilibrium_log_moles.rs:275)

**Проблема:** Структура содержит 22 поля, объединяющие термохимические данные, настройки решателя, состояние температурной серии, опубликованные решения, отчёты валидации и фазы. Это классический God Object. Новый API (`EquilibriumProblem` + `EquilibriumSolution`) уже решает эту проблему, но старый фасад остаётся в production.

**Важность:** Высокая. Затрудняет тестирование и поддержку.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.** Новый typed
pipeline уже служит source of truth для новых bridge workflows; legacy facade
будет сужаться по мере миграции, а не заменяться внезапно.

#### C.2 Две параллельные иерархии ошибок

**Где:** [`equilibrium_nonlinear.rs:55-72`](equilibrium_nonlinear.rs:55) и [`equilibrium_nonlinear.rs:89-188`](equilibrium_nonlinear.rs:89)

**Проблема:** `SolveError` (для численных решателей) и `ReactionExtentError` (для всего остального) — две пересекающиеся иерархии. `SolveError` мог бы быть вариантом `ReactionExtentError`, но они разделены. Это приводит к тому, что в некоторых местах ошибка оборачивается (`ReactionExtentError::SolveError`), а в некоторых — нет.

**Важность:** Средняя.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.**
`SolveError` остаётся внутренним численным уровнем, а `ReactionExtentError`
сохраняет доменный контекст и cascade trace. Слияние требует отдельного API
решения, не механического переименования.

#### C.3 `from_problem` — деструктуризация типобезопасности

**Где:** [`equilibrium_log_moles.rs:615-644`](equilibrium_log_moles.rs:615)

**Проблема:** Метод `from_problem` принимает типобезопасный `EquilibriumProblem`, но раскладывает его в 22 публичных поля `EquilibriumLogMoles`. Вся типобезопасность теряется. Это мост для миграции, но он должен быть временным.

**Важность:** Средняя.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.**
`from_problem` —
контролируемый migration bridge; его устранение возможно только после сужения
legacy facade и переноса оставшихся callers.

#### C.4 `solve_fixed_active_set_candidate` — клонирование всего solver'а

**Где:** [`equilibrium_workflows.rs:1439-1476`](equilibrium_workflows.rs:1439)

**Проблема:** Для решения с проекцией активного набора создаётся полная копия `EquilibriumLogMoles` (`let mut local = EquilibriumLogMoles::empty()`), в которую копируются поля по одному. Это 30+ строк ручного копирования. Любое новое поле в `EquilibriumLogMoles` нужно не забыть добавить и сюда.

**Важность:** Средняя. Хрупкий код.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.** Ручная
реконструкция active-set solver изолирована в одном private workflow и покрыта
projection regressions, но должна исчезнуть при окончательном переносе этого
пути на immutable prepared-problem snapshots.

### D. Пробелы в тестовом покрытии — статус

> **Статус:** 70 из 70 пунктов закрыты. D.2 закрыт отдельным serial-only
> `solve_for_T_range_with_phase_control` regression path.

#### D.1 ✅ `EquilibriumLogMoles::solve()` с Legacy бэкендами

**Статус:** 5 тестов в `equilibrium_log_moles_tests.rs` (`legacy_lm_solve_publishes_accepted_solution`, `legacy_nr_solve_publishes_accepted_solution`, `legacy_tr_solve_publishes_accepted_solution`, `legacy_solve_fails_gracefully_without_stoich_matrix`, `legacy_solve_publishes_solve_report_with_attempts`).

#### D.2 ✅ `solve_for_T_range` с phase control

**Статус:** `solve_for_T_range_with_phase_control()` выполняет bounded
`solve_with_phase_control()` для каждой точки последовательной температурной
серии, сохраняет accepted continuation seed/phase mask и публикует только
полностью согласованную sweep table. Offline NASA-gas regression:
`phase_controlled_temperature_sweep_publishes_each_accepted_point`.

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

**Статус (26.07.2026):** ✅ **Выполнен.** Все три значения получили
семантические имена рядом с numerical acceptance boundary.

#### E.2 Избыточное логирование

В [`solve_candidate_from_seed`](equilibrium_log_moles.rs:1908-1910) три `info!` подряд с полным дампом initial_guess, n0 и stoich_matrix. Для production это слишком подробно. Должно быть `debug!`.

**Статус (26.07.2026):** ✅ **Выполнен.** Полные дампы переведены на
`debug!`; штатный `info!` больше не засоряется внутренним состоянием solver.

#### E.3 Дублирование кода в `solve_for_T_range`, `solve_for_T_range_par`, `solve_for_T_range_par2`

Три реализации температурной серии содержат значительное дублирование логики обновления Gibbs, валидации и сохранения результатов. Можно было бы выделить общий шаг итерации.

**Статус (26.07.2026):** ❌ **Отклонён в исходной формулировке.** Обычный
serial sweep и serial phase-control sweep используют общие
`refresh_temperature_gibbs` и `publish_temperature_sweep`. Два parallel path
намеренно остаются отдельными: у них другой ownership/continuation contract,
который нельзя скрыть одной общей итерацией.

#### E.4 `PHASE_CONTROL_TRACE_MOLE_FLOOR = 1e-300`

**Где:** [`equilibrium_workflows.rs:31`](equilibrium_workflows.rs:31)

Значение `1e-300` находится на грани представимости в f64 (min positive normal ≈ 2.2e-308). `ln(1e-300) ≈ -690.8`, что далеко от -Inf, но при дальнейшем делении на phase totals может привести к underflow. Рекомендуется `1e-100` или документировать риск.

**Статус (26.07.2026):** ❌ **Отклонён.** `ln(1e-300)` конечен, а floor нужен
именно для inactive phase coordinates. Менять физически чувствительный порог
без stability benchmark нельзя; его назначение задокументировано и покрыто
phase-control regressions.

#### E.5 `DEFAULT_TRACE_MOLE_FLOOR = 1e-30`

**Где:** [`equilibrium_problem.rs:28`](equilibrium_problem.rs:28)

`ln(1e-30) ≈ -69.1`. Это безопасно для f64. Замечаний нет.

**Статус (26.07.2026):** ❌ **Отклонён как проблема.** Диагностика сама
подтверждает безопасность значения; это именованный, валидируемый seed policy.

#### E.6 Неконсистентность: `scaling_flag` в `EquilibriumSolverSettings` по умолчанию `false`

**Где:** [`equilibrium_log_moles.rs:153`](equilibrium_log_moles.rs:153)

Масштабирование выключено по умолчанию, хотя в архитектурном обзоре оно описано как важная фича. Рекомендуется включить по умолчанию.

**Статус (26.07.2026):** ❌ **Отклонён.** Default должен сохранять
историческую численную семантику. Масштабирование включается явным policy
параметром; изменение default требует comparative stability benchmark.

---

### Приоритеты для исправления

| Приоритет | Категория | Пункт | Итог |
|-----------|-----------|-------|------|
| P0 | Bug | A.4 | ❌ Отклонён как локальный патч; нужен отдельный typed migration этап |
| P0 | Design | C.1 | ❌ Отклонён как локальный патч; `EquilibriumProblem` остаётся source of truth |
| P1 | Bug | A.2 | ✅ Выполнен: RST readiness теперь требует полного контекста |
| P1 | Bug | A.1 | ✅ Выполнен: CODATA 2018 `R` |
| P1 | Overengineering | B.1 | ❌ Отклонён: global resource cap независим от per-attempt лимита |
| P1 | Overengineering | B.2 | ✅ Выполнен: мёртвый `VariableScalingContract` удалён |
| P1 | Test | D.1-D.70 | ✅ Выполнено: 70 из 70, включая serial phase-control sweep |
| P2 | Design | C.2 | ❌ Отклонён как локальный патч; требует отдельной error-model migration |
| P2 | Design | C.4 | ❌ Отклонён как локальный патч; часть prepared-problem migration |
| P2 | Code | E.1 | ✅ Выполнен: numerical constants именованы |
| P2 | Code | E.4 | ❌ Отклонён: trace floor нельзя менять без benchmark |
| P3 | Code | E.2 | ✅ Выполнен: подробный logging понижен до `debug!` |
| P3 | Code | E.3 | ❌ Отклонён в исходной форме; serial helpers выделены, parallel semantics отдельны |

Итог: SourceCraft Diagnostics зафиксирован как журнал принятых исправлений и
отклонённых локальных предложений; открытые архитектурные миграции ведутся
отдельными пунктами основного плана.
