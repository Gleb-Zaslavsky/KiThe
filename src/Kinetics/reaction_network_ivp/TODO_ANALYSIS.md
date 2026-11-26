# Reaction-Network IVP Implementation Plan

Status: public numerical foundation, broad correctness matrix, and release
performance baseline implemented. SolidState v1 is integrated into the
prepared vector IVP. All KiThe-owned checklist items below are closed;
upstream-dependent callback context and failure provenance remain separately
identified. The physical contract is in
`MODEL.md`. Follow
`KiThe_ARCHITECTURAL_GUIDELINE.md` throughout.

## Purpose And Scope

Build a public direct-kinetics calculator for user-defined reaction networks.
The user supplies species, reactions, rate laws, initial concentrations, an
integration interval, and a prescribed temperature program T(t). KiThe builds
the symbolic IVP and delegates numerical integration to RustedSciThe (RST).

This is a third domain workflow, not another numerical solver implementation:

- `SimpleReactorsIVP` owns the coupled spatial combustion model.
- `solid_state_kinetics_IVP` owns the existing scalar conversion workflow and
  the typed solid-state kinetic model catalog.
- `reaction_network_ivp` owns arbitrary reaction networks in physical time.

Initial physical scope: closed, homogeneous, constant-volume kinetics. No
energy balance, thermochemistry lookup, equation of state, diffusion, spatial
transport, or automatic reverse rates inferred from equilibrium constants.
Concentrations need not sum to one and are not mass or mole fractions.

Current foundation: validated public problem/rate/temperature types, symbolic
RHS assembly, BDF/Radau/Backward Euler/LSODE2 StiffAPI integration, prepared
initial-state, constant/linear-temperature, and complete reaction-coefficient
rebinding, cooperative cancellation, transactional publication, output
validation, and independently switchable native telemetry/logging. Stories
cover analytical linear/nonlinear networks, stiff nonlinear scale behavior,
prescribed heating, backend/frontend parity, rebinding, invariants, and
cancellation. KiThe-owned stage diagnostics report equation assembly, RST
preparation, integration, and combined output postprocessing. RN-09 has a
release-measured baseline and the supported matrix-layout/callback-policy
axes; no optimization claim is made from one machine's measurements.

### SolidState v1 Decisions

- **Reaction capacity and shared-resource policy:** v1 uses explicit per-stage
  extent capacity `q_j` and initial `alpha_j`; neither is inferred from
  concentrations. Reject a network with two stages consuming the same species
  whenever at least one is SolidState, regardless of whether the other stage
  is SolidState, MassAction, or CustomExpr. Identify both reaction IDs and the
  shared species in a typed construction error. A resource-depletion law is a
  separate future model; do not add an implicit depletion multiplier.
- **Initial resource feasibility:** for each species with net consumption
  `nu[j,i] < 0` by a SolidState stage, require
  `C_i(0) >= -nu[j,i] * q_j * (1-alpha_j0)` at construction. `C_i(0)` is
  the concentration at the start of this solve, and the right-hand side is
  the stage's remaining consumption capacity. This intentionally does not
  count production by other stages or future feeds in v1; reject a shortfall
  with a typed error containing the stage, species, required and available
  concentrations. This check does not replace runtime nonnegative-state
  validation.
- **Catalog endpoints:** SolidState conversion has the physical range
  `0 <= alpha <= 1`. Validate the solver-admissible domain of each
  `KineticModel` (including its parameters). An endpoint is supported only
  when the RHS and state Jacobian are finite and the stage does not drive
  conversion past one; otherwise an explicit analytic endpoint treatment
  would be needed. Reject an invalid initial value or out-of-domain solver state
  with a typed error identifying the stage and model; never epsilon-shift,
  clamp, or silently normalize alpha.
- **Conversion direction:** SolidState `alpha_j` is nondecreasing in v1:
  `alpha_j' = k_j(T)*f_j(alpha_j) >= 0`. Reject parameter ranges that make
  the catalog rate negative on its admitted domain; check evaluated rates
  too, with a typed stage/model error if a negative or nonfinite value occurs.
  Do not clamp a negative rate to zero. Model reverse conversion as a separate
  reaction, not by giving a SolidState stage a negative rate.
- **A1 naming:** keep A1 absent; do not alias it to F1. Add it only with an
  explicit formula and cited model convention. Existing F1 remains unchanged.
- **Solver pool:** intentionally limited to the four current `StiffAPI`
  families: LSODE2, standalone BDF, Radau, and Backward Euler. Do not add
  non-stiff methods unless a concrete use case justifies the additional
  settings/capability/telemetry matrix.

## RN-01: Physical And State Contracts

- [x] Document time in seconds, temperature in kelvin, concentrations in
  mol/m^3, and reaction rates in mol reaction/(m^3*s). Document activation
  energy in J/mol and rate-order-dependent units of pre-exponential factors.
- [x] Use the existing reaction-by-species stoichiometric convention:
  nu[j,i] = products[j,i] - reactants[j,i], dC[i]/dt = sum_j nu[j,i]*r[j].
- [x] Keep reaction order distinct from stoichiometry. Elementary MassAction
  uses reactant stoichiometry; any supported empirical orders must be explicit.
- [x] Pin the example A -> 2B; B + C -> D to A'=-r1, B'=2*r1-r2,
  C'=-r2, D'=r2. Never multiply a rate by product stoichiometry twice.
- [x] Define deterministic species/reaction ordering, stable IDs, and safe
  internal symbolic names distinct from user labels. Reject duplicate IDs and
  ambiguous references; permit distinct reaction IDs with the same equation.
- [x] Allow abstract labels such as A and B without requiring molecular
  formulas, molar masses, or thermodynamic records. Elemental conservation
  checks require separately supplied valid composition data.
- [x] Treat explicit forward/reverse reactions as independent stages unless
  a later typed reversible-rate contract is added. Do not infer kinetics from
  a reversible arrow without a defined reverse law.
- [x] Decide and document empty-network and inert-species behavior; reject
  empty species sets and unmatched initial-condition entries explicitly.
  Empty reaction sets and inert species are now covered by a BDF solve story.

- [x] Reuse `solid_state_kinetics_IVP::KineticModel` and `get_expr()`;
  do not duplicate formulas or launch a nested scalar IVP per reaction.
- [x] Use the approved stage-specific representation:
  alpha_j' = k_j(T)*f_j(alpha_j), r_j = q_j*alpha_j', where q_j is the stage
  capacity in mol reaction/m^3; publish alpha_j separately from concentrations.
- [x] Require explicit q_j and alpha_j0; never substitute
  alpha_j = 1-C_A/C_A0 or infer capacity from initial concentration.
- [x] Reject any pair of stages that consume the same species when at least
  one stage is SolidState, including a MassAction or CustomExpr counterpart;
  report both reaction IDs and the species in a typed construction error.
  A later coupled-resource model requires a separately agreed depletion law;
  do not clamp concentrations or invent one.
- [x] Require initial concentration to cover each SolidState stage's full
  remaining net consumption: `C_i(0) >= -nu[j,i]*q_j*(1-alpha_j0)` where
  `nu[j,i] < 0`. Do not credit later production or feed in v1; report
  available/required amounts in a typed construction error.
- [x] Choose the physical alpha/domain policy: `0 <= alpha <= 1`, with
  model- and parameter-specific solver-admissible domains and explicit
  endpoint limits only where implemented; no epsilon shift or clamp.
- [x] Make SolidState conversion nondecreasing in v1; reject model parameters
  admitting a negative rate, and reject negative/nonfinite evaluated rates
  with typed stage/model context. No rate clipping or implicit reverse stage.
- [x] Classify finite-real alpha domains and nonnegative-rate parameter
  restrictions for the catalog; validate initial/output states and the
  explicit KiThe RHS probe, with interior and endpoint regression tests.
  Native intermediate callbacks still need typed stage/model context (RST
  dependent, below); multi-point parameter sweeps remain a follow-up.
- [x] Keep A1 absent and do not alias it to F1 without an explicit formula
  and model convention.

## RN-02: Typed Inputs And Public Facade

- [x] Design validated problem/configuration types with one authoritative
  representation for species, reactions, initial state, and temperature.
- [x] Provide typed rate-law variants for MassAction, explicit empirical-order
  MassAction, and complete CustomExpr. `CustomExpr` is not supplemented by
  hidden concentrations or Arrhenius terms. SolidState uses a separate
  conversion coordinate and explicit stage capacity under RN-01.
- [x] Support constant/Arrhenius rate coefficients without making thermodata
  lookup a prerequisite. Validate numeric parameters according to their law.
- [x] Expose public prepare/solve/result access without exporting mutable
  internals that bypass validation or leave cached equations stale.
- [x] Store explicit reaction-by-species stoichiometry and keep empirical
  orders separate from it. `KinData` is intentionally not required for the
  direct user-defined network API; any future import adapter needs its own audit.
- [x] Keep direct user-defined mechanisms independent of database I/O; a
  later library-import adapter should preserve reaction/library provenance.

## RN-03: Prescribed Temperature And Symbolic Assembly

- [x] Accept a symbolic temperature program T(t), with constant and linear
  convenience constructors. T is prescribed, not an integrated thermal state.
- [x] Substitute T(t) structurally into rate expressions with Expr APIs; do
  not manipulate serialized expression strings to rename or bind symbols.
- [x] Validate `CustomExpr` symbols against `t`, `T`, and in-range canonical
  state names `c0..cN`; report reaction ID on invalid rate symbols. Arbitrary
  declared custom parameters and conversion states are not currently part of
  this API. Temperature expressions accept only `t`.
- [x] Assemble one symbolic RHS vector for the full state. RST prepares the
  state Jacobian at fixed independent time; KiThe does not duplicate it.
- [x] Validate finite positive temperature at construction endpoints and at
  every published sample. For symbolic T(t), endpoint checks do not prove
  positivity at internal solver trial states; RST callback failures remain
  native until it exposes richer evaluation context.
- [x] Do not add implicit clipping or domain-changing rewrites for fractional
  powers/logarithms. Invalid runtime evaluations are left visible to the
  numeric callback/backend boundary; richer expression-role context is open.
- [x] The public input is an already-built RST `Expr`, not a string parser, so
  this layer has no parse/panic boundary. Revisit if a text-expression API is
  introduced.
- [x] Defer piecewise/tabulated temperature programs: current `Constant`,
  `Linear`, and symbolic `T(t)` inputs are sufficient; no concrete use case
  currently requires a table. Revisit only with a real profile requirement.

## RN-04: RST Backend Integration And Prepared Lifecycle

- [x] Use the vector StiffAPI route, not the scalar `solve_1d()` preset.
- [x] Cover all current StiffAPI families: LSODE2, standalone BDF, Radau, and
  Backward Euler. Distinguish LSODE2's internal BDF/Adams controller from the
  standalone BDF solver.
- [x] Keep public solver support limited to the four StiffAPI families for
  now; non-stiff RST methods remain out of scope unless a concrete need arises.
- [x] Reuse the reactor adapter's StiffAPI preparation pattern, not its
  thermal-state offset,
  spatial coordinate conventions, or reactor-only stopping conditions.
- [x] Preserve supported shared controls and reject unsupported options:
  integration tolerances (where supported), time interval, output policy,
  frontend, matrix layout, callback policy, telemetry, and logging.
- [x] Cover ExprLegacy/AtomView and the native callback route currently exposed
  by StiffAPI; sequential is the conservative default. AOT is not represented
  by this API and should not be claimed as covered.
- [x] Reject unsupported combinations before integration, without silent
  substitution of solver family, matrix backend, or execution policy.
- [x] Separate preparation from integration. Reuse prepared residual/Jacobian
  machinery for constant/linear temperature and initial-state changes using
  RST `try_solve_with_parameter_values`.
- [x] Rebind numeric MassAction/EmpiricalMassAction rate coefficients without
  rebuilding symbolic callbacks. `RateCoefficientBinding` supplies a complete
  set keyed by reaction ID; ordering is irrelevant, while omissions, duplicate
  IDs, unknown IDs, kind mismatches, and invalid numeric values are rejected
  before RST is touched. Constant binds one value; Arrhenius binds `(A,Ea)`.
  Combined temperature+rate binding is available. Only accepted solves replace
  the accepted binding; rebind validation or integration failure preserves the
  last published result. Arbitrary symbolic T(t), custom-expression parameters,
  and rate-law shape changes still require rebuilding.
- [x] Define invalidation boundaries: schema/topology/rate-law changes rebuild;
  supported numeric parameter rebinding must not trigger symbolic rebuilding.
- [x] Forward cooperative cancellation through supported RST APIs; preserve
  native diagnostics when exposed. Cancellation is not
  convergence, acceptance, or permission to publish a partial final result.
- [x] Keep cancellation transactional: any native partial state is diagnostic
  only and must never replace `last_solution`; cancellation remains a typed
  failure. The last accepted solution is retained.
- [x] Characterize native mid-run cancellation/partial telemetry availability
  with an ignored wall-clock story. On the current RST 0.5.2 fixture BDF and
  Backward Euler expose telemetry and 32 partial samples, Radau exposes
  telemetry and its one-point failure snapshot, while LSODE2 finishes before
  the signal and reports telemetry for the completed attempt. This is measured
  behavior, not a deterministic cancellation-latency guarantee for LSODE2.

## RN-05: Rich Typed Errors From The Start

- [x] Design typed categories for input validation, species/reaction identity,
  stoichiometry/rate configuration, temperature, unsupported capability,
  backend preparation/integration, cancellation, and output validation.
  There is no KiThe expression-string parser or export layer in this API.
- [x] Carry structured context at KiThe-owned validation/publication
  boundaries: field path, reaction/species ID, output species/time/value, or
  sample index as applicable. Do not parse English error messages to classify
  failures.
- [x] Localize the native callback-context limit and specify the upstream
  request below. KiThe cannot truthfully add reaction/expression role or
  trial-state values to an opaque RST callback error.
- [x] Preserve native RST errors through `std::error::Error::source()` and
  typed variants. Avoid public `Result<_, String>` and lossy blanket conversion.
- [x] Preserve attempt-scoped telemetry and partial trajectory from RST reports
  in distinct diagnostics, never as an accepted solution. Reports and partial
  samples carry a `solve_id`; pre-start cancellation after a prior success is
  covered across all four families. Immediate callback failures can still
  return only the ID and typed error for BDF/LSODE2; this known upstream gap is
  recorded under Depends On RST Capability/API.
- [x] Separate malformed input, invalid runtime model domain, numerical
  non-convergence, and physically rejected candidate results.
- [x] Test local error categories/context and preserve native error sources
  across preparation and cancellation boundaries. Coverage includes invalid
  inputs/symbols/temperature, unsupported capabilities, cancellation, output
  values/times, and non-finite diagnostic RHS.

## RN-06: Opt-In Logging And Telemetry

- [x] Logging and telemetry are independent settings, both Off by default.
  Enabling either must not implicitly enable the other.
- [x] Use repository modes exactly: Off, Counters, Timings, Detailed. Off must
  avoid telemetry clock reads/counter updates; timings-only must not fabricate
  counters. Reject unsupported mode/backend combinations rather than aliasing.
  The facade forwarding matrix is implemented; native mode semantics remain
  owned by and must be verified against RST.
- [x] Forward the selected telemetry mode into the actual RST IVP solver.
  Symbolic preparation/evaluator detail remains part of the native family
  snapshot and is not flattened into fabricated KiThe counters.
- [x] Keep native family-tagged RST snapshots separate from KiThe stage data.
  KiThe currently reports equation assembly, native preparation, integration,
  and combined output postprocessing/validation.
- [x] Keep validated `NetworkProblem` construction outside prepared-session
  timings: it precedes the prepared lifecycle and no evidence shows it to be a
  material cost. Current totals distinguish equation assembly, backend
  preparation, integration, and output postprocessing. Split output extraction,
  materialization, validation, or publication only if measurement justifies it;
  native nested timings are not disjoint budget slices.
- [x] Distinguish not requested from measured zero; mode-specific absent
  counters/timings remain `None`. Native solver and symbolic evaluator
  counters describe different levels and need not be equal.
- [x] Preserve KiThe stage totals across integration and postprocessing errors
  on `PreparedNetworkIvp`; native partial snapshots are retained where RST
  exposes them. Constructor failure cannot return a prepared-session snapshot.
- [x] Provide opt-in native `SolverLogLevel` through the common StiffAPI
  preset, defaulting to `Off`; native capability/configuration remains owned by
  RST. KiThe does not initialize a logger, select a path, or print in solve.
  Effective level forwarding is asserted for all four native families.
- [x] Do not add a second KiThe iteration logger: native RST logs and the
  separate KiThe stage snapshots cover the current diagnostics contract.
  Revisit only if callers identify an event/context not exposed by either.
- [x] Test default-Off behavior, independent settings, all four modes, and
  forwarding to each StiffAPI family. Native mode semantics remain an RST
  contract. RST reports may omit telemetry when a callback fails before the
  native family records a failed status; KiThe preserves that absence.

## RN-07: Acceptance, Results, And Transactionality

- [x] Validate output shape, finite values, ordered time samples, terminal
  status/time coverage, and model-specific admissibility independently of
  solver success. Distinguish localized events from output truncation.
- [x] Keep the current strict concentration policy: reject every published
  value below zero; do not introduce a roundoff-negative corridor, clipping,
  or renormalization. Revisit only if a required case demonstrates a concrete
  roundoff-only rejection and a conservation-preserving correction is defined.
- [x] Do not make automatically discovered stoichiometric nullspace relations
  a generic acceptance gate. Keep known invariants as explicit story
  assertions. A future public invariant check must be opt-in with caller-typed
  weights and scale-aware tolerances; element/mass claims require composition
  data.
- [x] Publish named concentrations, time, prescribed temperature, and optional
  native telemetry. Conversion histories are blocked on conversion semantics;
  acceptance evidence and output interpolation metadata remain open.
  Separate native samples from requested/interpolated output.
- [x] Build and validate a candidate before atomically replacing the last
  accepted result. Failure/cancellation must retain the previous valid result.
- [x] Keep plotting, table formatting, and file export outside this numerical
  module; no such I/O occurs in integration or output materialization.

## RN-08: Correctness And Story Evidence

- [x] Test stoichiometry, ordering, rate semantics, temperature substitution,
  symbolic state derivatives against independent finite differences, and
  typed failure context for covered boundaries.
- [x] A -> 2B at constant T: analytical concentration curves and 2*A+B invariant.
- [x] A -> 2B; B+C -> D: componentwise RHS, integrated invariant C+D,
  non-negative output, and same-input StiffAPI family parity are covered.
- [x] Empty reaction set preserves inert species exactly on the default BDF
  route.
- [x] Explicit reversible first-order pair with analytical linear-system
  reference, a stiff consecutive first-order network, and a separate stiff
  nonlinear logistic reference across all four StiffAPI families.
- [x] Nonlinear prescribed-temperature law with a closed-form integrated
  trajectory, checked at all published samples for both frontends and all
  four solver families.
- [x] Add an analytical nonlinear autocatalytic network with ExprLegacy and
  AtomViewNative across all four families; validate each published sample
  against its logistic closed form and the conserved total.
- [x] Scalar F1: agreement with analytical alpha(t) and the existing scalar
  KineticModelIVP under identical model/temperature conventions.
- [x] Verify a mixed MassAction/SolidState/CustomExpr network and typed
  rejection of SolidState/CustomExpr, SolidState/MassAction, and two
  SolidState stages sharing a net-consumed species; cover catalog interior
  domains and representative singular endpoints.
- [x] Independent state-Jacobian finite-difference parity and a
  first-order-system scale characterization at `1e-3..1e3 mol/m^3` across
  all current frontend/family combinations.
- [x] Add a genuinely stiff nonlinear analytic network and concentration-scale
  coverage beyond linear first-order kinetics. Two independent nonlinear
  logistic subsystems have a 1000:1 rate-scale ratio; both are checked against
  their closed forms at `1e-3`, `1`, and `1e3` concentration scales. This is a
  reference fixture, not a universal nonlinear-stiffness guarantee.
- [x] Early cancellation and retention, unsupported capability, telemetry
  forwarding, malformed symbols, and default logging-disabled lifecycle are
  covered. A malformed middle reaction rejects whole problem construction.
  Attempt identity, pre-start cancellation, and native failure retention are
  covered; the availability matrix for immediate callback failures is pinned
  in Story 26, including the remaining BDF/LSODE2 RST limitation.
- [x] Add curated STORY_TESTS.md: name, question/hypothesis, exact commands,
  fixtures/options, numerical criteria/results, conclusions, and limitations.

## RN-09: Performance And Campaigns After Correctness

- [x] Add a Criterion baseline for 3/15/40-state linear chains and all current
  StiffAPI families, separating cold prepare+solve from repeated prepared
  solves. Detailed diagnostics are a separate telemetry-on pass.
- [x] Extend the baseline with ExprLegacy/AtomViewNative comparisons across
  sizes and all four families; matching analytical scale parity is asserted.
- [x] Extend across matrix layout and callback policy where each family
  supports the combination; add the parallel crossover only as an
  apple-to-apple axis. Unsupported combinations are reported explicitly.
- [x] Run the release campaign and record actual measurements before making
  performance claims. Results are a baseline, not a performance gate.
- [x] Emit a Tabled stage-attribution report for assembly, backend preparation,
  integration, output postprocessing, and unattributed wall time. Stage times
  are concrete KiThe phases; native RST detail remains a separate snapshot.
- [x] Keep complementary evidence: executable correctness assertions,
  Criterion benchmarks, machine-readable campaign records, and story findings.
- [x] Add a non-fail-fast PowerShell matrix with bounded build jobs, plan-only,
  suite selection, per-step exit codes, release profile, and aggregate status.
- [x] Store timestamped release output in `test_reports/<run>/reports` and raw
  Cargo/solver output in `technical`. No file-I/O performance gates.
- [x] Add an example through the public facade. GUI/task-document integration
  is a later explicit phase, not a prerequisite for the numerical API.

## Remaining Work And Ownership

The following are intentionally not implemented by inference. SolidState v1
choices and its basic numerical route are implemented; verify broader model
and backend coverage before expanding the model.

### Resolved Domain/API Contract

- SolidState v1 now has an explicit capacity contract: each stage supplies
  `q_j` and `alpha_j0`; there is no concentration-derived capacity or
  `alpha=1-C/C0` shortcut. Shared consumption by a SolidState stage and any
  other stage is rejected at construction, with both reaction IDs and the
  species in the error. A coupled resource-limited extension remains a later
  model decision requiring an explicit depletion closure. The physical
  conversion interval is `[0, 1]`, further restricted by each model's finite
  real evaluation domain; unsupported endpoints and out-of-domain trial
  states fail with typed stage/model context, never implicit clipping. Each
  net-consumed species must have enough initial stock for the SolidState
  stage's full remaining capacity, regardless of future production. SolidState
  conversion is nondecreasing; invalid negative or nonfinite stage rates are
  errors, not reverse conversion or zero-clipped rates.
- A1/F1 naming is resolved for now: no A1 alias; add an A1 model only after
  its exact expression and model convention are supplied. Scalar F1 parity
  and mixed SolidState-network stories cover the capacity/domain contract
  above (RN-08).
- Strict rejection of every negative concentration is the chosen publication
  policy; no clipping or renormalization. Automatically derived nullspace
  relations will not reject generic results. A future public invariant check
  must be explicitly requested with caller-declared weights/tolerances; no
  elemental/mass invariant is assumed without composition data.
- Public solver scope is resolved for now: only the four `StiffAPI` families;
  non-stiff RST methods may be reconsidered for a concrete use case.
- Piecewise/tabulated temperature is intentionally deferred because it is not
  currently needed. Existing `Constant`, `Linear`, and symbolic `T(t)` remain
  the supported set; if a use case appears, define knot/interpolation and
  discontinuity/restart semantics before implementation.

### Depends On RST Capability/API

- [ ] **Deferred: cancellable solve with rebound parameters.** `StiffAPI`
  currently exposes parameter rebinding and `CancellationToken` through
  separate solve paths. Do not emulate cancellation in KiThe with a token
  pre-check: that cannot interrupt an active native solve. The planned
  experimental-kinetics batch workflow may run its temperature conditions
  sequentially without mid-solve cancellation for now; typical solve duration
  is expected to be shorter than a user's reaction time. Revisit this when
  longer workloads or interactive cancellation make it valuable, and add an
  upstream `StiffAPI` operation that combines numeric rebinds with cooperative
  cancellation while preserving accepted-solution transactional semantics.
- [ ] **Deferred: changing the prepared solve interval.** Parameterized solves
  currently retain the time interval used during preparation. The initial
  batch workflow uses one common interval for every temperature condition,
  which is sufficient for the agreed scope. If a use case requires a different
  final time per condition (for example, a common final temperature under
  different heating rates), first add and validate an RST prepared-solver API
  for changing interval bounds without rebuilding the symbolic network. Until
  then, do not imply that each heating-rate run reaches the same final
  temperature; report the actual interval and resulting temperature range.

- Symbolic StiffAPI callbacks do not currently return KiThe's typed
  `SolidStateDomain` from intermediate trial states. KiThe validates initial
  and published states and its explicit RHS probe; preserve the native
  callback error until an upstream typed domain-callback hook exists.
- [x] RST 0.5.2 attempt reports provide a `solve_id` and cancellation
  disposition for all four StiffAPI families, with optional telemetry and
  diagnostic partial trajectory. KiThe consumes these reports whenever
  telemetry is enabled and retains their identity in
  `NetworkIvpAttemptDiagnostics`.
- [x] KiThe keeps the original allocation-light solve route for telemetry
  `Off`; it does not request an RST report or copy partial state in that mode.
- [x] Debug stories cover success followed by native callback failure and
  pre-start cancellation after a successful run. A failed/cancelled attempt
  cannot replace the accepted solution, and any partial trajectory remains a
  distinct diagnostic-only type with the current attempt ID.
- [ ] Native callback diagnostics still omit trial time/state and some native
  failure paths cannot identify the residual/Jacobian element or domain
  category. Preserve only RST-provided typed facts; do not infer them from
  display strings or substitute the last accepted state.
- [ ] A callback failure before the native solver records a failed status may
  have a valid attempt ID and typed error but no telemetry or partial state.
  Keep these fields absent rather than attaching the previous solve's snapshot;
  request upstream support if early callback failures need counters/state.
- [x] Repeated the wall-clock mid-run cancellation characterization on RST
  0.5.2 and recorded it in Story 25. Repeat it after future solver/API changes;
  LSODE2 may finish before the request, which is a reported disposition and not
  automatically a defect.
- Add reaction/expression role and trial-state context for native callback
  errors only when the evaluator boundary can supply that context.

### Safe Follow-Ups

- Extend the nonlinear scale sweep to strongly coupled networks if a concrete
  use case justifies the larger oracle and matrix. The current stiff fixture
  uses independent subsystems so its reference remains exact.
- Measure `NetworkProblem` validation/construction separately only if caller
  workloads show it is material. Current stage telemetry measures prepared
  assembly, backend preparation, integration, and output postprocessing.
- Add performance axes only when correctness-backed use cases justify them;
  the current release baseline does not itself call for optimization.
