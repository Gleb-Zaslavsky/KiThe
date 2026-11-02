# ReactorsBVP Technical TODO

Scope: BVP reactor implementation in `src/ReactorsBVP/` plus its structured
GUI adapter in `src/gui/`. This checklist treats performance as important,
but correctness and route-contract clarity come first.

Status reconciliation (2026-11-02): the lockfile uses RustedSciThe 0.4.31.
Several recorded stories and the 2026-10-30 release campaign predate that
version; they remain valid historical evidence, but any behavior affected by
the RST update must be rechecked before calling it a current-version guarantee.

## P1: Unify reactor solves on RustedSciThe `BVP_api`

Scope decision: migrate the reactor's fixed-parameter BVP workflows to the
public `RustedSciThe::numerical::BVP_api` contract. Supported solver families
for this work are Damped, Frozen/Naive, and BVP_sci. `SciWithFreeParameters` is
not part of the reactor problem and is explicitly out of scope. BVP and IVP
work now proceed independently against the locked RST 0.4.31.

Starting point before this migration (RST 0.4.28):
- The parser-driven production route still builds `NRBVP` directly from
  `DampedSolverOptions`; the typed RST damped-settings parser is already used.
- `BVPSolver::solve_BVPsci()` is a second direct-native route with hard-coded
  ExprLegacy, Dense, and telemetry Off; it returns `()` and logs failures.
- `ReactorBvpSolverConfig` already owns symbolic frontend, matrix backend,
  Lambdify/AOT, and AOT build/execution settings. Extend this facade rather
  than replacing it.
- `BVP_api` already provides validated `BvpProblem`, typed Damped/Frozen/Sci
  configurations, a common solution view, route-tagged telemetry snapshots,
  typed source-preserving errors, and common logging-level forwarding.
- RST 0.4.29 exposes `BvpSolverSettings::with_telemetry_mode` and tagged
  native snapshots. Timings-only is intentionally rejected for Damped/Frozen
  configurations that may select AOT; BVP_sci supports all four modes.
- The reactor has one boundary condition per state in canonical state order.
  Damped/Frozen require uniform legacy meshes; Damped initial guesses are
  interval-major by columns while Frozen/Sci guesses are node-major by columns.
  Mesh and guess conversion must be explicit and tested.
- The GUI's `strategy = Damped/Frozen/Naive` is distinct from the explicit
  `solver_family = damped/frozen/sci` route selector. Bounds, adaptive-grid
  controls, and solver options must only be shown for families that accept
  them.
- `ReactorError` currently erases backend errors into strings in solve paths.

### Implementation plan

- [x] Define a capability matrix for Damped, Frozen/Naive, and BVP_sci.
  Record which common and family-specific settings each route supports. Reject
  unsupported combinations explicitly; never silently discard bounds,
  tolerances, execution, or backend selections. Keep `SciWithFreeParameters`
  excluded. See the RST 0.4.29 route matrix below.
- [x] Confirm upstream construction/configuration contract with RST.
  Ensure public typed APIs can produce each supported `BvpSolverConfig` from
  solver settings and expose route telemetry after both successful and failed
  solves. Add upstream work only for an actual missing capability; do not
  duplicate task-document parsing in KiThe.
- [x] Keep KiThe's generated-callback policy separate from RST solver-family
  selection. `ReactorBvpSolverConfig` owns symbolic assembly, matrix layout,
  and Lambdify/AOT generation for reactor callbacks; the hybrid parser and
  `BvpSolverSettings` own Damped/Frozen/Sci selection and each route's options.
  Telemetry is forwarded through RST's route-aware API, and execution policies
  remain where each backend can honor them. This avoids a duplicate KiThe
  solver-family config while retaining the Damped + AtomView + Banded + Lambdify
  defaults.
- [x] Add one reactor-to-`BvpProblem` adapter. Build equations/state order,
  independent variable, explicit mesh, and indexed endpoint conditions from
  the canonical reactor snapshot. Validate the boundary-condition count and
  retain the current physical boundary semantics.
- [x] Normalize initial guesses for each family. Test interval-vs-node column
  conventions, `n_steps` meaning, endpoint inclusion, and uniform mesh
  requirements. Use explicit errors for route-incompatible nonuniform grids.
- [x] Route direct and parser-driven solves through `BvpSolver`. Keep physics,
  reactor initial guesses, and postprocessing in KiThe; keep solver-settings
  parsing and native solver configuration in RST. Replace the hard-coded
  `solve_BVPsci` configuration with the selected typed route or remove the
  obsolete parallel entry point once all callers are migrated.
- [x] Make solve publication transactional. Publish solution, mesh, and
  argument metadata only after successful convergence and finite/dimension
  validation; preserve the previous successful state on a failed subsequent
  solve.
- [x] Preserve typed failures. Add a source-bearing `ReactorError` variant for
  `BvpApiError`, with reactor stage/context where useful. Do not stringify
  native errors. Retain the route-native telemetry snapshot with the latest
  solve outcome, including failure outcomes when RST exposes partial data.
- [x] Connect KiThe telemetry modes to each route's native telemetry without
  flattening unlike snapshot fields. Verify `Off` truly disables collection;
  clarify that logging is independent, and map the common parallel policy only
  to capabilities supported by the selected family/backend. GUI/task settings
  retain their route-native telemetry; direct `BVPSolver` callers can set an
  optional facade-wide override, which RST validates and applies before solve.
- [x] Update task-document and GUI contracts. The hybrid document continues to
  split physics/initial guess (KiThe) from solver settings (RST). The structured
  GUI now selects Damped/Frozen/Sci explicitly, conditionally renders each
  family's supported settings, removes incompatible fields on route changes,
  and round-trips Sci Lambdify/AOT settings through the RST hybrid parser.
  Rendered hover/help stories cover pre-existing Damped/Frozen controls,
  physical-section fields, adaptive-grid controls, Sci conditional banded
  controls, and Damped/Sci AOT options. The structured editor is scrollable so
  controls below the viewport remain reachable. Help assertions check expected
  content rather than merely requiring a nonempty panel.
- [x] Complete rendered hover/help checks for the tested BVP GUI controls.
  Stories cover route selectors and settings, legacy physical-section headings
  and representative inputs, adaptive-grid controls, conditional Sci banded
  controls, and Damped/Sci AOT controls. A bounded scroll area keeps lower
  controls reachable in a constrained viewport.
- [x] Expose optional Damped AOT residual/Jacobian target chunk counts in the
  structured GUI. RST 0.4.31 hybrid-parser round-trip, positive-count
  validation, AOT-only visibility, hover help, and compact-viewport reachability
  are covered in the GUI tests and `STORY_TESTS.md`.
- [x] Add focused correctness tests for problem/BC mapping, route-specific
  mesh/guess dimensions, scalar Damped/Sci solves and Frozen boundary rejection, typed setup
  failures, transactional state retention, telemetry modes, parser round-trip,
  and GUI family/options round-trip. Stories record command, hypothesis,
  observations, and conclusion.
- [ ] Re-run and finish the compact physical-reactor acceptance matrix on the
  locked RST 0.4.31. The latest detailed HMX story (recorded on 0.4.30) reports
  Frozen as typed `BvpApiError::MissingResult`/`NoResult`; Naive publishes a
  finite profile with prescribed boundary values reconstructed. That story
  does not establish Naive's interior-profile accuracy against an independent
  physical/discretization reference, so do not describe either route as fully
  HMX-accepted yet. Recheck both behaviors on 0.4.31 and preserve transactional
  publication/telemetry assertions. The analytic Frozen/Naive matrix across
  ExprLegacy/AtomView and Dense/Sparse/Banded is solver-control evidence, not a
  substitute for HMX acceptance. Older text claiming a Naive endpoint error of
  about `1.81e-1` conflicts with the later boundary-reconstruction story and is
  superseded; do not use it as the current failure description.
- [x] Update examples and API documentation to use the unified facade and show
  route-aware telemetry/error handling.
- [x] Add release benchmark coverage for Damped, Frozen, Naive, and BVP_sci
  routes, including small/medium/large analytic controls, stage attribution,
  route-native telemetry, and a separately gated AOT/Lambdify comparison where
  supported. Frozen/Naive AOT is measured only for AtomView Sparse/Banded, the
  available compiled routes, and the benchmark asserts profile parity before
  recording timings. Their physical HMX acceptance remains a separate open
  item. Implementation:
  `benches/reactor_bvp_baseline.rs`.
- [x] Run the release baseline campaign and record its current results in
  `STORY_TESTS.md`. Run `reactor-bvp-all-release-20261030-142037-511Z` passed
  all correctness steps and the 66-case Lambdify Criterion matrix. No
  performance gates are added until repeated measurements establish stable
  thresholds.
- [ ] Rerun the optional AOT benchmark axis with the corrected TCC probe
  (`tcc -v`); the completed campaign skipped AOT timings despite TCC being
  installed, while dedicated Damped/Sci AOT correctness stories passed. The
  probe fix is implemented, but no release AOT timing result is recorded yet.
- [ ] Record release measurements for the Sequential/Auto/Forced-Parallel axis
  and the optional SuperLU route. Debug dispatch-policy stories pass, but the
  latest story still has pending Criterion results. Only label a run SuperLU
  when native factorization diagnostics confirm zero fallback; otherwise mark
  the route unavailable rather than reporting a Faer measurement as SuperLU.

### Current RST capability status

The typed route/parser capability matrix below is implemented; its telemetry
and settings work originated in RST 0.4.29 and has since been extended. The
locked-version validation status is tracked separately above so historical
version labels are not mistaken for fresh 0.4.31 test evidence.

| Capability | Damped | Frozen / Naive | BVP_sci |
| --- | --- | --- | --- |
| `solver_family` | `damped` | `frozen` | `sci` |
| Nonlinear/collocation strategy | Damped Newton | Frozen or Naive | Fixed collocation route |
| Core numerics | scheme, matrix method, abs tolerance, iteration limit | scheme, matrix method, tolerance, iteration limit | execution, assembly, matrix layout, tolerance and collocation limits |
| Backend-specific options | Lambdify/AOT, symbolic frontend, banded solver, refinement, strategy/adaptive-grid settings | Lambdify/AOT, symbolic frontend, strategy parameters; no reactor bounds/adaptive-grid sections | Lambdify/AOT, ExprLegacy/AtomView, dense/sparse/banded layouts, execution and output policies |
| Reactor guess shape | interval-major; `n_steps` columns | node-major; `n_steps + 1` columns | node-major; `n_steps + 1` columns |
| Telemetry | Off/Counters/Detailed; Timings requires Lambdify-only backend policy | Off/Counters/Timings/Detailed | Off/Counters/Timings/Detailed, including AOT |
| Logging lifecycle | common text-log level through BVP facade; defaults Off | same as Damped | common text-log level through BVP facade; defaults Off |

- [x] Frozen and Naive analytic story matrix covers ExprLegacy/AtomView and
  Dense/Sparse/Banded through the KiThe `solve_BVPfrozen` facade. Compact-HMX
  physical acceptance is still open: the latest recorded RST 0.4.30 story has
  Frozen `NoResult`; Naive restores its endpoints, but interior-profile accuracy
  has no independent acceptance bound. Recheck on locked RST 0.4.31 above.
- [x] Release baseline includes Frozen and Naive routes with separate native
  telemetry rows; current single-host timings are recorded in `STORY_TESTS.md`.

The GUI exposes these route-specific settings without conflating the family
selector with the `strategy` choice. Damped AOT plus Timings is visibly called
out and rejected by the RST parser; Sci AOT Timings remains available because
that route accepts it.

- [x] RST 0.4.29 exposes telemetry mapping through `BvpSolverSettings` and the
  hybrid parser. The structured GUI exposes the four telemetry choices for all
  three families; route-specific unsupported combinations are rejected by the
  parser.
- [x] Direct programmatic `BVPSolver` callers can select the common telemetry
  mode through `with_bvp_telemetry_mode`; no second telemetry collector or
  flattened snapshot was introduced.
- [x] The structured BVP GUI exposes `telemetry` as Off/Counters/Timings/
  Detailed for Damped, Frozen, and Sci, defaults to Off, and keeps it
  independent from Rust logging; GUI stories verify route
  serialization and Sci parser acceptance.
- [x] The RST 0.4.29 settings-only parser produces typed Damped, Frozen, and
  Sci `BvpSolverSettings` without owning equations, mesh, boundary values, or
  initial guesses. KiThe filters the hybrid task document to solver-owned
  sections and attaches its own initial guess through `with_initial_guess`.
- [x] Legacy documents without `solver_family` remain accepted at the KiThe
  boundary: `strategy: Frozen|Naive` infers the Frozen family, otherwise the
  legacy strategy-only form maps to Damped. New task templates declare
  `solver_family` explicitly.
- [x] Migrate logging to the common `loglevel` setting and remove
  `dont_save_log` after its removal in RST 0.4.31. The level defaults to Off;
  telemetry and bounded BVP decision-event logging remain independent.
- [x] Verify Damped/Frozen native telemetry after unsuccessful outcomes and
  preserve the status and available partial counters through the KiThe facade.
  RustedSciThe 0.4.30 adds `BvpSolveStatus` (`NotStarted`, `InProgress`,
  `Converged`, `NoResult`, `Failed`) to snapshots; its native error path keeps
  accumulated Counters and last-attempt timings when enabled. KiThe's compact
  Frozen story checks `NoResult`, retained Counters, and transactional solution
  publication; its successful Naive retry checks `Converged`. RST's BVP_api
  failure matrix covers native `Failed` status and source-preserving errors.

Progress (2026-10-29):
- The current direct and parser-driven Damped solve paths now construct a
  `BvpProblem` and run through `BVP_api::BvpSolver`; the old direct NRBVP
  constructor remains only as a compatibility/test helper for now.
- The separate fixed-parameter `solve_BVPsci()` route now also uses
  `BVP_api::BvpSolver` with a typed `BvpSolverConfig::Sci`. Its former warning-
  and-return failure behavior is now `Result<_, ReactorError>`, and it retains
  Sci telemetry and publishes the adaptive mesh/profile only after validation.
  `SciWithFreeParameters` remains explicitly out of scope.
- A public `solve_BVPfrozen(initial_guess, FrozenSolverOptions)` entry point
  now uses `BvpSolverConfig::Frozen`; it infers node count from the initial
  guess, validates node-major conventions at the boundary, and retains the
  Frozen telemetry snapshot. The GUI's `strategy` field is deliberately not
  used as a proxy for this solver-family choice.
- `ReactorError` now retains `BvpApiError` and RST hybrid settings errors as
  error sources; `BVPSolver` stores the latest route-tagged telemetry snapshot.
- A scalar BVP story verifies a real Damped API solve, node/solution layout,
  mesh publication, telemetry capture, and preservation of the last successful
  result after a rejected retry. A separate story asserts that an RST
  configuration failure remains a `BvpApiError` source rather than a string;
  a scalar collocation story verifies the API Sci route, adaptive mesh,
  analytic endpoint value, and route telemetry; a Frozen story verifies its
  fixed node mesh, analytic endpoint, and telemetry variant. All are recorded
  in `STORY_TESTS.md`.
- RST 0.4.29's hybrid settings parser is now the parser-driven production
  route. Damped settings still require per-state `bounds` and `rel_tolerance`;
  missing or malformed fields are preserved as typed settings errors. Direct
  API defaults remain governed by the direct BVP API contract.
- A parser story now exercises Damped, Frozen, and Sci selection through the
  KiThe-to-RST adapter, verifies legacy family inference and the distinct
  interval/node initial-guess dimensions, and checks that a Sci settings map
  does not need reactor-owned physics sections.
- The structured GUI now exposes the RST telemetry mode for its current
  Damped, Frozen, and Sci forms, defaults it to Off, and preserves the
  independent logging policy. The family selector independently chooses the
  route; route-specific controls and bounds/adaptive sections are conditional.
  Damped AOT plus Timings is surfaced as unavailable and rejected before solve,
  while Sci AOT Timings round-trips successfully through RST 0.4.29. The
  existing nonlinear `strategy` remains a separate choice within Damped/Frozen.
- Direct programmatic solves now accept an optional common telemetry override
  on `BVPSolver`. It is applied through RST's public
  `BvpSolverSettings::with_telemetry_mode`; absent an override, route-native
  parser settings remain intact. RST's typed unsupported-mode error is retained
  as `ReactorError::BvpTelemetry`. The Sci Counters story checks nonzero callback
  counts and the native Sci snapshot variant.
- The reporting split is verified through data-level tests: task, equation,
  reaction-rate, balance, estimate, postprocessing, and render snapshots are
  available independently, and the existing print/export/plot methods consume
  those snapshots rather than requiring tests to parse console output.
- Idiomatic read-only accessors now cover thermal conductivity, mean molar
  mass, state boundary conditions, reaction heat, and dimensionless-temperature
  tolerance/bounds. Public legacy fields remain for source compatibility and
  direct new callers toward these accessors.
- The new GUI family story covers Damped-to-Sci/Frozen switching, Sci
  Lambdify/AOT serialization, removal of AOT-only values on return to Lambdify,
  route-specific bounds, and parser acceptance. The GUI suite passes 45/45;
  the focused RST hybrid settings parser story also passes.
- BVP GUI hover stories now exercise actual hover paths across Damped/Frozen/Sci,
  legacy physical sections, adaptive refinement, and conditional AOT/banded
  controls. The editor has a bounded vertical scroll area; `Maximum grid
  refinements` now forwards help from both its label and numeric widget.
- The compact HMX test fixtures now use a non-overlapping product label, and
  the synthetic A/B/C policy fixture defines A as a valid custom group while
  retaining B and C as element symbols. Production validation remains strict.
- The API adapter now publishes solution matrices in KiThe's established
  node-major layout (`mesh_points x state_count`), matching postprocessing,
  plotting, and balance checks. The previous state-major publication caused
  real parser/postprocessing and plot failures.
- Compact physical-reactor stories exercise all three API families. Against
  RustedSciThe 0.4.30, KiThe converts its public full nodal Frozen seed to the
  reduced `n_steps` unknown vector by omitting each state's prescribed endpoint.
  RST then publishes `n_steps + 1` mesh/profile values with the boundary values
  restored. The scalar `y'=-y`, `y(0)=1` story now confirms the first sample is
  exactly the x=0 boundary; `0.9961089494` was the first interior value, not a
  boundary violation. On compact HMX, Frozen returns typed `MissingResult`
  (`NoResult` telemetry) and preserves the previous Damped profile; Naive now
  publishes a profile with restored endpoints (`Converged` telemetry). These
  stories verify shape, boundary, and lifecycle contracts, not an all-node
  discretization error bound.
- A second compact inlet composition (`product_fraction = 0.02`, versus the
  baseline `0.001`) passes endpoint and finite-profile checks through both
  Damped and BVP_sci. This broadens measured compact-case evidence but is not a
  universal model-family guarantee.
- Current Debug `ReactorsBVP` lib run: 131 passed, 0 failed. This includes the
  parser solve/postprocessing paths and compact generated-backend matrix. AOT
  compiler-gated cases returned early because `KITHE_RUN_BVP_AOT_TESTS` was not
  set; they are not evidence that AOT compilers were exercised.

## P0: correctness, crashes, and user-facing failures

- [x] Replace parser `unwrap`/`expect`/`assert` paths with typed `ReactorError`.
  Functions:
  - `task_parser_reactor_BVP.rs::set_reactor_params_from_hashmap`
  - `task_parser_reactor_BVP.rs::parse_toleranse_and_bounds`
  - `task_parser_reactor_BVP.rs::parse_basic_settings`
  - `task_parser_reactor_BVP.rs::set_initial_guess_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  malformed task files currently can panic instead of returning a recoverable diagnostic. This is user-facing input, so panics are not acceptable.
  Done when:
  all required document fields are extracted through helper functions returning `Result`, tests cover missing section, wrong scalar type, wrong vector length, missing initial guess keys, and the parser accepts both vector and scalar-pair `bounds` encodings for compatibility.

- [x] Stop silently ignoring failed setup steps.
  Functions:
  - `task_parser_reactor_BVP.rs::set_reactor_params_from_hashmap`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  calls such as `let _ = self.fast_react_set(...)` and `let _ = self.setup_bvp()` can leave the reactor half-initialized and later fail far away from the real cause.
  Done when:
  these methods return `Result`, propagate errors with context, and tests assert that invalid reactions or invalid physical parameters fail at the parsing/setup boundary.

- [x] Validate finite positive physical quantities before derived calculations.
  Functions:
  - `SimpleReactorBVP.rs::check_task`
  - `SimpleReactorBVP.rs::mean_molar_mass`
  - `SimpleReactorBVP.rs::transport_coefficients`
  - `SimpleReactorBVP.rs::peclet_numbers`
  - `SimpleReactorBVP.rs::ideal_gas_density`
  Why:
  `inf`, `NaN`, zero denominators, or empty mass-fraction sums can produce meaningless solver inputs that look numerically valid.
  Done when:
  `M`, `P`, `Tm`, `Cp`, `Lambda`, `m`, `L`, diffusion coefficients, boundary fractions, and molar masses are checked for finite valid ranges before use.

- [x] Remove remaining panics from computational runtime paths.
  Functions:
  - `SimpleReactorBVP.rs::transport_coefficients`
  - `createBVP.rs::create_bvp_equations`
  - `createBVP.rs::set_solver_BC`
  - `SimpleReactorBVP3.rs::check_energy_balance`
  - `SimpleReactorBVP3.rs::check_material_balance`
  - `SimpleReactorBVP3.rs::postprocessing`
  Why:
  reactor setup, solving, and postprocessing should report domain errors instead of crashing the process.
  Done when:
  no non-test `panic!`, `unwrap`, or `expect` remains on data-dependent paths in `ReactorsBVP`.

- [x] Make solver execution report convergence failure explicitly.
  Functions:
  - `SimpleReactorBVP.rs::BVPSolver::solve_NRBVP_impl`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  `bvp.solve()` is called, then `get_result()` is accepted without checking whether the solver converged or produced a finite solution.
  Done when:
  solver status/result is validated, non-finite solution matrices are rejected, and a test covers failed/invalid solve output if the backend exposes status hooks.

- [x] Keep explicit postprocessing flags stable on the solved HMX regression path.
  Functions:
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  the solved BVP path should keep honoring `return_to_dimension`, `plot`, `gnuplot`, `save`, `save_to_csv`, and `no_plots_in_terminal` without regressing into panics or silent no-op behavior.
  Done when:
  a known solvable HMX task still solves with explicit postprocessing toggles and the regression test stays green.

## P1: hot-path performance and allocation pressure

Current status: the remaining copy pressure is now concentrated in owned solver/backend handoffs, not in avoidable intermediate cloning. For the current roadmap this pass is treated as complete, and the next practical work lives in P2 and P3.

- [x] Reduce large clones before BVP solver construction.
  Functions:
  - `SimpleReactorBVP.rs::BVPSolver::solve_NRBVP_impl`
  - `SimpleReactorBVP.rs::BVPSolver::solve_BVPsci`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  `eq_system`, `unknowns`, boundary-condition maps, parser data, and meshes are cloned at the solver handoff. For large mechanisms and grids this is expensive before the numerical work even starts.
  Done when:
  unavoidable backend ownership requirements are isolated in one handoff helper, avoidable clones are removed, and tests verify the handoff still preserves equation/order semantics.
  Current status:
  `BVPSolver::backend_boundary_conditions()` now centralizes the boundary-condition shape conversion, and both the direct solver path and parser-driven handoff paths reuse it. `solve_from_map()` and `solve_from_parsed()` now read parsed settings from borrowed maps instead of cloning them into a temporary parser snapshot, but the remaining `eq_system` and `unknowns` clones are still owned handoff copies.
  The parser wrappers for postprocessing and initial guess extraction also now borrow the parsed cache instead of cloning it first.
  The parser handoff also now moves the produced `x_mesh` into reactor state instead of cloning it back out of the backend snapshot.
  At this point the remaining copies are concentrated at the owned solver/backend boundary, so the item is mostly about respecting ownership rather than shaving off more avoidable clones.
  `SimpleReactorBVP3::check_energy_balance()` and `SimpleReactorBVP3::check_material_balance()` now borrow the solved state instead of cloning temporary snapshots, and `create_bvp_equations()` now borrows the thermal-effects slice and `qm` from the normalized context instead of copying them while the parallel assembly path avoids an intermediate result buffer. The transport snapshot is also now normalized into an ordered `d_ro_values` vector once per assembly, so the remaining copy cost is concentrated in the owned solver/backend boundary.
  `transport_coefficients()` now refreshes the cached `D_ro_map` in place, and `peclet_numbers()` borrows that cached map instead of keeping an extra local copy during the Peclet calculation.
  `set_solver_BC()` now builds boundary conditions from canonical `substances` and `unknowns` order, with `map_of_equations` checked as a contract instead of being used as the traversal source. That removes one more unordered read path from the hot setup flow.
  `set_reactor_params_from_hashmap()` now preallocates the main boundary, diffusion, substances, and group maps up front, so the parser setup avoids repeated reallocations while it materializes the task snapshot.

- [x] Introduce a reusable solver handoff builder.
  Functions:
  - `SimpleReactorBVP.rs::BVPSolver::solve_NRBVP_impl`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  there are multiple almost-identical NRBVP setup paths. Duplication increases the chance that one path gets stale tolerances, bounds, argument name, or boundary conditions.
  Done when:
  one helper builds the NRBVP object from a validated `BVPSolver` snapshot and both parser-driven solve paths use it.
  Current status:
  `BVPSolver::build_nrbvp_backend()` now constructs the backend snapshot from a single `NrbvpHandoffConfig`, and `task_parser_reactor_BVP.rs::solve_and_store_nrbvp()` reuses it.

- [x] Avoid repeated `HashMap` cloning in equation and boundary-condition creation.
  Functions:
  - `createBVP.rs::create_bvp_equations`
  - `createBVP.rs::set_solver_BC`
  - `SimpleReactorBVP.rs::transport_coefficients`
  - `SimpleReactorBVP.rs::Le_number`
  Why:
  these methods clone maps/vectors that are then only read. This is small for toy tests but expensive for large species lists.
  Done when:
  read-only paths borrow data, output maps are moved into fields where possible, and tests cover unchanged equation ordering and boundary-condition mapping.
  Current status:
  `create_bvp_equations()` now borrows the normalized transport and kinetic snapshot instead of rebuilding it per species, `set_solver_BC()` reads canonical solver order instead of iterating unordered maps, and the transport cache is refreshed in place before Lewis/Peclet calculations borrow it. That leaves no meaningful P1 work in this branch for the current scope.

- [x] Cache reusable symbolic and numeric vectors during equation assembly.
  Functions:
  - `createBVP.rs::create_bvp_equations`
  Why:
  symbolic expressions are cloned and simplified inside nested reaction/species loops. Large mechanisms amplify this cost.
  Done when:
  rate expressions, molar masses, powers, and transport coefficients are prepared once, simplification is performed only where it changes expression complexity materially, and existing equation tests remain green.

- [x] Introduce policy-based symbolic RHS assembly.
  Functions:
  - `createBVP.rs::create_bvp_equations`
  - `SimpleReactorBVP.rs::SimpleReactorTask::set_symbolic_rhs_assembly_policy`
  - `SimpleReactorBVP.rs::SimpleReactorTask::with_symbolic_rhs_assembly_policy`
  Why:
  the assembly path should stay sequential for tiny systems and switch to parallel work only when the problem is large enough to benefit from Rayon overhead.
  Done when:
  `Auto` selects a sensible threshold, explicit sequential/parallel modes are available, and tests confirm both paths build the same symbolic system.

- [x] Remove unconditional debug output from hot paths.
  Functions:
  - `SimpleReactorBVP.rs::BVPSolver::debug_solution`
  - `SimpleReactorBVP.rs::setup_bvp`
  - `SimpleReactorBVP.rs::mean_molar_mass`
  - `SimpleReactorBVP.rs::peclet_numbers`
  - `task_parser_reactor_BVP.rs::set_initial_guess_from_map`
  Why:
  terminal I/O can dominate runtime and make batch runs unreadable.
  Done when:
  debug printing is behind `log` levels or explicit flags, and solver/setup methods are quiet by default.
  Current status:
  `setup_bvp()`, `Le_number()`, `solve_BVPsci()`, and the main NR solve path are now quiet, and the remaining explicit diagnostics use `log` output rather than raw stdout. The old debug helpers still exist, but they are no longer unconditional console noise.

## P1: RustedSciThe damped NRBVP backend migration

Current status:
RustedSciThe `0.4.8` exposes the modern damped BVP API through `NRBVP::new_with_options(..., DampedSolverOptions)`. KiThe now routes reactor NRBVP construction through a local facade in `solver_backend.rs`, while old positional methods remain as compatibility wrappers for scheme, strategy, tolerances, and bounds. Reactor BVPs always start from symbolic equations, so the practical user-facing matrix choices are `Banded` and `Sparse`; dense and numeric-closure routes are intentionally out of scope for this module.

- [x] Map the new RustedSciThe NRBVP API.
  RustedSciThe types and methods identified:
  - `NRBVP::new_with_options(eq_system, initial_guess, values, arg, border_conditions, t0, t_end, n_steps, options)`
  - `DampedSolverOptions::{sparse_damped, banded_damped}` as the relevant reactor routes
  - `DampedSolverOptions::{with_abs_tolerance, with_rel_tolerance, with_max_iterations, with_bounds, with_strategy_params, with_logging_level}`
  - `DampedSolverOptions::{with_scheme, forward_derivative, trapezoid_derivative, with_scheme_name}`
  - `DampedSolverOptions::{with_symbolic_assembly_backend, with_sparse_generated_backend_defaults, with_banded_generated_backend_defaults, with_banded_lambdify}`
  - `DampedSolverOptions::{with_banded_atomview_c_tcc, with_banded_atomview_c_gcc, with_banded_atomview_zig, with_banded_atomview_for_repeated_solves}`
  - `GeneratedBackendConfig` and `AotBuildPolicy`, `AotExecutionPolicy`, `AotChunkingPolicy`
  - `BvpSymbolicAssemblyBackend::{ExprLegacy, AtomView}`
  - `MatrixBackend::{Banded, SparseCol}` as the relevant generated matrix families
  - `BackendSelectionPolicy::{LambdifyOnly, AotOnly, PreferAotThenLambdify}` as the relevant symbolic execution policies
  Why:
  this API is no longer just a Newton-strategy wrapper; it now controls residual/Jacobian execution, symbolic assembly, matrix representation, AOT build policy, AOT compiler backend, and runtime chunking.
  Done when:
  the TODO records which RustedSciThe knobs must be represented in KiThe before implementation starts.

- [x] Introduce a KiThe solver facade module for damped NRBVP options.
  File:
  - `src/ReactorsBVP/solver_backend.rs`
  Public/semipublic KiThe types:
  - `ReactorBvpSolverConfig`
  - `ReactorBvpExecutionBackend` with at least `Lambdify`, `Aot`
  - `ReactorBvpMatrixBackend` with `Sparse`, `Banded`
  - `ReactorBvpSymbolicBackend` with `ExprLegacy`, `AtomView`
  - `ReactorBvpAotCompiler` with `CGcc`, `CTcc`, `Zig`
  - `ReactorBvpAotConfig` with compiler, build policy/profile, compile preset, execution policy, and chunk targets
  Why:
  KiThe should not expose every RustedSciThe internal type directly in user-facing reactor code. A local facade lets us keep stable names, document combustion-oriented defaults, and absorb future RustedSciThe API shifts in one place.
  Current status:
  the facade builds `DampedSolverOptions` for `Banded/Sparse + Lambdify` and for `Banded/Sparse + AOT` with `tcc`, `gcc`, or `zig`. It keeps `Lambdify + AtomView + Banded` as the ordinary default and leaves actual AOT compilation to RustedSciThe at solve time.
  KiThe's task parser now delegates solver-side settings parsing to RustedSciThe's native BVP parser and only keeps the reactor-physics side, `initial_guess`, and a small compatibility bridge for legacy grid-refinement spellings plus scalar short-hand bounds/tolerances.

- [x] Make the new KiThe default `Lambdify + AtomView + Banded`.
  Target default mapping:
  - start from `DampedSolverOptions::banded_damped()`
  - force lambdify execution with `with_banded_lambdify()` or equivalent `GeneratedBackendConfig::banded_lambdify_defaults()`
  - keep `BvpSymbolicAssemblyBackend::AtomView`
  - keep the existing derivative scheme until explicitly changed by user config
  Why:
  this matches the requested production default for reactor BVPs while avoiding AOT toolchain dependence in ordinary runs and CI.
  Current status:
  `NrbvpHandoffConfig::new()` uses `ReactorBvpSolverConfig::default_lambdify()`, direct reactor solves and parser-driven solves pass through that default, and tests inspect the produced RustedSciThe options without running a heavy solve.

- [x] Replace `NrbvpHandoffConfig` positional solver fields with one options object.
  Functions:
  - `SimpleReactorBVP.rs::NrbvpHandoffConfig`
  - `SimpleReactorBVP.rs::BVPSolver::build_nrbvp_backend`
  - `SimpleReactorBVP.rs::BVPSolver::solve_NRBVP_impl`
  - `task_parser_reactor_BVP.rs::solve_and_store_nrbvp`
  Implemented change:
  `NrbvpHandoffConfig` stores the KiThe facade config, and `build_nrbvp_backend()` calls `NRBVP::new_with_options(...)`.
  Compatibility:
  keep existing `solve_NRBVP*` positional methods initially, but make them construct the facade/default config internally.
  Current status:
  direct and parser-driven solver paths both construct NRBVP through `new_with_options`, and the snapshot test checks old positional inputs plus the new generated backend default.

- [x] Extend task-document parsing for the lambdify solver backend settings.
  Functions:
  - `task_parser_reactor_BVP.rs::parse_solver_settings_from_map`
  - task template in `task_parser_reactor_BVP.rs`
  Supported document keys:
  - `generated_backend`: `banded_lambdify`, `sparse_lambdify`, `banded_aot`, `sparse_aot`, `banded_aot_tcc`, `banded_aot_gcc`, `banded_aot_zig`, `sparse_aot_tcc`, `sparse_aot_gcc`, `sparse_aot_zig`
  - `symbolic_backend`: `AtomView` or `ExprLegacy`
  - `matrix_backend`: `banded`, `sparse`
  - `aot_compiler` / `aot_c_compiler`: `tcc`, `gcc`, `zig`
  - `aot_build_policy`: `use_if_available`, `build_if_missing`, `require_prebuilt`, `rebuild_always`
  - `aot_build_profile`: `release`, `debug`
  - `aot_compile_preset`: `production`, `fast_build`, `dev_fastest`
  - `aot_execution_policy`: `auto`, `sequential`, `parallel`
  - `aot_target_chunks`, `aot_residual_target_chunks`, `aot_jacobian_target_chunks`
  Why:
  the parser is currently limited to old fields such as `method: Sparse` and cannot express the modern backend matrix.
  Current status:
  the parser maps absent backend keys to the new default, maps explicit lambdify/AOT backend keys into the facade, supports RustedSciThe-like AOT lifecycle keys, and validates positive chunk counts before constructing solver options.

- [x] Wire the AOT branch behind the facade.
  Target:
  - `ReactorBvpExecutionBackend::Aot`
  - `ReactorBvpAotConfig`
  - RustedSciThe helpers such as `with_banded_atomview_c_tcc`, `with_banded_atomview_c_gcc`, and `with_banded_atomview_zig`
  Why:
  AOT must not silently fall back to lambdify when the user explicitly requests compiled callbacks.
  Current status:
  `banded_aot_*` and `sparse_aot_*` produce real `DampedSolverOptions`; `tcc/gcc` fill the C compiler field, `zig` selects the Zig codegen backend, and parser tests cover build policy, execution policy, compile preset, compiler aliasing, and chunking. Ordinary tests do not execute a compiler.

- [x] Normalize Lambdify and AOT task documents into mutually exclusive execution contracts.
  Contract:
  - Lambdify documents use `backend_policy: lambdify_only` and carry no AOT compiler, build, or execution-lifecycle fields.
  - AOT documents use `backend_policy: aot_only`; matrix and compiler choices remain explicit.
  Compatibility:
  contradictory older documents are normalized before the native RustedSciThe parser sees them, while unknown backend names remain untouched so RustedSciThe can report a typed validation error.
  Current status:
  GUI controls, task templates, parser normalization, and story/regression tests all enforce this contract. The Lambdify parser-driven solve no longer triggers AOT DLL generation through a stale `build_if_missing` policy.

- [x] Add gated AOT solve/integration tests.
  Target:
  - compact combustion fixture using `AOT + tcc + AtomView + Banded`
  - optional `gcc`, `zig`, and sparse AOT comparisons when toolchains are present
  Policy:
  these tests must be ignored or gated by explicit environment/toolchain probes, because normal CI and user machines may not have `tcc`, `gcc`, or `zig` available.
  Done when:
  unavailable toolchains are skipped with a clear reason, available toolchains solve the compact task, and solutions are compared against the lambdify default with finite-profile assertions.
  Current status:
  `reactor_bvp_matrix_tests.rs` now contains gated AOT smoke tests for `tcc`, `gcc`, and `zig` on both banded and sparse matrix routes; they skip cleanly when `KITHE_RUN_BVP_AOT_TESTS` is not set or the toolchain is unavailable. The tests now compare numeric profiles against the reference backend instead of checking only shapes.

- [x] Move old combustion solve tests to the new default.
  Targets:
  - `simple_reactor_bvp_tests.rs` combustion/HMX tests that actually run NRBVP
  - `simple_reactor_bvp_tests2.rs` direct backend-construction tests
  Required default:
  `Lambdify + AtomView + Banded`.
  Done when:
  existing real combustion solve tests pass through the facade default and still validate finite solution, mesh, balance/postprocessing behavior, and stable solver state.
  Current status:
  `simple_reactor_bvp_tests.rs` now routes the HMX combustion smoke cases through the new default facade (`Banded` with lambdify execution and AtomView symbolic assembly), and `simple_reactor_bvp_tests2.rs` now builds the direct NRBVP smoke cases on the same default path as well.

- [x] Add backend matrix tests using the same compact combustion task.
  Test variants:
  - default: `Lambdify + AtomView + Banded`
  - compatibility baseline: `Lambdify + ExprLegacy + Sparse`
  - AOT leader: `AOT + C/tcc + AtomView + Banded`
  - optional AOT comparisons: `C/gcc`, `Zig`, and sparse AOT when toolchains are available
  Test policy:
  - normal CI tests should skip unavailable AOT toolchains gracefully or mark strict AOT tests as ignored
  - compare solutions against the default with tolerances on profiles, not exact matrices
  - always assert finite solution and mesh
  Done when:
  one reusable compact combustion fixture can exercise all configured backends without duplicating reactor setup.
  Current status:
  `reactor_bvp_story_tests.rs` and `reactor_bvp_matrix_tests.rs` now provide the reusable compact HMX fixture, backend-shape checks, solve smoke tests for the default and sparse lambdify routes, a compatibility baseline for `ExprLegacy + Sparse`, and gated AOT solve/profile-comparison tests across both banded and sparse matrix routes.

- [x] Decide how much RustedSciThe task parser syntax to mirror.
  Reference:
  RustedSciThe already parses keys such as `generated_backend`, `symbolic_backend`, `matrix_backend`, `aot_codegen_backend`, `aot_c_compiler`, `aot_build_policy`, `aot_build_profile`, `aot_compile_preset`, `aot_execution_policy`, `banded_linear_solver`, and `refinement_steps`.
  Decision:
  do not mirror the full solver syntax in KiThe; let RustedSciThe own solver keys directly and keep only the narrow reactor-side aliases that belong to physics compatibility.
  Boundary rule:
  KiThe keeps physics and reactor-specific initial-condition strategy, including `process_conditions`, `boundary_conditions`, and `initial_guess`.
  RustedSciThe owns solver-engine knobs such as tolerances, bounds, symbolic backend choice, matrix backend choice, AOT policy, compiler selection, chunking, solve execution settings, and native solver-side task parsing.
  Done when:
  task templates and examples use one canonical spelling for solver settings, solver aliases are tested natively in RustedSciThe, and KiThe no longer needs to duplicate backend-specific key names beyond the physics-side compatibility layer.

- [x] Finish the solver-settings handoff in `task_parser_reactor_BVP.rs`.
  Local KiThe-only responsibilities after the next RustedSciThe pass:
  - keep `set_reactor_params_from_hashmap`
  - keep `parse_physics_run_settings_from_map`
  - keep `set_initial_guess_from_map_data`
  - keep reactor-side post-solve output helpers only if they are still clearly reactor-specific
  Remaining local bridge:
  - `normalize_reactor_physics_task_map` only for physics-side shorthand cleanup (`grid_refinement` spelling and `C/J/Teta/q` family expansion)
  Removed local solver-settings code:
  - `parse_solver_tolerance_and_bounds_from_map`
  - `parse_solver_backend_config_from_map`
  - `apply_native_solver_settings`
  - the temporary solver-side compatibility branch inside `solve_and_store_nrbvp`
  Migration rule:
  move the typed solver settings into RustedSciThe first, then simplify KiThe solve paths to a single handoff call that only passes physics + initial guess + raw solver-settings document.
  Current status:
  KiThe now parses the solver contract through typed RustedSciThe entry points and no longer performs local solver-settings translation on the hot solve path:
  - `task_parser_bvp::parse_bvp_solver_settings_from_document()` for solver-selection and generated-backend settings
  - `task_parser_damped::parse_bvp_damped_solver_settings_from_document()` for damped solver scalars, bounds, and tolerances
  The remaining local bridge is now only about reactor-side compatibility normalization, especially `grid_refinement` spelling cleanup and the `C/J/Teta/q` family expansion into full reactor-variable maps.
  The removed local helpers are `parse_solver_tolerance_and_bounds_from_map`, `apply_native_solver_settings`, `parse_solver_backend_config_from_spec`, and the old solver-facade bridge inside `solve_and_store_nrbvp`, so the manual solver-settings parsing path is gone.
  Done when:
  KiThe solve path no longer rewrites solver settings locally, and the only remaining solver-facing code here is a thin bridge to the native RustedSciThe parser or builder.

- [x] Native RustedSciThe parser handoff for the KiThe client.
  Target files:
  - `task_parser_bvp.rs`
  - `task_parser_damped.rs`
  Desired shape:
  - expose `parse_*_from_document` split entry points for `problem` vs `solver_settings`
  - keep a typed `*_Spec`/`*_SettingsSpec` result instead of mutating `NRBVP` directly during parsing
  - add `try_*` methods that return `Result<..., BvpTaskError>` or a solver-specific typed error
  - keep old `set_*`/`parse_*` names only as thin compatibility wrappers if needed
  - remove `panic!` on bad user input; unknown backend/matrix/grid-refinement names must become typed errors
  - make `set_params_from_hashmap` and `set_postpocessing_from_hashmap` bridge-free and fallible in the typed path
  - keep solver settings parsing isolated from physical problem parsing
  Concrete solver settings to own natively:
  - `generated_backend`
  - `matrix_backend`
  - `symbolic_backend`
  - `aot_*`
  - `strategy_params`
  - `adaptive_strategy`
  - `grid_refinement`
  - `bounds`
  - `rel_tolerance`
  - `scheme`, `strategy`, `method`, `linear_sys_method`, `abs_tolerance`, `max_iterations`, and common `loglevel`
  Specific panic sites to remove:
  - unknown grid refinement method dispatch
  - missing / malformed solver settings sections
  - postprocessing parser fallthroughs that should be typed errors instead of process aborts
  Current status:
  RST 0.4.29 provides settings-only BVP and Damped parser entry points with
  typed settings/spec results and fallible `try_*` application paths. Its hybrid
  parser accepts solver-owned settings without KiThe equations, mesh, boundary
  conditions, or initial guesses. KiThe calls those public entry points and
  keeps only physics-side normalization plus the handoff of its own problem and
  initial guess. Older mutating parser methods remain RST compatibility APIs;
  they are not used as the typed KiThe production handoff.

## P1: numerical contract and validation coverage

- [x] Strengthen initial guess validation.
  Functions:
  - `reactor_BVP_utils.rs::InitialTemplate::generate`
  - `reactor_BVP_utils.rs::InitialConfig::generate_initial_guess`
  - `reactor_BVP_utils.rs::InitialConfig::only_one_value_for_all_initial`
  Why:
  `n_steps < 2`, zero variables, or non-finite template parameters can create invalid matrices or divide by zero.
  Done when:
  initial guess creation rejects invalid grid sizes and non-finite values, with focused tests.
  Current status:
  `InitialTemplate::generate()` now rejects `n_steps < 2` and non-finite template parameters, `generate_initial_guess()` rejects empty unknown sets and short grids, and `only_one_value_for_all_initial()` validates finite inputs and non-zero dimensions.

- [x] Validate concentration and boundary-condition semantics.
  Functions:
  - `SimpleReactorBVP.rs::check_task`
  - `SimpleReactorBVP.rs::mean_molar_mass`
  - `createBVP.rs::set_solver_BC`
  Why:
  boundary mass fractions are used as physical weights. Missing, negative, non-finite, or zero-sum fractions can invalidate density, molar mass, and source terms.
  Done when:
  tests cover valid normalized fractions, zero-sum fractions, negative fractions, missing species, and optional flux/heat-flux defaults.

- [x] Validate dimensions of kinetic artifacts before equation generation.
  Functions:
  - `createBVP.rs::create_bvp_equations`
  Why:
  the function already checks some lengths, but `vec_of_molmasses`, `stecheo_reags`, `stecheo_matrx`, `K_sym_vec`, `D_ro_map`, and substances must all agree before indexing begins.
  Done when:
  all index-sensitive structures are checked before nested loops and mismatch tests fail with `ReactorError::InvalidConfiguration`.

- [x] Add invariant tests for the full `setup_bvp` pipeline.
  Functions:
  - `SimpleReactorBVP.rs::setup_bvp`
  - `createBVP.rs::create_bvp_equations`
  - `createBVP.rs::set_solver_BC`
  Why:
  this is the central state transition from user task to solver-ready system.
  Done when:
  tests assert `unknowns.len() == eq_system.len() == BorderConditions.len() == 2*n + 2`, finite `M`, finite `Pe_q`, finite `Pe_D`, and stable variable ordering.

## P2: API, UX, and maintainability

- [x] Split parser extraction helpers out of `task_parser_reactor_BVP.rs`.
  Functions:
  - `set_reactor_params_from_hashmap`
  - `parse_toleranse_and_bounds`
  - `parse_basic_settings`
  - `set_postpocessing_from_hashmap`
  - `set_initial_guess_from_map`
  Why:
  the same nested `DocumentMap` extraction pattern is repeated many times and is easy to get wrong.
  Done when:
  small helpers such as `required_f64`, `optional_string`, `required_vector`, and `required_section` exist and are covered by parser tests.
  Current status:
  `required_section()`, `required_values()`, `required_f64()`, `required_usize()`, `required_string()`, `optional_string()`, `required_vector_f64()`, and `required_bool()` now centralize the recurring `DocumentMap` extraction patterns, so the parser no longer hand-rolls the same nested lookup logic in each caller.

- [x] Normalize naming without breaking public callers abruptly.
  Targets:
  - `Lambda`, `Diffusion`, `BorderConditions`, `Teta`, `Pe_D`, `Pe_q`, `Q`, `M`
  Why:
  mixed Rust/non-Rust naming makes the API harder to use and obscures which fields are user inputs versus derived solver state.
  Done when:
  idiomatic accessors or aliases exist, old public names are documented as compatibility fields or deprecated in a planned migration.
  Current status:
  scientific names that are standard in the domain, such as `Pe_D`, `Pe_q`, `Q`, and `M`, should stay as-is; the cleanup here is about API clarity and compatibility wrappers, not forcing the physics vocabulary into snake_case. Read-only accessors now expose `Diffusion`, `boundary_condition`, `D_ro_map`, `Pe_D`, `Pe_q`, and the scaling triple through a cleaner surface without breaking old callers.
  `scaling_config()` exposes the grouped `ScalingConfig`; `thermal_conductivity()`, `mean_molar_mass_kg_per_mol()`, `state_boundary_conditions()`, `heat_of_reaction()`, `dimensionless_temperature_tolerance()`, and `dimensionless_temperature_bounds()` cover the remaining legacy field names. Existing public fields remain available for compatibility and are documented as such where applicable.

- [x] Separate pure reporting from printing and plotting side effects.
  Functions:
  - `SimpleReactorBVP2.rs::pretty_print_task`
  - `SimpleReactorBVP2.rs::pretty_print_equations`
  - `SimpleReactorBVP2.rs::pretty_print_reaction_rates`
  - `SimpleReactorBVP3.rs::plot`
  - `SimpleReactorBVP3.rs::gnuplot`
  - `SimpleReactorBVP3.rs::plot_in_terminal`
  - `SimpleReactorBVP3.rs::save_to_file`
  - `SimpleReactorBVP3.rs::save_to_csv`
  Why:
  batch and library users need data/results without terminal output, GUI popups, or file writes.
  Done when:
  pure data-returning helpers exist, side-effect methods call them, and tests can inspect reports without printing.
  Current status:
  `task_report()`, `equation_report_rows()`, `reaction_rate_report_rows()`, `balance_report()`, `estimate_values_report()`, and `solution_render_data()` now expose the report data without touching stdout, so tests can assert on canonical rows, cached balance metrics, quick estimates, and owned solution snapshots directly. The side-effecting pretty-printers and export methods still exist as thin wrappers.

- [x] Convert postprocessing methods to return `Result`.
  Functions:
  - `SimpleReactorBVP3.rs::check_balances`
  - `SimpleReactorBVP3.rs::check_energy_balance`
  - `SimpleReactorBVP3.rs::check_material_balance`
  - `SimpleReactorBVP3.rs::postprocessing`
  - `SimpleReactorBVP3.rs::estimate_values`
  Why:
  these methods assume solution, mesh, molar masses, and variable names exist. Missing solver output should be a typed error.
  Done when:
  postprocessing can be called safely before and after solve, with tests for missing solution and missing mesh.
  Current status:
  `check_balances()`, `check_energy_balance()`, `check_material_balance()`, `postprocessing()`, and `estimate_values()` now return typed errors, and `postprocessing_report()` gives tests a pure snapshot of the transformed solution before the solver state is rewritten.

## P3: cleanup and documentation

- [x] Extend the BVP performance matrix with BVP_sci SuperLU, callback execution
  policies, and KiThe-side setup cost.
  Current status: the benchmark conditionally adds the SuperLU route only when
  a configured bridge library exists, then checks native factorization use and
  fallback counters. The largest analytic case now compares Sequential, Auto,
  and Forced Parallel callback policies for both Damped and BVP_sci. A separate
  40-species/39-reaction fixture times KiThe `setup_bvp()` before RST handoff
  under Sequential/Auto/Parallel symbolic RHS assembly and asserts equation
  parity. The callback-policy story passes in Debug; refreshed release timings
  are pending the next campaign. See `STORY_TESTS.md`.

- [x] Document expected units and dimensional assumptions in API docs.
  Targets:
  - `FastElemReact`
  - `SimpleReactorTask::set_parameters`
  - `SimpleReactorTask::set_transport_properties`
  - `SimpleReactorTask::set_operating_conditions`
  - `ScalingConfig`
  Why:
  wrong units can produce plausible but physically meaningless BVP solutions.
  Current status:
  `FastElemReact`, `SimpleReactorTask::set_parameters`, `SimpleReactorTask::set_transport_properties`, `SimpleReactorTask::set_operating_conditions`, and `ScalingConfig` now carry explicit unit-oriented docs, and the example guides show the expected dimensioned setup flow.

- [x] Fix typos and user-facing wording.
  Targets:
  - `parse_toleranse_and_bounds`
  - `set_postpocessing_from_hashmap`
  - `from_mass_fractions_to_molar_conentration`
  - `energy_balane_error_abs`
  - `energy_balane_error_rel`
  - printed text such as `numver`, `int the task`
  Why:
  small naming issues make the API and diagnostics look less reliable.
  Compatibility aliases with corrected spelling now exist for `parse_tolerance_and_bounds()`, `set_postprocessing_from_hashmap()`, and `from_mass_fractions_to_molar_concentration()`, while the legacy spellings remain as wrappers so older callers do not break immediately.
  Current status:
  The most visible user-facing wording issues in `SimpleReactorBVP3.rs`, `SimpleReactorBVP.rs`, and the parser docs were cleaned up, grid-refinement aliases now accept human-friendly spellings, and the remaining legacy spellings are explicitly marked as compatibility aliases.

- [x] Add example-guides for BVP workflows after P0/P1 stabilization.
  Suggested examples:
  - direct elementary reactions to `setup_bvp`
  - document parser to solver
  - custom tolerances and bounds
  - postprocessing without plotting
  Why:
  examples should reflect the safer API after parser and solver error handling are fixed.
  Current status:
  `examples/reactors_bvp_guides.md` now covers direct setup, parser-to-solver,
  custom tolerances and bounds, postprocessing snapshots, native route-tagged
  telemetry, and source-preserving solve/configuration errors.

## P4: Physical task-document lifecycle stories

- [x] Exercise an edit through an actual BVP Physics-tab reaction widget, then Save As, reopen, and
  solve. The story verifies the reaction parameter, typed solver snapshot, and unchanged source file.
- [x] Make malformed replacement loading transactional through the dirty-document confirmation
  path. The story verifies the active document, typed GUI snapshot, selected path, and dirty flag
  survive a parse failure.
- [x] Extend failed-load coverage to a document with an existing successful result and prove that
  result remains associated only with the preserved active document. The malformed-replacement
  story snapshots the successful plot and verifies it remains unchanged after the rejected load.
- [x] Add semantic document round-trip assertions for species/order, reactions, physical and
  boundary inputs, solver/backend settings, initial guess, and postprocessing. Include integral
  floats to guard the Float-versus-Integer serialization contract. The GUI lifecycle story compares
  the full typed solver snapshot and checks the persisted physical, thermal, boundary, and transport
  fields individually.
- [x] A rejected BVP rerun clears the previous plot before validation can fail, so stale output is
  not presented as the result of the edited task.
- [x] Test a numerical solve failure followed by correction and successful rerun. A one-iteration
  HMX run fails; restoring the iteration budget completes and publishes a finite Kelvin profile.
- [x] Compare the checked-in physical task through `solve_from_file` and the GUI route against the
  same physical acceptance envelope. The GUI story runs the same Save-As document by both routes
  and compares the complete postprocessed mesh and solution matrix.
- [x] Keep file lifecycle stories correctness-only; no file-I/O performance gates were added.
