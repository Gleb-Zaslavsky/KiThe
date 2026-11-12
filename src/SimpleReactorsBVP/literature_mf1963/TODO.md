# MF-1963 Independent Reference TODO

Scope: an independent SciPy reference for the first-order Merzhanov-Filonenko
flow-reactor model. This folder must not call KiThe's production BVP solver or
RST. The frozen profiles are numerical fixtures, not claims that KiThe has
already reproduced the model.

## MF-00: independent mathematical reference

- [x] Record the dimensional model, nondimensionalization, boundary conditions,
  and validity limits in `MF_1963_REFERENCE.md`.
- [x] Select weak/long-induction and intense/compact-front cases from the
  independently computed branch profiles, not arbitrary seed values.
- [x] Implement a standalone NumPy/SciPy `solve_bvp` generator with analytic
  source derivative, stable-manifold initialization, mesh/domain convergence
  studies, and deterministic CSV/JSON output; use Matplotlib only for plots.
- [x] Add local assertions for equation derivative, unit/parameter round trips,
  physical bounds, boundary conditions, and the integrated source identity.
- [x] Generate and inspect both frozen reference profiles and convergence
  metadata in the supported Python environment.
- [x] Record the generated solver versions, residuals, mesh sizes, continuous
  source maxima, integral checks, and low/high regime diagnostics in the
  reference document and parent BVP story log.
- [ ] Independently review every dimensional conversion before any later
  KiThe comparison; in particular verify `Q_mol = Q_mass * M` and the
  temperature-dependent density mapping.

## Follow-up, explicitly out of MF-00 scope

- [ ] A later MF-01 may map the frozen model into KiThe's `FastElemReact` and
  compare dimensional profiles. It must preserve the constant-property and
  first-order assumptions used by this reference.
- [ ] Do not add second-order kinetics, search for a critical ignition
  parameter, reconstruct digitized literature curves, or claim uniqueness in
  this stage.

## MF-01: KiThe integration verification

- [x] Adapt the dimensional inputs to the actual `SimpleReactorTask` API:
  `Q` is J/mol, `n` is the temperature exponent, density is fixed-composition
  ideal gas, and transport input must invert KiThe's `rhoD` normalization.
- [x] Add preparation assertions for `P/Tm/M`, recovered `rho(T0)`, `rhoD`,
  `Le`, `Pe_q=xi_max`, reaction topology/rate, and temperature/flux scaling.
- [x] Compare MF-HIGH through the full KiThe BVP_sci Sparse/ExprLegacy/Lambdify
  path from both the MF-00 profile seed and an independently constructed smooth
  seed.
- [x] Exercise both production routes, Damped with a nontrivial grid strategy
  and CorrectionNewton, against MF-HIGH and MF-LOW; retain independent physical
  profile acceptance rather than treating solver status as correctness.
- [x] Compare tuned Damped Pearson, Grcar-Smooke, and TwoPoint indicators
  against DoublePoints on both regimes using the same reference seed and
  refinement budget. None reaches the hot reference branch; see the full error
  table in the parent `STORY_TESTS.md`.
- [x] Record the current production-route outcomes and partial telemetry in the
  parent Story log. CorrectionNewton passes MF-HIGH; the other measured routes
  expose unresolved branch selection/convergence gaps.
- [ ] Determine why the MF-LOW BVP_sci long-domain/narrow-front solve fails;
  do not alter production equations or tolerances just to obtain PASS.
- [x] Recheck the MF-HIGH hot-root handoff after the BVP_api guess-forwarding
  fix in RST 0.4.34. The documented state-major seed now reaches and preserves
  the hot root; older cold-selection conclusions are superseded. MF-LOW and
  adaptive-refinement robustness remain open; see the 0.4.34 update below.
- [x] Compare Damped fixed-mesh behavior from the frozen reference seed, a
  smooth seed, and a BVP_sci hot-profile seed; independently verify the cold
  output against the MF trapezoid equations and first integral.
- [x] Test whether stale Jacobian reuse caused cold-branch selection by
  repeating the fixed-mesh BVP_sci hot-seed solve with `max_jac=1,2,5,50`;
  refreshing more often did not change the selected branch.
- [x] Repeat an MF-LOW fixed-mesh Damped solve with `max_jac=1` and `50`.
  Refreshing the Jacobian improves the independent MF defect and iteration
  count, but the returned profile remains on the cold branch.
- [x] Continue the MF-HIGH hot branch with BVP_sci as a control and separately
  continue Damped from its correctly published cold root. BVP_sci retains the
  hot branch through the tested `delta=4.143`; Damped tracks the cold branch
  from `delta=0.1` through `4.142`, then returns `NoResult` at `4.143`.
  Correctly packed hot seeds already select the cold root at `delta=0.1`.
  The measured Damped outcome bracket `(4.142,4.143)` is specific to this fixed
  mesh/configuration and is not evidence of a mathematical fold.
- [x] Continue the Damped cold root from `delta=3` with BVP_sci: the cold
  solution persists through `delta=4`, while the ordinary BVP_sci continuation
  fails at `delta=4.5`. This does not prove a fold; the continuation method may
  also have reached its convergence limit.
- [x] Correct the Damped initial-guess matrix construction. For mixed endpoint
  boundary conditions, each matrix row must contain that state's own free-node
  values; node-major flattening of a full profile is not a valid reduced seed.
  The earlier warm Damped continuation obtained from the incorrectly reshaped
  seed is withdrawn. See the correction appended to the parent
  `STORY_TESTS.md`.
- [x] Refine the Damped cold-continuation failure bracket from `(4.1,4.2)` to
  `(4.142,4.143)` with `0.001` parameter steps. Failure at `4.143` reports
  `NoResult`; one rerun reported 59 iterations, 346 damping trials, and 289
  rejections. Counters vary between runs.
- [x] Distinguish the independent discrete cold-branch fold from Damped's
  observed failure. Pseudo-arclength continuation finds a fold near
  `delta=4.142789` on 301/601/1201-interval trapezoid systems. Damped's
  `(4.142,4.143)` outcome bracket is consistent with that discrete fold, but
  remains configuration-specific and is not a continuous-BVP claim.
  Continuous collocation reaches `delta=4.1427` and exceeds its mesh-node cap
  at `4.1428`; the continuous fold remains unresolved.
- [x] Determine current Damped decision telemetry semantics. In RST 0.4.32
  Legacy mode, acceptance compares undamped correction norms
  (`s1 < 1 || s1 < s0`), not residual decrease. `DampingTrial` records alpha
  and NaN; rejection events retain the next alpha and correction norm, but
  accepted events and the complete comparison operands are absent.
- [ ] Request upstream Damped decision telemetry with typed acceptance and
  rejection details. Do not redefine the existing numeric event fields from
  KiThe.
- [ ] Find a CorrectionNewton configuration/initialization that reaches the
  MF-LOW reference, or document its measured lower capability boundary.
- [ ] Run and classify Frozen independently for both regimes; it is a
  comparison route, not one of the two production routes requested here.
- [x] Define independent dimensional-profile gates for temperature, species
  fractions, heat flux, finiteness, Kelvin positivity, and species closure; log
  the current route matrix in `STORY_TESTS.md`.
- [ ] Reconcile the remaining MF-01 requirements when the full instruction is
  available; the current pasted attachment truncates partway through section 8.

## MF-02B diagnostic completion (2026-11-10)

- [x] Verify the corrected hot-seed mapping and first-step audit.
- [x] Build and test an independent exact hot root for the 301-interval
  trapezoid system, map it to the six KiThe physical states, and pass that seed
  directly to Damped. Damped still publishes the cold root.
- [x] Run independent pseudo-arclength continuation and fixed-scaled sparse
  Jacobian spectra at 301, 601, and 1201 intervals.
- [x] Preserve machine-readable continuation data, plots, and the full
  diagnostic report in this folder.
- [x] Run the complete MF test module: 14 passed, 1 ignored.
- [ ] Upstream RST telemetry enhancement and a public/testable full production
  residual/Jacobian evaluation hook remain follow-up work. No production
  solver algorithm or acceptance criteria were changed in MF-02B.

## MF-02B.1 correction and MF-03 handoff (2026-11-10)

- [x] Correct the row layout of the independent trapezoid `F_delta`; test the
  full vector against centered finite differences on two grids and multiple
  nonuniform states. Relative errors are `3.56e-11` and `4.43e-11`.
- [x] Update pseudo-arclength tangent through the bordered system, refine the
  actual tangent sign-change interval to arc step `1e-5`, and record original
  plus extended-corrector residuals.
- [x] Repeat 301/601/1201 interval continuation. Corrected tangent-zero fold
  estimates: `4.142800580301762`, `4.142791345592427`, `4.142789031353778`;
  old coarse estimates differ by at most `9.2e-9`.
- [x] Regenerate JSON/CSV and plot artifacts; run the Python reference suite
  (`8 passed`).
- [ ] MF-03 production residual/Jacobian evaluation still requires a safe RST
  API. BVP_Damp's typed accepted-iteration trace is now available through
  `BvpSolver::damped_decision_trace_snapshot()`; do not duplicate callback
  assembly to obtain full-state residual/Jacobian values.

## RST 0.4.34 guess-forwarding regression (2026-11-11)

- [x] Convert the numerical-only Damped guess-packing diagnostic into an
  assertion gate using public `BVP_api`, ExprLegacy, and the exact 301-interval
  six-state trapezoid profile. The independent full-system residual is
  `4.10e-13`; the documented state-major guess reaches the same hot root with
  initial residual `4.10e-13` and `max(u0)=100`.
- [x] Replay the MF-1963 group after the facade correction. The MF-HIGH
  fixed-mesh, hot-seed, seed-layout, and Jacobian-refresh checks now select the
  hot branch and pass their physical-reference gates. The prior cold-branch
  reports and continuation conclusions based on the malformed facade handoff
  are superseded.
- [ ] Investigate MF-HIGH continuation accuracy: at `delta=0.125`, Damped stays
  on the hot branch but the independent `p` trapezoid defect is `1.445e-4`,
  above the existing `1e-5` gate. Do not relax the gate without numerical
  analysis.
- [ ] Investigate MF-LOW Damped `MissingResult`: the 2049-interval fixed-mesh
  run ended after 37 iterations, 282 damping trials, and 247 rejections.
- [ ] Investigate MF-LOW adaptive DoublePoints `MissingResult`: the 2-refinement
  run ended after 300 iterations, 2270 damping trials, and 1983 rejections.
- [x] Full debug MF group on RST 0.4.34: `11 passed, 3 failed, 1 ignored`;
  the three failures are the two MF-LOW cases and the MF-HIGH continuation
  residual gate above. The new state-major facade gate passes independently.

## MF-04 parity and robustness investigation (2026-11-11)

Keep new work split by question. Rust tests belong in topical child modules;
Python parity probes and shared numeric fixtures live beside, not inside, the
large KiThe integration test module. Every matrix row must run independently
so one failure cannot hide later strategy/backend outcomes.

### MF-04A Damped adaptive-grid matrix

- [x] Move the adaptive strategy matrix into `adaptive_refinement.rs`; every
  case/strategy row runs after earlier solve failures.
- [x] Run each MF regime x strategy combination independently: DoublePoints,
  Pearson `(1e-2, 1.4)`, Grcar-Smooke `(1e-2, 5e-2, 1.4)`, and TwoPoint
  `(2e-2, 2e-2, 1.4)`, with the existing two-refinement budget and common
  physical reference gates.
- [x] Record final node counts, refinement counts, solve status, profile
  errors, and independent defects. Do not treat additional nodes as evidence
  of improved accuracy; compare on a common dense reference mesh.
- [x] Ensure a failed row is reported but does not prevent the remaining rows
  from running. Keep the current strict numerical/physical acceptance values.
- [x] Refresh the MF-04 Story table after all eight rows completed.
- [x] Rerun the final eight-row matrix after modularization and expected-status
  assertions: `1 passed` in `394.93 s`; all four MF-LOW cases still fail to
  publish, while only Pearson/Grcar pass the MF-HIGH independent accuracy gate.
- [ ] Capture per-pass indicator values from upstream telemetry if exposed;
  current public snapshot reports refinement totals, not each strategy's raw
  indicator values.

### MF-04B SciPy / BVP_sci parity probes

- [x] Add a dedicated `scipy_parity` module, standalone Python solver, and a
  shared 17-digit two-state mesh/seed fixture. No side crate or test body in
  the large `mod.rs`.
- [x] Compare the MF-HIGH two-state numerical BVP (`theta'=p`,
  `p'=p-source(theta)`) using identical equations, BCs, mesh, seed, and solver
  tolerances in SciPy `solve_bvp` and RST `BVP_api::Sci` with ExprLegacy.
- [x] Use the same reference-derived hot seed generated once and serialized for both
  solvers. The existing Python reference does not use an arbitrary guess: it
  constructs a nonuniform seed by integrating from the stable manifold with
  DOP853. Preserve that distinction in reports.
- [x] Pair MF-LOW on both a 2049-node uniform mesh and a 2049-node clustered
  mesh, with identical equations, boundaries, tolerances, and a 6000-node cap
  in SciPy and RST. Both hit the node cap on the uniform fixture; both converge
  and agree within `3.84e-13` in theta and `1.99e-13` in p on the clustered
  fixture. This localizes the tested mismatch to mesh distribution/resolution,
  not an isolated SciPy-vs-RST discrepancy.
- [x] Extend the same standalone numerical comparison to the six-state
  MF-HIGH system. Assert the standalone RHS against KiThe's assembled RHS at a
  nontrivial state, then compare all six solved components with the frozen
  SciPy profile on shared coordinates.
- [x] Extend six-state parity to MF-LOW on the shared clustered 2049-node
  reference-derived seed and add a same-size uniform-mesh control. At tolerance
  `2e-7` and cap 6000, both implementations converge on the clustered mesh
  (SciPy 3817 nodes, RST 3985; max component delta `2.86e-8`) and both hit the
  node cap on the uniform mesh. Fixtures use production coordinate `z=xi/Pe`.
- [x] Probe whether the same nonuniform mesh can be passed directly to public
  Damped or CorrectionNewton. Both reject it before solve with typed
  `UnsupportedRoute` on the then-current RST release. This was a historical
  capability boundary, superseded by MF-05 / RST 0.4.35 explicit-mesh support.
- [x] Prepare a paste-ready RST request and a RST-public-API reproducer using
  the shared full six-state MF-LOW seed, successful clustered mesh, SciPy output
  fixture, and uniform control. The request specifies local-step residual/
  Jacobian handling and adaptive Damped preservation of the supplied mesh.
- [ ] After RST adds explicit-mesh support, run Damped and CorrectionNewton on
  the exact clustered fixture; compare both with the already-passing BVP_sci
  control and the same-size uniform mesh. (RST support and initial KiThe
  production transfer are now complete; MF-05 continues the broader matrix.)
- [x] Extend the matched meshes and reference-derived seeds to the full
  six-state system for MF-HIGH and MF-LOW; assert production RHS parity,
  boundary satisfaction, and componentwise SciPy/RST profiles.
- [ ] Add independent six-state Jacobian/ODE-defect checks before treating this
  parity fixture as a broader full-system accuracy gate.
- [ ] Add several separately named seed cases (stable-manifold, reference
  profile, and simple smooth guesses). Treat each as a basin experiment, not a
  universal arbitrary-initial-guess convergence promise.
- [x] A Rust/Python disagreement is actionable only after fixtures, ordering,
  boundary rows, coordinates, and numeric parameters are byte-for-byte
  matched. The paired MF-HIGH and MF-LOW fixtures now meet this condition.
  Keep both standalone reproducers self-contained and numerical.

### MF-04C Damped continuation accuracy

- [x] Isolate the `delta=0.1 -> 0.125` continuation pair in
  `continuation_accuracy.rs`.
  The matched 301-node run preserves the hot branch but reproduces the
  independent `p` defect `1.445325e-4` (gate `1e-5`), versus `6.509993e-9`
  at `delta=0.1`; profile change is `89.87534 K`.
- [ ] Continue localizing whether the defect is due to
  the stopping criterion, mesh resolution, or seed continuation. Preserve the
  current `p`-defect gate `1e-5`.
- [ ] Record the accepted solution, correction norm/threshold, available
  detailed Damped decision trace, raw independent trapezoid defects, and
  comparison to BVP_sci on the same mesh.
- [ ] Determine whether the `1.445e-4` defect comes from the Damped stopping
  criterion, discretization/state scaling, or continuation seed handling.
  Do not change solver tolerances or equations before this localization.

### MF-04D Remaining route boundaries

- [x] Rerun all four MF-LOW adaptive Damped strategies independently in
  MF-04A; all return `MissingResult`, with distinct work counts recorded in
  `STORY_TESTS.md`. Confirmed the public Damped route rejects arbitrary
  nonuniform input mesh before solving. Fixed-mesh and CorrectionNewton need a
  fresh post-0.4.34 side-by-side run.
- [ ] Classify Frozen/Naive separately on MF-LOW and MF-HIGH only after the
  paired baseline is stable. It is a comparison route, not a production
  correctness substitute.
- [ ] Update obsolete historical conclusions in a clearly marked current
  section; retain old measurements as history and never mix pre-0.4.34 seed
  forwarding results into current acceptance tables.

## MF-05 RST 0.4.35 solver-selection integration (2026-11-11)

Goal: transfer the new standalone RST MF-1963 solver-selection evidence into
KiThe in two distinct layers: numerical RST regression gates and physical
`SimpleReactorTask` acceptance stories. Do not treat the standalone RST matrix
as proof that the production physics facade has already been validated.

Primary RST inputs:

- `RustedSciThe/src/numerical/BVP_api/KITHE_MF1963_SOLVER_SELECTION.md`
- `RustedSciThe/src/numerical/BVP_api/mf1963_nonuniform_story_tests.rs`
- RustedSciThe version locked by this work: `0.4.35`.

### MF-05A Reconcile the current API and test baseline

- [x] Run a focused compile/test inventory for the MF folder against RST
  `0.4.35`; classify compile failures separately from behavior changes.
- [x] Audit direct `BVP_api` call sites and `SimpleReactorTask` handoff shapes:
  domain endpoints, interval/node count, mesh layout, initial-guess layout,
  mixed-boundary packing, and the Damped adaptive-mesh lifecycle.
- [x] Replace the old assumption that Damped and CorrectionNewton reject every
  explicit nonuniform mesh. The positive RST capability is now part of the
  expected contract; do not simply delete the old tests without an equivalent
  assertion.
- [x] Keep historical RST 0.4.32-0.4.34 observations as dated history. Mark
  the 0.4.35 results as the current baseline only after the focused gates pass.

### MF-05B Add the RST solver-selection regression module

- [x] Copy/adapt the upstream `mf1963_nonuniform_story_tests.rs` into a
  separately named KiThe MF diagnostic module, retaining its numerical-only
  character: public RST BVP API, symbolic equations, explicit mesh, and
  node-aligned guesses; no KiThe physics objects in this layer.
- [x] Bring over only the mesh/seed/reference fixtures required by those gates.
  Reuse byte-identical fixtures already present in this folder where their
  coordinates, ordering, and precision match; record fixture provenance.
- [x] Convert the upstream experiment outcomes into explicit regression
  expectations: successful configurations assert convergence and numerical
  checks; known unsuccessful configurations assert the typed failure/status
  that is intentionally being characterized. Never make expected failures
  look like successful physical solutions.
- [x] Cover MF-LOW and MF-HIGH in the RST-only regression module with the tested
  route/configuration set:
  BVP_sci/Sparse; Damped with TwoPoint, Pearson, Grcar-Smooke, and DoublePoints;
  CorrectionNewton/BandedAuto; Frozen default and measured refresh policies.
- [x] Keep seed provenance explicit: KiThe physical-profile seed, independent
  smooth seed, SciPy/reference seed, and solver-generated cross-seed are
  separate experiment classes. Do not call a solver-generated seed an
  ordinary caller guess.
- [x] Assert the appropriate numerical contract for each route. In particular,
  distinguish nonlinear/discrete convergence from agreement with the SciPy
  collocation profile; compare fixed-mesh trapezoid routes at matched meshes
  and report discretization error separately.
- [x] Preserve common-coordinate componentwise comparisons, endpoint checks,
  finiteness, relevant invariants, and independent ODE/trapezoid defect checks.
  Additional nodes alone are never an accuracy assertion.
- [x] Run each matrix row independently and print a compact result row so one
  route failure cannot suppress the remaining route/configuration evidence.

### MF-05C Extend the production `SimpleReactorTask` API for explicit meshes

- [x] Add a supported way for the production reactor workflow to receive an
  explicit finite, strictly increasing mesh in the solver's independent
  coordinate and a profile aligned to that
  mesh. Preserve the existing uniform `n_steps` path as the ordinary default.
- [x] Define and document one unambiguous KiThe-side layout contract for the
  explicit profile, including state/node ordering and whether boundary nodes
  are included. Validate mesh endpoints, dimensions, finite values, and
  agreement between supplied boundary values and the physical boundary
  conditions before invoking RST.
- [x] Carry the mesh/profile through the normal reactor preparation, BVP
  construction, solver handoff, transactional publication, and dimensional
  postprocessing. Do not create a test-only path that bypasses the production
  equation builder or postprocessor.
- [x] Map the MF-LOW clustered fixture's `z in [0,1]` coordinates directly to
  KiThe's configured solver coordinate. Do not apply `x=L*z` here: the
  production BVP is already formulated on normalized `z`; `L` enters the
  prepared coefficients and dimensional postprocessing. Retain all six state
  profiles and mixed endpoint conditions.
- [x] Verify Damped's semantics explicitly: the supplied mesh is its starting
  adaptive mesh, after which the solver owns refinement. Verify that
  CorrectionNewton and Frozen use the supplied mesh without silently
  resampling it. Verify BVP_sci behavior separately because it owns its
  collocation adaptation.
- [x] Add focused API contract tests for nonuniform mesh validation, aligned
  guess validation, boundary consistency, and preservation of the previous
  published solution/mesh if preparation or solve fails.
- [x] Keep the existing uniform-grid callers and task-document format
  unchanged unless the API audit demonstrates that a document-level mesh is
  required for the production use case. An internal/programmatic task option
  is preferable to expanding text syntax without need.

### MF-05D Build the physical HIGH/LOW solver matrix through KiThe

- [x] Using real `SimpleReactorTask` setup, execute the measured configurations
  through the explicit-mesh facade and physical postprocessing: clustered
  MF-LOW for BVP_sci, CorrectionNewton, Frozen refresh-on-correction, and all
  four Damped refinement policies; MF-HIGH for CorrectionNewton, all four
  Damped policies, and Frozen default/refresh-large/periodic-refresh; plus the
  same-size MF-LOW uniform BVP_sci node-cap control.
- [x] Use the RST report's measured settings as candidates, not universal
  defaults: Damped TwoPoint `(0.01,0.01,2.0)`, Pearson `(1e-3,1.4)`,
  Grcar-Smooke `(1e-3,1e-2,1.4)`, and DoublePoints; Frozen
  `at_high_morm=1e-6` and `every_m=1` as explicit alternatives; and the
  measured BandedAuto CorrectionNewton configuration. Revalidate every value
  through KiThe's adapter rather than assuming forwarding is complete.
- [x] Separate these outcomes in the initial production table: backend returned a
  solution; solver's own discrete stopping contract passed; independent
  physical/reference acceptance passed. A fixed-mesh Newton root may be a
  valid discrete solution while remaining less accurate than the refined
  collocation reference.
- [x] Check dimensional temperature, species fractions, heat flux, endpoint
  conditions, finiteness, species closure, and independent model residuals
  using the existing MF-00/MF-01 acceptance definitions. Do not relax gates to
  make a route green; if a route misses a gate, report which metric and retain
  it as an explicit limitation. The initial production stories assert finite
  dimensional profiles, Kelvin positivity, species closure, reference metrics,
  route status, and independent defects for the fixed/adaptive Newton roots.
- [ ] Include a small, intentional set of production-facade seed comparisons
  where RST findings depend on basin/seed provenance. RST-only cross-seed
  evidence exists, but the production LOW transfer currently uses the common
  clustered profile seed; keep caller and solver-generated seeds distinct.
- [x] State conclusions narrowly: characterize the tested model, mesh, seed,
  solver, options, and acceptance thresholds. Do not claim universal solver
  capability from these two regimes.

### MF-05E Remove only proven duplicate numerical reproducers

- [x] Compare every test in `six_state_parity.rs`,
  `rst_damped_hot_branch_reproducer.rs`, and `scipy_parity.rs` against the copied RST gate
  and the new production matrix.
- [x] Remove the obsolete rejection-only `nonuniform_mesh_api.rs` and
  `damped_mesh_control.rs`. Positive mesh acceptance is now covered by
  `rst_solver_selection.rs` and the production facade tests in
  `production_mesh_matrix.rs`.
- [x] Remove a numerical-only test/fixture only when the new RST gate covers
  the same equations, coordinates, boundary ordering, seed, configuration,
  and assertion. This audit found no additional exact duplicate: `six_state_parity.rs`
  owns KiThe RHS/standalone SciPy parity; `rst_damped_hot_branch_reproducer.rs`
  owns the state-major versus compact seed-packing regression; and
  `scipy_parity.rs` owns the separate reduced two-state mesh-distribution control.
- [x] Preserve or migrate unique KiThe guarantees, especially production RHS
  parity, dimensional mapping, full-profile handoff, mixed-boundary packing,
  physical postprocessing, and independent physical acceptance. Similar
  numbers from an RST-only reproducer are not a replacement for these checks.
- [x] Keep the two-state SciPy probe: it is an intentionally reduced numerical
  control, while the six-state fixtures test the complete KiThe-mapped system.
- [ ] Split or reduce oversized test modules only along clear ownership
  boundaries: RST numerical reproductions, KiThe physical mapping, and
  production solver acceptance. Avoid scattering one matrix across unrelated
  files.

### MF-05F Documentation and completion evidence

- [x] Add the current RST 0.4.35 production-facade snapshot to `STORY_TESTS.md` with the full
  story-test structure: name, exact commands, hypothesis/question, setup and
  options, observed numerical/physical outcomes, interpretation, and limits.
- [x] Distinguish current evidence from historical measurements in
  `MF_1963_KITHE_MAPPING.md` and `TODO.md`; do not overwrite earlier results
  that explain how the investigation evolved.
- [x] Update the MF-04 open items superseded by explicit-mesh support and the
  new parity evidence. Leave unrelated open items (for example continuation
  accuracy) open until their own assertion gates are resolved.
- [x] Run focused debug correctness gates and the complete MF ignored matrix
  serially. The full `literature_mf1963` debug run passed 50/50 tests, including
  ignored tests, after converting the two stale LOW unconditional-success
  assertions into typed/status characterizations.
- [ ] Run release correctness after the debug matrix is stable;
  performance benchmarking is a separate later decision and is not a gate for
  this correctness/API integration.
- [ ] Final review: every claimed route/configuration has a reproducible test;
  every removed reproducer names its replacement; test output clearly separates
  solver convergence from physical acceptance; no tolerances were weakened.

### MF-06 Production-facade parity with standalone solver-selection matrix

- [x] Compare compact standalone RST equations with the exact equation vector
  emitted by the KiThe production builder at every node of the common MF-LOW
  clustered seed. Record the maximum scaled discrepancy and verify it is
  expression-evaluation roundoff, not a missing equation or normalization.
- [x] Reconcile the apparently conflicting MF-LOW solver tables: the RST
  matrix records numerical convergence (and some rows' error against its own
  SciPy fixture); the KiThe facade matrix additionally applies the independent
  MF-00 physical-profile gate. Common Damped node counts and `theta` errors
  match after converting `theta` error to Kelvin.
- [x] Add an assertion diagnostic comparing the two separately generated
  MF-LOW SciPy references. They differ by maxima `1.177788894e-3 theta`,
  `1.008443983e-3 p`, and `1.177788894e-5 Y_A`; the theta delta explains the
  `0.0176668 K` BVP_sci error against the MF-00 reference, not a facade error.
- [x] For BVP_sci, CorrectionNewton, Damped TwoPoint/Pearson/Grcar-Smooke/
  DoublePoints, and Frozen, solve both the compact system and the exact
  KiThe-generated equations directly through public RST `BVP_api`; compare
  those results with the same configured routes through `SimpleReactorTask`.
- [x] Hold the 2049-node clustered mesh, six-state seed, solver options, and
  physical acceptance thresholds constant. Compare state vectors on shared
  input coordinates; do not require BVP_sci's adaptive final node count to be
  identical between equivalent expression forms.
- [x] Assert direct-RST/facade physical classifications match for all seven
  routes and facade handoff does not worsen the exact-equation direct-RST
  result. No acceptance tolerance was relaxed.
- [x] Record canonical-expression deltas separately from actual facade parity
  in `STORY_TESTS.md` MF-06. Exact-equation direct-RST/facade state discrepancy
  is at most `6.4e-12` scaled; all physical classifications match the
  standalone baseline.
- [ ] Run release correctness after the Debug parity gate; performance
  benchmarking remains a separate later step.
- [ ] Reconcile the independent SciPy reference artifacts
  `mf1963_mf_low.csv` and `mf1963_mf_low_six_state_scipy_solution.csv` at the
  generation/accuracy level. Keep both explicitly named until their measured
  drift is understood; do not silently substitute one oracle for the other.

### MF-07 Consolidated correctness decision

- [x] Give explicit binary answers for production-facade reference correctness
  and measured direct-RST parity, with their tested scope and caveats.
- [x] Consolidate HIGH/LOW solver settings, seed, numerical convergence,
  output mesh size, reference errors, and physical-gate status in
  `STORY_TESTS.md` MF-07. Include unsuccessful and RST-only diagnostic rows so
  convergence is not confused with physical acceptance.
- [x] Mark earlier MF mapping claims superseded where MF-05/MF-07 changed the
  observed LOW result; retain the old stories as dated diagnostic history.
- [x] Add a matched HIGH direct-RST-versus-production-facade parity test. The
  ignored Debug matrix covers nine MF-HIGH configurations using the exact
  KiThe-generated equations, the same 301-node mesh, caller seed, and options.
  Direct RST and facade agree to at most `1.6e-15` scaled state error and have
  identical MF-00 gate classifications; all nine pairs pass.
- [x] Sweep MF-LOW Damped tolerances and refinement thresholds for TwoPoint,
  Pearson, and Grcar-Smooke. Tightening nonlinear tolerance from `2e-7` to
  `1e-9` did not materially improve profile error. Stricter tested indicators
  improved each baseline, but none reached the physical gate; see MF-09 in
  `STORY_TESTS.md`.
- [ ] Investigate whether further adaptive-strategy tuning can approach the
  LOW physical gate without explosive mesh growth. The tested Pearson/Grcar
  threshold produced 31,792 nodes yet still had `6.90 K` maximum temperature
  error, versus the unchanged gate of `0.5 K`; uniform DoublePoints remains the
  only tested Damped route that passes LOW.
- [x] Isolate gradient-indicator coefficients (`g`) from value-jump
  coefficients (`d`) instead of changing both together. TwoPoint is unchanged
  when only `g` is reduced 100x; Grcar-Smooke adds nodes, but even `g=1e-4`
  yields only `11.35 K` error at 21,312 nodes. Pearson has no separate gradient
  coefficient: reducing its jump threshold to `5e-5` yields `4.59 K` at
  62,521 nodes. A trial at `d=1e-5` was stopped at about `12.6 GB` resident
  memory before completion. Details and commands are in MF-10.
- [ ] Do not continue lowering Pearson `d` blindly. If adaptive LOW work
  resumes, first inspect refinement-mark distribution/decision telemetry and
  compare interval placement against the high-gradient reference region; set
  an explicit node/memory budget before another sweep.

**Current production matrix:** MF-LOW physical-reference gates pass for BVP_sci
and Damped/DoublePoints with four refinements. CorrectionNewton,
Frozen/`at_high_morm=1e-6`, and Damped TwoPoint/Pearson/Grcar-Smooke return finite
discrete roots but miss the common physical gate. MF-HIGH passes for
CorrectionNewton, all four tested Damped strategies, and Frozen default,
`at_high_morm=1e-6`, and `every_m=1`. Exact evidence is recorded in
`STORY_TESTS.md` MF-05. Production cross-seed coverage and the complete release
correctness run remain open.

The prior unconditional LOW success tests have been converted to explicit
characterizations: uniform BVP_sci reaches the typed node cap at both 6000
nodes (production explicit-mesh control) and 30000 nodes (legacy convenience
path); fixed-mesh Damped on the uniform seed reports `NoResult` for
`max_jac=1/50`, while the adaptive clustered DoublePoints configuration passes
the physical gate. These expected outcomes remain distinct from the measured
HIGH/LOW solver-option matrix and do not relax acceptance.
