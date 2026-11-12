# Reactor BVP Story Tests

## Story: Hybrid Task Settings Route to the Selected RST Family

**Name:** KiThe parses solver-only BVP settings through the RST 0.4.31 hybrid parser.

**Command:**

```powershell
cargo test --lib --no-default-features hybrid_bvp_parser_routes_damped_frozen_and_sci_settings -- --nocapture
```

**Question / hypothesis:** Can the hybrid task boundary delegate solver-owned
settings to RST for Damped, Frozen, and Sci without passing reactor physics or
initial guesses into that parser, while retaining compatibility with older
strategy-only task documents?

**What is checked:** Damped, Frozen, and Sci maps are parsed into their typed
`BvpSolverSettings` variants. The test checks Damped's interval-based initial
guess columns versus node-based Frozen/Sci columns, confirms legacy inference
without `solver_family`, confirms the common `loglevel` setting is owned by
the RST settings parser, and verifies the removed `dont_save_log` key is
rejected as an unknown setting.

**Observed result:** Pass on RustedSciThe 0.4.31. All three settings variants
are produced, the legacy Damped task is accepted, and the family-specific
initial-guess dimensions are `n_steps` for Damped and `n_steps + 1` for
Frozen/Sci.

**Conclusion:** The parser boundary now supports the three fixed-parameter
families without making RST own the physical problem. This verifies settings
selection and handoff shape, not full reactor solves for each task-document
family; the direct API solve stories cover numerical route behavior separately.

## Story: BVP GUI Telemetry Is Independent from Logging

**Name:** Structured BVP GUI exposes the RST telemetry mode and defaults to Off.

**Command:**

```powershell
cargo test --lib --no-default-features bvp_solver_telemetry_choice_is_serialized_without_enabling_logging -- --nocapture
```

**Question / hypothesis:** Does the BVP GUI let a user select RST's telemetry
collection mode without coupling that choice to Rust logging or log-file policy?

**What is checked:** The seeded task starts with `telemetry: off` and
`loglevel: off`. The GUI changes telemetry to `counters` without changing the
logging level.

**Observed result:** Pass. The selected telemetry value is serialized as
`counters`, while `loglevel` remains `off`.

**Conclusion:** The task-document path exposes RST telemetry with collection
disabled by default and keeps telemetry distinct from logging. Route-specific
unsupported combinations are still rejected by RST's typed settings parser.

## Story: BVP GUI Selects Solver Family

**Name:** Structured BVP GUI emits route-specific Damped, Frozen, and Sci settings.

**Command:**

```powershell
cargo test --lib --no-default-features bvp_sci_family_controls_round_trip_through_rst_settings_parser -- --nocapture
cargo test --lib --no-default-features bvp_damped_aot_timings_are_explained_and_rejected_by_rst -- --nocapture
```

**Question / hypothesis:** Can a user select a BVP solver family independently
of its internal strategy and produce a settings document accepted by the RST
0.4.29 hybrid parser, including switching between Sci Lambdify and AOT?

**What is checked:** The GUI switches from Damped to Sci and Frozen, renders
route-specific controls, omits Damped-only bounds on Frozen, parses Sci settings
for both execution modes, and removes AOT-only fields when returning to
Lambdify. Sci AOT accepts Timings; Damped AOT shows an unavailable-action
warning for Timings and its serialized settings are rejected by RST. Hovering
every Sci control and conditional AOT control reveals its corresponding help
content. The existing Damped route remains covered by the GUI settings tests.

**Observed result:** Pass. Sci Lambdify, Sci AOT (including Timings), and Frozen
settings parse into their corresponding `BvpSolverSettings` variants. Damped
AOT plus Timings is visibly identified as unsupported and rejected by RST; the
GUI cleans incompatible settings while changing solver routes.

**Conclusion:** Solver-family selection is now an explicit GUI contract rather
than an inference from nonlinear `strategy`. This verifies settings routing and
parser compatibility, not full end-to-end reactor solves for every family.

## Story: Damped BVP Uses the Shared API Transactionally

**Name:** Scalar Damped BVP solves through `BVP_api` and preserves the last accepted result.

**Command:**

```powershell
cargo test --lib --no-default-features story_scalar_bvp_solves_through_bvp_api_and_retains_native_telemetry -- --nocapture
```

**Question / hypothesis:** Does the reactor-facing handoff preserve the Damped
route's interval-count initial-guess convention while publishing the shared API
solution in KiThe's node-major layout? If a later request is rejected before
solve, does the previously accepted solution and mesh remain intact? Does the
direct facade's common Counters request reach the native Damped telemetry?

**What is checked:** A scalar symbolic problem is solved over 32 intervals.
The test checks a finite `33 x 1` profile, 33 published mesh points including
both domain endpoints, a Damped snapshot with nonzero residual calls under Counters,
and unchanged solution and mesh after an invalid one-interval retry.

**Observed result:** Pass. The API result is published as a `33 x 1` KiThe
solution with a 33-point mesh; Counters reaches the native snapshot, and the
rejected retry leaves solution and mesh unchanged.

**Conclusion:** The Damped adapter's mesh/solution handoff, direct telemetry
override, and pre-publication failure behavior are covered. This is an adapter
lifecycle story, not a physical reactor validation or a performance measurement.

## Story: Compact Reactor Damped Route Meets Boundary Conditions

**Name:** Compact HMX reactor publishes a finite node-major profile through the Damped BVP API.

**Command:**

```powershell
cargo test --lib --no-default-features story_compact_reactor_solves_through_damped_bvp_api_with_boundary_conditions -- --nocapture
```

**Question / hypothesis:** Does the production reactor setup map its generated
equations and endpoint conditions into the shared Damped API and publish a
finite profile with the requested physical boundary values?

**What is checked:** The compact HMX reactor is prepared through `setup_bvp`,
solved on 32 intervals, and checked for a `33 x state_count` node-major profile,
finite values, and endpoint residuals below `1e-5` for every state.

**Observed result:** Pass. The facade publishes all 33 mesh points and every
configured endpoint value meets the `1e-5` absolute residual bound.

**Conclusion:** This is an end-to-end physical-reactor correctness gate for the
Damped family. It is not a performance measurement and does not establish
Frozen or Sci behavior.

## Story: Compact Reactor BVP_sci Route Meets Boundary Conditions

**Name:** Compact HMX reactor publishes a finite adaptive profile through BVP_sci.

**Command:**

```powershell
cargo test --lib --no-default-features story_compact_reactor_solves_through_sci_bvp_api_with_boundary_conditions -- --nocapture
```

**Question / hypothesis:** Can BVP_sci solve the same production reactor model
through the shared API while preserving its adaptive mesh and endpoint values?

**What is checked:** The reactor is solved from a 32-node initial profile with
an adaptive node cap of 128. The test checks a finite node-major profile, an
adaptive mesh no smaller than the initial mesh, and endpoint residuals below
`1e-5` for every state.

**Observed result:** Pass. BVP_sci publishes a finite adaptive profile and all
configured endpoint values meet the `1e-5` absolute residual bound.

**Conclusion:** This is an end-to-end physical-reactor correctness gate for
the fixed-parameter Sci route; it does not cover Sci AOT or free parameters.

## Story: Compact Reactor Frozen and Naive Profile Publication

**Name:** Frozen and Naive publish full finite profiles with prescribed endpoint values.

**Command:**

```powershell
cargo test --lib --no-default-features story_compact_reactor_frozen_and_naive_publish_full_nodal_profiles -- --nocapture
```

**Question / hypothesis:** Can the compact production reactor publish usable
full nodal profiles through both Frozen and Naive, with telemetry and endpoint
values consistent with each route's result?

**What is checked:** The production Damped route first publishes the compact
HMX profile. Frozen and Naive each receive that 33-node profile as their
initial guess. Both routes publish finite 33-node profiles and retain
route-tagged Counters snapshots with `Converged` status. Every Frozen boundary
value is checked within `1e-5`; Naive's reconstructed endpoints are checked
against their prescribed values.

**Observed result (Debug, current RST):** Both Frozen and Naive publish finite
profiles with the expected shape, converged telemetry, and endpoint values.
This test checks the public profile and boundary contracts, not internal
discretization accuracy at every node.

**Conclusion:** Frozen's reduced internal seed and full public profile are
distinct representations. KiThe removes the prescribed endpoint from a
caller-provided full seed before handing it to RST; RST restores that endpoint
when publishing the full mesh/profile. `NoResult` remains a typed, transactional
failure, distinct from a converged candidate rejected by KiThe validation.

## Story: RST Configuration Errors Retain Their Source

**Name:** Missing required Damped bounds and relative tolerances remain a typed RST error.

**Command:**

```powershell
cargo test --lib --no-default-features story_damped_api_preserves_rst_configuration_error_source -- --nocapture
```

**Question / hypothesis:** When RST rejects an invalid Damped configuration,
does KiThe retain the `BvpApiError` in the error chain instead of flattening it
to a message string?

**What is checked:** A Damped request omits the required bounds and relative
tolerance maps. The test checks the `ReactorError::BvpApi` variant, a reachable
native error source, and that no solution is published.

**Observed result:** Pass. RST's configuration rejection is preserved as a
typed source and no result is stored.

**Conclusion:** The current Damped route has a source-preserving API error
boundary. Parser-settings errors and other solver families need their own
coverage as those routes are migrated.

## Story: Adaptive Collocation Uses the Shared API

**Name:** Scalar BVP_sci solve uses the common BVP facade and publishes its adaptive mesh.

**Command:**

```powershell
cargo test --lib --no-default-features story_scalar_sci_bvp_uses_shared_api_and_publishes_adaptive_mesh -- --nocapture
```

**Question / hypothesis:** Can KiThe route the fixed-parameter symbolic
collocation problem through `BVP_api::BvpSolver` while preserving the Sci
solver's adaptive mesh and solution values, while a direct programmatic caller
can request native telemetry through one route-independent facade setting?

**What is checked:** A scalar `y' = -y` problem with a left endpoint condition
is solved over `[0, 1]`. The test checks a finite published profile, mesh
endpoints and refinement-compatible shape, the analytic right-end value, a
Counters run with nonzero RHS calls, and a second explicit Off run with zero
RHS/Jacobian calls.

**Observed result:** Pass. The API publishes a finite profile on the adaptive
mesh, the right endpoint agrees with `exp(-1)` within `2e-3`. Counters reaches
the RST Sci snapshot with nonzero RHS calls; Off reports zero RHS and Jacobian
calls while still solving successfully.

**Conclusion:** The direct facade override feeds RST's existing telemetry
policy and retains its route-native snapshot; it does not create a second
collector. Parser-owned telemetry remains unchanged when no override is set.
This covers fixed-parameter Sci only, not `SciWithFreeParameters`, which is
outside the reactor problem, or GUI/task-document family selection.

## Story: Frozen Publishes a Boundary-Aligned Scalar Profile

**Name:** The reactor facade rejects Frozen output outside the configured endpoint tolerance.

**Command:**

```powershell
cargo test --lib --no-default-features story_scalar_frozen_bvp_publishes_boundary_aligned_profile -- --nocapture
```

**Question / hypothesis:** Does the KiThe facade translate a full nodal seed
to RST's reduced Frozen unknowns and publish the reconstructed profile with
the left endpoint aligned to `x=0`?

**What is checked:** The scalar `y' = -y` problem has `y(0) = 1` and is passed
to Frozen as a full profile on 257 nodes with tolerance `1e-6`. The test checks
the full node-major profile and mesh shape, `mesh[0] = 0`, the published value
`y(0) = 1`, finite values, and a Frozen Counters snapshot with `Converged`
status.

**Observed result:** Pass on RustedSciThe 0.4.30. The first public profile
sample is the prescribed boundary value at `x=0`; the former `0.9961089494`
sample is an interior point, not the left endpoint.

**Conclusion:** The prior mismatch was between the reduced Frozen unknown
vector and the full published mesh, not evidence that the solver violated the
left boundary condition. This story verifies the publication contract; it
does not claim a global error bound for the internal discretization.

## Story: BVP GUI Hover Coverage Across Solver Families

**Name:** Visible BVP controls show their own contextual help, including controls below the fold.

**Commands:**

```powershell
cargo test --lib --no-default-features bvp_damped_and_frozen_controls_render_content_aware_hover_help -- --nocapture
cargo test --lib --no-default-features bvp_damped_aot_timings_are_explained_and_rejected_by_rst -- --nocapture
cargo test --lib --no-default-features bvp_sci_family_controls_round_trip_through_rst_settings_parser -- --nocapture
```

**Question / hypothesis:** Do real pointer-hover paths display the corresponding help for the pre-existing BVP controls and the route-specific Damped, Frozen, Sci, and conditional AOT/banded controls?

**What is checked:** The tests hover solver-family selectors, core settings, telemetry/logging controls, symbolic and generated backends, AOT options, Sci matrix bandwidth controls, adaptive-grid settings, legacy section headers, and representative physical inputs. Help assertions check content rather than merely requiring a nonempty panel. `Maximum grid refinements` is tested through both its label and numeric control.

**Observed result:** Pass. The structured BVP editor now scrolls independently, so lower AOT controls can be reached and hovered in a constrained viewport. The maximum-refinement control also updates the contextual help from its registered description.

**Conclusion:** The tested BVP GUI routes expose contextual help on their rendered controls and conditional states; this is not inferred from static help-map entries alone.

## Story: Compact Reactor Inlet Composition Across Damped and Sci

**Name:** Damped and BVP_sci satisfy the endpoint contract for a second compact inlet composition.

**Command:**

```powershell
cargo test --lib --no-default-features story_compact_reactor_inlet_composition_is_accepted_by_damped_and_sci -- --nocapture
```

**Question / hypothesis:** Do both accepted API families retain physical endpoint values when the compact HMX fixture's inlet product fraction changes from `0.001` to `0.02`?

**What is checked:** The alternate fixture is prepared through the production reactor setup, then solved independently through Damped and BVP_sci. Each published profile must be finite, have a mesh-consistent node-major shape, and satisfy every configured endpoint value within `1e-5`.

**Observed result:** Pass for Damped and BVP_sci. Both routes publish finite profiles satisfying the endpoint checks for the alternate composition.

**Conclusion:** The compact correctness evidence now covers two inlet compositions for the supported Damped and Sci routes. This does not close Frozen/Naive acceptance or establish a universal physical-model guarantee.

## Story: BVP Symbolic Frontend and Matrix Layout Correctness Matrix

**Name:** Compact HMX profiles remain equivalent across the supported Damped and BVP_sci numerical routes.

**Commands:**

```powershell
cargo test --release --lib --no-default-features matrix_compact_hmx_damped_frontend_and_matrix_matrix_matches_reference -- --nocapture
cargo test --release --lib --no-default-features matrix_compact_hmx_sci_frontend_and_layout_matrix_matches_boundaries -- --nocapture
cargo test --release --lib --no-default-features matrix_compact_hmx_large_grid_remains_correct_on_banded_and_sparse_routes -- --nocapture
cargo test --release --lib --no-default-features matrix_compact_hmx_sci_large_grid_remains_bounded -- --nocapture
```

**Question / hypothesis:** For a real reactor-generated BVP, do ExprLegacy and AtomView preserve the published profile and endpoint contract across every supported matrix layout, including a substantially refined mesh?

**What is checked:** Damped compares Banded and Sparse with both symbolic frontends against its reference profile. BVP_sci compares Dense, Sparse, and BandedAuto with both frontends and verifies the actual boundary conditions. Separate 96-interval runs check profile dimensions, finite values, and cross-layout equivalence/bounded output. Frozen/Naive's real-HMX endpoint limitation remains separately documented; the analytic solver-family matrix below now exercises these API routes independently.

**Observed result:** The Debug stories and the 2026-11-08 Release campaign pass for the 32-interval route matrices and both 96-interval cases. The exact run IDs and combined route/performance evidence are in `Story: Release BVP Correctness and Performance Campaign` below.

**Conclusion:** This is a physical-model correctness matrix, not a performance claim. The larger grid provides a second size point without substituting a synthetic equation system for the reactor fixture.

## Story: BVP AOT and Lambdify Produce Equivalent Reactor Profiles

**Name:** TCC-generated BVP callbacks preserve the HMX profile for matched frontend and matrix routes.

**Commands:**

```powershell
$env:KITHE_RUN_BVP_AOT_TESTS = '1'
cargo test --release --lib --no-default-features matrix_compact_hmx_sci_aot_matches_lambdify_across_frontends_and_layouts -- --nocapture
cargo test --release --lib --no-default-features matrix_compact_hmx_damped_aot_matches_lambdify_across_frontends_and_matrices -- --nocapture
```

**Question / hypothesis:** Does AOT execution agree with Lambdify when solver family, symbolic frontend, and matrix layout are held fixed?

**What is checked:** The Damped matrix covers two frontends by Banded/Sparse; BVP_sci covers two frontends by Dense/Sparse/BandedAuto. Each generated profile is compared componentwise with its matched Lambdify profile. Frozen and Naive AOT are included in the release benchmark only for their supported AtomView + Sparse/Banded routes; the benchmark asserts profile parity before recording timings. These gates require TCC and are deliberately opt-in.

**Observed result:** Both targeted HMX AOT tests pass locally with TCC; all ten matched Damped/Sci route pairs satisfy the configured profile tolerances. Frozen/Naive AOT parity is pending the opt-in expanded benchmark run.

**Conclusion:** AOT is checked as an execution axis, not conflated with a solver or frontend change. The correctness tests do not imply a speedup; release timings come from `reactor_bvp_baseline`.

## Story: BVP Solver Performance Baseline Across Problem Sizes

**Name:** Supported BVP families expose comparable wall-time, stage, and native telemetry records from small to large systems.

**Command:**

```powershell
cargo bench --bench reactor_bvp_baseline --no-default-features -j 2 -- --noplot
```

Optional TCC axis:

```powershell
$env:KITHE_RUN_BVP_AOT_BENCH = '1'
cargo bench --bench reactor_bvp_baseline --no-default-features -j 2 -- --noplot
```

**Question / hypothesis:** How do Damped, Frozen, Naive, and BVP_sci costs change from a one-state, coarse problem through a 24-state, 128-interval problem, and where do preparation, validation, native solve, and backend callback/factorization costs differ?

**What is checked/measured:** A deterministic analytic family `y_i' = -k_i y_i`, `y_i(0)=1` is used at 1/24, 8/64, and 24/128 state/interval sizes. The fixed route matrix is Damped Banded/Sparse, Frozen and Naive Dense/Sparse/Banded, and BVP_sci Dense/Sparse/BandedAuto, each with ExprLegacy/AtomView. The benchmark prints separate macro, KiThe/API stage, and route-native telemetry tables; Criterion supplies release wall-time estimates. AOT-TCC is a separately matched optional axis: Damped/Sci on their supported routes and Frozen/Naive on AtomView Sparse/Banded. The benchmark asserts AOT/Lambdify profile parity before recording comparison timings. The analytic fixture is a numerical performance control; physical validity is covered by the HMX stories above.

**Observed result:** Release campaign `reactor-bvp-all-release-20261030-142037-511Z` passed all six steps in `6284.471 s` with build jobs `1` (revision `92ff18d5e06381c3c34daeb9106f0a4e1430fe85`, dirty working tree). Correctness logs report 11 physical stories, 15 HMX matrix stories, one Frozen/Naive analytic matrix, and one story for each Damped and Sci AOT route: all passed. Criterion completed all 66 Lambdify scenario/routes. The optional AOT benchmark axis was unexpectedly skipped: the benchmark probed `tcc --version`, while the installed TinyCC accepts `tcc -v`. The detector is fixed; AOT benchmark parity/timings still need a follow-up run. Full campaign dashboard: `test_reports/reactor-bvp-all-release-20261030-142037-511Z/reports/summary.md`.

Large analytic case, Criterion medians (`24 states / 128 intervals`, ms):

| Family | Frontend | Layout | Median ms |
| --- | --- | --- | ---: |
| Damped | ExprLegacy | Banded | 90.00 |
| Damped | ExprLegacy | Sparse | 85.68 |
| Damped | AtomView | Banded | 92.16 |
| Damped | AtomView | Sparse | 86.75 |
| BVP_sci | ExprLegacy | Dense | 6.84 |
| BVP_sci | ExprLegacy | Sparse | 0.32 |
| BVP_sci | ExprLegacy | BandedAuto | 0.35 |
| BVP_sci | AtomView | Dense | 7.13 |
| BVP_sci | AtomView | Sparse | 0.41 |
| BVP_sci | AtomView | BandedAuto | 0.55 |
| Frozen | ExprLegacy | Dense | 560.11 |
| Frozen | ExprLegacy | Sparse | 46.90 |
| Frozen | ExprLegacy | Banded | 35.66 |
| Frozen | AtomView | Dense | 253.35 |
| Frozen | AtomView | Sparse | 34.19 |
| Frozen | AtomView | Banded | 37.57 |
| Naive | ExprLegacy | Dense | 581.01 |
| Naive | ExprLegacy | Sparse | 50.12 |
| Naive | ExprLegacy | Banded | 37.31 |
| Naive | AtomView | Dense | 316.54 |
| Naive | AtomView | Sparse | 35.16 |
| Naive | AtomView | Banded | 38.92 |

**Conclusion:** For this large analytic control, Frozen/Naive Dense is several times slower than their Sparse/Banded routes; Sci Dense is also much slower than Sparse/BandedAuto. AtomView notably reduces Dense Frozen/Naive time here, but is not uniformly faster for every family/layout. Treat these as single-host baselines, not performance gates. Native telemetry fields have route-specific timing semantics and must be compared within a solver family, not summed across families. No AOT speed conclusion is available from this campaign.

## Story: Frozen and Naive Frontend/Layout Correctness Matrix

**Name:** Frozen and Naive publish the same bounded profile across both symbolic assemblies and all supported matrix forms.

**Command:**

```powershell
cargo test --lib --no-default-features frozen_and_naive_routes_match_across_frontends_and_matrix_layouts -- --nocapture --test-threads=1
```

**Question / hypothesis:** After the RST Frozen/Naive repair and AtomView conversion fix, do both
strategies solve a bounded BVP through the public KiThe reactor facade for ExprLegacy and AtomView with Dense,
Sparse, and Banded matrices?

**What is checked:** The analytic `y' = -0.7 y`, `y(0)=1` profile is solved for all 12 strategy,
frontend, and layout combinations. The initial guess is deliberately perturbed from the analytic
solution. Each route must publish finite values in `[0,1]`, stay within `4e-2` of the analytic
profile, retain Frozen-family telemetry, and match the other routes within the same profile envelope.

**Observed result:** Pass in Debug on the current RST lockfile through the KiThe facade. All 12
strategy/frontend/layout combinations published matching profiles with zero measured route delta.
This analytic control is deliberately separate from compact HMX. The earlier
compact Frozen `NoResult` observation is historical; the current compact story
now asserts successful Frozen and Naive profile publication.

## Story: Largest-Case Callback Policy and KiThe Preparation Attribution

**Name:** The largest BVP workload preserves its solution across callback policies, reports actual dispatch choices, and exposes KiThe setup time before RST handoff.

**Commands:**

```powershell
cargo test --lib --no-default-features story_large_bvp_callback_policies_match_and_report_dispatches -- --nocapture --test-threads=1
cargo test --lib --no-default-features story_kithe_pre_rst_symbolic_policies_preserve_large_equation_system -- --nocapture --test-threads=1
cargo test --lib --no-default-features matrix_compact_hmx_sci_superlu_matches_faer_without_fallback -- --nocapture --test-threads=1
cargo bench --bench reactor_bvp_baseline --no-default-features -j 1 -- --noplot
```

Optional SuperLU route: set `SUPERLU_LIBRARY` to the RST SuperLU bridge library before running the campaign. The route is otherwise omitted, not silently represented as a Faer result.

**Question / hypothesis:** For the largest current analytic matrix (24 states, 128 intervals), do Damped and BVP_sci produce the same profile under Sequential, Auto, and Forced Parallel callback policies? Does native telemetry confirm that the selected policy was actually used? Separately, how much time does KiThe spend preparing a large reactor task before passing equations to RST?

**What is checked/measured:** The correctness story runs both solver families with AtomView/Sparse and all three callback policies. It compares the full published profile and requires at least one parallel dispatch for Forced Parallel. The detailed benchmark prints dispatch counters and callback timing alongside solve time. A 40-species/39-reaction setup fixture measures `SimpleReactorTask::setup_bvp()` under KiThe's Sequential/Auto/Parallel symbolic RHS assembly policies; the generated equation strings must be identical across policies. Criterion records the same preparation-policy axis independently.

**Observed result:** Debug correctness story passes. For Damped, Sequential recorded `0 parallel / 5 sequential` dispatches, Auto and Forced Parallel each recorded `5 / 0`. For BVP_sci, Sequential recorded `0 / 513`, Auto selected sequential (`0 / 513`), and Forced Parallel recorded `513 / 0`. All six runs published matching profiles within `2e-5`. The Release benchmark now records callback-policy timings and KiThe setup cost; see the 2026-11-08 campaign report below. `Auto` is a policy decision, not a synonym for forced parallel, and differs by solver family.

**Conclusion:** Keep Sequential, Auto, and Forced Parallel as distinct measurement axes. Do not infer whole-solver parallelism from callback dispatch counters: this comparison isolates residual/Jacobian callback execution. SuperLU must be verified through native factorization counters and zero fallback count before its timing is treated as a SuperLU measurement.

**Conclusion:** The repaired routes are eligible for the general BVP solver matrix and baseline;
this does not yet establish physical HMX acceptance for Frozen/Naive.

## Story: Checked-In Physical BVP Task Documents

**Name:** Readable HMX BVP task files load and solve through the production file route.

**Command:**

```powershell
cargo test --lib --no-default-features checked_in_physical_bvp_documents_load_and_solve -- --nocapture --test-threads=1
```

**Question / hypothesis:** Do checked-in, user-readable BVP documents preserve their physical
composition, kinetics, and solver settings when parsed by the standard `solve_from_file` path?

**What is checked:** Both compact HMX documents are solved from disk. The story checks finite
profiles, physical temperature in `0..3000 K`, and all concentration profiles in `0..1`. The pair
uses the same HMX decomposition model with distinct inlet composition, exercising more than one
physical input while keeping the fixture small and reviewable.

**Observed result:** Debug story passes for both files. The compact HMX profile spans approximately
`600..868.98 K`; the alternate-inlet profile spans approximately `600..863.85 K`. Each contains
33 points. File I/O is correctness coverage only, not a performance gate.

## Story: BVP GUI File Lifecycle

**Name:** A physical BVP task can be edited, saved, reopened, and solved in the GUI.

**Command:**

```powershell
cargo test --lib --no-default-features physical_bvp_task_doc_gui_load_edit_save_reload_and_run -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the structured BVP GUI preserve physical inputs through an
edit/save/reopen cycle and run the reopened document through its normal calculation action?

**What is checked:** The test opens the alternate-inlet HMX task, edits the reaction pre-exponential
factor through its Physics-tab widget, saves to a temporary task file, and reopens it. It verifies
species order, reaction data, thermal effects, diffusion coefficients, boundary conditions, physical
inputs, integral-valued Float types, and the full typed GUI solver snapshot (including backend,
initial guess, and postprocessing settings). The checked-in source file must remain unchanged.
The same saved task is solved through direct `solve_from_file` and the GUI calculation path. With
`return_to_dimension: true`, the published coordinate must end at physical `L`, temperature values
must be plausible kelvin values, and the GUI plot must match the direct physical solution
componentwise within `1e-8 * max(|reference|, 1)` (mesh `1e-12`). The test then runs with
`return_to_dimension: false` and checks that the GUI keeps the dimensionless endpoint `z=1` and labels
the temperature as `Teta`, without applying a second conversion.

**Observed result:** Debug GUI story passes for both dimensional and dimensionless publication. The
actual GUI widget edit survives Save As and reload, while the source file remains unchanged. The
physical assertions prevent two routes from agreeing merely because both applied the same extra
scale transform. This also catches a prior serialization defect where integral-valued `Float`s were
written without a decimal point and reparsed as `Integer`, causing otherwise valid solver inputs to
fail validation. File I/O is not performance-measured.

The focused scale contracts are also covered by:

```powershell
cargo test --lib --no-default-features test_set_solver_bc_scales_physical_heat_and_species_fluxes -- --nocapture --test-threads=1
cargo test --lib --no-default-features test_scaling_config_ -- --nocapture --test-threads=1
```

They verify conversion of physical outlet flux inputs to solver coordinates and rejection of
non-finite or non-positive `dT`, `L`, and `T_scale`.

## Story: BVP Document Failure and Result Provenance

**Name:** Rejected file replacements preserve the active result; rejected reruns clear stale output.

**Commands:**

```powershell
cargo test --lib --no-default-features malformed_bvp_replacement_preserves_dirty_active_document_transactionally -- --nocapture --test-threads=1
cargo test --lib --no-default-features rejected_bvp_rerun_clears_previous_plot_result -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does a failed document replacement leave the current task and its
successful result intact, while a new run request invalidates old output before new-input
validation can fail?

**What is checked:** After a successful HMX solve, a malformed replacement is requested while the
document has unsaved edits. The active document, typed GUI snapshot, path, dirty flag, and complete
plot payload remain unchanged. In the second story, an invalid solver backend is submitted after a
successful run; the request fails visibly and the previous plot is cleared.

**Observed result:** Both Debug stories pass. File-load failure is transactional with respect to the
existing result; failed rerun cannot leave old plot data displayed as if it belonged to the rejected
inputs.

## Story: BVP Numerical Failure Recovery

**Name:** A failed HMX BVP solve can be corrected and rerun without stale output.

**Command:**

```powershell
cargo test --lib --no-default-features bvp_numerical_failure_can_be_repaired_and_rerun -- --nocapture --test-threads=1
```

**Question / hypothesis:** When the nonlinear iteration budget is insufficient for the physical
task, does the GUI report a numerical failure, then produce a fresh physical profile after the user
repairs that setting?

**What is checked:** The checked-in HMX task is run with `max_iterations=1`, which must fail without
a plot. The task is repaired to `max_iterations=100` and rerun. The successful result must contain a
finite profile and a physically bounded Kelvin temperature column.

**Observed result:** Debug story passes: the constrained solve fails, the corrected solve completes,
and only the recovered run publishes a new plot.

## Story: BVP Solver-Family Physical Readiness in the GUI

**Name:** The selected solver family reports its measured compact-HMX acceptance status.

**Command:**

```powershell
cargo test --lib --no-default-features bvp_family_switch_exposes_capability_specific_controls_in_compact_window -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can users distinguish physical-route validation status without losing
the ability to inspect or select a solver family?

**What is checked:** In the compact BVP window, the solver panel reports compact-HMX physical
acceptance for Damped, BVP_sci, and Frozen/Naive as the family selector changes. Frozen/Naive remain
selectable; the GUI does not imply physical acceptance where the endpoint contract is not met.
This is presentation of measured status, not a change to solver policy or the acceptance gate.

**Observed result:** Debug GUI story passes. Damped and BVP_sci display accepted status, while
Frozen/Naive display that compact-HMX boundary acceptance remains open. The route selection and
typed solver configuration are unchanged.

## Story: BVP Conditional Controls in a Compact Viewport

**Name:** Every currently supported conditional solver control is reachable in a compact BVP window.

**Command:**

```powershell
cargo test --lib --no-default-features bvp_conditional_solver_controls_are_reachable_in_compact_viewport -- --nocapture --test-threads=1
```

**Question / hypothesis:** In an 820x620 window, can a user scroll each solver-family panel to every
control revealed by the selected execution mode, matrix layout, and callback policy?

**What is checked:** The story checks control bounds after scrolling for Damped Lambdify sequential
and parallel-threshold settings, Damped AOT build controls, and BVP_sci Lambdify dense/sparse/banded
and AOT controls. Conditional sparse factorization, callback threshold, bandwidth, and AOT build
fields are included. The test seeds valid document settings to exercise each rendered branch;
selector interactions and typed setting round-trips are covered by separate stories.

**Observed result:** Debug story passes. Each checked control is brought into the compact viewport
and its widget rectangle intersects the visible screen. This is UI reachability coverage only; it
does not claim solver correctness or extend route capability.

## Story: BVP_sci AOT Document Controls

**Name:** BVP_sci AOT artifact and chunk settings remain document-driven and parser-valid.

**Command:**

```powershell
cargo test --lib --no-default-features bvp_sci_aot_artifact_and_chunk_controls_roundtrip_into_hybrid_api -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the GUI expose the artifact-lifecycle and chunking settings that
the locked RST hybrid parser supports, and does the resulting task document preserve them as typed
solver configuration?

**What is checked:** The AOT BVP_sci panel renders output/handoff paths, crate/module overrides,
and residual/sparse/banded chunk controls. The production hybrid parser reconstructs each supplied
value. The story also verifies that a build profile is omitted for `require_prebuilt`, where RST
rejects that otherwise-inapplicable field.

**Observed result:** Debug story passes against RustedSciThe 0.4.30. All artifact and chunking values
survive the GUI document path; no compiler is invoked. Damped target chunk counts were not part of
this BVP_sci story and are covered separately against the RST 0.4.31 Damped parser below.
Run it with:

```powershell
cargo test --lib --no-default-features malformed_bvp_replacement_preserves_dirty_active_document_transactionally -- --nocapture --test-threads=1
```

## Story: Damped AOT Target Chunk Counts in the GUI

**Name:** Optional Damped AOT chunk-count overrides reach the RST hybrid settings API.

**Command:**

```powershell
cargo test --lib --no-default-features bvp_damped_aot_chunk_targets_roundtrip_into_hybrid_api -- --nocapture --test-threads=1
cargo test --lib --no-default-features bvp_aot_chunk_target_counts_must_be_positive -- --nocapture --test-threads=1
cargo test --lib --no-default-features bvp_conditional_solver_controls_are_reachable_in_compact_viewport -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can users optionally override RST's Damped AOT residual and Jacobian
target chunk counts, while leaving defaults untouched when overrides are disabled and avoiding
irrelevant controls for Lambdify?

**What is checked:** Both controls render on the Damped AOT path, their values are normalized to the
integer representation expected by the production hybrid parser, and the parsed generated-backend
configuration contains the requested `ByTargetChunkCount` policies. Switching to Lambdify removes
the AOT-only values and widgets. Zero is rejected before solver handoff. Hover/help content and
reachability in an 820x620 viewport are covered as well.

**Observed result:** Debug stories pass on RustedSciThe 0.4.31. The GUI forwards residual target
count 5 and Jacobian target count 3 as typed RST policies; zero is rejected, and Lambdify emits no
chunk overrides.

**Conclusion:** The document-driven GUI now exposes the RST 0.4.31 Damped AOT chunk controls without
changing RST defaults or adding these AOT-only settings to the Lambdify route. This verifies settings
construction, not runtime performance or generated-code execution.

## Story: BVP-01 Analytical Transport Reference

**Name:** The production BVP setup and supported solver routes reproduce an independent source-free transport solution.

**Command:**

```powershell
cargo test --lib --no-default-features bvp01_ -- --nocapture --test-threads=1
```

**Question / hypothesis:** Do the KiThe BVP equations, physical endpoint conversions, and dimensional
publication match a closed-form solution for coupled temperature and two-species transport when the
reaction rate and heat source are exactly zero?

**What is checked:** A real `X0 => X1` reaction record with `A=n=E=Q=0` passes kinetic preparation
and symbolic equation assembly. The pseudo-species share a `C2H4O2` atom inventory and the test
checks prepared molar masses, mean molar mass, reference density, `K_i`, Peclet numbers, state order,
source-free RHS values, and converted boundary conditions. An independent stable `F(Pe,z)` reference
is checked at `Pe=0`, `1`, and `30`. The test solves an exact-zero-flux constant-profile case,
moderate equal-transport and unequal-transport cases, and a `Pe=30` outlet-layer case. Damped and
Frozen are checked against the same moderate case; BVP_sci covers all nontrivial transport cases.
Published `x`, `T`, `q`, `J`, and mass fractions are compared componentwise to the physical oracle.
The unequal-diffusivity case explicitly expects `X0+X1 > 1` at the outlet; unit-sum normalization is
not a valid invariant for these independent Fick equations.

**Acceptance distinction:** A returned solver profile is not called analytically accurate solely
because the route reports success. Damped and Frozen use 256 fixed intervals for the `Pe=1` check;
the fixed-grid `Pe=30` Damped and Frozen studies must show decreasing analytic profile error at 32,
64, 128, 256, and 512 intervals. These BVP-01 Damped analytic runs explicitly use the `ExprLegacy`
symbolic frontend; AtomView's adaptive-grid panic is tracked separately in the BVP TODO. A separate
Damped BVP_api run starts at 32 intervals with
`SolverParams.adaptive = Some(AdaptiveGridConfig { version: 1, max_refinements: 2,
grid_method })`; the strategy matrix covers `DoublePoints`, `Pearson(1e-3, 1.4)`,
`GrcarSmooke(1e-3, 1e-2, 1.4)`, and `TwoPoint(0.02, 0.02, 2.5)`. Each route must report actual
refinements, grow the published mesh, and not worsen the physical analytic profile error versus the
same 32-interval seed. `DoublePoints` additionally must strictly improve that error. BVP_sci uses
adaptive collocation and is checked over the whole profile, including the outlet layer.
Production `Pe=0` is not claimed: task validation requires positive mass flux and positive transport
coefficients. The exact `F(0,z)` limit is tested in the independent reference, and a separate facade
case checks exact zero boundary fluxes at positive Peclet numbers.

**Observed result:** Debug correctness stories pass. The zero-rate reaction and formula-based mass
preparation succeed. At `Pe=1`, Damped and Frozen have the same measured maximum physical profile
error of about `7.17e-1` on 256 intervals; BVP_sci returns the expected outlet values for the equal-diffusivity case
(`T=361.362055881 K`, `X0=0.726424112`, `X1=0.273575888`), unequal-diffusivity case
(`X1=0.321306132`, total fraction `1.047730244`), and `Pe=30` case (`T=323.15 K`,
`X0=0.75`, `X1=0.25`; adaptive mesh 68 nodes). Both Damped and Frozen `Pe=30` profile error
decreases monotonically from about `31.13`, `17.99`, `9.79`, `5.14`, to `2.63` as the grid is refined
from 32 through 64, 128, 256, and 512 intervals; this records mesh convergence but also shows that
even 512 intervals leave a visible maximum profile error. The expanded adaptive BVP_api strategy
matrix passes with ExprLegacy: `DoublePoints` grows 33 nodes to 129 and reduces the max profile error
from `31.13` to `9.79`; Pearson `(1e-3, 1.4)` and GrcarSmooke `(1e-3, 1e-2, 1.4)` each grow to
786 nodes, while TwoPoint `(0.02, 0.02, 2.5)` grows to 46 nodes in the current Debug rerun. All three report two refinements
and do not worsen the max error, but each remains at approximately `31.13`; mesh growth alone
therefore does not establish improved global profile accuracy for these settings on this case. A
parameter sweep then varied estimator thresholds and the TwoPoint growth factor, starting from 33
nodes and allowing four refinements. Pearson `(0.01, 1.4)` and Grcar-Smooke `(0.01, 0.05, 1.4)`
grew to 120 nodes but retained the same `31.13` maximum error, located in heat flux near
`z=0.96875`. TwoPoint `(0.02, 0.02, 1.4)` reached 72 nodes in four passes and reduced maximum
error to `2.630348`, matching the 512-interval fixed-grid error with far fewer nodes. The test now
asserts this physical-error improvement and mesh budget; other sweep rows are comparative evidence,
not accuracy claims. The
normalized physical first-integral
diagnostics are near floating-point noise for the tested routes (thermal residual at most
`6.20e-11 W/m²`; species residual at most `4.54e-13 kg/(m² s)` in the moderate Damped/Frozen
comparison).

**Conclusion:** The proposed source-free model is representable by the existing production BVP
pipeline without changing physical equations. Exact zero reaction rate is supported. The test also
exposes a practical mesh-accuracy issue for fixed-grid Damped/Frozen at high Peclet: error decreases
with refinement, but the tested 512-interval grid does not yet match adaptive BVP_sci accuracy.
High-Peclet fixed-grid use therefore needs further refinement/acceptance work, and solver convergence
must not be equated with profile accuracy. These are correctness findings, not performance
measurements.

The fixed-mesh CorrectionNewton route is additionally compared with the same
moderate `Pe=1` physical oracle at 128 and 512 intervals. It must reduce the
whole-profile physical error on the finer grid and retain the RST-native
CorrectionNewton telemetry through KiThe's `last_telemetry`.

**Additional command:**

```powershell
cargo test --lib --no-default-features bvp01_correction_newton_matches_analytic_profile_and_retains_telemetry -- --nocapture
```

**Additional observed result (Debug, RST 0.4.32):** CorrectionNewton's maximum
physical profile error decreases from `1.871144e-3` at 128 intervals to
`1.169458e-4` at 512. The facade retains nonzero residual/Jacobian/linear-solve
counters and callback timing (`0.800 ms` in this run). The route is fixed-mesh;
this is a mesh-convergence check, not adaptive-grid evidence.

## Story: BVP-02 Analytical Multispecies Reaction-Diffusion

**Name:** The production multispecies BVP reproduces an independently assembled linear reaction-diffusion oracle.

**Commands:**

```powershell
cargo test --lib --no-default-features bvp02 -- --nocapture --test-threads=1
cargo test --lib --no-default-features bvp01_ -- --nocapture --test-threads=1
```

**Question / hypothesis:** Do the existing KiThe reaction preprocessing, mass-weighted source terms, transport scaling, interleaved state ordering, mixed boundary conditions, and three BVP solver routes agree with a matrix-exponential reference that is constructed independently of the symbolic KiThe RHS?

**What is checked:** BVP-02a uses the fixed 4-species/5-reaction network and independently asserts its source matrix
`W=[[-0.60,0,0,0],[0.25,-1.05,0,0],[0.35,0.60,-1,0],[0,0.45,1,0]]`.
BVP-02b uses the ordered 12-species/20-reaction network. The test builds `W` from the frozen reaction table,
forms the `2N x 2N` operator, solves the mixed-boundary linear system by LU, checks its residual and estimates
`cond2(E_JJ)` from singular values without explicitly inverting the block. It compares local production RHS
and symbolic Jacobian against the independent operator and finite differences. Damped uses `DoublePoints`
from a constant-inlet/zero-flux guess, Frozen uses a dense 1025-node guess, and BVP_sci uses a non-exact
17-node guess for the 12-species case. Physical `C_i`, `J_i`, `T`, and `q` errors are reported separately.
Mass balance and the integral `sum(J_i)-mL*sum(C_i)` are checked. Sum-of-fractions and sum-of-flux checks are
applied only for equal diffusion; BVP-02c intentionally reports those two criteria as not applicable.
Additional checks cover mesh convergence, ExprLegacy/AtomView parity on BVP-02a, unequal species diffusion,
nonnegative published fractions, and the `alpha={0,0.01,0.1,1}` zero-kinetics continuation.

**Observed result:** Debug tests pass. For BVP-02a, `cond2(E_JJ)=1.728`; Damped/DoublePoints at 257 nodes
has max errors `dC=2.067e-4`, `dJ=1.965e-6 kg/(m2 s)`, `dT=7.173e-2 K`, and `dq=7.173e-1 W/m2`.
Frozen at 1025 nodes reduces these to `5.170e-5`, `4.919e-7`, `1.796e-2 K`, and `1.796e-1 W/m2`;
BVP_sci uses 99 nodes and reaches `6.086e-12`, `1.034e-13`, `1.529e-9 K`, and `1.529e-8 W/m2`.
The BVP-02a ExprLegacy/AtomView Damped profiles differ by at most `1.137e-13`.

For BVP-02b, `cond2(E_JJ)=3.318`. Damped (257 nodes) reports
`dC=9.142e-5`, `dJ=7.172e-7`, `dT=7.173e-2 K`, `dq=7.173e-1 W/m2`; Frozen (1025 nodes)
reports `2.287e-5`, `1.796e-7`, `1.796e-2 K`, `1.796e-1 W/m2`; BVP_sci selects 94 nodes and reports
`4.914e-12`, `8.659e-14`, `1.752e-9 K`, `1.752e-8 W/m2`. Its debug solve takes about 105 seconds.
Equal-diffusion mass-fraction and total-flux invariants remain at roundoff. Frozen mesh refinement from
129 through 1025 nodes strictly decreases each field's reference error. BVP-02c passes with species-keyed
unequal diffusion without normalization. BVP-02d passes for all four kinetic scales and has exactly zero
source at `alpha=0`.

**Conclusion:** The nonzero-chemistry source construction and production RHS/Jacobian agree with the
independent model; all three solver families pass the small and large equal-diffusion cases. Damped and
Frozen show the same thermal discretization error already characterized by BVP-01, while their reactive
species errors decrease with denser meshes. BVP_sci is the most accurate tested route but expensive for
the 12-species debug case. Unequal diffusion does not trigger an invalid composition normalization.
These are correctness findings, not performance guarantees. Release performance observations are now available in the campaign dashboard below; they are single-host baselines, not gates.

CorrectionNewton is additionally tested on BVP-02a at 129 and 257 fixed nodes against the same
matrix-exponential oracle. It reduces all four profile errors from `5.155e-7 / 5.865e-9 /
1.842e-4 K / 1.842e-3 W/m2` to `1.299e-7 / 1.478e-9 / 4.642e-5 K / 4.642e-4 W/m2`
and retains nonzero residual/Jacobian/linear-solve counts and callback timings.

The small-network Damped refinement sweep starts from 33 nodes and allows
three passes. Pearson `(0.01, 1.4)` and Grcar-Smooke `(0.01, 0.05, 1.4)` each
produce 166 nodes and improve every physical-field error over the fixed seed:
`dY=3.756e-4`, `dJ=3.668e-6`, `dT=0.1300 K`, and `dq=1.373 W/m2`, versus
`8.013e-4`, `7.587e-6`, `0.2770 K`, and `2.770 W/m2`. TwoPoint
`(0.02, 0.02, 1.4)` uses 102 nodes and also improves all four, but less than
Pearson/Grcar. This is independent-oracle evidence, not merely a strategy
smoke test.

## Story: BVP-03 Analytical Thermochemical Reaction-Diffusion

**Name:** The production BVP reproduces a linear reaction-diffusion system with reaction heat and an independent energy invariant.

**Commands:**

```powershell
cargo test --lib --no-default-features bvp03 -- --nocapture --test-threads=1
cargo test --lib --no-default-features bvp01_ -- --nocapture --test-threads=1
cargo test --lib --no-default-features bvp02 -- --nocapture --test-threads=1
```

**Question / hypothesis:** With `Ea=n=0`, does the real KiThe model correctly couple the linear species source into heat release, preserve the exact energy first integral, and agree with a full matrix-exponential reference across solver routes?

**Independent model:** The BVP-02 ordered reaction networks and inlet compositions are reused unchanged. Species enthalpies in J/mol are `[0,-40,-25,-65,-15,-80,-35,-100,-50,-90,-60,-115]`; the small network uses the first four. The test forms the signed stoichiometric matrix `S`, computes `Q=-gamma*S^T*h`, reaction rates `r=R_C*C`, and mass sources `omega=D_M*S*R_C*C`, without reading KiThe's symbolic RHS to construct the oracle. It assembles the complete block matrix in `(Teta,q,C-vector,J-vector)` order, then verifies conversion to KiThe's interleaved order `(Teta,q,C0,J0,...)`. The mixed boundary system is solved by pivoted LU, checked for residual and boundary satisfaction, and condition-estimated by singular values without explicit inversion. The TW Padé-13 matrix-exponential action is reused from the test-only IVP validation code.

Production BVP properties remain those of BVP-02: `T(0)=298.15 K`, `q_phys(1)=1000 W/m2`, zero species fluxes at the right boundary, `L=0.01 m`, `m=0.01 kg/(m2 s)`, `rho=1 kg/m3`, `D_i=1e-4 m2/s` for equal diffusion, `lambda=0.1 W/(m K)`, `Cp=1000 J/(kg K)`, and `T_scale=100 K`. Explicit zero species fluxes are used; no default near-zero boundary values are relied on.

The synthetic species all use formula `C2H4O2` and the same molar mass, as in BVP-02. Thus mass and elemental conservation are explicitly checked for these isomer-like conversions; this fixture does not establish behavior for reactions between species with different molecular masses.

**Gamma calibration:** For the 4/5 and 12/20 systems, respectively:

| gamma | max chemical `delta T`, 4/5 | max chemical `delta T`, 12/20 |
|---:|---:|---:|
| 0 | 0 K | 0 K |
| 10 | 0.958 K | 0.320 K |
| 100 | 9.580 K | 3.203 K |
| 300 | 28.739 K | 9.608 K |
| 1,000 | 95.797 K | 32.027 K |
| 10,000 | 957.971 K | 320.271 K |

The production acceptance value is `gamma=300`: it gives a resolvable 5–30 K correction on both networks and `cond2(E_JJ)` estimates `5.838e5` (4/5) and `1.235e7` (12/20). At `gamma=1e4`, the 12/20 boundary condition estimate rises to `1.372e10`; this is a calibration stress point, not an accepted operating fixture. The exact analytic response is linear in gamma. The checked reaction-heat vector for the small network at gamma 300 is `[12000,-4500,12000,7500,7500] J/mol`. The ordered 12/20 vector is `[12000,-4500,12000,-15000,19500,-13500,19500,-15000,12000,-9000,16500,-34500,4500,12000,3000,10500,10500,3000,7500,4500] J/mol`. The first 12 reactions form a closed cycle whose net reaction heat is zero by Hess's law.

**What is checked:** All five/all twenty heats are derived from species enthalpies, not independently hand-entered. Reaction atom-count surrogate and mass-source column sums are asserted. Production RHS and every symbolic Jacobian entry are compared to the independent block operator at a nontrivial state; representative columns also use finite differences. Damped uses `DoublePoints` with two refinements; Frozen uses 1025 nodes for the principal route comparison; BVP_sci starts from 33 nodes (small) or 17 (large) and reports its adaptive result. Profiles are compared at the solver's actual physical mesh coordinates, not by matching adaptive node indices. The measured thermal chemical correction is checked separately from total-temperature error. Both internal and dimensional energy invariants use direct state/enthalpy formulas, with summed absolute term scales and nonzero normalization floors. Reaction-heat trapezoidal quadrature is checked against the independent species-flux identity and its mesh error is measured separately.

**Observed solver comparison at gamma 300:**

| Case | Solver | Nodes | max `dC` | max physical `dJ`, kg/(m2 s) | max `dT`, K | max `dq`, W/m2 | energy invariant relative | heat-integral relative |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 4/5 | Damped | 257 | `2.067e-4` | `1.965e-6` | `1.258e-2` | `1.999e-1` | `2.34e-16` | `5.11e-5` |
| 4/5 | Frozen | 1025 | `5.170e-5` | `4.919e-7` | `3.154e-3` | `5.002e-2` | `6.58e-16` | `1.28e-5` |
| 4/5 | BVP_sci | 76 | `1.015e-11` | `1.735e-13` | `1.423e-10` | `8.175e-10` | `2.02e-16` | `1.54e-7` |
| 12/20 | Damped | 257 | `9.142e-5` | `7.172e-7` | `5.048e-2` | `5.756e-1` | `2.06e-16` | `6.87e-4` |
| 12/20 | Frozen | 1025 | `2.287e-5` | `1.796e-7` | `1.264e-2` | `1.440e-1` | `5.20e-16` | `1.72e-4` |
| 12/20 | BVP_sci | 84 | `5.821e-12` | `1.086e-13` | `7.003e-10` | `2.023e-8` | `1.30e-16` | `1.70e-5` |

The thermochemical small-network Damped sweep starts from 33 nodes and allows
three adaptive passes. Pearson `(0.01, 1.4)` and Grcar-Smooke `(0.01, 0.05,
1.4)` each use 170 nodes and improve all four profile errors over the fixed
seed (`dC=3.698e-4`, `dJ=3.601e-6`, `dT=0.02235 K`, `dq=0.3752 W/m2`, versus
`8.013e-4`, `7.587e-6`, `0.04821 K`, `0.7737 W/m2`). TwoPoint
`(0.02, 0.02, 1.4)` uses 134 nodes with a slightly larger remaining error.
The integrated energy invariant remains near machine precision for all these
configurations.

CorrectionNewton also passes the BVP-03a thermochemical oracle at 129 and 257
fixed nodes. Every profile error decreases on the finer mesh
(`dC/dJ/dT/dq=1.299e-7 / 1.478e-9 / 1.818e-5 K / 5.548e-5 W/m2`), while the
relative energy-invariant defect stays below `2e-13`. The facade retains
nonzero native callback and linear-solve telemetry.

The normalized componentwise errors use floors `1e-6` for `C_i` and physical `J_i`; temperature error is normalized by `max(|T_ref-T_in|,1 K)`, and heat-flux error by `max(|q_ref|,1 W/m2)`. Chemical-temperature error is separately normalized by `max(|delta T_chem|,1 K)`. These ratios are diagnostics, while pass/fail profile gates remain field-specific absolute tolerances.

| Case | Solver | relative `dC` | relative `dJ` | relative `dT` | relative `dTchem` | relative `dq` |
|---|---|---:|---:|---:|---:|---:|
| 4/5 | Damped | `1.006e-3` | `4.789e-3` | `1.368e-4` | `4.376e-4` | `2.302e-4` |
| 4/5 | Frozen | `2.517e-4` | `1.205e-3` | `3.430e-5` | `1.097e-4` | `5.758e-5` |
| 4/5 | BVP_sci | `5.063e-11` | `4.359e-10` | `3.783e-12` | `6.499e-12` | `9.204e-13` |
| 12/20 | Damped | `1.716e-3` | `1.940e-2` | `6.932e-4` | `5.254e-3` | `1.037e-3` |
| 12/20 | Frozen | `4.288e-4` | `4.877e-3` | `1.736e-4` | `1.315e-3` | `2.596e-4` |
| 12/20 | BVP_sci | `6.304e-11` | `9.403e-10` | `6.464e-11` | `2.099e-10` | `3.883e-11` |

The largest relative error in the chemical temperature correction is `5.25e-3` for 12/20 Damped and `1.32e-3` for Frozen; BVP_sci is about `2.10e-10`. Max-absolute-error coordinates `(z_C,z_J,z_T,z_q)` are `(1.000,0.480,1.000,0.203)` for large Damped, `(1.000,0.480,1.000,0.204)` for large Frozen, and `(0.672,0.000,0.406,0.000)` for large BVP_sci. The exact energy-invariant residual maximum is also reported with its coordinate. ExprLegacy/AtomView Damped profiles for the 4/5 thermochemical case differ by at most `2.274e-13`.

**Metamorphic and invariant results:** With unequal diffusion `D_i=1e-4*(1+i/10)`, Damped passes with energy-invariant relative defect `1.07e-15`; the maximum departure of `sum(C_i)` from one is `2.805e-2`, confirming that the test did not renormalize composition. For the equal-diffusion 12/20 case, production profiles at `gamma=0` and `gamma=300` differ in species state by only `4.16e-17`. Across analytic references, species profiles are invariant and temperature/heat response is linear in gamma. The analytic energy-invariant defects are approximately `1.5e-13` relative internally and `8.6e-14` in dimensional form. Exact-reference heat quadrature relative error decreases from `7.70e-5` at 32 intervals to `4.81e-6` at 128 and `3.01e-7` at 512.

Damped and Frozen mesh refinements from 129 to 513/1025 nodes reduce all four profile errors monotonically. The BVP_sci 12/20 Debug solve takes about 99 seconds. All nine BVP-03 tests passed in the Release campaign together with the BVP-01/02 regression groups. These remain correctness results, not performance claims.

## Story: BVP-04 Temperature-Dependent Reaction-Diffusion

**Name:** Positive-activation Arrhenius rates in the production BVP match an independent nonautonomous linear reference.

**Commands:**

```powershell
cargo test --lib --no-default-features bvp04 -- --nocapture --test-threads=1
cargo test --release --lib --no-default-features bvp04c_large_nonautonomous_network_matches_frozen_dop853_profile -- --include-ignored --nocapture --test-threads=1
uv run --with numpy --with scipy python src/SimpleReactorsBVP/analytical_reaction_diffusion_tests/generate_bvp04_reference.py
```

**Question / hypothesis:** Does `SimpleReactorTask` evaluate the real Arrhenius law at the scaled physical temperature, and do its production RHS, symbolic Jacobian, and solver routes reproduce an independently integrated temperature-dependent reaction-diffusion model? The fixture sets `Q=0`, isolating temperature-dependent kinetics from reaction-heat coupling.

**Production contract checked:** The audited implementation uses `k(T)=A*T^n*exp(-E/(R*T))`, with `R=8.314` and `E` in J/mol; `T=T0+T_scale*Teta`; first-order molar rate `k*rho*C_reactant/M_reactant`; mass source `omega_i += M_i*nu_i*r`; and species equations `C_i'=J_i/(rho*D_i)`, `J_i'=m*L*J_i/(rho*D_i)-L^2*omega_i/(rho*D_i)`. For each reaction `A` is chosen so that its specified `k_ref` is recovered at 330 K. The imposed temperature field uses physical right heat flux `q_phys,R=1000 W/m2` and `Pe_q=1`:

```text
T(z)=T0 + q_phys,R*L/lambda
      * (exp(Pe_q*(z-1))-exp(-Pe_q))/Pe_q
```

Internally the boundary value is dimensionless, `qhat_R=q_phys,R*L/T_scale`; the dimensional `1000 W/m2` is not inserted directly into the dimensionless state. `T0=298.15 K`, `T_scale=100 K`, `L=0.01 m`, `rho=1 kg/m3`, `m=0.01 kg/(m2 s)`, `D_i=1e-4 m2/s`, `lambda=0.1 W/(m K)`, and `Cp=1000 J/(kg K)`. Thus `T(1)=361.36205588 K`.

**Independent reference:** A separate Python generator integrates the full fundamental matrix for the nonautonomous first-order system using SciPy DOP853. It solves the mixed boundary condition by a linear solve, not by KiThe's BVP routines. Four tolerance/maximum-step combinations establish convergence; the final two profile states differ by at most `9.992e-16` (4/5), `3.608e-16` (12/20), and `6.106e-16` (O3 case). Frozen CSVs contain 2001 uniformly spaced points. `metadata.json` records the SciPy version, initial fluxes, `cond2(Phi_JJ)`, mixed-boundary residual, convergence delta, and mass/element invariant drift. The measured condition estimates are `1.853` (4/5), `4.152` (12/20), and `1.577` (O3); mixed-boundary residuals are at most `1.32e-20`. Reference mass and elemental invariant relative drift is below `4.4e-15`.

**Coverage and results:** The ordinary Debug command completed with 7 passed, 0 failed, 1 ignored. Production local RHS and every symbolic Jacobian entry are compared against independent analytic blocks; representative derivatives also use finite differences. The isothermal `q_R=0` Arrhenius limit is checked against a constant-operator matrix exponential. The small and unequal-molar-mass O3/O2/O/N2 nonautonomous cases are checked across Damped, Frozen, and BVP_sci; Damped ExprLegacy/AtomView full-profile delta is `2.274e-13`.

| Case | Route | Published nodes | max `dC` | max `dJ` (kg/(m2 s)) | max `dT` (K) | max `dq` (W/m2) | mass invariant | element invariant |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 4/5 | BVP_sci | 155 | `2.015e-8` | `7.379e-10` | `1.422e-10` | `1.422e-9` | `2.711e-16` | `2.408e-16` |
| 4/5 | Frozen | 2049 | `7.814e-5` | `1.212e-6` | `8.980e-3` | `8.980e-2` | `2.033e-15` | `1.807e-15` |
| 4/5 | Damped | 1025 | `1.562e-4` | `2.424e-6` | `1.796e-2` | `1.796e-1` | `2.711e-16` | `3.614e-16` |
| O3/O2/O/N2 | Damped | 1025 | `1.065e-4` | `1.680e-6` | `1.796e-2` | `1.796e-1` | `0` | `3.913e-16` |
| O3/O2/O/N2 | Frozen | 2049 | `5.329e-5` | `8.405e-7` | `8.980e-3` | `8.980e-2` | `1.355e-15` | `1.565e-15` |
| O3/O2/O/N2 | BVP_sci | 145 | `1.495e-8` | `5.731e-10` | `3.782e-10` | `3.782e-9` | `0` | `3.912e-16` |

The errors compare published profiles with the independent DOP853 reference; the species-flux column is converted from internal `J` to physical units. Mass and elemental columns measure solution-level first-integral drift, not reference-profile error. Distinct O3/O2/O/N2 molecular masses are checked against the production parser and stoichiometric mass/O/N balances. Frozen and Damped show fixed-grid thermal discretization error; BVP_sci is much closer to the continuous reference on the completed cases.

The 4/5 nonautonomous case also has a same-seed Damped refinement parameter
sweep (65 nodes, three allowed passes) against the frozen DOP853 profile.
Grcar-Smooke `(0.01, 0.05, 1.4)` produced 173 nodes and improved every field
versus the unrefined 65-node result: `dC=9.981e-4`, `dJ=1.618e-5`,
`dT=0.1175 K`, `dq=1.215 W/m2` (seed: `1.223e-3`, `1.895e-5`, `0.1410 K`,
`1.410 W/m2`). Mass and elemental invariants stayed below `3e-16`. Pearson
`(0.01, 1.4)` was close but slightly worse than Grcar; TwoPoint `(0.02, 0.02,
1.4)` used fewer nodes (126) but did not improve the seed error on all fields.
The active assertion selects the measured Grcar configuration and checks both
physical improvement and invariant preservation.

**Open result:** The 12/20 BVP_sci reference gate remains ignored in ordinary Debug and in the 2026-11-08 default Release campaign. A 128-node cap previously produced `MaxNodesExceeded`; a Debug attempt with 512 nodes did not complete within four minutes and was stopped. It must still be explicitly measured in Release before claiming full large-network BVP-04 coverage. The ordinary BVP-04 and broader BVP Release regression filters passed; these are correctness results, not performance claims.

## Story: Adaptive Grid Strategy Parameterization Across Physical References

**Name:** Tuned Damped grid-refinement parameters improve physical accuracy on representative transport and reaction-diffusion BVPs.

**Commands:**

```powershell
cargo test --lib --no-default-features bvp01_damped_pe30_refinement_parameter_sweep_reports_physical_accuracy -- --nocapture --test-threads=1
cargo test --lib --no-default-features bvp02a_damped_refinement_parameter_sweep_measures_physical_accuracy -- --nocapture --test-threads=1
cargo test --lib --no-default-features bvp03a_damped_refinement_parameter_sweep_measures_thermochemical_accuracy -- --nocapture --test-threads=1
cargo test --lib --no-default-features bvp04b_damped_refinement_parameter_sweep_measures_nonautonomous_accuracy -- --nocapture --test-threads=1
cargo test --lib --no-default-features zz1962b_damped_grid_parameter_sweep_measures_branch_accuracy -- --nocapture --test-threads=1
```

**Question / hypothesis:** Do estimator thresholds and mesh-growth parameters materially change the accuracy/cost tradeoff of production Damped refinement, and which settings improve an independent physical profile rather than merely adding nodes?

**Common setup:** All routes use the production KiThe `BVP_api` Damped path with the `ExprLegacy` frontend. Each comparison uses the same physical case and initial mesh for its seed and candidates; errors are computed against an independent analytical, matrix-exponential, or frozen DOP853 reference. These are correctness/mesh-quality stories, not runtime benchmarks.

| Physical case | Seed and refinement budget | Measured strategy and parameters | Published mesh | Physical result |
|---|---|---|---:|---|
| BVP-01, source-free transport, `Pe=30` outlet layer | 33 nodes, up to 4 passes | `TwoPoint(0.02, 0.02, 1.4)` | 72 nodes | Max profile error falls from `31.13085` to `2.630348`, matching the fixed 512-interval error. |
| BVP-01, same case | 33 nodes, up to 4 passes | `Pearson(0.01, 1.4)`; `GrcarSmooke(0.01, 0.05, 1.4)` | 120 nodes each | Max error stays `31.13085`; worst error remains heat flux near `z=0.96875`. Refinement count/node growth alone would have been misleading. |
| BVP-02a, 4 species / 5 reactions | 33 nodes, up to 3 passes | `Pearson(0.01, 1.4)` or `GrcarSmooke(0.01, 0.05, 1.4)` | 166 nodes each | All four errors improve: `dY=3.756e-4`, `dJ=3.668e-6`, `dT=0.1300 K`, `dq=1.373 W/m2`; seed values are `8.013e-4`, `7.587e-6`, `0.2770 K`, `2.770 W/m2`. |
| BVP-02a, same case | 33 nodes, up to 3 passes | `TwoPoint(0.02, 0.02, 1.4)` | 102 nodes | Improves all four fields, but less than Pearson/Grcar. |
| BVP-03a, thermochemical 4/5 network | 33 nodes, up to 3 passes | `Pearson(0.01, 1.4)` or `GrcarSmooke(0.01, 0.05, 1.4)` | 170 nodes each | All fields improve: `dC=3.698e-4`, `dJ=3.601e-6`, `dT=0.02235 K`, `dq=0.3752 W/m2`; energy-invariant error remains near machine precision. |
| BVP-03a, same case | 33 nodes, up to 3 passes | `TwoPoint(0.02, 0.02, 1.4)` | 134 nodes | Improves all fields, but leaves slightly larger error than Pearson/Grcar. |
| BVP-04b, temperature-dependent 4/5 network | 65 nodes, up to 3 passes | `GrcarSmooke(0.01, 0.05, 1.4)` | 173 nodes | All fields improve: `dC=9.981e-4`, `dJ=1.618e-5`, `dT=0.1175 K`, `dq=1.215 W/m2`; mass/element invariants stay below `3e-16`. |
| ZZ-1962-B, low branch at `eta_max=10` | 65 nodes, up to 3 passes | `GrcarSmooke(0.01, 0.05, 1.4)` | 207 nodes | All branch-profile errors improve: `dT=0.2346`, `dY=7.821e-5`, `dJ=2.637e-7`, `drate=8.110e-6`; peak-rate coordinate error is unchanged at zero. |
| ZZ-1962-B, same branch | 65 nodes, up to 3 passes | `Pearson(0.01, 1.4)` | 197 nodes | Improves all four fields, but less than Grcar. |
| ZZ-1962-B, same branch | 65 nodes, up to 3 passes | `TwoPoint(0.02, 0.02, 1.4)` | 117 nodes | Worse than the fixed seed in all four compared fields; fewer nodes do not make it a better result. |

**Conclusion:** There is no universally best refinement strategy. The tuned TwoPoint configuration is the clear choice for the BVP-01 outlet layer; Pearson/Grcar are more effective on the tested reaction-diffusion cases; Grcar is best among these candidates for BVP-04b and the ZZ-1962-B low branch. Keep these settings attached to their validated cases rather than promoting them as global defaults. The full solver/refinement behavior remains bounded by the tested meshes and fixtures; no universal convergence guarantee is implied.

## Story: CorrectionNewton Across Physical Reference Cases

**Name:** Fixed-mesh CorrectionNewton reproduces independent transport, reaction-diffusion, thermochemical, and literature references while retaining its native telemetry.

**Commands:**

```powershell
cargo test --lib --no-default-features bvp01_correction_newton_matches_analytic_profile_and_retains_telemetry -- --nocapture
cargo test --lib --no-default-features bvp02a_correction_newton_matches_matrix_exponential_and_retains_telemetry -- --nocapture
cargo test --lib --no-default-features bvp03a_correction_newton_matches_thermochemical_oracle_and_energy_invariant -- --nocapture
cargo test --lib --no-default-features zz1962b_correction_newton_matches_low_branch_and_retains_telemetry -- --nocapture
```

**Question / hypothesis:** Does the new fixed-mesh RST route work through KiThe's public BVP facade on more than the source-free transport fixture, including independent species/energy references and a nonlinear literature branch?

**Observed results (Debug, RST 0.4.32):** Every gate uses Sparse CorrectionNewton, explicit per-state tolerances and bounds, and Detailed telemetry. BVP-02a at 129→257 nodes reduces `dY/dJ/dT/dq` from `5.155e-7 / 5.865e-9 / 1.842e-4 K / 1.842e-3 W/m2` to `1.299e-7 / 1.478e-9 / 4.642e-5 K / 4.642e-4 W/m2`. BVP-03a reduces `dC/dJ/dT/dq` from `5.155e-7 / 5.865e-9 / 7.215e-5 K / 2.202e-4 W/m2` to `1.299e-7 / 1.478e-9 / 1.818e-5 K / 5.548e-5 W/m2`; energy-invariant relative error is `1.5e-13` at 257 nodes. ZZ-1962-B low branch reduces `dT/dY/dJ/drate` from `8.312e-4 / 2.771e-7 / 1.962e-9 / 2.729e-8` at 129 nodes to `2.100e-4 / 7.000e-8 / 5.012e-10 / 6.893e-9` at 257 nodes, with `dq=1.504e-3 W/m2` and physical invariants passing. Each gate asserts nonzero residual/Jacobian/linear-solve counts and callback timings.

**Conclusion:** CorrectionNewton is verified on four distinct physical/reference families, not only a toy or transport case. This is fixed-mesh correctness evidence; it does not imply adaptive refinement. GUI handoff, partial failed-attempt telemetry, and task-file lifecycle are covered by the later stories below.

## Story: Release BVP Correctness and Performance Campaign

**Name:** The RST 0.4.32 reactor routes preserve the physical reference gates, and the release dashboard reports comparable stage and backend timings across the supported test matrix.

**Commands:**

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorBvpCampaign.ps1 -Suite correctness -BuildJobs 1
powershell -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorBvpCampaign.ps1 -Suite performance -BuildJobs 1
```

**Question / hypothesis:** After integrating CorrectionNewton, tuning adaptive-grid parameters, and updating the RST 0.4.32 routes, do the independent physical gates still pass in Release? Does the benchmark expose meaningful differences between solver families, symbolic frontends, matrix layouts, and Lambdify/AOT execution without treating unlike routes as interchangeable?

**Observed correctness result (Release):** The correctness campaign had 13 successful steps: one shared test-binary build and 12 focused filters, with no failures or empty filters. BVP-01/02/03/04 analytical groups, the ZZ-1962-B physical stories, HMX symbolic and Frozen/Naive matrices, callback-policy dispatch, KiThe-side preparation parity, and both opt-in AOT correctness stories passed. The selected Release gates include the CorrectionNewton physical references and the tuned-grid assertions recorded above. The intentionally ignored large BVP-04 Sci reference remains outside this ordinary campaign.

**Observed performance result (Release):** Criterion completed 93 cases, including scalar-small, eight-state-medium, and twenty-four-state-large routes, all configured symbolic/matrix axes, and matched AOT-TCC cases where supported. The dashboard also contains detailed route-stage, native telemetry, callback-policy, KiThe preparation, and AOT/Lambdify tables. At twenty-four states, CorrectionNewton ExprLegacy/Sparse measured `1.384 ms`, BandedAuto `2.218 ms`, and Dense `3463 ms`; BVP_sci ExprLegacy/Sparse measured `0.335 ms`. Damped ExprLegacy/Sparse measured `88.15 ms` with Lambdify and `133.95 ms` with AOT-TCC. These are single-host observations, not performance gates or universal rankings; notably, Dense is a poor layout for this large CorrectionNewton case, while AOT is not automatically faster for these Damped runs.

**Campaign records:** RST 0.4.32, `rustc 1.98.1`, Windows 10.0.26200, build jobs `1`. Correctness report: `test_reports/reactor-bvp-correctness-release-20261108-203755-471Z/reports/summary.md`. Performance report: `test_reports/reactor-bvp-performance-release-20261108-205631-718Z/reports/summary.md` (93 Criterion rows and detailed diagnostic tables). The correctness campaign duration includes a cold Release LTO test build (`1035.397 s`); subsequent filtered test execution totaled under 40 seconds. The performance campaign took `3922.565 s`, including its cold benchmark/dependency build.

**Conclusion:** The tested production and diagnostic routes pass the selected Release correctness envelope, and the current dashboard provides actionable route/layout comparisons. Continue to treat solver-specific telemetry according to its own semantics; do not sum callback, native solver, and facade timings as if they were disjoint wall-clock partitions.

## Story: ZZ-1962-B Fixed-Inlet Literature Benchmark

**Name:** The KiThe fixed-inlet-concentration BVP agrees with an independent reduced scalar model and preserves its thermochemical invariants.

**Commands:**

```powershell
cargo test --lib --no-default-features zz1962b -- --nocapture --test-threads=1
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_zz1962b/generate_reduced_reference.py
```

**Question / hypothesis:** Does the existing `SimpleReactorTask` represent the requested first-order, heat-releasing `A -> B` problem consistently, and do its BVP_sci, Damped, and Frozen solutions agree with a separately solved scalar boundary-value problem on the same steady branch? This is a correctness story, not a performance measurement.

**Literature scope:** The motivating source is R. M. Zaydel and Ya. B. Zel'dovich, “On possible regimes of stationary combustion,” *Journal of Applied Mechanics and Technical Physics*, 1962, vol. 3, issue 4, pp. 27–32. The bibliographic record is indexed by [Math-Net.Ru](https://www.mathnet.ru/php/contents.phtml?issue=4&jrnid=pmtf&option_lang=rus&series=0&volume=3&wshow=issue&year=1962) and the journal issue listing [Siberian Branch RAS](https://sibran.ru/journals/issue.php?ARTICLE_ID=160090&ID=9122&PAGEN_1=189). The checked case is deliberately the task's **variant B**, not a claim of exact article boundary-condition reproduction.

**Dimensional setup:** `A -> B`; both species have formula `C2H4O2` and molar mass `0.060052 kg/mol`. Values are `T0=300 K`, `rho=1 kg/m3`, `m=0.01 kg/(m2 s)`, `Cp=1000 J/(kg K)`, `lambda=0.1 W/(m K)`, `D_A=D_B=1e-4 m2/s`, `tau=0.5 s`, `A0=2 s-1`, `E=7 R T0=17459.4 J/mol`, `Q_mass=10 Cp T0=3.0e6 J/kg`, `Q_molar=Q_mass*M=180156 J/mol`, and `L=0.05 m`. The Arrhenius exponent is `n=0`.

**Transport-reference adaptation:** The current `SimpleReactorTask::transport_coefficients()` derives the reference density from molar mass, pressure, and a reference temperature, and scales the transport coefficient with that temperature. To realize the requested constant `rho=1`, `D=1e-4`, and `Le=1` under that existing contract, the fixture sets the transport reference to `Tm=298.15 K` and `P=rho*R*298.15/M`; the physical left boundary remains `T(0)=300 K`. This is a test-fixture accommodation to the current API, not a change to the physical boundary temperature or production model. The asserted groups are `alpha=D/((m/rho)^2 tau)=2`, `beta=E/(R T0)=7`, `delta=Q_mass/(Cp T0)=10`, `Le=lambda/(rho D Cp)=1`, and `eta_max=L/((m/rho)tau)=10`.

**Mapping and boundary conditions:** KiThe uses `z=x/L` and physical `T=T0+3000*Teta`; the article coordinate is `eta=x/(u tau)`, with `u=m/rho`, so `eta=10 z` for this case. The explicit KiThe boundary values are `T(0)=T0`, `Y_A(0)=1`, `Y_B(0)=0`, and at `z=1`: `q=0`, `J_A=0`, `J_B=0`. The reduced variable is `theta=(T-T0)/3000=Y_B=1-Y_A`, which yields:

```text
theta' = alpha*theta'' + (1-theta)*exp[-beta/(1+delta*theta)]
theta(0)=0, theta'(eta_max)=0
```

The article's inlet condition is instead the total incoming A flux, `m*Y_A(0)-rho*D*Y_A'(0)=m`; this fixture pins `Y_A(0)=1`. The resulting `epsilon_inlet=(rho*D/m)*abs(Y_A'(0))` is reported as a measure of the neglected inlet diffusive contribution. Therefore the computed profiles characterize the stated variant B and are not presented as exact article curves.

**Independent reference and branch evidence:** `generate_reduced_reference.py` solves the scalar equation with SciPy `solve_bvp`, independently of KiThe's symbolic RHS and BVP solver. It stores 2001-point CSV profiles and metadata including SciPy version, collocation residuals, BC residuals, and tolerance/mesh refinement. At `eta_max=10`, multistart found three distinct solutions (low/middle/high) with respectively `Y_A(out)=0.98815695 / 0.85095026 / 0.14798157`, `T(out)=335.529 / 747.149 / 2856.055 K`, and `epsilon_inlet=2.13703e-3 / 5.88181e-3 / 1.03743e-1`. Fine-refinement profile deltas are at most `5.3e-11`, `3.9e-10`, and `2.3e-10` for those branches. Multistart outcomes are evidence of multiple roots, not a proof that all roots have been found.

The independent high branch at finite domains gives:

| `eta_max` | `Y_A(out)` | `T(out)`, K | `epsilon_inlet` | peak-rate `eta` |
|---:|---:|---:|---:|---:|
| 10 | 0.14798157 | 2856.055 | `1.037427e-1` | 4.555 |
| 15 | 0.03016502 | 3209.505 | `1.041308e-1` | 4.545 |
| 20 | 0.00606029 | 3281.819 | `1.041322e-1` | 4.540 |

At `eta_max=15`, one tested seed failed because SciPy exceeded its adaptive node limit; at `eta_max=20`, three seeds failed (node limit or singular collocation Jacobian). These are recorded as failed initial guesses, not as evidence that the corresponding branches cease to exist. The high branch approaches strong reactant depletion as the finite domain grows, but zero outlet gradient alone is not treated as proof of complete conversion.

**Production-route results:** All three production routes were initialized from the independent low branch at `eta_max=10`. Each matched that branch, passed the first integral, species sum, flux sum, and BVP-03 energy-invariant checks. Maximum profile errors compare against the frozen independent scalar profile at the solver's actual physical mesh coordinates.

| Route | Published nodes | `Y_A(L)` | `T(L)`, K | max `dT`, K | max `dY` | max `dJ` | max `d(rate)` | `epsilon_inlet` | energy-invariant relative drift |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BVP_sci | 219 | 0.98815695 | 335.529149 | `1.908e-6` | `6.36e-10` | `2.133e-12` | `6.608e-11` | `2.137031e-3` | `4.387e-16` |
| Damped / DoublePoints | 1025 | 0.98816747 | 335.497593 | `3.929e-2` | `1.31e-5` | `5.487e-8` | `1.332e-6` | `2.136482e-3` | `1.316e-15` |
| Frozen | 2049 | 0.98816221 | 335.513359 | `1.967e-2` | `6.555e-6` | `2.746e-8` | `6.668e-7` | `2.136756e-3` | `1.974e-15` |

The directly evaluated first-integral defect is at most `4.3e-14 K`; species-sum and flux-sum defects remain below their asserted `2e-5` bounds. Damped uses two `DoublePoints` refinements from 256 points; Frozen uses 2049 points; BVP_sci starts at 101 points with tolerance `1e-9`. These route-specific discretizations explain why their profile errors differ; no performance ranking is implied.

The tuned-grid story additionally finds Grcar-Smooke `(0.01, 0.05, 1.4)` as
the best measured Damped setting on the low branch: from the same 65-node
seed it reaches 207 nodes and improves `dT/dY/dJ/drate` to
`0.2346 / 7.821e-5 / 2.637e-7 / 8.110e-6`; Pearson `(0.01, 1.4)` is slightly
worse, while TwoPoint `(0.02, 0.02, 1.4)` increases all four errors.

CorrectionNewton independently solves the same low branch at 129 and 257
fixed nodes. At 257 nodes, `dT/dY/dJ/drate` is
`2.100e-4 / 7.000e-8 / 5.012e-10 / 6.893e-9`, heat-flux error is
`1.504e-3 W/m2`, and the temperature/species/energy invariants pass. All four
profile errors decrease from the 129-node run, and native telemetry counters
and callback timings are retained.

**Conclusion and remaining boundary:** At the time of this baseline, the production-route gate covered the low branch at `eta_max=10`, while multibranch behavior and finite-domain references were characterized independently in the reduced model. The later **High-Branch Validation** story below extends the full BVP_sci evidence to the high branch and `eta_max=15/20`. The inlet BC still differs from the article and the high-branch inlet correction is about 10.4%, so quantitative coincidence with published curves is not claimed. The BVP-01 through BVP-04 stories remain unchanged.

## Story: ZZ-1962-B High-Branch Validation

**Name:** The full production BVP can reproduce the high-temperature branch at the nominal domain, and each route's branch selection is measured explicitly.

**Commands:**

```powershell
cargo test --release --lib --no-default-features zz1962b_full_bvp_sci_reproduces_high_temperature_branch_eta10 -- --include-ignored --nocapture --test-threads=1
cargo test --release --lib --no-default-features zz1962b_full_bvp_sci_reaches_high_branch_from_smooth_nonreference_seed -- --include-ignored --nocapture --test-threads=1
cargo test --release --lib --no-default-features zz1962b_full_bvp_damped_reports_branch_selection_eta10 -- --include-ignored --nocapture --test-threads=1
cargo test --release --lib --no-default-features zz1962b_full_bvp_frozen_reports_branch_selection_eta10 -- --include-ignored --nocapture --test-threads=1
cargo test --release --lib --no-default-features zz1962b_high_branch_bvp_sci_continues_domain_10_to_15_to_20 -- --include-ignored --nocapture --test-threads=1
```

These are addressable, ignored-by-default high-branch checks because Debug BVP_sci can take many minutes. They are correctness tests, not performance gates. The `eta=10` strict gate uses BVP_sci with 301 initial nodes, `max_nodes=1200`, and tolerance `1e-9`; it leaves all benchmark physical parameters and boundary conditions unchanged.

**Acceptance A — high branch at `eta_max=10`: PASS (Release).** Starting from the frozen independent high-branch profile, production BVP_sci returned 1094 nodes and matched the full high-branch profile. Published outlet values were `YA=0.147981571`, `T=2856.055286 K`; `epsilon_inlet=1.0374271e-1`. The reaction-rate maximum is at approximately `eta=4.555` (`x=0.022775 m`); the measured width above 10% of peak is `0.0447083 m`. The direct high-branch reference comparison reported max errors for `T,Y_A/J_A,normalized rate` within the asserted profile gates.

**Acceptance B — non-reference seed: PASS (Release).** A smooth seed was generated without reading a frozen profile:

```text
theta_guess(eta)=0.85*(1-exp(-0.2*eta))/(1-exp(-0.2*eta_max))
theta'_guess(eta)=0.85*0.2*exp(-0.2*eta)/(1-exp(-0.2*eta_max))
```

The state matrix is formed in the actual order `Teta,q,C_A,J_A,C_B,J_B`; the derivatives are converted to the internal `z=x/L` flux variables. BVP_sci reached the same high branch with 1094 published nodes. Against the high reference, the measured maxima were `dT=2.941e-4 K`, `dY=9.804e-8`, `dJ=1.212e-9`, normalized-rate error `1.099e-7`; peak-location error was `3.333e-3` in `eta`. This satisfies the same strict eta=10 high-branch gate as the reference-profile seed.

**Damped/Frozen branch-selection diagnostic:** Both routes were given the high-branch frozen reference as their initial state, but both converged to the low branch. This is a measured branch-selection/basin result; it is not labelled as a production-equation failure.

| Route | Nodes | Nearest branch | `YA(L)` | `T(L)`, K | max `dT` to low branch | max `dY` | max `dJ` | max `dq`, W/m2 | `epsilon_inlet` | energy-invariant relative drift |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| Damped / DoublePoints | 1025 | low | 0.988167469 | 335.497593 | `3.929e-2 K` | `1.310e-5` | `5.487e-8` | `1.646e-1` | `2.136482e-3` | `2.851e-15` |
| Frozen | 2049 | low | 0.988162214 | 335.513359 | `1.967e-2 K` | `6.555e-6` | `2.746e-8` | `8.239e-2` | `2.136756e-3` | `2.851e-15` |

`dJ` is the maximum physical species-flux error; `dq` is the physical heat-flux error against the nearest independent scalar branch. Both retain positive temperature and bounded mass fractions; thermal first-integral, species-sum, flux-sum, and energy-invariant checks pass. Their inability to select high branch from this seed remains visible as a branch-selection result, not hidden by comparing only solver status.

**Acceptance C — domain continuation: PASS (Release).** Continuation starts from the production eta=10 high-branch solution. Each extension maps the previous published solution through physical `x`, appends a constant outlet tail, then reinitializes the model at the new length so Peclet numbers and flux scales are recomputed. The independent frozen reference is specific to each `eta_max`.

| `eta_max` | `L`, m | nodes | `YA(L)` | `T(L)`, K | max `dT`, K | max `dY` | max `dJ` | max `dq`, W/m2 | `epsilon_inlet` | 10%-rate-zone width, m | max scaled ODE defect | abs defect at max, internal units | row / `z` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 10 | 0.050 | 1094 | 0.147981571 | 2856.055286 | `2.941e-4` | `9.804e-8` | `1.212e-9` | `3.637e-3` | `1.0374271e-1` | `4.470833e-2` | `4.082e-6` | `4.082e-6` | 4 / 1.000000 |
| 15 | 0.075 | 1524 | 0.030165016 | 3209.504953 | `7.218e-4` | `2.406e-7` | `3.031e-9` | `9.092e-3` | `1.0413081e-1` | `6.970833e-2` | `3.144e-6` | `5.526e-6` | 1 / 0.098333 |
| 20 | 0.100 | 1818 | 0.006060286 | 3281.819143 | `1.299e-3` | `4.331e-7` | `5.435e-9` | `1.630e-2` | `1.0413222e-1` | `6.458333e-2` | `4.186e-6` | `9.809e-6` | 1 / 0.074583 |

Domain-extension gates are `dT<2e-3 K`, `dY<1e-6`, `dJ<2e-8`, `dq<0.05 W/m2`, and normalized-rate error `<2e-5`; the stricter nominal eta=10 branch-selection gate remains unchanged. Peak physical rate is about `5.18632 mol/(m3 s)` near `x=0.0228 m` for all three domains. The 10%-of-peak width is a mesh-based diagnostic, not a universal physical flame thickness.

The BVP_sci high-branch checks also assert physical bounds, all endpoint conditions, `T=T0+3000Y_B`, `Y_A+Y_B=1`, `J_A+J_B=0`, the independent BVP-03 energy invariant, the full physical heat-flux profile, reaction-rate profile/peak, and 10%-of-peak zone width. They independently reconstruct the production ODE residual from the published physical profile using a nonuniform three-point finite-difference derivative, transform it back to internal coordinates, and compare against `eq_system`; the maximum rowwise defect is scaled by `1+|RHS|` and required below `1e-4`. The measured maximum over eta=10/15/20 was `4.186e-6`, providing a 24x margin. This is a mesh-based diagnostic, not the native collocation residual. Full experimental stability of the branch is outside this steady-BVP benchmark.
# Story: Larger Cantera Subset Multi-Route Capability Probe

**Name:** A connected multispecies Cantera subset is validated through `SubsData`, then the Damped fixed-grid outcome is recorded without treating `Ok(())` as convergence.

**Commands:**

```powershell
cargo test --lib --no-default-features larger_cantera_subset -- --nocapture
cargo test --lib --no-default-features larger_cantera_subset_bvp_sci_solves_the_same_physical_case -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features larger_cantera_subset_frozen_dense_65_node_diagnostic -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features larger_cantera_subset_frozen_banded_257_publishes_physical_profile -- --include-ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Can Damped produce a converged baseline for a larger, thermochemically complete gas-phase reaction set, so fixed and refined grids can later be compared against physical and independent solution-quality checks?

**Fixture:** Five irreversible elementary reactions (Cantera/GRI30 records 296, 297, 299, 300, 301), twelve species including N2, 1200 K inlet, 0.01 m reactor, 0.03 kg/m3 inlet density, 0.2 kg/(m2 s) mass flux, and 2e-5 m2/s species diffusivity. The test resolves Cp, enthalpy, and molar mass through `SubsData`/`calculate_molar_mass`; reaction element balances are asserted. Thermochemistry comes from `NASA_gas` and `nuig_thermo`. Unsupported reversible, falloff, and third-body reactions are excluded because `FastElemReact` does not model those forms.

**Observed result (Debug, RST 0.4.31; historical):** Preparation and balance story passes. The sparse Lambdify / ExprLegacy Damped route with trapezoid discretization, 33 requested nodes, `abs_tolerance=1e-5`, `max_jac=1`, and up to 12 damping trials returned telemetry status `NoResult` after 3 iterations, 2 damping trials, and 5 linear solves. This status is superseded by the current-version measurements below.

The same physical input converges through BVP_sci from the inlet-profile seed. It publishes 78 nodes, temperatures from 1200.000 to 1200.906 K, maximum species change `2.473e-3`, and maximum mass-fraction-sum error `1.11e-16`; the story also asserts finite bounded species fractions. The Debug run took about 76 s, so the test is ignored by default and run explicitly. This supplies a solver-route comparison and a useful numerical reference, but it is not an independent analytic oracle and cannot be compared to a Damped profile until Damped converges.

On RST 0.4.31, Frozen was checked separately on two fixed grids/routes. `dense_naive` at 65 nodes returned a finite profile with `T=[1200.000, 1201.480] K`, mass-fraction closure error `2.22e-16`, and max species change `2.475e-3`, but minimum mass fraction was `-3.385e-4`; it failed physical non-negativity. The 257-node `banded_frozen` attempt failed in the generated residual callback with `NonFiniteCallbackValue` (row 2289, value `inf`). These are historical measurements and are superseded for the 257-node banded route by the current result below.

| Route | Grid / configuration | Numerical outcome | Physical gate |
|---|---|---|---|
| BVP_sci | inlet seed, 33 nodes, adaptive collocation, max 256 | converged; 78 published nodes | pass: bounded finite fractions, nonzero conversion, closure `1.11e-16` |
| Frozen | `dense_naive`, 65 nodes | finite profile returned | fail: minimum mass fraction `-3.385e-4` |
| Frozen | `banded_frozen`, 257 nodes (RST 0.4.31) | typed `CallbackExecutionFailed` wrapping `NonFiniteCallbackValue(inf)` | fail (historical) |
| Damped | Sparse Lambdify / ExprLegacy, trapezoid, fixed 33 nodes | `NoResult` after 3 iterations | no profile to assess |
| Damped | same initial mesh, `DoublePoints`, one refinement | RST panics while reconstructing failed solve | no profile to assess |

The ignored one-refinement reproducer historically panicked in RST `BVP_utils.rs:824` while constructing the failed solve result. RST 0.4.32 fixes that failure-path panic; current positive refinement assertions are recorded below.

**Historical conclusion (RST 0.4.31):** The earlier comparison could not isolate refinement strategy because the fixed 33-node Damped solve did not converge. RST 0.4.32 changes that result; the updated strategy comparison follows.

## Story: Large Cantera Subset Grid and Solver Routes (RST 0.4.32)

**Name:** On the 12-species Cantera subset, Damped grid refinement improves both physical admissibility and agreement with the converged BVP_sci profile; CorrectionNewton requires a denser initial mesh than 33 nodes.

**Commands (Debug; intentionally excluded from the ordinary fast suite):**

```powershell
cargo test --lib --no-default-features larger_cantera_subset -- --nocapture
cargo test --lib --no-default-features larger_cantera_subset_damped_grid_strategy_matrix_reports_physical_quality -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features larger_cantera_subset_correction_newton_mesh_diagnostic -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features larger_cantera_subset_bvp_sci_solves_the_same_physical_case -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features larger_cantera_subset_frozen_dense_65_node_diagnostic -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features larger_cantera_subset_frozen_banded_257_publishes_physical_profile -- --include-ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Do the RST 0.4.32 Damped fixes make the previously failing large mechanism solvable, and does adaptive refinement improve the published physical profile rather than merely increasing node count? How does fixed-mesh CorrectionNewton behave under a controlled initial-mesh sweep?

**Shared physical fixture and reference:** The same five irreversible elementary Cantera/GRI30 reactions (records 296, 297, 299, 300, 301), twelve species, 1200 K inlet, 0.01 m reactor and thermochemistry resolved via `SubsData` are used across routes. BVP_sci from the same 33-node inlet seed is the numerical comparison profile, not an independent analytic oracle; it publishes 78 nodes, `T=1200.000..1200.906 K`, maximum species change `2.473e-3`, and closure error `1.11e-16`.

**Damped conditions:** ExprLegacy/Lambdify, Sparse, trapezoid, `abs_tolerance=1e-5`, `max_jac=1`, at most 12 damping trials, 33-node inlet seed, and at most two refinement passes. Tuned strategies are DoublePoints; Pearson `(1e-2, 1.4)`; Grcar-Smooke `(1e-2, 5e-2, 1.4)`; and TwoPoint `(2e-2, 2e-2, 1.4)`. The test now asserts convergence, finite plausible temperatures, closure `<1e-3`, `min(Y)>=-1e-12`, and strict improvement in all four interpolated profile errors (`T`, species fractions, species fluxes, heat flux) over the fixed-grid Damped baseline.

| Damped route | Published nodes | `min(Y)` | max closure error | max `dT`, K | max `dY` | max `dJ` | max `dq`, W/m2 | refinements / iterations / linear solves |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Fixed grid baseline | 34 | `-1.296e-4` | `2.22e-16` | `2.433e-1` | `4.268e-4` | `1.880e-5` | `1.231e2` | 0 / 4 / 12 |
| DoublePoints | 133 | `-1.061e-17` | `5.940e-14` | `1.812e-2` | `3.195e-5` | `1.562e-6` | `9.655` | 2 / 6 / 16 |
| Pearson `(0.01,1.4)` | 294 | `0` | `2.949e-13` | `1.280e-3` | `1.431e-6` | `2.575e-7` | `4.631e-1` | 2 / 6 / 16 |
| Grcar-Smooke `(0.01,0.05,1.4)` | 296 | `-4.105e-18` | `2.949e-13` | `1.280e-3` | `1.431e-6` | `2.575e-7` | `4.631e-1` | 2 / 6 / 16 |
| TwoPoint `(0.02,0.02,1.4)` | 70 | `0` | `3.253e-14` | `1.812e-2` | `3.195e-5` | `1.562e-6` | `9.655` | 2 / 6 / 16 |

**CorrectionNewton mesh sweep:** Sparse CorrectionNewton, ExprLegacy/Lambdify, explicit per-state tolerances and nonnegative species bounds, Detailed telemetry. The 33-node seed fails after one accepted step and three bound-limited steps (4 residual evaluations, 2 Jacobian evaluations, 4 linear solves); the error retains this partial telemetry. At 65 and 129 requested nodes it converges and publishes physically admissible profiles: 66/130 published nodes, respectively, `T` range `1200.000..1200.902/1200.903 K`, `min(Y)=0`, closure `2.22e-16/4.44e-16`; each requires 6 residual evaluations, 1 Jacobian evaluation and 6 linear solves. The diagnostic asserts at least one successful mesh and physical checks on every published result; mesh-dependent failure remains visible instead of being hidden by an all-routes-must-pass gate.

**Current Frozen result (Debug, RST 0.4.33):** `banded_frozen` at 257 nodes converges in 7 iterations and publishes a physical profile: `T=[1200.000,1200.911] K`, `min(Y)=0`, `max(Y)=0.7656`, species closure error `2.22e-16`, and maximum species change `2.472e-3`. This replaces the earlier 0.4.31 callback-failure characterization. The dense-naive 65-node result remains a separate historical physical-bounds failure.

**Conclusion:** This is direct evidence that on this specific large mechanism refinement improves Damped solution quality, and the current banded Frozen route is also physically admissible at 257 nodes. Pearson/Grcar give substantially tighter Damped agreement with BVP_sci than DoublePoints/TwoPoint at the tested two-pass budget. CorrectionNewton is viable on this fixture from 65/129-node seeds, but not the 33-node seed. All current measurements are Debug correctness diagnostics; no performance claim is made.

## RST-Native 0.4.32 Regression

**Name:** RST 0.4.32 preserves the large-system BVP route fixes using only RST `Expr` and `BVP_api`, without constructing a KiThe reactor or thermochemical data.

**Commands:**

```powershell
cargo test --lib --no-default-features rst_native_ -- --nocapture --test-threads=1
cargo test --lib --no-default-features rst_native_sci_route_accepts_the_frozen_equation_system -- --ignored --nocapture --test-threads=1
cargo test --lib --no-default-features rst_damped_adaptive_grid_strategy_diagnostic_matrix -- --nocapture --test-threads=1
cargo test --lib --no-default-features correction_newton_group_bounds_and_tolerances_expand_to_reactor_states -- --nocapture
```

**Question / hypothesis:** Is the observed behavior produced by KiThe's reactor/formulation layer, or can the same explicit first-order system trigger it through the public RST symbolic BVP interface alone?

**Fixture and setup:** [rst_native_large_bvp_reproducer.rs](rst_native_large_bvp_reproducer.rs) fixes the 26 RHS expressions, abstract state ordering `u0..u25`, mesh, boundary conditions, and initial profiles. Equations are parsed with RST `Expr::parse_expression`; all solves use `BvpProblem` and `BvpSolver` directly. The only preprocessing is a token-for-token rename from exported coordinate labels to `u0..u25`; no KiThe model, reaction parser, thermochemistry, scaling facade, or postprocessing is called. Sci uses a 33-node seed, ExprLegacy/Lambdify/Dense, tolerance `1e-5`, and a 256-node cap. Frozen uses `dense_naive`, 65 nodes, and tolerance `1e-5`. Damped uses ExprLegacy/Lambdify/Sparse, trapezoid discretization, 33 intervals, tolerance `1e-5`, `max_jac=1`, and at most 12 damping trials; the adaptive variant enables one `DoublePoints` refinement.

**Observed result (Debug, RST 0.4.32):** Fixed-grid Damped converges with unchanged fixture/configuration: 4 iterations, 8 damping trials, 12 linear solves. One-pass `DoublePoints` also converges (5 iterations, one refinement) without panic. Damped tests verify every published state/node lies within its configured inclusive bounds. Frozen `dense_naive` still returns a finite 65-node profile with minimum even-index state `-3.385104e-4`; its group-sum invariant remains `2.22e-16`, so numerical convergence does not establish admissibility. The separate six-state adaptive matrix passes all eight ExprLegacy/AtomView x Sequential/Parallel x Pearson/Grcar combinations, each with two refinements and a 786-node result. The Sci route test remains ignored in Debug due its long runtime; refresh its current-version profile when running the full overnight matrix.

**Conclusion:** RST 0.4.32 fixes the tested Damped trapezoid Jacobian/RHS-coordinate defect, AtomView adaptive coefficient-merge failure, and `DoublePoints` failed-result panic. Retain these as positive assertions without loosening tolerances. Frozen's negative profile is a distinct unresolved numerical/acceptance result. The fixture isolates RST routes from KiThe model construction.

## CorrectionNewton Settings Handoff

**Name:** The hybrid reactor parser recognizes RST `correction_newton` settings and expands grouped tolerances/bounds into per-state keys.

**Command:**

```powershell
cargo test --lib --no-default-features correction_newton_group_bounds_and_tolerances_expand_to_reactor_states -- --nocapture
```

**Question / hypothesis:** Does the RST 0.4.32 hybrid settings route retain native CorrectionNewton options while translating KiThe's grouped `C/J/Teta/q` fields into actual unknown names?

**Observed result (Debug, RST 0.4.32):** The parser-level assertion passes. Group fields expand to each substance state plus `Teta` and `q`; explicit `C0` settings override group defaults. Initial-guess width is checked as one column per interval.

**Conclusion:** This verifies settings translation, grouped bounds/tolerances, telemetry selection, and logging-level forwarding. The end-to-end file and failure lifecycle are covered by the next story; GUI family selection and option handoff are covered below.

## Story: CorrectionNewton Task-File and Failure Lifecycle

**Name:** A physical BVP task file can select CorrectionNewton, publish a bounded profile, and preserve the last success after a failed retry.

**Command:**

```powershell
cargo test --lib --no-default-features correction_newton_physical_task_file_solves_and_reports_telemetry -- --nocapture --test-threads=1
cargo test --lib --no-default-features correction_newton_group_bounds_and_tolerances_expand_to_reactor_states -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the CorrectionNewton path work through the user-facing BVP task-file lifecycle, preserve native telemetry and typed errors, and avoid replacing a previously published result when a later solve fails?

**What is checked:** The test derives a temporary CorrectionNewton variant of the checked-in compact HMX task, changing only the solver settings, and executes it through `solve_from_file`. The physical profile and Detailed residual/Jacobian/linear-solve counters are checked. A second one-iteration attempt must fail with a typed BVP error, retain partial counter telemetry, and leave the prior solution matrix and mesh unchanged. The parser test independently verifies that `loglevel: debug` and `telemetry: counters` are forwarded into native CorrectionNewton options.

**Observed result (Debug and targeted Release, RustedSciThe 0.4.32):** The task-file solve publishes 33 points with `T=600.000..868.984 K`; successful telemetry reports 16 residual evaluations, 2 Jacobian evaluations, and 16 linear solves. The deliberately constrained retry fails after 2 residual evaluations; the preceding physical profile and mesh remain intact. The lifecycle story passes in Release; the parser-forwarding story passes in Debug.

**Conclusion:** The CorrectionNewton family is integrated through task-file parsing, physical solve, success/failure telemetry, typed error propagation, configured logging-level forwarding, and transactional publication. This is correctness evidence, not a release performance measurement.

## Story: CorrectionNewton Selection in the BVP GUI

**Name:** The structured BVP GUI selects fixed-mesh CorrectionNewton and preserves its settings through the RST hybrid parser.

**Command:**

```powershell
cargo test --lib --no-default-features correction_newton_family_roundtrips_solver_options_and_aot_request -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can the user select CorrectionNewton in the Solver tab, configure its frontend/execution/layout and solver limits, and produce settings accepted by the RST hybrid parser, including a typed AOT toolchain request?

**What is checked:** The GUI selects `correction_newton`, renders fixed-mesh controls, serializes Sparse/ExprLegacy/Lambdify settings and per-state tolerances/bounds, then round-trips them through `parse_hybrid_bvp_solver_settings_from_document`. The test switches to AtomView/AOT and verifies parsed settings plus the typed AOT request. Adaptive-grid controls are intentionally absent for this fixed-mesh solver.

**Observed result (Debug, RustedSciThe 0.4.32):** The story passes for both Lambdify defaults and AtomView/AOT. The latter requests the `correction_newton_aot` generated backend.

**Conclusion:** CorrectionNewton is selectable, configured, and parser-validated in the GUI. Its physical task-file solving and failed-rerun publication behavior are covered by the separate lifecycle story.

## Story: Large-Mechanism Physical Solver Benchmark Axis

**Name:** Criterion measures tuned Damped refinement strategies and fixed-mesh CorrectionNewton using the correctness-story Cantera fixture.

**Commands:**

```powershell
cargo test --lib --no-default-features explicit_species_universe_keeps_inert_species_and_is_transactional -- --nocapture --test-threads=1
cargo check --bench reactor_bvp_baseline --no-default-features -j 1
cargo bench --bench reactor_bvp_baseline --no-default-features -j 1 -- reactor_bvp_large_mechanism_physical --noplot
```

**Question / hypothesis:** Does the performance harness use the same physical reaction/thermochemistry setup as the large-mechanism stories, and does it expose the cost of each tuned grid strategy separately from fixed-mesh CorrectionNewton?

**What is checked:** `KinData::set_species_universe` validates and transactionally installs the explicit 12-species universe including inert N2, invalidating composition-derived caches before normal analysis rebuilds them. A shared support fixture then supplies the five Cantera reactions, resolved thermochemistry, reactor properties, and inlet state to correctness tests and Criterion. Benchmark cases are Damped with DoublePoints, Pearson `(1e-2, 1.4)`, Grcar-Smooke `(1e-2, 5e-2, 1.4)`, and TwoPoint `(2e-2, 2e-2, 1.4)`, each with two refinements from 33 intervals; CorrectionNewton uses a fixed 65-interval mesh. Correctness assertions remain in dedicated stories rather than being treated as timing output.

**Observed result (Release, RustedSciThe 0.4.32):** Campaign `reactor-bvp-performance-release-20261109-005551-800Z` completed successfully with TCC enabled. Criterion medians for the large physical fixture were `960.98 ms` (DoublePoints), `2147.3 ms` (Pearson), `2150.1 ms` (Grcar-Smooke), `735.05 ms` (TwoPoint), and `39.869 ms` (fixed-65 CorrectionNewton). These routes have different meshes and error behavior; timing alone is not a ranking of physical quality. The complete campaign report is `test_reports/reactor-bvp-performance-release-20261109-005551-800Z/reports/summary.md`.

**Conclusion:** The physical performance axis is measured in Release and shares its fixture with correctness stories. Use recorded physical errors/mesh evidence alongside these timings; no performance gates are set from this single-host run.

## Story: Release AOT, Callback-Policy, and Physical-Route Measurements

**Name:** The repaired BVP routes complete a Release AOT/Lambdify comparison, expose callback-policy overhead, and report timings for the large physical mechanism.

**Command:**

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorBvpCampaign.ps1 -Suite performance -BuildJobs 1
```

**Question / hypothesis:** With solver, symbolic frontend, matrix layout, and scenario held fixed, does TCC AOT improve end-to-end time over Lambdify? Does callback parallelism pay off on the largest analytic case? How do the tuned physical grid strategies compare in elapsed time?

**What is checked/measured:** The campaign sets `KITHE_RUN_BVP_AOT_BENCH=1`; every AOT result must match its Lambdify profile within the existing `2e-5` gate before timing is reported. Detailed rows separate setup/preparation and native solve; the 24-state/128-interval callback-policy axis records dispatch counters. Criterion separately measures the Cantera-derived physical mechanism with two adaptive refinements and fixed-65 CorrectionNewton. SuperLU is included only if an RST bridge path is configured and native telemetry proves nonzero SuperLU factorization with zero fallback.

**Observed result (Release, 2026-11-09, RST 0.4.32):** Campaign `reactor-bvp-performance-release-20261109-005551-800Z` passed in `2090.273 s`; TCC was available and all matched AOT profiles passed parity. On the 24-state case, selected end-to-end medians were:

| Route | Lambdify ms | AOT-TCC ms | AOT delta |
| --- | ---: | ---: | ---: |
| CorrectionNewton / ExprLegacy / Sparse | 1.8181 | 1.9346 | +6.4% |
| Damped / ExprLegacy / Sparse | 98.6367 | 153.6052 | +55.7% |
| BVP_sci / ExprLegacy / Sparse | 0.4149 | 0.5776 | +39.2% |
| Frozen / AtomView / Sparse | 37.7872 | 72.2022 | +91.1% |
| Naive / AtomView / Sparse | 31.7666 | 72.0425 | +126.8% |

The totals include the measured AOT preparation/build path; they are not callback-only timings. Most sampled large routes favored Lambdify end-to-end, so AOT remains an alternative execution route, not a default speed recommendation.

For the largest callback-policy case, all profiles matched exactly at the printed precision:

| Solver | Policy | Callback dispatches (parallel/sequential) | Callback ms | Solve ms |
| --- | --- | ---: | ---: | ---: |
| Damped | Sequential | 0 / 3 | 0.4407 | 109.8312 |
| Damped | Auto | 3 / 0 | 0.5845 | 106.9692 |
| Damped | Forced Parallel | 3 / 0 | 0.5830 | 110.8377 |
| BVP_sci | Sequential | 0 / 513 | 0.0795 | 0.1862 |
| BVP_sci | Auto | 0 / 513 | 0.0818 | 0.1832 |
| BVP_sci | Forced Parallel | 513 / 0 | 6.9152 | 7.0627 |

Auto selected sequential callbacks for BVP_sci. Forced Parallel took about `38x` the Sequential solve time on this case; for Damped, all three totals were close and callbacks themselves remained below `0.6 ms`. KiThe setup for the 40-species/39-reaction preparation fixture measured `33.2315 ms` Sequential, `27.0049 ms` Auto, and `29.3614 ms` Parallel.

The physical Criterion medians are recorded above. No SuperLU route was present: neither `SUPERLU_LIBRARY` nor `SUPERLU_EXPERIMENTAL_LIBRARY` was configured, so no SuperLU timing is claimed. Full macro/stage/native tables, route medians, confidence intervals, environment, and raw output are in `test_reports/reactor-bvp-performance-release-20261109-005551-800Z/reports/summary.md` and `technical/bvp-release-benchmark.log`.

**Conclusion:** Keep Auto/Sequential as the practical BVP_sci callback choice for this workload; do not force parallelism based on equation count alone. AOT parity is established, but these end-to-end measurements do not show a general speed advantage. SuperLU remains unmeasured pending a bridge-enabled host. These are baselines, not performance gates.

## Story: State Bounds and Positive Kelvin Temperature

**Name:** The KiThe BVP facade applies state bounds and never publishes a nonpositive absolute temperature.

**Command:**

```powershell
cargo test --lib --no-default-features test_bvp_profile_bounds_reject_interior_violation_with_state_and_node -- --nocapture --test-threads=1
cargo test --lib --no-default-features test_positive_kelvin_temperature_becomes_a_teta_bound -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the reactor's affine temperature map turn the physical requirement `T[K] > 0` into an admissible solver bound, while retaining a final dimensional check for routes without bounds support?

**What is checked:** With `T = T_scale*Teta + dT`, the computed floor lies just above `-dT/T_scale`; the physical Kelvin value at the floor is positive. The floor tightens the Damped/CorrectionNewton `Teta` bounds without weakening a stricter user lower bound, and an incompatible upper bound is rejected. A valid three-node state profile passes; a negative `C0` at interior node 1 yields `ReactorError::BvpStateBoundViolation` with state, node, value, and interval. A separate dimensional-profile assertion rejects nonpositive Kelvin temperature. Production validation occurs before assigning the new solution and mesh. A small floating-point slack avoids rejecting roundoff at inclusive state bounds.

**Observed result (Debug and targeted Release):** Both regressions pass: user bounds are merged with the physical temperature floor, incompatible bounds fail clearly, and invalid interior composition/temperature values are rejected with typed context.

**Conclusion:** Damped and CorrectionNewton cannot step below the positive-Kelvin floor; every configured reactor route also has a publication-time Kelvin check. Generic `BVPSolver` instances not attached to a `SimpleReactorTask` have no physical temperature transform and retain their generic state semantics.

## Story: ZZ-1962-B Smooth-Seed Branch Selection

**Name:** Production Damped and Frozen routes converge to a physically admissible branch from the same smooth, non-reference initial profile.

**Commands:**

```powershell
cargo test --lib --no-default-features zz1962b_full_bvp_damped_reports_selection_from_smooth_seed_eta10 -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features zz1962b_full_bvp_frozen_reports_selection_from_smooth_seed_eta10 -- --include-ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Can the production Damped or Frozen route reach the ZZ-1962-B high-temperature branch from a smooth seed that is not sampled from the high-branch reference, or does it return to the already observed low branch?

**What is checked:** Both routes start from `smooth_high_branch_guess(eta_max=10, theta_R=0.85, a=0.2)`. Damped starts with 256 nodes and applies its default two DoublePoints refinements; Frozen uses 129 nodes. Each test validates finite values, composition bounds, first integral, species/flux sums, and energy invariant. The report scores all three reference branches and asserts that the nearest is the high branch, so a change in basin selection is an explicit correctness regression.

**Observed result (Debug, current RST):** Both routes select the high branch from the smooth seed. Damped publishes 1025 nodes (`T_out=2855.188128 K`, `Y_A,out=0.148270624`; high-branch score `1.027e1`). Frozen publishes 129 nodes (`T_out=2848.741768 K`, `Y_A,out=0.150419411`; high-branch score `8.340e1`). Their nearest branch is high by a wide margin over middle/low; first-integral drift is below `9.1e-13`, mass-fraction sum error below `8.9e-16`, flux-sum error below `4.4e-19`, and relative energy drift below `9.0e-16`.

**Conclusion:** With the current facade/RST handoff, both production routes reach the high branch from this smooth non-reference seed. The branch-selection stories now assert this measured result while retaining independent physical invariants; this does not establish universal high-branch convergence for arbitrary seeds.

## Story: Density Model Closure

**Name:** The default fixed-composition ideal-gas density changes chemical sources with temperature while preserving the configured constant `rho*D` transport closure.

**Command:**

```powershell
cargo test --lib --no-default-features SimpleReactorsBVP::density_model_tests -- --nocapture --test-threads=1
cargo test --lib --no-default-features analytical_transport_tests -- --nocapture --test-threads=1
cargo test --lib --no-default-features analytical_reaction_diffusion_tests -- --nocapture --test-threads=1
cargo test --lib --no-default-features story_compact_reactor_solves_through_damped_bvp_api_with_boundary_conditions -- --nocapture --test-threads=1
cargo test --lib --no-default-features matrix_compact_hmx_sci_frontend_and_layout_matrix_matches_boundaries -- --nocapture --test-threads=1
cargo test --lib --no-default-features correction_newton_physical_task_file_solves_and_reports_telemetry -- --nocapture --test-threads=1
cargo test --lib --no-default-features checked_in_physical_bvp_documents_load_and_solve -- --nocapture --test-threads=1
cargo test --lib --no-default-features SimpleReactorsBVP:: -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the new default implement `rho(T)=rho_ref*Tm/T` through the symbolic temperature state, preserve the existing constant-density mode when explicitly selected, and leave transport snapshots unchanged?

**What is checked:** A compact HMX task with `Tm=1500 K`, `dT=600 K`, and `T_scale=600 K` exercises a nontrivial affine temperature map. The tests assert the inverse-temperature law at `Tm/2`, `Tm`, and `2*Tm`, compare the symbolic derivative with its analytic value and central finite differences, and check first- and second-order kinetic source ratios `(Tm/T)` and `(Tm/T)^2`. They also assert that `SimpleReactorTask::new()` selects `IdealGasFixedComposition`, that explicit `Constant` remains temperature-independent, and that `D_ro_map`, `Pe_D`, `Pe_q`, `Cp`, `Lambda`, and the species concentration transport equation are identical between density modes. BVP-01, BVP-02/03, BVP-04, and ZZ-1962-B fixtures explicitly select `Constant` to retain their analytical/reference equations.

**Observed result (Debug):** Seven density-model regressions pass, including full compact-HMX BVP setup, both kinetic orders, the symbolic density and Arrhenius-source derivatives, and the unchanged-transport assertion. The full ordinary `SimpleReactorsBVP` lib suite passed 202 tests with 14 intentionally ignored diagnostics. This includes the BVP-01 analytical transport group, BVP-02/03/04 and ZZ-1962-B references, compact-HMX Damped/BVP_sci/CorrectionNewton routes, task parser, and checked-in task-document solves. The Damped solve used the existing positive-Kelvin bounds.

**Conclusion:** Density variation is confined to reaction concentration/source terms. The model keeps constant `K_i=rho D_i`; it does not claim composition-dependent density or a full gas-kinetic diffusion law. The literary MF-1963 reproduction remains a separate follow-up.

## Story: MF-1963 Independent Mathematical Reference

**Name:** An independent SciPy collocation reference resolves both a weak long-induction state and an intense compact-front state of the first-order Merzhanov-Filonenko flow model.

**Commands:**

```powershell
uv run --with numpy --with scipy --with matplotlib python src/SimpleReactorsBVP/literature_mf1963/generate_reference.py
uv run --with numpy --with scipy --with matplotlib python -m unittest discover -s src/SimpleReactorsBVP/literature_mf1963 -p 'test_*.py' -v
```

**Question / hypothesis:** Can the dimensional first-order model be independently solved without KiThe/RST, with dimensionless equations, finite-domain conditions, source integration, and two clearly distinct physical regimes all checked numerically?

**What is checked:** The reference keeps the original Arrhenius law, fixed `cp`, `lambda`, and `rho*D`, `Le=1`, and `rho(T)=rho0*T0/T`. It uses `theta(0)=0`, `p(xi_max)=0`; the outlet temperature is not prescribed. SciPy `solve_bvp` uses the analytic source Jacobian and a stable-manifold seed. Independent checks cover derivative finite differences, dimensional/dimensionless and delta/Arrhenius round trips, finite physical profiles, boundary conditions, adaptive source quadrature versus the integrated ODE identity, domain extension, collocation tolerance, and a low-regime first-order convective approximation. Each case has 17-digit CSV profiles, continuous-peak diagnostics, metadata, and temperature/source SVG plots. The script imports neither KiThe nor RustedSciThe and does not call a production reactor solver.

**Observed result (Python Debug, NumPy 2.4.6 / SciPy 1.17.1):** All 7 tests pass. Both cases converge with `solve_bvp` status 0 at tolerance `2e-7`; the maximum RMS collocation residuals are `1.991e-7` and `1.993e-7`.

| Case | `delta` | Domain `x_max` | Adaptive nodes | `lambda*dT/dx` at inlet | Source peak `x` | `S_max` | Outlet `T` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MF-LOW | `1e-4` | `120.173 m` | 4514 | `0.0150 W/m^2` | `119.963 m` | `38.2031` | `1800 K` |
| MF-HIGH | `0.1` | `0.005427 m` | 855 | `2.672e5 W/m^2` | `4.358e-4 m` | `3.82031e4` | `1800 K` |

The independently quadrature-computed source integral differs from the ODE identity by `1.86e-8` (MF-LOW) and `1.80e-7` (MF-HIGH). Extending the domain tail from `1e-8` to `1e-10` changes the common-domain theta/slope profiles by at most `2.42e-7 / 2.57e-7` for MF-LOW and `5.52e-8 / 1.23e-6` for MF-HIGH. The loose-to-reference mesh-tolerance comparison (`2e-5` vs `2e-7`) differs by at most `4.96e-4 / 5.08e-4` and `9.36e-6 / 5.13e-4`, respectively. The MF-LOW first-order convective estimate places `theta=90` at `xi=12085.38`, about `0.73%` downstream of the full BVP value.

**Conclusion:** MF-00 now has a reproducible, independent mathematical baseline for a weak long-induction case and an intense near-inlet case. These results validate only the stated reduced model and finite-domain fixtures; they do not establish KiThe parity, a critical delta, uniqueness, or a quantitative digitization of the 1963 figures. The dimensional-to-`FastElemReact` mapping remains a separate review gate before any future KiThe comparison.

## Story: MF-1963 KiThe Mapping, Production Solvers, and Grid Strategy Matrix

**Name:** The MF-1963 dimensional model maps into the existing SimpleReactorTask inputs, and the MF-HIGH profile is reproduced by the KiThe BVP_sci Sparse route from independent seeds.

**Commands:**

```powershell
cargo test --lib --no-default-features mf1963_actual_kithe_inputs_recover_transport_and_scaling_contract -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_bvp_sci_matches_mf_high_full_physical_profile -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_bvp_sci_reaches_mf_high_from_independent_smooth_seed -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_bvp_sci_matches_mf_low_full_physical_profile -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_production_ -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can the actual KiThe facade reproduce the independent model after correctly converting molar heat, diffusivity input, temperature coordinates, heat flux, and species fluxes? How do the production Damped and CorrectionNewton routes compare with the independent physical profile on both regimes, and is a BVP_sci failure tied to a poor MF-00 seed?

**What is checked:** The formal species `A` and `B` have the same `C2H4O2` formula and molar mass `0.060052 kg/mol`; `A=>B` is first order by stoichiometry while `FastElemReact.n=0` is the temperature exponent. Heat is converted from `1.5e6 J/kg` to `90078 J/mol`. The test checks KiThe's reference density, the actual `Diffusion -> rhoD` formula, `Le=1`, `Pe_q=xi_max`, prepared reaction rate with temperature-dependent density, and the six-state seed mapping. MF-HIGH is solved via BVP_api using ExprLegacy/Lambdify/Sparse/StandardFaer from both the independent frozen profile and a smooth analytic seed. Production Damped uses ExprLegacy/Lambdify/Sparse, the same reference seed per regime, and a positive-Kelvin Teta lower bound. Its adaptive strategy matrix compares DoublePoints with measured Pearson `(1e-2,1.4)`, Grcar-Smooke `(1e-2,5e-2,1.4)`, and TwoPoint `(2e-2,2e-2,1.4)`, each capped at two refinements. CorrectionNewton uses Sparse, bounds, and fixed meshes. Both production solvers are exercised on both MF regimes. Independent T, composition, flux, and closure gates are applied after dimensional postprocessing; backend convergence status alone is not counted as physical acceptance.

**Observed result (Debug):** Input mapping passed for both cases: `L_LOW=120.173246 m`, `L_HIGH=0.00542708 m`, both recover `rhoD=1e-4 kg/(m s)` and `Le=1`. MF-HIGH passed from both seeds with 520 output nodes; versus MF-00, maximum temperature error was `1.397e-2 K`, maximum component fraction error `9.315e-6`, heat-flux error `4.764 W/m^2` on a `2.672e5 W/m^2` inlet scale, and mass-fraction closure error zero. These updated values use Arrhenius parameters constructed with KiThe's `R_G=8.314`. The test explicitly uses Sparse; the compatibility convenience wrapper hard-codes Dense.

**Production route matrix (Debug):**

| Regime / route | Configuration | Solver outcome | Independent physical acceptance |
| --- | --- | --- | --- |
| MF-HIGH / CorrectionNewton | Sparse, ExprLegacy/Lambdify, 301-node initial profile, `1e-5` tolerances, bounds | Converged; publishes 302 nodes | PASS: max `T` error `0.2604 K`, max mass-fraction error `1.736e-4`, heat-flux error `65.68 W/m2`, closure `1.33e-15` |
| MF-HIGH / Damped | Sparse, ExprLegacy/Lambdify, DoublePoints x2, 301-node reference seed, abs/rel `1e-8`, positive-Kelvin bound | Reports Converged; 5 iterations, 3 Jacobians, 2 refinements; publishes 1205 nodes | FAIL: selects cold branch (`T=300..300.188 K`), max temperature-profile error `1.500e3 K` |
| MF-LOW / CorrectionNewton | Sparse, ExprLegacy/Lambdify, 2049-node reference seed, `1e-5` tolerances, bounds | Nonconvergence; observed 10 iterations, 56 residual evaluations, 5 Jacobians, 46 damping trials (41 rejected) | FAIL: no published solution |
| MF-LOW / Damped | Sparse, ExprLegacy/Lambdify, DoublePoints x2, 2049-node reference seed, abs/rel `1e-8`, positive-Kelvin bound | Reports Converged; 70 iterations, 4 Jacobians, 2 refinements; publishes 8197 nodes | FAIL: max temperature-profile error `1.355e3 K` |

**Damped adaptive strategy sweep (Debug):** Same two-pass cap, reference seed, backend, tolerances, and bounds in every row. The three tuned indicators are borrowed as starting candidates from earlier KiThe BVP fixtures, not assumed to be universally optimal.

| Regime | Strategy | Output nodes / refinements | Max `T` error (K) | Max mass-fraction error | Relative heat-flux error | Closure error |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| MF-LOW | DoublePoints | 8197 / 2 | `1.355e3` | `9.036e-1` | `9.144e-1` | `0` |
| MF-LOW | Pearson `(0.01,1.4)` | 2130 / 2 | `1.356e3` | `9.037e-1` | `9.851e-1` | `0` |
| MF-LOW | Grcar-Smooke `(0.01,0.05,1.4)` | 2129 / 2 | `1.356e3` | `9.037e-1` | `9.851e-1` | `0` |
| MF-LOW | TwoPoint `(0.02,0.02,1.4)` | 2069 / 2 | `1.355e3` | `9.036e-1` | `9.144e-1` | `0` |
| MF-HIGH | DoublePoints | 1205 / 2 | `1.500e3` | `9.999e-1` | `1.000` | `0` |
| MF-HIGH | Pearson `(0.01,1.4)` | 302 / 0 | `1.500e3` | `9.999e-1` | `9.999e-1` | `0` |
| MF-HIGH | Grcar-Smooke `(0.01,0.05,1.4)` | 302 / 0 | `1.500e3` | `9.999e-1` | `9.999e-1` | `0` |
| MF-HIGH | TwoPoint `(0.02,0.02,1.4)` | 302 / 0 | `1.500e3` | `9.999e-1` | `9.999e-1` | `0` |

No tested grid strategy changes branch selection: MF-LOW refines twice for all four methods but remains far from the reference; on MF-HIGH, only DoublePoints refines while Pearson/Grcar/TwoPoint stop at 302 nodes and still return the cold profile. Thus this sweep does not support the hypothesis that the earlier failure was merely caused by choosing DoublePoints. It instead narrows the problem to Damped's convergence/branch selection or its acceptance behavior for this nonlinear thermal front. Both production routes are exercised, but only CorrectionNewton on MF-HIGH passes the physical profile gate. The MF-LOW BVP_sci attempt with 2049 initial nodes and max 30000 nodes also fails with `BVP_sci route failed`; the MF-HIGH BVP_sci route passes from both seeds. These are measured route limitations, not reasons to weaken the physical gate. MF-00 uses `R=8.31446261815324`, whereas KiThe uses `R_G=8.314`; MF-01 constructs exact requested parameters under `R_G`, which evaluated under MF-00's R gives `beta=0.0500027822`, `gamma=1.000055643271e-2`, and `delta` higher by `0.106%`.

**Historical conclusion, superseded by MF-05/MF-07:** At the time of this
initial mapping story, MF-LOW had no passing production route and the tested
Damped strategies selected the cold branch. Later explicit clustered-mesh
tests changed that result: MF-LOW BVP_sci and Damped/DoublePoints now pass the
physical gate. Retain the original diagnostics as history, but use MF-07 for
the current solver-selection conclusion. The supplied MF-01 attachment still
truncates during section 8; missing acceptance criteria are not inferred.

## Story: MF-HIGH Damped fixed-mesh branch and seed-basin diagnosis

**Name:** Damped converges to the cold MF-HIGH root from reference, smooth, and BVP_sci hot-profile seeds.

**Commands:**

```powershell
cargo test --lib --no-default-features mf1963_damped_fixed_mesh_trace_localizes_cold_branch_selection -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_fixed_mesh_seed_basin_matrix -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_hot_seed_jacobian_refresh_matrix -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_mf_low_fixed_mesh_root_diagnostic -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_corrector_newton_seed_layout_diagnostic -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_hot_continuation_01_to_0125_reports_independent_defect -- --nocapture --test-threads=1
```

**Question / hypothesis:** Was the earlier Damped failure to match the MF-HIGH hot profile caused by adaptive refinement, a malformed/transposed seed, or convergence to another solution branch? Does the same fixed-mesh Damped route preserve a hot profile already found by BVP_sci?

**What is checked:** The first test disables adaptive refinement and prints the state-major seed dimensions/order, physical seed/profile ranges, Damped decision events, native telemetry, and published physical profile. A second test runs the same fixed-mesh Sparse/ExprLegacy/Lambdify Damped setup from (1) the frozen reference seed, (2) an independent smooth seed, and (3) a hot profile returned by BVP_sci, converted back into Damped's dimensionless state coordinates. All three use 301 intervals, `abs/rel=1e-8`, the same bounds, and no refinement. The independently derived MF trapezoid defects and `eta=gamma*theta` first integral are evaluated on each published solution. CorrectionNewton is run on the same fixed mesh from the same reference seed as a production-method control.

A third test holds the BVP_sci hot seed and all other inputs fixed while varying Damped's Jacobian refresh interval over `max_jac=1,2,5,50`.

**Observed result (Debug):** The input is passed in the documented Damped shape `6 x 301`, state order `Teta,q,C0,J0,C1,J1`; this equals `state_count x (mesh.len()-1)` for the 302-point mesh. Every seed spans `300..1800 K`. All three Damped runs report `Converged` after 3 iterations, 6 residual requests, 1 Jacobian request, 3 damping trials, 0 rejected trials, and 0 refinements. Despite distinct seed provenance, all publish the same cold profile, approximately `300..300.1876 K`, with maximum hot-reference errors `1.499863e3 K`, `0.9999086` in mass fraction, and `0.9999207` relative heat flux. Independent defects remain small: `max|r_theta| <= 3.59e-15`, `max|r_p| <= 1.90e-10`, and `max|eta-gamma*theta| <= 4.80e-16`. Thus the output is a numerically valid cold root of the measured MF discretization, not merely a nonphysical profile accepted by a weak facade gate.

CorrectionNewton from the same reference seed passes the hot-profile gate on the same 302-node fixed mesh: `dT=0.2604 K`, `dY=1.736e-4`, relative heat-flux error `2.432e-4`. BVP_sci independently passes at 520 nodes with `dT=1.397e-2 K`, `dY=9.315e-6`, and relative heat-flux error about `1.78e-5`.

Jacobian-refresh sweep: all four values again select the cold branch with the same profile errors. `max_jac=1` uses 2 Jacobian requests; `max_jac=2,5,50` each use 1. Refresh frequency therefore does not explain the branch selection on this fixture.

MF-LOW fixed-mesh check (2049 intervals, reference seed, no refinement): both `max_jac=1` and `50` report convergence to the same cold profile (`T=300..454.046 K`, max `dT=1.353896e3 K`, `dY=0.9025971`, relative heat-flux error `0.3016510`). With `max_jac=1`, Damped uses 10 iterations, 21 residual requests, 5 Jacobian requests, and its independent MF defects are approximately `(1.4e-13, 7.8e-9, 7e-17)` for `(r_theta,r_p,eta-gamma*theta)`. At `max_jac=50`, it uses 66 iterations, 325 residual requests, 2 Jacobian requests, and defects approximately `(1.4e-13, 5.4e-8, 7e-17)`. More frequent refresh improves the defect and work count substantially, but does not select the hot branch. The default-refresh MF-LOW result is therefore only characterized as a near-root at this discretization/tolerance, not asserted to be an exact discrete root.

**Interpretation:** This rules out grid refinement and seed shape/transposition as the MF-HIGH cause. The equations admit at least cold and hot solutions; on this setup Damped's Newton/damping path is attracted to the cold basin even when initialized from the BVP_sci hot solution. Its `Converged` label is legitimate for that cold root, but insufficient as application-level acceptance for a task requesting the hot combustion branch. The RST Damped legacy policy checks reduction/convergence of Newton corrections; the independent residual check here separately confirms the cold root. RST's detailed `DampingTrial` event currently reports `value_b=NaN`, so it does not provide the trial residual needed for a complete iteration trace.

**Conclusion:** Do not change MF equations, relax physical tolerances, or blame refinement/Jacobian reuse for MF-HIGH branch selection. Jacobian reuse does affect MF-LOW convergence quality and cost, but its tested refresh settings still select the cold profile. The next investigation should target branch continuation / basin tracking for Damped, then test whether a continuation path can keep it on the hot branch. The result does not establish that Damped cannot solve the hot branch under every configuration.

### Correction: Damped reduced-seed layout and continuation (2026-11-10)

The earlier Damped continuation above used an invalid seed-matrix layout and its apparent warm Damped results are withdrawn. With mixed endpoint boundary conditions, the public `Damped` initial guess is a `state_count x interval_count` matrix: each row contains that state's free-node values, omitting its own fixed boundary node. Flattening a complete profile node-major and reshaping it into that matrix interleaves state values at the endpoints. The mistake was visible when the next continuation seed was rejected before solver entry (`Teta=138.43` outside `[-19.9,120]`); it was not a Damped convergence failure.

The regression `mf1963_damped_corrector_newton_seed_layout_diagnostic` compares the old full-grid-as-reduced input with the corrected per-state reduced matrix. Both converge to the same cold MF-HIGH root (`Tmax=300.188 K`, `dT=1.499863e3 K`, `dY=0.9999086` from the hot oracle); both independent residual and first-integral checks pass. The corrected `mf1963_damped_fixed_mesh_seed_basin_matrix` repeats from frozen-reference, smooth, and BVP_sci hot profiles; all three select this cold root. CorrectionNewton still passes the hot profile at the same 301-interval fixed mesh (`dT=0.2604 K`, `dY=1.736e-4`), and BVP_sci remains an independent hot-branch control.

The historical reaction-parameter sweep gave Damped a correctly packed CorrectionNewton hot seed at `delta=0.1`, then continued from each Damped result while BVP_sci independently continued the hot branch. With fixed 301 intervals, Sparse/ExprLegacy/Lambdify, `max_jac=1`, and no refinement, the recorded run selected the cold root at `delta=0.1`, continued through `delta=4.142`, and returned `NoResult` at `4.143`; BVP_sci retained its hot branch. The long sweep was removed when the known `delta=0.125` accuracy issue was isolated into a short focused test. Its high-delta values above remain historical measurements, not a currently rerunnable gate.

This brackets the observed Damped cold-continuation failure for this configuration in `(4.142,4.143)`. It does **not** establish a mathematical fold or prove that the cold branch disappears; current Damped decision telemetry lacks trial residual/correction values needed to distinguish branch geometry from solver stagnation. The earlier `(4.1,4.2)` warm-seed account is superseded. The complete debug module run on 2026-11-10 passed: 12 passed, 1 ignored.

### Reaction-parameter continuation

The final test starts from the MF-HIGH hot BVP_sci solution at `delta=0.1`, continues that solution through `delta=10`, and at each point passes the continued hot profile to fixed-mesh Damped (`301` intervals, `max_jac=1`, no refinement). The tested parameters are `0.1, 0.125, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 3.5, 4, 4.5, 5, 10`. BVP_sci retains the hot branch at every point (`T_max=1800 K`; 600 nodes at `0.125`, 621 at `0.15`, 637 at `0.2`, 634 at `0.3`, 727 at `0.5`, 687 at `0.75`, 746 at `1`, 803 at `1.5`, 786 at `2`, 770 at `3`, 767 at `3.5`, 849 at `4`, 868 at `4.5`, 863 at `5`, and 922 at `10`).

Damped returns a cold root at every tested parameter through `delta=4.1`: its `T_max` rises from `300.188 K` at `0.1` to `308.484 K` at `3`, `311.364 K` at `3.5`, `316.617 K` at `4`, and `319.036 K` at `4.1`. The independent discrete defects remain small; at `delta=4`, `(max|r_theta|, max|r_p|, max|eta-gamma*theta|)=(3.566e-15, 1.212e-11, 6.939e-17)`. Starting at `delta=4.2`, Damped returns `NoResult` at `4.2`, `4.3`, `4.4`, `4.5`, `5`, and `10`, rather than publishing a profile. At `4.2`, telemetry reports 99 iterations, 818 damping trials, 721 rejections, and 50 Jacobian requests; at `4.3`, 45 iterations, 378 trials, 335 rejections, and 23 Jacobian requests; at `4.4`, 55 iterations, 472 trials, 420 rejections, and 28 Jacobian requests. At `4.5`, the corresponding counts are 63, 553, 493, and 32. These counts vary between runs, so they characterize an observed run rather than exact thresholds. RST's current `DampingTrial.value_b` remains NaN, preventing reconstruction of the rejected trial residual history.

A separate cold-seed Damped continuation starts at the published cold `delta=4.1` result and retries Damped at `delta=4.2`. It also returns `NoResult` (115 iterations, 957 damping trials, 846 rejections, 58 Jacobian requests in the observed run). Thus the `4.2` failure is not explained solely by using a hot BVP_sci profile as the Damped seed. It still does not prove the cold branch ceases to exist: the run establishes a Damped solve failure from both hot- and cold-branch seeds, not a mathematical fold.

To check whether the cold root itself disappears in the same interval, the test also takes Damped's cold `delta=3` profile and continues it with BVP_sci. This cold continuation succeeds at `delta=3.5` (`T_max=311.364 K`) and `4` (`316.617 K`), then BVP_sci fails at `4.5`. That failure does not prove a saddle-node/fold; it may be a limitation of ordinary parameter continuation or the solve configuration. The supported conclusion is narrower: on this fixed mesh and Damped configuration, the hot branch is outside the observed corrector basin through `delta=4.1`; Damped transitions from a published cold root at `4.1` to `NoResult` at `4.2`, while BVP_sci continues the hot branch throughout. This brackets Damped's observed outcome change in `(4.1,4.2)` but does not establish branch loss versus solver stagnation, nor a universal parameter threshold.

**Final correction to the earlier continuation subsection:** its claims that the hot BVP_sci profile was supplied to Damped at every parameter and that Damped stayed warm through `delta=4.1` came from the invalid node-major seed reshape described above. Do not use its warm-branch numbers or `(4.1,4.2)` bracket. The corrected historical sweep values are retained in the paragraph above, but that sweep is no longer an executable test. The current executable continuation test covers only `delta=0.1 -> 0.125`; rerunning the high-delta cold-branch bracket is separate follow-up work.

### Story: MF-02B Independent Hot Root and Cold Discrete Fold

**Name:** Independent trapezoid roots separate Damped branch selection from the cold discrete fold

**Commands:**

```powershell
cargo test --lib --no-default-features mf1963_independent_trapezoid_hot_root_is_probed_through_damped -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_hot_seed_roundtrip_and_first_scalar_newton_audit -- --nocapture --test-threads=1
cargo test --lib --no-default-features literature_mf1963 -- --nocapture --test-threads=1
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/diagnose_cold_branch.py
uv run --with numpy --with scipy python -m unittest test_reference -v
uv run --with matplotlib python src/SimpleReactorsBVP/literature_mf1963/plot_cold_branch.py
```

**Question / hypothesis:** Is Damped being given a malformed hot profile, or does it leave even an independently solved root of the same scalar trapezoid discretization? Separately, does the cold branch itself end near the observed Damped `NoResult` at `delta=4.143`?

**Hot-seed result:** The six-state physical mapping of the independent 301-interval hot trapezoid root has independent defects `(r_theta, r_p, eta-gamma*theta) = (2.815e-14, 2.998e-12, 2.220e-16)` and reduced-seed packing error `0`. Its scalar residual L2 is `7.321e-12`. The Damped facade receives this seed on the same 301-interval mesh with fixed settings and no refinement, but publishes the cold profile: `Tmax=300.187596 K` instead of `1800 K`, maximum profile temperature difference `1.499863e3 K`, and maximum species-fraction difference `0.9999086`. Damped reports three iterations, six residual calls, two Jacobian calls, three damping trials, and zero rejections. The event records show `alpha=1` and `value_b=NaN` for each trial.

The CorrectionNewton hot profile is not itself an exact root of the independent scalar trapezoid equations: its initial scalar residual L2 is `4.330685e-3`; the independently evaluated first Newton correction has norm `4.361342e-2`, and the full-step scalar trial residual is `2.281375e-8`. Directional finite-difference checks of the independent scalar Jacobian are between `3.0e-8` and `6.7e-8` relative error. These checks validate the independent diagnostic equations/Jacobian, not the unavailable production six-state callbacks.

**Cold-branch result, corrected by MF-02B.1:** The original parameter derivative had an indexing defect: it populated the left BC derivative row and omitted the final p-balance row. After correcting it, testing the full vector against centered finite differences, updating tangents through the bordered system, and locally refining the tangent sign-change interval to arc step `1e-5`, fold estimates are `4.142800580301762` (301 intervals), `4.142791345592427` (601), and `4.142789031353778` (1201). These differ from the old coarse estimates by `9.17e-9`, `1.15e-9`, and `6.25e-9`, respectively. The maximum infinity residuals of the extended corrector/nearest branch roots are `6.94e-11`, `2.14e-12`, and `5.03e-12`. Fixed-scaled `sigma_min` estimates at the nearest fold roots are `1.24e-10`, `2.12e-11`, and `4.89e-12`. The hot trapezoid branch remains near `1800 K` through `delta=4.143` on all three grids. This establishes a fold of these discrete cold branches, not of the continuous BVP. SciPy collocation reaches `delta=4.1427` but hits its mesh-node cap at `4.1428`, so the continuous fold remains unresolved.

**Telemetry finding:** Existing KiThe/RST snapshots are insufficient to reconstruct trial decisions: `DampingTrial` has no explicit accepted/rejected result or decision metrics, and the current `value_b` is `NaN`. RST 0.4.32 Legacy acceptance compares undamped correction norms (`s1 < 1 || s1 < s0`), not residual decrease. The exact accepted-step trajectory and iteration at which the profile leaves the hot branch therefore remain unknown. A direct call to the exposed RST `fun/jac` callbacks is not a supported full-state evaluation path: those callbacks operate on internal chunks and fail when invoked with the packed full state. No production solver, equations, tolerances, or acceptance rules were changed.

**Conclusion:** The seed-layout defect is ruled out. The independently solved hot root also fails to remain hot under Damped, so the issue is not explained by the CorrectionNewton-versus-trapezoid seed mismatch alone. The present evidence does not distinguish a Damped iteration/globalization issue from a production callback/discretization difference, because the full production residual and Jacobian cannot be evaluated independently through the current public API. Damped's `NoResult` at `4.143` is consistent with the independently measured discrete cold fold; a continuous physical fold is not established. Full commands, scaling definitions, raw continuation data, plot, limitations, and an upstream telemetry request are in `literature_mf1963/MF_02B_DIAGNOSTIC_REPORT.md`.

**Observed verification:** Debug MF module, 2026-11-10: `14 passed, 1 ignored`. The detailed run is diagnostic, not a performance measurement.

## Story: MF-1963 Damped Guess Forwarding on RST 0.4.34

**Name:** The BVP facade forwards the documented Damped state-major initial guess without changing its basin.

**Commands:**

```powershell
cargo test --lib --no-default-features rst_damped_state_major_guess_preserves_hot_root -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_ -- --nocapture --test-threads=1
```

**Question / hypothesis:** Did the RST 0.4.34 BVP_api guess-forwarding fix restore the public `state_count x interval_count` state-major contract for mixed left/right boundary conditions? Do the previous MF-1963 Damped branch and strategy stories still reproduce their earlier outcomes?

**Facade gate:** The standalone test uses only symbolic `Expr`, public `BVP_api`, six generic states, and a numeric fixture. The 301-interval full-state trapezoid residual is independently `4.0987e-13`. With the documented state-major input, RST reports initial residual `4.0956e-13`, converges in one iteration, and publishes the same hot root (`max(u0)=99.99999999981`). A control using the internal compact permutation instead has initial residual `5.5167e4` and publishes the small root (`max(u0)=0.0125064`). This locks the public facade behavior; it does not assert that arbitrary same-shaped state matrices encode the same profile.

**Replayed MF results:** On RST 0.4.34, the MF-HIGH fixed-mesh Damped check retains the hot reference branch and passes the existing physical profile gate (`max dT=0.2647 K`, `max dY=1.764e-4`, relative flux error `2.482e-4`). The independent exact-hot-root Damped probe also passes with output temperature difference `7.05e-12 K` and fraction difference `5.00e-15`. The corrected seed-layout comparison passes for both seed constructions (`dT~4.28e-3 K`, `dY~2.85e-6`, relative flux error `~5.03e-6`). Jacobian refresh settings `1/2/5/50` retain the same hot profile.

The full `mf1963_` group finished with `11 passed, 3 failed, 1 ignored` in `242.43 s`. The remaining failures are not the old MF-HIGH cold-branch assertion:

| Case | RST 0.4.34 observation | Classification |
| --- | --- | --- |
| MF-HIGH Damped continuation | Hot branch is retained at `delta=0.1` and `0.125`, but at `0.125` the independent `p` trapezoid defect is `1.445e-4`, above the existing `1e-5` gate | Open accuracy issue; tolerance unchanged |
| MF-LOW fixed mesh, 2049 intervals | `MissingResult` after 37 iterations, 282 damping trials, 247 rejections | Open convergence issue |
| MF-LOW adaptive DoublePoints | `MissingResult` after 2 refinements, 300 iterations, 2270 trials, 1983 rejections | Open refinement/convergence issue |

**Decision:** The RST 0.4.34 facade guess-forwarding defect is closed by the new assertion gate. Earlier conclusions that a correctly supplied MF-HIGH hot seed necessarily switches to the cold branch are superseded. Do not infer that the entire MF-1963 family is solved: MF-HIGH continuation accuracy at `delta=0.125` and MF-LOW fixed/adaptive convergence remain open. Keep their numerical assertions strict and investigate them separately.

## MF-04: Adaptive Damped Strategies and Matched SciPy / RST Probe (2026-11-11)

**Name:** Each MF regime/refinement strategy is measured independently, and the matched MF-HIGH two-state problem produces the same profile in standalone SciPy and RST BVP_sci.

**Commands:**

```powershell
cargo test --lib --no-default-features mf1963_damped_adaptive_strategy_matrix_reports_every_case -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_hot_continuation_01_to_0125_reports_independent_defect -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_two_state_rst_sci_ -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_six_state_rst_sci_matches_standalone_scipy_fixture -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_explicit_mesh -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mfhigh_uniform_mesh_runs_through_production_facade -- --nocapture --test-threads=1
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/scipy_parity/solve_scipy.py
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/six_state_parity/prepare_seed.py
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/six_state_parity/solve_scipy.py
```

**Question / hypothesis:** Did the earlier adaptive strategy story fail before reaching MF-HIGH? Is a discrepancy between SciPy and BVP_sci present when both receive the same numerical equations, boundary conditions, mesh, and serialized seed?

**Adaptive matrix setup:** `adaptive_refinement.rs` runs all eight regime/strategy combinations even when a solve returns an error. Damped uses ExprLegacy/Lambdify/Sparse, `abs/rel=1e-8`, max 150 iterations, `max_jac=50`, and at most two refinements. MF-LOW starts with 2049 nodes; MF-HIGH with 301. The classifier requires the existing physical profile limits plus independent trapezoid defects; additional nodes alone are not counted as accuracy.

| Regime | Strategy | Outcome | Final nodes / refinements | Profile errors `(max dT K, max dY, rel dq)` | Independent defects `(theta,p,integral)` |
| --- | --- | --- | ---: | ---: | ---: |
| MF-LOW | DoublePoints | MissingResult | - / 2 | - | - |
| MF-LOW | Pearson `(1e-2,1.4)` | MissingResult | - / 2 | - | - |
| MF-LOW | Grcar-Smooke `(1e-2,5e-2,1.4)` | MissingResult | - / 2 | - | - |
| MF-LOW | TwoPoint `(2e-2,2e-2,1.4)` | MissingResult | - / 2 | - | - |
| MF-HIGH | DoublePoints | Accuracy fail | 1205 / 2 | `(1.654e-2,1.103e-5,1.517e-5)` | `(3.71e-14,5.86e-7,3.66e-14)` |
| MF-HIGH | Pearson `(1e-2,1.4)` | PASS | 419 / 2 | `(2.048e-2,1.366e-5,3.136e-5)` | `(3.06e-14,1.29e-8,4.44e-16)` |
| MF-HIGH | Grcar-Smooke `(1e-2,5e-2,1.4)` | PASS | 419 / 2 | `(2.048e-2,1.366e-5,3.136e-5)` | `(3.06e-14,1.29e-8,4.44e-16)` |
| MF-HIGH | TwoPoint `(2e-2,2e-2,1.4)` | Accuracy fail | 340 / 1 | `(6.945e-2,4.630e-5,6.415e-5)` | `(2.95e-14,1.84e-5,4.44e-16)` |

MF-LOW failure counters from the final rerun differ by strategy: DoublePoints `294 iterations / 2481 residual requests / 16 Jacobians / 2187 trials / 1906 rejected`; Pearson `157 / 1357 / 15 / 1200 / 1058`; Grcar `164 / 1441 / 15 / 1277 / 1128`; TwoPoint `291 / 2535 / 15 / 2244 / 1966`. All four reached two refinements. RST preparation ranged from `0.535 s` to `2.519 s` in this Debug run; these timings are diagnostic only. Pearson and Grcar both meet the current MF-HIGH profile and independent-defect gates on 419 nodes. DoublePoints and TwoPoint publish finite hot profiles but fail the unchanged independent `p`-defect gate. MF-LOW remains unresolved for every tested strategy.

**Matched SciPy/RST fixtures:** The two-state equations are `theta'=p`, `p'=p-source(theta)` with `beta=0.05`, `gamma=0.01`, `theta(0)=0`, and `p(xi_max)=0`. Both routes consume the same 17-digit mesh/seed and frozen SciPy output per case. SciPy uses analytic `fun_jac`, tolerance `2e-7`; RST uses `BVP_api::Sci`, ExprLegacy/Lambdify/Sparse, the same tolerance and boundaries. Python 3 / SciPy 1.18.1 and RST agree componentwise on MF-HIGH to `1.421e-14` in theta and `2.274e-13` in p (both publish 933 nodes). On MF-LOW's clustered 2049-node mesh, both publish 2551 nodes and agree to `3.837e-13` in theta and `1.990e-13` in p. The SciPy MF-LOW solve takes 3 iterations, max RMS collocation residual `1.993e-7`, and sampled max ODE defect `3.779e-7`.

**Six-state parity:** `six_state_parity.rs` solves the full numerical system
`[Teta,q,C0,J0,C1,J1]` directly through RST `BVP_api::Sci`, not through the
KiThe solver handoff. The fixture coordinates are `z=xi/Pe` on `[0,1]`, matching
the production RHS coordinate; state values are lifted from the same frozen
two-state reference. At a nontrivial state, standalone RHS values are checked
against equations assembled by `SimpleReactorTask` for both regimes. SciPy
1.18.1 and RST use identical 17-digit mesh/seed, equations, boundaries,
`2e-7` tolerance, and a 6000-node cap for MF-LOW.

| Regime / mesh | SciPy | RST BVP_sci | Max component delta `[Teta,q,C0,J0,C1,J1]` |
|---|---|---|---|
| MF-HIGH, 301-node reference seed | 823 nodes / 5 iterations | 823 nodes / 5 iterations, 4 refinements | `[2.84e-14,4.26e-14,1.11e-16,2.17e-19,1.11e-16,2.17e-19]` |
| MF-LOW, 2049-node clustered reference seed | 3817 nodes / 9 iterations | 3985 nodes / 25 iterations, 7 refinements | `[4.60e-11,2.86e-8,4.60e-13,2.86e-13,4.60e-13,2.86e-13]` |
| MF-LOW, 2049-node uniform control | `MaxNodesExceeded` at 6000 | `MaxNodesExceeded` at 6000 | no solution |

This confirms SciPy/RST parity on the full six-state model for the tested
clustered MF-LOW seed. The uniform-grid control fails in both implementations;
the parity test does not establish that every seed or mesh distribution will
converge, and it does not test Damped's adaptive strategy.

Reproduce/update the independent six-state fixture with:

```powershell
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/six_state_parity/prepare_seed.py
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/six_state_parity/solve_scipy.py
cargo test --lib --no-default-features mf1963_six_state_rst_sci_matches_standalone_scipy_fixture -- --nocapture --test-threads=1
```

**MF-LOW mesh-distribution control:** The same MF-LOW equations, 2049 initial values interpolated from one independent profile, tolerance `2e-7`, and 6000-node cap were also tested on a uniform mesh. SciPy and RST both hit their mesh-node cap (`MaxNodesExceeded` in RST). On the 2049-node nonuniform mesh clustered like the independent reference, both converge. This supports a mesh-resolution/distribution diagnosis for the tested failure; it does not prove every Damped refinement policy will succeed. KiThe production Damped's four tested adaptive indicators still return `MissingResult` after two refinements.

**Fixed-mesh route capability (historical, superseded by RST 0.4.35):** Older releases rejected this mesh before their numerical loops. RST 0.4.35 now accepts explicit meshes and full node-aligned profiles for Damped and CorrectionNewton. The current production-facade behavior and its qualified MF-LOW physical outcomes are recorded in MF-05 below; retain this paragraph only as the dated reason for the original limitation.

**Focused continuation accuracy:** `continuation_accuracy.rs` keeps the strict independent `p`-defect threshold `1e-5`. On the same 301-node hot branch, Damped defect increases from `6.509993e-9` at `delta=0.1` to `1.445325e-4` at `delta=0.125`, while the branch remains hot (`Tmax=1800 K`). The profile changes by `89.87534 K`. The characterization currently asserts that the known defect remains above the strict gate; it does not count the `delta=0.125` point as an accepted physical solution.

**Current boundary at the time of MF-04:** No SciPy-vs-RST discrepancy was found in either regime when mesh and seed match. MF-LOW's uniform-grid failure is reproduced in both solvers; an equally sized clustered mesh passes in both, identifying mesh distribution as the demonstrated factor. The then-current public Damped route could not accept that custom mesh. RST 0.4.35 explicit-mesh support and current production behavior supersede that limitation; the continuation-accuracy and Frozen questions remain separate.

## MF-05: RST 0.4.35 explicit meshes and production-facade solver outcomes (2026-11-11)

**Name:** The explicit clustered MF-LOW mesh is transferred through the production reactor facade without resampling, and solver convergence is classified separately from physical-reference acceptance.

**Commands:**

```powershell
cargo test --lib --no-default-features rst_solver_selection -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_explicit_mesh -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_clustered_mesh_runs_through_production_facade_for_selected_solvers -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_uniform_mesh_control_reaches_node_cap_through_production_facade -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mfhigh_uniform_mesh_runs_through_production_facade -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mfhigh_damped_grid_strategy_matrix_matches_reference_gate -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mfhigh_frozen_refresh_policy_matrix_matches_reference_gate -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_frozen_refresh_on_correction_runs_through_production_facade -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_bvp_sci_uniform_low_mesh_reports_node_limit_or_passes_physical_gate -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_mf_low_fixed_mesh_classifies_solver_and_physical_outcomes -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_explicit_mesh_rejection_preserves_published_solution -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does RST 0.4.35's node-aligned nonuniform-mesh support work through KiThe's real `SimpleReactorTask` equation builder, boundary mapping, solver handoff, and dimensional postprocessing? For each route, did the backend merely return a converged discrete root, or did that root also satisfy the independent MF-1963 physical profile gate?

**Setup:** `solve_with_bvp_api_on_mesh` accepts solver-coordinate mesh `z` and a full state-major profile of shape `6 x node_count`, including both boundary nodes. It validates finite, strictly increasing mesh coordinates, configured endpoints, profile dimensions/finite values, and profile values at the configured boundary nodes before invoking RST. The normal transactional solve, telemetry capture, physical checks, publication, and dimensional postprocessing remain in use. Damped treats the supplied mesh as its starting adaptive grid; BVP_sci can adapt its collocation mesh; CorrectionNewton remains fixed-mesh.

**MF-LOW clustered production matrix:** The source mesh has 2049 nodes on `z=[0,1]`; all routes receive the same full six-state profile. BVP_sci uses ExprLegacy/Lambdify/Sparse at `2e-7`; CorrectionNewton uses BandedAuto at `2e-7` absolute/per-state relative tolerances; Frozen uses Banded/Lambdify/trapezoid and `at_high_morm=1e-6`; Damped uses Banded/Lambdify/trapezoid, `max_jac=3`, damping factor `0.5`, and `2e-7` absolute/relative tolerances. Damped refinement limits are two except DoublePoints, which uses four.

| Route | Solver result | Published nodes | Max dT | Max dY | Relative heat-flux error | Independent trapezoid defect | Physical reference gate |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| BVP_sci / Sparse | Converged | 3807-3985 | `1.767e-2 K` | `1.178e-5` | `3.433e-5` | Collocation residual gate | PASS; tiny negative roundoff in depleted fraction is bounded by `1e-12` |
| CorrectionNewton / BandedAuto | Converged discrete root | 2049 | `1.254e2 K` | `8.360e-2` | `1.468e-1` | max `2.075e-8` | FAIL; alternate discrete branch and `min(Y_A)=-2.646e-8` |
| Frozen / `at_high_morm=1e-6` | Converged discrete root | 2049 | `1.254e2 K` | `8.360e-2` | `1.468e-1` | max `1.249e-10` | FAIL; same measured alternate branch as CorrectionNewton, `min(Y_A)=-2.646e-8` |
| Damped / TwoPoint `(0.01,0.01,2.0)` | Converged discrete root after refinement | 2207 | `2.223e1 K` | `1.482e-2` | `2.608e-2` | max `1.135e-7` | FAIL; `min(Y_A)=-3.570e-8` |
| Damped / Pearson `(1e-3,1.4)` | Converged discrete root after refinement | 4574 | `1.188e1 K` | `7.919e-3` | `1.394e-2` | max `2.508e-8` | FAIL; `min(Y_A)=-1.329e-10` |
| Damped / Grcar-Smooke `(1e-3,1e-2,1.4)` | Converged discrete root after refinement | 4574 | `1.188e1 K` | `7.919e-3` | `1.394e-2` | max `2.508e-8` | FAIL; effectively same measured branch/profile as Pearson |
| Damped / DoublePoints, 4 refinements | Converged discrete root after refinement | 32769 | `4.130e-1 K` | `2.754e-4` | `4.857e-4` | max `6.377e-7` | **PASS**; nonnegative within strict gate, closure `2.22e-15` |

The discrete state-closure error is at most `2.22e-15` in this matrix. The explicit `at_high_morm` Frozen attempt returns a solution but does not pass the physical-reference gate; its result matches the CorrectionNewton alternate branch, so it is not counted as a LOW success. DoublePoints passes the LOW physical gate, at the cost of a 32769-node final mesh and a substantially longer solve. The test asserts each measured classification instead of weakening shared acceptance thresholds.

**MF-HIGH production matrix:** CorrectionNewton/BandedAuto and Damped with TwoPoint, Pearson, Grcar-Smooke, and DoublePoints all pass the physical-reference gate from the same 301-node uniform seed. Damped settings are TwoPoint `(0.01,0.01,2.0)`, Pearson `(1e-3,1.4)`, Grcar-Smooke `(1e-3,1e-2,1.4)` (each capped at two refinements), and DoublePoints (four refinements). Measured final node counts/errors are: TwoPoint `424`, `dT=1.799e-2 K`, `dY=1.199e-5`, relative `dq=3.197e-5`; Pearson and Grcar-Smooke `1810`, `dT=1.747e-2 K`, `dY=1.165e-5`, relative `dq=1.993e-5`; DoublePoints `4801`, `dT=1.632e-2 K`, `dY=1.088e-5`, relative `dq=1.906e-5`. CorrectionNewton and Frozen's default, `at_high_morm=1e-6`, and `every_m=1` all pass on the fixed 301-node mesh; all Frozen profiles match the CorrectionNewton root within `1.1e-8`. CorrectionNewton's measured errors are `dT=0.2693 K`, `dY=1.795e-4`, relative `dq=2.471e-4`, closure zero at published precision.

**MF-LOW uniform-mesh control:** Production BVP_sci/Sparse on a 2049-node uniform mesh with the same `2e-7` tolerance and `max_nodes=6000` returns typed `MaxNodesExceeded { max_nodes: 6000 }` and publishes no solution. The clustered-mesh BVP_sci route passes. This reproduces the mesh-distribution limitation through the KiThe facade.

The older convenience-path uniform control, with 2049 input nodes, tolerance `1e-6`, and a much larger `max_nodes=30000`, also returns the typed Sci mesh-limit error `MaxNodesExceeded { max_nodes: 30000 }`. Its formerly unconditional success assertion is now a characterization: if a future RST version solves this uniform case, it must pass the same physical-reference thresholds before the test is green.

**MF-LOW fixed-mesh Damped characterization:** With the old uniform 2049-node seed, no grid refinement, and `max_jac=1` or `50`, Damped returns `MissingResult` / telemetry status `NoResult` after 37 or 77 iterations respectively. This is distinct from adaptive Damped on the explicit clustered seed: DoublePoints with four refinements passes the physical gate, while the tested TwoPoint/Pearson/Grcar settings converge to rejected alternate profiles. The fixed-mesh failure test retains counters and requires any future successful result to pass the physical gate and independent defects.

**API rejection/transaction checks:** Nonmonotone mesh and boundary-inconsistent profile are rejected before a solver run. A failed replacement guess leaves the previously published solution and physical mesh unchanged. The standalone numerical regression separately exercises the RST route/options matrix without KiThe physics objects.

**Interpretation and limits:** RST 0.4.35 explicit-mesh support is integrated into the production facade and demonstrated on MF-LOW/MF-HIGH. On MF-LOW, BVP_sci and Damped/DoublePoints pass the physical-reference gate. CorrectionNewton, Frozen/`at_high_morm`, and Damped TwoPoint/Pearson/Grcar-Smooke return discrete roots but miss the same physical gate; for CorrectionNewton and Frozen the measured error bands identify the same alternate branch. The MF-HIGH matrix passes for all tested production methods/policies. These are narrow results for the stated equations, mesh, seed, and settings, not universal solver rankings. The production API transfer is covered, but release correctness and later release performance campaigns remain outstanding.

**Integrated Debug verification:**

```powershell
cargo test --lib --no-default-features literature_mf1963 -- --include-ignored --nocapture --test-threads=1
```

Result on the current tree: `50 passed, 0 failed, 0 ignored`, finished in `703.27 s`. The run includes the full MF ignored matrix. Two formerly red LOW tests now classify the observed uniform-grid Sci node cap and fixed-mesh Damped `NoResult`; neither treats those outcomes as physical success. Release correctness and performance campaigns have not yet been run.

## MF-06: Production-Facade Parity Against Standalone RST (2026-11-12)

**Name:** The MF-LOW production facade preserves generated equations, scaling,
initial mesh, seed, and selected RST route; it does not worsen the standalone
route's physical acceptance result.

**Command:**

```powershell
cargo test --lib --no-default-features mf1963_mflow_reference_fixtures_report_cross_workflow_drift -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_production_facade_matches_standalone_rst_solver_outputs -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mfhigh_production_facade_matches_standalone_rst_solver_outputs -- --ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_damped_tolerance_and_grid_threshold_sweep -- --ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Do the outcomes measured in the RST-only
`KITHE_MF1963_SOLVER_SELECTION` / `production_mesh_matrix` setup survive the
full `SimpleReactorTask` path? Does KiThe preserve equations and scaling, the
nonuniform mesh, node-aligned initial guess, and route-specific settings? The
decisive control is direct public RST `BVP_api` using the exact equation
expressions generated by KiThe, not only independently re-derived equations.

**Setup:** Seven MF-LOW configurations are compared: BVP_sci/Sparse;
CorrectionNewton/BandedAuto; Damped/TwoPoint `(0.01,0.01,2.0)`;
Damped/Pearson `(1e-3,1.4)`; Damped/Grcar-Smooke `(1e-3,1e-2,1.4)`;
Damped/DoublePoints with four refinements; and Frozen with
`at_high_morm=1e-6`. Each route is run using (1) compact standalone RST
equations, (2) direct RST with the exact equation vector from the KiThe reactor
builder, and (3) the production facade. Route, options, 2049-node clustered
input mesh, and full six-state seed are held fixed. Published profiles are
mapped back to `[Teta,q,C0,J0,C1,J1]` for comparison. The unchanged physical
gate is `dT < 0.5 K`, `dY < 2e-3`, relative heat-flux error `< 2e-3`, closure
`< 2e-8`, and minimum species fraction `>= -1e-12`.

**Equation parity:** The KiThe-generated and compact standalone equations are
mathematically equivalent but differently parenthesized symbolic forms.
Across the shared seed their maximum scaled row discrepancies are
`[0, 1.2463e-11, 4.435e-16, 1.337e-13, 4.435e-16, 1.337e-13]`; the largest
is in the thermal row. This is a small floating-point evaluation difference,
not a missing equation or scale factor. Direct RST with the exact generated
equations isolates the facade handoff from this expression-form difference.

**Two independent LOW reference fixtures:** The RST-only report evaluates its
BVP_sci profile against `fixtures/mf1963_mf_low_six_state_scipy_solution.csv`;
the KiThe physical gate evaluates against the independent
`fixtures/mf1963_mf_low.csv`. These are distinct SciPy artifacts, not one
shared oracle. A new test maps both to common `xi` coordinates and measures
maximum absolute differences of `1.177788894e-3` in `theta`, `1.008443983e-3`
in `p`, and `1.177788894e-5` in `Y_A`; the largest theta difference is at
`xi=11994.060988`. Since `T=300+15*theta`, this is `0.0176668 K`, explaining
the KiThe BVP_sci error against the MF-00 reference. Direct RST using the exact
KiThe-generated equations and the production facade agree on shared seed
nodes to the tight tolerances below, so this delta is not caused by the facade.
The reason the two independently generated SciPy fixtures differ remains open.

| Route | RST-only outcome | Direct RST / facade final nodes | Max scaled state delta on shared input nodes | MF-00 physical gate, direct canonical / facade |
|---|---|---:|---:|---|
| BVP_sci / Sparse | Converged; 3807 nodes; max error vs RST SciPy fixture `1.29e-8` | `3807 / 3985` | `3.37e-12` | PASS / PASS |
| CorrectionNewton / BandedAuto | Converged fixed-grid trapezoid system; profile differs from collocation | `2049 / 2049` | `2.22e-15` | FAIL / FAIL |
| Damped / TwoPoint | Converged, 2 refinements | `2207 / 2207` | `5.39e-12` | FAIL / FAIL |
| Damped / Pearson | Converged, 2 refinements | `4574 / 4574` | `5.67e-12` | FAIL / FAIL |
| Damped / Grcar-Smooke | Converged, 2 refinements | `4574 / 4574` | `2.22e-15` | FAIL / FAIL |
| Damped / DoublePoints | Converged, 4 refinements; closest tested Damped profile to RST SciPy fixture | `32769 / 32769` | `2.06e-13` | PASS / PASS |
| Frozen / `at_high_morm=1e-6` | Converged in 8 iterations / 8 Jacobian refreshes | `2049 / 2049` | `6.37e-12` | FAIL / FAIL |

Direct RST using the exact KiThe-generated equation vector and the facade
solution differ by at most `6.4e-12` scaled state error across all seven
routes; the facade retains all 2049 input nodes. BVP_sci can choose slightly
different adaptive final meshes (`3807` versus `3985` here), so final node
count equality is not an acceptance condition. Damped mesh counts `2207`,
`4574`, and `32769` match the standalone matrix exactly.

RST reports Damped and CorrectionNewton errors in dimensionless `theta`, while
the KiThe gate reports Kelvin: one `theta` unit is `15 K`. The values align:
CorrectionNewton `theta=8.36` maps to `125.4 K`; TwoPoint `1.48` to `22.2 K`;
Pearson/Grcar `0.792` to `11.9 K`; DoublePoints `0.0275` to `0.413 K`. Thus
RST's “Passed” for TwoPoint/Pearson/Grcar means the numerical solve converged;
it does not imply they pass the stricter KiThe MF-00 profile gate. This
resolves the apparent status mismatch without relaxing either criterion.

The RST-only 4097-node CorrectionNewton refinement, its SciPy-profile seed,
and Frozen default/`every_m=1`/CorrectionNewton cross-seeds are additional
experiments, not rows in this seven-route production-facade parity matrix. The
facade parity gate deliberately fixes one caller-owned 2049-node physical seed
for every route; those alternate seed/mesh experiments still need separate
production-facade coverage before claiming parity for them.

**Conclusion:** No evidence was found of equation loss, incorrect scaling,
changed initial guess, or route-setting loss in the production handoff. The
standalone solver outcomes are reproduced: direct-RST and facade
classifications agree, and exact-equation profiles match closely. Against the
KiThe MF-00 physical gate, BVP_sci and Damped/DoublePoints pass; the other
listed routes converge but fail that stricter profile gate. The drift between
the two SciPy reference artifacts remains a separate open item. This is a
Debug correctness/parity story, not a performance result or release signoff.

## MF-07: HIGH/LOW Production Capability and RST Parity Summary (2026-11-12)

**Name:** The production reactor facade can reproduce the independent MF-00
reference in both regimes, and its measured LOW routes agree with direct RST
when equations, seed, mesh, and options are matched.

**Commands:**

```powershell
cargo test --lib --no-default-features mf1963_bvp_sci_matches_mf_high_full_physical_profile -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_bvp_sci_reaches_mf_high_from_independent_smooth_seed -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_production_correction_newton_matches_mf_high_profile -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mfhigh_damped_grid_strategy_matrix_matches_reference_gate -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mfhigh_frozen_refresh_policy_matrix_matches_reference_gate -- --include-ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_reference_fixtures_report_cross_workflow_drift -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_production_facade_matches_standalone_rst_solver_outputs -- --include-ignored --nocapture --test-threads=1
```

**Answers:**

1. **Does the production facade reproduce the semi-analytical reference with
   at least some solvers? Yes.** MF-HIGH passes through BVP_sci, CorrectionNewton,
   all four tested Damped strategies, and all three Frozen refresh policies.
   MF-LOW passes through BVP_sci and Damped/DoublePoints with four refinements.
   “Pass” here means the unchanged independent MF-00 physical-profile gate,
   not just that the nonlinear solver returned a vector.
2. **Does the production facade agree with standalone RST? Yes, for the tested
   apples-to-apples MF-LOW matrix.** Seven routes use the same 2049-node seed,
   settings, and equation graph in direct public RST and the facade. On common
   seed nodes the maximum scaled state difference is `6.4e-12`; every route has
   the same MF-00 physical-gate classification. The matched MF-HIGH parity
   matrix is now also complete; see MF-08. These results are limited to the
   explicitly listed routes, inputs, meshes, and settings.
3. **Which routes solve each regime?** The detailed matrix follows. It keeps
   solver convergence and physical-reference acceptance in separate columns.

**Common reference and gates:** MF-LOW uses `delta=1e-4`, `Pe=12017.324603855763`,
and the 2049-node clustered physical-profile seed. MF-HIGH uses `delta=0.1`,
`Pe=0.5427079731331537`, and a 301-node uniform physical-profile seed; the
BVP_sci control was also started from an independently reconstructed smooth
seed. The production gate compares dimensional profiles with the MF-00
independent SciPy reference and requires `max dT < 0.5 K`, `max dY < 2e-3`,
relative heat-flux error `< 2e-3`, closure `< 2e-8`, and no species fraction
below `-1e-12`. MF-00 is an independent numerical/semi-analytical reference,
not a closed-form solution or a digitization of the 1963 paper.

All production BVP_sci rows use ExprLegacy/Lambdify/Sparse with
StandardFaer. The MF-HIGH BVP_sci facade stories use tolerance `1e-6`, a 20,000
node cap, and 301 input nodes; the separate RST-only HIGH matrix uses tighter
`2e-7` tolerance and obtains 823 nodes, so those mesh counts are not an
apples-to-apples discrepancy. The paired MF-LOW test uses tolerance `2e-7`,
boundary tolerance `2e-7`, and a 6000-node cap.

CorrectionNewton production rows use BandedAuto, absolute and per-state
relative tolerances `2e-7`, max 100 iterations, and a lower Teta bound just
above zero Kelvin. Damped rows use Banded/Lambdify/Trapezoid, absolute and
per-state relative tolerances `2e-7`, max 100 iterations, `max_jac=3`,
`max_damp_iter=10`, damping factor `0.5`, the same positive-Kelvin Teta bound,
and `[-1e6,1e6]` bounds on other states. Frozen rows use Banded/Lambdify/
Trapezoid, tolerance `2e-7`, max 100 iterations, and the same input bounds.

| Regime | Solver / settings | Initial guess | Solver outcome / output nodes | MF-00 reference result |
|---|---|---|---|---|
| LOW | BVP_sci / Sparse, ExprLegacy+Lambdify, `tol=2e-7` | Clustered physical profile, 2049 nodes | Converged; direct RST 3807 nodes, facade 3985 | **PASS**; `dT=1.767e-2 K`, `dY=1.178e-5`, relative `dq=3.433e-5`; RST's own SciPy fixture comparison reports max state error `1.29e-8` |
| LOW | CorrectionNewton / BandedAuto, abs/per-state rel `2e-7` | Same clustered profile, 2049 nodes | Converged fixed-grid root; 2049 nodes | **FAIL**; `dT=1.254e2 K`, `dY=8.360e-2`, relative `dq=1.468e-1`; alternate discrete branch |
| LOW | Damped / TwoPoint `(0.01,0.01,2.0)`, 2 refinements | Same clustered profile, 2049 nodes | Converged; 2207 nodes | **FAIL**; `dT=2.223e1 K`, `dY=1.482e-2`, relative `dq=2.608e-2` |
| LOW | Damped / Pearson `(1e-3,1.4)`, 2 refinements | Same clustered profile, 2049 nodes | Converged; 4574 nodes | **FAIL**; `dT=1.188e1 K`, `dY=7.919e-3`, relative `dq=1.394e-2` |
| LOW | Damped / Grcar-Smooke `(1e-3,1e-2,1.4)`, 2 refinements | Same clustered profile, 2049 nodes | Converged; 4574 nodes | **FAIL**; same measured profile/gate band as Pearson |
| LOW | Damped / DoublePoints, 4 refinements | Same clustered profile, 2049 nodes | Converged; 32769 nodes | **PASS**; `dT=4.130e-1 K`, `dY=2.754e-4`, relative `dq=4.857e-4` |
| LOW | Frozen / default `Frozen_naive` | Same clustered profile, 2049 nodes | RST-only run fails safely on non-finite residual; no solution | Not accepted; this default case is not in the seven-route facade parity test |
| LOW | Frozen / `at_high_morm=1e-6` | Same clustered profile, 2049 nodes | Converged; 2049 nodes; 8 iterations / 8 Jacobian refreshes | **FAIL**; `dT=1.254e2 K`, `dY=8.360e-2`, relative `dq=1.468e-1`; alternate discrete branch |
| LOW | Frozen / `every_m=1` | Same clustered profile, 2049 nodes | RST-only run does not converge within 100 iterations | No accepted profile; not in the seven-route facade parity test |
| HIGH | BVP_sci / Sparse, ExprLegacy+Lambdify, `tol=1e-6` | Uniform physical profile, 301 nodes; also smooth seed | Converged; 520 nodes | **PASS** from both seeds; max `dT=1.397e-2 K`, `dY=9.315e-6`, heat-flux error `4.764 W/m^2` |
| HIGH | CorrectionNewton / BandedAuto, abs/per-state rel `2e-7` | Uniform physical profile, 301 nodes | Converged fixed mesh | **PASS**; `dT=2.693e-1 K`, `dY=1.795e-4`, relative `dq=2.471e-4` |
| HIGH | Damped / TwoPoint `(0.01,0.01,2.0)`, 2 refinements | Same uniform profile, 301 nodes | Converged; 424 nodes | **PASS**; `dT=1.799e-2 K`, `dY=1.199e-5`, relative `dq=3.197e-5` |
| HIGH | Damped / Pearson `(1e-3,1.4)`, 2 refinements | Same uniform profile, 301 nodes | Converged; 1810 nodes | **PASS**; `dT=1.747e-2 K`, `dY=1.165e-5`, relative `dq=1.993e-5` |
| HIGH | Damped / Grcar-Smooke `(1e-3,1e-2,1.4)`, 2 refinements | Same uniform profile, 301 nodes | Converged; 1810 nodes | **PASS**; same measured error band as Pearson |
| HIGH | Damped / DoublePoints, 4 refinements | Same uniform profile, 301 nodes | Converged; 4801 nodes | **PASS**; `dT=1.632e-2 K`, `dY=1.088e-5`, relative `dq=1.906e-5` |
| HIGH | Frozen / default, `at_high_morm=1e-6`, or `every_m=1` | Same uniform profile, 301 nodes | All converge in 3 iterations on fixed mesh | **PASS** for all three; profiles agree with CorrectionNewton within `1.1e-8` |

**Interpretation:** RST's LOW table labels its successful discrete solves
“Passed”; the KiThe production table adds the stricter MF-00 profile gate.
For example, RST's `theta` errors map to KiThe temperature errors by `15 K`
per theta unit: TwoPoint `1.48 -> 22.2 K`, Pearson/Grcar `0.792 -> 11.9 K`,
and DoublePoints `0.0275 -> 0.413 K`. These are consistent outcomes, not
contradictory solver results. The two MF-LOW SciPy reference fixtures do have
a measured drift and remain separately named pending reconciliation.

The standalone-only CorrectionNewton 4097-node run, its SciPy-profile seed,
and Frozen cross-seed experiments are not claimed as production-facade rows.
HIGH direct-pair coverage is recorded separately in MF-08; the seven-row LOW
matrix remains the corresponding LOW evidence.

## MF-08: MF-HIGH Direct-RST / Production-Facade Parity (2026-11-12)

**Name:** The production facade preserves the MF-HIGH equation graph, mesh,
initial profile, selected solver options, and resulting physical profile when
compared with direct public RST `BVP_api` execution.

**Command:**

```powershell
cargo test --lib --no-default-features mf1963_mfhigh_production_facade_matches_standalone_rst_solver_outputs -- --ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Does MF-HIGH's successful production-facade matrix
also agree with direct RST if both receive the exact KiThe-generated equations
and identical mesh, node-aligned seed, and route settings? A facade/direct
mismatch would indicate altered equations, scaling, guess, or configuration.

**Setup:** `delta=0.1`, `Pe=0.5427079731331537`, 301-node uniform mesh, and
the KiThe physical-profile seed. For each route, the test builds equations
through `SimpleReactorTask::setup_bvp()`, verifies equation parity, then
solves once through direct RST `BVP_api` and once through the production
facade. Routes are BVP_sci/Sparse (`tol=1e-6`), CorrectionNewton/BandedAuto,
Damped TwoPoint `(0.01,0.01,2.0)`/2 refinements, Pearson `(1e-3,1.4)`/2,
Grcar-Smooke `(1e-3,1e-2,1.4)`/2, DoublePoints/4, and Frozen default,
`at_high_morm=1e-6`, and `every_m=1`. All other settings match the MF-07
production cases. Both solutions are checked against the unchanged MF-00
physical gate.

**Results:** KiThe equations and the compact RST equations agree at the probe
points with maximum scaled row discrepancies
`[0, 6.43e-15, 3.57e-16, 3.95e-17, 3.57e-16, 3.95e-17]`. All nine direct-RST
and facade pairs have identical final node counts and pass the MF-00 gate.
Maximum scaled componentwise state deltas on common mesh nodes:

| Route | Direct / facade nodes | Max scaled state delta | Direct / facade physical gate |
|---|---:|---:|---|
| BVP_sci / Sparse, `tol=1e-6` | `520 / 520` | `1.46e-15` | PASS / PASS |
| CorrectionNewton / BandedAuto | `301 / 301` | `1.57e-15` | PASS / PASS |
| Damped / TwoPoint `(0.01,0.01,2.0)`, 2 | `424 / 424` | `9.72e-16` | PASS / PASS |
| Damped / Pearson `(1e-3,1.4)`, 2 | `1810 / 1810` | `1.15e-15` | PASS / PASS |
| Damped / Grcar-Smooke `(1e-3,1e-2,1.4)`, 2 | `1810 / 1810` | `3.33e-16` | PASS / PASS |
| Damped / DoublePoints, 4 | `4801 / 4801` | `3.58e-16` | PASS / PASS |
| Frozen / default | `301 / 301` | `4.21e-16` | PASS / PASS |
| Frozen / `at_high_morm=1e-6` | `301 / 301` | `1.40e-15` | PASS / PASS |
| Frozen / `every_m=1` | `301 / 301` | `4.21e-16` | PASS / PASS |

**Conclusion:** The binary answer for this matched MF-HIGH matrix is **yes**:
the production facade reproduces direct RST outcomes without measurable
degradation. Together with MF-06's seven LOW route pairs, direct-RST/facade
parity is now demonstrated for both regimes. This does not claim parity for
unlisted seeds, meshes, solver options, or other physical models.

## MF-09: MF-LOW Damped Tolerance and Grid-Indicator Sweep (2026-11-12)

**Name:** On the fixed MF-LOW clustered seed, tightening nonlinear tolerances
is distinguished from changing adaptive-grid indicators; any improvement is
measured against the same physical reference and acceptance gate.

**Command:**

```powershell
cargo test --lib --no-default-features mf1963_mflow_damped_tolerance_and_grid_threshold_sweep -- --ignored --nocapture --test-threads=1
```

**Question / hypothesis:** The RST-only matrix called LOW TwoPoint/Pearson/
Grcar results converged, while the independent physical-profile gate showed
errors of `22.2 K` and `11.9 K`. Can tighter Damped `atol/rtol`, or more
selective grid refinement, bring these adaptive strategies into the physical
acceptance envelope?

**Setup:** Production facade, MF-LOW (`delta=1e-4`,
`Pe=12017.324603855763`), same 2049-node clustered physical seed, Banded /
Lambdify / Trapezoid, bounds unchanged, `max_jac=3`, `max_damp_iter=10`,
`damp_factor=0.5`, and up to four refinements. For each strategy, run the
baseline coefficients at `atol=rtol=2e-7`, the same coefficients at
`atol=rtol=1e-9`, and stricter indicator coefficients at `2e-7`:
TwoPoint `(0.005,0.005,1.4)`, Pearson `(1e-4,1.4)`, and Grcar-Smooke
`(1e-4,1e-3,1.4)`. The MF-00 gate is unchanged: `dT<0.5 K`, `dY<2e-3`,
relative heat-flux error `<2e-3`, closure `<2e-8`, and nonnegative species
fractions within `1e-12`.

**Results:** All runs produced finite profiles and passed the independent
discrete trapezoid-defect checks. Tightening nonlinear tolerance from `2e-7`
to `1e-9` kept the same mesh sizes and changed reference errors only in the
last reported digits. Stricter indicator coefficients improved each
strategy's normalized worst-gate score, but none passed the physical gate:

| Strategy | Variant | `atol=rtol` | Final nodes | max `dT` (K) | max `dY` | relative `dq` | MF-00 gate |
|---|---|---:|---:|---:|---:|---:|---|
| TwoPoint `(0.01,0.01,2.0)` | baseline | `2e-7` | `2209` | `21.9878` | `1.46586e-2` | `2.58024e-2` | FAIL |
| TwoPoint `(0.01,0.01,2.0)` | tighter tolerance | `1e-9` | `2209` | `21.9879` | `1.46586e-2` | `2.58025e-2` | FAIL |
| TwoPoint `(0.005,0.005,1.4)` | stricter indicator | `2e-7` | `2462` | `14.5352` | `9.69015e-3` | `1.70575e-2` | FAIL |
| Pearson `(1e-3,1.4)` | baseline | `2e-7` | `4690` | `11.4906` | `7.66043e-3` | `1.34853e-2` | FAIL |
| Pearson `(1e-3,1.4)` | tighter tolerance | `1e-9` | `4690` | `11.4906` | `7.66042e-3` | `1.34853e-2` | FAIL |
| Pearson `(1e-4,1.4)` | stricter indicator | `2e-7` | `31792` | `6.89875` | `4.59916e-3` | `8.09640e-3` | FAIL |
| Grcar-Smooke `(1e-3,1e-2,1.4)` | baseline | `2e-7` | `4690` | `11.4906` | `7.66043e-3` | `1.34853e-2` | FAIL |
| Grcar-Smooke `(1e-3,1e-2,1.4)` | tighter tolerance | `1e-9` | `4690` | `11.4906` | `7.66042e-3` | `1.34853e-2` | FAIL |
| Grcar-Smooke `(1e-4,1e-3,1.4)` | stricter indicator | `2e-7` | `31792` | `6.89875` | `4.59916e-3` | `8.09640e-3` | FAIL |

For context, the existing Damped/DoublePoints LOW story passes at `32769`
nodes with `dT=0.413 K`, `dY=2.754e-4`, and relative `dq=4.857e-4`. Thus
node count alone does not explain the adaptive-strategy gap: Pearson/Grcar
used a comparable number of nodes but remained outside the physical gate.

**Conclusion:** The candidate adaptive indicators are directionally better,
but not “almost passing” under KiThe's unchanged physical gate; their best
tested temperature error remains about 13.8 times the limit and the species/
flux errors also exceed their limits. Nonlinear tolerance tightening is not
the useful control in this experiment. Keep DoublePoints as the demonstrated
LOW Damped pass; further adaptive-strategy tuning remains open and must be
judged by physical profile error, not solver convergence or grid size alone.

## MF-10: Isolated MF-LOW Gradient-Refinement Coefficients (2026-11-12)

**Name:** Tightening the gradient indicator adds points in steep-profile
regions only when that criterion is active; value-jump and gradient
coefficients must be varied independently to identify which one improves MF-LOW.

**Commands:**

```powershell
cargo test --lib --no-default-features mf1963_mflow_damped_twopoint_coefficient_sweep -- --ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_damped_grcar_coefficient_sweep -- --ignored --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_mflow_damped_pearson_coefficient_sweep -- --ignored --nocapture --test-threads=1
```

**Question / hypothesis:** After nonlinear `atol/rtol` proved ineffective,
does a stricter gradient-change threshold `g` direct Damped refinement toward
the MF-LOW reaction layer? For Pearson, which has no separate `g` coefficient,
does lowering its value-jump threshold `d` have a similar effect?

**Setup:** Same MF-LOW production facade, 2049-node clustered physical seed,
Banded/Lambdify/Trapezoid, `atol=rtol=2e-7`, four maximum refinements, and
unchanged physical bounds. Each coefficient is isolated against its baseline:
TwoPoint holds `d=0.01` and changes only `g`; Grcar-Smooke holds `d=1e-3`
and changes only `g`; Pearson changes only `d` and keeps its buffer parameter
`C=1.4`. `transition_nodes` counts published nodes whose `C0` lies in
`[0.05,0.95]`, a diagnostic proxy for the species transition layer. The same
MF-00 physical gate remains in force.

**Results:**

| Strategy | Coefficient change | Nodes | Transition nodes | max `dT` (K) | max `dY` | relative `dq` | Gate |
|---|---|---:|---:|---:|---:|---:|---|
| TwoPoint | baseline `d=.01,g=.01` | `2209` | `1008` | `21.9878` | `1.46586e-2` | `2.58024e-2` | FAIL |
| TwoPoint | `g=.001` (10x tighter), same `d` | `2209` | `1008` | `21.9878` | `1.46586e-2` | `2.58024e-2` | FAIL |
| TwoPoint | `g=.0001` (100x tighter), same `d` | `2209` | `1008` | `21.9878` | `1.46586e-2` | `2.58024e-2` | FAIL |
| Grcar-Smooke | baseline `d=1e-3,g=1e-2` | `4690` | `3141` | `11.4906` | `7.66043e-3` | `1.34853e-2` | FAIL |
| Grcar-Smooke | `g=1e-3`, same `d` | `5322` | `3766` | `11.4580` | `7.63869e-3` | `1.34471e-2` | FAIL |
| Grcar-Smooke | `g=1e-4`, same `d` | `21312` | `18623` | `11.3521` | `7.56809e-3` | `1.33228e-2` | FAIL |
| Pearson | baseline `d=1e-3,C=1.4` | `4690` | `3141` | `11.4906` | `7.66043e-3` | `1.34853e-2` | FAIL |
| Pearson | `d=5e-5`, same `C` | `62521` | `54475` | `4.59165` | `3.06110e-3` | `5.38891e-3` | FAIL |

An additional Pearson trial at `d=1e-5` did not complete: the Debug test
process reached approximately `12.6 GB` resident memory, so it was stopped
before a result was published. No accuracy claim is made for that run.

**Conclusion:** The isolated `g` coefficient is not a useful control for
TwoPoint on this fixture. Grcar-Smooke responds to lower `g` with denser
transition sampling, but the profile error improves only slightly even at
21,312 nodes. Pearson's value-jump threshold gives the largest measured gain,
but needs 62,521 nodes and still misses the physical gate by a wide margin.
Do not lower thresholds blindly: inspect refinement marks and impose explicit
node/memory budgets before further sweeps. Solver tolerance remains unchanged.

