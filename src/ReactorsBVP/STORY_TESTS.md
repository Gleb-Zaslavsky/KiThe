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

## Story: Compact Reactor Frozen No-Result and Boundary Reconstruction

**Name:** Frozen preserves the prior profile on `NoResult`; Naive publishes a full profile with restored endpoints.

**Command:**

```powershell
cargo test --lib --no-default-features story_compact_reactor_frozen_no_result_and_naive_boundary_reconstruction -- --nocapture
```

**Question / hypothesis:** Does the updated Frozen contract distinguish a
native `NoResult` from a published profile, preserve the last accepted Damped
solution on that failure, and reconstruct prescribed endpoints when Naive
does publish a profile?

**What is checked:** The production Damped route first publishes the compact
HMX profile. Frozen and Naive each receive that converged 33-node profile as
their initial guess. The Frozen strategy returns typed
`BvpApiError::MissingResult`; its Counters snapshot reports `NoResult`, and the
previous Damped solution and mesh remain published. Naive publishes a finite
full nodal profile whose configured endpoints equal their prescribed values;
its route-tagged Counters snapshot reports `Converged`.

**Observed result:** Pass on RustedSciThe 0.4.30. Frozen reports no result and
does not replace the Damped publication. Naive's reconstructed endpoints pass
the facade's boundary gate. This test checks the public profile shape and
endpoint contract, not internal discretization accuracy at every node.

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

**Observed result:** Pending release campaign. Debug targeted runs pass for the 32-interval route matrices and both 96-interval cases.

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
This analytic control is deliberately separate from compact HMX,
where Frozen can report `NoResult` and Naive's endpoint accuracy remains outside the production
acceptance tolerance.

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

**Observed result:** Debug story passes. For Damped, Sequential recorded `0 parallel / 5 sequential` dispatches, Auto and Forced Parallel each recorded `5 / 0`. For BVP_sci, Sequential recorded `0 / 513`, Auto selected sequential (`0 / 513`), and Forced Parallel recorded `513 / 0`. All six runs published matching profiles within `2e-5`. These counts confirm that `Auto` is a policy decision, not a synonym for forced parallel, and can differ by solver family. Benchmark timings for KiThe setup and callback-policy wall time await the next release campaign.

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
The same saved task is then solved through direct `solve_from_file` and the GUI calculation path;
the complete postprocessed mesh and solution matrix handed to the plot window must agree componentwise
within `1e-8 * max(|reference|, 1)` (mesh `1e-12`).

**Observed result:** Debug GUI story passes and both routes produce matching physical profiles. The
actual GUI widget edit survives Save As and reload, while the source file remains unchanged. This catches a prior
serialization defect where integral-valued `Float`s were written without a decimal point and
reparsed as `Integer`, causing otherwise valid solver inputs to fail validation. File I/O is not
performance-measured.

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
