# SimpleReactorsIVP Story Tests

## Story: StiffAPI Family Routing for Condensed IVP

**Question / hypothesis:** Can the condensed reactor build its physical RHS in KiThe and solve
through the RST StiffAPI using LSODE2, standalone BDF, or Radau, without confusing the LSODE2
controller with standalone solver families? Do the new document route and telemetry reporting retain
the selected family and native diagnostics?

**Command:**

```powershell
cargo test --lib --no-default-features -j 2 SimpleReactorsIVP:: -- --nocapture
```

**Fixture and route:** The canonical small H2 -> 2H condensed-reactor fixture has the same physical
parameters and initial state across the solver-family stories. The LSODE2 route retains its
Auto/BDF/Adams controller selection. Standalone BDF and Radau use the same KiThe-built symbolic RHS
through `StiffSolverRequest`. A separate hybrid-document story parses only solver settings through
RST, applies KiThe physical inputs locally, and prepares the resulting caller-built BDF request with
`HybridStiffSettings::prepare`.

**Acceptance checks:** Published trajectories have the expected family and state shape, finite
values, monotonic coordinate, requested interval endpoints, and the expected H2/H conversion trend.
Detailed BDF and Radau telemetry is exposed in report rows with native counters and timing names.
The public config stores LSODE2 controller/native execution, BDF step budget, and Radau step/Newton
limits in disjoint tagged variants. GUI and task-document round-trips preserve the selected variant
and the distinct RST `max_iterations` interpretation. Radau Dense, Sparse, and BandedAuto are
exercised by the matrix story below for both symbolic frontends.

**Observed result:** Current targeted run: 73 SimpleReactorsIVP tests passed, 0 failed, 0 ignored. The hybrid BDF document route, standalone
BDF and Dense Radau solve stories, family-aware task round-trips, and native telemetry row assertions
all passed. The GUI integration filter passed 30/30 tests after changing the demo fixture from
placeholder species `A`/`B` to valid `H2`/`H` chemistry.

**Conclusion:** The StiffAPI family handoff is working for the tested LSODE2, standalone BDF, and
Dense Radau paths. This is correctness evidence for the canonical small fixture, not a claim that
the families have equivalent adaptive meshes, performance, or broad physical parity.

**Follow-up limitations:** Solver capability rejection, StiffAPI preparation errors, and native
integration errors now retain distinct types and source chains; parser, postprocessing, and
publication error boundaries still need a broader audit. Sparse Radau correctness is established
for the bounded story fixture only; broader mechanism coverage remains open.

## Story: Callback Execution Policy Parity

**Question / hypothesis:** On the largest common bounded IVP fixture that completes successfully
for all three StiffAPI families, do `Sequential`, `Auto`, and forced `Parallel` callback policies
preserve the published endpoint? Where RST exposes dispatch counters, does forced parallelism
actually dispatch parallel callback work?

**Commands:**

```powershell
cargo test --lib --no-default-features callback_execution_policies_preserve_largest_solved_ivp_fixture_and_report_dispatches -- --nocapture --test-threads=1
cargo bench --bench reactor_ivp_baseline --no-default-features -j 1 -- --noplot
```

**Fixture and comparison:** Four-state H2 -> 2H condensed IVP, 900 K, domain `[0, 0.02]`,
`max_step = 1e-3`, `rtol = 1e-6`, `atol = 1e-8`, ExprLegacy frontend and Dense matrix layout.
The test runs LSODE2, standalone BDF, and Radau with all three policies. The benchmark has a
separate Criterion group for the same nine solver/policy pairs; telemetry is Off in timed samples.
Detailed policy counters/timings are reported separately and are not mixed into wall-time samples.

**Acceptance checks:** Every solve finishes with finite output and matching state shape. The
endpoint must agree with the Sequential reference within `2e-6 * (1 + |reference|)`. Where the
telemetry exposes the configured policy, it must match the request; otherwise dispatch behavior is
checked directly. Where callback dispatch counters are available, Sequential must report no
parallel dispatch, and forced Parallel must report at least one.

**Observed result:** Refreshed focused run on RustedSciThe 0.4.31 passed. LSODE2 reported
Sequential/Auto/forced-Parallel dispatches as `0/64`, `0/64`, and `64/0` (parallel/sequential),
with `64` residual and `40` Jacobian callbacks. Standalone BDF reported `0/43`, `0/43`, and
`43/0`, with `42` residual and `1` Jacobian callback for each policy; BDF native solver counters
also report `42` residual calls and `1` Jacobian call. Radau reported `0/143`, `0/143`, and
`143/0`, with native counters of `142` residual calls and `1` Jacobian call. Radau exposes its own
native solver counters rather than the symbolic residual/Jacobian callback fields. Endpoint parity
passed for all nine routes.

**Conclusion:** The 0.4.31 post-solve telemetry snapshot reaches the KiThe facade for standalone
BDF: its callback counters are nonzero and reflect the selected dispatch policy. Radau's distinct
telemetry schema is not treated as missing BDF-style callback attribution. This is correctness and
telemetry evidence only; release performance timings remain for the scheduled benchmark campaign.

### Historical Sparse Radau Jacobian Shape Diagnostic (Superseded)

**Command:**

```powershell
cargo test --lib --no-default-features radau_matrix_layouts_match_for_both_symbolic_frontends -- --nocapture --test-threads=1
```

The prior 0.4.29 failure and RST shape analysis are retained as historical evidence only. The old
ignored reproducer was removed after its expected failure stopped representing the current route.
On the current RST, the production guard is removed and Sparse Radau is covered by a bounded
physical-profile matrix against Dense and BandedAuto for ExprLegacy and AtomView. See the active
story below; do not rerun the old expected-failure assertion.

## Story: Failed Solve Retains Partial Telemetry

**Question / hypothesis:** When an instrumented standalone BDF run exhausts its one-step budget,
does the returned typed native failure retain the RST telemetry gathered before failure, while
keeping the incomplete trajectory unpublished?

**Command:**

```powershell
cargo test --lib --no-default-features -j 2 failed_standalone_bdf_solve_retains_partial_native_telemetry -- --nocapture
```

**Fixture and route:** Canonical H2 -> 2H condensed fixture, standalone BDF, `max_steps = 1`, and
`Detailed` telemetry. The interval and physical inputs are unchanged; only the intentionally
insufficient step budget forces a native runtime failure.

**Acceptance checks:** The result is `IvpError::SolverExecution { error, telemetry }` with native solver kind BDF and a
`StiffTelemetrySnapshot::Bdf` in Detailed mode with at least one residual callback. No solve
snapshot is published.

**Observed result:** The forced failure returned the expected typed error and partial telemetry;
the task retained no successful snapshot.

**Conclusion:** KiThe can expose partial attempt telemetry when the selected RST solver provides it.
Telemetry remains `None` when disabled or when the native route does not make a failed-attempt
snapshot available.

## Story: Typed StiffAPI Failure Boundaries

**Question / hypothesis:** Does KiThe preserve the distinction between an unsupported solver
capability, StiffAPI preparation failure, and native integration failure instead of flattening all
three into `CalculationError(String)`? Does a capability rejection happen before integration?

**Command:**

```powershell
cargo test --lib --no-default-features -j 2 SimpleReactorsIVP::reactor_ivp_story_tests -- --nocapture
```

**Acceptance checks:** StiffAPI preparation and native solve failures are carried by
`IvpError::SolverPreparation(StiffApiError)` and `IvpError::SolverExecution { error, telemetry }`;
`Error::source()` exposes the original RST error. No successful trajectory is published after an
error.

**Observed result:** The focused story suite passes. GUI status can still display the error via
`Display`, while callers inspecting the returned `IvpError` retain access to its source chain.

**Conclusion:** Native StiffAPI preparation and solve failures no longer lose their RST error type
at the primary facade boundary. This does not complete typed errors for parser, postprocessing, or
legacy compatibility APIs.

## Story: Radau Matrix and Frontend Parity

**Name:** Bounded Radau trajectories agree across Dense, Sparse, and BandedAuto for both symbolic frontends.

**Command:**

```powershell
cargo test --lib --no-default-features radau_matrix_layouts_match_for_both_symbolic_frontends -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the repaired Sparse Radau route now integrate a bounded KiThe-built
reactor system, and do its published trajectories match the Dense reference when frontend and
physical inputs are held fixed?

**What is checked:** The low-rate H2 -> 2H scenario uses the same initial state, thermal inputs,
interval, and tolerances across ExprLegacy/AtomView and Dense/Sparse/BandedAuto. Every trajectory
must finish, remain finite, keep physical temperature in `[0, 3000] K`, keep species fractions in
`[0, 1]`, show measurable reaction progress, and conserve the H2/H fraction sum. Full profile
values are compared componentwise to the first route with a `2e-6` maximum absolute tolerance.

**Observed result:** Pass in the current Debug run. All six routes completed, met the reaction
progress/conservation checks, and matched the reference profile within `2e-6`; the former Sparse
callback-shape error did not recur.

**Conclusion:** Sparse Radau is enabled in the KiThe adapter for the tested bounded scenario.
This does not imply every Radau workload or every sparse pattern is guaranteed; the release matrix
benchmark adds a performance observation axis, while this story remains the correctness gate.

## Story: Multi-Step AtomView Frontend Parity

**Name:** AtomView and ExprLegacy produce matching bounded profiles for a consecutive two-reaction network.

**Command:**

```powershell
cargo test --lib --no-default-features multistep_oxidation_expr_and_atomview_publish_matching_bounded_profiles -- --nocapture --test-threads=1
```

**Question / hypothesis:** After the Expr-to-Atom conversion fix, does the larger multi-step symbolic
expression path that previously needed manual diagnostics complete a real integration under both
frontends and preserve bounded physical state?

**What is checked:** The same two-step synthetic oxidation network, mixture, thermal inputs, interval,
and tolerances are solved through LSODE2 using ExprLegacy and AtomView. Both runs must finish, publish
finite profiles, keep temperature in `[0, 3000] K` and species fractions in `[0, 1]`, retain KiThe
stage telemetry, and match componentwise within `1e-7` absolute error.

**Observed result:** Pass in Debug on the current RST lockfile. Both integrations finished with
bounded state; the maximum absolute profile delta was `0` at the published output points.

**Conclusion:** This closes the former preparation/conversion diagnostic as an executable correctness
story rather than logging caught failures and continuing.

## Story: KiThe IVP Stage Attribution

**Question / hypothesis:** Does IVP telemetry distinguish KiThe-side setup and publication work
from RST request preparation and numerical integration, without reading the clock when telemetry
is disabled or mixing counters with timings in the single-channel modes?

**Commands:**

```powershell
cargo test --lib --no-default-features -j 2 stage_telemetry_modes_keep_counter_and_clock_costs_separate -- --nocapture
cargo test --lib --no-default-features -j 2 test_canonical_condensed_story_fixture_solves_with_standalone_bdf -- --nocapture
```

**Fixture and route:** The stage-mode contract uses a zero-work operation to verify `Off`,
`Counters`, and `Timings` semantics. The integration check uses the canonical condensed H2 -> 2H
fixture through standalone BDF with `Detailed` telemetry.

**Acceptance checks:** `Off` has no stage records; `Counters` has call counts but no elapsed values;
`Timings` has elapsed values but no counters. A detailed BDF result reports KiThe stages
`contract_validation`, `scaling_processing`, `kinetic_processing`, `peclet_numbers`,
`equation_assembly`, `pre_solve_validation`, `solver_request_build`, `rst_solver_preparation`,
`native_integration`, `trajectory_materialization`, `candidate_validation`, `conservation`,
`solution_publication`, `postprocessing`, and end-to-end `total_workflow`, separately from RST's
native timing rows.

**Observed result:** Both focused tests passed. The existing canonical BDF solve remains accepted;
the report exposes separate KiThe stage rows and RST native timing/counter rows.

**Conclusion:** The previous opaque preparation/integration boundary is now split for one-shot IVP
reports, with mode semantics matching the shared telemetry policy. This is not yet a release
performance baseline: cold versus repeated preparation and end-to-end wall-time reconciliation
remain open. Failed-stage partial attribution is also not yet attached to errors.

## Story: Multistep Network Assembly Coverage

**Question / hypothesis:** Can the condensed IVP setup assemble a small atom-balanced consecutive
reaction network with an intermediate, rather than only a one-step conversion?

**Command:**

```powershell
cargo test --lib --no-default-features -j 2 multistep_oxidation_coverage_fixture_assembles_balanced_network -- --nocapture
```

**Fixture and route:** Synthetic methane oxidation coverage network:
`2CH4 + 3O2 -> 2CO + 4H2O`, followed by `2CO + O2 -> 2CO2`. The reaction
stoichiometry is atom-balanced; rate and heat parameters are synthetic and are not an
experimental mechanism. The test exercises KiThe parsing, kinetic preprocessing, and RHS assembly.

**Acceptance checks:** Setup succeeds with seven canonical states (`Teta`, `q`, and five species),
seven equations, and two reaction-rate expressions.

**Observed result:** The assembly story passed. The full-solve diagnostic localized the AtomView
panic to LSODE2 lazy native callback setup, before the first residual/Jacobian callback. KiThe setup,
request construction, and StiffAPI preparation pass for both frontends. The failing row is the
heat-flux RHS (`rhs[1]`); its second-reaction heat term contains the minimal failing subtree:

```text
-1790971310706.9194 * C1 * (48197.07247411639 * C2)^2
```

Each factor converts to AtomView independently; combining them panics at
`Coefficient::reduce` line 114 in RustedSciThe 0.4.29. The isolated diagnostic command is:

```powershell
cargo test --lib --no-default-features -j 2 multistep_fixture_localizes_atomview_conversion_overflow -- --include-ignored --nocapture --test-threads=1
```

An ExprLegacy end-to-end attempt with the original high-order synthetic rates was stopped after
taking too long, so this story makes no full-trajectory acceptance claim.

**Conclusion and limits:** This closes multi-step assembly coverage and localizes an upstream
AtomView coefficient-reduction defect for a generated heat-flux expression. It does not demonstrate
a converged physical trajectory, solver parity, validated kinetic data, or a performance workload.
An exploratory peroxide-decomposition profile was rejected as a correctness fixture because the
published trajectory contained negative species values; it is not included in the solver matrix and
does not establish a solver regression. A second validated multi-step physical oracle is still needed.

## Story: Bounded IVP Frontend Parity

**Question / hypothesis:** For two bounded low-rate condensed H2 -> 2H inlet scenarios, do
ExprLegacy and AtomView publish matching endpoints through LSODE2, standalone BDF, and Radau when
the physical inputs, initial state, integration interval, tolerance, and Dense matrix layout are
held fixed? Do the published temperature and species fractions stay within basic physical bounds?

**Command:**

```powershell
cargo test --lib --no-default-features -j 2 symbolic_frontends_match_endpoint_for_each_solver_family -- --nocapture
```

**Fixture and route:** Synthetic elementary H2 -> 2H reaction with Arrhenius parameters
`A=1e-3`, `n=0`, `E=8e3`, zero thermal effect, `q=0`. The two matched inlet states are
`520 K / (H2=0.85, H=0.15)` and `900 K / (H2=0.70, H=0.30)`. The common solver configuration
is Dense, `x=[0,0.02]`, max step `1e-3`, `rtol=1e-6`, `atol=1e-8`. These are numerical coverage
parameters, not experimental kinetics.

**Acceptance checks:** All endpoint components are finite and bounded by 10 in solver coordinates;
the dimensional temperature profile is in `[0,3000] K`, every species mass fraction across the
trajectory is in `[0,1]`, both species evolve measurably, and their mass fractions sum to one within
`1e-12`. Within each of the six scenario/family pairs, every endpoint component must match between
frontends within `1e-8 * (1 + abs(reference_component))`.

**Observed result:** Both scenarios passed for LSODE2, BDF, and Radau. The maximum absolute
endpoint delta was `0` for each ExprLegacy/AtomView pair in this run. The independent physical
smoke test below also passed the separate two-step network. These are bounded synthetic controls,
not cross-family equivalence or an independent physical oracle.

## Story: IVP Physical Plausibility Smoke

**Question / hypothesis:** Before using additional IVP scenarios in solver/performance coverage,
does the complete published trajectory remain within conservative physical bounds: temperature
from 0 to 3000 K and each species mass fraction from 0 to 1?

**Command:**

```powershell
cargo test --lib --no-default-features -j 2 physical_smoke_scenarios_keep_temperature_and_mass_fractions_bounded -- --nocapture --test-threads=1
```

**Fixture and route:** The bounded H2 -> 2H baseline at 520 K and a two-step synthetic network,
`2CH4 + 3O2 -> 2CO + 4H2O`, followed by `2CO + O2 -> 2CO2`, with low synthetic rates, positive
initial fractions, LSODE2/ExprLegacy/Dense, and a short interval. These rate constants are smoke
parameters only; this is not an experimentally validated mechanism.

**Acceptance checks:** Every dimensional temperature sample is finite and in `[0,3000] K`. Every
species sample is finite and in `[0,1]`. Both solves must finish successfully. The test intentionally
does not promote this fixture to an external physical oracle or broad solver-family gate.

**Observed result:** Both cases passed; the observed temperature range was `520.000..520.000 K`
for each run. The two-step trajectory published seven state variables and five output points.

**Conclusion and limits:** These assertions reject obvious nonphysical blow-up before benchmarking.
They do not prove kinetic accuracy or elemental-conservation quality for the synthetic network.

## Story: Multistep Elemental Balance Across IVP Families

**Question / hypothesis:** Does the balanced two-step oxidation network preserve its independently
computed C/H/O elemental mass inventory over the full trajectory when run with LSODE2, standalone
BDF, or Radau and either ExprLegacy or AtomView?

**Command:**

```powershell
cargo test --lib --no-default-features -j 1 multistep_oxidation_preserves_elemental_inventory_across_solver_families_and_frontends -- --nocapture --test-threads=1
```

**Fixture and route:** Synthetic atom-balanced reactions `2CH4 + 3O2 -> 2CO + 4H2O` and
`2CO + O2 -> 2CO2`, with the same initial mixture, thermal inputs, short interval, Dense layout,
and tolerances for all six family/frontend combinations. Deliberately low rates keep this a bounded
conservation/adapter fixture, not a validated combustion mechanism.

**Acceptance checks:** Each run finishes; dimensional temperature stays in `[0,3000] K`; species
mass fractions remain finite and in `[0,1]`; and C/H/O elemental mass fractions are recomputed
from the published species states and molecular weights and remain within `2e-7` absolute/relative
scale of their initial values at every native trajectory point.

**Observed result:** All six combinations passed. LSODE2 returned 5 points, BDF 12, and Radau 11
for each frontend; temperature remained at 520 K to six printed decimals. The elemental check is
performed per route, so this also checks conservation despite each family using its own adaptive
mesh.

**Conclusion and limits:** The multistep physical-consistency matrix is broader and now checks an
independent elemental invariant, not only finiteness and bounds. Kinetic rates remain synthetic and
near-quiescent; an independently validated medium multi-step kinetic oracle is still required for
claims about reaction accuracy or thermal feedback.

## Story: Atomic IVP Result Publication

**Question / hypothesis:** If conversion of a validated solver trajectory into dimensional output
fails, are the previously published mesh, solution, and solve metadata left untouched?

**Command:**

```powershell
cargo test --lib --no-default-features -j 1 failed_candidate_postprocessing_leaves_published_state_unchanged -- --nocapture --test-threads=1
```

**Fixture and route:** A previously published two-point state is retained while a candidate
temperature conversion intentionally overflows during the fallible postprocessing preparation
step. The test exercises the same pure candidate-postprocessing helper used before the production
commit boundary, including the solve's conservation/postprocessing preparation stage.

**Acceptance checks:** Postprocessing returns an error and leaves the stored mesh, solution, latest
solve snapshot, and conservation report unchanged. Production solve now computes conservation and
dimensional output from local candidate values and only mutates published result fields after
postprocessing succeeds.

**Observed result:** The regression passed. No solver or benchmark path was changed by this story.

**Conclusion and limits:** The result publication boundary is atomic with respect to fallible
trajectory postprocessing. GUI callers still deliberately clear old results when edited inputs are
submitted, as documented by the GUI lifecycle contract.

## Story: PERF-IVP-01 Release Baseline and Lifecycle Matrix

**Question / hypothesis:** Across two bounded inlet conditions, how do LSODE2, standalone BDF, and
Radau compare across ExprLegacy and AtomView on identical Dense routes? How much is spent in KiThe
setup versus RST integration, and what is saved when the prepared RST solver is reused? For LSODE2,
how does AOT compare with Lambdify for Auto/BDF/Adams controllers and both symbolic frontends?

**Command:**

```powershell
cargo bench --bench reactor_ivp_baseline --no-default-features -j 2 -- --noplot --sample-size 20 --warm-up-time 0.2 --measurement-time 0.5
```

For the optional compiler-backed axis, run with `KITHE_RUN_IVP_AOT_BENCH=1` and `tcc` on `PATH`.

Release profile with the repository's configured LTO. Criterion measures with telemetry Off;
three separate Detailed runs produce the stage and RST diagnostic tables. The measured solve call
starts after construction of `SimpleReactorTask` and includes KiThe setup, request preparation,
integration, validation, publication, and postprocessing.

**Fixture and route:** Two bounded low-rate synthetic H2 -> 2H scenarios (`520 K, H2/H=0.85/0.15`
and `900 K, H2/H=0.70/0.30`), four state variables, common interval and tolerances. Twenty Criterion
case/routes run: both scenarios x LSODE2/BDF Dense and Radau Dense/Sparse/Banded x ExprLegacy/AtomView.
Separate lifecycle rows compare same-initial-state repeats and prepared initial-composition sweeps
for LSODE2 and standalone BDF. The AOT axis is LSODE2-only because the current KiThe adapter
explicitly rejects AOT for standalone BDF and Radau; it covers Auto/BDF/Adams controller modes,
both symbolic frontends, and both bounded scenarios. AOT rows include setup, problem construction,
RST preparation, integration, native detailed telemetry, physical bounds, and endpoint parity with
Lambdify. Unsupported routes are not timed.

**Acceptance checks:** Each detailed preflight run must return a finished status, finite values,
temperature in `[0,3000] K`, and all species fractions in `[0,1]`.
The paired frontend correctness gate is the preceding story; this benchmark is measurement evidence,
not a performance gate. Criterion intervals below are the observed ranges from 20 samples and a
0.5 s target measurement per route; short runs and machine load limit their precision.

**Historical Dense-only Criterion medians (release, 20 samples):**

| Inlet | Solver / frontend | 95% estimate interval | Median |
| --- | --- | ---: | ---: |
| 520 K | LSODE2 / ExprLegacy | `1.1155..1.1368 ms` | `1.1244 ms` |
| 520 K | LSODE2 / AtomView | `1.3407..1.4037 ms` | `1.3603 ms` |
| 520 K | BDF / ExprLegacy | `1.2784..1.2855 ms` | `1.2814 ms` |
| 520 K | BDF / AtomView | `1.2610..1.3739 ms` | `1.2903 ms` |
| 520 K | Radau / ExprLegacy | `1.1730..1.3178 ms` | `1.1976 ms` |
| 520 K | Radau / AtomView | `1.2357..1.3300 ms` | `1.2489 ms` |
| 900 K | LSODE2 / ExprLegacy | `1.4643..1.4659 ms` | `1.4648 ms` |
| 900 K | LSODE2 / AtomView | `1.4694..1.5638 ms` | `1.5364 ms` |
| 900 K | BDF / ExprLegacy | `1.2638..1.2653 ms` | `1.2642 ms` |
| 900 K | BDF / AtomView | `1.2870..1.3307 ms` | `1.3201 ms` |
| 900 K | Radau / ExprLegacy | `1.3221..1.3239 ms` | `1.3233 ms` |
| 900 K | Radau / AtomView | `1.2771..1.3395 ms` | `1.2841 ms` |

These are end-to-end telemetry-Off Criterion estimates from the earlier Dense-only matrix. The
expanded layout matrix has not yet been run. The short measurement window and host load produce
noticeable run-to-run movement; keep these as a dated baseline, not a performance gate.

**Stage findings:** In the initial detailed sample, KiThe `kinetic_processing` was the largest single
stage at roughly `0.8..1.1 ms` (`~60..83%` of workflow). These shares are independent medians per
stage divided by the median workflow duration; they are attribution estimates, not an additive
partition. The harness now chooses the diagnostic repetition closest to median workflow time for
paired stage/native breakdowns in future runs. This is cold full-solve behavior: each
call re-analyzes kinetics and rebuilds the symbolic rate expressions. Native integration ranged
from about `0.06 ms` for BDF/ExprLegacy to `0.62 ms` for LSODE2/AtomView in the diagnostic sample.
RST preparation was generally below `0.25 ms`, but was frontend-sensitive for this tiny problem.
The detailed stage percentages are per-stage summaries and must not be summed across independently
aggregated repeats. Native RST timings are nested within the KiThe integration stage and are printed
separately. The current report exposes detailed native rows for BDF and Radau; LSODE2 has only the
KiThe `native_integration` duration in this report.

**Prepared lifecycle measurements:**

| Solver / lifecycle | One-time prepare | First solve | Repeated solve median |
| --- | ---: | ---: | ---: |
| LSODE2 / same initial state | `0.0069 ms` | `0.2287 ms` | `0.0638 ms` |
| LSODE2 / prepared inlet-composition sweep | `0.0050 ms` | `0.1858 ms` | `0.0629 ms` |
| BDF / same initial state | `0.0470 ms` | `0.0259 ms` | `0.0197 ms` |
| BDF / prepared inlet-composition sweep | `0.0414 ms` | `0.0240 ms` | `0.0192 ms` |

The prepare values are StiffAPI problem preparation only; they exclude KiThe setup and are not
added to first/repeat solve time. The first LSODE2 run pays initial driver/callback work that the
retained prepared object avoids on subsequent runs. This is direct evidence that prepared reuse is
useful for repeated ranges, even on this small fixture.

**AOT versus Lambdify, LSODE2 only:** Run with `KITHE_RUN_IVP_AOT_BENCH=1` and `tcc` available.
Every pair used identical bounded inputs, Dense layout, and one of Auto/BDF/Adams controllers.
`prepare / integration / total` values below are medians over three Detailed runs; the same
per-route temporary artifact directory is reused across the three runs. Native `IvpTelemetry`
reports cold symbolic/AOT stages and warm callbacks separately. AOT/Lambdify endpoint parity and
the same temperature/fraction bounds passed for all 12 AOT combinations. Times are instrumented
diagnostics, not the telemetry-Off Criterion data above.

| Inlet / controller / frontend | Lambdify prepare / solve / total (ms) | AOT prepare / solve / total (ms) |
| --- | ---: | ---: |
| 520 K / Auto / ExprLegacy | `0.0499 / 0.1942 / 1.0200` | `0.0845 / 0.2154 / 1.1689` |
| 520 K / Auto / AtomView | `0.1627 / 0.2793 / 1.2050` | `0.1715 / 0.2253 / 1.3402` |
| 520 K / BDF / ExprLegacy | `0.0574 / 0.0698 / 0.6521` | `0.0619 / 0.0538 / 0.6413` |
| 520 K / BDF / AtomView | `0.1691 / 0.1249 / 0.8617` | `0.0986 / 0.0694 / 0.7555` |
| 520 K / Adams / ExprLegacy | `0.0463 / 0.0953 / 0.7341` | `0.0542 / 0.1069 / 0.8066` |
| 520 K / Adams / AtomView | `0.1164 / 0.1218 / 0.8267` | `0.1153 / 0.1338 / 0.8837` |
| 900 K / Auto / ExprLegacy | `0.0479 / 0.1910 / 0.7792` | `0.0537 / 0.1869 / 0.8190` |
| 900 K / Auto / AtomView | `0.1393 / 0.1992 / 0.8083` | `0.0855 / 0.1938 / 0.7472` |
| 900 K / BDF / ExprLegacy | `0.0505 / 0.0576 / 0.6235` | `0.0637 / 0.0558 / 0.7103` |
| 900 K / BDF / AtomView | `0.1523 / 0.1234 / 0.7718` | `0.1025 / 0.0672 / 0.6676` |
| 900 K / Adams / ExprLegacy | `0.0475 / 0.0981 / 0.7278` | `0.0546 / 0.1030 / 0.6313` |
| 900 K / Adams / AtomView | `0.1389 / 0.1535 / 0.7647` | `0.0992 / 0.1112 / 0.7898` |

The AOT result is mixed: some AtomView/BDF and warm-inlet routes improve, while several other
controller/frontend pairs are slower. It is not an across-the-board single-solve win. AOT remains
an LSODE2-only capability; standalone BDF and Radau reject that execution backend before solving.

**Conclusion:** Cold `kinetic_processing` remains the largest KiThe stage for ExprLegacy routes,
while LSODE2/AtomView can spend a substantial share inside native integration. Prepared reuse
materially reduces repeated solve cost. AOT benefits are route-dependent, and neither AOT nor
frontend selection should be generalized from this four-state control. Larger independently
validated workloads are still needed before prioritizing performance work. A previous exploratory run using
the older high-rate/default-domain fixture was discarded: it produced enormous finite states and
was not a valid performance/correctness baseline.

## Story: Checked-In Physical IVP Task Documents

**Name:** Physical IVP task files parse and produce bounded reactor solutions.

**Command:**

```powershell
cargo test --lib --no-default-features checked_in_physical_ivp_documents_load_and_solve -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can a user keep a readable physical IVP task as a normal file, load it
through the application document path, and solve it without constructing the mechanism in code?

**What is checked:** The H2 dissociation, two-step oxidation, and heat-releasing H2 documents are
loaded by `ReactorIvpApp`, parsed into the normal task configuration, and run through the production
IVP solve path. The test verifies finite states, temperature in `0..3000 K`, and species fractions
in `0..1`. The heat-releasing case must show measurable temperature and species evolution; its
elemental and mass-fraction balances must remain within `5%`, and its energy residual must close
within `1 J` and `0.01%`. An independent trapezoidal integration of the published profiles must
also close within `1 J`. Mechanisms use named species and indexed reaction fields rather than an
opaque JSON blob.

**Observed result:** Debug story passes for all three checked-in documents. The bounded smoke cases
retain their `520 K` inlet profiles; the coupled heat-release case publishes `520..522 K` and
asserts nonzero hydrogen conversion. The postprocessor energy residual is `-4.14e-2 J`
(`-3.85e-4%`); an independent trapezoidal integration gives `+5.90e-1 J`. Investigation found
that `simpsons()` skipped the second interval whenever it fell back on a trapezoid for nonuniform
mesh spacing. The corrected iterator and irregular-grid unit tests are the regression gate. No
solver or physics tolerance was relaxed, and no performance claim is made; file I/O is intentionally
outside performance measurement.

## Story: IVP GUI File Lifecycle

**Name:** Editing a physical IVP file preserves the task and remains runnable after reload.

**Command:**

```powershell
cargo test --lib --no-default-features physical_task_doc_gui_edit_save_reload_preview_and_solve -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can the GUI open a checked-in physical task, edit and save it as a new
file, reopen it, Preview and solve it, while rejecting an invalid replacement without partially
overwriting the active task?

**What is checked:** The test edits the heat-releasing H2 inlet temperature through the GUI numeric widget,
performs Save As/reopen and compares every typed config field, including the edited value, complete
reaction payload, all solver-route options, telemetry, logging, and postprocessing. It exercises
Preview and Run Calculation, then rejects a malformed indexed Arrhenius payload with
reaction context, verifies the active task remains unchanged, repairs the file, and completes a
second Save As/reopen/solve cycle. The checked-in source file must remain unchanged.

The reopened task is separately solved through the direct production API; final temperature and
species fractions and conservation diagnostics from the GUI worker must agree within `1e-8`.
Before the successful Save As, an intentionally unavailable destination is rejected; the active
source path, dirty state, and document fingerprint remain unchanged.

**Observed result:** Debug GUI lifecycle story passes. The widget edit survives Save As and reload;
the reopened task completes with finite solution data and a Kelvin temperature result. The
invalid-load attempt preserves the active document. Preview output is informational; there are no
file-I/O performance assertions.

## Story: Malformed IVP Mechanism Diagnostics

**Name:** Invalid indexed reactions fail with field context and do not replace the active task.

**Command:**

```powershell
cargo test --lib --no-default-features malformed_indexed_mechanisms_fail_with_context_without_replacing_active_task -- --nocapture --test-threads=1
```

**Question / hypothesis:** Do malformed readable mechanisms report the relevant reaction/field and
leave the currently loaded physical task intact?

**What is checked:** Four malformed documents cover an invalid reaction equation, a non-contiguous
reaction index, a product species absent from the declared mechanism, and a reaction/thermal-effect
count mismatch. Each failed load must expose diagnostic context while preserving the active
document fingerprint and path. The separate GUI lifecycle story covers a malformed Arrhenius
vector, repair, Save As/reopen, and successful solve.

**Observed result:** All four malformed variants are rejected transactionally in Debug. No partial
task replacement is observed; repair-and-run is covered by the physical GUI lifecycle story.

## Story: Failed IVP Rerun Result Invalidation

**Name:** A failed rerun cannot leave the previous solution presented as the result for edited inputs.

**Command:**

```powershell
cargo test --lib --no-default-features failed_ivp_rerun_discards_previous_solution_and_plot -- --nocapture --test-threads=1
```

**Question / hypothesis:** If a successful IVP result exists and a subsequent run uses invalid
solver settings, does the GUI discard the earlier published result rather than showing it as current?

**What is checked:** The story first completes a physical solve and confirms the solution snapshot
and plot exist. It then sets an invalid absolute tolerance and requests another run. Before the
worker result arrives, the prior publication must already be cleared. After failure, no solution
snapshot, conservation report, renderable solver state, or plot may remain.

**Observed result:** Debug story passes. The failed rerun reports the configuration error without
retaining the previous solution as current output.

## Story: IVP Solver Selection by Keyboard

**Name:** Keyboard focus and arrow navigation operate every solver settings selector.

**Command:**

```powershell
cargo test --lib --no-default-features solver_family_combo_opens_and_selects_with_keyboard -- --nocapture --test-threads=1
cargo test --lib --no-default-features all_solver_combos_support_keyboard_arrow_selection -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can a keyboard user focus the IVP solver selectors, change their values
with arrow keys, and use Enter/Space to open or close the option menu without relying on pointer
input?

**What is checked:** The story navigates the solver family, symbolic frontend, LSODE2 method and
matrix, execution backend, callback policy, telemetry, logging, and conditional AOT toolchain/profile
controls using AccessKit focus followed by keyboard input. It asserts the resulting typed GUI config.
The focused control uses egui's visible widget highlight; no separate hidden focus state is maintained.

**Observed result:** Debug story passes. Arrow navigation updates the typed selections across all
tested route states, including controls that appear only after selecting AOT. Existing mouse-driven
stories continue to cover opening menus and selecting individual entries.

## Story: TW-01 Analytical Multistep IVP Control

**Name:** The production IVP routes reproduce an independently constructed analytical profile for a
12-species, 20-reaction linear network.

**Commands:**

```powershell
cargo test --lib --no-default-features tw01_ -- --nocapture --test-threads=1
cargo test --lib --no-default-features expm_action_passes_independent_closed_form_checks -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the current KiThe mechanism parser, stoichiometry assembly, symbolic
RHS construction, scaling, and production StiffAPI route agree with a separate matrix-exponential
solution, rather than only agreeing across solver backends?

**What is checked:** TW-01 defines twelve pseudo-species with the same C2H4O2 composition and twenty
first-order reaction rows with distinct rates. The reference uses the steady convective species
balance `dY/dz=(L/m) diag(M) S^T r` with reactants negative and products positive. The input
substance ordering is reversed before production setup, so comparisons map concentration states back
by species name rather than relying on accidental index agreement. The test independently
constructs the linear generator from the reaction table and checks its two hand-derived rows,
column-sum invariant, positive profile, and thermal closed form. KiThe's parsed stoichiometry,
elemental matrix, Arrhenius constants, and zero heat release are checked against the fixed input.
The frozen CSV contains the analytical profile at
`z = 0.0, 0.1, ..., 1.0`; its ordinary test only reads it. LSODE2, standalone BDF, and Radau are each
compared to the analytical solution at every point on their own native adaptive mesh. Dimensional
publication is checked separately (`x=L*z`, `T=300+100*Teta`, published `q=T_scale*q_internal/L`),
without comparing different solver meshes by row index. A second LSODE2 run verifies that tighter
requested tolerances reduce profile error. The independent Padé-13 implementation has separate
zero, diagonal, nilpotent, 2x2 rotation, semigroup, and conserved-invariant checks.

**Observed result:** All five TW-01 tests and the standalone matrix-exponential test pass in Debug.
Maximum absolute species errors on the tested native meshes were `5.65e-10` (LSODE2), `1.60e-9`
(BDF), and `3.79e-13` (Radau); maximum total mass-fraction drift was below `1.89e-15`, and the
maximum C/H/O elemental mass-fraction drift was `1.11e-15`. Temperature errors were below
`2.95e-10 K`. Tightening LSODE2 tolerances from `rtol=1e-5, atol=1e-7` to
`rtol=1e-9, atol=1e-11` reduced maximum absolute state error from `5.85e-6` to `5.65e-10`.

**Scope:** This is a synthetic analytical numerical control, not an experimental or externally
validated kinetic mechanism. The solver integrates `z=x/L` and postprocessing reports physical
`x=L*z`. It is evidence for the documented steady convective model, not for other transport
assumptions or general kinetic accuracy.

## Story: TW-03 Non-Isothermal Arrhenius Kinetics With Zero Reaction Heat

**Name:** LSODE2, standalone BDF, and Radau reproduce three independently derived Arrhenius
profiles on the existing 12-species, 20-reaction control mechanism.

**Commands:**

```powershell
cargo test --lib --no-default-features tw03_ -- --nocapture --test-threads=1
cargo test --lib --no-default-features tw01_ -- --nocapture --test-threads=1
cargo test --lib --no-default-features tw02_ -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the production kinetic path use the absolute temperature and the
same activation-energy convention as an independent Arrhenius evaluator, while preserving the
current dimensionless mass-fraction RHS under an analytically prescribed temperature field?

**What is checked:** The audited KiThe contract is `k=A*T^n*exp(-E/(8.314*T))`, with `E` in
`J/mol`, `R=8.314 J/(mol K)`, and `T=300+100*Teta` in kelvin. The prefactors are independently
calibrated so every reaction retains its TW-01 rate at 400 K. The test evaluates all 20 production
symbolic constants at 400, 425, and 450 K, checks the analytic Arrhenius derivative, compares
production RHS rows to a separately assembled matrix, and confirms zero chemical heat.

For the distinct-Ea scenarios, the computed prefactors are:

| Reaction | `E` (J/mol) | `A` (1/s) |
| --- | ---: | ---: |
| R01 | 20000 | 2.863679222677e2 |
| R02 | 45000 | 8.279458321950e5 |
| R03 | 30000 | 4.137228398592e3 |
| R04 | 65000 | 4.926693562928e8 |
| R05 | 40000 | 1.338883051493e5 |
| R06 | 55000 | 1.979090962880e7 |
| R07 | 25000 | 1.103911254408e3 |
| R08 | 70000 | 2.631150866939e9 |
| R09 | 35000 | 3.349180734245e4 |
| R10 | 60000 | 9.585328940386e7 |
| R11 | 25000 | 1.379889068009e3 |
| R12 | 70000 | 1.661779494909e9 |
| R13 | 55000 | 6.850699486891e6 |
| R14 | 35000 | 3.535246330592e4 |
| R15 | 75000 | 1.058760312013e10 |
| R16 | 20000 | 2.250033674961e2 |
| R17 | 60000 | 7.188996705290e7 |
| R18 | 45000 | 1.354820452683e6 |
| R19 | 80000 | 1.820617322849e10 |
| R20 | 30000 | 1.241168519578e4 |

TW-03a uses distinct nonzero activation energies at constant 400 K and the TW-01 matrix exponential.
TW-03b uses a varying 400 to 450.025 K field and a common 40 kJ/mol activation energy; composite
Simpson quadrature reduces the system to a scalar reparameterization followed by a matrix
exponential. Its `N=256` / `N=1024` sigma difference is `1.684e-12`, and normalized commutator is
at roundoff. TW-03c uses the fixed distinct-Ea table, verifies a nonzero commutator, and builds an
independent fourth-order Magnus profile. The measured `N=640/1280/2560/5120` L2 refinements are
`2.495e-9 / 1.561e-10 / 9.758e-12`, with observed order `3.999`; the frozen reference uses 2560
steps per unit interval. Benchmark solver requests use `rtol=1e-11, atol=1e-13` so the dimensional
temperature criterion is met despite the deliberately large internal `qhat` state.

All three cases retain zero reaction heat, check positive finite fractions and total mass balance,
and compare each solver at its own native output points (no cross-mesh row pairing). The frozen
11-point profiles in [`tw_03_profiles.csv`](combustion_validation/fixtures/tw_03_profiles.csv) are
read-only in ordinary tests; method and convergence provenance is recorded in
[`tw_03_reference_metadata.toml`](combustion_validation/fixtures/tw_03_reference_metadata.toml).
The test converts internal initial `qhat=5000` back to dimensional `q_physical` expected by the
task API. This preserves the requested temperature field while keeping TW-01's actual
`rho,m,Lambda,L` and chemical residence scale. The alternative parameter values and
`-rho*L^2*...` species matrix in the supplied draft do not match the current KiThe production model
and were not used.

**Observed result:** All seven TW-03 tests pass in Debug, including all nine solver/scenario
combinations. In the full serial campaign, maximum species profile errors were `9.07e-12`
(LSODE2), `3.70e-11` (BDF), and `6.15e-12` (Radau); maximum temperature discrepancy was
`2.47e-9 K`; total mass-fraction drift was below `5.6e-15`. TW-03c's minimum mass fraction
remained positive at approximately
`6.50e-4`.

**Scope:** This is a synthetic numerical validation of Arrhenius evaluation, non-autonomous RHS
assembly, and solver routes. It does not validate an experimental mechanism or independently
resolve the separate physical sign convention of the thermal balance.

## Story: TW-02 Coupled Thermochemical Analytical IVP Control

**Name:** All supported stiff solver families reproduce the analytical coupled 14-state TW-02
profile while preserving the TW-01 chemistry trajectory.

**Commands:**

```powershell
cargo test --lib --no-default-features tw02_ -- --nocapture --test-threads=1
cargo test --lib --no-default-features tw01_ -- --nocapture --test-threads=1
```

**Question / hypothesis:** When nonzero reaction heats are derived consistently from species
enthalpies and coupled into the heat equation, do the existing KiThe reaction-rate and thermal
source assembly paths agree with a separately assembled constant-coefficient 14x14 system and
its matrix-exponential solution?

**What is checked:** TW-02 reuses the exact TW-01 species, 20 reactions, kinetic constants, initial
composition, and physical parameter values. It assigns synthetic molar enthalpies
`[0,-40,-25,-65,-15,-80,-35,-100,-50,-90,-60,-115] J/mol` and derives each reaction heat as
`Q_j = h_reactant - h_product`. The test asserts all twenty resulting heats, two independent
stoichiometric cycle sums (Hess-law consistency), and invariance to a common enthalpy offset.
At five points, it checks each production reaction rate and its individual `Q_j * rate_j`
contribution against the independently assembled rate map, the total heat source, and the heat RHS
against a separately built 14x14 matrix.

The matrix exponential supplies the reference at every point of each solver's own adaptive mesh.
The frozen CSV stores an 11-point analytical profile: the common species profile and `Teta,q` for
`gamma = +1,+10,-1`. Assertions cover componentwise species/thermal profile error, dimensional
temperature, finite positive species values, total mass fraction, C/H/O elemental mass fractions,
the local heat RHS computed from independently reconstructed reaction heat, and local heat-source
agreement. Metamorphic checks require `gamma=0` to recover TW-01, chemistry to remain invariant
under gamma, and thermal changes to scale linearly with gamma. The reference metadata and fixture
are read-only during normal tests.

**Observed result:** All five TW-02 tests pass in Debug. LSODE2, standalone BDF, and Radau pass for
each gamma value. Across the nine solver/gamma combinations, maximum absolute mass-fraction error
is `1.59e-9`, maximum `Teta` error is `1.74e-11`, maximum q error is `1.32e-9`, maximum
dimensional temperature error is `1.74e-9 K`, and maximum local energy-RHS discrepancy is
`1.51e-15`. The largest C/H/O mass-fraction balance drift is `1.22e-15`; all temperature profiles
remain within `400.000..400.001 K`.

**Model reference:** The physical assumptions, units, signed stoichiometry, dimensional balances,
and nondimensionalization are specified in
[`MATHEMATICAL_MODEL.md`](MATHEMATICAL_MODEL.md). TW-02 remains a synthetic analytical control;
it verifies the coupled RHS and solver routes against a separate matrix exponential, not against
experimental combustion data.

## Story: TW-04 Fully Coupled Nonlinear Thermochemical IVP

**Name:** All supported StiffAPI families reproduce an independently integrated 14-state
Arrhenius/thermochemical IVP across calibrated weak, moderate, and strong heat-feedback regimes.

**Commands:**

```powershell
cargo test --lib --no-default-features tw04_ -- --nocapture --test-threads=1
Push-Location src/SimpleReactorsIVP/combustion_validation/external_reference
py -3 -m pytest tw_04/test_reference.py -q
py -3 -m tw_04.calibrate
Pop-Location
```

`calibrate` is an explicit research command. It writes the probe and three-level
SciPy convergence matrix to `external_reference/tw_04/calibration.json`; ordinary
Rust tests only read the frozen profile and provenance JSON.

**Question / hypothesis:** Does the production pipeline correctly couple all 20
first-order Arrhenius reactions to the temperature and heat-flux equations when
reaction heats follow Hess's law? Do LSODE2, standalone BDF, and Radau match an
independently assembled SciPy reference while conserving mass, C/H/O inventory,
and the exact algebraic energy invariant?

**Fixture and route:** Reuses the exact 12 equal-formula-mass TW-01 pseudo-species,
20 reactions, initial composition, TW-02 component enthalpies, and TW-03c
activation energies and reference-rate calibration. The physical configuration is
the established `MATHEMATICAL_MODEL.md` setup (`rho=1000 kg/m^3`, `m=1 kg/(m^2 s)`,
`Lambda=10000 W/(m K)`, `L=0.01 m`, `Cp=1000 J/(kg K)`, `T=300+100*Teta`).
The task API receives dimensional initial `q=5e7 W/m^2`, which the production
builder maps to internal `qhat=5000`; dimensional heat flux is derived as
`qhat*T_scale/L`.

The independent Python implementation builds its own stoichiometric matrix,
concentrations, Arrhenius constants, reaction rates, mass sources, heat release,
and RHS. DOP853 is the frozen route; Radau independently cross-checks it. Both
are run at three tolerances for calibration, and the frozen metadata captures the
full mechanism, enthalpies, derived reaction heats, parameters, grid, versions,
solver options, method deltas, mass error, and energy-invariant error. Each frozen
profile uses a common 0.001 grid with extra points around maximum `|Qvol|`; its
RHS derivatives provide cubic-Hermite interpolation onto each Rust backend's own
adaptive output points.

**Calibrated scenarios:**

| Scenario | `gamma` | max `|T_gamma-T_0|` | max `|Y_gamma-Y_0|` | departure from weak linear extrapolation | min `Y` |
| --- | ---: | ---: | ---: | ---: | ---: |
| TW-04a, zero-feedback control | 0 | 0 K | 0 | 0 K | `6.496e-4` |
| TW-04b, weak | `1e3` | `1.874e-2 K` | `8.788e-5` | reference | `6.498e-4` |
| TW-04c, moderate | `1e5` | `2.063 K` | `9.679e-3` | `3.149e-1 K` vs `gamma=1` tangent | `6.827e-4` |
| TW-04d, strong nonlinear | `5e5` | `2.190e1 K` | `9.196e-2` | `1.848e1 K` vs `gamma=1` tangent | `2.001e-3` |

The `gamma=0` endpoint temperature is `450.025008 K`. The strong case stays within
`400..428.669 K`; all selected reference trajectories have strictly positive mass
fractions. It is not labelled “strong” merely because gamma is numerically large:
the profile departs by `18.48 K` from the small-gamma linear extrapolation.

**Acceptance checks:** At control states, compare all 20 production Arrhenius
constants/rates, 12 species RHS components, `Teta'`, `qhat'`, and total `Qvol`
componentwise with the independent equations. Across each complete trajectory,
check all species, `Teta`, `qhat`, dimensional temperature, and dimensional heat
flux; finite positive fractions; total mass and C/H/O mass inventories; and the
raw/normalized energy invariant. No clipping or renormalization is applied.
Metamorphic assertions cover `gamma=0` against TW-03c, `Ea=0` constant-rate
behavior, linear scaling of Hess-derived heats with gamma, invariance to a common
enthalpy offset, the weak limit, and a measurable strong nonlinear response.

Comparison budgets are `1e-9 + 1e-7*|Y_i|` for each mass fraction,
`1e-9 + 1e-7*|Teta|`, `1e-5 + 1e-7*|qhat|`, and
`2e-6 + 1e-7*|T_K| K` for dimensional temperature. The derived physical heat
flux uses `2e-2 + 1e-9*|q_physical| W/m^2`. Conservation uses `2e-9` absolute
mass/element error; energy uses normalized error `1e-8` and raw error `1e-3` in
internal heat-flux units. The frozen `rtol=1e-12`, `max_step=2e-4` DOP853/Radau
cross-check is at least three orders tighter than the species acceptance budget
on the strong profile.

**Observed result:** Debug Rust run: 4 TW-04 tests passed, including the local
production-RHS gate and all 12 scenario/backend combinations. Across the solver
matrix the largest mass-fraction difference from the frozen reference was
`4.92e-11`, largest dimensional temperature difference `1.16e-8 K`, and largest
internal `qhat` difference `1.51e-6`. Minimum published fraction was
`6.496e-4`; mass and C/H/O conservation stayed within `2e-9`; the independently
reconstructed normalized energy invariant stayed below `1e-8`. Python tests:
6 passed. The calibrated DOP853/Radau full-profile differences were at most
`1.21e-9 K` and `1.04e-9` in mass fraction for the strong scenario; tighter
`rtol=1e-12`, `max_step=2e-4` frozen-reference cross-checks were smaller still.

**Conclusion:** The current nonlinear production RHS and all three solver-family
routes pass this synthetic, independently coded, thermochemically coupled case.
This is mathematical/code-path validation for the declared simplified model, not
experimental combustion validation or an independent physical audit of the sign
convention for positive heat flux. Debug and Python correctness runs passed;
the attempted Release compile was stopped after prolonged compilation near
12 GB RAM, so Release status remains unverified.
