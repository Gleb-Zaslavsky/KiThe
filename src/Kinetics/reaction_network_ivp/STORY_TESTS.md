# Reaction-Network IVP Story Tests

These records describe executable correctness tests. They are not performance
measurements. The named tests live in `tests.rs` and use the public network
facade with RST 0.5.0 from the current lockfile.

## Story 1: Four Native Families Solve An Analytical Network

**Name:** `A -> 2B` with prescribed constant temperature.

**Question / hypothesis:** Do BDF, Radau, Backward Euler, and LSODE2 receive
the same symbolic reaction network, preserve species ordering, and publish a
physically admissible solution matching `A(t)=exp(-0.5t)`,
`B(t)=2(1-exp(-0.5t))` over `0 <= t <= 2 s`?

**Commands:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features four_stiffapi_families_solve_the_same_analytical_network -- --nocapture
cargo test --lib --no-default-features bdf_facade_solves_analytical_a_to_two_b_without_telemetry_by_default -- --nocapture
```

**Fixture and acceptance:** `A(0)=1`, `B(0)=0 mol/m^3`, `k=0.5 s^-1`,
`T=700 K`. Each family has `max_step=0.1 s`; default BDF uses
`rtol=1e-8`, `atol=1e-10`. The cross-family test accepts final `A` error
`<2e-2`, final `B` error `<4e-2`, and invariant `|2A+B-2|<2e-2`.
The tighter BDF test checks all returned samples against the analytical
curves (`A <2e-5`, `B <4e-5`) and invariant (`<2e-5`). The facade separately
requires finite non-negative concentrations, increasing time, and terminal
coverage.

**Result:** Both tests passed in Debug with RST 0.5.0. Default telemetry is
Off and the published BDF solution has no telemetry snapshot.

**Conclusion / limitation:** This proves the route and elementary physical
contract on one small network. It does not establish equally tight accuracy
for every family, large/stiff networks, or solid-state stages.

## Story 2: Prepared Rerun And Transactional Publication

**Name:** Same prepared network, changed initial concentrations.

**Question / hypothesis:** Can each StiffAPI family reuse its prepared
residual/Jacobian callbacks for a fresh solve on the same time interval,
without treating the previous endpoint as a continuation seed? Does an
invalid subsequent request leave the last accepted solution intact?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features prepared_session_reuses_all_four_families_and_retains_last_good_result -- --nocapture
```

**Fixture and acceptance:** First run: `A(0)=1`, `B(0)=0`. Second run:
`A(0)=2`, `B(0)=0.5 mol/m^3`. At `t=2 s`, the second reference is
`A=2exp(-1)`, `B=0.5+2(2-2exp(-1))`. For each family, final absolute
errors must be `<2e-2` for `A` and `<4e-2` for `B`. Counters telemetry is
requested and must yield a native snapshot. A third request with `A(0)=-1`
must return a typed input error and leave the second accepted endpoint
unchanged.

**Result:** Passed in Debug for BDF, Radau, Backward Euler, and LSODE2 with
RST 0.5.0. The session uses `try_solve_with_initial_state()` on a previously
prepared native solver; equations, temperature law, and interval are unchanged.

**Conclusion / limitation:** Same-formulation initial-state reuse is working.
Numeric parameter rebinding, changed intervals, and retention after an actual
native integration failure remain separate tests/work items.

## Story 3: Cancellation Never Publishes A Partial Result

**Name:** Cancelled prepared rerun with a previously accepted solution.

**Question / hypothesis:** For each of the four StiffAPI families, does an
already-cancelled token produce a typed native error without replacing the
first accepted trajectory?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features cancellation_does_not_replace_an_accepted_solution -- --nocapture
```

**Fixture and acceptance:** Prepare the same `A -> 2B` network as Stories 1
and 2 with Counters telemetry. Accept the first run from `A(0)=1`, `B(0)=0`.
Request a second run from `A(0)=2`, `B(0)=0.5` with an already-cancelled
`CancellationToken`. The rerun must return
`NetworkIvpError::CancellationAwareIntegration` with a native error source;
the previous published endpoint must remain byte-for-byte unchanged.

**Result:** Passed in Debug for BDF, Radau, Backward Euler, and LSODE2 with
RST 0.5.0.

**Conclusion / limitation:** Cooperative cancellation is wired through the
common facade and publication is transactional for this early-cancel case.
The test does not characterize mid-integration cancellation latency or the
availability of a partial trajectory/telemetry snapshot after such a stop.

## Story 4: Numeric Temperature Rebind Without Symbolic Repreparation

**Name:** Constant-temperature and linear-heating parameter loops.

**Question / hypothesis:** Can one prepared StiffAPI solver reuse the same
symbolic residual/Jacobian after changing the isothermal temperature or the
linear-heating rate, with optional simultaneous change of initial
concentrations? Are the published temperature and concentrations consistent
with the *new* numeric binding rather than the original one?

**Commands:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features prepared_isothermal_arrhenius_rebinds_temperature_and_initial_state -- --nocapture
cargo test --lib --no-default-features prepared_linear_heating_rate_rebinds_without_rebuilding_equations -- --nocapture
cargo test --lib --no-default-features parameterized_temperature_rebind_has_atom_view_expr_legacy_parity -- --nocapture
cargo test --lib --no-default-features invalid_temperature_rebind_is_typed_and_does_not_replace_last_solution -- --nocapture
```

**Fixture and acceptance:** Both examples solve `A -> B` over `0..2 s`.
The Arrhenius case uses `A=1 s^-1`, `E_a=R*500 J/mol`, first `T=500 K`,
then `T=1000 K` and `C(0)=[2,0.5]`. Its analytical reference is
`C_A(2)=C_A(0)*exp(-2*exp(-E_a/(R*T)))`. The linear case uses complete
custom rate `r=C_A*T/500`, first `T=500+100*t`, then `T=500+200*t` and
`C(0)=[2,0.5]`; its independent reference is
`C_A(t)=C_A(0)*exp(-T0*t/500-beta*t^2/1000)`.

The temperature parameter prefix is `[p_temperature]` or
`[p_temperature,p_heating_rate]`; parameterized rate coefficients follow in
stable reaction order (`p_rate_{reaction_index}_k` for Constant and
`p_rate_{reaction_index}_a/e` for Arrhenius). BDF, Radau, Backward Euler, and
LSODE2 must solve both conditions in a single `PreparedNetworkIvp` session.
Continuous
reference errors must be `<1e-2` for `A` (Backward Euler Arrhenius rebind:
`<2e-2`) and `<2e-2` for `B`; the Backward Euler result must additionally
match its independent per-step implicit-Euler recurrence to `<1e-8`.
Parameterized AtomView/ExprLegacy BDF endpoints must agree to `<1e-8`.
Mismatched law shape and negative Kelvin must return typed errors without
changing the accepted solution or binding.

**Result:** All four tests passed in Debug with RST 0.5.0. In the
isothermal example, Backward Euler's second `A` endpoint was
`0.605389157555665`; its discrete recurrence gives
`0.6053891575556639`. The continuous reference is `0.5945715963705381`;
the difference is the expected first-order time-discretization error at
`max_step=0.05 s`, not a failed parameter bind.

**Conclusion / limitation:** Numeric T/linear-beta rebinding is functional
without rebuilding the KiThe Expr graph or creating a second native solver.
This is a correctness story, not a measured speedup. Kinetic rate parameter
rebinding and arbitrary symbolic-temperature shape changes remain outside
this contract.

## Story 5: Explicit Forward And Reverse Stages

**Name:** Reversible first-order pair represented by two user reactions.

**Question / hypothesis:** Are forward and reverse laws treated as independent
stages, preserving the declared species ordering and the analytical reversible
solution without inferring a reverse rate?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features explicit_reverse_reactions_preserve_species_order_and_approach_equilibrium -- --nocapture
```

**Fixture and acceptance:** `A -> B` with `k_f=2 s^-1` and `B -> A` with
`k_r=1 s^-1`, `A(0)=1`, `B(0)=0`, at constant `600 K`. At `t=5 s`,
`A=1/3+(2/3)exp(-3t)` and `B=(2/3)(1-exp(-3t))`. Both endpoint errors
must be `<2e-5`; every sampled state must preserve `A+B=1` within `2e-6`.

**Result:** Passed in Debug with RST 0.5.0. Stoichiometry is two ordered
rows, one for each explicitly supplied reaction.

**Conclusion / limitation:** Reverse kinetics remain user-specified. No
equilibrium constant or thermodynamic reverse law is inferred.

## Story 6: Public-Facade Example Compiles

**Name:** Arrhenius plus custom-rate chain with linear prescribed heating.

**Question / hypothesis:** Can an application define abstract species,
heterogeneous rate-law variants, and a non-isothermal input using only the
public facade and compile it as an ordinary crate example?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo check --example reaction_network_ivp --no-default-features
```

**Result:** Passed in Debug. The example constructs `A -> 2B` with an
Arrhenius coefficient and `B -> C` with a complete custom expression, then
prints the named final concentrations and prescribed temperature.

**Conclusion / limitation:** This is API/compilation evidence, not a physical
acceptance or performance gate.

## Story 7: RHS Domain Error Context

**Name:** Non-finite rate identifies its output component.

**Question / hypothesis:** Does the independent RHS diagnostic return a typed
domain error naming the affected species and time, and does a malformed
middle reaction reject the entire network construction with that reaction's
identity rather than yielding a partial mechanism?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features rhs_domain_error_identifies_the_nonfinite_state_component -- --nocapture
cargo test --lib --no-default-features validation_keeps_reaction_context_and_rejects_unknown_symbols -- --nocapture
```

**Result:** Passed. A singular custom expression reports time `0.25 s`,
component index `1`, and species ID `B`. A network with valid reactions before
and after an invalid middle custom expression returns the typed
`RateSymbol` error for reaction `malformed-middle`; no `NetworkProblem` is
published.

**Limitation:** During integrated RST callbacks, native evaluation errors are
still preserved at the RST boundary; KiThe does not yet enrich them with a
per-reaction expression path.

## Story 8: Native And Calculator Telemetry Mode Matrix

**Name:** Off/Counters/Timings/Detailed across all StiffAPI families and
KiThe-owned preparation/integration stages.

**Question / hypothesis:** Does the facade forward each independent telemetry
mode to BDF, Radau, Backward Euler, and LSODE2, while preserving the distinct
mode semantics in calculator-owned stage diagnostics?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features telemetry_modes_are_forwarded_independently_to_all_stiffapi_families -- --nocapture
cargo test --lib --no-default-features prepared_session_reuses_all_four_families_and_retains_last_good_result -- --nocapture
```

**Fixture and acceptance:** Identical `A -> 2B` network, same tolerances,
temperature, and time interval for all routes. The test crosses all four
families with `Off`, `Counters`, `Timings`, and `Detailed`. Off must publish
no snapshots; every enabled mode must publish a native snapshot and four
KiThe stages: equation assembly, backend preparation, integration, and output
postprocessing. `Counters` contains stage calls without durations; `Timings`
contains durations without fabricated calls; `Detailed` contains both. A
prepared two-run story verifies assembly/preparation are counted once while
integration/postprocessing are counted twice.

**Result:** All 16 combinations passed in Debug with RST 0.5.0. Native
snapshots were absent only for Off. Calculator stages were absent in Off;
Counters contained calls but no durations; Timings contained durations but no
calls; Detailed contained both. The prepared two-run test reported one
assembly/preparation and two integration/postprocessing calls.

**Conclusion / limitation:** Both native mode forwarding and KiThe-owned
stage collection are covered. This does not assert native counter/timing
values or equal semantics between families; those remain native RST contracts.
Output postprocessing is currently one macro-stage, not a split between
trajectory extraction, row materialization, and physical validation. Failed
prepared solves retain KiThe stage totals on the runner, but not every RST
family exposes a failed native snapshot.

## Story 9: Coupled Reaction Network And Stoichiometric Invariant

**Name:** `A -> 2B; B + C -> D` across all StiffAPI families.

**Question / hypothesis:** Does the assembled network integrate both reaction
rates with the declared stoichiometry, preserve the independent invariant
`C+D=1 mol/m^3`, and produce comparable endpoints across BDF, Radau, Backward
Euler, and LSODE2?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features coupled_two_reaction_network_preserves_invariant_across_stiffapi_families -- --nocapture
```

**Fixture and acceptance:** `A(0)=1`, `C(0)=1`, `B(0)=D(0)=0 mol/m^3`,
`k1=0.5`, `k2=0.1`, constant `T=700 K`, interval `0..2 s`, and
`max_step=0.005 s`. Every sample must be finite, nonnegative, and preserve
`C+D` within `2e-6`; all endpoint components must agree with BDF within
`2e-3`.

**Result:** Passed in Debug with RST 0.5.0 for BDF, Radau, Backward Euler,
and LSODE2. All sampled states were finite/nonnegative, `C+D` stayed within
the `2e-6` gate, and all cross-family endpoint differences stayed below
`2e-3`.

**Conclusion / limitation:** This is an internal invariant for the declared
abstract stoichiometry; it does not establish elemental conservation without
composition metadata or a broad large/stiff-network guarantee.

## Story 10: Independently Checked State Jacobian

**Name:** Mixed symbolic rates, derivative matrix versus central differences.

**Question / hypothesis:** Does differentiating the assembled full-state RHS
with respect to each concentration at fixed time agree with central finite
differences of the independently probed numeric RHS, including Arrhenius,
empirical-order, and complete custom-rate stages under linear prescribed
temperature?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features assembled_state_jacobian_matches_independent_rhs_finite_differences -- --nocapture
```

**Fixture and acceptance:** Three species and three reactions are evaluated
at `t=0.37 s` and positive concentrations `[1.2, 0.9, 0.55] mol/m^3`.
Every symbolic Jacobian element must agree with its centered finite-difference
estimate to `2e-7 * max(|symbolic|, |finite_difference|, 1)`.

**Result:** Passed in Debug with RST 0.5.0.

**Conclusion / limitation:** This validates KiThe's symbolic RHS derivative
semantics on the tested expression families. RST still owns callback
compilation and integration Jacobian execution; this test is not an AOT or
solver-level derivative benchmark.

## Story 11: Inert Species And Empty Reaction Set

**Name:** Valid zero-reaction network.

**Question / hypothesis:** Does the public API allow an inert network and keep
each concentration exactly constant through the requested time interval?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features empty_reaction_set_integrates_inert_species_without_changing_state -- --nocapture
```

**Fixture and acceptance:** Two species start at `[1.25, 0] mol/m^3`, with no
reactions, constant `400 K`, and interval `0..3 s`. Species order must be
preserved, at least both interval endpoints must be published, and every row
must equal the initial state exactly.

**Result:** Passed in Debug with RST 0.5.0. The default BDF route published
the initial state unchanged at every returned sample.

**Conclusion / limitation:** This confirms inert-state semantics on the
default BDF route only; it is not a full-family zero-RHS guarantee.

## Story 12: Release Performance Baseline Definition

**Name:** Linear-chain cold setup versus prepared rerun.

**Question / hypothesis:** How do cold prepare+solve and repeated prepared
solves scale with network size (3, 15, 40 states) across BDF, Radau, Backward
Euler, and LSODE2, and which KiThe preparation/integration stages account for
the measured prepared workflow?

**Command:**

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactionNetworkIvpCampaign.ps1 -Suite performance -BuildJobs 1
```

**Fixture and measurement:** Each fixture is a first-order linear chain at
`700 K`, `0..2 s`, with `max_step=0.02 s`. Criterion crosses
ExprLegacy/AtomViewNative, 3/15/40 states, and the four solver families.
Matrix layout is Dense and callbacks are sequential. Criterion samples use
telemetry Off.
The separate Tabled diagnostic pass reports prepared-workflow wall time,
equation assembly, native preparation, integration, output postprocessing,
and signed unattributed time. Rows failing preflight are shown as failures and
not benchmarked. The benchmark is a baseline, not a performance gate.

**Result:** Release campaign passed with RST 0.5.0 on 2026-11-21. Run ID:
`reaction-network-ivp-performance-release-20261121-175402-396Z`. Criterion
medians, confidence intervals, and detailed KiThe stage attribution are in
`test_reports/reaction-network-ivp-performance-release-20261121-175402-396Z/reports/summary.md`;
the sibling `technical` directory preserves raw output.

Representative 40-state cold / prepared Criterion medians (ms):

| Family | ExprLegacy cold | ExprLegacy prepared | AtomViewNative cold | AtomViewNative prepared |
|---|---:|---:|---:|---:|
| BDF | 0.9284 | 0.5172 | 1.0171 | 0.4814 |
| Radau | 4.5618 | 3.7220 | 4.3243 | 3.6210 |
| Backward Euler | 2.0704 | 1.5329 | 2.0167 | 1.4189 |
| LSODE2 | 2.8126 | 2.1926 | 3.0153 | 2.1044 |

For 40-state Radau, the diagnostic pass attributes approximately
`1.37-3.85 ms` to integration depending on frontend/layout, with preparation
around `0.29-0.73 ms`; exact combinations are in the report. Criterion and
one-pass stage attribution are distinct measurements and must not be summed.

**Conclusion / limits:** This is a machine- and fixture-specific baseline,
not a performance gate or universal ranking. The stage-attribution pass had
one `ExprLegacy/BDF/Auto` first-call outlier (`2154.862 ms` preparation),
while its Criterion median was `0.9237 ms`; the cause is not isolated, so the
diagnostic outlier must not be treated as steady-state solver cost. Forced
Parallel was substantially slower for this 40-state workload; see Story 17.
No speedup is claimed.

## Story 13: Frontend Parity Across Scales And Solver Families

**Name:** ExprLegacy versus AtomViewNative on a first-order analytic network.

**Question / hypothesis:** Do both symbolic frontends publish the same
first-order solution for BDF, Radau, Backward Euler, and LSODE2 at small,
unit, and large concentration scales, without changing tolerances to force
agreement?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features frontend_and_family_matrix_preserves_first_order_solution_across_scales -- --nocapture
```

**Fixture and acceptance:** `A -> B`, `k=0.5 s^-1`, `700 K`, interval
`0..2 s`, and `C_A(0) = 1e-3, 1, 1e3 mol/m^3`. Each family/frontend result
must match normalized analytic endpoints `A/A0=exp(-1)` and
`B/A0=1-exp(-1)` within `2e-2`; corresponding frontend endpoints must agree
within `2e-5 * A0`.

**Result:** Passed in Debug with RST 0.5.0 for all 24 family/frontend/scale
combinations.

**Conclusion / limitation:** This demonstrates tested scale invariance over
six orders of magnitude for a linear first-order system. It does not establish
general concentration-scale robustness for nonlinear or stiff networks.

## Story 14: Stiff Consecutive First-Order Network

**Name:** Fast `A -> B` followed by slow `B -> C`.

**Question / hypothesis:** Can all current StiffAPI families solve a stiff
consecutive reaction network and remain consistent with its independent
closed-form trajectory without publishing negative concentrations?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features stiff_consecutive_first_order_network_matches_closed_form_for_all_families -- --nocapture
```

**Fixture and acceptance:** `A(0)=1`, `B(0)=C(0)=0 mol/m^3`, `k_fast=100 s^-1`,
`k_slow=0.2 s^-1`, constant `T=700 K`, interval `0..0.1 s`, and
`max_step=0.005 s`. The rate ratio is 500:1. With `t=0.1 s`, the reference is
`A=exp(-k_fast*t)`,
`B=k_fast/(k_fast-k_slow)*(exp(-k_slow*t)-exp(-k_fast*t))`, and
`C=1-A-B`. All endpoint absolute errors must be below `2e-3`; every returned
sample must be finite, non-negative, and preserve `A+B+C=1` within `2e-5`.

**Result:** Passed in Debug for BDF, Radau, Backward Euler, and LSODE2 with
RST 0.5.0.

**Conclusion / limitation:** This adds a genuinely stiff *linear* network gate.
It does not establish robustness for nonlinear stiff kinetics, discontinuous
temperature inputs, or arbitrary concentration scales.

## Story 15: Prepared Rate-Coefficient Rebinding

**Name:** Full reaction-ID binding set on a prepared network.

**Question / hypothesis:** Can a prepared symbolic network accept a complete
replacement set of kinetic coefficients keyed by reaction ID, independent of
input ordering, and solve the resulting system with all StiffAPI families?
Do incomplete, duplicate, or type-incompatible sets fail before integration
and preserve the last accepted result?

**Commands:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features prepared_rate_coefficients_rebind_by_reaction_id_across_families -- --nocapture
cargo test --lib --no-default-features invalid_rate_rebind_is_typed_and_preserves_last_solution -- --nocapture
```

**Fixture and acceptance:** `A -> B -> C`, initial state `[1,0,0]`, constant
`T=700 K`, interval `0..1 s`, and rebound rates `k1=1.2 s^-1`,
`k2=0.4 s^-1`. Bindings are intentionally supplied in reverse reaction order.
The endpoint oracle is `A=exp(-k1)`,
`B=k1/(k1-k2)*(exp(-k2)-exp(-k1))`, and `C=1-A-B`; each family must be within
`1e-2`. Validation stories require typed errors for missing entries,
duplicates, and a Constant/Arrhenius kind mismatch, while retaining the
previous accepted endpoint.

**Result:** Both tests passed in Debug with RST 0.5.0. BDF, Radau, Backward
Euler, and LSODE2 reproduced the analytical endpoint gate.

**Conclusion / limitation:** Rebinding does not rebuild symbolic callbacks and
is explicit rather than stateful: each call must provide every rebindable
reaction coefficient. It does not yet provide a rate parameter-loop report or
benchmark the savings versus cold preparation.

## Story 16: Localized Output Validation Errors

**Name:** Invalid output samples retain component and coordinate context.

**Question / hypothesis:** When a native trajectory contains a negative or
non-finite concentration, or a non-increasing time sample, does the KiThe
publication gate return a typed error with the relevant species/time/index
rather than a generic message?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features invalid_published_samples_report_component_and_time -- --nocapture
```

**Fixture and acceptance:** A synthetic two-row candidate for `A -> 2B` has
`A=-1e-12 mol/m^3` at `t=0.1 s`; it must return
`InvalidConcentration { species_id: "A", time: 0.1, value: -1e-12 }` without
publishing a partial result. A second candidate with duplicate time samples
must return `InvalidTimeSample` with row index 1 and the previous time.

**Result:** Passed in Debug with RST 0.5.0.

**Conclusion / limitation:** KiThe-owned output rejection now carries concrete
sample context. This does not add missing reaction/expression context to
failures raised inside native RST callbacks.

## Story 17: Release Matrix Layout And Callback Policy

**Name:** Supported 40-state StiffAPI layouts and callback policies.

**Question / hypothesis:** For the largest baseline chain, which matrix
layouts and callback policies are accepted by each family, and does forced
parallel evaluation help this callback workload? Unsupported combinations
must be visible rather than silently substituted.

**Command:**

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactionNetworkIvpCampaign.ps1 -Suite performance -BuildJobs 1
```

**Fixture and measurement:** 40-state/39-reaction first-order chain. Both
ExprLegacy and AtomViewNative are measured. BDF, Radau, and LSODE2 compare
Sequential, forced Parallel (`min_work=1`), and Auto (`min_work=1`) callbacks;
Backward Euler has no callback-policy control and is reported with `None`.
Dense, Sparse, and BandedAuto are attempted only where supported. Criterion
medians are separate from the Tabled one-pass stage attribution.

**Representative callback medians (ms, ExprLegacy):**

| Family / layout | Sequential | Parallel | Auto |
|---|---:|---:|---:|
| BDF / Dense | 0.9596 | 5.0635 | 0.9237 |
| Radau / Dense | 4.5143 | 22.8327 | 4.4973 |
| Radau / BandedAuto | 1.5165 | 12.3592 | 1.6088 |
| LSODE2 / Dense | 2.9649 | 13.3544 | 2.8907 |
| LSODE2 / BandedAuto | 1.8252 | 12.9507 | 1.8497 |

Sparse and all supported AtomViewNative combinations are detailed in the
campaign report. BDF and Backward Euler reject Sparse/BandedAuto; the report
marks these `UNSUPPORTED`, not as failures and not as fallback runs.

**Result:** The release campaign passed. Forced Parallel was substantially
slower than Sequential/Auto in every supported family/frontend/layout
combination in this fixture (roughly 3-10x). Auto tracked Sequential closely
in the warmed Criterion measurements. Parallel is therefore not recommended
by this fixture; larger or more expensive callbacks need their own crossover
measurement.

**Conclusion / limitation:** This is capability and performance evidence for
one linear chain at 40 states, not a general threshold for nonlinear
networks, other machines, or Rayon pool states. The diagnostic pass also
recorded a single BDF/Auto preparation outlier; see Story 12 and the raw
report before interpreting first-call latency.

## Story 18: Native Logging Level Forwarding

**Name:** Opt-in common solver logging without enabling telemetry.

**Question / hypothesis:** Does the public calculator default to silent native
logging and forward an explicitly selected level to each supported StiffAPI
family independently of telemetry?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features native_logging_level_is_forwarded_to_all_stiffapi_families -- --nocapture
```

**Fixture and acceptance:** Prepare the same `A -> 2B` network with each of
BDF, Radau, Backward Euler, and LSODE2, selecting `SolverLogLevel::Debug` and
leaving telemetry at `Off`. Inspect each prepared native solver's logging
level; no integration is run, so this tests configuration forwarding without
creating output or initializing a logger. The default facade options must use
`SolverLogLevel::Off`.

**Result:** Passed in Debug with RST 0.5.0. All four native families
reported `Debug`; telemetry remained `Off`.

**Conclusion / limitation:** The calculator exposes native opt-in logging and
does not initialize a process logger, choose a file, or emit duplicate
iteration events. Actual log sink behavior and output policy belong to RST.

## Story 19: Release Correctness Campaign

**Name:** Focused release verification of the reaction-network IVP module.

**Question / hypothesis:** Do the module's current correctness stories and
public example build pass under the release profile with RST 0.5.0?

**Command:**

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactionNetworkIvpCampaign.ps1 -Suite correctness -BuildJobs 1
```

**Result:** Campaign `reaction-network-ivp-correctness-release-20261121-180250-892Z`
passed. The release library-test binary built, the public example check passed,
and all 28 targeted module tests passed (`0 failed`, `0 ignored`). The report and raw
logs are under the timestamped directory in `test_reports/`.

**Conclusion / limitation:** This verifies the implementation state before
Stories 18, 20, 21, and 22 were added; all four additions have since passed
in Debug. It is a focused module campaign, not a full-crate release test run.

## Story 20: Nonlinear Custom-Rate Analytic Reference

**Name:** Autocatalytic `A -> B` logistic trajectory.

**Question / hypothesis:** Can the full reaction-network assembly integrate a
nonlinear custom rate law and reproduce a closed-form trajectory across both
symbolic frontends and all current StiffAPI families?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features nonlinear_custom_rate_network_matches_logistic_reference_across_frontends_and_families -- --nocapture
```

**Fixture and acceptance:** `A -> B` with `r=3*A*B`, initial
`A(0)=1`, `B(0)=0.02 mol/m^3`, `T=700 K`, and `t=0..0.1 s`. The total
`M=A+B=1.02` is conserved. The independent reference is
`B(t)=M/[1+(A0/B0)*exp(-3*M*t)]`, `A(t)=M-B(t)`. Every published point for
BDF, Radau, Backward Euler, and LSODE2 using ExprLegacy and AtomViewNative
must be within `1e-4 mol/m^3` of the oracle and preserve `A+B` within
`2e-8 mol/m^3`. Backward Euler uses its supported native tolerance contract;
the shared step cap is `2e-4 s`.

**Result:** Passed in Debug for all eight frontend/family combinations.

**Conclusion / limitation:** This closes an ordinary nonlinear analytical
reference gap. It is not a genuinely stiff nonlinear gate and does not
establish concentration-scale robustness for nonlinear kinetics.

## Story 21: Nonlinear Rate Under Prescribed Heating

**Name:** Quadratic temperature-dependent rate under linear heating.

**Question / hypothesis:** Does substituting a changing prescribed temperature
into a nonlinear rate expression preserve the analytical concentration
trajectory across both symbolic frontends and all current StiffAPI families?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features nonlinear_temperature_rate_matches_closed_form_across_frontends_and_families -- --nocapture
```

**Fixture and acceptance:** `A -> B` with
`r=0.4*C_A*(T/500)^2`, `A(0)=1`, `B(0)=0`, `T(t)=500+100t K`, and
`t=0..1 s`. The exact integrated coefficient is
`0.4*(t + 0.2*t^2 + 0.04*t^3/3)`, so
`A(t)=exp(-integrated_coefficient)` and `B(t)=1-A(t)`. Every published
sample for BDF, Radau, Backward Euler, and LSODE2 with ExprLegacy and
AtomViewNative must be within `1e-4 mol/m^3`; `A+B` is conserved within
`2e-8 mol/m^3`, and published temperature must equal the prescribed linear
program within `1e-10 K`. The shared step cap is `2e-4 s`; Backward Euler uses
its native tolerance contract.

**Result:** Passed in Debug for all eight frontend/family combinations.

**Conclusion / limitation:** This closes nonlinear rate dependence on a
prescribed temperature for a linearly varying profile. It does not cover
temperature discontinuities or arbitrary piecewise interpolation. Story 22
adds a separate nonlinear stiff reference.

## Story 22: Stiff Nonlinear Network Across Concentration Scales

**Name:** Fast and slow autocatalytic logistic subsystems.

**Question / hypothesis:** Can the full network solver resolve nonlinear
kinetics with separated time scales and preserve scale metamorphism across
both frontends and every current StiffAPI family?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features stiff_nonlinear_logistic_network_matches_closed_form_across_scales -- --nocapture
```

**Fixture and acceptance:** Two independent `A -> B` subsystems use
`r_fast=100*A_fast*B_fast` and
`r_slow=0.1*A_slow*B_slow`, with each initial pair `[1, 0.02]` times the
concentration scale and `t=0..0.1 s`. Their kinetic scales differ by `1000`.
Initially, `k*M` is `102` and `0.102 s^-1`, while the nonzero eigenvalues
of the two logistic Jacobian blocks are `98` and `0.098 s^-1`.
The exact normalized solution for each pair is
`B/M=1/[1+50*exp(-k*M*t)]`, `A/M=1-B/M`, with `M=1.02` and `k=100` or
`0.1 s^-1`. The fixture crosses scales `1e-3`, `1`, and `1e3 mol/m^3`, both
frontends, and all four solver families. Every published sample must remain
finite/non-negative, normalized component error must be below `2e-4`, and
each pair total must be within `2e-5` of `M`. Backward Euler uses
`max_step=2e-6 s`; other families use `2e-3 s`.

**Result:** Passed in Debug for all 24 scale/frontend/family combinations.

**Conclusion / limitation:** This supplies a genuinely nonlinear stiff
reference with an independent closed-form oracle and a three-order
concentration-scale sweep. The two subsystems are intentionally uncoupled to
keep that oracle exact; this does not establish robustness for strongly
coupled nonlinear stiff networks or arbitrary scales.

## Story 23: SolidState F1 Network Contract

**Name:** Typed SolidState extent in the prepared reaction-network IVP.

**Question / hypothesis:** Does a single F1 conversion stage reproduce its
closed form through both symbolic frontends and every StiffAPI family, keep its conversion coordinate
separate from species concentrations, and reject resource/domain violations
without publishing a bad prepared result?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features reaction_network_ivp::solid_state_tests -- --nocapture --test-threads=1
```

**Fixture and acceptance:** `A -> B`, `q=1 mol/m^3`, `alpha(0)=0.1`,
`k=0.5 s^-1`, `T=700 K`, `t=0..1 s`. Independent oracle:
`alpha(t)=1-0.9*exp(-0.5*t)`, `C_A(t)=1-(alpha(t)-0.1)`.
The test asserts the RHS component values, the final conversion and
concentration for BDF, Radau, Backward Euler, and LSODE2 with ExprLegacy and
AtomViewNative, the preserved
initial conversion on prepared rebind, and a noncompeting F1 + MassAction
network, plus a conservative three-law network with a CustomExpr stage.
Separate assertions require typed errors for shared consumption,
insufficient stock, singular D1 at `alpha=0`, negative model coefficient,
out-of-domain RHS evaluation, and invalid or decreasing published conversion.
The catalog interior probe covers all 21 named/parameterized model variants.
A failed
initial-stock rebind must retain the last accepted solution.

**Result:** Debug run passed, 8/8 SolidState tests; full reaction-network IVP
suite passed, 40/40 tests.

**Conclusion / limitation:** The production facade now integrates a typed
SolidState extent within the same symbolic vector IVP. This is not yet an
all-catalog or all-frontend guarantee; native intermediate callback failures
still lack KiThe stage/model context. No release performance conclusion.

## Story 24: Scalar F1 And Network F1 Parity

**Name:** Same F1 conversion law through scalar and network StiffAPI routes.

**Question / hypothesis:** Does the new network route reproduce the existing
scalar solid-state result when both use identical `alpha(0)`, Arrhenius
parameters, linear heating, interval, BDF family, and tolerances?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features f1_network_matches_existing_scalar_solid_state_route -- --nocapture
```

**Fixture and acceptance:** `F1`, `alpha(0)=1e-5`, `T(t)=700+t K`,
`A=0.5 s^-1`, `E=5000 J/mol`, `t=0..1 s`, BDF with `max_step=0.01 s`,
`rtol=1e-8`, `atol=1e-10`. The final conversion difference must be
`<1e-5`; network `A` concentration must equal `1-(alpha-alpha0)` within
`1e-6 mol/m^3`.

**Result:** Passed in Debug. The scalar route retains its native `alpha0=1e-5`
convention; the network route supplies that value explicitly.

**Conclusion / limitation:** This verifies F1 route parity for one BDF
heating fixture. It does not compare all catalog models or prove identical
solver trajectories at every sample.

## Story 25: Failed-Attempt Cancellation Diagnostics

**Name:** Native telemetry and bounded partial trajectory after cancellation.

**Question / hypothesis:** Which StiffAPI families expose current-attempt
telemetry and partial samples after a mid-run cancellation signal, without
turning those samples into an accepted solution?

**Command:**

```powershell
$env:CARGO_BUILD_JOBS='1'
cargo test --lib --no-default-features mid_run_cancellation_reports_native_diagnostic_availability -- --include-ignored --nocapture --test-threads=1
```

**Fixture and acceptance:** `A -> B`, `k=0.5 s^-1`, `t=0..10 s`, requested
`max_step=1e-5 s`, Detailed telemetry. A separate thread requests cancellation
after 10 ms. The diagnostic test prints each native outcome and available
snapshot sizes; if a solve fails it must not publish a solution. A separate
non-ignored test pre-cancels all four families after a prior accepted solve
and asserts that no old telemetry or trajectory is misattributed to the new
attempt.

| Family | Observed outcome in two Debug runs | Failure telemetry | Partial points retained | Accepted solution |
| --- | --- | --- | ---: | --- |
| BDF | typed failure | yes | 32 | no |
| Radau | typed failure | no | 0 | no |
| Backward Euler | typed failure | yes | 32 | no |
| LSODE2 | completed before cancellation | not a failed attempt | 0 | yes |

**Result:** Ignored characterization passed on the second run; the first run
revealed that asserting cancellation for every wall-clock schedule was
incorrect. The pre-cancel correctness gate passes for all four families.

**Conclusion / limitation:** KiThe retains bounded raw diagnostics where RST
exposes them, with no acceptance or physical-validity claim. Radau needs an
upstream failure snapshot with attempt identity. This test does not establish
LSODE2 mid-run cancellation latency, and its timing is not a performance gate.
