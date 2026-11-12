# MF-02B Diagnostic Report

Date: 2026-11-10  
Scope: diagnosis only. No production solver, equations, tolerances, bounds, or
acceptance rules were changed.

## Executive Summary

1. The reduced Damped hot-seed mapping is correct. Its round-trip error from a
   CorrectionNewton physical profile is `1.421085e-14` across
   `Teta,q,C0,J0,C1,J1`; the independently generated 301-interval hot
   trapezoid root packs with zero error.
2. Damped leaves both the hot CorrectionNewton profile and an independently
   solved hot root of the scalar trapezoid system. From the latter it reports
   three iterations and publishes the cold solution (`Tmax=300.187596 K`)
   instead of the hot solution (`Tmax=1800 K`). The seed defect is ruled out,
   but the exact mechanism cannot be localized to Newton globalization versus
   a difference in the production callback system: the current public path
   does not expose a supported full-state production residual/Jacobian
   evaluator or intermediate iterates.
3. MF-02B.1 found and fixed a row-indexing defect in the independent
   `F_delta`. After repeating continuation with a finite-difference-tested
   derivative and a bordered tangent update, the cold trapezoid fold remains
   near `delta=4.142789`. The old coarse estimates and new tangent-zero
   estimates differ by at most `9.2e-9`. This supports a fold of the discrete
   branches, not of the continuous BVP.
4. Existing Damped telemetry is insufficient to reconstruct accepted/rejected
   trial decisions. The observed `DampingTrial` events have `alpha=1` and
   `value_b=NaN`, without trial merits, explicit acceptance, or reason.

## Reproduction

The independent script imports only NumPy and SciPy, not KiThe or RST:

```powershell
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/diagnose_cold_branch.py
uv run --with matplotlib python src/SimpleReactorsBVP/literature_mf1963/plot_cold_branch.py
```

Focused Rust checks and complete MF module:

```powershell
cargo test --lib --no-default-features mf1963_damped_hot_seed_roundtrip_and_first_scalar_newton_audit -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_independent_trapezoid_hot_root_is_probed_through_damped -- --nocapture --test-threads=1
cargo test --lib --no-default-features mf1963_damped_corrector_newton_seed_layout_diagnostic -- --nocapture --test-threads=1
cargo test --lib --no-default-features literature_mf1963 -- --nocapture --test-threads=1
```

Result on 2026-11-10: `14 passed, 1 ignored`; the run took about 90 seconds.
These are Debug correctness/diagnostic runs, not performance measurements.

Files produced by the independent study:

- `diagnose_cold_branch.py`: independent collocation, trapezoid, continuation,
  and spectral calculations.
- `mf1963_mf_high_cold_branch_diagnostic.json`: full machine-readable results.
- `mf1963_mf_high_cold_branch_pseudo_arclength_all_grids.csv`: all
  pseudo-arclength points on 301/601/1201 grids.
- `mf1963_mf_high_cold_branch_pseudo_arclength.csv` and
  `mf1963_mf_high_cold_branch_collocation.csv`: compact branch data.
- `mf1963_mf_high_hot_trapezoid_301.csv`: independent hot discrete root used
  by the Rust Damped probe.
- `mf1963_mf_high_cold_branch_diagnostics.png`: rendered overview.

## A. Hot Seed and Damped Branch Selection

### Seed-layout audit

The full physical profile uses node-major storage (`node x state`), while the
Damped reduced initial state is state-major (`state x free_node`). For the
mixed boundary conditions, each state omits its own fixed endpoint: left-fixed
states omit node 0 and right-fixed states omit the final node. The corrected
helper normalizes physical mesh coordinates to `[0,1]` before interpolation.

The CorrectionNewton profile to Damped reduced-state round trip has maximum
absolute error `1.421085e-14`. A separate assertion checks each component's
free-node mapping; this rules out the prior reshape/packing error.

### CorrectionNewton hot profile on the independent trapezoid equations

The CorrectionNewton profile is a hot solution of its own route. Evaluated
against the independently implemented scalar trapezoid system at 301
intervals, it has:

| Quantity | Measured value |
|---|---:|
| Initial independent scalar residual L2 | `4.330684810e-3` |
| Initial residual infinity norm | `1.231944074e-3` |
| First independent Newton correction L2 | `4.361341533e-2` |
| Correction component range | `[-9.063850420e-3, 6.585935730e-3]` |
| Linear solve residual L2 | `6.166910237e-18` |
| Independent full-step trial residual L2 | `2.281375197e-8` |

Directional finite-difference checks of this independent scalar Jacobian gave
relative errors `3.018424e-8` (temperature state), `6.696247e-8` (slope), and
`5.860685e-8` (mixed direction). This validates the diagnostic scalar
Jacobian; it is not a parity test against the production six-state Jacobian.

### Exact independent hot trapezoid root passed to Damped

An independent Newton solve produced a hot root on the exact 301-interval
trapezoid scalar system. It was mapped to the six physical KiThe states and
passed to the fixed-mesh Damped facade with the prescribed configuration:
Sparse, ExprLegacy, Lambdify, `max_jac=1`, no adaptive refinement,
`abs/rel=1e-8`, and the existing bounds.

| Check | Result |
|---|---:|
| Independent scalar residual L2 | `7.321271866e-12` |
| Reconstructed scalar balance defects | `(2.814690e-14, 2.998490e-12)` |
| First-integral defect `max|eta-gamma*theta|` | `2.220446e-16` |
| Reduced-state seed packing error | `0` |
| Input profile `Tmax` | `1800 K` |
| Damped output `Tmax` | `300.187596 K` |
| Maximum profile temperature difference | `1.499863e3 K` |
| Maximum species-fraction difference | `0.9999086` |
| Damped iterations / residual calls / Jacobian calls | `3 / 6 / 2` |
| Damping trials / rejected trials / refinements | `3 / 0 / 0` |

Every observed trial event has `alpha=1` and `value_b=NaN`. The actual Legacy
acceptance metric in RST 0.4.32 is the **undamped Newton correction norm**
`s1 = ||Delta||_2`; a trial is accepted when `s1 < 1` or `s1 < s0`, where
`s0` is the preceding correction norm. It is not a residual-decrease merit
test. RST also computes/retains a residual norm for convergence reporting, but
that value is not present in these trial events. Zero rejection counters and
full-step trial events are consistent with the cold transition, but there is
no explicit accepted flag in the event record and no intermediate profile, so
the exact iteration where the hot branch is lost is unknown.

### Production-discretization parity limitation

The independent scalar model uses the fixed domain
`xi_max=0.5427079731331537`, `theta'=p`,
`p'=p-delta*exp(theta/(1+beta*theta))*(1-gamma*theta)/(1+beta*theta)`,
`theta(0)=0`, `p(xi_max)=0`, with trapezoid interval equations. The mapped
root also satisfies the inspected six-state closure/first-integral relations
to near machine precision.

However, a componentwise comparison
`F_Damped(u_hot) == F_CorrectionNewton(u_hot)` and corresponding Jacobian
comparison was **not** obtained. The public RST `fun/jac` references in the
constructed Damped object are internal chunk callbacks, not a full packed
state evaluator. Calling them directly with the packed full state is invalid
(type mismatch before conversion; after sparse conversion, chunk-local
indexing panics). This route is not a supported diagnostic API and is not used
as a regression. Consequently the report does not claim exact equality of
production residual definitions/scaling or six-state Jacobians.

The independent hot-root experiment does establish that seed layout and
CorrectionNewton-to-trapezoid discretization mismatch alone are insufficient
to explain Damped's cold publication: Damped also leaves a verified root of
the independent scalar trapezoid model. Remaining plausible causes include a
difference in production callback/discretization details or a globalization,
state reconstruction, and stopping interaction inside Damped. The available
evidence cannot distinguish these causes.

### MF-03 implementation boundary

The KiThe lockfile uses RustedSciThe `0.4.32`; its source is available here
only as Cargo's registry package. The separate local RustedSciThe checkout is
version `0.2.9`, so it is not a valid source tree for implementing or
validating an API change against the locked dependency. The only available
full-state-adjacent callbacks in `0.4.32` are chunk callbacks; invoking them
as if they were a stable full-state evaluator is invalid and can panic.

Implementing MF-03 correctly requires an upstream-compatible RST change that
exposes its existing internal full-state assembly in a typed, non-mutating
diagnostic/test entry point and adds opt-in accepted-step telemetry. KiThe
must not duplicate that assembly, and the Cargo registry cache must not be
patched as a substitute for an upstream source change. This is the point where
local work pauses pending the RST API change. The paste-ready request above is
the handoff; the Damped criterion remains unchanged.

## B. Cold Branch and Fold Evidence

### MF-02B.1 correction to parameter derivative

The old `F_delta` loop placed its first derivative value in flattened residual
row 0, which is the left boundary condition, and omitted the last `p`-balance
row. The corrected residual layout is interleaved: `r[0,0]` is the left
boundary condition; the `p` balances are `r[i,0]` for `i=1..N`, including the
last interval; theta balances and the right boundary condition occupy the
odd flattened rows. The corrected derivative is

```text
(F_delta)[0,0] = 0
(F_delta)[i,0] = h/2 * (g(theta[i-1]) + g(theta[i])) for i=1..N
(F_delta)[i,1] = 0 for i=0..N
```

It is implemented once in `trapezoid_parameter_derivative()` and reused for
the initial tangent, augmented Newton corrector, and tangent updates. A full
vector centered finite-difference test runs on two grids and nonuniform
temperature/slope states:

| Intervals | Relative `F_delta` error | Maximum absolute component error |
|---:|---:|---:|
| 7 | `3.558012e-11` | `3.771156e-7` |
| 13 | `4.425817e-11` | `3.746163e-7` |

The test explicitly verifies that the first and last boundary rows are zero,
all theta rows are zero, and both the first and last `p`-balance derivative
rows are populated.

The old tangent update solved `J du = -F_delta` at every point, which is not a
robust tangent construction near a singular `J`. The corrected continuation
updates tangent through the bordered system using the previous tangent as the
normalization row. When the tangent's parameter component changes sign, it
retries only that crossing interval with a halved arc step down to `1e-5`,
then interpolates the zero of `d(delta)/ds`.

### Old and corrected fold estimates

The old values below were the maximum sampled `delta` from the former coarse
arc trace; the corrected values interpolate the zero of the continuation
tangent component after local crossing-step refinement.

| Intervals | Old sampled `delta_fold` | Corrected tangent-zero estimate | Difference |
|---:|---:|---:|---:|
| 301 | `4.142800571129940` | `4.142800580301762` | `9.172e-9` |
| 601 | `4.142791344441408` | `4.142791345592427` | `1.151e-9` |
| 1201 | `4.142789025103465` | `4.142789031353778` | `6.250e-9` |

The corrected estimates preserve the same mesh-convergent branch picture; the
old `F_delta` defect was real, but did not materially move these estimates.
The final arc step bracketing the tangent sign change is at most `1e-5`. The
original residual and extended corrector infinity norms at the nearest
pre-fold roots are at most `6.94e-11` for 301 intervals, `2.14e-12` for 601,
and `5.03e-12` for 1201.

Independent pseudo-arclength continuation solves the trapezoid algebraic
system with fixed domain and state scaling. The following state/spectrum data
are from the nearest pre-fold corrector root (the fold parameter itself is
interpolated from the tangent sign change):

| Intervals | `delta_peak` | `Tmax` (K) | inlet `p(0)` | scaled `sigma_min(J)` | estimated condition |
|---:|---:|---:|---:|---:|---:|
| 301 | `4.142800580301762` | `322.253021` | `3.988537521` | `1.242381e-10` | `4.4642e11` |
| 601 | `4.142791345592427` | `322.253013` | `3.988524144` | `2.121441e-11` | `5.2201e12` |
| 1201 | `4.142789031353778` | `322.253011` | `3.988520776` | `4.889676e-12` | `4.5258e13` |

The fixed scaling is `J_scaled = diag(row_scale) J diag(column_scale)`.
Boundary rows have row scale 1; trapezoid ODE rows are scaled by `1/h`.
Temperature and slope state columns both use fixed scale `0.05`. The
singular-value estimates are sparse estimates from SciPy and are tied to this
explicit scaling; raw values and condition estimates are also retained in the
JSON report. At these nearest roots, normalized tangent parameter components
are `7.70e-9`, `1.76e-9`, and `4.25e-10`; corresponding scaled state
derivatives with respect to `delta` are approximately `1.30e8`, `5.67e8`, and
`2.35e9`. This is the expected sharp growth near the discrete turning point.

At `delta=4`, the same fixed-scaled minimum singular-value estimates are
`3.049424e-4`, `2.159639e-4`, and `1.528304e-4` on those grids. The fold
estimates converge to approximately `delta=4.142789`, and all three discrete
systems have passed their cold fold before `delta=4.143`. This supports a
**discrete cold-branch fold** and explains why ordinary cold continuation
becomes singular in the observed Damped bracket. It does not prove the
continuous BVP has a fold at the same value.

The independent hot trapezoid branch persists through `delta=4.143` on all
three grids with `Tmax` near `1800 K`. Hence the cold fold does not remove the
hot branch.

For the continuous SciPy collocation model, cold continuation reaches
`delta=4.1427` (`Tmax=322.099 K`, 19036 nodes). At `delta=4.1428`, `solve_bvp`
returns “maximum number of mesh nodes is exceeded.” This is a resource/mesh
limit, not a fold detection. Continuous-branch status at and beyond that point
remains unresolved.

## C. Damped Telemetry Assessment

For a successful exact-hot-seed run, Detailed logging currently exposes only
three `DampingTrial` entries:

| Iteration | Trial | `value_a` | `value_b` | Rejection event |
|---:|---:|---:|---:|---|
| 1 | 1 | `1.0` | `NaN` | no |
| 2 | 2 | `1.0` | `NaN` | no |
| 3 | 3 | `1.0` | `NaN` | no |

The counters say 3 trials and 0 rejections, but there is no explicit accepted
event, current/trial correction norm, the threshold comparison operands,
residual values already available to the solver, rejection reason, per-state
correction, Jacobian age/refresh event, or intermediate solution. Therefore
these data cannot explain the hot-to-cold transition. The `value_b` contract
should not be silently redefined at the KiThe call site.

### Paste-ready upstream RST request

> Please extend Damped's opt-in detailed decision telemetry without changing
> its numerical algorithm, merit function, damping schedule, Jacobian refresh
> policy, bounds behavior, or stopping criteria. For each Newton iteration,
> expose the current correction metric and previous comparison value, residual
> norm if already available, Jacobian age and refresh status, and accepted
> damping factor plus final convergence metrics. For each trial, expose the
> attempted alpha, correction metric actually used by the acceptance rule,
> that rule's threshold/operands and outcome, residual norm if already
> computed, accepted/rejected, and a typed reason for rejection (criterion not
> met, bounds, nonfinite evaluation, linear solve failure, or other). Also expose the iteration
> termination reason. Do not perform additional residual/Jacobian evaluations
> solely for telemetry; record values already available to the solver. Keep
> the event schema bounded and ensure Telemetry Off has no material extra
> work. Please document existing `DampingTrial.value_a/value_b` semantics
> before changing either field.

Separately, a test-only or public diagnostic entry point is needed to evaluate
the prepared full residual and Jacobian at a supplied full state. It must
perform the same unpacking, boundary insertion, chunk assembly, scaling, and
ordering as the production solve. The currently exposed chunk callbacks do
not fulfill that contract.

## Final Classification

| Question | Finding |
|---|---|
| A. Why does Damped leave the hot branch? | Not caused by the corrected seed packing. It also leaves an independent hot root of the scalar trapezoid equations. Exact cause remains unresolved without production residual/Jacobian evaluation and per-iteration accepted-state telemetry. |
| B. Does a cold solution exist at `delta=4.143`? | The independent 301/601/1201 trapezoid systems have a converged cold fold near `4.142789`, so their tracked cold branches turn before `4.143`. Existence/nonexistence for the continuous BVP is unresolved. |
| C. Is current telemetry sufficient? | No. It reports trial alpha/counters but not merits, acceptance, reasons, or intermediate states. |

No solver fix is proposed by this diagnostic. Any production algorithm change
should be a separate task after the RST callback/telemetry limitations are
addressed.
