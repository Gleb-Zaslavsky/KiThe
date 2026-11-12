# BVP-01: Analytical Transport Benchmark

This folder owns an independent analytical correctness fixture for the
production gas-phase BVP reactor. The benchmark must validate KiThe's existing
equation assembly, boundary conversion, solver routes, and dimensional
publication; it must not replace or alter production equations.

## Audit Against The Current Implementation

- [x] The canonical state order is `Teta, q, C0, J0, C1, J1, ...`. It is
  assembled in `createBVP.rs` as two thermal states followed by each species
  concentration and flux state.
- [x] With one inert/source-free reaction, the production equations reduce to
  the task's independent transport equations:

  ```text
  Teta' = q / Lambda
  q'    = Pe_q q
  C_i'  = J_i / K_i
  J_i'  = Pe_D_i J_i
  ```

  Here `q` and `J_i` are internal scaled flux states, `K_i = rho D_i`,
  `Pe_q = m Cp L / Lambda`, and `Pe_D_i = m L / K_i`. The production
  equation builder also includes heat/species source terms; this fixture must
  demonstrate that they are exactly zero for the chosen zero-rate reaction.
- [x] `set_solver_BC()` interprets `T` as inlet temperature, `q` as the
  physical outlet heat flux and `<species>_flux` as the physical outlet mass
  flux. It converts them using `Teta=(T-dT)/T_scale`,
  `q_hat=q_physical*L/T_scale`, and `J_hat=J_physical*L`. Missing flux fields
  default to `1e-10`, so cases that require exact zero flux must pass explicit
  zero values.
- [x] At `Tm=298.15 K`, transport preparation uses
  `rho_ref=M*P/(R_G*298.15)` and `K_i=D_i*rho_ref`; `R_G=8.314`. Choosing
  `P=rho*R_G*Tm/M` therefore realizes the requested reference density, subject
  to floating-point rounding and the fixture's parsed molar mass.
- [x] The solver publishes physical coordinates and fields by default:
  `x=L*z`, `T=dT+T_scale*Teta`, `q_physical=q_hat*T_scale/L`,
  `J_physical=J_hat/L`; mass fractions are unchanged. The test must inspect
  this production publication path, not manually apply the transform twice.
- [x] Existing tests live in flat `reactor_bvp_*_tests.rs` modules, while
  `reactor_bvp_test_support.rs` contains reusable physical fixtures. This
  benchmark gets an isolated directory so its analytic reference and staged
  diagnostics do not become production code or mix with HMX stories.
- [x] Confirm end-to-end that a zero-pre-exponential reaction (`A=0`,
  `n=0`, `E=0`, `Q=0`) survives `fast_react_set`, kinetic preparation,
  `check_task`, and symbolic BVP assembly and produces identically zero source
  expressions. `fast_react_set` checks finiteness rather than positivity of
  `A`, but this is not by itself proof that the full pipeline accepts the
  fixture. Do not substitute a small nonzero `A` if exact zero is rejected;
  report the concrete blocker and decide on a test-only construction that
  still exercises the production equation builder.
- [x] Confirm the proposed pseudo-species formula parses and obtains the
  intended molar mass from the actual mechanism/data path. The proposed
  identical formula for both pseudo-species makes this a transport fixture,
  not a chemical mass-conservation test.

## Independent Analytical Reference

Use `z=x/L` and define

```text
F(p,z) = exp(-p) * expm1(p*z) / p,  p != 0
F(0,z) = z.
```

For right-boundary internal flux `q_R` or `J_R_i`, and inlet values
`Teta_L`, `C_i,L`, the reference is

```text
q(z)       = q_R * exp(p_q*(z-1))
Teta(z)    = Teta_L + (q_R/Lambda) * F(p_q,z)
J_i(z)     = J_R_i * exp(p_D_i*(z-1))
C_i(z)     = C_i,L + (J_R_i/K_i) * F(p_D_i,z).
```

The test implementation must be independent of KiThe's symbolic expressions
and must handle `p=0` without division by zero. For larger positive `p`, use a
numerically stable equivalent that avoids forming `exp(p)`; test the actual
parameter range through `p=30`. Validate reference values at endpoints and
interior points before using them as solver oracles.

First-integral diagnostics for each computed profile are

```text
q(z)   - q(0)   - Lambda*Pe_q*(Teta(z)-Teta(0)) = 0
J_i(z) - J_i(0) - K_i*Pe_D_i*(C_i(z)-C_i(0))    = 0.
```

These are supplemental checks; componentwise profile error against the closed
form remains the primary oracle.

## Cases And Gates

| Case | Parameters / purpose | Analytic gate |
| --- | --- | --- |
| BVP-01a | `Pe=0` limiting reference plus a production constant-profile case with explicit zero right fluxes at `Pe=1` | The independent `F(0,z)` branch is tested directly. The facade test uses exact zero boundary fluxes and gets constant profiles; it does not mislabel the production case as `Pe=0`. |
| BVP-01b | `Pe_q=Pe_D0=Pe_D1=1`, `Lambda=0.1`, `D0=D1=1e-4`, physical `q_R=1000 W/m^2`, `J_R0=+0.002`, `J_R1=-0.002 kg/(m^2 s)` | `T(1)=361.3620558829 K`; `C0(1)=0.7264241118`, `C1(1)=0.2735758882`; `C0+C1=1` for this equal-transport setup. |
| BVP-01c | As 01b except `D1=2e-4`, hence `Pe_D1=0.5` | `C0(1)=0.7264241118`, `C1(1)=0.3213061319`; total fraction is about `1.0477302437`. Do not assert unit sum: independent Fick closures with unequal diffusivities do not preserve it. |
| BVP-01d | `Pe_q=Pe_D0=Pe_D1=30`, `Lambda=1/300`, `D0=D1=1/300000`, physical `q_R=250 W/m^2`, `J_R0=+0.0015`, `J_R1=-0.0015` | `T(1)=323.15 K`, `C0(1)=0.75`, `C1(1)=0.25`; compare full profiles and resolve the outlet boundary layers under mesh refinement. |

Parameter notes from the proposal: `M=0.060052 kg/mol`, `Tm=298.15 K`,
`Cp=1000 J/(kg K)`, `m=0.01 kg/(m^2 s)`, `L=0.01 m`, `dT=298.15 K`,
`T_scale=100 K`, `rho=1 kg/m^3`, and `P=rho*R_G*Tm/M` (approximately
`41277.877506 Pa`). The test must assert the actual derived density,
`K_i`, and Peclet numbers before solving; this catches inconsistencies from
formula parsing, unit conversion, or constants.

**Open contract for BVP-01a:** the current task validation requires positive
mass flux and strictly positive `Lambda` and diffusion coefficients. Thus
exactly zero Peclet numbers are not obviously representable through the public
physical task contract. Test the analytic `p=0` branch directly in the
independent reference helper. The production constant-profile case instead
uses positive Peclet numbers with exact zero right-boundary fluxes; it tests
the zero-flux physical solution without claiming exact-zero Peclet coverage.

## Solver And Mesh Matrix

- [x] Exercise currently supported Damped, Frozen/Naive, and BVP_sci routes
  with the same physical fixture, boundary conditions, and acceptance
  quantities. Record route errors as data; do not weaken a shared analytic
  tolerance to make every isolated route pass. Existing historical evidence
  says Frozen/Naive had a physical-acceptance gap on compact HMX, so this simple
  source-free case is useful evidence but must not be presumed to close that
  separate gap.
- [x] Use a basic non-analytic initial guess for the primary solver
  correctness gate. An analytic initial guess may be added as a separate
  diagnostic to distinguish solver failure from equation/coordinate errors;
  it must not be the only passing route.
- [x] For 01a-01c, establish a modest fixed-grid baseline and report max
  absolute/relative profile errors, endpoint errors, first-integral errors,
  boundary-condition residuals, solver termination, iterations, and final
  mesh. Set tolerances from observed mesh/solver convergence, not from the
  proposal's suggested values alone.
- [x] For 01d, compare at least two successively refined fixed meshes for
  Damped/Frozen where supported, and the adaptive BVP_sci mesh. Report error
  versus the analytic solution over the whole domain and separately near
  `z=1`; endpoint-only comparison is insufficient for a boundary-layer case.
- [x] Exercise Damped's own BVP_api adaptive-grid settings explicitly: choose
  `DoublePoints`, cap refinements, enable counters, and assert that the native
  telemetry reports actual refinements and lower analytic error than the
  fixed seed mesh. This analytical Damped fixture explicitly selects the
  `ExprLegacy` symbolic frontend to avoid the separately tracked RST AtomView
  adaptive-grid panic. Do not infer Damped adaptivity from a manual mesh sweep.
- [x] Confirm each route receives the intended grid and guess shape. Current
  adapter contracts differ: Damped uses interval-count columns, while
  Frozen/Sci use nodal columns; Frozen removes the prescribed endpoint from
  its internal seed. Do not compare nominal `n_steps` without recording actual
  mesh nodes/intervals.
- [x] Verify dimensional output (`x`, `T`, physical `q/J`) on at least one
  facade/file-to-solution story, as well as solver-space fields if exposed.
  Assert values and units, not just labels or finite output.

## Implementation And Reporting Checklist

- [x] Add an independent analytic reference helper and parameter constructors
  under this folder; keep the helper free of calls to production equation
  assembly.
- [x] Add a staged smoke test for parameters/transforms, then exact-zero
  source assembly, then one moderate (`Pe=1`) solver route, then the route
  matrix, then the `Pe=30` mesh study. Each stage should print enough context
  to identify the first divergence.
- [x] Add named correctness stories to `src/SimpleReactorsBVP/STORY_TESTS.md` with
  hypothesis, exact command, acceptance definition, observed route results,
  and limitations (especially BVP-01a representability and 01c non-unit
  fraction sum).
- [x] Keep this as a correctness/verification campaign; no performance gate
  is required. Solver timings may be included as diagnostics but are not an
  acceptance criterion.
- [x] Run focused tests first, then the relevant BVP correctness matrix. Do
  not change production equations, solver tolerances, or GUI behavior merely
  to force agreement with this reference; investigate any mismatch at the
  first failing stage.

## Current Status

Code-level equation audit, zero-rate kinetic preparation, pseudo-species mass
preparation, moderate-case Damped/Frozen/Sci correctness, unequal-diffusivity
Sci behavior, adaptive `Pe=30` Sci behavior, and explicit Damped BVP_api grid
refinement through `ExprLegacy` are covered by tests. The Damped adaptive story selects
`GridRefinementMethod::DoublePoints`, sets a two-refinement limit, and asserts
the reported refinement count and improvement over its fixed-grid seed. The
targeted Debug run and all 152 non-ignored `SimpleReactorsBVP` lib tests were green
before this final adaptive Damped story was added; rerun the focused test.
Exact production `Pe=0` remains unavailable because the task contract requires
strictly positive mass flux, thermal conductivity, and diffusivities; the
zero-limit helper and a separate exact-zero-flux facade case cover the
corresponding limiting behaviors without changing that claim.
