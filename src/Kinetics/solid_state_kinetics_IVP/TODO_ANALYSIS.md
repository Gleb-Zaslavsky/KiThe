# Solid-State Kinetics IVP TODO

This folder owns the solid-state kinetic model catalog and its scalar IVP
workflow. The public path remains `Kinetics::solid_state_kinetics_IVP` through
`mod.rs`.

## Scope and module boundaries

- [x] Move the former monolithic module into a folder while preserving its
  public module path.
- [x] Separate symbolic kinetic models from the IVP lifecycle and keep the
  existing tests grouped by concern.
- [x] Replace `ODE_api2::UniversalODESolver` as the execution backend with
  `numerical::StiffAPI`. The state dimension is one and the equations can be
  stiff; prepare requests with `StiffPreset::solve_1d()` and retain the
  selected native stiff solver's typed errors/results.
- [x] Route BDF, Radau, Backward Euler, and LSODE2 selections into their
  corresponding StiffAPI requests. Non-stiff RK45 choices were removed from
  the solid-state CLI/GUI selector rather than silently mapped to another
  method. The selector and constructor use native `StiffSolverKind`.
- [x] Exercise the solid-state GUI with `egui_kittest`: all 21 model names and
  formulas render within the catalog viewport without row overlap, all four
  native stiff solver choices are selectable, and the solver selector and Run
  action remain within a compact 820x620 window. Problem inputs are arranged
  compactly, and model actions have an unambiguous accessible label.
- [x] Replace string-keyed solver settings with
  `SolidStateIvpSolverOptions` throughout production callers.
- [x] Remove the string-keyed solver-parameter adapter and the `ODE_api2` solver
  selector compatibility boundary. Production callers now use
  `SolidStateIvpSolverOptions` and `StiffSolverKind` directly; non-stiff methods
  are unrepresentable in this scalar stiff-solver facade.
- [x] Expose opt-in StiffAPI telemetry through this facade, defaulting to Off,
  and preserve the native per-family snapshot rather than flattening unlike
  counters into a fabricated common record.
- [x] Keep `solver.rs` focused on configuration, validation, and dispatch;
  trajectory-to-`IVPSolution` materialization now lives in `solution.rs`,
  independent of the native solver result type.
- [ ] Add a typed domain error that preserves StiffAPI's source error. Migrate
  public `Result<_, String>` boundaries in a separate API-conscious step.
- [x] Invalidate cached solver/result state whenever problem, heating,
  Arrhenius, model, solver, or stopping configuration changes.
- [x] Remove duplicated model name/parameter state. `KineticModel` is now the
  authoritative representation; the RHS expression is derived from it and the
  current thermal/Arrhenius inputs when preparing each solve.
- [x] Avoid rebuilding/lambdifying the selected model once per output sample:
  `solution.rs` builds one evaluator from the stored model and reuses it across
  the trajectory.
- [x] Replace the arbitrary map stopping condition with a typed scalar
  conversion threshold and validate it before preparation. Current common
  output truncation is not event localization and can only retain the first
  accepted sample past the requested threshold; LSODE2 also receives its native
  stopping condition.
- [ ] Separate plotting/file output from numerical integration and cover its
  error behavior independently.

## Correctness audit

These are investigation items, not approved formula changes. Preserve current
formula behavior during structural/backend migration and add independent
numeric references before changing any expression.

- [x] Audit the Avrami-Erofeev `A2`, `A3`, and `A4` exponents against their
  documented rate laws and `get_g`; corrected the three negative exponents to
  the positive powers required by those independent integral forms.
- [x] Audit `P3`: its implemented `a^0.75` rate and `get_g = 4*a^0.25` are
  consistent. Corrected the stale `R3` formula label to `(1-a)^(2/3)`.
- [x] Verify the explicit `JMA` and `SB` derivatives against independently
  evaluated central finite differences; corrected both expressions.
- [x] Added numeric rate and derivative checks for all 21 model variants at an
  interior conversion, plus rate checks near `a=0` and `a=1`. The 15 models
  with analytical integrals also satisfy `g'(a) * f(a) = 1` numerically.
- [ ] Validate model parameter vector length and finite/domain constraints at
  the public boundary. `KineticModelNames::from_names_to_expr` must not panic
  on missing parameters; duplicate validation policy should be consistent
  with `KineticModel::from_name_and_params`.

## Lifecycle and robustness stories

- [x] Add a story changing time, heating, and Arrhenius settings after
  `set_model`; its trajectory matches a freshly configured solve. Published
  rates are checked against that same final configuration.
- [ ] Add stories proving a failed reconfiguration or solve cannot expose a
  stale result from an earlier successful solve.
- [ ] Check all scalar inputs for finiteness as well as sign; reject NaN and
  infinities with field-specific errors.
- [ ] Check published trajectory dimensions, finite values, and conversion
  range policy; document whether overshoot beyond `[0, 1]` is rejected,
  tolerated, or clipped (never silently clip without a stated contract).
- [ ] Add cross-solver correctness stories for representative stiff and
  non-stiff regimes using the same model, input, and acceptance criteria.
- [ ] Add an analytical/reference solution story for a simple first-order
  isothermal case before changing solver defaults or tolerances.

## Suggested implementation order

1. Establish the StiffAPI scalar solve path and typed solver selection while
   preserving existing physical equations.
2. Add focused backend-parity/correctness stories and update public callers.
3. Formula audit is now covered by independent value, derivative, and integral
   consistency assertions; preserve these gates for future model edits.
4. Consolidate mutable configuration and result invalidation.
5. Introduce typed errors and finish output/result separation.
