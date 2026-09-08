# Postprocessing and Logging

This directory owns read-only observation of accepted equilibrium results. It
collects optional timing and diagnostic evidence, renders it for logs or user
interfaces, and derives interpolated temperature-series views without changing
the solver state.

## Contents

- `equilibrium_timing.rs`: opt-in stage timings and immutable timing reports;
- `equilibrium_diagnostics.rs`: bounded typed lifecycle/TPD/backend events and
  optional live sinks;
- `equilibrium_diagnostics_display.rs`: explicit text and `log` facade output
  for structured diagnostics, plus quiet accepted-solution execution summaries
  that retain phase transitions, extensive-normalization recovery, and P,H
  route fallback evidence;
- `equilibrium_presentation.rs`: generic row-oriented accepted-solution views;
- `equilibrium_display.rs`: display-only filtering, units, and number formats;
- `equilibrium_temperature_postprocessing.rs`: PCHIP resampling of stable
  temperature series;
- `equilibrium_range_presentation.rs`: P,T range projections that preserve
  phase-transition boundaries and forbid smoothing across them.

## Boundary

No module here changes a candidate, cache, phase set, continuation seed, or
publication decision. Diagnostics are disabled by default; timing avoids clock
reads unless explicitly enabled. `format_solution_execution_summary()` and
`format_ph_solution_execution_summary()` remain useful with verbose diagnostics
disabled: they read only immutable accepted reports and expose committed
phase-control history, extensive-normalization recovery, and P,H route choices.
Interpolation receives immutable accepted points and may only produce a
separate render/export grid.

`ph/equilibrium_ph_range_presentation.rs` remains with the P,H workflow
because it exposes P,H-specific continuation and enthalpy evidence. Execution
cancellation/progress and reproducibility capsules remain parent-level runtime
contracts rather than presentation concerns.

## Stable Rust Paths

The parent module uses explicit `#[path]` declarations, so callers continue to
import `ChemEquilibrium::equilibrium_timing`,
`ChemEquilibrium::equilibrium_diagnostics`, and
`ChemEquilibrium::equilibrium_temperature_postprocessing` unchanged.
