# Fixed-Pressure, Fixed-Enthalpy Workflow

This directory groups the `P,H = const` capability: its typed constraint
boundary, monolithic and nested formulations, thermochemistry preparation,
continuation, and presentation reports.

## Contents

- `equilibrium_constraints.rs`: typed enthalpy and temperature constraints;
- `equilibrium_ph_thermochemistry.rs`: format-agnostic thermochemistry
  capabilities and provenance;
- `equilibrium_ph_formulation.rs` and `equilibrium_ph_monolithic.rs`: coupled
  composition-temperature formulation and runner;
- `equilibrium_ph_nested.rs`: safeguarded scalar reference/recovery route;
- `equilibrium_ph_options.rs`: route and acceptance options;
- `equilibrium_ph_range.rs` and `_presentation.rs`: transactional continuation
  across nearby `P,H` targets;
- `equilibrium_ph_workflow.rs`: public resolved-data facade.

## Boundary

The fixed-`P,T` solver, phase-control lifecycle, and generic temperature-range
workflow remain in the parent directory. The P,H workflow consumes those
canonical capabilities; it must not fork the production phase or nonlinear
solver orchestration.

## Stable Rust Paths

The parent module uses explicit `#[path]` declarations. Existing imports such
as `ChemEquilibrium::equilibrium_ph_workflow` therefore remain unchanged.
