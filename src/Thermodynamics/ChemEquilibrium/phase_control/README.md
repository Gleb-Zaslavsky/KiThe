# Phase Control

This directory owns the bounded outer lifecycle around a fixed-set canonical
equilibrium solve.

## Contents

- `equilibrium_active_set.rs`: immutable projection into one active phase set;
- `equilibrium_phase_stability.rs`: TPD and constrained phase-stability
  primitives;
- `equilibrium_workflows.rs`: retained mutable phase-control orchestration and
  compatibility helpers;
- `prepared_phase_control_runner.rs`: immutable prepared-runner path;
- `equilibrium_workflow_tests.rs`: lifecycle, hysteresis, rollback, and legacy
  compatibility characterization.

## Boundary

The directory decides whether a candidate phase set can be activated or
deactivated, but it does not own resolved thermochemistry lookup or the public
`ResolvedPhaseSystem` bridge. The bridge builds immutable canonical problems;
nonlinear backend dispatch remains in `../nonlinear_solvers/`.

Rust module names remain stable through explicit `#[path]` declarations in
the parent module.
