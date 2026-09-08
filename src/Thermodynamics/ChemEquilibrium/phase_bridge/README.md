# Resolved-Phase Bridge

This directory is the typed boundary between the phase/data subsystem and the
canonical chemical-equilibrium engine.

## Contents

- `phase_equilibrium_problem.rs`: converts a resolved phase system into a
  canonical equilibrium problem with provenance and timing;
- `phase_equilibrium_solution.rs`: immutable phase-qualified accepted result;
- `phase_equilibrium_workflow.rs`: narrow public fixed-`P,T` facade;
- `phase_equilibrium_problem_tests.rs` and `equilibrium_phase_bridge_tests.rs`:
  resolved-data boundary contracts;
- `equilibrium_multiphase_story_tests.rs`: offline public-facade stories.

## Boundary

This bridge never owns phase activation rules or nonlinear backend algorithms.
It receives already-resolved data, prepares the canonical problem, delegates
to phase control and the solver policy, then publishes one immutable result.

The original Rust module paths are retained through `#[path]` declarations in
`ChemEquilibrium.rs`.
