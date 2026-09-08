# Equilibrium-Constant Validation

This directory owns the independent equilibrium-constant (`K_eq`) path used
to validate accepted canonical equilibrium solutions on its intentionally
small supported domain.

## Contents

- `equilibrium_constant_problem.rs`: reaction-extent problem and activity
  model;
- `equilibrium_constant_solver.rs`: independent extent solver;
- `equilibrium_constant_validation.rs`: validation modes, tolerances, and
  typed evidence;
- `equilibrium_constant_cross_validation.rs`: comparison and status reporting
  for canonical versus independent results;
- `*_tests.rs`: solver/domain regressions.

## Boundary

This is a validation subsystem, not a second production multiphase solver.
The canonical log-mole/phase-control path remains the only production source
of equilibrium solutions. `K_eq` is run only when the selected problem falls
inside its explicit applicability contract.

`equilibrium_reaction_basis.rs` remains in the parent directory because it is
a general typed reaction-basis contract shared by the canonical engine and
this independent validator.

## Stable Rust Paths

`ChemEquilibrium.rs` uses explicit `#[path]` declarations, so existing imports
such as `ChemEquilibrium::equilibrium_constant_problem` and test filters keep
their established names while the source files are grouped here.
