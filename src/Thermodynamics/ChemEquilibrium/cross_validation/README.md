# Cross-Validation Modules

This directory groups the independent validation infrastructure for pure-phase
P,T and P,H scenarios:

- `phase_boundary_validation.rs` owns independent boundary mathematics;
- `phase_boundary_production_adapter.rs` converts immutable production
  evidence into validation inputs;
- `real_pure_phase_fixtures.rs` resolves shared offline chemistry;
- `pure_phase_*` modules contain P,T/P,H validators and their lifecycle/live
  regression stories;
- `real_boudouard_ph_validation.rs` and its tests own the stricter I4
  Boudouard `P,H` reference-state chain without turning JANAF I5 data into an
  external P,H reference table;
- `phase_boundary_*_tests.rs` covers the lower-level boundary contracts.
The independent K_eq validation subsystem has its own
`../equilibrium_constants/` directory. It is separate because it validates a
different mathematical formulation, whereas this directory owns pure-phase
boundary and lifecycle evidence.

The module names remain re-exported from `ChemEquilibrium.rs` through explicit
`#[path]` declarations. This preserves established internal paths while
keeping the filesystem organized around the validation responsibility.
