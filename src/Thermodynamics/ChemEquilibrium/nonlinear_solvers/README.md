# Nonlinear Solver Backends

This directory contains concrete nonlinear solver implementations and the
mechanical infrastructure that dispatches, budgets, validates, and reports
their attempts.

## Contents

- `equilibrium_nonlinear.rs`: retained hand-written LM, NR, and trust-region
  implementations plus shared solver errors;
- `NR_Legacy.rs`: older compatibility solver implementation;
- `equilibrium_legacy_backend.rs`: adapter for the retained legacy backends;
- `equilibrium_rst_backend.rs`: RustedSciThe wrapper and symbolic prepared
  problems;
- `equilibrium_backend_adapter.rs`: backend-neutral request/result contract;
- `equilibrium_solver_policy.rs`: ordered fallback policy, budgets, and
  attempt reports;
- `equilibrium_validation.rs`: backend-independent acceptance gate;
- `*_tests.rs`: RST, matrix, and acceptance regressions.

## Boundary

This directory does not construct thermochemical problems or decide phase
lifecycle. `equilibrium_log_moles`, prepared runners, and phase-control
workflows remain in the parent directory as orchestration code. Keeping the
boundary explicit prevents a numerical backend from acquiring ownership of
physical state or publication.

## Stable Rust Paths

The parent module declares these files with `#[path]`, preserving public and
crate-private paths such as `ChemEquilibrium::equilibrium_rst_backend` during
the filesystem refactor.
