# ChemEquilibrium Behavior Inventory

This file is a current snapshot of the public and test-facing surface of the
`Thermodynamics/ChemEquilibrium` subsystem. It is intentionally descriptive,
not aspirational. The goal is to help future refactors preserve the right
behavior and retire the right legacy paths.

## Canonical modules

- `equilibrium_calculator.rs`
- `equilibrium_problem.rs`
- `equilibrium_prepared_runner.rs`
- `equilibrium_element_inventory.rs`
- `equilibrium_multiphase_domain.rs`
- `equilibrium_validation.rs`
- `equilibrium_solver_policy.rs`
- `equilibrium_rst_backend.rs`
- `equilibrium_backend_adapter.rs`
- `phase_bridge/phase_equilibrium_workflow.rs`
- `phase_control/prepared_phase_control_runner.rs`
- `ph/equilibrium_ph_workflow.rs`
- `equilibrium_temperature_range.rs`
- `equilibrium_temperature_postprocessing.rs`
- `equilibrium_nonlinear.rs`
- `equilibrium_ids.rs`

## Legacy or transitional modules

- `equilibrium_legacy_backend.rs`
- `easy_equilibrium.rs`
- `NR_Legacy.rs`

## Public entry points currently exposed by the canonical stack

- `equilibrium_calculator.rs`
  - application-facing builder for explicit species, elemental inventory,
    `P,T`/`P,H` point and range requests
- `phase_bridge/phase_equilibrium_workflow.rs`
  - resolved-system `P,T` facade, including fixed declared phases, bounded
    phase control, elemental inventory, and transactional normalization
    recovery
- `ph/equilibrium_ph_workflow.rs`
  - resolved-system `P,H` facade with monolithic, nested, and `Auto` route
    selection; monolithic is attempted first by default
- `equilibrium_ids.rs`
  - typed identifiers for species, elements, phases, and reactions
- `equilibrium_problem.rs`
  - problem, prepared problem, solution, scaling, and initial-guess contracts
- `equilibrium_solver_policy.rs`
  - solver backend selection, attempt outcomes, budgets, and solve reports
- `equilibrium_validation.rs`
  - backend-independent candidate validation and candidate comparison helpers
- `equilibrium_rst_backend.rs`
  - typed RST solve contract, opaque prepared RST problem, and RST adapters
- `phase_control/prepared_phase_control_runner.rs`
  - immutable active-set lifecycle, stability checks, continuation, and
    accepted-only publication
- `equilibrium_workflows.rs`
  - phase-management primitives and retained mutable compatibility helpers
- `equilibrium_log_moles.rs`
  - mutable compatibility host, solver settings, sweep compatibility, and
    pure residual/Jacobian helpers
- `equilibrium_nonlinear.rs`
  - reaction basis, legacy numerical solvers, and temporary nonlinear support

## Test families

### Canonical formulation and problem tests

- `equilibrium_calculator.rs` tests
- `equilibrium_problem_tests.rs`
- `equilibrium_validation_tests.rs`
- `equilibrium_log_moles_tests.rs`
- `equilibrium_workflow_tests.rs`
- `equilibrium_public_api_tests.rs`
- phase-bridge tests and stories
- P,H workflow and range tests

### Solver policy and backend tests

- `equilibrium_solver_policy.rs` tests
- `equilibrium_backend_adapter.rs` tests
- `equilibrium_rst_backend_tests.rs`
- `equilibrium_rst_matrix_tests.rs`

### Low-level numerical support tests

- `equilibrium_nonlinear.rs` tests

### Transitional / historical regression suites

- `easy_equilibrium.rs` internal tests
- `NR_Legacy.rs` internal tests

### Independent and frozen evidence

- `cross_validation/` reaction-extent and pure-phase reference checks
- `frozen_reference/` source-backed thermochemistry and lifecycle stories
- `test_suites/live_data/` offline real-data production stories

### Classification notes

- Canonical regression and solver-policy tests:
  - `equilibrium_problem_tests.rs`
  - `equilibrium_validation_tests.rs`
  - `equilibrium_log_moles_tests.rs`
  - `equilibrium_workflow_tests.rs`
  - `equilibrium_rst_backend_tests.rs`
  - `equilibrium_rst_matrix_tests.rs`
  - `equilibrium_solver_policy.rs` tests
- Low-level numerical characterization:
  - `equilibrium_nonlinear.rs` tests
- Migration and legacy regression coverage:
  - `NR_Legacy.rs` tests
  - `easy_equilibrium.rs` tests
- Debug / experimental studies:
  - none currently retained as separate standalone source files

## Current characterization notes

- `equilibrium_log_moles.rs` is the main orchestration entry point for the
  retained mutable compatibility path, not the production facade.
- `equilibrium_calculator.rs` is the application-facing production facade.
  It resolves the requested species universe and delegates `P,T` and `P,H`
  work to the typed workflows.
- `equilibrium_prepared_runner.rs` owns one prepared fixed-formulation solve;
  `phase_control/prepared_phase_control_runner.rs` owns the immutable bounded
  lifecycle around those solves.
- The production default is an explicit RST-first cascade with the configured
  legacy numerical family retained as fallback. A caller may select a single
  backend or an explicit ordered policy.
- `equilibrium_backend_adapter.rs` is the mechanical boundary between the
  domain layer and concrete numerical backends.
- `equilibrium_rst_backend.rs` owns the RST-specific symbolic bridge and
  should not leak RST internals upward.
- `equilibrium_validation.rs` is the backend-independent acceptance gate and
  should remain the only place where candidate evidence is compared against
  conservation and tolerance thresholds.
- `equilibrium_nonlinear.rs` and `equilibrium_legacy_backend.rs` contain the
  active handwritten Legacy LM/NR/TR fallback family. They are distinct from
  the deprecated standalone `NR_Legacy.rs` compatibility module.
- Legacy and classical modules remain part of the tree for migration and
  characterization, but they are no longer the preferred place for new
  equilibrium behavior.

## Known audit targets for later passes

- remaining public mutation sequences that still need to collapse into one
  validated construction plus one solve call
- remaining `panic!`, `expect`, unchecked indexing, and direct console output
  are concentrated in the legacy/classical modules; canonical `equilibrium_*`
  modules are already much better but still need a final sweep
- measured repeated P,H property evaluation before changing ownership or
  caching; no performance redesign is implied without measurements
- legacy module deletion criteria once the regression matrix is fully covered

## Known difficult or fallback-prone fixtures

- `equilibrium_rst_matrix_tests::high_temperature_n2_dissociation_falls_back_after_rejected_nielsen_candidate`
  is intentionally hard and must remain a fallback-prone regression fixture.
  It exists to prove that the solver can reject one candidate and continue to a
  better backend without mutating published state.
- `equilibrium_rst_matrix_tests::every_selected_rst_backend_accepts_the_o2_dissociation_fixture`
  is a simple acceptance fixture and should stay in the fast regression set.
- Currently retained canonical equilibrium failures are recorded in TODO as
  soon as they are discovered; there is no reason to soften a failing assertion
  into a passing one just to keep the suite green.
