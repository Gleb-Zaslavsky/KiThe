# Chemical Equilibrium Refactoring Plan

## Scope and decisions

This checklist covers `Thermodynamics/ChemEquilibrium` only. Its goal is one
maintainable equilibrium engine with explicit numerical policies, typed failure
reporting, independent validation, and deterministic tests.

## P8 - Independent Pure-Phase Boundary Cross-Validation

The original validation kernel is scoped to synthetic ideal-gas reactions
with one pure condensed candidate at fixed `P,T`. It remains an independent
harness, not a second production phase-control path. Later P8-P10 layers reuse
that kernel with immutable offline real-data fixtures and P,H adapters while
still excluding network NIST, Cantera, and non-ideal activity models.

- [x] Let `PurePhaseBoundaryProblem` accept independently supplied elemental
  compositions and reject a non-conserving supplied reaction instead of
  silently projecting it onto a conservative basis.
- [x] Add an opt-in strict-family structural contract based on an independent
  SVD rank calculation: one full reaction direction and no gas-only reaction
  direction.
- [x] Cover valid structure, non-conservation, excessive full reaction-space
  dimension, residual gas-only chemistry, and explicit opt-in behavior with
  focused synthetic tests.
- [x] Add a controlled analytic-temperature synthetic fixture with a known
  boundary `T*` and an independent `ln(Q)-ln(K)` bisection reference.
- [x] Bridge that fixture to canonical phase-control evidence and compare
  `TPD_candidate` against `R*T*(ln(Q)-ln(K))/nu_candidate` before hysteresis.
- [x] Keep outer-loop topology, appearance/disappearance, hysteresis, cycle,
  and rollback tests separate from pure-phase boundary algebra. They remain in
  `prepared_phase_control_runner.rs` and live lifecycle matrices rather than
  being duplicated by the independent `K_eq` validator.
- [x] Compare finite equilibrium from the independent `K_eq` extent solver
  against canonical final equilibrium for the strict family.
- [x] Prove that the independent finite-equilibrium state is invariant under
  synthetic reaction-coordinate scalings `0.5*nu`, `nu`, and `2*nu`.
- [x] Repeat the finite-equilibrium comparison under stoichiometric scalings
  `0.5*nu`, `nu`, and `2*nu`; equilibrium state must remain invariant.

### P8.1 - Pure-Phase Boundary Lifecycle Evidence

This follow-up keeps the three validation layers explicit: boundary algebra,
fixed-topology Gibbs equilibrium, and the production active-set lifecycle.

- [x] Classify `Boundary` as `HysteresisDependent` rather than implicitly
  inactive; report topology agreement as `Option<bool>`.
- [x] Split cross-validation evidence into thermodynamic, topology, and
  composition agreement fields. The aggregate remains only a convenience.
- [x] Make the one-point TPD test compare TPD and independent driving force
  directly, without borrowing finite K_eq composition as fake evidence.
- [x] Add an analytic three-temperature sweep proving both sign and magnitude
  agreement of independent `ln(Q)-ln(K)` and canonical TPD.
- [x] Extend scaling tests with reaction-coordinate quantities, inverse extent
  scaling, and invariant per-candidate driving force.
- [x] Add a production `PhaseManager` hysteresis matrix for activate, inactive
  band retention, active hold, and safe deactivation.
- [x] Add an end-to-end `PreparedPhaseControlRunner` synthetic story:
  inactive candidate -> negative evaluated TPD -> activation record -> restart
  -> accepted two-phase state matching independent finite K_eq equilibrium.
- [x] Add the symmetric end-to-end stable-inactive story: positive TPD, no
  activation record, and no finite independent two-phase solution.
- [x] Assert lifecycle provenance on the activation story: the transition must
  retain the evaluated negative TPD and the TPD-derived restart seed, rather
  than merely widening a numerical phase mask.
- [x] Publish post-seed phase totals in activation records, so audit evidence
  describes the physical restart state rather than the pre-activation trace.

### P8.2 - Synthetic Boundary Regression Completion

Keep the independent validator, canonical thermodynamic evidence, and
production lifecycle as three separately named test layers. Do not expand this
synthetic work to real databases, Boudouard chemistry, P,H, or multi-component
candidate phases.

- [x] Move canonical TPD and fixed-topology Gibbs comparisons out of
  `phase_boundary_validation.rs`; that module must test only the independent
  validator, while `phase_boundary_cross_validation_tests.rs` owns cross-system
  evidence.
- [x] Add test-local canonical-TPD bisection and prove that its mathematical
  root agrees with both the independent `ln(Q)-ln(K)` root and known `T*`.
- [x] Lock exact `PhaseManager` hysteresis inequalities: equality with
  `dg_create` must not activate and equality with `dg_keep` must not deactivate.
- [x] Compare lifecycle TPD values numerically, not only by sign, against the
  independent per-candidate driving force in both appearance and stable-absence
  stories.
- [x] Extend scaling evidence with invariant gas mole fractions.
- [x] Add deliberately broken canonical evidence fixtures proving that the
  cross-validation report independently localizes thermodynamic, topology, and
  composition disagreement.
- [x] Add a second, non-collinear strict synthetic phase-forming family with
  three gas species and `nu_candidate != 1`, and repeat the key TPD, lifecycle,
  and finite-equilibrium cross-validation stories.

### P8.3 - Pure-Phase Lifecycle Completion

Close the remaining synthetic lifecycle evidence without modifying production
hysteresis semantics, diagnostic-report structure, or the public API. This
stage remains limited to one-component pure candidate phases and deliberately
excludes real databases, non-ideal activity models, P,H, and candidate
solutions.

- [x] Make the test-local canonical-TPD boundary bisection fail explicitly on
  exhausted iterations instead of returning an unverified bracket midpoint.
- [x] Assert that each controlled appearance fixture makes exactly one
  `Activate` transition, while stable-inactive fixtures make none.
- [x] Add an end-to-end `PreparedPhaseControlRunner` disappearance story.
  For a one-component pure phase with no positive interior solution, the
  canonical route is active fixed-set failure -> validated boundary recovery
  on the gas-only set -> positive TPD -> `Deactivate` -> accepted restart.
  The test locks that physical route, conservation, trace reseeding, and the
  exactly-one-transition contract.
- [x] Assess a history-dependent end-to-end hysteresis story using the current
  continuation API. Do not add a second continuation or orchestration API just
  to make this synthetic test possible. Assessment correction: the existing
  runner already exposes crate-local numeric retargeting and accepted
  `PhaseSet` continuation, and the typed temperature-range facade uses both.
  The actual history-aware regression belongs to P8.7 below.

### P8.4 - Cross-Validation Identity and Evidence Completeness

The comparator must prove that the independent result, canonical result, and
declared boundary problem describe the same physical case. Vector lengths and
the absence of a failed check are not sufficient evidence.

- [x] Introduce a deterministic boundary-case identity/fingerprint containing
  conditions, ordered gas-component identities, candidate identity,
  stoichiometry, and the gas-only boundary inventory. Either retain this
  identity in `PurePhaseValidationResult` or make the high-level comparator run
  the independent validator directly from the supplied problem.
- [x] Add ordered component identities to `CanonicalPurePhaseEvidence`.
  Composition comparison must align by identity or reject a mismatched layout;
  it must not assume that two unnamed `Vec<f64>` values use the same order.
- [x] Keep the current low-level comparison helper crate-private/test-only if
  needed, and expose one canonical cross-validation entry point that cannot be
  called with an independent result from a different problem.
- [x] Separate agreement from evidence completeness. Add an explicit coverage
  or status contract such as `Complete`, `ConsistentButPartial`,
  `InsufficientEvidence`, and `Disagreed`; absence of applicable checks must
  never produce a fully validated result.
- [x] Preserve the three independent diagnostic axes: thermodynamic,
  topology, and composition. The aggregate status may summarize them but must
  not erase `None`/not-applicable evidence.
- [x] Add contract tests:
  - [x] `cross_validation_rejects_independent_result_from_another_problem`;
  - [x] `cross_validation_rejects_independent_result_at_another_temperature`;
  - [x] `canonical_evidence_matches_components_by_identity_not_vector_position`;
  - [x] `consistent_species_permutation_preserves_cross_validation`;
  - [x] `zero_applicable_checks_cannot_produce_validated_status`;
  - [x] `partial_evidence_is_consistent_but_not_complete`.

### P8.5 - Production Evidence Adapters

Connect the independent mathematics to evidence that is actually published by
the production outer loop. Avoid test-only manual reconstruction of canonical
TPD and final compositions wherever a typed production report already owns the
same facts.

- [x] Add a typed adapter for a stable inactive candidate using the final
  accepted `PhaseStabilityReport` plus the accepted component layout.
- [x] Add a typed adapter for phase appearance using the activation
  `PhaseTransitionRecord` as boundary TPD evidence and the final accepted
  solution as composition/topology evidence. A final active-phase TPD report
  must not be substituted for the pre-activation gas-only boundary TPD.
- [x] Make both adapters verify candidate phase identity, one-component pure
  topology, component ordering, finite TPD, and finite feasible mole values.
- [x] Route the existing activation, stable-inactive, and disappearance
  stories through the high-level comparator after their lifecycle assertions.
  - [x] Activation and stable-inactive now use the immutable production adapter
    plus independent local NASA `H2O(g) <=> H2O(l)` validation. Activation
    compares TPD at the pre-activation gas boundary reconstructed from the
    transition restart seed, while composition remains the final accepted
    two-phase state.
  - [x] Disappearance now has a symmetric immutable adapter. It accepts only a
    recorded `BoundaryUnstableActivePhase` transition from active to inactive,
    uses the transition's reduced-boundary TPD/restart state, and never infers
    physical absence from final trace moles. The local NASA water/liquid story
    passes this evidence through the same independent `ln(Q)-ln(K)` comparator.
    The existing synthetic story remains the focused lower-level runner proof.
- [x] Add adapter tests:
  - [x] `stable_inactive_outcome_builds_complete_canonical_evidence`;
  - [x] `activation_transition_builds_boundary_and_final_composition_evidence`;
  - [x] `adapter_rejects_wrong_candidate_phase`;
  - [x] `adapter_rejects_multicomponent_candidate` (the gas assemblage cannot
    masquerade as a pure condensed candidate);
  - [x] `adapter_rejects_missing_or_nonfinite_tpd` (wrong lifecycle route is
    rejected before evidence publication);
  - [x] `adapter_rejects_component_layout_mismatch` through immutable
    solution/layout bounds and named phase lookup.
  - [x] Disappearance rejects both a still-active candidate and an inactive
    candidate without a matching deactivation transition.

### P8.6 - Metamorphic Thermodynamic Matrix

Exercise transformations whose expected physical effect is known analytically.
These tests provide more independent information than adding another reaction
family with the same pressure, amount scale, and component ordering.

- [x] Add combined absolute-plus-relative composition tolerances. Absolute
  tolerances alone cannot compare otherwise equivalent systems spanning trace
  inventories through very large mole counts. The independent scalar extent
  bisection now operates in normalized extent coordinates, and the scaling
  matrix uses scale-aware absolute-plus-relative assertions.
- [x] Verify ideal-gas pressure dependence:
  `delta ln(Q) = delta_nu_gas * ln(P2/P1)`, including matching canonical TPD.
- [x] Verify invariance when pressure and reference pressure are multiplied by
  the same positive factor, preserving `P/P0`.
- [x] Add an inert gas with zero reaction stoichiometry and an independent
  elemental row. Its Gibbs contribution is zero, but dilution must alter
  reacting-species mole fractions, `ln(Q)`, and TPD consistently.
- [x] Scale the complete initial inventory by
  `1e-9, 1e-3, 1, 1e3, 1e9`. Boundary prediction and mole fractions must be
  invariant; equilibrium extent and every physical phase amount must scale
  linearly.
- [x] Permute gas species, elemental-composition rows, canonical phase
  component indices, and reported output consistently. Boundary root, TPD,
  topology, and identity-aligned composition must remain invariant.
- [x] Lock the independent classification tolerance with residuals at
  `-1.01*tol`, `-0.99*tol`, `0`, `0.99*tol`, and `1.01*tol`. Keep this
  tolerance separate from production `dg_create`/`dg_keep` hysteresis.
- [ ] Add tests:
  - [x] `pressure_dependence_matches_delta_nu_gas_ln_pressure`;
  - [x] `joint_pressure_reference_scaling_is_invariant`;
  - [x] `inert_dilution_changes_q_and_matches_canonical_tpd`;
  - [x] `global_inventory_scaling_preserves_boundary_and_scales_equilibrium`;
  - [x] `species_permutation_preserves_boundary_root_and_final_equilibrium`;
  - [x] `boundary_classification_respects_independent_tolerance_edges`.

### P8.7 - History-Aware End-to-End Lifecycle

Use the existing `PreparedPhaseControlRunner::retarget_numeric`, continuation
phase-set support, and typed temperature-range path. Do not emulate history by
calling `PhaseManager` manually.

- [x] Construct a controlled temperature-dependent synthetic boundary and
  solve a point where the pure candidate is unambiguously active. Retarget the
  same prepared runner to a nearby point inside the hysteresis band while
  carrying only the previously accepted seed and `PhaseSet`.
- [x] Solve the same second point from a fresh inactive history. Prove that
  discrete topology follows the declared hysteresis policy while continuous
  TPD evidence remains identical.
- [x] Run ascending and descending grids across the analytic boundary. Record
  transition temperatures, reasons, counts, continuation provenance, and
  final complementarity evidence.
- [x] Characterize the log-moles boundary explicitly. A one-component pure
  phase with no positive interior root may require validated boundary recovery
  rather than an ordinary active-set solve. If an active in-band history cannot
  be represented, expose a typed lifecycle limitation or fix the canonical
  boundary representation; do not loosen nonlinear acceptance tolerances.
  - Fixed: boundary recovery now examines the accepted continuation seed,
    rather than construction-time `initial_moles`. A phase created at an
    earlier temperature can therefore be considered for later disappearance.
- [x] Require every continuation seed to come from an accepted point. Inject a
  failed intermediate solve and prove transactional restoration of both seed
  and phase set before the next attempt.
  - [x] Owner-level `PreparedPhaseControlRunner` tests inject a failure after
    continuation consumption and assert restoration of seed, phase set, and
    streamed rollback diagnostics. This is deliberately not duplicated in the
    independent cross-validation suite.
- [ ] Add tests:
  - [x] `same_in_band_point_uses_previous_active_set_history`;
  - [x] `ascending_and_descending_boundary_sweeps_expose_hysteresis`;
  - [x] `failed_boundary_point_does_not_poison_continuation` (covered by
    `failed_lifecycle_attempt_restores_accepted_continuation_state`);
  - `boundary_history_preserves_element_inventory_and_trace_contract`.

### P8.8 - Independent Root and Error Matrix

Lock all typed failure modes of the independent implementation. These tests
belong in a separate `phase_boundary_validation_tests.rs` module so production
code and the cross-system suite remain readable.

- [x] Cover temperature-root endpoint acceptance, invalid/unordered brackets,
  missing sign changes, problem-factory temperature mismatch, non-finite
  factory output, and explicit iteration-budget exhaustion.
  - [x] Endpoint acceptance, unordered/missing-sign/mismatched-factory
    rejection, and explicit `MaxIterations` are covered in the dedicated
    `phase_boundary_validation_tests` module.
  - [x] Non-finite factory Gibbs output preserves the candidate component and
    requested-temperature context.
- [x] Cover finite-extent roots close to zero and close to a gas-species
  positivity boundary without evaluating a non-positive activity.
- [x] Add a favorable-at-zero case with no finite interior two-phase root and
  require typed `ValidationNotApplicable`, documenting that exact gas-species
  disappearance is outside this validator's current scope.
- [ ] Cover invalid structural, boundary, cross-validation, and scalar-solver
  tolerances, including NaN and infinity.
- [x] Cover non-finite gas and candidate Gibbs closures with typed errors that
  preserve the failing temperature and component context.
  - [x] Boundary and scalar-solver invalid settings are covered; structural and
    cross-validation tolerance matrices remain separate work.
  - [x] Candidate Gibbs now uses `InvalidDG0` with deterministic candidate
    component index `gas_species.len()`, matching gas-closure failure typing.
- [ ] Add tests:
  - [x] `temperature_root_accepts_each_bracket_endpoint`;
  - [x] `temperature_root_rejects_unbracketed_and_mismatched_factory_cases`;
  - [x] `temperature_root_fails_explicitly_after_iteration_budget`;
  - [x] `finite_extent_root_is_robust_near_each_feasibility_boundary`;
  - [x] `favorable_boundary_without_interior_root_is_not_applicable`;
  - [x] `invalid_settings_and_nonfinite_gibbs_return_typed_errors`.

### P8.9 - Validation Coverage Inventory

Use the following levels consistently when describing evidence:

- **I1 - independent thermodynamics:** direct `ln(Q)-ln(K)` or scalar extent
  mathematics that does not reuse canonical Gibbs residual equations;
- **I2 - independent numerical route:** a separate scalar bisection/root solve
  compared with canonical TPD or fixed-topology minimization;
- **I3 - production lifecycle:** the real prepared active-set outer loop,
  continuation, hysteresis, rollback, and immutable reports;
- **I4 - offline real data:** local repository resolution and real
  thermochemical closures, with no network and no JSON mutation.

Passing at I3 or I4 does not imply I1 independence. Conversely, an I1
synthetic fixture does not prove the repository-to-solver production path.

| Physical / numerical scenario | I1-I2 independent evidence | I3 production lifecycle | I4 offline real data | Remaining gap |
|---|---|---|---|---|
| Pure phase should appear | Two synthetic reaction families; canonical TPD magnitude/sign; scalar finite extent | Prepared activation transition and restart composition | Shared local NASA water/liquid, water/ice, and Boudouard/graphite fixtures pass the complete comparator | Current I1+I3+I4 coverage includes two chemical families and is sufficient |
| Pure phase remains absent | Positive independent boundary residual and canonical TPD | Stable-inactive prepared solve | High-temperature water/liquid and hot carbon stories | Current coverage is sufficient; more same-family tests would duplicate it |
| Pure phase disappears | Independent positive TPD exists at the reduced boundary | Synthetic boundary recovery, continuation-aware active-to-inactive sweep, immutable disappearance adapter | Local NASA water/liquid disappearance passes the high-level comparator | Current I1+I3+I4 coverage is sufficient |
| Exact boundary | Independent endpoint/root bisection and tolerance-edge matrix | History-dependent topology is explicitly not overclassified | No stable real-data exact-root fixture | Optional real boundary fixture; do not hard-code a database-dependent exact temperature |
| Hysteresis / path dependence | Independent continuous TPD is fixed while topology differs | In-band previous-active vs fresh-inactive; ascending/descending sweep | Real marginal liquid retention | Add no duplicate test unless P,H continuation exhibits a distinct contract |
| Pressure and reference pressure | Analytical `delta_nu*ln(P2/P1)` plus canonical TPD | Same canonical activity implementation | Covered indirectly by real P,T cases | Sufficient for ideal-gas physics |
| Inert dilution | Independent quotient change plus canonical TPD | No separate lifecycle story | No dedicated real-data inert fixture | Optional; current I1-I2 evidence already isolates the equation |
| Inventory / reaction-coordinate scaling | `1e-9..1e9`, stoichiometric scaling, near-zero and near-positivity roots | Conservation and trace contracts in phase-control tests | Large real-data scaling matrices elsewhere | Generated corpus still needs explicit rank/comparator status per case |
| Species / element permutation | Species, element-row/column, identity-aligned evidence permutation | Canonical layout ordering tested elsewhere | Repository provenance preserves named components | Sufficient; avoid additional positional-vector tests |
| Failure and rollback | Typed scalar/root/settings/Gibbs errors | Injected continuation rollback, budget and cycle evidence | Real water budget rollback | Sufficient for this stage |
| P,H phase lifecycle | Local water stable absence/disappearance compare independent `ln(Q)-ln(K)` with canonical TPD; water, ice, and Boudouard points add independent scalar P,H roots | Monolithic/nested P,H lifecycle, accepted-only continuation, rollback, and hysteresis are covered | Shared local NASA water/liquid, water/ice, and graphite fixture families are resolved read-only | Core pure-phase I1-I4 matrix is complete; non-NASA payload remains a separate format-diversity gap |
| Non-ideal or multi-component candidate phase | Not applicable to current pure-phase scalar contract | Unsupported models rejected explicitly | None | Future physics, not a missing test for the current validator |
| Thermochemistry format diversity | Formula-agnostic closure contract | Same solver API for all closures | NASA gas/condensed only for stable local fixtures | Local non-NASA fixture remains a fundamental evidence gap |

Prioritized gaps after the inventory:

1. [x] Complete immutable disappearance evidence and high-level comparison.
2. [ ] Add rank/conservation/comparator-completeness assertions to each fixed
   generated fixture; do not grow the corpus until those axes are complete.
3. [ ] Add a stable local non-NASA fixture when one exists in the repository.
4. [ ] Keep exact real P,T boundary temperature deferred, but develop the
   independent synthetic P,H boundary validator under P9 before adding another
   real-data fixture.

Evidence ownership (use this map before adding another test):

- `phase_boundary_validation_tests.rs` owns scalar-root contracts, malformed
  inputs, feasibility-boundary behavior, and typed independent-validator
  failures (I1-I2).
- `phase_boundary_cross_validation_tests.rs` owns independent K_eq/extent vs
  canonical TPD comparisons, metamorphic invariants, deterministic generated
  fixtures, and synthetic prepared-runner lifecycle stories (I1-I3).
- `prepared_phase_control_runner.rs` owns transactional continuation rollback
  and streamed rejection/rollback diagnostics (I3). Cross-validation tests
  should not duplicate these unless they also compare independent physics.
- `pure_phase_pt_live_data_tests.rs` owns shared-fixture P,T appearance,
  stable-absence, and disappearance comparisons against independent
  `ln(Q)-ln(K)` boundary problems (I1-I4). It also owns the read-only JSON
  contract for those stories.
- `equilibrium_multiphase_story_tests.rs` owns general immutable fixed-P,T
  result/facade stories. It consumes shared fixtures where phase-pair data are
  needed but no longer owns a second hand-built boundary-validation fixture.
- `equilibrium_live_data_tests.rs` owns repository-to-solver water/ice,
  water/liquid, carbon, hysteresis, budget rollback, and release evidence
  using local thermochemical libraries (I4).

Before accepting a new validation test, require it to add at least one missing
cell in the matrix, raise an existing scenario to a stronger independence
level, or reproduce a distinct regression. A new temperature, tolerance, or
backend alone is not additional physical evidence.

### P8.10 - Deterministic Generated and Offline Evidence

Add breadth only after identity, completeness, adapter, and metamorphic
contracts are stable.

- [x] Build a deterministic generated-fixture matrix (fixed seed, no flaky
  randomness) for strict one-reaction families. Vary stoichiometric scale,
  inventory scale, pressure ratio, species permutation, and boundary sign.
  - [x] The initial fixed corpus carries a stable fixture id and exercises
    boundary sign, stoichiometric scale, inventory scale, and pressure ratio;
    the species/element permutation leg remains the named P8.6 fixture.
- [ ] For every generated case, check independent conservation/rank evidence,
  finite-extent acceptance when applicable, canonical TPD magnitude/sign, and
  comparator coverage status.
- [x] Keep generated failures reproducible by printing the fixture seed and
  all physical inputs in the assertion message.
- [x] Add one or two stable offline real-data pure-phase stories only after the
  synthetic suite is complete. Reuse immutable local thermochemistry and prove
  that tests do not modify JSON libraries; do not make network NIST access part
  of this suite. Shared water/liquid P,T now covers appearance, stable absence,
  and disappearance; shared water/ice and Boudouard fixtures add independent
  solid-phase appearance stories. Every story snapshots the local JSON files.
- [x] Do not duplicate general backend cascade, cycle, cancellation, or
  rollback matrices already owned by their production workflow modules. Add a
  cross-validation story only when it contributes independent K_eq evidence.
  The former hand-built water adapter story was removed after the shared P,T
  matrix superseded it.

### P8.11 - Boundary Module Engineering Hygiene

- [ ] Repair mojibake in thermodynamic formulas and documentation (`Delta G`,
  `Sigma`, `nu`, `xi`, degree superscripts) using the repository's established
  UTF-8 encoding; add a source-text audit preventing known corrupted sequences
  from returning.
- [ ] Move the independent validator unit tests out of
  `phase_boundary_validation.rs` into
  `phase_boundary_validation_tests.rs`. Keep cross-system tests in
  `phase_boundary_cross_validation_tests.rs`.
- [ ] If the production module remains difficult to navigate after test
  extraction, split only cohesive responsibilities: problem/structure,
  independent roots, evidence comparison, and report formatting. Do not create
  another equilibrium orchestration facade.
- [ ] Use one canonical high-precision molar gas constant within
  `ChemEquilibrium`. Algorithmic independence must come from different
  equations and solvers, not from duplicated physical constants that can
  silently drift.

- [x] Treat the logarithmic `Chem_eq_K_eq*` formulation as the candidate
  source of truth for equilibrium calculations.
- [x] Treat the legacy classical equilibrium stack as a migration/reference path, not as a
  second production architecture where new equilibrium logic should be added.
- [x] Keep hand-written nonlinear solvers only as explicit temporary fallback
  backends. They must not silently define the public equilibrium contract.
- [x] Use the reaction-equilibrium-constant approach as an independent
  validation solver for systems within its documented applicability range.
- [x] Defer the final `phase_*` bridge and GUI until the core problem, result,
  validation, and solver-policy contracts are stable.
- [x] Replace the misleading `Chem_eq_K_eq*` module family with names that
  reflect its actual responsibilities: log-moles formulation, workflows, and
  nonlinear support.

## P9 - Independent Pure-Phase `P,H` Cross-Validation

### Scope and evidence boundary

Build a second mathematical route for one ideal-gas phase plus one pure,
one-component condensed candidate at fixed pressure and total enthalpy. This
is a validation subsystem, not a replacement P,H engine and not a second
phase-control orchestration stack.

The evidence ladder is intentionally staged:

- **I1 - independent thermodynamics:** direct `ln(Q)-ln(K)` and additive
  enthalpy equations in reaction extent/temperature coordinates;
- **I2 - independent numerical route:** nested scalar extent root inside a
  safeguarded scalar temperature root, without canonical P,H residuals,
  Jacobians, backend policy, or prepared P,T workflow;
- **I3 - production lifecycle:** only after I1-I2, compare immutable
  production P,H reports, active-set transitions, and final states;
- **I4 - offline real data:** only after the synthetic I1-I3 suite is stable;
  local resolved records only, with no network or JSON mutation.

Do not describe I1-I2 as independent thermochemical data: initially they may
use the same synthetic functions as production comparison fixtures. Their
independence is the equations and numerical route.

### P9.1 - Keep the P,H validator separate and typed

- [x] Create `pure_phase_ph_validation.rs`; do not turn the fixed-P,T
  `PurePhaseBoundaryProblem` into a broad P,T/P,H union type.
- [x] Define a typed `PurePhasePhProblem` for one strict phase-forming
  reaction with:
  - ordered gas identities, initial gas moles, gas stoichiometry, and one
    positive candidate stoichiometry;
  - physical initial candidate amount, candidate identity, pressure, and
    reference pressure;
  - target total enthalpy and a finite ordered temperature bracket;
  - fallible standard Gibbs and molar enthalpy capabilities. Reuse
    `MolarThermoFunction` where its typed `Result` contract fits instead of
    creating a parallel closure alias merely for this validator;
  - optional independent elemental composition for structural checks.
- [x] Preserve the strict independent-family contract whenever element data is
  supplied: conservative reaction, full reaction-space rank one, and gas-only
  reaction-space rank zero. Reuse only small independent SVD/rank helpers from
  P8, never the canonical reaction-basis builder.
- [x] Keep the scalar bisection local because the P,T extent solver is
  intrinsically tied to a zero-candidate boundary state, while P,H needs a
  general signed physical extent interval. The only shared code is the small
  independent SVD/rank helper; no broad P,T problem type leaks into P,H.
- [x] Return typed `ReactionExtentError` variants for invalid dimensions,
  non-finite functions, infeasible extents, invalid brackets, and exhausted
  budgets. The validator must not fall back silently to `solve_resolved_ph`.

### P9.2 - Independent nested scalar route

- [x] At fixed temperature solve the inner chemical equation
  `ln(Q(xi))-ln(K(T)) = 0` on the physically feasible extent interval.
  Keep the phase model explicit: ideal-gas activities and unit pure-condensed
  activity only.
- [x] Define `H_eq(T) = H(xi_eq(T), T)` from the additive physical amounts and
  independent molar enthalpy functions. Do not call canonical P,H residual,
  Jacobian, enthalpy model, or P,T solve objects from this path.
- [x] Solve `H_eq(T)-H_target = 0` with a bracketed/safeguarded scalar method.
  Require a verified sign bracket, finite inner state/enthalpy, a positive
  finite temperature, and explicit non-convergence failure; never accept an
  unverified midpoint.
- [x] Publish a compact `PurePhasePhEquilibriumResult`: accepted temperature,
  extent, physical gas/candidate moles, chemical log residual, enthalpy
  residual, outer iterations, total inner solves, and total inner iterations.
  Do not expose a sprawling per-iteration diagnostic API in the first pass.

### P9.3 - Deterministic synthetic I1-I2 corpus

- [x] Build a named analytic fixture from predeclared `T_star` and `xi_star`,
  then derive `H_target = H(xi_star, T_star)`. The expected point must follow
  from fixture construction, never from the validator under test.
- [x] Add the base I1-I2 recovery test: recover `T_star`, `xi_star`, physical
  gas/candidate moles, conservation, `abs(lnQ-lnK)`, and enthalpy residual.
- [x] Add reaction-coordinate scaling cases (`0.5*nu`, `nu`, `2*nu`). Check
  inverse extent scaling while temperature, physical amounts, mole fractions,
  total enthalpy, and per-candidate driving force stay invariant.
- [x] Add monotone target-enthalpy perturbations around the exact fixture;
  assert the fixture-specific temperature ordering instead of claiming a
  universal sign for `dT/dH`.
- [x] Add the typed failure matrix: reversed/non-positive/non-finite
  temperature brackets, unbracketed enthalpy target, unavailable inner extent
  root, non-finite Gibbs/enthalpy callback, infeasible extent, and exhausted
  inner/outer iteration budgets.
- [x] Label every test comment with its evidence layer (`I1`, `I2`, later
  `I3`/`I4`) so a passing synthetic root cannot be misread as production
  lifecycle coverage.

### P9.4 - Fixed-topology P,H comparison

- [x] After the independent synthetic route is green, build the same physical
  fixture through the canonical fixed-declared-phase P,H path only. Do not use
  bounded phase control in this comparison.
- [x] Introduce `PurePhasePhCrossValidationReport` rather than forcing P,H
  semantics into the P,T report. It must retain temperature, thermodynamic,
  topology, composition, and enthalpy agreement as `Option<bool>` where
  `None` means genuinely not applicable, never implicit success.
- [x] Compare accepted temperature, physical gas/candidate moles, gas mole
  fractions, additive total enthalpy, conservation, and residual contracts.
  The comparator must keep independent and canonical identity/order checks
  explicit.
- [x] Add deliberately mismatched synthetic evidence to prove that the report
  localizes temperature, thermodynamic, composition, and enthalpy disagreement
  independently.

### P9.5 - Production P,H lifecycle only after fixed topology

- [x] Add a clearly favorable inactive-to-active pure-phase story, far outside
  the hysteresis band. Require real TPD/transition evidence; a candidate in a
  numerical recovery mask is not physical activation.
- [x] Add the stable-inactive counterpart and route its final stability
  evidence through the P,H comparator.
  - The I3 fixture now compares the dimensionless I1 `ln(Q)-ln(K)` driving
    force against the canonical pure-phase TPD through `TPD = R*T*residual`;
    both must be positive before an inactive phase can remain absent.
- [x] Add active-to-inactive recovery. Accept the physically meaningful
  `BoundaryUnstableActivePhase` route when no positive interior state exists;
  do not force an artificial `n_phase < phase_eps` mechanism.
- [x] Only then add an accepted-target continuation sweep and compare every
  accepted P,H point against the independent nested validator. Failed trial
  coordinates must never become continuation seeds.
- [x] Keep an in-band hysteresis scenario deferred until the single-point
  appearance/disappearance contracts are established.
  - The I3 history test now establishes an active phase at a neighboring
    favorable P,H point, retargets to an in-band TPD, and contrasts that
    continuation with a fresh gas-only solve of the identical target. The
    independent scalar I2 route is intentionally not asserted for this
    near-boundary point because its fixed-topology contract requires a
    positive inner extent throughout the complete temperature bracket.

### P9.6 - Offline real-data evidence

The fixed-topology local-water I4 bridge now proves that resolved real closures
can be shared by canonical and independent routes without network access or
library mutation. Phase-controlled water/ice stories prove production
lifecycle behavior, but do not automatically become I1/I2 evidence: the
independent scalar route must admit the same interior topology. Keep that
distinction explicit in reports and test comments.

- [x] Add a pinned P9 evidence inventory before growing the corpus. A passing
  production test sharing the same local database is I3/I4 evidence, not
  automatically independent I1/I2 evidence.

  | Scenario | I1 | I2 | I3 | I4 | Owner / status |
  |---|---:|---:|---:|---:|---|
  | `H2O(g)/O2(g) -> H2O(l)`, 350 K fixed topology | yes | yes | fixed topology | yes | `pure_phase_ph_live_data_tests::p9_i4_local_water_fixed_topology_scalar_route_matches_canonical_ph` |
  | `H2O(g)/O2(g)`, 550 K, liquid initially absent | yes | n/a | yes | yes | `p9_i1_i3_local_hot_water_ph_keeps_liquid_inactive_and_matches_boundary_tpd` |
  | `H2O(g)/O2(g)`, 550 K, liquid initially active | yes | boundary-only | yes | yes | `p9_i1_i3_i4_local_hot_water_ph_deactivates_initial_liquid_transactionally` |
  | `2 CO(g) <=> CO2(g) + C(gr)`, 700 K fixed topology | yes | yes | fixed topology | yes | `p9_i1_i2_i4_boudouard_independent_reference_matches_fixed_topology_ph` |
  | `H2O(g) -> H2O(s)`, 250 K phase control | yes | yes | yes | yes | `p9_i4_local_ice_phase_control_activates_solid_and_matches_independent_ph` |
- [x] Add a pinned local `H2O(g)/O2(g)` plus `H2O(l)` preflight fixture,
  subject to the actual overlap of locally resolved native intervals. It
  resolves through the ordinary immutable repository/capability path with
  explicit `NASA_gas`/`NASA_cond` provenance, NIST disabled, and byte-for-byte
  library immutability before and after the solve.
  - The initial I4 point is intentionally bracketed within `345..355 K`.
    The same local fixture has no liquid-containing inner root at `375 K`, so
    a broad bracket would make the test physically inapplicable rather than
    more robust.
- [x] For that fixture, add an I4 fixed-topology point safely inside the
  liquid-containing region. The test feeds resolved local `G_i(T)` and
  `H_i(T)` capabilities to the independent I1/I2 scalar validator and compares
  its accepted temperature, gas composition, liquid moles, additive enthalpy,
  chemical residual, and elemental conservation with canonical
  fixed-topology `P,H`.
- [x] Add an I3/I4 real appearance story only where the independent scalar
  route has a genuine interior enthalpy bracket for the same lifecycle target.
  - The first local water/liquid candidate is not suitable: its P,H outer
    bracket rejects a topology jump rather than manufacturing a liquid branch.
  - The real water/ice 250 K lifecycle now supplies complete I1-I4 evidence.
    Its former missing bracket was a fixture defect: canonical target enthalpy
    represented `0.5 mol` total water, while the independent problem silently
    added `0.1 mol` initial ice and therefore solved a `0.6 mol` inventory.
    Using one gas-only scenario for P,T target construction, independent P,H,
    and canonical P,H restores the honest bracket without relaxing tolerances.
- [x] Add the I1/I3/I4 hot-water stable-absence story. At the recovered
  gas-only boundary, compare positive independent `ln(Q)-ln(K)` with positive
  canonical TPD using `TPD = R*T*(ln(Q)-ln(K))/nu_candidate`; do not require a
  finite liquid-containing I2 root on the inactive side. The shared local
  `H2O(g)/O2(g)/H2O(l)` test derives its target enthalpy at 550 K, keeps liquid
  inactive, and compares final TPD evidence against independent boundary
  thermodynamics.
- [x] Add the I3/I4 water disappearance story from positive liquid inventory.
  Accept either documented physical route: `VanishingUnstableActivePhase` for
  a finite interior phase driven below the destruction threshold, or
  `BoundaryUnstableActivePhase` when no positive fixed-topology root remains.
  The final result must publish one deactivation, valid P,H acceptance, and
  boundary driving-force agreement when boundary recovery was used. The local
  550 K water story takes `BoundaryUnstableActivePhase`; independent boundary
  thermodynamics is deliberately evaluated at the recorded reduced-boundary
  restart composition, not at the original two-phase inventory.
- [ ] Only after isolated water points are stable, add a short local enthalpy
  sweep through appearance/retention and reverse disappearance. Compare every
  state admitting a positive independent I2 root with that route; use the
  boundary-driving-force comparison for gas-only endpoints. Do not certify a
  database-dependent exact transition temperature.
- [ ] Add real history-dependent hysteresis only if an observed local point
  lies inside the configured band: previously accepted liquid-active history
  must retain it while a fresh inactive history must not create it. The
  independent route supplies continuous driving force only; it does not choose
  a topology inside an artificial hysteresis band.
- [x] Every current P9 I4 test proves no network access and byte-for-byte JSON
  library immutability before and after execution. Do not add online NIST,
  live Cantera/CEA, non-ideal, or multicomponent-solution physics to this
  validation stage.
- [x] Add a second local family, `2 CO <=> CO2 + C(graphite)`, after water is
  stable. Its 700 K fixed-topology P,H story matches an independent scalar
  solution on temperature, composition, enthalpy, chemical residual, and the
  canonical solver's accepted element-balance scale. It adds a different
  gas-phase reaction geometry without pretending to prove format diversity.

### P9.7 - Real Boudouard `P,H` Phase-Control Evidence

Keep the existing frozen JANAF thermochemistry and `P,T` boundary modules
unchanged: they remain I5 reaction/pressure characterization, not an external
`P,H` table. This pass reuses only `RealPurePhaseFamily::BoudouardCarbon` and
the generic independent `PurePhasePhProblem`/validator to complete the missing
local I1-I4 lifecycle evidence at one explicit convention:
`P = p0 = 100000 Pa`.

- [x] Replace the current Boudouard fixed-topology test's canonical-`P,T`
  source point with an independent scalar extent root at a deliberately
  selected interior temperature. Materialize `H_target` solely from the local
  additive `CO/CO2/C(gr)` enthalpy closures at that independent state.
  - [x] Add a small typed independent state materialization API so a test can
    obtain feasible physical moles and additive enthalpy from an accepted
    extent without copying private closure arithmetic or calling canonical
    `P,H` residual code.
  - [x] Recover that target through the independent nested `P,H` root and the
    canonical fixed-topology `P,H` route; compare temperature, extent-scaled
    physical moles, gas fractions, total enthalpy, chemical residual, and
    scale-aware elemental conservation.
- [x] Add a Boudouard I3 appearance story from gas-only inventory only after
  the interior state passes. Its target must come from an independently
  materialized positive-graphite state; require an activation transition,
  negative pre-activation TPD, accepted conservation, and final agreement
  with the independent nested `P,H` state.
  - [x] Review `PhSolveMode::Auto` for a monolithic all-active recovery-probe
    rejection. Deliberately keep the current answer: `Auto` falls back only
    after a classified retryable *numerical* failure; a phase-stability
    rejection is physical evidence and must not be bypassed silently. The I3
    topology-changing Boudouard story therefore requests explicit guarded
    `NestedTemperature`, consistent with the existing real water/ice Auto
    probe-rejection contract.
- [x] Add stable-inactive and active-to-inactive Boudouard `P,H` stories.
  For an inactive final state, compare positive independent `ln(Q)-ln(K)` to
  positive canonical TPD at the accepted gas-only boundary. If no positive
  graphite I2 root exists, treat that as boundary evidence, not a failed
  interior-root test.
- [x] Add a short forward/reverse accepted-target enthalpy sweep only after
  the three isolated lifecycle stories are stable. Reuse only the last
  accepted production state, compare interior points to freshly materialized
  independent I2 problems, and use boundary evidence for gas-only endpoints.
  Keep real in-band hysteresis deferred unless a natural local point occurs.
  - [x] Make the independent scalar P,H validator branch-aware: a positive
    condensed-phase extent may exist only on a strict subinterval of the
    declared common temperature range. Scan only for adjacent valid interior
    roots and bracket `H-H_target` inside that continuous branch; do not turn
    a physical endpoint disappearance into a solver or data failure.
- [x] Add an ignored compact diagnostic table for the selected interior state
  and each lifecycle route: `P`, `p0`, `H_target`, independent/canonical/
  production temperature and moles, phase status/transitions, TPD or chemical
  driving force, and validation errors. Snapshot both local JSON libraries
  and frozen JANAF rows around every Boudouard test; no network NIST, external
  `P,H` reference table, CEA/Cantera, non-ideal activity, or large CHON sweep
  belongs in this stage. The P,H range portion prints both final lifecycle
  transitions and the larger total of nested-temperature trial phase events;
  they are intentionally not conflated.

### P10 - Shared real pure-phase fixture layer for P,T and P,H

The real-data suite must describe a resolved physical family once, then adapt
that same immutable chemistry to independent P,T, independent P,H, and
canonical requests. It must not duplicate water/carbon record selection or
reaction vectors across test modules.

- [x] Create a crate-private real pure-phase fixture module. Its base fixture
  owns resolved phase data, canonical component identities, selected-record
  provenance, explicit reaction stoichiometry, elemental composition, and the
  common thermochemistry interval. It contains no target enthalpy, expected
  phase transition, or other scenario assertion. `real_pure_phase_fixtures`
  now holds `RealPurePhaseFamily`, immutable resolved payloads, and shared
  inventory data used by both P,T and P,H scenarios.
- [x] Resolve and report an offline inventory through ordinary KiThe APIs for
  four deliberately restricted validation families: `H2O(g)/H2O(l)`,
  `H2O(g)/H2O(s)`, `2 CO(g) <=> CO2(g) + C(condensed)`, and
  `CH4(g) <=> 2 H2(g) + C(condensed)`. Pin provenance and record identities;
  do not assume a library-specific phase or polymorph name without checking
  the resolved report. The ignored inventory diagnostic currently finds all
  four locally: water/liquid `[273.15, 600] K`, water/ice `[200, 273.15] K`,
  and both carbon families `[200, 5000] K`.
- [ ] Have each constructor validate: exact requested states, no NIST use,
  immutable JSON libraries, non-empty common temperature range, G/H capability,
  Cp availability when monolithic P,H is requested, elemental conservation,
  `full reaction dimension == 1`, and `gas-only reaction dimension == 0`.
  Incompatible real families must return a typed not-applicable reason rather
  than silently widening their species universe or weakening the scalar test.
- [x] The resolved offline fixture rejects enabled NIST fallback and validates
  common G/H/Cp capability, phase physical state/model/component ordering, and
  non-empty local Thermo provenance from the expected NASA library. Its P,T/P,H
  adapter matrix checks elemental conservation and strict reaction dimensions
  for all four declared families, while byte snapshots prove that resolving
  and materializing every family leaves JSON libraries unchanged.
- [ ] Add narrow adapters from one fixture to `PurePhaseBoundaryProblem`,
  `PurePhasePhProblem`, and existing typed canonical requests. Reuse the
  existing `ResolvedThermochemistry` capabilities directly; never copy
  coefficients into a second closure layer.
  - [x] `PurePhaseBoundaryProblem` and `PurePhasePhProblem` adapters now
    materialize phase-qualified component order, zero-stoichiometry inerts,
    and independent element matrices from one resolved payload.
  - [ ] Add a narrow typed canonical-request adapter once at least two shared
    stories would otherwise repeat `MultiphaseInitialComposition` construction.
- [x] Keep P,T/P,H scenario inputs separate from fixture chemistry. Only add
  `RealPtScenario` or `RealPhScenario` types if repeated tests demonstrate
  actual duplication; do not build a framework merely to name a temperature.
  `RealPurePhaseGasScenario` and `RealPurePhaseInventory` are shared by
  independent and canonical P,T/P,H routes. Inerts and candidate inventory
  remain explicit scenario inputs rather than fixture chemistry.
- [x] Add inventory/capability tests before lifecycle tests. The report must
  say which families support P,T and P,H and why a family is unavailable.
  It is test/diagnostic output only, never default production logging.
  The diagnostic is ignored and prints availability/range/component layout;
  normal tests verify that all inventoried families construct strict P,T/P,H
  independent problems.
- [x] Add a dedicated shared-fixture P,T I1-I4 matrix instead of retaining the
  historical hand-built water pair in the general multiphase story module.
  `pure_phase_pt_live_data_tests` now compares immutable production evidence
  with independently materialized `ln(Q)-ln(K)` problems for water/liquid
  appearance, stable absence and disappearance, water/ice appearance, and
  Boudouard/graphite appearance. The suite also proves byte-for-byte JSON
  immutability.
- [x] After inventory review, select one family with a wide interior interval
  for the first shared real P,T/P,H cross-validation story. Expand to
  appearance/disappearance/hysteresis only one scenario at a time and only
  where the independent scalar contract is applicable. Dedicated P,T stories
  now cover water/liquid appearance, stable absence and disappearance plus
  water/ice and Boudouard/graphite appearance. P,H stories cover fixed topology,
  stable absence, disappearance, ice appearance, and Boudouard chemistry. Every
  route derives its canonical and independent inputs from the same validated
  inventory instead of reconstructing dense vectors by hand.

### P10.1 - Validation Integrity and Scale Contracts

The shared real-data fixtures make it possible to compare P,T and P,H routes
honestly. The comparators themselves must now reject cross-case evidence and
remain meaningful across trace and large inventories.

- [x] Give `PurePhasePhProblem`, `PurePhasePhEquilibriumResult`, and canonical
  P,H evidence one typed case identity containing component identities,
  physical inventory, stoichiometry, pressure/reference pressure, target
  enthalpy, and temperature bracket. The comparator must reject an independent
  or canonical record originating from another case before evaluating numeric
  deltas.
- [x] Publish independent accepted-state elemental conservation evidence in
  the P,H scalar result and compare it explicitly with canonical conservation
  evidence. A canonical balance below tolerance alone is not a two-route
  conservation comparison.
- [x] Replace absolute-only molecule/enthalpy comparison thresholds with
  combined absolute-plus-relative contracts for P,T and P,H comparator paths.
  Keep residual and TPD tolerances dimensionally explicit; do not hide them in
  a generic floating-point epsilon.
- [x] Complete a table-driven invalid-settings matrix for structural, P,T
  boundary/cross-validation, and P,H solver/cross-validation tolerances using
  zero, negative, `NaN`, and infinity inputs. Every rejection must be typed.
- [x] Upgrade the deterministic generated P,T corpus so every case supplies
  independent element/rank evidence and a complete high-level comparator
  result, not only a boundary residual and hand-computed TPD.
- [x] Add symmetric P,H metamorphic checks: full inventory plus target-H
  scaling, joint pressure/reference-pressure scaling, and identity-preserving
  component permutation. Add inert dilution only after its P,H expected
  enthalpy contract is stated explicitly.
- [ ] Add short shared-fixture P,T temperature and P,H enthalpy sweeps after
  the single-point contracts are stable. Reuse accepted continuation only;
  compare each point with the independent scalar route where a positive
  interior root exists and use boundary evidence for gas-only endpoints.
  - [x] P,T: the offline water/liquid `350 -> 450 -> 550 K` story uses the
    production range API, verifies accepted-only continuation and phase-set
    reuse, then independently validates the appearance boundary and final
    stable-inactive boundary without mutating local JSON libraries.
  - [ ] P,H: add the matching shared-fixture target-enthalpy sweep. Each
    accepted range point must be compared with a freshly materialized
    independent P,H scalar problem carrying that point's exact target H.
- [ ] Keep real in-band hysteresis and an exact database transition temperature
  deferred unless an observed local point supplies evidence naturally.

### P11 - I5 frozen external reference evidence

I5 is an independent evidence layer, not an extension of the local repository
path. I4 resolves real thermochemistry from KiThe's offline databases and runs
the production solver. I5 compares KiThe results with small authoritative
external numerical tables frozen in the test tree. Passing either level does
not imply passing the other.

- [x] Add a test-only, read-only frozen-reference loader that does not call the
  NASA/NIST handlers, access the network, mutate local JSON libraries, or enter
  the runtime dependency graph. The chosen first-pass format is a strict JSON
  pair: one provenance document and one typed row document. It reuses existing
  `serde`/`serde_json`; no dependency was added.
- [x] Separate numerical rows from typed provenance. Metadata records stable
  dataset identity, evidence kind, expected data filename, source
  organization/name/version, citation/table/stable identifier, manual
  transcription statement, column meanings, explicit units, and optional
  source precision/uncertainty.
- [x] Make synthetic infrastructure fixtures impossible to confuse with real
  I5 evidence through `FrozenReferenceEvidenceKind`. The initial three-row
  temperature/pressure fixture exists only to exercise parsing and validation.
- [x] Add a generic typed-row contract without introducing a universal
  `HashMap<String, f64>` schema. The first narrow record is
  `TemperaturePressureReference`; future phase-boundary, species-
  thermochemistry, and full-equilibrium tables must define their own row types
  and dimensional schemas.
- [x] Validate non-empty datasets, finite/positive values where required by the
  concrete row type, strictly increasing temperatures, duplicate temperatures,
  required non-empty provenance, unique columns, exact typed schema and units,
  optional uncertainty, dataset-id agreement, and metadata/data filename
  agreement. Parsing is strict and never skips malformed or unknown fields.
- [x] Add negative tests for empty, duplicate, unordered, non-positive,
  malformed, missing-field, unknown-field, non-finite textual/numeric,
  mismatched-unit, dataset-id, and filename cases. Prove byte-for-byte that a
  successful load does not modify either fixture file.
- [x] Preserve row-level provenance context for future assertion failures while
  keeping comparison tolerances outside the loader. Published uncertainty and
  test/model acceptance policy are different contracts.
- [x] Review and freeze the first authoritative external dataset: ten IAPWS
  `H2O(g) <=> H2O(l)` low-pressure saturation points at 275..320 K. The JSON
  provenance pins `IAPWS SR1-86(1992)`, equation (1), the explicit correlation
  table identity, SI units, and the statement that values were calculated once
  from the published correlation before being frozen. Tests never implement or
  call IAPWS at runtime.
- [x] Add the semantic `WaterSaturationPressureReference` row type. Its strict
  validation owns water-specific invariants: temperatures above the triple
  point, positive and strictly increasing temperatures, and strictly
  increasing saturation pressures. These constraints do not leak into generic
  reference rows.
- [x] Add an ignored IAPWS P,T diagnostic that finds one pressure root from
  the independent `ln(Q)-ln(K)` boundary and another from canonical TPD,
  records both against every frozen IAPWS row, prints external and internal
  relative errors, and proves frozen JSON immutability. The I1/I2-versus-TPD
  root has a strict internal contract; KiThe-versus-IAPWS remains diagnostic
  until the observed model discrepancy is reviewed. The diagnostic now uses
  direct `H2O(g)/H2O(l)` with no zero-stoichiometry carrier; the independent
  and canonical TPD roots agree to the strict internal tolerance.
- [x] Move IAPWS row evidence, summary statistics, and table rendering into a
  typed `WaterSaturationComparisonReport`; retain the pressure-root search in
  the physical diagnostic. The report keeps I1/I2 and canonical TPD roots,
  residuals, iterations, external deltas, and strict internal-root validation
  separate, so a printed table is no longer the only evidence carrier.
- [x] Define the first source/model-specific external comparison contract as
  explicit `CharacterizationOnly`. It pins the IAPWS dataset/source identity
  and current NASA-gas/NASA-condensed ideal-phase scope, but deliberately does
  not turn the observed model discrepancy into a generic solver tolerance.
  A later physics review must choose any enforced external acceptance bounds.
- [x] Add `dataset_format_version` to the strict frozen-data schema and record
  the IAPWS correlation-evaluation/rounding provenance plus source precision.
  Missing or zero format versions are rejected; the loader remains read-only
  and deliberately has no refresh or download path.
- [x] Generalize the canonical fixed-P,T bridge to reduce dependent elemental
  constraints deterministically for physically valid low-dimensional systems.
  The direct `H2O(g) + H2O(l)` fixed-declared path now retains complete H/O
  conservation evidence in its bridge report while providing only a stable
  independent H basis to the square log-moles formulation. Real local NASA
  gas/condensed water solves without an O2 carrier, and a boundary regression
  proves that a zero-stoichiometry O2 carrier does not change `p_H2O`.
- [x] Extend the same direct low-dimensional support through bounded TPD
  phase-control. The failure was an active-set/TPD contract bug, not bad NASA
  data and not a backend-specific convergence defect: after liquid activation
  away from saturation, the all-active pure-water equations correctly have no
  interior coexistence root, but the liquid-only recovery could not prove gas
  absence because every ideal-gas phase was unconditionally labelled
  `FixedGasAssemblage`. An inactive gas phase is now evaluated against the
  condensed reference only when no gas assemblage is active. Existing active
  gas assemblages still retain the one-gas-phase policy. The production path
  now performs the validated `gas -> liquid` replacement, conserves inventory,
  and the direct IAPWS TPD boundary no longer needs `O2` or relaxed acceptance.
- [x] Complete the two-sided gas/condensed regression contract exposed by the
  direct-water lifecycle fix. This is internal phase-stability hardening, not
  a new external dataset or a real-fluid model:
  - [x] Keep an active ideal-gas reference as `FixedGasAssemblage` with no gas
    TPD while the inactive liquid remains an evaluated candidate.
  - [x] Evaluate an inactive ideal gas against a condensed-only reference and
    retain a finite TPD instead of `FixedGasAssemblage`; the focused synthetic
    unit contract is present.
  - [x] On real local water data, require positive inactive-gas TPD above the
    saturation boundary and negative inactive-gas TPD below it.
  - [x] Find the canonical gas-side root
    `TPD(liquid | gas reference) = 0` and the independent condensed-side root
    `TPD(gas | liquid reference) = 0` using test-local log-pressure bisection.
    Require both roots and the independent I1/I2 root to agree under one strict
    internal tolerance. Keep this search out of the production API.
  - [x] Lock the high-pressure complete-condensation lifecycle: liquid
    activation, failed positive-interior coexistence solve, validated
    liquid-only recovery, positive missing-gas TPD, gas deactivation, complete
    inventory conservation, and transactional publication.
  - [x] Add the reverse low-pressure lifecycle from accepted liquid-only
    history: negative gas TPD must permit evaporation and publish a
    gas-containing accepted state. Assert thermodynamic evidence and topology,
    not one incidental sequence of nonlinear restarts.
  - [x] Once the two-sided roots are stable, extend the IAPWS-specific typed
    report with separate liquid-from-gas and gas-from-liquid root deltas only
    if this improves auditability without mixing lifecycle details into the
    generic frozen-reference loader. The ten-point IAPWS diagnostic reports
    zero at printed precision for every I1/I2-to-liquid-TPD,
    I1/I2-to-gas-TPD, and liquid-TPD-to-gas-TPD delta.
- [x] Keep the external IAPWS result in `CharacterizationOnly`: retain visible
  max/RMS external errors but do not freeze the current 1.09--1.47 percent
  discrepancy as an acceptance target. Likely contributors include consistency
  of local NASA gas/condensed standard Gibbs functions, reference-state and
  polynomial approximation effects, and only to a lesser degree at these low
  pressures ideal-gas versus real-fluid behavior. Quantitative attribution is
  a separate validation study.
- [ ] Add a second frozen IAPWS I5 benchmark for `H2O(s, ice Ih) <=> H2O(g)`.
  It complements, rather than dilutes, the liquid-water characterization:
  - [x] Freeze the eight reviewed IAPWS R14-08(2011), section 4 equation (6)
    sublimation rows at 200..270 K with source identity and transcription
    provenance in read-only JSON. Regression tests must not implement the
    IAPWS correlation dynamically.
  - [x] Give ice sublimation a separate typed row contract: finite positive
    values, increasing temperature/pressure, and the official 50..273.16 K
    domain. It must not inherit liquid-above-triple-point semantics merely
    because both tables contain `(T, p)`.
  - [x] Add an IAPWS-specific typed report with `I1/I2`, `TPD(ice|gas)`, and
    `TPD(gas|ice)` roots; strict three-way internal equality; relative, RMS,
    and mean-signed external diagnostics. Keep it distinct from the generic
    frozen loader and from the liquid-only report.
  - [x] Add the ignored real-data diagnostic using log-pressure root search.
    Before solving, intersect the local `H2O(g)`/`H2O(s)` validity interval
    with frozen rows, explicitly report excluded rows, and never extrapolate
    a local polynomial to retain a reference point. At 200 K prove trace
    seeds, log-moles, activities, and `ln(P/P0)` remain finite.
  - [x] Lock direct `gas -> ice deposition` above, and `ice -> gas
    sublimation` below, the *KiThe internal* 250 K boundary. Lifecycle
    pressure offsets must be relative to that internal root so the test
    isolates phase-control from the external model discrepancy.
  - [ ] Record release output in `STORY_TESTS.md`; keep IAPWS comparison
    `CharacterizationOnly` until a source/model-specific discrepancy review
    justifies an external acceptance contract.
- [ ] Add the first Boudouard I4/I5 reaction-thermochemistry layer before any
  Boudouard pressure-boundary, TPD, or lifecycle story:
  - [x] Freeze eleven primary NIST-JANAF rows at 500..1500 K for CO(g) C-093
    and CO2(g) C-095. Pin C(ref), graphite C-002 in provenance and retain the
    primary `Delta_f G` and `log10 Kf` columns rather than a pre-combined Kp.
    The JSON is read-only and tests never access JANAF online.
  - [x] Add a separate typed row schema and a Boudouard-specific comparison
    report. Derive `Delta_r G = Delta_f G(CO2) - 2 Delta_f G(CO)` and derive
    `log10 Kp` independently from frozen Gibbs and frozen log-Kf columns.
    Validate those two JANAF representations under a rounding-aware internal
    contract, not machine precision.
  - [x] Resolve only the existing `RealPurePhaseFamily::BoudouardCarbon`
    fixture, assert exact phase-qualified layout `CO`, `CO2`, `C(gr)`, local
    NASA-gas/NASA-cond provenance, full reaction dimension one, and gas-only
    reaction dimension zero. A graphite identity mismatch must fail rather
    than silently selecting a generic carbon record.
  - [x] Make the JANAF `p0 = 100000 Pa` comparison convention explicit. This
    first layer compares only local standard `G0(T)` closures, therefore it
    applies no hidden ideal-gas pressure correction or 1-bar/1-atm tolerance.
  - [x] Expose the local thermochemistry standard-state pressure as typed
    component provenance. Resolution now reads only explicit JSON fields
    (`standard_state_pressure_pa` or `reference_pressure_pa`) and reports
    `Undeclared` otherwise; it never guesses 1 bar or 1 atm from a NASA/NIST
    polynomial type. The current bundled NASA Boudouard records are explicitly
    observed as `Undeclared`.
  - [ ] Before promoting the separate JANAF pressure-boundary
    characterization into an external acceptance contract, add reviewed
    standard-state-pressure metadata for every selected local record or
    library family. The current `G0(T)` and boundary diagnostics deliberately
    make no conversion or 1-bar/1-atm inference from an undocumented
    convention.
  - [x] Require the local common validity interval to cover every frozen row;
    otherwise report exact excluded temperatures and refuse extrapolation.
  - [x] Add an ignored diagnostic table with JANAF/KiThe reaction Gibbs,
    `log10 Kp`, signed/RMS/max deltas, frozen-file immutability, and external
    `CharacterizationOnly` policy. Debug evidence shows a smooth local bias
    of 41--62 J/mol and no basis to alter records or relax the solver.
  - [ ] Record the release diagnostic in `STORY_TESTS.md`. The separate
    `janaf_boudouard_boundary` module now owns the JANAF-derived `P,T`
    boundary, TPD, and lifecycle characterization; this thermochemistry-only
    layer must remain free of those concerns.
- [x] Add a separate Boudouard `P,T` I5 pressure-boundary and lifecycle layer.
  The existing JANAF module remains thermochemistry-only; this layer may reuse
  its frozen primary CO/CO2/C(gr) rows and derived `Kp(T)`, but must own its
  pressure roots, canonical TPD evidence, and production lifecycle tests.
  - [x] Derive the external 50/50-gas analytical oracle directly from the
    frozen JANAF rows: `P_boundary = 2 * 100000 Pa / Kp(T)`. Do not create a
    second hand-transcribed Kp table or obtain this boundary from KiThe.
  - [x] Use only the existing 800/900/1000 K frozen rows on the first pass;
    validate finite positive composition, Kp, and pressure values, and retain
    `2 CO -> CO2 + C(gr)` orientation with `nu_C > 0`.
  - [x] Build both local roots with `reference_pressure = 100000 Pa`: the
    independent I1/I2 `ln(Q)-ln(K)` root and the canonical `TPD(C | gas)`
    root. Search in log-pressure with typed unbracketed/non-finite/budget
    failure, then require their strict internal agreement.
  - [x] Keep JANAF-to-local boundary deltas `CharacterizationOnly`. The local
    NASA records currently publish `ThermochemistryStandardStatePressure::Undeclared`;
    print that fact and do not disguise it with a 1-bar/1-atm correction.
  - [x] Repeat the structural contract at every boundary fixture: gas-only
    reaction dimension zero and full reaction dimension one. Do not copy the
    water `TPD(gas | condensed)` symmetry story: graphite-only inventory cannot
    represent the same C/O elemental inventory as CO/CO2 gas.
  - [x] Add one 900 K local-lifecycle story after the two roots agree: at
    `3 * P_local`, initially inactive graphite must activate with negative
    pre-activation TPD; at `P_local / 3` it must remain inactive with positive
    final TPD. Both routes must retain conservation, immutable transition
    evidence, and byte-for-byte local/frozen JSON snapshots.
  - [x] Add the direct pressure metamorphic invariant for the fixed 50/50 gas:
    `ln Q(P2) - ln Q(P1) = -ln(P2/P1)`. This protects reaction orientation,
    ideal-gas pressure exponent, and reference-pressure plumbing without
    duplicating a production TPD test.
  - [x] Add a typed table/report and ignored release diagnostic containing
    JANAF Kp/boundary, local I1/I2 root, local TPD root, external errors,
    internal agreement, and source-pressure provenance. Record the release
    result in `STORY_TESTS.md` before treating this evidence as complete.
- [x] Keep current multi-gas semantics explicit and outside the direct-water
  fix. While any ideal-gas phase is active, declared gas phases belong to one
  fixed gas assemblage for stability purposes. Independent immiscible gas
  phases require a future physical/activity-normalization decision.
- [ ] Review the recorded IAPWS diagnostic and set a source/model-justified
  external comparison contract. It must distinguish IAPWS accuracy from the
  present local NASA thermochemistry, ideal-gas activity, and pure-condensed
  liquid model; it must not inherit generic solver tolerances from the loader.
- [ ] After the first dataset format is reviewed, add separate typed schemas as
  evidence requires them: ATcT species thermochemistry and published NASA CEA
  full-equilibrium cases. The JANAF Boudouard reaction-thermochemistry schema
  now demonstrates why heterogeneous tables must not be forced into the
  temperature/pressure row.
  - [x] Define a semantic frozen-row type for the first complete CEA output:
    it carries the HP inputs, final temperature, total `kg-mol/kg`, and named
    species amounts rather than a positional composition vector. Require the
    exact published eleven-component H/O universe, with duplicate, missing,
    unknown, non-finite, and negative entries rejected by the frozen loader.
  - [x] Freeze the NASA CEA Tutorial H2/O2 HP I5 source case from NASA document
    `20240016039` (Leader et al., AIAA SciTech 2025). Keep the external source
    table read-only and explicitly mark its rounded presentation values as
    `CharacterizationOnly`, not a machine-precision acceptance oracle.
  - [ ] Resolve exactly the declared gas `H/H2/H2O/H2O2/HO2/O/O2/O3/OH` and
    the pinned IAPWS-fixture water records `H2O(L)`/`H2O(s)` offline. Record
    component identity, library, record key, physical phase, G/H/Cp capability,
    common temperature interval, and standard-state pressure provenance. A
    missing exact record must return `ValidationNotApplicable`, never select a
    similar species or polymorph.
    - [x] Add the strict offline preflight and prove it refuses the current
      local data rather than silently dropping condensed candidates: the common
      `H2O(L)=[273.15, 600] K` and `H2O(s)=[200, 273.15] K`, so their common
      interval with the otherwise compatible `NASA_gas=[200, 6000] K` records
      collapses to `273.15 K`. It excludes both the 2000 K reactant state and
      the 3181.23 K CEA result. This is a local-data capability gap, not a
      solver failure. It is not grounds for claiming a full eleven-component
      solve. The executable I5 layer may compare the explicitly named nine-gas
      subsystem because both external condensed rows are exactly zero, but it
      must label them `ExternallyAbsentExcluded` and must not claim local TPD
      or phase-stability evidence. Revisit the full universe only with reviewed
      high-temperature condensed-record coverage or an explicit, physically
      justified policy for candidates that are inapplicable before TPD
      construction.
  - [x] Reconstruct the 1 kg H2/O2 reactant mixture from CEA's O/F *mass*
    ratio using the same local molar-mass machinery as phase resolution. Print
    mass, kmol, H/O totals, and the local `H_target` built at 2000 K; keep
    physical `P = 101325 Pa` separate from every record's standard-state `p0`.
  - [ ] Run the ordinary production general `P,H` phase-control workflow over
    exactly that declared universe, initially physical H2/O2 only with trace
    seeds for other gases and inactive condensed water. Report component/
    element counts, matrix rank, reaction-space dimension, actual P,H route,
    accepted phase lifecycle, and strict internal conservation.
    - [x] Add the explicitly scoped gas-only executable layer: resolve exactly
      the nine CEA gas identities offline, reconstruct physical H2/O2 input,
      run canonical `P,H Auto`, and preserve its actual route/fallback evidence.
      This is a fixed one-gas-phase characterization, not the deferred
      eleven-component phase-control proof.
  - [ ] Add a typed identity-based CEA comparison report: temperature and
    total amounts; major/minor/trace species classes; absolute/relative errors
    and `delta_log10` where both amounts are positive; separate condensed
    topology/stability evidence; and distinct external rounded-table versus
    internal strict-conservation diagnostics.
    - [x] The gas-only report compares all nine local amounts by exact CEA
      identity, classifies major/minor/trace rows, reports absolute/relative and
      log-space differences, and retains the two zero condensed rows as typed
      exclusions with reasons. Their topology/stability evidence remains open.
  - [x] Add an ignored offline diagnostic table plus immutable snapshots of
    frozen and local JSON files. Do not establish external numeric tolerances,
    add CEA/Cantera dependencies, permit network access, or auto-expand the
    declared universe on this first complex case.

#### P11.1 - NASA CEA H2/O2 HP gas-only route evidence

The frozen CEA table remains an eleven-row external transcription while the
first executable local problem deliberately contains only its nine gas
components. `H2O(L)` and `H2O(cr)` are `ExternallyAbsentExcluded`, not locally
TPD-validated inactive phases: both external amounts are exactly zero and the
pinned condensed records are outside the required temperature domain. Do not
extrapolate, alter record intervals, introduce pseudo-records, or expand the
published gas universe.

- [x] Keep identity availability separate from thermochemical applicability in
  eleven-row preflight. The full fixture must return
  `ValidationNotApplicable`, never a generic solver error, when any exact
  record does not cover the required temperature.
- [x] Enforce the gas-only eligibility invariant at the executable boundary:
  every locally excluded condensed CEA row must have an exactly zero external
  amount. A future positive external condensed amount must reject gas-only
  construction with typed `ValidationNotApplicable`; it must never be compared
  against a missing local component.
- [x] Extend the typed CEA comparison report with diagnostic-only aggregates:
  major-species max/RMS/mean-signed relative errors; minor/trace max/RMS
  `delta_log10`; absolute and relative temperature/total-amount differences.
  Excluded zero rows are outside all numerical aggregates.
- [x] Add one ignored route-matrix diagnostic for the same exact nine-gas
  fixture. Run `NestedTemperature`, `Monolithic`, and `Auto` separately and
  retain success/failure, accepted path, temperature, residual, elemental
  balance, component comparison, typed failure family, and compact backend
  attempt records. `Auto` must not stand in for an explicit monolithic result.
- [x] If direct monolithic remains rejected, characterize rather than tune:
  use physically meaningful temperature seeds (reactant 2000 K, 3000 K,
  3500 K, and external 3181.23 K only in the ignored diagnostic) plus a small
  trace-floor sweep. Preserve the production trace policy and do not use CEA's
  final composition as a production seed. A broad element-conserving local
  gas seed may be an additional diagnostic probe only.
- [x] When both direct routes accept, compare their nine-component state,
  total amount, temperature, enthalpy acceptance, residual, and balances under
  a strict *internal* tolerance independent of external CEA rounding.
- [ ] Record release route-matrix evidence in `STORY_TESTS.md`. Keep all NASA
  CEA deltas `CharacterizationOnly` until more than one complex reviewed case
  supports a source/model-specific external contract.
- [x] Design a **generic** monolithic P,H temperature-seed recovery policy from
  this evidence. The exact H2/O2 case converges from 3000--3500 K to the same
  accepted state as `NestedTemperature`, while 2000 K exhausts the current
  LM/NR/TR cascade regardless of trace floor or broad local composition seed.
  Do not hard-code a CEA temperature or a case-specific 3000 K retry. First
  define reusable candidate seeds, acceptance/transaction semantics, bounded
  retry budget, route reporting, and cross-fixture regressions.
  - Implemented as `PhMonolithicSeedPolicy`: the caller seed is always first;
    retryable numerical rejection may start fresh transactions at the bounded
    midpoint and quarter points. Candidates are interval-derived,
    duplicate-free, and stop at the first accepted state. Reports retain each
    physical temperature, failure class, backend work, and selected attempt.
    `InitialOnly` preserves strict diagnostics. The CEA H2/O2 case now recovers
    from 2000 K at the generic 3100 K midpoint and remains monolithic.
  - Monolithic recovery now obeys the same declared backend-attempt,
    nonlinear-iteration, wall-time, cancellation, and phase-transition budgets
    as nested P,H. The gate runs on both success and typed error paths before
    `Auto` may select fallback, so an exhausted seed cascade cannot silently
    spend beyond the caller's limit and then start another formulation.
  - `PhMonolithicSeedRecoveryFailed` preserves every seed's original typed
    backend cause. Recursive work counters and presentation rows keep those
    attempts visible to budget enforcement, CLI reports, and GUI diagnostics.
- [x] Document the manual maintainer workflow for adding or revising a frozen
  dataset. Automatic download/update, runtime network access, generated KiThe
  golden snapshots, solver timings, residual dumps, and lifecycle logs remain
  explicitly outside I5.

#### P11.2 - Argonne/STANJAN CHON fixed-`P,T` characterization

This is the first I5 case for the general multicomponent production `P,T`
formulation rather than a pure-phase boundary, scalar reaction, or `P,H`
workflow. The frozen Table 4 source output retains all sixteen published
identities. KiThe deliberately solves an exact 15-component NASA-gas universe:
the source's `C5H12 = 0` row is an external zero reactant not solved locally,
not a missing species and not part of numeric error metrics.

- [x] Freeze Argonne National Laboratory / S. M. Aithal Table 4 provenance,
  `T=2500 K`, `P=35 atm=3546375 Pa`, original pentane-methane-air molecular
  feed, explicit C/H/O/N totals, and all sixteen STANJAN mole fractions as
  read-only typed external evidence.
- [x] Prove in code that the published `C5H12 + CH4 + O2 + N2` feed and the
  local `4 CH4 + 2 CO2 + 8 O2 + 37.6 N2` feed both reconstruct exactly
  `C=6, H=16, O=20, N=75.2`. This proof must not require a local pentane
  thermochemistry record.
- [x] Resolve precisely `CH4/O2/CO2/H2O/N2/N/O/NO/OH/H/N2O/CO/H2/NO2/HO2`
  from offline `NASA_gas`, without NIST fallback or automatic candidate
  expansion; preflight each record's key, physical state, interval, and
  `G(2500 K)`.
- [x] Assert the actual structure is `15 components / 4 elements / rank 4 /
  11 reaction directions`, then run the normal production fixed-`P,T` path
  with zero final phase transitions.
- [x] Add identity-based major/minor/trace diagnostics plus frozen/local JSON
  immutability checks. Preserve STANJAN source rounding (`sum=1.000008107`)
  without renormalizing it and keep all external deltas
  `CharacterizationOnly`.
- [ ] Record the optimized-profile characterization in `STORY_TESTS.md`; do
  not turn the first one-point STANJAN comparison into external acceptance
  tolerances. A future second CHON point or independently reviewed
  thermochemistry source is needed before choosing such bounds.

#### P11.3 - NASA TP-1907 Table 11.3E CHON + graphite multiphase characterization

This is the first I5 benchmark which couples a general CHON+Ar gas reaction
space to a real pure condensed candidate through the production bounded active
set. It is not a scalar Boudouard proxy: the reviewed local universe is 17
NASA-gas components plus `C(gr)`, with five conserved elements and thirteen
reaction directions.

- [x] Freeze the selected NASA TP-1907 Table 11.3E rows at 680, 700, 720, and
  740 K for `H/C=2.000`, `F/A=0.084535`, `ER=1.250`, dry air, and exactly
  101325 Pa. Keep gas and condensed rows as separate typed values; all source
  data, metadata, and normalisation evidence remain read-only and offline.
- [x] State and validate the first source-normalisation contract explicitly:
  the printed gas plus condensed values sum to one to five-decimal rounding
  only when `C(gr)` is included, so the selected Table 11.3E rows are stored as
  `SystemTotalMoleFraction`, never as silently gas-normalised values.
- [x] Resolve the reviewed 17-species NASA gas universe plus a separate
  inactive `C(gr)` NASA-condensed phase with NIST disabled; retain record keys,
  provenance, common temperature coverage, element rank, and reaction-space
  dimension as preflight evidence. `H2O(s)` and `H2O(l)` are typed external
  zero rows and are not extrapolated into the 700--720 K local system.
- [x] Build the source-faithful ordinary-molecule feed from a conceptual `CH2`
  basis, the ordinary `ER=1.25`, and TP-1906 dry air
  (`O2 + 3.727587 N2 + 0.0447068 Ar + 0.0015228 CO2`). `F/A` and chemical ER
  are independently checked rounded diagnostics, never inverse inputs; exact
  C/H/O/N/Ar equality is proven before solving.
- [x] Add an ignored production `P,T` characterization using normal bounded
  phase control with gas initially active and graphite initially inactive. The
  debug route currently gives the expected topology: 700 K activates graphite
  from negative TPD, while 720 K keeps it inactive with positive TPD. It also
  prints source/local system-fraction rows and preserves frozen/local JSON.
- [x] Extend the ignored source-audit diagnostic to all frozen
  680/700/720/740 K rows. It hard-checks only the external topology
  (active, active, inactive, inactive), prints gas and graphite values under
  system-total normalisation, and keeps graphite quantity as characterization.
- [x] Add forward and reverse accepted-state continuation diagnostics across
  all four TP-1907 rows. They distinguish continuation seeds and accepted
  transitions from trial events; both routes reproduce the external active /
  inactive topology order.
- [x] Prove that a one-point typed range is numerically equivalent to an
  independent bounded `P,T` solve at all four frozen temperatures. The initial
  range point has no physical continuation state, and component moles, phase
  status, and graphite system fraction agree before multi-point continuation
  is characterized.
- [x] Locate the internal canonical `TPD(C(gr) | gas)=0` temperature through a
  finite, explicitly bracketed bisection. The first source-faithful result is
  about `705.69 K`, inside the NASA external `[700, 720] K` topology bracket;
  it is explicitly not presented as an interpolated NASA boundary.

#### P11.4 - Frozen-reference software-regression envelopes

Frozen external values remain immutable source evidence. The following guards
are deliberately broad bounds on previously characterized **KiThe software
behavior**, not source uncertainty estimates or new physical acceptance
criteria. Internal I1/I2/TPD, conservation, topology, and file-immutability
contracts remain substantially stricter and independent.

- [x] Keep benchmark-specific envelopes beside their comparator/tests; do not
  serialize KiThe thresholds or expected KiThe numbers into frozen JSON and do
  not parse `STORY_TESTS.md` at runtime.
- [x] Add IAPWS liquid-water and ice-Ih boundary envelopes for maximum and RMS
  external relative error while retaining strict three-route internal roots.
- [x] Add JANAF Boudouard thermochemistry and pressure-boundary envelopes for
  `|delta G|`, `|delta log10 K|`, and external boundary error; do not constrain
  signed bias or weaken the strict local I1/I2-to-TPD agreement.
- [x] Add NASA CEA H2/O2 P,H and Argonne/STANJAN CHON envelopes using existing
  identity-aware aggregate reports. Exclude explicitly out-of-domain condensed
  rows and external zero/reactant rows from numeric aggregates.
- [x] Add TP-1907 guards with hard topology/TPD/root-bracket contracts, broad
  gas aggregate quality bounds, a separate 680 K graphite relative guard, and
  a 700 K graphite absolute system-fraction guard. Never use one ill-conditioned
  near-boundary relative threshold for both temperatures.
- [x] Make every envelope failure report dataset, temperature/species where
  applicable, observed metric, and reviewed guard. Do not freeze accepted
  backend, iteration count, or timing as a quality metric.
- [x] Update `STORY_TESTS.md` to distinguish immutable external source evidence
  from reviewed software-regression envelopes after each guard is implemented.

#### P11.5 - NIST ThermoML benzene/toluene multicomponent candidate phase

This is the first I5 case whose inactive candidate phase has an unknown binary
composition. It must exercise the normal `IdealGas` plus `IdealSolution` path,
not a benchmark-specific VLE model. Experimental ThermoML P-x data characterize
the ideal model; strict algorithmic evidence comes from independent Raoult
arithmetic versus canonical TPD composition minimization.

- [x] Inventory exact offline records at 353.15 K: `NASA_gas:C6H6`,
  `NASA_gas:C7H8`, `NASA_cond:C6H6(L)`, and `NASA_cond:C7H8(L)` all cover the
  target temperature. Keep `nuig_thermo` as an explicit future fallback only if
  a required native NASA-condensed record becomes unavailable.
- [x] Freeze the reviewed NIST ThermoML P-x subset and provenance locally. It
  contains only published liquid benzene composition and total pressure; do not
  invent an experimental vapour composition or add network access to tests.
- [x] Add the independent Raoult bubble/dew oracle from the frozen pure
  endpoints, including exact algebraic round-trip tests.
- [x] Resolve `C6H6`/`C7H8` gas and `C6H6(L)`/`C7H8(L)` liquid states from
  local `SubsData` with explicit offline state/model declarations and no NIST
  fallback. `nuig_thermo` remains a future fallback only if the reviewed NASA
  condensed records become unavailable.
- [x] Add canonical gas and liquid candidate TPD-minimization tests at three
  interior ThermoML compositions. Require both minimum TPD and recovered
  argmin composition; use absolute simplex errors, not relative errors.
- [x] Prove the pressure sign and bounded phase activation around the local
  bubble pressure: gas remains stable below it and the binary liquid appears
  above it from the TPD minimizer composition.
- [x] Add an ignored release characterization table for local NASA Raoult P-x
  versus frozen NIST P-x. It is source-comparison evidence, not a strict
  regression tolerance or a replacement for I1/I3 proof.
- [x] Add two-phase chemical-potential equality, phase-qualified identity, and
  per-species/element conservation stories after liquid activation. The
  molecular-species check is phase-aware (gas plus its declared liquid peer),
  not an ambiguous aggregate over arbitrary same-named records.

#### P11.6 - Ternary VLE 2D-simplex preflight

The first ternary candidate is `toluene + ethylbenzene + chlorobenzene`, with
one ideal gas phase and one three-component ideal liquid solution. This pass is
strictly a capability and source-data audit. Do not create a ternary
phase-control benchmark, add activity coefficients, interpolate source rows,
or substitute chemically similar compounds until the preflight reaches a clear
Outcome A.

- [x] Inventory the six exact local gas/liquid records with NIST fallback
  disabled: canonical identity, selected library/key, physical state, `G(T)`/
  `H(T)` support, standard-pressure provenance, and each temperature interval.
  The executable probe finds only `NASA_gas:C7H8`,
  `NASA_cond:C7H8(L)`, and `NASA_gas:C8H10,ethylbenz`. Exact
  `NASA_cond:C8H10(L),ethylbenz`, `NASA_gas:C6H5Cl`, and
  `NASA_cond:C6H5Cl(L)` are absent.
- [x] Reject the candidate with `ValidationNotApplicable` if any exact state is
  missing. The typed preflight refuses substitutions and proves that no local
  fallback path enables online NIST lookup.
- [x] Compute and report the common intersection of all six local intervals;
  no experimental condition may be selected before this intersection is known.
  It is intentionally unavailable because the six-state inventory is incomplete.
- [x] Audit the official machine-readable NIST ThermoML payload for the exact
  ternary: compound order, property type, published `T/P/x/y` variables,
  uncertainties, and number of genuine ternary interior rows. DOI
  `10.1021/je020186c` has 48 ternary isobaric `T-x` rows at 26.66, 53.33,
  79.99, and 101.32 kPa; it publishes liquid `x` and temperature uncertainty,
  but no experimental vapor composition.
- [x] Identify official pure-component saturation support for the later
  independent ternary Raoult oracle. The reviewed payload itself provides no
  pure-saturation table, so a future executable fixture must freeze a separate
  authoritative source.
- [x] If any requirement fails, record the exact missing capability and stop
  before TPD/lifecycle implementation. The primary candidate is not executable
  with the current repository; no frozen ternary rows or artificial simplex
  benchmark are added.
- [ ] After adding exact offline liquid ethylbenzene plus gas/liquid
  chlorobenzene records, rerun the preflight, require a nonempty common interval,
  select four to six genuine ThermoML interior rows, freeze only published
  `T/P/x` values, and then implement the 3-component / 2D-simplex TPD story.
- [ ] If the missing local records are not added, inventory actual offline
  gas/liquid pairs first and choose a different ternary NIST VLE source from
  that repository-driven shortlist. Do not select a literature case before the
  local six-state inventory succeeds.

#### P11.7 - Explicit test-only frozen thermochemistry for the ternary VLE candidate

The production-only P11.6 preflight is a negative capability contract and must
remain so.  A later test-only thermochemistry universe may complement it, but
must never enter `SubsData` search, mutate JSON libraries, or turn the
production preflight from `ValidationNotApplicable` into `Ready`.

- [ ] First make a reviewed source inventory for the three missing exact
  states: liquid ethylbenzene, gas chlorobenzene, and liquid chlorobenzene.
  For every candidate source record, retain CAS/formula/state, primary source
  and table, representation, temperature range, `Cp/H/S/G` capabilities,
  formation/reference convention, and standard pressure.  Prefer reviewed
  NASA/CEA, NIST/TRC, JANAF/equivalent, or the original peer-reviewed source;
  reject undocumented aggregators.  Do not fit a closure to the ternary
  ThermoML VLE rows.
- [x] Freeze the approved NIST WebBook ethylbenzene numerical seeds separately:
  liquid/gas formation-enthalpy anchors, liquid entropy and Cp anchors,
  vaporization-enthalpy points, the bounded Majer-Svoboda correlation, and the
  bounded Antoine pressure oracle.  This is I5 source evidence only; it does
  not yet claim a gas/liquid `G0(T)` closure or alter the production preflight.
- [ ] Decide whether production NASA records may be combined with external
  frozen records only after an explicit reference-state alignment audit:
  pressure, enthalpy/entropy zero, formation convention, and elemental
  reference state.  Present NASA records report an *undeclared* machine
  pressure, so this decision must not be inferred silently.  Any `1 atm` to
  `1 bar` correction must be explicit, documented, and unit-tested; otherwise
  use a compatible external gas/liquid pair or retain `ValidationNotApplicable`.
- [ ] Add a separate, read-only frozen dataset only after the inventory passes.
  Give it semantic external provenance (never a forged `NASA_gas` or
  `NASA_cond` identity), metadata beside the numerical representation, exact
  state labels, temperature bounds, and immutable-byte coverage.  It is
  test-only and must be selected by an explicit adapter, never a global
  fallback.
- [ ] Build the adapter around the existing generic
  `ResolvedThermochemistry::from_functions` capability contract:
  `G0(T)`, optional `H(T)`/`Cp(T)`, bounds, state, reference pressure, and
  provenance.  It must contain no compound-specific branches and must return a
  typed out-of-range error instead of extrapolating.
- [ ] Before any ternary lifecycle claim, validate each gas/liquid pair with
  `Delta G_vap(T) = G0_gas(T) - G0_liquid(T)` and independently compare the
  implied pure saturation pressure with authoritative pure-component vapour
  pressure data.  This evidence route must be independent from ThermoML DOI
  `10.1021/je020186c`.
- [ ] Add the frozen-capability preflight beside the preserved production-only
  negative preflight.  It must require all six exact state-qualified
  components, three liquid components / a two-dimensional simplex, finite
  `G0` (and `H` where supplied), shared temperature coverage for selected
  source rows, preserved provenance/reference pressure, and no library/data
  mutation.
- [ ] Only after those gates pass, freeze four to six published interior
  ThermoML `T/P/x` rows as characterization evidence and add the normal
  gas-only / boundary / two-phase TPD lifecycle stories.  Do not claim an
  external VLE error envelope before the independent pure-state checks are
  established.

#### P11.8 - Ternary `P,T` Antoine-gauge fixture

The first strict ternary algorithmic test need not wait for an absolute
`Cp/H/S` library.  With only toluene, ethylbenzene, chlorobenzene, and their
gas/liquid copies, the molecular element matrix has rank three and the
six-component system has three transfer directions.  An explicit test-only
gauge therefore supplies the required standard-Gibbs differences from
independent pure `Psat(T)` while remaining prohibited from `P,H` and
production lookup.

- [x] Freeze generic NIST Antoine records for all three identities with exact
  CAS/formula/elemental composition, `log10(P/bar)` coefficients, source route,
  and validity ranges.  Derive, rather than hard-code, their common interval
  `335.19..384.66 K`; use explicit `p0 = 100000 Pa` throughout.
- [x] Implement bounded generic `G0_gas=0` / `G0_liquid=RT ln(Psat/p0)` gauge
  functions plus the pure pressure round-trip.  Test gas/liquid sign semantics,
  extrapolation rejection, source immutability, molecule-matrix rank three,
  and exactly three duplicated phase-transfer directions.
- [x] Audit genuine ThermoML interior rows in the derived common interval and
  freeze only published `T/P/x` values. Four source rows at the central and
  ethylbenzene-rich compositions, each at `26.66` and `53.33 kPa`, lie inside
  the common window. The Raoult-derived `y` remains a solver-side independent
  result, never a synthetic source column.
- [x] Add the scalar ideal-Raoult bubble-temperature root and its derived
  vapour composition as I1/I2 evidence before invoking canonical TPD.
- [x] Run canonical liquid and vapor candidate TPD at the same boundary:
  require a three-component candidate, simplex dimension two, `TPD_min ~= 0`,
  and recovered argmin rather than an expected-composition seed. This is a
  direct test-only `IdealTpdProblem` route over the frozen gauge, not a fake
  `ResolvedPhaseSystem` or a legacy mutable workflow.
- [x] Add gas-to-liquid and liquid-to-gas active-set decision stories around
  the pressure boundary, preserving the complete incipient ternary
  composition. They exercise canonical stability reporting and `PhaseManager`
  transition classification without falsely packaging test-only gauge data as
  a production-resolved system.
- [x] Characterize, but do not use as an acceptance envelope, `T_NIST -
  T_Raoult` for selected experimental rows.  The ignored story prints every
  source point plus RMS/max temperature deltas. Keep this external comparison
  distinct from strict Antoine-gauge I1/I2 and canonical I3 evidence.

#### P11.9 - Full canonical lifecycle over the ternary Antoine gauge

The gauge is deliberately a `P,T`-only test universe, but it must now prove
that the normal immutable phase-control runner carries a multicomponent
candidate through activation, accepted re-solve, continuation, and reduction
of topology. `PreparedPhaseControlRunner` over an explicit raw
`EquilibriumProblem` is the canonical test boundary here; never substitute a
fake `ResolvedPhaseSystem` or the legacy mutable workflow.

- [x] Add an independent generic ternary Rachford-Rice oracle with typed
  `TwoPhase`, `AllLiquid`, and `AllVapor` outcomes. Validate exact simplex,
  reconstructed bulk inventory, endpoint classifications, and an interior
  `beta` reference built from `beta*y + (1-beta)*x`.
- [x] Build one central `P = 53.33 kPa` benchmark inventory from the
  independently derived Raoult `x/y` and a strictly interior vapor fraction.
  Keep `x`, `y`, and bulk `z` as distinct quantities.
- [x] Run gas-only -> liquid activation through `PreparedPhaseControlRunner`:
  require a negative liquid TPD, a three-component incipient composition,
  an accepted two-phase output, molecular conservation, and cross-phase
  chemical-potential equality. Compare `beta/x/y` with Rachford-Rice.
- [x] Mirror liquid-only -> gas activation with the same accepted-state and
  independent-flash contracts.
- [x] Demonstrate at least one accepted multicomponent disappearance with
  finite final absence TPD evidence. Permit either supported production
  disappearance reason; do not force one numerical branch.
- [ ] Add a short forward/reverse fixed-pressure temperature story using only
  accepted continuation state. The first accepted two-phase -> gas-only
  continuation handoff, topology reduction, and conservation are now covered;
  reverse sweep and no-chatter evidence remain.
- [ ] Locate and exercise one genuine in-band hysteresis point from active and
  inactive histories. Outside the band, require history-independent topology.
- [ ] Add an ignored compact release diagnostic table with phase statuses,
  continuation use, flash classification, `beta`, TPDs, transition/trial
  counts, maximum composition/conservation deltas, and separate forward /
  reverse transition locations.
- [ ] Keep the established external ThermoML temperature characterization as
  a conservative software-regression guard only (`RMS < 1.5 K`,
  `max |delta| < 2.0 K`), never as the lifecycle truth.

#### P11.10 - Fixed-inventory full ternary production lifecycle

The earlier P11.9 story proves the runner can reach an independently checked
two-phase split from either one-phase start. This stricter follow-up fixes the
physical inventory to `z = [0.334, 0.333, 0.333] mol` at `P = 53_330 Pa` and
uses independently recomputed Antoine/Raoult bubble and dew boundaries. It is
still test-only `P,T` ideal physics: no Antoine data, phase-control policy, or
production code may be adjusted merely to obtain an expected topology.

- [x] Add an independent bounded dew-temperature oracle alongside the existing
  bubble and Rachford-Rice routes. Recompute and characterize `Tb`, `Td`, the
  endpoint incipient compositions, and the physical two-phase interval from
  frozen source data rather than treating approved diagnostic numbers as truth.
- [x] Add the exact seven-point fixed-inventory grid (`Tb-3`, `Tb-0.5`,
  `Tb+0.5`, midpoint, `Td-0.5`, `Td+0.5`, `Td+3` K) and prove all independent
  Rachford-Rice classifications before production lifecycle is invoked.
- [x] Run the full forward production sequence liquid-only -> two-phase ->
  gas-only. The first point must be independent; later points may receive only
  the previous accepted log-mole/phase-set continuation. Require appearance
  between points 2/3, disappearance between points 5/6, no accepted chatter,
  molecular conservation, and final complementarity.
- [x] Compare every accepted two-phase forward point with independent flash:
  `beta`, phase-qualified `x/y`, per-component conservation, and cross-phase
  chemical-potential equality. The midpoint is the primary strict comparator.
- [x] Run a new independent reverse sequence gas-only -> two-phase ->
  liquid-only. Require liquid appearance, gas disappearance, and equality of
  forward/reverse interior physical states; transition locations may differ
  only within the declared hysteresis band.
- [x] Find a genuine continuous-TDP point strictly inside the actual production
  hysteresis band. At the same `P,T,z`, prove inactive and active accepted
  histories retain different allowed topologies while their TPD values agree;
  add outside-band history-independent controls. A dew-side analogue is
  optional if it adds no distinct contract.
- [x] Add an ignored release table containing forward/reverse point topology,
  oracle regime, continuation origin, phase fraction, TPDs, transition counts,
  composition and conservation deltas, chemical-potential mismatch, and
  aggregate no-chatter metrics. The release execution is recorded in
  `STORY_TESTS`.
- [x] Keep frozen Antoine and ThermoML immutability plus the existing external
  ThermoML `RMS < 1.5 K`, `max |delta| < 2.0 K` software-regression guard.

#### P11.11 - NASA TP-1906/1907 CHON + graphite `P,H` lifecycle

This is the first heterogeneous external `P,H` story. NASA TP-1906 supplies
the target **specific equilibrium-mixture enthalpy**; TP-1907 separately owns
composition and graphite-topology evidence. The production request receives
only pressure, extensive `H_target`, and the closed element inventory. The
published temperature is comparison evidence, never the answer or a
benchmark-specific seed.

- [x] Freeze the four TP-1906 Table 11.3E heterogeneous `H [J/g]` rows at
  680/700/720/740 K with separate provenance and explicit source units. Do
  not merge them into TP-1907 composition JSON or reinterpret them as molar,
  gas-only, or frozen-composition enthalpy.
- [x] Reuse the reviewed TP-1906 dry-air / TP-1907 executable feed and prove
  semantic joins by temperature, pressure, H/C, ER, chemical ER, and dry-air
  convention. Derive each total `H_target [J]` from the exact executable
  inventory mass, retaining source `J/g`, mass `g`, and total `J` in reports.
- [x] Add a fixed-`P,T` enthalpy-reference preflight at every source row.
  Characterize local `h [J/g] - h_TP1906 [J/g]` before interpreting any
  recovered `P,H` temperature. Never apply an empirical enthalpy offset.
- [x] Add isolated canonical `P,H` cases and compare solved temperature,
  topology, graphite amount, system-normalized gas composition, enthalpy
  residual, conservation, and a fixed-`P,T` witness at recovered temperature.
- [x] Add forward `H680 -> H740` and fresh reverse `H740 -> H680` ranges.
  Only accepted states may continue; retain accepted transitions separately
  from nested trial events and require the expected single graphite
  disappearance/appearance with no accepted chatter.
- [x] Add route/seed evidence for nested, monolithic when eligible, and Auto.
  The diagnostic uses a common 900 K numerical seed, records the accepted
  route, and never provides a source temperature as an answer or seed. Its
  monolithic finding is retained as the separate blocker below.
- [x] Repair bounded monolithic `P,H` phase control before promoting `Auto` for
  near-boundary condensed phases. The defect was in the shared row-scaling
  contract, not TPD: the dimensionless reaction residual was divided by a
  dimensional `max(|Delta G0|, R*T*||nu||)` factor, so a raw affinity error of
  order one could pass a `1e-6` acceptance gate. Reaction rows now remove only
  the arbitrary reaction-basis norm. The ordinary H700 regression requires
  monolithic `P,H` to retain finite graphite and agree with an independently
  solved canonical `P,T` witness to `1e-6` relative moles.
- [x] Add an ignored release characterization table and only then review
  conservative software-regression envelopes. Snapshot frozen TP-1906,
  TP-1907, and local library files around all offline stories.

#### P11.12 - Frozen-reference numerical formulation hardening

The reviewed datasets now exercise production pathways. Use them to lock down
coordinate and lifecycle invariants that synthetic fixtures cannot prove.

- [x] Add a real TP-1907 fixed-active reaction-basis metamorphic regression:
  permute and rescale non-zero reaction columns while preserving `A^T*N = 0`,
  then require the same accepted mole vector, affinity, and element balance.
- [x] Add a real bounded TP-1907 trace-floor regression on both sides of the
  graphite transition. Floors used only for log coordinates must not change
  accepted topology or materially change component amounts.
- [x] Promote the four TP-1906/1907 P,H source rows into a normal,
  non-printing regression: every isolated target and both accepted-state
  continuation directions must retain the P,T witness manifold, graphite
  topology, enthalpy contract, and frozen/local file immutability.
- [ ] Promote selected ignored external-characterization envelopes into small
  ordinary assertions where their contracts are now stable: NASA CEA H2/O2
  P,H, Argonne/STANJAN CHON, IAPWS liquid/ice, and the ternary external
  temperature envelope. Keep detailed release tables ignored.
- [x] Add a `FrozenReferenceCatalog` that verifies every dataset identifier,
  evidence kind, schema version, expected row count, and source file is unique
  and referenced by its intended adapter. Reject unsupported positive schema
  versions instead of accepting every non-zero value.
- [x] Add real-data extensive-inventory-scaling metamorphic coverage for both
  TP-1907 P,T and TP-1906 P,H. Scaling the full closed inventory and the
  extensive P,H target must retain temperature/topology and scale every amount.
  The fixed-`P,T` reference covers `10^-3..10^3`; bounded nested `P,H` covers
  the controlled two-order `10^-1..10^1` probe with a tightened scalar/inner
  contract. Extreme `P,H` brackets remain a separate cascade-conditioning
  concern, not a reason to weaken the accepted-state invariant.
- [x] Expand the TP-1907 CHON + graphite extensive-scaling evidence without
  tuning production thresholds to the fixture.
  - [x] Add a lower-level reaction-row normalization matrix (`10^-8..10^8`
    and a sign reversal) that proves raw affinities and their row scales change
    together while the dimensionless residual remains invariant.
  - [x] Characterize fresh bounded P,T solves at 680/700/720 K over
    `10^-4..10^4` inventory factors. Compare topology, phase totals, gas
    composition, TPD sign/value, conservation, and component moles by typed
    identity. Print the active-phase amount next to the numerical trace floor
    and `phase_eps` so a threshold interaction is visible rather than hidden.
    The ignored matrix is now a strict production regression: invariant
    accepted states cover every factor through `10^4`, with active phase
    totals still far above `phase_eps`. The original physical-coordinate
    failure at `10^4` remains an explicit negative witness when
    `ExtensiveNormalizationPolicy::Disabled` is selected. A
    scale-aware trace seed was explicitly tried and does not remove the
    `10^4` failure. Fixed-active diagnostics show that a fresh `10^4`
    inventory seed fails while `accepted_log_moles + ln(10^4)` accepts the
    identical scaled formulation. This rules out the trace floor and
    phase-lifecycle thresholds as the direct cause, but is only evidence of a
    fresh-seed conditioning/basin gap: it does not yet localize one defect
    across legacy and RST backends. A first solver-coordinate-origin
    experiment did not recover the fresh case robustly and was reverted; no
    production coordinate contract changed as a result.
    - [x] Isolate the `10^4` fresh-seed failure with minimal per-backend
      reproductions, including initial residual, finite bounds, step-control,
      and candidate-rejection evidence. Only after that may we reconsider a
      solver-internal normalized log-mole coordinate origin. The resulting
      production route proved recovery of both fixed-active and bounded
      TP-1907 cases while preserving public physical `ln(n)`/mole results.
    - [x] Diagnose the real TP-1907 `10^4` fresh-inventory gap before changing
      any production numerical policy.
      - [x] Add a diagnostic-only classification with mutually exclusive
        results: `StrictScaleInvariant`, `ThresholdLimited`,
        `FreshSeedBasinLimited`, and `PhysicalScaleRegression`. Do not assign
        `FreshSeedBasinLimited` until a transformed accepted seed proves the
        scaled physical state exists. The 680/700/720 K matrix now assigns
        `FreshSeedBasinLimited` only after its temperature-appropriate
        transformed physical oracle has been accepted.
      - [x] Extend the fixed-`P,T` diagnostic to 680, 700, and 720 K: the
        first two retain graphite, while 720 K excludes it. For each point
        compare (A) ordinary fresh production seed, (B) exact `ln(10^4)`
        translation of an accepted unit-scale state, and (C) an
        answer-independent scale-aware candidate seed. Route B is an oracle
        only and must never become a production solving path. At 720 K the
        oracle is correctly reduced to gas-only: the graphite coordinate stays
        numerical trace, the accepted topology is gas-only, normalized gas
        moles/fractions agree with unit scale, and `TPD(C(gr)) > 0` remains
        intensive (within the diagnostic nonlinear-route envelope).
        - [x] The fixed-active 680/700 K pass now proves
          `FreshSeedBasinLimited`: ordinary fresh and both input-only trace
          variants fail, while the transformed accepted seed preserves moles,
          reaction affinity, and balances. Lowering the absolute trace floor
          to `1e-36` and raising it relatively both fail, so a trace-floor
          policy change is not a justified production fix.
        - [x] Complete 720 K through the physically correct reduced gas-only
          projection. The all-active fixed formulation is intentionally not a
          valid oracle there because graphite is inactive; its iteration-limit
          failure must not be counted as the scale-conditioning result. The
          transformed `10^4` gas-only seed now accepts with gas active and
          graphite inactive, preserves normalized gas moles and fractions,
          and retains positive intensive graphite TPD. Ordinary fresh bounded
          production with recovery disabled still fails, while the default
          typed recovery accepts through a physical-coordinate retry. Thus
          720 K independently confirms `FreshSeedBasinLimited` without a
          graphite-activation artifact.
      - [x] Measure the actual seed at every relevant boundary before inventing
        a new policy. `LogMolesInitialGuess::from_moles_with_policy` already
        translates positive physical initial amounts under uniform scaling;
        trace coordinates deliberately remain governed by their explicit
        floor. Record where that covariance is lost, if at all: global/reduced
        active seed, phase-control restart, or backend iteration.
        - [x] Component-level 680/700 K instrumentation confirms that every
          positive physical input coordinate shifts by `ln(10^4)` to machine
          precision. The proposed input-total-normalized Route E has scale
          exactly one and is bitwise-equivalent to ordinary Route A, so it
          cannot recover the case and must not be promoted to production.
          The accepted oracle instead has a strongly different reactive
          composition (many input-zero species become major and some input
          reactants become trace). The current explanation is therefore
          relative-composition/basin conditioning, not lost global extensive
          scale or trace-floor magnitude.
      - [x] Print structured provenance for every A/B/C route: initial
        log-mole min/max/mean/range, trace-coordinate count, raw and scaled
        reaction/element residual blocks, finite bounds, backend attempts,
        iteration counters, terminal typed error, and component-level
        recovered-mole mismatch. The release-only diagnostic
        now prints seed and component tables plus A/B residual blocks and
        terminal solver evidence. The 720 K row explicitly marks C/D as not
        applicable because 680/700 already falsify the trace-floor hypothesis;
        it records the reduced gas-only B route and graphite TPD instead. Do
        not tighten a global error tolerance until a new witness identifies a
        different component-level mismatch.
      - [x] Establish the production representation boundary for exact
        extensive normalization and select an explicit on-failure policy.
        - [x] Introduce a typed answer-independent `ExtensiveNormalization`
          built only from physical input inventory. It must translate component
          moles, element/phase totals, total enthalpy, extensive residuals,
          and physical absolute thresholds explicitly rather than scattering
          raw `/ scale` operations. It is public through the equilibrium
          prelude but contains no solve or routing policy.
        - [x] Publish a concise audited classification: extensive physical
          values transform; T/P/p0, mole fractions, activities, G/H/S standard
          states, reaction affinity, chemical potentials, TPD, hysteresis, and
          topology remain intensive/invariant. Record separate semantics for
          trace floor, `phase_eps`, absolute acceptance tolerances, and finite
          log-mole floors. The audited table lives with `ExtensiveNormalization`.
        - [x] Preserve physical absolute-mole threshold meaning in normalized
          solver space: map `trace_floor` and `phase_eps` to `/ scale`, keep
          relative trace fractions unchanged while scaling their absolute cap,
          and provide the same explicit conversion for absolute energy/mole
          tolerances. Do not scale TPD create/keep thresholds. The typed phase
          policy adapter changes only `phase_eps`; P,H test requests convert
          absolute enthalpy tolerance explicitly.
        - [x] Keep normalized diagnostics internal and reconstruct physical
          extensive balance/enthalpy evidence before publication. Discovery
          diagnostics are disabled so live sinks cannot receive normalized
          mole units. `ExtensiveNormalizationRecoveryEvidence` retains the
          scale, original typed failure, discovery backend, optional physical
          retry backend/failure, and whether physical publication required
          reconstruction. Frozen P,H evidence prints normalized and
          reconstructed physical enthalpy errors as distinct values.
        - [x] Define the continuation boundary before routing: accepted phase
          topology and physical composition remain physical state; a normalized
          numerical iterate alone must never become a continuation seed. The
          mapping documents this rule and does not alter continuation behavior.
        - [x] Add focused algebraic tests for round trips, scale invariance,
          enthalpy conversion, physical threshold conversion, and non-scaling
          of intensive values, then migrate frozen diagnostics to this shared
          abstraction without changing default solve paths. Unit coverage also
          rejects invalid scales/inventories and proves trace/phase-epsilon
          semantics; P,T/P,H ignored recovery matrices use the shared mapping.
      - [x] Test exact extensive normalization as a separate test-only
        formulation recovery. Derive `s = sum(n0 > 0)` from the current
        request, rebuild the same fixed-`P,T` problem with `n0 / s`, use its
        ordinary fresh seed, and reconstruct physical moles only after
        acceptance. The `10^4` TP-1907 matrix now accepts at 680/700 K with
        gas+graphite and at 720 K with the correct gas-only reduced set:
        reconstructed moles agree with B within `4.6e-9`, gas fractions within
        `2.6e-15`, and the 720 K graphite TPD remains positive/intensive
        (relative difference `2.5e-6`). This is evidence for absolute
          extensive-scale conditioning. The later production policy uses this
          exact request-derived transform only after a classified numerical
          failure; it does not depend on this oracle.
        - [x] Add the full `10^-4,10^-2,1,10^2,10^4` P,T normalization
          metamorphic matrix for 680/700/720 K. Each internal request must
          have an order-one inventory, reconstruct the appropriate B oracle,
          preserve the 720 K gas-only topology and positive graphite TPD, and
          retain the existing F/I/B negative evidence at `10^4`. The ignored
          matrix now passes for all fifteen points: the largest reconstructed
          mole mismatch is `1.7e-8` (720 K at `10^-4`), gas fractions agree
          within `4.1e-10`, and the intensive graphite-TPD difference remains
          `2.5e-6`.
        - [x] Only after the P,T scale matrix is complete, characterize nested
          P,H normalization on H700/H720 by applying the same `n0 / s` and
          `H_target / s` transformation. Reports must distinguish normalized
          internal enthalpy residual from reconstructed physical Joules; no
          P,H production retry may be introduced in this diagnostic pass. At
          physical factor `10^4`, ordinary fresh nested P,H still fails while
          normalized ordinary nested P,H accepts both rows. H700 reconstructs
          the graphite-present base state with max mole error `1.5e-9`; H720
          reconstructs the graphite-absent state with `1.3e-11`; both recover
          exactly the base temperature within the printed precision. The
          diagnostic publishes internal (`~4e-5 J`) and reconstructed physical
          (`~2..3 J`) enthalpy errors separately, so normalized Joules cannot
          leak into physical evidence.
        - [ ] Inspect why failed backend attempts currently publish
          `iterations=0` with `metrics=None`: distinguish initialization,
          Jacobian/linear-solve, bounds/step, and iteration-limit failure.
          This is solver observability work independent of normalization.
      - [x] Classify the `P,H` `10^4` failure as an inner bounded-`P,T`
        conditioning failure rather than a missing outer enthalpy root. The
        canonical nested route now recovers its numerical inner trials through
        the shared P,T normalization boundary and accepts both H700/H720 while
        preserving temperature, topology, physical moles, and enthalpy.
      - [x] Reject the conditional Route-C seed branch: answer-independent
        interior seeds did not succeed while A failed. The selected recovery
        is exact extensive representation normalization, not a new guessed
        equilibrium composition.
      - [x] Only if answer-independent Route C succeeds while A fails, design
        a deterministic input-derived active-coordinate seed/recovery policy.
        It must preserve `10^-4..10^2` behavior, leave numerical trace floors
        explicit, avoid magic shifts, and never depend on a previously solved
        scale-one reference.
        - [x] Since direct seed replacement was rejected, investigate a
          chemically feasible answer-independent *relative-composition* seed:
          it must respect the current closed element inventory and active phase
          set without borrowing a unit-scale equilibrium, external reference
          fractions, or a temperature-range continuation result. Compare its
          initial residual and element-manifold distance against A/B before
          adding any production recovery path.
        - [x] Prototype the recovery as test-only direct `A^T n = b`
          interiorization, not as a second equilibrium solver. Start with a
          deterministic bounded affine projection of a uniform active-species
          target: keep every inactive phase at numerical trace, enforce a
          dimensionless positive floor only for the current active species,
          and iteratively bind coordinates that would violate that floor.
          Normalize by active input inventory before SVD solves. No Gibbs,
          `K_eq`, accepted state, external composition, or reaction-basis
          coordinates may enter this prototype. The diagnostic helper now also
          supports a floor-injected input-anchored target; the injection occurs
          in normalized coordinates, so raw zero feed components cannot make
          bounded-affine active-set selection scale-sensitive.
        - [x] Before any nonlinear retry, add direct seed evidence for
          680/700/720 K and scale factors `10^-4`, `1`, `10^4`: finite active
          coordinates, exact element conservation, inactive graphite retained
          at 720 K, deterministic output, and seed-scale covariance. Record
          requested versus achieved dimensionless interior fraction and a
          typed diagnostic outcome when a full positive interior is not
          structurally available. Both uniform and input-anchored direct seeds
          now have a full `1e-6` interior at all three points, preserve the
          active/inactive graphite scope, and scale through `10^-4..10^4`.
        - [x] Compare ordinary fresh (F), element-feasible interior (I), and
          transformed accepted oracle (B) at `10^4`. Print raw/scaled residual
          blocks and log-space geometry. I is evidence only until it reaches
          the same physical topology/composition/TPD contract as B with the
          unchanged strict LM configuration. At 680/700 K the input-anchored
          I lowers the raw reaction block from about `65` to `56..59`, while
          uniform I remains about `63..65`; both preserve element balance near
          `1e-11` and nevertheless make LM/NR/TR hit their iteration limit.
          At 720 K the gas-only I likewise fails while gas-only B succeeds and
          retains positive graphite TPD. Therefore direct manifold feasibility
          and the numerical trace floor are not the missing recovery mechanism.
        - [x] Reject further answer-independent relative-composition seed
          invention for this failure after preserving the negative evidence.
          Exact request-derived extensive normalization solves the conditioning
          problem without tuning `phase_eps`, trace floors, acceptance
          tolerances, or borrowing B's accepted composition.
      - [x] Promote the
        full `10^-4..10^4` 680/700/720 K P,T matrix to strict metamorphic
        evidence, then rerun scaled H700/H720 `P,H` cases. Do not change
        `phase_eps`, trace floors, hysteresis thresholds, thermochemistry,
        element totals, or acceptance tolerances merely to admit `10^4`.
  - [x] Add the corresponding P,H
    matrix for the graphite-present and graphite-absent target rows. Scale the
    full inventory and `H_target` together; carry only accepted seeds and
    report route selection separately from the physical invariant. The ignored
    strict 700/720 K nested-P,H matrix now accepts every factor through `10^4`;
    it asserts recovery provenance at `10^4` and distinguishes physical retry
    from reconstructed physical publication.
  - [x] Add typed `ExtensiveNormalizationPolicy::{Disabled,
    OnNumericalFailure}` to the canonical P,T solve options and reproducibility
    snapshot. Recovery is attempted only for classified numerical failures and
    only when the request-derived inventory scale differs materially from one.
  - [x] Prefer a final physical-coordinate retry after normalized basin
    discovery. If the unchanged physical formulation remains ill-conditioned,
    reconstruct every public extensive field through one audited boundary and
    retain the retry failure in typed recovery evidence.
  - [x] Turn the P,T and P,H scale characterizers into strict regressions:
    a failed point now fails the test instead of merely printing `FAILED`.
    The default route accepts all current `10^-4..10^4` rows; an explicit
    disabled-policy assertion preserves the historical `10^4` failure witness.
  - [x] Extend the same representation policy to the prepared typed P,T
    temperature-range runner. The cached formulation remains the fast path;
    only a failed point enters the canonical recovery transaction with its
    current accepted log-mole seed and, for bounded solves, accepted
    `PhaseSet`. `RecoveryFormulation` makes the extra build auditable, range
    progress is emitted exactly once per point, and a failed recovery publishes
    neither a point nor a partial `TemperatureRangeSolution`. Recovered-point
    wall time includes both the rejected prepared attempt and recovery. The strict real
    TP-1907 `680/700/720 K` range now recovers factor `10^4` at its first point
    and then continues through the prepared physical formulation.
  - [x] Add a retained physical-unit recovery event/report for opt-in
    diagnostics. Internal discovery events are suppressed, while
    `ExtensiveNormalizationRecoveryAccepted` records scale, trigger, discovery
    backend, optional physical retry, and reconstruction. The event respects
    diagnostics mode/event limits, reaches the live sink, renders in CLI and
    GUI diagnostics, and is backed by the full typed recovery evidence.
  - [ ] Record the strict production P,T and nested P,H matrices in release
    mode in `STORY_TESTS.md`; the development regressions are green.
  - [ ] Keep `10^-6..10^6` as an ignored release characterization and
    `10^-9..10^9` as diagnostic-only evidence. A failure in those bands must
    be classified as a threshold/conditioning boundary, never "fixed" by
    weakening phase-control or acceptance contracts.
- [ ] Add real-data layout-permutation metamorphic coverage for both TP-1907
  P,T and TP-1906 P,H. Reorder phase-qualified components only through typed
  layouts, then compare results by component identity rather than vector index.
- [x] Add a compact route/backend evidence matrix on reviewed fixtures:
  default RST versus legacy NR for a fixed P,T case, and nested/monolithic/Auto
  P,H where each route is applicable. Compare accepted physical state, not
  iteration count.
- [x] Extract a test-only accepted-solution assertion helper for finite values,
  conservation, complementarity, and file immutability. Keep scenario-specific
  topology and external assertions explicit rather than hiding them behind a
  generic fixture framework.
- [ ] Split `frozen_reference/mod.rs` by ownership (loader/catalog metadata,
  typed schemas, and shared assertions) before it grows further; preserve the
  current public module paths through re-exports.
- [ ] Enrich frozen metadata only from reviewed sources. In particular, do not
  invent an Argonne stable identifier/version or add byte hashes before a
  repository-wide line-ending policy exists.

### P9.5a - Synthetic lifecycle evidence hardening

The following is deliberately narrow. It closes evidence gaps identified by a
review of the completed I1-I3 suite; it must not create another P,H solver or
duplicate general phase-control regressions already owned by the runner.

- [x] Keep the fixed-topology I1/I2-to-canonical comparison separate from I3
  phase selection. It proves equal equations under one declared topology, not
  correct activation/deactivation.
- [x] Keep accepted-target continuation transactional. The lifecycle story
  injects a rejected target and proves that the next accepted solve starts from
  the preceding accepted phase set; no mutable inspection API is warranted only
  to expose an internal numerical seed.
- [x] Keep one in-band history scenario: the same continuous thermodynamic
  point retains a previously active candidate but does not create it from fresh
  inactive history. More synthetic in-band fixtures add no new contract unless
  they exercise a different route.
- [x] The fixed-topology comparator already proves complete evidence and
  localizes temperature, thermodynamic, composition, enthalpy, conservation,
  and identity disagreements.
- [x] Add one explicit comparator case with matching identities but unavailable
  canonical chemical or conservation evidence. Its corresponding axis must be
  `None` and `is_complete_match()` must be false, proving that `None` means
  unavailable/not applicable rather than implicit success.
- [x] Tighten the controlled synthetic disappearance story. Determine and
  document which route its constructed fixture actually follows. If it is the
  intended no-interior-root case, require exactly one
  `BoundaryUnstableActivePhase` transition rather than accepting either reason;
  retain the broader production contract because both routes are legitimate.
- [x] At the recovered synthetic gas-only boundary, add explicit I1-to-I3
  evidence: independently compute `ln(Q)-ln(K) > 0`, require canonical
  `TPD > 0`, and compare them with the reaction-normalized
  `TPD = R*T*(ln(Q)-ln(K))/nu_candidate` relation. Do not invent a positive
  I2 root after the candidate has physically disappeared.
- [x] Make simple topology expectations explicit: one activation for the
  favorable inactive fixture, zero transitions for stable absence, and one
  deactivation for controlled boundary disappearance. Do not impose those
  counts on continuation or hysteresis histories.
- [x] Complete the compact I1/I2 typed-failure matrix only where a distinct
  contract is still absent: infeasible extent interval and exhausted inner
  bisection budget. The existing invalid bracket/target, unbracketed outer and
  inner roots, non-finite thermochemistry, and exhausted outer budget are
  already sufficient; do not add duplicate error variants.

### P9 completion gate

- [x] The synthetic I1-I3 P,H milestone is complete: the deterministic scalar
  fixture recovers its constructed point, fixed-topology comparison and
  coordinate scaling agree, and the production lifecycle covers appearance,
  stable absence, boundary disappearance, accepted-only continuation, and
  history-dependent in-band retention.
- [x] P9 core pure-phase validation is complete through offline local-data I4:
  fixed topology, stable absence, appearance, disappearance, and a second
  carbon chemistry family all compare independent evidence with canonical
  accepted states. Non-NASA format diversity and optional real hysteresis
  discovery remain separately tracked evidence extensions, not blockers for
  this pure-phase milestone.

## Non-negotiable contracts

A nonlinear backend reporting convergence is not sufficient for accepting an
equilibrium result. Every accepted result must pass the same backend-independent
validation gate.

- [x] Require finite positive temperature and pressure and finite input data.
- [x] Require a dimensionally consistent species list, initial composition,
  elemental-composition matrix, phase assignment, and thermochemical vector.
- [x] Define the contract for zero and trace species in log space. Do not hide
  clipping, floors, or exponent saturation inside residual evaluation.
  The phase-control helpers now use a named trace-floor constant instead of a
  hidden numeric literal, and the log-mole seed path stays centralized in the
  typed initial-guess policy.
- [x] Reject NaN/Inf residuals, Jacobians, iterates, and reconstructed mole
  numbers with typed errors.
- [x] Check non-negative mole numbers, elemental conservation, normalized
  nonlinear residuals, and equilibrium optimality before publishing a result.
  - [x] The candidate report now exposes both L2 and RMS residual norms, so
    callers can reason about normalized residual evidence without recomputing
    it from scratch.
  - [x] Reaction-affinity optimality now has its own tolerance in the
    backend-independent acceptance gate, instead of being only reported
    passively after the solve.
  - [x] The validation gate now rejects candidates whose reaction-affinity
    block exceeds its tolerance, even when the backend-facing residual is
    small.
- [x] Keep species and reaction ordering explicit and deterministic. The
  validated problem and prepared snapshot now preserve canonical species,
  phase, and reaction ordering, and the golden fixtures plus solver-report
  tests lock that behavior down.
- [x] Never publish partial state from a failed setup, solve, fallback attempt,
  temperature sweep, or phase-control restart.
  - [x] Legacy matrix/system setup now validates the numerical settings before
    staging, so malformed tolerances or budgets cannot publish a partly rebuilt
    reaction basis.
  - [x] Single-temperature solve now stages automatic seeds and reconstructed
    mole/result bundles locally, publishing them only after an accepted solve.
  - [x] Sequential and parallel sweep paths now validate their complete
    column-oriented mole table before replacing any published sweep results.
  - [x] A failed solve after a previously accepted solve keeps the published
    solution, mole state, validation report, and solve report intact.

## P0 - Establish a safe canonical core

### P0.1 Characterize the current behavior

- [x] Inventory the public entry points in `Chem_eq_K_eq*`,
  `easy_equilibrium.rs`, and the retired classical equilibrium stack.
  A current behavior snapshot now lives in `behavior_inventory.md` and
  separates canonical modules, transitional modules, and their test families.
- [x] Record which current tests are physical regression fixtures, solver unit
  tests, debug experiments, or duplicates.
  Canonical regression and solver-policy tests live in the `equilibrium_*`
  modules, low-level numerical characterization lives in
  `equilibrium_nonlinear.rs`, legacy/classical coverage stays in the
  the retired classical stack and `NR_Legacy.rs` suites, while retained
  `debug_*` studies remain explicitly experimental. `Untitled-1.rs` has been
  removed.
  - [x] Audit the retired classical test family before deletion.
    The migration notes now separate adopted canonical contracts, P4
    phase-bridge scenarios, and tests tied only to obsolete state.
- [x] Preserve representative currently solvable systems as golden regression
  fixtures before changing numerical backends. The representative O2/O and
  N2/N fixtures now live in `equilibrium_golden_fixtures_tests.rs` and the
  RST matrix tests, so backend changes immediately trip a stable regression.
- [x] Record difficult and currently failing systems separately. A known
  failure must not be converted into a weak passing assertion.
  - [x] Hard fallback-prone fixtures are listed in `behavior_inventory.md` so
    future refactors can keep them distinct from easy acceptance regressions.
  - [x] No canonical equilibrium failure is being softened into a passing
    assertion; retained hard cases remain explicit regression fixtures.
- [x] Identify all `panic!`, `assert!`, `expect`, unchecked indexing, and
  production `unwrap` paths in setup, solving, phase control, and output.
  The next audit pass should focus on this boundary now that the public
  surface inventory is explicit.
  Current hotspots are concentrated in the retained legacy/classical modules;
  the canonical `equilibrium_*` path is already much more typed and its
  remaining `unwrap` / `assert` hits are in tests or legacy-only code.

### P0.2 Typed domain boundary

- [x] Introduce a validated `EquilibriumProblem` containing species identity,
  conditions, initial composition, elemental constraints, thermochemical data,
  phase/activity information, and numerical formulation options.
- [x] Introduce an immutable `PreparedEquilibriumProblem` for matrices,
  reaction bases, element totals, and phase-derived data.
- [x] Add residual scaling to the prepared snapshot as an explicit typed
  contract once its units and backend contract are defined. The prepared
  boundary now owns `ResidualScalingContract`; variable scaling remains a
  separate future design because it is a different coordinate transform.
- [x] Introduce `EquilibriumSolution`; accepted legacy solves expose an
  immutable snapshot instead of requiring callers to combine mutable fields.
- [x] Introduce `EquilibriumSolveReport` with policy, ordered attempts, and
  accepted backend details. RST attempts additionally expose termination,
  iteration/evaluation counters, linear solves, and elapsed time; comparable
  legacy counters remain unavailable until that temporary adapter is replaced.
  The report now also exposes summary helpers so UI/tests can reason about
  fallbacks and skipped attempts without hand-walking the raw attempt vector.
- [x] Introduce typed identifiers for species, elements, phases, and reactions
  where raw indices can be confused. Typed boundary accessors now exist on the
  validated problem and prepared snapshot; the internal ordered arrays are
  still the storage format for now.
- [x] Separate domain settings from backend-specific solver settings.
  `EquilibriumLogMoles` now owns a dedicated `EquilibriumSolverSettings`
  bundle, while the domain/problem boundary stays focused on chemistry and
  composition data.
- [x] Replace public mutation sequences such as `set_problem -> create_* ->
  solve` with validated construction plus one explicit solve operation.
  - [x] The canonical `solve_problem_with` configuration hook now receives
    only `EquilibriumSolverSettings`, so post-validation tuning cannot mutate
    chemistry, phase layout, conditions, or published state.
  - [x] Public workflow constructors that still need a mutable sweep engine now
    route through `EquilibriumProblem` + `EquilibriumLogMoles::from_problem`
    instead of rebuilding the solver by staging mutable fields step by step.

### P0.3 Typed errors

- [x] Split `ReactionExtentError` into meaningful categories: invalid input,
  thermochemical lookup, formulation, residual/Jacobian evaluation, backend
  failure, non-convergence, invalid candidate solution, validation mismatch,
  unsupported model, and all-backends-failed. RST linear solve, singular
  Jacobian, and numerical-breakdown errors now have a typed backend-failure
  category instead of being mislabelled as legacy evaluation errors. Canonical
  Gibbs-cache, residual-generator dimension, and condition failures no longer
  use the catch-all `Other` variant. The `Other` variant has been removed from
  the public error enum; the remaining work is to make every category carry
  enough domain identity for a user to locate the bad input.
  - [x] `ReactionExtentErrorKind` now exposes a machine-readable top-level
    classification for the requested failure families, so tests and policy code
    can branch on categories instead of parsing display text.
- [x] Preserve error sources instead of converting them to debug strings in
  every backend adapter and legacy compatibility path.
  - [x] `ReactionExtentError` and `SolveError` now implement stable `Display`
    and `std::error::Error` contracts. Nested `SubsDataError`/`SolveError`
    values remain available through `source()`.
  - [x] Solver-cascade attempt reports now retain stable `Display` diagnostics
    rather than private `Debug` formatting.
  - [x] Recoverable attempt failures carry a typed report category so callers
    do not need to infer retry diagnostics by parsing a message.
- [x] Include species/reaction identifiers and offending values in diagnostics.
  - [x] Canonical residual/Jacobian condition checks now identify the exact
    invalid scalar; non-finite Gibbs values identify the species index and
    temperature; invalid phase totals identify the phase index and value.
  - [x] Invalid equilibrium candidates now carry the offending field name in
    the error itself, so setup/acceptance failures are no longer just opaque
    "candidate rejected" messages.
- [x] Make invalid input and unsupported physics non-retryable. The retry
  boundary now has typed helpers on `ReactionExtentError`, and the cascade
  helper explicitly classifies only numerical backend failures as retryable.
- [x] Return every backend attempt in `AllBackendsFailed`, rather than only the
  final hand-written solver error.
  - [x] A cascade that reaches a non-retryable error after earlier attempts now
    returns `CascadeAborted { attempts, cause }`; ordinary exhaustion remains
    `AllBackendsFailed`, so neither case loses its diagnostic trace.

### P0.4 One logarithmic formulation

- [x] Reject non-square systems before closure construction: independent
  reaction equations plus elemental-balance equations must equal the number of
  log-mole unknowns.
- [x] Move residual, Jacobian, scaling, and mole reconstruction into pure
  formulation functions over `PreparedEquilibriumProblem`. The canonical
  residual and Jacobian now have pure `evaluate_*` entry points, while legacy
  closure/wrapper APIs delegate to the same formulas. Mole reconstruction,
  immutable accepted-result packaging, and explicit row scaling now also live
  on the prepared boundary. Variable scaling remains a separate coordinate
  transformation and is intentionally not implied by row scaling.
- [x] Ensure every numerical backend receives exactly the same residual and
  Jacobian contract.
  - [x] The offline formulation fixture now checks every analytic Jacobian
    entry against the derivative of the corresponding symbolic residual.
- [x] Remove duplicate closure construction from single solve, parallel solve,
  temperature sweep, and phase-control paths.
  - [x] Single-solve and parallel-solve Gibbs closure construction now share
    one ordered lookup helper, so the closure ordering contract is no longer
    duplicated in two separate code paths.
  - [x] Temperature-sweep and workflow constructors now reuse the same shared
    Gibbs-function assembly path instead of hand-building another closure
    ordering layer.
- [x] Validate the analytic Jacobian against finite differences across normal,
  trace-species, high-temperature, and low-temperature states. The regression
  matrix now covers several fixed temperature/state pairs on the canonical
  typed problem; richer real thermochemistry temperature fixtures remain a
  separate future layer.
- [x] Document the meaning and units of every residual block and scale factor.
- [x] Guard `exp(log_n)` reconstruction against overflow and underflow without
  silently accepting a distorted solution.

### P0.5 Backend-independent acceptance gate

- [x] Add a pure `validate_equilibrium_candidate` function.
  The acceptance gate now takes typed residual data plus typed acceptance
  criteria, rather than a bundle of unstructured vectors and bare tolerances.
- [x] Report absolute and scaled residual norms separately.
- [x] Report maximum elemental-balance error by element.
- [x] Report minimum mole number and all non-finite values.
- [x] Reject truncated or empty raw/acceptance residual vectors before their
  norms can be interpreted as convergence evidence.
- [x] Check reaction affinity or chemical-potential optimality appropriate to
  the implemented model.
  The acceptance gate now reports the reaction-affinity block separately from
  the elemental-balance block, so the physical optimality evidence is explicit
  and can be compared without reusing the wrong residual rows.
- [x] Distinguish `BackendConverged` from `SolutionAccepted` in reports.
- [x] Publish solution/state only after this gate succeeds.

### P0.6 Isolate legacy nonlinear solvers

- [x] Move `LMSolver`, `NRSolver`, and `TrustRegionSolver` behind one temporary
  legacy adapter. `equilibrium_legacy_backend` now owns their mechanical
  dispatch; policy, budgets, validation, retry, and publication remain in the
  canonical orchestration path.
- [x] Remove panics, assertions on user data, global logger initialization, and
  direct console output from legacy production paths. The canonical log-moles
  engine no longer initializes a process-global logger, and the legacy trust
  region solver now exposes a singular LU step as `SolveError::SingularMatrix`
  instead of silently substituting a zero step. The public Newton step-limit
  helper now returns typed input errors instead of asserting vector lengths.
  Remaining console output and legacy solver diagnostics still require
  isolation.
  - [x] Canonical sweep-table construction is now pure (`moles_table`); the
    solver and workflow APIs no longer print directly to stdout.
  - [x] Correct the Trust Region quadratic-model reduction to use `J^T f` for
    its linear term. The previous `f^T p` expression was mathematically wrong
    and caused ill-conditioned systems to shrink the trust radius to machine
    noise; the regression fixture now converges deterministically.
  - [x] Canonical `equilibrium_*` runtime paths are now free of production
    `panic!/unwrap/expect` sites; the remaining occurrences are confined to
    tests, doc examples, or the legacy/classical layer that is explicitly out
    of the canonical path.
- [x] Do not select legacy solvers implicitly outside an explicit solver policy.
  The default no longer falls back to the handwritten cascade; legacy
  backends remain available only through an explicit policy.
- [x] Add characterization tests before changing their behavior.
  Canonical regression, backend-matrix, and fallback-story fixtures now pin
  the current behavior before more structural changes land.
- [x] Define the legacy-backend retention rule: the handwritten backend stays
  as an explicit compatibility/fallback path, even after the RST matrix covers
  the main regression set. The remaining work is to keep its policy boundary
  strict and its diagnostics typed.
  - [x] Review the retained handwritten solvers as production-quality
    fallbacks, not as disposable demo code. `LMSolver`, `NRSolver`, and
    `TrustRegionSolver` should keep explicit step control, bounded failure
    modes, and deterministic test coverage even though they are legacy.
  - [x] Remove `panic!` / `assert!` / silent `None` publication paths from
    `NR_Legacy.rs` where a typed failure or an explicit best-effort step can
    be returned instead.
  - [x] Decide whether the two currently failing legacy regression tests are
    strict contract tests or characterization tests. If the solver is allowed
    to stop without convergence for a given fixture, encode that explicitly in
    the test instead of relying on an implicit panic or unwrap.
  - [x] Audit the inner step control in `LMSolver`, `NRSolver`, and
    `TrustRegionSolver` for finite-value checks, dimension checks, and clear
    iteration-exhaustion reports.

## P1 - Reliable solver architecture

### P1.1 Common backend adapter

- [x] Define an internal `EquilibriumNonlinearBackend` interface receiving a
  prepared problem, initial iterate, residual, Jacobian, and attempt budget.
  The adapter now owns one mechanical request/result boundary for both legacy
  and RST backends.
- [x] Return a typed `SolverAttemptReport` with backend name, termination
  reason, iteration/evaluation counts, norms, elapsed time, and candidate.
- [x] Keep RustedSciThe types inside adapter modules so domain types do not
  depend on one numerical library's API. The prepared symbolic problem is now
  carried as an opaque adapter-level wrapper instead of leaking the concrete
  RST problem type into the domain layer.
- [x] Extend the typed problem boundary with immutable symbolic thermochemical
  expressions (or a dedicated symbolic formulation snapshot) before making an
  RST-symbolic policy the default for every `EquilibriumProblem`. The current
  typed input intentionally stores numeric Gibbs closures only. The adapter
  now builds a dedicated immutable symbolic thermochemistry snapshot before it
  assembles the RST problem, so the lookup pipeline is no longer intertwined
  with the symbolic bridge itself.
- [x] Support deterministic fake backends for cascade tests.

### P1.2 RustedSciThe migration

- [x] Audit the current RustedSciThe 0.4.12 `Nonlinear_systems` API through its
  `NonlinearProblem`, `JacobianProvider`, `SolverEngine`, and typed result APIs.
- [x] Implement symbolic RST adapters for ordinary LM, MINPACK LM, Nielsen LM,
  trust-region LM, Powell dogleg, and damped Newton. The canonical adapter
  passes residual `Expr` plus temperature/settings to
  `SymbolicNonlinearProblem`; RST owns Lambdify/Jacobian preparation.
- [x] Benchmark and characterize candidates before fixing the default order.
  Do not choose a default solely from method names or one easy fixture.
  - [x] Current RST default order is now characterized by a dedicated matrix
    test, so the chosen order is explicit and regression-protected even before
    any future performance-driven reordering.
- [x] Map RST termination reasons without losing details.
- [x] Verify scaling, bounds, and stopping tolerances have the same meaning at
  the domain and backend boundaries.
  - [x] RST solve-contract construction is now centralized as a typed
    `RustedSciTheSolveContract`, so tolerance and iteration budget are no
    longer scattered across call sites. Feasibility bounds remain a
    formulation-side concern and are not silently implied by the backend
    contract.
- [x] Use the RST symbolic Lambdify backend for the first migration; keep AOT
  out until a measurement shows a material evaluation bottleneck.

### P1.3 Explicit solver cascade

- [x] Replace the loop inside `solver_impl` with a typed `SolverPolicy`.
- [x] Support at least `Single(backend)` and `Cascade(Vec<backend>)`; provide a
  documented `Auto` policy only after the backend matrix is measured.
- [x] Give every attempt an explicit iteration budget and enforce a global
  cascade budget. `SolverCascadeBudget` caps started backends, per-backend
  iterations, and total allocated iterations; unstarted backends are retained
  as ordered `Skipped` entries in the report. Evaluation budgets still depend
  on support from the numerical backend API.
- [x] Restart each backend from the original validated initial guess unless a
  named warm-start policy explicitly permits reuse of a previous candidate.
  The cascade now documents and tests this invariant; no implicit candidate
  reuse is allowed.
- [x] Retry only recoverable numerical failures such as non-convergence,
  singular steps, or rejected candidate solutions.
- [x] Never retry invalid input, missing thermochemical data, dimension errors,
  or unsupported activity/phase models.
- [x] Run the common acceptance gate after every backend success. If the
  candidate is invalid, record why and continue the cascade.
- [x] Preserve deterministic backend order and expose all attempts in the final
  solve report.
- [x] If several valid candidates are retained, choose between them only by a
  documented physical/numerical criterion, never by "first finite vector".
  A dedicated comparator now ranks accepted candidates by residual norm,
  raw residual norm, balance error, and minimum mole evidence. The current
  cascade still stops on the first accepted backend, so this criterion is now
  explicit and testable for any future multi-candidate retention path.
- [x] Make logging a view of the report, not the only record of fallback.
  - [x] Solve reports now expose typed summary accessors for accepted,
    fallback, skipped, and ordered attempt views.
  - [x] The main solve path now logs the typed solve report summary instead of
    a separate ad hoc fallback message.

### P1.4 Initialization, scaling, and continuation

- [x] Centralize construction and validation of log-mole initial guesses. The
  typed constructor, solver settings, and legacy mutable setter now share the
  same finite-value and dimension validation; temperature-sweep/phase-control
  seed paths now reuse the same trace-seed policy instead of re-encoding the
  floor contract in multiple places.
- [x] Define a deterministic trace-species floor contract. A named absolute
  coordinate default now replaces hidden duplicated literals, and a typed
  system-scale-relative policy is available when callers want scale-aware
  seeding explicitly. The production default remains the absolute floor until
  we finish the numerical characterization for a relative default.
- [x] Test row scaling and variable scaling independently. Row scaling now has
  typed dimension/scale validation and direct residual/Jacobian tests.
  Variable scaling now has its own typed coordinate-scaling contract and
  tests, but it remains intentionally separate from the solver path; the
  remaining work is integration, not conceptual separation.
- [x] Treat temperature continuation as an explicit policy with transactional
  per-temperature results and a clear warm-start contract. Sequential and
  chunked sweeps now use `ContinuationSeedPolicy`: either reuse the previous
  accepted solution or restart every point from the configured seed. The fully
  parallel `par2` path is explicitly independent because no ordered prior
  point exists while tasks run concurrently. All sweep paths publish
  `TemperatureSolveSnapshot` evidence for accepted points.
- [x] Make legacy mutable `set_problem` and explicit-mole-map updates validate
  before publishing state. They now reject invalid pressure, non-finite or
  negative amounts, out-of-range/duplicate phase indices, and unknown mapped
  substances without clearing a previously published solution.
- [x] On continuation failure, record the failed temperature and attempts;
  never silently leave vectors with different lengths or stale state. Sweeps
  publish `TemperatureSolveFailure` records with the temperature, diagnostic
  message, and cascade trace when available, alongside separate accepted-point
  snapshots.
- [x] Do not add multi-start without regression evidence justifying its
  complexity. The deterministic backend cascade remains the production policy;
  multi-start is an explicitly deferred extension, not a release blocker.

## P2 - Independent equilibrium-constant validation

The K-equilibrium validator is valuable precisely because it is a second
formulation. It may share immutable thermochemical input and common final
invariant checks, but it must not reuse the main residual/Jacobian assembly.

### P2.1 Reaction-basis contract

- [x] Make `ReactionBasis` typed and validate its dimensions, rank, nullspace
  accuracy, species ordering, and elemental conservation.
  - The independent validator now owns `ValidatedReactionBasis`, which binds
    the species order to a finite species-by-reaction matrix and validates
    dimensions, declared rank, reaction count, and `A^T * N`. Replacing the
    canonical solver's older public `ReactionBasis` remains separate work.
- [x] Detect underdetermined, overdetermined, and numerically rank-ambiguous
  bases with typed errors.
- [x] Make basis normalization/sign conventions deterministic so comparisons
  and fixtures are stable.
  - Explicit validator bases now normalize each reaction by its first
    significant species coefficient and orient that coefficient as a
    reactant. Multi-dimensional SVD basis canonicalization is still open.
- [x] Test basis invariance under species and element permutations.

### P2.2 Independent K_eq solver

- [x] Implement a distinct reaction-extent/K_eq problem and solver path for
  small systems with a known independent reaction basis.
  The first pass now solves one-reaction systems in extent space with a
  safeguarded Newton/bracketing loop and returns its own typed solve report
  plus an independent validation report.
- [x] Compute `ln(K)` and reaction quotients through code independent from the
  main log-moles residual builder.
  `EquilibriumConstantProblem` evaluates standard reaction Gibbs energies,
  dimensionless ideal-gas activities, and `ln(Q) - ln(K)` without calling the
  canonical residual/Jacobian implementation.
- [x] Define supported activity models and phase combinations explicitly.
  The first contract supports one ideal-gas phase only; condensed and mixed
  phases must become explicit variants rather than implicit corrections.
- [x] Return `ValidationNotApplicable` for unsupported or excessively large
  systems; do not report a false solver failure.
  The first-pass solver refuses multi-reaction systems explicitly instead of
  pretending that they failed numerically.
- [x] Give the validation solver its own typed report and numerical tolerances.
  Candidate assessment now returns ordered per-reaction `ln(Q)`, `ln(K)`, and
  residual evidence plus the maximum residual and acceptance decision. The
  future extent solver will add iteration/backend evidence to this report.
- [x] Keep it optional in production (`ValidationMode::Off/WhenApplicable/
  Required`) and enabled broadly in tests.
  The new solver already exposes an explicit on/off/applicable policy, and
  `solve_if_applicable` now returns `None` immediately in `Off` mode while the
  solver tests cover applicable, `Off`, and `Required`-mode paths. The
  top-level solve boundary now threads the validation mode and publishes an
  optional independent K_eq status alongside the main accepted solution.

### P2.3 Cross-validation report

- [x] Compare main and K_eq solutions by species moles/mole fractions with
  absolute, relative, and trace-species-aware tolerances.
  The new cross-validation report compares deterministic species moles and
  mole fractions and records the largest disagreement.
- [x] Compare elemental balances independently for both solutions.
  Canonical balance evidence remains in the canonical candidate report, while
  the independent validator remains logically separate.
- [x] Evaluate `ln(Q) - ln(K)` for every independent reaction.
  The independent validation report carries the per-reaction residuals, and
  the comparison layer exposes the largest one as part of the bridge report.
- [x] Compare total Gibbs energy or another documented equilibrium objective.
  The bridge report now computes a total ideal-gas Gibbs comparison directly
  from the thermochemical inputs and mole numbers.
- [x] Report the largest disagreement with species/reaction identity.
  Each comparison row carries the species name, index, canonical value, K_eq
  value, and both absolute deltas.
- [x] Distinguish main-solver failure, validator failure, not-applicable, and
  genuine cross-validation mismatch.
  The report layer is in place, and a typed status helper now classifies
  canonical failure, validator failure, not-applicable, and compared-result
  branches separately. The top-level solve path now routes and publishes the
  resulting status instead of collapsing those cases into one branch.

### P2.4 Validation fixture matrix

- [x] Add small analytically tractable dissociation/association systems.
  - Offline contracts now cover `A2 <=> 2A` with pressure dependence and
    `2NO <=> N2 + O2` with the closed-form extent inherited from the old test
    module, now exercised through the independent extent solver.
    The inherited `2N2O <=> 2N2 + O2` and `2NO2 <=> N2 + 2O2` cases now also
    cover three-species reconstruction and the nonzero pressure exponent.
    Thermochemical-library association fixtures remain open.
- [x] Add water-gas and hydrogen/oxygen equilibrium fixtures.
  - The cross-validation matrix now includes a water-gas-shift fixture in
    addition to the earlier diatomic and nitric-oxide contracts.
- [x] Add methane/air lean, stoichiometric, and rich fixtures.
  - The methane-combustion baseline, lean, and rich cases are now covered in
    the independent validation matrix.
- [x] Add inert dilution and pressure-variation fixtures.
  - A hydrogen/oxygen fixture with inert `N2` now exercises both pressure
    dependence and diluent handling in the independent validator bridge.
- [x] Add low-, medium-, and high-temperature sweeps.
  - The methane-combustion cross-validation matrix now exercises a
    low/medium/high temperature sweep at fixed composition.
- [x] Compare the canonical solver and K_eq validator for every applicable
    fixture, not only for one final scalar.
- [x] Keep all validation fixtures offline and deterministic.

## P3 - Test and reliability matrix

### P3.1 Formulation unit tests

- [x] Test residual blocks, Jacobian blocks, scaling, log/mole conversion, and
  phase totals independently.
- [x] Test analytic versus finite-difference Jacobians over a state matrix.
- [x] Test NaN/Inf, overflow, underflow, zero totals, duplicate species,
  unassigned species, singular composition matrices, and malformed phases.
  The public phase-index boundary now returns a dimension error instead of
  indexing past the phase map; broaden the malformed-phase matrix further.
- [x] Test exact dimensions and ordering at every boundary.

### P3.2 Backend matrix

- [x] Run every supported RST backend against the same curated fixture set.
  The first `O2/O` case is covered, high-temperature `N2/N` provides a real
  rejected-candidate -> RST fallback story, and the dilute inert-diluent
  matrix now exercises the same accepted-backend reporting contract.
- [x] Record convergence, acceptance, iterations, residuals, balance errors,
  and runtime without making runtime assertions flaky.
- [x] Keep legacy backend parity tests as characterization while the handwritten
  LM/NR/TR implementations remain supported fallback backends. Their deletion
  criteria are intentionally deferred; this does not make them a second
  production orchestration path.
- [x] Include difficult cases that require fallback, not only easy systems on
  which every solver succeeds immediately.

### P3.3 Cascade story tests

- [x] First backend fails recoverably; second succeeds and is accepted.
- [x] First backend returns a nominal success with an invalid candidate; the
  acceptance gate rejects it and the next backend succeeds.
- [x] Invalid input produces zero backend attempts.
- [x] Every backend fails; the error contains every attempt in exact order.
- [x] Global budget exhaustion stops the cascade deterministically.
- [x] A failed attempt cannot mutate the original problem, initial guess,
  published solution, or temperature-sweep results. The RST rejected-candidate
  fallback story now verifies that the original seed and initial composition
  remain unchanged; failed numerical attempts and temperature sweeps remain.
- [x] Warm-start behavior occurs only when explicitly selected.

### P3.4 Physical and metamorphic tests

- [x] Check non-negative concentrations and elemental conservation for all
  accepted fixtures.
- [x] Check invariance under species, element, and reaction reordering.
- [x] Check consistent results when all initial mole totals are scaled.
- [x] Check robustness to reasonable initial-guess perturbations.
- [x] Check inert-species addition and removal where the model predicts it.
- [x] Check temperature/pressure trends against known qualitative behavior.

### P3.5 Integration and regression tests

  - [x] Build problems from real offline `SubsData` thermochemistry and retain
    provenance in reports.
    Local `SubsData` thermochemistry is the normal path; NIST fallback should
    remain a rare and explicitly reported fallback path, not the default.
  - [x] Cover mixed-library thermochemical inputs without network access.
- [x] Add transactional tests for setup, solve, phase-control restarts, and
  serial/parallel temperature sweeps.
- [x] Replace debug binaries/scripts with assertions in named test modules or
  documented examples.
  The standalone debug comparison file has been removed, and no separate
  debug study source files remain in this subsystem.
- [x] Separate fast unit tests, deterministic integration tests, and explicitly
  ignored expensive studies.

### P3.6 Retained ideas from the classical implementation

The retired classical equilibrium modules are not a compatibility target. The
following list records the useful contracts discovered during their deletion
audit so the implementation can be removed without losing sound ideas.

Technical work that is independent from the future phase bridge:

- [x] Add a pure typed `SpeciesCapacityReport` derived from elemental totals:
  `n_i_max = min(b_e / a_ie)` over elements present in species `i`.
  - Do not copy the legacy `1.2` safety multiplier or string-keyed solver
    bounds.
  - Use the exact capacity for candidate diagnostics, initial-guess checks,
    impossible-composition detection, and composition-ordering tests.
  - Keep log-moles positivity as the production feasibility mechanism; this
    report is evidence and validation, not a second hidden solver policy.
- [x] `EquilibriumProblem::validate` now rejects species rows without any
  positive elemental support, so impossible compositions are caught before
  preparation.
- [x] Add optional, pure `FormulationDiagnostics` for prepared problems.
  Include numerical rank, singular values, a documented condition estimate,
  and suspicious null directions where available.
  - Compute expensive SVD diagnostics only during preparation or by explicit
    request, never on every nonlinear iteration.
  - Return typed data for reports/tests; do not print directly or let a
    diagnostic heuristic silently reject an otherwise valid solution.
  - The prepared boundary now exposes `formulation_diagnostics(tolerance)` and
    `preview_with_diagnostics(tolerance)` without changing solve behavior.
- [x] Add a typed `EquilibriumProblemPreview` (or equivalent report rows) for
  CLI, GUI, examples, and snapshot tests. It should expose conditions, species
  and element ordering, element matrix, reaction basis, initial inventory,
  scaling contract, solver policy, thermochemical provenance, and optional
  formulation diagnostics without owning any solving behavior.
  - The current preview exposes the validated problem state, reaction basis
    size, exact element totals, and the species-capacity evidence that can be
    derived without the future phase bridge.
  - The preview also now exposes stable summary rows and `Display` output for
    CLI/snapshot consumers without adding any solving behavior.
  - `SpeciesCapacityReport` and `FormulationDiagnostics` now also have stable
    `Display` output for human-readable reports and regression snapshots.
- [x] Give the independent K_eq cross-validation report the same typed
  summary/`Display` treatment so comparisons can be inspected without ad hoc
  string assembly.

Ideas intentionally deferred until P4:

- [x] Validate named initial compositions transactionally at the
  `ResolvedPhaseSystem -> EquilibriumProblem` adapter: every component must be
  known, represented exactly once by its typed component id, and aligned with
  deterministic solver ordering. `new_with_sparse_initial_composition` now
  enforces this contract. Do not restore the legacy
  `HashMap<Option<String>, ...>` contract.
- [x] Defer a small independent Lagrange-stationarity validator unless
  multi-reaction regression evidence shows a real validation gap after the
  phase bridge exists.
  - Its possible value is mathematical independence from the reaction-basis
    residual and applicability beyond the one-reaction extent validator.
  - Do not migrate the legacy mutable `Solver`, auxiliary `Np` unknowns, raw
    mole variables, or unscaled equations. Any future experiment must use
    log-moles, typed phase/component ids, scaled residuals, and a typed report.

Classical ideas already superseded by the canonical engine:

- [x] SVD reaction-basis discovery and elemental-conservation validation.
- [x] Temperature-band discovery and thermochemical coefficient refresh.
- [x] Explicit sequential warm-start/continuation policy.
- [x] Ordered nonlinear backend cascade, acceptance gate, budgets, and attempt
  reports.
- [x] Independent reaction `delta G`, `ln(K)`, quotient, and extent
  reconstruction through the K_eq validation subsystem.

Classical behavior that must not be migrated:

- [x] Global Lagrange/Newton interpolation followed by value clamping. Global
  polynomial oscillation and clamping can hide invalid composition and break
  elemental balances; future plotting interpolation belongs to a separate
  postprocessing layer with explicit error bounds.
- [x] Parallel mutable maps, placeholder zero closures during clone, string
  variable names such as `N0`/`Lambda0`, direct `println!`, and unchecked
  `unwrap`-driven output.
- [x] Auxiliary phase-total unknowns when phase totals can be derived from the
  ordered species mole vector.
- [x] The classical solver cascade and residual-norm-only acceptance rule.

## P4 - Unify the two equilibrium worlds

This phase begins only after P0-P3 contracts are stable.

The first production target is a closed reacting system at fixed pressure and
temperature. The supported phase-model slice is deliberately narrow:

- one ideal-gas phase containing any number of gas components;
- any number of one-component pure liquid/solid/condensed phases;
- no ideal/non-ideal solution phase and no multiple gas phases until those
  activity models have explicit equations and validation fixtures.

This is not merely an adapter task. `ResolvedPhaseSystem` already has the
correct phase-qualified component identity, while `EquilibriumProblem` still
uses unique bare strings and integer-only phase records. The domain boundary
must be corrected before the two systems are connected.

### Current production top-level priorities

The canonical fixed-`P,T` equations, backend cascade, independent K_eq
validation, bounded phase control, immutable result, and resolved-phase bridge
are already implemented. The remaining work must now be judged by whether it
makes `solve_resolved_pt` the single reliable production entry point.

#### Production-readiness audit (2026-07-26)

**Verdict:** the supported fixed-`P,T` slice is a strong production candidate,
but it is not yet the only internally canonical and externally unambiguous
workflow. The first production milestone remains deliberately limited to one
ideal-gas phase plus any number of one-component pure condensed phases.
Multiple gas phases and condensed solution models are separate physical-model
features, not defects in this milestone.

Release blockers, in execution order:

1. [x] Remove `EquilibriumLogMoles` as the mutable orchestration host from the
   canonical `solve_resolved_pt` path. Both fixed and bounded modes now use
   immutable prepared inner problems. `PreparedPhaseControlRunner` owns only
   active-set transitions, restart seeds, and reports; the historical
   `solve_with_phase_control` remains available solely for compatibility
   callers and is no longer reached by the resolved-data facade.
2. [x] Honor independent-validation policy in bounded phase control.
   `solve_fixed_active_set_candidate` currently forces
   `keq_validation_mode = Off`, and `from_phase_control_parts` publishes
   `keq_validation_status = None`. In particular, `Required` must never be
   silently downgraded. Either validate the final applicable active set or
   return an explicit typed `NotApplicable`/error according to the requested
   mode.
3. [x] Separate physical published composition from internal positive
   log-coordinate trace floors. An inactive phase currently remains visible
   through `component_moles()` and `phase_total()` as a tiny positive amount
   while its status is `Inactive`. Define and test one unambiguous result
   contract: physical public amounts are zero for inactive phases, while
   numerical trace coordinates remain available only through explicit
   `numerical_*` diagnostic accessors. Bounded mixed-phase stories now cover
   excluded phase publication and trace-floor separation.
4. [x] Finish the independent typed policy boundary. `EquilibriumSolveOptions`
   and `PhaseControlPolicy` now expose the supported configuration through
   validated typed builders. There is one source of truth for trace policy;
   `with_solver_backend` clears an explicit conflicting policy; scalar phase
   limits, explicit-set duplicates, and resolved phase-index bounds are
   rejected before the bounded loop. Raw settings/manager constructors remain
   crate-private migration ingress only, so the public prelude does not leak
   mutable backend orchestration objects.
5. [x] Add a phase-qualified sparse/named initial-composition constructor to
   `PhaseEquilibriumPipelineRequest`. The lower layer already provides
   `MultiphaseInitialComposition::from_sparse`; the top-level pipeline should
   not require users to supply an anonymous `Vec<f64>` in an ordering that is
   only known after resolution. `new_with_sparse_initial_composition` now
   provides this contract and has dense-equivalence and duplicate-entry tests.

Mandatory reliability evidence before the production label:

- [x] Record live offline continued-lifecycle cycle detection as a deferred
  evidence gap. The
  supported fixed-P,T phase models do not naturally generate a physical cycle;
  add this only when a real cyclic fixture exists rather than manufacturing
  one by corrupting the phase policy. The existing synthetic cycle tests remain
  the sharper unit layer.
  - [x] Real disappearance, hysteresis retention, budget termination, and
    failed-transition rollback are covered by `equilibrium_live_data_tests`.
  The zero-inventory `AllCandidatePhases` boundary is now normalized before
  the positive log-moles solve and covered by a live water regression. A
  separate positive-inventory boundary-recovery scenario is also covered: the
  live high-temperature water fixture deactivates liquid water through the
  typed `BoundaryUnstableActivePhase` transition instead of publishing a
  backend NaN/`AllBackendsFailed` error.
- [x] Add inventory-scale metamorphic tests over several orders of magnitude
  on the public resolved-phase facade. The acceptance gate now applies the
  explicit contract `error <= absolute_tolerance +
  relative_tolerance * abs(original_element_total)` and the live H/O fixture
  exercises both small and large inventory scales.
- [x] Add explicit tests for bounded `K_eq` modes (`Off`, `WhenApplicable`,
  `Required`) and for the typed solver/phase-policy precedence rules covered
  by the public facade. The default `Off` path is now asserted not to publish
  a validator report; `WhenApplicable` publishes either comparison evidence or
  typed `ValidatorNotApplicable`, while `Required` rejects an unsupported
  mixed-phase request.
- [x] Make retained fallback backends library-quiet by default. The reachable
  legacy NR and shared LM/NR/TR fallback internals now route iteration-level
  diagnostics through `debug`; direct production `println!` calls were removed.
  Test-only demonstrations may still print their intermediate values.
- [x] Run the ignored 20/50/100/200-species benchmark after the immutable
  orchestration path is final. The baseline is recorded in
  `PHASE_CONTROL_BENCHMARK.md`. The original debug profile grew from `2.14 ms`
  at 20 species to `8.81 s` at 200 species; a later release-build run measured
  `67.2 us` and `70.5061 ms` respectively. Both profiles are retained as
  machine/build-specific characterization, and the release run is the useful
  deployment baseline. Cache/allocation work still requires a realistic
  transition benchmark.

Required migration and documentation closure:

- [x] Migrate `chem_equilibrium_gas_example.rs` and the legacy gas section of
  `chem_equilibrium_guides.md` to the typed resolved-phase pipeline facade.
  Retained internal consumers still need a separate migration audit.
- [x] Mark `easy_equilibrium` and broad mutable workflow entry points as
  compatibility/experimental APIs. `EasyEquilibrium`, the `gas_solver*`
  helpers, and `solve_with_phase_control` now carry Rust deprecation messages. Their
  unchecked indexing, `unwrap`, and direct console output remain outside the
  production prelude; the APIs stay available for characterization and
  migration.
- [x] Update architecture documents after orchestration migration so the
  dependency graph names exactly one production entry point and clearly
  distinguishes retained numerical fallback implementations from legacy
  orchestration. The fixed-P,T boundary and compatibility status are recorded
  in `ARCHITECTURE_RU.md`.

The following do **not** block the first production milestone: GUI work,
temperature-range orchestration, multiple ideal-gas phases, condensed solution
activity models, projection caching, speculative multi-start, and deletion of
the handwritten LM/NR/TR fallback implementations.

#### P0 - Remove the remaining canonical-path correctness risks

- [x] Replace the internal mutable `EquilibriumLogMoles` host used by
  `PhaseEquilibriumProblemBundle::solve_with_bounded_phase_control` with
  canonical orchestration over immutable `EquilibriumProblem` /
  `PreparedEquilibriumProblem`, active-set projections, backend policies, and
  acceptance reports. Retain the handwritten legacy nonlinear backends only as
  explicit fallback implementations behind the common backend contract.
- [x] Build an independent element-constraint basis for every active-set
  projection. A physically admissible reduced phase set must not be rejected
  merely because the full declared element matrix contains columns that become
  dependent in that active set. Preserve the mapping back to full element
  labels and validate conservation against the original full totals. The active
  set now carries the retained independent element rank explicitly, and the
  projection accepts physically valid rank-deficient cases such as single
  `H2O(g)` H/O inventories and permuted element columns.
- [x] Add regression tests for rank-deficient but physically valid active sets,
  including the one-component `H2O(g)` H/O case, element permutations, and
  transactional failure when the reduced basis is genuinely infeasible. The
  active-set test matrix now covers acceptance of valid rank-deficient
  projections, element-order permutations, and infeasible total vectors.

#### P1 - Finish the production request and data boundary

- [x] Replace the broad legacy-facing `EquilibriumSolverSettings` and mutable
  public-field `PhaseManager` at the facade boundary with small validated
  `EquilibriumSolveOptions` and `PhaseControlPolicy` types. Backend order,
  global budget, acceptance tolerances, K_eq validation, trace policy, and
  phase hysteresis must each have one source of truth.
  - The typed wrappers now exist and are wired into the high-level
    pipeline facade. Legacy setters remain as a migration shim, but the public
    request objects no longer have to expose the raw backend controller types
    as their only shape. The old raw-settings setters are now marked
    deprecated so the typed path is visually and semantically primary, and the
    main story consumers have already moved to `EquilibriumSolveOptions` and
    `PhaseControlPolicy`, including the bounded phase-control story paths. The
    request objects now also expose an explicit reset helper for the safe fixed
    path, while the enum-based mode setter is transition-only.
- [x] Add a narrow public re-export/prelude for the supported production path:
  phase-system specification/resolution, conditions, initial composition,
  solve options, `solve_resolved_pt`, immutable solution, and typed reports.
  Internal formulation and compatibility modules must not be required imports
  for normal users. A typed `ChemEquilibrium::prelude` now re-exports the
  production path while leaving legacy modules opt-in.
- [x] Connect the existing `SubstanceSystemFactory` /
  `ThermoRepository` resolution workflow to a high-level equilibrium builder
  so callers can request `PhaseSpec -> resolve -> solve` without manually
  constructing per-phase `SubsData` maps. Keep resolution separately callable
  for users who need to inspect or reuse the immutable resolved system. A
  typed `PhaseEquilibriumPipelineRequest` now resolves the spec and can return
  both the immutable resolved system and the accepted solution as one
  transaction.
- [x] Preserve lookup provenance and NIST fallback decisions from the
  repository through build, solve, and final result reports. The resolved
  report now survives into build and solution bundles, and the public pipeline
  outcome forwards the original lookup report directly.
- [x] Decide and test the safe multiphase default explicitly:
  `FixedDeclaredPhases` and `BoundedPhaseControl` must never be selected by an
  accidental or undocumented default. The facade now exposes explicit helper
  constructors, and the default is tested to remain the fixed-declared-phase
  path.

#### P2 - Prove the top-level workflow and retire duplicate orchestration

- [x] Defer a live offline resolved-system story for cycle detection until a
  physically credible fixed-P,T fixture exists. The
  real disappearance, hysteresis, budget-termination, and failed-transition
  rollback stories are now covered; a physical fixed-P,T cycle fixture is
  intentionally not manufactured for a checklist.
- [x] Assert that the live integration layer leaves all thermochemical library
  files byte-for-byte unchanged.
- [x] Migrate the retained user-facing examples and production entry point to
  `solve_resolved_pt`; old/new parity remains confined to characterization
  tests and compatibility modules.
- [x] Deprecate the broad mutable equilibrium facade and old phase-control
  helpers after the active-set matrix passed. The explicit legacy nonlinear
  fallback backends remain supported.
- [x] Update the architecture documents and examples so only the typed
  resolved-phase facade is described as the production path.

#### P3 - Optimize only from measured evidence

- [x] Provide final canonical characterization at the supported scales and
  record projection, Jacobian, nonlinear, outer-loop, and total timings.
  Synthetic 20/50/100/200 active-set timings are recorded in
  `PHASE_CONTROL_BENCHMARK.md`; real exact-element NASA coverage exercises
  20/50/100 candidates from the local NASA gas catalog, using the allowed
  element alphabet `{C,H,O}`. The strict all-three-element subset contains
  only about 20 records and is covered separately.
  - [x] Add an ignored release story over 20/50/100 real C/H/O-limited NASA
    candidates using the same legacy-NR contract and a byte-for-byte library
    immutability check. Release evidence was captured on 2026-07-28 and is
    recorded in `PHASE_CONTROL_BENCHMARK.md`.
    The first run exposed a fixture-contract error rather than a solver limit:
    strict all-three-element matching has only about 20 local NASA gas records,
    while the allowed-element (`SubsetOf`) search supplies 100+ candidates.
    The live regression now pins both semantics and the low-level search keeps
    element sets isolated per `(library, substance)`.
- [x] Cache `ActiveSetProjection` and reduced numeric formulations for repeated
  `A -> B -> A` phase sets. The bounded range runner now retains one numeric
  and, for RST policies, one symbolic preparation per active mask; cache
  cardinality is published in the typed report and covered by live stories.

Explicitly deferred beyond this production engine milestone: GUI, multiple
ideal-gas phases, ideal/non-ideal condensed solutions, and speculative
multi-start. The typed fixed-`P,T` and temperature-range facades are in scope;
physical extensions of the supported phase models remain separate work.

### P4.0 Freeze the fixed-P,T multiphase contract

- [x] Document the thermodynamic ensemble as closed-system Gibbs minimization
  at fixed `P`, `T`, and conserved elemental totals.
- [x] Define the supported-model matrix in code and documentation:
  `IdealGas` maps to the ideal-gas activity law and `PureCondensed` maps to
  unit activity.
- [x] Require a `PureCondensed` phase to contain exactly one component. A
  multi-component condensed phase is a solution and must return
  `UnsupportedModel`, not silently reuse the ideal-solution equation.
- [x] Reject more than one ideal-gas phase until a physically meaningful
  immiscible/multiple-gas model exists.
- [x] Keep physical state separate from activity model. `Liquid`, `Solid`, and
  `Condensed` select records and describe output; `PhaseModel` selects the
  chemical-potential equation.
- [x] State explicitly that zero initial amount means "candidate phase absent
  initially", not "phase excluded from equilibrium". Exclusion must be a
  separate typed input choice.
  `equilibrium_multiphase_domain` now makes this distinction at the physical
  input boundary; a later `PhaseSet` bridge owns explicit exclusion.

### P4.1 Make phase-qualified identity canonical in the equilibrium problem

- [x] Introduce an equilibrium component descriptor that retains
  `PhaseComponentId`, bare substance name, physical state, phase model, and a
  stable display label. The bare name may repeat across phases; the qualified
  component id may not.
  `equilibrium_component::EquilibriumComponentDescriptor` now derives the
  bare name and stable label from the qualified id and carries the exact
  physical-state, phase-model, and solver activity-model mapping.
- [x] Replace `EquilibriumProblem::species: Vec<String>` as the canonical
  identity with the ordered component descriptors. Keep labels as a derived
  view only.
  `EquilibriumProblem` now owns ordered `EquilibriumComponentDescriptor`
  values; `species()` is a derived label view retained only for numerical
  compatibility while the dense legacy formulation is migrated.
- [x] Replace or rename the integer-only equilibrium `PhaseId` so it cannot be
  confused with the phase subsystem's semantic `PhaseId`. Use a typed phase
  index internally and preserve the semantic id at the boundary.
  The dense solver id is now `equilibrium_ids::PhaseIndex`; bridge descriptors
  retain `phase_layout::PhaseId` as the semantic identity.
- [x] Store an ordered equilibrium phase descriptor instead of the current
  bare `Phase { kind, species: Vec<usize> }`; derive index ranges/maps once
  from `SystemLayout`.
  - [x] `EquilibriumPhaseDescriptor` now lives beside the canonical component
    descriptor and is the source of truth in both bridge metadata and
    `EquilibriumProblem`. The legacy numeric `Phase` vector is derived once
    from descriptor ranges for existing residual/Jacobian backends.
  - [x] Retain the derived legacy numeric `Phase` projection only inside the
    numerical compatibility payload. It is not canonical identity and is not
    exposed by the production facade; removing it is tied to a future backend
    payload rewrite and is not useful release work by itself.
- [x] Permit the same chemical substance in multiple phases while continuing
  to reject duplicate `PhaseComponentId` values.
- [x] Carry a layout fingerprint/revision in the bridge result so a solution
  cannot be applied to a different resolved component order.
  `MultiphaseInitialComposition` is now fingerprint-bound to its layout; the
  future solver bridge/result uses the same fingerprint.

Required tests:

- [x] `gas::H2O` and `liquid::H2O` become two solver unknowns with one shared
  molecular formula and two distinct thermochemical records.
- [x] Permuting input map insertion order does not change component, phase,
  matrix, residual, or result order.
  The bridge regression builds the same local NASA gas/condensed system from
  opposite `HashMap` insertion orders and compares metadata, descriptors,
  labels, element matrix, initial coordinates, and `G0(T)` ordering.
- [x] Duplicate qualified ids, duplicate phase ids, and mismatched phase ranges
  fail before thermochemistry or a nonlinear backend is invoked.
  `EquilibriumProblem::new_with_phase_descriptors` now rejects duplicate
  semantic phase ids and non-contiguous/overlapping descriptor ranges at its
  typed input boundary; tests assert that these failures occur before solver
  setup.

### P4.2 Add a typed initial-composition boundary

- [x] Introduce `MultiphaseInitialComposition` aligned to `SystemLayout`, with
  finite non-negative physical mole numbers and at least one positive amount.
- [x] Provide a strict sparse constructor keyed by `PhaseComponentId`. Unknown
  keys are errors; omitted known components become physical zero only through
  an explicit sparse-input policy.
- [x] Reject bare-name maps when one substance occurs in more than one phase.
  Never guess whether `H2O` means gas or liquid.
  The typed boundary accepts only dense layout-order input or sparse
  `PhaseComponentId` input; there is deliberately no bare-name constructor.
- [x] Compute conserved element totals from physical initial amounts before
  trace floors are introduced. Numerical log-mole seeds must not create mass.
- [x] Convert physical zeroes to log-coordinate seeds through the existing
  typed trace policy only after the physical inventory has been validated.
  The bridge computes its element totals from physical input before creating
  `LogMolesInitialGuess`; the local NASA-gas solve verifies that a trace-seeded
  zero `H2O` leaves the `H=4`, `O=2` inventory exactly unchanged.
- [x] Make request construction transactional: failed composition validation
  must not publish a partially prepared problem, cache, or lookup report.
  `PhaseEquilibriumBuildRequest::new` validates the resolved layout,
  composition fingerprint, supported-model policy, and immutable metadata
  before returning the request; it mutates neither subsystem.

Required tests:

- [x] Dense and sparse constructors produce the same ordered vector.
- [x] Ambiguous bare names, negative/NaN/Inf amounts, unknown components, and
  an all-zero inventory return typed errors.
- [x] Trace seeding leaves computed elemental totals bit-for-bit unchanged.

### P4.3 Build one `ResolvedPhaseSystem -> EquilibriumProblem` adapter

- [x] Add the adapter inside `ChemEquilibrium` (for example
  `phase_equilibrium_problem.rs`). The phase subsystem supplies data and typed
  layout; it must not know solver policies, residuals, or backend types.
  - [x] The foundational `phase_equilibrium_problem` module now owns the typed
    request and immutable structural metadata.
  - [x] Add standard-state extraction and final `EquilibriumProblem`
    construction without leaking solver policy into `phase_*`.
- [x] Define a single request containing `&ResolvedPhaseSystem`, fixed
  `EquilibriumConditions`, typed initial composition, trace-seed policy, and
  supported-model policy.
- [x] Return a typed bridge bundle containing the canonical
  `EquilibriumProblem`, `SystemLayout`, component/phase index maps, immutable
  lookup provenance, and a build report. Do not return parallel unrelated
  vectors that callers must keep aligned manually.
  - [x] Introduce immutable `PhaseEquilibriumMetadata` with the canonical
    layout, fingerprint, component/phase descriptors, index maps, and retained
    `ResolvedPhaseSystemReport`.
  - [x] Add the numerical `EquilibriumProblem` and standard-state build report
    only after all thermochemical extraction has succeeded.
- [x] Build the element-composition matrix and element labels in exact
  `SystemLayout` component order. Validate that the same substance resolved in
  two phases has identical molecular composition even when its thermochemical
  records differ.
- [x] Extract one standard-state `G0(T)` function/expression per qualified
  component from its own phase-local `SubsData` record. Do not pass the
  composition-corrected Gibbs value from `evaluate_gibbs`: the canonical
  equilibrium residual already applies mixing and pressure activity terms.
- [x] Evaluate every standard-state model once at the requested temperature
  during preparation and reject missing, out-of-range, NaN, or infinite data
  with the offending `PhaseComponentId` and source provenance.
- [x] Map `IdealGas` and one-component `PureCondensed` to canonical activity
  models; return a typed unsupported-model error for every other combination.
- [x] Preserve `ResolvedPhaseSystemReport` provenance by component, including
  library, record key, physical-state match, and NIST fallback evidence.
- [x] Add a stable preview/report surface showing conditions, ordered
  components, phase models, initial amounts, element totals, thermochemical
  sources, and solver-facing labels without running a solver.

Required tests:

- [x] Adapter output agrees exactly with direct `SubsData` standard-state Gibbs
  and elemental-composition calculations for each component.
- [x] A regression proves that ideal-gas mixing/pressure terms are applied once,
  not once in `SubsData` and again in the equilibrium residual.
- [x] Missing data in the last phase leaves no partially built bridge result.
- [x] Local NASA gas plus NASA condensed records build fully offline and keep
  per-component provenance.

### P4.4 Solve a fixed active phase set through the canonical engine

- [x] Make the pure numeric and symbolic formulations consume the same typed
  phase/activity descriptors; remove `PhaseKind` switches duplicated across
  residual, Jacobian, symbolic generation, validation, and K_eq support.
  `Phase` now stores the canonical `PhaseActivityModel` directly. The numeric
  residual, RST symbolic generator, phase-stability logic, and K_eq guard read
  that same value. `PhaseKind` is only a temporary type alias for older test
  fixture literals, not a second runtime representation.
- [x] Implement and test the chemical-potential/activity term for each
  supported model:
  - ideal gas: `ln(x_i * P / P0)`;
  - one-component pure condensed phase: zero activity correction.
  - `equilibrium_activity` now owns the common `ln(a_i)` contract and the
    pressure offset used by both numeric and symbolic residual construction.
- [x] Build an immutable active-set projection that maps global
  `PhaseComponentId` values to local nonlinear coordinates and can scatter an
  accepted local solution back to the full `SystemLayout`.
  - The ChemEquilibrium-local `ActiveSetProjection` currently uses typed
    `SpeciesId`/`PhaseId`; upgrading its global identity to `PhaseComponentId`
    belongs to the later phase-subsystem bridge.
- [x] Solve the all-declared-phases-active case first through the existing
  backend cascade and backend-independent acceptance gate.
  `PhaseEquilibriumProblemBundle::solve_with` now consumes the prepared
  bridge bundle, configures only `EquilibriumSolverSettings`, and returns one
  immutable solution bundle with the accepted snapshot, build provenance, and
  backend trace. The offline NASA-gas contract covers this end-to-end path.
- [x] Extend the RST symbolic thermochemistry snapshot to build from the typed
  bridge bundle rather than from the legacy single `solver.subs_data` field.
  The bridge now owns phase-qualified symbolic `G0(T)` expressions alongside
  numeric closures. It injects them into the canonical solver, and the RST
  adapter prefers that immutable ordered snapshot before falling back to the
  historical mutable `SubsData` facade. The local NASA-gas bridge regression
  asserts that the accepted backend is RustedSciThe.
- [x] Ensure legacy nonlinear backends can consume the same prepared active-set
  problem as an explicit fallback, without creating a second multiphase
  formulation.
- [x] Keep fixed-active-set solve publication transactional: only an accepted,
  globally re-expanded solution and its complete reports become visible.
  The non-phase-controlled bridge solve is also consuming and transactional:
  invalid backend settings return no accepted result and leave the resolved
  phase data untouched.

Required tests:

- [x] Numeric residual, analytic Jacobian, and RST symbolic residual agree for
  ideal-gas plus pure-condensed fixtures.
  - The mixed fixture runs at `P != P0` and also caught a corrected Jacobian
    defect: phase-total derivatives apply to non-reacting species in the same
    phase as a reaction participant.
- [x] Every configured RST backend sees identical equations and component
  order; fallback attempts preserve the same active-set projection.
  The bridge matrix configures every RST backend independently against one
  local NASA-gas bundle. Each reaches the same bridge-owned symbolic contract;
  methods that do not accept the fixture report a normal one-attempt backend
  failure rather than a setup error. A Nielsen -> LM cascade preserves the
  exact qualified component order and accepts the second attempt.
- [x] A failed backend cascade cannot overwrite a previous accepted
  multiphase result.
  Bridge solves publish only immutable `PhaseEquilibriumSolutionBundle` values.
  A regression retains an accepted NASA-gas snapshot, forces a later Nielsen
  failure, and verifies both accepted moles and retained source provenance are
  unchanged.

### P4.5 Replace legacy phase control with a bounded active-set algorithm

This subsection is about the current helper-level API surface in
`equilibrium_workflows.rs`: `solve_with_phase_control`,
`compute_phase_creation_dg`, `activate_phase`,
`deactivate_phases`, and `deactivate_phases_mass_conserving`.
The complaint is not the phase subsystem interface itself; it is the behavior
and contract of these old helpers.

- [x] Preserve the updated log-mole seed before every outer-loop `continue`
  so the next nonlinear attempt starts from the new phase-adjusted state.
- [x] Add a hard outer-loop iteration cap and a typed
  `PhaseControlDidNotConverge { iterations }` error for phase-control
  non-convergence.
- [x] Retire the current `deactivate_phases_mass_conserving` behavior from
  production paths. Phase removal should seed the deactivated phase to trace
  floor only, not transfer its moles to an arbitrary receiver species.
- [x] Replace the current phase-creation criterion with an explicit phase
  stability contract. Production should support only the physically justified
  one-component condensed-phase case first; multicomponent ideal-solution
  stability should return `ValidationNotApplicable` until tangent-plane
  minimization exists.
- [x] Make `activate_phase` a seed-only helper: mark the phase active, seed it
  with a small positive mole number, and let the next full nonlinear solve
  restore the element balance.
- [x] Keep the active-phase set as persistent state across restarts instead of
  recreating it opportunistically inside one solve call.
- [x] Encode hysteresis explicitly with separate create/keep thresholds and a
  typed transition plan so phase flip-flopping cannot consume the outer loop.
- [x] Do not promote the current `PhaseManager`, `compute_phase_creation_dg`,
  or in-place `activate_phase`/`deactivate_phases` helpers to the canonical
  path. Characterize any useful behavior, then replace their unbounded mutable
  loop with typed orchestration.
  - The raw creation-score and compatibility activation/deactivation helpers
    have been removed. The canonical boundary now exposes only physical
    stability reports and typed seed-only transitions.
- [x] Introduce `PhaseSet`/`PhaseStatus` for declared, active, inactive,
  excluded, appeared, and disappeared phases. Every transition must retain the
  semantic `PhaseId`.
  - `InitialPhaseSet` supports inventory-derived, all-candidate, and explicit
    active/excluded policies. Excluded phases do not enter stability
    classification, while the nonlinear backend sees only the immutable active
    mask for one fixed-set attempt.
- [x] Solve only active components, then compute a typed phase-stability report
  for every inactive candidate. For the first supported slice, implement the
  mathematically justified stability criterion for one-component pure
  condensed phases against the accepted gas/condensed state.
  - Each outer iteration now builds a reduced species/phase projection,
    recomputes its reaction basis, solves it against the unchanged original
    element totals, and expands only the accepted candidate back into declared
    ordering.
- [x] Derive any elemental potentials/reduced chemical potentials required by
  the stability test from the accepted state with a checked linear solve and a
  reported residual; do not infer phase stability from raw `sum(G0)`.
- [x] Remove an active phase only when its amount is below the destruction
  threshold and the resulting inactive phase satisfies the stability
  inequality. Never transfer its moles to an arbitrary first species.
- [x] Add one most-unstable inactive phase per restart, apply only a numerical
  trace seed, and recompute the reduced active problem against the original
  physical element totals. The seed is not claimed to conserve elements; every
  solved candidate must pass the common element-balance gate before it can
  trigger another transition or be published.
- [x] Use separate create/keep thresholds (hysteresis), a maximum restart
  count, and visited-phase-set cycle detection. Report oscillation and budget
  exhaustion as typed non-convergence errors.
- [x] Stage the whole outer solve locally. If any restart, stability test, or
  backend cascade fails, retain the previously published solution unchanged.
  - `solve_candidate_from_seed` now returns a validated unpublished bundle;
    `solve_with_phase_control` publishes it only after the transition
    classifier reaches a fixed point.
- [x] Accept the final result only when the nonlinear candidate gate passes,
  every active phase is internally valid, and every inactive candidate meets
  the phase-stability tolerance.
  - The outer loop refuses to publish material in an inactive phase, rejects
    unsupported multicomponent solution models, rejects missing/rank-deficient
    elemental-potential references, and checks the elemental-potential fit
    residual before using a phase driving force.

Required tests:

- [x] A stable absent condensed phase stays absent.
- [x] An unstable absent condensed phase appears and converges after a bounded
  restart.
- [x] A vanishing unstable/stable boundary fixture exercises hysteresis without
  cycling.
- [x] Deliberately alternating phase sets terminate with a typed cycle/budget
  report instead of looping forever.
- [x] A failure after an accepted fixed-set candidate but before outer-loop
  convergence leaves the previously published solution and reports unchanged.
- [x] Every intermediate and final state conserves each element within the
  common acceptance tolerance.
  - Transition records retain the candidate validation evidence that admitted
    each fixed-set solve; trace seeds are restart coordinates, not published
    thermodynamic states.

### P4.5a Production hardening and regression evidence

The bounded active-set architecture is now canonical, but production readiness
requires proof that every reduced solve preserves the original closed-system
inventory and that near-boundary numerical behavior cannot manufacture phase
transitions.

- [x] Assert in `solve_fixed_active_set_candidate` and its regression tests
  that element totals always originate from the full physical problem, never
  from a reduced seed or a trace-seeded inactive component vector.
  The local solver now receives totals recomputed from full `n0` and the full
  element matrix; the regression compares the accepted phase-control moles
  with that original inventory.
- [x] Validate active-set feasibility before invoking a nonlinear backend:
  the original element-total vector must belong to the column space of the
  active species element matrix. Reject an impossible active set with typed
  `InvalidProblem { field: "phase_active_set", .. }`.
  `ActiveSetProjection` now checks the reduced SVD feasibility residual before
  any backend call. Its pre-existing square/rank gate is now classified as the
  same typed invalid active-set error instead of `ValidationNotApplicable`.
- [x] Document and test the strict reduced-system contract: the active element
  matrix must retain every conserved-element direction required by the
  physical inventory, and the reduced reaction/elements equation count must
  be square in the present log-moles formulation. Explain why a phase set
  without a carbon-bearing component cannot solve a carbon-containing system.
  The projection documentation and its impossible H/C active-set regression
  now make this boundary explicit.
- [x] Define the inactive ideal-gas contract. Until a true gas-phase split or
  tangent-plane criterion exists, trace-floor gas species must not generate
  an artificial stability driving force or activate a second gas phase.
  `compute_phase_stability_reports` classifies every gas phase as
  `FixedIdealGas` with no driving force; a regression uses an extreme Gibbs
  value on an inactive trace gas and proves it cannot enter the transition
  plan.
- [x] Complete the activity-law audit with a source-level regression: numeric
  residual, analytic Jacobian, symbolic residual, and phase-stability code
  must delegate to `equilibrium_activity`; no duplicate `ln(a_i)` formula may
  be added outside that module.
  The deprecated residual path now also routes phase offsets through the same
  helper, and a regression test locks the canonical and deprecated residual
  behavior together.
 - [x] Expand `ActiveSetProjection` tests for non-consecutive global indices,
   multiple active phases, local phase remapping, scatter/project round trips,
   preserved ordering, inactive trace floors, reduced basis dimensions,
   dependent element columns, and impossible active-set rejection.
  - [x] Sparse non-consecutive active species now round-trip through
    `project_log_moles` and `scatter_log_moles` while preserving global
    ordering and inactive floors.
  - [x] Projected local values now round-trip through `scatter_log_moles` and
    `project_log_moles` without reordering.
  - [x] The reduced basis reports the expected reaction dimension for a
    nontrivial active set, so the projection contract is pinned instead of
    inferred.
  - [x] Impossible sparse active sets now fail before backend publication
    instead of being silently coerced into a reducible phase inventory.
- [x] Expand phase-control stories: appearance, disappearance, appearance
  followed by disappearance, disappearance followed by reappearance, several
  sequential transitions, maximum-restart termination, A -> B -> A cycle,
  and transactional preservation after a failed intermediate restart.
  - [x] A first transactional story now covers appearance followed by a
    failed restaging attempt, and proves that the previously published
    solution and moles remain intact.
  - [x] The intermediate publication also survives a second solve attempt
    that fails before changing the published active set or mole state.
  - [x] Maximum-restart termination is covered by a typed story that confirms
    the solver stops cleanly after one transition budget is exhausted.
  - [x] The activation/deactivation seed helpers now cover a full
    appearance -> disappearance -> reappearance cycle without relying on a
    brittle full-solve stress case.
  - [x] Several sequential transition scenarios are now covered through
    fresh fixed-fixture solves so the publication contract remains stable
    across appearance, disappearance, and reappearance cases.
- [x] Add deterministic numerical stress fixtures covering trace species,
  large Gibbs-energy scales, nearly rank-deficient reaction bases, several
  pure condensed candidates, and low/high temperatures. Assert finite values,
  accepted element conservation, and typed failure rather than NaN/Inf.
  - [x] A first stress fixture now exercises trace species against extreme
    Gibbs-energy scales across low, nominal, and high temperatures, while
    checking finite accepted moles and exact element conservation.
  - [x] A second stress fixture now covers multiple pure condensed
    candidates and confirms that the phase-control solve remains finite and
    element-conserving.
- [x] Run at least LM and Trust-Region policies through the same phase-control
  fixtures and compare accepted phase sets, mole vectors, and validation
  reports within documented tolerances.
  - [x] LM and TR now match on both a pure-condensed appearance fixture and
    a stable gas-only fixture, with checked phase sets, mole vectors, and
    element-balance tolerances.
- [x] Make hysteresis temperature-aware. Store dimensionless create/keep
  coefficients or derive thresholds from the solve temperature immediately
  before phase control; remove the current accidental 298 K default from the
  physical policy contract. Phase control now resolves the hysteresis band
  from the current solve temperature instead of baking in a fixed 298 K
  constant, and the policy itself can be expressed either explicitly or as
  dimensionless RT factors.
- [x] Add ignored, opt-in performance characterizations for 50--100 species
  and several phases. Record outer iterations, nonlinear iterations,
  transitions, and elapsed time, but do not impose flaky wall-clock limits in
  the default test suite.
  - [x] A first ignored characterization scaffold now runs a larger
    multi-phase synthetic system and prints iterations, transitions, and
    elapsed time without affecting the default suite.
- [x] Establish an always-offline long-term regression matrix for O2 <-> 2O,
  N2 <-> 2N, N/O mixtures, condensed appearance/disappearance, trace species,
  and deterministic generated small systems. Each fixture must state a
  physical invariant, not merely retain a historical iterate.
  - [x] A first offline matrix now covers O2/O dissociation, N2/N
    dissociation, a diluted N/O gas mixture, and synthetic condensed
    appearance/disappearance, with finite accepted states and bounded
    element-balance drift.
  - [x] A second offline matrix now covers additional deterministic small gas
    systems, including NO/N2/O2 and diluted O2/O/N2 cases, without network
    access.

### P4.6 Publish a phase-aware immutable solution and reports

- [x] Add `MultiphaseEquilibriumSolution` containing fixed conditions, layout
  revision/fingerprint, ordered component moles, phase totals, phase-local mole
  fractions, active/inactive status, and the accepted canonical solution.
  - [x] The fixed-active bridge result is now published as an immutable
    `MultiphaseEquilibriumSolution` with fingerprint checks, phase totals,
    local mole fractions, qualified lookups, provenance, and backend evidence.
  - [x] Bridge-backed bounded phase control now publishes its accepted
    `PhaseControlledSolveReport`, final `PhaseSet` statuses, and complementarity
    acceptance report through the same immutable solution type. The numerical
    inner loop is now owned by `PreparedPhaseControlRunner`, which constructs
    immutable reduced problems from the bridge-owned data and symbolic `G0`
    snapshot; `EquilibriumLogMoles` is retained only for compatibility paths.
- [x] Provide lookups by `PhaseComponentId` and `PhaseId`; expose aggregate
  totals by bare substance only as an explicit derived view.
- [x] Retain build/lookup provenance and the complete nested backend solve
  report in the result bundle.
  - [x] When fixed-phase K_eq validation is enabled, its typed status is also
    retained and emitted as a stable result-summary row instead of being lost
    inside the mutable solver host.
- [x] Add `PhaseTransitionReport` entries with previous/new phase sets, reason,
  seed, stability metric, backend outcome, and elemental-balance evidence.
  - The immutable `PhaseControlledSolveReport` now retains typed initial/final
    phase sets, ordered transitions, explicit physical reasons, exact restart
    seeds, phase totals, driving forces, per-candidate validation evidence, the
    final validation report, and every nonlinear backend report. The remaining
    report work is integration into the final multiphase solution bundle.
- [x] Add a final `MultiphaseAcceptanceReport` combining canonical residual and
  element checks with phase stability/complementarity checks.
  - [x] The final bundle now exists and combines phase-control evidence,
    canonical validation, phase-stability reports, and a complementarity
    summary with stable rows and `Display` output.
- [x] Provide stable summary rows and `Display` output for CLI, future GUI, and
  snapshot tests; core code must not print directly.
  - [x] `PhaseControlledSolveReport` now exposes stable summary rows and a
    `Display` implementation for CLI and snapshot-friendly output.
- [x] Reject result reconstruction when the result layout does not match the
  resolved-system fingerprint.
  - [x] `MultiphaseInitialComposition` already refuses reconstruction against
    a foreign layout fingerprint, and the regression tests now pin that
    contract down explicitly.
  - [x] `MultiphaseEquilibriumSolution` now independently rejects metadata and
    build provenance originating from distinct valid resolved layouts before
    publishing a queryable result.

### P4.7 Complete the multiphase test matrix

- [x] Add a dedicated `equilibrium_phase_bridge_tests.rs` for identity,
  ordering, adapter validation, provenance, and transactional failures.
- [x] Add a dedicated `equilibrium_multiphase_story_tests.rs` with module-level
  documentation describing each physical hypothesis and expected result.
  - [x] Cover the remaining live transition scenarios supported by the current
    fixed-P,T model contract:
  - [x] ideal-gas-only parity with the current canonical solver;
  - [x] gas plus one stable pure solid;
  - [x] gas plus one stable pure liquid;
  - [x] the same molecule represented in gas and condensed phases;
  - [x] two independent pure condensed candidate phases;
  - [x] unsupported multi-component condensed solution;
  - [x] phase appearance on live resolved thermochemical data;
  - [x] disappearance within one continued solve on live resolved
    thermochemical data, including positive-inventory liquid boundary recovery;
  - [x] hysteresis retention on live resolved thermochemical data;
  - [x] cycle detection remains covered by deterministic synthetic state-machine
    regressions; a live fixture is explicitly deferred until credible physics
    produces one;
  - [x] offline mixed NASA gas/NASA condensed lookup provenance.
  - [x] The dedicated equilibrium_phase_bridge_tests.rs module now covers real offline NASA gas/condensed provenance, a real NASA gas plus solid H2O(s) fixture, the same molecule as two phase-qualified components, bounded inactive-liquid startup, and transactional no-mutation-on-failure publication.
  - [x] The new `equilibrium_live_data_tests.rs` module adds a separate live-data layer over the bundled local thermochemistry repository: explicit offline lookup policy, real NASA gas resolution, one stable real gas solve, one bounded real multiphase solve, and transactional non-mutation of the resolved source view.
  - [x] The live-data layer now also includes a real water gas/liquid
    temperature-shift regression: 350 K carries a non-trace liquid inventory,
    550 K shifts the inventory toward vapor, both solves stay inside the
    `H2O(L)` record range, and layout identity remains unchanged.
  - [x] The live bounded pipeline now snapshots the canonical substance-base,
    address-catalog, and element-composition JSON files before solve and proves
    that all three remain byte-for-byte unchanged afterward.
  - [x] A live `NASA_gas` H2/O2/H2O reaction fixture now solves through the
    public pipeline with element conservation and an accepted independent
    `K_eq` cross-validation report. `EquilibriumSolveOptions` exposes this as
    an explicit typed opt-in rather than requiring callers to mutate backend
    settings directly.
  - [x] The same reactive fixture pins stable fingerprints of its three source
    NASA7 records. Updating local thermochemical data now makes the fixture
    fail loudly until its physical contract is deliberately reviewed.
  - [x] A live explicit-state lookup now resolves `gas::H2O` and
    `solid::H2O(s)` as distinct phase-qualified records from the local
    NASA_gas/NASA_cond catalog, with NIST fallback disabled.
  - [x] Fix phase-qualified live lookup and thermochemistry preparation.
    `PhaseSpec::physical_state` is now installed as a typed `SubsData` lookup
    constraint, and bridge construction selects the NASA coefficient interval
    at the requested temperature before snapshotting numeric/symbolic `G0(T)`.
    Previously liquid water could resolve as `NASA_gas:H2O`, while freshly
    parsed NASA calculators silently published zero-coefficient `G0(T)=0`.
  - [x] Replace the invalid ice diagnostic built on zero `G0`. The real
    `H2O(g)/H2O(s)` fixture at 250 K resolves
    `NASA_gas:H2O` plus `NASA_cond:H2O(s)`, activates the initially absent ice
    phase, and transfers almost all water inventory into the solid while
    preserving the gas oxygen inventory.
  - [x] Defer Fe phase-lifecycle fixtures until defining an exact-record
    selection policy for condensed polymorphs (`Fe(a)`, `Fe(c)`, `Fe(d)`).
    The earlier zero-driving-force observation is invalid because it was made
    before the zero-`G0` bridge defect was fixed; do not retain it as physical
    evidence or let a state-only resolver guess a polymorph.
  - [x] Promote a physically scaled Boudouard fixture
    `2 CO(g) <=> CO2(g) + C(gr)` into the live matrix. A finite CO/CO2 initial
    mixture avoids trace-floor chemical potentials, activates real
    `NASA_cond:C(gr)` at 700 K, and shows lower graphite inventory at 1400 K.
    The deliberately harsher all-CO boundary case remains a separate
    initialization benchmark rather than the production lifecycle fixture.
  - [x] The local NASA_gas H2O/O2 plus NASA_cond H2O fixture now reaches
    the bounded bridge with the zero-inventory liquid phase initially inactive.
    The outer loop derives its first phase set from physical `n0`, not the
    positive numerical trace coordinate, and solves the square gas-only
    projection without relaxing acceptance tolerances.
  - [x] `AllCandidatePhases` now has the same physical boundary contract:
    inventory-free candidates are normalized to inactive before a positive
    log-moles solve, while remaining eligible for stability-driven activation.
    The live high-temperature water story exercises this path without a
    synthetic trace inventory.
  - [x] The public water phase-pair story now shows temperature-driven vapor /
    condensed dominance on the real `solve_resolved_pt` facade rather than only
    in synthetic workflow tests.
  - [x] Accept that cycle detection during one continued live run is pinned
    primarily by synthetic workflow regressions. Positive-inventory
    disappearance, hysteresis retention, budget termination, and rollback are
    now covered by real local thermochemistry; these must not be conflated
    with zero-inventory normalization.
    Ice, liquid water, and graphite now provide stable offline appearance and
    temperature-shift fixtures for extending that matrix without inventing
    thermochemistry.
  - [x] Keep rank-deficient-but-representable active sets valid. The active-set
    projection derives the independent reaction/element rank and separately
    verifies that the closed-system element totals lie in the range of
    `A_active^T`; it no longer equates the number of declared element columns
    with the number of independent constraints.
- [x] Define an explicit element-candidate search contract for real-data
  equilibrium fixtures. `ThermoData` and `SubsData` now expose
  `ElementSearchMode::{AnyRequested, SubsetOf, ExactSet}` and a dedicated
  `search_by_exact_elements` method. The historical `search_by_elements_only`
  remains the broader subset-of-elements mode; it is not silently redefined.
- [x] Build a production candidate-selection layer on top of element search.
  It must intersect the chosen element mode with ordered thermochemical
  library preference, physical-state/phase policy, temperature-interval
  support, and per-record provenance before constructing an equilibrium
  problem. `EquilibriumCandidateSelector` now performs this as a read-only
  deterministic transaction and returns selected records plus explicit
  rejection reasons. It chooses one acceptable record per substance according
  to library preference, preserves physical-state evidence and record keys,
  screens recognizable coefficient intervals, and leaves unknown interval
  schemas visible as `Unknown` rather than silently rejecting them.
  - [x] Add `EquilibriumCandidatePhasePlan` and the
    `PhaseEquilibriumPipelineRequest::from_candidate_selection` boundary.
    Every selected exact `record_key` must be assigned once to an explicit
    phase/model; the selected library is pinned through explicit lookup
    instructions. The builder rejects omitted, duplicate, or unknown records
    transactionally and never guesses a phase model from element data.
- [x] For every applicable small fixed-phase set, compare the accepted result
  with the independent equilibrium-constant solver. Report non-applicability
  explicitly for phase-appearance decisions that the K_eq problem does not
  model. The gas-only fixed-phase bridge now retains a `Compared(report)` K_eq
  validation status in the immutable summary, and the acceptance evidence is
  checked in the dedicated story test.
  - [x] Keep `WhenApplicable` validation observational: failure of the
    independent K_eq solver is retained as `ValidatorFailed` evidence and does
    not reject an otherwise accepted canonical solution. `Required` remains
    the explicit fail-closed mode.
  - [x] Scale the default total-Gibbs comparison tolerance for live problems.
    The validator still requires close species amounts and fractions, while
    avoiding false rejection from a few microjoules of absolute objective
    drift on megajoule-scale solutions.
- [x] Add exact component-order and numeric/closure/symbolic equivalence
  tests, then run the existing RST backend/fallback matrix over at least one
  physical multiphase fixture. The real NASA gas + condensed-water bridge test
  now compares prepared residual/Jacobian, legacy closure/Jacobian, and
  symbolic residual/Jacobian on one real multiphase fixture, while the
  bounded phase-control matrix still pins the explicit legacy-vs-RST policy
  parity on the same physical system.
  - [x] A physical bounded phase-control fixture now compares explicit legacy
    and explicit RST solver policies on the same resolved multiphase system and
    confirms that the accepted component order and moles remain aligned.
  - [x] Make the implicit production solver policy RST-first while retaining
    LM/NR/TR legacy fallbacks. Explicit single-backend policies remain strict
    and never acquire an undeclared fallback.
- [x] Adopt physically meaningful fixtures from the retired classical stack
  only after restating their expected invariants. Oxygen dissociation,
  reactive H/O, water phase pairs, ice, and graphite now test conservation,
  residual, provenance, and lifecycle contracts; fixtures tied only to legacy
  mutable state or one historical iterate were not retained.
- [x] Keep all default tests offline and leave thermochemical library files
  byte-for-byte unchanged. The dedicated live-data regression snapshots every
  canonical JSON file used by normal thermo lookup around a real bounded solve.

### P4.8 Migrate workflows and retire the duplicate engine

- [x] Add one public one-shot facade such as `solve_resolved_pt` accepting the
  resolved phase system, typed initial composition, conditions, solver policy,
  and phase-control policy.
  - [x] `ResolvedPhaseEquilibriumRequest` and `solve_resolved_pt` now own the
    fixed-declared-phase bridge transaction, including typed numerical
    settings and immutable result publication.
  - [x] `PhaseEquilibriumSolveMode::BoundedPhaseControl(PhaseManager)` now
    routes through the same bridge-owned thermochemistry and publishes one
    `MultiphaseEquilibriumSolution` with phase-control and complementarity
    evidence. The mutable solver remains an internal compatibility host, not
    an alternative public workflow.
- [x] Migrate useful legacy equilibrium workflows and examples one
  physical scenario at a time onto that facade.
  - [x] The local NASA ideal-gas guide now uses `ResolvedPhaseSystem`, typed
    initial composition, and `solve_resolved_pt`; it is the first retained
    user-facing scenario that does not call legacy `gas_solver` directly.
  - [x] The resolved-phase example now uses the repository-backed typed
    pipeline as well; callers no longer need to assemble `SubsData` maps or
    `ResolvedPhaseSystem` manually for the standard guide scenario.
  - [x] A source audit found no non-test consumer outside the compatibility
    module that still calls `gas_solver*`, `solve`, or
    `solve_with_phase_control`; remaining calls are characterization coverage.
    The mutable temperature-range methods remain deprecated now that the typed
    range facade is available.
- [x] During migration, compare both implementations only in characterization
  tests; do not expose two production APIs as equivalent long-term choices.
- [x] Add a typed temperature-range facade above the fixed-`P,T` production
  solve. Resolve the real-data candidate set once, reuse the element matrix,
  reaction basis, and active-set projection while the layout is unchanged,
  refresh temperature-dependent standard-state data for each point, use the
  previous accepted physical solution as the continuation seed, and publish
  each point transactionally. Rebuild only after an accepted phase transition.
  - [x] Make the range request own an immutable resolved layout and a typed
    point grid; do not route the canonical range path through the mutable
    `EquilibriumLogMoles` sweep. `TemperatureRangeRequest` and
    `PhaseEquilibriumPipelineRequest::solve_temperature_range` now provide
    this fixed-declared-phase path.
  - [x] Add per-point timing records for repository lookup, coefficient
    refresh, numeric closure refresh/construction, symbolic refresh or reuse,
    equation/Jacobian preparation, nonlinear solve, phase control, validation,
    and postprocessing. Keep a sweep summary with total, mean, median, and
    worst-point durations. The fixed-declared range now publishes the full
    typed point timing report produced by the solve and a
    `TemperatureRangeDurationSummary`; bounded points additionally publish
    transition counts, phase-set reuse, and projection-cache cardinality.
  - [x] Distinguish one-time setup from per-point work in the report. The
    initial baseline of approximately 445.9 ms (NR), 457.9 ms (LM), and
    466.2 ms (TR) for five points was diagnosed as an accidental RST
    promotion: the compatibility wrapper left `solver_policy` unset even
    after receiving a legacy `Solvers` selector. After the boundary fix,
    the same debug run measured approximately 5.6 ms (NR), 9.2 ms (TR), and
    16.2 ms (LM) for five points. The typed report now separates
    `initial_formulation_timing` from accepted-point timing and exposes
    total/mean/median/worst point durations. The release baseline is recorded
    below for the retained real-data characterization story.
    - [x] Release characterization is now recorded for the real 20-species,
      three-point story: legacy NR took 494 us ascending / 409.4 us descending
      across accepted points; default RST took 106.5217 ms ascending. The RST
      run reused symbolic preparation for 2 of 3 points; the first point
      crosses a real NASA coefficient interval boundary. These are measured
      point totals, not a claim that the two backend families have equivalent
      setup costs.
  - [x] Add a regression comparing one fixed-`P,T` solve with the equivalent
    one-point range solve. The real local fixture compares accepted moles under
    the same legacy backend and keeps any setup cost visible in the timing
    report.
  - [x] Reuse the element matrix, phase/component layout, reaction basis, and
    allocation buffers when the active phase set is unchanged. Add counters
    for reused versus rebuilt objects so a faster result cannot be caused by
    accidentally dropping required work. The fixed range reports one
    formulation build and point-level formulation reuse counters.
  - [x] Refresh only temperature-dependent standard-state data at each point;
    measure numeric closure rebuild separately from symbolic expression reuse.
    Do not assume symbolic expressions are reusable if their captured values
    are temperature-specific; the current reusable RST path updates the
    shared `T` parameter and records that decision explicitly. Coefficient
    interval changes still require the symbolic-refresh branch below.
  - [x] Use the previous accepted physical solution as the continuation seed
    and test both ascending and descending grids. Record seed source and
    rejected/failed points without publishing a partial sweep as complete.
    The ignored real-data release story covers both directions and verifies
    that a failed point cannot publish a partial `TemperatureRangeSolution`.
  - [x] Add a real-data backend matrix over the same temperature grid for RST
    LM/Minpack/Nielsen/trust-region/Powell/Newton and legacy LM/NR/TR. The
    ignored release story records success/failure, per-point residual,
    conservation error, and timing; every run uses `SolverPolicy::Single`, so
    a cascade cannot hide a backend-specific failure. Target-machine release
    measurements remain an operator-run evidence step.
  - [x] Add phase-control range stories with an accepted phase transition and
    verify that projection/layout rebuild counters increase only after that
    transition. A temperature point without a transition must not rebuild the
    full active-set problem. The real offline water/ice fixture now covers the
    transition path, while the bounded live story covers the no-transition
    reuse path and reports both projection and reduced-formulation caches.
  - [x] Extend bounded-range reuse to the RST symbolic backend itself. Each
    active-set cache entry retains one `RstPreparedProblem`; unchanged Gibbs
    expressions update only its shared `T` parameter, while a new active mask
    or symbolic snapshot creates a new entry. The range report exposes the
    retained symbolic-entry count and live stories assert it.
  - [x] Add a scale-aware conservation and residual contract for every bounded
    range point, plus a byte-for-byte no-file-mutation check for the live local
    repository.
  - [x] Add the release characterization in release mode over at least 20,
  50, and 100 real species where exact element search supplies enough local
  records; retain the smaller five-point story as the fast regression
  fixture. The printed timing baseline still has to be collected on the
  target machine.
- [x] Make the retained `gas_solver_for_T_range` compatibility boundary honor
  its explicit legacy `Solvers` argument. It now installs
  `SolverPolicy::legacy_default` instead of allowing symbolic context to
  select the RST production cascade implicitly. The real T-range story pins
  that every accepted point reports a legacy backend.
- [x] Add an opt-in typed timing report to the canonical fixed-`P,T` solve.
  `EquilibriumTimingMode::Enabled` publishes repository lookup,
  thermochemistry, numeric-closure, symbolic, equation, numerical-preparation,
  projection, nonlinear-solve, phase-control, validation, postprocessing, and
  total durations without logging or clock reads on the default disabled path.
   - [x] Extend the timing report to the typed temperature-range facade. Record
  per-point durations plus cache/layout/projection reuse and rebuild counters;
  distinguish coefficient refresh from a genuine active-phase transition so
     range benchmarks expose the cost of changing temperature rather than only
     the cost of one fixed-`P,T` solve. Fixed-declared and bounded-range timing,
     formulation-reuse counters, active-set cache cardinality, symbolic RST
     cache cardinality, and total/mean/median/worst point summaries now exist.
- [x] Add a typed postprocessing adapter from `TemperatureRangeSolution` to the
  immutable interpolation/report layer. It preserves phase-qualified component
  labels and physical component amounts, supports descending solver grids by
  sorting a copied view, and keeps raw-row helpers as a lower-level API.
- [x] Move the direct mutable orchestration modules behind crate-private
  implementation boundaries and expose their retained compatibility surface
  only through `ChemEquilibrium::legacy`; LM/NR/TR numerical fallback code is
  intentionally retained there. A compatibility test now resolves
  `legacy::gas_solver` through that namespace, and the guide no longer presents
  the historical module paths as peer production APIs.
- [x] Migrate production equilibrium consumers away from the broad
  `ThermodynamicsCalculatorTrait`, raw Gibbs closures, and nested legacy phase
  maps to the narrow typed bridge.
  A source audit found no non-test production consumer of those representations;
  they remain only in the explicitly namespaced compatibility subsystem.
- [x] Deprecate `solve_with_phase_control` and the old mutable phase helpers
  after the active-set story matrix passes. The deprecation is advisory: the
  legacy implementation remains available as a compatibility host, while the
  resolved facade never calls it.
- [x] Remove duplicate mutable orchestration from the production public surface.
  The old workflow remains available only under `ChemEquilibrium::legacy` by
  explicit compatibility policy. Keep the handwritten legacy nonlinear
  backends as explicit policy-selected fallback implementations.
- [x] Update the phase-subsystem architecture documents with the final
  dependency direction and add a fixed-P,T multiphase usage example. The
  retained legacy dataflow is labeled historical compatibility rather than
  presented as a second production architecture.

### Recommended P4 implementation passes

1. [x] Component/phase identity refactor in `EquilibriumProblem` plus tests.
2. [x] Typed initial composition and pure `ResolvedPhaseSystem` adapter.
3. [x] Standard-state thermochemistry/provenance bridge and preview report.
4. [x] Fixed-active-set numeric, symbolic, RST, and legacy-fallback parity.
5. [x] Bounded phase-stability active-set orchestration.
6. [x] Phase-aware solution/report boundary and full acceptance gate.
7. [x] Initial offline physical story matrix and independent K_eq
   cross-validation. The remaining live transition matrix is tracked in the
   production P2 queue above.
8. [x] Classical workflow migration and production-boundary cleanup. Retained
   compatibility workflows are isolated under `ChemEquilibrium::legacy`.

### P4 definition of done

- [x] One typed request can solve a resolved closed multiphase system at fixed
  `P,T` without flattening phase-qualified identity.
- [x] Standard-state thermochemistry is evaluated once per component and every
  record retains lookup provenance.
- [x] All supported activity terms are explicit and tested; unsupported phase
  models fail before a backend starts.
- [x] Phase appearance/disappearance is bounded, transactional, conservative,
  and included in the acceptance report.
- [x] The result can be queried unambiguously by phase and component and is
  tied to the exact resolved layout that produced it.
- [x] Numeric, analytic-Jacobian, and symbolic formulations agree; the RST
  backend matrix and explicit legacy fallback pass the same physical fixture.
- [x] Offline gas/condensed stories pass elemental, residual, stability, and
  applicable independent K_eq validation gates.
- [x] No retained production workflow requires the duplicate classical
  or mutable equilibrium orchestration. The handwritten legacy nonlinear
  backends remain supported only as explicit fallbacks.

### Post-P4 optimization and review notes

These items are useful, but they are not blocking the core correctness work
above. Keep them in mind once the phase-control and multiphase acceptance
paths are stable.

- [x] Verify the phase-driving-force contract directly from chemical
  potentials, and add one regression that makes the formula visible at the
  phase-boundary level rather than only through downstream transition plans.
  The new regression checks a pure-condensed candidate against an active gas
  reference and pins the `g0 + RT ln(a)` activity contribution in the
  reported driving force.
- [x] Add a dedicated hysteresis story that exercises
  `appear -> keep -> disappear -> reappear` behavior explicitly and pins the
  transition thresholds against accidental regression.
  The regression now exercises the public phase-classification path with an
  explicit hysteresis band and checks appearance, keep, disappearance, and
  reappearance decisions in order.
- [x] Cache `ActiveSetProjection` and surrounding fixed-layout preparations for
  repeated active sets, including `A -> B -> A`. Bounded temperature sweeps
  expose projection, prepared-formulation, and symbolic cache cardinalities.
- [x] Add a separate benchmark suite for larger synthetic systems
  (20, 50, 100, 200 species) that reports outer iterations, nonlinear
  iterations, projection build time, Jacobian time, and total solve time.
  The ignored benchmark sweep now covers 20, 50, 100, and 200 species
  synthetic inventories and prints projection-build and solve timing together
  with iteration counts.
- [x] Do not rewrite allocation-heavy paths without measured evidence. Current
  characterization identifies nonlinear solve/backend setup as the dominant
  real-data cost; further allocation work is explicitly performance-driven.
- [x] After the phase* bridge was connected, add live resolved-phase-system
  regression fixtures and gradually migrate selected synthetic stories to
  those real inputs. Water gas/liquid/ice, graphite, reactive H/O, lookup
  provenance, conservation, lifecycle, and transactionality are now covered;
  synthetic cases remain as the fast unit layer.
- [x] Add an opt-in real large-system characterization: select 20
  C/H/O-limited candidates from the local element catalog, resolve them from
  `NASA_gas`,
  solve the fixed-`P,T` system, and print the stage timing report. This is a
  performance/diagnostic fixture, not a replacement for the smaller stable
  physical regressions.
- [x] Run the same real 20-species problem through each concrete RST and
  legacy backend as isolated `Single` policies. Keep backend failures visible
  instead of hiding them behind the production cascade, and require every
  accepted backend result to satisfy finite-mole, residual, and scale-aware
  elemental-balance checks.
- [x] Add a real 20-species, five-point temperature-range characterization for
  legacy LM/NR/TR. It verifies ordered accepted points, no failed temperatures,
  finite positive moles, preserves per-point backend reports, rejects hidden
  RST promotion, and records one sweep duration per backend. This is
  deliberately labeled compatibility coverage until the typed range facade
  replaces the mutable sweep.
- [x] Add the corresponding real-data temperature-range characterization after
  the typed range facade is finalized; report per-point coefficient refresh,
  continuation, rebuild, and solve timings. The ignored live story
  `live_large_element_limited_typed_temperature_range_story` covers real local
  NASA data in both directions and checks one formulation build, continuation,
  symbolic reuse, and point reports. The release run now passes and prints the
  measured point timing baseline; a broader 20/50/100-species release matrix
  remains separate characterization work.
- [x] Add an ignored release characterization over 100 real local NASA
  C/H/O-limited species and 50 temperature points, running every supported RST
  and retained legacy backend as an isolated `Single` policy. The test reports
  per-backend total/mean/median/worst timings, formulation builds and reuses,
  symbolic updates, residuals, conservation error, continuation usage, and
  preserves the JSON-library snapshot. Rebuilds are allowed at genuine NASA
  coefficient-interval boundaries; the contract requires complete build/reuse
  accounting rather than assuming one formulation for the whole grid. The
  target-machine release run is now recorded in `PHASE_CONTROL_BENCHMARK.md`:
  six backends succeeded and three reported explicit first-point failures.
- [x] Diagnose backend-specific failures and rare slow points from the 100 x 50
  release matrix before changing the production default. The range error
  boundary now preserves the boxed `ReactionExtentError` instead of flattening
  the backend trace into a string; the ignored 100 x 50 story prints the three
  slowest accepted points for each successful backend with seed-to-result
  log-coordinate movement, detailed stage timings, and RST counters, plus a
  separate typed attempt table for every failure. The strict monolithic P,H
  matrix additionally records failure point and `ReactionExtentErrorKind`.

  A concrete defect has now been corrected in the canonical fixed-`P,T` RST
  path: symbolic methods previously received unbounded log-mole coordinates,
  so a trial step could make `exp(y)` overflow before the common candidate
  validator ran. `PreparedEquilibriumProblem` now derives finite per-species
  box bounds from exact elemental capacities and the caller's trace seed, and
  passes them through `RustedSciTheSolveContract` into RST. The focused live
  20-species C/H/O backend matrix now accepts Nielsen LM, Powell Dogleg, and
  legacy TR; the old Nielsen non-finite residual no longer reproduces there.
  Damped Newton instead exposes an independent stagnation/high-residual
  rejection, which remains a valid method-specific outcome.

  The enhanced target-machine report isolated the common 1010 K anomaly:
  RST residual/Jacobian callbacks and the nonlinear solve took only about
  10--20 ms, while KiThe rebuilt `SymbolicNonlinearProblem` after a NASA
  coefficient-interval change. The canonical fixed-`P,T` graph now carries
  `G0_i` as checked equation parameters. A range point updates `[T, G0...]`
  in place, so an unchanged active set may never rebuild its symbolic graph;
  the same rule now applies to the cached bounded phase-control runner.
  Initial graph construction is recorded in `initial_formulation_timing`, and
  live 20-species stories assert zero per-point formulation build across an
  unchanged layout / active set. The 100 x 50 operator table now prints
  initial setup, initial symbolic construction, and initial numerical-problem
  preparation separately from later per-point formulation builds. The bridge
  builder now also records its outer `total`; a focused local-NASA test guards
  the invariant that setup total encloses every recorded setup stage. A second
  local-NASA contract compares residual and Jacobian entries of the
  parameterized RST graph directly against `PreparedEquilibriumProblem`, so a
  performance change cannot silently alter the fixed-`P,T` equations.

  - [x] Re-run the 100 x 50 release matrix after parameterizing `G0_i`.
    The recorded release table confirms 49 reuses, 49 parameter updates, and
    zero later `formulation_build` for every successful RST backend. The old
    9--11 second coefficient-boundary point is gone; the remaining roughly
    4.2--4.3 second cost is the explicitly reported one-time symbolic setup.
    A test-only baked-vs-parameterized 20-species Damped-Newton comparison
    rejects both graphs under the same strict gate, while an entrywise
    residual/Jacobian contract proves the reusable graph matches the canonical
    real NASA formulation. Rebuilding the historical baked 100-species graph
    is itself too expensive for a practical regression (it exceeded two
    release minutes), which is precisely why it is not retained as a fallback.
  - [ ] Classify the remaining strict-single-backend outcomes without weakening
    acceptance: the latest 100 x 50 release matrix shows Nielsen LM stopping
    at the first point with a finite residual after `MaxIterations`, Damped
    Newton stagnating after two iterations with a high residual, and legacy TR
    exhausting its iteration budget. The Damped Newton 20-species real-data
    baked-vs-parameterized regression already rules out a reusable-graph
    regression for that method. Add minimal regression only for behaviour that
    changes unexpectedly (especially renewed non-finite residual), rather
    than forcing every method to be a production default.

  RST attempt metrics split residual-callback time, Jacobian-callback time,
  and remaining engine overhead; the enhanced report prints that split for
  each slow accepted point and typed failure attempt. Legacy methods
  deliberately leave this optional evidence absent rather than inventing
  incomparable counters. Do not weaken the common acceptance gate or select a
  universal default from one inventory.

- [x] Add an explicit fixed-`P,T` multi-start recovery policy. The typed
  `ResolvedPhaseEquilibriumRequest` and high-level pipeline now accept an
  ordered list of `LogMolesInitialGuess` values, reuse the prepared formulation
  and one symbolic RST problem across seeds, compare all accepted candidates
  through the common validation ordering, and publish selected-seed evidence in
  both solution bundles. Multi-start is intentionally rejected for bounded
  phase control and temperature ranges until those workflows have a separate
  lifecycle-aware policy; it is not silently applied to an incompatible outer
  loop.
- [x] Keep the live-cycle evidence boundary explicit. Real local water/ice and
  liquid fixtures cover appearance, disappearance, hysteresis retention,
  rollback, and budget termination. A physically credible continued run that
  actually revisits an active phase set has not been observed; synthetic cycle
  state-machine tests remain the mandatory correctness layer, while a live
  cycle remains deferred rather than being manufactured with an invalid policy.

## P5 - Public API, GUI, and cleanup

The original P5 checklist below is retained as historical planning context;
the fixed-P,T public API items that were completed during P4 are synchronized
with their current status here. GUI and deferred physical models remain open.

- [x] Provide a small public facade: validated problem builder, solver policy,
  solve method, solution, solve report, and optional validation report.
  A dedicated `equilibrium_public_api_tests` module imports only
  `ChemEquilibrium::prelude` and covers fixed-`P,T`, one-point range parity,
  typed phase-qualified lookup, duplicate sparse-inventory rejection, backend
  policy selection, element-selection policy, and physical phase assignment.
- [x] Expose backend selection and cascade diagnostics in GUI only through typed
  controls; never require users to type internal enum names. The editor uses
  `GuiSolverBackend` and ordered typed cascade controls; egui tests cover every
  concrete backend and duplicate rejection.
- [x] Show conservation, residual, fallback-attempt, and K_eq validation status
  as first-class result sections. The result view reads the immutable accepted
  report and does not reconstruct solver diagnostics in the GUI.
- [x] Complete GUI story tests for success, fallback success, all-backends-
  failed, invalid input, validation mismatch, and save/load roundtrip.
  - [x] Success, invalid editor input, worker failure/rollback, lifecycle, and
    document roundtrip stories are covered.
  - [x] The ignored local fallback story covers accepted backend attempts and
    renders the diagnostic sections.
  - [x] The ignored local H2/O2/H2O P,H stories cover inner fallback and a
    true engine `AllBackendsFailed` publication failure; a generic worker
    error is not used as a substitute.
  - [x] Add a deterministic accepted-result layout-mismatch fixture. The
    ignored GUI result-layer story solves two real NASA systems, combines their
    accepted snapshots as one range payload, and verifies transactional
    rejection before any partial table can be published.
- [x] Split the implementation by responsibility: `domain`, `formulation`,
  `validation`, `backend`, `solver_policy`, `keq_validation`, and `report`.
- [x] Review `easy_equilibrium.rs` as the future facade: it is rejected as the
  production facade because it models only one reaction and has no phase-aware
  acceptance contract; it is retained as a deprecated compatibility helper.
- [x] Rename the canonical modules and their tests so public paths no longer
  imply that the main formulation is only an equilibrium-constant solver.
- [x] Keep the production facade free of direct console output and duplicate
  temperature-sweep orchestration. `Untitled-1.rs` is deleted; timing output is
  confined to opt-in ignored characterization tests. The retained
  single-reaction legacy helper no longer emits `println!`/`dbg!` output.

## P6 - Production fixed-pressure, fixed-enthalpy equilibrium

This stage adds the closed-system `P,H = const` formulation after the fixed
`P,T` engine has become the canonical production core. Temperature is an
unknown result, while pressure, total elemental inventory, and total enthalpy
are prescribed.

The first implemented path is a bounded scalar temperature solve around the
existing accepted `P,T` solve:

```text
target H, pressure, temperature bracket
                    |
                    v
        evaluate F(T) = H_eq(P, T) - H_target
                    |
                    v
        canonical solve_resolved_pt at trial T
                    |
                    v
      accepted composition + phase-control report
```

This path remains valuable as an independent reference and safeguarded
fallback. It must not, however, be presented as the final production
architecture: the canonical target is a coupled composition-temperature
system with an analytic Jacobian.

### P6.A Architectural correction: monolithic P,H is the production target

Current dependency map (30.07.2026):

```text
ResolvedPhaseEnthalpyRequest
  |
  +-- ResolvedThermochemistry / EnthalpyModel
  |
  +-- solve_resolved_ph
        |
        +-- fixed declared phases:
        |     PreparedMonolithicPhRunner          [enabled, RST symbolic + analytic legacy backends]
        |
        +-- nested/reference path:
        |     solve_bracketed_temperature_from
        |       |
        |       +-- PhTrialEvaluator::evaluate(T)
        |             |
        |             +-- PreparedPhaseEquilibriumTemplate::solve_at(T)
        |
        +-- bounded monolithic path:
              PreparedPhaseControlTemplate
                |
                +-- PreparedPhaseControlRunner
                      |
                      +-- solve_monolithic_active_set_candidate
```

The fixed-`P,T` equations that must be reused rather than copied are:

- `PreparedEquilibriumProblem::residual` and
  `evaluate_equilibrium_logmole_residual` for the canonical reaction and
  element rows;
- `PreparedEquilibriumProblem::jacobian` and
  `evaluate_equilibrium_logmole_jacobian` for the analytic log-mole block;
- `ResidualScalingContract` for matching residual/Jacobian row scaling;
- `PreparedEquilibriumRunner` and `EquilibriumNonlinearBackend` as the current
  backend-cascade boundary. The dedicated P,H RST payload now carries the
  coupled `(ln(n), theta_T)` residual rather than reusing a fixed-P,T problem;
  its exact symbolic route is deliberately limited to one native coefficient
  interval per component. A range crossing stays on the analytic P,H path
  until a piecewise-symbolic contract is designed explicitly.

The active-set boundary is `PreparedPhaseControlRunner::solve`: phase
activation/deactivation, hysteresis, projection caches, and transition reports
remain outside the nonlinear fixed-active-set solve. P,H must provide one
monolithic fixed-set candidate solver to this outer loop; phase-control state
must never enter the P,H residual.

The existing P,H code is split by ownership as follows:

- `equilibrium_ph_thermochemistry.rs`:
  `ThermochemistryProvenance`, `ResolvedThermochemistry`,
  `MolarThermoFunction`, `EnthalpyModel`, `EnthalpyEvaluation`, and private
  property/capability construction;
- `equilibrium_ph_formulation.rs`:
  monolithic unknown/row layouts, bounded temperature transform, P,H
  residual/Jacobian evaluation, row scaling, and dimension validation;
- `equilibrium_ph_monolithic.rs`:
  prepared fixed-active-set P,H problem, backend cascade, candidate
  reconstruction, and monolithic diagnostics;
- `equilibrium_ph_nested.rs`:
  the current safeguarded bracket, trial evaluator, nested budgets,
  monotonicity policy, and nested evidence;
- `equilibrium_ph_workflow.rs`:
  public request/result facade, `PhSolveMode`, classified fallback policy,
  phase-control integration, and final publication validation.

Reference tests retained during migration:

- `nonreacting_sensible_heat_fixture_recovers_the_analytic_root`;
- `inventory_scaling_preserves_temperature_and_mole_fractions`;
- `sampled_non_monotone_branch_is_rejected_by_default_but_explicit_compatibility_allows_it`;
- `ph_temperature_seed_is_a_real_trial_and_can_accept_the_root`;
- `bracketed_solver_rejects_unreachable_target`;
- live `live_reactive_pt_to_h_to_ph_recovers_temperature_and_composition`;
- live `live_bounded_water_pt_to_h_to_ph_preserves_phase_evidence_and_json`;
- live `live_reactive_gas_ph_backend_matrix`.

Migration contract:

- [x] Introduce `PhSolveMode::{Monolithic, NestedTemperature, Auto}`. The
  resolved-thermochemistry facade now defaults to `Auto`: bounded monolithic
  seed recovery is attempted first and the independent nested route remains a
  final classified fallback. The generic closure constructor remains explicitly
  nested because it cannot provide the Gibbs/Cp bundle required by the
  coupled formulation; `NestedTemperature` remains available as an
  independent reference route and `Auto` as the classified recovery mode.
- [x] Make `Auto` run monolithic first and fall back solely for classified
  numerical failures. Input, thermochemistry, dimension, and
  unsupported-physics errors never trigger fallback; the immutable report
  retains the typed fallback reason and route decision, so the production
  default cannot silently hide which formulation accepted the state.
- [x] Add unknown vector `[ln(n_0), ..., ln(n_{m-1}), theta_T]` with a smooth
  bounded temperature transform. Do not clip temperature inside residual
  evaluation. `equilibrium_ph_formulation` now owns deterministic unknown and
  row layouts plus a tested logistic transform that rejects floating-point
  saturation instead of clipping.
- [x] Add residual rows `[reaction equilibrium, element balances, enthalpy]`
  while preserving the exact P,T row ordering and one explicit energy scale.
  The prepared formulation is solver-agnostic and is published through the
  fixed-declared-phase facade only after the common candidate gate accepts it.
- [x] Reuse the analytic P,T log-mole Jacobian as the upper-left block and add
  the temperature column plus enthalpy row analytically.
- [x] Require component-aligned `Cp(T)` for the production monolithic path.
  Any finite-difference fallback must be a typed thermochemistry capability
  with provenance, not an implicit numerical trick in the solver. The prepared
  monolithic formulation rejects a missing component capability explicitly.
- [x] Keep resolved `G0(T)` synchronized with the coefficient interval used by
  the canonical `P,T` bridge. `ResolvedThermochemistry` now owns a private
  per-phase, temperature-keyed Gibbs snapshot: it selects all phase
  coefficients and constructs the phase `G0` functions once per temperature,
  then shares them across component residual rows. This fixed the live
  `P,T -> H -> monolithic P,H` mismatch without mutating resolved data or JSON
  libraries.
- [x] Extend the activity contract with `d ln(a) / dT`. The current ideal-gas
  and ideal-solution models have zero derivative at fixed pressure; future
  non-ideal models must provide or explicitly decline this capability.
- [x] Compare the full analytic P,H Jacobian against central finite
  differences block by block, including the temperature chain rule.
  The initial test covers reaction, element, and enthalpy rows against all
  log-mole and temperature columns; live thermochemistry fixtures remain part
  of backend integration.
- [x] Run one monolithic solve per fixed active set and carry its accepted
  composition and temperature into the next phase-control transition.
  - [x] `equilibrium_ph_monolithic::PreparedMonolithicPhRunner` now performs
    the coupled fixed-active-set solve through the common legacy backend
    cascade, reuses the shared P,T candidate gate, and rejects the existing
    P,T-only RST symbolic payload explicitly instead of faking compatibility.
  - [x] The fixed declared public facade bridges the accepted monolithic
    snapshot into the ordinary immutable phase-aware result with the solved
    temperature, matching lookup report, and backend evidence. No mutable
    `EquilibriumLogMoles` compatibility object participates in that path.
  - [x] Let `PreparedPhaseControlRunner` invoke this fixed-set runner after
    active-set transitions. This is implemented through the existing injected
    candidate callback, not a second orchestration facade. A monolithic
    candidate may probe a trace-seeded wider mask when the current reduced
    branch has no enthalpy root; the runner records that mask expansion as an
    ordinary activation transition before final publication.
    - [x] Extract the shared active-set lifecycle behind an injected fixed-set
      candidate callback. The existing P,T runner uses that path now, while
      hysteresis, boundary recovery, cycle detection, transition reports, and
      publication remain owned by `PreparedPhaseControlRunner`.
    - [x] Make phase-stability evaluation candidate-local: it now consumes the
      accepted candidate's Gibbs capabilities and conditions rather than the
      runner construction temperature. This is required before a P,H candidate
      with solved temperature can enter the same loop.
    - [x] Add the reduced monolithic P,H candidate adapter. It projects
      both `PreparedEquilibriumProblem` and `ResolvedThermochemistry` by the
      active mask, solve `[ln(n_active), theta_T]`, scatter back to the full
      layout, and provides full-layout `G0(T_solution)` to stability checks.
- [x] Return one result type for monolithic, nested, and fallback paths, with
  the selected path, fallback reason, backend attempts, residual blocks,
  energy error, conservation evidence, and timing stated explicitly.
  - [x] The existing `FixedPressureEnthalpySolution` now exposes the actual
    `PhSolvePath` through its common report; nested and monolithic live paths
    are both asserted, including `MonolithicPhaseControl`. A classified
    fallback reason and `Auto` policy remain open.

### P6.0 Freeze the physical and dimensional contract

- [x] Introduce a typed top-level equilibrium constraint instead of boolean
  flags. Keep fixed pressure/reference pressure in the existing condition
  types and distinguish:
  - prescribed `P,T`;
  - prescribed `P,H`;
  - the `P,H` temperature seed or bracket hint;
  - the accepted equilibrium temperature;
  - the prescribed total extensive enthalpy.
  `equilibrium_constraints::EquilibriumConstraint` now provides this boundary
  and is consumed by the typed P,H facade.
- [x] Make total enthalpy in joules the canonical engine input. The total is
  tied to the supplied closed-system inventory; any GUI molar, mass-specific,
  or mixture-normalized input is an explicit conversion before entering the
  engine. `TotalEnthalpyJoules` is now the typed `PH` boundary, while the raw
  constructor remains as a transition convenience.
- [x] Document the initial supported physical scope: ideal-gas phases and pure
  condensed phases at fixed pressure, using the current standard-state
  thermochemistry. Do not imply support for non-ideal solution enthalpy,
  pressure-dependent real-fluid enthalpy, excess enthalpy, kinetic/potential
  energy, or unmodelled heat/work terms. The public module documentation now
  makes this boundary explicit before any GUI exposes the workflow.
- [x] Define the enthalpy reference convention. All component enthalpies in one
  solve must come from thermochemically compatible records and reference
  states; `H_target` must use the same convention. The public contract states
  this explicitly and `ResolvedThermochemistry` already pins every component
  to its selected record and provenance.
- [x] State precisely when `P,H` represents an adiabatic isobaric calculation:
  no unaccounted heat transfer, no shaft/electrical work, and no omitted
  kinetic or potential energy. This is now part of the public module contract,
  not a hidden assumption of the scalar solver.
- [x] Do not require source compatibility at any cost. Preserve `P,T`
  numerical behavior and provide a clear migration path, but remove or
  deprecate an old constructor if retaining it would duplicate the canonical
  condition model. Compatibility constructors that duplicate the canonical
  resolved-thermochemistry boundary are now explicitly deprecated; the narrow
  `from_resolved_thermochemistry(...)` builder remains the production entry.

### P6.1 Extend the resolved thermochemistry capability

- [x] Add the first resolved-system `dH(T)` bridge. `EnthalpyModel::from_resolved_system`
  aligns capabilities to `SystemLayout` and rebuilds temperature-dependent
  property readers in private `SubsData` copies.
- [x] Add a typed per-component thermochemistry capability aligned exactly to
  `SystemLayout`, carrying at least:
  - standard Gibbs free energy `g0(T)`;
  - molar enthalpy `h(T)`;
  - optional heat capacity `cp(T)`;
  - valid temperature interval;
  - source library, record key, phase/state evidence, and lookup provenance.
  `ResolvedThermochemistry` carries this bundle and the accepted P,H result
  retains it for diagnostics.
- [x] Reuse the existing `SubsData` `dH(T)` and `Cp(T)` calculators and
  closures. Do not reconstruct enthalpy by differentiating Gibbs free energy,
  and do not copy NASA/NIST formula implementations into `ChemEquilibrium`.
  The resolved P,H bundle is format-agnostic: it consumes the common
  `ThermoCalculator` capability API after lookup and never branches on a
  library name or polynomial representation.
- [x] Build Gibbs, enthalpy, and optional heat-capacity capabilities from the
  same selected thermochemical record. Reject mixed-record property bundles
  unless an explicit, validated reconciliation policy is introduced later.
- [x] Preserve the current read-only repository transaction: resolving or
  evaluating `P,H` thermochemistry must never mutate the JSON libraries.
- [x] Compute the admissible temperature domain as the intersection of all
  selected records' valid intervals. Reject an empty intersection before any
  nonlinear or scalar solver starts.
- [ ] Replace the current single interval/envelope representation with an
  exact ordered set of disjoint valid segments if a supported NASA/NIST
  record can contain gaps. Until then, runtime property evaluation must still
  reject a temperature that falls inside an envelope gap; the envelope is not
  permission to extrapolate across missing coefficients.
- [x] Treat enthalpy as mandatory for the first production outer-temperature
  workflow. Keep `Cp` optional there; it becomes mandatory only for an
  analytic derivative/Newton acceleration or the monolithic formulation.
- [x] Expose a pure additive derivative-ready primitive for
  `dH/dT|n = sum_i n_i Cp_i`. Its result is explicitly named a partial
  derivative and does not pretend to include the implicit equilibrium term
  `sum_i h_i d(n_i)/dT`.
- [ ] Define a future extension point for phase-model enthalpy contributions.
  The current additive contract
  `H = sum_i n_i h_i(T)` must fail explicitly for any model requiring an
  unimplemented excess/mixing enthalpy term.

### P6.2 Add the canonical outer-temperature workflow

- [x] Add the first closure-backed bracketed workflow. The
  `equilibrium_ph_workflow` module provides `EnthalpyModel`, the
  `ResolvedThermochemistry` bundle, validated scalar controls, trial/report
  types, and `solve_resolved_ph` over the canonical `solve_resolved_pt`
  facade. The bundle-backed constructor is now the resolved-data boundary;
  continuation, budgets, and phase-transition policy remain open below.
- [x] Promote the prototype into the final narrow typed `P,H` request above
  the resolved-phase facade. The production request must own:
  - immutable `ResolvedPhaseSystem`;
  - typed initial composition;
  - pressure and reference pressure;
  - target total enthalpy;
  - validated temperature bounds and optional seed;
  - fixed-`P,T` solve options and phase-control policy;
  - scalar temperature-solver policy, budgets, cancellation, and timing mode.
  `ResolvedPhaseEnthalpyRequest` now captures an owned immutable
  `ResolvedPhaseSystem` snapshot, typed composition/constraint/bounds, inner
  solver and phase policy, scalar policy, budgets, cancellation, and timing.
- [x] Define a pure trial evaluation:
  `F(T) = H(solution_of_P_T(T), T) - H_target`.
  A trial is usable only after the inner `P,T` candidate passes the existing
  backend-independent acceptance gate. `PhTrialEvaluator` now creates the
  local immutable outcome while the outer workflow alone owns budgets,
  progress, continuation, and final publication.
- [x] Use a proven safeguarded bracketed scalar method (Brent-style or
  equivalent bisection/interpolation hybrid) as the production default.
  The current default accepts a secant proposal only inside a guarded interior
  portion of the valid sign bracket and otherwise uses bisection. Trial
  reports retain `LowerBound`, `UpperBound`, `Seed`, `Interpolation`, or
  `Bisection`, so the scalar path is auditable. No unguarded Newton step is
  enabled.
- [x] Add deterministic bracket construction and validation:
  - explicit user bracket takes precedence;
  - the required `P,H` temperature seed is evaluated as an additional trial
    strictly inside, and never instead of, user bounds; it deterministically
    narrows an already valid bracket when possible;
  - [x] endpoint and midpoint trial failures retain the zero-based trial
    index, temperature, typed inner cause, source chain, and retryability
    classification instead of becoming an unqualified scalar-solver error;
  - an unbracketed or unreachable target returns a typed error; a seed that
    exposes two sign-change intervals returns `enthalpy_multiple_brackets`
    instead of selecting a root by incidental evaluation order.
- [x] Reuse expensive immutable preparation across fixed-declared-phase trial
  temperatures: `PreparedPhaseEquilibriumTemplate` now retains the resolved
  layout, element matrix, reaction basis, provenance, numeric closures, and
  optional RST symbolic problem. The P,H workflow creates it lazily at the
  first real trial, then retargets only temperature-dependent
  thermochemistry while preserving deterministic multi-start recovery. The
  report publishes one build and per-trial reuse evidence. Bounded phase
  control remains intentionally separate because its active-set lifecycle
  cannot be reused across the non-monotone temperature order of a scalar
  bracket without changing the physical contract.
- [x] Use the nearest accepted fixed-declared-phase composition as an
  additional typed log-mole continuation seed without publishing trial state.
  Bounded phase control intentionally keeps its ordinary inventory seed until
  active-set-aware continuation is implemented; no active set is assumed to
  remain unchanged across a phase transition.
- [x] Add one global scalar temperature-evaluation budget and enforce it in
  the bracket helper itself, so endpoint and interior evaluations share the
  same limit rather than receiving a fresh budget per trial.
- [x] Forward `EquilibriumExecutionControl` into the outer P,H workflow and
  every inner P,T request, with explicit temperature-trial progress stages.
- [x] Enforce an optional wall-time budget at the scalar evaluation gate. The
  budget is global to the outer operation; inner cooperative cancellation is
  still required for interruption during a long backend call. The deadline is
  now checked after evaluator return as well, so a trial that finishes after
  the deadline cannot be accepted retroactively.
- [x] Make the complete `P,H` solve transactional. Failed trial solves,
  rejected phase transitions, cancellation, exhausted budgets, or a failed
  final acceptance check must not publish a partial result or mutate the
  caller's resolved system.
- [x] Define nested resource budgets explicitly: scalar evaluations, total
  inner backend attempts, total nonlinear iterations/evaluations, wall time,
  and phase-control transitions. A per-trial budget must not accidentally
  multiply into an unbounded outer budget. The scalar evaluation budget is
  global, and `PhTemperatureSolveOptions` now adds optional global
  `max_inner_backend_attempts`, `max_inner_nonlinear_iterations`, and
  `max_phase_control_transitions` limits.
  The outer transaction accounts for every started inner backend attempt and
  every reported nonlinear iteration before publishing its trial result,
  including all work belonging to every continuation multi-start seed. Phase
  transitions are counted from the immutable solution report and are subject
  to the same global outer budget.
  The request owns an immutable resolved snapshot, while accepted trials and
  the final result remain local until the scalar bracket and enthalpy contract
  succeed.
- [x] Integrate cooperative cancellation and progress events at both levels:
  bracket preparation, temperature trial start/accept/reject, inner backend
  attempt, phase transition, and final publication. Trial-level cancellation
  and progress are wired, including a typed rejected-trial event for inner
  P,T or enthalpy-evaluation failures. The P,H report now also emits inner
  backend start/finish, accepted transition, and publication start/finish
  events; cancellation is checked at each publication boundary.

### P6.3 Handle phase transitions and non-smooth enthalpy honestly

- [ ] Do not assume `H_eq(P,T)` is globally smooth or monotone. Phase
  appearance/disappearance can create derivative discontinuities, and
  coexistence can create very flat or discontinuous-looking numerical
  branches.
  - [x] The safeguarded outer solver now has an explicit
  `PhMonotonicityPolicy`. Its default rejects an observed reversal in the
  sampled enthalpy branch with a typed error; the legacy sign-bracket
  behavior is available only through an explicit compatibility policy and
  is recorded in the immutable solve report. This is sampled evidence, not
  a claim of global monotonicity.
  - [x] Every accepted temperature trial now retains the exact phase-id to
    lifecycle-status snapshot accepted by its inner P,T solve. Consumers can
    therefore compare adjacent active sets instead of inferring a transition
    only from aggregate counters.
- [ ] Define how the scalar solver treats active-set changes inside a bracket.
  Retain transition evidence per trial and prevent interpolation steps from
  treating values from incompatible failed branches as a smooth derivative.
  - [x] The bounded phase-control P,H path now forces bisection. Safeguarded
    secant interpolation remains available only for fixed declared phases,
    where the phase layout is invariant. This prevents a numerical slope from
    being inferred across two potentially different active sets.
  - [x] Every published trial now exposes a typed preparation reason:
    `FixedFormulationInitial`, `FixedFormulationReused`, or
    `BoundedPhaseControlIsolated`. The latter makes the absence of active-set
    reuse explicit rather than looking like a missed performance metric.
- [ ] Detect and report multiple brackets/multiple roots when sampled evidence
  reveals them. The first production policy must be deterministic (for
  example, nearest valid root to the requested seed), not dependent on hash or
  backend iteration order.
  - [x] The mandatory interior seed already detects the first concrete
    multiple-bracket evidence: equal-sign endpoints with an opposite-sign seed
    are rejected as `enthalpy_multiple_brackets`. Broader scan-based root
    enumeration remains intentionally deferred until phase-boundary semantics
    are finalized.
- [ ] Define boundary behavior for latent-heat/coexistence cases. If the
  requested enthalpy is achieved by changing phase fractions at nearly fixed
  temperature, the accepted result still must satisfy composition,
  complementarity, conservation, and enthalpy tolerances.
  - [x] Add a release-only real gas/ice `P,T -> H -> P,H` boundary fixture.
    The monolithic active-set candidate exhausts its backend cascade, but
    explicit `Auto` deterministically recovers through bounded nested `P,T`
    phase control, retains `AllBackendsFailed` as fallback evidence, activates
    the solid phase, preserves balances/enthalpy, and leaves local libraries
    unchanged.
  - [ ] Define a true phase-fraction/coexistence formulation only after the
    physical model and complementarity contract are specified. The successful
    ice recovery is evidence for one point solution, not proof that every
    latent-heat plateau is represented correctly.
- [x] Reuse the bounded phase-control hysteresis and cycle/budget protections.
  The bounded monolithic P,H adapter now runs through
  `PreparedPhaseControlRunner`, so hysteresis, transition limits, cycle
  detection, rollback, and typed lifecycle reports are shared with the P,T
  path. A scalar P,H trial cannot hide an inner phase-control cycle as a
  generic scalar-function failure; the accepted water fixture also verifies
  trace-seeded activation of a zero-inventory liquid phase.

### P6.4 Publish an immutable `P,H` solution and diagnostics

- [x] Return an immutable result that contains the accepted temperature,
  pressure, component moles, phase states, target/calculated total enthalpy,
  raw enthalpy error in joules, relative/scaled error, and the normal
  fixed-`P,T` acceptance reports. Bundle-backed results also retain the
  component-aligned thermochemistry provenance and common domain.
- [x] Add a validated enthalpy acceptance contract with both absolute and
  scale-aware relative tolerances. The scalar solver accepts only when
  `|H-H_target| <= max(abs_tol, relative_tol * scale)`; near-zero targets
  therefore retain a meaningful absolute floor.
- [x] Define one positive finite enthalpy scale from the immutable request and
  initial state, not from arbitrary solver iterates. Record the scale and both
  tolerance components in the report so acceptance is reproducible.
- [x] Retain complete nested evidence:
  - bracket endpoints and accepted root;
  - every temperature trial and its status;
  - inner backend attempts and solver metrics;
  - phase transitions and active-set changes;
  - rebuild/reuse reasons;
  - thermochemistry provenance;
  - timing by repository, property refresh, preparation, nonlinear solve,
    phase control, enthalpy evaluation, scalar orchestration, and
    postprocessing.
  Every accepted trial now retains its complete inner backend cascade,
  optional multi-start comparison, phase-control lifecycle, and final
  acceptance snapshots alongside compact counters and timing. Rejected scalar
  trials remain typed errors rather than partial result rows. Each published
  trial now also carries an immutable `PhTrialTimingReport` separating trial
  wall time, nested `P,T` time, and additive enthalpy evaluation. The timing
  report is opt-in and covered by the P,H unit evidence test. Each inner
  `PhaseTransitionRecord` now also carries a measured control-pass duration;
  live activation/deactivation stories verify that this timing is published.
  - [x] Publish a point-level `formulation_build` duration for typed T-range
    points. It accounts for reduced formulation and RST symbolic construction
    performed during that point and is covered by the live ice-transition
    range story.
   - [x] Preserve separate stopwatch intervals for every individual formulation
     cache entry when one point creates multiple active-set formulations. The
     public point duration remains an aggregate, while each accepted range
     point now carries a deterministic active-mask/cache-entry timing snapshot.
  - [x] Per-trial and aggregate backend-attempt, nonlinear-iteration, and
    phase-transition counters are now published together with nested timing;
    trial/evidence cardinality is validated before report publication.
  - [x] Fixed-declared-phase P,H reports now state the number of immutable
  formulation builds and retarget reuses. Bounded phase control deliberately
  reports no such reuse because rejected scalar trials cannot safely mutate
  or seed its active-set lifecycle.
  - [x] `Auto` now preserves both typed failure trees when its monolithic
    attempt and nested recovery both fail. Successful recovery retains the
    classified fallback reason in the immutable P,H report; failed recovery
    returns `PhAutoFallbackFailed { monolithic, nested }` rather than erasing
    the primary active-set evidence.
  - [x] Keep the default production cascade free of unsolicited console
    output. RustedSciThe Powell Dogleg remains an explicit benchmark/backend
    option, but is excluded from the default cascade because version 0.4.12
    prints internal `beta` diagnostics and rejected the large real-data case.
- [x] Keep independent `K_eq` validation scoped to applicable inner
  fixed-phase chemical-equilibrium candidates. It does not independently
  validate the outer enthalpy root and must not be presented as doing so. The
  P,H module documentation now explicitly separates this observational inner
  evidence from enthalpy-root acceptance.
- [x] Add a compact stable public facade and prelude exports only after the
  request/result/error contracts are validated. Keep raw prepared state and
  scalar orchestration internals crate-private. `solve_resolved_ph`, its typed
  request/result/report types, and constrained solver controls are exported
  from `ChemEquilibrium::prelude`; prepared runner/template internals remain
  crate-private.

### P6.5 Verification and release evidence

- [x] Add a typed fixed-pressure, fixed-enthalpy target-range facade with a
  strictly monotone ascending/descending enthalpy grid, continuation from the
  previous accepted physical composition and temperature, per-point timing,
  phase-transition evidence, and transactional publication. `PhRangeRequest`
  is the canonical engine boundary; it does not publish a partially solved
  batch after a point failure.
- [x] Cover the continuation contract with unit tests for grid/error/report
  semantics and an ignored live NASA H/O story that checks ascending and
  descending targets, accepted temperatures, provenance, conservation, and
  the seed hand-off between adjacent points.
- [x] Reuse prepared monolithic P,H formulation state across fixed-phase target
  points. Reaction basis, element totals, phase projection, analytic row
  structure, and RST symbolic expressions are built once; accepted
  `[log-moles, T]` state plus target/scale parameters are retargeted between
  points. The range report exposes formulation build/reuse counters.
- [x] Extend the same prepared-state contract to nested P,H batches without
  conflating it with the monolithic formulation. The nested route now shares
  only its fixed-P,T prepared template across target batches; every target
  retains an independent scalar bracket, while bounded phase-control remains
  isolated whenever active-set lifecycle can change. The live continuation
  story asserts one template build across the nested target range.
- [x] Add a live backend matrix for the fixed-phase monolithic P,H target
  range with strict single-backend policies and per-point timing. The ignored
  matrix keeps backend failure visible instead of hiding it behind a cascade.
- [x] Extend the live route matrix with nested/Auto batches, a real
  phase-transition range, and rollback after an unreachable target; these are
  separate lifecycle stories and must not be implied by the fixed-phase
  matrix. The ignored real-data route matrix now prints all three rows and
  keeps the rollback point index typed. Its release evidence is recorded in
  `STORY_TESTS.md`.
- [x] Add a large real-data nested P,H continuation story: twenty exact
  element-limited C/H/O species, nine interior target enthalpies, accepted
  temperature/composition hand-off, one fixed-P,T template build, conservation
  and enthalpy acceptance at every point.
- [x] Add a denser real water/ice Auto P,H story: five target enthalpies from
  250 to 270 K, route-dependent monolithic/nested acceptance, retained
  fallback reason, phase-transition evidence, and transactional output.
- [x] Record release evidence for the nested/Auto route matrix, the large
  twenty-species nested story, and the dense water/ice Auto story in
  `STORY_TESTS.md`. Keep these large live tests out of the default suite.
- [x] Implement and document the strict fixed-phase monolithic P,H target-range
  backend matrix. It remains separate because backend failures are deliberate
  characterization evidence rather than a route-lifecycle failure; the release
  command and aligned backend table are recorded in `STORY_TESTS.md` for the
  target-machine evidence refresh.
  - [x] The matrix now emits one aligned row per backend with point count,
    formulation builds/reuses, solve and wall timing, maximum residual,
    maximum element-balance error, and the complete typed failure reason.
    A debug characterization currently shows \`legacy_tr\` completing all
    three points while the other rows expose backend-specific acceptance or
    iteration-limit failures; this is evidence, not a production-default
    decision.

- [x] Keep every existing fixed-`P,T` regression green and add direct parity
  between the resolve/build pipeline request and the resolved `P,T` facade.
  The parity story fixes one backend policy, compares component amounts,
  residual quality, and layout fingerprint, and therefore checks orchestration
  equivalence without conflating it with backend-cascade selection.
- [x] Add the primary inverse story:
  solve a stable `P,T` problem at `T*`, calculate `H*`, then solve the same
  inventory at `P,H*` from a different seed/bracket and recover temperature,
  composition, conservation, and enthalpy within explicit tolerances.
  The ignored live NASA H2/O2/H2O story now exercises this contract through
  `ResolvedThermochemistry::from_resolved_system` and `solve_resolved_ph`.
- [x] Cover synthetic analytic fixtures where `h(T)` and the expected
  temperature/root are known exactly. Include a non-reacting sensible-heat
  case so failures in chemistry cannot mask errors in the energy equation;
  the same unit layer also checks inventory scaling invariance.
- [ ] Add real offline cross-format gas and gas-plus-pure-condensed stories
  with pinned lookup policy and provenance. Snapshot all canonical JSON files
  before and after the tests. The solver contract is library/polynomial-format
  agnostic, so this must include at least one stable non-NASA fixture rather
  than treating NASA as the production format.
  - [x] Ignored release stories now cover NASA-gas reactive P,T -> H -> P,H
    inversion and NASA-gas plus NASA-condensed bounded water inversion. Both
    retain provenance, verify the P,H reuse/phase evidence appropriate to
    their solve mode, and compare the canonical library files byte-for-byte.
  - [ ] Add a pinned local NIST P,T -> H -> P,H inverse fixture once the
    read-only repository actually contains a stable NIST record set. The
    current `all_keys_substance.json` advertises historical NIST addresses,
    but the canonical local repository does not contain the corresponding
    H2/O2/H2O records, so an offline test would only pin stale-index failure.
    Repair the library/index consistency first; do not turn this into a
    network-dependent regression. This fixture is evidence for the
    format-agnostic capability boundary, not an alternate solver path. See
    **F1** for the required atomic data-release and full offline matrix.
  - [x] Make the fallback contract explicit before adding that fixture:
    local data is authoritative; an enabled canonical fallback queries NIST
    only for the requested gas/liquid/solid state; an unspecified state fails
    closed; and parser/network failures are preserved as failures. The old
    unconstrained gas fallback remains compatibility-only and deprecated.
  - [x] Add an ignored online NIST state-matrix smoke test for the parser
    boundary. It may validate gas/liquid/solid payload availability and finite
    Shomate data, but it must not be used as equilibrium evidence or mutate
    local JSON. Offline NIST parity remains blocked until real payloads are
    checked into the repository.
    - [x] The current H2O run reports complete gas/liquid payloads and an
      incomplete solid navigation payload without Cp intervals. The latter is
      recorded as unavailable solid thermochemistry, never as a phase fallback.
  - [x] Add a read-only catalog consistency report before repairing that
    fixture. It canonicalizes library aliases and reports missing payload,
    orphan payload, and duplicate index pairs without modifying JSON. The
    live release diagnostic currently reports `5519` indexed pairs, `5472`
    payload pairs, and `47` indexed records without payload; this is now an
    explicit data-release blocker rather than a hidden solver failure.
- [x] Add phase-transition stories for water vapor/liquid/ice and another
  physically credible condensed system. Check latent-heat/coexistence
  behavior, hysteresis, rollback, and final complementarity.
   - [x] Real P,T water liquid/ice appearance and disappearance stories are
     covered; the ignored P,H gas/ice inverse story verifies bounded nested
     recovery after a rejected monolithic candidate. The release rerun
     recovered `T=250 K`, `solid_moles=4.998093e-1`, one transition, and zero
     scaled enthalpy error through `Auto`.
   - [x] The ignored real phase-transition release matrix now reports water
     gas/ice, water gas/liquid, hot-water disappearance, and graphite
     appearance/disappearance in one table. Every row checks finite
     non-negative amounts, conservation, complementarity, transition evidence,
     and byte-for-byte JSON immutability. The release run passed for all five
     scenarios, with per-case solve times from roughly 1.1 to 5.0 ms.
   - [x] The real ice temperature-range story now retains one deterministic
     build-duration snapshot per prepared active-set cache entry, instead of
     exposing only the aggregate point duration.
- [x] Test invalid and unreachable inputs: non-finite enthalpy, invalid
  pressure/temperature bounds, missing enthalpy capability, incompatible
  record intervals, unsupported phase enthalpy model, unbracketed target,
  inner all-backends-failed, cancellation, and exhausted global budget.
  - [x] The P,H unit layer already covers non-finite targets, invalid scalar
    budgets, missing/misaligned `Cp`, wrong enthalpy component counts,
    incompatible bundle domains, unbracketed targets, observed multiple
    brackets, endpoint/midpoint inner failures, cancellation before and after
    inner start, and global evaluation/wall-time budgets. These paths return
    typed errors and never publish a partial `FixedPressureEnthalpySolution`.
  - [ ] Add the remaining resolved-data cases only when their physical source
    exists: an actual record without `dH`, an unsupported non-additive phase
    model, and a repaired local NIST fixture. Do not fabricate malformed
    production records merely to make this checklist look complete.
  - [x] The ignored real-data validation matrix covers non-finite target,
    reversed bounds, non-positive pressure, unreachable target rollback, and
    JSON immutability.
- [x] Add metamorphic inventory-scaling tests. Multiplying every initial mole
  and total target enthalpy by the same factor should preserve equilibrium
  temperature and mole fractions while scaling extensive amounts. The
  analytic outer-solver fixture now pins this contract independently of
  chemistry.
- [x] Add a backend/cascade matrix for the inner `P,T` solves and a release
  characterization over small, medium, and large real systems. Report outer
  evaluations separately from measured inner nonlinear solve time; do not
  select a production default from wall time alone.
  - [x] The ignored live H2/O2/H2O `P,H` backend matrix runs each RST and
    legacy backend under a strict `Single` policy against one independently
    constructed `H*`. Its rounded table reports outer trials, total/wall and
    inner nonlinear timings, formulation build/reuse counts, inner work,
    final temperature, energy error, residual, balance, and an explicit error
    row for failed methods. It snapshots canonical JSON libraries before and
  after the complete matrix.
   - [x] The strict monolithic P,H target-range matrix now emits aligned
     per-backend timing, attempts, accepted-backend, enthalpy, residual,
     conservation, failure-point, failure-kind, and full typed-error evidence.
     Release execution remains the final characterization step.
   - [x] `live_nested_ph_inner_pt_cascade_release_matrix` now covers 5, 20,
     and 100 real local NASA candidates over three target enthalpies. It reports
     outer evaluations, inner backend attempts, inner solve time, formulation
     reuse, conservation, and immutable-library evidence. The debug baseline
     passes; the release command is recorded in `STORY_TESTS.md`.
- [x] Verify deterministic ascending/descending enthalpy sweeps in the typed
  P,H batch facade. Reuse accepted neighboring solutions transactionally,
  while keeping the semantic distinction explicit: one P,H target has one
  solved temperature.

### P6.6 Evaluate and promote the monolithic formulation

- [x] Only after the outer-temperature workflow was accepted, prototype a
  coupled unknown layout with log-moles and bounded/logarithmic temperature.
  The canonical reaction and element rows remain embedded in the prepared
  problem; no second multiplier-based formulation was introduced. Keep the
  layout typed; do not spread arithmetic such as `n_species + n_elements`
  through residual code.
- [x] Reuse the existing chemical and balance blocks and append exactly one
  scaled enthalpy equation. Do not create a second copy of the `P,T`
  formulation.
- [x] Derive and test the full temperature column of every chemical residual,
  including `g0(T)/(RT)`, pressure/activity terms, and any phase-model
  temperature dependency. A zero temperature column is invalid.
- [x] For additive ideal enthalpy, verify the analytic energy derivatives:
  `dH/dln(n_i) = n_i h_i(T)` and
  `dH/dln(T) = T sum_i n_i cp_i(T)`.
  Extend these formulas before enabling any excess-enthalpy model.
- [ ] If analytic `Cp` or activity derivatives are unavailable, isolate a
  scale-aware finite-difference fallback and compare the complete Jacobian
  against central differences away from phase boundaries. Do not replace the
  established analytic composition Jacobian wholesale.
- [x] Bound temperature through a validated transform or safeguarded step
  policy over the common thermochemistry interval; `ln(T)` alone enforces
  positivity but does not enforce the upper/lower data bounds.
- [ ] Compare monolithic and outer-temperature formulations on convergence
  basin, phase transitions, backend fallback behavior, residual quality,
  thermochemistry evaluations, and release timing. The real reactive inverse
  story now compares both paths. Resolved requests default to the auditable
  `Auto` policy: monolithic with bounded fresh temperature starts first, then
  nested recovery only after a classified numerical failure. A real
  water/liquid phase-activation story now directly
  compares explicit monolithic and nested routes with the same `P,H` target,
  phase lifecycle evidence, accepted temperature, composition, residual, and
  conservation. The real gas/ice story now additionally proves the hard-route
  contract: direct monolithic returns `AllBackendsFailed`, explicit nested
  accepts, and `Auto` publishes the nested-equivalent state while retaining the
  same fallback reason. The broader phase-transition/release matrix remains
  before declaring promotion complete.

### P6.7 GUI integration gate

- [x] The fixed-phase ideal/pure-condensed engine facade, immutable reports,
  typed errors, cancellation, budgets, and release backend evidence are now
  sufficient to begin GUI integration within the documented P6.0 scope. GUI
  must not advertise non-ideal phases, latent-heat coexistence, or an
  enthalpy-temperature range sweep as implemented capabilities.
- [x] Expose total enthalpy in explicit joules, temperature seed and bracket,
  inner solver policy, progress, cancellation, immutable energy evidence, and
  route-dependent diagnostics through typed GUI controls. Nested routes expose
  scalar trial rows; monolithic routes expose backend attempts, accepted
  backend, acceptance, and phase-control evidence without fabricated trials.
  The GUI selects the canonical monolithic route when the resolved
  thermochemistry supports one exact symbolic interval; a bracket crossing a
  native coefficient switch remains on the nested numeric reference route.
- [x] Do not offer a temperature-range control in `P,H` mode. A future batch
  operation should sweep target enthalpy and report the solved temperature at
  each point.
- [x] Add the basic GUI stories: valid local P,H solve, immutable energy and
  outer-trial diagnostics, lifecycle staleness, cancellation boundary, and
  document roundtrip.
- [x] Add GUI stories for unreachable target, inner fallback success,
  all-backends-failed, cancellation/rollback during an active P,H worker, and
  a real phase-transition publication using stable local fixtures.
  - [x] The ignored local GUI story now covers an unreachable finite P,H target:
    the worker fails visibly, publishes no partial snapshot, and leaves the
    canonical JSON files unchanged.
  - [x] Existing local GUI stories cover fallback acceptance, transactional
    range failure, cancellation/stale-result rejection, and real water/ice
    phase publication.
  - [x] All-backends-failed is a separate ignored local fixture and asserts
    the typed engine error plus absence of a partial snapshot; it does not
    manufacture a generic worker failure solely to satisfy the checklist.

## P6.8 Release-oriented monolithic P,H architecture hardening

The external architecture review is accepted as a follow-up plan, with the
following corrections to its scope. The coupled P,H mathematics and the
nested reference route already exist; this work must preserve the P,T facade,
the phase-control boundary, and the existing typed result contract. The GUI
now covers the basic P,H editor, worker, immutable energy/outer diagnostics,
trial table, cancellation boundary, and document lifecycle. The remaining
items below are engine architecture and release evidence, not a reason to
rebuild the GUI.

### P6.8.1 Dependency map and module boundaries

- [x] Produce a checked dependency map for `equilibrium_ph_workflow.rs`,
  `equilibrium_ph_formulation.rs`, `equilibrium_ph_monolithic.rs`, and the
  nested implementation. Record every shared type, thermochemistry adapter,
  `GibbsFn` construction site, `Cp` capability check, option/report type, and
  fallback decision before moving code. The map above is now synchronized
  with the extracted nested/thermochemistry modules and the shared typed
  nonlinear-system adapter.
- [x] Move `ResolvedThermochemistry`, its provenance, molar property
  function types, interval intersection, and Gibbs/enthalpy/Cp evaluation
  into a lower-level thermochemistry module. Formulation and runners may
  depend on this module; it must not depend on workflow orchestration,
  publication, phase lifecycle, GUI progress, or scalar bracket policy. The
  workflow retains only compatibility re-exports and the nested enthalpy
  adapter.
- [x] Extract the nested scalar-temperature algorithm and its trial types,
  bracket policy, monotonicity policy, continuation, budgets, and
  nested-specific report into `equilibrium_ph_nested.rs`. Keep the workflow
  as a facade, mode selector, fallback coordinator, and common publication
  boundary.
  - [x] Extract the stateless safeguarded interpolation, wall-time guard,
    sampled-monotonicity helpers, and the scalar bracket loop with unit tests.
    The module also owns `PhTemperatureTrial`, phase-state snapshots, timing,
    and inner-solve evidence. The workflow re-exports those types, adapts
    scalar records, attaches phase evidence, and publishes the final common
    result.
- [x] Do not introduce a second thermochemistry bundle or a second public
  result merely to perform this extraction. Reuse the existing immutable
  bundle and expose accessors where the current result already owns the data.
  The nested module owns scalar trial/evidence snapshots only; the common
  `FixedPressureEnthalpySolution` remains the single published result.

### P6.8.2 Typed thermochemistry evaluation

- [x] Remove every `unwrap_or(f64::NAN)` conversion from the new monolithic
  P,H path. An out-of-domain evaluation, missing coefficient, database error,
  or property evaluation failure must retain its typed cause through
  preparation/residual evaluation and must not be reported as an anonymous
  numerical backend breakdown. The remaining infallible `GibbsFn` crossing is
  an explicit, finite snapshot adapter at the accepted temperature.
- [x] Choose one explicit fallible boundary for P,H: either a fallible
  prepared-property/evaluation context or residual helpers that consume
  already evaluated Gibbs/enthalpy/Cp vectors. Do not change the established
  P,T closure API wholesale solely to introduce a new closure alias. The
  monolithic residual already consumes the typed evaluated context; phase
  control remains a documented compatibility boundary until its callback is
  made fallible.
- [ ] Replace the snapshot compatibility adapter with a fallible stability
  evaluator once phase-control no longer requires `GibbsFn = Fn(T) -> f64`.
- [x] Move mandatory `Cp(T)` capability validation to formulation preparation.
  Validate component count/order, common temperature domain, finite initial
  Gibbs/enthalpy/Cp values, and the selected capability policy before a
  nonlinear backend starts. Missing Cp is now an `InvalidProblem` before
  backend execution; a deliberately enabled finite-difference policy remains
  future work.
- [ ] If finite-difference Cp is supported, add explicit analytic/central,
  forward, and backward source diagnostics, boundary-aware step selection,
  and a Jacobian comparison test. Do not enforce an unconditional positive-Cp
  rule without first defining its physical scope.

### P6.8.3 Separate solve policies and diagnostics

- [x] Separate common physical/acceptance/execution options from monolithic
  nonlinear options and nested scalar/bracket options. Nested-only controls
  must not affect monolithic solves, and monolithic-only controls must not be
  silently applied to nested recovery.
  - [x] `PhAcceptanceOptions`, `PhMonolithicOptions`, and `PhNestedOptions`
    are typed route contracts; `PhTemperatureSolveOptions` exposes validated
    projections into them while retaining the compatibility builder surface.
  - [x] The fixed-active-set monolithic runner now receives only an immutable
    enthalpy acceptance contract; scalar bracket limits, phase-transition
    budgets, and progress controls remain in the workflow. The nested scalar
    engine now receives its own smaller `NestedBracketOptions` contract.
- [x] Separate monolithic diagnostics from nested diagnostics. A monolithic
  Newton iteration must not be represented as an outer temperature trial;
  nested reports must retain trials, bracket evidence, inner P,T work, and
  scalar residuals. Shared counters may remain in a common immutable result
  only when their semantics are identical. Monolithic backend and phase
  lifecycle evidence now lives in `PhMonolithicEvidence`, while `trials` is
  populated only by the nested scalar route.
- [x] Isolate the compatibility constructor
  `ResolvedPhaseEnthalpyRequest::new(...)` from the canonical API. It is now
  deprecated, explicitly documented as nested/reference-only, and the
  production `from_resolved_thermochemistry(...)` builder no longer calls it.
- [x] Audit `Auto`: retry nested only for classified retryable numerical
  failures; never retry for invalid input, missing capability, layout/data
  errors, cancellation, unsupported phase models, or invariant failures.
  Preserve the original monolithic error tree and test the classification
  table explicitly. `AllBackendsFailed` is retryable only when its non-empty
  trace consists solely of started numerical failures; rejected candidates,
  skipped attempts, and empty traces fail fast.
- [x] Keep phase-control lifecycle outside the monolithic residual and keep
  the nested route as a reference/fallback. Do not move hysteresis or phase
  transitions into the coupled equation itself without a new physical
  complementarity contract.

### P6.8.4 Backend and publication contracts

- [x] Remove the monolithic runner's dependence on `EquilibriumLogMoles` as a
  mutable orchestration host. Reuse neutral residual/backend helpers only;
  preparation, iteration, phase lifecycle, and result publication must use
  immutable prepared state and local candidates. The backend cascade is now a
  stateless adapter function; the old associated method is only a legacy
  wrapper.
- [x] Define a domain-independent prepared nonlinear-system contract that can
  represent both fixed P,T and P,H residual/Jacobian dimensions. RST must
  receive residual/Jacobian capabilities without knowing which coordinate is
  temperature or which row is enthalpy. The adapter now carries
  `PreparedNonlinearSystem` with dimension, residual, Jacobian, feasibility,
  and optional backend payload without attaching thermodynamic meaning to any
  coordinate.
  - [x] Add a dedicated symbolic monolithic `P,H` payload for RST. It uses the
    exact resolved component order, substitutes the bounded `T(theta)`
    expression into phase-qualified `G0(T)` and `H(T)`, appends the scaled
    enthalpy row, and gives RST the complete `(N + 1)` system. The analytic
    formulation remains an independent legacy/fallback route; no fixed-P,T
    payload or fake zero temperature column is reused.
  - [x] Keep symbolic thermochemistry optional in `ResolvedThermochemistry`.
    Numeric resolved data remain valid for legacy and nested P,H solving.
    Real symbolic capabilities are materialised lazily from private
    phase-local `SubsData` copies for the requested P,H bounds and are
    accepted only when every component remains inside one native NASA/NIST
    coefficient interval. Selecting RST without `G0/H` expressions, or across
    a native polynomial boundary, fails with a typed unsupported-capability
    error instead of fitting/extrapolating silently or changing lookup
    semantics.
- [x] Retain one physical P,H publication gate for both paths: finite and
  in-range temperature, finite log-moles/reconstructed moles, reaction and
  element residuals, enthalpy error, conservation, phase stability, and
  active-set consistency. Diagnostics may differ, acceptance meaning may not.
- [x] Refresh the immutable build report at the accepted monolithic
  temperature before publication. A fixed-active P,H solve is prepared at its
  seed temperature but may converge elsewhere; the report and accepted
  solution must carry the same `(P,T)` conditions. The live GUI monolithic
  story covers this publication contract.
- [x] Add a typed capability error that names the requested incompatible
  backend and the supported alternatives. The monolithic RST guard now uses
  `UnsupportedBackendCapability` both for absent symbolic capabilities and
  for a requested symbolic range crossing a native coefficient boundary. The
  error is non-retryable and therefore cannot be hidden by fallback.

### P6.8.5 Architecture tests and release characterization

- [x] Add architecture tests proving that monolithic mode performs no outer
  nested temperature trials, nested mode reports real scalar trials, and
  `Auto` records monolithic-then-nested fallback only after a retryable
  failure. Include an invalid-input case proving nested is not started.
  - [x] The real-data P,H story now asserts that monolithic reports no scalar
    trials and exposes `PhMonolithicEvidence`; the unit contract also records
    that invalid monolithic input emits no `TemperatureTrialStarted` event.
- [x] Add a route-dependent GUI diagnostics snapshot. The monolithic GUI view
  now presents backend-attempt summaries, accepted backend, acceptance rows,
  phase-control rows, and inner timing; it does not render a zero-valued
  temperature-trial section. Nested routes retain their real trial table.
- [ ] Extend the analytic Jacobian matrix with missing-Cp, typed property
  error, analytic/finite-difference Cp, zero activity-temperature derivative,
  non-zero activity-temperature derivative, lower/upper interval boundary,
  scaling, and every P,H block-column/row comparison.
  - [x] Add a synthetic RST P,H regression that verifies the symbolic
    `(ln(n), theta_T)` system converges through the normal candidate gate and
    never reuses the fixed-P,T symbolic payload.
  - [x] Add an ignored real-data NASA-gas parity case for RST-monolithic P,H
    against the independent analytic/legacy route. The 1900..2900 K fixture
    stays within the common native interval and compares accepted temperature,
    composition, conservation, and enthalpy rather than backend iteration
    counts. A companion real-data guard proves that 900..1100 K is rejected
    for the RST symbolic path when it crosses a coefficient switch, while the
    JSON repository remains unchanged.
   - [x] Add an ignored real-data NASA-gas Jacobian central-difference matrix.
     The fixture uses five locally resolved C/H/O records and points immediately
     inside both common temperature boundaries plus an interior point. It
     compares every analytic P,H block entry and verifies that the JSON source
     files remain unchanged. The test intentionally keeps ideal-model activity
     temperature derivatives at their physical zero value.
    - [x] Extend that real matrix across four inventory scales from `1e-6` to
      `1e3`, while retaining both native-interval boundaries and the complete
      residual/Jacobian block comparison.
    - [x] Release rerun confirms the same real Jacobian matrix: 432 entries,
      four inventory scales, three boundary/interior temperatures, and a
      six-dimensional coupled formulation.
  - [ ] Add the corresponding ignored NIST P,H parity fixture after a stable
    locally bundled NIST record set is selected. Do not weaken the exact
    single-native-interval contract to make this fixture pass.
- [x] Add one release characterization over the same real problem for
  monolithic and nested paths. Report backend attempts, residual/Jacobian and
  thermochemistry evaluations, active-set work, temperature trials, total
  inner P,T solves, reuse, and wall time. Keep performance ratios out of unit
  test assertions; use an ignored release benchmark instead. The ignored live
  P,H story now prints a table for nested, analytic monolithic, symbolic
  RST-LM monolithic, and Auto with backend metrics plus
  thermochemistry-preparation timing.
- [x] Add the remaining GUI stories when their engine fixtures exist:
  unreachable target, P,H inner fallback, all-backends-failed,
  cancellation/rollback, and a real phase-transition publication. The stories
  now exist in `equilibrium_gui_tests`; catalog-dependent cases remain
  explicitly ignored and preserve JSON immutability checks. The mixed
  gas/condensed phase-publication story selects Legacy NR explicitly because
  the symbolic monolithic route requires one native coefficient interval per
  component. The basic valid P,H, diagnostics, lifecycle, and
  document-roundtrip stories are covered.

### P6.8.6 Documentation and explicit non-goals

- [x] Update module rustdoc after extraction: workflow as facade, monolithic
  as fixed-active-set coupled solve, nested as reference/recovery algorithm,
  thermochemistry as capability/provenance layer, and formulation as the
  complete `[F_rxn, F_elem, F_H]` system with its block Jacobian. The module
  headers now describe these ownership boundaries and explicitly distinguish
  outer scalar trials from monolithic coupled evidence.
- [x] Define and implement typed P,H temperature-range batch solving with
  accepted-state continuation, transactional publication, and per-point
  evidence. Export/import, latent-heat coexistence modeling, and non-ideal
  phase physics remain outside this technical hardening pass until their
  separate contracts are defined.
- [x] Re-run the full P,T regression suite after every logical extraction;
  the current full `ChemEquilibrium` library run passes after the option and
  prepared-system extractions, and no P,T compatibility workaround became a
  hidden dependency of the monolithic P,H path.

## P7 - Replace phase stability with reacting-system TPD

This is a correctness replacement for the current ideal phase-control path,
not a future non-ideal extension. The contracts in `task.md` and
`arch_contract.md` were compared with the implementation on 24.08.2026.

### P7 review verdict

- [x] **Confirmed:** the present implementation already evaluates canonical
  chemical potentials as `mu_i = g_i^0(T) + R*T*ln(a_i)` through
  `PhaseActivityModel::log_activity`. Preserve and reuse that boundary; do not
  introduce a second activity implementation inside phase stability.
- [x] **Confirmed:** the current stability criterion is not general. It
  special-cases one-component condensed phases, rejects a multicomponent
  `IdealSolution`, and rejects a rank-deficient active elemental assemblage.
  Those restrictions are implementation artifacts rather than the required
  reacting-system stability contract.
- [x] **Confirmed:** `PhaseStabilityModel`, `driving_force`, and the current
  per-species `PhaseSeedPolicy` encode the old `q = 1` workflow. In particular,
  phase activation does not seed the candidate with the composition that
  minimizes its tangent-plane distance.
- [x] **Confirmed:** stability mathematics currently lives in
  `equilibrium_workflows.rs`, while both the prepared `P,T` runner and the
  monolithic `P,H` route consume its limited report and seed semantics.
- [x] **Retain:** active-set projection, bounded outer iterations, hysteresis,
  cycle detection, rollback, and transactional publication remain valid
  orchestration mechanisms. Rewire them to typed TPD evidence rather than
  reimplementing them.
- [x] **Qualified:** a one-component phase may retain an analytical fast path,
  but only as an internal optimization proven equivalent to the general
  `q = 1` TPD problem. It must not remain a separate public physical criterion
  or a hidden fallback.
- [x] **Qualified:** the solver core has `PhaseActivityModel::IdealSolution`,
  while the resolved phase domain currently exposes only `IdealGas` and
  `PureCondensed`. The migration must introduce an explicit ideal-solution
  capability at the phase-domain boundary before the top-level facade claims
  support for a multicomponent condensed solution. This does not imply any
  non-ideal excess-Gibbs model.

### P7.0 Freeze the contract and dependency graph (P0)

- [x] Record every producer and consumer of `compute_phase_stability_reports`,
  `PhaseStabilityReport`, `PhaseStabilityModel`, `ElementPotentialReport`,
  `driving_force`, `dg_create`, `dg_keep`, `seed_activated_phase`,
  `PhaseSeedPolicy`, `PhaseManager`, `PreparedActiveSetCandidate`,
  `PreparedPhaseControlRunner`, monolithic `P,H` active-set restarts,
  complementarity reports, GUI snapshots, presentation/export helpers, and
  reproducibility capsules.
- [x] Freeze the mathematical units and signs: `minimum_tpd`, `dg_create`, and
  `dg_keep` are molar Gibbs quantities in `J/mol`; a sufficiently negative
  minimum requests activation, zero is coexistence within tolerance, and a
  positive minimum is stable against creation of that candidate phase.
- [x] Define phase capability explicitly. `Excluded` means intentionally not
  evaluated and must remain distinct from an unsupported activity model, an
  elementally infeasible candidate, and a failed TPD minimization.
- [x] Add `IdealSolution` to the resolved phase-model contract with validation
  of its allowed physical states and ordered components. Keep
  `PureCondensed` as a one-component convenience constructor or normalize it
  to the same ideal-solution activity implementation; do not infer ideal
  mixing merely from a liquid/solid phase name.
  - [x] **Corrective audit closed:** `PhaseModel::IdealSolution` is now an
    explicit condensed ideal-mixing declaration, while `PureCondensed` rejects
    every component count other than one. `MultiphaseEquilibriumLayout` accepts
    the new semantic model and the bridge maps it explicitly to
    `PhaseActivityModel::IdealSolution`.
- [x] Treat this as an intentional API replacement. Do not add adapters that
  fabricate old `driving_force` or `PureCondensedSpecies` semantics from a new
  TPD result. Migrate all retained consumers in one dependency-directed pass.
  - [x] The audit leaves exactly one mutable owner: `PreparedPhaseControlRunner`.
    The pure `equilibrium_phase_stability` service receives immutable state and
    returns evidence only; workflows classify that evidence and no production
    caller retains obsolete driving-force or per-species seed semantics.
  - [x] `PhaseStabilityStatus` distinguishes policy exclusion,
    `FixedGasAssemblage`, absence of an independent reference assemblage, and
    a genuine evaluated TPD result. Typed errors cover failed construction or
    infeasible minimization rather than being encoded as a numeric status.

### P7.0a Corrective semantic-boundary audit (P0)

The post-refactor architecture review confirms that the constrained TPD kernel
is the right foundation and must not be replaced by the historical pure-phase
criterion. The unfinished work is at the semantic construction boundary and
in lifecycle evidence around it.

- [x] **Accept the diagnostic's central finding.** The current chain is
  inconsistent: `PhaseModel` has only `IdealGas` and `PureCondensed`, while
  `PhaseActivityModel` already has `IdealSolution`; `activity_model_for`
  currently maps `PureCondensed` to that numerical law, and the production
  layout enforces one component for the semantic condensed model. Unit-level
  multicomponent TPD tests therefore do not prove production support.
- [x] Introduce an explicit semantic `PhaseModel::IdealSolution` and validate
  the physical states for which the present ideal-mixing contract is defined.
  Keep `PureCondensed` semantically distinct and one-component, even though its
  numerical activity law is the `q = 1` special case of `IdealSolution`.
- [x] Migrate the complete definition chain without a compatibility adapter:
  `PhaseSpec -> ResolvedPhaseSystem -> MultiphaseEquilibriumLayout ->
  PhaseEquilibriumMetadata -> EquilibriumPhaseDescriptor ->
  PhaseActivityModel`. Update canonical ordering/fingerprints, factories,
  candidate selection, public prelude types, and all exhaustive matches.
- [x] Add a production bridge regression that constructs a multicomponent
  `IdealSolution` through `SubstanceSystemFactory` and the local repository,
  then through `ResolvedPhaseSystem` and `PhaseEquilibriumBuildRequest`, and
  proves that every ordered component,
  thermochemistry capability, element row, phase range, and
  `PhaseActivityModel::IdealSolution` reaches the resulting problem bundle.
  No branch on this path may require `components.len() == 1` for
  `IdealSolution`. The local NASA-condensed fixture proves two thermochemistry
  records, lookup provenance, element rows, descriptors, activity laws, and
  the full component range reach `PhaseEquilibriumProblemBundle` unchanged.
- [ ] Add an end-to-end fixed-`P,T` lifecycle story in which an initially
  inactive multicomponent ideal solution has negative accepted `minimum_tpd`,
  is seeded with the reported `incipient_composition`, is re-solved, and is
  published only after residual, conservation, feasibility, KKT, and
  complementarity validation pass.
- [x] Redesign `SupportedPhaseModelPolicy` after its dependency audit. The
  retained versioned capability contract is now `IdealPhaseModelsV1`; its
  documentation names ideal gas, pure condensed, multicomponent ideal
  solution, and canonical constrained TPD. The obsolete
  `FixedPressureTemperatureV1` spelling has no compatibility alias.
- [x] Extend the GUI phase-model enum and validation after the engine
  semantic bridge is canonical. The GUI must serialize/construct
  `IdealSolution` explicitly rather than disguising it as `PureCondensed`.
- [x] Audit reference-assemblage construction for inactive candidates, active
  retention, a single active phase, multiple active condensed phases, and a
  rank-deficient active elemental space. The audit retained the canonical rule:
  every candidate is compared only against other active phases, a single active
  condensed phase reports `ActiveWithoutReferenceAssemblage`, policy-excluded
  phases are never assigned a numeric TPD, and rank-deficient reference rows
  use the SVD geometry rather than an accidental full-rank assumption. Focused
  orchestration regressions cover exclusion, the missing-reference state,
  multiple active condensed phases, and a rank-one H/O reference space.
- [ ] Add the remaining monolithic `P,H` recovery matrix with typed evidence:
  reduced active-set success without a probe; reduced failure followed by a
  successful bounded all-active probe; a probed but TPD-stable phase remaining
  inactive; a TPD-unstable phase activated from `x*` rather than the neutral
  probe composition; stability evaluated at accepted `T*`; and complete
  rollback when both reduced and recovery solves fail. The recovery probe is
  numerical branch discovery, never an activation criterion.
  - [x] The runner now treats all-active probe occupancy as non-publishable.
    It restores inactive probe phases to trace and permits one deterministic
    transition only when `minimum_tpd < dg_create`, seeded from `x*`; after an
    activation it re-solves the new fixed set before acceptance. Real water
    stories prove that strict `Monolithic` and `Auto` reject a TPD-stable
    liquid probe rather than manufacturing appearance. Deterministic
    reduced-success/recovery/restart/rollback fixtures now protect the
    lifecycle itself; the remaining observability gap is typed probe counters
    in report evidence.
  - [x] **Regression review, retained evidence:** real strict/`Auto` water
    stories prove that a wider numerical probe with `minimum_tpd >= dg_create`
    is rejected rather than published; the dense water/ice `P,H` range proves
    every published TPD report carries the accepted point temperature; generic
    runner rollback preserves accepted continuation; direct multicomponent
    deactivation uses total phase amount; and capsules preserve
    `CanonicalTpdV1` plus named `IdealSolution` semantics.
  - [x] **Deterministic lifecycle harness:** crate-private injected-candidate
    stories now force, independently: (a) a reduced fixed-set success with no
    wider candidate; (b) a numerically wider TPD-stable candidate rejection;
    (c) a TPD-unstable multicomponent `IdealSolution` activation from a
    visibly non-uniform `x*`, followed by a new fixed-set solve; and (d)
    rollback when that restart fails. They assert masks, transition reason,
    `minimum_tpd`, `incipient_composition`, seed ratios, accepted candidate
    temperature, validation evidence, and unchanged accepted continuation.
    The real P,H stories remain responsible for enthalpy residuals; neither
    suite fixes optimizer iteration counts or SVD-call order.
  - [ ] Add compact typed numerical-probe evidence to the P,H lifecycle report
    for the preceding stories: attempted wider masks, outcome class
    (unneeded/failed/TPD-rejected/TPD-activated), and restart outcome. It must
    describe numerical branch discovery without presenting a probe candidate
    as a physical phase transition or adding a public mutable API.
- [x] Strengthen multicomponent deactivation regressions. A phase with one
  trace-small component but substantial `sum_i n_i` must remain active; a
  phase may be deactivated only from its total amount plus the accepted
  retention TPD criterion. The direct regression covers both a trace-small
  component with substantial phase total and a truly vanishing total phase.
- [x] Replace stale production wording such as `inactive pure phase` and
  `Detect phases that must be created (Delta G < 0)` with candidate-phase,
  minimum-TPD, and phase-stability terminology. Retain legitimate historical
  `q = 1` explanations and schema tests that deliberately reject the old
  `driving_force` representation.
- [x] Reassess reproducibility schema compatibility after the serialized
  semantic `PhaseModel` change. `EquilibriumPhaseSpecSnapshot` stores explicit
  model names rather than enum ordinals, so prior `IdealGas` and
  `PureCondensed` capsules retain their exact meaning while current capsules
  can record `IdealSolution`. `CanonicalTpdV1` and capsule schema v2 remain
  correct; the regression round-trips the new symbolic model value. GUI
  documents have the same named-enum property, so their schema remains v1 and
  older documents need no migration.
- [x] **Already satisfied:** canonical state construction, element-potential
  reconstruction, exact elemental feasibility, constrained ideal TPD,
  TPD-derived activation composition, transactional publication, typed
  stability evidence, presentation-only diagnostics, and rejection of stale
  `driving_force` fields are retained. No parallel pure-phase fast path or
  replacement minimizer is required by this audit.
- [x] **Deferred physical decision:** `FixedGasAssemblage` remains explicit
  while any ideal-gas assemblage is active, until the project decides whether
  multiple gas declarations represent one shared mixture, separate
  compartments, or competing phases. The narrower gas-vs-condensed boundary
  is now supported: when no gas assemblage is active, an inactive ideal gas is
  evaluated against the condensed reference so complete condensation and
  reverse evaporation are representable. Generic competing-gas TPD still
  depends on the unresolved activity-normalization contract.
- [x] **Out of scope:** non-ideal activity coefficients, fugacity/EOS phase
  split, and global non-convex multi-start TPD remain in the future-physics
  section. The semantic ideal-solution bridge must not imply these models.

### P7.1 Extract the stability domain and canonical inputs (P0)

- [x] Create `equilibrium_phase_stability.rs`. It owns chemical-potential
  snapshots, element-potential reconstruction, elemental-direction
  feasibility, candidate TPD problems, minimization, and typed stability
  reports. It must not mutate `PhaseSet`, choose transitions, publish a solve,
  or own cycle/hysteresis policy.
- [x] Keep `equilibrium_workflows.rs` and `PreparedPhaseControlRunner`
  responsible only for orchestration: fixed-active-set solve, request a pure
  stability analysis, classify the returned evidence, seed/reproject, detect
  cycles/budgets, and transactionally publish an accepted result.
- [x] Build one immutable accepted-state input containing solver-order moles,
  phase totals, canonical `g_i^0`, canonical activities/chemical potentials,
  conditions, element matrix, active/candidate masks, and layout identity.
  Validate dimensions, finiteness, positivity requirements, and snapshot
  alignment once at construction.
- [x] Reuse `PhaseActivityModel` for `ln(a_i)` and expose only the additional
  composition-level operation needed by TPD. The residual, Jacobian, and
  stability paths must share the same standard-state and pressure convention.
  - [x] First extraction: `equilibrium_phase_stability.rs` now owns an
    immutable canonical-state snapshot. The retained pure-phase workflow
    obtains `mu_i` from that snapshot rather than recomputing its own activity
    convention.

### P7.2 Reconstruct element potentials without a full-rank shortcut (P0)

- [x] Solve `A_active * lambda ~= mu_active` with a documented SVD tolerance.
  Publish `lambda`, numerical rank, singular-value/tolerance evidence, maximum
  absolute residual, scaled residual, and the solver species used in the fit.
- [x] Validate representability of accepted active chemical potentials with
  an absolute-plus-relative residual contract. A poor fit is a typed physical
  validation failure; it must not silently produce phase decisions.
- [x] Do not reject `rank < element_count`. Element potentials may be
  non-unique; stability must remain invariant to null-space-equivalent choices
  of `lambda` for an elementally feasible candidate direction.
- [x] Construct the active elemental range/null space once per unchanged
  active set and cache it in the prepared phase-control state. Rebuild it only
  after an active-set/layout change and report the rebuild/reuse decision.
  - [x] SVD fit now records rank, tolerance, residual, and reference species,
    accepts a rank-deficient valid fit, and the phase-control workflow no
    longer reintroduces a `rank < element_count` veto after reconstruction.
    A null-space-equivalence regression protects the resulting TPD value.
  - [x] `PreparedPhaseControlRunner` now owns a runner-scoped cache of immutable
    reference-assemblage geometry. One cache entry retains the SVD fit factors
    and full `A_active^T A_active` null space; unchanged reference species
    reuse it for elemental-potential fitting, feasibility, and constrained TPD
    constraints. `TemperatureRangeSolveReport` publishes entry/build/reuse
    counters, while a unit regression proves one build plus one reuse leaves
    the fit unchanged. Temperature, activities, chemical potentials, and TPD
    values remain uncached candidate-state data.

### P7.3 Enforce candidate elemental-direction feasibility (P0)

- [x] For candidate composition `x`, evaluate
  `c_beta(x) = A_beta^T * x` and require it to lie in
  `Range(A_active^T)`. Implement the check through an SVD/null-space basis with
  explicit absolute-plus-relative tolerance, not through a full-rank guard.
- [x] Include simplex constraints `x_i >= 0` and `sum(x_i) = 1` in the
  candidate problem itself. Never normalize an unconstrained answer after the
  optimizer and call it a constrained minimum.
- [x] Return a typed infeasible-candidate result when the admissible set is
  empty. Do not evaluate or classify a TPD value for a physically infeasible
  trial composition.
- [x] Report the candidate elemental composition and feasibility residual for
  every evaluated phase so complementarity and release evidence can audit the
  decision.
  - [x] The canonical ideal TPD path constructs elemental null-space rows
    `B*x = 0`, solves them as part of the candidate simplex problem, and then
    independently verifies `A_active^T*y ~= A_candidate^T*x`. No full-rank
    shortcut remains in the connected phase-control path.
  - [x] `CanonicalPhaseState` is the only accepted-state ingress for the
    connected TPD workflow and delegates activity evaluation to the shared
    `PhaseActivityModel`. `ElementalFeasibilityReport` retains candidate
    elemental totals, residual, and tolerance in every evaluated
    `PhaseStabilityReport`; an empty feasible simplex returns a typed error
    before any TPD classification or transition publication.

### P7.4 Implement the general ideal TPD minimum (P0)

- [x] Define the single canonical objective
  `TPD_beta(x) = sum_i x_i * (mu_i^beta(x,T,P) - a_i*lambda)` over the
  admissible simplex. Store the minimum and all thresholds in `J/mol`; no
  arbitrary trace amount belongs to the TPD definition.
- [x] Cover `IdealSolution` candidates with the convex analytical/specialized
  structure: both the full-rank closed-form case and the rank-deficient
  constrained case use deterministic minimization with independently checked
  KKT/feasibility evidence.
- [x] Return both `minimum_tpd` and normalized `incipient_composition = argmin
  TPD` in declared phase-component order. Validate finite objective values,
  finite/nonnegative composition, unit sum, elemental feasibility, and
  deterministic tie handling.
- [x] Handle boundary optima without evaluating `ln(0)`. Zero composition is
  valid in the mathematical minimizer; a positive numerical floor is allowed
  only when constructing log-mole restart coordinates.
- [x] If minimization, feasibility, or validation fails, return a typed error
  with phase/context evidence and abort the transition transaction. Do not
  fall back to uniform composition, a random start, the old pure-phase
  formula, or a penalty-only answer.
- [x] Keep the TPD problem interface ready for a future non-convex activity
  model and multi-start/global policy, but do not implement or imply a fake
  non-ideal model in this stage.
  - [x] Interior pass: `IdealTpdProblem` now supplies an analytical softmax
    minimum for full-rank inputs and a deterministic dual-Newton minimizer for
    rank-deficient interior simplexes. It validates independent elemental
    constraints, feasibility, and KKT residuals; the connected workflow uses
    it for both pure and multicomponent `IdealSolution` phases. Unit tests
    cover `q = 1`, the analytic multicomponent result, rank-deficient interior
    feasibility, and invariance to non-unique lambda representatives.
  - [x] Boundary pass: a deterministic phase-I linear program identifies the
    relative-interior feasible support, then dual Newton runs only on that
    face. Exact zero components remain in `incipient_composition`; the
    positive floor is introduced only by the separate log-mole seed builder.
    Regression tests cover a boundary-only candidate and seed-total recovery.
  - [ ] **Physical decision gate: separately declared `IdealGas` phases.** The
    present lifecycle treats them as one fixed gas assemblage whenever one gas
    phase is active, reported as `FixedGasAssemblage`, rather than as competing
    candidates. An inactive gas is evaluated only when the accepted reference
    contains no active gas, which closes the gas/condensed replacement case
    without deciding multi-gas semantics. Before generic competing-gas TPD,
    decide whether multiple declarations mean one shared mixture, separate
    compartments, or genuinely competing phases. The answer determines
    activity normalization and is therefore not an engineering cleanup.

### P7.5 Replace reports, hysteresis, and activation seeding (P1)

- [x] Replace the old physical-model/driving-force report with a typed result
  containing phase identity, active state, evaluation status,
  `minimum_tpd`, `incipient_composition`, element-potential evidence,
  feasibility evidence, minimizer diagnostics, conditions, and layout
  identity. Excluded/not-evaluated phases carry no fabricated numeric value.
  - [x] `PhaseStabilityReport` now has `PhaseStabilityStatus` plus
    `minimum_tpd`; obsolete `PhaseStabilityModel` and the misleading
    `driving_force` field are removed. Transition and complementarity evidence
    now use the same TPD terminology. Skipped phases report an explicit status
    and no numeric minimum.
  - [x] Evaluated reports retain conditions, canonical ordered layout evidence,
    elemental-feasibility residuals, and constrained-minimizer/KKT diagnostics.
    Skipped reports retain the same identity and conditions but no fabricated
    numerical stability evidence.
- [x] Apply `dg_create` and `dg_keep` only after the mathematical minimum has
  been accepted. Hysteresis is a transition policy over `minimum_tpd`; it must
  not alter the objective or define phase stability.
  - [x] `PhaseManager` classifies only accepted `minimum_tpd` values;
    thresholding is downstream of TPD construction and minimization.
- [x] Replace per-species seed semantics with a total phase-seed amount policy.
  For an activated phase construct `n_i = n_phase_seed*x_i_star`, apply a
  strictly positive floor only for log coordinates, and renormalize so the
  seeded component sum equals the requested total seed.
  - [x] `PhaseSeedPolicy` and the uniform `seed_activated_phase` helper have
    been removed after their last production caller migrated. The monolithic
    all-active recovery uses the same total-seed API with a documented neutral
    numerical probe; physical phase appearance always uses the TPD minimizer.
- [x] Preserve transactional mass handling: failed seed construction,
  projection, reduced solve, candidate validation, or complementarity must
  leave the last accepted composition, temperature, active set, caches, and
  published report unchanged.
  - [x] Phase transition records now retain the chosen `incipient_composition`
    independently of the floored `restart_seed`; both prepared and mutable
    orchestration paths seed an activated phase through
    `seed_activated_phase_with_composition`.
  - [x] A failed prepared-runner lifecycle now restores the earlier accepted
    continuation seed and phase set. The injected-failure regression proves a
    rejected solve cannot consume continuation intended for the next point.
    Projection/formulation/TPD geometry caches are immutable-layout
    memoization, never accepted-solution publication; every numerical use is
    retargeted before solving.
- [x] Update boundary recovery, transition records, cycle fingerprints, and
  complementarity checks to consume `minimum_tpd` and the minimizer-derived
  seed. Record the chosen seed composition in transition evidence.
  - [x] Transition and complementarity paths use `minimum_tpd`; activated
    transitions retain the exact TPD minimizer independently of the floored
    restart seed. Cycle fingerprints deliberately remain phase-set based until
    a composition-sensitive cycle policy has a physical contract.

### P7.6 Integrate the same criterion into P,T and P,H (P1)

- [x] Wire the prepared fixed-`P,T` phase-control runner to the extracted TPD
  service without rebuilding invariant matrices or activity metadata on every
  outer iteration when the active set is unchanged.
- [x] Run monolithic `P,H` stability analysis at the accepted solution
  temperature `T*`, never at the initial temperature seed. Do not create a
  separate enthalpy-specific phase criterion.
- [x] After a `P,H` transition, seed the phase with `x*`, retain accepted `T*`
  as the continuation temperature, rebuild the active-set formulation, and
  solve again against the same target enthalpy transactionally.
- [x] Propagate typed TPD failure and transition evidence through point,
  temperature-range, and enthalpy-range reports. A failed point must not seed
  the next continuation point.
  - [x] The prepared runner is the sole lifecycle owner for fixed `P,T` and
    bounded monolithic `P,H`: it evaluates TPD from the candidate's accepted
    conditions, seeds transitions from `x*`, and publishes the same immutable
    acceptance report through each point/range solution. `PhRangeRequest`
    advances composition and temperature only after a point is accepted; its
    real water/ice release story retains transition and phase-control evidence
    on each accepted point.

### P7.7 Migrate diagnostics and public consumers (P1)

- [x] Update `MultiphaseAcceptanceReport`, phase-control traces, presentation
  rows, GUI route-specific snapshots, tables, plotting boundaries,
  reproducibility capsules, and story-test output to use TPD terminology and
  retain minimizer/feasibility evidence where appropriate.
  - [x] Acceptance summary rows and the GUI phase-evidence table now retain
    per-phase status, accepted `minimum_tpd`, elemental-feasibility residual,
    and KKT residual when a phase was evaluated. Aggregate complementarity is
    displayed separately, so it cannot hide a weak individual minimization.
- [x] Update stable facade/prelude exports to expose the new immutable report
  types while keeping minimizer implementation details crate-private.
  - [x] `ChemEquilibrium::prelude` now exposes the immutable phase-control,
    transition, TPD, feasibility, minimizer, validation, and typed-index
    evidence. The canonical-state and minimizer implementation remains
    crate-private; external consumers can inspect accepted reports but cannot
    invoke a second mutable orchestration path.
- [x] Add a schema/version decision for serialized reproducibility evidence;
  reject stale old stability fields explicitly rather than interpreting them
  under new semantics.
  - [x] Reproducibility capsule schema v2 records `CanonicalTpdV1` explicitly.
    Its typed JSON loader rejects every older schema and recursively rejects
    obsolete `driving_force`/`PhaseStabilityModel` fields rather than silently
    reinterpreting them as TPD evidence.
- [x] Remove console-only diagnostics from minimization. Timings, iterations,
  termination, KKT/feasibility residuals, and cache reuse belong in typed
  reports and deterministic story tables.
  - [x] The canonical TPD/minimizer, prepared phase-control runner, and
    fixed-`P,T` range route contain no direct `print!`/`println!` calls.
    Runtime evidence is carried by immutable reports and story-table formatters;
    legacy backend debug logging is outside this minimization contract.

### P7.8 Required test matrix (P0/P1)

- [x] Pure `q = 1` regression: general TPD equals
  `g_beta^0(T) - a_beta*lambda`, including sign, units, and minimizer `[1]`.
- [x] Stable and unstable multicomponent ideal-solution fixtures: verify the
  sign and value of `minimum_tpd`, normalized `x*`, and transition decision.
- [x] Full-rank analytical ideal minimum: compare both `minimum_tpd` and `x*`
  against the closed-form reference, not only the final active set.
- [x] Seed regression: component seed ratios match `x*`, every log-mole seed is
  positive, and physical seeded moles sum to the requested phase amount after
  floor-and-renormalize handling.
- [x] Boundary optimum: permit exact zero components in `x*`, avoid `ln(0)`,
  and prove that the log-coordinate floor does not change the reported TPD
  minimizer.
- [x] Rank-deficient active assemblage: an elementally feasible candidate is
  accepted, and TPD is invariant under admissible null-space changes to
  `lambda` within tolerance.
- [x] Elementally infeasible candidate: return explicit infeasibility with no
  TPD classification, activation, or published partial state.
- [x] Hysteresis story: activation below `dg_create` is retained when the next
  accepted minimum lies between `dg_create` and `dg_keep`; test the reverse
  direction separately.
- [x] P,H regression: stability uses accepted `T*`; activation retains `T*` as
  the next continuation seed and preserves the target enthalpy contract.
- [x] Property/metamorphic tests: normalization, non-negativity, finite TPD,
  feasibility residual, element-potential residual, permutation invariance,
  inventory scaling, mass conservation, rollback, cycle/budget termination,
  and cache reuse across repeated active sets. The mathematical unit matrix
  includes direct component-permutation invariance; live reactive-gas stories
  cover inventory scaling, while phase-control regressions cover rollback,
  cycle/budget termination, and runner-scoped geometry-cache reuse.
- [ ] Re-run all retained `P,T`, `P,H`, fixed-phase, phase-transition, GUI, and
  release story matrices. Add offline real-data water/ice, water/liquid, and
  carbon/graphite evidence without mutating JSON libraries; synthetic tests
  remain as precise mathematical fixtures.
  - [x] The current release ledger already records the real TPD
    water/ice, water/liquid, hot-water, graphite, and hot-carbon transition
    matrix in `STORY_TESTS.md`, including conservation, complementarity,
    transition count, timing, and JSON immutability. This is retained evidence,
    not a substitute for the final whole-suite release campaign.
  - [x] The paired release ice T-range story records `tpd_geometry=1/1/3`:
    one runner-scoped active-element geometry build and three reuses while the
    real range crosses a phase transition. It simultaneously verifies the
    projection/prepared/RST cache layout, continuation, conservation, and JSON
    immutability.

### P7.9 Cleanup and definition of done (P1)

- [x] Delete obsolete `PhaseStabilityModel` variants/enum, old
  `driving_force` data paths, the multicomponent capability veto,
  full-element-rank rejection, uniform per-species activation seeding, and
  tests asserting the old limitation. The remaining `ValidationNotApplicable`
  uses are generic typed applicability errors, not a solution-model veto.
  `PureCondensed` retains its intentional one-component semantic constraint;
  the explicit `IdealSolution` model carries multicomponent mixtures. The
  legacy indexed numerical API now labels an undifferentiated ideal-solution
  activity law as `IdealSolution`, never falsely as `PureCondensed`.
- [ ] Repair the stale/corrupted phase-stability subsection in
  `ARCHITECTURE_RU.md`. It still describes scalar `Delta G` driving force and
  says multicomponent solutions are unsupported; rewrite it around constrained
  `minimum_tpd`, `incipient_composition`, and `FixedGasAssemblage` without
  altering the executable stability kernel.
- [x] Remove duplicated pure-phase stability helpers and any internal
  compatibility translation. An optional `q = 1` fast path must be tested
  against the general formulation and remain invisible to consumers.
- [ ] Mark P7 complete only when the canonical dataflow is:
  `accepted fixed-set state -> canonical mu -> element potentials -> feasible
  TPD minimization -> minimum_tpd + x* -> hysteresis + seed -> rebuilt
  active-set solve -> transactional acceptance` for both `P,T` and `P,H`.
  The kernel and current pure-condensed production path satisfy this flow, but
  P7 remains open until a multicomponent semantic `IdealSolution` traverses
  the same production facade and the bounded `P,H` recovery matrix proves that
  a neutral all-active probe cannot leak into physical activation semantics.
- [x] Document the remaining boundary honestly: P7 completes general phase
  stability for the currently supported ideal activity models. Non-ideal
  excess-Gibbs/fugacity models and global multi-start policy remain in F2.

## Deferred Foundation: Local NIST Fixture and New Physics

The current production path is deliberately scoped to ideal, fixed-pressure
equilibrium with local data. The items in this section are not small cleanup
tasks. They require an explicit data release or a new physical contract, and
must not be quietly approximated by the existing NASA-gas fixtures.

### F1. Local NIST fixture is a data-release blocker

- [ ] **Create a stable, read-only local NIST thermochemistry fixture.**
  The NIST parser is an online acquisition boundary, while production and
  regression equilibrium calculations must remain offline and deterministic.
  The current local index advertises historical NIST addresses but has no
  matching payload records; that is not an acceptable cross-format fixture.

  Required deliverables:
  - choose a small physically meaningful record set with pinned canonical
    names, library keys, physical states, temperature intervals, elemental
    composition, and NIST/Shomate payloads;
  - write the payload, address/index, and elemental-composition entries as one
    reviewed atomic data update, plus a manifest/fingerprint that identifies
    the data release used by the tests;
  - make the repository consistency report show every selected NIST index
    entry paired with an exact payload, without orphan aliases;
  - add offline `P,T` point/range and `P,H` inverse/range stories that retain
    NIST provenance, verify conservation and acceptance, and prove that no
    network call or JSON mutation occurred;
  - add at least one mixed local NASA/NIST story only after the individual
    NIST records are stable, with explicit per-component provenance and a
    documented temperature-domain intersection.

  Do not use live NIST access, a parser mock, or an incomplete historical
  index as a substitute. This fixture proves the format-agnostic contract; it
  is not a second solver implementation.

### F2. Future physical extensions: explicitly outside the current model

- [ ] **Non-ideal phase thermodynamics.** Introduce activity/fugacity models
  with typed parameters, temperature/pressure validity ranges, and analytic or
  independently validated derivative contracts. Ideal-gas and pure-condensed
  activity semantics must remain selectable baselines, not implicit fallbacks.
- [ ] **Liquid/solid solution models.** Support real solution excess Gibbs
  models, composition-dependent chemical potentials, and phase-qualified
  parameter provenance. Reject unsupported solution models before a nonlinear
  backend starts; never silently treat them as ideal.
- [ ] **Phase coexistence and latent heat.** Define the physical contract for
  liquid-solid, liquid-vapor, and solid-solid coexistence, phase fractions,
  and latent contributions in `P,H` solves. Existing phase-control appearance
  evidence is not a general coexistence model.
- [ ] **Non-ideal/global TPD extension.** After P7 establishes the general TPD
  contract for the current ideal activity models, add model-specific
  non-convex minimization, deterministic multi-start/global policy, and
  independently validated chemical-potential derivatives. Do not weaken P7
  or reuse an ideal closed-form minimum for a non-ideal phase.
- [ ] **Additional constraints and state variables.** Design separate typed
  workflows for `P,V`, `U,V`, pressure sweeps, and reactive-flash problems.
  They must not be encoded as ad-hoc switches inside the fixed-pressure
  `P,T`/`P,H` request types.
- [ ] **Electrolyte/charged-species support.** Add electroneutrality,
  reference-state conventions, ionic-strength/activity models, and a clear
  contract for aqueous phases only when suitable local data and validation
  fixtures exist.

### F3. Quality-of-life and observability backlog

These items are safe to implement independently of new physics. They must
consume immutable solution/report snapshots and must not reopen `SubsData` or
change accepted numerical results.

- [x] **Human-readable report views.** `equilibrium_presentation` now projects
  an accepted immutable solution into deterministic phase, component/provenance,
  backend-attempt, and timing rows for GUI, CLI, or export consumers. Its
  compact ASCII renderer is a convenience view only; structured rows remain the
  canonical presentation contract and retain failure evidence.
- [x] **Range presentation and plotting series.** `TemperatureRangePresentationReport`
  now exposes stable ordered `P,T` component-mole and phase-total series,
  point-level continuation/timing/validation rows, and explicit phase-transition
  boundaries. `PhRangePresentationReport` now mirrors this with a correctly
  typed target-enthalpy axis plus solved-temperature/component/phase raw series
  and route/fallback/continuation evidence. Both range views also publish raw
  residual/balance/per-point-time columns and phase-local component mole
  fractions. `EquilibriumDisplayPolicy::visible_amount_columns()` returns
  display-column indices without rebuilding or dropping the canonical raw
  series, so plotting/UI code can apply trace filtering safely.
- [x] **Phase-aware postprocessing.** PCHIP remains an optional presentation
  layer, but typed `P,T` range resampling is now forbidden when accepted
  phase-control transitions occurred; callers keep raw points or explicitly
  split phase-stable ranges. Raw and resampled rows remain distinct, and the
  existing policy still makes the linear/log interpolation space explicit.
  Future range presentation must extend the same guard to failed gaps and any
  future layout-changing workflow.
- [x] **Run comparison report.** `EquilibriumComparisonReport` compares two
  accepted single-point solutions only after exact phase-qualified layout
  validation. It reports component/phase mole deltas, conditions, backend,
  balance/residual evidence, and lookup provenance changes. A future range
  comparison may aggregate these immutable point comparisons without altering
  the strict layout contract.
- [x] **Reproducibility capsule.** `equilibrium_reproducibility` now exports
  a JSON-compatible immutable snapshot of solved conditions, phase specs,
  selected per-component record identities/provenance, effective backend
  cascade and numerical options, acceptance evidence, optional candidate
  selection, catalog consistency, and a caller-supplied data-release label.
  The record fingerprint is deliberately an identity fingerprint and the
  catalog value is structural evidence, not a false claim of a payload-content
  hash; **F1** remains responsible for a versioned local-data manifest. This
  is provenance/report export, not mutable project import/export.
- [x] **Failure-focused diagnostics.** `backend_attempt_rows_from_error()`
  exposes the complete ordered trace retained by `AllBackendsFailed` and
  `CascadeAborted`, including termination, iterations, callback timing, and
  typed cause. Range-point/seed context remains owned by the range reports and
  can be paired with these rows without fabricating a successful solution.
- [x] **Display thresholds and units policy.** `equilibrium_display` now
  provides validated trace-row filtering, scientific/fixed/engineering-SI
  number styles, and explicit fraction-versus-percent rendering. It only
  borrows or projects immutable presentation rows: thresholding never alters
  solver inputs, conservation checks, or raw data retained for export.
- [x] **Opt-in structured phase-lifecycle diagnostics.**
  `equilibrium_diagnostics` now records bounded typed events for the canonical
  `PreparedPhaseControlRunner`: solve/outer-iteration boundaries, accepted
  fixed-set candidates, constrained TPD evidence, hysteresis holds, accepted
  activation/deactivation, boundary recovery, and final publication. The
  default is disabled; enabled reports are attached immutably to the accepted
  solution and an optional live sink can observe the same events. The separate
  `equilibrium_diagnostics_display` adapter resolves phase/component labels and
  can emit a human-readable tree through the existing `log` facade without
  putting formatting or logging side effects in the physical solver.
- [x] **Complete engine diagnostics coverage for rejected paths and long
  ranges.** Typed rejected-candidate, boundary-probe, rollback, outer-budget,
  cycle, and `P,H` `Auto` route-fallback events stream through the live sink.
  Bounded P,T and P,H ranges support endpoints, every-N, every-point, and
  `TransitionsOnly` policies. The last policy defers the sink until an
  accepted point actually changes the phase set, then retains/replays that
  point trace only. Accepted `P,H Auto` results retain an immutable
  `PhRouteDecision::AutoFallback` independently of scalar trials.
- [x] **Expose diagnostics in the equilibrium GUI.** The document now maps an
  explicit lifecycle detail level and range retention policy to the canonical
  engine options. Accepted point snapshots project the immutable typed report
  into a phase-qualified collapsible decision tree; P,H also exposes its
  route decisions independently from scalar trials. While a worker runs, the
  same opt-in engine sink streams bounded transient events through the
  ticket/fingerprint gate. GUI formatting remains a consumer of engine
  evidence, never a second source of phase-control decisions.

## Recommended implementation passes

1. [x] Characterize current fixtures and introduce typed problem/result/error
   boundaries without changing numerical behavior.
2. [x] Extract the pure log-moles formulation and common acceptance gate.
3. [x] Characterize the RST backend matrix, then complete attempt counters,
   termination mapping, and a measured production default.
4. [x] Isolate legacy solvers and keep them behind an explicit fallback
   policy with typed diagnostics.
5. [x] Build the independent K_eq validator and cross-validation fixtures.
6. [x] Bridge `ResolvedPhaseSystem`, migrate retained classical workflows, and
   remove duplicate orchestration from the production surface.
7. [ ] Build the typed GUI on the stable facade and reports.

## Historical definition of done

The original checklist below is retained for traceability and synchronized
with the current fixed-`P,T` production contract. Legacy orchestration removal
is intentionally deferred until compatibility consumers have migrated; the
handwritten LM/NR/TR implementations remain supported numerical fallbacks.

- [x] There is one canonical equilibrium problem and solution model.
- [x] Every accepted solution has passed backend-independent numerical and
  physical validation.
- [x] Fallback order, retry rules, budgets, and attempts are explicit and fully
  reported.
- [x] RustedSciThe provides the production nonlinear backends.
- [x] The K_eq path independently validates every applicable small-system
  fixture and clearly reports when it is not applicable.
- [x] Tests cover formulation, every backend, cascade behavior, physical
  invariants, cross-validation, transactionality, and offline integration.
- [x] Legacy mutable orchestration is absent from the production facade. Its
  deprecated compatibility entry points live under `ChemEquilibrium::legacy`,
  while legacy numerical fallback backends remain intentionally supported.

---

## SourceCraft Diagnostics

Результаты ревизии кода модуля `ChemEquilibrium`, проведённой 23.07.2026.

### Проверка диагностики (26.07.2026)

Каждый пункт ниже имеет итоговый бинарный статус: `✅ выполнен` либо
`❌ отклонён`. Формулировка «отклонён как локальное исправление» означает, что
сама архитектурная тема признана реальной, но предложенный SourceCraft патч
небезопасен или преждевременен; связанная миграция остаётся в основном плане.

- [x] **A.1 подтверждён и исправлен.** `R` теперь использует полное значение
  CODATA 2018 `8.314_462_618_153_24 J/(mol K)`; есть прямой regression-test.
- [x] **A.2 подтверждён и исправлен.** Автоматический RST путь требует либо
  полный `gibbs_sym` snapshot в solver order, либо успешный thermochemistry
  lookup для каждого requested substance. Частично заполненные кэши больше не
  считаются символическим контекстом.
- [x] **A.3 подтверждён и исправлен.** Мутирующий compatibility method переименован
  в crate-private `publish_reconstructed_moles`; чистое преобразование остаётся
  свободной функцией `compute_species_moles`.
- [x] **B.2 подтверждён и исправлен.** Неиспользуемый
  `VariableScalingContract` и тесты его изолированной арифметики удалены. Если
  variable scaling понадобится, его следует вводить только вместе с реальной
  передачей масштаба во все backend contracts.
- [x] **E.1/E.2 подтверждены и исправлены.** Численные допуски получили имена,
  а подробные дампы initial guess, inventory и stoichiometric matrix понижены
  с `info!` до `debug!`.
- [x] **A.4/C.1/C.2/C.3/C.4 отклонены как локальные SourceCraft-патчи.** Полная
  миграция mutable legacy facade к `EquilibriumProblem`/`PreparedEquilibriumProblem`
  и устранение ручной active-set reconstruction остаются отдельным этапом, а
  не безопасным локальным патчем.
- [x] **D.2 подтверждён и исправлен.** Добавлен serial-only
  `solve_for_T_range_with_phase_control`: он разделяет canonical обновление
  Gibbs/publish sweep с ordinary serial path, переносит accepted seed и phase
  mask между точками и покрыт offline NASA-gas regression test. Parallel
  phase-control sweep остаётся намеренно неподдержанным, поскольку ему нужна
  отдельная independent semantics.
- [x] **B.1, B.3, B.4, B.5, E.4-E.6 отклонены как дефекты.** Это осознанные
  контракты или policy choices: глобальный budget независим от per-attempt
  лимита, `Required` честно сообщает о неприменимости K_eq validation,
  explicit hysteresis имеет размерную семантику, а trace floors/default scaling
  нельзя менять без отдельных stability benchmarks.

### A. Потенциальные ошибки (bugs)

#### A.1 `R = 8.314` — неверное значение газовой постоянной

**Где:** [`equilibrium_log_moles.rs:100`](equilibrium_log_moles.rs:100)

**Проблема:** `pub const R: f64 = 8.314;` — это значение не является стандартной молярной газовой постоянной. Правильное значение: `8.314462618`. Ошибка в 0.00046 J/(mol·K) даёт систематическое смещение в логарифме константы равновесия: `ln(K) = -ΔG/(RT)`. При T=1000K, ΔG=-100 kJ/mol: `ln(K)_true = 100000/(8.31446*1000) = 12.027`, `ln(K)_current = 100000/(8.314*1000) = 12.028`. Ошибка мала (~0.01%), но накапливается в температурных сериях.

**Важность:** Средняя. Для инженерных расчётов ошибка незначительна, но для научных публикаций — недопустима.

**Статус (26.07.2026):** ✅ **Выполнен.** Константа заменена на полное
значение CODATA 2018; прямой regression-test защищает значение.

#### A.2 `has_rst_symbolic_context()` — некорректная эвристика

**Где:** [`equilibrium_log_moles.rs:490-496`](equilibrium_log_moles.rs:490)

**Проблема:** Метод проверяет наличие символьного контекста через `gibbs_sym.len() == substances.len() && !gibbs_sym.is_empty()`, но также возвращает `true`, если `search_results`, `search_states`, `therm_map_of_sym` или `therm_map_of_fun` непусты. Последние четыре условия не гарантируют, что символьные Gibbs-функции действительно построены для всех веществ. Это может привести к тому, что `SolverPolicy::rusted_scithe_default()` будет выбран, но RST-адаптер не сможет построить SymbolicNonlinearProblem.

**Важность:** Средняя. Может проявляться как трудноотлавливаемая ошибка "RST backend failed" при определённых последовательностях вызовов.

**Статус (26.07.2026):** ✅ **Выполнен.** Эвристика требует либо полного
`gibbs_sym` snapshot, либо успешного thermo lookup для каждого вещества;
частичные caches больше не включают RST автоматически.

#### A.3 `compute_species_moles` — публичная функция с побочным эффектом

**Где:** [`equilibrium_log_moles.rs:2173-2178`](equilibrium_log_moles.rs:2173)

**Проблема:** Публичный метод `compute_species_moles(&mut self, sol: Vec<f64>)` мутирует `self.moles` и `self.map_of_moles_for_each_substance`. Название предполагает чистую функцию (как у свободной функции `compute_species_moles(log_moles: &[f64])`), но на самом деле это метод с побочным эффектом. Это нарушает принцип наименьшего удивления.

**Важность:** Низкая. Косметическая проблема именования.

**Статус (26.07.2026):** ✅ **Выполнен.** Мутирующий метод стал
crate-private `publish_reconstructed_moles`; чистая свободная функция сохранила
имя `compute_species_moles`.

#### A.4 `EquilibriumLogMoles` — 22 публичных поля

**Где:** [`equilibrium_log_moles.rs:275-320`](equilibrium_log_moles.rs:275)

**Проблема:** Структура имеет 22 публичных поля, которые можно изменять извне в любом порядке. Нет гарантии, что после изменения одного поля (например, `T`) остальные поля (например, `gibbs`) остаются консистентными. Это прямой путь к багам "забыл обновить gibbs после смены T".

**Важность:** Высокая. Основной источник производственных ошибок.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.** Это реальный
долг legacy facade, но безопасный путь — дальнейшая миграция callers на
`EquilibriumProblem`/`PreparedEquilibriumProblem`, а не приватизация полей
одним разрушающим патчем.

### B. Overengineering

#### B.1 `SolverCascadeBudget` — избыточная сложность для текущих потребностей

**Где:** [`equilibrium_solver_policy.rs:38-46`](equilibrium_solver_policy.rs:38)

**Проблема:** Бюджет каскада содержит три независимых лимита (`max_attempts`, `max_iterations_per_attempt`, `max_total_iterations`), которые проверяются в [`solve_backend_cascade`](equilibrium_log_moles.rs:2064-2073). При этом `max_total_iterations` вычисляется как `max_iter * backends.len()` — то есть всегда пропорционален двум другим. Три лимита вместо одного — это overengineering. Достаточно двух: `max_attempts` и `max_iterations_per_attempt`.

**Важность:** Средняя. Усложняет понимание без реальной выгоды.

**Статус (26.07.2026):** ❌ **Отклонён.** `max_total_iterations` не обязан
равняться произведению двух остальных лимитов: он задаёт независимый глобальный
resource cap и уже проверяется regression-тестами.

#### B.2 `VariableScalingContract` — объявлен, но нигде не используется

**Где:** [`equilibrium_problem.rs:313-384`](equilibrium_problem.rs:313)

**Проблема:** `VariableScalingContract` полностью реализован (с `apply_iterate`, `unscale_iterate`, валидацией), но не используется ни в одном бэкенде. Row scaling (`ResidualScalingContract`) используется, а variable scaling — нет. Это мёртвый код.

**Важность:** Средняя. Увеличивает когнитивную нагрузку при чтении.

**Статус (26.07.2026):** ✅ **Выполнен.** Неиспользуемый
`VariableScalingContract` удалён вместе с изолированными тестами; при реальной
потребности его следует вводить сразу через backend contracts.

#### B.3 `EquilibriumConstantSolverMode::Required` — нереализуемый контракт

**Где:** [`equilibrium_constant_solver.rs:74-75`](equilibrium_constant_solver.rs:74)

**Проблема:** Режим `Required` требует, чтобы K_eq валидация была выполнена для каждого решения. Но K_eq решатель работает только для однореакционных систем. Для многореакционных систем `Required` гарантированно упадёт с ошибкой. Этот режим нельзя использовать в production, только в тестах.

**Важность:** Низкая. Но вводит в заблуждение.

**Статус (26.07.2026):** ❌ **Отклонён.** `Required` — намеренный строгий
контракт для workflows, где независимая validation обязательна. Неприменимость
к многореакционной задаче возвращается как typed error, а не маскируется.

#### B.4 `TemperaturePostprocessingPolicy` — избыточная гибкость

**Где:** [`equilibrium_temperature_postprocessing.rs`](postprocessing_and_logging/equilibrium_temperature_postprocessing.rs)

**Проблема:** Политика постобработки поддерживает три типа сеток (RawOnly, Uniform, Explicit) и два типа интерполяции (Linear, Log). При этом в коде нет ни одного вызова с Explicit или Log. Вся гибкость существует только для тестов.

**Важность:** Низкая. YAGNI-нарушение.

**Статус (26.07.2026):** ❌ **Отклонён.** Это публичная policy граница для
postprocessing и GUI/CLI consumers; отсутствие текущих production callers не
делает корректные modes мёртвым кодом.

#### B.5 `PhaseHysteresisPolicy::Explicit` — дублирование TemperatureScaled

**Где:** [`equilibrium_workflows.rs:1017-1025`](equilibrium_workflows.rs:1017)

**Проблема:** `Explicit { dg_create, dg_keep }` — это фактически `TemperatureScaled` с замороженной температурой. При T=const они эквивалентны. Можно было бы обойтись одним вариантом.

**Важность:** Низкая.

**Статус (26.07.2026):** ❌ **Отклонён.** `TemperatureScaled` хранит
безразмерные множители `R*T`, а `Explicit` — физические пороги `ΔG` в J/mol.
Они намеренно различаются при температурных сериях.

### C. Недостатки дизайна

#### C.1 `EquilibriumLogMoles` — God Object

**Где:** [`equilibrium_log_moles.rs:275-320`](equilibrium_log_moles.rs:275)

**Проблема:** Структура содержит 22 поля, объединяющие термохимические данные, настройки решателя, состояние температурной серии, опубликованные решения, отчёты валидации и фазы. Это классический God Object. Новый API (`EquilibriumProblem` + `EquilibriumSolution`) уже решает эту проблему, но старый фасад остаётся в production.

**Важность:** Высокая. Затрудняет тестирование и поддержку.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.** Новый typed
pipeline уже служит source of truth для новых bridge workflows; legacy facade
будет сужаться по мере миграции, а не заменяться внезапно.

#### C.2 Две параллельные иерархии ошибок

**Где:** [`equilibrium_nonlinear.rs:55-72`](equilibrium_nonlinear.rs:55) и [`equilibrium_nonlinear.rs:89-188`](equilibrium_nonlinear.rs:89)

**Проблема:** `SolveError` (для численных решателей) и `ReactionExtentError` (для всего остального) — две пересекающиеся иерархии. `SolveError` мог бы быть вариантом `ReactionExtentError`, но они разделены. Это приводит к тому, что в некоторых местах ошибка оборачивается (`ReactionExtentError::SolveError`), а в некоторых — нет.

**Важность:** Средняя.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.**
`SolveError` остаётся внутренним численным уровнем, а `ReactionExtentError`
сохраняет доменный контекст и cascade trace. Слияние требует отдельного API
решения, не механического переименования.

#### C.3 `from_problem` — деструктуризация типобезопасности

**Где:** [`equilibrium_log_moles.rs:615-644`](equilibrium_log_moles.rs:615)

**Проблема:** Метод `from_problem` принимает типобезопасный `EquilibriumProblem`, но раскладывает его в 22 публичных поля `EquilibriumLogMoles`. Вся типобезопасность теряется. Это мост для миграции, но он должен быть временным.

**Важность:** Средняя.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.**
`from_problem` —
контролируемый migration bridge; его устранение возможно только после сужения
legacy facade и переноса оставшихся callers.

#### C.4 `solve_fixed_active_set_candidate` — клонирование всего solver'а

**Где:** [`equilibrium_workflows.rs:1439-1476`](equilibrium_workflows.rs:1439)

**Проблема:** Для решения с проекцией активного набора создаётся полная копия `EquilibriumLogMoles` (`let mut local = EquilibriumLogMoles::empty()`), в которую копируются поля по одному. Это 30+ строк ручного копирования. Любое новое поле в `EquilibriumLogMoles` нужно не забыть добавить и сюда.

**Важность:** Средняя. Хрупкий код.

**Статус (26.07.2026):** ❌ **Отклонён как локальное исправление.** Ручная
реконструкция active-set solver изолирована в одном private workflow и покрыта
projection regressions, но должна исчезнуть при окончательном переносе этого
пути на immutable prepared-problem snapshots.

### D. Пробелы в тестовом покрытии — статус

> **Статус:** 70 из 70 пунктов закрыты. D.2 закрыт отдельным serial-only
> `solve_for_T_range_with_phase_control` regression path.

#### D.1 ✅ `EquilibriumLogMoles::solve()` с Legacy бэкендами

**Статус:** 5 тестов в `equilibrium_log_moles_tests.rs` (`legacy_lm_solve_publishes_accepted_solution`, `legacy_nr_solve_publishes_accepted_solution`, `legacy_tr_solve_publishes_accepted_solution`, `legacy_solve_fails_gracefully_without_stoich_matrix`, `legacy_solve_publishes_solve_report_with_attempts`).

#### D.2 ✅ `solve_for_T_range` с phase control

**Статус:** `solve_for_T_range_with_phase_control()` выполняет bounded
`solve_with_phase_control()` для каждой точки последовательной температурной
серии, сохраняет accepted continuation seed/phase mask и публикует только
полностью согласованную sweep table. Offline NASA-gas regression:
`phase_controlled_temperature_sweep_publishes_each_accepted_point`.

#### D.3 ✅ `EquilibriumLogMoles::from_problem()`

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`from_problem_preserves_all_problem_data`, `from_problem_rejects_invalid_problem`).

#### D.4 ✅ `accepted_solution()` после `solve_with_phase_control()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`accepted_solution_after_phase_control_returns_ok_for_single_phase`, `accepted_solution_after_phase_control_preserves_element_totals`, `accepted_solution_after_phase_control_fails_without_solve`).

#### D.5 ✅ `current_problem_snapshot()`

**Статус:** Косвенно через `from_problem` тесты (D.3). Прямой вызов невозможен — метод приватный.

#### D.6 ✅ `run_keq_cross_validation()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`solve_with_keq_validation_enabled_does_not_crash`).

#### D.7 ✅ `EquilibriumSolverSettings::validate()`

**Статус:** 11 тестов в `equilibrium_log_moles_tests.rs` (`solver_settings_validate_accepts_default_settings`, `solver_settings_validate_rejects_zero_max_iter`, `solver_settings_validate_rejects_non_finite_tol`, `solver_settings_validate_rejects_zero_tol`, `solver_settings_validate_rejects_negative_lambda`, `solver_settings_validate_rejects_non_finite_alpha_min`, `solver_settings_validate_rejects_delta_max_less_than_delta_init`, `solver_settings_validate_rejects_eta_out_of_range`, `solver_settings_validate_accepts_single_backend_policy`, `solver_settings_validate_rejects_zero_budget_limits`, `solver_settings_validate_rejects_non_finite_keq_tolerances`).

#### D.8 ✅ `compute_phase_totals()`

**Статус:** 3 теста в `equilibrium_log_moles_tests.rs` (`compute_phase_totals_sums_moles_by_phase`, `compute_phase_totals_handles_empty_input`, `compute_phase_totals_handles_single_phase`).

#### D.9 ✅ `initial_phase_activity()`

**Статус:** 5 тестов в `equilibrium_log_moles_tests.rs` (`initial_phase_activity_marks_active_phases`, `initial_phase_activity_rejects_dimension_mismatch`, `initial_phase_activity_rejects_invalid_phase_eps`, `initial_phase_activity_rejects_out_of_bounds_phase`, `initial_phase_activity_rejects_non_finite_moles`).

#### D.10 ✅ `build_multiphase_acceptance_report()`

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`build_multiphase_acceptance_report_rejects_dimension_mismatch`, `build_multiphase_acceptance_report_accepts_matching_dimensions`).

#### D.11 ✅ `PhaseManager::detect_phase_destruction()`

**Статус:** 4 теста в `equilibrium_log_moles_tests2.rs` (`detect_phase_destruction_returns_empty_when_no_phase_below_threshold`, `detect_phase_destruction_ignores_inactive_phases_below_threshold`, `detect_phase_destruction_finds_active_phase_below_threshold`, `detect_phase_destruction_handles_empty_input`).

#### D.12 ✅ `PhaseManager::classify_phases()` (без температуры)

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`classify_phases_returns_no_transition_when_all_phases_stable`, `classify_phases_rejects_temperature_scaled_hysteresis`, `classify_phases_rejects_dimension_mismatch`).

#### D.13 ✅ `reject_repeated_phase_set()` с пустым множеством

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`reject_repeated_phase_set_accepts_first_occurrence`, `reject_repeated_phase_set_rejects_duplicate`).

#### D.14 ✅ `validate_phase_set_candidate()` с корректным набором

**Статус:** 4 теста в `equilibrium_log_moles_tests.rs` (`validate_phase_set_candidate_accepts_valid_set`, `validate_phase_set_candidate_rejects_dimension_mismatch`, `validate_phase_set_candidate_rejects_negative_total`, `validate_phase_set_candidate_rejects_inactive_phase_with_moles`).

#### D.15 ✅ `seed_activated_phase()` с разными `PhaseSeedPolicy`

**Статус:** 7 тестов в `equilibrium_log_moles_tests2.rs` (`seed_activated_phase_trace_floor_sets_log_moles_to_floor`, `seed_activated_phase_absolute_per_species_sets_positive_moles`, `seed_activated_phase_relative_to_system_total_scales_with_inventory`, `seed_activated_phase_rejects_phase_with_no_species`, `seed_activated_phase_rejects_dimension_mismatch`, `seed_activated_phase_absolute_rejects_non_positive_moles`, `seed_activated_phase_relative_rejects_invalid_fraction`).

#### D.16 ✅ `deactivate_phases_seed_only()` с несколькими фазами

**Статус:** 5 тестов в `equilibrium_log_moles_tests2.rs` (`deactivate_phases_seed_only_sets_species_to_trace_floor`, `deactivate_phases_seed_only_handles_multiple_phases`, `deactivate_phases_seed_only_rejects_non_positive_trace_floor`, `deactivate_phases_seed_only_rejects_non_finite_trace_floor`, `deactivate_phases_seed_only_handles_empty_deactivation_list`).

#### D.17 ✅ `PhaseEquilibriumProblemBundle::solve_with()`

**Статус:** Тесты существуют в `phase_equilibrium_problem_tests.rs` (например, `local_phase_data_solves_through_one_accepted_bridge_bundle`).

#### D.18 ✅ `PhaseEquilibriumSolutionBundle::into_multiphase_solution()`

**Статус:** Тесты существуют в `phase_equilibrium_problem_tests.rs` и `equilibrium_multiphase_story_tests.rs`.

#### D.19 ✅ `MultiphaseEquilibriumSolution::moles_for()` и `mole_fraction_for()`

**Статус:** Тесты существуют в `equilibrium_multiphase_story_tests.rs` (`accepted_fixed_phase_solution_exposes_qualified_amounts_totals_and_summary`).

#### D.20 ✅ `MultiphaseEquilibriumSolution::aggregate_moles_by_substance()`

**Статус:** Тесты существуют в `equilibrium_multiphase_story_tests.rs`.

#### D.21 ✅ `TemperatureSolveSnapshot` и `TemperatureSolveFailure`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`temperature_solve_snapshot_holds_expected_fields`, `temperature_solve_failure_holds_error_message`, `temperature_solve_failure_carries_backend_attempts`).

#### D.22 ✅ `TemperatureWorkerSeed::apply()`

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`temperature_worker_seed_apply_copies_state_to_local_solver`, `temperature_worker_seed_apply_clears_previous_solution`).

#### D.23 ✅ `continuation_seed_for_point()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`continuation_seed_for_point_uses_previous_accepted_when_policy_says_so`, `continuation_seed_for_point_uses_configured_seed_when_independent`, `continuation_seed_for_point_returns_none_when_no_seed_available`).

#### D.24 ✅ `build_temperature_ranges()` и `build_temperature_gibbs_cache()`

**Статус:** 4 теста в `equilibrium_log_moles_tests2.rs` (`build_temperature_ranges_returns_single_range_on_empty_solver`, `build_temperature_gibbs_cache_returns_empty_on_empty_solver`, `build_temperature_point_index_returns_empty_for_empty_ranges`, `build_temperature_point_index_maps_temperatures_to_range_ids`).

#### D.25 ✅ `solve_temperature_point_from_seed()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`solve_temperature_point_from_seed_fails_on_empty_seed`).

#### D.26 ✅ `snapshot_published_solution()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`snapshot_published_solution_fails_without_publication`).

#### D.27 ✅ `ordered_gibbs_functions()` и `build_gibbs_functions()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`ordered_gibbs_functions_returns_error_for_missing_substance`, `ordered_gibbs_functions_returns_empty_for_empty_input`, `build_gibbs_functions_returns_empty_on_empty_solver`).

#### D.28 ✅ `build_parallel_gibbs_functions()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`build_parallel_gibbs_functions_returns_empty_on_empty_solver`).

#### D.29 ✅ `reconstructed_mole_state()`

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`reconstructed_mole_state_rejects_dimension_mismatch`, `reconstructed_mole_state_rejects_non_finite_log_moles`).

#### D.30 ✅ `check_task()` с некорректными данными

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`check_task_fails_on_empty_solver_without_settings`).

#### D.31 ✅ `has_rst_symbolic_context()`

**Статус:** 2 теста: 1 в `equilibrium_log_moles_tests.rs` (`has_rst_symbolic_context_returns_false_for_empty_solver`), 1 в `equilibrium_log_moles_tests2.rs` (`has_rst_symbolic_context_returns_false_for_empty_solver`).

#### D.32 ✅ `validate_temperature_range()` с граничными значениями

**Статус:** 4 прямых теста в `equilibrium_log_moles_tests2.rs` + 7 косвенных через `solve_for_T_range`/`par`/`par2` rejection тесты.

#### D.33 ✅ `EquilibriumLogMoles::new()` и `EquilibriumLogMoles::empty()`

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`empty_constructor_creates_default_solver`, `new_constructor_creates_default_solver`).

#### D.34 ✅ `set_initial_guess()`

**Статус:** 1 тест в `equilibrium_log_moles_tests.rs` (`set_initial_guess_validates_length`).

#### D.35 ✅ `clear_published_state()`

**Статус:** 2 теста: 1 в `equilibrium_log_moles_tests.rs` (`clear_published_state_resets_all_publication_fields`), 1 в `equilibrium_log_moles_tests2.rs` (`clear_published_state_clears_all_fields`).

#### D.36 ✅ `validate_mutable_problem_shape()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`validate_mutable_problem_shape_accepts_valid_input`, `validate_mutable_problem_shape_rejects_non_finite_moles`, `validate_mutable_problem_shape_rejects_negative_moles`).

#### D.37 ✅ `stage_equilibrium_system()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`stage_equilibrium_system_panics_on_empty_solver`, `#[should_panic]`).

#### D.38 ✅ `resolved_initial_guess()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`resolved_initial_guess_derives_seed_from_n0_when_no_explicit_guess`, `resolved_initial_guess_uses_explicit_seed_when_provided`, `resolved_initial_guess_rejects_dimension_mismatch`).

#### D.39 ✅ `build_solve_contract()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`build_solve_contract_fails_on_empty_solver`).

#### D.40 ✅ `temperature_worker_seed()`

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`temperature_worker_seed_getters_return_expected_values`).

#### D.41 ✅ `publish_solve_candidate()` с частичными данными

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`publish_solve_candidate_updates_all_publication_fields`, `publish_solve_candidate_accepts_missing_keq_validation`).

#### D.42 ✅ `solver_impl()` с разными политиками

**Статус:** 1 тест в `equilibrium_log_moles_tests2.rs` (`solver_impl_fails_on_empty_solver_with_empty_policy`).

#### D.43 ✅ `solve_backend_cascade()` с пустым списком бэкендов

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`solve_backend_cascade_rejects_empty_backends`, `solve_backend_cascade_rejects_zero_budget`).

#### D.44 ✅ `recoverable_backend_failure_kind()` со всеми вариантами ошибок

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`recoverable_backend_failure_kind_returns_none_for_invalid_problem`, `recoverable_backend_failure_kind_returns_solver_for_solve_error`).

#### D.45 ✅ `compute_species_moles()` (свободная функция)

**Статус:** 3 теста в `equilibrium_log_moles_tests.rs` (`compute_species_moles_rejects_non_finite_log_moles`, `compute_species_moles_rejects_underflow_to_zero`, `compute_species_moles_round_trips_positive_moles`).

#### D.46 ✅ `compute_element_totals()`

**Статус:** 3 теста в `equilibrium_log_moles_tests.rs` (`compute_element_totals_matches_manual_calculation`, `compute_element_totals_rejects_dimension_mismatch`, `compute_element_totals_rejects_non_finite_input`).

#### D.47 ✅ `reaction_standard_gibbs()`

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`reaction_standard_gibbs_computes_weighted_sum`, `reaction_standard_gibbs_handles_multiple_reactions`).

#### D.48 ✅ `equilibrium_scaling()`

**Статус:** 2 теста в `equilibrium_log_moles_tests.rs` (`equilibrium_scaling_returns_positive_factors`, `equilibrium_scaling_rejects_invalid_temperature`).

#### D.49 ✅ `species_to_phase_map()`

**Статус:** 3 теста в `equilibrium_log_moles_tests.rs` (`species_to_phase_map_assigns_each_species_to_its_phase`, `species_to_phase_map_rejects_out_of_range_species_index`, `species_to_phase_map_rejects_unassigned_species`).

#### D.50 ✅ `reaction_phase_stoichiometry()`

**Статус:** 1 тест в `equilibrium_log_moles_tests.rs` (`reaction_phase_stoichiometry_aggregates_by_phase`).

#### D.51 ✅ `validate_logmole_system_dimensions()`

**Статус:** 5 тестов в `equilibrium_log_moles_tests2.rs` (`validate_logmole_system_dimensions_rejects_element_matrix_mismatch`, `validate_logmole_system_dimensions_rejects_species_phase_mismatch`, `validate_logmole_system_dimensions_rejects_out_of_bounds_phase`, `validate_logmole_system_dimensions_rejects_delta_n_mismatch`, `validate_logmole_system_dimensions_accepts_valid_dimensions`).

#### D.52 ✅ `validate_residual_conditions()`

**Статус:** 4 теста в `equilibrium_log_moles_tests2.rs` (`validate_residual_conditions_accepts_valid_conditions`, `validate_residual_conditions_rejects_non_finite_temperature`, `validate_residual_conditions_rejects_zero_pressure`, `validate_residual_conditions_rejects_negative_reference_pressure`).

#### D.53 ✅ `scaled_residual()` и `scaled_jacobian()`

**Статус:** 2 теста в `equilibrium_log_moles_tests2.rs` (`scaled_residual_applies_row_scaling`, `scaled_jacobian_applies_row_scaling`).

#### D.54 ✅ `scale_residual_rows()` и `scale_jacobian_rows()` с некорректными scale

**Статус:** 7 тестов в `equilibrium_log_moles_tests2.rs` (`scale_residual_rows_rejects_dimension_mismatch`, `scale_residual_rows_rejects_non_positive_scale`, `scale_residual_rows_rejects_non_finite_scale`, `scale_residual_rows_applies_correct_division`, `scale_jacobian_rows_rejects_dimension_mismatch`, `scale_jacobian_rows_rejects_non_positive_scale`, `scale_jacobian_rows_applies_correct_division`).

#### D.55 ✅ `evaluate_equilibrium_logmole_residual()` с вырожденными входными данными

**Статус:** 4 теста в `equilibrium_log_moles_tests2.rs` (`evaluate_equilibrium_logmole_residual_rejects_dimension_mismatch`, `evaluate_equilibrium_logmole_residual_rejects_invalid_temperature`, `evaluate_equilibrium_logmole_residual_rejects_non_finite_log_moles`, `evaluate_equilibrium_logmole_residual_returns_correct_length`).

#### D.56 ✅ `evaluate_equilibrium_logmole_jacobian()` с вырожденными входными данными

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`evaluate_equilibrium_logmole_jacobian_rejects_dimension_mismatch`, `evaluate_equilibrium_logmole_jacobian_returns_correct_shape`, `evaluate_equilibrium_logmole_jacobian_rejects_non_finite_log_moles`).

#### D.57 ✅ `equilibrium_logmole_residual2()` и `equilibrium_logmole_jacobian2()`

**Статус:** 3 теста в `equilibrium_log_moles_tests2.rs` (`equilibrium_logmole_residual2_creates_closure`, `equilibrium_logmole_residual2_rejects_wrong_log_moles_length`, `equilibrium_logmole_jacobian2_returns_correct_shape`).

#### D.58 ✅ `PhaseActivityModel::log_activity()` с разными моделями

**Статус:** Тесты существуют в `equilibrium_activity.rs` (`gas_and_solution_activity_offsets_follow_their_contracts`, `pure_condensed_phase_has_unit_activity`).

#### D.59 ✅ `EquilibriumComponentDescriptor` с пустым веществом

**Статус:** Тесты существуют в `equilibrium_problem_tests.rs`.

#### D.60 ✅ `EquilibriumConstantProblem` с некорректными данными

**Статус:** Тесты существуют в `equilibrium_constant_tests.rs` и `equilibrium_constant_solver_tests.rs`.

#### D.61 ✅ `EquilibriumConstantSolver::solve()` с некорректными границами

**Статус:** Тесты существуют в `equilibrium_constant_solver_tests.rs`.

#### D.62 ✅ `EquilibriumConstantSolver::solve()` с неограниченной задачей

**Статус:** Тесты существуют в `equilibrium_constant_solver_tests.rs`.

#### D.63 ✅ `compare_equilibrium_constant_solutions()` с идентичными решениями

**Статус:** Тесты существуют в `equilibrium_constant_cross_validation.rs`.

#### D.64 ✅ `classify_equilibrium_constant_cross_validation()` со всеми статусами

**Статус:** Тесты существуют в `equilibrium_constant_cross_validation.rs`.

#### D.65 ✅ `MultiphaseEquilibriumLayout` с некорректными PhaseSpec

**Статус:** Тесты существуют в `equilibrium_multiphase_domain_tests.rs`.

#### D.66 ✅ `MultiphaseInitialComposition::from_sparse()`

**Статус:** Тесты существуют в `equilibrium_multiphase_domain_tests.rs`.

#### D.67 ✅ `MultiphaseInitialComposition::element_totals()`

**Статус:** Тесты существуют в `equilibrium_multiphase_domain_tests.rs`.

#### D.68 ✅ `PhaseEquilibriumMetadata::from_resolved()` с пустой фазовой системой

**Статус:** Тесты существуют в `phase_equilibrium_problem_tests.rs`.

#### D.69 ✅ `build_phase_equilibrium_problem()` с некорректными условиями

**Статус:** Тесты существуют в `phase_equilibrium_problem_tests.rs`.

#### D.70 ✅ `solve_resolved_pt()` с разными `PhaseEquilibriumSolveMode`

**Статус:** Тесты существуют в `equilibrium_multiphase_story_tests.rs`.

### E. Прочие замечания

#### E.1 Магические числа

- `1e-8` в [`solve_with_phase_control`](equilibrium_workflows.rs:1617) — seed для активируемой фазы. Не вынесено в константу.
- `1e-12` в [`solve_candidate_from_seed`](equilibrium_log_moles.rs:1921) — допуск для отрицательных молей в feasibility check. Не вынесено в константу.
- `10.0` в [`solve_candidate_from_seed`](equilibrium_log_moles.rs:1937) — множитель для element_balance_tolerance. Не вынесено в константу.

**Статус (26.07.2026):** ✅ **Выполнен.** Все три значения получили
семантические имена рядом с numerical acceptance boundary.

#### E.2 Избыточное логирование

В [`solve_candidate_from_seed`](equilibrium_log_moles.rs:1908-1910) три `info!` подряд с полным дампом initial_guess, n0 и stoich_matrix. Для production это слишком подробно. Должно быть `debug!`.

**Статус (26.07.2026):** ✅ **Выполнен.** Полные дампы переведены на
`debug!`; штатный `info!` больше не засоряется внутренним состоянием solver.

#### E.3 Дублирование кода в `solve_for_T_range`, `solve_for_T_range_par`, `solve_for_T_range_par2`

Три реализации температурной серии содержат значительное дублирование логики обновления Gibbs, валидации и сохранения результатов. Можно было бы выделить общий шаг итерации.

**Статус (26.07.2026):** ❌ **Отклонён в исходной формулировке.** Обычный
serial sweep и serial phase-control sweep используют общие
`refresh_temperature_gibbs` и `publish_temperature_sweep`. Два parallel path
намеренно остаются отдельными: у них другой ownership/continuation contract,
который нельзя скрыть одной общей итерацией.

#### E.4 `PHASE_CONTROL_TRACE_MOLE_FLOOR = 1e-300`

**Где:** [`equilibrium_workflows.rs:31`](equilibrium_workflows.rs:31)

Значение `1e-300` находится на грани представимости в f64 (min positive normal ≈ 2.2e-308). `ln(1e-300) ≈ -690.8`, что далеко от -Inf, но при дальнейшем делении на phase totals может привести к underflow. Рекомендуется `1e-100` или документировать риск.

**Статус (26.07.2026):** ❌ **Отклонён.** `ln(1e-300)` конечен, а floor нужен
именно для inactive phase coordinates. Менять физически чувствительный порог
без stability benchmark нельзя; его назначение задокументировано и покрыто
phase-control regressions.

#### E.5 `DEFAULT_TRACE_MOLE_FLOOR = 1e-30`

**Где:** [`equilibrium_problem.rs:28`](equilibrium_problem.rs:28)

`ln(1e-30) ≈ -69.1`. Это безопасно для f64. Замечаний нет.

**Статус (26.07.2026):** ❌ **Отклонён как проблема.** Диагностика сама
подтверждает безопасность значения; это именованный, валидируемый seed policy.

#### E.6 Неконсистентность: `scaling_flag` в `EquilibriumSolverSettings` по умолчанию `false`

**Где:** [`equilibrium_log_moles.rs:153`](equilibrium_log_moles.rs:153)

Масштабирование выключено по умолчанию, хотя в архитектурном обзоре оно описано как важная фича. Рекомендуется включить по умолчанию.

**Статус (26.07.2026):** ❌ **Отклонён.** Default должен сохранять
историческую численную семантику. Масштабирование включается явным policy
параметром; изменение default требует comparative stability benchmark.

---

### Приоритеты для исправления

| Приоритет | Категория | Пункт | Итог |
|-----------|-----------|-------|------|
| P0 | Bug | A.4 | ❌ Отклонён как локальный патч; нужен отдельный typed migration этап |
| P0 | Design | C.1 | ❌ Отклонён как локальный патч; `EquilibriumProblem` остаётся source of truth |
| P1 | Bug | A.2 | ✅ Выполнен: RST readiness теперь требует полного контекста |
| P1 | Bug | A.1 | ✅ Выполнен: CODATA 2018 `R` |
| P1 | Overengineering | B.1 | ❌ Отклонён: global resource cap независим от per-attempt лимита |
| P1 | Overengineering | B.2 | ✅ Выполнен: мёртвый `VariableScalingContract` удалён |
| P1 | Test | D.1-D.70 | ✅ Выполнено: 70 из 70, включая serial phase-control sweep |
| P2 | Design | C.2 | ❌ Отклонён как локальный патч; требует отдельной error-model migration |
| P2 | Design | C.4 | ❌ Отклонён как локальный патч; часть prepared-problem migration |
| P2 | Code | E.1 | ✅ Выполнен: numerical constants именованы |
| P2 | Code | E.4 | ❌ Отклонён: trace floor нельзя менять без benchmark |
| P3 | Code | E.2 | ✅ Выполнен: подробный logging понижен до `debug!` |
| P3 | Code | E.3 | ❌ Отклонён в исходной форме; serial helpers выделены, parallel semantics отдельны |

Итог: SourceCraft Diagnostics зафиксирован как журнал принятых исправлений и
отклонённых локальных предложений; открытые архитектурные миграции ведутся
отдельными пунктами основного плана.

---

## SourceCraft Diagnostics (дополнение 27.08.2026)

### F. Потенциальные ошибки (bugs) — новые модули P8/P9

#### F.1 Monolithic P,H runner не сходится на синтетическом фикстуре i3

**Где:** [`pure_phase_ph_validation_tests.rs:225`](pure_phase_ph_validation_tests.rs:225)

Тест
`i3_fixed_topology_monolithic_ph_matches_independent_nested_scalar_solution`
падает с `AllBackendsFailed`:

```text
Legacy(LM): nonlinear solver reached its iteration limit
Legacy(NR): nonlinear solver reached its iteration limit
Legacy(TR): nonlinear solver encountered a singular matrix
```

**Причина:** production monolithic `P,H` runner
(`PreparedMonolithicPhRunner::solve_from_temperature_seed`) не сходится на
синтетическом газ+конденсированная фаза фикстуре при `tol = 1e-11`,
`max_iter = 300` c legacy каскадом. Сходимость связана с численной
чувствительностью монолитного `[ln(n), theta_T]` формулирования на данном
дата-наборе, а не с тестом документации.

**Статус:** ❌ **Открыто.** Требуется диагностика:
1. исследовать conditioning монолитного Якобиана на фикстуре;
2. проверить стартовый seed и масштабирование enthalpy residual;
3. рассмотреть усиление каскада через RST бэкенды либо отдельный
   continuation/масштабный ход перед требованием `1e-11`.
---

## Extensive-normalization production-boundary audit (2026-09-01)

This post-implementation audit is limited to correctness and evidence of the
existing `ExtensiveNormalization` recovery. It does not add a solver
algorithm, change tolerances, or widen the recovery trigger.

| Area | Status | Evidence / remaining work |
| --- | --- | --- |
| Public physical-unit result boundary | PASS | Public solution accessors publish physical extensive quantities; normalized values are explicitly diagnostic-only. |
| Numerical versus physical phase totals | PASS | `phase_totals` remains separate from `numerical_phase_totals`; frozen recovery tests verify physical totals and topology. |
| Recovery provenance | PASS | Evidence retains the original failure, discovery backend, physical retry outcome, and reconstruction flag. |
| Disabled policy | PASS | Frozen TP-1907 `P,T` and TP-1906 `P,H` matrices preserve the historical failure when recovery is disabled. |
| Recovery trigger classification | PASS | The workflow matrix test excludes invalid input, cancellation, candidate rejection, domain, and capability errors. |
| Near-unit guard | PASS | Named and tested boundary: recovery requires `scale < 0.1` or `scale > 10.0`. |
| Small-scale symmetry | PASS | Pure transformation tests cover `1e-6`; the real matrix covers `1e-4`. |
| P,T / P,H reconstruction | PASS | Frozen matrices verify physical amounts, conservation, topology, and enthalpy after recovery. |
| Continuation and range progress | PASS | Prepared TP-1907 range evidence records one recovery formulation, two physical reuses, and exact point lifecycle events. Release output is recorded in `STORY_TESTS.md` section 46. |
| Diagnostic leakage | PASS | Formatter explains physical failure, normalized basin discovery, and physical retry/reconstruction without printing normalized mole fields. |
| Options snapshot | PASS | Both recovery policies round-trip through `serde_json`. |
| Timing accounting | PASS | Prepared-range timing includes the rejected attempt and recovery transaction; build/reuse counters remain distinct. |
| Extreme arithmetic | PASS | Unit tests cover `1e-250` and `1e250`; unsafe operations return typed errors. |
| Backend metrics (`iterations=0`, `metrics=None`) | DEFERRED | No normalization defect was found; missing-metrics conventions need a dedicated backend-reporting audit. |
| Recovery cancellation and global budgets | PARTIAL | Ordinary execution-control tests pass; a dedicated cancellation-during-recovery matrix remains useful before release. |
| Unsupported model capability | PASS | Recovery is reached only after the typed numerical-failure gate. |
| Immutable source data / independent solves | PASS | Frozen-reference immutability and independent-solve coverage remain unchanged. |

### Audit conclusion

The reviewed normalization boundary is production-safe for the current ideal
phase-model scope. Remaining items are observability evidence gaps, not known
public-unit or recovery-correctness defects. Future non-ideal models require a
separate extensivity proof before this recovery policy is enabled for them.

---

## P,H target-range transactional extensive-normalization recovery (planned)

**Scope.** Extend the established production recovery boundary from one
fixed-`P,T` transaction to a fixed-pressure `P,H` target range. This is a
numerical-conditioning feature for the current ideal-model scope, not a new
physical model, a new `P,H` solver, or permission to use an accepted
reference solution as a seed.

The existing isolated TP-1906/1907 `P,H` normalization matrix proves the
physical transformation for individual points. The missing work is range
orchestration: a recovered point must be published in physical units and the
next target must continue from that physical accepted state only.

### Contract

- [x] Keep `PhRangeRequest` as the sole owner of target ordering, public
  progress events, accepted continuation, and all-or-nothing range
  publication. The point-level recovery helper must not create a second P,H
  range or publish an internal normalized point.
- [x] Add one private P,H point recovery boundary in `equilibrium_ph_workflow`.
  It must run only after the ordinary physical P,H route has returned a
  classified numerical error and `ExtensiveNormalizationPolicy` is
  `OnNumericalFailure`. Invalid input, cancellation, capability/domain errors,
  and rejected physical candidates remain ordinary errors.
- [x] Reuse `ExtensiveNormalization` and `ExtensiveNormalizationPolicy` rather
  than implementing local `/ scale` and `* scale` arithmetic. Normalize only
  extensive quantities: initial moles, element totals through the existing
  composition transform, total enthalpy target, physical absolute enthalpy
  tolerance, absolute trace/phase thresholds, and any explicitly extensive
  error values. Preserve `P`, `p0`, temperature bounds, temperature itself,
  Gibbs/enthalpy/Cp data, TPD, hysteresis, and dimensionless tolerances.
- [x] Solve the exact normalized P,H equivalent with recovery recursively
  disabled and with a local normalized prepared state. A normalized trial must
  never mutate the physical prepared template, physical phase history, or the
  range's currently accepted continuation state.
- [x] After normalized acceptance, prefer a physical-coordinate retry seeded
  from the reconstructed normalized answer. If that retry is still numerically
  ill-conditioned, publish only through an audited reconstruction boundary.
  Both routes must preserve the original physical failure, normalized discovery
  backend, retry outcome when present, and the chosen publication route.
- [ ] Add/extend immutable P,H recovery evidence. It must distinguish internal
  normalized enthalpy error from reconstructed physical-Joule error, preserve
  the physical inventory scale and physical/normalized targets for labelled
  diagnostics, and never expose normalized quantities through ordinary public
  solution accessors.

### Continuation and transactions

- [x] A recovered first target remains one physical `Initial` range point even
  when it used a normalized internal formulation. The next target must receive
  the recovered physical component moles, physical phase set/history, and
  accepted physical temperature.
- [x] Subsequent successful targets remain `Continued`; do not overload this
  classification to mean "normalized recovery". Add a separate recovery/formula
  axis to point evidence if the existing report cannot represent both facts
  without ambiguity.
- [x] Guard against inherited/double normalization: a following target may use
  recovery only after *its own* classified physical failure. The predecessor's
  recovery provenance must not become a solver mode flag.
- [ ] Count a rejected physical attempt plus normalized discovery/physical
  retry in timing and formulation evidence without increasing the public
  physical point count. Preserve existing backend-attempt, iteration, global
  budget, cancellation, and range rollback semantics.
- [ ] Public range progress for four accepted physical targets must remain
  exactly four `PointStarted` and four `PointAccepted` events. Backend attempts,
  scalar temperature trials, and normalized discovery are internal diagnostics,
  never additional range points.

### Frozen NASA regression matrix

- [x] Add an ignored release story using frozen NASA TP-1906/1907 CHON +
  graphite at `P = 101325 Pa`, physical inventory factor `1e4`, and the
  established scaled targets `H680 -> H700 -> H720 -> H740`.
- [ ] Assert source/frozen/local-JSON byte snapshots are unchanged. Start H680
  from a genuine large physical request; do not seed it from a unit-scale
  solution. The implementation must record actual recovery locations rather
  than hard-coding an expected count.
- [x] Assert physical continuation semantics: first point `Initial`, later
  points `Continued`; every continuation seed equals the preceding accepted
  physical solution; no partial suffix is published if a point and its recovery
  both fail.
- [ ] Compare every accepted range point with an independently solved unit-scale
  P,H baseline and a same-inventory P,T witness at recovered temperature.
  Check temperature, component moles after scaling, element balances, target
  enthalpy/residual in physical units, graphite amount, and phase topology.
- [ ] Require the forward topology `gas+graphite, gas+graphite, gas, gas` and
  an accepted graphite disappearance between H700 and H720. Keep nested trial
  phase events separate from accepted range topology transitions.
- [ ] Add an independent reverse range `H740 -> H720 -> H700 -> H680`. Its
  first point is `Initial`; compare matching physical targets forward/reverse
  within the established production envelope while retaining intentional
  hysteresis semantics near the boundary.
- [ ] Add a companion `ExtensiveNormalizationPolicy::Disabled` witness. It
  must retain the historical typed numerical failure at the difficult physical
  point, emit no `PointAccepted` for it, and publish no partial range.
- [ ] Add a focused cancellation-during-recovery test if the existing execution
  hooks can inject it without a production-only test hook. Otherwise leave this
  explicitly deferred with the current ordinary cancellation coverage noted.

### Current implementation note

The point recovery boundary and the frozen two-point story are now implemented.
The range commits continuation only after a physical publication succeeds; a
failed prepared monolithic point also restores its mutable RST parameter buffer.
The remaining unchecked items above are deliberate follow-up work: richer
normalized-versus-physical error fields in the immutable report, global budget
accounting across all recovery attempts, the full four-point forward/reverse
release matrix, and cancellation injected during recovery.

### Observability and release evidence

- [ ] Extend P,H range reports/presentation/diagnostics so a recovered point
  explains: physical failure, normalized-equivalent acceptance, physical retry
  or reconstruction, scale, labelled internal versus physical enthalpy error,
  and retained typed backend evidence. Do not print unlabeled normalized moles.
- [ ] Add a compact ignored release table with target/source temperature,
  preparation, recovered temperature, topology, accepted transitions, nested
  trial events, physical and normalized enthalpy errors, recovery provenance,
  build/reuse counts, and public progress summary. Record its release output in
  `STORY_TESTS.md`.
- [ ] Treat missing backend iteration metrics as separate reporting debt. This
  work must preserve existing typed attempts but must not redesign backend
  metrics merely to implement P,H normalization recovery.

---

## I5: competing water phase candidates and declaration-order invariance

This stage is a feasibility-first validation of the existing production phase
control machinery. It must not change phase-control policy, thresholds, local
thermochemistry, or candidate sorting before a reproducible order-dependent
failure is observed.

### Feasibility and data preflight

- [x] Identify the exact local records and phase-qualified identifiers for
  `Ar(g)`, `H2O(g)`, `H2O(l)`, and ordinary ice Ih (`H2O(s)`).
- [x] Record local provenance: `NASA_gas::Ar`, `NASA_gas::H2O`,
  `NASA_cond::H2O(L)`, and `NASA_cond::H2O(s)`. The resolved common local
  temperature domain is the singleton endpoint `273.150000 K`; it is not a
  usable below-triple-point interval.
- [x] Confirm the frozen IAPWS row coverage: liquid saturation rows begin at
  `275 K`, while ice-Ih sublimation rows end at `270 K`. There is no shared
  frozen-IAPWS row near `273.15 K`.
- [x] Record the available reference-pressure metadata: all four selected
  local records currently publish `Undeclared` standard-state pressure
  metadata. The equilibrium request therefore continues to require an
  explicit reference pressure; no pressure convention is inferred from the
  lookup report.
- [ ] Expose full per-record validity intervals in the characterization
  report before selecting a lifecycle temperature; the current preflight only
  publishes the resolved common interval.
- [ ] Reuse the existing frozen IAPWS liquid-saturation and ice-sublimation
  datasets/helpers; do not create duplicate external data.
- [x] Stop the first feasibility pass before phase-control: the local common
  domain is only a singleton and the external frozen row sets do not overlap.
  Never silently extrapolate a thermochemical record. A lifecycle test now
  requires a separately documented source-faithful bridge or a revised local
  record domain.
- [x] Add a test-only IAPWS gauge at `273.15 K`: the independent liquid and
  ice boundaries are `611.212846 Pa` and `611.153475 Pa`; the resulting Ar/H2O
  gas-only state has both condensed candidates below `dg_create`, with ice
  preferred by `2.205e-1 J/mol`. This is an external characterization gauge,
  not a claim that the local NASA records provide a below-triple-point model.
- [x] Run the production active-set runner for `[gas, liquid, ice]` and
  `[gas, ice, liquid]`. Both permutations expose negative initial canonical
  TPD candidates and converge to the same physical fixed point `gas + ice`;
  comparison is by component identity and conservation, not transition order.
- [ ] Add gas+liquid and gas+ice initial-history permutations after the
  supported public API for constructing accepted active histories is selected.

Preflight output (debug, 2026-09-07):

```text
Ar + H2O competing-phase local preflight
common local temperature bounds: 273.150000..273.150000 K
NASA_gas::Ar       record=Ar    standard-state pressure=Undeclared
NASA_gas::H2O      record=H2O   standard-state pressure=Undeclared
NASA_cond::H2O(s)  record=H2O(s) standard-state pressure=Undeclared
NASA_cond::H2O(L)  record=H2O(L) standard-state pressure=Undeclared
```

---

### Controlled condensed-phase Gibbs degeneracy

- [x] Add a test-only controlled Gibbs family in which only
  `G0_B - G0_A = delta_G` changes; keep it separate from frozen IAPWS
  chemistry and production APIs.
- [x] Characterize `delta_G` from `1e0` through `1e-6 J/mol` plus exact zero.
  For all positive splittings both declaration orders select the same winner,
  preserve the physical state, and reproduce the analytic losing TPD.
- [x] Treat `delta_G=0` as a mathematically non-unique representative: do
  not require a phase label, while requiring finite termination,
  complementarity, conservation, and equal physical state after permutation.
- [x] Add a separate observed-TPD-scatter measurement for repeated and
  permuted well-resolved runs. The measured value is characterization for this
  fixture only; it is not promoted to a production threshold. Debug and
  release both measured zero scatter for the current deterministic fixture.
- [ ] Classify the first tolerance-limited splitting using that evidence. Do
  not infer a universal solver resolution from this one fixture.
- [x] Add the requested `1e-6..1e-15` plus exact-zero resolution matrix,
  including requested versus representable `f64` splitting, both declaration
  orders, measured losing TPD, relative error, physical-state delta,
  complementarity, transitions, and explicit floating-point-collapse labels.
- [ ] Review the recorded matrix output and document the empirical boundaries
  between quantitative resolution, order-independent winner selection,
  tolerance-limited behavior, and floating-point collapse. Debug and release
  agree: quantitative and unique-winner boundary `1e-7`, first
  tolerance-limited row `1e-8`, and f64 collapse from `1e-13`. These boundaries
  remain fixture-specific characterization, not production policy.

### I5: vanishing species inside one active phase

- [x] Add an analytic synthetic ideal-gas fixture with five species sharing
  one conserved element and exact Gibbs-derived target weights.
- [x] Cover dynamic ranges through `1e-32` and a dedicated vanishing-species
  sweep through `1e-40`. Check finite positive moles, element balance, and
  scale-aware log errors rather than one absolute mole tolerance.
- [x] Confirm that tiny species do not trigger phase disappearance: every
  accepted result retains one active gas phase.
- [x] Record release evidence for the three baseline tests. The `1e-40`
  vanishing-species row remains finite and tracks the analytical amount with
  balance error below `3e-15`.
- [x] Add seed and species-permutation matrices. Both preserve the canonical
  physical composition and one active gas phase.
- [x] Add a solver-independent Gibbs-weight representability matrix through
  ratio `1e-100`; the requested ratios remain representable in `f64`, so this
  fixture does not confuse a floating-point collapse with a solver failure.
- [x] Assert chemical-potential equality directly from the returned log-mole
  state, in addition to mole-ratio and conservation checks. The observed
  deviations are below `1e-6 J/mol` for the tested dynamic ranges.
- [x] Add an extensive-scaling characterization for `N=1e-8, 1, 1e8 mol`.
  The small-inventory row is finite and conservative but currently exposes a
  scale-limited fraction error of about `8.5e-5`; do not conflate this with a
  species floor or silently label it invariant.
- [ ] Extend the supported-backend matrix and measure a solver-specific
  numerical floor. The current evidence reaches `1e-40` in the vanishing
  sweep and `1e-100` in fixture representability, but the extensive `1e8 mol`
  case still reports `AllBackendsFailed`; this remains characterization, not a
  production-policy change.

## I5: ideal Ar/H2O/CO2 three-phase lifecycle

This proposed `gas -> gas + solid -> gas + ice + CO2(s) -> ... -> gas`
story is blocked at the data-feasibility boundary, not at phase control. It
must remain an offline, source-faithful fixture: online NIST lookup and guessed
species names are not substitutes for a local dry-ice thermochemistry record.

- [x] Inspect the bundled catalog before constructing a reaction basis or
  touching phase-control policy. `CO2(s)` is present only in synthetic
  candidate-selection tests, not in the runtime catalog.
- [x] Add a negative offline lookup regression: the declared local
  `NASA_gas::Ar/H2O/CO2`, `NASA_cond::H2O(s)`, `NASA_cond::CO2(s)` system
  deterministically returns `SubstanceNotFound("CO2(s)")` with NIST fallback
  disabled.
- [x] Verify the broader index: local `CO2` entries are gas records in
  `Cantera_nasa_base_gas`, `CEA`, and `nuig_thermo`; no indexed
  `NASA_cond`/other condensed `CO2` record is available. Do not relabel a gas
  record as dry ice.
- [x] Identify a second independent feasibility blocker: the local
  `NASA_gas::Ar/H2O` plus `NASA_cond::H2O(s)` interval begins at `200 K`, while
  the supplied frozen CO2(s) Antoine correlation ends at `195.89 K`. Adding
  dry ice alone would therefore not establish a common interval for the
  initially proposed external oracle.
- [x] Review the official NIST WebBook CO2 phase-change evidence. It records
  `T_triple = 216.58 K`, `P_triple = 5.185 bar`, and
  `Delta_sub H(207 K) = 26.1 kJ/mol` (the latter is based on `198..216 K`
  data), but its published Antoine correlation remains only `154.26..195.89
  K`. A triple-point datum and one sublimation-enthalpy anchor are not a
  source-complete `G0_CO2(s, T)` closure or a 200--216 K pressure oracle.
- [x] Confirm the runtime gas-side constraint independently: the selected
  local `NASA_gas::CO2` record is declared for `200..6000 K`, so it does not
  overlap the Giauque--Egan/NIST `154.26..195.89 K` sublimation-pressure range
  either.
  A physically closed 154--195 K CO2(s) characterization can still be useful
  I5 evidence, but cannot be connected to the current runtime gas record.
- [x] Establish the reference-pressure contract: production activities use
  the explicit `EquilibriumConditions::reference_pressure()` supplied with a
  solve; there is no safe hidden global `p0`. Any future frozen or runtime
  CO2(s) closure must therefore carry a declared source pressure and be
  converted only through `RT ln(p/p0)` with the request's explicit `p0`.
- [ ] Build a complete, source-faithful CO2 gas/solid closure before adding a
  runtime dry-ice record or revisiting the Ar/H2O/CO2 lifecycle. This replaces
  the earlier provisional `CalorimetricSolid` / sublimation-reconstruction
  plan: a primary Gibbs EoS is available for the solid, so independently
  fitting `Cp`, `H`, `S`, and `G` would make a weaker and internally
  inconsistent production representation.
  - [x] Freeze the reviewed low-temperature rows from the official
    natural-abundance `CO2-total.txt` table from
    Tashkun--Harvey, JPCRD 54, 023102 (2025), DOI `10.1063/5.0276615`,
    NIST dataset DOI `10.18434/mds2-2364`. It supplies `Cp`, `S`, `H`, and
    uncertainties on `1..6000 K` at 1 K increments. Record its SHA-256,
    license/provenance, natural-isotope composition, SI units, and the
    declared standard pressure `1e5 Pa`. The small source-faithful window is
    deliberate: the full 1..6000 K table remains recoverable from the pinned
    NIST URL/hash rather than being duplicated in the test repository.
  - [ ] Inspect the Tashkun--Harvey enthalpy zero convention before joining it
    to KiThe. Compare its `Cp/S/H/G` with the selected local `CO2(g)` record
    at 200 K and 298.15 K; a constant formation-enthalpy offset is not a
    harmless detail because gas/solid chemical-potential differences require
    one common standard-state datum.
  - [ ] Keep exactly one production `CO2(g)` identity. Add a reviewed
    low-temperature segment, either as direct tabulated interpolation or as a
    `Cp(T)`-only fit analytically integrated to `H/S`, with integration
    constants fixed by the existing canonical record at `T_join = 200 K`.
    Never independently fit `Cp`, `H`, `S`, and `G`, and never extrapolate the
    existing NASA segment below its declared 200 K lower bound.
  - [ ] Add a splice characterization at `195, 198, 199, 200, 201, 205 K`.
    It must report `Cp/H/S/G` and left/right discontinuities. `H/S/G` must be
    continuous by construction; any `Cp` difference is source evidence, not
    a value to hide by altering the anchor.
  - [ ] Freeze the NIST-JANAF `CO2(g)` 100 K and 200 K rows as independent
    anchors (`p0 = 0.1 MPa`, `Tr = 298.15 K`), and characterize rather than
    demand bitwise agreement with the newer 2025 evaluation.
  - [ ] Implement a test-only, source-faithful adapter for the primary Siah,
    Campestrini, Stringari dry-ice-I Gibbs EoS, JCED 70, 2890--2905 (2025),
    DOI `10.1021/acs.jced.5c00260`. Freeze the reviewed equation, parameter
    values, units, standard-pressure convention, validity domain
    (`T <= 400 K`, `p <= 1.2 GPa`), and source digest; do not manufacture a
    solid NASA polynomial from phase-change anchors. Source review confirmed
    the article's model scope and the existence of analytic derivatives in its
    supporting information, but the actual equation/parameter tables must be
    acquired as a reviewed frozen source artifact before code is written.
  - [ ] Derive every solid standard-state function from the one Gibbs EoS at
    its declared `p0`: `G0 = g(T,p0)`, `S0 = -dg/dT`, `H0 = G0 + T*S0`, and
    `Cp0 = dH/dT`. Prefer the source's analytic derivatives from supporting
    information; validate identities numerically only as a secondary check.
  - [ ] Quantify `g_s(T,p_actual) - g_s(T,p0)` over the intended low-pressure
    fixture range. The ideal-pure-condensed adapter is permitted only when
    that omitted pressure correction is demonstrably small relative to the
    phase-stability/TPD scales being asserted; otherwise this is a missing
    production physics capability, not a tolerance problem.
  - [ ] Retain Giauque--Egan (1937) and NIST WebBook Antoine/sublimation data
    as independent frozen evidence, never as fitted input for the new solid
    closure. On `160, 170, 180, 185, 190, 194, 195 K`, reconstruct
    `p_sub = p0 * exp(-(G0_g - G0_s)/(R*T))` and separately compare
    `H0_g - H0_s` at 194.67 K against `6030 cal/mol` (about 25.23 kJ/mol).
  - [ ] Add one explicit closure report and verdict taxonomy:
    `CO2GasLowTemperatureClosed`, `CO2SolidPhaseClosed`, and
    `CO2GasSolidClosureValidated`; otherwise emit a specific gas splice,
    solid-Gibbs, sublimation, or reference-pressure mismatch. The report must
    include gas/solid `Cp/H/S`, reconstructed and external sublimation
    pressure, `Delta_H_sub`, and all source provenance.
  - [ ] Only after all of the above passes, promote the adapter through the
    thermochemistry/repository boundary as a reviewed local `CO2(s)` record,
    then resume the Ar/H2O/CO2 lifecycle feasibility study. The lifecycle is
    not a substitute for validating the CO2 standard-state closure.
- [ ] Once all five local records resolve, perform the required independent
  scalar preflight before any production lifecycle test: choose a finite
  one-solid/two-solid temperature interval using only local chemical-potential
  boundaries, then compare canonical TPD and phase-order permutations.
- [ ] Keep the final lifecycle as a phase-order-invariance regression only
  after that feasibility gate is passed. Candidate enumeration may affect the
  transition history but must not affect the physical fixed point.

### Test-only IAPWS/NIST sublimation-gauge lifecycle

This is a separate, explicitly non-production P,T validation route. It tests
the universal multiphase lifecycle `1 -> 2 -> 3 -> 2 -> 1` while the complete
CO2 gas/solid production closure remains deferred above. The gauge must never
enter `SubsData`, `ThermoRepository`, JSON libraries, or P,H workflows.

- [x] Reuse the official IAPWS ice-Ih sublimation equation on `50..273.16 K`
  and freeze/reuse the NIST CO2 Antoine sublimation relation on
  `154.26..195.89 K`; use the shared explicit test convention `p0 = 100000 Pa`.
- [x] Define only relative standard Gibbs functions: all gas `G0 = 0`,
  `G0_ice = R*T*ln(p_sub_H2O/p0)`, and
  `G0_CO2_s = R*T*ln(p_sub_CO2/p0)`. Document that this is unsuitable for
  P,H and has no formation-energy interpretation.
- [x] Prove gauge algebra independently: each gas/solid equality reproduces
  its source sublimation pressure before invoking TPD or a nonlinear solver.
- [x] Implement a scalar ideal oracle for gas-only, one-solid, and
  gas+ice+dry-ice states. It must use only `T/P`, inventory, and both external
  correlations, and determine a non-degenerate `160..194 K` design with a
  finite three-phase interval. The initial design is `P=100000 Pa`,
  `n_Ar=1 mol`, `n_H2O=1e-7 mol`, `n_CO2=1 mol`: it yields gas-only at
  194 K, gas+dry-ice at 185 K, and gas+ice+dry-ice at 170 K with both
  condensed amounts materially above their respective inventory floors.
- [x] Run the ordinary production prepared phase-control runner over forward
  cooling and reverse heating. The accepted continuation path proves topology
  changes `1 -> 2 -> 3` and `3 -> 2 -> 1`, with point-by-point agreement with
  the independent identity-aligned scalar oracle.
- [x] Prove that the three-phase physical fixed point is invariant under both
  `[gas, ice, dry_ice]` and `[gas, dry_ice, ice]` declaration orders. This
  compares physical species amounts and maps active masks back to phase
  identity; it deliberately does not constrain transition history.
- [ ] Add direct chemical-potential equality assertions and canonical TPD
  boundary values against the independent sublimation equations.
- [x] Assert finite canonical TPD evidence, zero empty/chattering transition
  records, and fresh-versus-continued agreement at the selected interior
  points. The test also requires exactly two accepted transitions in each
  direction.
- [ ] Keep a failure taxonomy (`ThirdPhaseActivationFailure`, fixed-point,
  premature disappearance, chatter, order dependence, continuation dependence,
  scalar-oracle mismatch). Do not change phase thresholds, sorting, or solver
  policy without a separately recorded negative characterization.

---

## Transactional rollback during failed multicomponent continuation

This stage is separate from invalid-input/range rollback. Its target is a
real TP-1906/1907 CHON + graphite continuation in which point `S1` has already
been accepted, point `S2` has entered meaningful fixed-active/phase-control
work, and a numerical failure is followed by a supported successful route.

### Feasibility before production changes

- [x] Reuse the existing real TP-1906/1907 fixture and avoid a synthetic
  multicomponent surrogate.
- [x] Add an ignored iteration-budget matrix that distinguishes a failure at
  point zero from a failure after accepted continuation points.
- [x] Run the debug feasibility matrix. Budgets `1, 2, 4, 8` fail at point
  zero; budgets `16, 32, 64, 100` accept all three points (`700, 720, 740 K`).
  No natural post-acceptance failure was found, so this matrix is not itself
  rollback evidence.
- [x] Repeat the feasibility matrix in release before finalizing the evidence
  classification. The release result matches debug: no natural
  post-acceptance failure. Do not call an initial-point failure a rollback
  scenario.
- [x] If no natural failure can be reproduced, document the tested budgets,
  backend policy, boundary locations, and conclude whether a minimal
  `cfg(test)` failpoint is justified. The real transition-before-commit
  boundary now has a one-shot test-only failpoint; production builds do not
  contain it.

Debug/release characterization output (2026-09-07):

```text
NASA TP-1906/1907 continuation rollback feasibility
budget  status       failing_point  accepted_points  detail
     1  FAILED                   0               0  post_acceptance=false
     2  FAILED                   0               0  post_acceptance=false
     4  FAILED                   0               0  post_acceptance=false
     8  FAILED                   0               0  post_acceptance=false
    16  OK           -                            3  complete
    32  OK           -                            3  complete
    64  OK           -                            3  complete
   100  OK           -                            3  complete
natural post-acceptance failure found: false
```

The release run completed successfully with the same rows and verdict. Since
the natural failure search is exhausted, the next evidence layer may use only
a minimal `cfg(test)` failpoint at the real transition-before-commit boundary;
it must not be a generic callback failure at solve entry.

### Transactional evidence once a failure exists

- [x] Exercise failure after a real local phase transition has been assembled
  and diagnosed but before its restart seed is committed. The focused runner
  regression confirms that the transition event remains observable while the
  accepted continuation seed and phase set are restored.
- [x] Compare clean, failed-then-recovered, and fresh routes at `S1` and `S2`
  on the real TP-1906/1907 CHON + graphite range; compare physical state, not
  storage indices or diagnostic counters.
- [x] Verify the real range retry is not poisoned by the rejected trial. The
  one-shot failure is consumed, the retry publishes all points, and its
  physical results match both clean and fresh routes.
- [ ] Extend the real route comparison to explicit hysteresis-history and
  prepared-cache fingerprints once those internals have a stable report-level
  representation; the current public contract compares accepted physical
  state and typed range preparation evidence.
- [ ] Preserve rejected attempt diagnostics and trial events without publishing
  rejected transitions as accepted physical history.
- [ ] Verify repeated failed transactions are physically idempotent.
- [x] Add focused release story evidence and record the verdict
`TransactionalRollbackInvariant` in `STORY_TESTS.md`. The recorded release
run accepted the `A/B/C` matrix in `8.43 s`.

## NASA CEA RP-1311 Example 3 H,P frozen regression

This is a test-only frozen external regression. The source layer is complete;
the production comparison remains deliberately separate from source loading.

- [x] Identify the source problem as NASA CEA RP-1311 Example 3: explicit CEA
  Air, liquid fuel split `0.4 C7H8(L) + 0.6 C8H18(L)`, `O/F=17`, reactants at
  `700 K` and `298.15 K`, and pressure points `100/10/1 bar`.
- [x] Count the supplied non-trace source universe: 40 species. Every one has
  an exact spelling match in local `NASA_gas` (`Ar` through `OH`) with NASA7
  Gibbs, enthalpy, and heat-capacity coefficients. Their common interval is
  `200..6000 K`, which covers all three published equilibrium temperatures
  `2418.660/2390.593/2338.840 K`.
- [x] Resolve the two liquid reactants offline. `C7H8(L)` is an exact local
  `NASA_cond` key with interval `178.15..500 K`; n-octane is present as the
  explicit local key `C8H18(L),n-octa` with interval `220..300 K`. The latter
  is a `ReviewedAlias` for the CEA identity `C8H18(L)`, not an invented
  fallback, and it covers `298.15 K`.
- [x] Record the exact CEA Air elemental convention without replacing it by
  generic dry air: `N=1.561680`, `O=0.419590`, `Ar=0.009365`, `C=0.000319`
  per nominal mole of Air.
- [x] Confirm that all above-threshold gas records have the thermodynamic
  capabilities required by the planned ideal-gas H,P route. The current raw
  NASA records do not carry an explicit standard-state-pressure field; the
  future fixture must therefore state and test its pressure convention
  explicitly rather than infer one from absent metadata.
- [x] Obtain and freeze the exact published CEA mole fractions for all three
  pressure rows. The task description supplies temperatures and total
  enthalpy, but not the numerical 40-species composition table; no frozen
  rows or composition regression may be fabricated from that summary.
- [x] Inspect and classify the CEA trace-species universe below `1e-15` once
  the complete source output is available. Missing trace-only records may be
  documented as negligible only after their source abundance is known.
- [x] Add the typed immutable dataset adapter and register its metadata pair
  in the repository frozen-reference catalog. The source preflight passes
  with three rows, 40 non-trace identities, and 39 trace-only identities.
- [x] Reconstruct the physical feed independently on a `1 kg fuel + 17 kg
  explicit CEA Air` basis. The resulting total mass is `18 kg`, elemental
  totals are positive for `C/H/O/N/Ar`, and the frozen source enthalpy gives
  `H_target=5.721084e6 J` on this extensive basis.
- [x] Resolve the two pure condensed fuel records offline and evaluate their
  local enthalpies at `298.15 K`: `C7H8(L)=1.217868e4 J/mol` and
  `C8H18(L),n-octa=-2.502669e5 J/mol`. Their mass-weighted one-kilogram fuel
  contribution is `-1.261649e6 J`.
- [x] Add the test-only atom-balanced CEA Air decomposition
  `N2=0.780840`, `O2=0.209476`, `Ar=0.009365`, `CO2=0.000319`. It reproduces
  the four published elemental constraints and sums to one nominal mole.
- [x] Resolve those four ordinary NASA gas records offline, with no NIST
  fallback, and evaluate their enthalpies at `700 K` with explicit local
  provenance.
- [x] Characterize the complete source-side reactant enthalpy. The local
  result is `317813.0918 J/kg` versus CEA `317838.0000 J/kg`, delta
  `-24.9082 J/kg`, relative `7.83675e-5`; classify this as
  `ReactantEnthalpyConventionAligned`.
- [x] Add an ignored production-path H,P characterization over all three
  frozen pressure rows. It uses the canonical `solve_resolved_ph` route,
  local `NASA_gas` only, source CEA enthalpy as the frozen target, and checks
  finite temperature, normalized published composition, residual, and element
  balance. The test-side initial seed uses an element-equivalent `CO + H2`
  representation because the local gas catalog has no octane gas record; it is
  not a production Air/fuel alias.
- [x] Correct the fixture standard-state pressure to the fixed NASA gas value
  `P0=100000 Pa` and assert `P/P0 = 100, 10, 1` for the three rows.
- [x] Correct the liquid-fuel unit error. `reconstruct_feed()` now uses SI
  `kg/mol`; the production fixture resolves the local liquid record masses,
  explicitly converts `g/mol -> kg/mol`, and reconstructs `1 kg` fuel,
  `17 kg` Air, and `18 kg` total independently. It prints `n_C7H8`,
  `n_C8H18`, fuel C/H atoms, nominal Air, elemental totals, and H target.
  The CO/H2 computational seed is checked against the locally reconstructed
  C/H/O/N/Ar totals before solving.
- [x] Re-run all three rows after both fixture corrections. KiThe gives
  `2419.501/2391.538/2339.773 K` versus CEA
  `2418.660/2390.593/2338.840 K`, or `+0.841/+0.945/+0.933 K`. Maximum
  major-species relative error is below `9.54e-2`; minor/trace maximum
  `|delta_log10|` is below `9.97e-1`; release residuals are below
  `6.93e-14` and balances below `9.7e-12`.

Current source verdict: `ExternalFixtureSourceComplete`; production route
verdict: `ReactantEnthalpyConventionAligned`, `ProductionPHCharacterized`.
The early constant-temperature result is classified
`ExternalFixtureStandardStatePressureError`; the later `632.79 K` result is
classified `ExternalFixtureFuelAmountUnitError`. After both corrections the
fixture agrees with the three CEA temperatures to under `1 K`; no external
composition mismatch remains to investigate from this evidence.
No production algorithm, thermochemistry, tolerance, scaling, or solver policy
was changed.

### Independent boundary and candidate evidence

- [x] Add a characterization-only catalog preflight for
  `H2O(g) <=> H2O(l)` and `H2O(g) <=> H2O(ice Ih)` independent of the
  production TPD minimizer. The local NASA domain is only the singleton
  `273.15 K`, and the frozen liquid/ice tables do not overlap; no local
  boundary extrapolation is performed.
- [x] Select `T_test=273.15 K` only for the independent IAPWS gauge after
  checking numerical separation: `p_liquid=611.212846 Pa`,
  `p_ice=611.153475 Pa`, and the candidate Gibbs-force separation is
  `2.205e-1 J/mol`, well above the diagnostic threshold used by the test.
- [x] Construct a gas-only Ar/H2O inventory with both absent condensed phases
  satisfying `TPD_liquid < dg_create` and `TPD_ice < dg_create`; print both
  values and their separation using the test-only IAPWS gauge.
- [x] Add the independent ideal scalar `gas + ice` oracle from the ice boundary,
  including gas water amount, ice amount, composition, and element balance.
- [x] Compare the independent candidate driving forces with the canonical
  initial TPD reports using explicit signs and units.

### Metamorphic phase-order matrix

- [x] Run physically identical production bounded phase-control requests with
  phase declarations `gas, liquid, ice` and `gas, ice, liquid`, fully
  permuting all coupled indices and arrays.
- [x] Compare final results by physical component identity, not storage index:
  expected
  topology is `gas + ice`, with liquid inactive. Check moles, gas fractions,
  conservation, complementarity, and the independent scalar candidate oracle.
- [x] Do not require identical transition histories or transition counts. Store
  histories as characterization and assert only fixed-point invariance.
- [x] Add initial-history cases `gas only`, `gas + liquid`, and `gas + ice`
  using positive physical phase inventories, without synthetic invalid
  accepted states. All six history/order routes converge to `gas + ice`.
- [x] Record a clear verdict: `PhaseCandidateOrderInvariant` or
  `PhaseCandidateOrderDependent`, with numerical max deltas and full failure
  evidence. The current verdict is `PhaseCandidateOrderInvariant`; no
  production fix was needed.

### Frozen/release evidence

- [x] Keep the IAPWS gauge, strict canonical TPD/production assertions, and
  catalog feasibility characterization as separate evidence channels.
- [x] Add the focused phase permutation/unpermutation regression and physical
  identity matching. A release-only story is not required for this tiny
  deterministic fixture; the debug story is strict and reproducible.
- [x] Record exact identifiers, provenance, domains, selected conditions,
  candidate TPDs, oracle values, transition histories, and the command in
  `STORY_TESTS.md`. Release characterization was completed and recorded; the
  release result agrees with the debug result.
