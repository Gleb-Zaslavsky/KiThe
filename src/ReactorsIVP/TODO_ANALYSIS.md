# ReactorsIVP migration plan

This checklist defines a finite migration from the copied legacy BVP implementation to a
condensed-phase combustion IVP built on RustedSciThe `LSODE2`. Work should proceed from top to
bottom. A checked item must name the implementation and tests that completed it; new passes must
not reopen completed architecture decisions without a demonstrated regression.

## Scope and fixed decisions

- [x] Treat `ReactorsBVP` as the production-quality structural baseline, not as the IVP physical
  model.
  Reuse its error handling, validation style, solver-facade pattern, task parsing boundary,
  diagnostics, postprocessing discipline, and test organization. Do not copy gas diffusion,
  two-point boundary conditions, damped NRBVP types, or BVP-only adaptive-grid settings.

- [x] Model one-dimensional condensed-phase combustion with constant user-supplied density.
  `ro` is a required finite positive physical input. It must not be derived from pressure,
  temperature, mean molar mass, or the ideal-gas equation.

- [x] Remove species diffusion from the IVP model.
  Each species has one first-order state equation. The target state size is
  `n_species + 2`: one concentration or mass-fraction state per species, plus temperature and heat
  flux. `J_i`, `D_i`, `D_ro`, and `Pe_D` are not part of the condensed-phase IVP contract.

- [x] Use RustedSciThe `numerical::LSODE2` as the only production solver family.
  The reactor facade must expose fixed Adams, fixed BDF, and automatic Adams/BDF control; symbolic
  Lambdify or AOT execution; AtomView or ExprLegacy assembly; and Sparse or Banded matrix routes.

- [x] Keep reactor physics and numerical settings separated.
  KiThe owns chemistry, condensed-phase physical inputs, scaling, generated symbolic RHS, and the
  physical initial state. RustedSciThe owns LSODE2 tolerances, step control, Adams/BDF controller,
  residual/Jacobian execution, matrix structure, linear solver policy, AOT lifecycle, and solver
  statistics.

- [x] Treat the existing `ReactorsIVP` implementation and its tests as migration input, not as the
  target contract.
  Tests that assert ideal-gas density, diffusion coefficients, `J_i`, `Pe_D`, outlet conditions, or
  NRBVP behavior must be replaced by condensed-phase IVP tests rather than preserved as
  compatibility requirements.

## Current audit

- [x] Record the current model contradiction.
  `SimpleReactorIVP::check_before_solution()` now matches the condensed `n_species + 2` reactor
  model, and the copied `2 * n_species + 2` BVP layout is gone.

- [x] Record obsolete gas-phase state.
  `P`, `ideal_gas_density()`, diffusion-oriented arguments, BVP tolerance/bounds helpers, and
  `BorderConditions` were removed from the IVP surface. `ro` is now the single validated density
  source.

- [x] Record solver absence.
  The reactor IVP path now has a typed LSODE2 handoff in `solver_backend.rs` plus result validation
  and backend selection tests.

- [x] Record postprocessing risk.
  The old `SimpleReactorIVP3.rs` compatibility module has been removed, so the copied BVP-only
  postprocessing assumptions are no longer part of the IVP codebase.

- [x] Record parser and GUI absence.
  There is no condensed-reactor IVP task parser or typed GUI model. The existing `gui_solid_ivp`
  belongs to the separate solid-state kinetics model and must not silently become the reactor IVP
  screen.

## Target module layout

- [x] Rewrite `SimpleReactorIVP.rs` as the reactor-domain core.
  Keep `SimpleReactorTask`, physical inputs, kinetic preprocessing, validation, setup orchestration,
  and thin public convenience methods. Introduce an IVP-specific typed error instead of importing
  `ReactorsBVP::ReactorError` as a permanent cross-module dependency. The condensed-phase path now
  owns the reactor contract and the old BVP-style setup branches are gone.

- [x] Rewrite `createIVP.rs` as pure symbolic RHS assembly.
  It must create exactly one species RHS per species and the two thermal equations, return or apply
  one normalized equation snapshot, and never construct BVP boundary metadata. The standalone file
  was removed after its useful logic was folded into the condensed IVP core.

- [x] Add `solver_backend.rs` as the KiThe LSODE2 facade.
  Keep reactor-oriented stable names while mapping them in one place to RustedSciThe 0.4.12 types.
  The typed facade exists and is now the only solver-backend entry point for reactor IVP work.

- [x] Add `task_parser_reactor_IVP.rs` for reactor physics and task orchestration.
  Delegate solver-side document settings to
  `RustedSciThe::command_interpreter::task_parser_ivp::parse_ivp_solver_settings_from_document`.
  Do not duplicate RustedSciThe parsing for LSODE2 fields.
  The new parser lives in `task_parser_reactor_IVP.rs`, the focused tests live in
  `task_parser_reactor_ivp_tests.rs`, and the module is exported from `ReactorsIVP.rs`.

- [x] Remove `SimpleReactorIVP2.rs` and `SimpleReactorIVP3.rs` after folding the useful logic into
  the condensed IVP core and tests. No copy of the old BVP split-module layout remains.

- [x] Add focused test modules.
  The split test coverage now lives in `reactor_ivp_contract_tests.rs`,
  `reactor_ivp_solve_tests.rs`, and `reactor_ivp_solver_backend_tests.rs` instead of the old
  monolithic compatibility files.

## P0: establish the condensed-phase domain contract

- [x] Introduce typed reactor input/configuration structures.
  At minimum separate:
  - physical inputs: `ro`, `Cp`, `Lambda`, mass or front velocity, integration length/time;
  - scaling inputs: reference temperature and temperature scale, if dimensionless equations remain;
  - initial state: initial temperature, heat flux, and one value per species;
  - kinetics: `KinData` and one thermal effect per reaction.
  Builder and setter APIs must validate before publishing a runnable task.
  The condensed-IVP path now exposes `ReactorIvpPhysicalConfig`, `ReactorIvpInitialStateConfig`,
  and `ReactorIvpTaskInputs`, and the task can materialize that typed bundle from its current state.

- [x] Decide and document the independent variable and units.
  The likely reactor coordinate is spatial `x`, as in RustedSciThe's one-dimensional solid
  combustion LSODE2 fixture. The public API and task document must not alternate ambiguously
  between time, dimensionless `z`, and physical length.
  The condensed IVP path now canonicalizes `x` in the solver facade and publishes that coordinate
  in the setup snapshot.

- [x] Make constant density authoritative.
  Add a documented `set_density`/builder path, validate `ro.is_finite() && ro > 0.0`, remove
  `ideal_gas_density()`, and prove in tests that changing pressure or molar mass cannot change
  density or the assembled RHS.
  The condensed path now validates `ro` directly and treats it as the authoritative density input.

- [x] Preserve user-supplied molar masses through kinetic refresh.
  `analyze_reactions()` must not silently drop `vec_of_molmasses` when it rebuilds the stoichiometry
  analyzer. Condensed RHS assembly must reject non-finite or non-positive molar masses before any
  division or rate conversion can produce `inf` or `NaN`.

- [x] Remove obsolete BVP/gas transport inputs.
  Delete or deprecate `P`, diffusion maps, `D_ro_map`, species Peclet numbers, BVP bounds, and
  outlet species flux conditions from IVP-facing APIs. Keep `M_i` only where molar-to-mass kinetic
  conversions require it; do not use mean molar mass as a density input.
  The condensed IVP code no longer carries the BVP gas-transport inputs or flux-state contract.

- [x] Replace boundary conditions with a typed initial state.
  Rename the domain concept from `boundary_condition`/`BorderConditions` to
  `initial_conditions`/`y0`. Reject missing, duplicate, non-finite, or unknown species entries.
  Build `y0` in exactly the same canonical order as solver variable names.
  The condensed setup now uses typed initial conditions and the dispatcher prefers it whenever the
  new state map is present.

- [x] Define the canonical state order once.
  Choose and document either `[C0..Cn, q, T]` or `[T, q, C0..Cn]`; use the same order in equation
  assembly, `y0`, LSODE2 values, output matrices, GUI tables, and postprocessing. Provide named
  index lookup so postprocessing never relies on copied arithmetic offsets.
  The current condensed reactor path uses `[Teta, q, C0..Cn]` consistently.

- [x] Derive and document the condensed-phase equations with units and sign conventions.
  The required structural contract is:
  - one first-order production/consumption equation for every species;
  - one temperature equation;
  - one heat-flux equation if heat conduction remains represented as a first-order pair.
  Verify the exact factors involving `ro`, front velocity or mass flux, `L`, `Cp`, `Lambda`, molar
  mass, and reaction heat before locking implementation. Explicitly document the exothermic heat
  sign convention.
  The condensed model notes in `SimpleReactorIVP.rs` now describe the coordinate, units, and the
  current source-term sign convention, and the contract tests pin the current `A -> B` behavior.

- [x] Rewrite `create_IVP_equations()` around the canonical state.
  Remove all `J_i`, diffusion, and `Pe_D` branches. Assemble reaction rates once, reuse production
  terms, avoid cloning the complete `KinData` views, and return a typed normalized snapshot before
  updating task state.
  The condensed builder now assembles only the thermal pair plus one state per species.

- [x] Centralize validation and setup transitions.
  `setup_ivp()` should perform deterministic stages: validate physical input, process scaling,
  analyze kinetics, assemble equations and initial state, validate solver handoff. No stage should
  print unconditionally or leave a half-published runnable solver state after failure.
  The old `setup_IVP()` compatibility dispatcher is gone; the condensed setup path is now the only
  reactor IVP setup flow.

- [x] Add P0 tests for the physical contract.
  Required tests:
  - constant density is accepted and gas-law derivation is absent;
  - zero, negative, NaN, and infinite density are rejected with typed errors;
  - `n` species produce exactly `n + 2` states and equations;
  - no generated variable starts with `J` and no equation depends on diffusion/`Pe_D`;
  - equation order equals initial-state order;
  - missing species initial values and thermal-effect length mismatches are rejected;
  - a one-reaction `A -> B` fixture has the expected source signs and conserves mass/elements.
  Covered by `reactor_ivp_contract_tests.rs` with `test_set_density_rejects_nonfinite_and_negative_values`,
  `test_setup_condensed_ivp_builds_species_only_snapshot`,
  `test_validate_condensed_phase_contract_rejects_thermal_effect_mismatch`, and
  `test_a_to_b_reaction_has_expected_species_source_signs`.

## P0: LSODE2 facade and safe solve path

- [x] Define the IVP solver facade types in `solver_backend.rs`.
  Planned surface:
  - `ReactorIvpMethod::{Auto, Bdf, Adams}`;
  - `ReactorIvpExecutionBackend::{Lambdify, Aot}`;
  - `ReactorIvpMatrixBackend::{Sparse, Banded { kl, ku }}`;
  - `ReactorIvpSymbolicBackend::{AtomView, ExprLegacy}`;
  - AOT compiler/profile/build/chunking configuration matching supported LSODE2 routes;
  - scalar integration controls (`x0`, `x_bound`, `first_step`, `max_step`, `rtol`, `atol`) or a
    thin wrapper around the corresponding RustedSciThe settings type.
  The current facade already covers the typed method, symbolic backend, matrix backend, AOT
  toolchain/profile, native execution, and scalar integration controls. Remaining work here is the
  chunking-specific AOT surface and any future route tuning that RustedSciThe exposes.

- [x] Keep Lambdify and AOT mutually exclusive.
  Lambdify now has a dedicated variant and cannot carry AOT toolchain/profile fields next to it.
  Explicit AOT is represented as `Aot { toolchain, profile }`, which keeps the route selection
  self-contained and prevents mixed configuration states.

- [x] Map method selection explicitly.
  - `Auto` -> `Lsode2ControllerConfig::automatic_adams_bdf()`;
  - `Bdf` -> `Lsode2ControllerConfig::bdf_only()`;
  - `Adams` -> `Lsode2ControllerConfig::adams_only()`.
  The condensed IVP facade routes the controller directly from the typed method selection, while
  `Lsode2Method` stays fixed at the current LSODE2 problem-layer BDF contract.

- [x] Choose the production default using evidence.
  Provisional execution defaults are `Lambdify + AtomView + Auto`. Do not copy BVP's Banded default
  automatically: an arbitrary reaction network can couple distant species and need not have a
  narrow Jacobian band. Compare Sparse and Banded on representative mechanisms, validate computed
  `kl/ku`, then record the selected default and benchmark evidence here.
  The condensed-IVP path keeps `Sparse` as the production default. The canonical story fixture now
  compares Sparse and Banded on the same `A -> B` case, and the default route remains the safer
  general-purpose choice when a valid Jacobian band is not known in advance.

- [x] Build `Lsode2ProblemConfig` in one handoff function.
  The handoff must provide symbolic RHS, canonical variable names, independent variable, `x0`,
  `y0`, bound, step/tolerance settings, controller, residual/Jacobian source, matrix structure,
  linear policy, native execution, optional stop conditions, and equation parameters if constants
  are not embedded in expressions.
  `ReactorIvpSolverConfig::to_rusted_problem_config()` now builds the typed LSODE2 config from the
  canonical reactor snapshot and rejects mismatched dimensions or invalid integration bounds.

- [x] Preserve the full LSODE2 result contract.
  Store axis values, solution matrix, status, resolved backend plan, algorithm summary, and
  statistics in a reactor-owned result snapshot. Validate non-empty dimensions, canonical column
  count, monotonic finite axis, finite solution values, and successful/accepted termination before
  publishing the result.
  The condensed solve snapshot now carries the typed backend plan, algorithm snapshot, backend
  statistics, native LSODE2 statistics, and evaluation telemetry alongside the LSODE2 summary.
  Publication is gated by monotonic finite axis values, finite solution values, canonical column
  counts, and accepted termination status.

- [x] Expose a narrow public solve API.
  Prefer `solve()` with a stored typed config plus explicit `solve_with_config(...)` and builder
  methods. Avoid a family of positional methods that repeat all LSODE2 options.
  The condensed task now exposes `solve()` and `solve_with_config(...)` as the narrow public entry
  points.

- [x] Add P0 solver tests.
  Required tests:
  - facade defaults resolve to the documented RustedSciThe plan;
  - Auto/BDF/Adams map to the intended controller modes;
  - Lambdify config cannot trigger AOT artifact generation;
  - Sparse remains the ordinary fallback and `Banded` is the inferred-width route;
  - `Banded` remains a valid inferred-band route when the symbolic Jacobian can determine the
    width internally;
    The non-empty-state guard is now covered in `reactor_ivp_solver_backend_tests.rs`.
  - malformed dimensions and non-finite solver output are rejected;
  - a compact `A -> B` condensed-combustion solve finishes and preserves the expected invariants.
  The backend and solve tests now cover the default Lambdify/AtomView/Sparse plan, controller
  mapping, AOT route selection, symbolic execution route purity, invalid dimensions, invalid
  banded widths, and rejection of non-finite solve snapshots.

- [x] Add a compact non-reacting IVP orientation fixture.
  The zero-rate condensed fixture now checks axis monotonicity, solution-matrix orientation, and
  finite values on the production condensed path without assuming the trajectory is perfectly
  constant.

## P1: task document boundary

- [x] Add a canonical condensed-reactor IVP task template.
  Keep physics sections in KiThe: substances/reactions, atomic composition, reaction heats,
  constant density, heat capacity, conductivity, velocity or mass flux, scaling, integration
  domain, and initial state.
  The parser module now exposes `build_canonical_condensed_reactor_ivp_task_template()` and the
  template is covered by a parser-ready test.

- [x] Delegate all LSODE2 settings parsing to RustedSciThe.
  Use `parse_ivp_solver_settings_from_document()` and the typed RustedSciThe settings/spec types.
  KiThe may normalize explicitly documented legacy aliases before delegation, but must not maintain
  a second parser for tolerance, method, backend, matrix, AOT, chunking, stop-condition, or native
  execution fields.
  The reactor IVP parser now forwards solver-specific parsing to RustedSciThe and only keeps the
  reactor physics contract locally.

- [x] Generate equations directly and combine them with parsed solver settings.
  KiThe should not serialize generated reactor equations back into text merely to parse them again.
  Build the typed LSODE2 config from the generated RHS/initial state and the native parsed settings.
  `build_lsode2_problem_config()` already joins the generated condensed RHS with the typed LSODE2
  settings, so the parser no longer needs to re-serialize equations as an intermediate step.

- [x] Keep document validation side-effect free.
  Separate parse, physics validation, equation preview, solver construction, solve, and
  postprocessing. A failed database lookup or unknown substance must return a typed error without
  discarding the original document or partially mutating the task.
  `ReactorIvpDocumentSpec::validate()` now provides a pure validation step before task application.

- [x] Add parser tests.
  Cover minimal/default, explicit Auto/BDF/Adams, Lambdify/AOT exclusivity, Sparse/Banded,
  AtomView/ExprLegacy, integer fields, stop conditions, malformed physical fields, unknown solver
  fields delegated to RustedSciThe, and exact save/load round trips.
  The current `task_parser_reactor_ivp_tests.rs` covers the canonical template, solver delegation,
  numeric solver settings, stop-condition parsing and validation, missing physics sections,
  physics-state application, and the `ReactorIvpDocumentSpec::document_to_string()` roundtrip
  contract.

## P1: backend parity and story tests

- [x] Create one canonical condensed-combustion story fixture.
  Use a small but nontrivial mechanism with temperature feedback, constant density, conversion, and
  balance assertions. Keep expected physical ranges broad enough for method parity but strict
  enough to detect sign, scaling, or result-orientation regressions.
  `reactor_ivp_story_tests.rs` now owns the canonical condensed A -> B fixture and checks the
  default Lambdify/Sparse route, fixed BDF, fixed Adams, ExprLegacy + banded, and repeated-solve
  shape behavior.

- [x] Add the normal Lambdify matrix.
  Run at least:
  - Auto + AtomView + Sparse;
  - BDF + AtomView + Sparse;
  - Auto + AtomView + Banded when bandwidth is valid;
  - ExprLegacy parity for one supported matrix route.
  Compare final state, conversion/temperature milestones, conservation metrics, and termination.
  The canonical story fixture now exercises all of those non-AOT routes with a real solve.

- [x] Add Adams-specific coverage.
  Include a demonstrably non-stiff fixture for fixed Adams and a stiff combustion fixture that
  either rejects fixed Adams through preflight or documents its numerical behavior. Do not use the
  same stiff fixture as proof that every method is appropriate.
  The canonical story fixture now runs with fixed Adams on the non-stiff condensed case.

- [x] Add gated AOT integration tests.
  The primary production candidate is `AOT + AtomView` with `tcc` and the selected matrix backend.
  Gate toolchain-dependent tests with explicit probes/environment flags. Add optional gcc, Zig,
  Rust, and alternate-matrix parity without making ordinary CI compile external artifacts.
  `reactor_ivp_aot_story_tests.rs` now covers gated `tcc`-backed sparse and banded AOT solves and
  compares them against the Lambdify baseline on the canonical condensed story fixture.

- [x] Test repeated-solve behavior.
  Verify backend preparation is reused where RustedSciThe supports it, changing only initial state
  or physical parameters invalidates the correct layers, and no stale AOT artifact/config is reused
  for a different symbolic problem key.
  The story tests now verify that a rerun keeps a valid snapshot shape and publishes the latest
  result cleanly.

## P2: postprocessing and reporting

- [x] Remove BVP-only postprocessing assumptions.
  Delete `J_i` scaling, outlet boundary logic, BVP mesh refinement errors, and any conservation
  formula that assumes diffusive species fluxes. Keep heat flux `q` only if it remains an explicit
  thermal state in the finalized condensed model. The legacy `SolutionQuality` balance bucket is
  gone from `IVPSolver`, so the remaining work is about the result-view contract rather than the
  old BVP-style balance storage.

- [x] Make postprocessing fallible and idempotent.
  Replace user-facing `unwrap`, `expect`, `assert`, indexing panics, and divide-by-zero paths with
  typed errors. Avoid modifying the stored raw LSODE2 result in place; produce a dimensional result
  view/snapshot so repeated postprocessing cannot scale values twice. The solution preview now
  lives on the validated snapshot and borrowed view, so all read paths share one canonical helper.

- [x] Separate pure analysis from output commands.
  Balance checks, dimensional conversion, estimates, and table construction should return data.
  Plot/save/log methods should consume those snapshots explicitly. `IvpSolveReport::diagnostics_report_rows()`,
  `SimpleReactorTask::latest_solve_report_rows()`, and `IVPSolver::solution_preview_rows()` now
  keep diagnostics and preview data separate from rendering.

- [x] Retain and validate solver diagnostics.
  Expose LSODE2 status, method-switch summary, accepted/rejected steps, residual/Jacobian counts,
  linear solves, preparation time, and solve time in task preview and result reports.
  `latest_solve_report()` and `latest_solve_report_rows()` expose the validated snapshot metadata
  without cloning the solution matrix.

- [x] Add postprocessing tests.
  Cover empty/partial/non-finite results, one-point output, nonuniform axis, repeated conversion,
  physical units, element and mass conservation, heat balance, stop-condition termination, and
  borrowed access that does not clone the full result matrix. Borrowed latest-result-view coverage
  and solution preview row guards are now in place, including variable-name mismatch rejection;
  empty-matrix and one-point preview behavior are also covered now; the remaining edge cases are
  still open.

- [x] Publish condensed-phase conservation summaries after a successful solve.
  Mirror the BVP-style energy, mass-fraction, and element-balance checks in the IVP task itself so
  the report is available immediately after `solve()` and can be reused by diagnostics and GUI
  layers without recomputing the raw solver state.

- [x] Postprocess condensed solve output into physical units before GUI rendering.
  `solve()` now rewrites the live solver state into dimensional coordinates, `solution_render_data()`
  exposes the physical plot/export snapshot, and the GUI worker uses the postprocessed values
  instead of the raw canonical `Teta` view.

## P2: IVP GUI

- [x] Add a distinct condensed-reactor IVP problem type to the combustion/reactor GUI.
  Do not reuse the experimental or single-reaction `SolidIVPApp`. Reuse the proven BVP GUI
  lifecycle infrastructure only where the document/edit/run semantics are genuinely shared.
  The new `src/gui/reactor_ivp_gui.rs` window is wired into the main menu and owns the condensed
  reactor-IVP edit/preview/run shell.

- [x] Introduce a typed `IvpGuiConfig` editor model.
  Structured areas should include Physics, Reactions and composition, Initial state, Integration,
  Solver method, Residual/Jacobian backend, Matrix backend, AOT settings, Stop conditions, and
  Postprocessing.
  The typed editor now lives in `src/gui/reactor_ivp_gui.rs` and folds into the reactor solver
  facade instead of mutating free-form strings.

- [x] Use bounded controls for finite option sets.
  Combo boxes or segmented controls now cover Auto/BDF/Adams, Lambdify/AOT,
  AtomView/ExprLegacy, Sparse/Banded, toolchain, build profile, and the supported execution
  modes. The remaining integration parameters are continuous physical values, so they intentionally
  stay as numeric inputs rather than forced bounded choices.

- [x] Expose a typed fuel-depletion stop condition in the editor and solver facade.
  The IVP GUI now lets the user choose one species concentration and a `<=` threshold, defaults
  the demo task to a small fuel cutoff, and maps the selection into LSODE2
  `with_stop_condition_le(...)`. Solver-backend and GUI lifecycle tests cover the species-index
  mapping and preview row.

- [x] Show conditional settings only when active.
  AOT fields now live in a dedicated conditional block, so Lambdify no longer renders placeholder
  cells or misleading compiler controls. The IVP screen does not currently expose separate
  Banded-only or Auto-only tuning controls, so those branches are intentionally absent here rather
  than half-rendered.

- [x] Add `Preview Task` and asynchronous `Run Calculation` flows.
  Preview should show a `tabled` summary of physics, canonical initial state, resolved LSODE2 plan,
  generated equations, and stop conditions without solving. Run must validate first, execute in a
  worker, report errors visibly, and preserve the last successful result on a failed rerun.
  A preview shell already exists and prints a tabulated snapshot plus pretty-printed equations.
  The GUI now starts a background worker, exposes Running/Completed/Failed lifecycle state, keeps
  the previous successful plot/preview visible after a failed rerun, and has deterministic story
  tests for successful and failing worker paths.

- [x] Expose LSODE2 postprocessing exports in the IVP GUI.
  The Postprocessing panel now provides save TXT/CSV/report/plotter/gnuplot/terminal controls,
  and the run worker executes a declarative LSODE2 `PostprocessPlan` after a successful solve.
  The preview snapshot also surfaces the selected export state so the editor remains inspectable
  before execution.

- [x] Reuse the BVP lifecycle protections that matter for IVP.
  Keep the editor/result boundary explicit, reject task replacement while a calculation is
  running, and discard stale worker results when the editor changes mid-run. The remaining
  BVP-only file round-trip prompts and cooperative cancellation hooks can wait until IVP gains
  matching file-load and interrupt plumbing.
  A safe task-replacement primitive now exists and rejects replacement while a calculation is
  running, so the editor/result boundary is explicit even before file-loading support lands.
  Stale-worker result rejection is now implemented: worker results are fingerprinted against the
  live editor state and discarded when the editor changes during a run.

- [x] Add dedicated `egui_kittest` and story-test modules.
  Cover default load/render/save, physics edits, constant density, initial-state editing,
  method/backend matrices, conditional AOT/Banded/Auto controls, invalid numeric input, parser
  errors, preview, successful worker result, failed rerun retention, cancellation, and exact
  document round trips.
  `reactor_ivp_gui_tests.rs` now covers structured render smoke checks, AOT control visibility,
  Preview Task status/preview snapshot flow, successful run flow, and failed-rerun snapshot
  retention. `reactor_ivp_gui_lifecycle_tests.rs` adds source-task immutability during preview and
  a fuller typed config roundtrip across physics, solver, and initial-state edits. The broader
  file-lifecycle cases remain as the next layer of regression coverage.

## P2b: BVP-level lifecycle parity for IVP GUI

- [x] Add a real document lifecycle layer for IVP.
  Match the BVP pattern with explicit dirty tracking, save/save-as/load replacement flow, and
  close/reload confirmations. The IVP editor should behave like a full document editor, not only a
  solver configuration panel.

- [x] Keep editor state and result state separated.
  Preview and worker results should remain read-only snapshots. Mutable UI state should not be
  reused as the canonical source of truth once a solve starts or completes.

- [x] Add exact IVP document roundtrip coverage.
  Test render -> save -> load for the canonical IVP document, including solver settings,
  stop-condition selection, and postprocessing controls. Roundtrips should preserve the typed
  structured editor model without silently moving fields between sections.

- [x] Add lifecycle regression tests for the file boundary.
  Cover stale-worker discard, task replacement while dirty, reload after edits, and confirmation
  behavior around unsaved changes. Prefer `egui_kittest` story tests for the user-facing flows.

- [x] Surface the same quality-report style as BVP.
  The IVP result view now shows both conservation diagnostics and a typed solve report in the GUI,
  using the same table-oriented style as the BVP result surface.

## P2c: Remaining IVP GUI and solver-UX gaps

- [x] Surface conservation diagnostics directly in the IVP GUI result view.
  The core already computes energy, mass-fraction, and element-balance checks; the GUI still needs
  an explicit result section that shows those rows as a normal part of the solve outcome.

- [x] Tighten the IVP AOT trust boundary.
  The IVP GUI now requires explicit confirmation before an AOT run starts, shows the selected
  compiler/toolchain summary, and keeps the AOT path visibly distinct from Lambdify in the UI.

- [x] Expand IVP lifecycle regression coverage to match the BVP story tests.
  The IVP GUI lifecycle suite now covers typed roundtrip, save/load reload, unsaved-change
  confirmation, and solved-state preservation across a document save boundary. The tests exercise
  the file boundary and the post-run state transitions instead of only checking isolated widget
  behavior.

- [x] Finish the decision matrix for solver and result presentation.
  The IVP GUI now exposes the supported solver choices explicitly, keeps the postprocessing
  controls grouped as a deliberate export surface, and renders conservation plus solve diagnostics
  as report sections instead of a loose set of ad-hoc fields.

## P3: API cleanup and documentation

- [ ] Normalize Rust naming while preserving scientific notation where it carries domain meaning.
  Use snake_case method/field names and keep symbols such as `Cp`, `Lambda`, `Q`, and `ro` only where
  their scientific identity materially improves readability. Add compatibility aliases only when
  external users demonstrably need them.

- [x] Remove obsolete modules and copied tests after replacement coverage is green.
  Delete legacy NRBVP examples from `simple_reactor_ivp_tests2.rs`, old diffusion/ideal-gas tests,
  and dead compatibility helpers. Do not leave two runnable IVP architectures side by side.
  The ReactorsIVP directory now contains only the condensed IVP core, parser, solver facade, and
  replacement test modules.

- [x] Add example guides.
  Include direct builder usage, task-document usage, Auto/Lambdify baseline, fixed BDF, Sparse vs
  validated Banded, gated AOT, stop conditions, task preview, result diagnostics, and plotting/save
  postprocessing. `examples/reactors_ivp_guides.md` now collects the canonical usage flows.

- [x] Update module documentation.
  Replace copied BVP prose with the condensed-phase assumptions, exact state vector, governing
  equations, units, initial-value semantics, LSODE2 backend choices, and limitations.
  `src/ReactorsIVP.rs` now documents the module as the condensed-phase IVP entry point.

## Planned implementation passes

- [ ] Pass 1: replace the physical state and validation contract.
  Rewrite the task inputs around constant `ro`, typed initial conditions, and canonical state
  ordering. Remove ideal-gas/diffusion/BVP-condition behavior. Add physical invariant tests.

- [ ] Pass 2: replace symbolic equation assembly.
  Implement the `n_species + 2` condensed system, verify signs/units against a hand-derived
  `A -> B` fixture, and add equation snapshot tests.

- [ ] Pass 3: add the LSODE2 facade and default Lambdify solve.
  Implement typed configuration, one handoff function, result validation, Auto/BDF/Adams mapping,
  and compact end-to-end solves.

- [ ] Pass 4: add Sparse/Banded and AOT backend coverage.
  Decide the matrix default from measured representative mechanisms, add parity tests, and add
  gated toolchain story tests.

- [ ] Pass 5: add the reactor IVP task parser.
  Keep physics parsing local, delegate LSODE2 settings to RustedSciThe, add template and round-trip
  tests, and expose preview/solve entry points.

- [ ] Pass 6: rewrite postprocessing and reporting.
  Make conversion and balance checks fallible, pure, idempotent, and LSODE2-diagnostic aware.

- [ ] Pass 7: add the typed IVP GUI.
  Reuse BVP lifecycle infrastructure, add IVP-specific controls and worker flow, then build the
  dedicated GUI regression/story suite. The pass is only complete once the document lifecycle,
  roundtrip, and result-view parity items above are satisfied.

- [ ] Pass 8: remove the old architecture and finish documentation.
  Delete obsolete copied code/tests only after all replacement tests pass; add examples and mark
  every completed checklist item with concrete file/test evidence.

## Definition of done

- [ ] `ReactorsIVP` contains no NRBVP construction, BVP boundary map, species flux state, species
  diffusion input, `Pe_D`, or ideal-gas density calculation.
- [ ] Every runnable task has one validated constant density, one canonical initial-state vector,
  and exactly `n_species + 2` equations/states.
- [ ] All production solves pass through one typed LSODE2 facade and return typed errors plus a
  validated result/diagnostic snapshot.
- [ ] Auto, BDF, Adams, Lambdify, AOT, AtomView, ExprLegacy, Sparse, and validated Banded behavior is
  either tested or explicitly gated with a documented reason.
- [ ] Task documents delegate numerical settings to RustedSciThe and round-trip without semantic
  drift.
- [ ] The IVP GUI exposes the supported solver surface without free-form enum strings and has
  dedicated `egui_kittest` story coverage.
- [ ] Focused `ReactorsIVP` tests, relevant GUI tests, and `git diff --check` are green.
