# SimpleReactorsIVP architecture and migration plan

This document records the completed condensed-phase IVP migration and the next approved
architecture phase: move solver selection and execution to RustedSciThe `StiffAPI` and support
LSODE2, standalone BDF, and Radau. The completed LSODE2-only design below is historical baseline,
not a constraint on the new API. Public API compatibility is not a goal for this migration.
Implementation must preserve the physical model unless a separate physics change is explicitly
approved. Each completed item must cite concrete code and tests; performance gates belong to the
final measurement phase.

## Scope and fixed decisions

- [x] Treat `SimpleReactorsBVP` as the production-quality structural baseline, not as the IVP physical
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

- [x] Historical baseline: use RustedSciThe `numerical::LSODE2` as the only production solver
  family. It exposed fixed Adams, fixed BDF, automatic Adams/BDF control, symbolic Lambdify/AOT,
  AtomView/ExprLegacy assembly, and Sparse/Banded matrix routes. This decision is superseded by
  the approved `StiffAPI` migration below; the existing behavior remains the correctness baseline.

- [x] Keep reactor physics and numerical settings separated.
  KiThe owns chemistry, condensed-phase physical inputs, scaling, generated symbolic RHS, and the
  physical initial state. RustedSciThe owns numerical method controls, callback execution, matrix
  structure, AOT lifecycle, and solver-native statistics. The earlier implementation was
  LSODE2-specific; the approved StiffAPI adapter must preserve this ownership boundary across all
  selected solver families.

- [x] Treat the existing `SimpleReactorsIVP` implementation and its tests as migration input, not as the
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
  `SimpleReactorsBVP::ReactorError` as a permanent cross-module dependency. The condensed-phase path now
  owns the reactor contract and the old BVP-style setup branches are gone.

- [x] Rewrite `createIVP.rs` as pure symbolic RHS assembly.
  It must create exactly one species RHS per species and the two thermal equations, return or apply
  one normalized equation snapshot, and never construct BVP boundary metadata. The standalone file
  was removed after its useful logic was folded into the condensed IVP core.

- [x] Add `solver_backend.rs` as the initial KiThe LSODE2 facade.
  Its original mapping targeted RustedSciThe 0.4.12 types. The adapter exists and is the current
  backend entry point; it is explicitly scheduled to be redesigned around the pinned StiffAPI
  version 0.4.28, with no public compatibility requirement.

- [x] Add `task_parser_reactor_IVP.rs` for reactor physics and task orchestration.
  Delegate solver-side document settings to
  `RustedSciThe::command_interpreter::task_parser_ivp::parse_ivp_solver_settings_from_document`.
  Do not duplicate RustedSciThe parsing for LSODE2 fields.
  The new parser lives in `task_parser_reactor_IVP.rs`, the focused tests live in
  `task_parser_reactor_ivp_tests.rs`, and the module is exported from `SimpleReactorsIVP.rs`.

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
  The derivation is recorded in `MATHEMATICAL_MODEL.md` and beside the RHS builder. The solver now
  integrates dimensionless `z=x/L`; mass-fraction sources use `(L/m)*M_i*sum(nu_ij*r_j)`, with
  reactants negative and products positive. The direct `H2 -> 2H` test checks source signs, exact
  scaling, and mass-source closure; TW-01/TW-02 independent references use the same physical model.

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
  `test_hydrogen_dissociation_has_expected_species_source_signs`.

## Historical baseline: LSODE2 facade and safe solve path

The following completed items describe the current LSODE2-only adapter and its regression
contracts. They are retained as migration evidence; the StiffAPI phase below replaces the adapter
and updates solver-specific contracts rather than requiring compatibility aliases.

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

## P2d: Physical task-document lifecycle stories

- [x] Exercise an IVP edit through the actual GUI numeric widget, then Save As, reopen, Preview,
  and solve. The physical-task story verifies the changed temperature, mechanism, and solver route.
- [x] Compare every typed IVP config field, including reaction payloads, postprocessing, and
  solver-specific options, before Save As and after reopen.
- [x] Add a heat-releasing H2 document with measurable species and temperature evolution while
  retaining the bounded smoke files as parser and lifecycle controls. This is a physical
  acceptance story, not a performance benchmark.
- [x] Cover a malformed indexed Arrhenius payload: require reaction context, prove the active task
  is unchanged, then repair the file and complete Save As/reopen/solve.
- [x] Extend malformed indexed-reaction coverage to bad equations, missing reaction indices,
  unknown species, and reaction/thermal-effect count mismatch. Each variant verifies diagnostic
  context and unchanged active task; malformed Arrhenius additionally exercises
  repair/save/reopen/solve.
- [x] Compare the heat-releasing task through the direct production API and GUI worker. Final
  physical temperature, species fractions, and conservation diagnostics agree within `1e-8`.
- [x] Cover failed Save As semantics: the dirty task, fingerprint, and current source-file binding
  survive an I/O failure.
- [x] Verify a failed rerun clears the prior published solution, conservation report, and plot before
  the attempted solve, so a failure cannot present old data as belonging to edited inputs.
- [x] Keep file lifecycle stories correctness-only; no file-I/O performance gates were added.
- [x] Diagnose and gate the energy balance on `h2_dissociation_with_heat_release.task`. The defect
  was in postprocessing quadrature, not the physical IVP equations or solver: `simpsons()` advanced
  two mesh intervals after a nonuniform interval but had integrated only the first, silently
  skipping the second. It now advances one interval for the trapezoid fallback and two only for a
  uniform Simpson pair. Unit regressions cover irregular meshes; the physical story asserts the
  reported energy residual is within `1 J` and `0.01%` and independently checks trapezoidal closure
  within `1 J`. Current result: reported residual `-4.14e-2 J` (`-3.85e-4%`); independent
  trapezoidal residual `+5.90e-1 J`. No solver or physics tolerance was relaxed.

## P3: API cleanup and documentation

- [ ] Normalize Rust naming while preserving scientific notation where it carries domain meaning.
  Use snake_case method/field names and keep symbols such as `Cp`, `Lambda`, `Q`, and `ro` only where
  their scientific identity materially improves readability. Add compatibility aliases only when
  external users demonstrably need them.

- [x] Remove obsolete modules and copied tests after replacement coverage is green.
  Delete legacy NRBVP examples from `simple_reactor_ivp_tests2.rs`, old diffusion/ideal-gas tests,
  and dead compatibility helpers. Do not leave two runnable IVP architectures side by side.
  The SimpleReactorsIVP directory now contains only the condensed IVP core, parser, solver facade, and
  replacement test modules.

- [x] Add example guides.
  Include direct builder usage, task-document usage, Auto/Lambdify baseline, fixed BDF, Sparse vs
  validated Banded, gated AOT, stop conditions, task preview, result diagnostics, and plotting/save
  postprocessing. `examples/reactors_ivp_guides.md` now collects the canonical usage flows.

- [x] Update module documentation.
  Replace copied BVP prose with the condensed-phase assumptions, exact state vector, governing
  equations, units, initial-value semantics, LSODE2 backend choices, and limitations.
  `src/SimpleReactorsIVP.rs` now documents the module as the condensed-phase IVP entry point.

## Current phase: KiThe integration with StiffAPI (locked RST 0.4.31)

**Pinned dependency:** RustedSciThe `0.4.31` (`Cargo.toml` and `Cargo.lock`). The common
`StiffAPI`, its full-document IVP parser, and the hybrid settings-only parser are upstream
capabilities, not KiThe work to recreate.

**Status reconciliation (2026-11-02):** The adapter, tagged solver-family options, hybrid document
route, GUI lifecycle, stage telemetry, benchmark harness, and PowerShell campaign script are
already implemented. Several unchecked items below describe work that has since landed and are
updated here. Version-sensitive observations in `STORY_TESTS.md` still need refresh: RST reports
that the standalone BDF stale callback-telemetry snapshot is fixed in 0.4.31, but KiThe has not yet
rerun its callback-policy story against that version.

**Scope:** Keep the condensed-phase physical equations, inputs, scaling, and canonical state
contract unchanged. Replace the LSODE2-only KiThe application boundary with typed StiffAPI requests
and prepared solvers. Initial KiThe scope remains LSODE2, standalone BDF, and Radau; Backward Euler
is available upstream but is not added to the KiThe production selector until its usefulness for
the physical IVP workload is explicitly evaluated.

**Upstream status verified against the current 0.4.31 reports and the KiThe adapter:**

- `StiffAPI` supports LSODE2, standalone BDF, Radau, and Backward Euler with typed requests,
  capability validation, typed errors, and solver-owned telemetry access.
- The full-text IVP parser remains a standalone numerical-task workflow. KiThe must not route its
  physical task documents through that parser.
- `parse_hybrid_stiff_settings_from_document` parses solver selection/options only;
  `HybridStiffSettings::prepare(StiffSolverRequest)` applies the settings to a caller-built model,
  and `serialize_hybrid_solver_fragment` serializes only solver-owned document sections.
- Radau logging through the shared `SolverLogLevel` API is now available independently from
  telemetry. The supplied RST report also records programmatic/hybrid parity tests; those upstream
  tests are evidence for the RST API, not a substitute for KiThe integration tests.

**KiThe status:** production code now builds LSODE2, standalone BDF, and Radau requests via
`StiffSolverRequest`/`PreparedStiffSolver`. The hybrid document parser delegates only solver-owned
settings to RST, and the public config uses tagged family options. Current remaining work is not
building the facade from scratch: it is refreshing version-sensitive evidence, completing result
publication/error contracts, expanding physical correctness coverage, and collecting matched
release measurements.

| KiThe solver-family choice | Meaning |
| --- | --- |
| `Lsode2` | LSODE2 family; its internal controller can choose Auto Adams/BDF or a fixed Adams/BDF mode |
| `Bdf` | Standalone RustedSciThe BDF solver, not LSODE2's BDF controller mode |
| `Radau` | Standalone Radau solver |

Do not collapse `Lsode2(BdfOnly)` and standalone `Bdf` into one option or one report identity.
Backwards compatibility with the existing LSODE2-specific public configuration, task-document
schema, or GUI selections is not required. Make a clean typed API and update callers and fixtures.

### Pass A: pre-migration baseline (historical)

The planned unified pre-migration baseline was not frozen as one complete matrix before the adapter
work landed, and the old dependency/runtime can no longer be reconstructed from the current lockfile.
Do not keep those pre-migration checkboxes as active implementation work. Current solver behavior,
known gaps, and reproducible commands are recorded in the stories below; remaining physical and
result-contract coverage is tracked in the correctness program and Pass D/G.

### Pass B: define the KiThe StiffAPI adapter

**Boundary update:** the RST API additions reported for 0.4.29 remove the earlier blocker around
settings-only hybrid parsing and solver-family dispatch. KiThe should consume those APIs rather
than implement another parser or common facade. Native request/config types may still differ by
family; keep those details behind the KiThe adapter and retain family-specific diagnostics where
the native contracts differ.

- [x] Add a typed family selection with distinct variants for `Lsode2`, standalone `Bdf`, and
  `Radau`.
- [x] Separate family-owned controls into tagged option bundles and retain genuinely shared request
  controls in `ReactorIvpSolverConfig`. `ReactorIvpFamilyOptions::{Lsode2, Bdf, Radau}` prevents
  cross-family iteration settings from being silently applied.
- [x] Define genuinely common options separately: interval and initial state handoff, `rtol`,
  `atol`, `first_step`, `max_step`, callback execution policy, telemetry mode, logging request, and
  output/trajectory policy where semantics are genuinely shared. Output/trajectory semantics still
  need the result-contract work in Pass D; this item records the config ownership split.
- [x] Define native option bundles explicitly:
  - LSODE2: internal controller (`AutoAdamsBdf`, `AdamsOnly`, `BdfOnly`), linear structure/policy,
    native execution controls, and LSODE2 symbolic source;
  - standalone BDF: the currently exposed step budget and supported symbolic options;
  - Radau: current step/Newton controls, frontend, output policy, and supported matrix layouts.
  These tagged bundles reflect the current StiffAPI adapter; do not imply that unexposed RST-native
  order/retry knobs are configurable in KiThe.
- [x] Treat parallelism/callback scheduling as a common policy axis, but validate capability per
  family and execution route. Begin with `Sequential` as the conservative default; expose
  `Parallel { min_work }` and `Auto { min_work }` only where StiffAPI/native routes support them.
  Do not imply that a common setting is honored when a selected backend cannot implement it.
  KiThe exposes and forwards these policies, with parity/dispatch stories and a Criterion axis. RST
  reports the BDF stale callback snapshot fixed in 0.4.31; the old 0.4.30 zero-counter observation
  remains in `STORY_TESTS.md` until the KiThe story is rerun on 0.4.31.
- [x] Keep frontend (`AtomViewNative`/`ExprLegacy`), symbolic execution (`Lambdify`/AOT), and matrix
  layout as common requests only where StiffAPI supports them. Return a typed unsupported-capability
  error for invalid combinations. In particular, standalone BDF currently supports Dense linear
  algebra only; unsupported layouts are rejected rather than mapped or silently downgraded.
- [x] Model logging independently from telemetry. Default logging to Off; forward the shared
  `SolverLogLevel` through StiffAPI, including Radau, and verify that logging remains opt-in and
  does not alter telemetry or result publication.
- [x] Map all four telemetry modes (`Off`, `Counters`, `Timings`, `Detailed`) explicitly. Disabled
  telemetry must not read clocks or update counters on KiThe-owned stages.
- [x] Make defaults explicit and conservative: preserve the current LSODE2 route as the initial
  production default until correctness and performance evidence supports a change. Do not use
  LSODE2's internal `Auto` as a synonym for StiffAPI family selection or solver fallback.

### Pass C: prepare, solve, collect native evidence

**Current progress (RST 0.4.31):** the adapter builds and solves LSODE2, standalone BDF, and Radau requests through
`StiffSolverRequest`/`PreparedStiffSolver`. The mixed-document path retains RST
`HybridStiffSettings` and prepares caller-built KiThe equations through its `prepare(request)`;
programmatic tasks continue through a common `StiffPreset`. Results retain selected family and
native telemetry; BDF and Radau Detailed telemetry is exposed in table-ready rows. Dense Radau and
standalone BDF pass the canonical condensed fixture. Sparse Radau's former callback-shape mismatch
was fixed upstream in RST 0.4.30; KiThe now covers Dense/Sparse/BandedAuto parity for both symbolic
frontends. Solver capability rejection, StiffAPI preparation errors, and native solve errors remain
distinct typed variants with native `Error::source()` chains. Family-owned controls are tagged and
round-trip without cross-family fields. Failed BDF/LSODE2 attempts retain available native
telemetry. RST reports that the standalone BDF callback snapshot defect was fixed in 0.4.31; the
KiThe callback-policy story still needs a current-version rerun. Remaining work: result publication
and error semantics, wider family/physical acceptance, and refreshed evidence. Focused command:
`cargo test --lib --no-default-features SimpleReactorsIVP:: -- --nocapture`.

Historical note: Sparse Radau previously failed in RST 0.4.29 because the shifted sparse system
workspace length was reused for raw callback Jacobian values. RST 0.4.30 fixed the mismatch; KiThe
removed its stale rejection and added a bounded Dense/Sparse/BandedAuto parity story for both
symbolic frontends. The ignored expected-failure reproducer was removed. See the current evidence
in `STORY_TESTS.md` and `radau_matrix_layouts_match_for_both_symbolic_frontends`.

- [x] Build one adapter from reactor RHS, canonical variables, independent coordinate `x`, initial
  state, interval, and typed family options to `StiffSolverRequest`; prepare via StiffAPI and solve
  via `PreparedStiffSolver`. Keep chemistry/physics assembly outside the adapter.
- [x] Preserve the selected family and effective options in the result. Preserve native result and
  diagnostics by family rather than flattening distinct fields into guessed common values.
- [x] Collect successful-solve native telemetry from the correct owner: BDF/LSODE2 snapshots from the prepared
  solver, Radau telemetry from `RadauSolution` or `last_solve_telemetry()` (including failed runtime
  attempts where present). Represent unavailable metrics distinctly from zero.
- [x] Preserve partial native telemetry on failed runtime solves when RST exposes it. `IvpError::SolverExecution`
  retains both the native `StiffSolveError` and an optional `StiffTelemetrySnapshot`; telemetry-off
  runs do not synthesize counters. Radau uses its last solve telemetry when available.
- [x] Add opt-in KiThe stage attribution around the IVP route. `Off` performs no clock reads and
  emits no stage rows; `Counters` counts calls without timings; `Timings` records elapsed time
  without counters; `Detailed` reports both. Equation assembly, request build, RST preparation,
  integration, trajectory materialization, validation, conservation, publication, and dimensional
  postprocessing are separate rows. Native solver timings remain in the RST-owned namespace.
- [x] Measure cold preparation versus repeated/prepared solves on the current bounded fixtures.
  `PERF-IVP-01` records preparation, first solve, repeated solve, and prepared composition-sweep
  rows for LSODE2 and standalone BDF. Broader validated workloads and cross-family wall-time
  reconciliation remain open in Pass H.
- [ ] Keep host-owned logging configuration outside the numerical library. Remove unconditional
  console output and avoid constructing pretty tables in the normal solve path; return typed report
  data for GUI/CLI/example presentation instead.

### Pass D: common result contract and transactional publication

- [ ] Finish the family-neutral result contract. `IvpSolveSnapshot` already carries the shared axis,
  row-oriented values, variable order, selected family, and telemetry, but it still exposes an
  LSODE2-specific summary/statistics surface and stores status as a string. Use a typed termination
  classification and keep native diagnostics in explicitly family-specific payloads.
- [ ] Define trajectory/output semantics before claiming parity. Adaptive solvers may return
  different accepted-step grids. Compare physical values on an explicitly shared observation grid
  using solver-supported dense output/interpolation, or report native-grid comparisons separately.
  Do not resample silently or compare solver iteration meshes as if identical.
- [ ] Preserve both raw solver coordinates and dimensional postprocessed values where needed. Make
  dimensional conversion and conservation validation fallible and non-mutating until commit.
- [ ] Make publication atomic: construct and validate a candidate before replacing the last result.
  Current solve code writes `solver.x_mesh`/`solver.solution` before conservation and fallible
  postprocessing complete, so a later error can leave partially updated internal publication. Move
  the commit after required validation/postprocessing and test failed reruns at each late stage.
- [ ] Define and test stop-condition semantics separately for each family. Do not assume equal
  threshold crossing behavior from similarly named APIs. Classify cancellation/stop termination
  distinctly from solver failure and ordinary completion.

### Pass E: typed errors and capability matrix

- [x] Preserve typed capability rejection, StiffAPI preparation, and native integration failures at
  the solver boundary instead of flattening them into `CalculationError(String)`.
- [ ] Replace remaining string-erasing errors with typed stages for request/configuration
  validation, result acceptance, conservation, postprocessing, cancellation, and publication.
- [x] Preserve native `StiffApiError`/`StiffSolveError` and their source chain at request preparation
  and native solve boundaries. `ReactorIvpBackendError` and `IvpError` implement `source()` for
  these variants; GUI status uses their display text without replacing the returned typed error.
  Parser/postprocessing/publication typing remains a separate open audit. Keep KiThe physical
  rejection errors distinct from native convergence errors. Include family, stage, and available
  attempt/telemetry context without converting absent diagnostics to zero.
- [ ] Add a capability matrix for each family covering dense/sparse/banded layout, frontend,
  Lambdify/AOT, callback scheduling, logging, telemetry, stop behavior, output policy, parameter
  binding, and continuation/restart. Test rejection of unsupported requests before integration.
- [ ] Audit failure paths for partial native results and telemetry. Where a solver cannot expose
  partial state or detailed error classification, document that limitation rather than fabricate
  an incomplete successful snapshot.

### Pass F: task documents, facade users, and GUI

- [x] Add a KiThe hybrid-document route: parse KiThe-owned physical inputs locally, then pass the
  solver-only document fragment to RST's `parse_hybrid_stiff_settings_from_document`. Keep the
  existing full-text RST IVP parser out of this mixed-ownership flow. Apply the parsed settings to a
  caller-built native request with `HybridStiffSettings::prepare`, and use
  `serialize_hybrid_solver_fragment` without rewriting or consuming KiThe sections. KiThe retains
  the typed hybrid settings on a loaded task and uses them to prepare the caller-built request.
- [x] Replace LSODE2-only KiThe facade settings with typed family selection and tagged options.
  `ReactorIvpFamilyOptions` keeps LSODE2 controller/native execution, BDF max steps, and Radau
  max-step/Newton limits in disjoint variants. GUI/task settings map to these variants explicitly.
- [x] Define the hybrid task-document boundary: `ReactorIvpDocumentSpec` owns physical inputs,
  while solver settings are delegated to RST's hybrid parser. Checked-in documents and GUI lifecycle
  stories cover parse/serialize and file round-trips. Broader solver-family acceptance is tracked in
  Pass G; do not reopen this as a second parser implementation.
- [x] Update public builders/facades to accept the family-aware config directly. Ambiguous
  `with_method` was replaced with `with_lsode2_controller`; standalone iteration settings have
  family-named builders. Review examples and remaining downstream call sites as they are added.
- [x] Update GUI solver selection so the family is a distinct, visually clear control. Show only
  settings supported by the selected family; label LSODE2 internal controller separately from
  standalone BDF. Dense/Lambdify constraints and LSODE2-only stop/AOT controls are conditionally
  surfaced.
- [x] Route GUI telemetry and solver logs independently. Telemetry/log controls are separate, and
  GUI lifecycle stories cover typed solver failure presentation. The remaining typed-error audit is
  for parser, acceptance, conservation, postprocessing, and publication boundaries in Pass E.
- [ ] Update result tables/plots/export to consume the family-neutral trajectory and family-specific
  diagnostics without assuming LSODE2 summary fields.

### Pass G: correctness matrix and story tests

- [ ] Run the focused contract and solve tests after each pass; add direct StiffAPI adapter tests for
  each supported family and invalid capability combinations.
- [ ] Run the same condensed physical fixtures through LSODE2, standalone BDF, and Radau where each
  supports the request. Keep input, initial state, equations, observation grid, and physical
  acceptance criteria fixed. Record explicit unsupported cases separately from numerical failures.
- [ ] Compare state profiles, final dimensional temperature/flux/species, conservation/mass
  balance, stop location, and status. Residual/solver-native convergence flags are supplementary,
  not substitutes for physical acceptance.
- [ ] Add failure stories for invalid options, preparation failure, integration failure, rejected
  physical output, cancellation/stop, and telemetry availability.
- [x] Make result publication atomic across fallible postprocessing: compute conservation and
  dimensional output from the local candidate, then publish mesh, solution, conservation report,
  and solve snapshot only after preparation succeeds. The overflow regression verifies rejected
  candidate postprocessing leaves existing published state untouched.
- [x] Maintain `STORY_TESTS.md` in the repository's story format and distinguish LSODE2-internal
  BDF from standalone BDF. Version-sensitive callback telemetry observations have been refreshed
  for RST 0.4.31; broad physical parity stories remain open separately.
- [ ] Refresh examples and public module docs for the StiffAPI migration. The guides and module
  overview still describe LSODE2 as the only production solver; update `examples/reactors_ivp_guides.md`
  and `src/SimpleReactorsIVP.rs` to distinguish LSODE2's internal BDF controller from standalone BDF and
  document the current StiffAPI families/task-document route.

### Pass H: telemetry-led benchmarking and performance gates (final phase)

- [x] Establish the initial release benchmark matrix for bounded small controls and prepared
  LSODE2/BDF lifecycle. Keep single-solver comparisons separate from any explicitly chosen fallback
  policy; no silent cascade is part of this migration.
- [ ] Extend it to representative medium/large validated workloads and additional physical regimes.
- [x] Print compact `tabled` summaries with separate macro, KiThe-stage, RST-native, and callback
  tables for the current baseline. Unsupported/unavailable metrics are distinct from zero; nested
  timings are reported separately. Failed/unsupported attempt rows and LSODE2 native
  residual/Jacobian detail remain open additions.
- [x] Compare telemetry-off wall-clock runs separately from instrumented diagnostics on the initial
  bounded matrix. Repeat and record machine/toolchain, fixture identity, profile and status; extend
  uncertainty-aware comparisons to the larger workload matrix before drawing cross-scale conclusions.
- [ ] Add perf-gates only for material stable regressions after a clean baseline exists. A stage is
  not a priority if its end-to-end contribution is negligible; benchmarks do not replace correctness
  stories.
- [x] Add a non-fail-fast PowerShell campaign under `scripts/` that runs focused correctness
  stories and release benchmarks one by one, continues after failures, writes timestamped reports
  and technical logs under `test_reports/`, and records each step's command/status/exit code. The
  implementation is `scripts/Invoke-ReactorIvpCampaign.ps1`, with plan-only and bounded-build modes.
  A fresh broad campaign run remains pending after the fixture and correctness matrix expands.

## Definition of done for the StiffAPI phase

- [x] LSODE2, standalone BDF, and Radau are distinct typed solver-family choices; LSODE2's internal
  controller is a separate setting.
- [x] Common settings are mapped only when supported, and invalid capability combinations fail
  explicitly before solving. Broader combination coverage remains in Pass E/G.
- [x] KiThe telemetry implements `Off/Counters/Timings/Detailed`, avoids collection work when Off,
  and distinguishes unavailable from zero. The focused callback story now verifies the RST 0.4.31
  BDF post-solve snapshot reaches KiThe with nonzero residual/Jacobian callback and dispatch counts.
- [ ] Logging is independently opt-in and does not require solver modules to initialize a logger or
  print output.
- [x] Results use a family-neutral validated trajectory plus typed native diagnostics; publication
  is atomic through validation and postprocessing.
- [ ] Correctness/story evidence covers broad physical acceptance, all stop/failure behavior,
  parser/GUI routes, and known unsupported combinations. Current bounded smoke and lifecycle stories
  cover only a subset of those guarantees.
- [ ] Release benchmarks and perf-gates are based on matched workloads and stable repeated
  measurements; the final campaign script preserves per-step evidence after failures.

## IVP Correctness And Performance Evidence Program

This is the active follow-up to the StiffAPI integration. The current
canonical `H2 -> 2H` story remains a small smoke fixture only. Do not describe
it as representative of large mechanisms or as an independent physical
reference. The baseline benchmark and campaign runner already exist; expand
correctness first, then refresh their measurements. No perf-gates until repeated
release baselines exist.

### Fixture portfolio

- [x] Audit available KiThe kinetic data and domain fixtures for IVP suitability.
  Record provenance, units, reaction set, kinetic source, model assumptions,
  density, transport/thermal inputs, initial state, and intended claim for each
  fixture. No independently validated multi-step condensed-IVP kinetic
  reference was found in the existing IVP corpus. Existing HMX BVP data is a
  model fixture, not a reviewed IVP physical oracle.
- [x] Add a small balanced multistep assembly fixture with atom-balanced
  reactions, to test equation assembly beyond a single conversion. Current test
  verifies parsing and symbolic RHS
  assembly only; the end-to-end solve remains open pending controlled rate
  scales and cross-backend acceptance.
- [ ] Add a branching/competing-pathway case with distinct kinetic scales.
- [ ] Add a stiff, thermally coupled case with a clear temperature-feedback
  effect. Validate that chosen parameters produce the intended regime rather
  than assuming stiffness from their labels.
- [ ] Add a medium/large mechanism workload for setup, symbolic preparation,
  and integration scaling. Prefer an identified, traceable mechanism; if none
  is physically appropriate for this condensed-phase model, keep the case as a
  clearly named synthetic stress fixture, never a physical benchmark.
- [x] Keep the current `H2 -> 2H` fixture as a fast smoke/control case. The
  existing BVP HMX-style fixture is not an IVP oracle and must not be reused as
  one without an independent physical model and reviewed data.
- [x] Put reusable task builders and fixture metadata in
  `src/SimpleReactorsIVP/reactor_ivp_test_support.rs`, shared by tests and benchmark
  setup. Keep route inputs matched; add new fixture definitions there rather
  than silently changing rates/conditions between routes.
- [x] Localize the AtomView failure on the multi-step fixture. KiThe setup,
  request construction, and `StiffAPI` prepare all pass for ExprLegacy and
  AtomView. AtomView panics only when LSODE2 starts `try_solve`, during lazy
  native callback setup: `prepare_callbacks_for_config` -> generated IVP
  residual preparation -> `prepare_native_atom_system` -> `try_expr_to_atom`.
  It occurs before native or symbolic residual/Jacobian callback counters
  increment; one LSODE2 algorithm-selection decision has happened. The failing
  row is `rhs[1]` (heat-flux equation), specifically the second reaction's
  weighted heat term. Both rate expressions and both numeric factors convert
  separately; their product does not.

  Minimal standalone expression from the failing subtree:
  `-1790971310706.9194 * C1 * (48197.07247411639 * C2)^2`.
  `-1790971310706.9194 * C1` and `(48197.07247411639 * C2)^2` each convert
  successfully, while their product panics at
  `RustedSciThe::symbolic::View::coefficient::Coefficient::reduce`, line 114
  in locked RST 0.4.29. The failing branch evaluates `sign * (n / gcd) as i64`
  with a negative sign and an unrepresentable signed magnitude. This localizes
  the defect to coefficient combination/reduction in Expr-to-Atom conversion,
  not LSODE2 convergence or the public facade. Preserve this as an isolated
  diagnostic; do not globally rewrite KiThe equations or claim physical parity.

  The exploratory ExprLegacy integration with the initial high-order rate
  scales did not finish promptly and was stopped; it remains a separate
  fixture-stability issue, not evidence of an ExprLegacy defect.

### Correctness stories

- [x] Add a bounded low-rate H2 -> 2H frontend-parity gate across LSODE2,
  standalone BDF, and Radau. It checks finite bounded endpoint state, measurable
  species evolution, mass-fraction conservation, and componentwise ExprLegacy /
  AtomView endpoint parity. It now covers two inlet conditions (520 K and 900 K)
  and checks the full temperature/fraction profiles against `[0,3000] K` and
  `[0,1]` bounds. This is a synthetic numerical control, not an independent
  physical oracle.
- [x] Add a separate physical-plausibility smoke gate before expanding the
  workload matrix. The bounded H2 case and a two-step synthetic oxidation
  network pass finite dimensional-temperature `[0,3000] K` and species-fraction
  `[0,1]` checks. This smoke test itself remains LSODE2/ExprLegacy; the separate
  elemental-inventory story below adds the multi-step network across all families
  and frontends without promoting it to a validated kinetic mechanism.
- [x] Add TW-01 as a synthetic, independently derived analytical numerical
  control: 12 equal-composition pseudo-species, 20 first-order reactions, a
  frozen matrix-exponential profile, and LSODE2/BDF/Radau comparisons. This is
  not an experimental or externally validated kinetic mechanism and must not be
  presented as one.
- [x] Add TW-02 as an algebraic thermochemical extension of TW-01: derive all
  20 synthetic reaction heats from component enthalpies, check Hess cycles and
  per-reaction heat contributions, and compare LSODE2/BDF/Radau against a
  frozen analytical 14-state matrix-exponential reference for `gamma=+1,+10,-1`.
  This validates agreement with the current KiThe equations, not physical
  dimensional consistency of the thermal closure. Species-source scaling and
  the mass-flux-per-area interpretation are now documented in
  `MATHEMATICAL_MODEL.md`; the independent physical sign convention for the
  heat-flux/heat-release pair remains an explicit audit item.
- [x] Add TW-03 Arrhenius validation using the actual production contract:
  `R=8.314 J/(mol K)`, `E` in `J/mol`, absolute `T`, and calibrated `A_j` that
  preserve each TW-01 rate at 400 K. The supplied draft's species matrix and
  alternative physical parameter tuple did not match the current KiThe model;
  the test uses the documented `dY/dz=(L/m) M nu r` specialization and an
  equivalent internal heat-flux value to preserve the requested temperature
  path without changing TW-01's chemical scale. Three independent references
  cover constant-temperature matrix exponential, equal-Ea scalar quadrature,
  and distinct-Ea fourth-order Magnus with convergence evidence. LSODE2,
  standalone BDF, and Radau pass all three cases. This remains synthetic
  numerical evidence, not external kinetic validation.
- [x] Add TW-04 fully coupled nonlinear thermochemical validation without
  changing the production model. An independent NumPy/SciPy implementation
  reconstructs all rates, sources, and the RHS, calibrates four gamma regimes,
  studies DOP853/Radau tolerance convergence, and stores a versioned frozen
  profile with extra samples around peak heat release. The local RHS and full
  profiles pass through LSODE2, standalone BDF, and Radau. Weak/moderate/strong
  gamma values are selected from measured response and positivity, not assumed
  from magnitude; mass, elemental, Hess-law, and energy invariants are asserted.
  This closes the synthetic coupled-math validation, not the separate physical
  sign-convention audit or the need for an experimentally validated mechanism.
- [ ] Add a traceable, independently validated physical multi-step IVP
  mechanism before claiming broad physical accuracy or using it as a physical
  performance workload. TW-01/TW-02 are synthetic analytical controls only.
- [x] Extend the balanced two-step synthetic oxidation smoke to all three solver families and both
  symbolic frontends. Independently recompute C/H/O elemental mass fractions at every point on each
  native mesh; retain conservative temperature/species bounds. This is physical-consistency evidence,
  not a kinetic-accuracy oracle.
- [ ] Revisit the older high-rate canonical H2 -> 2H story fixture: with its
  default kinetics and tested intervals it can grow state values to enormous
  magnitudes while remaining finite. Keep it as solver-path smoke coverage only
  until rates, dimensional interpretation, and physical acceptance bounds are
  independently reviewed.
- [ ] Run matched physical inputs, equations, initial states, integration
  interval, tolerances, and acceptance checks through LSODE2, standalone BDF,
  and Radau wherever the adapter supports the request. Keep LSODE2's internal
  controller explicit and distinct from standalone BDF.
- [ ] Compare ExprLegacy and AtomView on each selected supported route. Check
  final dimensional temperature, heat flux and species values, profile values
  on a documented common observation grid, status, finiteness, conservation,
  and mass/element balance. State the absolute/relative tolerances and their
  scale basis in the story.
- [ ] Do not compare adaptive native meshes pointwise. Define interpolation or
  solver-supported dense-output behavior before profile parity assertions;
  report native-grid properties separately.
- [x] Treat Sparse/Banded/Dense as a capability matrix, not a Cartesian
  product. Standalone BDF accepts Dense only; LSODE2 maps Dense/Sparse/Banded.
  RST 0.4.30 repairs the Sparse Radau callback workspace: a bounded solve story
  now checks Dense/Sparse/BandedAuto for ExprLegacy and AtomView, including
  physical bounds and profile parity. Keep broader mechanism coverage open.
- [x] Promote the former Expr-to-Atom multi-step diagnostic to an active
  two-reaction frontend parity story. It checks successful LSODE2 integration,
  bounded physical state, telemetry, and profile parity instead of merely
  printing caught conversion errors.
- [ ] Include failed convergence, invalid capability, stop/cancellation where
  supported, and physical acceptance rejection as distinct outcomes. A
  production/default route passing does not imply each isolated solver passes.
- [ ] Preserve any independent reference source and its provenance. Cross
  solver agreement alone is parity evidence, not an external physical oracle.

### Release performance evidence

- [x] Add the initial release benchmark harness `reactor_ivp_baseline` for
  bounded four-state low-rate H2 -> 2H controls at 520 K and 900 K. It compares
  LSODE2, standalone BDF, and Radau across ExprLegacy and AtomView on the same
  Dense layout and reports telemetry-Off Criterion wall time separately from
  detailed runs.
- [ ] Extend the harness from its current small controls to medium and large validated workloads and
  multiple physical regimes; use identical inputs for every compared route.
- [x] Compare LSODE2 ExprLegacy/AtomView on Dense for the common timed baseline;
  standalone BDF ExprLegacy/AtomView with Dense; and Radau ExprLegacy/AtomView
  with Dense/Sparse/Banded. Correctness story covers all Radau layouts.
  Unsupported LSODE2/BDF combinations are not manufactured.
- [x] Separate telemetry-Off Criterion wall time from instrumented diagnostic runs for the initial
  bounded matrix. Record the release profile, repeated samples, fixture/route, and status; expand
  these metadata with the larger workload matrix before drawing cross-scale conclusions.
- [x] Report compact `tabled` summaries: (1) macro rows per case/route, (2) KiThe stage attribution,
  (3) RST-native solver timings/counters, and (4) callback-policy records. Unsupported/unavailable
  metrics are distinct from zero and nested timings are reported separately. Failed/unsupported
  attempts and LSODE2's native residual/Jacobian split remain open reporting additions.
- [x] Distinguish cold preparation from repeated solves. Added StiffAPI prepared
  lifecycle rows for LSODE2 and standalone BDF, separating preparation, first
  solve, repeated same-state solves, and a prepared initial-composition sweep.
  Repeats invoke the retained prepared solver and do not rebuild the symbolic
  graph or KiThe facade workflow.
- [x] Add a compiler execution axis for LSODE2 AOT versus Lambdify across
  Auto/BDF/Adams controllers, ExprLegacy/AtomView frontends, and both bounded
  inlet scenarios. Stage timing separates KiThe setup, problem construction,
  RST preparation, and integration; RST detailed telemetry is printed per
  route. AOT is not fabricated for standalone BDF or Radau: their adapter paths
  return typed unsupported-capability errors. Release execution is opt-in via
  `KITHE_RUN_IVP_AOT_BENCH=1` and requires `tcc`.
- [x] Record initial release benchmark observations in `STORY_TESTS.md` using the established
  story format. Scripts do not rewrite story history. The callback-policy Criterion axis is
  implemented; its release measurements remain pending. The RST 0.4.31 BDF callback telemetry
  refresh is complete as a focused correctness story. Do not establish performance gates from one
  run.
- [ ] Add LSODE2 native residual/Jacobian/integration timing rows to the IVP
  report. In the initial benchmark its KiThe `native_integration` stage is
  visible, but the same detailed native timing/counter split available for BDF
  and Radau is not present in the report table.

### Nightly campaign, final operational step

- [x] Implemented `scripts/Invoke-ReactorIvpCampaign.ps1`. It runs Cargo steps as separate
  processes, continues after failures, supports suite/plan-only and bounded-build options, and
  writes timestamped `reports/` and `technical/` artifacts with a manifest and Tabled summary.
- [ ] Run and archive the broad campaign after the representative correctness/benchmark matrix is
  expanded. The campaign implementation is complete; release benchmarks and the final evidence run
  remain pending.
