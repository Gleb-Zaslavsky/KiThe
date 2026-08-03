# GUI BVP Technical TODO

Scope: `src/gui/` with emphasis on the combustion BVP menu and its tests.
The current BVP screen still behaves like a generic document editor, while the reactor
solver backend has grown into a much richer typed configuration surface.

## P0: expose the new BVP solver surface in the GUI

- [x] Split the BVP menu into clear sections for physics, solver backend, initial guess, and postprocessing.
  Why:
  the current UI mixes reactor physics and solver-engine knobs in one raw document editor, which makes the new API hard to discover.

- [x] Add first-class controls for the modern solver backend settings.
  Required fields:
  - `generated_backend`
  - `matrix_backend`
  - `symbolic_backend`
  - `aot_c_compiler`
  - `aot_build_policy`
  - `aot_build_profile`
  - `aot_compile_preset`
  - `aot_execution_policy`
  - `banded_linear_solver`
  - `refinement_steps`
  Why:
  the new combustion solver API supports lambdify and AOT execution, but the GUI still only exposes the old `scheme/method/strategy`-style fields.

- [x] Make the default visible and explicit in the UI.
  Target default:
  `Lambdify + AtomView + Banded`
  Why:
  users should see the production default immediately instead of guessing from template comments.

- [x] Remove outdated solver wording from the BVP help text.
  Current stale concepts:
  - `Dense` as a solver method
  - solver settings described only as tolerances and bounds
  - backend choices not shown in the help map
  Why:
  the UI should describe the real solver contract, not the old simplified one.

## P1: keep the GUI and the typed parser aligned

- [x] Update the BVP template to use the new canonical solver keys where appropriate.
  Why:
  the template is still the main entry point for many users, so it should mirror the typed solver API instead of hiding it.

- [x] Keep legacy documents working through compatibility aliases only.
  Why:
  old task files still exist, but the GUI should not teach stale names as the primary path.

- [x] Add regression tests for the BVP menu document generation.
  Cover:
  - default template contents
  - round-trip parsing of the new solver backend fields
  - legacy solver documents still loading
  - `gui_plot` / postprocessing behavior for solved BVP runs
  - `egui_kittest`-driven interaction checks for the combustion/BVP menu so we do not rely on manual clicking

- [x] Add story tests for the combustion screen presets and backend mode switching.
  Cover:
  - default BVP preset
  - AOT preset
  - sparse compatibility preset
  - validation of missing or malformed backend fields

- [x] Show species composition sections as first-class BVP UI instead of treating them as legacy raw document content.
  Why:
  per-species atomic composition is part of the reactor physics contract, not a compatibility-only edge case.

- [x] Surface `RUN CALCULATION` feedback directly in the GUI.
  Why:
  solver failures and validation problems should be visible to the user instead of disappearing into logs.

- [x] Split the solver backend UI into a mutually exclusive `Lambdify` / `AOT` selector.
  Why:
  the backend choice is not a free-form text field; the GUI should expose the actual finite set of supported execution modes.

- [x] Enforce the `Lambdify` / `AOT` choice as one coherent solver contract.
  Implemented behavior:
  - `Lambdify` sets `backend_policy: lambdify_only` and removes every AOT artifact-lifecycle field.
  - `AOT` sets `backend_policy: aot_only` and creates the compiler/build defaults only for that mode.
  Why:
  RustedSciThe supports execution selection and AOT artifact management as independent low-level settings, but the KiThe GUI deliberately presents them as mutually exclusive user modes. A Lambdify run must never compile or load an AOT artifact as a side effect.

- [x] Normalize `refinement_steps` to `usize` before handing the document to the solver.
  Why:
  the solver contract rejects floating-point values for this field, so the GUI must coerce the value into the expected typed form.

- [x] Expose grid refinement as a finite strategy selector instead of raw method keys.
  Covered strategies:
  - `pearson`
  - `grcarsmooke`
  - `twopnt`
  - `easy`
  - `doubleoints`
  Why:
  RustedSciThe accepts a short parser-supported strategy list, so the GUI should prevent invalid free-form method names and keep only one active refinement method.

- [x] Add a BVP save/read/save roundtrip regression test.
  Why:
  the GUI workflow relies on saving a configured task file and reading it back later; the generated document must therefore be parser-stable and must not emit duplicate keys or empty parser-invalid sections.

- [x] Keep advanced solver sections next to the main solver backend panel.
  Why:
  `rel_tolerance`, `strategy_params`, `adaptive_strategy`, and `grid_refinement` are solver settings, not reactor physics or species-composition data.

- [x] Replace the generic `linear_sys_method` optional editor with a typed override selector.
  Why:
  the normal value is `None`/Auto; the generic optional editor could create meaningless `Some(Float(0.0))` payloads for a field that RustedSciThe treats as `Option<String>`.

- [x] Use `refinement_steps = 5` as the KiThe BVP GUI/template default.
  Why:
  `0` remains a valid opt-out, but a small nonzero default is a better practical starting point for the banded solver route.

- [x] Replace the legacy `strategy_params.adaptive` editor with one complete adaptive-grid toggle.
  Why:
  RustedSciThe enables adaptation from the `adaptive_strategy` and `grid_refinement` sections. The GUI now creates or removes that complete contract atomically and seeds the only supported production version (`1`).

- [x] Collapse `method` and `matrix_backend` into one visible linear-algebra backend control.
  Why:
  Both fields select the same Sparse/Banded route for reactor BVP tasks. KiThe keeps both lower-level aliases synchronized for RustedSciThe compatibility without asking the user to configure the same choice twice.

## P2: polish and maintainability

- [ ] Reduce the amount of raw document editing needed for common BVP workflows.
- [x] Separate user-facing physics fields from advanced solver fields visually.
- [ ] Keep helper text and template comments in sync with the solver backend API.
- [ ] Remove dead wording that talks about the old backend model once the new controls exist.
- [x] Normalize `max_iterations` to `usize` before handing the document to the solver.

# Transport and Thermochemistry GUI Technical TODO

Scope: `transport_gui.rs`, `thermochemistry_gui.rs`, and their dedicated
`egui_kittest` regression modules. The database handlers and `SubsData` already
provide the typed calculation boundary; these screens should be thin, state-safe
clients of that boundary rather than independent implementations of the same
business logic.

Priority definitions:

- **P0**: realistic application crash or loss of the active GUI workflow.
- **P1**: potentially incorrect property source/result or misleading visible state.
- **P2**: architecture and validation debt likely to cause future regressions.
- **P3**: incomplete UX, wording, and maintainability work.

## P0: fail safely when the thermodynamic catalog cannot be loaded

- [x] Replace `ThermoData::new()` in both GUI constructors with the fallible
  `ThermoData::try_new()` path.
  Affected paths: `TransportApp::default` and `ThermochemistryApp::default`.
  Required behavior: a missing, malformed, or schema-invalid catalog produces a
  visible typed startup error and an inert calculation surface; it must not panic
  and terminate the entire GUI application.
  Done when: constructor tests inject missing/invalid catalog fixtures and prove
  that both screens remain renderable while exposing the load failure.
  Implemented: both GUI screens now build from `ThermoData::try_new()` and fall
  back to an inert empty catalog plus a visible startup error banner when loading
  fails. The dedicated GUI tests inject a synthetic catalog error and confirm that
  the screens stay alive and refuse to calculate.

## P1: use the canonical property aggregation boundary

- [x] Move transport calculation orchestration from raw handlers to `SubsData`.
  Affected paths: `perform_transport_calculation` and
  `calculate_heat_capacity_for_transport`.
  Required behavior: library priority, explicit search instructions, aliases,
  Cp lookup, molar mass, pressure, and transport calculation follow the same typed
  contract used by reactor code. The GUI must not choose the first matching thermo
  library or derive molar mass directly from an arbitrary display name.
  Done when: the GUI calculation path reads numeric values through canonical
  `SubsData` accessors and tests cover competing thermo sources, aliases, missing
  Cp, and a library name that is not itself a parseable chemical formula.

- [x] Move thermochemistry calculation orchestration to the same `SubsData`
  search and derived-value contract.
  Required behavior: explicit selection of a thermo library remains possible, but
  lookup, calculator initialization, typed errors, and numeric output no longer
  duplicate low-level handler code inside the GUI.
  Done when: NASA and NIST GUI stories use the aggregator and agree with direct
  `SubsData` results at the same temperature and units.

- [x] Invalidate dependent GUI state after every input-source mutation.
  Affected state: `selected_substance`, `calculated_*`, `search_results`, and
  `plot_window`.
  Required behavior: changing library, substance, temperature, pressure, or units
  cannot leave a result that appears to belong to the new inputs. Failed
  calculations clear the previous successful numeric snapshot.
  Done when: lifecycle tests cover success followed by library change, substance
  change, invalid input, unsupported temperature, and unit change.

- [x] Make the thermochemistry screen's visible contract match its calculations.
  Decision taken: remove pressure and Gibbs wording from this screen instead of
  pretending it computes `dG(T, P)`.
  Result:
  the visible contract now matches the actual canonical calculation path, which
  computes only `Cp`, `dH`, and `dS` from temperature.
  Done when: the chosen contract is explicit in the controls, output, and tests;
  no accepted field is silently ignored.

- [x] Fix the temperature-range window lifecycle and user-visible errors.
  Closing the window now persists, invalid ranges and calculation failures are
  rendered inside the GUI, and opening plots without a selected substance
  produces actionable feedback.
  Done when: `egui_kittest` covers open, calculate, failure, close, and reopen.
  Covered so far: missing-substance feedback is visible, and valid range
  calculation builds a ready plot snapshot instead of silently doing nothing.

- [x] Repair corrupted unit labels and settle on one encoding-safe notation.
  Affected labels include thermal conductivity, viscosity, Cp, and entropy units.
  Preferred fallback: ASCII forms such as `W/(m*K)`, `uPa*s`, and `J/(mol*K)` if
  source encoding cannot be guaranteed consistently.
  Done when: source scans contain no mojibake and widget/output tests assert the
  exact visible labels.

## P2: typed input and deterministic editor state

- [x] Introduce one shared typed condition parser for temperature, pressure, and
  temperature ranges.
  Required behavior: reject non-finite, non-positive, overflowing, and malformed
  values before entering a calculator; require `T0 < Tend`; return field-specific
  visible errors.
  Done when: table-driven tests cover zero, negative values, `NaN`, infinity,
  malformed strings, reversed ranges, and valid boundary values on both screens.

- [x] Separate editable inputs from read-only result snapshots.
  Required behavior: raw database JSON and calculated reports remain selectable or
  copyable but cannot be edited as if they were active application state.
  Done when: rendering tests distinguish input widgets from read-only result views.

- [x] Make substance lists deterministic and selection-aware.
  Required behavior: library keys are sorted, the selected substance is visibly
  selected, and changing a filter does not silently change calculation ownership.
  Done when: repeated construction produces the same ordering and interaction tests
  verify the selected row state.

- [x] Remove duplicate/dead search paths after the typed calculation migration.
  Affected code: `search_substance`, `search_substance_by_name`, raw JSON lookup,
  and the commented low-level usage blocks at the top of both modules.
  Done when: each screen has one search transition and one calculation transition,
  both documented with their state effects.
  Progress: the dead free-text `search_substance` methods and the old low-level
  example comments have been removed; the remaining live path is the selection-
  driven snapshot update through `search_substance_by_name`.
  Completed: both calculation paths now require an explicit selected substance;
  typed selection remains the only way to enter a calculation snapshot.

## P3: complete or remove decorative actions

- [ ] Implement or remove `Export Data` on both screens.
  Status: intentionally deferred for a future UX pass. We keep it visible as a
  placeholder for now, but it is not part of the current hardening work.
  Required behavior: an enabled export action writes a typed snapshot and keeps a
  persistent success/error status; an unavailable feature is visibly disabled.

- [ ] Implement or remove transport `Load from File`.
  Status: intentionally deferred for a future UX pass. We keep it visible as a
  placeholder for now, but it is not part of the current hardening work.
  Required behavior: loading uses a documented typed format and never appends a
  placeholder message that looks like completed work.

- [ ] Review module-level feature documentation after the functional migration.
  Remove claims that are not represented by working controls and keep comments
  aligned with the `SubsData` ownership boundary.

## Required transport and thermochemistry GUI test matrix

- [x] Transport picker exposes CEA and Aramco transport backends.
- [x] CEA and Aramco happy paths return positive finite transport values.
- [x] Thermochemistry picker excludes transport-only CEA records.
- [x] NASA happy path returns finite thermochemical values.
- [x] Drive library selection, substance selection, and calculation through real
  `egui_kittest` widget interactions instead of calling calculation methods directly.
- [x] Cover missing/malformed catalog startup without a panic.
- [x] Cover state invalidation after library, substance, condition, and unit changes.
- [x] Cover all non-finite and non-positive numeric input classes.
- [ ] Compare GUI results with canonical `SubsData` results for NASA, NIST, CEA,
  and Aramco scenarios.
  Note: NASA, CEA, and Aramco are covered by the new canonical comparison tests.
  NIST still needs a stable calculable fixture in the local catalog before we
  can assert the same GUI-vs-SubsData contract.
- [x] Cover plot-window lifecycle and visible range-calculation failures.
- [x] Assert exact encoding-safe unit labels and result text.

# Chemical Equilibrium GUI Technical TODO

This section defines the GUI boundary for the production chemical-equilibrium
API. The editor must consume only
`Thermodynamics::ChemEquilibrium::prelude`; legacy mutable equilibrium
orchestrators are not a GUI dependency.

The first implementation target is ideal fixed-pressure equilibrium with
explicitly declared phases at `P,T = const`. Bounded phase control and
temperature ranges are supported extensions of the same typed model. The
fixed-pressure, fixed-total-enthalpy (`P,H`) engine is now also connected
through a separate typed request path; the GUI must never simulate it by
silently substituting a temperature.

## Architectural decisions

- [x] Start the feature as focused modules instead of one large application
  file:
  - [x] `equilibrium_gui.rs`: first egui editor/application surface;
  - [x] `equilibrium_gui_model.rs`: serializable editor/document model and pure
    validation;
  - [x] `equilibrium_gui_request.rs`: the only conversion boundary from validated
    GUI values to production `ChemEquilibrium::prelude` requests;
  - [x] `equilibrium_gui_execution.rs`: transactional run gate and
    stale-result policy; background resolve/solve worker and
    immutable outcomes;
  - `equilibrium_gui_plot.rs`: equilibrium result series and adapters for both
    the existing embedded plot window and KiThePlot;
  - separate config, kittest, lifecycle, story, and plot-adapter test modules.

- [x] Introduce a versioned `EquilibriumGuiDocument` that owns only editable,
  serializable values. Keep transient catalog previews, repository handles,
  workers, plots, and accepted solver results outside the document.

- [x] Make `EquilibriumGuiConfig` a sum of typed policy blocks rather than a
  collection of booleans:

  ```text
  EquilibriumGuiConfig
  +-- problem: EquilibriumProblemSpec
  +-- inventory: EquilibriumInventorySpec
  +-- lookup: EquilibriumLookupPolicy
  +-- phase_mode: EquilibriumPhaseMode
  +-- solver: EquilibriumSolverPolicy
  +-- diagnostics: EquilibriumDiagnosticsPolicy
  +-- postprocessing: EquilibriumPostprocessingPolicy
  ```

- [x] Store physical values in canonical SI units after parsing. Editable text
  and display-unit choices belong to field drafts; validated requests must not
  contain unit-dependent strings.

- [x] Key every initial amount and result column by a phase-qualified component
  identity (phase plus exact record/substance key), never by a bare substance
  name. `gas::H2O` and `liquid::H2O` are distinct components even when their
  molecular compositions are equal.

- [x] Publish one immutable GUI-owned result snapshot per successful run.
  `EquilibriumGuiResultSnapshot` is built transactionally from accepted point
  or range outcomes, retains the source solution through `Arc`, and rejects
  empty or layout-inconsistent range payloads before publication.
  Point and range snapshots must retain conditions, canonical component order,
  physical moles, mole fractions, phase status, lookup provenance,
  conservation/acceptance evidence, backend attempt reports, optional
  equilibrium-constant validation, and optional timing. Plot windows receive
  `Arc` snapshots and never borrow mutable solver state.

- [x] Keep production defaults in the engine. The GUI expresses
  `ProductionDefault` explicitly and only construct detailed tolerances,
  cascades, or phase-control parameters when the user selects an advanced
  override. The remaining advanced override surface is tracked separately.

## Engine/API gaps discovered during GUI planning

- [x] Re-export every type required by the GUI request builder through
  `ChemEquilibrium::prelude`. In particular,
  `EquilibriumConstantValidationMode` is accepted by
  `EquilibriumSolveOptions` and is now exported through the production
  prelude; the GUI does not import an internal equilibrium path to reach it.

- [x] Give production point/range workflows a cloneable cooperative execution
  control with an optional typed progress sink. Cancellation is checked before
  lookup/formulation, between range points, before/after backend attempts, and
  from residual evaluation where the backend exposes that callback. A
  cancelled request returns typed `ReactionExtentError::Cancelled`; the GUI
  keeps the worker receiver until its terminal message and never publishes a
  partial range.

- [x] Report progress only at stable engine boundaries: repository lookup,
  formulation preparation, point start, and accepted point. The GUI renders
  the last typed event and does not invent percentages from elapsed time.
  Backend-attempt progress remains a solver-report concern until the backend
  contract exposes safe attempt callbacks.

## P0: typed document and request boundary

### P0.1 Thermodynamic problem

- [x] Model the thermodynamic constraint as an enum:
  - `FixedPT { pressure, reference_pressure, temperature }`;
  - `FixedPH { pressure, reference_pressure, target_total_enthalpy_j,
    temperature_bounds, seed }`.
  P,H uses extensive joules at the GUI boundary and validates the scalar
  bracket/seed before repository access. Its request resolves thermochemistry
  in the worker and publishes a read-only energy-contract report.

- [x] Model fixed-`P,T` temperature input as:
  - `Point { temperature }`;
  - `Range { start, end, point_count }`.
  Validate finite positive temperatures, finite positive pressures,
  `point_count >= 2`, non-equal endpoints, and both ascending and descending
  grids. Use the production `TemperatureGrid` for final validation.

- [x] Reserve an advanced postprocessing range policy independent from the
  solver grid: no resampling, or PCHIP resampling with a user-selected output
  grid/count and interpolation space. Resampling must never masquerade as
  additional solved equilibrium points. The display adapter now implements
  linear/log PCHIP, clamping, descending grids, and explicit rejection of
  single-point or non-positive log data.

### P0.2 Inventory modes

- [x] Model substance input as an enum, not as two simultaneously active
  panels:
  - `ExplicitPhases`: phase rows with typed physical state, activity model,
    exact component keys, and initial moles;
  - `ElementCandidates`: requested elements plus a typed candidate-selection
    policy.

- [x] Provide a simple explicit-species preset that creates one ideal-gas phase,
  while retaining an advanced phase editor. The simple mode is only a builder
  for the same canonical phase specification; it must not have a separate
  solving path. The typed document constructor and editor action now use the
  same explicit-phase request path.

- [x] Define the element-mode data contract as a two-stage workflow:
  1. preview deterministic candidates with selected/rejected reasons,
     temperature support, physical state, library, record key, and provenance;
  2. confirm candidate-to-phase assignments and edit initial component moles.
  Candidate discovery cannot directly launch a solve because the solver
  requires an explicit initial inventory. The GUI now runs local-catalog
  discovery in a background worker and renders an immutable audit projection;
  assignment controls now exist and preserve the selected record's library
  provenance; the candidate table now also edits assigned initial amounts.
  Manual inclusion toggles beyond explicit assignment remain.

- [x] Keep the candidate preview as derived state tied to the document
  fingerprint. Any change to elements, library policy, state filters,
  temperature interval, or candidate limit makes the preview stale; the
  fingerprint gate prevents it from being shown as current.

- [x] Expose candidate selection controls already supported by the engine:
  exact element set versus subset-of set, allowed physical states,
  temperature-coverage filter, deterministic candidate limit, and explicit
  candidate-to-phase assignments. Exact/subset, state filters, temperature
  bounds, element rows, and candidate limits are now editable in the GUI;
  candidate-to-phase assignment is now wired to the document and pins the
  selected library into the production lookup instructions. Manual inclusion
  toggles remain; assignment amounts and pinned provenance are editable and
  visible. Do not infer phase activity models from catalog metadata.

- [x] Validate unique phase IDs, unique phase-qualified components, finite
  non-negative initial moles, at least one positive amount, and complete
  candidate phase assignment before constructing a production request.

### P0.3 Library lookup

- [x] Model lookup as:
  - `DefaultPolicy`, using canonical repository defaults;
  - `ExplicitPolicy`, with ordered priority libraries, permitted libraries,
    explicit record instructions, and NIST fallback policy.

- [x] Populate library choices through the repository/API rather than a
  hard-coded GUI list. Preserve library order where it represents preference,
  but sort unordered catalog display lists deterministically. The editor now
  loads the local catalog in a background worker, normalizes the choices, and
  uses them to seed ordered priority/permitted entries; manual canonical-name
  editing remains available when the catalog is not loaded.

- [x] Treat library selection and physical phase as separate dimensions.
  Selecting `NASA_cond` does not itself mean “solid”, and requesting a liquid
  phase does not silently rewrite a substance name.

- [x] Display lookup/build provenance and failures per phase-qualified
  component. A partially resolved system must not become runnable.

### P0.4 Solver and phase-control policies

- [x] Model nonlinear backend selection as:
  - `ProductionDefault`;
  - `Single(SolverBackendChoice)`;
  - an advanced ordered custom cascade, if exposed.
  Include all supported RST backends and retained legacy LM/NR/TR numerical
  fallbacks, using stable display names mapped in one place.

- [x] Add typed advanced overrides for tolerance, maximum iterations, scaling,
  cascade budget, and trace-species seed policy. Defaults remain unset so the
  engine owns canonical production values.
  Progress: tolerance, iteration, scaling, trace-species seed overrides, and a
  user-owned ordered backend cascade are typed and validated. The GUI exposes
  add/remove/reorder controls and the request boundary maps the cascade to the
  engine policy, including the optional validated cascade budget. Scaling now defaults to the engine's stable `false` value;
  enabling it remains an explicit advanced choice after the real water/ice
  story exposed its sensitivity. The budget remains absent unless explicitly
  enabled by the user.

- [x] Model phase behavior as:
  - fixed declared phases;
  - bounded phase control with `PhaseControlPolicy`.
  Advanced controls include initial active phases, phase epsilon, outer-loop
  budget, and explicit or temperature-scaled hysteresis.

- [x] Enforce capability constraints in validation:
  multi-start is available only for the fixed active-set `P,T` workflow;
  range continuation and phase-control policies use their production typed
  contracts; unsupported combinations are disabled with a reason rather than
  silently downgraded. The GUI currently exposes no multi-start control, so it
  cannot be selected accidentally; the engine request boundary rejects it for
  ranges and bounded phase control. `P,H` uses a separate outer scalar solve,
  the same inner solver policy, and the same explicit phase-mode boundary. A
  future multi-start UI must add capability metadata before exposing that option.

### P0.5 Diagnostics and logging

- [x] Separate three concepts currently easy to conflate:
  - timing collection (`EquilibriumTimingMode`);
  - solver/report detail retained in the result;
  - application logging visibility.
  The default production solve remains quiet and timing-free.

- [x] Add diagnostics policy controls for timing, backend-attempt details,
  conservation/acceptance details, phase-transition reports, and independent
  equilibrium-constant validation when applicable.

- [x] Keep equilibrium-constant validation off or “when applicable” by default.
  The GUI default is `WhenApplicable`; it is an independent validator for
  suitable systems, not a mandatory production gate for arbitrary multiphase
  problems. The default is covered by a model contract test.

- [x] Never change process-global logger configuration from an equilibrium
  document. GUI “logging” means retaining and displaying this run's typed
  diagnostic events; the document and request builder contain no logger
  initialization or global logger mutation.

### P0.6 Pure validation and conversion

- [x] Implement pure, field-addressable validation returning
  `EquilibriumGuiValidationReport`, with errors and non-blocking warnings.
  Validation must not open databases, mutate libraries, start threads, or
  construct solver state.

- [x] Implement one request builder that consumes only a validated config and
  either:
  - builds `PhaseEquilibriumPipelineRequest` for a single point;
  - resolves once and builds `TemperatureRangeRequest` for a range;
  - builds a deferred P,H request which resolves one immutable system and
    `ResolvedThermochemistry` bundle inside the worker.
  No view code may assemble engine requests directly. The first implementation
  lives in `equilibrium_gui_request.rs`; repository resolution and candidate
  discovery remain execution-layer work.

- [x] Add initial unit tests for enum transitions and invalid combinations,
  including P,H bracket/seed validation, duplicate components, zero inventory,
  stale element previews, descending grids, phase-control options, and
  single-backend selection. The remaining request-builder-specific cases stay
  pending with the production conversion boundary.

- [x] Implement the first production request boundary in
  `equilibrium_gui_request.rs`: map validated phases and sparse initial moles
  to the PT pipeline or P,H request, map point/range temperatures to the
  canonical condition/grid types, and map GUI solver/phase/diagnostic policies
  to engine policies. Repository resolution and actual candidate discovery are
  intentionally still execution-layer work.

## P1: execution lifecycle and editor

### P1.1 Background execution

- [x] Define a background-worker-facing lifecycle gate for candidate preview,
  repository resolution, point solve, and range solve. Long solves must not
  block egui. Cancel requests cooperative engine cancellation and the worker
  remains owned until its terminal event is drained. `EquilibriumApp` moves a
  canonical point/range request to a background thread and drains progress plus
  one terminal event on the egui thread.

- [x] Tag every worker message with a monotonically increasing run ID and an
  input fingerprint. Discard stale previews/results after any owning input
  changes, following the reactor IVP lifecycle contract.

- [x] Model lifecycle explicitly:
  `Idle -> Resolving -> Solving -> Cancelling -> Completed/Failed`.
  Range execution also reports accepted point count and current temperature
  when the engine exposes progress. Closing the window must not publish an
  orphaned result into a later document.

- [x] Publish results transactionally. A failed point or range leaves the last
  accepted result clearly marked as belonging to the previous fingerprint; it
  never mixes new inputs with old output. The app also rejects late worker
  messages after an editor fingerprint change, and result snapshots validate
  homogeneous phase/component layouts before publication.

### P1.2 Editor workflow

- [x] Build one work-focused editor with stable sections: Problem, Components,
  Lookup, Phase policy, Solver, Diagnostics, Run status, Results. Avoid a
  marketing/landing screen and avoid nested cards. Result panels remain pending
  until immutable solve snapshots are connected.

- [x] Use segmented controls for mutually exclusive modes, checkboxes/toggles
  for binary diagnostics, numeric fields for physical values, and menus for
  solver/library policies. Advanced controls remain collapsed until selected.

- [x] Provide add/remove/reorder controls for phases, components, elements, and
  custom solver cascades. All four editors are present and invalidate prepared
  requests. Phase, component, and element rows use semantic egui scopes with a
  deterministic fallback while a row is incomplete, so reorder does not bind
  text-editor state to the vector index. The phase identity behavior is covered
  by egui-kittest; persisted row IDs are intentionally not duplicated in the
  document schema because phase IDs and catalog keys already provide the
  canonical identity.

- [x] Show candidate-preview rows as an auditable table with inclusion state,
  exact record key, phase/state, library, temperature support, and rejection
  reason. Selection and phase assignment are visually explicit; manual
  inclusion overrides beyond the engine's selected/rejected decision remain a
  separate policy question.

- [x] Keep accepted results read-only. Show point composition or range summary,
  phase statuses/transitions, residual and balance contracts, selected backend,
  fallback attempts, lookup provenance, timing, and validation status. The
  first read-only result table now shows phase-qualified component values,
  phase status, solver-attempt count, and range reuse evidence; the result view
  now also exposes first-class read-only sections for conservation, residuals,
  fallback attempts, K_eq validation, phase/acceptance evidence, and per-
  component lookup provenance. The sections read the immutable engine report;
  the GUI does not recompute acceptance metrics.

### P1.3 Main GUI integration

- [x] Register the equilibrium modules in `src/gui.rs`.

- [x] Add an independently owned `EquilibriumApp` window to `MainApp`; opening,
  closing, and reopening must preserve the current document without sharing
  mutable state with thermochemistry or reactor windows. The menu ownership
  helper, native child-window kittest, and persistence regression test now
  exercise this boundary.

- [x] Add a clearly named main-menu action. Ownership/lifecycle interaction
  tests for the `MainApp` menu owner and native child-window rendering are
  covered.

## P2: result postprocessing and both plot systems

- [x] Define one `EquilibriumPlotSeries` domain model containing temperatures,
  phase-qualified labels, units, moles, mole fractions, phase totals, and
  active/inactive masks. The role is provided by the immutable
  `EquilibriumGuiResultSnapshot` plus the shared `EquilibriumGuiPlotData`
  projection: basis-specific units and lifecycle masks are aligned with every
  column, and both plot backends consume the same metadata-bearing adapter.
  The public type name remains GUI-specific until a cross-application plotting
  domain is actually needed.

- [x] Add an initial adapter from `EquilibriumGuiPlotData` to the existing
  embedded `gui_plot::PlotWindow`/`egui_plot` path. PNG/export wiring remains
  part of the display-policy pass and does not run another solve. The editor
  now exposes an explicit Open embedded plot action for an accepted snapshot.

- [x] Add an initial KiThePlot `kithe_plot::DataSource` adapter for
  `EquilibriumGuiPlotData` without depending on experimental-kinetics
  `SampledColumns` or `PlotModel`; general display-policy integration remains.

- [x] Let the user choose embedded plot, KiThePlot, both, or neither. This is a
  display policy and must not trigger another solve; both actions now consume
  the same accepted immutable plot adapter.

- [x] Support result basis selection: component moles, component mole
  fractions, and phase totals. Use phase-qualified legend labels and preserve
  exact solver-grid values alongside optional resampled display series.

- [x] Add filtering for visible components/phases and optional logarithmic
  display handling that never rewrites stored zeros or accepted physical
  values.
  Progress: component/phase series can now be hidden through transient GUI
  display state and the common adapter rejects an empty visible selection.
  Stored results are not rewritten. Lifecycle masks are preserved through
  filtering and display resampling; `log10` is a shared display projection
  with explicit positive-value validation, so accepted zeros are never faked.

- [x] Verify point-result table behavior separately from range plots. A
  one-point solve remains one exact point in the shared plot adapter and is
  never sent through PCHIP resampling.

## P3: persistence, hardening, and deferred features

- [x] Add versioned save/load roundtrip tests for `EquilibriumGuiDocument`.
  Solver results, repository instances, workers, and open plot windows are not
  serialized; `EquilibriumApp::load_document_json` starts in a clean idle
  state and drops accepted results, workers, prepared requests, and plots.

- [x] Add `egui_kittest` coverage for:
  direct-species versus element modes; simple gas versus multiphase editing;
  default versus explicit libraries; point versus ascending/descending range;
  production default versus every concrete backend; diagnostics toggles;
  P,H controls and energy-report visibility; phase-control capability constraints; visible validation
  errors; candidate preview and invalidation.
  Progress: element mode, P,H controls, phase hysteresis controls, visible
  validation errors, candidate invalidation, custom cascade editing, and all
  nine concrete backend request preparations are covered. Range direction and
  display-resampling controls, explicit lookup policy, and diagnostics
  sections are now exercised too; postprocessing is kept
  outside the solver fingerprint. MainApp ownership and native child-window
  rendering are covered; OS-level native window behavior remains outside
  this egui test boundary.

- [x] Add lifecycle tests for cancellation, stale preview/result discard,
  success followed by editor changes, close/reopen, failed fallback cascades,
  and transactional range failure.
  Progress: cancellation, stale publication, previous-result retention after a
  current failure, document load reset, and close/reopen ownership are covered
  by direct lifecycle and egui tests. A real repository-backed range failure
  now verifies visible failure, no partial snapshot publication, and catalog
  immutability; real H2/O2/H2O P,H stories now verify inner fallback,
  all-backends-failed publication, cancellation, and rendered backend-attempt
  diagnostics. The complete ignored `equilibrium_gui_tests` story set passed
  in the release profile with serialized execution.

- [x] Complete the real-worker diagnostic matrix with a deterministic
  accepted-candidate validation-mismatch story. The true engine
  `AllBackendsFailed` story and P,H inner-fallback story are covered by the
  ignored local-catalog suite; the regular suite covers the UI and lifecycle
  contracts. The new ignored result-layer story mutates only a test copy of a
  real accepted candidate's validation payload and verifies that the GUI
  rejects it transactionally. Repository-backed stories remain explicitly
  ignored and require the local thermochemical catalogs.

- [x] Add offline story tests using stable local thermochemical fixtures for:
  ideal-gas `P,T`; a real multiphase system; element-defined candidate preview;
  same molecule in two phases; temperature continuation with and without phase
  transitions; provenance, conservation, timing, and backend-attempt display.
  Snapshot the relevant library files before/after tests and assert that GUI
  workflows never mutate them. Progress: ignored local-catalog `H2O` point,
  gas/ice multiphase point, water temperature-range, and exact H/O
  element-candidate preview stories now run through the GUI worker or the
  repository-backed preview path. The candidate story renders the assignment
  audit and snapshots the canonical substance, key, and element-index files
  before and after the workflow; it also assigns one live candidate and
  prepares the canonical request through the solver boundary. The point/range
  stories verify immutable
  provenance, phase-qualified layout, backend-attempt evidence, continuation
  reuse, and non-zero point timing. Both transition and no-transition
  continuation stories are covered, and a separate out-of-coverage range story
  verifies transactional failure without partial result publication.
  The range fixture stays at 250--270 K because the bundled `NASA_cond`
  `H2O(s)` record ends at 273.15 K; the GUI transition story stays inside
  that phase-specific coefficient coverage. The element-defined
  candidate preview now uses the live local catalog and verifies file
  immutability; the GUI water range also asserts a real phase-control
  transition. A separate no-transition continuation story remains useful.

- [x] Add plot-adapter tests asserting identical labels, x grids, values, and
  units for embedded and KiThePlot paths, plus empty/single-point/filtered data
  behavior. Embedded/KiThePlot parity, single-point PCHIP rejection,
  filtered-column behavior, unit metadata, mask preservation, and the shared
  positive-only `log10` display projection are covered. Physical zeros remain
  in the accepted snapshot and are rejected explicitly by the log projection;
  they are not rewritten to a fake finite value.

- [x] Connect `P,H = const` only after the engine-level typed workflow exists.
  The GUI uses total joules, explicit temperature bounds, worker-side
  thermochemistry resolution, and an immutable result energy contract. The
  accepted snapshot also retains the outer route, trial/iteration counts,
  fallback reason, phase transitions, formulation reuse, acceptance limit,
  and optional wall-time evidence; Results renders these diagnostics and the
  complete outer temperature-trial table in a dedicated P,H section. The
  editor fingerprint invalidates a prepared P,H
  request after target, bound, seed, or pressure changes. The
  ignored local-catalog story
  `offline_local_h2o_ph_story_publishes_energy_contract` is the real-data gate.

- [ ] Defer generic export/import formats until the point/range result schema
  stabilizes. When implemented, export immutable typed snapshots rather than
  screen text or mutable solver objects.

## First implementation sequence

1. Create the typed document/config enums and pure validation tests.
2. Implement the production request builder and engine-boundary unit tests.
3. Add worker lifecycle with injected executor fakes and stale-result tests.
4. Build the explicit-species point editor and result diagnostics.
5. Add element candidate preview and phase assignment.
6. Add typed temperature ranges and continuation reports.
7. Add the shared plot-series model and both plotting adapters.
8. Integrate `EquilibriumApp` into `MainApp`, then complete offline story and
   `egui_kittest` matrices.
