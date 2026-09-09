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

## Application facade migration and ergonomic equilibrium UI

This is a separate architectural/UI plan for the equilibrium calculator. The
goal is not to wrap the existing GUI request types indefinitely. After the
migration, `EquilibriumCalculator` is the single application request builder;
the GUI owns editable state and worker lifecycle, while the engine owns
resolution, solving, fallback, and immutable result evidence.

### Target ownership

- [x] Define the target dependency direction and keep it one-way:

  ```text
  EquilibriumGuiDocument (editable form state)
              |
              v
  ValidatedEquilibriumGuiConfig
              |
              v
  EquilibriumCalculator::builder()
              |
              v
  canonical P,T/P,H workflows -> immutable calculator outcome
              |
              v
  EquilibriumGuiResultSnapshot / presentation / plots
  ```

- [x] Add one GUI-to-facade adapter that maps validated GUI values to
  `EquilibriumCalculatorBuilder`, including phase specs, phase-qualified
  inventory, repository/lookup policy, pressure and fixed standard-state
  pressure, P,T/P,H point or range, solver policy, phase control, timing,
  diagnostics, and P,T postprocessing.
- [x] Make the adapter return the facade builder or a typed adapter error; it
  must not resolve the repository, mutate the document, start a worker, or
  duplicate engine validation.
- [x] Preserve worker-owned `EquilibriumExecutionControl`, cancellation,
  diagnostic sinks, stale-ticket rejection, and transactional publication.
  Attach these concerns through facade options rather than reintroducing a
  second solve implementation in the GUI.
- [x] Convert facade outcomes to the existing GUI result snapshot only at one
  boundary. Point and P,T/P,H range outcomes now cross this boundary without
  rebuilding accepted numerical results; snapshots retain raw points and
  range-level continuation/recovery/timing evidence. The visible UI rendering
  of the additional route-specific fields remains a separate presentation task.

### Removal strategy: debloat, do not layer wrappers

- [x] Inventory every caller of `EquilibriumGuiSolveRequest`,
  `EquilibriumGuiPhRequest`, and the manual request assembly helpers in
  `equilibrium_gui_request.rs`. Production GUI now uses only `Facade`; the
  no callers remain outside the facade boundary; the manual builder and
  `EquilibriumGuiPhRequest` have been removed.
- [x] Migrate the production preparation path for P,T point and P,T range to
  the facade, with parity tests against the
  current accepted result, conservation, provenance, backend policy, timing,
  and continuation reuse.
- [x] Migrate P,H point and P,H range second, preserving nested/monolithic
  route selection and GUI cancellation/fallback evidence. The production app
  now prepares every P,T/P,H point or range through `Facade`; compatibility
  request variants remain only for focused legacy/parity tests.
- [x] Once all production GUI callers use the facade, delete the duplicate
  `EquilibriumGuiSolveRequest` variants and manual P,T/P,H assembly. The
  worker request boundary now has exactly one `Facade` variant; the manual
  P,T/P,H builders, `EquilibriumGuiPhRequest`, and their compatibility tests
  were removed rather than retained as wrappers.
- [x] Keep only narrowly scoped GUI adapters that perform view-model mapping;
  remove engine policy, layout construction, thermochemistry resolution, and
  solver invocation from GUI-specific request structs. The GUI adapter now
  maps validated state only; resolution and orchestration are facade-owned.
- [x] Retain low-level request construction only in explicitly labelled
  diagnostic/architecture examples or tests that exercise the engine boundary
  itself. The GUI tree no longer contains manual request assembly or legacy
  request variants; every prepared GUI request is `Facade`.

### Primary user-facing tab

- [x] Make the first visible tab `Problem` (or `Setup`) and keep it sufficient
  for the normal calculation path:
  calculation mode (`P,T` or `P,H`), point/range selector, pressure,
  temperature or target enthalpy, bounds/grid, inventory mode, amounts, and
  the default lookup policy. `Problem` and `Components and phases` are open
  on the initial Setup viewport; secondary sections remain collapsed.
- [x] Make the first tab usable without opening advanced settings. The default
  document builds the production facade request, keeps lookup local/offline,
  uses the production solver cascade, and leaves timing, lifecycle tracing,
  plotting, and resampling disabled. A GUI regression test pins that contract.
- [x] Show only controls relevant to the selected mode. The Problem tab now
  switches between P,T temperature controls and P,H enthalpy/bounds controls;
  range-only controls remain in the P,T range branch.
- [x] Keep validation errors adjacent to the relevant field and show a compact
  readiness summary before starting the worker: mode, pressure/P0, inventory
  source, phase count, and requested points. Section-local validation feedback
  covers problem, inventory, lookup, phase policy, solver, diagnostics, and
  postprocessing fields; the command area retains only cross-tab fallback
  issues.

### Secondary thematic tabs

- [x] `Setup` owns the normal input path: explicit substances versus
  element-defined candidate search, phase physical states/models,
  phase-qualified amounts, candidate preview, and assignment. This keeps the
  minimum runnable problem on the first tab rather than splitting it across
  advanced settings.
- [x] `Libraries`: default versus explicit library selection, permitted
  libraries, offline mode, exact-state NIST fallback, shared repository status,
  and per-component provenance preview. Network access must be visibly
  opt-in. Default/explicit policy, local catalog selection, closed permitted
  sets, explicit NIST opt-in, catalog load status/failure, and phase-qualified
  declared lookup instructions are implemented. Actual resolved provenance
  remains result-side evidence rather than draft metadata; it cannot be known
  honestly before a request is resolved.
- [x] `Solver`: default production cascade versus one backend/custom cascade,
  solver budgets, tolerances, scaling, trace-seed policy, and a deliberate
  P,H route-mode policy. `Auto` remains the serialized production default;
  explicit `Monolithic` and `Nested temperature` are exposed only for P,H and
  map directly to the calculator facade. Keep this collapsed/secondary because
  most users should not tune it.
- [x] `Phase control`: bounded lifecycle, TPD create/keep thresholds,
  hysteresis, transition/cycle budgets, and initial phase history. Fixed
  solves hide every bounded-only setting. Bounded solves expose the engine's
  outer-loop budget and the two layout-independent initial policies: positive
  inventory or all declared candidate phases. Index-based explicit exclusions
  remain intentionally outside the GUI because they are not stable document
  semantics; cycle detection itself is automatic engine safety logic.
- [x] `Diagnostics and output`: timing, lifecycle logging, retained event
  limit/range policy, display thresholds, units, table density, and optional
  P,T resampling/interpolation. These settings affect observation and
  presentation, not physical equations. Existing output controls cover timing,
  diagnostic trace retention, explicit bounded retained-event override, plot target/basis/scale, PCHIP, and detailed/compact result tables; removed
  legacy report-retention toggles that did not reach the engine or alter an
  accepted snapshot. Resampling
  is intentionally visible only for P,T ranges. Lifecycle trace controls are
  likewise visible only for bounded phase control, with range retention only
  for P,T ranges. SI labels are explicit at every result field. The optional
  mole-fraction display cutoff is presentation-only and defaults to zero;
  future unit conversion must retain source units and never alter snapshots.
- [x] `Results` remains a read-only view rather than an input tab. It now
  exposes a compact P,H range continuation summary (direction, accepted
  points, reuse, normalization recovery, phase transitions, and timing) in
  addition to the existing point diagnostics. It must
  select the correct route-specific snapshot: monolithic P,H evidence must not
  be represented as zero temperature trials, and failed runs must not replace
  the last accepted snapshot.

### UI ergonomics and regression gates

- [x] Replace the current long single-surface settings layout with a compact
  tabbed surface. The implemented surface provides Setup, Phase control,
  Libraries, Numerics, Output, and Results tabs; detailed tab-specific
  accessibility and layout tests remain below.
  tab bar and stable per-tab sections. Avoid nested cards and avoid showing
  advanced controls as permanent visual noise.
- [x] Add `egui_kittest` coverage for default first-viewport setup, tab
  switching, conditional control visibility, invalid-input focus, and facade
  request preparation for point/range P,T/P,H configurations.
- [x] Do not retain a pre-migration GUI path solely for parity stories. The
  old path has been removed after facade execution, lifecycle, snapshot, and
  real-worker stories covered the canonical route. Future evidence compares
  the GUI facade directly with canonical engine/frozen-reference outcomes.
- [x] Add lifecycle stories for facade-backed worker cancellation, stale result
  discard, failed solve rollback, diagnostic sink delivery, and plot snapshot
  publication. Cancellation, stale-result discard, failed publication rollback,
  and bounded live diagnostic delivery are covered by GUI/unit stories; the
  real offline P,T range story also constructs an embedded plot strictly from
  its accepted snapshot, without a second solve.
- [x] Add a UI smoke assertion that the default view contains the essential
  problem controls and does not expose solver/diagnostics/library internals in
  the first tab.
- [x] Update GUI documentation and examples after migration. The equilibrium
  guide now describes the facade-backed GUI tabs and explicitly identifies the
  low-level reactive-gas example as explanatory rather than an application
  integration pattern.

## Remaining GUI quality-of-life work

These items are intentionally limited to the equilibrium calculator UI. They
must not introduce a second request path or duplicate validation owned by the
calculator facade.

### Story-level GUI regression coverage

- [x] Add `egui_kittest` story tests for the complete default point workflow:
  open Setup, enter a valid explicit inventory, validate, prepare the canonical
  facade request, run it, switch to Results, and inspect the accepted snapshot.
- [x] Add a matching element-defined candidate workflow, including candidate
  preview, phase assignment, validation feedback, and successful preparation.
  The ignored local-catalog story `offline_local_element_candidate_story_keeps_catalogs_unchanged_and_renders_audit`
  covers the full real-data UI path; deterministic preparation constraints are
  covered in `equilibrium_gui_story_tests`.
- [x] Cover both P,T and P,H point routes at the facade boundary, including the visible route-specific
  result evidence. P,H must verify that monolithic evidence is shown as
  monolithic evidence and is not represented by fabricated temperature trials.
- [x] Cover P,T temperature range at the facade boundary.
- [ ] Add P,H target-range stories: ascending and descending input where
  supported, continuation summary, formulation reuse, and immutable
  publication of accepted points.
- [x] Add a negative story for invalid pressure. Missing inventory,
  invalid phase/model combinations, unsupported P,H route settings, and
  unavailable lookup data remain to be covered with the same field-local
  error contract.
- [ ] Add presentation stories for tab isolation, conditional controls,
  diagnostics trace visibility, compact/detailed result tables, display cutoff,
  PCHIP visibility, and plot creation from an accepted snapshot without a
  second solve.
- [ ] Keep expensive local-catalog worker stories ignored/release-oriented;
  keep deterministic document, validation, routing, and rendering stories in
  the ordinary test suite. Every story must assert behavior, not merely print
  the widget tree.

### Per-tab contextual help

- [x] Add a small help model keyed by stable tab/control identifiers rather than
  by translated display text. Help lookup must be pure and must not depend on
  repository state or a solver worker.
- [x] Store equilibrium-calculator help content in `src/assets/`, extending the
  existing `help_eng.*` and `help_rus.*` resources. Keep English and Russian
  resources semantically equivalent and document the selected resource format.
- [x] Provide contextual help for `Setup`: P,T versus P,H, point versus range,
  explicit substances versus element-defined search, phase state/model, and
  initial amounts.
- [x] Provide contextual help for `Phase control`: fixed versus bounded mode,
  creation/keep driving forces, hysteresis, iteration budget, and initial phase
  set. Explain that these are numerical lifecycle policies, not new physical
  models.
- [x] Provide contextual help for `Libraries`: default/explicit priorities,
  permitted libraries, offline mode, NIST opt-in, catalog status, and the
  distinction between declared lookup instructions and resolved provenance.
- [x] Provide contextual help for `Numerics`: production cascade, explicit
  backend selection, P,H route mode, iteration budgets, scaling, and trace
  seed policy. State clearly that advanced values are optional overrides.
- [x] Provide contextual help for `Output`: timing, lifecycle retention, event
  limits, result basis, table density, display cutoff, plotting, and PCHIP.
  Explain that these settings affect observation/presentation only and do not
  alter accepted equilibrium values.
- [x] Provide contextual help for `Results`: accepted snapshots, route-specific
  evidence, conservation/residual diagnostics, provenance, phase transitions,
  continuation reuse, and fallback attempts.
- [x] Add `egui_kittest` checks that every visible primary/secondary tab exposes a help
  entry, that help is reachable without changing the document, and that
  missing/unknown keys fail safely without blocking calculation.
- [x] Add resource-loading and language-fallback tests. Missing optional help
  text must degrade to a concise stable label; it must never invalidate or
  modify a prepared facade request.

## Remaining equilibrium GUI story matrix

This is the release-oriented gap list for user-facing calculator behavior.
Existing model, accessibility, and ignored worker tests do not replace these
stories: each story must assert an observable contract and must not merely
print the widget tree.

### Results and presentation

- [x] Add a deterministic accepted-snapshot story covering component amounts,
  mole fractions, phase totals, conservation/residual evidence, provenance,
  backend attempts, and route-specific P,H diagnostics.
- [x] Prove that changing the editable document after publication cannot mutate
  the accepted snapshot or cause it to be replaced by an invalid result.
- [x] Add Results table stories for detailed and compact density. Verify that
  compact mode changes only visible columns and that the snapshot retains all
  physical quantities.
- [x] Add display-cutoff stories proving that a cutoff hides rows only in the
  presentation and does not alter component arrays, totals, conservation, or
  plotting data.
- [x] Add plotting stories for `None`, `Embedded`, `KiThePlot`, and `Both`, for
  point and P,T range snapshots. Verify plot creation consumes the accepted
  snapshot and does not start another solve.
- [x] Cover PCHIP, log scale, result basis, hidden series, repeated plot opening,
  and empty/no-result plot errors.

### Complete user workflows

- [x] Add a full element-defined UI story: choose elements, inspect candidate
  preview, assign a candidate to a phase, validate, prepare, run, and inspect
  the accepted result.
- [x] Add a bounded phase-control UI story covering fixed/bounded switching,
  hysteresis, initial phase policy, lifecycle budget, phase transition, and
  rendered lifecycle trace.
- [x] Add a Libraries UI story covering explicit priority, closed permitted set,
  offline mode, explicit NIST opt-in, catalog loading/failure/recovery, and
  accepted per-component provenance.
- [x] Add a Diagnostics UI story covering timing, lifecycle trace, range trace
  retention, event limit/truncation, fallback evidence, and phase transitions.

### Invalid-input and route matrix

- [x] Cover missing substances, negative/NaN/infinite amounts, duplicate phase
  IDs, duplicate phase-qualified components, invalid candidate policies, and
  empty/invalid library names.
- [x] Cover invalid P,T ranges, invalid P,H bounds/seeds, unsupported
  state/model combinations, unsupported route settings, invalid display cutoff,
  and invalid PCHIP point counts.
- [x] For every invalid story, assert field-local feedback, no prepared request,
  no worker start, and no mutation of the previous accepted snapshot.

### Help and future range support

- [ ] Add content-aware `egui_kittest` checks for every tab help section, not
  only Help widget presence. Pure resource tests check key terms today, but the
  rendered story verifies only that the Help widget exists.
- [ ] Test incomplete localized help through an injected pure lookup fixture and
  verify English fallback without changing the document or request. The current
  Russian `Setup` assertion exercises an existing section, not the fallback.
- [ ] Add P,H target-range stories after the canonical P,H batch API exists:
  continuation, formulation reuse, normalization recovery, ascending/descending
  grids, per-point diagnostics, and transactional publication.

## Equilibrium GUI hover tooltips

The combustion GUI provides short hover explanations for labels, widgets, and
buttons. Add the same low-friction affordance to the equilibrium calculator,
while keeping the full tab Help sections as the detailed reference.

- [x] Add the initial stable tooltip catalog keyed by control identifiers, not
  translated display text, and connect it to the main tabs, actions, fields,
  and list-management controls.
- [ ] Explain each control's purpose, units, accepted value/domain, default
  behavior, and whether changing it invalidates a prepared request or only
  changes presentation. Include action buttons such as candidate refresh,
  prepare, run, cancel, clear, and plot.
- [x] Provide English and Russian tooltip resources with the same keys and
  equivalent meaning. Missing localized entries must use the existing English
  fallback and must never block rendering or calculation.
- [x] Keep tooltips concise and complementary to the tab Help text: the tooltip
  should answer "what is this control?", while Help explains the workflow and
  the reason for the policy.
- [x] Use egui's native hover-tooltip/accessibility path so the text appears on
  mouse hover and remains available to keyboard/accessibility inspection.
- [ ] Add `egui_kittest` coverage that every visible interactive control has a
  non-empty tooltip, language switching changes tooltip text without mutating
  the document, and unknown/missing keys degrade safely.
- [x] Add story assertions for conditional controls: PCHIP only for P,T ranges,
  bounded phase controls only in bounded mode, and route-specific P,H controls.
- [x] Record tooltip coverage in GUI documentation and keep the catalog in sync
  when a new control is introduced or an old control is removed.

### Remaining tooltip implementation pass

- [ ] Extend the catalog and UI wiring to every ComboBox, checkbox, selectable
  option, and conditional control in the calculator. The coverage inventory
  must include temperature point/range editors, candidate policies, physical
  state/model selectors, trace-seed options, PCHIP/log/clamp controls, and all
  clear/refresh/preview/assignment actions.
- [x] Pass the selected `EquilibriumHelpLanguage` through shared field-rendering
  helpers so text-field tooltips are localized as well as tab Help text.
- [ ] Add `egui_kittest` hover assertions for the actual visible controls rather
  than only checking that catalog entries are non-empty. A missing tooltip for
  a newly rendered interactive control must fail the story.
- [x] Add conditional-visibility stories proving that PCHIP appears only for
  P,T ranges, bounded phase controls only in bounded mode, P,H route controls
  only for P,H, and range lifecycle controls only for ranges.
- [x] Add a pure coverage check mapping rendered control IDs to tooltip keys,
  with a safe explicit allow-list for intentionally tooltip-free decorative
  labels. Keep this check synchronized with the GUI control inventory.
- [x] Document the tooltip authoring rule: every new interactive control gets
  an English entry, a Russian entry or deliberate fallback, a concise purpose
  statement, and a UI/story assertion.

### Hover/help audit (2026-09-09)

Current baseline is green (`10` pure help/catalog tests and `39` equilibrium
GUI story tests), but it proves only the registered subset. The former
`RENDERED_TOOLTIP_CONTROL_IDS` list was removed because it was maintained next
to the catalog and could not independently detect a newly rendered widget that
was omitted from both lists. Rendered-widget coverage must be established by an
independent widget census plus actual hover stories.

#### P0 correctness gaps

- [x] Pass the selected `EquilibriumHelpLanguage` into the `Point` and `Range`
  controls. `render_temperature` currently ignores its `language` argument and
  requests English tooltips explicitly; add a real Russian hover assertion.
- [x] Attach field help to the editable `TextEdit` response as well as its
  label. At present hovering the actual input box provides no explanation even
  though the adjacent label does.
- [x] Introduce a strict catalog lookup for tests (`Option`/missing-key error)
  while retaining a safe generic fallback in production. The current generic
  fallback makes an unknown key look valid and allows spelling mistakes or
  missing translations to pass non-empty assertions.
- [x] Add stable field identifiers for dynamic labels. `Priority 1` and
  `Permitted 1`, plus `Element`, `Trace floor`, `Trace fraction`, and
  `Minimum trace floor`, no longer fall through to the generic tooltip. The
  full typed catalog refactor remains tracked in P2.

#### P1 missing UI wiring

- [x] Complete Setup wiring: language choices, element add/remove actions,
  every physical-state candidate checkbox, candidate target ComboBox options,
  and the simple ideal-gas preset must expose localized hover text on the
  interactive widget itself.
- [x] Complete Libraries wiring: catalog selector and its options, `Load local
  library choices`, `Engine default`, `Explicit policy`, and per-row remove
  actions need stable tooltip keys and hover assertions.
- [x] Complete Numerics wiring: P,H route options, production-default selector,
  concrete backend selector and backend options, `Use custom cascade`, custom
  cascade selectors, and both trace-seed strategy options need tooltips. The
  text must distinguish production defaults, diagnostic single-backend use,
  fallback ordering, and request invalidation.
- [x] Complete Phase-control wiring: add option-level help for both initial
  phase-set choices and state clearly that these alter the initial active-set
  policy rather than the physical inventory.
- [x] Complete Diagnostics wiring: add option-level help for phase/range trace
  retention and K_eq validation modes, including memory cost and applicability.
- [x] Complete Output wiring: plot-target options, result-table density and its
  options, and dynamic series visibility need complete catalog entries and
  rendered hover stories. Explain that all are presentation-only and do not
  invalidate the prepared request.

#### P1 help-content and test gaps

- [x] Reconcile Help terminology with visible labels (`Reference pressure`,
  `Model`, `Engine default`, and `Production default cascade`) so users can
  find the described control. A pure bilingual resource test now pins these
  rendered labels rather than relying on approximate terminology.
- [x] Add a compact control-contract table for every tab covering purpose,
  units/domain, effective default, visibility condition, and mutation class:
  request-invalidating, preparation-only, runtime diagnostics, or display-only.
  Both complete Help resources now carry mirrored tables and a pure test makes
  their localized structural headings part of the contract.
- [x] Audit Russian tooltip completeness explicitly. Every catalog key must
  either resolve to text distinct from English or appear in the deliberate
  `RUSSIAN_TOOLTIP_FALLBACK_KEYS` allow-list. The initial audit translated all
  seven implicit fallbacks, so the allow-list is currently empty.
- [x] Replace the claimed injected missing-localization story with a real pure
  fixture. `section_from_document` now permits an incomplete synthetic
  localized Markdown resource to prove English-section fallback and unknown-key
  behavior without changing a GUI document, prepared request, or worker state.
- [x] Build an independent widget census for each rendered tab/state and compare
  it with the tooltip contract. Cover default P,T, P,T range, P,H, element
  search, explicit lookup, bounded phases, custom cascade, diagnostics, PCHIP,
  accepted Results, running, failed, and cancelled states. The new baseline
  census covers every tab selector and the default Setup viewport using actual
  accessibility labels, independently from tooltip IDs. The conditional census
  now covers P,T range, P,H, element search, bounded lifecycle, custom cascade,
  and PCHIP. A Results census renders `Solving`, `Cancelling`, and `Failed`
  states; accepted snapshots and stale-publication behavior are covered by the
  local worker tests, including rendered diagnostic sections. The census is
  independent of tooltip IDs.
- [ ] Add actual hover assertions for labels, text editors, checkboxes,
  ComboBox headers and options, enabled/disabled actions, and dynamic rows.
  Also assert that language switching and all hover operations leave the
  editable document, prepared request, worker state, and accepted snapshot
  unchanged. First add a deterministic test helper for nested
  `CollapsingHeader` state: the current accesskit harness cannot reliably open
  every conditional subsection through repeated pointer/accessibility clicks.
  The baseline now verifies read-only hover behavior for a text editor, mode
  selector, action button, and bounded checkbox. ComboBox option and
  dynamic-row coverage are now present. The current renderer omits unavailable
  actions instead of rendering disabled widgets, so disabled-action hover
  requires an explicit UI decision and remains open.

#### P2 maintainability

- [x] Visually group the mutually exclusive P,T/P,H calculation modes and the
  explicit/elements inventory modes in the Setup tab. Keep the selected
  `selectable_label` state visible inside the grouped frames, and keep the
  controls request-invalidating.
- [x] Update bilingual Setup Help and mode/inventory hover resources to explain
  the facade boundary, physical elemental inventory, candidate-only search,
  and the effect of the selected calculation mode. Widget-census coverage now
  pins the visible `Calculation mode` and `Inventory input` groups.

- [x] Replace the duplicated declared/rendered string inventories with the
  canonical `EquilibriumTooltipDescriptor` catalog. Catalog lookup, strict
  validation, duplicate checks, and inventory iteration now derive from one
  source; stable string IDs remain only at the resource/UI boundary. Do not
  treat this catalog contract as evidence that every widget was rendered or
  hovered: that belongs to story tests.
- [x] Make `EquilibriumTooltipDescriptor` store a typed `EquilibriumHelpKey`
  and bind catalog iteration/strict lookup to that key. Strings are now exposed
  only through `as_str()` at the UI/resource boundary; the incremental call-site
  migration can remain ergonomic without a fragile giant enum.
- [x] Move the core action, route, field, lookup, candidate, solver, phase
  policy, and diagnostics short tooltip translations out of the large Rust
  `match` and into paired structured resources beside the full Help documents.
  The resource keys are checked for parity, catalog membership, and non-empty
  localized text at test time; lookup remains compile-time/offline.
- [x] Migrate every catalogued short tooltip translation from the Rust `match`
  into structured resources. The last five field entries now live in both
  Markdown resources, and a catalog test requires English and Russian resource
  entries for every registered key. Canonical rendering is resource-only; the
  old match is isolated behind a deprecated compatibility shim and has no
  internal callers.
- [ ] Remove the deprecated `legacy_tooltip` shim after the external-consumer
  compatibility audit. It must not be used as a fallback by the calculator.
- [x] Split tooltip tests into catalog-contract tests and rendered-widget story
  tests. `equilibrium_gui_help.rs` owns pure catalog/resource/fallback
  assertions; `equilibrium_gui_story_tests.rs` owns rendered calculator
  workflows. A green catalog test is explicitly not evidence that every widget
  was rendered or hovered.
