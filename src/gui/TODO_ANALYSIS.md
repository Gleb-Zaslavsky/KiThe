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
