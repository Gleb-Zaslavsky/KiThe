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
