# BVP GUI Refactoring Checklist

Scope: the reactor BVP GUI implemented in `combustion.rs`, its plot handoff in
`gui_plot.rs`, and the BVP-specific GUI/story tests. This checklist is based on
the current RustedSciThe `0.4.12` contract.

Priority definitions:

- **P0**: security issue or realistic user-data loss.
- **P1**: wrong solver configuration, unusable long-running workflow, or invalid typed data.
- **P2**: architectural debt that makes future correctness regressions likely.
- **P3**: diagnostics, wording, and maintainability improvements.

## P0: security and data safety

- [x] Reject untrusted `aot_c_compiler` values before solver execution.
  Affected paths: `render_solver_backend_choice_row`,
  `render_bvp_solver_backend_panel`, and `run_calculation`.
  Required behavior: loaded task documents may select only the explicitly supported
  C compilers. An arbitrary executable name or path must never reach RustedSciThe's
  `Command::new(...)` path.
  Done when: invalid compiler values produce a visible validation error and a test
  proves that a document-provided executable path is rejected without execution.
  Implemented: the GUI accepts only `tcc` and `gcc`; both request and execution
  paths validate the field before solver
  handoff. The story test uses `calc.exe` and verifies rejection before confirmation.

- [x] Require explicit user confirmation before the first AOT toolchain invocation.
  Required behavior: the confirmation identifies the selected codegen backend,
  compiler, and artifact build policy. Lambdify must never show this confirmation.
  Done when: GUI story tests cover accept, cancel, and Lambdify paths.
  Implemented: AOT requests are represented by a stable toolchain snapshot. New or
  changed snapshots open a modal confirmation, cancellation never enters the solver,
  and Lambdify bypasses the dialog. Story tests cover all three paths.

- [x] Prevent stale file paths from overwriting a previous task after problem changes.
  Affected path: problem-type switch inside `CombustionApp::show`.
  Required behavior: replacing the active problem resets `current_file_path`, plot,
  run status, edit buffers, and dirty state unless the entire application state is
  intentionally restored.
  Done when: a regression test loads file A, changes problem type, invokes Save, and
  proves that file A was not overwritten.
  Implemented: problem replacement reconstructs the complete `CombustionApp` state,
  clearing the file binding, plots, status, edit buffers, and AOT confirmation state.
  The regression test verifies that the previous file contents remain unchanged.

- [x] Prevent silent replacement of existing sections and fields.
  Affected paths: `Create Section` and all `Add Field` handlers.
  Required behavior: duplicate names are rejected or require explicit replacement
  confirmation.
  Done when: regression tests cover duplicate section and duplicate field paths.
  Implemented: every structured, species, and raw editor path uses shared no-replace
  insertion helpers. Duplicate attempts return an error and preserve the old payload.

- [x] Require confirmation before deleting required BVP sections or fields.
  Required behavior: deleting required physics/solver data is a distinct pending action;
  cancellation preserves the document and confirmation clears the intended entry only.
  Done when: egui interaction tests cover required-field, required-section, optional-data,
  confirm, and cancel paths.
  Implemented: `render_section` reports deletion intents instead of mutating the map.
  Structured BVP editors queue a typed request for the confirmation modal, while raw
  optional data retains a direct deletion path. Story tests cover both modal outcomes
  and optional-data deletion.

## P1: solver correctness and execution lifecycle

- [x] Replace string-fragment AOT detection with one typed execution-backend enum.
  Affected functions: `synchronize_exclusive_execution_backend`,
  `sync_generated_backend_from_controls`, and `render_generated_backend_mode_row`.
  Required behavior: all RustedSciThe presets, including `banded_aot` and
  `sparse_aot`, retain their AOT meaning. Unknown presets must be reported rather
  than displayed as Lambdify.
  Done when: a table-driven test covers every supported Lambdify/AOT preset and
  round-trips each preset without changing its execution mode.
  Implemented: `GuiExecutionBackend` recognizes the explicit RST Lambdify/AOT
  families, including unsuffixed `banded_aot` and `sparse_aot`. Unsupported presets
  remain explicit instead of being displayed as Lambdify.

- [x] Model AOT codegen backend and C compiler as separate controls.
  Required controls:
  - codegen backend: `C`, `Zig`, `Rust`;
  - C compiler, visible only for C: `tcc`, `gcc` and any other explicitly approved C compiler;
  - no `aot_c_compiler` field for Zig or Rust.
  Done when: generated documents parsed by RustedSciThe produce exactly the selected
  `AotCodegenBackend`, and matrix tests cover C/tcc, C/gcc, Zig, and Rust.
  Implemented: the codegen ComboBox exposes C/Zig/Rust, the compiler ComboBox exposes
  tcc/gcc only for C, and non-C selections remove `aot_c_compiler`. The tests inspect
  the built RustedSciThe `DampedSolverOptions`, not only document strings.

- [x] Preserve self-contained AOT preset semantics during GUI rendering.
  Affected functions: `ensure_bvp_solver_backend_defaults` and
  `sync_generated_backend_from_controls`.
  Required behavior: loading `banded_aot_gcc` must not seed tcc and rewrite the
  preset; loading `banded_aot_zig` must not override Zig with explicit C settings.
  Done when: load-render-save tests preserve all canonical AOT presets.
  Implemented: rendering no longer recomputes a preset unless the matrix control
  actually disagrees with its prefix. Preset-only gcc/zig settings are inferred into
  typed controls without rewriting the preset.

- [x] Move BVP solving off the egui event/render thread.
  Affected function: `CombustionApp::run_calculation`.
  Introduce an explicit state such as `Idle`, `Running`, `Completed`, `Failed`, and
  `Cancelling`, with a worker result channel and repaint notification.
  Required behavior: the GUI remains responsive, displays Running before the solve
  finishes, blocks duplicate runs, and supports cancellation when the backend permits it.
  Done when: story tests use a controllable fake worker to verify every state transition.
  Implemented: the GUI launches one named worker thread, polls an mpsc channel without
  blocking, schedules periodic repaint, disables duplicate Run, and owns explicit
  Idle/Running/Cancelling/Completed/Failed states. A typed worker payload transfers only
  plot-ready owned data back to egui. Deterministic fake-worker stories cover success,
  failure, duplicate Run, and cancellation.

- [ ] Add true cooperative cancellation after RustedSciThe exposes an interrupt hook.
  Current limitation: cancellation keeps the GUI responsive and discards the result, but
  an already-running damped solve continues until the backend returns because version
  `0.4.12` has no cancellation token in this solve path.

- [x] Make integer conversion strict and non-lossy.
  Affected functions: `value_as_usize`, `normalize_usize_field`,
  `render_value`, and `render_solver_backend_usize_row`.
  Required behavior: negative, non-finite, fractional, and overflowing values are
  rejected visibly. Values such as `1.9` must never become `1`, and negative editor
  values must never wrap to a huge `usize`.
  Done when: boundary tests cover zero, maximum accepted value, negative values,
  fractions, NaN, infinity, strings, and overflow.
  Implemented by the shared `task_value_conversion` boundary, visible GUI errors,
  canonical `Value::Usize` storage, and parser-side rejection tests.
 
- [x] Remove the generic `Optional -> Some(Float(0.0))` behavior from typed BVP fields.
  Immediate affected field: `solver_settings.loglevel`.
  Required behavior: optional string fields use a typed selector and restoring Some
  creates the correct inner type. `strategy` and other finite option sets should also
  use dropdowns instead of free-form text.
  Done when: None-to-Some interaction tests produce parser-valid RustedSciThe settings.
  `loglevel` and `linear_sys_method` now use typed `Option<String>` selectors;
  the generic renderer never guesses the missing inner type.

- [x] Validate the complete task before starting a potentially long solve.
  Required behavior: physics and solver parser errors are collected into a visible
  validation report; Run remains disabled for known-invalid state.
  Done when: malformed physics, invalid backend, incompatible AOT settings, and bad
  numeric values all fail before worker creation.
  Implemented: parser-owned preflight independently checks reactor physics, run
  settings, initial guess, postprocessing, tolerances, bounds, and native RST solver
  options. GUI-owned checks add strict scalar normalization and AOT security. The
  report is cached by a stable document fingerprint, displayed above Run, and keeps
  the button disabled until the known errors are corrected.

## P2: model and state architecture

- [x] Introduce a typed `BvpGuiConfig` as the structured editor model.
  It contains typed physics, solver backend, initial guess, adaptive-grid,
  advanced-solver, and postprocessing settings. Unknown compatibility sections are
  preserved separately.
  Required conversions: `try_from_document`, `to_document`, and an explicit
  compatibility/migration report.
  Implemented: a typed snapshot layer exists in `src/gui/bvp_gui_config.rs`,
  `CombustionApp` keeps a cached snapshot and migration report, the snapshot apply
  path has idempotence coverage, and the structured BVP editor now keeps the typed
  snapshot synced alongside the legacy raw panes while the raw document remains a
  compatibility layer.

- [x] Make rendering pure with respect to loaded task data.
  Affected functions: `ensure_bvp_solver_backend_defaults`,
  `normalize_bvp_adaptive_grid_settings`, and `normalize_grid_refinement_section`.
  Required behavior: migrations run once during load or explicit user actions, not on
  every frame. Unknown or future RustedSciThe settings are preserved and reported.
  Done when: render-only tests prove that an untouched document is byte-stable after
  load-render-save.
  Implemented: the solver-backend and adaptive-grid normalizers now run during
  load/seeding only. The render path no longer rewrites those sections, and the new
  canonical render regression keeps the loaded task byte-stable.

- [x] Stop auto-creating canonical structured sections during UI rendering.
  Required behavior: `render_bvp_section_group` and `render_bvp_dynamic_section_group`
  must report missing canonical sections instead of materializing them as a side
  effect of paint code.
  Done when: missing sections stay missing after render, and the UI path that repairs
  them lives in load/migration or an explicit recovery action.
  Implemented: the structured BVP renderers now show a visible warning for missing
  canonical sections and leave the document unchanged during paint. The structured
  solver UI now also uses the typed snapshot as its editor state instead of creating
  sections opportunistically during render.

- [x] Centralize backend normalization in one shared typed implementation.
  Remove duplicated interpretations from GUI helpers and reactor task normalization.
  Required behavior: GUI display, saved document, parser handoff, and direct solver
  facade all derive from the same execution/matrix/symbolic backend model.
  Done when: one backend matrix test is reused across GUI and parser boundaries.
  Implemented: the raw renderer and the reactor-task parser both delegate backend
  synchronization to the shared helper in `bvp_gui_config.rs`. That helper now also
  reconstructs backend state from legacy `backend_policy`-only documents, and the
  AOT toolchain request/validation helpers live in the same shared module.

- [x] Add explicit document lifecycle state.
  Suggested state: current path, dirty flag, load/save status, last validation report,
  and pending destructive action.
  Current status: `DocumentLifecycleState` now stores the active path, dirty flag,
  last file-operation result, and a clean document fingerprint snapshot. Save/load
  paths update that lifecycle state, the footer keeps the dirty flag synced to the
  current document fingerprint, and the GUI shows the result visibly. Problem
  switching now also uses an explicit confirmation path when the document is dirty
  or a calculation is active. Read Task and window closing now also go through
  explicit confirm paths before overwriting dirty work or discarding an active
  task. The destructive-action entry points now resync dirty state from the current
  fingerprint before deciding whether to prompt, and the edit paths sync the typed
  snapshot plus dirty flag immediately after UI mutations land. Remaining overwrite
  prompts now stay within the same explicit lifecycle flow instead of bypassing it.
  Required behavior: successful edits mark the document dirty; close, reload, problem
  switch, and overwrite paths prompt when unsaved work exists.

- [x] Return typed results from GUI file operations.
  Affected function: `save_document` and the Read Task handler.
  Required behavior: I/O and parse failures update visible GUI status and do not change
  `current_file_path` on failure. Console output is diagnostic only.
  Implemented: `save_document` and `load_document_from_path` now return typed
  `GuiFileOperationResult` values, update the lifecycle state, and show visible status
  messages. Tests cover successful save/load, save failure, and load failure without
  rebinding the active path.

- [x] Validate solver-to-plot handoff dimensions explicitly.
  Required contract: mesh length equals solution rows and unknown count equals solution
  columns. Mismatches must be reported instead of being silently truncated by `zip`.
  Done when: plot regression tests cover empty, mismatched, and valid result snapshots.
  Implemented: `PlotWindow` now checks mesh and column alignment before rendering,
  `save_plot_to_png` rejects mismatched series lengths, and dedicated regression tests
  cover mesh/solution and label/column mismatches.

## P3: UX, diagnostics, and maintainability

- [x] Replace remaining `println!` and `eprintln!` user feedback with GUI status plus
  `info!`, `warn!`, or `error!` diagnostics.
  Implemented: the remaining GUI-side stub/status prints now route through `log`
  macros, including the experimental kinetics panels and the GUI regression helper.

- [x] Update help text to describe only the exclusive KiThe backend contract.
  Remove unsupported numeric/fallback backend policies from user-facing help and stop
  describing Zig as a C compiler.
  Implemented: the BVP help map now describes the exclusive Lambdify/AOT split,
  keeps the matrix backend wording canonical, and limits AOT help to the actual
  visible backend choices.

- [x] Split `combustion.rs` into focused modules.
  Suggested boundaries: typed document model, solver backend controls, physics editor,
  file lifecycle, calculation worker, and postprocessing/plot handoff.
  Implemented: the document lifecycle and confirmation flow now live in
  `combustion_lifecycle.rs`, while calculation execution, validation, preview,
  and worker lifecycle moved to `combustion_execution.rs`. The main file now keeps
  the shared render/editor surface plus the type definitions that still belong to
  the GUI root.

- [x] Add stable widget IDs and deterministic field ordering for raw compatibility data.
  This prevents focus/state jumps caused by rendering `HashMap` entries in arbitrary order.
  Raw compatibility sections now render in sorted order, and each field row uses a
  stable `push_id` keyed by section plus field name.

- [x] Add persistent success/error feedback for plot export and task save operations.
  A label created only in the click frame is not sufficient feedback.
  Task save/load status already persists in the document lifecycle bar, and plot
  export now keeps a stable status message in the plot window itself.

- [x] Add a `Preview Task` button next to `Run` that prints a `Tabled` summary of the
  full BVP document, including `solver_settings`, and pretty-prints every equation via
  `Expr::pretty_print()` without mutating task state.
  Implemented: the BVP screen now exposes a read-only preview snapshot, prints the
  normalized document plus assembled equations through `tabled`, and updates the GUI
  status without touching the loaded task data.

## Required regression and story-test matrix

- [x] Add real ComboBox interaction tests for `Lambdify -> AOT -> Lambdify`.
- [x] Add real interaction tests for `C/tcc`, `C/gcc`, `Zig`, and `Rust`.
- [x] Add load-render-save tests for every canonical RustedSciThe generated-backend preset.
  Implemented: the regression now covers banded and sparse Lambdify, AOT tcc/gcc/zig,
  and the Rust codegen variants through real GUI render/save/load roundtrips.
- [x] Add malicious/invalid task tests for arbitrary compiler paths and unknown backends.
- [x] Add strict numeric editor tests for all solver counters and tolerances.
  Implemented: the GUI tests now cover `process_conditions.n_steps`,
  `solver_settings.max_iterations`, `solver_settings.refinement_steps`,
  `solver_settings.abs_tolerance`, `strategy_params.max_jac`,
  `strategy_params.max_damp_iter`, and `adaptive_strategy.max_refinements`,
  including fractional rejection and lossless normalization paths.
- [x] Add worker lifecycle tests for running, duplicate Run, cancellation, success, and failure.
- [x] Add successful solved-BVP GUI handoff tests, including plot enabled and disabled.
  Implemented: the regression suite now covers both post-solve plot handoff paths,
  one with `gui_plot` enabled and one with it disabled, so the solver result reaches
  the GUI in both cases without regressing the final window state.
- [x] Add `Preview Task` story coverage for the table snapshot and equation
  pretty-print path.
  Implemented: regression tests cover the preview snapshot builder, the equation
  pretty-print path, and the egui click flow for the new button.
- [ ] Keep the existing default, adaptive-grid, compatibility, and round-trip tests green.

Current baseline at checklist creation:

- `combustion_gui_tests`: 20 passed;
- `combustion_story_tests`: 15 passed;
- `gui_plot` tests: 2 passed.
