# ReactorsBVP Technical TODO

Scope: `src/ReactorsBVP/` only. This checklist intentionally treats performance issues as high priority because BVP reactor runs can be long CPU-bound workloads.

## P0: correctness, crashes, and user-facing failures

- [x] Replace parser `unwrap`/`expect`/`assert` paths with typed `ReactorError`.
  Functions:
  - `task_parser_reactor_BVP.rs::set_reactor_params_from_hashmap`
  - `task_parser_reactor_BVP.rs::parse_toleranse_and_bounds`
  - `task_parser_reactor_BVP.rs::parse_basic_settings`
  - `task_parser_reactor_BVP.rs::set_initial_guess_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  malformed task files currently can panic instead of returning a recoverable diagnostic. This is user-facing input, so panics are not acceptable.
  Done when:
  all required document fields are extracted through helper functions returning `Result`, tests cover missing section, wrong scalar type, wrong vector length, missing initial guess keys, and the parser accepts both vector and scalar-pair `bounds` encodings for compatibility.

- [x] Stop silently ignoring failed setup steps.
  Functions:
  - `task_parser_reactor_BVP.rs::set_reactor_params_from_hashmap`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  calls such as `let _ = self.fast_react_set(...)` and `let _ = self.setup_bvp()` can leave the reactor half-initialized and later fail far away from the real cause.
  Done when:
  these methods return `Result`, propagate errors with context, and tests assert that invalid reactions or invalid physical parameters fail at the parsing/setup boundary.

- [x] Validate finite positive physical quantities before derived calculations.
  Functions:
  - `SimpleReactorBVP.rs::check_task`
  - `SimpleReactorBVP.rs::mean_molar_mass`
  - `SimpleReactorBVP.rs::transport_coefficients`
  - `SimpleReactorBVP.rs::peclet_numbers`
  - `SimpleReactorBVP.rs::ideal_gas_density`
  Why:
  `inf`, `NaN`, zero denominators, or empty mass-fraction sums can produce meaningless solver inputs that look numerically valid.
  Done when:
  `M`, `P`, `Tm`, `Cp`, `Lambda`, `m`, `L`, diffusion coefficients, boundary fractions, and molar masses are checked for finite valid ranges before use.

- [x] Remove remaining panics from computational runtime paths.
  Functions:
  - `SimpleReactorBVP.rs::transport_coefficients`
  - `createBVP.rs::create_bvp_equations`
  - `createBVP.rs::set_solver_BC`
  - `SimpleReactorBVP3.rs::check_energy_balance`
  - `SimpleReactorBVP3.rs::check_material_balance`
  - `SimpleReactorBVP3.rs::postprocessing`
  Why:
  reactor setup, solving, and postprocessing should report domain errors instead of crashing the process.
  Done when:
  no non-test `panic!`, `unwrap`, or `expect` remains on data-dependent paths in `ReactorsBVP`.

- [x] Make solver execution report convergence failure explicitly.
  Functions:
  - `SimpleReactorBVP.rs::BVPSolver::solve_NRBVP_impl`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  `bvp.solve()` is called, then `get_result()` is accepted without checking whether the solver converged or produced a finite solution.
  Done when:
  solver status/result is validated, non-finite solution matrices are rejected, and a test covers failed/invalid solve output if the backend exposes status hooks.

- [x] Keep explicit postprocessing flags stable on the solved HMX regression path.
  Functions:
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  the solved BVP path should keep honoring `return_to_dimension`, `plot`, `gnuplot`, `save`, `save_to_csv`, and `no_plots_in_terminal` without regressing into panics or silent no-op behavior.
  Done when:
  a known solvable HMX task still solves with explicit postprocessing toggles and the regression test stays green.

## P1: hot-path performance and allocation pressure

Current status: the remaining copy pressure is now concentrated in owned solver/backend handoffs, not in avoidable intermediate cloning. For the current roadmap this pass is treated as complete, and the next practical work lives in P2 and P3.

- [x] Reduce large clones before BVP solver construction.
  Functions:
  - `SimpleReactorBVP.rs::BVPSolver::solve_NRBVP_impl`
  - `SimpleReactorBVP.rs::BVPSolver::solve_BVPsci`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  `eq_system`, `unknowns`, boundary-condition maps, parser data, and meshes are cloned at the solver handoff. For large mechanisms and grids this is expensive before the numerical work even starts.
  Done when:
  unavoidable backend ownership requirements are isolated in one handoff helper, avoidable clones are removed, and tests verify the handoff still preserves equation/order semantics.
  Current status:
  `BVPSolver::backend_boundary_conditions()` now centralizes the boundary-condition shape conversion, and both the direct solver path and parser-driven handoff paths reuse it. `solve_from_map()` and `solve_from_parsed()` now read parsed settings from borrowed maps instead of cloning them into a temporary parser snapshot, but the remaining `eq_system` and `unknowns` clones are still owned handoff copies.
  The parser wrappers for postprocessing and initial guess extraction also now borrow the parsed cache instead of cloning it first.
  The parser handoff also now moves the produced `x_mesh` into reactor state instead of cloning it back out of the backend snapshot.
  At this point the remaining copies are concentrated at the owned solver/backend boundary, so the item is mostly about respecting ownership rather than shaving off more avoidable clones.
  `SimpleReactorBVP3::check_energy_balance()` and `SimpleReactorBVP3::check_material_balance()` now borrow the solved state instead of cloning temporary snapshots, and `create_bvp_equations()` now borrows the thermal-effects slice and `qm` from the normalized context instead of copying them while the parallel assembly path avoids an intermediate result buffer. The transport snapshot is also now normalized into an ordered `d_ro_values` vector once per assembly, so the remaining copy cost is concentrated in the owned solver/backend boundary.
  `transport_coefficients()` now refreshes the cached `D_ro_map` in place, and `peclet_numbers()` borrows that cached map instead of keeping an extra local copy during the Peclet calculation.
  `set_solver_BC()` now builds boundary conditions from canonical `substances` and `unknowns` order, with `map_of_equations` checked as a contract instead of being used as the traversal source. That removes one more unordered read path from the hot setup flow.
  `set_reactor_params_from_hashmap()` now preallocates the main boundary, diffusion, substances, and group maps up front, so the parser setup avoids repeated reallocations while it materializes the task snapshot.

- [x] Introduce a reusable solver handoff builder.
  Functions:
  - `SimpleReactorBVP.rs::BVPSolver::solve_NRBVP_impl`
  - `task_parser_reactor_BVP.rs::solve_from_map`
  - `task_parser_reactor_BVP.rs::solve_from_parsed`
  Why:
  there are multiple almost-identical NRBVP setup paths. Duplication increases the chance that one path gets stale tolerances, bounds, argument name, or boundary conditions.
  Done when:
  one helper builds the NRBVP object from a validated `BVPSolver` snapshot and both parser-driven solve paths use it.
  Current status:
  `BVPSolver::build_nrbvp_backend()` now constructs the backend snapshot from a single `NrbvpHandoffConfig`, and `task_parser_reactor_BVP.rs::solve_and_store_nrbvp()` reuses it.

- [x] Avoid repeated `HashMap` cloning in equation and boundary-condition creation.
  Functions:
  - `createBVP.rs::create_bvp_equations`
  - `createBVP.rs::set_solver_BC`
  - `SimpleReactorBVP.rs::transport_coefficients`
  - `SimpleReactorBVP.rs::Le_number`
  Why:
  these methods clone maps/vectors that are then only read. This is small for toy tests but expensive for large species lists.
  Done when:
  read-only paths borrow data, output maps are moved into fields where possible, and tests cover unchanged equation ordering and boundary-condition mapping.
  Current status:
  `create_bvp_equations()` now borrows the normalized transport and kinetic snapshot instead of rebuilding it per species, `set_solver_BC()` reads canonical solver order instead of iterating unordered maps, and the transport cache is refreshed in place before Lewis/Peclet calculations borrow it. That leaves no meaningful P1 work in this branch for the current scope.

- [x] Cache reusable symbolic and numeric vectors during equation assembly.
  Functions:
  - `createBVP.rs::create_bvp_equations`
  Why:
  symbolic expressions are cloned and simplified inside nested reaction/species loops. Large mechanisms amplify this cost.
  Done when:
  rate expressions, molar masses, powers, and transport coefficients are prepared once, simplification is performed only where it changes expression complexity materially, and existing equation tests remain green.

- [x] Introduce policy-based symbolic RHS assembly.
  Functions:
  - `createBVP.rs::create_bvp_equations`
  - `SimpleReactorBVP.rs::SimpleReactorTask::set_symbolic_rhs_assembly_policy`
  - `SimpleReactorBVP.rs::SimpleReactorTask::with_symbolic_rhs_assembly_policy`
  Why:
  the assembly path should stay sequential for tiny systems and switch to parallel work only when the problem is large enough to benefit from Rayon overhead.
  Done when:
  `Auto` selects a sensible threshold, explicit sequential/parallel modes are available, and tests confirm both paths build the same symbolic system.

- [x] Remove unconditional debug output from hot paths.
  Functions:
  - `SimpleReactorBVP.rs::BVPSolver::debug_solution`
  - `SimpleReactorBVP.rs::setup_bvp`
  - `SimpleReactorBVP.rs::mean_molar_mass`
  - `SimpleReactorBVP.rs::peclet_numbers`
  - `task_parser_reactor_BVP.rs::set_initial_guess_from_map`
  Why:
  terminal I/O can dominate runtime and make batch runs unreadable.
  Done when:
  debug printing is behind `log` levels or explicit flags, and solver/setup methods are quiet by default.
  Current status:
  `setup_bvp()`, `Le_number()`, `solve_BVPsci()`, and the main NR solve path are now quiet, and the remaining explicit diagnostics use `log` output rather than raw stdout. The old debug helpers still exist, but they are no longer unconditional console noise.

## P1: RustedSciThe damped NRBVP backend migration

Current status:
RustedSciThe `0.4.8` exposes the modern damped BVP API through `NRBVP::new_with_options(..., DampedSolverOptions)`. KiThe now routes reactor NRBVP construction through a local facade in `solver_backend.rs`, while old positional methods remain as compatibility wrappers for scheme, strategy, tolerances, and bounds. Reactor BVPs always start from symbolic equations, so the practical user-facing matrix choices are `Banded` and `Sparse`; dense and numeric-closure routes are intentionally out of scope for this module.

- [x] Map the new RustedSciThe NRBVP API.
  RustedSciThe types and methods identified:
  - `NRBVP::new_with_options(eq_system, initial_guess, values, arg, border_conditions, t0, t_end, n_steps, options)`
  - `DampedSolverOptions::{sparse_damped, banded_damped}` as the relevant reactor routes
  - `DampedSolverOptions::{with_abs_tolerance, with_rel_tolerance, with_max_iterations, with_bounds, with_strategy_params, with_loglevel}`
  - `DampedSolverOptions::{with_scheme, forward_derivative, trapezoid_derivative, with_scheme_name}`
  - `DampedSolverOptions::{with_symbolic_assembly_backend, with_sparse_generated_backend_defaults, with_banded_generated_backend_defaults, with_banded_lambdify}`
  - `DampedSolverOptions::{with_banded_atomview_c_tcc, with_banded_atomview_c_gcc, with_banded_atomview_zig, with_banded_atomview_for_repeated_solves}`
  - `GeneratedBackendConfig` and `AotBuildPolicy`, `AotExecutionPolicy`, `AotChunkingPolicy`
  - `BvpSymbolicAssemblyBackend::{ExprLegacy, AtomView}`
  - `MatrixBackend::{Banded, SparseCol}` as the relevant generated matrix families
  - `BackendSelectionPolicy::{LambdifyOnly, AotOnly, PreferAotThenLambdify}` as the relevant symbolic execution policies
  Why:
  this API is no longer just a Newton-strategy wrapper; it now controls residual/Jacobian execution, symbolic assembly, matrix representation, AOT build policy, AOT compiler backend, and runtime chunking.
  Done when:
  the TODO records which RustedSciThe knobs must be represented in KiThe before implementation starts.

- [x] Introduce a KiThe solver facade module for damped NRBVP options.
  File:
  - `src/ReactorsBVP/solver_backend.rs`
  Public/semipublic KiThe types:
  - `ReactorBvpSolverConfig`
  - `ReactorBvpExecutionBackend` with at least `Lambdify`, `Aot`
  - `ReactorBvpMatrixBackend` with `Sparse`, `Banded`
  - `ReactorBvpSymbolicBackend` with `ExprLegacy`, `AtomView`
  - `ReactorBvpAotCompiler` with `CGcc`, `CTcc`, `Zig`
  - `ReactorBvpAotConfig` with compiler, build policy/profile, compile preset, execution policy, and chunk targets
  Why:
  KiThe should not expose every RustedSciThe internal type directly in user-facing reactor code. A local facade lets us keep stable names, document combustion-oriented defaults, and absorb future RustedSciThe API shifts in one place.
  Current status:
  the facade builds `DampedSolverOptions` for `Banded/Sparse + Lambdify` and for `Banded/Sparse + AOT` with `tcc`, `gcc`, or `zig`. It keeps `Lambdify + AtomView + Banded` as the ordinary default and leaves actual AOT compilation to RustedSciThe at solve time.
  KiThe's task parser now delegates solver-side settings parsing to RustedSciThe's native BVP parser and only keeps the reactor-physics side, `initial_guess`, and a small compatibility bridge for legacy grid-refinement spellings plus scalar short-hand bounds/tolerances.

- [x] Make the new KiThe default `Lambdify + AtomView + Banded`.
  Target default mapping:
  - start from `DampedSolverOptions::banded_damped()`
  - force lambdify execution with `with_banded_lambdify()` or equivalent `GeneratedBackendConfig::banded_lambdify_defaults()`
  - keep `BvpSymbolicAssemblyBackend::AtomView`
  - keep the existing derivative scheme until explicitly changed by user config
  Why:
  this matches the requested production default for reactor BVPs while avoiding AOT toolchain dependence in ordinary runs and CI.
  Current status:
  `NrbvpHandoffConfig::new()` uses `ReactorBvpSolverConfig::default_lambdify()`, direct reactor solves and parser-driven solves pass through that default, and tests inspect the produced RustedSciThe options without running a heavy solve.

- [x] Replace `NrbvpHandoffConfig` positional solver fields with one options object.
  Functions:
  - `SimpleReactorBVP.rs::NrbvpHandoffConfig`
  - `SimpleReactorBVP.rs::BVPSolver::build_nrbvp_backend`
  - `SimpleReactorBVP.rs::BVPSolver::solve_NRBVP_impl`
  - `task_parser_reactor_BVP.rs::solve_and_store_nrbvp`
  Implemented change:
  `NrbvpHandoffConfig` stores the KiThe facade config, and `build_nrbvp_backend()` calls `NRBVP::new_with_options(...)`.
  Compatibility:
  keep existing `solve_NRBVP*` positional methods initially, but make them construct the facade/default config internally.
  Current status:
  direct and parser-driven solver paths both construct NRBVP through `new_with_options`, and the snapshot test checks old positional inputs plus the new generated backend default.

- [x] Extend task-document parsing for the lambdify solver backend settings.
  Functions:
  - `task_parser_reactor_BVP.rs::parse_solver_settings_from_map`
  - task template in `task_parser_reactor_BVP.rs`
  Supported document keys:
  - `generated_backend`: `banded_lambdify`, `sparse_lambdify`, `banded_aot`, `sparse_aot`, `banded_aot_tcc`, `banded_aot_gcc`, `banded_aot_zig`, `sparse_aot_tcc`, `sparse_aot_gcc`, `sparse_aot_zig`
  - `symbolic_backend`: `AtomView` or `ExprLegacy`
  - `matrix_backend`: `banded`, `sparse`
  - `aot_compiler` / `aot_c_compiler`: `tcc`, `gcc`, `zig`
  - `aot_build_policy`: `use_if_available`, `build_if_missing`, `require_prebuilt`, `rebuild_always`
  - `aot_build_profile`: `release`, `debug`
  - `aot_compile_preset`: `production`, `fast_build`, `dev_fastest`
  - `aot_execution_policy`: `auto`, `sequential`, `parallel`
  - `aot_target_chunks`, `aot_residual_target_chunks`, `aot_jacobian_target_chunks`
  Why:
  the parser is currently limited to old fields such as `method: Sparse` and cannot express the modern backend matrix.
  Current status:
  the parser maps absent backend keys to the new default, maps explicit lambdify/AOT backend keys into the facade, supports RustedSciThe-like AOT lifecycle keys, and validates positive chunk counts before constructing solver options.

- [x] Wire the AOT branch behind the facade.
  Target:
  - `ReactorBvpExecutionBackend::Aot`
  - `ReactorBvpAotConfig`
  - RustedSciThe helpers such as `with_banded_atomview_c_tcc`, `with_banded_atomview_c_gcc`, and `with_banded_atomview_zig`
  Why:
  AOT must not silently fall back to lambdify when the user explicitly requests compiled callbacks.
  Current status:
  `banded_aot_*` and `sparse_aot_*` produce real `DampedSolverOptions`; `tcc/gcc` fill the C compiler field, `zig` selects the Zig codegen backend, and parser tests cover build policy, execution policy, compile preset, compiler aliasing, and chunking. Ordinary tests do not execute a compiler.

- [x] Normalize Lambdify and AOT task documents into mutually exclusive execution contracts.
  Contract:
  - Lambdify documents use `backend_policy: lambdify_only` and carry no AOT compiler, build, or execution-lifecycle fields.
  - AOT documents use `backend_policy: aot_only`; matrix and compiler choices remain explicit.
  Compatibility:
  contradictory older documents are normalized before the native RustedSciThe parser sees them, while unknown backend names remain untouched so RustedSciThe can report a typed validation error.
  Current status:
  GUI controls, task templates, parser normalization, and story/regression tests all enforce this contract. The Lambdify parser-driven solve no longer triggers AOT DLL generation through a stale `build_if_missing` policy.

- [x] Add gated AOT solve/integration tests.
  Target:
  - compact combustion fixture using `AOT + tcc + AtomView + Banded`
  - optional `gcc`, `zig`, and sparse AOT comparisons when toolchains are present
  Policy:
  these tests must be ignored or gated by explicit environment/toolchain probes, because normal CI and user machines may not have `tcc`, `gcc`, or `zig` available.
  Done when:
  unavailable toolchains are skipped with a clear reason, available toolchains solve the compact task, and solutions are compared against the lambdify default with finite-profile assertions.
  Current status:
  `reactor_bvp_matrix_tests.rs` now contains gated AOT smoke tests for `tcc`, `gcc`, and `zig` on both banded and sparse matrix routes; they skip cleanly when `KITHE_RUN_BVP_AOT_TESTS` is not set or the toolchain is unavailable. The tests now compare numeric profiles against the reference backend instead of checking only shapes.

- [x] Move old combustion solve tests to the new default.
  Targets:
  - `simple_reactor_bvp_tests.rs` combustion/HMX tests that actually run NRBVP
  - `simple_reactor_bvp_tests2.rs` direct backend-construction tests
  Required default:
  `Lambdify + AtomView + Banded`.
  Done when:
  existing real combustion solve tests pass through the facade default and still validate finite solution, mesh, balance/postprocessing behavior, and stable solver state.
  Current status:
  `simple_reactor_bvp_tests.rs` now routes the HMX combustion smoke cases through the new default facade (`Banded` with lambdify execution and AtomView symbolic assembly), and `simple_reactor_bvp_tests2.rs` now builds the direct NRBVP smoke cases on the same default path as well.

- [x] Add backend matrix tests using the same compact combustion task.
  Test variants:
  - default: `Lambdify + AtomView + Banded`
  - compatibility baseline: `Lambdify + ExprLegacy + Sparse`
  - AOT leader: `AOT + C/tcc + AtomView + Banded`
  - optional AOT comparisons: `C/gcc`, `Zig`, and sparse AOT when toolchains are available
  Test policy:
  - normal CI tests should skip unavailable AOT toolchains gracefully or mark strict AOT tests as ignored
  - compare solutions against the default with tolerances on profiles, not exact matrices
  - always assert finite solution and mesh
  Done when:
  one reusable compact combustion fixture can exercise all configured backends without duplicating reactor setup.
  Current status:
  `reactor_bvp_story_tests.rs` and `reactor_bvp_matrix_tests.rs` now provide the reusable compact HMX fixture, backend-shape checks, solve smoke tests for the default and sparse lambdify routes, a compatibility baseline for `ExprLegacy + Sparse`, and gated AOT solve/profile-comparison tests across both banded and sparse matrix routes.

- [x] Decide how much RustedSciThe task parser syntax to mirror.
  Reference:
  RustedSciThe already parses keys such as `generated_backend`, `symbolic_backend`, `matrix_backend`, `aot_codegen_backend`, `aot_c_compiler`, `aot_build_policy`, `aot_build_profile`, `aot_compile_preset`, `aot_execution_policy`, `banded_linear_solver`, and `refinement_steps`.
  Decision:
  do not mirror the full solver syntax in KiThe; let RustedSciThe own solver keys directly and keep only the narrow reactor-side aliases that belong to physics compatibility.
  Boundary rule:
  KiThe keeps physics and reactor-specific initial-condition strategy, including `process_conditions`, `boundary_conditions`, and `initial_guess`.
  RustedSciThe owns solver-engine knobs such as tolerances, bounds, symbolic backend choice, matrix backend choice, AOT policy, compiler selection, chunking, solve execution settings, and native solver-side task parsing.
  Done when:
  task templates and examples use one canonical spelling for solver settings, solver aliases are tested natively in RustedSciThe, and KiThe no longer needs to duplicate backend-specific key names beyond the physics-side compatibility layer.

- [x] Finish the solver-settings handoff in `task_parser_reactor_BVP.rs`.
  Local KiThe-only responsibilities after the next RustedSciThe pass:
  - keep `set_reactor_params_from_hashmap`
  - keep `parse_physics_run_settings_from_map`
  - keep `set_initial_guess_from_map_data`
  - keep reactor-side post-solve output helpers only if they are still clearly reactor-specific
  Remaining local bridge:
  - `normalize_reactor_physics_task_map` only for physics-side shorthand cleanup (`grid_refinement` spelling and `C/J/Teta/q` family expansion)
  Removed local solver-settings code:
  - `parse_solver_tolerance_and_bounds_from_map`
  - `parse_solver_backend_config_from_map`
  - `apply_native_solver_settings`
  - the temporary solver-side compatibility branch inside `solve_and_store_nrbvp`
  Migration rule:
  move the typed solver settings into RustedSciThe first, then simplify KiThe solve paths to a single handoff call that only passes physics + initial guess + raw solver-settings document.
  Current status:
  KiThe now parses the solver contract through typed RustedSciThe entry points and no longer performs local solver-settings translation on the hot solve path:
  - `task_parser_bvp::parse_bvp_solver_settings_from_document()` for solver-selection and generated-backend settings
  - `task_parser_damped::parse_bvp_damped_solver_settings_from_document()` for damped solver scalars, bounds, and tolerances
  The remaining local bridge is now only about reactor-side compatibility normalization, especially `grid_refinement` spelling cleanup and the `C/J/Teta/q` family expansion into full reactor-variable maps.
  The removed local helpers are `parse_solver_tolerance_and_bounds_from_map`, `apply_native_solver_settings`, `parse_solver_backend_config_from_spec`, and the old solver-facade bridge inside `solve_and_store_nrbvp`, so the manual solver-settings parsing path is gone.
  Done when:
  KiThe solve path no longer rewrites solver settings locally, and the only remaining solver-facing code here is a thin bridge to the native RustedSciThe parser or builder.

- [ ] Native RustedSciThe parser notes for the next pass.
  Target files:
  - `task_parser_bvp.rs`
  - `task_parser_damped.rs`
  Desired shape:
  - expose `parse_*_from_document` split entry points for `problem` vs `solver_settings`
  - keep a typed `*_Spec`/`*_SettingsSpec` result instead of mutating `NRBVP` directly during parsing
  - add `try_*` methods that return `Result<..., BvpTaskError>` or a solver-specific typed error
  - keep old `set_*`/`parse_*` names only as thin compatibility wrappers if needed
  - remove `panic!` on bad user input; unknown backend/matrix/grid-refinement names must become typed errors
  - make `set_params_from_hashmap` and `set_postpocessing_from_hashmap` bridge-free and fallible in the typed path
  - keep solver settings parsing isolated from physical problem parsing
  Concrete solver settings to own natively:
  - `generated_backend`
  - `matrix_backend`
  - `symbolic_backend`
  - `aot_*`
  - `strategy_params`
  - `adaptive_strategy`
  - `grid_refinement`
  - `bounds`
  - `rel_tolerance`
  - `scheme`, `strategy`, `method`, `linear_sys_method`, `abs_tolerance`, `max_iterations`, `loglevel`, `dont_save_log`
  Specific panic sites to remove:
  - unknown grid refinement method dispatch
  - missing / malformed solver settings sections
  - postprocessing parser fallthroughs that should be typed errors instead of process aborts
  End state:
  KiThe no longer carries solver-settings compatibility helpers except for a very small migration bridge, and RustedSciThe exposes a typed parser that the KiThe reactor facade can call directly.

## P1: numerical contract and validation coverage

- [x] Strengthen initial guess validation.
  Functions:
  - `reactor_BVP_utils.rs::InitialTemplate::generate`
  - `reactor_BVP_utils.rs::InitialConfig::generate_initial_guess`
  - `reactor_BVP_utils.rs::InitialConfig::only_one_value_for_all_initial`
  Why:
  `n_steps < 2`, zero variables, or non-finite template parameters can create invalid matrices or divide by zero.
  Done when:
  initial guess creation rejects invalid grid sizes and non-finite values, with focused tests.
  Current status:
  `InitialTemplate::generate()` now rejects `n_steps < 2` and non-finite template parameters, `generate_initial_guess()` rejects empty unknown sets and short grids, and `only_one_value_for_all_initial()` validates finite inputs and non-zero dimensions.

- [x] Validate concentration and boundary-condition semantics.
  Functions:
  - `SimpleReactorBVP.rs::check_task`
  - `SimpleReactorBVP.rs::mean_molar_mass`
  - `createBVP.rs::set_solver_BC`
  Why:
  boundary mass fractions are used as physical weights. Missing, negative, non-finite, or zero-sum fractions can invalidate density, molar mass, and source terms.
  Done when:
  tests cover valid normalized fractions, zero-sum fractions, negative fractions, missing species, and optional flux/heat-flux defaults.

- [x] Validate dimensions of kinetic artifacts before equation generation.
  Functions:
  - `createBVP.rs::create_bvp_equations`
  Why:
  the function already checks some lengths, but `vec_of_molmasses`, `stecheo_reags`, `stecheo_matrx`, `K_sym_vec`, `D_ro_map`, and substances must all agree before indexing begins.
  Done when:
  all index-sensitive structures are checked before nested loops and mismatch tests fail with `ReactorError::InvalidConfiguration`.

- [x] Add invariant tests for the full `setup_bvp` pipeline.
  Functions:
  - `SimpleReactorBVP.rs::setup_bvp`
  - `createBVP.rs::create_bvp_equations`
  - `createBVP.rs::set_solver_BC`
  Why:
  this is the central state transition from user task to solver-ready system.
  Done when:
  tests assert `unknowns.len() == eq_system.len() == BorderConditions.len() == 2*n + 2`, finite `M`, finite `Pe_q`, finite `Pe_D`, and stable variable ordering.

## P2: API, UX, and maintainability

- [x] Split parser extraction helpers out of `task_parser_reactor_BVP.rs`.
  Functions:
  - `set_reactor_params_from_hashmap`
  - `parse_toleranse_and_bounds`
  - `parse_basic_settings`
  - `set_postpocessing_from_hashmap`
  - `set_initial_guess_from_map`
  Why:
  the same nested `DocumentMap` extraction pattern is repeated many times and is easy to get wrong.
  Done when:
  small helpers such as `required_f64`, `optional_string`, `required_vector`, and `required_section` exist and are covered by parser tests.
  Current status:
  `required_section()`, `required_values()`, `required_f64()`, `required_usize()`, `required_string()`, `optional_string()`, `required_vector_f64()`, and `required_bool()` now centralize the recurring `DocumentMap` extraction patterns, so the parser no longer hand-rolls the same nested lookup logic in each caller.

- [ ] Normalize naming without breaking public callers abruptly.
  Targets:
  - `Lambda`, `Diffusion`, `BorderConditions`, `Teta`, `Pe_D`, `Pe_q`, `Q`, `M`
  Why:
  mixed Rust/non-Rust naming makes the API harder to use and obscures which fields are user inputs versus derived solver state.
  Done when:
  idiomatic accessors or aliases exist, old public names are documented as compatibility fields or deprecated in a planned migration.
  Current status:
  scientific names that are standard in the domain, such as `Pe_D`, `Pe_q`, `Q`, and `M`, should stay as-is; the cleanup here is about API clarity and compatibility wrappers, not forcing the physics vocabulary into snake_case. Read-only accessors now expose `Diffusion`, `boundary_condition`, `D_ro_map`, `Pe_D`, `Pe_q`, and the scaling triple through a cleaner surface without breaking old callers.
  `scaling_config()` now exposes the grouped `ScalingConfig` snapshot directly, so callers no longer need to rebuild the dimensionless setup from `dT`, `L`, and `T_scale` by hand.

- [ ] Separate pure reporting from printing and plotting side effects.
  Functions:
  - `SimpleReactorBVP2.rs::pretty_print_task`
  - `SimpleReactorBVP2.rs::pretty_print_equations`
  - `SimpleReactorBVP2.rs::pretty_print_reaction_rates`
  - `SimpleReactorBVP3.rs::plot`
  - `SimpleReactorBVP3.rs::gnuplot`
  - `SimpleReactorBVP3.rs::plot_in_terminal`
  - `SimpleReactorBVP3.rs::save_to_file`
  - `SimpleReactorBVP3.rs::save_to_csv`
  Why:
  batch and library users need data/results without terminal output, GUI popups, or file writes.
  Done when:
  pure data-returning helpers exist, side-effect methods call them, and tests can inspect reports without printing.
  Current status:
  `task_report()`, `equation_report_rows()`, `reaction_rate_report_rows()`, `balance_report()`, `estimate_values_report()`, and `solution_render_data()` now expose the report data without touching stdout, so tests can assert on canonical rows, cached balance metrics, quick estimates, and owned solution snapshots directly. The side-effecting pretty-printers and export methods still exist as thin wrappers.

- [x] Convert postprocessing methods to return `Result`.
  Functions:
  - `SimpleReactorBVP3.rs::check_balances`
  - `SimpleReactorBVP3.rs::check_energy_balance`
  - `SimpleReactorBVP3.rs::check_material_balance`
  - `SimpleReactorBVP3.rs::postprocessing`
  - `SimpleReactorBVP3.rs::estimate_values`
  Why:
  these methods assume solution, mesh, molar masses, and variable names exist. Missing solver output should be a typed error.
  Done when:
  postprocessing can be called safely before and after solve, with tests for missing solution and missing mesh.
  Current status:
  `check_balances()`, `check_energy_balance()`, `check_material_balance()`, `postprocessing()`, and `estimate_values()` now return typed errors, and `postprocessing_report()` gives tests a pure snapshot of the transformed solution before the solver state is rewritten.

## P3: cleanup and documentation

- [x] Document expected units and dimensional assumptions in API docs.
  Targets:
  - `FastElemReact`
  - `SimpleReactorTask::set_parameters`
  - `SimpleReactorTask::set_transport_properties`
  - `SimpleReactorTask::set_operating_conditions`
  - `ScalingConfig`
  Why:
  wrong units can produce plausible but physically meaningless BVP solutions.
  Current status:
  `FastElemReact`, `SimpleReactorTask::set_parameters`, `SimpleReactorTask::set_transport_properties`, `SimpleReactorTask::set_operating_conditions`, and `ScalingConfig` now carry explicit unit-oriented docs, and the example guides show the expected dimensioned setup flow.

- [x] Fix typos and user-facing wording.
  Targets:
  - `parse_toleranse_and_bounds`
  - `set_postpocessing_from_hashmap`
  - `from_mass_fractions_to_molar_conentration`
  - `energy_balane_error_abs`
  - `energy_balane_error_rel`
  - printed text such as `numver`, `int the task`
  Why:
  small naming issues make the API and diagnostics look less reliable.
  Compatibility aliases with corrected spelling now exist for `parse_tolerance_and_bounds()`, `set_postprocessing_from_hashmap()`, and `from_mass_fractions_to_molar_concentration()`, while the legacy spellings remain as wrappers so older callers do not break immediately.
  Current status:
  The most visible user-facing wording issues in `SimpleReactorBVP3.rs`, `SimpleReactorBVP.rs`, and the parser docs were cleaned up, grid-refinement aliases now accept human-friendly spellings, and the remaining legacy spellings are explicitly marked as compatibility aliases.

- [x] Add example-guides for BVP workflows after P0/P1 stabilization.
  Suggested examples:
  - direct elementary reactions to `setup_bvp`
  - document parser to solver
  - custom tolerances and bounds
  - postprocessing without plotting
  Why:
  examples should reflect the safer API after parser and solver error handling are fixed.
  Current status:
  `examples/reactors_bvp_guides.md` now covers direct setup, parser-to-solver, custom tolerances and bounds, and postprocessing snapshots.
