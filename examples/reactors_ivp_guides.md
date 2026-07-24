# Reactor IVP Example Guides

These guides cover the condensed-phase reactor IVP workflow in `src/ReactorsIVP/`.
They intentionally avoid the BVP model and the solid-state kinetics code.

Run the preview example with:

```powershell
cargo run --example reactors_ivp_preview
```

## 1. Build a condensed reactor task in code

Use this when you already have the mechanism in memory and want to assemble a
runnable task directly.

```rust
use KiThe::ReactorsIVP::SimpleReactorIVP::SimpleReactorTask;
use KiThe::ReactorsIVP::solver_backend::{ReactorIvpMatrixBackend, ReactorIvpMethod, ReactorIvpSolverConfig};
use RustedSciThe::numerical::LSODE2::{Lsode2AotProfile, Lsode2AotToolchain, Lsode2NativeExecutionConfig};
use std::collections::HashMap;

let mut task = SimpleReactorTask::new();
task.set_problem_name("Condensed burn");
task.set_problem_description("Minimal IVP builder example");
task.set_density(1200.0)?;
task.set_transport_properties(0.25, 1000.0);
task.set_initial_conditions(HashMap::from([
    ("T".to_string(), 450.0),
    ("q".to_string(), 0.15),
    ("A".to_string(), 0.7),
    ("B".to_string(), 0.3),
]));

task.set_solver_backend_config(
    ReactorIvpSolverConfig::default()
        .with_method(ReactorIvpMethod::Auto)
        .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
        .with_native_execution(Lsode2NativeExecutionConfig::bridge_solve())
);
```

## 2. Parse a task document and preview it

Use this when the task is stored as a document and you want to validate the
physics layer before solving.

```rust
use KiThe::ReactorsIVP::SimpleReactorIVP::SimpleReactorTask;
use KiThe::ReactorsIVP::task_parser_reactor_IVP::{
    build_canonical_condensed_reactor_ivp_task_template,
    parse_reactor_ivp_document_from_str,
};

let template = build_canonical_condensed_reactor_ivp_task_template();
let spec = parse_reactor_ivp_document_from_str(&template)?;

let mut task = SimpleReactorTask::new();
spec.apply_to_task(&mut task)?;
task.kindata.substances = vec!["A".to_string(), "B".to_string()];
spec.validate(&task.kindata.substances)?;
```

## 3. Inspect a typed preview before solving

Use this to print the canonical task snapshot and the pretty-printed symbolic
equations without mutating the source task.

```rust
use KiThe::gui::reactor_ivp_gui::build_reactor_ivp_task_preview_snapshot;
use KiThe::gui::reactor_ivp_gui::IvpGuiConfig;

let config = IvpGuiConfig::from_task(&task);
let preview = build_reactor_ivp_task_preview_snapshot(&task, &config)?;
preview.print_to_console();
```

## 4. Solve with the default route

The default solver route is Lambdify + AtomView + Sparse with LSODE2 method
selection controlled by the typed facade.

```rust
let snapshot = task.solve()?;
let view = task.latest_result_view()?;
let rows = view.solution_preview_rows()?;
```

## 5. Switch solver policies explicitly

Use the typed solver config when you need fixed BDF, fixed Adams, or AOT.
The supported matrix routes are Sparse and validated Banded.

```rust
use KiThe::ReactorsIVP::solver_backend::{
    ReactorIvpExecutionBackend, ReactorIvpMatrixBackend, ReactorIvpMethod,
    ReactorIvpSolverConfig, ReactorIvpSymbolicBackend,
};

let config = ReactorIvpSolverConfig::default()
    .with_method(ReactorIvpMethod::Bdf)
    .with_symbolic_backend(ReactorIvpSymbolicBackend::AtomView)
    .with_matrix_backend(ReactorIvpMatrixBackend::Banded)
    .with_execution_backend(ReactorIvpExecutionBackend::Aot {
        toolchain: Lsode2AotToolchain::CTcc,
        profile: Lsode2AotProfile::Release,
    });
```

## 6. Read diagnostics after the solve

Use the validated result view when you want status rows, residual information,
or result previews without reaching back into the mutable task.

```rust
let report_rows = task.latest_solve_report_rows()?;
let summary = task.latest_solve_report()?;
```

## Notes

- `ReactorIvpApp` owns the GUI lifecycle and worker handling.
- `task_parser_reactor_IVP.rs` owns the condensed-phase document boundary.
- `solver_backend.rs` owns the typed LSODE2 route selection.
- The solid-state IVP and BVP paths are separate subsystems and should not be
  mixed into this guide.
