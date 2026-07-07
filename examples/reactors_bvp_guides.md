# ReactorsBVP Example Guides

These guides cover the ordinary `src/ReactorsBVP/` workflows and intentionally avoid
`src/experimental_kinetics/`.

The snippets below are short, test-derived walkthroughs. They are meant as
practical onboarding notes for the API surface, not as a deep numerical manual.

## 1. Direct reactor setup

Use this path when you already know the kinetic mechanism and want to assemble
the BVP object in memory.

```rust
use KiThe::ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;
use KiThe::ReactorsBVP::reactor_BVP_utils::ScalingConfig;

fn main() -> Result<(), KiThe::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let mut reactor = SimpleReactorTask::new();

    // Fill in the physical inputs first.
    reactor.set_problem_name("Direct BVP example");
    reactor.set_problem_description("Manual setup without the parser");
    reactor.set_scaling(ScalingConfig::new(100.0, 0.1, 100.0))?;
    reactor.set_operating_conditions(101_325.0, 450.0, 0.01);

    // Derived reports are available before solving.
    let report = reactor.task_report();
    println!("Problem: {:?}", report.problem_name);
    println!("Scaling dT: {}", reactor.scaling_config().dT);

    Ok(())
}
```

## 2. Parser to solver workflow

Use this path when the problem definition lives in a document file.

```rust
use KiThe::ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;

fn main() -> Result<(), KiThe::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let mut reactor = SimpleReactorTask::new();

    // Parse the task document and let the reactor build the solver state.
    // The document includes process_conditions, kinetics, bounds, and postprocessing.
    reactor.solve_from_file("path/to/reactor_task.txt".into())?;

    Ok(())
}
```

## 3. Custom tolerances and bounds

Use this path when you want to build solver constraints from smaller config maps.

```rust
use KiThe::ReactorsBVP::reactor_BVP_utils::{create_bounds_map, create_tolerance_map};
use KiThe::ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;
use std::collections::HashMap;

fn main() {
    let reactor = SimpleReactorTask::new();

    // Species keys are expanded to the current mechanism automatically.
    let mut tolerance = HashMap::new();
    tolerance.insert("C".to_string(), 1e-7);
    tolerance.insert("J".to_string(), 1e-7);

    let mut bounds = HashMap::new();
    bounds.insert("C".to_string(), (-10.0, 10.0));
    bounds.insert("J".to_string(), (-1e20, 1e20));

    let tolerance_map = create_tolerance_map(tolerance, &reactor.kindata.substances);
    let bounds_map = create_bounds_map(bounds, &reactor.kindata.substances);

    println!("Tolerance keys: {:?}", tolerance_map.keys().collect::<Vec<_>>());
    println!("Bounds keys: {:?}", bounds_map.keys().collect::<Vec<_>>());
}
```

## 4. Postprocessing without plotting

Use this path when you want data snapshots but do not need terminal or file output.

```rust
use KiThe::ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;

fn main() {
    let reactor = SimpleReactorTask::new();

    // The report helpers are pure and easy to inspect in tests.
    let task = reactor.task_report();
    let balances = reactor.balance_report();
    println!("Task has {} reactions", task.reactions.len());
    println!("Energy balance error: {}", balances.energy_balane_error_abs);
}
```

## Notes

- The high-level helpers in `SimpleReactorTask` are thin wrappers around the
  canonical internal snapshots.
- The compatibility names remain available for older callers, but new code
  should prefer the read-only accessors and report snapshots.
- The naming of domain quantities such as `Pe_D`, `Pe_q`, `Q`, and `M` is
  intentionally preserved because it matches the scientific notation used by the
  model itself.
