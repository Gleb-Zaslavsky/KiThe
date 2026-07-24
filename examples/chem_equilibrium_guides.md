# ChemEquilibrium Guides

This page collects the main user-facing equilibrium workflows that are now
considered canonical in `Thermodynamics/ChemEquilibrium`.

## What lives where

- `equilibrium_log_moles.rs`
  - canonical log-moles equilibrium solve
  - solver policy, backend selection, accepted solution snapshots
  - temperature-range solve entry points
- `equilibrium_constant_problem.rs`
  - independent reaction-extent validation problem
- `equilibrium_constant_solver.rs`
  - small-system validator solver for `ln(Q) - ln(K)`
- `equilibrium_constant_validation.rs`
  - backend-independent validation reports
- `equilibrium_temperature_postprocessing.rs`
  - sweep resampling and report generation for plots / exports

## 1. Canonical gas-phase equilibrium solve

Use this when you want the production equilibrium path with solver fallback,
validation, and accepted-solution publication.

See:
- `src/Thermodynamics/ChemEquilibrium/equilibrium_workflows.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_log_moles.rs`
- `examples/chem_equilibrium_gas_example.rs`

```rust,ignore
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_workflows::gas_solver;

let mut solver = gas_solver(
    vec!["CO".to_string(), "CO2".to_string(), "O2".to_string()],
    1500.0,
    101325.0,
    Solvers::LM,
    Some("info"),
    true,
)
.expect("failed to prepare equilibrium solver");

solver.solve().expect("equilibrium solve failed");
let accepted = solver.accepted_solution().expect("no accepted solution");
println!("{accepted:?}");
println!("{}", solver.moles_table());
```

## 2. Independent equilibrium-constant validation

Use this as a second opinion for small systems. It is intentionally narrower
than the canonical solver and operates in reaction-extent space.

See:
- `src/Thermodynamics/ChemEquilibrium/equilibrium_constant_problem.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_constant_solver.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_constant_validation.rs`
- `examples/chem_equilibrium_constant_validation_example.rs`

```rust,ignore
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::{
    EquilibriumConstantActivityModel, EquilibriumConstantProblem, MOLAR_GAS_CONSTANT,
};
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolver;
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::{
    EquilibriumConstantValidationTolerances, validate_equilibrium_constants,
};
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_reaction_basis::{
    ReactionBasisTolerances, ValidatedReactionBasis,
};

let basis = ValidatedReactionBasis::new(
    vec!["A2".to_string(), "A".to_string()],
    &nalgebra::DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
    nalgebra::DMatrix::from_column_slice(2, 1, &[-1.0, 2.0]),
    1,
    ReactionBasisTolerances::default(),
)
.unwrap();

let problem = EquilibriumConstantProblem::new(
    basis,
    vec![1.0, 0.0],
    vec![std::rc::Rc::new(|_| -50_000.0), std::rc::Rc::new(|_| 0.0)],
    EquilibriumConditions::new(1500.0, 2.0 * 101325.0, 101325.0).unwrap(),
    EquilibriumConstantActivityModel::IdealGas,
)
.unwrap();

let solution = EquilibriumConstantSolver::default().solve(&problem).unwrap();
let validation = validate_equilibrium_constants(
    &problem,
    &solution.moles,
    EquilibriumConstantValidationTolerances::default(),
)
.unwrap();
println!("accepted = {}", validation.accepted);
```

## 3. Temperature-sweep postprocessing

Use this when a solved temperature sweep needs a smoother export grid or a
compact textual report for logs and previews.

See:
- `src/Thermodynamics/ChemEquilibrium/equilibrium_temperature_postprocessing.rs`
- `examples/chem_equilibrium_temperature_postprocessing_example.rs`

```rust,ignore
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
    postprocess_temperature_series, TemperatureInterpolationPolicy,
    TemperatureInterpolationSpace, TemperaturePostprocessingPolicy,
    TemperatureResamplingGrid,
};

let rows = vec![
    (1000.0, vec![0.80, 0.15, 0.05]),
    (1300.0, vec![0.72, 0.20, 0.08]),
    (1600.0, vec![0.60, 0.26, 0.14]),
];

let policy = TemperaturePostprocessingPolicy {
    grid: TemperatureResamplingGrid::Uniform { points: 8 },
    interpolation: TemperatureInterpolationPolicy {
        space: TemperatureInterpolationSpace::Log,
        clamp: true,
    },
};

let report = postprocess_temperature_series(
    vec!["CO".to_string(), "CO2".to_string(), "O2".to_string()],
    &rows,
    &policy,
)
.unwrap();

println!("{}", report.render_table());
```

## 4. Solver policy and fallback order

The canonical engine keeps solver ordering explicit. The policy is part of the
contract and should be chosen intentionally, not inferred from hidden defaults.

```rust,ignore
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;

let policy = SolverPolicy::rusted_scithe_default();
println!("{:?}", policy.ordered_backends());
```

## Recommended tests

- `src/Thermodynamics/ChemEquilibrium/equilibrium_log_moles_tests.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_constant_tests.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_constant_solver_tests.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_temperature_postprocessing.rs`

