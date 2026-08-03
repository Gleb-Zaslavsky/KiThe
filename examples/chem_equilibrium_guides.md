# ChemEquilibrium Guides

This page collects the main user-facing equilibrium workflows that are now
considered canonical in `Thermodynamics/ChemEquilibrium`.

## Runnable guides

The following examples use the production prelude and bundled local NASA
thermochemistry. They are intentionally small, but each runs through the
same typed lookup, preparation, validation, and publication path as a larger
application.

| Scenario | Runnable example | Main result |
| --- | --- | --- |
| One `P,T` point | `chem_equilibrium_pt_point_guide.rs` | Accepted composition with lookup provenance. |
| `P,T` temperature range | `chem_equilibrium_pt_temperature_range_guide.rs` | Transactional range with accepted-state continuation. |
| One `P,H` point | `chem_equilibrium_ph_point_guide.rs` | Equilibrium composition and solved temperature. |
| `P,H` enthalpy range | `chem_equilibrium_ph_enthalpy_range_guide.rs` | Continued target-enthalpy sweep with per-point evidence. |

Run an example with `cargo run --example <example_name>`, for example:

```text
cargo run --example chem_equilibrium_ph_enthalpy_range_guide --no-default-features
```

## What lives where

- `ChemEquilibrium::legacy`
  - namespaced compatibility access to the prepared log-moles host and
    historical mutable workflows retained for migration
- `equilibrium_constant_problem.rs`
  - independent reaction-extent validation problem
- `equilibrium_constant_solver.rs`
  - small-system validator solver for `ln(Q) - ln(K)`
- `equilibrium_constant_validation.rs`
  - backend-independent validation reports
- `equilibrium_temperature_postprocessing.rs`
  - sweep resampling and report generation for plots / exports
- `phase_equilibrium_workflow.rs`
  - typed `ResolvedPhaseSystem -> MultiphaseEquilibriumSolution` facade for
    fixed-`P,T` systems with phase-qualified component identity
- `phase_equilibrium_workflow.rs` owns the production orchestration; the
  mutable phase-control helpers remain behind `ChemEquilibrium::legacy` and
  are not a second production entry point

## 1. Canonical gas-phase equilibrium solve

Use this when you want the production equilibrium path with solver fallback,
validation, and accepted-solution publication.

See:
- `src/Thermodynamics/ChemEquilibrium/phase_equilibrium_workflow.rs`
- `src/Thermodynamics/ChemEquilibrium/phase_equilibrium_problem.rs`
- `examples/chem_equilibrium_gas_example.rs`

```rust,ignore
use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumConditions, EquilibriumSolveOptions, PhaseEquilibriumPipelineRequest,
    SubstanceSystemSpecBuilder, SubstancesContainer,
};

let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
    "CO".to_string(), "CO2".to_string(), "O2".to_string(),
]))
.with_library_priorities(vec!["NASA_gas".to_string()])
.with_search_in_nist(false)
.build()?;

let outcome = PhaseEquilibriumPipelineRequest::new(
    spec,
    vec![0.25, 0.25, 0.5],
    EquilibriumConditions::new(1500.0, 101325.0, 101325.0)?,
)
.with_solve_options(EquilibriumSolveOptions::new().with_production_cascade())
.solve()?;

println!("lookup report: {:?}", outcome.lookup_report());
println!("{}", outcome.solution());
```

The older `ChemEquilibrium::legacy::gas_solver` helper remains available only
as a compatibility workflow. It owns mutable orchestration state and should
not be used by new production code.

Explicit backend policies do not require internal module imports:

```rust,ignore
use KiThe::Thermodynamics::ChemEquilibrium::prelude::{
    EquilibriumSolveOptions, LegacyEquilibriumSolver, RustedSciTheSolver,
    SolverBackend, SolverPolicy,
};

let rst_only = EquilibriumSolveOptions::new().with_solver_policy(
    SolverPolicy::Single(SolverBackend::RustedSciThe(
        RustedSciTheSolver::MinpackLevenbergMarquardt,
    )),
)?;

let legacy_fallback = SolverBackend::Legacy(LegacyEquilibriumSolver::NR);
```

The historical standalone solver and one-reaction helper live under
`ChemEquilibrium::legacy::nr` and `ChemEquilibrium::legacy::single_reaction`.

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

## 5. Resolved phase-system equilibrium

Use this boundary when thermochemistry has already been resolved through
`SubsData` and the calculation must retain phase-qualified identity and lookup
provenance. The one-shot facade is intentionally transactional: it neither
mutates `ResolvedPhaseSystem` nor exposes the historical mutable solver.

See:
- `src/Thermodynamics/ChemEquilibrium/phase_equilibrium_workflow.rs`
- `src/Thermodynamics/ChemEquilibrium/phase_equilibrium_solution.rs`
- `examples/chem_equilibrium_resolved_phase_example.rs`

```rust,ignore
let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())?;
let initial = MultiphaseInitialComposition::from_dense(&layout, vec![2.0, 1.0, 0.0])?;
let conditions = EquilibriumConditions::new(1200.0, 101_325.0, 101_325.0)?;

let solution = solve_resolved_pt(ResolvedPhaseEquilibriumRequest::new(
    &resolved,
    conditions,
    initial,
))?;
println!("{solution}");
```

## Recommended tests

- `src/Thermodynamics/ChemEquilibrium/equilibrium_log_moles_tests.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_constant_tests.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_constant_solver_tests.rs`
- `src/Thermodynamics/ChemEquilibrium/equilibrium_temperature_postprocessing.rs`
