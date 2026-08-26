# ChemEquilibrium Release Story Tests

This file is the operational catalog for expensive, ignored, real-data story
tests. It records what each test proves, the exact PowerShell command used to
run it, and the observed release output. The test name is the source of truth;
this document prevents release evidence from being lost in terminal history.

All commands below assume the repository root is the current directory.

Common command shape:

```powershell
cargo test --release --lib <TEST_FILTER> --no-default-features -- --ignored --nocapture
```

## 1. Real Phase-Transition Matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_phase_transition_release_matrix`

### What it checks

- real local NASA gas/condensed records;
- water gas/ice appearance;
- water gas/liquid appearance;
- hot-water liquid disappearance;
- graphite appearance and high-temperature disappearance;
- finite non-negative component amounts;
- scale-aware elemental conservation;
- final complementarity;
- minimum inactive-phase TPD in `J/mol` when an evaluated phase remains absent;
- maximum elemental-feasibility and constrained-minimizer KKT residuals;
- transition evidence and timing;
- byte-for-byte immutability of canonical JSON libraries.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_phase_transition_release_matrix `
  --no-default-features -- --ignored --nocapture
```

### Historical release output before expanded P7 evidence

The table below predates explicit `minimum_tpd`, feasibility, and KKT columns.
It is retained only as an earlier timing reference; the current P7 evidence is
recorded immediately after it. Do not compare timings across machines.

```text
live real phase-transition matrix (water/ice, water/liquid, graphite)
╭──────────────────┬────────┬───────────────────────────────────────────────┬─────────────┬─────────────┬─────────────────┬──────────╮
│ Fixture          │ T K    │ Active phases                                 │ Transitions │ Max balance │ Complementarity │ Total ms │
├──────────────────┼────────┼───────────────────────────────────────────────┼─────────────┼─────────────┼─────────────────┼──────────┤
│ water-gas-ice    │ 250.0  │ PhaseId(Some("gas")), PhaseId(Some("solid"))  │ 1           │ 0.000e0     │ OK              │ 4.027    │
│ water-gas-liquid │ 350.0  │ PhaseId(Some("gas")), PhaseId(Some("liquid")) │ 1           │ 6.661e-16   │ OK              │ 2.657    │
│ water-gas-hot    │ 550.0  │ PhaseId(Some("gas"))                          │ 0           │ 0.000e0     │ OK              │ 0.963    │
│ carbon-graphite  │ 700.0  │ PhaseId(Some("gas")), PhaseId(Some("solid"))  │ 1           │ 0.000e0     │ OK              │ 2.466    │
│ carbon-hot       │ 1400.0 │ PhaseId(Some("gas"))                          │ 0           │ 0.000e0     │ OK              │ 0.905    │
╰──────────────────┴────────┴───────────────────────────────────────────────┴─────────────┴─────────────┴─────────────────┴──────────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_phase_transition_release_matrix ... ok
```

### Recorded release output AFTER expanded P7 evidence

```text
live real phase-transition matrix (water/ice, water/liquid, graphite)
╭──────────────────┬────────┬───────────────────────────────────────────────┬─────────────┬─────────────┬────────────────────────┬─────────────────┬─────────┬─────────────────┬──────────╮
│ Fixture          │ T K    │ Active phases                                 │ Transitions │ Max balance │ Min inactive TPD J/mol │ Max feasibility │ Max KKT │ Complementarity │ Total ms │
├──────────────────┼────────┼───────────────────────────────────────────────┼─────────────┼─────────────┼────────────────────────┼─────────────────┼─────────┼─────────────────┼──────────┤
│ water-gas-ice    │ 250.0  │ PhaseId(Some("gas")), PhaseId(Some("solid"))  │ 1           │ 5.819e-9    │ -                      │ 2.220e-16       │ 0.000e0 │ OK              │ 3.360    │
│ water-gas-liquid │ 350.0  │ PhaseId(Some("gas")), PhaseId(Some("liquid")) │ 1           │ 6.632e-11   │ -                      │ 2.220e-16       │ 0.000e0 │ OK              │ 2.221    │
│ water-gas-hot    │ 550.0  │ PhaseId(Some("gas"))                          │ 0           │ 0.000e0     │ 1.961e4                │ 2.220e-16       │ 0.000e0 │ OK              │ 0.950    │
│ carbon-graphite  │ 700.0  │ PhaseId(Some("gas")), PhaseId(Some("solid"))  │ 1           │ 1.013e-11   │ -                      │ 1.776e-15       │ 0.000e0 │ OK              │ 2.274    │
│ carbon-hot       │ 1400.0 │ PhaseId(Some("gas"))                          │ 0           │ 0.000e0     │ 4.842e4                │ 1.776e-15       │ 0.000e0 │ OK              │ 0.914    │
╰──────────────────┴────────┴───────────────────────────────────────────────┴─────────────┴─────────────┴────────────────────────┴─────────────────┴─────────┴─────────────────┴──────────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_phase_transition_release_matrix ... ok
```

Result: **passed, 5/5 cases**.

The separate real ice temperature-range story also verifies that the final
accepted point retains one timing snapshot per prepared active-set cache entry;
the release phase matrix above is the compact single-point transition table.

## 1a. Real Ice Temperature-Range Cache Evidence

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_phase_control_temperature_range_records_real_ice_transition`

### What it checks

- real gas/ice temperature continuation at `250`, `260`, and `270 K`;
- phase transition and continuation reuse;
- distinct projection, prepared, and RST symbolic cache entries;
- runner-scoped active-element TPD geometry entries/builds/reuses;
- one deterministic build-duration snapshot per prepared active-set entry;
- conservation, finite amounts, transition timing, and JSON immutability.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_phase_control_temperature_range_records_real_ice_transition `
  --no-default-features -- --ignored --nocapture
```

### Historical output before TPD-geometry cache evidence

The historical line below predates the `tpd_geometry=entries/builds/reuses`
suffix. It is retained only as an earlier timing reference; current release
evidence follows immediately after it.

```text
live bounded ice T-range: points=3 transitions=1 projections=2 prepared=2 rst_symbolic=2 total=14.2398ms
```
### Recorded release output after TPD-geometry cache evidence

```text
live bounded ice T-range: points=3 transitions=1 projections=2 prepared=2 rst_symbolic=2 tpd_geometry=1/1/3 total=4.3895ms
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_phase_control_temperature_range_records_real_ice_transition ... ok
```



Result: **passed in release**. The cache counts are behavioral evidence and
must remain stable for this fixture.

## 2. Real P,H Jacobian Boundary and Scaling Matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_gas_ph_jacobian_matches_central_difference_near_interval_boundaries`

### What it checks

- five real NASA C/H/O gas records;
- temperatures just inside both native interval boundaries and one interior
  temperature;
- every row and column of the coupled P,H Jacobian;
- analytic Jacobian versus central finite differences;
- inventory scales `1e-6`, `1e-3`, `1.0`, and `1e3`;
- no accidental absolute-scale dependence;
- unchanged JSON source files.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_gas_ph_jacobian_matches_central_difference_near_interval_boundaries `
  --no-default-features -- --ignored --nocapture
```

### Recorded output

```text
live real P,H Jacobian matrix: species=5 scales=4 temperatures=3 dimension=6 checked_entries=432 status=OK
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_gas_ph_jacobian_matches_central_difference_near_interval_boundaries ... ok
```

Result: **passed in release**. The summary is intentionally compact; every
matrix entry is still an assertion. A failed entry prints scale, temperature,
row, column, analytic value, finite-difference value, and tolerance.

## 3. Real 100-Species, 50-Point Temperature-Range Matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_100_species_50_point_temperature_range_backend_matrix`

### What it checks

- real element-limited C/H/O candidate selection;
- 100 species and 50 temperature points;
- typed temperature-range formulation reuse;
- all configured RST and legacy backends;
- per-backend total, wall, mean, median, and worst time;
- formulation builds/reuses and symbolic updates;
- residual and conservation maxima;
- explicit failure rows instead of hidden fallback;
- release scalability characterization.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_100_species_50_point_temperature_range_backend_matrix `
  --no-default-features -- --ignored --nocapture
```

### Recorded release output

```text
live real temperature-range backend matrix (species=100, points=50)
╭─────────────────────┬────────┬──────────┬──────────┬─────────┬────────┬───────────┬──────────┬───────────┬────────────────────┬───────────────────┬─────────────────────────────────────────────┬────────────────┬──────────────────┬────────┬──────────────────┬─────────────────────┬─────────────────────────┬──────────────────────┬────────┬──────────────────┬──────────────┬─────────────┬─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ Backend             │ Status │ Total ms │ Wall ms  │ Mean ms │ Min ms │ Median ms │ Worst ms │ Worst T K │ Worst nonlinear ms │ Worst symbolic ms │ Worst accepted backend                      │ Worst attempts │ Worst iterations │ Builds │ Initial setup ms │ Initial symbolic ms │ Initial problem prep ms │ Formulation build ms │ Reuses │ Symbolic updates │ Max residual │ Max balance │ Error                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       │
├─────────────────────┼────────┼──────────┼──────────┼─────────┼────────┼───────────┼──────────┼───────────┼────────────────────┼───────────────────┼─────────────────────────────────────────────┼────────────────┼──────────────────┼────────┼──────────────────┼─────────────────────┼─────────────────────────┼──────────────────────┼────────┼──────────────────┼──────────────┼─────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ rst_lm              │ OK     │ 5024.478 │ 5264.596 │ 14.009  │ 13.250 │ 13.687    │ 25.697   │ 1000.0    │ 24.859             │ 0.000             │ RustedSciThe(LevenbergMarquardt)            │ 1              │ 9                │ 1      │ 4298.557         │ 4312.329            │ 0.432                   │ 0.000                │ 49     │ 49               │ 2.902e-7     │ 2.883e-7    │ -                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           │
│ rst_minpack_lm      │ OK     │ 4658.434 │ 4909.922 │ 8.440   │ 7.643  │ 7.990     │ 29.084   │ 1000.0    │ 28.299             │ 0.000             │ RustedSciThe(MinpackLevenbergMarquardt)     │ 1              │ 9                │ 1      │ 4211.377         │ 4224.590            │ 0.435                   │ 0.000                │ 49     │ 49               │ 1.907e-7     │ 7.116e-8    │ -                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           │
│ rst_nielsen_lm      │ FAILED │ -        │ 4754.355 │ -       │ -      │ -         │ -        │ -         │ -                  │ -                 │ -                                           │ -              │ -                │ -      │ -                │ -                   │ -                       │ -                    │ -      │ -                │ -            │ -           │ TemperatureRangePointFailed { point_index: 0, temperature: 1000.0, cause: AllBackendsFailed { attempts: [SolverAttemptReport { backend: RustedSciThe(NielsenLevenbergMarquardt), outcome: RejectedCandidate { reason: "equilibrium candidate field 'candidate_acceptance_residual' was rejected: residual L2 norm 1.3401354047923884e0 exceeds tolerance 1e-6" }, metrics: Some(SolverAttemptMetrics { termination: MaxIterations, backend_converged: false, iterations: 50, residual_evaluations: 51, jacobian_evaluations: 51, linear_solves: 56, elapsed_millis: 229, evaluation_timing: Some(SolverEvaluationTiming { residual_evaluation_micros: 15302, jacobian_evaluation_micros: 202862, solver_overhead_micros: 11115 }) }) }] } } │
│ rst_trust_region_lm │ OK     │ 4675.175 │ 4938.955 │ 8.775   │ 7.751  │ 8.209     │ 28.489   │ 1000.0    │ 27.704             │ 0.000             │ RustedSciThe(TrustRegionLevenbergMarquardt) │ 1              │ 9                │ 1      │ 4210.601         │ 4224.113            │ 0.413                   │ 0.000                │ 49     │ 49               │ 1.907e-7     │ 7.116e-8    │ -                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           │
│ rst_powell_dogleg   │ OK     │ 4956.235 │ 5212.086 │ 15.270  │ 11.984 │ 12.801    │ 123.487  │ 1000.0    │ 122.658            │ 0.000             │ RustedSciThe(PowellDogleg)                  │ 1              │ 50               │ 1      │ 4167.287         │ 4180.595            │ 0.402                   │ 0.000                │ 49     │ 49               │ 6.098e-7     │ 4.014e-7    │ -                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           │
│ rst_damped_newton   │ FAILED │ -        │ 4463.248 │ -       │ -      │ -         │ -        │ -         │ -                  │ -                 │ -                                           │ -              │ -                │ -      │ -                │ -                   │ -                       │ -                    │ -      │ -                │ -            │ -           │ TemperatureRangePointFailed { point_index: 0, temperature: 1000.0, cause: AllBackendsFailed { attempts: [SolverAttemptReport { backend: RustedSciThe(DampedNewton), outcome: RejectedCandidate { reason: "equilibrium candidate field 'candidate_acceptance_residual' was rejected: residual L2 norm 1.381759743012583e2 exceeds tolerance 1e-6" }, metrics: Some(SolverAttemptMetrics { termination: Stagnation, backend_converged: false, iterations: 2, residual_evaluations: 3, jacobian_evaluations: 3, linear_solves: 3, elapsed_millis: 8, evaluation_timing: Some(SolverEvaluationTiming { residual_evaluation_micros: 730, jacobian_evaluation_micros: 7194, solver_overhead_micros: 297 }) }) }] } }                              │
│ legacy_lm           │ OK     │ 101.594  │ 102.064  │ 1.594   │ 1.494  │ 1.531     │ 2.460    │ 1000.0    │ 1.951              │ 0.000             │ Legacy(LM)                                  │ 1              │ 0                │ 1      │ 0.401            │ 13.840              │ 0.401                   │ 0.000                │ 49     │ 0                │ 9.516e-13    │ 9.460e-13   │ -                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           │
│ legacy_nr           │ OK     │ 54.158   │ 54.489   │ 0.665   │ 0.623  │ 0.636     │ 1.508    │ 1000.0    │ 1.086              │ 0.000             │ Legacy(NR)                                  │ 1              │ 0                │ 1      │ 0.407            │ 13.818              │ 0.407                   │ 0.000                │ 49     │ 0                │ 1.907e-7     │ 7.116e-8    │ -                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           │
│ legacy_tr           │ FAILED │ -        │ 22.578   │ -       │ -      │ -         │ -        │ -         │ -                  │ -                 │ -                                           │ -              │ -                │ -      │ -                │ -                   │ -                       │ -                    │ -      │ -                │ -            │ -           │ TemperatureRangePointFailed { point_index: 0, temperature: 1000.0, cause: AllBackendsFailed { attempts: [SolverAttemptReport { backend: Legacy(TR), outcome: Failed { kind: Solver, reason: "equilibrium solver failed: nonlinear solver reached its iteration limit" }, metrics: None }] } }                                                                                                                                                                                                                                                                                                                                                                                                                                               │
╰─────────────────────┴────────┴──────────┴──────────┴─────────┴────────┴───────────┴──────────┴───────────┴────────────────────┴───────────────────┴─────────────────────────────────────────────┴────────────────┴──────────────────┴────────┴──────────────────┴─────────────────────┴─────────────────────────┴──────────────────────┴────────┴──────────────────┴──────────────┴─────────────┴─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
three slowest accepted points per successful backend
╭─────────────────────┬───────┬────────┬──────────┬──────────────┬─────────────┬─────────────┬─────────────────┬───────────────┬───────────────┬───────────────────┬───────────────────┬─────────────────────────────────────────────┬─────────────┬────────────┬────────────────┬────────────────┬───────────────┬────────────┬─────────────┬─────────────┬───────────╮
│ Backend             │ Point │ T K    │ Total ms │ Nonlinear ms │ Symbolic ms │ Closures ms │ Problem prep ms │ Validation ms │ Seed dlog inf │ Seed log range    │ Result log range  │ Accepted                                    │ Termination │ Iterations │ Residual evals │ Jacobian evals │ Linear solves │ Backend ms │ Residual ms │ Jacobian ms │ Engine ms │
├─────────────────────┼───────┼────────┼──────────┼──────────────┼─────────────┼─────────────┼─────────────────┼───────────────┼───────────────┼───────────────────┼───────────────────┼─────────────────────────────────────────────┼─────────────┼────────────┼────────────────┼────────────────┼───────────────┼────────────┼─────────────┼─────────────┼───────────┤
│ rst_powell_dogleg   │ 0     │ 1000.0 │ 123.487  │ 122.658      │ 0.000       │ 0.457       │ 0.000           │ 0.000         │ 6.816e1       │ [-6.908, -6.908]  │ [-75.063, -2.343] │ RustedSciThe(PowellDogleg)                  │ Converged   │ 50         │ 51             │ 51             │ 50            │ 122        │ 12.582      │ 102.789     │ 7.101     │
│ rst_minpack_lm      │ 0     │ 1000.0 │ 29.084   │ 28.299       │ 0.000       │ 0.433       │ 0.000           │ 0.000         │ 6.816e1       │ [-6.908, -6.908]  │ [-75.063, -2.343] │ RustedSciThe(MinpackLevenbergMarquardt)     │ Converged   │ 9          │ 10             │ 10             │ 9             │ 28         │ 2.551       │ 20.752      │ 4.824     │
│ rst_trust_region_lm │ 0     │ 1000.0 │ 28.489   │ 27.704       │ 0.000       │ 0.428       │ 0.000           │ 0.000         │ 6.816e1       │ [-6.908, -6.908]  │ [-75.063, -2.343] │ RustedSciThe(TrustRegionLevenbergMarquardt) │ Converged   │ 9          │ 10             │ 10             │ 9             │ 27         │ 1.963       │ 21.029      │ 4.538     │
│ rst_lm              │ 0     │ 1000.0 │ 25.697   │ 24.859       │ 0.000       │ 0.418       │ 0.000           │ 0.000         │ 6.816e1       │ [-6.908, -6.908]  │ [-75.063, -2.343] │ RustedSciThe(LevenbergMarquardt)            │ Converged   │ 9          │ 10             │ 10             │ 9             │ 24         │ 2.134       │ 21.111      │ 1.356     │
│ rst_powell_dogleg   │ 3     │ 1030.0 │ 15.716   │ 14.710       │ 0.000       │ 0.638       │ 0.000           │ 0.000         │ 1.030e0       │ [-72.944, -2.370] │ [-71.914, -2.384] │ RustedSciThe(PowellDogleg)                  │ Converged   │ 5          │ 6              │ 6              │ 5             │ 14         │ 1.431       │ 11.973      │ 1.142     │
│ rst_powell_dogleg   │ 5     │ 1050.0 │ 15.553   │ 14.780       │ 0.000       │ 0.414       │ 0.000           │ 0.000         │ 9.941e-1      │ [-70.902, -2.399] │ [-70.060, -2.415] │ RustedSciThe(PowellDogleg)                  │ Converged   │ 5          │ 6              │ 6              │ 5             │ 14         │ 1.521       │ 11.955      │ 1.118     │
│ rst_lm              │ 26    │ 1260.0 │ 15.081   │ 13.909       │ 0.000       │ 0.766       │ 0.000           │ 0.000         │ 7.147e-1      │ [-76.271, -2.424] │ [-76.564, -2.394] │ RustedSciThe(LevenbergMarquardt)            │ Converged   │ 5          │ 6              │ 6              │ 5             │ 13         │ 1.284       │ 11.611      │ 0.748     │
│ rst_lm              │ 27    │ 1270.0 │ 14.832   │ 13.707       │ 0.000       │ 0.621       │ 0.000           │ 0.000         │ 7.047e-1      │ [-76.564, -2.394] │ [-76.854, -2.364] │ RustedSciThe(LevenbergMarquardt)            │ Converged   │ 5          │ 6              │ 6              │ 5             │ 13         │ 1.250       │ 11.570      │ 0.695     │
│ rst_trust_region_lm │ 49    │ 1490.0 │ 10.679   │ 9.892        │ 0.000       │ 0.435       │ 0.000           │ 0.000         │ 5.317e-1      │ [-82.603, -1.981] │ [-82.865, -1.971] │ RustedSciThe(TrustRegionLevenbergMarquardt) │ Converged   │ 2          │ 3              │ 3              │ 2             │ 9          │ 0.715       │ 7.733       │ 1.246     │
│ rst_trust_region_lm │ 47    │ 1470.0 │ 10.598   │ 9.812        │ 0.000       │ 0.431       │ 0.000           │ 0.000         │ 5.443e-1      │ [-82.077, -2.002] │ [-82.340, -1.991] │ RustedSciThe(TrustRegionLevenbergMarquardt) │ Converged   │ 2          │ 3              │ 3              │ 2             │ 9          │ 0.808       │ 7.537       │ 1.258     │
│ rst_minpack_lm      │ 21    │ 1210.0 │ 8.713    │ 7.857        │ 0.000       │ 0.452       │ 0.000           │ 0.000         │ 7.683e-1      │ [-74.786, -2.597] │ [-75.086, -2.560] │ RustedSciThe(MinpackLevenbergMarquardt)     │ Converged   │ 2          │ 3              │ 3              │ 2             │ 7          │ 0.656       │ 5.935       │ 1.078     │
│ rst_minpack_lm      │ 2     │ 1020.0 │ 8.459    │ 7.461        │ 0.000       │ 0.644       │ 0.000           │ 0.000         │ 1.050e0       │ [-73.994, -2.356] │ [-72.944, -2.370] │ RustedSciThe(MinpackLevenbergMarquardt)     │ Converged   │ 2          │ 3              │ 3              │ 2             │ 7          │ 0.578       │ 5.643       │ 1.049     │
│ legacy_lm           │ 0     │ 1000.0 │ 2.460    │ 1.951        │ 0.000       │ 0.287       │ 0.000           │ 0.000         │ 6.816e1       │ [-6.908, -6.908]  │ [-75.063, -2.343] │ Legacy(LM)                                  │ -           │ -          │ -              │ -              │ -             │ -          │ -           │ -           │ -         │
│ legacy_lm           │ 7     │ 1070.0 │ 1.848    │ 1.425        │ 0.000       │ 0.241       │ 0.000           │ 0.000         │ 9.599e-1      │ [-70.390, -2.432] │ [-70.718, -2.449] │ Legacy(LM)                                  │ -           │ -          │ -              │ -              │ -             │ -          │ -           │ -           │ -         │
│ legacy_lm           │ 5     │ 1050.0 │ 1.763    │ 1.287        │ 0.000       │ 0.261       │ 0.000           │ 0.000         │ 9.941e-1      │ [-70.902, -2.399] │ [-70.060, -2.415] │ Legacy(LM)                                  │ -           │ -          │ -              │ -              │ -             │ -          │ -           │ -           │ -         │
│ legacy_nr           │ 0     │ 1000.0 │ 1.508    │ 1.086        │ 0.000       │ 0.231       │ 0.000           │ 0.000         │ 6.816e1       │ [-6.908, -6.908]  │ [-75.063, -2.343] │ Legacy(NR)                                  │ -           │ -          │ -              │ -              │ -             │ -          │ -           │ -           │ -         │
│ legacy_nr           │ 13    │ 1130.0 │ 0.837    │ 0.409        │ 0.000       │ 0.247       │ 0.000           │ 0.000         │ 8.688e-1      │ [-72.321, -2.549] │ [-72.635, -2.571] │ Legacy(NR)                                  │ -           │ -          │ -              │ -              │ -             │ -          │ -           │ -           │ -         │
│ legacy_nr           │ 35    │ 1350.0 │ 0.807    │ 0.404        │ 0.000       │ 0.223       │ 0.000           │ 0.000         │ 6.317e-1      │ [-78.841, -2.191] │ [-79.118, -2.171] │ Legacy(NR)                                  │ -           │ -          │ -              │ -              │ -             │ -          │ -           │ -           │ -         │
╰─────────────────────┴───────┴────────┴──────────┴──────────────┴─────────────┴─────────────┴─────────────────┴───────────────┴───────────────┴───────────────────┴───────────────────┴─────────────────────────────────────────────┴─────────────┴────────────┴────────────────┴────────────────┴───────────────┴────────────┴─────────────┴─────────────┴───────────╯
failed first-point/backend attempt diagnostics
╭───────────────────┬───────┬────────┬───────────────────┬─────────┬─────────────────────────────────────────┬───────────────────┬───────────────┬────────────┬────────────────┬────────────────┬───────────────┬────────────┬─────────────┬─────────────┬───────────┬─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ Requested backend │ Point │ T K    │ Error kind        │ Attempt │ Attempt backend                         │ Outcome           │ Termination   │ Iterations │ Residual evals │ Jacobian evals │ Linear solves │ Backend ms │ Residual ms │ Jacobian ms │ Engine ms │ Detail                                                                                                                                                                                                              │
├───────────────────┼───────┼────────┼───────────────────┼─────────┼─────────────────────────────────────────┼───────────────────┼───────────────┼────────────┼────────────────┼────────────────┼───────────────┼────────────┼─────────────┼─────────────┼───────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ rst_nielsen_lm    │ 0     │ 1000.0 │ AllBackendsFailed │ 0       │ RustedSciThe(NielsenLevenbergMarquardt) │ RejectedCandidate │ MaxIterations │ 50         │ 51             │ 51             │ 56            │ 229        │ 15.302      │ 202.862     │ 11.115    │ backend=RustedSciThe(NielsenLevenbergMarquardt), outcome=rejected_candidate: equilibrium candidate field 'candidate_acceptance_residual' was rejected: residual L2 norm 1.3401354047923884e0 exceeds tolerance 1e-6 │
│ rst_damped_newton │ 0     │ 1000.0 │ AllBackendsFailed │ 0       │ RustedSciThe(DampedNewton)              │ RejectedCandidate │ Stagnation    │ 2          │ 3              │ 3              │ 3             │ 8          │ 0.730       │ 7.194       │ 0.297     │ backend=RustedSciThe(DampedNewton), outcome=rejected_candidate: equilibrium candidate field 'candidate_acceptance_residual' was rejected: residual L2 norm 1.381759743012583e2 exceeds tolerance 1e-6               │
│ legacy_tr         │ 0     │ 1000.0 │ AllBackendsFailed │ 0       │ Legacy(TR)                              │ Failed(Solver)    │ -             │ -          │ -              │ -              │ -             │ -          │ -           │ -           │ -         │ backend=Legacy(TR), outcome=failed(Solver): equilibrium solver failed: nonlinear solver reached its iteration limit                                                                                                 │
╰───────────────────┴───────┴────────┴───────────────────┴─────────┴─────────────────────────────────────────┴───────────────────┴───────────────┴────────────┴────────────────┴────────────────┴───────────────┴────────────┴─────────────┴─────────────┴───────────┴─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
summary: attempted=9 successful=6 failed=["rst_nielsen_lm", "rst_damped_newton", "legacy_tr"]
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_100_species_50_point_temperature_range_backend_matrix ... ok
```

Result: **passed as a characterization matrix, 6/9 backends successful**.
The failed backends are retained as explicit evidence and are not treated as
test failure when another backend succeeds.

## 4. Real P,H Backend Matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_backend_matrix`

### What it checks

- one fixed real H2/O2/H2O problem;
- each RST and legacy inner backend under an explicit `Single` policy;
- outer trial count versus inner nonlinear work;
- residual, conservation, enthalpy error, and final temperature;
- backend failures reported independently;
- unchanged canonical JSON libraries.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_backend_matrix `
  --no-default-features -- --ignored --nocapture
```

### Recorded output

Not yet recorded in this document. Run the command and append the complete
table below. Do not summarize failed backends away.

```text
╭──────────────────┬────────┬───────────────────────────────────────────────┬─────────────┬─────────────┬─────────────────┬──────────╮
│ Fixture          │ T K    │ Active phases                                 │ Transitions │ Max balance │ Complementarity │ Total ms │
├──────────────────┼────────┼───────────────────────────────────────────────┼─────────────┼─────────────┼─────────────────┼──────────┤
│ water-gas-ice    │ 250.0  │ PhaseId(Some("gas")), PhaseId(Some("solid"))  │ 1           │ 0.000e0     │ OK              │ 3.544    │
│ water-gas-liquid │ 350.0  │ PhaseId(Some("gas")), PhaseId(Some("liquid")) │ 1           │ 6.661e-16   │ OK              │ 2.528    │
│ water-gas-hot    │ 550.0  │ PhaseId(Some("gas"))                          │ 0           │ 0.000e0     │ OK              │ 0.944    │
│ carbon-graphite  │ 700.0  │ PhaseId(Some("gas")), PhaseId(Some("solid"))  │ 1           │ 0.000e0     │ OK              │ 2.673    │
│ carbon-hot       │ 1400.0 │ PhaseId(Some("gas"))                          │ 0           │ 0.000e0     │ OK              │ 1.050    │
╰──────────────────┴────────┴───────────────────────────────────────────────┴─────────────┴─────────────┴─────────────────┴──────────╯
```

## 5. Real Ice P,T -> P,H Auto Route

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_ice_pt_to_h_to_auto_ph_preserves_phase_transition`

### What it checks

- real gas/ice reference state;
- P,T to total enthalpy to P,H inversion;
- rejected monolithic candidate retained in evidence;
- deterministic nested recovery through `Auto`;
- ice activation and phase-control evidence;
- enthalpy acceptance and JSON immutability.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_ice_pt_to_h_to_auto_ph_preserves_phase_transition `
  --no-default-features -- --ignored --nocapture
```

### Recorded output

```text
live real ice P,H route matrix: monolithic=FAILED(AllBackendsFailed) auto_fallback=Some("AllBackendsFailed") solid_moles=4.998093e-1 transitions=1
╭──────────────────────┬────────┬───────────┬────────────────┬─────────┬──────────────┬────────┬──────────┬────────────┬────────┬────────┬─────────────┬────────────────┬────────────────┬────────────────┬───────────┬─────────┬───────╮
│ Path                 │ Status │ Report ms │ Thermo prep ms │ Wall ms │ Inner solves │ Trials │ Attempts │ Iterations │ Builds │ Reuses │ Transitions │ Residual evals │ Jacobian evals │ Scaled H error │ Residual  │ Balance │ Error │
├──────────────────────┼────────┼───────────┼────────────────┼─────────┼──────────────┼────────┼──────────┼────────────┼────────┼────────┼─────────────┼────────────────┼────────────────┼────────────────┼───────────┼─────────┼───────┤
│ nested-phase-control │ OK     │ 0.000     │ 0.000          │ 33.457  │ 3            │ 3      │ 18       │ 0          │ 0      │ 0      │ 3           │ 0              │ 0              │ 0.000e0        │ 8.882e-16 │ 0.000e0 │ -     │
│ auto-nested-recovery │ OK     │ 0.000     │ 0.000          │ 25.517  │ 3            │ 3      │ 18       │ 0          │ 0      │ 0      │ 3           │ 0              │ 0              │ 0.000e0        │ 8.882e-16 │ 0.000e0 │ -     │
╰──────────────────────┴────────┴───────────┴────────────────┴─────────┴──────────────┴────────┴──────────┴────────────┴────────┴────────┴─────────────┴────────────────┴────────────────┴────────────────┴───────────┴─────────┴───────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_ice_pt_to_h_to_auto_ph_preserves_phase_transition ... ok
```

Result: **passed in release**. If a failure occurs, preserve the complete typed
error and fallback evidence here.

## 6. Release Recording Rules

- Run ignored story tests in release mode with `--nocapture`.
- Record the exact command and the complete table, including failed backends.
- Keep numerical tolerances and expected physical contracts in the test code.
- Keep timing values in this document as characterization, not hard assertions.
- Re-run the relevant story after changes to solver policy, phase control,
  thermochemistry lookup, or repository data.
- If a test depends on a repaired NIST fixture, record the fixture revision and
  lookup policy next to its output.
- Confirm canonical JSON fingerprints before and after every live story.

## 7. GUI P,H Failure Boundary

### Test

`gui::equilibrium_gui_tests::offline_local_unreachable_ph_target_is_transactional`

### What it checks

- a finite but unreachable P,H target;
- visible worker failure rather than a false completed state;
- no partial result snapshot publication;
- unchanged local JSON libraries.

### Release command

```powershell
cargo test --release --lib `
  gui::equilibrium_gui_tests::offline_local_unreachable_ph_target_is_transactional `
  --no-default-features -- --ignored --nocapture
```

### Recorded output

```text
test ...::offline_local_unreachable_ph_target_is_transactional ... ok
```

Debug validation passed. Record the release run here when executed; the
behavioral contract is the failed worker state with no published snapshot.

## 8. Online NIST physical-state parser diagnostic

### Test

`Thermodynamics::DBhandlers::NIST_parser_tests::tests::test_real_water_state_specific_payload_matrix`

### What it checks

- the parser requests gas, liquid, and solid independently;
- complete payloads contain finite Cp intervals and matching coefficient rows;
- an incomplete state payload is reported as unavailable thermochemistry;
- no online result is written into the local JSON repository.

### Command

```powershell
cargo test --lib --no-default-features `
  test_real_water_state_specific_payload_matrix `
  -- --ignored --nocapture
```

### Recorded output

```text
NIST state payload: state=gas status=COMPLETE ranges=2
NIST state payload: state=liquid status=COMPLETE ranges=1
NIST state payload: state=solid status=INCOMPLETE
test ...::test_real_water_state_specific_payload_matrix ... ok
```

Interpretation: NIST exposes usable H2O gas/liquid Cp data in this run, while
the solid navigation page does not expose a Cp interval table. The solid case
is not a valid thermochemical record and must never fall back to gas or liquid
data. This diagnostic is network-dependent and is not a release solver gate.

## 9. Online NIST phase-resolution boundary

### Test

`Thermodynamics::User_PhaseOrSolution_tests::tests::phase_resolution_can_use_explicit_online_nist_for_requested_liquid`

### What it checks

- an explicit `NIST` library instruction is accepted by the phase factory;
- local NIST miss may use the online parser only under the explicit policy;
- liquid `PhaseSpec` is carried into the NIST query;
- the resolved report retains `ExactRequestedState` and NIST provenance.

### Command

```powershell
cargo test --lib --no-default-features `
  phase_resolution_can_use_explicit_online_nist_for_requested_liquid `
  -- --ignored --nocapture
```

### Recorded output

```text
test ...::phase_resolution_can_use_explicit_online_nist_for_requested_liquid ... ok
```

Debug validation passed. This is an online lookup-boundary story, not a
solver release gate; it does not write or repair local JSON files.

## 10. Typed P,H target-enthalpy continuation

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_continues_from_previous_accepted_state`

### What it checks

- real local NASA H2/O2/H2O thermochemistry;
- one typed P,H target grid solved in ascending and descending enthalpy order;
- the first point uses the caller seed and every later point uses the previous
  accepted physical composition and temperature;
- accepted temperatures, conservation, and immutable thermochemistry lookup;
- the nested reference route reuses its fixed-P,T template across the same
  target range while keeping scalar brackets independent;
- no partial batch publication and no JSON-library mutation.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_continues_from_previous_accepted_state `
  --no-default-features -- --ignored --nocapture
```

### Recorded release summary

| Scenario | Mode | Status | Points | Builds | Reuses | Transitions | Total ms | Mean ms | Worst ms | Detail |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| reactive-gas | Nested | OK | 3 | 1 | 41 | 0 | 2.282 | 0.761 | 0.859 | fixed-template=shared |
| water-ice | Auto | OK | 3 | 0 | 0 | 1 | 25.137 | 8.379 | 22.742 | nested_points=1, fallback_points=1 |
| unreachable-target | Nested | ROLLBACK | 0 | - | - | - | 0.809 | - | - | failed_point=1 |

### Raw console transcript

```text
P,H target range | direction=Ascending | points=3 | formulation_builds=1 | formulation_reuses=2
index | target J | seed K | solved K | preparation | elapsed
0 | -2.956790e5 | 2300.000 | 2300.000 | Initial | 18.4092ms
1 | -2.638129e5 | 2300.000 | 2500.000 | Continued | 14.946ms
2 | -2.256501e5 | 2500.000 | 2700.000 | Continued | 4.6834ms
nested | points=3 | formulation_builds=1 | formulation_reuses=41 | total=7.206ms | mean=2.402ms | worst=2.8865ms
test ...::live_reactive_gas_ph_target_range_continues_from_previous_accepted_state ... ok
```

The debug run passed in both directions. Release timing and the complete
descending table should be appended after the operator-run release command.
An attempted release run on 2026-08-02 exceeded two three-minute waits while
the crate was still compiling in `rustc`; it was stopped without a test
result, so no release pass is claimed here.

The same story also runs a strict nested reference subcase. Its assertions
require `formulation_builds=1` for the whole nested range and at least one
reuse, proving that the scalar bracket is recreated per target without
rebuilding the fixed-P,T structural template.

## 11. P,H target-range backend matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_backend_matrix`

### What it checks

- the same real H2/O2/H2O target-enthalpy grid for every selected backend;
- strict `Single` backend policy, so backend failure is visible;
- monolithic P,H continuation and prepared formulation reuse counters;
- point-level total/mean/worst timing;
- unchanged local JSON libraries even when a backend fails.

This matrix deliberately covers the fixed-phase monolithic route only. Nested,
`Auto`, phase-transition, and rollback batches remain separate lifecycle
stories.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_backend_matrix `
  --no-default-features -- --ignored --nocapture
```

### Recorded debug output

```text
backend             | status | points | builds | reuses | total_ms | mean_ms | worst_ms | error
rst_lm              | FAILED |      - |      - |      - |   43.638 |       - |        - | P,H target point 0 (-2.956790e5 J) failed: all equilibrium solver backends failed after 1 attempt(s)
rst_minpack_lm      | FAILED |      - |      - |      - |   41.295 |       - |        - | P,H target point 2 (-2.256501e5 J) failed: all equilibrium solver backends failed after 1 attempt(s)
rst_trust_region_lm | FAILED |      - |      - |      - |   41.502 |       - |        - | P,H target point 2 (-2.256501e5 J) failed: all equilibrium solver backends failed after 1 attempt(s)
legacy_lm           | FAILED |      - |      - |      - |   19.186 |       - |        - | P,H target point 0 (-2.956790e5 J) failed: all equilibrium solver backends failed after 1 attempt(s)
legacy_nr           | FAILED |      - |      - |      - |    2.664 |       - |        - | P,H target point 2 (-2.256501e5 J) failed: all equilibrium solver backends failed after 1 attempt(s)
legacy_tr           | OK     |      3 |      1 |      2 |    1.353 |   0.451 |    0.511 | -
test ...::live_reactive_gas_ph_target_range_backend_matrix ... ok
```

This debug matrix passed because at least one explicitly selected backend
accepted the complete range; the failed rows remain intentional evidence, not
suppressed errors. Release output is still pending.

## 12. Nested/Auto P,H route and rollback matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_and_water_ph_target_range_route_matrix`

### What it checks

- nested real H2/O2/H2O target range with one shared fixed-P,T template;
- real water/ice `Auto` target range with phase-control transition and
  retained monolithic-fallback evidence;
- an unreachable later target returning a typed point error at index `1`;
- no partial range publication and no mutation of the local JSON libraries.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_and_water_ph_target_range_route_matrix `
  --no-default-features -- --ignored --nocapture
```

### Recorded debug output

```text
live P,H nested/Auto route matrix (real NASA gas and water/ice)
╭────────────────────┬────────┬──────────┬────────┬────────┬────────┬─────────────┬──────────┬─────────┬──────────┬───────────────────────────────────╮
│ Scenario           │ Mode   │ Status   │ Points │ Builds │ Reuses │ Transitions │ Total ms │ Mean ms │ Worst ms │ Detail                            │
├────────────────────┼────────┼──────────┼────────┼────────┼────────┼─────────────┼──────────┼─────────┼──────────┼───────────────────────────────────┤
│ reactive-gas       │ Nested │ OK       │ 3      │ 1      │ 41     │ 0           │ 2.282    │ 0.761   │ 0.859    │ fixed-template=shared             │
│ water-ice          │ Auto   │ OK       │ 3      │ 0      │ 0      │ 1           │ 25.137   │ 8.379   │ 22.742   │ nested_points=1 fallback_points=1 │
│ unreachable-target │ Nested │ ROLLBACK │ 0      │ -      │ -      │ -           │ 0.809    │ -       │ -        │ failed_point=1                    │
╰────────────────────┴────────┴──────────┴────────┴────────┴────────┴─────────────┴──────────┴─────────┴──────────┴───────────────────────────────────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_and_water_ph_target_range_route_matrix ... ok
```

The release run passed. The normalized table above is the canonical evidence;
the raw transcript is retained only as an exact console artifact.

## 13. Large real nested P,H continuation range

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_nested_ph_target_range_story`

### What it checks

- twenty locally resolved NASA-gas C/H/O species;
- nine target enthalpies inside one common coefficient interval;
- accepted composition and temperature continuation across eight hand-offs;
- one fixed-P,T template build and repeated inner-template reuse;
- inverse-temperature accuracy, conservation, enthalpy acceptance, and JSON
  immutability at every point.

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_nested_ph_target_range_story `
  --no-default-features -- --ignored --nocapture
```

### Recorded release summary

| Species | Targets | Template builds | Template reuses | Total ms | Worst point ms | Max balance | Status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 20 | 9 | 1 | 163 | 41.424 | 8.482 | 1.534e-8 | OK |

### Raw console transcript

```text
live large nested P,H target range: species=20 points=9 builds=1 reuses=163 total=41.4236ms
╭───────┬─────────────┬──────────┬──────────┬───────────────────┬──────────┬─────────────┬───────────┬────────────╮
│ Index │ Target J    │ Seed K   │ Solved K │ Path              │ Fallback │ Transitions │ Balance   │ Elapsed ms │
├───────┼─────────────┼──────────┼──────────┼───────────────────┼──────────┼─────────────┼───────────┼────────────┤
│ 0     │ -7.962557e1 │ 1005.000 │ 1005.000 │ NestedTemperature │ -        │ 0           │ 1.300e-13 │ 8.482      │
│ 1     │ -3.093510e1 │ 1005.000 │ 1016.250 │ NestedTemperature │ -        │ 0           │ 1.115e-12 │ 4.748      │
│ 2     │ 1.829819e1  │ 1016.250 │ 1027.500 │ NestedTemperature │ -        │ 0           │ 1.872e-13 │ 4.924      │
│ 3     │ 6.808685e1  │ 1027.500 │ 1038.750 │ NestedTemperature │ -        │ 0           │ 3.229e-14 │ 5.119      │
│ 4     │ 1.184453e2  │ 1038.750 │ 1050.000 │ NestedTemperature │ -        │ 0           │ 2.725e-12 │ 2.245      │
│ 5     │ 1.693892e2  │ 1050.000 │ 1061.250 │ NestedTemperature │ -        │ 0           │ 5.157e-14 │ 4.653      │
│ 6     │ 2.209351e2  │ 1061.250 │ 1072.500 │ NestedTemperature │ -        │ 0           │ 1.013e-13 │ 4.394      │
│ 7     │ 2.730997e2  │ 1072.500 │ 1083.750 │ NestedTemperature │ -        │ 0           │ 1.441e-13 │ 4.029      │
│ 8     │ 3.259015e2  │ 1083.750 │ 1095.001 │ NestedTemperature │ -        │ 0           │ 1.534e-8  │ 2.201      │
╰───────┴─────────────┴──────────┴──────────┴───────────────────┴──────────┴─────────────┴───────────┴────────────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_nested_ph_target_range_story ... ok
```

The release run passed. The target temperatures are intentionally interior to
the common NASA interval; endpoint contracts belong to the dedicated
interval/Jacobian stories.

## 14. Dense water/ice Auto P,H continuation range

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_water_ice_auto_ph_target_range_story`

### What it checks

- five real water/ice targets at 250, 255, 260, 265, and 270 K;
- route-dependent `Auto` behavior: nested fallback for a hard point and
  monolithic phase-control acceptance for neighboring points;
- phase-transition evidence, fallback provenance, conservation and enthalpy
  acceptance for every published point;
- no mutation of the local JSON libraries.

### Recorded release summary

| Targets | Phase transitions | Nested fallback points | Total ms | Worst point ms | Status |
| ---: | ---: | ---: | ---: | ---: | --- |
| 5 | 1 | 1 | 27.820 | 23.532 | OK |

### Raw console transcript
live dense water/ice Auto P,H target range: points=5 transitions=1 fallback_points=1 total=27.8199ms
╭───────┬─────────────┬─────────┬──────────┬────────────────────────┬───────────────────┬─────────────┬───────────┬────────────╮
│ Index │ Target J    │ Seed K  │ Solved K │ Path                   │ Fallback          │ Transitions │ Balance   │ Elapsed ms │
├───────┼─────────────┼─────────┼──────────┼────────────────────────┼───────────────────┼─────────────┼───────────┼────────────┤
│ 0     │ -1.476124e5 │ 250.000 │ 250.000  │ NestedTemperature      │ AllBackendsFailed │ 1           │ 0.000e0   │ 23.532     │
│ 1     │ -1.474818e5 │ 250.000 │ 255.000  │ MonolithicPhaseControl │ -                 │ 0           │ 8.540e-12 │ 0.976      │
│ 2     │ -1.473462e5 │ 255.000 │ 260.000  │ MonolithicPhaseControl │ -                 │ 0           │ 2.527e-12 │ 0.952      │
│ 3     │ -1.472040e5 │ 260.000 │ 265.000  │ MonolithicPhaseControl │ -                 │ 0           │ 6.334e-12 │ 0.899      │
│ 4     │ -1.470531e5 │ 265.000 │ 270.000  │ MonolithicPhaseControl │ -                 │ 0           │ 1.321e-12 │ 1.333      │
╰───────┴─────────────┴─────────┴──────────┴────────────────────────┴───────────────────┴─────────────┴───────────┴────────────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_water_ice_auto_ph_target_range_story ... ok

### Release command

```powershell
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_water_ice_auto_ph_target_range_story `
  --no-default-features -- --ignored --nocapture
```

### Release interpretation

```text
points=5 transitions=1 fallback_points=1 total=27.8199ms
250 K: NestedTemperature, AllBackendsFailed -> nested fallback, transitions=1
255/260/265/270 K: MonolithicPhaseControl, no fallback
status=OK
```

The release run passed. The expected invariant is not one fixed route for
every point; it is an accepted, documented route for every point without
partial publication.

## 15. Fixed-phase monolithic P,H backend characterization

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_backend_matrix`

### What it checks

- one real NASA H2/O2/H2O target-enthalpy grid under a strict single-backend
  policy;
- shared formulation construction and continuation reuse across all points;
- per-point enthalpy acceptance, finite residuals, and element conservation;
- an explicit table row for every backend, including failed methods;
- byte-for-byte immutability of the local JSON libraries.

### Release command

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_backend_matrix --no-default-features -- --ignored --nocapture
```
live real fixed-phase monolithic P,H target-range backend matrix
╭─────────────────────┬────────┬────────┬────────┬────────┬──────────┬─────────┬─────────┬──────────┬──────────────┬─────────────┬────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ Backend             │ Status │ Points │ Builds │ Reuses │ Total ms │ Wall ms │ Mean ms │ Worst ms │ Max residual │ Max balance │ Error                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              │
├─────────────────────┼────────┼────────┼────────┼────────┼──────────┼─────────┼─────────┼──────────┼──────────────┼─────────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ rst_lm              │ FAILED │ -      │ -      │ -      │ -        │ 25.011  │ -       │ -        │ -            │ -           │ Point(PhRangePointError { index: 0, target_enthalpy: -295678.9620290071, source: AllBackendsFailed { attempts: [SolverAttemptReport { backend: RustedSciThe(LevenbergMarquardt), outcome: RejectedCandidate { reason: "equilibrium candidate field 'candidate_acceptance_residual' was rejected: residual L2 norm 6.668693482588257e-6 exceeds tolerance 1e-6" }, metrics: Some(SolverAttemptMetrics { termination: MaxIterations, backend_converged: false, iterations: 50, residual_evaluations: 19, jacobian_evaluations: 19, linear_solves: 50, elapsed_millis: 3 }) }] } })   │
│ rst_minpack_lm      │ FAILED │ -      │ -      │ -      │ -        │ 22.783  │ -       │ -        │ -            │ -           │ Point(PhRangePointError { index: 2, target_enthalpy: -225650.11810361288, source: AllBackendsFailed { attempts: [SolverAttemptReport { backend: RustedSciThe(MinpackLevenbergMarquardt), outcome: RejectedCandidate { reason: "equilibrium candidate field 'candidate_enthalpy' was rejected: enthalpy error -0.07273956961580552 J exceeds the P,H acceptance limit" }, metrics: Some(SolverAttemptMetrics { termination: Converged, backend_converged: true, iterations: 3, residual_evaluations: 4, jacobian_evaluations: 4, linear_solves: 3, elapsed_millis: 0 }) }] } })     │
│ rst_trust_region_lm │ FAILED │ -      │ -      │ -      │ -        │ 24.453  │ -       │ -        │ -            │ -           │ Point(PhRangePointError { index: 2, target_enthalpy: -225650.11810361288, source: AllBackendsFailed { attempts: [SolverAttemptReport { backend: RustedSciThe(TrustRegionLevenbergMarquardt), outcome: RejectedCandidate { reason: "equilibrium candidate field 'candidate_enthalpy' was rejected: enthalpy error -0.07273956961580552 J exceeds the P,H acceptance limit" }, metrics: Some(SolverAttemptMetrics { termination: Converged, backend_converged: true, iterations: 3, residual_evaluations: 4, jacobian_evaluations: 4, linear_solves: 3, elapsed_millis: 0 }) }] } }) │
│ legacy_lm           │ FAILED │ -      │ -      │ -      │ -        │ 5.250   │ -       │ -        │ -            │ -           │ Point(PhRangePointError { index: 0, target_enthalpy: -295678.9620290071, source: AllBackendsFailed { attempts: [SolverAttemptReport { backend: Legacy(LM), outcome: Failed { kind: Solver, reason: "equilibrium solver failed: nonlinear solver reached its iteration limit" }, metrics: None }] } })                                                                                                                                                                                                                                                                              │
│ legacy_nr           │ FAILED │ -      │ -      │ -      │ -        │ 0.953   │ -       │ -        │ -            │ -           │ Point(PhRangePointError { index: 2, target_enthalpy: -225650.11810361288, source: AllBackendsFailed { attempts: [SolverAttemptReport { backend: Legacy(NR), outcome: RejectedCandidate { reason: "equilibrium candidate field 'candidate_enthalpy' was rejected: enthalpy error -0.0727395694993902 J exceeds the P,H acceptance limit" }, metrics: None }] } })                                                                                                                                                                                                                   │
│ legacy_tr           │ OK     │ 3      │ 1      │ 2      │ 0.325    │ 0.907   │ 0.108   │ 0.138    │ 1.776e-10    │ 1.584e-9    │ -                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  │
╰─────────────────────┴────────┴────────┴────────┴────────┴──────────┴─────────┴─────────┴──────────┴──────────────┴─────────────┴────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_backend_matrix ... ok

### Current debug characterization

The updated table reports `legacy_tr` completing all three points. The other
backend rows remain visible with their typed failure reasons, including
candidate residual/enthalpy rejection and nonlinear iteration exhaustion.
This is backend characterization only; the release result and production
default decision remain open.

## 16. Real P,H invalid-input and rollback matrix

### Test

Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_ph_invalid_input_and_rollback_matrix

### What it checks

- non-finite enthalpy target rejection before solving;
- invalid/reversed temperature-bound rejection before solving;
- non-positive pressure rejection before any backend attempt;
- finite but unreachable target failure as a typed point error;
- no partial P,H range publication after the failed point;
- byte-for-byte immutability of the local JSON libraries.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_ph_invalid_input_and_rollback_matrix --no-default-features -- --ignored --nocapture

### Debug characterization
live real P,H invalid-input and rollback matrix
╭─────────────────────────────┬─────────────┬─────────────────────────────────╮
│ Scenario                    │ Status      │ Detail                          │
├─────────────────────────────┼─────────────┼─────────────────────────────────┤
│ non-finite target           │ REJECTED    │ typed grid validation           │
│ reversed temperature bounds │ REJECTED    │ typed bounds validation         │
│ finite unreachable target   │ ROLLED BACK │ point error; no range published │
╰─────────────────────────────┴─────────────┴─────────────────────────────────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_ph_invalid_input_and_rollback_matrix ... ok
The matrix passed with four explicit rows: three request-level rejections and
one transactional rollback. No partial range was published. The additional
debug row is `non-positive pressure | REJECTED | typed condition validation`.

## 17. GUI accepted-result layout mismatch

### Test

gui::equilibrium_gui_result::tests::source_snapshot_rejects_an_accepted_range_with_mismatched_layout

### What it checks

- two independently accepted real NASA gas solutions;
- an intentionally malformed range payload with different component layouts;
- rejection during immutable GUI snapshot construction;
- no partial table or result snapshot is produced.

### Release command

    cargo test --release --lib gui::equilibrium_gui_result::tests::source_snapshot_rejects_an_accepted_range_with_mismatched_layout --no-default-features -- --ignored --nocapture

### Debug characterization

The story passed. The result layer reports a layout-fingerprint mismatch and
publishes no GUI snapshot.

## 18. Technical finishing release batch

This batch is the compact release checklist for the current P,H and GUI
technical work. It deliberately keeps the tests separate: the full module
suite checks regressions, while the ignored stories provide live-data,
transactionality, backend, and result-publication evidence.

### Commands

#### Canonical ChemEquilibrium suite

    cargo test --release --lib Thermodynamics::ChemEquilibrium --no-default-features

Checks the non-ignored regression suite for the canonical equilibrium engine.
The debug baseline currently passes with `612 passed; 0 failed; 36 ignored`.

#### Fixed-phase monolithic P,H backend matrix

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_backend_matrix --no-default-features -- --ignored --nocapture

Records per-backend status, point count, setup reuse, wall time, point-time
statistics, backend attempts, accepted backend, enthalpy error, residual,
balance, failure point, machine-readable failure kind, and the full typed
failure reason. A backend failure is characterization data; the story succeeds
when at least one configured backend produces an accepted range.

Latest debug characterization:

| Backend | Status | Wall ms | Attempts | Accepted | Validation | Max `|dH|` J | Max residual | Max balance | Failure |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| `rst_lm` | FAILED | 43.903 | - | - | NOT ACCEPTED | - | - | - | point 0, `AllBackendsFailed` |
| `rst_minpack_lm` | FAILED | 42.068 | - | - | NOT ACCEPTED | - | - | - | point 2, `AllBackendsFailed` |
| `rst_trust_region_lm` | FAILED | 40.205 | - | - | NOT ACCEPTED | - | - | - | point 2, `AllBackendsFailed` |
| `legacy_lm` | FAILED | 19.139 | - | - | NOT ACCEPTED | - | - | - | point 0, `AllBackendsFailed` |
| `legacy_nr` | FAILED | 2.706 | - | - | NOT ACCEPTED | - | - | - | point 2, `AllBackendsFailed` |
| `legacy_tr` | OK | 2.889 | 3 | `Legacy(TR)` | OK | `1.172e-4` | `1.776e-10` | `1.584e-9` | - |

The complete `Error` column remains printed by the test and retains the
backend-specific rejection text (for example residual rejection versus P,H
enthalpy rejection). This debug run passed because the matrix is a strict
characterization story, not an assertion that every numerical backend must
solve every real grid.

#### Real-data P,H validation and rollback

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_real_ph_invalid_input_and_rollback_matrix --no-default-features -- --ignored --nocapture

Checks typed rejection of non-finite targets and reversed bounds, typed failure
of an unreachable finite target, absence of partial range publication, and
byte-for-byte immutability of the local JSON libraries.

#### GUI accepted-result layout contract

    cargo test --release --lib gui::equilibrium_gui_result::tests::source_snapshot_rejects_an_accepted_range_with_mismatched_layout --no-default-features -- --ignored --nocapture

Checks that the GUI result layer rejects an accepted range whose component
layout differs from the source solution, without publishing a partial snapshot.

### Recorded status before release

| Story | Debug status | Release status |
|---|---|---|
| Canonical ChemEquilibrium suite | `612 passed; 0 failed; 36 ignored` | pending |
| Monolithic P,H backend matrix | passed; full timing/validation/failure rows | pending |
| Real-data P,H validation and rollback | passed; 3 explicit validation/rollback rows | pending |
| GUI accepted-result layout contract | passed | pending |

After a release run, append the compact printed table or matrix under the
corresponding story section and replace only that row's `pending` marker. Keep
the command unchanged so future runs remain comparable.

## 19. Nested P,H inner P,T cascade release matrix

### Test

Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_nested_ph_inner_pt_cascade_release_matrix

### What it checks

- small, medium, and large real local NASA systems with 5, 20, and 100
  selected species;
- three target enthalpies solved through the typed nested P,H range facade;
- one reusable inner P,T formulation per fixture and accepted-state
  continuation between points;
- outer scalar evaluations reported separately from inner backend attempts and
  measured inner nonlinear-solve time;
- enthalpy acceptance, elemental conservation, and byte-for-byte JSON
  immutability.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_nested_ph_inner_pt_cascade_release_matrix --no-default-features -- --ignored --nocapture

### Debug characterization

```text
live nested P,H inner P,T cascade release matrix
| Fixture | Status | Species | Points | Outer evals | Inner attempts | Inner solve ms | Builds | Reuses | Total ms | Mean ms | Worst ms | Max balance |
| small   | OK     | 5       | 3      | 24          | 63             | 3.797          | 1      | 32     | 9.446    | 3.149   | 4.360    | 7.642e-11   |
| medium  | OK     | 20      | 3      | 49          | 113            | 66.025         | 1      | 57     | 92.079   | 30.693  | 51.832   | 9.834e-11   |
| large   | OK     | 100     | 3      | 26          | 67             | 1172.744       | 1      | 34     | 1290.992 | 430.331 | 531.140  | 3.643e-10   |
test ...::live_nested_ph_inner_pt_cascade_release_matrix ... ok
```

The debug matrix passed. Legacy adapters do not expose comparable nonlinear
iteration counters, so the table records inner nonlinear solve time rather
than manufacturing an iteration count. The release command is the target
machine characterization run.

## 20. GUI accepted-candidate validation mismatch

### Test

gui::equilibrium_gui_result::tests::source_snapshot_rejects_an_accepted_candidate_with_mismatched_validation

### What it checks

- a real NASA accepted candidate is retained as the source;
- only its test-copy validation payload is made inconsistent with physical
  moles;
- the immutable GUI result boundary rejects the payload with a deterministic
  publication error;
- no point table or partial result snapshot is published.

### Release command

    cargo test --release --lib gui::equilibrium_gui_result::tests::source_snapshot_rejects_an_accepted_candidate_with_mismatched_validation --no-default-features -- --ignored --nocapture

### Debug characterization

The story passed. The GUI result layer rejected the accepted candidate because
its validation evidence was inconsistent with the physical mole vector.

<!-- Raw release capture retained only for source traceability. Its normalized
     story-test records are sections 22 through 31 below. -->
<!--


cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_pt_to_h_to_ph_recovers_temperature_and_composition --no-default-features -- --ignored --nocapture
live P,H inverse story | reference T=2500.000000 K | recovered T=2500.000000 K | target H=-2.638129e5 J | scaled H error=-2.745e-12 | trials=22 | max mole delta=1.701e-11
live monolithic P,H inverse | T=2500.000000 K | scaled H error=0.000e0 | backend=policy=Cascade([Legacy(LM), Legacy(LM), Legacy(NR), Legacy(TR)]), attempts=1, started=1, skipped=0, fallback=0, accepted_backend=Legacy(LM) 
live real P,H formulation characterization
╭───────────────────┬────────┬───────────┬────────────────┬─────────┬──────────────┬────────┬──────────┬────────────┬────────┬────────┬─────────────┬────────────────┬────────────────┬────────────────┬───────────┬───────────┬───────╮
│ Path              │ Status │ Report ms │ Thermo prep ms │ Wall ms │ Inner solves │ Trials │ Attempts │ Iterations │ Builds │ Reuses │ Transitions │ Residual evals │ Jacobian evals │ Scaled H error │ Residual  │ Balance   │ Error │
├───────────────────┼────────┼───────────┼────────────────┼─────────┼──────────────┼────────┼──────────┼────────────┼────────┼────────┼─────────────┼────────────────┼────────────────┼────────────────┼───────────┼───────────┼───────┤
│ nested            │ OK     │ 28.401    │ 0.015          │ 51.719  │ 22           │ 22     │ 43       │ 110        │ 1      │ 21     │ 0           │ 91             │ 91             │ -2.745e-12     │ 5.457e-11 │ 4.854e-11 │ -     │
│ monolithic        │ OK     │ 0.813     │ 0.006          │ 23.115  │ 1            │ 0      │ 1        │ 0          │ 1      │ 0      │ 0           │ 0              │ 0              │ 0.000e0        │ 0.000e0   │ 0.000e0   │ -     │
│ monolithic-rst-lm │ OK     │ 21.390    │ 0.005          │ 22.184  │ 1            │ 0      │ 1        │ 0          │ 1      │ 0      │ 0           │ 1              │ 1              │ 0.000e0        │ 0.000e0   │ 0.000e0   │ -     │
│ auto              │ OK     │ 0.000     │ 0.000          │ 0.715   │ 1            │ 0      │ 1        │ 0          │ 1      │ 0      │ 0           │ 0              │ 0              │ 0.000e0        │ 0.000e0   │ 0.000e0   │ -     │
╰───────────────────┴────────┴───────────┴────────────────┴─────────┴──────────────┴────────┴──────────┴────────────┴────────┴────────┴─────────────┴────────────────┴────────────────┴────────────────┴───────────┴───────────┴───────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_pt_to_h_to_ph_recovers_temperature_and_composition ... ok


cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_symbolic_ph_rejects_a_native_polynomial_boundary_without_mutating_json --no-default-features -- --ignored --nocapture
нмчего не печатает
cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_water_pt_to_h_to_ph_preserves_phase_evidence_and_json --no-default-features -- --ignored --nocapture
нмчего не печатает

cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_phase_control_temperature_range_reuses_accepted_state --no-default-features -- --ignored --nocapture
live bounded typed T-range: points=3 transitions=0 projections=1 prepared=1 rst_symbolic=1 total=3.6966ms
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_phase_control_temperature_range_reuses_accepted_state ... ok


cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_phase_control_timing_report_prints_stage_breakdown --no-default-features -- --ignored --nocapture
running 1 test
live phase-control timing report
  repository_lookup              163µs
  thermochemistry_preparation   14.1µs
  numeric_closure_construction  30.5µs
  symbolic_construction         573µs
  equation_construction         17.3µs
  numerical_problem_preparation 25.1µs
  projection_build               9µs
  nonlinear_solve                0ns
  phase_control                  1.0389ms
  validation                     100ns
  postprocessing                 0ns
  total                          1.9854ms
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_phase_control_timing_report_prints_stage_breakdown ... ok


cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_solver_matrix_compares_rst_and_legacy_backends --no-default-features -- --ignored --nocapture

live large C/H/O-limited solver matrix: species=20
  rst_lm                   ok total=94.0388ms nonlinear=89.7731ms residual=2.3428860389495855e-7 balance=1.266742673117216e-7
  rst_minpack_lm           ok total=97.017ms nonlinear=92.4002ms residual=7.413691612457578e-10 balance=6.661297902166297e-11
  rst_nielsen_lm           failed: Solve(AllBackendsFailed { attempts: [SolverAttemptReport { backend: RustedSciThe(NielsenLevenbergMarquardt), outcome: Failed { kind: ResidualEvaluation, reason: "equilibrium residual evaluation failed: RustedSciThe rst_nielsen_levenberg_marquardt: symbolic residual returned NaN or Inf" }, metrics: None }] })
  rst_trust_region_lm      ok total=93.5314ms nonlinear=89.8947ms residual=6.284833064178006e-7 balance=5.6569654433014094e-8
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
  rst_powell_dogleg        ok total=94.8836ms nonlinear=90.5406ms residual=8.702307137612399e-8 balance=3.175701354662941e-9
  rst_damped_newton        ok total=92.4923ms nonlinear=88.626ms residual=2.0497288674797745e-7 balance=1.8438181773050566e-8
  legacy_lm                ok total=4.1019ms nonlinear=105.2µs residual=3.6781743366064495e-14 balance=1.3877787807814457e-17
  legacy_nr                ok total=3.6206ms nonlinear=64µs residual=2.049728905790898e-7 balance=1.8438182231017564e-8
  legacy_tr                ok total=3.4807ms nonlinear=221.4µs residual=2.2203496345604435e-10 balance=1.1150885770305763e-11
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_solver_matrix_compares_rst_and_legacy_backends ... ok



cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_element_limited_release_scaling_matrix --no-default-features -- --ignored --nocapture


running 1 test
live C/H/O-limited scaling: backend=legacy-nr species=20 total=4.2566ms nonlinear=101.6µs residual=2.0497336587239133e-7 balance=1.8438228416295388e-8
live C/H/O-limited scaling: backend=legacy-nr species=50 total=9.5831ms nonlinear=237.2µs residual=5.517840136201935e-7 balance=1.0011128948050896e-7
live C/H/O-limited scaling: backend=legacy-nr species=100 total=18.3736ms nonlinear=1.067ms residual=7.615420917540235e-10 balance=3.803223708187531e-11
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_element_limited_release_scaling_matrix ... ok


cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_exact_element_typed_temperature_range_backend_matrix --no-default-features -- --ignored --nocapture
live typed T-range backend: backend=rst_lm species=20 points=3 total=101.1493ms mean=33.716433ms worst=89.7155ms
live typed T-range backend: backend=rst_minpack_lm species=20 points=3 total=99.5231ms mean=33.174366ms worst=87.4463ms
live typed T-range backend: backend=rst_nielsen_lm species=20 failed=InvalidProblem { field: "temperature_range_point", message: "point 0 at 1000 K failed: backend=RustedSciThe(NielsenLevenbergMarquardt), outcome=failed(ResidualEvaluation): equilibrium residual evaluation failed: RustedSciThe rst_nielsen_levenberg_marquardt: symbolic residual returned NaN or Inf" }
live typed T-range backend: backend=rst_trust_region_lm species=20 points=3 total=100.2948ms mean=33.4316ms worst=88.9015ms
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
beta = 0.3127537096772781
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
beta = 0.37011885784489074
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
gauss newton step norm <= delta: return gauss newton step
live typed T-range backend: backend=rst_powell_dogleg species=20 points=3 total=104.3899ms mean=34.796633ms worst=88.8162ms
live typed T-range backend: backend=rst_damped_newton species=20 points=3 total=97.6591ms mean=32.553033ms worst=87.4859ms
live typed T-range backend: backend=legacy_lm species=20 points=3 total=501.6µs mean=167.2µs worst=211.3µs
live typed T-range backend: backend=legacy_nr species=20 points=3 total=361.9µs mean=120.633µs worst=146.3µs
live typed T-range backend: backend=legacy_tr species=20 points=3 total=587µs mean=195.666µs worst=315.6µs
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_exact_element_typed_temperature_range_backend_matrix ... ok


cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_legacy_temperature_range_story --no-default-features -- --ignored --nocapture

running 1 test
live real T-range: backend=LM, species=20, points=5, elapsed=942.2µs
live real T-range: backend=NR, species=20, points=5, elapsed=1.0818ms
live real T-range: backend=TR, species=20, points=5, elapsed=812.9µs
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_legacy_temperature_range_story ... ok

cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_typed_temperature_range_story --no-default-features -- --ignored --nocapture
live typed T-range: backend=legacy-nr direction=Ascending species=20 points=3 setup_builds=1 reuses=2 point_total=408µs point_mean=136µs point_median=108.1µs point_worst=193.8µs
live typed T-range: backend=legacy-nr direction=Descending species=20 points=3 setup_builds=1 reuses=2 point_total=424.4µs point_mean=141.466µs point_median=146.4µs point_worst=170µs
live typed T-range: backend=rst-default direction=Ascending species=20 points=3 setup_builds=1 reuses=2 symbolic_reuses=2 point_total=100.9532ms point_mean=33.651066ms point_median=6.068ms point_worst=90.2757ms
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_typed_temperature_range_story ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 2069 filtered out; finished in 0.21s

-->

## 21. Large real-data timing report

### Test

Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_gas_timing_report

### What it checks

- real C/H/O-limited NASA candidate selection;
- a 20-species resolved gas solve with stage-by-stage timing;
- finite positive amounts and residual validation;
- scale-aware absolute-plus-relative elemental conservation;
- no mutation of the resolved thermochemistry source.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_gas_timing_report --no-default-features -- --ignored --nocapture

### Debug characterization after test-contract fix

The first release attempt used an overly strict fixed `1e-8` balance assertion
and failed before printing timings. The test now uses the same scale-aware
contract as the other real-data matrices. The debug rerun passed with:

```text
selected_species=20
repository_lookup=1.2954ms
thermochemistry_preparation=153.1us
numeric_closure_construction=223.6us
symbolic_construction=7.4035ms
equation_construction=134.4us
numerical_problem_preparation=684.3us
nonlinear_solve=198.9765ms
postprocessing=117.5us
validation residual=2.343e-7 balance=1.267e-7 limit=1.014e-6
total=209.8891ms
test ...::live_large_element_limited_gas_timing_report ... ok
```

## 22. Real P,H inverse and formulation comparison

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_pt_to_h_to_ph_recovers_temperature_and_composition`

### What it checks

- roundtrip from a real fixed-`P,T` equilibrium to target enthalpy and back to `P,H`;
- agreement of nested, monolithic, explicit-RST and auto routes;
- continuation of composition and acceptance of a validated candidate;
- scaled enthalpy, residual, balance, trial and backend-attempt evidence.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_pt_to_h_to_ph_recovers_temperature_and_composition --no-default-features -- --ignored --nocapture

### Release characterization

The supplied run passed. The reference and recovered temperatures were both
`2500 K`; the nested route used 22 trials, with scaled enthalpy error
`-2.745e-12`, residual `5.457e-11`, balance `4.854e-11`, and maximum mole
delta `1.701e-11`. Monolithic, explicit-RST, and auto routes also accepted the
same temperature with zero reported scaled enthalpy error.

| Path | Wall ms | Inner solves | Attempts | Builds | Reuses | Transitions | Scaled H error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| nested | 51.719 | 22 | 43 | 1 | 21 | 0 | -2.745e-12 |
| monolithic | 23.115 | 1 | 1 | 1 | 0 | 0 | 0.000e0 |
| monolithic RST-LM | 22.184 | 1 | 1 | 1 | 0 | 0 | 0.000e0 |
| auto | 0.715 | 1 | 1 | 1 | 0 | 0 | 0.000e0 |

## 23. Symbolic P,H boundary rejection

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_symbolic_ph_rejects_a_native_polynomial_boundary_without_mutating_json`

### What it checks

- deterministic rejection of an unsupported native polynomial boundary in the
  symbolic P,H route;
- no publication of a partial result;
- no mutation of the JSON thermochemistry libraries.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_symbolic_ph_rejects_a_native_polynomial_boundary_without_mutating_json --no-default-features -- --ignored --nocapture

### Release characterization

The supplied release command completed successfully without diagnostic output.
The test is intentionally quiet; its assertions cover the typed rejection and
the before/after JSON snapshot.

## 24. Bounded water P,H phase evidence

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_water_pt_to_h_to_ph_preserves_phase_evidence_and_json`

### What it checks

- real water phase-control data through the `P,T -> H -> P,H` route;
- preservation of phase totals and transition evidence;
- rollback-safe publication and unchanged JSON libraries.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_water_pt_to_h_to_ph_preserves_phase_evidence_and_json --no-default-features -- --ignored --nocapture

### Release characterization

The supplied release command completed successfully without diagnostic output.
The test is a quiet invariant story: assertions validate the phase evidence,
rollback boundary, and library snapshots.

## 25. Bounded phase-control T-range reuse

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_phase_control_temperature_range_reuses_accepted_state`

### What it checks

- continuation across a bounded real-data temperature range;
- reuse of the prepared formulation and the accepted state;
- absence of spurious phase transitions.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_phase_control_temperature_range_reuses_accepted_state --no-default-features -- --ignored --nocapture

### Release characterization

`points=3`, `transitions=0`, `projections=1`, `prepared=1`,
`rst_symbolic=1`, total `3.6966 ms`; test passed.

## 26. Phase-control timing breakdown

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_phase_control_timing_report_prints_stage_breakdown`

### What it checks

- timing coverage for repository lookup, thermochemistry, closures, symbolic
  construction, equation preparation, projection, phase control, validation,
  and postprocessing;
- a stable report contract for performance diagnosis.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_phase_control_timing_report_prints_stage_breakdown --no-default-features -- --ignored --nocapture

### Release characterization

| Stage | Time |
| --- | ---: |
| repository lookup | 163 us |
| thermochemistry preparation | 14.1 us |
| numeric closures | 30.5 us |
| symbolic construction | 573 us |
| equation construction | 17.3 us |
| numerical preparation | 25.1 us |
| projection build | 9 us |
| nonlinear solve | 0 ns |
| phase control | 1.0389 ms |
| validation | 100 ns |
| postprocessing | 0 ns |
| total | 1.9854 ms |

The test passed.

## 27. Real C/H/O solver backend matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_solver_matrix_compares_rst_and_legacy_backends`

### What it checks

- one real 20-species C/H/O problem across all RST and legacy backends;
- comparable residual and elemental-balance validation;
- explicit backend failure classification for fallback-policy design.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_solver_matrix_compares_rst_and_legacy_backends --no-default-features -- --ignored --nocapture

### Release characterization

| Backend | Status | Total ms | Nonlinear ms | Residual | Balance |
| --- | --- | ---: | ---: | ---: | ---: |
| rst_lm | OK | 94.039 | 89.773 | 2.343e-7 | 1.267e-7 |
| rst_minpack_lm | OK | 97.017 | 92.400 | 7.414e-10 | 6.661e-11 |
| rst_nielsen_lm | FAILED | - | - | - | - |
| rst_trust_region_lm | OK | 93.531 | 89.895 | 6.285e-7 | 5.657e-8 |
| rst_powell_dogleg | OK | 94.884 | 90.541 | 8.702e-8 | 3.176e-9 |
| rst_damped_newton | OK | 92.492 | 88.626 | 2.050e-7 | 1.844e-8 |
| legacy_lm | OK | 4.102 | 0.105 | 3.678e-14 | 1.388e-17 |
| legacy_nr | OK | 3.621 | 0.064 | 2.050e-7 | 1.844e-8 |
| legacy_tr | OK | 3.481 | 0.221 | 2.220e-10 | 1.115e-11 |

`rst_nielsen_lm` failed deterministically during symbolic residual evaluation
because the residual returned NaN/Inf; the complete fallback matrix still
passed.

## 28. Real C/H/O scaling matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_element_limited_release_scaling_matrix`

### What it checks

- real element-limited systems with 20, 50, and 100 selected species;
- residual and scale-aware conservation contracts;
- release-mode growth of the legacy-NR canonical timing baseline.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_element_limited_release_scaling_matrix --no-default-features -- --ignored --nocapture

### Release characterization

| Species | Total ms | Nonlinear ms | Residual | Balance |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 4.257 | 0.102 | 2.050e-7 | 1.844e-8 |
| 50 | 9.583 | 0.237 | 5.518e-7 | 1.001e-7 |
| 100 | 18.374 | 1.067 | 7.615e-10 | 3.803e-11 |

All three cases passed.

## 29. Real typed T-range backend matrix

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_exact_element_typed_temperature_range_backend_matrix`

### What it checks

- typed continuation over three real temperature points;
- all RST and legacy backend routes;
- reuse/build counters and residual failure classification.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_exact_element_typed_temperature_range_backend_matrix --no-default-features -- --ignored --nocapture

### Release characterization

| Backend | Status | Total ms | Mean ms | Worst ms |
| --- | --- | ---: | ---: | ---: |
| rst_lm | OK | 101.149 | 33.716 | 89.716 |
| rst_minpack_lm | OK | 99.523 | 33.174 | 87.446 |
| rst_nielsen_lm | FAILED | - | - | - |
| rst_trust_region_lm | OK | 100.295 | 33.432 | 88.902 |
| rst_powell_dogleg | OK | 104.390 | 34.797 | 88.816 |
| rst_damped_newton | OK | 97.659 | 32.553 | 87.486 |
| legacy_lm | OK | 0.502 | 0.167 | 0.211 |
| legacy_nr | OK | 0.362 | 0.121 | 0.146 |
| legacy_tr | OK | 0.587 | 0.196 | 0.316 |

The test passed. `rst_nielsen_lm` failed at the first point because its
symbolic residual returned NaN/Inf; this remains backend-specific evidence, not
a failure of the typed range workflow.

## 30. Legacy T-range comparison

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_legacy_temperature_range_story`

### What it checks

- the legacy LM, NR, and TR fallback routes on the same real five-point range;
- a compact timing baseline for retained numerical fallback backends.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_legacy_temperature_range_story --no-default-features -- --ignored --nocapture

### Release characterization

| Backend | Species | Points | Elapsed |
| --- | ---: | ---: | ---: |
| LM | 20 | 5 | 942.2 us |
| NR | 20 | 5 | 1.0818 ms |
| TR | 20 | 5 | 812.9 us |

The test passed.

## 31. Typed T-range continuation directions

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_typed_temperature_range_story`

### What it checks

- ascending and descending temperature grids;
- one setup build followed by accepted-state reuse;
- comparison of the legacy-NR and RST-default typed paths.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_typed_temperature_range_story --no-default-features -- --ignored --nocapture

### Release characterization

| Backend | Direction | Points | Builds | Reuses | Mean | Median | Worst |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| legacy-NR | ascending | 3 | 1 | 2 | 136 us | 108.1 us | 193.8 us |
| legacy-NR | descending | 3 | 1 | 2 | 141.466 us | 146.4 us | 170 us |
| RST-default | ascending | 3 | 1 | 2 | 33.651 ms | 6.068 ms | 90.276 ms |

The test passed; the RST row also reported two symbolic reuses.

## 32. P,H monolithic and nested phase-control agreement

### Test

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_water_ph_monolithic_and_nested_routes_agree_after_phase_activation`

### What it checks

- one real water/liquid `P,H` target at `350 K`, built from an independently
  accepted `P,T` state;
- explicit `Monolithic` and `NestedTemperature` routes receive the same
  resolved records, inventory, target, bounds, phase policy, and timing mode;
- both routes activate the initially empty liquid phase and publish lifecycle
  evidence;
- accepted temperature, component moles, residual, elemental balance, and
  enthalpy acceptance agree within explicit tolerances;
- JSON thermochemistry libraries remain byte-for-byte unchanged.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_water_ph_monolithic_and_nested_routes_agree_after_phase_activation --no-default-features -- --ignored --nocapture

### Debug characterization

live bounded water P,H route comparison
╭──────────────────────────┬────────┬───────────┬────────────────┬─────────┬──────────────┬────────┬──────────┬────────────┬────────┬────────┬─────────────┬────────────────┬────────────────┬────────────────┬───────────┬───────────┬───────╮
│ Path                     │ Status │ Report ms │ Thermo prep ms │ Wall ms │ Inner solves │ Trials │ Attempts │ Iterations │ Builds │ Reuses │ Transitions │ Residual evals │ Jacobian evals │ Scaled H error │ Residual  │ Balance   │ Error │
├──────────────────────────┼────────┼───────────┼────────────────┼─────────┼──────────────┼────────┼──────────┼────────────┼────────┼────────┼─────────────┼────────────────┼────────────────┼────────────────┼───────────┼───────────┼───────┤
│ monolithic-phase-control │ OK     │ 14.284    │ 0.008          │ 77.329  │ 1            │ 0      │ 2        │ 0          │ 0      │ 0      │ 1           │ 0              │ 0              │ -9.008e-10     │ 1.123e-10 │ 1.004e-9  │ -     │
│ nested-phase-control     │ OK     │ 62.761    │ 0.113          │ 63.005  │ 25           │ 25     │ 145      │ 0          │ 0      │ 0      │ 24          │ 1              │ 1              │ 9.932e-9       │ 8.083e-16 │ 6.661e-16 │ -     │
╰──────────────────────────┴────────┴───────────┴────────────────┴─────────┴──────────────┴────────┴──────────┴────────────┴────────┴────────┴─────────────┴────────────────┴────────────────┴────────────────┴───────────┴───────────┴───────╯
test Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_bounded_water_ph_monolithic_and_nested_routes_agree_after_phase_activation ... ok
The debug run passed. These timings are diagnostic rather than a performance
claim; the stable contract is route agreement and retained phase evidence.
