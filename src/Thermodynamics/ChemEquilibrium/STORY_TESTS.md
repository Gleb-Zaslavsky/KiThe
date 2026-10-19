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

## 3a. Production default at upper ElementInventory scale

### Story name

Real TP-1907 P,T production default must retain a scale-robust backend at
inventory scale `1e8`.

### Tests

`Thermodynamics::ChemEquilibrium::equilibrium_multiphase_story_tests::production_default_retains_a_scale_robust_backend_on_real_tp1907_story`

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_is_backend_invariant`

### What it checks

- the real offline NASA TP-1907 CHON + graphite system at inventory scale
  `1e8`;
- the production default cascade publishes an accepted result from the
  scale-robust backend set;
- an explicit `Single` policy cannot hide backend-specific outcomes behind a
  fallback;
- accepted states preserve topology, finite physical moles, element totals,
  residuals, relative conservation, and the molecular oracle;
- default solver selection remains protected against promoting a currently
  fragile isolated RST method to the only nonlinear backend.

### Hypothesis / question

Does the production default cascade still publish a valid physical result for
a real `ElementInventory` problem at `1e8`, even when an individual nonlinear
backend is not robust at that scale?

### Acceptance criteria

- the accepted backend belongs to the qualified robust set;
- the published physical state remains finite and conserves the elemental
  inventory within the scale-aware contract;
- an isolated backend failure does not invalidate the complete production
  cascade.

### Release commands

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_multiphase_story_tests::production_default_retains_a_scale_robust_backend_on_real_tp1907_story `
  --no-default-features -- --include-ignored --exact --nocapture

cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_is_backend_invariant `
  --no-default-features -- --include-ignored --exact --nocapture
```

### Recorded result

```text
TP-1907 backend/default story: test ...::production_default_retains_a_scale_robust_backend_on_real_tp1907_story ... ok

TP-1907 backend matrix at scale=1e8:
  rst-lm: FAILED, normalized basin discovery, candidate residual 5.12501771773085e1
  rst-minpack-lm: FAILED, normalized basin discovery, candidate residual 1.6277711838251514e2
  rst-nielsen-lm: FAILED, normalized basin discovery, candidate residual 1.6277711838251489e2
  rst-trust-region-lm: residual=2.712e-14, balance=1.192e-7
  rst-damped-newton: FAILED, normalized basin discovery, candidate residual 8.859999596823666e1
  legacy-lm: residual=1.557e-14, balance=1.192e-7
  legacy-nr: residual=1.693e-14, balance=3.532e-6
  legacy-tr: residual=1.474e-14, balance=1.192e-7
  accepted: rst-trust-region-lm, legacy-lm, legacy-nr, legacy-tr
  failed in isolated Single mode: rst-lm, rst-minpack-lm, rst-nielsen-lm, rst-damped-newton
test ...::production_default_retains_a_scale_robust_backend_on_real_tp1907_story ... ok
test ...::tp1907_upper_inventory_is_backend_invariant ... ok
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 2559 filtered out
```

The four isolated RST failures are strict normalized-recovery outcomes, not
evidence that the thermodynamic equations or the complete production cascade
are invalid. RST Trust-Region and all retained legacy methods return a
candidate, while Legacy NR reports a physical absolute balance of about
`3.53e-6` at this scale and is accuracy-failing in the strict matrix. The
production default must therefore remain a cascade, not an isolated RST
LM/Nielsen/Damped Newton selection.

### Current facade acceptance matrix: local and RST backends

#### Story name

TP-1907 upper `ElementInventory` facade must retain at least one physically
accurate backend in each solver family.

#### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_frozen_facade_upper_inventory_matches_the_same_oracle`

#### Command

```powershell
cargo test --lib `
  Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_frozen_facade_upper_inventory_matches_the_same_oracle `
  --no-default-features -- --exact --nocapture
```

#### Question / hypothesis

Does the real TP-1907 upper-scale `ElementInventory` problem have at least
one acceptable RST method and at least one acceptable local Legacy method when
each method is exercised through the normal typed facade with an explicit
`Single` policy? A solver is not accepted merely because it returns or reports
convergence: the published physical state must pass the same oracle and
conservation checks.

#### Classification

- `PASS`: scaled and raw residuals, physical elemental balance, component
  amounts, mole fractions, topology, finiteness, and non-negativity all pass;
- `ACCURACY_FAIL`: a candidate was returned or rejected by the acceptance
  gate, but it did not meet the physical/oracle tolerances;
- `HARD_FAIL`: the backend produced no candidate because of a typed backend,
  singularity, evaluation, or iteration failure without an accepted-candidate
  rejection.

The test prints every backend row and detailed attempt trace. It fails only
when no RST backend passes or when no Legacy backend passes. Thus a weak
isolated method remains visible without hiding the fact that the family still
has a usable method.

#### Recorded result

```text
TP-1907 upper inventory facade backend acceptance matrix
| backend | family | status | scaled residual | raw residual | max physical balance | max component rel | max fraction error | topology |
|---|---|---|---:|---:|---:|---:|---:|---|
| rst-lm | RST | ACCURACY_FAIL | - | - | - | - | - | false |
| rst-minpack-lm | RST | ACCURACY_FAIL | - | - | - | - | - | false |
| rst-nielsen-lm | RST | ACCURACY_FAIL | - | - | - | - | - | false |
| rst-trust-region-lm | RST | PASS | 2.712e-14 | 2.712e-14 | 1.192e-7 | 6.742e-14 | 1.665e-16 | true |
| rst-damped-newton | RST | ACCURACY_FAIL | - | - | - | - | - | false |
| legacy-lm | Legacy | PASS | 1.557e-14 | 1.557e-14 | 1.192e-7 | 3.091e-13 | 3.053e-16 | true |
| legacy-nr | Legacy | ACCURACY_FAIL | 1.693e-14 | 1.693e-14 | 3.532e-6 | 6.325e-12 | 3.469e-16 | true |
| legacy-tr | Legacy | PASS | 1.474e-14 | 1.474e-14 | 1.192e-7 | 1.616e-13 | 1.388e-16 | true |
```

Summary: RST `1/5` and Legacy `2/3` methods passed the complete physical
gate. No method was a `HARD_FAIL` in this run. RST LM, Minpack LM, Nielsen LM,
and Damped Newton returned candidates rejected by the strict acceptance
residual; Legacy NR returned a numerically small residual but missed the
physical balance limit (`3.532e-6 > 1e-6`). The family-level assertions passed
because Trust-Region LM, Legacy LM, and Legacy TR matched the independent
molecular oracle.

This matrix deliberately differs from the production-default test above:
the production cascade may accept RST LM under its ordinary solver tolerance
before the stronger physical balance gate is applied. That behavior remains a
separate correctness issue and is now visible rather than being confused with
the family-level availability result.

### Conclusion

The hypothesis is supported. The production default remains valid at `1e8`
because it retains a scale-robust backend in its cascade. The explicit matrix
also records which isolated methods are currently unsuitable as the sole
default backend.

## 3a.1. Paste-ready RST witness for the upper-scale failure

### Story name

TP-1907 upper `ElementInventory` scale must have a standalone symbolic RST
reproducer.

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_emits_paste_ready_rst_symbolic_problem`

### What it checks

- prepares the real TP-1907 fixed-`P,T` graph at inventory scale `1e8`;
- exports the exact symbolic residual expressions, variable order, numeric
  equation parameters (`T`, `g0_i`), initial log-moles, and production bounds;
- runs every recommended RST method plus explicit Powell Dogleg in `Single`
  mode, without phase control, normalization, or facade fallback logic;
- prints the full physical contract, rank evidence, raw/scaled residuals,
  molecular oracle, termination, callback counters, and typed backend errors
  for direct transfer to the RustedSciThe test suite.

The `1e8` qualifier is intentional. This witness does not claim that every
inventory scale is inaccurate: the neighboring facade scale tests cover
`1e-3` through `1e6`, while this reproducer isolates the upper stress point
that currently exposes the conditioning failure.

### Hypothesis / question

Can the RST convergence problem observed through the TP-1907 facade be reduced
to the prepared symbolic residual/Jacobian problem and reproduced without any
KiThe thermochemistry or orchestration code?

### Command

```powershell
cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_emits_paste_ready_rst_symbolic_problem `
  --no-default-features -- --include-ignored --exact --nocapture *> tp1907_rst_reproducer.txt
```

The generated `tp1907_rst_reproducer.txt` contains a Rust `vec![...]` payload
that can be copied into an RST-side test. The payload is intentionally emitted
by an ignored diagnostic test and is not part of the ordinary correctness
gate.

### Full physical contract emitted by the witness

The output is now a complete transfer contract, not only a symbolic equation
dump:

- physical point: `T = 700 K`, `P = 101325 Pa`, inventory multiplier `1e8`;
- element order: `Ar, C, H, N, O`;
- species order: `gas::Ar, gas::CH4, gas::CO, gas::CO2, gas::H, gas::HO2,
  gas::H2, gas::H2O, gas::N, gas::NH3, gas::NO, gas::NO2, gas::N2,
  gas::N2O, gas::O, gas::OH, gas::O2, graphite::C(gr)`;
- the complete `A[element, species]` matrix, `b`, `species_phase`, phase
  membership, bounds, seed, and all numeric Gibbs parameters are printed in
  `tp1907_rst_reproducer.txt`;
- the physical balance is exactly `balance[e] = sum_i A[e,i] * n[i] - b[e]`;
- phase order is `gas` with species `0..17` and `graphite` with species `17`.

The printed `residual_row_scale_w` name is historical. It is the positive
denominator vector `D`, not the diagonal of `W`:

```text
F_raw(y) = [reaction-affinity rows; A*n - b]
W = diag(1 / D)
F_scaled(y) = W * F_raw(y)
F_scaled[i] = F_raw[i] / D[i]
J_scaled[i,:] = J_raw[i,:] / D[i]
```

For reaction rows `D[i] = ||nu_i||`; for elemental rows
`D[r + e] = abs(b[e])`, or `1` when `b[e] == 0`. Thus the elemental scaled
row is `(A_e*n - b_e) / abs(b_e)` and is never a replacement for the physical
balance check.

For elemental input, the production normalization rule is also explicit:
`s = sum(b) = 1.440533904e9`, `b_hat = b/s`, `n_hat = n/s`, and
`y_hat = log(n_hat) = log(n) - log(s)`. The emitted witness deliberately uses
raw physical `b` and `y`; row scaling `W` is printed separately and cannot
replace physical balance checks.

The direct witness coordinates are always `y = log(n)`. Plain mole coordinates
`n` are not passed to RST. The normalized recovery representation uses
`y_hat = y - log(s)`, evaluates `n_hat = exp(y_hat)`, and reconstructs the
published physical vector as `n = s * n_hat`. Under the default facade policy,
the physical solve is attempted first and this normalized representation is a
recovery route after a classified numerical failure, not an unconditional first
step.

There is an important current adapter distinction. Legacy backends consume the
runner's scaled residual/Jacobian closures. RST currently receives the separate
`RstPreparedProblem`, whose symbolic callbacks are raw; the generic scaled `f/j`
closures are not passed into that RST call. Consequently the direct RST matrix
in this witness is a raw-system characterization. The common validation gate
then computes both views, and its `residual_l2_norm` is the scaled view when
scaling is enabled. Making RST consume a parameterized scaled symbolic problem
is a separate integration task, not something this witness silently assumes.

The molecular oracle is solved at the actual `1e8` physical scale through the
production `EquilibriumCalculator` molecular-input facade with an explicit
`SolverPolicy::Single(SolverBackend::Legacy(Solvers::TR))` policy. It is
independent of the direct standalone RST callback invocation and independent
of the RST algorithms under test, while still exercising the real facade,
phase-control, validation, and publication boundaries.
The current debug run reports `raw_residual_l2 = 1.24957e-6`,
`scaled_residual_l2 = 2.42217e-14`, and
`max_abs(A*n-b) = 1.19209e-7`, so it passes the proposed physical thresholds.
Both `gas` and `graphite` are active; the final graphite minimum TPD is
`+1.87993e-9`, accepted by the phase-control stability tolerance. Phase
control is therefore an active-set/complementarity condition, not an extra
equality that can be added to this fixed-topology graph without redesigning
the formulation.

The ideal phase-stability value is
`TPD(x) = sum_i x_i * (g_reduced[i] + R*T*ln(x_i))`, with zero-support terms
omitted, and is expressed in `J/mol`. At `700 K` the default temperature-scaled
thresholds are `dg_create = -5.820123832707268e-3 J/mol` and
`dg_keep = 5.820123832707268e-5 J/mol`. An inactive candidate is activated when
`minimum_tpd < dg_create`; an active phase is removed only when its total is
below `phase_eps` and `minimum_tpd > dg_keep`. The default `phase_eps` is
`1e-30 mol`. A stable active phase therefore need not have strictly negative
TPD; it must be inside the hysteresis keep threshold.

The singular values contain one exact zero at the stated numerical threshold,
so the `18 x 18` reduced Jacobian has rank `17`. This remains an intentional
rank-deficiency safety regression. It is not yet the full physical RST gate:
that gate must wait for a full-rank formulation carrying normalization and
the outer phase-control acceptance contract.

### Recorded debug characterization

The current direct, unnormalized `18 x 18` graph reports:

| RST method | Result |
|---|---|
| `rst_levenberg_marquardt` | typed `SingularJacobian` error |
| `rst_minpack_levenberg_marquardt` | `Stagnation`, residual `1.2486e2` |
| `rst_nielsen_levenberg_marquardt` | `StepTooSmall`, residual `1.2486e2` |
| `rst_trust_region_levenberg_marquardt` | `MaxIterations` at 500, residual `5.2978e1` |
| `rst_damped_newton` | typed `SingularJacobian` error |
| `rst_powell_dogleg` | typed `NumericalBreakdown/SingularMatrix` error |

This is a solver-development witness, not a new equilibrium benchmark. Its
purpose is to let RST investigate conditioning and scaling while preserving
the exact TP-1907 symbolic payload that triggered the KiThe facade failure.

### Full acceptance gate status

The future non-ignored gate is specified but intentionally not marked green.
For each of RST Minpack LM, Nielsen LM, and Trust-Region LM it must require:

- `raw_residual_l2 < 1e-5` and separately reported `scaled_residual_l2`;
- `max_abs(A*n-b) <= 1e-6`;
- component relative error `<= 1e-5` or absolute error `<= 1e-10`;
- mole-fraction error `<= 1e-8`;
- finite, non-negative published moles, exact component/phase ordering, and
  oracle topology/phase-control acceptance.

The current isolated methods do not satisfy this gate because the transferred
reduced graph is rank deficient. The gate must be enabled only after RST can
solve the augmented physical formulation; weakening tolerances or accepting a
scaled residual alone is explicitly out of scope.

The distinction from the current common KiThe gate is deliberate: when row
scaling is enabled, `EquilibriumCandidateReport.residual_l2_norm` is the norm
of `F_scaled` and is the value currently compared with the backend tolerance;
`raw_residual_l2_norm` is recorded diagnostic evidence. Physical element
balances are independently recomputed as `A*n-b` and currently use the
configured absolute-plus-relative tolerance. The strict raw and absolute
physical thresholds above belong to this future TP-1907 transfer gate.

The production-correct pipeline as implemented today is therefore:
physical input -> physical fixed-active problem build -> raw RST symbolic
preparation or scaled legacy closures -> nonlinear backend cascade -> candidate
validation (`F_raw` and `F_scaled`, plus physical `A*n-b`) -> phase-control TPD
and active-set transition/restart -> publication. If the physical attempt has a
classified numerical failure at a material extensive scale, normalized discovery
(`b_hat`, `n_hat`, `y_hat`) is inserted, followed by physical retry or the typed
physical publication boundary. A future RST gate must state explicitly whether
its symbolic callback is raw or scaled; it must not infer this from the name
`residual_row_scale_w`.

### 3a.1a. Strict KiThe/RST parity export from the accepted facade route

#### Story name

TP-1907 upper `ElementInventory` must export a lossless enough accepted KiThe
state for componentwise RST residual/Jacobian parity.

#### Test and command

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_emits_paste_ready_rst_parity_fixture`

```powershell
cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_emits_paste_ready_rst_parity_fixture `
  --no-default-features -- --include-ignored --exact --nocapture *> tp1907_rst_parity_fixture.txt
```

#### Question / hypothesis

Can the already accepted KiThe upper-inventory facade result be transferred to
RustedSciThe as a deterministic numerical fixture, including its exact
coordinates, equations, scaling, phase state, derivatives, physical checks,
and the actual recovery solve history? The export is diagnostic only and does
not alter solver order or acceptance behavior.

#### Current result

**PASS: export generated and test passed.** The request uses an explicit
`Single(RST Trust-Region LM)` policy. The actual accepted route is:

```text
physical attempt (MaxIterations)
-> extensive normalization
-> normalized RST discovery (Converged)
-> physical retry (StepTooSmall)
-> reconstructed physical publication
```

The accepted backend is `RustedSciThe(TrustRegionLevenbergMarquardt)`, the
inventory normalization scale is `1.44053390400000000e9`, and the captured
normalized discovery reports `18` iterations, `34` residual calls, `16`
Jacobian calls, `18` linear solves, `15` accepted steps, and `3` rejected
steps. Its denormalized final coordinate agrees with facade `y_final` at
`0.00000000000000000e0` in the observed run.

The fixture recomputes the accepted physical state on the same `18 x 18`
prepared RST graph. Observed values were:

```text
recomputed raw residual L2       = 1.19255846588116014e-6
recomputed scaled residual L2    = 2.81236124142565315e-14
exp(y_final) max balance error    = 1.19209289550781250e-6
published moles max balance error = 1.19209289550781250e-7
final graphite minimum TPD        = 1.90902937902137637e-9
```

Both `raw_jacobian` and `scaled_jacobian` are printed as full matrices, along
with `raw_residual`, `scaled_residual`, `jt_f_raw`, `jt_f_scaled`, `A`, `b`,
species/phase ordering, bounds, parameter values, and full-precision
`y_final`, `exp(y_final)`, and published moles. The published-moles balance is
the facade publication check; the separate `exp(y_final)` balance is retained
so RST can see any floating-point log/exp round-trip effect directly.
In this run the `exp(y_final)` value is slightly above the proposed strict
`1e-6` transfer limit (`1.19209289550781250e-6`), while the published vector
passes. That is deliberately exposed as a parity diagnostic, not hidden by
substituting published moles for `exp(y_final)`.

#### Trace contract

The actual normalized discovery `SolveResult.history` is captured from the
facade call, not from a substitute hand-built solve. RST `0.4.17` exposes only
iteration, post-step residual norm, step norm, and accepted flag. Therefore the
fixture prints `lambda`, regularization form, trial `alpha`, exact internal
trial residual-before/after, and per-iteration termination as `None`, and
prints the final termination separately. Filling those fields is a separate
RST telemetry TODO; inferring them from aggregate counters would invalidate
parity evidence.

### 3a.1b. Strict symbolic string round-trip gate

#### Story name

The equations exported by the accepted KiThe TP-1907 route must be a faithful
public RST string fixture, not merely a human-readable rederivation.

#### Test and command

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_symbolic_string_round_trip_matches_original_graph`

```powershell
cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_symbolic_string_round_trip_matches_original_graph `
  --no-default-features -- --include-ignored --exact --nocapture
```

#### Question / hypothesis

Does the exact `equations`, ordered variables, and parameter binding returned by
`symbolic_payload_for_test()` produce the same residual and Jacobian when
reconstructed through the public RST
`SymbolicNonlinearProblem::from_strings_with_options(...)` path at the exact
accepted `y_final`? The gate requires componentwise residual and Jacobian drift
no greater than `1e-12`, and prints both L2 norms, the first differing row or
Jacobian element, and an FNV-1a-64 hash of the complete equation list.

#### Current result

**RED, localized before numerical comparison.** RST `0.4.17` panics while
parsing the first exported equation with `Invalid expression format`. The
complete equation-list hash is `FNV-1a-64 = 0x6dc514792152d82c`. The
offending fragment is:

```text
... - -0.12027235504272604 * (...)
```

The serializer emits subtraction of a negative coefficient as `a - -b`, while
the current RST parser does not accept that unary-minus form. The test therefore
does not yet reach residual or Jacobian evaluation; no numerical parity claim
can be made from this run. The failure is deterministic and belongs to the
serializer/parser contract, not to the TP-1907 nonlinear convergence issue.

The test remains ignored until the RST string path accepts a semantically
faithful canonical form. Its assertion is still strict and must not be replaced
by tolerance relaxation or an ad-hoc equation rewrite. The emitted diagnostic
also prints the complete first failing equation and its stable equation hash so
the RST test can reproduce the parser defect exactly.

## 3a.2. TP-1907 analytic versus symbolic Jacobian contract

### Story name

The real TP-1907 upper-scale failure must not be caused by a disagreement
between KiThe's analytical formulation and RST's lambdified symbolic graph.

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_analytic_and_lambdified_symbolic_contracts_agree`

### What it checks

- evaluates the canonical analytical residual and Jacobian;
- evaluates the identical prepared parameterized graph through RST's
  lambdified residual and symbolic Jacobian;
- compares both representations at the `1e8` initial point and a deterministic
  nearby perturbation;
- uses relative error assertions, so large conservation rows are not hidden by
  a small absolute tolerance.

### Command

```powershell
cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_analytic_and_lambdified_symbolic_contracts_agree `
  --no-default-features -- --exact --nocapture
```

### Recorded debug result

```text
TP-1907 analytic-vs-symbolic upper scale
max_residual_abs=2.842170943040401e-14
max_residual_relative=9.770375732637176e-16
max_jacobian_abs=2.632681394231215e-15
max_jacobian_relative=2.632681394231215e-15
test result: ok; 1 passed; 0 failed
```

The residual and Jacobian representations agree to floating-point roundoff.
The upper-scale RST issue therefore remains a nonlinear conditioning/step
acceptance problem, not a symbolic differentiation or lambdification defect.
The backend matrix still supplies the LM comparison (`legacy-lm` versus RST
LM); RST currently has no plain undamped NR method, so `rst-damped-newton` is
only a Newton-family comparison, not an algorithm-identical `NR versus NR`
pair.

## 3b. ElementInventory route equivalence at production scale

### Story name

The molecular, elemental, and formal routes must describe the same physical
equilibrium, including through the public facade and under extensive P,H
scaling.

### Tests

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_frozen_pt_molecular_elemental_and_formal_routes_match`

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_frozen_ph_molecular_elemental_and_formal_routes_reuse_target`

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together`

`Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::elemental_calculator_ph_is_extensive_under_inventory_and_enthalpy_scaling`

### Hypothesis / question

Does `ElementInventory` preserve the molecular equilibrium state through the
P,T and P,H facades, and does joint scaling of inventory and total enthalpy
preserve the intensive P,H solution?

### Acceptance criteria

- `ElementInventory` and formal elemental input reach the same physical state
  as the molecular TP-1907/TP-1906 route;
- the public P,T facade preserves that equivalence on the regular production
  scales; the upper-scale per-backend gate is recorded separately in section
  3a;
- P,H preserves the recovered temperature and target enthalpy when inventory
  and total enthalpy are scaled together;
- route comparisons are made against the same molecular oracle, with residual,
  conservation, topology, and target-enthalpy assertions.

### Release commands

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_frozen_pt_molecular_elemental_and_formal_routes_match --no-default-features -- --include-ignored --exact --nocapture

cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_frozen_ph_molecular_elemental_and_formal_routes_reuse_target --no-default-features -- --include-ignored --exact --nocapture

cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together --no-default-features -- --include-ignored --exact --nocapture
```

### Recorded release result

```text
TP-1907 P,T elemental/formal routes at T=700 K:
  b=[0.05364815999999999, 1.00182736, 2.0, 8.9462088, 2.4036547199999996]
  max external absolute discrepancy=4.187077292968633e-4
  residual=9.793148882444655e-14 balance=9.547918011776346e-14
  test result: ok; 1 passed; 0 failed; finished in 0.27s

TP-1907 upper facade backend acceptance matrix: see section 3a,
  RST Trust-Region LM PASS, Legacy LM PASS, Legacy TR PASS;
  RST LM/Minpack LM/Nielsen LM/Damped Newton ACCURACY_FAIL,
  Legacy NR ACCURACY_FAIL.

TP-1906 P,H elemental/formal routes:
  recovered_T=699.508902431 K H_target=-4.200392761e5 J
  H_error=-8.119e-4 residual=1.842e-8 balance=6.242e-9
  test result: ok; 1 passed; 0 failed; finished in 2.72s

TP-1906 P,H extensive inventory and target scaling:
  test result: ok; 1 passed; 0 failed; finished in 2.93s
```

### Conclusion

The hypothesis is supported. All four route/facade stories passed in release;
the elemental and formal routes remain pinned to the molecular oracle, and the
P,H state remains extensive under joint inventory and target scaling.

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

### Recorded release characterization
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

## Frozen-reference software-regression policy

Frozen datasets are immutable source evidence, not a claim that every rounded
external value is a source-uncertainty acceptance standard for KiThe. The
benchmark-local tests therefore keep two independent layers:

- strict internal I1/I2/TPD, conservation, topology, normalization, and
  file-immutability contracts;
- conservative software-regression envelopes around repeatedly observed KiThe
  agreement with the external source.

The latter are deliberately broader than the recorded values and live only in
the relevant test modules, never in frozen JSON or in this document as runtime
input. Exceeding one requires review of the implementation; it does not prove
that the external source or the local physical model is inaccurate. The current
reviewed guards are: IAPWS liquid `max/rms < 2%`; IAPWS ice `max/rms < 2.5%`;
JANAF Boudouard `max |delta G| < 100 J/mol`, `max |delta log10 K| < 0.01`, and
boundary `max/rms < 1.5%`; NASA CEA `|delta T| < 10 K`, total `< 2%`, major
`max/rms < 10%/7.5%`, minor/trace `|dlog10| < 0.1`; STANJAN major `max/rms <
7.5%/5%`, minor/trace `|dlog10| < 0.1/0.15`; and TP-1907 has its separate
topology, TPD-bracket, gas-composition, and graphite guards. Solver identity,
iteration count, and wall time remain recorded characterization only.

Every ignored frozen-reference story also retains at least one executable
assertion. Tables are for inspection; they cannot be the only mechanism that
detects a behavioral regression.

## 33. Frozen IAPWS low-pressure water saturation boundary

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_iapws_tests::i5_iapws_water_liquid_low_pressure_boundary_diagnostic`

### What it checks

- ten frozen external IAPWS SR1-86(1992) liquid/vapor saturation-pressure
  points at 275..320 K;
- an independent `ln(Q)-ln(K)` pressure root, canonical
  `TPD(liquid | gas)` root, and canonical `TPD(gas | liquid)` root from the
  same local resolved water/liquid records;
- strict three-way internal root agreement separately from the external IAPWS
  comparison;
- byte-for-byte immutability of the frozen metadata and numerical rows.

The story uses direct `H2O(g)/H2O(l)` with no carrier. Above saturation the
bounded outer loop validates a liquid-only reduced state and deactivates the
gas phase; below saturation it retains the gas-only state. This also guards
the inactive-gas TPD path needed for complete condensation and its reverse
evaporation route from a liquid-only reference.

### Debug command

    cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_iapws_tests::i5_iapws_water_liquid_low_pressure_boundary_diagnostic --no-default-features -- --ignored --nocapture

### Recorded characterization

| T K | IAPWS Pa | KiThe I1/I2 pH2O Pa | KiThe TPD pH2O Pa | I1 error % | TPD error % | Internal % |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 275 | 698.436 | 708.715 | 708.715 | 1.47178 | 1.47178 | 0.00000 |
| 280 | 991.759 | 1006.182 | 1006.182 | 1.45429 | 1.45429 | 0.00000 |
| 285 | 1388.931 | 1408.765 | 1408.765 | 1.42800 | 1.42800 | 0.00000 |
| 290 | 1919.877 | 1946.655 | 1946.655 | 1.39477 | 1.39477 | 0.00000 |
| 295 | 2621.115 | 2656.653 | 2656.653 | 1.35583 | 1.35583 | 0.00000 |
| 300 | 3536.718 | 3583.115 | 3583.115 | 1.31188 | 1.31188 | 0.00000 |
| 305 | 4719.327 | 4778.942 | 4778.942 | 1.26322 | 1.26322 | 0.00000 |
| 310 | 6231.204 | 6306.588 | 6306.588 | 1.20979 | 1.20979 | 0.00000 |
| 315 | 8145.307 | 8239.078 | 8239.078 | 1.15124 | 1.15124 | 0.00000 |
| 320 | 10546.384 | 10661.025 | 10661.025 | 1.08701 | 1.08701 | 0.00000 |

Maximum external relative error: `1.471778e-2`; RMS external relative error:
`1.318758e-2`; maximum I1/I2-vs-TPD relative discrepancy: `0.0` at printed
precision. The second, two-sided canonical table has a maximum three-way
internal relative discrepancy of `0.0` at printed precision. The external
difference remains diagnostic for IAPWS-versus-local-model review, while the
established `max/rms < 2%` software-regression envelope prevents an unreviewed
implementation drift.


cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_iapws_tests::tests::i5_iapws_water_liquid_low_pressure_boundary_diagnostic --no-default-features -- --ignored --nocapture

running 1 test
T [K]      IAPWS [Pa]    I1/I2 pH2O      TPD pH2O   I1 err [%]  TPD err [%] internal [%]
275.0         698.436       708.715       708.715      1.47178      1.47178      0.00000
280.0         991.759      1006.182      1006.182      1.45429      1.45429      0.00000
285.0        1388.931      1408.765      1408.765      1.42800      1.42800      0.00000
290.0        1919.877      1946.655      1946.655      1.39477      1.39477      0.00000
295.0        2621.115      2656.653      2656.653      1.35583      1.35583      0.00000
300.0        3536.718      3583.115      3583.115      1.31188      1.31188      0.00000
305.0        4719.327      4778.942      4778.942      1.26322      1.26322      0.00000
310.0        6231.204      6306.588      6306.588      1.20979      1.20979      0.00000
315.0        8145.307      8239.078      8239.078      1.15124      1.15124      0.00000
320.0       10546.384     10661.025     10661.025      1.08701      1.08701      0.00000
summary: max external relative error=1.471778e-2, rms I1/I2 external=1.318758e-2, rms TPD external=1.318758e-2, max I1/I2-vs-TPD=0.000000e0
T [K]      I1/I2 pH2O  liq|gas pH2O  gas|liq pH2O    I1-liq [%]    I1-gas [%]   liq-gas [%]
275.0         708.715       708.715       708.715       0.00000       0.00000       0.00000
280.0        1006.182      1006.182      1006.182       0.00000       0.00000       0.00000
285.0        1408.765      1408.765      1408.765       0.00000       0.00000       0.00000
290.0        1946.655      1946.655      1946.655       0.00000       0.00000       0.00000
295.0        2656.653      2656.653      2656.653       0.00000       0.00000       0.00000
300.0        3583.115      3583.115      3583.115       0.00000       0.00000       0.00000
305.0        4778.942      4778.942      4778.942       0.00000       0.00000       0.00000
310.0        6306.588      6306.588      6306.588       0.00000       0.00000       0.00000
315.0        8239.078      8239.078      8239.078       0.00000       0.00000       0.00000
320.0       10661.025     10661.025     10661.025       0.00000       0.00000       0.00000
summary: max three-way internal relative error=0.000000e0
test Thermodynamics::ChemEquilibrium::frozen_reference_iapws_tests::tests::i5_iapws_water_liquid_low_pressure_boundary_diagnostic ... ok

## 34. Frozen IAPWS ice-Ih low-pressure sublimation boundary

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_iapws_ice_tests::tests::i5_iapws_water_ice_ih_low_pressure_boundary_diagnostic`

### What it checks

- eight frozen IAPWS R14-08(2011) section 4 equation (6) sublimation-pressure
  rows for `H2O(s, ice Ih) <=> H2O(g)` at 200..270 K;
- local `H2O(g)` and `H2O(s)` validity before solving, without extending a
  polynomial beyond its declared interval;
- independent I1/I2, `TPD(ice | gas)`, and `TPD(gas | ice)` roots found in
  logarithmic pressure coordinates;
- strict three-way internal agreement, plus an external relative-error,
  RMS-error, and signed-bias characterization;
- byte-for-byte immutability of both frozen JSON files.

The companion non-ignored lifecycle test works relative to the *local*
250 K boundary: three times that pressure must deposit ice; one third must
sublimate it. That separation keeps phase-control evidence independent of the
external IAPWS/model discrepancy.

### Debug command

    cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_iapws_ice_tests::tests::i5_iapws_water_ice_ih_low_pressure_boundary_diagnostic --no-default-features -- --ignored --nocapture

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_iapws_ice_tests::tests::i5_iapws_water_ice_ih_low_pressure_boundary_diagnostic --no-default-features -- --ignored --nocapture

### Recorded debug characterization

| T K | IAPWS Pa | KiThe I1/I2 Pa | TPD ice\|gas Pa | TPD gas\|ice Pa | External % | Internal % |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 | 0.162604 | 0.165326 | 0.165326 | 0.165326 | 1.67407 | 0.00000 |
| 210 | 0.701727 | 0.713380 | 0.713380 | 0.713380 | 1.66060 | 0.00000 |
| 220 | 2.654072 | 2.697812 | 2.697812 | 2.697812 | 1.64803 | 0.00000 |
| 230 | 8.947353 | 9.093507 | 9.093507 | 9.093507 | 1.63349 | 0.00000 |
| 240 | 27.266844 | 27.707033 | 27.707033 | 27.707033 | 1.61437 | 0.00000 |
| 250 | 76.012670 | 77.220056 | 77.220056 | 77.220056 | 1.58840 | 0.00000 |
| 260 | 195.801674 | 198.843580 | 198.843580 | 198.843580 | 1.55356 | 0.00000 |
| 270 | 470.061878 | 477.150577 | 477.150577 | 477.150577 | 1.50804 | 0.00000 |

Debug summary: maximum external relative error `1.674069e-2`, RMS external
relative error `1.610953e-2`, mean signed external error `1.610071e-2`, and
maximum three-way internal relative discrepancy `0.0` at printed precision.
The positive, smoothly decreasing external bias is recorded for model review;
it is not a source-accuracy acceptance threshold. The test also applies a
separate conservative KiThe software-regression envelope of `max/rms < 2.5%`.
Release output is pending.

## 35. Frozen JANAF Boudouard reaction thermochemistry

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_tests::tests::i5_janaf_boudouard_reaction_thermochemistry_diagnostic`

### What it checks

- eleven frozen NIST-JANAF primary species rows at 500..1500 K: CO(g) C-093,
  CO2(g) C-095, and provenance for graphite reference-state C-002;
- JANAF reaction Gibbs energy and `log10(Kp)` derived twice from the frozen
  columns: formation Gibbs energies and tabulated formation `log Kf`;
- the local offline `BoudouardCarbon` fixture, exact `CO/CO2/C(gr)` identity,
  its common validity range, and reaction-space dimensions `gas-only=0`,
  `full=1`;
- local standard `G0(T)` versus JANAF under an explicit `p0 = 100000 Pa`
  convention, with no pressure-boundary or phase-control calculation;
- typed source-pressure provenance: a record may declare
  `standard_state_pressure_pa` or `reference_pressure_pa`; the current local
  NASA records declare neither and are printed as `undeclared`, never silently
  treated as 1 bar or 1 atm;
- byte-for-byte immutability of both frozen JSON files.

The scope is intentionally thermochemistry only. Any later JANAF-derived
Boudouard pressure boundary, TPD, activation/deactivation, or P,H scenario is
a separate validation layer and must not be folded into this diagnosis.

### Debug command

    cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_tests::tests::i5_janaf_boudouard_reaction_thermochemistry_diagnostic --no-default-features -- --ignored --nocapture

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_tests::tests::i5_janaf_boudouard_reaction_thermochemistry_diagnostic --no-default-features -- --ignored --nocapture

### Recorded debug characterization

| T K | JANAF dGr kJ/mol | KiThe dGr kJ/mol | dG error J/mol | JANAF log10K | KiThe log10K | dlog10K |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 500 | -84.111000 | -84.069486 | 41.514 | 8.786844 | 8.782507 | -0.004337 |
| 600 | -66.210000 | -66.164903 | 45.097 | 5.763980 | 5.760054 | -0.003926 |
| 700 | -48.362000 | -48.313393 | 48.607 | 3.608746 | 3.605119 | -0.003627 |
| 800 | -30.592000 | -30.544400 | 47.600 | 1.997414 | 1.994306 | -0.003108 |
| 900 | -12.916000 | -12.868233 | 47.767 | 0.749610 | 0.746838 | -0.002772 |
| 1000 | 4.664000 | 4.712912 | 48.912 | -0.243618 | -0.246172 | -0.002555 |
| 1100 | 22.149000 | 22.198822 | 49.822 | -1.051748 | -1.054113 | -0.002366 |
| 1200 | 39.540000 | 39.592352 | 52.352 | -1.721098 | -1.723377 | -0.002279 |
| 1300 | 56.841000 | 56.897148 | 56.148 | -2.283855 | -2.286111 | -0.002256 |
| 1400 | 74.058000 | 74.116687 | 58.687 | -2.763084 | -2.765273 | -0.002190 |
| 1500 | 91.192000 | 91.254283 | 62.283 | -3.175526 | -3.177694 | -0.002169 |

Debug summary: JANAF's two frozen derivation routes differ by at most
`1.389514e-3 log10K`, compatible with printed-table rounding. KiThe's maximum
energy difference is `62.283 J/mol`; RMS is `51.131 J/mol`; mean signed error
is `50.799 J/mol`. The maximum `|dlog10K|` is `0.004337`. The smooth signed
bias is visible evidence for a future thermochemistry review, not a reason to
modify records or claim source-accuracy acceptance. The test separately applies
the conservative KiThe envelope `|delta G| < 100 J/mol` and `|delta log10 K| <
0.01`. Release output is pending. The diagnostic additionally reports `undeclared` standard-state
pressure provenance for `CO`, `CO2`, and `C(gr)`; this blocks a future JANAF
pressure-boundary comparison until reviewed local source metadata is supplied.

## 36. Frozen JANAF Boudouard P,T pressure boundary and graphite lifecycle

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_boundary_tests::tests::i5_janaf_boudouard_pressure_boundary_diagnostic`

### What it checks

- an external 50/50 `CO/CO2` analytical Boudouard pressure boundary derived
  directly from the frozen JANAF primary rows at 800, 900, and 1000 K:
  `P = 2 * 100000 Pa / Kp(T)`;
- two deliberately independent local boundary routes under the same explicit
  `reference_pressure = 100000 Pa`: scalar I1/I2 `ln(Q)-ln(K)=0` and
  canonical `TPD(C(gr) | gas)=0`, searched in `ln(P)`;
- the structural invariant `gas-only reaction dimension = 0`, full dimension
  `= 1`, and the pressure metamorphic law for `delta nu_gas = -1`;
- local graphite lifecycle at 900 K: `3 * P_local` activates graphite with a
  negative pre-activation TPD, while `P_local / 3` retains inactive graphite
  with a positive stability TPD and accepted scale-aware conservation;
- byte-for-byte immutability of the frozen JANAF JSON and local JSON
  libraries.

JANAF-to-local values are characterization only. The local NASA payloads for
`CO`, `CO2`, and `C(gr)` expose standard-state pressure as `undeclared`, so the
diagnostic prints that fact and does not apply a hidden 1-bar or 1-atm shift.

### Debug command

    cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_boundary_tests::tests::i5_janaf_boudouard_pressure_boundary_diagnostic --no-default-features -- --ignored --nocapture

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_boundary_tests::tests::i5_janaf_boudouard_pressure_boundary_diagnostic --no-default-features -- --ignored --nocapture

### Recorded debug characterization

| T K | JANAF log10K | JANAF P Pa | I1/I2 P Pa | TPD P Pa | I1/I2 external error | TPD external error | Internal error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 800 | 1.997414 | 2011.9462 | 2026.3957 | 2026.3957 | 0.71819% | 0.71819% | 0.00000% |
| 900 | 0.749610 | 35597.5009 | 35825.4580 | 35825.4580 | 0.64037% | 0.64037% | 0.00000% |
| 1000 | -0.243618 | 350467.3784 | 352535.1659 | 352535.1659 | 0.59001% | 0.59001% | 0.00000% |

Debug summary: maximum JANAF-to-local relative difference `7.181865e-3`; RMS
I1/I2 difference `6.516596e-3`; signed mean `6.495230e-3`; maximum local
I1/I2-versus-TPD discrepancy `0.0` at printed precision. This is reproducible
evidence for thermochemistry/provenance review, not an acceptance threshold.
### Recorded release characterization

The same read-only command passed in the optimized profile. The local roots
remained equal at printed precision; the small last-digit difference at 1000 K
is normal floating-point variation between profiles.

| T K | JANAF P Pa | I1/I2 P Pa | TPD P Pa | I1/I2 external error | TPD external error | Internal error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 800 | 2011.9462 | 2026.3957 | 2026.3957 | 0.71819% | 0.71819% | 0.00000% |
| 900 | 35597.5009 | 35825.4580 | 35825.4580 | 0.64037% | 0.64037% | 0.00000% |
| 1000 | 350467.3784 | 352535.1984 | 352535.1984 | 0.59002% | 0.59002% | 0.00000% |

Release summary: maximum external relative difference `7.181865e-3`; RMS
I1/I2 external difference `6.516624e-3`; mean signed difference
`6.495261e-3`; maximum internal root discrepancy `0.0` at printed precision.
The source remains characterization evidence; the test also enforces a separate
conservative KiThe boundary envelope of `max/rms < 1.5%`.

## 37. Local Boudouard P,H graphite lifecycle

### Test

`Thermodynamics::ChemEquilibrium::real_boudouard_ph_tests::p9_i1_i2_i3_i4_boudouard_ph_lifecycle_diagnostic`

### What it checks

- one independent I1/I2 interior Boudouard state for
  `2 CO(g) <=> CO2(g) + C(gr)` at the explicit convention
  `P = p0 = 100000 Pa`, with `H_target` assembled only from local molar
  enthalpy closures;
- I3 graphite appearance from a gas-only inventory, stable high-temperature
  graphite absence, and active-to-inactive graphite disappearance;
- I1/TPD evidence for the gas-only endpoints and I2 composition/enthalpy
  agreement for interior graphite-bearing states;
- a five-point forward continuation sweep using only previous accepted
  production states, including final lifecycle transitions separately from
  all nested-temperature trial phase events;
- byte-for-byte immutability of both local library JSON inputs and frozen
  JANAF Boudouard provenance rows.

The external JANAF data remain I5 thermochemistry characterization only. This
story does not claim a published external P,H equilibrium table.

### Debug command

    cargo test --lib Thermodynamics::ChemEquilibrium::real_boudouard_ph_tests::p9_i1_i2_i3_i4_boudouard_ph_lifecycle_diagnostic --no-default-features -- --ignored --nocapture

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::real_boudouard_ph_tests::p9_i1_i2_i3_i4_boudouard_ph_lifecycle_diagnostic --no-default-features -- --ignored --nocapture

### Recorded debug characterization

Boudouard P,H lifecycle diagnostic | P=p0=100000 Pa
route              | H target J      | T K       | CO mol       | CO2 mol      | C(gr) mol    | transitions | boundary TPD
independent I2     | -2.212100e5 |  700.0000 | +9.453186e-3 | +5.952734e-1 | +5.902734e-1 | -           | chemical=-9.592e-13
appearance I3      | -2.212100e5 |  700.0000 | +9.453186e-3 | +5.952734e-1 | +5.902734e-1 |           1 | -8.013828e4
stable absence     | -9.119399e4 | 1400.0000 | +1.190000e0 | +5.000000e-3 | +0.000000e0 |           0 | +1.046691e4
disappearance      | -9.119399e4 | 1400.0009 | +1.190000e0 | +5.000108e-3 | +0.000000e0 |           1 | +1.046727e4
sweep index | H target J      | seed K    | solved K  | C(gr) mol    | final transitions | trial phase events | preparation
          0 | -1.209431e5 |  650.0000 | 1084.1760 | +9.894344e-2 |                 1 |                 28 | Initial
          1 | -1.190897e5 | 1084.1760 | 1091.2023 | +8.972952e-2 |                 0 |                  1 | Continued
          2 | -1.114502e5 | 1091.2023 | 1126.8952 | +5.331764e-2 |                 0 |                  1 | Continued
          3 | -1.035063e5 | 1126.8952 | 1186.4436 | +2.076303e-2 |                 0 |                  1 | Continued
          4 | -9.119399e4 | 1186.4436 | 1400.0008 | +0.000000e0 |                 1 |                 21 | Continued
test Thermodynamics::ChemEquilibrium::real_boudouard_ph_tests::p9_i1_i2_i3_i4_boudouard_ph_lifecycle_diagnostic ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 2257 filtered out; finished in 0.44s


The forward sweep begins with graphite present and ends carbon-free. Its
reported `trial phase events` include all nested P,H temperature trials; they
must not be interpreted as the final accepted lifecycle-transition count.
Release output is pending.

## 38. Frozen NASA CEA H2/O2 gas-only P,H equilibrium characterization

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_h2_o2_hp_tests::tests::i5_nasa_cea_hp_gas_only_equilibrium_characterization`

### What it checks

- reconstructs the one-kilogram H2/O2 reactant inventory from CEA's O/F mass
  ratio using local molar masses and constructs the local extensive enthalpy
  target at 2000 K;
- resolves exactly the nine gas identities printed by CEA from the offline
  `NASA_gas` library, with no network fallback or species expansion;
- solves the fixed-pressure, fixed-enthalpy problem through production
  `P,H Auto` and retains the route decision, strict residual, elemental
  balance, and immutable local/frozen-file snapshots;
- compares KiThe and the rounded external CEA table by semantic species
  identity, including major/minor/trace classification and absolute, relative,
  and logarithmic amount differences;
- retains `H2O(L)` and `H2O(cr)` as explicit
  `ExternallyAbsentExcluded` rows instead of silently deleting them.

Only condensed water was excluded from the local solve. Gas-phase `H2O`
remains a normal solved component. The external CEA table reports both
condensed rows as exactly zero, but that alone is not the local exclusion
criterion: the bundled `H2O(L)` record ends at 600 K and `H2O(s)` at 273.15 K,
so neither can be evaluated near the 3181 K equilibrium state. Consequently,
this story characterizes the nine-species gas equilibrium and **does not claim
that KiThe locally evaluated TPD and rejected liquid water or ice**. A full
eleven-component phase-stability result remains deferred until reviewed
high-temperature condensed thermochemistry or another explicit applicability
policy exists.

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_h2_o2_hp_tests::tests::i5_nasa_cea_hp_gas_only_equilibrium_characterization --no-default-features -- --ignored --nocapture

### Recorded release characterization

`Auto` first attempted the monolithic route. Its three configured backends
failed, after which the accepted nested-temperature route converged. This is
retained as route evidence rather than hidden behind the successful result.

- dataset: `nasa_cea.h2_o2.hp.scitech_2025.v1`;
- accepted route: `NestedTemperature`;
- fallback: `AllBackendsFailed` after three monolithic attempts;
- pressure and activity reference pressure: `101325 Pa`;
- reactant temperature: `2000 K`;
- local target enthalpy: `3.296826777e6 J`;
- CEA temperature: `3181.230 K`;
- KiThe temperature: `3183.221 K`;
- signed temperature difference: `+1.991 K`;
- CEA total amount: `5.310000000e-2 kgmol/kg`;
- KiThe total amount: `5.313595257e-2 kgmol/kg`;
- residual L2 norm: `6.186e-8`;
- maximum elemental-balance error: `5.193e-8`;
- final phase transitions and nested trial phase events: `0` and `0`.

| CEA identity | Class | CEA kgmol/kg | KiThe kgmol/kg | Absolute difference | Relative difference | dlog10 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| H | Major | 3.5000000e-3 | 3.5704255e-3 | +7.0425509e-5 | +2.012e-2 | +0.009 |
| H2 | Major | 3.2000000e-3 | 3.3277471e-3 | +1.2774709e-4 | +3.992e-2 | +0.017 |
| H2O | Major | 2.0000000e-2 | 2.0175312e-2 | +1.7531231e-4 | +8.766e-3 | +0.004 |
| H2O2 | Trace | 2.2000000e-7 | 2.2530587e-7 | +5.3058691e-9 | +2.412e-2 | +0.010 |
| HO2 | Minor | 6.6400000e-6 | 6.5835776e-6 | -5.6422421e-8 | -8.497e-3 | -0.004 |
| O | Major | 5.2000000e-3 | 5.2895651e-3 | +8.9565125e-5 | +1.722e-2 | +0.007 |
| O2 | Major | 1.2400000e-2 | 1.2554998e-2 | +1.5499845e-4 | +1.250e-2 | +0.005 |
| O3 | Trace | 1.2100000e-8 | 1.2007209e-8 | -9.2790923e-11 | -7.669e-3 | -0.003 |
| OH | Major | 8.7500000e-3 | 8.2110832e-3 | -5.3891681e-4 | -6.159e-2 | -0.028 |
| H2O(L) | ExternallyAbsentExcluded | 0 | excluded | - | - | - |
| H2O(cr) | ExternallyAbsentExcluded | 0 | excluded | - | - | - |

The temperature and gas composition are close to the rounded CEA result, while
strict internal conservation remains independently visible. These differences
are not source-derived acceptance tolerances; the benchmark nevertheless has a
separate conservative KiThe software-regression envelope documented above.

## 39. Frozen Argonne/STANJAN CHON fixed-P,T characterization

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_argonne_stanjan_chon_tests::tests::i5_argonne_stanjan_chon_fixed_pt_equilibrium_characterization`

### What it checks

- retains the complete 16-row published STANJAN output for the
  pentane-methane-air case at `2500 K`, `35 atm`;
- proves that the original `C5H12 + CH4 + O2 + N2` feed and the local
  `4 CH4 + 2 CO2 + 8 O2 + 37.6 N2` feed have the same closed C/H/O/N totals;
- resolves exactly the local 15-species `NASA_gas` universe without network
  fallback or automatic species expansion, and checks each record at `2500 K`;
- runs the normal general fixed-`P,T` production Gibbs formulation with its
  `15 components / 4 elements / rank 4 / 11 reaction directions` structure;
- compares every local component by identity, retains source `C5H12 = 0` as
  `ExternallyZeroReactantNotSolved`, and excludes it from numeric metrics;
- proves local-library and frozen-reference JSON files are byte-for-byte
  unchanged.

The STANJAN values are rounded published evidence. They are deliberately not
renormalized and remains **CharacterizationOnly** with respect to source
accuracy. It nevertheless has a separate conservative KiThe
software-regression envelope; it neither converts one source table into an
external acceptance tolerance nor silently treats absent local pentane
thermochemistry as an incomplete solver universe.

### Debug command

    cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_argonne_stanjan_chon_tests::tests::i5_argonne_stanjan_chon_fixed_pt_equilibrium_characterization --no-default-features -- --ignored --nocapture

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_argonne_stanjan_chon_tests::tests::i5_argonne_stanjan_chon_fixed_pt_equilibrium_characterization --no-default-features -- --ignored --nocapture

### Recorded debug characterization

- dataset: `argonne.stanjan.pentane_methane_air.tp.2500k_35atm.v1`;
- accepted backend: `RustedSciThe(LevenbergMarquardt)`;
- physical pressure / activity reference: `3546375 Pa / 101325 Pa`;
- final phase transitions: `0`;
- local mole-fraction sum: `1.000000000`;
- frozen STANJAN sum: `1.000008107`, retained as printed source rounding;
- residual L2 norm: `1.123e-7`;
- maximum element-balance error: `1.087e-7`;
- major-species maximum relative delta: `3.743e-2` (`NO`);
- minor maximum `|delta log10|`: `4.001e-2`;
- trace maximum `|delta log10|`: `6.029e-2`.

The exact species table is intentionally printed by the ignored command. Its
current differences are useful provenance/model characterization, not
source-accuracy thresholds. The benchmark additionally has the separate
conservative KiThe software-regression envelope documented above.

## 40. Frozen NASA TP-1907 CHON + graphite multiphase topology characterization

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_nasa_tp1907_chon_graphite_source_audit_and_sweep`

### What it checks

- preserves NASA TP-1907 Table 11.3E rows at 680, 700, 720, and 740 K with
  distinct gas and condensed system-total fraction fields;
- resolves an offline reviewed `17 NASA_gas + C(gr) NASA_cond` universe with
  NIST disabled, `5` independent elements, and `13` reaction directions;
- reconstructs the source feed as `CH2 + 1.2 * (O2 + 3.727587 N2 +
  0.0447068 Ar + 0.0015228 CO2)`; published F/A and chemical ER are reported
  as diagnostics rather than inputs;
- starts with gas only and runs ordinary bounded phase control for all four
  rows, so graphite must arise or remain inactive through canonical TPD and
  active-set lifecycle evidence;
- records the external topology sequence active/active/inactive/inactive,
  source/local system fractions, TPD evidence, accepted-state continuation in
  both directions, internal TPD root bracket, and local/frozen JSON immutability;
- keeps zero `H2O(s)` and `H2O(l)` rows as explicit external exclusions rather
  than extrapolating their local thermochemistry to this graphite case.

The selected source supports **topology and characterization**, not a blanket
source-accuracy oracle. The TP-1906 dry-air convention is now explicit; the
conceptual `CH2` basis remains a transparent local inventory model, not a NASA
molecular fuel record. A separate non-ignored regression enforces active /
active / inactive / inactive topology, TPD signs, a `700..720 K` internal root
bracket, major-gas `max/rms < 10%/5%`, graphite relative error `< 30%` at
680 K, and graphite absolute system-fraction error `< 1e-3` at the
near-boundary 700 K point. Those are conservative KiThe regression guards,
not a statement of external source uncertainty. It also proves that a one-point
typed range matches an independent bounded `P,T` solve component-by-component
at every printed temperature. The earlier apparent mismatch compared range
**moles** with direct **system fractions**; both continuation directions below
now print the same system-fraction normalization.

### Debug command

    cargo test --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_nasa_tp1907_chon_graphite_source_audit_and_sweep --no-default-features -- --ignored --nocapture

### Release command

    cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_nasa_tp1907_chon_graphite_source_audit_and_sweep --no-default-features -- --ignored --nocapture

### Recorded source-faithful debug characterization

| T K | NASA C(gr), system | KiThe C(gr), system | NASA topology | KiThe topology | Boundary TPD | Transitions |
| ---: | ---: | ---: | --- | --- | ---: | ---: |
| 680 | 3.420e-3 | 3.866e-3 | active | active | -6.320e2 | 1 |
| 700 | 5.500e-4 | 9.212e-4 | active | active | -1.630e2 | 1 |
| 720 | 0 | 0 | inactive | inactive | +4.849e2 | 0 |
| 740 | 0 | 0 | inactive | inactive | +1.393e3 | 0 |

The isolated-row source-faithful graphite error is about `+13%` at 680 K and
`+67%` at near-boundary 700 K. Major gas values remain close; the small NH3
column has larger relative error. The canonical internal TPD root is
`705.691406 K`, within the external topology-only 700--720 K bracket. The
continued sweep has the expected phase topology in both directions. Its first
point is now independently proven equal to the corresponding direct `P,T`
solution; no numerical discrepancy was present after using a common
system-fraction basis.

### Recorded accepted-state continuation (debug)

| Direction | T K | Graphite status | KiThe C(gr), system | Continuation seed | Final transitions | Nonlinear iterations |
| --- | ---: | --- | ---: | --- | ---: | ---: |
| forward | 680 | active | 3.866045e-3 | no | 1 | 20 |
| forward | 700 | active | 9.212074e-4 | yes | 0 | 5 |
| forward | 720 | inactive | 0 | yes | 1 | 0 |
| forward | 740 | inactive | 0 | yes | 0 | 4 |
| reverse | 740 | inactive | 0 | no | 0 | 9 |
| reverse | 720 | inactive | 0 | yes | 0 | 4 |
| reverse | 700 | active | 9.212078e-4 | yes | 1 | 17 |
| reverse | 680 | active | 3.866092e-3 | yes | 0 | 6 |

The different iteration counts are expected continuation behavior: only an
already accepted state may seed the next point. The active/inactive topology
and normalized graphite amount remain consistent in both traversal directions.
### Release

NASA TP-1906/1907 source reconstruction
  H/C=2.000000 ER=1.250000 alpha=1.20000000
  dry air / O2: N2=3.7275870 Ar=0.0447068 CO2=0.0015228
  diagnostics: F/A=0.084535289 (published 0.084535), chemical ER=1.2500000 (published 1.2496)
NASA TP-1907 Table 11.3E | source-faithful CH2 + dry-air P=101325 Pa
T K    NASA C(gr)    KiThe C(gr)    NASA topology  KiThe topology  TPD              transitions
680.0  3.420000e-3     3.866045e-3     active          Active -6.320353e2      1
  accepted backend=RustedSciThe(LevenbergMarquardt) nonlinear_iterations=20
  component     NASA(system)   KiThe(system)  relative       dlog10       note
  Ar           8.430000e-3    8.429045e-3    -1.133e-4   -4.919e-5   -
  CH4          1.318000e-2    1.284739e-2    -2.524e-2   -1.110e-2   -
  CO           3.720000e-3    3.757883e-3    +1.018e-2   +4.400e-3   -
  CO2          1.371900e-1    1.369329e-1    -1.874e-3   -8.146e-4   -
  H2           3.097000e-2    3.130526e-2    +1.083e-2   +4.676e-3   -
  H2O          9.981000e-2    1.000316e-1    +2.220e-3   +9.632e-4   -
  NH3          6.000000e-5    5.698829e-5    -5.020e-2   -2.237e-2   -
  N2           7.032200e-1    7.027729e-1    -6.359e-4   -2.762e-4   -
  O2           0.000000e0     5.919523e-32   NaN   NaN   -
  C(gr)        3.420000e-3    3.866045e-3    +1.304e-1   +5.324e-2   -
  H2O(s)       0.000000e0     excluded       -            -            published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed
  H2O(l)       0.000000e0     excluded       -            -            published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed
700.0  5.500000e-4     9.212078e-4     active          Active -1.630329e2      1
  accepted backend=RustedSciThe(LevenbergMarquardt) nonlinear_iterations=17
  component     NASA(system)   KiThe(system)  relative       dlog10       note
  Ar           8.420000e-3    8.417517e-3    -2.949e-4   -1.281e-4   -
  CH4          1.245000e-2    1.214559e-2    -2.445e-2   -1.075e-2   -
  CO           5.800000e-3    5.855963e-3    +9.649e-3   +4.170e-3   -
  CO2          1.384800e-1    1.382662e-1    -1.544e-3   -6.710e-4   -
  H2           3.740000e-2    3.777407e-2    +1.000e-2   +4.322e-3   -
  H2O          9.460000e-2    9.475045e-2    +1.590e-3   +6.902e-4   -
  NH3          6.000000e-5    5.770285e-5    -3.829e-2   -1.695e-2   -
  N2           7.022300e-1    7.018113e-1    -5.963e-4   -2.590e-4   -
  O2           0.000000e0     4.375242e-31   NaN   NaN   -
  C(gr)        5.500000e-4    9.212078e-4    +6.749e-1   +2.240e-1   -
  H2O(s)       0.000000e0     excluded       -            -            published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed
  H2O(l)       0.000000e0     excluded       -            -            published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed
720.0  0.000000e0      0.000000e0      inactive        Inactive 4.848804e2       0
  accepted backend=RustedSciThe(LevenbergMarquardt) nonlinear_iterations=9
  component     NASA(system)   KiThe(system)  relative       dlog10       note
  Ar           8.390000e-3    8.386112e-3    -4.635e-4   -2.013e-4   -
  CH4          1.038000e-2    1.023516e-2    -1.395e-2   -6.103e-3   -
  CO           8.360000e-3    8.488646e-3    +1.539e-2   +6.632e-3   -
  CO2          1.379000e-1    1.378787e-1    -1.545e-4   -6.709e-5   -
  H2           4.383000e-2    4.427570e-2    +1.017e-2   +4.394e-3   -
  H2O          9.167000e-2    9.148570e-2    -2.010e-3   -8.740e-4   -
  NH3          6.000000e-5    5.675731e-5    -5.404e-2   -2.413e-2   -
  N2           6.994000e-1    6.991932e-1    -2.956e-4   -1.284e-4   -
  O2           0.000000e0     3.101130e-30   NaN   NaN   -
  C(gr)        0.000000e0     0.000000e0     NaN   NaN   -
  H2O(s)       0.000000e0     excluded       -            -            published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed
  H2O(l)       0.000000e0     excluded       -            -            published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed
740.0  0.000000e0      0.000000e0      inactive        Inactive 1.393323e3       0
  accepted backend=RustedSciThe(LevenbergMarquardt) nonlinear_iterations=9
  component     NASA(system)   KiThe(system)  relative       dlog10       note
  Ar           8.350000e-3    8.346148e-3    -4.613e-4   -2.004e-4   -
  CH4          7.960000e-3    7.804825e-3    -1.949e-2   -8.550e-3   -
  CO           1.145000e-2    1.159900e-2    +1.301e-2   +5.615e-3   -
  CO2          1.365000e-1    1.364524e-1    -3.487e-4   -1.515e-4   -
  H2           5.002000e-2    5.044357e-2    +8.468e-3   +3.662e-3   -
  H2O          8.961000e-2    8.943742e-2    -1.926e-3   -8.372e-4   -
  NH3          5.000000e-5    5.420311e-5    +8.406e-2   +3.505e-2   -
  N2           6.960700e-1    6.958624e-1    -2.982e-4   -1.295e-4   -
  O2           0.000000e0     2.102217e-29   NaN   NaN   -
  C(gr)        0.000000e0     0.000000e0     NaN   NaN   -
  H2O(s)       0.000000e0     excluded       -            -            published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed
  H2O(l)       0.000000e0     excluded       -            -            published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed
accepted-state continuation
  direction=forward
    T= 680.0 K status=Active C(gr)/system=3.866045e-3 continuation=false transitions=1 backend=RustedSciThe(LevenbergMarquardt) iterations=20
    T= 700.0 K status=Active C(gr)/system=9.212074e-4 continuation=true transitions=0 backend=RustedSciThe(LevenbergMarquardt) iterations=5
    T= 720.0 K status=Inactive C(gr)/system=0.000000e0 continuation=true transitions=1 backend=RustedSciThe(LevenbergMarquardt) iterations=0
    T= 740.0 K status=Inactive C(gr)/system=0.000000e0 continuation=true transitions=0 backend=RustedSciThe(LevenbergMarquardt) iterations=4
  direction=reverse
    T= 740.0 K status=Inactive C(gr)/system=0.000000e0 continuation=false transitions=0 backend=RustedSciThe(LevenbergMarquardt) iterations=9
    T= 720.0 K status=Inactive C(gr)/system=0.000000e0 continuation=true transitions=0 backend=RustedSciThe(LevenbergMarquardt) iterations=4
    T= 700.0 K status=Active C(gr)/system=9.212078e-4 continuation=true transitions=1 backend=RustedSciThe(LevenbergMarquardt) iterations=17
    T= 680.0 K status=Active C(gr)/system=3.866092e-3 continuation=true transitions=0 backend=RustedSciThe(LevenbergMarquardt) iterations=6
canonical graphite TPD boundary: T=705.691406 K tpd=5.688e-6 iterations=26 NASA bracket=[700, 720] K
test Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_nasa_tp1907_chon_graphite_source_audit_and_sweep ... ok

## 41. Frozen NIST ThermoML benzene/toluene binary candidate solution

### Tests

`frozen_reference_nist_benzene_toluene_vle_tests::{offline_benzene_toluene_fixture_preserves_state_specific_nasa_provenance, canonical_tpd_recovers_local_raoult_binary_liquid_minimum, bounded_phase_control_activates_the_binary_liquid_only_above_local_bubble_pressure}`

### What it checks

- a reviewed five-row NIST ThermoML P-x subset at `353.15 K`, with no invented
  experimental vapor-composition values;
- local, offline `NASA_gas:C6H6/C7H8` and `NASA_cond:C6H6(L)/C7H8(L)` lookup
  provenance, explicitly with NIST fallback disabled;
- exact independent local-Raoult versus canonical ideal-solution TPD agreement
  for three interior liquid compositions, including the two-component argmin;
- physical pressure sign: gas is stable below its local bubble pressure and
  bounded phase control creates the liquid above it;
- phase-qualified transfer chemical-potential equality plus molecular and
  elemental conservation after activation.

The separate P-x printout compares local NASA ideal-Raoult pressures with the
published NIST values. It is characterization of source/model differences, not
an external accuracy acceptance threshold.

### Release command

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nist_benzene_toluene_vle_tests::i5_nist_benzene_toluene_px_characterization --no-default-features -- --ignored --nocapture
```

### Recorded release characterization

The release run completed successfully. Local NASA ideal-Raoult pressure is
systematically below the frozen NIST P-x values by roughly `7.3..12.2%` over
this selected isotherm. That is source/model characterization, not a TPD or
active-set failure: the strict local I1-to-I3 argmin and lifecycle tests are
separate and pass.

```text
NIST ThermoML benzene/toluene P-x characterization at 353.15 K
source=P-x only; no experimental vapor composition is inferred
    x_B     NIST P kPa  frozen Raoult     local NASA  NASA-NIST %
  0.000        38.5000        38.5000        33.7900     -12.2337
  0.207        50.5000        51.3340        46.0798      -8.7529
  0.500        69.5000        69.5000        63.4754      -8.6684
  0.702        82.7000        82.0240        75.4683      -8.7444
  1.000       100.5000       100.5000        93.1609      -7.3026
test Thermodynamics::ChemEquilibrium::frozen_reference_nist_benzene_toluene_vle_tests::i5_nist_benzene_toluene_px_characterization ... ok
```

## 42. Frozen NIST ThermoML ternary Antoine-gauge VLE candidate

### Tests

`frozen_reference_nist_ternary_vle_antoine_gauge_tests::{canonical_tpd_recovers_both_three_component_boundary_minima_from_the_antoine_gauge, canonical_phase_control_decision_tracks_both_ternary_boundary_directions, i5_nist_ternary_vle_antoine_raoult_temperature_characterization}` and `frozen_reference_nist_ternary_vle_lifecycle_tests::{production_gas_to_liquid_activation_reaches_accepted_ternary_flash_state, production_liquid_to_gas_activation_reaches_accepted_ternary_flash_state, production_liquid_to_gas_lifecycle_deactivates_the_exhausted_liquid_phase, production_temperature_continuation_uses_only_the_previously_accepted_state, i5_nist_ternary_antoine_gauge_production_lifecycle_characterization}`

### What it checks

- three independent frozen NIST WebBook Antoine records, with a common valid
  temperature interval of `335.19..384.66 K` and one explicit `p0 = 1 bar`;
- four genuine interior NIST ThermoML `T/P/x` rows for toluene, ethylbenzene,
  and chlorobenzene. ThermoML provides no vapor composition; `y` is derived
  exclusively by the independent scalar Raoult root;
- a test-only `G0_gas = 0`, `G0_liquid = RT ln(Psat/p0)` gauge that is never
  sent to `SubsData`, production lookup, or the `P,H` path;
- exact local I1/I3 agreement: gas-to-liquid and liquid-to-gas canonical TPD
  minima are zero at the same Raoult boundary, both recover a three-component
  (two-dimensional-simplex) argmin, and `PhaseManager` selects the correct
  activation direction just across the boundary;
- byte-for-byte immutability of the frozen Antoine and ThermoML source files.

The final temperature table is external characterization only. It compares
independent pure-component Antoine/Raoult calculations with the published
mixture boiling data; it is not an external acceptance tolerance and does not
weaken the exact I1/I3 checks above.

### Release command

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_antoine_gauge_tests::tests::i5_nist_ternary_vle_antoine_raoult_temperature_characterization --no-default-features -- --ignored --nocapture
```

### Recorded characterization

The recorded release characterization finds KiThe's independent Antoine and
ideal-Raoult bubble temperatures `0.524..0.854 K` below the four selected
NIST ThermoML points. The RMS temperature difference is `0.7073 K`, or about
`0.14..0.24%` per row. This is good agreement for the deliberately simple,
untuned ideal-mixture characterization, while remaining visibly non-zero.

```text
NIST ThermoML ternary VLE Antoine/Raoult temperature characterization
source=T/P/x only; vapor composition is derived, never read from ThermoML
  x_tol   x_ethyl   x_chl      P kPa   T_NIST K  T_Raoult K  delta K  relative %
  0.126   0.748   0.126     26.660    361.870      361.223   -0.647    -0.1788
  0.126   0.748   0.126     53.330    382.790      382.266   -0.524    -0.1369
  0.334   0.333   0.333     26.660    355.690      354.836   -0.854    -0.2401
  0.334   0.333   0.333     53.330    376.390      375.630   -0.760    -0.2020
summary: rows=4 RMS delta=0.7073 K max abs delta=0.8540 K
test Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_antoine_gauge_tests::tests::i5_nist_ternary_vle_antoine_raoult_temperature_characterization ... ok
```

### Production lifecycle release command

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_lifecycle_tests::i5_nist_ternary_antoine_gauge_production_lifecycle_characterization --no-default-features -- --ignored --nocapture
```

### Recorded production lifecycle characterization

This is a separate `P,T` lifecycle story over the same frozen gauge. It uses
the ordinary bounded `PreparedPhaseControlRunner`, not a ternary-specific
solver. Rachford-Rice remains an independent phase-regime and composition
oracle. The central `53.33 kPa` inventory reaches the same two-phase state
from gas-only and liquid-only starts; at the independently selected dew-side
point, gas activation is followed by liquid disappearance. The liquid TPD is
near zero for the accepted two-phase state and positive once the stable
gas-only state is reached.

```text
ternary Antoine-gauge production lifecycle characterization
  story                         T K      flash      initial final  beta       TPD_gas    TPD_liquid transitions
  gas -> liquid activation      375.130  TwoPhase gas     gas+liq  0.37107          -     1.22e-12           1
  liquid -> gas activation      375.130  TwoPhase liquid  gas+liq  0.37107          -    -1.51e-12           1
  liquid -> gas disappearance   377.430  AllVapor liquid  gas     1.00000          -       2.96e0           2
test Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_lifecycle_tests::i5_nist_ternary_antoine_gauge_production_lifecycle_characterization ... ok
```

### Fixed-inventory full lifecycle release command

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_full_lifecycle_tests::i5_nist_ternary_fixed_inventory_full_lifecycle_characterization --no-default-features -- --ignored --nocapture
```

### Recorded fixed-inventory lifecycle characterization

This stricter P11.10 story fixes `P = 53.33 kPa` and
`z = [0.334, 0.333, 0.333] mol`.  The production runner is compared point by
point to an independent scalar Antoine/Rachford-Rice oracle.  Both directions
enter and leave the same two-phase interval with exactly two accepted
transitions and no accepted chatter.  At interior points, the largest oracle
differences are `4.146e-12` for vapor fraction, `5.969e-13` for liquid
composition, and `9.150e-13` for vapor composition.  Maximum conservation
error is `1.341e-13 mol`; maximum cross-phase chemical-potential mismatch is
`5.639e-11 J/mol`.

The final line is genuine production hysteresis evidence, not an artificial
threshold: at one physical `P,T,z`, the gas-candidate TPD lies inside the
default create/keep band.  Thus a fresh inactive gas candidate stays absent,
whereas accepted active history retains the two-phase topology.  Outside that
band, both histories select the same topology.

```text
ternary Antoine-gauge fixed-inventory full lifecycle
P=53330 Pa, z=[0.334, 0.333, 0.333] mol, p0=100000 Pa
forward:
  T K      route       oracle     topology  beta      transitions TPD_gas   TPD_liquid
  372.630  initial    AllLiquid liquid     0.00000           0    2.93e2           -
  375.130  continued  AllLiquid liquid     0.00000           0    4.87e1           -
  376.130  continued  TwoPhase gas+liquid  0.12450           1         -    3.34e-12
   component        x_RR       x_KiThe      y_RR       y_KiThe
   toluene         0.313649   0.313649   0.477115   0.477115
   ethylbenzene    0.345522   0.345522   0.244937   0.244937
   chlorobenzene   0.340828   0.340828   0.277948   0.277948
  377.352  continued  TwoPhase gas+liquid  0.44756           0         -    4.31e-12
   component        x_RR       x_KiThe      y_RR       y_KiThe
   toluene         0.265416   0.265416   0.418655   0.418655
   ethylbenzene    0.377276   0.377276   0.278349   0.278349
   chlorobenzene   0.357308   0.357308   0.302996   0.302996
  378.574  continued  TwoPhase gas+liquid  0.82215           0         -    6.37e-12
   component        x_RR       x_KiThe      y_RR       y_KiThe
   toluene         0.219420   0.219420   0.358786   0.358786
   ethylbenzene    0.411643   0.411643   0.315988   0.315988
   chlorobenzene   0.368937   0.368937   0.325226   0.325226
  379.574  continued  AllVapor gas        1.00000           1         -      4.94e1
  382.074  continued  AllVapor gas        1.00000           0         -      2.96e2
reverse:
  T K      route       oracle     topology  beta      transitions TPD_gas   TPD_liquid
  382.074  initial    AllVapor gas        1.00000           0         -      2.96e2
  379.574  continued  AllVapor gas        1.00000           0         -      4.94e1
  378.574  continued  TwoPhase gas+liquid  0.82215           1         -    2.90e-12
   component        x_RR       x_KiThe      y_RR       y_KiThe
   toluene         0.219420   0.219420   0.358786   0.358786
   ethylbenzene    0.411643   0.411643   0.315988   0.315988
   chlorobenzene   0.368937   0.368937   0.325226   0.325226
  377.352  continued  TwoPhase gas+liquid  0.44756           0         -    3.97e-12
   component        x_RR       x_KiThe      y_RR       y_KiThe
   toluene         0.265416   0.265416   0.418655   0.418655
   ethylbenzene    0.377276   0.377276   0.278349   0.278349
   chlorobenzene   0.357308   0.357308   0.302996   0.302996
  376.130  continued  TwoPhase gas+liquid  0.12450           0         -    4.23e-12
   component        x_RR       x_KiThe      y_RR       y_KiThe
   toluene         0.313649   0.313649   0.477115   0.477115
   ethylbenzene    0.345522   0.345522   0.244937   0.244937
   chlorobenzene   0.340828   0.340828   0.277948   0.277948
  375.130  continued  AllLiquid liquid     0.00000           1    4.87e1           -
  372.630  continued  AllLiquid liquid     0.00000           0    2.93e2           -
hysteresis: T=375.629568833 K gas-TPD=-1.545963e-3 J/mol band=(-3.123158e-3, 3.123158e-5) inactive=[false, true] active-history=[true, true]
summary: max|delta beta|=4.146e-12 max|delta x|=5.969e-13 max|delta y|=9.150e-13 max conservation=1.341e-13 max|delta mu|=5.639e-11 J/mol
accepted transitions: forward=2 reverse=2 accepted chatter=0
test Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_full_lifecycle_tests::i5_nist_ternary_fixed_inventory_full_lifecycle_characterization ... ok
```

## 43. Frozen NASA TP-1906/1907 CHON + graphite P,H lifecycle

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_nasa_tp1906_tp1907_chon_graphite_ph_production_lifecycle_characterization`

### What it checks

- a source-faithful TP-1906/1907 CHON plus graphite inventory at one
  atmosphere, with a fixed total inventory mass of `179.957703644 g`;
- P,T preflight enthalpy characterization at `680..740 K`, using the local
  solution enthalpy only to construct independent P,H targets;
- isolated P,H solves from the deliberately distant common `900 K` seed;
- graphite topology on both sides of the `700..720 K` boundary;
- forward and reverse P,H continuation, where every next point may use only
  the previously accepted physical state.

The P,T enthalpy differences of about `1.0..1.5 J/g` are source/model
characterization, not P,H residuals. The P,H targets themselves close to
`|scaled H error| < 5.4e-9`; KiThe recovers the expected graphite topology at
all four points and retains it under both continuation directions.

### Release command

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_nasa_tp1906_tp1907_chon_graphite_ph_production_lifecycle_characterization --no-default-features -- --ignored --nocapture
```

### Recorded release characterization

```text
NASA TP-1906/1907 CHON+graphite P,H characterization
P=101325 Pa inventory_mass_g=179.957703644
P,T enthalpy-reference preflight:
  T source K    h source J/g    h local J/g    delta J/g   graphite expected/local
     680.000    -2375.500000   -2374.503837     0.996163     true/true
     700.000    -2334.100000   -2333.036758     1.063242     true/true
     720.000    -2291.400000   -2289.896706     1.503294    false/false
     740.000    -2247.100000   -2245.621263     1.478737    false/false
isolated P,H cases (common numerical seed=900 K):
  T source  H target J     T KiThe  delta T  graphite expected/local  H scaled  route/trials
  680.000   -4.274895e5    679.498   -0.502     true/true     -5.334e-9 NestedTemperature/26
  700.000   -4.200393e5    699.509   -0.491     true/true      3.590e-9 NestedTemperature/28
  720.000   -4.123551e5    719.309   -0.691    false/false     1.235e-9 NestedTemperature/28
  740.000   -4.043830e5    739.339   -0.661    false/false    -4.645e-9 NestedTemperature/28
forward accepted P,H continuation:
  H target J       seed/continued  T KiThe  graphite  accepted transitions  trial events
    -4.274895e5 Initial  679.498     true         23            23
    -4.200393e5 Continued  699.509     true          8             8
    -4.123551e5 Continued  719.309    false         28            28
    -4.043830e5 Continued  739.339    false          1             1
reverse accepted P,H continuation:
  H target J       seed/continued  T KiThe  graphite  accepted transitions  trial events
    -4.043830e5 Initial  739.339    false          2             2
    -4.123551e5 Continued  719.309    false          5             5
    -4.200393e5 Continued  699.509     true         23            23
    -4.274895e5 Continued  679.498     true          1             1
test Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_nasa_tp1906_tp1907_chon_graphite_ph_production_lifecycle_characterization ... ok
```

## 44. Test-only exact extensive-normalization recovery

### Tests

`frozen_reference_nasa_tp1907_chon_graphite_tests::{i5_tp1907_large_inventory_extensive_normalization_matrix, i5_tp1907_extensive_normalization_scale_metamorphism}` and `frozen_reference_nasa_tp1906_chon_graphite_ph_tests::i5_tp1906_extensive_ph_normalization_recovery_matrix`

### What it checks

- a large physical input is transformed only as `n0_tilde = n0 / s`, where
  `s = sum(n0_i > 0)` is derived from that request itself; the normalized
  internal inventory is one mole;
- `P`, `T`, reference pressure, standard-state thermochemistry, nonlinear
  backend settings, TPD thresholds, and hysteresis policy are not modified;
  the public physical trace floor and `phase_eps` are instead represented as
  their explicit normalized values `/s`, while a relative trace fraction stays
  dimensionless and only its physical absolute cap is transformed;
- ordinary fresh solver initialization is used for the normalized request;
  no accepted reference solution, external equilibrium composition, or
  continuation seed participates in its construction;
- on acceptance, component moles are reconstructed as `n = s * n_tilde` and
  compared only afterwards with the diagnostic transformed-answer oracle B;
- for `P,H`, total target enthalpy is transformed by the same factor,
  `H_target_tilde = H_target / s`, and internal normalized-Joule error is
  reported separately from reconstructed physical-Joule error.

### Statement established by these tests

Within the current ideal NASA TP-1906/1907 CHON + graphite formulation,
uniform extensive scaling is an exact physical symmetry, but large absolute
inventory can place the ordinary fresh nonlinear formulation outside its
numerical basin. The same physical problem, rebuilt with an order-one input
inventory and solved through the ordinary fresh path, reconstructs the
accepted large-scale physical state.

This is evidence for **`AbsoluteExtensiveScaleConditioning`**. It proves a
different failure mechanism from the already rejected trace-floor and
element-feasible-interior hypotheses. These tests established the semantic
boundary before routing was selected. The later production route uses the same
request-derived transform only after a classified numerical failure and keeps
its original failure plus publication route in typed evidence.

### Release commands

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_tp1907_large_inventory_extensive_normalization_matrix --no-default-features -- --ignored --exact --nocapture

cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_tp1907_extensive_normalization_scale_metamorphism --no-default-features -- --ignored --exact --nocapture

cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_extensive_ph_normalization_recovery_matrix --no-default-features -- --ignored --exact --nocapture
```

### Recorded development characterization

The following is a debug characterization. The same commands above still need
to be recorded in release mode before this becomes release evidence.

For the `P,T` `10^4` recovery matrix, ordinary large fresh solves remain
rejected while normalized ordinary fresh solves accept. The reconstructed state
matches oracle B at all three physical topologies, including gas-only `720 K`.

```text
TP-1907 exact extensive-normalization matrix
T K    input factor  route  internal total  topology       status   max n/B error  max gas-x error  TPD error
680.0         1.0e4  N          1.000000e0  gas+graphite  OK           3.205e-10        5.274e-16  -
700.0         1.0e4  N          1.000000e0  gas+graphite  OK           3.167e-11        6.939e-16  -
720.0         1.0e4  N          1.000000e0  gas            OK          2.122e-14        1.110e-16   2.486e-6
```

The `P,T` metamorphic matrix keeps the internal inventory at one mole over
five physical scales. Under the shared representation boundary, the largest
reconstructed-mole error is `2.907e-8` and the largest gas-fraction error is
`4.047e-10`. The graphite TPD comparison remains `2.486e-6` relative at all
720 K rows and stays positive.

```text
TP-1907 extensive-normalization scale metamorphism
T K    physical factor  internal total  topology       max n/B error  max gas-x error  TPD error
680.0          1.0e-4      1.000000e0  gas+graphite     3.942e-10        2.220e-16  -
680.0          1.0e-2      1.000000e0  gas+graphite     6.609e-14        2.914e-16  -
680.0           1.0e0      1.000000e0  gas+graphite     3.581e-12        3.469e-16  -
680.0           1.0e2      1.000000e0  gas+graphite      1.060e-9        1.110e-16  -
680.0           1.0e4      1.000000e0  gas+graphite     3.205e-10        5.274e-16  -
700.0          1.0e-4      1.000000e0  gas+graphite     1.313e-12        1.804e-16  -
700.0          1.0e-2      1.000000e0  gas+graphite     1.405e-12        9.714e-17  -
700.0           1.0e0      1.000000e0  gas+graphite      2.907e-8        8.327e-17  -
700.0           1.0e2      1.000000e0  gas+graphite     1.097e-10        1.110e-16  -
700.0           1.0e4      1.000000e0  gas+graphite     3.167e-11        6.939e-16  -
720.0          1.0e-4      1.000000e0  gas               1.676e-8        4.047e-10   2.486e-6
720.0          1.0e-2      1.000000e0  gas               1.144e-10        2.759e-12   2.486e-6
720.0           1.0e0      1.000000e0  gas               1.634e-13        4.122e-15   2.486e-6
720.0           1.0e2      1.000000e0  gas               1.289e-14        2.220e-16   2.486e-6
720.0           1.0e4      1.000000e0  gas               2.122e-14        1.110e-16   2.486e-6
```

The `P,H` story uses the ordinary nested temperature route. The original
physical `10^4` fresh request remains a required failure witness when recovery
is explicitly disabled; the normalized request recovers the base temperature
and graphite topology.

```text
NASA TP-1906/1907 P,H exact extensive-normalization recovery
T source K  physical factor  internal total  T normalized  delta T K  graphite  max n/base err  H internal J  H physical J  status
   700.000            1.0e4      1.000000e0    699.508908    0.000e0      true        1.462e-9      4.395e-5       2.737e0  OK
   720.000            1.0e4      1.000000e0    719.309220    0.000e0     false       1.263e-11     -3.590e-5      -2.236e0  OK
```

### Recorded release characterization

The corrected tests were rerun in the release profile. The P,T matrix
accepted eight of nine explicitly selected backends; the single rejected
backend remains an honest method-specific diagnostic and does not invalidate
the matrix because the test requires at least one accepted backend.

```text
P,T typed T-range matrix: species=20, points=3
Backend              Status    Total ms  Mean ms  Worst ms
rst_lm               OK          6.1794   2.0598    2.1178
rst_minpack_lm       OK         12.4264   4.1421   10.3735
rst_nielsen_lm       OK        124.3239  41.4413   59.9004
rst_trust_region_lm  OK          7.9978   2.6659    6.1642
rst_powell_dogleg    OK         11.3534   3.7845    6.1882
rst_damped_newton    FAILED             --        --
legacy_lm            OK          0.5295   0.1765    0.2206
legacy_nr            OK          0.3752   0.1251    0.1606
legacy_tr             OK          0.6023   0.2008    0.3305
summary: attempted=9 successful=8 failed=[rst_damped_newton]
```

`rst_damped_newton` failed at the first point after one iteration with a
`Stagnation` termination and could not discover a normalized basin. This is
retained as backend evidence, not converted into a false success.

The P,H matrix also passed in release. It preserved the recovered temperature
and graphite topology at both probes; the normalized internal inventory was
exactly one mole.

```text
P,H extensive-normalization matrix: factor=1.0e4
T source K  T normalized K  delta T K  graphite  max n/base err  H internal J  H physical J  status
700.000      699.508908       0.000e0    true        1.462e-9       4.395e-5       2.737e0  OK
720.000      719.309220       0.000e0   false        1.263e-11      -3.590e-5      -2.236e0  OK
test Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_extensive_ph_normalization_recovery_matrix ... ok

test result: ok; 1 passed; 0 failed; finished in 3.14s
```

## 45. Production extensive-normalization recovery

### Tests

`frozen_reference_nasa_tp1907_chon_graphite_tests::i5_tp1907_extensive_pt_scale_matrix_characterization`
and
`frozen_reference_nasa_tp1906_chon_graphite_ph_tests::i5_tp1906_extensive_ph_scale_matrix_characterization`

### What they check

- the canonical public P,T route accepts the real TP-1907 matrix at
  `10^-4, 10^-2, 1, 10^2, 10^4` without changing thermochemistry, solver
  tolerances, trace floors, `phase_eps`, or TPD hysteresis;
- a direct physical solve remains the first attempt; exact normalization is
  authorized only by a typed numerical failure and a materially non-unit
  request-derived inventory scale;
- `ExtensiveNormalizationPolicy::Disabled` preserves the historical physical
  `10^4` `AllBackendsFailed` witness;
- the accepted result publishes physical component moles, phase totals,
  element-balance evidence, topology, and intensive TPD;
- `ExtensiveNormalizationRecoveryEvidence` distinguishes a successful final
  physical retry from publication reconstructed through the audited physical
  boundary and retains the original/retry failures;
- opt-in diagnostics retain and render a typed recovery event without exposing
  the suppressed normalized-discovery lifecycle in solver-coordinate moles;
- nested P,H obtains the same behavior through its canonical inner P,T route,
  scaling the total enthalpy target together with inventory and preserving
  temperature, graphite topology, and physical enthalpy closure;
- every matrix row is now an assertion: printing `FAILED` can no longer leave
  either story test green.

### Statement established by these tests

For the current ideal multiphase formulation, the former `10^4` failures are
an extensive-coordinate conditioning problem, not missing physical roots.
Exact request-derived normalization is a production recovery representation,
not an answer-derived seed. Public results remain in physical units and expose
when recovery was used. At `680/700 K` the unchanged physical retry remains
ill-conditioned, so publication is reconstructed explicitly; at gas-only
`720 K` the recovered basin is accepted again in physical coordinates.

### Release commands

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_tp1907_extensive_pt_scale_matrix_characterization --no-default-features -- --ignored --exact --nocapture

cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_extensive_ph_scale_matrix_characterization --no-default-features -- --ignored --exact --nocapture
```

### Recorded development characterization

Release output has not yet been recorded. The strict debug P,T matrix accepted
all fifteen rows. The `10^4` recovery route was `reconstructed` at 680 and
700 K and `physical-retry` at 720 K. Maximum recovered-mole relative error was
`4.501e-5`; maximum gas-composition error was `3.305e-14`; the largest physical
absolute element-balance error was `4.491e-7 mol` and remained inside the
declared absolute-plus-relative acceptance contract.

The strict nested P,H matrix accepted all ten rows. At factor `10^4`, H700 used
the reconstructed physical boundary and H720 completed the physical retry.
Both retained their unit-scale temperatures to printed precision, their
graphite-present/graphite-absent topologies, and scaled enthalpy errors below
`7e-10`.

## 46. Prepared P,T range normalization recovery

### Test

`frozen_reference_nasa_tp1907_chon_graphite_tests::i5_tp1907_extensive_temperature_range_recovers_large_inventory_transactionally`

### What it checks

- the prepared bounded T-range path remains the primary fast path and invokes
  canonical exact normalization only for a numerically rejected point;
- recovery receives the current accepted log-mole seed and phase set, so it
  does not reset active-phase history or hysteresis semantics;
- every accepted point is published in physical units and retains the same
  topology, graphite stability sign, component amounts divided by scale, and
  element totals as the unit-inventory range;
- the report distinguishes `RecoveryFormulation` from ordinary prepared reuse
  and counts its additional formulation build; its point total includes both
  the rejected prepared attempt and the recovery transaction;
- internal recovery emits no duplicate point progress: the successful range
  has exactly three `PointStarted` and three `PointAccepted` events;
- with normalization explicitly disabled, point zero retains the historical
  `AllBackendsFailed` cause and emits no `PointAccepted`, proving transactional
  failure rather than partial range publication.

### Statement established by this test

Prepared continuation and extensive recovery are complementary rather than
competing routes. At physical scale `10^4`, the first `680 K` point discovers
and publishes a valid physical basin through normalization. The accepted state
then seeds `700 K` and `720 K`, which reuse the physical prepared formulation
without stale recovery evidence. The graphite phase remains present at
`680/700 K` and disappears at `720 K`, exactly as in the unit-scale range.

### Release command

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_tp1907_extensive_temperature_range_recovers_large_inventory_transactionally --no-default-features -- --ignored --exact --nocapture
```

### Recorded development characterization

```text
NASA TP-1907 prepared T-range extensive recovery
  T K    preparation           topology              max n/f err  graphite TPD  route
   680.0  RecoveryFormulation  gas+graphite             6.351e-13     1.6171e-9  recovery
   700.0  ReusedFormulation    gas+graphite             2.531e-13     1.7671e-9  direct
   720.0  ReusedFormulation    gas                      2.402e-14      4.8488e2  direct
```

The compact output above is development characterization; the release result
is recorded below.

### Recorded release characterization

The release run completed successfully in `1.59 s`:

```text
NASA TP-1907 prepared T-range extensive recovery
  T K    preparation           topology              max n/f err  graphite TPD  route
   680.0  RecoveryFormulation  PhaseId(Some("gas"))+PhaseId(Some("graphite"))    6.456e-13     1.6025e-9  recovery
   700.0  ReusedFormulation    PhaseId(Some("gas"))+PhaseId(Some("graphite"))    4.004e-14     1.6325e-9  direct
   720.0  ReusedFormulation    PhaseId(Some("gas"))                              1.738e-14      4.8488e2  direct
test Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_tp1907_extensive_temperature_range_recovers_large_inventory_transactionally ... ok

test result: ok; finished in 1.59s
```

The release evidence confirms that recovery is paid only at the first point:
the normalized discovery is published back in physical units, the following
points use the accepted physical continuation state, and graphite disappears
at `720 K` without a second recovery formulation.

## 47. Transactional P,H target-range extensive recovery

### Test

`frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_extensive_ph_target_range_recovery_is_transactional`

### What it proves

- A factor `10^4` physical P,H range can recover its first numerical basin
  through an exactly normalized equivalent.
- The accepted normalized state is reconstructed before publication; public
  moles, enthalpy, topology, and provenance remain physical.
- A physical-coordinate retry is attempted after normalized discovery, and its
  outcome is retained in recovery evidence when reconstruction is required.
- The first point remains `Initial`; the next point is `Continued` from the
  preceding accepted physical state. Recovery is never inherited as a solver
  mode by the next point.
- Range publication remains transactional: a failed recovery cannot publish a
  partial suffix or seed a later target.

### Development characterization

```text
NASA TP-1906/1907 transactional P,H range normalization recovery
factor=1.0e4 points=2 recoveries=1
    -4.200393e9 Initial T=699.508908 recovery=true transitions=25
    -4.123551e9 Continued T=719.309221 recovery=false transitions=28
```

The story is intentionally ignored and must be recorded in release because it
uses the local NASA fixture and the large-inventory numerical basin.

### Release command

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_extensive_ph_target_range_recovery_is_transactional --no-default-features -- --ignored --exact --nocapture
```

### Recorded release characterization

The release run completed successfully in `8.92 s`:

```text
NASA TP-1906/1907 transactional P,H range normalization recovery
factor=1.0e4 points=2 recoveries=1
    -4.200393e9 Initial T=699.508908 recovery=true transitions=25
    -4.123551e9 Continued T=719.309221 recovery=false transitions=28
test Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_extensive_ph_target_range_recovery_is_transactional ... ok

test result: ok; finished in 8.92s
```

This confirms the same transactional claim under the release build: one
large-inventory point used recovery, while the next point continued from the
accepted physical state and did not inherit normalized mode or evidence.

## 48. Multicomponent continuation rollback feasibility

### Test

`frozen_reference_nasa_tp1906_chon_graphite_rollback_tests::tests::i5_tp1906_multicomponent_continuation_rollback_feasibility_matrix`

### What it proves

This is a feasibility characterization, not yet rollback evidence. It uses
the real TP-1906/1907 CHON + graphite system and checks whether a restricted
iteration budget can fail after at least one accepted continuation point.
Failures at point zero are deliberately classified separately.

### Recorded release characterization

```text
NASA TP-1906/1907 continuation rollback feasibility
budget  status       failing_point  accepted_points  detail
     1  FAILED                   0               0  post_acceptance=false
     2  FAILED                   0               0  post_acceptance=false
     4  FAILED                   0               0  post_acceptance=false
     8  FAILED                   0               0  post_acceptance=false
    16  OK           -                            3  complete
    32  OK           -                            3  complete
    64  OK           -                            3  complete
   100  OK           -                            3  complete
natural post-acceptance failure found: false
```

Verdict: no natural post-acceptance failure was found. The next stage is a
minimal test-only failpoint at the actual transition-before-commit boundary,
followed by clean/failed/fresh route comparison. No production behavior is
changed by this feasibility matrix.

### Release command

```powershell
cargo test --release --lib --no-default-features rollback_feasibility_matrix -- --ignored --nocapture
```

## 49. Transition-boundary rollback evidence

### Test

`prepared_phase_control_runner::tests::transition_failpoint_rolls_back_after_transition_diagnostics`

### What it proves

The test-only failpoint is placed after a genuine TPD-driven phase activation
has changed the local phase set, emitted `TransitionAccepted` diagnostics, and
created a restart seed, but before the restart can be committed. The outer
transaction then restores the previously accepted continuation seed and phase
set. This checks the important distinction between retained diagnostic
evidence and published physical state.

The failpoint is compiled only for tests and is one-shot. It does not change
the production phase-control policy or create a production failure path.

### Local characterization

```text
test prepared_phase_control_runner::tests::transition_failpoint_rolls_back_after_transition_diagnostics ... ok
```

### Command

```powershell
cargo test --lib --no-default-features transition_failpoint_rolls_back_after_transition_diagnostics -- --nocapture
```

Verdict: `TransitionBoundaryRollbackInvariant` for the focused runner
transaction.

## 50. Real multicomponent continuation rollback route matrix

### Test

`frozen_reference_nasa_tp1906_chon_graphite_rollback_tests::tests::i5_tp1906_multicomponent_continuation_rollback_real_route_matrix`

### What it proves

This is the real-data A/B/C continuation check:

- Route A solves the three-point TP-1906/1907 CHON + graphite P,H range cleanly.
- Route B injects one failure after all phase transitions of the first
  accepted point, fails at the second point, then retries from a fresh
  transaction.
- Route C is an independent clean range.

The retry and fresh route are compared with Route A by target enthalpy,
temperature, and every physical component mole. Internal transition counts
and storage order are not used as equality criteria.

### Debug characterization

```text
NASA TP-1906/1907 real continuation rollback route matrix
route  status                     points  first_transitions  failure_point
A      CLEAN                          3                 25  -
B      FAILED@S1 -> RECOVERED        3                 25  1
C      FRESH CLEAN                   3                 25  -
verdict: TransactionalRollbackInvariant
```

The release test completed in `8.43 s`.

The failure is injected only after the first point has been accepted and its
25 phase transitions have completed. The failed second point is not published
as a partial `PhRangeSolution`; the subsequent retry publishes all three
points and agrees with both clean routes.

### Release command

```powershell
cargo test --release --lib --no-default-features i5_tp1906_multicomponent_continuation_rollback_real_route_matrix -- --ignored --nocapture
```

## 51. Test-only Ar/H2O/CO2 sublimation-gauge lifecycle

### Tests

`frozen_reference_ar_water_co2_sublimation_gauge::tests::production_runner_continues_the_full_gas_dry_ice_ice_lifecycle_both_directions`

`frozen_reference_ar_water_co2_sublimation_gauge::tests::three_phase_fixed_point_is_invariant_to_condensed_phase_declaration_order`

### What they prove

These tests use only a test-local relative-Gibbs gauge: IAPWS ice-Ih
sublimation pressure and the frozen NIST CO2 sublimation Antoine relation.
They do not introduce a runtime `CO2(s)` record and do not exercise P,H.

The ordinary prepared phase-control runner agrees with the independent scalar
oracle while cooling through:

```text
gas -> gas + CO2(s) -> gas + CO2(s) + H2O(ice)
```

and while heating back through the reverse lifecycle. The accepted paths have
exactly two transitions in each direction, finite TPD evidence, no empty
transition records, and accepted complementarity.

The three-phase fixed point is unchanged when the complete component layout is
permuted from `[gas, ice, dry_ice]` to `[gas, dry_ice, ice]`. Comparison is
performed by physical species identity, not by storage index or transition
history.

### Command

```powershell
cargo test --lib --no-default-features frozen_reference_ar_water_co2_sublimation_gauge -- --nocapture
```

Local result: `5 passed`. Release execution is still useful for final evidence;
the fixture itself is deterministic and has no network or database lookup.

## 52. IAPWS Ar/H2O exclusive competing liquid/ice candidates

### Test

`frozen_reference_iapws_exclusive_competing_candidates::tests::i5_iapws_water_has_two_unstable_condensed_candidates_and_order_invariant_fixed_point`

### What it proves

The test uses an independent test-only IAPWS gauge because the local NASA
gas/condensed records meet only at the singleton endpoint `273.15 K`, while
the frozen liquid and ice tables do not overlap. At `273.15 K` and `700 Pa`,
with `n(Ar)=0.1 mol` and `n(H2O)=1.0 mol`, the independent boundaries are:

```text
p_liquid = 611.212846 Pa
p_ice    = 611.153475 Pa
TPD_liquid = -9.158170e1 J/mol
TPD_ice    = -9.180231e1 J/mol
```

Both absent condensed phases are therefore genuine creation candidates and ice
is preferred by the independent gauge. The production phase-control runner
also reports both initial canonical candidates below `dg_create` and reaches
`gas + ice` for both declaration orders `[gas, liquid, ice]` and
`[gas, ice, liquid]`. Final physical moles, complementarity, and conservation
are compared after unpermuting component storage. Transition histories are
characterization only and are deliberately not required to be identical.

This test does not change production thresholds, candidate sorting,
thermochemistry, or reference-pressure policy.

### Command

```powershell
cargo test --lib --no-default-features i5_iapws_water_has_two_unstable_condensed_candidates_and_order_invariant_fixed_point -- --nocapture
```

Debug result: `passed`.

Release characterization:

```text
water candidate preflight: T=273.15 p_liquid=611.212845939855 p_ice=611.153475056703 delta_p=0.059370883152
Ar + H2O exclusive candidates: T=273.15 K P=700.0 Pa p_liquid=611.212846 Pa p_ice=611.153475 TPD_liquid=-9.158170e1 TPD_ice=-9.180231e1 final=gas+ice order_invariant=true
test ...::i5_iapws_water_has_two_unstable_condensed_candidates_and_order_invariant_fixed_point ... ok
test result: ok
```

Release command:

```powershell
cargo test --release --lib --no-default-features i5_iapws_water_has_two_unstable_condensed_candidates_and_order_invariant_fixed_point -- --nocapture
```

The final inactive-liquid TPD is also checked against the independent oracle
`RT ln(p_liquid / p_ice)` and is required to be positive. This demonstrates
that liquid is absent because it is thermodynamically losing, rather than only
because it was omitted from the active set.

### Metastable replacement and retention matrix

`frozen_reference_iapws_exclusive_competing_candidates::tests::i5_iapws_water_metastable_replacement_and_stable_ice_retention_matrix`

The six routes combine initial histories `gas`, `gas+liquid`, and `gas+ice`
with both declarations `[gas, liquid, ice]` and `[gas, ice, liquid]`. Every
route terminates at `gas+ice`, preserves the same physical mole state, and
passes complementarity. The verdict is:

```text
ExclusiveCompetitionHistoryInvariant
```

Debug characterization:

```text
history=gas        order=L->I final=gas+ice transitions=1
history=gas        order=I->L final=gas+ice transitions=1
history=gas+liquid order=L->I final=gas+ice transitions=2
history=gas+liquid order=I->L final=gas+ice transitions=2
history=gas+ice    order=L->I final=gas+ice transitions=0
history=gas+ice    order=I->L final=gas+ice transitions=0
verdict: ExclusiveCompetitionHistoryInvariant
```

The liquid amount is at the trace floor in every accepted state; the ice
amount and gas state agree across all six routes. Different transition counts
are history characterization, while the physical fixed point is a strict
invariant.

### Controlled Gibbs degeneracy matrix

`frozen_reference_iapws_exclusive_competing_candidates::tests::i5_controlled_condensed_gibbs_split_matrix_has_analytic_losing_tpd`

This is deliberately a synthetic numerical stress test, not frozen chemistry:
the only changed quantity is `G0_B - G0_A = delta_G`. For
`delta_G > 0`, phase A is the unique thermodynamic winner, both declaration
orders converge to the same physical state, and the inactive phase B has
`TPD_B = delta_G` within numerical tolerance. The exact `delta_G=0` row is
accepted as a non-unique representative and does not require a phase label.

Debug characterization:

```text
delta_G=1.0e0  transitions=1/1
delta_G=1.0e-1 transitions=1/1
delta_G=1.0e-2 transitions=1/1
delta_G=1.0e-3 transitions=1/1
delta_G=1.0e-4 transitions=1/1
delta_G=1.0e-5 transitions=1/1
delta_G=1.0e-6 transitions=1/1
delta_G=0.0e0  transitions=1/1
verdict: ResolvedCompetitionInvariant; exact delta_G=0 classified as non-unique representative
```

The companion scatter test repeats a well-resolved `delta_G=1e-6 J/mol` case
under both layouts. It records the observed minimum, maximum, and range of
the losing-phase TPD without treating that range as a global tolerance.

Release characterization:

```text
controlled Gibbs TPD scatter: delta_G=1.000e-6 samples=6 min=9.999985195464e-7 max=9.999985195464e-7 scatter=0.000e0
test ...::i5_controlled_condensed_gibbs_split_measures_observed_tpd_scatter ... ok
test result: ok
```

Release command:

```powershell
cargo test --release --lib --no-default-features i5_controlled_condensed_gibbs_split_measures_observed_tpd_scatter -- --nocapture
```

### Gibbs splitting resolution matrix

`frozen_reference_iapws_exclusive_competing_candidates::tests::i5_controlled_condensed_gibbs_split_resolution_matrix`

The matrix probes `delta_G` from `1e-6` down to `1e-15 J/mol` and exact zero.
It prints both the requested and actually representable `f64` perturbation,
both measured losing-phase TPD values, relative errors, winner labels, state
delta, transition counts, and a separate `FloatingPointCollapsedDegeneracy`
classification when the perturbation rounds away. It intentionally does not
use `abs(TPD-delta_G) < 1e-6` as a classification rule.

Debug characterization currently shows:

```text
smallest quantitatively resolved splitting = 1e-6 J/mol
first tolerance-limited splitting         = 1e-8 J/mol
floating-point collapse begins            = 1e-13 J/mol
exact degeneracy                           = finite termination, non-unique representative
```

The exact-zero row compares gas composition and total condensed amount, not
the arbitrary distribution between exactly degenerate phase labels. These
boundaries are fixture-specific and require release confirmation.

Release characterization confirmed the same matrix and summary:

```text
smallest_representable=1e-12
smallest_quantitative=1e-7
smallest_unique_winner=1e-7
first_tolerance_limited=1e-8
floating_point_collapse=1e-13
exact=ExactDegeneracyNonUnique
test result: ok
```

No production threshold or tie-breaker was changed.

## 53. Vanishing species inside one active phase

### Tests

`frozen_reference_extreme_dynamic_range::tests::analytic_single_phase_species_ratios_match_gibbs_weights`

`frozen_reference_extreme_dynamic_range::tests::i5_extreme_composition_dynamic_range_matrix`

`frozen_reference_extreme_dynamic_range::tests::i5_vanishing_species_resolution_sweep`

### What they prove

This is a synthetic analytical fixture, not external chemistry. Five species
share one abstract conserved element and one ideal gas phase. The exact
equilibrium fractions are constructed from Gibbs-derived weights, so the
tests independently verify mole fractions, log-ratios, positivity, and
conservation. Dynamic-range cases reach `1e-32`; the dedicated trace sweep
reaches `1e-40` and still reports the requested abundance rather than the
`1e-30` seed floor. The phase topology remains one active gas phase in every
case.

Current verdict:

```text
ExtremeCompositionDynamicRangeInvariant
VanishingSpeciesLogResolved
```

The seed, permutation, representability, and scaling matrices are separate
evidence and are intentionally not implied by the three baseline tests.

Release characterization: all three tests passed. The release run reproduced
the debug evidence, including the `epsilon=1e-40` row with
`actual_over_expected=1.000000e0`, finite positive mole amount, and balance
error `2.887e-15`. No species was pinned to the `1e-30` seed floor and the
single gas phase remained active throughout.

The subsequent seed and permutation extensions also pass. The first extensive
scaling probe at `N=1e-8 mol` is intentionally classified separately:
absolute balance error remains finite, but its relative composition error is
about `8.5e-5`; this is scale-limited characterization of the current numeric
contract, not evidence of a trace floor or a species-order defect.

The expanded debug run contains seven passing tests. It additionally proves:

- the returned state satisfies the analytic chemical-potential equality, with
  observed pairwise deviations below `1e-6 J/mol`;
- species order and initial seed do not change the physical composition;
- the Gibbs-weight fixture itself remains representable through ratio
  `1e-100`, so no floating-point collapse was observed in that matrix;
- extensive scaling is exact at `N=1 mol`, scale-limited at `N=1e-8 mol`,
  and the `N=1e8 mol` fresh solve is explicitly reported as
  `AllBackendsFailed`, rather than being counted as a false invariant.

The last point is an open numerical characterization item, not a hidden test
failure and not a reason to alter production tolerances without a separate
conditioning investigation.

## 54. NASA CEA RP-1311 Example 3 source preflight

### Test

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_rp1311_example3::tests::rp1311_example3_source_dataset_has_complete_semantic_universe`

### What it proves

The immutable source layer is complete before the production solver is used:

- three `H,P` rows exist at `100`, `10`, and `1 bar`;
- all 40 published non-trace species are retained in exact source order;
- 39 complete-calculation trace-only species are recorded categorically,
  without invented numerical mole fractions;
- the dataset id, metadata, rows, and frozen catalog entry agree.

### Debug command

```powershell
cargo test --lib --no-default-features rp1311_example3_source_dataset_has_complete_semantic_universe -- --nocapture
```

### Recorded debug output

```text
RP-1311 Example 3 preflight: rows=3 non_trace_species=40 trace_only_species=39 verdict=ExternalFixtureSourceComplete
test result: ok. 1 passed; 0 failed
```

The production H,P regression is intentionally not claimed by this test yet.
Feed reconstruction and enthalpy/reference-state validation are the next
stage.

### Feed reconstruction preflight

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_rp1311_example3::tests::rp1311_example3_feed_reconstruction_preserves_source_conventions`

The independent one-basis reconstruction passes:

```text
total_mass=18.000000 kg
O/F=17.000000
H_target=5.721084e6 J
elements={C,H,O,N,Ar: all positive}
```

This proves source feed bookkeeping only. It does not yet claim agreement
between local reactant enthalpies and the CEA H,P convention.

### Liquid fuel enthalpy preflight

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_rp1311_example3::tests::rp1311_example3_liquid_fuel_enthalpy_preflight_is_offline_and_explicit`

The two pure condensed records resolve without NIST fallback and cover
`298.15 K`:

```text
h_C7H8(L)             =  1.217868e4 J/mol
h_C8H18(L),n-octa     = -2.502669e5 J/mol
fuel contribution     = -1.261649e6 J
```

The CEA Air contribution remains explicitly deferred. CEA defines Air by an
elemental composition, so this preflight does not replace it with an assumed
dry-air thermochemical record.

### Production P,H comparison across the three CEA pressure rows

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_rp1311_example3::tests::i5_nasa_cea_rp1311_example3_production_ph_comparison`

Release/debug command:

```text
cargo test --release --lib --no-default-features i5_nasa_cea_rp1311_example3_production_ph_comparison -- --ignored --nocapture
```

The test invokes the canonical production `solve_resolved_ph` route for the
100, 10, and 1 bar frozen rows. The local gas catalog has no octane gas record,
so the initial computational seed uses an element-equivalent `CO + H2`
representation; the CEA source enthalpy remains the frozen P,H target. This
preserves the source C/H/O/N/Ar inventory without inventing a production
`Air` alias or silently enabling NIST fallback.

The first run exposed a fixture defect: it passed the current system pressure
as both `P` and `P0`. The fixture now uses fixed NASA gas standard pressure
`P0=100000 Pa` and asserts the expected ratios `100`, `10`, and `1`.

A second fixture defect was then found: liquid-fuel molecular masses expressed
in `g/mol` were used as `kg/mol`, reducing the fuel amount by `1000`. The
production fixture now resolves `C7H8(L)` and `C8H18(L),n-octa` locally,
converts `g/mol -> kg/mol` explicitly, and independently reconstructs the
physical `1 kg fuel + 17 kg Air = 18 kg` basis. The element-equivalent
`CO+H2` seed is checked against that local C/H/O/N/Ar inventory.

Corrected characterization:

```text
feed audit:
n_C7H8=4.341172768 mol  n_C8H18=5.252468660 mol
fuel_C=72.40795866 mol  fuel_H=129.2738180 mol
nominal_Air=586.9067991 mol  reconstructed_mass=18.00000000 kg
H_target=5.721084000e6 J

pressure_Pa | P/P0 | CEA_T_K | KiThe_T_K | delta_T_K | major_rel | minor_trace_max_log10 | residual | balance | backend
10000000    | 100  | 2418.660 | 2419.501 | +0.841 | 9.424e-2 | 9.622e-1 | 4.555e-14 | 2.274e-13 | Legacy(LM)
1000000     | 10   | 2390.593 | 2391.538 | +0.945 | 9.427e-2 | 9.735e-1 | 6.927e-14 | 2.274e-13 | Legacy(LM)
100000      | 1    | 2338.840 | 2339.773 | +0.933 | 9.537e-2 | 9.961e-1 | 5.603e-14 | 9.663e-12 | Legacy(LM)
```

What this proves:

- all three rows execute through the production P,H API;
- the local solve is finite, normalized, and element-conservative;
- the solver residual and balance contracts are satisfied;
- fixing `P0` restores the intended ideal-gas pressure convention;
- fixing the fuel `g/mol` versus `kg/mol` error restores the physical source
  inventory and CEA-scale equilibrium temperatures;
- all three temperatures agree with CEA to under `1 K` without a production
  algorithm, tolerance, seed-policy, or thermochemistry change;
- major species agree within about `9.6%`, while minor/trace evidence is
  retained as a separate logarithmic characterization metric;
- the two discarded outcomes are explicitly classified
  `ExternalFixtureStandardStatePressureError` and
  `ExternalFixtureFuelAmountUnitError`.

This is an external characterization result, not a claim that the local NASA
model reproduces the CEA pressure dependence. No production algorithm or
tolerance was changed to obtain it.

### Complete source-side reactant enthalpy characterization

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_rp1311_example3::tests::rp1311_example3_cea_air_reactant_enthalpy_characterization`

The test-side CEA Air decomposition is:

```text
N2  = 0.780840
O2  = 0.209476
Ar  = 0.009365
CO2 = 0.000319
```

It is atom-balanced against the published pseudo-reactant and sums to one
nominal mole. All four records resolve from local `NASA_gas` without NIST
fallback. The complete source-side enthalpy comparison is:

```text
Air molar mass       = 2.896541670e-2 kg/mol
nominal Air moles    = 5.869067991e2
H_air                = 6.982284751e6 J
H_fuel               = -1.261649098e6 J
H_local              = 5.720635653e6 J
H_CEA                = 5.721084000e6 J
h_local              = 3.178130918e5 J/kg
h_CEA                = 3.178380000e5 J/kg
delta                = -2.490816388e1 J/kg
relative             = 7.836748243e-5
verdict              = ReactantEnthalpyConventionAligned
```

This decomposition remains test-only and does not create a generic
production `Air` alias.

Command:

```powershell
cargo test --lib --no-default-features frozen_reference_extreme_dynamic_range -- --nocapture
```

Command:

```powershell
cargo test --lib --no-default-features i5_controlled_condensed_gibbs_split_resolution_matrix -- --nocapture
```

Command:

```powershell
cargo test --lib --no-default-features i5_controlled_condensed_gibbs_split_matrix_has_analytic_losing_tpd -- --nocapture
```
## 55. Frozen ElementInventory facade witnesses and independent K_eq validation

### Tests

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_h2_o2_hp_tests::tests::i5_nasa_cea_hp_element_inventory_facade_matches_frozen_reference`

`Thermodynamics::ChemEquilibrium::frozen_reference_argonne_stanjan_chon_tests::tests::i5_argonne_stanjan_chon_element_inventory_facade_matches_frozen_reference`

`Thermodynamics::ChemEquilibrium::cross_validation::real_boudouard_ph_tests::tests::p9_i1_i2_i4_boudouard_independent_reference_matches_fixed_topology_ph`

### Hypothesis / question

Does the public `EquilibriumCalculator` facade preserve independent frozen
NASA CEA and Argonne/STANJAN results when feeds are supplied as
`ElementInventory`, and does the real Boudouard reaction-quotient route agree
with the production result through an independent `K_eq` calculation?

### What it checks

- NASA CEA H2/O2: exact nine-species named `gas` phase, fixed-P,H
  `NestedTemperature` route, facade input provenance, external species envelope,
  and unchanged frozen/library snapshots;
- Argonne/STANJAN CHON: exact 15-species named `gas` phase, fixed-P,T
  production cascade, C/H/O/N atom totals `6 / 16 / 20 / 75.2`, species identity
  comparison, and unchanged snapshots;
- real Boudouard: independent reaction quotient versus evaluated `K_eq`, plus
  enthalpy and elemental-balance gates.

### Commands

    cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_cea_h2_o2_hp_tests::tests::i5_nasa_cea_hp_element_inventory_facade_matches_frozen_reference --no-default-features -- --ignored --nocapture
    cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_argonne_stanjan_chon_tests::tests::i5_argonne_stanjan_chon_element_inventory_facade_matches_frozen_reference --no-default-features -- --ignored --nocapture
    cargo test --release --lib Thermodynamics::ChemEquilibrium::cross_validation::real_boudouard_ph_tests::tests::p9_i1_i2_i4_boudouard_independent_reference_matches_fixed_topology_ph --no-default-features -- --nocapture

### Recorded result

The two new facade witnesses passed: CEA `|delta T|=1.991 K`, residual
`9.336e-9`, balance `3.815e-9`; Argonne residual `1.482e-11`, balance
`1.152e-11`. The existing Boudouard story remains the independent `K_eq`
source of truth; no artificial frozen benchmark was added.

## 56. Monolithic P,H source-property evaluation baseline

### Test name

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_monolithic_ph_property_evaluation_characterization`

`Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::facade_attaches_opt_in_instrumentation_to_a_ph_outcome`

### Command

```powershell
cargo test --release --lib `
  i5_tp1906_monolithic_ph_property_evaluation_characterization `
  --no-default-features -- --include-ignored --nocapture
```

### Hypothesis / question

How often does one source-backed monolithic `P,H` solve evaluate `G(T)`,
`H(T)`, and `Cp(T)`, and does the shared `evaluate_context` request property
families that are unused by the current residual or Jacobian caller?

### What it checks

- the real frozen NASA TP-1906/1907 CHON + graphite fixture reaches the
  accepted monolithic phase-control route;
- context, residual, Jacobian, and candidate-snapshot evaluation counts;
- component-level Gibbs, enthalpy, and heat-capacity calls;
- coefficient refreshes and source property calculation calls;
- counter consistency without imposing a machine-dependent performance limit.

### Baseline release output before context separation

```text
TP-1906 monolithic P,H property evaluation: species=18 context=176 residual=87 jacobian=36 snapshots=140 G_vectors=182 H_vectors=179 Cp_vectors=178 G_component_calls=3213 H_component_calls=3196 Cp_component_calls=3142 G_coefficient_refreshes=90 property_coefficient_selections=6338 property_calculations=6338
test Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_monolithic_ph_property_evaluation_characterization ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 2529 filtered out; finished in 0.02s
```

### Release output after context separation

```text
TP-1906 monolithic P,H property evaluation: species=18 context=176 residual_context=140 jacobian_context=36 residual=87 jacobian=36 snapshots=140 G_vectors=146 H_vectors=179 Cp_vectors=38 G_component_calls=2578 H_component_calls=3196 Cp_component_calls=670 G_coefficient_refreshes=90 property_coefficient_selections=3866 property_calculations=3866
test Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_monolithic_ph_property_evaluation_characterization ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 2529 filtered out; finished in 0.02s
```

### Result and decision

The baseline hypothesis was confirmed. After separating the contexts, the
accepted route and all thermochemistry contracts remain green while the same
release story changes `G_vectors` from `182` to `146`, `Cp_vectors` from `178`
to `38`, and source property calculations from `6338` to `3866`. The remaining
Gibbs and heat-capacity calls belong to preparation, stability, or Jacobian and
residual work respectively; no machine-dependent timing threshold is imposed.
The story now asserts that residual-side contexts are exactly the candidate
snapshot contexts and Jacobian-side contexts are exactly the analytic Jacobian
contexts.

## 57. ElementInventory semantic-to-numeric boundary

### Test name

`Thermodynamics::ChemEquilibrium::phase_bridge::phase_equilibrium_problem_tests::elemental_and_species_feeds_with_the_same_b_share_the_fixed_pt_solution`

### Command

```powershell
cargo test phase_equilibrium_problem_tests --lib --no-default-features
```

### Hypothesis / question

Can the bridge preserve one physical conserved vector `b` while keeping
semantic `ElementInventory`, phase-qualified physical composition, and the
independent solver-column representation distinct?

### What it checks

- an explicit molecular feed and an `ElementInventory` feed with the same
  physical H/O totals build the same component matrix and conserved totals;
- reaction-basis rank and reaction matrix remain identical;
- the report preserves different input provenance and seed sources;
- accepted fixed-P,T solutions agree component-by-component;
- the shared dense-moles-to-element-totals calculation does not alter the
  existing validation contract.

### Recorded result

```text
test Thermodynamics::ChemEquilibrium::phase_bridge::phase_equilibrium_problem_tests::elemental_and_species_feeds_with_the_same_b_share_the_fixed_pt_solution ... ok

test result: ok. 35 passed; 0 failed; 0 ignored; 0 measured; 2495 filtered out; finished in 0.19s
```

The consolidation decision is deliberate: `ElementInventory` remains the
semantic authority for closed-system `b`; dense phase composition remains a
physical seed value; `EquilibriumProblem` receives the numeric independent
column vector required by the solver. The public representation boundaries are
therefore retained, while duplicate totals arithmetic and duplicate seed
reconstruction are removed.

## 58. Backend seam preserves deterministic cascade policy

### Test names

`Thermodynamics::ChemEquilibrium::equilibrium_backend_adapter::tests::fake_backend_can_be_used_as_a_deterministic_test_double`

`Thermodynamics::ChemEquilibrium::equilibrium_backend_adapter::tests::fake_backend_can_force_a_repeatable_failure`

`Thermodynamics::ChemEquilibrium::equilibrium_log_moles::solver_cascade_story_tests::controlled_cascade_reports_each_started_backend_attempt`

`Thermodynamics::ChemEquilibrium::equilibrium_log_moles::solver_cascade_story_tests::first_backend_fails_recoverably_second_succeeds_and_is_accepted`

`Thermodynamics::ChemEquilibrium::equilibrium_log_moles::solver_cascade_story_tests::first_backend_returns_a_rejected_candidate_and_the_next_one_succeeds`

`Thermodynamics::ChemEquilibrium::equilibrium_log_moles::solver_cascade_story_tests::every_backend_fails_and_the_report_keeps_the_declared_order`

`Thermodynamics::ChemEquilibrium::equilibrium_log_moles::solver_budget_tests::every_started_backend_receives_the_original_validated_seed`

`Thermodynamics::ChemEquilibrium::equilibrium_log_moles::solver_cascade_story_tests::global_budget_exhaustion_stops_later_backends_deterministically`

`Thermodynamics::ChemEquilibrium::phase_bridge::phase_equilibrium_problem_tests::prepared_runner_honours_cancellation_before_starting_a_backend`

### Hypothesis / question

Can `EquilibriumNonlinearBackend` remain as the internal backend seam without
being a redundant abstraction, while preserving deterministic fallback,
failure evidence, retry classification, attempt ordering, progress, and
cancellation behavior across the prepared `P,T`, historical host, and
monolithic `P,H` routes?

### What it checks

- the uniform `PreparedNonlinearSystem` contract rejects invalid dimensions;
- a fake backend can deterministically succeed or produce a typed backend
  failure;
- a recoverable first attempt allows the next backend to succeed;
- rejected candidates and exhaustive failures retain ordered attempt evidence;
- progress reports expose both started and finished backend attempts;
- cancellation is returned as a typed error before a backend starts;
- all started backends receive the original validated seed through the shared
  cascade policy.

### Commands

```powershell
cargo test equilibrium_backend_adapter --lib --no-default-features
cargo test controlled_cascade_reports_each_started_backend_attempt --lib --no-default-features
cargo test first_backend_fails_recoverably_second_succeeds_and_is_accepted --lib --no-default-features
cargo test first_backend_returns_a_rejected_candidate_and_the_next_one_succeeds --lib --no-default-features
cargo test every_backend_fails_and_the_report_keeps_the_declared_order --lib --no-default-features
cargo test every_started_backend_receives_the_original_validated_seed --lib --no-default-features
cargo test global_budget_exhaustion_stops_later_backends_deterministically --lib --no-default-features
cargo test prepared_runner_honours_cancellation_before_starting_a_backend --lib --no-default-features
```

### Recorded result

The adapter module passed `4` tests. Each of the seven cascade, budget, seed,
and cancellation commands passed its focused test (`1 passed; 0 failed`). The
seam is therefore retained as the owner of backend dispatch and attempt
evidence. Cancellation remains cooperative at the boundaries before and after
a backend attempt; interrupting a solver internally is outside this contract
and was not changed.

## 59. P,H callback lifetime and nested template ownership

### Test names

`Thermodynamics::ChemEquilibrium::ph::equilibrium_ph_workflow::tests::legacy_request_accepts_a_borrowed_enthalpy_callback`

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_continues_from_previous_accepted_state`

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_nested_ph_target_range_story`

### Hypothesis / question

Can the P,H request preserve borrowed callback compatibility without making the
canonical resolved-data route depend on caller-owned mutable state, while the
nested target-range runner reuses only topology-stable preparation and keeps
target-specific scalar and lifecycle state independent?

### What it checks

- deprecated arbitrary-callback construction accepts a callback borrowing a
  local value, so the compatibility lifetime is not accidentally narrowed to
  `'static`;
- canonical resolved thermochemistry remains owned by clone-cheap `Arc`
  capabilities and the request owns its resolved-system snapshot;
- nested P,H ranges hand the accepted temperature and composition to the next
  target;
- fixed declared phases build one reusable inner P,T template, while every
  target retains an independent scalar bracket;
- rejected scalar trials and bounded active-set lifecycle state are not shared
  as continuation state.

### Commands

```powershell
cargo test legacy_request_accepts_a_borrowed_enthalpy_callback --lib --no-default-features
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_continues_from_previous_accepted_state `
  --no-default-features -- --ignored --nocapture
cargo test --release --lib `
  Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_nested_ph_target_range_story `
  --no-default-features -- --ignored --nocapture
```

### Recorded result

The focused borrowed-callback test passed (`1 passed; 0 failed`). The recorded
three-target continuation story passes in debug with one template build and
accepted-state hand-off; its release rerun remains explicitly pending in the
earlier Story Test 10. The large nested story has a recorded release pass with
nine targets, one template build, and `163` template reuses. The ownership
decision is therefore to keep the generic callback lifetime for the
deprecated/reference route and keep the nested continuation template as a
narrow mutable cache, without introducing a new shared P,H state object.

## 60. P,H resolved-data copying and lock ownership

### Test names

`Thermodynamics::ChemEquilibrium::ph::equilibrium_ph_thermochemistry::tests::thermochemistry_bundle_is_send_sync_for_shared_capabilities`

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_reactive_gas_ph_target_range_continues_from_previous_accepted_state`

`Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_large_element_limited_nested_ph_target_range_story`

### Hypothesis / question

Which P,H data are deep-copied, which are cheap `Arc` capability clones, and
which synchronization is required by the current thermochemistry ownership
contract? Can that be documented without changing snapshot isolation or
assuming parallel P,H range execution that the production routes do not use?

### What it checks

- `ResolvedThermochemistry` and its public capability function type satisfy
  `Send + Sync` at compile time;
- `ResolvedPhaseSystem::Clone` remains understood as a deep snapshot over
  `SubsData`, not as a cheap shared pointer operation;
- numeric G/H/Cp callbacks share immutable capability handles while private
  calculator/cache mutation remains behind `Arc<Mutex<...>>`;
- retained symbolic sources remain isolated from numeric evaluation state;
- the nested range stories continue to reuse only topology-stable preparation,
  with accepted composition and temperature hand-off between targets.

### Commands

```powershell
cargo test equilibrium_ph_thermochemistry --lib --no-default-features
cargo test equilibrium_ph_range --lib --no-default-features
```

The real-data continuation commands and their release evidence remain in
Story Tests 10, 13, and 19; they are not duplicated here.

### Recorded result

The thermochemistry suite passed `3` tests and the range/presentation suite
passed `6` tests. The audit confirms that replacing the resolved snapshot with
a blanket `Arc`, or removing the locks, would alter isolation and/or mutable
coefficient-cache semantics without a demonstrated production need. Runtime
clone allocation and lock-contention benchmarks remain a separate future
performance task, not an assumption of this story.

## 61. Common opt-in equilibrium instrumentation

### Test name

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_tp1906_monolithic_ph_property_evaluation_characterization`

### Command

```powershell
cargo test i5_tp1906_monolithic_ph_property_evaluation_characterization --lib --no-default-features -- --ignored --nocapture
cargo test facade_attaches_opt_in_instrumentation_to_a_ph_outcome --lib --no-default-features
cargo check --release --features equilibrium-instrumentation --lib
```

### Hypothesis / question

Can fine-grained thermochemistry evaluation counts live in the common
postprocessing/diagnostic layer instead of in P,H thermochemistry, while the
ordinary build remains disabled by default and release benchmarking has an
explicit opt-in feature? Are the measured counts unchanged by the move?

### What it checks

- the former `PhThermochemistryEvaluationCounters` and every `record_ph_*`
  caller are owned by `equilibrium_instrumentation`;
- one reusable `capture` call collects the monolithic P,H evaluation contract;
- the accepted TP-1906 solve and all characterization assertions remain the
  same (`176` contexts, `140` residual contexts, `36` Jacobian contexts,
  `146/179/38` G/H/Cp vectors, and `3866` property calculations);
- stage durations continue to come from the existing enabled
  `EquilibriumTimingReport`, while the instrumentation report supplies the
  outer capture duration and fine-grained counts;
- no cross-thread aggregation is claimed until the workflow explicitly passes
  a collector handle to worker closures.

### Recorded result

The focused characterization passed (`1 passed; 0 failed`) with the same
evaluation counts as Story Test 56. The common instrumentation unit tests
passed (`2 passed; 0 failed`). The facade attachment test passed in both the
default test build and the feature-enabled calculator suite (`22 passed; 0
failed`). Release library checking with
`equilibrium-instrumentation` passed. The default build keeps hooks compiled
out. The facade attachment test now proves that an explicitly enabled mode
publishes the snapshot on a typed P,H outcome and keeps stage timing in the
existing inner timing report. Explicit worker-thread propagation remains a
follow-up TODO item.

## 62. ElementInventory range preparation does not duplicate the bridge

### Test name

`Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::facade_solves_pt_then_ph_without_manual_resolved_request_assembly`

`Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::facade_accepts_element_inventory_for_pt_temperature_range`

### Command

```powershell
cargo test facade_solves_pt_then_ph_without_manual_resolved_request_assembly --lib --no-default-features -- --nocapture
cargo test facade_accepts_element_inventory_for_pt_temperature_range --lib --no-default-features -- --nocapture
```

### Hypothesis / question

Can elemental range constructors derive their first numerical composition from
the resolved element matrix without performing a full thermochemical bridge
build, while the first real solve still publishes the same physical inventory,
provenance, continuation, and prepared-formulation evidence?

### What it checks

- the `P,T` and `P,H` elemental constructors use matrix-only feasible-seed
  preparation;
- the canonical bridge is built exactly once per range solve, rather than once
  during construction and once again for the first point;
- the range still accepts all points transactionally and reuses continuation
  state/formulation;
- every accepted point still reports `ElementInventory` with unchanged `b`;
- the instrumentation counter measures bridge preparation, not numerical
  iterations, and remains opt-in at the facade boundary.

### Recorded result

Both focused tests passed (`1 passed; 0 failed` each). The captured output was:

```text
ElementInventory range optimization P,H bridge_builds=1 points=2 continuation_points=1 formulation_builds=1
ElementInventory range optimization P,T bridge_builds=1 points=3 formulation_builds=1 formulation_reuses=2
```

The pre-change call graph required two bridge builds for each elemental range
(`constructor seed derivation` plus `first real solve`). The optimized path
retains one canonical bridge build and moves only answer-independent element
matrix/seed construction to the lightweight domain helper. Existing warnings
are unrelated deprecated APIs and unused legacy fields; no numerical failure
was observed.

## 63. Explicit multi-start failure retains every typed backend trace

### Test name

`Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::error_contract_tests::multi_start_failure_preserves_each_typed_cascade_and_aggregates_work`

`Thermodynamics::ChemEquilibrium::equilibrium_presentation::tests::multi_start_seed_failure_flattens_every_preserved_backend_trace`

### Command

```powershell
cargo test multi_start_failure_preserves_each_typed_cascade_and_aggregates_work --lib --no-default-features -- --nocapture
cargo test multi_start_seed_failure_flattens_every_preserved_backend_trace --lib --no-default-features -- --nocapture
cargo check --lib --no-default-features
```

### Hypothesis / question

When every explicit initial log-mole seed is rejected, does the prepared
runner preserve the complete work of every seed and its backend cascade, while
keeping the established solver-budget semantics? Can diagnostic presentation
consume that nested evidence without reconstructing it from a summary string?

### What it checks

- `SolverCascadeBudget` is documented and retained as a per-seed inner-cascade
  budget for explicit multi-start recovery;
- the runner's all-seeds-failed path uses typed
  `MultiStartSeedRecoveryFailed` causes instead of `InvalidProblem` plus only
  the first error text;
- every seed retains its original `AllBackendsFailed` trace, seed ordering and
  typed classification;
- aggregate started-backend and nonlinear-iteration values come from the
  recursive `ReactionExtentError` API, not duplicated local counters;
- presentation and frozen diagnostic traversal expose both nested backend
  traces when no candidate exists.

### Recorded result

Both focused tests passed (`1 passed; 0 failed` each), and `cargo check`
passed. The error retains two ordered seed causes, reports
`started_backend_attempts() = 2`, remains classified as
`AllBackendsFailed`, and the presentation projection returns both backend
rows. No solver equations or backend order were changed.

## 64. Reaction-basis nullspace and validation invariants

### Test name

`Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::tests_reaction_basis::rectangular_rank_deficient_basis_has_complete_independent_nullspace`

`Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::tests_reaction_basis::reaction_basis_rejects_invalid_matrix_and_rank_tolerances`

`Thermodynamics::ChemEquilibrium::equilibrium_constant_tests::tests::reaction_basis_rejects_dependent_reaction_columns`

### Command

```powershell
cargo test tests_reaction_basis --lib --no-default-features -- --nocapture
cargo test equilibrium_constant_tests --lib --no-default-features -- --nocapture
```

### Hypothesis / question

Does the SVD-derived reaction basis always contain the complete independent
nullspace implied by `m - rank(A)`, including rectangular rank-deficient
matrices? Do conservation and normalization remain valid, and does the
independent K_eq basis reject a dependent explicit reaction matrix instead of
trusting only its declared element rank?

### What it checks

- finite positive rank tolerance and finite positive-dimensional element
  matrices are required;
- near-rank-deficient singular values follow the documented absolute
  threshold;
- thin-SVD nullspace vectors and Gram-Schmidt completion are mutually
  independent and conserve every element;
- no post-hoc coefficient truncation changes the conservation equation;
- `ValidatedReactionBasis` retains species ordering and deterministic
  normalization while rejecting dependent reaction columns.

### Recorded result

All six direct reaction-basis tests and all twelve equilibrium-constant
contract tests passed. The previously failing rectangular rank-deficient
case now returns three independent nullspace columns with conservation below
the test tolerance. The implementation keeps `ReactionBasis` and
`ValidatedReactionBasis` as separate contracts; only their verified numerical
invariants were strengthened.

## 65. TP-1906 P,H extensive scaling envelope

### Test name

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together`

### Command

```powershell
cargo test tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together --lib --no-default-features -- --nocapture
cargo test elemental_calculator_ph_is_extensive_under_inventory_and_enthalpy_scaling --lib --no-default-features -- --nocapture
```

### Hypothesis / question

If the physical inventory is multiplied by `s` and the total enthalpy target
is multiplied by the same `s`, does the canonical nested P,H path preserve
the intensive temperature and phase topology while scaling every published
component amount by `s`? Separately, does the public elemental calculator
facade preserve the same law at its existing factor-10 contract?

### What it checks

- the real frozen TP-1906/1907 CHON + graphite fixture is solved at the
  baseline target and at `s = 1e-3, 1e-1, 1e1, 1e3, 1e4, 1e6, 1e8` through
  the canonical nested P,H request with `s*n_initial` and `s*H`;
- the public elemental calculator route independently checks `s = 10` with
  `s*b` and `s*H`, using the ordinary facade and no test-only solver;
- accepted temperature is intensive, graphite active/inactive topology is
  unchanged, enthalpy/residual/balance acceptance remains bounded, and
  component moles divided by `s` agree with the baseline;
- the check is a metamorphic scaling contract, not a new equilibrium oracle;
  the frozen chemical source and target remain unchanged apart from the
  declared extensive factor.

### Recorded result

The frozen nested P,H story passed for all seven factors (`1 passed; 0
failed`) in `87.17 s` in the debug profile. The independent public elemental
calculator story passed at factor `10` (`1 passed; 0 failed`). A TP-1906 elemental probe
at `1e-3` failed the intensive-temperature assertion with a reproducible
`1.19e-2 K` drift, while `1e-1` showed `1.49e-5 K`; the lower elemental
boundary is therefore open numerical debt, not a relaxed assertion. The
frozen nested contract is verified through `1e8`; behavior beyond `1e8` is
uncharacterized.

## 66. Production public surface and phase-interface compatibility

### Test name

`Thermodynamics::ChemEquilibrium::equilibrium_public_api_tests::tests::*`

`Thermodynamics::User_PhaseOrSolution_tests::tests::test_one_and_many_phase_facades_implement_the_same_narrow_roles`

### Command

```powershell
cargo test equilibrium_public_api_tests --lib --no-default-features -- --nocapture
cargo test test_one_and_many_phase_facades_implement_the_same_narrow_roles --lib --no-default-features -- --nocapture
```

### Hypothesis / question

Can the supported production workflow be consumed through the curated
`ChemEquilibrium::prelude` without naming implementation modules, while the
narrow phase interfaces remain usable by both one-phase and many-phase
facades? Is there evidence strong enough to remove or reorganize the public
surface now?

### What it checks

- the prelude names and exercises the application facade, elemental inventory,
  numerical seed, range seed constructors, solver policy, immutable reports,
  stability evidence, and transactional input validation;
- the retained legacy solver policy is explicit and the deprecated namespace
  is not confused with the production entry point;
- one-phase and many-phase facades implement the shared narrow layout,
  property, preparation, symbolic, and equilibrium-assembly contracts;
- no public removal or module reorganization is inferred without downstream
  consumer evidence.

### Recorded result

The public facade contract passed `6 passed; 0 failed`, and the phase-interface
compatibility witness passed `1 passed; 0 failed`. The current public surface
is therefore retained; removal of legacy exports remains a breaking-release
decision rather than an internal cleanup.

## 67. Public elemental P,H scale envelope and lower-bound characterization

### Story name

The public elemental P,H facade must expose its scale behavior honestly: joint
inventory and enthalpy scaling is an extensive transformation, but a passing
solve outside the measured envelope must not be promoted into a supported
equivalence contract.

### Test name

`Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::elemental_calculator_ph_scale_envelope_characterization`

`Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::elemental_calculator_ph_recovers_a_large_inventory_point_transactionally`

### Command

```powershell
cargo test elemental_calculator_ph_scale_envelope_characterization --lib --no-default-features -- --include-ignored --nocapture
cargo test elemental_calculator_ph_recovers_a_large_inventory_point_transactionally --lib --no-default-features -- --nocapture
```

### Hypothesis / question

Is the lower elemental P,H boundary caused by an invalid extensive
transformation, or by numerical acceptance and inner-basin conditioning? Does
the public point facade remain usable through the already measured upper
window, without claiming anything beyond `1e8`?

### What it checks

- the ordinary `EquilibriumCalculator` facade is used for every scale;
- both elemental inventory and total enthalpy are multiplied by the same `s`;
- accepted cases compare intensive temperature, scaled component moles,
  residual, and relative elemental balance against the `s=1` result;
- failed cases retain their typed error instead of being hidden by a relaxed
  assertion;
- the measurement includes lower probes `1e-6..1e-1`, the existing pair `1`
  and `10`, and upper probes through `1e8`.
- the dedicated point transaction starts with the same physical facade solve,
  then recovers a classified large-inventory numerical failure in exactly
  normalized coordinates; it asserts physical publication provenance,
  intensive temperature, and component amounts divided by the scale.

### Recorded result

```text
baseline T=2499.999998428 K
1e-6   error: enthalpy_non_monotone
1e-5   ok, delta_T=4.277201e-4 K
1e-4   ok, delta_T=4.292847e-4 K
1e-3   ok, delta_T=4.292847e-4 K
1e-2   error: enthalpy_acceptance
1e-1   ok, delta_T=4.102910e-4 K
1      ok, delta_T=0
1e1    ok, delta_T=5.343509e-7 K
1e2    ok, delta_T=1.858272e-6 K
1e3..1e8 ok, delta_T=4.046777e-5 K after typed normalized recovery
```

The dedicated `1e4` transaction story passed (`1 passed; 0 failed`). It
retained normalized-recovery evidence at the public physical boundary and
matched the unit-scale temperature within `1e-4 K` while every component
amount scaled by `1e4`.

### Conclusion / next action

The extensive transformation itself is valid where a solution is accepted,
but the current strict lower contract is not established below `s=1`: small
scales either show a reproducible temperature drift or fail typed P,H
acceptance. Direct facade P,H points now share the typed, failure-triggered
normalization transaction with the canonical point workflow, restoring the
measured upper `1e3..1e8` window without changing the default solver order or
acceptance thresholds. This remains evidence rather than a contract beyond
`1e8`. Improving the lower boundary needs a separately measured
preconditioning policy because those physical solves are accepted before
failure recovery can be authorized.

## 68. TP-1907 normalization regressions after exact reaction-basis restoration

### Story name

Restoring the exact conservation-preserving reaction basis must not regress the
independent extensive-normalization evidence or turn a fresh-seed basin
diagnostic into a false physical failure.

### Test names

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_tp1907_extensive_normalization_scale_metamorphism`

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_tp1907_large_inventory_extensive_normalization_matrix`

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::i5_tp1907_large_inventory_seed_basin_matrix_diagnostic`

### Command

```powershell
cargo test i5_tp1907_extensive_normalization_scale_metamorphism --lib --no-default-features -- --include-ignored --nocapture
cargo test i5_tp1907_large_inventory_extensive_normalization_matrix --lib --no-default-features -- --include-ignored --nocapture
cargo test i5_tp1907_large_inventory_seed_basin_matrix_diagnostic --lib --no-default-features -- --include-ignored --nocapture
```

### Hypothesis / question

Did the reaction-basis hardening break the exact normalized formulation, or is
the observed failure limited to an under-seeded fresh nonlinear basin? Does
the normalized route still reconstruct the same physical TP-1907 solution?

### What it checks

- exact `A^T*N = 0` reaction-basis conservation is retained without post-hoc
  coefficient truncation;
- the normalized inventory remains one mole and is derived from the physical
  request, not from an accepted oracle;
- 680/700 K preserve gas plus graphite, while 720 K uses the correct gas-only
  projection and positive graphite TPD;
- physical reconstruction preserves component amounts, gas fractions,
  residual, balance, and topology;
- the seed-basin diagnostic distinguishes fresh-seed failure from a valid
  transformed accepted seed and records the original backend trace;
- the metamorphic reference solve uses the existing canonical row-scaling
  contract, so arbitrary reaction-column norms do not become a numerical
  artifact of the comparison.

### Recorded result

All three tests passed (`1 passed; 0 failed` each). The scale metamorphism
accepted all fifteen rows for factors `1e-4, 1e-2, 1, 1e2, 1e4` across 680,
700, and 720 K. The large-inventory matrix accepted its three `1e4` rows with
the largest displayed reconstructed mole error `3.450e-9` and normalized
720 K gas-only topology. The basin matrix classified 680/700 K as
`FreshSeedBasinLimited` and 720 K as `StrictScaleInvariant`; the transformed
accepted seed remained valid in every applicable comparison.

The diagnostic fixed-active fixture uses `1e-20` only as a seed-only basin
hint. The production recovery rule is independent of that fixture: when an
absolute physical trace floor collapses after normalization, it uses a
relative normalized seed and retains the transformed absolute floor as its
minimum bound. No frozen value, equation, or acceptance tolerance was
changed.

The full ChemEquilibrium run also exposed and closed a test-only conditioning
gap: the exact-basis metamorphic reference cascade failed when its deliberately
rescaled reaction rows were left unscaled. Enabling `scaling_flag` only for
that reference solve restores deterministic comparison of the same physical
root; production defaults and backend order remain unchanged.

### Conclusion

The failure was a fresh-seed conditioning issue exposed by the exact reaction
basis, not a broken extensive-normalization identity. The exact basis remains
in production, and the three previously failing regressions are green. The
remaining direct elemental P,H point limitation is recorded separately in
Story 67 and must not be hidden by changing solver tolerances.

## 69. Production facade performance baseline: P,T, P,H, and ElementInventory

### Story name

The canonical `EquilibriumCalculator` facade must have a repeatable performance
baseline across the main point and range workflows before optimization work
begins. Timing and evaluation counts must come from the shared internal
`postprocessing_and_logging` instrumentation, not benchmark-specific wrappers
around solver internals.

### Benchmark target and command

`benches/equilibrium_baseline.rs`

```powershell
$env:CARGO_PROFILE_BENCH_LTO='false'
$env:CARGO_PROFILE_BENCH_CODEGEN_UNITS='16'
$env:CARGO_PROFILE_RELEASE_LTO='false'
$env:CARGO_PROFILE_RELEASE_CODEGEN_UNITS='16'
cargo bench --bench equilibrium_baseline --features equilibrium-instrumentation -- --noplot
```

The recorded run used the optimized `bench` profile with LTO disabled only for
practical build time, `10` fixed-repeat instrumentation observations per row,
and `10` Criterion samples per scenario with a `1 s` warm-up and `2 s`
measurement window. The repository's ordinary default build remains
uninstrumented; these settings are part of the measurement identity.

### Hypothesis / question

Which stages and evaluation counters dominate the production facade for the
main fixed-`P,T`, fixed-`P,H`, and continuation-range tasks? What is the
baseline cost of replacing explicit molecular initial amounts with a physical
`ElementInventory` while keeping the same declared H2/O2/H2O universe and
NASA-gas repository?

### What it checks

- all cases enter through `EquilibriumCalculator`, `production_cascade()`,
  offline local NASA data, and the public facade outcome types;
- molecular and elemental inventory inputs are measured separately for `P,T`
  and `P,H` points;
- elemental `P,T` and `P,H` ranges measure continuation and accepted-point
  work, with range stage values reported per accepted point;
- the table uses `EquilibriumTimingReport` and
  `EquilibriumInstrumentationReport` from the shared diagnostic layer;
- Criterion supplies wall-time distributions while the internal table retains
  stage/counter evidence; neither timing nor counters are used as correctness
  assertions.

### Status

This first compact baseline is superseded. Its wide table mixed facade wall
time, attached inner reports, inclusive P,H fields, and counters from an older
instrumentation contract. The current benchmark and its authoritative
macro/stage tables are recorded in Story 72; no timing values from this story
should be used for optimization comparisons.

## 70. Expanded component-count performance baseline and P,H attribution

### Story name

The performance baseline must include real high-dimensional frozen universes,
not only the three-species smoke system. The measurement must expose whether
the observed P,H wall time belongs to the reported nonlinear stages or to
facade preparation that is currently outside the timing report.

### Benchmark target and command

`benches/equilibrium_baseline.rs`

```powershell
$env:CARGO_PROFILE_BENCH_LTO='false'
$env:CARGO_PROFILE_BENCH_CODEGEN_UNITS='16'
$env:CARGO_PROFILE_RELEASE_LTO='false'
$env:CARGO_PROFILE_RELEASE_CODEGEN_UNITS='16'
cargo bench --bench equilibrium_baseline --features equilibrium-instrumentation -- --noplot
```

### Systems and question

- `small3`: H2/O2/H2O control, including molecular and ElementInventory P,T/P,H
  points and the two continuation ranges;
- `stanjan15`: the exact 15-species Argonne/STANJAN CHON gas universe;
- `tp1907_18`: the exact 17-species NASA gas universe plus solid `C(gr)`;
- `rp1311_40`: the exact 40-species NASA CEA RP-1311 gas universe.

For the three large systems the baseline measures the production
`ElementInventory` P,T and P,H facade routes. TP-1907 P,H uses the frozen
nested-temperature route and its source-backed target, because the monolithic
route is intentionally rejected across the graphite coefficient boundary.
The 40-component P,H case uses nested temperature as the stable accepted
route for this baseline. These route choices are printed through accepted
backend/attempt/iteration evidence and are not silently mixed.

### Hypothesis / question

Does the cost grow materially with component count, and which measured stage
should be optimized first? In particular, can the previous 22--25 ms small
P,H facade interval be attributed to inner equilibrium work, or is it setup,
symbolic preparation, recovery, or result assembly outside the current timing
report?

### Status

The intermediate wide matrix, legacy-control matrix, and focused P,H rerun
are removed. They mixed facade, workflow, and inclusive child timings from
different report contracts, so their values are not valid baseline evidence.
The current benchmark command and the authoritative macro and detailed stage
views are recorded in Story 72.

One compatibility finding is retained: a fixed legacy LM policy did not
complete the TP-1907 nested P,H fixture during normalized basin discovery,
while the production cascade accepted it with legacy TR. This is a
convergence guard, not a successful timing row or an apples-to-apples
performance comparison.

## 71. Optional legacy nonlinear-solver diagnostic table

### Name

`legacy_solver_diagnostics_are_scoped_to_calculator_capture` plus the
solver-diagnostic table printed by `equilibrium_baseline`.

### Command

```text
cargo test --lib Thermodynamics::ChemEquilibrium::equilibrium_backend_adapter::tests::legacy_solver_diagnostics_are_scoped_to_calculator_capture --no-default-features -- --exact --nocapture

$env:CARGO_PROFILE_BENCH_LTO='false'
$env:CARGO_PROFILE_BENCH_CODEGEN_UNITS='16'
$env:CARGO_PROFILE_RELEASE_LTO='false'
$env:CARGO_PROFILE_RELEASE_CODEGEN_UNITS='16'
cargo bench --bench equilibrium_baseline --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" -- __solver_diagnostics__ --noplot
```

### Question and hypothesis

Do the retained legacy LM, NR, and TR solvers expose enough internal timing to
decide later whether residual/Jacobian or linear-system work is worth
parallelizing, without adding runtime cost to ordinary calculations? The
hypothesis is that a calculator-level instrumentation capture can own one
short-lived per-attempt scope, while the solver callbacks report calls and
durations upward through the existing `SolverAttemptMetrics` contract.

### What is checked

- Legacy LM, NR, and TR residual callbacks are counted and timed separately.
- Legacy analytic Jacobian callbacks are counted and timed separately.
- Legacy dense linear-system solves are counted and timed separately.
- The adapter reports the snapshot for every completed instrumented legacy
  attempt; this benchmark table selects the accepted backend. P,H nested rows
  aggregate accepted inner-attempt diagnostics over their temperature trials.
- The scope is enabled only while the calculator's existing `capture` is
  active. A direct solver call remains uninstrumented and returns no fabricated
  zero-valued metrics.
- The benchmark emits a separate `tabled` solver table alongside the
  hierarchical facade/stage report. It contains backend wall clock, iterations,
  residual/Jacobian/
  linear calls, absolute microseconds, and each measured time as a percentage
  of backend wall time.
- RST rows retain their available callback/overhead metrics but show `n/a` for
  linear-system time and percentage because the RST adapter does not expose
  that internal operation.

### Current report

The earlier compact solver table is removed because it was a historical
snapshot of the diagnostic output. Current residual, Jacobian, linear-system,
callback, backend wall-time, and failed-attempt evidence is emitted by the
separate backend-attempt and accepted-solver tables in Story 72.

### Result and decision

The hypothesis is supported. The legacy instruments are opt-in, scoped to the
calculator capture, and the payload reaches the accepted backend metrics
without changing solver equations or order. The first measurement also shows
why parallelism must not be assumed: in the tiny legacy P,T cases callback and
linear timings are below useful resolution and solver overhead dominates,
while the larger nested TP-1907 P,H case has measurable Jacobian and residual
work. Do not parallelize yet. First collect the same table at larger stable
component counts and compare instrumentation overhead; then evaluate only a
measured solver-level hypothesis.

### Limitations

Failed legacy attempts now publish a partial residual/Jacobian/linear-system
snapshot together with their backend-call wall time in the full attempt table,
so their fallback cost is visible at callback level. `wall-only` is retained
only for backends that do not expose the corresponding internal split, such as
RST linear algebra.
RST linear algebra remains `n/a` until its backend exposes a corresponding
hook.

## 72. Hierarchical performance-budget report without hidden percentages

### Story name

The performance report must explain every large timing budget at the level
where an optimization can act. A single wide table is not sufficient when it
mixes complete-call partitions with inclusive parent/child timings.

### Benchmark command

```powershell
$env:CARGO_PROFILE_BENCH_LTO='false'
$env:CARGO_PROFILE_BENCH_CODEGEN_UNITS='16'
$env:CARGO_PROFILE_RELEASE_LTO='false'
$env:CARGO_PROFILE_RELEASE_CODEGEN_UNITS='16'
cargo bench --bench equilibrium_baseline --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" -- --noplot
```

`cargo bench` uses Cargo's optimized benchmark profile. The fixed-repeat
pretty tables are emitted before Criterion's statistical benchmark output and
are intended for attribution, not as assertions about a universal machine
time. Without `equilibrium-benchmark-prettyprint`, the benchmark now prints
only the compact macro columns; the obsolete wide diagnostic row is no longer
part of either benchmark interface.

### Question and hypothesis

When a macro row reports tens of percent in the workflow, can that budget be
followed into facade, P,H orchestration, fixed-P,T stages, the nonlinear
envelope, and every backend attempt without counting a parent and its child as
independent work? The hypothesis is yes, provided each table declares its own
scope and the nonlinear remainder is shown explicitly.

### Report contract

The pretty report now emits these independent views:

1. **Macro budget:** complete facade wall split only into mutually exclusive
   `facade outside workflow` and `workflow`; these two columns sum to 100%.
2. **Facade detail:** resolution, request preparation, route selection,
   normalization, phase control, thermochemistry, and result assembly as
   facade-owned observations.
3. **P,H workflow detail:** scalar/enthalpy orchestration, trial totals and
   preparation, bracket control, monolithic preparation/execution/publication,
   with parent/child overlap explicitly labelled.
4. **Inner fixed-P,T detail:** thermochemistry, numeric closure, normalization,
   recovery, symbolic and equation preparation, numerical problem preparation,
   active-set formulation build, RST preparation, projection build, nonlinear
   envelope, phase control, validation, postprocessing, repository lookup, and
   result assembly. It also includes `outside nonlinear envelope` as an
   arithmetic complement. Monolithic P,H rows do not fabricate an inner P,T
   table; their payload is shown in the P,H table and backend table.
5. **Nonlinear envelope:** parent scope/wall, all attempt wall including
   wall-only failures, known active-set/RST/projection preparation, and the
   remaining runner/orchestration budget. It also prints backend coverage,
   total attempt count, `no metrics` count, and `unknown wall` count. This
   distinguishes a real orchestration gap from a backend whose internal split
   was not published.
6. **Backend attempt detail:** one row per backend attempt, including accepted,
   failed, rejected-candidate, skipped, fallback, and normalization-recovery
   scopes. It reports outcome, samples, wall, iterations, residual/Jacobian/
   linear-solve/solver-overhead times and percentages relative to that
   attempt's own wall.
7. **Accepted solver callback summary:** retained compact view for comparing
   callback cost across accepted legacy and RST attempts.

All tables use milliseconds for wall time and percentages rounded to one
decimal place. Detailed stage rows are inclusive diagnostics and are not
supposed to sum to 100%. A failed attempt without published counters is shown
as `wall-only`: its real elapsed call time is retained, while iterations and
callback splits remain `n/a`. `unknown wall` is a separate count and is zero
for the release capture; skipped backends are not started and therefore do not
consume wall time. The remaining debt is to retain partial legacy callback
timing when an interrupted legacy attempt never publishes a metrics snapshot.

### Decision use

This report is now the authoritative format for performance work. For every
`ElementInventory`, molecular, P,T, nested P,H, monolithic P,H, legacy, and
RST row, first use the macro partition to locate the owning workflow, then use
the matching lower-level table and its denominator. Do not sum the inclusive
facade/P,H/inner rows, and do not treat the runner/orchestration remainder as
a residual or Jacobian cost without an attempt-level measurement. A change is
performance evidence only when the same case, route, backend order, profile,
and correctness story remain fixed.

### Runtime verification of the rendered hierarchy

The optimized run on 2026-09-12 rendered all seven tables and passed the
benchmark's per-scenario coverage assertions for all 23 scenarios. The fixed-
repeat dashboard itself completed; Criterion later returned non-zero because
its statistical comparison flagged a `ph_element_inventory_range` regression,
which is a performance signal rather than a dashboard or correctness failure.
The assertions require every scenario to have macro, facade, nonlinear, backend
attempt, and accepted-solver rows; P,T scenarios must have inner detail, while
P,H scenarios must expose either nested trials plus inner detail or monolithic
execution detail. Criterion completed normally in this capture. Representative
macro rows from the latest capture were:

| case | wall | facade outside workflow | workflow | partition |
| --- | ---: | ---: | ---: | ---: |
| `stanjan15_ph_element_inventory_point` | 571.58 ms | 1.01 ms / 0.2% | 570.57 ms / 99.8% | 100.0% |
| `tp1907_18_ph_element_inventory_point` | 146.75 ms | 2.18 ms / 1.5% | 144.58 ms / 98.5% | 100.0% |
| `rp1311_40_pt_element_inventory_point` | 255.15 ms | 1.18 ms / 0.5% | 253.97 ms / 99.5% | 100.0% |
| `rp1311_40_ph_element_inventory_point` | 381.71 ms | 16.72 ms / 4.4% | 364.99 ms / 95.6% | 100.0% |

This is the macro table only. It is intentionally the one additive view: the
two partition columns account for the complete facade wall time. The detailed
facade, P,H workflow, inner fixed-P,T, nonlinear-envelope, backend-attempt,
and accepted-solver tables are separate rendered views with their own parent
scope and denominator; their percentages must not be added to this table or
to one another.

### Current detailed dashboard capture

The following rows are copied from the same ten-repeat optimized `bench`
profile capture on 2026-09-12. They are deliberately limited to the large
production cases and the small `ElementInventory` control, so the evidence is
readable while the command above still renders all 23 scenarios. Every table
below has an explicit scope. The rows are not a replacement for the complete
machine output; they are the retained optimization dashboard for the cases
that own the large budgets.

#### Facade-owned stages

| case | resolution | request preparation | route selection | normalization | phase control | thermochemistry | result assembly |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `ph_element_inventory_point` | 0.05 ms / 0.3% | 0.03 ms / 0.1% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.02 ms / 0.1% | 0.00 ms / 0.0% |
| `stanjan15_pt_element_inventory_point` | 0.21 ms / 1.2% | 0.26 ms / 1.5% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% |
| `stanjan15_ph_element_inventory_point` | 0.26 ms / 0.0% | 0.13 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.11 ms / 0.0% | 0.00 ms / 0.0% |
| `tp1907_18_pt_element_inventory_point` | 0.30 ms / 0.3% | 0.39 ms / 0.4% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% |
| `tp1907_18_ph_element_inventory_point` | 0.36 ms / 0.2% | 0.20 ms / 0.1% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.16 ms / 0.1% | 0.00 ms / 0.0% |
| `rp1311_40_pt_element_inventory_point` | 0.60 ms / 0.2% | 0.57 ms / 0.2% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% |
| `rp1311_40_ph_element_inventory_point` | 0.58 ms / 0.2% | 0.28 ms / 0.1% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.00 ms / 0.0% | 0.25 ms / 0.1% | 0.00 ms / 0.0% |

#### P,H workflow stages

These are inclusive observations. For nested P,H, `trial total` contains
`trial inner execution`; for monolithic P,H, `monolithic execution` contains
backend execution. The rows therefore explain the workflow without pretending
that all percentages form a second additive partition.

| case | stage | time | % of workflow | optimization meaning |
| --- | --- | ---: | ---: | --- |
| `ph_element_inventory_point` | monolithic RST preparation | 19.35 ms | 93.2% | RST payload/preparation boundary |
| `ph_element_inventory_point` | monolithic backend execution | 1.02 ms | 4.9% | accepted backend callbacks |
| `ph_element_inventory_point` | monolithic validation | 0.03 ms | 0.2% | acceptance checks |
| `ph_element_inventory_point` | monolithic execution | 20.69 ms | 99.6% | inclusive parent, not additive |
| `stanjan15_ph_element_inventory_point` | monolithic preparation | 0.36 ms | 0.1% | outer preparation |
| `stanjan15_ph_element_inventory_point` | monolithic RST preparation | 554.19 ms | 97.1% | primary performance target |
| `stanjan15_ph_element_inventory_point` | monolithic backend execution | 6.22 ms | 1.1% | accepted backend callbacks |
| `stanjan15_ph_element_inventory_point` | monolithic validation | 0.14 ms | 0.0% | not performance-critical |
| `stanjan15_ph_element_inventory_point` | monolithic execution | 570.18 ms | 99.9% | inclusive parent, not additive |
| `tp1907_18_ph_element_inventory_point` | P,H scalar orchestration | 144.51 ms | 100.0% | inclusive outer parent |
| `tp1907_18_ph_element_inventory_point` | trial total | 144.37 ms | 99.9% | contains trial children |
| `tp1907_18_ph_element_inventory_point` | trial preparation | 20.33 ms | 14.1% | template/prepared-state construction |
| `tp1907_18_ph_element_inventory_point` | trial inner execution | 123.92 ms | 85.7% | fixed-P,T work inside trials |
| `tp1907_18_ph_element_inventory_point` | bracket orchestration | 0.14 ms | 0.1% | not a target |
| `rp1311_40_ph_element_inventory_point` | P,H scalar orchestration | 364.88 ms | 100.0% | inclusive outer parent |
| `rp1311_40_ph_element_inventory_point` | trial total | 364.45 ms | 99.9% | contains trial children |
| `rp1311_40_ph_element_inventory_point` | trial preparation | 230.85 ms | 63.3% | primary preparation target |
| `rp1311_40_ph_element_inventory_point` | trial inner execution | 133.09 ms | 36.5% | fixed-P,T work inside trials |
| `rp1311_40_ph_element_inventory_point` | bracket orchestration | 0.43 ms | 0.1% | not a target |

#### Inner stages and nonlinear ownership

`outside nonlinear envelope` is an exclusive arithmetic complement. Other
inner rows are inclusive measurements and are intentionally read as ownership
evidence, not summed percentages.

| case | inner scope | outside envelope | RST preparation | numeric closure | nonlinear solve | validation |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `stanjan15_pt_element_inventory_point` | fixed P,T | 0.23 ms / 1.4% | 12.62 ms / 76.3% | 0.04 ms / 0.3% | 16.32 ms / 98.6% | 0.01 ms / 0.0% |
| `tp1907_18_pt_element_inventory_point` | fixed P,T | 77.21 ms / 78.2% | 20.24 ms / 20.5% | 0.06 ms / 0.1% | 21.55 ms / 21.8% | 0.01 ms / 0.0% |
| `tp1907_18_ph_element_inventory_point` | fixed P,T trial | 0.94 ms / 0.8% | 0.00 ms / 0.0% | 0.70 ms / 0.6% | 122.81 ms / 99.2% | 0.16 ms / 0.1% |
| `rp1311_40_pt_element_inventory_point` | fixed P,T | 0.58 ms / 0.2% | 233.85 ms / 92.1% | 0.10 ms / 0.0% | 253.39 ms / 99.8% | 0.02 ms / 0.0% |
| `rp1311_40_ph_element_inventory_point` | fixed P,T trial | 3.88 ms / 2.9% | 0.00 ms / 0.0% | 2.94 ms / 2.2% | 128.64 ms / 97.1% | 0.51 ms / 0.4% |

#### Nonlinear envelope and all-attempt accounting

This is the table that answers what lies inside the large workflow numbers.
`all attempt wall` includes accepted, rejected, failed, fallback, and recovery
attempts. `known preparation` is the measured active-set/RST/projection part;
`runner remainder` is the still-unattributed orchestration inside the parent,
not silently assigned to the solver. `no metrics` is different from
`unknown wall`: the current capture has elapsed wall time for every row.

| case | parent scope | parent wall | envelope | attempts | all attempt wall | known preparation | runner remainder | no metrics | unknown wall |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `stanjan15_ph_element_inventory_point` | monolithic backend execution | 6.22 ms | 6.22 ms / 100.0% | 2 | 6.07 ms / 97.6% | 0.00 ms / 0.0% | 0.15 ms / 2.4% | 1/2 | 0 |
| `tp1907_18_pt_element_inventory_point` | inner recovery transaction | 98.76 ms | 98.76 ms / 100.0% | 15 | 33.22 ms / 33.6% | 20.24 ms / 20.5% | 45.30 ms / 45.9% | 5/15 | 0 |
| `tp1907_18_ph_element_inventory_point` | inner fixed P,T | 123.75 ms | 122.81 ms / 99.2% | 98 | 122.31 ms / 98.8% | 0.00 ms / 0.0% | 0.50 ms / 0.4% | 49/98 | 0 |
| `rp1311_40_pt_element_inventory_point` | inner fixed P,T | 253.97 ms | 253.39 ms / 99.8% | 1 | 4.44 ms / 1.7% | 233.85 ms / 92.1% | 15.10 ms / 5.9% | 0 | 0 |
| `rp1311_40_ph_element_inventory_point` | inner fixed P,T | 132.52 ms | 128.64 ms / 97.1% | 49 | 127.70 ms / 96.4% | 0.00 ms / 0.0% | 0.94 ms / 0.7% | 0 | 0 |

#### Nonlinear solver callback table

This separate table is the solver-level evidence. It reports the accepted
backend callback work; the full rendered attempt table additionally retains
each failed and rejected backend row. RST does not expose its internal linear
solve, so `n/a` is an interface limitation, not zero work.

| case | backend | backend wall | iterations | residual calls / time | Jacobian calls / time | linear solves / time | overhead |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `stanjan15_ph_element_inventory_point` | RST | 2.95 ms | 10 | 10 / 1,516.5 us (51.5%) | 10 / 1,255.9 us (42.6%) | 10 / n/a | 173.8 us (5.9%) |
| `tp1907_18_ph_element_inventory_point` | legacy TR | 3.43 ms | 207 | 404 / 925.5 us (27.0%) | 212 / 1,485.3 us (43.3%) | 207 / 194.6 us (5.7%) | 820.1 us (23.9%) |
| `rp1311_40_pt_element_inventory_point` | RST | 4.44 ms | 10 | 8 / 928.1 us (20.9%) | 8 / 3,299.1 us (74.3%) | 10 / n/a | 210.9 us (4.8%) |
| `rp1311_40_ph_element_inventory_point` | RST | 69.32 ms | 174 | 160 / 10,693.1 us (15.4%) | 160 / 55,816.8 us (80.5%) | 174 / n/a | 2,784.7 us (4.0%) |
| `ph_element_inventory_point_legacy_lm` | legacy LM | 0.14 ms | 8 | 17 / 26.8 us (19.1%) | 8 / 16.4 us (11.7%) | 8 / 0.0 us (0.0%) | 96.9 us (69.2%) |

#### Representative backend-attempt ledger

The production output contains one row for every started attempt. These rows
show the important transitions from the current capture: failed and rejected
attempts are retained with their wall time, and the accepted attempt is shown
with callback detail. The full table includes the remaining seed/trial rows.

| case | scope | backend | outcome | metrics | wall | residual / Jacobian / linear / overhead |
| --- | --- | --- | --- | --- | ---: | --- |
| `stanjan15_ph_element_inventory_point` | monolithic | RST LM | failed | wall-only | 3.12 ms | n/a / n/a / n/a / n/a |
| `stanjan15_ph_element_inventory_point` | monolithic | RST Minpack LM | accepted | reported+callbacks | 2.95 ms | 1.52 / 1.26 / n/a / 0.17 ms |
| `tp1907_18_pt_element_inventory_point` | recovery/physical-failure | RST LM | failed | wall-only | 0.52 ms | n/a / n/a / n/a / n/a |
| `tp1907_18_pt_element_inventory_point` | recovery/normalized-discovery | Legacy LM | accepted | reported+callbacks | 0.13 ms | 0.02 / 0.00 / 0.01 / 0.10 ms |
| `tp1907_18_pt_element_inventory_point` | recovery/physical-retry | RST LM | accepted | reported+callbacks | 0.20 ms | 0.10 / 0.10 / n/a / 0.00 ms |
| `tp1907_18_ph_element_inventory_point` | trial-0 | RST LM | failed | wall-only | 0.20 ms | n/a / n/a / n/a / n/a |
| `tp1907_18_ph_element_inventory_point` | trial-0 | RST Minpack LM | rejected candidate | reported+callbacks | 4.84 ms | 2.24 / 2.02 / n/a / 0.58 ms |
| `tp1907_18_ph_element_inventory_point` | trial-0 | Legacy TR | accepted | reported+callbacks | 0.15 ms | 0.00 / 0.00 / 0.03 / 0.12 ms |
| `rp1311_40_pt_element_inventory_point` | point | RST LM | accepted | reported+callbacks | 4.44 ms | 0.93 / 3.30 / n/a / 0.21 ms |
| `rp1311_40_ph_element_inventory_point` | trial-0 | RST LM | accepted | reported+callbacks | 4.09 ms | 0.79 / 3.12 / n/a / 0.18 ms |

The accepted callback table explains solver work, while the attempt ledger
explains solver selection cost. For example, TP-1907 P,H has `98` started
attempts and `49/98` without structured metrics, but zero unknown wall time;
those calls remain visible in the full attempt table as measured `wall-only`
rows. Therefore the remaining optimization questions are concrete: RST
preparation for STANJAN/RP P,T, P,H trial preparation for RP-1311, and the
callback implementation itself only after the separate nonlinear-solver
optimization task is opened.

The lower-level tables rendered failed calls as `wall-only` with their actual
elapsed time and `n/a` only for unavailable internal counters. `unknown wall`
was zero across the capture. TP-1907 P,T recovery is explicitly visible as a
recovery transaction; RP-1311/40 P,T exposes its large RST preparation stage
and the remaining runner/orchestration remainder. Thus the report identifies
the next optimization owner instead of leaving a generic unexplained budget.

The post-fix focused check
`cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium --quiet`
completed with `979 passed; 0 failed; 71 ignored`. The three frozen cases that
were reported as failures in the parallel all-library run each pass when run
individually; their parallel failure is the pre-existing relative-path/global
library snapshot sensitivity, not an equilibrium solve or timing failure. The
same all-library run also reported seven unrelated failures in kinetics, GUI,
settings, and the legacy phase test.

## 73. PERF-RST-01: isolated RST solver comparison [HISTORICAL]

**Status:** Historical pre-prepared-adapter capture from 2026-09-12. This
section is retained for the original solver-compatibility evidence only. It
predates the parameterized prepared RST adapter, the production Sequential
Lambdify policy, and the addition of Backtracking LM to the benchmark matrix.
Do not use its wall times or four-method list as the current baseline; use
Stories 103, 104, 105, and 108 plus the current benchmark source.

### Story name

The RST preparation benchmark must show not only where the calculator spends
time, but also which RST method can produce an accepted equilibrium for the
same facade request. `LevenbergMarquardt`, `MinpackLevenbergMarquardt`,
`NielsenLevenbergMarquardt`, and `DampedNewton` are separate methods and must
not be represented by a production cascade row.

### Command

```powershell
$env:CARGO_PROFILE_BENCH_LTO='false'
$env:CARGO_PROFILE_BENCH_CODEGEN_UNITS='16'
cargo bench --bench equilibrium_rst_preparation --features equilibrium-instrumentation -- --noplot
```

The capture was made on 2026-09-12 in the optimized Cargo benchmark profile.
The diagnostic dashboard uses three fixed repetitions per case; Criterion
uses ten measurements per benchmark ID. `equilibrium_rst_preparation` uses
the normal `EquilibriumCalculator` facade and
`postprocessing_and_logging` instrumentation. Every row uses
`SolverPolicy::Single(RustedSciThe(...))`; legacy fallback is disabled.

### Question and hypothesis

Can the RST methods be compared on both acceptance and cost as the system grows
from 3 to 15, 18, and 40 components? The answer must be split into two
dimensions: accepted-solve time and solver compatibility. A short rejected
attempt is useful diagnostic evidence, but is not a faster equilibrium solve.

### Matrix

The benchmark contains 44 fixed-policy cases. All four RST methods run the
same molecular P,T size ladder (3, 15, 18, 40 components), small elemental
P,T, small monolithic elemental P,H, and 15-component elemental P,T. LM and
Minpack additionally cover molecular and elemental P,H, TP-1907/18 and
RP-1311/40 elemental points, and the three-point elemental continuation range.

### Result: molecular P,T solver comparison

The numbers below are Criterion median wall times. `failed` means the isolated
RST call did not publish an accepted solution; its measured rejection time is
not a successful-solve performance result.

| RST method | 3 components | 15 components | 18 components | 40 components |
| --- | ---: | ---: | ---: | ---: |
| LM | 1.057 ms, accepted | 15.268 ms, accepted | 22.776 ms, failed | 253.39 ms, failed |
| Minpack LM | 1.100 ms, accepted | 15.003 ms, accepted | 29.654 ms, failed | 252.47 ms, accepted |
| Nielsen LM | 1.426 ms, accepted | 34.532 ms, failed | 30.794 ms, failed | 259.03 ms, failed |
| Damped Newton | 1.092 ms, accepted | 28.625 ms, failed | 22.514 ms, failed | 251.84 ms, failed |

The 40-component molecular row identifies Minpack LM as the only accepted
method in this isolated slice. Damped Newton has a low rejection time there,
but that is not evidence that it is the best solver. LM remains the best
accepted choice on the 3/15-component molecular rows in this capture. Nielsen
LM and Damped Newton require further compatibility work before a production
default decision can be based on speed.

### Result: elemental and P,H controls

| Case | LM | Minpack LM | Nielsen LM | Damped Newton |
| --- | ---: | ---: | ---: | ---: |
| small elemental P,T | failed | 1.35 ms, accepted | 2.33 ms, accepted | failed |
| small elemental P,H | 21.89 ms, accepted | 21.98 ms, accepted | failed | 22.17 ms, accepted |
| STANJAN-15 elemental P,T | 16.65 ms, accepted | 16.26 ms, accepted | 17.98 ms, accepted | failed |

The P,H and elemental controls confirm that solver ranking is route- and
fixture-dependent. No method is promoted or demoted from this benchmark
alone; the isolated acceptance result must first pass the ordinary and frozen
correctness gates.

### Detailed timing evidence

The rendered report contains four separate views rather than one misleading
sum of overlapping percentages:

1. macro facade wall partition into `outside workflow` and `workflow`;
2. concrete facade and route stages, including thermochemistry, normalization,
   preparation, projection, validation, and publication;
3. every RST attempt, including rejected candidates and failed attempts;
4. execution/reuse counters and RST callback timing.

Representative accepted callback rows from the same capture were:

| Case | RST method | backend wall | iterations | residual | Jacobian | overhead |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| small molecular P,T | LM | 0.39 ms | 3 | 0.17 ms / about 43% | 0.20 ms / about 50% | about 0.02 ms |
| STANJAN-15 molecular P,T | LM | about 1 ms | 7 | about 0.4 ms | about 0.6 ms | about 0.05 ms |
| RP-1311/40 molecular P,T | Minpack LM | about 3.5 ms | 5 | about 0.8 ms | about 2.4 ms | about 0.32 ms |

RST does not currently publish a separate internal linear-system duration;
the report prints that field as `n/a`, never as zero. Failed and rejected
attempts remain in the ledger with their wall time and callback counts where
RST published metrics.

### Verification and decision

`cargo check --bench equilibrium_rst_preparation --no-default-features`
passed, `cargo fmt -- benches/equilibrium_rst_preparation.rs` completed, and
`git diff --check` reported no whitespace errors. The dashboard completed and
rendered all matrix tables. Criterion returned exit code 1 because its stored
statistical baseline flagged performance regressions in some existing IDs;
there was no benchmark panic or hidden fallback. This is a performance signal,
not a correctness pass.

The correctness gate remains open: run the focused calculator suite, matching
frozen references, molecular/elemental equivalence stories, and phase
lifecycle tests before using any solver ranking for optimization or defaults.

## 74. PERF-PREP-01: elemental P,T bridge reuse

### Story name

The elemental P,T facade must not prepare the same canonical
`PhaseEquilibriumProblemBundle` twice merely because it needs the bundle's
initial numerical seed. One bridge build is the ownership invariant for a
one-shot elemental P,T point; normalized recovery is allowed to build its own
transformed retry after a physical numerical failure.

### Commands

```powershell
cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::facade_element_inventory_pt_point_builds_canonical_bridge_once --quiet
cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests --quiet
$env:CARGO_PROFILE_BENCH_LTO='false'
$env:CARGO_PROFILE_BENCH_CODEGEN_UNITS='16'
cargo bench --bench equilibrium_rst_preparation --features equilibrium-instrumentation -- rp1311_40_pt_element_inventory_rst_minpack_levenberg_marquardt --noplot
```

### Question and hypothesis

Does passing the seed-producing bridge into the first physical transaction
remove the duplicate preparation without changing the accepted equilibrium,
recovery behavior, or diagnostic stage ownership? The expected invariant is
`phase_equilibrium_builds = 1` for the public elemental P,T point. Any wall-time
improvement is expected to be modest here because the dominant RP-1311/40 cost
is RST preparation, not bridge construction.

### Correctness result

The focused point test passed and printed no additional backend route. It
asserts `phase_equilibrium_builds() == 1` and conservation error below
`1e-7`. The full `equilibrium_calculator::tests` module passed: `24 passed,
0 failed, 1 ignored`. The existing elemental P,T range test also passed and
continues to assert one canonical bridge build and formulation reuse.

### Performance result

The optimized RP-1311/40 elemental P,T Minpack diagnostic observations were
`251.39/252.20/254.98 ms` before and `247.34/244.79/244.73 ms` after this
change. The concrete post-change stage sample was `request preparation
0.43 ms`, `RST preparation 228.06 ms`, and `nonlinear solve 246.20 ms`.
Criterion measured `[246.88 ms, 248.84 ms, 250.79 ms]` and reported
`No change in performance detected` (`p = 0.44`). The comparison therefore
confirms removal of duplicate bridge work and no correctness regression, but
does not claim a material benchmark speedup.

### Decision and next action

`PERF-PREP-01` is accepted as a correctness-preserving preparation cleanup.
It does not close the performance target for large cases. Continue with cold
versus repeated/prepared RST graph measurements; the next optimization must
attack the measured `RST preparation` budget while retaining the same isolated
solver policy and detailed attempt ledger.

## 75. PERF-PREP-02: prepared fixed-P,T session

### Story name

Repeated fixed-phase P,T calculations over one unchanged species universe must
reuse the prepared formulation and parameterized RST graph. The public facade
must expose this as a typed session so that preparation time is not paid again
for every temperature point.

### Commands

```powershell
cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::prepared_pt_session_reuses_rst_graph_across_sequential_points --quiet
cargo check --bench equilibrium_prepared_session --no-default-features --quiet
$env:CARGO_PROFILE_BENCH_LTO='false'
$env:CARGO_PROFILE_BENCH_CODEGEN_UNITS='16'
cargo bench --bench equilibrium_prepared_session --features equilibrium-instrumentation -- --noplot
```

### Question and hypothesis

Does a prepared session remove repeated symbolic RST construction while
preserving ordinary accepted-solution publication and continuation from the
last accepted seed? The expected evidence is one retained graph, a large
construction budget on the first point, and no symbolic construction on the
later point.

### Matrix and scope

`equilibrium_prepared_session` uses the normal `EquilibriumCalculator` facade
and the internal timing/instrumentation reports. The capture covers molecular
3-, 15-, and 40-component gas systems, elemental 3- and 15-component systems,
and the 18-component TP-1907 gas-plus-graphite ElementInventory fixture. Each
row compares preparation plus the first point with a later temperature point
on the same fixed-phase session. P,H and bounded phase-control routes are
intentionally excluded because their rollback and active-set lifecycles need a
separate prepared boundary.

### Correctness result

The focused sequential-session test passed. It asserts that the RST graph is
prepared, solves two temperatures through the session, and checks both
published solution temperatures and enabled timing on the second result. The
ordinary one-shot facade path was not changed by this experiment.

### Performance result

A representative diagnostic capture produced the following rows:

| Case | Prepared construction | Cold operation (construction + first point) | Hot point | Saved versus cold |
| --- | ---: | ---: | ---: | ---: |
| TP-1907/18 elemental P,T | 23.095 ms | 35.776 ms | 0.682 ms | 98.1% |
| RP-1311/40 molecular P,T | 267.132 ms | 275.111 ms | 2.881 ms | 99.0% |
| STANJAN-15 elemental P,T | 13.171 ms | 16.881 ms | 1.024 ms | 93.9% |

The concrete stage table attributes RP-1311/40 preparation to
`266.490 ms` of symbolic construction; its hot point has about `2.620 ms` in
the nonlinear solve and no new symbolic construction. The macro and concrete
tables are separate because the concrete stages are inclusive and must not be
summed as independent percentages.

A separate equal two-point capture compared three public
facade workflows. These are diagnostic wall samples, not Criterion medians:

| Case | One-shot loop | Prepared session | Built-in range | Session saving | Range saving |
| --- | ---: | ---: | ---: | ---: | ---: |
| small3 molecular P,T | 2.284 ms | 1.506 ms | 1.534 ms | 34.1% | 32.8% |
| STANJAN-15 molecular P,T | 30.835 ms | 16.369 ms | 15.839 ms | 46.9% | 48.6% |
| TP-1907/18 elemental P,T | 69.474 ms | 33.835 ms | 34.533 ms | 51.3% | 50.3% |
| RP-1311/40 molecular P,T | 495.888 ms | 231.799 ms | 245.918 ms | 53.3% | 50.4% |
| small3 elemental P,T | 2.984 ms | 2.225 ms | 2.217 ms | 25.4% | 25.7% |
| STANJAN-15 elemental P,T | 31.178 ms | 15.782 ms | 16.202 ms | 49.4% | 48.0% |

The result answers two different questions. Against a one-shot loop, both
prepared session and built-in range remove repeated symbolic construction and
save about `25-53%` for this two-point workload. Against each other, they are
close; the new session is therefore a reusable public ownership boundary, not
a claim that range internals were previously uncached.

### Decision and next action

`PERF-PREP-02` is accepted as a measured ownership/API optimization for fixed
P,T. The large result confirms that prepared formulation ownership is a real
performance target. The next work is a correctness matrix against frozen and
ElementInventory stories, followed by allocation/shared-buffer measurements;
P,H prepared ownership remains open and must not be inferred from this result.

## 76. PERF-PH-01: focused P,H preparation and continuation baseline

### Story name

P,H performance work must identify the owner of the large wall-time budget for
both monolithic and nested routes. The report must distinguish formulation and
trial preparation from inner execution and bracket control, and must retain the
backend-attempt evidence needed to interpret fallback cost.

### Commands

```powershell
cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::facade_ph_monolithic_range_reuses_prepared_formulation_without_changing_acceptance -- --nocapture
cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::facade_solves_pt_then_ph_without_manual_resolved_request_assembly -- --nocapture
cargo check --bench equilibrium_ph_preparation --no-default-features --quiet
cargo bench --bench equilibrium_ph_preparation --no-default-features --features equilibrium-instrumentation -- --noplot
```

The optimized capture was executed from the generated release benchmark
binary with the same `--noplot` arguments after the benchmark target was
built. The benchmark uses the public `EquilibriumCalculator` facade,
`EquilibriumTimingMode::Enabled`, shared calculator instrumentation, and
`tabled`; it does not construct equations or time benchmark-local solver
wrappers.

### Question and hypothesis

Where does cold P,H time go after the facade and workflow reports are made
complete? The hypothesis was that monolithic cases are dominated by one-time
RST graph preparation, while nested cases are dominated by repeated trial
template preparation or inner execution. Route selection and bracket control
were expected to be small. P,H range reuse was expected to show one prepared
formulation build followed by reuse without changing accepted results.

### Matrix

| Fixture | Components | Input | Route |
| --- | ---: | --- | --- |
| `small3_molecular` | 3 | molecular | monolithic |
| `stanjan15_molecular` | 15 | molecular | monolithic |
| `stanjan15_elemental` | 15 | `ElementInventory` | monolithic |
| `tp1907_18_elemental` | 18 | `ElementInventory` + graphite | nested |
| `rp1311_40_elemental` | 40 | `ElementInventory` | nested |

### Correctness result

Both focused facade tests passed. The monolithic range published three
accepted points with `formulation_builds=1`, `continuation_points=2`, and
`formulation_reuses >= 2`; every target had an accepted enthalpy error and a
finite temperature. The existing ElementInventory P,H point/range test also
passed with one canonical bridge build, two points, and one continuation point.
No production equation, solver order, fallback classification, or rollback
contract was changed by this baseline.

### Diagnostic result

The following is one optimized diagnostic capture. The detailed stage columns
are inclusive and therefore are not percentages to be summed together. The
only additive boundary is facade wall time versus the complete P,H workflow;
the nested trial and inner columns explain the workflow budget at the next
level.

| Case | Route | Wall ms | P,H total ms | Mono prep ms | Mono RST ms | Mono backend ms | Trial prep ms | Trial template ms | Trial inner ms | Bracket ms | Builds | Reuses | Attempts |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `small3_molecular` | mono | 20.384 | 20.168 | 0.139 | 18.565 | 1.132 | 0 | 0 | 0 | 0 | 1 | 0 | 2 |
| `stanjan15_molecular` | mono | 532.658 | 531.790 | 0.443 | 492.215 | 28.844 | 0 | 0 | 0 | 0 | 1 | 0 | 6 |
| `stanjan15_elemental` | mono | 515.889 | 514.857 | 0.613 | 497.955 | 5.999 | 0 | 0 | 0 | 0 | 1 | 0 | 2 |
| `tp1907_18_elemental` | nested | 162.752 | 159.879 | 0 | 0 | 0 | 19.949 | 19.917 | 139.457 | 0.272 | 1 | 10 | 98 |
| `rp1311_40_elemental` | nested | 271.188 | 252.112 | 0 | 0 | 0 | 225.594 | 225.587 | 26.258 | 0.050 | 1 | 2 | 6 |

The attempt table printed by the same command contains every accepted,
rejected, and failed inner backend attempt. Representative evidence is:

| Case | Attempt evidence | Accepted backend |
| --- | --- | --- |
| small molecular | RST LM rejected at `0.539 ms`; RST Minpack accepted at `0.518 ms` | RST Minpack LM |
| STANJAN-15 molecular | RST LM failed; Minpack and Trust Region were rejected; Nielsen and Damped Newton failed | Legacy LM |
| STANJAN-15 elemental | RST LM failed at `2.974 ms`; RST Minpack accepted at `2.887 ms` | RST Minpack LM |
| TP-1907/18 nested | `11` trials and `98` backend attempts retained, including RST and legacy fallbacks | shown per trial in attempt table |
| RP-1311/40 nested | `3` trials and `6` backend attempts retained | shown per trial in attempt table |

### Range evidence

| Case | Targets | Total ms | Mean point ms | Formulation builds | Formulation reuses | Continuation points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `small3_molecular` | 3 | 21.725 | 0.735 | 1 | 2 | 2 |
| `stanjan15_elemental` | 3 | 522.552 | 5.351 | 1 | 2 | 2 |
| `rp1311_40_elemental` | 2 | 483.658 | 231.887 | 1 | 42 | 1 |

The RP-1311/40 nested reuse count is aggregated over nested trial
formulations, so it is not a point count; it must not be compared directly to
the `2` point targets. The one-build result and continuation point are the
relevant range ownership evidence.

### Decision and next action

`PERF-PH-01` is accepted as a diagnostic baseline and correctness gate. The
large monolithic budgets are in RST preparation (`18.821/20.353 ms` for the
small case and `496.815/534.766 ms` for STANJAN-15 molecular), while nested
TP/RP budgets are respectively trial inner execution and trial-template
preparation. Bracket control is not a target at these sizes. This is not yet a
performance win and does not authorize nonlinear solver changes.

The next experiment described by this baseline was executed as `PERF-PH-02`
below. The legacy fallback attempt table remains part of this baseline and was
not removed or averaged away.

## 77. PERF-PH-02: prepared monolithic P,H session

### Story name

Repeated fixed-phase monolithic P,H targets should pay for thermochemistry,
the fixed formulation, and the parameterized RST graph once. Later targets
must retarget only scalar enthalpy parameters and the last accepted numerical
seed, while preserving accepted physical results and transactional continuation.

### Commands

```powershell
cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::prepared_ph_session_reuses_monolithic_formulation_across_targets -- --nocapture
cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::prepared_ph_element_inventory_session_preserves_acceptance -- --nocapture
cargo test --lib --no-default-features Thermodynamics::ChemEquilibrium::equilibrium_calculator::tests::facade_ph_monolithic_range_reuses_prepared_formulation_without_changing_acceptance -- --nocapture
cargo check --bench equilibrium_ph_preparation --no-default-features --quiet
cargo bench --bench equilibrium_ph_preparation --no-default-features --features equilibrium-instrumentation -- --noplot prepared_session_two_points
```

The benchmark uses the public facade and the shared opt-in instrumentation.
The capture prints the complete P,H stage/attempt/range tables and a separate
prepared-session table; Criterion also executes the three session benchmark
IDs: `one_shot_two_points`, `prepared_session_two_points`, and
`prepared_session_hot_point`.

### Question and hypothesis

Does a public prepared P,H ownership boundary remove repeated fixed-formulation
and RST graph construction? The hypothesis is that the first prepared point
has approximately the ordinary one-shot cost, the session preparation pays
that cost once, and a later point is reduced to hot numeric execution. The
accepted backend and physical acceptance contract must remain unchanged.

### Implementation

`PreparedEquilibriumPhSession` is exported from the production prelude. It is
implemented on top of the existing `PreparedPhContinuationState`; it does not
introduce a second equation system or alter solver order. The initial scope is
deliberately fixed declared phases with explicit `PhSolveMode::Monolithic`.
Bounded phase control, `Auto` nested fallback, nested P,H, and cross-identity
cache sharing remain outside this session.

The session retargets its owned request in place and advances its composition
and temperature seed only after a successful accepted point. Failed solves do
not replace the continuation seed. Both molecular and ElementInventory
correctness tests passed.

### Performance result

One optimized diagnostic capture produced these two-target rows:

| Case | One-shot two points | Session preparation | Session first point | Session hot point | Session total | Reduction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `small3_molecular` | 42.789 ms | 19.817 ms | 1.454 ms | 0.407 ms | 21.678 ms | 49.3% |
| `stanjan15_molecular` | 1185.978 ms | 542.876 ms | 29.699 ms | 4.505 ms | 577.080 ms | 51.3% |
| `stanjan15_elemental` | 1127.899 ms | 545.127 ms | 6.485 ms | 6.018 ms | 557.630 ms | 50.6% |

The session first and second points each report `fixed_formulation_builds=0`
and `fixed_formulation_reuses=1`. The session total is the sum of preparation,
first point, and hot point for this capture; the reduction is against the
two-point one-shot loop. These are diagnostic wall samples, not a claim about
Criterion confidence intervals.

### Correctness result

The two session tests passed. They assert accepted enthalpy error, finite
temperature, and formulation reuse for molecular input; the molecular session
also tries an unreachable target first and then solves the original target,
proving that a failed point does not replace the continuation seed. The
ElementInventory test additionally asserts element balance. The monolithic
range regression also passed with one build and two continuation reuses. The
calculator test module remained green: `28 passed, 1 ignored`.

### Decision and next action

`PERF-PH-02` is accepted as a correctness-preserving performance/API
optimization. It demonstrates approximately `49-51%` reduction for this
two-target workload by removing the second one-time preparation. It does not
optimize nonlinear residual/Jacobian/linear-system execution.

Next, measure the internal split of the one-time session preparation and the
remaining nested P,H trial-template cost. Only shared immutable payloads or
large buffers with a material measured effect should be changed; the existing
legacy attempt ledger, fallback ordering, scale envelope, and frozen lifecycle
tests remain mandatory gates.

## 78. PERF-PH-03: P,H preparation and nested trial detail

### Story name

The P,H performance dashboard must explain cold prepared-session construction
and every nested temperature trial without hiding work inside one preparation
or workflow number.

### Command

```powershell
cargo bench --bench equilibrium_ph_preparation --no-default-features --features equilibrium-instrumentation -- --noplot prepared_session_two_points
```

The benchmark uses the public facade and the opt-in shared instrumentation.
It prints the existing macro/stage/attempt/session tables, a cold monolithic
preparation breakdown, and one row for every nested TP-1907/18 and RP-1311/40
temperature trial. The capture below was run from the optimized benchmark
executable after the bench target was linked; the package-level Cargo command
also started rebuilding the unrelated GUI binary.

### Question and hypothesis

Which part of the one-time prepared P,H state is responsible for its large
cold wall time, and does nested P,H spend its time in template construction,
scalar orchestration, or inner equilibrium execution? The hypothesis was that
symbolic expression construction plus RST payload materialization dominate
cold monolithic preparation, while nested trial rows expose whether the
one-time template is amortized and whether the remaining cost belongs to the
inner nonlinear solve.

### Cold monolithic preparation output

| case | facade setup | state total | thermochemistry | normalization | equation | symbolic construction | P,H formulation | RST preparation | state remainder |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `small3_molecular` | 0.123 ms | 25.005 ms | 0.693 ms | 0.000 ms | 0.003 ms | 5.710 ms | 0.020 ms | 18.497 ms | 0.069 ms |
| `stanjan15_molecular` | 0.336 ms | 756.420 ms | 5.757 ms | 0.000 ms | 0.005 ms | 517.481 ms | 0.070 ms | 232.854 ms | 0.191 ms |
| `stanjan15_elemental` | 0.876 ms | 751.531 ms | 5.284 ms | 0.000 ms | 0.008 ms | 513.946 ms | 0.128 ms | 231.744 ms | 0.326 ms |

The stage values are measured inside one preparation invocation. The facade
setup is the measured wall outside the prepared state; the state remainder is
the explicit complement after the named exclusive preparation fields.
`symbolic construction` covers the transformed Gibbs and enthalpy expression
graph; `RST preparation` covers RST symbolic-problem and lambdify/backend
materialization. The values are ownership evidence, not a promise that
independently printed wall samples will be identical.

### Nested P,H output

The existing macro row gave this aggregate ownership view for the same
optimized executable:

| case | trials | template preparation | trial preparation | inner execution | bracket orchestration |
| --- | ---: | ---: | ---: | ---: | ---: |
| `tp1907_18_elemental` | 11 | 20.391 ms | 20.424 ms | 138.731 ms | 0.210 ms |
| `rp1311_40_elemental` | 3 | 232.878 ms | 232.887 ms | 25.556 ms | 0.072 ms |

The separate `PERF-PH-03 nested P,H trial detail` table contains all 14
trial rows and, for each row, prints temperature, scalar step, preparation
kind, template-build marker, trial total, inner execution/equilibrium,
enthalpy, postprocessing, and inner equation/RST/nonlinear/validation times.
Representative rows from that table were:

| case/trial | preparation | template build | trial total | inner execution | nonlinear solve | enthalpy |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `tp1907_18_elemental/0` | `FixedFormulationInitial` | true | 12.806 ms | 12.780 ms | 12.622 ms | 0.023 ms |
| `tp1907_18_elemental/10` | `FixedFormulationReused` | false | 12.773 ms | 12.759 ms | 12.665 ms | 0.012 ms |
| `rp1311_40_elemental/0` | `FixedFormulationInitial` | true | 4.346 ms | 4.301 ms | 4.034 ms | 0.043 ms |
| `rp1311_40_elemental/1` | `FixedFormulationReused` | false | 12.918 ms | 12.884 ms | 12.673 ms | 0.031 ms |

### Correctness result

The calculator module remained green: `28 passed, 1 ignored`. The prepared
molecular P,H session, its failed-target retry/rollback check, and the
ElementInventory acceptance test all passed. The benchmark's three prepared
session Criterion IDs (`small3_molecular`, `stanjan15_molecular`, and
`stanjan15_elemental`) completed successfully.

### Decision and next action

`PERF-PH-03` closes the diagnostic gap for the next optimization decision.
The large cold monolithic cost is not in normalization, equation assembly, or
P,H formulation. It is split between symbolic expression construction and RST
payload materialization. On STANJAN-15, symbolic construction accounts for
about `68.4%` of state preparation and RST materialization for about `30.8%`;
on the 3-component case the corresponding shares are about `22.8%` and
`74.0%`. Thermochemistry remains small but is measurable at roughly `5-6 ms`
on STANJAN-15.
For nested TP-1907/18, the template is paid once and subsequent rows are
nonlinear-execution dominated. For RP-1311/40, the first template build is
the dominant cost and is already reused by later trials. Bracket control and
enthalpy evaluation are not performance targets.

The next implementation experiment should first isolate immutable symbolic
expression construction from request-specific numeric state, then measure RST
backend materialization separately. A shared payload across P,T and P,H is
only valid if the differing equation dimensions and temperature/enthalpy
parameters are represented explicitly; the current P,H graph cannot be
reused by identity with a P,T graph. No callback parallelization or
solver-order change is justified by this story alone.

## 79. PERF-PH-04: Remove redundant P,H symbolic simplification

### Story name

Monolithic P,H RST preparation must not run a second full symbolic
simplification pass over already simplified physical residuals or over the
scaled enthalpy row.

### Question and hypothesis

Does removing the redundant `.simplify()` work at the P,H builder boundary
reduce cold preparation and end-to-end P,H time without changing equations,
accepted backend, attempt order, or physical acceptance? The hypothesis was
that the large symbolic trees make this repeated pass material, especially for
15-component systems.

### Change

`prepare_rst_symbolic_ph_problem_with_timing` still builds exactly the same
scaled residual and enthalpy equations. It now applies the validated constant
row scale and enthalpy scale directly, leaving the existing per-residual
generator simplification and RST differentiation unchanged. No nonlinear
solver, cascade order, tolerance, or fallback behavior was changed.

### Command and paired release capture

```powershell
cargo build --release --bench equilibrium_ph_preparation `
  --no-default-features --features equilibrium-instrumentation
```

The old and new release bench executables were run three times with
`prepared_session_two_points`; the table reports the median `ph_total_ms` from
the `PERF-PH-01` macro row. The same case matrix and public facade were used
on both sides.

| case | before median | after median | change |
| --- | ---: | ---: | ---: |
| `stanjan15_molecular` | 545.701 ms | 383.452 ms | -29.7% |
| `stanjan15_elemental` | 522.100 ms | 353.739 ms | -32.2% |

The paired samples retained the same accepted routes: legacy LM for the
molecular case and RST Minpack LM for the ElementInventory case. The RST
rejections/failures and fallback order remained visible in the attempt table.
The measurement is a wall-time experiment rather than a claim that every
machine will reproduce the same absolute milliseconds.

### Correctness result

The ordinary calculator module passed `28 passed, 1 ignored`. TP-1906 frozen
P,H molecular/elemental/formal route equivalence passed, including the same
recovered temperature, enthalpy error, residual, and element balance; the
nested/Auto physical-state agreement test also passed.

### Decision and next action

`PERF-PH-04` is accepted as a correctness-preserving preparation optimization:
the paired large-case medians show a material reduction, while nonlinear
callback execution was not modified. The next experiment remains immutable
symbolic expression reuse across compatible prepared sessions. It must keep
P,T and P,H graph dimensions separate unless an explicit shared RST payload
contract is introduced and measured.

## 80. PERF-PH-05: Skip intermediate P,H residual simplification

### Story name

The monolithic P,H RST route should pass its already-structured residuals
directly to RST instead of simplifying every physical residual before RST
differentiates and materializes the final payload.

### Command and question

```powershell
cargo build --release --bench equilibrium_ph_preparation `
  --no-default-features --features equilibrium-instrumentation
```

The release executable was run three times with
`--noplot prepared_session_two_points`. Does removing only the intermediate
P,H residual-generator simplification preserve the accepted routes and
reduce the cold preparation budget? The comparison control is the
post-PERF-PH-04 release capture, before this P,H generator change. The table
uses the median of three current runs; absolute wall values are machine and
load dependent.

### Hypothesis and implementation

The P,H RST path now uses a private raw generator mode. Fixed-P,T continues
to use the canonical simplified generator. The P,H path still performs the
final row scaling, adds the enthalpy equation, and lets RST own the
subsequent symbolic differentiation. Equation dimensions, variables,
parameters, tolerances, solver order, and fallback behavior are unchanged.

### Cold preparation evidence

| case | PERF-PH-04 control state total | PERF-PH-05 median state total | symbolic median | RST median | state change |
| --- | ---: | ---: | ---: | ---: | ---: |
| `stanjan15_molecular` | 383.825 ms | 233.075 ms | 10.763 ms | 216.483 ms | -39.3% |
| `stanjan15_elemental` | 461.680 ms | 227.570 ms | 10.673 ms | 210.912 ms | -50.7% |

The three current runs were:

| case | run 1 state / symbolic / RST | run 2 state / symbolic / RST | run 3 state / symbolic / RST |
| --- | --- | --- | --- |
| `stanjan15_molecular` | 223.500 / 6.956 / 211.898 ms | 235.466 / 11.339 / 218.069 ms | 233.075 / 10.763 / 216.483 ms |
| `stanjan15_elemental` | 227.570 / 10.673 / 210.912 ms | 224.970 / 11.147 / 208.026 ms | 230.591 / 10.353 / 214.509 ms |

RST materialization remains the dominant current owner at roughly
`~210-218 ms`; this experiment removes the former symbolic-construction
spike rather than claiming that the whole P,H preparation is solved.

### Correctness and route invariants

The following gates passed after the change:

```text
calculator tests: 28 passed, 1 ignored
TP-1906 frozen P,H molecular/elemental/formal route reuse: passed
TP-1906 nested monolithic versus Auto physical state: passed
TP-1907 frozen P,T molecular/elemental/formal route match: passed
numeric/legacy closure versus symbolic residual contract: passed
```

The prepared-session benchmark completed successfully for `small3_molecular`,
`stanjan15_molecular`, and `stanjan15_elemental`. Accepted routes and backend
attempt evidence remained unchanged: legacy LM for the molecular case and
RST Minpack LM for the ElementInventory case.

The full repository command `cargo test --lib --no-default-features` was also
attempted, but ended with `2387 passed, 10 failed, 154 ignored`. The failures
were outside this change (kinetics, GUI, settings, phase-spec, and concurrent
local-library fixture access). The three affected frozen tests pass when run
serially, including the P,H reuse and lifecycle cases listed above.

### Decision and next action

`PERF-PH-05` is accepted as a localized, correctness-preserving P,H
preparation optimization. It does not alter the fixed-P,T symbolic route or
the nonlinear solver itself. The next performance target is RST payload
materialization and immutable expression reuse at that boundary; any cache
must keep P,T and P,H graph identities separate and must be validated against
the same frozen and route-invariant gates.

## 81. PERF-PH-06: Keep one production RST symbolic payload

### Story name

The production RST adapter should retain only the callable symbolic payload
needed by the nonlinear solver. It should not retain a second copy of the
equations and symbolic Jacobian that is useful to the dependency's broader
diagnostic/AOT API but is not read by KiThe during a solve.

### Command and question

```powershell
cargo build --release --bench equilibrium_ph_preparation `
  --no-default-features --features equilibrium-instrumentation
```

Run the resulting executable three times with
`--noplot prepared_session_two_points`. Does the RST preparation boundary
remove the outer duplicate symbolic ownership while preserving the same
residual/Jacobian closures, solver cascade, accepted backend, and physical
result?

### Hypothesis and implementation

`RstSymbolicProblem` now owns one RST `Jacobian`, its lambdified residual and
Jacobian closures, variable names, and parameter values. It invokes the same
`Jacobian::calc_jacobian` and parameterized/non-parameterized lambdification
methods as the previous `SymbolicNonlinearProblem` route. The production path
therefore changes ownership and materialization overhead only; it does not
change equations, differentiation, scaling, variables, tolerances, solver
order, fallback behavior, or the separate P,T/P,H graph contracts.

### Release evidence

The table shows `facade_setup / state_total / thermochemistry / normalization /
equation / symbolic / P,H formulation / RST / state_unattributed`, in ms.
The last column is only the remainder after the named stages and is not part
of the RST share.

| case | run 1 | run 2 | run 3 |
| --- | --- | --- | --- |
| `stanjan15_molecular` | 0.342 / 95.692 / 3.173 / 0.000 / 0.011 / 6.961 / 0.142 / 85.067 / 0.275 | 0.317 / 106.540 / 3.262 / 0.000 / 0.005 / 6.836 / 0.069 / 96.096 / 0.210 | 0.377 / 97.766 / 3.400 / 0.000 / 0.010 / 6.935 / 0.185 / 86.755 / 0.320 |
| `stanjan15_elemental` | 0.438 / 109.782 / 3.221 / 0.000 / 0.007 / 7.458 / 0.121 / 98.592 / 0.281 | 0.453 / 101.216 / 3.103 / 0.000 / 0.006 / 6.850 / 0.088 / 90.791 / 0.313 | 0.450 / 93.542 / 3.530 / 0.000 / 0.006 / 7.038 / 0.231 / 82.448 / 0.223 |

Median comparison with the accepted `PERF-PH-05` control:

| case | PERF-PH-05 control state | PERF-PH-06 median state | PERF-PH-05 RST | PERF-PH-06 RST | state change |
| --- | ---: | ---: | ---: | ---: | ---: |
| `stanjan15_molecular` | 233.075 ms | 97.766 ms | 216.483 ms | 86.755 ms | -58.1% |
| `stanjan15_elemental` | 227.570 ms | 101.216 ms | 210.912 ms | 90.791 ms | -55.5% |

The RST preparation stage remains the largest named owner, but its median
falls by `59.9%` for the molecular case and `57.0%` for ElementInventory.
This is evidence for a materialization/ownership cost, not evidence that the
nonlinear solve became faster: solver execution is outside this cold
preparation measurement.

### Correctness and route invariants

The following focused gates passed after the adapter change:

```text
cargo test --lib --no-default-features equilibrium_rst_backend_tests -- --quiet
3 passed

cargo test --lib --no-default-features \
  prepared_ph_session_reuses_monolithic_formulation_across_targets -- --quiet
1 passed

cargo test --lib --no-default-features \
  tp1906_frozen_ph_molecular_elemental_and_formal_routes_reuse_target -- --quiet
1 passed

cargo test --lib --no-default-features \
  tp1907_frozen_pt_molecular_elemental_and_formal_routes_match -- --quiet
1 passed
```

The prepared-session benchmark reported `Success` for `small3_molecular`,
`stanjan15_molecular`, and `stanjan15_elemental`. The accepted backend routes
and cascade evidence remained unchanged: legacy LM for the molecular case and
RST Minpack LM for the ElementInventory case.

### Decision and next action

`PERF-PH-06` is accepted as a correctness-preserving RST preparation
optimization. It reduces duplicate dependency-side symbolic ownership without
altering the production equations or nonlinear solver. The next experiment is
safe reuse of an already materialized immutable RST `Jacobian`/expression
payload across compatible prepared sessions, only where graph identity,
parameter ownership, and P,T versus P,H dimensions are explicit. Nonlinear
callback optimization remains a separate benchmark family.

## 82. PERF-PH-07: Share immutable P,H retarget payloads

### Story name

Prepared monolithic P,H sessions should not deep-copy large immutable
thermochemistry capabilities or row-scale buffers when a new enthalpy target
is retargeted. Only target-local numerical state and the scalar energy
parameters should change.

### Command and question

```powershell
cargo build --release --bench equilibrium_ph_preparation `
  --no-default-features --features equilibrium-instrumentation
```

Run the resulting executable three times with
`--noplot prepared_session_two_points` and inspect the
`PERF-PH-01 prepared monolithic P,H session` table. Question: does sharing the
immutable payload preserve the prepared-session contract and reduce hot-point
retarget overhead without changing the nonlinear result?

### Hypothesis and implementation

`PreparedPhFormulation` now stores `Arc<ResolvedThermochemistry>` and
`Arc<Vec<f64>>` for the P,T row scales. `retarget` and the immutable runner
clone increment reference counts instead of cloning the complete capability
bundle and scale vector. The fixed-active numerical problem remains owned by
the retargeted formulation because it carries the current conditions, Gibbs
closures, and continuation seed. No residual, Jacobian, temperature transform,
enthalpy scaling, solver order, or publication rule changed.

### Release evidence

The prepared-session columns are `one_shot_two_points / preparation /
first_point / second_point`, in ms. The before values are the available single
PERF-PH-06 release control capture; the after values are medians of three fresh
PERF-PH-07 captures. The control is therefore directional evidence rather than
a statistically matched benchmark series.

| case | PERF-PH-06 control | PERF-PH-07 median | preparation change | first point change | second point change |
| --- | --- | --- | ---: | ---: | ---: |
| `stanjan15_molecular` | 264.617 / 93.219 / 35.088 / 5.713 | 257.339 / 85.618 / 28.758 / 4.559 | -8.2% | -18.0% | -20.2% |
| `stanjan15_elemental` | 229.317 / 105.706 / 6.531 / 8.934 | 208.498 / 95.541 / 6.175 / 4.547 | -9.6% | -5.5% | -49.1% |

All three post-change captures completed the prepared-session benchmark with
`Success` for `small3_molecular`, `stanjan15_molecular`, and
`stanjan15_elemental`. The cold preparation and first-point numbers remain
solver/load-sensitive; the meaningful architectural evidence is that the
retarget path now shares the large immutable bundle, while the second-point
hot path did not regress in either case.

### Correctness and ownership evidence

```text
cargo test --lib --no-default-features equilibrium_ph_formulation -- --quiet
12 passed

cargo test --lib --no-default-features \
  prepared_ph_session_reuses_monolithic_formulation_across_targets -- --quiet
1 passed

cargo test --lib --no-default-features \
  tp1906_frozen_ph_molecular_elemental_and_formal_routes_reuse_target -- --quiet
1 passed

cargo test --lib --no-default-features \
  tp1907_frozen_pt_molecular_elemental_and_formal_routes_match -- --quiet
1 passed
```

The formulation test also asserts pointer identity for thermochemistry and
row-scale storage before and after retarget. This protects the intended
ownership property directly; numerical frozen tests protect the physical
contract.

### Decision and next action

`PERF-PH-07` is retained as a narrow ownership optimization, with the caveat
that its before/after timing control is not a matched multi-run experiment.
The next measurement should target the remaining cloned structural
`PreparedEquilibriumProblem` in P,H retargeting. It must first prove that this
copy is material for 15/18/40-component cases and preserve transactional
rollback before any ownership change. Nonlinear callback work remains a
separate benchmark family.

## 83. PERF-PH-08: Retarget structural P,H state (rejected)

### Story name

Remove the per-target clone of `PreparedEquilibriumProblem` from the prepared
monolithic P,H session, while preserving transactional rollback and the exact
fixed-active physical contract.

### Command and question

```powershell
cargo build --release --bench equilibrium_ph_preparation `
  --no-default-features --features equilibrium-instrumentation
```

Run the refreshed executable three times with
`--noplot prepared_session_two_points`. Is the structural retarget clone a
material hot-path cost for the current 15-component molecular and
ElementInventory cases?

### Hypothesis and attempted implementation

The attempted change moved the fixed-active `PreparedEquilibriumProblem` into
an `Arc` owned by `PreparedPhFormulation`, retargeting only target-local P,H
values. The formulation, workflow, and frozen tests were kept unchanged in
their equations and solver behavior. The change was reverted after the
measurement because the result did not establish a repeatable gain.

### Release evidence

The prepared-session columns are `one_shot_two_points / preparation /
first_point / second_point`, in ms. The PERF-PH-07 values are the preceding
three-run control; PERF-PH-08 values are three runs with the structural-share
attempt.

| case | PERF-PH-07 control median | PERF-PH-08 attempt median | preparation change | first point change | second point change |
| --- | --- | --- | ---: | ---: | ---: |
| `stanjan15_molecular` | 257.339 / 85.618 / 28.758 / 4.559 | 261.069 / 95.152 / 29.633 / 4.738 | +11.1% | +3.0% | +3.9% |
| `stanjan15_elemental` | 208.498 / 95.541 / 6.175 / 4.547 | 214.751 / 97.203 / 6.189 / 4.320 | +1.7% | +0.2% | -5.0% |

The attempt did not improve cold preparation or first-point execution. The
ElementInventory second point improved modestly, while the molecular second
point regressed; this is insufficient evidence against run-to-run noise and
does not isolate the structural copy from solver variance.

### Correctness and rollback

During the attempt the focused gates passed: formulation `12/12`, workflow
`35/35`, RST adapter `3/3`, prepared P,H `1/1`, TP-1906 P,H `1/1`, and TP-1907
P,T `1/1`. The attempt was nevertheless rolled back because the performance
acceptance criterion was not met. The rejected measurement remains recorded;
the production code retains the accepted PERF-PH-07 ownership change.

### Decision and next action

`PERF-PH-08` is rejected and reverted. Before revisiting this copy, add a
matched before/after allocation or profiler measurement for 15, 18, and 40
components, with structural-copy time separated from nonlinear execution.
Do not change the transactional ownership model based on the current wall
table alone. Nonlinear callback work remains a separate benchmark family.

## 84. PERF-PH-09: Isolate prepared P,H structural retarget cost

### Story name

Measure the remaining `PreparedEquilibriumProblem::retarget_with_gibbs`
structural clone and validation inside a reusable monolithic P,H session,
without attributing that work to the nonlinear solver.

### Command and question

```powershell
cargo bench --bench equilibrium_ph_preparation `
  --no-default-features --features equilibrium-instrumentation -- --noplot
```

The benchmark now prints a dedicated retarget table with
`one_shot_two_points / preparation / first_point / second_point /
first_retarget / second_retarget` in ms, plus an explicit route-availability
table. The question is whether the remaining structural copy is large enough
to justify another transactional ownership experiment.

### Hypothesis and scope

The measurement must cover small, 15-component, and 40-component prepared
monolithic routes. TP-1907 is an 18-component ElementInventory fixture, but
its current production P,H route is nested; forcing fixed-monolithic mode is a
rejected numerical route, so it is retained in the nested trial table rather
than reported as a misleading zero-cost retarget row.

### Capture evidence

The release timing capture completed for `small3_molecular`,
`stanjan15_molecular`, `stanjan15_elemental`, and `rp1311_40_elemental`; the
new retarget columns were non-zero for both prepared points. The table below
is `one-shot-two-points / preparation / first-point / second-point /
first-retarget / second-retarget`, all in ms:

| case | one-shot two points | preparation | first point | second point | first retarget | second retarget |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `small3_molecular` | 17.599 | 7.249 | 1.164 | 0.434 | 0.005 | 0.003 |
| `stanjan15_molecular` | 244.457 | 81.450 | 29.139 | 4.150 | 0.009 | 0.006 |
| `stanjan15_elemental` | 204.823 | 99.116 | 6.404 | 4.453 | 0.009 | 0.006 |
| `rp1311_40_elemental` | 2495.210 | 1146.921 | 51.618 | 10.131 | 0.014 | 0.018 |

The separate allocation-enabled release capture reported the following
allocation calls/bytes for the first and second accepted points:

| case | first calls | first bytes | second calls | second bytes |
| --- | ---: | ---: | ---: | ---: |
| `small3_molecular` | 1686 | 200591 | 1073 | 123814 |
| `stanjan15_molecular` | 38242 | 6487512 | 21319 | 3098061 |
| `stanjan15_elemental` | 6311 | 1232544 | 21320 | 3099587 |
| `rp1311_40_elemental` | 16263 | 8518743 | 15135 | 2718725 |

These allocation totals include thermochemistry, nonlinear execution, and
publication; they are evidence for the size of the point-wide allocation
budget, not an attribution of every allocation to `retarget_with_gibbs`.
The retarget-specific timing is currently tens of microseconds, while cold
preparation and nonlinear work are orders of magnitude larger. This does not
justify another ownership experiment.

### Correctness gate

```powershell
cargo test --lib --no-default-features equilibrium_ph_workflow -- --quiet
cargo test --lib --no-default-features equilibrium_ph_formulation -- --quiet
cargo check --bench equilibrium_ph_preparation --no-default-features `
  --features equilibrium-instrumentation --quiet

cargo bench --bench equilibrium_ph_preparation `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-allocation" `
  -- --noplot prepared_session_two_points
```

The prepared-session test also asserts that enabled reports publish a
non-zero retarget duration for both the first and second accepted points.
Timing remains opt-in: disabled execution does not read the clock for this
stage.

### Decision and next action

`PERF-PH-09` is accepted as measurement infrastructure and a negative
materiality result for the structural-copy hypothesis. Do not repeat the
rejected `Arc<PreparedEquilibriumProblem>` experiment. The 18-component
TP-1907 result belongs to the nested-route baseline and is a separate
stability/solver case; its fixed-monolithic `N/A` status is intentional.

## 85. PERF-PH-10: Reuse the accepted Gibbs snapshot

### Story name

Reuse the immutable Gibbs closure snapshot produced for the previous accepted
prepared P,H point when the next point starts at that same continuation
temperature.

### Command and question

```powershell
cargo bench --bench equilibrium_ph_preparation `
  --no-default-features --features equilibrium-instrumentation `
  -- --noplot prepared_session_two_points
```

The prepared-session table now includes `first_gibbs_ms`, `second_gibbs_ms`,
and `first_gibbs_reused / second_gibbs_reused`. The question is whether the
Gibbs evaluation and closure materialization can be removed from the hot
continuation point without changing the accepted physical result.

### Hypothesis and implementation

After an accepted point, the workflow already evaluates a constant Gibbs
snapshot for publication. The prepared continuation now retains that snapshot
only at the accepted-result commit boundary. A later point reuses it when its
initial temperature is bitwise identical to the cached temperature; otherwise
it evaluates a fresh snapshot. Failed targets cannot warm or replace this
cache, preserving transactional continuation semantics.

### Release evidence

The filtered release capture reported the following values in ms. The Gibbs
columns are `first-point preparation / second-point cache clone`; the boolean
columns prove the lifecycle transition:

| case | first point | second point | first Gibbs | second Gibbs | first reused | second reused |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `small3_molecular` | 1.066 | 0.506 | 0.024 | 0.000 | false | true |
| `stanjan15_molecular` | 29.563 | 4.634 | 0.080 | 0.000 | false | true |
| `stanjan15_elemental` | 5.994 | 4.863 | 0.128 | 0.000 | false | true |
| `rp1311_40_elemental` | 52.656 | 10.437 | 0.215 | 0.000 | false | true |

The same capture retained successful prepared-session routes for all four
cases. The failed-target-then-success test also proves that a rejected target
does not make the first accepted point report a false cache hit. This is a
narrow preparation optimization; it does not parallelize or alter nonlinear
callbacks, and it does not cache target-dependent enthalpy scaling.

### Correctness gate

```powershell
cargo test --lib --no-default-features `
  prepared_ph_session_reuses_monolithic_formulation_across_targets -- --quiet
cargo test --lib --no-default-features equilibrium_ph_workflow -- --quiet
cargo check --bench equilibrium_ph_preparation --no-default-features `
  --features equilibrium-instrumentation --quiet
```

### Decision and next action

`PERF-PH-10` is accepted as a cache-localization optimization. The next
measurement should determine whether repeated P,H request validation and
initial enthalpy-vector evaluation are material; any cache must respect the
changing continuation temperature and failed-point rollback boundary.

## 86. PERF-PH-11: Remove repeated prepared P,H validation and enthalpy traversal

### Story name

Measure and minimize request validation and initial molar-enthalpy preparation
at the prepared P,H continuation boundary without reusing thermochemistry at a
different seed temperature.

### Command and question

```powershell
cargo bench --bench equilibrium_ph_preparation `
  --no-default-features --features equilibrium-instrumentation `
  -- --noplot prepared_session_two_points
```

The prepared-session table now reports
`first/second_request_validation_ms`,
`first/second_initial_enthalpy_ms`, and their reuse flags. The question is
whether request validation or construction of the initial molar enthalpy
vector is large enough to justify a continuation cache.

### Hypothesis and implementation

The old prepared-session boundary performed full request validation before the
prepared workflow, then validated again while preparing the nonlinear point.
The retarget operation now validates only structural fields and restores the
previous request if that check fails. The prepared workflow evaluates or clones
one initial enthalpy vector, validates against that vector, and commits its
temperature-keyed cache only after an accepted result.

Exact temperature matching is required. A changed continuation temperature
must evaluate a fresh vector, and a failed target must not update the cache.
The same transaction rule as the Gibbs cache therefore remains in force.

### Release evidence

The filtered release capture reported the following values in ms. A `true`
reuse flag means the vector was cloned from the prepared continuation cache;
`0.000` validation means it was below the displayed millisecond precision.

| case | first point | second point | first validation | second validation | first enthalpy | second enthalpy | first reused | second reused |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `small3_molecular` | 1.046 | 0.449 | 0.000 | 0.000 | 0.000 | 0.001 | true | false |
| `stanjan15_molecular` | 30.171 | 5.040 | 0.000 | 0.000 | 0.000 | 0.006 | true | false |
| `stanjan15_elemental` | 6.649 | 4.760 | 0.000 | 0.000 | 0.000 | 0.007 | true | false |
| `rp1311_40_elemental` | 51.864 | 10.556 | 0.000 | 0.000 | 0.000 | 0.020 | true | false |

The first point begins at the prepared session temperature and reuses the
vector captured during cold preparation. The second point changes its seed
temperature in these fixtures and correctly performs a fresh evaluation. This
is a correctness signal, not an invitation to reuse a vector across
temperatures.

### Correctness gate

```powershell
cargo test --lib --no-default-features `
  prepared_ph_session_reuses_monolithic_formulation_across_targets -- --quiet
cargo test --lib --no-default-features equilibrium_ph_workflow -- --quiet
cargo check --bench equilibrium_ph_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-allocation" --quiet
```

Observed gates:

```text
prepared_ph_session_reuses_monolithic_formulation_across_targets: 1 passed
equilibrium_ph_workflow: 35 passed
equilibrium_calculator: 28 passed, 1 ignored
nasa_tp1907_chon_graphite_tests: 22 passed, 6 ignored
nasa_tp1906_chon_graphite_ph_tests: 6 passed, 6 ignored
bench allocation configuration: cargo check passed
```

### Decision and next action

`PERF-PH-11` is accepted as a de-duplication and diagnostics improvement. The
measured validation path is below materiality, and initial enthalpy traversal
is only `0.001–0.020 ms` in the 3/15/40-component fixtures. Do not add a more
general thermochemistry workspace cache based on this evidence. The next
performance work remains the separately instrumented nonlinear callback path;
prepared P,H orchestration is no longer the leading target in these cases.

## 87. PERF-LEGACY-01: Isolated Legacy callback budget

### Story name

As a performance engineer, I need residual, analytic Jacobian, linear-system,
and solver-overhead timing for each retained Legacy nonlinear backend so that
the next optimization targets measured work rather than an opaque facade
`workflow` bucket.

### Command

```powershell
cargo bench --bench equilibrium_legacy_callbacks `
  --no-default-features `
  --features equilibrium-instrumentation `
  -- --noplot
```

`cargo bench` uses Cargo's optimized benchmark profile. The benchmark selects
exactly one `Legacy(LM)`, `Legacy(NR)`, or `Legacy(TR)` backend per case through
the public equilibrium facade. It prints two independent `tabled` tables:

1. `accepted callback breakdown`: callback counts, callback time, shares of
   backend wall time, iterations, and solver overhead for accepted attempts;
2. `complete attempt ledger`: every started attempt, including rejected and
   failed attempts. `wall-only` means that the selected backend did not expose
   structured callback metrics; failed Legacy attempts now do expose them.

### Question and hypothesis

Question: where does time go inside Legacy LM, NR, and TR for P,T, monolithic
P,H, and nested P,H, and does ElementInventory change the answer?

Hypothesis: small controls are dominated by solver overhead, while larger
universes make residual, Jacobian, and linear-system work material. Therefore
callback parallelism must not be selected from call counts alone.

### Coverage

The matrix covers public-facade molecular and ElementInventory controls over
3, 15, 18, and 40-component universes, with separate `P,H monolithic` and
`P,H nested` labels. The 18-component TP-1907 and 40-component RP-1311
fixtures are the large phase-aware controls. Cases where a deliberately
isolated Legacy backend cannot converge remain visible in the attempt ledger;
they are not silently replaced by the production cascade.

### Release evidence

Representative rows from the optimized dashboard are below. Percentages are
within the selected backend attempt and therefore add up to the reported
callback/overhead budget; they are not facade percentages.

| case | backend | wall ms | iterations | residual | Jacobian | linear | overhead |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `rp1311_40_pt_element_inventory` | Legacy LM | 0.178 | 10 | 29.2% | 7.3% | 20.2% | 43.3% |
| `rp1311_40_pt_element_inventory` | Legacy NR | 0.146 | 11 | 36.3% | 10.3% | 33.6% | 19.9% |
| `rp1311_40_pt_element_inventory` | Legacy TR | 0.683 | 42 | 24.6% | 6.6% | 25.6% | 43.2% |
| `stanjan15_ph_element_inventory` | Legacy LM | 6.530 | 35 | 13.7% | 6.4% | 0.0% | 79.9% |
| `stanjan15_ph_element_inventory` | Legacy TR | 0.965 | 14 | 24.7% | 14.0% | 0.2% | 61.1% |
| `tp1907_18_ph_element_inventory` | Legacy TR, nested | 0.167 | 32 | 1.2% | 0.6% | 19.8% | 78.4% |

The full command prints all matrix rows. Failed Legacy attempts retain their
callback snapshot even when the nonlinear call returns an error; an
unavailable backend-specific split remains explicitly marked `wall-only`.
The first capture also showed
intentional isolated-backend failures: LM failed on one TP-1907 nested P,H
case, NR failed on TP-1907 and both Stanjan P,H cases in this seed setup, and
TR rejected one small P,H candidate and failed the TP-1907 P,T case. This is
useful solver-selection evidence, not a correctness regression: production
policy is not changed by this benchmark.

### Correctness gate

```powershell
cargo test --lib --no-default-features equilibrium_nonlinear -- --quiet
cargo test --lib --no-default-features equilibrium_rst_backend_tests -- --quiet
cargo test --lib --no-default-features equilibrium_backend_adapter -- --quiet
cargo check --bench equilibrium_legacy_callbacks --no-default-features `
  --features equilibrium-instrumentation --quiet
```

Observed gate:

```text
equilibrium_nonlinear: 36 passed
equilibrium_rst_backend_tests: 6 passed
equilibrium_backend_adapter: 7 passed
tp1907_upper_inventory_is_backend_invariant: 1 passed
production_default_retains_a_scale_robust_backend_on_real_tp1907_story: 1 passed
failed_legacy_attempt_preserves_callback_metrics: 1 passed
cargo check equilibrium_legacy_callbacks: passed
full Criterion matrix (36 case/backend functions): passed
```

### Decision and next action

`PERF-LEGACY-01` is accepted as measurement infrastructure. The immediate
optimization target is to split the currently unaccounted solver remainder in
large P,H and TR cases before calling it orchestration overhead. It includes
iteration bookkeeping, norms, trial-vector construction, feasibility and
line-search/trust-region decisions, matrix products or copies outside the
timed linear solve, and any other Legacy solver work not covered by the three
callback timers. Residual, Jacobian, and linear-system work is material in the
40-component cases, but the data does not yet justify parallel callbacks
because callback shares vary strongly by backend and route. No equations,
solver order, or production default changed. Failed Legacy callback metrics
are now part of the internal attempt ledger. Allocation counters and a
before/after optimization capture remain open TODO items.

## 88. PERF-LEGACY-02: Remove diagnostics from the ordinary Legacy path

### Story name

As a production user, I need Legacy solver callbacks and linear solves to run
without diagnostic-scope checks when instrumentation is not compiled in, so
the optional benchmark tooling does not become ordinary solver orchestration
overhead.

### Commands

```powershell
cargo check --lib --no-default-features --quiet
cargo check --bench equilibrium_legacy_callbacks --no-default-features --quiet
cargo check --bench equilibrium_legacy_callbacks --no-default-features `
  --features equilibrium-instrumentation --quiet
cargo test --lib --no-default-features equilibrium_nonlinear -- --quiet
cargo test --lib --no-default-features equilibrium_backend_adapter -- --quiet
```

### Question and hypothesis

Question: does the normal library build pay for the Legacy diagnostic timer
scope on every residual, Jacobian, and linear-system call?

Hypothesis: compile-time selection can make the disabled path a direct callback
or solve invocation. Instrumented builds retain the same per-attempt metrics,
while ordinary builds contain no thread-local diagnostic lookup or `Instant`
creation in these hot calls.

### Change under test

`timed_legacy_residual`, `timed_legacy_jacobian`, and the linear-system helper
now have a direct `#[cfg(not(...))]` implementation. The timing implementation
and the diagnostic counter functions exist only for tests or the
`equilibrium-instrumentation` feature. No equations, tolerances, iteration
order, acceptance rule, or backend policy changed.

The same Legacy pass retains only ownership-preserving copy elision already
visible in the solver code: LM uses the normal-equation product as its mutable
system, and LM/NR/TR move owned residual buffers where the value is no longer
needed. The mixed Criterion deltas do not justify claiming a general solver
speedup, so these changes remain correctness-preserving candidates rather
than a new performance baseline.

### Correctness and build evidence

```text
equilibrium_nonlinear: 36 passed
equilibrium_backend_adapter: 7 passed
cargo check --lib --no-default-features: passed
cargo check equilibrium_legacy_callbacks without instrumentation: passed
cargo check equilibrium_legacy_callbacks with instrumentation: passed
```

The instrumented release callback benchmark still publishes the same accepted
callback table and complete attempt ledger. Its repeated TR measurement did
not establish a statistically significant wall-time change, so this story
does not claim a numeric speedup for the diagnostic build. A separate attempt
to precompute J^T J and remove the TR Jacobian clone was rejected after a
`1.655 ms` result versus the previous `1.590 ms`-class measurement; that
experiment was reverted.

The matching release run without `equilibrium-instrumentation` completed the
same `rp1311_40_pt_element_inventory_legacy_TR` case at `1.533 ms`. The
preceding instrumented run was `1.643 ms`, a `-6.98%` on/off difference in
this large TR control by the Criterion comparison. The displayed medians
themselves imply about `-6.7%`; the distinction is the stored Criterion
baseline denominator. This measures the cost of optional diagnostics, not an
algorithmic TR speedup; a clean historical no-feature baseline is still
needed before generalizing the percentage.

### Decision and next action

Accept the compile-time cleanup as a production-path optimization and keep
instrumentation opt-in. Do not label the broad `solver_overhead` remainder as
orchestration yet: the next Legacy optimization must add a measured split for
matrix products, norms, trial-vector construction, feasibility, and
line-search/trust-region work, then target only a material bucket. The
instrumented benchmark remains the regression gate for preserving callback
diagnostics.

## 89. PERF-LEGACY-03: Legacy P,H end-to-end optimization gate

### Story name

Before changing Legacy P,H orchestration, measure the complete chain from the
public facade to the accepted Legacy backend. The decision must distinguish
facade cost, P,H orchestration, trial work, nonlinear solve time, and the
backend-attempt wall clock instead of treating the whole `workflow` bucket as
solver time.

### Command

```powershell
$env:CARGO_PROFILE_BENCH_LTO='false'
$env:CARGO_PROFILE_BENCH_CODEGEN_UNITS='16'
cargo bench --bench equilibrium_baseline `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  -- --noplot
```

`cargo bench` uses Cargo's optimized `bench` profile. The table below is the
fixed-repeat diagnostic capture from the same release-style run; Criterion's
statistical output is intentionally not used as the stage attribution.

### Question / hypothesis

Question: is Legacy orchestration large enough to justify a dedicated
optimization, or is the nonlinear backend itself the material part of Legacy
P,H wall time?

Hypothesis: orchestration should be deferred when nonlinear work is only a few
percent of facade time. Large P,H cases must be measured through the complete
chain before applying that rule, because nested and monolithic routes have
different timing scopes.

### What the columns mean

- `profile` identifies the execution graph: `P,T accepted`, `P,H monolithic`,
  `P,H nested`, or `P,H recovery`. Rows from different profiles are not a
  single homogeneous population and must not be ranked by percentage alone.
- `facade total` is the complete public-facade wall measurement.
- `P,H workflow` is the outer P,H workflow timing, shown as a share of facade
  total.
- `trial total` is the accepted nested scalar-trial aggregate, shown as a
  share of P,H workflow. It is intentionally `n/a (monolithic)` for the
  monolithic route.
- `nonlinear solve` is the nested inner nonlinear aggregate as a share of
  trial total. For monolithic P,H it is the measured monolithic backend
  execution as a share of P,H workflow.
- `Legacy backend wall` is the accepted Legacy backend metric wall time. It is
  printed without a percentage because nested accepted backend metrics may
  aggregate several trial/attempt scopes and therefore do not share one
  denominator with the outer P,H timing.

The percentages are parent-relative. They are not additive across columns:
`P,H workflow` contains `trial total`, and `trial total` contains the nested
nonlinear solve. `Legacy backend wall` is a separate backend-attempt
observation, not another segment to add to the chain. Only the macro table's
`facade outside workflow + solve transaction` partition is additive. The
workflow and diagnostic detail tables are inclusive views; their `profile` and
`scope` columns are the comparison key.

### Release evidence

| case | backend | P,H route | profile | facade total | P,H workflow | trial total | nonlinear solve | Legacy backend wall |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `ph_element_inventory_point_legacy_lm` | Legacy LM | monolithic | P,H monolithic | 0.28 ms | 0.19 ms / 66.7% | n/a | 0.13 ms / 69.2% | 0.13 ms |
| `stanjan15_ph_element_inventory_point_legacy_lm` | Legacy LM | monolithic | P,H monolithic | 7.95 ms | 7.31 ms / 91.9% | n/a | 6.97 ms / 95.4% | 6.97 ms |
| `stanjan15_ph_element_inventory_point_legacy_tr` | Legacy TR | monolithic | P,H monolithic | 1.99 ms | 1.38 ms / 69.1% | n/a | 1.05 ms / 76.1% | 1.04 ms |
| `tp1907_18_ph_element_inventory_point_legacy_tr` | Legacy TR | nested | P,H nested | 3.86 ms | 3.26 ms / 84.4% | 3.12 ms / 95.8% | 2.13 ms / 68.2% | 0.46 ms |
| `rp1311_40_ph_element_inventory_point_legacy_lm` | Legacy LM | nested | P,H nested | 12.34 ms | 10.72 ms / 86.9% | 10.25 ms / 95.6% | 6.12 ms / 59.7% | 3.48 ms |
| `rp1311_40_ph_element_inventory_point_legacy_tr` | Legacy TR | nested | P,H recovery | 16.55 ms | 15.56 ms / 94.0% | 1.44 ms / 9.2% | 0.12 ms / 8.4% | 0.26 ms |

### Result and decision

The measurement rejects a global statement that Legacy orchestration is the
next optimization target. In the small control it is not material to the
complete call. In STANJAN-15 monolithic P,H, Legacy LM spends about `6.97 ms`
of `7.95 ms` facade time in the monolithic nonlinear execution, so
orchestration is not the dominant opportunity. Legacy TR has a smaller
nonlinear share, but its total is also much lower.

For nested TP-1907/18 and RP-1311/40, the accepted trial chain is visible and
the nonlinear share is substantial for LM (`6.12 ms` of `10.25 ms` trial),
while Legacy TR's accepted backend wall is small. The previously observed
`15.81 ms` outside-workflow remainder was a reporting defect: the normalized
recovery result retained the inner recovery timing but not the complete outer
P,H transaction wall. The corrected RP-1311/40 Legacy TR capture attributes
`15.56 ms` to the complete `P,H recovery` workflow and leaves only `0.99 ms`
outside the workflow for facade-owned work. The recovery parent and its
physical-failure, normalized-discovery, and physical-retry attempt rows are
now printed separately.

Therefore Legacy orchestration optimization is deferred until a material
bucket is demonstrated inside a matching profile and parent scope. The next
valid Legacy experiment should use the attempt ledger and solver callback
table to split the nonlinear remainder into concrete operations. A change is
worthwhile only when its owning bucket is a material fraction of the same
complete-facade case and all correctness gates remain green.

### Correctness and coverage gate

```powershell
cargo check --bench equilibrium_baseline --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" --quiet
cargo test --lib --no-default-features equilibrium_ph_workflow -- --quiet
cargo test --lib --no-default-features equilibrium_nonlinear -- --quiet
cargo test --lib --no-default-features equilibrium_backend_adapter -- --quiet
cargo test --lib --no-default-features nasa_tp1907_chon_graphite_tests -- --quiet
```

The benchmark's pretty-coverage assertions require every baseline scenario to
have macro, facade, nonlinear, attempt, and accepted-solver rows. The Legacy
P,H envelope additionally requires one row for every non-production Legacy
P,H case, including both monolithic and nested routes.

## 90. PERF-DASHBOARD-01: Heterogeneous stage evidence is explicitly scoped

### Story name

The equilibrium performance dashboard must make every large wall-time bucket
actionable without pretending that different execution graphs have the same
denominator. Shared benchmark fixtures must also have one source of truth.

### Command

```powershell
cargo build --release --bench equilibrium_baseline `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint"
$exe = Get-ChildItem target\release\deps\equilibrium_baseline-*.exe |
  Sort-Object LastWriteTime -Descending | Select-Object -First 1
& $exe.FullName --noplot
```

### Question / hypothesis

Question: can a reader identify the owner of a large facade wall time and
compare Legacy with RST without adding overlapping parent and child rows?

Hypothesis: a macro partition plus separate scoped detail tables is sufficient
if every row carries its execution `profile` and its parent `scope`.

### Output contract

- **Macro time budget:** only `facade outside workflow + solve transaction`
  is additive and must equal 100% of facade wall time.
- **Facade detail:** resolution, request preparation, route selection,
  normalization, phase control, thermochemistry, and result assembly.
- **P,H workflow detail:** scalar orchestration, enthalpy evaluation, trial
  preparation/execution, monolithic preparation/execution, validation,
  publication, and the complete recovery transaction where applicable.
- **Inner/recovery detail:** fixed-P,T stages or the explicitly named P,H
  recovery transaction; recovery is not mislabeled as fixed-P,T.
- **Nonlinear envelope:** parent wall, all attempt wall, known preparation,
  runner/orchestration remainder, backend coverage, missing metrics, and
  missing wall time.
- **Attempt ledger:** every accepted, failed, rejected, fallback, and
  normalization-recovery attempt, including Legacy callback diagnostics.
- **Accepted solver table:** residual, Jacobian, linear-system, and solver
  overhead timings for the accepted backend only.

Every table includes `profile`: `P,T accepted`, `P,H monolithic`, `P,H
nested`, or `P,H recovery`. Percentages are valid only relative to the row's
own parent `scope` and profile. Backend wall is a separate observation and is
never added to facade or workflow wall. This prevents an RST preparation
parent, a Legacy callback attempt, and a nested scalar trial from being read as
three additive slices of one clock.

### Representative output

The corrected 40-component Legacy TR recovery point printed:

```text
profile=P,H recovery
facade total       16.55 ms
P,H workflow       15.56 ms / 94.0%
recovery parent    15.56 ms / 100.0% of P,H workflow
facade outside      0.99 ms / 6.0%
nonlinear solve     0.12 ms / 8.4% of accepted nested trial
Legacy backend      0.26 ms (separate attempt aggregate)
```

The prior `15.81 ms` outside-workflow value no longer appears. It was caused
by a normalized recovery result carrying inner recovery timing while its outer
P,H timing still described only the normalized solve. The facade now restores
the complete physical transaction wall; the remaining outside-workflow time is
ordinary facade-owned resolution/request/thermochemistry work.

### Result

The dashboard is decision-capable for optimization: additive macro budget,
inclusive detail, backend-attempt ledger, and callback diagnostics are
separate views with explicit comparison keys. Common frozen species lists,
pressure constants, phase-component construction, and the TP-1907 elemental
seed are shared through `benches/support/mod.rs`; case-specific builders and
timing contracts remain local to each benchmark target.

## 91. PERF-LEGACY-04: Legacy solver remainder is decomposed by operation

### Story name

The retained Legacy LM, NR, and TR backends must expose enough per-attempt
timing evidence to decide whether their `solver overhead` is worth optimizing.
The benchmark must preserve failed attempts and nested P,H attempts, because a
single accepted row cannot explain retry cost.

### Command

```powershell
cargo bench --bench equilibrium_legacy_overhead `
  --no-default-features `
  --features equilibrium-instrumentation `
  -- --noplot
```

This is a release-style Criterion benchmark. The dashboard is printed once
before Criterion repetitions. For a quicker local rebuild, the same target
can be built with LTO disabled, but before/after comparisons must use the same
profile settings.

### Question / hypothesis

Question: on large Legacy calculations, is the remaining backend time caused
by matrix algebra, norms, trial-vector construction, feasibility, globalization
logic, or an unmeasured remainder? Is any one named bucket large enough to
justify an optimization?

Hypothesis: small cases are dominated by fixed solver overhead, while large
40-component cases should reveal a material matrix/linear or residual backend
bucket. Norms, trial-vector construction, feasibility, and globalization are
not optimization targets unless their measured parent-relative share is
material.

### Output contract

The benchmark prints two intentionally different tables.

- **Macro facade/backend comparison:** `facade_ms`, accepted Legacy attempt
  wall and its facade share, all-attempt wall and its facade share, accepted
  overhead share, structured detail total, and the accepted unattributed
  remainder. The accepted row is one attempt; `all_attempts_ms` includes every
  measured attempt in the route.
- **Internal Legacy attempt detail:** one row per measured accepted, failed, or
  rejected attempt with callback, linear system, matrix operations, norms,
  trial-vector, feasibility, globalization, and unattributed milliseconds.
  Percentages are relative to that attempt's backend wall. P,H nested trial
  and recovery attempts remain separate rows.

The internal columns are mutually exclusive by construction: callback and
linear time come from the existing callback instrument, the five new fields
come from the Legacy solver stage instrument, and `unattributed` is the
remaining `solver_overhead` after subtracting those five fields. This table is
therefore suitable for optimization triage; it is not an inclusive workflow
table.

### Release evidence

Representative output from the focused release-style run:

| profile | backend | facade | accepted / facade | all attempts | all / facade | overhead share | material named buckets |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 40-component P,T ElementInventory | Legacy LM | 1.448 ms | 0.186 ms / 12.8% | 0.186 ms | 12.8% | 53.2% | matrix 0.030 ms; unattributed 0.069 ms |
| 40-component P,T ElementInventory | Legacy NR | 1.252 ms | 0.150 ms / 12.0% | 0.150 ms | 12.0% | 27.3% | linear 0.056 ms; unattributed 0.041 ms |
| 40-component P,T ElementInventory | Legacy TR | 1.871 ms | 0.699 ms / 37.4% | 0.699 ms | 37.4% | 48.5% | matrix 0.129 ms; linear 0.210 ms; unattributed 0.210 ms |
| 40-component nested P,H ElementInventory | Legacy LM | 3.846 ms | 0.192 ms / 5.0% | 0.542 ms | 14.1% | 50.5% | matrix 0.031 ms; unattributed 0.066 ms |
| 40-component nested P,H ElementInventory | Legacy NR | 3.316 ms | 0.141 ms / 4.3% | 0.419 ms | 12.6% | 31.9% | linear 0.050 ms; unattributed 0.045 ms |
| 40-component nested P,H ElementInventory | Legacy TR | 4.823 ms | 0.266 ms / 5.5% | 1.141 ms | 23.7% | 42.9% | matrix 0.049 ms; linear 0.079 ms; unattributed 0.065 ms |

In the same capture, norms, trial-vector construction, feasibility, and
globalization were below `0.001 ms` at display precision in the 40-component
accepted rows. The 15-component monolithic P,H rows are a separate profile:
feasibility measured about `0.624 ms` of `1.025 ms` accepted LM wall and
`1.339 ms` of `2.462 ms` accepted TR wall. Those figures are not pooled with
the 40-component elemental rows.

### Result and decision

The benchmark answers the previous `solver overhead` question: it is not one
single mysterious operation. On large elemental cases, linear systems and
matrix algebra are measurable, while a material remainder still exists after
those named stages. On the 15-component monolithic P,H route, feasibility is
the largest named non-callback bucket. The very small norm, trial-vector, and
globalization measurements do not justify standalone optimization.

No algorithmic Legacy optimization is accepted by this story. The next
experiment should target one matched 40-component matrix/linear/unattributed
bucket, retain the same correctness matrix, and compare accepted plus failed
attempt paths. The benchmark executable and the ordinary library build keep
the diagnostic capture opt-in.

### Correctness and coverage gate

```powershell
cargo check --bench equilibrium_legacy_overhead `
  --no-default-features `
  --features equilibrium-instrumentation --quiet
cargo test --lib --no-default-features equilibrium_nonlinear -- --quiet
cargo test --lib --no-default-features equilibrium_backend_adapter -- --quiet
cargo test --lib --no-default-features equilibrium_solver_policy -- --quiet
cargo test --lib --no-default-features nasa_tp1907_chon_graphite_tests -- --quiet
```

The benchmark's detail table is the required evidence gate for the next
Legacy solver optimization. A change is accepted only when correctness,
attempt classification, and the relevant profile's facade/backend wall all
remain stable enough for a meaningful before/after comparison.

## 92. PERF-LEGACY-05: Trust-region scalar model reduction avoids dense J^T J

### Story name

The Legacy Trust Region solver must calculate the scalar quadratic-model term
without materializing a dense `J^T J` matrix when the algorithm only needs
`p^T J^T J p`. The replacement must preserve the equation, step acceptance,
backend order, and recovery behavior.

### Command

```powershell
cargo rustc --release --bench equilibrium_legacy_overhead `
  --no-default-features `
  --features equilibrium-instrumentation `
  --config 'profile.release.lto="off"' `
  --config 'profile.release.codegen-units=16' `
  -- --noplot
```

The focused dashboard is run directly after this build for repeatable local
inspection. Before/after captures use the same release profile and machine.

### Question / hypothesis

Question: does the scalar-only Trust Region model-reduction expression spend
avoidable time constructing a dense `J^T J` matrix?

Hypothesis: replacing that construction with `Jp` followed by `(Jp)^T(Jp)`
will remove the matrix-operation bucket while preserving the accepted
solution. The wall-time win is expected to be case- and noise-dependent, so
the bucket collapse is stronger evidence than one dashboard run.

### Change under test

The production expression now evaluates
`p^T J^T J p = (Jp)^T(Jp)` directly. It does not change the Trust Region
linear solve, the model equation, the trust-radius update, or the solver
ordering. The stage instrument records the `Jp` product as matrix work and
the two scalar products as norm work.

### Release evidence

The matched 40-component Legacy TR observations were:

| profile | facade | accepted backend | matrix operations | linear system | unattributed solver remainder |
| --- | ---: | ---: | ---: | ---: | ---: |
| P,T ElementInventory, before | 1.871 ms | 0.699 ms | 0.129 ms | 0.210 ms | 0.210 ms |
| P,T ElementInventory, after run 1 | 2.017 ms | 0.631 ms | 0.001 ms | 0.196 ms | 0.224 ms |
| P,T ElementInventory, after run 2 | 1.780 ms | 0.627 ms | 0.001 ms | 0.196 ms | 0.238 ms |
| P,T ElementInventory, after run 3 | 1.714 ms | 0.589 ms | 0.001 ms | 0.196 ms | 0.237 ms |
| nested P,H ElementInventory, before | 4.823 ms | 0.266 ms | 0.049 ms | 0.079 ms | 0.065 ms |
| nested P,H ElementInventory, after run 1 | 5.171 ms | 0.237 ms | 0.001 ms | 0.076 ms | 0.074 ms |
| nested P,H ElementInventory, after run 2 | 4.593 ms | 0.232 ms | 0.001 ms | 0.076 ms | 0.080 ms |
| nested P,H ElementInventory, after run 3 | 4.644 ms | 0.223 ms | 0.001 ms | 0.076 ms | 0.082 ms |

The matrix bucket collapsed to display noise in both large profiles. The
after wall times suggest a useful improvement, but these dashboard runs are
not a statistically controlled Criterion before/after series; no exact speedup
percentage is claimed. The previous rejected experiment that cloned/reused a
full dense transformation remains rejected and is not part of this result.

### Correctness gate and decision

```powershell
cargo test --lib --no-default-features equilibrium_nonlinear -- --quiet
cargo test --lib --no-default-features nasa_tp1907_chon_graphite_tests -- --quiet
cargo test --lib --no-default-features nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together -- --quiet
cargo check --bench equilibrium_legacy_overhead `
  --no-default-features `
  --features equilibrium-instrumentation --quiet
```

The gates passed: 28 nonlinear tests, 22 TP-1907 frozen tests, and the
isolated TP-1906 extensivity test passed; the latter took 82.16 s in the debug
test profile because it solves seven physical scales. The narrow optimization
is accepted as correctness-clean and structurally justified. Broad Legacy
orchestration remains deferred until a controlled benchmark shows a material
facade-level gain; the next Legacy work must target a named material bucket,
not generic `solver overhead`.

## 93. PERF-LEGACY-06: Trust-region boundary Cauchy steps skip unused LU

### Story name

The Legacy Trust Region solver must not compute the Newton/LU step when the
Cauchy step already reaches the trust-region boundary. The dogleg result,
singular-Jacobian behavior when Newton data is required, and accepted solution
must remain correct.

### Question / hypothesis

Question: how much of the measured Legacy TR linear-system bucket is spent on
Newton data that is discarded because `||p_cauchy|| >= delta`?

Hypothesis: calculating the Cauchy step first and delaying LU until the full
Newton step or dogleg interpolation is needed will reduce linear-system wall
time on boundary iterations without changing the selected step.

### Change under test

The solver now computes Cauchy data before Newton data. When Cauchy already
lies on or beyond the trust boundary, it scales that step directly. Otherwise
it performs the existing LU solve and retains the existing dogleg selection.
The zero-Jacobian case still reaches LU and remains a singular-Jacobian error.

### Command

```powershell
cargo rustc --release --bench equilibrium_legacy_overhead `
  --no-default-features `
  --features equilibrium-instrumentation `
  --config "profile.release.lto='off'" `
  --config "profile.release.codegen-units=16"
```

The resulting executable is run with `--noplot`; before/after captures use
the same release profile and machine.

### Release evidence

| profile | facade | accepted TR wall | linear system | matrix | unattributed |
| --- | ---: | ---: | ---: | ---: | ---: |
| P,T ElementInventory, before | 1.871 ms | 0.699 ms | 0.210 ms | 0.129 ms | 0.210 ms |
| P,T ElementInventory, after A | 2.152 ms | 0.572 ms | 0.130 ms | 0.003 ms | 0.220 ms |
| P,T ElementInventory, after B | 1.961 ms | 0.566 ms | 0.125 ms | 0.003 ms | 0.219 ms |
| P,T ElementInventory, after C | 2.009 ms | 0.705 ms | 0.135 ms | 0.003 ms | 0.352 ms |
| nested P,H ElementInventory, before | 4.823 ms | 0.266 ms | 0.079 ms | 0.049 ms | 0.065 ms |
| nested P,H ElementInventory, after A | 5.241 ms | 0.227 ms | 0.060 ms | 0.003 ms | 0.080 ms |
| nested P,H ElementInventory, after B | 4.627 ms | 0.221 ms | 0.060 ms | 0.001 ms | 0.080 ms |
| nested P,H ElementInventory, after C | 4.699 ms | 0.224 ms | 0.061 ms | 0.005 ms | 0.077 ms |

The new ordering reduces the observed linear bucket in the large profiles,
while the dashboard remains noisy and the unattributed bucket is not reduced
by this change. These are focused dashboard observations, not a controlled
Criterion before/after series; no exact facade speedup is claimed.

### Correctness gate and decision

```powershell
cargo test --lib --no-default-features equilibrium_nonlinear -- --quiet
cargo test --lib --no-default-features nasa_tp1907_chon_graphite_tests -- --quiet
cargo test --lib --no-default-features nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together -- --quiet
cargo check --bench equilibrium_legacy_overhead `
  --no-default-features `
  --features equilibrium-instrumentation --quiet
```

The gates passed: 36 nonlinear tests, 22 TP-1907 frozen tests, and the
isolated TP-1906 extensivity test passed. The optimization is accepted as a
narrow TR linear-work reduction. Broad Legacy orchestration and parallel
callback work remain deferred; the next candidate must be selected from the
remaining detailed dashboard buckets.

## 94. PERF-LEGACY-07: Large Legacy P,H non-solver budget is decomposed

### Story name

The large Legacy P,H ElementInventory cases must expose every material stage
outside the nonlinear solver. A `workflow` or `recovery` parent row is useful
for partitioning, but it is not an optimization target until its child timing
and backend-attempt evidence are visible.

### Question / hypothesis

Question: for Legacy LM/NR/TR on a 40-component ElementInventory case, how much
time remains outside nonlinear solving, and is it ordinary nested orchestration
or extensive-normalization recovery?

Hypothesis: ordinary nested P,H LM/NR contains a measurable numeric-closure
cost, while TR's large P,H wall is caused by a recovery transaction and must be
analyzed as separate recovery children rather than pooled with normal trials.

### Command

Release dashboard, ten repetitions per scenario, with the public facade and
all diagnostic tables enabled:

```powershell
cargo rustc --release --bench equilibrium_baseline `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  --config "profile.release.lto='off'" `
  --config "profile.release.codegen-units=16"

target\release\deps\equilibrium_baseline-<fresh-hash>.exe --noplot
```

### Release evidence

The following is a fresh release capture. The macro table is additive only in
the facade partition; nested child rows are inclusive observations and must
not be summed.

| case | facade | facade outside workflow | P,H workflow | trial / recovery child | inner nonlinear | inner outside nonlinear | named non-solver observation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| RP-1311/40 P,H Legacy LM | 13.78 ms | 1.79 ms | 11.99 ms | trial 11.52 ms | 7.32 ms | 2.87 ms | numeric closure 2.21 ms |
| RP-1311/40 P,H Legacy NR | 10.81 ms | 1.74 ms | 9.07 ms | trial 8.60 ms | 4.47 ms | 2.81 ms | numeric closure 2.17 ms |
| RP-1311/40 P,H Legacy TR | 15.74 ms | 1.15 ms | 14.59 ms | recovery 14.59 ms | 0.12 ms* | 14.47 ms* | recovery children below |

The TR recovery transaction is decomposed by its own complete recovery wall:

| recovery child | time | share |
| --- | ---: | ---: |
| physical failure | 1.89 ms | 12.9% |
| normalization | 0.11 ms | 0.8% |
| normalized discovery | 10.75 ms | 73.6% |
| physical retry | 1.63 ms | 11.2% |
| reconstruction | 0.00 ms | 0.0% |
| recovery remainder | 0.21 ms | 1.5% |

`recovery remainder` now contains only control work not assigned to a named
child. The initial failed physical attempt has its own wall-time row and its
backend-attempt table remains the authoritative place to inspect individual
failed calls. Recovery child rows include their own solver work; therefore the
starred TR complement is a conservative published nonlinear-envelope
complement, not a claim that all normalized-discovery wall is non-solver work.

The RP-1311/40 P,T LM/NR/TR controls are present as well (`1.35`, `1.30`, and
`1.58 ms` in this capture). They are below the project materiality gate of
`10 ms` and are retained for coverage, not optimization priority.

### Decision

Do not start a generic Legacy-orchestration refactor from the macro percentage.
For the only ordinary large cases above the gate, numeric closure is visible at
`2.21/2.17 ms` for LM/NR but remains below the current `10 ms` materiality
threshold. For TR, the failed physical attempt is now separated and the true
unassigned recovery remainder is only `0.21 ms`. No Legacy orchestration work
is justified at this stage; facade resolution, request preparation,
validation, and every other sub-`10 ms` aggregate remain out of scope.

The benchmark now exposes the required macro, facade, P,H, inner/recovery,
nonlinear-envelope, backend-attempt, and accepted-callback tables for every
Legacy scenario in the matrix. This story is a measurement gate, not an
algorithmic optimization acceptance.

## 95. PERF-LEGACY-08/09: Numeric-closure clone elision and solver trial-buffer probe

### Story name

The numeric closure must not deep-copy a large thermochemistry payload twice
when the phase bridge already owns a private working copy. A Legacy solver
allocation change is worthwhile only when trial-vector construction is a
material part of the measured nonlinear wall.

### Question / hypothesis

Question 1: does the phase bridge call the historical public Gibbs-closure
factory in a way that repeats a large `SubsData` clone?

Hypothesis 1: a crate-private in-place entry point, used only after the bridge
has created its private working copy, removes the redundant clone while the
public copy-on-entry contract remains unchanged.

Question 2: can LM/NR/TR reuse one `Vec<f64>` for trial states without changing
accepted candidates or line-search behavior?

Hypothesis 2: allocation reuse might reduce Legacy solver overhead on the
40-component profiles, but it must be rejected if the stage remains below
measurement resolution or wall-time captures do not improve.

### Change under test

`PERF-LEGACY-08` factors the Gibbs-closure loop into a private local helper.
`SubsData::calculate_dG0_fun_one_phase()` keeps its historical copy-on-entry
behavior; the phase bridge calls a crate-private
`calculate_dG0_fun_one_phase_in_place()` only on its already-cloned working
payload. No equations, callback contract, solver order, or accepted-state
rules changed.

`PERF-LEGACY-09` temporarily reused a trial `Vec<f64>` in Legacy LM, NR, and
TR. The probe was reverted after the detailed dashboard showed trial-vector
work at `0.000 ms` for the large accepted rows and repeated facade/backend
captures moved within ordinary run noise.

### Command

```powershell
cargo rustc --release --bench equilibrium_baseline `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  --config "profile.release.lto='off'" `
  --config "profile.release.codegen-units=16"

target\release\deps\equilibrium_baseline-<fresh-hash>.exe --noplot
```

The Legacy detail harness was also used for the trial-buffer decision:

```powershell
cargo bench --bench equilibrium_legacy_overhead `
  --no-default-features `
  --features equilibrium-instrumentation -- --noplot
```

### Release evidence

Three fresh release captures after `PERF-LEGACY-08` gave the following
RP-1311/40 nested P,H observations. Numeric-closure timing is an inclusive
stage row; it must not be added to the nonlinear row or to another inclusive
stage. Captures 1-3 were made while the temporary `PERF-LEGACY-09` probe was
present; that probe does not touch closure construction.

| capture | LM numeric closure | NR numeric closure | LM backend wall | NR backend wall |
| --- | ---: | ---: | ---: | ---: |
| 1 | 1.06 ms | 0.93 ms | 10.14 ms | 7.16 ms |
| 2 | 1.03 ms | 0.94 ms | 10.15 ms | 7.20 ms |
| 3 | 1.12 ms | 1.06 ms | 11.00 ms | 7.84 ms |
| median | 1.06 ms | 0.94 ms | 10.15 ms | 7.20 ms |

After removing the rejected trial-buffer probe, a current release smoke run
reported `0.81 ms` numeric closure for both LM and NR, with facade walls
`12.12 ms` and `9.09 ms` respectively and `dashboard_exit=0`. This is a
current-tree consistency check, not an additional before/after sample.

The previous retained Story 94 baseline named numeric closure at `2.21 ms`
for LM and `2.17 ms` for NR. This supports removal of a repeated large
payload clone and a substantial reduction of that named preparation stage;
the release captures are not a controlled Criterion before/after series, so
they do not claim the same percentage reduction for complete facade wall.

The detailed Legacy harness measured trial-vector construction at `0.000 ms`
at display precision in the large accepted rows. The corresponding LM/NR/TR
facade and backend walls varied materially more than that stage across
captures. The trial-buffer probe therefore has no accepted performance claim
and its code was removed.

### Correctness gate and decision

```powershell
cargo test --lib --no-default-features equilibrium_nonlinear -- --quiet
cargo test --lib --no-default-features equilibrium_backend_adapter -- --quiet
cargo test --lib --no-default-features phase_equilibrium_problem_tests -- --quiet
cargo test --lib --no-default-features equilibrium_calculator -- --quiet
cargo test --lib --no-default-features nasa_tp1907_chon_graphite_tests -- --quiet
cargo test --lib --no-default-features nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together -- --quiet
```

The relevant gates passed: nonlinear `36/36`, backend adapter `7/7`, phase
bridge `35/35`, calculator `28/28`, frozen TP-1907 `22 passed + 6 ignored`,
and the TP-1906 extensive scale test `1/1`. A complete lib run also exposed
three unrelated pre-existing Kinetics failures; none is in the changed
ChemEquilibrium path.

`PERF-LEGACY-08` is accepted narrowly as large-payload clone elision. The
public API remains protected by its historical copy semantics. `PERF-LEGACY-09`
is rejected as below the measurement/materiality threshold. The next useful
nonlinear experiment is a shared numeric evaluation workspace for residual
and Jacobian state, with matched allocation/callback evidence; callback
parallelism remains deferred until callback work itself is a material budget.

## 96. PERF-PH-12/13: Changing-temperature P,H Gibbs ownership and trial snapshot reuse

### Story name

Changing-temperature P,H evaluation must not clone an owned thermochemistry
payload at the Gibbs refresh boundary, and a feasible trial must not be fully
evaluated again when the nonlinear backend immediately requests its residual.

### Command

Release preparation benchmark with optional counters and allocation evidence:

```powershell
cargo rustc --release --bench equilibrium_ph_preparation `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-allocation" `
  --config "profile.release.lto='off'" `
  --config "profile.release.codegen-units=16"

target\release\deps\equilibrium_ph_preparation-<fresh-hash>.exe --noplot
```

Correctness gates:

```powershell
cargo test --lib --no-default-features equilibrium_ph_thermochemistry -- --quiet
cargo test --lib --no-default-features equilibrium_ph_formulation -- --quiet
cargo test --lib --no-default-features equilibrium_ph_monolithic -- --quiet
cargo test --lib --no-default-features equilibrium_nonlinear -- --quiet
cargo test --lib --no-default-features phase_equilibrium_problem_tests -- --quiet
cargo test --lib --no-default-features nasa_tp1907_chon_graphite_tests -- --quiet
cargo test --lib --no-default-features nasa_tp1906_chon_graphite_ph_tests -- --quiet
```

### What we check

Hypothesis 1: `PhaseGibbsSource::evaluate` already owns its private
`SubsData`; calling the public `calculate_dG0_fun_one_phase()` adds a second
copy-on-entry. The crate-private in-place factory should preserve values and
errors while removing that redundant large-payload copy. The public API must
retain its historical copy-on-entry behavior.

Hypothesis 2: `is_feasible()` and `residual()` commonly receive the same exact
trial vector. Reusing the complete `PhCandidateSnapshot` should remove the
second Gibbs/enthalpy/mole/residual evaluation without accepting invalid
trials, publishing rejected state, or memoizing a different vector.

### Implementation

`PhaseGibbsSource::evaluate` now calls
`calculate_dG0_fun_one_phase_in_place()` on its mutex-protected, source-owned
payload. The public `SubsData::calculate_dG0_fun_one_phase()` remains
copy-on-entry for external callers.

`PreparedMonolithicPhRunner` now keeps one local exact-vector cache. A feasible
snapshot is stored only after the existing full candidate evaluation and mole
checks pass. The following residual consumes it only when every coordinate is
exactly equal. A failed feasibility evaluation, an invalid trial, a different
vector, or any consumed entry falls back to the original evaluation path.

### Release changing-target evidence

The benchmark prepares a session, solves the declared first target, then
solves a different second total-enthalpy target from that accepted state. It
also prints allocation mode, point wall, exact candidate snapshot counts, and
target-local allocated bytes. One current release capture:

| case | preparation | first point | second point | first snapshots | second snapshots | first allocated bytes | second allocated bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| small3 molecular | 10.187 ms | 1.511 ms | 0.679 ms | 3 | 2 | 122720 | 95886 |
| STANJAN-15 molecular | 377.227 ms | 35.710 ms | 5.228 ms | 32 | 20 | 3696563 | 1515919 |
| STANJAN-15 elemental | 376.579 ms | 7.920 ms | 5.265 ms | 2 | 20 | 930799 | 1517445 |
| RP-1311/40 elemental | 4716.048 ms | 51.761 ms | 9.913 ms | 2 | 3 | 7764103 | 1994519 |

The second target is genuinely different, but it is not a refresh-only
microbenchmark: prepared symbolic state and backend attempts can dominate it.
The `0.000 ms` second Gibbs-snapshot rows in the detailed table mean that the
prepared formulation reused its existing Gibbs capability; they do not prove
that every numeric `PhaseGibbsSource` refresh is free. A matched old/new
refresh-only capture is still needed for a percentage claim about PERF-PH-12.

The snapshot counter is the useful structural evidence for PERF-PH-13. It is
reported beside the changing-target wall and allocation numbers so a matched
before/after run can distinguish fewer full chemical evaluations from solver
iteration noise. This single capture is not a claim that the whole point wall
has decreased monotonically: first/second wall values move with backend route,
iterations, and basin selection.

### Matched PERF-PH-13 comparison

The same release profile, cases, public facade, prepared-session operation, and
three process repetitions were run with `trial_cache_mode=disabled-baseline`
and `trial_cache_mode=enabled`. The baseline was selected only through the
internal `equilibrium-benchmark-no-trial-cache` feature; production builds use
the enabled path.

| case | baseline second point median | enabled second point median | enabled wall delta | baseline snapshots | enabled snapshots | enabled allocation delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| small3 molecular | 0.592 ms | 0.650 ms | +9.8% | 2 | 2 | 0.0% |
| STANJAN-15 molecular | 5.244 ms | 5.515 ms | +5.2% | 26 | 20 | -4.8% |
| STANJAN-15 elemental | 5.273 ms | 5.100 ms | -3.3% | 26 | 20 | -4.8% |
| RP-1311/40 elemental | 10.347 ms | 10.037 ms | -3.0% | 3 | 3 | 0.0% |

The wall result is mixed and remains within ordinary solver/run variation; it
does not justify a percentage speedup claim. The structural result is clear
for the STANJAN cases: six complete candidate snapshots and about `4.8%` of
target-local allocation bytes disappear. RP-1311/40 does not exercise a
material feasibility-to-residual reuse opportunity in this prepared route.
Keep PERF-PH-13 as a narrow repeated-work reduction, and do not add a more
general cache or optimize its small key vector without a new workload showing
material benefit.

### Correctness result

All listed focused gates passed: thermochemistry `3/3`, formulation `12/12`,
monolithic runner `5/5`, nonlinear `36/36`, phase bridge `35/35`, frozen
TP-1907 `22 passed + 6 ignored`, and TP-1906 P,H frozen tests `22 passed + 6
ignored`. The prepared-session benchmark completed every listed case and
Criterion operation successfully.

### Decision

Accept PERF-PH-12 as a narrow ownership correction and PERF-PH-13 as a
localized full-snapshot reuse. Neither change changes equations, solver order,
acceptance tolerances, rollback, or the public copy-on-entry contract. The
matched series confirms structural work reduction but no uniform wall-time
speedup, so no facade percentage is reported. Further work should target a
refresh-only P,H measurement and interval/failure gates, not a broader cache.

## Story 97: PERF-PH-14 native-interval Gibbs closure reuse

### Command

```text
cargo test --lib --no-default-features tp1907_gibbs_reuses_native_intervals_between_temperature_points -- --nocapture
cargo test --lib --no-default-features nasa_tp1906_chon_graphite_ph_tests -- --quiet --test-threads=1
cargo check --lib --no-default-features
```

### Question / hypothesis

Does the private resolved P,H Gibbs source rebuild the complete set of callable
closures only when the selected native coefficient interval changes, while
preserving exact Gibbs values and the existing failure boundary? The useful
performance unit here is a full closure rebuild, not the cost of copying a
small temperature value.

### What changed

`PhaseGibbsSource` now stores the actual selected `(T_min, T_max)` interval for
each component after a successful refresh. A new temperature reuses all
closures only when it is contained in every stored interval. A crossing or an
out-of-range temperature follows the existing `extract_all_thermal_coeffs`
and in-place Gibbs factory path. The `ThermoCalculator` capability API and
`SubsData` gained a read-only current-interval query so the cache does not
guess from a fixture breakpoint or from a rounded temperature. The public
`from_functions` capability path is unchanged.

The common optional instrumentation gained
`gibbs_interval_reuses()`, separate from `gibbs_coefficient_refreshes()`.
The counter is disabled with the rest of equilibrium instrumentation in an
ordinary library build.

### Result

The focused TP-1907 frozen probe evaluated one complete Gibbs vector at
`700.0 K`, then at `700.5 K`, then in reverse at `699.5 K`. All returned
values were finite. The result was:

```text
refreshes=2  reuses=4
```

The two refreshes are the independent `gas` and `graphite` phase sources on
the first vector. Each later vector reuses one closure set per phase, hence
four reuse events. This is an assertion over the observed counters, not a
wall-time claim.

The sequential frozen P,H module completed `7 passed; 0 failed; 6 ignored`.
The library compile check passed. Existing warnings are unrelated deprecation
and unused-code warnings.

### Decision and remaining gates

Accept the narrow interval-aware reuse as a production change: it removes
complete closure construction from repeated in-interval temperatures without
changing equations, solver order, public ownership, or acceptance behavior.
Do not generalize it into a multi-temperature closure cache yet. Add the
remaining correctness matrix before claiming full boundary coverage:
crossing in both directions, components with different native intervals,
out-of-range and invalid-source errors, recovery after a failed refresh, and
agreement with the canonical P,T Gibbs/H/Cp values and provenance.
### Story 98: PERF-PH-14 native Gibbs interval reuse pays for itself in a matched release probe

**Name:** `PERF-PH-14 native Gibbs interval reuse`

**Command:**

```text
cargo bench --bench equilibrium_gibbs_interval --no-default-features --features "equilibrium-instrumentation" -- --noplot
cargo bench --bench equilibrium_gibbs_interval --no-default-features --features "equilibrium-instrumentation equilibrium-benchmark-no-gibbs-interval-reuse" -- --noplot
```

**Question / hypothesis:** When one resolved TP-1907 thermochemistry bundle is
evaluated at many temperatures that remain inside the selected native NASA
intervals, can the calculator reuse the already built Gibbs closures? The
control build forces the pre-change behavior: refresh closures for every new
temperature. The setup uses the public `EquilibriumCalculator` facade and a
frozen TP-1907 gas-only molecular feed; P,T basin discovery is outside the
timed sequence.

**What is measured:** 384 Gibbs vectors, in the sequence
`[700.0, 700.5, 699.5] K` repeated 128 times. The table reports wall time for
the sequence, mean time per vector, and instrumentation evidence for closure
refreshes and interval reuses. This is a focused thermochemistry benchmark,
not a P,H or nonlinear-solver benchmark.

**Result:** Five direct runs of each optimized binary produced the following
range and mean:

| mode | vectors | wall range | mean wall | mean/vector | refreshes | interval reuses |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| native interval reuse | 384 | 0.291--0.345 ms | 0.307 ms | 0.80 us | 2 | 766 |
| refresh every temperature (control) | 384 | 6.492--7.315 ms | 6.805 ms | 17.72 us | 768 | 0 |

The matched local speedup is approximately `22.2x`, or `95.5%` lower wall
time for this Gibbs-evaluation sequence. The counters prove that the current
path builds one closure set per phase on the first vector and reuses both
phase sets for every subsequent temperature in the interval.

**Decision:** Accept `PERF-PH-14` as a measured local optimization. Do not
extrapolate the percentage to complete P,H wall time: nonlinear solve and
other P,H stages are intentionally absent from this probe. Remaining
correctness gates are interval crossing in both directions, heterogeneous
component intervals, out-of-range and invalid-source errors, and failed
refresh recovery. The benchmark is ready to serve as a regression baseline
for later P,H changing-temperature measurements.

### Story 99: PERF-PH-15 production P,H path reaches Gibbs interval reuse

**Name:** `PERF-PH-15 changing-temperature P,H production comparison`

**Command:**

```text
cargo bench --bench equilibrium_gibbs_interval --no-default-features --features "equilibrium-instrumentation" -- --noplot
cargo bench --bench equilibrium_gibbs_interval --no-default-features --features "equilibrium-instrumentation equilibrium-benchmark-no-gibbs-interval-reuse" -- --noplot
```

**Question / hypothesis:** Is `PhaseGibbsSource` actually used by a
production changing-temperature P,H route, and does its interval-aware reuse
reduce complete facade wall time? The benchmark uses the public facade and
compares the current native path with the matched control build that forces a
Gibbs refresh for every changed temperature.

**Cases:**

| case | route | points | purpose |
| --- | --- | ---: | --- |
| TP-1907 | nested P,H | 3 | negative control: compatibility nested route does not consume `PhaseGibbsSource` |
| STANJAN-15 | monolithic P,H | 3 | production caller of `ResolvedThermochemistry::evaluate_gibbs` |

**Result:** The STANJAN-15 monolithic row consistently reported `1` Gibbs
refresh and `41` interval reuses in the current build, versus `42` refreshes
and `0` reuses in the control. One representative run was `136.293 ms`
current versus `146.162 ms` control. However, ten alternating wall samples
gave `138.836 ms` current versus `137.651 ms` control, a difference of
`-0.86%` for the current path, which is within runtime noise. The nested
TP-1907 row reported `0/0` in both builds, confirming that it is not a caller
of this cache.

**Decision:** Confirmed: the new Gibbs closure construction/reuse is on the
production monolithic P,H path. The local closure optimization is accepted,
but no complete P,H speedup is claimed: the 3-point range is dominated by
nonlinear solve, orchestration, and validation. Any further performance work
must target those larger stages or use a longer P,H range to amplify this
small contribution. Correctness remained green for all three accepted points
in both matched builds.

### Story 100: PERF-NL-01 Legacy callback same-iterate measurement

**Name:** `PERF-NL-01 shared numeric evaluation workspace gate`

**Command:**

```text
cargo bench --bench equilibrium_legacy_callbacks --no-default-features --features "equilibrium-instrumentation" -- --noplot
```

The dashboard was executed from the optimized benchmark binary after Cargo
finished its release/LTO build. The benchmark uses the public facade and one
selected Legacy backend per case.

**Question / hypothesis:** Do residual and Jacobian callbacks repeatedly
receive the identical trial iterate, and is the duplicated numeric state
construction large enough to justify a solver-attempt-local workspace? The
measurement includes Legacy LM, NR, and TR; P,T and P,H; small, STANJAN-15,
TP-1907/18, and RP-1311/40 cases. The accepted table and complete attempt
ledger retain accepted, rejected, and failed attempts separately.

**What was added:** The optional Legacy diagnostic scope now records adjacent
residual/Jacobian callback pairs with the same floating-point bit signature.
The report exposes `legacy_same_iterate_callback_pairs`; RST rows remain
`n/a`. The production build without `equilibrium-instrumentation` does not
compute this signature.

**Representative release result:**

| case | backend | route | backend wall | residual + Jacobian | same-iterate pairs | pair share |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| RP-1311/40, 40 components | Legacy LM | P,T | `0.175 ms` | `0.043 ms` | `10` | `100.0%` |
| RP-1311/40, 40 components | Legacy TR | P,H nested | `0.504 ms` | `0.213 ms` | `42` | `100.0%` |
| STANJAN-15, 15 components | Legacy LM | P,H monolithic | `1.370 ms` | `0.497 ms` | `35` | `100.0%` |
| TP-1907/18, 18 components | Legacy TR | P,H nested | `0.155 ms` | `<0.001 ms` | `32` | `100.0%` |

The exact pair count matched `min(residual_calls, jacobian_calls)` in the
representative accepted rows, confirming the duplicate-evaluation hypothesis.
However, even the largest measured callback budget was below `0.5 ms`, and
all visible started-attempt rows remained below the project materiality gate
of `10 ms` or `1%` of the complete request.

**Correctness:** The focused diagnostic unit test passed. The benchmark
completed all 36 configured facade operations with the existing accepted or
expected failure outcomes. No equation, tolerance, backend order, callback
contract, or continuation state was changed.

**Decision:** Reject the shared numeric workspace implementation for the
current fixture scale. Keep the counter as evidence for future scaling work;
reopen the optimization only when a valid larger or sustained-throughput
workload makes repeated callback state construction material. Do not infer
that a `100%` pair share alone warrants a cache.

### Story 101: PERF-RST-02 RST worker scheduling and Jacobian lock gate [HISTORICAL]

**Status:** Historical pre-migration worker-scaling capture. The current
`equilibrium_rst_workers` source now uses the prepared adapter and includes
Backtracking LM; the tables below cover only the earlier four-method matrix.
They remain useful for locating the old worker-scaling signal, but are not a
complete current solver-policy or wall-time baseline. Current policy crossover
evidence is Story 108, and current prepared-route medians are Story 103.

**Name:** `PERF-RST-02 RST callback scheduling/materiality gate`

**Command:**

```text
$env:RAYON_NUM_THREADS="1"
$env:RST_WORKERS_FILTER="rp1311_40"
$env:RST_WORKERS_REPORT_REPETITIONS="1"
cargo bench --bench equilibrium_rst_workers --no-default-features --features "equilibrium-instrumentation" -- --noplot rp1311_40
```

The same optimized benchmark executable was then run in fresh processes with
`RAYON_NUM_THREADS=2`, `4`, and unset (`default`, reported as 24 workers).
The benchmark uses the normal `EquilibriumCalculator` facade and retains a
separate stage and attempt tables. The stage table reports every existing
facade, P,T/P,H workflow, and inner equilibrium timing against its local
scope denominator; the attempt table reports every measured RST attempt,
including rejected attempts. The focused target is registered in `Cargo.toml` as
`equilibrium_rst_workers`; it does not alter solver policy or production
closure selection.

**Question / hypothesis:** Does the Rayon worker scheduler and the shared
per-entry Jacobian write lock consume a material part of real RST runtime, so
that lock-free indexed assembly or a sequential small-system closure is
warranted? The test compares 15-component monolithic elemental P,H and
18-/40-component elemental P,T cases across RST LM, Minpack LM, Nielsen LM,
and Damped Newton. The macro table keeps facade/workflow/backend scopes
compact. The stage table expands facade, P,T or P,H transaction, and inner
equilibrium scopes into concrete rows with `stage ms`, `scope ms`, and `% of
scope`. The attempt table keeps backend wall, residual, Jacobian, and solver
overhead in separate columns and includes call counts.

**Representative macro result:**

| case | RST policy | workers | status | facade ms | workflow ms | RST backend ms |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| RP-1311/40 elemental P,T | LM | 1 | accepted | `541.953` | `540.956` | `5.808` |
| RP-1311/40 elemental P,T | LM | 2 | accepted | `340.495` | `339.489` | `2.885` |
| RP-1311/40 elemental P,T | LM | 4 | accepted | `233.448` | `232.469` | `2.164` |
| RP-1311/40 elemental P,T | LM | 24 | accepted | `174.988` | `173.945` | `4.070` |
| STANJAN-15 elemental monolithic P,H | Minpack LM | 1 | accepted | `565.556` | `564.781` | `5.569` |
| STANJAN-15 elemental monolithic P,H | Minpack LM | 2 | accepted | `316.834` | `315.994` | `3.062` |
| STANJAN-15 elemental monolithic P,H | Minpack LM | 4 | accepted | `183.576` | `182.723` | `1.937` |
| STANJAN-15 elemental monolithic P,H | Minpack LM | 24 | accepted | `116.966` | `116.058` | `3.714` |

**Worker-attributed stage result (one diagnostic sample per fresh process):**

| case | workers | scope | stage | stage ms | scope ms | % of scope |
| --- | ---: | --- | --- | ---: | ---: | ---: |
| RP-1311/40 elemental P,T | 1 | P,T workflow | RST preparation | `532.925` | `544.919` | `97.8%` |
| RP-1311/40 elemental P,T | 1 | P,T workflow | nonlinear solve | `544.742` | `544.919` | `100.0%` |
| RP-1311/40 elemental P,T | 2 | P,T workflow | RST preparation | `330.260` | `341.599` | `96.7%` |
| RP-1311/40 elemental P,T | 2 | P,T workflow | nonlinear solve | `341.430` | `341.599` | `100.0%` |
| RP-1311/40 elemental P,T | 4 | P,T workflow | RST preparation | `221.997` | `232.451` | `95.5%` |
| RP-1311/40 elemental P,T | 4 | P,T workflow | nonlinear solve | `232.303` | `232.451` | `99.9%` |
| RP-1311/40 elemental P,T | 24 | P,T workflow | RST preparation | `158.131` | `175.223` | `90.2%` |
| RP-1311/40 elemental P,T | 24 | P,T workflow | nonlinear solve | `175.082` | `175.223` | `99.9%` |
| STANJAN-15 elemental monolithic P,H | 1 | P,H transaction | monolithic RST preparation | `581.693` | `592.207` | `98.2%` |
| STANJAN-15 elemental monolithic P,H | 1 | P,H transaction | monolithic backend execution | `5.947` | `592.207` | `1.0%` |
| STANJAN-15 elemental monolithic P,H | 4 | P,H transaction | monolithic RST preparation | `177.807` | `184.215` | `96.5%` |
| STANJAN-15 elemental monolithic P,H | 4 | P,H transaction | monolithic backend execution | `1.923` | `184.215` | `1.0%` |

The executable prints all stage rows, including zero-valued stages. The table
above keeps the material rows readable in this story record. `RST
preparation` is the measured symbolic/problem preparation portion recorded by
the RST adapter; `nonlinear solve` and `monolithic execution` are inclusive
wall-clock scopes around the runner. Therefore these rows are intentionally
not additive: the RST backend attempt table is the separate source for
residual/Jacobian/solver-overhead time.

**Attempt-level result:** RP-1311/40 LM at one worker spent `0.841 ms`
(`14.5%`) in residual, `4.821 ms` (`83.0%`) in Jacobian, and `0.145 ms`
(`2.5%`) in RST solver overhead across eight residual and eight Jacobian
calls. At four workers the split was `0.543 ms` (`25.1%`), `1.425 ms`
(`65.9%`), and `0.195 ms` (`9.0%`). STANJAN-15 Minpack LM at one worker
spent `2.091 ms` residual, `3.228 ms` Jacobian, and `0.249 ms` overhead;
all ten residual and ten Jacobian calls were retained in the attempt ledger.
Nielsen LM and Damped Newton rejected on the representative fixtures; their
attempt rows remain evidence of solver behavior, not accepted performance
comparisons.

**Correctness:** The focused release executable completed its configured
facade operations with the existing accepted/rejected outcomes. No equations,
seeds, tolerances, backend order, closure route, or solution publication
logic changed. `cargo check --bench equilibrium_rst_workers --no-default-features
--features equilibrium-instrumentation` passed.

**Decision:** The RST callback and its remaining solver overhead are below the
project `10 ms` optimization threshold. Reject lock-removal and sequential
closure implementation for now. The request-level worker scaling is real,
but the accepted RST attempt is only a few milliseconds while the facade and
workflow are hundreds of milliseconds; this points to a separate
orchestration/normalization or repeated-preparation budget. The worker-
attributed stage table is now the instrument for the next runs: its rows are
inclusive, so a parent such as `nonlinear solve` contains callback work and
must not be added to its children. A future RST change requires comparing the
same stage rows across fresh 1/2/4/default-worker processes and repeated
samples first.

### Story 102: PERF-RST-03 RST preparation decomposition [SUPERSEDED]

**Name:** `PERF-RST-03 RST preparation substage attribution (historical)`

**Status:** Superseded on 2026-09-17 by the prepared RST adapter migration.

**Command:**

```text
$env:RAYON_NUM_THREADS="1"
$env:RST_WORKERS_FILTER="rp1311_40_pt_element_inventory_rst_lm"
$env:RST_WORKERS_REPORT_REPETITIONS="1"
cargo bench --bench equilibrium_rst_workers --no-default-features --features "equilibrium-instrumentation" -- --noplot rp1311_40_pt_element_inventory_rst_lm
```

The command and table below belong to the pre-migration KiThe-side timing
wrappers. They are retained as historical evidence only and must not be used
as the current performance baseline. The prepared adapter now delegates graph
construction and lambdification to RustedSciThe's
`PreparedSymbolicNonlinearProblem`. RustedSciThe 0.4.18 exposes one aggregate
`preparation_duration` for this direct Lambdify path, not separate public
durations for symbolic differentiation, residual lambdification, and Jacobian
lambdification. The current dashboard therefore prints one explicitly named
`RST backend preparation (aggregate; internal split unavailable)` row. It no
longer prints zero-valued child rows, which would falsely imply measured
zero-cost work.

**Historical question / hypothesis:** Which concrete operation owned the
hundreds of milliseconds previously reported only as `RST preparation`, and
did that owner change with P,T versus monolithic P,H or with the Rayon worker
count?

**Historical release diagnostic result:** A fresh optimized executable
captured the following pre-migration split for RP-1311/40 P,T:

| workers | RST preparation | equation generation | symbolic differentiation | residual lambdification | Jacobian lambdification | problem materialization |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `539.271 ms` | `104.253 ms` | `415.322 ms` | `1.741 ms` | `17.926 ms` | `434.995 ms` |
| 2 | `323.302 ms` | `102.742 ms` | `205.815 ms` | `1.158 ms` | `13.558 ms` | `220.537 ms` |
| 4 | `229.110 ms` | `103.997 ms` | `112.770 ms` | `1.246 ms` | `11.067 ms` | `125.089 ms` |
| 24 | `156.739 ms` | `103.518 ms` | `37.481 ms` | `1.426 ms` | `14.279 ms` | `53.193 ms` |

For STANJAN-15 monolithic P,H at one worker, the historical RST preparation
scope was
`555.608 ms`: symbolic differentiation `536.590 ms`, Jacobian
lambdification `6.980 ms`, equation generation `6.004 ms`, row scaling
`0.768 ms`, and residual lambdification `1.825 ms`. The corresponding parent
`problem materialization` was `545.398 ms`; it includes differentiation and
lambdification and must not be added to them. The release samples show that
equation generation is nearly worker-invariant, while symbolic differentiation
is the stage that benefits from additional workers.

**Historical correctness:** `cargo check --lib --no-default-features` and
`cargo check --bench equilibrium_rst_workers --no-default-features
--features equilibrium-instrumentation` passed for the old instrumentation.
Those checks do not validate the current prepared adapter.

**Decision:** The old table is closed as a historical measurement. The current
optimization boundary is the aggregate RST preparation cost plus the native
RST callback table. Further differentiation/lambdification optimization
requires a RustedSciThe telemetry/API change; KiThe must not infer those
sub-stages from aggregate timing.

### Story 104: PERF-RST-04 prepared-backend release baseline [HISTORICAL]

**Status:** Historical post-migration, parallel-default capture from
2026-09-17. The architecture is current, but the numbers predate the measured
Sequential production-default decision and omit the new crossover evidence.
Use Story 103 for the current cold/warm Sequential baseline and Story 108 for
the Sequential/Parallel policy measurement. Rerun this method matrix before
using its medians to choose an RST solver.

**Name:** `PERF-RST-04 prepared RST route baseline after migration`

**Commands:**

```text
cargo bench --bench equilibrium_rst_preparation --no-default-features \
  --features equilibrium-instrumentation -- --noplot
```

The older `equilibrium_rst_workers` command belongs to the historical
pre-policy worker capture in Story 101. It is not a current Sequential versus
Parallel experiment because the production adapter now selects Sequential
explicitly.

**Question / hypothesis:** After moving the production adapter to the
parameterized prepared RST problem, what is the current cold-start wall time
for representative P,T and P,H element-inventory cases, and does the
diagnostic dashboard identify all available work without fabricating missing
RST sub-stages?

**Current release capture:** Criterion medians from the optimized benchmark
run on 2026-09-17 are shown below. They are representative route/method
baselines, not a claim that one method is accepted for every case.

| case | RST method | route | benchmark outcome | Criterion median |
| --- | --- | --- | --- | ---: |
| `stanjan15_pt_molecular` | LM | P,T | returned | `11.018 ms` |
| `tp1907_18_pt_molecular` | LM | P,T | returned | `16.080 ms` |
| `rp1311_40_pt_molecular` | LM | P,T | returned | `190.710 ms` |
| `stanjan15_ph_elemental` | Minpack LM | P,H monolithic | returned | `125.710 ms` |
| `tp1907_18_ph_elemental` | LM | P,H nested | returned | `33.242 ms` |
| `rp1311_40_ph_elemental` | LM | P,H nested | returned | `261.720 ms` |

The same run also measured the RST alternatives. For P,T TP-1907/18, the
medians were LM `16.080 ms`, Minpack LM `22.408 ms`, Nielsen LM `18.804 ms`,
and Damped Newton `16.063 ms`. For P,T RP-1311/40 they were respectively
`190.710 ms`, `187.370 ms`, `188.300 ms`, and `180.510 ms`. For large P,H,
Minpack LM was `125.710 ms` on STANJAN-15 elemental, while TP-1907/18
elemental was `35.524 ms` with Minpack LM and RP-1311/40 elemental was
`272.400 ms`.

The dashboard now has separate inclusive macro-stage and RST-attempt tables.
The RST preparation row is explicitly aggregate when the dependency does not
expose an internal split. Native callback telemetry remains separate and
reports residual/Jacobian/linear-solve call counts, durations, termination,
and solver overhead. Inclusive parent rows must not be added to their child
rows.

For the focused RP-1311/40 P,T worker run with one worker, the current tables
reported `facade 543.353 ms`, `P,T workflow 542.422 ms`, and an RST
preparation scope of `524.725 ms`. Its available rows were:

| scope | stage or backend metric | time | local share |
| --- | --- | ---: | ---: |
| P,T RST preparation | equation generation | `101.066 ms` | `19.3%` |
| P,T RST preparation | RST backend preparation aggregate | `423.628 ms` | `80.7%` |
| RST attempt | residual callbacks (`18` calls) | `1.004 ms` | `12.6%` of backend |
| RST attempt | Jacobian callbacks (`8` calls) | `6.810 ms` | `85.6%` of backend |
| RST attempt | solver overhead | `0.082 ms` | `1.0%` of backend |

This is intentionally not presented as a complete decomposition of the
`423.628 ms`: the current RST dependency does not expose that internal split.
It is, however, enough to distinguish KiThe/RST preparation cost from the
actual nonlinear callbacks. The old symbolic-differentiation numbers in Story
102 cannot be substituted into this table.

**Correctness:** Focused RST adapter tests (`11 passed`),
`cargo check --bench equilibrium_rst_workers --no-default-features
--features equilibrium-instrumentation`, and `cargo fmt -- --check` passed.
The TP-1907 matrix remains the acceptance gate; it characterizes per-method
pass, accuracy-fail, and backend-fail outcomes rather than requiring every
RST method to solve every topology.

**Decision:** The prepared architecture and available telemetry are ready for
performance work. Current KiThe measurements can optimize aggregate preparation
and facade stages; a detailed split inside RST preparation is blocked on a
future RustedSciThe preparation report. The next benchmark refresh must rerun
the prepared-session cold/warm cases and update Story 103 before using its old
numbers for decisions.

### Story 103: PERF-PREP-03 final KiThe cold/warm payload experiment

**Name:** `PERF-PREP-03 parameterized RST payload removes repeated symbolic preparation`

**Commands:**

```text
cargo bench --bench equilibrium_prepared_session --no-default-features \
  --features equilibrium-instrumentation -- --noplot
cargo bench --bench equilibrium_ph_preparation --no-default-features \
  --features equilibrium-instrumentation -- --noplot prepared_session_two_points
```

**Question / hypothesis:** Does the public KiThe prepared-session boundary
reuse the parameterized residual/Jacobian payload when temperature changes in
P,T and when `target_enthalpy` changes in monolithic P,H? The expected result
is one symbolic payload build, numeric parameter refreshes on later points,
unchanged accepted physical solutions, and a material reduction for repeated
large requests.

**Historical release capture (superseded by the 2026-09-17 refresh below):**
The earlier optimized executables were run after the detailed RST preparation
instrumentation was added. The P,T capture used
the public `EquilibriumCalculator` facade and compared a two-point one-shot
loop with one `PreparedEquilibriumPtSession`:

| case | components | preparation | cold operation | hot point | saved vs cold |
| --- | ---: | ---: | ---: | ---: | ---: |
| `tp1907_18_elemental_pt` | 18 | 14.772 ms | 27.992 ms | 0.567 ms | 98.0% |
| `rp1311_40_molecular_pt` | 40 | 155.155 ms | 162.602 ms | 2.804 ms | 98.3% |
| `stanjan15_elemental_pt` | 15 | 8.583 ms | 13.607 ms | 0.651 ms | 95.2% |

The equal two-point wall comparison was:

| case | one-shot loop | prepared session | built-in range | session saving |
| --- | ---: | ---: | ---: | ---: |
| `tp1907_18_elemental_pt` | 57.747 ms | 27.561 ms | 29.739 ms | 52.3% |
| `rp1311_40_molecular_pt` | 354.454 ms | 166.646 ms | 179.454 ms | 53.0% |
| `stanjan15_elemental_pt` | 25.663 ms | 14.174 ms | 13.485 ms | 44.8% |

For monolithic P,H, the second point used a different enthalpy target on the
same prepared session:

| case | one-shot two points | preparation | first point | second point | prepared total | saving |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `small3_molecular` | 18.044 ms | 7.236 ms | 1.037 ms | 0.463 ms | 8.736 ms | 51.6% |
| `stanjan15_molecular` | 242.563 ms | 78.783 ms | 28.243 ms | 4.006 ms | 111.032 ms | 54.2% |
| `stanjan15_elemental` | 210.935 ms | 77.749 ms | 5.401 ms | 3.641 ms | 86.791 ms | 58.9% |
| `rp1311_40_elemental` | 2421.277 ms | 1185.289 ms | 49.382 ms | 9.224 ms | 1243.895 ms | 48.6% |

The P,H capture reports `first_builds=0`, `first_reuses=1`,
`second_builds=0`, and `second_reuses=1` for each prepared session. The
`preparation_ms` value is paid once; the later target changes only numeric
state and continuation seed. The large RP-1311/40 case therefore removes
approximately `1177 ms` from this two-target workflow without changing the
equations or solver policy.

**2026-09-17 prepared-path refresh (historical parallel-default capture):**
The first post-migration capture used the then-current parallel lambdification
policy. Its numbers remain useful for comparison, but they are no longer the
production baseline after the measured sequential-policy decision below. The
P,T two-point macro capture was:

| case | one-shot loop | prepared session | built-in range | session saving |
| --- | ---: | ---: | ---: | ---: |
| `small3_molecular_pt` | `1.245 ms` | `1.013 ms` | `0.999 ms` | `18.6%` |
| `stanjan15_molecular_pt` | `22.625 ms` | `11.354 ms` | `11.804 ms` | `49.8%` |
| `tp1907_18_elemental_pt` | `46.774 ms` | `23.979 ms` | `24.050 ms` | `48.7%` |
| `rp1311_40_molecular_pt` | `371.540 ms` | `186.840 ms` | `186.320 ms` | `49.7%` |
| `small3_elemental_pt` | `1.889 ms` | `1.561 ms` | `1.544 ms` | `17.4%` |
| `stanjan15_elemental_pt` | `24.194 ms` | `12.381 ms` | `12.385 ms` | `48.8%` |

For the public P,H prepared-session two-point Criterion group, the current
medians are:

| case | prepared two-point median | current Criterion change | stage sample: preparation / first / second |
| --- | ---: | ---: | ---: |
| `small3_molecular` | `9.722 ms` | `-1.2%` | `8.884 / 1.188 / 0.442 ms` |
| `stanjan15_molecular` | `188.660 ms` | `+24.4%` | `97.195 / 54.217 / 9.878 ms` |
| `stanjan15_elemental` | `139.820 ms` | `+21.0%` | `84.900 / 9.265 / 10.505 ms` |
| `rp1311_40_elemental` | `1499.700 ms` | `+15.7%` | `1225.797 / 205.694 / 42.787 ms` |

The P,H dashboard identifies the accepted methods and stage ownership. For
example, the current STANJAN-15 molecular point accepted RST Backtracking LM;
the large P,H cases remain preparation-dominated, while the changing-target
solve is still much cheaper than rebuilding the formulation. The positive
Criterion changes are a performance regression against the stored benchmark
history, not a correctness failure; they require a separate cascade/backend
attribution pass before more optimization is attempted.

**2026-09-17 sequential production-default refresh:** The prepared RST
adapter now constructs its symbolic payload with
`LambdifyExecutionPolicy::Sequential`. The parallel policy remains available
as an explicit RST option, but the current KiThe matrix is below the point
where Rayon scheduling pays for itself. The fresh P,T two-point capture is:

| case | one-shot loop | prepared session | built-in range | session saving |
| --- | ---: | ---: | ---: | ---: |
| `small3_molecular_pt` | `0.952 ms` | `0.380 ms` | `0.382 ms` | `60.1%` |
| `stanjan15_molecular_pt` | `19.590 ms` | `9.549 ms` | `10.054 ms` | `51.3%` |
| `tp1907_18_elemental_pt` | `40.260 ms` | `19.225 ms` | `19.797 ms` | `52.2%` |
| `rp1311_40_molecular_pt` | `382.669 ms` | `177.238 ms` | `192.212 ms` | `53.7%` |
| `small3_elemental_pt` | `0.819 ms` | `0.300 ms` | `0.338 ms` | `63.4%` |
| `stanjan15_elemental_pt` | `21.239 ms` | `9.812 ms` | `11.898 ms` | `53.8%` |

The fresh P,H Criterion medians for the prepared two-target workflow are:

| case | prepared two-point median | stage sample: preparation / first / second |
| --- | ---: | ---: |
| `small3_molecular` | `8.123 ms` | `7.646 / 0.173 / 0.085 ms` |
| `stanjan15_molecular` | `173.892 ms` | `82.654 / 47.431 / 8.918 ms` |
| `stanjan15_elemental` | `132.829 ms` | `89.128 / 7.647 / 8.441 ms` |
| `rp1311_40_elemental` | `1570.050 ms` | `1136.992 / 257.841 / 59.347 ms` |

The stage columns are one diagnostic capture and are intentionally not added
to the Criterion median: they are separate observations and the RP case has
visible preparation variance. The important range result is stable: for
medium and large P,T cases a prepared session cuts the two-point workflow by
roughly `51--54%`, while the second P,H target reuses the prepared payload and
is much cheaper than a fresh construction. The small cases also benefit from
sequential execution because they avoid Rayon scheduling overhead.

This capture makes Sequential the current KiThe production default. Before
changing that decision, add explicit larger-system measurements (at least
80--100 equations) comparing Sequential and Parallel; the present evidence
does not claim that Sequential wins for arbitrarily large RST systems.

**Correctness:** Both release benchmark commands completed successfully. The
existing prepared P,T and prepared P,H tests remain the correctness gate and
assert accepted temperatures/enthalpy, element balance where applicable,
formulation reuse, and failed-target seed rollback. The complete post-
migration correctness result is recorded in Story 106; the wider library
suite was not globally green because of 18 unrelated failures.

**Decision:** The KiThe-side experiment is complete. Parameterized payload
reuse is a material optimization for repeated P,T/P,H calculations, but it
does not remove the cold-start cost. Further improvements to the symbolic
payload belong in RustedSciThe: reusable immutable expression/Jacobian/
lambdified structures, efficient numeric parameter refresh, and a common
opt-in diagnostic contract for all nonlinear backends. That work is recorded
as a separate RST handoff rather than mixed into the KiThe facade. The KiThe
default execution policy is Sequential for the current problem-size matrix;
Parallel remains an explicit experiment/escape hatch until the larger-system
threshold is measured.

### Story 105: RST backtracking LM added to the production cascade

**Name:** `PERF-RST-05 backtracking Levenberg-Marquardt cascade coverage`

**Question / hypothesis:** Does the new RustedSciThe
`BacktrackingLevenbergMarquardtMethod` work through KiThe's prepared symbolic
adapter under the same bounds, residual/Jacobian payload, acceptance gate,
and telemetry contract as the existing RST methods? It is placed after the
established RST methods and before Legacy fallback as the closest RST
counterpart to the robust legacy LM route; its position remains provisional
until the matrix benchmark is complete.

**Implementation contract:** The adapter maps the method in both production
and test-history solve paths, exposes the stable name
`rst_backtracking_levenberg_marquardt`, includes it in
`RustedSciTheSolver::recommended_cascade()` after the established RST methods,
and includes it in the RST
preparation and worker-attribution benchmark matrices. It uses the same
prepared `residual/Jacobian` binding and native RST statistics as every other
RST backend.

**Correctness command:**

```text
cargo test equilibrium_rst_backend --lib --no-default-features -- --nocapture
cargo test equilibrium_rst_matrix_tests --lib --no-default-features -- --nocapture
cargo test rusted_scithe_default_policy_has_a_stable_documented_backend_order --lib --no-default-features -- --nocapture
```

**Required current benchmark commands:**

```text
cargo bench --bench equilibrium_rst_preparation --no-default-features \
  --features equilibrium-instrumentation -- --noplot
$env:RST_CROSSOVER_PROBE_ONLY = "1"
cargo bench --bench equilibrium_rst_policy_crossover --no-default-features \
  -- --noplot
Remove-Item Env:RST_CROSSOVER_PROBE_ONLY
```

`equilibrium_rst_workers` is retained as historical request-level evidence in
Story 101. It is not a current Backtracking or Sequential-versus-Parallel
acceptance command because the production adapter selects `Sequential`
explicitly.

**Acceptance:** The method must be judged by the same physical validation as
the existing methods: finite published moles, physical element-balance
error, residual tolerance, and topology. A backend that terminates but fails
the common acceptance gate is recorded as `ACCURACY_FAIL`, not as a successful
solver. Its residual/Jacobian calls, linear solves, and solver overhead must
appear in the native RST attempt table.

**Robustness role:** Backtracking LM is retained as a late rescue method, not
as a first-choice performance method. Its line-search step reduction can
accept a point after the earlier RST methods reject a candidate or encounter
an ill-conditioned Jacobian. Therefore a few additional milliseconds can be
the correct trade for turning a failed cascade into an accepted physical
solution. The PERF-PH-14 A/B capture confirms a measurable wall-time cost,
but does not justify removing the method merely from that cost.

**Initial TP-1907 result:** On the upper-inventory `scale=1e8` frozen case,
the method was reached as the last RST fallback and returned `ACCURACY_FAIL`
after `500` iterations. The normalized route reported residuals around
`8.75e1` and backend wall time around `326 ms`; the existing Trust-Region LM
still passed before this fallback was reached. This is a correctness and
ordering characterization, not evidence that backtracking LM should become
the default first method.

**Current TP-1907 release acceptance snapshot:** The same normalized problem,
seed, bounds, and physical acceptance gate were used for every isolated
backend. The facade test passes because at least one RST method and at least
one Legacy method pass the physical gate; individual method failure is retained
as diagnostic evidence.

| backend | status | residual evidence | max physical balance | interpretation |
| --- | --- | ---: | ---: | --- |
| RST LM | `ACCURACY_FAIL` | rejected: `5.074e-8` best | n/a | above the `1e-12` candidate gate |
| RST Backtracking LM | `ACCURACY_FAIL` | rejected: `8.726e1` best | n/a | stalled after the fallback budget |
| RST Minpack LM | `ACCURACY_FAIL` | rejected: `1.787e-11` best | n/a | close, but above `1e-12` |
| RST Nielsen LM | `ACCURACY_FAIL` | rejected: `9.435e-7` best | n/a | insufficient final residual |
| RST Trust-Region LM | `PASS` | `2.712e-14` / `2.712e-14` | `1.192e-7` | accepted RST method |
| RST Damped Newton | `ACCURACY_FAIL` | rejected: `6.312e1` best | n/a | failed to reach the physical candidate |
| Legacy LM | `PASS` | `1.557e-14` / `1.557e-14` | `1.192e-7` | accepted Legacy method |
| Legacy NR | `ACCURACY_FAIL` | `1.693e-14` / `1.693e-14` | `3.532e-6` | formal residual passes, balance fails |
| Legacy TR | `PASS` | `1.474e-14` / `1.474e-14` | `1.192e-7` | accepted Legacy method |

This confirms that solver termination or a small formal residual is not alone
an acceptance proof: `Legacy NR` fails the physical balance gate, while RST
Trust-Region LM supplies the passing RST route.

### Story 106: post-migration full correctness matrix

**Name:** `CORR-RST-06 prepared RST migration full-suite audit`

**Command:**

```text
cargo test --lib --no-default-features -- --include-ignored --nocapture
```

**Question / hypothesis:** Does moving the production RST adapter to the
prepared parameterized problem preserve the calculator and frozen-reference
correctness contracts, especially the P,H lifecycle and graphite topology
routes?

**Observed full-suite result on 2026-09-17:** `2544 passed, 18 failed,
0 ignored, 0 measured`. The complete library suite is therefore not globally
green. The 18 failures are outside the equilibrium/RST migration scope:
Kinetics, live-data, NIST fitting, User_PhaseOrSolution, GUI lifecycle, and
settings tests.

**Equilibrium-specific follow-up:** The two P,H frozen lifecycle tests listed
by the full run were executed directly after the suite:

```text
cargo test --lib --no-default-features -- --include-ignored --exact \
  Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::i5_nasa_tp1906_tp1907_chon_graphite_ph_production_lifecycle_characterization \
  --nocapture
cargo test --lib --no-default-features -- --exact \
  Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_four_point_ph_lifecycle_preserves_pt_manifold_and_graphite_topology \
  --nocapture
```

Both passed. The ignored characterization retained the expected graphite
topology transitions, accepted P,H continuation in both directions, and
finite scaled enthalpy errors. The focused RST adapter, TP-1907 backend
matrix, and facade oracle gate also pass. Therefore this migration has no
identified equilibrium correctness regression; the unrelated 18 failures
remain separate technical debt and must not be counted as RST evidence.

**Decision:** The prepared RST migration clears the equilibrium correctness
gate. Continue with the prepared-session cold/warm benchmark refresh; do not
claim a globally green library until the unrelated failures are handled.

### Story 107: controlled P,H cascade attribution for Backtracking LM

**Name:** `PERF-PH-14 repeated Backtracking LM contribution to prepared P,H wall time`

**Command:**

```text
cargo bench --bench equilibrium_ph_preparation --no-default-features \
  --features equilibrium-instrumentation -- --noplot prepared_session_two_points
```

**Question / hypothesis:** Does the new Backtracking LM account for the fresh
P,H regression, or is the dominant change in prepared formulation/RST
preparation or in the other backend attempts? Use the same public facade,
prepared session, targets, seed, bounds, acceptance gate, and Legacy fallback;
only remove Backtracking LM from the comparison cascade.

**A/B diagnostic capture on 2026-09-17:**

| case | policy | preparation | first point | second point | total | first backend | second backend |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| `small3_molecular` | production + Backtracking | `7.552 ms` | `1.256 ms` | `0.511 ms` | `9.319 ms` | RST Minpack LM | RST LM |
| `small3_molecular` | production without Backtracking | `7.902 ms` | `1.269 ms` | `0.480 ms` | `9.651 ms` | RST Minpack LM | RST LM |
| `stanjan15_molecular` | production + Backtracking | `112.231 ms` | `57.966 ms` | `9.891 ms` | `180.088 ms` | RST Backtracking LM | Legacy LM |
| `stanjan15_molecular` | production without Backtracking | `110.938 ms` | `51.475 ms` | `8.189 ms` | `170.602 ms` | Legacy LM | Legacy LM |
| `stanjan15_elemental` | production + Backtracking | `119.274 ms` | `10.758 ms` | `10.448 ms` | `140.480 ms` | RST Minpack LM | Legacy LM |
| `stanjan15_elemental` | production without Backtracking | `116.755 ms` | `10.579 ms` | `8.862 ms` | `136.195 ms` | RST Minpack LM | Legacy LM |
| `rp1311_40_elemental` | production + Backtracking | `1297.024 ms` | `190.742 ms` | `51.368 ms` | `1539.134 ms` | RST Minpack LM | RST Minpack LM |
| `rp1311_40_elemental` | production without Backtracking | `1095.244 ms` | `235.139 ms` | `44.363 ms` | `1374.746 ms` | RST Minpack LM | RST Minpack LM |

All rows were accepted. Removing Backtracking reduced the observed total by
about `5.6%` for STANJAN molecular, `3.1%` for STANJAN elemental, and `10.7%`
for RP-1311/40 in this single diagnostic capture. The RP difference is not a
clean solver-only measurement because preparation varied by roughly `202 ms`.
The result nevertheless proves that Backtracking is a measurable late-cascade
cost and that its presence alone does not explain the full historical P,H
regression. The accepted backend can also change between points, so the
attempt table remains part of the evidence rather than just the final method
name.

**Superseding three-repeat capture on 2026-09-17:** The same public-facade
case matrix, prepared-session operation, target pair, fallback family, and
acceptance gate were run three times with and without Backtracking LM. The
values below are means of the three complete workflows; the earlier table is
retained as a single-capture historical observation.

| case | production + Backtracking | production without Backtracking | delta |
| --- | ---: | ---: | ---: |
| `small3_molecular` | `8.593 ms` | `8.307 ms` | `+3.44%` |
| `stanjan15_molecular` | `168.657 ms` | `161.942 ms` | `+4.15%` |
| `stanjan15_elemental` | `127.178 ms` | `124.562 ms` | `+2.10%` |
| `rp1311_40_elemental` | `1485.063 ms` | `1495.577 ms` | `-0.70%` |

The repeated result confirms a modest Backtracking cost on the two STANJAN
routes, while the RP-1311/40 difference is within preparation/attempt noise.
It does not explain the old P,H regression and does not justify moving the
rescue method earlier. The current benchmark source still includes
Backtracking in both the production cascade and the full RST method matrices.

**Decision:** Keep Backtracking LM in the production cascade for correctness
and robustness coverage, but do not promote it earlier. The repeated A/B gate
is complete; future changes to its implementation must preserve the accepted
backend sequence and rerun this same three-repeat comparison.

### Story 108: PERF-RST-06 Sequential versus Parallel Lambdify crossover

**Name:** `PERF-RST-06 prepared RST callback execution-policy crossover`

**Command:**

```powershell
$env:RST_CROSSOVER_PROBE_ONLY = "1"
cargo bench --bench equilibrium_rst_policy_crossover --no-default-features -- --noplot
Remove-Item Env:RST_CROSSOVER_PROBE_ONLY
```

**Question / hypothesis:** At what prepared-system size does RST's parallel
Lambdify callback policy repay Rayon scheduling overhead? The fixture keeps the
same variables, parameter binding, equations, and `min_work=1`; only
`LambdifyExecutionPolicy::Sequential` versus `Parallel` changes. Preparation
and repeated runtime callbacks are reported separately because preparation
policy timing does not answer the solver-loop question.

**Fixture and measurement:** The probe uses independent parameterized symbolic
equations with sizes `3, 15, 18, 40, 64, 80, 100, 128, 256, 512, 1024`.
Preparation values are the mean of three constructions. Callback values are
the mean per `residual_into + jacobian_into` pair over three samples of forty
repetitions. The callback ratio is `parallel / sequential`; values below `1`
mean that Parallel is faster.

**Optimized release probe result on 2026-09-17:**

| equations | prep Sequential / Parallel | callbacks Sequential / Parallel | callback ratio | callback winner |
| ---: | ---: | ---: | ---: | --- |
| 3 | `0.292 / 0.137 ms` | `0.2 / 42.3 us` | `181.946x` | Sequential |
| 15 | `0.382 / 0.258 ms` | `1.0 / 33.7 us` | `33.920x` | Sequential |
| 18 | `0.333 / 0.340 ms` | `1.6 / 37.1 us` | `22.680x` | Sequential |
| 40 | `0.628 / 0.494 ms` | `4.2 / 34.8 us` | `8.260x` | Sequential |
| 64 | `0.891 / 0.781 ms` | `10.5 / 45.9 us` | `4.379x` | Sequential |
| 80 | `0.995 / 1.011 ms` | `13.5 / 57.0 us` | `4.235x` | Sequential |
| 100 | `1.374 / 1.220 ms` | `20.3 / 59.9 us` | `2.959x` | Sequential |
| 128 | `1.744 / 1.706 ms` | `30.8 / 121.8 us` | `3.959x` | Sequential |
| 256 | `4.488 / 4.082 ms` | `115.3 / 302.1 us` | `2.620x` | Sequential |
| 512 | `13.121 / 12.910 ms` | `455.4 / 551.5 us` | `1.211x` | Sequential |
| 1024 | `47.982 / 47.800 ms` | `1862.1 / 1544.5 us` | `0.829x` | Parallel |

**Conclusion:** For this workload the crossover is bracketed by `512` and
`1024` equations. Parallel preparation is roughly comparable and therefore is
not the decision point; the crossover comes from callback arithmetic finally
becoming large enough to amortize Rayon. Sequential remains the correct
default for the current KiThe production matrix (`15--40` components and the
tested large cases), while Parallel becomes a candidate for genuinely
thousand-equation systems.

This is a policy threshold measurement, not a chemical correctness test. The
fixture has independent cubic scalar expressions and does not model phase
control, dense coupling, or a real TP/RP graph. Before changing the default,
repeat the same experiment with replicated real symbolic payloads and record
the solver iteration wall time, not only isolated callback time.

### Story 109: PERF-RST-07 current five-method prepared route matrix

**Name:** `PERF-RST-07 current Sequential RST method suitability matrix`

**Command:**

```powershell
$env:RST_PREPARATION_REPORT_ONLY = "1"
cargo bench --bench equilibrium_rst_preparation --no-default-features `
  --features equilibrium-instrumentation -- --noplot
Remove-Item Env:RST_PREPARATION_REPORT_ONLY
```

For a terminal-sized solver summary, set the compact switch as well:

```powershell
$env:RST_PREPARATION_REPORT_ONLY = "1"
$env:RST_PREPARATION_COMPACT = "1"
cargo bench --bench equilibrium_rst_preparation --no-default-features `
  --features equilibrium-instrumentation -- --noplot
Remove-Item Env:RST_PREPARATION_REPORT_ONLY
Remove-Item Env:RST_PREPARATION_COMPACT
```

To retain every stage and every backend attempt without terminal-width loss,
redirect the full mode to a UTF-8 text file and inspect its sections separately:

```powershell
$env:RST_PREPARATION_REPORT_ONLY = "1"
cargo bench --bench equilibrium_rst_preparation --no-default-features `
  --features equilibrium-instrumentation -- --noplot `
  *> .\target\rst-preparation-full.txt
Remove-Item Env:RST_PREPARATION_REPORT_ONLY
```

The compact table is the current matrix summary; the file remains the
authoritative detailed stage/attempt dashboard.

**Question / hypothesis:** After the prepared RST migration and the
Sequential production-policy decision, which of the five available RST
methods can solve the same facade routes, and what is the current cold-route
wall time? The matrix must distinguish physical acceptance from rejection
time; a fast rejected method is not a successful solver.

**Scope and measurement:** The report-only run used the normal
`EquilibriumCalculator` facade, fixed `SolverPolicy::Single` routes, the same
seeds/bounds/acceptance rules for each method, and the production prepared
RST payload. It covered 51 fixed-policy routes, five methods, and three
repetitions per route: 153 observations. Values below are medians of the three
wall samples. `FAIL` means the isolated method did not publish an accepted
solution; its wall time is retained only as rejection-cost evidence.

**Current release summary (2026-09-17, latest compact rerun):**

| route | RST LM | Backtracking LM | Minpack LM | Nielsen LM | Damped Newton |
| --- | ---: | ---: | ---: | ---: | ---: |
| small3 molecular P,T | `0.34 P` | `0.35 P` | `0.42 P` | `0.55 P` | `0.51 P` |
| small3 elemental P,T | `0.52 F` | `0.41 P` | `0.31 P` | `0.47 P` | `0.45 F` |
| small3 elemental P,H | `9.73 P` | `10.25 P` | `7.82 P` | `8.41 F` | `7.94 P` |
| STANJAN-15 molecular P,T | `10.67 P` | `42.77 P` | `10.88 P` | `27.24 F` | `23.60 F` |
| STANJAN-15 elemental P,T | `12.65 P` | `28.10 F` | `12.19 P` | `10.92 P` | `23.46 F` |
| STANJAN-15 molecular P,H | `748.98 F` | n/a | `413.24 P` | n/a | n/a |
| STANJAN-15 elemental P,H | `519.88 F` | n/a | `122.65 P` | n/a | n/a |
| TP-1907/18 molecular P,T | `16.21 F` | `19.50 F` | `24.12 F` | `24.39 F` | `15.96 F` |
| TP-1907/18 elemental P,T | `38.94 F` | n/a | `40.33 F` | n/a | n/a |
| TP-1907/18 elemental P,H | `39.51 F` | n/a | `47.56 F` | n/a | n/a |
| RP-1311/40 molecular P,T | `208.13 F` | `320.83 F` | `224.61 P` | `208.28 F` | `201.53 F` |
| RP-1311/40 elemental P,T | `187.19 P` | n/a | `267.25 P` | n/a | n/a |
| RP-1311/40 elemental P,H | `534.81 P` | n/a | `492.86 P` | n/a | n/a |

`P` means physically accepted; `F` means `failed:AllBackendsFailed` for the
isolated policy. The P,H and larger ElementInventory rows intentionally use
only the methods supported by the benchmark's fixed-policy route matrix.
For rejected rows the compact report now prints `n/a (rejected)` for workflow
time, because the facade never publishes a successful workflow timing after
an early isolated-method rejection. The wall value remains rejection-cost
evidence and must not be compared with an accepted workflow as if it were a
successful solve.

**Interpretation:** There is no universal fastest RST method. Minpack LM is
the only accepted method in the 40-component molecular P,T row and is also
the clear accepted method for STANJAN-15 P,H. Trust-Region LM remains the
known accepted method for the TP-1907 upper-inventory gate. Nielsen LM is
useful on small and STANJAN elemental P,T but rejects several molecular and
P,H routes. Backtracking LM is accepted on small routes but is slower and
does not rescue the TP-1907/RP-1311 molecular rows in this isolated matrix.

The matrix therefore supports keeping the heterogeneous cascade and
Backtracking as a late robustness candidate. It does not justify reordering
the production cascade solely by median wall time. Preparation and callback
policy are separate questions; use Story 108 for Sequential versus Parallel,
and use the attempt table from this run for residual/Jacobian/solver-overhead
attribution.

### Story 110: PERF-PH-16 share resolved P,H recovery ownership

**Name:** `PERF-PH-16 P,H recovery request shares resolved-system snapshot`

**Command:**

```powershell
cargo test --lib --no-default-features equilibrium_ph_workflow -- --nocapture
cargo test --lib --no-default-features nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_frozen_ph_molecular_elemental_and_formal_routes_reuse_target -- --nocapture
cargo bench --bench equilibrium_ph_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  -- --noplot prepared_session_two_points
```

**Question / hypothesis:** Does P,H route fallback/recovery copy an entire
resolved phase system merely to retain the original physical request? The
correctness requirement is unchanged equations, route, accepted state,
element balances, and rollback. The performance target is only the request
ownership boundary; symbolic preparation and nonlinear execution are measured
separately and are outside this change.

**Change:** `ResolvedPhaseEnthalpyRequest` now stores its immutable
`ResolvedPhaseSystem` behind `Arc`. `cloned_with_ph_solve_mode()` and the
normalized recovery request therefore clone a shared handle. The request's
target, composition seed, options, and recovery transaction semantics remain
value-local. `ResolvedThermochemistry` and result publication were deliberately
left unchanged for a separate measurement.

**Correctness result, 2026-09-17:**

| gate | result |
| --- | --- |
| P,H workflow unit suite | `35 passed, 0 failed` |
| TP-1906 molecular, elemental, formal target-reuse route | `1 passed, 0 failed` |
| equations, solver order, acceptance and balance contract | unchanged |

**Focused release benchmark result:** Criterion measured the existing public
facade prepared-session two-target workflow. The reported interval is the
Criterion estimate; the change column is against its stored local baseline and
is not a matched clean checkout A/B series.

| case | current median | Criterion change | decision |
| --- | ---: | ---: | --- |
| `small3_molecular` | `9.0525 ms` | `-4.16%, p=0.01` | measurable local win |
| `stanjan15_molecular` | `189.93 ms` | `-6.06%, p=0.03` | modest win |
| `stanjan15_elemental` | `152.83 ms` | `+3.94%, p=0.03` | no material conclusion |
| `rp1311_40_elemental` | `1.8074 s` | `-1.18%, p=0.63` | no detected change |

The result is consistent with the hypothesis: sharing the resolved snapshot
can matter on small and medium routes, but it cannot explain the large
RP-1311/40 budget. That budget remains dominated by RST preparation and
nonlinear/backend execution. This is a partial ownership optimization, not a
completion of copy-on-failure.

**Decision:** Keep the change. The next experiment must either share the
large `ResolvedThermochemistry` payload or redesign the attempt boundary so a
successful physical solve does not need a recovery-preserving request copy.
Both require a new correctness gate before benchmarking; do not infer their
value from the current wall table alone.

### Story 112: PERF-PH-18 explain nested fixed-template preparation

**Name:** `PERF-PH-18 nested P,H template preparation child attribution`

**Command:**

```powershell
cargo check --bench equilibrium_ph_preparation --no-default-features `
  --features equilibrium-instrumentation
cargo test --lib --no-default-features equilibrium_ph_workflow -- --nocapture
cargo test --lib --no-default-features equilibrium_ph_range -- --nocapture
cargo bench --bench equilibrium_ph_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  -- --noplot point
```

**Question / hypothesis:** What is hidden inside the large nested P,H
`trial_preparation` macro-stage? The first fixed-layout template is built once
and must be distinguished from later reused trials. Only a child stage large
enough to matter should become an optimization target.

**Instrumentation change:** `PhTrialTimingReport::template_preparation()` now
retains the `EquilibriumTimingReport` produced while the fixed template is
created. The benchmark prints a separate `PERF-PH-18` table with the template
total, request, thermochemistry, equation, symbolic, numeric closure,
numerical problem, RST, active-set, projection, and unattributed portions.
Reused trials publish zero template cost. No numerical equations, solver order,
seed, phase lifecycle, rollback, or acceptance rule changed.

**Correctness result, 2026-09-18:**

| gate | result |
| --- | --- |
| P,H workflow unit suite | `35 passed, 0 failed` |
| P,H range unit suite | `6 passed, 0 failed` |
| benchmark compilation | passed |
| equations, solver order, acceptance and balance contract | unchanged |

**Release capture result:**

| case | point | template ms | symbolic-construction ms | other measured child stages | unattributed ms |
| --- | ---: | ---: | ---: | ---: | ---: |
| `tp1907_18_elemental` | `0` | `14.737` | `14.556` | each `<=0.021` | `0.118` |
| `rp1311_40_elemental` | `0` | `259.279` | `258.780` | each `<=0.076` | `0.298` |

All subsequent fixed-formulation trials in both cases reported
`template_ms=0.000`. The symbolic-construction column is the KiThe timing
around the RST prepared-problem call. Because the current RST release reports
`backend_stage_breakdown_available=false`, it is an aggregate RST preparation
measurement, not proof that differentiation, lambdification, or one specific
RST subroutine owns the full time.

**Decision:** The clone and all other KiThe child stages are not candidates for
the next optimization. The next measurement belongs inside the RST prepared
symbolic path: split graph construction, symbolic differentiation, residual
and Jacobian lambdification, and problem materialization while keeping
telemetry disabled by default. Only after that split should sequential versus
parallel RST preparation or another RST optimization be attempted. The
one-time cost is already amortized by prepared P,H ranges, so the benchmark
must report both cold template cost and reused-point cost.

### Story 111: PERF-PH-17 measure recovery request clone cost before optimizing

**Name:** `PERF-PH-17 P,H recovery-preserving request clone attribution`

**Command:**

```powershell
cargo check --bench equilibrium_ph_preparation --no-default-features `
  --features equilibrium-instrumentation
cargo test --lib --no-default-features equilibrium_ph_workflow -- --nocapture
cargo test --lib --no-default-features equilibrium_ph_range -- --nocapture
cargo test --lib --no-default-features nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_frozen_ph_molecular_elemental_and_formal_routes_reuse_target -- --nocapture
cargo bench --bench equilibrium_ph_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  -- --noplot point
```

**Question / hypothesis:** Is the request clone retained before every ordinary
P,H attempt large enough to justify a copy-on-failure or shared
`ResolvedThermochemistry` refactor? The measurement must be separate from
`trial_preparation`, `scalar_orchestration`, and nonlinear execution. A stage
below roughly one percent of a large case is not a useful optimization target
unless a later structural change makes it grow.

**Instrumentation change:** `PhTemperatureTimingReport` now exposes
`recovery_request_clone()`. The single-point transaction measures the clone
used as the working request; the range driver measures the clone retained for
classified recovery. Timing is taken only when the existing equilibrium
instrumentation is enabled and is zero otherwise. No solver, seed, bounds,
rollback, recovery, or acceptance behavior changed.

**Correctness result, 2026-09-18:**

| gate | result |
| --- | --- |
| P,H workflow unit suite | `35 passed, 0 failed` |
| P,H range unit suite | `6 passed, 0 failed` |
| TP-1906 molecular, elemental, formal target-reuse route | `1 passed, 0 failed` |
| benchmark compilation | passed |
| equations, solver order, acceptance and balance contract | unchanged |

**Release capture result:** The public-facade point capture measured one
ordinary point for each representative route. `wall_ms` is the external
point wall clock, `ph_total_ms` is the instrumented P,H transaction, and
`recovery_clone_ms` is the isolated request-copy cost.

| case | route | wall ms | P,H ms | clone ms | clone/P,H | clone/wall |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `small3_molecular` | monolithic | 8.317 | 8.088 | 0.021 | 0.3% | 0.3% |
| `stanjan15_molecular` | monolithic | 160.007 | 159.297 | 0.029 | 0.0% | 0.0% |
| `stanjan15_elemental` | monolithic | 112.686 | 111.874 | 0.039 | 0.0% | 0.0% |
| `tp1907_18_elemental` | nested | 69.189 | 67.198 | 0.058 | 0.1% | 0.1% |
| `rp1311_40_elemental` | nested | 214.829 | 200.537 | 0.104 | 0.1% | 0.0% |

**Decision:** Do not optimize this clone now. Even on the largest measured
case it is about one tenth of a millisecond and one tenth of one percent of
the P,H budget. The earlier `Arc<ResolvedPhaseSystem>` ownership change stays
in place, but another deep ownership experiment for
`ResolvedThermochemistry` is not justified by this evidence. Continue with a
stage that is already large in the telemetry; for monolithic STANJAN-15 that
is RST symbolic preparation, while for nested RP-1311/40 it is trial
preparation. If either macro-stage remains too broad, add a child timer there
before changing code.

### Story 113: PERF-RST-08 complete preparation evidence boundary

**Name:** `PERF-RST-08 RST preparation evidence is detailed without fabricated zeros`

**Commands:**

```powershell
cargo fmt --all -- --check
cargo check --lib --no-default-features
cargo check --bench equilibrium_rst_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint"
cargo test --lib --no-default-features equilibrium_rst_backend_tests::rst_preparation_telemetry_reports_aggregate_when_backend_split_is_unavailable -- --nocapture
```

**Question / hypothesis:** Can one RST preparation report identify every
measured KiThe child stage and the exact RST lifecycle boundary, without
mistaking an unavailable backend sub-stage for zero cost? This is required
before selecting another performance target.

**Instrumentation and table contract:** `RstPreparationTimingReport` now
retains RST's public effective backend, artifact action, optional AOT build
duration, and optional generated residual/Jacobian job counts. The benchmark
`equilibrium_rst_preparation` prints three complementary levels:

| level | denominator / meaning |
| --- | --- |
| macro partition | external wall versus facade/workflow; orientation only |
| concrete stage ownership | measured KiThe stages and the aggregate RST preparation call |
| RST preparation evidence | backend metadata, build/jobs, and availability of an internal RST split |

The concrete preparation rows include equation generation, temperature
substitution, row scaling, enthalpy equation, variable materialization,
parameter evaluation, log-mole bounds, the aggregate RST preparation call,
and optional artifact build. RST differentiation/lambdification rows are
printed only when the dependency reports them. With RustedSciThe 0.4.18
direct Lambdify they are explicitly `unavailable`, not `0.00 ms`.

**Correctness result, 2026-09-18:**

| gate | result |
| --- | --- |
| library check | passed |
| benchmark check with instrumentation and prettyprint | passed |
| RST preparation telemetry regression | `1 passed, 0 failed` |
| production equations, solver order, and acceptance behavior | unchanged |

**Decision:** The KiThe-side preparation budget is now sufficiently detailed
for the next optimization decision without another broad release run. The
remaining unknown is specifically inside RST's aggregate prepared-symbolic
call. It must be split in RustedSciThe into graph construction, symbolic
differentiation, residual/Jacobian lambdification, and materialization before
optimizing that call or drawing conclusions about its internal bottleneck.

**Nested-template routing correction:** The release dashboard first exposed a
second instrumentation gap: nested P,H preparation stores its cold RST build
under `trial.template_preparation`, not under `trial.inner_timing`. The RST
benchmark now records that scope separately. A fast dev-profile runtime check
after the correction produced RP-1311/40 `trial-0 template RST preparation`
rows with equation-generation and aggregate RST preparation values; these are
sanity evidence only, not a release performance baseline. The release command
above remains the authoritative performance capture command.

### Story 114: PERF-RST-09 KiThe symbolic equation-builder decomposition

**Name:** `PERF-RST-09 KiThe-owned equation preparation is separated from RST preparation`

**Question / hypothesis:** Which part of the preparation budget belongs to
KiThe, and which part begins only after KiThe has handed symbolic expressions
to RustedSciThe? The hypothesis is that
`multiphase_equilibrium_residual_generator_sym*` is KiThe work, while
expression differentiation, lambdification, and prepared-problem
materialization belong to RST.

**Command:**

```powershell
$env:RST_PREPARATION_REPORT_ONLY='1'
cargo bench --bench equilibrium_rst_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  -- --noplot
```

**Ownership contract:**

| dashboard row | owner | meaning |
| --- | --- | --- |
| `KiThe symbolic input validation` | KiThe | dimensions and input-shape checks |
| `KiThe symbolic phase topology` | KiThe | species-to-phase map and reaction phase stoichiometry |
| `KiThe symbolic variable expansion` | KiThe | construction of `y` and `exp(y)` expression nodes |
| `KiThe symbolic phase totals` | KiThe | symbolic phase mole totals |
| `KiThe symbolic activity offsets` | KiThe | phase-model pressure/reference-pressure terms |
| `KiThe symbolic reaction rows` | KiThe | reaction-affinity expression assembly |
| `KiThe symbolic element-balance rows` | KiThe | physical elemental-balance expression assembly |
| `KiThe symbolic residual simplification` | KiThe | optional `.simplify()` pass, when enabled |
| `equation generation` | KiThe | inclusive parent for the equation-generator call |
| `RST backend preparation (aggregate)` | RST boundary | public prepared-symbolic construction call |

The `equation generation` row is inclusive and must not be added to its
KiThe child rows. The RST aggregate must likewise not be added to the KiThe
rows. This prevents double counting while making the ownership boundary
explicit.

**Implementation result, 2026-09-18:** KiThe now measures all sub-stages
available before the RST call for fixed P,T, monolithic P,H, nested P,H cold
templates, and range preparation. The diagnostic collector is disabled by
default and has no timing calls on the normal path. Existing correctness tests
remain green, including five RST backend contract tests and the parity-test
compilation.

**Release result, 2026-09-18:** The optimized benchmark completed successfully
(`Finished bench profile [optimized]`, exit code `0`). The report printed the
macro table, concrete stage table, KiThe equation-builder child rows, RST
metadata table, every backend attempt, and execution/reuse evidence. The
release capture confirms the current RST public report exposes `lambdify`,
`not-applicable` artifact action, no AOT job counts, and no internal
differentiation/lambdification split. Therefore all unallocated preparation
time after the KiThe child rows is currently an RST-side question, not an
unmeasured KiThe stage.

**Representative release medians:** The following values are medians of the
three printed repetitions in the same optimized run. They use the named RST
preparation scope as denominator; `equation generation` is inclusive and is
not added to its child rows.

| case / scope | KiThe reaction rows | KiThe simplification | KiThe equation generation | RST aggregate |
| --- | ---: | ---: | ---: | ---: |
| STANJAN-15 molecular P,T | `6.22 ms` | `5.95 ms` | `6.37 ms` | `2.74 ms` |
| STANJAN-15 elemental P,T | `8.32 ms` | `7.93 ms` | `8.48 ms` | `3.80 ms` |
| RP-1311/40 molecular P,T | `115.21 ms` | `109.21 ms` | `116.38 ms` | `62.81 ms` |
| RP-1311/40 elemental P,T | `149.56 ms` | `142.14 ms` | `151.41 ms` | `81.96 ms` |
| RP-1311/40 elemental P,H `trial-0 template` | `176.39 ms` | `166.83 ms` | `177.59 ms` | `66.79 ms` |

The remaining KiThe child rows for the large cases are small: elemental
balance assembly is about `1.0--1.7 ms`, while phase topology, variable
expansion, phase totals, and activity offsets are below the displayed
`0.1 ms` resolution. The immediate KiThe-side candidate is therefore the
symbolic reaction-row builder and its simplification pass, not bounds,
parameter evaluation, or elemental-balance assembly. The RST-side request
remains necessary only to split its `62--82 ms` aggregate into graph
construction, differentiation, lambdification, and materialization.

**Route-specific qualification:** `symbolic residual simplification` is not a
universal RST stage:

| route | simplification behavior | release observation |
| --- | --- | --- |
| fixed P,T cold preparation | enabled | dominant KiThe cost for large systems |
| nested P,H cold template | enabled in the fixed-temperature inner template | dominant KiThe child cost for RP-1311/40 |
| monolithic P,H | intentionally disabled; raw expressions go directly to RST | `0.00 ms`; RST aggregate is about `100 ms` for STANJAN-15 |
| prepared/range warm point | graph is reused and parameters are rebound | no repeated equation generation/simplification |
| Legacy backend | no RST symbolic preparation | not applicable |

This distinction is part of the optimization contract. Any future attempt to
remove or reduce simplification must target the fixed-P,T/nested-template
route, preserve the existing raw monolithic P,H route, and pass symbolic
parity plus physical correctness tests before comparing release timings.

**Decision:** KiThe-side preparation is now sufficiently decomposed for a
performance decision. The next request belongs to RustedSciThe: expose child
timings for graph construction, differentiation, residual lambdification,
Jacobian lambdification, and prepared-problem materialization. Until then,
optimize only a KiThe child row that is visibly large in the detailed table;
do not infer an RST bottleneck from the aggregate label.

### Story 115: PERF-RST-10 fixed-P,T raw-expression A/B experiment

**Name:** `PERF-RST-10 removing fixed-P,T symbolic simplification is a gated experiment`

**Question / hypothesis:** Can fixed-P,T send the same KiThe symbolic
expression tree directly to RST, as monolithic P,H already does, without
losing residual/Jacobian parity, solver acceptance, or physical correctness?

**Experimental switch:** The production default is now the raw expression
route. The old simplified route remains available for a matched baseline with:

```powershell
$env:KITHE_RST_SYMBOLIC_SIMPLIFICATION='simplified'
cargo bench --profile release --bench equilibrium_rst_preparation `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint equilibrium-benchmark-simplified-symbolic" `
  -- --noplot
```

The corresponding production benchmark prints
`rst_symbolic_simplification=Raw (production)`. The environment variable
`KITHE_RST_SYMBOLIC_SIMPLIFICATION=simplified` is a diagnostic baseline
override; `raw`/`off` explicitly select production behavior. Monolithic P,H
already used the raw generator and is unchanged.

**Correctness gates:**

| gate | result |
| --- | --- |
| TP-1907 upper-inventory acceptance matrix | passed; 1 test |
| TP-1907 frozen group | `23 passed, 9 ignored, 0 failed` |
| phase-bridge correctness group | `35 passed, 0 failed` |
| ChemEquilibrium production gate | `992 passed, 0 failed, 74 ignored` |

The broader crate run still has unrelated failures outside ChemEquilibrium,
plus two TP-1906 tests blocked by a missing local thermochemistry library
file. The relevant TP-1907, phase-bridge, and full ChemEquilibrium gates pass
with the raw graph.

**Release A/B result, three repetitions:**

| case / scope | simplified total | raw total | raw simplification |
| --- | ---: | ---: | ---: |
| STANJAN-15 molecular P,T | `6.37 + 2.74 = 9.11 ms` | `0.36 + 2.80 = 3.16 ms` | `0.00 ms` |
| RP-1311/40 elemental P,T | `151.41 + 81.96 = 233.37 ms` | `6.48 + 75.46 = 81.94 ms` | `0.00 ms` |
| RP-1311/40 elemental P,H `trial-0 template` | `177.59 + 66.79 = 244.38 ms` | `6.61 + 75.92 = 82.53 ms` | `0.00 ms` |

The raw route reduces cold preparation by roughly `65--66%` on the large
cases. The RST aggregate does not increase materially and is slightly lower
in this matched capture. This proves that the large KiThe simplification cost
is avoidable for these cases; the full ChemEquilibrium gate confirms that the
promotion preserves the current route's correctness contract.

**Decision:** Promote raw fixed-P,T to the production default. Retain
`equilibrium-benchmark-simplified-symbolic` and the environment override as a
reproducible pre-change baseline. This is a preparation-path optimization,
not an equation change: both modes start from the same KiThe graph and use the
same RST solve/acceptance pipeline.

### Story 116: PERF-RST-11 RST native telemetry is attached to the calculator report

**Name:** `PERF-RST-11 RST 0.4.19 preparation and solver telemetry reaches the common equilibrium report`

**Question / hypothesis:** Does the prepared RST adapter preserve the detailed
telemetry now exposed by RustedSciThe, or does the calculator still collapse it
to one opaque preparation row and broad solver wall time? The required result is
an honest split: RST preparation children under the aggregate preparation row,
and native nonlinear-solver counters under the common backend-attempt timing
object. A stage that RST did not execute must remain unavailable rather than
appearing as a zero-cost stage.

**Commands:**

```powershell
cargo test --lib --no-default-features `
  Thermodynamics::ChemEquilibrium::equilibrium_rst_backend_tests::rst_preparation_telemetry_reports_detailed_backend_stages `
  -- --nocapture

cargo test --lib --no-default-features `
  Thermodynamics::ChemEquilibrium::equilibrium_rst_backend_tests::every_recommended_rst_strategy_solves_the_same_symbolic_problem `
  -- --nocapture

cargo check --benches --no-default-features
```

**Checks:**

| area | assertion |
| --- | --- |
| RST preparation | `backend_stage_breakdown_available=true` for the enabled 0.4.19 preparation report |
| RST preparation | differentiation, residual callback preparation, Jacobian callback preparation, input validation, prepared assembly, and unattributed time are forwarded when reported |
| unavailable stages | `FinalValidation` and other non-applicable records remain `None`; they are not fabricated zeros |
| RST solve | every recommended RST strategy exposes state/trial callback counters, linear factorization/system-solve timing, accepted/rejected steps, and per-iteration telemetry |
| Legacy solve | RST-only fields remain `None`; existing Legacy callback timing is unchanged |
| production history | RST iteration history remains disabled; only statistics needed for the common attempt report are retained |
| benchmark integration | both preparation and worker dashboards compile with the detailed RST child rows |

**Implementation result:** The adapter now opts into
`PreparationTelemetryMode::Collect` only when `EquilibriumTimingMode::Enabled`
is requested. It maps the RST records into `RstPreparationTimingReport`, keeps
the public aggregate duration intact, and forwards `SolveStatistics` plus its
`SolveAttemptStatistics` into `SolverEvaluationTiming`. The focused telemetry
tests and bench compilation pass.

**Performance follow-up:** No new release timing numbers are claimed by this
Story. The next release benchmark must print the detailed rows before any RST
preparation optimization is selected; old aggregate-only Story 114 numbers are
not a baseline for the new child-stage attribution.

### Story 117: PERF-RST-12 release attribution of RST preparation and solver work

**Name:** `PERF-RST-12 release dashboard explains RP-1311/40 RST preparation and native solver time`

**Question / hypothesis:** Does the new telemetry explain the large release
budget in RP-1311/40, and does the native RST solver itself dominate the work?
The relevant threshold is `10 ms`: sub-millisecond rows are recorded for
completeness but are not optimization targets.

**Commands:**

```powershell
$env:RST_PREPARATION_REPORT_ONLY='1'
cargo bench --profile release --bench equilibrium_rst_preparation `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  -- --noplot

$env:RST_WORKERS_FILTER='rp1311_40'
$env:RST_WORKERS_REPORT_REPETITIONS='3'
cargo bench --profile release --bench equilibrium_rst_workers `
  --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-prettyprint" `
  -- --noplot
```

**Release observations, 2026-10-06:**

| case / scope | aggregate | dominant child evidence |
| --- | ---: | --- |
| RP-1311/40 P,T elemental inventory preparation | `78.75--81.66 ms` | symbolic differentiation `56.47--57.91 ms`; Jacobian lambdification `7.95--11.41 ms`; input validation `6.36--6.89 ms`; prepared assembly `5.50--6.43 ms` |
| RP-1311/40 nested P,H trial-0 template preparation | `77.14--78.65 ms` | symbolic differentiation `56.14--57.33 ms`; Jacobian lambdification `7.66--8.21 ms`; input validation `6.26--6.62 ms`; prepared assembly `5.52--6.37 ms` |
| RP-1311/40 preparation residual | `0.25--0.30 ms` | unattributed time; not a bottleneck |

**Native solver evidence for RP-1311/40 P,T element inventory:**

| method | facade median | status | solver wall / callback split |
| --- | ---: | --- | --- |
| RST LM | `106.53 ms` | accepted | about `5.3--5.8 ms`; Jacobian `4.4--5.0 ms`, residual `0.74--0.76 ms`, overhead `0.08 ms` |
| RST Minpack LM | `105.06 ms` | accepted | about `5.9--6.6 ms`; Jacobian `4.6--5.2 ms`, residual `0.71--0.77 ms`, overhead `0.04--0.06 ms` |
| RST Backtracking LM | `232.35 ms` | rejected | two attempts; `49--67` residual calls per attempt and `710--901` rejected steps in native trace rows |
| RST Nielsen LM | `223.42 ms` | rejected | two attempts; Jacobian callbacks dominate each attempt at about `90%` of backend wall |
| RST Damped Newton | `200.89 ms` | rejected | two attempts; callback wall is small per attempt, but the method does not pass acceptance |

The native telemetry table confirms that linear algebra is not the large
budget: for accepted LM, factorization is about `0.04 ms` and system solve
about `0.006 ms`. The dominant accepted-solver operation is Jacobian callback
evaluation, but the dominant end-to-end cold cost remains RST symbolic
preparation, especially symbolic differentiation. Failed methods are expensive
because of repeated attempts and rejected steps, not because of preparation
attribution ambiguity.

**Decision:** The next optimization candidate is RST symbolic differentiation,
with Jacobian lambdification second. Do not optimize linear solves or generic
solver overhead yet. The current telemetry makes those paths measurable, but
their observed cost is far below the `10 ms` threshold on the accepted large
case.

### Story 118: TP-1907 equilibrium output parity across RST symbolic frontends

**Name:** `TP-1907 facade A/B compares ExprLegacy and AtomView equilibrium outputs`

**Question / hypothesis:** Callback parity at one seed does not guarantee that
the nonlinear solve publishes the same equilibrium. With the same physical
problem, inventory, numerical seed, production cascade, and acceptance checks,
do ExprLegacy and AtomView return matching component fractions, temperature,
balance, and active topology? At the upper inventory limit, does each frontend
still have an RST backend that passes physical acceptance?

**Commands:**

```powershell
cargo test --lib --no-default-features `
  Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_facade_atom_view_and_expr_legacy_match_published_equilibrium `
  -- --nocapture

cargo test --lib --no-default-features `
  Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_upper_inventory_records_atom_view_acceptance_gap `
  -- --nocapture

cargo test --lib --no-default-features `
  Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_atom_view_and_expr_legacy_have_componentwise_callback_parity_at_the_seed `
  -- --nocapture
```

**Matched accepted-solution result, TP-1907 at 700 K, inventory scale `1e-3`:**

| measure | ExprLegacy | AtomView | difference |
| --- | ---: | ---: | ---: |
| accepted RST backend | LM | LM | same |
| temperature | `700.000000000 K` | `700.000000000 K` | `0 K` |
| max absolute component mole-fraction difference | - | - | `2.288253e-12` |
| max relative component-mole difference | - | - | `3.866215e-8` |
| max physical element-balance error | `5.020315e-7` | `5.020315e-7` | same |
| phase topology | gas + graphite | gas + graphite | same |

Both routes pass through the normal facade and production cascade using the
same elemental inventory and numerical seed. The per-component fractions are
printed by the test; the maximum absolute difference is far below the existing
`1e-8` fraction acceptance threshold.

**Upper-inventory acceptance result, scale `1e8`:**

The following table is the original RST 0.4.25 result; it is superseded by
the RST 0.4.26 debug retest in Story 119 below.

| frontend | RST backends passing physical acceptance | Legacy controls passing |
| --- | ---: | ---: |
| ExprLegacy | `1` (Trust-Region LM) | `2` |
| AtomView | `0` | `2` |

The callback-only parity test at the same upper-scale seed reports maximum
relative residual drift `1.321652e-9` and maximum relative Jacobian drift
`1.298984e-21` (global-magnitude normalized). Despite that close local parity,
the full upper-inventory solve does not pass acceptance with any tested RST
strategy through AtomView. Since AtomView has no accepted equilibrium at this
scale, there are no valid final fractions to compare there. The matching
Legacy controls show that the shared fixture and facade remain solvable.

**Historical conclusion on RST 0.4.25:** For accepted solutions at the small
inventory, frontend output differences are negligible. At `1e8`, AtomView
callbacks were close at sampled states but no tested RST backend passed the
facade acceptance matrix. Do not loosen physical tolerances to make this
comparison green.

### Story 119: TP-1907 staged AtomView failure localization

**Name:** `TP-1907 staged frontend diagnostics separate callback parity from solver-path failure`

**Question / hypothesis:** At inventory scale `1e8`, does AtomView first diverge
while preparing/evaluating the symbolic callbacks, or only after a nonlinear
method takes trial steps? Compare identical prepared equations at progressively
more informative states, then run the same Trust-Region LM directly with RST
Detailed telemetry on both physical and normalized fixed-active formulations.

**Command:**

```powershell
cargo test --lib --no-default-features `
  tp1907_atom_view_failure_is_localized_by_stage_and_native_telemetry `
  -- --nocapture
```

The diagnostic is isolated in
`frozen_reference/nasa_tp1907_frontend_diagnostics.rs`. It compares callbacks
at the physical seed, midpoint to ExprLegacy's accepted facade state, that
accepted state, and the normalized seed. Direct solve traces include iteration
residual norms, accept/reject status, callback counts/durations, and
linear-operation telemetry; only the first and last five iterations print.

**Observed callback parity at scale `1e8`:**

| common evaluation point | residual L2, ExprLegacy / AtomView | max residual delta | max Jacobian delta | worst residual row |
| --- | ---: | ---: | ---: | --- |
| physical seed | `1.248644e2 / 1.248644e2` | `1.192093e-7` | `1.162098e-12` | 16 |
| midpoint to accepted ExprLegacy state | `1.626005e8 / 1.626005e8` | `1.192093e-7` | `1.225270e-12` | 16 |
| accepted ExprLegacy facade state | `1.250011e-6 / 1.364169e-6` | `1.192093e-7` | `1.153244e-12` | 16 |
| normalized seed | `8.960091e1 / 8.960091e1` | `1.302993e-10` | `1.162098e-12` | 2 |

The largest callback discrepancy is in residual row 16, within rows 13--17
used for element balances in the exported TP-1907 system. The sampled states
show no large callback or Jacobian mismatch.

**Direct Trust-Region LM replay:**

| formulation / frontend | termination | iterations | final residual L2 | key evidence |
| --- | --- | ---: | ---: | --- |
| physical / ExprLegacy | MaxIterations | 500 | `5.297682e1` | 915 residual and 415 Jacobian calls |
| physical / AtomView | callback error | 11 residual calls | - | the solver passed an all-NaN 18-vector to the callback |
| normalized / ExprLegacy | StepTooSmall | 32 | `8.583697e-4` | 61 residual and 29 Jacobian calls |
| normalized / AtomView | StepTooSmall | 53 | `8.583697e-4` | 83 residual and 30 Jacobian calls |

At the normalized direct-replay endpoints, residual vectors differ by at most
`3.504897e-7`; `y17` (`graphite::C(gr)`) is `-27.9621` for ExprLegacy and
`-157.4445` for AtomView. This is a large coordinate difference, but both
coordinates exponentiate to amounts effectively indistinguishable from zero
at the system scale. It is therefore not evidence by itself of a meaningful
published-composition difference or a distinct physical state. The direct
fixed-active replay remains a diagnostic slice, not the facade's full
normalization/recovery workflow; its `StepTooSmall` outcome does not contradict
the accepted ExprLegacy facade route. Frontend comparisons should prioritize
published mole fractions, temperature, balances, and topology, with a stated
floor for trace amounts, rather than raw log-moles alone.

**Captured first failing trial:** A test-only recording wrapper retained the
actual inputs to each RST callback. AtomView's first failing input has NaN in
all 18 coordinates; evaluating that same input with ExprLegacy also returns a
non-finite residual. The immediately preceding finite state has residual L2
`7.682054e1`; its graphite coordinate is at the lower log-mole bound,
`y17=-708.3964`. At that shared finite state both frontends report Jacobian
maximum magnitude `8.946209e8`, estimated condition number `2.071418e9`, and
gradient L2 `6.933567e6`. The finite residual/Jacobian callback values agree to
the previously observed precision. Therefore AtomView did not produce the
first NaN: a Trust-Region LM internal step/subproblem produced an all-NaN trial
state, and the callback error is where it becomes visible. The exact internal
operation that first makes the step non-finite (linear subproblem, damping, or
trial assembly) still needs RST-side telemetry.

**Conclusion and RST recommendation:** The first bad object is the Trust-Region
LM trial vector, not an AtomView residual value computed from finite inputs.
Add an explicit finite check immediately after the regularized linear solve,
after constructing the nonlinear step, and before residual evaluation of
`trial_x`. On failure, return a typed numerical-breakdown result containing
iteration, current `x`/residual norm, Jacobian/gradient norms, damping and
trust-region radius, linear-solve status, and the first non-finite step/trial
coordinate. Keep this opt-in under Detailed telemetry. Also inspect why the
method reaches this transition after the graphite coordinate saturates its
lower bound. Do not label it an AtomView callback defect unless a finite input
produces a frontend-only mismatch. Normalized fixed-active runs both stall at
`8.6e-4`; their raw graphite log-coordinates differ, but both imply amounts
effectively indistinguishable from zero at this scale. This is not evidence of
a meaningful published-composition difference. No solver tolerance or
physical acceptance criterion changed.

**RST 0.4.26 debug retest: production facade at `1e8`:**

```powershell
cargo test -j 1 --lib --no-default-features `
  tp1907_upper_inventory_facade_atom_view_and_expr_legacy_match_published_equilibrium `
  -- --nocapture --test-threads=1

cargo test -j 1 --lib --no-default-features `
  tp1907_upper_inventory_records_frontend_acceptance_parity `
  -- --nocapture --test-threads=1
```

The production facade now accepts both frontends with the same RST LM backend.
At the common physical seed, residual relative drift is `1.321652e-9` and
Jacobian relative drift is `2.942790e-24`. The public facade outputs agree:

| Metric | ExprLegacy | AtomView | Difference |
| --- | ---: | ---: | ---: |
| recovered T | `700 K` | `700 K` | `0 K` |
| maximum mole-fraction delta | - | - | `8.327e-16` |
| maximum component relative delta | - | - | `5.223e-13` |
| maximum physical balance error | `4.832e-5` | `4.834e-5` | both within `1e-12 * sum(b)` |
| phase topology | gas + graphite | gas + graphite | equal |

The isolated fixed-active backend matrix also accepts `rst-trust-region-lm`
for both frontends (scaled residuals `1.730e-14` ExprLegacy and `3.079e-14`
AtomView); all other tested RST methods retain their prior rejection status.
The old `tp1907_upper_inventory_records_atom_view_acceptance_gap` expectation
was stale and is replaced by an acceptance-parity regression.

This closes the **published production-facade AtomView-vs-ExprLegacy parity**
question at `1e8`, not every direct-solver diagnostic. The raw fixed-active
Trust-Region replay still produces an all-NaN physical AtomView trial; the
normalized direct replays still stop at residual `8.583697e-4`, with a large
graphite log-coordinate difference representing a negligible amount. That
low-level numerical edge case remains separately open and must not be confused
with the accepted facade result.

### Story 120: Public ExprLegacy versus AtomView facade benchmark

**Name:** `Public frontend selection is honored; AtomView parity and speed depend on the solve route`

**Question / hypothesis:** When selected through the public calculator builder,
does AtomView preserve the ExprLegacy published result and improve end-to-end
time on representative normal-scale cases? The `1e8` stress case remains a
separate numerical failure investigation.

**Command (optimized release profile):**

```powershell
$env:RST_FRONTEND_AB_ONLY='1'
cargo bench --bench equilibrium_rst_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-rst-frontends" `
  -- equilibrium_rst_frontend_ab --noplot
```

Each pair selects the frontend through the public facade builder and holds the
case, RST Minpack LM, seed, bounds, and acceptance policy fixed. Criterion's
central facade estimate and the separate one-shot diagnostic preparation
measurement are not additive. Accepted pairs are compared by recovered
temperature, phase count, and published mole fractions; gates are
`|delta T| <= 1e-3 K`, maximum absolute mole-fraction difference `<= 1e-8`, and
equal phase count. AtomView is exercised through the production route, not a
test-only scoped override.

**Historical release results (RST 0.4.25):** central Criterion estimate;
preparation is a separate single diagnostic observation. The STANJAN-15 P,H
failure below was fixed in RST 0.4.26 and is superseded by the retest below;
the other measurements remain historical performance data, not a current
0.4.26 baseline.

| normal-scale case | status | facade estimate ExprLegacy / AtomView | RST preparation ExprLegacy / AtomView | max `|delta x|` | max `|delta T|` |
| --- | --- | ---: | ---: | ---: | ---: |
| Small3 P,T molecular | both pass | `0.560 / 0.452 ms` | `0.310 / 0.420 ms` | `3.23e-14` | `0 K` |
| STANJAN-15 P,T element inventory | both pass | `10.059 / 8.533 ms` | `5.077 / 3.521 ms` | `3.67e-13` | `0 K` |
| RP-1311/40 P,T element inventory | both pass | `197.13 / 132.14 ms` | `100.032 / 37.059 ms` | `2.10e-12` | `0 K` |
| STANJAN-15 P,H molecular, monolithic | Expr passes; AtomView fails | `378.57 / failed` | `108.338 / failed` | n/a | n/a |
| RP-1311/40 P,H element inventory | both pass | `341.52 / 251.82 ms` | `n/a / n/a` | `2.07e-12` | `1.54e-9 K` |

On RST 0.4.25, the AtomView STANJAN-15 P,H attempt failed at normalized basin discovery:
Minpack LM reaches its 50-iteration limit (88 residual and 38 Jacobian calls);
the candidate residual norm is `1.5655`, above the `1e-6` acceptance
tolerance. This is a real route-level failure, not a parity pass. The
RP-1311/40 P,H pair passes output comparisons, but its preparation timing is
unavailable (`n/a`, not zero). Among P,T cases, RP-1311/40 is the strongest
measured gain: facade estimate `197.13 -> 132.14 ms`; the independent
preparation observation is `100.032 -> 37.059 ms`. The one-shot preparation
diagnostics and Criterion estimates are separate measurements and should not
be combined arithmetically.

**Measurement correction:** the earlier Story 120 capture claimed all five
pairs passed and described an internal scoped override. That selector did not
reach monolithic P,H graph construction, so the earlier STANJAN-15 P,H row
was not a valid frontend A/B comparison. The public-selector capture above
supersedes its measurements and conclusion. P,T results were also recaptured
through the public API.

**Historical conclusion on RST 0.4.25:** AtomView was output-equivalent on the
measured P,T cases and RP-1311/40 P,H, but failed the STANJAN-15 monolithic P,H
acceptance gate. These normal-scale observations did not explain the separate
upper-scale `1e8` failure from Story 119.

**Historical targeted failure trace (RST 0.4.25, optimized release):**

```powershell
$env:RST_FRONTEND_AB_ONLY='1'
$env:RST_FRONTEND_CASE='stanjan15_ph_molecular'
$env:RST_FRONTEND_REPORT_ONLY='1'
cargo bench --bench equilibrium_rst_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-rst-frontends" `
  -- equilibrium_rst_frontend_ab --noplot
```

This runs one common facade case for both frontends and prints detailed
per-iteration RST counters; it skips Criterion sampling. ExprLegacy was
accepted (`T=2499.999999638 K`). AtomView produced two rejected candidates in
the physical and normalized-recovery attempts:

| AtomView attempt | RST termination | iterations | accepted / rejected steps | candidate residual L2 | acceptance limit |
| --- | --- | ---: | ---: | ---: | ---: |
| physical | `MaxIterations` | 50 | `37 / 13` | `1.5655113` | `1e-6` |
| normalized recovery | `Converged` | 49 | `38 / 11` | `2.8657492e-2` | `1e-6` |

Most iterations accepted a trial; rejected steps were intermittent rather than
a sustained rejection storm. The second attempt's backend termination label
does not mean the equilibrium candidate passed: the facade correctly applies
its stricter common acceptance gate and rejects the residual. Detailed RST
telemetry exposes per-iteration callback counts, accepted/rejected steps, and
stage durations, but not residual norms or iterate coordinates. Therefore
this trace confirmed a convergence/acceptance gap on RST 0.4.25. The callback
comparison and retest below localize it to the AtomView small-constant defect;
no KiThe tolerance or solver setting was changed.

**RST 0.4.26 retest and regression gates:**

```powershell
$env:RST_FRONTEND_AB_ONLY='1'
$env:RST_FRONTEND_CASE='stanjan15_ph_molecular'
$env:RST_FRONTEND_REPORT_ONLY='1'
cargo bench -j 1 --bench equilibrium_rst_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-rst-frontends" `
  -- equilibrium_rst_frontend_ab --noplot
```

The report-only production-facade A/B now accepts both frontends with the same
RST Minpack LM, seed, bounds, and unchanged acceptance gate:

| frontend | status | one-shot facade time | RST preparation | recovered T | phase count | max mole-fraction delta vs ExprLegacy |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| ExprLegacy | accepted | `349.250 ms` | `93.762 ms` | `2499.999999638 K` | 1 | reference |
| AtomView | accepted | `130.454 ms` | `40.974 ms` | `2499.999999638 K` | 1 | `1.804112e-16` |

These are the latest one-shot diagnostic timings from a repeated release
capture, not Criterion estimates. An earlier release capture measured
`330.645 / 127.589 ms` facade time and `95.891 / 39.594 ms` preparation for
ExprLegacy / AtomView respectively; both captures agree on acceptance and
published output. The normal
`lib` regression `stanjan15_tiny_constant_and_public_ph_facade_are_frontend_equivalent`
asserts that the parsed coefficient `-5.09024475e-16` and the amplified term
`-5.09024475e-16 * T^4` at `T=2500` survive AtomView Lambdification and agree
with ExprLegacy within eight `f64` epsilons. It also gates STANJAN initial
residual/Jacobian parity at `1e-12` and requires both public P,H facade routes
to be accepted, have matching phase topology, recover temperature within
`1e-3 K`, and match mole fractions within `1e-8`. The benchmark now separately
asserts that both frontends are accepted for this case.

On the full STANJAN symbolic callbacks, the initial-state max differences are
`1.42e-14` for residual and `1.22e-15` for Jacobian; at the AtomView diagnostic
endpoint they are `2.32e-14` and `8.53e-13`. The old low-level iteration-trace
helper in `stanjan15_atom_view_diagnostics.rs` uses a manually assembled
`solve_resolved_ph` request and its own 50-iteration path; it is retained only
for characterization and is not the public-facade acceptance gate. The
production-facade A/B and the new regression above are the authoritative
checks for the former failure.

**Conclusion on RST 0.4.26:** the normal-scale STANJAN-15 P,H AtomView failure
is resolved. The prior residual/Jacobian discrepancy came from AtomView
dropping very small nonzero constants; the minimal probe and STANJAN callback
parity now pass. The production-facade AtomView-vs-ExprLegacy comparison at
`1e8` is covered by the debug retest in Story 119: both frontends are accepted
and their published outputs agree within the stated tolerances. The isolated
raw Trust-Region replay remains an open diagnostic edge case, not a
production-facade acceptance failure.

**Current three-pass release refresh (RST 0.4.26, 2026-10-20):**

```powershell
$env:RST_FRONTEND_AB_ONLY='1'
cargo bench --bench equilibrium_rst_preparation --no-default-features `
  --features "equilibrium-instrumentation equilibrium-benchmark-rst-frontends" `
  -- --noplot
```

The release bench uses three explicitly separated evidence passes:

1. End-to-end public-facade Criterion A/B: fixed RST Minpack LM, same facade
   inputs, seeds, bounds, and acceptance policy; only ExprLegacy vs AtomView
   changes. A Criterion estimate is emitted only if both routes are accepted
   and their published outputs pass the parity gate.
2. One-shot preparation attribution: KiThe equation/formulation stages plus
   RST symbolic differentiation, residual/Jacobian lambdification, validation,
   and other reported preparation children. Aggregate and child times are
   nested and must not be added. Missing telemetry is printed as unavailable.
3. KiThe-native Legacy LM control on the same six facade requests. This is a
   whole-route control, not an algorithm-identical RST comparison; it does not
   attribute time to RST symbolic preparation.

Criterion used 10 samples, 1 s warmup, and a 2 s target measurement per entry.
The values below are estimates with Criterion confidence intervals. The
native column is a different solver implementation, so interpret it as a
facade-route reference rather than a frontend-only speed comparison.

| Case | RST Minpack + ExprLegacy | RST Minpack + AtomView | KiThe Legacy LM | AtomView reduction vs ExprLegacy |
| --- | ---: | ---: | ---: | ---: |
| Small3 P,T molecular | `0.549 ms [0.540, 0.558]` | `0.479 ms [0.473, 0.488]` | `0.400 ms [0.396, 0.405]` | `12.8%` |
| STANJAN-15 P,T ElementInventory | `11.579 ms [11.465, 11.705]` | `10.315 ms [10.170, 10.495]` | `6.355 ms [6.298, 6.414]` | `10.9%` |
| RP-1311/40 P,T ElementInventory | `216.614 ms [215.246, 218.078]` | `157.658 ms [151.384, 164.908]` | `106.568 ms [104.165, 109.652]` | `27.2%` |
| STANJAN-15 monolithic P,H molecular | `408.353 ms [403.297, 413.563]` | `155.899 ms [147.352, 166.798]` | `6.396 ms [6.346, 6.451]` | `61.8%` |
| RP-1311/40 nested P,H ElementInventory | `395.928 ms [392.006, 400.057]` | `290.647 ms [288.158, 293.960]` | `116.335 ms [114.189, 119.419]` | `26.6%` |

All five measured frontend pairs passed the same output gate: equal phase
count, `|delta T| <= 1e-3 K`, and maximum absolute mole-fraction delta
`<= 1e-8`. Report-only captures showed matching temperatures and maximum
mole-fraction deltas between `1.11e-16` and `1.80e-16` for these pairs.

**Unpaired 18-component diagnostic:** TP-1907/18 P,T ElementInventory is not
a performance result. RST Minpack + ExprLegacy failed acceptance: physical
candidate residual `5.1315e-3` exceeded `1e-6`, and normalized recovery
residual `3.7149e-4` exceeded `1e-12`. RST Minpack + AtomView and KiThe Legacy
LM were accepted, but there is no ExprLegacy reference with which to assert
output parity. The final run therefore emitted no Criterion estimate for
this case. A one-sided AtomView estimate from an interrupted exploratory run
is explicitly excluded.

**Preparation-stage evidence (separate one-shot release observations):**

| Case / stage, ms (ExprLegacy -> AtomView) | ExprLegacy | AtomView | Change |
| --- | ---: | ---: | ---: |
| STANJAN-15 P,T RST preparation aggregate | `3.669` | `2.575` | `-29.8%` |
| STANJAN-15 P,T symbolic differentiation | `2.531` | `0.694` | `-72.6%` |
| STANJAN-15 P,T residual / Jacobian lambdification | `0.056 / 0.305` | `0.077 / 0.293` | `+36.1% / -4.2%` |
| STANJAN-15 P,H RST preparation aggregate | `83.487` | `31.338` | `-62.5%` |
| STANJAN-15 P,H symbolic differentiation | `72.902` | `15.156` | `-79.2%` |
| STANJAN-15 P,H residual / Jacobian lambdification | `1.021 / 2.719` | `0.483 / 1.242` | `-52.7% / -54.3%` |
| RP-1311/40 P,T RST preparation aggregate | `84.928` | `31.037` | `-63.5%` |
| RP-1311/40 P,T symbolic differentiation | `61.359` | `11.885` | `-80.6%` |
| RP-1311/40 P,T residual / Jacobian lambdification | `0.686 / 9.033` | `0.589 / 5.691` | `-14.1% / -37.0%` |

RST input validation was effectively unchanged on the large P,T cases
(`6.319 / 6.222 ms` for RP-1311/40). RP-1311/40 P,H preparation attribution
was unavailable on both routes; its facade and Criterion measurements remain
valid, but no internal-stage conclusion is drawn. These single-shot stage
observations are not Criterion estimates and are not summed with the full
facade times.

**Conclusion:** AtomView is output-equivalent and faster than ExprLegacy on
the five matched workloads. The main measured preparation gain is symbolic
differentiation; the largest end-to-end gains are STANJAN-15 P,H and
RP-1311/40 P,T. This run does not establish a universal AtomView default:
Minpack's unpaired TP-1907/18 ExprLegacy failure needs a separate correctness
investigation, and the native LM figures use a different nonlinear method.
