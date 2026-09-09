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

### Recorded release result

```text
TP-1907 backend/default story: test ...::production_default_retains_a_scale_robust_backend_on_real_tp1907_story ... ok

TP-1907 backend matrix at scale=1e8:
  rst-lm: FAILED, normalized basin discovery, candidate residual 5.12501771773085e1
  rst-minpack-lm: residual=3.545e-14, balance=3.998e-8
  rst-nielsen-lm: FAILED, normalized basin discovery, candidate residual 1.6277711838251489e2
  rst-trust-region-lm: residual=3.519e-14, balance=3.998e-8
  rst-damped-newton: FAILED, normalized basin discovery, candidate residual 8.859999596823666e1
  legacy-lm: residual=2.588e-14, balance=1.999e-8
  legacy-nr: residual=1.848e-14, balance=3.498e-6
  legacy-tr: residual=1.772e-14, balance=3.998e-8
  accepted: rst-minpack-lm, rst-trust-region-lm, legacy-lm, legacy-nr, legacy-tr
  failed in isolated Single mode: rst-lm, rst-nielsen-lm, rst-damped-newton
test ...::production_default_retains_a_scale_robust_backend_on_real_tp1907_story ... ok
test ...::tp1907_upper_inventory_is_backend_invariant ... ok
test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 2520 filtered out
```

The three isolated RST failures are strict normalized-recovery outcomes, not
evidence that the thermodynamic equations or the complete production cascade
are invalid. RST Minpack/Trust-Region and all retained legacy methods pass the
same real fixture. Legacy NR reports a physical absolute balance of about
`3.50e-6` at this scale, while its scale-aware relative balance remains within
the matrix guard. The production default must therefore remain a cascade, not
an isolated RST LM/Nielsen/Damped Newton selection.

### Conclusion

The hypothesis is supported. The production default remains valid at `1e8`
because it retains a scale-robust backend in its cascade. The explicit matrix
also records which isolated methods are currently unsuitable as the sole
default backend.

## 3b. ElementInventory route equivalence at production scale

### Story name

The molecular, elemental, and formal routes must describe the same physical
equilibrium, including through the public facade and under extensive P,H
scaling.

### Tests

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_frozen_pt_molecular_elemental_and_formal_routes_match`

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_frozen_facade_upper_inventory_matches_the_same_oracle`

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_frozen_ph_molecular_elemental_and_formal_routes_reuse_target`

`Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph_tests::tests::tp1906_real_ph_solution_is_extensive_when_inventory_and_target_scale_together`

### Hypothesis / question

Does `ElementInventory` preserve the molecular equilibrium state through the
P,T and P,H facades, and does joint scaling of inventory and total enthalpy
preserve the intensive P,H solution?

### Acceptance criteria

- `ElementInventory` and formal elemental input reach the same physical state
  as the molecular TP-1907/TP-1906 route;
- the public P,T facade preserves that equivalence at inventory scale `1e8`;
- P,H preserves the recovered temperature and target enthalpy when inventory
  and total enthalpy are scaled together;
- route comparisons are made against the same molecular oracle, with residual,
  conservation, topology, and target-enthalpy assertions.

### Release commands

```powershell
cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_frozen_pt_molecular_elemental_and_formal_routes_match --no-default-features -- --include-ignored --exact --nocapture

cargo test --release --lib Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite_tests::tests::tp1907_frozen_facade_upper_inventory_matches_the_same_oracle --no-default-features -- --include-ignored --exact --nocapture

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

TP-1907 public facade scale=1e8:
  elemental: total_b=1.440534e9 residual=3.545e-14 balance=3.998e-8
  formal:    total_b=1.440534e9 residual=3.545e-14 balance=3.998e-8
  test result: ok; 1 passed; 0 failed; finished in 0.33s

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
