# Phase-Control Benchmark Baseline

This is a characterization baseline for the retained fixed-`P,T` phase-control
workflow. It is deliberately not a performance assertion: timings depend on
the machine, compiler, BLAS configuration, and background load. The benchmark
must remain opt-in with:

```text
cargo test --lib Thermodynamics::ChemEquilibrium::equilibrium_workflow_tests::tests::phase_control_benchmark_sweep_across_larger_inventories -- --ignored --nocapture
```

## Run

Date: 2026-07-26  
Fixture: synthetic closed phase-control inventory  
Outer iterations: 1 for every case  
Transitions: 0 for every case  
Reaction count: 0 for every case

| Total species | Phases | Active species | Projection build | Full solve |
|---:|---:|---:|---:|---:|
| 20 | 4 | 20 | 2.1379 ms | 3.9443 ms |
| 50 | 5 | 50 | 42.1631 ms | 50.0253 ms |
| 100 | 10 | 100 | 525.2114 ms | 582.3572 ms |
| 200 | 20 | 200 | 8.8057523 s | 8.8991402 s |

## Release build run

The following second baseline was captured by the application in a release
build on 2026-07-26. It is kept alongside the earlier debug-oriented run
because the two profiles answer different questions and must not be mixed.

| Total species | Phases | Active species | Projection build | Full solve |
|---:|---:|---:|---:|---:|
| 20 | 4 | 20 | 67.2 us | 257.8 us |
| 50 | 5 | 50 | 507.4 us | 640.6 us |
| 100 | 10 | 100 | 5.378 ms | 6.0646 ms |
| 200 | 20 | 200 | 70.5061 ms | 67.0024 ms |

This release run also had one outer iteration, zero transitions, and zero
reaction degrees of freedom in every case. It is therefore a useful fixed-set
orchestration baseline, not a complete production temperature-range benchmark.
At the measured full-solve cost, 100 repeated fixed-set calls would be roughly
64 ms for 50 species and 606 ms for 100 species. A real temperature sweep must
measure coefficient extraction, Gibbs refresh, residual/Jacobian construction,
backend iterations, and any active-set transitions separately.

## Interpretation

The debug profile was dominated by `ActiveSetProjection::build`, especially at
200 components, while the release profile is much smaller and should be the
more relevant deployment baseline. The projection depends on layout and active
set, not on the temperature-dependent standard-state Gibbs values. A future
range workflow should therefore reuse it while the active set is unchanged and
rebuild it only after a phase transition. The fixture has no transitions and no
reaction degrees of freedom, so it measures projection and fixed-set
orchestration overhead rather than chemical stiffness or a realistic
phase-appearance workload.

## Real local NASA release scaling

Command:

```text
cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_element_limited_release_scaling_matrix -- --ignored --nocapture
```

Fixture: 20/50/100 real `NASA_gas` records selected from the local catalog by
the allowed element alphabet `{C,H,O}` (`SubsetOf` semantics), with exact
offline lookup and legacy NR as an explicit single backend.

| Species | Total | Nonlinear solve | Residual L2 | Max element balance |
|---:|---:|---:|---:|---:|
| 20 | 3.9209 ms | 141.9 us | 2.0497117e-7 | 1.8438041e-8 |
| 50 | 8.9967 ms | 490.9 us | 5.5178320e-7 | 1.0011114e-7 |
| 100 | 20.0289 ms | 2.4965 ms | 7.6152139e-10 | 3.8032251e-11 |

All three points passed finite-positive-mole, residual, conservation, and
byte-for-byte library immutability checks. These values are a target-machine
release characterization, not a universal performance guarantee.

## Real 100-species, 50-point typed temperature range

Command:

```text
cargo test --release --lib Thermodynamics::ChemEquilibrium::equilibrium_live_data_tests::live_100_species_50_point_temperature_range_backend_matrix -- --ignored --nocapture
```

Fixture: 100 real local `NASA_gas` records from the C/H/O-limited catalog,
50 points from 1000 K to 1490 K, and one isolated `Single` policy per backend.
The test checks positive finite moles, residuals, scale-aware elemental
conservation, continuation, complete formulation build/reuse accounting, and
byte-for-byte source-library immutability.

Latest release evidence captured on 2026-07-28:

| Backend | Status | Total | Mean | Min | Median | Worst | Builds | Reuses | Symbolic updates | Max residual | Max balance |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| RST LM | OK | 21.368 s | 244.824 ms | 32.767 ms | 33.654 ms | 10.502 s | 1 | 49 | 49 | 8.456e-7 | 8.326e-7 |
| RST Minpack LM | OK | 21.357 s | 223.265 ms | 25.261 ms | 25.775 ms | 9.806 s | 1 | 49 | 49 | 1.907e-7 | 7.116e-8 |
| RST Nielsen LM | FAILED | - | - | - | - | - | - | - | - | - | - |
| RST Trust-region LM | OK | 20.476 s | 229.525 ms | 24.872 ms | 25.529 ms | 10.143 s | 1 | 49 | 49 | 1.907e-7 | 7.116e-8 |
| RST Powell Dogleg | FAILED | - | - | - | - | - | - | - | - | - | - |
| RST Damped Newton | OK | 20.981 s | 242.131 ms | 24.375 ms | 24.774 ms | 10.820 s | 1 | 49 | 49 | 1.907e-7 | 7.116e-8 |
| Legacy LM | OK | 175.791 ms | 3.018 ms | 2.355 ms | 3.202 ms | 3.917 ms | 1 | 49 | 0 | 9.514e-13 | 9.460e-13 |
| Legacy NR | OK | 97.803 ms | 1.356 ms | 1.285 ms | 1.298 ms | 3.950 ms | 1 | 49 | 0 | 1.907e-7 | 7.116e-8 |
| Legacy TR | FAILED | - | - | - | - | - | - | - | - | - | - |

Six of nine backends completed the full range. Nielsen LM failed because its
symbolic residual produced `NaN/Inf`; Powell Dogleg returned a candidate with
residual `0.8886`, rejected by the common acceptance gate; legacy TR reached
its nonlinear iteration limit. These are distinct failure classes and must not
be collapsed into one generic solver failure. The successful RST methods also
show rare 9.8-10.8 s points that dominate their total range time. `legacy_nr`
was the fastest successful backend in this run; that observation requires more
fixtures before it can drive a production default.

The next performance pass should add a second benchmark with a stable real
thermochemical fixture and at least one active-set transition. It should report
projection build, Jacobian, nonlinear, outer-loop, and total timings separately
before changing the data structures.
