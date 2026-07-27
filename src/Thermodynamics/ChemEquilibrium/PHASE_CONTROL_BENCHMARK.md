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

The next performance pass should add a second benchmark with a stable real
thermochemical fixture and at least one active-set transition. It should report
projection build, Jacobian, nonlinear, outer-loop, and total timings separately
before changing the data structures.
