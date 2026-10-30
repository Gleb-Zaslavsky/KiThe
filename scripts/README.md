# Equilibrium campaign runner

`Invoke-EquilibriumCampaign.ps1` runs the maintained equilibrium story-test and
benchmark matrix as individual Cargo processes. It continues after a failed
step and keeps each step's Cargo output in its own log.

Run from the repository root in Windows PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-EquilibriumCampaign.ps1 -Suite smoke
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-EquilibriumCampaign.ps1 -Suite correctness -BuildJobs 1
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-EquilibriumCampaign.ps1 -Suite performance -BuildJobs 1
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-EquilibriumCampaign.ps1 -Suite all -BuildJobs 1 -FailOnAny
```

Use `-PlanOnly` to print the exact command queue without invoking Cargo or
creating report artifacts. `smoke` runs fast integrity/API checks;
`correctness` adds the selected backend, multiphase, GUI-range, and wide-range
story tests; `performance` runs the optimized Cargo bench profile with the
instrumentation features needed by each benchmark. `all` runs both groups.

Each executed campaign creates a unique directory under `test_reports/` with
`reports/campaign.json`, `reports/summary.md`, and per-step logs under
`technical/`. The summary is rendered with `tabled`; the JSON manifest retains
commands, revisions, toolchain, timings, statuses, and log paths. Story-test
documents are intentionally not modified by the runner.

## Kinetics tests

Run every `Kinetics::` library test, including ignored stress stories, while
excluding `experimental_kinetics`:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-KineticsTests.ps1
```

The runner uses one Cargo build job and one test thread by default. Pass
`-BuildJobs N` to change build parallelism; test execution remains serial. The
stress stories cover full Cantera payload parsing and recursive Cantera mechanism
construction.

## ReactorIVP release matrix

`Invoke-ReactorIvpCampaign.ps1` runs the IVP correctness stories and release
performance matrix sequentially. Each Cargo command has a separate raw log;
failures do not prevent later steps from running. The default `all` suite runs
the release story tests, gated AOT correctness (when `tcc` is installed), and
the Criterion plus Detailed benchmark with the LSODE2 AOT/Lambdify axis enabled.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorIvpCampaign.ps1 -PlanOnly
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorIvpCampaign.ps1 -Suite all -BuildJobs 1
```

Use `-Suite correctness` or `-Suite performance` to run one half of the
campaign, and `-FailOnAny` when the caller should receive a nonzero exit code
if a step fails. `-PlanOnly` prints the queue without invoking Cargo. Results
are stored under a unique `test_reports/reactor-ivp-*/` directory: the Tabled
summary includes test status, Criterion medians/intervals, and the Tabled stage,
prepared-lifecycle, native-telemetry, and AOT/Lambdify matrices. Full Cargo
output and renderer diagnostics remain under `technical/`.

## ReactorBVP release matrix

`Invoke-ReactorBvpCampaign.ps1` runs the physical compact-HMX stories, the
frontend/layout and enlarged-grid matrix, compiler-gated AOT correctness, and
the small/medium/large Criterion baseline. Damped and BVP_sci are included;
Frozen/Naive is intentionally omitted from comparative acceptance/performance
routes while its known physical boundary failure remains documented by its
negative correctness story.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorBvpCampaign.ps1 -PlanOnly
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorBvpCampaign.ps1 -Suite correctness -BuildJobs 1
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorBvpCampaign.ps1 -Suite performance -BuildJobs 1
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\Invoke-ReactorBvpCampaign.ps1 -Suite all -BuildJobs 1 -FailOnAny
```

The `all` campaign enables AOT comparisons if TCC is installed; otherwise it
records those steps as not applicable and still completes the Lambdify matrix.
Each run writes a unique `test_reports/reactor-bvp-*/` directory containing a
Tabled summary, JSON manifest, and separate raw logs under `technical/`.
Use `-PlanOnly` to inspect the command queue without running Cargo or creating
reports.
