# Experimental Kinetics GUI Story Tests

The story suite lives in `tests/` and is registered by the GUI module under
`cfg(test)`. It drives the real `egui` app with `egui_kittest`; deterministic
CSV files are generated from `AdvancedTGAConfig`/`VirtualTGA`. The stories check
both user-visible workflow outcomes and domain state. Simulation stories print
lines prefixed with `[STORY]` so the numerical record can be extracted from
`cargo test -- --nocapture` output.

Run only these scenarios with:

```text
cargo test --lib gui::experimental_kinetics_gui::story_tests -- --nocapture --test-threads=1
```

## What The Filter Metrics Mean

Each filter report prints all five values from `FilterQualityReport`:

- `rough1` and `rough2` are filtered/raw standard-deviation ratios of first and
  second finite differences. A value below one means the respective difference
  signal became less rough; this alone does not prove that the physical signal
  became more accurate.
- `rmse`, `nrmse`, and `corr` compare the filtered result with its noisy input.
  They describe how much the filter changed that input, not its error against
  the unknown physical truth.
- `rmse_to_truth` is available because the tests have a paired noise-free
  synthetic run with the same kinetic parameters and seed. The `reduction`
  percentage compares noisy-input RMSE and filtered-output RMSE against that
  known signal. Positive is an improvement; negative is a degradation.

This distinction matters: a filter can reduce roughness while also biasing a
curve. The tests therefore require the mass Savitzky-Golay pass to reduce
noise-free-truth RMSE and require the final rate pass to reduce both roughness
ratios, but they report conversion-filter fidelity without claiming that every
second smoothing pass is beneficial.

## Three-Stage Quantitative Matrix

The `tests/matrix/` suite separates synthetic generation, GUI filter handling,
and GUI kinetic-method execution. It prints `tabled` tables, including the
fixed seed, so a run can be copied from `cargo test ... -- --nocapture` and
compared after changes. The matrix is diagnostic: it deliberately does not
assert that every filter or kinetic method must improve the answer.

### Generator Stage

`generation.rs` pairs a noisy series with the same-seed noise-free series for
isothermal and non-isothermal conditions. Profiles are Gaussian, uniform,
drift, impulse, and flicker-like noise. The generator table reports the seed,
number of runs and points, and mass RMSE to the paired clean signal; it also
checks deterministic output for equal seeds.

### GUI Filter Stage

`filters.rs` runs the selected-curve GUI model handlers for SG, Hampel/MAD,
rolling mean, LOWESS, linear spline resampling, and LSQ-spline resampling. Each
case reports all five `FilterQualityReport` metrics and RMSE against the
noise-free truth, on both thermal scenarios and all five noise profiles. The
same seed is used for noisy and clean runs; filter outputs and provenance are
checked.

`golden_pipeline.rs` clicks the front-panel Golden Pipeline action and the
dialog's `Calculate` button for both thermal scenarios, all five noise profiles,
and SG windows 11/31/51 (30 runs; 360 source points resampled to 240, then two
derivative-edge rows removed). It prints all five `FilterQualityReport` values
comparing final pipeline mass with the raw input interpolated to the output time grid,
plus mass RMSE against the paired noise-free truth before and after the
complete pipeline. The final
experiment, output bindings, finite metrics, and dialog close are checked.

On success the Golden Pipeline dialog now remains open and displays the SG
mass-filter provenance report: both roughness ratios, RMSE, normalized RMSE,
correlation, and the compared column names. It explicitly explains that these
compare the filtered signal with its input and do not establish truth accuracy.
On failure the dialog still displays the error text. The LSQ resampling stage
does not have a directly comparable same-grid quality report because it changes
the sampling grid; a visible explanation for that limitation remains a TODO.
The table below includes paired-truth RMSE from the test fixture, which is not
available for real experimental data and is not shown in the application.

Representative Golden Pipeline outcomes (seed 96000-series; the table printed
by the test contains all 30 cases):

| Case | SG window | `rough1` / `rough2` | Truth RMSE before → after (mg) | Change |
|---|---:|---:|---:|---:|
| Isothermal, Gaussian | 11 | 0.045 / 0.003 | 0.0971 → 0.0260 | 73.2% better |
| Isothermal, impulse | 31 | 0.033 / 0.003 | 0.1349 → 0.0423 | 68.6% better |
| Isothermal, drift | 11 | 1.000 / 9.316 | 0.0207 → 0.0207 | 0.1% worse |
| Isothermal, drift | 51 | 1.384 / 169.024 | 0.0207 → 0.0236 | 13.8% worse |
| Non-isothermal, Gaussian | 31 | 0.166 / 0.004 | 0.0878 → 0.0241 | 72.5% better |
| Non-isothermal, flicker | 31 | 0.422 / 0.019 | 0.2278 → 0.2214 | 2.8% better |

This is a useful caution, not a universal ranking: roughness ratios can be very
small while truth accuracy barely changes, and filtering a drift-dominated
signal with a wider SG window can make it worse. The full table and fixed seeds
are available with the matrix test's `--nocapture` output.

Example, isothermal Gaussian noise (seed 82000, 360 points): SG reduces truth
RMSE from 0.1241 to 0.0376 mg (69.7%), LOWESS to 0.0372 mg (70.0%), and LSQ
spline to 0.0291 mg (76.6%). Hampel changes this Gaussian trace only slightly
(1.0% reduction). Rolling mean lowers roughness ratios but increases truth
RMSE to 0.2339 mg (-87.7%), exposing phase/shape bias on the decaying trend.
Linear resampling reproduces the same signal on a new grid and has 0% truth
error reduction, as expected; it is a resampling operation, not a denoiser.
The matrix also shows that filter rankings depend strongly on noise type: a
method that is useful for impulses need not improve drift or flicker.

### GUI Kinetic Stage

`kinetics.rs` applies the selected mass filter through the GUI model handler,
then uses the visible Fit Model / Isoconversional windows and their Calculate
buttons. The isothermal branch fits DecExp independently at three temperatures,
regresses `ln(k)` against `1/T`, and also runs Friedman Integral. The
non-isothermal branch runs OFW, KAS, Starink, Differential Friedman, and
Vyazovkin; it additionally reports an independent first-order Arrhenius check
from the processed rate. Isoconversional methods estimate Ea, not the
pre-exponential factor A, so their A cells are intentionally `n/a`.

Observed representative results (true `Ea=85 kJ/mol`, `A=2e5 1/s`):

| Thermal/noise/filter | GUI method | R2 | Ea (kJ/mol) | Ea error | A (1/s) | A error |
|---|---|---:|---:|---:|---:|---:|
| Isothermal Gaussian / SG 31/3 | DecExp + Arrhenius | 1.00000 | 84.964 | 0.0% | 1.9848e5 | 0.8% |
| Isothermal Gaussian / SG 31/3 | Friedman Integral | 0.99972 | 80.022 | 5.9% | n/a | n/a |
| Isothermal impulse / Hampel 11/3 | DecExp + Arrhenius | 1.00000 | 84.995 | 0.0% | 1.9972e5 | 0.1% |
| Non-isothermal Gaussian / SG 31/3 | Differential Friedman | 0.99879 | 84.195 | 0.9% | n/a | n/a |
| Non-isothermal Gaussian / SG 31/3 | Vyazovkin | 0.00000 | 83.579 | 1.7% | n/a | n/a |
| Non-isothermal impulse / Hampel 11/3 | Differential Friedman | 0.99694 | 84.483 | 0.6% | n/a | n/a |
| Non-isothermal impulse / Hampel 11/3 | Vyazovkin | 0.00000 | 83.874 | 1.3% | n/a | n/a |
| Non-isothermal Gaussian / SG 31/3 | OFW / KAS / Starink | 0.99548 / 0.99477 / 0.99480 | 127.553 / 125.350 / 125.603 | 47.5-50.1% | n/a | n/a |

The last row is a notable scientific diagnostic: three integral methods give a
large common bias on this single-step synthetic fixture while Differential
Friedman and Vyazovkin return Ea close to the configured value. The integral
methods' high regression R2 does not remove that bias. Vyazovkin currently
exports `R2=0` in this result path; this is visible in the matrix and needs
separate interpretation rather than being presented as a high regression
quality. The table is evidence for follow-up validation of model assumptions,
temperature/rate units, and integral/Vyazovkin reporting; it is not a reason to
change expected values or weaken assertions without a scientific explanation.
Flicker-like noise also exposes parameter sensitivity: LOWESS recovers Ea
reasonably here but its fitted A is about 30-37% low in both thermal branches.

The matrix currently runs three representative noise/filter pairs in the
kinetic stage (Gaussian/SG, impulse/Hampel, flicker/LOWESS). The full five-noise
by six-filter coverage belongs to the filter stage, keeping the more expensive
kinetic-method matrix broad enough to compare methods but fast enough for a
regular test run.

## Scenarios And Fixtures

### Continue A Run From A Second File

Two CSVs each contain three samples with `t`, `T`, and `m`. The second file has
a reset clock and is appended through the GUI continuation dialog. The expected
inferred step is 1 s, the combined boundary is 3 s, and the result has six
aligned rows while the original experiment remains. An out-of-order preview is
rejected atomically. The test prints the row counts and outcome.

### Reconstruct A Missing Interval

Seed `7401`; 2,400 points at 1 s; one-step Arrhenius mass (`m0=100`,
`k0=2e5`, `E=82 kJ/mol`, `R=8.314`); non-isothermal fixture starts at 520 K
with beta `0.02 K/min`; Gaussian mass noise sigma `0.04 mg`. Samples 950..1149
are absent, not represented by nulls or zeros. The GUI fits a linear `T(t)` and
DecExp `m(t)` using points on both sides, evaluates 200 shared-grid points, and
commits a separate reconstructed experiment. The observed fit R2 values are
1.0000000 for temperature and 0.9999895 for mass. An intentionally unsuitable
linear mass fit reports R2 about 0.90266 and is not accepted; source rows remain
unchanged.

### Three Isothermal Runs

Seed `91743`; 600 samples at 2 s; temperatures 500, 530, and 560 K; one-step
Arrhenius truth `m0=100`, `k0=2e5`, `E=85 kJ/mol`; Gaussian mass noise sigma
`0.12 mg`, no temperature noise. The processing sequence is calibration,
Savitzky-Golay mass smoothing (`window=31`, `poly=3`), time conversion seconds
to hours, conversion with a 0.01 h baseline, SG conversion smoothing
(`31/3`), differentiation, joint removal of derivative-edge rows, and SG rate
cleanup (`11/3`). The point windows cover nominal spans of 60 s and 20 s at the
original 2 s sampling interval.

Observed mass RMSE against the noise-free signal falls from 0.124418 to
0.033645 mg (73.0%) at 500 K, 0.123824 to 0.037366 mg (69.8%) at 530 K, and
0.117424 to 0.052741 mg (55.1%) at 560 K. The rate cleanup lowers both
difference roughness ratios in all three runs. Conversion smoothing is less
uniform: its truth RMSE improves 6.2% at 500 K but worsens 7.1% and 37.9% at
530 K and 560 K. This is an explicit diagnostic result, not hidden as a pass.

GUI DecExp fits use the active smoothed mass and recover per-hour rate constants
with R2 above 0.9999. Arrhenius regression gives `E=84,961.37 J/mol` and
`A=1.983320e5`, against `85,000 J/mol` and `2.0e5`; the fit R2 is 1.000000.
GUI Friedman Integral runs on a clean CSV
round-trip; over `eta=0.1..0.9`, its mean Ea is `78,109.88 J/mol` (8.1% low),
which is currently within the story's 12% gate and remains visible in the
printed report.

### Three Non-Isothermal Runs

Seed `44921`; 360 samples at 2 s; start temperature 500 K; heating rates 2, 5,
and 10 K/min; the same `E=85 kJ/mol`, `A=2e5`, and mass model; Gaussian mass
noise sigma `0.03 mg`, no temperature noise. Processing uses SG mass (`31/3`),
SG conversion (`11/3`), central derivative with shared endpoint cleanup, then
SG rate cleanup (`11/3`).

Mass truth RMSE falls 46.7%, 9.0%, and 35.0% for 2, 5, and 10 K/min. The final
rate filter lowers both roughness ratios for every heating rate. In contrast,
the second SG pass on conversion slightly worsens truth RMSE in all three runs
(5.0%, 8.1%, and 10.1%); this supports retaining a user-selectable branch and
checking fidelity rather than treating smoother appearance as proof of
improvement.

GUI Differential Friedman returns mean Ea `84,573.58 J/mol` for
`eta=0.1..0.9` (0.5% low; individual values 80,847..87,080 J/mol). An
independent regression of `ln((deta/dt)/(1-eta))` against `1/T` gives
`E=84,998.78 J/mol`, `A=1.981528e5`, and R2 `0.996149`. The method is rerun
after exporting only `t`, `T`, `sg_mass`, `eta`, and `deta_dt` and importing the
CSV again. Missing heating-rate metadata is covered as an explicit failure.

## Interpreting The Current Results

On these fixtures, the broad mass SG pass has clear evidence of denoising and
preserves the kinetic truth well. A second SG pass over already-smoothed
conversion is not consistently helpful. The rate cleanup reduces finite-
difference roughness, but its report's RMSE is relative to the pre-cleanup rate,
not a direct physical-error estimate. These are fixture-specific measurements,
not universal recommendations for window sizes; the printed seed, sampling
interval, window, and truth comparison make future changes auditable.

### Integral Isoconversional Recovery Gate

The original short non-isothermal matrix fixture starts at 500 K, uses 360
samples at 2 s, and has heating rates 2, 5, and 10 K/min. Domain-only replay
with seed 94000 reproduces the GUI result: OFW reports 127.55 kJ/mol, KAS
125.35 kJ/mol, and Starink 125.60 kJ/mol for a truth of 85 kJ/mol. Removing
noise and all filtering changes these values negligibly. Differential Friedman
reports 84.20 kJ/mol and independent Vyazovkin reports 83.58 kJ/mol on the same
short data. Increasing the hot-start run to 5,000 samples reduces but does not
eliminate the integral-method bias (about 101-105 kJ/mol).

The short fixture's shared conversion range is limited by the slowest heating
run: its conversion reaches only about 0.265, while the fastest reaches 0.928.
The synthetic model also starts at 500 K with conversion reset to zero even
though its Arrhenius rate is already nonzero there. This is consistent with a
finite lower-temperature-bound error in the integral approximations, which
assume negligible reaction before the recorded range; it is not yet a general
validation of the approximation for every physical system. No formula was
changed based on this fixture.

The release gate uses a separate noise-free cold-start fixture (350 K, 3,000
samples at 2 s, same three heating rates and kinetic parameters), then exercises
the GUI's SG mass, time-unit, conversion, SG conversion, derivative, and rate
cleanup path. The gate requires all applicable non-isothermal methods to recover
the known Ea within 8% over the common 0.1-0.9 conversion interval. Observed
mean results: OFW 88.775, KAS 84.342, Starink 84.636, differential Friedman
84.734, and Vyazovkin 84.422 kJ/mol. Friedman Integral is excluded because it
requires isothermal data. The broader noisy matrix remains diagnostic: when
noise leaves conversion non-monotone, it records the explicit method
precondition failure instead of treating the numerical result as a recovery
pass.
