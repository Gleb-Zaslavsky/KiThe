## Typical Workflows

### One experiment

#### Import and inspect

1. Open **File Manager > Import Data** for a raw TXT/CSV file with a header row, or use **New Experiment** when you also want to enter experiment metadata. Select the file and bind time (`t`), mass/instrument signal (`m`), and temperature (`T`) when present. Check the units before importing.
2. For an isothermal experiment without a measured temperature column, enter the known temperature in the experiment metadata. A `T` column is not required. If the instrument recorded a nearly constant temperature, bind that column; you can calculate its average later.
3. Click **Import experiment** / **Select**. The standard `m vs t` plot and, when temperature is bound, `T vs t` plot are created. Inspect them before transforming the data.

#### Clean and prepare

1. Remove startup or shutdown portions that are outside the experiment. Select the unwanted interval and use **Cut selected**, or set a known boundary with **Cut before time**. These operations remove rows; they do not fill missing measurements.
2. Use **Move time to zero** when the experiment should start at `t = 0`. Check the plot after each row or time transformation.
3. If the instrument column is voltage, select it and use **From mV to mg**. Configure the calibration line first. This creates a calibrated mass column; keep the voltage source when you want to compare the calibration result with the original signal.
4. Calculate relative mass/conversion from the mass column you intend to analyze. Enter a nonzero baseline time in **Input value**: the mean mass up to that time is used as initial mass. A value of `0` is rejected because it would select only the first sample. Relative mass and conversion are dimensionless and use opposite direction conventions.
5. Apply filters and smoothers to the intended source column. Savitzky-Golay (SG) is useful for local polynomial smoothing; LSQ spline resampling is useful for a smooth curve on a chosen time grid. Prefer a separate output column to preserve the source for comparison. For SG and other window methods, inspect the actual time span covered by the point window; for splines, inspect knot spacing and the resampled time range.
6. In **Column Manager**, inspect sample values, units, active bindings, transformation provenance, and quality comparisons. Confirm the output column contains numeric values and smoothing has not erased the trend or distorted the result.
7. Save the prepared experiment as typed CSV. It preserves column metadata and bindings, unlike a raw text export. Reopen it when you need to verify the round trip.

#### Isothermal versus non-isothermal

- **Isothermal:** if there is no temperature trace, store the known temperature in experiment metadata. If a temperature trace was measured, bind it and optionally use **Average on column** over the full range to inspect its mean. Use isothermal kinetic methods; do not assign a heating rate to a constant-temperature run.
- **Non-isothermal:** bind the measured temperature trace and record the heating rate in K/min. Convert Celsius to Kelvin when necessary before kinetic calculations. Inspect `T vs t` and its units; a non-isothermal method requires a valid temperature history, not just a metadata setpoint.

### Series of experiments

1. Import each run as a separate experiment with a unique ID identifying its condition, for example `iso_500K` or `ramp_2Kmin`. Bind columns and record the isothermal temperature or heating rate for each file. Do not combine separate runs into one experiment just to create a series.
2. Process runs one at a time: trim invalid beginning/end regions, shift time if appropriate, calibrate voltage to mass if required, calculate relative mass/conversion, then apply the chosen smoothing or resampling. Keep column meaning and units consistent across runs being compared.
3. In **Column Manager**, inspect each experiment's bindings, samples, and provenance. Use compare-samples to compare selected columns across runs; plotting does not require forcing adaptive or resampled grids to equal lengths.
4. Save the typed series after processing. Before a series kinetic method, verify each experiment has the required inputs: isothermal temperature metadata for isothermal methods, or a bound temperature trace and heating-rate metadata for non-isothermal methods.

#### Isothermal and non-isothermal series

- **Isothermal series:** create one experiment per temperature; temperature metadata is sufficient and a `T` column is optional. Independently calibrate mass, calculate conversion, and smooth each run. Check IDs and temperatures, then run an isothermal series method.
- **Non-isothermal series:** create one experiment per heating rate, bind its temperature trace, and record the rate in K/min. Apply the same preprocessing sequence to each run, checking the valid conversion range and temperature units. Then run a non-isothermal series method.
- Keep regimes separate. A missing `T` trace is valid for an explicitly isothermal experiment with a known temperature, but it cannot replace the temperature program of a non-isothermal experiment.
