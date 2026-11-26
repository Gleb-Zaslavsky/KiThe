## Experimental Kinetics: Help

## File Manager
The "File Manager" menu provides tools for loading, saving, and managing experiment data.

# Import Data
Imports raw TXT/CSV files after column binding and typed single-experiment KiThe CSV files. A typed series can be opened as the current document or added to it. Import never sets the destination used by Save. The file browser searches names, sorts directories/files, and filters to CSV/TXT by default; choose All files for other extensions.

# Open Document
Opens a typed KiThe series CSV as the current document, replacing the current series. If the current document has unsaved changes, choose Save, Discard, or Cancel first. Raw measurement files are imported, not opened as the document.

# New experiment
Opens a dialog for creating a new experiment. The window consists of two panels: left – file browser, right – metadata.

In the file browser you can navigate disks and folders, select a data file. The selected file is shown in the "File path" field. In metadata you can specify:
- Experiment identifier (unique name)
- Heating rate (K/min) – for non‑isothermal experiments
- Isothermal temperature (K) – for isothermal ones
- Comment (free text)

Below are column binding blocks:
- "bind m" – mass column binding: enter column name and choose unit (default mV). The program will interpret that column as mass.
- "bind T" – temperature column binding (units: °C, K)
- "bind t" – time column binding (units: s, h)

After filling, click the "Select" button – data will be loaded, metadata saved, bindings applied.

# Append to experiment
Adds a second file to an existing experiment after binding its columns. The resulting experiment is added to the current series, or replaces the source experiment when that option is selected.

# Save / Save As
Save writes the complete current series to a typed KiThe CSV, including all experiments, columns, metadata, bindings, provenance, and operation history. The first Save asks for a path; later Saves reuse that document path. Save As chooses a new path for the complete document. Ctrl/Cmd+S runs Save.

# Export CSV / Export TXT
Exports the whole series, selected experiments, or one experiment with selected columns. Column choices show their classification, unit, and active binding. Export settings and paths are independent of the current document; exporting never changes what Save writes. CSV preserves KiThe metadata, while TXT is plain tab-delimited data for one experiment.

# Manage Data
Opens the data manager where you can view, edit, and delete loaded experiments.

# Manage Plot
Opens the plot configuration dialog: choose displayed curves, adjust colors, labels, and other visual parameters.

# Plot Editor
Opens the print-ready plot editor for journal/report figures. It does not change live plots or experiment data.

## Math
The "Math" menu contains tools for mathematical data processing: filtering, smoothing, interpolation.

# Hampel filter
Detects and handles outliers in data. Configurable window size, threshold (sigma), and replacement strategy (median, NaN, drop).

The filter window contains:
- Column selection (col_name) – dropdown list of columns of the currently selected plot.
- Parameter window – size of the sliding window for median and MAD calculation.
- Parameter sigma – threshold value for outlier detection.
- Strategy selection (ReplaceWithMedian, ReplaceWithNaN, Drop) – how to treat outliers.
- Field out_col (optional) – name of a new column for the result; if omitted, the result replaces the original column.
- Buttons "Apply" and "Close".

# Savitzky Golay
Applies polynomial smoothing to reduce noise. Parameters: window size, polynomial order, derivative.

The smoothing window contains:
- Column selection (col) – dropdown list of columns.
- Parameter window – window size (odd number).
- Parameter poly_order – polynomial degree.
- Parameter deriv – derivative order (0 – smoothing, 1 – first derivative, etc.).
- Parameter delta – step along the X axis for derivative calculation.
- Field out_col (optional) – name of a new column for the result.
- Buttons "Apply" and "Close".

# Rolling mean
Computes a moving average over a given window. Simplifies data, reducing high‑frequency noise.

The window contains:
- Column selection (col_name) – dropdown list of columns.
- Parameter window – window size.
- Field out_col (optional) – name of a new column.
- Buttons "Apply" and "Close".

# Splines
Interpolates data using splines (linear, cosine, Bézier). Allows resampling data on a new time grid.

The window contains:
- Field new_time_col – name of the new time column (default "t_spline").
- Parameter n_points – number of points in the new grid.
- Spline type selection (Linear, Cosine, Bezier).
- For Bezier an additional parameter bezier_tension (tension).
- Field out_col (optional) – column name for interpolated values.
- Buttons "Apply" and "Close".

# LSQ Splines
Approximates data with least‑squares B‑splines. Configurable spline degree, number of points, and internal knots.

The window contains:
- Parameter degree – B‑spline degree.
- Parameter n_points – number of points in the new grid.
- Parameter n_internal_knots – number of internal knots.
- Field out_col (optional) – column name for the result.
- Buttons "Apply" and "Close".

# LOWESS
Locally weighted smoothing (LOESS) for non‑parametric regression. Configurable data fraction, number of iterations, delta parameter.

The window contains:
- Time column selection (time_col) – dropdown list of columns used as independent variable.
- List of columns to smooth – checkboxes for each column (except time). For each selected column you can specify an output column name (optional).
- Parameter fraction – fraction of data used in each local regression window.
- Parameter iterations – number of robust smoothing iterations.
- Parameter delta – parameter that speeds up computation.
- Checkbox parallel – use parallel computation.
- Buttons "Apply" and "Close".

## Kinetic Methods
The "Kinetic Methods" menu provides a set of methods for thermal decomposition kinetics analysis.

# Isoconversional
Computes effective activation energy as a function of conversion degree. Available methods: Ozawa‑Flynn‑Wall (OFW), Kissinger‑Akahira‑Sunose (KAS), Friedman, etc. Allows selecting experiments for analysis (all or only selected). The result is saved as a new experiment.

The analysis window contains:
- Method selection (Method) – dropdown list of available methods.
- Checkbox "Use selected only" – if checked, you can pick specific experiments from a list.
- Field "Output ID" – name under which the result will be saved in the experiment list.
- Button "Calculate" – start calculation.
- Field "Status" – shows progress or errors.

# Kissinger
Determines activation energy from DTG peak temperatures. Uses the dependence of peak temperature on heating rate. Results are displayed in the window (Ea, R²).

The window contains:
- Result display fields: Ea (kJ/mol) and R².
- Field "Kinetic expression" for entering an expression (under development).
- Buttons "Calculate kinetic expression" and "Calculate pre‑exponential factor" (under development).
- Button "Calculate" – start calculation.
- Field "System messages" – execution status.

# Combined Kinetics Analysis
Jointly fits model parameters (reaction order n, m) and activation energy by minimizing residuals across all experiments. Allows setting search ranges for n and m, number of steps, conversion degree α range.

The window contains:
- Search parameters: n_min, n_max, n_steps; m_min, m_max, m_steps; alpha_min, alpha_max, refinement_steps.
- Result display fields: n, m, Ea (kJ/mol), R².
- Button "Calculate" – start analysis.
- Field "System messages" – status.

# Criado Master Curve
Builds a master curve to check how well data match a particular model function f(α). (Under development)

# Fit Model
Interactive fitting of a kinetic model to experimental data. (Under development)

# Is this a sublimation?
Analyses data for conformity with a sublimation model (zero order). Outputs a verdict with confidence (Yes/Maybe/No) and average values of Ea, k.

The window contains:
- Result display fields: Mean Ea (kJ/mol), Mean k, Verdict.
- Button "Calculate" – start analysis.
- Field "System messages" – status.

# Golden Pipeline
Automated pipeline for sequential application of several kinetic methods with report generation.

## Direct Problem
The "Direct Problem" menu solves direct kinetic problems: simulating curves from given parameters.

# Solve IVP
Solves an initial value problem (IVP) for a system of ordinary differential equations describing kinetics. Allows obtaining theoretical mass, rate, etc. curves. (Under development)

# Solve BVP
Solves a boundary value problem (BVP) for more complex models, e.g., diffusion processes. (Under development)

## Test Options
The "Test Options" menu is intended for testing and debugging the program.

# Run Tests
Executes built‑in unit tests to verify correct operation.

# Load testing data
Loads predefined test data for demonstration.

# Create synthetic data
Generates artificial TGA data with given parameters: mode (isothermal/non‑isothermal), kinetic model, noise, outliers. Useful for method validation. The configuration window lets you set number of points, time step, initial temperature, heating rates, model parameters (m0, k0, E, R), and add noise (Gaussian, uniform, drift) and outliers (spikes).

# Show and redact settings
Opens the program settings window: calibration line (coefficients k, b), number of plot points, logging level. Changes are saved to tgasettings.json.

# Show Logs
Displays the program event log (logs) for diagnostics.

# Help
Shows this reference guide.

## Quick Action Panel
The Quick Action Panel contains buttons for frequently used operations on selected data.

# Manage Plots
Opens the plot management dialog (same as the menu item). In the dialog you can create new plots by selecting experiment, X column, and Y column. For each plot you can define a unique configuration. Buttons "Add Plot" add new configurations, "Create All" creates all defined plots.

# Column Manager
Opens the data column management dialog. Allows renaming, deleting, creating new columns in the data. You can also change data types and units.

# Golden Pipeline
Quick access to the golden pipeline.

# Clear Selected
Clears the selection on the plot.

# Reset View
Restores the original plot bounds.

# Refresh
Forces a plot redraw.

# Cut before time
Removes data before the specified time (in the "Input value" field).

# Cut selected
Trims data according to the selected range along the X or Y axis (first select X/Y with the button).

# Zoom To Selection
Zooms the plot to the selected region.

# Move time to zero
Subtracts the minimum time from the time column so that time starts at zero.

# Delete Plot
Deletes the selected curve from the plot.

# Calculate relative mass
Recalculates mass into relative mass (fraction of initial mass). Requires specifying the initial mass in the "Input value" field.

# Calculate conversion
Computes the conversion degree α (0 to 1). Requires specifying initial and final mass.

# From mV to mg
Calibrates mass from voltage (mV) using calibration coefficients. Creates a new column with mass in mg.

# From s to h
Converts time from seconds to hours.

# From C° to K
Converts temperature from degrees Celsius to Kelvin.

# Average on column
Computes the average value of the selected column over the whole range or over the selected interval.

# Conversion rate
Computes the derivative of conversion degree with respect to time (dα/dt).

# Dimensionless mass rate
Computes the derivative of relative mass with respect to time.

# Mass rate
Computes the derivative of mass with respect to time (dm/dt).

# T rate
Computes the derivative of temperature with respect to time (dT/dt).

# Arithmetic operations
Buttons "Add", "Sub", "Mul", "Divide", "ln", "exp" apply the corresponding mathematical operation to the selected column using the value from the "Input value" field. First you need to select a column (X or Y) and optionally check "apply only to the selected chart".

## Typical Workflows

### One experiment

#### Import and inspect

1. Open **File Manager > Import Data** for a raw TXT/CSV file with a header row, or use **New Experiment** when you also want to enter experiment metadata. Select the file and bind the time column (`t`), mass/signal column (`m`), and temperature column (`T`) when present. Check the displayed units before importing.
2. For an isothermal experiment without a measured temperature column, enter its known temperature in the isothermal metadata. A temperature column is not required. If the instrument recorded a nearly constant temperature, bind that column instead; you can calculate its average later.
3. Click **Import experiment** / **Select**. The standard `m vs t` and, when temperature is bound, `T vs t` plots are created. Inspect the curves before transforming data.

#### Clean and prepare

1. Remove startup or shutdown portions that are not part of the experiment. Select the unwanted plot interval and use **Cut selected**, or use **Cut before time** for a known boundary. These operations remove rows from the experiment; they do not fill missing measurements.
2. Use **Move time to zero** when the experiment should start at `t = 0`. Check the plot again after each row or time transformation.
3. If the instrument column is voltage, select that signal and use **From mV to mg**. The calibration coefficients must already be configured. This creates a calibrated mass column; retain the voltage column when you want to compare the original signal with the calibration result.
4. Calculate conversion from the mass column you intend to analyze. Choose a nonzero baseline interval in the `Input value` field: the mean mass up to that time is used as the initial mass. A value of `0` is rejected because it would select only the first sample. Relative mass and conversion are dimensionless quantities; they have opposite direction conventions.
5. Apply filters and smoothers to the intended source column. Savitzky-Golay (SG) is useful for local polynomial smoothing; LSQ spline resampling is useful when a smooth curve on a chosen time grid is needed. Prefer a separate output column so the source remains available for comparison. For SG and rolling-window methods, inspect the time span covered by the chosen point window; for splines, inspect knot spacing and the resampled time range.
6. Open **Column Manager** and inspect the sample values, units, active bindings, transformation provenance, and quality comparison. Confirm that the new column contains numeric values and that the smoothing did not erase meaningful trend or distort the result.
7. Save the prepared experiment as the typed CSV format. It preserves column metadata and bindings, unlike a raw text export. Reopen the saved file when you need to verify the round trip.

#### Isothermal versus non-isothermal

- **Isothermal:** set the known temperature in experiment metadata if no temperature trace exists. If a temperature trace was measured, bind it and optionally use **Average on column** over the selected whole range to inspect or record its representative value. Use isothermal kinetic methods; do not infer a heating rate from a constant-temperature experiment.
- **Non-isothermal:** bind the measured temperature trace and record the heating rate in K/min. Convert Celsius to Kelvin when necessary before kinetic calculations. Inspect `T vs t` and its units; non-isothermal kinetic methods require a valid temperature history, not only a metadata setpoint.

### Series of experiments

1. Import each run as a separate experiment with a unique ID that identifies its condition, for example `iso_500K` or `ramp_2Kmin`. Bind each file's columns and record its isothermal temperature or heating rate. Do not concatenate distinct runs into one experiment merely to make them look like a series.
2. Process runs one at a time: trim invalid beginning/end regions, shift time if appropriate, calibrate voltage to mass if required, calculate relative mass/conversion, then apply the chosen smoothing or resampling. Keep the same physical meaning, units, and transformation choices across runs when comparing them.
3. In **Column Manager**, inspect each experiment's bindings, samples, and provenance. Use the compare-samples view to compare selected columns across runs; different adaptive/resampled grids do not need to be forced to equal lengths for plotting.
4. Save the typed experiment series after processing. Before running a series kinetic method, verify that every experiment has the required inputs: isothermal temperature metadata for isothermal methods, or bound temperature data and heating-rate metadata for non-isothermal methods.

#### Isothermal and non-isothermal series

- **Isothermal series:** create one experiment per temperature, with the temperature stored in metadata (a `T` column is optional). Apply mass calibration, conversion, and smoothing independently to each run. Confirm IDs and temperatures, then run an isothermal series method.
- **Non-isothermal series:** create one experiment per heating rate, bind its temperature trace, and store the heating rate in K/min. Apply the same preprocessing sequence to each run, checking the valid conversion range and temperature units. Then run a non-isothermal series method.
- Keep the two regimes separate. A missing `T` trace is valid for an explicitly isothermal run, but it is not a substitute for a temperature program in a non-isothermal run.
