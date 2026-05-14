## Experimental Kinetics: Help

## File Manager
The "File Manager" menu provides tools for loading, saving, and managing experiment data.

# Import Data
Loads data from CSV, TXT, or Excel files. Various column formats are supported. After import, data appears in the table and on the plot.

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

# Append data to experiment
Adds new data to an existing experiment, merging them into a single dataset.

# Export Data
Saves current experiment data to a file of chosen format (CSV, TXT, Excel).

# Manage Data
Opens the data manager where you can view, edit, and delete loaded experiments.

# Manage Plot
Opens the plot configuration dialog: choose displayed curves, adjust colors, labels, and other visual parameters.

# Save As CSV/TXT/Excel/Image
Exports data in a specific format. "Save as image" saves the current plot to an image file (PNG, JPEG, etc.) via the KiThe Plot editor window.

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