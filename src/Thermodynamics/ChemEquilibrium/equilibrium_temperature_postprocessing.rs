//! Postprocessing for temperature-range equilibrium sweeps.
//!
//! # Purpose
//!
//! This module provides **postprocessing** for temperature-range equilibrium
//! calculations. After [`EquilibriumLogMoles::solve_for_T_range`](super::equilibrium_log_moles::EquilibriumLogMoles::solve_for_T_range)
//! produces a series of solutions at discrete temperature points, this module
//! can:
//!
//! - Preserve the raw solved points as-is.
//! - Build a **smoother render/export grid** via interpolation (PCHIP).
//! - Apply logarithmic interpolation for strictly positive quantities.
//! - Generate formatted tables for display or export.
//!
//! The module is intentionally **separate from the solver** — it takes solved
//! temperature series as immutable input and produces output without modifying
//! the solver state.
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`TemperaturePostprocessingPolicy`] | Controls resampling and interpolation strategy |
//! | [`TemperatureResamplingGrid`] | How the output grid is constructed |
//! | [`TemperatureInterpolationPolicy`] | Per-value interpolation (linear or log) |
//! | [`TemperatureInterpolationSpace`] | Enum: Linear or Log interpolation |
//! | [`TemperaturePostprocessingResult`] | Output container with interpolated values |
//!
//! # Key Functions
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`postprocess_temperature_series`] | Main entry point — processes a temperature series |
//!
//! # Dataflow
//!
//! ```text
//!   EquilibriumLogMoles::solve_for_T_range() produces:
//!     temperatures: Vec<f64>
//!     solution_values: Vec<Vec<f64>>  (one per species)
//!     │
//!     v
//!   TemperaturePostprocessingPolicy
//!     ├── grid: TemperatureResamplingGrid
//!     │     ├── RawOnly → keep original points
//!     │     ├── Uniform { points } → build uniform grid
//!     │     └── Explicit(vec) → use user-provided grid
//!     │
//!     ├── interpolation: HashMap<String, TemperatureInterpolationPolicy>
//!     │     ├── "linear" → TemperatureInterpolationSpace::Linear
//!     │     └── "log" → TemperatureInterpolationSpace::Log
//!     │
//!     └── clamp: bool → clamp interpolated values to input range
//!     │
//!     v
//!   postprocess_temperature_series(temperatures, values, policy)
//!     ├── Build PCHIP interpolant for each species
//!     ├── Evaluate on output grid
//!     ├── Apply log/linear interpolation per policy
//!     └── Return TemperaturePostprocessingResult
//!     │
//!     v
//!   TemperaturePostprocessingResult
//!     ├── output_temperatures: Vec<f64>
//!     ├── output_values: Vec<Vec<f64>>
//!     └── species_labels: Vec<String>
//! ```
//!
//! # Examples
//!
//! ```rust, ignore
//! use equilibrium_temperature_postprocessing::*;
//!
//! let policy = TemperaturePostprocessingPolicy {
//!     grid: TemperatureResamplingGrid::Uniform { points: 100 },
//!     interpolation: HashMap::new(), // defaults to linear
//!     clamp: true,
//! };
//!
//! let result = postprocess_temperature_series(
//!     &temperatures, &values, &policy
//! ).unwrap();
//! ```
//!
//! # Non-obvious Details
//!
//! - Uses **PCHIP** (Piecewise Cubic Hermite Interpolating Polynomial) from
//!   RustedSciThe, which preserves monotonicity and avoids overshoot.
//! - **Log-space interpolation** (`TemperatureInterpolationSpace::Log`) is
//!   recommended for mole numbers and mole fractions, which are strictly
//!   positive and often span several orders of magnitude.
//! - The `clamp` option prevents extrapolation beyond the solved temperature
//!   range, which would be physically meaningless.
//!
//! # Related Modules
//!
//! - [`equilibrium_log_moles`](super::equilibrium_log_moles) — solver that produces the raw series
//! - [`equilibrium_workflows`](super::equilibrium_workflows) — convenience wrappers for T-range solves
//!

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use RustedSciThe::numerical::optimization::inter_n_extrapolate::{
    InterpolationSpace as PchipSpace, Pchip,
};
use prettytable::{Cell, Row, Table};
use std::fmt;

/// How the output temperature grid should be constructed.
#[derive(Debug, Clone, PartialEq)]
pub enum TemperatureResamplingGrid {
    /// Keep the original solved points only.
    RawOnly,
    /// Build a uniform grid with the requested number of points.
    Uniform { points: usize },
    /// Use an explicit user-provided temperature grid.
    Explicit(Vec<f64>),
}

/// Interpolation policy for one postprocessing pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TemperatureInterpolationSpace {
    /// Interpolate the values directly.
    Linear,
    /// Interpolate `ln(value)` and exponentiate the result.
    ///
    /// This is appropriate for strictly positive quantities such as mole
    /// numbers and mole fractions.
    Log,
}

/// Value interpolation controls for one postprocessing pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TemperatureInterpolationPolicy {
    /// Interpolation space: `Linear` (direct PCHIP) or `Log` (PCHIP in log-space, then exponentiate).
    pub space: TemperatureInterpolationSpace,
    /// Whether to clamp interpolated values to the range of the original data.
    /// Prevents overshoot artifacts from PCHIP at sharp transitions.
    pub clamp: bool,
}

impl Default for TemperatureInterpolationPolicy {
    fn default() -> Self {
        Self {
            space: TemperatureInterpolationSpace::Linear,
            clamp: true,
        }
    }
}

/// Typed postprocessing policy for one temperature sweep.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperaturePostprocessingPolicy {
    pub grid: TemperatureResamplingGrid,
    pub interpolation: TemperatureInterpolationPolicy,
}

impl Default for TemperaturePostprocessingPolicy {
    fn default() -> Self {
        Self {
            grid: TemperatureResamplingGrid::RawOnly,
            interpolation: TemperatureInterpolationPolicy::default(),
        }
    }
}

/// One solved temperature sweep in row-major form.
///
/// `rows[i][j]` is the value of the `j`-th labelled series at temperature
/// `temperatures[i]`.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureSweepSeries {
    /// Column labels for each series (e.g., species names or "T").
    labels: Vec<String>,
    /// Temperature points in Kelvin, strictly increasing.
    temperatures: Vec<f64>,
    /// Data rows: `rows[i]` contains values at `temperatures[i]` for each label.
    rows: Vec<Vec<f64>>,
}

impl TemperatureSweepSeries {
    /// Builds one validated raw temperature series from solver rows.
    pub fn from_rows(
        labels: Vec<String>,
        rows: &[(f64, Vec<f64>)],
    ) -> Result<Self, ReactionExtentError> {
        if labels.is_empty() {
            return Err(invalid_series(
                "at least one series label is required for temperature postprocessing",
            ));
        }
        if labels.iter().any(|label| label.trim().is_empty()) {
            return Err(invalid_series("series labels must not be empty"));
        }
        if rows.is_empty() {
            return Err(invalid_series(
                "at least one solved temperature point is required",
            ));
        }

        let mut temperatures = Vec::with_capacity(rows.len());
        let mut previous_temperature = f64::NEG_INFINITY;

        for (index, (temperature, values)) in rows.iter().enumerate() {
            if !temperature.is_finite() {
                return Err(invalid_series(format!(
                    "temperature[{index}] must be finite"
                )));
            }
            if *temperature <= previous_temperature {
                return Err(invalid_series(
                    "temperature grid must be strictly increasing for resampling",
                ));
            }
            if values.len() != labels.len() {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "temperature row {index} has {} values but {} series labels were supplied",
                    values.len(),
                    labels.len()
                )));
            }
            for (column, value) in values.iter().enumerate() {
                if !value.is_finite() {
                    return Err(invalid_series(format!(
                        "series value at row {index}, column {column} must be finite"
                    )));
                }
            }
            temperatures.push(*temperature);
            previous_temperature = *temperature;
        }

        Ok(Self {
            labels,
            temperatures,
            rows: rows.iter().map(|(_, values)| values.clone()).collect(),
        })
    }

    /// Creates one raw series from an accepted temperature sweep.
    pub fn from_temperature_rows(
        labels: Vec<String>,
        rows: &[(f64, Vec<f64>)],
    ) -> Result<Self, ReactionExtentError> {
        Self::from_rows(labels, rows)
    }

    /// Number of sampled temperatures.
    pub fn point_count(&self) -> usize {
        self.temperatures.len()
    }

    /// Number of tracked series.
    pub fn series_count(&self) -> usize {
        self.labels.len()
    }

    /// Borrow the temperature grid.
    pub fn temperatures(&self) -> &[f64] {
        &self.temperatures
    }

    /// Borrow the series labels.
    pub fn labels(&self) -> &[String] {
        &self.labels
    }

    /// Borrow the raw row-major data.
    pub fn rows(&self) -> &[Vec<f64>] {
        &self.rows
    }

    /// Returns the data as stable labeled columns.
    pub fn columns(&self) -> Vec<(String, Vec<f64>)> {
        self.labels
            .iter()
            .enumerate()
            .map(|(column, label)| {
                (
                    label.clone(),
                    self.rows.iter().map(|row| row[column]).collect(),
                )
            })
            .collect()
    }

    /// Builds a stable summary for preview panes and logs.
    pub fn summary_rows(&self) -> Vec<TemperaturePostprocessingRow> {
        let mut rows = vec![
            TemperaturePostprocessingRow {
                section: "temperature_series",
                label: "point_count".to_string(),
                value: self.point_count().to_string(),
            },
            TemperaturePostprocessingRow {
                section: "temperature_series",
                label: "series_count".to_string(),
                value: self.series_count().to_string(),
            },
            TemperaturePostprocessingRow {
                section: "temperature_series",
                label: "temperature_min".to_string(),
                value: format!("{:.6}", self.temperatures.first().copied().unwrap_or(0.0)),
            },
            TemperaturePostprocessingRow {
                section: "temperature_series",
                label: "temperature_max".to_string(),
                value: format!("{:.6}", self.temperatures.last().copied().unwrap_or(0.0)),
            },
        ];

        for label in &self.labels {
            rows.push(TemperaturePostprocessingRow {
                section: "series",
                label: label.clone(),
                value: "available".to_string(),
            });
        }

        rows
    }

    fn validate_resampling_targets(
        &self,
        grid: &TemperatureResamplingGrid,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        match grid {
            TemperatureResamplingGrid::RawOnly => Ok(self.temperatures.clone()),
            TemperatureResamplingGrid::Uniform { points } => {
                if *points < 2 {
                    return Err(invalid_series(
                        "uniform resampling requires at least two points",
                    ));
                }
                if self.temperatures.len() < 2 {
                    return Err(invalid_series(
                        "resampling requires at least two source temperatures",
                    ));
                }
                let start = self.temperatures[0];
                let end = *self.temperatures.last().unwrap();
                let step = (end - start) / (*points as f64 - 1.0);
                Ok((0..*points)
                    .map(|index| start + step * index as f64)
                    .collect())
            }
            TemperatureResamplingGrid::Explicit(grid) => {
                if grid.is_empty() {
                    return Err(invalid_series(
                        "explicit temperature resampling grid must not be empty",
                    ));
                }
                let mut previous = f64::NEG_INFINITY;
                for (index, value) in grid.iter().enumerate() {
                    if !value.is_finite() {
                        return Err(invalid_series(format!(
                            "explicit temperature grid entry {index} must be finite"
                        )));
                    }
                    if *value <= previous {
                        return Err(invalid_series(
                            "explicit temperature grid must be strictly increasing",
                        ));
                    }
                    previous = *value;
                }
                Ok(grid.clone())
            }
        }
    }

    /// Resamples every series onto the requested temperature grid.
    pub fn resample(
        &self,
        policy: &TemperaturePostprocessingPolicy,
    ) -> Result<Option<Self>, ReactionExtentError> {
        let new_temps = self.validate_resampling_targets(&policy.grid)?;
        if matches!(policy.grid, TemperatureResamplingGrid::RawOnly) {
            return Ok(None);
        }
        if new_temps == self.temperatures {
            return Ok(Some(self.clone()));
        }

        let pchip_space = match policy.interpolation.space {
            TemperatureInterpolationSpace::Linear => PchipSpace::Linear,
            TemperatureInterpolationSpace::Log => PchipSpace::Log,
        };

        let mut resampled_rows = Vec::with_capacity(new_temps.len());
        let columns = self.columns();
        let interpolators: Vec<Pchip> = columns
            .iter()
            .map(|(_, values)| validate_and_build_pchip(&self.temperatures, values, pchip_space))
            .collect::<Result<_, _>>()?;

        for temperature in &new_temps {
            let mut row = Vec::with_capacity(self.series_count());
            for interpolator in &interpolators {
                row.push(interpolator.eval(*temperature, policy.interpolation.clamp));
            }
            resampled_rows.push(row);
        }

        Ok(Some(Self {
            labels: self.labels.clone(),
            temperatures: new_temps,
            rows: resampled_rows,
        }))
    }
}

/// One stable summary row for preview and logs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TemperaturePostprocessingRow {
    pub section: &'static str,
    pub label: String,
    pub value: String,
}

impl fmt::Display for TemperaturePostprocessingRow {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {} = {}", self.section, self.label, self.value)
    }
}

/// Full postprocessing result: raw solved data plus optional resampled data.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperaturePostprocessingResult {
    pub raw: TemperatureSweepSeries,
    pub resampled: Option<TemperatureSweepSeries>,
}

impl TemperaturePostprocessingResult {
    /// Stable preview rows for CLI/GUI/debug display.
    pub fn summary_rows(&self) -> Vec<TemperaturePostprocessingRow> {
        let mut rows = self.raw.summary_rows();
        rows.push(TemperaturePostprocessingRow {
            section: "postprocessing",
            label: "resampled".to_string(),
            value: self.resampled.is_some().to_string(),
        });
        if let Some(resampled) = &self.resampled {
            rows.push(TemperaturePostprocessingRow {
                section: "postprocessing",
                label: "resampled_point_count".to_string(),
                value: resampled.point_count().to_string(),
            });
        }
        rows
    }

    /// Renders the report into a compact table suitable for logs and tests.
    pub fn render_table(&self) -> String {
        let mut table = Table::new();
        table.add_row(Row::new(vec![
            Cell::new("section"),
            Cell::new("label"),
            Cell::new("value"),
        ]));

        for row in self.summary_rows() {
            table.add_row(Row::new(vec![
                Cell::new(row.section),
                Cell::new(&row.label),
                Cell::new(&row.value),
            ]));
        }

        table.to_string()
    }
}

impl fmt::Display for TemperaturePostprocessingResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for row in self.summary_rows() {
            writeln!(f, "{row}")?;
        }
        Ok(())
    }
}

/// Converts solved temperature rows into a postprocessed result.
pub fn postprocess_temperature_series(
    labels: Vec<String>,
    rows: &[(f64, Vec<f64>)],
    policy: &TemperaturePostprocessingPolicy,
) -> Result<TemperaturePostprocessingResult, ReactionExtentError> {
    let raw = TemperatureSweepSeries::from_rows(labels, rows)?;
    let resampled = raw.resample(policy)?;
    Ok(TemperaturePostprocessingResult { raw, resampled })
}

fn validate_and_build_pchip(
    x: &[f64],
    y: &[f64],
    space: PchipSpace,
) -> Result<Pchip, ReactionExtentError> {
    if x.len() != y.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "interpolation has {} x-values and {} y-values",
            x.len(),
            y.len()
        )));
    }
    if x.len() < 2 {
        return Err(invalid_series(
            "interpolation requires at least two source temperatures",
        ));
    }
    for (index, temperature) in x.iter().enumerate() {
        if !temperature.is_finite() {
            return Err(invalid_series(format!(
                "source temperature {index} must be finite"
            )));
        }
        if index > 0 && *temperature <= x[index - 1] {
            return Err(invalid_series(
                "source temperatures must be strictly increasing",
            ));
        }
    }
    if matches!(space, PchipSpace::Log) && y.iter().any(|value| *value <= 0.0) {
        return Err(invalid_series(
            "log-space postprocessing requires strictly positive values",
        ));
    }
    if y.iter().any(|value| !value.is_finite()) {
        return Err(invalid_series("series values must be finite"));
    }

    Ok(Pchip::new(x, y, space))
}

fn invalid_series(message: impl Into<String>) -> ReactionExtentError {
    ReactionExtentError::InvalidProblem {
        field: "temperature_postprocessing",
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_rows() -> Vec<(f64, Vec<f64>)> {
        vec![
            (300.0, vec![1.0, 0.1]),
            (400.0, vec![2.0, 0.2]),
            (500.0, vec![4.0, 0.4]),
        ]
    }

    #[test]
    fn raw_series_preserves_row_and_column_order() {
        let series = TemperatureSweepSeries::from_rows(
            vec!["A".to_string(), "B".to_string()],
            &sample_rows(),
        )
        .unwrap();

        assert_eq!(series.point_count(), 3);
        assert_eq!(series.series_count(), 2);
        assert_eq!(series.temperatures(), &[300.0, 400.0, 500.0]);
        assert_eq!(series.rows()[1], vec![2.0, 0.2]);
        let columns = series.columns();
        assert_eq!(columns[0].0, "A");
        assert_eq!(columns[0].1, vec![1.0, 2.0, 4.0]);
        assert_eq!(columns[1].0, "B");
        assert_eq!(columns[1].1, vec![0.1, 0.2, 0.4]);
    }

    #[test]
    fn uniform_resampling_builds_a_strictly_increasing_grid() {
        let series = TemperatureSweepSeries::from_rows(
            vec!["A".to_string(), "B".to_string()],
            &sample_rows(),
        )
        .unwrap();
        let policy = TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::Uniform { points: 5 },
            interpolation: TemperatureInterpolationPolicy {
                space: TemperatureInterpolationSpace::Linear,
                clamp: true,
            },
        };
        let resampled = series.resample(&policy).unwrap().unwrap();

        assert_eq!(
            resampled.temperatures(),
            &[300.0, 350.0, 400.0, 450.0, 500.0]
        );
        assert_eq!(resampled.series_count(), 2);
        assert!(
            resampled
                .rows()
                .iter()
                .all(|row| row.iter().all(|v| v.is_finite()))
        );
    }

    #[test]
    fn log_space_resampling_keeps_positive_series_positive() {
        let series = TemperatureSweepSeries::from_rows(
            vec!["A".to_string()],
            &[
                (300.0, vec![1e-6]),
                (400.0, vec![1e-4]),
                (500.0, vec![1e-2]),
            ],
        )
        .unwrap();
        let policy = TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::Uniform { points: 9 },
            interpolation: TemperatureInterpolationPolicy {
                space: TemperatureInterpolationSpace::Log,
                clamp: true,
            },
        };
        let resampled = series.resample(&policy).unwrap().unwrap();

        assert!(resampled.rows().iter().all(|row| row[0] > 0.0));
        assert!(
            resampled
                .rows()
                .windows(2)
                .all(|window| window[0][0] <= window[1][0])
        );
    }

    #[test]
    fn explicit_grid_respects_requested_temperatures() {
        let series = TemperatureSweepSeries::from_rows(
            vec!["A".to_string(), "B".to_string()],
            &sample_rows(),
        )
        .unwrap();
        let policy = TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::Explicit(vec![300.0, 450.0, 500.0]),
            interpolation: TemperatureInterpolationPolicy::default(),
        };
        let resampled = series.resample(&policy).unwrap().unwrap();

        assert_eq!(resampled.temperatures(), &[300.0, 450.0, 500.0]);
        assert_eq!(resampled.rows().len(), 3);
    }

    #[test]
    fn raw_only_policy_keeps_resampling_optional() {
        let series = TemperatureSweepSeries::from_rows(
            vec!["A".to_string(), "B".to_string()],
            &sample_rows(),
        )
        .unwrap();
        let policy = TemperaturePostprocessingPolicy::default();

        assert!(series.resample(&policy).unwrap().is_none());
    }

    #[test]
    fn postprocess_temperature_series_exposes_raw_and_resampled_views() {
        let policy = TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::Uniform { points: 4 },
            interpolation: TemperatureInterpolationPolicy {
                space: TemperatureInterpolationSpace::Linear,
                clamp: true,
            },
        };
        let result = postprocess_temperature_series(
            vec!["A".to_string(), "B".to_string()],
            &sample_rows(),
            &policy,
        )
        .unwrap();

        assert_eq!(result.raw.point_count(), 3);
        assert_eq!(result.resampled.as_ref().unwrap().point_count(), 4);
        assert!(
            result
                .summary_rows()
                .iter()
                .any(|row| row.section == "postprocessing" && row.label == "resampled")
        );
    }

    #[test]
    fn postprocessing_report_renders_a_table_for_logs_and_reports() {
        let policy = TemperaturePostprocessingPolicy {
            grid: TemperatureResamplingGrid::Uniform { points: 4 },
            interpolation: TemperatureInterpolationPolicy {
                space: TemperatureInterpolationSpace::Linear,
                clamp: true,
            },
        };
        let result = postprocess_temperature_series(
            vec!["A".to_string()],
            &[(300.0, vec![1.0]), (500.0, vec![3.0])],
            &policy,
        )
        .unwrap();

        let rendered = result.render_table();
        assert!(rendered.contains("section"));
        assert!(rendered.contains("temperature_series"));
        assert!(rendered.contains("postprocessing"));
        assert!(rendered.contains("resampled_point_count"));
        assert!(format!("{result}").contains("resampled"));
    }
}
