//! Preprocessing adapters for experimental-kinetics datasets.
//!
//! KiThe owns column/schema validation, output naming, units, provenance, and
//! error context. RustedSciThe owns the numerical filter and spline primitives;
//! their behavior is guarded here by integration regression tests.

use crate::Kinetics::experimental_kinetics::column_provenance::ColumnTransformKind;
use crate::Kinetics::experimental_kinetics::math::ndarray_statistics::filter_quality_from_arrays;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    AffectedColumns, TGADataset, TGADomainError,
};
use RustedSciThe::numerical::data_processing::LSQSplines::{
    SolverKind, make_lsq_univariate_spline,
};
use RustedSciThe::numerical::data_processing::lowess_wrapper::{
    LowessConfig, lowess_smooth_values,
};
use RustedSciThe::numerical::data_processing::savgol2::SG_filter_dyn;
use RustedSciThe::numerical::data_processing::splines::{
    SplineKind, spline_eval, spline_keys, spline_resample, uniform_grid_from,
};
use log::{info, warn};
use ndarray::{Array1, s};
use polars::prelude::*;
use std::any::Any;
use std::collections::HashSet;
use std::f64;
use std::time::Instant;
/* TODO!:
SG как альтернатива derive_rate

unit-инференс для SG-deriv

SciPy-like SG (#2) как отдельный mode

lazy-реализация rolling_mean_many


*/
//====================================================================
//SPLINES
//===========================================================================
pub fn extract_f64_column(df: &DataFrame, col: &str) -> Result<Vec<f64>, TGADomainError> {
    let s = df
        .column(col)
        .map_err(|_| TGADomainError::ColumnNotFound(col.into()))?;

    let actual_type = s.dtype().to_string();
    let s = s
        .cast(&DataType::Float64)
        .map_err(|_| TGADomainError::NonNumericColumn {
            column: col.into(),
            actual_type,
        })?;
    let s = s.f64().map_err(TGADomainError::PolarsError)?;

    if s.null_count() > 0 {
        return Err(TGADomainError::InvalidOperation(format!(
            "Operation does not support nulls in column '{}'",
            col
        )));
    }

    Ok(s.into_no_null_iter().collect())
}

fn extract_optional_f64_column(
    df: &DataFrame,
    col: &str,
) -> Result<Vec<Option<f64>>, TGADomainError> {
    let series = df
        .column(col)
        .map_err(|_| TGADomainError::ColumnNotFound(col.into()))?;
    let actual_type = series.dtype().to_string();
    let casted = series
        .cast(&DataType::Float64)
        .map_err(|_| TGADomainError::NonNumericColumn {
            column: col.into(),
            actual_type,
        })?;
    let values = casted.f64().map_err(TGADomainError::PolarsError)?;
    Ok(values.iter().collect())
}
pub fn resample_all_columns(
    df: &DataFrame,
    time_col: &str,
    x_new: &[f64],
    kind: SplineKind,
) -> Result<DataFrame, TGADomainError> {
    validate_spline_kind(kind)?;
    if x_new.len() < 2
        || x_new.iter().any(|value| !value.is_finite())
        || x_new.windows(2).any(|pair| pair[1] <= pair[0])
    {
        return Err(TGADomainError::InvalidFilterParameters {
            filter: "spline resampling".into(),
            details:
                "explicit output grid must contain at least two finite, strictly increasing points"
                    .into(),
        });
    }

    let x_old = extract_f64_column(df, time_col)?;
    if x_old.len() < 2
        || x_old.iter().any(|value| !value.is_finite())
        || x_old.windows(2).any(|pair| pair[1] <= pair[0])
    {
        return Err(TGADomainError::InvalidColumnValues {
            column: time_col.into(),
            details: "spline resampling requires at least two finite, strictly increasing X values"
                .into(),
            indices: Vec::new(),
        });
    }
    if x_new[0] < x_old[0] || x_new[x_new.len() - 1] > x_old[x_old.len() - 1] {
        return Err(TGADomainError::OutOfRange);
    }

    let mut cols: Vec<Column> = Vec::with_capacity(df.width());

    for col in df.columns() {
        let name = col.name();
        let y_new = if name.as_str() == time_col {
            x_new.to_vec()
        } else {
            let y_old = extract_f64_column(df, name)?;
            if y_old.len() != x_old.len() {
                return Err(TGADomainError::ExperimentDataLengthMismatch {
                    experiment_id: "<dataframe>".into(),
                    column: name.to_string(),
                    expected: x_old.len(),
                    actual: y_old.len(),
                });
            }
            let keys = spline_keys(&x_old, &y_old, kind).map_err(|error| {
                contextual_backend_error(
                    "resample_all_columns",
                    name,
                    format!("spline construction failed: {error}"),
                )
            })?;
            let interior: Vec<f64> = x_new
                .iter()
                .copied()
                .filter(|value| *value > x_old[0] && *value < x_old[x_old.len() - 1])
                .collect();
            let interior_values = if interior.is_empty() {
                Vec::new()
            } else {
                spline_eval(keys, &interior).map_err(|error| {
                    contextual_backend_error(
                        "resample_all_columns",
                        name,
                        format!(
                            "spline evaluation failed (kind={kind:?}, points={}): {error}",
                            x_new.len()
                        ),
                    )
                })?
            };
            let mut interior_values = interior_values.into_iter();
            x_new
                .iter()
                .map(|value| {
                    if *value == x_old[0] {
                        Ok(y_old[0])
                    } else if *value == x_old[x_old.len() - 1] {
                        Ok(y_old[y_old.len() - 1])
                    } else {
                        interior_values.next().ok_or_else(|| {
                            contextual_backend_error(
                                "resample_all_columns",
                                name,
                                "spline returned fewer values than requested",
                            )
                        })
                    }
                })
                .collect::<Result<Vec<_>, _>>()?
        };
        validate_backend_output("resample_all_columns", name, &y_new, x_new.len())?;
        cols.push(Series::new(name.clone(), y_new).into());
    }

    Ok(DataFrame::new(x_new.len(), cols)?)
}

//=======================================================================
//SMOOTHING API
//=======================================================================

#[derive(Clone, Copy, Debug)]
pub enum SGMode {
    FirPadding,
    // SciPyLike, // intentionally not implemented yet
}
#[derive(Debug, Clone)]
pub enum SmoothStrategy {
    RollingMean {
        window: usize,
    },

    Hampel {
        window: usize,
        n_sigma: f64,
        strategy: HampelStrategy,
    },

    SavitzkyGolay {
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
        mode: SGMode,
    },
    // --- placeholders ---
    Lowess {
        frac: f64,
    },

    Spline {
        knots: usize,
    },
}
#[derive(Debug, Clone)]
pub enum HampelStrategy {
    ReplaceWithMedian,
    ReplaceWithNaN,
    Drop,
}

fn hampel_median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 0 {
        v[n / 2 - 1] / 2.0 + v[n / 2] / 2.0
    } else {
        v[n / 2]
    }
}
fn median_abs_deviation(v: &[f64], med: f64) -> f64 {
    let devs: Vec<f64> = v.iter().map(|x| (x - med).abs()).collect();
    hampel_median(devs)
}

fn default_output_name(prefix: &str, input_col: &str, output_col: Option<&str>) -> String {
    match output_col.map(str::trim) {
        Some(name) if !name.is_empty() => name.to_string(),
        _ => format!("{}{}", prefix, input_col),
    }
}

fn validate_filter_output(
    dataset: &TGADataset,
    input: &str,
    output: &str,
) -> Result<(), TGADomainError> {
    validate_output_mapping(dataset, &[input], &[output.to_string()])
}

fn validate_output_mapping(
    dataset: &TGADataset,
    inputs: &[&str],
    outputs: &[String],
) -> Result<(), TGADomainError> {
    if inputs.len() != outputs.len() {
        return Err(TGADomainError::InvalidOperation(format!(
            "Expected one output name per input column ({} inputs, {} outputs)",
            inputs.len(),
            outputs.len()
        )));
    }

    let frame_schema = dataset.frame.clone().collect_schema()?;
    let mut seen = HashSet::with_capacity(outputs.len());
    for (index, (&input, output)) in inputs.iter().zip(outputs).enumerate() {
        if output.trim().is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "Output column names cannot be empty".into(),
            ));
        }
        if !seen.insert(output.as_str()) {
            return Err(TGADomainError::ColumnAlreadyExists(output.clone()));
        }
        if inputs
            .iter()
            .enumerate()
            .any(|(other_index, other_input)| other_index != index && *other_input == output)
        {
            return Err(TGADomainError::ColumnAlreadyExists(output.clone()));
        }
        if output != input
            && (dataset.has_column_name(output) || frame_schema.get(output).is_some())
        {
            return Err(TGADomainError::ColumnAlreadyExists(output.clone()));
        }
    }
    Ok(())
}

fn validate_hampel_parameters(window: usize, n_sigma: f64) -> Result<(), TGADomainError> {
    if window < 3 || window % 2 == 0 || !n_sigma.is_finite() || n_sigma < 0.0 {
        return Err(TGADomainError::InvalidFilterParameters {
            filter: "Hampel".into(),
            details: "window must be odd and >= 3; n_sigma must be finite and >= 0".into(),
        });
    }
    Ok(())
}

fn validate_spline_kind(kind: SplineKind) -> Result<(), TGADomainError> {
    if matches!(kind, SplineKind::Bezier(tension) if !tension.is_finite()) {
        return Err(TGADomainError::InvalidFilterParameters {
            filter: "spline resampling".into(),
            details: "Bezier tension must be finite".into(),
        });
    }
    Ok(())
}

fn contextual_backend_error(
    operation: &str,
    column: &str,
    details: impl Into<String>,
) -> TGADomainError {
    TGADomainError::InvalidOperation(format!("RustedSciThe: {}", details.into())).with_context(
        operation,
        None,
        Some(column),
    )
}

fn validate_backend_output(
    operation: &str,
    column: &str,
    values: &[f64],
    expected_len: usize,
) -> Result<(), TGADomainError> {
    if values.len() != expected_len || values.iter().any(|value| !value.is_finite()) {
        return Err(TGADomainError::NoConvergence {
            method: operation.into(),
            details: format!(
                "backend returned {} values (expected {}) or non-finite output for column '{}'",
                values.len(),
                expected_len,
                column
            ),
        }
        .with_context(operation, None, Some(column)));
    }
    Ok(())
}

fn validate_lowess_inputs(
    operation: &str,
    column: &str,
    x_column: &str,
    x: &[f64],
    y: &[f64],
    config: &LowessConfig,
) -> Result<(), TGADomainError> {
    let invalid_parameters = |details: String| {
        TGADomainError::InvalidFilterParameters {
            filter: "LOWESS".into(),
            details: format!("{details} (X='{x_column}', Y='{column}')"),
        }
        .with_context(operation, None, Some(column))
    };

    config.validate().map_err(invalid_parameters)?;
    if !config.delta.is_finite() {
        return Err(invalid_parameters("delta must be finite".into()));
    }
    if x.len() != y.len() {
        return Err(invalid_parameters(format!(
            "x/y length mismatch: {} vs {}",
            x.len(),
            y.len()
        )));
    }
    if x.len() < 2 {
        return Err(invalid_parameters(
            "at least two paired samples are required".into(),
        ));
    }

    let invalid_x: Vec<_> = x
        .iter()
        .enumerate()
        .filter_map(|(index, value)| (!value.is_finite()).then_some(index))
        .collect();
    let invalid_y: Vec<_> = y
        .iter()
        .enumerate()
        .filter_map(|(index, value)| (!value.is_finite()).then_some(index))
        .collect();
    if !invalid_x.is_empty() || !invalid_y.is_empty() {
        let mut indices = invalid_x;
        indices.extend(invalid_y);
        indices.sort_unstable();
        indices.dedup();
        return Err(TGADomainError::InvalidColumnValues {
            column: column.into(),
            details: "LOWESS requires finite x and y values".into(),
            indices,
        }
        .with_context(operation, None, Some(column)));
    }
    Ok(())
}

fn panic_payload_message(payload: &(dyn Any + Send)) -> String {
    if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else if let Some(message) = payload.downcast_ref::<&'static str>() {
        (*message).to_string()
    } else {
        "backend panicked with a non-string payload".into()
    }
}

fn resolved_output_names(
    prefix: &str,
    input_cols: &[&str],
    output_cols: &[Option<&str>],
) -> Result<Vec<String>, TGADomainError> {
    if output_cols.is_empty() {
        return Ok(input_cols
            .iter()
            .map(|col| default_output_name(prefix, col, None))
            .collect());
    }

    if output_cols.len() != input_cols.len() {
        return Err(TGADomainError::InvalidOperation(format!(
            "Expected {} output column names, got {}",
            input_cols.len(),
            output_cols.len()
        )));
    }

    Ok(input_cols
        .iter()
        .zip(output_cols.iter())
        .map(|(&input, &output)| default_output_name(prefix, input, output))
        .collect())
}

fn maybe_filter_quality(
    raw: &[f64],
    filtered: &[f64],
    raw_name: &str,
    filtered_name: &str,
) -> Option<crate::Kinetics::experimental_kinetics::math::ndarray_statistics::FilterQualityReport> {
    filter_quality_from_arrays(
        &Array1::from_vec(raw.to_vec()),
        &Array1::from_vec(filtered.to_vec()),
        raw_name,
        filtered_name,
    )
    .ok()
}

fn log_filter_quality(
    method: &str,
    raw_name: &str,
    filtered_name: &str,
    quality: Option<
        &crate::Kinetics::experimental_kinetics::math::ndarray_statistics::FilterQualityReport,
    >,
) {
    match quality {
        Some(report) => {
            info!(
                "{} quality for '{}' -> '{}': {}",
                method, raw_name, filtered_name, report
            );
        }
        None => {
            warn!(
                "{} quality for '{}' -> '{}' could not be computed",
                method, raw_name, filtered_name
            );
        }
    }
}

impl TGADataset {
    pub fn smooth_columns(
        mut self,
        cols: &[&str],
        strategy: SmoothStrategy,
    ) -> Result<Self, TGADomainError> {
        for &col in cols {
            self = match &strategy {
                SmoothStrategy::RollingMean { window } => {
                    self.rolling_mean_as(col, *window, None)?
                }

                SmoothStrategy::Hampel {
                    window,
                    n_sigma,
                    strategy,
                } => self.hampel_filter_as(col, *window, *n_sigma, strategy.clone(), None)?,

                SmoothStrategy::SavitzkyGolay {
                    window,
                    poly_order,
                    deriv,
                    delta,
                    mode,
                } => match mode {
                    SGMode::FirPadding => {
                        self.sg_filter_column_as(col, *window, *poly_order, *deriv, *delta, None)?
                    }
                },

                SmoothStrategy::Lowess { frac } => {
                    self.lowess_filter_column_as(col, *frac, None)?
                }

                SmoothStrategy::Spline { .. } => {
                    return Err(TGADomainError::NotImplemented);
                }
            };
        }
        Ok(self)
    }
    /// Apply a trailing fixed-window mean without changing frame height.
    /// The first `window - 1` outputs and windows containing nulls are null.
    /// Input nulls are accepted; unit and nature are inherited by the output.
    pub fn rolling_mean_as(
        mut self,
        col_name: &str,
        window: usize,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        if window == 0 {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "rolling mean".into(),
                details: "window size must be greater than zero".into(),
            });
        }
        let out_col = default_output_name("rolling_", col_name, out_col);
        if out_col != col_name {
            self.ensure_column_name_is_available(&out_col)?;
        }
        let input_df = self.frame.clone().collect()?;
        if window > input_df.height() {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "rolling mean".into(),
                details: format!(
                    "window size {window} exceeds frame height {}",
                    input_df.height()
                ),
            });
        }
        let raw_values = extract_optional_f64_column(&input_df, col_name)?;

        let new_frame = self.frame.clone().with_column(
            col(col_name)
                .rolling_mean(RollingOptionsFixedWindow {
                    window_size: window,
                    min_periods: window,
                    ..Default::default()
                })
                .alias(&out_col),
        );
        let output_df = new_frame.clone().collect()?;
        let output_series = output_df.column(&out_col)?.f64()?;
        let mut raw_valid = Vec::new();
        let mut filtered_valid = Vec::new();
        for (index, raw) in raw_values.iter().enumerate() {
            if let (Some(raw), Some(filtered)) = (
                raw.and_then(|value| value.is_finite().then_some(value)),
                output_series.get(index),
            ) {
                if filtered.is_finite() {
                    raw_valid.push(raw);
                    filtered_valid.push(filtered);
                }
            }
        }

        self.push_undo_snapshot();

        info!(
            "new column with rolling averaged data {:?} has been created",
            out_col
        );
        self.schema
            .update_schema_for_numeric_transform(col_name, &out_col);
        let quality = maybe_filter_quality(&raw_valid, &filtered_valid, col_name, &out_col);
        log_filter_quality("rolling_mean", col_name, &out_col, quality.as_ref());
        self.record_column_transform(
            col_name,
            &out_col,
            ColumnTransformKind::RollingMean { window },
            true,
            quality,
        );
        self.log_operation(
            "rolling_mean",
            AffectedColumns::Specific(vec![out_col.clone()]),
            None,
            format!(
                "Rolling mean on '{}' -> '{}' (window={})",
                col_name, out_col, window
            ),
            true,
        );
        self.frame = output_df.lazy();

        Ok(self)
    }
    pub fn rolling_mean(self, col_name: &str, window: usize) -> Result<Self, TGADomainError> {
        self.rolling_mean_as(col_name, window, Some(col_name))
    }
    /// Apply a centered Hampel filter to a null-free column.
    /// Edge samples are copied unchanged; use `hampel_filter_null_safe` for
    /// nullable input.
    pub fn hampel_filter_as(
        mut self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: HampelStrategy,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let out_col = default_output_name("hampel_", col, out_col);
        validate_hampel_parameters(window, n_sigma)?;
        validate_filter_output(&self, col, &out_col)?;
        let mut df = self.frame.clone().collect()?;
        let s = df.column(col)?.f64()?;
        if s.null_count() > 0 {
            return Err(TGADomainError::InvalidOperation(format!(
                "Hampel filter '{}' requires a null-free input; use hampel_filter_null_safe for nullable data",
                col
            )));
        }

        let values: Vec<f64> = s.into_no_null_iter().collect();
        let invalid_indices: Vec<usize> = values
            .iter()
            .enumerate()
            .filter_map(|(index, value)| (!value.is_finite()).then_some(index))
            .collect();
        if !invalid_indices.is_empty() {
            return Err(TGADomainError::InvalidColumnValues {
                column: col.into(),
                details: "Hampel requires finite values".into(),
                indices: invalid_indices,
            });
        }
        if window > values.len() {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "Hampel".into(),
                details: format!("window {window} exceeds column length {}", values.len()),
            });
        }
        self.push_undo_snapshot();
        let is_drop = matches!(&strategy, HampelStrategy::Drop);
        let n = values.len();
        let k = window / 2;

        let mut mask = vec![true; n];
        let mut out = values.clone();

        for i in k..(n - k) {
            let w = &values[(i - k)..=(i + k)];
            let med = hampel_median(w.to_vec());
            let mad = median_abs_deviation(w, med);
            if !mad.is_finite() {
                return Err(TGADomainError::InvalidColumnValues {
                    column: col.into(),
                    details: "Hampel deviation overflowed".into(),
                    indices: vec![i],
                });
            }
            let thresh = n_sigma * 1.4826 * mad;

            if (values[i] - med).abs() > thresh {
                match strategy {
                    HampelStrategy::ReplaceWithMedian => out[i] = med,
                    HampelStrategy::ReplaceWithNaN => out[i] = f64::NAN,
                    HampelStrategy::Drop => mask[i] = false,
                }
            }
        }

        let mut new_df = df
            .with_column(Series::new(out_col.clone().into(), out).into())?
            .clone();
        if is_drop {
            new_df = new_df.filter(&BooleanChunked::from_slice("mask".into(), &mask))?;
        }
        let filtered_values = if is_drop {
            None
        } else {
            Some(extract_f64_column(&new_df, &out_col)?)
        };
        info!(
            "new column with Hampel filtered data named {:?} has been created",
            out_col
        );
        let list_of_cols = &new_df.get_column_names();
        info!("now list of columns is {:?}", list_of_cols);
        self.schema
            .update_schema_for_numeric_transform(col, &out_col);
        let quality = filtered_values
            .as_deref()
            .and_then(|filtered| maybe_filter_quality(&values, filtered, col, &out_col));
        log_filter_quality("hampel", col, &out_col, quality.as_ref());
        self.record_column_transform(
            col,
            &out_col,
            ColumnTransformKind::Hampel {
                window,
                n_sigma,
                strategy: format!("{:?}", &strategy),
            },
            true,
            quality,
        );
        self.log_operation(
            "hampel_filter",
            if is_drop {
                AffectedColumns::All
            } else {
                AffectedColumns::Specific(vec![out_col.clone()])
            },
            None,
            format!(
                "Hampel filter on '{}' -> '{}' (window={}, n_sigma={}, strategy={:?})",
                col, out_col, window, n_sigma, strategy
            ),
            true,
        );
        self.frame = new_df.clone().lazy();
        Ok(self)
    }
    pub fn hampel_filter(
        self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: HampelStrategy,
    ) -> Result<Self, TGADomainError> {
        self.hampel_filter_as(col, window, n_sigma, strategy, Some(col))
    }

    pub fn hampel_filter_null_safe_as(
        mut self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: HampelStrategy,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let out_col = default_output_name("hampel_", col, out_col);
        validate_hampel_parameters(window, n_sigma)?;
        validate_filter_output(&self, col, &out_col)?;
        let mut df = self.frame.clone().collect()?;
        let s = df.column(col)?.f64()?;

        let values: Vec<Option<f64>> = s.iter().collect();
        let invalid_indices: Vec<usize> = values
            .iter()
            .enumerate()
            .filter_map(|(index, value)| value.filter(|value| !value.is_finite()).map(|_| index))
            .collect();
        if !invalid_indices.is_empty() {
            return Err(TGADomainError::InvalidColumnValues {
                column: col.into(),
                details: "Hampel requires finite non-null values".into(),
                indices: invalid_indices,
            });
        }
        if window > values.len() {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "Hampel".into(),
                details: format!("window {window} exceeds column length {}", values.len()),
            });
        }
        self.push_undo_snapshot();
        let is_drop = matches!(&strategy, HampelStrategy::Drop);
        let n = values.len();
        let k = window / 2;

        let mut mask = vec![true; n];
        let mut out = values.clone();

        for i in k..(n - k) {
            let Some(xi) = values[i] else {
                continue;
            };

            let window_vals: Vec<f64> = values[(i - k)..=(i + k)]
                .iter()
                .filter_map(|v| *v)
                .collect();

            if window_vals.len() < 3 {
                continue;
            }

            let med = hampel_median(window_vals.clone());
            let mad = median_abs_deviation(&window_vals, med);

            if !mad.is_finite() {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Hampel deviation overflowed in column '{}' near row {}",
                    col, i
                )));
            }

            if mad == 0.0 {
                continue;
            }

            let thresh = n_sigma * 1.4826 * mad;

            if (xi - med).abs() > thresh {
                match strategy {
                    HampelStrategy::ReplaceWithMedian => out[i] = Some(med),
                    HampelStrategy::ReplaceWithNaN => out[i] = None,
                    HampelStrategy::Drop => mask[i] = false,
                }
            }
        }

        let new_df = df
            .with_column(Series::new(out_col.clone().into(), out).into())?
            .clone();
        let new_df = if is_drop {
            let mask = BooleanChunked::from_slice("mask".into(), &mask);
            new_df.filter(&mask)?
        } else {
            new_df
        };
        let filtered_values = if is_drop {
            None
        } else {
            Some(extract_optional_f64_column(&new_df, &out_col)?)
        };
        info!(
            "new column with Hampel filtered data named {:?} has been created",
            out_col
        );
        let list_of_cols = &new_df.get_column_names();
        info!("now list of columns is {:?}", list_of_cols);
        self.schema
            .update_schema_for_numeric_transform(col, &out_col);
        let quality = filtered_values.as_ref().and_then(|filtered| {
            let mut raw_valid = Vec::new();
            let mut filtered_valid = Vec::new();
            for (raw, output) in values.iter().zip(filtered) {
                if let (Some(raw), Some(output)) = (raw, output) {
                    if raw.is_finite() && output.is_finite() {
                        raw_valid.push(*raw);
                        filtered_valid.push(*output);
                    }
                }
            }
            maybe_filter_quality(&raw_valid, &filtered_valid, col, &out_col)
        });
        log_filter_quality("hampel", col, &out_col, quality.as_ref());
        self.record_column_transform(
            col,
            &out_col,
            ColumnTransformKind::Hampel {
                window,
                n_sigma,
                strategy: format!("{:?}", &strategy),
            },
            true,
            quality,
        );
        self.log_operation(
            "hampel_filter_null_safe",
            if is_drop {
                AffectedColumns::All
            } else {
                AffectedColumns::Specific(vec![out_col.clone()])
            },
            None,
            format!(
                "Null-safe Hampel filter on '{}' -> '{}' (window={}, n_sigma={}, strategy={:?})",
                col, out_col, window, n_sigma, strategy
            ),
            true,
        );
        self.frame = new_df.clone().lazy();
        Ok(self)
    }
    pub fn hampel_filter_null_safe(
        self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: HampelStrategy,
    ) -> Result<Self, TGADomainError> {
        self.hampel_filter_null_safe_as(col, window, n_sigma, strategy, Some(col))
    }
    /*
    pub fn ewm_mean(
        mut self,
        col_name: &str,
        alpha: f64,
    ) -> Self {
        self.frame = self.frame.with_column(
            col(col_name)
                .ewm_mean(EWMOptions {
                    alpha,
                    adjust: true,
                    min_periods: 1,
                    ..Default::default()
                })
                .alias(col_name)
        );
        self
    }
    */
    /// Apply a Savitzky-Golay filter using RustedSciThe's nearest-value edge padding.
    /// The result keeps the input length and has no nulls; null and non-finite input
    /// values are rejected before calling the backend.
    pub fn sg_filter_column_as(
        mut self,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let out_col = default_output_name("sg_", col, out_col);
        if window < 3
            || window % 2 == 0
            || poly_order > window.saturating_sub(2)
            || deriv > poly_order
            || !delta.is_finite()
            || delta <= 0.0
        {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "Savitzky-Golay".into(),
                details: "window must be odd and >= poly_order + 2; deriv <= poly_order; delta must be finite and > 0".into(),
            }
            .with_context("sg_filter_column", None, Some(col)));
        }
        validate_filter_output(&self, col, &out_col)?;
        let mut df = self.frame.clone().collect()?;

        let values = extract_f64_column(&df, col)
            .map_err(|error| error.with_context("sg_filter_column", None, Some(col)))?;
        let invalid_indices: Vec<usize> = values
            .iter()
            .enumerate()
            .filter_map(|(index, value)| (!value.is_finite()).then_some(index))
            .collect();
        if !invalid_indices.is_empty() {
            return Err(TGADomainError::InvalidColumnValues {
                column: col.into(),
                details: "Savitzky-Golay requires finite values".into(),
                indices: invalid_indices,
            }
            .with_context("sg_filter_column", None, Some(col)));
        }
        if window > values.len() {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "Savitzky-Golay".into(),
                details: format!("window {window} exceeds column length {}", values.len()),
            }
            .with_context("sg_filter_column", None, Some(col)));
        }
        let raw_values = values.clone();

        // The backend exposes a panic-based legacy API. Contain that failure at
        // this boundary and return it through the normal GUI-visible error path.
        let filtered = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            SG_filter_dyn(values.iter(), window, poly_order, Some(deriv), Some(delta))
        }))
        .map_err(|payload| {
            contextual_backend_error(
                "sg_filter_column",
                col,
                format!(
                    "Savitzky-Golay failed (window={window}, poly_order={poly_order}, deriv={deriv}, delta={delta}): {}",
                    panic_payload_message(payload.as_ref())
                ),
            )
        })?;
        validate_backend_output("sg_filter_column", col, &filtered, raw_values.len())?;
        let values = Column::new(out_col.clone().into(), filtered);
        let new_df = df
            .with_column(values)
            .map_err(TGADomainError::PolarsError)?;
        let filtered_values = extract_f64_column(&new_df, &out_col)?;

        self.push_undo_snapshot();
        info!(
            "new column with SG filtered data named {:?} has been created",
            out_col
        );
        let list_of_cols = &new_df.get_column_names();
        info!("now list of columns is {:?}", list_of_cols);
        self.frame = new_df.clone().lazy();
        self.schema
            .update_schema_for_numeric_transform(col, &out_col);
        let quality = maybe_filter_quality(&raw_values, &filtered_values, col, &out_col);
        log_filter_quality("savitzky_golay", col, &out_col, quality.as_ref());
        self.record_column_transform(
            col,
            &out_col,
            ColumnTransformKind::SavitzkyGolay {
                window,
                poly_order,
                deriv,
                delta,
                mode: "FirPadding".to_string(),
            },
            true,
            quality,
        );
        self.log_operation(
            "sg_filter_column",
            AffectedColumns::Specific(vec![out_col.clone()]),
            None,
            format!(
                "Savitzky-Golay filter on '{}' -> '{}' (window={}, order={}, deriv={})",
                col, out_col, window, poly_order, deriv
            ),
            true,
        );
        Ok(self)
    }
    pub fn sg_filter_column(
        self,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
    ) -> Result<Self, TGADomainError> {
        let start = Instant::now();
        let S = self.sg_filter_column_as(col, window, poly_order, deriv, delta, Some(col))?;
        info!(
            "SG filter completed in {:?} ms",
            start.elapsed().as_millis()
        );
        Ok(S)
    }
    //===============================================================================
    fn lowess_reference_x(&self, df: &DataFrame) -> Result<Vec<f64>, TGADomainError> {
        if let Some(time_col) = self.schema.time.as_deref() {
            if df
                .get_column_names()
                .iter()
                .any(|name| name.as_str() == time_col)
            {
                return extract_f64_column(df, time_col);
            }
        }
        Ok((0..df.height()).map(|i| i as f64).collect())
    }

    /// Apply LOWESS against bound time (or row indices when time is unbound).
    /// It preserves frame height and rejects null/non-finite inputs.
    pub fn lowess_filter_column(self, col: &str, frac: f64) -> Result<Self, TGADomainError> {
        self.lowess_filter_column_as(col, frac, Some(col))
    }

    pub fn lowess_filter_column_as(
        mut self,
        col: &str,
        frac: f64,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let out_col = default_output_name("lowess_", col, out_col);
        let mut df = self.frame.clone().collect()?;
        let x = self.lowess_reference_x(&df)?;
        let y = extract_f64_column(&df, col)?;
        let raw_values = y.clone();
        let cfg = LowessConfig {
            fraction: frac,
            ..LowessConfig::default()
        };
        let x_column = self.schema.time.as_deref().unwrap_or("<row index>");
        validate_lowess_inputs("lowess_filter_column", col, x_column, &x, &y, &cfg)?;
        validate_filter_output(&self, col, &out_col)?;

        let filtered = lowess_smooth_values(&x, &y, &cfg).map_err(|e| {
            contextual_backend_error(
                "lowess_filter_column",
                col,
                format!("LOWESS fit failed (config={cfg:?}): {e:?}"),
            )
        })?;
        validate_backend_output("lowess_filter_column", col, &filtered, y.len())?;
        self.push_undo_snapshot();
        df.with_column(Column::new(out_col.clone().into(), filtered))?;
        let filtered_values = extract_f64_column(&df, &out_col)?;
        self.schema
            .update_schema_for_numeric_transform(col, &out_col);
        let quality = maybe_filter_quality(&raw_values, &filtered_values, col, &out_col);
        log_filter_quality("lowess", col, &out_col, quality.as_ref());
        self.frame = df.lazy();
        self.record_column_transform(
            col,
            &out_col,
            ColumnTransformKind::Lowess { frac },
            true,
            quality,
        );
        self.log_operation(
            "lowess_filter",
            AffectedColumns::Specific(vec![out_col.clone()]),
            None,
            format!(
                "LOWESS filter on '{}' -> '{}' (fraction={})",
                col, out_col, frac
            ),
            true,
        );
        Ok(self)
    }

    pub fn lowess_smooth_column(
        self,
        time_col: &str,
        column: &str,
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        self.lowess_smooth_column_as(time_col, column, Some(column), config)
    }

    pub fn lowess_smooth_column_as(
        self,
        time_col: &str,
        column: &str,
        out_column: Option<&str>,
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        self.lowess_smooth_columns_as(time_col, &[column], &[out_column], config)
    }

    pub fn lowess_smooth_columns(
        self,
        time_col: &str,
        columns: &[&str],
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        let out_cols: Vec<Option<&str>> = columns.iter().map(|&col| Some(col)).collect();
        self.lowess_smooth_columns_as(time_col, columns, &out_cols, config)
    }

    /// Smooth multiple columns against an explicit X column in one operation.
    /// Frame height is preserved; null and non-finite inputs are rejected.
    pub fn lowess_smooth_columns_as(
        mut self,
        time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        let mut df = self.frame.clone().collect()?;
        if columns.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "LOWESS requires at least one data column".into(),
            ));
        }
        let out_names = resolved_output_names("lowess_", columns, out_columns)?;
        let x = extract_f64_column(&df, time_col)?;
        validate_output_mapping(&self, columns, &out_names)?;
        let mut results = Vec::with_capacity(columns.len());
        for (&col, out_name) in columns.iter().zip(out_names.iter()) {
            let y = extract_f64_column(&df, col)?;
            validate_lowess_inputs("lowess_smooth_columns", col, time_col, &x, &y, &config)?;
            let smoothed = lowess_smooth_values(&x, &y, &config).map_err(|e| {
                contextual_backend_error(
                    "lowess_smooth_columns",
                    col,
                    format!("LOWESS fit failed (config={config:?}): {e:?}"),
                )
            })?;
            validate_backend_output("lowess_smooth_columns", col, &smoothed, y.len())?;
            results.push((col, out_name, y, smoothed));
        }

        self.push_undo_snapshot();
        for (col, out_name, raw_values, smoothed) in results {
            let quality = maybe_filter_quality(&raw_values, &smoothed, col, out_name);
            log_filter_quality("lowess", col, out_name, quality.as_ref());

            df.with_column(Column::new(out_name.clone().into(), smoothed))?;
            self.schema
                .update_schema_for_numeric_transform(col, out_name);
            self.record_column_transform(
                col,
                out_name,
                ColumnTransformKind::Lowess {
                    frac: config.fraction,
                },
                true,
                quality,
            );
        }

        self.frame = df.lazy();
        self.log_operation(
            "lowess_smooth_columns",
            AffectedColumns::Specific(out_names.clone()),
            None,
            format!(
                "LOWESS smoothing of columns {:?} using '{}' (fraction={})",
                columns, time_col, config.fraction
            ),
            true,
        );
        Ok(self)
    }

    //===============================================================================

    // SPLINES

    /// Resample selected values on a new grid while preserving the frame's row count.
    ///
    /// Existing columns remain in the frame. Resampled columns and the new grid are
    /// padded with trailing nulls when `n_points` is smaller than the source height.
    /// With explicit output names, source columns remain unchanged; omitted names
    /// overwrite the selected columns. Output metadata inherits unit and nature,
    /// records a numeric-derived origin, and the active binding follows the output.
    pub fn spline_resample_columns(
        self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let out_cols: Vec<Option<&str>> = columns.iter().map(|&col| Some(col)).collect();
        self.spline_resample_columns_as(time_col, new_time_col, columns, &out_cols, n_points, kind)
    }

    pub fn spline_resample_columns_as(
        mut self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        validate_spline_kind(kind)?;
        let mut df = self.frame.clone().collect()?;
        let n = df.height();
        if n_points < 2 || n_points > n {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "spline resampling".into(),
                details: format!(
                    "n_points must be in 2..=frame height (requested {}, rows {})",
                    n_points, n
                ),
            }
            .with_context("spline_resample", None, Some(time_col)));
        }
        if columns.is_empty() || columns.contains(&time_col) {
            return Err(TGADomainError::InvalidOperation(
                "Spline resampling requires at least one non-time data column".into(),
            ));
        }
        let out_names = resolved_output_names("spline_", columns, out_columns)?;

        // --- old time ---
        let t_old = extract_f64_column(&df, time_col)?;
        if t_old.iter().any(|value| !value.is_finite())
            || t_old.windows(2).any(|pair| pair[1] <= pair[0])
        {
            return Err(TGADomainError::InvalidOperation(
                "Spline resampling requires finite, strictly increasing time values".into(),
            ));
        }
        let mut all_inputs = columns.to_vec();
        all_inputs.push(time_col);
        let mut all_outputs = out_names.clone();
        all_outputs.push(new_time_col.to_string());
        validate_output_mapping(&self, &all_inputs, &all_outputs)?;
        self.push_undo_snapshot();

        // --- resample requested columns ---
        for (&col, out_name) in columns.iter().zip(out_names.iter()) {
            let y_old = extract_f64_column(&df, col)?;
            if y_old.iter().any(|value| !value.is_finite()) {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Spline resampling requires finite values in column '{}'",
                    col
                )));
            }

            let (_, y_new) = spline_resample(&t_old, &y_old, n_points, kind).map_err(|error| {
                contextual_backend_error(
                    "spline_resample",
                    col,
                    format!("spline evaluation failed (points={n_points}, kind={kind:?}): {error}"),
                )
            })?;
            validate_backend_output("spline_resample", col, &y_new, n_points)?;
            let y_series = Column::new(out_name.clone().into(), y_new);
            let n_new = y_series.len();
            let y_series = y_series.extend_constant(AnyValue::Null, n - n_new)?;
            df = df.with_column(y_series)?.clone();
            self.schema
                .update_schema_for_numeric_transform(col, &out_name);
            self.record_column_transform(
                col,
                out_name,
                ColumnTransformKind::SplineResample {
                    time_col: time_col.to_string(),
                    new_time_col: new_time_col.to_string(),
                    n_points,
                    kind: format!("{:?}", kind),
                },
                true,
                None,
            );
        }
        info!("updating ");
        // --- new grid ---
        let t_new = uniform_grid_from(&t_old, n_points, kind).map_err(|error| {
            contextual_backend_error(
                "spline_resample",
                time_col,
                format!("grid generation failed (points={n_points}, kind={kind:?}): {error}"),
            )
        })?;
        validate_backend_output("spline_resample", time_col, &t_new, n_points)?;

        // --- add new time column ---
        let time_column = Column::new(new_time_col.into(), t_new);
        let n_new = time_column.len();
        let time_column = time_column.extend_constant(AnyValue::Null, n - n_new)?;
        df = df.clone().with_column(time_column)?.clone();
        self.schema
            .update_schema_for_numeric_transform(time_col, new_time_col);
        self.record_column_transform(
            time_col,
            new_time_col,
            ColumnTransformKind::SplineResample {
                time_col: time_col.to_string(),
                new_time_col: new_time_col.to_string(),
                n_points,
                kind: format!("{:?}", kind),
            },
            true,
            None,
        );
        let list_of_cols = &df.get_column_names();
        info!("now list of columns is {:?}", list_of_cols);
        self.frame = df.lazy();
        self.log_operation(
            "spline_resample",
            AffectedColumns::All,
            None,
            format!(
                "Spline resampling of {:?} and time '{}' -> '{}' (points={}, kind={:?})",
                columns, time_col, new_time_col, n_points, kind
            ),
            true,
        );
        info!(
            "new columns with splined data named {:?} have been created",
            out_columns
        );

        Ok(self)
    }

    pub fn lsq_spline_resample_column(
        self,
        time_col: &str,
        new_time_col: &str,
        column: &str,
        n_points: usize,
    ) -> Result<Self, TGADomainError> {
        self.lsq_spline_resample_column_as(
            time_col,
            new_time_col,
            column,
            Some(column),
            n_points,
            3,
            24,
            SolverKind::Banded,
        )
    }

    pub fn lsq_spline_resample_column_as(
        self,
        time_col: &str,
        new_time_col: &str,
        column: &str,
        out_column: Option<&str>,
        n_points: usize,

        degree: usize,
        n_internal_knots: usize,
        solver: SolverKind,
    ) -> Result<Self, TGADomainError> {
        self.lsq_spline_resample_columns_as(
            time_col,
            new_time_col,
            &[column],
            &[out_column],
            n_points,
            degree,
            n_internal_knots,
            solver,
        )
    }

    pub fn lsq_spline_resample_columns(
        self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        n_points: usize,
    ) -> Result<Self, TGADomainError> {
        let out_cols: Vec<Option<&str>> = columns.iter().map(|&col| Some(col)).collect();
        self.lsq_spline_resample_columns_as(
            time_col,
            new_time_col,
            columns,
            &out_cols,
            n_points,
            3,
            24,
            SolverKind::Banded,
        )
    }

    pub fn lsq_spline_resample_columns_as(
        mut self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        n_points: usize,
        degree: usize,
        n_internal_knots: usize,
        solver: SolverKind,
    ) -> Result<Self, TGADomainError> {
        let start = Instant::now();
        if n_points < 2 {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "LSQ spline resampling".into(),
                details: "n_points must be >= 2".into(),
            }
            .with_context("lsq_spline_resample", None, Some(time_col)));
        }
        if degree == 0 {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "LSQ spline resampling".into(),
                details: "degree must be >= 1".into(),
            }
            .with_context("lsq_spline_resample", None, Some(time_col)));
        }

        let mut df = self.frame.clone().collect()?;
        let n = df.height();
        if n_points > n {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "LSQ spline resampling".into(),
                details: format!(
                    "cannot preserve row alignment when n_points ({}) exceeds frame height ({})",
                    n_points, n
                ),
            }
            .with_context("lsq_spline_resample", None, Some(time_col)));
        }
        if columns.is_empty() || columns.contains(&time_col) {
            return Err(TGADomainError::InvalidOperation(
                "LSQ spline resampling requires at least one non-time data column".into(),
            ));
        }
        let out_names = resolved_output_names("lsq_spline_", columns, out_columns)?;
        let t_old = extract_f64_column(&df, time_col)?;
        if t_old.len() < 2 {
            return Err(TGADomainError::InvalidOperation(
                "LSQ spline resampling requires at least 2 points".into(),
            ));
        }
        if t_old.iter().any(|value| !value.is_finite()) {
            return Err(TGADomainError::InvalidOperation(
                "LSQ spline resampling requires finite time values".into(),
            ));
        }
        if degree >= t_old.len() {
            return Err(TGADomainError::InvalidOperation(format!(
                "LSQ spline degree {} must be smaller than the {} input points",
                degree,
                t_old.len()
            )));
        }
        for i in 1..t_old.len() {
            if t_old[i] <= t_old[i - 1] {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Time column must be strictly monotonic for LSQ spline resampling:
                     t[{}] {} is smaller than t[{}] = {} ",
                    i,
                    t_old[i],
                    i - 1,
                    t_old[i - 1]
                )));
            }
        }

        let t_new = uniform_grid_from(&t_old, n_points, SplineKind::Linear).map_err(|error| {
            contextual_backend_error(
                "lsq_spline_resample",
                time_col,
                format!("grid generation failed (points={n_points}, kind=Linear): {error}"),
            )
        })?;
        validate_backend_output("lsq_spline_resample", time_col, &t_new, n_points)?;
        if t_new.windows(2).any(|pair| pair[1] <= pair[0]) {
            warn!(
                "LSQ spline generated a non-monotonic grid for '{}' -> '{}'",
                time_col, new_time_col
            );
            return Err(contextual_backend_error(
                "lsq_spline_resample",
                time_col,
                "generated a non-monotonic time grid",
            ));
        }
        let mut all_inputs = columns.to_vec();
        all_inputs.push(time_col);
        let mut all_outputs = out_names.clone();
        all_outputs.push(new_time_col.to_string());
        validate_output_mapping(&self, &all_inputs, &all_outputs)?;
        let max_internal_knots = t_old.len() - degree - 1;
        if n_internal_knots > max_internal_knots {
            return Err(TGADomainError::InvalidFilterParameters {
                filter: "LSQ spline resampling".into(),
                details: format!(
                    "requested {} internal knots, but at most {} are valid for {} points and degree {}",
                    n_internal_knots, max_internal_knots, t_old.len(), degree
                ),
            }
            .with_context("lsq_spline_resample", None, Some(time_col)));
        }
        self.push_undo_snapshot();
        let solver_label = match solver {
            SolverKind::DenseQR => "DenseQR".to_string(),
            SolverKind::Banded => "Banded".to_string(),
        };
        info!(
            "Starting LSQ spline resampling '{}' -> '{}' for columns {:?} | input_rows={} target_points={} degree={} internal_knots={} solver={} | generated monotonic grid points={} spanning [{:.6}, {:.6}]",
            time_col,
            new_time_col,
            columns,
            n,
            n_points,
            degree,
            n_internal_knots,
            solver_label,
            t_new.len(),
            t_new.first().copied().unwrap_or(f64::NAN),
            t_new.last().copied().unwrap_or(f64::NAN)
        );
        let t_old_arr = Array1::from_vec(t_old.clone());
        let t_new_arr = Array1::from_vec(t_new.clone());

        let used_internal_knots = n_internal_knots;
        let xmin = t_old[0];
        let xmax = t_old[t_old.len() - 1];
        let internal_knots = if used_internal_knots == 0 {
            Array1::<f64>::zeros(0)
        } else {
            Array1::linspace(xmin, xmax, used_internal_knots + 2)
                .slice(s![1..used_internal_knots + 1])
                .to_owned()
        };

        for (&col, out_name) in columns.iter().zip(out_names.iter()) {
            let y_old = extract_f64_column(&df, col)?;
            if y_old.iter().any(|value| !value.is_finite()) {
                return Err(TGADomainError::InvalidOperation(format!(
                    "LSQ spline resampling requires finite values in column '{}'",
                    col
                )));
            }
            let y_old_arr = Array1::from_vec(y_old);

            let spline =
                make_lsq_univariate_spline(&t_old_arr, &y_old_arr, &internal_knots, degree, solver)
                    .map_err(|error| {
                        contextual_backend_error(
                            "lsq_spline_resample",
                            col,
                            format!(
                                "spline fit failed (degree={degree}, internal_knots={n_internal_knots}, solver={solver_label}): {error:?}"
                            ),
                        )
                    })?;

            let y_new = spline.evaluate_batch_array(&t_new_arr).to_vec();
            validate_backend_output("lsq_spline_resample", col, &y_new, n_points)?;
            let y_series = Column::new(out_name.clone().into(), y_new);
            let n_new = y_series.len();
            let y_series = y_series.extend_constant(AnyValue::Null, n - n_new)?;
            df = df.with_column(y_series)?.clone();
            self.schema
                .update_schema_for_numeric_transform(col, out_name);
            self.record_column_transform(
                col,
                out_name,
                ColumnTransformKind::LsqSplineResample {
                    time_col: time_col.to_string(),
                    new_time_col: new_time_col.to_string(),
                    n_points,
                    degree,
                    n_internal_knots,
                    solver: solver_label.clone(),
                },
                true,
                None,
            );
        }

        let time_column = Column::new(new_time_col.into(), t_new);
        let n_new = time_column.len();
        let time_column = time_column.extend_constant(AnyValue::Null, n - n_new)?;
        df = df.with_column(time_column)?.clone();
        warn!(
            "LSQ spline resampling does not compute a direct filter-quality report because the output grid differs from the input and the frame is padded to preserve alignment"
        );
        self.schema
            .update_schema_for_numeric_transform(time_col, new_time_col);
        self.record_column_transform(
            time_col,
            new_time_col,
            ColumnTransformKind::LsqSplineResample {
                time_col: time_col.to_string(),
                new_time_col: new_time_col.to_string(),
                n_points,
                degree,
                n_internal_knots,
                solver: solver_label.clone(),
            },
            true,
            None,
        );
        self.frame = df.lazy();
        self.log_operation(
            "lsq_spline_resample",
            AffectedColumns::All,
            None,
            format!(
                "LSQ spline resampling of {:?} and time '{}' -> '{}' (points={}, degree={}, internal_knots={}, solver={})",
                columns, time_col, new_time_col, n_points, degree, n_internal_knots, solver_label
            ),
            true,
        );
        info!(
            "LSQ spline resampling completed in {:?} ms",
            start.elapsed().as_millis()
        );
        Ok(self)
    }
}
