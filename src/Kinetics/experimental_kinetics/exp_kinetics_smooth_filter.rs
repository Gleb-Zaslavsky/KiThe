use crate::Kinetics::experimental_kinetics::LSQSplines::{SolverKind, make_lsq_univariate_spline};
use crate::Kinetics::experimental_kinetics::lowess_wrapper::{LowessConfig, lowess_smooth_values};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{TGADataset, TGADomainError};
use crate::Kinetics::experimental_kinetics::savgol2::SG_filter_dyn;
use crate::Kinetics::experimental_kinetics::splines::{
    SplineKind, spline_resample, uniform_grid_from,
};
use log::info;
use ndarray::{Array1, s};
use polars::prelude::*;
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

    let s = s
        .cast(&DataType::Float64)
        .map_err(|_| TGADomainError::InvalidColumnType(col.into()))?;
    let s = s.f64().unwrap();

    if s.null_count() > 0 {
        return Err(TGADomainError::InvalidOperation(format!(
            "Operation does not support nulls in column '{}'",
            col
        )));
    }

    Ok(s.into_no_null_iter().collect())
}
pub fn resample_all_columns(
    df: &DataFrame,
    time_col: &str,
    x_new: &[f64],
    kind: SplineKind,
) -> Result<DataFrame, TGADomainError> {
    let mut cols: Vec<Column> = Vec::new();

    for col in df.columns() {
        let name = col.name();

        let y_old = extract_f64_column(df, name)?;
        let x_old = extract_f64_column(df, time_col)?;

        let (_, y_new) = spline_resample(&x_old, &y_old, x_new.len(), kind)
            .map_err(|e| TGADomainError::InvalidOperation(e))?;

        cols.push(Series::new(name.clone(), y_new).into());
    }

    let height = cols.iter().map(|c| c.len()).max().unwrap_or(0);
    Ok(DataFrame::new(height, cols)?)
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
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 0 {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
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

impl TGADataset {
    pub fn smooth_columns(
        mut self,
        cols: &[&str],
        strategy: SmoothStrategy,
    ) -> Result<Self, TGADomainError> {
        for &col in cols {
            self = match &strategy {
                SmoothStrategy::RollingMean { window } => self.rolling_mean(col, *window),

                SmoothStrategy::Hampel {
                    window,
                    n_sigma,
                    strategy,
                } => self.hampel_filter(col, *window, *n_sigma, strategy.clone())?,

                SmoothStrategy::SavitzkyGolay {
                    window,
                    poly_order,
                    deriv,
                    delta,
                    mode,
                } => match mode {
                    SGMode::FirPadding => {
                        self.sg_filter_column(col, *window, *poly_order, *deriv, *delta)?
                    }
                },

                SmoothStrategy::Lowess { frac } => self.lowess_filter_column(col, *frac)?,

                SmoothStrategy::Spline { .. } => {
                    return Err(TGADomainError::NotImplemented);
                }
            };
        }
        Ok(self)
    }
    /// Apply a simple moving average filter to smooth the TGA data.
    //=======================================================
    // attention: feature "rolling_window" must be enabled
    ///  длина таблицы не меняется
    ///  первые window-1 значений → null
    /// metadata (unit, origin) не меняются
    pub fn rolling_mean_as(mut self, col_name: &str, window: usize, out_col: Option<&str>) -> Self {
        let out_col = default_output_name("rolling_", col_name, out_col);
        self.frame = self.frame.with_column(
            col(col_name)
                .rolling_mean(RollingOptionsFixedWindow {
                    window_size: window,
                    min_periods: window,
                    ..Default::default()
                })
                .alias(&out_col),
        );

        info!(
            "new column with rolling averaged data {:?} has been created",
            out_col
        );
        self.schema.update_schema(col_name, &out_col);

        self
    }
    pub fn rolling_mean(self, col_name: &str, window: usize) -> Self {
        self.rolling_mean_as(col_name, window, Some(col_name))
    }
    /// Hampel применяется к колонке без null
    /// Рекомендуется вызывать после trim_edges
    pub fn hampel_filter_as(
        mut self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: HampelStrategy,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let out_col = default_output_name("hampel_", col, out_col);
        let mut df = self.frame.clone().collect()?;
        let s = df.column(col)?.f64()?;

        let values: Vec<f64> = s.into_no_null_iter().collect();
        let n = values.len();
        let k = window / 2;

        let mut mask = vec![true; n];
        let mut out = values.clone();

        for i in k..(n - k) {
            let w = &values[(i - k)..=(i + k)];
            let med = hampel_median(w.to_vec());
            let mad = median_abs_deviation(w, med);
            let thresh = n_sigma * 1.4826 * mad;

            if (values[i] - med).abs() > thresh {
                match strategy {
                    HampelStrategy::ReplaceWithMedian => out[i] = med,
                    HampelStrategy::ReplaceWithNaN => out[i] = f64::NAN,
                    HampelStrategy::Drop => mask[i] = false,
                }
            }
        }

        let new_df = df.with_column(Series::new(out_col.clone().into(), out).into())?;

        let new_df = if let HampelStrategy::Drop = strategy {
            &mut new_df.filter(&BooleanChunked::from_slice("mask".into(), &mask))?
        } else {
            new_df
        };
        info!(
            "new column with Hampel filtered data named {:?} has been created",
            out_col
        );
        let list_of_cols = &new_df.get_column_names();
        info!("now list of columns is {:?}", list_of_cols);
        self.schema.update_schema(col, &out_col);
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
        let mut df = self.frame.clone().collect()?;
        let s = df.column(col)?.f64()?;

        let values: Vec<Option<f64>> = s.into_iter().collect();
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

        let new_df = df.with_column(Series::new(out_col.clone().into(), out).into())?;

        let new_df = if let HampelStrategy::Drop = strategy {
            let mask = BooleanChunked::from_slice("mask".into(), &mask);
            &mut new_df.filter(&mask)?
        } else {
            new_df
        };
        info!(
            "new column with Hampel filtered data named {:?} has been created",
            out_col
        );
        let list_of_cols = &new_df.get_column_names();
        info!("now list of columns is {:?}", list_of_cols);
        self.schema.update_schema(col, &out_col);
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
    /// Savitzky–Golay smoothing (FIR-based).
    ///
    /// Properties:
    /// - output length == input length
    /// - does NOT introduce nulls
    /// - edge handling via nearest-value padding
    /// - suitable for smoothing and low-order derivatives
    ///
    /// Note:
    /// If edge effects are undesirable, call `trim_edges` explicitly.
    /// rolling_mean НЕ уменьшает число строк

    ///результат всегда той же длины

    /// если min_periods = window, то:

    /// первые window - 1 значений становятся NULL

    /// в конце — NULL не появляется
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
        let mut df = self.frame.clone().collect()?;

        let s = df
            .column(col)
            .map_err(TGADomainError::PolarsError)?
            .f64()
            .map_err(TGADomainError::PolarsError)?;

        if s.null_count() > 0 {
            return Err(TGADomainError::InvalidOperation(format!(
                "SG filter does not support nulls in column '{}'",
                col
            )));
        }

        let values: Vec<f64> = s.into_no_null_iter().collect();

        let filtered = SG_filter_dyn(values.iter(), window, poly_order, Some(deriv), Some(delta));
        let values = Column::new(out_col.clone().into(), filtered);
        let new_df = df
            .with_column(values)
            .map_err(TGADomainError::PolarsError)?;

        info!(
            "new column with SG filtered data named {:?} has been created",
            out_col
        );
        let list_of_cols = &new_df.get_column_names();
        info!("now list of columns is {:?}", list_of_cols);
        self.frame = new_df.clone().lazy();
        self.schema.update_schema(col, &out_col);
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

        if x.len() != y.len() {
            return Err(TGADomainError::InvalidOperation(format!(
                "LOWESS input length mismatch: x has {}, y has {} for column '{}'",
                x.len(),
                y.len(),
                col
            )));
        }

        let cfg = LowessConfig {
            fraction: frac,
            ..LowessConfig::default()
        };

        let filtered = lowess_smooth_values(&x, &y, &cfg).map_err(|e| {
            TGADomainError::InvalidOperation(format!("LOWESS fit failed for '{}': {:?}", col, e))
        })?;

        df.with_column(Column::new(out_col.clone().into(), filtered))?;
        self.schema.update_schema(col, &out_col);
        self.frame = df.lazy();
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

    pub fn lowess_smooth_columns_as(
        mut self,
        time_col: &str,
        columns: &[&str],
        out_columns: &[Option<&str>],
        config: LowessConfig,
    ) -> Result<Self, TGADomainError> {
        let mut df = self.frame.clone().collect()?;
        let out_names = resolved_output_names("lowess_", columns, out_columns)?;
        let x = extract_f64_column(&df, time_col)?;

        for (&col, out_name) in columns.iter().zip(out_names.iter()) {
            let y = extract_f64_column(&df, col)?;
            if x.len() != y.len() {
                return Err(TGADomainError::InvalidOperation(format!(
                    "LOWESS input length mismatch for column '{}': x has {}, y has {}",
                    col,
                    x.len(),
                    y.len()
                )));
            }

            let smoothed = lowess_smooth_values(&x, &y, &config).map_err(|e| {
                TGADomainError::InvalidOperation(format!(
                    "LOWESS fit failed for column '{}': {:?}",
                    col, e
                ))
            })?;

            df.with_column(Column::new(out_name.clone().into(), smoothed))?;
            self.schema.update_schema(col, out_name);
        }

        self.frame = df.lazy();
        Ok(self)
    }

    //===============================================================================

    // SPLINES

    /// Smooth a column using spline interpolation over time axis.
    ///
    /// # Preconditions
    /// - time column must be bound
    /// - or you can use another column as x argument (but it must be present)
    /// - target column must exist
    /// - no nulls in time or target column
    /// - time must be strictly monotonic
    ///
    /// # Guarantees
    /// - preserves row count
    /// - does NOT introduce nulls
    /// - overwrites the column in-place
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
        let mut df = self.frame.clone().collect()?;
        let n = df.height();
        let out_names = resolved_output_names("spline_", columns, out_columns)?;

        // --- old time ---
        let t_old = extract_f64_column(&df, time_col)?;

        // --- resample requested columns ---
        for (&col, out_name) in columns.iter().zip(out_names.iter()) {
            let y_old = extract_f64_column(&df, col)?;

            let (_, y_new) = spline_resample(&t_old, &y_old, n_points, kind)
                .map_err(|e| TGADomainError::InvalidOperation(e))?;
            let y_series = Column::new(out_name.clone().into(), y_new);
            let n_new = y_series.len();
            let y_series = y_series.extend_constant(AnyValue::Null, n - n_new)?;
            df = df.with_column(y_series)?.clone();
            self.schema.update_schema(col, &out_name);
        }
        info!("updating ");
        // --- new grid ---
        let t_new = uniform_grid_from(&t_old, n_points, kind)
            .map_err(|e| TGADomainError::InvalidOperation(e))?;

        // --- add new time column ---
        let time_column = Column::new(new_time_col.into(), t_new);
        let n_new = time_column.len();
        let time_column = time_column.extend_constant(AnyValue::Null, n - n_new)?;
        df = df.clone().with_column(time_column)?.clone();
        self.schema.update_schema(time_col, new_time_col);
        let list_of_cols = &df.get_column_names();
        info!("now list of columns is {:?}", list_of_cols);
        self.frame = df.lazy();
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
            return Err(TGADomainError::InvalidOperation(
                "n_points must be >= 2".into(),
            ));
        }
        if degree == 0 {
            return Err(TGADomainError::InvalidOperation(
                "degree must be >= 1".into(),
            ));
        }

        let mut df = self.frame.clone().collect()?;
        let n = df.height();
        let out_names = resolved_output_names("lsq_spline_", columns, out_columns)?;

        let t_old = extract_f64_column(&df, time_col)?;
        if t_old.len() < 2 {
            return Err(TGADomainError::InvalidOperation(
                "LSQ spline resampling requires at least 2 points".into(),
            ));
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

        let t_new = uniform_grid_from(&t_old, n_points, SplineKind::Linear)
            .map_err(|e| TGADomainError::InvalidOperation(e))?;
        let t_old_arr = Array1::from_vec(t_old.clone());
        let t_new_arr = Array1::from_vec(t_new.clone());

        let max_internal_knots = t_old.len().saturating_sub(degree + 1);
        let used_internal_knots = n_internal_knots.min(max_internal_knots);
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
            let y_old_arr = Array1::from_vec(y_old);

            let spline =
                make_lsq_univariate_spline(&t_old_arr, &y_old_arr, &internal_knots, degree, solver)
                    .map_err(|e| {
                        TGADomainError::InvalidOperation(format!(
                            "LSQ spline fit failed for column '{}': {:?}",
                            col, e
                        ))
                    })?;

            let y_new = spline.evaluate_batch_array(&t_new_arr).to_vec();
            let y_series = Column::new(out_name.clone().into(), y_new);
            let n_new = y_series.len();
            let y_series = y_series.extend_constant(AnyValue::Null, n - n_new)?;
            df = df.with_column(y_series)?.clone();
            self.schema.update_schema(col, out_name);
        }

        let time_column = Column::new(new_time_col.into(), t_new);
        let n_new = time_column.len();
        let time_column = time_column.extend_constant(AnyValue::Null, n - n_new)?;
        df = df.with_column(time_column)?.clone();
        self.schema.update_schema(time_col, new_time_col);
        self.frame = df.lazy();
        info!(
            "LSQ spline resampling completed in {:?} ms",
            start.elapsed().as_millis()
        );
        Ok(self)
    }
}
