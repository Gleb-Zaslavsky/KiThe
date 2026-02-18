use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{TGADataset, TGADomainError};
use crate::Kinetics::experimental_kinetics::savgol2::SG_filter_dyn;
use crate::Kinetics::experimental_kinetics::splines::{
    SplineKind, spline_resample, uniform_grid_from,
};
use polars::prelude::*;
use std::f64;
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

    for col in df.get_columns() {
        let name = col.name();

        let y_old = extract_f64_column(df, name)?;
        let x_old = extract_f64_column(df, time_col)?;

        let (_, y_new) = spline_resample(&x_old, &y_old, x_new.len(), kind)
            .map_err(|e| TGADomainError::InvalidOperation(e))?;

        cols.push(Series::new(name.clone(), y_new).into());
    }

    Ok(DataFrame::new(cols)?)
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

                SmoothStrategy::Lowess { .. } | SmoothStrategy::Spline { .. } => {
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
    pub fn rolling_mean(mut self, col_name: &str, window: usize) -> Self {
        self.frame = self.frame.with_column(
            col(col_name)
                .rolling_mean(RollingOptionsFixedWindow {
                    window_size: window,
                    min_periods: window,
                    ..Default::default()
                })
                .alias(col_name),
        );

        self
    }
    /// Hampel применяется к колонке без null
    /// Рекомендуется вызывать после trim_edges
    pub fn hampel_filter(
        mut self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: HampelStrategy,
    ) -> Result<Self, TGADomainError> {
        let df = self.frame.clone().collect()?;
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

        let mut new_df = df.clone();
        new_df.replace(col, Series::new(col.into(), out))?;

        if let HampelStrategy::Drop = strategy {
            new_df = new_df.filter(&BooleanChunked::from_slice("mask".into(), &mask))?;
        }

        self.frame = new_df.lazy();
        Ok(self)
    }

    pub fn hampel_filter_null_safe(
        mut self,
        col: &str,
        window: usize,
        n_sigma: f64,
        strategy: HampelStrategy,
    ) -> Result<Self, TGADomainError> {
        let df = self.frame.clone().collect()?;
        let s = df.column(col)?.f64()?;

        let values: Vec<Option<f64>> = s.into_iter().collect();
        let n = values.len();
        let k = window / 2;

        let mut mask = vec![true; n];
        let mut out = values.clone();

        for i in k..(n - k) {
            let Some(xi) = values[i] else {
                continue; // null не фильтруем
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

        let mut new_df = df.clone();
        let new_series = Series::new(col.into(), out);
        new_df.replace(col, new_series)?;

        if let HampelStrategy::Drop = strategy {
            let mask = BooleanChunked::from_slice("mask".into(), &mask);
            new_df = new_df.filter(&mask)?;
        }

        self.frame = new_df.lazy();
        Ok(self)
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
    pub fn sg_filter_column(
        mut self,
        col: &str,
        window: usize,
        poly_order: usize,
        deriv: usize,
        delta: f64,
    ) -> Result<Self, TGADomainError> {
        let df = self.frame.clone().collect()?;

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

        let mut new_df = df.clone();
        new_df.replace(col, Series::new(col.into(), filtered))?;

        self.frame = new_df.lazy();
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
        mut self,
        time_col: &str,
        new_time_col: &str,
        columns: &[&str],
        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let df = self.frame.clone().collect()?;

        // --- old time ---
        let t_old = extract_f64_column(&df, time_col)?;

        // --- new grid ---
        let t_new = uniform_grid_from(&t_old, n_points, kind)
            .map_err(|e| TGADomainError::InvalidOperation(e))?;

        let mut series: Vec<Column> = Vec::new();

        // --- resample requested columns ---
        for &col in columns {
            let y_old = extract_f64_column(&df, col)?;

            let (_, y_new) = spline_resample(&t_old, &y_old, n_points, kind)
                .map_err(|e| TGADomainError::InvalidOperation(e))?;

            series.push(Series::new(col.into(), y_new).into());
        }

        // --- add new time column ---
        series.push(Series::new(new_time_col.into(), t_new).into());

        let new_df = DataFrame::new(series)?;

        // --- update schema ---
        if self.schema.time.as_deref() == Some(time_col) {
            self.schema.time = Some(new_time_col.into());
        }

        self.frame = new_df.lazy();
        Ok(self)
    }
    /*
    pub fn spline_smooth_column(
        mut self,
        col: &str,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let time_col = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();

        // --- Collect ---
        let df = self.frame.clone().collect()?;

        let time = df.column(&time_col)?.f64().map_err(|_| {
            TGADomainError::InvalidOperation(format!("Time column '{}' must be f64", time_col))
        })?;

        let values = df.column(col)?.f64().map_err(|_| {
            TGADomainError::InvalidOperation(format!("Column '{}' must be f64", col))
        })?;

        if time.null_count() > 0 || values.null_count() > 0 {
            return Err(TGADomainError::InvalidOperation(format!(
                "Spline smoothing does not support nulls (column '{}')",
                col
            )));
        }

        let t: Vec<f64> = time.into_no_null_iter().collect();
        let y: Vec<f64> = values.into_no_null_iter().collect();

        if t.len() < 4 {
            return Err(TGADomainError::InvalidOperation(
                "Spline smoothing requires at least 4 points".into(),
            ));
        }

        // --- Check monotonicity ---
        for i in 1..t.len() {
            if t[i] <= t[i - 1] {
                return Err(TGADomainError::InvalidOperation(
                    "Time column must be strictly monotonic for spline smoothing".into(),
                ));
            }
        }

        // --- Build spline ---
        let interp = kind.to_interpolation();

        let keys = t
            .iter()
            .zip(y.iter())
            .map(|(&x, &v)| Key::<f64, f64>::new(x , v , interp))
            .collect();

        let spline = Spline::<f64, f64>::from_vec(keys);

        // --- Evaluate spline at original time points ---
        let mut smoothed = Vec::with_capacity(t.len());
        for &x in &t {
            match spline.sample(x ) {
                Some(v) => smoothed.push(v as f64),
                None => {
                    return Err(TGADomainError::InvalidOperation(
                        "Spline evaluation failed (out of bounds)".into(),
                    ));
                }
            }
        }

        // --- Write back ---
        let mut new_df = df.clone();
        new_df.replace(col, Series::new(col.into(), smoothed))?;

        self.frame = new_df.lazy();

        Ok(self)
    }

     /// onto a uniform time grid of size `n_points`.
    ///
    /// This operation:
    /// - rebuilds time column
    /// - resamples all f64 columns
    /// - preserves schema bindings
    ///
    pub fn spline_resample(
        mut self,
        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let time_col = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();

        let df = self.frame.clone().collect()?;

        let time_series = df.column(&time_col)?;
        let time = if matches!(time_series.dtype(), DataType::Float64) {
            let ts = time_series.f64()?;
            ts.clone()
        } else {
            let cast_series = time_series.cast(&DataType::Float64)?;
           let ts =  cast_series.f64()?;
           ts.clone()
        };
        if time.null_count() > 0 {
            return Err(TGADomainError::InvalidOperation(
                "Resampling does not support nulls in time column".into(),
            ));
        }

        let t: Vec<f64> = time.into_no_null_iter().collect();
        let t_min = *t.first().unwrap();
        let t_max = *t.last().unwrap();

        if n_points < 2 {
            return Err(TGADomainError::InvalidOperation(
                "n_points must be >= 2".into(),
            ));
        }

        if t.len() < 2 {
            return Err(TGADomainError::InvalidOperation(
                "Spline resampling requires at least 2 points".into(),
            ));
        }

        if matches!(kind, SplineKind::Cubic) && t.len() < 4 {
            return Err(TGADomainError::InvalidOperation(
                "Cubic spline resampling requires at least 4 points".into(),
            ));
        }

        for i in 1..t.len() {
            if t[i] <= t[i - 1] {
                return Err(TGADomainError::InvalidOperation(
                    "Time column must be strictly monotonic for spline resampling".into(),
                ));
            }
        }

        // New uniform grid
        let new_time: Vec<f64> = (0..n_points)
            .map(|i| {
                t_min + (t_max - t_min) * (i as f64) / ((n_points - 1) as f64)
            })
            .collect();

        let interp = kind.to_interpolation();

        let t_f64 = t;
        let new_time_f64: Vec<f64> = new_time.clone();

        let mut new_columns: Vec<Column> = Vec::new();

        for col in df.get_columns() {
            if *col.name() == time_col {
                new_columns.push(Series::new(col.name().clone(), &new_time).into());
                continue;
            }

            if let Ok(values) = col.f64() {
                if values.null_count() > 0 {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Resampling does not support nulls (column '{}')",
                        col.name()
                    )));
                }

                let y: Vec<f64> = values.into_no_null_iter().collect();
                let y_f64: Vec<f64> = y;

                let keys = t_f64
                    .iter()
                    .zip(y_f64.iter())
                    .map(|(&x, &v)| Key::new(x, v, interp))
                    .collect();

                let spline = Spline::from_vec(keys);

                let mut resampled = Vec::with_capacity(new_time_f64.len());
                for &x in &new_time_f64 {
                    match spline.sample(x) {
                        Some(v) => resampled.push(v as f64),
                        None => return Err(TGADomainError::InvalidOperation(
                            "Spline evaluation failed".into(),
                        )),
                    }
                }

                new_columns.push(Series::new(col.name().clone(), resampled).into());
            } else {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Resampling only supports f64 columns (column '{}')",
                    col.name()
                )));
            }
        }

        let new_df = DataFrame::new(new_columns)?;
        self.frame = new_df.lazy();

        Ok(self)
    }


     pub fn spline_resample2(
        mut self,
        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let time_col = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();

        let df = self.frame.clone().collect()?;

        let time = df.column(&time_col)?.f64()?;
        if time.null_count() > 0 {
            return Err(TGADomainError::InvalidOperation(
                "Resampling does not support nulls in time column".into(),
            ));
        }

        let t: Vec<f64> = time.into_no_null_iter().collect();
        let t_min = *t.first().unwrap();
        let t_max = *t.last().unwrap();

        if n_points < 2 {
            return Err(TGADomainError::InvalidOperation(
                "n_points must be >= 2".into(),
            ));
        }

        // New uniform grid
        let new_time: Vec<f64> = (0..n_points)
            .map(|i| {
                t_min + (t_max - t_min) * (i as f64) / ((n_points - 1) as f64)
            })
            .collect();

        let interp = kind.to_interpolation();

        let mut new_columns: Vec<Column> = Vec::new();
        for col in df.get_columns() {
            if *col.name() == time_col {
              new_columns.push(Series::new(col.name().clone(), &new_time).into());
                continue;
            }

            if let Ok(values) = col.f64() {
                if values.null_count() > 0 {
                    return Err(TGADomainError::InvalidOperation(format!(
                        "Resampling does not support nulls (column '{}')",
                        col.name()
                    )));
                }

                let y: Vec<f64> = values.into_no_null_iter().collect();

                let keys: Vec<Key<f64, f64>> = t
                    .iter()
                    .zip(y.iter())
                    .map(|(&x, &v)| Key::new(x, v, interp))
                    .collect();

                let spline = Spline::from_vec(keys);

                let resampled: Vec<f64> = new_time
                    .iter()
                    .map(|&x| {
                        spline.sample(x).ok_or_else(|| {
                            TGADomainError::InvalidOperation(
                                "Spline evaluation failed".into(),
                            )
                        })
                    })
                    .collect::<Result<_, _>>()?;

                new_columns.push(Series::new(col.name().clone(), resampled).into());
            } else {
             return Err(TGADomainError::InvalidOperation(format!(
                    "Resampling only supports f64 columns (column '{}')",
                    col.name()
                )));
            }
        }

        let new_df = DataFrame::new(new_columns)?;
        self.frame = new_df.lazy();

        Ok(self)
    }
     */
}
