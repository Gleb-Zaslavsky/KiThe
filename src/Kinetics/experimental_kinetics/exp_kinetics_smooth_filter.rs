use crate::Kinetics::experimental_kinetics::exp_kinetics_main::{TGADataset, TGADomainError};
use crate::Kinetics::experimental_kinetics::savgol2:: SG_filter_dyn;
use polars::prelude::*;
use std::f64;

/* TODO!:
SG как альтернатива derive_rate

unit-инференс для SG-deriv

SciPy-like SG (#2) как отдельный mode

lazy-реализация rolling_mean_many


*/



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
                    SGMode::FirPadding => self.sg_filter_column(
                        col,
                        *window,
                        *poly_order,
                        *deriv,
                        *delta,
                    )?,
                },

             
                | SmoothStrategy::Lowess { .. }
                | SmoothStrategy::Spline { .. } => {
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
}
