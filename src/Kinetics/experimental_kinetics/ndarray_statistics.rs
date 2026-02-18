//! # Statistical Analysis and ndarray Integration
//!
//! This module provides bridge between Polars-based TGADataset and ndarray for statistical operations.
//!
//! ## Main Data Structures
//!
//! - `ColumnStats`: Statistical summary (mean, variance, std_dev, min, max) for a single column
//! - `CorrelationMatrix`: Pearson correlation matrix for multiple columns
//! - `StatisticalAnalysis`: Container aggregating column statistics, correlation, and covariance
//!
//! ## Mathematical Operations
//!
//! - **Mean**: μ = (1/n) Σ xᵢ
//! - **Variance**: σ² = (1/n) Σ (xᵢ - μ)²
//! - **Standard Deviation**: σ = √(σ²)
//! - **Pearson Correlation**: r = Cov(X,Y) / (σₓ σᵧ)
//! - **Covariance Matrix**: Cov = (1/n) XᵀX where X is centered
//!
//! ## Usage
//!
//! ```rust,ignore
//! // Extract single column as ndarray
//! let mass_array = dataset.column_as_array1("mass")?;
//!
//! // Compute statistics for multiple columns
//! let stats = dataset.compute_statistics(&["mass", "temperature"], true)?;
//! stats.pretty_print();
//!
//! // Pearson correlation between two columns
//! let corr = dataset.pearson_correlation("mass", "temperature")?;
//!
//! // Individual statistics
//! let mean = dataset.column_mean("mass")?;
//! let std_dev = dataset.column_std_dev("mass")?;
//! ```

use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{TGADataset, TGADomainError};

use ndarray::{Array1, Array2, Axis};
use ndarray_stats::CorrelationExt;
use ndarray_stats::interpolate::Interpolate;
use polars::prelude::*;
use std::collections::HashMap;
impl TGADataset {
    /// Extract a column as Array1<f64>
    /// Contract:
    /// - column must exist
    /// - column must be f64
    /// - column must contain NO nulls
    pub fn column_as_array1(&self, col: &str) -> Result<Array1<f64>, TGADomainError> {
        let df = self.frame.clone().collect()?;

        let series = df
            .column(col)
            .map_err(|_| TGADomainError::ColumnNotFound(col.into()))?;

        let s = series.f64().map_err(TGADomainError::PolarsError)?;

        if s.null_count() > 0 {
            return Err(TGADomainError::InvalidOperation(format!(
                "Column '{}' contains nulls",
                col
            )));
        }

        let v: Vec<f64> = s.into_no_null_iter().collect();
        Ok(Array1::from_vec(v))
    }
    /// Extract two columns as Array2<f64> with shape (n_rows, 2)
    pub fn columns_as_array2(
        &self,
        col_x: &str,
        col_y: &str,
    ) -> Result<Array2<f64>, TGADomainError> {
        let x = self.column_as_array1(col_x)?;
        let y = self.column_as_array1(col_y)?;

        if x.len() != y.len() {
            return Err(TGADomainError::InvalidOperation(
                "Columns have different lengths".into(),
            ));
        }

        let n = x.len();
        let mut arr = Array2::<f64>::zeros((n, 2));

        for i in 0..n {
            arr[[i, 0]] = x[i];
            arr[[i, 1]] = y[i];
        }

        Ok(arr)
    }

    /// Insert or replace a column from Array1<f64>
    pub fn with_array1_column(
        mut self,
        col_name: &str,
        data: Array1<f64>,
    ) -> Result<Self, TGADomainError> {
        let series = Series::new(col_name.into(), data.to_vec());
        self.frame = self.frame.with_column(lit(series));
        Ok(self)
    }

    /// Extract multiple columns as Array2<f64> with shape (n_rows, n_cols)
    pub fn columns_as_array2_multi(&self, cols: &[&str]) -> Result<Array2<f64>, TGADomainError> {
        if cols.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "No columns specified".into(),
            ));
        }

        let arrays: Result<Vec<_>, _> =
            cols.iter().map(|&col| self.column_as_array1(col)).collect();
        let arrays = arrays?;

        let n_rows = arrays[0].len();
        let n_cols = arrays.len();

        for arr in &arrays {
            if arr.len() != n_rows {
                return Err(TGADomainError::InvalidOperation(
                    "Columns have different lengths".into(),
                ));
            }
        }

        let mut result = Array2::<f64>::zeros((n_rows, n_cols));
        for (j, arr) in arrays.iter().enumerate() {
            for i in 0..n_rows {
                result[[i, j]] = arr[i];
            }
        }

        Ok(result)
    }
}

//=================================================================
// STATISTICAL DATA STRUCTURES
//=================================================================

/// Statistical summary for a single column
#[derive(Debug, Clone)]
pub struct ColumnStats {
    pub name: String,
    pub mean: f64,
    pub variance: f64,
    pub std_dev: f64,
    pub min: f64,
    pub max: f64,
    pub count: usize,
}

/// Correlation matrix result
#[derive(Debug, Clone)]
pub struct CorrelationMatrix {
    pub columns: Vec<String>,
    pub matrix: Array2<f64>,
}

/// Statistical analysis container
#[derive(Debug, Clone)]
pub struct StatisticalAnalysis {
    pub column_stats: HashMap<String, ColumnStats>,
    pub correlation: Option<CorrelationMatrix>,
    pub covariance: Option<Array2<f64>>,
}

//=================================================================
// STATISTICAL OPERATIONS
//=================================================================

impl ColumnStats {
    pub fn from_array(name: String, data: &Array1<f64>) -> Self {
        let n = data.len();
        let mean = data.sum() / n as f64;
        let variance = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
        let std_dev = variance.sqrt();
        let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        Self {
            name,
            mean,
            variance,
            std_dev,
            min,
            max,
            count: n,
        }
    }
}

impl StatisticalAnalysis {
    pub fn new() -> Self {
        Self {
            column_stats: HashMap::new(),
            correlation: None,
            covariance: None,
        }
    }

    pub fn compute_from_dataset(
        dataset: &TGADataset,
        cols: &[&str],
        compute_correlation: bool,
    ) -> Result<Self, TGADomainError> {
        let mut analysis = Self::new();

        // Compute individual column statistics
        for &col in cols {
            let arr = dataset.column_as_array1(col)?;
            let stats = ColumnStats::from_array(col.to_string(), &arr);
            analysis.column_stats.insert(col.to_string(), stats);
        }

        // Compute correlation and covariance if requested
        if compute_correlation && cols.len() > 1 {
            let data = dataset.columns_as_array2_multi(cols)?;

            // Pearson correlation
            let corr = data.t().pearson_correlation().map_err(|e| {
                TGADomainError::InvalidOperation(format!("Correlation failed: {:?}", e))
            })?;

            analysis.correlation = Some(CorrelationMatrix {
                columns: cols.iter().map(|s| s.to_string()).collect(),
                matrix: corr,
            });

            // Covariance
            let cov = compute_covariance(&data);
            analysis.covariance = Some(cov);
        }

        Ok(analysis)
    }

    pub fn pretty_print(&self) {
        println!("\n=== Statistical Analysis ===\n");

        for (name, stats) in &self.column_stats {
            println!("Column: {}", name);
            println!("  Mean:     {:.6}", stats.mean);
            println!("  Std Dev:  {:.6}", stats.std_dev);
            println!("  Variance: {:.6}", stats.variance);
            println!("  Min:      {:.6}", stats.min);
            println!("  Max:      {:.6}", stats.max);
            println!("  Count:    {}", stats.count);
            println!();
        }

        if let Some(corr) = &self.correlation {
            println!("Pearson Correlation Matrix:");
            println!("Columns: {:?}", corr.columns);
            println!("{:.4}", corr.matrix);
            println!();
        }

        if let Some(cov) = &self.covariance {
            println!("Covariance Matrix:");
            println!("{:.4}", cov);
        }
    }
}

/// Compute covariance matrix from data matrix (n_rows, n_cols)
fn compute_covariance(data: &Array2<f64>) -> Array2<f64> {
    let n_rows = data.nrows() as f64;
    let n_cols = data.ncols();

    // Compute means
    let means = data.mean_axis(Axis(0)).unwrap();

    // Center the data
    let mut centered = data.clone();
    for i in 0..data.nrows() {
        for j in 0..n_cols {
            centered[[i, j]] -= means[j];
        }
    }

    // Compute covariance: (X^T * X) / n
    let cov = centered.t().dot(&centered) / n_rows;
    cov
}

//=================================================================
// CONVENIENCE METHODS ON TGADataset
//=================================================================

impl TGADataset {
    /// Compute statistics for specified columns
    pub fn compute_statistics(
        &self,
        cols: &[&str],
        include_correlation: bool,
    ) -> Result<StatisticalAnalysis, TGADomainError> {
        StatisticalAnalysis::compute_from_dataset(self, cols, include_correlation)
    }

    /// Compute Pearson correlation between two columns
    pub fn pearson_correlation(&self, col_x: &str, col_y: &str) -> Result<f64, TGADomainError> {
        let data = self.columns_as_array2(col_x, col_y)?;
        let corr = data.t().pearson_correlation().map_err(|e| {
            TGADomainError::InvalidOperation(format!("Correlation failed: {:?}", e))
        })?;
        Ok(corr[[0, 1]])
    }

    /// Compute mean of a column
    pub fn column_mean(&self, col: &str) -> Result<f64, TGADomainError> {
        let arr = self.column_as_array1(col)?;
        Ok(arr.sum() / arr.len() as f64)
    }

    /// Compute variance of a column
    pub fn column_variance(&self, col: &str) -> Result<f64, TGADomainError> {
        let arr = self.column_as_array1(col)?;
        let mean = arr.sum() / arr.len() as f64;
        let variance = arr.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / arr.len() as f64;
        Ok(variance)
    }

    /// Compute standard deviation of a column
    pub fn column_std_dev(&self, col: &str) -> Result<f64, TGADomainError> {
        Ok(self.column_variance(col)?.sqrt())
    }

    /// Linear interpolation: given x_old, y_old, and x_new, compute y_new
    /// Linear interpolation: given x_old, y_old, and x_new, compute y_new.
    /// This implements a small, dependency-free linear interpolator with
    /// simple extrapolation using the first/last interval.
    pub fn interpolate_linear(
        &self,
        x_col: &str,
        y_col: &str,
        x_new: &Array1<f64>,
    ) -> Result<Array1<f64>, TGADomainError> {
        let x_old = self.column_as_array1(x_col)?;
        let y_old = self.column_as_array1(y_col)?;

        if x_old.len() != y_old.len() {
            return Err(TGADomainError::InvalidOperation(
                "x and y columns have different lengths".into(),
            ));
        }

        let n = x_old.len();
        if n == 0 {
            return Err(TGADomainError::InvalidOperation(
                "Empty data for interpolation".into(),
            ));
        }

        let xs: Vec<f64> = x_old.to_vec();
        let ys: Vec<f64> = y_old.to_vec();

        let mut out = Array1::<f64>::zeros(x_new.len());

        for (idx, &xq) in x_new.iter().enumerate() {
            // Before first point -> extrapolate using first interval or constant if single point
            if xq <= xs[0] {
                if n == 1 {
                    out[idx] = ys[0];
                } else {
                    let (x0, x1) = (xs[0], xs[1]);
                    let (y0, y1) = (ys[0], ys[1]);
                    let t = (xq - x0) / (x1 - x0);
                    out[idx] = y0 + t * (y1 - y0);
                }
                continue;
            }

            // After last point -> extrapolate using last interval or constant if single point
            if xq >= xs[n - 1] {
                if n == 1 {
                    out[idx] = ys[0];
                } else {
                    let (x0, x1) = (xs[n - 2], xs[n - 1]);
                    let (y0, y1) = (ys[n - 2], ys[n - 1]);
                    let t = (xq - x0) / (x1 - x0);
                    out[idx] = y0 + t * (y1 - y0);
                }
                continue;
            }

            // Binary search for interval [lo, lo+1] where xs[lo] <= xq < xs[lo+1]
            let mut lo = 0usize;
            let mut hi = n - 1;
            while lo + 1 < hi {
                let mid = (lo + hi) / 2;
                if xs[mid] <= xq {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }

            let (x0, x1) = (xs[lo], xs[lo + 1]);
            let (y0, y1) = (ys[lo], ys[lo + 1]);
            let t = (xq - x0) / (x1 - x0);
            out[idx] = y0 + t * (y1 - y0);
        }

        Ok(out)
    }
}

//=================================================================
// TESTS
//=================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset_test::tests::ds_from_csv;
    use crate::Kinetics::experimental_kinetics::testing_mod::{
        NoiseModel, SpikeModel, VirtualTGA, VirtualTGAConfig,
    };
    fn make_test_dataset() -> TGADataset {
        let cfg = VirtualTGAConfig {
            n_points: 1000,
            dt: 1.0,
            temperature: 400.0,
            temp_noise: NoiseModel { sigma: 2.0 },
            m0: 100.0,
            k: 1e-3,
            mass_noise: NoiseModel { sigma: 0.5 },
            spikes: Some(SpikeModel {
                probability: 0.002,
                amplitude: 20.0,
            }),
            seed: 42,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();

        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();

        let ds = ds_from_csv(&csv);
        // Prevent the NamedTempFile from being deleted — the dataset uses a lazy
        // CSV reader that needs the file to remain available during collection.
        std::mem::forget(csv);

        // --- Full pipeline ---
        let ds = ds
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap();
        ds
    }

    #[test]
    fn test_column_as_array1() {
        let ds = make_test_dataset();
        let mass = ds.column_as_array1("mass").unwrap();
        assert_eq!(mass.len(), 1000);
        assert!(mass.iter().all(|&x| x.is_finite()));
    }

    #[test]
    fn test_columns_as_array2() {
        let ds = make_test_dataset();
        let arr = ds.columns_as_array2("time", "mass").unwrap();
        assert_eq!(arr.shape(), &[1000, 2]);
    }

    #[test]
    fn test_columns_as_array2_multi() {
        let ds = make_test_dataset();
        let arr = ds
            .columns_as_array2_multi(&["time", "temperature", "mass"])
            .unwrap();
        assert_eq!(arr.shape(), &[1000, 3]);
    }

    #[test]
    fn test_column_stats() {
        let ds = make_test_dataset();
        let mean = ds.column_mean("mass").unwrap();
        let variance = ds.column_variance("mass").unwrap();
        let std_dev = ds.column_std_dev("mass").unwrap();

        assert!(mean > 0.0 && mean < 200.0);
        assert!(variance > 0.0);
        assert!((std_dev - variance.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_pearson_correlation() {
        let ds = make_test_dataset();
        let corr = ds.pearson_correlation("time", "mass").unwrap();
        assert!(corr >= -1.0 && corr <= 1.0);
        assert!(corr < -0.5); // mass decreases with time
    }

    #[test]
    fn test_compute_statistics() {
        let ds = make_test_dataset();
        let stats = ds
            .compute_statistics(&["time", "mass", "temperature"], true)
            .unwrap();

        assert_eq!(stats.column_stats.len(), 3);
        assert!(stats.correlation.is_some());
        assert!(stats.covariance.is_some());

        let corr_matrix = stats.correlation.unwrap();
        assert_eq!(corr_matrix.matrix.shape(), &[3, 3]);

        // Diagonal should be 1.0
        for i in 0..3 {
            assert!((corr_matrix.matrix[[i, i]] - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_with_array1_column() {
        let ds = make_test_dataset();
        let new_data = Array1::from_vec(vec![1.0; 1000]);
        let ds2 = ds.with_array1_column("test_col", new_data).unwrap();

        let retrieved = ds2.column_as_array1("test_col").unwrap();
        assert_eq!(retrieved.len(), 1000);
        assert!(retrieved.iter().all(|&x| (x - 1.0).abs() < 1e-10));
    }

    #[test]
    fn test_interpolate_linear() {
        let ds = make_test_dataset();
        let x_new = Array1::linspace(0.0, 999.0, 100);
        let y_new = ds.interpolate_linear("time", "mass", &x_new).unwrap();

        assert_eq!(y_new.len(), 100);
        assert!(y_new.iter().all(|&x| x.is_finite()));
    }

    #[test]
    fn test_interpolate_linear_extrapolate() {
        let ds = make_test_dataset();
        let n = ds.column_as_array1("time").unwrap().len();

        // Create exact linear mass: mass = 2*time + 1
        let new_mass = Array1::from_vec((0..n).map(|i| 2.0 * (i as f64) + 1.0).collect());
        let ds2 = ds.with_array1_column("mass", new_mass).unwrap();

        let x_new = Array1::from_vec(vec![-1.0, 0.0, 0.5, (n - 1) as f64, n as f64]);
        let y_new = ds2.interpolate_linear("time", "mass", &x_new).unwrap();

        let expected = vec![
            2.0 * -1.0 + 1.0,
            2.0 * 0.0 + 1.0,
            2.0 * 0.5 + 1.0,
            2.0 * ((n - 1) as f64) + 1.0,
            2.0 * (n as f64) + 1.0,
        ];
        for (a, b) in y_new.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn test_interpolate_linear_single_point() {
        let ds = make_test_dataset();
        let n = ds.column_as_array1("time").unwrap().len();

        // Constant mass column
        let new_mass = Array1::from_vec(vec![42.0; n]);
        let ds2 = ds.with_array1_column("mass", new_mass).unwrap();

        let x_new = Array1::from_vec(vec![0.0, 10.0, 1000.0]);
        let y_new = ds2.interpolate_linear("time", "mass", &x_new).unwrap();

        assert!(y_new.iter().all(|&v| (v - 42.0).abs() < 1e-10));
    }

    #[test]
    fn test_statistical_analysis_pretty_print() {
        let ds = make_test_dataset();
        let stats = ds
            .compute_statistics(&["mass", "temperature"], true)
            .unwrap();
        stats.pretty_print(); // Should not panic
    }
}
