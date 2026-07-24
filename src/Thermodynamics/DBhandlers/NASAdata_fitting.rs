use crate::Thermodynamics::DBhandlers::NASAdata::{
    Coeffs, Cp, Cp_sym, NASAError, NASAdata, dh, dh_sym, ds, ds_sym,
};
use crate::Thermodynamics::DBhandlers::fitting_reports::compare_property_series;

use RustedSciThe::numerical::optimization::fitting_features::SewTwoFunctions;
use RustedSciThe::numerical::optimization::sym_fitting::Fitting;
use RustedSciThe::symbolic::symbolic_engine::Expr;

use prettytable::{Cell, Row, Table};
use std::collections::HashMap;

use nalgebra::{DMatrix, DVector};
////////////////////////////////////////STABLE LINEAR LEAST SQUARES FITTING FOR NASA-7 COEFFICIENTS//////////////////////////////////////
use std::f64::EPSILON;

/// Gas constant J / mol / K
const R: f64 = 8.314_462_618_15324_f64;
/// Result of a fit: coefficients and diagnostics
pub struct FitResult {
    /// coefficients in NASA form: a1..a7 (and optionally a8 as cp offset in Cp/R units)
    pub coeffs: Vec<f64>,
    pub stats: FitStats,
}

/// Diagnostics for fit quality
#[derive(Debug, Clone)]
pub struct FitStats {
    pub r2_cp: f64,
    pub r2_h: f64,
    pub r2_s: f64,
    pub rmse_cp: f64,
    pub rmse_h: f64,
    pub rmse_s: f64,
    pub maxabs_cp: f64,
    pub maxabs_h: f64,
    pub maxabs_s: f64,
    /// condition number of the weighted & column-scaled design matrix
    pub cond_number: f64,
    /// singular values (useful for debugging)
    pub singular_values: Vec<f64>,
}

/// High-level fit function:
/// - temps_K: temperatures (K)
/// - cp_si: Cp (J / mol / K)
/// - h_si: H (J / mol)
/// - s_si: S (J / mol / K)
/// - allow_cp_offset: if true returns 8 coefficients (last is Cp offset in Cp/R units)
/// - lambda: Tikhonov regularization parameter (recommended: 1e-10 .. 1e-6)
/// - returns FitResult
pub fn fit_nasa7_stable(
    temps_K: &[f64],
    cp_si: &[f64],
    h_si: &[f64],
    s_si: &[f64],
    allow_cp_offset: bool,
    lambda: f64,
) -> Result<FitResult, NASAError> {
    validate_fit_inputs("NASA", temps_K, cp_si, h_si, s_si)?;
    let n = temps_K.len();
    let cols = if allow_cp_offset { 8 } else { 7 };
    let rows = 3 * n;

    // --- Targets converted to NASA dimensionless form ---
    // cp_tgt = Cp / R
    // h_tgt  = H / (R * T)
    // s_tgt  = S / R
    let mut cp_tgt = vec![0.0; n];
    let mut h_tgt = vec![0.0; n];
    let mut s_tgt = vec![0.0; n];
    for i in 0..n {
        let t = temps_K[i];
        cp_tgt[i] = cp_si[i] / R;
        h_tgt[i] = h_si[i] / (R * t);
        s_tgt[i] = s_si[i] / R;
    }

    // --- row weights (1/std) with safety fallback ---
    let w_cp = 1.0 / (std2(&cp_tgt).max(1e-16));
    let w_h = 1.0 / (std2(&h_tgt).max(1e-16));
    let w_s = 1.0 / (std2(&s_tgt).max(1e-16));

    // --- Build design matrix X (rows: 3*n; cols: 7 or 8) and y ---
    // Row order per sample i:
    // row 3*i:   Cp/R  -> [1, T, T^2, T^3, T^4, (a6=0), (a7=0), (a8?)]
    // row 3*i+1: H/(R T) -> [1, T/2, T^2/3, T^3/4, T^4/5, 1/T, 0, 0]
    // row 3*i+2: S/R -> [ln T, T, T^2/2, T^3/3, T^4/4, 0, 1, 0]
    let mut X = DMatrix::<f64>::zeros(rows, cols);
    let mut y = DVector::<f64>::zeros(rows);

    for (i, &T) in temps_K.iter().enumerate() {
        let t = T;
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;

        let row_cp = 3 * i;
        let row_h = row_cp + 1;
        let row_s = row_cp + 2;

        // Cp/R row
        X[(row_cp, 0)] = 1.0; // a1
        X[(row_cp, 1)] = t; // a2 * T
        X[(row_cp, 2)] = t2; // a3 * T^2
        X[(row_cp, 3)] = t3; // a4 * T^3
        X[(row_cp, 4)] = t4; // a5 * T^4
        // a6 (1/T) not in Cp row
        // a7 (S constant) not in Cp row
        if allow_cp_offset {
            X[(row_cp, 7)] = 1.0; // a8 Cp offset
        }
        y[row_cp] = cp_tgt[i];

        // H/(R T) row
        X[(row_h, 0)] = 1.0; // a1
        X[(row_h, 1)] = t / 2.0; // a2 * T/2
        X[(row_h, 2)] = t2 / 3.0; // a3 * T^2/3
        X[(row_h, 3)] = t3 / 4.0; // a4 * T^3/4
        X[(row_h, 4)] = t4 / 5.0; // a5 * T^4/5
        if cols >= 6 {
            X[(row_h, 5)] = 1.0 / t;
        } // a6 / T term
        y[row_h] = h_tgt[i];

        // S/R row
        X[(row_s, 0)] = t.ln(); // a1 * ln T
        X[(row_s, 1)] = t; // a2 * T
        X[(row_s, 2)] = t2 / 2.0; // a3 * T^2/2
        X[(row_s, 3)] = t3 / 3.0; // a4 * T^3/3
        X[(row_s, 4)] = t4 / 4.0; // a5 * T^4/4
        if cols >= 7 {
            X[(row_s, 6)] = 1.0;
        } // a7 constant for S
        y[row_s] = s_tgt[i];
    }

    // --- Apply row weights: multiply rows of X and entries of y ---
    for i in 0..n {
        scale_row2(&mut X, &mut y, 3 * i, w_cp);
        scale_row2(&mut X, &mut y, 3 * i + 1, w_h);
        scale_row2(&mut X, &mut y, 3 * i + 2, w_s);
    }

    // --- Column scaling (RMS per column) to improve conditioning ---
    let mut col_scale = vec![1.0f64; cols];
    for j in 0..cols {
        // compute RMS of column j
        let mut sumsq = 0.0;
        for i in 0..rows {
            let v = X[(i, j)];
            sumsq += v * v;
        }
        let rms = ((sumsq / rows as f64).sqrt()).max(EPSILON);
        col_scale[j] = rms;
        // scale column
        if rms != 0.0 {
            for i in 0..rows {
                X[(i, j)] /= rms;
            }
        }
    }

    // --- SVD of scaled & weighted design matrix ---
    let svd = X.svd(true, true);
    let sv = svd.singular_values.clone();
    let mut singular_values: Vec<f64> = Vec::with_capacity(sv.len());
    for k in 0..sv.len() {
        singular_values.push(sv[k]);
    }

    // condition number
    let cond_number = if singular_values.len() > 0 {
        let max = singular_values[0].abs();
        let min = singular_values[singular_values.len() - 1]
            .abs()
            .max(EPSILON);
        max / min
    } else {
        f64::INFINITY
    };

    // --- Tikhonov regularized solution via SVD ---
    // X = U * S * V^T
    // theta_scaled = V * diag( s_i / (s_i^2 + lambda) ) * U^T * y
    let u_opt = svd
        .u
        .as_ref()
        .ok_or_else(|| NASAError::FittingError("SVD did not compute U".to_string()))?;
    let v_t_opt = svd
        .v_t
        .as_ref()
        .ok_or_else(|| NASAError::FittingError("SVD did not compute V^T".to_string()))?;

    // compute u^T * y  (u is rows x cols? but u has size rows x r, we use full U with same cols)
    // nalgebra's U has dimensions rows x min(rows,cols); Vt is min(rows,cols) x cols
    let u = u_opt;
    let v_t = v_t_opt;

    // Project y onto U: u_t_y = U^T * y  (length = r = sv.len())
    let r = sv.len();
    let mut u_t_y = DVector::<f64>::zeros(r);
    for i in 0..r {
        // column i of U is U.column(i)
        let mut accum = 0.0;
        for row in 0..rows {
            accum += u[(row, i)] * y[row];
        }
        u_t_y[i] = accum;
    }

    // build filtered coefficients in SVD space: d_i = s_i / (s_i^2 + lambda)
    let mut d = vec![0.0f64; r];
    for i in 0..r {
        let s_i = sv[i];
        d[i] = (s_i) / (s_i * s_i + lambda.max(0.0));
    }

    // compute z = d * (U^T y) elementwise
    let mut z = DVector::<f64>::zeros(r);
    for i in 0..r {
        z[i] = d[i] * u_t_y[i];
    }

    // theta_scaled = V * z  ; V has size cols x r (note v_t is r x cols)
    // so compute theta_scaled_j = sum_i V[j,i] * z_i  where V = (Vt)^T
    let mut theta_scaled = DVector::<f64>::zeros(cols);
    for j in 0..cols {
        let mut accum = 0.0;
        for i in 0..r {
            accum += v_t[(i, j)] * z[i];
        }
        theta_scaled[j] = accum;
    }

    // --- unscale coefficients (column scaling) ---
    let mut theta = vec![0.0f64; cols];
    for j in 0..cols {
        let scl = col_scale[j].max(EPSILON);
        theta[j] = theta_scaled[j] / scl;
    }

    // --- Prepare output coefficients (length 7 or 8) ---
    let coeffs = theta.clone();

    // --- Diagnostics: compute predictions in original SI units and stats ---
    let mut cp_pred = vec![0.0f64; n];
    let mut h_pred = vec![0.0f64; n];
    let mut s_pred = vec![0.0f64; n];
    for (i, &T) in temps_K.iter().enumerate() {
        cp_pred[i] = nasa7_cp_si(T, &coeffs);
        h_pred[i] = nasa7_h_si(T, &coeffs);
        s_pred[i] = nasa7_s_si(T, &coeffs);
    }

    let stats = FitStats {
        r2_cp: r2(cp_si, &cp_pred),
        r2_h: r2(h_si, &h_pred),
        r2_s: r2(s_si, &s_pred),
        rmse_cp: rmse(cp_si, &cp_pred),
        rmse_h: rmse(h_si, &h_pred),
        rmse_s: rmse(s_si, &s_pred),
        maxabs_cp: max_abs_err(cp_si, &cp_pred),
        maxabs_h: max_abs_err(h_si, &h_pred),
        maxabs_s: max_abs_err(s_si, &s_pred),
        cond_number,
        singular_values,
    };

    if coeffs.iter().any(|value| !value.is_finite())
        || !stats.r2_cp.is_finite()
        || !stats.r2_h.is_finite()
        || !stats.r2_s.is_finite()
        || !stats.rmse_cp.is_finite()
        || !stats.rmse_h.is_finite()
        || !stats.rmse_s.is_finite()
        || !stats.maxabs_cp.is_finite()
        || !stats.maxabs_h.is_finite()
        || !stats.maxabs_s.is_finite()
        || !stats.cond_number.is_finite()
        || stats.singular_values.iter().any(|value| !value.is_finite())
    {
        return Err(NASAError::FittingError(
            "NASA fitting produced non-finite coefficients or diagnostics".to_string(),
        ));
    }

    Ok(FitResult { coeffs, stats })
}

/// Check whether fit is acceptable using thresholds (adjust as desired)
pub fn fit_ok(stats: &FitStats) -> bool {
    // Strict defaults — may be relaxed for larger T ranges
    let r2_min = 0.999;
    let rmse_cp_max = 1e-3; // J/mol/K
    let rmse_h_max = 1e-2; // J/mol (enthalpy is larger absolute scale; use small value)
    let rmse_s_max = 1e-3; // J/mol/K

    let max_cp = 3e-3;
    let max_h = 5e-2;
    let max_s = 5e-3;

    stats.r2_cp >= r2_min
        && stats.r2_h >= r2_min
        && stats.r2_s >= r2_min
        && stats.rmse_cp <= rmse_cp_max
        && stats.rmse_h <= rmse_h_max
        && stats.rmse_s <= rmse_s_max
        && stats.maxabs_cp <= max_cp
        && stats.maxabs_h <= max_h
        && stats.maxabs_s <= max_s
}

/// --- NASA formula helpers returning SI units (J/mol or J/mol/K) ---

/// Cp in J / mol / K from NASA coefficients (a1..a7) ; if a8 present, it's Cp offset in Cp/R
pub fn nasa7_cp_si(T: f64, a: &[f64]) -> f64 {
    let mut cp_r = 0.0;
    cp_r += a[0];
    cp_r += a[1] * T;
    cp_r += a[2] * T * T;
    cp_r += a[3] * T * T * T;
    cp_r += a[4] * T * T * T * T;
    if a.len() >= 8 {
        cp_r += a[7]; // cp offset in Cp/R
    }
    cp_r * R
}

/// H in J / mol
pub fn nasa7_h_si(T: f64, a: &[f64]) -> f64 {
    let t = T;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let h_rt =
        a[0] + a[1] * t / 2.0 + a[2] * t2 / 3.0 + a[3] * t3 / 4.0 + a[4] * t4 / 5.0 + a[5] / t;
    h_rt * R * t
}

/// S in J / mol / K
pub fn nasa7_s_si(T: f64, a: &[f64]) -> f64 {
    let t = T;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let s_r = a[0] * t.ln() + a[1] * t + a[2] * t2 / 2.0 + a[3] * t3 / 3.0 + a[4] * t4 / 4.0 + a[6];
    s_r * R
}

/// ----------------- small utility functions -----------------

fn std2(v: &[f64]) -> f64 {
    let n = v.len() as f64;
    if n <= 1.0 {
        return 0.0;
    }
    let mean = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    var.sqrt()
}

fn scale_row2(X: &mut DMatrix<f64>, y: &mut DVector<f64>, row: usize, w: f64) {
    let cols = X.ncols();
    for j in 0..cols {
        X[(row, j)] *= w;
    }
    y[row] *= w;
}

fn r2(obs: &[f64], pred: &[f64]) -> f64 {
    let n = obs.len();
    if n == 0 {
        return 0.0;
    }
    let mean = obs.iter().sum::<f64>() / n as f64;
    let ss_tot: f64 = obs.iter().map(|v| (v - mean).powi(2)).sum();
    let ss_res: f64 = obs
        .iter()
        .zip(pred.iter())
        .map(|(o, p)| (o - p).powi(2))
        .sum();
    if ss_tot.abs() < 1e-20 {
        return 1.0;
    }
    1.0 - ss_res / ss_tot
}

fn rmse(obs: &[f64], pred: &[f64]) -> f64 {
    let n = obs.len() as f64;
    let mse = obs
        .iter()
        .zip(pred.iter())
        .map(|(o, p)| (o - p).powi(2))
        .sum::<f64>()
        / n;
    mse.sqrt()
}

fn max_abs_err(obs: &[f64], pred: &[f64]) -> f64 {
    obs.iter()
        .zip(pred.iter())
        .map(|(o, p)| (o - p).abs())
        .fold(0.0, f64::max)
}

fn validate_fit_inputs(
    label: &str,
    temps_K: &[f64],
    cp_si: &[f64],
    h_si: &[f64],
    s_si: &[f64],
) -> Result<(), NASAError> {
    if temps_K.len() != cp_si.len() || temps_K.len() != h_si.len() || temps_K.len() != s_si.len() {
        return Err(NASAError::FittingError(format!(
            "{} fitting input lengths must match: temps={}, cp={}, h={}, s={}",
            label,
            temps_K.len(),
            cp_si.len(),
            h_si.len(),
            s_si.len()
        )));
    }

    let validate_slice = |name: &str, values: &[f64]| -> Result<(), NASAError> {
        if let Some((idx, value)) = values
            .iter()
            .enumerate()
            .find(|(_, value)| !value.is_finite())
        {
            return Err(NASAError::FittingError(format!(
                "{} fitting input contains non-finite {} at index {}: {}",
                label, name, idx, value
            )));
        }
        Ok(())
    };

    validate_slice("temperature", temps_K)?;
    validate_slice("Cp", cp_si)?;
    validate_slice("H", h_si)?;
    validate_slice("S", s_si)?;
    Ok(())
}

fn validate_temperature_bounds(label: &str, T_min: f64, T_max: f64) -> Result<(), NASAError> {
    if !T_min.is_finite() || !T_max.is_finite() {
        return Err(NASAError::FittingError(format!(
            "{} fitting bounds must be finite, got T_min={}, T_max={}",
            label, T_min, T_max
        )));
    }
    if T_min >= T_max {
        return Err(NASAError::FittingError(format!(
            "{} fitting bounds must satisfy T_min < T_max, got T_min={}, T_max={}",
            label, T_min, T_max
        )));
    }
    Ok(())
}

fn validate_coeffs_finite(label: &str, coeffs: &[f64]) -> Result<(), NASAError> {
    if let Some((idx, value)) = coeffs
        .iter()
        .enumerate()
        .find(|(_, value)| !value.is_finite())
    {
        return Err(NASAError::FittingError(format!(
            "{} coefficients must be finite, coefficient {} = {}",
            label, idx, value
        )));
    }
    Ok(())
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fit NASA-7 coefficients (a1..a7) using stacked linear least squares.
///
/// Inputs:
///  - temps_K: temperatures in Kelvin
///  - cp_si: Cp (J / mol / K)
///  - h_si:  H  (J / mol)
///  - s_si:  S  (J / mol / K)
///  - allow_cp_offset: if true, fit an 8th coefficient (cp offset) added to Cp only
///
/// Returns: Vec<f64> of length 7 (or 8 if allow_cp_offset==true)
pub fn fit_nasa7(
    temps_K: &[f64],
    cp_si: &[f64],
    h_si: &[f64],
    s_si: &[f64],
    allow_cp_offset: bool,
) -> Result<Vec<f64>, NASAError> {
    validate_fit_inputs("NASA", temps_K, cp_si, h_si, s_si)?;

    let n = temps_K.len();
    let p = if allow_cp_offset { 8 } else { 7 }; // columns

    // Gas constant (J / mol / K)
    const R: f64 = 1.987;

    // Convert targets to NASA dimensionless form:
    // cp_tgt = Cp / R
    // h_tgt  = H / (R * T)
    // s_tgt  = S / R
    let mut cp_tgt = Vec::with_capacity(n);
    let mut h_tgt = Vec::with_capacity(n);
    let mut s_tgt = Vec::with_capacity(n);

    for i in 0..n {
        let t = temps_K[i];
        cp_tgt.push(cp_si[i] / R);
        // avoid dividing by zero (T should always be > 0)
        h_tgt.push(h_si[i] / (R * t));
        s_tgt.push(s_si[i] / R);
    }

    // Weight determination: 1 / stddev per quantity
    let w_cp = 1.0 / std(&cp_tgt);
    let w_h = 0.1 / std(&h_tgt);
    let w_s = 1.0 / std(&s_tgt);

    let rows = 3 * n;
    let cols = p;
    let mut X = DMatrix::<f64>::zeros(rows, cols);
    let mut y = DVector::<f64>::zeros(rows);

    for (i, &T) in temps_K.iter().enumerate() {
        let t = T;
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;

        let row_cp = 3 * i;
        let row_h = 3 * i + 1;
        let row_s = 3 * i + 2;

        // Cp/R row: a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
        // if allow_cp_offset: add column a8 (last column) with coefficient 1
        X[(row_cp, 0)] = 1.0; // a1
        X[(row_cp, 1)] = t; // a2
        X[(row_cp, 2)] = t2; // a3
        X[(row_cp, 3)] = t3; // a4
        X[(row_cp, 4)] = t4; // a5
        if allow_cp_offset {
            // columns: [a1,a2,a3,a4,a5,a6,a7,a8]
            X[(row_cp, 7)] = 1.0; // a8 (Cp-only offset)
        }
        // a6 (enthalpy 1/T term) and a7 (S constant) are zero in Cp row
        // Fill zeros for missing columns implicitly (matrix was zeroed).
        y[row_cp] = cp_tgt[i];

        // H/(R T) row:
        // a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6 / T
        X[(row_h, 0)] = 1.0; // a1
        X[(row_h, 1)] = t / 2.0; // a2 * T/2
        X[(row_h, 2)] = t2 / 3.0; // a3 * T^2/3
        X[(row_h, 3)] = t3 / 4.0; // a4 * T^3/4
        X[(row_h, 4)] = t4 / 5.0; // a5 * T^4/5
        X[(row_h, 5)] = 1.0 / t; // a6 * (1/T)
        // a7, a8 (if present) are zero for H row
        y[row_h] = h_tgt[i];

        // S/R row:
        // a1*ln T + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
        X[(row_s, 0)] = t.ln(); // a1 * ln T
        X[(row_s, 1)] = t; // a2 * T
        X[(row_s, 2)] = t2 / 2.0; // a3 * T^2/2
        X[(row_s, 3)] = t3 / 3.0; // a4 * T^3/3
        X[(row_s, 4)] = t4 / 4.0; // a5 * T^4/4
        // a6 appears only in H row
        // a7 is the additive constant for S:
        X[(row_s, 6)] = 1.0; // a7
        y[row_s] = s_tgt[i];
    }

    // Scale rows by weights
    for i in 0..n {
        scale_row(&mut X, &mut y, 3 * i, w_cp);
        scale_row(&mut X, &mut y, 3 * i + 1, w_h);
        scale_row(&mut X, &mut y, 3 * i + 2, w_s);
    }

    // Solve via SVD
    let svd = X.svd(true, true);
    let theta = svd
        .solve(&y, 1e-12)
        .map_err(|e| NASAError::FittingError(format!("SVD solve failed: {}", e)))?;

    Ok(theta.data.as_vec().clone())
}

/// compute standard deviation (population std)
fn std(v: &[f64]) -> f64 {
    let n = v.len() as f64;
    let mean = v.iter().copied().sum::<f64>() / n;
    let var = v.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    var.sqrt()
}

fn scale_row(X: &mut DMatrix<f64>, y: &mut DVector<f64>, row: usize, w: f64) {
    if w == 1.0 {
        return;
    }
    for c in 0..X.ncols() {
        X[(row, c)] *= w;
    }
    y[row] *= w;
}

pub(crate) type Nasa7Coefficients = (f64, f64, f64, f64, f64, f64, f64);

/// Returns the shared NASA-7 vector when adjacent ranges describe the same polynomial.
///
/// Re-fitting an already global polynomial is both unnecessary and numerically harmful:
/// the NASA basis is ill-conditioned in raw Kelvin powers, so an optimizer may move far
/// in coefficient space while changing the sampled property values only slightly.
fn identical_adjacent_coefficients(
    coeffs_min: &Coeffs,
    coeffs_max: &Coeffs,
) -> Option<Nasa7Coefficients> {
    (coeffs_min.coeff == coeffs_max.coeff).then_some(coeffs_min.coeff)
}

/// Extracts a finite, complete result from RustedSciThe's LM fitting wrapper.
///
/// `Fitting::solve()` records an unsuccessful termination by leaving its result
/// empty. Convert that state into the typed thermodynamics error expected by the
/// public fitting methods instead of panicking on an `Option::unwrap()`.
fn lm_solution(
    fitting: &Fitting,
    stage: &str,
    expected_coefficients: &[&str],
) -> Result<(HashMap<String, f64>, f64), NASAError> {
    let solution = fitting
        .solution_map()
        .ok_or_else(|| NASAError::FittingError(format!("{stage} LM solver did not converge")))?;
    let r_squared = fitting.r_squared().ok_or_else(|| {
        NASAError::FittingError(format!("{stage} LM solver did not produce R squared"))
    })?;

    if !r_squared.is_finite() {
        return Err(NASAError::FittingError(format!(
            "{stage} LM solver produced non-finite R squared: {r_squared}"
        )));
    }
    for &name in expected_coefficients {
        let value = solution.get(name).ok_or_else(|| {
            NASAError::FittingError(format!(
                "{stage} LM solution is missing coefficient `{name}`"
            ))
        })?;
        if !value.is_finite() {
            return Err(NASAError::FittingError(format!(
                "{stage} LM solution produced non-finite coefficient `{name}`: {value}"
            )));
        }
    }

    Ok((solution, r_squared))
}

impl NASAdata {
    /// Builds quality reports for a known NASA-7 vector without running another fit.
    pub(crate) fn reports_for_coefficients(
        &self,
        temps: &[f64],
        cp_vals: &[f64],
        h_vals: &[f64],
        s_vals: &[f64],
        coeffs: Nasa7Coefficients,
    ) -> Result<HashMap<String, FittingReport>, NASAError> {
        let (a, b, c, d, e, f, g) = coeffs;
        let mut reports = HashMap::new();
        reports.insert(
            "Cp fitting report".to_string(),
            self.direct_compare2(cp_vals, temps, Cp_sym(a, b, c, d, e))?,
        );
        reports.insert(
            "dh fitting report".to_string(),
            self.direct_compare2(h_vals, temps, dh_sym(a, b, c, d, e, f))?,
        );
        reports.insert(
            "ds fitting report".to_string(),
            self.direct_compare2(s_vals, temps, ds_sym(a, b, c, d, e, g))?,
        );
        Ok(reports)
    }

    /// Generates thermodynamic data from two adjacent NASA coefficient ranges
    ///
    /// Creates temperature points and calculates Cp, H, S values using original coefficients
    /// from two adjacent ranges. Used as reference data for fitting new coefficients.
    ///
    /// # Arguments
    /// * `coeffs_min` - NASA coefficients for lower temperature range
    /// * `coeffs_max` - NASA coefficients for upper temperature range  
    /// * `T_min` - Minimum temperature for data generation [K]
    /// * `T_max` - Maximum temperature for data generation [K]
    ///
    /// # Returns
    /// Tuple of (temperatures, Cp_values, H_values, S_values)
    pub fn prepare_fit_data_for_2_range_fitting(
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), NASAError> {
        validate_temperature_bounds("NASA adjacent", T_min, T_max)?;
        let (_, T_center) = (coeffs_min.T.0, coeffs_min.T.1);
        let (T_center_check, _) = (coeffs_max.T.0, coeffs_max.T.1);
        if T_center != T_center_check {
            return Err(NASAError::InvalidTemperatureRange);
        }
        if coeffs_min.T.0 >= coeffs_min.T.1 || coeffs_max.T.0 >= coeffs_max.T.1 {
            return Err(NASAError::FittingError(
                "NASA interval bounds must be strictly ordered".to_string(),
            ));
        }
        validate_coeffs_finite(
            "NASA lower interval",
            &[
                coeffs_min.coeff.0,
                coeffs_min.coeff.1,
                coeffs_min.coeff.2,
                coeffs_min.coeff.3,
                coeffs_min.coeff.4,
                coeffs_min.coeff.5,
                coeffs_min.coeff.6,
            ],
        )?;
        validate_coeffs_finite(
            "NASA upper interval",
            &[
                coeffs_max.coeff.0,
                coeffs_max.coeff.1,
                coeffs_max.coeff.2,
                coeffs_max.coeff.3,
                coeffs_max.coeff.4,
                coeffs_max.coeff.5,
                coeffs_max.coeff.6,
            ],
        )?;

        let a = coeffs_min.coeff.0;
        let b = coeffs_min.coeff.1;
        let c = coeffs_min.coeff.2;
        let d = coeffs_min.coeff.3;
        let e = coeffs_min.coeff.4;
        let f = coeffs_min.coeff.5;
        let g = coeffs_min.coeff.6;

        let mut temps = (T_min as usize..=T_center as usize)
            .step_by(5 as usize)
            .map(|t| t as f64)
            .collect::<Vec<_>>();
        let mut cp_vals = temps
            .iter()
            .map(|&T| Cp(T, a, b, c, d, e))
            .collect::<Vec<_>>();
        let mut h_vals = temps
            .iter()
            .map(|&T| dh(T, a, b, c, d, e, f))
            .collect::<Vec<_>>();
        let mut s_vals = temps
            .iter()
            .map(|&T| ds(T, a, b, c, d, e, g))
            .collect::<Vec<_>>();
        let a = coeffs_max.coeff.0;
        let b = coeffs_max.coeff.1;
        let c = coeffs_max.coeff.2;
        let d = coeffs_max.coeff.3;
        let e = coeffs_max.coeff.4;
        let f = coeffs_max.coeff.5;
        let g = coeffs_max.coeff.6;

        let temps2 = (T_center as usize + 5..=T_max as usize)
            .step_by(5 as usize)
            .map(|t| t as f64)
            .collect::<Vec<_>>();

        temps.extend(temps2.clone());
        let cp_vals2 = temps2
            .clone()
            .iter()
            .map(|&T| Cp(T, a, b, c, d, e))
            .collect::<Vec<_>>();
        let h_vals2 = temps2
            .clone()
            .iter()
            .map(|&T| dh(T, a, b, c, d, e, f))
            .collect::<Vec<_>>();
        let s_vals2 = temps2
            .iter()
            .map(|&T| ds(T, a, b, c, d, e, g))
            .collect::<Vec<_>>();
        cp_vals.extend(cp_vals2);
        h_vals.extend(h_vals2);
        s_vals.extend(s_vals2);
        // println!("temps = {:?}", temps);
        // println!("cp_vals = {:?}", cp_vals);
        // println!("h_vals = {:?}", h_vals);
        // println!("s_vals = {:?}", s_vals);

        assert_eq!(temps.len(), cp_vals.len());
        assert_eq!(temps.len(), h_vals.len());
        assert_eq!(temps.len(), s_vals.len());
        Ok((temps, cp_vals, h_vals, s_vals))
    }

    /// Generates thermodynamic data from multiple non-adjacent NASA coefficient ranges
    ///
    /// Combines data from multiple temperature ranges that may have gaps between them.
    /// Only includes data points within the specified T_min to T_max bounds.
    ///
    /// # Arguments
    /// * `coeffs` - Vector of NASA coefficients for different temperature ranges
    /// * `T_min` - Minimum temperature bound [K]
    /// * `T_max` - Maximum temperature bound [K]
    ///
    /// # Returns
    /// Tuple of (temperatures, Cp_values, H_values, S_values)
    fn prepare_fit_data_for_non_adjacent_fitting(
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), NASAError> {
        validate_temperature_bounds("NASA non-adjacent", T_min, T_max)?;
        let mut Tvec: Vec<f64> = Vec::new();
        let mut cp_vec: Vec<f64> = Vec::new();
        let mut h_vec: Vec<f64> = Vec::new();
        let mut s_vec: Vec<f64> = Vec::new();

        for coeff in coeffs.iter() {
            if coeff.T.0 >= coeff.T.1 {
                return Err(NASAError::FittingError(format!(
                    "NASA interval bounds must be strictly ordered, got {:?}",
                    coeff.T
                )));
            }
            validate_coeffs_finite(
                "NASA interval",
                &[
                    coeff.coeff.0,
                    coeff.coeff.1,
                    coeff.coeff.2,
                    coeff.coeff.3,
                    coeff.coeff.4,
                    coeff.coeff.5,
                    coeff.coeff.6,
                ],
            )?;
            let (T_min_r, T_max_r) = (coeff.T.0, coeff.T.1);
            let mut Tleft = T_min_r;
            let mut Tright = T_max_r;

            // Skip if range is completely outside T_min..T_max
            if T_min_r < T_min && T_max_r <= T_min {
                continue;
            }
            if T_min_r >= T_max && T_max_r > T_max {
                continue;
            }

            // Adjust range boundaries to fit within T_min..T_max
            if T_min_r <= T_min && T_min <= T_max_r {
                Tleft = T_min;
                Tright = T_max_r;
            }
            if T_min_r <= T_max && T_max <= T_max_r {
                Tright = T_max;
                Tleft = T_min_r;
            }
            if T_min_r >= T_min && T_max >= T_max_r {
                Tleft = T_min_r;
                Tright = T_max_r;
            }

            let (a, b, c, d, e, f, g) = coeff.coeff;
            let temps = (Tleft as usize + 5..=Tright as usize)
                .step_by(5)
                .map(|t| t as f64)
                .collect::<Vec<_>>();
            let cp_vals = temps
                .iter()
                .map(|&T| Cp(T, a, b, c, d, e))
                .collect::<Vec<_>>();
            let h_vals = temps
                .iter()
                .map(|&T| dh(T, a, b, c, d, e, f))
                .collect::<Vec<_>>();
            let s_vals = temps
                .iter()
                .map(|&T| ds(T, a, b, c, d, e, g))
                .collect::<Vec<_>>();

            Tvec.extend(temps);
            cp_vec.extend(cp_vals);
            h_vec.extend(h_vals);
            s_vec.extend(s_vals);
        }

        Ok((Tvec, cp_vec, h_vec, s_vec))
    }
    /// Fits NASA-7 coefficients using simultaneous linear least squares
    ///
    /// Uses the `fit_nasa7` function to simultaneously fit all 7 NASA coefficients
    /// by solving a stacked linear system for Cp, H, and S equations.
    ///
    /// # Arguments
    /// * `coeffs_min` - NASA coefficients for lower temperature range
    /// * `coeffs_max` - NASA coefficients for upper temperature range
    /// * `T_min` - Minimum fitting temperature [K]
    /// * `T_max` - Maximum fitting temperature [K]
    ///
    /// # Returns
    /// Ok(()) if fitting succeeds, stores coefficients in self.coeffs
    pub fn fit_cp_dh_ds(
        &mut self,
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<(), NASAError> {
        let preserved_coeffs = identical_adjacent_coefficients(&coeffs_min, &coeffs_max);
        let (temps, cp_vals, h_vals, s_vals) =
            Self::prepare_fit_data_for_2_range_fitting(coeffs_min, coeffs_max, T_min, T_max)?;

        // Identical source polynomials are already the exact global solution. Preserve
        // their coefficients bit-for-bit instead of sending them through an SVD solve.
        if let Some(coeffs) = preserved_coeffs {
            self.coeffs = coeffs;
            let report_map =
                self.reports_for_coefficients(&temps, &cp_vals, &h_vals, &s_vals, coeffs)?;
            Self::print_quality_table(&report_map);
            return Ok(());
        }

        let coeffs = fit_nasa7(&temps, &cp_vals, &h_vals, &s_vals, false)?;
        if coeffs.len() != 7 {
            return Err(NASAError::FittingError(format!(
                "NASA solver returned {} coefficients instead of 7",
                coeffs.len()
            )));
        }
        let a = coeffs[0];
        let b = coeffs[1];
        let c = coeffs[2];
        let d = coeffs[3];
        let e = coeffs[4];
        let f = coeffs[5];
        let g = coeffs[6];

        let coeffs = (a, b, c, d, e, f, g);
        self.coeffs = coeffs;
        let report_map =
            self.reports_for_coefficients(&temps, &cp_vals, &h_vals, &s_vals, coeffs)?;
        Self::print_quality_table(&report_map);
        Ok(())
    }

    pub fn fitting_cp_dh_ds_non_adjacent(
        &mut self,
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NASAError> {
        let mut Tvec: Vec<f64> = Vec::new();
        let mut cp_vec: Vec<f64> = Vec::new();
        let mut h_vec: Vec<f64> = Vec::new();
        let mut s_vec: Vec<f64> = Vec::new();

        for coeffs in coeffs.iter() {
            let (T_min_r, T_max_r) = (coeffs.T.0, coeffs.T.1);
            let mut Tleft = T_min_r;
            let mut Tright = T_max_r;

            // Skip if range is completely outside T_min..T_max
            if T_min_r < T_min && T_max_r <= T_min {
                continue;
            }
            if T_min_r >= T_max && T_max_r > T_max {
                continue;
            }

            // Adjust range boundaries to fit within T_min..T_max
            if T_min_r <= T_min && T_min <= T_max_r {
                Tleft = T_min;
                Tright = T_max_r;
            }
            if T_min_r <= T_max && T_max <= T_max_r {
                Tright = T_max;
                Tleft = T_min_r;
            }
            if T_min_r >= T_min && T_max >= T_max_r {
                Tleft = T_min_r;
                Tright = T_max_r;
            }

            let a = coeffs.coeff.0;
            let b = coeffs.coeff.1;
            let c = coeffs.coeff.2;
            let d = coeffs.coeff.3;
            let e = coeffs.coeff.4;
            let f = coeffs.coeff.5;
            let g = coeffs.coeff.6;

            let temps = (Tleft as usize + 5..=Tright as usize)
                .step_by(5 as usize)
                .map(|t| t as f64)
                .collect::<Vec<_>>();
            let cp_vals = temps
                .iter()
                .map(|&T| Cp(T, a, b, c, d, e))
                .collect::<Vec<_>>();
            let h_vals = temps
                .iter()
                .map(|&T| dh(T, a, b, c, d, e, f))
                .collect::<Vec<_>>();
            let s_vals = temps
                .iter()
                .map(|&T| ds(T / 1000.0, a, b, c, d, e, g))
                .collect::<Vec<_>>();
            Tvec.extend(temps);
            cp_vec.extend(cp_vals);
            h_vec.extend(h_vals);
            s_vec.extend(s_vals);
        }
        if Tvec.len() != cp_vec.len() || Tvec.len() != h_vec.len() || Tvec.len() != s_vec.len() {
            return Err(NASAError::FittingError(
                "NASA non-adjacent fitting assembled inconsistent sample lengths".to_string(),
            ));
        }

        //  let coeffs = fit_nasa7(&Tvec, &cp_vec, &h_vec, &s_vec, false);
        let coeffs = fit_nasa7_stable(&Tvec, &cp_vec, &h_vec, &s_vec, false, 1e-9)?.coeffs;
        if coeffs.len() != 7 {
            return Err(NASAError::FittingError(format!(
                "NASA solver returned {} coefficients instead of 7",
                coeffs.len()
            )));
        }
        let a = coeffs[0];
        let b = coeffs[1];
        let c = coeffs[2];
        let d = coeffs[3];
        let e = coeffs[4];
        let f = coeffs[5];
        let g = coeffs[6];

        let coeffs = (a, b, c, d, e, f, g);
        self.coeffs = coeffs;
        let report_map = self.reports_for_coefficients(&Tvec, &cp_vec, &h_vec, &s_vec, coeffs)?;
        Self::print_quality_table(&report_map);

        Ok(report_map)
    }

    ///////////////////////////////////////////////////////////////////

    /// Stacked weighted LM fitting for simultaneous optimization
    ///
    /// Fits one NASA-7 vector against separate `Cp`, `dH`, and `dS` residual
    /// families. Each family is scaled to a comparable numerical magnitude, but
    /// it remains an independent residual block in the LM objective.
    ///
    /// Keeping the residuals separate gives LM direct information about each
    /// property instead of only their weighted scalar sum. The existing
    /// property-wise reports already rejected inaccurate final fits; the stacked
    /// objective improves the optimization formulation itself.
    ///
    /// # Arguments
    /// * `coeffs_min` - Coefficients for lower temperature range
    /// * `coeffs_max` - Coefficients for upper temperature range
    /// * `T_min` - Lower bound of fitting temperature range
    /// * `T_max` - Upper bound of fitting temperature range
    ///
    /// # Returns
    /// * `Ok(HashMap<String, FittingReport>)` - Quality report for combined fitting
    /// * `Err(NASAError)` - If temperature ranges are not adjacent or fitting fails
    pub fn fitting_adjacent3(
        &mut self,
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NASAError> {
        let preserved_coeffs = identical_adjacent_coefficients(&coeffs_min, &coeffs_max);
        let (_, T_center) = (coeffs_min.T.0, coeffs_min.T.1);
        let (T_center_check, _) = (coeffs_max.T.0, coeffs_max.T.1);
        if T_center != T_center_check {
            return Err(NASAError::InvalidTemperatureRange);
        }
        let (temps, cp_vals, h_vals, s_vals) = Self::prepare_fit_data_for_2_range_fitting(
            coeffs_min.clone(),
            coeffs_max.clone(),
            T_min,
            T_max,
        )?;
        if let Some(coeffs) = preserved_coeffs {
            self.coeffs = coeffs;
            let report_map =
                self.reports_for_coefficients(&temps, &cp_vals, &h_vals, &s_vals, coeffs)?;
            Self::print_quality_table(&report_map);
            return Ok(report_map);
        }
        // Scale the three residual families without collapsing them into one value.
        let w_cp = 1.0;
        let w_h = 1e-3;
        let w_s = 1e-1;
        let Cp_func_to_fit =
            Expr::parse_expression("1.987* (a + b * t + c * t^2 + d * t^3 + e * t^4)");
        let dh_func_to_fit = Expr::parse_expression(
            "1.987 *t* (a + b * t/2 + (c/3) * t^2 + (d/4) * t^3 + (e/5) * t^4 + f/t)",
        );
        let ds_func_to_fit = Expr::parse_expression(
            "1.987 * (a* ln( t ) + b * t + c * t^2 / 2 + d * t^3 / 3.0 + e * t^4 / 4.0 + g)",
        );
        let weighted_equations = vec![
            Expr::Const(w_cp) * Cp_func_to_fit.clone(),
            Expr::Const(w_h) * dh_func_to_fit.clone(),
            Expr::Const(w_s) * ds_func_to_fit.clone(),
        ];
        // RustedSciThe stores multi-equation targets in equation-major order.
        let mut weighted_targets = cp_vals.iter().map(|value| w_cp * value).collect::<Vec<_>>();
        weighted_targets.extend(h_vals.iter().map(|value| w_h * value));
        weighted_targets.extend(s_vals.iter().map(|value| w_s * value));
        let initial_guess = vec![
            (coeffs_min.coeff.0 + coeffs_max.coeff.0) / 2.0,
            (coeffs_min.coeff.1 + coeffs_max.coeff.1) / 2.0,
            (coeffs_min.coeff.2 + coeffs_max.coeff.2) / 2.0,
            (coeffs_min.coeff.3 + coeffs_max.coeff.3) / 2.0,
            (coeffs_min.coeff.4 + coeffs_max.coeff.4) / 2.0,
            (coeffs_min.coeff.5 + coeffs_max.coeff.5) / 2.0,
            (coeffs_min.coeff.6 + coeffs_max.coeff.6) / 2.0,
        ];
        let unknowns = vec![
            "a".to_string(),
            "b".to_string(),
            "c".to_string(),
            "d".to_string(),
            "e".to_string(),
            "f".to_string(),
            "g".to_string(),
        ];
        let fitting = Fitting::new()
            .with_data(temps.clone(), weighted_targets)
            .with_equations(weighted_equations)
            .with_unknowns(unknowns)
            .with_arg("t".to_string())
            .with_initial_guess(initial_guess)
            .with_tolerance(1e-10)
            .with_f_tolerance(1e-12)
            .with_g_tolerance(1e-12)
            .with_max_iterations(1000)
            .build();
        let (map_of_solutions, r_squared) = lm_solution(
            &fitting,
            "NASA adjacent simultaneous",
            &["a", "b", "c", "d", "e", "f", "g"],
        )?;
        println!("r_squared for stacked fitting: {r_squared}");
        let a = map_of_solutions["a"];
        let b = map_of_solutions["b"];
        let c = map_of_solutions["c"];
        let d = map_of_solutions["d"];
        let e = map_of_solutions["e"];
        let f = map_of_solutions["f"];
        let g = map_of_solutions["g"];

        let coeffs = (a, b, c, d, e, f, g);
        self.coeffs = coeffs;
        let report_map =
            self.reports_for_coefficients(&temps, &cp_vals, &h_vals, &s_vals, coeffs)?;
        Self::print_quality_table(&report_map);

        Ok(report_map)
    }

    /// Sequential fitting for multiple non-adjacent temperature ranges
    ///
    /// Uses a 2-stage approach: first fits entropy (dS) to determine shared coefficients (a,b,c,d,e,g),
    /// then fits enthalpy coefficient f separately. Handles multiple temperature ranges with gaps.
    ///
    /// # Arguments
    /// * `coeffs` - Vector of NASA coefficients from different temperature ranges
    /// * `T_min` - Minimum fitting temperature [K]
    /// * `T_max` - Maximum fitting temperature [K]
    ///
    /// # Returns
    /// HashMap with fitting quality reports for each thermodynamic property
    pub fn fitting_non_adjacent(
        &mut self,
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NASAError> {
        let (temps, cp_vals, h_vals, s_vals) =
            Self::prepare_fit_data_for_non_adjacent_fitting(coeffs.clone(), T_min, T_max)?;

        /////////////////////////////fitting ds//////////////////////////////////

        let ds_func_to_fit = Expr::parse_expression(
            "1.987 * (a* ln( t ) + b * t + c * t^2 / 2 + d * t^3 / 3.0 + e * t^4 / 4.0 + g)",
        );

        let mut sew = SewTwoFunctions {
            f1: Expr::Const(0.0),
            f2: Expr::Const(0.0),
            x_left: 0.0,
            x_central: 0.0,
            n_points: 0,
            x_right: 0.0,
            fitting_data: (temps.clone(), s_vals.clone()),
            fitting: Fitting::new(),
        };

        let mut initial_guess = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        for coeff in coeffs.iter() {
            initial_guess[0] += coeff.coeff.0;
            initial_guess[1] += coeff.coeff.1;
            initial_guess[2] += coeff.coeff.2;
            initial_guess[3] += coeff.coeff.3;
            initial_guess[4] += coeff.coeff.4;

            initial_guess[5] += coeff.coeff.6;
        }
        initial_guess = initial_guess
            .iter()
            .map(|x| x / coeffs.len() as f64)
            .collect();
        sew.fit(
            ds_func_to_fit.clone(),
            Some(vec![
                "a".to_string(),
                "b".to_string(),
                "c".to_string(),
                "d".to_string(),
                "e".to_string(),
                "g".to_string(),
            ]),
            "t".to_string(),
            initial_guess,
            None,
            None,
            None,
            None,
            None,
        );

        let (map_of_solutions_6_coeffs, r_squared) = lm_solution(
            &sew.fitting,
            "NASA non-adjacent entropy",
            &["a", "b", "c", "d", "e", "g"],
        )?;
        println!("r_squared for s: {r_squared}");
        if 1.0 - r_squared >= 5.0e-2 {
            return Err(NASAError::FittingError(format!(
                "Poor entropy fitting quality: R² = {:.6}",
                r_squared
            )));
        }
        let a = map_of_solutions_6_coeffs["a"];
        let b = map_of_solutions_6_coeffs["b"];
        let c = map_of_solutions_6_coeffs["c"];
        let d = map_of_solutions_6_coeffs["d"];
        let e = map_of_solutions_6_coeffs["e"];
        let g = map_of_solutions_6_coeffs["g"];
        /////////////////// fitting h ///////////////////////
        let dh_func_to_fit = Expr::parse_expression(
            "1.987 *t* (a + b * t/2 + (c/3) * t^2 + (d/4) * t^3 + (e/5) * t^4 + f/t)",
        );
        let dh_func_to_fit = dh_func_to_fit
            .clone()
            .set_variable_from_map(&map_of_solutions_6_coeffs);
        let mut sew2 = SewTwoFunctions {
            f1: Expr::Const(0.0),
            f2: Expr::Const(0.0),
            x_left: 0.0,
            x_central: 0.0,
            n_points: 0,
            x_right: 0.0,
            fitting_data: (temps.clone(), h_vals.clone()),
            fitting: Fitting::new(),
        };
        let mut initial_guess = vec![0.0; 1];
        for coeff in coeffs.iter() {
            initial_guess[0] += coeff.coeff.5;
        }
        initial_guess = initial_guess
            .iter()
            .map(|x| x / coeffs.len() as f64)
            .collect();
        sew2.fit(
            dh_func_to_fit.clone(),
            Some(vec!["f".to_string()]),
            "t".to_string(),
            initial_guess,
            None,
            None,
            None,
            None,
            None,
        );
        let (f_map, r_squared) = lm_solution(&sew2.fitting, "NASA non-adjacent enthalpy", &["f"])?;
        println!("r_squared for h: {r_squared}");
        if 1.0 - r_squared >= 5.0e-2 {
            return Err(NASAError::FittingError(format!(
                "Poor enthalpy fitting quality: R² = {:.6}",
                r_squared
            )));
        }
        let f = f_map["f"];
        let coeffs = (a, b, c, d, e, f, g);
        self.coeffs = coeffs;
        // Validate the published coefficients against every raw property family.
        let report_map =
            self.reports_for_coefficients(&temps, &cp_vals, &h_vals, &s_vals, coeffs)?;
        Self::print_quality_table(&report_map);

        Ok(report_map.clone())
    }

    /// Sequential fitting for adjacent temperature ranges using entropy-first approach
    ///
    /// 3-stage fitting process:
    /// 1. Fit entropy (dS) to determine shared coefficients (a,b,c,d,e,g)
    /// 2. Fit enthalpy coefficient f using fixed shared coefficients
    /// 3. Validate heat capacity fitting with final coefficients
    ///
    /// # Arguments
    /// * `coeffs_min` - NASA coefficients for lower temperature range
    /// * `coeffs_max` - NASA coefficients for upper temperature range
    /// * `T_min` - Minimum fitting temperature [K]
    /// * `T_max` - Maximum fitting temperature [K]
    ///
    /// # Returns
    /// HashMap with fitting quality reports for each thermodynamic property
    pub fn fitting_adjacent2(
        &mut self,
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NASAError> {
        let preserved_coeffs = identical_adjacent_coefficients(&coeffs_min, &coeffs_max);
        let (_, T_center) = (coeffs_min.T.0, coeffs_min.T.1);
        let (T_center_check, _) = (coeffs_max.T.0, coeffs_max.T.1);
        if T_center != T_center_check {
            return Err(NASAError::InvalidTemperatureRange);
        }
        let (temps, cp_vals, h_vals, s_vals) = Self::prepare_fit_data_for_2_range_fitting(
            coeffs_min.clone(),
            coeffs_max.clone(),
            T_min,
            T_max,
        )?;
        if let Some(coeffs) = preserved_coeffs {
            self.coeffs = coeffs;
            let report_map =
                self.reports_for_coefficients(&temps, &cp_vals, &h_vals, &s_vals, coeffs)?;
            Self::print_quality_table(&report_map);
            return Ok(report_map);
        }
        /////////////////////////////fitting ds//////////////////////////////////

        let ds_func_to_fit = Expr::parse_expression(
            "1.987 * (a* ln( t ) + b * t + c * t^2 / 2 + d * t^3 / 3.0 + e * t^4 / 4.0 + g)",
        );

        let mut sew = SewTwoFunctions {
            f1: Expr::Const(0.0),
            f2: Expr::Const(0.0),
            x_left: 0.0,
            x_central: 0.0,
            n_points: 0,
            x_right: 0.0,
            fitting_data: (temps.clone(), s_vals.clone()),
            fitting: Fitting::new(),
        };

        let initial_guess = vec![
            (coeffs_min.coeff.0 + coeffs_max.coeff.0) / 2.0,
            (coeffs_min.coeff.1 + coeffs_max.coeff.1) / 2.0,
            (coeffs_min.coeff.2 + coeffs_max.coeff.2) / 2.0,
            (coeffs_min.coeff.3 + coeffs_max.coeff.3) / 2.0,
            (coeffs_min.coeff.4 + coeffs_max.coeff.4) / 2.0,
            (coeffs_min.coeff.6 + coeffs_max.coeff.6) / 2.0,
        ];
        sew.fit(
            ds_func_to_fit.clone(),
            Some(vec![
                "a".to_string(),
                "b".to_string(),
                "c".to_string(),
                "d".to_string(),
                "e".to_string(),
                "g".to_string(),
            ]),
            "t".to_string(),
            initial_guess,
            None,
            None,
            None,
            None,
            None,
        );

        let (map_of_solutions_6_coeffs, r_squared) = lm_solution(
            &sew.fitting,
            "NASA adjacent sequential entropy",
            &["a", "b", "c", "d", "e", "g"],
        )?;
        println!("r_squared for s: {r_squared}");
        if 1.0 - r_squared >= 5.0e-2 {
            return Err(NASAError::FittingError(format!(
                "Poor entropy fitting quality: R² = {:.6}",
                r_squared
            )));
        }
        let a = map_of_solutions_6_coeffs["a"];
        let b = map_of_solutions_6_coeffs["b"];
        let c = map_of_solutions_6_coeffs["c"];
        let d = map_of_solutions_6_coeffs["d"];
        let e = map_of_solutions_6_coeffs["e"];
        let g = map_of_solutions_6_coeffs["g"];
        /////////////////// fitting h ///////////////////////
        let dh_func_to_fit = Expr::parse_expression(
            "1.987 *t* (a + b * t/2 + (c/3) * t^2 + (d/4) * t^3 + (e/5) * t^4 + f/t)",
        );
        let dh_func_to_fit = dh_func_to_fit
            .clone()
            .set_variable_from_map(&map_of_solutions_6_coeffs);
        let mut sew2 = SewTwoFunctions {
            f1: Expr::Const(0.0),
            f2: Expr::Const(0.0),
            x_left: 0.0,
            x_central: 0.0,
            n_points: 0,
            x_right: 0.0,
            fitting_data: (temps.clone(), h_vals.clone()),
            fitting: Fitting::new(),
        };
        sew2.fit(
            dh_func_to_fit.clone(),
            Some(vec!["f".to_string()]),
            "t".to_string(),
            vec![(coeffs_min.coeff.5 + coeffs_max.coeff.5) / 2.0],
            None,
            None,
            None,
            None,
            None,
        );
        let (f_map, r_squared) =
            lm_solution(&sew2.fitting, "NASA adjacent sequential enthalpy", &["f"])?;
        println!("r_squared for h: {r_squared}");
        if 1.0 - r_squared >= 5.0e-2 {
            return Err(NASAError::FittingError(format!(
                "Poor enthalpy fitting quality: R² = {:.6}",
                r_squared
            )));
        }
        let f = f_map["f"];
        let coeffs = (a, b, c, d, e, f, g);
        self.coeffs = coeffs;
        // Validate the published coefficients against every raw property family.
        let report_map =
            self.reports_for_coefficients(&temps, &cp_vals, &h_vals, &s_vals, coeffs)?;

        // Print fitting quality report table
        Self::print_quality_table(&report_map);

        Ok(report_map)
    }
    /// Compares fitted function against reference data and generates quality metrics
    ///
    /// Evaluates the fitted symbolic function at given temperature points and compares
    /// with reference data to compute L1/L2 norms, maximum relative error, and identify
    /// points with significant deviations (>1% relative error).
    ///
    /// # Arguments
    /// * `y_data` - Reference data values (Cp, H, or S)
    /// * `T_vec` - Temperature points [K]
    /// * `fit_function` - Fitted symbolic function to evaluate
    ///
    /// # Returns
    /// FittingReport with error norms and deviation statistics
    pub fn direct_compare2(
        &self,
        y_data: &[f64],
        T_vec: &[f64],
        fit_function: Expr,
    ) -> Result<FittingReport, NASAError> {
        let report = compare_property_series(y_data, T_vec, fit_function)
            .map_err(|err| NASAError::FittingError(err.to_string()))?;
        Ok(FittingReport {
            l2_norm: report.l2_norm,
            max_norm: report.max_norm,
            number_of_significant_points: report.number_of_significant_points,
            T_range: report.t_range,
        })
    }

    /// Prints a formatted table of fitting quality metrics
    ///
    /// Displays L1 norm, L2 norm, max norm, number of significant deviation points,
    /// and temperature range where significant deviations occur for each fitted property.
    ///
    /// # Arguments
    /// * `report_map` - HashMap containing fitting reports for different properties
    pub fn print_quality_table(report_map: &HashMap<String, FittingReport>) {
        // Print fitting quality report table
        let mut table = Table::new();
        table.add_row(Row::new(vec![
            Cell::new("Property"),
            Cell::new("L2 Norm"),
            Cell::new("Max Norm"),
            Cell::new("Points"),
            Cell::new("T Range (K)"),
        ]));

        for (property, report) in report_map {
            let t_range_str = match report.T_range {
                Some((min, max)) => format!("{:.1}-{:.1}", min, max),
                None => "N/A".to_string(),
            };

            table.add_row(Row::new(vec![
                Cell::new(property),
                Cell::new(&format!("{:.2e}", report.l2_norm)),
                Cell::new(&format!("{:.2e}", report.max_norm)),
                Cell::new(&report.number_of_significant_points.to_string()),
                Cell::new(&t_range_str),
            ]));
        }

        println!("\nFitting Quality Report:");
        table.printstd();
    }
}

/// Quality metrics for NASA coefficient fitting
///
/// Contains error norms and statistics about fitting quality,
/// including identification of temperature ranges with poor fits.
#[derive(Debug, Clone)]
pub struct FittingReport {
    /// L2 (RMS) norm of relative errors
    pub l2_norm: f64,
    /// Maximum relative error encountered
    pub max_norm: f64,
    /// Number of points with >1% relative error
    pub number_of_significant_points: usize,
    /// Temperature range where significant deviations occur [K]
    pub T_range: Option<(f64, f64)>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_fit_nasa7_rejects_bad_inputs() {
        let temps = [300.0, 400.0, 500.0];
        let cp = [10.0, 11.0];
        let h = [1.0, 2.0, 3.0];
        let s = [4.0, 5.0, 6.0];

        let err = fit_nasa7(&temps, &cp, &h, &s, false).unwrap_err();
        assert!(matches!(err, NASAError::FittingError(_)));

        let cp_bad = [10.0, f64::NAN, 12.0];
        let err = fit_nasa7(&temps, &cp_bad, &h, &s, false).unwrap_err();
        assert!(matches!(err, NASAError::FittingError(_)));
    }

    #[test]
    fn test_fit_nasa7_stable_recovers_known_coeffs() {
        let coeffs = [2.5, 1.0e-3, -2.0e-6, 3.0e-9, -4.0e-12, 8.0, 1.2];
        let temps = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0];
        let cp: Vec<f64> = temps.iter().map(|&t| nasa7_cp_si(t, &coeffs)).collect();
        let h: Vec<f64> = temps.iter().map(|&t| nasa7_h_si(t, &coeffs)).collect();
        let s: Vec<f64> = temps.iter().map(|&t| nasa7_s_si(t, &coeffs)).collect();

        let fitted =
            fit_nasa7_stable(&temps, &cp, &h, &s, false, 1e-10).expect("stable fit should succeed");
        assert_eq!(fitted.coeffs.len(), 7);
        assert!(fit_ok(&fitted.stats));
        for (expected, actual) in coeffs.iter().zip(fitted.coeffs.iter()) {
            assert_relative_eq!(*expected, *actual, epsilon = 1e-2);
        }
    }

    #[test]
    fn test_prepare_fit_data_for_2_range_fitting_keeps_boundary_once() {
        let coeffs_min = Coeffs {
            T: (300.0, 500.0),
            coeff: (2.0, 1.0e-3, 0.0, 0.0, 0.0, 5.0, 1.0),
        };
        let coeffs_max = Coeffs {
            T: (500.0, 900.0),
            coeff: (2.0, 1.0e-3, 0.0, 0.0, 0.0, 5.0, 1.0),
        };

        let (temps, cp, h, s) =
            NASAdata::prepare_fit_data_for_2_range_fitting(coeffs_min, coeffs_max, 300.0, 700.0)
                .expect("boundary-aware preparation should succeed");

        assert_eq!(temps.len(), cp.len());
        assert_eq!(temps.len(), h.len());
        assert_eq!(temps.len(), s.len());
        assert_eq!(temps.iter().filter(|&&t| t == 500.0).count(), 1);
        assert!(temps.windows(2).all(|window| window[0] < window[1]));
    }

    #[test]
    fn test_lm_solution_converts_missing_backend_result_into_typed_error() {
        let fitting = Fitting::new();
        let error = lm_solution(&fitting, "test stage", &["a"])
            .expect_err("an unsolved fitting backend must not look successful");

        assert!(
            matches!(error, NASAError::FittingError(message) if message.contains("did not converge"))
        );
    }
}
