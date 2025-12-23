//! # NIST Thermodynamic Data Fitting Module
//!
//! This module provides advanced fitting algorithms for NIST thermodynamic polynomial coefficients.
//! It handles the complex mathematical relationships between heat capacity (Cp), enthalpy (dH),
//! and entropy (dS) using symbolic mathematics and sequential optimization.
//!
//! ## Two Main Fitting Algorithms
//!
//! ### 1. Adjacent Range Fitting (`fitting_adjacent_weighted`)
//! - **Use case**: Two adjacent temperature ranges that share a boundary point
//! - **Method**: Sequential two-stage optimization to avoid conflicting constraints
//!   - Stage 1: Fit Cp + weighted dS to determine shared coefficients (a,b,c,d,e,g)
//!   - Stage 2: Fit dH separately to determine enthalpy-specific coefficients (f,h)
//! - **Advantage**: Respects thermodynamic relationships, avoids optimization conflicts
//!
//! ### 2. Non-Adjacent Range Fitting (`fitting_non_adjacent`)
//! - **Use case**: Multiple temperature ranges with gaps between them
//! - **Method**: Same sequential approach as adjacent fitting but handles multiple ranges
//!   - Combines data from all overlapping ranges within target temperature bounds
//!   - Uses identical weighting scheme for consistency
//! - **Advantage**: Handles complex multi-range scenarios while maintaining accuracy
//!
//! ## NIST Polynomial Forms
//! The fitting targets these standard NIST equations:
//! - **Cp(T)**: `a + b*t + c*t² + d*t³ + e*t⁻²` (where t = T/1000)
//! - **dH(T)**: `a*t + b*t²/2 + c*t³/3 + d*t⁴/4 - e/t + f - h`
//! - **dS(T)**: `a*ln(t) + b*t + c*t²/2 + d*t³/3 - e/(2*t²) + g`
//!
//! ## Key Features
//! - Symbolic mathematics for exact coefficient relationships
//! - Weighted optimization to balance different thermodynamic properties
//! - Comprehensive error reporting and quality metrics
//! - Support for both simultaneous and sequential fitting approaches

use crate::Thermodynamics::DBhandlers::NISTdata::*;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};
use prettytable::{Cell, Row, Table};
use std::collections::HashMap;

/// Compute global NIST-form coefficients using SVD least squares fitting.
///
/// This function simultaneously fits all three thermodynamic properties (Cp, dH, dS)
/// using a weighted least squares approach with SVD decomposition.
///
/// # Arguments
/// * `temps_K` - Temperature sampling points in Kelvin
/// * `cp_orig` - Heat capacity values [J/mol·K]
/// * `h_orig` - Enthalpy values [kJ/mol]
/// * `s_orig` - Entropy values [J/mol·K]
///
/// # Returns
/// Vector of 8 NIST coefficients [a,b,c,d,e,f,g,h]
///
/// # Mathematical Approach
/// Constructs a system of equations where each temperature point contributes 3 rows:
/// - Cp equation: `a + b*t + c*t² + d*t³ + e*t⁻² + h = Cp`
/// - dH equation: `a*t + b*t²/2 + c*t³/3 + d*t⁴/4 - e/t + f = dH`
/// - dS equation: `a*ln(t) + b*t + c*t²/2 + d*t³/3 - e/(2*t²) + g = dS`
///
/// Uses inverse standard deviation weighting to balance the different property scales.
pub fn fit_nist_scaled8(
    temps_K: &[f64],
    cp_orig: &[f64],
    h_orig: &[f64],
    s_orig: &[f64],
) -> Vec<f64> {
    let n = temps_K.len();
    assert!(cp_orig.len() == n && h_orig.len() == n && s_orig.len() == n);

    // Calculate inverse standard deviation weights to normalize different property scales
    // This ensures Cp [~30 J/mol·K], dH [~-100 kJ/mol], and dS [~200 J/mol·K] contribute equally
    let w_cp = 1.0 / std(cp_orig);
    let w_h = 1.0 / std(h_orig);
    let w_s = 1.0 / std(s_orig);

    let rows = 3 * n;
    let cols = 8; // <-- now eight coefficients
    let mut X = DMatrix::<f64>::zeros(rows, cols);
    let mut y = DVector::<f64>::zeros(rows);

    for (i, &T) in temps_K.iter().enumerate() {
        let t = T / 1000.0;
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;

        let row_cp = 3 * i;
        let row_h = 3 * i + 1;
        let row_s = 3 * i + 2;

        //--------------------------------------------------------
        // Cp row: Cp = h + a + b t + c t^2 + d t^3 + e t^(-2)
        //--------------------------------------------------------
        X[(row_cp, 0)] = 1.0; // a
        X[(row_cp, 1)] = t; // b
        X[(row_cp, 2)] = t2; // c
        X[(row_cp, 3)] = t3; // d
        X[(row_cp, 4)] = 1.0 / t2; // e term (t^-2)
        X[(row_cp, 5)] = 0.0; // f
        X[(row_cp, 6)] = 0.0; // g
        X[(row_cp, 7)] = 1.0; // h (Cp-only constant)
        y[row_cp] = cp_orig[i];

        //--------------------------------------------------------
        // H row
        //--------------------------------------------------------
        X[(row_h, 0)] = t;
        X[(row_h, 1)] = t2 / 2.0;
        X[(row_h, 2)] = t3 / 3.0;
        X[(row_h, 3)] = t4 / 4.0;
        X[(row_h, 4)] = -1.0 / t;
        X[(row_h, 5)] = 1.0; // f constant
        X[(row_h, 6)] = 0.0; // g
        X[(row_h, 7)] = 0.0; // h
        y[row_h] = h_orig[i];

        //--------------------------------------------------------
        // S row
        //--------------------------------------------------------
        X[(row_s, 0)] = t.ln();
        X[(row_s, 1)] = t;
        X[(row_s, 2)] = t2 / 2.0;
        X[(row_s, 3)] = t3 / 3.0;
        X[(row_s, 4)] = -0.5 / t2;
        X[(row_s, 5)] = 0.0; // f
        X[(row_s, 6)] = 1.0; // g constant
        X[(row_s, 7)] = 0.0; // h
        y[row_s] = s_orig[i];
    }

    //--- apply weights ---
    for i in 0..n {
        scale_row(&mut X, &mut y, 3 * i, w_cp);
        scale_row(&mut X, &mut y, 3 * i + 1, w_h);
        scale_row(&mut X, &mut y, 3 * i + 2, w_s);
    }

    //--- SVD LS solve ---
    let svd = X.svd(true, true);
    let theta = svd.solve(&y, 1e-12).expect("svd solve failed");

    theta.data.as_vec().clone()
}

//------------------------- Helper Functions -------------------------------

/// Calculate standard deviation for weighting purposes
fn std(v: &[f64]) -> f64 {
    let m = v.iter().copied().sum::<f64>() / v.len() as f64;
    let var = v.iter().map(|&x| (x - m).powi(2)).sum::<f64>() / v.len() as f64;
    var.sqrt()
}

/// Apply weight to a specific row in the least squares system
fn scale_row(X: &mut DMatrix<f64>, y: &mut DVector<f64>, row: usize, w: f64) {
    let cols = X.ncols();
    for c in 0..cols {
        X[(row, c)] *= w;
    }
    y[row] *= w;
}
impl NISTdata {
    /// Prepare thermodynamic data from two adjacent temperature ranges for fitting.
    ///
    /// This function extracts Cp, dH, and dS values from two coefficient sets that
    /// share a common boundary temperature, then combines them into unified datasets.
    ///
    /// # Arguments
    /// * `coeffs_min` - Lower temperature range coefficients
    /// * `coeffs_max` - Higher temperature range coefficients  
    /// * `T_min` - Minimum temperature for fitting range [K]
    /// * `T_max` - Maximum temperature for fitting range [K]
    ///
    /// # Returns
    /// Tuple of (temperatures, cp_values, h_values, s_values)
    pub fn prepare_fit_data_for_2_range_fitting(
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), NISTError> {
        let (_, T_center) = (coeffs_min.T.0, coeffs_min.T.1);
        let (T_center_check, _) = (coeffs_max.T.0, coeffs_max.T.1);
        if T_center != T_center_check {
            return Err(NISTError::InvalidTemperatureRange);
        }

        let a = coeffs_min.coeff.0;
        let b = coeffs_min.coeff.1;
        let c = coeffs_min.coeff.2;
        let d = coeffs_min.coeff.3;
        let e = coeffs_min.coeff.4;
        let f = coeffs_min.coeff.5;
        let g = coeffs_min.coeff.6;
        let h = coeffs_min.coeff.7;

        let mut temps = (T_min as usize..=T_center as usize)
            .step_by(5 as usize)
            .map(|t| t as f64)
            .collect::<Vec<_>>();
        let mut cp_vals = temps
            .iter()
            .map(|&T| calculate_cp(T / 1000.0, a, b, c, d, e))
            .collect::<Vec<_>>();
        let mut h_vals = temps
            .iter()
            .map(|&T| calculate_dh(T / 1000.0, a, b, c, d, e, f, g, h))
            .collect::<Vec<_>>();
        let mut s_vals = temps
            .iter()
            .map(|&T| calculate_s(T / 1000.0, a, b, c, d, e, f, g, h))
            .collect::<Vec<_>>();
        let a = coeffs_max.coeff.0;
        let b = coeffs_max.coeff.1;
        let c = coeffs_max.coeff.2;
        let d = coeffs_max.coeff.3;
        let e = coeffs_max.coeff.4;
        let f = coeffs_max.coeff.5;
        let g = coeffs_max.coeff.6;
        let h = coeffs_max.coeff.7;

        let temps2 = (T_center as usize + 5..=T_max as usize)
            .step_by(5 as usize)
            .map(|t| t as f64)
            .collect::<Vec<_>>();

        temps.extend(temps2.clone());
        let cp_vals2 = temps2
            .clone()
            .iter()
            .map(|&T| calculate_cp(T / 1000.0, a, b, c, d, e))
            .collect::<Vec<_>>();
        let h_vals2 = temps2
            .clone()
            .iter()
            .map(|&T| calculate_dh(T / 1000.0, a, b, c, d, e, f, g, h))
            .collect::<Vec<_>>();
        let s_vals2 = temps2
            .iter()
            .map(|&T| calculate_s(T / 1000.0, a, b, c, d, e, f, g, h))
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
    pub fn fit_cp_dh_ds(
        &mut self,
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<(), NISTError> {
        let mut report_map: HashMap<String, FittingReport> = HashMap::new();
        let (temps, cp_vals, h_vals, s_vals) =
            Self::prepare_fit_data_for_2_range_fitting(coeffs_min, coeffs_max, T_min, T_max)?;

        let coeffs = fit_nist_scaled8(&temps, &cp_vals, &h_vals, &s_vals);
        println!("a..g = {:?}", coeffs);
        assert_eq!(coeffs.len(), 8);
        let a = coeffs[0];
        let b = coeffs[1];
        let c = coeffs[2];
        let d = coeffs[3];
        let e = coeffs[4];
        let f = coeffs[5];
        let g = coeffs[6];
        let h = coeffs[7];
        let Cp_sym = calculate_cp_sym(a, b, c, d, e);
        let dh_sym = calculate_dh_sym(a, b, c, d, e, f, g, h) / Expr::Const(1000.0);
        let ds_sym = calculate_s_sym(a, b, c, d, e, f, g, h);
        let c_fitting_report = self.direct_compare2(&cp_vals, &temps, Cp_sym);
        report_map.insert("Cp fitting report".to_string(), c_fitting_report);
        let dh_fitting_report = self.direct_compare2(&h_vals, &temps, dh_sym);
        report_map.insert("dh fitting report".to_string(), dh_fitting_report);
        let ds_fitting_report = self.direct_compare2(&s_vals, &temps, ds_sym);
        report_map.insert("ds fitting report".to_string(), ds_fitting_report);

        let coeffs = (a, b, c, d, e, f, g, h);
        self.coeffs = Some(coeffs);
        Self::print_quality_table(&report_map);
        Ok(())
    }

    pub fn fitting_cp_dh_ds_non_adjacent(
        &mut self,
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NISTError> {
        let mut report_map: HashMap<String, FittingReport> = HashMap::new();
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
            let h = coeffs.coeff.7;

            let temps = (Tleft as usize + 5..=Tright as usize)
                .step_by(5 as usize)
                .map(|t| t as f64)
                .collect::<Vec<_>>();
            let cp_vals = temps
                .iter()
                .map(|&T| calculate_cp(T / 1000.0, a, b, c, d, e))
                .collect::<Vec<_>>();
            let h_vals = temps
                .iter()
                .map(|&T| calculate_dh(T / 1000.0, a, b, c, d, e, f, g, h))
                .collect::<Vec<_>>();
            let s_vals = temps
                .iter()
                .map(|&T| calculate_s(T / 1000.0, a, b, c, d, e, f, g, h))
                .collect::<Vec<_>>();
            Tvec.extend(temps);
            cp_vec.extend(cp_vals);
            h_vec.extend(h_vals);
            s_vec.extend(s_vals);
        }
        assert_eq!(Tvec.len(), cp_vec.len());
        assert_eq!(Tvec.len(), h_vec.len());
        assert_eq!(Tvec.len(), s_vec.len());

        let coeffs = fit_nist_scaled8(&Tvec, &cp_vec, &h_vec, &s_vec);
        println!("a..g = {:?}", coeffs);
        assert_eq!(coeffs.len(), 8);
        let a = coeffs[0];
        let b = coeffs[1];
        let c = coeffs[2];
        let d = coeffs[3];
        let e = coeffs[4];
        let f = coeffs[5];
        let g = coeffs[6];
        let h = coeffs[7];
        let Cp_sym = calculate_cp_sym(a, b, c, d, e);
        let dh_sym = calculate_dh_sym(a, b, c, d, e, f, g, h) / Expr::Const(1000.0);
        let ds_sym = calculate_s_sym(a, b, c, d, e, f, g, h);
        let c_fitting_report = self.direct_compare2(&cp_vec, &Tvec, Cp_sym);
        report_map.insert("Cp fitting report".to_string(), c_fitting_report);
        let dh_fitting_report = self.direct_compare2(&h_vec, &Tvec, dh_sym);
        report_map.insert("dh fitting report".to_string(), dh_fitting_report);
        let ds_fitting_report = self.direct_compare2(&s_vec, &Tvec, ds_sym);
        report_map.insert("ds fitting report".to_string(), ds_fitting_report);

        let coeffs = (a, b, c, d, e, f, g, h);
        self.coeffs = Some(coeffs);
        Self::print_quality_table(&report_map);

        Ok(report_map)
    }

    //////////////////////////////////////////////////////////////////////////////////////////////

    /// Sequential fitting method for adjacent temperature ranges.
    ///
    /// This is the recommended approach that avoids conflicting constraints between
    /// dH and dS by using a two-stage sequential optimization:
    ///
    /// **Stage 1**: Fit Cp + weighted dS to determine shared coefficients (a,b,c,d,e,g)
    /// **Stage 2**: Fit dH separately using fixed shared coefficients to find (f,h)
    ///
    /// # Arguments
    /// * `coeffs_min` - Lower temperature range coefficients
    /// * `coeffs_max` - Higher temperature range coefficients
    /// * `T_min` - Minimum fitting temperature [K]
    /// * `T_max` - Maximum fitting temperature [K]
    ///
    /// # Mathematical Approach
    /// Uses symbolic expressions and weighted optimization to respect the
    /// thermodynamic relationships while avoiding optimization conflicts.
    pub fn fitting_adjacent_weighted(
        &mut self,
        coeffs_min: Coeffs,
        coeffs_max: Coeffs,
        T_min: f64,
        T_max: f64,
    ) -> Result<(), NISTError> {
        use RustedSciThe::numerical::optimization::fitting_features::SewTwoFunctions;
        use RustedSciThe::numerical::optimization::sym_fitting::Fitting;

        let mut report_map: HashMap<String, FittingReport> = HashMap::new();
        let (temps, cp_vals, h_vals, s_vals) = Self::prepare_fit_data_for_2_range_fitting(
            coeffs_min.clone(),
            coeffs_max.clone(),
            T_min,
            T_max,
        )?;
        //  let w_cp = 1.0 / std(&cp_vals);
        // let w_h = 1.0 / std(&h_vals);
        //  let w_s = 1.0 / std(&s_vals);

        // Use fixed weights instead of std-based weights for better control
        // w_cp = 1.0: Heat capacity baseline weight
        // w_s = 10.0: Higher weight for entropy to improve fitting quality
        let w_cp = 1.0;
        let w_s = 10.0;

        // Define symbolic expressions for NIST thermodynamic equations
        // These represent the mathematical relationships between coefficients and properties
        let Cp_sym = Expr::parse_expression("(a + b*t + c*t^2 + d*t^3 + e/t^2)");
        let dh_sym =
            Expr::parse_expression("(a*t + (b*t^2)/2 + (c*t^3)/3 + (d*t^4)/4 - e/t + f - h)");
        let ds_sym =
            Expr::parse_expression("(a*ln( t ) + b*t + (c*t^2)/2 + (d*t^3)/3 - e/(2*t^2) + g)");

        // Create weighted combination of Cp and dS for Stage 1 fitting
        let combined_func_to_fit = Expr::Const(w_cp) * Cp_sym + Expr::Const(w_s) * ds_sym;
        // Substitute t = T/1000 to convert from Kelvin to the NIST standard form
        let subst = Expr::Var("T".to_string()) / Expr::Const(1000.0);
        let combined_func_to_fit = combined_func_to_fit.substitute_variable("t", &subst);
        let y_data =
            w_cp * DVector::from_vec(cp_vals.clone()) + w_s * DVector::from_vec(s_vals.clone());

        let mut sew = SewTwoFunctions {
            f1: Expr::Const(0.0),
            f2: Expr::Const(0.0),
            x_left: 0.0,
            x_central: 0.0,
            n_points: 0,
            x_right: 0.0,
            fitting_data: (temps.clone(), y_data.data.as_vec().clone()),
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
            combined_func_to_fit.clone(),
            Some(vec![
                "a".to_string(),
                "b".to_string(),
                "c".to_string(),
                "d".to_string(),
                "e".to_string(),
                "g".to_string(),
            ]),
            "T".to_string(),
            initial_guess,
            Some(1e-6),
            Some(1e-6),
            Some(1e-6),
            None,
            Some(50),
        );

        let r_squared = sew.get_r_ssquared().unwrap();
        println!("r_squared for combined fitting: {}", r_squared);
        assert!(
            1.0 - r_squared < 5e-2,
            "Combined fitting failed with R² = {}",
            r_squared
        );
        let coeffs = sew.get_result().unwrap();
        assert_eq!(coeffs.len(), 6);
        let a = coeffs[0];
        let b = coeffs[1];
        let c = coeffs[2];
        let d = coeffs[3];
        let e = coeffs[4];

        let g = coeffs[5];

        let mut map_of_coeffs = HashMap::from([
            ("a".to_string(), a),
            ("b".to_string(), b),
            ("c".to_string(), c),
            ("d".to_string(), d),
            ("e".to_string(), e),
            ("g".to_string(), g),
        ]);

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
        let initial_guess = vec![
            (coeffs_min.coeff.5 + coeffs_max.coeff.5) / 2.0,
            (coeffs_min.coeff.7 + coeffs_max.coeff.7) / 2.0,
        ];
        let dh_sym = dh_sym
            .substitute_variable("t", &subst)
            .set_variable_from_map(&map_of_coeffs);

        sew2.fit(
            dh_sym.clone(),
            Some(vec!["f".to_string(), "h".to_string()]),
            "T".to_string(),
            initial_guess,
            Some(1e-6),
            Some(1e-6),
            Some(1e-6),
            None,
            Some(50),
        );
        let coeffs = sew2.get_result().unwrap();
        assert_eq!(coeffs.len(), 2);
        let f = coeffs[0];
        let h = coeffs[1];
        map_of_coeffs.insert("f".to_string(), f);
        map_of_coeffs.insert("h".to_string(), h);
        let Cp_sym = calculate_cp_sym(a, b, c, d, e);
        let dh_sym = calculate_dh_sym(a, b, c, d, e, f, g, h) / Expr::Const(1000.0);
        let ds_sym = calculate_s_sym(a, b, c, d, e, f, g, h);
        let combined_func_to_fit = combined_func_to_fit.set_variable_from_map(&map_of_coeffs);
        let comb_func_fitting_report =
            self.direct_compare2(&y_data.data.as_vec().clone(), &temps, combined_func_to_fit);
        report_map.insert("weighted function".to_string(), comb_func_fitting_report);
        let c_fitting_report = self.direct_compare2(&cp_vals, &temps, Cp_sym);
        report_map.insert("Cp fitting report".to_string(), c_fitting_report);
        let dh_fitting_report = self.direct_compare2(&h_vals, &temps, dh_sym);
        report_map.insert("dh fitting report".to_string(), dh_fitting_report);
        let ds_fitting_report = self.direct_compare2(&s_vals, &temps, ds_sym);
        report_map.insert("ds fitting report".to_string(), ds_fitting_report);

        self.coeffs = Some((a, b, c, d, e, f, g, h));

        Self::print_quality_table(&report_map);
        Ok(())
    }

    /// Non-adjacent fitting method for multiple temperature intervals
    pub fn fitting_non_adjacent(
        &mut self,
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NISTError> {
        use RustedSciThe::numerical::optimization::fitting_features::SewTwoFunctions;
        use RustedSciThe::numerical::optimization::sym_fitting::Fitting;

        let mut report_map: HashMap<String, FittingReport> = HashMap::new();
        let (temps, cp_vals, h_vals, s_vals) =
            Self::prepare_fit_data_for_non_adjacent_fitting(coeffs.clone(), T_min, T_max)?;

        let w_cp = 1.0;
        let w_s = 10.0;
        // Step 1: Fit Cp + dS (they work well together)
        let Cp_sym = Expr::parse_expression("(a + b*t + c*t^2 + d*t^3 + e/t^2)");
        let ds_sym =
            Expr::parse_expression("(a*ln( t ) + b*t + (c*t^2)/2 + (d*t^3)/3 - e/(2*t^2) + g)");
        let combined_func_to_fit = Expr::Const(w_cp) * Cp_sym + Expr::Const(w_s) * ds_sym;
        let subst = Expr::Var("T".to_string()) / Expr::Const(1000.0);
        let combined_func_to_fit = combined_func_to_fit.substitute_variable("t", &subst);

        let y_data =
            w_cp * DVector::from_vec(cp_vals.clone()) + w_s * DVector::from_vec(s_vals.clone());

        let mut sew = SewTwoFunctions {
            f1: Expr::Const(0.0),
            f2: Expr::Const(0.0),
            x_left: 0.0,
            x_central: 0.0,
            n_points: 0,
            x_right: 0.0,
            fitting_data: (temps.clone(), y_data.data.as_vec().clone()),
            fitting: Fitting::new(),
        };

        // Average initial guess from all coefficient sets
        let mut initial_guess = vec![0.0; 6];
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
            combined_func_to_fit,
            Some(vec![
                "a".to_string(),
                "b".to_string(),
                "c".to_string(),
                "d".to_string(),
                "e".to_string(),
                "g".to_string(),
            ]),
            "T".to_string(),
            initial_guess,
            None,
            None,
            None,
            None,
            None,
        );

        let coeffs_step1 = sew.get_result().unwrap();
        let (a, b, c, d, e, g) = (
            coeffs_step1[0],
            coeffs_step1[1],
            coeffs_step1[2],
            coeffs_step1[3],
            coeffs_step1[4],
            coeffs_step1[5],
        );

        // Step 2: Fit f and h for dH using fixed a,b,c,d,e
        let dh_sym =
            Expr::parse_expression("(a*t + (b*t^2)/2 + (c*t^3)/3 + (d*t^4)/4 - e/t + f - h)");
        let dh_func_fixed = dh_sym
            .substitute_variable("a", &Expr::Const(a))
            .substitute_variable("b", &Expr::Const(b))
            .substitute_variable("c", &Expr::Const(c))
            .substitute_variable("d", &Expr::Const(d))
            .substitute_variable("e", &Expr::Const(e))
            .substitute_variable("t", &subst);

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

        let mut initial_guess_fh = vec![0.0; 2];
        for coeff in coeffs.iter() {
            initial_guess_fh[0] += coeff.coeff.5;
            initial_guess_fh[1] += coeff.coeff.7;
        }
        initial_guess_fh = initial_guess_fh
            .iter()
            .map(|x| x / coeffs.len() as f64)
            .collect();

        sew2.fit(
            dh_func_fixed,
            Some(vec!["f".to_string(), "h".to_string()]),
            "T".to_string(),
            initial_guess_fh,
            None,
            None,
            None,
            None,
            None,
        );

        let coeffs_step2 = sew2.get_result().unwrap();
        let (f, h) = (coeffs_step2[0], coeffs_step2[1]);

        // Final coefficients and reports
        let Cp_sym = calculate_cp_sym(a, b, c, d, e);
        let dh_sym = calculate_dh_sym(a, b, c, d, e, f, g, h) / Expr::Const(1000.0);
        let ds_sym = calculate_s_sym(a, b, c, d, e, f, g, h);

        let c_fitting_report = self.direct_compare2(&cp_vals, &temps, Cp_sym);
        report_map.insert("Cp fitting report".to_string(), c_fitting_report);
        let dh_fitting_report = self.direct_compare2(&h_vals, &temps, dh_sym);
        report_map.insert("dh fitting report".to_string(), dh_fitting_report);
        let ds_fitting_report = self.direct_compare2(&s_vals, &temps, ds_sym);
        report_map.insert("ds fitting report".to_string(), ds_fitting_report);

        self.coeffs = Some((a, b, c, d, e, f, g, h));
        Self::print_quality_table(&report_map);
        Ok(report_map)
    }

    fn prepare_fit_data_for_non_adjacent_fitting(
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), NISTError> {
        let mut Tvec: Vec<f64> = Vec::new();
        let mut cp_vec: Vec<f64> = Vec::new();
        let mut h_vec: Vec<f64> = Vec::new();
        let mut s_vec: Vec<f64> = Vec::new();

        for coeff in coeffs.iter() {
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

            let (a, b, c, d, e, f, g, h) = coeff.coeff;
            let temps = (Tleft as usize + 5..=Tright as usize)
                .step_by(5)
                .map(|t| t as f64)
                .collect::<Vec<_>>();
            let cp_vals = temps
                .iter()
                .map(|&T| calculate_cp(T / 1000.0, a, b, c, d, e))
                .collect::<Vec<_>>();
            let h_vals = temps
                .iter()
                .map(|&T| calculate_dh(T / 1000.0, a, b, c, d, e, f, g, h))
                .collect::<Vec<_>>();
            let s_vals = temps
                .iter()
                .map(|&T| calculate_s(T / 1000.0, a, b, c, d, e, f, g, h))
                .collect::<Vec<_>>();

            Tvec.extend(temps);
            cp_vec.extend(cp_vals);
            h_vec.extend(h_vals);
            s_vec.extend(s_vals);
        }

        Ok((Tvec, cp_vec, h_vec, s_vec))
    }

    //////////////////////////////////MISC/////////////////////////////////////
    /// Print quality table for fitting reports
    pub fn print_quality_table(report_map: &HashMap<String, FittingReport>) {
        let mut table = Table::new();
        table.add_row(Row::new(vec![
            Cell::new("Report Type"),
            Cell::new("L2 Norm"),
            Cell::new("Max Norm"),
            Cell::new("Significant Points"),
            Cell::new("T Range"),
        ]));

        for (report_name, report) in report_map {
            let t_range_str = match report.T_range {
                Some((min_t, max_t)) => format!("{:.1}-{:.1}", min_t, max_t),
                None => "None".to_string(),
            };
            table.add_row(Row::new(vec![
                Cell::new(report_name),
                Cell::new(&format!("{:.6}", report.l2_norm)),
                Cell::new(&format!("{:.6}", report.max_norm)),
                Cell::new(&report.number_of_significant_points.to_string()),
                Cell::new(&t_range_str),
            ]));
        }
        table.printstd();
    }

    pub fn direct_compare2(
        &self,
        y_data: &[f64],
        T_vec: &[f64],
        fit_function: Expr,
    ) -> FittingReport {
        let n = y_data.len();

        let fitted_data_func = fit_function.lambdify1D();
        let y_data_fitted = T_vec
            .iter()
            .map(|&T| fitted_data_func(T))
            .collect::<Vec<f64>>();
        assert_eq!(y_data.len(), y_data_fitted.len());
        let mut l1_norm = 0.0;
        let mut l2_norm = 0.0;
        let mut max_norm: f64 = 0.0;
        let mut map_of_residuals: HashMap<usize, (f64, f64)> = HashMap::new();

        for (i, (y, y_fit)) in y_data.iter().zip(y_data_fitted.iter()).enumerate() {
            let diff = (y - y_fit).abs();
            l1_norm += diff;
            l2_norm += diff * diff;
            let rel_diff = diff / y;
            if rel_diff.abs() > 10e-2 {
                map_of_residuals.insert(i, (T_vec[i], rel_diff));
            }
            max_norm = max_norm.max((diff / y).abs());
        }

        l1_norm = l1_norm / (2.0 * n as f64);
        l2_norm = l2_norm.sqrt() / (2.0 * n as f64);

        let mut T_range: Option<(f64, f64)> = None;
        if !map_of_residuals.is_empty() {
            let max_T = *map_of_residuals
                .iter()
                .map(|(_, (T, _))| T)
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap();
            let min_T = *map_of_residuals
                .iter()
                .map(|(_, (T, _))| T)
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap();
            T_range = Some((min_T, max_T));
        }

        FittingReport {
            l2_norm,
            max_norm,
            number_of_significant_points: map_of_residuals.len(),
            T_range,
        }
    }
}

///////////////////////////TESTS//////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::DBhandlers::NIST_parser::{Phase, SearchType};
    use approx::assert_relative_eq;
    use std::{thread, time};
    #[allow(non_upper_case_globals)]
    const smtime: time::Duration = time::Duration::from_secs(5);
    use thread::sleep;

    #[test]
    fn test_fit_cp_dh_ds_adjacent() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO2".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if nist.coeffs_map.len() >= 2 {
            let coeffs_min = nist.coeffs_map.get(&0).unwrap().clone();
            let coeffs_max = nist.coeffs_map.get(&1).unwrap().clone();

            if coeffs_min.T.1 == coeffs_max.T.0 {
                let T_center = coeffs_min.T.1;
                let T_min = T_center - 100.0;
                let T_max = T_center + 100.0;
                nist.set_T_interval(T_min, T_max);

                let result =
                    nist.fit_cp_dh_ds(coeffs_min.clone(), coeffs_max.clone(), T_min, T_max);
                assert!(result.is_ok());
                assert!(nist.coeffs.is_some());

                let _ = nist.calculate_cp_dh_ds(T_center);
                assert!(nist.Cp > 0.0);
                assert!(nist.dh != 0.0);
                assert!(nist.ds != 0.0);
            }
        }
        sleep(smtime);
    }

    #[test]
    fn test_fit_cp_dh_ds_consistency() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if nist.coeffs_map.len() >= 2 {
            let coeffs_min = nist.coeffs_map.get(&0).unwrap().clone();
            let coeffs_max = nist.coeffs_map.get(&1).unwrap().clone();

            if coeffs_min.T.1 == coeffs_max.T.0 {
                let T_center = coeffs_min.T.1;
                let T_min = T_center - 50.0;
                let T_max = T_center + 50.0;
                nist.set_T_interval(T_min, T_max);
                let result = nist.fit_cp_dh_ds(coeffs_min, coeffs_max, T_min, T_max);
                assert!(result.is_ok());

                let _ = nist.create_closure_cp_dh_ds();
                let _ = nist.calculate_cp_dh_ds(T_center);

                let cp_calc = nist.Cp;
                let dh_calc = nist.dh;
                let ds_calc = nist.ds;

                let cp_closure = (nist.C_fun)(T_center);
                let dh_closure = (nist.dh_fun)(T_center);
                let ds_closure = (nist.ds_fun)(T_center);

                assert_relative_eq!(cp_closure, cp_calc, epsilon = 1e-6);
                assert_relative_eq!(dh_closure, dh_calc, epsilon = 1e-6);
                assert_relative_eq!(ds_closure, ds_calc, epsilon = 1e-6);
            }
        }
        sleep(smtime);
    }

    #[test]
    fn test_fit_cp_dh_ds_error_handling() {
        let mut nist = NISTdata::new();

        let coeffs_min = Coeffs {
            T: (300.0, 600.0),
            coeff: (33.0, -5.0, 1.0, -0.1, 0.01, -1000.0, 200.0, 0.0),
        };
        let coeffs_max = Coeffs {
            T: (700.0, 1000.0),
            coeff: (35.0, -4.0, 0.8, -0.08, 0.008, -1200.0, 220.0, 0.0),
        };

        let result = nist.fit_cp_dh_ds(coeffs_min, coeffs_max, 400.0, 800.0);
        assert!(result.is_err());
        sleep(smtime);
    }

    #[test]
    fn test_fit_cp_dh_ds_identical_coeffs() {
        let mut nist = NISTdata::new();
        let identical_coeffs = (33.0, -5.0, 1.0, -0.1, 0.01, -1000.0, 200.0, 0.0);

        let coeffs_min = Coeffs {
            T: (300.0, 600.0),
            coeff: identical_coeffs,
        };
        let coeffs_max = Coeffs {
            T: (600.0, 900.0),
            coeff: identical_coeffs,
        };

        let result = nist.fit_cp_dh_ds(coeffs_min, coeffs_max, 400.0, 800.0);
        assert!(result.is_ok());
        assert!(nist.coeffs.is_some());

        let fitted_coeffs = nist.coeffs.unwrap();
        assert_relative_eq!(fitted_coeffs.0, identical_coeffs.0, epsilon = 1e-2);
        assert_relative_eq!(fitted_coeffs.1, identical_coeffs.1, epsilon = 1e-2);
        sleep(smtime);
    }
    ////////////////////////////NON ADJACENT//////////////////////////////////////
    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_basic() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO2".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if nist.coeffs_map.len() >= 3 {
            let coeffs_vec = vec![
                nist.coeffs_map.get(&0).unwrap().clone(),
                nist.coeffs_map.get(&1).unwrap().clone(),
                nist.coeffs_map.get(&2).unwrap().clone(),
            ];

            let T_min = 400.0;
            let T_max = 1200.0;
            nist.set_T_interval(T_min, T_max);
            let result = nist.fitting_cp_dh_ds_non_adjacent(coeffs_vec, T_min, T_max);
            assert!(result.is_ok());
            assert!(nist.coeffs.is_some());

            let reports = result.unwrap();
            assert!(reports.contains_key("Cp fitting report"));
            assert!(reports.contains_key("dh fitting report"));
            assert!(reports.contains_key("ds fitting report"));
        }
        sleep(smtime);
    }

    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_single_range() {
        let mut nist = NISTdata::new();
        let coeffs_vec = vec![Coeffs {
            T: (300.0, 800.0),
            coeff: (33.0, -5.0, 1.0, -0.1, 0.01, -1000.0, 200.0, 0.0),
        }];
        nist.set_T_interval(300.0, 800.0);
        let result = nist.fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 700.0);
        assert!(result.is_ok());
        assert!(nist.coeffs.is_some());
        sleep(smtime);
    }

    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_overlapping_ranges() {
        let mut nist = NISTdata::new();
        let coeffs_vec = vec![
            Coeffs {
                T: (300.0, 600.0),
                coeff: (33.0, -5.0, 1.0, -0.1, 0.01, -1000.0, 200.0, 0.0),
            },
            Coeffs {
                T: (500.0, 900.0),
                coeff: (35.0, -4.0, 0.8, -0.08, 0.008, -1200.0, 220.0, 0.0),
            },
        ];
        nist.set_T_interval(400.0, 800.0);
        let result = nist.fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 800.0);
        assert!(result.is_ok());
        assert!(nist.coeffs.is_some());
        sleep(smtime);
    }

    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_partial_overlap() {
        let mut nist = NISTdata::new();
        let coeffs_vec = vec![
            Coeffs {
                T: (200.0, 500.0),
                coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
            },
            Coeffs {
                T: (700.0, 1200.0),
                coeff: (40.0, -6.0, 1.5, -0.15, 0.015, -1500.0, 250.0, 0.0),
            },
        ];

        let result = nist.fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 1000.0);
        assert!(result.is_ok());
        assert!(nist.coeffs.is_some());
        sleep(smtime);
    }

    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_uses_both_intervals() {
        let mut nist1 = NISTdata::new();
        let mut nist2 = NISTdata::new();

        let coeffs1 = Coeffs {
            T: (300.0, 600.0),
            coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
        };
        let coeffs2 = Coeffs {
            T: (700.0, 1000.0),
            coeff: (40.0, -6.0, 1.5, -0.15, 0.015, -1500.0, 250.0, 0.0),
        };

        let result1 = nist1.fitting_cp_dh_ds_non_adjacent(vec![coeffs1.clone()], 400.0, 900.0);
        assert!(result1.is_ok());
        let coeffs_single = nist1.coeffs.unwrap();

        let result2 = nist2.fitting_cp_dh_ds_non_adjacent(vec![coeffs1, coeffs2], 400.0, 900.0);
        assert!(result2.is_ok());
        let coeffs_both = nist2.coeffs.unwrap();

        let diff_a = (coeffs_single.0 - coeffs_both.0).abs();
        let diff_b = (coeffs_single.1 - coeffs_both.1).abs();

        assert!(
            diff_a > 1e-6 || diff_b > 1e-6,
            "Coefficients should differ when using both intervals"
        );
        sleep(smtime);
    }

    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_consistency() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("H2O".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if nist.coeffs_map.len() >= 2 {
            let coeffs_vec = vec![
                nist.coeffs_map.get(&0).unwrap().clone(),
                nist.coeffs_map.get(&1).unwrap().clone(),
            ];

            let T_min = 1000.0;
            let T_max = 1900.0;

            let result = nist.fitting_cp_dh_ds_non_adjacent(coeffs_vec, T_min, T_max);
            assert!(result.is_ok());

            let _ = nist.calculate_cp_dh_ds(1750.0);
            assert!(nist.Cp > 0.0);
            assert!(nist.dh != 0.0);
            assert!(nist.ds != 0.0);
        }
        sleep(smtime);
    }
    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_debug_data_points() {
        let mut nist = NISTdata::new();

        let coeffs_vec = vec![
            Coeffs {
                T: (300.0, 500.0),
                coeff: (25.0, -2.0, 0.3, -0.03, 0.003, -500.0, 150.0, 0.0),
            },
            Coeffs {
                T: (600.0, 900.0),
                coeff: (35.0, -4.0, 0.8, -0.08, 0.008, -1000.0, 200.0, 0.0),
            },
        ];

        let result = nist.fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 800.0);
        assert!(result.is_ok());

        let _ = nist.calculate_cp_dh_ds(450.0);
        let cp1 = nist.Cp;
        let _ = nist.calculate_cp_dh_ds(700.0);
        let cp2 = nist.Cp;

        assert!(cp1 > 0.0);
        assert!(cp2 > 0.0);
        assert_ne!(
            cp1, cp2,
            "Cp values should differ at different temperatures"
        );
        sleep(smtime);
    }
    ////////////////////////////////////////////////////////////////////////
    #[test]
    fn test_fit_nist_scaled8_function() {
        let temps = vec![400.0, 500.0, 600.0, 700.0, 800.0];
        let cp_vals = vec![29.1, 30.2, 31.5, 32.8, 34.1];
        let h_vals = vec![-110.5, -107.2, -103.8, -100.3, -96.7];
        let s_vals = vec![197.7, 206.2, 213.8, 220.6, 226.7];

        let coeffs = fit_nist_scaled8(&temps, &cp_vals, &h_vals, &s_vals);
        assert_eq!(coeffs.len(), 8);

        assert!(coeffs[0] > 0.0);
        sleep(smtime);
    }

    #[test]
    fn test_fit_cp_dh_ds_multiple_substances() {
        let substances = vec!["H2O", "NO2", "O2"];

        for substance in substances {
            let mut nist = NISTdata::new();
            let _ = nist.get_data_from_NIST(substance.to_owned(), SearchType::All, Phase::Gas);
            nist.parse_coefficients().unwrap();

            if nist.coeffs_map.len() >= 2 {
                let coeffs_min = nist.coeffs_map.get(&0).unwrap().clone();
                let coeffs_max = nist.coeffs_map.get(&1).unwrap().clone();

                if coeffs_min.T.1 == coeffs_max.T.0 {
                    let T_center = coeffs_min.T.1;
                    let T_min = T_center - 50.0;
                    let T_max = T_center + 50.0;

                    let result = nist.fit_cp_dh_ds(coeffs_min, coeffs_max, T_min, T_max);
                    assert!(result.is_ok(), "Fitting failed for {}", substance);
                    assert!(nist.coeffs.is_some());
                }
            }
        }
        sleep(smtime);
    }
    ///////////////////////////////////////////////////////////////////////////
    ///
    #[test]
    fn test_fitting_adjacent_weighted_identical_coeffs2() {
        use crate::Thermodynamics::DBhandlers::NISTdata::Coeffs;

        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO2".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();
        let identical_coeffs = nist.coeffs_map.get(&0).unwrap().coeff.clone();

        let coeffs_min = Coeffs {
            T: (400.0, 700.0),
            coeff: identical_coeffs,
        };
        let coeffs_max = Coeffs {
            T: (700.0, 1000.0),
            coeff: identical_coeffs,
        };
        nist.set_T_interval(600.0, 800.0);
        let result = nist.fitting_adjacent_weighted(coeffs_min, coeffs_max, 600.0, 800.0);
        assert!(result.is_ok());

        let fitted_coeffs = nist.coeffs.unwrap();
        assert_relative_eq!(fitted_coeffs.0, identical_coeffs.0, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.1, identical_coeffs.1, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.2, identical_coeffs.2, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.3, identical_coeffs.3, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.4, identical_coeffs.4, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.5, identical_coeffs.5, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.6, identical_coeffs.6, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.7, identical_coeffs.7, epsilon = 1e-3);
        sleep(smtime)
    }

    #[test]
    fn test_fitting_adjacent_weighted_identical_coeffs() {
        use crate::Thermodynamics::DBhandlers::NISTdata::Coeffs;

        let mut nist = NISTdata::new();
        let identical_coeffs = (35.0, -4.0, 0.8, -0.08, 0.008, -1200.0, 220.0, 0.0);

        let coeffs_min = Coeffs {
            T: (400.0, 700.0),
            coeff: identical_coeffs,
        };
        let coeffs_max = Coeffs {
            T: (700.0, 1000.0),
            coeff: identical_coeffs,
        };
        nist.set_T_interval(600.0, 800.0);
        let result = nist.fitting_adjacent_weighted(coeffs_min, coeffs_max, 600.0, 800.0);
        assert!(result.is_ok());

        let fitted_coeffs = nist.coeffs.unwrap();
        assert_relative_eq!(fitted_coeffs.0, identical_coeffs.0, epsilon = 1e-1);
        assert_relative_eq!(fitted_coeffs.1, identical_coeffs.1, epsilon = 1e-1);
        assert_relative_eq!(fitted_coeffs.2, identical_coeffs.2, epsilon = 1e-1);
        assert_relative_eq!(fitted_coeffs.3, identical_coeffs.3, epsilon = 1e-1);
        assert_relative_eq!(fitted_coeffs.4, identical_coeffs.4, epsilon = 1e-1);
        assert_relative_eq!(fitted_coeffs.5, identical_coeffs.5, epsilon = 1e-1);
        assert_relative_eq!(fitted_coeffs.6, identical_coeffs.6, epsilon = 1e-1);
        assert_relative_eq!(fitted_coeffs.7, identical_coeffs.7, epsilon = 1e-1);
        sleep(smtime)
    }

    #[test]
    fn test_fitting_adjacent_methods2() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO2".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();
        let _ = nist.extract_coefficients(400.0);
        let _ = nist.create_closure_cp_dh_ds();
        if nist.coeffs_map.len() >= 2 {
            let coeffs_min = nist.coeffs_map.get(&0).unwrap().clone();
            let coeffs_max = nist.coeffs_map.get(&1).unwrap().clone();

            if coeffs_min.T.1 == coeffs_max.T.0 {
                let T_center = coeffs_min.T.1;
                let T_min = T_center - 100.0;
                let T_max = T_center + 200.0;

                nist.set_T_interval(T_min, T_max);
                let result = nist.fitting_adjacent_weighted(coeffs_min, coeffs_max, T_min, T_max);
                assert!(result.is_ok());
                assert!(nist.coeffs.is_some());
                let _ = nist.calculate_cp_dh_ds(T_center);
                println!("Cp: {}", nist.Cp);
                println!("dh {}", nist.dh);
                println!("ds {}", nist.ds);
                let Cp_closure = nist.C_fun;
                let dh_closure = nist.dh_fun;
                let ds_closure = nist.ds_fun;
                assert_relative_eq!(nist.Cp, Cp_closure(T_center), epsilon = 5.0);
                assert_relative_eq!(nist.dh, dh_closure(T_center), epsilon = 500.0);
                assert_relative_eq!(nist.ds, ds_closure(T_center), epsilon = 5.0);
                println!("Cp from closure: {}", Cp_closure(T_center));
                println!("dh from closure: {}", dh_closure(T_center));
                println!("ds from closure: {}", ds_closure(T_center));
            }
        }
        sleep(smtime)
    }

    ////////////////////////////FITTING NON ADJACENT//////////////////////////////////////
    #[test]
    fn test_fitting_non_adjacent_basic() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO2".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();
        let _ = nist.extract_coefficients(400.0);
        let _ = nist.create_closure_cp_dh_ds();
        if nist.coeffs_map.len() >= 2 {
            let coeffs_min = nist.coeffs_map.get(&0).unwrap().clone();
            let coeffs_max = nist.coeffs_map.get(&1).unwrap().clone();

            if coeffs_min.T.1 == coeffs_max.T.0 {
                let T_center = coeffs_min.T.1;
                let T_min = T_center - 100.0;
                let T_max = T_center + 200.0;

                nist.set_T_interval(T_min, T_max);
                let result = nist.fitting_non_adjacent(vec![coeffs_min, coeffs_max], T_min, T_max);
                assert!(result.is_ok());
                assert!(nist.coeffs.is_some());
                let _ = nist.calculate_cp_dh_ds(T_center);
                println!("Cp: {}", nist.Cp);
                println!("dh {}", nist.dh);
                println!("ds {}", nist.ds);
                let Cp_closure = nist.C_fun;
                let dh_closure = nist.dh_fun;
                let ds_closure = nist.ds_fun;
                assert_relative_eq!(nist.Cp, Cp_closure(T_center), epsilon = 5.0);
                assert_relative_eq!(nist.dh, dh_closure(T_center), epsilon = 500.0);
                assert_relative_eq!(nist.ds, ds_closure(T_center), epsilon = 5.0);
                println!("Cp from closure: {}", Cp_closure(T_center));
                println!("dh from closure: {}", dh_closure(T_center));
                println!("ds from closure: {}", ds_closure(T_center));
            }
        }
        sleep(smtime)
    }

    #[test]
    fn test_fitting_non_adjacent_sequential_vs_simultaneous() {
        let mut nist1 = NISTdata::new();
        let mut nist2 = NISTdata::new();

        let coeffs_vec = vec![
            Coeffs {
                T: (300.0, 600.0),
                coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
            },
            Coeffs {
                T: (700.0, 1000.0),
                coeff: (40.0, -6.0, 1.5, -0.15, 0.015, -1500.0, 250.0, 0.0),
            },
        ];

        let result1 = nist1.fitting_non_adjacent(coeffs_vec.clone(), 400.0, 900.0);
        assert!(result1.is_ok());
        let coeffs_sequential = nist1.coeffs.unwrap();

        let result2 = nist2.fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 900.0);
        assert!(result2.is_ok());
        let coeffs_simultaneous = nist2.coeffs.unwrap();

        // Sequential fitting should give different results than simultaneous
        let diff_a = (coeffs_sequential.0 - coeffs_simultaneous.0).abs();
        assert!(
            diff_a > 1e-6,
            "Sequential and simultaneous fitting should differ"
        );
        sleep(smtime);
    }

    #[test]
    fn test_fitting_non_adjacent_consistency() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("H2O".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if nist.coeffs_map.len() >= 2 {
            let coeffs_vec = vec![
                nist.coeffs_map.get(&0).unwrap().clone(),
                nist.coeffs_map.get(&1).unwrap().clone(),
            ];

            let T_min = 500.0;
            let T_max = 1000.0;

            let result = nist.fitting_non_adjacent(coeffs_vec, T_min, T_max);
            assert!(result.is_ok());

            let _ = nist.calculate_cp_dh_ds(750.0);
            assert!(nist.Cp > 0.0);
            assert!(nist.dh != 0.0);
            assert!(nist.ds != 0.0);
        }
        sleep(smtime);
    }

    #[test]
    fn test_prepare_fit_data_for_non_adjacent_fitting_basic() {
        let coeffs_vec = vec![
            Coeffs {
                T: (300.0, 600.0),
                coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
            },
            Coeffs {
                T: (700.0, 1000.0),
                coeff: (40.0, -6.0, 1.5, -0.15, 0.015, -1500.0, 250.0, 0.0),
            },
        ];

        let result = NISTdata::prepare_fit_data_for_non_adjacent_fitting(coeffs_vec, 400.0, 900.0);
        assert!(result.is_ok());

        let (temps, cp_vals, h_vals, s_vals) = result.unwrap();
        assert!(!temps.is_empty());
        assert_eq!(temps.len(), cp_vals.len());
        assert_eq!(temps.len(), h_vals.len());
        assert_eq!(temps.len(), s_vals.len());

        // Should have data from both ranges
        assert!(temps.iter().any(|&t| t >= 400.0 && t <= 600.0));
        assert!(temps.iter().any(|&t| t >= 700.0 && t <= 900.0));
        sleep(smtime);
    }

    #[test]
    fn test_prepare_fit_data_for_non_adjacent_fitting_no_overlap() {
        let coeffs_vec = vec![
            Coeffs {
                T: (100.0, 200.0),
                coeff: (25.0, -2.0, 0.3, -0.03, 0.003, -500.0, 150.0, 0.0),
            },
            Coeffs {
                T: (1500.0, 2000.0),
                coeff: (45.0, -8.0, 2.0, -0.2, 0.02, -2000.0, 300.0, 0.0),
            },
        ];

        let result = NISTdata::prepare_fit_data_for_non_adjacent_fitting(coeffs_vec, 400.0, 800.0);
        assert!(result.is_ok());

        let (temps, cp_vals, h_vals, s_vals) = result.unwrap();
        // Should be empty since no ranges overlap with 400-800K
        assert!(temps.is_empty());
        assert!(cp_vals.is_empty());
        assert!(h_vals.is_empty());
        assert!(s_vals.is_empty());
        sleep(smtime);
    }

    #[test]
    fn test_prepare_fit_data_for_non_adjacent_fitting_partial_overlap() {
        let coeffs_vec = vec![
            Coeffs {
                T: (200.0, 500.0),
                coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
            },
            Coeffs {
                T: (700.0, 1200.0),
                coeff: (40.0, -6.0, 1.5, -0.15, 0.015, -1500.0, 250.0, 0.0),
            },
        ];

        let result = NISTdata::prepare_fit_data_for_non_adjacent_fitting(coeffs_vec, 400.0, 1000.0);
        assert!(result.is_ok());

        let (temps, _, _, _) = result.unwrap();
        assert!(!temps.is_empty());

        // Should have data from 400-500K and 700-1000K
        let min_temp = temps
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap()
            .clone();

        let max_temp = temps
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap()
            .clone();

        assert!(min_temp >= 400.0);
        assert!(max_temp <= 1000.0);
        sleep(smtime);
    }
}
