//! # NIST Thermodynamic Data Fitting Module
//!
//! This module fits NIST thermodynamic polynomial coefficients while preserving
//! the linked heat-capacity (Cp), enthalpy (dH), and entropy (dS) model.
//!
//! ## Two Main Fitting Algorithms
//!
//! ### Production stitching
//! Adjacent and non-adjacent stitching use one deterministic linear SVD fit over
//! samples from all selected intervals. The model has seven identifiable
//! coordinates: the eighth published coefficient is a preserved enthalpy gauge.
//!
//! ### Compatibility fitting
//! The older symbolic sequential optimizers remain available through explicitly
//! named compatibility methods; production NIST stitching does not depend on
//! nonlinear fitting for this coefficient-linear model.
//!
//! ## NIST Polynomial Forms
//! The fitting targets these standard NIST equations:
//! - **Cp(T)**: `a + b*t + c*t^2 + d*t^3 + e*t^-2` (where t = T/1000)
//! - **dH(T)**: `a*t + b*t^2/2 + c*t^3/3 + d*t^4/4 - e/t + f - h`
//! - **dS(T)**: `a*ln(t) + b*t + c*t^2/2 + d*t^3/3 - e/(2*t^2) + g`
//! ## Key Features
//!
//! - Shared coefficient equations for all three properties
//! - Deterministic sampling, interval ordering, and join-boundary ownership
//! - Property-level validation before fitted coefficients are published

use crate::Thermodynamics::DBhandlers::NISTdata::*;
use crate::Thermodynamics::DBhandlers::fitting_reports::{
    compare_property_series, compare_property_series_with_absolute_relative_tolerance,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};
use prettytable::{Cell, Row, Table};
use std::collections::HashMap;
const MAX_PROPERTY_RELATIVE_ERROR: f64 = 2.0e-2;
const NIST_CP_ABSOLUTE_TOLERANCE: f64 = 0.05;
const NIST_H_ABSOLUTE_TOLERANCE_KJ_PER_MOL: f64 = 0.05;
const NIST_S_ABSOLUTE_TOLERANCE: f64 = 0.05;
const MIN_LM_R_SQUARED: f64 = 0.95;
const FIT_SAMPLE_STEP_K: f64 = 5.0;
const MAX_FIT_SAMPLES: usize = 250_000;

/// Compute global NIST-form coefficients using SVD least squares fitting.
///
/// This function simultaneously fits all three thermodynamic properties (Cp, dH, dS)
/// using a weighted least squares approach with SVD decomposition.
///
/// # Arguments
/// * `temps_K` - Temperature sampling points in Kelvin
/// * `cp_orig` - Heat capacity values [J/molВ·K]
/// * `h_orig` - Enthalpy values [kJ/mol]
/// * `s_orig` - Entropy values [J/molВ·K]
///
/// # Returns
/// Vector of 8 NIST coefficients `[A,B,C,D,E,F,G,H]`.
///
/// Property samples identify only `F-H`. This low-level fit uses the canonical
/// gauge `H = 0` and returns `F = F-H`. A metadata-owning caller may select
/// another `H` and add it to `F` without changing any fitted property.
///
/// # Mathematical Approach
/// - Cp equation: `a + b*t + c*t^2 + d*t^3 + e*t^-2 = Cp`
/// - dH equation: `a*t + b*t^2/2 + c*t^3/3 + d*t^4/4 - e/t + (f-h) = dH`
/// - dS equation: `a*ln(t) + b*t + c*t^2/2 + d*t^3/3 - e/(2*t^2) + g = dS`
/// Uses inverse RMS weighting to balance the property scales without dividing
/// by zero for a constant property series.
pub fn fit_nist_scaled8(
    temps_K: &[f64],
    cp_orig: &[f64],
    h_orig: &[f64],
    s_orig: &[f64],
) -> Result<Vec<f64>, NISTError> {
    let n = temps_K.len();
    validate_fit_inputs(temps_K, cp_orig, h_orig, s_orig)?;

    // Calculate inverse standard deviation weights to normalize different property scales
    // This ensures Cp [~30 J/molВ·K], dH [~-100 kJ/mol], and dS [~200 J/molВ·K] contribute equally
    let w_cp = 1.0 / property_scale(cp_orig);
    let w_h = 1.0 / property_scale(h_orig);
    let w_s = 1.0 / property_scale(s_orig);

    let rows = 3 * n;
    // Independent coordinates are A..E, Q=(F-H), and G.
    let cols = 7;
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
        X[(row_cp, 5)] = 0.0; // F-H
        X[(row_cp, 6)] = 0.0; // G
        y[row_cp] = cp_orig[i];

        //--------------------------------------------------------
        // H row
        //--------------------------------------------------------
        X[(row_h, 0)] = t;
        X[(row_h, 1)] = t2 / 2.0;
        X[(row_h, 2)] = t3 / 3.0;
        X[(row_h, 3)] = t4 / 4.0;
        X[(row_h, 4)] = -1.0 / t;
        X[(row_h, 5)] = 1.0; // F-H
        X[(row_h, 6)] = 0.0; // G
        y[row_h] = h_orig[i];

        //--------------------------------------------------------
        // S row
        //--------------------------------------------------------
        X[(row_s, 0)] = t.ln();
        X[(row_s, 1)] = t;
        X[(row_s, 2)] = t2 / 2.0;
        X[(row_s, 3)] = t3 / 3.0;
        X[(row_s, 4)] = -0.5 / t2;
        X[(row_s, 5)] = 0.0; // F-H
        X[(row_s, 6)] = 1.0; // G
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
    let largest_singular_value = svd.singular_values.iter().copied().fold(0.0_f64, f64::max);
    let rank_tolerance = largest_singular_value * 1e-12;
    let rank = svd
        .singular_values
        .iter()
        .filter(|&&value| value > rank_tolerance)
        .count();
    if rank < cols {
        return Err(NISTError::Fit(NISTFitError::Numerical(
            NISTFitNumericalError::RankDeficient {
                rank,
                required: cols,
            },
        )));
    }
    let theta = svd.solve(&y, rank_tolerance).map_err(|error| {
        NISTError::Fit(NISTFitError::Numerical(
            NISTFitNumericalError::LinearSolveFailed(error.to_string()),
        ))
    })?;

    let theta = theta.data.as_vec();
    let coeffs = vec![
        theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6], 0.0,
    ];
    if coeffs.iter().any(|value| !value.is_finite()) {
        return Err(NISTError::Fit(NISTFitError::Numerical(
            NISTFitNumericalError::NonFiniteCoefficients,
        )));
    }

    Ok(coeffs)
}
//------------------------- Helper Functions -------------------------------

fn property_scale(values: &[f64]) -> f64 {
    let rms = (values.iter().map(|value| value * value).sum::<f64>() / values.len() as f64).sqrt();
    rms.max(1.0)
}

fn reanchor_enthalpy_gauge(coeffs: &mut [f64], source_h: f64) {
    coeffs[5] += source_h;
    coeffs[7] = source_h;
}

#[cfg(test)]
mod offline_tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_fit_nist_scaled8_function() -> Result<(), NISTError> {
        let temps = vec![400.0, 500.0, 600.0, 700.0, 800.0];
        let cp_vals = vec![29.1, 30.2, 31.5, 32.8, 34.1];
        let h_vals = vec![-110.5, -107.2, -103.8, -100.3, -96.7];
        let s_vals = vec![197.7, 206.2, 213.8, 220.6, 226.7];

        let coeffs = fit_nist_scaled8(&temps, &cp_vals, &h_vals, &s_vals)?;
        assert_eq!(coeffs.len(), 8);
        assert!(coeffs[0] > 0.0);
        Ok(())
    }

    #[test]
    fn test_fit_nist_scaled8_rejects_bad_inputs() {
        let temps = [300.0, 400.0, 500.0];
        let cp = [10.0, 11.0];
        let h = [1.0, 2.0, 3.0];
        let s = [4.0, 5.0, 6.0];

        let err = fit_nist_scaled8(&temps, &cp, &h, &s).unwrap_err();
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Input(
                NISTFitInputError::LengthMismatch { .. }
            ))
        ));

        let cp_bad = [10.0, f64::INFINITY, 12.0];
        let err = fit_nist_scaled8(&temps, &cp_bad, &h, &s).unwrap_err();
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Input(
                NISTFitInputError::NonFiniteValue { .. }
            ))
        ));
    }

    #[test]
    fn test_fit_nist_scaled8_rejects_rank_deficient_design() {
        let temperatures = [500.0; 5];
        let values = [1.0; 5];
        let error = fit_nist_scaled8(&temperatures, &values, &values, &values)
            .expect_err("repeated temperatures cannot identify seven fit coordinates");

        assert!(matches!(
            error,
            NISTError::Fit(NISTFitError::Numerical(
                NISTFitNumericalError::RankDeficient { .. }
            ))
        ));
    }

    #[test]
    fn typed_fit_error_survives_thermo_facade_conversion() {
        let original = NISTError::Fit(NISTFitError::Input(NISTFitInputError::EmptySamples));
        assert!(std::error::Error::source(&original).is_some());

        let facade_error =
            crate::Thermodynamics::DBhandlers::thermo_api::ThermoError::from(original);
        assert!(matches!(
            facade_error,
            crate::Thermodynamics::DBhandlers::thermo_api::ThermoError::NISTFit(
                NISTFitError::Input(NISTFitInputError::EmptySamples)
            )
        ));
        assert!(std::error::Error::source(&facade_error).is_some());
    }

    #[test]
    fn test_direct_compare2_rejects_length_mismatch() {
        let nist = NISTdata::new();
        let err = nist
            .direct_compare2(&[1.0, 2.0], &[300.0], Expr::Const(1.0))
            .expect_err("length mismatch must be rejected explicitly");
        assert!(matches!(err, NISTError::FittingError(message) if message.contains("same length")));
    }

    #[test]
    fn test_direct_compare2_accepts_zero_reference_values_without_panicking() {
        let nist = NISTdata::new();
        let report = nist
            .direct_compare2(&[0.0, 0.0], &[300.0, 400.0], Expr::Const(0.0))
            .expect("zero-valued reference data should still produce a report");
        assert_eq!(report.number_of_significant_points, 0);
        assert!(report.l2_norm.is_finite());
        assert!(report.max_norm.is_finite());
    }

    #[test]
    fn test_validate_property_reports_accepts_good_reports() {
        let report_map = HashMap::from([
            (
                "Cp fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
            (
                "dh fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
            (
                "ds fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
        ]);

        assert!(
            NISTdata::validate_property_reports(
                "unit-test",
                &report_map,
                &[
                    "Cp fitting report",
                    "dh fitting report",
                    "ds fitting report"
                ],
            )
            .is_ok()
        );
    }

    #[test]
    fn test_validate_property_reports_rejects_large_error() {
        let report_map = HashMap::from([
            (
                "Cp fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 3.0e-2,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
            (
                "dh fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
            (
                "ds fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
        ]);

        let err = NISTdata::validate_property_reports(
            "unit-test",
            &report_map,
            &[
                "Cp fitting report",
                "dh fitting report",
                "ds fitting report",
            ],
        )
        .expect_err("large property error must be rejected");
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Quality(NISTFitQualityError::Rejected { property, .. }))
                if property.contains("Cp fitting report")
        ));
    }

    #[test]
    fn test_direct_compare2_rejects_empty_inputs() {
        let nist = NISTdata::new();
        let err = nist
            .direct_compare2(&[], &[], Expr::Const(0.0))
            .expect_err("empty vectors must be rejected explicitly");
        assert!(matches!(
            err,
            NISTError::FittingError(message) if message.contains("must not be empty")
        ));
    }

    #[test]
    fn test_fit_nist_scaled8_rejects_empty_vectors() {
        let err = fit_nist_scaled8(&[], &[], &[], &[])
            .expect_err("empty vectors must be rejected explicitly");
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Input(NISTFitInputError::EmptySamples))
        ));
    }

    #[test]
    fn test_prepare_fit_data_for_non_adjacent_fitting_rejects_empty_selection() {
        let error = NISTdata::prepare_fit_data_for_non_adjacent_fitting(
            vec![Coeffs {
                T: (1000.0, 1100.0),
                coeff: (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            }],
            10.0,
            20.0,
        )
        .expect_err("intervals outside the requested range must fail explicitly");

        assert!(matches!(
            error,
            NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::NoSamples { .. }
            ))
        ));
    }

    #[test]
    fn test_fit_nist_scaled8_reconstructs_one_exact_shomate_polynomial() {
        let source = (30.0, -3.0, 0.5, -0.05, 0.005, -8.0, 180.0, 5.0);
        let temperatures = (300..=900).step_by(25).map(f64::from).collect::<Vec<_>>();
        let cp = temperatures
            .iter()
            .map(|temperature| {
                let t = temperature / 1000.0;
                calculate_cp(t, source.0, source.1, source.2, source.3, source.4)
            })
            .collect::<Vec<_>>();
        let h = temperatures
            .iter()
            .map(|temperature| {
                let t = temperature / 1000.0;
                calculate_dh(
                    t, source.0, source.1, source.2, source.3, source.4, source.5, source.6,
                    source.7,
                )
            })
            .collect::<Vec<_>>();
        let s = temperatures
            .iter()
            .map(|temperature| {
                let t = temperature / 1000.0;
                calculate_s(
                    t, source.0, source.1, source.2, source.3, source.4, source.5, source.6,
                    source.7,
                )
            })
            .collect::<Vec<_>>();

        let fitted = fit_nist_scaled8(&temperatures, &cp, &h, &s)
            .expect("an exact Shomate polynomial must have an exact linear fit");

        assert_eq!(fitted.len(), 8);
        assert_eq!(fitted[7], 0.0, "the low-level fit uses the H=0 gauge");
        for (index, temperature) in temperatures.iter().copied().enumerate() {
            let t = temperature / 1000.0;
            let fitted_cp = calculate_cp(t, fitted[0], fitted[1], fitted[2], fitted[3], fitted[4]);
            let fitted_h = calculate_dh(
                t, fitted[0], fitted[1], fitted[2], fitted[3], fitted[4], fitted[5], fitted[6],
                fitted[7],
            );
            let fitted_s = calculate_s(
                t, fitted[0], fitted[1], fitted[2], fitted[3], fitted[4], fitted[5], fitted[6],
                fitted[7],
            );
            assert!((fitted_cp - cp[index]).abs() <= 1e-9);
            assert!((fitted_h - h[index]).abs() <= 1e-9);
            assert!((fitted_s - s[index]).abs() <= 1e-9);
        }
    }

    #[test]
    fn test_adjacent_svd_preserves_enthalpy_gauge_and_publishes_atomically() {
        let source = (30.0, -3.0, 0.5, -0.05, 0.005, -8.0, 180.0, 5.0);
        let mut nist = NISTdata::new();
        nist.fit_cp_dh_ds(
            Coeffs {
                T: (300.0, 600.0),
                coeff: source,
            },
            Coeffs {
                T: (600.0, 900.0),
                coeff: source,
            },
            350.0,
            850.0,
        )
        .expect("identical adjacent polynomials must fit exactly");

        let fitted = nist.coeffs.expect("accepted fit must be published");
        assert!((fitted.7 - source.7).abs() <= 1e-12);
        assert!((fitted.5 - fitted.7 - (source.5 - source.7)).abs() <= 1e-9);
    }

    #[test]
    fn adjacent_svd_two_range_fixture_reproduces_cp_h_s_across_join() {
        let lower = Coeffs {
            T: (300.0, 600.0),
            coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -8.0, 180.0, 5.0),
        };
        let upper = Coeffs {
            T: (600.0, 900.0),
            coeff: (
                30.02, -3.002, 0.5002, -0.05002, 0.005002, -3.98, 180.02, 9.0,
            ),
        };
        let mut nist = NISTdata::new();
        nist.fit_cp_dh_ds(lower.clone(), upper.clone(), 350.0, 850.0)
            .expect("nearby Shomate intervals should be accepted by the SVD stitcher");
        let fitted = nist
            .coeffs
            .expect("accepted coefficients must be published");

        let mut maximum_absolute_error = [0.0_f64; 3];
        for temperature in (350..=850).map(f64::from) {
            let source = if temperature <= 600.0 {
                lower.coeff
            } else {
                upper.coeff
            };
            let t = temperature / 1000.0;
            let reference = [
                calculate_cp(t, source.0, source.1, source.2, source.3, source.4),
                calculate_dh(
                    t, source.0, source.1, source.2, source.3, source.4, source.5, source.6,
                    source.7,
                ),
                calculate_s(
                    t, source.0, source.1, source.2, source.3, source.4, source.5, source.6,
                    source.7,
                ),
            ];
            let actual = [
                calculate_cp(t, fitted.0, fitted.1, fitted.2, fitted.3, fitted.4),
                calculate_dh(
                    t, fitted.0, fitted.1, fitted.2, fitted.3, fitted.4, fitted.5, fitted.6,
                    fitted.7,
                ),
                calculate_s(
                    t, fitted.0, fitted.1, fitted.2, fitted.3, fitted.4, fitted.5, fitted.6,
                    fitted.7,
                ),
            ];
            let absolute_tolerances = [0.05, 0.05, 0.05];
            for property in 0..3 {
                let error = (actual[property] - reference[property]).abs();
                maximum_absolute_error[property] = maximum_absolute_error[property].max(error);
                let allowed = absolute_tolerances[property] + 2.0e-3 * reference[property].abs();
                assert!(
                    error <= allowed,
                    "property {property} at {temperature} K: fitted={}, source={}, abs_error={}, allowed={}",
                    actual[property],
                    reference[property],
                    error,
                    allowed
                );
            }
        }

        assert!(maximum_absolute_error.iter().all(|error| *error > 0.0));
    }

    #[test]
    fn shomate_svd_and_rst_lm_use_the_same_weighted_three_property_problem() {
        use RustedSciThe::numerical::optimization::sym_fitting::Fitting;

        let lower = Coeffs {
            T: (300.0, 600.0),
            coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -8.0, 180.0, 5.0),
        };
        let upper = Coeffs {
            T: (600.0, 900.0),
            coeff: (
                30.02, -3.002, 0.5002, -0.05002, 0.005002, -3.98, 180.02, 9.0,
            ),
        };
        let (temperatures, cp, h, s) = NISTdata::prepare_fit_data_for_2_range_fitting(
            lower.clone(),
            upper.clone(),
            350.0,
            850.0,
        )
        .unwrap();
        let weights = [
            property_scale(&cp).recip(),
            property_scale(&h).recip(),
            property_scale(&s).recip(),
        ];

        let svd = fit_nist_scaled8(&temperatures, &cp, &h, &s).unwrap();

        let t = Expr::Var("T".to_string()) / Expr::Const(1000.0);
        let equations = [
            Expr::parse_expression("a + b*t + c*t^2 + d*t^3 + e/t^2"),
            Expr::parse_expression("a*t + b*t^2/2 + c*t^3/3 + d*t^4/4 - e/t + q"),
            Expr::parse_expression("a*ln(t) + b*t + c*t^2/2 + d*t^3/3 - e/(2*t^2) + g"),
        ]
        .into_iter()
        .map(|equation| equation.substitute_variable("t", &t))
        .zip(weights)
        .map(|(equation, weight)| Expr::Const(weight) * equation)
        .collect::<Vec<_>>();
        let mut targets = Vec::with_capacity(temperatures.len() * 3);
        targets.extend(cp.iter().map(|value| weights[0] * value));
        targets.extend(h.iter().map(|value| weights[1] * value));
        targets.extend(s.iter().map(|value| weights[2] * value));

        let initial = vec![
            0.5 * (lower.coeff.0 + upper.coeff.0),
            0.5 * (lower.coeff.1 + upper.coeff.1),
            0.5 * (lower.coeff.2 + upper.coeff.2),
            0.5 * (lower.coeff.3 + upper.coeff.3),
            0.5 * (lower.coeff.4 + upper.coeff.4),
            0.5 * ((lower.coeff.5 - lower.coeff.7) + (upper.coeff.5 - upper.coeff.7)),
            0.5 * (lower.coeff.6 + upper.coeff.6),
        ];
        let mut rst = Fitting::new()
            .with_data(temperatures.clone(), targets)
            .with_equations(equations)
            .with_unknowns(
                ["a", "b", "c", "d", "e", "q", "g"]
                    .map(str::to_owned)
                    .to_vec(),
            )
            .with_arg("T".to_string())
            .with_initial_guess(initial)
            .with_tolerance(1e-10)
            .with_f_tolerance(1e-12)
            .with_g_tolerance(1e-12)
            .with_max_iterations(1000);
        rst.solve()
            .expect("RST LM should solve the same weighted Shomate least-squares problem");
        let rst = rst
            .get_map_of_solutions()
            .expect("RST LM should publish fitted parameters");
        let rst = ["a", "b", "c", "d", "e", "q", "g"]
            .map(|name| rst[name])
            .to_vec();

        let mut max_prediction_delta = [0.0_f64; 3];
        let mut svd_objective = 0.0;
        let mut rst_objective = 0.0;
        for (index, temperature) in temperatures.iter().copied().enumerate() {
            let x = temperature / 1000.0;
            let svd_values = [
                calculate_cp(x, svd[0], svd[1], svd[2], svd[3], svd[4]),
                calculate_dh(
                    x, svd[0], svd[1], svd[2], svd[3], svd[4], svd[5], svd[6], 0.0,
                ),
                calculate_s(
                    x, svd[0], svd[1], svd[2], svd[3], svd[4], svd[5], svd[6], 0.0,
                ),
            ];
            let rst_values = [
                calculate_cp(x, rst[0], rst[1], rst[2], rst[3], rst[4]),
                calculate_dh(x, rst[0], rst[1], rst[2], rst[3], rst[4], rst[5], 0.0, 0.0),
                calculate_s(x, rst[0], rst[1], rst[2], rst[3], rst[4], 0.0, rst[6], 0.0),
            ];
            let source_values = [cp[index], h[index], s[index]];
            for property in 0..3 {
                let allowed = 0.05 + 0.02 * source_values[property].abs();
                assert!(
                    (svd_values[property] - source_values[property]).abs() <= allowed,
                    "SVD property {property} at {temperature} K exceeds shared acceptance"
                );
                assert!(
                    (rst_values[property] - source_values[property]).abs() <= allowed,
                    "RST property {property} at {temperature} K exceeds shared acceptance"
                );
                max_prediction_delta[property] = max_prediction_delta[property]
                    .max((svd_values[property] - rst_values[property]).abs());
                svd_objective +=
                    (weights[property] * (svd_values[property] - source_values[property])).powi(2);
                rst_objective +=
                    (weights[property] * (rst_values[property] - source_values[property])).powi(2);
            }
        }

        assert!(
            rst_objective <= svd_objective * (1.0 + 1e-5) + 1e-16,
            "same weighted objective: SVD={svd_objective:e}, RST LM={rst_objective:e}"
        );
        assert!(
            max_prediction_delta.iter().all(|delta| *delta < 1e-3),
            "componentwise property deltas (Cp,H,S) = {max_prediction_delta:?}"
        );
    }

    #[test]
    fn test_adjacent_svd_quality_rejection_preserves_previous_coefficients() {
        let previous = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);
        let mut nist = NISTdata::new();
        nist.coeffs = Some(previous);

        let error = nist
            .fit_cp_dh_ds(
                Coeffs {
                    T: (300.0, 600.0),
                    coeff: (10.0, 0.0, 0.0, 0.0, 0.0, -10.0, 50.0, 0.0),
                },
                Coeffs {
                    T: (600.0, 900.0),
                    coeff: (1000.0, 0.0, 0.0, 0.0, 0.0, 500.0, 500.0, 0.0),
                },
                350.0,
                850.0,
            )
            .expect_err("a discontinuous synthetic pair must fail the common quality gate");

        assert!(matches!(
            error,
            NISTError::Fit(NISTFitError::Quality(NISTFitQualityError::Rejected { .. }))
        ));
        assert_eq!(nist.coeffs, Some(previous));
    }

    #[test]
    fn test_prepare_fit_data_for_2_range_fitting_keeps_boundary_once() {
        let coeffs_min = Coeffs {
            T: (300.0, 500.0),
            coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
        };
        let coeffs_max = Coeffs {
            T: (500.0, 900.0),
            coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
        };

        let (temps, cp, h, s) =
            NISTdata::prepare_fit_data_for_2_range_fitting(coeffs_min, coeffs_max, 300.0, 700.0)
                .expect("boundary-aware preparation should succeed");

        assert_eq!(temps.len(), cp.len());
        assert_eq!(temps.len(), h.len());
        assert_eq!(temps.len(), s.len());
        assert_eq!(temps.iter().filter(|&&t| t == 500.0).count(), 1);
        assert!(temps.windows(2).all(|window| window[0] < window[1]));
    }

    #[test]
    fn test_validate_lm_fit_result_accepts_finite_solver_output() {
        let coeffs =
            NISTdata::validate_lm_fit_result("unit-test", Some(vec![1.0, 2.0, 3.0]), 3, Some(0.99))
                .expect("well-formed LM output should be accepted");

        assert_eq!(coeffs, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_validate_lm_fit_result_rejects_non_finite_or_weak_output() {
        let err =
            NISTdata::validate_lm_fit_result("unit-test", Some(vec![1.0, f64::NAN]), 2, Some(0.99))
                .expect_err("non-finite coefficients must be rejected");
        assert!(matches!(
            err,
            NISTError::FittingError(message) if message.contains("non-finite coefficient")
        ));

        let err = NISTdata::validate_lm_fit_result("unit-test", Some(vec![1.0, 2.0]), 2, Some(0.9))
            .expect_err("weak solver termination must be rejected");
        assert!(matches!(
            err,
            NISTError::FittingError(message) if message.contains("terminated too poorly")
        ));
    }
}
/// Apply weight to a specific row in the least squares system
fn scale_row(X: &mut DMatrix<f64>, y: &mut DVector<f64>, row: usize, w: f64) {
    let cols = X.ncols();
    for c in 0..cols {
        X[(row, c)] *= w;
    }
    y[row] *= w;
}

fn validate_fit_inputs(
    temps_K: &[f64],
    cp_orig: &[f64],
    h_orig: &[f64],
    s_orig: &[f64],
) -> Result<(), NISTError> {
    if temps_K.is_empty() {
        return Err(NISTError::Fit(NISTFitError::Input(
            NISTFitInputError::EmptySamples,
        )));
    }
    if temps_K.len() != cp_orig.len()
        || temps_K.len() != h_orig.len()
        || temps_K.len() != s_orig.len()
    {
        return Err(NISTError::Fit(NISTFitError::Input(
            NISTFitInputError::LengthMismatch {
                temperatures: temps_K.len(),
                cp: cp_orig.len(),
                h: h_orig.len(),
                s: s_orig.len(),
            },
        )));
    }

    let validate_slice = |name: &'static str, values: &[f64]| -> Result<(), NISTError> {
        if let Some((idx, value)) = values
            .iter()
            .enumerate()
            .find(|(_, value)| !value.is_finite())
        {
            return Err(NISTError::Fit(NISTFitError::Input(
                NISTFitInputError::NonFiniteValue {
                    field: name,
                    index: idx,
                    value: *value,
                },
            )));
        }
        Ok(())
    };

    validate_slice("temperature", temps_K)?;
    if let Some((idx, value)) = temps_K.iter().enumerate().find(|(_, value)| **value <= 0.0) {
        return Err(NISTError::Fit(NISTFitError::Input(
            NISTFitInputError::NonPositiveTemperature {
                index: idx,
                value: *value,
            },
        )));
    }
    validate_slice("Cp", cp_orig)?;
    validate_slice("H", h_orig)?;
    validate_slice("S", s_orig)?;
    Ok(())
}

fn validate_temperature_bounds(T_min: f64, T_max: f64) -> Result<(), NISTError> {
    if !T_min.is_finite() || !T_max.is_finite() || T_min <= 0.0 || T_max <= 0.0 {
        return Err(NISTError::Fit(NISTFitError::Preparation(
            NISTFitPreparationError::InvalidRange {
                minimum: T_min,
                maximum: T_max,
            },
        )));
    }
    if T_min >= T_max {
        return Err(NISTError::Fit(NISTFitError::Preparation(
            NISTFitPreparationError::InvalidRange {
                minimum: T_min,
                maximum: T_max,
            },
        )));
    }
    Ok(())
}

fn temperature_samples(start: f64, end: f64, include_start: bool) -> Result<Vec<f64>, NISTError> {
    if !start.is_finite() || !end.is_finite() || start <= 0.0 || end <= 0.0 || start > end {
        return Err(NISTError::Fit(NISTFitError::Preparation(
            NISTFitPreparationError::InvalidRange {
                minimum: start,
                maximum: end,
            },
        )));
    }

    let span = end - start;
    let estimated_samples = (span / FIT_SAMPLE_STEP_K).ceil() + 2.0;
    if !estimated_samples.is_finite() || estimated_samples > MAX_FIT_SAMPLES as f64 {
        return Err(NISTError::Fit(NISTFitError::Preparation(
            NISTFitPreparationError::SampleLimitExceeded {
                requested: estimated_samples,
                maximum: MAX_FIT_SAMPLES,
            },
        )));
    }
    let estimated_samples = estimated_samples as usize;

    if start == end {
        return Ok(if include_start {
            vec![start]
        } else {
            Vec::new()
        });
    }

    let mut samples = Vec::with_capacity(estimated_samples);
    if include_start {
        samples.push(start);
    }

    let tolerance = 8.0 * f64::EPSILON * start.abs().max(end.abs()).max(1.0);
    let mut offset = FIT_SAMPLE_STEP_K;
    while start + offset < end - tolerance {
        samples.push(start + offset);
        offset += FIT_SAMPLE_STEP_K;
    }

    if samples.last().is_some_and(|last| end - *last <= tolerance) {
        *samples.last_mut().expect("last sample was checked") = end;
    } else {
        samples.push(end);
    }
    Ok(samples)
}

fn validate_coeffs_finite(label: &str, coeffs: &[f64]) -> Result<(), NISTError> {
    if let Some((idx, value)) = coeffs
        .iter()
        .enumerate()
        .find(|(_, value)| !value.is_finite())
    {
        return Err(NISTError::Fit(NISTFitError::Input(
            NISTFitInputError::NonFiniteValue {
                field: if label.contains("interval") {
                    "source coefficient"
                } else {
                    "coefficient"
                },
                index: idx,
                value: *value,
            },
        )));
    }
    Ok(())
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
        validate_temperature_bounds(T_min, T_max)?;
        let (_, T_center) = (coeffs_min.T.0, coeffs_min.T.1);
        let (T_center_check, _) = (coeffs_max.T.0, coeffs_max.T.1);
        if T_center != T_center_check {
            return Err(NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::InvalidRange {
                    minimum: T_center,
                    maximum: T_center_check,
                },
            )));
        }
        if !coeffs_min.T.0.is_finite()
            || !coeffs_min.T.1.is_finite()
            || !coeffs_max.T.0.is_finite()
            || !coeffs_max.T.1.is_finite()
            || coeffs_min.T.0 <= 0.0
            || coeffs_max.T.0 <= 0.0
            || coeffs_min.T.0 >= coeffs_min.T.1
            || coeffs_max.T.0 >= coeffs_max.T.1
        {
            return Err(NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::InvalidRange {
                    minimum: coeffs_min.T.0,
                    maximum: coeffs_max.T.1,
                },
            )));
        }
        if coeffs_min.T.0 <= 0.0
            || coeffs_max.T.0 <= 0.0
            || T_min < coeffs_min.T.0
            || T_max > coeffs_max.T.1
            || T_min > T_center
            || T_max < T_center
        {
            return Err(NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::InvalidRange {
                    minimum: T_min,
                    maximum: T_max,
                },
            )));
        }
        validate_coeffs_finite(
            "NIST lower interval",
            &[
                coeffs_min.coeff.0,
                coeffs_min.coeff.1,
                coeffs_min.coeff.2,
                coeffs_min.coeff.3,
                coeffs_min.coeff.4,
                coeffs_min.coeff.5,
                coeffs_min.coeff.6,
                coeffs_min.coeff.7,
            ],
        )?;
        validate_coeffs_finite(
            "NIST upper interval",
            &[
                coeffs_max.coeff.0,
                coeffs_max.coeff.1,
                coeffs_max.coeff.2,
                coeffs_max.coeff.3,
                coeffs_max.coeff.4,
                coeffs_max.coeff.5,
                coeffs_max.coeff.6,
                coeffs_max.coeff.7,
            ],
        )?;

        let a = coeffs_min.coeff.0;
        let b = coeffs_min.coeff.1;
        let c = coeffs_min.coeff.2;
        let d = coeffs_min.coeff.3;
        let e = coeffs_min.coeff.4;
        let f = coeffs_min.coeff.5;
        let g = coeffs_min.coeff.6;
        let h = coeffs_min.coeff.7;

        let mut temps = temperature_samples(T_min, T_center, true)?;
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

        let upper_temps = temperature_samples(T_center, T_max, false)?;
        let cp_upper = upper_temps
            .iter()
            .map(|&temperature| calculate_cp(temperature / 1000.0, a, b, c, d, e))
            .collect::<Vec<_>>();
        let h_upper = upper_temps
            .iter()
            .map(|&temperature| calculate_dh(temperature / 1000.0, a, b, c, d, e, f, g, h))
            .collect::<Vec<_>>();
        let s_upper = upper_temps
            .iter()
            .map(|&temperature| calculate_s(temperature / 1000.0, a, b, c, d, e, f, g, h))
            .collect::<Vec<_>>();
        temps.extend(upper_temps);
        cp_vals.extend(cp_upper);
        h_vals.extend(h_upper);
        s_vals.extend(s_upper);

        if temps.len() != cp_vals.len()
            || temps.len() != h_vals.len()
            || temps.len() != s_vals.len()
        {
            return Err(NISTError::FittingError(
                "adjacent fit data generation produced mismatched vector lengths".to_string(),
            ));
        }
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
        let source_h = coeffs_min.coeff.7;
        let (temps, cp_vals, h_vals, s_vals) =
            Self::prepare_fit_data_for_2_range_fitting(coeffs_min, coeffs_max, T_min, T_max)?;

        let mut coeffs = fit_nist_scaled8(&temps, &cp_vals, &h_vals, &s_vals)?;
        reanchor_enthalpy_gauge(&mut coeffs, source_h);
        if coeffs.len() != 8 {
            return Err(NISTError::FittingError(format!(
                "NIST solver returned {} coefficients instead of 8",
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
        let h = coeffs[7];
        let Cp_sym = calculate_cp_sym(a, b, c, d, e);
        let dh_sym = calculate_dh_sym(a, b, c, d, e, f, g, h);
        let ds_sym = calculate_s_sym(a, b, c, d, e, f, g, h);
        let c_fitting_report =
            Self::compare_nist_property(&cp_vals, &temps, Cp_sym, NIST_CP_ABSOLUTE_TOLERANCE)?;
        report_map.insert("Cp fitting report".to_string(), c_fitting_report);
        let dh_fitting_report = Self::compare_nist_property(
            &h_vals,
            &temps,
            dh_sym,
            NIST_H_ABSOLUTE_TOLERANCE_KJ_PER_MOL,
        )?;
        report_map.insert("dh fitting report".to_string(), dh_fitting_report);
        let ds_fitting_report =
            Self::compare_nist_property(&s_vals, &temps, ds_sym, NIST_S_ABSOLUTE_TOLERANCE)?;
        report_map.insert("ds fitting report".to_string(), ds_fitting_report);

        Self::validate_property_reports(
            "NIST adjacent SVD",
            &report_map,
            &[
                "Cp fitting report",
                "dh fitting report",
                "ds fitting report",
            ],
        )?;
        self.coeffs = Some((a, b, c, d, e, f, g, h));
        Ok(())
    }

    pub fn fitting_cp_dh_ds_non_adjacent(
        &mut self,
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NISTError> {
        let mut report_map: HashMap<String, FittingReport> = HashMap::new();
        validate_temperature_bounds(T_min, T_max)?;
        let (Tvec, cp_vec, h_vec, s_vec) =
            Self::prepare_fit_data_for_non_adjacent_fitting(coeffs.clone(), T_min, T_max)?;
        let source_h = coeffs
            .iter()
            .filter(|coefficients| coefficients.T.1 >= T_min && coefficients.T.0 <= T_max)
            .min_by(|left, right| left.T.0.total_cmp(&right.T.0))
            .map(|coefficients| coefficients.coeff.7)
            .ok_or_else(|| {
                NISTError::Fit(NISTFitError::Preparation(
                    NISTFitPreparationError::NoSamples {
                        minimum: T_min,
                        maximum: T_max,
                    },
                ))
            })?;

        let mut coeffs = fit_nist_scaled8(&Tvec, &cp_vec, &h_vec, &s_vec)?;
        reanchor_enthalpy_gauge(&mut coeffs, source_h);
        if coeffs.len() != 8 {
            return Err(NISTError::FittingError(format!(
                "NIST solver returned {} coefficients instead of 8",
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
        let h = coeffs[7];
        let Cp_sym = calculate_cp_sym(a, b, c, d, e);
        let dh_sym = calculate_dh_sym(a, b, c, d, e, f, g, h);
        let ds_sym = calculate_s_sym(a, b, c, d, e, f, g, h);
        let c_fitting_report =
            Self::compare_nist_property(&cp_vec, &Tvec, Cp_sym, NIST_CP_ABSOLUTE_TOLERANCE)?;
        report_map.insert("Cp fitting report".to_string(), c_fitting_report);
        let dh_fitting_report = Self::compare_nist_property(
            &h_vec,
            &Tvec,
            dh_sym,
            NIST_H_ABSOLUTE_TOLERANCE_KJ_PER_MOL,
        )?;
        report_map.insert("dh fitting report".to_string(), dh_fitting_report);
        let ds_fitting_report =
            Self::compare_nist_property(&s_vec, &Tvec, ds_sym, NIST_S_ABSOLUTE_TOLERANCE)?;
        report_map.insert("ds fitting report".to_string(), ds_fitting_report);

        Self::validate_property_reports(
            "NIST non-adjacent SVD",
            &report_map,
            &[
                "Cp fitting report",
                "dh fitting report",
                "ds fitting report",
            ],
        )?;
        self.coeffs = Some((a, b, c, d, e, f, g, h));
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
    /// Compatibility-only symbolic fit. Production stitching uses
    /// [`Self::fit_cp_dh_ds`], because the Shomate model is coefficient-linear.
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
        )
        .map_err(|source| NISTError::FittingBackend {
            stage: "adjacent Cp/S",
            source,
        })?;

        let coeffs = Self::validate_lm_fit_result(
            "combined fitting",
            sew.get_result(),
            6,
            sew.get_r_ssquared(),
        )?;
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
        )
        .map_err(|source| NISTError::FittingBackend {
            stage: "adjacent enthalpy",
            source,
        })?;
        let coeffs = Self::validate_lm_fit_result(
            "enthalpy stage",
            sew2.get_result(),
            2,
            sew2.get_r_ssquared(),
        )?;
        let f = coeffs[0];
        let h = coeffs[1];
        map_of_coeffs.insert("f".to_string(), f);
        map_of_coeffs.insert("h".to_string(), h);
        let Cp_sym = calculate_cp_sym(a, b, c, d, e);
        let dh_sym = calculate_dh_sym(a, b, c, d, e, f, g, h) / Expr::Const(1000.0);
        let ds_sym = calculate_s_sym(a, b, c, d, e, f, g, h);
        let combined_func_to_fit = combined_func_to_fit.set_variable_from_map(&map_of_coeffs);
        let comb_func_fitting_report =
            self.direct_compare2(&y_data.data.as_vec().clone(), &temps, combined_func_to_fit)?;
        report_map.insert("weighted function".to_string(), comb_func_fitting_report);
        let c_fitting_report = self.direct_compare2(&cp_vals, &temps, Cp_sym)?;
        report_map.insert("Cp fitting report".to_string(), c_fitting_report);
        let dh_fitting_report = self.direct_compare2(&h_vals, &temps, dh_sym)?;
        report_map.insert("dh fitting report".to_string(), dh_fitting_report);
        let ds_fitting_report = self.direct_compare2(&s_vals, &temps, ds_sym)?;
        report_map.insert("ds fitting report".to_string(), ds_fitting_report);

        self.coeffs = Some((a, b, c, d, e, f, g, h));

        Ok(())
    }

    /// Non-adjacent fitting method for multiple temperature intervals
    /// Compatibility-only symbolic fit. Production stitching uses
    /// [`Self::fitting_cp_dh_ds_non_adjacent`].
    pub fn fitting_non_adjacent(
        &mut self,
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<HashMap<String, FittingReport>, NISTError> {
        use RustedSciThe::numerical::optimization::fitting_features::SewTwoFunctions;
        use RustedSciThe::numerical::optimization::sym_fitting::Fitting;

        if coeffs.is_empty() {
            return Err(NISTError::FittingError(
                "NIST non-adjacent fitting requires at least one interval".to_string(),
            ));
        }

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
        )
        .map_err(|source| NISTError::FittingBackend {
            stage: "non-adjacent Cp/S",
            source,
        })?;

        let coeffs_step1 = Self::validate_lm_fit_result(
            "NIST LM stage",
            sew.get_result(),
            6,
            sew.get_r_ssquared(),
        )?;
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
        )
        .map_err(|source| NISTError::FittingBackend {
            stage: "non-adjacent enthalpy",
            source,
        })?;

        let coeffs_step2 = Self::validate_lm_fit_result(
            "NIST enthalpy stage",
            sew2.get_result(),
            2,
            sew2.get_r_ssquared(),
        )?;
        let (f, h) = (coeffs_step2[0], coeffs_step2[1]);

        // Final coefficients and reports
        let Cp_sym = calculate_cp_sym(a, b, c, d, e);
        let dh_sym = calculate_dh_sym(a, b, c, d, e, f, g, h) / Expr::Const(1000.0);
        let ds_sym = calculate_s_sym(a, b, c, d, e, f, g, h);

        let c_fitting_report = self.direct_compare2(&cp_vals, &temps, Cp_sym)?;
        report_map.insert("Cp fitting report".to_string(), c_fitting_report);
        let dh_fitting_report = self.direct_compare2(&h_vals, &temps, dh_sym)?;
        report_map.insert("dh fitting report".to_string(), dh_fitting_report);
        let ds_fitting_report = self.direct_compare2(&s_vals, &temps, ds_sym)?;
        report_map.insert("ds fitting report".to_string(), ds_fitting_report);

        self.coeffs = Some((a, b, c, d, e, f, g, h));
        Ok(report_map)
    }

    fn prepare_fit_data_for_non_adjacent_fitting(
        coeffs: Vec<Coeffs>,
        T_min: f64,
        T_max: f64,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), NISTError> {
        validate_temperature_bounds(T_min, T_max)?;
        let mut coeffs = coeffs;
        for coeff in &coeffs {
            if !coeff.T.0.is_finite()
                || !coeff.T.1.is_finite()
                || coeff.T.0 <= 0.0
                || coeff.T.0 >= coeff.T.1
            {
                return Err(NISTError::Fit(NISTFitError::Preparation(
                    NISTFitPreparationError::InvalidRange {
                        minimum: coeff.T.0,
                        maximum: coeff.T.1,
                    },
                )));
            }
            validate_coeffs_finite(
                "NIST interval",
                &[
                    coeff.coeff.0,
                    coeff.coeff.1,
                    coeff.coeff.2,
                    coeff.coeff.3,
                    coeff.coeff.4,
                    coeff.coeff.5,
                    coeff.coeff.6,
                    coeff.coeff.7,
                ],
            )?;
        }
        coeffs.sort_by(|left, right| left.T.0.total_cmp(&right.T.0));

        let range_tolerance = 8.0 * f64::EPSILON * T_max.abs().max(1.0);
        for pair in coeffs.windows(2) {
            if pair[1].T.0 < pair[0].T.1 - range_tolerance {
                return Err(NISTError::Fit(NISTFitError::Preparation(
                    NISTFitPreparationError::OverlappingIntervals {
                        left: pair[0].T,
                        right: pair[1].T,
                    },
                )));
            }
        }

        let mut Tvec: Vec<f64> = Vec::new();
        let mut cp_vec: Vec<f64> = Vec::new();
        let mut h_vec: Vec<f64> = Vec::new();
        let mut s_vec: Vec<f64> = Vec::new();
        let mut previous_overlap_end = None;

        for coeff in &coeffs {
            let sample_start = coeff.T.0.max(T_min);
            let sample_end = coeff.T.1.min(T_max);
            if sample_start > sample_end {
                continue;
            }

            let owns_start = previous_overlap_end.is_none_or(|previous_end: f64| {
                (sample_start - previous_end).abs() > range_tolerance
            });
            let temps = temperature_samples(sample_start, sample_end, owns_start)?;
            previous_overlap_end = Some(sample_end);

            let (a, b, c, d, e, f, g, h) = coeff.coeff;
            Tvec.extend(temps.iter().copied());
            cp_vec.extend(
                temps
                    .iter()
                    .map(|&temperature| calculate_cp(temperature / 1000.0, a, b, c, d, e)),
            );
            h_vec.extend(
                temps
                    .iter()
                    .map(|&temperature| calculate_dh(temperature / 1000.0, a, b, c, d, e, f, g, h)),
            );
            s_vec.extend(
                temps
                    .iter()
                    .map(|&temperature| calculate_s(temperature / 1000.0, a, b, c, d, e, f, g, h)),
            );
        }

        if Tvec.is_empty() {
            return Err(NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::NoSamples {
                    minimum: T_min,
                    maximum: T_max,
                },
            )));
        }
        if Tvec.len() != cp_vec.len() || Tvec.len() != h_vec.len() || Tvec.len() != s_vec.len() {
            return Err(NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::InconsistentSampleLengths,
            )));
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
    ) -> Result<FittingReport, NISTError> {
        let report = compare_property_series(y_data, T_vec, fit_function)
            .map_err(|err| NISTError::FittingError(err.to_string()))?;
        Ok(FittingReport {
            l2_norm: report.l2_norm,
            max_norm: report.max_norm,
            number_of_significant_points: report.number_of_significant_points,
            T_range: report.t_range,
        })
    }

    fn compare_nist_property(
        y_data: &[f64],
        T_vec: &[f64],
        fit_function: Expr,
        absolute_tolerance: f64,
    ) -> Result<FittingReport, NISTError> {
        let report = compare_property_series_with_absolute_relative_tolerance(
            y_data,
            T_vec,
            fit_function,
            absolute_tolerance,
            MAX_PROPERTY_RELATIVE_ERROR,
        )
        .map_err(|err| {
            NISTError::Fit(NISTFitError::Input(NISTFitInputError::ComparisonFailed(
                err.to_string(),
            )))
        })?;
        Ok(FittingReport {
            l2_norm: report.l2_norm,
            max_norm: report.max_norm,
            number_of_significant_points: report.number_of_significant_points,
            T_range: report.t_range,
        })
    }

    fn validate_lm_fit_result<C>(
        context: &str,
        coeffs: Option<C>,
        expected_len: usize,
        r_squared: Option<f64>,
    ) -> Result<Vec<f64>, NISTError>
    where
        C: AsRef<[f64]>,
    {
        let r_squared = r_squared.ok_or_else(|| {
            NISTError::FittingError(format!("{context} did not publish an R^2 value"))
        })?;
        if !r_squared.is_finite() {
            return Err(NISTError::FittingError(format!(
                "{context} produced a non-finite R^2 value"
            )));
        }
        if r_squared < MIN_LM_R_SQUARED {
            return Err(NISTError::FittingError(format!(
                "{context} terminated too poorly: R^2 = {:.6} < {:.6}",
                r_squared, MIN_LM_R_SQUARED
            )));
        }

        let coeffs = coeffs.ok_or_else(|| {
            NISTError::FittingError(format!("{context} did not publish fitted coefficients"))
        })?;
        let coeffs = coeffs.as_ref();
        if coeffs.len() != expected_len {
            return Err(NISTError::FittingError(format!(
                "{context} returned {} coefficients instead of {}",
                coeffs.len(),
                expected_len
            )));
        }
        if let Some((idx, value)) = coeffs
            .iter()
            .enumerate()
            .find(|(_, value)| !value.is_finite())
        {
            return Err(NISTError::FittingError(format!(
                "{context} produced a non-finite coefficient at index {}: {}",
                idx, value
            )));
        }

        Ok(coeffs.to_vec())
    }
    fn validate_property_reports(
        context: &str,
        report_map: &HashMap<String, FittingReport>,
        required_reports: &[&str],
    ) -> Result<(), NISTError> {
        for report_name in required_reports {
            let report = report_map.get(*report_name).ok_or_else(|| {
                NISTError::Fit(NISTFitError::Quality(NISTFitQualityError::MissingReport {
                    property: (*report_name).to_string(),
                }))
            })?;

            if !report.l2_norm.is_finite() || !report.max_norm.is_finite() {
                return Err(NISTError::Fit(NISTFitError::Quality(
                    NISTFitQualityError::NonFiniteMetrics {
                        property: (*report_name).to_string(),
                    },
                )));
            }

            if report.max_norm > MAX_PROPERTY_RELATIVE_ERROR {
                return Err(NISTError::Fit(NISTFitError::Quality(
                    NISTFitQualityError::Rejected {
                        property: format!("{context} {report_name}"),
                        max_error: report.max_norm,
                        tolerance: MAX_PROPERTY_RELATIVE_ERROR,
                    },
                )));
            }

            if report.number_of_significant_points != 0 {
                return Err(NISTError::Fit(NISTFitError::Quality(
                    NISTFitQualityError::SignificantResiduals {
                        property: format!("{context} {report_name}"),
                        significant_points: report.number_of_significant_points,
                    },
                )));
            }
        }

        Ok(())
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
    #[ignore]
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
    #[ignore]
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
    #[ignore]
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
    #[ignore]
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
    #[ignore]
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
    #[ignore]
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
        let error = nist
            .fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 800.0)
            .expect_err("overlapping source intervals are ambiguous");
        assert!(matches!(
            error,
            NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::OverlappingIntervals { .. }
            ))
        ));
    }

    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_partial_overlap() {
        let mut nist = NISTdata::new();
        let source = (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0);
        let coeffs_vec = vec![
            Coeffs {
                T: (200.0, 500.0),
                coeff: source,
            },
            Coeffs {
                T: (700.0, 1200.0),
                coeff: source,
            },
        ];

        let reports = nist
            .fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 1000.0)
            .expect("compatible disjoint intervals clipped by the requested range should fit");
        assert!(nist.coeffs.is_some());
        assert!(reports.values().all(|report| report.max_norm <= 1e-9));
    }

    #[test]
    fn test_fitting_cp_dh_ds_non_adjacent_uses_both_intervals() {
        let source = (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0);
        let coeffs1 = Coeffs {
            T: (300.0, 600.0),
            coeff: source,
        };
        let coeffs2 = Coeffs {
            T: (700.0, 1000.0),
            coeff: source,
        };

        let mut nist = NISTdata::new();
        let reports = nist
            .fitting_cp_dh_ds_non_adjacent(vec![coeffs1, coeffs2], 400.0, 900.0)
            .expect("both disjoint intervals share one exact source polynomial");
        assert!(nist.coeffs.is_some());
        assert!(reports.values().all(|report| report.max_norm <= 1e-9));

        let (temps, _, _, _) = NISTdata::prepare_fit_data_for_non_adjacent_fitting(
            vec![
                Coeffs {
                    T: (300.0, 600.0),
                    coeff: source,
                },
                Coeffs {
                    T: (700.0, 1000.0),
                    coeff: source,
                },
            ],
            400.0,
            900.0,
        )
        .expect("each interval contributes samples in the requested range");
        assert!(temps.iter().any(|&temperature| temperature <= 600.0));
        assert!(temps.iter().any(|&temperature| temperature >= 700.0));
    }

    #[test]
    #[ignore]
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

        let source = (25.0, -2.0, 0.3, -0.03, 0.003, -500.0, 150.0, 0.0);
        let coeffs_vec = vec![
            Coeffs {
                T: (300.0, 500.0),
                coeff: source,
            },
            Coeffs {
                T: (600.0, 900.0),
                coeff: source,
            },
        ];

        let reports = nist
            .fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 800.0)
            .expect("compatible disjoint intervals should produce an accepted fit");
        assert!(reports.values().all(|report| report.max_norm <= 1e-9));

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
    }
    ////////////////////////////////////////////////////////////////////////
    #[test]
    #[ignore]
    fn test_fit_nist_scaled8_function() -> Result<(), NISTError> {
        let temps = vec![400.0, 500.0, 600.0, 700.0, 800.0];
        let cp_vals = vec![29.1, 30.2, 31.5, 32.8, 34.1];
        let h_vals = vec![-110.5, -107.2, -103.8, -100.3, -96.7];
        let s_vals = vec![197.7, 206.2, 213.8, 220.6, 226.7];

        let coeffs = fit_nist_scaled8(&temps, &cp_vals, &h_vals, &s_vals)?;
        assert_eq!(coeffs.len(), 8);

        assert!(coeffs[0] > 0.0);
        sleep(smtime);
        Ok(())
    }

    #[test]
    #[ignore]
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
    #[ignore]
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
    #[ignore]
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
    #[ignore]
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
    #[ignore]
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
    #[ignore]
    fn test_fitting_non_adjacent_sequential_vs_simultaneous() {
        // Keep the historical test name so existing filtered commands continue
        // to work. The production contract is now the deterministic linear
        // route; algorithm disagreement is not a scientific acceptance metric.
        let mut nist = NISTdata::new();
        let source = (30.0, -3.0, 0.5, -0.05, 0.005, -8.0, 180.0, 5.0);
        let coeffs_vec = vec![
            Coeffs {
                T: (300.0, 600.0),
                coeff: source,
            },
            Coeffs {
                T: (700.0, 1000.0),
                coeff: source,
            },
        ];

        let reports = nist
            .fitting_cp_dh_ds_non_adjacent(coeffs_vec, 400.0, 900.0)
            .expect("one exact polynomial over disjoint ranges must fit linearly");
        for property in [
            "Cp fitting report",
            "dh fitting report",
            "ds fitting report",
        ] {
            let report = reports.get(property).expect("missing property report");
            assert!(report.max_norm <= 1e-9, "{property}: {report:?}");
            assert_eq!(report.number_of_significant_points, 0);
        }

        let fitted = nist.coeffs.expect("accepted fit must publish coefficients");
        assert!((fitted.5 - fitted.7 - (source.5 - source.7)).abs() <= 1e-9);
    }

    #[test]
    #[ignore]
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
    #[ignore]
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
    fn non_adjacent_sampling_is_ordered_deduplicated_and_keeps_fractional_bounds() {
        let lower = Coeffs {
            T: (300.0, 700.0),
            coeff: (30.0, 0.0, 0.0, 0.0, 0.0, -800.0, 180.0, 0.0),
        };
        let upper = Coeffs {
            T: (700.0, 900.0),
            coeff: (40.0, 0.0, 0.0, 0.0, 0.0, -1500.0, 250.0, 0.0),
        };

        let (temps, cp, h, s) = NISTdata::prepare_fit_data_for_non_adjacent_fitting(
            vec![upper, lower.clone()],
            450.25,
            750.5,
        )
        .expect("valid disjoint intervals should produce fit samples");

        assert_eq!(temps.first(), Some(&450.25));
        assert_eq!(temps.last(), Some(&750.5));
        assert!(temps.windows(2).all(|pair| pair[0] < pair[1]));
        assert_eq!(
            temps
                .iter()
                .filter(|&&temperature| temperature == 700.0)
                .count(),
            1
        );
        assert_eq!(cp.len(), temps.len());
        assert_eq!(h.len(), temps.len());
        assert_eq!(s.len(), temps.len());
        let boundary_index = temps
            .iter()
            .position(|&temperature| temperature == 700.0)
            .unwrap();
        assert_eq!(
            cp[boundary_index],
            calculate_cp(0.7, 30.0, 0.0, 0.0, 0.0, 0.0)
        );
    }

    #[test]
    fn non_adjacent_sampling_rejects_overlapping_source_intervals() {
        let intervals = vec![
            Coeffs {
                T: (300.0, 710.0),
                coeff: (30.0, 0.0, 0.0, 0.0, 0.0, -800.0, 180.0, 0.0),
            },
            Coeffs {
                T: (700.0, 900.0),
                coeff: (40.0, 0.0, 0.0, 0.0, 0.0, -1500.0, 250.0, 0.0),
            },
        ];

        let error = NISTdata::prepare_fit_data_for_non_adjacent_fitting(intervals, 400.0, 800.0)
            .expect_err("ambiguous overlapping source intervals must be rejected");
        assert!(matches!(
            error,
            NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::OverlappingIntervals { .. }
            ))
        ));
    }

    #[test]
    fn nist_fit_rejects_non_positive_sample_temperature() {
        let error = fit_nist_scaled8(&[-1.0, 300.0], &[1.0, 1.0], &[1.0, 1.0], &[1.0, 1.0])
            .expect_err("Shomate fitting cannot use non-positive temperatures");
        assert!(matches!(
            error,
            NISTError::Fit(NISTFitError::Input(
                NISTFitInputError::NonPositiveTemperature { .. }
            ))
        ));
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

        let error = NISTdata::prepare_fit_data_for_non_adjacent_fitting(coeffs_vec, 400.0, 800.0)
            .expect_err("an empty selected range is not a valid fit input");
        assert!(matches!(
            error,
            NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::NoSamples { .. }
            ))
        ));
    }

    #[test]
    #[ignore]
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

    #[test]
    fn test_fit_nist_scaled8_rejects_bad_inputs() {
        let temps = [300.0, 400.0, 500.0];
        let cp = [10.0, 11.0];
        let h = [1.0, 2.0, 3.0];
        let s = [4.0, 5.0, 6.0];

        let err = fit_nist_scaled8(&temps, &cp, &h, &s).unwrap_err();
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Input(
                NISTFitInputError::LengthMismatch { .. }
            ))
        ));

        let cp_bad = [10.0, f64::INFINITY, 12.0];
        let err = fit_nist_scaled8(&temps, &cp_bad, &h, &s).unwrap_err();
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Input(
                NISTFitInputError::NonFiniteValue { .. }
            ))
        ));
    }

    #[test]
    #[ignore]
    fn test_direct_compare2_rejects_length_mismatch() {
        let nist = NISTdata::new();
        let err = nist
            .direct_compare2(&[1.0, 2.0], &[300.0], Expr::Const(1.0))
            .expect_err("length mismatch must be rejected explicitly");
        assert!(matches!(err, NISTError::FittingError(message) if message.contains("same length")));
    }

    #[test]
    #[ignore]
    fn test_direct_compare2_accepts_zero_reference_values_without_panicking() {
        let nist = NISTdata::new();
        let report = nist
            .direct_compare2(&[0.0, 0.0], &[300.0, 400.0], Expr::Const(0.0))
            .expect("zero-valued reference data should still produce a report");
        assert_eq!(report.number_of_significant_points, 0);
        assert!(report.l2_norm.is_finite());
        assert!(report.max_norm.is_finite());
    }

    #[test]
    #[ignore]
    fn test_validate_property_reports_accepts_good_reports() {
        let report_map = HashMap::from([
            (
                "Cp fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
            (
                "dh fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
            (
                "ds fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
        ]);

        assert!(
            NISTdata::validate_property_reports(
                "unit-test",
                &report_map,
                &[
                    "Cp fitting report",
                    "dh fitting report",
                    "ds fitting report"
                ],
            )
            .is_ok()
        );
    }

    #[test]
    #[ignore]
    fn test_validate_property_reports_rejects_large_error() {
        let report_map = HashMap::from([
            (
                "Cp fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 3.0e-2,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
            (
                "dh fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
            (
                "ds fitting report".to_string(),
                FittingReport {
                    l2_norm: 1.0e-6,
                    max_norm: 1.0e-4,
                    number_of_significant_points: 0,
                    T_range: Some((300.0, 1000.0)),
                },
            ),
        ]);

        let err = NISTdata::validate_property_reports(
            "unit-test",
            &report_map,
            &[
                "Cp fitting report",
                "dh fitting report",
                "ds fitting report",
            ],
        )
        .expect_err("large property error must be rejected");
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Quality(NISTFitQualityError::Rejected { property, .. }))
                if property.contains("Cp fitting report")
        ));
    }

    #[test]
    #[ignore]
    fn test_direct_compare2_rejects_empty_inputs() {
        let nist = NISTdata::new();
        let err = nist
            .direct_compare2(&[], &[], Expr::Const(0.0))
            .expect_err("empty vectors must be rejected explicitly");
        assert!(matches!(
            err,
            NISTError::FittingError(message) if message.contains("must not be empty")
        ));
    }

    #[test]
    #[ignore]
    fn test_fit_nist_scaled8_rejects_empty_vectors() {
        let err = fit_nist_scaled8(&[], &[], &[], &[])
            .expect_err("empty vectors must be rejected explicitly");
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Input(NISTFitInputError::EmptySamples))
        ));
    }

    #[test]
    #[ignore]
    fn test_prepare_fit_data_for_non_adjacent_fitting_rejects_empty_selection() {
        let err = NISTdata::prepare_fit_data_for_non_adjacent_fitting(
            vec![Coeffs {
                T: (1000.0, 1100.0),
                coeff: (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            }],
            10.0,
            20.0,
        )
        .expect_err("intervals outside the requested range must fail explicitly");
        assert!(matches!(
            err,
            NISTError::Fit(NISTFitError::Preparation(
                NISTFitPreparationError::NoSamples { .. }
            ))
        ));
    }

    #[test]
    #[ignore]
    fn test_prepare_fit_data_for_2_range_fitting_keeps_boundary_once() {
        let coeffs_min = Coeffs {
            T: (300.0, 500.0),
            coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
        };
        let coeffs_max = Coeffs {
            T: (500.0, 900.0),
            coeff: (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0),
        };

        let (temps, cp, h, s) =
            NISTdata::prepare_fit_data_for_2_range_fitting(coeffs_min, coeffs_max, 300.0, 700.0)
                .expect("boundary-aware preparation should succeed");

        assert_eq!(temps.len(), cp.len());
        assert_eq!(temps.len(), h.len());
        assert_eq!(temps.len(), s.len());
        assert_eq!(temps.iter().filter(|&&t| t == 500.0).count(), 1);
        assert!(temps.windows(2).all(|window| window[0] < window[1]));
    }

    #[test]
    #[ignore]
    fn test_validate_lm_fit_result_accepts_finite_solver_output() {
        let coeffs =
            NISTdata::validate_lm_fit_result("unit-test", Some(vec![1.0, 2.0, 3.0]), 3, Some(0.99))
                .expect("well-formed LM output should be accepted");

        assert_eq!(coeffs, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    #[ignore]
    fn test_validate_lm_fit_result_rejects_non_finite_or_weak_output() {
        let err =
            NISTdata::validate_lm_fit_result("unit-test", Some(vec![1.0, f64::NAN]), 2, Some(0.99))
                .expect_err("non-finite coefficients must be rejected");
        assert!(matches!(
            err,
            NISTError::FittingError(message) if message.contains("non-finite coefficient")
        ));

        let err = NISTdata::validate_lm_fit_result("unit-test", Some(vec![1.0, 2.0]), 2, Some(0.9))
            .expect_err("weak solver termination must be rejected");
        assert!(matches!(
            err,
            NISTError::FittingError(message) if message.contains("terminated too poorly")
        ));
    }
}
