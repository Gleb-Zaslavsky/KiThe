use RustedSciThe::symbolic::symbolic_engine::Expr;
use thiserror::Error;

/// Shared comparison metrics for property-wise fitting validation.
///
/// The helper keeps the contract stable across NASA and NIST fitting paths:
/// compare raw property samples against the published symbolic expression,
/// reject malformed inputs early, and report relative-error diagnostics in one
/// place.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PropertyComparisonReport {
    pub(crate) l2_norm: f64,
    pub(crate) max_norm: f64,
    pub(crate) number_of_significant_points: usize,
    pub(crate) t_range: Option<(f64, f64)>,
}

/// Errors produced while validating or comparing property series.
#[derive(Debug, Error, PartialEq)]
pub(crate) enum PropertyComparisonError {
    #[error("comparison vectors must have the same length: y_data={y_data}, T_vec={t_vec}")]
    LengthMismatch { y_data: usize, t_vec: usize },
    #[error("comparison vectors must not be empty")]
    EmptyInput,
    #[error("comparison vectors must contain only finite values")]
    NonFiniteInput,
    #[error("symbolic comparison produced non-finite fitted values")]
    NonFiniteOutput,
}

const SIGNIFICANT_RESIDUAL_THRESHOLD: f64 = 0.10;
const ZERO_REFERENCE_FLOOR: f64 = 1e-12;

/// Compare raw property samples against a symbolic fit expression.
///
/// The contract is intentionally strict:
/// - equal-length, non-empty, finite input vectors only;
/// - relative error is measured against a small reference floor to avoid
///   division by zero when the raw property is near zero;
/// - the returned range only covers points that cross the significant-residual
///   threshold.
pub(crate) fn compare_property_series(
    y_data: &[f64],
    t_vec: &[f64],
    fit_function: Expr,
) -> Result<PropertyComparisonReport, PropertyComparisonError> {
    if y_data.len() != t_vec.len() {
        return Err(PropertyComparisonError::LengthMismatch {
            y_data: y_data.len(),
            t_vec: t_vec.len(),
        });
    }
    if y_data.is_empty() {
        return Err(PropertyComparisonError::EmptyInput);
    }
    if y_data.iter().any(|value| !value.is_finite()) || t_vec.iter().any(|value| !value.is_finite())
    {
        return Err(PropertyComparisonError::NonFiniteInput);
    }

    let fitted = fit_function.lambdify1D();
    let fitted_values = t_vec.iter().map(|&t| fitted(t)).collect::<Vec<f64>>();
    if fitted_values.iter().any(|value| !value.is_finite()) {
        return Err(PropertyComparisonError::NonFiniteOutput);
    }

    let mut squared_relative_sum = 0.0;
    let mut max_norm: f64 = 0.0;
    let mut significant_points: Vec<(f64, f64)> = Vec::new();

    for ((&reference, &temperature), &fitted_value) in
        y_data.iter().zip(t_vec.iter()).zip(fitted_values.iter())
    {
        let abs_error = (reference - fitted_value).abs();
        let reference_scale = reference.abs().max(ZERO_REFERENCE_FLOOR);
        let relative_error = abs_error / reference_scale;

        squared_relative_sum += relative_error * relative_error;
        max_norm = max_norm.max(relative_error);

        if relative_error > SIGNIFICANT_RESIDUAL_THRESHOLD {
            significant_points.push((temperature, relative_error));
        }
    }

    let l2_norm = (squared_relative_sum / y_data.len() as f64).sqrt();
    let t_range = if significant_points.is_empty() {
        None
    } else {
        let min_t = significant_points
            .iter()
            .map(|(temperature, _)| *temperature)
            .fold(f64::INFINITY, f64::min);
        let max_t = significant_points
            .iter()
            .map(|(temperature, _)| *temperature)
            .fold(f64::NEG_INFINITY, f64::max);
        Some((min_t, max_t))
    };

    Ok(PropertyComparisonReport {
        l2_norm,
        max_norm,
        number_of_significant_points: significant_points.len(),
        t_range,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compare_property_series_rejects_empty_input() {
        let err = compare_property_series(&[], &[], Expr::Const(0.0)).unwrap_err();
        assert!(matches!(err, PropertyComparisonError::EmptyInput));
    }

    #[test]
    fn compare_property_series_rejects_length_mismatch() {
        let err = compare_property_series(&[1.0], &[1.0, 2.0], Expr::Const(1.0)).unwrap_err();
        assert!(matches!(
            err,
            PropertyComparisonError::LengthMismatch { .. }
        ));
    }

    #[test]
    fn compare_property_series_handles_zero_references() {
        let report = compare_property_series(&[0.0, 0.0], &[300.0, 400.0], Expr::Const(0.0))
            .expect("zero-valued references should be valid");
        assert!(report.l2_norm.is_finite());
        assert!(report.max_norm.is_finite());
        assert_eq!(report.number_of_significant_points, 0);
    }
}
