use crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError;
use ndarray::Array1;

/// Summary of an ordinary least-squares fit to paired numeric samples.
#[derive(Clone, Debug)]
pub struct LinearRegressionResult {
    pub slope: f64,
    pub intercept: f64,
    pub r: f64,
    pub r2: f64,
    /// Residual standard error is unavailable when there are only two points.
    pub stderr: Option<f64>,
}

impl Default for LinearRegressionResult {
    fn default() -> Self {
        Self {
            slope: 0.0,
            intercept: 0.0,
            r: 0.0,
            r2: 0.0,
            stderr: None,
        }
    }
}

/// Fits `y = slope * x + intercept` and rejects invalid or degenerate inputs.
pub fn linear_regression(
    x: &Array1<f64>,
    y: &Array1<f64>,
) -> Result<LinearRegressionResult, TGADomainError> {
    validate_pair(x, y, "linear regression")?;
    let n = x.len() as f64;
    let x_mean = x.sum() / n;
    let y_mean = y.sum() / n;
    let dx = x - x_mean;
    let dy = y - y_mean;
    let sxy = dx.dot(&dy);
    let sxx = dx.dot(&dx);
    let syy = dy.dot(&dy);

    if !sxx.is_finite() || sxx <= 0.0 || !syy.is_finite() || syy <= 0.0 {
        return Err(TGADomainError::InvalidOperation(
            "Linear regression requires non-constant finite X and Y values".into(),
        ));
    }

    let slope = sxy / sxx;
    let intercept = y_mean - slope * x_mean;
    let r = (sxy / (sxx.sqrt() * syy.sqrt())).clamp(-1.0, 1.0);
    let r2 = r * r;
    let stderr = (n > 2.0).then(|| ((syy - slope * sxy).max(0.0) / (n - 2.0)).sqrt());

    if !slope.is_finite() || !intercept.is_finite() || !r.is_finite() {
        return Err(TGADomainError::InvalidOperation(
            "Linear regression produced non-finite coefficients".into(),
        ));
    }

    Ok(LinearRegressionResult {
        slope,
        intercept,
        r,
        r2,
        stderr,
    })
}

/// Computes Pearson correlation for equal-length, non-constant finite samples.
pub fn pearson(x: &Array1<f64>, y: &Array1<f64>) -> Result<f64, TGADomainError> {
    validate_pair(x, y, "Pearson correlation")?;
    let n = x.len() as f64;
    let dx = x - (x.sum() / n);
    let dy = y - (y.sum() / n);
    let sxx = dx.dot(&dx);
    let syy = dy.dot(&dy);
    if !sxx.is_finite() || !syy.is_finite() || sxx <= 0.0 || syy <= 0.0 {
        return Err(TGADomainError::InvalidOperation(
            "Pearson correlation requires finite, non-constant data".into(),
        ));
    }
    let correlation = dx.dot(&dy) / (sxx.sqrt() * syy.sqrt());
    if !correlation.is_finite() {
        return Err(TGADomainError::InvalidOperation(
            "Pearson correlation produced a non-finite value".into(),
        ));
    }
    Ok(correlation.clamp(-1.0, 1.0))
}

/// Computes the unbiased sample covariance of equal-length finite samples.
pub fn covariance(x: &Array1<f64>, y: &Array1<f64>) -> Result<f64, TGADomainError> {
    validate_pair(x, y, "covariance")?;
    let n = x.len() as f64;
    let dx = x - (x.sum() / n);
    let dy = y - (y.sum() / n);
    let result = dx.dot(&dy) / (n - 1.0);
    if !result.is_finite() {
        return Err(TGADomainError::InvalidOperation(
            "Covariance produced a non-finite value".into(),
        ));
    }
    Ok(result)
}

fn validate_pair(x: &Array1<f64>, y: &Array1<f64>, operation: &str) -> Result<(), TGADomainError> {
    if x.len() != y.len() {
        return Err(TGADomainError::InvalidOperation(format!(
            "{operation} requires equal input lengths (X={}, Y={})",
            x.len(),
            y.len()
        )));
    }
    if x.len() < 2 {
        return Err(TGADomainError::InvalidOperation(format!(
            "{operation} requires at least two paired values"
        )));
    }
    if x.iter().chain(y.iter()).any(|value| !value.is_finite()) {
        return Err(TGADomainError::InvalidOperation(format!(
            "{operation} inputs must be finite"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn linear_regression_exact_line() {
        let x = Array1::from_vec(vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        let y = Array1::from_vec(vec![1.0, 3.0, 5.0, 7.0, 9.0]);
        let out = linear_regression(&x, &y).unwrap();
        assert_relative_eq!(out.slope, 2.0, epsilon = 1e-12);
        assert_relative_eq!(out.intercept, 1.0, epsilon = 1e-12);
        assert_relative_eq!(out.r, 1.0, epsilon = 1e-12);
        assert_relative_eq!(out.r2, 1.0, epsilon = 1e-12);
        assert_eq!(out.stderr, Some(0.0));
    }

    #[test]
    fn pearson_matches_perfect_negative_correlation() {
        let x = Array1::from_vec(vec![0.0, 1.0, 2.0, 3.0]);
        let y = Array1::from_vec(vec![0.0, -1.0, -2.0, -3.0]);
        assert_relative_eq!(pearson(&x, &y).unwrap(), -1.0, epsilon = 1e-12);
    }

    #[test]
    fn covariance_matches_expected() {
        let x = Array1::from_vec(vec![1.0, 2.0, 3.0]);
        let y = Array1::from_vec(vec![1.0, 5.0, 7.0]);
        assert_relative_eq!(covariance(&x, &y).unwrap(), 3.0, epsilon = 1e-12);
        assert_relative_eq!(covariance(&y, &x).unwrap(), 3.0, epsilon = 1e-12);
    }

    #[test]
    fn regression_helpers_reject_invalid_inputs() {
        let constant = Array1::from_vec(vec![1.0, 1.0]);
        let varying = Array1::from_vec(vec![2.0, 3.0]);
        assert!(linear_regression(&constant, &varying).is_err());
        assert!(pearson(&constant, &varying).is_err());

        let x = Array1::from_vec(vec![1.0, 2.0]);
        let short_y = Array1::from_vec(vec![2.0]);
        assert!(covariance(&x, &short_y).is_err());
    }

    #[test]
    fn two_point_fit_marks_standard_error_unavailable() {
        let x = Array1::from_vec(vec![1.0, 2.0]);
        let y = Array1::from_vec(vec![2.0, 4.0]);
        assert_eq!(linear_regression(&x, &y).unwrap().stderr, None);
    }
}
