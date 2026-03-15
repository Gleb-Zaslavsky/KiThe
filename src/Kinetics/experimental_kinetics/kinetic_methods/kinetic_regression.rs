use ndarray::Array1;
use ndarray_stats::CorrelationExt;
// Структура результата регрессии
#[derive(Clone, Debug)]
pub struct LinearRegressionResult {
    pub slope: f64,
    pub intercept: f64,

    pub r: f64,
    pub r2: f64,

    pub stderr: f64,
}
impl Default for LinearRegressionResult {
    fn default() -> Self {
        LinearRegressionResult {
            slope: 0.0,
            intercept: 0.0,
            r: 0.0,
            r2: 0.0,
            stderr: 0.0,
        }
    }
}
// Линейная регрессия (OLS)

//Используем ndarray.

pub fn linear_regression(x: &Array1<f64>, y: &Array1<f64>) -> LinearRegressionResult {
    let n = x.len() as f64;

    let x_mean = x.mean().unwrap();
    let y_mean = y.mean().unwrap();

    let dx = x - x_mean;
    let dy = y - y_mean;

    let sxy = dx.dot(&dy);
    let sxx = dx.dot(&dx);

    let slope = sxy / sxx;

    let intercept = y_mean - slope * x_mean;

    let r = sxy / (sxx.sqrt() * dy.dot(&dy).sqrt());

    let r2 = r * r;

    let stderr = ((dy.dot(&dy) - slope * sxy) / (n - 2.0)).sqrt();

    LinearRegressionResult {
        slope,
        intercept,
        r,
        r2,
        stderr,
    }
}

// Pearson correlation (через ndarray_stats)

pub fn pearson(x: &Array1<f64>, y: &Array1<f64>) -> f64 {
    let mut m = ndarray::stack(ndarray::Axis(0), &[x.view(), y.view()]).unwrap();

    let corr = m.pearson_correlation().unwrap();

    corr[[0, 1]]
}

// Covariance
pub fn covariance(x: &Array1<f64>, y: &Array1<f64>) -> f64 {
    let xm = x.mean().unwrap();
    let ym = y.mean().unwrap();

    let dx = x - xm;
    let dy = y - ym;

    dx.dot(&dy) / (x.len() as f64 - 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array1;

    #[test]
    fn linear_regression_exact_line() {
        let x = Array1::from_vec(vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        let y = Array1::from_vec(vec![1.0, 3.0, 5.0, 7.0, 9.0]);

        let out = linear_regression(&x, &y);

        assert_relative_eq!(out.slope, 2.0, epsilon = 1e-12);
        assert_relative_eq!(out.intercept, 1.0, epsilon = 1e-12);
        assert_relative_eq!(out.r, 1.0, epsilon = 1e-12);
        assert_relative_eq!(out.r2, 1.0, epsilon = 1e-12);
        assert_relative_eq!(out.stderr, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn pearson_matches_perfect_negative_correlation() {
        let x = Array1::from_vec(vec![0.0, 1.0, 2.0, 3.0]);
        let y = Array1::from_vec(vec![0.0, -1.0, -2.0, -3.0]);

        let r = pearson(&x, &y);

        assert_relative_eq!(r, -1.0, epsilon = 1e-12);
    }

    #[test]
    fn covariance_matches_expected() {
        let x = Array1::from_vec(vec![1.0, 2.0, 3.0]);
        let y = Array1::from_vec(vec![1.0, 5.0, 7.0]);

        let cov = covariance(&x, &y);
        assert_relative_eq!(cov, 3.0, epsilon = 1e-12);

        let cov_rev = covariance(&y, &x);
        assert_relative_eq!(cov_rev, 3.0, epsilon = 1e-12);
    }
}
