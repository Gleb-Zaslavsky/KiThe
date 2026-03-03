use ndarray::{Array1, Array2};
use ndarray_linalg::QR;
use ndarray_linalg::Solve;
use thiserror::Error;
// use ndarray_linalg::QRInto;
// bspline.rs


#[derive(Clone, Debug)]
pub struct BSpline {
    pub knots: Array1<f64>,
    pub coeffs: Array1<f64>,
    pub degree: usize,
}
impl BSpline {
    pub fn evaluate(&self, x: f64) -> f64 {
            deboor(
            x,
            &self.knots,
            &self.coeffs,
            self.degree,
        )
    }
}
//=================================================================================
//          Basis + de Boor (production safe)
//======================================================================================
pub fn find_span(n: usize, k: usize, x: f64, t: &Array1<f64>) -> usize {
    if x == t[n + 1] {
        return n;
    }

    let mut low = k;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;

    while x < t[mid] || x >= t[mid + 1] {
        if x < t[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }

    mid
}
pub fn basis_functions(
    span: usize,
    x: f64,
    k: usize,
    t: &Array1<f64>,
) -> Vec<f64> {
    let mut left = vec![0.0; k + 1];
    let mut right = vec![0.0; k + 1];
    let mut N = vec![0.0; k + 1];

    N[0] = 1.0;

    for j in 1..=k {
        left[j] = x - t[span + 1 - j];
        right[j] = t[span + j] - x;

        let mut saved = 0.0;

        for r in 0..j {
            let temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }

        N[j] = saved;
    }

    N
}
pub fn deboor(
    x: f64,
    t: &Array1<f64>,
    c: &Array1<f64>,
    k: usize,
) -> f64 {
    let n = c.len() - 1;
    let span = find_span(n, k, x, t);

    let mut d = vec![0.0; k + 1];

    for j in 0..=k {
        d[j] = c[span - k + j];
    }

    for r in 1..=k {
        for j in (r..=k).rev() {
            let i = span - k + j;
            let alpha = (x - t[i]) / (t[i + k + 1 - r] - t[i]);
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
        }
    }

    d[k]
}
//=======================================================================================
//             Design Matrix (Dense)
//=======================================================================================
pub fn design_matrix_dense(
    x: &Array1<f64>,
    knots: &Array1<f64>,
    degree: usize,
) -> Array2<f64> {
    let m = x.len();
    let n = knots.len() - degree - 1;

    let mut A = Array2::<f64>::zeros((m, n));

    for (i, &xi) in x.iter().enumerate() {
        let span = find_span(n - 1, degree, xi, knots);
        let basis = basis_functions(span, xi, degree, knots);

        for j in 0..=degree {
            A[[i, span - degree + j]] = basis[j];
        }
    }

    A
}
//============================================================================
// Solver enum (production API)
//==============================================================================

#[derive(Debug, Error)]
pub enum LssError {
    #[error("Linear algebra failure")]
    Linalg,
    #[error("Invalid spline input")]
InvalidInput,
}

#[derive(Clone, Copy)]
pub enum SolverKind {
    DenseQR,
    Banded, // future
}
pub fn solve_lsq(
    A: Array2<f64>,
    y: &Array1<f64>,
    w: &Array1<f64>,
    solver: SolverKind,
) -> Result<Array1<f64>, LssError> {
    match solver {
        SolverKind::DenseQR => dense_normal_eq(A, y, w),
        SolverKind::Banded => {
            Err(LssError::Linalg) // placeholder
        }
    }
}
fn dense_qr(
    mut A: Array2<f64>,
    y: &Array1<f64>,
    w: &Array1<f64>,
) -> Result<Array1<f64>, LssError> {

    let mut yw = y.clone();

    // apply weights
    for i in 0..w.len() {
        let wi = w[i];
        A.row_mut(i).mapv_inplace(|v| v * wi);
        yw[i] *= wi;
    }

    let (Q, R) = A.qr().map_err(|_| LssError::Linalg)?;

    let qt_y = Q.t().dot(&yw);

    let coeffs = R.solve_into(qt_y)
        .map_err(|_| LssError::Linalg)?;

    Ok(coeffs)
}
//1️⃣ Проверка Schoenberg–Whitney

fn dense_normal_eq(
    mut A: Array2<f64>,
    y: &Array1<f64>,
    w: &Array1<f64>,
) -> Result<Array1<f64>, LssError> {

    let mut yw = y.clone();

    for i in 0..w.len() {
        let wi = w[i];
        A.row_mut(i).mapv_inplace(|v| v * wi);
        yw[i] *= wi;
    }

    let ata = A.t().dot(&A);
    let aty = A.t().dot(&yw);

    let coeffs = ata.solve_into(aty)
        .map_err(|_| LssError::Linalg)?;

    Ok(coeffs)
}

//Условие:
//Для j = 0..(n-k-2):
//t[j] < x[i_j] < t[j+k+1]
//Мы реализуем проверку существования хотя бы одного x в каждом интервале.

pub fn check_schoenberg_whitney(
    x: &Array1<f64>,
    knots: &Array1<f64>,
    degree: usize,
) -> Result<(), &'static str> {

    let n = knots.len() - degree - 1;

    for j in 0..(n - degree - 1) {
        let left = knots[j];
        let right = knots[j + degree + 1];

        let exists = x.iter().any(|&xi| xi > left && xi < right);

        if !exists {
            return Err("Schoenberg–Whitney condition violated");
        }
    }

    Ok(())
}
//==============================================================================
//                High-level API
//====================================================================================

pub fn make_lsq_spline(
    x: &Array1<f64>,
    y: &Array1<f64>,
    knots: &Array1<f64>,
    degree: usize,
    weights: Option<&Array1<f64>>,
    solver: SolverKind,
) -> Result<BSpline, LssError> {

    let w = weights
        .cloned()
        .unwrap_or_else(|| Array1::ones(x.len()));

    let A = design_matrix_dense(x, knots, degree);

    let coeffs = solve_lsq(A, y, &w, solver)?;

    Ok(BSpline {
        knots: knots.clone(),
        coeffs,
        degree,
    })
}


pub fn make_lsq_univariate_spline(
    x: &Array1<f64>,
    y: &Array1<f64>,
    internal_knots: &Array1<f64>,
    degree: usize,
    solver: SolverKind,
) -> Result<BSpline, LssError> {

    let xmin = x[0];
    let xmax = x[x.len() - 1];

    // build full knot vector
    let mut knots = Vec::new();

    for _ in 0..=degree {
        knots.push(xmin);
    }

    for &k in internal_knots {
        knots.push(k);
    }

    for _ in 0..=degree {
        knots.push(xmax);
    }

    let knots = Array1::from(knots);

    check_schoenberg_whitney(x, &knots, degree)
        .map_err(|_| LssError::InvalidInput)?;

    make_lsq_spline(
        x,
        y,
        &knots,
        degree,
        None,
        solver,
    )
}
//=================================================================================================
// TESTING
//=================================================================================================
#[cfg(test)]
mod tests {
    use ndarray::ArrayBase;

    use super::*;
 //   ✅ Тест 1: Аппроксимация exp(-x²)
#[test]
fn test_lsq_basic_fit() {
    use ndarray::Array1;

    let x:ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>, f64> = Array1::linspace(-3., 3., 200);
    let y = x.mapv(|v | (-v * v).exp());

    let internal = Array1::linspace(-2., 2., 10);

    let spline = make_lsq_univariate_spline(
        &x,
        &y,
        &internal,
        3,
        SolverKind::DenseQR,
    ).unwrap();

    let xs = Array1::linspace(-3., 3., 100);

    let ys = xs.mapv(|v| spline.evaluate(v));

    let mse = (&ys - &xs.mapv(|v| (-v * v).exp()))
        .mapv(|v| v * v)
        .mean()
        .unwrap();

    assert!(mse < 1e-3);
}

//✅ Тест 2: Проверка Schoenberg–Whitney
#[test]
fn test_sw_violation() {
    use ndarray::array;

    let x = array![0., 1., 2., 3.];
    let y = x.clone();

    let internal = array![10.]; // outside

    let result = make_lsq_univariate_spline(
        &x,
        &y,
        &internal,
        3,
        SolverKind::DenseQR,
    );

    assert!(result.is_err());
}
   


    
}