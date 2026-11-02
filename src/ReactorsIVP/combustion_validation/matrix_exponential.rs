use nalgebra::{DMatrix, DVector};

const PADE13_THETA: f64 = 5.371_920_351_148_152;
const PADE13_COEFFICIENTS: [f64; 14] = [
    64_764_752_532_480_000.0,
    32_382_376_266_240_000.0,
    7_771_770_303_897_600.0,
    1_187_353_796_428_800.0,
    129_060_195_264_000.0,
    10_559_470_521_600.0,
    670_442_572_800.0,
    33_522_128_640.0,
    1_323_241_920.0,
    40_840_800.0,
    960_960.0,
    16_380.0,
    182.0,
    1.0,
];

pub(super) fn expm_action(matrix: &DMatrix<f64>, vector: &DVector<f64>) -> DVector<f64> {
    assert_eq!(
        matrix.nrows(),
        matrix.ncols(),
        "expm requires a square matrix"
    );
    assert_eq!(
        matrix.ncols(),
        vector.len(),
        "matrix/vector dimensions differ"
    );

    let norm = matrix_one_norm(matrix);
    let squarings = if norm <= PADE13_THETA || norm == 0.0 {
        0
    } else {
        (norm / PADE13_THETA).log2().ceil().max(0.0) as u32
    };
    let scaled = matrix / 2.0_f64.powi(squarings as i32);
    let mut exponential = pade13(&scaled);
    for _ in 0..squarings {
        exponential = &exponential * &exponential;
    }
    exponential * vector
}

fn matrix_one_norm(matrix: &DMatrix<f64>) -> f64 {
    (0..matrix.ncols())
        .map(|column| matrix.column(column).iter().map(|value| value.abs()).sum())
        .fold(0.0, f64::max)
}

fn pade13(matrix: &DMatrix<f64>) -> DMatrix<f64> {
    let n = matrix.nrows();
    let identity = DMatrix::identity(n, n);
    let a2 = matrix * matrix;
    let a4 = &a2 * &a2;
    let a6 = &a4 * &a2;
    let b = &PADE13_COEFFICIENTS;

    let u_inner = &a6 * (b[13] * &a6 + b[11] * &a4 + b[9] * &a2)
        + b[7] * &a6
        + b[5] * &a4
        + b[3] * &a2
        + b[1] * &identity;
    let u = matrix * u_inner;
    let v = &a6 * (b[12] * &a6 + b[10] * &a4 + b[8] * &a2)
        + b[6] * &a6
        + b[4] * &a4
        + b[2] * &a2
        + b[0] * identity;
    let denominator = &v - &u;
    denominator
        .lu()
        .solve(&(&v + &u))
        .expect("Pade-13 denominator should be nonsingular for finite matrices")
}

#[cfg(test)]
mod tests {
    use super::expm_action;
    use nalgebra::{DMatrix, DVector};

    fn assert_vector_close(actual: &DVector<f64>, expected: &DVector<f64>, tolerance: f64) {
        assert_eq!(actual.len(), expected.len());
        for index in 0..actual.len() {
            assert!(
                (actual[index] - expected[index]).abs() <= tolerance,
                "index={index} actual={} expected={} error={}",
                actual[index],
                expected[index],
                (actual[index] - expected[index]).abs()
            );
        }
    }

    #[test]
    fn expm_action_passes_independent_closed_form_checks() {
        let vector = DVector::from_vec(vec![1.0, 2.0]);
        assert_vector_close(&expm_action(&DMatrix::zeros(2, 2), &vector), &vector, 1e-14);

        let diagonal = DMatrix::from_diagonal(&DVector::from_vec(vec![0.5, -0.25]));
        let diagonal_result = expm_action(&diagonal, &vector);
        assert_vector_close(
            &diagonal_result,
            &DVector::from_vec(vec![0.5_f64.exp(), 2.0 * (-0.25_f64).exp()]),
            2e-14,
        );

        let nilpotent = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 0.0, 0.0]);
        assert_vector_close(
            &expm_action(&nilpotent, &vector),
            &DVector::from_vec(vec![3.0, 2.0]),
            2e-14,
        );

        let rotation = DMatrix::from_row_slice(2, 2, &[0.0, -1.0, 1.0, 0.0]);
        assert_vector_close(
            &expm_action(&rotation, &vector),
            &DVector::from_vec(vec![
                1.0_f64.cos() - 2.0 * 1.0_f64.sin(),
                1.0_f64.sin() + 2.0 * 1.0_f64.cos(),
            ]),
            2e-14,
        );

        let a = 0.4 * &rotation;
        let b = 0.7 * &rotation;
        let initial = DVector::from_vec(vec![1.0, 0.0]);
        let sequential = expm_action(&a, &expm_action(&b, &initial));
        let combined = expm_action(&(a + b), &initial);
        assert_vector_close(&sequential, &combined, 3e-14);
    }
}
