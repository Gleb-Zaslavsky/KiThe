use crate::Kinetics::experimental_kinetics::kinetic_methods::kinetic_regression::{
    LinearRegressionResult, linear_regression,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::{
    KineticDataView, KineticMethod, KineticRequirements, TGADomainError, check_requirements,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::ColumnNature;
use ndarray::Array1;
#[derive(Debug, Clone)]
pub struct CKAResult {
    pub n: f64,
    pub m: f64,

    pub ea: f64,

    pub regression: LinearRegressionResult,
}
#[derive(Debug, Clone)]
pub struct CombinedKineticAnalysis {
    pub n_min: f64,
    pub n_max: f64,
    pub n_steps: usize,

    pub m_min: f64,
    pub m_max: f64,
    pub m_steps: usize,

    pub eta_min: f64,
    pub eta_max: f64,

    pub refinement_steps: usize,
}
impl Default for CombinedKineticAnalysis {
    fn default() -> Self {
        Self {
            n_min: -2.0,
            n_max: 5.0,
            n_steps: 40,

            m_min: -2.0,
            m_max: 5.0,
            m_steps: 40,

            eta_min: 0.02,
            eta_max: 0.98,
            refinement_steps: 10,
        }
    }
}
impl CombinedKineticAnalysis {
    pub fn solve(&self, data: &KineticDataView) -> Result<CKAResult, TGADomainError> {
        let r = 8.314462618;

        let (x, a, b, c) = self.collect_data(data);

        if x.len() < 20 {
            return Err(TGADomainError::InvalidOperation(
                "Not enough valid data points".into(),
            ));
        }

        let pc = precompute(&x, &a, &b, &c);

        let (n0, m0) = self.grid_search_nm(&pc);

        let (n, m) = self.refine_nm(&pc, n0, m0);

        let reg = final_regression(&x, &a, &b, &c, n, m);

        let ea = -reg.slope * r;

        Ok(CKAResult {
            n,
            m,
            ea,

            regression: reg,
        })
    }
    fn collect_data(&self, data: &KineticDataView) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut x = Vec::new();
        let mut a = Vec::new();
        let mut b = Vec::new();
        let mut c = Vec::new();

        for exp in &data.experiments {
            for i in 0..exp.conversion.len() {
                let eta = exp.conversion[i];
                let rate = exp.conversion_rate[i];

                if eta < self.eta_min || eta > self.eta_max {
                    continue;
                }

                if rate <= 0.0 {
                    continue;
                }

                let T = exp.temperature[i];

                x.push(1.0 / T);
                a.push(rate.ln());
                b.push((1.0 - eta).ln());
                c.push(eta.ln());
            }
        }

        (x, a, b, c)
    }

    fn grid_search_nm(&self, pc: &CKAPrecomputed) -> (f64, f64) {
        let mut best_n = 0.0;
        let mut best_m = 0.0;
        let mut best_corr = -1.0;

        let dn = (self.n_max - self.n_min) / self.n_steps as f64;
        let dm = (self.m_max - self.m_min) / self.m_steps as f64;

        for i in 0..=self.n_steps {
            let n = self.n_min + i as f64 * dn;

            for j in 0..=self.m_steps {
                let m = self.m_min + j as f64 * dm;

                let corr = correlation_nm(pc, n, m).abs();

                if corr > best_corr {
                    best_corr = corr;

                    best_n = n;
                    best_m = m;
                }
            }
        }

        (best_n, best_m)
    }
    fn refine_nm(&self, pc: &CKAPrecomputed, n0: f64, m0: f64) -> (f64, f64) {
        let step_n = (self.n_max - self.n_min) / self.n_steps as f64 / 10.0;

        let step_m = (self.m_max - self.m_min) / self.m_steps as f64 / 10.0;

        let mut best_n = n0;
        let mut best_m = m0;

        let mut best_corr = correlation_nm(pc, n0, m0).abs();

        for i in -(self.refinement_steps as i32)..=(self.refinement_steps as i32) {
            for j in -(self.refinement_steps as i32)..=(self.refinement_steps as i32) {
                let n = n0 + i as f64 * step_n;
                let m = m0 + j as f64 * step_m;

                let corr = correlation_nm(pc, n, m).abs();

                if corr > best_corr {
                    best_corr = corr;

                    best_n = n;
                    best_m = m;
                }
            }
        }

        (best_n, best_m)
    }
}

//==================================================================================================
// KineticMethod impl (wrapper)
//==================================================================================================
impl KineticMethod for CombinedKineticAnalysis {
    type Output = CKAResult;

    fn name(&self) -> &'static str {
        "Combined kinetic analysis"
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        self.solve(data)
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 1,
            needs_conversion: true,
            needs_conversion_rate: true,
            needs_temperature: false,
            needs_heating_rate: false,
        }
    }
    fn check_input(&self, data: &KineticDataView) -> Result<(), TGADomainError> {
        check_requirements(data, &self.requirements())
    }
    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        vec![
            ColumnNature::Conversion,
            ColumnNature::Temperature,
            ColumnNature::ConversionRate,
        ]
    }
}
#[derive(Debug, Clone)]
struct CKAPrecomputed {
    cov_xa: f64,
    cov_xb: f64,
    cov_xc: f64,

    var_x: f64,

    var_a: f64,
    var_b: f64,
    var_c: f64,

    cov_ab: f64,
    cov_ac: f64,
    cov_bc: f64,
}

fn correlation_nm(pc: &CKAPrecomputed, n: f64, m: f64) -> f64 {
    let cov_xy = pc.cov_xa - n * pc.cov_xb - m * pc.cov_xc;

    let var_y =
        pc.var_a + n * n * pc.var_b + m * m * pc.var_c - 2.0 * n * pc.cov_ab - 2.0 * m * pc.cov_ac
            + 2.0 * n * m * pc.cov_bc;

    if var_y <= 0.0 {
        return 0.0;
    }

    cov_xy / (pc.var_x * var_y).sqrt()
}

fn final_regression(
    x: &[f64],
    a: &[f64],
    b: &[f64],
    c: &[f64],
    n: f64,
    m: f64,
) -> LinearRegressionResult {
    let mut y = Vec::with_capacity(x.len());

    for i in 0..x.len() {
        y.push(a[i] - n * b[i] - m * c[i]);
    }

    linear_regression(&Array1::from_vec(x.to_vec()), &Array1::from_vec(y))
}
fn precompute(x: &[f64], a: &[f64], b: &[f64], c: &[f64]) -> CKAPrecomputed {
    CKAPrecomputed {
        cov_xa: covariance(x, a),
        cov_xb: covariance(x, b),
        cov_xc: covariance(x, c),

        var_x: covariance(x, x),

        var_a: covariance(a, a),
        var_b: covariance(b, b),
        var_c: covariance(c, c),

        cov_ab: covariance(a, b),
        cov_ac: covariance(a, c),
        cov_bc: covariance(b, c),
    }
}
fn covariance(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;

    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;

    let mut s = 0.0;

    for (&xi, &yi) in x.iter().zip(y) {
        s += (xi - mx) * (yi - my);
    }

    s / n
}

//========================================================================================================
// TESTS
//=================================================================================================
#[cfg(test)]
mod tests_differential {
    use super::*;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion_tests::tests::simulate_tga_first_order2;
    use crate::Kinetics::experimental_kinetics::kinetic_methods_tests::tests::build_view_from_cfg_exact_m0;
    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::base_advanced_config_non_isothermal;
    use std::time::Instant;

    #[test]
    fn friedman_differential_compute_with_mock_non_isothermal_data() {
        let cfg = base_advanced_config_non_isothermal(
            700.0,
            1e5,
            80_000.0,
            0.1,
            10_000,
            vec![0.5, 3.0, 5.0],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let result = CombinedKineticAnalysis::default().solve(&view).unwrap();
        println!("{:?}", result);
    }

    #[test]
    fn friedman_differential_compute_with_mock_non_isothermal_data2() {
        let now = Instant::now();
        let cfg = base_advanced_config_non_isothermal(
            420.0,
            1e6,
            100_000.0,
            1.0,
            30_000,
            vec![0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let result = CombinedKineticAnalysis::default().solve(&view).unwrap();
        println!("{:?}", result);
    }

    #[test]
    fn friedman_differential_compute_with_mock_non_isothermal_data3() {
        let cfg = base_advanced_config_non_isothermal(
            600.0,
            1e7,
            150_000.0,
            0.5,
            50_000,
            vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let result = CombinedKineticAnalysis::default().solve(&view).unwrap();
        println!("{:?}", result);
    }
}
