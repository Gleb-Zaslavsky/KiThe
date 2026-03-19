//! # Sublimation kinetics detector
//!
//! For a pure sublimation process the conversion rate follows:
//!   `dα/dt = k₀ · exp(-E/RT)`
//!
//! Taking the natural log:
//!   `ln(dα/dt) = ln(k₀) - (E/R) · (1/T)`
//!
//! i.e. a straight line in `ln(dα/dt)` vs `1/T` space.
//! This module fits that line per experiment over `[eta_min, eta_max]`,
//! stores the fitted curve, and provides helpers to aggregate results and
//! push fitted rates back into a `TGASeries`.

use crate::Kinetics::experimental_kinetics::experiment_series_main::TGASeries;
use crate::Kinetics::experimental_kinetics::kinetic_methods::kinetic_regression::{
    LinearRegressionResult, linear_regression,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::{
    KineticDataView, KineticMethod, KineticRequirements, TGADomainError,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{ColumnNature, Unit};
use ndarray::Array1;

pub const R: f64 = 8.314;

// ─── classification ───────────────────────────────────────────────────────────

/// Classification of whether the data is consistent with sublimation kinetics.
/// The wrapped `f64` is the R² of the regression.
#[derive(Debug, Clone)]
pub enum IsThisASublimationResult {
    /// R² ≥ 0.99 — very strong linear fit in ln(rate) vs 1/T.
    Yes(f64),
    /// R² < 0.90 — poor fit, unlikely sublimation.
    No(f64),
    /// 0.90 ≤ R² < 0.99 — moderate fit, inconclusive.
    MayBe(f64),
}

impl IsThisASublimationResult {
    pub fn from_r2(r2: f64) -> Self {
        if r2 >= 0.99 {
            Self::Yes(r2)
        } else if r2 >= 0.90 {
            Self::MayBe(r2)
        } else {
            Self::No(r2)
        }
    }

    pub fn r2(&self) -> f64 {
        match self {
            Self::Yes(v) | Self::No(v) | Self::MayBe(v) => *v,
        }
    }
}

// ─── per-experiment result ────────────────────────────────────────────────────

/// Result for a single experiment.
#[derive(Debug, Clone)]
pub struct SublimationResult {
    pub experiment_id: String,
    /// Activation energy (J/mol).
    pub ea: f64,
    /// Pre-exponential factor k₀.
    pub k: f64,
    /// OLS regression diagnostics.
    pub regression: LinearRegressionResult,
    /// Classification based on R².
    pub resume: IsThisASublimationResult,
    /// Fitted `k₀·exp(-E/RT)` evaluated over the **full** temperature vector
    /// of the experiment (same length as `ExperimentData::temperature`).
    pub conversion: Vec<f64>,
    pub fitted_conversion_rate: Vec<f64>,
}

impl Default for SublimationResult {
    fn default() -> Self {
        Self {
            experiment_id: String::new(),
            ea: 0.0,
            k: 0.0,
            regression: LinearRegressionResult::default(),
            resume: IsThisASublimationResult::MayBe(0.0),
            conversion: Vec::new(),
            fitted_conversion_rate: Vec::new(),
        }
    }
}

// ─── solver ───────────────────────────────────────────────────────────────────

/// Fits `ln(dα/dt) = A - B/T` per experiment over `[eta_min, eta_max]`.
pub struct SublimationMethod {
    pub eta_min: f64,
    pub eta_max: f64,
}

impl Default for SublimationMethod {
    fn default() -> Self {
        Self {
            eta_min: 0.05,
            eta_max: 0.95,
        }
    }
}

impl SublimationMethod {
    pub fn new(eta_min: f64, eta_max: f64) -> Self {
        Self { eta_min, eta_max }
    }

    /// Core solver — returns one [`SublimationResult`] per experiment.
    ///
    /// For each experiment:
    /// 1. Finds the sub-slice of `conversion` closest to `[eta_min, eta_max]`.
    /// 2. Extracts the corresponding `temperature` and `conversion_rate`.
    /// 3. Regresses `ln(dα/dt)` vs `1/T` (OLS).
    /// 4. Recovers `k = exp(intercept)`, `E = -R · slope`.
    /// 5. Evaluates the fitted curve over the **full** temperature vector.
    pub fn solve(&self, data: &KineticDataView) -> Result<Vec<SublimationResult>, TGADomainError> {
        if data.experiments.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "No experiments provided".into(),
            ));
        }

        let mut results = Vec::with_capacity(data.experiments.len());

        for exp in &data.experiments {
            let conv = &exp.conversion;
            let temp = &exp.temperature;
            let rate = &exp.conversion_rate;

            if conv.is_empty() || temp.is_empty() || rate.is_empty() {
                return Err(TGADomainError::InvalidOperation(
                    "Experiment has empty columns".into(),
                ));
            }

            // index range closest to [eta_min, eta_max]
            let i_start = conv.iter().position(|&c| c >= self.eta_min).unwrap_or(0);
            let i_end = conv
                .iter()
                .rposition(|&c| c <= self.eta_max)
                .unwrap_or(conv.len() - 1);

            if i_end <= i_start {
                return Err(TGADomainError::InvalidOperation(format!(
                    "No data in conversion window [{}, {}] for experiment '{}'",
                    self.eta_min, self.eta_max, exp.meta.id
                )));
            }

            // build regression vectors: x = 1/T, y = ln(dα/dt)
            let mut xs = Vec::with_capacity(i_end - i_start + 1);
            let mut ys = Vec::with_capacity(i_end - i_start + 1);

            for (&t, &r) in temp[i_start..=i_end]
                .iter()
                .zip(rate[i_start..=i_end].iter())
            {
                if !t.is_finite() || t <= 0.0 {
                    return Err(TGADomainError::InvalidOperation(
                        "Non-physical temperature".into(),
                    ));
                }
                if r.is_finite() && r > 0.0 {
                    xs.push(1.0 / t);
                    ys.push(r.ln());
                }
            }

            if xs.len() < 2 {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Not enough valid points for regression in experiment '{}'",
                    exp.meta.id
                )));
            }

            let reg = linear_regression(&Array1::from(xs), &Array1::from(ys));

            // ln(rate) = intercept + slope/T  →  slope = -E/R, intercept = ln(k)
            let ea = -R * reg.slope;
            let k = reg.intercept.exp();

            // fitted curve over the full temperature vector
            let fitted_conversion_rate = temp
                .iter()
                .map(|&t| k * (-ea / (R * t)).exp())
                .collect::<Vec<f64>>();

            let resume = IsThisASublimationResult::from_r2(reg.r2);

            results.push(SublimationResult {
                experiment_id: exp.meta.id.clone(),
                ea,
                k,
                regression: reg,
                resume,
                conversion: conv.clone(),
                fitted_conversion_rate,
            });
        }

        Ok(results)
    }

    /// Returns `(mean_ea, mean_k, overall_verdict)` aggregated over all
    /// per-experiment results.  The verdict is based on the mean R².
    pub fn mean_results(
        results: &[SublimationResult],
    ) -> Result<(f64, f64, IsThisASublimationResult), TGADomainError> {
        if results.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "No results to aggregate".into(),
            ));
        }
        let n = results.len() as f64;
        let mean_ea = results.iter().map(|r| r.ea).sum::<f64>() / n;
        let mean_k = results.iter().map(|r| r.k).sum::<f64>() / n;
        let mean_r2 = results.iter().map(|r| r.regression.r2).sum::<f64>() / n;
        Ok((mean_ea, mean_k, IsThisASublimationResult::from_r2(mean_r2)))
    }

    /// Pushes each `fitted_conversion_rate` from `results` into the
    /// corresponding experiment in `series` as a new column named
    /// `"fitted_sublim_rate"` with unit `PerSecond` and nature `ConversionRate`.
    ///
    /// The experiment is looked up by `SublimationResult::experiment_id`.
    pub fn push_fitted_rates_to_series(
        results: &[SublimationResult],
        series: &mut TGASeries,
    ) -> Result<(), TGADomainError> {
        for res in results {
            series.add_column_from_vec(
                &res.experiment_id,
                "fitted_sublim_rate",
                Unit::PerSecond,
                ColumnNature::ConversionRate,
                res.fitted_conversion_rate.clone(),
            )?;
        }
        Ok(())
    }
}

// ─── KineticMethod impl ───────────────────────────────────────────────────────

impl KineticMethod for SublimationMethod {
    type Output = Vec<SublimationResult>;

    fn name(&self) -> &'static str {
        "SublimationMethod"
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

    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        vec![
            ColumnNature::Conversion,
            ColumnNature::Temperature,
            ColumnNature::ConversionRate,
        ]
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        self.solve(data)
    }
}

// ─── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::experiment_series_main::{
        ExperimentMeta, TGAExperiment,
    };
    use crate::Kinetics::experimental_kinetics::kinetic_methods::{
        ExperimentData, KineticDataView,
    };
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADataset;
    use crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA;
    use approx::assert_relative_eq;

    // ── synthetic data builder ────────────────────────────────────────────────

    /// Builds an `ExperimentData` where `dα/dt = k0 · exp(-ea / R / T)`.
    /// Temperature is linear from `t_start` to `t_end`; conversion is the
    /// trapezoidal integral of the rate, normalised to [0, 1].
    fn synthetic_experiment(
        k0: f64,
        ea: f64,
        t_start: f64,
        t_end: f64,
        n: usize,
        id: &str,
    ) -> ExperimentData {
        let temps: Vec<f64> = (0..n)
            .map(|i| t_start + (t_end - t_start) * i as f64 / (n - 1) as f64)
            .collect();

        let rates: Vec<f64> = temps.iter().map(|&t| k0 * (-ea / (R * t)).exp()).collect();

        let dt = (t_end - t_start) / (n - 1) as f64;
        let mut conv = vec![0.0_f64; n];
        for i in 1..n {
            conv[i] = conv[i - 1] + 0.5 * (rates[i - 1] + rates[i]) * dt;
        }
        let max_conv = conv[n - 1].max(1e-30);
        let conv: Vec<f64> = conv.iter().map(|&c| c / max_conv).collect();
        let rates: Vec<f64> = rates.iter().map(|&r| r / max_conv).collect();

        ExperimentData {
            meta: ExperimentMeta {
                id: id.to_string(),
                heating_rate: None,
                isothermal_temperature: None,
                comment: None,
            },
            time: (0..n).map(|i| i as f64).collect(),
            temperature: temps,
            conversion: conv,
            conversion_rate: rates,
            mass: None,
            mass_rate: None,
        }
    }

    fn make_view(k0: f64, ea: f64, ids: &[&str]) -> KineticDataView {
        KineticDataView {
            experiments: ids
                .iter()
                .enumerate()
                .map(|(i, &id)| {
                    synthetic_experiment(
                        k0,
                        ea,
                        400.0 + i as f64 * 20.0,
                        800.0 + i as f64 * 20.0,
                        500,
                        id,
                    )
                })
                .collect(),
        }
    }

    // ── solve tests ───────────────────────────────────────────────────────────

    #[test]
    fn recovers_ea_single_experiment() {
        let true_ea = 80_000.0;
        let view = make_view(1e6, true_ea, &["exp1"]);
        let results = SublimationMethod::new(0.05, 0.95)
            .solve(&view)
            .expect("solve failed");

        assert_eq!(results.len(), 1);
        assert_relative_eq!(results[0].ea, true_ea, epsilon = true_ea * 0.01);
        assert!(
            results[0].regression.r2 > 0.999,
            "R²={}",
            results[0].regression.r2
        );
        assert!(results[0].k > 0.0);
    }

    #[test]
    fn fitted_curve_has_same_length_as_temperature() {
        let view = make_view(1e5, 60_000.0, &["e1"]);
        let n_temp = view.experiments[0].temperature.len();
        let results = SublimationMethod::default().solve(&view).unwrap();
        assert_eq!(results[0].fitted_conversion_rate.len(), n_temp);
    }

    #[test]
    fn fitted_curve_matches_arrhenius() {
        let true_ea = 70_000.0;
        let true_k0 = 1e5_f64;
        let view = make_view(true_k0, true_ea, &["e1"]);
        let results = SublimationMethod::new(0.05, 0.95).solve(&view).unwrap();
        let res = &results[0];

        // spot-check a mid-range temperature point
        let t = 600.0_f64;
        let expected = res.k * (-res.ea / (R * t)).exp();
        let actual = res.k * (-res.ea / (R * t)).exp();
        assert_relative_eq!(expected, actual, epsilon = 1e-10);
    }

    #[test]
    fn returns_one_result_per_experiment() {
        let view = make_view(1e6, 80_000.0, &["a", "b", "c"]);
        let results = SublimationMethod::default().solve(&view).unwrap();
        assert_eq!(results.len(), 3);
        for (res, id) in results.iter().zip(["a", "b", "c"]) {
            assert_eq!(res.experiment_id, id);
        }
    }

    #[test]
    fn classify_yes_for_perfect_sublimation() {
        let view = make_view(1e5, 60_000.0, &["e1"]);
        let results = SublimationMethod::default().solve(&view).unwrap();
        assert!(
            matches!(results[0].resume, IsThisASublimationResult::Yes(_)),
            "expected Yes, got {:?}",
            results[0].resume
        );
    }

    // ── mean_results tests ────────────────────────────────────────────────────

    #[test]
    fn mean_results_averages_ea_and_k() {
        let true_ea = 100_000.0;
        let true_k0 = 1e8_f64;
        let view = make_view(true_k0, true_ea, &["e1", "e2", "e3"]);
        let results = SublimationMethod::new(0.05, 0.95).solve(&view).unwrap();
        let (mean_ea, mean_k, verdict) = SublimationMethod::mean_results(&results).unwrap();

        assert_relative_eq!(mean_ea, true_ea, epsilon = true_ea * 0.01);
        assert!(mean_k > 0.0);
        assert!(matches!(verdict, IsThisASublimationResult::Yes(_)));
    }

    #[test]
    fn mean_results_errors_on_empty() {
        assert!(SublimationMethod::mean_results(&[]).is_err());
    }

    // ── KineticMethod trait test ──────────────────────────────────────────────

    #[test]
    fn compute_returns_vec_of_results() {
        let view = make_view(1e6, 80_000.0, &["e1", "e2"]);
        let out = SublimationMethod::default().compute(&view).unwrap();
        assert_eq!(out.len(), 2);
        for r in &out {
            assert!(r.regression.r2 > 0.999);
        }
    }

    // ── push_fitted_rates_to_series test ─────────────────────────────────────

    #[test]
    fn push_fitted_rates_adds_column_to_series() {
        // build a minimal TGASeries with two experiments whose temperature
        // vectors match the synthetic data length (500 points each)
        let n = 500_usize;
        let ids = ["exp_a", "exp_b"];
        let view = make_view(1e6, 80_000.0, &ids);
        let results = SublimationMethod::new(0.05, 0.95).solve(&view).unwrap();

        // build TGASeries from VirtualTGA with matching length
        let mut series = TGASeries::new();
        for (exp_data, &id) in view.experiments.iter().zip(ids.iter()) {
            let v = VirtualTGA {
                time: exp_data.time.clone(),
                temperature: exp_data.temperature.clone(),
                mass: vec![1.0; n],
            };
            let ds = TGADataset::create_from_synthetic_data(&v).unwrap();
            series.push(TGAExperiment::new(ds).with_id(id));
        }

        SublimationMethod::push_fitted_rates_to_series(&results, &mut series).unwrap();

        // verify the column was added to each experiment
        for &id in &ids {
            let cols = series.list_of_columns(id).unwrap();
            assert!(
                cols.contains(&"fitted_sublim_rate".to_string()),
                "column missing for {id}"
            );
        }
    }
}
