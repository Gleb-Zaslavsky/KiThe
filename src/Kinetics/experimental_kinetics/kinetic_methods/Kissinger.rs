//! # Kissinger method
//!
//! ## What this module does
//! Implements the **Kissinger method** for estimating the global activation
//! energy Eₐ from a series of non-isothermal TGA experiments run at different
//! heating rates β.
//!
//! Unlike isoconversional methods (OFW, KAS, Starink, Friedman) which produce
//! Eₐ(α) curves, the Kissinger method yields a **single** Eₐ value by
//! exploiting the shift of the reaction-rate peak temperature Tₘ with β.
//! It is therefore best suited for single-step reactions where Eₐ is
//! approximately constant over the whole conversion range.
//!
//! ## Pipeline
//! ```text
//! KineticDataView  (one ExperimentData per heating rate)
//!      │
//!      ▼
//! ConversionGridBuilder::build_nonisothermal
//!      │
//!      ▼
//! extract_reaction_peaks  →  Vec<ReactionPeak>  (one Tₘ per experiment)
//!      │
//!      ▼
//! KissingerSolver::solve  →  linear regression of ln(β/Tₘ²) vs 1/Tₘ
//!      │
//!      ▼
//! IsoconversionalResult  (single IsoLayerResult, eta = NaN)
//! ```
//!
//! ## Main data structures
//! - [`KissingerSolver`] — stateless solver; takes a `ConversionGrid`,
//!   extracts peaks, and runs the regression.
//! - [`Kissinger`] — zero-size struct implementing `KineticMethod`; owns the
//!   grid-building step.
//! - [`ReactionPeak`] — temperature, conversion, rate, and time at the
//!   maximum d α/dt for one experiment.
//!
//! ## Math
//! At the peak of the reaction rate d α/dt the derivative d²α/dt² = 0.
//! For a first-order Arrhenius reaction this gives:
//!
//! ```text
//! ln(β / Tₘ²) = ln(A·R / Eₐ) − (Eₐ / R) · (1 / Tₘ)
//! ```
//!
//! A linear regression of `ln(β / Tₘ²)` vs `1/Tₘ` across all heating rates
//! yields slope = `−Eₐ/R`, so `Eₐ = −slope · R`.
//!
//! ## Non-trivial technique
//! Peak detection is a simple argmax over the `conversion_rate` row of the
//! `ConversionGrid`.  Because the grid is already resampled onto a uniform α
//! axis, the argmax corresponds to the α layer with the highest rate, which
//! is a robust proxy for the true kinetic peak even when the raw signal is
//! slightly noisy.  No smoothing or derivative computation is required.

use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::IsoLayerResult;
use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::IsoconversionalResult;
use crate::Kinetics::experimental_kinetics::kinetic_methods::kinetic_regression::linear_regression;
use crate::Kinetics::experimental_kinetics::kinetic_methods::*;

/// Stateless solver for the Kissinger method.
///
/// Call [`KissingerSolver::solve`] with a pre-built [`ConversionGrid`] to
/// obtain a single [`IsoLayerResult`] carrying the global Eₐ.
pub struct KissingerSolver;

impl KissingerSolver {
    /// Extracts the reaction-rate peak temperature Tₘ for each experiment,
    /// then regresses `ln(β/Tₘ²)` vs `1/Tₘ` to recover Eₐ.
    ///
    /// Returns an [`IsoLayerResult`] with `eta = NaN` (the Kissinger method
    /// does not resolve Eₐ as a function of conversion).
    pub fn solve(grid: &ConversionGrid) -> Result<IsoLayerResult, TGADomainError> {
        let r = 8.314462618;

        let peaks = extract_reaction_peaks(grid);

        let mut x = Vec::with_capacity(peaks.len());
        let mut y = Vec::with_capacity(peaks.len());

        for (i, peak) in peaks.iter().enumerate() {
            let beta = grid.meta[i]
                .heating_rate
                .ok_or_else(|| TGADomainError::InvalidOperation("Heating rate missing".into()))?;

            let Tm = peak.temperature;

            x.push(1.0 / Tm);
            y.push((beta / (Tm * Tm)).ln());
        }

        let reg = linear_regression(&Array1::from_vec(x), &Array1::from_vec(y));

        let ea = -reg.slope * r;

        Ok(IsoLayerResult {
            eta: f64::NAN,
            ea,
            k: None,
            regression: reg,
        })
    }
}

/// `KineticMethod` entry point for the Kissinger method.
///
/// Builds a non-isothermal `ConversionGrid` and delegates to
/// [`KissingerSolver::solve`].  Requires at least 2 experiments with
/// distinct heating rates and non-empty conversion-rate columns.
pub struct Kissinger;

impl KineticMethod for Kissinger {
    type Output = IsoLayerResult;

    fn name(&self) -> &'static str {
        "Kissinger"
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 2,
            needs_conversion: false,
            needs_conversion_rate: true,
            needs_temperature: true,
            needs_heating_rate: true,
        }
    }

    /// Validates the view against `requirements` before computing.
    fn check_input(&self, data: &KineticDataView) -> Result<(), TGADomainError> {
        check_requirements(data, &self.requirements())
    }

    /// Builds the grid, extracts peaks, runs the Kissinger regression, and
    /// wraps the single result in an [`IsoconversionalResult`].
    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        let grid = ConversionGridBuilder::new().build_nonisothermal(data)?;

        let layer = KissingerSolver::solve(&grid)?;

        Ok(layer)
    }
}

/// Temperature, conversion, rate, and time at the d α/dt maximum for one
/// non-isothermal experiment.
pub struct ReactionPeak {
    /// Temperature at the rate peak (K).
    pub temperature: f64,
    /// Conversion α at the rate peak.
    pub conversion: f64,
    /// Maximum d α/dt value (s⁻¹).
    pub rate: f64,
    /// Time at the rate peak (s).
    pub time: f64,
}

/// Finds the d α/dt peak for each experiment row in `grid`.
///
/// Uses a simple argmax over the `conversion_rate` matrix row.  Because the
/// grid is already on a uniform α axis, no smoothing is needed — the argmax
/// is a robust proxy for the true kinetic peak temperature Tₘ.
pub fn extract_reaction_peaks(grid: &ConversionGrid) -> Vec<ReactionPeak> {
    let n_exp = grid.temperature.nrows();

    let mut peaks = Vec::with_capacity(n_exp);

    for i in 0..n_exp {
        let rate_row = grid.conversion_rate.row(i);
        let temp_row = grid.temperature.row(i);
        let time_row = grid.time.row(i);

        let mut max_idx = 0;
        let mut max_rate = rate_row[0];

        for j in 1..rate_row.len() {
            if rate_row[j] > max_rate {
                max_rate = rate_row[j];
                max_idx = j;
            }
        }

        peaks.push(ReactionPeak {
            temperature: temp_row[max_idx],
            conversion: grid.eta[max_idx],
            rate: rate_row[max_idx],
            time: time_row[max_idx],
        });
    }

    peaks
}

/// Estimates the pre-exponential factor A for each experiment from the
/// Kissinger peak data.
///
/// Uses the relation derived from the peak condition:
/// `A = −β · Eₐ / (R · Tₘ² · f'(αₘ)) · exp(Eₐ / (R · Tₘ))`
///
/// where `fprime` is the derivative of the reaction model function f(α)
/// evaluated at the peak conversion αₘ.
pub fn estimate_preexponential(
    peaks: &[ReactionPeak],
    ea: f64,
    beta: &[f64],
    fprime: fn(f64) -> f64,
) -> Vec<f64> {
    let r = 8.314462618;

    peaks
        .iter()
        .zip(beta)
        .map(|(p, &b)| {
            let Tm = p.temperature;
            let fp = fprime(p.conversion);
            (-b * ea / (r * Tm * Tm * fp)) * (ea / (r * Tm)).exp()
        })
        .collect()
}

//================================================================================================================================
/*
найти максимум функции y(x)
вернуть:

x_peak
y_peak
index
ыполнить параболическую интерполяцию по трём точкам
(i-1, i, i+1)
Это стандартный численный приём.
 Параболическая интерполяция
Если точки
(x1,y1)
(x2,y2)
(x3,y3)
то вершина параболы:
x_peak = x2 + 0.5*(y1 - y3)/(y1 - 2y2 + y3)*(x3-x2)
Использование:
let peak =
    PeakDetector::detect_max(
        &exp.temperature,
        &exp.conversion_rate,
    )
    .ok_or_else(|| TGADomainError::InvalidOperation(
        "Peak not found".into(),
    ))?;

let tm = peak.x;
*/

pub struct Peak {
    pub index: usize,

    pub x: f64,
    pub y: f64,
}

pub struct PeakDetector;

impl PeakDetector {
    pub fn detect_max(x: &[f64], y: &[f64]) -> Option<Peak> {
        if x.len() < 3 {
            return None;
        }

        let mut idx = 0;
        let mut ymax = y[0];

        for i in 1..y.len() {
            if y[i] > ymax {
                ymax = y[i];
                idx = i;
            }
        }

        // если максимум на границе
        if idx == 0 || idx == y.len() - 1 {
            return Some(Peak {
                index: idx,
                x: x[idx],
                y: y[idx],
            });
        }

        // соседние точки
        let x1 = x[idx - 1];
        let x2 = x[idx];
        let x3 = x[idx + 1];

        let y1 = y[idx - 1];
        let y2 = y[idx];
        let y3 = y[idx + 1];

        let denom = y1 - 2.0 * y2 + y3;

        if denom.abs() < 1e-12 {
            return Some(Peak {
                index: idx,
                x: x2,
                y: y2,
            });
        }

        let dx = x3 - x2;

        let delta = 0.5 * (y1 - y3) / denom * dx;

        let x_peak = x2 + delta;

        Some(Peak {
            index: idx,
            x: x_peak,
            y: y2,
        })
    }
}

pub struct PeakWindow {
    pub start_idx: usize,

    pub peak_idx: usize,

    pub end_idx: usize,
}
/*
Детектор пиков

Алгоритм:

1️⃣ найти максимум dα/dt
2️⃣ двигаться влево пока скорость > threshold
3️⃣ двигаться вправо пока скорость > threshold

let peaks = extract_reaction_peaks(grid);

for (i, peak) in peaks.iter().enumerate() {

    let beta = grid.meta[i].heating_rate.unwrap();

    let Tm = peak.temperature;

    x.push(1.0 / Tm);

    y.push((beta / (Tm * Tm)).ln());
}
*/
pub struct ReactionPeak2 {
    pub temperature: f64,

    pub conversion: f64,

    pub rate: f64,

    pub time: f64,

    pub window: PeakWindow,
}

pub fn detect_peak_window(rate: &[f64]) -> PeakWindow {
    let mut peak_idx = 0;
    let mut peak_val = rate[0];

    for i in 1..rate.len() {
        if rate[i] > peak_val {
            peak_val = rate[i];
            peak_idx = i;
        }
    }

    let threshold = 0.01 * peak_val;

    let mut start = peak_idx;

    while start > 0 && rate[start] > threshold {
        start -= 1;
    }

    let mut end = peak_idx;

    while end + 1 < rate.len() && rate[end] > threshold {
        end += 1;
    }

    PeakWindow {
        start_idx: start,

        peak_idx,

        end_idx: end,
    }
}

pub fn extract_reaction_peaks2(grid: &ConversionGrid) -> Vec<ReactionPeak2> {
    let n_exp = grid.temperature.nrows();

    let mut peaks = Vec::with_capacity(n_exp);

    for i in 0..n_exp {
        let rate = grid.conversion_rate.row(i);

        let window = detect_peak_window(rate.as_slice().unwrap());

        let idx = window.peak_idx;

        peaks.push(ReactionPeak2 {
            temperature: grid.temperature[[i, idx]],

            conversion: grid.eta[idx],

            rate: grid.conversion_rate[[i, idx]],

            time: grid.time[[i, idx]],

            window,
        });
    }

    peaks
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion_tests::tests::simulate_tga_first_order2;
    use crate::Kinetics::experimental_kinetics::kinetic_methods_tests::tests::build_view_from_cfg_exact_m0;
    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::base_advanced_config_non_isothermal;
    use std::time::Instant;

    // ── grid-level tests (simulate_tga_first_order2) ─────────────────────────

    #[test]
    fn kissinger_grid_recovers_activation_energy() {
        let e_true = 80_000.0;
        let grid = simulate_tga_first_order2(
            500.0,
            1e5,
            e_true,
            &[2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            0.5,
            50_000,
            Some((5..90).map(|x| x as f64 / 100.0).collect()),
        );
        grid.report();
        let layer = KissingerSolver::solve(&grid).unwrap();
        assert!(layer.ea.is_finite());
        assert!(
            layer.regression.r2 > 0.99,
            "R² = {} is below 0.99",
            layer.regression.r2
        );
        println!(
            "Kissinger Ea = {:.1} J/mol  (true = {:.1})",
            layer.ea, e_true
        );
    }

    #[test]
    fn kissinger_grid_recovers_activation_energy2() {
        let e_true = 140_000.0;
        let grid = simulate_tga_first_order2(
            600.0,
            1e8,
            e_true,
            &[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0],
            1.0,
            175_000,
            Some((5..90).map(|x| x as f64 / 100.0).collect()),
        );
        grid.report();
        let layer = KissingerSolver::solve(&grid).unwrap();
        assert!(layer.ea.is_finite());
        assert!(
            layer.regression.r2 > 0.99,
            "R² = {} is below 0.99",
            layer.regression.r2
        );
        println!(
            "Kissinger Ea = {:.1} J/mol  (true = {:.1})",
            layer.ea, e_true
        );
    }

    #[test]
    fn kissinger_grid_no_alpha_grid() {
        let now = Instant::now();
        let grid = simulate_tga_first_order2(
            500.0,
            1e5,
            85_000.0,
            &[2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            0.5,
            50_000,
            None,
        );
        println!("grid built in {} ms", now.elapsed().as_millis());
        grid.report();
        let now = Instant::now();
        let layer = KissingerSolver::solve(&grid).unwrap();
        println!("solver completed in {} ms", now.elapsed().as_millis());
        assert!(layer.ea.is_finite());
        assert!(
            layer.regression.r2 > 0.99,
            "R² = {} is below 0.99",
            layer.regression.r2
        );
    }

    #[test]
    fn kissinger_grid_wider_heating_rates() {
        let grid = simulate_tga_first_order2(
            550.0,
            5e7,
            150_000.0,
            &[0.5, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            1.0,
            80_000,
            Some((5..90).map(|x| x as f64 / 100.0).collect()),
        );
        grid.report();
        let layer = KissingerSolver::solve(&grid).unwrap();
        assert!(layer.ea.is_finite());
        assert!(
            layer.regression.r2 > 0.99,
            "R² = {} is below 0.99",
            layer.regression.r2
        );
    }

    // ── full-pipeline tests (build_view_from_cfg_exact_m0) ───────────────────

    #[test]
    fn kissinger_compute_with_mock_non_isothermal_data() {
        let cfg = base_advanced_config_non_isothermal(
            400.0,
            1e5,
            90_000.0,
            0.1,
            30_000,
            vec![0.5, 3.0, 5.0],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();

        let result = Kissinger.compute(&view);
        match result {
            Ok(result) => {
                println!(
                    "Kissinger result: Ea = {:.1} J/mol, R² = {:.4}",
                    result.ea, result.regression.r2
                );
                let layer = &result;
                assert!(layer.ea.is_finite());
                assert!((0.0..=1.0).contains(&layer.regression.r2));
            }
            Err(e) => {
                panic!("Kissinger compute failed: {:?}", e);
            }
        }
    }

    #[test]
    fn kissinger_compute_with_mock_non_isothermal_data2() {
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
        println!("view built in {} ms", now.elapsed().as_millis());

        let now = Instant::now();
        let result = Kissinger.compute(&view).unwrap();
        println!("elapsed: {} ms", now.elapsed().as_millis());

        let layer = &result;
        assert!(layer.ea.is_finite());
        assert!(
            layer.regression.r2 > 0.99,
            "R² = {} is below 0.99",
            layer.regression.r2
        );
        println!("Kissinger Ea = {:.1} J/mol", layer.ea);
    }

    #[test]
    fn kissinger_compute_with_mock_non_isothermal_data3() {
        let cfg = base_advanced_config_non_isothermal(
            600.0,
            1e7,
            150_000.0,
            0.5,
            50_000,
            vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let result = Kissinger.compute(&view).unwrap();

        let layer = &result;
        assert!(layer.ea.is_finite());
        assert!(
            layer.regression.r2 > 0.99,
            "R² = {} is below 0.99",
            layer.regression.r2
        );
    }

    // ── structural / edge-case tests ─────────────────────────────────────────

    #[test]
    fn kissinger_result_has_single_layer_with_nan_eta() {
        let grid = simulate_tga_first_order2(
            500.0,
            1e5,
            80_000.0,
            &[2.0, 4.0, 6.0],
            0.5,
            20_000,
            Some((5..90).map(|x| x as f64 / 100.0).collect()),
        );
        let layer = KissingerSolver::solve(&grid).unwrap();
        // Kissinger produces one global result, not per-α layers
        assert!(layer.eta.is_nan());
    }

    #[test]
    fn extract_peaks_returns_one_peak_per_experiment() {
        let grid = simulate_tga_first_order2(
            500.0,
            1e5,
            80_000.0,
            &[1.0, 2.0, 3.0, 4.0],
            0.5,
            20_000,
            Some((5..90).map(|x| x as f64 / 100.0).collect()),
        );
        let peaks = extract_reaction_peaks(&grid);
        assert_eq!(peaks.len(), 4);
        for peak in &peaks {
            assert!(peak.temperature.is_finite() && peak.temperature > 0.0);
            assert!(peak.rate > 0.0);
            assert!(peak.time >= 0.0);
        }
    }
}
