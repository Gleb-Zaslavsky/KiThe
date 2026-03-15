use super::kinetic_regression::linear_regression;
use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::{
    IsoLayerResult, IsoconversionalResult,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::{
    ConversionGrid, ConversionGridBuilder, GridInterpolation, KineticDataView, KineticMethod,
    KineticRequirements, TGADomainError,
};
use ndarray::Array1;
//============================================================================================================
//          DIFFERENTIAL FRIEDMAN
//===========================================================================================================
pub struct DifferentialFriedmanSolver;

impl Default for DifferentialFriedmanSolver {
    fn default() -> Self {
        Self
    }
}

impl DifferentialFriedmanSolver {
    pub fn solve(&self, grid: &ConversionGrid) -> Result<IsoconversionalResult, TGADomainError> {
        let r = 8.314462618;

        let n_eta = grid.eta.len();
        let n_exp = grid.inv_temperature.nrows();

        let mut layers = Vec::with_capacity(n_eta);

        for k in 0..n_eta {
            let inv_t = grid.inv_temperature.column(k);
            let rate = grid.conversion_rate.column(k);

            let mut x = Vec::with_capacity(n_exp);
            let mut y = Vec::with_capacity(n_exp);

            for (&it, &r) in inv_t.iter().zip(rate.iter()) {
                if r <= 0.0 {
                    continue;
                }

                x.push(it);
                y.push(r.ln());
            }

            if x.len() < 2 {
                continue;
            }

            let reg = linear_regression(&Array1::from_vec(x), &Array1::from_vec(y));

            let ea = -reg.slope * r;

            layers.push(IsoLayerResult {
                eta: grid.eta[k],

                ea,

                k: None,

                regression: reg,
            });
        }

        Ok(IsoconversionalResult {
            method: "Friedman differential",

            layers,
        })
    }
}

pub struct DifferentialFriedman;

impl KineticMethod for DifferentialFriedman {
    type Output = IsoconversionalResult;

    fn name(&self) -> &'static str {
        "Friedman differential"
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        let grid = ConversionGridBuilder::new()
            .interpolation(GridInterpolation::Linear)
            .build(data)?;
        DifferentialFriedmanSolver.solve(&grid)
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 3,
            needs_conversion: true,
            needs_conversion_rate: true,
            needs_temperature: true,
            needs_heating_rate: false,
        }
    }
}
//================================================================================================================================
//              INTEGRAL FRIEDMAN
//=======================================================================================================================
pub struct FriedmanIntegralSolver;

impl FriedmanIntegralSolver {
    pub fn solve(&self, grid: &ConversionGrid) -> Result<IsoconversionalResult, TGADomainError> {
        let r = 8.314462618;

        let n_eta = grid.eta.len();
        let n_exp = grid.inv_temperature.nrows();

        let mut layers = Vec::with_capacity(n_eta);

        for k in 0..n_eta {
            let inv_t = grid.inv_temperature.column(k);
            let t = grid.time.column(k);

            let mut x = Vec::with_capacity(n_exp);
            let mut y = Vec::with_capacity(n_exp);

            for (&it, &time) in inv_t.iter().zip(t.iter()) {
                if time <= 0.0 {
                    continue;
                }

                x.push(it);
                y.push(-time.ln());
            }

            if x.len() < 2 {
                continue;
            }

            let reg = linear_regression(&Array1::from_vec(x), &Array1::from_vec(y));

            let ea = -reg.slope * r;

            layers.push(IsoLayerResult {
                eta: grid.eta[k],

                ea,

                k: None,

                regression: reg,
            });
        }

        Ok(IsoconversionalResult {
            method: "Friedman integral",

            layers,
        })
    }
}

pub struct FriedmanIntegral;

impl KineticMethod for FriedmanIntegral {
    type Output = IsoconversionalResult;

    fn name(&self) -> &'static str {
        "Friedman integral"
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        let grid = ConversionGridBuilder::new()
            .interpolation(GridInterpolation::Linear)
            .build_isothermal(data)?;
        FriedmanIntegralSolver.solve(&grid)
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 3,
            needs_conversion: true,
            needs_conversion_rate: false,
            needs_temperature: true,
            needs_heating_rate: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion_tests::tests::simulate_tga_first_order_isothermal;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::tests::{
         build_view_from_cfg_exact_m0,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::base_advanced_config_isothermal;
    use std::time::Instant;

    // ── helpers ──────────────────────────────────────────────────────────────

    /// Temperatures used across most tests (K).
    fn iso_temps() -> Vec<f64> {
        vec![520.0, 540.0, 560.0, 580.0, 600.0, 620.0]
    }

    // ── grid-level tests (simulate_tga_first_order_isothermal) ───────────────

    #[test]
    fn friedman_integral_grid_recovers_activation_energy() {
        let e_true = 80_000.0;
        let grid = simulate_tga_first_order_isothermal(
            0.0, // t0 unused for isothermal
            1e5,
            e_true,
            &iso_temps(),
            0.5,
            50_000,
            Some((5..90).map(|x| x as f64 / 100.0).collect()),
        );
        grid.report();
        let result = FriedmanIntegralSolver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.90, 50, Some(0.99));
    }

    #[test]
    fn friedman_integral_grid_recovers_activation_energy2() {
        let e_true = 120_000.0;
        let grid = simulate_tga_first_order_isothermal(
            0.0,
            1e8,
            e_true,
            &[550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0],
            0.3,
            150_000,
            Some((5..90).map(|x| x as f64 / 100.0).collect()),
        );
        grid.report();
        let result = FriedmanIntegralSolver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.10, 0.90, 50, Some(0.99));
    }

    #[test]
    fn friedman_integral_grid_no_alpha_grid() {
        let now = Instant::now();
        let e_true = 85_000.0;
        let grid =
            simulate_tga_first_order_isothermal(0.0, 1e5, e_true, &iso_temps(), 0.5, 50_000, None);
        println!("grid built in {} ms", now.elapsed().as_millis());
        grid.report();
        let now = Instant::now();
        let result = FriedmanIntegralSolver.solve(&grid).unwrap();
        println!("solver completed in {} ms", now.elapsed().as_millis());
        result.pretty_print_and_assert(0.05, 0.90, 50, Some(0.99));
    }

    #[test]
    fn friedman_integral_grid_wider_temperature_range() {
        let grid = simulate_tga_first_order_isothermal(
            0.0,
            5e6,
            115_000.0,
            &[500.0, 530.0, 560.0, 590.0, 620.0, 650.0, 680.0],
            1.0,
            80_000,
            Some((5..90).map(|x| x as f64 / 100.0).collect()),
        );
        grid.report();
        let result = FriedmanIntegralSolver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.10, 0.90, 50, Some(0.99));
    }

    // ── full-pipeline tests (build_view_from_cfg_exact_m0) ───────────────────

    #[test]
    fn friedman_integral_compute_with_mock_isothermal_data() {
        let cfg = base_advanced_config_isothermal(1e5, 80_000.0, iso_temps(), 0.1, 10_000);
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let result = FriedmanIntegral.compute(&view).unwrap();
        assert!(!result.layers.is_empty());
        for layer in &result.layers {
            assert!(layer.ea.is_finite());
            assert!((0.0..=1.0).contains(&layer.regression.r2));
        }
    }

    #[test]
    fn friedman_integral_compute_with_mock_isothermal_data2() {
        let now = Instant::now();
        let cfg = base_advanced_config_isothermal(
            1e6,
            100_000.0,
            vec![520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0],
            1.0,
            30_000,
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        println!("view built in {} ms", now.elapsed().as_millis());

        let now = Instant::now();
        let result = FriedmanIntegral.compute(&view).unwrap();
        println!("elapsed: {} ms", now.elapsed().as_millis());
        result.pretty_print_and_assert(0.05, 0.90, 100, Some(0.99));
    }

    #[test]
    fn friedman_integral_compute_with_mock_isothermal_data3() {
        let cfg = base_advanced_config_isothermal(
            1e7,
            150_000.0,
            vec![600.0, 625.0, 650.0, 675.0, 700.0, 725.0],
            0.5,
            50_000,
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let result = FriedmanIntegral.compute(&view).unwrap();
        result.pretty_print_and_assert(0.05, 0.90, 100, Some(0.99));
    }
}

#[cfg(test)]
mod tests_differential {
    use super::*;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion_tests::tests::simulate_tga_first_order2;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::tests::build_view_from_cfg_exact_m0;
    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::base_advanced_config_non_isothermal;
    use std::time::Instant;

    // ── grid-level tests (simulate_tga_first_order2) ─────────────────────────

    #[test]
    fn friedman_differential_grid_recovers_activation_energy() {
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
        let result = DifferentialFriedmanSolver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.90, 50, Some(0.99));
    }

    #[test]
    fn friedman_differential_grid_recovers_activation_energy2() {
        let e_true = 120_000.0;
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
        let result = DifferentialFriedmanSolver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.90, 50, Some(0.99));
    }

    #[test]
    fn friedman_differential_grid_no_alpha_grid() {
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
        let result = DifferentialFriedmanSolver.solve(&grid).unwrap();
        println!("solver completed in {} ms", now.elapsed().as_millis());
        result.pretty_print_and_assert(0.05, 0.90, 50, Some(0.99));
    }

    #[test]
    fn friedman_differential_grid_wider_heating_rates() {
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
        let result = DifferentialFriedmanSolver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.90, 50, Some(0.99));
    }

    // ── full-pipeline tests (build_view_from_cfg_exact_m0) ───────────────────

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
        let result = DifferentialFriedman.compute(&view).unwrap();
        assert!(!result.layers.is_empty());
        for layer in &result.layers {
            let eta = layer.eta;
            if eta > 0.05 && eta < 0.95 {
                assert!(layer.ea.is_finite());
                assert!((0.0..=1.0).contains(&layer.regression.r2));
            }
        }
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
        println!("view built in {} ms", now.elapsed().as_millis());

        let now = Instant::now();
        let result = DifferentialFriedman.compute(&view).unwrap();
        println!("elapsed: {} ms", now.elapsed().as_millis());
        result.pretty_print_and_assert(0.05, 0.90, 100, Some(0.99));
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
        let result = DifferentialFriedman.compute(&view).unwrap();
        result.pretty_print_and_assert(0.05, 0.90, 100, Some(0.99));
    }
}
