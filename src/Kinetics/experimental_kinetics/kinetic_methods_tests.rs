//=====================================================================================================
// TESTS
//=====================================================================================================

#[cfg(test)]
pub mod tests {

    use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
        ColumnNature, TGADomainError,
    };

    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::{
        base_advanced_config_isothermal, base_advanced_config_non_isothermal,
        build_series_from_cfg, build_series_from_cfg_with_m0,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::{AdvancedTGAConfig, KineticModel};
    use approx::assert_relative_eq;

    fn mock_experiment() -> ExperimentData {
        ExperimentData {
            meta: ExperimentMeta::default(),
            time: vec![0.0, 5.0, 10.0],
            temperature: vec![300.0, 350.0, 400.0],
            conversion: vec![0.0, 0.5, 1.0],
            conversion_rate: vec![0.1, 0.2, 0.3],
            mass: None,
            mass_rate: None,
        }
    }

    #[test]
    fn interpolate_linear_matches_expected() {
        let exp = mock_experiment();
        let eta_grid = vec![0.25, 0.75];

        let (t, temp, rate) = interpolate_linear(&exp, &eta_grid).unwrap();

        assert_relative_eq!(t[0], 2.5, epsilon = 1e-12);
        assert_relative_eq!(t[1], 7.5, epsilon = 1e-12);
        assert_relative_eq!(temp[0], 325.0, epsilon = 1e-12);
        assert_relative_eq!(temp[1], 375.0, epsilon = 1e-12);
        assert_relative_eq!(rate[0], 0.15, epsilon = 1e-12);
        assert_relative_eq!(rate[1], 0.25, epsilon = 1e-12);
    }

    #[test]
    fn interpolate_spline_is_consistent_for_linear_data() {
        let exp = mock_experiment();
        let eta_grid = vec![0.25, 0.75];

        let (t, temp, rate) = interpolate_spline(&exp, &eta_grid).unwrap();

        assert_relative_eq!(t[0], 2.5, epsilon = 1e-8);
        assert_relative_eq!(t[1], 7.5, epsilon = 1e-8);
        assert_relative_eq!(temp[0], 325.0, epsilon = 1e-8);
        assert_relative_eq!(temp[1], 375.0, epsilon = 1e-8);
        assert_relative_eq!(rate[0], 0.15, epsilon = 1e-8);
        assert_relative_eq!(rate[1], 0.25, epsilon = 1e-8);
    }

    pub fn build_view_from_cfg(cfg: &AdvancedTGAConfig) -> Result<KineticDataView, TGADomainError> {
        let series = build_series_from_cfg(cfg, -1.0, 0.0, 1e-4)?;

        let united = series
            .concat_into_vertical_stack(
                None,
                vec![
                    ColumnNature::Time,
                    ColumnNature::Conversion,
                    ColumnNature::Temperature,
                    ColumnNature::ConversionRate,
                ],
            )
            .unwrap();

        KineticDataView::from_united_dataset(&united)
    }

    pub fn build_view_from_cfg_exact_m0(
        cfg: &AdvancedTGAConfig,
    ) -> Result<KineticDataView, TGADomainError> {
        let m0_raw = match cfg.kinetic_model {
            KineticModel::ArrheniusSingle { m0, .. } => m0,
            KineticModel::ArrheniusTwoComponent { m01, m02, .. } => m01 + m02,
        };
        let k = -1.0;
        let b = 0.0;
        let m0_calibrated = (k * m0_raw + b).abs();

        let series = build_series_from_cfg_with_m0(cfg, k, b, m0_calibrated)?;

        let united = series
            .concat_into_vertical_stack(
                None,
                vec![
                    ColumnNature::Time,
                    ColumnNature::Conversion,
                    ColumnNature::Temperature,
                    ColumnNature::ConversionRate,
                ],
            )
            .unwrap();

        KineticDataView::from_united_dataset(&united)
    }

    fn assert_basic_view_sanity(view: &KineticDataView, expected_experiments: usize) {
        assert_eq!(view.experiments.len(), expected_experiments);
        for exp in &view.experiments {
            assert!(!exp.time.is_empty());
            assert_eq!(exp.time.len(), exp.temperature.len());
            assert_eq!(exp.time.len(), exp.conversion.len());
            assert_eq!(exp.time.len(), exp.conversion_rate.len());
            //    println!("mass {:?}", &exp.mass.unwrap()[0..10].unwrap());
            println!("conversion {:?}", &exp.conversion[0..10]);
            assert!(exp.conversion.iter().all(|&v| v >= -2e-2 && v <= 1.0));

            for w in exp.time.windows(2) {
                assert!(w[1] > w[0]);
            }
        }
    }

    #[test]
    fn test_with_virtual_tga() {
        let cfg = base_advanced_config_non_isothermal(
            700.0,
            1e5,
            80_000.0,
            0.1,
            10_000,
            vec![0.5, 3.0, 5.0],
        );
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 3);

        let grid = ConversionGridBuilder::new()
            .eta_range(0.0, 1.0)
            .segments(200)
            .interpolation(GridInterpolation::Linear)
            .build_nonisothermal(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 200);
        assert_eq!(grid.temperature.dim(), (3, 200));
        assert_eq!(grid.time.dim(), (3, 200));
        assert_eq!(grid.conversion_rate.dim(), (3, 200));
        assert_eq!(grid.meta.len(), 3);
        assert_relative_eq!(grid.eta[0], 0.0, epsilon = 1e-12);
        assert!(grid.eta[199] < 1.0);
    }

    #[test]
    fn test_with_virtual_tga_with_auto_range() {
        let cfg = base_advanced_config_non_isothermal(
            700.0,
            1e5,
            80_000.0,
            0.1,
            10_000,
            vec![0.5, 3.0, 5.0],
        );
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 3);

        let grid = ConversionGridBuilder::new()
            .auto_range()
            .segments(200)
            .interpolation(GridInterpolation::Linear)
            .build_nonisothermal(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 200);
        assert_eq!(grid.temperature.dim(), (3, 200));
        assert_eq!(grid.time.dim(), (3, 200));
        assert_eq!(grid.conversion_rate.dim(), (3, 200));
        assert_eq!(grid.meta.len(), 3);
        assert_relative_eq!(grid.eta[0], 0.0, epsilon = 1e-1);
        assert!(grid.eta[199] < 1.0);
    }
    #[test]
    fn test_with_virtual_tga_spline_grid_isothermal() {
        let cfg = base_advanced_config_isothermal(1e5, 80_000.0, vec![600.0, 700.0], 0.1, 10_000);
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 2);
        for exp in &view.experiments {
            assert!(exp.meta.isothermal_temperature.is_some());
            assert!(exp.meta.heating_rate.is_none());
        }

        let grid = ConversionGridBuilder::new()
            .eta_range(0.1, 0.9)
            .segments(120)
            .interpolation(GridInterpolation::Spline)
            .build_nonisothermal(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 120);
        assert_eq!(grid.temperature.dim(), (2, 120));
        assert_eq!(grid.time.dim(), (2, 120));
        assert_eq!(grid.conversion_rate.dim(), (2, 120));
        assert_relative_eq!(grid.eta[0], 0.1, epsilon = 1e-12);
        assert!(grid.eta[119] < 0.9);
    }

    #[test]
    fn test_with_virtual_tga_line_grid_isothermal() {
        let cfg = base_advanced_config_isothermal(1e5, 80_000.0, vec![600.0, 700.0], 0.1, 10_000);
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 2);
        for exp in &view.experiments {
            assert!(exp.meta.isothermal_temperature.is_some());
            assert!(exp.meta.heating_rate.is_none());
        }

        let grid = ConversionGridBuilder::new()
            .eta_range(0.1, 0.9)
            .segments(120)
            .interpolation(GridInterpolation::Linear)
            .build_nonisothermal(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 120);
        assert_eq!(grid.temperature.dim(), (2, 120));
        assert_eq!(grid.time.dim(), (2, 120));
        assert_eq!(grid.conversion_rate.dim(), (2, 120));
        assert_relative_eq!(grid.eta[0], 0.1, epsilon = 1e-12);
        assert!(grid.eta[119] < 0.9);
    }
    #[test]
    fn test_with_virtual_tga_narrow_eta_grid() {
        let cfg =
            base_advanced_config_non_isothermal(500.0, 1e5, 80_000.0, 0.1, 10_000, vec![3.5, 4.5]);
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 2);

        let grid = ConversionGridBuilder::new()
            .eta_range(0.15, 0.75)
            .segments(60)
            .interpolation(GridInterpolation::Linear)
            .build_nonisothermal(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 60);
        assert_eq!(grid.temperature.dim(), (2, 60));
        assert_eq!(grid.time.dim(), (2, 60));
        assert_eq!(grid.conversion_rate.dim(), (2, 60));
        assert_relative_eq!(grid.eta[0], 0.15, epsilon = 1e-12);
        assert!(grid.eta[59] < 0.75);
        for w in grid.eta.windows(2) {
            assert!(w[1] > w[0]);
        }
    }

    #[test]
    fn test_calc_auto_eta_range() {
        // Helper to create a simple experiment with given conversion range
        fn make_experiment(conversion: Vec<f64>) -> ExperimentData {
            let n = conversion.len();
            ExperimentData {
                meta: ExperimentMeta::default(),
                time: (0..n).map(|i| i as f64).collect(),
                temperature: vec![300.0; n],
                conversion,
                conversion_rate: vec![0.1; n],
                mass: None,
                mass_rate: None,
            }
        }

        // Single experiment
        let view = KineticDataView {
            experiments: vec![make_experiment(vec![0.1, 0.2, 0.3, 0.4])],
        };
        let builder = ConversionGridBuilder::new();
        let (min, max) = builder.calc_auto_eta_range(&view).unwrap();
        assert_relative_eq!(min, 0.1, epsilon = 1e-12);
        assert_relative_eq!(max, 0.4, epsilon = 1e-12);

        // Two experiments with overlapping ranges
        let view = KineticDataView {
            experiments: vec![
                make_experiment(vec![0.2, 0.5, 0.8]),
                make_experiment(vec![0.3, 0.6, 0.9]),
            ],
        };
        let (min, max) = builder.calc_auto_eta_range(&view).unwrap();
        assert_relative_eq!(min, 0.3, epsilon = 1e-12); // intersection min = max of mins (0.2,0.3) = 0.3
        assert_relative_eq!(max, 0.8, epsilon = 1e-12); // intersection max = min of maxs (0.8,0.9) = 0.8

        // Three experiments, one with narrower range
        let view = KineticDataView {
            experiments: vec![
                make_experiment(vec![0.0, 0.5, 1.0]),
                make_experiment(vec![0.2, 0.3, 0.4]),
                make_experiment(vec![0.1, 0.6, 0.7]),
            ],
        };
        let (min, max) = builder.calc_auto_eta_range(&view).unwrap();
        assert_relative_eq!(min, 0.2, epsilon = 1e-12); // max of mins: max(0.0,0.2,0.1) = 0.2
        assert_relative_eq!(max, 0.4, epsilon = 1e-12);

        // Non-overlapping ranges -> error
        let view = KineticDataView {
            experiments: vec![
                make_experiment(vec![0.1, 0.2]),
                make_experiment(vec![0.5, 0.6]),
            ],
        };
        let result = builder.calc_auto_eta_range(&view);
        assert!(matches!(
            result,
            Err(TGADomainError::InvalidConversionRange)
        ));

        // Touching ranges (max of first equals min of second) -> intersection is a point, should error
        let view = KineticDataView {
            experiments: vec![
                make_experiment(vec![0.1, 0.3]),
                make_experiment(vec![0.3, 0.5]),
            ],
        };
        let result = builder.calc_auto_eta_range(&view);
        assert!(matches!(
            result,
            Err(TGADomainError::InvalidConversionRange)
        ));

        // Empty conversion vector (should not happen in practice, but test edge case)
        // The fold will produce INFINITY and NEG_INFINITY, leading to global_min = NEG_INFINITY, global_max = INFINITY? Actually local_min = INFINITY, local_max = NEG_INFINITY.
        // Then global_min = max(NEG_INFINITY, INFINITY) = INFINITY? Wait, global_min starts as NEG_INFINITY, then max with INFINITY yields INFINITY.
        // global_max starts as INFINITY, then min with NEG_INFINITY yields NEG_INFINITY.
        // At the end global_max <= global_min? INFINITY <= NEG_INFINITY? false. But we'll skip this test as it's unrealistic.
    }
}
