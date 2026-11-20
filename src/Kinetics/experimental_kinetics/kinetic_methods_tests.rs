//=====================================================================================================
// TESTS
//=====================================================================================================

#[cfg(test)]
pub mod tests {

    use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
    use crate::Kinetics::experimental_kinetics::experiment_series2::{
        KineticColumnUnits, UnitedDataset,
    };
    use crate::Kinetics::experimental_kinetics::kinetic_methods::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
        ColumnNature, TGADomainError, Unit,
    };

    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::{
        base_advanced_config_isothermal, base_advanced_config_non_isothermal,
        build_series_from_cfg, build_series_from_cfg_with_m0,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::{AdvancedTGAConfig, KineticModel};
    use approx::assert_relative_eq;
    use polars::df;
    use polars::prelude::IntoLazy;

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

    fn partial_conversion_experiment() -> ExperimentData {
        let mut exp = mock_experiment();
        exp.conversion = vec![0.2, 0.5, 0.8];
        exp
    }

    #[test]
    fn conversion_grid_builders_share_no_extrapolation_contract() {
        let data = KineticDataView {
            experiments: vec![partial_conversion_experiment()],
        };

        // Both the allocating and in-place builder accept a range whose
        // requested endpoints exactly match the measured interval.
        for grid in [
            ConversionGridBuilder::new()
                .eta_range(0.2, 0.8)
                .segments(8)
                .build_withcopy(&data),
            ConversionGridBuilder::new()
                .eta_range(0.2, 0.8)
                .segments(8)
                .build_nonisothermal(&data),
        ] {
            assert!(grid.is_ok(), "exact measured boundaries must be accepted");
        }

        // A requested range extending beyond either measured endpoint must
        // fail identically, even when the generated half-open grid might not
        // happen to sample the out-of-range upper endpoint.
        for (min, max) in [(0.1, 0.5), (0.5, 0.9)] {
            let allocating = ConversionGridBuilder::new()
                .eta_range(min, max)
                .segments(8)
                .build_withcopy(&data);
            let in_place = ConversionGridBuilder::new()
                .eta_range(min, max)
                .segments(8)
                .build_nonisothermal(&data);
            assert!(
                allocating.is_err(),
                "allocating path accepted [{min}, {max}]"
            );
            assert!(in_place.is_err(), "in-place path accepted [{min}, {max}]");
            match in_place {
                Err(TGADomainError::ConversionRangeNotCovered {
                    requested_min,
                    requested_max,
                    available_min,
                    available_max,
                    ..
                }) => {
                    assert_eq!(requested_min, min);
                    assert_eq!(requested_max, max);
                    assert_eq!(available_min, 0.2);
                    assert_eq!(available_max, 0.8);
                }
                _ => panic!("expected typed conversion range error"),
            }
        }
    }

    #[test]
    fn conversion_validation_handles_long_terminal_plateau_in_linear_time() {
        let plateau_len = 100_000;
        let mut exp = partial_conversion_experiment();
        exp.conversion = vec![0.2, 0.5];
        exp.conversion
            .extend(std::iter::repeat(1.0).take(plateau_len));
        exp.time = (0..exp.conversion.len()).map(|idx| idx as f64).collect();
        exp.temperature = vec![350.0; exp.conversion.len()];
        exp.conversion_rate = vec![0.1; exp.conversion.len()];

        let data = KineticDataView {
            experiments: vec![exp],
        };
        assert_eq!(
            ConversionGridBuilder::new()
                .calc_auto_eta_range(&data)
                .unwrap(),
            (0.2, 1.0)
        );
    }

    #[test]
    fn linear_interpolators_reject_extrapolation_but_accept_endpoints() {
        let exp = partial_conversion_experiment();

        for eta_grid in [&[0.1][..], &[0.9][..]] {
            assert!(interpolate_linear(&exp, eta_grid).is_err());

            let mut t = ndarray::Array1::zeros(eta_grid.len());
            let mut temperature = ndarray::Array1::zeros(eta_grid.len());
            let mut rate = ndarray::Array1::zeros(eta_grid.len());
            assert!(
                interpolate_linear_into(
                    &exp,
                    eta_grid,
                    &mut t.view_mut(),
                    &mut temperature.view_mut(),
                    &mut rate.view_mut(),
                )
                .is_err()
            );
        }

        let eta_grid = [0.2, 0.8];
        let (t, temperature, rate) = interpolate_linear(&exp, &eta_grid).unwrap();
        assert_eq!(t, vec![0.0, 10.0]);
        assert_eq!(temperature, vec![300.0, 400.0]);
        assert_eq!(rate, vec![0.1, 0.3]);

        let mut t_into = ndarray::Array1::zeros(eta_grid.len());
        let mut temperature_into = ndarray::Array1::zeros(eta_grid.len());
        let mut rate_into = ndarray::Array1::zeros(eta_grid.len());
        interpolate_linear_into(
            &exp,
            &eta_grid,
            &mut t_into.view_mut(),
            &mut temperature_into.view_mut(),
            &mut rate_into.view_mut(),
        )
        .unwrap();
        assert_eq!(t_into.to_vec(), t);
        assert_eq!(temperature_into.to_vec(), temperature);
        assert_eq!(rate_into.to_vec(), rate);
    }

    #[test]
    fn dt_matrix_keeps_forward_intervals_for_nonuniform_grid() {
        let time = ndarray::array![[0.0, 1.0, 4.0, 10.0]];
        let dt = ConversionGridBuilder::compute_dt_matrix(&time).unwrap();
        assert_eq!(dt.row(0).to_vec(), vec![1.0, 3.0, 6.0, 6.0]);
        // A solver evaluating conversion layer k must consume dt[k - 1].
        assert_eq!(dt[[0, 0]], time[[0, 1]] - time[[0, 0]]);
        assert_eq!(dt[[0, 1]], time[[0, 2]] - time[[0, 1]]);
    }

    #[test]
    fn dt_matrix_rejects_short_or_non_increasing_time_grids() {
        for time in [
            ndarray::array![[0.0]],
            ndarray::array![[0.0, 0.0]],
            ndarray::array![[1.0, 0.0]],
            ndarray::array![[0.0, f64::NAN]],
        ] {
            assert!(ConversionGridBuilder::compute_dt_matrix(&time).is_err());
        }
    }

    #[test]
    fn conversion_grid_rejects_invalid_conversion_data_and_ranges() {
        let mut cases = Vec::new();

        let mut duplicate = mock_experiment();
        duplicate.conversion = vec![0.0, 0.5, 0.5];
        duplicate.time.truncate(3);
        duplicate.temperature.truncate(3);
        duplicate.conversion_rate.truncate(3);
        cases.push(duplicate);

        let mut decreasing = mock_experiment();
        decreasing.conversion = vec![0.0, 0.7, 0.6];
        decreasing.time.truncate(3);
        decreasing.temperature.truncate(3);
        decreasing.conversion_rate.truncate(3);
        cases.push(decreasing);

        let mut not_finite = mock_experiment();
        not_finite.conversion[1] = f64::NAN;
        cases.push(not_finite);

        let mut misaligned = mock_experiment();
        misaligned.time.pop();
        cases.push(misaligned);

        for experiment in cases {
            let data = KineticDataView {
                experiments: vec![experiment],
            };
            assert!(
                ConversionGridBuilder::new()
                    .segments(10)
                    .build_nonisothermal(&data)
                    .is_err()
            );
            assert!(
                ConversionGridBuilder::new()
                    .calc_auto_eta_range(&data)
                    .is_err()
            );
        }

        let data = KineticDataView {
            experiments: vec![mock_experiment()],
        };
        for (min, max) in [(-0.1, 0.8), (0.2, 1.1), (0.7, 0.6), (f64::NAN, 0.8)] {
            assert!(matches!(
                ConversionGridBuilder::new()
                    .eta_range(min, max)
                    .segments(10)
                    .build_nonisothermal(&data),
                Err(TGADomainError::InvalidConversionRange)
            ));
        }
        assert!(
            ConversionGridBuilder::new()
                .safety_margin(-0.01, 0.0)
                .segments(10)
                .build_nonisothermal(&data)
                .is_err()
        );

        let mut noisy_conversion = mock_experiment();
        noisy_conversion.conversion = vec![-0.01, 0.5, 1.01];
        let noisy_data = KineticDataView {
            experiments: vec![noisy_conversion],
        };
        let grid = ConversionGridBuilder::new()
            .segments(10)
            .build_nonisothermal(&noisy_data)
            .expect("small measured overshoots should not make the physical grid invalid");
        assert!(grid.eta.iter().all(|&eta| (0.0..=1.0).contains(&eta)));
    }

    #[test]
    fn kinetic_view_normalizes_minute_celsius_and_rate_units() {
        let frame = df!(
            "exp_id" => ["sample", "sample"],
            "time" => [0.0, 1.0],
            "temperature" => [25.0, 35.0],
            "conversion" => [0.1, 0.2],
            "conversion_rate" => [6.0, 12.0],
            "mass" => [0.01, 0.02],
            "mass_rate" => [0.6, 1.2]
        )
        .unwrap()
        .lazy();
        let meta = ExperimentMeta {
            id: "sample".into(),
            isothermal_temperature: Some(298.15),
            ..ExperimentMeta::default()
        };
        let units = KineticColumnUnits {
            time: Unit::Minute,
            temperature: Unit::Celsius,
            conversion_rate: Unit::PerMinute,
            mass: Unit::Gram,
            mass_rate: Unit::MilligramPerMinute,
        };
        let united = UnitedDataset::new(frame, vec![meta])
            .with_units(std::collections::HashMap::from([("sample".into(), units)]));

        let view = KineticDataView::from_united_dataset_by_nature(
            &united,
            vec![
                ColumnNature::Time,
                ColumnNature::Temperature,
                ColumnNature::Conversion,
                ColumnNature::ConversionRate,
                ColumnNature::Mass,
                ColumnNature::MassRate,
            ],
        )
        .unwrap();
        let exp = &view.experiments[0];
        assert_eq!(exp.time, vec![0.0, 60.0]);
        assert_eq!(exp.temperature, vec![298.15, 308.15]);
        assert_eq!(exp.conversion_rate, vec![0.1, 0.2]);
        assert_eq!(exp.mass.as_deref(), Some(&[10.0, 20.0][..]));
        assert_eq!(exp.mass_rate.as_deref(), Some(&[0.01, 0.02][..]));
        assert!((exp.meta.isothermal_temperature.unwrap() - 298.15).abs() < 1e-12);

        let required_view = KineticDataView::from_united_dataset(&united).unwrap();
        assert_eq!(required_view.experiments[0].time, vec![0.0, 60.0]);
        assert_eq!(
            required_view.experiments[0].temperature,
            vec![298.15, 308.15]
        );
        assert_eq!(required_view.experiments[0].conversion_rate, vec![0.1, 0.2]);
    }

    #[test]
    fn isothermal_kinetic_view_uses_metadata_temperature_without_measured_column() {
        let frame = df!(
            "exp_id" => ["iso", "iso", "iso"],
            "time" => [0.0, 1.0, 2.0],
            "conversion" => [0.0, 0.2, 0.4],
            "conversion_rate" => [0.1, 0.1, 0.1]
        )
        .unwrap()
        .lazy();
        let meta = ExperimentMeta {
            id: "iso".into(),
            isothermal_temperature: Some(523.15),
            thermal_mode: crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentThermalMode::Isothermal,
            ..ExperimentMeta::default()
        };
        let united = UnitedDataset::new(frame, vec![meta]);

        let view = KineticDataView::from_united_dataset_by_nature(
            &united,
            vec![ColumnNature::Time, ColumnNature::Conversion],
        )
        .unwrap();
        let experiment = &view.experiments[0];
        assert_eq!(experiment.temperature, vec![523.15; 3]);
        assert_eq!(experiment.time, vec![0.0, 1.0, 2.0]);
    }

    #[test]
    fn kinetic_view_normalizes_each_experiment_using_its_own_units() {
        let frame = df!(
            "exp_id" => ["minutes", "minutes", "seconds", "seconds"],
            "time" => [0.0, 1.0, 0.0, 60.0],
            "temperature" => [25.0, 35.0, 298.15, 308.15],
            "conversion" => [0.1, 0.2, 0.1, 0.2],
            "conversion_rate" => [6.0, 12.0, 0.1, 0.2],
            "mass" => [0.01, 0.02, 10.0, 20.0],
            "mass_rate" => [0.6, 1.2, 0.01, 0.02]
        )
        .unwrap()
        .lazy();
        let metas = ["minutes", "seconds"]
            .into_iter()
            .map(|id| ExperimentMeta {
                id: id.into(),
                ..ExperimentMeta::default()
            })
            .collect();
        let minute_units = KineticColumnUnits {
            time: Unit::Minute,
            temperature: Unit::Celsius,
            conversion_rate: Unit::PerMinute,
            mass: Unit::Gram,
            mass_rate: Unit::MilligramPerMinute,
        };
        let second_units = KineticColumnUnits {
            time: Unit::Second,
            temperature: Unit::Kelvin,
            conversion_rate: Unit::PerSecond,
            mass: Unit::Milligram,
            mass_rate: Unit::MilligramPerSecond,
        };
        let united =
            UnitedDataset::new(frame, metas).with_units(std::collections::HashMap::from([
                ("minutes".into(), minute_units),
                ("seconds".into(), second_units),
            ]));
        let view = KineticDataView::from_united_dataset_by_nature(
            &united,
            vec![
                ColumnNature::Time,
                ColumnNature::Temperature,
                ColumnNature::Conversion,
                ColumnNature::ConversionRate,
                ColumnNature::Mass,
                ColumnNature::MassRate,
            ],
        )
        .unwrap();

        assert_eq!(view.experiments[0].time, view.experiments[1].time);
        assert_eq!(
            view.experiments[0].temperature,
            view.experiments[1].temperature
        );
        assert_eq!(
            view.experiments[0].conversion_rate,
            view.experiments[1].conversion_rate
        );
        assert_eq!(view.experiments[0].mass, view.experiments[1].mass);
        assert_eq!(view.experiments[0].mass_rate, view.experiments[1].mass_rate);
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

        let mut view = KineticDataView::from_united_dataset(&united)?;
        for exp in &mut view.experiments {
            while exp.conversion.len() > 2 {
                let n = exp.conversion.len();
                let saturated_tail =
                    exp.conversion[n - 1] >= 1.0 - 1e-12 && exp.conversion[n - 2] >= 1.0 - 1e-12;
                if !saturated_tail {
                    break;
                }
                exp.conversion.pop();
                exp.time.pop();
                exp.temperature.pop();
                exp.conversion_rate.pop();
            }
        }
        Ok(view)
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
    fn kinetic_view_rejects_misaligned_optional_series_before_method_checks() {
        let mut experiment = mock_experiment();
        experiment.meta.id = "misaligned-run".to_string();
        experiment.mass_rate = Some(vec![0.1, 0.2]);
        let view = KineticDataView {
            experiments: vec![experiment],
        };

        let err = view.validate_lengths().unwrap_err();
        assert!(matches!(
            err,
            TGADomainError::ExperimentDataLengthMismatch {
                experiment_id,
                column,
                expected: 3,
                actual: 2,
            } if experiment_id == "misaligned-run" && column == "mass_rate"
        ));

        assert!(matches!(
            check_requirements(&view, &KineticRequirements::default()),
            Err(TGADomainError::ExperimentDataLengthMismatch { .. })
        ));
    }

    #[test]
    fn kinetic_view_allows_omitted_optional_columns() {
        let view = KineticDataView {
            experiments: vec![mock_experiment()],
        };
        assert!(view.validate_lengths().is_ok());
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

        let builder = ConversionGridBuilder::new();
        let (eta_min, eta_max) = builder.calc_auto_eta_range(&view).unwrap();

        let grid = builder
            .eta_range(eta_min, eta_max)
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
