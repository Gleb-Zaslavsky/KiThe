#[cfg(test)]
mod tests {

    use crate::Kinetics::experimental_kinetics::exp_engine_api::{
        GoldenPipelineConfig, SavGolConfig, SplineConfig,
    };
    use crate::Kinetics::experimental_kinetics::exp_kinetics_column_manipulation::column_stats;
    use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::*;
    use crate::Kinetics::experimental_kinetics::ndarray_statistics::{
        FilterQualityReport, filter_quality_from_arrays,
    };
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{TGADataset, Unit};
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset_test::tests::{
        ds_from_csv, make_csv,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::{
        AdvancedTGAConfig, ExperimentMode, KineticModel, NoiseConfig, NoiseKind, NoiseModel,
        PipelineInvariantTest, SpikeModel, VirtualTGA, VirtualTGAConfig,
    };
    use RustedSciThe::numerical::data_processing::lowess_wrapper::LowessConfig;
    use RustedSciThe::numerical::data_processing::splines::SplineKind;
    use ndarray::Array1;
    use polars::prelude::*;

    fn build_advanced_synthetic_pair(
        mass_noise: Option<NoiseConfig>,
        temp_noise: Option<NoiseConfig>,
        spikes: Option<SpikeModel>,
        seed: u64,
    ) -> (TGADataset, TGADataset) {
        let base_cfg = AdvancedTGAConfig {
            n_points: 1_800,
            dt: 0.5,
            experiment_mode: ExperimentMode::NonIsothermal {
                t0: 360.0,
                heating_rates: vec![6.0],
            },
            kinetic_model: KineticModel::ArrheniusSingle {
                m0: 100.0,
                k0: 1.0e5,
                e: 82_000.0,
                r: 8.314,
            },
            mass_noise: None,
            temp_noise: None,
            spikes: None,
            seed,
        };

        let clean_simulated = VirtualTGA::generate_series(&base_cfg);
        let mut noisy_cfg = base_cfg.clone();
        noisy_cfg.mass_noise = mass_noise;
        noisy_cfg.temp_noise = temp_noise;
        noisy_cfg.spikes = spikes;

        let noisy_simulated = VirtualTGA::generate_series(&noisy_cfg);

        let clean =
            TGADataset::create_from_synthetic_data(&clean_simulated[0].virtual_tga).unwrap();
        let noisy =
            TGADataset::create_from_synthetic_data(&noisy_simulated[0].virtual_tga).unwrap();
        (clean, noisy)
    }

    fn filter_quality_for_columns(
        reference: &TGADataset,
        candidate: &TGADataset,
        reference_col: &str,
        candidate_col: &str,
    ) -> FilterQualityReport {
        let reference_df = reference.frame.clone().collect().unwrap();
        let candidate_df = candidate.frame.clone().collect().unwrap();

        let reference_values = reference_df.column(reference_col).unwrap().f64().unwrap();
        let candidate_values = candidate_df.column(candidate_col).unwrap().f64().unwrap();

        let mut raw = Vec::new();
        let mut filtered = Vec::new();

        for (reference_value, candidate_value) in reference_values
            .into_iter()
            .zip(candidate_values.into_iter())
        {
            if let (Some(raw_value), Some(filtered_value)) = (reference_value, candidate_value) {
                raw.push(raw_value);
                filtered.push(filtered_value);
            }
        }

        filter_quality_from_arrays(
            &Array1::from_vec(raw),
            &Array1::from_vec(filtered),
            reference_col,
            candidate_col,
        )
        .unwrap()
    }

    fn quality_score(report: &FilterQualityReport) -> f64 {
        report.normalized_rmse
            + report.roughness_ratio_1st
            + report.roughness_ratio_2nd
            + (1.0 - report.correlation)
    }

    fn shuffled_csv_from_make_csv(n_points: usize, seed: u64) -> tempfile::NamedTempFile {
        let source = make_csv(n_points, seed);
        let content = std::fs::read_to_string(source.path()).unwrap();
        let mut lines = content.lines();
        let header = lines.next().unwrap().to_string();
        let mut rows: Vec<String> = lines.map(|line| line.to_string()).collect();
        if rows.len() > 12 {
            rows.swap(10, 11);
            rows.swap(100, 101);
        }

        let shuffled = tempfile::NamedTempFile::new().unwrap();
        let mut text = String::new();
        text.push_str(&header);
        text.push('\n');
        for row in rows {
            text.push_str(&row);
            text.push('\n');
        }
        std::fs::write(shuffled.path(), text).unwrap();
        shuffled
    }

    fn assert_filter_improves_quality(
        label: &str,
        reference: &TGADataset,
        noisy: &TGADataset,
        source_col: &str,
        filtered_col: &str,
        transform: fn(
            TGADataset,
        ) -> Result<
            TGADataset,
            crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
        >,
    ) {
        let raw_report = filter_quality_for_columns(reference, noisy, source_col, source_col);
        let reference_filtered = transform(reference.clone()).unwrap();
        let noisy_filtered = transform(noisy.clone()).unwrap();
        let filtered_report = filter_quality_for_columns(
            &reference_filtered,
            &noisy_filtered,
            filtered_col,
            filtered_col,
        );

        let raw_score = quality_score(&raw_report);
        let filtered_score = quality_score(&filtered_report);

        assert!(
            filtered_score < raw_score,
            "{}: filtered score {:.6} should be smaller than raw score {:.6}",
            label,
            filtered_score,
            raw_score
        );
        assert!(
            filtered_report.roughness_ratio_1st <= raw_report.roughness_ratio_1st * 1.05,
            "{}: first roughness ratio did not improve enough (raw {:.6}, filtered {:.6})",
            label,
            raw_report.roughness_ratio_1st,
            filtered_report.roughness_ratio_1st
        );
    }

    //================================================================================================
    // ROLLING MEAN
    #[test]
    fn rolling_mean_smoothes_signal() {
        let csv = make_csv(100_000, 9876);
        let ds = ds_from_csv(&csv);
        // Use a relatively large window to ensure smoothing effect
        let ds_sm = ds.clone().rolling_mean("mass", 125);
        let df0 = ds.frame.collect().unwrap();
        let df_smooth = ds_sm.frame.collect().unwrap();

        let m0: Vec<f64> = df0
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let m_smooth: Vec<f64> = df_smooth
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        // lengths equal
        assert_eq!(m0.len(), m_smooth.len());
        // variance should not increase (rough heuristic)
        let mean0 = m0.iter().sum::<f64>() / m0.len() as f64;
        let var0 = m0.iter().map(|v| (v - mean0).powi(2)).sum::<f64>() / m0.len() as f64;
        let mean_smooth = m_smooth.iter().sum::<f64>() / m_smooth.len() as f64;
        let var_smooth = m_smooth
            .iter()
            .map(|v| (v - mean_smooth).powi(2))
            .sum::<f64>()
            / m_smooth.len() as f64;
        println!(
            "mean0 = {}, var0 = {}, mean1 = {}, var1 = {}",
            mean0, var0, mean_smooth, var_smooth
        );
        assert!(var_smooth <= var0 * 1.05); // allow small numerical tolerance
    }

    #[test]
    fn smooth_columns_defaults_to_new_column_and_records_provenance() {
        let csv = make_csv(1200, 2024);
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .smooth_columns(&["mass"], SmoothStrategy::RollingMean { window: 21 })
            .unwrap();

        let df = ds2.frame.clone().collect().unwrap();
        assert!(df.column("mass").is_ok());
        assert!(df.column("rolling_mass").is_ok());

        let provenance = ds2.column_provenance_text("rolling_mass").unwrap();
        assert!(provenance.contains("raw("));
        assert!(provenance.contains("rolling_mean(window=21)"));
    }

    #[test]
    fn provenance_chain_keeps_previous_steps_when_transforming_again() {
        let csv = make_csv(1200, 2025);
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .smooth_columns(&["mass"], SmoothStrategy::RollingMean { window: 15 })
            .unwrap()
            .trim_null_edges();
        let ds3 = ds2
            .lowess_filter_column_as("rolling_mass", 0.2, Some("mass_lowess"))
            .unwrap();

        let provenance = ds3.column_provenance_text("mass_lowess").unwrap();
        assert!(provenance.contains("raw("));
        assert!(provenance.contains("rolling_mean(window=15)"));
        assert!(provenance.contains("lowess(frac=0.2)"));
    }
    //================================================================================================
    //              HAMPEL/MAD
    #[test]
    fn hampel_filter_replace_with_median() {
        use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy;
        let csv = make_csv(1000, 777);
        let mut ds = ds_from_csv(&csv);

        // Add artificial outliers
        let df = ds.frame.clone().collect().unwrap();
        let mut mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        mass[100] = mass[100] * 10.0; // outlier
        mass[500] = mass[500] * 0.1; // outlier

        ds.frame = ds.frame.with_column(Series::new("mass".into(), mass).lit());

        let ds_filtered = ds
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap();
        let df_filtered = ds_filtered.frame.collect().unwrap();

        assert_eq!(df_filtered.height(), df.height());
        assert!(df_filtered.column("mass").is_ok());
    }

    #[test]
    fn hampel_filter_drop_outliers() {
        use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy;
        let csv = make_csv(200, 888);
        let mut ds = ds_from_csv(&csv);

        let df = ds.frame.clone().collect().unwrap();
        let mut mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        mass[50] = mass[50] * 20.0; // strong outlier

        ds.frame = ds.frame.with_column(Series::new("mass".into(), mass).lit());

        let ds_filtered = ds
            .hampel_filter("mass", 7, 2.0, HampelStrategy::Drop)
            .unwrap();
        let df_filtered = ds_filtered.frame.collect().unwrap();

        assert!(df_filtered.height() < df.height());
    }

    #[test]
    fn tga_hampel_replace_keeps_length_and_damps_outliers() {
        let csv = make_csv(1000, 777);
        let mut ds = ds_from_csv(&csv);

        // eager modification (SAFE)
        let mut df = ds.frame.clone().collect().unwrap();
        let mut mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        let original = mass[100];
        mass[100] *= 10.0;

        df.replace("mass", Column::new("mass".into(), mass))
            .unwrap();
        ds.frame = df.clone().lazy();

        let ds2 = ds
            .hampel_filter_null_safe("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap();

        let df2 = ds2.frame.collect().unwrap();

        assert_eq!(df2.height(), df.height());

        let filtered = df2.column("mass").unwrap().f64().unwrap().get(100).unwrap();

        assert!(filtered < original * 2.0);
    }

    #[test]
    fn tga_hampel_drop_reduces_rows_consistently() {
        let csv = make_csv(300, 888);
        let mut ds = ds_from_csv(&csv);

        let mut df = ds.frame.clone().collect().unwrap();
        let mut mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        mass[50] *= 20.0;

        df.replace("mass", Column::new("mass".into(), mass))
            .unwrap();
        ds.frame = df.clone().lazy();

        let ds2 = ds
            .hampel_filter_null_safe("mass", 7, 2.0, HampelStrategy::Drop)
            .unwrap();

        let df2 = ds2.frame.collect().unwrap();

        assert!(df2.height() < df.height());

        // column alignment invariant
        assert_eq!(
            df2.column("time").unwrap().len(),
            df2.column("mass").unwrap().len()
        );
    }
    //==========================================================================
    // SAVITSKY-GOLAY
    #[test]
    fn tga_sg_filter_smooths_mass() {
        let csv = make_csv(500, 42);
        let ds = ds_from_csv(&csv);

        let ds2 = ds.sg_filter_column("mass", 11, 3, 0, 1.0).unwrap();

        let df = ds2.frame.collect().unwrap();

        let mass = df.column("mass").unwrap().f64().unwrap();

        let roughness: f64 = mass
            .into_no_null_iter()
            .collect::<Vec<_>>()
            .windows(2)
            .map(|w| (w[1] - w[0]).abs())
            .sum();

        assert!(roughness.is_finite());
    }

    #[test]
    fn lowess_filter_column_overwrites_column() {
        let csv = make_csv(800, 903);
        let ds = ds_from_csv(&csv);

        let df0 = ds.frame.clone().collect().unwrap();
        let m0: Vec<f64> = df0
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let rough0: f64 = m0.windows(2).map(|w| (w[1] - w[0]).abs()).sum();

        let ds2 = ds.lowess_filter_column("mass", 0.2).unwrap();
        let df2 = ds2.frame.collect().unwrap();

        assert_eq!(df2.height(), df0.height());
        assert_eq!(df2.column("mass").unwrap().null_count(), 0);

        let m1: Vec<f64> = df2
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let rough1: f64 = m1.windows(2).map(|w| (w[1] - w[0]).abs()).sum();
        assert!(rough1 <= rough0 * 1.05);
    }

    #[test]
    fn lowess_smooth_columns_as_adds_named_columns() {
        let csv = make_csv(700, 91011);
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .lowess_smooth_columns_as(
                "time",
                &["mass", "temperature"],
                &[Some("mass_lowess"), Some("temperature_lowess")],
                LowessConfig {
                    fraction: 0.25,
                    ..LowessConfig::default()
                },
            )
            .unwrap();

        let df = ds2.frame.collect().unwrap();
        assert!(df.column("mass_lowess").is_ok());
        assert!(df.column("temperature_lowess").is_ok());
        assert_eq!(df.column("mass_lowess").unwrap().len(), df.height());
        assert_eq!(df.column("temperature_lowess").unwrap().len(), df.height());
        assert_eq!(df.column("mass_lowess").unwrap().null_count(), 0);
        assert_eq!(df.column("temperature_lowess").unwrap().null_count(), 0);
    }

    #[test]
    fn smooth_strategy_lowess_path_works() {
        let csv = make_csv(600, 314159);
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .smooth_columns(&["mass"], SmoothStrategy::Lowess { frac: 0.3 })
            .unwrap();

        let df = ds2.frame.collect().unwrap();
        assert!(df.column("mass").is_ok());
        assert_eq!(df.column("mass").unwrap().null_count(), 0);
    }

    fn rolling_mean_mass(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        Ok(ds.rolling_mean_as("mass", 21, Some("rolling_mass")))
    }

    fn hampel_mass(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        ds.hampel_filter_null_safe_as(
            "mass",
            11,
            3.0,
            HampelStrategy::ReplaceWithMedian,
            Some("hampel_mass"),
        )
    }

    fn sg_mass(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        ds.sg_filter_column_as("mass", 15, 3, 0, 1.0, Some("sg_mass"))
    }

    fn lowess_mass(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        ds.lowess_filter_column_as("mass", 0.2, Some("lowess_mass"))
    }

    fn rolling_mean_temperature(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        Ok(ds.rolling_mean_as("temperature", 21, Some("rolling_temperature")))
    }

    fn hampel_temperature(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        ds.hampel_filter_null_safe_as(
            "temperature",
            11,
            3.0,
            HampelStrategy::ReplaceWithMedian,
            Some("hampel_temperature"),
        )
    }

    fn sg_temperature(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        ds.sg_filter_column_as("temperature", 15, 3, 0, 1.0, Some("sg_temperature"))
    }

    fn lowess_temperature(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        ds.lowess_filter_column_as("temperature", 0.2, Some("lowess_temperature"))
    }

    #[test]
    fn smoothing_filters_improve_mass_quality_on_synthetic_noise() {
        let (clean, noisy) = build_advanced_synthetic_pair(
            Some(NoiseConfig {
                kind: NoiseKind::Gaussian { sigma: 0.45 },
            }),
            Some(NoiseConfig {
                kind: NoiseKind::Drift { slope: 0.0009 },
            }),
            Some(SpikeModel {
                probability: 0.02,
                amplitude: 1.8,
            }),
            20240621,
        );

        let scenarios: &[(
            &str,
            &str,
            fn(
                TGADataset,
            ) -> Result<
                TGADataset,
                crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
            >,
        )] = &[
            ("rolling mean", "rolling_mass", rolling_mean_mass),
            ("hampel", "hampel_mass", hampel_mass),
            ("savitzky-golay", "sg_mass", sg_mass),
            ("lowess", "lowess_mass", lowess_mass),
        ];

        for (label, filtered_col, transform) in scenarios {
            assert_filter_improves_quality(label, &clean, &noisy, "mass", filtered_col, *transform);
        }
    }

    #[test]
    fn smoothing_filters_improve_temperature_quality_on_synthetic_noise() {
        let (clean, noisy) = build_advanced_synthetic_pair(
            Some(NoiseConfig {
                kind: NoiseKind::Impulse {
                    probability: 0.01,
                    amplitude: 0.9,
                },
            }),
            Some(NoiseConfig {
                kind: NoiseKind::Flicker { amplitude: 0.08 },
            }),
            None,
            20240622,
        );

        let scenarios: &[(
            &str,
            &str,
            fn(
                TGADataset,
            ) -> Result<
                TGADataset,
                crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
            >,
        )] = &[
            (
                "rolling mean",
                "rolling_temperature",
                rolling_mean_temperature,
            ),
            ("hampel", "hampel_temperature", hampel_temperature),
            ("savitzky-golay", "sg_temperature", sg_temperature),
            ("lowess", "lowess_temperature", lowess_temperature),
        ];

        for (label, filtered_col, transform) in scenarios {
            assert_filter_improves_quality(
                label,
                &clean,
                &noisy,
                "temperature",
                filtered_col,
                *transform,
            );
        }
    }

    fn spline_mass(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        ds.spline_resample_columns_as(
            "time",
            "time_splined",
            &["mass"],
            &[Some("mass_splined")],
            300,
            SplineKind::Cosine,
        )
    }

    fn lsq_spline_mass(
        ds: TGADataset,
    ) -> Result<
        TGADataset,
        crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError,
    > {
        ds.lsq_spline_resample_columns_as(
            "time",
            "time_lsq",
            &["mass"],
            &[Some("mass_lsq")],
            300,
            3,
            20,
            RustedSciThe::numerical::data_processing::LSQSplines::SolverKind::Banded,
        )
    }

    #[test]
    fn resampling_filters_improve_mass_quality_and_keep_monotonic_time() {
        let (clean, noisy) = build_advanced_synthetic_pair(
            Some(NoiseConfig {
                kind: NoiseKind::Gaussian { sigma: 0.45 },
            }),
            Some(NoiseConfig {
                kind: NoiseKind::Gaussian { sigma: 0.12 },
            }),
            Some(SpikeModel {
                probability: 0.01,
                amplitude: 1.2,
            }),
            20240623,
        );

        let spline_clean = spline_mass(clean.clone()).unwrap();
        let spline_noisy = spline_mass(noisy.clone()).unwrap();
        let lsq_clean = lsq_spline_mass(clean.clone()).unwrap();
        let lsq_noisy = lsq_spline_mass(noisy.clone()).unwrap();

        let raw_report = filter_quality_for_columns(&clean, &noisy, "mass", "mass");

        let spline_report = filter_quality_for_columns(
            &spline_clean,
            &spline_noisy,
            "mass_splined",
            "mass_splined",
        );
        let lsq_report = filter_quality_for_columns(&lsq_clean, &lsq_noisy, "mass_lsq", "mass_lsq");

        assert!(quality_score(&spline_report) < quality_score(&raw_report));
        assert!(quality_score(&lsq_report) < quality_score(&raw_report));

        let spline_df = spline_noisy.frame.collect().unwrap();
        let spline_time: Vec<f64> = spline_df
            .column("time_splined")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert!(!spline_time.is_empty());
        assert!(spline_time.iter().all(|v| v.is_finite()));

        let lsq_df = lsq_noisy.frame.collect().unwrap();
        let lsq_time: Vec<f64> = lsq_df
            .column("time_lsq")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert!(!lsq_time.is_empty());
        assert!(lsq_time.iter().all(|v| v.is_finite()));
    }
    //==============================================================================
    // FULL PIPELINE

    #[test]
    fn pipeline_smoke_test() {
        let csv = make_csv(1000, 123);
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .trim_edges(5, 5)
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            .assert_no_nulls("hampel")
            .unwrap()
            .rolling_mean("mass", 9)
            .trim_null_edges()
            .assert_no_nulls("rolling_mean")
            .unwrap()
            .sg_filter_column("mass", 11, 3, 0, 1.0)
            .unwrap()
            .derive_mass_rate("dm_dt")
            .unwrap();

        let df = ds2.frame.collect().unwrap();

        assert!(df.column("dm_dt").is_ok());
        assert!(df.height() > 10);
    }

    #[test]
    fn golden_tga_pipeline_produces_physical_dm_dt() {
        let cfg = VirtualTGAConfig {
            n_points: 10_0,
            dt: 0.1,
            temperature: 600.0,
            temp_noise: NoiseModel { sigma: 0.2 },
            m0: 10.0,
            k: 0.01,
            mass_noise: NoiseModel { sigma: 0.02 },
            spikes: Some(SpikeModel {
                probability: 0.01,
                amplitude: 1.0,
            }),
            seed: 2024,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();
        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .trim_edges(10, 10)
            .hampel_filter("mass", 15, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            .assert_no_nulls("hampel")
            .unwrap()
            .sg_filter_column("mass", 15, 3, 0, 1.0)
            .unwrap()
            .check_nulls_for_operation("SG")
            .derive_mass_rate("dm_dt")
            .unwrap();

        let df = ds2.frame.collect().unwrap();
        let dm_dt = df.column("dm_dt").unwrap().f64().unwrap();
        let dm_dt_vec: Vec<f64> = dm_dt.clone().into_no_null_iter().collect();
        println!("{:?}", dm_dt_vec);
        // физический инвариант: dm/dt < 0
        let neg_fraction =
            dm_dt.into_no_null_iter().filter(|v| *v < 0.0).count() as f64 / df.height() as f64;

        assert!(neg_fraction > 0.9);
    }

    #[test]
    fn golden_pipeline_smoke_test() {
        let cfg = VirtualTGAConfig {
            n_points: 5000,
            dt: 1.0,
            temperature: 400.0,
            temp_noise: NoiseModel { sigma: 2.0 },
            m0: 100.0,
            k: 1e-3,
            mass_noise: NoiseModel { sigma: 0.5 },
            spikes: Some(SpikeModel {
                probability: 0.002,
                amplitude: 20.0,
            }),
            seed: 42,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();
        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();
        let ds = ds_from_csv(&csv);

        let ds = ds
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap()
            // калибровка
            .calibrate_mass(1.0, 0.1, "mass")
            .unwrap()
            // отсечь начальный мусор
            .trim_edges(50, 0)
            // единицы
            .celsius_to_kelvin()
            .seconds_to_hours()
            // выбросы
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            // сглаживание
            .sg_filter_column("mass", 11, 3, 0, 1.0)
            .unwrap()
            .rolling_mean("mass", 9)
            // синхронизация
            .trim_null_edges()
            // степень превращения
            .dimensionless_mass(0.0, 10., "alpha")
            .unwrap()
            .derive_mass_rate("dalpha_dt")
            .unwrap()
            // финальная чистка
            .trim_null_edges();

        let df = ds.frame.collect().unwrap();

        // --- Инварианты ---
        assert!(df.height() > 1000);
        for col in df.columns() {
            assert_eq!(col.null_count(), 0);
        }
    }

    #[test]
    fn golden_pipeline_test2() {
        // --- Virtual experiment ---

        let cfg = VirtualTGAConfig {
            n_points: 5000,
            dt: 1.0,
            temperature: 400.0,
            temp_noise: NoiseModel { sigma: 2.0 },
            m0: 100.0,
            k: 1e-3,
            mass_noise: NoiseModel { sigma: 0.5 },
            spikes: Some(SpikeModel {
                probability: 0.002,
                amplitude: 20.0,
            }),
            seed: 42,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();

        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();

        let ds = ds_from_csv(&csv);

        // --- Full pipeline ---
        let ds = ds
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap()
            // calibration
            .calibrate_mass(1.0, 0.1, "mass")
            .unwrap()
            // trim startup
            .trim_edges(50, 0)
            // units
            .celsius_to_kelvin()
            .seconds_to_hours()
            // outliers
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            // smoothing
            .sg_filter_column("mass", 11, 3, 0, 1.0)
            .unwrap()
            .rolling_mean("mass", 9)
            // synchronization after smoothing
            .trim_null_edges()
            // kinetics
            .dimensionless_mass(0.0, 10.0, "alpha")
            .unwrap()
            .derive_mass_rate("dalpha_dt")
            .unwrap()
            // final cleanup
            .trim_null_edges();
        assert!(ds.get_dalpha_dt().is_ok());
        assert!(ds.get_alpha().is_ok());
        // --- Invariants ---
        let inv = PipelineInvariantTest::from_lazy(&ds.frame).unwrap();

        inv.same_length(&["time", "mass", "temperature", "alpha", "dalpha_dt"]);
        inv.no_nulls(&["time", "mass", "temperature", "alpha", "dalpha_dt"]);
        inv.monotonic_increasing("time");

        assert!(inv.df.height() > 1000);
    }

    #[test]
    fn apply_golden_pipeline_adds_expected_columns_and_reports_them() {
        let csv = make_csv(1200, 20260304);
        let ds = ds_from_csv(&csv)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap();

        let cfg = GoldenPipelineConfig {
            k: 1.0,
            b: 0.0,
            time_cut_before: None,
            time_cut_after: None,
            sav_gol_config: SavGolConfig {
                window_size: 11,
                poly_degree: 3,
                deriv: 0,
                delta: 1.0,
            },
            averaging_time: 1.0,
            spline_config: SplineConfig {
                n_points: 300,
                n_internal_points: 200,
                degree: 3,
            },
            save_to_new_experiment: false,
            del_old_experiment: false,
        };

        let (ds2, new_columns) = ds.apply_golden_pipeline(cfg).unwrap();
        let df = ds2.frame.collect().unwrap();

        for col in [
            "alpha",
            "eta",
            "dT/dt",
            "dalpha_dt",
            "deta_dt",
            "time_splined",
        ] {
            assert!(df.column(col).is_ok(), "missing output column {}", col);
            assert!(
                new_columns.iter().any(|c| c == col),
                "not reported: {}",
                col
            );
        }

        assert_eq!(
            ds2.schema.columns.get("dT/dt").unwrap().unit,
            Unit::KelvinPerHour
        );
        assert_eq!(
            ds2.schema.columns.get("dalpha_dt").unwrap().unit,
            Unit::PerHour
        );
        assert_eq!(
            ds2.schema.columns.get("deta_dt").unwrap().unit,
            Unit::PerHour
        );

        let time_splined = df.column("time_splined").unwrap().f64().unwrap();
        let non_null: Vec<f64> = time_splined.into_iter().flatten().collect();
        assert!(!non_null.is_empty());
        assert!(non_null.iter().all(|v| v.is_finite()));
        assert!(non_null.windows(2).all(|pair| pair[1] > pair[0]));
    }

    #[test]
    fn apply_golden_pipeline_repairs_non_monotonic_time_before_processing() {
        let csv = shuffled_csv_from_make_csv(1200, 20260305);
        let ds = ds_from_csv(&csv)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap();

        let cfg = GoldenPipelineConfig {
            k: 1.0,
            b: 0.0,
            time_cut_before: None,
            time_cut_after: None,
            sav_gol_config: SavGolConfig {
                window_size: 11,
                poly_degree: 3,
                deriv: 0,
                delta: 1.0,
            },
            averaging_time: 1.0,
            spline_config: SplineConfig {
                n_points: 300,
                n_internal_points: 200,
                degree: 3,
            },
            save_to_new_experiment: false,
            del_old_experiment: false,
        };

        let (ds2, _) = ds.apply_golden_pipeline(cfg).unwrap();
        let df = ds2.frame.collect().unwrap();
        let time_splined = df.column("time_splined").unwrap().f64().unwrap();
        let non_null: Vec<f64> = time_splined.into_iter().flatten().collect();
        assert!(non_null.windows(2).all(|pair| pair[1] > pair[0]));
    }

    #[test]
    fn spline_smooth_removes_spikes() {
        let cfg = VirtualTGAConfig {
            n_points: 50_000,
            dt: 1.0,
            temperature: 500.0,
            temp_noise: NoiseModel { sigma: 0.5 },
            m0: 10.0,
            k: 15e-5,
            mass_noise: NoiseModel { sigma: 0.05 },
            spikes: Some(SpikeModel {
                probability: 0.02,
                amplitude: 2.0,
            }),
            seed: 42,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();

        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();

        let ds = ds_from_csv(&csv);

        let ds = ds
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap()
            .trim_edges(10, 10)
            .hampel_filter_null_safe("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            .trim_null_edges()
            .assert_no_nulls("hampel")
            .unwrap()
            .spline_resample_columns("time", "new_time", &["mass"], 300, SplineKind::Cosine)
            .unwrap();
        let _inv = PipelineInvariantTest::from_lazy(&ds.frame).unwrap();
        let mass = ds.get_mass().unwrap();
        println!("mass {:?}", mass);
        let df = ds.frame.collect().unwrap();
        let stats = column_stats(&df, "mass").unwrap();
        println!("columstats {:?}", stats);
        let mass_vals: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_iter()
            .flatten()
            .collect();
        let time_vals: Vec<f64> = df
            .column("new_time")
            .unwrap()
            .f64()
            .unwrap()
            .into_iter()
            .flatten()
            .collect();
        assert!(!mass_vals.is_empty());
        assert!(!time_vals.is_empty());
        assert!(mass_vals.iter().all(|v| v.is_finite()));
        assert!(time_vals.windows(2).all(|pair| pair[1] > pair[0]));
        // --- Invariants ---
    }

    #[test]
    fn lsq_spline_resample_columns_as_adds_padded_columns() {
        let csv = make_csv(1000, 31415);
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .lsq_spline_resample_columns_as(
                "time",
                "time_lsq",
                &["mass"],
                &[Some("mass_lsq")],
                300,
                3,
                20,
                RustedSciThe::numerical::data_processing::LSQSplines::SolverKind::Banded,
            )
            .unwrap();

        let df = ds2.frame.collect().unwrap();
        assert_eq!(df.height(), 1000);

        let mass_lsq = df.column("mass_lsq").unwrap();
        let time_lsq = df.column("time_lsq").unwrap();
        assert_eq!(mass_lsq.null_count(), 700);
        assert_eq!(time_lsq.null_count(), 700);

        let mass_vals: Vec<f64> = mass_lsq.f64().unwrap().into_iter().flatten().collect();
        let time_vals: Vec<f64> = time_lsq.f64().unwrap().into_iter().flatten().collect();
        assert_eq!(mass_vals.len(), 300);
        assert_eq!(time_vals.len(), 300);
        assert!(mass_vals.iter().all(|v| v.is_finite()));
        for i in 1..time_vals.len() {
            assert!(time_vals[i] > time_vals[i - 1]);
        }
    }

    #[test]
    fn lsq_spline_resample_column_default_overwrites_column() {
        let csv = make_csv(600, 2718);
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .lsq_spline_resample_column("time", "time_lsq_default", "mass", 200)
            .unwrap();

        let df = ds2.frame.collect().unwrap();
        assert_eq!(df.height(), 600);
        assert_eq!(df.column("mass").unwrap().null_count(), 400);
        assert_eq!(df.column("time_lsq_default").unwrap().null_count(), 400);
    }
}
