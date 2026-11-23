#[cfg(test)]
mod tests {
    use crate::Kinetics::experimental_kinetics::exp_engine_api::{
        GoldenPipelineConfig, SavGolConfig, SplineConfig,
    };
    use crate::Kinetics::experimental_kinetics::experiment_series_main::*;
    use crate::Kinetics::experimental_kinetics::math::exp_kinetics_smooth_filter::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADataset;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset_test::tests::{
        ds_from_csv, make_csv,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::*;
    use polars::prelude::*;

    #[test]
    fn failed_multi_column_lsq_does_not_commit_partial_changes_to_series() {
        let csv = make_csv(80, 783);
        let mut ds = ds_from_csv(&csv);
        let mut df = ds.frame.clone().collect().unwrap();
        let mut temperature: Vec<f64> = df
            .column("temperature")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        temperature[35] = f64::NAN;
        df.with_column(Series::new("temperature".into(), temperature).into())
            .unwrap();
        ds.frame = df.lazy();

        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(ds).with_id("atomic"))
            .unwrap();
        let before = series.get_experiment_by_id("atomic").unwrap();
        let before_frame = before.dataset.frame.clone().collect().unwrap();
        let before_columns = before.dataset.schema.columns.clone();
        let before_schema_bindings = (
            before.dataset.schema.time.clone(),
            before.dataset.schema.temperature.clone(),
            before.dataset.schema.mass.clone(),
            before.dataset.schema.alpha.clone(),
            before.dataset.schema.eta.clone(),
        );
        let before_provenance = before_columns
            .iter()
            .map(|(name, meta)| (name.clone(), format!("{:?}", meta.provenance)))
            .collect::<std::collections::HashMap<_, _>>();
        let before_history = before.dataset.history_of_operations.feed_lines();
        let before_changed_columns = before.dataset.history_of_operations.columns_has_changed();
        let before_undo_depth = before.dataset.undo_depth();
        let before_series_len = series.len();

        let result = series.lsq_spline_resample_columns_as(
            "atomic",
            "time",
            "time_lsq",
            &["mass", "temperature"],
            &[Some("mass_lsq"), Some("temperature_lsq")],
            40,
            3,
            10,
            RustedSciThe::numerical::data_processing::LSQSplines::SolverKind::Banded,
        );

        assert!(result.is_err());
        let after = series.get_experiment_by_id("atomic").unwrap();
        let after_frame = after.dataset.frame.clone().collect().unwrap();
        assert_eq!(after_frame, before_frame);
        assert_eq!(after.dataset.schema.columns.len(), before_columns.len());
        assert_eq!(after.dataset.schema.time, before_schema_bindings.0);
        assert_eq!(after.dataset.schema.temperature, before_schema_bindings.1);
        assert_eq!(after.dataset.schema.mass, before_schema_bindings.2);
        assert_eq!(after.dataset.schema.alpha, before_schema_bindings.3);
        assert_eq!(after.dataset.schema.eta, before_schema_bindings.4);
        let after_provenance = after
            .dataset
            .schema
            .columns
            .iter()
            .map(|(name, meta)| (name.clone(), format!("{:?}", meta.provenance)))
            .collect::<std::collections::HashMap<_, _>>();
        assert_eq!(after_provenance, before_provenance);
        assert_eq!(
            after.dataset.history_of_operations.feed_lines(),
            before_history
        );
        assert_eq!(
            after.dataset.history_of_operations.columns_has_changed(),
            before_changed_columns
        );
        assert_eq!(after.dataset.undo_depth(), before_undo_depth);
        assert_eq!(series.len(), before_series_len);
        assert!(!after.dataset.schema.columns.contains_key("mass_lsq"));
        assert!(!after.dataset.schema.columns.contains_key("temperature_lsq"));
        assert!(!after.dataset.schema.columns.contains_key("time_lsq"));
    }
    #[test]
    fn golden_pipeline_experiment() {
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

        let virtual_tga = VirtualTGA::generate(&cfg).unwrap();
        let txt = virtual_tga.write_txt().unwrap();

        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();

        let ds = ds_from_csv(&csv);

        let exp = TGAExperiment::new(ds)
            .with_heating_rate(10.0)
            .with_comment("synthetic noisy test");

        let exp = exp // <TGADataset as Clone>::clone(&exp)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap()
            .calibrate_mass(1.0, 0.1, "mass")
            .unwrap()
            .trim_edges(50, 0)
            .unwrap()
            .celsius_to_kelvin()
            .unwrap()
            .seconds_to_hours()
            .unwrap()
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            .sg_filter_column("mass", 11, 3, 0, 1.0)
            .unwrap()
            .rolling_mean("mass", 9)
            .unwrap()
            .trim_null_edges()
            .unwrap()
            .dimensionless_mass(0.0, 10.0, "alpha")
            .unwrap()
            .derive_mass_rate("dalpha_dt")
            .unwrap()
            .trim_null_edges()
            .unwrap();

        let df = exp.dataset.frame.collect().unwrap();

        assert!(df.height() > 1000);
        for col in df.columns() {
            assert_eq!(col.null_count(), 0);
        }
    }

    #[test]
    fn experiment_wrapper_mean_functions() {
        let csv = make_csv(1000, 456);
        let ds = ds_from_csv(&csv);
        let exp = TGAExperiment::new(ds);
        let _ = exp.mean_on_interval("mass", "time", 5.0, 10.0).unwrap();
        let _ = exp
            .mean_on_interval_on_own_range("mass", 5.0, 10.0)
            .unwrap();
        let _ = exp.mean_on_column("mass").unwrap();
    }

    #[test]
    fn series_wrapper_mean_functions() {
        let csv = make_csv(500, 789);
        let ds = ds_from_csv(&csv);
        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(ds.clone()).with_id("a"))
            .unwrap();
        let _ = series
            .mean_on_interval("a", "mass", "time", 0.0, 100.0)
            .unwrap();
        let _ = series
            .mean_on_interval_on_own_range("a", "mass", 0.0, 100.0)
            .unwrap();
        let _ = series.mean_on_column("a", "mass").unwrap();
    }

    #[test]
    fn undo_last_via_series_restores_renamed_column() {
        let csv = make_csv(400, 790);
        let ds = ds_from_csv(&csv);
        let mut series = TGASeries::new();
        series.push(TGAExperiment::new(ds).with_id("a")).unwrap();

        series
            .rename_column("a", "temperature", "temperature_renamed")
            .unwrap();
        assert!(
            series
                .get_experiment_by_id("a")
                .unwrap()
                .dataset
                .schema
                .columns
                .contains_key("temperature_renamed")
        );

        series.undo_last("a").unwrap();
        let exp = series.get_experiment_by_id("a").unwrap();
        assert!(exp.dataset.schema.columns.contains_key("temperature"));
        assert!(
            !exp.dataset
                .schema
                .columns
                .contains_key("temperature_renamed")
        );
        assert!(
            exp.dataset
                .history_of_operations
                .feed_text()
                .contains("undo_last")
        );
    }

    #[test]
    fn series_mutation_rejects_duplicate_id_without_corrupting_index() {
        let csv = make_csv(80, 793);
        let ds = ds_from_csv(&csv);
        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(ds.clone()).with_id("first"))
            .unwrap();
        series
            .push(TGAExperiment::new(ds).with_id("second"))
            .unwrap();
        let original_comment = series
            .get_experiment_by_id("first")
            .unwrap()
            .meta
            .comment
            .clone();

        let err = series.try_mutate_by_id("first", |exp| {
            exp.meta.comment = Some("partial mutation".to_string());
            exp.meta.id = "second".to_string();
            Ok(())
        });
        assert!(matches!(err, Err(TGADomainError::DuplicateExperimentId(id)) if id == "second"));
        assert_eq!(
            series.get_experiment_by_id("first").unwrap().meta.id,
            "first"
        );
        assert_eq!(
            series.get_experiment_by_id("second").unwrap().meta.id,
            "second"
        );
        assert_eq!(
            series.get_experiment_by_id("first").unwrap().meta.comment,
            original_comment
        );
    }

    #[test]
    fn series_index_stays_consistent_after_push_rename_and_removal() {
        let csv = make_csv(20, 795);
        let ds = ds_from_csv(&csv);
        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(ds.clone()).with_id("first"))
            .unwrap();
        series
            .push(TGAExperiment::new(ds).with_id("second"))
            .unwrap();

        assert_eq!(series.index_by_id("first").unwrap(), 0);
        assert_eq!(series.index_by_id("second").unwrap(), 1);

        series.set_experiment_id(1, "renamed").unwrap();
        assert!(matches!(
            series.index_by_id("second"),
            Err(TGADomainError::ExperimentNotFound(_))
        ));
        assert_eq!(series.index_by_id("renamed").unwrap(), 1);

        series.remove_experiment_by_id("first").unwrap();
        assert_eq!(series.index_by_id("renamed").unwrap(), 0);
        assert_eq!(series.get_experiment(0).unwrap().meta.id, "renamed");

        series.drop_experiment(0).unwrap();
        assert!(series.is_empty());
        assert!(matches!(
            series.index_by_id("renamed"),
            Err(TGADomainError::ExperimentNotFound(_))
        ));
    }

    #[test]
    fn failed_try_mutate_by_id_discards_all_candidate_changes() {
        let csv = make_csv(24, 794);
        let ds = ds_from_csv(&csv);
        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(ds).with_id("atomic"))
            .unwrap();
        let original_comment = series
            .get_experiment_by_id("atomic")
            .unwrap()
            .meta
            .comment
            .clone();
        let original_history = series
            .get_experiment_by_id("atomic")
            .unwrap()
            .dataset
            .history_of_operations
            .vector_of_operations
            .len();

        let result = series.try_mutate_by_id("atomic", |exp| {
            exp.meta.comment = Some("must not commit".to_string());
            exp.dataset.frame = exp
                .dataset
                .frame
                .clone()
                .with_columns([lit(1.0).alias("partial_column")]);
            Err::<(), TGADomainError>(TGADomainError::InvalidOperation("injected failure".into()))
        });

        assert!(matches!(result, Err(TGADomainError::InvalidOperation(_))));
        let unchanged = series.get_experiment_by_id("atomic").unwrap();
        assert_eq!(unchanged.meta.comment, original_comment);
        assert_eq!(
            unchanged
                .dataset
                .history_of_operations
                .vector_of_operations
                .len(),
            original_history
        );
        assert!(
            unchanged
                .dataset
                .frame
                .clone()
                .collect()
                .unwrap()
                .column("partial_column")
                .is_err()
        );
    }

    #[test]
    fn series_bindings_update_in_place_and_keep_undo_history() {
        let csv = make_csv(24, 796);
        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(ds_from_csv(&csv)).with_id("binding"))
            .unwrap();
        let initial_history_len = series
            .get_experiment_by_id("binding")
            .unwrap()
            .dataset
            .history_of_operations
            .vector_of_operations
            .len();

        series.bind_time("binding", "time", Unit::Second).unwrap();
        series
            .bind_temperature("binding", "temperature", Unit::Celsius)
            .unwrap();
        series
            .bind_mass("binding", "mass", Unit::Milligram)
            .unwrap();
        series
            .bind_column("binding", ColumnRole::Time, "time", Unit::Second)
            .unwrap();

        let dataset = &series.get_experiment_by_id("binding").unwrap().dataset;
        assert_eq!(dataset.schema.time.as_deref(), Some("time"));
        assert_eq!(dataset.schema.temperature.as_deref(), Some("temperature"));
        assert_eq!(dataset.schema.mass.as_deref(), Some("mass"));
        assert_eq!(dataset.undo_depth(), 4);
        assert_eq!(
            dataset.history_of_operations.vector_of_operations.len(),
            initial_history_len + 4
        );
        assert_eq!(
            dataset
                .history_of_operations
                .vector_of_operations
                .last()
                .unwrap()
                .operation_name,
            "bind_column"
        );
    }

    #[test]
    fn failed_in_place_binding_leaves_dataset_state_unchanged() {
        let csv = make_csv(24, 797);
        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(ds_from_csv(&csv)).with_id("binding"))
            .unwrap();

        let before = series.get_experiment_by_id("binding").unwrap();
        let old_time = before.dataset.schema.time.clone();
        let old_history_len = before
            .dataset
            .history_of_operations
            .vector_of_operations
            .len();
        let old_undo_depth = before.dataset.undo_depth();

        let result = series.bind_time("binding", "time", Unit::Milligram);

        assert!(matches!(
            result,
            Err(TGADomainError::InvalidUnitForBinding { .. })
        ));
        let after = series.get_experiment_by_id("binding").unwrap();
        assert_eq!(after.dataset.schema.time, old_time);
        assert_eq!(
            after
                .dataset
                .history_of_operations
                .vector_of_operations
                .len(),
            old_history_len
        );
        assert_eq!(after.dataset.undo_depth(), old_undo_depth);
        assert_eq!(after.dataset.schema.columns["time"].unit, Unit::Second);
    }

    #[test]
    fn in_place_series_transforms_preserve_results_history_and_undo() {
        let csv = make_csv(24, 798);
        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(ds_from_csv(&csv)).with_id("ops"))
            .unwrap();

        let initial = series
            .get_experiment_by_id("ops")
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let initial_height = initial.height();
        let initial_mass = initial
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .get(0)
            .unwrap();

        series.scale_column("ops", "mass", 2.0).unwrap();
        let scaled = series
            .get_experiment_by_id("ops")
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        assert!(
            (scaled
                .column("mass")
                .unwrap()
                .f64()
                .unwrap()
                .get(0)
                .unwrap()
                - initial_mass * 2.0)
                .abs()
                < 1e-12
        );
        assert!(
            series
                .get_experiment_by_id("ops")
                .unwrap()
                .dataset
                .can_undo()
        );
        series.undo_last("ops").unwrap();

        series.filter_rows("ops", col("time").lt(lit(0.5))).unwrap();
        let filtered_height = series
            .get_experiment_by_id("ops")
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap()
            .height();
        assert!(filtered_height > 0 && filtered_height < initial_height);
        series.undo_last("ops").unwrap();
        assert_eq!(
            series
                .get_experiment_by_id("ops")
                .unwrap()
                .dataset
                .frame
                .clone()
                .collect()
                .unwrap()
                .height(),
            initial_height
        );
    }

    fn golden_cfg(save_to_new_experiment: bool, del_old_experiment: bool) -> GoldenPipelineConfig {
        GoldenPipelineConfig {
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
            isothermal_temperature: None,
            spline_config: SplineConfig {
                n_points: 300,
                n_internal_points: 200,
                degree: 3,
            },
            save_to_new_experiment,
            del_old_experiment,
        }
    }

    #[test]
    fn isothermal_golden_pipeline_accepts_dataset_without_temperature_column() {
        let file = tempfile::NamedTempFile::new().unwrap();
        let mut csv = String::from("time,mass\n");
        for index in 0..1200 {
            let time = index as f64;
            let mass = 100.0 * (-0.0002 * time).exp();
            csv.push_str(&format!("{time:.1},{mass:.8}\n"));
        }
        std::fs::write(file.path(), csv).unwrap();

        let dataset = TGADataset::from_csv_raw(file.path())
            .unwrap()
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap();
        let experiment = TGAExperiment::new(dataset.clone())
            .with_id("isothermal_without_sensor")
            .with_isothermal_temperature(523.15);
        let mut config = golden_cfg(false, false);
        config.spline_config.n_points = 250;
        config.spline_config.n_internal_points = 1000;

        let missing_mode = TGAExperiment::new(dataset.clone())
            .with_id("missing_thermal_metadata")
            .apply_golden_pipeline(config.clone());
        assert!(matches!(
            missing_mode,
            Err(TGADomainError::TemperatureNotBound)
        ));

        let (processed, _) = experiment.apply_golden_pipeline(config).unwrap();
        let columns = processed.dataset.list_of_columns();
        assert!(columns.iter().any(|column| column == "dalpha_dt"));
        assert!(columns.iter().any(|column| column == "deta_dt"));
        assert!(!columns.iter().any(|column| column == "dT/dt"));
        assert!(processed.dataset.schema.temperature.is_none());
        assert_eq!(
            processed.meta.resolved_thermal_mode(),
            ExperimentThermalMode::Isothermal
        );
    }

    #[test]
    fn isothermal_golden_pipeline_preserves_measured_temperature_without_dT_dt() {
        let file = tempfile::NamedTempFile::new().unwrap();
        let mut csv = String::from("time,temperature,mass\n");
        for index in 0..1200 {
            let time = index as f64;
            let temperature = 523.15 + 0.02 * (time / 60.0).sin();
            let mass = 100.0 * (-0.0002 * time).exp();
            csv.push_str(&format!("{time:.1},{temperature:.5},{mass:.8}\n"));
        }
        std::fs::write(file.path(), csv).unwrap();

        let dataset = TGADataset::from_csv_raw(file.path())
            .unwrap()
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_temperature("temperature", Unit::Kelvin)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap();
        let experiment = TGAExperiment::new(dataset)
            .with_id("isothermal_with_sensor")
            .with_isothermal_temperature(523.15);
        let mut config = golden_cfg(false, false);
        config.spline_config.n_points = 250;
        config.spline_config.n_internal_points = 1000;

        let (processed, _) = experiment.apply_golden_pipeline(config).unwrap();
        let columns = processed.dataset.list_of_columns();
        assert!(columns.iter().any(|column| column == "temperature"));
        assert!(processed.dataset.schema.temperature.is_some());
        assert!(!columns.iter().any(|column| column == "dT/dt"));
        assert!(columns.iter().any(|column| column == "deta_dt"));
    }

    #[test]
    fn series_apply_golden_pipeline_keeps_source_when_configured() {
        let csv = make_csv(1200, 555001);
        let ds = ds_from_csv(&csv);
        let exp = TGAExperiment::new(ds).with_id("exp1");
        let mut series = TGASeries::new();
        series.push(exp).unwrap();

        series
            .apply_golden_pipeline("exp1", golden_cfg(true, false))
            .unwrap();

        assert_eq!(series.len(), 2);
        assert!(series.index_by_id("exp1").is_ok());
        assert!(series.index_by_id("proceeded_exp1").is_ok());
        let source = series.get_experiment_by_id("exp1").unwrap();
        assert!(source.dataset.schema.eta.is_none());
        assert!(source.dataset.schema.deta_dt.is_none());
    }

    #[test]
    fn series_apply_golden_pipeline_drops_source_when_configured() {
        let csv = make_csv(1200, 555002);
        let ds = ds_from_csv(&csv);
        let exp = TGAExperiment::new(ds).with_id("exp1");
        let mut series = TGASeries::new();
        series.push(exp).unwrap();

        series
            .apply_golden_pipeline("exp1", golden_cfg(true, true))
            .unwrap();

        assert_eq!(series.len(), 1);
        assert!(series.index_by_id("exp1").is_err());
        assert!(series.index_by_id("proceeded_exp1").is_ok());
    }

    #[test]
    fn golden_pipeline_output_id_collision_does_not_mutate_source() {
        let csv = make_csv(1200, 555004);
        let source = TGAExperiment::new(ds_from_csv(&csv)).with_id("exp1");
        let occupied = TGAExperiment::new(ds_from_csv(&csv)).with_id("proceeded_exp1");
        let mut series = TGASeries::new();
        series.push(source).unwrap();
        series.push(occupied).unwrap();

        let source_columns = series.list_of_columns("exp1").unwrap();
        let source_history = series
            .get_experiment_by_id("exp1")
            .unwrap()
            .dataset
            .history_of_operations
            .vector_of_operations
            .len();
        let result = series.apply_golden_pipeline("exp1", golden_cfg(true, false));

        assert!(matches!(
            result,
            Err(TGADomainError::DuplicateExperimentId(_))
        ));
        assert_eq!(series.len(), 2);
        assert_eq!(series.list_of_columns("exp1").unwrap(), source_columns);
        assert_eq!(
            series
                .get_experiment_by_id("exp1")
                .unwrap()
                .dataset
                .history_of_operations
                .vector_of_operations
                .len(),
            source_history
        );
        assert!(
            series
                .get_experiment_by_id("exp1")
                .unwrap()
                .dataset
                .schema
                .eta
                .is_none()
        );
    }

    #[test]
    fn series_csv_roundtrip_preserves_minute_unknown_and_dimensionless_units() {
        let csv = make_csv(200, 555003);
        let mut ds = ds_from_csv(&csv)
            .bind_time("time", Unit::Minute)
            .unwrap()
            .add_numeric_column(
                "aux_unknown",
                Unit::Unknown,
                vec![1.0; 200].into(),
                ColumnNature::Unknown,
            )
            .unwrap()
            .add_numeric_column(
                "aux_dimensionless",
                Unit::Dimensionless,
                vec![2.0; 200].into(),
                ColumnNature::Unknown,
            )
            .unwrap()
            .add_numeric_column(
                "fit_curve",
                Unit::Milligram,
                vec![3.0; 200].into(),
                ColumnNature::FittedCurve,
            )
            .unwrap()
            .add_numeric_column(
                "fit_r2",
                Unit::Dimensionless,
                vec![0.99; 200].into(),
                ColumnNature::R2,
            )
            .unwrap()
            .add_numeric_column(
                "concentration_A",
                Unit::MolePerCubicMeter,
                vec![0.5; 200].into(),
                ColumnNature::Concentration,
            )
            .unwrap();
        ds.schema.R2 = Some("fit_r2".to_string());
        ds.schema.columns.get_mut("concentration_A").unwrap().origin = ColumnOrigin::Simulated;

        let mut provenance =
            crate::Kinetics::experimental_kinetics::column_provenance::ColumnProvenance::raw(
                "mass",
            );
        provenance.push_step(
            vec!["time".to_string(), "mass".to_string()],
            "fit_curve",
            crate::Kinetics::experimental_kinetics::column_provenance::ColumnTransformKind::Fitting {
                model: "Linear".to_string(),
                x_col: "time".to_string(),
                y_col: "mass".to_string(),
                method: "LM".to_string(),
                parameters: vec![("a".to_string(), 1.25), ("b".to_string(), 0.5)],
                r2: 0.98,
                tolerance: 1e-8,
                max_iter: 120,
            },
            Some(
                crate::Kinetics::experimental_kinetics::column_provenance::TransformQuality::computed(
                    crate::Kinetics::experimental_kinetics::math::ndarray_statistics::FilterQualityReport {
                        column_raw: "mass".to_string(),
                        column_filtered: "fit_curve".to_string(),
                        roughness_ratio_1st: 0.7,
                        roughness_ratio_2nd: 0.6,
                        rmse: 0.2,
                        normalized_rmse: 0.03,
                        correlation: 0.99,
                    },
                    crate::Kinetics::experimental_kinetics::column_provenance::QualityReference::NamedColumn(
                        "mass".to_string(),
                    ),
                ),
            ),
            Some(9),
            true,
        );
        ds.schema.columns.get_mut("fit_curve").unwrap().provenance = provenance;

        let mut series = TGASeries::new();
        let mut experiment = TGAExperiment::new(ds).with_id("exp1");
        experiment.meta.heating_rate = None;
        experiment.meta.isothermal_temperature = Some(523.15);
        experiment.meta.thermal_mode = ExperimentThermalMode::Isothermal;
        series.push(experiment).unwrap();

        let out = tempfile::NamedTempFile::new().unwrap();
        series.to_csv_series(out.path()).unwrap();

        let content = std::fs::read_to_string(out.path()).unwrap();
        assert!(content.contains("# format_version=2"));
        assert!(content.contains("Minute"));
        assert!(content.contains("Unknown"));
        assert!(content.contains("Dimensionless"));
        assert!(content.contains("isothermal"));

        let restored = TGASeries::from_csv_series(out.path()).unwrap();
        let exp = &restored.experiments()[0];
        assert_eq!(exp.meta.thermal_mode, ExperimentThermalMode::Isothermal);
        assert_eq!(
            exp.meta.resolved_thermal_mode(),
            ExperimentThermalMode::Isothermal
        );
        assert!(
            exp.dataset
                .history_of_operations
                .feed_text()
                .contains("from_csv_series")
        );
        assert_eq!(
            exp.dataset.schema.columns.get("time").unwrap().unit,
            Unit::Minute
        );
        assert_eq!(
            exp.dataset.schema.columns.get("aux_unknown").unwrap().unit,
            Unit::Unknown
        );
        assert_eq!(
            exp.dataset
                .schema
                .columns
                .get("aux_dimensionless")
                .unwrap()
                .unit,
            Unit::Dimensionless
        );
        assert_eq!(
            exp.dataset.schema.columns["time"].nature,
            ColumnNature::Time
        );
        assert_eq!(
            exp.dataset.schema.columns["temperature"].nature,
            ColumnNature::Temperature
        );
        assert_eq!(
            exp.dataset.schema.columns["mass"].nature,
            ColumnNature::Mass
        );
        assert_eq!(
            exp.dataset.schema.columns["fit_curve"].nature,
            ColumnNature::FittedCurve
        );
        assert_eq!(
            exp.dataset.schema.columns["fit_curve"].origin,
            ColumnOrigin::NumericDerived
        );
        assert_eq!(exp.dataset.schema.R2.as_deref(), Some("fit_r2"));
        assert_eq!(
            exp.dataset.schema.columns["fit_r2"].nature,
            ColumnNature::R2
        );
        assert_eq!(
            exp.dataset.schema.columns["concentration_A"].unit,
            Unit::MolePerCubicMeter
        );
        assert_eq!(
            exp.dataset.schema.columns["concentration_A"].nature,
            ColumnNature::Concentration
        );
        assert_eq!(
            exp.dataset.schema.columns["concentration_A"].origin,
            ColumnOrigin::Simulated
        );
        let fit_provenance = &exp.dataset.schema.columns["fit_curve"].provenance;
        assert_eq!(fit_provenance.steps.len(), 2);
        assert_eq!(fit_provenance.steps[0].output_column, "mass");
        assert_eq!(fit_provenance.steps[1].output_column, "fit_curve");
        assert_eq!(fit_provenance.steps[1].operation_id, Some(9));
        let quality = fit_provenance.steps[1].quality.as_ref().unwrap();
        assert_eq!(quality.report.as_ref().unwrap().normalized_rmse, 0.03);
        assert!(matches!(
            quality.compared_to,
            crate::Kinetics::experimental_kinetics::column_provenance::QualityReference::NamedColumn(_)
        ));

        let legacy = content
            .lines()
            .filter_map(|line| {
                if line.starts_with("#META\tPROV\t") {
                    return None;
                }
                if line.starts_with("#META\tEXP\t") {
                    let parts = line.split('\t').collect::<Vec<_>>();
                    return Some(
                        parts
                            .iter()
                            .take(5)
                            .chain(parts.iter().skip(6))
                            .copied()
                            .collect::<Vec<_>>()
                            .join("\t"),
                    );
                }
                if line == "# format_version=2" {
                    return Some("# format_version=1".to_string());
                }
                if line.starts_with("#META\tCOL\t") {
                    Some(line.split('\t').take(7).collect::<Vec<_>>().join("\t"))
                } else {
                    Some(line.to_string())
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let legacy_file = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(legacy_file.path(), legacy).unwrap();
        let restored_legacy = TGASeries::from_csv_series(legacy_file.path()).unwrap();
        assert_eq!(
            restored_legacy.experiments()[0]
                .meta
                .resolved_thermal_mode(),
            ExperimentThermalMode::Isothermal
        );
        let legacy_dataset = &restored_legacy.experiments()[0].dataset;
        assert_eq!(
            legacy_dataset.schema.columns["time"].nature,
            ColumnNature::Time
        );
        assert_eq!(
            legacy_dataset.schema.columns["fit_curve"].nature,
            ColumnNature::Unknown
        );
        assert_eq!(
            legacy_dataset.schema.columns["fit_curve"].origin,
            ColumnOrigin::NumericDerived
        );
        assert!(
            legacy_dataset.schema.columns["fit_curve"]
                .provenance
                .feed_text()
                .contains("legacy series csv")
        );

        let no_version = content
            .lines()
            .filter(|line| !line.starts_with("# format_version="))
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(legacy_file.path(), no_version).unwrap();
        assert_eq!(
            TGASeries::from_csv_series(legacy_file.path())
                .unwrap()
                .len(),
            1
        );
    }

    #[test]
    fn series_csv_roundtrip_escapes_metadata_and_supports_empty_series() {
        let mut series = TGASeries::new();
        let dataset = TGADataset::default();
        let id = "experiment\\with\ttabs\nand\nlines";
        let mut experiment = TGAExperiment::new(dataset).with_id(id);
        experiment.meta.comment = Some("first\tline\nsecond\\tail\r".to_string());
        series.push(experiment).unwrap();

        let out = tempfile::NamedTempFile::new().unwrap();
        series.to_csv_series(out.path()).unwrap();
        let restored = TGASeries::from_csv_series(out.path()).unwrap();
        let restored_experiment = restored.get_experiment_by_id(id).unwrap();
        assert_eq!(
            restored_experiment.meta.comment.as_deref(),
            Some("first\tline\nsecond\\tail\r")
        );

        let empty_path = tempfile::NamedTempFile::new().unwrap();
        TGASeries::new().to_csv_series(empty_path.path()).unwrap();
        let empty = TGASeries::from_csv_series(empty_path.path()).unwrap();
        assert!(empty.is_empty());

        std::fs::write(empty_path.path(), "# format_version=999\n").unwrap();
        assert!(matches!(
            TGASeries::from_csv_series(empty_path.path()),
            Err(TGADomainError::InvalidOperation(message)) if message.contains("Unsupported series CSV format version")
        ));

        std::fs::write(
            empty_path.path(),
            "# format_version=1\n# format_version=1\n",
        )
        .unwrap();
        assert!(matches!(
            TGASeries::from_csv_series(empty_path.path()),
            Err(TGADomainError::InvalidOperation(message)) if message.contains("more than once")
        ));

        std::fs::write(empty_path.path(), "# format_version=not-a-number\n").unwrap();
        assert!(matches!(
            TGASeries::from_csv_series(empty_path.path()),
            Err(TGADomainError::InvalidOperation(message)) if message.contains("Invalid series CSV format version")
        ));
    }

    #[test]
    fn selected_series_csv_preserves_only_selected_experiment_columns_and_metadata() {
        let mut series = TGASeries::new();
        let selected_source = make_csv(24, 801);
        let omitted_source = make_csv(19, 802);
        let mut selected = TGAExperiment::new(ds_from_csv(&selected_source)).with_id("selected");
        selected.meta.isothermal_temperature = Some(523.15);
        selected.meta.thermal_mode = ExperimentThermalMode::Isothermal;
        selected.meta.comment = Some("kept metadata".to_string());
        series.push(selected).unwrap();
        series
            .push(TGAExperiment::new(ds_from_csv(&omitted_source)).with_id("omitted"))
            .unwrap();

        let path = tempfile::NamedTempFile::new().unwrap();
        let experiment_ids = vec!["selected".to_string()];
        let columns = std::collections::HashMap::from([(
            "selected".to_string(),
            vec!["time".to_string(), "mass".to_string()],
        )]);
        series
            .to_csv_series_selected(path.path(), &experiment_ids, &columns)
            .unwrap();

        let restored = TGASeries::from_csv_series(path.path()).unwrap();
        assert_eq!(
            restored
                .experiments()
                .iter()
                .map(|experiment| experiment.meta.id.as_str())
                .collect::<Vec<_>>(),
            vec!["selected"]
        );
        let experiment = restored.get_experiment_by_id("selected").unwrap();
        assert_eq!(experiment.meta.isothermal_temperature, Some(523.15));
        assert_eq!(
            experiment.meta.thermal_mode,
            ExperimentThermalMode::Isothermal
        );
        assert_eq!(experiment.meta.comment.as_deref(), Some("kept metadata"));
        assert_eq!(experiment.dataset.schema.time.as_deref(), Some("time"));
        assert_eq!(experiment.dataset.schema.mass.as_deref(), Some("mass"));
        assert_eq!(experiment.dataset.schema.temperature, None);
        let frame = experiment.dataset.frame.clone().collect().unwrap();
        assert_eq!(frame.width(), 2);
        assert_eq!(frame.height(), 24);
    }

    #[test]
    fn plain_txt_export_is_tab_delimited_and_rejects_empty_column_selection() {
        let mut series = TGASeries::new();
        let source = make_csv(12, 803);
        series
            .push(TGAExperiment::new(ds_from_csv(&source)).with_id("plain"))
            .unwrap();

        let path = tempfile::NamedTempFile::new().unwrap();
        let selected = vec!["time".to_string(), "mass".to_string()];
        series
            .to_plain_txt_experiment("plain", path.path(), &selected)
            .unwrap();
        let text = std::fs::read_to_string(path.path()).unwrap();
        assert!(text.lines().next().unwrap().contains("time\tmass"));
        assert!(!text.contains("#META"));
        assert!(!text.contains("KiThe TGA"));
        assert_eq!(text.lines().count(), 13);

        let empty_path = tempfile::tempdir().unwrap().path().join("empty.txt");
        assert!(matches!(
            series.to_plain_txt_experiment("plain", &empty_path, &[]),
            Err(TGADomainError::InvalidOperation(message)) if message.contains("at least one column")
        ));
        assert!(!empty_path.exists());
    }

    #[test]
    fn series_csv_roundtrip_preserves_all_column_metadata_variants() {
        let units = [
            Unit::Second,
            Unit::Minute,
            Unit::Hour,
            Unit::Kelvin,
            Unit::Celsius,
            Unit::MilliVolt,
            Unit::Milligram,
            Unit::MilligramPerSecond,
            Unit::MilligramPerMinute,
            Unit::MilligramPerHour,
            Unit::KelvinPerSecond,
            Unit::KelvinPerMinute,
            Unit::KelvinPerHour,
            Unit::CelsiusPerSecond,
            Unit::CelsiusPerMinute,
            Unit::CelsiusPerHour,
            Unit::PerSecond,
            Unit::PerMinute,
            Unit::PerHour,
            Unit::Gram,
            Unit::Dimensionless,
            Unit::Unknown,
        ];
        let natures = [
            ColumnNature::Time,
            ColumnNature::Temperature,
            ColumnNature::Mass,
            ColumnNature::Conversion,
            ColumnNature::DimensionlessMass,
            ColumnNature::MassRate,
            ColumnNature::ConversionRate,
            ColumnNature::TemperatureRate,
            ColumnNature::DimensionlessMassRate,
            ColumnNature::Unknown,
            ColumnNature::ActivationEnergy,
            ColumnNature::PredexFactor,
            ColumnNature::R2,
            ColumnNature::FittedCurve,
        ];
        let origins = [
            ColumnOrigin::Raw,
            ColumnOrigin::PolarsDerived,
            ColumnOrigin::NumericDerived,
            ColumnOrigin::Imported,
        ];

        let synthetic = VirtualTGA {
            time: vec![0.0, 1.0, 2.0],
            temperature: vec![300.0, 301.0, 302.0],
            mass: vec![10.0, 9.9, 9.8],
        };
        let mut dataset = TGADataset::create_from_synthetic_data(&synthetic).unwrap();
        let mut column_names = Vec::with_capacity(units.len());
        for (index, unit) in units.iter().enumerate() {
            let name = format!("metadata_{index}");
            dataset = dataset
                .add_numeric_column(
                    &name,
                    *unit,
                    vec![index as f64; 3].into(),
                    natures[index % natures.len()],
                )
                .unwrap();
            dataset.schema.columns.get_mut(&name).unwrap().origin = origins[index % origins.len()];
            column_names.push(name);
        }

        let mut series = TGASeries::new();
        series
            .push(TGAExperiment::new(dataset).with_id("metadata"))
            .unwrap();
        let path = tempfile::NamedTempFile::new().unwrap();
        series.to_csv_series(path.path()).unwrap();
        let restored = TGASeries::from_csv_series(path.path()).unwrap();
        let restored_dataset = &restored.get_experiment_by_id("metadata").unwrap().dataset;

        for (index, name) in column_names.iter().enumerate() {
            let metadata = &restored_dataset.schema.columns[name];
            assert_eq!(metadata.unit, units[index], "unit for {name}");
            assert_eq!(
                metadata.nature,
                natures[index % natures.len()],
                "nature for {name}"
            );
            assert_eq!(
                metadata.origin,
                origins[index % origins.len()],
                "origin for {name}"
            );
        }
    }
}
