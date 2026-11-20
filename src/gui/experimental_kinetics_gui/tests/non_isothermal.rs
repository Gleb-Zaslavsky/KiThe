use super::support::{
    deterministic_config, empty_app, generated, import_csv, read_column, with_app,
    write_simulated_csv,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::isoconversion::IsoconversionalMethod;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
use crate::Kinetics::experimental_kinetics::testing_mod::ExperimentMode;
use egui_kittest::kittest::Queryable;
use tempfile::tempdir;

#[test]
fn story_non_isothermal_three_rate_workflow_runs_friedman_and_roundtrips_export() {
    const E_TRUE: f64 = 85_000.0;
    let mut config = deterministic_config(false, 44_921);
    config.n_points = 360;
    config.experiment_mode = ExperimentMode::NonIsothermal {
        t0: 500.0,
        heating_rates: vec![2.0, 5.0, 10.0],
    };
    config.mass_noise = Some(
        crate::Kinetics::experimental_kinetics::testing_mod::NoiseConfig {
            kind: crate::Kinetics::experimental_kinetics::testing_mod::NoiseKind::Gaussian {
                sigma: 0.03,
            },
        },
    );
    config.temp_noise = None;
    let generated_runs = generated(&config);
    let mut truth_config = config.clone();
    truth_config.mass_noise = None;
    let truth_runs = generated(&truth_config);
    let dir = tempdir().unwrap();
    let mut inputs = Vec::new();
    for run in &generated_runs {
        let beta = run.meta.heating_rate.unwrap();
        let path = dir.path().join(format!("ramp_{beta:.0}.csv"));
        write_simulated_csv(&path, run, true);
        inputs.push((path, beta));
    }

    let mut app = empty_app();
    app.new_experiment_dialog.current_dir = dir.path().to_path_buf();
    let export_dir = tempdir().unwrap();
    let export_path = export_dir.path().join("series_export.csv");
    with_app(app, |harness| {
        for (path, beta) in &inputs {
            import_csv(harness, path, &format!("ramp_{beta:.0}"), None, Some(*beta));
        }
        let ids = harness.state().app.model.list_of_experiments();
        assert_eq!(ids.len(), 3);
        let mut arrhenius_x = Vec::new();
        let mut arrhenius_y = Vec::new();

        for (index, id) in ids.iter().enumerate() {
            let model = &mut harness.state_mut().app.model;
            model.settings.set_calibration_coeffs(-1.0, 40.0);
            model.series.bind_time(id, "t", Unit::Second).unwrap();
            model
                .series
                .bind_temperature(id, "T", Unit::Kelvin)
                .unwrap();
            model.series.bind_mass(id, "m", Unit::MilliVolt).unwrap();
            model
                .calibrate_mass_from_voltage_new_column(id, "mass_mg")
                .expect("calibrate voltage signal");
            model
                .series
                .bind_mass(id, "mass_mg", Unit::Milligram)
                .unwrap();
            model
                .series
                .sg_filter_column_as(id, "mass_mg", 31, 3, 0, 1.0, Some("sg_mass"))
                .expect("smooth mass");
            let quality = model
                .series
                .get_experiment_by_id(id)
                .unwrap()
                .dataset
                .filter_quality_report("mass_mg", "sg_mass")
                .unwrap();
            assert!(
                quality.roughness_ratio_1st < 1.0,
                "mass smoothing ineffective: {quality}"
            );
            let raw_mass = read_column(model, id, "mass_mg");
            let smoothed_mass = read_column(model, id, "sg_mass");
            let clean_mass = &truth_runs[index].virtual_tga.mass;
            let raw_truth_rmse = super::support::rmse(&raw_mass, clean_mass);
            let smooth_truth_rmse = super::support::rmse(&smoothed_mass, clean_mass);
            assert!(
                smooth_truth_rmse < raw_truth_rmse,
                "SG increased mass error against known noise-free signal for {id}: raw={raw_truth_rmse}, filtered={smooth_truth_rmse}"
            );
            super::support::print_filter_report(
                config.seed,
                id,
                "calibrated mass",
                31,
                3,
                &quality,
                Some((raw_truth_rmse, smooth_truth_rmse)),
            );
            model
                .series
                .bind_mass(id, "sg_mass", Unit::Milligram)
                .unwrap();
            model.from_s_to_h_for_experiment(id).unwrap();
            model
                .conversion_for_experiment(id.clone(), 0.005, "eta")
                .unwrap();
            let raw_eta = read_column(model, id, "eta");
            model
                .series
                .sg_filter_column_as(id, "eta", 11, 3, 0, 1.0, Some("eta"))
                .expect("smooth conversion");
            let smoothed_eta = read_column(model, id, "eta");
            let eta_report = super::support::quality_from_slices(
                &raw_eta,
                &smoothed_eta,
                "eta_before_sg",
                "eta_after_sg",
            );
            let truth_eta = super::support::ideal_conversion(
                &truth_runs[index].virtual_tga.time,
                &truth_runs[index].virtual_tga.mass,
                0.005 * 3600.0,
            );
            let raw_eta_truth_rmse = super::support::rmse(&raw_eta, &truth_eta);
            let smooth_eta_truth_rmse = super::support::rmse(&smoothed_eta, &truth_eta);
            super::support::print_filter_report(
                config.seed,
                id,
                "conversion",
                11,
                3,
                &eta_report,
                Some((raw_eta_truth_rmse, smooth_eta_truth_rmse)),
            );
            model.derive_deta_dt(id, Some("deta_dt")).unwrap();
            model
                .series
                .drop_nulls(id)
                .expect("remove derivative edge rows together");
            let raw_rate = read_column(model, id, "deta_dt");
            model
                .series
                .sg_filter_column_as(id, "deta_dt", 11, 3, 0, 1.0, Some("deta_dt"))
                .expect("rate cleanup");
            let smoothed_rate = read_column(model, id, "deta_dt");
            let rate_report = super::support::quality_from_slices(
                &raw_rate,
                &smoothed_rate,
                "deta_dt_before_sg",
                "deta_dt_after_sg",
            );
            assert!(
                rate_report.roughness_ratio_1st < 1.0 && rate_report.roughness_ratio_2nd < 1.0,
                "final SG did not reduce rate roughness for {id}: {rate_report}"
            );
            super::support::print_filter_report(
                config.seed,
                id,
                "conversion rate",
                11,
                3,
                &rate_report,
                None,
            );

            let exp = model.series.get_experiment_by_id(id).unwrap();
            assert!(exp.meta.heating_rate.is_some());
            assert_eq!(exp.dataset.schema.columns["t"].unit, Unit::Hour);
            assert_eq!(exp.dataset.schema.columns["T"].unit, Unit::Kelvin);
            assert_eq!(exp.dataset.schema.columns["mass_mg"].unit, Unit::Milligram);
            assert!(exp.dataset.column_provenance_text("deta_dt").is_some());
            let frame = exp.dataset.frame.clone().collect().unwrap();
            for column in ["t", "T", "mass_mg", "sg_mass", "eta", "deta_dt"] {
                assert_eq!(
                    frame.column(column).unwrap().null_count(),
                    0,
                    "null in {id}.{column}"
                );
            }
            assert!(
                read_column(model, id, "eta")
                    .iter()
                    .all(|value| value.is_finite())
            );
            let eta = read_column(model, id, "eta");
            let temperature = read_column(model, id, "T");
            let rate_per_hour = read_column(model, id, "deta_dt");
            for ((conversion, temp), rate_h) in eta.iter().zip(temperature).zip(rate_per_hour) {
                if (0.1..=0.9).contains(conversion) && rate_h > 0.0 {
                    // This fixture is first order: d(alpha)/dt / (1-alpha) = k(T).
                    let k_per_second = rate_h / 3600.0 / (1.0 - conversion);
                    if k_per_second.is_finite() && k_per_second > 0.0 {
                        arrhenius_x.push(1.0 / temp);
                        arrhenius_y.push(k_per_second.ln());
                    }
                }
            }
        }

        // The kinetic action is driven through the visible method window.
        harness
            .state_mut()
            .app
            .kinetic_methods_state
            .isoconversional_method = IsoconversionalMethod::FriedmanDifferential;
        harness
            .state_mut()
            .app
            .kinetic_methods_state
            .isoconversional_output_id = "friedman_story_result".into();
        harness.run();
        super::support::click_menu_entry(harness, "Kinetic Methods", "Isoconversional");
        harness.get_by_label("Calculate").click();
        harness.run();
        let status = &harness
            .state()
            .app
            .kinetic_methods_state
            .isoconversional_status;
        assert!(
            status.starts_with("Done."),
            "Differential Friedman failed: {status}"
        );
        let result_id = "friedman_story_result";
        assert!(
            harness
                .state()
                .app
                .model
                .list_of_experiments()
                .contains(&result_id.to_string())
        );
        let ea = read_column(&harness.state().app.model, result_id, "Ea");
        let eta = read_column(&harness.state().app.model, result_id, "eta");
        let finite_ea = ea
            .into_iter()
            .zip(eta)
            .filter_map(|(value, conversion)| {
                (value.is_finite() && (0.1..=0.9).contains(&conversion)).then_some(value)
            })
            .collect::<Vec<_>>();
        assert!(!finite_ea.is_empty());
        let mean_ea = finite_ea.iter().sum::<f64>() / finite_ea.len() as f64;
        println!(
            "[STORY][non-isothermal-method] seed={} method=FriedmanDifferential beta=[2,5,10]K/min finite_Ea_points={} mean_Ea={mean_ea:.2}J/mol range=[{:.2},{:.2}]J/mol E_true={E_TRUE:.2}J/mol",
            config.seed,
            finite_ea.len(),
            finite_ea.iter().copied().fold(f64::INFINITY, f64::min),
            finite_ea.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        );
        assert!(
            (mean_ea - E_TRUE).abs() / E_TRUE < 0.25,
            "Differential Friedman Ea={mean_ea}; truth={E_TRUE}"
        );
        let arrhenius = crate::Kinetics::experimental_kinetics::kinetic_methods::kinetic_regression::linear_regression(
            &ndarray::Array1::from_vec(arrhenius_x),
            &ndarray::Array1::from_vec(arrhenius_y),
        )
        .expect("fit Arrhenius parameters recovered from the non-isothermal series");
        let recovered_e = -arrhenius.slope * 8.314;
        let recovered_a = arrhenius.intercept.exp();
        println!(
            "[STORY][non-isothermal-arrhenius] seed={} rate=eta_rate/(1-eta) E_true={E_TRUE:.2}J/mol E_fit={recovered_e:.2}J/mol A_true=2.000000e5 A_fit={recovered_a:.6e} R2={:.6}",
            config.seed, arrhenius.r2
        );
        assert!(
            (recovered_e - E_TRUE).abs() / E_TRUE < 0.15,
            "rate-derived E={recovered_e}; truth={E_TRUE}"
        );
        assert!(
            (recovered_a / 2.0e5).ln().abs() < 0.7,
            "rate-derived A={recovered_a}; truth=2e5"
        );

        let clean_series =
            super::support::make_clean_export_series(&harness.state().app.model, &ids);
        harness.state_mut().app.model.series = clean_series;
        harness.state_mut().app.new_experiment_dialog.current_dir = export_dir.path().to_path_buf();
        harness.run();
        super::support::click_menu_entry(harness, "File Manager", "Save As CSV");
        harness.get_by_label("Save").click();
        harness.run();
        assert!(export_path.is_file());
        let mut roundtrip = crate::gui::experimental_kinetics_gui::model::PlotModel::new();
        roundtrip.from_csv_series(&export_path).unwrap();
        assert_eq!(roundtrip.list_of_experiments().len(), 3);
        for id in &ids {
            let columns = roundtrip.list_of_columns(&id).unwrap();
            assert_eq!(columns.len(), 5);
            for name in ["t", "T", "sg_mass", "eta", "deta_dt"] {
                assert!(columns.contains(&name.to_string()));
            }
            assert!(columns.contains(&"eta".into()));
            assert!(columns.contains(&"deta_dt".into()));
        }
    });

    let mut reloaded = crate::gui::experimental_kinetics_gui::model::PlotModel::new();
    reloaded.from_csv_series(&export_path).unwrap();
    let mut reloaded_app = empty_app();
    reloaded_app.model = reloaded;
    with_app(reloaded_app, |harness| {
        harness
            .state_mut()
            .app
            .kinetic_methods_state
            .isoconversional_method = IsoconversionalMethod::FriedmanDifferential;
        harness.run();
        super::support::click_menu_entry(harness, "Kinetic Methods", "Isoconversional");
        harness.get_by_label("Calculate").click();
        harness.run();
        assert!(
            harness
                .state()
                .app
                .kinetic_methods_state
                .isoconversional_status
                .starts_with("Done.")
        );
    });
}

#[test]
fn story_non_isothermal_method_rejects_missing_heating_rate_metadata() {
    let mut config = deterministic_config(false, 98);
    config.experiment_mode = ExperimentMode::NonIsothermal {
        t0: 420.0,
        heating_rates: vec![2.0, 5.0, 10.0],
    };
    let mut model = crate::gui::experimental_kinetics_gui::model::PlotModel::new();
    model.generate_synthetic_data_from_config(&config).unwrap();
    let mut no_beta_series =
        crate::Kinetics::experimental_kinetics::experiment_series_main::TGASeries::new();
    for (index, run) in
        crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA::generate_series(&config)
            .unwrap()
            .into_iter()
            .enumerate()
    {
        no_beta_series
            .create_from_synthetic_data(
                &run.virtual_tga,
                crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta {
                    id: format!("no_beta_{index}"),
                    heating_rate: None,
                    isothermal_temperature: None,
                    thermal_mode: crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentThermalMode::Unspecified,
                    comment: None,
                },
            )
            .unwrap();
    }
    model.series = no_beta_series;
    for id in model.list_of_experiments() {
        model.conversion_for_experiment(id, 0.005, "eta").unwrap();
    }

    let mut app = empty_app();
    app.model = model;
    with_app(app, |harness| {
        harness
            .state_mut()
            .app
            .kinetic_methods_state
            .isoconversional_method = IsoconversionalMethod::OFW;
        harness.run();
        super::support::click_menu_entry(harness, "Kinetic Methods", "Isoconversional");
        harness.get_by_label("Calculate").click();
        harness.run();
        let state = harness.state();
        assert!(
            state
                .app
                .kinetic_methods_state
                .isoconversional_status
                .contains("Error")
        );
        assert!(
            !state
                .app
                .model
                .list_of_experiments()
                .contains(&"isoconversional_result".to_string())
        );
    });
}
