use super::support::{
    click_menu_entry, deterministic_config, empty_app, generated, import_csv, read_column,
    with_app, write_simulated_csv,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::isoconversion::IsoconversionalMethod;
use crate::Kinetics::experimental_kinetics::math::fitting::FittingModelName;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
use crate::Kinetics::experimental_kinetics::testing_mod::ExperimentMode;
use eframe::egui;
use egui_kittest::kittest::Queryable;
use ndarray::Array1;
use tempfile::tempdir;

fn fit_k_for_experiment(
    harness: &mut egui_kittest::Harness<super::support::StoryState>,
    id: &str,
) -> f64 {
    if !harness.state().app.kinetic_methods_state.fit_model_open {
        super::support::click_menu_entry(harness, "Kinetic Methods", "Fit Model");
    }
    {
        let fit = &mut harness.state_mut().app.kinetic_methods_state;
        fit.fit_model_experiment_id = id.to_string();
        fit.fit_model_x_col = "t".to_string();
        fit.fit_model_y_col = "sg_mass".to_string();
        fit.fit_model_output_col = format!("fit_{id}");
        fit.fit_model_model = FittingModelName::DecExp;
        fit.fit_model_initial_guess.clear();
    }
    harness.run();
    harness.get_by_label("Fit").click();
    harness.run();
    let result = harness
        .state()
        .app
        .kinetic_methods_state
        .fit_model_result
        .as_ref()
        .expect("GUI did not produce a DecExp fit");
    assert!(result.r2 > 0.98, "{} fit R2={}", id, result.r2);
    let k = result
        .coefficients
        .iter()
        .find_map(|(name, value)| (name == "k").then_some(*value))
        .expect("DecExp fit did not report k");
    println!(
        "[STORY][isothermal-fit] experiment={id} model=DecExp k_per_hour={k:.6} R2={:.6}",
        result.r2
    );
    k
}

#[test]
fn story_three_noisy_isothermal_runs_recover_arrhenius_parameters_and_export() {
    const E_TRUE: f64 = 85_000.0;
    const A_TRUE: f64 = 2.0e5;
    const R: f64 = 8.314;
    let mut config = deterministic_config(true, 91_743);
    config.mass_noise = Some(
        crate::Kinetics::experimental_kinetics::testing_mod::NoiseConfig {
            kind: crate::Kinetics::experimental_kinetics::testing_mod::NoiseKind::Gaussian {
                sigma: 0.12,
            },
        },
    );
    config.temp_noise = None;
    let generated_runs = generated(&config);
    let mut truth_config = config.clone();
    truth_config.mass_noise = None;
    let truth_runs = generated(&truth_config);
    let dir = tempdir().unwrap();
    let files = generated_runs
        .iter()
        .map(|run| {
            let temp = run.meta.isothermal_temperature.unwrap();
            let path = dir.path().join(format!("isothermal_{temp:.0}.csv"));
            write_simulated_csv(&path, run, true);
            (path, temp)
        })
        .collect::<Vec<_>>();

    let mut app = empty_app();
    app.new_experiment_dialog.current_dir = dir.path().to_path_buf();
    let export_dir = tempdir().unwrap();
    let export_path = export_dir.path().join("series_export.csv");
    with_app(app, |harness| {
        for (path, temperature) in &files {
            let id = format!("iso_{temperature:.0}");
            import_csv(harness, path, &id, Some(*temperature), None);
        }

        let ids = harness.state().app.model.list_of_experiments();
        assert_eq!(ids.len(), 3);
        let mut estimated_k = Vec::new();
        let mut reciprocal_temperature = Vec::new();
        for ((id, temperature), truth) in ids.iter().zip([500.0, 530.0, 560.0]).zip(&truth_runs) {
            let model = &mut harness.state_mut().app.model;
            model.settings.set_calibration_coeffs(-1.0, 40.0);
            model
                .series
                .bind_time(id, "t", Unit::Second)
                .expect("bind seconds");
            model
                .series
                .bind_temperature(id, "T", Unit::Kelvin)
                .expect("bind kelvin");
            model
                .series
                .bind_mass(id, "m", Unit::MilliVolt)
                .expect("bind raw voltage");
            model
                .calibrate_mass_from_voltage_new_column(id, "mass_mg")
                .expect("calibrate raw voltage to mass");
            model
                .series
                .bind_mass(id, "mass_mg", Unit::Milligram)
                .expect("bind calibrated mass");
            let quality = {
                model
                    .series
                    .sg_filter_column_as(id, "mass_mg", 31, 3, 0, 1.0, Some("sg_mass"))
                    .expect("smooth calibrated mass");
                model
                    .series
                    .get_experiment_by_id(id)
                    .unwrap()
                    .dataset
                    .filter_quality_report("mass_mg", "sg_mass")
                    .expect("compute smoothing quality")
            };
            assert!(
                quality.roughness_ratio_1st < 1.0 && quality.roughness_ratio_2nd < 1.0,
                "SG did not reduce high-frequency noise: {quality}"
            );
            let raw_mass = read_column(model, id, "mass_mg");
            let smoothed_mass = read_column(model, id, "sg_mass");
            let clean_mass = &truth.virtual_tga.mass;
            let raw_truth_rmse = super::support::rmse(&raw_mass, clean_mass);
            let smooth_truth_rmse = super::support::rmse(&smoothed_mass, clean_mass);
            assert!(
                smooth_truth_rmse < raw_truth_rmse,
                "SG increased mass error against the known noise-free signal: raw={raw_truth_rmse}, filtered={smooth_truth_rmse}"
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
                .expect("use smoothed mass as active source");
            model.from_s_to_h_for_experiment(id).unwrap();
            model
                .conversion_for_experiment(id.to_string(), 0.01, "eta")
                .expect("derive conversion");
            let raw_eta = read_column(model, id, "eta");
            model
                .series
                .sg_filter_column_as(id, "eta", 31, 3, 0, 1.0, Some("eta"))
                .expect("smooth conversion before differentiation");
            let smoothed_eta = read_column(model, id, "eta");
            let eta_report = super::support::quality_from_slices(
                &raw_eta,
                &smoothed_eta,
                "eta_before_sg",
                "eta_after_sg",
            );
            let truth_eta = super::support::ideal_conversion(
                &truth.virtual_tga.time,
                &truth.virtual_tga.mass,
                0.01 * 3600.0,
            );
            let raw_eta_truth_rmse = super::support::rmse(&raw_eta, &truth_eta);
            let smooth_eta_truth_rmse = super::support::rmse(&smoothed_eta, &truth_eta);
            super::support::print_filter_report(
                config.seed,
                id,
                "conversion",
                31,
                3,
                &eta_report,
                Some((raw_eta_truth_rmse, smooth_eta_truth_rmse)),
            );
            model
                .derive_deta_dt(id, Some("deta_dt"))
                .expect("derive conversion rate");
            model
                .series
                .drop_nulls(id)
                .expect("remove derivative edge rows together");
            let raw_rate = read_column(model, id, "deta_dt");
            model
                .series
                .sg_filter_column_as(id, "deta_dt", 11, 3, 0, 1.0, Some("deta_dt"))
                .expect("clean rate noise");
            let smooth_rate = read_column(model, id, "deta_dt");
            let rate_report = super::support::quality_from_slices(
                &raw_rate,
                &smooth_rate,
                "deta_dt_before_sg",
                "deta_dt_after_sg",
            );
            assert!(
                rate_report.roughness_ratio_1st < 1.0 && rate_report.roughness_ratio_2nd < 1.0,
                "final SG did not reduce rate roughness: {rate_report}"
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

            let dataset = &model.series.get_experiment_by_id(id).unwrap().dataset;
            assert_eq!(dataset.schema.columns["t"].unit, Unit::Hour);
            assert_eq!(dataset.schema.columns["mass_mg"].unit, Unit::Milligram);
            assert!(dataset.schema.columns.contains_key("deta_dt"));
            assert!(dataset.column_provenance_text("deta_dt").is_some());
            for name in ["t", "T", "mass_mg", "sg_mass", "eta", "deta_dt"] {
                assert_eq!(
                    dataset
                        .frame
                        .clone()
                        .collect()
                        .unwrap()
                        .column(name)
                        .unwrap()
                        .null_count(),
                    0,
                    "unexpected null in {id}.{name}"
                );
            }
            let raw_times = read_column(model, id, "t");
            assert!(raw_times.last().unwrap() > &0.3);
            estimated_k.push(fit_k_for_experiment(harness, id) / 3600.0);
            reciprocal_temperature.push(1.0 / temperature);
        }

        let fit = crate::Kinetics::experimental_kinetics::kinetic_methods::kinetic_regression::linear_regression(
            &Array1::from_vec(reciprocal_temperature.clone()),
            &Array1::from_vec(estimated_k.iter().map(|k| k.ln()).collect()),
        )
        .expect("fit Arrhenius parameters recovered from the isothermal series");
        let e_estimated = -fit.slope * R;
        let a_estimated = fit.intercept.exp();
        println!(
            "[STORY][isothermal-arrhenius] seed={} model=ArrheniusSingle E_true={E_TRUE:.2}J/mol E_fit={e_estimated:.2}J/mol A_true={A_TRUE:.6e} A_fit={a_estimated:.6e} R2={:.6}",
            config.seed, fit.r2
        );
        assert!(
            (e_estimated - E_TRUE).abs() / E_TRUE < 0.12,
            "E={e_estimated} J/mol; truth={E_TRUE}"
        );
        assert!(
            (a_estimated / A_TRUE).ln().abs() < 0.65,
            "A={a_estimated}; truth={A_TRUE}"
        );

        // Run the app's explicit isothermal method against the prepared series.
        {
            let kinetic = &mut harness.state_mut().app.kinetic_methods_state;
            kinetic.isoconversional_method = IsoconversionalMethod::FriedmanIntegral;
            kinetic.isoconversional_output_id = "isothermal_friedman_result".into();
        }
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
            "Friedman Integral failed: {status}"
        );
        harness
            .state()
            .app
            .model
            .series
            .get_experiment_by_id("isothermal_friedman_result")
            .expect("GUI kinetic method did not publish its result");
        let ea = read_column(
            &harness.state().app.model,
            "isothermal_friedman_result",
            "Ea",
        );
        let result_eta = read_column(
            &harness.state().app.model,
            "isothermal_friedman_result",
            "eta",
        );
        let finite_ea = ea
            .into_iter()
            .zip(result_eta)
            .filter_map(|(value, eta)| {
                (value.is_finite() && (0.1..=0.9).contains(&eta)).then_some(value)
            })
            .collect::<Vec<_>>();
        assert!(!finite_ea.is_empty());
        let mean_ea = finite_ea.iter().sum::<f64>() / finite_ea.len() as f64;
        assert!(
            (mean_ea - E_TRUE).abs() / E_TRUE < 0.12,
            "Friedman Integral mean Ea={mean_ea}; truth={E_TRUE}"
        );
        println!(
            "[STORY][isothermal-method] method=FriedmanIntegral eta_range=[0.1,0.9] finite_Ea_points={} mean_Ea={mean_ea:.2}J/mol range=[{:.2},{:.2}]J/mol",
            finite_ea.len(),
            finite_ea.iter().copied().fold(f64::INFINITY, f64::min),
            finite_ea.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        );

        let clean_series =
            super::support::make_clean_export_series(&harness.state().app.model, &ids);
        harness.state_mut().app.model.series = clean_series;
        harness.state_mut().app.new_experiment_dialog.current_dir = export_dir.path().to_path_buf();
        harness.run();
        super::support::click_menu_entry(harness, "File Manager", "Save As CSV");
        harness.get_by_label("Save").click();
        harness.run();
        assert!(
            export_path.is_file(),
            "CSV export did not create the expected file"
        );
        let mut roundtrip = crate::gui::experimental_kinetics_gui::model::PlotModel::new();
        roundtrip.from_csv_series(&export_path).unwrap();
        assert_eq!(roundtrip.list_of_experiments().len(), 3);
        for id in &ids {
            let columns = roundtrip.list_of_columns(id).unwrap();
            assert_eq!(columns.len(), 5);
            for name in ["t", "T", "sg_mass", "eta", "deta_dt"] {
                assert!(columns.contains(&name.to_string()));
            }
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
            .isoconversional_method = IsoconversionalMethod::FriedmanIntegral;
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
fn story_isothermal_golden_pipeline_runs_without_measured_temperature() {
    let config = deterministic_config(true, 29_104);
    let run = generated(&config).remove(0);
    let temperature = run.meta.isothermal_temperature.unwrap();
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("isothermal_t_m.csv");
    let mut csv = String::from("t,m\n");
    for (time, mass) in run.virtual_tga.time.iter().zip(&run.virtual_tga.mass) {
        csv.push_str(&format!("{time:.8},{mass:.8}\n"));
    }
    std::fs::write(&input_path, csv).unwrap();

    let mut app = empty_app();
    app.new_experiment_dialog.current_dir = dir.path().to_path_buf();
    app.model.settings.set_calibration_coeffs(1.0, 0.0);
    with_app(app, |harness| {
        import_csv(
            harness,
            &input_path,
            "isothermal_t_m",
            Some(temperature),
            None,
        );
        harness.set_size(eframe::egui::vec2(1300.0, 1000.0));
        harness.run();
        harness.get_by_label("golden pipeline").click();
        harness.run();
        {
            let state = &mut harness
                .state_mut()
                .app
                .quick_actions_state
                .golden_pipeline_dialog;
            state.config.isothermal_temperature = Some(temperature);
            state.config.sav_gol_config.window_size = 11;
            state.config.spline_config.n_points = 240;
            state.config.spline_config.n_internal_points = 12;
            state.config.save_to_new_experiment = true;
            state.config.del_old_experiment = true;
        }
        harness.run();
        harness.get_by_label("Calculate").click();
        harness.run();

        let model = &harness.state().app.model;
        assert!(
            harness
                .state()
                .app
                .quick_actions_state
                .golden_pipeline_dialog
                .last_error
                .is_none(),
            "Golden Pipeline failed: {:?}; status: {}",
            harness
                .state()
                .app
                .quick_actions_state
                .golden_pipeline_dialog
                .last_error,
            model.message
        );
        let output = model
            .series
            .get_experiment_by_id("proceeded_isothermal_t_m")
            .unwrap_or_else(|error| {
                panic!(
                    "Golden Pipeline should publish its output; ids={:?}: {error}",
                    model
                        .series
                        .experiments()
                        .iter()
                        .map(|exp| exp.meta.id.as_str())
                        .collect::<Vec<_>>()
                )
            });
        assert!(output.dataset.schema.temperature.is_none());
        assert!(output.dataset.schema.dT_dt.is_none());
        assert!(output.dataset.schema.time.is_some());
        assert!(output.dataset.schema.mass.is_some());
        assert_eq!(output.meta.resolved_thermal_mode(), crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentThermalMode::Isothermal);
        assert!(model.message.contains("Golden Pipeline completed"));
        assert!(
            harness
                .state()
                .app
                .quick_actions_state
                .golden_pipeline_dialog
                .quality_report
                .is_some()
        );
        assert!(model.series.get_experiment_by_id("isothermal_t_m").is_err());
    });
}

#[test]
fn gui_import_accepts_isothermal_t_m_file_without_measured_temperature() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("isothermal_t_m.csv");
    std::fs::write(&path, "t,m\n0.0,10.0\n1.0,9.8\n2.0,9.5\n").unwrap();
    let mut app = empty_app();
    app.new_experiment_dialog.current_dir = dir.path().to_path_buf();

    with_app(app, |harness| {
        import_csv(harness, &path, "iso_no_sensor", Some(523.15), None);
        let experiment = harness
            .state()
            .app
            .model
            .series
            .get_experiment_by_id("iso_no_sensor")
            .unwrap();
        assert_eq!(experiment.meta.resolved_thermal_mode(), crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentThermalMode::Isothermal);
        assert_eq!(experiment.meta.isothermal_temperature, Some(523.15));
        assert!(experiment.dataset.schema.temperature.is_none());
        assert_eq!(
            harness
                .state()
                .app
                .model
                .plots
                .iter()
                .filter(|curve| curve.experiment_id == "iso_no_sensor")
                .count(),
            1
        );
    });
}

#[test]
fn gui_import_can_adopt_mean_of_measured_celsius_temperature_as_isothermal_setpoint() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("isothermal_sensor.csv");
    std::fs::write(&path, "t,T,m\n0.0,25.0,10.0\n1.0,26.0,9.8\n2.0,27.0,9.5\n").unwrap();
    let mut app = empty_app();
    app.new_experiment_dialog.current_dir = dir.path().to_path_buf();

    with_app(app, |harness| {
        click_menu_entry(harness, "File Manager", "Import Data");
        harness.get_by_label("[F] isothermal_sensor.csv").click();
        harness.run();
        let id_input = harness.get_by_label("Experiment ID:");
        id_input.focus();
        harness.key_press_modifiers(egui::Modifiers::COMMAND, egui::Key::A);
        harness.key_press(egui::Key::Backspace);
        id_input.type_text("iso_sensor_mean");
        harness.get_by_label("Isothermal").click();
        harness.run();
        harness
            .get_by_label("Use mean of the bound measured-temperature column")
            .click();
        harness.run();
        harness.get_by_label("Import experiment").click();
        harness.run();

        assert!(
            harness
                .state()
                .app
                .new_experiment_dialog
                .last_error
                .is_none(),
            "GUI import failed: {:?}",
            harness.state().app.new_experiment_dialog.last_error
        );
        let experiment = harness
            .state()
            .app
            .model
            .series
            .get_experiment_by_id("iso_sensor_mean")
            .unwrap();
        assert_eq!(experiment.meta.thermal_mode, crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentThermalMode::Isothermal);
        assert!((experiment.meta.isothermal_temperature.unwrap() - 299.15).abs() < 1e-10);
        assert_eq!(experiment.dataset.schema.temperature.as_deref(), Some("T"));
        assert!(harness.state().app.model.message.contains("299.150 K"));
    });
}

#[test]
fn story_isothermal_method_rejects_series_without_isothermal_metadata() {
    let mut config = deterministic_config(false, 771);
    config.experiment_mode = ExperimentMode::NonIsothermal {
        t0: 300.0,
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
            .isoconversional_method = IsoconversionalMethod::FriedmanIntegral;
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
                .contains("Error")
        );
        assert!(
            !harness
                .state()
                .app
                .model
                .list_of_experiments()
                .contains(&"isoconversional_result".to_string())
        );
    });
}

#[test]
fn story_nonisothermal_method_rejects_declared_isothermal_series() {
    let config = deterministic_config(true, 817);
    let mut model = crate::gui::experimental_kinetics_gui::model::PlotModel::new();
    model.generate_synthetic_data_from_config(&config).unwrap();
    let mut app = empty_app();
    app.model = model;

    with_app(app, |harness| {
        harness
            .state_mut()
            .app
            .kinetic_methods_state
            .isoconversional_method = IsoconversionalMethod::KAS;
        harness.run();
        super::support::click_menu_entry(harness, "Kinetic Methods", "Isoconversional");
        harness.get_by_label("Calculate").click();
        harness.run();
        let status = &harness
            .state()
            .app
            .kinetic_methods_state
            .isoconversional_status;
        assert!(status.contains("isothermal"), "unexpected status: {status}");
        assert_eq!(harness.state().app.model.series.len(), 3);
    });
}
