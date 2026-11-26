use super::support::{empty_app, with_app};
use crate::Kinetics::experimental_kinetics::math::fitting::{FitColumnRequest, FittingModelName};
use crate::Kinetics::experimental_kinetics::testing_mod::{
    AdvancedTGAConfig, ExperimentMode, KineticModel, NoiseConfig, NoiseKind, VirtualTGA,
};
use egui_kittest::kittest::Queryable;
use tempfile::tempdir;

fn noisy_gap_csv() -> (tempfile::TempDir, std::path::PathBuf, f64, f64) {
    let cfg = AdvancedTGAConfig {
        n_points: 2400,
        dt: 1.0,
        experiment_mode: ExperimentMode::NonIsothermal {
            t0: 520.0,
            heating_rates: vec![0.02],
        },
        kinetic_model: KineticModel::ArrheniusSingle {
            m0: 100.0,
            k0: 2.0e5,
            e: 82_000.0,
            r: 8.314,
        },
        mass_noise: Some(NoiseConfig {
            kind: NoiseKind::Gaussian { sigma: 0.04 },
        }),
        temp_noise: None,
        spikes: None,
        seed: 7401,
    };
    let source = VirtualTGA::generate_series(&cfg)
        .unwrap()
        .remove(0)
        .virtual_tga;
    let (gap_start, gap_end) = (950, 1150);
    let mut csv = String::from("t,T,m\n");
    for index in 0..source.time.len() {
        if (gap_start..gap_end).contains(&index) {
            continue;
        }
        csv.push_str(&format!(
            "{:.8},{:.8},{:.8}\n",
            source.time[index], source.temperature[index], source.mass[index]
        ));
    }
    let dir = tempdir().unwrap();
    let path = dir.path().join("gapped_isothermal.csv");
    std::fs::write(&path, csv).unwrap();
    (dir, path, gap_start as f64, gap_end as f64)
}

#[test]
fn story_fit_both_sides_of_missing_interval_and_commit_aligned_reconstruction() {
    let (dir, path, gap_start, gap_end) = noisy_gap_csv();
    let mut app = empty_app();
    // This story tests a real isothermal acquisition with known metadata.
    app.new_experiment_dialog.current_dir = dir.path().to_path_buf();
    with_app(app, |harness| {
        super::support::import_csv(harness, &path, "missing_interval", None, Some(0.02));
        let id = "missing_interval";

        let sample_count = harness.state().app.model.settings.n_points();
        let source_before_failed_fit = harness
            .state()
            .app
            .model
            .series
            .get_experiment_by_id(id)
            .unwrap();
        let history_before_failed_fit = source_before_failed_fit
            .dataset
            .history_of_operations
            .feed_lines();
        let undo_before_failed_fit = source_before_failed_fit.dataset.undo_depth();
        let plots_before_failed_fit = harness
            .state()
            .app
            .model
            .plots
            .iter()
            .map(|plot| {
                (
                    plot.plot_short_name.clone(),
                    plot.points.clone(),
                    plot.shown,
                    plot.selected,
                )
            })
            .collect::<Vec<_>>();
        harness.state_mut().app.model.settings.set_n_points(None);
        let failed_plot_fit = harness.state_mut().app.model.fit_column_and_plot(
            FitColumnRequest::new(id, "t", "T", FittingModelName::Linear)
                .with_output_col("unplotted_fit"),
        );
        assert!(failed_plot_fit.is_err());
        let unchanged = harness
            .state()
            .app
            .model
            .series
            .get_experiment_by_id(id)
            .unwrap();
        assert!(
            !unchanged
                .dataset
                .schema
                .columns
                .contains_key("unplotted_fit")
        );
        assert_eq!(
            unchanged.dataset.history_of_operations.feed_lines(),
            history_before_failed_fit
        );
        assert_eq!(unchanged.dataset.undo_depth(), undo_before_failed_fit);
        let plots_after_failed_fit = harness
            .state()
            .app
            .model
            .plots
            .iter()
            .map(|plot| {
                (
                    plot.plot_short_name.clone(),
                    plot.points.clone(),
                    plot.shown,
                    plot.selected,
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(plots_after_failed_fit, plots_before_failed_fit);
        harness
            .state_mut()
            .app
            .model
            .settings
            .set_n_points(sample_count);

        let source_rows = super::support::read_column(&harness.state().app.model, id, "t").len();
        let source_times = super::support::read_column(&harness.state().app.model, id, "t");

        super::support::click_menu_entry(harness, "Kinetic Methods", "Fit Model");
        {
            let fit = &mut harness.state_mut().app.kinetic_methods_state;
            fit.fit_model_experiment_id = id.into();
            fit.fit_model_x_col = "t".into();
            fit.fit_model_y_col = "T".into();
            fit.fit_model_output_col = "temperature_rebuilt".into();
            fit.fit_model_model =
                crate::Kinetics::experimental_kinetics::math::fitting::FittingModelName::Linear;
            fit.fit_gap_enabled = true;
            fit.fit_gap_start = gap_start;
            fit.fit_gap_end = gap_end;
            fit.fit_evaluation_points = 200;
        }
        harness.run();
        harness.get_by_label("Fit").click();
        harness.run();
        let temperature_r2 = harness
            .state()
            .app
            .kinetic_methods_state
            .fit_model_result
            .as_ref()
            .unwrap()
            .r2;
        assert!(temperature_r2 > 0.9999);

        {
            let fit = &mut harness.state_mut().app.kinetic_methods_state;
            fit.fit_model_y_col = "m".into();
            fit.fit_model_output_col = "mass_rebuilt".into();
            fit.fit_model_model =
                crate::Kinetics::experimental_kinetics::math::fitting::FittingModelName::DecExp;
        }
        harness.run();
        harness.get_by_label("Fit").click();
        harness.run();
        let mass_fit = harness
            .state()
            .app
            .kinetic_methods_state
            .fit_model_result
            .as_ref()
            .unwrap();
        assert!(mass_fit.r2 > 0.995, "mass fit R2={}", mass_fit.r2);
        println!(
            "[STORY][gap-reconstruction] seed=7401 noise=Gaussian(sigma=0.04mg) gap=t[{gap_start:.0},{gap_end:.0})s missing_rows=200 observed_rows={source_rows} eval_points=200 shared_x_rows={} T_model=Linear R2_T={temperature_r2:.7} mass_model=DecExp R2_mass={:.7}",
            source_rows + 200,
            mass_fit.r2
        );
        assert_eq!(
            harness.state().app.kinetic_methods_state.fit_sessions.len(),
            2
        );

        harness
            .get_by_label("Evaluate fitted columns")
            .scroll_to_me();
        harness.run();
        harness.get_by_label("Evaluate fitted columns").click();
        harness.run();
        let batch = harness
            .state()
            .app
            .kinetic_methods_state
            .fit_pending_batch
            .as_ref()
            .unwrap();
        assert_eq!(batch.columns.len(), 2);
        assert_eq!(batch.evaluation_x.len(), 200);
        assert!(
            batch
                .columns
                .iter()
                .all(|column| column.values.len() == batch.x_values.len())
        );
        assert_eq!(batch.x_values.len(), source_rows + 200);
        assert_eq!(
            super::support::read_column(&harness.state().app.model, id, "t").len(),
            source_rows,
            "fit/evaluation must not mutate the acquired source"
        );

        let previous_sample_count = harness.state().app.model.settings.n_points();
        let plots_before_failed_commit = harness
            .state()
            .app
            .model
            .plots
            .iter()
            .map(|plot| {
                (
                    plot.plot_short_name.clone(),
                    plot.points.clone(),
                    plot.shown,
                    plot.selected,
                )
            })
            .collect::<Vec<_>>();
        let mut stale_reconstruction_curve =
            harness.state().app.model.plots.first().unwrap().clone();
        stale_reconstruction_curve.experiment_id = "failed_plot_commit".into();
        stale_reconstruction_curve.experiment_index = 99;
        stale_reconstruction_curve.x_name = "t".into();
        stale_reconstruction_curve.y_name = "temperature_rebuilt".into();
        stale_reconstruction_curve.plot_short_name =
            "failed_plot_commit_t_temperature_rebuilt".into();
        stale_reconstruction_curve.points = vec![[0.0, 1.0], [1.0, 2.0]];
        stale_reconstruction_curve.shown = false;
        harness
            .state_mut()
            .app
            .model
            .plots
            .insert(1, stale_reconstruction_curve);
        let previous_view = (
            harness.state().app.model.interaction.view_range,
            harness.state().app.model.interaction.view_y_range,
        );
        harness.state_mut().app.model.settings.set_n_points(None);
        let failed_commit = {
            let state = harness.state_mut();
            let batch = state
                .app
                .kinetic_methods_state
                .fit_pending_batch
                .as_ref()
                .unwrap();
            state
                .app
                .model
                .push_fitted_evaluation_batch_and_plot(batch, "failed_plot_commit")
        };
        assert!(matches!(
            failed_commit,
            Err(crate::gui::experimental_kinetics_gui::model::TGAGUIError::PostCommitWarning(_))
        ));
        assert_eq!(harness.state().app.model.series.len(), 2);
        assert!(
            harness
                .state()
                .app
                .model
                .series
                .get_experiment_by_id("failed_plot_commit")
                .is_ok()
        );
        assert!(
            harness
                .state()
                .app
                .model
                .plots
                .iter()
                .all(|plot| plot.experiment_id != "failed_plot_commit")
        );
        assert!(
            harness
                .state()
                .app
                .model
                .series
                .get_experiment_by_id("missing_interval")
                .is_ok()
        );
        assert_eq!(
            (
                harness.state().app.model.interaction.view_range,
                harness.state().app.model.interaction.view_y_range,
            ),
            previous_view
        );
        let remaining_plots = harness
            .state()
            .app
            .model
            .plots
            .iter()
            .map(|plot| {
                (
                    plot.plot_short_name.clone(),
                    plot.points.clone(),
                    plot.shown,
                    plot.selected,
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(remaining_plots, plots_before_failed_commit);
        assert!(
            harness
                .state()
                .app
                .kinetic_methods_state
                .fit_pending_batch
                .is_some()
        );
        harness
            .state_mut()
            .app
            .model
            .settings
            .set_n_points(previous_sample_count);

        harness
            .state_mut()
            .app
            .new_experiment_dialog
            .column_manager_open = true;
        harness.run();
        harness.get_by_label("Add reconstructed experiment").click();
        harness.run();
        let reconstructed = harness
            .state()
            .app
            .model
            .series
            .get_experiment_by_id("missing_interval_reconstructed")
            .unwrap();
        let frame = reconstructed.dataset.frame.clone().collect().unwrap();
        assert_eq!(frame.height(), source_rows + 200);
        assert!(
            frame
                .columns()
                .iter()
                .all(|column| column.null_count() == 0)
        );
        let rebuilt_time = super::support::read_column(
            &harness.state().app.model,
            "missing_interval_reconstructed",
            "t",
        );
        assert!(rebuilt_time.windows(2).all(|pair| pair[1] > pair[0]));
        for rebuilt in ["temperature_rebuilt", "mass_rebuilt"] {
            let values = super::support::read_column(
                &harness.state().app.model,
                "missing_interval_reconstructed",
                rebuilt,
            );
            assert_eq!(values.len(), rebuilt_time.len());
            assert!(
                reconstructed
                    .dataset
                    .column_provenance_text(rebuilt)
                    .is_some()
            );
        }
        assert_eq!(source_times.len(), source_rows);
        assert_eq!(harness.state().app.model.series.len(), 3);
    });
}

#[test]
fn story_unsuitable_fit_reports_low_quality_without_touching_source() {
    let (dir, path, gap_start, gap_end) = noisy_gap_csv();
    let mut app = empty_app();
    app.new_experiment_dialog.current_dir = dir.path().to_path_buf();
    with_app(app, |harness| {
        super::support::import_csv(harness, &path, "poor_fit", Some(520.0), None);
        super::support::click_menu_entry(harness, "Kinetic Methods", "Fit Model");
        {
            let state = &mut harness.state_mut().app.kinetic_methods_state;
            state.fit_model_experiment_id = "poor_fit".into();
            state.fit_model_x_col = "t".into();
            state.fit_model_y_col = "m".into();
            state.fit_model_output_col = "bad_fit".into();
            state.fit_model_model =
                crate::Kinetics::experimental_kinetics::math::fitting::FittingModelName::Linear;
            state.fit_gap_enabled = true;
            state.fit_gap_start = gap_start;
            state.fit_gap_end = gap_end;
        }
        harness.run();
        let before = super::support::read_column(&harness.state().app.model, "poor_fit", "m");
        harness.get_by_label("Fit").click();
        harness.run();
        let state = harness.state();
        let fit = state
            .app
            .kinetic_methods_state
            .fit_model_result
            .as_ref()
            .unwrap();
        assert!(
            fit.r2 < 0.95,
            "an unsuitable linear mass fit unexpectedly passed: R2={}",
            fit.r2
        );
        println!(
            "[STORY][gap-reconstruction-negative] seed=7401 model=Linear applied_to=exponential_mass R2={:.7} low_quality=true decision_left_to_user=true source_rows={} measured_values_unchanged=true",
            fit.r2,
            before.len()
        );
        assert!(
            state
                .app
                .kinetic_methods_state
                .fit_model_status
                .contains("R2=")
        );
        assert!(
            state
                .app
                .model
                .series
                .get_experiment_by_id("poor_fit")
                .unwrap()
                .dataset
                .schema
                .columns
                .contains_key("bad_fit")
        );
        assert_eq!(
            super::support::read_column(&state.app.model, "poor_fit", "m"),
            before
        );
    });
}
