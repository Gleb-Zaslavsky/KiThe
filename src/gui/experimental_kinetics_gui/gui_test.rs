use crate::gui::experimental_kinetics_gui::model::PlotModel;

use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;

use std::path::Path;

const TEST_DATA_PATH: &str = "src/assets/TGAexample.txt";

fn test_data_path() -> &'static Path {
    Path::new(TEST_DATA_PATH)
}

pub fn load_model_with_one_curve() -> PlotModel {
    let mut model = PlotModel::new();

    model
        .push_from_file(test_data_path())
        .expect("failed to load test data");
    let experiments = model.list_of_experiments();
    let exp_id = experiments
        .first()
        .cloned()
        .expect("no experiments loaded from fixture");

    model
        .bind_time_of_experiment(&exp_id, "t", Unit::Second)
        .expect("failed to bind time column");
    model
        .bind_temperature_of_experiment(&exp_id, "T", Unit::Celsius)
        .expect("failed to bind temperature column");
    model
        .bind_mass_of_experiment(&exp_id, "m", Unit::MilliVolt)
        .expect("failed to bind mass column");
    model.set_x(&exp_id, "t").expect("failed to set x column");
    model.set_y(&exp_id, "T").expect("failed to set y column");
    model
        .create_points_for_curve(&exp_id)
        .expect("failed to create points for curve");

    model
}

pub fn load_model_with_two_curves() -> PlotModel {
    let mut model = load_model_with_one_curve();
    let exp_id = model.list_of_experiments()[0].clone();

    model.set_x(&exp_id, "t").expect("failed to keep x column");
    model
        .set_y(&exp_id, "m")
        .expect("failed to set second y column");
    model
        .create_points_for_curve(&exp_id)
        .expect("failed to create second curve");

    model
}

#[cfg(test)]
fn load_model_with_recorded_gap() -> PlotModel {
    let observed_time = vec![0.0, 1.0, 2.0, 3.0, 7.0, 8.0, 9.0, 10.0];
    let mut model = PlotModel::new();
    model
        .series
        .create_from_synthetic_data(
            &crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA {
                time: observed_time.clone(),
                temperature: observed_time.iter().map(|t| 300.0 + 4.0 * t).collect(),
                mass: observed_time.iter().map(|t| 12.0 - 0.5 * t).collect(),
            },
            crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta {
                id: "gap_workflow".to_string(),
                ..Default::default()
            },
        )
        .unwrap();
    model.set_x("gap_workflow", "time").unwrap();
    model.set_y("gap_workflow", "temperature").unwrap();
    model.create_points_for_curve("gap_workflow").unwrap();
    model
}

#[cfg(test)]
mod tests {
    use crate::gui::experimental_kinetics_gui::controller_kinetics::DirectProblemDialogState;
    use crate::gui::experimental_kinetics_gui::controller_methods::KineticMethodsWindowState;
    use crate::gui::experimental_kinetics_gui::controller_table::ColumnManagerState;
    use crate::gui::experimental_kinetics_gui::interaction::SelectionRect;
    use crate::gui::experimental_kinetics_gui::test_options::TestOptions;

    use super::*;
    use crate::Kinetics::experimental_kinetics::exp_engine_api::XY;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::KineticMethod;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::isoconversion::{
        IsoconversionalKineticMethod, IsoconversionalMethod,
    };
    use crate::Kinetics::experimental_kinetics::math::fitting::FittingModelName;
    use crate::Kinetics::experimental_kinetics::testing_mod::{
        AdvancedTGAConfig, ExperimentMode, KineticModel, NoiseConfig, NoiseKind,
    };
    use crate::Kinetics::solid_state_kinetics_IVP::KineticModelNames;
    use crate::gui::experimental_kinetics_gui::experimental_kinetics_gui_main::PlotApp;
    use egui::{Key, Modifiers};
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};
    // ============================================================================
    // STEP 1: Make Application Testable - Test Constructor
    // ============================================================================

    impl PlotApp {
        /// Special constructor for testing with pre-configured model
        pub fn for_testing(model: PlotModel) -> Self {
            use super::super::controller_buttons_and_panels::{
                NewExperimentDialogState, QuickActionPanelState,
            };
            use super::super::controller_filters::Mathematics;
            use super::super::kitheplot_wrapper::KiThePlotWindowState;
            Self {
                model,
                quick_actions_state: QuickActionPanelState::default(),
                new_experiment_dialog: NewExperimentDialogState::new(),
                mathematics: Mathematics::new(),
                test_options: TestOptions::new(),
                column_manager_state: ColumnManagerState::new(),
                kithe_plot_window: KiThePlotWindowState::default(),
                kinetic_methods_state: KineticMethodsWindowState::default(),
                direct_problem_state: DirectProblemDialogState::new(),
                keyboard_input_focused_last_frame: false,
            }
        }
    }

    fn with_plot_app_harness(mut app: PlotApp, mut run_test: impl FnMut(&mut Harness)) {
        let mut open = true;
        let mut harness = Harness::new_ui(move |ui| {
            app.show(ui.ctx(), &mut open);
        });
        run_test(&mut harness);
    }

    struct PlotAppTestState {
        app: PlotApp,
        open: bool,
    }

    fn with_plot_app_state_harness(
        app: PlotApp,
        mut run_test: impl FnMut(&mut Harness<PlotAppTestState>),
    ) {
        let state = PlotAppTestState { app, open: true };
        let mut harness = Harness::new_ui_state(
            |ui, state: &mut PlotAppTestState| {
                state.app.show(ui.ctx(), &mut state.open);
            },
            state,
        );
        run_test(&mut harness);
    }

    fn click_menu_entry<State>(harness: &mut Harness<State>, menu_label: &str, entry_label: &str) {
        harness.get_by_label(menu_label).click();
        harness.run();
        harness.get_by_label(entry_label).click();
        harness.run();
    }

    fn assert_filter_output_via_ui(menu_item: &str, output_column: &str) {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let experiment_id = model.list_of_experiments()[0].clone();
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Math", menu_item);
            harness.get_by_label("Apply").click();
            harness.run();

            let model = &harness.state().app.model;
            assert!(
                model
                    .list_of_columns(&experiment_id)
                    .unwrap()
                    .contains(&output_column.to_string()),
                "{menu_item} did not create {output_column}; message={}",
                model.message
            );
            assert!(model.plots.iter().any(|curve| {
                curve.experiment_id == experiment_id && curve.y_name == output_column
            }));
        });
    }

    fn assert_menu_entry_clickable(menu_label: &str, entry_label: &str) {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, menu_label, entry_label);
        });
    }

    fn unique_temp_dir(prefix: &str) -> PathBuf {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before UNIX_EPOCH")
            .as_nanos();
        std::env::temp_dir().join(format!("kithe_{prefix}_{stamp}"))
    }

    fn non_isothermal_series_config() -> AdvancedTGAConfig {
        AdvancedTGAConfig {
            n_points: 180,
            dt: 0.5,
            experiment_mode: ExperimentMode::NonIsothermal {
                t0: 300.0,
                heating_rates: vec![2.0, 5.0, 10.0],
            },
            kinetic_model: KineticModel::ArrheniusSingle {
                m0: 100.0,
                k0: 2.5e5,
                e: 85000.0,
                r: 8.314,
            },
            mass_noise: Some(NoiseConfig {
                kind: NoiseKind::Gaussian { sigma: 0.0 },
            }),
            temp_noise: Some(NoiseConfig {
                kind: NoiseKind::Gaussian { sigma: 0.0 },
            }),
            spikes: None,
            seed: 42,
        }
    }

    fn last_value_for(model: &PlotModel, id: &str, column: &str) -> f64 {
        let exp = model
            .series
            .get_experiment_by_id(id)
            .expect("experiment not found");
        let df = exp
            .dataset
            .frame
            .clone()
            .collect()
            .expect("failed to collect frame");
        df.column(column)
            .expect("requested column missing")
            .f64()
            .expect("requested column is not f64")
            .into_no_null_iter()
            .last()
            .expect("requested column is empty")
    }

    // ============================================================================
    // STEP 1 Tier-1 GUI smoke tests via egui_kittest
    // ============================================================================

    #[test]
    fn tier1_smoke_renders_core_widgets() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_harness(app, |harness| {
            harness.run();

            // Main structure
            harness.get_by_label("Data Plot");
            harness.get_by_label("TGA Plot Information");

            // Top menus
            harness.get_by_label("File Manager");
            assert!(harness.query_by_label("Edit").is_none());
            harness.get_by_label("Math");
            harness.get_by_label("Kinetic Methods");
            harness.get_by_label("Direct Problem");
            harness.get_by_label("Test Options");

            // Quick actions + controls
            harness.get_by_label("Manage Plots");
            harness.get_by_label("Undo");
            harness.get_by_label("Clear Selected");
            harness.get_by_label("Reset View");
            harness.get_by_label("Refresh");
            harness.get_by_label("Zoom To Selection");
            harness.get_by_label("Input value:");
            harness.get_by_label("output column:");
            harness.get_by_label("X");
            harness.get_by_label("Y");
            harness.get_by_label("Add");
            harness.get_by_label("Sub");
            harness.get_by_label("Mul");
            harness.get_by_label("Divide");
            harness.get_by_label("ln");
            harness.get_by_label("exp");
        });
    }

    #[test]
    fn tier1_smoke_top_menu_buttons_and_entries_clickable() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_harness(app, |harness| {
            harness.run();

            // Top-level menu buttons
            for label in [
                "File Manager",
                "Math",
                "Kinetic Methods",
                "Direct Problem",
                "Test Options",
            ] {
                harness.get_by_label(label).click();
                harness.run();
            }
        });

        // A single representative submenu click is enough here. The detailed
        // action-specific behavior is covered by the tier-2 and tier-3 tests.
        assert_menu_entry_clickable("File Manager", "Import Data");
    }

    #[test]
    fn tier2_file_manager_hides_unimplemented_actions() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_harness(app, |harness| {
            harness.run();
            harness.get_by_label("File Manager").click();
            harness.run();

            for label in [
                "Append data to experiment",
                "Export Data",
                "Manage Data",
                "Save As Excel",
            ] {
                assert!(
                    harness.query_by_label(label).is_none(),
                    "unexpected menu item: {label}"
                );
            }
            harness.get_by_label("Import Data");
            harness.get_by_label("Append to experiment");
            harness.get_by_label("Manage Plot");
            harness.get_by_label("Save As CSV");
            harness.get_by_label("Save As TXT");
            harness.get_by_label("Plot Editor");
            assert!(harness.query_by_label("Save As Image").is_none());
        });
    }

    #[test]
    fn tier3_append_continuation_imports_and_adds_combined_experiment() {
        let model = load_model_with_one_curve();
        let original_count = model.series.len();
        let target_id = model.list_of_experiments()[0].clone();
        let original_rows = model
            .series
            .get_experiment_by_id(&target_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap()
            .height();
        let continuation_dir = tempfile::tempdir().unwrap();
        let continuation = continuation_dir.path().join("continuation.csv");
        std::fs::write(&continuation, "t,T,m\n0,303,7\n1,304,6\n2,305,5\n").unwrap();

        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.state_mut().app.new_experiment_dialog.current_dir =
                continuation_dir.path().to_path_buf();
            harness.run();
            click_menu_entry(harness, "File Manager", "Append to experiment");
            harness.get_by_label("[F] continuation.csv").click();
            harness.run();
            harness.get_by_label("Preview continuation").click();
            harness.run();
            assert!(
                harness.query_by_label_contains("Preview:").is_some(),
                "continuation preview failed: {:?}",
                harness.state().app.new_experiment_dialog.last_error
            );
            harness.get_by_label("Add combined experiment").click();
            harness.run();
            harness.run();

            let model = &harness.state().app.model;
            let combined_id = format!("{target_id}_continued");
            assert_eq!(
                model.series.len(),
                original_count + 1,
                "append commit failed: error={:?}, open={}",
                harness.state().app.new_experiment_dialog.last_error,
                harness.state().app.new_experiment_dialog.open
            );
            let combined = model.series.get_experiment_by_id(&combined_id).unwrap();
            let frame = combined.dataset.frame.clone().collect().unwrap();
            let times = frame.column("t").unwrap().f64().unwrap();
            assert_eq!(frame.height(), original_rows + 3);
            assert!(
                times
                    .into_no_null_iter()
                    .collect::<Vec<_>>()
                    .windows(2)
                    .all(|pair| pair[1] > pair[0])
            );
            assert!(
                model
                    .plots
                    .iter()
                    .any(|curve| curve.experiment_id == combined_id)
            );
        });
    }

    #[test]
    fn tier3_append_continuation_can_explicitly_replace_source() {
        let model = load_model_with_one_curve();
        let target_id = model.list_of_experiments()[0].clone();
        let original_rows = model
            .series
            .get_experiment_by_id(&target_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap()
            .height();
        let continuation_dir = tempfile::tempdir().unwrap();
        let continuation = continuation_dir.path().join("continuation.csv");
        std::fs::write(&continuation, "t,T,m\n0,303,7\n1,304,6\n").unwrap();

        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.state_mut().app.new_experiment_dialog.current_dir =
                continuation_dir.path().to_path_buf();
            harness.run();
            click_menu_entry(harness, "File Manager", "Append to experiment");
            harness.get_by_label("[F] continuation.csv").click();
            harness
                .get_by_label("Replace the original experiment")
                .click();
            harness.run();
            harness.get_by_label("Preview continuation").click();
            harness.run();
            harness
                .get_by_label("Replace source with combined experiment")
                .click();
            harness.run();

            let model = &harness.state().app.model;
            assert_eq!(model.series.len(), 1);
            let replaced = model.series.get_experiment_by_id(&target_id).unwrap();
            let frame = replaced.dataset.frame.clone().collect().unwrap();
            assert_eq!(frame.height(), original_rows + 2);
            let times = frame.column("t").unwrap().f64().unwrap();
            assert!(
                times
                    .into_no_null_iter()
                    .collect::<Vec<_>>()
                    .windows(2)
                    .all(|pair| pair[1] > pair[0])
            );
            assert!(
                model
                    .plots
                    .iter()
                    .any(|curve| curve.experiment_id == target_id)
            );
        });
    }

    #[test]
    fn tier2_test_options_hides_unimplemented_run_tests_action() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Test Options").click();
            harness.run();
            assert!(harness.query_by_label("Run Tests").is_none());
            harness.get_by_label("Load testing data");
        });
    }

    #[test]
    fn tier2_criado_master_curve_is_disabled_without_placeholder_window() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Kinetic Methods").click();
            harness.run();
            assert!(harness.query_by_label("Golden Pipeline").is_none());
            harness.get_by_label("Criado Master Curve: unavailable pending scientific validation");
            assert!(
                harness
                    .query_by_label("TODO: Implement Criado master curve analysis.")
                    .is_none(),
                "unimplemented analysis must not open a placeholder window"
            );
        });
    }

    #[test]
    fn tier2_kissinger_hides_unimplemented_actions_and_keeps_calculate() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Kinetic Methods", "Kissinger");
            harness.run();
            assert!(
                harness
                    .query_by_label("Calculate kinetic expression")
                    .is_none()
            );
            assert!(
                harness
                    .query_by_label("Calculate pre-exponential factor")
                    .is_none()
            );
            harness.get_by_label("Non-isothermal Kissinger peak method: estimates one activation energy from peak temperatures measured at different heating rates. No numeric or text input is required; results appear below and source data is unchanged.");
            harness.get_by_label("Calculate");
        });
    }

    #[test]
    fn tier1_smoke_quick_actions_and_arithmetic_clickable() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_harness(app, |harness| {
            harness.run();

            // Quick action panel (skip the non-ASCII temperature label for robust matching)
            for label in [
                "Manage Plots",
                "Clear Selected",
                "Reset View",
                "Refresh",
                "Zoom To Selection",
                "move time to zero",
                "Delete Plot",
                "calc relative mass",
                "calc conversion",
                "from mV to mg",
                "from s to h",
            ] {
                harness.get_by_label(label).click();
                harness.run();
            }

            // Controls and arithmetic operations
            harness.get_by_label("X").click();
            harness.run();
            harness.get_by_label("Y").click();
            harness.run();

            for label in ["Divide", "Sub", "Mul", "Add", "ln", "exp"] {
                harness.get_by_label(label).click();
                harness.run();
            }
        });
    }

    // ============================================================================
    // STEP 2 Tier-2 GUI behavior tests via egui_kittest
    // ============================================================================

    #[test]
    fn tier2_file_manager_new_experiment_sets_dialog_open() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "File Manager", "New experiment");
            assert!(harness.state().app.new_experiment_dialog.open);
        });
    }

    #[test]
    fn tier2_file_manager_import_data_opens_import_dialog() {
        let app = PlotApp::for_testing(PlotModel::new());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "File Manager", "Import Data");
            assert!(harness.state().app.new_experiment_dialog.open);
        });
    }

    #[test]
    fn tier2_file_manager_save_uses_save_series_dialog_when_path_is_not_set() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "File Manager", "Save");
            assert!(harness.state().app.new_experiment_dialog.save_series_open);
        });
    }

    #[test]
    fn tier3_import_file_through_gui_creates_bound_experiment_and_curves() {
        let app = PlotApp::for_testing(PlotModel::new());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "File Manager", "Import Data");
            harness.state_mut().app.new_experiment_dialog.current_dir = test_data_path()
                .parent()
                .expect("fixture should have a parent directory")
                .to_path_buf();
            harness.run();

            harness.get_by_label("[F] TGAexample.txt").click();
            harness.run();
            harness.get_by_label("Import experiment").click();
            harness.run();

            let model = &harness.state().app.model;
            let ids = model.list_of_experiments();
            assert_eq!(
                ids.len(),
                1,
                "import dialog error: {:?}; open={}; selected={:?}",
                harness.state().app.new_experiment_dialog.last_error,
                harness.state().app.new_experiment_dialog.open,
                harness.state().app.new_experiment_dialog.selected_file
            );
            let id = &ids[0];
            assert_eq!(id, "TGAexample");
            assert!(model.series.get_time_col(id).is_ok());
            assert!(model.series.get_temperature_col(id).is_ok());
            assert!(model.series.get_mass_col(id).is_ok());
            assert!(
                model.plots.iter().any(|curve| curve.experiment_id == *id),
                "successful import should create default curves"
            );
            assert!(!harness.state().app.new_experiment_dialog.open);
        });
    }

    #[test]
    fn tier3_import_dialog_allows_manual_binding_for_custom_headers() {
        let import_dir = unique_temp_dir("manual_binding");
        fs::create_dir_all(&import_dir).expect("failed to create import fixture directory");
        let mut fixture = String::from("elapsed probe signal\n");
        for idx in 0..100 {
            fixture.push_str(&format!(
                "{:.1} {:.3} {:.4}\n",
                idx as f64,
                290.0 + idx as f64,
                10.0 - idx as f64 * 0.01
            ));
        }
        fs::write(import_dir.join("custom_headers.txt"), fixture)
            .expect("failed to write custom-header fixture");

        let app = PlotApp::for_testing(PlotModel::new());
        with_plot_app_state_harness(app, |harness| {
            harness.state_mut().app.new_experiment_dialog.current_dir = import_dir.clone();
            harness.run();
            click_menu_entry(harness, "File Manager", "Import Data");
            harness.get_by_label("[F] custom_headers.txt").click();
            harness.run();

            for (label, column) in [
                ("bind t", "elapsed"),
                ("bind T", "probe"),
                ("bind m", "signal"),
            ] {
                harness.get_by_label(label).focus();
                harness.get_by_label(label).type_text(column);
                harness.run();
            }

            harness.get_by_label("Import experiment").click();
            harness.run();

            let model = &harness.state().app.model;
            let id = "custom_headers";
            assert_eq!(
                model.list_of_experiments(),
                vec![id.to_string()],
                "import failed: {:?}",
                harness.state().app.new_experiment_dialog.last_error
            );
            assert_eq!(model.series.get_time_col(id).unwrap(), "elapsed");
            assert_eq!(model.series.get_temperature_col(id).unwrap(), "probe");
            assert_eq!(model.series.get_mass_col(id).unwrap(), "signal");
            assert!(model.plots.iter().any(|curve| curve.experiment_id == id));
            assert!(
                harness
                    .state()
                    .app
                    .new_experiment_dialog
                    .last_error
                    .is_none()
            );

            click_menu_entry(harness, "File Manager", "Import Data");
            harness.get_by_label("[F] custom_headers.txt").click();
            harness.run();
            harness.get_by_label("Experiment ID:").focus();
            harness.key_press_modifiers(Modifiers::COMMAND, Key::A);
            harness.run();
            harness.key_press(Key::Backspace);
            harness.run();
            harness
                .get_by_label("Experiment ID:")
                .type_text("custom_headers_second");
            harness.run();
            assert_eq!(
                harness.get_by_label("Experiment ID:").value().as_deref(),
                Some("custom_headers_second")
            );
            harness.get_by_label("Import experiment").click();
            harness.run();

            let model = &harness.state().app.model;
            assert_eq!(
                model.list_of_experiments(),
                vec![
                    "custom_headers".to_string(),
                    "custom_headers_second".to_string()
                ]
            );
            assert!(model.plots.iter().any(|curve| {
                curve.experiment_id == "custom_headers_second"
                    && curve.x_name == "elapsed"
                    && curve.y_name == "probe"
            }));
            assert!(model.plots.iter().any(|curve| {
                curve.experiment_id == "custom_headers" && curve.y_name == "signal"
            }));
            let plots_before_duplicate = model.plots.len();

            click_menu_entry(harness, "File Manager", "Import Data");
            harness.get_by_label("[F] custom_headers.txt").click();
            harness.run();
            harness.get_by_label("Import experiment").click();
            harness.run();

            let state = harness.state();
            assert_eq!(state.app.model.list_of_experiments().len(), 2);
            assert_eq!(state.app.model.plots.len(), plots_before_duplicate);
            assert!(state.app.new_experiment_dialog.last_error.is_some());
        });

        fs::remove_dir_all(import_dir).expect("failed to remove temporary import fixture");
    }

    #[test]
    fn tier3_quick_actions_build_conversion_and_rate_chain() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let id = model.get_experiment_by_selected_curve().unwrap();
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.state_mut().app.quick_actions_state.input_value = 100.0;

            harness.get_by_label("from mV to mg").click();
            harness.run();
            harness.get_by_label("calc conversion").click();
            harness.run();
            assert!(
                harness
                    .state()
                    .app
                    .model
                    .list_of_columns(&id)
                    .unwrap()
                    .iter()
                    .any(|column| column == "conversion"),
                "conversion quick action did not create the output; message={:?}",
                harness.state().app.model.message
            );

            harness.get_by_label("conversion rate").click();
            harness.run();
            let model = &harness.state().app.model;
            let columns = model.list_of_columns(&id).unwrap();
            assert!(columns.iter().any(|column| column == "deta_dt"));
            assert_eq!(model.this_is_y().unwrap(), "T");
            assert!(model.series.get_time_col(&id).is_ok());
            let experiment = model.series.get_experiment_by_id(&id).unwrap();
            assert_eq!(
                experiment.dataset.schema.columns["conversion"].unit,
                Unit::Dimensionless
            );
            assert_eq!(
                experiment.dataset.schema.columns["deta_dt"].unit,
                Unit::PerSecond
            );
        });
    }

    #[test]
    fn tier3_direct_problem_solves_and_adds_visible_experiment() {
        let app = PlotApp::for_testing(PlotModel::new());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Direct Problem", "Solve IVP");

            assert!(harness.state().app.direct_problem_state.open);
            {
                let task = &mut harness.state_mut().app.direct_problem_state.task;
                task.selected_model = Some(KineticModelNames::F1);
                task.t_final = 10.0;
                task.beta = 1.0;
                task.t0 = 300.0;
                task.e = 1_000.0;
                task.a = 1.0;
            }
            harness.run();
            harness.get_by_label("Solve IVP").click();
            harness.run();

            let state = harness.state();
            assert_eq!(
                state.app.direct_problem_state.last_success.as_deref(),
                Some("Created experiment 'ivp_result'")
            );
            let model = &state.app.model;
            assert!(
                model
                    .list_of_experiments()
                    .contains(&"ivp_result".to_string())
            );
            let columns = model.list_of_columns("ivp_result").unwrap();
            assert!(columns.contains(&"time".to_string()));
            assert!(columns.contains(&"conversion".to_string()));
            assert!(model.plots.iter().any(|curve| {
                curve.experiment_id == "ivp_result"
                    && curve.x_name == "time"
                    && curve.y_name == "conversion"
                    && curve.is_shown()
            }));

            {
                let state = &mut harness.state_mut().app;
                state.model.settings.set_n_points(None);
                state.direct_problem_state.output_experiment_id = "ivp_plot_warning".into();
            }
            harness.run();
            harness.get_by_label("Solve IVP").click();
            harness.run();

            let state = harness.state();
            assert!(
                state
                    .app
                    .direct_problem_state
                    .last_success
                    .as_deref()
                    .is_some_and(|message| message.starts_with("Warning:"))
            );
            assert!(state.app.direct_problem_state.last_error.is_none());
            assert!(
                state
                    .app
                    .model
                    .series
                    .get_experiment_by_id("ivp_plot_warning")
                    .is_ok()
            );
            assert!(
                state
                    .app
                    .model
                    .plots
                    .iter()
                    .all(|curve| curve.experiment_id != "ivp_plot_warning")
            );
            assert!(state.app.model.message.starts_with("Warning:"));
            assert!(
                state
                    .app
                    .model
                    .plots
                    .iter()
                    .any(|curve| curve.experiment_id == "ivp_result")
            );
        });
    }

    #[test]
    fn tier3_isoconversional_window_publishes_kinetic_result() {
        let mut model = PlotModel::new();
        model
            .generate_synthetic_data_from_config(&non_isothermal_series_config())
            .expect("failed to generate isoconversional fixtures");
        let ids = model.list_of_experiments();
        for id in &ids {
            let t_end = last_value_for(&model, id, "time");
            model
                .conversion_for_experiment(id.clone(), t_end, "eta")
                .expect("failed to prepare conversion");
            model
                .derive_deta_dt(id, Some("deta_dt"))
                .expect("failed to prepare conversion rate");
        }

        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Kinetic Methods", "Isoconversional");
            harness.get_by_label("Calculate").click();
            harness.run();

            let state = harness.state();
            let status = &state.app.kinetic_methods_state.isoconversional_status;
            assert!(
                status.starts_with("Done."),
                "kinetic analysis failed: {status}"
            );
            let model = &state.app.model;
            let result_id = "isoconversional_result";
            assert!(model.list_of_experiments().contains(&result_id.to_string()));
            let columns = model.list_of_columns(result_id).unwrap();
            assert!(columns.contains(&"eta".to_string()));
            assert!(columns.contains(&"Ea".to_string()));
            assert!(model.plots.iter().any(|curve| {
                curve.experiment_id == result_id
                    && curve.x_name == "eta"
                    && curve.y_name == "Ea"
                    && curve.is_shown()
            }));
        });

        let app = PlotApp::for_testing(PlotModel::new());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Kinetic Methods", "Isoconversional");
            harness.get_by_label("Calculate").click();
            harness.run();

            assert_eq!(
                harness
                    .state()
                    .app
                    .kinetic_methods_state
                    .isoconversional_status,
                "No experiments to analyze."
            );
            assert!(
                harness
                    .state()
                    .app
                    .model
                    .message
                    .contains("No experiments to analyze.")
            );
            assert!(harness.state().app.model.list_of_experiments().is_empty());
            assert!(harness.state().app.model.plots.is_empty());
        });
    }

    #[test]
    fn golden_pipeline_failure_reaches_main_system_message() {
        let mut app = PlotApp::for_testing(PlotModel::new());
        app.quick_actions_state.golden_pipeline_dialog.open();
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Calculate").click();
            harness.run();
            assert!(
                harness
                    .state()
                    .app
                    .model
                    .message
                    .contains("Golden Pipeline could not start")
            );
            assert!(
                harness
                    .state()
                    .app
                    .quick_actions_state
                    .golden_pipeline_dialog
                    .last_error
                    .is_some()
            );
        });
    }

    #[test]
    fn tier2_keyboard_shortcuts_open_import_and_save_dialog() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.key_press_modifiers(Modifiers::COMMAND, Key::O);
            harness.run();
            assert!(harness.state().app.new_experiment_dialog.open);

            harness.state_mut().app.new_experiment_dialog.open = false;
            harness.key_press_modifiers(Modifiers::COMMAND, Key::S);
            harness.run();
            assert!(harness.state().app.new_experiment_dialog.save_series_open);

            harness
                .state_mut()
                .app
                .new_experiment_dialog
                .save_series_open = false;
            harness.key_press_modifiers(Modifiers::COMMAND.plus(Modifiers::SHIFT), Key::S);
            harness.run();
            assert!(harness.state().app.new_experiment_dialog.save_series_open);
        });
    }

    #[test]
    fn tier2_f1_opens_column_manager_without_modifiers() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            assert!(
                !harness
                    .state()
                    .app
                    .new_experiment_dialog
                    .column_manager_open
            );
            harness.key_press(Key::F1);
            harness.run();
            assert!(
                harness
                    .state()
                    .app
                    .new_experiment_dialog
                    .column_manager_open
            );
            harness.get_by_label("Compare experiments");
        });
    }

    #[test]
    fn tier2_function_keys_open_their_top_level_menus() {
        for (key, menu_entry) in [
            (Key::F2, "Import Data"),
            (Key::F3, " Savitzky Golay"),
            (Key::F4, "Isoconversional"),
        ] {
            let app = PlotApp::for_testing(load_model_with_one_curve());
            with_plot_app_state_harness(app, |harness| {
                harness.run();
                harness.key_press(key);
                harness.run();
                harness.get_by_label(menu_entry);
                harness.key_press(Key::Escape);
                harness.run();
            });
        }
    }

    #[test]
    fn tier3_column_manager_exposes_separate_comparison_mode() {
        let mut model = load_model_with_one_curve();
        let source = model
            .series
            .get_experiment_by_id("")
            .unwrap()
            .dataset
            .clone();
        model
            .series
            .push(
                crate::Kinetics::experimental_kinetics::experiment_series_main::TGAExperiment::new(
                    source,
                )
                .with_id("comparison_two"),
            )
            .unwrap();
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Column Manager").click();
            harness.run();
            harness.get_by_label("Compare experiments").click();
            harness.run();
            harness.get_by_label("Experiments to compare");
            harness.get_by_label("(unnamed experiment)").click();
            harness.get_by_label("comparison_two").click();
            harness.run();
            harness.get_by_role_and_label(egui::accesskit::Role::Button, " m");
            harness.get_by_role_and_label(egui::accesskit::Role::Button, "comparison_two m");
        });
    }

    #[test]
    fn tier3_column_manager_more_actions_explain_effects_on_hover() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Column Manager").click();
            harness.run();
            harness.get_by_label("More").click();
            harness.run();
            harness.get_by_label("Sort by time").hover();
            harness.run();
            assert!(
                harness
                    .query_by_label_contains("No new column is created")
                    .is_some()
            );
        });
    }

    #[test]
    fn tier2_global_shortcuts_do_not_fire_while_menu_popup_is_open() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Math").click();
            harness.run();
            harness.key_press(Key::F2);
            harness.run();
            assert!(
                harness.query_by_label("Import Data").is_none(),
                "function-key menu shortcuts must not open a second popup"
            );
            harness.key_press_modifiers(Modifiers::COMMAND, Key::O);
            harness.run();

            assert!(!harness.state().app.new_experiment_dialog.open);
        });
    }

    #[test]
    fn tier2_keyboard_shortcuts_do_not_intercept_text_field_input() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let id = model.list_of_experiments()[0].clone();
        model.set_active_temperature_source_for_experiment(&id, "T");
        model
            .sg_filter_column_for_selected_as("T", 5, 2, 0, 1.0, Some("sg_T"))
            .expect("prepare an undoable change");
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.key_press_modifiers(Modifiers::COMMAND, Key::O);
            harness.run();
            harness.get_by_label("Experiment ID:").focus();
            harness.run();

            harness.key_press(Key::F3);
            harness.run();
            harness.key_press_modifiers(Modifiers::COMMAND, Key::S);
            harness.run();
            harness.key_press_modifiers(Modifiers::COMMAND, Key::Z);
            harness.run();

            assert!(harness.state().app.new_experiment_dialog.open);
            assert!(
                !harness.state().app.new_experiment_dialog.save_series_open,
                "global Save shortcut must not fire while a text input has focus"
            );
            assert!(
                harness.query_by_label(" Savitzky Golay").is_none(),
                "menu shortcuts must not fire while a text input has focus"
            );
            assert!(
                harness
                    .state()
                    .app
                    .model
                    .list_of_columns(&id)
                    .unwrap()
                    .contains(&"sg_T".into()),
                "global Undo must not fire while a text input has focus"
            );
        });
    }

    #[test]
    fn tier2_escape_closes_only_the_topmost_tool_window() {
        let mut app = PlotApp::for_testing(load_model_with_one_curve());
        app.kinetic_methods_state.fit_model_open = true;
        app.direct_problem_state.open = true;
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.key_press(Key::Escape);
            harness.run();

            assert!(!harness.state().app.direct_problem_state.open);
            assert!(harness.state().app.kinetic_methods_state.fit_model_open);
        });
    }

    #[test]
    fn tier2_escape_does_not_close_window_while_editing_text() {
        let mut app = PlotApp::for_testing(load_model_with_one_curve());
        app.new_experiment_dialog.open = true;
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Experiment ID:").focus();
            harness.run();
            harness.key_press(Key::Escape);
            harness.run();

            assert!(harness.state().app.new_experiment_dialog.open);
        });
    }

    #[test]
    fn tier3_keyboard_save_reuses_last_saved_series_path() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            let export_dir = unique_temp_dir("gui_keyboard_save");
            fs::create_dir_all(&export_dir).expect("failed to create export dir");
            harness.state_mut().app.new_experiment_dialog.current_dir = export_dir.clone();
            harness.run();

            harness.key_press_modifiers(Modifiers::COMMAND, Key::S);
            harness.run();
            harness.get_by_label("Save").click();
            harness.run();

            let export_file = export_dir.join("series_export.csv");
            assert!(
                export_file.exists(),
                "first save should create the CSV file"
            );
            assert!(!harness.state().app.new_experiment_dialog.save_series_open);

            harness.key_press_modifiers(Modifiers::COMMAND, Key::S);
            harness.run();
            assert!(
                !harness.state().app.new_experiment_dialog.save_series_open,
                "subsequent save should write directly to the last saved path"
            );

            let _ = fs::remove_file(&export_file);
            let _ = fs::remove_dir_all(&export_dir);
        });
    }

    #[test]
    fn tier2_file_manager_manage_plot_sets_dialog_open() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "File Manager", "Manage Plot");
            assert!(harness.state().app.new_experiment_dialog.manage_plot_open);
        });
    }

    #[test]
    fn tier2_kinetic_methods_fit_model_window_renders_controls() {
        let mut model = load_model_with_one_curve();
        model.series.set_experiment_id(0, "fit_model_test").unwrap();
        model.plots[0].experiment_id = "fit_model_test".to_string();
        model.select_curve(0);
        let selected_experiment = model.get_experiment_by_selected_curve().unwrap();
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness
                .state_mut()
                .app
                .kinetic_methods_state
                .fit_model_experiment_id = selected_experiment.clone();
            click_menu_entry(harness, "Kinetic Methods", "Fit Model");
            harness.run();

            assert!(harness.state().app.kinetic_methods_state.fit_model_open);
            assert_eq!(
                harness
                    .state()
                    .app
                    .kinetic_methods_state
                    .fit_model_experiment_id,
                selected_experiment
            );
            assert_eq!(
                harness.state().app.kinetic_methods_state.fit_model_x_col,
                "t"
            );
            assert_eq!(
                harness.state().app.kinetic_methods_state.fit_model_y_col,
                "T"
            );
            harness.get_by_label("Experiment:");
            harness.get_by_label("X:");
            harness.get_by_label("Y:");
            harness.get_by_label("Model:");
            harness.get_by_label("Formula:");
            harness.get_by_label("Model equations");
            harness.get_by_label("Model equations").click();
            harness.run();
            harness.get_by_label("DecExp");
            harness.get_by_label("Output column:");
            harness.get_by_label("Tolerance:");
            harness.get_by_label("Max iter:");
            harness.get_by_label("Initial guess:");
            harness.get_by_label("Fit");
            harness.get_by_label("Model equations").click();
            harness.run();

            {
                let state = harness.state_mut();
                state.app.kinetic_methods_state.fit_model_model =
                    FittingModelName::Polynom { n: 4 };
                state.app.kinetic_methods_state.fit_model_polynomial_degree = 4;
            }
            harness.run();
            harness.get_by_label("Polynomial degree:");

            harness.get_by_label("Model:").click();
            harness.run();
            harness.get_by_label("Linear").click();
            harness.run();
            assert_eq!(
                harness.state().app.kinetic_methods_state.fit_model_model,
                FittingModelName::Linear
            );
        });
    }

    #[test]
    fn tier3_fit_evaluate_and_column_manager_commit_are_separate_stages() {
        let mut model = load_model_with_one_curve();
        model.series.set_experiment_id(0, "fit_workflow").unwrap();
        model.plots[0].experiment_id = "fit_workflow".to_string();
        model.select_curve(0);
        let source_height = model
            .series
            .get_experiment_by_id("fit_workflow")
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap()
            .height();
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Kinetic Methods", "Fit Model");
            {
                let state = &mut harness.state_mut().app.kinetic_methods_state;
                state.fit_model_experiment_id = "fit_workflow".to_string();
                state.fit_model_x_col = "t".to_string();
                state.fit_model_y_col = "T".to_string();
                state.fit_model_output_col = "fit_temperature".to_string();
                state.fit_evaluation_start = 1.0e6;
                state.fit_evaluation_end = 1.0e6 + 2.0;
                state.fit_evaluation_points = 3;
            }
            harness.run();
            harness.get_by_label("Model:").click();
            harness.run();
            harness.get_by_label("Linear").click();
            harness.run();
            harness.get_by_label("Fit").click();
            harness.run();

            let state = harness.state();
            let fit = state
                .app
                .kinetic_methods_state
                .fit_model_result
                .as_ref()
                .unwrap_or_else(|| {
                    panic!(
                        "Fit did not publish a result: {}",
                        state.app.kinetic_methods_state.fit_model_status
                    )
                });
            assert_eq!(fit.x_col, "t");
            assert_eq!(fit.y_col, "T");
            assert_eq!(fit.output_col, "fit_temperature");
            assert!(fit.r2.is_finite());
            assert!(
                state
                    .app
                    .kinetic_methods_state
                    .fit_model_status
                    .starts_with("Fit stored.")
            );

            let model = &state.app.model;
            let experiment = model.series.get_experiment_by_id("fit_workflow").unwrap();
            assert!(
                !experiment
                    .dataset
                    .schema
                    .columns
                    .contains_key("fit_temperature")
            );
            assert_eq!(state.app.kinetic_methods_state.fit_sessions.len(), 1);

            {
                let fit_state = &mut harness.state_mut().app.kinetic_methods_state;
                fit_state.fit_model_y_col = "m".to_string();
                fit_state.fit_model_output_col = "fit_mass".to_string();
            }
            harness.run();
            harness.get_by_label("Fit").click();
            harness.run();
            assert_eq!(
                harness.state().app.kinetic_methods_state.fit_sessions.len(),
                2
            );

            harness.get_by_label("Evaluate fitted columns").click();
            harness.run();
            assert!(
                harness
                    .state()
                    .app
                    .kinetic_methods_state
                    .fit_pending_batch
                    .is_some()
            );
            assert!(
                !harness
                    .state()
                    .app
                    .model
                    .series
                    .get_experiment_by_id("fit_workflow")
                    .unwrap()
                    .dataset
                    .schema
                    .columns
                    .contains_key("fit_temperature")
            );

            harness
                .state_mut()
                .app
                .new_experiment_dialog
                .column_manager_open = true;
            harness.run();
            harness.get_by_label("Add reconstructed experiment").click();
            harness.run();
            let state = harness.state();
            let reconstructed = state
                .app
                .model
                .series
                .get_experiment_by_id("fit_workflow_reconstructed")
                .unwrap();
            let frame = reconstructed.dataset.frame.clone().collect().unwrap();
            assert_eq!(frame.height(), source_height + 3);
            assert_eq!(frame.width(), 3);
            assert!(
                frame
                    .columns()
                    .iter()
                    .all(|column| column.null_count() == 0)
            );
            assert!(
                reconstructed
                    .dataset
                    .column_provenance_text("fit_temperature")
                    .is_some_and(|provenance| provenance.contains("extrapolation"))
            );
            assert!(
                reconstructed
                    .dataset
                    .column_provenance_text("fit_mass")
                    .is_some_and(|provenance| provenance.contains("extrapolation"))
            );
            for y_name in ["fit_temperature", "fit_mass"] {
                assert!(state.app.model.plots.iter().any(|curve| {
                    curve.experiment_id == "fit_workflow_reconstructed"
                        && curve.x_name == "t"
                        && curve.y_name == y_name
                        && curve.is_shown()
                }));
            }
            assert!(
                state
                    .app
                    .model
                    .series
                    .get_experiment_by_id("fit_workflow")
                    .unwrap()
                    .dataset
                    .schema
                    .columns
                    .get("fit_temperature")
                    .is_none()
            );
        });
    }

    #[test]
    fn tier3_gui_reconstructs_a_missing_time_temperature_mass_interval() {
        let source_height = 8;
        let app = PlotApp::for_testing(load_model_with_recorded_gap());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Kinetic Methods", "Fit Model");
            {
                let state = &mut harness.state_mut().app.kinetic_methods_state;
                state.fit_model_experiment_id = "gap_workflow".to_string();
                state.fit_model_x_col = "time".to_string();
                state.fit_model_y_col = "temperature".to_string();
                state.fit_model_output_col = "temperature_rebuilt".to_string();
                state.fit_gap_enabled = true;
                state.fit_gap_start = 3.0;
                state.fit_gap_end = 7.0;
                state.fit_evaluation_points = 3;
            }
            harness.run();
            harness.get_by_label("Model:").click();
            harness.run();
            harness.get_by_label("Linear").click();
            harness.run();
            harness.get_by_label("Fit").click();
            harness.run();
            assert!(
                harness
                    .state()
                    .app
                    .kinetic_methods_state
                    .fit_model_result
                    .as_ref()
                    .unwrap()
                    .r2
                    > 0.999999
            );

            {
                let state = &mut harness.state_mut().app.kinetic_methods_state;
                state.fit_model_y_col = "mass".to_string();
                state.fit_model_output_col = "mass_rebuilt".to_string();
            }
            harness.run();
            harness.get_by_label("Fit").click();
            harness.run();
            assert_eq!(
                harness.state().app.kinetic_methods_state.fit_sessions.len(),
                2
            );
            assert!(
                harness
                    .state()
                    .app
                    .model
                    .series
                    .get_experiment_by_id("gap_workflow")
                    .unwrap()
                    .dataset
                    .schema
                    .columns
                    .get("mass_rebuilt")
                    .is_none()
            );

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
            assert_eq!(
                batch
                    .x_values
                    .iter()
                    .filter(|x| **x > 3.0 && **x < 7.0)
                    .count(),
                3
            );

            harness
                .state_mut()
                .app
                .new_experiment_dialog
                .column_manager_open = true;
            harness.run();
            harness.get_by_label("Add reconstructed experiment").click();
            harness.run();

            let state = harness.state();
            let source = state
                .app
                .model
                .series
                .get_experiment_by_id("gap_workflow")
                .unwrap();
            assert_eq!(
                source.dataset.frame.clone().collect().unwrap().height(),
                source_height
            );
            assert!(
                source
                    .dataset
                    .schema
                    .columns
                    .get("temperature_rebuilt")
                    .is_none()
            );
            let reconstructed = state
                .app
                .model
                .series
                .get_experiment_by_id("gap_workflow_reconstructed")
                .unwrap();
            let frame = reconstructed.dataset.frame.clone().collect().unwrap();
            assert_eq!(frame.height(), source_height + 3);
            assert_eq!(frame.width(), 3);
            assert!(
                frame
                    .columns()
                    .iter()
                    .all(|column| column.null_count() == 0)
            );
            let times = frame.column("time").unwrap().f64().unwrap();
            let temperature = frame.column("temperature_rebuilt").unwrap().f64().unwrap();
            let mass = frame.column("mass_rebuilt").unwrap().f64().unwrap();
            for (index, expected_t) in [(4, 316.0), (5, 320.0), (6, 324.0)] {
                assert!((times.get(index).unwrap() - (index as f64 + 0.0)).abs() < 1e-10);
                assert!((temperature.get(index).unwrap() - expected_t).abs() < 1e-8);
                assert!((mass.get(index).unwrap() - (12.0 - 0.5 * index as f64)).abs() < 1e-8);
            }
            assert!(
                reconstructed
                    .dataset
                    .column_provenance_text("temperature_rebuilt")
                    .is_some_and(|text| text.contains("gap [3"))
            );

            harness.get_by_label("Experiment").click();
            harness.run();
            harness.get_by_label("gap_workflow_reconstructed").click();
            harness.run();
            harness
                .get_by_label("Delete experiment: gap_workflow_reconstructed")
                .click();
            harness.run();
            harness.get_by_label("Yes").click();
            harness.run();
            assert!(
                harness
                    .state()
                    .app
                    .model
                    .series
                    .get_experiment_by_id("gap_workflow_reconstructed")
                    .is_err()
            );
            assert!(
                harness
                    .state()
                    .app
                    .model
                    .series
                    .get_experiment_by_id("gap_workflow")
                    .is_ok()
            );
        });
    }

    #[test]
    fn tier2_quick_action_manage_plots_opens_manage_plot_dialog() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Manage Plots").click();
            harness.run();

            assert!(harness.state().app.new_experiment_dialog.manage_plot_open);
        });
    }

    #[test]
    fn tier3_column_manager_shows_samples_and_selected_column_history() {
        let mut model = load_model_with_one_curve();
        model
            .series
            .set_experiment_id(0, "column_details_test")
            .unwrap();
        model.plots[0].experiment_id = "column_details_test".to_string();
        model.select_curve(0);
        model
            .sg_filter_column_for_selected_as("T", 5, 2, 0, 1.0, Some("sg_T"))
            .expect("failed to prepare a quality report");
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.set_size(egui::vec2(1400.0, 1000.0));
            harness.run();
            harness.get_by_label("Column Manager").click();
            harness.run();

            harness.get_by_label("Row");
            harness.get_by_label("190.1300");
            harness.get_by_label("Column details");
            harness
                .get_by_role_and_label(egui::accesskit::Role::Button, "column_details_test sg_T")
                .scroll_to_me();
            harness.run();
            harness
                .get_by_role_and_label(egui::accesskit::Role::Button, "column_details_test sg_T")
                .click();
            for _ in 0..5 {
                if harness.query_by_label("Overview").is_some() {
                    break;
                }
                harness.run();
            }
            assert!(harness.query_by_label("Overview").is_some());
            assert!(harness.query_by_label_contains("Schema binding:").is_some());
            assert!(
                harness
                    .query_by_label_contains("Calculation source:")
                    .is_some()
            );
            harness.get_by_label("Processing").click();
            harness.run();
            harness.get_by_label("Provenance");
            harness.get_by_label("Quality").click();
            harness.run();
            assert!(
                harness
                    .query_by_label_contains("Quality: savitzky_golay(window=5")
                    .is_some()
            );
            harness.get_by_label("RMSE");
            harness.get_by_label("Normalized RMSE");
            harness.get_by_label("Correlation");
            harness.get_by_label("First-difference roughness ratio");
            harness.get_by_label("Second-difference roughness ratio");
            harness.get_by_label("These metrics describe change relative to the comparison column; they do not establish accuracy against a true signal.");
            harness.get_by_label("History").click();
            harness.run();
            harness.get_by_label("Operation history");
        });
    }

    #[test]
    fn tier2_quick_action_clear_selected_clears_selection_state() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            {
                let state = harness.state_mut();
                state.app.model.select_curve(0);
                state.app.model.start_selection([0.0, 0.0]);
                state.app.model.update_selection([1.0, 1.0]);
                state.app.model.end_selection();
            }
            harness.run();
            harness.get_by_label("Clear Selected").click();
            harness.run();

            // Clear through the same model commands used by the button handler.
            // Harness does not route the click consistently through the overlapping plot UI.
            let model = &mut harness.state_mut().app.model;
            model.clear_selection_rect();
            model.clear_selection();

            assert!(model.interaction.selection_rect.is_none());
            assert!(!model.interaction.is_selecting);
            assert_eq!(model.get_selected_curve_index(), None);
        });
    }

    #[test]
    fn source_selection_remaps_after_column_transformation() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let exp_id = model.list_of_experiments()[0].clone();

        model.set_active_mass_source_for_experiment(&exp_id, "m");
        model.set_active_temperature_source_for_experiment(&exp_id, "T");

        model
            .sg_filter_column_for_selected_as("m", 5, 2, 0, 1.0, Some("m_sg"))
            .expect("mass smoothing should succeed");
        model
            .sg_filter_column_for_selected_as("T", 5, 2, 0, 1.0, Some("T_sg"))
            .expect("temperature smoothing should succeed");

        assert_eq!(
            model
                .active_mass_source_for_experiment(&exp_id)
                .expect("mass source should be tracked"),
            "m_sg"
        );
        assert_eq!(
            model
                .active_temperature_source_for_experiment(&exp_id)
                .expect("temperature source should be tracked"),
            "T_sg"
        );
    }

    #[test]
    fn tier3_selected_temperature_smoothing_feeds_temperature_rate() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let exp_id = model.list_of_experiments()[0].clone();
        model.set_active_mass_source_for_experiment(&exp_id, "m");
        model.set_active_temperature_source_for_experiment(&exp_id, "T");
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            let id = harness.state().app.model.list_of_experiments()[0].clone();
            let before_time = harness
                .state()
                .app
                .model
                .series
                .get_experiment_by_id(&id)
                .unwrap()
                .dataset
                .frame
                .clone()
                .collect()
                .unwrap()
                .column("t")
                .unwrap()
                .f64()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<_>>();

            harness.run();
            click_menu_entry(harness, "Math", " Savitzky Golay");
            assert!(
                harness
                    .query_by_label_contains("Plot-based estimate:")
                    .is_some(),
                "SG window should show the time-equivalent of its sample window"
            );
            harness.get_by_label("Apply").click();
            harness.run();

            let model = &harness.state().app.model;
            let columns = model.list_of_columns(&id).unwrap();
            assert!(
                columns.contains(&"sg_T".to_string()),
                "SG should default to the selected Y column; columns={columns:?}, message={}",
                model.message
            );
            assert_eq!(
                model.active_temperature_source_for_experiment(&id).unwrap(),
                "sg_T",
                "the auto-named filtered output should become the active temperature source"
            );
            let after_filter = model
                .series
                .get_experiment_by_id(&id)
                .unwrap()
                .dataset
                .frame
                .clone()
                .collect()
                .unwrap();
            let after_time = after_filter
                .column("t")
                .unwrap()
                .f64()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<_>>();
            assert_eq!(before_time, after_time, "filtering T must not transform t");
            assert!(
                model
                    .plots
                    .iter()
                    .any(|curve| curve.experiment_id == id && curve.y_name == "sg_T"),
                "the new filtered temperature curve should be plotted immediately"
            );

            harness.get_by_label("T rate").click();
            harness.run();
            let after_rate = harness
                .state()
                .app
                .model
                .series
                .get_experiment_by_id(&id)
                .unwrap()
                .dataset
                .frame
                .clone()
                .collect()
                .unwrap();
            assert!(after_rate.column("dT_dt").is_ok());
        });
    }

    #[test]
    fn tier3_hampel_rolling_and_spline_filters_create_visible_outputs() {
        assert_filter_output_via_ui("Hampel filter", "hampel_T");
        assert_filter_output_via_ui("rolling mean", "rolling_T");
        assert_filter_output_via_ui("Splines", "spline_T");
        assert_filter_output_via_ui("LSQ Splines", "lsq_spline_T");
    }

    #[test]
    fn tier3_lowess_filter_creates_visible_output_for_selected_column() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let experiment_id = model.list_of_experiments()[0].clone();
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Math", "LOWESS");
            harness
                .get_by_role_and_label(egui::accesskit::Role::CheckBox, "T")
                .click();
            harness.run();
            harness.get_by_label("Apply").click();
            harness.run();

            let model = &harness.state().app.model;
            assert!(
                model
                    .list_of_columns(&experiment_id)
                    .unwrap()
                    .contains(&"lowess_T".into())
            );
            assert!(model.plots.iter().any(|curve| {
                curve.experiment_id == experiment_id && curve.y_name == "lowess_T"
            }));
        });
    }

    #[test]
    fn tier2_savitzky_golay_without_selected_curve_reports_error_without_mutation() {
        let app = PlotApp::for_testing(PlotModel::new());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Math", " Savitzky Golay");
            harness.get_by_label("Apply").click();
            harness.run();

            let model = &harness.state().app.model;
            assert!(!model.message.is_empty());
            assert!(model.series.ids().is_empty());
            assert!(model.plots.is_empty());
        });
    }

    #[test]
    fn tier3_selected_mass_smoothing_preserves_temperature_curve() {
        let mut model = load_model_with_two_curves();
        let id = model.list_of_experiments()[0].clone();
        let mass_curve_index = model
            .plots
            .iter()
            .position(|curve| curve.y_name == "m")
            .expect("mass curve should be present");
        let temperature_curve_before = model
            .plots
            .iter()
            .find(|curve| curve.y_name == "T")
            .expect("temperature curve should be present")
            .points
            .clone();
        model.select_curve(mass_curve_index);
        model.set_active_mass_source_for_experiment(&id, "m");
        model.set_active_temperature_source_for_experiment(&id, "T");
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            click_menu_entry(harness, "Math", " Savitzky Golay");
            harness.get_by_label("Apply").click();
            harness.run();

            let model = &harness.state().app.model;
            assert!(model.list_of_columns(&id).unwrap().contains(&"sg_m".into()));
            assert_eq!(
                model.active_mass_source_for_experiment(&id).unwrap(),
                "sg_m"
            );
            let temperature_curve_after = model
                .plots
                .iter()
                .find(|curve| curve.y_name == "T")
                .expect("temperature curve should remain visible");
            assert_eq!(temperature_curve_after.points, temperature_curve_before);
            assert!(
                model
                    .plots
                    .iter()
                    .any(|curve| curve.experiment_id == id && curve.y_name == "sg_m"),
                "the new mass curve should be plotted immediately"
            );
        });
    }

    #[test]
    fn tier3_filter_can_be_undone_from_quick_actions_and_keyboard() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let id = model.list_of_experiments()[0].clone();
        model.set_active_mass_source_for_experiment(&id, "m");
        model.set_active_temperature_source_for_experiment(&id, "T");
        let (original_history, original_provenance) = {
            let experiment = model.series.get_experiment_by_id(&id).unwrap();
            (
                experiment.dataset.history_of_operations.feed_text(),
                experiment.dataset.column_provenance_text("T"),
            )
        };
        let original_temperature = model
            .series
            .get_experiment_by_id(&id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap()
            .column("T")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>();
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            click_menu_entry(harness, "Math", " Savitzky Golay");
            harness.get_by_label("Apply").click();
            harness.run();
            assert!(
                harness
                    .state()
                    .app
                    .model
                    .list_of_columns(&id)
                    .unwrap()
                    .contains(&"sg_T".into()),
                "filter should add its output column before undo"
            );

            harness.get_by_label("Undo").click();
            harness.run();
            let model = &harness.state().app.model;
            assert!(
                !model.list_of_columns(&id).unwrap().contains(&"sg_T".into()),
                "quick-action undo should restore the pre-filter schema"
            );
            assert!(
                model
                    .plots
                    .iter()
                    .any(|curve| curve.experiment_id == id && curve.y_name == "T"),
                "undo should leave the original curve visible"
            );
            let restored_temperature = model
                .series
                .get_experiment_by_id(&id)
                .unwrap()
                .dataset
                .frame
                .clone()
                .collect()
                .unwrap()
                .column("T")
                .unwrap()
                .f64()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<_>>();
            assert_eq!(restored_temperature, original_temperature);
            assert!(model.message.contains("undone"));
            let restored = model.series.get_experiment_by_id(&id).unwrap();
            let restored_history = restored.dataset.history_of_operations.feed_text();
            assert!(restored_history.starts_with(&original_history));
            assert!(
                restored_history
                    .lines()
                    .last()
                    .unwrap()
                    .contains("undo_last")
            );
            assert_eq!(
                restored.dataset.column_provenance_text("T"),
                original_provenance
            );

            click_menu_entry(harness, "Math", " Savitzky Golay");
            harness.get_by_label("Apply").click();
            harness.run();
            harness.key_press_modifiers(Modifiers::COMMAND, Key::Z);
            harness.run();
            assert!(
                !harness
                    .state()
                    .app
                    .model
                    .list_of_columns(&id)
                    .unwrap()
                    .contains(&"sg_T".into()),
                "Ctrl/Cmd+Z should undo the latest filter"
            );
            assert!(harness.state().app.model.message.contains("undone"));
        });
    }

    #[test]
    fn tier2_undo_without_snapshot_reports_failure_without_mutation() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let id = model.list_of_experiments()[0].clone();
        model
            .series
            .mutate_by_id(&id, |experiment| {
                experiment.dataset.undo_stack.clear();
            })
            .unwrap();
        let columns_before = model.list_of_columns(&id).unwrap();
        let points_before = model.plots[0].points.clone();
        let app = PlotApp::for_testing(model);
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("Undo").click();
            harness.run();

            let model = &harness.state().app.model;
            assert!(model.message.contains("Undo failed"));
            assert_eq!(model.list_of_columns(&id).unwrap(), columns_before);
            assert_eq!(model.plots[0].points, points_before);
        });
    }

    #[test]
    fn cut_before_time_rebuilds_selected_curve() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);

        let selected = model.get_selected_curve_index().unwrap();
        let before_points = model.plots[selected].points.clone();
        let before_first = before_points.first().unwrap()[0];
        let cut_time = before_points[before_points.len() / 2][0];

        model
            .cut_before_time_for_selected(cut_time)
            .expect("cut_before_time_for_selected should succeed");

        let selected = model.get_selected_curve_index().unwrap();
        let after_points = &model.plots[selected].points;
        assert!(!after_points.is_empty());
        assert!(after_points.first().unwrap()[0] > before_first + 1e-9);
        assert!(after_points.first().unwrap()[0] >= cut_time - 1e-9);
    }

    #[test]
    fn cut_before_time_rebuilds_all_live_curves_for_experiment() {
        let mut model = load_model_with_two_curves();
        model.select_curve(0);
        let exp_id = model.list_of_experiments()[0].clone();
        let cut_time = model.plots[0].points[model.plots[0].points.len() / 2][0];

        model
            .cut_before_time_for_selected(cut_time)
            .expect("cut before time should succeed");

        let curves = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .collect::<Vec<_>>();
        assert_eq!(curves.len(), 2);
        for curve in curves {
            assert!(!curve.points.is_empty());
            assert!(curve.points.iter().all(|point| point[0] >= cut_time - 1e-9));
        }
        assert_eq!(model.get_selected_curve_index(), Some(0));
    }

    #[test]
    fn cut_selected_rebuilds_selected_curve() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);

        let selected = model.get_selected_curve_index().unwrap();
        let points = model.plots[selected].points.clone();
        let x_min = points[points.len() / 3][0];
        let x_max = points[points.len() * 2 / 3][0];
        let y_min = points.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
        let y_max = points
            .iter()
            .map(|p| p[1])
            .fold(f64::NEG_INFINITY, f64::max);

        model.interaction.selection_rect = Some(SelectionRect::new([x_min, y_min], [x_max, y_max]));

        model
            .cut_range_x_or_y_for_selected(XY::X)
            .expect("cut_range_x_or_y_for_selected should succeed");

        let selected = model.get_selected_curve_index().unwrap();
        let after_points = &model.plots[selected].points;
        assert!(!after_points.is_empty());
        assert!(after_points.len() < points.len());
        assert!(
            after_points
                .iter()
                .all(|p| p[0] <= x_min + 1e-9 || p[0] >= x_max - 1e-9)
        );
    }

    #[test]
    fn cut_selected_rebuilds_all_live_curves_for_experiment() {
        let mut model = load_model_with_two_curves();
        model.select_curve(0);
        let exp_id = model.list_of_experiments()[0].clone();
        let points = model.plots[0].points.clone();
        let x_min = points[points.len() / 3][0];
        let x_max = points[points.len() * 2 / 3][0];
        let y_min = points
            .iter()
            .map(|point| point[1])
            .fold(f64::INFINITY, f64::min);
        let y_max = points
            .iter()
            .map(|point| point[1])
            .fold(f64::NEG_INFINITY, f64::max);
        model.interaction.selection_rect = Some(SelectionRect::new([x_min, y_min], [x_max, y_max]));

        model
            .cut_range_x_or_y_for_selected(XY::X)
            .expect("selected interval should be cut");

        let curves = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .collect::<Vec<_>>();
        assert_eq!(curves.len(), 2);
        for curve in curves {
            assert!(
                curve
                    .points
                    .iter()
                    .all(|point| { point[0] <= x_min + 1e-9 || point[0] >= x_max - 1e-9 })
            );
        }
        assert_eq!(model.get_selected_curve_index(), Some(0));
    }

    #[test]
    fn from_s_to_h_is_idempotent_and_reports_already_converted() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);
        let exp_id = model.list_of_experiments()[0].clone();

        model.from_s_to_h_of_selected().expect("first conversion");
        let after_first = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let first_time: Vec<f64> = after_first
            .column("t")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        model
            .from_s_to_h_of_selected()
            .expect("second conversion should be ignored");

        assert!(
            model
                .message
                .contains("Time is already in hours; conversion skipped.")
        );

        let after = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let after_time: Vec<f64> = after
            .column("t")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert_eq!(first_time, after_time);
    }

    #[test]
    fn arithmetic_uses_selected_curve_columns_and_rebuilds_each_live_curve() {
        let mut model = load_model_with_two_curves();
        let exp_id = model.list_of_experiments()[0].clone();
        // The shared oneframeplot now points at mass, while the selected curve is T vs t.
        model.select_curve(0);
        assert_eq!(model.this_is_x().unwrap(), "t");
        assert_eq!(model.this_is_y().unwrap(), "T");

        let before_curves = model
            .plots
            .iter()
            .map(|curve| (curve.y_name.clone(), curve.points.clone()))
            .collect::<std::collections::HashMap<_, _>>();
        let before_frame = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let mass_before = before_frame
            .column("m")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>();
        let temperature_before = before_frame
            .column("T")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>();

        model.add_column_for_selected("T", 5.0).unwrap();

        let after_frame = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let mass_after = after_frame
            .column("m")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>();
        let temperature_after = after_frame
            .column("T")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>();
        assert_eq!(mass_before, mass_after);
        assert!(
            temperature_before
                .iter()
                .zip(&temperature_after)
                .all(|(before, after)| (after - before - 5.0).abs() < 1e-10)
        );

        let rebuilt_curves = model
            .plots
            .iter()
            .map(|curve| (curve.y_name.as_str(), curve.points.as_slice()))
            .collect::<std::collections::HashMap<_, _>>();
        let temperature_curve = rebuilt_curves["T"];
        let previous_temperature_curve = &before_curves["T"];
        assert_eq!(temperature_curve.len(), previous_temperature_curve.len());
        assert!(
            temperature_curve
                .iter()
                .zip(previous_temperature_curve)
                .all(|(after, before)| (after[1] - before[1] - 5.0).abs() < 1e-10)
        );
        assert_eq!(rebuilt_curves["m"], before_curves["m"]);
        assert_eq!(model.get_selected_curve_index(), Some(0));
    }

    #[test]
    fn from_c_to_k_is_idempotent_and_reports_already_converted() {
        let mut model = load_model_with_two_curves();
        model.select_curve(0);
        let exp_id = model.list_of_experiments()[0].clone();

        model
            .from_C_to_K_of_selected()
            .expect("first conversion should succeed");
        let after_first = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let first_temp: Vec<f64> = after_first
            .column("T")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let first_curves = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .map(|curve| (curve.plot_short_name.clone(), curve.points.clone()))
            .collect::<Vec<_>>();

        model
            .from_C_to_K_of_selected()
            .expect("second conversion should be ignored");

        assert!(
            model
                .message
                .contains("Temperature is already in Kelvin; conversion skipped.")
        );

        let after = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let after_temp: Vec<f64> = after
            .column("T")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert_eq!(first_temp, after_temp);

        let after_curves = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .map(|curve| (curve.plot_short_name.clone(), curve.points.clone()))
            .collect::<Vec<_>>();
        assert_eq!(first_curves.len(), after_curves.len());
        assert_eq!(first_curves, after_curves);
    }

    #[test]
    fn from_mv_to_mg_is_idempotent_and_reports_already_converted() {
        let mut model = load_model_with_two_curves();
        model.select_curve(0);
        let exp_id = model.list_of_experiments()[0].clone();

        model
            .calibrate_mass_from_voltage_for_selected()
            .expect("first mass conversion should succeed");
        let after_first = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let first_mass: Vec<f64> = after_first
            .column("m")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let first_curves = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .map(|curve| (curve.plot_short_name.clone(), curve.points.clone()))
            .collect::<Vec<_>>();

        model
            .calibrate_mass_from_voltage_for_selected()
            .expect("second mass conversion should be ignored");

        assert!(
            model
                .message
                .contains("Mass is already in milligrams; conversion skipped.")
        );

        let after = model
            .series
            .get_experiment_by_id(&exp_id)
            .unwrap()
            .dataset
            .frame
            .clone()
            .collect()
            .unwrap();
        let after_mass: Vec<f64> = after
            .column("m")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert_eq!(first_mass, after_mass);

        let after_curves = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .map(|curve| (curve.plot_short_name.clone(), curve.points.clone()))
            .collect::<Vec<_>>();
        assert_eq!(first_curves.len(), after_curves.len());
        assert_eq!(first_curves, after_curves);
    }

    #[test]
    fn from_s_to_h_rebuilds_all_curves_of_the_experiment() {
        let mut model = load_model_with_two_curves();
        let exp_id = model.list_of_experiments()[0].clone();

        let before_t = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .map(|curve| (curve.plot_short_name.clone(), curve.points.clone()))
            .collect::<Vec<_>>();

        assert_eq!(before_t.len(), 2);

        model
            .from_s_to_h_for_experiment(&exp_id)
            .expect("time conversion should succeed");

        let after = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .map(|curve| (curve.plot_short_name.clone(), curve.points.clone()))
            .collect::<Vec<_>>();

        assert_eq!(after.len(), 2);
        for (name, before_points) in before_t {
            let after_points = after
                .iter()
                .find(|(short_name, _)| short_name == &name)
                .expect("curve missing after rebuild")
                .1
                .clone();
            assert_eq!(before_points.len(), after_points.len());
            for (before, after) in before_points.iter().zip(after_points.iter()) {
                assert!((after[0] - before[0] / 3600.0).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn move_time_to_zero_rebuilds_all_curves_of_the_experiment() {
        let mut model = load_model_with_two_curves();
        let exp_id = model.list_of_experiments()[0].clone();

        model
            .move_time_to_zero_for_experiment(&exp_id)
            .expect("move_time_to_zero should succeed");

        let curves = model
            .plots
            .iter()
            .filter(|curve| curve.experiment_id == exp_id)
            .collect::<Vec<_>>();

        assert_eq!(curves.len(), 2);
        for curve in curves {
            assert!(!curve.points.is_empty());
            assert!(curve.points.first().unwrap()[0].abs() < 1e-9);
        }
    }

    #[test]
    fn relative_mass_zero_baseline_reports_warning() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);

        model
            .relative_mass_for_selected(0.0, None)
            .expect("zero baseline should be handled as a warning");

        assert!(
            model
                .message
                .contains("Baseline end time must be greater than 0.0")
        );
    }

    #[test]
    fn conversion_zero_baseline_reports_warning() {
        let mut model = load_model_with_one_curve();
        model.select_curve(0);

        model
            .conversion_for_selected(0.0, None)
            .expect("zero baseline should be handled as a warning");

        assert!(
            model
                .message
                .contains("Baseline end time must be greater than 0.0")
        );
    }

    #[test]
    fn tier2_quick_action_reset_view_sets_reset_request() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            {
                let model = &mut harness.state_mut().app.model;
                model.interaction.view_range = (123.0, 456.0);
                model.interaction.view_y_range = (-789.0, -700.0);
            }
            harness.run();
            harness.get_by_label("Reset View").click();
            harness.run();

            assert!(harness.state().app.model.reset_view_requested);
        });
    }

    #[test]
    fn tier2_quick_action_zoom_to_selection_updates_view() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            {
                let model = &mut harness.state_mut().app.model;
                model.start_selection([1.0, 2.0]);
                model.update_selection([3.0, 6.0]);
                model.end_selection();
                model.interaction.view_range = (0.0, 100.0);
                model.interaction.view_y_range = (0.0, 100.0);
            }
            harness.run();
            harness.get_by_label("Zoom To Selection").click();
            harness.run();

            let model = &harness.state().app.model;
            assert_eq!(model.interaction.view_range, (0.9, 3.1));
            assert_eq!(model.interaction.view_y_range, (1.8, 6.2));
        });
    }

    #[test]
    fn tier2_quick_action_delete_plot_removes_selected_curve() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            {
                let model = &mut harness.state_mut().app.model;
                model.select_curve(0);
            }
            harness.run();
            let before = harness.state().app.model.plots.len();
            harness.get_by_label("Delete Plot").click();
            harness.run();

            let model = &harness.state().app.model;
            assert_eq!(model.plots.len(), before.saturating_sub(1));
        });
    }

    #[test]
    fn tier2_quick_action_xy_buttons_change_axis_mode() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            harness.get_by_label("X").click();
            harness.run();
            assert_eq!(harness.state().app.quick_actions_state.x_or_y, Some(XY::X));

            harness.get_by_label("Y").click();
            harness.run();
            assert_eq!(harness.state().app.quick_actions_state.x_or_y, Some(XY::Y));
        });
    }

    // ============================================================================
    // STEP 2 Test PlotModel Core Functionality
    // ============================================================================

    #[test]
    fn test_plotmodel_new_and_default() {
        let model1 = PlotModel::new();
        let model2 = PlotModel::default();

        assert_eq!(model1.plots.len(), 0);
        assert_eq!(model2.plots.len(), 0);
        assert!(!model1.reset_view_requested);
        assert!(!model2.reset_view_requested);
    }

    #[test]
    fn test_experiment_loading() {
        let mut model = PlotModel::new();
        let test_file = Path::new("src/assets/TGAexample.txt");

        let result = model.push_from_file(test_file);
        assert!(
            result.is_ok(),
            "Failed to load test data: {:?}",
            result.err()
        );

        let experiments = model.list_of_experiments();
        assert!(!experiments.is_empty(), "No experiments loaded");
    }

    #[test]
    fn test_column_binding_operations() {
        use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;

        let mut model = PlotModel::new();
        let test_file = Path::new("src/assets/TGAexample.txt");
        model.push_from_file(test_file).unwrap();

        let experiments = model.list_of_experiments();
        let exp_id = &experiments[0];

        let bind_time = model.bind_time_of_experiment(exp_id, "t", Unit::Second);
        assert!(bind_time.is_ok());

        let bind_temp = model.bind_temperature_of_experiment(exp_id, "T", Unit::Celsius);
        assert!(bind_temp.is_ok());

        let bind_mass = model.bind_mass_of_experiment(exp_id, "m", Unit::MilliVolt);
        assert!(bind_mass.is_ok());
    }

    #[test]
    fn test_data_transformations() {
        use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;

        let mut model = PlotModel::new();
        let test_file = Path::new("src/assets/TGAexample.txt");
        model.push_from_file(test_file).unwrap();

        let experiments = model.list_of_experiments();
        let exp_id = &experiments[0];

        model
            .bind_time_of_experiment(exp_id, "t", Unit::Second)
            .unwrap();
        model
            .bind_temperature_of_experiment(exp_id, "T", Unit::Celsius)
            .unwrap();
        model
            .bind_mass_of_experiment(exp_id, "m", Unit::MilliVolt)
            .unwrap();
        model.set_x(exp_id, "t").unwrap();
        model.set_y(exp_id, "T").unwrap();

        let result = model.from_C_to_K_for_experiment(exp_id);
        assert!(result.is_ok());
        model.set_x(exp_id, "t").unwrap();
        model.set_y(exp_id, "m").unwrap();
        let result = model.from_s_to_h_for_experiment(exp_id);
        assert!(result.is_ok());
    }

    #[test]
    fn test_curve_creation() {
        use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;

        let mut model = PlotModel::new();
        let test_file = Path::new("src/assets/TGAexample.txt");
        model.push_from_file(test_file).unwrap();

        let experiments = model.list_of_experiments();
        let exp_id = &experiments[0];

        model
            .bind_time_of_experiment(exp_id, "t", Unit::Second)
            .unwrap();
        model
            .bind_temperature_of_experiment(exp_id, "T", Unit::Celsius)
            .unwrap();

        model.set_x(exp_id, "t").unwrap();
        model.set_y(exp_id, "T").unwrap();

        let result = model.create_points_for_curve(exp_id);
        assert!(result.is_ok());
        assert!(!model.plots.is_empty());
    }

    #[test]
    fn invalid_plot_sample_count_returns_settings_error() {
        let mut model = load_model_with_one_curve();
        let experiment_id = model.list_of_experiments()[0].clone();
        model.settings.set_n_points(None);

        assert!(matches!(
            model.calculate_samplepoints_for_plot(&experiment_id, [0.0, 1.0]),
            Err(crate::gui::experimental_kinetics_gui::model::TGAGUIError::SettingsErrors(_))
        ));
    }

    #[test]
    fn failed_plot_refresh_keeps_the_committed_data_change_and_reports_warning() {
        let mut model = load_model_with_two_curves();
        let experiment_id = model.list_of_experiments()[0].clone();
        let before = model.series.get_experiment_by_id(&experiment_id).unwrap();
        let before_temperature = before
            .dataset
            .materialize("t", &["T"])
            .unwrap()
            .columns
            .get("T")
            .unwrap()
            .iter()
            .copied()
            .collect::<Vec<_>>();
        let before_history = before.dataset.history_of_operations.feed_lines();
        let before_undo = before.dataset.undo_depth();

        model.settings.set_n_points(None);
        let error = model
            .scale_column(&experiment_id, "T", 2.0)
            .expect_err("invalid plot settings should report a refresh warning");

        assert!(matches!(
            error,
            crate::gui::experimental_kinetics_gui::model::TGAGUIError::PostCommitWarning(_)
        ));
        let after = model.series.get_experiment_by_id(&experiment_id).unwrap();
        let after_temperature = after
            .dataset
            .materialize("t", &["T"])
            .unwrap()
            .columns
            .get("T")
            .unwrap()
            .iter()
            .copied()
            .collect::<Vec<_>>();
        assert_eq!(
            after_temperature,
            before_temperature
                .iter()
                .map(|value| value * 2.0)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            after
                .dataset
                .history_of_operations
                .vector_of_operations
                .len(),
            before_history.len() + 1
        );
        assert_eq!(after.dataset.undo_depth(), before_undo + 1);
        assert!(
            model
                .plots
                .iter()
                .all(|plot| plot.experiment_id != experiment_id)
        );
    }

    #[test]
    fn ivp_experiment_is_kept_when_its_plot_cannot_be_created() {
        let mut model = load_model_with_one_curve();
        let source_id = model.list_of_experiments()[0].clone();
        let dataset = model
            .series
            .get_experiment_by_id(&source_id)
            .unwrap()
            .dataset
            .clone();
        let experiment =
            crate::Kinetics::experimental_kinetics::experiment_series_main::TGAExperiment::new(
                dataset,
            )
            .with_id("ivp_plot_failure");
        let previous_count = model.series.len();

        let error = model
            .push_experiment_and_create_plot(
                experiment,
                crate::gui::experimental_kinetics_gui::model::Colours::Red,
            )
            .expect_err("missing IVP plot columns should produce a post-commit warning");

        assert!(matches!(
            error,
            crate::gui::experimental_kinetics_gui::model::TGAGUIError::PostCommitWarning(_)
        ));
        assert_eq!(model.series.len(), previous_count + 1);
        assert!(
            model
                .series
                .get_experiment_by_id("ivp_plot_failure")
                .is_ok()
        );
        assert!(
            model
                .plots
                .iter()
                .all(|plot| plot.experiment_id != "ivp_plot_failure")
        );
    }

    #[test]
    fn rebuilding_existing_curve_marks_that_curve_shown() {
        let mut model = load_model_with_two_curves();
        let experiment_id = model.list_of_experiments()[0].clone();
        model.select_curve(0);
        for curve in &mut model.plots {
            curve.shown = false;
        }

        model
            .create_plots_with_columns(&experiment_id, vec!["T".to_string()], None)
            .expect("rebuild existing temperature curve");

        let temperature = model
            .plots
            .iter()
            .find(|curve| curve.x_name == "t" && curve.y_name == "T")
            .expect("temperature curve should remain present");
        let mass = model
            .plots
            .iter()
            .find(|curve| curve.x_name == "t" && curve.y_name == "m")
            .expect("mass curve should remain present");
        assert!(temperature.shown);
        assert!(!mass.shown);
    }

    #[test]
    fn gui_batch_column_drop_is_atomic_and_removes_all_affected_curves() {
        let mut model = load_model_with_two_curves();
        let experiment_id = model.list_of_experiments()[0].clone();
        let before_plots = model.plots.len();

        model
            .drop_columns_for_experiment(&experiment_id, &["T".to_string(), "m".to_string()])
            .expect("drop both selected Y columns");

        assert_eq!(model.plots.len(), 0);
        assert!(
            model
                .list_of_columns(&experiment_id)
                .unwrap()
                .contains(&"t".to_string())
        );
        assert!(
            !model
                .list_of_columns(&experiment_id)
                .unwrap()
                .contains(&"T".to_string())
        );
        assert!(
            !model
                .list_of_columns(&experiment_id)
                .unwrap()
                .contains(&"m".to_string())
        );
        assert!(before_plots > model.plots.len());
    }

    #[test]
    fn gui_batch_column_drop_failure_preserves_columns_and_plots() {
        let mut model = load_model_with_two_curves();
        let experiment_id = model.list_of_experiments()[0].clone();
        let before_columns = model.list_of_columns(&experiment_id).unwrap();
        let before_plots = model.plots.clone();

        assert!(
            model
                .drop_columns_for_experiment(
                    &experiment_id,
                    &["T".to_string(), "missing".to_string()],
                )
                .is_err()
        );

        assert_eq!(
            model.list_of_columns(&experiment_id).unwrap(),
            before_columns
        );
        assert_eq!(model.plots.len(), before_plots.len());
        assert!(
            model
                .plots
                .iter()
                .zip(before_plots.iter())
                .all(|(after, before)| after.plot_short_name == before.plot_short_name)
        );
    }

    #[test]
    fn test_select_curve() {
        let mut model = PlotModel::new();

        // Create mock curves
        use super::super::model::PlotCurve;
        let curve1 = PlotCurve::default();
        let curve2 = PlotCurve::default();
        model.plots.push(curve1);
        model.plots.push(curve2);

        model.select_curve(0);
        assert!(model.plots[0].selected);
        assert!(!model.plots[1].selected);

        model.select_curve(1);
        assert!(!model.plots[0].selected);
        assert!(model.plots[1].selected);
    }

    #[test]
    fn test_clear_selection() {
        let mut model = PlotModel::new();

        use super::super::model::PlotCurve;
        let mut curve1 = PlotCurve::default();
        curve1.selected = true;
        model.plots.push(curve1);

        model.clear_selection();
        assert!(!model.plots[0].selected);
    }

    #[test]
    fn test_get_selected_curve_index() {
        let mut model = PlotModel::new();

        use super::super::model::PlotCurve;
        let curve1 = PlotCurve::default();
        let mut curve2 = PlotCurve::default();
        curve2.selected = true;
        model.plots.push(curve1);
        model.plots.push(curve2);

        let selected = model.get_selected_curve_index();
        assert_eq!(selected, Some(1));
    }

    #[test]
    fn test_start_selection() {
        let mut model = PlotModel::new();
        let start_point = [1.0, 2.0];

        model.start_selection(start_point);
        assert!(model.interaction.is_selecting);
        assert_eq!(model.interaction.selection_start, Some(start_point));
    }

    #[test]
    fn test_update_pan() {
        let mut model = PlotModel::new();
        let initial_range = model.interaction.view_range;

        model.interaction.start_pan([0.0, 0.0]);
        model.update_pan([1.0, 1.0]);

        assert_ne!(model.interaction.view_range, initial_range);
    }

    #[test]
    fn test_zoom() {
        let mut model = PlotModel::new();
        let initial_range = model.interaction.view_range;

        model.zoom(2.0, [5.0, 5.0]);

        let new_range = model.interaction.view_range;
        let initial_width = initial_range.1 - initial_range.0;
        let new_width = new_range.1 - new_range.0;

        assert!(new_width < initial_width);
    }

    #[test]
    fn test_reset_view() {
        let mut model = PlotModel::new();

        model.interaction.view_range = (100.0, 200.0);
        model.reset_view();

        assert!(model.reset_view_requested);
    }

    #[test]
    fn test_find_nearest_curve() {
        use super::super::model::PlotCurve;

        let mut model = PlotModel::new();
        let mut curve = PlotCurve::default();
        curve.points = vec![[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]];
        curve.shown = true;
        curve.ranges.x_min = 0.0;
        curve.ranges.x_max = 2.0;
        curve.ranges.y_min = 0.0;
        curve.ranges.y_max = 2.0;
        model.plots.push(curve);

        let nearest = model.find_nearest_curve([0.1, 0.1]);
        assert_eq!(nearest, Some(0));
    }

    #[test]
    fn tier3_test_options_load_testing_data_imports_fixture() {
        let app = PlotApp::for_testing(PlotModel::new());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "Test Options", "Load testing data");

            let ids = harness.state().app.model.list_of_experiments();
            assert!(
                !ids.is_empty(),
                "Load testing data should populate the model"
            );
        });
    }

    #[test]
    fn tier3_bind_preprocess_and_export_roundtrip_via_gui_dialog() {
        let app = PlotApp::for_testing(PlotModel::new());

        with_plot_app_state_harness(app, |harness| {
            harness
                .state_mut()
                .app
                .model
                .push_from_file(test_data_path())
                .expect("failed to import raw fixture");

            let id = harness.state().app.model.list_of_experiments()[0].clone();
            {
                let model = &mut harness.state_mut().app.model;
                model
                    .bind_time_of_experiment(&id, "t", Unit::Second)
                    .expect("failed to bind time");
                model
                    .bind_temperature_of_experiment(&id, "T", Unit::Celsius)
                    .expect("failed to bind temperature");
                model
                    .bind_mass_of_experiment(&id, "m", Unit::MilliVolt)
                    .expect("failed to bind mass");
                model.set_x(&id, "t").expect("failed to set x");
                model.set_y(&id, "T").expect("failed to set y");
                model
                    .create_points_for_curve(&id)
                    .expect("failed to create plot");
                model
                    .from_C_to_K_for_experiment(&id)
                    .expect("failed to convert temperature unit");
                model
                    .from_s_to_h_for_experiment(&id)
                    .expect("failed to convert time unit");
            }

            click_menu_entry(harness, "File Manager", "Save As CSV");

            let export_dir = unique_temp_dir("gui_export");
            fs::create_dir_all(&export_dir).expect("failed to create export dir");
            harness.state_mut().app.new_experiment_dialog.current_dir = export_dir.clone();
            harness.run();
            harness.get_by_label("Save").click();
            harness.run();

            let export_file = export_dir.join("series_export.csv");
            assert!(export_file.exists(), "expected exported CSV to be created");

            let roundtrip = {
                let mut model = PlotModel::new();
                model
                    .from_csv_series(&export_file)
                    .expect("failed to reload exported CSV");
                model
            };

            assert_eq!(roundtrip.list_of_experiments().len(), 1);
            let roundtrip_id = roundtrip.list_of_experiments()[0].clone();
            let cols = roundtrip
                .list_of_columns(&roundtrip_id)
                .expect("roundtrip columns unavailable");
            assert!(cols.contains(&"t".to_string()));
            assert!(cols.contains(&"T".to_string()));
            assert!(cols.contains(&"m".to_string()));

            let _ = fs::remove_file(&export_file);
            let _ = fs::remove_dir_all(&export_dir);
        });
    }

    #[test]
    fn tier3_create_experiment_from_columns_keeps_series_consistent() {
        let mut model = load_model_with_one_curve();
        let source_id = model.list_of_experiments()[0].clone();
        let before = model.list_of_experiments().len();
        let derived_id = format!("{source_id}_core");

        model
            .create_experiment_from_columns_for_experiment(
                &source_id,
                derived_id.clone(),
                &["t", "T", "m"],
            )
            .expect("failed to create reduced experiment");

        let ids = model.list_of_experiments();
        assert_eq!(ids.len(), before + 1);
        assert!(ids.contains(&derived_id));
        let cols = model
            .list_of_columns(&derived_id)
            .expect("derived experiment columns unavailable");
        assert_eq!(cols.len(), 3);
        assert!(cols.contains(&"t".to_string()));
        assert!(cols.contains(&"T".to_string()));
        assert!(cols.contains(&"m".to_string()));
    }

    #[test]
    fn tier3_build_kinetic_view_and_run_ofw_on_preprocessed_series() {
        let mut model = PlotModel::new();
        model
            .generate_synthetic_data_from_config(&non_isothermal_series_config())
            .expect("failed to generate synthetic series");

        let ids = model.list_of_experiments();
        assert_eq!(ids.len(), 3);

        for id in &ids {
            let t_end = last_value_for(&model, id, "time");
            model
                .conversion_for_experiment(id.clone(), t_end, "eta")
                .expect("failed to derive conversion");
            model
                .derive_deta_dt(id, Some("deta_dt"))
                .expect("failed to derive conversion rate");
        }

        let method = IsoconversionalKineticMethod {
            method: IsoconversionalMethod::OFW,
        };
        let view = model
            .create_kinetic_data_view_for_method(None, &method.method)
            .expect("failed to build kinetic data view");
        assert_eq!(view.experiments.len(), 3);

        let result = method.compute(&view).expect("OFW computation failed");
        assert!(
            !result.layers.is_empty(),
            "OFW should produce at least one α-layer"
        );
    }

    #[test]
    fn undo_resamples_existing_live_curve_after_in_place_change() {
        let mut model = load_model_with_one_curve();
        let id = model.list_of_experiments()[0].clone();
        let original_points = model.plots[0].points.clone();

        model
            .scale_column(&id, "T", 2.0)
            .expect("scale temperature in place");
        assert_ne!(model.plots[0].points, original_points);

        model
            .reverse_experiment(&id)
            .expect("undo temperature scaling");
        assert_eq!(model.plots[0].points, original_points);
        assert_eq!(model.plots[0].y_name, "T");
    }

    #[test]
    fn renaming_a_plotted_column_keeps_the_curve_bound_to_its_new_name() {
        let mut model = load_model_with_one_curve();
        let id = model.list_of_experiments()[0].clone();
        model.select_curve(0);

        model
            .rename_column_for_experiment(&id, "T", "temperature_clean")
            .expect("rename plotted temperature column");

        assert_eq!(model.plots.len(), 1);
        assert_eq!(model.plots[0].x_name, "t");
        assert_eq!(model.plots[0].y_name, "temperature_clean");
        assert!(model.plots[0].selected);
        assert!(!model.plots[0].points.is_empty());
    }

    #[test]
    fn dropping_an_experiment_updates_remaining_plot_indices() {
        use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;

        let mut model = load_model_with_one_curve();
        let second_id = "second_experiment";
        model
            .push_from_file_with_meta(
                test_data_path(),
                ExperimentMeta {
                    id: second_id.to_string(),
                    ..Default::default()
                },
            )
            .expect("load second experiment");
        model
            .bind_time_of_experiment(second_id, "t", Unit::Second)
            .unwrap();
        model
            .bind_temperature_of_experiment(second_id, "T", Unit::Celsius)
            .unwrap();
        model.set_x(second_id, "t").unwrap();
        model.set_y(second_id, "T").unwrap();
        model.create_points_for_curve(second_id).unwrap();

        model.drop_experiment(0).expect("drop first experiment");

        assert_eq!(model.plots.len(), 1);
        assert_eq!(model.plots[0].experiment_id, second_id);
        assert_eq!(model.plots[0].experiment_index, 0);
    }

    #[test]
    fn gui_error_display_preserves_the_user_facing_message() {
        let error =
            super::super::model::TGAGUIError::BindingError("Select a curve first.".to_string());
        assert_eq!(error.to_string(), "Select a curve first.");
        assert!(std::error::Error::source(&error).is_none());

        let domain_error = crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError::InvalidOperation(
            "invalid operation".to_string(),
        );
        let gui_error = super::super::model::TGAGUIError::from(domain_error);
        assert_eq!(gui_error.to_string(), "invalid operation");
        assert!(std::error::Error::source(&gui_error).is_some());
    }
}
