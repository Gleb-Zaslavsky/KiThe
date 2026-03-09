use crate::gui::experimental_kinetics_gui::model::PlotModel;

use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
use crate::gui::experimental_kinetics_gui::controller_table::ColumnManagerState;
use crate::gui::experimental_kinetics_gui::test_options::TestOptions;

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

#[cfg(test)]
mod tests {

    use super::*;
    use crate::Kinetics::experimental_kinetics::exp_engine_api::XY;
    use crate::gui::experimental_kinetics_gui::experimental_kinetics_gui_main::PlotApp;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
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
        let mut harness = Harness::new_state(
            |ctx, state: &mut PlotAppTestState| {
                state.app.show(ctx, &mut state.open);
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
            harness.get_by_label("Math");
            harness.get_by_label("Kinetic Methods");
            harness.get_by_label("Direct Problem");
            harness.get_by_label("Test Options");

            // Quick actions + controls
            harness.get_by_label("Manage Plots");
            harness.get_by_label("Clear Selected");
            harness.get_by_label("Reset View");
            harness.get_by_label("Refresh");
            harness.get_by_label("Zoom To Selection");
            harness.get_by_label("Input value:");
            harness.get_by_label("new column:");
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

            // File Manager
            click_menu_entry(harness, "File Manager", "Import Data");
            click_menu_entry(harness, "File Manager", "Export Data");
            click_menu_entry(harness, "File Manager", "Manage Data");
            click_menu_entry(harness, "File Manager", "Save As Image");

            // Math
            click_menu_entry(harness, "Math", "Filter Data");
            click_menu_entry(harness, "Math", "Differentiate");
            click_menu_entry(harness, "Math", "Smooth");
            click_menu_entry(harness, "Math", "Average Data");
            click_menu_entry(harness, "Math", "Splines");

            // Kinetic Methods
            click_menu_entry(harness, "Kinetic Methods", "Estimate Rates");
            click_menu_entry(harness, "Kinetic Methods", "Fit Mechanism");

            // Direct Problem
            click_menu_entry(harness, "Direct Problem", "Solve IVP");
            click_menu_entry(harness, "Direct Problem", "Solve BVP");

            // Test Options
            click_menu_entry(harness, "Test Options", "Run Tests");
            click_menu_entry(harness, "Test Options", "Show Logs");
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
    fn tier2_file_manager_manage_plot_sets_dialog_open() {
        let app = PlotApp::for_testing(load_model_with_one_curve());
        with_plot_app_state_harness(app, |harness| {
            harness.run();
            click_menu_entry(harness, "File Manager", "Manage Plot");
            assert!(harness.state().app.new_experiment_dialog.manage_plot_open);
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

            let model = &harness.state().app.model;
            assert_eq!(model.get_selected_curve_index(), None);
            assert!(model.interaction.selection_rect.is_none());
        });
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
}
