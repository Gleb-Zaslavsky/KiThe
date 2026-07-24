//! Story-style regression tests for the thermochemistry GUI.
//!
//! Hypotheses:
//! - Transport-only CEA records stay out of the thermochemistry library picker.
//! - The thermochemistry screen still renders the expected controls for NASA/NIST flows.
//! - A canonical NASA entry can still produce nonzero Cp, dH, and dS values.
//!
//! Expected result:
//! the thermochemistry screen should keep a clean backend boundary and remain calculable.

#[cfg(test)]
mod tests {
    use super::super::thermochemistry_gui::ThermochemistryApp;
    use crate::Thermodynamics::DBhandlers::thermo_api::EnergyUnit;
    use crate::Thermodynamics::thermo_lib_api::ThermoLibraryError;
    use eframe::egui;
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::{NodeT, Queryable};
    use std::cell::RefCell;
    use std::rc::Rc;

    fn parse_report_value(report: &str, prefix: &str) -> Option<f64> {
        report
            .lines()
            .find(|line| line.starts_with(prefix))
            .and_then(|line| line.split_once("= "))
            .and_then(|(_, rest)| rest.split_whitespace().next())
            .and_then(|value| value.parse::<f64>().ok())
    }

    fn assert_thermochemistry_matches_canonical(library: &str, substance: &str, temperature: f64) {
        let mut app = ThermochemistryApp::new();
        app.selected_library = library.to_string();
        app.update_substances_for_library();
        app.selected_substance = substance.to_string();
        app.temperature = temperature.to_string();
        app.set_energy_unit(EnergyUnit::J);

        let mut canonical = app
            .build_calculation_subs_data(substance)
            .expect("thermochemistry GUI should build canonical SubsData for comparison");
        canonical
            .extract_thermal_coeffs(substance, temperature)
            .expect("canonical thermochemistry calculation should extract coefficients");
        let (expected_cp, expected_dh, expected_ds) = canonical
            .calculate_thermo_properties(substance, temperature)
            .expect("canonical thermochemistry calculation should succeed");

        app.calculate_properties();

        let report_cp = parse_report_value(&app.search_results, "Cp")
            .expect("thermochemistry report should contain Cp");
        let report_dh = parse_report_value(&app.search_results, "dH")
            .expect("thermochemistry report should contain dH");
        let report_ds = parse_report_value(&app.search_results, "dS")
            .expect("thermochemistry report should contain dS");
        let expected_cp_display = format!("{:.3}", expected_cp)
            .parse::<f64>()
            .expect("rounded Cp should parse");
        let expected_dh_display = format!("{:.3}", expected_dh / 1000.0)
            .parse::<f64>()
            .expect("rounded dH should parse");
        let expected_ds_display = format!("{:.3}", expected_ds)
            .parse::<f64>()
            .expect("rounded dS should parse");

        assert_eq!(app.calculated_cp, Some(expected_cp));
        assert_eq!(app.calculated_dh, Some(expected_dh));
        assert_eq!(app.calculated_ds, Some(expected_ds));
        assert_eq!(report_cp, expected_cp_display);
        assert_eq!(report_dh, expected_dh_display);
        assert_eq!(report_ds, expected_ds_display);
        assert!(app.search_results.contains("J/(mol*K)"));
        assert!(app.search_results.contains("kJ/mol"));
    }

    fn find_calculable_thermochemistry_case(library: &str) -> (String, f64) {
        let mut app = ThermochemistryApp::new();
        app.selected_library = library.to_string();
        app.update_substances_for_library();

        let candidate_temperatures = [298.15, 400.0, 500.0, 700.0, 1000.0, 1500.0];

        for temperature in candidate_temperatures {
            for substance in &app.available_substances {
                let calculable = match app.build_calculation_subs_data(substance) {
                    Ok(mut subs_data) => {
                        subs_data
                            .extract_thermal_coeffs(substance, temperature)
                            .is_ok()
                            && subs_data
                                .calculate_thermo_properties(substance, temperature)
                                .is_ok()
                    }
                    Err(_) => false,
                };

                if calculable {
                    return (substance.clone(), temperature);
                }
            }
        }

        panic!("library should expose at least one calculable substance at the test temperature");
    }

    #[test]
    fn thermochemistry_gui_does_not_list_cea_as_a_thermochemistry_backend() {
        let app = ThermochemistryApp::new();

        assert!(
            !app.thermo_data
                .thermo_libs
                .as_ref()
                .iter()
                .any(|library| library == "CEA")
        );
        assert!(
            app.thermo_data
                .transport_libs
                .as_ref()
                .iter()
                .any(|library| library == "CEA")
        );
        assert_ne!(app.selected_library, "CEA");
    }

    #[test]
    fn thermochemistry_gui_renders_calculation_controls() {
        let app = Rc::new(RefCell::new(ThermochemistryApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_role_and_label(Role::ComboBox, "Library Source");
        harness.get_by_role_and_label(Role::Button, "Calculate Properties");
        harness.get_by_role_and_label(Role::Button, "View Plots");
        harness.get_by_role_and_label(Role::Button, "Search in NIST");
        assert!(
            harness.query_all_by_value("Pressure (Pa)").next().is_none(),
            "thermochemistry screen should not expose a pressure control when it does not use pressure"
        );
    }

    #[test]
    fn thermochemistry_gui_nasa_branch_returns_nonzero_properties() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();

        app.calculate_properties();

        assert!(app.calculated_cp.is_some_and(|value| value > 0.0));
        assert!(app.calculated_dh.is_some());
        assert!(app.calculated_ds.is_some_and(|value| value.is_finite()));
        assert!(app.search_results.contains("Thermodynamic properties"));
        assert!(app.search_results.contains("J/(mol*K)"));
        assert!(app.search_results.contains("kJ/mol"));
    }

    #[test]
    fn thermochemistry_gui_matches_canonical_subsdata_results_for_nasa() {
        let (nasa_substance, nasa_temperature) = find_calculable_thermochemistry_case("NASA_gas");
        assert_thermochemistry_matches_canonical("NASA_gas", &nasa_substance, nasa_temperature);
    }

    #[test]
    fn thermochemistry_gui_builds_canonical_subsdata_snapshot_for_calculation() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();

        let mut subs_data = app
            .build_calculation_subs_data("CO")
            .expect("thermochemistry GUI should build a canonical SubsData snapshot");

        let state = subs_data
            .get_substance_search_state("CO")
            .expect("substance should be present in canonical search state");
        assert!(state.thermo().is_resolved());
        subs_data
            .extract_thermal_coeffs("CO", 500.0)
            .expect("canonical SubsData should extract thermochemistry coefficients");
        let (cp, dh, ds) = subs_data
            .calculate_thermo_properties("CO", 500.0)
            .expect("canonical SubsData should calculate thermochemistry");
        assert!(cp.is_finite() && cp > 0.0);
        assert!(dh.is_finite());
        assert!(ds.is_finite());
    }

    #[test]
    fn thermochemistry_gui_clears_stale_state_when_library_changes() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();
        app.calculate_properties();

        assert!(app.calculated_cp.is_some());
        assert!(app.search_results.contains("Thermodynamic properties"));

        app.selected_library = app
            .thermo_data
            .thermo_libs
            .as_ref()
            .iter()
            .find(|library| *library != "NASA_gas")
            .cloned()
            .unwrap_or_else(|| "NASA_gas".to_string());
        app.update_substances_for_library();

        assert!(app.selected_substance.is_empty());
        assert!(app.calculated_cp.is_none());
        assert!(app.calculated_dh.is_none());
        assert!(app.calculated_ds.is_none());
        assert_eq!(app.search_results, "No search performed yet");
        assert!(app.plot_window.is_none());
        assert!(!app.show_plots_window);
    }

    #[test]
    fn thermochemistry_gui_invalid_calculation_clears_previous_snapshot() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();
        app.calculate_properties();

        assert!(app.calculated_cp.is_some());

        app.temperature = "invalid".to_string();
        app.calculate_properties();

        assert_eq!(app.search_results, "Temperature must be a valid number");
        assert!(app.calculated_cp.is_none());
        assert!(app.calculated_dh.is_none());
        assert!(app.calculated_ds.is_none());
    }

    #[test]
    fn thermochemistry_gui_rejects_non_finite_and_non_positive_temperature_inputs() {
        let cases = [
            ("NaN", "Temperature must be finite"),
            ("inf", "Temperature must be finite"),
            ("0", "Temperature must be positive"),
            ("-1", "Temperature must be positive"),
        ];

        for (temperature, expected_message) in cases {
            let mut app = ThermochemistryApp::new();
            app.selected_library = "NASA_gas".to_string();
            app.update_substances_for_library();
            app.selected_substance = "CO".to_string();
            app.temperature = temperature.to_string();
            app.calculate_properties();

            assert_eq!(app.search_results, expected_message);
            assert!(app.calculated_cp.is_none());
            assert!(app.calculated_dh.is_none());
            assert!(app.calculated_ds.is_none());
        }
    }

    #[test]
    fn thermochemistry_gui_requires_an_explicit_selection_for_calculation() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();
        app.substance_input = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();

        app.calculate_properties();

        assert_eq!(app.search_results, "Please select a substance first");
        assert!(app.calculated_cp.is_none());
        assert!(app.calculated_dh.is_none());
        assert!(app.calculated_ds.is_none());
    }

    #[test]
    fn thermochemistry_gui_result_snapshot_is_read_only() {
        let app = Rc::new(RefCell::new(ThermochemistryApp::new()));
        app.borrow_mut().search_results = "Frozen thermochemistry report".to_string();

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        let snapshot = harness
            .query_all_by_value("Frozen thermochemistry report")
            .find(|node| node.accesskit_node().role() == Role::MultilineTextInput)
            .expect("thermochemistry snapshot should be present as a multiline text input");
        snapshot.type_text(" edited");
        harness.run();

        assert_eq!(app.borrow().search_results, "Frozen thermochemistry report");
    }

    #[test]
    fn thermochemistry_gui_manual_invalidation_clears_plot_state() {
        let mut app = ThermochemistryApp::new();
        app.show_plots_window = true;
        app.plot_window = Some(crate::gui::gui_plot::PlotWindow::new(
            "Thermochemistry".to_string(),
            vec!["Cp".to_string()],
            nalgebra::DVector::from_vec(vec![1.0]),
            nalgebra::DMatrix::from_vec(1, 1, vec![1.0]),
        ));
        app.calculated_cp = Some(1.0);
        app.calculated_dh = Some(2.0);
        app.calculated_ds = Some(3.0);
        app.search_results = "old".to_string();

        app.invalidate_calculation_snapshot();

        assert_eq!(app.search_results, "No search performed yet");
        assert!(app.calculated_cp.is_none());
        assert!(app.calculated_dh.is_none());
        assert!(app.calculated_ds.is_none());
        assert!(app.plot_window.is_none());
        assert!(!app.show_plots_window);
    }

    #[test]
    fn thermochemistry_gui_unit_change_invalidates_stale_snapshot() {
        let mut app = ThermochemistryApp::new();
        app.calculated_cp = Some(1.0);
        app.calculated_dh = Some(2.0);
        app.calculated_ds = Some(3.0);
        app.search_results = "old".to_string();
        app.plot_status = "old plot".to_string();

        app.set_energy_unit(EnergyUnit::Cal);

        assert_eq!(app.search_results, "No search performed yet");
        assert!(app.calculated_cp.is_none());
        assert!(app.calculated_dh.is_none());
        assert!(app.calculated_ds.is_none());
        assert_eq!(app.plot_status, "No temperature range calculated yet");
        assert!(app.plot_window.is_none());
        assert!(!app.show_plots_window);
    }

    #[test]
    fn thermochemistry_gui_substance_list_is_deterministic_and_selection_aware() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();

        let mut sorted = app.available_substances.clone();
        sorted.sort();
        assert_eq!(app.available_substances, sorted);

        let first = app.available_substances.first().cloned().unwrap();
        app.search_substance_by_name(&first);
        assert_eq!(app.selected_substance, first);
    }

    #[test]
    fn thermochemistry_gui_plot_window_reports_invalid_temperature_range() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.t0 = "1000.0".to_string();
        app.tend = "500.0".to_string();
        app.show_plots_window = true;

        app.calculate_range_data();

        assert!(app.plot_status.contains("Invalid temperature range"));
        assert!(app.plot_window.is_none());
    }

    #[test]
    fn thermochemistry_gui_plot_window_clears_stale_snapshot_after_failure() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.t0 = "300.0".to_string();
        app.tend = "350.0".to_string();
        app.show_plots_window = true;

        app.calculate_range_data();
        assert!(app.plot_window.is_some());

        app.t0 = "400.0".to_string();
        app.tend = "350.0".to_string();
        app.calculate_range_data();

        assert!(app.plot_status.contains("Invalid temperature range"));
        assert!(app.plot_window.is_none());
        assert!(app.show_plots_window);
    }

    #[test]
    fn thermochemistry_gui_plot_window_reports_missing_substance_selection() {
        let mut app = ThermochemistryApp::new();
        app.show_plots_window = true;
        app.t0 = "300.0".to_string();
        app.tend = "400.0".to_string();

        app.calculate_range_data();

        assert!(app.plot_status.contains("Please select a substance first"));
        assert!(app.plot_window.is_none());
    }

    #[test]
    fn thermochemistry_gui_plot_window_builds_ready_state_for_valid_range() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.t0 = "300.0".to_string();
        app.tend = "350.0".to_string();
        app.show_plots_window = true;

        app.calculate_range_data();

        let plot_window = app
            .plot_window
            .as_ref()
            .expect("valid range should create a plot window snapshot");
        assert!(app.plot_status.contains("Temperature range ready"));
        assert!(plot_window.visible);
        assert_eq!(plot_window.arg, "Temperature (K)");
        assert_eq!(plot_window.values.len(), 3);
        assert_eq!(plot_window.t_result.len(), 100);
        assert_eq!(plot_window.y_result.ncols(), 3);
        assert_eq!(plot_window.y_result.nrows(), 100);
    }

    #[test]
    fn thermochemistry_gui_plot_window_is_removed_after_close() {
        let mut app = ThermochemistryApp::new();
        app.selected_library = "NASA_gas".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.t0 = "300.0".to_string();
        app.tend = "350.0".to_string();
        app.show_plots_window = true;

        app.calculate_range_data();
        let plot_window = app
            .plot_window
            .as_mut()
            .expect("valid range should create a plot window snapshot");
        plot_window.visible = false;

        let ctx = egui::Context::default();
        let mut open = true;
        let _ = ctx.run_ui(egui::RawInput::default(), |ui| {
            app.show(ui.ctx(), &mut open);
        });

        assert!(app.plot_window.is_none());
    }

    #[test]
    fn thermochemistry_gui_catalog_failure_becomes_visible_and_inert() {
        let mut app =
            ThermochemistryApp::from_catalog_result(Err(ThermoLibraryError::EmptyLibraryCatalog {
                path: "missing-thermo-catalog.json".to_string(),
            }));

        assert!(app.startup_error.is_some());
        assert!(app.available_substances.is_empty());
        assert!(app.calculated_cp.is_none());
        assert!(app.calculated_dh.is_none());
        assert!(app.calculated_ds.is_none());

        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();
        app.calculate_properties();

        assert!(app.search_results.contains("Thermo catalog failed to load"));
        assert!(app.calculated_cp.is_none());
        assert!(app.calculated_dh.is_none());
        assert!(app.calculated_ds.is_none());
    }
}
