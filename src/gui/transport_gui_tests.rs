//! Story-style regression tests for the transport GUI.
//!
//! Hypotheses:
//! - CEA is exposed through the transport backend picker, not through the thermochemistry picker.
//! - The CEA transport path returns real values instead of the old zero placeholder.
//! - The Aramco transport path can derive the missing supporting inputs from existing thermo data.
//!
//! Expected result:
//! the transport screen should behave like a real calculation surface, not a dead-end form.

#[cfg(test)]
mod tests {
    use super::super::transport_gui::TransportApp;
    use crate::Thermodynamics::DBhandlers::transport_api::{LambdaUnit, ViscosityUnit};
    use crate::Thermodynamics::thermo_lib_api::ThermoLibraryError;
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

    fn assert_transport_matches_canonical(
        library: &str,
        substance: &str,
        temperature: f64,
        pressure: f64,
    ) {
        let mut app = TransportApp::new();
        app.selected_library = library.to_string();
        app.update_substances_for_library();
        app.selected_substance = substance.to_string();
        app.temperature = temperature.to_string();
        app.pressure = pressure.to_string();
        app.set_lambda_unit(LambdaUnit::WPerMK);
        app.set_viscosity_unit(ViscosityUnit::PaS);

        let mut canonical = app
            .build_calculation_subs_data(substance, pressure)
            .expect("transport GUI should build canonical SubsData for comparison");
        canonical
            .extract_thermal_coeffs(substance, temperature)
            .expect("canonical transport calculation should extract thermochemistry coefficients");
        let (cp, _, _) = canonical
            .calculate_thermo_properties(substance, temperature)
            .expect("canonical transport calculation should resolve Cp");
        canonical
            .extract_transport_coeffs(substance, temperature)
            .expect("canonical transport calculation should extract transport coefficients");
        let (expected_lambda, expected_viscosity) = canonical
            .calculate_transport_properties(substance, temperature, Some(cp), None)
            .expect("canonical transport calculation should succeed");

        app.calculate_properties();

        let report_lambda = parse_report_value(&app.search_results, "Thermal conductivity")
            .expect("transport report should contain thermal conductivity");
        let report_viscosity = parse_report_value(&app.search_results, "Viscosity")
            .expect("transport report should contain viscosity");
        let expected_lambda_display = format!("{:.5}", expected_lambda)
            .parse::<f64>()
            .expect("rounded transport lambda should parse");
        let expected_viscosity_display = format!("{:.8}", expected_viscosity)
            .parse::<f64>()
            .expect("rounded transport viscosity should parse");

        assert_eq!(app.calculated_lambda, Some(expected_lambda));
        assert_eq!(app.calculated_viscosity, Some(expected_viscosity));
        assert_eq!(report_lambda, expected_lambda_display);
        assert_eq!(report_viscosity, expected_viscosity_display);
        assert!(app.search_results.contains("W/m*K"));
        assert!(app.search_results.contains("Pa*s"));
    }

    #[test]
    fn transport_gui_dropdown_exposes_transport_backends() {
        let app = Rc::new(RefCell::new(TransportApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Library Source")
            .click_accesskit();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "CEA");
        harness.get_by_role_and_label(Role::Button, "Aramco_transport");
    }

    #[test]
    fn transport_gui_cea_branch_returns_nonzero_properties() {
        let mut app = TransportApp::new();
        app.selected_library = "CEA".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();

        app.calculate_properties();

        assert!(
            app.calculated_lambda.is_some_and(|value| value > 0.0),
            "{}",
            app.search_results
        );
        assert!(
            app.calculated_viscosity.is_some_and(|value| value > 0.0),
            "{}",
            app.search_results
        );
        assert!(app.search_results.contains("Thermal conductivity"));
        assert!(app.search_results.contains("Viscosity"));
        assert!(app.search_results.contains("uW/m*K"));
        assert!(app.search_results.contains("uPa*s"));
    }

    #[test]
    fn transport_gui_aramco_branch_uses_derived_inputs_and_returns_nonzero_values() {
        let mut app = TransportApp::new();
        app.selected_library = "Aramco_transport".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();

        app.calculate_properties();

        assert!(
            app.calculated_lambda.is_some_and(|value| value > 0.0),
            "{}",
            app.search_results
        );
        assert!(
            app.calculated_viscosity.is_some_and(|value| value > 0.0),
            "{}",
            app.search_results
        );
        assert!(app.search_results.contains("Transport properties"));
        assert!(app.search_results.contains("uW/m*K"));
        assert!(app.search_results.contains("uPa*s"));
    }

    #[test]
    fn transport_gui_matches_canonical_subsdata_results_for_cea_and_aramco() {
        assert_transport_matches_canonical("CEA", "CO", 500.0, 101325.0);
        assert_transport_matches_canonical("Aramco_transport", "CO", 500.0, 101325.0);
    }

    #[test]
    fn transport_gui_builds_canonical_subsdata_snapshot_for_calculation() {
        let mut app = TransportApp::new();
        app.selected_library = "Aramco_transport".to_string();
        app.update_substances_for_library();

        let subs_data = app
            .build_calculation_subs_data("CO", 101325.0)
            .expect("transport GUI should build a canonical SubsData snapshot");

        let state = subs_data
            .get_substance_search_state("CO")
            .expect("substance should be present in canonical search state");
        assert!(state.transport().is_resolved());
        assert!(state.thermo().is_resolved());
        assert_eq!(subs_data.pressure(), Some(101325.0));
        assert_eq!(subs_data.molar_mass_unit(), Some("g/mol"));
        assert!(subs_data.molar_mass_by_substance.get("CO").is_some());
    }

    #[test]
    fn transport_gui_clears_stale_state_when_library_changes() {
        let mut app = TransportApp::new();
        app.selected_library = "Aramco_transport".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();
        app.calculate_properties();

        assert!(app.calculated_lambda.is_some());
        assert!(app.search_results.contains("Transport properties"));

        app.selected_library = app
            .transport_data
            .transport_libs
            .as_ref()
            .iter()
            .find(|library| *library != "Aramco_transport")
            .cloned()
            .unwrap_or_else(|| "Aramco_transport".to_string());
        app.update_substances_for_library();

        assert!(app.selected_substance.is_empty());
        assert!(app.calculated_lambda.is_none());
        assert!(app.calculated_viscosity.is_none());
        assert_eq!(app.search_results, "No search performed yet");
    }

    #[test]
    fn transport_gui_invalid_calculation_clears_previous_snapshot() {
        let mut app = TransportApp::new();
        app.selected_library = "Aramco_transport".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();
        app.calculate_properties();

        assert!(app.calculated_lambda.is_some());

        app.temperature = "invalid".to_string();
        app.calculate_properties();

        assert_eq!(app.search_results, "Temperature must be a valid number");
        assert!(app.calculated_lambda.is_none());
        assert!(app.calculated_viscosity.is_none());
    }

    #[test]
    fn transport_gui_rejects_non_finite_and_non_positive_numeric_inputs() {
        let cases = [
            ("NaN", "101325.0", "Temperature must be finite"),
            ("inf", "101325.0", "Temperature must be finite"),
            ("0", "101325.0", "Temperature must be positive"),
            ("-1", "101325.0", "Temperature must be positive"),
            ("500.0", "NaN", "Pressure must be finite"),
            ("500.0", "inf", "Pressure must be finite"),
            ("500.0", "0", "Pressure must be positive"),
            ("500.0", "-1", "Pressure must be positive"),
        ];

        for (temperature, pressure, expected_message) in cases {
            let mut app = TransportApp::new();
            app.selected_library = "Aramco_transport".to_string();
            app.update_substances_for_library();
            app.selected_substance = "CO".to_string();
            app.temperature = temperature.to_string();
            app.pressure = pressure.to_string();
            app.calculate_properties();

            assert_eq!(app.search_results, expected_message);
            assert!(app.calculated_lambda.is_none());
            assert!(app.calculated_viscosity.is_none());
        }
    }

    #[test]
    fn transport_gui_requires_an_explicit_selection_for_calculation() {
        let mut app = TransportApp::new();
        app.selected_library = "Aramco_transport".to_string();
        app.update_substances_for_library();
        app.substance_input = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();

        app.calculate_properties();

        assert_eq!(app.search_results, "Please select a substance first");
        assert!(app.calculated_lambda.is_none());
        assert!(app.calculated_viscosity.is_none());
    }

    #[test]
    fn transport_gui_result_snapshot_is_read_only() {
        let app = Rc::new(RefCell::new(TransportApp::new()));
        app.borrow_mut().search_results = "Frozen transport report".to_string();

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        let snapshot = harness
            .query_all_by_value("Frozen transport report")
            .find(|node| node.accesskit_node().role() == Role::MultilineTextInput)
            .expect("transport snapshot should be present as a multiline text input");
        snapshot.type_text(" edited");
        harness.run();

        assert_eq!(app.borrow().search_results, "Frozen transport report");
    }

    #[test]
    fn transport_gui_manual_invalidation_clears_result_snapshot() {
        let mut app = TransportApp::new();
        app.calculated_lambda = Some(1.0);
        app.calculated_viscosity = Some(2.0);
        app.search_results = "old".to_string();

        app.invalidate_calculation_snapshot();

        assert_eq!(app.search_results, "No search performed yet");
        assert!(app.calculated_lambda.is_none());
        assert!(app.calculated_viscosity.is_none());
    }

    #[test]
    fn transport_gui_unit_change_updates_displayed_values() {
        let mut app = TransportApp::new();
        app.selected_library = "Aramco_transport".to_string();
        app.update_substances_for_library();
        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();
        app.set_lambda_unit(crate::Thermodynamics::DBhandlers::transport_api::LambdaUnit::WPerMK);
        app.set_viscosity_unit(
            crate::Thermodynamics::DBhandlers::transport_api::ViscosityUnit::PaS,
        );
        app.calculate_properties();

        let raw_lambda = app
            .calculated_lambda
            .expect("successful run should cache the raw lambda value");
        let raw_viscosity = app
            .calculated_viscosity
            .expect("successful run should cache the raw viscosity value");
        let lambda_before = parse_report_value(&app.search_results, "Thermal conductivity")
            .expect("report should contain thermal conductivity");
        let viscosity_before = parse_report_value(&app.search_results, "Viscosity")
            .expect("report should contain viscosity");

        app.set_lambda_unit(crate::Thermodynamics::DBhandlers::transport_api::LambdaUnit::MWPerMK);
        app.set_viscosity_unit(
            crate::Thermodynamics::DBhandlers::transport_api::ViscosityUnit::MKPaS,
        );

        let lambda_after = parse_report_value(&app.search_results, "Thermal conductivity")
            .expect("updated report should contain thermal conductivity");
        let viscosity_after = parse_report_value(&app.search_results, "Viscosity")
            .expect("updated report should contain viscosity");

        assert!(
            lambda_after > lambda_before * 900.0 && lambda_after < lambda_before * 1_100.0,
            "lambda should be rescaled when the display unit changes: before={lambda_before}, after={lambda_after}, report={}",
            app.search_results
        );
        assert!(
            viscosity_after > viscosity_before * 900_000.0
                && viscosity_after < viscosity_before * 1_100_000.0,
            "viscosity should be rescaled when the display unit changes: before={viscosity_before}, after={viscosity_after}, report={}",
            app.search_results
        );
        assert_eq!(app.calculated_lambda, Some(raw_lambda));
        assert_eq!(app.calculated_viscosity, Some(raw_viscosity));
        assert!(app.search_results.contains("mW/m*K"));
        assert!(app.search_results.contains("uPa*s"));
    }

    #[test]
    fn transport_gui_substance_list_is_deterministic_and_selection_aware() {
        let mut app = TransportApp::new();
        app.selected_library = "Aramco_transport".to_string();
        app.update_substances_for_library();

        let mut sorted = app.available_substances.clone();
        sorted.sort();
        assert_eq!(app.available_substances, sorted);

        let first = app.available_substances.first().cloned().unwrap();
        app.search_substance_by_name(&first);
        assert_eq!(app.selected_substance, first);
    }

    #[test]
    fn transport_gui_catalog_failure_becomes_visible_and_inert() {
        let mut app =
            TransportApp::from_catalog_result(Err(ThermoLibraryError::EmptyLibraryCatalog {
                path: "missing-transport-catalog.json".to_string(),
            }));

        assert!(app.startup_error.is_some());
        assert!(app.available_substances.is_empty());
        assert!(app.calculated_lambda.is_none());
        assert!(app.calculated_viscosity.is_none());

        app.selected_substance = "CO".to_string();
        app.temperature = "500.0".to_string();
        app.pressure = "101325.0".to_string();
        app.calculate_properties();

        assert!(app.search_results.contains("Thermo catalog failed to load"));
        assert!(app.calculated_lambda.is_none());
        assert!(app.calculated_viscosity.is_none());
    }
}
