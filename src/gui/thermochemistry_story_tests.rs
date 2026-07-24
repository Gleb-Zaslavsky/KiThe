//! Story-style regression tests for the thermochemistry plot workflow.
//!
//! Hypotheses:
//! - Clicking `View Plots` should open the temperature-range controls window.
//! - Clicking `Calculate Range` from that window should build a ready plot snapshot.
//! - Invalid temperature ranges should surface a visible error and not leave stale plot data behind.
//! - Real widget interaction should still be enough to select a substance and run a calculation.
//!
//! Expected result:
//! the plot workflow remains a real user-facing path rather than a hidden backend-only contract.

#[cfg(test)]
mod tests {
    use super::super::thermochemistry_gui::ThermochemistryApp;
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::cell::RefCell;
    use std::rc::Rc;

    #[test]
    fn thermochemistry_gui_story_view_plots_opens_the_range_window() {
        let app = Rc::new(RefCell::new(ThermochemistryApp::new()));
        {
            let mut app_mut = app.borrow_mut();
            app_mut.selected_library = "NASA_gas".to_string();
            app_mut.update_substances_for_library();
            app_mut.selected_substance = "CO".to_string();
            app_mut.t0 = "300.0".to_string();
            app_mut.tend = "350.0".to_string();
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "View Plots")
            .click();
        harness.run();

        assert!(app.borrow().show_plots_window);
        harness.get_by_role_and_label(Role::Button, "Calculate Range");
    }

    #[test]
    fn thermochemistry_gui_story_calculates_plot_snapshot_through_the_ui() {
        let app = Rc::new(RefCell::new(ThermochemistryApp::new()));
        {
            let mut app_mut = app.borrow_mut();
            app_mut.selected_library = "NASA_gas".to_string();
            app_mut.update_substances_for_library();
            app_mut.selected_substance = "CO".to_string();
            app_mut.t0 = "300.0".to_string();
            app_mut.tend = "350.0".to_string();
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "View Plots")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Calculate Range")
            .click();
        harness.run();

        let app_ref = app.borrow();
        assert!(app_ref.show_plots_window);
        assert!(app_ref.plot_status.contains("Temperature range ready"));
        let plot_window = app_ref
            .plot_window
            .as_ref()
            .expect("successful range calculation should build a plot snapshot");
        assert!(plot_window.visible);
        assert_eq!(plot_window.values.len(), 3);
    }

    #[test]
    fn thermochemistry_gui_story_reports_invalid_plot_range_through_the_ui() {
        let app = Rc::new(RefCell::new(ThermochemistryApp::new()));
        {
            let mut app_mut = app.borrow_mut();
            app_mut.selected_library = "NASA_gas".to_string();
            app_mut.update_substances_for_library();
            app_mut.selected_substance = "CO".to_string();
            app_mut.t0 = "1000.0".to_string();
            app_mut.tend = "350.0".to_string();
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "View Plots")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Calculate Range")
            .click();
        harness.run();

        let app_ref = app.borrow();
        assert!(app_ref.show_plots_window);
        assert!(app_ref.plot_status.contains("Invalid temperature range"));
        assert!(app_ref.plot_window.is_none());
    }

    #[test]
    fn thermochemistry_gui_story_selects_backend_and_calculates_through_the_ui() {
        let app = Rc::new(RefCell::new(ThermochemistryApp::new()));
        {
            let mut app_mut = app.borrow_mut();
            app_mut.temperature = "500.0".to_string();
        }

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
        harness
            .get_by_role_and_label(Role::Button, "NASA_gas")
            .click_accesskit();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "CO")
            .click_accesskit();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Calculate Properties")
            .click_accesskit();
        harness.run();

        let app_ref = app.borrow();
        assert_eq!(app_ref.selected_library, "NASA_gas");
        assert_eq!(app_ref.selected_substance, "CO");
        assert!(app_ref.search_results.contains("Thermodynamic properties"));
        assert!(app_ref.search_results.contains("Cp"));
        assert!(app_ref.search_results.contains("dH"));
        assert!(app_ref.search_results.contains("dS"));
    }
}
