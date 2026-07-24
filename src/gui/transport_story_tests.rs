//! Story-style regression tests for the transport GUI.
//!
//! Hypotheses:
//! - The transport backend picker should be usable through real widget interaction.
//! - Selecting a backend and substance through the UI should drive a successful calculation.
//! - The visible transport report should update from the same button click path a user would use.
//!
//! Expected result:
//! the transport screen should behave like a usable calculation workflow, not a set of hidden state setters.

#[cfg(test)]
mod tests {
    use super::super::transport_gui::TransportApp;
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::cell::RefCell;
    use std::rc::Rc;

    #[test]
    fn transport_gui_story_selects_backend_and_calculates_through_the_ui() {
        let app = Rc::new(RefCell::new(TransportApp::new()));
        {
            let mut app_mut = app.borrow_mut();
            app_mut.temperature = "500.0".to_string();
            app_mut.pressure = "101325.0".to_string();
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
            .get_by_role_and_label(Role::Button, "Aramco_transport")
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
        assert_eq!(app_ref.selected_library, "Aramco_transport");
        assert_eq!(app_ref.selected_substance, "CO");
        assert!(app_ref.search_results.contains("Transport properties"));
        assert!(app_ref.search_results.contains("Thermal conductivity"));
        assert!(app_ref.search_results.contains("Viscosity"));
    }
}
