//! Quality-regression tests for the reactor-IVP GUI.
//!
//! Hypotheses:
//! - the stop-condition editor must roundtrip into a solver configuration;
//! - a preview must render the selected species threshold explicitly;
//! - a real background run must surface a stop-condition termination instead of
//!   silently ignoring the selected fuel cutoff.
//!
//! Expected result:
//! the GUI stays aligned with the solver contract and exposes the stop rule as
//! a first-class user-facing control.

#[cfg(test)]
mod tests {
    use super::super::reactor_ivp_gui::{IvpGuiConfig, ReactorIvpApp};
    use RustedSciThe::numerical::LSODE2::Lsode2NativeExecutionConfig;
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::cell::RefCell;
    use std::rc::Rc;
    use std::time::Duration;

    #[test]
    fn stop_condition_roundtrips_through_preview_and_background_run() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        {
            let mut app_ref = app.borrow_mut();
            app_ref.config.stop_condition_enabled = true;
            app_ref.config.stop_condition_species_index = 0;
            app_ref.config.stop_condition_threshold = 0.2;
            app_ref.config.native_execution =
                Lsode2NativeExecutionConfig::faithful_bdf_solve(200_000, 200_000);
            app_ref
                .config
                .initial_state
                .species
                .insert("A".to_string(), 0.1);
            app_ref
                .config
                .initial_state
                .species
                .insert("B".to_string(), 0.9);
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Preview Task")
            .click_accesskit();
        harness.run();

        {
            let app_ref = app.borrow();
            let preview = app_ref
                .last_preview()
                .expect("preview should be stored after previewing the task");
            assert!(preview.summary_rows.iter().any(|row| {
                row.section == "solver"
                    && row.field == "stop_condition"
                    && row.value.contains("<=")
                    && row.value.contains("0.2")
            }));
        }

        harness
            .get_by_role_and_label(Role::Button, "Run Calculation")
            .click_accesskit();
        for _ in 0..20 {
            harness.step();
            if app
                .borrow()
                .last_status()
                .is_some_and(|status| status.contains("stopped_by_condition"))
            {
                break;
            }
            std::thread::sleep(Duration::from_millis(10));
        }

        let app_ref = app.borrow();
        assert!(
            app_ref
                .last_status()
                .is_some_and(|status| status.contains("stopped_by_condition")),
            "unexpected stop-condition status: {:?}, error: {:?}",
            app_ref.last_status(),
            app_ref.last_error()
        );
        assert!(app_ref.last_error().is_none());
    }

    #[test]
    fn gui_config_stop_condition_roundtrip_is_lossless() {
        let app = ReactorIvpApp::new();
        let mut config = IvpGuiConfig::from_task(app.task());
        config.stop_condition_enabled = true;
        config.stop_condition_species_index = 0;
        config.stop_condition_threshold = 0.125;

        let mut task = app.task.clone();
        config
            .apply_to_task(&mut task)
            .expect("typed GUI config should apply");

        let rebuilt = IvpGuiConfig::from_task(&task);
        assert!(rebuilt.stop_condition_enabled);
        assert_eq!(rebuilt.stop_condition_species_index, 0);
        assert_eq!(rebuilt.stop_condition_threshold, 0.125);
        assert!(task.solver_backend_config.stop_condition.is_some());
    }
}
