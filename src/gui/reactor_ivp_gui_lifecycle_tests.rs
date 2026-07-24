//! Lifecycle-style regression tests for the condensed reactor-IVP GUI.
//!
//! Hypotheses:
//! - preview must read from the typed editor state without mutating the source task;
//! - the typed GUI config must roundtrip through task application without losing
//!   physics, initial state, or solver backend choices;
//! - AOT-specific controls should stay attached to the typed config rather than
//!   leaking into the source task during preview.
//!
//! Expected result:
//! the GUI behaves like a structured editor with a safe preview path and a lossless
//! typed bridge to the reactor task.

#[cfg(test)]
mod tests {
    use super::super::reactor_ivp_gui::{
        IvpGuiConfig, ReactorIvpApp, ReactorIvpCalculationState, ReactorIvpGuiExecutionMode,
    };
    use crate::ReactorsIVP::solver_backend::{
        ReactorIvpExecutionBackend, ReactorIvpMatrixBackend, ReactorIvpMethod,
        ReactorIvpSymbolicBackend,
    };
    use RustedSciThe::numerical::LSODE2::{Lsode2AotProfile, Lsode2AotToolchain};
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::cell::RefCell;
    use std::rc::Rc;
    use std::time::Duration;
    use tempfile::tempdir;

    #[test]
    fn preview_does_not_mutate_source_task_and_reports_configured_route() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let original_problem_name = app.borrow().task().problem_name.clone();
        let original_execution_backend =
            app.borrow().task().solver_backend_config.execution_backend;
        let original_matrix_backend = app.borrow().task().solver_backend_config.matrix_backend;

        {
            let mut app_ref = app.borrow_mut();
            app_ref.config.problem_name = "Preview-only edit".to_string();
            app_ref.config.problem_description = "Preview should stay local".to_string();
            app_ref.config.execution_mode = ReactorIvpGuiExecutionMode::Aot;
            app_ref.config.aot_toolchain = Lsode2AotToolchain::CTcc;
            app_ref.config.aot_profile = Lsode2AotProfile::Debug;
            app_ref.config.method = ReactorIvpMethod::Adams;
            app_ref.config.symbolic_backend = ReactorIvpSymbolicBackend::ExprLegacy;
            app_ref.config.matrix_backend = ReactorIvpMatrixBackend::Banded;
            app_ref.config.x_bound = 0.02;
            app_ref.config.rtol = 1e-7;
            app_ref.config.atol = 1e-9;
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

        let app_ref = app.borrow();
        assert_eq!(app_ref.task().problem_name, original_problem_name);
        assert_eq!(
            app_ref.task().solver_backend_config.execution_backend,
            original_execution_backend
        );
        assert_eq!(
            app_ref.task().solver_backend_config.matrix_backend,
            original_matrix_backend
        );
        assert_eq!(
            app_ref.last_status(),
            Some("Task preview printed to console.")
        );
        let preview = app_ref
            .last_preview()
            .expect("preview snapshot should be stored after Preview Task");
        assert!(
            preview
                .summary_rows
                .iter()
                .any(|row| row.section == "problem"
                    && row.field == "name"
                    && row.value == "Preview-only edit")
        );
        assert!(
            preview
                .summary_rows
                .iter()
                .any(|row| row.section == "solver"
                    && row.field == "execution_backend"
                    && row.value.contains("AOT"))
        );
        assert!(
            preview
                .summary_rows
                .iter()
                .any(|row| row.section == "solver"
                    && row.field == "matrix_backend"
                    && row.value == "Banded")
        );
        assert!(
            preview
                .summary_rows
                .iter()
                .any(|row| row.section == "solver"
                    && row.field == "stop_condition"
                    && row.value.contains("<="))
        );
    }

    #[test]
    fn gui_config_roundtrip_preserves_full_typed_contract_after_multiple_edits() {
        let app = ReactorIvpApp::new();
        let mut config = IvpGuiConfig::from_task(app.task());

        config.problem_name = "Typed roundtrip".to_string();
        config.problem_description = "reactor IVP lifecycle test".to_string();
        config.physical.ro = 980.0;
        config.physical.Cp = 1110.0;
        config.physical.Lambda = 0.31;
        config.physical.m = 0.018;
        config.physical.L = 0.027;
        config.scaling.dT = 92.0;
        config.scaling.T_scale = 85.0;
        config.scaling.L = 0.027;
        config.initial_state.temperature = 612.0;
        config.initial_state.heat_flux = 0.09;
        config.initial_state.species.insert("A".to_string(), 0.78);
        config.initial_state.species.insert("B".to_string(), 0.22);
        config.thermal_effects = vec![-1.45e5];
        config.stop_condition_enabled = true;
        config.stop_condition_species_index = 0;
        config.stop_condition_threshold = 1.0e-4;
        config.method = ReactorIvpMethod::Bdf;
        config.symbolic_backend = ReactorIvpSymbolicBackend::ExprLegacy;
        config.matrix_backend = ReactorIvpMatrixBackend::Banded;
        config.execution_mode = ReactorIvpGuiExecutionMode::Aot;
        config.aot_toolchain = Lsode2AotToolchain::CTcc;
        config.aot_profile = Lsode2AotProfile::Release;
        config.x0 = 0.0;
        config.x_bound = 0.015;
        config.first_step = Some(1.0e-6);
        config.max_step = 7.5e-4;
        config.rtol = 1.0e-7;
        config.atol = 1.0e-9;

        let mut task = app.task().clone();
        config
            .apply_to_task(&mut task)
            .expect("typed GUI config should apply");

        let rebuilt = IvpGuiConfig::from_task(&task);
        assert_eq!(rebuilt.problem_name, "Typed roundtrip");
        assert_eq!(rebuilt.problem_description, "reactor IVP lifecycle test");
        assert_eq!(rebuilt.physical.ro, 980.0);
        assert_eq!(rebuilt.physical.Cp, 1110.0);
        assert_eq!(rebuilt.physical.Lambda, 0.31);
        assert_eq!(rebuilt.physical.m, 0.018);
        assert_eq!(rebuilt.physical.L, 0.027);
        assert_eq!(rebuilt.scaling.dT, 92.0);
        assert_eq!(rebuilt.scaling.T_scale, 85.0);
        assert_eq!(rebuilt.scaling.L, 0.027);
        assert_eq!(rebuilt.initial_state.temperature, 612.0);
        assert_eq!(rebuilt.initial_state.heat_flux, 0.09);
        assert_eq!(rebuilt.initial_state.species["A"], 0.78);
        assert_eq!(rebuilt.initial_state.species["B"], 0.22);
        assert_eq!(rebuilt.thermal_effects, vec![-1.45e5]);
        assert_eq!(rebuilt.method, ReactorIvpMethod::Bdf);
        assert_eq!(
            rebuilt.symbolic_backend,
            ReactorIvpSymbolicBackend::ExprLegacy
        );
        assert_eq!(rebuilt.matrix_backend, ReactorIvpMatrixBackend::Banded);
        assert_eq!(rebuilt.execution_mode, ReactorIvpGuiExecutionMode::Aot);
        assert_eq!(rebuilt.aot_toolchain, Lsode2AotToolchain::CTcc);
        assert_eq!(rebuilt.aot_profile, Lsode2AotProfile::Release);
        assert_eq!(rebuilt.x0, 0.0);
        assert_eq!(rebuilt.x_bound, 0.015);
        assert_eq!(rebuilt.first_step, Some(1.0e-6));
        assert_eq!(rebuilt.max_step, 7.5e-4);
        assert_eq!(rebuilt.rtol, 1.0e-7);
        assert_eq!(rebuilt.atol, 1.0e-9);
        assert!(rebuilt.stop_condition_enabled);
        assert_eq!(rebuilt.stop_condition_species_index, 0);
        assert_eq!(rebuilt.stop_condition_threshold, 1.0e-4);
        assert_eq!(
            task.solver_backend_config.execution_backend,
            ReactorIvpExecutionBackend::Aot {
                toolchain: Lsode2AotToolchain::CTcc,
                profile: Lsode2AotProfile::Release,
            }
        );
        assert_eq!(
            task.solver_backend_config.stop_condition,
            Some(crate::ReactorsIVP::solver_backend::ReactorIvpStopCondition::new(0, 1.0e-4))
        );
    }

    #[test]
    fn document_lifecycle_save_load_roundtrip_preserves_typed_editor_state() {
        let mut app = ReactorIvpApp::new();
        app.config.problem_name = "Lifecycle Roundtrip".to_string();
        app.config.problem_description = "IVP document lifecycle".to_string();
        app.config.physical.ro = 990.0;
        app.config.physical.Cp = 1120.0;
        app.config.physical.Lambda = 0.29;
        app.config.physical.m = 0.019;
        app.config.physical.L = 0.031;
        app.config.scaling.dT = 88.0;
        app.config.scaling.T_scale = 71.0;
        app.config.scaling.L = 0.031;
        app.config.initial_state.temperature = 640.0;
        app.config.initial_state.heat_flux = 0.11;
        app.config
            .initial_state
            .species
            .insert("A".to_string(), 0.74);
        app.config
            .initial_state
            .species
            .insert("B".to_string(), 0.26);
        app.config.thermal_effects = vec![-1.65e5];
        app.config.stop_condition_enabled = true;
        app.config.stop_condition_species_index = 0;
        app.config.stop_condition_threshold = 1.0e-4;
        app.config.method = ReactorIvpMethod::Adams;
        app.config.symbolic_backend = ReactorIvpSymbolicBackend::ExprLegacy;
        app.config.matrix_backend = ReactorIvpMatrixBackend::Banded;
        app.config.execution_mode = ReactorIvpGuiExecutionMode::Aot;
        app.config.aot_toolchain = Lsode2AotToolchain::CTcc;
        app.config.aot_profile = Lsode2AotProfile::Debug;
        app.config.x0 = 0.0;
        app.config.x_bound = 0.013;
        app.config.first_step = Some(1.0e-6);
        app.config.max_step = 6.0e-4;
        app.config.rtol = 1.0e-7;
        app.config.atol = 1.0e-9;

        app.refresh_document_lifecycle_dirty_state();
        assert!(
            app.document_lifecycle.dirty,
            "editing the config should mark the document dirty"
        );

        let tempdir = tempdir().expect("temporary directory should be created");
        let save_path = tempdir.path().join("reactor_ivp_document.txt");
        let save_result = app.save_document(save_path.clone());
        assert!(matches!(
            save_result,
            crate::gui::document_lifecycle::GuiFileOperationResult::Saved { .. }
        ));
        assert_eq!(app.current_file_path.as_deref(), Some(save_path.as_path()));
        assert!(!app.document_lifecycle.dirty);

        let saved_text = std::fs::read_to_string(&save_path).expect("saved document should exist");
        assert!(saved_text.contains("process_conditions"));
        assert!(saved_text.contains("solver_options"));
        assert!(saved_text.contains("stop_condition"));

        let mut loaded = ReactorIvpApp::new();
        let load_result = loaded.load_document_from_path(save_path.clone());
        assert!(matches!(
            load_result,
            crate::gui::document_lifecycle::GuiFileOperationResult::Loaded { .. }
        ));
        assert_eq!(
            loaded.current_file_path.as_deref(),
            Some(save_path.as_path())
        );
        assert_eq!(loaded.config.problem_name, "Lifecycle Roundtrip");
        assert_eq!(loaded.config.problem_description, "IVP document lifecycle");
        assert_eq!(loaded.config.physical.ro, 990.0);
        assert_eq!(loaded.config.physical.Cp, 1120.0);
        assert_eq!(loaded.config.physical.Lambda, 0.29);
        assert_eq!(loaded.config.physical.m, 0.019);
        assert_eq!(loaded.config.physical.L, 0.031);
        assert_eq!(loaded.config.scaling.dT, 88.0);
        assert_eq!(loaded.config.scaling.T_scale, 71.0);
        assert_eq!(loaded.config.scaling.L, 0.031);
        assert_eq!(loaded.config.initial_state.temperature, 640.0);
        assert_eq!(loaded.config.initial_state.heat_flux, 0.11);
        assert_eq!(loaded.config.initial_state.species["A"], 0.74);
        assert_eq!(loaded.config.initial_state.species["B"], 0.26);
        assert_eq!(loaded.config.thermal_effects, vec![-1.65e5]);
        assert!(loaded.config.stop_condition_enabled);
        assert_eq!(loaded.config.stop_condition_species_index, 0);
        assert_eq!(loaded.config.stop_condition_threshold, 1.0e-4);
        assert_eq!(loaded.config.method, ReactorIvpMethod::Adams);
        assert_eq!(
            loaded.config.symbolic_backend,
            ReactorIvpSymbolicBackend::ExprLegacy
        );
        assert_eq!(
            loaded.config.matrix_backend,
            ReactorIvpMatrixBackend::Banded
        );
        assert_eq!(
            loaded.config.execution_mode,
            ReactorIvpGuiExecutionMode::Aot
        );
        assert_eq!(loaded.config.aot_toolchain, Lsode2AotToolchain::CTcc);
        assert_eq!(loaded.config.aot_profile, Lsode2AotProfile::Debug);
        assert_eq!(loaded.config.x_bound, 0.013);
        assert_eq!(loaded.config.first_step, Some(1.0e-6));
        assert_eq!(loaded.config.max_step, 6.0e-4);
        assert_eq!(loaded.config.rtol, 1.0e-7);
        assert_eq!(loaded.config.atol, 1.0e-9);
        assert!(!loaded.document_lifecycle.dirty);
    }

    #[test]
    fn document_lifecycle_queues_confirmation_when_replacing_dirty_work() {
        let mut app = ReactorIvpApp::new();
        let tempdir = tempdir().expect("temporary directory should be created");
        let save_path = tempdir.path().join("reactor_ivp_document.txt");

        let _ = app.save_document(save_path.clone());
        app.config.problem_name = "Unsaved change".to_string();
        app.refresh_document_lifecycle_dirty_state();
        assert!(app.document_lifecycle.dirty);

        let queued = app.request_document_load(save_path.clone());
        assert!(!queued);
        assert_eq!(app.pending_document_load.as_ref(), Some(&save_path));
    }

    #[test]
    fn solved_state_survives_save_but_reload_starts_clean() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        app.borrow_mut()
            .start_test_calculation_worker(Duration::from_millis(25), true);
        harness.step();

        for _ in 0..20 {
            if !matches!(
                app.borrow().calculation_state(),
                ReactorIvpCalculationState::Running { .. }
            ) {
                break;
            }
            std::thread::sleep(Duration::from_millis(10));
            harness.step();
        }

        let solved_preview_rows = app
            .borrow()
            .last_preview()
            .expect("successful solve should publish a preview")
            .summary_rows
            .len();
        let solved_conservation_rows = app
            .borrow()
            .task()
            .latest_conservation_report()
            .expect("successful solve should publish a conservation report")
            .sum_of_mass_fractions
            .len();

        let tempdir = tempdir().expect("temporary directory should be created");
        let save_path = tempdir.path().join("reactor_ivp_solved_roundtrip.txt");
        let save_result = app.borrow_mut().save_document(save_path.clone());
        assert!(matches!(
            save_result,
            crate::gui::document_lifecycle::GuiFileOperationResult::Saved { .. }
        ));

        let app_ref = app.borrow();
        assert!(matches!(
            app_ref.calculation_state(),
            ReactorIvpCalculationState::Completed
        ));
        assert_eq!(
            app_ref
                .last_preview()
                .expect("saved solved app should keep the preview snapshot")
                .summary_rows
                .len(),
            solved_preview_rows
        );
        assert_eq!(
            app_ref
                .task()
                .latest_conservation_report()
                .expect("saved solved app should keep the conservation report")
                .sum_of_mass_fractions
                .len(),
            solved_conservation_rows
        );
        drop(app_ref);

        let mut loaded = ReactorIvpApp::new();
        let load_result = loaded.load_document_from_path(save_path.clone());
        assert!(matches!(
            load_result,
            crate::gui::document_lifecycle::GuiFileOperationResult::Loaded { .. }
        ));
        assert!(matches!(
            loaded.calculation_state(),
            ReactorIvpCalculationState::Idle
        ));
        assert!(loaded.last_preview().is_none());
        assert!(loaded.task().last_solve_snapshot.is_none());
        assert!(loaded.task().last_conservation_report.is_none());
        assert_eq!(
            loaded.current_file_path.as_deref(),
            Some(save_path.as_path())
        );
    }

    #[test]
    fn worker_lifecycle_reaches_completed_state_after_background_run() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        app.borrow_mut()
            .start_test_calculation_worker(Duration::from_millis(25), true);
        harness.step();

        assert!(matches!(
            app.borrow().calculation_state(),
            ReactorIvpCalculationState::Running { .. }
        ));

        for _ in 0..20 {
            if !matches!(
                app.borrow().calculation_state(),
                ReactorIvpCalculationState::Running { .. }
            ) {
                break;
            }
            std::thread::sleep(Duration::from_millis(10));
            harness.step();
        }

        let app_ref = app.borrow();
        assert!(matches!(
            app_ref.calculation_state(),
            ReactorIvpCalculationState::Completed
        ));
        assert!(
            app_ref
                .last_status()
                .is_some_and(|status| status.contains("Calculation completed"))
        );
        assert!(app_ref.last_error().is_none());
        assert!(app_ref.last_preview().is_some_and(
            |snapshot| !snapshot.summary_rows.is_empty() && !snapshot.equation_rows.is_empty()
        ));
    }

    #[test]
    fn worker_lifecycle_failure_keeps_previous_preview_and_reports_error() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
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

        let prior_preview_rows = app
            .borrow()
            .last_preview()
            .expect("preview should exist before failure")
            .summary_rows
            .len();

        app.borrow_mut()
            .start_test_calculation_worker(Duration::from_millis(25), false);
        harness.step();

        for _ in 0..20 {
            if !matches!(
                app.borrow().calculation_state(),
                ReactorIvpCalculationState::Running { .. }
            ) {
                break;
            }
            std::thread::sleep(Duration::from_millis(10));
            harness.step();
        }

        let app_ref = app.borrow();
        assert!(matches!(
            app_ref.calculation_state(),
            ReactorIvpCalculationState::Failed
        ));
        assert!(
            app_ref
                .last_error()
                .is_some_and(|error| error.contains("synthetic reactor IVP worker failure"))
        );
        assert_eq!(
            app_ref
                .last_preview()
                .expect("previous preview should be retained after worker failure")
                .summary_rows
                .len(),
            prior_preview_rows
        );
    }

    #[test]
    fn worker_lifecycle_discards_stale_results_after_editor_changes() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
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

        let prior_preview_rows = app
            .borrow()
            .last_preview()
            .expect("preview should exist before stale result test")
            .summary_rows
            .len();

        app.borrow_mut()
            .start_test_calculation_worker(Duration::from_millis(25), true);
        app.borrow_mut().config.problem_name = "Changed while running".to_string();
        harness.step();

        for _ in 0..20 {
            if !matches!(
                app.borrow().calculation_state(),
                ReactorIvpCalculationState::Running { .. }
            ) {
                break;
            }
            std::thread::sleep(Duration::from_millis(10));
            harness.step();
        }

        let app_ref = app.borrow();
        assert!(matches!(
            app_ref.calculation_state(),
            ReactorIvpCalculationState::Failed
        ));
        assert!(
            app_ref
                .last_status()
                .is_some_and(|status| status.contains("ignored because the editor changed")),
            "unexpected stale-result status: {:?}, error: {:?}",
            app_ref.last_status(),
            app_ref.last_error()
        );
        assert!(
            app_ref
                .last_error()
                .is_some_and(|error| error.contains("Stale calculation result discarded"))
        );
        assert_eq!(
            app_ref
                .last_preview()
                .expect("previous preview should be retained after stale result rejection")
                .summary_rows
                .len(),
            prior_preview_rows
        );
    }

    #[test]
    fn replace_task_resets_editor_state_and_preserves_a_safe_switch_boundary() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
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

        let mut replacement = app.borrow().task().clone();
        replacement.problem_name = Some("Replacement problem".to_string());
        replacement.problem_description = Some("New task snapshot".to_string());

        app.borrow_mut()
            .replace_task(replacement)
            .expect("task replacement should succeed when idle");

        let app_ref = app.borrow();
        assert_eq!(app_ref.last_status(), Some("Task replaced."));
        assert!(app_ref.last_error().is_none());
        assert!(app_ref.last_preview().is_none());
        assert_eq!(
            app_ref.config.problem_name,
            "Replacement problem".to_string()
        );
        assert_eq!(
            app_ref.config.problem_description,
            "New task snapshot".to_string()
        );
        assert_eq!(
            app_ref.task().problem_name.as_deref(),
            Some("Replacement problem")
        );
        assert_eq!(
            app_ref.task().problem_description.as_deref(),
            Some("New task snapshot")
        );
    }

    #[test]
    fn replace_task_rejects_active_calculation() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        app.borrow_mut()
            .start_test_calculation_worker(Duration::from_millis(25), true);
        harness.step();

        let mut replacement = app.borrow().task().clone();
        replacement.problem_name = Some("Should not replace while busy".to_string());

        let error = app
            .borrow_mut()
            .replace_task(replacement)
            .expect_err("active calculation should reject task replacement");

        assert!(
            error
                .to_string()
                .contains("cannot replace the active reactor IVP task")
        );
        assert!(matches!(
            app.borrow().calculation_state(),
            ReactorIvpCalculationState::Running { .. }
        ));
    }
}
