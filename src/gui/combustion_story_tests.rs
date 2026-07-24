//! Story tests for the combustion BVP GUI.
//!
//! These tests treat the BVP screen as a user-facing workflow rather than a
//! bag of controls. Each scenario checks a concrete hypothesis:
//!
//! - A fresh BVP preset should expose the structured layout immediately.
//! - Canonical solver defaults should be visible and stable on first render.
//! - AOT-oriented backend values should expose the AOT-specific controls.
//! - Lambdify and AOT should remain mutually exclusive document contracts.
//! - AOT should require explicit consent before an external toolchain can run.
//! - A task document must not be able to select an arbitrary executable as a compiler.
//! - Structured BVP data deletion should require confirmation, while optional raw data remains editable.
//! - Legacy documents should still render when older spellings are present.
//! - A solved BVP should still hand postprocessing state through to the GUI.
//! - The plot window should appear only when `gui_plot` is enabled.
//!
//! Expected result:
//! the BVP screen should behave like a structured form with a compatibility
//! fallback, not like a raw document editor with a few extra labels on top.

#[cfg(test)]
mod tests {
    use super::super::combustion::{
        CalculationState, CombustionApp, DocumentDeleteRequest, GuiFileOperationResult,
        ProblemsEnum, apply_document_deletion,
    };
    use RustedSciThe::command_interpreter::task_parser::{DocumentMap, Value};
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::cell::RefCell;
    use std::path::PathBuf;
    use std::rc::Rc;
    use std::time::Duration;

    fn first_value<'a>(document: &'a DocumentMap, section: &str, field: &str) -> &'a Value {
        document
            .get(section)
            .and_then(|section_map| section_map.get(field))
            .and_then(|slot| slot.as_ref())
            .and_then(|values| values.first())
            .expect("expected BVP field to exist")
    }

    fn as_string(value: &Value) -> &str {
        match value {
            Value::String(value) => value,
            other => panic!("expected Value::String, got {other:?}"),
        }
    }

    fn configure_aot(app: &mut CombustionApp, compiler: &str) {
        let solver_settings = app
            .document
            .get_mut("solver_settings")
            .expect("solver_settings section must exist in the BVP template");
        solver_settings.insert(
            "generated_backend".to_string(),
            Some(vec![Value::String("banded_aot_tcc".to_string())]),
        );
        solver_settings.insert(
            "backend_policy".to_string(),
            Some(vec![Value::String("aot_only".to_string())]),
        );
        solver_settings.insert(
            "aot_codegen_backend".to_string(),
            Some(vec![Value::String("C".to_string())]),
        );
        solver_settings.insert(
            "aot_c_compiler".to_string(),
            Some(vec![Value::String(compiler.to_string())]),
        );
        solver_settings.insert(
            "aot_build_policy".to_string(),
            Some(vec![Value::String("build_if_missing".to_string())]),
        );
    }

    /// Toggle the plot handoff flag without mutating the rest of the BVP task.
    fn set_gui_plot(app: &mut CombustionApp, enabled: bool) {
        let postprocessing = app
            .document
            .get_mut("postprocessing")
            .expect("postprocessing section must exist in the BVP template");
        postprocessing.insert("gui_plot".to_string(), Some(vec![Value::Boolean(enabled)]));
    }

    /// Load the known BVP story document used by the GUI regression suite.
    fn load_problem_test_document(app: &mut CombustionApp) {
        let result = app.load_document_from_path(PathBuf::from("problem_test.txt"));
        assert!(
            matches!(result, GuiFileOperationResult::Loaded { .. }),
            "problem_test.txt should load as a valid BVP document"
        );
    }

    fn wait_for_calculation_with_limit(app: &mut CombustionApp, attempts: usize, sleep: Duration) {
        let ctx = egui::Context::default();
        for _ in 0..attempts {
            app.poll_calculation_worker(&ctx);
            if !matches!(
                app.calculation_state,
                CalculationState::Running { .. } | CalculationState::Cancelling { .. }
            ) {
                return;
            }
            std::thread::sleep(sleep);
        }
        panic!("calculation worker did not finish within the test timeout");
    }

    fn wait_for_calculation(app: &mut CombustionApp) {
        wait_for_calculation_with_limit(app, 200, Duration::from_millis(5));
    }

    #[test]
    fn default_bvp_preset_exposes_structured_screen_and_canonical_solver_defaults() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();

        for label in [
            "Physics",
            "Solver backend",
            "Initial guess",
            "Postprocessing",
            "Legacy raw document",
            "Execution and assembly",
            "Residuals and Jacobian Backend",
        ] {
            harness.get_by_label(label);
        }

        let app_ref = app.borrow();
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "generated_backend"
            )),
            "banded_lambdify"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "symbolic_backend"
            )),
            "AtomView"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "matrix_backend"
            )),
            "Banded"
        );
    }

    #[test]
    fn structured_sections_are_reported_but_not_recreated_by_render() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            app_ref.document.remove("initial_guess");
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("Missing canonical section: initial_guess");

        let app_ref = app.borrow();
        assert!(
            !app_ref.document.contains_key("initial_guess"),
            "rendering should not recreate missing canonical sections"
        );
    }

    #[test]
    fn aot_oriented_backend_values_survive_rendering() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            let solver_settings = app_ref
                .document
                .get_mut("solver_settings")
                .expect("solver_settings section must exist in the BVP template");

            solver_settings.insert(
                "generated_backend".to_string(),
                Some(vec![Value::String("banded_aot_tcc".to_string())]),
            );
            solver_settings.insert(
                "backend_policy".to_string(),
                Some(vec![Value::String("prefer_aot_then_lambdify".to_string())]),
            );
            solver_settings.insert(
                "aot_compile_preset".to_string(),
                Some(vec![Value::String("production".to_string())]),
            );
            solver_settings.insert(
                "aot_execution_policy".to_string(),
                Some(vec![Value::String("parallel".to_string())]),
            );
            solver_settings.insert(
                "aot_build_policy".to_string(),
                Some(vec![Value::String("require_prebuilt".to_string())]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Residuals and Jacobian Backend");
        harness.get_by_label("AOT codegen backend");
        harness.get_by_label("C compiler");
        harness.get_by_label("AOT build policy");
        harness.get_by_label("AOT build profile");
        harness.get_by_label("AOT compile preset");
        harness.get_by_label("AOT execution policy");

        let app_ref = app.borrow();
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "generated_backend"
            )),
            "banded_aot_tcc"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "backend_policy"
            )),
            "aot_only"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "aot_execution_policy"
            )),
            "parallel"
        );
    }

    #[test]
    fn sparse_compatibility_backend_values_survive_rendering() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            let solver_settings = app_ref
                .document
                .get_mut("solver_settings")
                .expect("solver_settings section must exist in the BVP template");

            solver_settings.insert(
                "generated_backend".to_string(),
                Some(vec![Value::String("sparse_lambdify".to_string())]),
            );
            solver_settings.insert(
                "method".to_string(),
                Some(vec![Value::String("Sparse".to_string())]),
            );
            solver_settings.insert(
                "matrix_backend".to_string(),
                Some(vec![Value::String("Sparse".to_string())]),
            );
            solver_settings.insert(
                "backend_policy".to_string(),
                Some(vec![Value::String("prefer_aot_then_numeric".to_string())]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Residuals and Jacobian Backend");
        harness.get_by_label("Linear algebra backend");

        let app_ref = app.borrow();
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "generated_backend"
            )),
            "sparse_lambdify"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "matrix_backend"
            )),
            "Sparse"
        );
        assert_eq!(
            as_string(first_value(&app_ref.document, "solver_settings", "method")),
            "Sparse"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "backend_policy"
            )),
            "lambdify_only"
        );
        assert!(
            !app_ref
                .document
                .get("solver_settings")
                .is_some_and(|section| section.contains_key("aot_build_policy"))
        );
    }

    #[test]
    fn legacy_backend_document_still_renders_with_raw_compatibility_data() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            let solver_settings = app_ref
                .document
                .get_mut("solver_settings")
                .expect("solver_settings section must exist in the BVP template");
            solver_settings.remove("dont_save_log");
            solver_settings.insert(
                "dont_save_logs".to_string(),
                Some(vec![Value::Boolean(true)]),
            );
            solver_settings.remove("generated_backend");
            solver_settings.remove("matrix_backend");
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Solver backend");
        harness.get_by_label("Legacy raw document");

        let app_ref = app.borrow();
        assert!(
            app_ref
                .document
                .get("solver_settings")
                .and_then(|section| section.get("dont_save_logs"))
                .is_some(),
            "legacy alias should stay available in the raw compatibility layer"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "generated_backend"
            )),
            "banded_lambdify"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "matrix_backend"
            )),
            "Banded"
        );
    }

    #[test]
    fn untrusted_aot_compiler_is_blocked_before_solver_execution() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        configure_aot(&mut app, "calc.exe");

        app.request_calculation();

        assert!(app.last_run_is_error);
        assert!(app.pending_aot_confirmation.is_none());
        let message = app
            .last_run_message
            .as_deref()
            .expect("blocked AOT request should produce a visible message");
        assert!(message.contains("not allowed by the GUI"));
        assert!(message.contains("calc.exe"));
    }

    #[test]
    fn unknown_solver_backend_is_rejected_before_worker_creation() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document
            .get_mut("solver_settings")
            .expect("solver_settings should exist")
            .insert(
                "generated_backend".to_string(),
                Some(vec![Value::String("mystery_backend".to_string())]),
            );

        app.request_calculation();

        assert_eq!(app.calculation_state, CalculationState::Failed);
        assert!(app.pending_aot_confirmation.is_none());
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("solver settings"))
        );
    }

    #[test]
    fn incompatible_aot_codegen_and_compiler_are_rejected_before_confirmation() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        configure_aot(&mut app, "tcc");
        app.document
            .get_mut("solver_settings")
            .expect("solver_settings should exist")
            .insert(
                "aot_codegen_backend".to_string(),
                Some(vec![Value::String("Zig".to_string())]),
            );

        app.request_calculation();

        assert_eq!(app.calculation_state, CalculationState::Failed);
        assert!(app.pending_aot_confirmation.is_none());
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("cannot be used with the Zig"))
        );
    }

    #[test]
    fn aot_run_requires_confirmation_and_can_be_cancelled() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        configure_aot(&mut app.borrow_mut(), "tcc");
        app.borrow_mut().request_calculation();
        assert!(app.borrow().pending_aot_confirmation.is_some());

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Confirm AOT toolchain execution");
        harness.get_by_label("Cancel AOT run").click();
        harness.run();

        let app_ref = app.borrow();
        assert!(app_ref.pending_aot_confirmation.is_none());
        assert!(!app_ref.last_run_is_error);
        assert_eq!(
            app_ref.last_run_message.as_deref(),
            Some("AOT calculation cancelled before execution.")
        );
    }

    #[test]
    fn confirmed_aot_run_continues_into_normal_task_validation() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            configure_aot(&mut app_ref, "tcc");
            // Keep this story deterministic: approval should cross the trust
            // boundary, then fail before any compiler is invoked.
            app_ref.document.remove("process_conditions");
            app_ref.request_calculation();
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Run AOT calculation").click();
        // Approval starts a repainting worker state. Advance one frame here;
        // the polling helper below owns completion waiting.
        harness.step();

        wait_for_calculation(&mut app.borrow_mut());

        let app_ref = app.borrow();
        assert!(app_ref.pending_aot_confirmation.is_none());
        assert!(app_ref.last_run_is_error);
        assert!(
            app_ref
                .last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("Calculation failed"))
        );
    }

    #[test]
    fn lambdify_run_does_not_request_aot_confirmation() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document.remove("process_conditions");

        app.request_calculation();

        assert!(app.pending_aot_confirmation.is_none());
        assert_eq!(app.calculation_state, CalculationState::Failed);
        assert!(app.last_run_is_error);
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("task validation failed"))
        );
    }

    #[test]
    fn invalid_task_shows_validation_report_and_disables_run_button() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        app.borrow_mut().document.remove("process_conditions");

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(1400.0, 1600.0))
            .build_ui(move |ui| {
                app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
            });
        harness.run();
        harness.get_by_label("Task validation");
        let run_button = harness.get_by_role_and_label(Role::Button, "🚀 RUN CALCULATION!");
        run_button.scroll_to_me();
        run_button.click();
        harness.step();

        let app_ref = app.borrow();
        assert_eq!(app_ref.calculation_state, CalculationState::Idle);
        assert!(!app_ref.validation_report.is_valid());
    }

    #[test]
    fn cancelling_required_field_deletion_preserves_document_data() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        app.borrow_mut().pending_document_deletion = Some(DocumentDeleteRequest::Field {
            section: "process_conditions".to_string(),
            field: "Tm".to_string(),
        });

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Confirm BVP data deletion");
        harness.get_by_label("Cancel deletion").click();
        harness.run();

        let app_ref = app.borrow();
        assert!(app_ref.pending_document_deletion.is_none());
        assert!(
            app_ref
                .document
                .get("process_conditions")
                .is_some_and(|section| section.contains_key("Tm"))
        );
        assert_eq!(
            app_ref.last_run_message.as_deref(),
            Some("Deletion cancelled; document data was preserved.")
        );
    }

    #[test]
    fn confirming_required_section_deletion_removes_its_payload() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        app.borrow_mut().pending_document_deletion = Some(DocumentDeleteRequest::Section {
            section: "reactions".to_string(),
        });

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Confirm deletion").click();
        harness.run();

        let app_ref = app.borrow();
        assert!(app_ref.pending_document_deletion.is_none());
        assert!(
            app_ref
                .document
                .get("reactions")
                .is_none_or(|section| section.is_empty())
        );
    }

    #[test]
    fn optional_raw_field_deletion_does_not_create_pending_confirmation() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document.insert(
            "custom_metadata".to_string(),
            [(
                "note".to_string(),
                Some(vec![Value::String("keep me".to_string())]),
            )]
            .into_iter()
            .collect(),
        );
        let request = DocumentDeleteRequest::Field {
            section: "custom_metadata".to_string(),
            field: "note".to_string(),
        };

        assert!(apply_document_deletion(&mut app.document, &request));
        assert!(app.pending_document_deletion.is_none());
        assert!(
            app.document
                .get("custom_metadata")
                .is_some_and(|section| !section.contains_key("note"))
        );
    }

    #[test]
    fn backend_comboboxes_drive_lambdify_c_zig_rust_roundtrip() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();

        harness
            .get_by_role_and_label(Role::ComboBox, "Residuals and Jacobian Backend")
            .click();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "AOT").click();
        harness.run();
        {
            let app_ref = app.borrow();
            assert_eq!(
                as_string(first_value(
                    &app_ref.document,
                    "solver_settings",
                    "generated_backend"
                )),
                "banded_aot_tcc"
            );
        }

        harness
            .get_by_role_and_label(Role::ComboBox, "C compiler")
            .click();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "gcc").click();
        harness.run();
        assert_eq!(
            as_string(first_value(
                &app.borrow().document,
                "solver_settings",
                "generated_backend"
            )),
            "banded_aot_gcc"
        );

        harness
            .get_by_role_and_label(Role::ComboBox, "AOT codegen backend")
            .click();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "Zig").click();
        harness.run();
        {
            let app_ref = app.borrow();
            assert_eq!(
                as_string(first_value(
                    &app_ref.document,
                    "solver_settings",
                    "generated_backend"
                )),
                "banded_aot_zig"
            );
            assert!(
                !app_ref
                    .document
                    .get("solver_settings")
                    .is_some_and(|section| section.contains_key("aot_c_compiler"))
            );
        }

        harness
            .get_by_role_and_label(Role::ComboBox, "AOT codegen backend")
            .click();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "Rust").click();
        harness.run();
        {
            let app_ref = app.borrow();
            assert_eq!(
                as_string(first_value(
                    &app_ref.document,
                    "solver_settings",
                    "generated_backend"
                )),
                "banded_aot"
            );
            assert_eq!(
                as_string(first_value(
                    &app_ref.document,
                    "solver_settings",
                    "aot_codegen_backend"
                )),
                "Rust"
            );
        }

        harness
            .get_by_role_and_label(Role::ComboBox, "Residuals and Jacobian Backend")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Lambdify")
            .click();
        harness.run();
        let app_ref = app.borrow();
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "generated_backend"
            )),
            "banded_lambdify"
        );
        assert!(
            !app_ref
                .document
                .get("solver_settings")
                .is_some_and(|section| section.contains_key("aot_codegen_backend"))
        );
    }

    #[test]
    fn worker_lifecycle_blocks_duplicate_runs_and_completes() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.start_test_calculation_worker(Duration::from_millis(20), Ok(()));
        let first_run_id = match app.calculation_state {
            CalculationState::Running { run_id } => run_id,
            other => panic!("expected Running state, got {other:?}"),
        };

        app.start_test_calculation_worker(Duration::from_millis(1), Ok(()));
        assert_eq!(
            app.calculation_state,
            CalculationState::Running {
                run_id: first_run_id
            }
        );
        assert_eq!(
            app.last_run_message.as_deref(),
            Some("A calculation is already running.")
        );

        wait_for_calculation(&mut app);
        assert_eq!(app.calculation_state, CalculationState::Completed);
        assert_eq!(
            app.last_run_message.as_deref(),
            Some("Calculation completed successfully.")
        );
        assert!(!app.last_run_is_error);
    }

    #[test]
    fn worker_lifecycle_cancellation_discards_pending_result() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        app.borrow_mut()
            .start_test_calculation_worker(Duration::from_secs(1), Ok(()));

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(1400.0, 1600.0))
            .build_ui(move |ui| {
                app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
            });
        harness.step();
        harness.get_by_label("Running calculation...");
        harness
            .get_by_role_and_label(Role::Button, "Cancel calculation")
            .scroll_to_me();
        harness.step();
        harness
            .get_by_role_and_label(Role::Button, "Cancel calculation")
            .click();
        harness.step();
        assert!(
            matches!(
                app.borrow().calculation_state,
                CalculationState::Cancelling { .. }
            ),
            "unexpected state after cancel click: {:?}",
            app.borrow().calculation_state
        );

        wait_for_calculation(&mut app.borrow_mut());

        let app_ref = app.borrow();
        assert_eq!(app_ref.calculation_state, CalculationState::Idle);
        assert_eq!(
            app_ref.last_run_message.as_deref(),
            Some("Calculation cancelled.")
        );
        assert!(app_ref.plot_window.is_none());
    }

    #[test]
    fn worker_lifecycle_surfaces_background_failure() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.start_test_calculation_worker(
            Duration::from_millis(5),
            Err("synthetic worker failure".to_string()),
        );

        wait_for_calculation(&mut app);

        assert_eq!(app.calculation_state, CalculationState::Failed);
        assert_eq!(
            app.last_run_message.as_deref(),
            Some("synthetic worker failure")
        );
        assert!(app.last_run_is_error);
    }

    #[test]
    fn successful_bvp_run_hands_plot_window_back_to_the_gui() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        load_problem_test_document(&mut app);
        set_gui_plot(&mut app, true);

        app.request_calculation();
        wait_for_calculation_with_limit(&mut app, 2000, Duration::from_millis(5));

        assert_eq!(app.calculation_state, CalculationState::Completed);
        assert_eq!(
            app.last_run_message.as_deref(),
            Some("Calculation completed successfully.")
        );
        assert!(!app.last_run_is_error);

        let plot_window = app
            .plot_window
            .as_ref()
            .expect("successful BVP run should create a plot window when gui_plot is enabled");
        assert!(plot_window.visible);
    }

    #[test]
    fn successful_bvp_run_skips_plot_window_when_gui_plot_is_disabled() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        load_problem_test_document(&mut app);
        set_gui_plot(&mut app, false);

        app.request_calculation();
        wait_for_calculation_with_limit(&mut app, 2000, Duration::from_millis(5));

        assert_eq!(app.calculation_state, CalculationState::Completed);
        assert_eq!(
            app.last_run_message.as_deref(),
            Some("Calculation completed successfully.")
        );
        assert!(!app.last_run_is_error);
        assert!(app.plot_window.is_none());
    }
}
