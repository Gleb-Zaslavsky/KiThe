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
//! - A checked-in physical BVP document should survive GUI edit/save/reopen and run.
//!
//! Expected result:
//! the BVP screen should behave like a structured form with a compatibility
//! fallback, not like a raw document editor with a few extra labels on top.

#[cfg(test)]
mod tests {
    use super::super::combustion::{
        BvpEditorTab, CalculationState, CombustionApp, DocumentDeleteRequest,
        GuiFileOperationResult, ProblemsEnum, apply_document_deletion,
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

    fn assert_control_reachable(harness: &mut Harness<'_>, label: &str, screen_rect: egui::Rect) {
        harness.get_by_label(label).scroll_to_me();
        harness.run();
        harness.run();
        let rect = harness.get_by_label(label).rect();
        assert!(
            rect.intersects(screen_rect),
            "control `{label}` did not scroll into the compact viewport: {rect:?} outside {screen_rect:?}"
        );
    }

    fn set_solver_setting(app: &Rc<RefCell<CombustionApp>>, field: &str, value: &str) {
        let mut app = app.borrow_mut();
        let section = app
            .document
            .get_mut("solver_settings")
            .expect("solver_settings section must exist in the BVP template");
        section.insert(
            field.to_string(),
            Some(vec![Value::String(value.to_string())]),
        );
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

    fn set_return_to_dimension(app: &mut CombustionApp, enabled: bool) {
        let postprocessing = app
            .document
            .get_mut("postprocessing")
            .expect("postprocessing section must exist in the BVP template");
        postprocessing.insert(
            "return_to_dimension".to_string(),
            Some(vec![Value::Boolean(enabled)]),
        );
    }

    /// Load the known BVP story document used by the GUI regression suite.
    fn load_problem_test_document(app: &mut CombustionApp) {
        let result =
            app.load_document_from_path(PathBuf::from("src/gui/fixtures/bvp_problem.fixture"));
        assert!(
            matches!(result, GuiFileOperationResult::Loaded { .. }),
            "tracked BVP GUI fixture should load as a valid document"
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
    fn physical_bvp_task_doc_gui_load_edit_save_reload_and_run() {
        let source = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("task_docs/reactors/bvp/hmx_compact_two_percent_product.task");
        let source_before = std::fs::read_to_string(&source).expect("source task should exist");
        let mut editor = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        assert!(matches!(
            editor.load_document_from_path(source.clone()),
            GuiFileOperationResult::Loaded { .. }
        ));
        assert_eq!(
            as_string(first_value(
                &editor.document,
                "process_conditions",
                "problem_name"
            )),
            "Compact HMX decomposition - alternate inlet"
        );

        let editor = Rc::new(RefCell::new(editor));
        let mut open = true;
        let editor_for_ui = Rc::clone(&editor);
        let mut harness = Harness::new_ui(move |ui| {
            editor_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("A").scroll_to_me();
        harness.run();
        harness.get_by_label("A").click();
        harness.run();
        harness.key_combination_modifiers(egui::Modifiers::CTRL, &[egui::Key::A]);
        harness.get_by_label("A").type_text("140000");
        harness.key_press(egui::Key::Enter);
        harness.run();
        assert_eq!(
            first_value(
                &editor.borrow().document,
                "reactions",
                "HMX=>DecompositionProduct"
            )
            .as_vector()
            .expect("reaction parameters should remain a vector")[0],
            140000.0
        );
        set_gui_plot(&mut editor.borrow_mut(), true);
        set_return_to_dimension(&mut editor.borrow_mut(), true);
        editor.borrow_mut().sync_document_state_after_edit();
        let gui_config_before_save = editor.borrow().bvp_gui_config.clone();

        let tempdir = tempfile::tempdir().expect("temporary directory should be created");
        let edited_path = tempdir.path().join("hmx_edited.task");
        assert!(matches!(
            editor.borrow_mut().save_document(edited_path.clone()),
            GuiFileOperationResult::Saved { .. }
        ));
        assert_eq!(
            std::fs::read_to_string(&source).expect("source task should remain readable"),
            source_before,
            "Save As must not overwrite the checked-in source task"
        );

        let mut reopened = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        assert!(matches!(
            reopened.load_document_from_path(edited_path.clone()),
            GuiFileOperationResult::Loaded { .. }
        ));
        assert!(matches!(
            first_value(&reopened.document, "process_conditions", "Tm"),
            Value::Float(1500.0)
        ));
        assert!(matches!(
            first_value(&reopened.document, "process_conditions", "t0"),
            Value::Float(0.0)
        ));
        assert!(matches!(
            first_value(&reopened.document, "bounds", "Teta"),
            Value::Float(-100.0)
        ));
        let reopened_species = reopened.document["process_conditions"]["substances"]
            .as_ref()
            .expect("species field should remain populated")
            .iter()
            .map(as_string)
            .collect::<Vec<_>>();
        assert_eq!(reopened_species, ["HMX", "DecompositionProduct"]);
        for (field, expected) in [
            ("L", 9.0e-4),
            ("dT", 600.0),
            ("T_scale", 600.0),
            ("P", 1.0e6),
            ("Cp", 1464.4),
            ("Lambda", 0.07),
            ("m", 0.00431),
            ("M", 0.0342),
        ] {
            assert_eq!(
                first_value(&reopened.document, "process_conditions", field),
                &Value::Float(expected),
                "physical task field {field} changed during Save As/reopen"
            );
        }
        assert_eq!(
            first_value(&reopened.document, "process_conditions", "thermal_effects")
                .as_vector()
                .expect("thermal effects should remain a vector"),
            &[30000.0]
        );
        for (field, expected) in [("HMX", 0.98), ("DecompositionProduct", 0.02), ("T", 600.0)] {
            assert_eq!(
                first_value(&reopened.document, "boundary_condition", field),
                &Value::Float(expected),
                "boundary input {field} changed during Save As/reopen"
            );
        }
        for species in ["HMX", "DecompositionProduct"] {
            assert_eq!(
                first_value(&reopened.document, "diffusion_coefficients", species),
                &Value::Float(6.97e-6),
                "diffusion coefficient for {species} changed during Save As/reopen"
            );
        }
        assert_eq!(reopened.bvp_gui_config, gui_config_before_save);
        assert_eq!(
            first_value(&reopened.document, "reactions", "HMX=>DecompositionProduct")
                .as_vector()
                .expect("reaction parameters should remain a vector")[0],
            140000.0
        );

        let mut direct = crate::SimpleReactorsBVP::SimpleReactorBVP::SimpleReactorTask::new();
        direct
            .solve_from_file(edited_path.clone())
            .expect("the saved BVP task should solve through the direct file route");
        let direct_mesh = direct
            .solver
            .x_mesh
            .as_ref()
            .expect("direct file solve should publish its converted mesh");
        let direct_solution = direct
            .solver
            .solution
            .as_ref()
            .expect("direct file solve should publish its converted solution");
        let length = direct.characteristic_length();
        assert!(direct_mesh[0].abs() <= 1e-12);
        assert!(
            (direct_mesh[direct_mesh.len() - 1] - length).abs() <= 1e-12,
            "published endpoint {} should match physical length {}",
            direct_mesh[direct_mesh.len() - 1],
            length
        );
        let temperature_index = direct
            .solver
            .unknowns
            .iter()
            .position(|name| name == "Teta")
            .expect("temperature state should exist");
        for temperature in direct_solution.column(temperature_index).iter() {
            assert!(
                (0.0..3000.0).contains(temperature),
                "published temperature should be in kelvin, got {temperature}"
            );
        }

        reopened.request_calculation();
        wait_for_calculation(&mut reopened);
        assert_eq!(
            reopened.calculation_state,
            CalculationState::Completed,
            "reopened physical BVP task failed: {:?}",
            reopened.last_run_message
        );
        assert!(
            !reopened.last_run_is_error,
            "reopened physical BVP task failed: {:?}",
            reopened.last_run_message
        );
        let plot = reopened
            .plot_window
            .as_ref()
            .expect("gui_plot=true should hand the completed solution to the plot window");
        assert_eq!(plot.t_result.len(), direct_mesh.len());
        assert_eq!(plot.y_result.shape(), direct_solution.shape());
        assert_eq!(plot.values.len(), direct.solver.unknowns.len());
        for (row, (&gui_x, &direct_x)) in plot.t_result.iter().zip(direct_mesh.iter()).enumerate() {
            assert!(
                (gui_x - direct_x).abs() <= 1e-12,
                "GUI/direct mesh differs at row {row}: {gui_x} vs {direct_x}"
            );
            for column in 0..plot.y_result.ncols() {
                let gui_value = plot.y_result[(row, column)];
                let direct_value = direct_solution[(row, column)];
                let tolerance = 1e-8 * direct_value.abs().max(1.0);
                assert!(
                    (gui_value - direct_value).abs() <= tolerance,
                    "GUI/direct solution differs at ({row}, {column}): {gui_value} vs {direct_value}"
                );
            }
        }
        let expected_plot_values = direct
            .solver
            .unknowns
            .iter()
            .map(|name| match name.as_str() {
                "Teta" => "Temperature (K)".to_string(),
                "q" => "q (W/m^2)".to_string(),
                name if name.starts_with('J') => format!("{name} (kg/(m^2 s))"),
                name if name.starts_with('C') => format!("{name} (mass fraction)"),
                _ => name.clone(),
            })
            .collect::<Vec<_>>();
        assert_eq!(plot.values, expected_plot_values);

        set_return_to_dimension(&mut reopened, false);
        reopened.sync_document_state_after_edit();
        reopened.request_calculation();
        wait_for_calculation(&mut reopened);
        assert_eq!(reopened.calculation_state, CalculationState::Completed);
        let dimensionless_plot = reopened
            .plot_window
            .as_ref()
            .expect("dimensionless result should also be plotted");
        assert!(dimensionless_plot.arg.contains("/ L"));
        assert!(dimensionless_plot.t_result[0].abs() <= 1e-12);
        assert!(
            (dimensionless_plot.t_result[dimensionless_plot.t_result.len() - 1] - 1.0).abs()
                <= 1e-12
        );
        assert!(
            dimensionless_plot
                .values
                .iter()
                .any(|value| value == "Teta (dimensionless)")
        );
        assert!(
            dimensionless_plot
                .values
                .iter()
                .any(|value| value == "q_hat (scaled)")
        );
        assert!(
            dimensionless_plot
                .values
                .iter()
                .any(|value| value == "J_hat0 (scaled)")
        );
    }

    #[test]
    fn malformed_bvp_replacement_preserves_dirty_active_document_transactionally() {
        let tempdir = tempfile::tempdir().expect("temporary directory should be created");
        let active_path = tempdir.path().join("active.task");
        let fixture_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("task_docs/reactors/bvp/hmx_compact.task");
        let mut template = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        assert!(matches!(
            template.load_document_from_path(fixture_path),
            GuiFileOperationResult::Loaded { .. }
        ));
        set_gui_plot(&mut template, true);
        template.sync_document_state_after_edit();
        std::fs::write(&active_path, template.document_to_string())
            .expect("valid task should be written");

        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        assert!(matches!(
            app.load_document_from_path(active_path.clone()),
            GuiFileOperationResult::Loaded { .. }
        ));
        app.request_calculation();
        wait_for_calculation(&mut app);
        assert_eq!(app.calculation_state, CalculationState::Completed);
        let successful_plot = app
            .plot_window
            .as_ref()
            .expect("active physical task should publish its successful plot");
        let successful_plot_data = (
            successful_plot.arg.clone(),
            successful_plot.values.clone(),
            successful_plot.t_result.clone(),
            successful_plot.y_result.clone(),
        );
        app.document
            .entry("user_notes".to_string())
            .or_default()
            .insert(
                "draft".to_string(),
                Some(vec![Value::String("unsaved work".to_string())]),
            );
        app.sync_document_state_after_edit();
        assert!(app.document_lifecycle.dirty);

        let document_before = app.document_to_string();
        let gui_config_before = app.bvp_gui_config.clone();
        let path_before = app.current_file_path.clone();
        let malformed_path = tempdir.path().join("malformed.task");
        std::fs::write(
            &malformed_path,
            "process_conditions\nmalformed line without colon\n",
        )
        .expect("malformed task should be written");

        assert!(!app.request_document_load(malformed_path.clone()));
        assert!(!app.confirm_pending_document_load());
        assert_eq!(app.document_to_string(), document_before);
        assert_eq!(app.bvp_gui_config, gui_config_before);
        assert_eq!(app.current_file_path, path_before);
        assert!(app.document_lifecycle.dirty);
        let preserved_plot = app
            .plot_window
            .as_ref()
            .expect("failed replacement must preserve the prior successful result");
        assert_eq!(preserved_plot.arg, successful_plot_data.0);
        assert_eq!(preserved_plot.values, successful_plot_data.1);
        assert_eq!(preserved_plot.t_result, successful_plot_data.2);
        assert_eq!(preserved_plot.y_result, successful_plot_data.3);
        match app.document_lifecycle.last_file_operation.as_ref() {
            Some(GuiFileOperationResult::Failed {
                kind: super::super::combustion::GuiFileOperationKind::Load,
                message,
                ..
            }) => assert!(
                !message.is_empty(),
                "load failure should explain the parse error"
            ),
            other => panic!("expected a typed load failure, got {other:?}"),
        }
    }

    #[test]
    fn rejected_bvp_rerun_clears_previous_plot_result() {
        let fixture_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("task_docs/reactors/bvp/hmx_compact.task");
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        assert!(matches!(
            app.load_document_from_path(fixture_path),
            GuiFileOperationResult::Loaded { .. }
        ));
        set_gui_plot(&mut app, true);
        app.sync_document_state_after_edit();
        app.request_calculation();
        wait_for_calculation(&mut app);
        assert_eq!(app.calculation_state, CalculationState::Completed);
        assert!(app.plot_window.is_some());

        app.document
            .get_mut("solver_settings")
            .expect("solver settings should exist")
            .insert(
                "generated_backend".to_string(),
                Some(vec![Value::String("not_a_backend".to_string())]),
            );
        app.sync_document_state_after_edit();
        app.request_calculation();

        assert_eq!(app.calculation_state, CalculationState::Failed);
        assert!(app.plot_window.is_none());
        assert!(app.last_run_is_error);
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("solver settings"))
        );
    }

    #[test]
    fn bvp_numerical_failure_can_be_repaired_and_rerun() {
        let fixture_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("task_docs/reactors/bvp/hmx_compact.task");
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        assert!(matches!(
            app.load_document_from_path(fixture_path),
            GuiFileOperationResult::Loaded { .. }
        ));
        set_gui_plot(&mut app, true);
        set_return_to_dimension(&mut app, true);
        app.document
            .get_mut("solver_settings")
            .expect("solver settings should exist")
            .insert("max_iterations".to_string(), Some(vec![Value::Integer(1)]));
        app.sync_document_state_after_edit();

        app.request_calculation();
        wait_for_calculation(&mut app);
        assert_eq!(
            app.calculation_state,
            CalculationState::Failed,
            "one nonlinear iteration should not satisfy the HMX BVP: {:?}",
            app.last_run_message
        );
        assert!(app.plot_window.is_none());

        app.document
            .get_mut("solver_settings")
            .expect("solver settings should exist")
            .insert(
                "max_iterations".to_string(),
                Some(vec![Value::Integer(100)]),
            );
        app.sync_document_state_after_edit();
        app.request_calculation();
        wait_for_calculation(&mut app);
        assert_eq!(
            app.calculation_state,
            CalculationState::Completed,
            "corrected BVP task did not recover: {:?}",
            app.last_run_message
        );
        assert!(!app.last_run_is_error);
        let plot = app
            .plot_window
            .as_ref()
            .expect("successful recovery should publish a new plot");
        assert!(!plot.t_result.is_empty());
        assert!(plot.y_result.iter().all(|value| value.is_finite()));
        let temperature_index = plot
            .values
            .iter()
            .position(|name| name == "Temperature (K)")
            .unwrap_or_else(|| {
                panic!(
                    "recovered result should include physical temperature; columns={:?}",
                    plot.values
                )
            });
        assert!(
            plot.y_result
                .column(temperature_index)
                .iter()
                .all(|temperature| (200.0..5000.0).contains(temperature))
        );
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

        harness.get_by_role_and_label(Role::Button, "Physics");
        for label in ["Legacy raw document"] {
            harness.get_by_label(label);
        }
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
        harness.run();
        for label in [
            "Solver backend",
            "Initial guess",
            "Postprocessing",
            "Equation assembly",
            "Assembly parallelism",
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
        app.borrow_mut().active_bvp_tab = BvpEditorTab::Solver;
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
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
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
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
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
            solver_settings.insert(
                "legacy_custom_field".to_string(),
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
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
        harness.run();
        harness.get_by_label("Solver backend");
        harness
            .get_by_role_and_label(Role::Button, "Physics")
            .click();
        harness.run();
        harness.get_by_label("Legacy raw document");

        let app_ref = app.borrow();
        assert!(
            app_ref
                .document
                .get("solver_settings")
                .and_then(|section| section.get("legacy_custom_field"))
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
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(820.0, 620.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
        harness.run();

        harness
            .get_by_role_and_label(Role::ComboBox, "Residuals and Jacobian Backend")
            .scroll_to_me();
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
            .scroll_to_me();
        harness.run();
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
    fn bvp_family_switch_exposes_capability_specific_controls_in_compact_window() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(820.0, 620.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
        harness.run();
        harness.get_by_label("Compact HMX physical acceptance: passed for Damped.");
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "sci")
            .click_accesskit();
        harness.run();

        assert_eq!(
            as_string(first_value(
                &app.borrow().document,
                "solver_settings",
                "solver_family"
            )),
            "sci"
        );
        harness.get_by_label("Compact HMX physical acceptance: passed for BVP_sci.");
        for label in ["Execution", "Symbolic assembly", "Matrix layout"] {
            harness.get_by_label(label);
        }

        harness
            .get_by_role_and_label(Role::ComboBox, "Execution")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Execution")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "aot")
            .click_accesskit();
        harness.run();
        harness.get_by_label("AOT codegen backend");
        harness.get_by_label("C compiler");

        harness
            .get_by_role_and_label(Role::ComboBox, "AOT codegen backend")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "AOT codegen backend")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Zig")
            .click_accesskit();
        harness.run();
        assert!(
            !app.borrow()
                .document
                .get("solver_settings")
                .is_some_and(|section| section.contains_key("aot_c_compiler")),
            "non-C AOT codegen must remove the irrelevant C compiler setting"
        );

        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "frozen")
            .click_accesskit();
        harness.run();
        assert_eq!(
            as_string(first_value(
                &app.borrow().document,
                "solver_settings",
                "solver_family"
            )),
            "frozen"
        );
        harness.get_by_label(
            "Compact HMX physical acceptance is open: Frozen/Naive do not currently meet the boundary contract.",
        );
        harness.get_by_label("Core solve");
        harness.get_by_label("Linear algebra backend");
    }

    #[test]
    fn bvp_conditional_solver_controls_are_reachable_in_compact_viewport() {
        let size = egui::Vec2::new(820.0, 620.0);
        let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, size);
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(size)
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
        harness.run();

        for label in [
            "Solver family",
            "scheme",
            "Linear algebra backend",
            "Nonlinear strategy",
            "Linear system method override",
            "abs_tolerance",
            "max_iterations",
            "Solver telemetry",
            "loglevel",
            "Residuals and Jacobian Backend",
            "Symbolic backend",
            "Lambdify callback policy",
            "Banded linear solver",
            "refinement_steps",
            "Assembly parallelism",
        ] {
            assert_control_reachable(&mut harness, label, screen_rect);
        }

        set_solver_setting(&app, "lambdify_execution_policy", "parallel");
        harness.run();
        assert_control_reachable(&mut harness, "Minimum callback work", screen_rect);

        set_solver_setting(&app, "generated_backend", "banded_aot_tcc");
        set_solver_setting(&app, "backend_policy", "aot_only");
        harness.run();
        for label in [
            "AOT codegen backend",
            "C compiler",
            "AOT build policy",
            "AOT build profile",
            "AOT compile preset",
            "AOT execution policy",
            "Residual target chunks",
            "Jacobian target chunks",
        ] {
            assert_control_reachable(&mut harness, label, screen_rect);
        }

        set_solver_setting(&app, "solver_family", "sci");
        set_solver_setting(&app, "execution", "lambdify");
        set_solver_setting(&app, "matrix_layout", "sparse");
        set_solver_setting(&app, "execution_policy", "parallel");
        harness.run();
        for label in [
            "Execution",
            "Symbolic assembly",
            "Matrix layout",
            "Tolerance",
            "Boundary tolerance",
            "Maximum mesh nodes",
            "Maximum Newton iterations",
            "Maximum Jacobian refreshes",
            "Maximum mesh refinements",
            "Maximum backtracking steps",
            "Solver telemetry",
            "Logging level",
            "Solution output",
        ] {
            assert_control_reachable(&mut harness, label, screen_rect);
        }

        assert_control_reachable(&mut harness, "Sparse factorization", screen_rect);
        assert_control_reachable(&mut harness, "Minimum callback work", screen_rect);

        set_solver_setting(&app, "matrix_layout", "banded");
        harness.run();
        for label in ["Lower bandwidth", "Upper bandwidth"] {
            assert_control_reachable(&mut harness, label, screen_rect);
        }

        set_solver_setting(&app, "matrix_layout", "sparse");
        set_solver_setting(&app, "execution", "aot");
        set_solver_setting(&app, "aot_residual_chunking", "output_count");
        set_solver_setting(&app, "aot_sparse_jacobian_chunking", "rows");
        set_solver_setting(&app, "aot_sparse_jacobian_chunk_size", "4");
        harness.run();
        for label in [
            "AOT codegen backend",
            "AOT build policy",
            "AOT build profile",
            "AOT output directory",
            "AOT handoff path",
            "Generated crate name",
            "Generated module name",
            "Residual chunking",
            "Residual chunk size",
            "Sparse Jacobian chunking",
            "Sparse Jacobian chunk size",
        ] {
            assert_control_reachable(&mut harness, label, screen_rect);
        }

        set_solver_setting(&app, "matrix_layout", "banded");
        set_solver_setting(&app, "aot_banded_jacobian_chunking", "diagonals");
        set_solver_setting(&app, "aot_banded_jacobian_chunk_size", "4");
        harness.run();
        for label in ["Banded Jacobian chunking", "Banded Jacobian chunk size"] {
            assert_control_reachable(&mut harness, label, screen_rect);
        }
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
        assert_eq!(
            app.calculation_state,
            CalculationState::Completed,
            "test worker did not complete: {:?}",
            app.last_run_message
        );
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
        set_return_to_dimension(&mut app, true);

        app.request_calculation();
        wait_for_calculation_with_limit(&mut app, 2000, Duration::from_millis(5));

        assert_eq!(
            app.calculation_state,
            CalculationState::Completed,
            "BVP run did not complete: {:?}",
            app.last_run_message
        );
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
        let temperature_column = plot_window
            .values
            .iter()
            .position(|name| name == "Temperature (K)")
            .unwrap_or_else(|| {
                panic!(
                    "the GUI plot should label the dimensional temperature in kelvins; columns={:?}",
                    plot_window.values
                )
            });
        let temperatures = plot_window.y_result.column(temperature_column);
        assert!(
            temperatures.iter().all(|temperature| temperature.is_finite()
                && (200.0..5000.0).contains(temperature)),
            "plotted temperature must be a physical Kelvin profile, got {temperatures:?}"
        );
    }

    #[test]
    fn successful_bvp_run_skips_plot_window_when_gui_plot_is_disabled() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        load_problem_test_document(&mut app);
        set_gui_plot(&mut app, false);

        app.request_calculation();
        wait_for_calculation_with_limit(&mut app, 2000, Duration::from_millis(5));

        assert_eq!(
            app.calculation_state,
            CalculationState::Completed,
            "BVP run did not complete: {:?}",
            app.last_run_message
        );
        assert_eq!(
            app.last_run_message.as_deref(),
            Some("Calculation completed successfully.")
        );
        assert!(!app.last_run_is_error);
        assert!(app.plot_window.is_none());
    }
}
