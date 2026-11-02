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

    fn assert_ivp_config_equal(expected: &IvpGuiConfig, actual: &IvpGuiConfig) {
        macro_rules! assert_fields {
            ($($field:ident),+ $(,)?) => {
                $(
                    assert_eq!(
                        expected.$field,
                        actual.$field,
                        "IVP config field `{}` changed during document roundtrip",
                        stringify!($field)
                    );
                )+
            };
        }

        assert_fields!(
            problem_name,
            problem_description,
            physical,
            initial_state,
            substances,
            thermal_effects,
            stop_condition_enabled,
            stop_condition_species_index,
            stop_condition_threshold,
            solver_family,
            method,
            symbolic_backend,
            matrix_backend,
            execution_mode,
            callback_policy,
            callback_min_work,
            aot_toolchain,
            aot_profile,
            bdf_max_steps,
            radau_max_steps,
            radau_max_newton_iterations,
            x0,
            x_bound,
            first_step,
            max_step,
            rtol,
            atol,
            telemetry,
            logging,
            postprocessing,
            native_execution,
        );
        assert_eq!(expected.scaling.dT, actual.scaling.dT);
        assert_eq!(expected.scaling.L, actual.scaling.L);
        assert_eq!(expected.scaling.T_scale, actual.scaling.T_scale);
        assert_eq!(
            serde_json::to_value(&expected.reactions).unwrap(),
            serde_json::to_value(&actual.reactions).unwrap(),
            "IVP reaction payload changed during document roundtrip"
        );
    }

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
        config.initial_state.species.insert("H2".to_string(), 0.78);
        config.initial_state.species.insert("H".to_string(), 0.22);
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
        assert_eq!(rebuilt.initial_state.species["H2"], 0.78);
        assert_eq!(rebuilt.initial_state.species["H"], 0.22);
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
            .insert("H2".to_string(), 0.74);
        app.config
            .initial_state
            .species
            .insert("H".to_string(), 0.26);
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
        assert_eq!(loaded.config.initial_state.species["H2"], 0.74);
        assert_eq!(loaded.config.initial_state.species["H"], 0.26);
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
    fn physical_task_doc_gui_edit_save_reload_preview_and_solve() {
        let source_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("task_docs/reactors/ivp/h2_dissociation_with_heat_release.task");
        let source_before =
            std::fs::read_to_string(&source_path).expect("source task should exist");
        let mut editor = ReactorIvpApp::new();
        assert!(matches!(
            editor.load_document_from_path(source_path.clone()),
            crate::gui::document_lifecycle::GuiFileOperationResult::Loaded { .. }
        ));
        assert_eq!(editor.config.substances, vec!["H2", "H"]);
        assert_eq!(editor.config.reactions[0].eq, "H2=>2H");

        let editor = Rc::new(RefCell::new(editor));
        let mut open = true;
        let editor_for_ui = Rc::clone(&editor);
        let mut edit_harness = Harness::new_ui(move |ui| {
            editor_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        edit_harness.run();
        edit_harness.get_by_label("Temperature").scroll_to_me();
        edit_harness.run();
        edit_harness.get_by_label("Temperature").click();
        edit_harness.run();
        edit_harness.key_combination_modifiers(egui::Modifiers::CTRL, &[egui::Key::A]);
        edit_harness.get_by_label("Temperature").type_text("530");
        edit_harness.key_press(egui::Key::Enter);
        edit_harness.run();
        assert_eq!(editor.borrow().config.initial_state.temperature, 530.0);
        assert!(editor.borrow().document_lifecycle.dirty);
        let config_before_save = editor.borrow().config.clone();

        let tempdir = tempdir().expect("temporary directory should be created");
        let edited_path = tempdir.path().join("h2_dissociation_edited.task");
        let failed_save_path = tempdir.path().join("missing-parent").join("failed.task");
        let dirty_fingerprint = editor.borrow().current_document_fingerprint().unwrap();
        assert!(matches!(
            editor.borrow_mut().save_document(failed_save_path),
            crate::gui::document_lifecycle::GuiFileOperationResult::Failed { .. }
        ));
        assert_eq!(
            editor.borrow().current_file_path.as_ref(),
            Some(&source_path)
        );
        assert!(editor.borrow().document_lifecycle.dirty);
        assert_eq!(
            editor.borrow().current_document_fingerprint().unwrap(),
            dirty_fingerprint,
            "failed Save As must not mutate the active task"
        );
        assert!(matches!(
            editor.borrow_mut().save_document(edited_path.clone()),
            crate::gui::document_lifecycle::GuiFileOperationResult::Saved { .. }
        ));
        assert_eq!(
            std::fs::read_to_string(&source_path).expect("source task should remain readable"),
            source_before,
            "Save As must not overwrite the checked-in source task"
        );

        let mut reloaded = ReactorIvpApp::new();
        assert!(matches!(
            reloaded.load_document_from_path(edited_path.clone()),
            crate::gui::document_lifecycle::GuiFileOperationResult::Loaded { .. }
        ));
        assert_eq!(reloaded.config.initial_state.temperature, 530.0);
        assert_eq!(reloaded.config.substances, vec!["H2", "H"]);
        assert_eq!(reloaded.config.reactions[0].eq, "H2=>2H");
        assert_ivp_config_equal(&config_before_save, &reloaded.config);
        let mut direct_task = reloaded.task().clone();
        reloaded
            .config
            .apply_to_task(&mut direct_task)
            .expect("reopened document should apply to the direct production task");
        let direct_snapshot = direct_task
            .solve()
            .expect("reopened document should solve through the direct production API");
        let direct_balance = direct_task
            .latest_conservation_report()
            .expect("direct API solve should publish conservation diagnostics")
            .clone();

        let invalid_path = tempdir.path().join("invalid.task");
        let edited_document = std::fs::read_to_string(&edited_path)
            .expect("saved edited task should remain readable");
        let malformed_mechanism = edited_document.replace(
            "reaction_000_arrhenius: [0.01, 0, 8000]",
            "reaction_000_arrhenius: [1e-3, 0.0]",
        );
        assert_ne!(malformed_mechanism, edited_document);
        std::fs::write(&invalid_path, malformed_mechanism)
            .expect("malformed mechanism task should be written");
        let before_failed_load = reloaded.current_document_fingerprint().unwrap();
        let load_result = reloaded.load_document_from_path(invalid_path.clone());
        match load_result {
            crate::gui::document_lifecycle::GuiFileOperationResult::Failed { message, .. } => {
                assert!(
                    message.contains("kinetics reaction 000"),
                    "mechanism parse error should identify the reaction: {message}"
                );
            }
            other => panic!("malformed mechanism should be rejected, got {other:?}"),
        }
        assert_eq!(
            reloaded.current_document_fingerprint().unwrap(),
            before_failed_load,
            "a rejected file must not partially replace the active editor state"
        );

        let repaired_path = tempdir.path().join("repaired.task");
        let repaired_document = std::fs::read_to_string(&invalid_path)
            .expect("malformed task should remain readable")
            .replace(
                "reaction_000_arrhenius: [1e-3, 0.0]",
                "reaction_000_arrhenius: [0.01, 0, 8000]",
            );
        assert_ne!(
            repaired_document,
            std::fs::read_to_string(&invalid_path).unwrap()
        );
        std::fs::write(&repaired_path, repaired_document).expect("repaired task should be written");
        assert!(matches!(
            reloaded.load_document_from_path(repaired_path),
            crate::gui::document_lifecycle::GuiFileOperationResult::Loaded { .. }
        ));
        assert_eq!(reloaded.config.initial_state.temperature, 530.0);
        assert_eq!(reloaded.config.reactions[0].eq, "H2=>2H");
        let recovered_path = tempdir.path().join("recovered.task");
        assert!(matches!(
            reloaded.save_document(recovered_path.clone()),
            crate::gui::document_lifecycle::GuiFileOperationResult::Saved { .. }
        ));
        let mut reloaded = ReactorIvpApp::new();
        assert!(matches!(
            reloaded.load_document_from_path(recovered_path),
            crate::gui::document_lifecycle::GuiFileOperationResult::Loaded { .. }
        ));
        assert_eq!(reloaded.config.initial_state.temperature, 530.0);

        let app = Rc::new(RefCell::new(reloaded));
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
        assert!(app.borrow().last_preview().is_some());
        harness
            .get_by_role_and_label(Role::Button, "Run Calculation")
            .click_accesskit();
        for _ in 0..250 {
            harness.step();
            if matches!(
                app.borrow().calculation_state(),
                ReactorIvpCalculationState::Completed | ReactorIvpCalculationState::Failed
            ) {
                break;
            }
            std::thread::sleep(Duration::from_millis(5));
        }
        assert!(
            matches!(
                app.borrow().calculation_state(),
                ReactorIvpCalculationState::Completed
            ),
            "reloaded physical document did not solve: {:?}",
            app.borrow().last_error()
        );
        let app_ref = app.borrow();
        let render = app_ref
            .task()
            .solution_render_data()
            .expect("reloaded task should publish physical solution data");
        assert!(render.solution.iter().all(|value| value.is_finite()));
        assert_eq!(
            render.unknowns.first().map(String::as_str),
            Some("Temperature (K)")
        );
        let gui_final = render.solution.row(render.solution.nrows() - 1);
        let direct_final = direct_snapshot
            .values
            .row(direct_snapshot.values.nrows() - 1);
        let direct_temperature_index = direct_snapshot
            .variable_names
            .iter()
            .position(|name| name == "Teta")
            .expect("direct output should include temperature");
        let direct_temperature = direct_final[direct_temperature_index]
            * direct_task.scaling.T_scale
            + direct_task.scaling.dT;
        assert!(
            (gui_final[0] - direct_temperature).abs() <= 1e-8,
            "GUI/direct final temperatures differ: {} vs {direct_temperature}",
            gui_final[0]
        );
        for species_index in 0..2 {
            let state_name = format!("C{species_index}");
            let gui_index = render
                .unknowns
                .iter()
                .position(|name| name == &state_name)
                .expect("GUI output should include each species fraction");
            let direct_index = direct_snapshot
                .variable_names
                .iter()
                .position(|name| name == &state_name)
                .expect("direct output should include each species");
            assert!(
                (gui_final[gui_index] - direct_final[direct_index]).abs() <= 1e-8,
                "GUI/direct final {state_name} fractions differ"
            );
        }
        let gui_balance = app_ref
            .task()
            .latest_conservation_report()
            .expect("GUI solve should publish conservation diagnostics");
        assert!(
            (gui_balance.energy_balance_error_abs - direct_balance.energy_balance_error_abs).abs()
                <= 1e-8
        );
        assert!(
            (gui_balance.energy_balance_error_rel - direct_balance.energy_balance_error_rel).abs()
                <= 1e-8
        );
        assert_eq!(
            gui_balance.sum_of_mass_fractions,
            direct_balance.sum_of_mass_fractions
        );
        assert_eq!(
            gui_balance.atomic_mass_balance_error,
            direct_balance.atomic_mass_balance_error
        );
    }

    #[test]
    fn malformed_indexed_mechanisms_fail_with_context_without_replacing_active_task() {
        let source_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("task_docs/reactors/ivp/h2_dissociation_with_heat_release.task");
        let source = std::fs::read_to_string(&source_path).expect("fixture should be readable");
        let cases = [
            (
                "bad equation",
                source.replace(
                    "reaction_000_equation: H2=>2H",
                    "reaction_000_equation: H2=>2H+",
                ),
                "reaction",
            ),
            (
                "missing reaction index",
                source
                    .replace("reaction_000_equation", "reaction_001_equation")
                    .replace("reaction_000_arrhenius", "reaction_001_arrhenius"),
                "contiguous",
            ),
            (
                "unknown species",
                source.replace(
                    "reaction_000_equation: H2=>2H",
                    "reaction_000_equation: H2=>2X",
                ),
                "species",
            ),
            (
                "reaction/thermal-effect count mismatch",
                source.replace(
                    "thermal_effects: [-150000.0]",
                    "thermal_effects: [-150000.0, -1.0]",
                ),
                "reaction count",
            ),
        ];

        let mut app = ReactorIvpApp::new();
        assert!(matches!(
            app.load_document_from_path(source_path.clone()),
            crate::gui::document_lifecycle::GuiFileOperationResult::Loaded { .. }
        ));
        let original_fingerprint = app.current_document_fingerprint().unwrap();
        let original_path = app.current_file_path.clone();
        let tempdir = tempdir().expect("temporary directory should be created");

        for (index, (label, malformed, expected_context)) in cases.into_iter().enumerate() {
            assert_ne!(malformed, source, "case {label} must alter the fixture");
            let path = tempdir.path().join(format!("malformed_{index}.task"));
            std::fs::write(&path, malformed).expect("malformed task should be written");
            match app.load_document_from_path(path) {
                crate::gui::document_lifecycle::GuiFileOperationResult::Failed {
                    message, ..
                } => assert!(
                    message.to_lowercase().contains(expected_context),
                    "{label} error lacks expected context `{expected_context}`: {message}"
                ),
                other => panic!("{label} should be rejected without replacing the task: {other:?}"),
            }
            assert_eq!(
                app.current_document_fingerprint().unwrap(),
                original_fingerprint,
                "{label} failure partially replaced the active task"
            );
            assert_eq!(app.current_file_path, original_path, "{label}");
        }
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
    fn failed_ivp_rerun_discards_previous_solution_and_plot() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        app.borrow_mut()
            .start_test_calculation_worker(Duration::from_millis(10), true);
        for _ in 0..40 {
            harness.step();
            if matches!(
                app.borrow().calculation_state(),
                ReactorIvpCalculationState::Completed | ReactorIvpCalculationState::Failed
            ) {
                break;
            }
            std::thread::sleep(Duration::from_millis(10));
        }
        assert!(matches!(
            app.borrow().calculation_state(),
            ReactorIvpCalculationState::Completed
        ));
        assert!(app.borrow().task().last_solve_snapshot.is_some());
        assert!(app.borrow().has_plot_window());

        app.borrow_mut().config.atol = -1.0;
        app.borrow_mut().request_calculation();
        assert!(app.borrow().task().last_solve_snapshot.is_none());
        assert!(!app.borrow().has_plot_window());
        for _ in 0..40 {
            harness.step();
            if matches!(
                app.borrow().calculation_state(),
                ReactorIvpCalculationState::Completed | ReactorIvpCalculationState::Failed
            ) {
                break;
            }
            std::thread::sleep(Duration::from_millis(10));
        }

        let app_ref = app.borrow();
        assert!(matches!(
            app_ref.calculation_state(),
            ReactorIvpCalculationState::Failed
        ));
        assert!(app_ref.task().last_solve_snapshot.is_none());
        assert!(app_ref.task().last_conservation_report.is_none());
        assert!(app_ref.task().solution_render_data().is_err());
        assert!(!app_ref.has_plot_window());
        assert!(app_ref.last_error().is_some());
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
