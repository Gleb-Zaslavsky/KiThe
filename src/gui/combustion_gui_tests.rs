#[cfg(test)]
mod tests {
    use super::super::combustion::{CombustionApp, ProblemsEnum};
    use crate::ReactorsBVP::task_parser_reactor_BVP::normalize_reactor_physics_task_map;
    use egui_kittest::kittest::Queryable;
    use egui_kittest::Harness;
    use RustedSciThe::command_interpreter::task_parser::{DocumentMap, DocumentParser, Value};
    use RustedSciThe::numerical::BVP_Damp::task_parser_damped;
    use std::cell::RefCell;
    use std::rc::Rc;
    use tempfile::tempdir;

    fn first_value<'a>(document: &'a DocumentMap, section: &str, field: &str) -> &'a Value {
        document
            .get(section)
            .and_then(|section_map| section_map.get(field))
            .and_then(|slot| slot.as_ref())
            .and_then(|values| values.first())
            .expect("expected BVP solver field to exist")
    }

    fn as_string(value: &Value) -> &str {
        match value {
            Value::String(value) => value,
            other => panic!("expected Value::String, got {other:?}"),
        }
    }

    fn as_bool(value: &Value) -> bool {
        match value {
            Value::Boolean(value) => *value,
            other => panic!("expected Value::Boolean, got {other:?}"),
        }
    }

    fn as_usize(value: &Value) -> usize {
        match value {
            Value::Usize(value) => *value,
            other => panic!("expected Value::Usize, got {other:?}"),
        }
    }

    fn as_vector(value: &Value) -> &[f64] {
        match value {
            Value::Vector(value) => value,
            other => panic!("expected Value::Vector, got {other:?}"),
        }
    }

    fn load_problem_test_document() -> DocumentMap {
        let content = std::fs::read_to_string("problem_test.txt")
            .expect("problem_test.txt must be readable from the crate root");
        let mut parser = DocumentParser::new(content);
        parser
            .parse_document()
            .expect("problem_test.txt must parse as a valid BVP document");
        parser
            .get_result()
            .cloned()
            .expect("parsed problem_test.txt must yield a document map")
    }

    #[test]
    fn bvp_solver_backend_defaults_are_seeded() {
        let app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let solver_settings = app
            .document
            .get("solver_settings")
            .expect("solver_settings section must be seeded for the BVP template");

        assert_eq!(as_string(first_value(&app.document, "solver_settings", "scheme")), "forward");
        assert_eq!(as_string(first_value(&app.document, "solver_settings", "method")), "Banded");
        assert_eq!(
            as_string(first_value(
                &app.document,
                "solver_settings",
                "generated_backend"
            )),
            "banded_lambdify"
        );
        assert_eq!(
            as_string(first_value(&app.document, "solver_settings", "matrix_backend")),
            "Banded"
        );
        assert_eq!(
            as_usize(first_value(
                &app.document,
                "solver_settings",
                "refinement_steps"
            )),
            5
        );
        assert_eq!(
            as_string(first_value(&app.document, "solver_settings", "symbolic_backend")),
            "AtomView"
        );
        assert!(matches!(
            first_value(&app.document, "solver_settings", "linear_sys_method"),
            Value::Optional(None)
        ));
        assert!(as_bool(first_value(
            &app.document,
            "solver_settings",
            "dont_save_log"
        )));
        assert_eq!(
            as_string(first_value(
                &app.document,
                "solver_settings",
                "backend_policy"
            )),
            "lambdify_only"
        );
        assert!(!solver_settings.contains_key("aot_c_compiler"));
        assert!(!solver_settings.contains_key("aot_build_policy"));
        assert!(!solver_settings.contains_key("aot_execution_policy"));
    }

    #[test]
    fn bvp_solver_backend_panel_renders_new_controls() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let mut open = true;
        let mut harness = Harness::new_ui(move |ui| {
            app.show(ui.ctx(), &mut open);
        });

        harness.run();

        for label in [
            "Solver backend",
            "Core solve",
            "Execution and assembly",
            "Residuals and Jacobian Backend",
            "scheme",
            "Linear algebra backend",
            "Linear system method override",
            "Symbolic backend",
            "Banded linear solver",
            "refinement_steps",
        ] {
            harness.get_by_label(label);
        }

    }

    #[test]
    fn bvp_solver_backend_keeps_linear_system_override_typed() {
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
                "linear_sys_method".to_string(),
                Some(vec![Value::Optional(Some(Box::new(Value::Float(0.0))))]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("Linear system method override");

        let app_ref = app.borrow();
        assert!(matches!(
            first_value(&app_ref.document, "solver_settings", "linear_sys_method"),
            Value::Optional(None)
        ));
    }

    #[test]
    fn bvp_screen_renders_structured_sections() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let mut open = true;
        let mut harness = Harness::new_ui(move |ui| {
            app.show(ui.ctx(), &mut open);
        });

        harness.run();

        for label in [
            "Physics",
            "Initial guess",
            "Postprocessing",
            "Legacy raw document",
            "process_conditions",
            "boundary_condition",
            "initial_guess",
            "postprocessing",
        ] {
            harness.get_by_label(label);
        }
    }

    #[test]
    fn bvp_screen_renders_advanced_solver_sections_separately_from_physics() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let mut open = true;
        let mut harness = Harness::new_ui(move |ui| {
            app.show(ui.ctx(), &mut open);
        });

        harness.run();

        for label in [
            "Advanced solver settings",
            "rel_tolerance",
            "strategy_params",
            "Adaptive grid refinement",
            "Adaptive strategy version: 1",
            "Grid refinement strategy",
        ] {
            harness.get_by_label(label);
        }
    }

    #[test]
    fn bvp_adaptive_checkbox_builds_complete_native_rst_contract() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        app.borrow_mut().document = load_problem_test_document();

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        assert!(!app.borrow().document.contains_key("adaptive_strategy"));
        assert!(!app.borrow().document.contains_key("grid_refinement"));

        harness.get_by_label("Adaptive grid refinement").click();
        harness.run();

        let app_ref = app.borrow();
        assert_eq!(
            as_usize(first_value(
                &app_ref.document,
                "adaptive_strategy",
                "version"
            )),
            1
        );
        assert_eq!(
            as_usize(first_value(
                &app_ref.document,
                "adaptive_strategy",
                "max_refinements"
            )),
            3
        );
        assert!(app_ref
            .document
            .get("grid_refinement")
            .is_some_and(|section| section.len() == 1));
        assert!(!app_ref
            .document
            .get("strategy_params")
            .is_some_and(|section| section.contains_key("adaptive")));

        let solver_document = normalize_reactor_physics_task_map(&app_ref.document);
        task_parser_damped::parse_bvp_damped_solver_settings_from_document(&solver_document)
            .expect("the GUI must emit a complete adaptive-grid contract accepted by RST");
    }

    #[test]
    fn bvp_matrix_backend_aliases_collapse_to_one_user_choice() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            let solver_settings = app_ref
                .document
                .get_mut("solver_settings")
                .expect("solver_settings section must exist");
            solver_settings.insert(
                "method".to_string(),
                Some(vec![Value::String("Sparse".to_string())]),
            );
            solver_settings.insert(
                "matrix_backend".to_string(),
                Some(vec![Value::String("Banded".to_string())]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Linear algebra backend");

        let app_ref = app.borrow();
        assert_eq!(
            as_string(first_value(&app_ref.document, "solver_settings", "method")),
            "Sparse"
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
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "generated_backend"
            )),
            "sparse_lambdify"
        );
    }

    #[test]
    fn bvp_screen_preserves_custom_backend_values_and_seeds_missing_defaults() {
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
                Some(vec![Value::String("custom_backend".to_string())]),
            );
            solver_settings.remove("matrix_backend");
            solver_settings.remove("aot_build_policy");
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("Solver backend");

        let app_ref = app.borrow();
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "generated_backend"
            )),
            "custom_backend"
        );
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "matrix_backend"
            )),
            "Banded"
        );
        assert!(!app_ref
            .document
            .get("solver_settings")
            .is_some_and(|section| section.contains_key("aot_build_policy")));
    }

    #[test]
    fn bvp_screen_recreates_missing_structured_sections() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            app_ref.document.remove("initial_guess");
            app_ref.document.remove("postprocessing");
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("Initial guess");
        harness.get_by_label("Postprocessing");

        let app_ref = app.borrow();
        assert!(app_ref.document.contains_key("initial_guess"));
        assert!(app_ref.document.contains_key("postprocessing"));
    }

    #[test]
    fn bvp_screen_canonicalizes_legacy_grid_refinement_aliases() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            let grid_refinement = app_ref
                .document
                .get_mut("grid_refinement")
                .expect("grid_refinement section must exist in the BVP template");
            grid_refinement.remove("grcar_smooke");
            grid_refinement.insert(
                "grcarsmooke".to_string(),
                Some(vec![Value::Float(0.05), Value::Float(0.05), Value::Float(1.25)]),
            );
            grid_refinement.insert(
                "double_points".to_string(),
                Some(vec![Value::Vector(vec![])]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("Physics");
        harness.get_by_label("grid_refinement");
        harness.get_by_label("Grid refinement strategy");

        let app_ref = app.borrow();
        let grid_refinement = app_ref
            .document
            .get("grid_refinement")
            .expect("grid_refinement section must still be present");
        assert_eq!(
            grid_refinement.len(),
            1,
            "the structured editor should keep one active refinement strategy"
        );
        assert!(grid_refinement.contains_key("grcarsmooke"));
        assert!(!grid_refinement.contains_key("grcar_smooke"));
        assert!(!grid_refinement.contains_key("double_points"));
        assert_eq!(
            as_vector(first_value(&app_ref.document, "grid_refinement", "grcarsmooke")),
            &[0.05, 0.05, 1.25]
        );
    }

    #[test]
    fn bvp_grid_refinement_defaults_to_single_dropdown_strategy() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("Grid refinement strategy");

        let app_ref = app.borrow();
        let grid_refinement = app_ref
            .document
            .get("grid_refinement")
            .expect("grid_refinement section should exist");
        assert_eq!(
            grid_refinement.len(),
            1,
            "RST accepts one active grid refinement method per task"
        );
    }

    #[test]
    fn bvp_screen_renders_atomic_composition_as_structured_sections() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document = load_problem_test_document();

        let mut open = true;
        let mut harness = Harness::new_ui(move |ui| {
            app.show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("Atomic composition");
        harness.get_by_label("HMX");
        harness.get_by_label("HMXprod");
    }

    #[test]
    fn run_calculation_reports_validation_failure_instead_of_silence() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document.remove("process_conditions");

        app.run_calculation();

        assert!(app.last_run_is_error);
        let message = app
            .last_run_message
            .as_deref()
            .expect("run_calculation should surface an error message");
        assert!(message.contains("Calculation failed"));
    }

    #[test]
    fn bvp_solver_backend_normalizes_refinement_steps_to_usize() {
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
                "refinement_steps".to_string(),
                Some(vec![Value::Float(3.0)]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();

        let app_ref = app.borrow();
        let value = first_value(&app_ref.document, "solver_settings", "refinement_steps");
        match value {
            Value::Usize(step) => assert_eq!(*step, 3),
            other => panic!("refinement_steps should be normalized to Value::Usize, got {other:?}"),
        }
    }

    #[test]
    fn bvp_solver_backend_normalizes_max_iterations_to_usize() {
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
                "max_iterations".to_string(),
                Some(vec![Value::Integer(77)]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();

        let app_ref = app.borrow();
        let value = first_value(&app_ref.document, "solver_settings", "max_iterations");
        match value {
            Value::Usize(step) => assert_eq!(*step, 77),
            other => panic!("max_iterations should be normalized to Value::Usize, got {other:?}"),
        }
    }

    #[test]
    fn bvp_solver_backend_normalizes_string_encoded_max_iterations_to_usize() {
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
                "max_iterations".to_string(),
                Some(vec![Value::String("88".to_string())]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();

        let app_ref = app.borrow();
        let value = first_value(&app_ref.document, "solver_settings", "max_iterations");
        match value {
            Value::Usize(step) => assert_eq!(*step, 88),
            other => panic!("max_iterations should accept string input and normalize to Value::Usize, got {other:?}"),
        }
    }

    #[test]
    fn bvp_document_save_read_save_roundtrip_is_stable() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        {
            let mut app_ref = app.borrow_mut();
            let grid_refinement = app_ref
                .document
                .get_mut("grid_refinement")
                .expect("grid_refinement section should exist");
            grid_refinement.clear();
            grid_refinement.insert(
                "twopnt".to_string(),
                Some(vec![Value::Vector(vec![0.02, 0.03, 1.4])]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label("Grid refinement strategy");

        let temp_dir = tempdir().expect("temp dir should be available");
        let path = temp_dir.path().join("bvp_roundtrip.txt");
        app.borrow().save_document(path.clone());

        let saved = std::fs::read_to_string(&path).expect("saved BVP document should be readable");
        let mut parser = DocumentParser::new(saved.clone());
        parser
            .parse_document()
            .expect("saved BVP document should parse again");

        let mut reloaded = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        reloaded.document = parser
            .get_result()
            .cloned()
            .expect("roundtrip parser should return a document");

        assert_eq!(
            saved,
            reloaded.document_to_string(),
            "save -> read -> save should be a stable reflection of the same BVP task"
        );
    }
}
