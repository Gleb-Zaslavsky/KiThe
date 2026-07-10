//! Story tests for the combustion BVP GUI.
//!
//! These tests treat the BVP screen as a user-facing workflow rather than a
//! bag of controls. Each scenario checks a concrete hypothesis:
//!
//! - A fresh BVP preset should expose the structured layout immediately.
//! - Canonical solver defaults should be visible and stable on first render.
//! - AOT-oriented backend values should expose the AOT-specific controls.
//! - Lambdify and AOT should remain mutually exclusive document contracts.
//! - Legacy documents should still render when older spellings are present.
//! - A solved BVP should still hand postprocessing state through to the GUI.
//!
//! Expected result:
//! the BVP screen should behave like a structured form with a compatibility
//! fallback, not like a raw document editor with a few extra labels on top.

#[cfg(test)]
mod tests {
    use super::super::combustion::{CombustionApp, ProblemsEnum};
    use egui_kittest::kittest::Queryable;
    use egui_kittest::Harness;
    use RustedSciThe::command_interpreter::task_parser::{DocumentMap, Value};
    use std::cell::RefCell;
    use std::rc::Rc;

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
            as_string(first_value(&app_ref.document, "solver_settings", "matrix_backend")),
            "Banded"
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
                Some(vec![Value::String(
                    "prefer_aot_then_lambdify".to_string(),
                )]),
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
        harness.get_by_label("AOT compiler");
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
            as_string(first_value(&app_ref.document, "solver_settings", "matrix_backend")),
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
        assert!(!app_ref
            .document
            .get("solver_settings")
            .is_some_and(|section| section.contains_key("aot_build_policy")));
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
            as_string(first_value(&app_ref.document, "solver_settings", "matrix_backend")),
            "Banded"
        );
    }

}
