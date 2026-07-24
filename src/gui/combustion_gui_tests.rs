#[cfg(test)]
mod tests {
    use super::super::combustion::{
        CombustionApp, GuiFileOperationKind, GuiFileOperationResult, ProblemsEnum,
        insert_field_if_absent, insert_section_if_absent, normalize_bvp_usize_fields,
        sorted_document_section_names, sorted_section_field_names,
    };
    use crate::ReactorsBVP::task_parser_reactor_BVP::normalize_reactor_physics_task_map;
    use RustedSciThe::command_interpreter::task_parser::{DocumentMap, DocumentParser, Value};
    use RustedSciThe::numerical::BVP_Damp::task_parser_damped;
    use RustedSciThe::symbolic::codegen::codegen_aot_driver::AotCodegenBackend;
    use RustedSciThe::symbolic::codegen::codegen_backend_selection::BackendSelectionPolicy;
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::cell::RefCell;
    use std::collections::HashMap;
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

        assert_eq!(
            as_string(first_value(&app.document, "solver_settings", "scheme")),
            "forward"
        );
        assert_eq!(
            as_string(first_value(&app.document, "solver_settings", "method")),
            "Banded"
        );
        assert_eq!(
            as_string(first_value(
                &app.document,
                "solver_settings",
                "generated_backend"
            )),
            "banded_lambdify"
        );
        assert_eq!(
            as_string(first_value(
                &app.document,
                "solver_settings",
                "matrix_backend"
            )),
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
            as_string(first_value(
                &app.document,
                "solver_settings",
                "symbolic_backend"
            )),
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
    fn bvp_canonical_render_keeps_structured_document_stable() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        let document_before = app.borrow().document_to_string();
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("Solver backend");
        harness.get_by_label("Adaptive grid refinement");

        assert_eq!(app.borrow().document_to_string(), document_before);
    }

    #[test]
    fn raw_compatibility_fields_are_sorted_for_stable_rendering() {
        let mut section: HashMap<String, Option<Vec<Value>>> = HashMap::new();
        section.insert(
            "zeta".to_string(),
            Some(vec![Value::String("late".to_string())]),
        );
        section.insert(
            "alpha".to_string(),
            Some(vec![Value::String("first".to_string())]),
        );
        section.insert("middle".to_string(), None);

        assert_eq!(
            sorted_section_field_names(&section),
            vec!["alpha", "middle", "zeta"]
        );
    }

    #[test]
    fn raw_compatibility_sections_are_sorted_for_stable_rendering() {
        let mut document: DocumentMap = HashMap::new();
        document.insert("zeta".to_string(), HashMap::new());
        document.insert("alpha".to_string(), HashMap::new());
        document.insert("middle".to_string(), HashMap::new());

        assert_eq!(
            sorted_document_section_names(&document),
            vec!["alpha", "middle", "zeta"]
        );
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
        assert!(
            app_ref
                .document
                .get("grid_refinement")
                .is_some_and(|section| section.len() == 1)
        );
        assert!(
            !app_ref
                .document
                .get("strategy_params")
                .is_some_and(|section| section.contains_key("adaptive"))
        );

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
        assert!(
            !app_ref
                .document
                .get("solver_settings")
                .is_some_and(|section| section.contains_key("aot_build_policy"))
        );
    }

    #[test]
    fn bvp_aot_presets_preserve_typed_codegen_and_compiler_contracts() {
        struct Case {
            preset: &'static str,
            matrix: &'static str,
            explicit_codegen: Option<&'static str>,
            expected_codegen: AotCodegenBackend,
            expected_compiler: Option<&'static str>,
        }

        fn canonical_generated_backend(case: &Case) -> String {
            let matrix_prefix = match case.matrix {
                "Sparse" => "sparse",
                _ => "banded",
            };

            match case.expected_codegen {
                AotCodegenBackend::C => format!(
                    "{}_aot_{}",
                    matrix_prefix,
                    case.expected_compiler.unwrap_or("tcc")
                ),
                AotCodegenBackend::Zig => format!("{}_aot_zig", matrix_prefix),
                AotCodegenBackend::Rust => format!("{}_aot", matrix_prefix),
            }
        }

        let cases = [
            Case {
                preset: "banded_aot",
                matrix: "Banded",
                explicit_codegen: None,
                expected_codegen: AotCodegenBackend::C,
                expected_compiler: Some("tcc"),
            },
            Case {
                preset: "sparse_aot",
                matrix: "Sparse",
                explicit_codegen: None,
                expected_codegen: AotCodegenBackend::C,
                expected_compiler: Some("tcc"),
            },
            Case {
                preset: "banded_aot_gcc",
                matrix: "Banded",
                explicit_codegen: None,
                expected_codegen: AotCodegenBackend::C,
                expected_compiler: Some("gcc"),
            },
            Case {
                preset: "sparse_aot_tcc",
                matrix: "Sparse",
                explicit_codegen: None,
                expected_codegen: AotCodegenBackend::C,
                expected_compiler: Some("tcc"),
            },
            Case {
                preset: "sparse_aot_gcc",
                matrix: "Sparse",
                explicit_codegen: None,
                expected_codegen: AotCodegenBackend::C,
                expected_compiler: Some("gcc"),
            },
            Case {
                preset: "banded_aot_zig",
                matrix: "Banded",
                explicit_codegen: None,
                expected_codegen: AotCodegenBackend::Zig,
                expected_compiler: None,
            },
            Case {
                preset: "sparse_aot_zig",
                matrix: "Sparse",
                explicit_codegen: None,
                expected_codegen: AotCodegenBackend::Zig,
                expected_compiler: None,
            },
            Case {
                preset: "banded_aot",
                matrix: "Banded",
                explicit_codegen: Some("Rust"),
                expected_codegen: AotCodegenBackend::Rust,
                expected_compiler: None,
            },
            Case {
                preset: "sparse_aot",
                matrix: "Sparse",
                explicit_codegen: Some("Rust"),
                expected_codegen: AotCodegenBackend::Rust,
                expected_compiler: None,
            },
        ];

        for case in cases {
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
                    "generated_backend".to_string(),
                    Some(vec![Value::String(case.preset.to_string())]),
                );
                solver_settings.insert(
                    "method".to_string(),
                    Some(vec![Value::String(case.matrix.to_string())]),
                );
                solver_settings.insert(
                    "matrix_backend".to_string(),
                    Some(vec![Value::String(case.matrix.to_string())]),
                );
                solver_settings.remove("aot_codegen_backend");
                solver_settings.remove("aot_c_compiler");
                if let Some(codegen) = case.explicit_codegen {
                    solver_settings.insert(
                        "aot_codegen_backend".to_string(),
                        Some(vec![Value::String(codegen.to_string())]),
                    );
                }
            }

            let mut open = true;
            let app_for_ui = Rc::clone(&app);
            let mut harness = Harness::new_ui(move |ui| {
                app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
            });
            harness.run();

            let app_ref = app.borrow();
            let expected_generated_backend = canonical_generated_backend(&case);
            assert_eq!(
                as_string(first_value(
                    &app_ref.document,
                    "solver_settings",
                    "generated_backend"
                )),
                expected_generated_backend,
                "rendering changed canonical backend for {}",
                case.preset
            );
            let normalized = normalize_reactor_physics_task_map(&app_ref.document);
            let spec =
                task_parser_damped::parse_bvp_damped_solver_settings_from_document(&normalized)
                    .unwrap_or_else(|error| {
                        panic!("preset {} failed RST parsing: {error}", case.preset)
                    });
            let options = spec
                .build_solver_options()
                .unwrap_or_else(|error| panic!("preset {} failed RST build: {error}", case.preset));
            assert_eq!(
                options.generated_backend_config.backend_policy_override,
                Some(BackendSelectionPolicy::AotOnly),
                "wrong execution policy for {}",
                case.preset
            );
            assert_eq!(
                options.generated_backend_config.aot_codegen_backend, case.expected_codegen,
                "wrong codegen backend for {}",
                case.preset
            );
            assert_eq!(
                options.generated_backend_config.aot_c_compiler.as_deref(),
                case.expected_compiler,
                "wrong C compiler for {}",
                case.preset
            );

            let serialized = app_ref.document_to_string();
            let mut parser = DocumentParser::new(serialized);
            parser
                .parse_document()
                .unwrap_or_else(|error| panic!("preset {} failed roundtrip: {error}", case.preset));
            let reparsed = parser
                .get_result()
                .expect("roundtrip parser should return a document");
            assert_eq!(
                as_string(first_value(
                    reparsed,
                    "solver_settings",
                    "generated_backend"
                )),
                expected_generated_backend,
                "serialized roundtrip changed canonical backend for {}",
                case.preset
            );
        }
    }

    #[test]
    fn bvp_generated_backend_presets_roundtrip_through_load_render_save() {
        struct Case {
            preset: &'static str,
            matrix: &'static str,
            expected_codegen: Option<&'static str>,
            expected_compiler: Option<&'static str>,
        }

        fn seed_generated_backend_case(app: &mut CombustionApp, case: &Case) {
            let solver_settings = app
                .document
                .get_mut("solver_settings")
                .expect("solver_settings section must exist");
            solver_settings.insert(
                "generated_backend".to_string(),
                Some(vec![Value::String(case.preset.to_string())]),
            );
            solver_settings.insert(
                "method".to_string(),
                Some(vec![Value::String(case.matrix.to_string())]),
            );
            solver_settings.insert(
                "matrix_backend".to_string(),
                Some(vec![Value::String(case.matrix.to_string())]),
            );
            solver_settings.insert(
                "symbolic_backend".to_string(),
                Some(vec![Value::String("AtomView".to_string())]),
            );
            solver_settings.remove("aot_codegen_backend");
            solver_settings.remove("aot_c_compiler");

            if let Some(codegen) = case.expected_codegen {
                solver_settings.insert(
                    "aot_codegen_backend".to_string(),
                    Some(vec![Value::String(codegen.to_string())]),
                );
            }
            if let Some(compiler) = case.expected_compiler {
                solver_settings.insert(
                    "aot_c_compiler".to_string(),
                    Some(vec![Value::String(compiler.to_string())]),
                );
            }
        }

        fn render_bvp_frame(app: Rc<RefCell<CombustionApp>>) {
            let mut open = true;
            let app_for_ui = Rc::clone(&app);
            let mut harness = Harness::new_ui(move |ui| {
                app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
            });

            harness.run();
            harness.get_by_label("Residuals and Jacobian Backend");
            harness.get_by_label("Execution and assembly");
        }

        let cases = [
            Case {
                preset: "banded_lambdify",
                matrix: "Banded",
                expected_codegen: None,
                expected_compiler: None,
            },
            Case {
                preset: "sparse_lambdify",
                matrix: "Sparse",
                expected_codegen: None,
                expected_compiler: None,
            },
            Case {
                preset: "banded_aot_tcc",
                matrix: "Banded",
                expected_codegen: Some("C"),
                expected_compiler: Some("tcc"),
            },
            Case {
                preset: "sparse_aot_tcc",
                matrix: "Sparse",
                expected_codegen: Some("C"),
                expected_compiler: Some("tcc"),
            },
            Case {
                preset: "banded_aot_gcc",
                matrix: "Banded",
                expected_codegen: Some("C"),
                expected_compiler: Some("gcc"),
            },
            Case {
                preset: "sparse_aot_gcc",
                matrix: "Sparse",
                expected_codegen: Some("C"),
                expected_compiler: Some("gcc"),
            },
            Case {
                preset: "banded_aot_zig",
                matrix: "Banded",
                expected_codegen: Some("Zig"),
                expected_compiler: None,
            },
            Case {
                preset: "sparse_aot_zig",
                matrix: "Sparse",
                expected_codegen: Some("Zig"),
                expected_compiler: None,
            },
            Case {
                preset: "banded_aot",
                matrix: "Banded",
                expected_codegen: Some("Rust"),
                expected_compiler: None,
            },
            Case {
                preset: "sparse_aot",
                matrix: "Sparse",
                expected_codegen: Some("Rust"),
                expected_compiler: None,
            },
        ];

        for case in cases {
            let temp_dir = tempdir().expect("temp dir should be available");
            let first_path = temp_dir.path().join(format!("{}_first.txt", case.preset));
            let second_path = temp_dir.path().join(format!("{}_second.txt", case.preset));

            let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
                ProblemsEnum::BVPSimple,
            )));
            {
                let mut app_ref = app.borrow_mut();
                seed_generated_backend_case(&mut app_ref, &case);
            }

            render_bvp_frame(Rc::clone(&app));

            match app.borrow_mut().save_document(first_path.clone()) {
                GuiFileOperationResult::Saved { path } => assert_eq!(path, first_path),
                other => panic!(
                    "expected first save to succeed for preset {}, got {other:?}",
                    case.preset
                ),
            }
            let first_saved =
                std::fs::read_to_string(&first_path).expect("first save should be readable");

            let reloaded = Rc::new(RefCell::new(CombustionApp::new_with_problem(
                ProblemsEnum::BVPSimple,
            )));
            match reloaded
                .borrow_mut()
                .load_document_from_path(first_path.clone())
            {
                GuiFileOperationResult::Loaded { path } => assert_eq!(path, first_path),
                other => panic!(
                    "expected reload to succeed for preset {}, got {other:?}",
                    case.preset
                ),
            }

            render_bvp_frame(Rc::clone(&reloaded));

            match reloaded.borrow_mut().save_document(second_path.clone()) {
                GuiFileOperationResult::Saved { path } => assert_eq!(path, second_path),
                other => panic!(
                    "expected second save to succeed for preset {}, got {other:?}",
                    case.preset
                ),
            }
            let second_saved =
                std::fs::read_to_string(&second_path).expect("second save should be readable");

            assert_eq!(
                first_saved, second_saved,
                "load-render-save changed the serialized document for preset {}",
                case.preset
            );

            let reloaded_ref = reloaded.borrow();
            assert_eq!(
                as_string(first_value(
                    &reloaded_ref.document,
                    "solver_settings",
                    "generated_backend"
                )),
                case.preset,
                "generated_backend drifted for preset {}",
                case.preset
            );
            assert_eq!(
                as_string(first_value(
                    &reloaded_ref.document,
                    "solver_settings",
                    "method"
                )),
                case.matrix,
                "matrix backend drifted for preset {}",
                case.preset
            );
            assert_eq!(
                as_string(first_value(
                    &reloaded_ref.document,
                    "solver_settings",
                    "matrix_backend"
                )),
                case.matrix,
                "matrix_backend drifted for preset {}",
                case.preset
            );
            assert_eq!(
                as_string(first_value(
                    &reloaded_ref.document,
                    "solver_settings",
                    "backend_policy"
                )),
                if case.preset.contains("lambdify") {
                    "lambdify_only"
                } else {
                    "aot_only"
                },
                "backend policy drifted for preset {}",
                case.preset
            );

            let solver_settings = reloaded_ref
                .document
                .get("solver_settings")
                .expect("solver_settings section should still exist");
            if let Some(codegen) = case.expected_codegen {
                assert_eq!(
                    as_string(first_value(
                        &reloaded_ref.document,
                        "solver_settings",
                        "aot_codegen_backend"
                    )),
                    codegen,
                    "aot_codegen_backend drifted for preset {}",
                    case.preset
                );
            } else {
                assert!(
                    !solver_settings.contains_key("aot_codegen_backend"),
                    "lambdify preset {} should not keep AOT codegen controls",
                    case.preset
                );
            }

            if let Some(compiler) = case.expected_compiler {
                assert_eq!(
                    as_string(first_value(
                        &reloaded_ref.document,
                        "solver_settings",
                        "aot_c_compiler"
                    )),
                    compiler,
                    "aot_c_compiler drifted for preset {}",
                    case.preset
                );
            } else {
                assert!(
                    !solver_settings.contains_key("aot_c_compiler"),
                    "preset {} should not keep a C compiler field",
                    case.preset
                );
            }
        }
    }

    #[test]
    fn bvp_screen_reports_missing_structured_sections_without_recreating_them() {
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
        harness.get_by_label("Missing canonical section: initial_guess");
        harness.get_by_label("Missing canonical section: postprocessing");

        let app_ref = app.borrow();
        assert!(!app_ref.document.contains_key("initial_guess"));
        assert!(!app_ref.document.contains_key("postprocessing"));
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
                Some(vec![
                    Value::Float(0.05),
                    Value::Float(0.05),
                    Value::Float(1.25),
                ]),
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
            as_vector(first_value(
                &app_ref.document,
                "grid_refinement",
                "grcarsmooke"
            )),
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
            solver_settings.insert("max_iterations".to_string(), Some(vec![Value::Integer(77)]));
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
            other => panic!(
                "max_iterations should accept string input and normalize to Value::Usize, got {other:?}"
            ),
        }
    }

    #[test]
    fn bvp_numeric_solver_fields_normalize_all_counter_slots_without_touching_tolerance() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        {
            let process_conditions = app
                .document
                .get_mut("process_conditions")
                .expect("process_conditions section must exist");
            process_conditions.insert(
                "n_steps".to_string(),
                Some(vec![Value::String("256".to_string())]),
            );

            let solver_settings = app
                .document
                .get_mut("solver_settings")
                .expect("solver_settings section must exist");
            solver_settings.insert(
                "max_iterations".to_string(),
                Some(vec![Value::Integer(150)]),
            );
            solver_settings.insert(
                "refinement_steps".to_string(),
                Some(vec![Value::Float(6.0)]),
            );
            solver_settings.insert(
                "abs_tolerance".to_string(),
                Some(vec![Value::Float(2.5e-6)]),
            );

            let strategy_params = app
                .document
                .get_mut("strategy_params")
                .expect("strategy_params section must exist");
            strategy_params.insert(
                "max_jac".to_string(),
                Some(vec![Value::Optional(Some(Box::new(Value::Float(4.0))))]),
            );
            strategy_params.insert(
                "max_damp_iter".to_string(),
                Some(vec![Value::Optional(Some(Box::new(Value::String(
                    "12".to_string(),
                ))))]),
            );
        }

        let adaptive_strategy = app
            .document
            .get_mut("adaptive_strategy")
            .expect("adaptive_strategy section must exist");
        adaptive_strategy.insert("version".to_string(), Some(vec![Value::Usize(1)]));
        adaptive_strategy.insert("max_refinements".to_string(), Some(vec![Value::Integer(8)]));

        normalize_bvp_usize_fields(&mut app.document)
            .expect("all numeric solver counters should normalize losslessly");

        assert_eq!(
            as_usize(first_value(&app.document, "process_conditions", "n_steps")),
            256
        );
        assert_eq!(
            as_usize(first_value(
                &app.document,
                "solver_settings",
                "max_iterations"
            )),
            150
        );
        assert_eq!(
            as_usize(first_value(
                &app.document,
                "solver_settings",
                "refinement_steps"
            )),
            6
        );
        assert_eq!(
            as_usize(first_value(&app.document, "strategy_params", "max_jac")),
            4
        );
        assert_eq!(
            as_usize(first_value(
                &app.document,
                "strategy_params",
                "max_damp_iter"
            )),
            12
        );
        assert_eq!(
            as_usize(first_value(
                &app.document,
                "adaptive_strategy",
                "max_refinements"
            )),
            8
        );
        assert!(matches!(
            first_value(&app.document, "solver_settings", "abs_tolerance"),
            Value::Float(value) if *value == 2.5e-6
        ));
    }

    #[test]
    fn bvp_process_conditions_n_steps_rejects_fractional_value_without_truncation() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document
            .get_mut("process_conditions")
            .expect("process_conditions section should exist")
            .insert("n_steps".to_string(), Some(vec![Value::Float(12.5)]));

        let error = normalize_bvp_usize_fields(&mut app.document)
            .expect_err("fractional grid counters must be rejected");

        assert!(error.contains("process_conditions.n_steps"));
        assert!(error.contains("without truncation"));
        assert_eq!(
            first_value(&app.document, "process_conditions", "n_steps"),
            &Value::Float(12.5),
            "failed conversion must preserve the original grid count"
        );
    }

    #[test]
    fn bvp_integer_handoff_rejects_fractional_counter_without_truncation() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document
            .get_mut("solver_settings")
            .expect("solver_settings should exist")
            .insert("max_iterations".to_string(), Some(vec![Value::Float(1.9)]));

        let error = normalize_bvp_usize_fields(&mut app.document)
            .expect_err("fractional solver counters must be rejected");

        assert!(error.contains("solver_settings.max_iterations"));
        assert!(error.contains("without truncation"));
        assert_eq!(
            first_value(&app.document, "solver_settings", "max_iterations"),
            &Value::Float(1.9),
            "failed conversion must preserve the user's original value"
        );
    }

    #[test]
    fn bvp_solver_backend_reports_invalid_counter_instead_of_showing_default() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        app.borrow_mut()
            .document
            .get_mut("solver_settings")
            .expect("solver_settings should exist")
            .insert(
                "refinement_steps".to_string(),
                Some(vec![Value::Float(f64::NAN)]),
            );

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();

        harness.get_by_label("non-finite float `NaN` cannot represent usize");
        assert!(matches!(
            first_value(&app.borrow().document, "solver_settings", "refinement_steps"),
            Value::Float(value) if value.is_nan()
        ));
    }

    #[test]
    fn bvp_loglevel_uses_typed_optional_string_shape() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        app.borrow_mut()
            .document
            .get_mut("solver_settings")
            .expect("solver_settings should exist")
            .insert(
                "loglevel".to_string(),
                Some(vec![Value::Optional(Some(Box::new(Value::Float(0.0))))]),
            );

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();

        assert!(matches!(
            first_value(&app.borrow().document, "solver_settings", "loglevel"),
            Value::Optional(None)
        ));
    }

    #[test]
    fn bvp_loglevel_selector_restores_some_with_a_string_payload() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        app.borrow_mut()
            .document
            .get_mut("solver_settings")
            .expect("solver_settings should exist")
            .insert("loglevel".to_string(), Some(vec![Value::Optional(None)]));

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "loglevel")
            .click();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "debug").click();
        harness.run();

        assert!(app.borrow().document_lifecycle.dirty);
        assert!(matches!(
            first_value(&app.borrow().document, "solver_settings", "loglevel"),
            Value::Optional(Some(inner))
                if matches!(inner.as_ref(), Value::String(value) if value == "debug")
        ));
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
        let save_result = app.borrow_mut().save_document(path.clone());
        match save_result {
            GuiFileOperationResult::Saved { path: saved_path } => assert_eq!(saved_path, path),
            other => panic!("expected save success, got {other:?}"),
        }
        assert!(!app.borrow().document_lifecycle.dirty);
        assert_eq!(
            app.borrow().document_lifecycle.last_status_text(),
            Some(format!("Saved document to {}.", path.display()))
        );

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

    #[test]
    fn bvp_document_save_failure_keeps_current_path_and_reports_status() {
        let temp_dir = tempdir().expect("temp dir should be available");
        let failing_target = temp_dir.path().to_path_buf();
        let previous_path = temp_dir.path().join("existing_task.txt");

        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.current_file_path = Some(previous_path.clone());
        app.document_lifecycle.current_path = Some(previous_path.clone());
        app.document_lifecycle.dirty = true;

        let result = app.save_document(failing_target.clone());

        match result {
            GuiFileOperationResult::Failed {
                kind: GuiFileOperationKind::Save,
                path: Some(path),
                ..
            } => assert_eq!(path, failing_target),
            other => panic!("expected save failure, got {other:?}"),
        }
        assert_eq!(app.current_file_path, Some(previous_path.clone()));
        assert_eq!(app.document_lifecycle.current_path, Some(previous_path));
        assert!(app.document_lifecycle.dirty);
        assert!(app.document_lifecycle.last_status_is_error());
    }

    #[test]
    fn bvp_document_load_failure_keeps_current_path_and_document() {
        let temp_dir = tempdir().expect("temp dir should be available");
        let missing_path = temp_dir.path().join("missing_task.txt");
        let previous_path = temp_dir.path().join("existing_task.txt");

        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let document_before = app.document_to_string();
        app.current_file_path = Some(previous_path.clone());
        app.document_lifecycle.current_path = Some(previous_path.clone());

        let result = app.load_document_from_path(missing_path.clone());

        match result {
            GuiFileOperationResult::Failed {
                kind: GuiFileOperationKind::Load,
                path: Some(path),
                ..
            } => assert_eq!(path, missing_path),
            other => panic!("expected load failure, got {other:?}"),
        }
        assert_eq!(app.current_file_path, Some(previous_path.clone()));
        assert_eq!(app.document_lifecycle.current_path, Some(previous_path));
        assert_eq!(app.document_to_string(), document_before);
        assert!(app.document_lifecycle.last_status_is_error());
    }

    #[test]
    fn bvp_document_load_success_updates_lifecycle_state() {
        let temp_dir = tempdir().expect("temp dir should be available");
        let path = temp_dir.path().join("valid_task.txt");
        let template = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        std::fs::write(&path, template.document_to_string())
            .expect("valid task should be writable");

        let mut app = CombustionApp::new_with_problem(ProblemsEnum::None);
        let result = app.load_document_from_path(path.clone());

        match result {
            GuiFileOperationResult::Loaded { path: loaded_path } => assert_eq!(loaded_path, path),
            other => panic!("expected load success, got {other:?}"),
        }
        assert_eq!(app.current_file_path, Some(path.clone()));
        assert_eq!(app.document_lifecycle.current_path, Some(path));
        assert!(!app.document_lifecycle.dirty);
        assert!(!app.document_lifecycle.last_status_is_error());
    }

    #[test]
    fn bvp_document_load_requests_confirmation_when_document_is_dirty() {
        let temp_dir = tempdir().expect("temp dir should be available");
        let path = temp_dir.path().join("valid_task.txt");
        let template = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        std::fs::write(&path, template.document_to_string())
            .expect("valid task should be writable");

        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let mut section = HashMap::new();
        section.insert(
            "note".to_string(),
            Some(vec![Value::String("dirty".to_string())]),
        );
        app.document
            .insert("ephemeral_ui_edit".to_string(), section);
        let expected_document = app.document_to_string();

        let requested = app.request_document_load(path.clone());

        assert!(!requested);
        assert_eq!(app.pending_document_load, Some(path.clone()));
        assert_eq!(app.document_to_string(), expected_document);
        assert!(app.current_file_path.is_none());
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("Confirm loading"))
        );
        assert!(!app.last_run_is_error);

        assert!(app.confirm_pending_document_load());
        assert!(app.pending_document_load.is_none());
        assert_eq!(app.current_file_path, Some(path.clone()));
        assert_eq!(app.document_lifecycle.current_path, Some(path));
        assert!(!app.document_lifecycle.dirty);
    }

    #[test]
    fn bvp_document_load_applies_immediately_when_document_is_clean() {
        let temp_dir = tempdir().expect("temp dir should be available");
        let path = temp_dir.path().join("valid_task.txt");
        let template = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        std::fs::write(&path, template.document_to_string())
            .expect("valid task should be writable");

        let mut app = CombustionApp::new_with_problem(ProblemsEnum::None);

        let requested = app.request_document_load(path.clone());

        assert!(requested);
        assert!(app.pending_document_load.is_none());
        assert_eq!(app.current_file_path, Some(path.clone()));
        assert_eq!(app.document_lifecycle.current_path, Some(path));
        assert!(!app.document_lifecycle.dirty);
        assert!(!app.document_lifecycle.last_status_is_error());
    }

    #[test]
    fn bvp_window_close_requests_confirmation_when_document_is_dirty() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document
            .entry("ephemeral_ui_edit".to_string())
            .or_insert_with(HashMap::new)
            .insert(
                "note".to_string(),
                Some(vec![Value::String("dirty".to_string())]),
            );

        let allowed = app.request_window_close();

        assert!(!allowed);
        assert!(app.pending_window_close);
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("Confirm closing"))
        );
        assert!(!app.last_run_is_error);

        app.cancel_pending_window_close();

        assert!(!app.pending_window_close);
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("Close cancelled"))
        );
        assert!(!app.last_run_is_error);
    }

    #[test]
    fn bvp_window_close_applies_immediately_when_document_is_clean() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::None);

        let allowed = app.request_window_close();

        assert!(allowed);
        assert!(!app.pending_window_close);
        assert!(app.last_run_message.is_none());
        assert!(!app.document_lifecycle.dirty);
    }

    #[test]
    fn bvp_window_close_detects_direct_document_mutations_without_render_sync() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let mut section = HashMap::new();
        section.insert(
            "note".to_string(),
            Some(vec![Value::String("dirty".to_string())]),
        );
        app.document
            .insert("ephemeral_ui_edit".to_string(), section);

        let allowed = app.request_window_close();

        assert!(!allowed);
        assert!(app.pending_window_close);
        assert!(app.document_lifecycle.dirty);
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("Confirm closing"))
        );
    }

    #[test]
    fn duplicate_sections_and_fields_preserve_existing_document_data() {
        let mut document = DocumentMap::new();
        insert_section_if_absent(&mut document, "physics".to_string())
            .expect("first section insertion should succeed");
        let section = document
            .get_mut("physics")
            .expect("new section should exist");
        insert_field_if_absent(section, "temperature".to_string(), Value::Float(900.0))
            .expect("first field insertion should succeed");

        let duplicate_field =
            insert_field_if_absent(section, "temperature".to_string(), Value::Float(1200.0));
        assert!(duplicate_field.is_err());
        assert_eq!(
            section
                .get("temperature")
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(&Value::Float(900.0))
        );

        let duplicate_section = insert_section_if_absent(&mut document, "physics".to_string());
        assert!(duplicate_section.is_err());
        assert!(
            document
                .get("physics")
                .is_some_and(|section| section.contains_key("temperature"))
        );
    }

    #[test]
    fn replacing_problem_clears_file_binding_without_overwriting_previous_task() {
        let temp_dir = tempdir().expect("temp dir should be available");
        let path = temp_dir.path().join("existing_bvp_task.txt");
        std::fs::write(&path, "original task contents").expect("test task should be writable");

        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.current_file_path = Some(path.clone());
        app.new_section_name = "pending section".to_string();
        app.new_field_names
            .insert("physics".to_string(), "pending field".to_string());
        app.last_run_message = Some("stale result".to_string());
        app.last_run_is_error = true;

        app.replace_problem(ProblemsEnum::None);

        assert_eq!(app.selected_problem, ProblemsEnum::None);
        assert!(app.document.is_empty());
        assert!(app.current_file_path.is_none());
        assert!(app.new_section_name.is_empty());
        assert!(app.new_field_names.is_empty());
        assert!(app.last_run_message.is_none());
        assert!(!app.last_run_is_error);
        assert_eq!(
            std::fs::read_to_string(path).expect("original task should remain readable"),
            "original task contents"
        );
    }

    #[test]
    fn request_problem_switch_defers_when_document_is_dirty() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let mut section = HashMap::new();
        section.insert(
            "note".to_string(),
            Some(vec![Value::String("dirty".to_string())]),
        );
        app.document
            .insert("ephemeral_ui_edit".to_string(), section);
        let before = app.document_to_string();

        let switched = app.request_problem_switch(ProblemsEnum::None);

        assert!(!switched);
        assert_eq!(app.selected_problem, ProblemsEnum::BVPSimple);
        assert_eq!(app.pending_problem_switch, Some(ProblemsEnum::None));
        assert_eq!(app.document_to_string(), before);
        assert!(
            app.last_run_message
                .as_deref()
                .is_some_and(|message| message.contains("Confirm switching"))
        );
        assert!(!app.last_run_is_error);
    }

    #[test]
    fn request_problem_switch_applies_immediately_when_document_is_clean() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);

        let switched = app.request_problem_switch(ProblemsEnum::None);

        assert!(switched);
        assert_eq!(app.selected_problem, ProblemsEnum::None);
        assert!(app.pending_problem_switch.is_none());
        assert!(app.document.is_empty());
        assert!(app.current_file_path.is_none());
        assert!(!app.document_lifecycle.dirty);
    }

    #[test]
    fn confirm_pending_problem_switch_replaces_the_active_problem() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let mut section = HashMap::new();
        section.insert(
            "note".to_string(),
            Some(vec![Value::String("dirty".to_string())]),
        );
        app.document
            .insert("ephemeral_ui_edit".to_string(), section);

        assert!(!app.request_problem_switch(ProblemsEnum::None));
        assert_eq!(app.pending_problem_switch, Some(ProblemsEnum::None));

        assert!(app.confirm_pending_problem_switch());
        assert_eq!(app.selected_problem, ProblemsEnum::None);
        assert!(app.pending_problem_switch.is_none());
        assert!(app.document.is_empty());
        assert!(app.current_file_path.is_none());
        assert!(!app.document_lifecycle.dirty);
    }
}
