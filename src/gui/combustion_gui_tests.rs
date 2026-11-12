#[cfg(test)]
mod tests {
    use super::super::bvp_gui_config::aot_toolchain_request;
    use super::super::bvp_gui_config::{BvpSolverFamily, set_bvp_solver_family};
    use super::super::combustion::{
        BvpEditorTab, CombustionApp, DocumentDeleteRequest, GuiFileOperationKind,
        GuiFileOperationResult, ProblemsEnum, add_bvp_reaction, add_bvp_species,
        apply_document_deletion, insert_field_if_absent, normalize_bvp_usize_fields,
        sorted_document_section_names, sorted_section_field_names,
    };
    use crate::SimpleReactorsBVP::SimpleReactorBVP::{DensityModel, SimpleReactorTask};
    use crate::SimpleReactorsBVP::task_parser_reactor_BVP::normalize_reactor_physics_task_map;
    use RustedSciThe::command_interpreter::task_parser::{DocumentMap, DocumentParser, Value};
    use RustedSciThe::command_interpreter::task_parser_bvp_hybrid::parse_hybrid_bvp_solver_settings_from_document;
    use RustedSciThe::numerical::BVP_CorrectionNewton::{
        CorrectionNewtonMatrixLayout, CorrectionNewtonSymbolicExecution,
        CorrectionNewtonSymbolicFrontend,
    };
    use RustedSciThe::numerical::BVP_Damp::task_parser_damped;
    use RustedSciThe::numerical::BVP_Damp::telemetry::BvpTelemetryMode;
    use RustedSciThe::numerical::BVP_api::BvpSolverSettings;
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
            .unwrap_or_else(|| panic!("expected BVP solver field {section}.{field} to exist"))
    }

    fn as_string(value: &Value) -> &str {
        match value {
            Value::String(value) => value,
            other => panic!("expected Value::String, got {other:?}"),
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

    fn prepare_solver_tab(app: &Rc<RefCell<CombustionApp>>) {
        app.borrow_mut().active_bvp_tab = BvpEditorTab::Solver;
    }

    fn load_problem_test_document() -> DocumentMap {
        let content = include_str!("fixtures/bvp_problem.fixture").to_string();
        let mut parser = DocumentParser::new(content);
        parser
            .parse_document()
            .expect("tracked BVP GUI fixture must parse as a valid document");
        parser
            .get_result()
            .cloned()
            .expect("parsed BVP GUI fixture must yield a document map")
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
        assert_eq!(
            as_string(first_value(&app.document, "solver_settings", "telemetry")),
            "off"
        );
        assert!(matches!(
            first_value(&app.document, "solver_settings", "loglevel"),
            Value::Optional(Some(value))
                if matches!(value.as_ref(), Value::String(level) if level == "off")
        ));
        assert!(matches!(
            first_value(&app.document, "solver_settings", "linear_sys_method"),
            Value::Optional(None)
        ));
        assert!(
            !app.document
                .get("solver_settings")
                .is_some_and(|section| section.contains_key("dont_save_log"))
        );
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
        app.active_bvp_tab = BvpEditorTab::Solver;
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
            "Solver telemetry",
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
    fn bvp_density_model_control_roundtrips_into_the_physics_parser() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        assert_eq!(
            as_string(first_value(
                &app.borrow().document,
                "process_conditions",
                "density_model"
            )),
            "ideal_gas_fixed_composition"
        );

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(900.0, 760.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Density model")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Density model")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Constant density")
            .click_accesskit();
        harness.run();

        let document = app.borrow().document.clone();
        assert_eq!(
            as_string(first_value(
                &document,
                "process_conditions",
                "density_model"
            )),
            "constant"
        );
        let mut reactor = SimpleReactorTask::new();
        reactor
            .set_reactor_params_from_hashmap(&document)
            .expect("the GUI-edited task should pass the production physics parser");
        assert_eq!(reactor.density_model, DensityModel::Constant);
    }

    #[test]
    fn correction_newton_family_roundtrips_solver_options_and_aot_request() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "correction_newton")
            .click();
        harness.run();

        for label in [
            "Fixed-mesh CorrectionNewton route",
            "Symbolic assembly",
            "Execution",
            "Matrix layout",
            "Absolute tolerance",
            "Default relative tolerance",
            "Maximum Newton iterations",
            "Maximum Jacobian age",
            "Maximum damping trials",
            "Solver telemetry",
            "Logging level",
        ] {
            harness.get_by_label(label);
        }

        let parse_settings = || {
            let app_ref = app.borrow();
            let normalized = normalize_reactor_physics_task_map(&app_ref.document);
            let mut solver_document = DocumentMap::new();
            for section_name in ["solver_settings", "bounds", "rel_tolerance"] {
                if let Some(section) = normalized.get(section_name) {
                    solver_document.insert(section_name.to_string(), section.clone());
                }
            }
            parse_hybrid_bvp_solver_settings_from_document(&solver_document)
        };

        let settings = parse_settings().expect("GUI CorrectionNewton settings should parse");
        let BvpSolverSettings::CorrectionNewton(options) = settings else {
            panic!("GUI should select the RST CorrectionNewton route");
        };
        assert_eq!(options.matrix_layout, CorrectionNewtonMatrixLayout::Sparse);
        assert_eq!(options.absolute_tolerance, 1e-5);
        assert_eq!(options.default_relative_tolerance, 1e-5);
        assert_eq!(options.max_iterations, 100);
        assert_eq!(options.max_jacobian_age, 5);
        assert_eq!(options.max_damping_trials, 7);
        assert!(options.bounds.contains_key("Teta"));
        assert!(options.relative_tolerances.contains_key("Teta"));
        drop(options);

        harness
            .get_by_role_and_label(Role::ComboBox, "Symbolic assembly")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "atom_view_native")
            .click();
        harness.run();
        assert!(matches!(
            parse_settings().expect("AtomView CorrectionNewton settings should parse"),
            BvpSolverSettings::CorrectionNewton(options)
                if options.symbolic_frontend == CorrectionNewtonSymbolicFrontend::AtomView
        ));

        harness
            .get_by_role_and_label(Role::ComboBox, "Execution")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Execution")
            .click();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "aot").click();
        harness.run();
        for label in [
            "AOT codegen backend",
            "AOT build policy",
            "AOT build profile",
        ] {
            harness.get_by_label(label);
        }
        assert!(matches!(
            parse_settings(),
            Ok(BvpSolverSettings::CorrectionNewton(options))
                if options.symbolic_execution == CorrectionNewtonSymbolicExecution::Aot
        ));
        let request = aot_toolchain_request(&app.borrow().document)
            .expect("CorrectionNewton AOT should request explicit toolchain confirmation");
        assert_eq!(request.generated_backend, "correction_newton_aot");
    }

    #[test]
    fn bvp_solver_telemetry_choice_is_serialized_without_enabling_logging() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver telemetry")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver telemetry")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "counters")
            .click();
        harness.run();

        let app_ref = app.borrow();
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "telemetry"
            )),
            "counters"
        );
        assert!(
            !app_ref
                .document
                .get("solver_settings")
                .is_some_and(|section| section.contains_key("dont_save_log"))
        );

        let mut normalized_document = normalize_reactor_physics_task_map(&app_ref.document);
        let mut solver_document = DocumentMap::new();
        for section_name in [
            "solver_settings",
            "strategy_params",
            "adaptive_strategy",
            "grid_refinement",
            "bounds",
            "rel_tolerance",
        ] {
            if let Some(section) = normalized_document.remove(section_name) {
                solver_document.insert(section_name.to_string(), section);
            }
        }
        let settings = parse_hybrid_bvp_solver_settings_from_document(&solver_document)
            .expect("GUI-produced solver settings should satisfy the RST hybrid contract");
        let BvpSolverSettings::Damped(options) = settings else {
            panic!("GUI default should remain on the Damped solver family");
        };
        assert_eq!(
            options.generated_backend_config.bvp_telemetry_mode,
            BvpTelemetryMode::Counters
        );
    }

    #[test]
    fn bvp_damped_aot_timings_are_explained_and_rejected_by_rst() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        {
            let mut app_ref = app.borrow_mut();
            let settings = app_ref
                .document
                .get_mut("solver_settings")
                .expect("solver_settings section exists");
            settings.insert(
                "generated_backend".into(),
                Some(vec![Value::String("banded_aot_tcc".into())]),
            );
            settings.insert(
                "telemetry".into(),
                Some(vec![Value::String("timings".into())]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness.get_by_label(
            "Timings telemetry is unavailable for the Damped AOT route; select Counters or Detailed.",
        );
        for (label, expected_help) in [
            ("Residuals and Jacobian Backend", "High-level solver preset"),
            ("AOT codegen backend", "AOT code generation backend."),
            ("C compiler", "C compiler choice"),
            ("AOT build policy", "AOT build lifecycle policy."),
            ("AOT build profile", "AOT build profile."),
            ("AOT compile preset", "AOT compile preset."),
            ("AOT execution policy", "Runtime execution policy"),
            (
                "Residual target chunks",
                "Optional positive target chunk count for Damped AOT residual generation.",
            ),
            (
                "Jacobian target chunks",
                "Optional positive target chunk count for Damped AOT Jacobian generation.",
            ),
            ("Symbolic backend", "Symbolic assembly backend."),
            ("Banded linear solver", "Banded linear solver policy."),
            ("refinement_steps", "Iterative refinement steps"),
            (
                "Assembly parallelism",
                "Auto selects sequential or parallel assembly by system size.",
            ),
        ] {
            app.borrow_mut().current_help_info.clear();
            harness.get_by_label(label).scroll_to_me();
            harness.run();
            harness.run();
            if matches!(label, "AOT compile preset" | "AOT execution policy") {
                harness.get_by_role_and_label(Role::ComboBox, label).hover();
            } else {
                harness.get_by_label(label).hover();
            }
            harness.run();
            harness.run();
            let actual_help = app.borrow().current_help_info.clone();
            assert!(
                actual_help.contains(expected_help),
                "hover help for Damped AOT control `{label}` should contain `{expected_help}`, got `{}`",
                actual_help
            );
        }

        let normalized = normalize_reactor_physics_task_map(&app.borrow().document);
        let mut solver_document = DocumentMap::new();
        for section_name in [
            "solver_settings",
            "strategy_params",
            "adaptive_strategy",
            "grid_refinement",
            "bounds",
            "rel_tolerance",
        ] {
            if let Some(section) = normalized.get(section_name) {
                solver_document.insert(section_name.to_string(), section.clone());
            }
        }
        assert!(parse_hybrid_bvp_solver_settings_from_document(&solver_document).is_err());
    }

    #[test]
    fn bvp_damped_and_frozen_controls_render_content_aware_hover_help() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        for (label, expected_help) in [
            ("Solver family", "BVP solver family:"),
            ("scheme", "Spatial discretization scheme."),
            ("Linear algebra backend", "Linear algebra backend."),
            ("Nonlinear strategy", "Convergence strategy."),
            (
                "Linear system method override",
                "Optional linear-system override.",
            ),
            ("abs_tolerance", "Absolute convergence tolerance."),
            ("max_iterations", "Maximum solver iterations."),
            ("Solver telemetry", "Internal solver telemetry."),
            ("loglevel", "Common BVP text-logging level."),
            ("Residuals and Jacobian Backend", "High-level solver preset"),
            ("Symbolic backend", "Symbolic assembly backend."),
            ("Banded linear solver", "Banded linear solver policy."),
            ("refinement_steps", "Iterative refinement steps"),
        ] {
            app.borrow_mut().current_help_info.clear();
            harness.get_by_label(label).scroll_to_me();
            harness.run();
            harness.get_by_label(label).hover();
            harness.run();
            assert!(
                app.borrow().current_help_info.contains(expected_help),
                "hover help for Damped control `{label}` should contain `{expected_help}`"
            );
        }

        harness
            .get_by_role_and_label(Role::Button, "Physics")
            .click();
        harness.run();
        for (label, expected_help) in [
            (
                "process_conditions",
                "Main physical and chemical parameters",
            ),
            ("Density model", "rho(T) = rho0*T0/T"),
            (
                "boundary_condition",
                "Physical inlet mass fractions and temperature [K]",
            ),
            ("diffusion_coefficients", "Molecular diffusion coefficients"),
            (
                "Reactions",
                "Chemical reactions with Arrhenius kinetic parameters",
            ),
            ("bounds", "Variable bounds for solver convergence"),
        ] {
            app.borrow_mut().current_help_info.clear();
            harness.get_by_label(label).scroll_to_me();
            harness.run();
            harness.get_by_label(label).hover();
            harness.run();
            let actual_help = app.borrow().current_help_info.clone();
            assert!(
                actual_help.contains(expected_help),
                "hover help for BVP control `{label}` should contain `{expected_help}`, got `{actual_help}`"
            );
        }

        harness.get_by_label("P").scroll_to_me();
        harness.run();
        harness.get_by_label("P").hover();
        harness.run();
        assert!(
            app.borrow()
                .current_help_info
                .contains("System pressure [Pa]."),
            "hovering the physical pressure field should show its units and meaning; got `{}`",
            app.borrow().current_help_info
        );

        harness.get_by_label("T").scroll_to_me();
        harness.run();
        harness.get_by_label("T").hover();
        harness.run();
        assert!(
            app.borrow()
                .current_help_info
                .contains("Initial/boundary temperature [K]."),
            "hovering the temperature boundary field should show its units and meaning"
        );

        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .scroll_to_me();
        harness.run();
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
            .click();
        harness.run();
        for (label, expected_help) in [
            ("Solver family", "BVP solver family:"),
            ("Nonlinear strategy", "Convergence strategy."),
            ("tolerance", "Route-specific convergence tolerance."),
            ("max_iterations", "Maximum solver iterations."),
            ("Solver telemetry", "Internal solver telemetry."),
            ("loglevel", "Common BVP text-logging level."),
            ("Residuals and Jacobian Backend", "High-level solver preset"),
            ("Symbolic backend", "Symbolic assembly backend."),
            ("Banded linear solver", "Banded linear solver policy."),
            ("refinement_steps", "Iterative refinement steps"),
        ] {
            app.borrow_mut().current_help_info.clear();
            harness.get_by_label(label).scroll_to_me();
            harness.run();
            harness.get_by_label(label).hover();
            harness.run();
            assert!(
                app.borrow().current_help_info.contains(expected_help),
                "hover help for Frozen control `{label}` should contain `{expected_help}`"
            );
        }
    }

    #[test]
    fn bvp_sci_family_controls_round_trip_through_rst_settings_parser() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .click();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "sci").click();
        harness.run();
        for label in [
            "BVP_sci collocation solver",
            "Execution",
            "Symbolic assembly",
            "Matrix layout",
            "Maximum mesh nodes",
            "Solver telemetry",
        ] {
            harness.get_by_label(label);
        }
        for (label, expected_help) in [
            ("Solver family", "BVP solver family:"),
            ("Execution", "callback execution route."),
            ("Symbolic assembly", "Symbolic expression frontend."),
            ("Matrix layout", "Linear-system matrix layout."),
            ("Tolerance", "Route-specific convergence tolerance."),
            ("Boundary tolerance", "boundary residual tolerance."),
            ("Maximum mesh nodes", "Maximum collocation mesh nodes"),
            ("Maximum Newton iterations", "Maximum Newton iterations"),
            ("Maximum Jacobian refreshes", "Maximum Jacobian refreshes"),
            (
                "Maximum mesh refinements",
                "Maximum adaptive mesh refinements",
            ),
            ("Maximum backtracking steps", "Maximum backtracking steps"),
            ("Solver telemetry", "Internal solver telemetry."),
            ("Logging level", "Common BVP text-logging level."),
        ] {
            app.borrow_mut().current_help_info.clear();
            harness.get_by_label(label).scroll_to_me();
            harness.run();
            harness.get_by_label(label).hover();
            harness.run();
            assert!(
                app.borrow().current_help_info.contains(expected_help),
                "hover help for `{label}` should contain `{expected_help}`"
            );
        }

        harness
            .get_by_role_and_label(Role::ComboBox, "Matrix layout")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Matrix layout")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "banded")
            .click();
        harness.run();
        for (label, expected_help) in [
            ("Lower bandwidth", "Lower bandwidth for BVP_sci"),
            ("Upper bandwidth", "Upper bandwidth for BVP_sci"),
        ] {
            app.borrow_mut().current_help_info.clear();
            harness.get_by_label(label).scroll_to_me();
            harness.run();
            harness.get_by_label(label).hover();
            harness.run();
            assert!(
                app.borrow().current_help_info.contains(expected_help),
                "hover help for conditional Sci control `{label}` should contain `{expected_help}`"
            );
        }

        let parse_gui_settings = || {
            let app_ref = app.borrow();
            let normalized = normalize_reactor_physics_task_map(&app_ref.document);
            let mut solver_document = DocumentMap::new();
            for section_name in [
                "solver_settings",
                "strategy_params",
                "adaptive_strategy",
                "grid_refinement",
                "bounds",
                "rel_tolerance",
            ] {
                if let Some(section) = normalized.get(section_name) {
                    solver_document.insert(section_name.to_string(), section.clone());
                }
            }
            parse_hybrid_bvp_solver_settings_from_document(&solver_document)
        };

        assert!(matches!(
            parse_gui_settings().expect("GUI-produced Sci settings should parse"),
            BvpSolverSettings::Sci { .. }
        ));

        harness
            .get_by_role_and_label(Role::ComboBox, "Execution")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Execution")
            .click();
        harness.run();
        harness.get_by_role_and_label(Role::Button, "aot").click();
        harness.run();
        harness.get_by_label("AOT build policy");
        for (label, expected_help) in [
            ("AOT build policy", "AOT build lifecycle policy."),
            ("AOT build profile", "AOT build profile."),
            ("AOT codegen backend", "AOT code generation backend."),
            ("C compiler", "C compiler choice"),
        ] {
            app.borrow_mut().current_help_info.clear();
            harness.get_by_label(label).scroll_to_me();
            harness.run();
            harness.get_by_label(label).hover();
            harness.run();
            assert!(
                app.borrow().current_help_info.contains(expected_help),
                "hover help for `{label}` should contain `{expected_help}`"
            );
        }
        assert!(matches!(
            parse_gui_settings().expect("Sci AOT settings should parse with telemetry Off"),
            BvpSolverSettings::Sci { .. }
        ));
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver telemetry")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "timings")
            .click();
        harness.run();
        assert!(matches!(
            parse_gui_settings().expect("Sci AOT Timings telemetry should parse"),
            BvpSolverSettings::Sci { .. }
        ));

        harness
            .get_by_role_and_label(Role::ComboBox, "Execution")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Execution")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "lambdify")
            .click();
        harness.run();
        assert!(matches!(
            parse_gui_settings().expect("Sci Lambdify settings should parse"),
            BvpSolverSettings::Sci { .. }
        ));
        assert!(
            !app.borrow()
                .document
                .get("solver_settings")
                .is_some_and(|section| section.contains_key("aot_build_policy"))
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
            .click();
        harness.run();
        harness.get_by_label("Fixed-mesh Frozen/Naive route");
        let app_ref = app.borrow();
        assert_eq!(
            as_string(first_value(
                &app_ref.document,
                "solver_settings",
                "solver_family"
            )),
            "frozen"
        );
        assert!(!app_ref.document.contains_key("bounds"));
        drop(app_ref);
        assert!(matches!(
            parse_gui_settings().expect("GUI-produced Frozen settings should parse"),
            BvpSolverSettings::Frozen(_)
        ));
    }

    #[test]
    fn bvp_solver_backend_keeps_linear_system_override_typed() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
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

        harness.get_by_role_and_label(Role::Button, "Physics");
        for label in [
            "Legacy raw document",
            "process_conditions",
            "boundary_condition",
        ] {
            harness.get_by_label(label);
        }
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
        harness.run();
        for label in [
            "Solver backend",
            "Initial guess",
            "initial_guess",
            "Postprocessing",
            "postprocessing",
        ] {
            harness.get_by_label(label);
        }
    }

    #[test]
    fn solver_tab_persists_kithe_symbolic_assembly_policy() {
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
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Assembly parallelism")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Assembly parallelism")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Parallel")
            .click();
        harness.run();

        assert_eq!(
            as_string(first_value(
                &app.borrow().document,
                "kithe_settings",
                "symbolic_rhs_assembly_policy"
            )),
            "parallel"
        );
        assert!(
            app.borrow()
                .document_to_string()
                .contains("symbolic_rhs_assembly_policy: parallel")
        );
        let serialized = app.borrow().document_to_string();
        let mut parser = DocumentParser::new(serialized);
        parser
            .parse_document()
            .expect("GUI settings document should round-trip through the task parser");
        assert_eq!(
            as_string(first_value(
                parser.get_result().expect("parsed map should exist"),
                "kithe_settings",
                "symbolic_rhs_assembly_policy"
            )),
            "parallel"
        );
    }

    #[test]
    fn bvp_screen_renders_advanced_solver_sections_separately_from_physics() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        let mut open = true;
        let mut harness = Harness::new_ui(move |ui| {
            app.show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
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
        harness
            .get_by_role_and_label(Role::Button, "Solver")
            .click();
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
        prepare_solver_tab(&app);
        app.borrow_mut().document = load_problem_test_document();

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        assert!(!app.borrow().document.contains_key("adaptive_strategy"));
        assert!(!app.borrow().document.contains_key("grid_refinement"));

        let adaptive_control = harness.get_by_label("Adaptive grid refinement");
        adaptive_control.scroll_to_me();
        harness.run();
        harness
            .get_by_label("Adaptive grid refinement")
            .click_accesskit();
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
        prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
            prepare_solver_tab(&app);
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
            if case.expected_codegen == AotCodegenBackend::C {
                assert_eq!(
                    options.generated_backend_config.aot_c_compiler.as_deref(),
                    case.expected_compiler,
                    "wrong C compiler for {}",
                    case.preset
                );
            } else {
                assert!(
                    !app_ref
                        .document
                        .get("solver_settings")
                        .is_some_and(|section| section.contains_key("aot_c_compiler")),
                    "non-C backend {} must not serialize a C compiler",
                    case.preset
                );
            }

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
            prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
        assert!(
            harness.query_all_by_label("HMX").count() >= 2,
            "HMX should be visible in both the physics and atomic-composition frames"
        );
        assert!(
            harness.query_all_by_label("DecompositionProduct").count() >= 2,
            "the species should appear in both the physics and atomic-composition frames"
        );
    }

    #[test]
    fn bvp_species_editor_adds_and_cascades_species_data() {
        let mut document = load_problem_test_document();
        add_bvp_species(&mut document, "NewSpecies").expect("species should be addable");
        assert!(
            document["process_conditions"]["substances"]
                .as_ref()
                .unwrap()
                .iter()
                .any(|value| value.as_string().is_some_and(|name| name == "NewSpecies"))
        );
        assert!(document["boundary_condition"].contains_key("NewSpecies"));
        assert!(document["diffusion_coefficients"].contains_key("NewSpecies"));
        assert!(document.contains_key("NewSpecies"));
        assert!(add_bvp_species(&mut document, "NewSpecies").is_err());

        document
            .get_mut("reactions")
            .unwrap()
            .insert("NewSpecies=>HMX".to_string(), Some(vec![Value::Float(1.0)]));
        assert!(apply_document_deletion(
            &mut document,
            &DocumentDeleteRequest::Species {
                species: "NewSpecies".to_string()
            }
        ));
        assert!(!document.contains_key("NewSpecies"));
        assert!(!document["boundary_condition"].contains_key("NewSpecies"));
        assert!(!document["diffusion_coefficients"].contains_key("NewSpecies"));
        assert!(!document["reactions"].contains_key("NewSpecies=>HMX"));
    }

    #[test]
    fn bvp_reaction_editor_keeps_distinct_reactions_and_typed_parameters() {
        let mut document = load_problem_test_document();
        add_bvp_reaction(&mut document, "A=>B").expect("first reaction should be addable");
        add_bvp_reaction(&mut document, "B=>C").expect("second reaction should be addable");
        assert!(add_bvp_reaction(&mut document, "A=>B").is_err());

        let reactions = &document["reactions"];
        assert_eq!(reactions.len(), 3);
        assert!(matches!(
            reactions["A=>B"].as_ref().unwrap().as_slice(),
            [Value::Vector(parameters)] if parameters == &[1.0, 0.0, 0.0, 0.0]
        ));
    }

    #[test]
    fn bvp_physics_uses_collection_editors_instead_of_generic_add_field() {
        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.document = load_problem_test_document();
        let mut open = true;
        let mut harness = Harness::builder()
            .with_accessibility_check(false)
            .build_ui(move |ui| app.show(ui.ctx(), &mut open));
        harness.run();

        harness.get_by_label("Substances");
        harness.get_by_label("Add substance");
        harness.get_by_label("Reactions");
        harness.get_by_label("Add reaction");
        assert!(harness.query_all_by_label("Add element").count() >= 2);
        assert!(
            harness.query_all_by_label("Add Field").next().is_none(),
            "structured physics collections must not use the obsolete generic field editor"
        );
    }

    #[test]
    fn bvp_window_has_no_problem_type_selector_and_physics_sections_are_visible_frames() {
        let mut app = CombustionApp::new();
        assert_eq!(app.selected_problem, ProblemsEnum::BVPSimple);
        let mut open = true;
        let mut harness = Harness::new_ui(move |ui| {
            app.show(ui.ctx(), &mut open);
        });

        harness.run();

        assert!(
            harness.query_all_by_label("Problem Type").next().is_none(),
            "the dedicated BVP window must not offer a BVP/None problem selector"
        );
        assert!(
            harness
                .query_all_by_label("Create Section")
                .next()
                .is_none(),
            "the structured BVP editor must not expose generic section creation"
        );
        assert!(
            harness
                .query_all_by_label("New section name:")
                .next()
                .is_none(),
            "the structured BVP editor must not expose a raw section-name field"
        );
        harness.get_by_role_and_label(Role::Button, "Physics");
        assert!(
            harness
                .query_all_by_label("process_conditions")
                .next()
                .is_some()
        );
        harness.get_by_role_and_label(Role::SpinButton, "P");
        assert!(
            harness
                .query_all_by_role_and_label(Role::Button, "process_conditions")
                .next()
                .is_none(),
            "physics sections should render as framed content, not collapsed menus"
        );
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
        prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
        prepare_solver_tab(&app);
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
            .scroll_to_me();
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
        prepare_solver_tab(&app);
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
    fn duplicate_fields_preserve_existing_document_data() {
        let mut section = HashMap::new();
        insert_field_if_absent(&mut section, "temperature".to_string(), Value::Float(900.0))
            .expect("first field insertion should succeed");

        let duplicate_field = insert_field_if_absent(
            &mut section,
            "temperature".to_string(),
            Value::Float(1200.0),
        );
        assert!(duplicate_field.is_err());
        assert_eq!(
            section
                .get("temperature")
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(&Value::Float(900.0))
        );

        assert!(section.contains_key("temperature"));
    }

    #[test]
    fn replacing_problem_clears_file_binding_without_overwriting_previous_task() {
        let temp_dir = tempdir().expect("temp dir should be available");
        let path = temp_dir.path().join("existing_bvp_task.txt");
        std::fs::write(&path, "original task contents").expect("test task should be writable");

        let mut app = CombustionApp::new_with_problem(ProblemsEnum::BVPSimple);
        app.current_file_path = Some(path.clone());
        app.new_field_names
            .insert("physics".to_string(), "pending field".to_string());
        app.last_run_message = Some("stale result".to_string());
        app.last_run_is_error = true;

        app.replace_problem(ProblemsEnum::None);

        assert_eq!(app.selected_problem, ProblemsEnum::None);
        assert!(app.document.is_empty());
        assert!(app.current_file_path.is_none());
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

    #[test]
    fn solver_route_specific_controls_roundtrip_into_the_hybrid_api() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(900.0, 760.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));
        harness.run();

        macro_rules! choose {
            ($label:literal, $choice:literal) => {{
                harness
                    .get_by_role_and_label(Role::ComboBox, $label)
                    .scroll_to_me();
                harness.run();
                harness
                    .get_by_role_and_label(Role::ComboBox, $label)
                    .click();
                harness.run();
                harness
                    .get_by_role_and_label(Role::Button, $choice)
                    .click_accesskit();
                harness.run();
            }};
        }

        choose!("Solver family", "sci");
        {
            let mut app = app.borrow_mut();
            let settings = app
                .document
                .get_mut("solver_settings")
                .expect("solver settings");
            settings.insert(
                "matrix_layout".into(),
                Some(vec![Value::String("sparse".into())]),
            );
            settings.insert(
                "execution_policy".into(),
                Some(vec![Value::String("parallel".into())]),
            );
            settings.insert("execution_min_work".into(), Some(vec![Value::Usize(128)]));
            settings.insert(
                "output_policy".into(),
                Some(vec![Value::String("final_only".into())]),
            );
            settings.insert(
                "sparse_factorization".into(),
                Some(vec![Value::String("superlu".into())]),
            );
        }
        harness.run();

        let document = app.borrow().document.clone();
        let settings = document.get("solver_settings").expect("solver settings");
        assert_eq!(
            settings
                .get("execution_min_work")
                .and_then(Option::as_ref)
                .and_then(|values| values.first()),
            Some(&Value::Usize(128))
        );

        let normalized = normalize_reactor_physics_task_map(&document);
        let solver_document = normalized
            .into_iter()
            .filter(|(section, _)| {
                [
                    "solver_settings",
                    "strategy_params",
                    "adaptive_strategy",
                    "grid_refinement",
                    "bounds",
                    "rel_tolerance",
                ]
                .contains(&section.as_str())
            })
            .collect::<DocumentMap>();
        let parsed = parse_hybrid_bvp_solver_settings_from_document(&solver_document)
            .expect("GUI-selected BVP_sci settings should parse through the production contract");
        let BvpSolverSettings::Sci {
            options,
            sparse_factorization,
            ..
        } = parsed
        else {
            panic!("GUI selected the BVP_sci route");
        };
        assert_eq!(
            options.execution_policy,
            RustedSciThe::numerical::BVP_sci::new::BvpSciExecutionPolicy::Parallel {
                min_work: 128
            }
        );
        assert_eq!(
            options.output_policy,
            RustedSciThe::numerical::BVP_sci::new::BvpSciOutputPolicy::FinalOnly
        );
        assert_eq!(
            sparse_factorization,
            RustedSciThe::numerical::BVP_api::BvpSciSparseFactorization::SuperLu
        );
    }

    #[test]
    fn bvp_sci_aot_artifact_and_chunk_controls_roundtrip_into_hybrid_api() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        {
            let mut app = app.borrow_mut();
            set_bvp_solver_family(&mut app.document, BvpSolverFamily::Sci);
            let settings = app
                .document
                .get_mut("solver_settings")
                .expect("solver settings");
            for (field, value) in [
                ("execution", "aot"),
                ("matrix_layout", "banded_auto"),
                ("assembly", "atom_view_native"),
                ("aot_build_policy", "build_if_missing"),
                ("aot_build_profile", "release"),
                ("aot_codegen_backend", "rust"),
                ("aot_output_dir", "artifacts/bvp-sci"),
                ("aot_handoff_path", "reports/bvp-handoff.json"),
                ("aot_crate_name", "bvp_generated"),
                ("aot_module_name", "generated_bvp"),
                ("aot_residual_chunking", "target_chunks"),
                ("aot_sparse_jacobian_chunking", "rows"),
                ("aot_banded_jacobian_chunking", "diagonals"),
            ] {
                settings.insert(
                    field.to_string(),
                    Some(vec![Value::String(value.to_string())]),
                );
            }
            for (field, value) in [
                ("aot_residual_chunk_size", 3),
                ("aot_sparse_jacobian_chunk_size", 5),
                ("aot_banded_jacobian_chunk_size", 2),
            ] {
                settings.insert(field.to_string(), Some(vec![Value::Integer(value)]));
            }
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(900.0, 760.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));
        harness.run();

        for label in [
            "AOT build policy",
            "AOT build profile",
            "AOT codegen backend",
            "AOT output directory",
            "AOT handoff path",
            "Generated crate name",
            "Generated module name",
            "Residual chunking",
            "Residual chunk size",
            "Sparse Jacobian chunking",
            "Sparse Jacobian chunk size",
            "Banded Jacobian chunking",
            "Banded Jacobian chunk size",
        ] {
            harness.get_by_label(label);
        }

        let document = app.borrow().document.clone();
        let normalized = normalize_reactor_physics_task_map(&document);
        let solver_document = normalized
            .into_iter()
            .filter(|(section, _)| {
                [
                    "solver_settings",
                    "strategy_params",
                    "adaptive_strategy",
                    "grid_refinement",
                    "bounds",
                    "rel_tolerance",
                ]
                .contains(&section.as_str())
            })
            .collect::<DocumentMap>();
        let BvpSolverSettings::Sci { aot_config, .. } =
            parse_hybrid_bvp_solver_settings_from_document(&solver_document)
                .expect("AOT GUI settings should satisfy the production hybrid parser")
        else {
            panic!("the GUI selected the BVP_sci route");
        };
        let aot = aot_config.expect("AOT execution should produce a generated config");
        assert_eq!(
            aot.output_parent_dir,
            Some(std::path::PathBuf::from("artifacts/bvp-sci"))
        );
        assert_eq!(
            aot.handoff_path,
            Some(std::path::PathBuf::from("reports/bvp-handoff.json"))
        );
        assert_eq!(aot.crate_name_override.as_deref(), Some("bvp_generated"));
        assert_eq!(aot.module_name_override.as_deref(), Some("generated_bvp"));
        assert_eq!(
            aot.residual_chunking_strategy,
            RustedSciThe::symbolic::codegen::codegen_runtime_api::ResidualChunkingStrategy::ByTargetChunkCount {
                target_chunks: 3,
            }
        );
        assert_eq!(
            aot.sparse_jacobian_chunking_strategy,
            RustedSciThe::symbolic::codegen::codegen_tasks::SparseChunkingStrategy::ByRowCount {
                rows_per_chunk: 5,
            }
        );
        assert_eq!(
            aot.banded_jacobian_chunking_strategy,
            RustedSciThe::symbolic::codegen::codegen_tasks::BandedChunkingStrategy::ByDiagonalCount {
                diagonals_per_chunk: 2,
            }
        );

        {
            let mut app = app.borrow_mut();
            let settings = app
                .document
                .get_mut("solver_settings")
                .expect("solver settings");
            settings.insert(
                "aot_build_policy".to_string(),
                Some(vec![Value::String("require_prebuilt".to_string())]),
            );
            settings.remove("aot_build_profile");
        }
        harness.run();
        assert!(
            !app.borrow()
                .document
                .get("solver_settings")
                .expect("solver settings")
                .contains_key("aot_build_profile"),
            "profile should be omitted for policies that do not compile"
        );
        let document = app.borrow().document.clone();
        let normalized = normalize_reactor_physics_task_map(&document);
        let solver_document = normalized
            .into_iter()
            .filter(|(section, _)| {
                [
                    "solver_settings",
                    "strategy_params",
                    "adaptive_strategy",
                    "grid_refinement",
                    "bounds",
                    "rel_tolerance",
                ]
                .contains(&section.as_str())
            })
            .collect::<DocumentMap>();
        assert!(parse_hybrid_bvp_solver_settings_from_document(&solver_document).is_ok());
    }

    #[test]
    fn bvp_damped_aot_chunk_targets_roundtrip_into_hybrid_api() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        {
            let mut app = app.borrow_mut();
            let settings = app
                .document
                .get_mut("solver_settings")
                .expect("solver settings");
            settings.insert(
                "generated_backend".into(),
                Some(vec![Value::String("banded_aot_tcc".into())]),
            );
            settings.insert(
                "backend_policy".into(),
                Some(vec![Value::String("aot_only".into())]),
            );
        }

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(900.0, 760.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));
        harness.run();

        let normalized = normalize_reactor_physics_task_map(&app.borrow().document);
        let solver_document = normalized
            .into_iter()
            .filter(|(section, _)| {
                [
                    "solver_settings",
                    "strategy_params",
                    "adaptive_strategy",
                    "grid_refinement",
                    "bounds",
                    "rel_tolerance",
                ]
                .contains(&section.as_str())
            })
            .collect::<DocumentMap>();
        let settings = parse_hybrid_bvp_solver_settings_from_document(&solver_document)
            .expect("unset chunk targets should preserve RST defaults");
        let BvpSolverSettings::Damped(options) = settings else {
            panic!("the GUI selected the Damped family");
        };
        assert_eq!(
            options
                .generated_backend_config
                .aot_chunking_policy
                .residual,
            None
        );
        assert_eq!(
            options
                .generated_backend_config
                .aot_chunking_policy
                .sparse_jacobian,
            None
        );

        {
            let mut app = app.borrow_mut();
            let settings = app
                .document
                .get_mut("solver_settings")
                .expect("solver settings");
            settings.insert(
                "aot_residual_target_chunks".into(),
                Some(vec![Value::Usize(5)]),
            );
            settings.insert(
                "aot_jacobian_target_chunks".into(),
                Some(vec![Value::Usize(3)]),
            );
        }
        harness.run();
        harness.get_by_role_and_label(Role::SpinButton, "Residual target chunks");
        harness.get_by_role_and_label(Role::SpinButton, "Jacobian target chunks");

        let normalized = normalize_reactor_physics_task_map(&app.borrow().document);
        let solver_document = normalized
            .into_iter()
            .filter(|(section, _)| {
                [
                    "solver_settings",
                    "strategy_params",
                    "adaptive_strategy",
                    "grid_refinement",
                    "bounds",
                    "rel_tolerance",
                ]
                .contains(&section.as_str())
            })
            .collect::<DocumentMap>();
        let settings = parse_hybrid_bvp_solver_settings_from_document(&solver_document)
            .expect("Damped AOT chunk targets should satisfy the RST parser");
        let BvpSolverSettings::Damped(options) = settings else {
            panic!("the GUI selected the Damped family");
        };
        assert_eq!(
            options.generated_backend_config.aot_chunking_policy.residual,
            Some(
                RustedSciThe::symbolic::codegen::codegen_runtime_api::ResidualChunkingStrategy::ByTargetChunkCount {
                    target_chunks: 5,
                }
            )
        );
        assert_eq!(
            options
                .generated_backend_config
                .aot_chunking_policy
                .sparse_jacobian,
            Some(
                RustedSciThe::symbolic::codegen::codegen_tasks::SparseChunkingStrategy::ByTargetChunkCount {
                    target_chunks: 3,
                }
            )
        );

        {
            let mut app = app.borrow_mut();
            let settings = app
                .document
                .get_mut("solver_settings")
                .expect("solver settings");
            settings.insert(
                "generated_backend".into(),
                Some(vec![Value::String("banded_lambdify".into())]),
            );
            settings.insert(
                "backend_policy".into(),
                Some(vec![Value::String("lambdify_only".into())]),
            );
        }
        harness.run();
        {
            let app = app.borrow();
            let settings = app
                .document
                .get("solver_settings")
                .expect("solver settings");
            assert!(!settings.contains_key("aot_residual_target_chunks"));
            assert!(!settings.contains_key("aot_jacobian_target_chunks"));
        }
        assert!(
            harness
                .query_all_by_role_and_label(Role::SpinButton, "Residual target chunks")
                .next()
                .is_none()
        );
        assert!(
            harness
                .query_all_by_role_and_label(Role::SpinButton, "Jacobian target chunks")
                .next()
                .is_none()
        );
    }

    #[test]
    fn bvp_aot_chunk_target_counts_must_be_positive() {
        for field in ["aot_residual_target_chunks", "aot_jacobian_target_chunks"] {
            let mut document = DocumentMap::new();
            let mut settings = HashMap::new();
            settings.insert(field.to_string(), Some(vec![Value::Usize(0)]));
            document.insert("solver_settings".to_string(), settings);

            let error = normalize_bvp_usize_fields(&mut document)
                .expect_err("zero is not a valid target chunk count");
            assert!(error.contains(field), "{error}");
            assert!(error.contains("positive chunk count"), "{error}");
        }
    }

    #[test]
    fn damped_lambdify_scheduling_and_event_logging_are_editable() {
        let app = Rc::new(RefCell::new(CombustionApp::new_with_problem(
            ProblemsEnum::BVPSimple,
        )));
        prepare_solver_tab(&app);
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();

        for (label, value) in [
            ("Lambdify callback policy", "auto"),
            ("Backend event logging", "detailed"),
        ] {
            harness
                .get_by_role_and_label(Role::ComboBox, label)
                .scroll_to_me();
            harness.run();
            harness.get_by_role_and_label(Role::ComboBox, label).click();
            harness.run();
            harness
                .get_by_role_and_label(Role::Button, value)
                .click_accesskit();
            harness.run();
        }

        let document = app.borrow().document.clone();
        let solver_settings = document.get("solver_settings").expect("solver settings");
        assert_eq!(
            solver_settings
                .get("lambdify_execution_policy")
                .and_then(Option::as_ref)
                .and_then(|values| values.first()),
            Some(&Value::String("auto".into()))
        );
        assert_eq!(
            solver_settings
                .get("lambdify_min_work")
                .and_then(Option::as_ref)
                .and_then(|values| values.first()),
            Some(&Value::Usize(128))
        );
        assert_eq!(
            solver_settings
                .get("bvp_logging_mode")
                .and_then(Option::as_ref)
                .and_then(|values| values.first()),
            Some(&Value::String("detailed".into()))
        );
        let normalized = normalize_reactor_physics_task_map(&document);
        let solver_document = normalized
            .into_iter()
            .filter(|(section, _)| {
                [
                    "solver_settings",
                    "strategy_params",
                    "adaptive_strategy",
                    "grid_refinement",
                    "bounds",
                    "rel_tolerance",
                ]
                .contains(&section.as_str())
            })
            .collect::<DocumentMap>();
        assert!(parse_hybrid_bvp_solver_settings_from_document(&solver_document).is_ok());
    }
}
