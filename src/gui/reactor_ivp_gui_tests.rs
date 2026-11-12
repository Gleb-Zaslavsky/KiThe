//! Story-style regression tests for the condensed reactor-IVP GUI.
//!
//! Hypotheses:
//! - The screen should expose the reactor-IVP sections immediately on first render.
//! - Preview Task should update the visible status and store a preview snapshot.
//! - AOT-specific controls should appear when the execution mode is switched to AOT.
//! - Postprocessing controls should expose the export surface used by LSODE2.
//!
//! Expected result:
//! the reactor-IVP GUI should behave like a structured editor with typed solver
//! controls, not like a raw text form.

#[cfg(test)]
mod tests {
    use super::super::reactor_ivp_gui::{
        IvpGuiConfig, ReactorIvpApp, ReactorIvpGuiCallbackPolicy, ReactorIvpGuiExecutionMode,
        build_reactor_ivp_task_preview_snapshot,
    };
    use crate::Kinetics::mechfinder_api::{ReactionData, ReactionKinetics};
    use crate::SimpleReactorsIVP::solver_backend::{
        ReactorIvpCallbackPolicy, ReactorIvpExecutionBackend, ReactorIvpFamilyOptions,
    };
    use RustedSciThe::numerical::SolverTelemetryMode;
    use RustedSciThe::numerical::logging::SolverLogLevel;
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::cell::RefCell;
    use std::panic::{self, AssertUnwindSafe};
    use std::rc::Rc;

    #[test]
    fn reactor_ivp_screen_exposes_structured_sections_and_actions() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();

        for label in [
            "Physics",
            "Substances",
            "Reactions",
            "Add substance",
            "Add reaction",
            "Stop condition",
            "Initial state",
            "Solver backend",
            "Integration",
            "Postprocessing",
            "Save TXT",
            "Save CSV",
            "Write report",
            "Plotters PNG",
            "Gnuplot PNG",
            "Terminal plot",
        ] {
            harness.get_by_label(label);
        }

        harness.get_by_role_and_label(Role::Button, "Preview Task");
        harness.get_by_role_and_label(Role::Button, "Run Calculation");
        for label in [
            "Solver family",
            "Symbolic frontend",
            "LSODE2 controller",
            "Matrix backend",
            "Execution backend",
            "Telemetry",
            "Solver logging",
        ] {
            harness.get_by_role_and_label(Role::ComboBox, label);
        }
    }

    #[test]
    fn solver_family_controls_are_capability_specific_and_reachable_in_compact_window() {
        for (choice, expected_family) in [
            (
                "BDF (standalone)",
                crate::SimpleReactorsIVP::solver_backend::ReactorIvpSolverFamily::Bdf,
            ),
            (
                "Radau",
                crate::SimpleReactorsIVP::solver_backend::ReactorIvpSolverFamily::Radau,
            ),
        ] {
            let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
            let mut open = true;
            let app_for_ui = Rc::clone(&app);
            let mut harness = Harness::builder()
                .with_size(egui::Vec2::new(820.0, 620.0))
                .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));

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
                .get_by_role_and_label(Role::Button, choice)
                .click_accesskit();
            harness.run();

            let app_ref = app.borrow();
            assert_eq!(app_ref.config.solver_family, expected_family);
            assert_eq!(
                app_ref.config.matrix_backend,
                crate::SimpleReactorsIVP::solver_backend::ReactorIvpMatrixBackend::Dense
            );
            assert_eq!(
                app_ref.config.execution_mode,
                ReactorIvpGuiExecutionMode::Lambdify
            );
            assert!(!app_ref.config.stop_condition_enabled);
            drop(app_ref);
            harness.get_by_label("Dense (required by selected family)");
        }
    }

    #[test]
    fn ivp_solver_settings_are_accessible_and_update_the_typed_configuration() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(820.0, 620.0))
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

        choose!("Symbolic frontend", "ExprLegacy");
        choose!("LSODE2 controller", "Adams only");
        choose!("Matrix backend", "Banded");
        choose!("Telemetry", "Detailed");
        choose!("Solver logging", "Debug");
        choose!("Execution backend", "AOT");
        choose!("AOT toolchain", "Zig");
        choose!("AOT profile", "Debug");

        let config = app.borrow().config.clone();
        assert_eq!(
            config.symbolic_backend,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpSymbolicBackend::ExprLegacy
        );
        assert_eq!(
            config.method,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpMethod::Adams
        );
        assert_eq!(
            config.matrix_backend,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpMatrixBackend::Banded
        );
        assert_eq!(format!("{:?}", config.telemetry), "Detailed");
        assert_eq!(format!("{:?}", config.logging), "Debug");
        assert_eq!(format!("{:?}", config.aot_toolchain), "Zig");
        assert_eq!(format!("{:?}", config.aot_profile), "Debug");
        assert_eq!(config.execution_mode, ReactorIvpGuiExecutionMode::Aot);
    }

    #[test]
    fn solver_family_combo_opens_and_selects_with_keyboard() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(820.0, 620.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));

        harness.run();
        let family = harness.get_by_role_and_label(Role::ComboBox, "Solver family");
        family.focus();
        harness.run();
        assert!(
            harness
                .get_by_role_and_label(Role::ComboBox, "Solver family")
                .is_focused(),
            "AccessKit focus should be reflected by egui's widget response"
        );

        harness.key_press(egui::Key::Enter);
        harness.run();
        harness.get_by_label("Radau");

        harness.key_press(egui::Key::Escape);
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .focus();
        harness.run();
        harness.key_press(egui::Key::Space);
        harness.run();
        harness.get_by_label("Radau");
        harness.key_press(egui::Key::Escape);
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .focus();
        harness.run();
        harness.key_press(egui::Key::ArrowDown);
        harness.run();
        assert_eq!(
            app.borrow().config.solver_family,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpSolverFamily::Bdf
        );
    }

    #[test]
    fn all_solver_combos_support_keyboard_arrow_selection() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(820.0, 620.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));
        harness.run();

        macro_rules! press_on {
            ($label:literal, $key:expr) => {{
                harness
                    .get_by_role_and_label(Role::ComboBox, $label)
                    .scroll_to_me();
                harness.run();
                harness
                    .get_by_role_and_label(Role::ComboBox, $label)
                    .focus();
                harness.run();
                harness.key_press($key);
                harness.run();
            }};
        }

        press_on!("Solver family", egui::Key::ArrowDown);
        assert_eq!(
            app.borrow().config.solver_family,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpSolverFamily::Bdf
        );
        press_on!("Solver family", egui::Key::ArrowUp);
        press_on!("Symbolic frontend", egui::Key::ArrowDown);
        assert_eq!(
            app.borrow().config.symbolic_backend,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpSymbolicBackend::ExprLegacy
        );
        press_on!("LSODE2 controller", egui::Key::ArrowDown);
        press_on!("Matrix backend", egui::Key::ArrowDown);
        press_on!("Execution backend", egui::Key::ArrowDown);
        harness.run();
        press_on!("AOT toolchain", egui::Key::ArrowDown);
        press_on!("AOT profile", egui::Key::ArrowDown);
        assert_eq!(
            app.borrow().config.method,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpMethod::Bdf
        );
        assert_eq!(
            app.borrow().config.matrix_backend,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpMatrixBackend::Banded
        );
        assert_eq!(
            app.borrow().config.execution_mode,
            ReactorIvpGuiExecutionMode::Aot
        );
        assert_eq!(
            app.borrow().config.symbolic_backend,
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpSymbolicBackend::ExprLegacy
        );
        assert_eq!(format!("{:?}", app.borrow().config.aot_toolchain), "Zig");
        assert_eq!(format!("{:?}", app.borrow().config.aot_profile), "Debug");

        press_on!("Execution backend", egui::Key::ArrowUp);
        press_on!("Callback policy", egui::Key::ArrowDown);
        press_on!("Telemetry", egui::Key::ArrowDown);
        press_on!("Solver logging", egui::Key::ArrowDown);
        let config = app.borrow().config.clone();
        assert_eq!(
            config.callback_policy,
            ReactorIvpGuiCallbackPolicy::Parallel
        );
        assert_eq!(config.telemetry, SolverTelemetryMode::Counters);
        assert_eq!(config.logging, SolverLogLevel::Error);
    }

    #[test]
    fn callback_policy_and_standalone_iteration_budgets_reach_solver_config() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(820.0, 620.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));
        harness.run();

        harness
            .get_by_role_and_label(Role::ComboBox, "Callback policy")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Callback policy")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Parallel")
            .click_accesskit();
        harness.run();
        harness.get_by_label("Minimum callback work");

        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .scroll_to_me();
        harness.run();
        harness
            .get_by_role_and_label(Role::ComboBox, "Solver family")
            .click();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "BDF (standalone)")
            .click_accesskit();
        harness.run();
        harness.get_by_label("BDF maximum steps");

        let mut config = app.borrow().config.clone();
        config.bdf_max_steps = 12_345;
        config.callback_min_work = 64;
        let solver = config.to_solver_config();
        assert_eq!(
            solver.callback_policy,
            Some(ReactorIvpCallbackPolicy::Parallel { min_work: 64 })
        );
        assert!(matches!(
            solver.family_options,
            ReactorIvpFamilyOptions::Bdf { max_steps: 12_345 }
        ));

        config.solver_family =
            crate::SimpleReactorsIVP::solver_backend::ReactorIvpSolverFamily::Radau;
        config.radau_max_steps = 23_456;
        config.radau_max_newton_iterations = 9;
        assert!(matches!(
            config.to_solver_config().family_options,
            ReactorIvpFamilyOptions::Radau {
                max_steps: 23_456,
                max_newton_iterations: 9
            }
        ));
        assert_eq!(
            config.callback_policy,
            ReactorIvpGuiCallbackPolicy::Parallel
        );
    }

    #[test]
    fn preview_task_updates_status_and_snapshot() {
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

        let app_ref = app.borrow();
        assert_eq!(
            app_ref.last_status(),
            Some("Task preview printed to console.")
        );
        assert!(app_ref.last_error().is_none());
        assert!(app_ref.last_preview().is_some_and(
            |snapshot| !snapshot.summary_rows.is_empty() && !snapshot.equation_rows.is_empty()
        ));
    }

    #[test]
    fn ivp_mechanism_edits_apply_with_aligned_thermal_effects_and_species_state() {
        let app = ReactorIvpApp::new();
        let mut config = app.config.clone();
        config.substances.push("He".to_string());
        config.initial_state.species.insert("He".to_string(), 0.0);
        config.reactions.push(
            ReactionData::try_new_elementary("H=>H2".to_string(), vec![250.0, 0.0, 0.0], None)
                .unwrap(),
        );
        config.thermal_effects.push(-150_000.0);

        let mut task = app.task().clone();
        config
            .apply_to_task(&mut task)
            .expect("edited mechanism should apply transactionally");
        assert_eq!(task.kindata.substances(), &["H2", "H", "He"]);
        assert_eq!(task.kindata.equations().len(), 2);
        assert_eq!(task.thermal_effects, vec![-2.0e5, -150_000.0]);
        assert_eq!(
            task.kindata.stoichiometric_analyzer().stecheo_matrx[0].len(),
            3
        );

        config.reactions.pop();
        config.thermal_effects.pop();
        config.substances.pop();
        config.initial_state.species.remove("He");
        config
            .apply_to_task(&mut task)
            .expect("removing a reaction and unused substance should keep state aligned");
        assert_eq!(task.kindata.substances(), &["H2", "H"]);
        assert_eq!(task.kindata.equations().len(), 1);
        assert_eq!(task.thermal_effects.len(), 1);
    }

    #[test]
    fn ivp_preview_report_has_console_ready_chemistry_and_equations() {
        let app = ReactorIvpApp::new();
        let snapshot = build_reactor_ivp_task_preview_snapshot(&app.task, &app.config)
            .expect("demo IVP task should preview");
        let report = snapshot.format_console_report();
        assert!(report.contains("Reactor IVP preview:"));
        assert!(report.contains("chemistry"));
        assert!(report.contains("H2=>2H"));
        assert!(report.contains("Reactor IVP equations:"));
    }

    #[test]
    fn ivp_reaction_parameters_accept_scientific_notation() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();

        for (label, value, expected) in [("A", "1e12", 1.0e12), ("E (J/mol)", "2.5e6", 2.5e6)] {
            harness
                .query_all_by_label(label)
                .next()
                .unwrap_or_else(|| panic!("missing reaction parameter input `{label}`"))
                .scroll_to_me();
            harness.run();
            harness.query_all_by_label(label).next().unwrap().click();
            harness.run();
            harness.key_combination_modifiers(egui::Modifiers::CTRL, &[egui::Key::A]);
            harness
                .query_all_by_label(label)
                .next()
                .unwrap()
                .type_text(value);
            harness.key_press(egui::Key::Enter);
            harness.run();

            let parameter_index = if label == "A" { 0 } else { 2 };
            let app_ref = app.borrow();
            let ReactionKinetics::Elementary(data) = &app_ref.config.reactions[0].data else {
                panic!("demo reaction should use elementary kinetics");
            };
            assert_eq!(data.Arrenius[parameter_index], expected, "input `{value}`");
        }
    }

    #[test]
    fn run_calculation_updates_status_and_result_preview() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Run Calculation")
            .click_accesskit();
        for _ in 0..20 {
            harness.step();
            if app
                .borrow()
                .last_status()
                .is_some_and(|status| status.contains("Calculation completed"))
            {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }

        let app_ref = app.borrow();
        assert!(
            app_ref
                .last_status()
                .is_some_and(|status| status.contains("Calculation completed")),
            "unexpected run status: {:?}, error: {:?}",
            app_ref.last_status(),
            app_ref.last_error()
        );
        assert!(
            app_ref.last_error().is_none(),
            "unexpected run error: {:?}",
            app_ref.last_error()
        );
        assert!(app_ref.last_preview().is_some_and(
            |snapshot| !snapshot.summary_rows.is_empty() && !snapshot.equation_rows.is_empty()
        ));
    }

    #[test]
    fn successful_run_exposes_conservation_report_section() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Run Calculation")
            .click_accesskit();
        for _ in 0..20 {
            harness.step();
            if app
                .borrow()
                .last_status()
                .is_some_and(|status| status.contains("Calculation completed"))
            {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }

        harness.run();
        harness.get_by_label("Conservation");

        let app_ref = app.borrow();
        let report = app_ref
            .task()
            .latest_conservation_report()
            .expect("successful solve should publish a conservation report");
        assert!(report.energy_balance_error_abs.is_finite());
        assert!(report.energy_balance_error_rel.is_finite());
        assert!(!report.rows().is_empty());
    }

    #[test]
    fn successful_run_exposes_solve_report_section() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Run Calculation")
            .click_accesskit();
        for _ in 0..20 {
            harness.step();
            if app
                .borrow()
                .last_status()
                .is_some_and(|status| status.contains("Calculation completed"))
            {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }

        harness.run();
        harness.get_by_label("Solve report");

        let app_ref = app.borrow();
        let rows = app_ref
            .task()
            .latest_solve_report_rows()
            .expect("successful solve should publish diagnostics rows");
        assert!(!rows.is_empty());
        assert!(rows.iter().any(|row| row.field == "status"));
        assert!(rows.iter().any(|row| row.field == "controller_mode"));
    }

    #[test]
    fn failed_run_keeps_previous_preview_state_visible() {
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

        {
            let mut app_ref = app.borrow_mut();
            app_ref.config.physical.ro = -1.0;
        }

        harness
            .get_by_role_and_label(Role::Button, "Run Calculation")
            .click_accesskit();
        for _ in 0..20 {
            harness.step();
            if app
                .borrow()
                .last_status()
                .is_some_and(|status| status.contains("Calculation failed"))
            {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }

        let app_ref = app.borrow();
        assert!(app_ref.last_error().is_some());
        assert!(
            app_ref
                .last_status()
                .is_some_and(|status| status.contains("Calculation failed")),
            "unexpected failure status: {:?}, error: {:?}",
            app_ref.last_status(),
            app_ref.last_error()
        );
        assert!(app_ref.last_preview().is_some_and(
            |snapshot| !snapshot.summary_rows.is_empty() && !snapshot.equation_rows.is_empty()
        ));
    }

    #[test]
    fn aot_execution_mode_reveals_aot_controls() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        app.borrow_mut().config.execution_mode = ReactorIvpGuiExecutionMode::Aot;

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness.get_by_label("AOT build");
        harness.get_by_role_and_label(Role::ComboBox, "AOT toolchain");
        harness.get_by_role_and_label(Role::ComboBox, "AOT profile");
        harness.get_by_label("AOT confirmation");
    }

    #[test]
    fn aot_confirmation_can_be_cancelled_from_keyboard_input() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        app.borrow_mut().config.execution_mode = ReactorIvpGuiExecutionMode::Aot;
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::builder()
            .with_size(egui::Vec2::new(820.0, 620.0))
            .build_ui(move |ui| app_for_ui.borrow_mut().show(ui.ctx(), &mut open));

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Run AOT calculation")
            .focus();
        harness.run();
        harness.key_press(egui::Key::Enter);
        harness.run();
        harness.get_by_role_and_label(Role::Button, "Cancel AOT run");
        harness
            .get_by_role_and_label(Role::Button, "Cancel AOT run")
            .focus();
        harness.run();
        harness.key_press(egui::Key::Enter);
        harness.run();

        assert_eq!(
            app.borrow().last_status(),
            Some("AOT calculation cancelled before execution.")
        );
    }

    #[test]
    fn aot_run_requires_confirmation_before_solver_start() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        app.borrow_mut().config.execution_mode = ReactorIvpGuiExecutionMode::Aot;

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Run AOT calculation")
            .click_accesskit();
        harness.run();

        let app_ref = app.borrow();
        assert!(matches!(
            app_ref.calculation_state(),
            crate::gui::reactor_ivp_gui::ReactorIvpCalculationState::Idle
        ));
        assert!(
            app_ref
                .last_status()
                .is_some_and(|status| status.contains("AOT run pending confirmation"))
        );
        assert!(app_ref.last_error().is_none());
        drop(app_ref);

        harness.get_by_role_and_label(Role::Button, "Confirm AOT run");
        harness.get_by_role_and_label(Role::Button, "Cancel AOT run");
    }

    #[test]
    fn aot_confirmation_can_be_confirmed_and_runs_the_solver() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        app.borrow_mut().config.execution_mode = ReactorIvpGuiExecutionMode::Aot;

        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Run AOT calculation")
            .click_accesskit();
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Confirm AOT run")
            .click_accesskit();
        let mut dialog_closed = false;
        for _ in 0..20 {
            harness.step();
            let dialog_result = panic::catch_unwind(AssertUnwindSafe(|| {
                harness.get_by_role_and_label(Role::Button, "Cancel AOT run");
            }));
            if dialog_result.is_err() {
                dialog_closed = true;
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }

        let app_ref = app.borrow();
        assert!(
            !app_ref
                .last_status()
                .is_some_and(|status| status.contains("pending confirmation"))
        );
        assert!(!matches!(
            app_ref.calculation_state(),
            crate::gui::reactor_ivp_gui::ReactorIvpCalculationState::Idle
        ));
        assert!(
            dialog_closed,
            "AOT confirmation dialog should close after confirmation"
        );
        drop(app_ref);
    }

    #[test]
    fn lambdify_mode_hides_aot_controls() {
        let app = Rc::new(RefCell::new(ReactorIvpApp::new()));
        let mut open = true;
        let app_for_ui = Rc::clone(&app);
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });

        harness.run();
        let toolchain_result = panic::catch_unwind(AssertUnwindSafe(|| {
            harness.get_by_label("toolchain");
        }));
        let profile_result = panic::catch_unwind(AssertUnwindSafe(|| {
            harness.get_by_label("profile");
        }));
        let section_result = panic::catch_unwind(AssertUnwindSafe(|| {
            harness.get_by_label("AOT settings");
        }));

        assert!(
            toolchain_result.is_err(),
            "AOT toolchain control should be hidden in Lambdify mode"
        );
        assert!(
            profile_result.is_err(),
            "AOT profile control should be hidden in Lambdify mode"
        );
        assert!(
            section_result.is_err(),
            "AOT section should be hidden in Lambdify mode"
        );
    }

    #[test]
    fn gui_config_roundtrip_preserves_physics_and_solver_contract() {
        let app = ReactorIvpApp::new();
        let mut config = IvpGuiConfig::from_task(&app.task);
        config.problem_name = "Roundtrip".to_string();
        config.execution_mode = ReactorIvpGuiExecutionMode::Aot;
        config.rtol = 1e-7;
        config.atol = 1e-9;
        config.x_bound = 2.0;

        let mut task = app.task.clone();
        config
            .apply_to_task(&mut task)
            .expect("typed GUI config should apply");

        let rebuilt = IvpGuiConfig::from_task(&task);
        assert_eq!(rebuilt.problem_name, "Roundtrip");
        assert_eq!(rebuilt.execution_mode, ReactorIvpGuiExecutionMode::Aot);
        assert_eq!(rebuilt.rtol, 1e-7);
        assert_eq!(rebuilt.atol, 1e-9);
        assert_eq!(rebuilt.x_bound, 2.0);
        assert_eq!(rebuilt.physical.ro, app.config.physical.ro);
        assert_eq!(
            task.solver_backend_config.execution_backend,
            ReactorIvpExecutionBackend::Aot {
                toolchain: app.config.aot_toolchain,
                profile: app.config.aot_profile,
            }
        );
    }

    #[test]
    fn preview_snapshot_includes_postprocessing_export_state() {
        let app = ReactorIvpApp::new();
        let mut config = IvpGuiConfig::from_task(&app.task);
        config.problem_name = "Preview exports".to_string();
        config.postprocessing.save_csv = true;
        config.postprocessing.csv_path = "preview_exports.csv".to_string();

        let snapshot = build_reactor_ivp_task_preview_snapshot(&app.task, &config)
            .expect("preview should build");
        assert!(snapshot.summary_rows.iter().any(|row| {
            row.section == "postprocessing" && row.field == "save_csv" && row.value == "true"
        }));
        assert!(snapshot.summary_rows.iter().any(|row| {
            row.section == "postprocessing"
                && row.field == "csv_path"
                && row.value == "preview_exports.csv"
        }));
    }

    #[test]
    fn postprocessing_plan_writes_csv_and_report_for_a_solved_task() {
        let app = ReactorIvpApp::new();
        let mut task = app.task.clone();
        let mut config = IvpGuiConfig::from_task(&task);
        config.problem_name = "Postprocess export".to_string();
        config.postprocessing.save_csv = true;
        config.postprocessing.csv_path = String::new();
        config.postprocessing.write_report = true;
        config.postprocessing.report_path = String::new();
        config.postprocessing.terminal_plot = false;

        task.solve().expect("demo IVP should solve");
        let dir = tempfile::tempdir().expect("tempdir should exist");
        let csv_path = dir.path().join("postprocess_export.csv");
        let report_path = dir.path().join("postprocess_export.md");
        config.postprocessing.csv_path = csv_path.to_string_lossy().to_string();
        config.postprocessing.report_path = report_path.to_string_lossy().to_string();

        let plan = config
            .postprocessing_plan(&task)
            .expect("postprocessing plan should build")
            .expect("configured export actions should produce a plan");
        let dataset = RustedSciThe::Utils::postprocessing::PostprocessDataset::new(
            task.solver.arg_name.clone(),
            task.solver.unknowns.clone(),
            task.solver
                .x_mesh
                .clone()
                .expect("solved task should have an axis"),
            task.solver
                .solution
                .clone()
                .expect("solved task should have a matrix"),
        )
        .expect("postprocess dataset should build");

        let report = plan
            .execute(&dataset)
            .expect("configured postprocessing should succeed");
        assert!(report.all_done());
        assert!(csv_path.exists(), "CSV export should be written");
        assert!(report_path.exists(), "report export should be written");
    }
}
