use super::*;
use crate::Kinetics::reaction_network_ivp::solve_network;
use egui::accesskit::Role;
use egui_kittest::{Harness, kittest::Queryable};
use std::{cell::RefCell, rc::Rc};

#[test]
fn default_document_builds_and_solves_through_public_facade() {
    let app = NetworkIvpApp::new();
    let problem = app.build_problem().unwrap();
    let options = app.build_options().unwrap();
    let result = solve_network(&problem, &options).unwrap();
    assert_eq!(result.species_ids, ["A", "B"]);
    assert!(result.concentrations.last().unwrap()[0] < 1.0);
}

#[test]
fn series_conditions_parse_heating_rates_from_kelvin_per_minute() {
    let mut app = NetworkIvpApp::new_for_experimental_kinetics();
    app.temperature_kind = TemperatureKind::Linear;
    app.temperature = "500".into();
    app.series_values = "1, 2.3".into();
    let bindings = app.parse_series_bindings().unwrap();
    assert_eq!(bindings.len(), 2);
    assert_eq!(
        bindings[0],
        TemperatureBinding::Linear {
            initial_kelvin: 500.0,
            kelvin_per_second: 1.0 / 60.0,
        }
    );
    assert_eq!(
        bindings[1],
        TemperatureBinding::Linear {
            initial_kelvin: 500.0,
            kelvin_per_second: 2.3 / 60.0,
        }
    );
}

#[test]
fn series_condition_parser_rejects_empty_non_finite_and_oversized_lists() {
    let mut app = NetworkIvpApp::new_for_experimental_kinetics();
    app.series_values = ", ;".into();
    assert!(app.parse_series_bindings().is_err());
    app.series_values = "400, NaN".into();
    assert!(app.parse_series_bindings().is_err());
    app.series_values = (0..101)
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",");
    assert!(app.parse_series_bindings().unwrap_err().contains("100"));
}

#[test]
fn series_condition_preflight_rejects_existing_output_id() {
    let mut app = NetworkIvpApp::new_for_experimental_kinetics();
    app.series_values = "400, 500".into();
    app.set_existing_experiment_ids(["network_500K".to_string()]);
    assert!(
        app.parse_series_bindings()
            .unwrap_err()
            .contains("network_500K")
    );
}

#[test]
fn invalid_document_reports_validation_without_replacing_last_result() {
    let mut app = NetworkIvpApp::new();
    app.species[0].initial = "-1".into();
    assert!(
        app.build_problem()
            .unwrap_err()
            .contains("initial_concentration")
    );
    app.species[0].initial = "1".into();
    app.reactions[0].kind = RateKind::Custom;
    app.reactions[0].custom = "bad_symbol".into();
    assert!(app.build_problem().unwrap_err().contains("bad_symbol"));
}

#[test]
fn solver_and_task_fields_are_independent() {
    let mut app = NetworkIvpApp::new();
    app.family = StiffSolverKind::Radau;
    app.temperature_kind = TemperatureKind::Linear;
    app.heating_rate = "2".into();
    assert!(matches!(
        app.build_problem().unwrap().temperature(),
        TemperatureProgram::Linear { .. }
    ));
    assert_eq!(app.build_options().unwrap().family, StiffSolverKind::Radau);
}

#[test]
fn custom_rate_does_not_parse_hidden_coefficient_fields() {
    let mut app = NetworkIvpApp::new();
    app.reactions[0].kind = RateKind::Custom;
    app.reactions[0].custom = "0.5*c0".into();
    app.reactions[0].prefactor = "not a number".into();
    app.reactions[0].energy = "not a number".into();
    assert!(app.build_problem().is_ok());
}

#[test]
fn backward_euler_hides_unsupported_common_controls_without_changing_the_task() {
    let mut app = NetworkIvpApp::new();
    app.family = StiffSolverKind::BackwardEuler;
    app.callback = CallbackKind::Parallel;
    app.rtol = "not used".into();
    app.atol = "not used".into();
    let options = app.build_options().unwrap();
    assert_eq!(options.callback_policy, None);
    assert_eq!(options.rtol, 1e-8);
    assert_eq!(options.atol, 1e-10);
    assert!(app.build_problem().is_ok());
}

#[test]
fn rendered_tabs_separate_task_from_solver_and_run_the_default_case() {
    struct State {
        app: NetworkIvpApp,
        open: bool,
    }
    let state = Rc::new(RefCell::new(State {
        app: NetworkIvpApp::new(),
        open: true,
    }));
    let for_ui = Rc::clone(&state);
    let mut harness = Harness::new_ui(move |ui| {
        let mut state = for_ui.borrow_mut();
        let State { app, open } = &mut *state;
        app.show(ui.ctx(), open);
    });
    harness.run();
    harness.get_by_role_and_label(Role::Window, "Reaction network IVP");
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness.run();
    assert_eq!(state.borrow().app.page, Page::Solver);
    harness
        .get_by_role_and_label(Role::Button, "Task")
        .click_accesskit();
    harness.run();
    assert_eq!(state.borrow().app.page, Page::Physics);
    harness
        .get_by_role_and_label(Role::Button, "Run")
        .click_accesskit();
    for _ in 0..50 {
        harness.run();
        if state.borrow().app.result.is_some() {
            break;
        }
        std::thread::sleep(Duration::from_millis(10));
    }
    let state = state.borrow();
    assert!(state.app.result.is_some(), "{:?}", state.app.message);
    assert_eq!(state.app.page, Page::Result);
}

#[test]
fn rendered_temperature_selector_exposes_only_the_selected_program_inputs() {
    let state = Rc::new(RefCell::new((NetworkIvpApp::new(), true)));
    let for_ui = Rc::clone(&state);
    let mut harness = Harness::builder()
        .with_size(egui::Vec2::new(1100.0, 820.0))
        .build_ui(move |ui| {
            let mut state = for_ui.borrow_mut();
            let (app, open) = &mut *state;
            app.show(ui.ctx(), open);
        });
    harness.run();

    harness.get_by_role_and_label(Role::ComboBox, "Temperature");
    state.borrow_mut().0.temperature_kind = TemperatureKind::Linear;
    harness.run();
    assert_eq!(state.borrow().0.temperature_kind, TemperatureKind::Linear);
    harness.get_by_role_and_label(Role::TextInput, "K/s");
    state.borrow_mut().0.temperature_kind = TemperatureKind::Symbolic;
    harness.run();
    assert_eq!(state.borrow().0.temperature_kind, TemperatureKind::Symbolic);
    harness.get_by_role_and_label(Role::TextInput, "T(t), K");
}

#[test]
fn rendered_solver_family_updates_conditional_controls() {
    let state = Rc::new(RefCell::new((NetworkIvpApp::new(), true)));
    let for_ui = Rc::clone(&state);
    let mut harness = Harness::builder()
        .with_size(egui::Vec2::new(1100.0, 820.0))
        .build_ui(move |ui| {
            let mut state = for_ui.borrow_mut();
            let (app, open) = &mut *state;
            app.show(ui.ctx(), open);
        });
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness.run();

    harness.get_by_role_and_label(Role::ComboBox, "Family");
    state.borrow_mut().0.family = StiffSolverKind::BackwardEuler;
    harness.run();
    assert_eq!(state.borrow().0.family, StiffSolverKind::BackwardEuler);
    harness.get_by_label("Unavailable");

    state.borrow_mut().0.family = StiffSolverKind::Radau;
    harness.run();
    assert_eq!(state.borrow().0.family, StiffSolverKind::Radau);
    harness.get_by_role_and_label(Role::ComboBox, "Callbacks");
}

#[test]
fn rendered_species_and_reaction_controls_add_and_remove_rows() {
    let state = Rc::new(RefCell::new((NetworkIvpApp::new(), true)));
    let for_ui = Rc::clone(&state);
    let mut harness = Harness::builder()
        .with_size(egui::Vec2::new(1200.0, 900.0))
        .build_ui(move |ui| {
            let mut state = for_ui.borrow_mut();
            let (app, open) = &mut *state;
            app.show(ui.ctx(), open);
        });
    harness.run();

    harness
        .get_by_role_and_label(Role::Button, "Add species")
        .click_accesskit();
    harness.run();
    assert_eq!(state.borrow().0.species.len(), 3);
    harness.get_by_role_and_label(Role::TextInput, "Species 3 ID");

    harness
        .get_by_role_and_label(Role::Button, "Add reaction")
        .click_accesskit();
    harness.run();
    assert_eq!(state.borrow().0.reactions.len(), 2);
}

#[test]
fn rendered_help_has_all_sections_and_switches_language() {
    let state = Rc::new(RefCell::new((NetworkIvpApp::new(), true)));
    let for_ui = Rc::clone(&state);
    let mut harness = Harness::builder()
        .with_size(egui::Vec2::new(1100.0, 820.0))
        .build_ui(move |ui| {
            let mut state = for_ui.borrow_mut();
            let (app, open) = &mut *state;
            app.show(ui.ctx(), open);
        });
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Help")
        .click_accesskit();
    harness.run();

    for section in [
        "Физическая модель",
        "Доступные солверы",
        "Основные настройки",
        "Типичный workflow",
    ] {
        harness.get_by_label(section);
    }
    harness
        .get_by_role_and_label(Role::Button, "English")
        .click_accesskit();
    harness.run();
    assert_eq!(state.borrow().0.help_language, HelpLanguage::English);
    for section in [
        "Physical model",
        "Available solvers",
        "Main settings",
        "Typical workflow",
    ] {
        harness.get_by_label(section);
    }
}
