//! Deterministic GUI-to-facade stories.
//!
//! These tests intentionally stop at request preparation. Local-catalog worker
//! stories remain in `equilibrium_gui_tests`; this module protects the public
//! document workflow without making ordinary GUI tests depend on filesystem
//! state or a nonlinear solve.

use super::equilibrium_gui::EquilibriumApp;
use super::equilibrium_gui_help::{EquilibriumHelpLanguage, text_with_fallback, tooltip};
use super::equilibrium_gui_model::{
    CandidatePolicyDraft, ComponentDraft, EquilibriumGuiDocument, EquilibriumInventoryDraft,
    EquilibriumLookupDraft, EquilibriumPhaseModeDraft, EquilibriumProblemDraft,
    EquilibriumSolverDraft, GuiElementSearchMode, GuiInterpolationSpace, GuiPhSolveMode,
    GuiPhaseModel, GuiPhysicalState, GuiResamplingDraft, GuiSolverBackend,
    PhTemperatureBoundsDraft, PhaseDraft, TemperatureDraft,
};
use super::equilibrium_gui_request::{
    EquilibriumGuiSolveRequest, build_equilibrium_facade_request,
};
use egui::accesskit::Role;
use egui_kittest::Harness;
use egui_kittest::kittest::Queryable;
use std::cell::RefCell;
use std::rc::Rc;

fn prepared_request(document: &EquilibriumGuiDocument) -> EquilibriumGuiSolveRequest {
    let validated = document
        .validate_for_run()
        .expect("story document must pass GUI validation");
    build_equilibrium_facade_request(validated, None)
        .expect("validated document must produce a facade request")
}

#[test]
fn default_point_story_prepares_the_canonical_facade_request() {
    let document = EquilibriumGuiDocument::new();

    assert!(matches!(
        prepared_request(&document),
        EquilibriumGuiSolveRequest::Facade(_)
    ));
}

#[test]
fn pt_range_story_preserves_range_shape_at_the_facade_boundary() {
    let mut document = EquilibriumGuiDocument::new();
    let EquilibriumProblemDraft::FixedPt { temperature, .. } = &mut document.config.problem else {
        panic!("default story document must use P,T");
    };
    *temperature = TemperatureDraft::Range {
        start_k: "300".into(),
        end_k: "1200".into(),
        point_count: "5".into(),
    };

    assert!(matches!(
        prepared_request(&document),
        EquilibriumGuiSolveRequest::Facade(_)
    ));
}

#[test]
fn ph_point_story_prepares_without_exposing_a_second_request_path() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "100000".into(),
        target_enthalpy_j: "1000".into(),
        temperature_bounds: PhTemperatureBoundsDraft {
            lower_k: "250".into(),
            upper_k: "2000".into(),
            seed_k: "1000".into(),
        },
    };

    assert!(matches!(
        prepared_request(&document),
        EquilibriumGuiSolveRequest::Facade(_)
    ));
}

#[test]
fn invalid_pressure_story_stops_before_facade_request_preparation() {
    let mut document = EquilibriumGuiDocument::new();
    let EquilibriumProblemDraft::FixedPt { pressure_pa, .. } = &mut document.config.problem else {
        panic!("default story document must use P,T");
    };
    *pressure_pa = "not-a-pressure".into();

    let report = document
        .validate_for_run()
        .expect_err("invalid pressure must remain a field-local GUI error");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "problem.pressure_pa")
    );
}

#[test]
fn element_candidate_story_requires_assignment_before_preparing() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.inventory = EquilibriumInventoryDraft::ElementCandidates {
        elements: vec!["H".into(), "O".into()],
        candidate_policy: CandidatePolicyDraft {
            element_mode: GuiElementSearchMode::Exact,
            ..CandidatePolicyDraft::default()
        },
        assignments: Vec::new(),
    };
    let report = document
        .validate_for_run()
        .expect_err("unassigned candidates must not produce a runnable request");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "inventory.assignments")
    );

    document.config.inventory = EquilibriumInventoryDraft::ElementCandidates {
        elements: vec!["H".into(), "O".into()],
        candidate_policy: CandidatePolicyDraft::default(),
        assignments: vec![PhaseDraft::default()],
    };
    assert!(matches!(
        prepared_request(&document),
        EquilibriumGuiSolveRequest::Facade(_)
    ));
}

#[test]
fn empty_inventory_story_is_rejected_before_request_construction() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.inventory = EquilibriumInventoryDraft::ExplicitPhases { phases: Vec::new() };
    let report = document
        .validate_for_run()
        .expect_err("an empty inventory cannot be prepared");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "inventory.phases")
    );
}

#[test]
fn incompatible_phase_state_and_model_story_is_field_local() {
    let mut document = EquilibriumGuiDocument::new();
    let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut document.config.inventory
    else {
        panic!("default inventory must be explicit");
    };
    phases[0].physical_state = GuiPhysicalState::Solid;
    phases[0].model = GuiPhaseModel::IdealGas;
    let report = document
        .validate_for_run()
        .expect_err("solid ideal gas must be rejected");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "inventory.phases[0].model")
    );
}

#[test]
fn ph_route_override_story_is_ignored_by_the_pt_request_contract() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.solver.ph_solve_mode = GuiPhSolveMode::Monolithic;
    document.config.solver.selection = EquilibriumSolverDraft::ProductionDefault;
    // The route selector is presentation-scoped to P,H. A stale serialized
    // value must not make an otherwise valid P,T document unrunnable.
    assert!(matches!(
        prepared_request(&document),
        EquilibriumGuiSolveRequest::Facade(_)
    ));
}

#[test]
fn default_setup_story_exposes_only_the_primary_calculator_controls() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_role_and_label(Role::Button, "Setup");
    harness.get_by_label("P,T = const");
    harness.get_by_label("Explicit species");
    assert!(harness.query_by_label("Concrete backend").is_none());
    assert!(harness.query_by_label("Collect timing").is_none());
    assert!(matches!(
        prepared_request(&app.borrow().document),
        EquilibriumGuiSolveRequest::Facade(_)
    ));
}

#[test]
fn default_widget_census_is_independent_from_the_tooltip_catalog() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();

    // These are rendered accessibility labels, deliberately not tooltip IDs.
    // If a primary control vanishes or is accidentally moved behind an advanced
    // tab, this fails even when the static tooltip catalog is still complete.
    for label in [
        "Setup",
        "Phase control",
        "Libraries",
        "Numerics",
        "Output",
        "Results",
        "English",
        "Русский",
        "Problem",
        "Calculation mode",
        "P,T = const",
        "P,H = const",
        "Pressure [Pa]",
        "Reference pressure [Pa]",
        "Temperature [K]",
        "Components and phases",
        "Inventory input",
        "Explicit species",
        "Search by elements",
        "Simple ideal-gas preset",
        "Validate document",
        "Prepare canonical request",
    ] {
        harness.get_by_label(label);
    }

    for advanced_only in ["Concrete backend", "Collect timing", "Library lookup"] {
        assert!(
            harness.query_by_label(advanced_only).is_none(),
            "{advanced_only} must not leak into the default Setup viewport"
        );
    }
}

#[test]
fn conditional_widget_census_covers_problem_inventory_solver_phase_and_output_states() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    // P,T range controls replace the single-temperature editor.
    if let EquilibriumProblemDraft::FixedPt { temperature, .. } =
        &mut app.borrow_mut().document.config.problem
    {
        *temperature = TemperatureDraft::Range {
            start_k: "300".into(),
            end_k: "800".into(),
            point_count: "4".into(),
        };
    }
    harness.run();
    for label in ["Start [K]", "End [K]", "Solved points"] {
        harness.get_by_label(label);
    }
    assert!(harness.query_by_label("Temperature [K]").is_none());

    // P,H replaces P,T temperature controls with the enthalpy contract.
    app.borrow_mut().document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "100000".into(),
        target_enthalpy_j: "1000".into(),
        temperature_bounds: PhTemperatureBoundsDraft::default(),
    };
    harness.run();
    for label in [
        "Target total enthalpy [J]",
        "Lower bound [K]",
        "Upper bound [K]",
        "Initial seed [K]",
    ] {
        harness.get_by_label(label);
    }
    assert!(harness.query_by_label("Start [K]").is_none());

    // Element search is a different inventory surface from direct components.
    app.borrow_mut().document.config.inventory = EquilibriumInventoryDraft::ElementCandidates {
        elements: vec!["H".into(), "O".into()],
        candidate_policy: CandidatePolicyDraft::default(),
        assignments: vec![PhaseDraft::default()],
    };
    harness.run();
    harness.get_by_label("Element matching");

    // Bounded control exposes lifecycle controls that fixed mode must not show.
    app.borrow_mut().document.config.phase_mode = EquilibriumPhaseModeDraft::Bounded {
        phase_epsilon: "1e-12".into(),
        dg_create: "-1e-6".into(),
        dg_keep: "1e-8".into(),
        max_phase_iterations: "20".into(),
        initial_phase_policy:
            super::equilibrium_gui_model::GuiInitialPhasePolicyDraft::FromInitialMoles,
    };
    harness
        .get_by_role_and_label(Role::Button, "Phase control")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase policy")
        .click_accesskit();
    harness.run();
    for label in [
        "Phase epsilon",
        "Creation driving force",
        "Keep driving force",
        "Initial phase set",
    ] {
        harness.get_by_label(label);
    }

    // A custom cascade has a dedicated editor rather than the default selector.
    app.borrow_mut().document.config.solver.selection = EquilibriumSolverDraft::CustomCascade {
        backends: Vec::new(),
    };
    harness
        .get_by_role_and_label(Role::Button, "Numerics")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Ordered fallback sequence");
    harness.get_by_label("Add backend");

    // PCHIP belongs only to accepted P,T-range presentation, never P,H.
    app.borrow_mut().document.config.problem = EquilibriumProblemDraft::FixedPt {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "100000".into(),
        temperature: TemperatureDraft::Range {
            start_k: "300".into(),
            end_k: "800".into(),
            point_count: "4".into(),
        },
    };
    app.borrow_mut().document.config.postprocessing.resampling = GuiResamplingDraft::Pchip {
        output_points: "20".into(),
        interpolation_space: GuiInterpolationSpace::Linear,
        clamp: true,
    };
    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Postprocessing and plots")
        .click_accesskit();
    harness.run();
    for label in [
        "PCHIP display resampling",
        "Display points",
        "Interpolation space",
    ] {
        harness.get_by_label(label);
    }
}

#[test]
fn results_census_renders_running_cancelling_and_failed_states() {
    fn assert_results_state(
        app: Rc<RefCell<EquilibriumApp>>,
        expected: &str,
        detail: Option<&str>,
    ) {
        let app_for_ui = Rc::clone(&app);
        let mut open = true;
        let mut harness = Harness::new_ui(move |ui| {
            app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
        });
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Results")
            .click_accesskit();
        harness.run();
        harness.get_by_label_contains(expected);
        if let Some(detail) = detail {
            harness.get_by_label_contains(detail);
        }
    }

    let running = Rc::new(RefCell::new(EquilibriumApp::new()));
    running
        .borrow_mut()
        .prepare_request()
        .expect("default request prepares");
    running
        .borrow_mut()
        .begin_prepared_run()
        .expect("running state is available");
    assert_results_state(Rc::clone(&running), "Run state: Solving", None);
    assert!(running.borrow().result_snapshot().is_none());

    let cancelling = Rc::new(RefCell::new(EquilibriumApp::new()));
    cancelling
        .borrow_mut()
        .prepare_request()
        .expect("cancelling request prepares");
    cancelling
        .borrow_mut()
        .begin_prepared_run()
        .expect("cancelling state is available");
    assert!(cancelling.borrow_mut().cancel_run());
    assert_results_state(Rc::clone(&cancelling), "Run state: Cancelling", None);
    assert!(cancelling.borrow().result_snapshot().is_none());

    let failed = Rc::new(RefCell::new(EquilibriumApp::new()));
    failed
        .borrow_mut()
        .prepare_request()
        .expect("failed-state request prepares");
    let ticket = failed
        .borrow_mut()
        .begin_prepared_run()
        .expect("failed-state ticket is available");
    assert_eq!(
        failed
            .borrow_mut()
            .publish_failure(ticket, "deterministic story failure"),
        super::equilibrium_gui_execution::EquilibriumGuiPublication::Accepted
    );
    assert_results_state(Rc::clone(&failed), "Run state: Failed", None);
    assert_eq!(
        failed.borrow().last_error(),
        Some("deterministic story failure")
    );
    assert!(failed.borrow().result_snapshot().is_none());
}

#[test]
fn hover_tooltips_cover_combo_options_and_dynamic_rows() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    app.borrow_mut().document.config.solver.selection = EquilibriumSolverDraft::SingleBackend {
        backend: GuiSolverBackend::RstLm,
    };
    let document_before = app.borrow().document.clone();
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Phase id").hover();
    harness.run();
    harness.get_by_label_contains("Stable phase identifier used to distinguish");
    harness.get_by_label("Substance").hover();
    harness.run();
    harness.get_by_label_contains("Thermochemical component identifier");

    harness
        .get_by_role_and_label(Role::Button, "Numerics")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Concrete backend").click_accesskit();
    harness.run();
    harness.get_by_label("Legacy LM fallback").hover();
    harness.run();
    harness.get_by_label_contains("one concrete nonlinear backend");

    assert_eq!(app.borrow().document, document_before);
    assert!(app.borrow().result_snapshot().is_none());
}

#[test]
fn secondary_tab_story_keeps_advanced_controls_out_of_setup() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Diagnostics")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Collect timing");
    assert!(harness.query_by_label("Explicit species").is_none());

    harness
        .get_by_role_and_label(Role::Button, "Libraries")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Library lookup");
    assert!(harness.query_by_label("Collect timing").is_none());
}

#[test]
fn route_conditional_controls_story_switches_between_pt_and_ph() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Temperature [K]");
    assert!(
        harness
            .query_by_label("Target total enthalpy [J]")
            .is_none()
    );

    app.borrow_mut().document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "100000".into(),
        target_enthalpy_j: "1000".into(),
        temperature_bounds: PhTemperatureBoundsDraft::default(),
    };
    harness.run();
    harness.get_by_label("Target total enthalpy [J]");
    assert!(harness.query_by_label("Temperature [K]").is_none());
    assert!(tooltip(EquilibriumHelpLanguage::English, "mode.ph").contains("enthalpy"));
    assert!(!tooltip(EquilibriumHelpLanguage::Russian, "mode.ph").is_empty());
}

#[test]
fn conditional_tooltip_story_has_distinct_range_and_bounded_explanations() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    assert!(tooltip(EquilibriumHelpLanguage::English, "output.pchip").contains("Resample"));
    assert!(tooltip(EquilibriumHelpLanguage::English, "phase.initial_set").contains("active"));
    assert!(tooltip(EquilibriumHelpLanguage::Russian, "output.pchip").contains("PCHIP"));
    assert!(app.borrow().document.validate().is_ok());
}

#[test]
fn hover_tooltip_story_renders_the_catalog_text_and_localizes_it() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Validate document").hover();
    harness.run();
    harness.get_by_label_contains("without starting a solve");

    assert!(!tooltip(EquilibriumHelpLanguage::Russian, "action.validate").is_empty());
    assert!(app.borrow().document.validate().is_ok());
}

#[test]
fn hover_operations_are_read_only_across_input_selector_action_and_checkbox() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let document_before = app.borrow().document.clone();
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Pressure [Pa]").hover();
    harness.run();
    harness.get_by_label_contains("System pressure in pascals");

    harness.get_by_label("P,T = const").hover();
    harness.run();
    harness.get_by_label_contains("fixed pressure and temperature");

    harness.get_by_label("Simple ideal-gas preset").hover();
    harness.run();
    harness.get_by_label_contains("small valid ideal-gas example");

    harness
        .get_by_role_and_label(Role::Button, "Phase control")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase policy")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Bounded phase control").hover();
    harness.run();
    harness.get_by_label_contains("activation and deactivation");

    let app = app.borrow();
    assert_eq!(app.document, document_before);
    assert!(app.result_snapshot().is_none());
    assert!(matches!(
        app.run_state(),
        super::equilibrium_gui_execution::EquilibriumGuiRunState::Idle
    ));
}

#[test]
fn hover_tooltip_story_covers_conditional_temperature_control() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Point").hover();
    harness.run();
    harness.get_by_label_contains("fixed temperature in kelvin");

    harness.get_by_label("Range").click_accesskit();
    harness.run();
    harness.get_by_label("Range").hover();
    harness.run();
    harness.get_by_label_contains("ordered temperature grid");
    assert!(app.borrow().document.validate().is_ok());
}

#[test]
fn hover_tooltip_story_localizes_temperature_mode_controls() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let before = app.borrow().document.clone();
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Русский").click_accesskit();
    harness.run();
    harness.get_by_label("Point").hover();
    harness.run();
    harness.get_by_label_contains("фиксированной температуре");
    assert_eq!(app.borrow().document, before);
}

#[test]
fn hover_tooltip_story_covers_phase_policy_switch() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase control")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase policy")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Bounded phase control").hover();
    harness.run();
    harness.get_by_label_contains("activation and deactivation");
    assert!(app.borrow().document.validate().is_ok());
}

#[test]
fn hover_tooltip_story_covers_solver_and_output_tabs() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Numerics").click_accesskit();
    harness.run();
    harness.get_by_label("Numerics").hover();
    harness.run();
    harness.get_by_label_contains("nonlinear backend cascade");

    harness.get_by_label("Output").click_accesskit();
    harness.run();
    harness.get_by_label("Output").hover();
    harness.run();
    harness.get_by_label_contains("diagnostics");
    assert!(app.borrow().document.validate().is_ok());
}

#[test]
fn hover_tooltip_story_covers_lookup_and_presentation_selectors() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Libraries").click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Library lookup")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Load local library choices").hover();
    harness.run();
    harness.get_by_label_contains("without resolving substances");

    harness.get_by_label("Output").click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Postprocessing and plots")
        .click_accesskit();
    harness.run();
    harness.get_by_label("None").hover();
    harness.run();
    harness.get_by_label_contains("Do not open a plot");
    assert!(app.borrow().document.validate().is_ok());
}

#[test]
fn ph_bounds_story_rejects_an_out_of_bracket_seed() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "100000".into(),
        target_enthalpy_j: "1000".into(),
        temperature_bounds: PhTemperatureBoundsDraft {
            lower_k: "300".into(),
            upper_k: "1000".into(),
            seed_k: "1200".into(),
        },
    };
    let report = document
        .validate_for_run()
        .expect_err("P,H seed outside the bracket must be rejected");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "problem.temperature_bounds.seed_k")
    );
}

#[test]
fn explicit_empty_lookup_policy_story_remains_explicit_until_resolution() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.lookup = EquilibriumLookupDraft::Explicit {
        priority_libraries: Vec::new(),
        permitted_libraries: Vec::new(),
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };
    // Empty lists are a valid explicit "no permitted catalog" policy. The
    // repository resolver, not GUI validation, reports that no record can be
    // resolved under it.
    assert!(matches!(
        prepared_request(&document),
        EquilibriumGuiSolveRequest::Facade(_)
    ));
}

#[test]
fn display_cutoff_story_rejects_negative_values_without_touching_solver_policy() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.postprocessing.display_fraction_cutoff = "-1e-6".into();
    let report = document
        .validate_for_run()
        .expect_err("display cutoff must be non-negative");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "postprocessing.display_fraction_cutoff")
    );
}

#[test]
fn negative_inventory_story_is_rejected_at_the_component_field() {
    let mut document = EquilibriumGuiDocument::new();
    let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut document.config.inventory
    else {
        panic!("default inventory must be explicit");
    };
    phases[0].components[0].initial_moles = "-1e-3".into();
    let report = document
        .validate_for_run()
        .expect_err("negative inventory must be rejected");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| { issue.field == "inventory.phases[0].components[0].initial_moles" })
    );
}

#[test]
fn duplicate_phase_story_is_rejected_before_request_construction() {
    let mut document = EquilibriumGuiDocument::new();
    let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut document.config.inventory
    else {
        panic!("default inventory must be explicit");
    };
    phases.push(phases[0].clone());
    let report = document
        .validate_for_run()
        .expect_err("duplicate phase IDs must be rejected");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "inventory.phases[1].id")
    );
}

#[test]
fn duplicate_phase_component_story_is_rejected_before_request_construction() {
    let mut document = EquilibriumGuiDocument::new();
    let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut document.config.inventory
    else {
        panic!("default inventory must be explicit");
    };
    let substance = phases[0].components[0].substance.clone();
    phases[0].components.push(ComponentDraft {
        substance,
        initial_moles: "1".into(),
        source_library: None,
    });
    let report = document
        .validate_for_run()
        .expect_err("duplicate phase-qualified components must be rejected");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| { issue.field == "inventory.phases[0].components[1]" })
    );
}

#[test]
fn non_finite_inventory_story_is_rejected_at_the_amount_field() {
    for value in ["NaN", "inf", "-inf"] {
        let mut document = EquilibriumGuiDocument::new();
        let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut document.config.inventory
        else {
            panic!("default inventory must be explicit");
        };
        phases[0].components[0].initial_moles = value.into();
        let report = document
            .validate_for_run()
            .expect_err("non-finite amount must be rejected");
        assert!(
            report
                .issues
                .iter()
                .any(|issue| { issue.field == "inventory.phases[0].components[0].initial_moles" })
        );
    }
}

#[test]
fn reversed_pt_range_story_is_rejected_before_request_construction() {
    let mut document = EquilibriumGuiDocument::new();
    let EquilibriumProblemDraft::FixedPt { temperature, .. } = &mut document.config.problem else {
        panic!("default story document must use P,T");
    };
    *temperature = TemperatureDraft::Range {
        start_k: "1000".into(),
        end_k: "1000".into(),
        point_count: "3".into(),
    };
    let report = document
        .validate_for_run()
        .expect_err("a zero-width range must be rejected");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field.contains("temperature")),
        "unexpected P,T range issues: {:?}",
        report.issues
    );
}

#[test]
fn blank_library_name_story_is_rejected_at_lookup_policy() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.lookup = EquilibriumLookupDraft::Explicit {
        priority_libraries: vec!["".into()],
        permitted_libraries: vec!["NASA_gas".into()],
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };
    let report = document
        .validate_for_run()
        .expect_err("blank library names must be rejected");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field.starts_with("lookup.")),
        "unexpected lookup issues: {:?}",
        report.issues
    );
}

#[test]
fn invalid_pchip_count_story_is_rejected_as_presentation_input() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.postprocessing.resampling = GuiResamplingDraft::Pchip {
        output_points: "1".into(),
        interpolation_space: GuiInterpolationSpace::Linear,
        clamp: false,
    };
    let report = document
        .validate_for_run()
        .expect_err("PCHIP needs at least two output points");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "postprocessing.resampling.output_points")
    );
}

#[test]
fn phase_control_story_exposes_bounded_policy_only_when_selected() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase control")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase policy")
        .click_accesskit();
    harness.run();
    assert!(harness.query_by_label("Creation driving force").is_none());

    app.borrow_mut().document.config.phase_mode = EquilibriumPhaseModeDraft::Bounded {
        phase_epsilon: "1e-12".into(),
        dg_create: "-1e-6".into(),
        dg_keep: "1e-8".into(),
        max_phase_iterations: "20".into(),
        initial_phase_policy:
            super::equilibrium_gui_model::GuiInitialPhasePolicyDraft::FromInitialMoles,
    };
    harness.run();
    harness.get_by_label("Creation driving force");
    harness.get_by_label("Keep driving force");
    harness.get_by_label("Initial phase set");
}

#[test]
fn output_story_exposes_presentation_controls_without_invalidating_the_document() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Postprocessing and plots")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Result table");
    harness.get_by_label("Hide below mole fraction (display only)");
    harness.get_by_label("Plot target");
    harness.get_by_label("Result basis");
    assert!(app.borrow().document.validate().is_ok());
}

#[test]
fn contextual_help_story_is_available_on_primary_and_secondary_tabs() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let before = app.borrow().document.clone();
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_role_and_label(Role::Button, "Help");
    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness.get_by_role_and_label(Role::Button, "Help");
    assert_eq!(app.borrow().document, before);
}

#[test]
fn pchip_story_is_visible_only_for_pt_ranges_and_is_display_only() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Postprocessing and plots")
        .click_accesskit();
    harness.run();
    assert!(harness.query_by_label("PCHIP display resampling").is_none());

    if let EquilibriumProblemDraft::FixedPt { temperature, .. } =
        &mut app.borrow_mut().document.config.problem
    {
        *temperature = TemperatureDraft::Range {
            start_k: "300".into(),
            end_k: "1000".into(),
            point_count: "5".into(),
        };
    }
    harness.run();
    app.borrow_mut()
        .prepare_request()
        .expect("range request prepares");
    let before_request = app.borrow().prepared_request_is_current();
    harness
        .get_by_label("PCHIP display resampling")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Display points");
    assert!(app.borrow().prepared_request_is_current() == before_request);
    assert!(app.borrow().document.validate().is_ok());
}

#[test]
fn help_language_story_switches_resources_without_changing_the_document() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let before = app.borrow().document.clone();
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Help language");
    harness.get_by_label("Русский").click_accesskit();
    harness.run();
    harness.get_by_label("Help");
    assert_eq!(app.borrow().document, before);
    assert!(!text_with_fallback(EquilibriumHelpLanguage::Russian, "Setup").is_empty());
}

#[test]
fn empty_results_story_is_read_only_and_plotting_fails_without_a_snapshot() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let before = app.borrow().document.clone();
    assert!(app.borrow_mut().open_embedded_plot().is_err());

    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Results")
        .click_accesskit();
    harness.run();
    harness.get_by_label("No accepted equilibrium result yet");
    assert_eq!(app.borrow().document, before);
}
