use crate::Kinetics::solid_state_kinetics_IVP::KineticModelNames;
use crate::gui::gui_solid_ivp::SolidIVPApp;
use crate::gui::ivp_common::{IvpTaskState, SolverTypeShotcut};
use egui::accesskit::Role;
use egui_kittest::Harness;
use egui_kittest::kittest::Queryable;
use strum::IntoEnumIterator;

#[test]
fn solid_state_model_catalog_renders_all_formulas_without_row_overlap() {
    let mut harness = Harness::new_ui_state(
        |ui, task: &mut IvpTaskState| task.show_model_table(ui),
        IvpTaskState::new(),
    );
    harness.set_size(egui::Vec2::new(1200.0, 900.0));
    harness.run();

    let mut rows = Vec::new();
    let models: Vec<_> = KineticModelNames::iter().collect();
    assert_eq!(models.len(), 21);

    for model in models {
        let model_label = format!("{model:?}");
        let name_rect = harness.get_by_label(&model_label).rect();
        let formula = model.map_of_names_and_formulas();
        let formula_rect = harness.get_by_label(&formula).rect();
        assert!(
            name_rect.height() > 0.0,
            "{model:?} name has no visible height"
        );
        assert!(
            formula_rect.height() > 0.0,
            "{model:?} formula has no visible height"
        );
        assert!(
            name_rect.min.x >= 0.0 && name_rect.max.x <= 1200.0,
            "{model:?} name is clipped horizontally: {name_rect:?}"
        );
        assert!(
            formula_rect.min.x >= 0.0 && formula_rect.max.x <= 1200.0,
            "{model:?} formula is clipped horizontally: {formula_rect:?}"
        );
        assert!(
            !name_rect.intersects(formula_rect),
            "{model:?} name overlaps its formula"
        );
        rows.push(name_rect);
    }

    rows.sort_by(|left, right| left.min.y.total_cmp(&right.min.y));
    for pair in rows.windows(2) {
        assert!(
            pair[1].min.y >= pair[0].max.y,
            "kinetic model rows overlap: {:?} and {:?}",
            pair[0],
            pair[1]
        );
    }
}

#[derive(Default)]
struct SolidIvpGuiState {
    app: SolidIVPApp,
    open: bool,
}

#[test]
fn solid_state_window_exposes_solver_choices_and_keeps_actions_reachable() {
    let state = SolidIvpGuiState {
        app: SolidIVPApp::new(),
        open: true,
    };
    let mut harness = Harness::new_ui_state(
        |ui, state: &mut SolidIvpGuiState| state.app.show(ui.ctx(), &mut state.open),
        state,
    );
    harness.set_size(egui::Vec2::new(820.0, 620.0));
    harness.run();

    let run_rect = harness.get_by_label("Run Model").rect();
    let solver_rect = harness.get_by_label("Select Solver").rect();
    for (label, rect) in [("Run Model", run_rect), ("Select Solver", solver_rect)] {
        assert!(
            rect.min.x >= 0.0 && rect.max.x <= 820.0,
            "{label} is outside the viewport horizontally: {rect:?}"
        );
        assert!(
            rect.min.y >= 0.0 && rect.max.y <= 620.0,
            "{label} is outside the viewport vertically: {rect:?}"
        );
    }

    for (label, expected) in [
        ("Radau", SolverTypeShotcut::Radau),
        ("LSODE2", SolverTypeShotcut::LSODE2),
        ("BDF", SolverTypeShotcut::BDF),
        ("Backward Euler", SolverTypeShotcut::BackwardEuler),
    ] {
        harness.get_by_label("Select Solver").click();
        harness.run();
        harness.get_by_label(label).click();
        harness.run();
        assert_eq!(harness.state().app.task.selected_solver, expected);
    }

    harness
        .get_all_by_role_and_label(Role::Button, "Use model")
        .next()
        .expect("first model action is visible")
        .click();
    harness.run();
    assert_eq!(
        harness.state().app.task.selected_model,
        Some(KineticModelNames::A2)
    );
}
