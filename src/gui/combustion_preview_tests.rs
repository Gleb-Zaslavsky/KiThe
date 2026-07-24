//! Story and regression tests for the BVP Preview Task action.
//!
//! Hypotheses:
//! - Preview Task should build a read-only normalized snapshot of the current BVP document.
//! - The snapshot should include solver_settings and keep the canonical section ordering stable.
//! - Equation preview rows should come from `Expr::pretty_print()` on the assembled local reactor.
//! - Clicking Preview Task should update the GUI status without mutating the loaded document.
//!
//! Expected result:
//! - The BVP GUI exposes a safe inspection path before running the solver.

#[cfg(test)]
mod tests {
    use super::super::combustion::{
        CombustionApp, ProblemsEnum, build_bvp_task_preview_snapshot,
        collect_bvp_preview_equation_rows,
    };
    use crate::ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;
    use RustedSciThe::command_interpreter::task_parser::{DocumentMap, DocumentParser};
    use egui::accesskit::Role;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;
    use std::cell::RefCell;
    use std::rc::Rc;

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
    fn preview_snapshot_keeps_document_read_only_and_lists_solver_settings() {
        let document = load_problem_test_document();
        let before = document.clone();

        let snapshot = build_bvp_task_preview_snapshot(&document, &ProblemsEnum::BVPSimple);

        assert_eq!(before, document);
        assert!(snapshot.validation_report.is_valid());
        assert!(!snapshot.parameter_rows.is_empty());
        assert!(
            snapshot
                .parameter_rows
                .iter()
                .any(|row| row.section == "solver_settings")
        );
        assert_eq!(
            snapshot
                .parameter_rows
                .first()
                .map(|row| row.section.as_str()),
            Some("solver_settings")
        );
        assert!(snapshot.equation_build_error.is_none());
        assert!(!snapshot.equation_rows.is_empty());
    }

    #[test]
    fn preview_equations_use_pretty_printed_solver_output() {
        let document = load_problem_test_document();
        let mut reactor = SimpleReactorTask::new();

        reactor
            .set_reactor_params_from_hashmap(&document)
            .expect("problem_test.txt should populate the reactor");
        reactor
            .setup_bvp()
            .expect("problem_test.txt should assemble equations");

        let rows = collect_bvp_preview_equation_rows(&reactor);
        assert_eq!(
            rows.len(),
            reactor
                .solver
                .eq_system
                .len()
                .min(reactor.solver.unknowns.len())
        );

        let first_row = rows
            .first()
            .expect("assembled equations should not be empty");
        assert_eq!(first_row.unknown, reactor.solver.unknowns[0]);
        assert_eq!(
            first_row.equation,
            reactor.solver.eq_system[0].pretty_print()
        );
    }

    #[test]
    fn preview_button_updates_status_without_mutating_document() {
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
        let baseline_document = app.borrow().document.clone();
        let preview_button = harness.get_by_role_and_label(Role::Button, "Preview Task");
        preview_button.scroll_to_me();
        // Let AccessKit scroll settle before sending the pointer click. The
        // preview button lives near the bottom of the BVP screen, so the
        // layout can shift after the scroll request.
        harness.run();
        harness
            .get_by_role_and_label(Role::Button, "Preview Task")
            .click_accesskit();
        harness.run();

        assert_eq!(baseline_document, app.borrow().document);
        assert_eq!(
            app.borrow().last_run_message.as_deref(),
            Some("Task preview printed to console.")
        );
        assert!(!app.borrow().last_run_is_error);
    }
}
