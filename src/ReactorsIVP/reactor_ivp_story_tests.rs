//! Story tests for the condensed-reactor IVP path.
//!
//! These tests exercise one canonical condensed-combustion fixture against
//! multiple solver-backend plans. The hypotheses are simple:
//! - the default typed plan should solve a small exothermic A -> B case;
//! - fixed BDF and fixed Adams should stay usable on a non-stiff story;
//! - the alternate symbolic assembly route should not break a supported matrix
//!   configuration.
//!
//! Expected outcome:
//! - solve status is accepted;
//! - the axis is monotonic and the solution matrix stays finite;
//! - the default backend remains Lambdify + AtomView + Sparse;
//! - the A -> B conversion trend survives backend changes.

#[cfg(test)]
mod tests {
    use super::super::SimpleReactorIVP::*;
    use super::super::solver_backend::{
        ReactorIvpExecutionBackend, ReactorIvpMatrixBackend, ReactorIvpMethod,
        ReactorIvpSolverConfig, ReactorIvpSymbolicBackend,
    };
    use crate::Kinetics::mechfinder_api::ReactionData;
    use approx::assert_abs_diff_eq;
    use std::collections::HashMap;

    fn build_canonical_story_task() -> SimpleReactorTask {
        let mut task = SimpleReactorTask::new();
        task.kindata
            .set_reaction_data_directly(
                vec![ReactionData::new_elementary(
                    "A=>B".to_string(),
                    vec![2.5e2, 0.0, 8.0e3],
                    None,
                )],
                None,
            )
            .expect("reaction setup should succeed");
        task.kindata.stecheodata.vec_of_molmasses = Some(vec![10.0, 20.0]);
        task.kindata.substances = vec!["A".to_string(), "B".to_string()];
        task.thermal_effects = vec![-1.5e5];
        task.ro = 1350.0;
        task.Cp = 1100.0;
        task.Lambda = 0.22;
        task.m = 0.012;
        task.scaling = crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig::new(80.0, 0.04, 90.0);
        task.initial_conditions = HashMap::from([
            ("T".to_string(), 520.0),
            ("q".to_string(), 0.05),
            ("A".to_string(), 0.85),
            ("B".to_string(), 0.15),
        ]);
        task
    }

    fn assert_story_solution(
        task: &SimpleReactorTask,
        snapshot: &IvpSolveSnapshot,
        expected_config: ReactorIvpSolverConfig,
    ) {
        assert_eq!(snapshot.backend_config(), expected_config);
        assert_eq!(snapshot.resolved_backend_plan(), expected_config);
        assert_eq!(snapshot.method(), expected_config.method);
        assert_eq!(snapshot.variable_count(), 4);
        assert_eq!(snapshot.summary.time_points, snapshot.time_points());
        assert_eq!(snapshot.summary.variable_count, snapshot.variable_count());
        assert!(snapshot.status().starts_with("finished"));
        assert!(snapshot.axis.iter().all(|value| value.is_finite()));
        assert!(snapshot.values.iter().all(|value| value.is_finite()));
        assert!(
            snapshot
                .axis
                .as_slice()
                .windows(2)
                .all(|pair| pair[1] >= pair[0])
        );
        assert_eq!(snapshot.method_label(), "bdf");
        assert_eq!(snapshot.values.ncols(), task.canonical_state_names().len());

        let initial_a = task.initial_conditions["A"];
        let initial_b = task.initial_conditions["B"];
        let final_row = snapshot.values.row(snapshot.values.nrows() - 1);

        assert!(
            (final_row[2] - initial_a).abs() > 1e-9,
            "A should change in the canonical story fixture"
        );
        assert!(
            (final_row[3] - initial_b).abs() > 1e-9,
            "B should change in the canonical story fixture"
        );

        assert_abs_diff_eq!(snapshot.axis[0], expected_config.x0, epsilon = 1e-12);
        assert_abs_diff_eq!(
            snapshot.axis[snapshot.axis.len() - 1],
            expected_config.x_bound,
            epsilon = 1e-10
        );
    }

    fn solve_story_with_config(
        config: ReactorIvpSolverConfig,
    ) -> (SimpleReactorTask, IvpSolveSnapshot) {
        let mut task = build_canonical_story_task();
        task.set_solver_backend_config(config);
        task.setup_condensed_ivp()
            .expect("canonical story setup should succeed");
        let snapshot = task.solve().expect("canonical story solve should succeed");
        (task, snapshot)
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_default_lambdify_sparse() {
        let config = ReactorIvpSolverConfig::default();
        let (task, snapshot) = solve_story_with_config(config);

        assert_eq!(
            task.solver_backend_config.execution_backend,
            ReactorIvpExecutionBackend::Lambdify
        );
        assert_eq!(
            task.solver_backend_config.symbolic_backend,
            ReactorIvpSymbolicBackend::AtomView
        );
        assert_eq!(
            task.solver_backend_config.matrix_backend,
            ReactorIvpMatrixBackend::Sparse
        );
        assert_eq!(task.solver_backend_config.method, ReactorIvpMethod::Auto);
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_fixed_bdf_sparse() {
        let config = ReactorIvpSolverConfig::default()
            .with_method(ReactorIvpMethod::Bdf)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let (task, snapshot) = solve_story_with_config(config);

        assert_eq!(task.solver_backend_config.method, ReactorIvpMethod::Bdf);
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_fixed_adams_sparse() {
        let config = ReactorIvpSolverConfig::default()
            .with_method(ReactorIvpMethod::Adams)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let (task, snapshot) = solve_story_with_config(config);

        assert_eq!(task.solver_backend_config.method, ReactorIvpMethod::Adams);
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_exprlegacy_and_banded_matrix() {
        let config = ReactorIvpSolverConfig::default()
            .with_method(ReactorIvpMethod::Auto)
            .with_symbolic_backend(ReactorIvpSymbolicBackend::ExprLegacy)
            .with_matrix_backend(ReactorIvpMatrixBackend::Banded)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let (task, snapshot) = solve_story_with_config(config);

        assert_eq!(
            task.solver_backend_config.symbolic_backend,
            ReactorIvpSymbolicBackend::ExprLegacy
        );
        assert!(task.solver_backend_config.matrix_backend.is_banded());
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_canonical_condensed_story_fixture_solves_with_inferred_banded_matrix() {
        let config = ReactorIvpSolverConfig::default()
            .with_method(ReactorIvpMethod::Auto)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_symbolic_backend(ReactorIvpSymbolicBackend::AtomView)
            .with_matrix_backend(ReactorIvpMatrixBackend::Banded)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let (task, snapshot) = solve_story_with_config(config);

        assert!(task.solver_backend_config.matrix_backend.is_banded());
        assert_story_solution(&task, &snapshot, config);
    }

    #[test]
    fn test_sparse_and_banded_routes_match_on_the_canonical_story_fixture() {
        let sparse_config = ReactorIvpSolverConfig::default()
            .with_method(ReactorIvpMethod::Auto)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_symbolic_backend(ReactorIvpSymbolicBackend::AtomView)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);
        let banded_config = ReactorIvpSolverConfig::default()
            .with_method(ReactorIvpMethod::Auto)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_symbolic_backend(ReactorIvpSymbolicBackend::AtomView)
            .with_matrix_backend(ReactorIvpMatrixBackend::Banded)
            .with_integration_domain(0.0, 0.02)
            .with_max_step(1.0e-3)
            .with_rtol(1.0e-6)
            .with_atol(1.0e-8);

        let (sparse_task, sparse_snapshot) = solve_story_with_config(sparse_config);
        let (banded_task, banded_snapshot) = solve_story_with_config(banded_config);

        assert_eq!(
            sparse_task.solver_backend_config.matrix_backend,
            ReactorIvpMatrixBackend::Sparse
        );
        assert!(banded_task.solver_backend_config.matrix_backend.is_banded());
        assert_eq!(sparse_snapshot.axis.len(), banded_snapshot.axis.len());
        assert_eq!(
            sparse_snapshot.values.ncols(),
            banded_snapshot.values.ncols()
        );
        for (lhs, rhs) in sparse_snapshot.axis.iter().zip(banded_snapshot.axis.iter()) {
            assert_abs_diff_eq!(lhs, rhs, epsilon = 1e-10);
        }
        for (lhs, rhs) in sparse_snapshot
            .values
            .row(sparse_snapshot.values.nrows() - 1)
            .iter()
            .zip(
                banded_snapshot
                    .values
                    .row(banded_snapshot.values.nrows() - 1)
                    .iter(),
            )
        {
            let scale = lhs.abs().max(rhs.abs()).max(1.0);
            let rel = (lhs - rhs).abs() / scale;
            assert!(
                rel <= 1.0e-8,
                "Sparse and banded routes diverged too much on the canonical story fixture: lhs={lhs:e}, rhs={rhs:e}, rel={rel:e}"
            );
        }
    }

    #[test]
    fn test_canonical_story_fixture_reuses_initial_state_shape_after_rerun() {
        let mut task = build_canonical_story_task();
        task.setup_condensed_ivp()
            .expect("canonical story setup should succeed");
        let first = task.solve().expect("first solve should succeed");

        let override_config = ReactorIvpSolverConfig::default()
            .with_method(ReactorIvpMethod::Bdf)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 1.0);
        let second = task
            .solve_with_config(override_config)
            .expect("second solve should succeed");

        assert_eq!(first.variable_count(), second.variable_count());
        assert!(second.time_points() > 0);
        assert_eq!(
            task.last_solve_snapshot.as_ref().unwrap().status(),
            second.status()
        );
    }
}
