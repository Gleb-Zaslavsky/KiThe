#[cfg(test)]
mod tests {
    use super::super::SimpleReactorIVP::*;
    use super::super::solver_backend::{
        ReactorIvpExecutionBackend, ReactorIvpMatrixBackend, ReactorIvpMethod,
        ReactorIvpSolverConfig,
    };
    use crate::Kinetics::mechfinder_api::ReactionData;
    use RustedSciThe::numerical::LSODE2::Lsode2AlgorithmSnapshot;
    use RustedSciThe::numerical::LSODE2::solver::{
        Lsode2EvaluationTelemetry, Lsode2TelemetryScope,
    };
    use approx::assert_abs_diff_eq;
    use nalgebra::{DMatrix, DVector};
    use std::collections::HashMap;

    fn build_solve_task() -> SimpleReactorTask {
        let mut task = SimpleReactorTask::new();
        task.kindata
            .set_reaction_data_directly(
                vec![ReactionData::new_elementary(
                    "A=>B".to_string(),
                    vec![1.0e2, 0.0, 1.0e4],
                    None,
                )],
                None,
            )
            .expect("reaction setup should succeed");
        task.kindata.stecheodata.vec_of_molmasses = Some(vec![10.0, 20.0]);
        task.kindata.substances = vec!["A".to_string(), "B".to_string()];
        task.thermal_effects = vec![-2.0e5];
        task.ro = 1200.0;
        task.Cp = 1000.0;
        task.Lambda = 0.25;
        task.m = 0.015;
        task.scaling =
            crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig::new(100.0, 0.05, 100.0);
        task.initial_conditions = HashMap::from([
            ("T".to_string(), 450.0),
            ("q".to_string(), 0.15),
            ("A".to_string(), 0.7),
            ("B".to_string(), 0.3),
        ]);
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default()
                .with_method(ReactorIvpMethod::Auto)
                .with_native_execution(
                    RustedSciThe::numerical::LSODE2::Lsode2NativeExecutionConfig::bridge_solve(),
                )
                .with_integration_domain(0.0, 0.01)
                .with_max_step(1.0e-3)
                .with_rtol(1.0e-6)
                .with_atol(1.0e-8),
        );
        task
    }

    fn build_zero_rate_task() -> SimpleReactorTask {
        let mut task = SimpleReactorTask::new();
        task.kindata
            .set_reaction_data_directly(
                vec![ReactionData::new_elementary(
                    "A=>B".to_string(),
                    vec![0.0, 0.0, 0.0],
                    None,
                )],
                None,
            )
            .expect("reaction setup should succeed");
        task.kindata.stecheodata.vec_of_molmasses = Some(vec![10.0, 20.0]);
        task.kindata.substances = vec!["A".to_string(), "B".to_string()];
        task.thermal_effects = vec![0.0];
        task.ro = 1200.0;
        task.Cp = 1000.0;
        task.Lambda = 0.25;
        task.m = 0.015;
        task.scaling =
            crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig::new(100.0, 0.05, 100.0);
        task.initial_conditions = HashMap::from([
            ("T".to_string(), 450.0),
            ("q".to_string(), 0.15),
            ("A".to_string(), 0.7),
            ("B".to_string(), 0.3),
        ]);
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default()
                .with_method(ReactorIvpMethod::Auto)
                .with_native_execution(
                    RustedSciThe::numerical::LSODE2::Lsode2NativeExecutionConfig::bridge_solve(),
                )
                .with_integration_domain(0.0, 0.01)
                .with_max_step(1.0e-3)
                .with_rtol(1.0e-6)
                .with_atol(1.0e-8),
        );
        task
    }

    fn assert_rhs_is_finite(task: &SimpleReactorTask) {
        let config = task
            .build_lsode2_problem_config()
            .expect("problem config should build");
        let names = config
            .values
            .iter()
            .map(|name| name.as_str())
            .collect::<Vec<_>>();
        let y0 = config.y0.as_slice();
        // The condensed IVP path must publish a numerically sane RHS snapshot.
        for (index, expr) in config.eq_system.iter().enumerate() {
            let ir = expr.lower_to_linear(&names);
            let value = ir.eval(y0);
            assert!(
                value.is_finite(),
                "RHS expression at index {} evaluated to a non-finite value: {}",
                index,
                value
            );
        }
    }

    fn fake_summary(
        status: &str,
        time_points: usize,
        variable_count: usize,
        final_t: Option<f64>,
        final_y: Option<DVector<f64>>,
    ) -> RustedSciThe::numerical::LSODE2::Lsode2SolveSummary {
        RustedSciThe::numerical::LSODE2::Lsode2SolveSummary {
            method: "bdf",
            jacobian_backend: "symbolic",
            linear_solver_backend: "banded",
            linear_solver_reason: "test",
            resolved_source: "symbolic:lambdify",
            resolved_structure: "sparse",
            status: status.to_string(),
            time_points,
            variable_count,
            final_t,
            final_y,
            max_abs_solution: 1.0,
            algorithm: Lsode2AlgorithmSnapshot {
                controller_mode: "automatic_adams_bdf",
                active_family: "bdf",
                mused_family: "bdf",
                mcur_family: "bdf",
                preferred_family: "bdf",
                executed_family: Some("bdf"),
                switch_reason: "test",
                switch_uses_fallback: false,
                method_switching_enabled: false,
                max_adams_order: 12,
                max_bdf_order: 5,
                stiffness_ratio_threshold: 1.0,
                method_switch_probe_steps: 0,
                switch_probe_countdown: 0,
                switch_probe_ready: false,
                tsw: None,
                last_handoff_jstart: None,
                bdf_current_order: None,
                bdf_max_order_cap: None,
                bdf_equal_step_count: None,
                note: "test",
            },
            statistics: Default::default(),
            native_statistics: Default::default(),
            evaluation_telemetry: Lsode2EvaluationTelemetry {
                scope: Lsode2TelemetryScope::None,
                ..Default::default()
            },
            native_step_probe: None,
            native_integration_preview: None,
            native_integration_solve: None,
            native_termination_kind: None,
        }
    }

    #[test]
    fn test_condensed_ivp_solve_publishes_valid_snapshot() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        assert_rhs_is_finite(&task);
        let snapshot = task.solve().expect("condensed solve should succeed");

        assert_eq!(snapshot.time_points(), snapshot.axis.len());
        assert_eq!(snapshot.variable_count(), 4);
        assert_eq!(snapshot.summary.time_points, snapshot.time_points());
        assert_eq!(snapshot.summary.variable_count, snapshot.variable_count());
        assert!(!snapshot.status().trim().is_empty());
        assert!(snapshot.status().starts_with("finished"));
        assert!(snapshot.axis.iter().all(|value| value.is_finite()));
        assert!(snapshot.values.iter().all(|value| value.is_finite()));
        assert_eq!(snapshot.backend_config(), task.solver_backend_config);
        assert_eq!(snapshot.resolved_backend_plan(), task.solver_backend_config);
        assert_eq!(snapshot.method(), ReactorIvpMethod::Auto);
        assert_eq!(snapshot.method_label(), "bdf");
        assert!(
            !snapshot
                .algorithm_snapshot()
                .controller_mode
                .trim()
                .is_empty()
        );
        assert!(snapshot.backend_statistics().solve_calls > 0);
        assert!(snapshot.native_statistics().native_step_attempts > 0);

        let axis = task
            .solver
            .x_mesh
            .as_ref()
            .expect("x mesh should be stored");
        let values = task
            .solver
            .solution
            .as_ref()
            .expect("solution should be stored");
        assert_eq!(axis.len(), snapshot.time_points());
        assert_eq!(values.nrows(), snapshot.time_points());
        assert_eq!(values.ncols(), snapshot.variable_count());
        assert_abs_diff_eq!(axis[0], snapshot.axis[0] * task.L, epsilon = 1e-12);
        assert_abs_diff_eq!(
            axis[axis.len() - 1],
            snapshot.axis[snapshot.axis.len() - 1] * task.L,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            values[(0, 0)],
            snapshot.values[(0, 0)] * task.scaling.T_scale + task.scaling.dT,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            values[(0, 1)],
            snapshot.values[(0, 1)] * task.scaling.T_scale / task.L,
            epsilon = 1e-12
        );
        assert!(task.last_solve_snapshot.is_some());
        assert_eq!(
            task.last_solve_snapshot.as_ref().unwrap().status(),
            snapshot.status()
        );

        let conservation = task
            .latest_conservation_report()
            .expect("conservation report should be published after a successful solve");
        assert!(conservation.energy_balance_error_abs.is_finite());
        assert!(conservation.energy_balance_error_rel.is_finite());
        assert!(
            conservation
                .sum_of_mass_fractions
                .iter()
                .all(|(_, value)| value.is_finite())
        );
        assert!(
            conservation
                .atomic_mass_balance_error
                .iter()
                .all(|(_, value)| value.is_finite())
        );
        assert_eq!(task.last_conservation_report.as_ref(), Some(conservation));
    }

    #[test]
    fn test_latest_solve_report_exposes_diagnostics_without_touching_solution_matrix() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        let snapshot = task.solve().expect("condensed solve should succeed");

        let report = task
            .latest_solve_report()
            .expect("latest solve report should exist");

        assert_eq!(report.status, snapshot.status());
        assert_eq!(report.method_label, snapshot.method_label());
        assert_eq!(report.time_points, snapshot.time_points());
        assert_eq!(report.variable_count, snapshot.variable_count());
        assert_eq!(report.backend_plan, snapshot.backend_config());
        assert_eq!(report.algorithm, *snapshot.algorithm_snapshot());
        assert_eq!(
            report.statistics.table_report(),
            snapshot.backend_statistics().table_report()
        );
        assert_eq!(report.native_statistics, *snapshot.native_statistics());
        assert_eq!(
            report.evaluation_telemetry,
            *snapshot.evaluation_telemetry()
        );
        assert_eq!(
            task.last_solve_snapshot.as_ref().map(|s| s.status()),
            Some(snapshot.status())
        );
    }

    #[test]
    fn test_solution_preview_rows_expose_data_without_printing() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        let snapshot = task.solve().expect("condensed solve should succeed");
        let rendered = task
            .solution_render_data()
            .expect("rendered solution should exist after solve");

        let rows = task
            .solver
            .solution_preview_rows()
            .expect("solution preview should exist");

        assert_eq!(rows.len(), snapshot.variable_count());
        assert_eq!(rows[0].index, 0);
        assert_eq!(rows[0].variable, "T");
        assert_abs_diff_eq!(rows[0].first, rendered.solution[(0, 0)], epsilon = 1e-12);
        assert_abs_diff_eq!(
            rows[0].middle,
            rendered.solution[(rendered.solution.nrows() / 2, 0)],
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            rows[0].last,
            rendered.solution[(rendered.solution.nrows() - 1, 0)],
            epsilon = 1e-12
        );
        assert!(rows.iter().all(|row| row.first.is_finite()));
        assert!(rows.iter().all(|row| row.middle.is_finite()));
        assert!(rows.iter().all(|row| row.last.is_finite()));
    }

    #[test]
    fn test_snapshot_solution_preview_rows_match_borrowed_view() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        let snapshot = task.solve().expect("condensed solve should succeed");

        let snapshot_rows = snapshot
            .solution_preview_rows()
            .expect("snapshot preview rows should exist");
        let borrowed_rows = snapshot
            .result_view()
            .solution_preview_rows()
            .expect("borrowed preview rows should exist");

        assert_eq!(snapshot_rows, borrowed_rows);
    }

    #[test]
    fn test_latest_solve_report_rows_expose_table_ready_diagnostics() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        let snapshot = task.solve().expect("condensed solve should succeed");

        let rows = task
            .latest_solve_report_rows()
            .expect("latest solve report rows should exist");

        assert!(rows.iter().any(|row| row.field == "status"));
        assert!(rows.iter().any(|row| row.field == "method_label"));
        assert!(rows.iter().any(|row| row.field == "controller_mode"));
        assert!(rows.iter().any(|row| row.field == "resolved_source"));
        assert!(rows.iter().any(|row| row.field == "resolved_structure"));
        assert!(rows.iter().any(|row| row.field == "accepted_steps"));

        let status = rows
            .iter()
            .find(|row| row.field == "status")
            .map(|row| row.value.as_str());
        assert_eq!(status, Some(snapshot.status()));

        let method_label = rows
            .iter()
            .find(|row| row.field == "method_label")
            .map(|row| row.value.as_str());
        assert_eq!(method_label, Some(snapshot.method_label()));
    }

    #[test]
    fn test_latest_solve_report_rows_reject_missing_snapshot() {
        let task = build_solve_task();
        let result = task.latest_solve_report_rows();
        assert!(matches!(result, Err(IvpError::MissingData(_))));
    }

    #[test]
    fn test_latest_result_view_exposes_solution_preview_rows() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        let snapshot = task.solve().expect("condensed solve should succeed");

        let view = task
            .latest_result_view()
            .expect("latest result view should exist");
        let rows = view
            .solution_preview_rows()
            .expect("preview rows should be available");

        assert_eq!(rows.len(), snapshot.variable_count());
        assert_eq!(rows[0].variable, "Teta");
        assert_abs_diff_eq!(rows[0].first, snapshot.values[(0, 0)], epsilon = 1e-12);
        assert_abs_diff_eq!(
            rows[0].middle,
            snapshot.values[(snapshot.values.nrows() / 2, 0)],
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            rows[0].last,
            snapshot.values[(snapshot.values.nrows() - 1, 0)],
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_latest_result_view_matches_latest_report_rows() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        task.solve().expect("condensed solve should succeed");

        let view = task
            .latest_result_view()
            .expect("latest result view should exist");
        let report_rows_from_view = view.report_rows();
        let report_rows_from_task = task
            .latest_solve_report_rows()
            .expect("latest solve report rows should exist");

        assert_eq!(report_rows_from_view, report_rows_from_task);
    }

    #[test]
    fn test_result_view_rejects_variable_name_mismatch() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        let mut snapshot = task.solve().expect("condensed solve should succeed");
        snapshot.variable_names.push("extra_variable".to_string());

        let view = snapshot.result_view();
        let result = view.solution_preview_rows();
        assert!(matches!(result, Err(IvpError::CalculationError(_))));
    }

    #[test]
    fn test_snapshot_solution_preview_rows_reject_empty_matrix() {
        let snapshot = IvpSolveSnapshot {
            axis: DVector::from_vec(vec![]),
            values: DMatrix::zeros(0, 4),
            variable_names: vec![
                "Teta".to_string(),
                "q".to_string(),
                "C0".to_string(),
                "C1".to_string(),
            ],
            summary: fake_summary("finished", 0, 4, None, None),
            backend_plan: ReactorIvpSolverConfig::default(),
            algorithm: fake_summary("finished", 0, 4, None, None).algorithm,
            statistics: Default::default(),
            native_statistics: Default::default(),
            evaluation_telemetry: Default::default(),
        };

        let result = snapshot.solution_preview_rows();
        assert!(matches!(result, Err(IvpError::MissingData(_))));
    }

    #[test]
    fn test_snapshot_solution_preview_rows_support_one_point_output() {
        let snapshot = IvpSolveSnapshot {
            axis: DVector::from_vec(vec![0.0]),
            values: DMatrix::from_row_slice(1, 4, &[10.0, 20.0, 0.7, 0.3]),
            variable_names: vec![
                "Teta".to_string(),
                "q".to_string(),
                "C0".to_string(),
                "C1".to_string(),
            ],
            summary: fake_summary(
                "finished",
                1,
                4,
                Some(0.0),
                Some(DVector::from_vec(vec![10.0, 20.0, 0.7, 0.3])),
            ),
            backend_plan: ReactorIvpSolverConfig::default(),
            algorithm: fake_summary(
                "finished",
                1,
                4,
                Some(0.0),
                Some(DVector::from_vec(vec![10.0, 20.0, 0.7, 0.3])),
            )
            .algorithm,
            statistics: Default::default(),
            native_statistics: Default::default(),
            evaluation_telemetry: Default::default(),
        };

        let rows = snapshot
            .solution_preview_rows()
            .expect("one-point preview should exist");

        assert_eq!(rows.len(), 4);
        assert_eq!(rows[0].first, 10.0);
        assert_eq!(rows[0].middle, 10.0);
        assert_eq!(rows[0].last, 10.0);
        assert_eq!(rows[1].first, 20.0);
        assert_eq!(rows[2].first, 0.7);
        assert_eq!(rows[3].first, 0.3);
    }

    #[test]
    fn test_latest_result_view_rejects_missing_snapshot() {
        let task = build_solve_task();
        let result = task.latest_result_view();
        assert!(matches!(result, Err(IvpError::MissingData(_))));
    }

    #[test]
    fn test_condensed_ivp_solve_keeps_canonical_snapshot_shape() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        assert_rhs_is_finite(&task);
        let snapshot = task
            .solve_condensed_ivp()
            .expect("condensed solve should succeed");

        let names = task.canonical_state_names();
        assert_eq!(names, vec!["Teta", "q", "C0", "C1"]);
        assert_eq!(snapshot.values.ncols(), names.len());
        assert_eq!(snapshot.values.nrows(), snapshot.axis.len());
        assert_eq!(snapshot.backend_config(), task.solver_backend_config);
        assert_eq!(snapshot.resolved_backend_plan(), task.solver_backend_config);
        assert!(
            snapshot
                .axis
                .iter()
                .zip(snapshot.axis.iter().skip(1))
                .all(|(current, next)| next >= current)
        );
        assert!(snapshot.evaluation_telemetry().accepted_steps > 0);
    }

    #[test]
    fn test_condensed_ivp_solve_with_config_uses_temporary_backend_override() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let override_config = ReactorIvpSolverConfig::default()
            .with_method(ReactorIvpMethod::Adams)
            .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
            .with_matrix_backend(ReactorIvpMatrixBackend::Sparse)
            .with_integration_domain(0.0, 0.02);

        let snapshot = task
            .solve_with_config(override_config)
            .expect("temporary backend override should solve");

        assert_eq!(snapshot.backend_config(), override_config);
        assert_eq!(snapshot.method(), ReactorIvpMethod::Adams);
        assert_eq!(task.solver_backend_config, override_config);
    }

    #[test]
    fn test_validate_condensed_solve_snapshot_rejects_nonfinite_values() {
        let mut task = build_solve_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let axis = DVector::from_vec(vec![0.0, f64::NAN]);
        let values =
            DMatrix::from_row_slice(2, 4, &[0.0, 0.0, 0.0, 0.0, 1.0, 2.0, f64::INFINITY, 4.0]);
        let summary = fake_summary(
            "finished",
            axis.len(),
            values.ncols(),
            Some(1.0),
            Some(DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0])),
        );

        let result = task.validate_condensed_solve_snapshot(&axis, &values, &summary);

        assert!(matches!(
            result,
            Err(IvpError::CalculationError(message))
            if message.contains("not finite")
        ));
    }

    #[test]
    fn test_zero_rate_condensed_ivp_preserves_axis_and_matrix_orientation() {
        let mut task = build_zero_rate_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let snapshot = task
            .solve()
            .expect("zero-rate condensed solve should succeed");

        assert_eq!(snapshot.axis.len(), snapshot.values.nrows());
        assert_eq!(snapshot.values.ncols(), task.canonical_state_names().len());
        assert!(
            snapshot
                .axis
                .as_slice()
                .windows(2)
                .all(|pair| pair[1] >= pair[0])
        );
        assert!(snapshot.values.iter().all(|value| value.is_finite()));
        assert_eq!(
            snapshot.values.row(0).len(),
            task.canonical_state_names().len()
        );
        assert_eq!(
            snapshot.values.row(snapshot.values.nrows() - 1).len(),
            task.canonical_state_names().len()
        );
    }
}
