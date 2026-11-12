use super::*;
use strum::IntoEnumIterator;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parameter_free_models() {
        let models = [
            KineticModelNames::A2,
            KineticModelNames::A3,
            KineticModelNames::A4,
            KineticModelNames::R2,
            KineticModelNames::R3,
            KineticModelNames::D1,
            KineticModelNames::D2,
            KineticModelNames::D3,
            KineticModelNames::D4,
            KineticModelNames::P2_3,
            KineticModelNames::P2,
            KineticModelNames::P3,
            KineticModelNames::F1,
            KineticModelNames::F2,
            KineticModelNames::F3,
        ];

        for model_name in models {
            let result = create_kinetic_model(model_name.clone(), vec![]);
            assert!(result.is_ok(), "Failed to create model {:?}", model_name);

            let (kinetic_model, jacobian) = result.unwrap();
            let expr = kinetic_model.get_expr();

            // Verify expressions are not empty/trivial
            assert_ne!(format!("{:?}", expr), "Const(0.0)");
            assert_ne!(format!("{:?}", jacobian), "Const(0.0)");
        }
    }

    #[test]
    fn test_parameterized_models() {
        // Test SB model
        let result = create_kinetic_model(KineticModelNames::SB, vec![1.0, 2.0, 0.5]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        let sbjac = SB_j(1.0, 2.0, 0.5);
        assert!(matches!(model, KineticModel::SB(1.0, 2.0, 0.5)));
        assert_eq!(jac, sbjac);

        // Test JMA model
        let result = create_kinetic_model(KineticModelNames::JMA, vec![2.0]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        let jmajac = JMA_jac(2.0);
        assert!(matches!(model, KineticModel::JMA(2.0)));
        assert_eq!(jac, jmajac);

        // Test Ac model
        let result = create_kinetic_model(KineticModelNames::Ac, vec![1.5]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        assert!(matches!(model, KineticModel::Ac(1.5)));
        let acjac = Ac_jac(1.5);
        assert_eq!(jac, acjac);

        // Test Dec model
        let result = create_kinetic_model(KineticModelNames::Dec, vec![3.0]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        let decjac = Dec_jac(3.0);
        assert!(matches!(model, KineticModel::Dec(3.0)));
        assert_eq!(jac, decjac);

        // Test PTe model
        let result = create_kinetic_model(KineticModelNames::PTe, vec![2.0, 1.0]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        let ptejac = PTe_jac(2.0, 1.0);
        assert!(matches!(model, KineticModel::PTe(2.0, 1.0)));
        assert_eq!(jac, ptejac);

        // Test SBtp model
        let result = create_kinetic_model(KineticModelNames::SBtp, vec![1.0, 2.0, 0.8]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        assert!(matches!(model, KineticModel::SBtp(1.0, 2.0, 0.8)));
        let sbtpjac = SBtp_jac(1.0, 2.0, 0.8);
        assert_eq!(jac, sbtpjac);
    }

    #[test]

    fn test_parameter_validation() {
        // Test wrong number of parameters
        assert!(create_kinetic_model(KineticModelNames::SB, vec![1.0]).is_err());
        assert!(create_kinetic_model(KineticModelNames::JMA, vec![]).is_err());
        assert!(create_kinetic_model(KineticModelNames::PTe, vec![1.0]).is_err());
        assert!(create_kinetic_model(KineticModelNames::SBtp, vec![1.0, 2.0]).is_err());
    }

    #[test]
    fn test_jacobian_consistency() {
        // Test that jacobians are properly generated for parameterized models
        let test_cases = [
            (KineticModelNames::SB, vec![1.0, 2.0, 0.5]),
            (KineticModelNames::JMA, vec![2.0]),
            (KineticModelNames::Ac, vec![1.5]),
            (KineticModelNames::Dec, vec![3.0]),
            (KineticModelNames::PTe, vec![2.0, 1.0]),
            (KineticModelNames::SBtp, vec![1.0, 2.0, 0.8]),
        ];

        for (model_name, params) in test_cases {
            let result = create_kinetic_model(model_name.clone(), params);
            assert!(result.is_ok(), "Failed to create model {:?}", model_name);

            let (kinetic_model, jacobian) = result.unwrap();

            // Verify jacobian is not trivial
            assert_ne!(format!("{:?}", jacobian), "Const(0.0)");

            // For models with explicit jacobian functions, verify they're used
            match kinetic_model {
                KineticModel::SB(_, _, _) => {
                    // Should use SB_j function, not automatic differentiation
                    assert!(
                        format!("{:?}", jacobian).contains("Pow")
                            || format!("{:?}", jacobian).contains("ln")
                    );
                }
                KineticModel::JMA(_) => {
                    assert!(
                        format!("{:?}", jacobian).contains("Pow")
                            || format!("{:?}", jacobian).contains("ln")
                    );
                }
                _ => {}
            }
        }
    }

    #[test]
    fn test_individual_functions() {
        // Test A2 function
        let expr = A2();
        println!("{:?}", expr);
        assert!(format!("{:?}", expr).contains("Pow"));
        assert!(format!("{:?}", expr).contains("Ln"));

        // Test parameterized functions
        let sb_expr = Sestak_Berggen(1.0, 2.0, 0.5);
        assert!(format!("{:?}", sb_expr).contains("Pow"));

        let jma_expr = Johnson_Mehl_Avrami(2.0);
        assert!(format!("{:?}", jma_expr).contains("Pow"));
        assert!(format!("{:?}", jma_expr).contains("Ln"));
    }

    #[test]
    fn test_model_names_mapping() {
        // Test that all model names have proper string mappings
        for model in KineticModelNames::iter() {
            let mapping = model.map_of_names_and_formulas();
            assert!(!mapping.is_empty(), "Empty mapping for {:?}", model);
            assert!(
                mapping.contains("a"),
                "Mapping should contain variable 'a' for {:?}",
                model
            );
        }
    }

    #[test]
    fn test_expression_evaluation() {
        // Test that expressions can be created and contain expected elements
        let models_with_params = [
            (KineticModelNames::A2, vec![]),
            (KineticModelNames::F1, vec![]),
            (KineticModelNames::SB, vec![1.0, 2.0, 0.5]),
            (KineticModelNames::JMA, vec![2.0]),
        ];

        for (model_name, params) in models_with_params {
            let result = create_kinetic_model(model_name.clone(), params);
            assert!(result.is_ok());

            let (kinetic_model, jacobian) = result.unwrap();
            let expr = kinetic_model.get_expr();

            // Verify expressions contain variable 'a'
            let expr_str = format!("{:?}", expr);
            let _jac_str = format!("{:?}", jacobian);

            assert!(
                expr_str.contains("Var(\"a\")"),
                "Expression should contain variable 'a' for {:?}",
                model_name
            );
            // Jacobian might be constant for some simple models, so we don't always require 'a'
        }
    }

    #[test]
    fn test_constants() {
        // Test that constants are properly defined
        assert_eq!(two, Expr::Const(2.0));
        assert_eq!(three, Expr::Const(3.0));
        assert_eq!(four, Expr::Const(4.0));
        assert_eq!(one, Expr::Const(1.0));
        assert_eq!(half, Expr::Const(0.5));
        assert_eq!(neg, Expr::Const(-1.0));
    }

    #[test]
    fn test_create_kinetic_model_error_handling() {
        // Test error handling in create_kinetic_model
        let result = create_kinetic_model(KineticModelNames::SB, vec![1.0, 2.0]); // Missing one parameter
        assert!(result.is_err());

        let result = create_kinetic_model(KineticModelNames::JMA, vec![1.0, 2.0]); // Too many parameters
        assert!(result.is_err());
    }
}

#[cfg(test)]
mod kinetic_model_ivp_tests {
    use super::*;
    use RustedSciThe::numerical::StiffAPI::StiffSolverKind as SolverType;
    use std::time::Instant;

    #[test]
    fn test_new() {
        let ivp = KineticModelIVP::new(SolverType::Bdf);

        assert!(ivp.solver.is_none());
        assert_eq!(ivp.t_final, 0.0);
        assert_eq!(ivp.beta, 0.0);
        assert_eq!(ivp.T0, 0.0);
        assert_eq!(ivp.E, 0.0);
        assert_eq!(ivp.A, 0.0);
        assert!(ivp.kinmodel.is_none());
        assert!(matches!(ivp.solvertype, SolverType::Bdf));

        assert_eq!(ivp.solver_options, SolidStateIvpSolverOptions::default());
    }

    #[test]
    fn test_new_lsode2_uses_lsode2_solver_type() {
        let ivp = KineticModelIVP::new_lsode2();

        assert!(matches!(ivp.solvertype, SolverType::Lsode2));
        assert!(ivp.solver.is_none());
    }

    #[test]
    fn test_default_uses_lsode2_solver_type() {
        let ivp = KineticModelIVP::default();

        assert!(matches!(ivp.solvertype, SolverType::Lsode2));
        assert!(ivp.solver.is_none());
    }

    #[test]
    fn test_set_problem() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);

        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        assert_eq!(ivp.t_final, 100.0);
        assert_eq!(ivp.beta, 10.0);
        assert_eq!(ivp.T0, 298.15);
        assert_eq!(ivp.E, 50000.0);
        assert_eq!(ivp.A, 1e6);
    }

    #[test]
    fn test_set_model() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        // Test parameter-free model
        let _ = ivp.set_model(KineticModelNames::A2, vec![]);

        assert!(ivp.kinmodel.is_some());
        assert!(matches!(ivp.kinmodel.as_ref().unwrap(), KineticModel::A2));

        // Test parameterized model
        let _ = ivp.set_model(KineticModelNames::JMA, vec![2.0]);

        assert!(matches!(
            ivp.kinmodel.as_ref().unwrap(),
            KineticModel::JMA(2.0)
        ));
        assert!(ivp.kinetic_expression().is_ok());
    }

    #[test]
    fn test_get_solution_vectors() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        ivp.set_problem(10.0, 5.0, 298.15, 50000.0, 1e6).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();
        ivp.check_task().unwrap();
        ivp.solve().unwrap();

        let sol = ivp.get_solution().unwrap();
        let n = sol.time.len();

        assert!(n > 2);
        assert_eq!(sol.temperature.len(), n);
        assert_eq!(sol.conversion.len(), n);
        assert_eq!(sol.conversion_rate.len(), n);

        for i in 0..n {
            let expected_t = ivp.T0 + ivp.beta * sol.time[i];
            assert!((sol.temperature[i] - expected_t).abs() < 1e-9);
            let arrhenius = ivp.A * (-ivp.E / (8.31446261815324 * expected_t)).exp();
            let expected_rate = arrhenius * (1.0 - sol.conversion[i]);
            assert!(
                (sol.conversion_rate[i] - expected_rate).abs()
                    <= 1e-12 * expected_rate.abs().max(1.0)
            );
        }
    }

    #[test]
    fn test_lsode2_smoke_solve_for_f1() {
        let mut ivp = KineticModelIVP::new_lsode2();
        ivp.set_problem(10.0, 5.0, 298.15, 50000.0, 1e6).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();
        ivp.check_task().unwrap();
        ivp.solve().unwrap();

        let sol = ivp.get_solution().unwrap();
        assert!(!sol.time.is_empty());
        assert_eq!(sol.temperature.len(), sol.time.len());
        assert_eq!(sol.conversion.len(), sol.time.len());
        assert_eq!(sol.conversion_rate.len(), sol.time.len());
    }

    /// Interpolates a solution value at `t` on a monotonic time grid.
    fn interpolate_solution(times: &[f64], values: &[f64], t: f64) -> Option<f64> {
        if times.len() != values.len() || times.len() < 2 {
            return None;
        }

        if t <= times[0] {
            return Some(values[0]);
        }

        let last_time = *times.last()?;
        if t >= last_time {
            return Some(*values.last()?);
        }

        let upper = times.partition_point(|&time| time < t);
        if upper == 0 || upper >= times.len() {
            return None;
        }

        let lower = upper - 1;
        let t0 = times[lower];
        let t1 = times[upper];
        let y0 = values[lower];
        let y1 = values[upper];
        let span = t1 - t0;
        if span.abs() < f64::EPSILON {
            return Some(y1);
        }

        let weight = (t - t0) / span;
        Some(y0 + weight * (y1 - y0))
    }

    #[test]
    fn test_lsode2_matches_bdf_for_f1() {
        let mut lsode2 = KineticModelIVP::new_lsode2();
        let mut bdf = KineticModelIVP::new(SolverType::Bdf);

        for ivp in [&mut lsode2, &mut bdf] {
            ivp.set_problem(10.0, 5.0, 298.15, 50000.0, 1e6).unwrap();
            ivp.set_model(KineticModelNames::F1, vec![]).unwrap();
            ivp.check_task().unwrap();
            ivp.solve().unwrap();
        }

        let lsode2_sol = lsode2.get_solution().unwrap();
        let bdf_sol = bdf.get_solution().unwrap();

        let shared_end = lsode2_sol
            .time
            .last()
            .copied()
            .unwrap_or_default()
            .min(bdf_sol.time.last().copied().unwrap_or_default());
        let sample_count = 10usize;
        let mut max_abs_diff = 0.0_f64;

        for sample_idx in 0..sample_count {
            let fraction = (sample_idx + 1) as f64 / sample_count as f64;
            let t = shared_end * fraction;
            let lsode2_value = interpolate_solution(&lsode2_sol.time, &lsode2_sol.conversion, t)
                .expect("LSODE2 solution should interpolate on the shared grid");
            let bdf_value = interpolate_solution(&bdf_sol.time, &bdf_sol.conversion, t)
                .expect("BDF solution should interpolate on the shared grid");
            max_abs_diff = max_abs_diff.max((lsode2_value - bdf_value).abs());
        }

        // Adaptive solvers may land on slightly different trajectories while still
        // agreeing on the chemistry in practice, so keep the envelope modest.
        assert!(
            max_abs_diff < 1e-2,
            "LSODE2 and BDF diverged too much: max abs diff = {}",
            max_abs_diff
        );
    }

    #[test]
    fn stiffapi_solve_1d_dispatches_all_supported_scalar_solver_families() {
        use RustedSciThe::numerical::StiffAPI::StiffSolverKind;

        let cases = [
            StiffSolverKind::Bdf,
            StiffSolverKind::Radau,
            StiffSolverKind::BackwardEuler,
            StiffSolverKind::Lsode2,
        ];
        let mut endpoints = Vec::new();

        for selection in cases {
            let mut ivp = KineticModelIVP::new(selection);
            ivp.set_problem(4.0, 2.0, 300.0, 35_000.0, 2.0e5).unwrap();
            ivp.set_model(KineticModelNames::F1, vec![]).unwrap();
            ivp.solve().unwrap();

            assert_eq!(ivp.solver.as_ref().unwrap().kind(), selection);
            let (time, state) = ivp.get_result().unwrap();
            assert_eq!(state.ncols(), 1);
            assert_eq!(state.nrows(), time.len());
            assert!(time.len() > 1);
            assert!(state.iter().all(|value| value.is_finite()));
            let final_conversion = state[(state.nrows() - 1, 0)];
            assert!((0.0..=1.0).contains(&final_conversion));
            endpoints.push(final_conversion);
        }

        let reference = endpoints[0];
        for value in endpoints.iter().skip(1) {
            assert!(
                (value - reference).abs() < 2e-2,
                "scalar stiff solver routes disagree: {endpoints:?}"
            );
        }
    }

    #[test]
    fn scalar_telemetry_is_off_by_default_and_native_when_enabled() {
        use RustedSciThe::numerical::SolverTelemetryMode;
        use RustedSciThe::numerical::StiffAPI::StiffTelemetrySnapshot;

        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        ivp.set_problem(2.0, 1.0, 300.0, 35_000.0, 2.0e5).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();
        ivp.solve().unwrap();
        assert!(ivp.telemetry_snapshot().is_none());

        let mut options = SolidStateIvpSolverOptions::default();
        options.telemetry = SolverTelemetryMode::Counters;
        ivp.set_solver_options(options);
        ivp.solve().unwrap();
        assert!(matches!(
            ivp.telemetry_snapshot(),
            Some(StiffTelemetrySnapshot::Bdf {
                mode: SolverTelemetryMode::Counters,
                ..
            })
        ));
    }

    #[test]
    fn invalid_typed_solver_options_fail_before_preparation() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        ivp.set_problem(2.0, 1.0, 300.0, 35_000.0, 2.0e5).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();

        let mut options = SolidStateIvpSolverOptions::default();
        options.max_steps = Some(0);
        ivp.set_solver_options(options);
        let error = ivp.solve().unwrap_err();
        assert!(error.contains("max_steps must be positive when specified"));
        assert!(ivp.get_result().is_err());
    }

    #[test]
    fn typed_solver_options_are_applied_directly() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        let options = SolidStateIvpSolverOptions {
            max_step: 0.025,
            rtol: 1e-5,
            max_steps: Some(2000),
            ..SolidStateIvpSolverOptions::default()
        };
        ivp.set_solver_options(options.clone());
        assert_eq!(ivp.solver_options, options);
    }

    #[test]
    fn changing_problem_after_solve_invalidates_published_trajectory() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        ivp.set_problem(2.0, 1.0, 300.0, 35_000.0, 2.0e5).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();
        ivp.solve().unwrap();
        assert!(ivp.get_result().is_ok());

        ivp.set_time(3.0).unwrap();
        assert!(ivp.get_result().is_err());
    }

    #[test]
    fn changing_thermal_and_arrhenius_inputs_after_model_selection_rebuilds_rhs() {
        let mut reconfigured = KineticModelIVP::new(SolverType::Bdf);
        reconfigured
            .set_problem(2.0, 1.0, 300.0, 35_000.0, 2.0e5)
            .unwrap();
        reconfigured
            .set_model(KineticModelNames::F1, vec![])
            .unwrap();
        reconfigured.set_time(3.0).unwrap();
        reconfigured.set_heating_program(315.0, 4.0).unwrap();
        reconfigured.set_arrhenius(45_000.0, 4.0e5).unwrap();
        reconfigured.solve().unwrap();

        let mut freshly_configured = KineticModelIVP::new(SolverType::Bdf);
        freshly_configured
            .set_problem(3.0, 4.0, 315.0, 45_000.0, 4.0e5)
            .unwrap();
        freshly_configured
            .set_model(KineticModelNames::F1, vec![])
            .unwrap();
        freshly_configured.solve().unwrap();

        let (changed_time, changed_state) = reconfigured.get_result().unwrap();
        let (fresh_time, fresh_state) = freshly_configured.get_result().unwrap();
        let changed_solution = reconfigured.get_solution().unwrap();
        let fresh_solution = freshly_configured.get_solution().unwrap();
        assert_eq!(changed_time.len(), fresh_time.len());
        assert_eq!(changed_state.shape(), fresh_state.shape());
        assert_eq!(
            changed_solution.conversion_rate.len(),
            fresh_solution.conversion_rate.len()
        );
        for (changed, fresh) in changed_time.iter().zip(fresh_time.iter()) {
            assert!((changed - fresh).abs() <= 1e-12);
        }
        for (changed, fresh) in changed_state.iter().zip(fresh_state.iter()) {
            assert!((changed - fresh).abs() <= 1e-12);
        }
        for (changed, fresh) in changed_solution
            .conversion_rate
            .iter()
            .zip(fresh_solution.conversion_rate.iter())
        {
            assert!((changed - fresh).abs() <= 1e-12);
        }
    }

    #[test]
    fn test_solve_returns_error_without_model() {
        let mut ivp = KineticModelIVP::new_lsode2();
        ivp.set_problem(10.0, 5.0, 298.15, 50000.0, 1e6).unwrap();

        let result = ivp.solve();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert_eq!(err, "kinetic model not set");
    }

    #[test]
    fn test_get_result_returns_error_without_solve() {
        let mut ivp = KineticModelIVP::new_lsode2();
        ivp.set_problem(10.0, 5.0, 298.15, 50000.0, 1e6).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();

        let result = ivp.get_result();
        assert!(result.is_err());
        match result {
            Err(err) => assert_eq!(err, "Solver not initialized. Call solve() first."),
            Ok(_) => panic!("get_result should require solve() first"),
        }
    }

    #[test]
    fn test_get_solution_returns_error_without_solve() {
        let mut ivp = KineticModelIVP::new_lsode2();
        ivp.set_problem(10.0, 5.0, 298.15, 50000.0, 1e6).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();

        let result = ivp.get_solution();
        assert!(result.is_err());
        match result {
            Err(err) => assert_eq!(err, "Solver not initialized. Call solve() first."),
            Ok(_) => panic!("get_solution should require solve() first"),
        }
    }

    #[test]
    fn test_plot_in_terminal_does_not_panic_before_solve() {
        let mut ivp = KineticModelIVP::new_lsode2();
        ivp.set_problem(10.0, 5.0, 298.15, 50000.0, 1e6).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();

        // Terminal plotting should become a no-op instead of panicking when no result exists.
        ivp.plot_in_terminal();
    }

    #[test]
    fn test_lsode2_handles_stiffer_d3_case() {
        let mut ivp = KineticModelIVP::new_lsode2();
        ivp.set_problem(10.0, 10.0, 298.15, 50000.0, 1e6).unwrap();
        ivp.set_model(KineticModelNames::D3, vec![]).unwrap();
        ivp.check_task().unwrap();
        ivp.solve().unwrap();

        let sol = ivp.get_solution().unwrap();
        assert!(sol.time.len() > 2);
        assert_eq!(sol.temperature.len(), sol.time.len());
        assert_eq!(sol.conversion.len(), sol.time.len());
        assert_eq!(sol.conversion_rate.len(), sol.time.len());
        assert!(sol.conversion.last().unwrap() >= sol.conversion.first().unwrap());
    }

    #[test]
    fn test_get_solution_conversion_rate_matches_derivative() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        ivp.set_problem(20.0, 5.0, 298.15, 40000.0, 5e5).unwrap();
        ivp.set_model(KineticModelNames::F1, vec![]).unwrap();
        ivp.check_task().unwrap();
        ivp.solve().unwrap();

        let sol = ivp.get_solution().unwrap();
        let n = sol.time.len();
        assert!(n > 2);

        let mut rel_err_sum = 0.0;
        let mut count = 0usize;
        for i in 1..(n - 1) {
            let dt = sol.time[i + 1] - sol.time[i - 1];
            if dt <= 0.0 {
                continue;
            }
            let da = sol.conversion[i + 1] - sol.conversion[i - 1];
            let deriv = da / dt;
            let rate = sol.conversion_rate[i];

            let denom = rate.abs().max(1e-12);
            let rel_err = (deriv - rate).abs() / denom;
            rel_err_sum += rel_err;
            count += 1;
        }

        let mean_rel_err = rel_err_sum / (count as f64);
        assert!(
            mean_rel_err < 0.2,
            "mean relative error too large: {}",
            mean_rel_err
        );
    }

    #[test]
    fn test_check_task_success() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);

        // Follow the correct order
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
        let _ = ivp.set_model(KineticModelNames::A2, vec![]);

        // Should not panic
        let _ = ivp.check_task();
    }

    #[test]

    fn test_check_task_fail_no_problem() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        let _ = ivp.set_model(KineticModelNames::A2, vec![]);

        // Should panic because problem not set
        let r = ivp.check_task();
        match r {
            Ok(_) => {}
            Err(e) => {
                assert_eq!(e, "t_final not set or invalid")
            }
        }
    }

    #[test]

    fn test_check_task_fail_no_model() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        // Should panic because model not set
        let r = ivp.check_task();
        match r {
            Ok(_) => {}
            Err(e) => {
                assert_eq!(e, "kinetic model not set")
            }
        }
    }

    #[test]
    fn test_complete_workflow() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);

        // Step 1: new() - already done
        assert!(ivp.kinmodel.is_none());

        // Step 2: set_problem()
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
        assert_eq!(ivp.t_final, 100.0);

        // Step 3: set_model()
        let _ = ivp.set_model(KineticModelNames::JMA, vec![2.0]);
        assert!(ivp.kinmodel.is_some());

        // Step 4: check_task()
        let _ = ivp.check_task(); // Should not panic
    }

    #[test]
    fn test_arrhenius_expression() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        let arrhenius = ivp.arrhenius();

        // Check that arrhenius expression contains expected components
        let expr_str = format!("{:?}", arrhenius);
        assert!(expr_str.contains("Exp"));
        assert!(expr_str.contains("Const(1000000.0)")); // A value
        assert!(expr_str.contains("Const(50000.0)")); // E value
    }

    #[test]
    fn test_different_solver_types() {
        let solver_types = [SolverType::Bdf, SolverType::BackwardEuler, SolverType::Bdf];

        for solver_type in solver_types {
            let mut ivp = KineticModelIVP::new(solver_type.clone());
            assert!(matches!(ivp.solvertype, ref _solver_type));

            let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
            let _ = ivp.set_model(KineticModelNames::F1, vec![]);
            let _ = ivp.check_task(); // Should work for all solver types
        }
    }

    #[test]
    fn test_multiple_models() {
        let test_cases = [
            (KineticModelNames::A2, vec![]),
            (KineticModelNames::F1, vec![]),
            (KineticModelNames::JMA, vec![2.0]),
            (KineticModelNames::SB, vec![1.0, 2.0, 0.5]),
            (KineticModelNames::PTe, vec![2.0, 1.0]),
        ];

        for (model_name, params) in test_cases {
            let mut ivp = KineticModelIVP::new(SolverType::Bdf);
            let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
            let _ = ivp.set_model(model_name.clone(), params);

            // Verify model was set correctly
            assert!(ivp.kinmodel.is_some());
            assert!(ivp.kinetic_expression().is_ok());

            // Should pass check_task
            let _ = ivp.check_task();
        }
    }

    #[test]
    fn test_parameter_validation_in_set_model() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6).unwrap();

        // This should return error because wrong number of parameters
        let result = ivp.set_model(KineticModelNames::JMA, vec![1.0, 2.0]); // JMA needs 1 param, not 2
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("JMA model requires 1 parameter")
        );
    }

    #[test]
    fn test_kinetic_expression_contains_arrhenius() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
        let _ = ivp.set_model(KineticModelNames::F1, vec![]);

        // The kinetic expression should be arrhenius * kinetic_model
        let expr_str = format!("{:?}", ivp.kinetic_expression().unwrap());

        // Should contain multiplication (arrhenius * model)
        assert!(expr_str.contains("Mul") || expr_str.contains("*"));
        // Should contain exponential (from arrhenius)
        assert!(expr_str.contains("Exp"));
        // Should contain variable 'a' (from kinetic model)
        assert!(expr_str.contains("Var(\"a\")"));
    }
    #[test]
    fn solve_with_radau() {
        let t = Instant::now();
        let mut ivp = KineticModelIVP::new(SolverType::Radau);
        let _ = ivp.set_problem(10.0, 10.0, 298.15, 50000.0, 1e6);

        // Test parameterized model
        let _ = ivp.set_model(KineticModelNames::D3, vec![]);
        let _ = ivp.check_task();
        let _ = ivp.solve();
        println!("Time elapsed: {:?}", t.elapsed().as_secs());
        let _ = ivp.plot();
    }

    #[test]
    fn solve_with_be() {
        let mut ivp = KineticModelIVP::new(SolverType::BackwardEuler);
        let _ = ivp.set_problem(10.0, 10.0, 298.15, 50000.0, 1e6);

        // Test parameterized model
        let _ = ivp.set_model(KineticModelNames::D3, vec![]);
        let _ = ivp.check_task();
        let _ = ivp.solve();
        let _ = ivp.plot();
    }
}

#[cfg(test)]
mod error_handling_tests {
    use std::time::Instant;

    use super::*;
    use RustedSciThe::numerical::StiffAPI::StiffSolverKind as SolverType;

    #[test]
    fn test_set_problem_validation() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);

        // Test negative t_final
        assert!(ivp.set_problem(-1.0, 10.0, 298.15, 50000.0, 1e6).is_err());

        // Test negative beta
        assert!(ivp.set_problem(100.0, -10.0, 298.15, 50000.0, 1e6).is_err());

        // Test negative T0
        assert!(ivp.set_problem(100.0, 10.0, -298.15, 50000.0, 1e6).is_err());

        // Test negative E
        assert!(ivp.set_problem(100.0, 10.0, 298.15, -50000.0, 1e6).is_err());

        // Test negative A
        assert!(ivp.set_problem(100.0, 10.0, 298.15, 50000.0, -1e6).is_err());

        // Test valid parameters
        assert!(ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6).is_ok());
    }

    #[test]
    fn test_set_model_error_propagation() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);
        ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6).unwrap();

        // Test invalid parameters
        let result = ivp.set_model(KineticModelNames::JMA, vec![1.0, 2.0]);
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("JMA model requires 1 parameter")
        );

        // Test valid parameters
        let result = ivp.set_model(KineticModelNames::JMA, vec![2.0]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_check_task_detailed_errors() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);

        // Test uninitialized state
        let result = ivp.check_task();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("t_final not set"));

        // Set problem partially
        ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6).unwrap();
        let result = ivp.check_task();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("kinetic model not set"));

        // Set model
        ivp.set_model(KineticModelNames::A2, vec![]).unwrap();
        let result = ivp.check_task();
        assert!(result.is_ok());
    }

    #[test]
    fn test_solve_without_setup() {
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);

        // Try to solve without setup
        let result = ivp.solve();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not set"));
    }

    #[test]
    fn test_gnuplot_without_solve() {
        let ivp = KineticModelIVP::new(SolverType::Bdf);

        let result = ivp.plot();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Solver not initialized"));
    }

    #[test]
    fn test_save_result_without_solve() {
        let ivp = KineticModelIVP::new(SolverType::Bdf);

        let result = ivp.save_result();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Solver not initialized"));
    }

    #[test]
    fn test_complete_error_free_workflow() {
        let t = Instant::now();
        let mut ivp = KineticModelIVP::new(SolverType::Radau);

        // All operations should return Ok
        assert!(ivp.set_problem(50.0, 10.0, 298.15, 50000.0, 1e6).is_ok());
        ivp.set_stop_condition(Some(0.99));
        assert!(ivp.set_model(KineticModelNames::F1, vec![]).is_ok());
        assert!(ivp.check_task().is_ok());
        assert!(ivp.solve().is_ok());
        assert!(ivp.plot().is_ok());
        assert!(ivp.save_result().is_ok());
        println!("Time taken: {:?}", t.elapsed().as_secs());
        ivp.plot_in_terminal()
    }

    #[test]
    fn test_complete_error_free_workflow2() {
        let t = Instant::now();
        let mut ivp = KineticModelIVP::new(SolverType::Bdf);

        // All operations should return Ok
        assert!(ivp.set_problem(50.0, 2.0, 298.15, 50000.0, 1e6).is_ok());
        ivp.set_stop_condition(Some(0.99));
        assert!(ivp.set_model(KineticModelNames::F1, vec![]).is_ok());
        assert!(ivp.check_task().is_ok());
        assert!(ivp.solve().is_ok());
        assert!(ivp.plot().is_ok());
        //  assert!(ivp.save_result().is_ok());
        println!("Time taken: {:?}", t.elapsed().as_secs());
        ivp.plot_in_terminal()
    }

    fn independent_rate(name: KineticModelNames, params: &[f64], a: f64) -> f64 {
        let remaining = 1.0 - a;
        let logarithmic_conversion = -remaining.ln();
        match name {
            KineticModelNames::A2 => 2.0 * remaining * logarithmic_conversion.powf(0.5),
            KineticModelNames::A3 => 3.0 * remaining * logarithmic_conversion.powf(2.0 / 3.0),
            KineticModelNames::A4 => 4.0 * remaining * logarithmic_conversion.powf(0.75),
            KineticModelNames::R2 => remaining.powf(0.5),
            KineticModelNames::R3 => remaining.powf(2.0 / 3.0),
            KineticModelNames::D1 => a.powf(-1.0),
            KineticModelNames::D2 => logarithmic_conversion.powf(-1.0),
            KineticModelNames::D3 => {
                1.5 * remaining.powf(2.0 / 3.0) / (1.0 - remaining.powf(1.0 / 3.0))
            }
            KineticModelNames::D4 => 1.5 / (remaining.powf(-1.0 / 3.0) - 1.0),
            KineticModelNames::P2_3 => a.powf(-0.5),
            KineticModelNames::P2 => a.powf(0.5),
            KineticModelNames::P3 => a.powf(0.75),
            KineticModelNames::F1 => remaining,
            KineticModelNames::F2 => remaining.powi(2),
            KineticModelNames::F3 => remaining.powi(3),
            KineticModelNames::SB => {
                remaining.powf(params[0])
                    * a.powf(params[1])
                    * logarithmic_conversion.powf(params[2])
            }
            KineticModelNames::JMA => {
                params[0] * remaining * logarithmic_conversion.powf(1.0 - 1.0 / params[0])
            }
            KineticModelNames::Ac => params[0] * a.powf(1.0 - 1.0 / params[0]),
            KineticModelNames::Dec => remaining.powf(params[0]),
            KineticModelNames::PTe => remaining.powf(params[0]) * a.powf(params[1]),
            KineticModelNames::SBtp => params[2] * remaining.powf(params[0]) * a.powf(params[1]),
        }
    }

    #[test]
    fn kinetic_rate_formulas_and_jacobians_match_independent_references() {
        let a = 0.37;
        let derivative_step = 1e-6;
        let parameters = |name| match name {
            KineticModelNames::SB => vec![1.7, 1.2, 0.6],
            KineticModelNames::JMA | KineticModelNames::Ac | KineticModelNames::Dec => vec![1.7],
            KineticModelNames::PTe => vec![1.7, 1.2],
            KineticModelNames::SBtp => vec![1.7, 1.2, 0.8],
            _ => vec![],
        };

        for name in KineticModelNames::iter() {
            let params = parameters(name);
            let (model, jacobian) = create_kinetic_model(name, params.clone()).unwrap();
            let rate = model.get_fn();
            let jacobian = jacobian.lambdify1D();
            let expected = independent_rate(name, &params, a);
            let actual = rate(a);
            assert!(
                (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
                "{name:?} rate mismatch at a={a}: actual={actual:.16e}, expected={expected:.16e}"
            );

            let expected_derivative = (independent_rate(name, &params, a + derivative_step)
                - independent_rate(name, &params, a - derivative_step))
                / (2.0 * derivative_step);
            let actual_derivative = jacobian(a);
            assert!(
                (actual_derivative - expected_derivative).abs()
                    <= 2e-5 * expected_derivative.abs().max(1.0),
                "{name:?} Jacobian mismatch at a={a}: actual={actual_derivative:.16e}, expected={expected_derivative:.16e}"
            );

            for edge_a in [1e-6, 1.0 - 1e-6] {
                let edge_actual = rate(edge_a);
                let edge_expected = independent_rate(name, &params, edge_a);
                assert!(
                    edge_actual.is_finite()
                        && (edge_actual - edge_expected).abs()
                            <= 1e-10 * edge_expected.abs().max(1.0),
                    "{name:?} edge value mismatch at a={edge_a}: actual={edge_actual:.16e}, expected={edge_expected:.16e}"
                );
            }
        }
    }

    #[test]
    fn analytical_integrals_are_consistent_with_rate_models() {
        let a = 0.37;
        let h = 1e-6;
        for name in KineticModelNames::iter() {
            let Ok(integral) = name.get_g() else {
                continue;
            };
            let model = name.from_names_to_expr(vec![]);
            let rate = model.get_fn()(a);
            let derivative = (integral(a + h) - integral(a - h)) / (2.0 * h);
            assert!(
                (derivative * rate - 1.0).abs() <= 2e-5,
                "{name:?}: g'(a)*f(a) should equal 1, got {:.16e} at a={a}",
                derivative * rate
            );
        }
    }
}
