//! Symbolic frontend parity and prepared-request boundary stories.

#[cfg(test)]
mod tests {
    use super::super::reactor_ivp_test_support::multistep_oxidation_coverage_task;
    use super::super::solver_backend::{
        ReactorIvpMatrixBackend, ReactorIvpSolverFamily, ReactorIvpSymbolicBackend,
    };
    use crate::Kinetics::mechfinder_api::ReactionData;
    use RustedSciThe::numerical::SolverTelemetryMode;
    use RustedSciThe::numerical::StiffAPI::StiffSolverRequest;
    use std::collections::HashMap;

    #[test]
    fn multistep_oxidation_expr_and_atomview_publish_matching_bounded_profiles() {
        let mut reference: Option<Vec<f64>> = None;
        for frontend in [
            ReactorIvpSymbolicBackend::ExprLegacy,
            ReactorIvpSymbolicBackend::AtomView,
        ] {
            let mut task = multistep_oxidation_coverage_task();
            task.kindata
                .set_reaction_data_directly(
                    vec![
                        ReactionData::new_elementary(
                            "2CH4+3O2=>2CO+4H2O".to_string(),
                            vec![2.5e2, 0.0, 8.0e3],
                            None,
                        ),
                        ReactionData::new_elementary(
                            "2CO+O2=>2CO2".to_string(),
                            vec![1.5e2, 0.0, 5.0e3],
                            None,
                        ),
                    ],
                    None,
                )
                .expect("two-step oxidation network should parse");
            task.initial_conditions = HashMap::from([
                ("T".into(), 520.0),
                ("q".into(), 0.0),
                ("CH4".into(), 0.2),
                ("O2".into(), 0.4),
                ("CO".into(), 0.1),
                ("H2O".into(), 0.1),
                ("CO2".into(), 0.2),
            ]);
            let config = task
                .solver_backend_config
                .with_family(ReactorIvpSolverFamily::Lsode2)
                .with_symbolic_backend(frontend)
                .with_matrix_backend(ReactorIvpMatrixBackend::Dense)
                .with_integration_domain(0.0, 1.0e-4)
                .with_max_step(1.0e-5)
                .with_rtol(1.0e-6)
                .with_atol(1.0e-9)
                .with_telemetry(SolverTelemetryMode::Detailed);
            task.set_solver_backend_config(config);
            let snapshot = task
                .solve()
                .unwrap_or_else(|error| panic!("{frontend:?} multistep solve failed: {error}"));
            assert!(snapshot.status().starts_with("finished"));
            let temperature = snapshot
                .values
                .column(0)
                .map(|value| value * task.scaling.T_scale + task.scaling.dT);
            assert!(
                temperature
                    .iter()
                    .all(|value| { value.is_finite() && (0.0..=3000.0).contains(value) })
            );
            for (column, name) in snapshot.variable_names.iter().enumerate() {
                if name == "Teta" || name == "q" {
                    continue;
                }
                assert!(
                    snapshot
                        .values
                        .column(column)
                        .iter()
                        .all(|fraction| { fraction.is_finite() && (0.0..=1.0).contains(fraction) })
                );
            }
            let profile = snapshot.values.iter().copied().collect::<Vec<_>>();
            if let Some(expected) = &reference {
                assert_eq!(profile.len(), expected.len());
                let max_delta = profile
                    .iter()
                    .zip(expected)
                    .map(|(actual, expected)| (actual - expected).abs())
                    .fold(0.0_f64, f64::max);
                println!("multistep ExprLegacy/AtomView max_abs_profile_delta={max_delta:.6e}");
                assert!(max_delta <= 1.0e-7, "frontend profile drift {max_delta:e}");
            } else {
                reference = Some(profile);
            }
            assert!(!snapshot.stage_telemetry.stages.is_empty());
        }
    }

    #[test]
    fn fixture_builder_constructs_distinct_route_request_type() {
        fn accepts_stiff_request(_: StiffSolverRequest) {}
        let mut task = multistep_oxidation_coverage_task();
        task.setup_condensed_ivp().unwrap();
        accepts_stiff_request(task.build_condensed_stiff_request().unwrap());
    }
}
