//! Quality-regression tests for the condensed reactor IVP path.
//!
//! Hypotheses:
//! - the solver must honor an explicit stop condition on a selected species
//!   when the supported native faithful route is enabled;
//! - postprocessed species concentrations should stay non-negative for a stable
//!   condensed-phase combustion fixture;
//! - the conservation report should stay within a few percent on the canonical
//!   one-reaction regression problem.
//!
//! Expected result:
//! the IVP reactor behaves like a numerically sane combustion model rather than
//! a loose copied BVP prototype.

#[cfg(test)]
mod tests {
    use super::super::SimpleReactorIVP::*;
    use super::super::solver_backend::{ReactorIvpMethod, ReactorIvpSolverConfig};
    use crate::Kinetics::mechfinder_api::ReactionData;
    use RustedSciThe::numerical::LSODE2::Lsode2NativeExecutionConfig;
    use approx::assert_abs_diff_eq;
    use std::collections::HashMap;

    fn build_balanced_reaction_task(
        initial_h2: f64,
        initial_o2: f64,
        initial_h2o2: f64,
    ) -> SimpleReactorTask {
        let mut task = SimpleReactorTask::new();
        task.kindata
            .set_reaction_data_directly(
                vec![ReactionData::new_elementary(
                    "H2+O2=>H2O2".to_string(),
                    // Keep the fixture chemically balanced but dynamically inert
                    // so the conservation contract is measured against a stable
                    // reference state rather than a fast-moving transient.
                    vec![0.0, 0.0, 8.0e3],
                    None,
                )],
                None,
            )
            .expect("reaction setup should succeed");
        task.kindata.stecheodata.vec_of_molmasses = Some(vec![2.016, 31.998, 34.014]);
        task.kindata.substances = vec!["H2".to_string(), "O2".to_string(), "H2O2".to_string()];
        task.thermal_effects = vec![0.0];
        task.ro = 1200.0;
        task.Cp = 1000.0;
        task.Lambda = 0.25;
        task.m = 0.015;
        task.scaling =
            crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig::new(100.0, 0.05, 100.0);
        task.initial_conditions = HashMap::from([
            ("T".to_string(), 450.0),
            ("q".to_string(), 0.0),
            ("H2".to_string(), initial_h2),
            ("O2".to_string(), initial_o2),
            ("H2O2".to_string(), initial_h2o2),
        ]);
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default()
                .with_method(ReactorIvpMethod::Auto)
                .with_native_execution(Lsode2NativeExecutionConfig::faithful_bdf_solve(
                    200_000, 200_000,
                ))
                .with_integration_domain(0.0, 0.001)
                .with_max_step(1.0e-3)
                .with_rtol(1.0e-5)
                .with_atol(1.0e-8),
        );
        task
    }

    fn build_baseline_task() -> SimpleReactorTask {
        build_balanced_reaction_task(0.7, 0.3, 0.0)
    }

    fn build_stop_condition_task() -> SimpleReactorTask {
        let mut task = build_balanced_reaction_task(0.1, 0.9, 0.0);
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default()
                .with_method(ReactorIvpMethod::Auto)
                .with_native_execution(Lsode2NativeExecutionConfig::faithful_bdf_solve(
                    200_000, 200_000,
                ))
                .with_stop_condition_le(0, 0.2)
                .with_integration_domain(0.0, 0.001)
                .with_max_step(1.0e-3)
                .with_rtol(1.0e-5)
                .with_atol(1.0e-8),
        );
        task
    }

    #[test]
    fn test_stop_condition_terminates_the_run_when_initial_state_is_below_threshold() {
        let mut task = build_stop_condition_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        let snapshot = task.solve().expect("condensed solve should succeed");

        assert!(
            snapshot
                .summary
                .native_termination_kind
                .as_deref()
                .is_some_and(|kind| kind == "reached_stop_condition")
                || snapshot.status().contains("stopped_by_condition"),
            "unexpected stop status: {}, native termination: {:?}",
            snapshot.status(),
            snapshot.summary.native_termination_kind
        );

        let values = task
            .solution_render_data()
            .expect("render data should exist after solve")
            .solution;
        let a_column = values.column(2).into_owned();
        assert!(
            a_column.iter().all(|value| *value >= -1.0e-10),
            "A concentration should not go negative: {:?}",
            a_column
        );
        assert!(
            a_column[a_column.len() - 1] <= 0.2 + 1.0e-8,
            "stop condition should cap A at the selected threshold, got {}",
            a_column[a_column.len() - 1]
        );
    }

    #[test]
    fn test_condensed_conservation_stays_within_a_few_percent() {
        let mut task = build_baseline_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        task.solve().expect("condensed solve should succeed");

        let report = task
            .latest_conservation_report()
            .expect("conservation report should be available after solve");

        assert!(
            report.energy_balance_error_rel.abs() <= 5.0,
            "energy balance relative error too large: {}",
            report.energy_balance_error_rel
        );
        assert!(
            report
                .atomic_mass_balance_error
                .iter()
                .all(|(_, value)| value.abs() <= 5.0e-2),
            "atomic balance errors should stay within a few percent: {:?}",
            report.atomic_mass_balance_error
        );
        assert!(
            report
                .sum_of_mass_fractions
                .iter()
                .all(|(_, value)| (*value - 1.0).abs() <= 5.0e-2),
            "mass-fraction sums should stay close to unity: {:?}",
            report.sum_of_mass_fractions
        );

        let rendered = task
            .solution_render_data()
            .expect("render data should exist after solve");
        assert_eq!(
            rendered.unknowns.first().map(String::as_str),
            Some("T"),
            "temperature should be exported with a dimensional label"
        );
        assert!(
            !rendered.unknowns.iter().any(|name| name == "Teta"),
            "legacy dimensionless temperature label should not leak into rendered output: {:?}",
            rendered.unknowns
        );
        for (index, col_index) in (2..rendered.solution.ncols()).enumerate() {
            let column = rendered.solution.column(col_index).into_owned();
            assert!(
                column.iter().all(|value| *value >= -1.0e-10),
                "species column {} went negative: {:?}",
                index,
                column
            );
        }
        assert_abs_diff_eq!(rendered.x_mesh[0], 0.0, epsilon = 1e-12);
    }
}
