//! Gated AOT story tests for the condensed-reactor IVP path.
//!
//! These tests verify the real solver backend contract on a small condensed
//! combustion story fixture. The hypotheses are:
//! - Lambdify remains the ordinary baseline route;
//! - AOT with `tcc` works on the same canonical condensed problem;
//! - sparse and banded matrix routes both survive the AOT handoff;
//! - AOT results stay numerically close to the Lambdify baseline.
//!
//! Expected outcome:
//! - tests are skipped unless `KITHE_RUN_IVP_AOT_TESTS=1` is set;
//! - `tcc` must be available in `PATH`;
//! - the solve finishes with accepted status;
//! - baseline and AOT solutions remain finite and close.

#[cfg(test)]
mod tests {
    use super::super::SimpleReactorIVP::*;
    use super::super::solver_backend::{
        ReactorIvpExecutionBackend, ReactorIvpMatrixBackend, ReactorIvpMethod,
        ReactorIvpSolverConfig, ReactorIvpSymbolicBackend,
    };
    use crate::Kinetics::mechfinder_api::ReactionData;
    use RustedSciThe::numerical::LSODE2::{
        Lsode2AotProfile, Lsode2AotToolchain, Lsode2BackendConfig, Lsode2ProblemConfig,
        solver::Lsode2Solver,
    };
    use approx::assert_abs_diff_eq;
    use std::collections::HashMap;
    use std::process::Command;

    fn command_available(command: &str) -> bool {
        let probe = if cfg!(windows) { "where" } else { "which" };
        Command::new(probe)
            .arg(command)
            .output()
            .map(|output| output.status.success())
            .unwrap_or(false)
    }

    fn should_run_aot_tests() -> bool {
        matches!(
            std::env::var("KITHE_RUN_IVP_AOT_TESTS").ok().as_deref(),
            Some("1") | Some("true") | Some("TRUE") | Some("yes") | Some("YES")
        ) && command_available("tcc")
    }

    fn build_canonical_aot_task(matrix_backend: ReactorIvpMatrixBackend) -> SimpleReactorTask {
        let mut task = SimpleReactorTask::new();
        task.kindata
            .set_reaction_data_directly(
                vec![ReactionData::new_elementary(
                    "A=>B".to_string(),
                    vec![1.0e3, 0.0, 1.0e4],
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
                .with_symbolic_backend(ReactorIvpSymbolicBackend::AtomView)
                .with_execution_backend(ReactorIvpExecutionBackend::Lambdify)
                .with_matrix_backend(matrix_backend)
                .with_integration_domain(0.0, 0.01)
                .with_max_step(1.0e-3)
                .with_rtol(1.0e-6)
                .with_atol(1.0e-8),
        );
        task
    }

    fn build_problem_config(
        matrix_backend: ReactorIvpMatrixBackend,
        execution_backend: ReactorIvpExecutionBackend,
    ) -> Lsode2ProblemConfig {
        let mut task = build_canonical_aot_task(matrix_backend);
        task.solver_backend_config.execution_backend = execution_backend;
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");
        task.build_lsode2_problem_config()
            .expect("problem config should build")
    }

    fn solve_problem_config(config: Lsode2ProblemConfig) -> (Vec<f64>, Vec<f64>, String) {
        let mut solver = Lsode2Solver::new(config).expect("solver should build");
        let summary = solver
            .solve_with_summary()
            .expect("LSODE2 solve should finish");
        let (axis, values) = solver.get_result();
        let axis = axis.as_slice().to_vec();
        let final_row = values
            .row(values.nrows() - 1)
            .iter()
            .copied()
            .collect::<Vec<_>>();
        (axis, final_row, summary.status)
    }

    fn assert_close_result(
        baseline: &(Vec<f64>, Vec<f64>, String),
        aot: &(Vec<f64>, Vec<f64>, String),
    ) {
        assert!(baseline.2.starts_with("finished"));
        assert!(aot.2.starts_with("finished"));
        assert_eq!(baseline.0.len(), aot.0.len());
        assert_eq!(baseline.1.len(), aot.1.len());
        for (lhs, rhs) in baseline.0.iter().zip(aot.0.iter()) {
            assert_abs_diff_eq!(lhs, rhs, epsilon = 1e-10);
        }
        for (lhs, rhs) in baseline.1.iter().zip(aot.1.iter()) {
            let scale = lhs.abs().max(rhs.abs()).max(1.0);
            let rel = (lhs - rhs).abs() / scale;
            assert!(
                rel <= 1.0e-8,
                "AOT result drifted too far from Lambdify baseline: lhs={lhs:e}, rhs={rhs:e}, rel={rel:e}"
            );
        }
    }

    #[test]
    fn aot_sparse_tcc_matches_lambdify_baseline() {
        if !should_run_aot_tests() {
            eprintln!(
                "Skipping IVP AOT sparse story test: set KITHE_RUN_IVP_AOT_TESTS=1 and ensure tcc is on PATH."
            );
            return;
        }

        let baseline = solve_problem_config(build_problem_config(
            ReactorIvpMatrixBackend::Sparse,
            ReactorIvpExecutionBackend::Lambdify,
        ));

        let output_dir = tempfile::tempdir().expect("temporary AOT output dir should be created");
        let aot_config = build_problem_config(
            ReactorIvpMatrixBackend::Sparse,
            ReactorIvpExecutionBackend::Aot {
                toolchain: Lsode2AotToolchain::CTcc,
                profile: Lsode2AotProfile::Release,
            },
        )
        .with_backend(Lsode2BackendConfig::native_sparse_faer_aot_c_tcc(
            output_dir.path(),
        ))
        .with_aot_parallel_chunking(2);

        let aot = solve_problem_config(aot_config);
        assert_close_result(&baseline, &aot);
    }

    #[test]
    fn aot_banded_tcc_matches_lambdify_baseline() {
        if !should_run_aot_tests() {
            eprintln!(
                "Skipping IVP AOT banded story test: set KITHE_RUN_IVP_AOT_TESTS=1 and ensure tcc is on PATH."
            );
            return;
        }

        let matrix_backend = ReactorIvpMatrixBackend::Banded;
        let baseline = solve_problem_config(build_problem_config(
            matrix_backend,
            ReactorIvpExecutionBackend::Lambdify,
        ));

        let output_dir = tempfile::tempdir().expect("temporary AOT output dir should be created");
        let aot_config = build_problem_config(
            matrix_backend,
            ReactorIvpExecutionBackend::Aot {
                toolchain: Lsode2AotToolchain::CTcc,
                profile: Lsode2AotProfile::Release,
            },
        )
        .with_backend(Lsode2BackendConfig::native_banded_faithful_aot_c_tcc(
            output_dir.path(),
        ))
        .with_aot_parallel_chunking(2);

        let aot = solve_problem_config(aot_config);
        assert_close_result(&baseline, &aot);
    }
}
