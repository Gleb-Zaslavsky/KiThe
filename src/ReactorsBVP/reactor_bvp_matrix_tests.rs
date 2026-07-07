#[cfg(test)]
mod tests {
    use super::super::SimpleReactorBVP::NrbvpHandoffConfig;
    use super::super::reactor_bvp_test_support::compact_hmx_reactor;
    use super::super::solver_backend::{
        ReactorBvpAotCompiler, ReactorBvpMatrixBackend, ReactorBvpSolverConfig,
        ReactorBvpSymbolicBackend,
    };
    use crate::ReactorsBVP::SimpleReactorBVP::ReactorError;
    use crate::ReactorsBVP::reactor_BVP_utils::{BoundsConfig, ToleranceConfig};
    use nalgebra::DMatrix;
    use std::process::Command;

    /// Build a compact initial guess for the tiny HMX combustion fixture.
    fn compact_initial_guess(unknowns: usize, n_steps: usize) -> DMatrix<f64> {
        DMatrix::from_element(unknowns, n_steps, 0.5)
    }

    /// Build the legacy NRBVP handoff wrapper with a selected solver facade.
    fn compact_hmx_handoff(
        initial_guess: DMatrix<f64>,
        solver_backend_config: ReactorBvpSolverConfig,
    ) -> NrbvpHandoffConfig {
        NrbvpHandoffConfig::new(
            initial_guess,
            0.0,
            1.0,
            32,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Sparse".to_string(),
            1e-8,
            Some(
                ToleranceConfig::new(1e-8, 1e-8, 1e-8, 1e-8)
                    .to_full_tolerance_map(&["HMX".to_string(), "HMXprod".to_string()]),
            ),
            80,
            Some(
                BoundsConfig::new((-0.1, 1.1), (-1e20, 1e20), (-100.0, 100.0), (-1e20, 1e20))
                    .to_full_bounds_map(&["HMX".to_string(), "HMXprod".to_string()]),
            ),
            Some("warn".to_string()),
            false,
        )
        .with_solver_backend_config(solver_backend_config)
    }

    /// Solve the compact combustion fixture through the NRBVP backend snapshot.
    fn solve_compact_backend(
        solver_backend_config: ReactorBvpSolverConfig,
    ) -> Result<(Vec<f64>, Vec<f64>), ReactorError> {
        let mut reactor = compact_hmx_reactor();
        reactor.setup_bvp()?;

        let initial_guess = compact_initial_guess(reactor.solver.unknowns.len(), 32);
        let mut backend = reactor
            .solver
            .build_nrbvp_backend(compact_hmx_handoff(initial_guess, solver_backend_config))?;

        backend.before_solve_preprocessing();
        let solve_result = backend
            .try_solve()
            .map_err(|err| ReactorError::CalculationError(format!("{err:?}")))?;
        if solve_result.is_none() {
            return Err(ReactorError::CalculationError(
                "compact backend solve did not converge".to_string(),
            ));
        }

        let solution = backend.get_result().ok_or_else(|| {
            ReactorError::CalculationError(
                "compact backend solve did not return a solution".to_string(),
            )
        })?;
        if solution.nrows() == 0 || solution.ncols() == 0 {
            return Err(ReactorError::CalculationError(
                "compact backend returned an empty solution".to_string(),
            ));
        }
        if solution.iter().any(|value| !value.is_finite()) {
            return Err(ReactorError::InvalidNumericValue(
                "compact backend returned non-finite values".to_string(),
            ));
        }
        if backend.x_mesh.is_empty() {
            return Err(ReactorError::CalculationError(
                "compact backend returned an empty mesh".to_string(),
            ));
        }
        if backend.x_mesh.iter().any(|value| !value.is_finite()) {
            return Err(ReactorError::InvalidNumericValue(
                "compact backend returned non-finite mesh values".to_string(),
            ));
        }

        Ok((
            solution.iter().copied().collect(),
            backend.x_mesh.iter().copied().collect(),
        ))
    }

    /// Assert that two compact-task solution profiles are numerically close.
    fn assert_profile_close(reference: &[f64], candidate: &[f64], atol: f64, label: &str) {
        assert_eq!(
            reference.len(),
            candidate.len(),
            "{label} length mismatch: reference={}, candidate={}",
            reference.len(),
            candidate.len()
        );
        for (idx, (expected, actual)) in reference.iter().zip(candidate.iter()).enumerate() {
            assert!(
                (*expected - *actual).abs() <= atol,
                "{label} mismatch at index {idx}: expected {expected}, got {actual}, atol={atol}"
            );
            assert!(actual.is_finite(), "{label} contains non-finite value at {idx}");
        }
    }

    /// Check whether the gated AOT integration tests are allowed to run on this machine.
    fn aot_tests_enabled(toolchain: &str) -> bool {
        if std::env::var("KITHE_RUN_BVP_AOT_TESTS").is_err() {
            eprintln!("skipping AOT smoke test: set KITHE_RUN_BVP_AOT_TESTS=1 to enable it");
            return false;
        }
        if Command::new(toolchain).arg("--version").output().is_err() {
            eprintln!("skipping AOT smoke test: {toolchain} is not available on this machine");
            return false;
        }
        true
    }

    /// Compare two compact-task solver configurations and optionally gate them on a toolchain.
    fn assert_config_matches_reference(
        reference_config: ReactorBvpSolverConfig,
        candidate_config: ReactorBvpSolverConfig,
        toolchain: Option<&str>,
        atol: f64,
    ) {
        if let Some(toolchain) = toolchain {
            if !aot_tests_enabled(toolchain) {
                return;
            }
        }

        let (reference_solution, reference_mesh) = solve_compact_backend(reference_config)
            .expect("reference backend solve should succeed");
        let (candidate_solution, candidate_mesh) = solve_compact_backend(candidate_config)
            .expect("candidate backend solve should succeed");

        assert_profile_close(&reference_solution, &candidate_solution, atol, "solution");
        assert_profile_close(&reference_mesh, &candidate_mesh, atol, "mesh");
    }

    /// Default and sparse lambdify routes should differ only where the matrix backend changes.
    #[test]
    fn matrix_default_and_sparse_lambdify_configs_are_distinct() {
        assert_eq!(
            ReactorBvpMatrixBackend::default(),
            ReactorBvpMatrixBackend::Banded
        );

        let default_options = ReactorBvpSolverConfig::default_lambdify()
            .to_rusted_options()
            .expect("default config should convert");
        let sparse_options = ReactorBvpSolverConfig::sparse_lambdify()
            .to_rusted_options()
            .expect("sparse config should convert");

        assert_eq!(
            default_options
                .generated_backend_config
                .symbolic_assembly_backend,
            RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            sparse_options
                .generated_backend_config
                .symbolic_assembly_backend,
            RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            default_options
                .generated_backend_config
                .matrix_backend_override,
            Some(RustedSciThe::symbolic::codegen::codegen_provider_api::MatrixBackend::Banded)
        );
        assert_eq!(
            sparse_options
                .generated_backend_config
                .matrix_backend_override,
            None
        );
    }

    /// Compact combustion should solve through the default lambdify facade.
    #[test]
    fn matrix_compact_hmx_default_lambdify_solves() {
        let (solution, mesh) = solve_compact_backend(ReactorBvpSolverConfig::default_lambdify())
            .expect("default lambdify solve should succeed");

        assert!(!solution.is_empty());
        assert!(!mesh.is_empty());
    }

    /// Sparse lambdify should stay usable on the same compact combustion fixture.
    #[test]
    fn matrix_compact_hmx_sparse_lambdify_solves() {
        let (solution, mesh) = solve_compact_backend(ReactorBvpSolverConfig::sparse_lambdify())
            .expect("sparse lambdify solve should succeed");

        assert!(!solution.is_empty());
        assert!(!mesh.is_empty());
    }

    /// ExprLegacy + Sparse should stay compatible with the modern default profile on the same task.
    #[test]
    fn matrix_compact_hmx_legacy_expr_sparse_matches_default_profile() {
        let legacy_sparse = ReactorBvpSolverConfig::sparse_lambdify()
            .with_symbolic_backend(ReactorBvpSymbolicBackend::ExprLegacy);

        assert_config_matches_reference(
            ReactorBvpSolverConfig::default_lambdify(),
            legacy_sparse,
            None,
            1e-6,
        );
    }

    /// AOT smoke tests are gated because toolchains like tcc are optional on user machines.
    #[test]
    fn matrix_compact_hmx_aot_tcc_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::default_lambdify(),
            ReactorBvpSolverConfig::banded_aot(ReactorBvpAotCompiler::CTcc),
            Some("tcc"),
            1e-6,
        );
    }

    /// GCC-backed AOT smoke test for the same compact combustion fixture.
    #[test]
    fn matrix_compact_hmx_aot_gcc_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::default_lambdify(),
            ReactorBvpSolverConfig::banded_aot(ReactorBvpAotCompiler::CGcc),
            Some("gcc"),
            1e-6,
        );
    }

    /// Zig-backed AOT smoke test for the same compact combustion fixture.
    #[test]
    fn matrix_compact_hmx_aot_zig_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::default_lambdify(),
            ReactorBvpSolverConfig::banded_aot(ReactorBvpAotCompiler::Zig),
            Some("zig"),
            1e-6,
        );
    }

    /// Sparse AOT should remain available as a compatibility matrix route when the toolchain exists.
    #[test]
    fn matrix_compact_hmx_sparse_aot_tcc_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::sparse_lambdify(),
            ReactorBvpSolverConfig::sparse_aot(ReactorBvpAotCompiler::CTcc),
            Some("tcc"),
            1e-6,
        );
    }

    /// Sparse AOT with GCC should keep the same compact-task smoke coverage.
    #[test]
    fn matrix_compact_hmx_sparse_aot_gcc_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::sparse_lambdify(),
            ReactorBvpSolverConfig::sparse_aot(ReactorBvpAotCompiler::CGcc),
            Some("gcc"),
            1e-6,
        );
    }

    /// Sparse AOT with Zig should keep the same compact-task smoke coverage.
    #[test]
    fn matrix_compact_hmx_sparse_aot_zig_is_gated() {
        assert_config_matches_reference(
            ReactorBvpSolverConfig::sparse_lambdify(),
            ReactorBvpSolverConfig::sparse_aot(ReactorBvpAotCompiler::Zig),
            Some("zig"),
            1e-6,
        );
    }
}
