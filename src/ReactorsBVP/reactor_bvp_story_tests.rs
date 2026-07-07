#[cfg(test)]
mod tests {
    use super::super::SimpleReactorBVP::NrbvpHandoffConfig;
    use super::super::reactor_bvp_test_support::compact_hmx_reactor;
    use super::super::solver_backend::ReactorBvpSolverConfig;
    use crate::Kinetics::User_reactions::KinDataState;
    use RustedSciThe::symbolic::codegen::codegen_provider_api::MatrixBackend;
    use RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend;
    use nalgebra::DMatrix;

    /// End-to-end story test for the modern default reactor BVP setup.
    #[test]
    fn story_compact_hmx_setup_uses_the_default_facade() {
        let mut reactor = compact_hmx_reactor();
        reactor
            .setup_bvp()
            .expect("compact HMX setup should succeed");
        let initial_guess = DMatrix::from_element(reactor.solver.unknowns.len(), 4, 0.5);
        let config = NrbvpHandoffConfig::new(
            initial_guess,
            0.0,
            1.0,
            4,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Sparse".to_string(),
            1e-8,
            None,
            50,
            None,
            Some("info".to_string()),
            false,
        )
        .with_solver_backend_config(ReactorBvpSolverConfig::default_lambdify());
        let backend = reactor
            .solver
            .build_nrbvp_backend(config)
            .expect("backend snapshot should build");

        assert_eq!(
            backend.generated_backend_config().symbolic_assembly_backend,
            BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            backend.generated_backend_config().matrix_backend_override,
            Some(MatrixBackend::Banded)
        );
        assert_eq!(
            reactor.solver.unknowns,
            vec!["Teta", "q", "C0", "J0", "C1", "J1"]
        );
        assert_eq!(reactor.kindata.state, KinDataState::Analyzed);
    }

    /// Story-level handoff test that keeps the solver facade and the backend snapshot aligned.
    #[test]
    fn story_build_nrbvp_backend_uses_facade_generated_options() {
        let mut reactor = compact_hmx_reactor();
        reactor
            .setup_bvp()
            .expect("compact HMX setup should succeed");
        let initial_guess = DMatrix::from_element(reactor.solver.unknowns.len(), 4, 0.5);
        let config = NrbvpHandoffConfig::new(
            initial_guess,
            0.0,
            1.0,
            4,
            "forward".to_string(),
            "Damped".to_string(),
            None,
            None,
            "Sparse".to_string(),
            1e-8,
            None,
            50,
            None,
            Some("info".to_string()),
            false,
        )
        .with_solver_backend_config(ReactorBvpSolverConfig::default_lambdify());

        let backend = reactor
            .solver
            .build_nrbvp_backend(config)
            .expect("backend snapshot should build");

        assert_eq!(
            backend.generated_backend_config().symbolic_assembly_backend,
            BvpSymbolicAssemblyBackend::AtomView
        );
        assert_eq!(
            backend.generated_backend_config().matrix_backend_override,
            Some(MatrixBackend::Banded)
        );
        assert_eq!(reactor.kindata.state, KinDataState::Analyzed);
    }
}
