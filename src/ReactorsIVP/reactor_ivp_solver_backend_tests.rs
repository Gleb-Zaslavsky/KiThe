#[cfg(test)]
mod tests {
    use super::super::SimpleReactorIVP::*;
    use super::super::solver_backend::{
        ReactorIvpExecutionBackend, ReactorIvpMatrixBackend, ReactorIvpMethod,
        ReactorIvpSolverConfig, ReactorIvpSymbolicBackend,
    };
    use crate::Kinetics::mechfinder_api::ReactionData;
    use RustedSciThe::numerical::LSODE2::{
        Lsode2AotProfile, Lsode2AotToolchain, Lsode2ControllerConfig, Lsode2LinearSolverPolicy,
        Lsode2LinearSystemStructure, Lsode2Method, Lsode2NativeExecutionConfig,
        Lsode2ResidualJacobianSource, Lsode2StopComparator, Lsode2SymbolicAssemblyBackend,
        Lsode2SymbolicExecutionMode,
    };
    use std::collections::HashMap;

    fn build_condensed_task() -> SimpleReactorTask {
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
        task
    }

    #[test]
    fn test_default_solver_config_builds_lambdify_atomview_sparse_problem_config() {
        let mut task = build_condensed_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let config = task
            .build_lsode2_problem_config()
            .expect("problem config should build");

        assert_eq!(config.method, Lsode2Method::Bdf);
        assert_eq!(
            config.controller,
            Lsode2ControllerConfig::automatic_adams_bdf()
        );
        assert_eq!(
            config.residual_jacobian_source,
            Lsode2ResidualJacobianSource::Symbolic {
                assembly: Lsode2SymbolicAssemblyBackend::AtomView,
                execution: Lsode2SymbolicExecutionMode::LambdifyExpr,
            }
        );
        assert_eq!(
            config.linear_system_structure,
            Lsode2LinearSystemStructure::Sparse
        );
        assert_eq!(config.linear_solver_policy, Lsode2LinearSolverPolicy::Auto);
        assert_eq!(
            config.native_execution,
            Lsode2NativeExecutionConfig::faithful_bdf_solve(200_000, 200_000)
        );
        assert_eq!(config.values, vec!["Teta", "q", "C0", "C1"]);
        assert_eq!(config.arg, "x");
        assert_eq!(config.y0.len(), 4);
        assert_eq!(
            task.solver_backend_config.execution_backend,
            ReactorIvpExecutionBackend::Lambdify
        );
    }

    #[test]
    fn test_method_selection_maps_to_controller_modes() {
        assert_eq!(
            ReactorIvpMethod::Auto.to_controller(),
            Lsode2ControllerConfig::automatic_adams_bdf()
        );
        assert_eq!(
            ReactorIvpMethod::Bdf.to_controller(),
            Lsode2ControllerConfig::bdf_only()
        );
        assert_eq!(
            ReactorIvpMethod::Adams.to_controller(),
            Lsode2ControllerConfig::adams_only()
        );
        assert_eq!(ReactorIvpMethod::Auto.to_lsode2_method(), Lsode2Method::Bdf);
        assert_eq!(ReactorIvpMethod::Bdf.to_lsode2_method(), Lsode2Method::Bdf);
        assert_eq!(
            ReactorIvpMethod::Adams.to_lsode2_method(),
            Lsode2Method::Bdf
        );
    }

    #[test]
    fn test_aot_solver_config_preserves_typed_route_selection() {
        let mut task = build_condensed_task();
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default()
                .with_method(ReactorIvpMethod::Adams)
                .with_symbolic_backend(ReactorIvpSymbolicBackend::ExprLegacy)
                .with_matrix_backend(ReactorIvpMatrixBackend::Banded)
                .with_aot_backend(Lsode2AotToolchain::CTcc, Lsode2AotProfile::Debug)
                .with_integration_domain(0.0, 0.5)
                .with_first_step(Some(1.0e-4))
                .with_max_step(1.0e-2)
                .with_rtol(1.0e-7)
                .with_atol(1.0e-9),
        );
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let config = task
            .build_lsode2_problem_config()
            .expect("problem config should build");

        assert_eq!(config.controller, Lsode2ControllerConfig::adams_only());
        assert_eq!(
            config.residual_jacobian_source,
            Lsode2ResidualJacobianSource::Symbolic {
                assembly: Lsode2SymbolicAssemblyBackend::ExprLegacy,
                execution: Lsode2SymbolicExecutionMode::Aot {
                    toolchain: Lsode2AotToolchain::CTcc,
                    profile: Lsode2AotProfile::Debug,
                },
            }
        );
        assert_eq!(
            config.linear_system_structure,
            Lsode2LinearSystemStructure::Banded { kl: 0, ku: 0 }
        );
        assert_eq!(
            config.native_execution,
            Lsode2NativeExecutionConfig::faithful_bdf_solve(200_000, 200_000)
        );
        assert_eq!(config.t0, 0.0);
        assert_eq!(config.t_bound, 0.5);
        assert_eq!(config.first_step, Some(1.0e-4));
        assert_eq!(config.max_step, 1.0e-2);
        assert_eq!(config.rtol, 1.0e-7);
        assert_eq!(config.atol, 1.0e-9);
    }

    #[test]
    fn test_lambdify_solver_config_stays_on_symbolic_execution_route() {
        let mut task = build_condensed_task();
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default()
                .with_method(ReactorIvpMethod::Auto)
                .with_execution_backend(ReactorIvpExecutionBackend::Lambdify),
        );
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let config = task
            .build_lsode2_problem_config()
            .expect("problem config should build");

        assert_eq!(
            config.residual_jacobian_source,
            Lsode2ResidualJacobianSource::Symbolic {
                assembly: Lsode2SymbolicAssemblyBackend::AtomView,
                execution: Lsode2SymbolicExecutionMode::LambdifyExpr,
            }
        );
        assert!(!matches!(
            config.residual_jacobian_source,
            Lsode2ResidualJacobianSource::Symbolic {
                execution: Lsode2SymbolicExecutionMode::Aot { .. },
                ..
            }
        ));
    }

    #[test]
    fn test_problem_config_rejects_mismatched_dimensions() {
        let task = build_condensed_task();
        let result = task.solver_backend_config.to_rusted_problem_config(
            vec![],
            vec!["x".to_string()],
            "t".to_string(),
            nalgebra::DVector::from_vec(vec![1.0]),
        );
        assert!(matches!(
            result,
            Err(super::super::solver_backend::ReactorIvpBackendError::InvalidConfiguration(_))
        ));
    }

    #[test]
    fn test_banded_backend_accepts_inferred_band_route_for_nonempty_state_count() {
        let result = ReactorIvpMatrixBackend::Banded.validate_for_dimension(2);
        assert!(result.is_ok());
    }

    #[test]
    fn test_banded_backend_maps_to_inferred_lsode2_band_structure() {
        let backend = ReactorIvpMatrixBackend::Banded;

        assert!(backend.validate_for_dimension(2).is_ok());
        assert!(backend.is_banded());
        assert_eq!(
            backend.to_rusted(),
            Lsode2LinearSystemStructure::Banded { kl: 0, ku: 0 }
        );
    }

    #[test]
    fn test_stop_condition_maps_to_species_concentration_le_rule() {
        let mut task = build_condensed_task();
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default().with_stop_condition_le(0, 1.0e-4),
        );
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let config = task
            .build_lsode2_problem_config()
            .expect("problem config should build");

        assert_eq!(config.stop_conditions.len(), 1);
        let stop = &config.stop_conditions[0];
        assert_eq!(stop.variable, "C0");
        assert_eq!(stop.target, 1.0e-4);
        assert_eq!(stop.comparator, Lsode2StopComparator::LessEqual);
    }

    #[test]
    fn test_stop_condition_rejects_out_of_range_species_index() {
        let mut task = build_condensed_task();
        task.set_solver_backend_config(
            ReactorIvpSolverConfig::default().with_stop_condition_le(99, 1.0e-4),
        );
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let result = task.build_lsode2_problem_config();
        assert!(result.is_err());
    }
}
