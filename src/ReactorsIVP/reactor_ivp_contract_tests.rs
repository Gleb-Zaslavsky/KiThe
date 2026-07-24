#[cfg(test)]
mod tests {
    use super::super::SimpleReactorIVP::*;
    use crate::Kinetics::mechfinder_api::ReactionData;
    use approx::assert_abs_diff_eq;
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
    fn test_condensed_setup_canonicalizes_spatial_coordinate() {
        let mut task = build_condensed_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        assert_eq!(task.solver.arg_name, "x");
        assert_eq!(
            task.solver.x_range,
            (
                task.solver_backend_config.x0,
                task.solver_backend_config.x_bound
            )
        );

        let inputs = task
            .condensed_inputs()
            .expect("typed condensed inputs should build");
        assert_eq!(inputs.physical.ro, 1200.0);
        assert_eq!(inputs.initial_state.temperature, 450.0);
        assert_eq!(inputs.initial_state.heat_flux, 0.15);
        assert_eq!(inputs.initial_state.species.get("A"), Some(&0.7));
        assert_eq!(inputs.initial_state.species.get("B"), Some(&0.3));
    }

    #[test]
    fn test_canonical_state_names_follow_thermal_then_species_order() {
        let task = build_condensed_task();
        assert_eq!(
            task.canonical_state_names(),
            vec![
                "Teta".to_string(),
                "q".to_string(),
                "C0".to_string(),
                "C1".to_string()
            ]
        );
    }

    #[test]
    fn test_set_density_rejects_nonfinite_and_negative_values() {
        let mut task = SimpleReactorTask::new();
        assert!(task.set_density(1000.0).is_ok());
        assert!(task.set_density(0.0).is_err());
        assert!(task.set_density(-1.0).is_err());
        assert!(task.set_density(f64::NAN).is_err());
        assert_abs_diff_eq!(task.ro, 1000.0);
    }

    #[test]
    fn test_build_initial_state_vector_uses_canonical_order() {
        let task = build_condensed_task();
        let y0 = task
            .build_initial_state_vector()
            .expect("initial state vector should build");

        assert_eq!(y0.len(), 4);
        assert_abs_diff_eq!(y0[0], (450.0 - 100.0) / 100.0, epsilon = 1e-12);
        assert_abs_diff_eq!(y0[1], 0.15, epsilon = 1e-12);
        assert_abs_diff_eq!(y0[2], 0.7, epsilon = 1e-12);
        assert_abs_diff_eq!(y0[3], 0.3, epsilon = 1e-12);
    }

    #[test]
    fn test_setup_condensed_ivp_builds_species_only_snapshot() {
        let mut task = build_condensed_task();
        let result = task.setup_condensed_ivp();
        assert!(result.is_ok(), "{result:?}");

        assert_eq!(
            task.solver.unknowns,
            vec![
                "Teta".to_string(),
                "q".to_string(),
                "C0".to_string(),
                "C1".to_string()
            ]
        );
        assert_eq!(task.solver.eq_system.len(), 4);
        assert_eq!(task.map_of_equations.len(), 4);
        assert!(task.map_of_equations.contains_key("Teta"));
        assert!(task.map_of_equations.contains_key("q"));
        assert!(task.map_of_equations.contains_key("A"));
        assert!(task.map_of_equations.contains_key("B"));
        assert!(
            !task
                .solver
                .unknowns
                .iter()
                .any(|name| name.starts_with('J'))
        );
        assert!(
            !task
                .map_of_equations
                .keys()
                .any(|key| key.ends_with("_flux"))
        );
        assert!(
            !task
                .map_of_equations
                .values()
                .any(|(_, expr)| expr.to_string().contains("Pe_D"))
        );
    }

    #[test]
    fn test_setup_condensed_ivp_publishes_canonical_solver_shape() {
        let mut task = build_condensed_task();
        let result = task.setup_condensed_ivp();
        assert!(result.is_ok(), "{result:?}");

        assert_eq!(task.solver.arg_name, "x");
        assert_eq!(task.solver.unknowns.len(), 4);
        assert_eq!(task.solver.unknowns[0], "Teta");
        assert_eq!(task.solver.unknowns[1], "q");
    }

    #[test]
    fn test_condensed_inputs_bundle_captures_physical_and_initial_state_contract() {
        let task = build_condensed_task();
        let inputs = task
            .condensed_inputs()
            .expect("typed condensed inputs should build");

        assert_eq!(inputs.physical.ro, 1200.0);
        assert_eq!(inputs.physical.Cp, 1000.0);
        assert_eq!(inputs.physical.Lambda, 0.25);
        assert_eq!(inputs.physical.m, 0.015);
        assert_eq!(inputs.physical.L, 0.05);
        assert_eq!(inputs.initial_state.temperature, 450.0);
        assert_eq!(inputs.initial_state.heat_flux, 0.15);
        assert_eq!(inputs.initial_state.species.get("A"), Some(&0.7));
        assert_eq!(inputs.initial_state.species.get("B"), Some(&0.3));
        assert_eq!(inputs.thermal_effects, vec![-2.0e5]);
    }

    #[test]
    fn test_typed_initial_state_rejects_missing_species_entry() {
        let task = build_condensed_task();
        let mut inputs = task
            .condensed_inputs()
            .expect("typed condensed inputs should build");
        inputs.initial_state.species.remove("B");

        let result = inputs.validate(&task.kindata.substances);
        assert!(matches!(result, Err(IvpError::MissingData(_))));
    }

    #[test]
    fn test_validate_condensed_phase_contract_rejects_thermal_effect_mismatch() {
        let mut task = build_condensed_task();
        task.thermal_effects = vec![-2.0e5, 1.0];

        let result = task.validate_condensed_phase_contract();

        assert!(matches!(result, Err(IvpError::InvalidConfiguration(_))));
    }

    #[test]
    fn test_a_to_b_reaction_has_expected_species_source_signs() {
        let mut task = build_condensed_task();
        task.setup_condensed_ivp()
            .expect("condensed setup should succeed");

        let config = task
            .build_lsode2_problem_config()
            .expect("problem config should build");
        let names = config
            .values
            .iter()
            .map(|name| name.as_str())
            .collect::<Vec<_>>();
        let y0 = config.y0.as_slice();

        let rhs_a = config.eq_system[2].lower_to_linear(&names).eval(y0);
        let rhs_b = config.eq_system[3].lower_to_linear(&names).eval(y0);

        assert!(
            rhs_a > 0.0,
            "A should follow the current source convention, got rhs = {rhs_a}"
        );
        assert!(
            rhs_b < 0.0,
            "B should follow the current source convention, got rhs = {rhs_b}"
        );
    }
}
