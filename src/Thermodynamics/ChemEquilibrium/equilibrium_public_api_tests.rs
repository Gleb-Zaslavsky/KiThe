//! Contract tests for the narrow production equilibrium facade.
//!
//! These tests deliberately import only [`super::prelude`]. They prevent a
//! public request, result, or error from exposing an otherwise unnameable
//! implementation type and exercise the offline top-level workflow as an
//! external consumer would.

#[cfg(test)]
mod tests {
    use super::super::prelude::*;

    const PRESSURE_PA: f64 = 101_325.0;
    const TEMPERATURE_K: f64 = 500.0;

    fn offline_gas_spec() -> SubstanceSystemSpec {
        SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec![
            "N2".to_string(),
            "O2".to_string(),
        ]))
        .with_library_priorities(vec!["NASA_gas".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("the public offline gas specification must validate")
    }

    fn legacy_nr_options() -> EquilibriumSolveOptions {
        EquilibriumSolveOptions::new()
            .with_solver_policy(SolverPolicy::Single(SolverBackend::Legacy(
                LegacyEquilibriumSolver::NR,
            )))
            .expect("the public single-backend policy must validate")
    }

    #[test]
    fn production_prelude_solves_fixed_and_one_point_range_consistently() {
        let conditions =
            EquilibriumConditions::new(TEMPERATURE_K, PRESSURE_PA, PRESSURE_PA).unwrap();
        let initial_moles = vec![0.79, 0.21];
        let options = legacy_nr_options();

        let fixed = PhaseEquilibriumPipelineRequest::new(
            offline_gas_spec(),
            initial_moles.clone(),
            conditions,
        )
        .with_solve_options(options.clone())
        .solve()
        .expect("the public fixed-P,T pipeline must solve");
        let range =
            PhaseEquilibriumPipelineRequest::new(offline_gas_spec(), initial_moles, conditions)
                .with_solve_options(options)
                .solve_temperature_range(TemperatureGrid::new(vec![TEMPERATURE_K]).unwrap())
                .expect("the public one-point range pipeline must solve");

        let fixed_solution: &MultiphaseEquilibriumSolution = fixed.solution();
        let ranged_solution = range.points()[0].solution();
        assert_eq!(
            fixed_solution.solve_report().accepted_backend,
            SolverBackend::Legacy(LegacyEquilibriumSolver::NR)
        );
        assert_eq!(
            fixed_solution.component_moles().len(),
            ranged_solution.component_moles().len()
        );
        for (&fixed_moles, &ranged_moles) in fixed_solution
            .component_moles()
            .iter()
            .zip(ranged_solution.component_moles())
        {
            let scale = fixed_moles.abs().max(ranged_moles.abs()).max(1e-30);
            assert!((fixed_moles - ranged_moles).abs() / scale <= 1e-10);
        }

        let nitrogen = PhaseComponentId::new(PhaseId::new(None), "N2");
        let nitrogen_moles = fixed_solution
            .moles_for(&nitrogen)
            .expect("phase-qualified public lookup must find N2");
        assert!(nitrogen_moles.is_finite() && nitrogen_moles > 0.0);
        assert!(
            fixed_solution
                .accepted_solution()
                .validation()
                .residual_l2_norm
                .is_finite()
        );
        assert!(
            fixed_solution
                .accepted_solution()
                .validation()
                .max_abs_element_balance_error
                <= 1e-10
        );
    }

    #[test]
    fn production_prelude_rejects_duplicate_sparse_inventory_transactionally() {
        let duplicate = PhaseComponentId::new(PhaseId::new(None), "N2");
        let request = PhaseEquilibriumPipelineRequest::new_with_sparse_initial_composition(
            offline_gas_spec(),
            vec![(duplicate.clone(), 0.79), (duplicate, 0.21)],
            EquilibriumConditions::new(TEMPERATURE_K, PRESSURE_PA, PRESSURE_PA).unwrap(),
        )
        .with_solve_options(legacy_nr_options());

        let error = request
            .solve()
            .expect_err("duplicate phase-qualified inventory must not publish a result");
        assert!(matches!(
            &error,
            PhaseEquilibriumPipelineError::Solve(ReactionExtentError::InvalidProblem {
                field: "initial_composition",
                ..
            })
        ));
        assert_eq!(
            match error {
                PhaseEquilibriumPipelineError::Solve(error) => error.kind(),
                _ => panic!("duplicate inventory must be rejected by typed solve validation"),
            },
            ReactionExtentErrorKind::Formulation
        );
    }

    #[test]
    fn production_prelude_can_describe_element_selection_and_physical_phases() {
        let policy = EquilibriumCandidatePolicy::new(ElementSearchMode::ExactSet);
        assert_eq!(policy.element_mode(), ElementSearchMode::ExactSet);

        let path = PhSolvePath::MonolithicFixedActiveSet;
        assert_eq!(path, PhSolvePath::MonolithicFixedActiveSet);

        let solid = EquilibriumCandidatePhaseAssignment::pure_condensed(
            PhaseId::new(Some("solid".to_string())),
            PhysicalState::Solid,
            vec!["C(gr)".to_string()],
        );
        assert_eq!(solid.physical_state(), PhysicalState::Solid);
        assert_eq!(solid.model(), PhaseModel::PureCondensed);
    }
}
