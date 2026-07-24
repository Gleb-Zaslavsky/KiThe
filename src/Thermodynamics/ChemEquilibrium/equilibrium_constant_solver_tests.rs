//! Regression tests for the independent equilibrium-constant extent solver.
//!
//! The first pass deliberately covers only one-reaction systems. These tests
//! keep the small-system validator honest without tying it to the canonical
//! log-moles residual path.

#[cfg(test)]
mod tests {
    use super::super::equilibrium_constant_problem::{
        EquilibriumConstantActivityModel, EquilibriumConstantProblem, MOLAR_GAS_CONSTANT,
    };
    use super::super::equilibrium_constant_solver::{
        EquilibriumConstantSolver, EquilibriumConstantSolverMode,
    };
    use super::super::equilibrium_log_moles::GibbsFn;
    use super::super::equilibrium_nonlinear::ReactionExtentError;
    use super::super::equilibrium_reaction_basis::{
        ReactionBasisTolerances, ValidatedReactionBasis,
    };
    use nalgebra::DMatrix;
    use std::rc::Rc;

    fn constant_gibbs(value: f64) -> GibbsFn {
        Rc::new(move |_| value)
    }

    fn dissociation_problem(target_k: f64) -> EquilibriumConstantProblem {
        let temperature = 1_500.0;
        let delta_g = -target_k.ln() * MOLAR_GAS_CONSTANT * temperature;
        let basis = ValidatedReactionBasis::new(
            vec!["A2".to_string(), "A".to_string()],
            &DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            DMatrix::from_column_slice(2, 1, &[-1.0, 2.0]),
            1,
            ReactionBasisTolerances::default(),
        )
        .unwrap();

        EquilibriumConstantProblem::new(
            basis,
            vec![1.0, 0.0],
            vec![constant_gibbs(-delta_g), constant_gibbs(0.0)],
            super::super::equilibrium_problem::EquilibriumConditions::new(
                temperature,
                2.0 * 101_325.0,
                101_325.0,
            )
            .unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap()
    }

    /// Independent algebraic reference for reactions of the form
    /// `2 A <=> a B + b C` with `(a, b) = (2, 1)` or `(1, 2)`.
    ///
    /// Starting from one mole of `A`, both variants have
    /// `Q = 4 (P/P_ref) eta^3 / ((1 + eta) (1 - 2 eta)^2)`.
    /// Bisection here deliberately avoids the production quotient and solver.
    fn ternary_dissociation_extent(target_k: f64, pressure_ratio: f64) -> f64 {
        let residual = |eta: f64| {
            4.0 * pressure_ratio * eta.powi(3) - target_k * (1.0 + eta) * (1.0 - 2.0 * eta).powi(2)
        };
        let (mut lower, mut upper) = (0.0, 0.5);
        for _ in 0..100 {
            let midpoint = 0.5 * (lower + upper);
            if residual(midpoint) > 0.0 {
                upper = midpoint;
            } else {
                lower = midpoint;
            }
        }
        0.5 * (lower + upper)
    }

    fn ternary_dissociation_problem(
        species: [&str; 3],
        product_coefficients: [f64; 2],
        target_k: f64,
        pressure_ratio: f64,
    ) -> EquilibriumConstantProblem {
        let [product_b, product_c] = product_coefficients;
        let basis = ValidatedReactionBasis::new(
            species.into_iter().map(str::to_string).collect(),
            // Two abstract elements make the conservation rank explicit:
            // A contains half of each product inventory.
            &DMatrix::from_row_slice(
                3,
                2,
                &[product_b / 2.0, product_c / 2.0, 1.0, 0.0, 0.0, 1.0],
            ),
            DMatrix::from_column_slice(3, 1, &[-2.0, product_b, product_c]),
            2,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let temperature = 2_000.0;
        let delta_g = -target_k.ln() * MOLAR_GAS_CONSTANT * temperature;

        EquilibriumConstantProblem::new(
            basis,
            vec![1.0, 0.0, 0.0],
            vec![
                constant_gibbs(-delta_g / 2.0),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
            ],
            super::super::equilibrium_problem::EquilibriumConditions::new(
                temperature,
                pressure_ratio * 101_325.0,
                101_325.0,
            )
            .unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap()
    }

    #[test]
    fn one_reaction_dissociation_solver_matches_the_analytic_fixture() {
        let target_k = 0.75_f64;
        let problem = dissociation_problem(target_k);
        let solver = EquilibriumConstantSolver::default();

        let solution = solver.solve(&problem).unwrap();
        let pressure_ratio =
            problem.conditions().pressure() / problem.conditions().reference_pressure();
        let expected_product_mole = (target_k / (target_k + 4.0 * pressure_ratio)).sqrt();
        let expected_extent = expected_product_mole;

        assert!((solution.extent - expected_extent).abs() < 1e-8);
        assert!((solution.moles[0] - (1.0 - expected_product_mole)).abs() < 1e-8);
        assert!((solution.moles[1] - 2.0 * expected_product_mole).abs() < 1e-8);
        assert!(solution.report.converged);
        assert!(solution.validation.accepted);
        assert!(solution.validation.max_abs_log_residual < 1e-8);
    }

    #[test]
    fn nitric_oxide_solver_matches_the_old_closed_form_extent() {
        let target_k = 1e-3_f64;
        let sqrt_k = target_k.sqrt();
        let extent = sqrt_k / (1.0 + 2.0 * sqrt_k);
        let basis = ValidatedReactionBasis::new(
            vec!["NO".to_string(), "N2".to_string(), "O2".to_string()],
            &DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 2.0, 0.0, 0.0, 2.0]),
            DMatrix::from_column_slice(3, 1, &[-2.0, 1.0, 1.0]),
            2,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let temperature = 4_500.0;
        let delta_g = -target_k.ln() * MOLAR_GAS_CONSTANT * temperature;
        let problem = EquilibriumConstantProblem::new(
            basis,
            vec![1.0, 0.0, 0.0],
            vec![
                constant_gibbs(-delta_g / 2.0),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
            ],
            super::super::equilibrium_problem::EquilibriumConditions::new(
                temperature,
                101_325.0,
                101_325.0,
            )
            .unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();

        let solution = EquilibriumConstantSolver::default()
            .solve(&problem)
            .unwrap();
        assert!((solution.extent - 2.0 * extent).abs() < 1e-8);
        assert!((solution.moles[0] - (1.0 - 2.0 * extent)).abs() < 1e-8);
        assert!((solution.moles[1] - extent).abs() < 1e-8);
        assert!((solution.moles[2] - extent).abs() < 1e-8);
        assert!(solution.validation.accepted);
    }

    #[test]
    fn ternary_dissociation_fixtures_preserve_pressure_and_stoichiometric_ordering() {
        // These two contracts are retained for regression coverage.
        // Their quotient algebra is identical, while the product ordering and
        // reconstructed mole vectors differ.
        for (species, coefficients) in [
            (["N2O", "N2", "O2"], [2.0, 1.0]),
            (["NO2", "N2", "O2"], [1.0, 2.0]),
        ] {
            for (target_k, pressure_ratio) in [(1e-4, 0.5), (0.2, 1.0), (5.0, 3.0)] {
                let problem =
                    ternary_dissociation_problem(species, coefficients, target_k, pressure_ratio);
                let solution = EquilibriumConstantSolver::default()
                    .solve(&problem)
                    .unwrap();
                let expected_extent = ternary_dissociation_extent(target_k, pressure_ratio);

                // Basis normalization maps [-2, a, b] to [-1, a/2, b/2],
                // so its coordinate is twice the physical reaction extent.
                assert!((solution.extent - 2.0 * expected_extent).abs() < 1e-8);
                assert!((solution.moles[0] - (1.0 - 2.0 * expected_extent)).abs() < 1e-8);
                assert!((solution.moles[1] - coefficients[0] * expected_extent).abs() < 1e-8);
                assert!((solution.moles[2] - coefficients[1] * expected_extent).abs() < 1e-8);
                assert!(solution.validation.accepted);
                assert!(solution.validation.max_abs_log_residual < 1e-8);
            }
        }
    }

    #[test]
    fn multi_reaction_problem_is_reported_as_not_applicable() {
        let basis = ValidatedReactionBasis::new(
            vec![
                "NO".to_string(),
                "N2".to_string(),
                "O2".to_string(),
                "NO2".to_string(),
            ],
            &DMatrix::from_row_slice(
                4,
                2,
                &[
                    1.0, 1.0, // NO
                    2.0, 0.0, // N2
                    0.0, 2.0, // O2
                    1.0, 2.0, // NO2
                ],
            ),
            DMatrix::from_column_slice(
                4,
                2,
                &[
                    -2.0, 1.0, 1.0, 0.0, // 2NO <=> N2 + O2
                    -4.0, 1.0, 0.0, 2.0, // 4NO <=> N2 + 2NO2
                ],
            ),
            2,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let problem = EquilibriumConstantProblem::new(
            basis,
            vec![1.0, 0.0, 0.0, 0.0],
            vec![
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
            ],
            super::super::equilibrium_problem::EquilibriumConditions::new(
                1_500.0, 101_325.0, 101_325.0,
            )
            .unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();

        let error = EquilibriumConstantSolver::default()
            .solve(&problem)
            .unwrap_err();
        assert!(matches!(
            error,
            ReactionExtentError::ValidationNotApplicable { .. }
        ));
    }

    #[test]
    fn solve_if_applicable_returns_some_for_supported_one_reaction_systems() {
        let problem = dissociation_problem(0.75);
        let solver = EquilibriumConstantSolver::default();

        let result = solver.solve_if_applicable(&problem).unwrap();

        assert!(result.is_some());
        let solution = result.unwrap();
        assert!(solution.report.converged);
        assert!(solution.validation.accepted);
    }

    #[test]
    fn solve_if_applicable_respects_required_mode_for_unsupported_systems() {
        let basis = ValidatedReactionBasis::new(
            vec![
                "NO".to_string(),
                "N2".to_string(),
                "O2".to_string(),
                "NO2".to_string(),
            ],
            &DMatrix::from_row_slice(
                4,
                2,
                &[
                    1.0, 1.0, // NO
                    2.0, 0.0, // N2
                    0.0, 2.0, // O2
                    1.0, 2.0, // NO2
                ],
            ),
            DMatrix::from_column_slice(
                4,
                2,
                &[
                    -2.0, 1.0, 1.0, 0.0, // 2NO <=> N2 + O2
                    -4.0, 1.0, 0.0, 2.0, // 4NO <=> N2 + 2NO2
                ],
            ),
            2,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let problem = EquilibriumConstantProblem::new(
            basis,
            vec![1.0, 0.0, 0.0, 0.0],
            vec![
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
                constant_gibbs(0.0),
            ],
            super::super::equilibrium_problem::EquilibriumConditions::new(
                1_500.0, 101_325.0, 101_325.0,
            )
            .unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();

        let solver = EquilibriumConstantSolver {
            mode: EquilibriumConstantSolverMode::Required,
            ..EquilibriumConstantSolver::default()
        };
        let error = solver.solve_if_applicable(&problem).unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::ValidationNotApplicable { .. }
        ));
    }

    #[test]
    fn solver_report_has_stable_summary_rows_and_display() {
        let problem = dissociation_problem(0.75);
        let solution = EquilibriumConstantSolver::default()
            .solve(&problem)
            .unwrap();
        let rows = solution.report.summary_rows();

        assert!(
            rows.iter()
                .any(|row| row.section == "solve" && row.label == "converged")
        );
        assert!(
            rows.iter()
                .any(|row| row.section == "solve" && row.label == "iterations")
        );
        assert!(
            rows.iter()
                .any(|row| row.section == "solve" && row.label == "bracketed")
        );
        let rendered = format!("{}", solution.report);
        assert!(rendered.contains("[solve] converged = true"));
        assert!(rendered.contains("[solve] iterations = "));
        assert!(rendered.contains("[solve] log_residual = "));
    }
}
