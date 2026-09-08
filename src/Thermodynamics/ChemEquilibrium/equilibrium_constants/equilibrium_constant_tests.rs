//! Contract tests for the independent equilibrium-constant foundation.
//!
//! The fixtures are adapted from earlier diatomic dissociation and
//! nitric-oxide decomposition studies: `A2 <=> 2A` and `2NO <=> N2 + O2`.
//! They verify reaction ordering, elemental conservation,
//! pressure-aware ideal-gas quotients, and the independent `ln(Q)-ln(K)` report.

#[cfg(test)]
mod tests {
    use super::super::equilibrium_constant_problem::{
        EquilibriumConstantActivityModel, EquilibriumConstantProblem, MOLAR_GAS_CONSTANT,
    };
    use super::super::equilibrium_constant_validation::{
        EquilibriumConstantValidationTolerances, validate_equilibrium_constants,
    };
    use super::super::equilibrium_log_moles::GibbsFn;
    use super::super::equilibrium_log_moles::{Phase, PhaseKind};
    use super::super::equilibrium_problem::{
        EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
    };
    use super::super::equilibrium_reaction_basis::{
        ReactionBasisTolerances, ValidatedReactionBasis,
    };
    use nalgebra::DMatrix;
    use std::rc::Rc;

    fn constant_gibbs(value: f64) -> GibbsFn {
        Rc::new(move |_| value)
    }

    fn dissociation_basis() -> ValidatedReactionBasis {
        // Species order [A2, A], one conserved element, reaction -A2 + 2A.
        ValidatedReactionBasis::new(
            vec!["A2".to_string(), "A".to_string()],
            &DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            DMatrix::from_column_slice(2, 1, &[-1.0, 2.0]),
            1,
            ReactionBasisTolerances::default(),
        )
        .unwrap()
    }

    #[test]
    fn reaction_basis_is_conservative_and_deterministically_oriented() {
        let basis = ValidatedReactionBasis::new(
            vec!["A2".to_string(), "A".to_string()],
            &DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            DMatrix::from_column_slice(2, 1, &[0.5, -1.0]),
            1,
            ReactionBasisTolerances::default(),
        )
        .unwrap();

        assert_eq!(basis.reactions().column(0).as_slice(), &[-1.0, 2.0]);
        assert_eq!(basis.reaction_count(), 1);
        assert!(basis.max_conservation_residual() <= 1e-12);
    }

    #[test]
    fn reaction_basis_rejects_a_nonconservative_fixture() {
        let error = ValidatedReactionBasis::new(
            vec!["A2".to_string(), "A".to_string()],
            &DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            DMatrix::from_column_slice(2, 1, &[-1.0, 1.0]),
            1,
            ReactionBasisTolerances::default(),
        )
        .unwrap_err();

        assert!(error.to_string().contains("elemental conservation"));
    }

    #[test]
    fn reaction_basis_rejects_an_overdetermined_fixture() {
        let error = ValidatedReactionBasis::new(
            vec!["A2".to_string(), "A".to_string()],
            &DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            DMatrix::from_column_slice(2, 2, &[-1.0, 2.0, -1.0, 2.0]),
            1,
            ReactionBasisTolerances::default(),
        )
        .unwrap_err();

        assert!(error.to_string().contains("expected 1"));
    }

    #[test]
    fn reaction_basis_rejects_a_numerically_zero_reaction_column() {
        let error = ValidatedReactionBasis::new(
            vec!["A2".to_string(), "A".to_string()],
            &DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            DMatrix::from_column_slice(2, 1, &[0.0, 0.0]),
            1,
            ReactionBasisTolerances::default(),
        )
        .unwrap_err();

        assert!(error.to_string().contains("numerically zero"));
    }

    #[test]
    fn reaction_basis_is_stable_under_species_permutation_up_to_reordering() {
        let original = ValidatedReactionBasis::new(
            vec!["A2".to_string(), "B2".to_string(), "AB".to_string()],
            &DMatrix::from_row_slice(3, 2, &[2.0, 0.0, 0.0, 2.0, 1.0, 1.0]),
            DMatrix::from_column_slice(3, 1, &[-1.0, -1.0, 2.0]),
            2,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let permuted = ValidatedReactionBasis::new(
            vec!["AB".to_string(), "A2".to_string(), "B2".to_string()],
            &DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 2.0, 0.0, 0.0, 2.0]),
            DMatrix::from_column_slice(3, 1, &[2.0, -1.0, -1.0]),
            2,
            ReactionBasisTolerances::default(),
        )
        .unwrap();

        let permuted_back = [
            permuted.reactions()[(1, 0)],
            permuted.reactions()[(2, 0)],
            permuted.reactions()[(0, 0)],
        ];
        let original_column = original.reactions().column(0);
        let scale = permuted_back[0] / original_column[0];

        for (lhs, rhs) in original_column.iter().zip(permuted_back.iter()) {
            assert!((*rhs - scale * *lhs).abs() < 1e-12);
        }
        assert!(original.max_conservation_residual() <= 1e-12);
        assert!(permuted.max_conservation_residual() <= 1e-12);
    }

    #[test]
    fn standard_gibbs_is_converted_to_ln_k_without_main_residual_code() {
        let temperature = 2_000.0;
        let problem = EquilibriumConstantProblem::new(
            dissociation_basis(),
            vec![1.0, 0.0],
            vec![constant_gibbs(100_000.0), constant_gibbs(30_000.0)],
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();

        let reaction = problem.basis().reaction_id(0).unwrap();
        let expected_delta_g = 2.0 * 30_000.0 - 100_000.0;
        let expected_ln_k = -expected_delta_g / (MOLAR_GAS_CONSTANT * temperature);
        assert!(
            (problem.standard_reaction_gibbs(reaction).unwrap() - expected_delta_g).abs() < 1e-12
        );
        assert!((problem.ln_equilibrium_constant(reaction).unwrap() - expected_ln_k).abs() < 1e-12);
    }

    #[test]
    fn diatomic_dissociation_fixture_satisfies_pressure_aware_reaction_quotient() {
        let pressure = 2.0 * 101_325.0;
        let reference_pressure = 101_325.0;
        let target_k = 0.75_f64;
        let pressure_ratio = pressure / reference_pressure;
        let alpha = (target_k / (target_k + 4.0 * pressure_ratio)).sqrt();
        let candidate = vec![1.0 - alpha, 2.0 * alpha];
        let temperature = 1_500.0;
        let delta_g = -target_k.ln() * MOLAR_GAS_CONSTANT * temperature;
        let problem = EquilibriumConstantProblem::new(
            dissociation_basis(),
            vec![1.0, 0.0],
            vec![constant_gibbs(-delta_g), constant_gibbs(0.0)],
            EquilibriumConditions::new(temperature, pressure, reference_pressure).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();

        let report = validate_equilibrium_constants(
            &problem,
            &candidate,
            EquilibriumConstantValidationTolerances::default(),
        )
        .unwrap();
        assert!(report.accepted);
        assert!(report.max_abs_log_residual < 1e-12);
    }

    #[test]
    fn nitric_oxide_fixture_matches_the_old_closed_form_extent() {
        // Species [NO, N2, O2], reaction -2NO + N2 + O2.
        let basis = ValidatedReactionBasis::new(
            vec!["NO".to_string(), "N2".to_string(), "O2".to_string()],
            &DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 2.0, 0.0, 0.0, 2.0]),
            DMatrix::from_column_slice(3, 1, &[-2.0, 1.0, 1.0]),
            2,
            ReactionBasisTolerances::default(),
        )
        .unwrap();
        let target_k = 1e-3_f64;
        let sqrt_k = target_k.sqrt();
        let extent = sqrt_k / (1.0 + 2.0 * sqrt_k);
        let candidate = vec![1.0 - 2.0 * extent, extent, extent];
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
            EquilibriumConditions::new(temperature, 101_325.0, 101_325.0).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();

        let report = validate_equilibrium_constants(
            &problem,
            &candidate,
            EquilibriumConstantValidationTolerances::default(),
        )
        .unwrap();
        assert!(report.accepted);
        assert!(report.max_abs_log_residual < 1e-12);
        // Delta-nu is zero for 2NO -> N2 + O2, so total moles stay constant.
        assert!((candidate.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn quotient_rejects_zero_candidate_species_instead_of_returning_infinity() {
        let problem = EquilibriumConstantProblem::new(
            dissociation_basis(),
            vec![1.0, 0.0],
            vec![constant_gibbs(0.0), constant_gibbs(0.0)],
            EquilibriumConditions::new(1_000.0, 101_325.0, 101_325.0).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();
        let reaction = problem.basis().reaction_id(0).unwrap();
        let error = problem
            .ln_reaction_quotient(reaction, &[1.0, 0.0])
            .unwrap_err();
        assert!(error.to_string().contains("strictly positive"));
    }

    #[test]
    fn prepared_problem_adapter_shares_domain_data_but_builds_an_independent_contract() {
        let initial_moles = vec![1.0, 0.0];
        let source = EquilibriumProblem::new(
            vec!["A2".to_string(), "A".to_string()],
            initial_moles.clone(),
            LogMolesInitialGuess::from_initial_moles(&initial_moles).unwrap(),
            DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
            vec![constant_gibbs(10_000.0), constant_gibbs(20_000.0)],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            }],
            EquilibriumConditions::new(2_000.0, 101_325.0, 101_325.0).unwrap(),
        )
        .unwrap();
        let prepared = PreparedEquilibriumProblem::new(source).unwrap();

        let independent = EquilibriumConstantProblem::from_prepared_ideal_gas(
            &prepared,
            ReactionBasisTolerances::default(),
        )
        .unwrap();

        assert_eq!(independent.basis().species(), &["A2", "A"]);
        assert_eq!(independent.initial_moles(), &[1.0, 0.0]);
        assert_eq!(independent.basis().reaction_count(), 1);
        assert!(independent.basis().max_conservation_residual() < 1e-10);
    }

    #[test]
    fn validation_report_has_stable_summary_rows_and_display() {
        let problem = EquilibriumConstantProblem::new(
            dissociation_basis(),
            vec![1.0, 0.0],
            vec![constant_gibbs(100_000.0), constant_gibbs(30_000.0)],
            EquilibriumConditions::new(2_000.0, 101_325.0, 101_325.0).unwrap(),
            EquilibriumConstantActivityModel::IdealGas,
        )
        .unwrap();
        let report = validate_equilibrium_constants(
            &problem,
            &[0.5, 1.0],
            EquilibriumConstantValidationTolerances::default(),
        )
        .unwrap();
        let rows = report.summary_rows();

        assert!(
            rows.iter()
                .any(|row| row.section == "validation" && row.label == "temperature")
        );
        assert!(
            rows.iter()
                .any(|row| row.section == "validation" && row.label == "accepted")
        );
        assert!(
            rows.iter()
                .any(|row| row.section == "reaction" && row.label == "0")
        );
        let rendered = format!("{report}");
        assert!(rendered.contains("[validation] temperature = 2000.000000"));
        assert!(rendered.contains("[reaction] 0 = lnQ="));
    }
}
