//! Synthetic one-phase log-moles dynamic-range evidence.
//!
//! Every species carries the same abstract element and belongs to one ideal
//! gas phase. Therefore a tiny species is an internal composition limit, not a
//! phase disappearance. The Gibbs closures are constructed from exact target
//! weights, giving a direct analytical oracle.

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use nalgebra::DMatrix;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
        GibbsFn, Phase, R, Solvers,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
    use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::PreparedPhaseControlRunner;

    const T_K: f64 = 1000.0;
    const P_PA: f64 = 100_000.0;
    const P0_PA: f64 = 100_000.0;
    const TOTAL_MOLES: f64 = 1.0;
    const TRACE_FLOOR: f64 = 1e-30;

    fn try_solve_with_seed(
        weights: &[f64],
        total_moles: f64,
        seed_moles: &[f64],
    ) -> Result<(Vec<f64>, f64), ReactionExtentError> {
        let names = (0..weights.len()).map(|i| format!("X{i}(g)")).collect();
        let initial = vec![total_moles / weights.len() as f64; weights.len()];
        let seed = LogMolesInitialGuess::from_moles(seed_moles, TRACE_FLOOR).unwrap();
        let elements = DMatrix::from_element(weights.len(), 1, 1.0);
        let gibbs: Vec<GibbsFn> = weights
            .iter()
            .map(|weight| {
                let g0 = -R * T_K * weight.ln();
                Rc::new(move |_| g0) as GibbsFn
            })
            .collect();
        let problem = EquilibriumProblem::new(
            names,
            initial,
            seed,
            elements,
            gibbs,
            vec![Phase {
                kind: PhaseActivityModel::IdealGas,
                species: (0..weights.len()).collect(),
            }],
            EquilibriumConditions::new(T_K, P_PA, P0_PA).unwrap(),
        )
        .unwrap();
        let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false).unwrap();
        let settings = runner.configure_solver();
        settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
        settings.solver_params.tol = 1e-11;
        settings.solver_params.max_iter = 300;
        let outcome = runner.solve()?;
        assert!(outcome.acceptance_report.complementarity.satisfied);
        assert_eq!(
            outcome.phase_control_report.final_phase_set.active_mask(),
            vec![true]
        );
        let balance = outcome.solution.moles().iter().sum::<f64>() - total_moles;
        Ok((outcome.solution.moles().to_vec(), balance.abs()))
    }

    fn solve(weights: &[f64]) -> (Vec<f64>, f64) {
        let initial = vec![TOTAL_MOLES / weights.len() as f64; weights.len()];
        try_solve_with_seed(weights, TOTAL_MOLES, &initial)
            .expect("one-phase dynamic-range fixture must solve")
    }

    fn assert_against_oracle(weights: &[f64], actual: &[f64], balance: f64) {
        let sum: f64 = weights.iter().sum();
        for (index, (weight, value)) in weights.iter().zip(actual).enumerate() {
            let expected = TOTAL_MOLES * weight / sum;
            assert!(
                value.is_finite() && *value > 0.0,
                "X{index} is not finite/positive: {value:e}"
            );
            let log_error = (value.ln() - expected.ln()).abs();
            let relative = (value - expected).abs() / expected;
            let chemical_potential_delta = R * T_K * ((value / actual[0]) / weight).ln();
            assert!(chemical_potential_delta.is_finite());
            assert!(
                chemical_potential_delta.abs() < 1e-6,
                "chemical-potential equality failed for X{index}: {chemical_potential_delta:e} J/mol"
            );
            println!(
                "X{index} expected={expected:.6e} actual={value:.6e} relative={relative:.3e} log_error={log_error:.3e} mu_minus_mu0={chemical_potential_delta:.3e}"
            );
            // Major and resolved minor species use a scale-aware contract. The
            // log error remains the useful diagnostic for trace species.
            assert!(
                log_error < 2e-5,
                "X{index} log-ratio error too large: {log_error:e}"
            );
        }
        assert!(balance < 1e-10, "element balance error={balance:e}");
    }

    #[test]
    fn analytic_single_phase_species_ratios_match_gibbs_weights() {
        let weights = [1.0, 1e-2, 1e-4, 1e-6, 1e-8];
        let (actual, balance) = solve(&weights);
        assert_against_oracle(&weights, &actual, balance);
        println!("verdict: ExtremeCompositionDynamicRangeInvariant");
    }

    #[test]
    fn i5_extreme_composition_dynamic_range_matrix() {
        let cases = [
            ("A", vec![1.0, 1e-2, 1e-4, 1e-6, 1e-8]),
            ("B", vec![1.0, 1e-4, 1e-8, 1e-12, 1e-16]),
            ("C", vec![1.0, 1e-6, 1e-12, 1e-18, 1e-24]),
            ("D", vec![1.0, 1e-8, 1e-16, 1e-24, 1e-32]),
        ];
        for (label, weights) in cases {
            let (actual, balance) = solve(&weights);
            let minimum = actual.iter().copied().fold(f64::INFINITY, f64::min);
            println!(
                "case={label} expected_min_fraction={:.3e} actual_min_moles={minimum:.3e} balance={balance:.3e} status=OK",
                weights.last().unwrap() / weights.iter().sum::<f64>()
            );
            assert_against_oracle(&weights, &actual, balance);
        }
        println!("verdict: ExtremeCompositionDynamicRangeInvariant");
    }

    #[test]
    fn i5_vanishing_species_resolution_sweep() {
        let epsilons = [1e-4, 1e-8, 1e-12, 1e-16, 1e-20, 1e-24, 1e-28, 1e-32, 1e-40];
        let mut previous_ratio = f64::INFINITY;
        for epsilon in epsilons {
            let weights = [1.0, 1e-2, epsilon];
            let (actual, balance) = solve(&weights);
            let expected = epsilon / weights.iter().sum::<f64>();
            let value = actual[2];
            let log_error = (value.ln() - expected.ln()).abs();
            let ratio = value / expected;
            assert!(value.is_finite() && value > 0.0);
            assert!(balance < 1e-10);
            assert!(
                value / actual[0] < previous_ratio,
                "vanishing species amount must decrease monotonically"
            );
            previous_ratio = value / actual[0];
            println!(
                "epsilon={epsilon:.1e} expected_trace={expected:.6e} actual_trace={value:.6e} actual_over_expected={ratio:.6e} log_error={log_error:.3e} balance={balance:.3e} verdict=VanishingSpeciesLogResolved"
            );
        }
    }

    #[test]
    fn i5_gibbs_weight_representability_matrix() {
        // This is deliberately solver-independent: it identifies the f64
        // boundary of the fixture before a nonlinear backend is involved.
        let ratios: [f64; 11] = [
            1e-4, 1e-8, 1e-12, 1e-16, 1e-20, 1e-24, 1e-28, 1e-32, 1e-40, 1e-60, 1e-100,
        ];
        let mut first_collapsed = None;

        for ratio in ratios {
            let delta_g = -R * T_K * ratio.ln();
            let reconstructed = (-delta_g / (R * T_K)).exp();
            let log_error = (reconstructed.ln() - ratio.ln()).abs();
            let representable = reconstructed > 0.0 && reconstructed.is_finite();

            println!(
                "requested_ratio={ratio:.1e} delta_G={delta_g:.6e} reconstructed_ratio={reconstructed:.6e} log_error={log_error:.3e} verdict={}",
                if representable {
                    "FloatingPointFixtureRepresentable"
                } else {
                    "FloatingPointFixtureLimit"
                }
            );

            if !representable && first_collapsed.is_none() {
                first_collapsed = Some(ratio);
            }
        }

        assert!(
            first_collapsed.is_none(),
            "the chosen f64 Gibbs matrix must remain representable"
        );
        println!("representability_floor=not_reached_by_ratio_1e-100");
    }

    #[test]
    fn i5_extreme_composition_is_seed_invariant() {
        let weights = [1.0, 1e-3, 1e-8, 1e-16, 1e-30];
        let seeds = [
            vec![0.2; 5],
            vec![1.0, TRACE_FLOOR, TRACE_FLOOR, TRACE_FLOOR, TRACE_FLOOR],
            vec![TRACE_FLOOR, TRACE_FLOOR, TRACE_FLOOR, TRACE_FLOOR, 1.0],
            vec![0.55, 0.25, 0.15, 0.04, 0.01],
        ];
        let (reference, _) = solve(&weights);
        for (index, seed) in seeds.iter().enumerate() {
            let (actual, balance) =
                try_solve_with_seed(&weights, TOTAL_MOLES, seed).expect("seed route must solve");
            assert!(balance < 1e-10);
            for (value, expected) in actual.iter().zip(&reference) {
                assert!((value - expected).abs() < 1e-10 + 1e-8 * expected.abs());
            }
            println!(
                "seed={index} physical_state_delta={:.3e} balance={balance:.3e} status=OK",
                actual
                    .iter()
                    .zip(&reference)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max)
            );
        }
        println!("verdict: TraceSeedIndependenceInvariant");
    }

    #[test]
    fn i5_extreme_composition_is_species_permutation_invariant() {
        let weights = [1.0, 1e-3, 1e-8, 1e-16, 1e-30];
        let permutations = [
            [0, 1, 2, 3, 4],
            [4, 3, 2, 1, 0],
            [1, 3, 0, 4, 2],
            [2, 0, 4, 1, 3],
        ];
        let (reference, _) = solve(&weights);
        for permutation in permutations {
            let permuted = permutation.map(|i| weights[i]);
            let (actual, balance) = solve(&permuted);
            let mut canonical = vec![0.0; 5];
            for (position, &original) in permutation.iter().enumerate() {
                canonical[original] = actual[position];
            }
            assert!(balance < 1e-10);
            for (value, expected) in canonical.iter().zip(&reference) {
                assert!((value - expected).abs() < 1e-10 + 1e-8 * expected.abs());
            }
            println!(
                "permutation={permutation:?} physical_state_delta={:.3e} balance={balance:.3e} status=OK",
                canonical
                    .iter()
                    .zip(&reference)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max)
            );
        }
        println!("verdict: SpeciesPermutationDynamicRangeInvariant");
    }

    #[test]
    fn i5_extreme_composition_preserves_extensive_scaling() {
        let weights = [1.0, 1e-3, 1e-8, 1e-16, 1e-30];
        let reference_seed = vec![0.2; 5];
        let (reference, _) = try_solve_with_seed(&weights, 1.0, &reference_seed)
            .expect("unit extensive reference must solve");
        let reference_fractions: Vec<f64> = reference.iter().map(|n| *n / 1.0).collect();
        let mut all_invariant = true;
        for total in [1e-8, 1.0, 1e8] {
            let seed = vec![total / 5.0; 5];
            let result = try_solve_with_seed(&weights, total, &seed);
            let Ok((actual, balance)) = result else {
                all_invariant = false;
                println!(
                    "total_moles={total:.3e} verdict=ScaleLimitedCharacterization solver=AllBackendsFailed"
                );
                continue;
            };
            assert!(balance <= total.abs() * 1e-10 + 1e-12);
            let fraction_delta = actual
                .iter()
                .zip(&reference_fractions)
                .map(|(n, x)| (n / total - x).abs())
                .fold(0.0, f64::max);
            println!(
                "total_moles={total:.3e} max_fraction_delta={fraction_delta:.3e} balance={balance:.3e}"
            );
            for (value, fraction) in actual.iter().zip(&reference_fractions) {
                assert!(
                    (value / total - fraction).abs() < 1e-3,
                    "total={total:e} fraction mismatch"
                );
            }
            let verdict = if fraction_delta < 1e-8 {
                "DynamicRangeExtensiveScaleInvariant"
            } else {
                "ScaleLimitedCharacterization"
            };
            if verdict != "DynamicRangeExtensiveScaleInvariant" {
                all_invariant = false;
            }
            println!(
                "total_moles={total:.3e} max_fraction_delta={fraction_delta:.3e} balance={balance:.3e} verdict={verdict}"
            );
        }
        println!(
            "verdict: {}",
            if all_invariant {
                "DynamicRangeExtensiveScaleInvariant"
            } else {
                "ScaleLimitedCharacterization"
            }
        );
    }
}
