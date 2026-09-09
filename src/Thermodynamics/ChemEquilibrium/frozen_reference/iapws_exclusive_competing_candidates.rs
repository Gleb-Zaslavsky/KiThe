//! Test-only Ar/H2O liquid-versus-ice candidate competition.
//!
//! The local NASA records have no useful interval below the triple point, so
//! this story uses independent IAPWS boundary equations only to construct a
//! relative standard-Gibbs gauge.  The active-set lifecycle, conservation,
//! acceptance, and phase-order comparison remain production code.

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use nalgebra::DMatrix;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
        GibbsFn, Phase, R, Solvers,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
    use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::{
        PreparedPhaseControlOutcome, PreparedPhaseControlRunner,
    };

    const P0_PA: f64 = 100_000.0;
    const T_K: f64 = 273.15;
    const P_TOTAL_PA: f64 = 700.0;
    const AR_MOLES: f64 = 0.1;
    const WATER_MOLES: f64 = 1.0;
    const TRACE_MOLES: f64 = 1e-30;

    // IAPWS-IF97 Region 4 saturation pressure, valid at the liquid side.
    fn liquid_boundary_pa(t: f64) -> f64 {
        const TC: f64 = 647.096;
        const N: [f64; 6] = [
            -7.859_517_83,
            1.844_082_59,
            -11.786_649_7,
            22.680_741_1,
            -15.961_871_9,
            1.801_225_02,
        ];
        const E: [f64; 6] = [1.0, 1.5, 3.0, 3.5, 4.0, 7.5];
        let theta = 1.0 - t / TC;
        let sum: f64 = N.iter().zip(E).map(|(a, p)| a * theta.powf(p)).sum();
        22.064e6 * (TC / t * sum).exp()
    }

    // The official IAPWS ice-Ih sublimation equation used by the existing
    // frozen sublimation evidence. Keeping it local makes this gauge
    // independent from the production thermochemistry records.
    fn ice_boundary_pa(t: f64) -> f64 {
        const TS: f64 = 273.16;
        const PS: f64 = 611.657;
        let theta = t / TS;
        let sum = -0.212_144_006e2 * theta.powf(0.333_333_333e-2)
            + 0.273_203_819e2 * theta.powf(0.120_666_667e1)
            - 0.610_598_130e1 * theta.powf(0.170_333_333e1);
        PS * (sum / theta).exp()
    }

    fn runner_with_initial(
        liquid_before_ice: bool,
        physical_initial: [f64; 4],
    ) -> PreparedPhaseControlRunner {
        let physical_names = ["Ar(g)", "H2O(g)", "H2O(l)", "H2O(s,ice Ih)"];
        let order = if liquid_before_ice {
            [0, 1, 2, 3]
        } else {
            [0, 1, 3, 2]
        };
        let names = order
            .iter()
            .map(|&i| physical_names[i].to_owned())
            .collect();
        let initial = order
            .iter()
            .map(|&i| physical_initial[i])
            .collect::<Vec<_>>();
        let seed = initial
            .iter()
            .map(|n| n.max(TRACE_MOLES))
            .collect::<Vec<_>>();
        let elements = DMatrix::from_row_slice(
            4,
            2,
            &order
                .iter()
                .flat_map(|&i| match i {
                    0 => [1.0, 0.0],
                    1 | 2 | 3 => [0.0, 1.0],
                    _ => unreachable!(),
                })
                .collect::<Vec<_>>(),
        );
        let gibbs = order
            .iter()
            .map(|&i| match i {
                0 | 1 => Rc::new(|_| 0.0) as GibbsFn,
                2 => Rc::new(|_| R * T_K * (liquid_boundary_pa(T_K) / P0_PA).ln()) as GibbsFn,
                3 => Rc::new(|_| R * T_K * (ice_boundary_pa(T_K) / P0_PA).ln()) as GibbsFn,
                _ => unreachable!(),
            })
            .collect();
        let problem = EquilibriumProblem::new(
            names,
            initial,
            LogMolesInitialGuess::from_moles(&seed, TRACE_MOLES).unwrap(),
            elements,
            gibbs,
            vec![
                Phase {
                    kind: PhaseActivityModel::IdealGas,
                    species: vec![0, 1],
                },
                Phase {
                    kind: PhaseActivityModel::IdealSolution,
                    species: vec![2],
                },
                Phase {
                    kind: PhaseActivityModel::IdealSolution,
                    species: vec![3],
                },
            ],
            EquilibriumConditions::new(T_K, P_TOTAL_PA, P0_PA).unwrap(),
        )
        .unwrap();
        let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false).unwrap();
        let settings = runner.configure_solver();
        settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
        settings.solver_params.tol = 1e-11;
        settings.solver_params.max_iter = 250;
        runner
    }

    fn runner(liquid_before_ice: bool) -> PreparedPhaseControlRunner {
        runner_with_initial(liquid_before_ice, [AR_MOLES, WATER_MOLES, 0.0, 0.0])
    }

    fn physical_moles(outcome: &PreparedPhaseControlOutcome, liquid_before_ice: bool) -> [f64; 4] {
        let m = outcome.solution.moles();
        if liquid_before_ice {
            [m[0], m[1], m[2], m[3]]
        } else {
            [m[0], m[1], m[3], m[2]]
        }
    }

    fn controlled_runner(delta_g: f64, a_before_b: bool) -> PreparedPhaseControlRunner {
        let base = R * T_K * (ice_boundary_pa(T_K) / P0_PA).ln();
        let order = if a_before_b {
            [0, 1, 2, 3]
        } else {
            [0, 1, 3, 2]
        };
        let names = order
            .iter()
            .map(|&i| ["Ar(g)", "H2O(g)", "condensed_a", "condensed_b"][i].to_owned())
            .collect();
        let physical_initial = [AR_MOLES, WATER_MOLES, 0.0, 0.0];
        let initial = order
            .iter()
            .map(|&i| physical_initial[i])
            .collect::<Vec<_>>();
        let seed = initial
            .iter()
            .map(|n| n.max(TRACE_MOLES))
            .collect::<Vec<_>>();
        let elements = DMatrix::from_row_slice(
            4,
            2,
            &order
                .iter()
                .flat_map(|&i| match i {
                    0 => [1.0, 0.0],
                    _ => [0.0, 1.0],
                })
                .collect::<Vec<_>>(),
        );
        let gibbs = order
            .iter()
            .map(|&i| match i {
                0 | 1 => Rc::new(|_| 0.0) as GibbsFn,
                2 => Rc::new(move |_| base) as GibbsFn,
                3 => Rc::new(move |_| base + delta_g) as GibbsFn,
                _ => unreachable!(),
            })
            .collect();
        let problem = EquilibriumProblem::new(
            names,
            initial,
            LogMolesInitialGuess::from_moles(&seed, TRACE_MOLES).unwrap(),
            elements,
            gibbs,
            vec![
                Phase {
                    kind: PhaseActivityModel::IdealGas,
                    species: vec![0, 1],
                },
                Phase {
                    kind: PhaseActivityModel::IdealSolution,
                    species: vec![2],
                },
                Phase {
                    kind: PhaseActivityModel::IdealSolution,
                    species: vec![3],
                },
            ],
            EquilibriumConditions::new(T_K, P_TOTAL_PA, P0_PA).unwrap(),
        )
        .unwrap();
        let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false).unwrap();
        let settings = runner.configure_solver();
        settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
        settings.solver_params.tol = 1e-11;
        settings.solver_params.max_iter = 250;
        runner
    }

    #[test]
    fn i5_iapws_water_has_two_unstable_condensed_candidates_and_order_invariant_fixed_point() {
        let p_liquid = liquid_boundary_pa(T_K);
        let p_ice = ice_boundary_pa(T_K);
        let p_water_initial = P_TOTAL_PA * WATER_MOLES / (AR_MOLES + WATER_MOLES);
        let tpd_liquid = R * T_K * (p_liquid / p_water_initial).ln();
        let tpd_ice = R * T_K * (p_ice / p_water_initial).ln();
        let dg_create = 1e-6;

        println!(
            "water candidate preflight: T={T_K} p_liquid={p_liquid:.12} p_ice={p_ice:.12} delta_p={:.12}",
            p_liquid - p_ice
        );
        assert!(
            p_ice < p_liquid,
            "IAPWS must place ice below liquid below triple point"
        );
        assert!(tpd_liquid < dg_create && tpd_ice < dg_create);
        assert!(
            (tpd_liquid - tpd_ice).abs() > 1e-4,
            "candidate ordering is below diagnostic resolution"
        );

        let mut a = runner(true);
        let mut b = runner(false);
        let out_a = a.solve().unwrap();
        let out_b = b.solve().unwrap();
        let ma = physical_moles(&out_a, true);
        let mb = physical_moles(&out_b, false);
        assert!(out_a.acceptance_report.complementarity.satisfied);
        assert!(out_b.acceptance_report.complementarity.satisfied);
        let initial_tpds_a = out_a
            .phase_control_report
            .transitions
            .first()
            .expect("candidate competition must produce an accepted transition")
            .minimum_tpds
            .clone();
        assert!(
            initial_tpds_a
                .get(1)
                .and_then(|v| *v)
                .is_some_and(|v| v < dg_create)
        );
        assert!(
            initial_tpds_a
                .get(2)
                .and_then(|v| *v)
                .is_some_and(|v| v < dg_create)
        );
        assert_eq!(
            out_a.phase_control_report.final_phase_set.active_mask(),
            vec![true, false, true]
        );
        assert_eq!(
            out_b.phase_control_report.final_phase_set.active_mask(),
            vec![true, true, false]
        );
        for (x, y) in ma.into_iter().zip(mb) {
            assert!((x - y).abs() < 1e-8);
        }
        assert!(ma[3] > 0.0 && ma[2] <= 1e-10);
        let expected_liquid_tpd = R * T_K * (p_liquid / p_ice).ln();
        let final_liquid_tpd = out_a.acceptance_report.phase_stability[1]
            .minimum_tpd
            .expect("final inactive liquid must have a canonical TPD");
        assert!(final_liquid_tpd > 0.0);
        assert!((final_liquid_tpd - expected_liquid_tpd).abs() < 1e-6);
        println!(
            "Ar + H2O exclusive candidates: T={T_K:.2} K P={P_TOTAL_PA:.1} Pa p_liquid={p_liquid:.6} Pa p_ice={p_ice:.6} TPD_liquid={tpd_liquid:.6e} TPD_ice={tpd_ice:.6e} final=gas+ice order_invariant=true"
        );
    }

    #[test]
    fn i5_iapws_water_metastable_replacement_and_stable_ice_retention_matrix() {
        let p_liquid = liquid_boundary_pa(T_K);
        let p_ice = ice_boundary_pa(T_K);
        let gas_liquid_water = AR_MOLES * p_liquid / (P_TOTAL_PA - p_liquid);
        let gas_ice_water = AR_MOLES * p_ice / (P_TOTAL_PA - p_ice);
        let histories = [
            ("gas", [AR_MOLES, WATER_MOLES, 0.0, 0.0]),
            (
                "gas+liquid",
                [
                    AR_MOLES,
                    gas_liquid_water,
                    WATER_MOLES - gas_liquid_water,
                    0.0,
                ],
            ),
            (
                "gas+ice",
                [AR_MOLES, gas_ice_water, 0.0, WATER_MOLES - gas_ice_water],
            ),
        ];
        let mut reference: Option<[f64; 4]> = None;
        let mut route_count = 0;
        for (history, initial) in histories {
            for liquid_before_ice in [true, false] {
                let mut runner = runner_with_initial(liquid_before_ice, initial);
                let outcome = runner.solve().unwrap_or_else(|error| {
                    panic!("history={history} liquid_before_ice={liquid_before_ice}: {error:?}")
                });
                let physical = physical_moles(&outcome, liquid_before_ice);
                assert_eq!(
                    outcome.phase_control_report.final_phase_set.active_mask(),
                    if liquid_before_ice {
                        vec![true, false, true]
                    } else {
                        vec![true, true, false]
                    }
                );
                assert!(outcome.acceptance_report.complementarity.satisfied);
                assert!(physical[3] > 0.0 && physical[2] <= 1e-10);
                if let Some(expected) = reference {
                    for (actual, expected) in physical.into_iter().zip(expected) {
                        assert!(
                            (actual - expected).abs() < 1e-8,
                            "history={history} order={liquid_before_ice} physical state differs: {actual:e} vs {expected:e}"
                        );
                    }
                } else {
                    reference = Some(physical);
                }
                route_count += 1;
                println!(
                    "history={history:<9} order={} final=gas+ice transitions={} liquid_moles={:.6e} ice_moles={:.6e}",
                    if liquid_before_ice { "L->I" } else { "I->L" },
                    outcome.phase_control_report.transitions.len(),
                    physical[2],
                    physical[3]
                );
            }
        }
        assert_eq!(route_count, 6);
        println!("verdict: ExclusiveCompetitionHistoryInvariant");
    }

    #[test]
    fn i5_controlled_condensed_gibbs_split_matrix_has_analytic_losing_tpd() {
        let deltas = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0.0];
        let initial_tpd_a = R
            * T_K
            * (ice_boundary_pa(T_K) / (P_TOTAL_PA * WATER_MOLES / (AR_MOLES + WATER_MOLES))).ln();
        assert!(initial_tpd_a < 0.0);
        let mut resolved = 0;
        for delta in deltas {
            let mut a_first = controlled_runner(delta, true);
            let mut b_first = controlled_runner(delta, false);
            let outcome_a = a_first
                .solve()
                .unwrap_or_else(|error| panic!("delta={delta:e} A->B: {error:?}"));
            let outcome_b = b_first
                .solve()
                .unwrap_or_else(|error| panic!("delta={delta:e} B->A: {error:?}"));
            assert!(outcome_a.acceptance_report.complementarity.satisfied);
            assert!(outcome_b.acceptance_report.complementarity.satisfied);
            let ma = outcome_a.solution.moles();
            let mb = outcome_b.solution.moles();
            let physical_b = [mb[0], mb[1], mb[3], mb[2]];
            assert!((ma[0] - physical_b[0]).abs() < 1e-8);
            assert!((ma[1] - physical_b[1]).abs() < 1e-8);
            assert!((ma[2] + ma[3] - physical_b[2] - physical_b[3]).abs() < 1e-8);
            if delta > 0.0 {
                assert_eq!(
                    outcome_a.phase_control_report.final_phase_set.active_mask(),
                    vec![true, true, false]
                );
                assert_eq!(
                    outcome_b.phase_control_report.final_phase_set.active_mask(),
                    vec![true, false, true]
                );
                let losing = outcome_a.acceptance_report.phase_stability[2]
                    .minimum_tpd
                    .expect("losing phase B must have TPD evidence");
                assert!(losing > 0.0);
                assert!((losing - delta).abs() < 1e-6);
                resolved += 1;
            } else {
                assert!(ma[2] + ma[3] > 0.0);
            }
            println!(
                "delta_G={delta:.1e} A/B physical_state_invariant=true transitions={}/{}",
                outcome_a.phase_control_report.transitions.len(),
                outcome_b.phase_control_report.transitions.len()
            );
        }
        assert_eq!(resolved, 7);
        println!(
            "verdict: ResolvedCompetitionInvariant; exact delta_G=0 classified as non-unique representative"
        );
    }

    #[test]
    fn i5_controlled_condensed_gibbs_split_measures_observed_tpd_scatter() {
        let delta_g = 1e-6;
        let mut samples = Vec::new();
        for a_before_b in [true, false, true, false, true, false] {
            let mut runner = controlled_runner(delta_g, a_before_b);
            let outcome = runner
                .solve()
                .expect("well-resolved controlled split must solve");
            assert!(outcome.acceptance_report.complementarity.satisfied);
            let loser_index = if a_before_b { 2 } else { 1 };
            let tpd = outcome.acceptance_report.phase_stability[loser_index]
                .minimum_tpd
                .expect("inactive losing candidate must expose TPD evidence");
            assert!(tpd.is_finite() && tpd > 0.0);
            samples.push(tpd);
        }
        let min = samples.iter().copied().fold(f64::INFINITY, f64::min);
        let max = samples.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let scatter = max - min;
        assert!(scatter.is_finite() && scatter >= 0.0);
        assert!(
            (min - delta_g).abs() <= 3.0e-12 && (max - delta_g).abs() <= 3.0e-12,
            "controlled TPD must retain the known Gibbs split"
        );
        assert!(
            scatter <= 1.0e-12,
            "permuted/repeated TPD evaluation unexpectedly scattered: {scatter:e}"
        );
        println!(
            "controlled Gibbs TPD scatter: delta_G={delta_g:.3e} samples={} min={min:.12e} max={max:.12e} scatter={scatter:.3e}",
            samples.len()
        );
    }

    #[test]
    fn i5_controlled_condensed_gibbs_split_resolution_matrix() {
        let deltas = [
            1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 0.0,
        ];
        let base = R * T_K * (ice_boundary_pa(T_K) / P0_PA).ln();
        let mut smallest_representable = None;
        let mut smallest_quantitative = None;
        let mut smallest_unique_winner = None;
        let mut first_tolerance_limited = None;
        let mut floating_point_collapse = None;

        println!("controlled Gibbs resolution matrix");
        println!(
            "requested_delta actual_delta tpd_AB tpd_BA rel_err_AB rel_err_BA winner_AB winner_BA state_delta transitions verdict"
        );
        for requested in deltas {
            let actual = (base + requested) - base;
            if actual != 0.0 {
                smallest_representable = Some(requested);
            }
            if actual == 0.0 && requested > 0.0 && floating_point_collapse.is_none() {
                floating_point_collapse = Some(requested);
            }
            let mut ab_runner = controlled_runner(requested, true);
            let mut ba_runner = controlled_runner(requested, false);
            let ab = ab_runner
                .solve()
                .expect("controlled AB matrix row must terminate");
            let ba = ba_runner
                .solve()
                .expect("controlled BA matrix row must terminate");
            assert!(ab.acceptance_report.complementarity.satisfied);
            assert!(ba.acceptance_report.complementarity.satisfied);

            let ab_moles = ab.solution.moles();
            let ba_moles = ba.solution.moles();
            let ba_physical = [ba_moles[0], ba_moles[1], ba_moles[3], ba_moles[2]];
            let state_delta = ab_moles
                .iter()
                .zip(ba_physical)
                .map(|(x, y)| (x - y).abs())
                .fold(0.0, f64::max);
            let winner_ab = ab.phase_control_report.final_phase_set.active_mask()[1];
            let winner_ba = ba.phase_control_report.final_phase_set.active_mask()[2];
            let tpd_ab = ab.acceptance_report.phase_stability[2].minimum_tpd;
            let tpd_ba = ba.acceptance_report.phase_stability[1].minimum_tpd;
            let rel_ab = tpd_ab.map(|tpd| (tpd - actual).abs() / requested.max(f64::MIN_POSITIVE));
            let rel_ba = tpd_ba.map(|tpd| (tpd - actual).abs() / requested.max(f64::MIN_POSITIVE));
            let both_positive = tpd_ab.is_some_and(|v| v > 0.0) && tpd_ba.is_some_and(|v| v > 0.0);
            let quantitatively_tracks = actual != 0.0
                && tpd_ab.is_some_and(|v| (v - actual).abs() <= actual.abs() * 1e-4 + 1e-14)
                && tpd_ba.is_some_and(|v| (v - actual).abs() <= actual.abs() * 1e-4 + 1e-14);
            let unique = requested > 0.0
                && winner_ab
                && !ab.phase_control_report.final_phase_set.active_mask()[2]
                && winner_ba
                && !ba.phase_control_report.final_phase_set.active_mask()[1];
            let verdict = if requested == 0.0 {
                "ExactDegeneracyNonUnique"
            } else if actual == 0.0 {
                "FloatingPointCollapsedDegeneracy"
            } else if quantitatively_tracks && unique {
                smallest_quantitative = Some(requested);
                smallest_unique_winner = Some(requested);
                "ResolvedCompetitionInvariant"
            } else {
                if first_tolerance_limited.is_none() {
                    first_tolerance_limited = Some(requested);
                }
                "ToleranceLimitedDegeneracy"
            };
            println!(
                "{requested:.3e} {actual:.12e} {:?} {:?} {:?} {:?} {winner_ab} {winner_ba} {state_delta:.3e} {}/{} {verdict}",
                tpd_ab,
                tpd_ba,
                rel_ab,
                rel_ba,
                ab.phase_control_report.transitions.len(),
                ba.phase_control_report.transitions.len()
            );
            assert!(state_delta.is_finite());
            if requested == 0.0 {
                assert!((ab_moles[0] - ba_physical[0]).abs() < 1e-8);
                assert!((ab_moles[1] - ba_physical[1]).abs() < 1e-8);
                assert!(
                    ((ab_moles[2] + ab_moles[3]) - (ba_physical[2] + ba_physical[3])).abs() < 1e-8
                );
                assert!(ab.phase_control_report.transitions.len() <= 2);
                assert!(ba.phase_control_report.transitions.len() <= 2);
            }
            assert!(
                both_positive
                    || requested == 0.0
                    || actual == 0.0
                    || verdict == "ToleranceLimitedDegeneracy"
            );
        }
        // These values are a fixture-level f64 characterization, not phase
        // creation/retention thresholds. They make a loss of numerical
        // resolution visible before it can be mistaken for physics.
        assert_eq!(smallest_representable, Some(1e-12));
        assert_eq!(smallest_quantitative, Some(1e-7));
        assert_eq!(smallest_unique_winner, Some(1e-7));
        assert_eq!(first_tolerance_limited, Some(1e-8));
        assert_eq!(floating_point_collapse, Some(1e-13));
        println!(
            "resolution summary: smallest_representable={smallest_representable:?} smallest_quantitative={smallest_quantitative:?} smallest_unique_winner={smallest_unique_winner:?} first_tolerance_limited={first_tolerance_limited:?} floating_point_collapse={floating_point_collapse:?} exact=ExactDegeneracyNonUnique"
        );
    }
}
