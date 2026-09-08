//! Test-only IAPWS/NIST standard-Gibbs gauge and independent ideal phase oracle.
//!
//! This module deliberately owns relative `G0` differences only. It is valid
//! for fixed `P,T` lifecycle tests and is not a thermochemistry-library record
//! or an enthalpy-capable `P,H` model.

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::R;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{GibbsFn, Phase, Solvers};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
    use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::{
        PreparedPhaseControlOutcome, PreparedPhaseControlRunner,
    };
    use nalgebra::DMatrix;

    const P0_PA: f64 = 100_000.0;
    const WATER_T_STAR_K: f64 = 273.16;
    const WATER_P_STAR_PA: f64 = 611.657;

    fn water_ice_sublimation_pa(temperature_k: f64) -> f64 {
        assert!((50.0..=273.16).contains(&temperature_k));
        let theta = temperature_k / WATER_T_STAR_K;
        let sum = -0.212_144_006e2 * theta.powf(0.333_333_333e-2)
            + 0.273_203_819e2 * theta.powf(0.120_666_667e1)
            - 0.610_598_130e1 * theta.powf(0.170_333_333e1);
        WATER_P_STAR_PA * (sum / theta).exp()
    }

    fn co2_sublimation_pa(temperature_k: f64) -> f64 {
        assert!((154.26..=195.89).contains(&temperature_k));
        100_000.0 * 10_f64.powf(6.812_28 - 1301.679 / (temperature_k - 3.494))
    }

    fn solid_gauge_gibbs(temperature_k: f64, sublimation_pressure_pa: f64) -> f64 {
        R * temperature_k * (sublimation_pressure_pa / P0_PA).ln()
    }

    #[derive(Debug, Clone, Copy, PartialEq)]
    struct OracleState {
        ice_active: bool,
        dry_ice_active: bool,
        water_gas_moles: f64,
        co2_gas_moles: f64,
        ice_moles: f64,
        dry_ice_moles: f64,
    }

    const TRACE_MOLES: f64 = 1e-30;

    /// Builds the real bounded phase-control runner. Only the standard Gibbs
    /// gauge is test-local; every active-set and nonlinear-solve operation is
    /// production machinery.
    fn gauge_gibbs() -> Vec<GibbsFn> {
        vec![
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(|temperature| {
                solid_gauge_gibbs(temperature, water_ice_sublimation_pa(temperature))
            }),
            Rc::new(|temperature| solid_gauge_gibbs(temperature, co2_sublimation_pa(temperature))),
        ]
    }

    fn runner_with_initial(
        temperature_k: f64,
        pressure_pa: f64,
        initial: Vec<f64>,
    ) -> PreparedPhaseControlRunner {
        runner_with_initial_phase_order(temperature_k, pressure_pa, initial, true)
    }

    fn runner_with_initial_phase_order(
        temperature_k: f64,
        pressure_pa: f64,
        physical_initial: Vec<f64>,
        ice_before_dry_ice: bool,
    ) -> PreparedPhaseControlRunner {
        // `EquilibriumProblem` deliberately requires every phase to occupy a
        // contiguous component range. A phase-order metamorphic test must
        // therefore permute the complete component layout, not only `Phase`s.
        let initial = if ice_before_dry_ice {
            physical_initial
        } else {
            vec![
                physical_initial[0],
                physical_initial[1],
                physical_initial[2],
                physical_initial[4],
                physical_initial[3],
            ]
        };
        let seed = initial
            .iter()
            .map(|moles| moles.max(TRACE_MOLES))
            .collect::<Vec<_>>();
        let physical_elements = [
            1.0, 0.0, 0.0, // Ar(g)
            0.0, 1.0, 0.0, // H2O(g)
            0.0, 0.0, 1.0, // CO2(g)
            0.0, 1.0, 0.0, // ice
            0.0, 0.0, 1.0, // dry ice
        ];
        let elements = DMatrix::from_row_slice(
            5,
            3,
            if ice_before_dry_ice {
                &physical_elements
            } else {
                &[
                    1.0, 0.0, 0.0, // Ar(g)
                    0.0, 1.0, 0.0, // H2O(g)
                    0.0, 0.0, 1.0, // CO2(g)
                    0.0, 0.0, 1.0, // dry ice
                    0.0, 1.0, 0.0, // ice
                ]
            },
        );
        let mut gibbs = gauge_gibbs();
        if !ice_before_dry_ice {
            gibbs.swap(3, 4);
        }
        let problem = EquilibriumProblem::new(
            if ice_before_dry_ice {
                vec!["Ar(g)", "H2O(g)", "CO2(g)", "H2O(s,ice)", "CO2(s)"]
            } else {
                vec!["Ar(g)", "H2O(g)", "CO2(g)", "CO2(s)", "H2O(s,ice)"]
            }
            .into_iter()
            .map(str::to_owned)
            .collect(),
            initial,
            LogMolesInitialGuess::from_moles(&seed, TRACE_MOLES).unwrap(),
            elements,
            gibbs,
            if ice_before_dry_ice {
                vec![
                    Phase {
                        kind: PhaseActivityModel::IdealGas,
                        species: vec![0, 1, 2],
                    },
                    Phase {
                        kind: PhaseActivityModel::IdealSolution,
                        species: vec![3],
                    },
                    Phase {
                        kind: PhaseActivityModel::IdealSolution,
                        species: vec![4],
                    },
                ]
            } else {
                vec![
                    Phase {
                        kind: PhaseActivityModel::IdealGas,
                        species: vec![0, 1, 2],
                    },
                    Phase {
                        kind: PhaseActivityModel::IdealSolution,
                        species: vec![3],
                    },
                    Phase {
                        kind: PhaseActivityModel::IdealSolution,
                        species: vec![4],
                    },
                ]
            },
            EquilibriumConditions::new(temperature_k, pressure_pa, P0_PA).unwrap(),
        )
        .unwrap();
        let mut runner = PreparedPhaseControlRunner::new(problem, Vec::new(), false).unwrap();
        let settings = runner.configure_solver();
        settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
        settings.solver_params.tol = 1e-11;
        settings.solver_params.max_iter = 250;
        runner
    }

    fn runner(
        temperature_k: f64,
        pressure_pa: f64,
        ar: f64,
        water: f64,
        co2: f64,
    ) -> PreparedPhaseControlRunner {
        runner_with_initial(temperature_k, pressure_pa, vec![ar, water, co2, 0.0, 0.0])
    }

    /// Retargets only numerical thermochemistry and carries a state that has
    /// already passed the preceding point's production acceptance contract.
    fn continue_at(
        runner: &mut PreparedPhaseControlRunner,
        previous: &PreparedPhaseControlOutcome,
        temperature_k: f64,
        pressure_pa: f64,
    ) {
        runner
            .retarget_numeric(
                EquilibriumConditions::new(temperature_k, pressure_pa, P0_PA).unwrap(),
                LogMolesInitialGuess::new(previous.solution.log_moles().to_vec()).unwrap(),
                gauge_gibbs(),
            )
            .unwrap();
        runner
            .set_continuation_phase_set(previous.phase_control_report.final_phase_set.clone())
            .unwrap();
    }

    fn assert_transition_evidence(outcome: &PreparedPhaseControlOutcome) {
        for transition in &outcome.phase_control_report.transitions {
            assert!(
                !transition.activated.is_empty() || !transition.deactivated.is_empty(),
                "accepted transition must change the active set"
            );
            assert!(
                transition
                    .minimum_tpds
                    .iter()
                    .flatten()
                    .all(|value| value.is_finite()),
                "transition contains a non-finite TPD"
            );
        }
        assert!(
            outcome
                .acceptance_report
                .phase_stability
                .iter()
                .filter_map(|report| report.minimum_tpd)
                .all(|value| value.is_finite())
        );
        assert!(outcome.acceptance_report.complementarity.satisfied);
    }

    fn assert_matches_oracle(
        outcome: &PreparedPhaseControlOutcome,
        expected: OracleState,
        ar: f64,
        _water: f64,
        _co2: f64,
    ) {
        let mask = outcome.phase_control_report.final_phase_set.active_mask();
        assert_eq!(
            mask,
            vec![true, expected.ice_active, expected.dry_ice_active]
        );
        let actual = outcome.solution.moles();
        for (actual, expected) in [
            (actual[0], ar),
            (actual[1], expected.water_gas_moles),
            (actual[2], expected.co2_gas_moles),
            (actual[3], expected.ice_moles),
            (actual[4], expected.dry_ice_moles),
        ] {
            assert!(
                (actual - expected).abs() <= 1e-12 + 2e-7 * expected.abs(),
                "actual={actual:e}, expected={expected:e}"
            );
        }
        assert!(outcome.acceptance_report.complementarity.satisfied);
    }

    fn assert_matches_oracle_with_phase_order(
        outcome: &PreparedPhaseControlOutcome,
        expected: OracleState,
        ar: f64,
        ice_before_dry_ice: bool,
    ) {
        let moles = outcome.solution.moles();
        let physical = if ice_before_dry_ice {
            moles.to_vec()
        } else {
            vec![moles[0], moles[1], moles[2], moles[4], moles[3]]
        };
        for (actual, expected) in [
            (physical[0], ar),
            (physical[1], expected.water_gas_moles),
            (physical[2], expected.co2_gas_moles),
            (physical[3], expected.ice_moles),
            (physical[4], expected.dry_ice_moles),
        ] {
            assert!(
                (actual - expected).abs() <= 1e-12 + 2e-7 * expected.abs(),
                "actual={actual:e}, expected={expected:e}"
            );
        }
        assert!(outcome.acceptance_report.complementarity.satisfied);
        let mask = outcome.phase_control_report.final_phase_set.active_mask();
        let expected_mask = if ice_before_dry_ice {
            vec![true, expected.ice_active, expected.dry_ice_active]
        } else {
            vec![true, expected.dry_ice_active, expected.ice_active]
        };
        assert_eq!(mask, expected_mask);
    }

    /// Enumerates the four physical pure-condensed masks algebraically. The
    /// returned state is the only mask satisfying positive condensed amounts
    /// and stable inactive partial pressures.
    fn oracle(temperature_k: f64, pressure_pa: f64, ar: f64, water: f64, co2: f64) -> OracleState {
        let pw = water_ice_sublimation_pa(temperature_k);
        let pc = co2_sublimation_pa(temperature_k);
        for ice_active in [false, true] {
            for dry_ice_active in [false, true] {
                let fixed = ar
                    + if ice_active { 0.0 } else { water }
                    + if dry_ice_active { 0.0 } else { co2 };
                let fractions = (if ice_active { pw / pressure_pa } else { 0.0 })
                    + if dry_ice_active {
                        pc / pressure_pa
                    } else {
                        0.0
                    };
                if fractions >= 1.0 {
                    continue;
                }
                let gas_total = fixed / (1.0 - fractions);
                let water_gas = if ice_active {
                    gas_total * pw / pressure_pa
                } else {
                    water
                };
                let co2_gas = if dry_ice_active {
                    gas_total * pc / pressure_pa
                } else {
                    co2
                };
                let ice = water - water_gas;
                let dry_ice = co2 - co2_gas;
                let water_partial = pressure_pa * water_gas / gas_total;
                let co2_partial = pressure_pa * co2_gas / gas_total;
                let valid = (!ice_active && water_partial <= pw) || (ice_active && ice > 0.0);
                let valid = valid
                    && ((!dry_ice_active && co2_partial <= pc)
                        || (dry_ice_active && dry_ice > 0.0));
                if valid {
                    return OracleState {
                        ice_active,
                        dry_ice_active,
                        water_gas_moles: water_gas,
                        co2_gas_moles: co2_gas,
                        ice_moles: ice.max(0.0),
                        dry_ice_moles: dry_ice.max(0.0),
                    };
                }
            }
        }
        panic!("no physical gas/ice/dry-ice oracle state at {temperature_k} K");
    }

    #[test]
    fn sublimation_gauge_reproduces_each_external_boundary_by_algebra() {
        for temperature_k in [160.0, 180.0, 194.0] {
            for pressure in [
                water_ice_sublimation_pa(temperature_k),
                co2_sublimation_pa(temperature_k),
            ] {
                let g0 = solid_gauge_gibbs(temperature_k, pressure);
                let recovered = P0_PA * (g0 / (R * temperature_k)).exp();
                assert!(((recovered - pressure) / pressure).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn scalar_oracle_has_a_non_degenerate_gas_one_solid_three_phase_interval() {
        let args = (100_000.0, 1.0, 1.0e-7, 1.0);
        let gas = oracle(194.0, args.0, args.1, args.2, args.3);
        let one_solid = oracle(185.0, args.0, args.1, args.2, args.3);
        let three_phase = oracle(170.0, args.0, args.1, args.2, args.3);
        assert!(!gas.ice_active && !gas.dry_ice_active);
        assert!(!one_solid.ice_active && one_solid.dry_ice_active);
        assert!(three_phase.ice_active && three_phase.dry_ice_active);
        assert!(three_phase.ice_moles / args.2 > 1e-3);
        assert!(three_phase.dry_ice_moles / args.3 > 1e-3);
        assert!(three_phase.water_gas_moles > 0.0 && three_phase.co2_gas_moles > 0.0);
    }

    #[test]
    fn production_runner_matches_gauge_oracle_at_each_topology_interior() {
        let (pressure, ar, water, co2) = (100_000.0, 1.0, 1e-7, 1.0);
        for temperature in [194.0, 185.0, 170.0] {
            let expected = oracle(temperature, pressure, ar, water, co2);
            let outcome = runner(temperature, pressure, ar, water, co2)
                .solve()
                .expect("production runner must solve the gauge topology interior");
            assert_matches_oracle(&outcome, expected, ar, water, co2);
        }
    }

    #[test]
    fn production_runner_continues_the_full_gas_dry_ice_ice_lifecycle_both_directions() {
        let (pressure, ar, water, co2) = (100_000.0, 1.0, 1e-7, 1.0);
        let cooling = [194.0, 185.0, 170.0];
        let heating = [170.0, 185.0, 194.0];

        // Fresh gas-only cooling: gas -> gas + dry ice -> gas + dry ice + ice.
        let mut forward = runner(cooling[0], pressure, ar, water, co2);
        let mut forward_outcome = forward.solve().unwrap();
        assert_transition_evidence(&forward_outcome);
        assert_matches_oracle(
            &forward_outcome,
            oracle(cooling[0], pressure, ar, water, co2),
            ar,
            water,
            co2,
        );
        let mut forward_transitions = forward_outcome.phase_control_report.transitions.len();
        for temperature in cooling.into_iter().skip(1) {
            continue_at(&mut forward, &forward_outcome, temperature, pressure);
            forward_outcome = forward.solve().unwrap();
            assert_transition_evidence(&forward_outcome);
            assert_matches_oracle(
                &forward_outcome,
                oracle(temperature, pressure, ar, water, co2),
                ar,
                water,
                co2,
            );
            forward_transitions += forward_outcome.phase_control_report.transitions.len();
        }

        // Reverse heating starts from the physical three-phase state, then
        // must deactivate ice before dry ice without inheriting a gas-only seed.
        let initial = oracle(heating[0], pressure, ar, water, co2);
        let mut reverse = runner_with_initial(
            heating[0],
            pressure,
            vec![
                ar,
                initial.water_gas_moles,
                initial.co2_gas_moles,
                initial.ice_moles,
                initial.dry_ice_moles,
            ],
        );
        let mut reverse_outcome = reverse.solve().unwrap();
        assert_transition_evidence(&reverse_outcome);
        assert_matches_oracle(&reverse_outcome, initial, ar, water, co2);
        let mut reverse_transitions = reverse_outcome.phase_control_report.transitions.len();
        for temperature in heating.into_iter().skip(1) {
            continue_at(&mut reverse, &reverse_outcome, temperature, pressure);
            reverse_outcome = reverse.solve().unwrap();
            assert_transition_evidence(&reverse_outcome);
            assert_matches_oracle(
                &reverse_outcome,
                oracle(temperature, pressure, ar, water, co2),
                ar,
                water,
                co2,
            );
            reverse_transitions += reverse_outcome.phase_control_report.transitions.len();
        }

        assert_eq!(
            forward_transitions, 2,
            "cooling must activate both solids once"
        );
        assert_eq!(
            reverse_transitions, 2,
            "heating must deactivate both solids once"
        );
    }

    #[test]
    fn three_phase_fixed_point_is_invariant_to_condensed_phase_declaration_order() {
        let (temperature, pressure, ar, water, co2) = (170.0, 100_000.0, 1.0, 1e-7, 1.0);
        let expected = oracle(temperature, pressure, ar, water, co2);
        let mut states = Vec::new();

        for ice_before_dry_ice in [true, false] {
            let outcome = runner_with_initial_phase_order(
                temperature,
                pressure,
                vec![ar, water, co2, 0.0, 0.0],
                ice_before_dry_ice,
            )
            .solve()
            .unwrap();
            assert_matches_oracle_with_phase_order(&outcome, expected, ar, ice_before_dry_ice);
            let moles = outcome.solution.moles();
            states.push(if ice_before_dry_ice {
                moles.to_vec()
            } else {
                vec![moles[0], moles[1], moles[2], moles[4], moles[3]]
            });
        }

        for (left, right) in states[0].iter().zip(&states[1]) {
            assert!(
                (left - right).abs() <= 1e-12 + 2e-7 * left.abs(),
                "physical state changed under phase declaration permutation: {left:e} vs {right:e}"
            );
        }
    }
}
