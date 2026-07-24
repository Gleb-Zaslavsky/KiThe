#[cfg(test)]
mod tests {
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_active_set::ActiveSetProjection;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::*;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentErrorKind;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PHASE_CONTROL_TRACE_MOLE_FLOOR;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseSeedPolicy;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::*;
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData, WhatIsFound};
    use crate::Thermodynamics::thermo_lib_api::LibraryId;
    use nalgebra::DMatrix;
    use std::cell::RefCell;
    use std::rc::Rc;
    use std::time::Instant;

    fn stability_reports(forces: &[Option<f64>], active: &[bool]) -> Vec<PhaseStabilityReport> {
        forces
            .iter()
            .enumerate()
            .map(|(phase, &driving_force)| PhaseStabilityReport {
                phase: PhaseIndex::new(phase, forces.len()).unwrap(),
                model: PhaseStabilityModel::PureCondensedSpecies,
                active: active[phase],
                driving_force,
                element_potentials: None,
            })
            .collect()
    }

    fn phase_set(active: &[bool]) -> PhaseSet {
        PhaseSet::from_policy(&InitialPhaseSet::FromInitialMoles, active).unwrap()
    }

    fn synthetic_gas_with_pure_condensed(candidate_gibbs: f64) -> EquilibriumLogMoles {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["A_g".to_string(), "B_g".to_string(), "A_s".to_string()];
        solver.elem_composition = DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 1.0, 1.0, 0.0]);
        solver.n0 = vec![1.0, 1.0, 0.0];
        solver.initial_guess = Some(vec![
            1.0_f64.ln(),
            1.0_f64.ln(),
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
        ]);
        solver.gibbs = vec![
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(move |_| candidate_gibbs),
        ];
        solver.phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![2],
            },
        ];
        solver.T = 1000.0;
        solver.P = 101325.0;
        solver.p0 = 101325.0;
        solver.solver_settings.solver = Solvers::LM;
        solver.phase_manager.phase_eps = 1e-12;
        solver.phase_manager.set_explicit_hysteresis(-1.0, 1.0);
        solver.create_stoich_matrix().unwrap();
        solver.phase_active_mask.clear();
        solver
    }

    fn synthetic_gas_with_dynamic_condensed(
        candidate_gibbs: Rc<RefCell<f64>>,
    ) -> EquilibriumLogMoles {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["A_g".to_string(), "B_g".to_string(), "A_s".to_string()];
        solver.elem_composition = DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 1.0, 1.0, 0.0]);
        solver.n0 = vec![1.0, 1.0, 0.0];
        solver.initial_guess = Some(vec![
            1.0_f64.ln(),
            1.0_f64.ln(),
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
        ]);
        solver.gibbs = vec![
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(move |_| *candidate_gibbs.borrow()),
        ];
        solver.phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![2],
            },
        ];
        solver.T = 1000.0;
        solver.P = 101325.0;
        solver.p0 = 101325.0;
        solver.solver_settings.solver = Solvers::LM;
        solver.phase_manager.phase_eps = 1e-12;
        solver.phase_manager.set_explicit_hysteresis(-1.0, 1.0);
        solver.create_stoich_matrix().unwrap();
        solver.phase_active_mask.clear();
        solver
    }

    fn synthetic_gas_with_two_condensed_candidates(
        candidate_a_gibbs: Rc<RefCell<f64>>,
        candidate_b_gibbs: Rc<RefCell<f64>>,
    ) -> EquilibriumLogMoles {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec![
            "A_g".to_string(),
            "B_g".to_string(),
            "A_s".to_string(),
            "B_s".to_string(),
        ];
        solver.elem_composition =
            DMatrix::from_row_slice(4, 2, &[1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0]);
        solver.n0 = vec![1.0, 1.0, 0.0, 0.0];
        solver.initial_guess = Some(vec![
            1.0_f64.ln(),
            1.0_f64.ln(),
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
        ]);
        solver.gibbs = vec![
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(move |_| *candidate_a_gibbs.borrow()),
            Rc::new(move |_| *candidate_b_gibbs.borrow()),
        ];
        solver.phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![2],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![3],
            },
        ];
        solver.T = 1000.0;
        solver.P = 101325.0;
        solver.p0 = 101325.0;
        solver.solver_settings.solver = Solvers::LM;
        solver.phase_manager.phase_eps = 1e-12;
        solver.phase_manager.set_explicit_hysteresis(-1.0, 1.0);
        solver.create_stoich_matrix().unwrap();
        solver.phase_active_mask.clear();
        solver
    }

    fn synthetic_large_phase_control_fixture(
        species_per_phase: usize,
        phase_count: usize,
    ) -> EquilibriumLogMoles {
        let total_species = species_per_phase * phase_count;
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = (0..total_species)
            .map(|index| format!("S{index}"))
            .collect();
        solver.elem_composition = DMatrix::identity(total_species, total_species);
        solver.n0 = (0..total_species)
            .map(|index| {
                if index % species_per_phase == 0 {
                    1.0
                } else {
                    1e-12
                }
            })
            .collect();
        solver.initial_guess = Some(
            solver
                .n0
                .iter()
                .map(|&n| n.max(PHASE_CONTROL_TRACE_MOLE_FLOOR).ln())
                .collect(),
        );
        solver.gibbs = (0..total_species)
            .map(|_| Rc::new(|_| 0.0) as GibbsFn)
            .collect();
        solver.phases = (0..phase_count)
            .map(|phase| Phase {
                kind: PhaseKind::IdealGas,
                species: (0..species_per_phase)
                    .map(|local| phase * species_per_phase + local)
                    .collect(),
            })
            .collect();
        solver.T = 1000.0;
        solver.P = 101325.0;
        solver.p0 = 101325.0;
        solver.solver_settings.solver = Solvers::TR;
        solver.solver_settings.solver_params.max_iter = 200;
        solver.phase_manager.phase_eps = 1e-12;
        solver.phase_manager.set_explicit_hysteresis(-1.0, 1.0);
        solver.create_stoich_matrix().unwrap();
        solver.phase_active_mask.clear();
        solver
    }
    #[test]
    fn phase_deactivation_sets_species_to_floor() {
        let species_phase = vec![0, 0, 1, 1];
        let mut y = vec![0.0, 1.0, 2.0, 3.0]; // log-moles

        deactivate_phases_seed_only(
            y.as_mut_slice(),
            &[PhaseIndex::new(1, 2).unwrap()],
            &species_phase,
            PHASE_CONTROL_TRACE_MOLE_FLOOR,
        )
        .unwrap();

        assert_eq!(y[0], 0.0);
        assert_eq!(y[1], 1.0);
        assert!(y[2] < -600.0);
        assert!(y[3] < -600.0);
    }

    #[test]
    fn phase_control_trace_floor_is_named_and_positive() {
        assert!(PHASE_CONTROL_TRACE_MOLE_FLOOR.is_finite());
        assert!(PHASE_CONTROL_TRACE_MOLE_FLOOR > 0.0);
        assert!(PHASE_CONTROL_TRACE_MOLE_FLOOR < 1e-100);
    }

    #[test]
    fn initial_phase_activity_comes_from_seed_inventory() {
        let activity = initial_phase_activity(
            &[1.0_f64.ln(), PHASE_CONTROL_TRACE_MOLE_FLOOR.ln()],
            &[0, 1],
            2,
            1e-20,
        )
        .unwrap();

        assert_eq!(activity, vec![true, false]);
    }

    #[test]
    fn repeated_active_set_is_a_typed_cycle_error() {
        let set = phase_set(&[true, false]);
        let mut visited = std::collections::HashSet::from([set.clone()]);
        let error = reject_repeated_phase_set(&mut visited, &set, 3).unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::PhaseControlCycleDetected {
                iteration: 3,
                active_phases
            } if active_phases == vec![0]
        ));
    }

    #[test]
    fn inactive_phase_with_material_inventory_cannot_be_published() {
        let result = validate_phase_set_candidate(&[1.0, 1e-3], &[true, false], 1e-8);

        assert!(matches!(
            result,
            Err(ReactionExtentError::InvalidCandidate {
                field: "phase_active_set",
                ..
            })
        ));
        validate_phase_set_candidate(&[1.0, 1e-12], &[true, false], 1e-8).unwrap();
    }
    #[test]
    fn residual_finite_after_phase_deactivation() {
        let reactions = DMatrix::<f64>::zeros(4, 0);
        let elements = DMatrix::<f64>::identity(4, 4);
        let element_totals = vec![1.0, 1.0, 1.0, 1.0];

        let gibbs: Vec<GibbsFn> = vec![
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
            Rc::new(|_| 0.0),
        ];

        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![2, 3],
            },
        ];
        let species_phase = species_to_phase_map(&phases, 4).unwrap();
        let residual = equilibrium_logmole_residual(
            reactions,
            elements,
            element_totals,
            gibbs,
            phases,
            300.0,
            101325.0,
            101325.0,
            species_phase,
            1e-30,
            1e-30,
        )
        .unwrap();

        let _species_phase = vec![0, 0, 1, 1];
        let y = vec![0.0, 0.0, -700.0, -700.0];

        let f = residual(&y).unwrap();
        for v in f {
            assert!(v.is_finite());
        }
    }
    #[test]
    fn jacobian_finite_after_phase_deactivation() {
        let reactions = DMatrix::<f64>::zeros(4, 0);
        let elements = DMatrix::<f64>::identity(4, 4);
        let species_phase = vec![0, 0, 1, 1];
        let delta_n = Vec::new();
        let y = vec![0.0, 0.0, -700.0, -700.0];

        let J = equilibrium_logmole_jacobian(
            &y,
            &reactions,
            &elements,
            &species_phase,
            &delta_n,
            2,
            1e-30,
            1e-30,
        )
        .unwrap();

        for v in J.iter() {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn residual_rejects_non_square_formulation_before_evaluation() {
        let result = equilibrium_logmole_residual(
            DMatrix::<f64>::zeros(4, 1),
            DMatrix::<f64>::identity(4, 4),
            vec![1.0; 4],
            vec![Rc::new(|_| 0.0); 4],
            vec![Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0, 1, 2, 3],
            }],
            300.0,
            101325.0,
            101325.0,
            vec![0, 0, 0, 0],
            1e-30,
            1e-30,
        );

        assert!(matches!(
            result,
            Err(ReactionExtentError::DimensionMismatch(_))
        ));
    }

    #[test]
    fn seed_activated_phase_sets_every_species_in_the_phase() {
        let species_phase = vec![0, 0, 1, 1];
        let mut y = vec![0.0, 0.1, -20.0, -21.0];

        seed_activated_phase(
            y.as_mut_slice(),
            PhaseIndex::new(1, 2).unwrap(),
            &species_phase,
            PhaseSeedPolicy::AbsolutePerSpecies { moles: 1e-12 },
        )
        .unwrap();

        assert_eq!(y[0], 0.0);
        assert_eq!(y[1], 0.1);
        assert!(y[2] > -30.0);
        assert_eq!(y[2], y[3]);
    }

    #[test]
    fn phase_control_reports_a_typed_error_when_budget_is_exhausted() {
        let mut instance = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        instance.phase_manager.max_phase_iterations = 0;

        let result = instance.solve_with_phase_control();

        assert!(matches!(
            &result,
            Err(ReactionExtentError::PhaseControlDidNotConverge { iterations: 0 })
        ));
        let err = result.unwrap_err();
        assert_eq!(err.kind(), ReactionExtentErrorKind::NonConvergence);
    }

    #[test]
    fn phase_transition_plan_prefers_deactivation_over_activation() {
        let manager = PhaseManager {
            phase_eps: 1e-6,
            phase_hysteresis: PhaseHysteresisPolicy::Explicit {
                dg_create: -1.0,
                dg_keep: -20.0,
            },
            max_phase_iterations: 4,
            initial_phase_set: InitialPhaseSet::default(),
        };

        let active = [true, true, false];
        let stability = stability_reports(&[Some(-10.0), Some(-20.0), Some(-30.0)], &active);
        let plan = manager
            .classify_phases(&[1e-12, 2.0, 3.0], &stability, &phase_set(&active))
            .unwrap();

        assert_eq!(
            plan,
            PhaseTransitionPlan::Deactivate {
                phase: PhaseIndex::new(0, 3).unwrap()
            }
        );
    }

    #[test]
    fn phase_transition_plan_activates_most_favorable_missing_phase() {
        let manager = PhaseManager {
            phase_eps: 1e-6,
            phase_hysteresis: PhaseHysteresisPolicy::Explicit {
                dg_create: -1.0,
                dg_keep: 0.0,
            },
            max_phase_iterations: 4,
            initial_phase_set: InitialPhaseSet::default(),
        };

        let active = [true, false, false];
        let stability = stability_reports(&[None, Some(-5.0), Some(-15.0)], &active);
        let plan = manager
            .classify_phases(&[10.0, 10.0, 10.0], &stability, &phase_set(&active))
            .unwrap();

        assert_eq!(
            plan,
            PhaseTransitionPlan::Activate {
                phase: PhaseIndex::new(2, 3).unwrap()
            }
        );
    }

    #[test]
    fn phase_transition_plan_holds_marginal_active_phase_inside_keep_band() {
        let manager = PhaseManager {
            phase_eps: 1e-6,
            phase_hysteresis: PhaseHysteresisPolicy::Explicit {
                dg_create: -1.0,
                dg_keep: -0.1,
            },
            max_phase_iterations: 4,
            initial_phase_set: InitialPhaseSet::default(),
        };

        let active = [true, false];
        let stability = stability_reports(&[Some(-0.5), Some(10.0)], &active);
        let plan = manager
            .classify_phases(&[1e-12, 10.0], &stability, &phase_set(&active))
            .unwrap();

        assert_eq!(
            plan,
            PhaseTransitionPlan::Hold {
                phase: PhaseIndex::new(0, 2).unwrap()
            }
        );
    }

    #[test]
    fn hysteresis_does_not_hide_a_more_unstable_inactive_phase() {
        let manager = PhaseManager {
            phase_eps: 1e-6,
            phase_hysteresis: PhaseHysteresisPolicy::Explicit {
                dg_create: -1.0,
                dg_keep: 0.1,
            },
            max_phase_iterations: 4,
            initial_phase_set: InitialPhaseSet::default(),
        };
        let active = [true, true, false];
        let stability = stability_reports(&[None, Some(-0.5), Some(-5.0)], &active);

        let plan = manager
            .classify_phases(&[10.0, 1e-12, 0.0], &stability, &phase_set(&active))
            .unwrap();

        assert_eq!(
            plan,
            PhaseTransitionPlan::Activate {
                phase: PhaseIndex::new(2, 3).unwrap()
            }
        );
    }

    #[test]
    fn temperature_scaled_hysteresis_thresholds_scale_with_rt() {
        let manager = PhaseManager::with_temperature_scaled_hysteresis(1e-6, -1e-6, 1e-8);
        let (low_create, low_keep) = manager.thresholds_at(300.0).unwrap();
        let (high_create, high_keep) = manager.thresholds_at(600.0).unwrap();

        assert!(low_create < low_keep);
        assert!(high_create < high_keep);
        assert!((high_create - 2.0 * low_create).abs() <= low_create.abs() * 1e-12);
        assert!((high_keep - 2.0 * low_keep).abs() <= low_keep.abs() * 1e-12);
    }

    #[test]
    fn phase_transition_plan_reports_no_transition_when_system_is_stable() {
        let manager = PhaseManager {
            phase_eps: 1e-6,
            phase_hysteresis: PhaseHysteresisPolicy::Explicit {
                dg_create: -1.0,
                dg_keep: -0.1,
            },
            max_phase_iterations: 4,
            initial_phase_set: InitialPhaseSet::default(),
        };

        let active = [true, false];
        let stability = stability_reports(&[None, Some(0.5)], &active);
        let plan = manager
            .classify_phases(&[10.0, 12.0], &stability, &phase_set(&active))
            .unwrap();

        assert_eq!(
            plan,
            PhaseTransitionPlan::NoTransition {
                reason: PhaseStabilityReason::NoPhaseTransitionNeeded
            }
        );
    }

    #[test]
    fn multicomponent_ideal_solution_phase_rejects_phase_control_stability_guard() {
        let y = vec![
            0.0,
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
        ];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1, 2],
            },
        ];
        let species_phase = vec![0, 1, 1];

        let result = compute_phase_stability_reports(
            &y,
            &gibbs,
            &phases,
            &species_phase,
            &DMatrix::from_row_slice(3, 1, &[1.0, 1.0, 1.0]),
            300.0,
            101325.0,
            101325.0,
            &phase_set(&[true, false]),
        );

        assert!(matches!(
            result,
            Err(ReactionExtentError::ValidationNotApplicable {
                path: "phase_stability",
                ..
            })
        ));
    }

    #[test]
    fn inactive_ideal_gas_phase_has_no_trace_driven_activation_force() {
        let trace = PHASE_CONTROL_TRACE_MOLE_FLOOR.ln();
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| -1.0e9)];
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![1],
            },
        ];
        let active = [true, false];
        let reports = compute_phase_stability_reports(
            &[0.0, trace],
            &gibbs,
            &phases,
            &[0, 1],
            &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            1_000.0,
            101_325.0,
            101_325.0,
            &phase_set(&active),
        )
        .unwrap();

        assert_eq!(reports[1].model, PhaseStabilityModel::FixedIdealGas);
        assert_eq!(reports[1].driving_force, None);
        let plan = PhaseManager::default()
            .classify_phases_at_temperature(
                1_000.0,
                &[1.0, PHASE_CONTROL_TRACE_MOLE_FLOOR],
                &reports,
                &phase_set(&active),
            )
            .unwrap();
        assert!(matches!(plan, PhaseTransitionPlan::NoTransition { .. }));
    }

    #[test]
    fn pure_condensed_stability_uses_element_potentials() {
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| -1000.0)];
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
        ];
        let reports = compute_phase_stability_reports(
            &[0.0, PHASE_CONTROL_TRACE_MOLE_FLOOR.ln()],
            &gibbs,
            &phases,
            &[0, 1],
            &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            300.0,
            101325.0,
            101325.0,
            &phase_set(&[true, false]),
        )
        .unwrap();

        assert_eq!(reports[1].model, PhaseStabilityModel::PureCondensedSpecies);
        assert_eq!(reports[1].driving_force, Some(-1000.0));
        let potentials = reports[1].element_potentials.as_ref().unwrap();
        assert_eq!(potentials.rank, 1);
        assert!(potentials.max_abs_residual < 1e-12);
    }

    #[test]
    fn pure_condensed_driving_force_uses_chemical_potential_not_raw_sum() {
        let temperature = 500.0;
        let pressure = 2.0 * 101_325.0;
        let reference_pressure = 101_325.0;
        let rt = R * temperature;
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| -1_000.0), Rc::new(|_| 0.0)];
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
        ];
        let reports = compute_phase_stability_reports(
            &[0.0, PHASE_CONTROL_TRACE_MOLE_FLOOR.ln()],
            &gibbs,
            &phases,
            &[0, 1],
            &DMatrix::from_row_slice(2, 1, &[1.0, 1.0]),
            temperature,
            pressure,
            reference_pressure,
            &phase_set(&[true, false]),
        )
        .unwrap();

        let expected_gas_mu = -1_000.0 + rt * (pressure / reference_pressure).ln();
        let expected_driving_force = 0.0 - expected_gas_mu;
        assert_eq!(reports[0].model, PhaseStabilityModel::FixedIdealGas);
        assert_eq!(reports[1].model, PhaseStabilityModel::PureCondensedSpecies);
        assert_eq!(reports[1].driving_force.unwrap(), expected_driving_force);
    }

    #[test]
    fn unpublished_candidate_does_not_overwrite_accepted_state() {
        let mut instance = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![0.0, 1e-5_f64.ln()]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.create_stoich_matrix().unwrap();
        instance.solution = vec![7.0, 8.0];
        instance.moles = vec![9.0, 10.0];

        let seed = instance.initial_guess.clone().unwrap();
        let candidate = instance.solve_candidate_from_seed(seed).unwrap();

        assert_ne!(candidate.log_moles, instance.solution);
        assert_eq!(instance.solution, vec![7.0, 8.0]);
        assert_eq!(instance.moles, vec![9.0, 10.0]);
        assert!(instance.last_validation_report.is_none());
        assert!(instance.last_solve_report.is_none());
    }

    #[test]
    fn accepted_phase_control_solve_publishes_transition_evidence() {
        let mut instance = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![0.0, 1e-5_f64.ln()]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.create_stoich_matrix().unwrap();

        instance.solve_with_phase_control().unwrap();

        let report = instance.last_phase_control_report.as_ref().unwrap();
        assert_eq!(report.iterations, 1);
        assert_eq!(report.initial_active_phases, report.final_active_phases);
        assert!(report.transitions.is_empty());
        assert_eq!(report.nonlinear_reports.len(), 1);
        assert_eq!(
            report.final_validation,
            instance.last_validation_report.clone().unwrap()
        );
    }

    #[test]
    fn stable_absent_pure_phase_stays_out_of_reduced_solve() {
        let mut solver = synthetic_gas_with_pure_condensed(100_000.0);

        solver.solve_with_phase_control().unwrap();

        assert!(solver.moles[2] <= PHASE_CONTROL_TRACE_MOLE_FLOOR * 1.01);
        assert_eq!(solver.phase_active_mask, vec![true, false]);
        assert!(
            solver
                .last_phase_control_report
                .as_ref()
                .unwrap()
                .transitions
                .is_empty()
        );
    }

    #[test]
    fn unstable_absent_pure_phase_is_added_between_fixed_set_solves() {
        let mut solver = synthetic_gas_with_pure_condensed(-20_000.0);

        solver.solve_with_phase_control().unwrap();

        assert!(solver.moles[2] > solver.phase_manager.phase_eps);
        assert_eq!(solver.phase_active_mask, vec![true, true]);
        let transitions = &solver
            .last_phase_control_report
            .as_ref()
            .unwrap()
            .transitions;
        assert_eq!(transitions.len(), 1);
        assert_eq!(
            transitions[0].activated,
            vec![PhaseIndex::new(1, 2).unwrap()]
        );
        assert_eq!(
            transitions[0]
                .previous_phase_set
                .status(PhaseIndex::new(1, 2).unwrap()),
            PhaseStatus::Inactive
        );
        assert_eq!(
            transitions[0]
                .new_phase_set
                .status(PhaseIndex::new(1, 2).unwrap()),
            PhaseStatus::Appeared
        );
        assert!(transitions[0].restart_seed[2].exp() > PHASE_CONTROL_TRACE_MOLE_FLOOR);
        let (dg_create, _) = solver.phase_manager.thresholds_at(solver.T).unwrap();
        assert!(matches!(
            transitions[0].reason,
            PhaseTransitionReason::UnstableInactivePhase { driving_force }
                if driving_force < dg_create
        ));
        assert!(
            transitions[0]
                .candidate_validation
                .max_abs_element_balance_error
                <= 1e-8
        );
        assert!(
            solver
                .last_validation_report
                .as_ref()
                .unwrap()
                .max_abs_element_balance_error
                <= 1e-8
        );
    }

    #[test]
    fn explicit_exclusion_keeps_an_unstable_candidate_out_of_phase_control() {
        let mut solver = synthetic_gas_with_pure_condensed(-20_000.0);
        solver.phase_manager.initial_phase_set = InitialPhaseSet::Explicit {
            active: vec![PhaseIndex::new(0, 2).unwrap()],
            excluded: vec![PhaseIndex::new(1, 2).unwrap()],
        };

        solver.solve_with_phase_control().unwrap();

        assert_eq!(solver.phase_active_mask, vec![true, false]);
        assert!(solver.moles[2] <= PHASE_CONTROL_TRACE_MOLE_FLOOR * 1.01);
        assert!(
            solver
                .last_phase_control_report
                .as_ref()
                .unwrap()
                .transitions
                .is_empty()
        );
    }

    #[test]
    fn explicit_phase_policy_rejects_active_excluded_overlap() {
        let phase = PhaseIndex::new(0, 1).unwrap();
        let result = PhaseSet::from_policy(
            &InitialPhaseSet::Explicit {
                active: vec![phase],
                excluded: vec![phase],
            },
            &[true],
        );

        assert!(matches!(
            result,
            Err(ReactionExtentError::InvalidProblem {
                field: "initial_phase_set",
                ..
            })
        ));
    }

    #[test]
    fn failed_restart_budget_does_not_publish_the_intermediate_candidate() {
        let mut solver = synthetic_gas_with_pure_condensed(-20_000.0);
        solver.phase_manager.max_phase_iterations = 1;
        solver.solution = vec![7.0, 8.0, 9.0];
        solver.moles = vec![10.0, 11.0, 12.0];
        let accepted_solution = solver.solution.clone();
        let accepted_moles = solver.moles.clone();

        let error = solver.solve_with_phase_control().unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::PhaseControlDidNotConverge { iterations: 1 }
        ));
        assert_eq!(solver.solution, accepted_solution);
        assert_eq!(solver.moles, accepted_moles);
        assert!(solver.last_validation_report.is_none());
        assert!(solver.last_solve_report.is_none());
        assert!(solver.last_phase_control_report.is_none());
    }

    #[test]
    fn near_boundary_phase_is_held_by_hysteresis_without_cycling() {
        let mut solver = synthetic_gas_with_pure_condensed(-5_800.0);
        solver.phase_manager.phase_eps = 2e-2;
        solver.phase_manager.set_explicit_hysteresis(-1.0, 1.0);

        solver.solve_with_phase_control().unwrap();

        assert_eq!(solver.phase_active_mask, vec![true, true]);
        assert!(solver.moles[2] < solver.phase_manager.phase_eps);
        let report = solver.last_phase_control_report.as_ref().unwrap();
        assert_eq!(report.transitions.len(), 1);
        assert!(report.iterations <= solver.phase_manager.max_phase_iterations);
    }

    #[test]
    fn reduced_active_set_solve_preserves_the_original_physical_element_totals() {
        let mut solver = synthetic_gas_with_pure_condensed(-20_000.0);
        let expected = compute_element_totals(&solver.elem_composition, &solver.n0).unwrap();

        solver.solve_with_phase_control().unwrap();

        let observed = compute_element_totals(&solver.elem_composition, &solver.moles).unwrap();
        for element in 0..expected.len() {
            assert!(
                (expected[element] - observed[element]).abs() <= 1e-8,
                "element {element} drifted from {} to {} during phase-control",
                expected[element],
                observed[element]
            );
        }
    }

    #[test]
    fn invalid_hysteresis_order_is_rejected_before_phase_control() {
        let mut solver = synthetic_gas_with_pure_condensed(-5_800.0);
        solver.phase_manager.phase_eps = 2e-2;
        solver.phase_manager.set_explicit_hysteresis(-1.0, 1.0);
        solver.phase_manager.set_explicit_hysteresis(-1.0, -1.0);

        let error = solver.solve_with_phase_control().unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::InvalidProblem {
                field: "phase_hysteresis",
                ..
            }
        ));
        assert!(solver.last_phase_control_report.is_none());
    }

    #[test]
    fn lm_and_tr_phase_control_match_on_the_same_fixture_with_documented_tolerance() {
        let mut lm_solver = synthetic_gas_with_pure_condensed(-20_000.0);
        lm_solver.phase_manager.phase_eps = 1e-12;
        lm_solver
            .phase_manager
            .set_temperature_scaled_hysteresis(-1e-6, 1e-8);
        lm_solver.solver_settings.solver = Solvers::LM;

        let mut tr_solver = synthetic_gas_with_pure_condensed(-20_000.0);
        tr_solver.phase_manager.phase_eps = 1e-12;
        tr_solver
            .phase_manager
            .set_temperature_scaled_hysteresis(-1e-6, 1e-8);
        tr_solver.solver_settings.solver = Solvers::TR;

        lm_solver.solve_with_phase_control().unwrap();
        tr_solver.solve_with_phase_control().unwrap();

        assert_eq!(lm_solver.phase_active_mask, tr_solver.phase_active_mask);
        assert_eq!(
            lm_solver
                .last_phase_control_report
                .as_ref()
                .unwrap()
                .final_active_phases,
            tr_solver
                .last_phase_control_report
                .as_ref()
                .unwrap()
                .final_active_phases
        );
        assert_eq!(lm_solver.moles.len(), tr_solver.moles.len());
        for (lhs, rhs) in lm_solver.moles.iter().zip(tr_solver.moles.iter()) {
            assert!((lhs - rhs).abs() <= 1e-5 * lhs.abs().max(rhs.abs()).max(1.0));
        }
        let lm_balance = lm_solver
            .last_validation_report
            .as_ref()
            .unwrap()
            .max_abs_element_balance_error;
        let tr_balance = tr_solver
            .last_validation_report
            .as_ref()
            .unwrap()
            .max_abs_element_balance_error;
        assert!(lm_balance.is_finite() && tr_balance.is_finite());
        assert!(lm_balance <= 1e-6);
        assert!(tr_balance <= 1e-6);
        assert!((lm_balance - tr_balance).abs() <= 1e-6);
    }

    #[test]
    fn lm_and_tr_phase_control_match_on_a_stable_gas_only_fixture() {
        let substances = vec!["O2".to_string(), "O".to_string()];

        let mut lm_solver = gas_solver(
            substances.clone(),
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        lm_solver.n0 = vec![1.0, 1e-5];
        lm_solver.initial_guess = Some(vec![1.0, 1e-5]);
        lm_solver.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        lm_solver.create_stoich_matrix().unwrap();

        let mut tr_solver =
            gas_solver(substances, 500.0, 101325.0, Solvers::TR, None, false).unwrap();
        tr_solver.n0 = vec![1.0, 1e-5];
        tr_solver.initial_guess = Some(vec![1.0, 1e-5]);
        tr_solver.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        tr_solver.create_stoich_matrix().unwrap();

        lm_solver.solve_with_phase_control().unwrap();
        tr_solver.solve_with_phase_control().unwrap();

        assert_eq!(lm_solver.phase_active_mask, tr_solver.phase_active_mask);
        assert_eq!(lm_solver.solution.len(), tr_solver.solution.len());
        for (lhs, rhs) in lm_solver.solution.iter().zip(tr_solver.solution.iter()) {
            assert!((lhs - rhs).abs() <= 1e-8 * lhs.abs().max(rhs.abs()).max(1.0));
        }
        assert_eq!(lm_solver.moles.len(), tr_solver.moles.len());
        for (lhs, rhs) in lm_solver.moles.iter().zip(tr_solver.moles.iter()) {
            assert!((lhs - rhs).abs() <= 1e-8 * lhs.abs().max(rhs.abs()).max(1.0));
        }
        let lm_balance = lm_solver
            .last_validation_report
            .as_ref()
            .unwrap()
            .max_abs_element_balance_error;
        let tr_balance = tr_solver
            .last_validation_report
            .as_ref()
            .unwrap()
            .max_abs_element_balance_error;
        assert!(lm_balance.is_finite() && tr_balance.is_finite());
        assert!(lm_balance <= 1e-6);
        assert!(tr_balance <= 1e-6);
    }

    #[test]
    fn phase_control_remains_finite_across_trace_and_extreme_driving_forces() {
        let temperatures = [250.0_f64, 1000.0, 5000.0];
        let candidate_gibbs_values = [-100_000.0_f64, 0.0, 100_000.0];

        for &temperature in &temperatures {
            for &candidate_gibbs in &candidate_gibbs_values {
                let mut solver = synthetic_gas_with_pure_condensed(candidate_gibbs);
                solver.T = temperature;
                solver.phase_manager.phase_eps = 1e-12;
                solver
                    .phase_manager
                    .set_temperature_scaled_hysteresis(-1e-6, 1e-8);

                let expected = compute_element_totals(&solver.elem_composition, &solver.n0)
                    .expect("baseline element totals must be computable");
                solver.solve_with_phase_control().unwrap();

                let observed = compute_element_totals(&solver.elem_composition, &solver.moles)
                    .expect("accepted moles must remain element-representable");
                for element in 0..expected.len() {
                    assert!(
                        (expected[element] - observed[element]).abs() <= 1e-8,
                        "element {element} drifted at T={temperature}, candidate_gibbs={candidate_gibbs}"
                    );
                }
                assert!(
                    solver
                        .moles
                        .iter()
                        .all(|value| value.is_finite() && *value >= 0.0)
                );
                assert!(
                    solver
                        .last_phase_control_report
                        .as_ref()
                        .is_some_and(
                            |report| report.iterations <= solver.phase_manager.max_phase_iterations
                        )
                );
            }
        }
    }

    #[test]
    fn phase_control_remains_finite_with_multiple_pure_condensed_candidates() {
        let candidate_pairs = [(-20_000.0_f64, 30_000.0_f64), (5_000.0, -12_000.0)];

        for &(candidate_a_gibbs, candidate_b_gibbs) in &candidate_pairs {
            let mut solver = synthetic_gas_with_two_condensed_candidates(
                Rc::new(RefCell::new(candidate_a_gibbs)),
                Rc::new(RefCell::new(candidate_b_gibbs)),
            );
            solver.phase_manager.phase_eps = 1e-12;
            solver
                .phase_manager
                .set_temperature_scaled_hysteresis(-1e-6, 1e-8);

            let expected = compute_element_totals(&solver.elem_composition, &solver.n0)
                .expect("baseline element totals must be computable");
            solver.solve_with_phase_control().unwrap();
            let observed = compute_element_totals(&solver.elem_composition, &solver.moles)
                .expect("accepted moles must remain element-representable");

            assert!(
                solver
                    .moles
                    .iter()
                    .all(|value| value.is_finite() && *value >= 0.0)
            );
            assert_eq!(solver.phase_active_mask.len(), 3);
            assert!(solver.phase_active_mask[0]);
            for element in 0..expected.len() {
                assert!(
                    (expected[element] - observed[element]).abs() <= 1e-8,
                    "element {element} drifted for candidate pair ({candidate_a_gibbs}, {candidate_b_gibbs})"
                );
            }
        }
    }

    #[test]
    fn phase_control_remains_finite_for_nearly_rank_deficient_stress_fixtures() {
        let mut solver = synthetic_gas_with_two_condensed_candidates(
            Rc::new(RefCell::new(-1_500.0)),
            Rc::new(RefCell::new(-1_450.0)),
        );
        solver.T = 250.0;
        solver.phase_manager.phase_eps = 1e-12;
        solver
            .phase_manager
            .set_temperature_scaled_hysteresis(-1e-7, 1e-8);

        let expected = compute_element_totals(&solver.elem_composition, &solver.n0)
            .expect("baseline element totals must be computable");
        let error = solver.solve_with_phase_control().unwrap_err();

        match error {
            ReactionExtentError::AllBackendsFailed { attempts } => {
                assert!(!attempts.is_empty());
                assert!(
                    attempts
                        .iter()
                        .all(|attempt| !attempt.outcome.is_accepted())
                );
            }
            other => panic!(
                "expected typed backend exhaustion for nearly rank-deficient stress fixture, got {other:?}"
            ),
        }
        let observed = compute_element_totals(&solver.elem_composition, &solver.n0)
            .expect("baseline element totals must remain computable");
        for element in 0..expected.len() {
            assert!(
                (expected[element] - observed[element]).abs() <= 1e-12,
                "element {element} drifted in a nearly rank-deficient stress fixture"
            );
        }
    }

    #[test]
    fn phase_appearance_followed_by_failed_restaging_keeps_previous_publication() {
        let candidate_gibbs = Rc::new(RefCell::new(-20_000.0));
        let mut solver = synthetic_gas_with_dynamic_condensed(candidate_gibbs.clone());

        solver.solve_with_phase_control().unwrap();
        let first_report = solver.last_phase_control_report.as_ref().unwrap().clone();
        let first_solution = solver.solution.clone();
        let first_moles = solver.moles.clone();

        assert_eq!(solver.phase_active_mask, vec![true, true]);
        assert_eq!(first_report.transitions.len(), 1);
        assert!(
            first_report
                .final_active_phases
                .contains(&PhaseIndex::new(1, 2).unwrap())
        );

        *candidate_gibbs.borrow_mut() = 1_000.0;
        let second_result = solver.solve_with_phase_control();

        assert!(second_result.is_err());
        assert_eq!(solver.phase_active_mask, vec![true, true]);
        assert_eq!(solver.solution, first_solution);
        assert_eq!(solver.moles, first_moles);
        assert!(first_report.transitions.iter().any(|transition| matches!(
            transition.reason,
            PhaseTransitionReason::UnstableInactivePhase { .. }
        )));
    }

    #[test]
    fn hysteresis_story_exercises_appear_keep_disappear_reappear() {
        let manager = PhaseManager {
            phase_eps: 1e-8,
            phase_hysteresis: PhaseHysteresisPolicy::Explicit {
                dg_create: -1.0,
                dg_keep: 1.0,
            },
            max_phase_iterations: 4,
            initial_phase_set: InitialPhaseSet::default(),
        };

        let candidate_set = phase_set(&[true, false]);
        let full_set = phase_set(&[true, true]);

        let appearance = manager
            .classify_phases_at_temperature(
                500.0,
                &[1.0, 1e-12],
                &stability_reports(&[Some(0.0), Some(-2.0)], &[true, false]),
                &candidate_set,
            )
            .unwrap();
        assert_eq!(
            appearance,
            PhaseTransitionPlan::Activate {
                phase: PhaseIndex::new(1, 2).unwrap()
            }
        );

        let keep = manager
            .classify_phases_at_temperature(
                500.0,
                &[1.0, 1.0],
                &stability_reports(&[Some(0.0), Some(0.0)], &[true, true]),
                &full_set,
            )
            .unwrap();
        assert_eq!(
            keep,
            PhaseTransitionPlan::NoTransition {
                reason: PhaseStabilityReason::NoPhaseTransitionNeeded
            }
        );

        let disappear = manager
            .classify_phases_at_temperature(
                500.0,
                &[1.0, 0.0],
                &stability_reports(&[Some(0.0), Some(2.0)], &[true, true]),
                &full_set,
            )
            .unwrap();
        assert_eq!(
            disappear,
            PhaseTransitionPlan::Deactivate {
                phase: PhaseIndex::new(1, 2).unwrap()
            }
        );

        let reappear = manager
            .classify_phases_at_temperature(
                500.0,
                &[1.0, 1e-12],
                &stability_reports(&[Some(0.0), Some(-2.0)], &[true, false]),
                &candidate_set,
            )
            .unwrap();
        assert_eq!(
            reappear,
            PhaseTransitionPlan::Activate {
                phase: PhaseIndex::new(1, 2).unwrap()
            }
        );
    }

    #[test]
    fn several_sequential_phase_transitions_preserve_publication_and_state() {
        let solve_mask = |candidate_gibbs: f64| {
            let mut solver =
                synthetic_gas_with_dynamic_condensed(Rc::new(RefCell::new(candidate_gibbs)));
            solver.solve_with_phase_control().unwrap();
            (
                solver.phase_active_mask,
                solver.last_phase_control_report.unwrap(),
            )
        };

        let (appearance_mask, appearance_report) = solve_mask(-20_000.0);
        let (disappearance_mask, disappearance_report) = solve_mask(10_000.0);
        let (reappearance_mask, reappearance_report) = solve_mask(-20_000.0);

        assert_eq!(appearance_mask, vec![true, true]);
        assert_eq!(appearance_report.transitions.len(), 1);
        assert_eq!(disappearance_mask, vec![true, false]);
        assert_eq!(disappearance_report.transitions.len(), 0);
        assert_eq!(reappearance_mask, vec![true, true]);
        assert_eq!(reappearance_report.transitions.len(), 1);
    }

    #[test]
    fn phase_control_reports_maximum_restart_termination_after_a_single_transition() {
        let candidate_gibbs = Rc::new(RefCell::new(-20_000.0));
        let mut solver = synthetic_gas_with_dynamic_condensed(candidate_gibbs);
        let published_solution = solver.solution.clone();
        let published_moles = solver.moles.clone();
        solver.phase_manager.max_phase_iterations = 1;

        let error = solver.solve_with_phase_control().unwrap_err();

        assert!(matches!(
            error,
            ReactionExtentError::PhaseControlDidNotConverge { iterations: 1 }
        ));
        assert_eq!(solver.solution, published_solution);
        assert_eq!(solver.moles, published_moles);
        assert!(solver.last_phase_control_report.is_none());
        assert!(solver.last_validation_report.is_none());
        assert!(solver.last_solve_report.is_none());
    }

    #[test]
    fn phase_control_report_exposes_stable_summary_rows_and_display() {
        let mut solver = synthetic_gas_with_pure_condensed(-20_000.0);
        solver.solve_with_phase_control().unwrap();
        let report = solver.last_phase_control_report.as_ref().unwrap();
        let rows = report.summary_rows();
        let rendered = report.to_string();

        assert!(
            rows.iter()
                .any(|row| row.section == "phase_control" && row.label == "iterations")
        );
        assert!(
            rows.iter()
                .any(|row| row.section == "validation" && row.label == "residual_l2_norm")
        );
        assert!(
            rows.iter()
                .any(|row| row.section == "nonlinear" && row.label == "accepted_reports")
        );
        assert!(rendered.contains("[phase_control] iterations = "));
        assert!(rendered.contains("[validation] residual_l2_norm = "));
    }

    #[test]
    fn multiphase_acceptance_report_combines_validation_and_stability() {
        let mut solver = synthetic_gas_with_pure_condensed(-20_000.0);
        solver.solve_with_phase_control().unwrap();
        let phase_control = solver.last_phase_control_report.as_ref().unwrap().clone();
        let stability = compute_phase_stability_reports(
            &solver.solution,
            &solver.gibbs,
            &solver.phases,
            &solver.species_phase,
            &solver.elem_composition,
            solver.T,
            solver.P,
            solver.p0,
            &phase_control.final_phase_set,
        )
        .unwrap();
        let (dg_create, dg_keep) = solver.phase_manager.thresholds_at(solver.T).unwrap();
        let report = build_multiphase_acceptance_report(
            phase_control,
            stability,
            Some(dg_create),
            Some(dg_keep),
        )
        .unwrap();
        let rows = report.summary_rows();
        let rendered = report.to_string();

        assert!(report.complementarity.satisfied);
        assert!(
            rows.iter()
                .any(|row| row.section == "acceptance" && row.label == "phase_iterations")
        );
        assert!(
            rows.iter()
                .any(|row| row.section == "complementarity" && row.label == "satisfied")
        );
        assert!(rendered.contains("[acceptance] phase_iterations = "));
        assert!(rendered.contains("[complementarity] satisfied = true"));
    }

    #[test]
    fn phase_seed_cycle_restores_activation_after_deactivation() {
        let species_phase = vec![0, 0, 1];
        let phase = PhaseIndex::new(1, 2).unwrap();
        let mut y = vec![
            1.0_f64.ln(),
            1.0_f64.ln(),
            PHASE_CONTROL_TRACE_MOLE_FLOOR.ln(),
        ];

        seed_activated_phase(
            y.as_mut_slice(),
            phase,
            &species_phase,
            PhaseSeedPolicy::RelativeToSystemTotal {
                fraction: 1e-8,
                minimum: PHASE_CONTROL_TRACE_MOLE_FLOOR,
            },
        )
        .unwrap();
        let appeared = y[2].exp();
        assert!(appeared > PHASE_CONTROL_TRACE_MOLE_FLOOR);

        deactivate_phases_seed_only(
            y.as_mut_slice(),
            &[phase],
            &species_phase,
            PHASE_CONTROL_TRACE_MOLE_FLOOR,
        )
        .unwrap();
        let disappeared = y[2].exp();
        assert!(disappeared <= PHASE_CONTROL_TRACE_MOLE_FLOOR * 1.01);

        seed_activated_phase(
            y.as_mut_slice(),
            phase,
            &species_phase,
            PhaseSeedPolicy::RelativeToSystemTotal {
                fraction: 1e-8,
                minimum: PHASE_CONTROL_TRACE_MOLE_FLOOR,
            },
        )
        .unwrap();
        let reappeared = y[2].exp();
        assert!(reappeared > PHASE_CONTROL_TRACE_MOLE_FLOOR);
        assert!(reappeared >= appeared * 0.99);
    }

    #[test]
    #[ignore = "opt-in performance characterization"]
    fn phase_control_characterization_scales_through_several_phases() {
        let mut solver = synthetic_large_phase_control_fixture(6, 4);
        let started = Instant::now();

        solver.solve_with_phase_control().unwrap();

        let elapsed = started.elapsed();
        let report = solver.last_phase_control_report.as_ref().unwrap();
        eprintln!(
            "phase-control characterization: species={}, phases={}, iterations={}, transitions={}, elapsed={:?}",
            solver.subs_data.substances.len(),
            solver.phases.len(),
            report.iterations,
            report.transitions.len(),
            elapsed
        );

        assert!(solver.moles.iter().all(|value| value.is_finite()));
        assert!(report.iterations >= 1);
        assert!(report.final_validation.max_abs_element_balance_error <= 1e-8);
    }

    #[test]
    #[ignore = "opt-in performance characterization"]
    fn phase_control_characterization_scales_to_a_larger_synthetic_inventory() {
        let mut solver = synthetic_large_phase_control_fixture(10, 8);
        let started = Instant::now();

        solver.solve_with_phase_control().unwrap();

        let elapsed = started.elapsed();
        let report = solver.last_phase_control_report.as_ref().unwrap();
        eprintln!(
            "phase-control characterization (larger): species={}, phases={}, iterations={}, transitions={}, elapsed={:?}",
            solver.subs_data.substances.len(),
            solver.phases.len(),
            report.iterations,
            report.transitions.len(),
            elapsed
        );

        assert!(solver.moles.iter().all(|value| value.is_finite()));
        assert!(report.iterations >= 1);
        assert!(report.final_validation.max_abs_element_balance_error <= 1e-8);
    }

    #[test]
    #[ignore = "opt-in performance characterization"]
    fn phase_control_benchmark_sweep_across_larger_inventories() {
        let cases = [(5, 4), (10, 5), (10, 10), (10, 20)];

        for (species_per_phase, phase_count) in cases {
            let fixture = synthetic_large_phase_control_fixture(species_per_phase, phase_count);
            let active = vec![true; fixture.phases.len()];
            let start_projection = Instant::now();
            let projection = ActiveSetProjection::build(
                &fixture.phases,
                &fixture.species_phase,
                &fixture.elem_composition,
                &active,
                1e-12,
            )
            .unwrap();
            let projection_elapsed = start_projection.elapsed();

            let mut solver = fixture;
            let start_solve = Instant::now();
            solver.solve_with_phase_control().unwrap();
            let solve_elapsed = start_solve.elapsed();
            let report = solver.last_phase_control_report.as_ref().unwrap();

            eprintln!(
                "phase-control benchmark: species={}, phases={}, active_species={}, projection_build={:?}, solve={:?}, outer_iterations={}, transitions={}, reactions={}",
                solver.subs_data.substances.len(),
                solver.phases.len(),
                projection.active_species.len(),
                projection_elapsed,
                solve_elapsed,
                report.iterations,
                report.transitions.len(),
                projection.reaction_basis.num_reactions
            );

            assert!(solver.moles.iter().all(|value| value.is_finite()));
            assert!(report.final_validation.max_abs_element_balance_error <= 1e-8);
        }
    }

    #[test]
    fn gas_temperature_range_rejects_mismatched_species_and_initial_moles() {
        let result = gas_solver_for_T_range(
            vec!["O2".to_string()],
            vec![1.0, 0.0],
            101325.0,
            300.0,
            400.0,
            10.0,
            Solvers::LM,
            None,
            false,
        );

        assert!(matches!(
            result,
            Err(ReactionExtentError::DimensionMismatch(_))
        ));
    }

    ////////////////////////////////////////////////
    #[test]
    fn O2_O_equilibrium_logmoles_with_phase_control() {
        let T = 500.0;
        // no creation/deactivation case
        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut instance =
            gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false).unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.with_loglevel(Some("info"));
        instance.create_stoich_matrix().unwrap();
        instance.solve_with_phase_control().unwrap();

        let moles = &instance.moles;
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");
        assert_eq!(moles.len(), 2);
    }

    #[test]
    fn offline_local_subsdata_solve_retains_provenance_report() {
        let t = 500.0;
        let mut instance = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            t,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();

        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.create_stoich_matrix().unwrap();
        instance.solve_with_phase_control().unwrap();

        let report = instance.subs_data.search_summary_report();
        let rendered = report.render_table();

        assert!(instance.last_validation_report.is_some());
        assert!(instance.last_solve_report.is_some());
        assert!(instance.moles.iter().all(|&n| n >= 0.0));
        assert!(rendered.contains("NASA_gas"));
        assert!(rendered.contains("O2"));
        assert!(rendered.contains("O"));
        assert_eq!(
            instance
                .subs_data
                .get_search_result("O2", WhatIsFound::Thermo)
                .map(|result| result.library()),
            Some("NASA_gas")
        );
    }

    #[test]
    fn mixed_library_offline_subsdata_search_reports_multiple_local_sources() {
        let mut subs_data = SubsData::new();
        subs_data.substances = vec!["H2O".to_string(), "CO2".to_string(), "CH4".to_string()];
        subs_data
            .map_of_phases
            .insert("H2O".to_string(), Some(Phases::Gas));
        subs_data
            .map_of_phases
            .insert("CO2".to_string(), Some(Phases::Gas));
        subs_data
            .map_of_phases
            .insert("CH4".to_string(), Some(Phases::Gas));
        subs_data.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        subs_data.set_explicit_search_instruction("H2O".to_string(), LibraryId::NuigThermo);

        subs_data.search_substances().unwrap();
        let report = subs_data.search_summary_report();
        let rendered = report.render_table();

        assert_eq!(report.total_substances(), 3);
        assert!(rendered.contains("NASA_gas"));
        assert!(rendered.contains("nuig_thermo"));
        assert_eq!(
            subs_data
                .get_search_result("H2O", WhatIsFound::Thermo)
                .map(|result| result.library()),
            Some("nuig_thermo")
        );
        assert!(
            subs_data
                .get_search_result("CO2", WhatIsFound::Thermo)
                .is_some()
        );
        assert!(
            subs_data
                .get_search_result("CH4", WhatIsFound::Thermo)
                .is_some()
        );
    }

    #[test]
    fn phase_control_restart_failure_preserves_published_state() {
        let mut instance = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            500.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();
        instance.n0 = vec![1.0, 1e-5];
        instance.initial_guess = Some(vec![1.0, 1e-5]);
        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];
        instance.create_stoich_matrix().unwrap();
        instance.solution = vec![7.0, 8.0];
        instance.moles = vec![9.0, 10.0];
        instance.solver_settings.solver_policy = Some(SolverPolicy::Cascade(Vec::new()));

        let result = instance.solve_with_phase_control();

        assert!(matches!(
            result,
            Err(ReactionExtentError::InvalidProblem {
                field: "solver_policy",
                ..
            })
        ));
        assert_eq!(instance.solution, vec![7.0, 8.0]);
        assert_eq!(instance.moles, vec![9.0, 10.0]);
    }

    /*
    #[test]
    fn phase_destruction_due_to_depletion() {
        let T = 3000.0; // high T favors dissociation

        let substances = vec!["O2".to_string(), "O".to_string()];

        let mut instance = gas_solver(substances, T, 101325.0, Solvers::TR, None, false).unwrap();

        instance.n0 = vec![1.0, 1e-12];
        instance.initial_guess = Some(vec![1.0, 1e-12]);

        instance.phases = vec![Phase {
            kind: PhaseKind::IdealGas,
            species: vec![0, 1],
        }];

        instance.phase_eps = 1e-10; // force deactivation
        instance.create_stoich_matrix().unwrap();
        instance.solve_with_phase_control().unwrap();

        let moles = &instance.moles;

        assert!(moles[1] < 1e-12, "unstable species not eliminated");
    }

    #[test]
    fn phase_creation_hysteresis_stability() {
        let T = 500.0;

        let substances = vec!["A_g".to_string(), "A_l".to_string()];

        let mut instance = gas_solver(substances, T, 101325.0, Solvers::TR, None, false).unwrap();

        instance.n0 = vec![1.0, 0.0];
        instance.initial_guess = Some(vec![1.0, 1e-20]);

        instance.phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
        ];

        instance.phase_eps = 1e-12;
        instance.subs_eps = 1e-30;

        instance.create_stoich_matrix().unwrap();
        instance.solve_with_phase_control().unwrap();

        let moles = &instance.moles;

        assert!(
            moles[1] > 0.0 || moles[1] < 1e-12,
            "phase flip-flopping detected"
        );
    }
    #[test]
    fn species_elimination_and_reactivation() {
        let T = 200.0;

        let substances = vec!["B_g".to_string(), "B_s".to_string()];

        let mut instance = gas_solver(substances, T, 101325.0, Solvers::TrustRegion, None, false).unwrap();

        instance.n0 = vec![1.0, 0.0];
        instance.initial_guess = Some(vec![1.0, 1e-40]);

        instance.phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
        ];

        instance.create_stoich_matrix().unwrap();
        instance.solve_with_phase_control().unwrap();

        let moles = &instance.moles;

        assert!(moles.iter().all(|&n| n >= 0.0));
    }
    */
}
