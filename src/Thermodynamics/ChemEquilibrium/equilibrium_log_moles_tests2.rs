#![allow(deprecated)]
#[cfg(test)]
mod tests {
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
        continuation_seed_for_point, equilibrium_logmole_jacobian2, equilibrium_logmole_residual2,
        evaluate_equilibrium_logmole_jacobian, evaluate_equilibrium_logmole_residual,
        recoverable_backend_failure_kind, scale_jacobian_rows, scale_residual_rows, scaled_jacobian,
        scaled_residual, temperature_failure, validate_logmole_system_dimensions,
        validate_residual_conditions, ContinuationSeedPolicy, EquilibriumLogMoles,
        EquilibriumSolveCandidate, GibbsFn, Phase, PhaseKind, Solvers,
        TemperatureSolveFailure, TemperatureSolveSnapshot,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::EquilibriumConstantCrossValidationStatus;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::EquilibriumConstantValidationMode;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{ReactionExtentError, SolveError};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::SolverParams;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        EquilibriumSolveReport, SolverAttemptFailureKind, SolverAttemptOutcome,
        SolverAttemptReport, SolverBackend, SolverCascadeBudget, SolverPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
        build_multiphase_acceptance_report, deactivate_phases_seed_only, seed_activated_phase,
        InitialPhaseSet, PhaseManager, PhaseSeedPolicy, PhaseSet, PhaseStabilityReport,
        PhaseTransitionPlan, PHASE_CONTROL_TRACE_MOLE_FLOOR,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseControlledSolveReport;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ids::PhaseIndex;
    use nalgebra::DMatrix;
    use std::collections::HashMap;
    use std::rc::Rc;

    // -----------------------------------------------------------------------
    // D.4 — accepted_solution() after solve_with_phase_control()
    // -----------------------------------------------------------------------
    #[test]
    fn accepted_solution_after_phase_control_returns_ok_for_single_phase() {
        // O2 dissociation, single gas phase — phase control should succeed
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.n0 = vec![1.0, 0.0];
        solver.T = 3000.0;
        solver.P = 101_325.0;
        solver.p0 = 101_325.0;
        solver.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        solver.gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        solver.phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        solver.species_phase = vec![0, 0];
        solver.solver_settings.solver = Solvers::LM;
        solver.create_stoich_matrix().unwrap();
        solver.solve_with_phase_control().unwrap();
        let solution = solver.accepted_solution();
        assert!(solution.is_ok());
        let sol = solution.unwrap();
        assert!(!sol.moles().is_empty());
        assert_eq!(sol.moles().len(), 2);
    }

    #[test]
    fn accepted_solution_after_phase_control_preserves_element_totals() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.n0 = vec![1.0, 0.0];
        solver.T = 3000.0;
        solver.P = 101_325.0;
        solver.p0 = 101_325.0;
        solver.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        solver.gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        solver.phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        solver.species_phase = vec![0, 0];
        solver.solver_settings.solver = Solvers::LM;
        solver.create_stoich_matrix().unwrap();
        solver.solve_with_phase_control().unwrap();
        let solution = solver.accepted_solution().unwrap();
        // total O = 2 * n_O2 + 1 * n_O should equal 2.0
        let total_o = 2.0 * solution.moles()[0] + 1.0 * solution.moles()[1];
        assert!((total_o - 2.0).abs() < 1e-6);
    }

    #[test]
    fn accepted_solution_after_phase_control_fails_without_solve() {
        let solver = EquilibriumLogMoles::empty();
        let result = solver.accepted_solution();
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.5 — current_problem_snapshot() — indirect test
    // -----------------------------------------------------------------------
    #[test]
    fn from_problem_creates_solver_that_solves_and_returns_accepted_solution() {
        // from_problem is the only public path that exercises snapshot internals.
        // We verify the round-trip: from_problem → solve → accepted_solution.
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.n0 = vec![1.0, 0.0];
        solver.T = 3000.0;
        solver.P = 101_325.0;
        solver.p0 = 101_325.0;
        solver.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        solver.gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        solver.phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        solver.species_phase = vec![0, 0];
        solver.solver_settings.solver = Solvers::LM;
        solver.create_stoich_matrix().unwrap();
        solver.solve().unwrap();
        let solution = solver.accepted_solution().unwrap();
        assert!(!solution.moles().is_empty());
        assert_eq!(solution.moles().len(), 2);
    }

    // -----------------------------------------------------------------------
    // D.6 — run_keq_cross_validation() — indirect test via enabling keq validation
    // -----------------------------------------------------------------------
    #[test]
    fn solve_with_keq_validation_enabled_does_not_crash() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.n0 = vec![1.0, 0.0];
        solver.T = 3000.0;
        solver.P = 101_325.0;
        solver.p0 = 101_325.0;
        solver.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        solver.gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        solver.phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        solver.species_phase = vec![0, 0];
        solver.solver_settings.solver = Solvers::LM;
        solver.solver_settings.keq_validation_mode = EquilibriumConstantValidationMode::WhenApplicable;
        solver.create_stoich_matrix().unwrap();
        solver.solve().unwrap();
        let solution = solver.accepted_solution().unwrap();
        assert!(!solution.moles().is_empty());
    }

    // -----------------------------------------------------------------------
    // D.10 — build_multiphase_acceptance_report()
    // -----------------------------------------------------------------------
    #[test]
    fn build_multiphase_acceptance_report_rejects_dimension_mismatch() {
        let phase_set = PhaseSet::from_policy(
            &InitialPhaseSet::AllCandidatePhases,
            &[true],
        )
        .unwrap();
        let report = PhaseControlledSolveReport {
            iterations: 1,
            initial_active_phases: vec![PhaseIndex::new(0, 1).unwrap()],
            final_active_phases: vec![PhaseIndex::new(0, 1).unwrap()],
            initial_phase_set: phase_set.clone(),
            final_phase_set: phase_set,
            transitions: vec![],
            final_validation: EquilibriumCandidateReport {
                residual_l2_norm: 0.0,
                residual_rms: 0.0,
                max_abs_residual: 0.0,
                raw_residual_l2_norm: 0.0,
                raw_residual_rms: 0.0,
                raw_max_abs_residual: 0.0,
                max_abs_element_balance_error: 0.0,
                reaction_affinity_l2_norm: 0.0,
                max_abs_reaction_affinity: 0.0,
                min_moles: 1.0,
            },
            nonlinear_reports: vec![],
        };
        // 1 phase in control report, but 2 stability reports → mismatch
        let stability = vec![
            PhaseStabilityReport {
                phase: PhaseIndex::new(0, 2).unwrap(),
                model: crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStabilityModel::FixedIdealGas,
                active: true,
                driving_force: None,
                element_potentials: None,
            },
            PhaseStabilityReport {
                phase: PhaseIndex::new(1, 2).unwrap(),
                model: crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStabilityModel::Unsupported,
                active: false,
                driving_force: None,
                element_potentials: None,
            },
        ];
        let result = build_multiphase_acceptance_report(report, stability, None, None);
        assert!(result.is_err());
    }

    #[test]
    fn build_multiphase_acceptance_report_accepts_matching_dimensions() {
        let phase_set = PhaseSet::from_policy(
            &InitialPhaseSet::AllCandidatePhases,
            &[true],
        )
        .unwrap();
        let stability = vec![
            PhaseStabilityReport {
                phase: PhaseIndex::new(0, 1).unwrap(),
                model: crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStabilityModel::FixedIdealGas,
                active: true,
                driving_force: None,
                element_potentials: None,
            },
        ];
        let report = PhaseControlledSolveReport {
            iterations: 1,
            initial_active_phases: vec![PhaseIndex::new(0, 1).unwrap()],
            final_active_phases: vec![PhaseIndex::new(0, 1).unwrap()],
            initial_phase_set: phase_set.clone(),
            final_phase_set: phase_set,
            transitions: vec![],
            final_validation: EquilibriumCandidateReport {
                residual_l2_norm: 0.0,
                residual_rms: 0.0,
                max_abs_residual: 0.0,
                raw_residual_l2_norm: 0.0,
                raw_residual_rms: 0.0,
                raw_max_abs_residual: 0.0,
                max_abs_element_balance_error: 0.0,
                reaction_affinity_l2_norm: 0.0,
                max_abs_reaction_affinity: 0.0,
                min_moles: 1.0,
            },
            nonlinear_reports: vec![],
        };
        let result = build_multiphase_acceptance_report(report, stability, Some(0.0), Some(0.0));
        assert!(result.is_ok());
        let acceptance = result.unwrap();
        assert_eq!(acceptance.phase_stability.len(), 1);
    }

    // -----------------------------------------------------------------------
    // D.11 — PhaseManager::detect_phase_destruction()
    // -----------------------------------------------------------------------
    #[test]
    fn detect_phase_destruction_returns_empty_when_no_phase_below_threshold() {
        let manager = PhaseManager::new(1e-3, 0.0, 0.0);
        let n_phase = vec![1.0, 0.5, 2.0];
        let active = vec![true, true, true];
        let destroyed = manager.detect_phase_destruction(&n_phase, &active);
        assert!(destroyed.is_empty());
    }

    #[test]
    fn detect_phase_destruction_ignores_inactive_phases_below_threshold() {
        let manager = PhaseManager::new(1e-3, 0.0, 0.0);
        let n_phase = vec![1.0, 1e-6, 2.0];
        let active = vec![true, false, true];
        let destroyed = manager.detect_phase_destruction(&n_phase, &active);
        // phase 1 is below threshold but inactive → not detected
        assert!(destroyed.is_empty());
    }

    #[test]
    fn detect_phase_destruction_finds_active_phase_below_threshold() {
        let manager = PhaseManager::new(1e-3, 0.0, 0.0);
        let n_phase = vec![1.0, 1e-6, 2.0];
        let active = vec![true, true, true];
        let destroyed = manager.detect_phase_destruction(&n_phase, &active);
        assert_eq!(destroyed, vec![1]);
    }

    #[test]
    fn detect_phase_destruction_handles_empty_input() {
        let manager = PhaseManager::new(1e-3, 0.0, 0.0);
        let destroyed = manager.detect_phase_destruction(&[], &[]);
        assert!(destroyed.is_empty());
    }

    // -----------------------------------------------------------------------
    // D.12 — PhaseManager::classify_phases() (без температуры)
    // -----------------------------------------------------------------------
    #[test]
    fn classify_phases_returns_no_transition_when_all_phases_stable() {
        let manager = PhaseManager::new(1e-30, 0.0, 0.0);
        let phase_set = PhaseSet::from_policy(
            &InitialPhaseSet::AllCandidatePhases,
            &[true],
        )
        .unwrap();
        let phase_totals = vec![1.0];
        let stability = vec![PhaseStabilityReport {
            phase: PhaseIndex::new(0, 1).unwrap(),
            model: crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStabilityModel::FixedIdealGas,
            active: true,
            driving_force: None,
            element_potentials: None,
        }];
        let plan = manager.classify_phases(&phase_totals, &stability, &phase_set).unwrap();
        assert!(matches!(plan, PhaseTransitionPlan::NoTransition { .. }));
    }

    #[test]
    fn classify_phases_rejects_temperature_scaled_hysteresis() {
        let manager = PhaseManager {
            phase_eps: 1e-30,
            phase_hysteresis: crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseHysteresisPolicy::TemperatureScaled {
                create_rt_factor: -1e-6,
                keep_rt_factor: 1e-8,
            },
            max_phase_iterations: 16,
            initial_phase_set: InitialPhaseSet::default(),
        };
        let phase_set = PhaseSet::from_policy(
            &InitialPhaseSet::AllCandidatePhases,
            &[true],
        )
        .unwrap();
        let phase_totals = vec![1.0];
        let stability = vec![PhaseStabilityReport {
            phase: PhaseIndex::new(0, 1).unwrap(),
            model: crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStabilityModel::FixedIdealGas,
            active: true,
            driving_force: None,
            element_potentials: None,
        }];
        let result = manager.classify_phases(&phase_totals, &stability, &phase_set);
        assert!(result.is_err());
    }

    #[test]
    fn classify_phases_rejects_dimension_mismatch() {
        let manager = PhaseManager::new(1e-30, 0.0, 0.0);
        let phase_set = PhaseSet::from_policy(
            &InitialPhaseSet::AllCandidatePhases,
            &[true, false],
        )
        .unwrap();
        // 2 phases in set, but only 1 total and 1 stability report
        let phase_totals = vec![1.0];
        let stability = vec![PhaseStabilityReport {
            phase: PhaseIndex::new(0, 1).unwrap(),
            model: crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStabilityModel::FixedIdealGas,
            active: true,
            driving_force: None,
            element_potentials: None,
        }];
        let result = manager.classify_phases(&phase_totals, &stability, &phase_set);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.15 — seed_activated_phase() с разными PhaseSeedPolicy
    // -----------------------------------------------------------------------
    #[test]
    fn seed_activated_phase_trace_floor_sets_log_moles_to_floor() {
        let mut log_moles = vec![1.0_f64.ln(), 2.0_f64.ln(), 3.0_f64.ln()];
        let species_phase = vec![0, 0, 1];
        seed_activated_phase(
            &mut log_moles,
            PhaseIndex::new(1, 2).unwrap(),
            &species_phase,
            PhaseSeedPolicy::TraceFloor,
        )
        .unwrap();
        // species 2 (phase 1) should be set to trace floor
        let expected_ln = PHASE_CONTROL_TRACE_MOLE_FLOOR.max(PHASE_CONTROL_TRACE_MOLE_FLOOR).ln();
        assert!((log_moles[2] - expected_ln).abs() < 1e-10);
        // species 0,1 (phase 0) should be unchanged
        assert!((log_moles[0] - 1.0_f64.ln()).abs() < 1e-10);
        assert!((log_moles[1] - 2.0_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn seed_activated_phase_absolute_per_species_sets_positive_moles() {
        let mut log_moles = vec![1.0_f64.ln(), 2.0_f64.ln()];
        let species_phase = vec![0, 0];
        seed_activated_phase(
            &mut log_moles,
            PhaseIndex::new(0, 1).unwrap(),
            &species_phase,
            PhaseSeedPolicy::AbsolutePerSpecies { moles: 0.1 },
        )
        .unwrap();
        let expected_ln = 0.1_f64.ln();
        assert!((log_moles[0] - expected_ln).abs() < 1e-10);
        assert!((log_moles[1] - expected_ln).abs() < 1e-10);
    }

    #[test]
    fn seed_activated_phase_relative_to_system_total_scales_with_inventory() {
        let mut log_moles = vec![1.0_f64.ln(), 0.0_f64.ln()];
        let species_phase = vec![0, 1];
        // total = 1.0 + 0.0 = 1.0, fraction = 0.1 → seed = 0.1
        seed_activated_phase(
            &mut log_moles,
            PhaseIndex::new(1, 2).unwrap(),
            &species_phase,
            PhaseSeedPolicy::RelativeToSystemTotal {
                fraction: 0.1,
                minimum: 1e-10,
            },
        )
        .unwrap();
        let expected_ln = 0.1_f64.ln();
        assert!((log_moles[1] - expected_ln).abs() < 1e-10);
    }

    #[test]
    fn seed_activated_phase_rejects_phase_with_no_species() {
        let mut log_moles = vec![0.0_f64.ln()];
        let species_phase = vec![0];
        // phase 1 has no species in species_phase
        let result = seed_activated_phase(
            &mut log_moles,
            PhaseIndex::new(1, 2).unwrap(),
            &species_phase,
            PhaseSeedPolicy::TraceFloor,
        );
        assert!(result.is_err());
    }

    #[test]
    fn seed_activated_phase_rejects_dimension_mismatch() {
        let mut log_moles = vec![0.0_f64.ln()];
        let species_phase = vec![0, 1]; // longer than log_moles
        let result = seed_activated_phase(
            &mut log_moles,
            PhaseIndex::new(0, 2).unwrap(),
            &species_phase,
            PhaseSeedPolicy::TraceFloor,
        );
        assert!(result.is_err());
    }

    #[test]
    fn seed_activated_phase_absolute_rejects_non_positive_moles() {
        let mut log_moles = vec![0.0_f64.ln()];
        let species_phase = vec![0];
        let result = seed_activated_phase(
            &mut log_moles,
            PhaseIndex::new(0, 1).unwrap(),
            &species_phase,
            PhaseSeedPolicy::AbsolutePerSpecies { moles: -1.0 },
        );
        assert!(result.is_err());
    }

    #[test]
    fn seed_activated_phase_relative_rejects_invalid_fraction() {
        let mut log_moles = vec![0.0_f64.ln()];
        let species_phase = vec![0];
        let result = seed_activated_phase(
            &mut log_moles,
            PhaseIndex::new(0, 1).unwrap(),
            &species_phase,
            PhaseSeedPolicy::RelativeToSystemTotal {
                fraction: 1.5, // > 1.0
                minimum: 1e-10,
            },
        );
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.16 — deactivate_phases_seed_only() с несколькими фазами
    // -----------------------------------------------------------------------
    #[test]
    fn deactivate_phases_seed_only_sets_species_to_trace_floor() {
        let mut y = vec![1.0_f64.ln(), 2.0_f64.ln(), 3.0_f64.ln()];
        let species_phase = vec![0, 1, 1];
        deactivate_phases_seed_only(
            &mut y,
            &[PhaseIndex::new(1, 2).unwrap()],
            &species_phase,
            1e-30,
        )
        .unwrap();
        let floor_ln = 1e-30_f64.ln();
        // species 0 (phase 0) unchanged
        assert!((y[0] - 1.0_f64.ln()).abs() < 1e-10);
        // species 1,2 (phase 1) set to floor
        assert!((y[1] - floor_ln).abs() < 1e-10);
        assert!((y[2] - floor_ln).abs() < 1e-10);
    }

    #[test]
    fn deactivate_phases_seed_only_handles_multiple_phases() {
        let mut y = vec![1.0_f64.ln(), 2.0_f64.ln(), 3.0_f64.ln(), 4.0_f64.ln()];
        let species_phase = vec![0, 0, 1, 2];
        deactivate_phases_seed_only(
            &mut y,
            &[
                PhaseIndex::new(0, 3).unwrap(),
                PhaseIndex::new(2, 3).unwrap(),
            ],
            &species_phase,
            1e-30,
        )
        .unwrap();
        let floor_ln = 1e-30_f64.ln();
        // phase 0 species (0,1) → floor
        assert!((y[0] - floor_ln).abs() < 1e-10);
        assert!((y[1] - floor_ln).abs() < 1e-10);
        // phase 1 species (2) → unchanged
        assert!((y[2] - 3.0_f64.ln()).abs() < 1e-10);
        // phase 2 species (3) → floor
        assert!((y[3] - floor_ln).abs() < 1e-10);
    }

    #[test]
    fn deactivate_phases_seed_only_rejects_non_positive_trace_floor() {
        let mut y = vec![1.0_f64.ln()];
        let species_phase = vec![0];
        let result = deactivate_phases_seed_only(
            &mut y,
            &[PhaseIndex::new(0, 1).unwrap()],
            &species_phase,
            -1.0,
        );
        assert!(result.is_err());
    }

    #[test]
    fn deactivate_phases_seed_only_rejects_non_finite_trace_floor() {
        let mut y = vec![1.0_f64.ln()];
        let species_phase = vec![0];
        let result = deactivate_phases_seed_only(
            &mut y,
            &[PhaseIndex::new(0, 1).unwrap()],
            &species_phase,
            f64::NAN,
        );
        assert!(result.is_err());
    }

    #[test]
    fn deactivate_phases_seed_only_handles_empty_deactivation_list() {
        let mut y = vec![1.0_f64.ln(), 2.0_f64.ln()];
        let species_phase = vec![0, 1];
        deactivate_phases_seed_only(&mut y, &[], &species_phase, 1e-30).unwrap();
        // all values unchanged
        assert!((y[0] - 1.0_f64.ln()).abs() < 1e-10);
        assert!((y[1] - 2.0_f64.ln()).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // D.21 — TemperatureSolveSnapshot and TemperatureSolveFailure
    // -----------------------------------------------------------------------
    #[test]
    fn temperature_solve_snapshot_holds_expected_fields() {
        let validation = EquilibriumCandidateReport {
            residual_l2_norm: 1e-8,
            residual_rms: 1e-9,
            max_abs_residual: 1e-7,
            raw_residual_l2_norm: 1e-6,
            raw_residual_rms: 1e-7,
            raw_max_abs_residual: 1e-5,
            max_abs_element_balance_error: 1e-10,
            reaction_affinity_l2_norm: 1e-8,
            max_abs_reaction_affinity: 1e-7,
            min_moles: 1e-3,
        };
        let solve_report = EquilibriumSolveReport {
            policy: SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM)),
            attempts: vec![SolverAttemptReport {
                backend: SolverBackend::Legacy(Solvers::LM),
                outcome: SolverAttemptOutcome::Accepted,
                metrics: None,
            }],
            accepted_backend: SolverBackend::Legacy(Solvers::LM),
        };
        let snapshot = TemperatureSolveSnapshot {
            temperature: 3000.0,
            log_moles: vec![-1.0, -2.0],
            moles: vec![0.5, 0.3],
            validation,
            solve_report,
        };
        assert!((snapshot.temperature - 3000.0).abs() < 1e-12);
        assert_eq!(snapshot.log_moles.len(), 2);
        assert_eq!(snapshot.moles.len(), 2);
        assert!((snapshot.moles[0] - 0.5).abs() < 1e-12);
        assert!(snapshot.validation.min_moles > 0.0);
        assert!(snapshot.solve_report.accepted_attempt().is_some());
    }

    #[test]
    fn temperature_solve_failure_holds_error_message() {
        let failure = TemperatureSolveFailure {
            temperature: 5000.0,
            message: "solver did not converge".to_string(),
            attempts: vec![],
        };
        assert!((failure.temperature - 5000.0).abs() < 1e-12);
        assert_eq!(failure.message, "solver did not converge");
        assert!(failure.attempts.is_empty());
    }

    #[test]
    fn temperature_solve_failure_carries_backend_attempts() {
        let failure = TemperatureSolveFailure {
            temperature: 6000.0,
            message: "all backends failed".to_string(),
            attempts: vec![
                SolverAttemptReport {
                    backend: SolverBackend::Legacy(Solvers::LM),
                    outcome: SolverAttemptOutcome::Failed {
                        kind: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptFailureKind::Solver,
                        reason: "max iterations".to_string(),
                    },
                    metrics: None,
                },
                SolverAttemptReport {
                    backend: SolverBackend::Legacy(Solvers::NR),
                    outcome: SolverAttemptOutcome::Failed {
                        kind: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptFailureKind::Solver,
                        reason: "singular".to_string(),
                    },
                    metrics: None,
                },
            ],
        };
        assert_eq!(failure.attempts.len(), 2);
        assert!(failure.attempts[0].is_started());
        assert!(!failure.attempts[0].is_accepted());
    }

    // -----------------------------------------------------------------------
    // D.22 — TemperatureWorkerSeed::apply()
    // -----------------------------------------------------------------------
    #[test]
    fn temperature_worker_seed_apply_copies_state_to_local_solver() {
        let mut source = EquilibriumLogMoles::empty();
        source.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        source.n0 = vec![1.0, 0.0];
        source.T = 3000.0;
        source.P = 101_325.0;
        source.p0 = 101_325.0;
        source.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        source.gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        source.phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        source.species_phase = vec![0, 0];
        source.solver_settings.solver = Solvers::LM;
        source.create_stoich_matrix().unwrap();

        let seed = source.temperature_worker_seed();
        let mut local = EquilibriumLogMoles::empty();
        let new_gibbs = vec![
            Rc::new(|_| 100.0) as GibbsFn,
            Rc::new(|_| 200.0) as GibbsFn,
        ];
        seed.apply(&mut local, 3500.0, new_gibbs, None);

        assert!((local.T - 3500.0).abs() < 1e-12);
        assert_eq!(local.subs_data.substances.len(), 2);
        assert_eq!(local.n0, vec![1.0, 0.0]);
        assert_eq!(local.solution.len(), 0); // cleared
        assert!(local.last_validation_report.is_none());
        assert!(local.last_solve_report.is_none());
    }

    #[test]
    fn temperature_worker_seed_apply_clears_previous_solution() {
        let mut source = EquilibriumLogMoles::empty();
        source.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        source.n0 = vec![1.0, 0.0];
        source.T = 3000.0;
        source.P = 101_325.0;
        source.p0 = 101_325.0;
        source.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        source.gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        source.phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        source.species_phase = vec![0, 0];
        source.solver_settings.solver = Solvers::LM;
        source.create_stoich_matrix().unwrap();

        let seed = source.temperature_worker_seed();
        let mut local = EquilibriumLogMoles::empty();
        // Pre-populate local with stale data
        local.solution = vec![1.0, 2.0];
        local.moles = vec![3.0, 4.0];
        local.last_validation_report = Some(EquilibriumCandidateReport {
            residual_l2_norm: 0.0,
            residual_rms: 0.0,
            max_abs_residual: 0.0,
            raw_residual_l2_norm: 0.0,
            raw_residual_rms: 0.0,
            raw_max_abs_residual: 0.0,
            max_abs_element_balance_error: 0.0,
            reaction_affinity_l2_norm: 0.0,
            max_abs_reaction_affinity: 0.0,
            min_moles: 1.0,
        });

        let new_gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        seed.apply(&mut local, 3000.0, new_gibbs, None);

        // Stale data must be cleared
        assert!(local.solution.is_empty());
        assert!(local.moles.is_empty());
        assert!(local.last_validation_report.is_none());
        assert!(local.last_solve_report.is_none());
    }

    // -----------------------------------------------------------------------
    // D.38 — resolved_initial_guess()
    // -----------------------------------------------------------------------
    #[test]
    fn resolved_initial_guess_derives_seed_from_n0_when_no_explicit_guess() {
        let solver = EquilibriumLogMoles::empty();
        // empty solver has n0 = [] → should succeed with empty seed
        let result = solver.resolved_initial_guess();
        assert!(result.is_ok());
        let (guess, derived) = result.unwrap();
        assert!(derived); // was derived from n0
        assert!(guess.is_empty());
    }

    #[test]
    fn resolved_initial_guess_uses_explicit_seed_when_provided() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.n0 = vec![1.0, 0.0];
        solver.initial_guess = Some(vec![-1.0, -5.0]);
        let result = solver.resolved_initial_guess();
        assert!(result.is_ok());
        let (guess, derived) = result.unwrap();
        assert!(!derived); // was NOT derived from n0
        assert!((guess[0] - (-1.0)).abs() < 1e-12);
        assert!((guess[1] - (-5.0)).abs() < 1e-12);
    }

    #[test]
    fn resolved_initial_guess_rejects_dimension_mismatch() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        solver.n0 = vec![1.0, 0.0];
        // 1 entry for 2 substances → mismatch
        solver.initial_guess = Some(vec![-1.0]);
        let result = solver.resolved_initial_guess();
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.40 — temperature_worker_seed()
    // -----------------------------------------------------------------------
    #[test]
    fn temperature_worker_seed_apply_preserves_solver_state() {
        // Verify that temperature_worker_seed() captures state by applying it
        // to a local solver and checking the result via public API.
        let mut source = EquilibriumLogMoles::empty();
        source.subs_data.substances = vec!["O2".to_string(), "O".to_string()];
        source.n0 = vec![1.0, 0.0];
        source.T = 3000.0;
        source.P = 101_325.0;
        source.p0 = 101_325.0;
        source.elem_composition = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        source.gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        source.phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        source.species_phase = vec![0, 0];
        source.solver_settings.solver = Solvers::LM;
        source.create_stoich_matrix().unwrap();

        let seed = source.temperature_worker_seed();
        let mut local = EquilibriumLogMoles::empty();
        let new_gibbs = vec![
            Rc::new(|_| 0.0) as GibbsFn,
            Rc::new(|_| 0.0) as GibbsFn,
        ];
        seed.apply(&mut local, 3000.0, new_gibbs, None);

        // After apply, local should have the captured state
        assert_eq!(local.subs_data.substances.len(), 2);
        assert_eq!(local.n0, vec![1.0, 0.0]);
        assert!((local.P - 101_325.0).abs() < 1e-6);
        assert!((local.p0 - 101_325.0).abs() < 1e-6);
        assert_eq!(local.phases.len(), 1);
        assert_eq!(local.species_phase, vec![0, 0]);
    }

    // -----------------------------------------------------------------------
    // D.41 — publish_solve_candidate()
    // -----------------------------------------------------------------------
    #[test]
    fn publish_solve_candidate_updates_all_publication_fields() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];

        let candidate = EquilibriumSolveCandidate {
            log_moles: vec![-1.0, -2.0],
            moles: vec![0.5, 0.3],
            mole_table: {
                let mut map = HashMap::new();
                map.insert("O2".to_string(), vec![0.5]);
                map.insert("O".to_string(), vec![0.3]);
                map
            },
            validation_report: EquilibriumCandidateReport {
                residual_l2_norm: 1e-8,
                residual_rms: 1e-9,
                max_abs_residual: 1e-7,
                raw_residual_l2_norm: 1e-6,
                raw_residual_rms: 1e-7,
                raw_max_abs_residual: 1e-5,
                max_abs_element_balance_error: 1e-10,
                reaction_affinity_l2_norm: 1e-8,
                max_abs_reaction_affinity: 1e-7,
                min_moles: 0.3,
            },
            solve_report: EquilibriumSolveReport {
                policy: SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM)),
                attempts: vec![SolverAttemptReport {
                    backend: SolverBackend::Legacy(Solvers::LM),
                    outcome: SolverAttemptOutcome::Accepted,
                    metrics: None,
                }],
                accepted_backend: SolverBackend::Legacy(Solvers::LM),
            },
            keq_validation_status: Some(
                EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable {
                    message: "single reaction system".to_string(),
                },
            ),
        };

        solver.publish_solve_candidate(candidate);

        assert_eq!(solver.solution, vec![-1.0, -2.0]);
        assert_eq!(solver.moles, vec![0.5, 0.3]);
        assert!(solver.last_validation_report.is_some());
        assert!(solver.last_solve_report.is_some());
        assert!(solver.last_keq_validation_status.is_some());
        assert_eq!(solver.map_of_moles_for_each_substance.len(), 2);
    }

    #[test]
    fn publish_solve_candidate_accepts_missing_keq_validation() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.subs_data.substances = vec!["O2".to_string(), "O".to_string()];

        let candidate = EquilibriumSolveCandidate {
            log_moles: vec![-1.0, -2.0],
            moles: vec![0.5, 0.3],
            mole_table: HashMap::new(),
            validation_report: EquilibriumCandidateReport {
                residual_l2_norm: 0.0,
                residual_rms: 0.0,
                max_abs_residual: 0.0,
                raw_residual_l2_norm: 0.0,
                raw_residual_rms: 0.0,
                raw_max_abs_residual: 0.0,
                max_abs_element_balance_error: 0.0,
                reaction_affinity_l2_norm: 0.0,
                max_abs_reaction_affinity: 0.0,
                min_moles: 0.3,
            },
            solve_report: EquilibriumSolveReport {
                policy: SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM)),
                attempts: vec![SolverAttemptReport {
                    backend: SolverBackend::Legacy(Solvers::LM),
                    outcome: SolverAttemptOutcome::Accepted,
                    metrics: None,
                }],
                accepted_backend: SolverBackend::Legacy(Solvers::LM),
            },
            keq_validation_status: None,
        };

        solver.publish_solve_candidate(candidate);
        assert!(solver.last_keq_validation_status.is_none());
    }

    // -----------------------------------------------------------------------
    // D.53 — scaled_residual() and scaled_jacobian()
    // -----------------------------------------------------------------------
    #[test]
    fn scaled_residual_applies_row_scaling() {
        let f = Box::new(|_xi: &[f64]| Ok(vec![10.0, 20.0, 30.0]));
        let scale = vec![2.0, 5.0, 10.0];
        let scaled = scaled_residual(f, scale);
        let result = scaled(&[1.0, 2.0, 3.0]).unwrap();
        assert!((result[0] - 5.0).abs() < 1e-12);  // 10/2
        assert!((result[1] - 4.0).abs() < 1e-12);  // 20/5
        assert!((result[2] - 3.0).abs() < 1e-12);  // 30/10
    }

    #[test]
    fn scaled_jacobian_applies_row_scaling() {
        let j = Box::new(|_xi: &[f64]| Ok(DMatrix::from_row_slice(2, 2, &[2.0, 4.0, 6.0, 8.0])));
        let scale = vec![2.0, 4.0];
        let scaled = scaled_jacobian(j, scale);
        let result = scaled(&[1.0, 2.0]).unwrap();
        assert!((result[(0, 0)] - 1.0).abs() < 1e-12); // 2/2
        assert!((result[(0, 1)] - 2.0).abs() < 1e-12); // 4/2
        assert!((result[(1, 0)] - 1.5).abs() < 1e-12); // 6/4
        assert!((result[(1, 1)] - 2.0).abs() < 1e-12); // 8/4
    }

    // -----------------------------------------------------------------------
    // D.54 — scale_residual_rows() and scale_jacobian_rows()
    // -----------------------------------------------------------------------
    #[test]
    fn scale_residual_rows_rejects_dimension_mismatch() {
        let result = scale_residual_rows(vec![1.0, 2.0], &[1.0]);
        assert!(result.is_err());
    }

    #[test]
    fn scale_residual_rows_rejects_non_positive_scale() {
        let result = scale_residual_rows(vec![1.0], &[0.0]);
        assert!(result.is_err());
        let result = scale_residual_rows(vec![1.0], &[-1.0]);
        assert!(result.is_err());
    }

    #[test]
    fn scale_residual_rows_rejects_non_finite_scale() {
        let result = scale_residual_rows(vec![1.0], &[f64::NAN]);
        assert!(result.is_err());
        let result = scale_residual_rows(vec![1.0], &[f64::INFINITY]);
        assert!(result.is_err());
    }

    #[test]
    fn scale_residual_rows_applies_correct_division() {
        let result = scale_residual_rows(vec![10.0, 20.0, 30.0], &[2.0, 4.0, 5.0]).unwrap();
        assert!((result[0] - 5.0).abs() < 1e-12);
        assert!((result[1] - 5.0).abs() < 1e-12);
        assert!((result[2] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn scale_jacobian_rows_rejects_dimension_mismatch() {
        let j = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let result = scale_jacobian_rows(j, &[1.0]);
        assert!(result.is_err());
    }

    #[test]
    fn scale_jacobian_rows_rejects_non_positive_scale() {
        let j = DMatrix::from_row_slice(1, 2, &[1.0, 2.0]);
        let result = scale_jacobian_rows(j, &[0.0]);
        assert!(result.is_err());
    }

    #[test]
    fn scale_jacobian_rows_applies_correct_division() {
        let j = DMatrix::from_row_slice(2, 3, &[2.0, 4.0, 6.0, 8.0, 10.0, 12.0]);
        let result = scale_jacobian_rows(j, &[2.0, 4.0]).unwrap();
        assert!((result[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((result[(0, 1)] - 2.0).abs() < 1e-12);
        assert!((result[(0, 2)] - 3.0).abs() < 1e-12);
        assert!((result[(1, 0)] - 2.0).abs() < 1e-12);
        assert!((result[(1, 1)] - 2.5).abs() < 1e-12);
        assert!((result[(1, 2)] - 3.0).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // D.55 — evaluate_equilibrium_logmole_residual()
    // -----------------------------------------------------------------------
    #[test]
    fn evaluate_equilibrium_logmole_residual_rejects_dimension_mismatch() {
        // 2 species, 1 reaction, 1 element
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]); // O2 ↔ 2O
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let element_totals = vec![2.0];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        let species_phase = vec![0, 0];
        let phase_stoich = vec![vec![-1.0]]; // 1 reaction, 1 phase

        // Wrong log_moles length
        let result = evaluate_equilibrium_logmole_residual(
            &[1.0], // 1 entry, expected 2
            &reactions,
            &elements,
            &element_totals,
            &gibbs,
            &phases,
            3000.0,
            101325.0,
            101325.0,
            &species_phase,
            &phase_stoich,
        );
        assert!(result.is_err());
    }

    #[test]
    fn evaluate_equilibrium_logmole_residual_rejects_invalid_temperature() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let element_totals = vec![2.0];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        let species_phase = vec![0, 0];
        let phase_stoich = vec![vec![-1.0]];

        let result = evaluate_equilibrium_logmole_residual(
            &[1.0, 2.0],
            &reactions,
            &elements,
            &element_totals,
            &gibbs,
            &phases,
            -1.0, // invalid temperature
            101325.0,
            101325.0,
            &species_phase,
            &phase_stoich,
        );
        assert!(result.is_err());
    }

    #[test]
    fn evaluate_equilibrium_logmole_residual_rejects_non_finite_log_moles() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let element_totals = vec![2.0];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        let species_phase = vec![0, 0];
        let phase_stoich = vec![vec![-1.0]];

        let result = evaluate_equilibrium_logmole_residual(
            &[f64::NAN, 2.0],
            &reactions,
            &elements,
            &element_totals,
            &gibbs,
            &phases,
            3000.0,
            101325.0,
            101325.0,
            &species_phase,
            &phase_stoich,
        );
        assert!(result.is_err());
    }

    #[test]
    fn evaluate_equilibrium_logmole_residual_returns_correct_length() {
        // 2 species, 1 reaction, 1 element → residual length = 1 + 1 = 2
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let element_totals = vec![2.0];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];
        let species_phase = vec![0, 0];
        let phase_stoich = vec![vec![-1.0]];

        let result = evaluate_equilibrium_logmole_residual(
            &[1.0_f64.ln(), 1.0_f64.ln()],
            &reactions,
            &elements,
            &element_totals,
            &gibbs,
            &phases,
            3000.0,
            101325.0,
            101325.0,
            &species_phase,
            &phase_stoich,
        );
        assert!(result.is_ok());
        let residual = result.unwrap();
        assert_eq!(residual.len(), 2); // r + e = 1 + 1
    }

    // -----------------------------------------------------------------------
    // D.56 — evaluate_equilibrium_logmole_jacobian()
    // -----------------------------------------------------------------------
    #[test]
    fn evaluate_equilibrium_logmole_jacobian_rejects_dimension_mismatch() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let species_phase = vec![0, 0];
        let phase_stoich = vec![vec![-1.0]];

        // Wrong log_moles length
        let result = evaluate_equilibrium_logmole_jacobian(
            &[1.0], // 1 entry, expected 2
            &reactions,
            &elements,
            &species_phase,
            &phase_stoich,
            1,
        );
        assert!(result.is_err());
    }

    #[test]
    fn evaluate_equilibrium_logmole_jacobian_returns_correct_shape() {
        // 2 species, 1 reaction, 1 element → Jacobian shape = 2×2
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let species_phase = vec![0, 0];
        let phase_stoich = vec![vec![-1.0]];

        let result = evaluate_equilibrium_logmole_jacobian(
            &[1.0_f64.ln(), 1.0_f64.ln()],
            &reactions,
            &elements,
            &species_phase,
            &phase_stoich,
            1,
        );
        assert!(result.is_ok());
        let jacobian = result.unwrap();
        assert_eq!(jacobian.nrows(), 2);
        assert_eq!(jacobian.ncols(), 2);
    }

    #[test]
    fn evaluate_equilibrium_logmole_jacobian_rejects_non_finite_log_moles() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let species_phase = vec![0, 0];
        let phase_stoich = vec![vec![-1.0]];

        let result = evaluate_equilibrium_logmole_jacobian(
            &[f64::NAN, 2.0],
            &reactions,
            &elements,
            &species_phase,
            &phase_stoich,
            1,
        );
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.57 — equilibrium_logmole_residual2() and equilibrium_logmole_jacobian2()
    // -----------------------------------------------------------------------
    #[test]
    fn equilibrium_logmole_residual2_creates_closure() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let element_totals = vec![2.0];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];

        let result = equilibrium_logmole_residual2(
            reactions,
            elements,
            element_totals,
            gibbs,
            phases,
            3000.0,
            101325.0,
            101325.0,
            1e-30,
            1e-30,
        );
        assert!(result.is_ok());
        let closure = result.unwrap();
        let residual = closure(&[1.0_f64.ln(), 1.0_f64.ln()]).unwrap();
        assert_eq!(residual.len(), 2);
    }

    #[test]
    fn equilibrium_logmole_residual2_rejects_wrong_log_moles_length() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let element_totals = vec![2.0];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let phases = vec![Phase { kind: PhaseKind::IdealGas, species: vec![0, 1] }];

        let closure = equilibrium_logmole_residual2(
            reactions,
            elements,
            element_totals,
            gibbs,
            phases,
            3000.0,
            101325.0,
            101325.0,
            1e-30,
            1e-30,
        )
        .unwrap();
        let result = closure(&[1.0_f64.ln()]); // 1 entry, expected 2
        assert!(result.is_err());
    }

    #[test]
    fn equilibrium_logmole_jacobian2_returns_correct_shape() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[2.0, 1.0]);
        let species_phase = vec![0, 0];
        let phase_stoich = vec![vec![-1.0]];

        let result = equilibrium_logmole_jacobian2(
            &[1.0_f64.ln(), 1.0_f64.ln()],
            &reactions,
            &elements,
            &species_phase,
            &phase_stoich,
            1,
            1e-30,
            1e-30,
        );
        assert!(result.is_ok());
        let jacobian = result.unwrap();
        assert_eq!(jacobian.nrows(), 2);
        assert_eq!(jacobian.ncols(), 2);
    }

    // Note: equilibrium_logmole_jacobian2 (deprecated) does not call
    // validate_logmole_system_dimensions, so dimension mismatches are not
    // caught as typed errors. This is a known limitation of the deprecated path.

    // -----------------------------------------------------------------------
    // D.32 — validate_temperature_range() — indirect via solve_for_T_range()
    // -----------------------------------------------------------------------
    #[test]
    fn solve_for_T_range_rejects_non_finite_T_start() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.solve_for_T_range(f64::NAN, 3000.0, 100.0);
        assert!(result.is_err());
    }

    #[test]
    fn solve_for_T_range_rejects_non_finite_T_end() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.solve_for_T_range(2000.0, f64::INFINITY, 100.0);
        assert!(result.is_err());
    }

    #[test]
    fn solve_for_T_range_rejects_non_positive_T_step() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.solve_for_T_range(2000.0, 3000.0, 0.0);
        assert!(result.is_err());
        let result = solver.solve_for_T_range(2000.0, 3000.0, -10.0);
        assert!(result.is_err());
    }

    #[test]
    fn solve_for_T_range_rejects_T_start_gte_T_end() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.solve_for_T_range(3000.0, 2000.0, 100.0);
        assert!(result.is_err());
        let result = solver.solve_for_T_range(3000.0, 3000.0, 100.0);
        assert!(result.is_err());
    }

    #[test]
    fn solve_for_T_range_rejects_non_finite_T_step() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.solve_for_T_range(2000.0, 3000.0, f64::NAN);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.23 — continuation_seed_for_point() — indirect via solve_for_T_range
    // -----------------------------------------------------------------------

    #[test]
    fn solve_for_T_range_par_rejects_invalid_temperature_range() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.solve_for_T_range_par(f64::NAN, 3000.0, 100.0);
        assert!(result.is_err());
    }

    #[test]
    fn solve_for_T_range_par2_rejects_invalid_temperature_range() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.solve_for_T_range_par2(2000.0, 2000.0, 100.0);
        assert!(result.is_err());
    }

    // ===================================================================
    // Direct tests for newly exposed pub(crate) methods
    // ===================================================================

    // -----------------------------------------------------------------------
    // D.23 — continuation_seed_for_point() — direct test
    // -----------------------------------------------------------------------
    #[test]
    fn continuation_seed_for_point_uses_previous_accepted_when_policy_says_so() {
        let seed = continuation_seed_for_point(
            ContinuationSeedPolicy::PreviousAccepted,
            &Some(vec![1.0, 2.0]),
            &Some(vec![3.0, 4.0]),
        );
        assert_eq!(seed, Some(vec![3.0, 4.0]));
    }

    #[test]
    fn continuation_seed_for_point_uses_configured_seed_when_independent() {
        let seed = continuation_seed_for_point(
            ContinuationSeedPolicy::IndependentPerPoint,
            &Some(vec![1.0, 2.0]),
            &Some(vec![3.0, 4.0]),
        );
        assert_eq!(seed, Some(vec![1.0, 2.0]));
    }

    #[test]
    fn continuation_seed_for_point_returns_none_when_no_seed_available() {
        let seed = continuation_seed_for_point(
            ContinuationSeedPolicy::PreviousAccepted,
            &None,
            &None,
        );
        assert!(seed.is_none());
    }

    // -----------------------------------------------------------------------
    // D.24 — build_temperature_ranges() and build_temperature_gibbs_cache()
    // -----------------------------------------------------------------------
    #[test]
    fn build_temperature_ranges_returns_single_range_on_empty_solver() {
        let mut solver = EquilibriumLogMoles::empty();
        let ranges = solver.build_temperature_ranges(300.0, 400.0, 50.0);
        assert!(ranges.is_ok());
        // Empty solver has no substances, so the while loop body executes once
        // and produces a single range covering the full interval
        let ranges = ranges.unwrap();
        assert_eq!(ranges.len(), 1);
        assert_eq!(ranges[0], (300.0, 400.0));
    }

    #[test]
    fn build_temperature_gibbs_cache_returns_empty_on_empty_solver() {
        let mut solver = EquilibriumLogMoles::empty();
        let cache = solver.build_temperature_gibbs_cache(&[]);
        assert!(cache.is_ok());
        assert!(cache.unwrap().is_empty());
    }

    #[test]
    fn build_temperature_point_index_returns_empty_for_empty_ranges() {
        let pairs = EquilibriumLogMoles::build_temperature_point_index(&[], 300.0, 400.0, 50.0);
        assert!(pairs.is_empty());
    }

    #[test]
    fn build_temperature_point_index_maps_temperatures_to_range_ids() {
        let ranges = vec![(300.0, 400.0), (400.0, 500.0)];
        let pairs = EquilibriumLogMoles::build_temperature_point_index(&ranges, 300.0, 500.0, 50.0);
        // 300→0, 350→0, 400→1, 450→1
        assert_eq!(pairs.len(), 4);
        assert_eq!(pairs[0], (300.0, 0));
        assert_eq!(pairs[1], (350.0, 0));
        assert_eq!(pairs[2], (400.0, 1));
        assert_eq!(pairs[3], (450.0, 1));
    }

    // -----------------------------------------------------------------------
    // D.25 — solve_temperature_point_from_seed()
    // -----------------------------------------------------------------------
    #[test]
    fn solve_temperature_point_from_seed_fails_on_empty_seed() {
        let solver = EquilibriumLogMoles::empty();
        let seed = solver.temperature_worker_seed();
        let result = EquilibriumLogMoles::solve_temperature_point_from_seed(
            &seed,
            300.0,
            vec![],
            None,
            0.0,
        );
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.26 — snapshot_published_solution()
    // -----------------------------------------------------------------------
    #[test]
    fn snapshot_published_solution_fails_without_publication() {
        let solver = EquilibriumLogMoles::empty();
        let result = solver.snapshot_published_solution(300.0);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.27 — ordered_gibbs_functions() and build_gibbs_functions()
    // -----------------------------------------------------------------------
    #[test]
    fn ordered_gibbs_functions_returns_error_for_missing_substance() {
        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut functions: HashMap<String, Box<dyn Fn(f64) -> f64 + Send + Sync>> =
            HashMap::new();
        functions.insert("O2".to_string(), Box::new(|_| 0.0));
        // "O" is missing
        let result = EquilibriumLogMoles::ordered_gibbs_functions(&substances, functions);
        assert!(result.is_err());
    }

    #[test]
    fn ordered_gibbs_functions_returns_empty_for_empty_input() {
        let result = EquilibriumLogMoles::ordered_gibbs_functions(&[], HashMap::new());
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    #[test]
    fn build_gibbs_functions_returns_empty_on_empty_solver() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.build_gibbs_functions();
        // Empty solver has no substances → empty Gibbs vec
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    // -----------------------------------------------------------------------
    // D.28 — build_parallel_gibbs_functions()
    // -----------------------------------------------------------------------
    #[test]
    fn build_parallel_gibbs_functions_returns_empty_on_empty_solver() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.build_parallel_gibbs_functions();
        // Empty solver has no substances → empty Gibbs vec
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    // -----------------------------------------------------------------------
    // D.29 — reconstructed_mole_state()
    // -----------------------------------------------------------------------
    #[test]
    fn reconstructed_mole_state_rejects_dimension_mismatch() {
        let solver = EquilibriumLogMoles::empty();
        // solver has 0 substances, but we pass a non-empty solution
        let result = solver.reconstructed_mole_state(&[1.0, 2.0]);
        assert!(result.is_err());
    }

    #[test]
    fn reconstructed_mole_state_rejects_non_finite_log_moles() {
        let solver = EquilibriumLogMoles::empty();
        let result = solver.reconstructed_mole_state(&[f64::NAN]);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.30 — check_task()
    // -----------------------------------------------------------------------
    #[test]
    fn check_task_fails_on_empty_solver_without_settings() {
        let solver = EquilibriumLogMoles::empty();
        // Default settings are valid, but n0.len() != substances.len() (0 == 0)
        // Actually empty solver has 0 substances and 0 n0, so it should pass
        let result = solver.check_task();
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // D.32 — validate_temperature_range() — direct test
    // -----------------------------------------------------------------------
    #[test]
    fn validate_temperature_range_accepts_valid_range() {
        let result = EquilibriumLogMoles::validate_temperature_range(300.0, 3000.0, 100.0);
        assert!(result.is_ok());
    }

    #[test]
    fn validate_temperature_range_rejects_non_finite_start() {
        let result = EquilibriumLogMoles::validate_temperature_range(f64::NAN, 3000.0, 100.0);
        assert!(result.is_err());
    }

    #[test]
    fn validate_temperature_range_rejects_non_positive_step() {
        let result = EquilibriumLogMoles::validate_temperature_range(300.0, 3000.0, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn validate_temperature_range_rejects_start_gte_end() {
        let result = EquilibriumLogMoles::validate_temperature_range(3000.0, 300.0, 100.0);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.35 — clear_published_state() — direct test
    // -----------------------------------------------------------------------
    #[test]
    fn clear_published_state_clears_all_fields() {
        let mut solver = EquilibriumLogMoles::empty();
        // Set some publication fields
        solver.solution = vec![1.0, 2.0];
        solver.moles = vec![3.0, 4.0];
        solver.clear_published_state();
        assert!(solver.solution.is_empty());
        assert!(solver.moles.is_empty());
    }

    // -----------------------------------------------------------------------
    // D.36 — validate_mutable_problem_shape()
    // -----------------------------------------------------------------------
    #[test]
    fn validate_mutable_problem_shape_accepts_valid_input() {
        let solver = EquilibriumLogMoles::empty();
        let result = solver.validate_mutable_problem_shape(&[]);
        assert!(result.is_ok());
    }

    #[test]
    fn validate_mutable_problem_shape_rejects_non_finite_moles() {
        let solver = EquilibriumLogMoles::empty();
        let result = solver.validate_mutable_problem_shape(&[f64::NAN]);
        assert!(result.is_err());
    }

    #[test]
    fn validate_mutable_problem_shape_rejects_negative_moles() {
        let solver = EquilibriumLogMoles::empty();
        let result = solver.validate_mutable_problem_shape(&[-1.0]);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.37 — stage_equilibrium_system()
    // -----------------------------------------------------------------------
    #[test]
    #[should_panic(expected = "Cannot compute the SVD of an empty matrix")]
    fn stage_equilibrium_system_panics_on_empty_solver() {
        let solver = EquilibriumLogMoles::empty();
        let _ = solver.stage_equilibrium_system(&[]);
    }

    // -----------------------------------------------------------------------
    // D.39 — build_solve_contract()
    // -----------------------------------------------------------------------
    #[test]
    fn build_solve_contract_fails_on_empty_solver() {
        let mut solver = EquilibriumLogMoles::empty();
        let result = solver.build_solve_contract();
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.40 — TemperatureWorkerSeed getters
    // -----------------------------------------------------------------------
    #[test]
    fn temperature_worker_seed_getters_return_expected_values() {
        let solver = EquilibriumLogMoles::empty();
        let seed = solver.temperature_worker_seed();
        assert_eq!(seed.species_count(), 0);
        assert!(seed.n0().is_empty());
        assert!(seed.pressure().is_finite());
        assert!(seed.reference_pressure().is_finite());
        assert!(seed.species_eps().is_finite());
        assert!(seed.substate_eps().is_finite());
    }

    // -----------------------------------------------------------------------
    // D.42 — solver_impl()
    // -----------------------------------------------------------------------
    #[test]
    fn solver_impl_fails_on_empty_solver_with_empty_policy() {
        let solver = EquilibriumLogMoles::empty();
        let f = |_: &[f64]| Ok(vec![0.0]);
        let feasible = |_: &[f64]| true;
        let validate = |_: &[f64]| {
            Ok(EquilibriumCandidateReport {
                residual_l2_norm: 0.0,
                residual_rms: 0.0,
                max_abs_residual: 0.0,
                raw_residual_l2_norm: 0.0,
                raw_residual_rms: 0.0,
                raw_max_abs_residual: 0.0,
                max_abs_element_balance_error: 0.0,
                reaction_affinity_l2_norm: 0.0,
                max_abs_reaction_affinity: 0.0,
                min_moles: 0.0,
            })
        };
        // Empty policy with no backends
        let policy = SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM));
        let result = solver.solver_impl(
            vec![0.0],
            &f,
            None,
            &feasible,
            &validate,
            policy,
            None,
        );
        // Should fail because solver has no stoich matrix etc.
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.43 — solve_backend_cascade()
    // -----------------------------------------------------------------------
    #[test]
    fn solve_backend_cascade_rejects_empty_backends() {
        let f = |_: &[f64]| Ok(vec![0.0]);
        let feasible = |_: &[f64]| true;
        let validate = |_: &[f64]| {
            Ok(EquilibriumCandidateReport {
                residual_l2_norm: 0.0,
                residual_rms: 0.0,
                max_abs_residual: 0.0,
                raw_residual_l2_norm: 0.0,
                raw_residual_rms: 0.0,
                raw_max_abs_residual: 0.0,
                max_abs_element_balance_error: 0.0,
                reaction_affinity_l2_norm: 0.0,
                max_abs_reaction_affinity: 0.0,
                min_moles: 0.0,
            })
        };
        let budget = SolverCascadeBudget::new(1, 100, 100);
        let params = SolverParams::default();
        let stoich = DMatrix::from_row_slice(0, 0, &[]);
        let result = EquilibriumLogMoles::solve_backend_cascade(
            &[],
            vec![0.0],
            &f,
            None,
            &feasible,
            &validate,
            SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM)),
            budget,
            &params,
            &[],
            &stoich,
            None,
        );
        assert!(result.is_err());
    }

    #[test]
    fn solve_backend_cascade_rejects_zero_budget() {
        let f = |_: &[f64]| Ok(vec![0.0]);
        let feasible = |_: &[f64]| true;
        let validate = |_: &[f64]| {
            Ok(EquilibriumCandidateReport {
                residual_l2_norm: 0.0,
                residual_rms: 0.0,
                max_abs_residual: 0.0,
                raw_residual_l2_norm: 0.0,
                raw_residual_rms: 0.0,
                raw_max_abs_residual: 0.0,
                max_abs_element_balance_error: 0.0,
                reaction_affinity_l2_norm: 0.0,
                max_abs_reaction_affinity: 0.0,
                min_moles: 0.0,
            })
        };
        let budget = SolverCascadeBudget::new(0, 0, 0);
        let params = SolverParams::default();
        let stoich = DMatrix::from_row_slice(0, 0, &[]);
        let result = EquilibriumLogMoles::solve_backend_cascade(
            &[],
            vec![0.0],
            &f,
            None,
            &feasible,
            &validate,
            SolverPolicy::Single(SolverBackend::Legacy(Solvers::LM)),
            budget,
            &params,
            &[],
            &stoich,
            None,
        );
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.44 — recoverable_backend_failure_kind()
    // -----------------------------------------------------------------------
    #[test]
    fn recoverable_backend_failure_kind_returns_none_for_invalid_problem() {
        let err = ReactionExtentError::InvalidProblem {
            field: "test",
            message: "test error".to_string(),
        };
        assert!(recoverable_backend_failure_kind(&err).is_none());
    }

    #[test]
    fn recoverable_backend_failure_kind_returns_solver_for_solve_error() {
        let err = ReactionExtentError::SolveError(SolveError::MaxIterations);
        let kind = recoverable_backend_failure_kind(&err);
        assert_eq!(kind, Some(SolverAttemptFailureKind::Solver));
    }

    // -----------------------------------------------------------------------
    // D.51 — validate_logmole_system_dimensions()
    // -----------------------------------------------------------------------
    #[test]
    fn validate_logmole_system_dimensions_rejects_element_matrix_mismatch() {
        // 2 species, 1 reaction, 1 element → 2 = 1 + 1 ✓
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        // Wrong: 1 species row instead of 2
        let elements = DMatrix::from_row_slice(1, 1, &[1.0]);
        let species_phase = vec![0, 0];
        let delta_n = vec![vec![1.0]];
        let result = validate_logmole_system_dimensions(
            &reactions, &elements, &species_phase, 1, &delta_n,
        );
        assert!(result.is_err());
    }

    #[test]
    fn validate_logmole_system_dimensions_rejects_species_phase_mismatch() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        // Only 1 entry for 2 species
        let species_phase = vec![0];
        let delta_n = vec![vec![1.0]];
        let result = validate_logmole_system_dimensions(
            &reactions, &elements, &species_phase, 1, &delta_n,
        );
        assert!(result.is_err());
    }

    #[test]
    fn validate_logmole_system_dimensions_rejects_out_of_bounds_phase() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        let species_phase = vec![0, 5]; // phase 5 > phase_count 1
        let delta_n = vec![vec![1.0]];
        let result = validate_logmole_system_dimensions(
            &reactions, &elements, &species_phase, 1, &delta_n,
        );
        assert!(result.is_err());
    }

    #[test]
    fn validate_logmole_system_dimensions_rejects_delta_n_mismatch() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        let species_phase = vec![0, 0];
        // delta_n has 2 reaction rows but only 1 reaction
        let delta_n = vec![vec![1.0], vec![0.0]];
        let result = validate_logmole_system_dimensions(
            &reactions, &elements, &species_phase, 1, &delta_n,
        );
        assert!(result.is_err());
    }

    #[test]
    fn validate_logmole_system_dimensions_accepts_valid_dimensions() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 2.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        let species_phase = vec![0, 0];
        let delta_n = vec![vec![1.0]];
        let result = validate_logmole_system_dimensions(
            &reactions, &elements, &species_phase, 1, &delta_n,
        );
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // D.52 — validate_residual_conditions()
    // -----------------------------------------------------------------------
    #[test]
    fn validate_residual_conditions_accepts_valid_conditions() {
        let result = validate_residual_conditions(300.0, 101325.0, 101325.0);
        assert!(result.is_ok());
    }

    #[test]
    fn validate_residual_conditions_rejects_non_finite_temperature() {
        let result = validate_residual_conditions(f64::NAN, 101325.0, 101325.0);
        assert!(result.is_err());
    }

    #[test]
    fn validate_residual_conditions_rejects_zero_pressure() {
        let result = validate_residual_conditions(300.0, 0.0, 101325.0);
        assert!(result.is_err());
    }

    #[test]
    fn validate_residual_conditions_rejects_negative_reference_pressure() {
        let result = validate_residual_conditions(300.0, 101325.0, -1.0);
        assert!(result.is_err());
    }

    // -----------------------------------------------------------------------
    // D.31 — has_rst_symbolic_context() — direct test
    // -----------------------------------------------------------------------
    #[test]
    fn has_rst_symbolic_context_returns_false_for_empty_solver() {
        let solver = EquilibriumLogMoles::empty();
        assert!(!solver.has_rst_symbolic_context());
    }

    // -----------------------------------------------------------------------
    // D.23 — temperature_failure() — direct test
    // -----------------------------------------------------------------------
    #[test]
    fn temperature_failure_creates_failure_from_solve_error() {
        let err = ReactionExtentError::SolveError(SolveError::MaxIterations);
        let failure = temperature_failure(3000.0, &err);
        assert_eq!(failure.temperature, 3000.0);
        assert!(failure.message.contains("iteration limit"));
        assert!(failure.attempts.is_empty());
    }

    #[test]
    fn temperature_failure_extracts_attempts_from_all_backends_failed() {
        let err = ReactionExtentError::AllBackendsFailed { attempts: vec![] };
        let failure = temperature_failure(3000.0, &err);
        assert!(failure.attempts.is_empty());
    }
}