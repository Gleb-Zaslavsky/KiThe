#[cfg(test)]
mod tests {
    use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq::*;
    use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq2::*;
    use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
    use nalgebra::{DMatrix, DVector};
    use std::collections::HashMap;
    use std::rc::Rc;
    #[test]
    fn phase_deactivation_sets_species_to_floor() {
        let species_phase = vec![0, 0, 1, 1];
        let mut y = vec![0.0, 1.0, 2.0, 3.0]; // log-moles

        deactivate_phases(&mut y, &[1], &species_phase, -700.0);

        assert_eq!(y[0], 0.0);
        assert_eq!(y[1], 1.0);
        assert!(y[2] < -600.0);
        assert!(y[3] < -600.0);
    }
    #[test]
    fn detect_inactive_phases_basic() {
        let n_phase = vec![1.0, 1e-12, 0.5];
        let inactive = detect_inactive_phases(&n_phase, 1e-8);

        assert_eq!(inactive, vec![1]);
    }

    #[test]
    fn residual_finite_after_phase_deactivation() {
        let reactions = DMatrix::<f64>::zeros(4, 1);
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

        let species_phase = vec![0, 0, 1, 1];
        let mut y = vec![0.0, 0.0, -700.0, -700.0];

        let f = residual(&y).unwrap();
        for v in f {
            assert!(v.is_finite());
        }
    }
    #[test]
    fn jacobian_finite_after_phase_deactivation() {
        let reactions = DMatrix::<f64>::zeros(4, 1);
        let elements = DMatrix::<f64>::identity(4, 4);
        let species_phase = vec![0, 0, 1, 1];
        let delta_n = vec![vec![2.0, 0.0]];
        let mut y = vec![0.0, 0.0, -700.0, -700.0];

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

    ////////////////////////////////////////////////
    #[test]
    fn O2_O_equilibrium_logmoles_with_phase_control() {
        let T = 500.0;
        // no creation/deactivation case
        let substances = vec!["O2".to_string(), "O".to_string()];
        let mut instance = gas_solver(substances, T, 101325.0, Solvers::LM, Some("info"), false);
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
        println!("moles {:?}", moles);
        assert!(moles[0] > 0.0 && moles[1] > 0.0, "Moles should be positive");
    }

    #[test]
    fn phase_creation_from_gas() {
        // dummy case we consider O2/O different
        let T = 400.0;

        let substances = vec![
            "N2".to_string(), // gas oxygen
            "N".to_string(),  // solid oxygen
        ];

        let mut instance = gas_solver(substances, T, 101325.0, Solvers::TR, None, false);

        instance.n0 = vec![1.0, 0.0];
        instance.initial_guess = Some(vec![1.0, 1e-30]);

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
        instance.with_loglevel(Some("info"));
        instance.create_stoich_matrix().unwrap();
        instance.solve_with_phase_control().unwrap();

        let moles = &instance.moles;
        println!("moles {:?}", moles);
        // --- Assertions ---
        assert!(moles[1] > 1e-6, "solid phase was not created");
        assert!(moles[0] < 1.0, "gas phase did not lose material");
    }

    #[test]
    fn phase_creation_from_gas2() {
        let mut user_subs = SubsData::new();
        let subs = vec!["O2".to_string(), "C".to_string()];
        user_subs.substances = subs;
        let mut map_of_phases = HashMap::new();
        map_of_phases.insert("O2".to_string(), Some(Phases::Gas));
        map_of_phases.insert("C".to_string(), Some(Phases::Condensed));
        user_subs.map_of_phases = map_of_phases;
        user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
        user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);
        let mut explicit_map = HashMap::new();
        explicit_map.insert("C".to_string(), "NASA_cond".to_string());
        user_subs.set_explicis_searh_instructions(explicit_map);
        user_subs.search_substances();
        println!("user_subs {:?}", user_subs);
    }
    /*
    #[test]
    fn phase_destruction_due_to_depletion() {
        let T = 3000.0; // high T favors dissociation

        let substances = vec!["O2".to_string(), "O".to_string()];

        let mut instance = gas_solver(substances, T, 101325.0, Solvers::TR, None, false);

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

        let mut instance = gas_solver(substances, T, 101325.0, Solvers::TR, None, false);

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

        let mut instance = gas_solver(substances, T, 101325.0, Solvers::TrustRegion, None, false);

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
