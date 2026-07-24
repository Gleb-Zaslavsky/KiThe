//! Offline long-term regression matrix for canonical equilibrium cases.
//!
//! The goal is not to preserve one historical iterate per solver backend.
//! Instead we keep a compact matrix of physically meaningful fixtures that
//! must remain finite, element-conserving, and offline:
//! - O2 <-> 2O dissociation;
//! - N2 <-> 2N dissociation;
//! - a diluted N/O gas mixture;
//! - condensed phase appearance and disappearance in a synthetic fixture.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::*;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    PHASE_CONTROL_TRACE_MOLE_FLOOR, gas_solver,
};
use nalgebra::DMatrix;
use std::rc::Rc;

fn gas_fixture(
    substances: Vec<&str>,
    T: f64,
    P: f64,
    n0: Vec<f64>,
    solver: Solvers,
) -> EquilibriumLogMoles {
    let mut solver = gas_solver(
        substances.into_iter().map(|s| s.to_string()).collect(),
        T,
        P,
        solver,
        None,
        false,
    )
    .unwrap();
    solver.n0 = n0;
    solver.initial_guess = Some(
        solver
            .n0
            .iter()
            .map(|&n| n.max(PHASE_CONTROL_TRACE_MOLE_FLOOR).ln())
            .collect(),
    );
    solver.phases = vec![Phase {
        kind: PhaseKind::IdealGas,
        species: (0..solver.n0.len()).collect(),
    }];
    solver.create_stoich_matrix().unwrap();
    solver
}

fn synthetic_condensed_fixture(candidate_gibbs: f64) -> EquilibriumLogMoles {
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
    solver.solver_settings.solver_params.max_iter = 200;
    solver.phase_manager.phase_eps = 1e-12;
    solver.phase_manager.set_explicit_hysteresis(-1.0, 1.0);
    solver.create_stoich_matrix().unwrap();
    solver.phase_active_mask.clear();
    solver
}

fn assert_balanced_and_finite(instance: &EquilibriumLogMoles) {
    assert!(
        instance
            .moles
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
    );
    let validation = instance
        .last_validation_report
        .as_ref()
        .expect("solve must publish validation");
    assert!(validation.residual_l2_norm.is_finite());
    assert!(validation.max_abs_element_balance_error <= 1e-5);
    let expected = compute_element_totals(&instance.elem_composition, &instance.n0).unwrap();
    let observed = compute_element_totals(&instance.elem_composition, &instance.moles).unwrap();
    for (lhs, rhs) in expected.iter().zip(observed.iter()) {
        assert!((lhs - rhs).abs() <= 1e-5);
    }
}

#[test]
fn offline_regression_matrix_covers_canonical_gas_and_condensed_cases() {
    let cases = [
        (
            "O2 <-> 2O",
            gas_fixture(
                vec!["O2", "O"],
                500.0,
                101325.0,
                vec![1.0, 1e-5],
                Solvers::LM,
            ),
        ),
        (
            "N2 <-> 2N",
            gas_fixture(
                vec!["N2", "N"],
                5500.0,
                101325.0,
                vec![1e-5, 1.0],
                Solvers::LM,
            ),
        ),
        (
            "N/O mixture",
            gas_fixture(
                vec!["N2", "N", "O2", "O"],
                5500.0,
                101325.0,
                vec![1.0, 1e-5, 1.0, 1e-5],
                Solvers::LM,
            ),
        ),
    ];

    for (label, mut solver) in cases {
        solver.solve_with_phase_control().unwrap_or_else(|error| {
            panic!("{label} should solve offline and stay finite: {error:?}")
        });
        assert_balanced_and_finite(&solver);
    }

    let mut condensed_appears = synthetic_condensed_fixture(-20_000.0);
    condensed_appears.solve_with_phase_control().unwrap();
    assert_balanced_and_finite(&condensed_appears);
    assert_eq!(condensed_appears.phase_active_mask, vec![true, true]);
    assert!(condensed_appears.moles[2] > condensed_appears.phase_manager.phase_eps);

    let mut condensed_disappears = synthetic_condensed_fixture(100_000.0);
    condensed_disappears.solve_with_phase_control().unwrap();
    assert_balanced_and_finite(&condensed_disappears);
    assert_eq!(condensed_disappears.phase_active_mask, vec![true, false]);
    assert!(condensed_disappears.moles[2] <= PHASE_CONTROL_TRACE_MOLE_FLOOR * 1.01);
}

#[test]
fn offline_regression_matrix_covers_deterministic_small_gas_systems() {
    let cases = [
        (
            "NO/N2/O2",
            gas_fixture(
                vec!["NO", "N2", "O2"],
                1800.0,
                101325.0,
                vec![1e-4, 1.0, 1.0],
                Solvers::LM,
            ),
        ),
        (
            "N/O/O2/N2",
            gas_fixture(
                vec!["N2", "N", "O2", "O"],
                3200.0,
                101325.0,
                vec![0.8, 1e-4, 0.6, 1e-4],
                Solvers::LM,
            ),
        ),
        (
            "diluted O2/O/N2",
            gas_fixture(
                vec!["O2", "O", "N2"],
                1000.0,
                202_650.0,
                vec![1.0, 1e-5, 5.0],
                Solvers::LM,
            ),
        ),
    ];

    for (label, mut solver) in cases {
        solver.solve_with_phase_control().unwrap_or_else(|error| {
            panic!("{label} should solve offline and stay finite: {error:?}")
        });
        assert_balanced_and_finite(&solver);
    }
}
