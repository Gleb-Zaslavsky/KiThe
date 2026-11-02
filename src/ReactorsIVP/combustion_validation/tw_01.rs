use super::matrix_exponential::expm_action;
use crate::Kinetics::mechfinder_api::ReactionData;
use crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig;
use crate::ReactorsIVP::SimpleReactorIVP::{IvpSolveSnapshot, SimpleReactorTask};
use crate::ReactorsIVP::solver_backend::{ReactorIvpSolverConfig, ReactorIvpSolverFamily};
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;

const FROZEN_PROFILE: &str = include_str!("fixtures/tw_01_profiles.csv");
pub(super) const SPECIES_COUNT: usize = 12;
pub(super) const RHO: f64 = 1_000.0;
const CP: f64 = 1_000.0;
pub(super) const LAMBDA: f64 = 1.0;
pub(super) const LENGTH: f64 = 0.01;
const MASS_FLOW: f64 = 1.0e-4;
pub(super) const DT: f64 = 300.0;
pub(super) const T_SCALE: f64 = 100.0;
pub(super) const INITIAL_T: f64 = 400.0;
pub(super) const INITIAL_Q: f64 = 0.05;
pub(super) const INITIAL_MASS_FRACTIONS: [f64; SPECIES_COUNT] = [
    0.12, 0.08, 0.07, 0.09, 0.10, 0.06, 0.11, 0.05, 0.07, 0.08, 0.09, 0.08,
];
pub(super) const REACTIONS: [(usize, usize, f64); 20] = [
    (0, 1, 0.70),
    (1, 2, 1.10),
    (2, 3, 0.50),
    (3, 4, 1.60),
    (4, 5, 0.80),
    (5, 6, 1.30),
    (6, 7, 0.60),
    (7, 8, 1.90),
    (8, 9, 0.90),
    (9, 10, 1.40),
    (10, 11, 0.75),
    (11, 0, 1.20),
    (0, 4, 0.45),
    (1, 5, 0.95),
    (2, 6, 1.70),
    (3, 7, 0.55),
    (4, 8, 1.05),
    (5, 9, 1.80),
    (6, 10, 0.65),
    (7, 11, 1.50),
];

pub(super) fn species_name(index: usize) -> String {
    format!("X{index}")
}

pub(super) fn tw01_task(family: ReactorIvpSolverFamily) -> SimpleReactorTask {
    let mut task = SimpleReactorTask::new();
    let groups = (0..SPECIES_COUNT)
        .map(|index| {
            (
                species_name(index),
                HashMap::from([
                    ("C".to_string(), 2),
                    ("H".to_string(), 4),
                    ("O".to_string(), 2),
                ]),
            )
        })
        .collect::<HashMap<_, _>>();
    let reactions = REACTIONS
        .iter()
        .map(|(from, to, rate)| {
            ReactionData::new_elementary(
                format!("{}=>{}", species_name(*from), species_name(*to)),
                vec![*rate, 0.0, 0.0],
                None,
            )
        })
        .collect();
    task.kindata
        .set_reaction_data_directly(reactions, Some(groups))
        .expect("TW-01 reaction rows and formula groups should parse");

    let canonical = (0..SPECIES_COUNT).map(species_name).collect::<Vec<_>>();
    task.kindata.substances = canonical.iter().rev().cloned().collect();
    task.kindata.stecheodata.vec_of_molmasses = Some(vec![60.052; SPECIES_COUNT]);
    task.thermal_effects = vec![0.0; REACTIONS.len()];
    task.ro = RHO;
    task.Cp = CP;
    task.Lambda = LAMBDA;
    task.m = MASS_FLOW;
    task.scaling = ScalingConfig::new(DT, LENGTH, T_SCALE);
    task.initial_conditions =
        HashMap::from([("T".to_string(), INITIAL_T), ("q".to_string(), INITIAL_Q)]);
    for (index, fraction) in INITIAL_MASS_FRACTIONS.iter().copied().enumerate() {
        task.initial_conditions
            .insert(species_name(index), fraction);
    }
    task.solver_backend_config = ReactorIvpSolverConfig::default()
        .with_family(family)
        .with_integration_domain(0.0, 1.0)
        .with_max_step(0.01)
        .with_rtol(1e-9)
        .with_atol(1e-11);
    task
}

pub(super) fn generator() -> DMatrix<f64> {
    let scale = RHO * LENGTH.powi(2);
    let mut matrix = DMatrix::zeros(SPECIES_COUNT, SPECIES_COUNT);
    for (from, to, rate) in REACTIONS {
        matrix[(from, from)] += scale * rate;
        matrix[(to, from)] -= scale * rate;
    }
    matrix
}

pub(super) fn initial_state() -> DVector<f64> {
    DVector::from_column_slice(&INITIAL_MASS_FRACTIONS)
}

pub(super) fn analytical_profile(z: f64) -> (DVector<f64>, f64, f64, f64) {
    let c = expm_action(&(z * generator()), &initial_state());
    let pe_q = LENGTH * MASS_FLOW * CP / LAMBDA;
    let q = INITIAL_Q * (pe_q * z).exp();
    let theta = (INITIAL_T - DT) / T_SCALE + INITIAL_Q / LAMBDA * (pe_q * z).exp_m1() / pe_q;
    let temperature = DT + T_SCALE * theta;
    (c, theta, q, temperature)
}

fn solution_column(snapshot: &IvpSolveSnapshot, name: &str) -> usize {
    snapshot
        .variable_names
        .iter()
        .position(|candidate| candidate == name)
        .unwrap_or_else(|| panic!("missing solver state {name}"))
}

fn species_index_from_state_name(name: &str, species_order: &[String]) -> usize {
    name.strip_prefix('C')
        .and_then(|index| index.parse::<usize>().ok())
        .and_then(|state_index| species_order.get(state_index))
        .and_then(|species| species.strip_prefix('X'))
        .and_then(|index| index.parse::<usize>().ok())
        .expect("concentration state should map to a named TW-01 species")
}

fn assert_solution_against_oracle(
    snapshot: &IvpSolveSnapshot,
    species_order: &[String],
    family: ReactorIvpSolverFamily,
) {
    let theta_index = solution_column(snapshot, "Teta");
    let q_index = solution_column(snapshot, "q");
    let mut max_species_error: f64 = 0.0;
    let mut max_theta_error: f64 = 0.0;
    let mut max_q_error: f64 = 0.0;
    let mut max_temperature_error: f64 = 0.0;
    let mut max_species_normalized: f64 = 0.0;
    let mut max_theta_normalized: f64 = 0.0;
    let mut max_q_normalized: f64 = 0.0;
    let mut minimum_fraction: f64 = f64::INFINITY;
    let mut maximum_mass_balance_error: f64 = 0.0;
    let mut maximum_element_balance_error: f64 = 0.0;
    let initial_total = INITIAL_MASS_FRACTIONS.iter().sum::<f64>();

    assert_eq!(snapshot.variable_count(), SPECIES_COUNT + 2);
    assert!(snapshot.status.starts_with("finished"));
    for (row, z) in snapshot.axis.iter().copied().enumerate() {
        let (expected_c, expected_theta, expected_q, expected_temperature) = analytical_profile(z);
        let actual_theta = snapshot.values[(row, theta_index)];
        let actual_q = snapshot.values[(row, q_index)];
        let actual_temperature = DT + T_SCALE * actual_theta;
        max_theta_error = max_theta_error.max((actual_theta - expected_theta).abs());
        max_q_error = max_q_error.max((actual_q - expected_q).abs());
        max_temperature_error =
            max_temperature_error.max((actual_temperature - expected_temperature).abs());

        let mut sum = 0.0;
        let mut elemental_mass_fractions = [0.0; 3];
        for (state_index, state_name) in snapshot.variable_names.iter().enumerate() {
            if !state_name.starts_with('C') {
                continue;
            }
            let canonical_index = species_index_from_state_name(state_name, species_order);
            let actual = snapshot.values[(row, state_index)];
            let expected = expected_c[canonical_index];
            let error = (actual - expected).abs();
            max_species_error = max_species_error.max(error);
            minimum_fraction = minimum_fraction.min(actual);
            sum += actual;
            for (element_index, (atom_count, atomic_mass)) in
                [(2.0, 12.011), (4.0, 1.008), (2.0, 15.999)]
                    .into_iter()
                    .enumerate()
            {
                elemental_mass_fractions[element_index] +=
                    actual * atom_count * atomic_mass / 60.052;
            }
            assert!(
                actual.is_finite() && actual > 0.0,
                "{family:?} {state_name} at z={z}: {actual}"
            );
            let normalized_error = error / (1e-9 + 1e-7 * expected.abs());
            max_species_normalized = max_species_normalized.max(normalized_error);
            assert!(
                normalized_error <= 1.0,
                "{family:?} {state_name} at z={z}: actual={actual:.16e}, reference={expected:.16e}, abs_error={error:.3e}, normalized_error={normalized_error:.3e}"
            );
        }
        let balance_error = (sum - initial_total).abs();
        maximum_mass_balance_error = maximum_mass_balance_error.max(balance_error);
        assert!(
            balance_error <= 2e-9,
            "{family:?} mass-fraction sum at z={z}: {sum:.16e}, inlet={initial_total:.16e}"
        );

        for (name, actual, expected) in [
            ("Teta", actual_theta, expected_theta),
            ("q", actual_q, expected_q),
        ] {
            let error = (actual - expected).abs();
            let normalized_error = error / (1e-9 + 1e-7 * expected.abs());
            if name == "Teta" {
                max_theta_normalized = max_theta_normalized.max(normalized_error);
            } else {
                max_q_normalized = max_q_normalized.max(normalized_error);
            }
            assert!(
                normalized_error <= 1.0,
                "{family:?} {name} at z={z}: actual={actual:.16e}, reference={expected:.16e}, abs_error={error:.3e}, normalized_error={normalized_error:.3e}"
            );
        }
        assert!(actual_temperature.is_finite() && (0.0..=3_000.0).contains(&actual_temperature));

        for (element_index, (atom_count, atomic_mass)) in
            [(2.0, 12.011), (4.0, 1.008), (2.0, 15.999)]
                .into_iter()
                .enumerate()
        {
            let expected = initial_total * atom_count * atomic_mass / 60.052;
            let error = (elemental_mass_fractions[element_index] - expected).abs();
            maximum_element_balance_error = maximum_element_balance_error.max(error);
            assert!(error <= 4e-11);
        }
    }
    println!(
        "TW-01 family={family:?} points={} max_species_abs={max_species_error:.3e} max_species_norm={max_species_normalized:.3e} max_theta_abs={max_theta_error:.3e} max_theta_norm={max_theta_normalized:.3e} max_q_abs={max_q_error:.3e} max_q_norm={max_q_normalized:.3e} max_T_K_abs={max_temperature_error:.3e} max_mass_balance={maximum_mass_balance_error:.3e} max_element_balance={maximum_element_balance_error:.3e} min_fraction={minimum_fraction:.6e}",
        snapshot.axis.len()
    );
}

fn assert_dimensional_publication(task: &SimpleReactorTask, snapshot: &IvpSolveSnapshot) {
    let axis = task
        .solver
        .x_mesh
        .as_ref()
        .expect("postprocessed physical axis");
    let values = task
        .solver
        .solution
        .as_ref()
        .expect("postprocessed physical states");
    assert_eq!(axis.len(), snapshot.axis.len());
    assert_eq!(values.shape(), snapshot.values.shape());
    let theta = solution_column(snapshot, "Teta");
    let q = solution_column(snapshot, "q");
    for row in 0..axis.len() {
        assert!((axis[row] - LENGTH * snapshot.axis[row]).abs() <= 2e-15);
        assert!(
            (values[(row, theta)] - (DT + T_SCALE * snapshot.values[(row, theta)])).abs() <= 2e-12
        );
        assert!((values[(row, q)] - (T_SCALE / LENGTH * snapshot.values[(row, q)])).abs() <= 2e-10);
        for column in 2..values.ncols() {
            assert_eq!(values[(row, column)], snapshot.values[(row, column)]);
        }
    }
}

fn max_absolute_profile_error(snapshot: &IvpSolveSnapshot, species_order: &[String]) -> f64 {
    let theta = solution_column(snapshot, "Teta");
    let q = solution_column(snapshot, "q");
    let mut maximum: f64 = 0.0;
    for (row, z) in snapshot.axis.iter().copied().enumerate() {
        let (expected_c, expected_theta, expected_q, _) = analytical_profile(z);
        maximum = maximum.max((snapshot.values[(row, theta)] - expected_theta).abs());
        maximum = maximum.max((snapshot.values[(row, q)] - expected_q).abs());
        for (column, name) in snapshot.variable_names.iter().enumerate() {
            if name.starts_with('C') {
                let index = species_index_from_state_name(name, species_order);
                maximum = maximum.max((snapshot.values[(row, column)] - expected_c[index]).abs());
            }
        }
    }
    maximum
}

#[test]
fn tw01_generator_and_pade_reference_are_self_consistent() {
    let matrix = generator();
    assert_eq!(matrix.shape(), (SPECIES_COUNT, SPECIES_COUNT));
    assert!(matrix.iter().all(|value| value.is_finite()));
    assert!((INITIAL_MASS_FRACTIONS.iter().sum::<f64>() - 1.0).abs() <= 2e-16);

    for column in 0..SPECIES_COUNT {
        let sum = (0..SPECIES_COUNT)
            .map(|row| matrix[(row, column)])
            .sum::<f64>();
        assert!(sum.abs() < 2e-16, "column={column} sum={sum}");
    }
    assert!((matrix[(0, 0)] - 0.115).abs() < 1e-15);
    assert!((matrix[(0, 11)] + 0.12).abs() < 1e-15);
    assert!((matrix[(4, 4)] - 0.185).abs() < 1e-15);
    assert!((matrix[(4, 3)] + 0.16).abs() < 1e-15);
    assert!((matrix[(4, 0)] + 0.045).abs() < 1e-15);
    assert!((matrix.transpose() * DVector::from_element(SPECIES_COUNT, 1.0)).norm() < 1e-14);

    for z in (0..=10).map(|index| index as f64 / 10.0) {
        let (c, theta, q, temperature) = analytical_profile(z);
        assert!((c.sum() - INITIAL_MASS_FRACTIONS.iter().sum::<f64>()).abs() < 2e-14);
        assert!(c.iter().all(|value| value.is_finite() && *value > 0.0));
        assert!(theta.is_finite() && q.is_finite() && temperature.is_finite());
    }
    let (_, theta_1, q_1, temperature_1) = analytical_profile(1.0);
    assert!((theta_1 - 1.050_025_008_335_42).abs() < 2e-14);
    assert!((q_1 - 0.050_050_025_008_335_4).abs() < 2e-16);
    assert!((temperature_1 - 405.002_500_833_542).abs() < 2e-12);
}

#[test]
fn tw01_production_assembly_matches_manual_reaction_table() {
    let mut task = tw01_task(ReactorIvpSolverFamily::Lsode2);
    task.setup_condensed_ivp()
        .expect("TW-01 should assemble through KiThe");

    assert_eq!(task.kindata.vec_of_equations.len(), REACTIONS.len());
    assert_eq!(
        task.kindata.substances,
        (0..SPECIES_COUNT)
            .rev()
            .map(species_name)
            .collect::<Vec<_>>()
    );
    assert_eq!(task.solver.unknowns.len(), SPECIES_COUNT + 2);
    assert_eq!(
        task.kindata.stecheodata.stecheo_matrx.len(),
        REACTIONS.len()
    );
    assert_eq!(
        task.kindata.stecheodata.vec_of_molmasses.as_ref().unwrap(),
        &vec![60.052; SPECIES_COUNT]
    );
    assert_eq!(
        task.kindata
            .stecheodata
            .unique_vec_of_elems
            .as_ref()
            .unwrap(),
        &["C", "H", "O"]
    );
    let elemental_matrix = task
        .kindata
        .stecheodata
        .matrix_of_elements
        .as_ref()
        .unwrap();
    assert_eq!(elemental_matrix.shape(), (SPECIES_COUNT, 3));
    for row in 0..SPECIES_COUNT {
        assert_eq!(elemental_matrix[(row, 0)], 2.0);
        assert_eq!(elemental_matrix[(row, 1)], 4.0);
        assert_eq!(elemental_matrix[(row, 2)], 2.0);
    }

    for (reaction_index, (from, to, rate)) in REACTIONS.iter().copied().enumerate() {
        let row = &task.kindata.stecheodata.stecheo_matrx[reaction_index];
        let mut expected = vec![0.0; SPECIES_COUNT];
        expected[from] = -1.0;
        expected[to] = 1.0;
        for (species_index, species_name) in task.kindata.substances.iter().enumerate() {
            let canonical_index = species_name
                .strip_prefix('X')
                .unwrap()
                .parse::<usize>()
                .unwrap();
            assert_eq!(row[species_index], expected[canonical_index]);
        }
        let k = &task.kindata.K_sym_vec.as_ref().unwrap()[reaction_index];
        for theta in [0.0, 1.0, 9.0] {
            let evaluated = k.eval_expression(&["Teta"], &[theta]);
            assert!((evaluated - rate).abs() <= rate.abs() * 1e-14);
        }
    }

    for (reaction_index, _) in REACTIONS.iter().enumerate() {
        let row = &task.kindata.stecheodata.stecheo_matrx[reaction_index];
        for element_index in 0..3 {
            let elemental_delta = (0..SPECIES_COUNT)
                .map(|species_index| {
                    elemental_matrix[(species_index, element_index)] * row[species_index]
                })
                .sum::<f64>();
            assert_eq!(
                elemental_delta, 0.0,
                "reaction {reaction_index} element {element_index}"
            );
        }
    }

    let initial = task
        .build_initial_state_vector()
        .expect("valid initial state");
    let variable_names = task
        .solver
        .unknowns
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    let zero_heat = task
        .heat_release
        .eval_expression(&variable_names, initial.as_slice());
    assert!(
        zero_heat.abs() <= 1e-14,
        "zero-heat TW-01 evaluated to {zero_heat}"
    );
    assert!((task.Pe_q - 0.001).abs() < 1e-15);
    assert!((initial[0] - 1.0).abs() < 1e-15);
    assert!((initial[1] - INITIAL_Q).abs() < 1e-15);

    let state_composition = task
        .kindata
        .substances
        .iter()
        .map(|species| {
            let index = species.strip_prefix('X').unwrap().parse::<usize>().unwrap();
            INITIAL_MASS_FRACTIONS[index]
        })
        .collect::<Vec<_>>();
    let matrix = generator();
    for state_index in 0..SPECIES_COUNT {
        let canonical_row = task.kindata.substances[state_index]
            .strip_prefix('X')
            .unwrap()
            .parse::<usize>()
            .unwrap();
        let expected = (0..SPECIES_COUNT)
            .map(|source_index| {
                let canonical_column = task.kindata.substances[source_index]
                    .strip_prefix('X')
                    .unwrap()
                    .parse::<usize>()
                    .unwrap();
                matrix[(canonical_row, canonical_column)] * state_composition[source_index]
            })
            .sum::<f64>();
        let actual = task.solver.eq_system[state_index + 2]
            .eval_expression(&variable_names, initial.as_slice());
        assert!((actual - expected).abs() <= 2e-13);
    }
}

#[test]
fn tw01_frozen_analytical_profile_matches_matrix_exponential() {
    let mut lines = FROZEN_PROFILE.lines();
    let header = lines.next().expect("frozen CSV header");
    assert!(header.starts_with("z,Teta,q,T_K,X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11"));
    let mut rows = 0;
    for line in lines {
        let values = line
            .split(',')
            .map(|value| value.parse::<f64>().expect("numeric frozen CSV cell"))
            .collect::<Vec<_>>();
        assert_eq!(values.len(), SPECIES_COUNT + 4);
        let z = values[0];
        let (species, theta, q, temperature) = analytical_profile(z);
        for (column, (actual, expected)) in
            [(values[1], theta), (values[2], q), (values[3], temperature)]
                .into_iter()
                .enumerate()
        {
            assert!(
                (actual - expected).abs() <= 2e-13,
                "frozen scalar column={} z={z} actual={actual:.17e} expected={expected:.17e}",
                column + 1
            );
        }
        for species_index in 0..SPECIES_COUNT {
            assert!(
                (values[species_index + 4] - species[species_index]).abs() <= 2e-13,
                "frozen X{species_index} at z={z} differs from matrix exponential"
            );
        }
        rows += 1;
    }
    assert_eq!(rows, 11);
}

#[test]
fn tw01_analytical_reference_accepts_all_stiff_api_families() {
    for family in [
        ReactorIvpSolverFamily::Lsode2,
        ReactorIvpSolverFamily::Bdf,
        ReactorIvpSolverFamily::Radau,
    ] {
        let mut task = tw01_task(family);
        let snapshot = task
            .solve_condensed_ivp()
            .unwrap_or_else(|error| panic!("TW-01 {family:?} solve failed: {error}"));
        assert_solution_against_oracle(&snapshot, &task.kindata.substances, family);
        assert_dimensional_publication(&task, &snapshot);
        assert_eq!(
            task.solver.x_mesh.as_ref().unwrap().len(),
            snapshot.axis.len()
        );
    }
}

#[test]
fn tw01_tighter_solver_tolerances_reduce_profile_error() {
    let mut coarse = tw01_task(ReactorIvpSolverFamily::Lsode2);
    coarse.solver_backend_config = coarse.solver_backend_config.with_rtol(1e-5).with_atol(1e-7);
    let coarse_snapshot = coarse
        .solve_condensed_ivp()
        .expect("coarse TW-01 LSODE2 solve should finish");
    let coarse_error = max_absolute_profile_error(&coarse_snapshot, &coarse.kindata.substances);

    let mut fine = tw01_task(ReactorIvpSolverFamily::Lsode2);
    let fine_snapshot = fine
        .solve_condensed_ivp()
        .expect("fine TW-01 LSODE2 solve should finish");
    let fine_error = max_absolute_profile_error(&fine_snapshot, &fine.kindata.substances);
    println!(
        "TW-01 tolerance refinement coarse_max_abs={coarse_error:.3e} fine_max_abs={fine_error:.3e}"
    );
    assert!(
        fine_error < coarse_error,
        "tightening tolerances did not improve error: coarse={coarse_error:.3e}, fine={fine_error:.3e}"
    );
}
