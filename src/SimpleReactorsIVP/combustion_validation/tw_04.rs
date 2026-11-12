//! Fully coupled thermochemical regression against an independent SciPy oracle.

use super::tw_01::{
    CP, DT, INITIAL_MASS_FRACTIONS, LAMBDA, LENGTH, MASS_FLUX, REACTIONS, RHO, SPECIES_COUNT,
    T_SCALE, species_name,
};
use super::tw_02::ENTHALPIES_J_PER_MOL;
use super::tw_03::{ACTIVATION_ENERGIES_J_PER_MOL, task_for_parameters};
use crate::SimpleReactorsIVP::SimpleReactorIVP::{IvpSolveSnapshot, SimpleReactorTask};
use crate::SimpleReactorsIVP::solver_backend::ReactorIvpSolverFamily;
const FROZEN_PROFILE: &str = include_str!("fixtures/tw_04_profiles.csv");
const TW03_FROZEN_PROFILE: &str = super::tw_03::FROZEN_PROFILE;
const FROZEN_METADATA: &str = include_str!("fixtures/tw_04_reference_metadata.json");
const MOLAR_MASS_KG_PER_MOL: f64 = 0.060_052;
const REFERENCE_TEMPERATURE_K: f64 = 400.0;
const INITIAL_QHAT: f64 = 5_000.0;
const WEAK_GAMMA: f64 = 1_000.0;
const MODERATE_GAMMA: f64 = 100_000.0;
const STRONG_GAMMA: f64 = 500_000.0;
const SCENARIOS: [(&str, f64); 4] = [
    ("tw_04a", 0.0),
    ("tw_04b", WEAK_GAMMA),
    ("tw_04c", MODERATE_GAMMA),
    ("tw_04d", STRONG_GAMMA),
];

#[derive(Debug, Clone)]
struct FrozenRow {
    scenario: String,
    gamma: f64,
    z: f64,
    state: [f64; 14],
    derivative: [f64; 14],
    q_volume: f64,
}

fn frozen_rows() -> Vec<FrozenRow> {
    let mut lines = FROZEN_PROFILE.lines();
    assert_eq!(
        lines.next(),
        Some(concat!(
            "scenario,gamma,z,theta,qhat,",
            "Y0,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,",
            "dstate0,dstate1,dstate2,dstate3,dstate4,dstate5,dstate6,",
            "dstate7,dstate8,dstate9,dstate10,dstate11,dstate12,dstate13,",
            "q_volume"
        ))
    );
    lines
        .map(|line| {
            let fields = line.split(',').collect::<Vec<_>>();
            assert_eq!(fields.len(), 32, "frozen row field count");
            let values = fields[1..]
                .iter()
                .map(|value| value.parse::<f64>().expect("numeric TW-04 cell"))
                .collect::<Vec<_>>();
            let mut state = [0.0; 14];
            state.copy_from_slice(&values[2..16]);
            let mut derivative = [0.0; 14];
            derivative.copy_from_slice(&values[16..30]);
            FrozenRow {
                scenario: fields[0].to_string(),
                gamma: values[0],
                z: values[1],
                state,
                derivative,
                q_volume: values[30],
            }
        })
        .collect()
}

fn scenario_rows<'a>(rows: &'a [FrozenRow], scenario: &str) -> Vec<&'a FrozenRow> {
    rows.iter().filter(|row| row.scenario == scenario).collect()
}

fn reference_state_at(rows: &[&FrozenRow], z: f64) -> [f64; 14] {
    let right = rows.partition_point(|row| row.z < z);
    if right == 0 {
        return rows[0].state;
    }
    if right == rows.len() {
        return rows[rows.len() - 1].state;
    }
    let left_row = rows[right - 1];
    let right_row = rows[right];
    let width = right_row.z - left_row.z;
    let s = (z - left_row.z) / width;
    let h00 = 2.0 * s.powi(3) - 3.0 * s.powi(2) + 1.0;
    let h10 = s.powi(3) - 2.0 * s.powi(2) + s;
    let h01 = -2.0 * s.powi(3) + 3.0 * s.powi(2);
    let h11 = s.powi(3) - s.powi(2);
    std::array::from_fn(|index| {
        h00 * left_row.state[index]
            + h10 * width * left_row.derivative[index]
            + h01 * right_row.state[index]
            + h11 * width * right_row.derivative[index]
    })
}

fn reaction_heats(gamma: f64) -> Vec<f64> {
    REACTIONS
        .iter()
        .map(|(from, to, _)| gamma * (ENTHALPIES_J_PER_MOL[*from] - ENTHALPIES_J_PER_MOL[*to]))
        .collect()
}

fn reaction_heats_from_enthalpy(enthalpy: &[f64; SPECIES_COUNT]) -> [f64; 20] {
    std::array::from_fn(|reaction| {
        let (from, to, _) = REACTIONS[reaction];
        enthalpy[from] - enthalpy[to]
    })
}

fn task_for(gamma: f64, family: ReactorIvpSolverFamily) -> SimpleReactorTask {
    let mut task = task_for_parameters(ACTIVATION_ENERGIES_J_PER_MOL, INITIAL_QHAT, family);
    task.thermal_effects = reaction_heats(gamma);
    task.solver_backend_config = task
        .solver_backend_config
        .with_max_step(1e-3)
        .with_rtol(1e-11)
        .with_atol(1e-13);
    assert!(
        (task.initial_conditions["q"] - INITIAL_QHAT * T_SCALE / LENGTH).abs() <= 1e-8,
        "task input keeps dimensional q in W/m^2"
    );
    task
}

fn zero_heat_arrhenius_task(family: ReactorIvpSolverFamily) -> SimpleReactorTask {
    let mut task = task_for_parameters([0.0; 20], INITIAL_QHAT, family);
    task.thermal_effects = reaction_heats(MODERATE_GAMMA);
    task
}

fn state_column(snapshot: &IvpSolveSnapshot, name: &str) -> usize {
    snapshot
        .variable_names
        .iter()
        .position(|candidate| candidate == name)
        .unwrap_or_else(|| panic!("missing TW-04 state {name}"))
}

fn canonical_species_index(state_name: &str, species_order: &[String]) -> usize {
    let state_index = state_name
        .strip_prefix('C')
        .and_then(|value| value.parse::<usize>().ok())
        .expect("species state name C<n>");
    species_order[state_index]
        .strip_prefix('X')
        .and_then(|value| value.parse::<usize>().ok())
        .expect("TW-04 canonical species name X<n>")
}

fn state_in_snapshot_order(
    snapshot: &IvpSolveSnapshot,
    row: usize,
    species_order: &[String],
) -> [f64; 14] {
    let mut state = [0.0; 14];
    state[0] = snapshot.values[(row, state_column(snapshot, "Teta"))];
    state[1] = snapshot.values[(row, state_column(snapshot, "q"))];
    for (column, name) in snapshot.variable_names.iter().enumerate() {
        if name.starts_with('C') {
            let canonical = canonical_species_index(name, species_order);
            state[canonical + 2] = snapshot.values[(row, column)];
        }
    }
    state
}

fn evaluate_production_rhs(task: &SimpleReactorTask, state: &[f64; 14]) -> Vec<f64> {
    let variables = task
        .solver
        .unknowns
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    task.solver
        .eq_system
        .iter()
        .map(|equation| equation.eval_expression(&variables, state))
        .collect()
}

fn independent_rates_and_heat(
    temperature: f64,
    mass_fractions: &[f64; SPECIES_COUNT],
    gamma: f64,
) -> (Vec<f64>, f64) {
    let pre_exponentials = super::tw_03::arrhenius_pre_exponentials(&ACTIVATION_ENERGIES_J_PER_MOL);
    let heats = reaction_heats(gamma);
    let rates = REACTIONS
        .iter()
        .enumerate()
        .map(|(reaction, (from, _, _))| {
            let k = pre_exponentials[reaction]
                * (-ACTIVATION_ENERGIES_J_PER_MOL[reaction] / (8.314 * temperature)).exp();
            k * RHO * mass_fractions[*from] / MOLAR_MASS_KG_PER_MOL
        })
        .collect::<Vec<_>>();
    let heat = rates
        .iter()
        .enumerate()
        .map(|(reaction, rate)| heats[reaction] * rate)
        .sum();
    (rates, heat)
}

fn assert_physical_invariants(snapshot: &IvpSolveSnapshot, species_order: &[String], gamma: f64) {
    let initial_y = INITIAL_MASS_FRACTIONS;
    let g = ENTHALPIES_J_PER_MOL.map(|enthalpy| enthalpy / MOLAR_MASS_KG_PER_MOL);
    let mut max_raw_energy_error: f64 = 0.0;
    let mut max_normalized_energy_error: f64 = 0.0;
    for row in 0..snapshot.axis.len() {
        let state = state_in_snapshot_order(snapshot, row, species_order);
        let mass = state[2..].iter().sum::<f64>();
        assert!(
            (mass - 1.0).abs() <= 2e-9,
            "mass sum at z={}: {mass}",
            snapshot.axis[row]
        );
        assert!(
            state[2..]
                .iter()
                .all(|fraction| fraction.is_finite() && *fraction > 0.0)
        );
        for element_fraction in [24.022 / 60.052, 4.032 / 60.052, 31.998 / 60.052] {
            let elemental_mass = state[2..]
                .iter()
                .map(|fraction| fraction * element_fraction)
                .sum::<f64>();
            assert!((elemental_mass - element_fraction).abs() <= 2e-9);
        }
        let theta_delta = state[0] - (REFERENCE_TEMPERATURE_K - DT) / T_SCALE;
        let energy = state[1]
            - INITIAL_QHAT
            - gamma * MASS_FLUX * LENGTH / T_SCALE
                * (0..SPECIES_COUNT)
                    .map(|species| g[species] * (state[species + 2] - initial_y[species]))
                    .sum::<f64>()
            - LAMBDA * (MASS_FLUX * CP * LENGTH / LAMBDA) * theta_delta;
        let scale = state[1].abs()
            + INITIAL_QHAT.abs()
            + (gamma * MASS_FLUX * LENGTH / T_SCALE).abs()
                * (0..SPECIES_COUNT)
                    .map(|species| (g[species] * (state[species + 2] - initial_y[species])).abs())
                    .sum::<f64>()
            + (LAMBDA * (MASS_FLUX * CP * LENGTH / LAMBDA) * theta_delta).abs()
            + 1.0;
        max_raw_energy_error = max_raw_energy_error.max(energy.abs());
        max_normalized_energy_error = max_normalized_energy_error.max(energy.abs() / scale);
    }
    assert!(
        max_normalized_energy_error <= 1e-8,
        "normalized energy error {max_normalized_energy_error:.3e}"
    );
    assert!(
        max_raw_energy_error <= 1e-3,
        "energy invariant raw error {max_raw_energy_error:.3e}"
    );
}

fn assert_snapshot_matches_reference(
    scenario: &str,
    family: ReactorIvpSolverFamily,
    gamma: f64,
    snapshot: &IvpSolveSnapshot,
    species_order: &[String],
    rows: &[&FrozenRow],
) {
    let theta_column = state_column(snapshot, "Teta");
    let q_column = state_column(snapshot, "q");
    let mut max_y_error: f64 = 0.0;
    let mut max_temperature_error: f64 = 0.0;
    let mut max_qhat_error: f64 = 0.0;
    let mut max_physical_heat_flux_error: f64 = 0.0;
    let mut min_fraction = f64::INFINITY;
    let mut max_heat_error: f64 = 0.0;
    for row in 0..snapshot.axis.len() {
        let z = snapshot.axis[row];
        let expected = reference_state_at(rows, z);
        let actual = state_in_snapshot_order(snapshot, row, species_order);
        let expected_temperature = DT + T_SCALE * expected[0];
        let actual_temperature = DT + T_SCALE * actual[0];
        let temperature_error = (actual_temperature - expected_temperature).abs();
        assert!(
            temperature_error <= 2e-6 + 1e-7 * expected_temperature.abs(),
            "{scenario} {family:?} T z={z}: actual={actual_temperature:.16e} K expected={expected_temperature:.16e} K"
        );
        assert!(
            (actual[0] - expected[0]).abs() <= 1e-9 + 1e-7 * expected[0].abs(),
            "{scenario} {family:?} Teta z={z}: actual={} expected={}",
            actual[0],
            expected[0]
        );
        assert!(
            (actual[1] - expected[1]).abs() <= 1e-5 + 1e-7 * expected[1].abs(),
            "{scenario} {family:?} qhat z={z}: actual={} expected={}",
            actual[1],
            expected[1]
        );
        let physical_heat_flux_error = (actual[1] - expected[1]).abs() * T_SCALE / LENGTH;
        let expected_physical_heat_flux = expected[1] * T_SCALE / LENGTH;
        assert!(
            physical_heat_flux_error <= 2e-2 + 1e-9 * expected_physical_heat_flux.abs(),
            "{scenario} {family:?} physical q z={z}: abs error={physical_heat_flux_error:.3e} W/m^2"
        );
        assert!((snapshot.values[(row, theta_column)] - actual[0]).abs() <= f64::EPSILON);
        assert!((snapshot.values[(row, q_column)] - actual[1]).abs() <= f64::EPSILON);
        for species in 0..SPECIES_COUNT {
            let expected_y = expected[species + 2];
            let actual_y = actual[species + 2];
            let error = (actual_y - expected_y).abs();
            assert!(
                error <= 1e-9 + 1e-7 * expected_y.abs(),
                "{scenario} {family:?} Y{species} z={z}: actual={actual_y:.16e} expected={expected_y:.16e} error={error:.3e}"
            );
            assert!(actual_y.is_finite() && actual_y > 0.0);
            min_fraction = min_fraction.min(actual_y);
            max_y_error = max_y_error.max(error);
        }
        let (_, expected_heat) = independent_rates_and_heat(
            expected_temperature,
            &expected[2..].try_into().unwrap(),
            gamma,
        );
        let (_, actual_heat) =
            independent_rates_and_heat(actual_temperature, &actual[2..].try_into().unwrap(), gamma);
        max_heat_error = max_heat_error.max((actual_heat - expected_heat).abs());
        max_temperature_error = max_temperature_error.max(temperature_error);
        max_qhat_error = max_qhat_error.max((actual[1] - expected[1]).abs());
        max_physical_heat_flux_error = max_physical_heat_flux_error.max(physical_heat_flux_error);
    }
    assert_physical_invariants(snapshot, species_order, gamma);
    println!(
        "{scenario} family={family:?} gamma={gamma:.0} points={} max_dY={max_y_error:.3e} max_dT_K={max_temperature_error:.3e} max_dqhat={max_qhat_error:.3e} max_dq_Wm2={max_physical_heat_flux_error:.3e} minY={min_fraction:.6e} heat_delta={max_heat_error:.3e}",
        snapshot.axis.len()
    );
}

#[test]
fn tw04_independent_fixture_has_frozen_scenarios_and_provenance() {
    let metadata: serde_json::Value =
        serde_json::from_str(FROZEN_METADATA).expect("TW-04 reference metadata JSON");
    assert_eq!(metadata["specification"], "TW-04");
    assert_eq!(metadata["reference_method"], "DOP853");
    assert_eq!(metadata["cross_check_method"], "Radau");
    assert!(metadata["python"].as_str().is_some());
    assert!(metadata["numpy"].as_str().is_some());
    assert!(metadata["scipy"].as_str().is_some());

    let rows = frozen_rows();
    for (scenario, gamma) in SCENARIOS {
        let selected = scenario_rows(&rows, scenario);
        assert!(
            selected.len() >= 1001,
            "{scenario} has dense fixed-grid coverage"
        );
        assert!(
            selected
                .iter()
                .all(|row| (row.gamma - gamma).abs() <= f64::EPSILON)
        );
        assert_eq!(selected.first().unwrap().z, 0.0);
        assert_eq!(selected.last().unwrap().z, 1.0);
    }
}

#[test]
fn tw04_production_rhs_matches_independent_rates_and_heat_at_control_states() {
    let rows = frozen_rows();
    for (scenario, gamma) in SCENARIOS {
        let selected = scenario_rows(&rows, scenario);
        let mut task = task_for(gamma, ReactorIvpSolverFamily::Lsode2);
        task.setup_condensed_ivp()
            .expect("TW-04 production equations assemble");
        assert_eq!(task.thermal_effects, reaction_heats(gamma));

        let mut observed = Vec::new();
        let heat_peak = selected
            .iter()
            .max_by(|left, right| left.q_volume.abs().total_cmp(&right.q_volume.abs()))
            .unwrap()
            .z;
        for target in [0.0, heat_peak, 0.4, 1.0] {
            observed.push(reference_state_at(&selected, target));
        }
        for state in observed {
            let expected_temperature = DT + T_SCALE * state[0];
            let (rates, heat) = independent_rates_and_heat(
                expected_temperature,
                &state[2..].try_into().unwrap(),
                gamma,
            );
            let mut solver_order_state = [0.0; 14];
            solver_order_state[..2].copy_from_slice(&state[..2]);
            for (solver_index, species_name) in task.kindata.substances.iter().enumerate() {
                let canonical = species_name
                    .strip_prefix('X')
                    .unwrap()
                    .parse::<usize>()
                    .unwrap();
                solver_order_state[solver_index + 2] = state[canonical + 2];
            }
            let variables = task
                .solver
                .unknowns
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>();
            let actual_heat = task
                .heat_release
                .eval_expression(&variables, &solver_order_state);
            assert!((actual_heat - heat).abs() <= 1e-10 + 1e-11 * heat.abs());
            for (reaction, expected_rate) in rates.iter().copied().enumerate() {
                let rate_constant = task
                    .kindata
                    .K_sym_vec
                    .as_ref()
                    .expect("reaction rate expressions")[reaction]
                    .eval_expression(&["Teta"], &[state[0]]);
                let reactant = REACTIONS[reaction].0;
                let reactant_state = task
                    .kindata
                    .substances
                    .iter()
                    .position(|name| name == &species_name(reactant))
                    .unwrap();
                let actual_rate = rate_constant * RHO * solver_order_state[reactant_state + 2]
                    / MOLAR_MASS_KG_PER_MOL;
                assert!(
                    (actual_rate - expected_rate).abs() <= 1e-9 + 1e-10 * expected_rate.abs(),
                    "{scenario} reaction={reaction} actual={actual_rate:.16e} expected={expected_rate:.16e}"
                );
            }
            let actual_rhs = evaluate_production_rhs(&task, &solver_order_state);
            let expected_theta_rhs = state[1] / LAMBDA;
            let expected_qhat_rhs =
                MASS_FLUX * CP * LENGTH / LAMBDA * state[1] - LENGTH.powi(2) / T_SCALE * heat;
            assert!((actual_rhs[0] - expected_theta_rhs).abs() <= 1e-12);
            assert!(
                (actual_rhs[1] - expected_qhat_rhs).abs()
                    <= 1e-11 + 1e-10 * expected_qhat_rhs.abs()
            );
            let mut expected_species_rhs = [0.0; SPECIES_COUNT];
            for (rate, (from, to, _)) in rates.iter().zip(REACTIONS) {
                expected_species_rhs[from] -= LENGTH / MASS_FLUX * rate * MOLAR_MASS_KG_PER_MOL;
                expected_species_rhs[to] += LENGTH / MASS_FLUX * rate * MOLAR_MASS_KG_PER_MOL;
            }
            for solver_index in 0..SPECIES_COUNT {
                let canonical = task.kindata.substances[solver_index]
                    .strip_prefix('X')
                    .unwrap()
                    .parse::<usize>()
                    .unwrap();
                assert!(
                    (actual_rhs[solver_index + 2] - expected_species_rhs[canonical]).abs()
                        <= 2e-10 + 2e-10 * expected_species_rhs[canonical].abs(),
                    "{scenario} Y{canonical} rhs actual={} expected={}",
                    actual_rhs[solver_index + 2],
                    expected_species_rhs[canonical]
                );
            }
        }
    }
}

#[test]
fn tw04_all_stiff_api_families_match_the_independent_frozen_profiles() {
    let rows = frozen_rows();
    for (scenario, gamma) in SCENARIOS {
        let expected = scenario_rows(&rows, scenario);
        for family in [
            ReactorIvpSolverFamily::Lsode2,
            ReactorIvpSolverFamily::Bdf,
            ReactorIvpSolverFamily::Radau,
        ] {
            let mut task = task_for(gamma, family);
            let snapshot = task.solve_condensed_ivp().unwrap_or_else(|error| {
                panic!("{scenario} {family:?} gamma={gamma:.0} solve failed: {error}")
            });
            assert_snapshot_matches_reference(
                scenario,
                family,
                gamma,
                &snapshot,
                &task.kindata.substances,
                &expected,
            );
        }
    }
}

#[test]
fn tw04_gamma_limits_and_nonlinear_response_are_observable() {
    let rows = frozen_rows();
    let zero = scenario_rows(&rows, "tw_04a");
    let weak = scenario_rows(&rows, "tw_04b");
    let moderate = scenario_rows(&rows, "tw_04c");
    let strong = scenario_rows(&rows, "tw_04d");
    let zero_end = zero.last().unwrap().state;
    let expected_end_theta = (REFERENCE_TEMPERATURE_K - DT) / T_SCALE
        + INITIAL_QHAT / LAMBDA * (MASS_FLUX * CP * LENGTH / LAMBDA).exp_m1()
            / (MASS_FLUX * CP * LENGTH / LAMBDA);
    assert!((zero_end[0] - expected_end_theta).abs() < 2e-12);
    let uniform = (0..=100)
        .map(|index| index as f64 / 100.0)
        .collect::<Vec<_>>();
    let moderate_temperature_shift = uniform
        .iter()
        .map(|z| {
            ((reference_state_at(&moderate, *z)[0] - reference_state_at(&zero, *z)[0]) * T_SCALE)
                .abs()
        })
        .fold(0.0_f64, f64::max);
    let weak_temperature_shift = uniform
        .iter()
        .map(|z| {
            ((reference_state_at(&weak, *z)[0] - reference_state_at(&zero, *z)[0]) * T_SCALE).abs()
        })
        .fold(0.0_f64, f64::max);
    let strong_temperature_shift = uniform
        .iter()
        .map(|z| {
            ((reference_state_at(&strong, *z)[0] - reference_state_at(&zero, *z)[0]) * T_SCALE)
                .abs()
        })
        .fold(0.0_f64, f64::max);
    let strong_composition_shift = uniform
        .iter()
        .flat_map(|z| {
            let strong_state = reference_state_at(&strong, *z);
            let zero_state = reference_state_at(&zero, *z);
            (2..14).map(move |state| (strong_state[state] - zero_state[state]).abs())
        })
        .fold(0.0_f64, f64::max);
    let strong_nonlinear_temperature_departure = uniform
        .iter()
        .map(|z| {
            let zero_t = reference_state_at(&zero, *z)[0] * T_SCALE;
            let weak_t = reference_state_at(&weak, *z)[0] * T_SCALE;
            let strong_t = reference_state_at(&strong, *z)[0] * T_SCALE;
            (strong_t - zero_t - STRONG_GAMMA / WEAK_GAMMA * (weak_t - zero_t)).abs()
        })
        .fold(0.0_f64, f64::max);
    assert!(weak_temperature_shift > 1e-3);
    assert!(moderate_temperature_shift > 1.0);
    assert!(strong_temperature_shift > 10.0);
    assert!(strong_composition_shift > 1e-2);
    assert!(strong_nonlinear_temperature_departure > 10.0);
    assert!(
        strong
            .iter()
            .all(|row| row.state[2..].iter().all(|fraction| *fraction > 0.0))
    );
    assert!(strong.iter().any(|row| row.q_volume.abs() > 0.0));

    let shifted_enthalpies = ENTHALPIES_J_PER_MOL.map(|enthalpy| enthalpy + 9876.5);
    assert_eq!(
        reaction_heats_from_enthalpy(&shifted_enthalpies),
        reaction_heats_from_enthalpy(&ENTHALPIES_J_PER_MOL)
    );
    for (gamma, _) in [
        (WEAK_GAMMA, "weak"),
        (MODERATE_GAMMA, "moderate"),
        (STRONG_GAMMA, "strong"),
    ] {
        let scaled = reaction_heats(gamma);
        let unit = reaction_heats(1.0);
        for (actual, base) in scaled.iter().zip(unit) {
            assert!((*actual - gamma * base).abs() <= 1e-10);
        }
    }

    let mut zero_activation = zero_heat_arrhenius_task(ReactorIvpSolverFamily::Lsode2);
    zero_activation
        .setup_condensed_ivp()
        .expect("zero-Ea production system assembles");
    for theta in [0.0, 1.0, 4.0] {
        for (reaction, (_, _, k_reference)) in REACTIONS.iter().enumerate() {
            let k = zero_activation.kindata.K_sym_vec.as_ref().unwrap()[reaction]
                .eval_expression(&["Teta"], &[theta]);
            assert!((k - k_reference).abs() <= 1e-14 * k_reference.abs());
        }
    }

    let mut tw03_lines = TW03_FROZEN_PROFILE.lines();
    assert!(
        tw03_lines
            .next()
            .unwrap()
            .starts_with("scenario,z,Teta,qhat,T_K,")
    );
    let mut tw03c_rows = Vec::new();
    for line in tw03_lines {
        let fields = line.split(',').collect::<Vec<_>>();
        if fields[0] != "TW-03c" {
            continue;
        }
        let values = fields[1..]
            .iter()
            .map(|value| value.parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        tw03c_rows.push(values);
    }
    assert_eq!(tw03c_rows.len(), 11);
    for row in tw03c_rows {
        let z = row[0];
        let tw04_state = reference_state_at(&zero, z);
        assert!(
            (tw04_state[0] - row[1]).abs() <= 2e-9,
            "TW-04a/TW-03c Teta z={z}"
        );
        assert!(
            (tw04_state[1] - row[2]).abs() <= 2e-6,
            "TW-04a/TW-03c qhat z={z}"
        );
        for species in 0..SPECIES_COUNT {
            assert!(
                (tw04_state[species + 2] - row[4 + species]).abs() <= 2e-8,
                "TW-04a/TW-03c Y{species} z={z}"
            );
        }
    }
}
