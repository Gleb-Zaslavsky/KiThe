use super::matrix_exponential::expm_action;
use super::tw_01::{
    CP, DT, INITIAL_MASS_FRACTIONS, LAMBDA, LENGTH, MASS_FLUX, REACTIONS, RHO, SPECIES_COUNT,
    T_SCALE, species_name, tw01_task,
};
use crate::Kinetics::mechfinder_api::ReactionData;
use crate::SimpleReactorsIVP::SimpleReactorIVP::{IvpSolveSnapshot, SimpleReactorTask};
use crate::SimpleReactorsIVP::solver_backend::ReactorIvpSolverFamily;
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;

pub(super) const FROZEN_PROFILE: &str = include_str!("fixtures/tw_03_profiles.csv");
const GAS_CONSTANT_J_PER_MOL_K: f64 = 8.314;
const REFERENCE_TEMPERATURE_K: f64 = 400.0;
const VARIABLE_QHAT: f64 = 5_000.0;
const PROFILE_GRID: [f64; 11] = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
pub(super) const ACTIVATION_ENERGIES_J_PER_MOL: [f64; 20] = [
    20_000.0, 45_000.0, 30_000.0, 65_000.0, 40_000.0, 55_000.0, 25_000.0, 70_000.0, 35_000.0,
    60_000.0, 25_000.0, 70_000.0, 55_000.0, 35_000.0, 75_000.0, 20_000.0, 60_000.0, 45_000.0,
    80_000.0, 30_000.0,
];

#[derive(Clone, Copy)]
enum Scenario {
    IsothermalDifferentEa,
    VariableEqualEa,
    VariableDifferentEa,
}

impl Scenario {
    fn activation_energies(self) -> [f64; 20] {
        match self {
            Self::IsothermalDifferentEa | Self::VariableDifferentEa => {
                ACTIVATION_ENERGIES_J_PER_MOL
            }
            Self::VariableEqualEa => [40_000.0; 20],
        }
    }

    fn initial_qhat(self) -> f64 {
        match self {
            Self::IsothermalDifferentEa => 0.0,
            Self::VariableEqualEa | Self::VariableDifferentEa => VARIABLE_QHAT,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::IsothermalDifferentEa => "TW-03a",
            Self::VariableEqualEa => "TW-03b",
            Self::VariableDifferentEa => "TW-03c",
        }
    }
}

fn reference_rates() -> [f64; 20] {
    REACTIONS.map(|(_, _, rate)| rate)
}

pub(super) fn arrhenius_pre_exponentials(activation_energies: &[f64; 20]) -> [f64; 20] {
    let rates = reference_rates();
    std::array::from_fn(|reaction| {
        rates[reaction]
            * (activation_energies[reaction] / (GAS_CONSTANT_J_PER_MOL_K * REFERENCE_TEMPERATURE_K))
                .exp()
    })
}

fn rate_constants(temperature_k: f64, activation_energies: &[f64; 20]) -> [f64; 20] {
    let pre_exponentials = arrhenius_pre_exponentials(activation_energies);
    std::array::from_fn(|reaction| {
        pre_exponentials[reaction]
            * (-activation_energies[reaction] / (GAS_CONSTANT_J_PER_MOL_K * temperature_k)).exp()
    })
}

fn temperature_field(z: f64, qhat_0: f64) -> (f64, f64, f64) {
    let pe_q = LENGTH * MASS_FLUX * CP / LAMBDA;
    let qhat = qhat_0 * (pe_q * z).exp();
    let theta_0 = (REFERENCE_TEMPERATURE_K - DT) / T_SCALE;
    let theta = theta_0 + qhat_0 / LAMBDA * (pe_q * z).exp_m1() / pe_q;
    let temperature_k = DT + T_SCALE * theta;
    (theta, qhat, temperature_k)
}

fn generator_for_rates(rates: &[f64; 20]) -> DMatrix<f64> {
    let scale = LENGTH * RHO / MASS_FLUX;
    let mut matrix = DMatrix::zeros(SPECIES_COUNT, SPECIES_COUNT);
    for (reaction, (from, to, _)) in REACTIONS.iter().copied().enumerate() {
        matrix[(from, from)] -= scale * rates[reaction];
        matrix[(to, from)] += scale * rates[reaction];
    }
    matrix
}

fn equal_activation_reference(z: f64, activation_energy: f64) -> DVector<f64> {
    let initial_temperature = REFERENCE_TEMPERATURE_K;
    let relative_rate = |position: f64| {
        let temperature = temperature_field(position, VARIABLE_QHAT).2;
        (activation_energy / GAS_CONSTANT_J_PER_MOL_K
            * (1.0 / initial_temperature - 1.0 / temperature))
            .exp()
    };
    let sigma = composite_simpson(relative_rate, z, 1024);
    expm_action(
        &(sigma * generator_for_rates(&reference_rates())),
        &DVector::from_column_slice(&INITIAL_MASS_FRACTIONS),
    )
}

fn composite_simpson(function: impl Fn(f64) -> f64, end: f64, intervals: usize) -> f64 {
    if end == 0.0 {
        return 0.0;
    }
    assert!(intervals >= 2 && intervals.is_multiple_of(2));
    let step = end / intervals as f64;
    let mut sum = function(0.0) + function(end);
    for index in 1..intervals {
        sum += if index.is_multiple_of(2) { 2.0 } else { 4.0 } * function(index as f64 * step);
    }
    sum * step / 3.0
}

fn commutator(left: &DMatrix<f64>, right: &DMatrix<f64>) -> DMatrix<f64> {
    left * right - right * left
}

fn magnus_step(
    state: &DVector<f64>,
    start: f64,
    step: f64,
    activation_energies: &[f64; 20],
) -> DVector<f64> {
    let root_three = 3.0_f64.sqrt();
    let first_position = start + step * (0.5 - root_three / 6.0);
    let second_position = start + step * (0.5 + root_three / 6.0);
    let first_rates = rate_constants(
        temperature_field(first_position, VARIABLE_QHAT).2,
        activation_energies,
    );
    let second_rates = rate_constants(
        temperature_field(second_position, VARIABLE_QHAT).2,
        activation_energies,
    );
    let first = generator_for_rates(&first_rates);
    let second = generator_for_rates(&second_rates);
    let omega = (step / 2.0) * (&first + &second)
        - (root_three * step.powi(2) / 12.0) * commutator(&first, &second);
    expm_action(&omega, state)
}

fn magnus_profile(z: f64, activation_energies: &[f64; 20], steps_per_unit: usize) -> DVector<f64> {
    let steps = (z * steps_per_unit as f64).ceil() as usize;
    if steps == 0 {
        return DVector::from_column_slice(&INITIAL_MASS_FRACTIONS);
    }
    let step = z / steps as f64;
    let mut state = DVector::from_column_slice(&INITIAL_MASS_FRACTIONS);
    for index in 0..steps {
        state = magnus_step(&state, index as f64 * step, step, activation_energies);
    }
    state
}

fn magnus_trajectory(
    axis: &[f64],
    activation_energies: &[f64; 20],
    steps_per_unit: usize,
) -> Vec<DVector<f64>> {
    let mut state = DVector::from_column_slice(&INITIAL_MASS_FRACTIONS);
    let mut previous = 0.0;
    let mut trajectory = Vec::with_capacity(axis.len());
    for &target in axis {
        assert!(target >= previous && target <= 1.0 + 1e-12);
        let interval = target - previous;
        let steps = (interval * steps_per_unit as f64).ceil() as usize;
        if steps > 0 {
            let step = interval / steps as f64;
            for index in 0..steps {
                state = magnus_step(
                    &state,
                    previous + index as f64 * step,
                    step,
                    activation_energies,
                );
            }
        }
        trajectory.push(state.clone());
        previous = target;
    }
    trajectory
}

fn independent_profile(scenario: Scenario, z: f64) -> DVector<f64> {
    let energies = scenario.activation_energies();
    match scenario {
        Scenario::IsothermalDifferentEa => expm_action(
            &(z * generator_for_rates(&reference_rates())),
            &DVector::from_column_slice(&INITIAL_MASS_FRACTIONS),
        ),
        Scenario::VariableEqualEa => equal_activation_reference(z, energies[0]),
        Scenario::VariableDifferentEa => magnus_profile(z, &energies, 2560),
    }
}

fn task_for(scenario: Scenario, family: ReactorIvpSolverFamily) -> SimpleReactorTask {
    task_for_parameters(
        scenario.activation_energies(),
        scenario.initial_qhat(),
        family,
    )
}

pub(super) fn task_for_parameters(
    energies: [f64; 20],
    initial_qhat: f64,
    family: ReactorIvpSolverFamily,
) -> SimpleReactorTask {
    let mut task = tw01_task(family);
    let pre_exponentials = arrhenius_pre_exponentials(&energies);
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
        .enumerate()
        .map(|(index, (from, to, _))| {
            ReactionData::new_elementary(
                format!("{}=>{}", species_name(*from), species_name(*to)),
                vec![pre_exponentials[index], 0.0, energies[index]],
                None,
            )
        })
        .collect();
    task.kindata
        .set_reaction_data_directly(reactions, Some(groups))
        .expect("TW-03 Arrhenius reactions should pass the production data boundary");
    task.thermal_effects = vec![0.0; REACTIONS.len()];
    task.initial_conditions
        .insert("q".to_string(), initial_qhat * T_SCALE / LENGTH);
    task.solver_backend_config = task.solver_backend_config.with_rtol(1e-11).with_atol(1e-13);
    task
}

fn state_column(snapshot: &IvpSolveSnapshot, name: &str) -> usize {
    snapshot
        .variable_names
        .iter()
        .position(|candidate| candidate == name)
        .unwrap_or_else(|| panic!("missing state {name}"))
}

fn species_index(state_name: &str, species_order: &[String]) -> usize {
    let state_index = state_name
        .strip_prefix('C')
        .and_then(|value| value.parse::<usize>().ok())
        .expect("species state name C<n>");
    species_order[state_index]
        .strip_prefix('X')
        .and_then(|value| value.parse::<usize>().ok())
        .expect("canonical TW-03 species name X<n>")
}

fn assert_solver_matches_reference(
    scenario: Scenario,
    family: ReactorIvpSolverFamily,
    snapshot: &IvpSolveSnapshot,
    species_order: &[String],
) {
    let theta_column = state_column(snapshot, "Teta");
    let q_column = state_column(snapshot, "q");
    let mut max_species_error: f64 = 0.0;
    let mut max_temperature_error: f64 = 0.0;
    let mut max_mass_error: f64 = 0.0;
    let mut max_element_error: f64 = 0.0;
    let mut minimum_fraction = f64::INFINITY;

    let magnus = if matches!(scenario, Scenario::VariableDifferentEa) {
        Some(magnus_trajectory(
            snapshot.axis.as_slice(),
            &scenario.activation_energies(),
            2560,
        ))
    } else {
        None
    };
    for (row, z) in snapshot.axis.iter().copied().enumerate() {
        let (theta, qhat, temperature) = temperature_field(z, scenario.initial_qhat());
        let composition = match scenario {
            Scenario::VariableDifferentEa => magnus.as_ref().unwrap()[row].clone(),
            _ => independent_profile(scenario, z),
        };
        let actual_temperature = DT + T_SCALE * snapshot.values[(row, theta_column)];
        let theta_error = (snapshot.values[(row, theta_column)] - theta).abs();
        let theta_tolerance = 1e-9 + 1e-7 * theta.abs();
        assert!(
            theta_error <= theta_tolerance,
            "{} {family:?} Teta z={z}: actual={:.16e} expected={theta:.16e}",
            scenario.name(),
            snapshot.values[(row, theta_column)]
        );
        let qhat_error = (snapshot.values[(row, q_column)] - qhat).abs();
        let qhat_tolerance = 1e-9 + 1e-7 * qhat.abs();
        assert!(
            qhat_error <= qhat_tolerance,
            "{} {family:?} qhat z={z}: actual={:.16e} expected={qhat:.16e}",
            scenario.name(),
            snapshot.values[(row, q_column)]
        );
        assert!(
            (actual_temperature - temperature).abs() <= 1e-6,
            "{} {family:?} T z={z}: actual={actual_temperature:.16e} expected={temperature:.16e}",
            scenario.name()
        );

        let mut mass_sum = 0.0;
        for (column, name) in snapshot.variable_names.iter().enumerate() {
            if !name.starts_with('C') {
                continue;
            }
            let canonical = species_index(name, species_order);
            let actual = snapshot.values[(row, column)];
            let expected = composition[canonical];
            assert!(
                actual.is_finite() && actual > 0.0,
                "{name} at z={z}: {actual}"
            );
            assert!(expected.is_finite() && expected > 0.0);
            let error = (actual - expected).abs();
            let normalized_error = error / (1e-9 + 1e-7 * expected.abs());
            assert!(
                normalized_error <= 1.0,
                "{} {family:?} {name} z={z}: actual={actual:.16e} expected={expected:.16e} error={error:.3e} normalized={normalized_error:.3e}",
                scenario.name()
            );
            max_species_error = max_species_error.max(error);
            minimum_fraction = minimum_fraction.min(actual);
            mass_sum += actual;
        }
        let mass_error = (mass_sum - 1.0).abs();
        assert!(mass_error <= 2e-9, "mass balance at z={z}: {mass_sum}");
        max_mass_error = max_mass_error.max(mass_error);
        for element_mass_fraction in [24.022 / 60.052, 4.032 / 60.052, 31.998 / 60.052] {
            let error = (mass_sum * element_mass_fraction - element_mass_fraction).abs();
            assert!(error <= 2e-9, "element balance at z={z}: {error:.3e}");
            max_element_error = max_element_error.max(error);
        }
        max_temperature_error = max_temperature_error.max((actual_temperature - temperature).abs());
    }
    println!(
        "{} family={family:?} points={} max_dY={max_species_error:.3e} max_dT_K={max_temperature_error:.3e} max_mass={max_mass_error:.3e} max_element={max_element_error:.3e} min_Y={minimum_fraction:.6e}",
        scenario.name(),
        PROFILE_GRID.len()
    );
}

#[test]
fn tw03_arrhenius_contract_uses_joules_per_mole_absolute_temperature_and_reference_rates() {
    for scenario in [
        Scenario::IsothermalDifferentEa,
        Scenario::VariableEqualEa,
        Scenario::VariableDifferentEa,
    ] {
        let energies = scenario.activation_energies();
        let mut task = task_for(scenario, ReactorIvpSolverFamily::Lsode2);
        task.setup_condensed_ivp()
            .expect("TW-03 should prepare through the ordinary KiThe IVP path");
        let symbolic = task.kindata.K_sym_vec.as_ref().unwrap();
        let variable_names = ["Teta"];
        for temperature in [400.0, 425.0, 450.0] {
            let rates = rate_constants(temperature, &energies);
            let theta = (temperature - DT) / T_SCALE;
            for (reaction, rate) in rates.into_iter().enumerate() {
                assert!(rate.is_finite() && rate > 0.0);
                if temperature == REFERENCE_TEMPERATURE_K {
                    assert!((rate - REACTIONS[reaction].2).abs() <= 1e-13);
                }
                let production_rate = symbolic[reaction].eval_expression(&variable_names, &[theta]);
                assert!(
                    (production_rate - rate).abs() <= 2e-12 * rate.abs(),
                    "{} reaction={} T={temperature} production={production_rate:.16e} oracle={rate:.16e}",
                    scenario.name(),
                    reaction + 1
                );
                let step = 1e-3;
                let plus = rate_constants(temperature + step, &energies)[reaction].ln();
                let minus = rate_constants(temperature - step, &energies)[reaction].ln();
                let finite_difference = (plus - minus) / (2.0 * step);
                let expected_derivative =
                    energies[reaction] / (GAS_CONSTANT_J_PER_MOL_K * temperature.powi(2));
                assert!((finite_difference - expected_derivative).abs() <= 2e-9);
            }
        }
    }

    let energies = Scenario::VariableDifferentEa.activation_energies();
    let pre_exponentials = arrhenius_pre_exponentials(&energies);
    println!("TW-03 different-Ea Arrhenius parameters (reaction, Ea J/mol, A 1/s):");
    for reaction in 0..REACTIONS.len() {
        println!(
            "R{:02}, {:.0}, {:.12e}",
            reaction + 1,
            energies[reaction],
            pre_exponentials[reaction]
        );
    }

    let zero_energy = [0.0; 20];
    let mut zero_energy_task =
        task_for_parameters(zero_energy, VARIABLE_QHAT, ReactorIvpSolverFamily::Lsode2);
    zero_energy_task
        .setup_condensed_ivp()
        .expect("zero-Ea control should assemble through the ordinary path");
    let zero_energy_constants = zero_energy_task.kindata.K_sym_vec.as_ref().unwrap();
    for temperature in [400.0, 450.0] {
        let theta = (temperature - DT) / T_SCALE;
        for (reaction, expression) in zero_energy_constants.iter().enumerate() {
            let actual = expression.eval_expression(&["Teta"], &[theta]);
            assert!((actual - REACTIONS[reaction].2).abs() <= 1e-14);
        }
    }
}

#[test]
fn tw03_temperature_field_matches_the_current_scaled_thermal_rhs() {
    let pe_q = LENGTH * MASS_FLUX * CP / LAMBDA;
    assert!((pe_q - 0.001).abs() <= 1e-15);
    for z in PROFILE_GRID {
        let (_, qhat, temperature) = temperature_field(z, VARIABLE_QHAT);
        let exact_temperature = 400.0 + 50.0 * (pe_q * z).exp_m1() / pe_q;
        assert!((temperature - exact_temperature).abs() <= 1e-12);
        assert!((qhat - VARIABLE_QHAT * (pe_q * z).exp()).abs() <= 1e-12);
    }
    let end = temperature_field(1.0, VARIABLE_QHAT).2;
    assert!((end - 450.025008335417).abs() <= 2e-10);
}

#[test]
fn tw03_equal_activation_reference_quadrature_is_stable_and_commuting() {
    let activation_energy = 40_000.0;
    let relative_rate = |z: f64| {
        let temperature = temperature_field(z, VARIABLE_QHAT).2;
        (activation_energy / GAS_CONSTANT_J_PER_MOL_K
            * (1.0 / REFERENCE_TEMPERATURE_K - 1.0 / temperature))
            .exp()
    };
    let coarse = composite_simpson(relative_rate, 1.0, 256);
    let fine = composite_simpson(relative_rate, 1.0, 1024);
    println!(
        "TW-03b quadrature N=256/1024 sigma={coarse:.16e}/{fine:.16e} delta={:.3e}",
        (fine - coarse).abs()
    );
    assert!((fine - coarse).abs() <= 2e-11);

    let first = generator_for_rates(&rate_constants(
        temperature_field(0.1, VARIABLE_QHAT).2,
        &[activation_energy; 20],
    ));
    let second = generator_for_rates(&rate_constants(
        temperature_field(0.9, VARIABLE_QHAT).2,
        &[activation_energy; 20],
    ));
    let scaled_commutator =
        commutator(&first, &second).norm() / (first.norm() * second.norm()).max(f64::MIN_POSITIVE);
    assert!(
        scaled_commutator <= 2e-16,
        "commutator={scaled_commutator:.3e}"
    );

    let energies = [activation_energy; 20];
    for z in [0.1, 0.5, 1.0] {
        let magnus = magnus_profile(z, &energies, 2560);
        let scalar_quadrature = equal_activation_reference(z, activation_energy);
        assert!(
            (&magnus - scalar_quadrature).norm() <= 2e-10,
            "TW-03b Magnus/quadrature mismatch at z={z}"
        );
    }
}

#[test]
fn tw03_magnus_fourth_order_reference_passes_exact_and_convergence_controls() {
    let constant = DMatrix::from_row_slice(2, 2, &[-0.7, 0.0, 0.7, 0.0]);
    let initial = DVector::from_vec(vec![1.0, 0.0]);
    let exact = expm_action(&constant, &initial);
    let mut actual = initial.clone();
    for _ in 0..20 {
        actual = magnus_matrix_step(&actual, &constant, &constant, 0.05);
    }
    assert!((&actual - exact).norm() <= 2e-13);

    let scalar_a = DMatrix::from_element(1, 1, 0.3);
    let scalar_b = DMatrix::from_element(1, 1, -0.2);
    let scalar_start = DVector::from_element(1, 2.0);
    let scalar_result = magnus_matrix_step(&scalar_start, &scalar_a, &scalar_b, 0.4);
    assert!((scalar_result[0] - 2.0 * 0.02_f64.exp()).abs() <= 2e-14);

    let energies = Scenario::VariableDifferentEa.activation_energies();
    let n640 = magnus_profile(1.0, &energies, 640);
    let n1280 = magnus_profile(1.0, &energies, 1280);
    let n2560 = magnus_profile(1.0, &energies, 2560);
    let n5120 = magnus_profile(1.0, &energies, 5120);
    let error_640_1280 = (&n640 - &n1280).norm();
    let error_1280_2560 = (&n1280 - &n2560).norm();
    let reference_refinement_error = (&n2560 - &n5120).norm();
    let observed_order = (error_640_1280 / error_1280_2560).log2();
    println!(
        "TW-03 Magnus convergence N=640/1280/2560/5120 errors={error_640_1280:.3e}/{error_1280_2560:.3e}/{reference_refinement_error:.3e} p={observed_order:.3}"
    );
    assert!((observed_order - 4.0).abs() <= 0.12, "p={observed_order}");
    assert!(reference_refinement_error <= 2e-10);

    let first = generator_for_rates(&rate_constants(
        temperature_field(0.0, VARIABLE_QHAT).2,
        &energies,
    ));
    let last = generator_for_rates(&rate_constants(
        temperature_field(1.0, VARIABLE_QHAT).2,
        &energies,
    ));
    let normalized_commutator =
        commutator(&first, &last).norm() / (first.norm() * last.norm()).max(f64::MIN_POSITIVE);
    assert!(
        normalized_commutator > 1e-5,
        "commutator={normalized_commutator:.3e}"
    );
}

fn magnus_matrix_step(
    state: &DVector<f64>,
    first: &DMatrix<f64>,
    second: &DMatrix<f64>,
    step: f64,
) -> DVector<f64> {
    let omega = (step / 2.0) * (first + second)
        - (3.0_f64.sqrt() * step.powi(2) / 12.0) * commutator(first, second);
    expm_action(&omega, state)
}

#[test]
fn tw03_production_rhs_matches_independent_reaction_matrix_at_control_states() {
    for scenario in [
        Scenario::IsothermalDifferentEa,
        Scenario::VariableEqualEa,
        Scenario::VariableDifferentEa,
    ] {
        let mut task = task_for(scenario, ReactorIvpSolverFamily::Lsode2);
        task.setup_condensed_ivp()
            .expect("TW-03 production symbolic problem should assemble");
        let elemental_matrix = task
            .kindata
            .stecheodata
            .matrix_of_elements
            .as_ref()
            .expect("TW-03 element composition matrix");
        for (reaction, stoichiometry) in task.kindata.stecheodata.stecheo_matrx.iter().enumerate() {
            for element in 0..elemental_matrix.ncols() {
                let elemental_delta = (0..SPECIES_COUNT)
                    .map(|species| elemental_matrix[(species, element)] * stoichiometry[species])
                    .sum::<f64>();
                assert_eq!(
                    elemental_delta,
                    0.0,
                    "{} reaction={reaction} element={element}",
                    scenario.name()
                );
            }
        }
        let energies = scenario.activation_energies();
        let variables = task
            .solver
            .unknowns
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        for z in [0.0, 0.4, 1.0] {
            let composition = independent_profile(scenario, z);
            let (theta, qhat, temperature) = temperature_field(z, scenario.initial_qhat());
            let mut state = vec![0.0; SPECIES_COUNT + 2];
            state[0] = theta;
            state[1] = qhat;
            for (index, name) in task.kindata.substances.iter().enumerate() {
                let canonical = name.strip_prefix('X').unwrap().parse::<usize>().unwrap();
                state[index + 2] = composition[canonical];
            }
            let rates = rate_constants(temperature, &energies);
            let expected_matrix = generator_for_rates(&rates);
            for (state_index, name) in task.kindata.substances.iter().enumerate() {
                let canonical_row = name.strip_prefix('X').unwrap().parse::<usize>().unwrap();
                let expected = (0..SPECIES_COUNT)
                    .map(|source| {
                        let source_name = &task.kindata.substances[source];
                        let canonical_column = source_name
                            .strip_prefix('X')
                            .unwrap()
                            .parse::<usize>()
                            .unwrap();
                        expected_matrix[(canonical_row, canonical_column)] * state[source + 2]
                    })
                    .sum::<f64>();
                let actual =
                    task.solver.eq_system[state_index + 2].eval_expression(&variables, &state);
                assert!(
                    (actual - expected).abs() <= 2e-10 + 2e-11 * expected.abs(),
                    "{} z={z} species={name} actual={actual:.16e} expected={expected:.16e}",
                    scenario.name()
                );
            }
            assert_eq!(task.heat_release.eval_expression(&variables, &state), 0.0);
        }
    }
}

#[test]
fn tw03_all_stiff_api_families_match_the_independent_scenario_references() {
    for scenario in [
        Scenario::IsothermalDifferentEa,
        Scenario::VariableEqualEa,
        Scenario::VariableDifferentEa,
    ] {
        for family in [
            ReactorIvpSolverFamily::Lsode2,
            ReactorIvpSolverFamily::Bdf,
            ReactorIvpSolverFamily::Radau,
        ] {
            let mut task = task_for(scenario, family);
            let snapshot = task.solve_condensed_ivp().unwrap_or_else(|error| {
                panic!("{} {family:?} solve failed: {error}", scenario.name())
            });
            assert_solver_matches_reference(scenario, family, &snapshot, &task.kindata.substances);
        }
    }
}

#[test]
fn tw03_frozen_profiles_match_independent_analytic_and_magnus_references() {
    let scenarios = [
        Scenario::IsothermalDifferentEa,
        Scenario::VariableEqualEa,
        Scenario::VariableDifferentEa,
    ];
    let mut lines = FROZEN_PROFILE.lines();
    assert_eq!(
        lines.next(),
        Some("scenario,z,Teta,qhat,T_K,X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11")
    );
    let mut row_count = 0;
    for scenario in scenarios {
        for expected_z in PROFILE_GRID {
            let values = lines
                .next()
                .expect("TW-03 frozen profile row")
                .split(',')
                .map(str::trim)
                .collect::<Vec<_>>();
            assert_eq!(values.len(), SPECIES_COUNT + 5);
            assert_eq!(values[0], scenario.name());
            let z = values[1].parse::<f64>().expect("frozen coordinate");
            assert!((z - expected_z).abs() <= f64::EPSILON);
            let (theta, qhat, temperature) = temperature_field(z, scenario.initial_qhat());
            for (column, expected) in [(2, theta), (3, qhat), (4, temperature)] {
                let actual = values[column].parse::<f64>().expect("frozen thermal value");
                assert!(
                    (actual - expected).abs() <= 2e-13,
                    "{} z={z} frozen column={column} actual={actual:.16e} expected={expected:.16e}",
                    scenario.name()
                );
            }
            let composition = independent_profile(scenario, z);
            for species in 0..SPECIES_COUNT {
                let actual = values[species + 5]
                    .parse::<f64>()
                    .expect("frozen species fraction");
                assert!(
                    (actual - composition[species]).abs() <= 2e-13,
                    "{} z={z} X{species} frozen={actual:.16e} oracle={:.16e}",
                    scenario.name(),
                    composition[species]
                );
            }
            row_count += 1;
        }
    }
    assert!(lines.next().is_none());
    assert_eq!(row_count, 3 * PROFILE_GRID.len());
}
