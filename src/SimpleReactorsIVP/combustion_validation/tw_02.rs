use super::matrix_exponential::expm_action;
use super::tw_01::{
    CP, DT, INITIAL_MASS_FRACTIONS, INITIAL_Q, LAMBDA, LENGTH, MASS_FLUX, REACTIONS, RHO,
    SPECIES_COUNT, T_SCALE, analytical_profile as tw01_profile, generator as tw01_generator,
    tw01_task,
};
use crate::SimpleReactorsIVP::SimpleReactorIVP::{IvpSolveSnapshot, SimpleReactorTask};
use crate::SimpleReactorsIVP::solver_backend::ReactorIvpSolverFamily;
use nalgebra::{DMatrix, DVector};

const FROZEN_PROFILE: &str = include_str!("fixtures/tw_02_profiles.csv");
const MOLAR_MASS_KG_PER_MOL: f64 = 0.060_052;
pub(super) const ENTHALPIES_J_PER_MOL: [f64; SPECIES_COUNT] = [
    0.0, -40.0, -25.0, -65.0, -15.0, -80.0, -35.0, -100.0, -50.0, -90.0, -60.0, -115.0,
];
const GAMMAS: [f64; 3] = [1.0, 10.0, -1.0];
const PROFILE_GRID: [f64; 11] = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

#[derive(Debug)]
struct Tw02Model {
    reaction_heats_j_per_mol: Vec<f64>,
    reaction_rate_map: DMatrix<f64>,
    full_system: DMatrix<f64>,
    initial_state: DVector<f64>,
}

impl Tw02Model {
    fn new(gamma: f64) -> Self {
        let reaction_heats_j_per_mol = REACTIONS
            .iter()
            .map(|(from, to, _)| ENTHALPIES_J_PER_MOL[*from] - ENTHALPIES_J_PER_MOL[*to])
            .collect::<Vec<_>>();
        let reaction_rate_map = reaction_rate_map();
        let beta = DVector::from_iterator(
            SPECIES_COUNT,
            (0..SPECIES_COUNT).map(|species| {
                -LENGTH.powi(2) / T_SCALE
                    * gamma
                    * (0..REACTIONS.len())
                        .map(|reaction| {
                            reaction_heats_j_per_mol[reaction]
                                * reaction_rate_map[(reaction, species)]
                        })
                        .sum::<f64>()
            }),
        );
        let chemistry = tw01_generator();
        let pe_q = LENGTH * MASS_FLUX * CP / LAMBDA;
        let mut full_system = DMatrix::zeros(SPECIES_COUNT + 2, SPECIES_COUNT + 2);
        full_system[(0, 1)] = 1.0 / LAMBDA;
        full_system[(1, 1)] = pe_q;
        for species in 0..SPECIES_COUNT {
            full_system[(1, species + 2)] = beta[species];
            for source in 0..SPECIES_COUNT {
                full_system[(species + 2, source + 2)] = chemistry[(species, source)];
            }
        }

        let mut initial_state = DVector::zeros(SPECIES_COUNT + 2);
        initial_state[0] = (400.0 - DT) / T_SCALE;
        initial_state[1] = INITIAL_Q;
        for species in 0..SPECIES_COUNT {
            initial_state[species + 2] = INITIAL_MASS_FRACTIONS[species];
        }
        Self {
            reaction_heats_j_per_mol,
            reaction_rate_map,
            full_system,
            initial_state,
        }
    }

    fn state(&self, z: f64) -> DVector<f64> {
        expm_action(&(z * &self.full_system), &self.initial_state)
    }

    fn chemical_heat_release(&self, composition: &DVector<f64>, gamma: f64) -> f64 {
        (0..REACTIONS.len())
            .map(|reaction| {
                let omega = (0..SPECIES_COUNT)
                    .map(|species| {
                        self.reaction_rate_map[(reaction, species)] * composition[species]
                    })
                    .sum::<f64>();
                gamma * self.reaction_heats_j_per_mol[reaction] * omega
            })
            .sum()
    }
}

fn reaction_rate_map() -> DMatrix<f64> {
    let mut map = DMatrix::zeros(REACTIONS.len(), SPECIES_COUNT);
    for (reaction, (reactant, _, pre_exponential)) in REACTIONS.iter().copied().enumerate() {
        map[(reaction, reactant)] = RHO * pre_exponential / MOLAR_MASS_KG_PER_MOL;
    }
    map
}

fn task_for(family: ReactorIvpSolverFamily, gamma: f64, model: &Tw02Model) -> SimpleReactorTask {
    let mut task = tw01_task(family);
    task.thermal_effects = model
        .reaction_heats_j_per_mol
        .iter()
        .map(|heat| gamma * heat)
        .collect();
    task
}

fn snapshot_species_index(state_name: &str, species_order: &[String]) -> usize {
    let state_index = state_name
        .strip_prefix('C')
        .and_then(|value| value.parse::<usize>().ok())
        .expect("concentration variable name C<n>");
    species_order[state_index]
        .strip_prefix('X')
        .and_then(|value| value.parse::<usize>().ok())
        .expect("TW-02 species name X<n>")
}

fn snapshot_species(
    snapshot: &IvpSolveSnapshot,
    species_order: &[String],
    row: usize,
) -> DVector<f64> {
    let mut composition = DVector::zeros(SPECIES_COUNT);
    for (column, state_name) in snapshot.variable_names.iter().enumerate() {
        if state_name.starts_with('C') {
            composition[snapshot_species_index(state_name, species_order)] =
                snapshot.values[(row, column)];
        }
    }
    composition
}

fn assert_profile(
    snapshot: &IvpSolveSnapshot,
    task: &SimpleReactorTask,
    model: &Tw02Model,
    gamma: f64,
    family: ReactorIvpSolverFamily,
) {
    let theta_column = snapshot
        .variable_names
        .iter()
        .position(|name| name == "Teta")
        .unwrap();
    let q_column = snapshot
        .variable_names
        .iter()
        .position(|name| name == "q")
        .unwrap();
    let variable_names = snapshot
        .variable_names
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    let mut errors = [0.0_f64; 5];
    let mut max_element_balance: f64 = 0.0;
    let mut max_energy_rhs_error: f64 = 0.0;
    let mut max_heat_source_error: f64 = 0.0;
    let mut min_temperature = f64::INFINITY;
    let mut max_temperature = f64::NEG_INFINITY;

    assert!(snapshot.status.starts_with("finished"));
    assert_eq!(snapshot.variable_count(), SPECIES_COUNT + 2);
    for (row, z) in snapshot.axis.iter().copied().enumerate() {
        let expected = model.state(z);
        let actual_composition = snapshot_species(snapshot, &task.kindata.substances, row);
        let expected_composition = expected.rows(2, SPECIES_COUNT).into_owned();
        for species in 0..SPECIES_COUNT {
            let error = (actual_composition[species] - expected_composition[species]).abs();
            errors[0] = errors[0].max(error);
            let normalized = error / (1e-9 + 1e-7 * expected_composition[species].abs());
            assert!(
                normalized <= 1.0,
                "{family:?}/gamma={gamma} X{species} z={z}: actual={} reference={} abs={error:.3e} normalized={normalized:.3e}",
                actual_composition[species],
                expected_composition[species]
            );
        }
        let actual_theta = snapshot.values[(row, theta_column)];
        let actual_q = snapshot.values[(row, q_column)];
        let theta_error = (actual_theta - expected[0]).abs();
        let q_error = (actual_q - expected[1]).abs();
        errors[1] = errors[1].max(theta_error);
        errors[2] = errors[2].max(q_error);
        assert!(
            theta_error / (1e-9 + 1e-7 * expected[0].abs()) <= 1.0,
            "{family:?}/gamma={gamma} Teta z={z} actual={actual_theta:.16e} reference={:.16e}",
            expected[0]
        );
        assert!(
            q_error / (1e-9 + 1e-7 * expected[1].abs()) <= 1.0,
            "{family:?}/gamma={gamma} q z={z} actual={actual_q:.16e} reference={:.16e}",
            expected[1]
        );
        let temperature = DT + T_SCALE * actual_theta;
        let reference_temperature = DT + T_SCALE * expected[0];
        let temperature_error = (temperature - reference_temperature).abs();
        errors[3] = errors[3].max(temperature_error);
        assert!(
            temperature_error <= 1e-6 + 1e-7 * reference_temperature.abs(),
            "{family:?}/gamma={gamma} T(K) z={z}: actual={temperature} reference={reference_temperature}"
        );
        min_temperature = min_temperature.min(temperature);
        max_temperature = max_temperature.max(temperature);
        assert!(temperature.is_finite() && temperature > 0.0);
        assert!(
            actual_composition
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
        );

        let species_sum_error =
            (actual_composition.sum() - INITIAL_MASS_FRACTIONS.iter().sum::<f64>()).abs();
        assert!(
            species_sum_error <= 2e-9,
            "{family:?}/gamma={gamma} mass-fraction sum error at z={z}: {species_sum_error:.3e}"
        );
        let element_mass_fraction = [
            actual_composition.sum() * 2.0 * 12.011 / 60.052,
            actual_composition.sum() * 4.0 * 1.008 / 60.052,
            actual_composition.sum() * 2.0 * 15.999 / 60.052,
        ];
        for (actual, atom_mass) in element_mass_fraction.iter().zip([24.022, 4.032, 31.998]) {
            let expected_inventory =
                INITIAL_MASS_FRACTIONS.iter().sum::<f64>() * atom_mass / 60.052;
            max_element_balance = max_element_balance.max((actual - expected_inventory).abs());
        }

        let expected_heat_source = model.chemical_heat_release(&actual_composition, gamma);
        let actual_heat_source = snapshot.values.row(row).iter().copied().collect::<Vec<_>>();
        let actual_heat_source = task
            .heat_release
            .eval_expression(&variable_names, &actual_heat_source);
        let heat_error = (actual_heat_source - expected_heat_source).abs();
        max_heat_source_error = max_heat_source_error.max(heat_error);
        assert!(
            heat_error <= 1e-8 + 1e-10 * expected_heat_source.abs(),
            "{family:?}/gamma={gamma} heat source at z={z}: actual={actual_heat_source:.16e} reference={expected_heat_source:.16e} error={heat_error:.3e}"
        );
        let expected_q_rhs = LENGTH * MASS_FLUX * CP / LAMBDA * actual_q
            - LENGTH.powi(2) / T_SCALE * expected_heat_source;
        let actual_state = snapshot.values.row(row).iter().copied().collect::<Vec<_>>();
        let actual_q_rhs = task.solver.eq_system[1].eval_expression(&variable_names, &actual_state);
        let q_rhs_error = (actual_q_rhs - expected_q_rhs).abs();
        max_energy_rhs_error = max_energy_rhs_error.max(q_rhs_error);
        assert!(
            q_rhs_error <= 1e-9 + 1e-12 * expected_q_rhs.abs(),
            "{family:?}/gamma={gamma} energy RHS at z={z}: actual={actual_q_rhs:.16e} reference={expected_q_rhs:.16e} error={q_rhs_error:.3e}"
        );
    }
    assert!(
        max_element_balance <= 1e-14,
        "element balance drift={max_element_balance:.3e}"
    );
    println!(
        "TW-02 family={family:?} gamma={gamma:+.0} points={} max_dC={:.3e} max_dTeta={:.3e} max_dq={:.3e} max_dT_K={:.3e} max_dQchem={max_heat_source_error:.3e} max_energy_rhs_error={max_energy_rhs_error:.3e} max_element_balance={max_element_balance:.3e} T_range={min_temperature:.3}..{max_temperature:.3} K",
        snapshot.axis.len(),
        errors[0],
        errors[1],
        errors[2],
        errors[3]
    );
}

#[test]
fn tw02_component_enthalpies_obey_hess_law_and_known_reaction_heats() {
    let model = Tw02Model::new(1.0);
    let expected = [
        40.0, -15.0, 40.0, -50.0, 65.0, -45.0, 65.0, -50.0, 40.0, -30.0, 55.0, -115.0, 15.0, 40.0,
        10.0, 35.0, 35.0, 10.0, 25.0, 15.0,
    ];
    assert_eq!(model.reaction_heats_j_per_mol.as_slice(), &expected);

    let mut ring_cycle = [0.0; 20];
    ring_cycle[..12].fill(1.0);
    let mut alternate_path_cycle = [0.0; 20];
    alternate_path_cycle[0..4].fill(1.0);
    alternate_path_cycle[12] = -1.0;
    for cycle in [ring_cycle, alternate_path_cycle] {
        let mut species_delta = [0.0; SPECIES_COUNT];
        for (reaction, (from, to, _)) in REACTIONS.iter().copied().enumerate() {
            species_delta[from] -= cycle[reaction];
            species_delta[to] += cycle[reaction];
        }
        assert!(species_delta.iter().all(|value| *value == 0.0));
        let cycle_heat = model
            .reaction_heats_j_per_mol
            .iter()
            .zip(cycle)
            .map(|(heat, flow)| heat * flow)
            .sum::<f64>();
        assert_eq!(cycle_heat, 0.0);
    }

    let shifted = ENTHALPIES_J_PER_MOL.map(|enthalpy| enthalpy + 1234.5);
    for (reaction, (from, to, _)) in REACTIONS.iter().copied().enumerate() {
        assert_eq!(
            shifted[from] - shifted[to],
            model.reaction_heats_j_per_mol[reaction]
        );
    }
}

#[test]
fn tw02_14_state_exponential_has_zero_heat_and_gamma_metamorphic_properties() {
    let baseline_matrix = Tw02Model::new(0.0);
    for z in PROFILE_GRID {
        let baseline = baseline_matrix.state(z);
        let (tw01_c, tw01_theta, tw01_q, _) = tw01_profile(z);
        assert!((baseline[0] - tw01_theta).abs() <= 2e-13);
        assert!((baseline[1] - tw01_q).abs() <= 2e-13);
        for species in 0..SPECIES_COUNT {
            assert!((baseline[species + 2] - tw01_c[species]).abs() <= 2e-13);
        }
    }

    let unit = Tw02Model::new(1.0);
    for gamma in [10.0, -1.0] {
        let scaled = Tw02Model::new(gamma);
        for z in PROFILE_GRID {
            let reference = baseline_matrix.state(z);
            let one = unit.state(z);
            let actual = scaled.state(z);
            for species in 0..SPECIES_COUNT {
                assert!((actual[species + 2] - reference[species + 2]).abs() <= 2e-13);
            }
            for state in [0, 1] {
                let predicted = reference[state] + gamma * (one[state] - reference[state]);
                assert!((actual[state] - predicted).abs() <= 2e-12);
            }
        }
    }
}

#[test]
fn tw02_kinetic_and_heat_source_assembly_matches_independent_matrix() {
    let model = Tw02Model::new(1.0);
    let mut task = task_for(ReactorIvpSolverFamily::Lsode2, 1.0, &model);
    task.setup_condensed_ivp()
        .expect("TW-02 should assemble through existing KiThe path");
    assert_eq!(task.thermal_effects, model.reaction_heats_j_per_mol);

    for z in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let state = model.state(z);
        let mut named_state = vec![0.0; SPECIES_COUNT + 2];
        named_state[0] = state[0];
        named_state[1] = state[1];
        for (state_species_index, species) in task.kindata.substances.iter().enumerate() {
            let canonical = species.strip_prefix('X').unwrap().parse::<usize>().unwrap();
            named_state[state_species_index + 2] = state[canonical + 2];
        }
        let variables = task
            .solver
            .unknowns
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>();
        for (reaction, (from, _, _)) in REACTIONS.iter().copied().enumerate() {
            let equation = &task.kindata.vec_of_equations[reaction];
            let actual_rate = task.map_eq_rate[equation].eval_expression(&variables, &named_state);
            let expected_rate = (0..SPECIES_COUNT)
                .map(|species| model.reaction_rate_map[(reaction, species)] * state[species + 2])
                .sum::<f64>();
            assert!(
                (actual_rate - expected_rate).abs() <= 1e-10 + 1e-12 * expected_rate.abs(),
                "reaction={} reactant={from} rate actual={actual_rate} expected={expected_rate}",
                reaction + 1
            );
            let actual_heat_contribution = task.thermal_effects[reaction] * actual_rate;
            let expected_heat_contribution =
                model.reaction_heats_j_per_mol[reaction] * expected_rate;
            assert!(
                (actual_heat_contribution - expected_heat_contribution).abs()
                    <= 1e-10 + 1e-12 * expected_heat_contribution.abs(),
                "reaction={} heat contribution actual={actual_heat_contribution} expected={expected_heat_contribution}",
                reaction + 1
            );
        }
        let actual_heat = task.heat_release.eval_expression(&variables, &named_state);
        let expected_heat =
            model.chemical_heat_release(&state.rows(2, SPECIES_COUNT).into_owned(), 1.0);
        assert!((actual_heat - expected_heat).abs() <= 1e-8 + 1e-10 * expected_heat.abs());
        let actual_q_rhs = task.solver.eq_system[1].eval_expression(&variables, &named_state);
        let expected_q_rhs = (0..SPECIES_COUNT + 2)
            .map(|column| model.full_system[(1, column)] * state[column])
            .sum::<f64>();
        assert!((actual_q_rhs - expected_q_rhs).abs() <= 1e-10 + 1e-12 * expected_q_rhs.abs());
    }
}

#[test]
fn tw02_all_gammas_match_analytic_solution_on_all_solver_families() {
    for gamma in GAMMAS {
        let model = Tw02Model::new(gamma);
        for family in [
            ReactorIvpSolverFamily::Lsode2,
            ReactorIvpSolverFamily::Bdf,
            ReactorIvpSolverFamily::Radau,
        ] {
            let mut task = task_for(family, gamma, &model);
            let snapshot = task.solve_condensed_ivp().unwrap_or_else(|error| {
                panic!("TW-02 family={family:?} gamma={gamma:+.0} failed: {error}")
            });
            assert_profile(&snapshot, &task, &model, gamma, family);
        }
    }
}

#[test]
fn tw02_frozen_analytical_profiles_match_matrix_exponential() {
    let mut lines = FROZEN_PROFILE.lines();
    let header = lines.next().expect("TW-02 frozen profile header");
    assert!(header.starts_with("z,X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11"));
    let mut row_count = 0;
    for (line, expected_z) in lines.zip(PROFILE_GRID) {
        let values = line
            .split(',')
            .map(|cell| cell.parse::<f64>().expect("numeric TW-02 fixture cell"))
            .collect::<Vec<_>>();
        assert_eq!(values.len(), 1 + SPECIES_COUNT + 2 * GAMMAS.len());
        assert!((values[0] - expected_z).abs() <= f64::EPSILON);
        for (gamma_index, gamma) in GAMMAS.iter().copied().enumerate() {
            let state = Tw02Model::new(gamma).state(expected_z);
            for species in 0..SPECIES_COUNT {
                assert!(
                    (values[species + 1] - state[species + 2]).abs() <= 2e-13,
                    "frozen X{species} gamma={gamma} z={expected_z}"
                );
            }
            let theta_column = 1 + SPECIES_COUNT + 2 * gamma_index;
            assert!((values[theta_column] - state[0]).abs() <= 2e-13);
            assert!((values[theta_column + 1] - state[1]).abs() <= 2e-13);
        }
        row_count += 1;
    }
    assert_eq!(row_count, PROFILE_GRID.len());
}
